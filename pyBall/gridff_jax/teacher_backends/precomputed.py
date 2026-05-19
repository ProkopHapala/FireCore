"""Precomputed teacher backend.

Reads ``(positions, pose_params, energies, forces)`` from a single ``.npz``
file produced offline. The user (or their student) generates this with any
DFT or QM code — VASP, Psi4, Gaussian, NWChem, ORCA, Q-Chem, whatever —
then drops the npz here. The backend never calls a calculator; it just
looks up each pose in the cached array and returns the stored values.

Why this exists
---------------
For substrate classes outside MAD-SURF's training distribution (NaCl, KBr,
oxides, …), there is no convenient runtime MLIP. The most general way to
inject ground truth is to pre-compute it once and supply the file.

NPZ format
----------
Required arrays:

  pose_params   (n_pose, 7)          float64  [u, v, z, qw, qx, qy, qz]
  energies      (n_pose,)            float64  eV  (interaction energy preferred)

Optional but recommended:

  forces        (n_pose, n_atom, 3)  float64  eV/Å (forces on molecule atoms)
                                              If absent, the backend returns
                                              zeros and sets metadata flag
                                              ``forces_available=False`` so the
                                              driver can auto-switch to
                                              ``--fit-mode energy``.

Also optional:

  positions          (n_pose, n_atom, 3) float64 Å — saved for traceability;
                                                   not used by the backend
                                                   (positions are derivable
                                                   from pose_params).
  adsorbate_symbols  (n_atom,)            str    — sanity check vs the pose batch
  metadata           dict (0-d object array)     — provenance info

Lookup is by *exact* match on the ``pose_params`` array (positions are not
used because they are deterministically derivable from pose_params + the
adsorbate geometry). Tolerance is configurable via ``metadata['atol']``
(default 1e-8). Missing-pose errors are raised with the offending row
indices so the user can debug their generator script.

Configuration
-------------
  kind            = "precomputed"
  metadata["npz_path"] = "/abs/path/to/teacher.npz"  (required)
  metadata["atol"]     = 1e-8  (optional)

(Reusing ``metadata`` keeps the existing ``TeacherBackendConfig`` schema
intact — no new fields needed.)
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from ..interfaces import TeacherBackend, TeacherResult


class PrecomputedTeacherBackend(TeacherBackend):
    def __init__(self, config, model_factory=None):
        self.config = config
        self.model_factory = model_factory   # unused but required by registry signature
        npz_path = (config.metadata or {}).get("npz_path")
        if not npz_path:
            raise ValueError(
                "precomputed teacher requires metadata['npz_path'] in the "
                "TeacherBackendConfig"
            )
        self.npz_path = Path(npz_path)
        if not self.npz_path.exists():
            raise FileNotFoundError(f"precomputed teacher npz not found: {self.npz_path}")
        self.atol = float((config.metadata or {}).get("atol", 1e-8))
        self._loaded = None

    def _load(self):
        if self._loaded is not None:
            return self._loaded
        d = np.load(self.npz_path, allow_pickle=True)
        required = ("pose_params", "energies")
        missing = [k for k in required if k not in d.files]
        if missing:
            raise KeyError(
                f"precomputed npz {self.npz_path} missing required arrays: {missing}. "
                f"Available: {list(d.files)}"
            )
        pose_params = np.asarray(d["pose_params"], dtype=float)
        energies = np.asarray(d["energies"], dtype=float)
        if pose_params.shape[0] != energies.shape[0]:
            raise ValueError(
                f"precomputed npz shape mismatch: pose_params {pose_params.shape}, "
                f"energies {energies.shape}"
            )

        # Forces are optional. If absent, fill with zeros and flag in metadata
        # so the driver can detect and switch to --fit-mode energy.
        forces_available = "forces" in d.files
        if forces_available:
            forces = np.asarray(d["forces"], dtype=float)
            if forces.shape[0] != pose_params.shape[0]:
                raise ValueError(
                    f"precomputed npz forces row count ({forces.shape[0]}) does "
                    f"not match energies/pose_params ({pose_params.shape[0]})"
                )
        else:
            # Need a sensible shape; infer n_atom from positions if present,
            # otherwise fall back to 1 atom (shape will be reshaped on demand
            # in evaluate_batch from poses.adsorbate.natoms).
            if "positions" in d.files:
                n_atom = int(np.asarray(d["positions"]).shape[1])
            else:
                n_atom = 1
            forces = np.zeros((pose_params.shape[0], n_atom, 3), dtype=float)
        # Sanity: adsorbate symbols if available
        ads_syms = None
        if "adsorbate_symbols" in d.files:
            ads_syms = [str(s) for s in d["adsorbate_symbols"].tolist()]
        self._loaded = {
            "pose_params": pose_params,
            "energies": energies,
            "forces": forces,
            "forces_available": forces_available,
            "adsorbate_symbols": ads_syms,
            "metadata": dict(d["metadata"].item()) if "metadata" in d.files else {},
        }
        return self._loaded

    def evaluate_batch(self, density, poses) -> TeacherResult:
        cache = self._load()
        ref_pose_params = cache["pose_params"]
        ref_energies = cache["energies"]
        ref_forces = cache["forces"]

        # Adsorbate sanity check
        ads_syms = cache["adsorbate_symbols"]
        if ads_syms is not None and ads_syms != list(poses.adsorbate.symbols):
            raise ValueError(
                f"adsorbate mismatch: npz has {ads_syms}, "
                f"pose batch has {list(poses.adsorbate.symbols)}"
            )

        query = np.asarray(poses.pose_params, dtype=float)
        if query.shape[1] != ref_pose_params.shape[1]:
            raise ValueError(
                f"pose_params width mismatch: query {query.shape[1]}, "
                f"npz {ref_pose_params.shape[1]}"
            )

        n_query = query.shape[0]
        n_atom = poses.adsorbate.natoms
        out_E = np.zeros(n_query, dtype=float)
        # Validate force shape vs the pose batch's adsorbate. Previously we
        # used ref_forces.shape[1:] blindly, which lets a (n_pose, K!=n_atom, 3)
        # NPZ slip through and feed wrong force tensors into the fit. Now we
        # require an exact (n_pose, n_atom, 3) shape match — the user gets a
        # clear error before any silent corruption.
        if cache["forces_available"]:
            if ref_forces.ndim != 3 or ref_forces.shape[1:] != (n_atom, 3):
                raise ValueError(
                    f"precomputed npz forces have shape {ref_forces.shape}, "
                    f"but the pose batch adsorbate has {n_atom} atoms (expected "
                    f"(n_pose, {n_atom}, 3)). Regenerate the npz with the right "
                    f"adsorbate or check that adsorbate_symbols matches."
                )
            out_F = np.zeros((n_query, n_atom, 3), dtype=float)
        else:
            out_F = np.zeros((n_query, n_atom, 3), dtype=float)
        missing = []
        # Vectorised exact-row lookup with tolerance
        for i in range(n_query):
            diff = np.abs(ref_pose_params - query[i][None, :])
            hits = np.where(np.all(diff <= self.atol, axis=1))[0]
            if hits.size == 0:
                missing.append(i)
                continue
            j = int(hits[0])
            out_E[i] = ref_energies[j]
            out_F[i] = ref_forces[j]

        if missing:
            raise LookupError(
                f"precomputed teacher: {len(missing)}/{n_query} poses not found "
                f"in {self.npz_path} (first missing query row indices: "
                f"{missing[:5]}). Regenerate the cache or increase atol "
                f"(currently {self.atol:g})."
            )

        return TeacherResult(
            energies=out_E,
            forces=out_F,
            metadata={
                "backend": "precomputed",
                "npz_path": str(self.npz_path),
                "atol": self.atol,
                "n_cached": int(len(ref_pose_params)),
                "n_queried": int(n_query),
                "forces_available": bool(cache["forces_available"]),
                "source_metadata": cache["metadata"],
            },
        )
