"""Generic ASE-calculator teacher backend.

Wraps any calculator with an ASE-compatible interface
(``calc.get_potential_energy()``, ``calc.get_forces()``) and exposes it as
a ``TeacherBackend``. Use this for any MLIP (MACE, ANI, NEP, AIMNet2,
schnetpack, …) and for runtime-DFT wrappers like ASE's VASP or Psi4
calculators — anything that ASE knows how to drive.

Configuration
-------------
Set in ``TeacherBackendConfig`` (or RunConfig.teacher_backend):

  kind             = "generic_calc"
  model_module     = "mace.calculators"     # importable module
  model_callable   = "MACECalculator"       # class or factory in that module
  model_path       = "/path/to/model.bin"   # passed as model_paths=... kwarg
  device           = "cuda"
  default_dtype    = "float64"
  interaction_energy = True                  # subtract slab + gas-phase mol
  teacher_tile     = "auto" or (NX, NY)
  metadata         = {                       # optional kwargs passed to ctor
      "ctor_kwargs": {"r_cutoff": 5.0, ...}
  }

The class instantiates the calculator once and caches it. Slab energies
and gas-phase molecule references are cached per tiling for efficiency.

Why this exists
---------------
``madsurf.py`` was previously the only teacher backend that actually ran
an ML potential, and it imported ``mace.calculators.MACECalculator``
unconditionally. That made the whole framework MACE-locked even though the
slab/tile/interaction-energy plumbing is generic ASE logic. This backend
extracts the generic plumbing so users can swap MLIPs without touching
the framework.
"""

from __future__ import annotations

import hashlib
import importlib
import math
import time

import numpy as np
from ase import Atoms

from ..interfaces import TeacherBackend, TeacherResult
from ..utils import quaternion_to_matrix


def _auto_tile(mol_positions: np.ndarray, cell: np.ndarray) -> tuple[int, int]:
    """Min tiling so the molecule fits without periodic-image self-interaction.

    Same convention as madsurf.py: ``ceil(extent / cell_length) + 1`` (the
    +1 buffer prevents image overlap from edge cases like rotated molecules).
    """
    cell = np.asarray(cell, dtype=float)
    mol_positions = np.asarray(mol_positions, dtype=float)
    tiles = []
    for axis in range(2):
        cell_length = float(np.linalg.norm(cell[axis]))
        mol_extent = float(mol_positions[:, axis].max() - mol_positions[:, axis].min())
        n = math.ceil(mol_extent / cell_length) + 1 if cell_length > 0.0 else 1
        tiles.append(max(n, 1))
    return (tiles[0], tiles[1])


class GenericASECalcBackend(TeacherBackend):
    """Generic ASE-calculator teacher backend.

    Subclass this to hardcode a specific calculator (see madsurf.py for an
    example using MACECalculator); or use it directly with
    ``model_module`` + ``model_callable`` set in the config.
    """

    def __init__(self, config, model_factory=None):
        self.config = config
        self.model_factory = model_factory
        self._callable = None
        self._calculator = None
        self._slab_energy = None
        self._slab_energy_tiled: dict[tuple[int, int], float] = {}
        self._molecule_reference_cache: dict[str, dict] = {}
        if model_factory is not None:
            self._callable = model_factory(config)
        elif getattr(config, "model_module", None) and getattr(config, "model_callable", None):
            module = importlib.import_module(config.model_module)
            self._callable = getattr(module, config.model_callable)

    # ----- override this to inject a specific calculator class -----
    def _build_calculator(self):
        """Instantiate the ASE calculator from the config.

        Default behaviour:
        1. If ``model_factory`` was supplied at __init__, use that as a
           callable that takes the config and returns a calculator.
        2. Otherwise, import ``config.model_callable`` from
           ``config.model_module`` and call it as
               cls(model_paths=config.model_path,
                   device=config.device,
                   default_dtype=config.default_dtype,
                   **(config.metadata.get("ctor_kwargs") or {}))
        """
        if self._callable is None:
            raise RuntimeError(
                "generic_calc backend needs either model_module+model_callable "
                "(in the TeacherBackendConfig) or an explicit model_factory."
            )
        ctor_kwargs = {}
        if getattr(self.config, "metadata", None):
            ctor_kwargs = dict(self.config.metadata.get("ctor_kwargs", {}))
        kwargs = {
            "model_paths": self.config.model_path,
            "device": self.config.device,
            "default_dtype": self.config.default_dtype,
            **ctor_kwargs,
        }
        # Construct the calculator. Do NOT swallow TypeError — that previously
        # masked real signature errors and could silently drop ``model_paths``,
        # instantiating the wrong calculator. If the calculator class needs a
        # different kwarg layout, the user must specify it explicitly via
        # ``metadata['ctor_kwargs']`` (which can include None to suppress a
        # default) or via a ``model_factory`` callable supplied at __init__.
        return self._callable(**kwargs)

    def _get_calculator(self):
        if self._calculator is None:
            self._calculator = self._build_calculator()
        return self._calculator

    @property
    def teacher_tile(self) -> tuple[int, int]:
        t = getattr(self.config, "teacher_tile", (1, 1))
        # Accept the documented 'auto' string (JSON configs use it) as well
        # as the canonical (0, 0) / (NX, NY) tuples. Maps to (0, 0) which
        # triggers automatic per-pose tile selection downstream.
        if isinstance(t, str):
            if t.strip().lower() == "auto":
                return (0, 0)
            raise ValueError(
                f"teacher_tile string must be 'auto' or a 'NX,NY' pair; got {t!r}"
            )
        return (int(t[0]), int(t[1]))

    # ----- slab / molecule reference plumbing -----
    def _slab_atoms(self, density) -> Atoms:
        return Atoms(
            symbols=density.symbols,
            positions=np.asarray(density.positions, dtype=float),
            cell=np.asarray(density.cell, dtype=float),
            pbc=[True, True, True],
        )

    def _tiled_slab_atoms(self, density, n_tile_x: int, n_tile_y: int) -> Atoms:
        cell = np.asarray(density.cell, dtype=float)
        positions = np.asarray(density.positions, dtype=float)
        symbols = list(density.symbols)
        tiled_positions = []
        tiled_symbols = []
        for ix in range(n_tile_x):
            for iy in range(n_tile_y):
                shift = ix * cell[0] + iy * cell[1]
                tiled_positions.append(positions + shift)
                tiled_symbols.extend(symbols)
        big_cell = cell.copy()
        big_cell[0] = big_cell[0] * n_tile_x
        big_cell[1] = big_cell[1] * n_tile_y
        return Atoms(
            symbols=tiled_symbols,
            positions=np.vstack(tiled_positions),
            cell=big_cell,
            pbc=[True, True, True],
        )

    def _resolve_tile(self, density, mol_positions: np.ndarray) -> tuple[int, int]:
        tx, ty = self.teacher_tile
        if tx <= 0 or ty <= 0:
            return _auto_tile(mol_positions, np.asarray(density.cell, dtype=float))
        return (tx, ty)

    def _get_slab_energy(self, density, calc, n_tile_x: int, n_tile_y: int) -> float:
        if n_tile_x == 1 and n_tile_y == 1:
            if self._slab_energy is None:
                slab = self._slab_atoms(density)
                slab.calc = calc
                self._slab_energy = float(slab.get_potential_energy())
            return self._slab_energy
        key = (n_tile_x, n_tile_y)
        if key not in self._slab_energy_tiled:
            slab = self._tiled_slab_atoms(density, n_tile_x, n_tile_y)
            slab.calc = calc
            self._slab_energy_tiled[key] = float(slab.get_potential_energy())
        return self._slab_energy_tiled[key]

    @staticmethod
    def _adsorbate_cache_key(adsorbate) -> str:
        """Cache key that fingerprints the adsorbate geometry + composition +
        charges so we never reuse a reference computed for a different
        conformer / protonation / charge state under the same name."""
        h = hashlib.sha1()
        h.update(str(adsorbate.name).encode())
        h.update(",".join(adsorbate.symbols).encode())
        h.update(np.ascontiguousarray(adsorbate.positions, dtype=np.float64).tobytes())
        if adsorbate.charges is not None:
            h.update(b"|q|")
            h.update(np.ascontiguousarray(adsorbate.charges, dtype=np.float64).tobytes())
        h.update(b"|anchor|")
        h.update(str(int(adsorbate.anchor_index)).encode())
        return h.hexdigest()[:24]

    def _molecule_reference(self, poses) -> dict:
        key = self._adsorbate_cache_key(poses.adsorbate)
        if key in self._molecule_reference_cache:
            return self._molecule_reference_cache[key]
        positions = np.asarray(poses.adsorbate.positions, dtype=float)
        positions = positions - positions.mean(axis=0, keepdims=True)
        cell = np.diag([40.0, 40.0, 40.0])
        kwargs = dict(
            symbols=poses.adsorbate.symbols,
            positions=positions + np.array([20.0, 20.0, 20.0]),
            cell=cell,
            pbc=[True, True, True],
        )
        # If the adsorbate carries fixed charges, pass them as initial_charges
        # — charge-aware calculators (DFT wrappers, polarizable models) use
        # this. Charge-blind calculators (MACE, ANI) ignore initial_charges.
        if poses.adsorbate.charges is not None:
            kwargs["charges"] = np.asarray(poses.adsorbate.charges, dtype=float)
        atoms = Atoms(**kwargs)
        calc = self._get_calculator()
        atoms.calc = calc
        self._molecule_reference_cache[key] = {
            "energy": float(atoms.get_potential_energy()),
            "forces": np.asarray(atoms.get_forces(), dtype=float),
            "ads_name": poses.adsorbate.name,
            "key": key,
        }
        return self._molecule_reference_cache[key]

    # ----- main entry point -----
    def evaluate_batch(self, density, poses) -> TeacherResult:
        if self._callable is not None and self._calculator is None:
            # If model_factory returned a *calculator*, use it directly.
            # If it returned a *class*, we'll lazily instantiate in _build_calculator.
            try:
                # Heuristic: if calling with no args fails, it's a class — let lazy path handle it.
                test = self._callable
                if hasattr(test, "get_potential_energy"):
                    self._calculator = test
            except Exception:
                pass

        calc = self._get_calculator()
        molecule_reference = self._molecule_reference(poses) if self.config.interaction_energy else None

        n_tile_x, n_tile_y = self._resolve_tile(density, poses.positions[0])
        use_tiling = n_tile_x > 1 or n_tile_y > 1

        if use_tiling:
            tiled_slab = self._tiled_slab_atoms(density, n_tile_x, n_tile_y)
            slab_symbols = list(tiled_slab.get_chemical_symbols())
            slab_positions = np.asarray(tiled_slab.get_positions(), dtype=float)
            slab_cell = np.asarray(tiled_slab.get_cell(), dtype=float)
            orig_cell = np.asarray(density.cell, dtype=float)
            centre_shift = 0.5 * (n_tile_x - 1) * orig_cell[0] + 0.5 * (n_tile_y - 1) * orig_cell[1]
        else:
            slab_symbols = list(density.symbols)
            slab_positions = np.asarray(density.positions, dtype=float)
            slab_cell = np.asarray(density.cell, dtype=float)
            centre_shift = np.zeros(3, dtype=float)

        slab_energy = self._get_slab_energy(density, calc, n_tile_x, n_tile_y)

        n_ads = poses.adsorbate.natoms
        energies: list[float] = []
        forces: list[np.ndarray] = []
        t0 = time.perf_counter()
        for ipose, positions in enumerate(poses.positions):
            mol_positions = np.asarray(positions, dtype=float) + centre_shift
            combined = Atoms(
                symbols=slab_symbols + poses.adsorbate.symbols,
                positions=np.vstack([slab_positions, mol_positions]),
                cell=slab_cell,
                pbc=[True, True, True],
            )
            combined.calc = calc
            total_energy = float(combined.get_potential_energy())
            total_forces = np.asarray(combined.get_forces(), dtype=float)
            if self.config.interaction_energy:
                if molecule_reference is None:
                    raise RuntimeError("interaction_energy=True but molecule reference was not built")
                energies.append(total_energy - slab_energy - molecule_reference["energy"])
            else:
                energies.append(total_energy)
            ads_forces = total_forces[-n_ads:, :]
            if self.config.interaction_energy:
                if molecule_reference is None:
                    raise RuntimeError("interaction_energy=True but molecule reference was not built")
                quat = np.asarray(poses.pose_params[ipose, 3:7], dtype=float)
                rot = quaternion_to_matrix(quat)
                molecule_forces = molecule_reference["forces"] @ rot.T
                ads_forces = ads_forces - molecule_forces
            forces.append(ads_forces)
        elapsed = time.perf_counter() - t0
        return TeacherResult(
            energies=np.asarray(energies, dtype=float),
            forces=np.asarray(forces, dtype=float),
            metadata={
                "backend": "generic_calc",
                "model_module": getattr(self.config, "model_module", None),
                "model_callable": getattr(self.config, "model_callable", None),
                "model_path": self.config.model_path,
                "device": self.config.device,
                "default_dtype": self.config.default_dtype,
                "interaction_energy": self.config.interaction_energy,
                "teacher_tile": [n_tile_x, n_tile_y],
                "teacher_eval_seconds": elapsed,
                "seconds_per_pose": elapsed / max(len(poses.positions), 1),
            },
        )
