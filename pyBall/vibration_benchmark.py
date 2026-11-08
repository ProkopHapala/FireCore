"""Export MMFF vibration Hessians as documented sparse .npz benchmarks for external solver repos."""
from __future__ import annotations

import json
import os
import time
from dataclasses import asdict, dataclass, field
from typing import Any

import numpy as np

from pyBall import FTIR, MMFF

try:
    import scipy.sparse as sp
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False


@dataclass
class BenchmarkMeta:
    name: str
    description: str
    source_xyz: str
    natoms: int
    ndof: int
    dx: float
    n_shells: int
    max_neigh: int
    rigid_shift: float
    relax_steps: int
    nnz: int
    n_negative_eigvals_raw: int
    n_vibrational_modes: int
    omega_min_vib: float
    omega_max: float
    eigh_seconds: float | None = None
    units: dict = field(default_factory=lambda: {
        "length": "Angstrom",
        "energy": "eV",
        "mass": "amu",
        "Hessian": "eV/Ang^2",
        "omega": "sqrt(eV/amu)/Ang — internal MMFF units; NOT cm^-1",
    })
    notes: list[str] = field(default_factory=list)

    def to_json(self) -> str:
        return json.dumps(asdict(self), indent=2)


def reconstruct_dense_from_blocks(neigh_idx: np.ndarray, neigh_count: np.ndarray, blocks: np.ndarray) -> np.ndarray:
    natoms, max_neigh = neigh_idx.shape
    dim = natoms * 3
    H = np.zeros((dim, dim), dtype=np.float64)
    for p in range(natoms):
        for j in range(int(neigh_count[p])):
            o = int(neigh_idx[p, j])
            if o < 0:
                continue
            H[o * 3:(o + 1) * 3, p * 3:(p + 1) * 3] = blocks[p, j]
    return 0.5 * (H + H.T)


def mass_weighted_spectrum(H: np.ndarray, M: np.ndarray, *, vib_floor: float = 1e-4, check_minimum: bool = False):
    """Return (eigvals_mw, omegas, n_negative_raw, omega_min_vib)."""
    n_nodes = M.shape[0] // 3
    m = np.diag(M).reshape(n_nodes, 3)[:, 0]
    m_sqrt = np.repeat(np.sqrt(m), 3)
    m_inv_sqrt = 1.0 / m_sqrt
    Hmw = (m_inv_sqrt[:, None] * H) * m_inv_sqrt[None, :]
    w = np.linalg.eigvalsh(Hmw)
    n_neg = int(np.sum(w < -1e-8))
    if check_minimum:
        from pyBall.FFfit_utils import assert_harmonic_spectrum_at_minimum
        om_cm = np.sign(w) * 521.5 * np.sqrt(np.abs(w))
        assert_harmonic_spectrum_at_minimum(om_cm, ctx="mass_weighted_spectrum: ")
    omegas = np.sqrt(np.maximum(w, 0.0))
    vib = omegas[omegas > vib_floor]
    omega_min_vib = float(vib.min()) if vib.size else 0.0
    return w, omegas, n_neg, omega_min_vib


def sparse_csr_arrays(A_csr) -> dict[str, np.ndarray]:
    if not _HAS_SCIPY:
        raise ImportError("scipy required for sparse_csr_arrays")
    A_csr = A_csr.tocsr()
    return {
        "K_csr_data": A_csr.data.astype(np.float64),
        "K_csr_indices": A_csr.indices.astype(np.int32),
        "K_csr_indptr": A_csr.indptr.astype(np.int32),
        "K_csr_shape": np.array(A_csr.shape, dtype=np.int32),
    }


def compute_mmff_hessian_bundle(
    xyz_path: str,
    *,
    dx: float = 1e-4,
    n_shells: int = 2,
    max_neigh: int = 64,
    rigid_shift: float = 1e6,
    relax: bool = True,
    relax_steps: int = 2000,
    vib_floor: float = 1e-4,
    cwd: str | None = None,
) -> dict[str, Any]:
    """Load XYZ/MOL2, optional MMFF relax, return dense+sparse Hessian data."""
    if not _HAS_SCIPY:
        raise ImportError("scipy is required for vibration benchmark export")
    orig = os.getcwd()
    if cwd:
        os.chdir(cwd)
    try:
        MMFF.clear()
        MMFF.setVerbosity(0, 0)
        MMFF.init(xyz_name=os.path.abspath(xyz_path), bEpairs=False, bMMFF=True, nPBC=(0, 0, 0))
        MMFF.setSwitches(NonBonded=-1, MMFF=+1, Angles=+1, SurfAtoms=-1, GridFF=-1, PiSigma=-1, PiPiI=-1)
        MMFF.getBuffs()
        if relax:
            MMFF.set_opt(dt_max=0.05, dt_min=0.01, damp_max=0.1, cvf_min=-0.1, cvf_max=0.1)
            MMFF.run(nstepMax=relax_steps, dt=0.05, Fconv=1e-4, damping=0.1, ialg=2, omp=False)
            MMFF.getBuffs()
        natoms = MMFF.natoms
        pos = MMFF.apos[:natoms].copy()
        inds = np.arange(natoms, dtype=np.int32)
        H_full = MMFF.getHessian3Nx3N(inds, dx=dx)
        H_full = 0.5 * (H_full + H_full.T)
        M = FTIR.get_mass_matrix(MMFF, natoms)
        neigh_idx, neigh_count = MMFF.buildNeighShells(
            MMFF.neighs.copy()[:natoms], n_shells=n_shells, max_neigh=max_neigh, include_self=True
        )
        blocks = MMFF.getHessianSparseBlocks(neigh_idx, neigh_count, dx=dx)
        if np.isnan(blocks).any() or np.isinf(blocks).any():
            raise ValueError("NaN/Inf in sparse Hessian blocks")
        H_recon = reconstruct_dense_from_blocks(neigh_idx, neigh_count, blocks)
        H_proj = FTIR.project_rigid_modes(H_full, M, pos, shift=rigid_shift)
        H_sparse_raw = FTIR.build_sparse_hessian_from_blocks(neigh_idx, neigh_count, blocks, symmetrize=True)
        H_sparse_proj = FTIR.prepare_sparse_hessian(neigh_idx, neigh_count, blocks, M, pos, shift=rigid_shift)
        w_raw, om_raw, n_neg_raw, _ = mass_weighted_spectrum(H_full, M, vib_floor=vib_floor, check_minimum=True)
        w_proj, om_proj, n_neg_proj, om_min_vib = mass_weighted_spectrum(H_proj, M, vib_floor=vib_floor)
        elements = []
        for i in range(natoms):
            tn = MMFF.getTypeName(i)
            elements.append((tn.split("_")[0] or "C"))
        mass_diag = np.diag(M).astype(np.float64)
        t0 = time.perf_counter()
        _, _, _, _ = mass_weighted_spectrum(H_proj, M, vib_floor=vib_floor)
        eigh_s = time.perf_counter() - t0
        recon_err = float(np.max(np.abs(H_recon - H_full)))
        return {
            "natoms": natoms,
            "ndof": natoms * 3,
            "pos": pos,
            "elements": np.array(elements, dtype=object),
            "mass_diag": mass_diag,
            "neigh_idx": neigh_idx.astype(np.int32),
            "neigh_count": neigh_count.astype(np.int32),
            "blocks": blocks.astype(np.float64),
            "H_full": H_full,
            "H_projected": H_proj,
            "H_recon_from_blocks": H_recon,
            "H_sparse_projected": H_sparse_proj,
            "recon_vs_full_max_abs": recon_err,
            "eigvals_mw_raw": w_raw,
            "eigvals_mw_projected": w_proj,
            "omegas_modes_raw": om_raw,
            "omegas_modes_projected": om_proj,
            "n_negative_raw": n_neg_raw,
            "n_negative_projected": n_neg_proj,
            "omega_min_vib": om_min_vib,
            "eigh_seconds": eigh_s,
            "dx": dx,
            "n_shells": n_shells,
            "max_neigh": max_neigh,
            "rigid_shift": rigid_shift,
            "relax_steps": relax_steps if relax else 0,
            "source_xyz": os.path.abspath(xyz_path),
        }
    finally:
        os.chdir(orig)


def export_benchmark_npz(
    out_path: str,
    bundle: dict[str, Any],
    meta: BenchmarkMeta,
    *,
    include_dense: bool = True,
    dense_max_ndof: int = 3000,
) -> None:
    """Write self-documented .npz for external sparse solver benchmarks."""
    if not _HAS_SCIPY:
        raise ImportError("scipy required")
    ndof = int(bundle["ndof"])
    payload: dict[str, Any] = {
        "meta_json": meta.to_json(),
        "name": meta.name,
        "natoms": np.int32(bundle["natoms"]),
        "ndof": np.int32(ndof),
        "pos": bundle["pos"],
        "elements": bundle["elements"],
        "mass_diag": bundle["mass_diag"],
        "neigh_idx": bundle["neigh_idx"],
        "neigh_count": bundle["neigh_count"],
        "blocks": bundle["blocks"],
        "dx": np.float64(bundle["dx"]),
        "n_shells": np.int32(bundle["n_shells"]),
        "max_neigh": np.int32(bundle["neigh_idx"].shape[1]),
        "rigid_shift": np.float64(bundle["rigid_shift"]),
        "omegas_modes_projected": bundle["omegas_modes_projected"],
        "eigvals_mw_projected": bundle["eigvals_mw_projected"],
        "n_negative_raw": np.int32(bundle["n_negative_raw"]),
        "n_negative_projected": np.int32(bundle["n_negative_projected"]),
        "omega_min_vib": np.float64(bundle["omega_min_vib"]),
        "recon_vs_full_max_abs": np.float64(bundle["recon_vs_full_max_abs"]),
        "eigh_seconds": np.float64(bundle["eigh_seconds"]),
    }
    payload.update(sparse_csr_arrays(bundle["H_sparse_projected"]))
    if include_dense and ndof <= dense_max_ndof:
        payload["H_dense_full"] = bundle["H_full"]
        payload["H_dense_projected"] = bundle["H_projected"]
        payload["H_dense_recon_blocks"] = bundle["H_recon_from_blocks"]
    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    np.savez_compressed(out_path, **payload)
    print(f"[vibration_benchmark] wrote {out_path} natoms={bundle['natoms']} ndof={ndof} nnz={meta.nnz} eigh={meta.eigh_seconds:.3f}s")
