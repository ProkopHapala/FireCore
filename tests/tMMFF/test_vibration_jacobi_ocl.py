#!/usr/bin/env python3
"""Stage 3 GPU Jacobi parity: one block-Jacobi iteration vs CPU (M-S3) on adamantane."""
import argparse
import os
import sys
from pathlib import Path

import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))

try:
    import scipy.sparse as sp
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

from pyBall import FTIR
from pyBall.OCL.VibSolver import VibSolver, solve_dynamic_stiffness_gpu

FIXTURE_ROOT = Path(__file__).resolve().parent.parent / "tSiNCs/fixtures/vibration_parallel"
ADAMANTANE_NPZ = FIXTURE_ROOT / "hessian_mmff/adamantane_mmff_dense.npz"


def load_adamantane():
    if not ADAMANTANE_NPZ.is_file():
        raise FileNotFoundError(f"Missing {ADAMANTANE_NPZ}; run bootstrap_vibration_parallel_fixtures.py")
    d = np.load(ADAMANTANE_NPZ, allow_pickle=True)
    return d["H"], d["M"], d["pos"]


def mass_weighted_mode_freqs(H, M):
    m_sqrt = np.sqrt(np.diag(M))
    m_inv_sqrt = 1.0 / m_sqrt
    K_mw = (m_inv_sqrt[:, None] * H) * m_inv_sqrt[None, :]
    w = np.linalg.eigvalsh(K_mw)
    return np.sqrt(np.maximum(w, 0.0))


def cpu_block_jacobi_step(A_csr, f, n_atoms):
    """One block-Jacobi iteration (zero initial guess, no momentum)."""
    rn, rb, rc, di = FTIR.pack_block_rows_from_sparse(A_csr, n_atoms)
    u = np.zeros(len(f), dtype=np.complex128)
    u_new = np.zeros(len(f), dtype=np.complex128)
    f = np.asarray(f, dtype=np.complex128).ravel()
    for i in range(n_atoms):
        si = slice(i * 3, (i + 1) * 3)
        acc = f[si].copy()
        for t in range(int(rc[i])):
            j = int(rn[i, t])
            sj = slice(j * 3, (j + 1) * 3)
            acc -= rb[i, t] @ u[sj]
        u_new[si] = di[i] @ acc
    return u_new


def main():
    ap = argparse.ArgumentParser(description="GPU vib_jacobi one-step parity (M-S3)")
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--shift", type=float, default=1e6)
    ap.add_argument("--mode-index", type=int, default=0)
    ap.add_argument("--omega-offset-frac", type=float, default=0.15, dest="omega_offset_frac")
    ap.add_argument("--solve-tol", type=float, default=1e-4, help="rel_err one-step GPU vs CPU")
    ap.add_argument("--stabilize", type=float, default=1e-6)
    args = ap.parse_args()

    if not _HAS_SCIPY:
        raise ImportError("scipy required")

    print("=== Stage 3: pyOpenCL vib_jacobi one-step vs CPU (adamantane) ===")
    H_raw, M, pos = load_adamantane()
    natoms = pos.shape[0]
    om_modes = mass_weighted_mode_freqs(FTIR.apply_rigid_mode_shift(H_raw, M, pos, shift=args.shift), M)
    omega_mode = float(om_modes[6 + args.mode_index])
    omega = omega_mode * (1.0 + args.omega_offset_frac)
    print(f"  test omega={omega:.6f} (mode {omega_mode:.6f} + offset {args.omega_offset_frac})")

    H_sparse = FTIR.apply_rigid_mode_shift(sp.csr_matrix(H_raw), M, pos, shift=args.shift)
    f = np.zeros(natoms * 3, dtype=np.float64)
    f[0] = 1.0
    A = FTIR._dynamic_stiffness_sparse(H_sparse, M, omega, eta=args.eta, stabilize=args.stabilize)

    u_cpu = cpu_block_jacobi_step(A, f, natoms)
    solver = VibSolver()
    solver.realloc(natoms)
    u_gpu = solve_dynamic_stiffness_gpu(
        H_sparse, M, omega, args.eta, f, n_iter=1, b_mix=0.0,
        stabilize=args.stabilize, solver=solver
    )
    rel_err = np.linalg.norm(u_gpu - u_cpu) / max(np.linalg.norm(u_cpu), 1e-30)
    print(f"  one-step rel_err GPU vs CPU = {rel_err:.3e} (tol={args.solve_tol})")
    if rel_err > args.solve_tol:
        raise AssertionError(f"M-S3 failed: rel_err={rel_err:.3e} > {args.solve_tol}")
    print("M-S3 PASSED")


if __name__ == "__main__":
    main()
