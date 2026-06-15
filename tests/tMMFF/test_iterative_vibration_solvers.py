#!/usr/bin/env python3
"""Compare iterative linear solvers on adamantane A(omega)u=f; save spectrum parity plots."""
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(__file__), "../.."))
from pyBall import FTIR

try:
    import scipy.sparse as sp
    from scipy.sparse.linalg import spsolve
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

FIXTURE = Path(__file__).resolve().parent.parent / "tSiNCs/fixtures/vibration_parallel/hessian_mmff/adamantane_mmff_dense.npz"
OUTDIR = Path(__file__).resolve().parent / "OUT_vibration_solver"


def symmetric_excitation(pos):
    com = np.mean(pos, axis=0)
    dir_x = np.zeros((pos.shape[0], 3))
    dir_x[:, 0] = pos[:, 0] - com[0]
    return dir_x, np.ones(pos.shape[0])


def benchmark_single_omega(A, f, u_ref):
    results = {}
    results['spsolve'] = (spsolve(A, f.astype(np.complex128)), 1)
    try:
        A_real = sp.csr_matrix(np.real(A.toarray()), dtype=np.float64)
        w = np.linalg.eigvalsh(A_real.toarray())
        if w[0] > 0:
            u, n, _ = FTIR.solve_cg_spd(A_real, np.real(f), tol=1e-9)
            results['cg'] = (u, n)
        else:
            results['cg'] = (None, f'not SPD min_eig={w[0]:.3f}')
    except Exception as e:
        results['cg'] = (None, str(e))
    try:
        u, n, _ = FTIR.solve_gmres(A, f, restart=30, max_cycles=30, tol=1e-9)
        results['gmres'] = (u, n)
    except Exception as e:
        results['gmres'] = (None, str(e))
    try:
        u, n, _ = FTIR.solve_cgne(A, f, tol=1e-10, max_iter=3000)
        results['cgne'] = (u, n)
    except Exception as e:
        results['cgne'] = (None, str(e))
    try:
        u, n, _ = FTIR.solve_preconditioned_richardson(A, f, tau='auto', max_iter=5000, tol=1e-4)
        results['richardson'] = (u, n)
    except Exception as e:
        results['richardson'] = (None, str(e))
    try:
        u, n, _ = FTIR.solve_linear_jacobi_momentum(A, f, method='block', max_iter=500, tol=1e-6)
        results['jacobi'] = (u, n)
    except Exception as e:
        results['jacobi'] = (None, str(e))
    print(f"  single-omega benchmark (ref norm={np.linalg.norm(u_ref):.3e}):")
    for name, val in results.items():
        if val[0] is None:
            print(f"    {name:12s}: FAILED {val[1]}")
        else:
            u, n = val
            err = np.linalg.norm(u - u_ref) / max(np.linalg.norm(u_ref), 1e-30)
            print(f"    {name:12s}: iters={n:5d}  rel_err={err:.3e}")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fmin", type=float, default=0.01)
    ap.add_argument("--fmax", type=float, default=3.5)
    ap.add_argument("--nfreq", type=int, default=400)
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--shift", type=float, default=1e6)
    ap.add_argument("--outdir", type=str, default=str(OUTDIR))
    ap.add_argument("--benchmark-only", action="store_true")
    args = ap.parse_args()

    if not _HAS_SCIPY:
        raise ImportError("scipy required")
    if not FIXTURE.is_file():
        raise FileNotFoundError(f"Missing {FIXTURE}")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    d = np.load(FIXTURE, allow_pickle=True)
    H_raw, M, pos = d["H"], d["M"], d["pos"]
    H = FTIR.apply_rigid_mode_shift(sp.csr_matrix(H_raw), M, pos, shift=args.shift)
    dir_x, charges = symmetric_excitation(pos)
    omegas = np.linspace(args.fmin, args.fmax, args.nfreq)

    # single-omega solver comparison
    omega_test = 0.915
    A_test = FTIR._dynamic_stiffness_sparse(H, M, omega_test, eta=args.eta)
    f_test = np.zeros(H.shape[0], dtype=np.complex128)
    f_test[0] = 1.0
    u_ref = spsolve(A_test, f_test)
    print(f"=== Single-omega solver benchmark at omega={omega_test} ===")
    benchmark_single_omega(A_test, f_test, u_ref)

    if args.benchmark_only:
        return

    print(f"=== Spectrum scan ({args.nfreq} frequencies) ===")
    res_ref = FTIR.mechanical_greens_probing_sparse(
        H, M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges, solver='spsolve')
    res_gmres = FTIR.mechanical_greens_probing_sparse(
        H, M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges,
        solver='gmres', gmres_restart=30, gmres_cycles=40, jacobi_tol=1e-7, fail_loud=True)
    # lossless CG case: eta=0, smaller grid for speed
    omegas_cg = np.linspace(0.3, 0.7, 80)
    Hsp = H
    res_cg = FTIR.mechanical_greens_probing_sparse(
        Hsp, M, omegas_cg, eta=0.0, direction_vec=dir_x, charges=charges, solver='cg', jacobi_tol=1e-8)
    res_cg_ref = FTIR.mechanical_greens_probing_sparse(
        Hsp, M, omegas_cg, eta=0.0, direction_vec=dir_x, charges=charges, solver='spsolve')

    rel_gmres = np.max(np.abs(res_gmres["energy"] - res_ref["energy"])) / (np.max(res_ref["energy"]) + 1e-12)
    rel_cg = np.max(np.abs(res_cg["energy"] - res_cg_ref["energy"])) / (np.max(res_cg_ref["energy"]) + 1e-12)
    print(f"  GMRES spectrum max rel diff vs spsolve: {rel_gmres:.3e}")
    print(f"  CG (eta=0) spectrum max rel diff vs spsolve: {rel_cg:.3e}")

    fig, axes = plt.subplots(2, 1, figsize=(10, 8))
    ax = axes[0]
    ax.plot(res_ref["omega"], res_ref["energy"], 'g-', lw=1.5, label='spsolve (ref)')
    ax.plot(res_gmres["omega"], res_gmres["energy"], 'r--', lw=1.2, alpha=0.85, label=f'GMRES (max rel diff {rel_gmres:.1e})')
    ax.set_xlabel(r'$\omega$'); ax.set_ylabel(r'$|u|^2$ response')
    ax.set_title(f'Complex A(ω) with η={args.eta}: GMRES vs direct sparse solve')
    ax.legend(); ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(res_cg_ref["omega"], res_cg_ref["energy"], 'g-', lw=1.5, label='spsolve')
    ax.plot(res_cg["omega"], res_cg["energy"], 'b--', lw=1.2, label=f'CG SPD η=0 (max rel diff {rel_cg:.1e})')
    ax.set_xlabel(r'$\omega$'); ax.set_ylabel(r'$|u|^2$ response')
    ax.set_title('Lossless probe (η=0, ω below first mode): CG vs spsolve')
    ax.legend(); ax.grid(True, alpha=0.3)

    fig.tight_layout()
    p = outdir / "spectrum_iterative_solvers.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print(f"  saved plot: {p}")


if __name__ == "__main__":
    main()
