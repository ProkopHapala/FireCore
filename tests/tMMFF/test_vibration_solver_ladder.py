#!/usr/bin/env python3
"""Staged validation ladder for sparse vibration linear solvers (M-S0..M-S2) + spectrum plots.

Fixtures: tests/tSiNCs/fixtures/vibration_parallel/ (gitignored; regenerate via bootstrap script).
Plots: tests/tMMFF/OUT_vibration_solver/ (use --plot).
"""
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
    from scipy.sparse.linalg import eigsh, spsolve
    _HAS_SCIPY = True
except ImportError:
    _HAS_SCIPY = False

FIXTURE_ROOT = Path(__file__).resolve().parent.parent / "tSiNCs/fixtures/vibration_parallel"
ADAMANTANE_NPZ = FIXTURE_ROOT / "hessian_mmff/adamantane_mmff_dense.npz"
ADAMANTANE_MODES = FIXTURE_ROOT / "expected/adamantane_omegas_modes.txt"
DIAMOND_BLOCKS = FIXTURE_ROOT / "hessian_mmff/diamond_nc_R6_sparse_blocks.npz"
DIAMOND_XYZ = FIXTURE_ROOT / "structures/diamond_nc_R6_relaxed.xyz"
DEFAULT_OUTDIR = Path(__file__).resolve().parent / "OUT_vibration_solver"


def _require_fixtures():
    missing = [p for p in (ADAMANTANE_NPZ, ADAMANTANE_MODES) if not p.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing vibration_parallel fixtures: "
            + ", ".join(str(p) for p in missing)
            + "\nRun: python3 tests/tSiNCs/bootstrap_vibration_parallel_fixtures.py"
        )


def load_adamantane_fixture():
    _require_fixtures()
    d = np.load(ADAMANTANE_NPZ, allow_pickle=True)
    ref_modes = np.loadtxt(ADAMANTANE_MODES)[:, 1] if ADAMANTANE_MODES.is_file() else d["omegas_modes"]
    return {
        "H": d["H"],
        "H_projected": d["H_projected"],
        "M": d["M"],
        "pos": d["pos"],
        "omegas_modes_ref": ref_modes,
        "natoms": int(d["natoms"]),
    }


def mass_weighted_mode_freqs(H, M):
    n = H.shape[0]
    m_sqrt = np.sqrt(np.diag(M))
    m_inv_sqrt = 1.0 / m_sqrt
    K_mw = (m_inv_sqrt[:, None] * H) * m_inv_sqrt[None, :]
    w = np.linalg.eigvalsh(K_mw)
    return np.sqrt(np.maximum(w, 0.0))


def symmetric_excitation(pos):
    charges = np.ones(pos.shape[0])
    com = np.mean(pos, axis=0)
    dir_x = np.zeros((pos.shape[0], 3))
    dir_x[:, 0] = pos[:, 0] - com[0]
    return dir_x, charges


def compare_mode_lists(modes_a, modes_b, n_compare=20, skip_rigid=6, tol=0.01):
    a = np.sort(modes_a)[skip_rigid:skip_rigid + n_compare]
    b = np.sort(modes_b)[skip_rigid:skip_rigid + n_compare]
    n = min(len(a), len(b))
    return np.max(np.abs(a[:n] - b[:n]))


def _draw_mode_lines(ax, omegas_modes, fmin, fmax, skip=6, color='gray', alpha=0.35):
    modes = np.sort(omegas_modes)[skip:]
    for m in modes:
        if fmin < m < fmax:
            ax.axvline(m, color=color, ls='--', lw=0.7, alpha=alpha)


def plot_spectrum_panel(outdir, omegas, spectra, labels, omegas_modes, title, fname, ylabel='|u|^2 response'):
    fig, ax = plt.subplots(figsize=(10, 5))
    for y, lab in zip(spectra, labels):
        ax.plot(omegas, y, lw=1.2, label=lab)
    _draw_mode_lines(ax, omegas_modes, omegas[0], omegas[-1])
    ax.set_xlabel(r'$\omega$')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    path = outdir / fname
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"  saved plot: {path}")
    return path


def plot_spectra_adamantane(args, fx=None):
    """Main spectrum comparison: sparse vs dense, with/without rigid projection."""
    if not _HAS_SCIPY:
        raise ImportError("scipy required for spectrum plots")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if fx is None:
        fx = load_adamantane_fixture()
    H, M, pos = fx["H_projected"], fx["M"], fx["pos"]
    H_raw = fx["H"]
    omegas_modes = mass_weighted_mode_freqs(H, M)
    omegas = np.linspace(args.fmin, args.fmax, args.nfreq)
    dir_x, charges = symmetric_excitation(pos)

    res_buggy = FTIR.mechanical_greens_probing_sparse(
        sp.csr_matrix(H_raw), M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges)
    res_sparse = FTIR.mechanical_greens_probing_sparse(
        sp.csr_matrix(H_raw), M, omegas, eta=args.eta, pos=pos, shift_rigid=args.shift,
        direction_vec=dir_x, charges=charges)
    res_dense = FTIR.mechanical_greens_probing(H, M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges)
    res_modes = FTIR.vibration_spectrum_from_modes(H, M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges)

    plot_spectrum_panel(
        outdir, omegas,
        [res_buggy["energy"], res_sparse["energy"], res_dense["energy"]],
        ["sparse (no rigid proj.)", "sparse + rigid proj.", "dense probe (ref.)"],
        omegas_modes,
        f"Adamantane vibration spectrum (eta={args.eta}) — sparse fix",
        "spectrum_sparse_vs_dense.png",
    )
    plot_spectrum_panel(
        outdir, omegas,
        [res_dense["energy"], res_modes["energy"]],
        ["dense per-omega solve", "mode-sum (eigh once)"],
        omegas_modes,
        "Dense probe vs mode-sum (note: different damping formulas)",
        "spectrum_dense_vs_modesum.png",
    )
    # zoom low-frequency
    mask = omegas < 1.5
    plot_spectrum_panel(
        outdir, omegas[mask],
        [res_buggy["energy"][mask], res_sparse["energy"][mask], res_dense["energy"][mask]],
        ["sparse no proj.", "sparse + proj.", "dense"],
        omegas_modes,
        "Low-frequency zoom (rigid-mode artifact visible without projection)",
        "spectrum_lowfreq_zoom.png",
    )
    # mode frequencies bar
    fig, ax = plt.subplots(figsize=(10, 4))
    modes_phys = np.sort(omegas_modes)[6:30]
    ax.bar(np.arange(len(modes_phys)), modes_phys, color='steelblue', alpha=0.8)
    ax.set_xlabel('mode index (physical)')
    ax.set_ylabel(r'$\omega_i$')
    ax.set_title('Adamantane projected-Hessian mode frequencies (eigh)')
    fig.tight_layout()
    p_modes = outdir / "mode_frequencies.png"
    fig.savefig(p_modes, dpi=160)
    plt.close(fig)
    print(f"  saved plot: {p_modes}")
    return outdir


def plot_jacobi_diagnostics(args, fx=None):
    """Diagnostic: block-dominance ratio and inertial-shift Jacobi convergence vs mu."""
    if not _HAS_SCIPY:
        return
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    if fx is None:
        fx = load_adamantane_fixture()
    H_raw, M, pos = fx["H"], fx["M"], fx["pos"]
    H = FTIR.apply_rigid_mode_shift(sp.csr_matrix(H_raw), M, pos, shift=args.shift)
    n_atoms = pos.shape[0]
    om_modes = mass_weighted_mode_freqs(H.toarray() if sp.issparse(H) else H, M)
    omega = float(om_modes[6]) * 1.15
    A = FTIR._dynamic_stiffness_sparse(H, M, omega, eta=args.eta, stabilize=args.stabilize)
    f = np.zeros(n_atoms * 3, dtype=np.complex128)
    f[0] = 1.0
    u_ref = spsolve(A, f)

    dom_A = FTIR.block_diag_dominance_ratio(A, n_atoms)
    dom_K = FTIR.block_diag_dominance_ratio(sp.csr_matrix(H_raw), n_atoms)
    print(f"  block dominance ratio: K={dom_K:.4f}  A(omega)={dom_A:.4f}  (PD needs >>1)")

    mus = np.logspace(0, 8, 13)
    errs_shifted, errs_true, iters, doms = [], [], [], []
    for mu in mus:
        As = A + FTIR.inertial_shift_matrix(M, n_atoms, mu)
        doms.append(FTIR.block_diag_dominance_ratio(As, n_atoms))
        u_true_ref = u_ref
        u_mu_ref = spsolve(As, f)
        try:
            u, nit, _ = FTIR.solve_block_jacobi_inertial(
                A, M, f, mu=mu, max_iter=2000, tol=1e-7, verbose=False)
            errs_shifted.append(np.linalg.norm(u - u_mu_ref) / max(np.linalg.norm(u_mu_ref), 1e-30))
            errs_true.append(np.linalg.norm(u - u_true_ref) / max(np.linalg.norm(u_true_ref), 1e-30))
            iters.append(nit)
            print(f"  mu={mu:8.1e} dom={doms[-1]:6.2f}  Jacobi iters={nit:4d}  "
                  f"err(shifted)={errs_shifted[-1]:.2e}  err(true A)={errs_true[-1]:.2e}")
        except Exception as e:
            errs_shifted.append(np.nan)
            errs_true.append(np.nan)
            iters.append(np.nan)
            print(f"  mu={mu:8.1e} dom={doms[-1]:6.2f}  FAILED ({e})")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4))
    ax1, ax2, ax3 = axes
    ax1.loglog(mus, doms, 'k-o', lw=1.5)
    ax1.axhline(1.0, color='r', ls='--', alpha=0.6, label='dominance=1')
    ax1.set_xlabel(r'inertial shift $\mu$')
    ax1.set_ylabel('block diag. dominance ratio')
    ax1.set_title('PD-like $A + \\mu M$ conditioning')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax2.loglog(mus, errs_shifted, 'b-o', label='vs spsolve(A+μM)')
    ax2.loglog(mus, errs_true, 'r-o', label='vs spsolve(A) true')
    ax2.set_xlabel(r'$\mu$')
    ax2.set_ylabel('relative error after Jacobi')
    ax2.set_title('Block Jacobi on shifted system')
    ax2.legend(fontsize=8)
    ax2.grid(True, alpha=0.3)
    ax3.semilogx(mus, iters, 'coral', marker='o', lw=1.5)
    ax3.set_xlabel(r'$\mu$')
    ax3.set_ylabel('iterations')
    ax3.set_title('Jacobi iteration count')
    ax3.grid(True, alpha=0.3)
    fig.suptitle(f'Adamantane A(ω) at ω={omega:.3f}: K dominance={dom_K:.3f}, A dominance={dom_A:.3f}', fontsize=10)
    fig.tight_layout()
    p = outdir / "jacobi_inertial_shift_diagnostic.png"
    fig.savefig(p, dpi=160)
    plt.close(fig)
    print(f"  saved plot: {p}")


def run_stage0(args):
    print("=== Stage 0: dense eigh reference (adamantane) ===")
    fx = load_adamantane_fixture()
    H, M = fx["H_projected"], fx["M"]
    omegas_modes = mass_weighted_mode_freqs(H, M)
    mode_diff = compare_mode_lists(omegas_modes, fx["omegas_modes_ref"], tol=args.mode_tol)
    print(f"  mode freq max diff vs fixture: {mode_diff:.6e} (tol={args.mode_tol})")
    if mode_diff > args.mode_tol:
        raise AssertionError(f"M-S0 failed: mode frequencies differ by {mode_diff:.6e} > {args.mode_tol}")
    if args.plot:
        plot_spectra_adamantane(args, fx)
    print("M-S0 PASSED")
    return {"omegas_modes": omegas_modes}


def run_stage1(args):
    if not _HAS_SCIPY:
        raise ImportError("scipy required for stage 1")
    print("=== Stage 1: scipy sparse probe with rigid projection (adamantane) ===")
    fx = load_adamantane_fixture()
    H, M, pos = fx["H_projected"], fx["M"], fx["pos"]
    H_raw = fx["H"]
    omegas = np.linspace(args.fmin, args.fmax, args.nfreq)
    dir_x, charges = symmetric_excitation(pos)

    res_buggy = FTIR.mechanical_greens_probing_sparse(
        sp.csr_matrix(H_raw), M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges)
    res_sparse = FTIR.mechanical_greens_probing_sparse(
        sp.csr_matrix(H_raw), M, omegas, eta=args.eta, pos=pos, shift_rigid=args.shift,
        fail_loud=args.fail_loud, direction_vec=dir_x, charges=charges)
    res_dense = FTIR.mechanical_greens_probing(H, M, omegas, eta=args.eta, direction_vec=dir_x, charges=charges)
    rel_env = np.max(np.abs(res_sparse["energy"] - res_dense["energy"])) / (np.max(res_dense["energy"]) + 1e-12)
    print(f"  buggy first peak omega ~ {res_buggy['omega'][np.argmax(res_buggy['energy'])]:.4f}")
    print(f"  projected sparse vs dense probe: max_rel_env={rel_env:.6e} (tol={args.env_tol})")
    if rel_env > args.env_tol:
        raise AssertionError(f"M-S1 failed: rel_env={rel_env:.6e} > {args.env_tol}")

    m_vec = np.diag(M)
    k = min(20, H.shape[0] - 7)
    w_eigsh, _ = eigsh(H, k=k, M=sp.diags(m_vec), which="SM")
    om_eigsh = np.sqrt(np.maximum(w_eigsh, 0.0))
    om_dense = np.sort(mass_weighted_mode_freqs(H, M))
    eigsh_diff = np.max(np.abs(np.sort(om_eigsh) - om_dense[:k]))
    print(f"  eigsh vs dense eigh: max diff={eigsh_diff:.6e}")
    if eigsh_diff > args.mode_tol:
        raise AssertionError(f"M-S1 eigsh check failed: {eigsh_diff:.6e}")

    if args.plot:
        plot_spectra_adamantane(args, fx)
        plot_jacobi_diagnostics(args, fx)

    print("M-S1 PASSED")


def run_stage2(args):
    if not _HAS_SCIPY:
        raise ImportError("scipy required for stage 2")
    print("=== Stage 2: CPU point-Jacobi vs spsolve (scalar chain) ===")
    n = 15
    K = sp.diags([-np.ones(n - 1), 2.0 * np.ones(n), -np.ones(n - 1)], offsets=[-1, 0, 1], format='csr', dtype=np.float64)
    f = np.ones(n, dtype=np.float64)
    u_direct = spsolve(K, f)
    u_iter, n_iter, rel_res = FTIR.solve_linear_jacobi_momentum(
        K, f, max_iter=args.jacobi_max_iter, tol=args.jacobi_tol, method='point', verbose=args.verbose)
    rel_err = np.linalg.norm(u_iter - u_direct) / max(np.linalg.norm(u_direct), 1e-12)
    print(f"  scalar chain: Jacobi iters={n_iter} rel_res={rel_res:.3e} rel_err vs spsolve={rel_err:.3e}")
    if rel_err > args.solve_tol:
        raise AssertionError(f"M-S2 failed: Jacobi rel_err={rel_err:.3e} > {args.solve_tol}")
    if args.plot:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        res_hist = []
        u = np.zeros(n)
        A = K.tocsr()
        dinv = 1.0 / A.diagonal()
        for k in range(min(n_iter + 5, 800)):
            u = dinv * (f - A @ u + A.diagonal() * u)
            res_hist.append(np.linalg.norm(f - A @ u) / np.linalg.norm(f))
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.semilogy(res_hist, 'b-', lw=1.2)
        ax.axhline(args.jacobi_tol, color='r', ls='--', label=f'tol={args.jacobi_tol}')
        ax.set_xlabel('Jacobi iteration')
        ax.set_ylabel('relative residual')
        ax.set_title('Scalar chain: point-Jacobi convergence')
        ax.legend()
        ax.grid(True, alpha=0.3)
        fig.tight_layout()
        p = outdir / "jacobi_scalar_convergence.png"
        fig.savefig(p, dpi=160)
        plt.close(fig)
        print(f"  saved plot: {p}")
    print("M-S2 PASSED")


def main():
    ap = argparse.ArgumentParser(description="Vibration solver validation ladder (M-S0..S2) + plots")
    ap.add_argument("--stage", type=int, default=None, choices=[0, 1, 2], help="Validation stage (omit with --plot-only)")
    ap.add_argument("--plot", action="store_true", help="Save spectrum/diagnostic PNGs")
    ap.add_argument("--plot-only", action="store_true", help="Only generate plots (no stage asserts)")
    ap.add_argument("--outdir", type=str, default=str(DEFAULT_OUTDIR))
    ap.add_argument("--fmin", type=float, default=0.01)
    ap.add_argument("--fmax", type=float, default=3.5)
    ap.add_argument("--nfreq", type=int, default=800)
    ap.add_argument("--eta", type=float, default=0.05)
    ap.add_argument("--shift", type=float, default=1e6)
    ap.add_argument("--mode-tol", type=float, default=0.01, dest="mode_tol")
    ap.add_argument("--env-tol", type=float, default=0.02, dest="env_tol")
    ap.add_argument("--solve-tol", type=float, default=1e-4, dest="solve_tol")
    ap.add_argument("--jacobi-max-iter", type=int, default=2000)
    ap.add_argument("--jacobi-tol", type=float, default=1e-8)
    ap.add_argument("--stabilize", type=float, default=1e-6)
    ap.add_argument("--fail-loud", action="store_true")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    if args.plot_only:
        args.plot = True
        print("=== Generating vibration solver plots (adamantane) ===")
        plot_spectra_adamantane(args)
        plot_jacobi_diagnostics(args)
        print(f"Plots in: {args.outdir}")
        return

    if args.stage is None:
        ap.error("specify --stage 0|1|2, or use --plot-only")

    if args.stage == 0:
        run_stage0(args)
    elif args.stage == 1:
        run_stage1(args)
    elif args.stage == 2:
        run_stage2(args)


if __name__ == "__main__":
    main()
