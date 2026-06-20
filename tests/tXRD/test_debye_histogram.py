#!/usr/bin/env python3
"""Test XRDDebye histogram + Debye transform on existing nanocrystal fixtures.

Run from repo root or tests/tXRD/:
    python tests/tXRD/test_debye_histogram.py

Expects fixtures at tests/tSiNCs/fixtures/vibration_parallel/structures/.
Outputs plots to OUT_XRD/.
"""
import os
import sys
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, REPO)

from pyBall.XRD import (
    XRDDebye, load_xyz_atoms, get_form_factor_table, compute_sigma_from_sparse_blocks
)


def compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist, bin_edges):
    """Helper: Debye transform + self-scattering term."""
    ff_table = get_form_factor_table(unique_Z, Q)
    I = engine.compute_spectrum(hist, bin_edges, Q, ff_table, unique_Z)
    n_species = len(unique_Z)
    elem_idx = atoms[:, 3].astype(np.int32)
    counts = np.bincount(elem_idx, minlength=n_species).astype(np.float32)
    I_self = np.sum(counts[:, None] * (ff_table ** 2), axis=0)
    return I + I_self


def test_diamond_nc():
    xyz_path = os.path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/structures/diamond_nc_R6_relaxed.xyz')
    hess_path = os.path.join(REPO, 'tests/tSiNCs/fixtures/vibration_parallel/hessian_mmff/diamond_nc_R6_sparse_blocks.npz')
    if not os.path.exists(xyz_path):
        raise FileNotFoundError(f"Fixture not found: {xyz_path}")

    print(f"Loading {xyz_path} ...")
    atoms, unique_Z = load_xyz_atoms(xyz_path)
    print(f"  natoms={atoms.shape[0]}  species={unique_Z}")

    engine = XRDDebye(preferred_vendor='nvidia')

    Q = np.linspace(0.5, 15.0, 800, dtype=np.float32)
    r_max = 20.0
    dr = 0.02

    # ── 1. Static histogram (no thermal broadening) ──
    print(f"\n[1] Static histogram (r_max={r_max}, dr={dr}) ...")
    hist_static, bin_edges = engine.compute_histogram(atoms, r_max=r_max, dr=dr)
    print(f"    hist shape={hist_static.shape}  total_pairs={hist_static.sum():.1f}")
    I_static = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_static, bin_edges)

    # ── 2. Constant thermal broadening ──
    sigma_const = 0.015  # Å, typical room-temperature C-C thermal width
    print(f"\n[2] Constant Gaussian broadening (sigma={sigma_const} Å) ...")
    hist_gconst, _ = engine.compute_histogram_gaussian(atoms, sigma_const, r_max=r_max, dr=dr)
    print(f"    total_pairs={hist_gconst.sum():.1f}")
    I_gconst = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_gconst, bin_edges)

    # ── 3. Hessian-based thermal broadening ──
    I_ghess = None
    if os.path.exists(hess_path):
        print(f"\n[3] Hessian-based Gaussian broadening ...")
        hess = np.load(hess_path)
        pos = atoms[:, :3].astype(np.float64)
        sigma_hess = compute_sigma_from_sparse_blocks(
            pos, hess['neigh_idx'], hess['neigh_count'], hess['blocks'],
            kBT=0.02585, default_sigma=0.02
        )
        # Report sigma statistics
        n_pairs = atoms.shape[0] * (atoms.shape[0] - 1) // 2
        triu = np.triu_indices(atoms.shape[0], k=1)
        sig_vals = sigma_hess[triu]
        print(f"    sigma range: {sig_vals.min():.5f} .. {sig_vals.max():.5f}  median={np.median(sig_vals):.5f}")
        hist_ghess, _ = engine.compute_histogram_gaussian(atoms, sigma_hess, r_max=r_max, dr=dr)
        print(f"    total_pairs={hist_ghess.sum():.1f}")
        I_ghess = compute_spectrum_with_self(engine, atoms, unique_Z, Q, hist_ghess, bin_edges)
    else:
        print(f"\n[3] Hessian file not found, skipping Hessian broadening: {hess_path}")

    # ── Save & Plot ──
    outdir = os.path.join(REPO, 'OUT_XRD')
    os.makedirs(outdir, exist_ok=True)

    npz_path = os.path.join(outdir, 'diamond_nc_R6_xrd.npz')
    save_dict = dict(Q=Q, I_static=I_static, I_gconst=I_gconst,
                     hist_static=hist_static, bin_edges=bin_edges, unique_Z=unique_Z)
    if I_ghess is not None:
        save_dict['I_ghess'] = I_ghess
    np.savez(npz_path, **save_dict)
    print(f"\nSaved NPZ: {npz_path}")

    try:
        import matplotlib.pyplot as plt

        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(3, 3)

        # Row 0: histograms
        ax = fig.add_subplot(gs[0, 0])
        ax.bar(bin_edges[:-1], hist_static[0], width=np.diff(bin_edges), align='edge', color='C0', alpha=0.5, label='static')
        ax.set_xlabel('r (Å)')
        ax.set_ylabel('pair count')
        ax.set_title('Histogram (C-C) — static')

        ax = fig.add_subplot(gs[0, 1])
        ax.bar(bin_edges[:-1], hist_gconst[0], width=np.diff(bin_edges), align='edge', color='C1', alpha=0.5, label='const σ')
        ax.set_xlabel('r (Å)')
        ax.set_ylabel('pair count')
        ax.set_title(f'Histogram — const σ={sigma_const} Å')

        if I_ghess is not None:
            ax = fig.add_subplot(gs[0, 2])
            ax.bar(bin_edges[:-1], hist_ghess[0], width=np.diff(bin_edges), align='edge', color='C2', alpha=0.5, label='Hessian σ')
            ax.set_xlabel('r (Å)')
            ax.set_ylabel('pair count')
            ax.set_title('Histogram — Hessian σ')

        # Row 1: I(Q) linear scale
        ax = fig.add_subplot(gs[1, :])
        ax.plot(Q, I_static, '-', lw=1.0, color='C0', label='static')
        ax.plot(Q, I_gconst, '-', lw=1.0, color='C1', label=f'const σ={sigma_const} Å')
        if I_ghess is not None:
            ax.plot(Q, I_ghess, '-', lw=1.0, color='C2', label='Hessian σ')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('I(Q)')
        ax.set_title('Powder XRD — diamond NC R=6Å (270 atoms)')
        ax.legend()

        # Row 2: I(Q) log scale + difference
        ax = fig.add_subplot(gs[2, 0])
        ax.semilogy(Q, I_static, '-', lw=1.0, color='C0', label='static')
        ax.semilogy(Q, I_gconst, '-', lw=1.0, color='C1', label='const σ')
        if I_ghess is not None:
            ax.semilogy(Q, I_ghess, '-', lw=1.0, color='C2', label='Hessian σ')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('log I(Q)')
        ax.set_title('Log scale')
        ax.legend()

        ax = fig.add_subplot(gs[2, 1])
        ax.plot(Q, I_gconst - I_static, '-', lw=1.0, color='C1', label='const σ − static')
        if I_ghess is not None:
            ax.plot(Q, I_ghess - I_static, '-', lw=1.0, color='C2', label='Hessian σ − static')
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('ΔI(Q)')
        ax.set_title('Difference vs static')
        ax.legend()

        ax = fig.add_subplot(gs[2, 2])
        if I_ghess is not None:
            ax.plot(Q, I_ghess - I_gconst, '-', lw=1.0, color='C3', label='Hessian − const')
            ax.set_title('Hessian vs const σ')
        else:
            ax.text(0.5, 0.5, 'Hessian data\nnot available', ha='center', va='center', transform=ax.transAxes)
        ax.set_xlabel('Q (Å⁻¹)')
        ax.set_ylabel('ΔI(Q)')
        ax.legend()

        fig.tight_layout()
        png_path = os.path.join(outdir, 'diamond_nc_R6_xrd.png')
        fig.savefig(png_path, dpi=200)
        print(f"Saved plot: {png_path}")
    except Exception as e:
        print(f"Plotting skipped: {e}")

    return I_static, Q


if __name__ == '__main__':
    test_diamond_nc()
