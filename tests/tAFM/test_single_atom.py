#!/usr/bin/env python3
"""
test_single_atom.py — Single-atom density projection test harness.
Thin top-level script calling reusable functions from pyBall.OCL.AFM.
Projects density for H and/or C atoms, measures integrated electron count and radial profile.
"""
import sys, os, argparse
import numpy as np

_THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT      = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL.AFM import project_single_atom, Z_TO_ZVAL, ELEM_Z

_FDATA_HCNO  = os.path.join(_ROOT, 'tests', 'Fireball', 'Fdata_HCNO')
_FDATA_BASIS = os.path.join(_FDATA_HCNO, 'basis')

ap = argparse.ArgumentParser()
ap.add_argument('--atom', choices=['H','C','both'], default='both')
ap.add_argument('--step', type=float, default=0.10, help='Grid step (A)')
ap.add_argument('--margin', type=float, default=6.0, help='Box half-size (A)')
ap.add_argument('--noplot', action='store_true')
ap.add_argument('--notiled', action='store_true')
args = ap.parse_args()

# ── Helper: radial profile from 3D grid ──
def radial_profile_grid(rho, step, n, rmax=5.0, nbins=200):
    """Spherical average of rho(r) around grid center."""
    c = n // 2
    R = int(np.ceil(rmax / step))
    x0, x1 = max(0, c-R), min(n, c+R+1)
    xs = (np.arange(x0, x1) - c) * step
    X, Y, Z = np.meshgrid(xs, xs, xs, indexing='ij')
    d = np.sqrt(X**2 + Y**2 + Z**2)
    v = rho[x0:x1, x0:x1, x0:x1].astype(np.float64)
    d = d.ravel(); v = v.ravel()
    m = d <= rmax
    d = d[m]; v = v[m]
    bins = np.linspace(0, rmax, nbins+1)
    ib = np.clip(np.searchsorted(bins, d, side='right')-1, 0, nbins-1)
    sum_v = np.bincount(ib, weights=v, minlength=nbins)
    cnt_v = np.bincount(ib, minlength=nbins)
    avg = sum_v / np.maximum(cnt_v, 1)
    rc = 0.5*(bins[:-1] + bins[1:])
    return rc, avg

# ── Main ──
atoms_to_test = []
if args.atom in ('H', 'both'): atoms_to_test.append('H')
if args.atom in ('C', 'both'): atoms_to_test.append('C')

for elem in atoms_to_test:
    Z = ELEM_Z[elem]
    nval = Z_TO_ZVAL[Z]
    print(f"\n{'='*60}")
    print(f"=== Single-atom test: {elem}  Z={Z}  N_val={nval} ===")
    print(f"{'='*60}")

    # Build diagonal density matrix with valence occupations
    rho_4x4 = np.zeros((4, 4), dtype=np.float32)
    if Z == 1:
        rho_4x4[0, 0] = 1.0
    elif Z == 6:
        rho_4x4[0, 0] = 2.0
        rho_4x4[1, 1] = 2.0 / 3.0
        rho_4x4[2, 2] = 2.0 / 3.0
        rho_4x4[3, 3] = 2.0 / 3.0

    print(f"\n--- Test A: hand-built diagonal rho (occ = valence) ---")
    rho_A, gs_A, int_A, projector = project_single_atom(
        Z, rho_4x4, args.step, args.margin, _FDATA_BASIS,
        use_tiled=(not args.notiled))
    print(f"  Expected {nval} e, got {int_A:.4f} e  (ratio = {int_A/nval:.4f})")

    # Compute radial overlap integrals from basis data
    wfs = projector.basis_data[Z]
    r_s = np.arange(wfs[0]['mesh']) * (wfs[0]['rcutoff'] / (wfs[0]['mesh'] - 1))
    R_s = wfs[0]['data']
    S_ss = np.trapz(R_s**2 * r_s**2, r_s)
    print(f"  Radial overlap S_ss = {S_ss:.6f}")
    if len(wfs) > 1:
        r_p = np.arange(wfs[1]['mesh']) * (wfs[1]['rcutoff'] / (wfs[1]['mesh'] - 1))
        R_p = wfs[1]['data']
        S_pp = np.trapz(R_p**2 * r_p**2, r_p)
        print(f"  Radial overlap S_pp = {S_pp:.6f}")
    else:
        S_pp = S_ss

    # Test B: rescaled density matrix
    print(f"\n--- Test B: rescaled rho (occ / S_ii) ---")
    rho_4x4_B = np.zeros((4, 4), dtype=np.float32)
    if Z == 1:
        rho_4x4_B[0, 0] = 1.0 / S_ss
    elif Z == 6:
        rho_4x4_B[0, 0] = 2.0 / S_ss
        rho_4x4_B[1, 1] = (2.0/3.0) / S_pp
        rho_4x4_B[2, 2] = (2.0/3.0) / S_pp
        rho_4x4_B[3, 3] = (2.0/3.0) / S_pp
    rho_B, gs_B, int_B, _ = project_single_atom(
        Z, rho_4x4_B, args.step, args.margin, _FDATA_BASIS,
        use_tiled=(not args.notiled))
    print(f"  Expected {nval} e, got {int_B:.4f} e  (ratio = {int_B/nval:.4f})")

    # Plots
    if not args.noplot:
        import matplotlib; matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        n = gs_A['ngrid'][0]
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle(f"Single-atom test: {elem} (Z={Z}, Nval={nval})", fontsize=11)

        rc_A, avg_A = radial_profile_grid(rho_A, args.step, n, rmax=5.0)
        rc_B, avg_B = radial_profile_grid(rho_B, args.step, n, rmax=5.0)
        axes[0].plot(rc_A, avg_A, label=f'Test A (int={int_A:.3f})', lw=1.2)
        axes[0].plot(rc_B, avg_B, label=f'Test B (int={int_B:.3f})', lw=1.2, ls='--')
        axes[0].set_xlabel('r (A)'); axes[0].set_ylabel('rho (e/A^3)')
        axes[0].set_title('Radial profile (linear)')
        axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.2)

        axes[1].semilogy(rc_A, np.maximum(avg_A, 1e-12), label='Test A', lw=1.2)
        axes[1].semilogy(rc_B, np.maximum(avg_B, 1e-12), label='Test B', lw=1.2, ls='--')
        axes[1].set_xlabel('r (A)'); axes[1].set_ylabel('rho (e/A^3)')
        axes[1].set_title('Radial profile (log)')
        axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.2)

        c = n // 2
        sl = rho_A[:, :, c].T
        im = axes[2].imshow(sl, origin='lower', cmap='magma', aspect='equal',
                            extent=[-args.margin, args.margin, -args.margin, args.margin])
        axes[2].set_title(f'XY slice (z=0) Test A  max={sl.max():.4f}')
        axes[2].set_xlabel('x (A)'); axes[2].set_ylabel('y (A)')
        plt.colorbar(im, ax=axes[2], shrink=0.8)

        plt.tight_layout()
        fn = f'single_atom_{elem}.png'
        plt.savefig(os.path.join(_THIS_DIR, fn), dpi=130, bbox_inches='tight')
        plt.close()
        print(f"  Saved {fn}")

print("\n=== Single-atom test complete ===")
