#!/usr/bin/env python3
"""
test_single_atom.py — Single-atom density projection test harness.
Projects density for a single H atom and a single C atom using:
  (a) Fireball SCF density matrix (from a real SCF calculation on the isolated atom)
  (b) Hand-built diagonal density matrix with known occupations
Measures integrated electron count and radial profile to verify normalization.
"""
import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR  = os.path.dirname(os.path.abspath(__file__))
_ROOT      = os.path.realpath(os.path.join(_THIS_DIR, '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

_FDATA_HCNO  = os.path.join(_ROOT, 'tests', 'Fireball', 'Fdata_HCNO')
_FDATA_BASIS = os.path.join(_FDATA_HCNO, 'basis')

ap = argparse.ArgumentParser()
ap.add_argument('--atom', choices=['H','C','both'], default='both')
ap.add_argument('--step', type=float, default=0.10, help='Grid step (Å)')
ap.add_argument('--margin', type=float, default=6.0, help='Box half-size (Å)')
ap.add_argument('--noplot', action='store_true')
ap.add_argument('--notiled', action='store_true')
args = ap.parse_args()

ELEM_Z  = {'H': 1, 'C': 6}
Z_NVAL  = {1: 1, 6: 4}
RCUT    = {1: 2.3, 6: 2.6}  # .wf rcutoff is in Bohr; abohr*rc + margin
BLOCK   = 8

# ── Helper: build grid spec centered at origin ──
def make_grid(margin, step):
    half = margin
    n = int(np.ceil(2*half / step / BLOCK)) * BLOCK
    origin = np.array([-n*step/2, -n*step/2, -n*step/2], dtype=np.float32)
    return {
        'origin': origin,
        'dA': [step, 0., 0.], 'dB': [0., step, 0.], 'dC': [0., 0., step],
        'ngrid': np.array([n, n, n], dtype=np.int32),
    }, n

# ── Helper: project with a given rho matrix ──
def project_one_atom(projector, Z, rho_4x4, step, margin, use_tiled):
    """Project density for a single atom at origin with given 4x4 density matrix block."""
    grid_spec, n = make_grid(margin, step)
    atomPos   = np.array([[0.0, 0.0, 0.0]], dtype=np.float64)
    atomTypes = np.array([Z], dtype=np.int32)
    atoms_dict = {
        'pos':  atomPos,
        'Rcut': np.array([RCUT[Z]], dtype=np.float64),
        'type': atomTypes,
    }
    # Build a minimal sparse rho: shape (natoms=1, neigh_max, numorb_max, numorb_max)
    # We need neigh_max and numorb_max from the projector's conventions.
    # For a single isolated atom, the only neighbor is itself.
    neigh_max   = 1  # minimal
    numorb_max  = 4
    rho_sparse  = np.zeros((1, neigh_max, numorb_max, numorb_max), dtype=np.float32)
    rho_sparse[0, 0, :, :] = rho_4x4
    # neigh_j: atom 0's neighbor #0 is atom 0+1=1 (1-based Fortran convention)
    neigh_j = np.zeros((1, neigh_max), dtype=np.int32)
    neigh_j[0, 0] = 1  # self-neighbor (1-based)

    class FakeNeighs:
        pass
    neighs = FakeNeighs()
    neighs.rho     = rho_sparse
    neighs.neigh_j = neigh_j.ravel()
    neighs.neigh_max  = neigh_max
    neighs.numorb_max = numorb_max

    print(f"  [project_one_atom] Z={Z} rho_diag={np.diag(rho_4x4)} grid={n}^3 step={step}")
    rho_grid = projector.project(rho_sparse, neighs, atoms_dict, grid_spec,
                                  nMaxAtom=64, use_tiled=use_tiled)
    dV = step**3
    integral = float(rho_grid.sum() * dV)
    print(f"  rho_grid shape={rho_grid.shape}  range=[{rho_grid.min():.6f}, {rho_grid.max():.6f}]")
    print(f"  Integrated electrons = {integral:.6f}")
    return rho_grid, grid_spec, integral

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
from pyBall.FireballOCL import Grid as ocl_grid

projector = ocl_grid.GridProjector(fdata_dir=_FDATA_BASIS)

atoms_to_test = []
if args.atom in ('H', 'both'):
    atoms_to_test.append('H')
if args.atom in ('C', 'both'):
    atoms_to_test.append('C')

for elem in atoms_to_test:
    Z = ELEM_Z[elem]
    nval = Z_NVAL[Z]
    print(f"\n{'='*60}")
    print(f"=== Single-atom test: {elem}  Z={Z}  N_val={nval} ===")
    print(f"{'='*60}")

    projector.load_basis([Z])

    # Build diagonal density matrix with valence occupations
    # For H (s only): diag = [n_s, 0, 0, 0]  with n_s = 1
    # For C (s + 3p): diag = [n_s, n_py, n_pz, n_px]  with n_s=2, n_p=2/3 each
    # Kernel convention: (s, py, pz, px) — matches Fireball Ortega convention
    rho_4x4 = np.zeros((4, 4), dtype=np.float32)
    if Z == 1:
        rho_4x4[0, 0] = 1.0   # 1 electron in s
    elif Z == 6:
        rho_4x4[0, 0] = 2.0         # 2 electrons in s
        rho_4x4[1, 1] = 2.0 / 3.0   # 2/3 in py
        rho_4x4[2, 2] = 2.0 / 3.0   # 2/3 in pz
        rho_4x4[3, 3] = 2.0 / 3.0   # 2/3 in px

    print(f"\n--- Test A: hand-built diagonal rho (occ = valence) ---")
    rho_A, gs_A, int_A = project_one_atom(projector, Z, rho_4x4, args.step, args.margin, not args.notiled)
    print(f"  Expected {nval} e, got {int_A:.4f} e  (ratio = {int_A/nval:.4f})")

    # If the integral is off, compute the overlap integral S_ss from the radial table
    # and try the rescaled density matrix rho_ij = n_i / S_ii
    bm = projector.basis_meta
    nz_map = bm['nz_map']
    sp_idx = nz_map[Z]
    n_nodes = bm['n_nodes']
    dr = bm['dr']
    # Read radial data for shell 0 (s)
    # The packed basis is [n_species][max_shells][n_nodes]
    # We stored it via load_basis — re-read from .wf files
    wfs = projector.basis_data[Z]
    r_s = np.arange(wfs[0]['mesh']) * (wfs[0]['rcutoff'] / (wfs[0]['mesh'] - 1))
    R_s = wfs[0]['data']
    S_ss = np.trapz(R_s**2 * r_s**2, r_s)  # ∫ R_s² r² dr (the Ylm integral is 1)
    print(f"  Radial overlap S_ss = ∫ R_s² r² dr = {S_ss:.6f}")
    if len(wfs) > 1:
        r_p = np.arange(wfs[1]['mesh']) * (wfs[1]['rcutoff'] / (wfs[1]['mesh'] - 1))
        R_p = wfs[1]['data']
        S_pp = np.trapz(R_p**2 * r_p**2, r_p)
        print(f"  Radial overlap S_pp = ∫ R_p² r² dr = {S_pp:.6f}")
    else:
        S_pp = S_ss

    # Test B: rescaled density matrix — ρ_ii = n_i / S_ii
    print(f"\n--- Test B: rescaled rho (occ / S_ii) ---")
    rho_4x4_B = np.zeros((4, 4), dtype=np.float32)
    if Z == 1:
        rho_4x4_B[0, 0] = 1.0 / S_ss
    elif Z == 6:
        rho_4x4_B[0, 0] = 2.0 / S_ss
        rho_4x4_B[1, 1] = (2.0/3.0) / S_pp
        rho_4x4_B[2, 2] = (2.0/3.0) / S_pp
        rho_4x4_B[3, 3] = (2.0/3.0) / S_pp
    rho_B, gs_B, int_B = project_one_atom(projector, Z, rho_4x4_B, args.step, args.margin, not args.notiled)
    print(f"  Expected {nval} e, got {int_B:.4f} e  (ratio = {int_B/nval:.4f})")

    # Radial profiles
    if not args.noplot:
        _, n = make_grid(args.margin, args.step)
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        fig.suptitle(f"Single-atom test: {elem} (Z={Z}, Nval={nval})", fontsize=11)

        rc_A, avg_A = radial_profile_grid(rho_A, args.step, n, rmax=5.0)
        rc_B, avg_B = radial_profile_grid(rho_B, args.step, n, rmax=5.0)
        axes[0].plot(rc_A, avg_A, label=f'Test A (int={int_A:.3f})', lw=1.2)
        axes[0].plot(rc_B, avg_B, label=f'Test B (int={int_B:.3f})', lw=1.2, ls='--')
        axes[0].set_xlabel('r (Å)'); axes[0].set_ylabel('⟨ρ⟩ (e/Å³)')
        axes[0].set_title('Radial profile (linear)')
        axes[0].legend(fontsize=8); axes[0].grid(True, alpha=0.2)

        axes[1].semilogy(rc_A, np.maximum(avg_A, 1e-12), label='Test A', lw=1.2)
        axes[1].semilogy(rc_B, np.maximum(avg_B, 1e-12), label='Test B', lw=1.2, ls='--')
        axes[1].set_xlabel('r (Å)'); axes[1].set_ylabel('⟨ρ⟩ (e/Å³)')
        axes[1].set_title('Radial profile (log)')
        axes[1].legend(fontsize=8); axes[1].grid(True, alpha=0.2)

        # 2D slice through center
        c = n // 2
        sl = rho_A[:, :, c].T
        im = axes[2].imshow(sl, origin='lower', cmap='magma', aspect='equal',
                            extent=[-args.margin, args.margin, -args.margin, args.margin])
        axes[2].set_title(f'XY slice (z=0) Test A  max={sl.max():.4f}')
        axes[2].set_xlabel('x (Å)'); axes[2].set_ylabel('y (Å)')
        plt.colorbar(im, ax=axes[2], shrink=0.8)

        plt.tight_layout()
        fn = f'single_atom_{elem}.png'
        plt.savefig(os.path.join(_THIS_DIR, fn), dpi=130, bbox_inches='tight')
        plt.close()
        print(f"  Saved {fn}")

print("\n=== Single-atom test complete ===")
