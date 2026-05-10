#!/usr/bin/env python3
"""
Test script for 2D Ewald electrostatic potential of ionic crystal surfaces.

Loads arbitrary substrate .xyz files (with lvs + charges in the 4th column)
via pyBall.AtomicSystem, delegates computation to pyBall.Ewald2D, and
produces 5 diagnostic figures.

Usage:
  python test_ewald_2d.py --xyz NaCl_1x1_L2.xyz [--n_harm 3] [--N_rep 30] [--grid 200]
  python test_ewald_2d.py --xyz CaF2_6L.xyz --n_harm 4

See:  doc/Topics/OnSurfaceAssembly/Ewald_2D.md  for the underlying derivation.
"""

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import argparse

sys.path.insert(0, os.path.expanduser("~/git/FireCore"))

from pyBall.AtomicSystem import AtomicSystem
from pyBall.Ewald2D import Ewald2D, generate_G_vectors, compute_w_per_ion, eval_potential_full_1d

# ============================================================
#  Plotting helpers  (thin wrappers around matplotlib)
# ============================================================

def plot_ions_xy(ax, rx, ry, q, a_vec, b_vec, Lx, Ly):
    """Scatter ions colored by charge sign, tiled 3×3 around cell."""
    colors = ['red' if qi > 0 else 'blue' for qi in q]
    sizes  = [80*max(abs(qi), 0.3) for qi in q]
    for dn in range(-1, 2):
        for dm in range(-1, 2):
            ax.scatter(rx + dn*a_vec[0] + dm*b_vec[0],
                       ry + dn*a_vec[1] + dm*b_vec[1],
                       c=colors, s=sizes, edgecolors='k', zorder=5)
    ax.set_xlim(0, Lx); ax.set_ylim(0, Ly)
    ax.set_aspect('equal')
    ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')


def plot_imshow_sym(ax, data, extent, title, cbar_label='value', cmap='RdBu_r'):
    """Symmetric imshow (vmin=−vmax) with colorbar."""
    vmax = np.max(np.abs(data))
    if vmax == 0: vmax = 1.0
    im = ax.imshow(data, extent=extent, origin='lower', cmap=cmap, vmin=-vmax, vmax=vmax)
    ax.set_aspect('equal'); ax.set_title(title)
    ax.set_xlabel('x (Å)'); ax.set_ylabel('y (Å)')
    plt.colorbar(im, ax=ax, label=cbar_label)
    return im


# ============================================================
#  MAIN
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="2D Ewald surface potential test")
    parser.add_argument("--xyz",    type=str, required=True,  help="Path to .xyz file (must have 'lvs' in comment and charges in 4th column)")
    parser.add_argument("--n_harm", type=int, default=3,      help="G-vector half-width |h|,|k| <= n_harm  (default: 3)")
    parser.add_argument("--N_rep",  type=int, default=30,     help="Number of periodic replicas for brute-force  (default: 30)")
    parser.add_argument("--grid",   type=int, default=200,    help="Grid resolution per dimension  (default: 200)")
    parser.add_argument("--noPlot", action="store_true",      help="Skip plt.show()")
    parser.add_argument("--prefix", type=str, default=None,   help="Output filename prefix (default: derived from xyz)")
    args = parser.parse_args()

    # --- Load system via AtomicSystem ---
    sys_at = AtomicSystem(fname=args.xyz)
    assert sys_at.lvec is not None, f"No lattice vectors found in {args.xyz}. Comment line must start with 'lvs'."
    assert sys_at.qs   is not None, f"No charges found in {args.xyz}. 4th column must contain charges."

    prefix = args.prefix or os.path.splitext(os.path.basename(args.xyz))[0]

    # --- Build Ewald2D from AtomicSystem ---
    ew = Ewald2D.from_AtomicSystem(sys_at, n_harm=args.n_harm)
    ew.print_info()

    Ng = args.grid
    X_xy, Y_xy = ew.make_xy_grid(Ng)
    extent_xy = [0, ew.Lx, 0, ew.Ly]

    # ================================================================
    #  FIGURE 1: Charge density projection
    # ================================================================
    rho = ew.charge_density_xy(X_xy, Y_xy)

    fig1, axes1 = plt.subplots(1, 2, figsize=(12, 5))
    plot_ions_xy(axes1[0], ew.rx, ew.ry, ew.q, ew.a_vec, ew.b_vec, ew.Lx, ew.Ly)
    axes1[0].set_title('Ion positions (red=+, blue=−)')
    plot_imshow_sym(axes1[1], rho, extent_xy, f'Charge density (N_G={ew.N_G})', cbar_label='ρ (e/Å²)')
    fig1.tight_layout()
    fig1.savefig(f"{prefix}_fig1_charge_density.png", dpi=150); print(f"Saved {prefix}_fig1_charge_density.png")

    # ================================================================
    #  FIGURE 2: XY potential slices at multiple heights
    # ================================================================
    z_heights = [ew.z_max + 0.5, ew.z_max + 1.0, ew.z_max + 2.0]
    fig2, axes2 = plt.subplots(1, len(z_heights), figsize=(5*len(z_heights), 5))
    if len(z_heights) == 1: axes2 = [axes2]
    for ax, zh in zip(axes2, z_heights):
        phi = ew.phi_vacuum_xy(X_xy, Y_xy, zh)
        plot_imshow_sym(ax, phi, extent_xy, f'φ at z = {zh:.1f} Å', cbar_label='φ (e/Å)')
    fig2.tight_layout()
    fig2.savefig(f"{prefix}_fig2_xy_slices.png", dpi=150); print(f"Saved {prefix}_fig2_xy_slices.png")

    # ================================================================
    #  FIGURE 3: XZ potential slice
    # ================================================================
    X_xz, Y_xz, Z_xz = ew.make_xz_grid(Ng)
    phi_xz = ew.phi_full_2d(X_xz, Y_xz, Z_xz)

    fig3, ax3 = plt.subplots(1, 1, figsize=(10, 6))
    vmax = np.percentile(np.abs(phi_xz), 98) or 1.0
    im = ax3.pcolormesh(X_xz, Z_xz, phi_xz, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
    # mark ions near y=0
    tol = 0.1 * ew.Ly
    mask = np.abs(ew.ry) < tol
    if np.any(mask):
        ion_colors = ['red' if qi > 0 else 'blue' for qi in ew.q[mask]]
        ax3.scatter(ew.rx[mask], ew.rz[mask], c=ion_colors, s=80, edgecolors='k', zorder=5)
    ax3.set_xlabel('x (Å)'); ax3.set_ylabel('z (Å)')
    ax3.set_title('Potential in XZ plane (y=0)')
    plt.colorbar(im, ax=ax3, label='φ (e/Å)')
    fig3.tight_layout()
    fig3.savefig(f"{prefix}_fig3_xz_slice.png", dpi=150); print(f"Saved {prefix}_fig3_xz_slice.png")

    # ================================================================
    #  FIGURE 4: 1D z-scan comparison (2D Ewald vs brute force)
    # ================================================================
    # Compare DIFFERENCES φ(point)−φ(ref) to eliminate the arbitrary
    # constant from conditional convergence of the 2D Coulomb sum.
    test_points = [
        (0.0,              0.0,             "above origin"),
        (ew.a_vec[0]/2,    0.0,             "above (a/2,0)"),
        (ew.a_vec[0]/4,    ew.b_vec[1]/4,   "hollow site"),
    ]
    ref_xy = (ew.a_vec[0]/4, ew.b_vec[1]/4)
    z_scan = np.linspace(ew.z_max + 0.3, ew.z_max + 5.0, 120)

    phi_2d_ref    = ew.phi_full_1d(ref_xy[0], ref_xy[1], z_scan)
    phi_brute_ref = ew.phi_brute_1d(ref_xy[0], ref_xy[1], z_scan, N_rep=args.N_rep)

    fig4, axes4 = plt.subplots(2, len(test_points), figsize=(5*len(test_points), 8))
    fig4.suptitle('Δφ = φ(point,z) − φ(hollow,z)', fontsize=12)

    for col, (x0, y0, label) in enumerate(test_points):
        print(f"  z-scan at ({x0:.2f}, {y0:.2f}) '{label}' ... ", end="", flush=True)
        dphi_2d    = ew.phi_full_1d(x0, y0, z_scan)   - phi_2d_ref
        dphi_brute = ew.phi_brute_1d(x0, y0, z_scan, N_rep=args.N_rep) - phi_brute_ref
        err = dphi_2d - dphi_brute
        print(f"max|err| = {np.max(np.abs(err)):.2e}")

        ax = axes4[0, col]
        ax.plot(z_scan, dphi_brute, 'k-',  lw=1.5, label='brute force')
        ax.plot(z_scan, dphi_2d,    'r--', lw=1.0, label=f'2D Ewald (N_G={ew.N_G})')
        ax.set_xlabel('z (Å)'); ax.set_ylabel('Δφ (e/Å)')
        ax.set_title(label); ax.legend(fontsize=8); ax.grid(True, alpha=0.3)

        ax = axes4[1, col]
        ax.plot(z_scan, err, 'b-', lw=1.0)
        ax.set_xlabel('z (Å)'); ax.set_ylabel('Δφ_2D − Δφ_brute (e/Å)')
        ax.set_title(f'Error  max|Δ|={np.max(np.abs(err)):.2e}')
        ax.grid(True, alpha=0.3); ax.axhline(0, color='k', lw=0.5)

    fig4.tight_layout()
    fig4.savefig(f"{prefix}_fig4_zscan.png", dpi=150); print(f"Saved {prefix}_fig4_zscan.png")

    # ================================================================
    #  FIGURE 5: Convergence study — error vs n_harm
    # ================================================================
    x0_A, y0_A = 0.0, 0.0
    x0_B, y0_B = ew.a_vec[0]/4, ew.b_vec[1]/4
    z_conv = ew.z_max + 1.0
    z1 = np.array([z_conv])

    # reference at large n_harm
    ew_ref = Ewald2D(ew.a_vec, ew.b_vec, ew.rx, ew.ry, ew.rz, ew.q, n_harm=12)
    dphi_ref = ew_ref.phi_full_1d(x0_A, y0_A, z1)[0] - ew_ref.phi_full_1d(x0_B, y0_B, z1)[0]
    print(f"\nConvergence study: Δφ = φ(origin) − φ(hollow) at z={z_conv:.1f}")
    print(f"  Reference (n_harm=12): Δφ_ref = {dphi_ref:.10f}")

    n_harm_list = [1, 2, 3, 4, 5, 6, 7, 8]
    errors_conv, ng_list = [], []
    for nh in n_harm_list:
        ew_t = Ewald2D(ew.a_vec, ew.b_vec, ew.rx, ew.ry, ew.rz, ew.q, n_harm=nh)
        dphi_t = ew_t.phi_full_1d(x0_A, y0_A, z1)[0] - ew_t.phi_full_1d(x0_B, y0_B, z1)[0]
        err_t = abs(dphi_t - dphi_ref)
        errors_conv.append(err_t); ng_list.append(ew_t.N_G)
        print(f"  n_harm={nh:2d}  N_G={ew_t.N_G:4d}  Δφ={dphi_t:.10f}  |err|={err_t:.2e}")

    fig5, ax5 = plt.subplots(1, 1, figsize=(7, 5))
    ax5.semilogy(n_harm_list, errors_conv, 'bo-', lw=1.5, ms=6)
    ax5.set_xlabel('n_harm (G-vector half-width)')
    ax5.set_ylabel('|Δφ(n_harm) − Δφ(ref)| (e/Å)')
    ax5.set_title(f'Convergence of Δφ at z={z_conv:.1f} Å')
    ax5.grid(True, alpha=0.3)
    for nh, ng, e in zip(n_harm_list, ng_list, errors_conv):
        ax5.annotate(f'N_G={ng}', (nh, e), textcoords="offset points", xytext=(5, 5), fontsize=7)
    fig5.tight_layout()
    fig5.savefig(f"{prefix}_fig5_convergence.png", dpi=150); print(f"Saved {prefix}_fig5_convergence.png")

    if not args.noPlot:
        plt.show()
    print("Done.")
