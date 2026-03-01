"""
SDF Sampling Grid — Fixed-Endpoint Isolines
=============================================

GEOMETRY (4 atoms, NaCl step edge cross-section):

        *         x=0,  z=a   — corner atom (top of step)
       o*o        x=-a, 0, +a, z=0  — lower terrace atoms

BOUNDARY LINES:
    L : x = -a   (above left lower-terrace atom)
    R : x = +a   (above right lower-terrace atom)

GRID CONSTRUCTION:
    For each SDF shell u_k:
      1. Walk the isoline SDF=u_k at high resolution (small arc-length steps)
         to get accurate arc-length parametrization
      2. Find the two crossing points with x=L and x=R (the "entry" and "exit")
      3. Extract the arc segment between L and R
      4. Distribute exactly n_v samples uniformly along that segment
         => every isoline has the SAME n_v samples, same indexing

    This gives a true [n_u × n_v] grid with integer (i,j) indexing:
      i = u-index (which shell, 0 = closest to surface)
      j = v-index (which column, 0 = left boundary, n_v-1 = right boundary)

    "Streamlines" = columns of the grid: samples[:,j] across all i values
    "Isolines"    = rows of the grid:    samples[i,:] across all j values

    The streamlines are NOT integrated — they are simply the column index j
    shared across all isolines. This avoids all branching issues.

PARAMETERS:
    --n_sub      : number of samples per isoline (subdivisions between L and R)
    --u_levels   : number of SDF shells
    --u_schedule : 'sinh' (dense near surface), 'linear', 'geom'
    --corner     : 'Na' or 'Cl' (which species is the corner atom)
    --walk_res   : arc-length step for isoline walking (Å, default 0.005)
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches

from sampling_lib import (
    build_grid_uv,
    compute_sdf_2d,
    make_atoms_nacl_step,
    make_u_levels,
    parse_u_list,
)


# ─────────────────────────────────────────────────────────────────────────────
# 1.  ARGS
# ─────────────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('--a',           type=float, default=2.82065)
    p.add_argument('--r_Na',        type=float, default=1.4)
    p.add_argument('--r_Cl',        type=float, default=2.0)
    p.add_argument('--corner',      type=str,   default='both', choices=['Na', 'Cl', 'both'])
    p.add_argument('--n_sub',       type=int,   default=9,     help='samples per isoline between L and R (n_v)')
    p.add_argument('--u_levels',    type=int,   default=20)
    p.add_argument('--u_min',       type=float, default=0.05)
    p.add_argument('--u_max',       type=float, default=5.0)
    p.add_argument('--u_schedule',  type=str,   default='sinh', choices=['linear', 'sinh', 'geom'])
    p.add_argument('--u_sinh_a',    type=float, default=2.5)
    p.add_argument('--u_list',      type=str,   default=None, help='comma-separated explicit distances (e.g., 0,0.1,0.2,0.4,0.8,2,4,8); overrides u_levels/schedule')
    p.add_argument('--combine_p',   type=float, default=4.0, help='p-order to aggregate per-atom distances (p→inf → min, p=1 → average)')
    p.add_argument('--walk_res',    type=float, default=0.005,  help='arc-length step for isoline walk (Å)')
    p.add_argument('--sdf_res',     type=int,   default=1000,   help='SDF grid resolution')
    try:
        return p.parse_args(argv)
    except SystemExit:
        return p.parse_args([])


# ─────────────────────────────────────────────────────────────────────────────
# 2.  ATOMS  (exactly 4, minimal step-edge section)
# ─────────────────────────────────────────────────────────────────────────────

def make_atoms(a, r_Na, r_Cl, corner_species='Na'):
    return make_atoms_nacl_step(a, r_Na, r_Cl, corner_species=corner_species)


# ─────────────────────────────────────────────────────────────────────────────
# 3.  SDF
# ─────────────────────────────────────────────────────────────────────────────

def compute_sdf(atoms, x_range, z_range, res, combine_p):
    return compute_sdf_2d(atoms, x_range, z_range, res, combine_p)


# ─────────────────────────────────────────────────────────────────────────────
# 4.  U-SCHEDULE
# ─────────────────────────────────────────────────────────────────────────────


# ─────────────────────────────────────────────────────────────────────────────
# 5.  ISOLINE WALKING  (the core of this design)
# ─────────────────────────────────────────────────────────────────────────────

def build_grid(sdf, x_lin, z_lin, u_levels, x_left, x_right, n_v, walk_res):
    return build_grid_uv(sdf, x_lin, z_lin, u_levels, x_left, x_right, n_v, walk_res)


# ─────────────────────────────────────────────────────────────────────────────
# 7.  PLOT
# ─────────────────────────────────────────────────────────────────────────────

def plot_grid(fig, ax_left, ax_right,atoms, X_grid, Z_grid, sdf, x_lin, z_lin,grid, u_levels_valid, arc_lengths,x_left, x_right, corner_species, args):

    a = args.a
    n_u, n_v = grid.shape[:2]

    def sdf_bg(ax):
        sdf_disp = np.clip(sdf, -1.8, args.u_max * 1.1)
        ax.imshow(sdf_disp,extent=[x_lin[0], x_lin[-1], z_lin[0], z_lin[-1]],origin='lower', cmap='RdYlGn', alpha=0.22,vmin=-1.6, vmax=args.u_max * 0.8,zorder=1, aspect='auto')
        ax.contour(X_grid, Z_grid, sdf, levels=u_levels_valid,colors='#1e3a5f', linewidths=0.4, alpha=0.25, zorder=2)
        ax.contour(X_grid, Z_grid, sdf, levels=[0.0],colors='#111827', linewidths=1.8, alpha=0.9,linestyles='--', zorder=3)

    def draw_atoms(ax):
        for at in atoms:
            color = '#3b82f6' if at['species'] == 'Na' else '#f97316'
            ec    = '#1e40af' if at['species'] == 'Na' else '#9a3412'
            alpha = 0.9
            ax.scatter(at['x'], at['z'], s=at['r']**2 * 60, c=color,edgecolors=ec, linewidth=0.9,zorder=8, alpha=alpha)
            ax.text(at['x'], at['z'],at['species'],ha='center', va='center',fontsize=6, color='white',fontweight='bold', zorder=9)

    def draw_boundary_lines(ax):
        zlo, zhi = z_lin[0], z_lin[-1]
        ax.axvline(x_left,  color='#7c3aed', lw=1.2, ls='-.',alpha=0.7, zorder=4, label=f'x={x_left/a:.0f}a')
        ax.axvline(x_right, color='#dc2626', lw=1.2, ls='-.',alpha=0.7, zorder=4, label=f'x={x_right/a:.0f}a')

    # ── Left panel: full grid ──────────────────────────────────────────────
    ax = ax_left
    sdf_bg(ax)
    draw_atoms(ax)
    draw_boundary_lines(ax)

    # Draw every isoline as a line through its samples
    for i in range(n_u):
        row = grid[i]   # (n_v, 2)
        ax.plot(row[:, 0], row[:, 1],color='#1d4ed8', lw=0.8, alpha=0.35, zorder=5)

    # Draw every streamline (column of grid)
    for j in range(n_v):
        col = grid[:, j]   # (n_u, 2)
        ax.plot(col[:, 0], col[:, 1], color='#0891b2', lw=0.6, alpha=0.3, zorder=5)

    # Scatter all grid points, colored by u_idx
    xs = grid[:, :, 0].ravel()
    zs = grid[:, :, 1].ravel()
    ui = np.repeat(np.arange(n_u), n_v)
    sc = ax.scatter(xs, zs, s=12, c=ui, cmap='plasma_r',alpha=0.8, zorder=10, edgecolors='none', vmin=0, vmax=n_u)
    plt.colorbar(sc, ax=ax, fraction=0.03, pad=0.02, label='u_idx (shell index)')

    # Annotate arc-length range
    ax.text(0.02, 0.98,f"Arc-length: {arc_lengths[0]:.2f}Å (inner) → {arc_lengths[-1]:.2f}Å (outer)",transform=ax.transAxes, va='top', fontsize=7,color='#374151', bbox=dict(fc='white', alpha=0.7, pad=2))

    ax.set_title(f"Corner: {corner_species}  |  Grid: {n_u}×{n_v}\n"f"u: {args.u_schedule}(a={args.u_sinh_a}), "f"v: {n_v} uniform samples per isoline",fontsize=9)

    # ── Right panel: selected 1D cuts ─────────────────────────────────────
    ax = ax_right
    sdf_bg(ax)
    draw_atoms(ax)
    draw_boundary_lines(ax)

    # Faint background dots
    ax.scatter(xs, zs, s=5, c='#fca5a5', alpha=0.2,
               zorder=6, edgecolors='none')

    # 3 isoline cuts (rows)
    iso_colors = ['#15803d', '#16a34a', '#4ade80']
    iso_idxs   = [1, n_u//3, 2*n_u//3]
    for i_row, (i, c) in enumerate(zip(iso_idxs, iso_colors)):
        row = grid[i]
        ax.plot(row[:, 0], row[:, 1],color=c, lw=2.2, marker='o',ms=4, markerfacecolor='w', markeredgecolor=c,zorder=11,label=f'isoline i={i}  u={u_levels_valid[i]:.2f}Å')
        # Show first and last point explicitly
        ax.scatter([row[0, 0],  row[-1, 0]],[row[0, 1],  row[-1, 1]],s=40, c=c, zorder=13, marker='D')

    # 4 streamline cuts (columns)
    stream_colors = ['#f59e0b', '#06b6d4', '#7c3aed', '#f97316']
    j_targets = [0, n_v//3, 2*n_v//3, n_v-1]
    j_labels  = ['j=0 (left boundary)', f'j={n_v//3}', f'j={2*n_v//3}',
                 f'j={n_v-1} (right boundary)']
    for j, c, lbl in zip(j_targets, stream_colors, j_labels):
        col = grid[:, j]
        ax.plot(col[:, 0], col[:, 1],color=c, lw=2.2, marker='s',ms=4, markerfacecolor='w', markeredgecolor=c,zorder=11, label=f'stream {lbl}')

    ax.legend(loc='upper right', fontsize=6.5, framealpha=0.92,
              ncol=1)
    ax.set_title(
        f"1D cuts  —  Corner: {corner_species}\n"
        f"■ = streamlines (cols)  ●  = isolines (rows)  ◆ = endpoints",
        fontsize=9)

    for ax in [ax_left, ax_right]:
        ax.set_xlabel('x (Å)', fontsize=9)
        ax.set_ylabel('z (Å)', fontsize=9)
        ax.set_xlim(x_left  - 0.8*a, x_right + 0.8*a)
        ax.set_ylim(z_lin[0], min(z_lin[-1], a + args.u_max + 0.5))
        ax.set_aspect('equal')

    na_patch = mpatches.Patch(color='#3b82f6', label='Na')
    cl_patch = mpatches.Patch(color='#f97316', label='Cl')
    ax_left.legend(handles=[na_patch, cl_patch] +
                   [plt.Line2D([0],[0], color='#7c3aed', ls='-.', lw=1.2, label=f'boundary x={x_left/a:.0f}a'),
                    plt.Line2D([0],[0], color='#dc2626', ls='-.', lw=1.2, label=f'boundary x={x_right/a:.0f}a')],
                   loc='upper left', fontsize=7.5, framealpha=0.9)


# ─────────────────────────────────────────────────────────────────────────────
# 8.  MAIN
# ─────────────────────────────────────────────────────────────────────────────

def run_case(args, corner_species):
    a = args.a
    print(f"\n{'─'*50}")
    print(f"  Corner species: {corner_species}")
    print(f"{'─'*50}")

    atoms = make_atoms(a, args.r_Na, args.r_Cl, corner_species)
    for at in atoms:
        print(f"    {at['species']:2s}  x={at['x']:+.3f}  z={at['z']:.3f}  r={at['r']}")

    # Boundary lines
    x_left  = 0   # above left lower-terrace atom
    x_right = 2*a   # above right lower-terrace atom

    # Domain: generous padding for SDF computation
    pad = args.u_max + a
    x_range = (x_left  - pad, x_right + pad)
    z_range = (-a - 0.3, a + args.u_max + 1.5)

    # SDF
    X_grid, Z_grid, sdf, x_lin, z_lin = compute_sdf(
        atoms, x_range, z_range, args.sdf_res, args.combine_p)

    # U levels (explicit list overrides schedule)
    if args.u_list:
        u_levels = parse_u_list(args.u_list)
    else:
        u_levels = make_u_levels(args.u_min, args.u_max, args.u_levels,
                                  args.u_schedule, args.u_sinh_a)
    print(f"  u levels: {u_levels[0]:.4f} ... {u_levels[-1]:.4f} Å")
    print(f"  n_v = {args.n_sub} samples per isoline")
    print(f"  walk_res = {args.walk_res} Å")
    print(f"  Building grid...")

    grid, u_vals, u_idx_valid, arc_lengths = build_grid( sdf, x_lin, z_lin, u_levels, x_left, x_right, args.n_sub, args.walk_res)

    if grid is None:
        print("  ERROR: no valid isolines found")
        return None

    n_u, n_v = grid.shape[:2]
    print(f"\n  ✓ Grid shape: [{n_u} × {n_v}]  ({n_u*n_v} total samples)")
    print(f"  Arc-length range: {arc_lengths[0]:.3f} → {arc_lengths[-1]:.3f} Å")

    return (atoms, X_grid, Z_grid, sdf, x_lin, z_lin, grid, u_vals, arc_lengths, x_left, x_right)


def main(argv=None):
    args = parse_args(argv)
    a = args.a
    print("=" * 55)
    print(f"SDF Fixed-Endpoint Grid  —  n_v={args.n_sub}  n_u={args.u_levels}")
    print(f"u-schedule: {args.u_schedule}  (sinh_a={args.u_sinh_a})")
    print(f"Boundary lines: x = -a and x = +a  (Δx = 2a = {2*a:.4f} Å)")
    print("=" * 55)

    cases = ['Na', 'Cl'] if args.corner == 'both' else [args.corner]
    n_cases = len(cases)

    fig, axes = plt.subplots(n_cases, 2,  figsize=(13, 7 * n_cases),squeeze=False)

    for row_idx, corner_species in enumerate(cases):
        result = run_case(args, corner_species)
        if result is None: continue
        (atoms, X_grid, Z_grid, sdf, x_lin, z_lin,grid, u_vals, arc_lengths, x_left, x_right) = result
        plot_grid(fig, axes[row_idx][0], axes[row_idx][1], atoms, X_grid, Z_grid, sdf, x_lin, z_lin,  grid, u_vals, arc_lengths, x_left, x_right, corner_species, args)

    fig.suptitle(
        f"NaCl Step-Edge: Fixed-Endpoint SDF Sampling Grid\n"
        f"Boundary lines at x=−a and x=+a  |  "
        f"Every isoline: exactly {args.n_sub} uniform samples  |  "
        f"True [n_u × n_v] grid with (i,j) indexing\n"
        f"u-schedule: {args.u_schedule}(sinh_a={args.u_sinh_a})  |  "
        f"walk_res={args.walk_res}Å",
        fontsize=10, y=1.005
    )
    plt.tight_layout()

    out = os.path.join(os.path.dirname(__file__), 'sdf_fixed_endpoint.png')
    fig.savefig(out, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"\nSaved → {out}")
    plt.show()


if __name__ == '__main__':
    main()
