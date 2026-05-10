#!/usr/bin/env python3
"""
plot_transition.py - Visualize the attractive->repulsive transition at multiple z-heights.

Produces:
1. Panel plot: rows=[Pauli, vdW, ES, Total], columns=[z distances]
2. 1D line plots: total energy vs x-position along molecule at each z

Usage:
    cd tests/tAFM/pyocl_fdbm
    python3 plot_transition.py --basis mio-1-1
    python3 plot_transition.py --basis 3ob-3-1
"""

import os, sys, argparse
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..')))


def atom_to_grid_idx(atom_pos, origin, step, ngrid):
    frac = (atom_pos - origin) / step
    idx = np.round(frac).astype(np.int32)
    idx = np.clip(idx, [0, 0, 0], np.array(ngrid) - 1)
    return idx


def plot_transition(src_dir, out_dir, basis, A, beta):
    os.makedirs(out_dir, exist_ok=True)

    # Load fields
    E_pauli = np.load(os.path.join(src_dir, 'E_Pauli_field.npy'))
    E_es = np.load(os.path.join(src_dir, 'E_es_field.npy'))
    E_vdw = np.load(os.path.join(src_dir, 'E_vdw_field.npy'))
    E_total = E_pauli + E_es + E_vdw

    origin = np.array([-11.36, -6.6, -4.2], dtype=np.float32)
    step = 0.15
    ngrid = E_total.shape

    # Atom positions
    atom_names, atom_pos = [], []
    with open(os.path.join(_THIS_DIR, 'pentacene.xyz'), 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    for line in lines[2:2+natoms]:
        p = line.split()
        atom_names.append(p[0])
        atom_pos.append([float(p[1]), float(p[2]), float(p[3])])
    atom_pos = np.array(atom_pos, dtype=np.float64)

    # Target atoms
    target_indices = [0, 1, 20, 21]
    z_list = [2.5, 2.8, 3.0, 3.2, 3.5]
    n_z = len(z_list)

    # z-grid values for extraction
    z_grid_vals = origin[2] + np.arange(ngrid[2]) * step

    # ── Panel plot: rows=fields, cols=z with colorbars ──────────────
    fields = [
        (E_pauli, 'Pauli', 'magma', False),
        (E_vdw, 'vdW', 'viridis', False),
        (E_es, 'Electrostatics', 'coolwarm', True),
        (E_total, 'Total', 'RdBu_r', True),
    ]
    n_rows = len(fields)
    fig, axes = plt.subplots(n_rows, n_z, figsize=(3.5 * n_z, 3.0 * n_rows))
    for irow, (field, title, cmap, sym) in enumerate(fields):
        for icol, z in enumerate(z_list):
            ax = axes[irow, icol]
            iz = int(np.clip(np.round((z - origin[2]) / step), 0, ngrid[2] - 1))
            sl = field[:, :, iz].T
            nx, ny = field.shape[:2]
            extent_xy = [float(origin[0]), float(origin[0]) + (nx - 1) * step,
                         float(origin[1]), float(origin[1]) + (ny - 1) * step]
            # Per-subplot vmin/vmax
            if sym:
                vabs = max(abs(float(sl.min())), abs(float(sl.max())), 1e-12)
                norm = plt.matplotlib.colors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
            else:
                vmin = float(sl.min()) if sl.min() < 0 else 0
                vmax = float(sl.max())
                norm = plt.matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
            im = ax.imshow(sl, origin='lower', cmap=cmap, norm=norm, extent=extent_xy, aspect='equal')
            # Colorbar on each subplot
            cb = fig.colorbar(im, ax=ax, shrink=0.6, pad=0.02)
            cb.ax.tick_params(labelsize=6)
            if irow == 0:
                ax.set_title(f'z={z:.1f}A', fontsize=9)
            if icol == 0:
                ax.set_ylabel(title, fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])
    fig.suptitle(f'{basis}: Pauli(A={A:.0f}, b={beta:.3f})  rows=fields, cols=z-height', fontsize=11)
    plt.tight_layout()
    fname = os.path.join(out_dir, 'transition_panels.png')
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")

    # ── Per-atom z-scan: Pauli, vdW, ES, Total vs z ───────────────────
    for target_idx in target_indices:
        tix, tiy, _ = atom_to_grid_idx(atom_pos[target_idx], origin, step, ngrid)
        z = z_grid_vals
        # Extract profiles at this grid column
        ep = E_pauli[tix, tiy, :]
        ev = E_vdw[tix, tiy, :]
        ee = E_es[tix, tiy, :]
        et = E_total[tix, tiy, :]

        fig, ax = plt.subplots(figsize=(8, 5))
        mask = (z >= 2.0) & (z <= 6.0)
        ax.plot(z[mask], ep[mask], '-', color='tab:red', lw=1.5, label='Pauli')
        ax.plot(z[mask], ev[mask], '-', color='tab:blue', lw=1.5, label='vdW')
        ax.plot(z[mask], ee[mask], '-', color='tab:green', lw=1.5, label='Electrostatics')
        ax.plot(z[mask], et[mask], '-', color='tab:purple', lw=2.0, label='Total')
        ax.axhline(0, color='k', lw=0.5, ls='--')
        ax.set_xlabel('z [Å]')
        ax.set_ylabel('Energy [eV]')
        ax.set_title(f'{basis}: atom {target_idx} ({atom_names[target_idx]}) at [{atom_pos[target_idx,0]:.3f}, {atom_pos[target_idx,1]:.3f}, {atom_pos[target_idx,2]:.3f}]')
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        ax.set_xlim([2.0, 6.0])
        ax.set_ylim([-0.5, 0.5])
        plt.tight_layout()
        fname = os.path.join(out_dir, f'zscan_atom{target_idx}.png')
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}")

    # ── Combined z-scan: all atoms on one plot ──────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    axes = axes.flatten()
    for iax, target_idx in enumerate(target_indices):
        ax = axes[iax]
        tix, tiy, _ = atom_to_grid_idx(atom_pos[target_idx], origin, step, ngrid)
        z = z_grid_vals
        ep = E_pauli[tix, tiy, :]
        ev = E_vdw[tix, tiy, :]
        ee = E_es[tix, tiy, :]
        et = E_total[tix, tiy, :]
        mask = (z >= 2.0) & (z <= 6.0)
        ax.plot(z[mask], ep[mask], '-', color='tab:red', lw=1.2, label='Pauli')
        ax.plot(z[mask], ev[mask], '-', color='tab:blue', lw=1.2, label='vdW')
        ax.plot(z[mask], ee[mask], '-', color='tab:green', lw=1.2, label='ES')
        ax.plot(z[mask], et[mask], '-', color='tab:purple', lw=2.0, label='Total')
        ax.axhline(0, color='k', lw=0.5, ls='--')
        ax.set_xlabel('z [Å]')
        ax.set_ylabel('Energy [eV]')
        ax.set_title(f'atom {target_idx} ({atom_names[target_idx]})')
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([2.0, 6.0])
        ax.set_ylim([-0.5, 0.5])
    fig.suptitle(f'{basis}: z-scan components at 4 peripheral carbons (A={A:.0f}, b={beta:.3f})', fontsize=12)
    plt.tight_layout()
    fname = os.path.join(out_dir, 'zscan_all_atoms.png')
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")

    # ── Max energy vs z (keep this) ─────────────────────────────────
    z_fine = np.arange(2.0, 5.0, 0.1)
    ny = E_total.shape[1]
    iy = ny // 2
    max_total, max_pauli, max_vdw = [], [], []
    for zv in z_fine:
        iz = int(np.clip(np.round((zv - origin[2]) / step), 0, E_total.shape[2] - 1))
        ix_min = int(np.clip(np.round((-6 - origin[0]) / step), 0, E_total.shape[0] - 1))
        ix_max = int(np.clip(np.round((6 - origin[0]) / step), 0, E_total.shape[0] - 1))
        region = slice(ix_min, ix_max + 1)
        max_total.append(E_total[region, iy, iz].max())
        max_pauli.append(E_pauli[region, iy, iz].max())
        max_vdw.append(E_vdw[region, iy, iz].min())

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.semilogy(z_fine, max_pauli, 'o-', color='tab:red', label='Max Pauli', markersize=4)
    ax.semilogy(z_fine, np.abs(max_vdw), 's-', color='tab:blue', label='|Max vdW|', markersize=4)
    ax.semilogy(z_fine, np.array(max_total).clip(1e-6, None), '^-', color='tab:green', label='Max Total', markersize=4)
    ax.axvline(2.8, color='gray', ls='--', alpha=0.5)
    ax.axvline(3.2, color='gray', ls='--', alpha=0.5)
    ax.set_xlabel('z [Å]')
    ax.set_ylabel('|Energy| [eV]')
    ax.set_title(f'{basis}: Max energy above molecule vs tip height')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fname = os.path.join(out_dir, 'maxE_vs_z.png')
    plt.savefig(fname, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {fname}")

    print(f"\nAll transition plots in: {out_dir}/")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'])
    args = parser.parse_args()

    FITTED = {'mio-1-1': {'A': 965.0, 'beta': 0.871}, '3ob-3-1': {'A': 643.0, 'beta': 0.796}}
    src_dir = os.path.join(_THIS_DIR, 'afm_fitted_' + args.basis.replace('-', '_'))
    out_dir = os.path.join(src_dir, 'transition')
    plot_transition(src_dir, out_dir, args.basis, FITTED[args.basis]['A'], FITTED[args.basis]['beta'])


if __name__ == "__main__":
    main()
