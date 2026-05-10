#!/usr/bin/env python3
"""
diagnostic_forcefield.py - Plot Pauli, vdW, total potential and Fz at relevant z-heights
where repulsive atoms should be visible.

Usage:
    cd tests/tAFM/pyocl_fdbm
    python3 diagnostic_forcefield.py --basis mio-1-1
    python3 diagnostic_forcefield.py --basis 3ob-3-1
"""

import os, sys, argparse
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))

sys.path.insert(0, os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..')))

def plot_field_slice(ax, field, origin, step, z, cmap='magma', title='', sym=False):
    nx, ny, nz = field.shape
    iz = int(np.clip(np.round((z - origin[2]) / step), 0, nz-1))
    slice_xy = field[:, :, iz].T
    extent_xy = [float(origin[0]), float(origin[0]) + (nx-1)*step,
                 float(origin[1]), float(origin[1]) + (ny-1)*step]
    if sym:
        vabs = max(abs(float(slice_xy.min())), abs(float(slice_xy.max())), 1e-12)
        norm = plt.matplotlib.colors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
    else:
        norm = None
    im = ax.imshow(slice_xy, origin='lower', cmap=cmap, norm=norm, extent=extent_xy, aspect='equal')
    ax.set_title(title)
    plt.colorbar(im, ax=ax, shrink=0.6)
    return iz

def diagnostic_plots(basis, A_pauli, beta_pauli):
    src_dir = os.path.join(_THIS_DIR, 'afm_fitted_' + basis.replace('-', '_'))
    out_dir = os.path.join(src_dir, 'diagnostic')
    os.makedirs(out_dir, exist_ok=True)

    # Load fields
    E_pauli = np.load(os.path.join(src_dir, 'E_Pauli_field.npy'))
    E_es = np.load(os.path.join(src_dir, 'E_es_field.npy'))
    E_vdw = np.load(os.path.join(src_dir, 'E_vdw_field.npy'))
    E_total = E_pauli + E_es + E_vdw

    # Grid spec
    origin = np.array([-11.36, -6.6, -4.2], dtype=np.float32)
    step = 0.15

    # Compute Fz = -dE_total/dz
    Fz = -np.gradient(E_total, step, axis=2)

    z_heights = [2.0, 2.5, 3.0]
    n_z = len(z_heights)

    print(f"Diagnostic plots for {basis} (A={A_pauli:.1f}, b={beta_pauli:.3f})")
    print(f"Grid: origin={origin}, step={step}")
    for z in z_heights:
        iz = int(np.clip(np.round((z - origin[2]) / step), 0, E_total.shape[2]-1))
        print(f"  z={z:.1f}A (iz={iz}): Pauli=[{E_pauli[:,:,iz].min():.3f}, {E_pauli[:,:,iz].max():.3f}] eV")
        print(f"    vdW=[{E_vdw[:,:,iz].min():.3f}, {E_vdw[:,:,iz].max():.3f}], total=[{E_total[:,:,iz].min():.3f}, {E_total[:,:,iz].max():.3f}]")
        print(f"    Fz=[{Fz[:,:,iz].min():.3f}, {Fz[:,:,iz].max():.3f}] eV/A")

    # ── Panel 1: Pauli potential ────────────────────────────────────────
    fig, axes = plt.subplots(1, n_z, figsize=(5*n_z, 4.5))
    if n_z == 1: axes = [axes]
    for ax, z in zip(axes, z_heights):
        iz = plot_field_slice(ax, E_pauli, origin, step, z, cmap='magma',
                              title=f'Pauli z={z:.1f}A')
    fig.suptitle(f'{basis}: Pauli Repulsion (A={A_pauli:.1f}, b={beta_pauli:.3f})', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'panel_pauli.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_dir}/panel_pauli.png")

    # ── Panel 2: vdW potential ────────────────────────────────────────
    fig, axes = plt.subplots(1, n_z, figsize=(5*n_z, 4.5))
    if n_z == 1: axes = [axes]
    for ax, z in zip(axes, z_heights):
        iz = plot_field_slice(ax, E_vdw, origin, step, z, cmap='viridis',
                              title=f'vdW z={z:.1f}A')
    fig.suptitle(f'{basis}: vdW Attraction', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'panel_vdw.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_dir}/panel_vdw.png")

    # ── Panel 3: Total potential (symmetric colormap) ─────────────────
    fig, axes = plt.subplots(1, n_z, figsize=(5*n_z, 4.5))
    if n_z == 1: axes = [axes]
    for ax, z in zip(axes, z_heights):
        iz = plot_field_slice(ax, E_total, origin, step, z, cmap='RdBu_r',
                              title=f'Total z={z:.1f}A', sym=True)
    fig.suptitle(f'{basis}: Total Potential (red=repulsive, blue=attractive)', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'panel_total.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_dir}/panel_total.png")

    # ── Panel 4: Fz (force in z, sym colormap) ────────────────────────
    fig, axes = plt.subplots(1, n_z, figsize=(5*n_z, 4.5))
    if n_z == 1: axes = [axes]
    for ax, z in zip(axes, z_heights):
        iz = plot_field_slice(ax, Fz, origin, step, z, cmap='RdBu_r',
                              title=f'Fz z={z:.1f}A', sym=True)
    fig.suptitle(f'{basis}: Fz (red=push-up/repulsive, blue=pull-down/attractive)', fontsize=12)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'panel_Fz.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_dir}/panel_Fz.png")

    # ── Combined comparison at z=2.5A (most informative) ─────────────
    z_target = 2.5
    iz = int(np.clip(np.round((z_target - origin[2]) / step), 0, E_total.shape[2]-1))
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, field, cmap, sym, title in [
        (axes[0,0], E_pauli, 'magma', False, 'Pauli'),
        (axes[0,1], E_vdw, 'viridis', False, 'vdW'),
        (axes[1,0], E_total, 'RdBu_r', True, 'Total'),
        (axes[1,1], Fz, 'RdBu_r', True, 'Fz'),
    ]:
        slice_xy = field[:, :, iz].T
        nx, ny, nz = field.shape
        extent_xy = [float(origin[0]), float(origin[0]) + (nx-1)*step,
                     float(origin[1]), float(origin[1]) + (ny-1)*step]
        if sym:
            vabs = max(abs(float(slice_xy.min())), abs(float(slice_xy.max())), 1e-12)
            norm = plt.matplotlib.colors.TwoSlopeNorm(vmin=-vabs, vcenter=0, vmax=vabs)
        else:
            norm = None
        im = ax.imshow(slice_xy, origin='lower', cmap=cmap, norm=norm, extent=extent_xy, aspect='equal')
        ax.set_title(f'{title} z={z_target:.1f}A')
        plt.colorbar(im, ax=ax, shrink=0.6)
    fig.suptitle(f'{basis}: Forcefield at z={z_target:.1f}A (A={A_pauli:.1f}, b={beta_pauli:.3f})', fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'combined_z2.5.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_dir}/combined_z2.5.png")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'])
    args = parser.parse_args()

    FITTED = {'mio-1-1': {'A': 965.0, 'beta': 0.871}, '3ob-3-1': {'A': 643.0, 'beta': 0.796}}
    diagnostic_plots(args.basis, FITTED[args.basis]['A'], FITTED[args.basis]['beta'])

if __name__ == "__main__":
    main()
