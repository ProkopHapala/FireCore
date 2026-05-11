#!/usr/bin/env python3
"""
diagnostic_forcefield.py - Plot Pauli, vdW, total potential and Fz at relevant z-heights.

Usage:
    cd tests/tAFM/pyocl_fdbm
    python3 diagnostic_forcefield.py --basis mio-1-1
    python3 diagnostic_forcefield.py --basis 3ob-3-1 --src_dir my_afm_output
"""

import os, sys, argparse
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL import AFM as afm
from pyBall import plotUtils as pu

def diagnostic_plots(src_dir, out_dir, basis, A_pauli, beta_pauli, step=0.15):
    os.makedirs(out_dir, exist_ok=True)

    # Load fields
    E_pauli = np.load(os.path.join(src_dir, 'E_Pauli_field.npy'))
    E_es = np.load(os.path.join(src_dir, 'E_es_field.npy'))
    E_vdw = np.load(os.path.join(src_dir, 'E_vdw_field.npy'))
    E_total = E_pauli + E_es + E_vdw

    # Try to read grid spec from src dir (step1_density/log.txt is sibling)
    fdbm_dir = os.path.dirname(src_dir) if os.path.basename(src_dir).startswith('afm_fitted_') else src_dir
    log_path = os.path.join(fdbm_dir, 'step1_density', 'log.txt')
    origin, ngrid, step_read = afm.read_grid_spec_from_log(log_path)
    if origin is None:
        # Fallback hardcoded for backward compat with existing outputs
        origin = np.array([-11.36, -6.6, -4.2], dtype=np.float32)
        ngrid = E_total.shape
        print(f"WARNING: Could not read grid from {log_path}, using fallback origin={origin}")
    if step_read is not None:
        step = step_read

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

    # Panels using shared plot helper
    for field, cmap, sym, prefix, title in [
        (E_pauli, 'magma', False, 'panel_pauli', f'{basis}: Pauli Repulsion (A={A_pauli:.1f}, b={beta_pauli:.3f})'),
        (E_vdw, 'viridis', False, 'panel_vdw', f'{basis}: vdW Attraction'),
        (E_total, 'RdBu_r', True, 'panel_total', f'{basis}: Total Potential (red=repulsive, blue=attractive)'),
        (Fz, 'RdBu_r', True, 'panel_Fz', f'{basis}: Fz (red=push-up/repulsive, blue=pull-down/attractive)'),
    ]:
        fig, axes = plt.subplots(1, n_z, figsize=(5*n_z, 4.5))
        if n_z == 1: axes = [axes]
        pu.plot_field_panel(axes, field, origin, step, z_heights, cmap=cmap, sym=sym)
        fig.suptitle(title, fontsize=12)
        plt.tight_layout()
        fname = os.path.join(out_dir, f'{prefix}.png')
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}")

    # Combined comparison at z=2.5A
    z_target = 2.5
    iz = int(np.clip(np.round((z_target - origin[2]) / step), 0, E_total.shape[2]-1))
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, field, cmap, sym, title in [
        (axes[0,0], E_pauli, 'magma', False, 'Pauli'),
        (axes[0,1], E_vdw, 'viridis', False, 'vdW'),
        (axes[1,0], E_total, 'RdBu_r', True, 'Total'),
        (axes[1,1], Fz, 'RdBu_r', True, 'Fz'),
    ]:
        pu.plot_field_slice(ax, field, origin, step, z_target, cmap=cmap, sym=sym,
                           title=f'{title} z={z_target:.1f}A')
    fig.suptitle(f'{basis}: Forcefield at z={z_target:.1f}A (A={A_pauli:.1f}, b={beta_pauli:.3f})', fontsize=13)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'combined_z2.5.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_dir}/combined_z2.5.png")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'])
    parser.add_argument('--src_dir', type=str, default=None, help='AFM fitted output directory')
    args = parser.parse_args()

    if args.src_dir:
        src_dir = os.path.join(_THIS_DIR, args.src_dir)
    else:
        src_dir = os.path.join(_THIS_DIR, f'afm_fitted_{args.basis.replace("-", "_")}')
    out_dir = os.path.join(src_dir, 'diagnostic')

    A = afm.PAULI_FITTED_DEFAULTS[args.basis]['A']
    beta = afm.PAULI_FITTED_DEFAULTS[args.basis]['beta']
    diagnostic_plots(src_dir, out_dir, args.basis, A, beta)

if __name__ == "__main__":
    main()
