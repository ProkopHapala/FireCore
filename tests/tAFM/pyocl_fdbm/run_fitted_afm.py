#!/usr/bin/env python3
"""
run_fitted_afm.py - Rerun AFM simulation (steps 3-6) using fitted Pauli parameters.

Loads existing density grids from FDBM pipeline output, recomputes Pauli with
fitted A and beta, then generates AFM images.

Usage:
    cd tests/tAFM/pyocl_fdbm
    python3 run_fitted_afm.py --basis mio-1-1
    python3 run_fitted_afm.py --basis 3ob-3-1 --xyz mymol.xyz --src_dir my_fdbm_output
"""

import os, sys, argparse, time
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall import atomicUtils as au
from pyBall.OCL import AFM as afm
from pyBall import plotUtils as pu

def run_afm_with_fitted_params(src_dir, out_dir, A_pauli, beta_pauli, atomPos, xyz_path, step=0.15):
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"Fitted AFM: src={src_dir} -> out={out_dir}")
    print(f"A_pauli={A_pauli:.2f}, beta_pauli={beta_pauli:.4f}")
    print(f"{'='*60}")

    # Load existing grids
    rho_grid = np.load(os.path.join(src_dir, 'step1_density', 'rho_grid.npy'))
    V_ES = np.load(os.path.join(src_dir, 'step2_electrostatics', 'V_ES.npy'))
    co_rho_total = np.load(os.path.join(src_dir, 'co_tip', 'co_rho_total.npy'))
    co_rho_delta = np.load(os.path.join(src_dir, 'co_tip', 'co_rho_delta.npy'))

    # Grid spec from FDBM log
    origin, ngrid, step_read = afm.read_grid_spec_from_log(os.path.join(src_dir, 'step1_density', 'log.txt'))
    if origin is None:
        raise FileNotFoundError(f"Could not read grid spec from {src_dir}/step1_density/log.txt")
    if step_read is not None:
        step = step_read
    print(f"Loaded grids: rho={rho_grid.shape}, V_ES={V_ES.shape}")
    print(f"Grid: origin={origin}, ngrid={ngrid}, step={step}")

    # Step 3: Pauli with fitted parameters
    t0 = time.time()
    E_pauli_field, grads_pauli = afm.compute_pauli_field(rho_grid, co_rho_total, step, A_pauli=A_pauli, beta_pauli=beta_pauli)
    np.save(os.path.join(out_dir, 'E_Pauli_field.npy'), E_pauli_field)
    print(f"  Pauli field range: [{E_pauli_field.min():.4f}, {E_pauli_field.max():.4f}] eV")
    print(f"  Step 3 (Pauli): {time.time()-t0:.2f}s")

    # Step 4: ES convolution
    t0 = time.time()
    E_es_field, grads_es = afm.compute_es_conv_field(V_ES, co_rho_delta, step)
    np.save(os.path.join(out_dir, 'E_es_field.npy'), E_es_field)
    print(f"  Step 4 (ES): {time.time()-t0:.2f}s")

    # Step 5: vdW
    t0 = time.time()
    # Use C6 table from run_fitted_afm.py original for exact backward compat
    C6_table = {1: 6.5, 6: 24.0, 7: 20.0, 8: 15.0}
    elem_z = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
    atomTypes = np.array([elem_z.get(e, 6) for e in au.load_xyz(xyz_path)[2]])
    E_vdw, grads_vdw = afm.compute_vdw_field(atomPos, atomTypes, origin, step, ngrid, C6_table=C6_table, C6_CO=30.0, RA=1.5)
    np.save(os.path.join(out_dir, 'E_vdw_field.npy'), E_vdw)
    print(f"  Step 5 (vdW): {time.time()-t0:.2f}s")

    # Step 6: AFM images
    t0 = time.time()
    scan_xs = np.linspace(atomPos[:,0].min()-3, atomPos[:,0].max()+3, 100)
    scan_ys = np.linspace(atomPos[:,1].min()-3, atomPos[:,1].max()+3, 100)
    heights = np.arange(3.0, 5.5, 0.1)

    nx, ny, nz = ngrid
    xs = origin[0] + np.arange(nx) * step
    ys = origin[1] + np.arange(ny) * step
    zs = origin[2] + np.arange(nz) * step
    F_total = grads_pauli + grads_es + grads_vdw

    from scipy.ndimage import map_coordinates
    def force_func(positions):
        query = np.stack([
            (positions[:,0] - xs[0]) / step,
            (positions[:,1] - ys[0]) / step,
            (positions[:,2] - zs[0]) / step
        ], axis=1)
        fx = map_coordinates(F_total[:,:,:,0], query.T, order=1, mode='constant', cval=0.0)
        fy = map_coordinates(F_total[:,:,:,1], query.T, order=1, mode='constant', cval=0.0)
        fz = map_coordinates(F_total[:,:,:,2], query.T, order=1, mode='constant', cval=0.0)
        return np.stack([fx, fy, fz], axis=1)

    FEs_relax = afm.pp_relax_2d(force_func, scan_xs, scan_ys, heights, mol_z=atomPos[:,2].max(), K_LAT=0.5, N_RELAX=50, step=step)
    Fz_relax = FEs_relax[:,:,:,2]
    df = afm.compute_df(Fz_relax, heights[1]-heights[0])
    np.save(os.path.join(out_dir, 'df.npy'), df)
    print(f"  df shape: {df.shape}, range: [{df.min():.4f}, {df.max():.4f}]")

    afm.save_afm_images(df, scan_xs, scan_ys, heights, out_dir, prefix='df')

    # Pauli field slices using shared plot helper
    z_slice = 2.0
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Pauli Energy (A={A_pauli:.1f}, b={beta_pauli:.3f})")
    pu.plot_field_slice(axes[0], E_pauli_field, origin, step, z_slice, cmap='magma', title=f"XY z={z_slice:.1f}A")
    pu.plot_field_slice(axes[1], E_pauli_field, origin, step, z_slice, cmap='magma', title="XZ y-center")
    pu.plot_field_slice(axes[2], E_pauli_field, origin, step, z_slice, cmap='magma', title="YZ x-center")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'pauli_slices.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {os.path.join(out_dir, 'pauli_slices.png')}")

    print(f"\nAll outputs in: {out_dir}/")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'])
    parser.add_argument('--xyz', type=str, default='pentacene.xyz', help='Molecule XYZ file')
    parser.add_argument('--src_dir', type=str, default=None, help='FDBM pipeline output directory')
    parser.add_argument('--A_pauli', type=float, default=None)
    parser.add_argument('--beta_pauli', type=float, default=None)
    args = parser.parse_args()

    A_pauli = args.A_pauli if args.A_pauli is not None else afm.PAULI_FITTED_DEFAULTS[args.basis]['A']
    beta_pauli = args.beta_pauli if args.beta_pauli is not None else afm.PAULI_FITTED_DEFAULTS[args.basis]['beta']

    if args.src_dir:
        src_dir = os.path.join(_THIS_DIR, args.src_dir)
    else:
        src_dir = os.path.join(_THIS_DIR, f'debug_dftb_{args.basis.replace("-", "_")}')
    out_dir = os.path.join(_THIS_DIR, f'afm_fitted_{args.basis.replace("-", "_")}')

    xyz_path = os.path.join(_THIS_DIR, args.xyz)
    mol_pos, _, mol_names, _, _ = au.load_xyz(xyz_path)
    atomPos = np.array(mol_pos, dtype=np.float64)
    run_afm_with_fitted_params(src_dir, out_dir, A_pauli, beta_pauli, atomPos, xyz_path)

if __name__ == "__main__":
    main()
