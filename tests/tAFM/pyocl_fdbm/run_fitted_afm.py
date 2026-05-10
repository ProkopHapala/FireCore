#!/usr/bin/env python3
"""
run_fitted_afm.py - Rerun AFM simulation (steps 3-6) using fitted Pauli parameters.

Loads existing density grids from FDBM pipeline output, recomputes Pauli with
fitted A and beta, then generates AFM images.

Usage:
    cd tests/tAFM/pyocl_fdbm
    python3 run_fitted_afm.py --basis mio-1-1
    python3 run_fitted_afm.py --basis 3ob-3-1
"""

import os, sys, argparse, time
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_XYZ_PATH = os.path.join(_THIS_DIR, 'pentacene.xyz')

sys.path.insert(0, os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..')))
from pyBall.OCL import AFM as afm
from pyBall.OCL import clUtils as clu

def load_xyz(fname):
    with open(fname, 'r') as f:
        lines = f.readlines()
    natoms = int(lines[0].strip())
    enames, apos = [], []
    for line in lines[2:2+natoms]:
        p = line.split()
        enames.append(p[0])
        apos.append([float(p[1]), float(p[2]), float(p[3])])
    return np.array(enames), np.array(apos, dtype=np.float64)

def run_afm_with_fitted_params(src_dir, out_dir, A_pauli, beta_pauli, atomPos, step=0.15):
    os.makedirs(out_dir, exist_ok=True)
    print(f"\n{'='*60}")
    print(f"Fitted AFM: src={src_dir} -> out={out_dir}")
    print(f"A_pauli={A_pauli:.2f}, beta_pauli={beta_pauli:.4f}")
    print(f"{'='*60}")

    # Load existing grids
    rho_grid = np.load(os.path.join(src_dir, 'step1_density', 'rho_grid.npy'))
    rho_na_grid = np.load(os.path.join(src_dir, 'step1_density', 'rho_na_grid.npy'))
    rho_diff = np.load(os.path.join(src_dir, 'step1_density', 'rho_diff.npy'))
    V_ES = np.load(os.path.join(src_dir, 'step2_electrostatics', 'V_ES.npy'))
    co_rho_total = np.load(os.path.join(src_dir, 'co_tip', 'co_rho_total.npy'))
    co_rho_delta = np.load(os.path.join(src_dir, 'co_tip', 'co_rho_delta.npy'))

    # Grid spec from step1
    with open(os.path.join(src_dir, 'step1_density', 'log.txt'), 'r') as f:
        for line in f:
            if 'Grid:' in line:
                # Parse: "Grid: origin=[-11.36  -6.6   -4.2 ] ngrid=[152  88  96]"
                parts = line.split('ngrid=')
                origin_str = parts[0].split('origin=')[1].strip()
                ngrid_str = parts[1].strip()
                origin = np.array([float(x) for x in origin_str.strip('[]').split()])
                ngrid = np.array([int(x) for x in ngrid_str.strip('[]').split()])
                break

    print(f"Loaded grids: rho={rho_grid.shape}, V_ES={V_ES.shape}")
    print(f"Grid: origin={origin}, ngrid={ngrid}")

    # Step 3: Pauli with fitted parameters
    t0 = time.time()
    dV = step**3
    rho_tip_total = co_rho_total
    nx_t, ny_t, nz_t = rho_tip_total.shape
    tip_kernel = np.roll(np.roll(np.roll(rho_tip_total[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
    overlap_raw = dV * np.real(np.fft.ifftn(np.fft.fftn(rho_grid) * np.fft.fftn(tip_kernel))).astype(np.float32)
    overlap_safe = np.clip(overlap_raw, 1e-30, None)
    E_pauli_field = A_pauli * (overlap_safe ** beta_pauli)
    grads_pauli = np.stack([np.gradient(E_pauli_field, step, axis=i) for i in range(3)], axis=-1)
    np.save(os.path.join(out_dir, 'E_Pauli_field.npy'), E_pauli_field)
    print(f"  Pauli field range: [{E_pauli_field.min():.4f}, {E_pauli_field.max():.4f}] eV")
    print(f"  Step 3 (Pauli): {time.time()-t0:.2f}s")

    # Step 4: ES
    t0 = time.time()
    rho_tip_delta = co_rho_delta
    nx_t, ny_t, nz_t = rho_tip_delta.shape
    tip_kernel = np.roll(np.roll(np.roll(rho_tip_delta[::-1,::-1,::-1], -(nx_t//2), axis=0), -(ny_t//2), axis=1), -(nz_t//2), axis=2)
    E_es_field = dV * np.real(np.fft.ifftn(np.fft.fftn(V_ES) * np.fft.fftn(tip_kernel))).astype(np.float32)
    grads_es = np.stack([np.gradient(E_es_field, step, axis=i) for i in range(3)], axis=-1)
    np.save(os.path.join(out_dir, 'E_es_field.npy'), E_es_field)
    print(f"  Step 4 (ES): {time.time()-t0:.2f}s")

    # Step 5: vdW
    t0 = time.time()
    # Hardcoded C6 from original script
    C6_table = {1: 6.5, 6: 24.0, 7: 20.0, 8: 15.0}
    elem_z = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
    atomTypes, _ = load_xyz(_XYZ_PATH)
    atomPos_ang = atomPos
    E_vdw = np.zeros_like(rho_grid)
    for iat in range(len(atomTypes)):
        Z = C6_table.get(elem_z[atomTypes[iat]], 10.0)
        dx = np.arange(ngrid[0]) * step + origin[0] - atomPos_ang[iat, 0]
        dy = np.arange(ngrid[1]) * step + origin[1] - atomPos_ang[iat, 1]
        dz = np.arange(ngrid[2]) * step + origin[2] - atomPos_ang[iat, 2]
        DX, DY, DZ = np.meshgrid(dx, dy, dz, indexing='ij')
        r = np.sqrt(DX**2 + DY**2 + DZ**2)
        r_safe = np.clip(r, 2.0, None)
        E_vdw += -Z / (r_safe**6)
    grads_vdw = np.stack([np.gradient(E_vdw, step, axis=i) for i in range(3)], axis=-1)
    np.save(os.path.join(out_dir, 'E_vdw_field.npy'), E_vdw)
    print(f"  Step 5 (vdW): {time.time()-t0:.2f}s")

    # Step 6: AFM images
    t0 = time.time()
    scan_xs = np.linspace(atomPos[:,0].min()-3, atomPos[:,0].max()+3, 100)
    scan_ys = np.linspace(atomPos[:,1].min()-3, atomPos[:,1].max()+3, 100)
    heights = np.arange(3.0, 5.5, 0.1)

    # Interpolated force field
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

    # Save AFM images at different heights
    for i in [0, len(heights)//2, -1]:
        h = heights[i]
        plt.figure(figsize=(6,5))
        plt.imshow(df[:,:,i].T, origin='lower', extent=[scan_xs[0], scan_xs[-1], scan_ys[0], scan_ys[-1]], cmap='afmhot')
        plt.title(f"df at h={h:.1f} A")
        plt.colorbar(label="df [Hz]")
        fname = os.path.join(out_dir, f"df_h{h:.1f}.png")
        plt.savefig(fname, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  Saved: {fname}")

    # Also save Pauli field slices
    z_slice = 2.0
    iz = int(np.clip(np.round((z_slice - origin[2]) / step), 0, nz-1))
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(f"Pauli Energy (A={A_pauli:.1f}, b={beta_pauli:.3f})")
    # XY slice
    extent_xy = [origin[0], origin[0]+(nx-1)*step, origin[1], origin[1]+(ny-1)*step]
    im = axes[0].imshow(E_pauli_field[:,:,iz].T, origin='lower', cmap='magma', extent=extent_xy)
    axes[0].set_title(f"XY z={z_slice:.1f}A")
    fig.colorbar(im, ax=axes[0])
    # XZ center
    extent_xz = [origin[0], origin[0]+(nx-1)*step, origin[2], origin[2]+(nz-1)*step]
    im = axes[1].imshow(E_pauli_field[:,ny//2,:].T, origin='lower', cmap='magma', extent=extent_xz)
    axes[1].set_title("XZ y-center")
    fig.colorbar(im, ax=axes[1])
    # YZ center
    extent_yz = [origin[1], origin[1]+(ny-1)*step, origin[2], origin[2]+(nz-1)*step]
    im = axes[2].imshow(E_pauli_field[nx//2,:,:].T, origin='lower', cmap='magma', extent=extent_yz)
    axes[2].set_title("YZ x-center")
    fig.colorbar(im, ax=axes[2])
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'pauli_slices.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {os.path.join(out_dir, 'pauli_slices.png')}")

    print(f"\nAll outputs in: {out_dir}/")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--basis', type=str, default='mio-1-1', choices=['mio-1-1', '3ob-3-1'])
    parser.add_argument('--A_pauli', type=float, default=None)
    parser.add_argument('--beta_pauli', type=float, default=None)
    args = parser.parse_args()

    FITTED_DEFAULTS = {
        'mio-1-1': {'A': 965.0, 'beta': 0.871},
        '3ob-3-1': {'A': 643.0, 'beta': 0.796},
    }
    A_pauli = args.A_pauli if args.A_pauli is not None else FITTED_DEFAULTS[args.basis]['A']
    beta_pauli = args.beta_pauli if args.beta_pauli is not None else FITTED_DEFAULTS[args.basis]['beta']

    src_dir = os.path.join(_THIS_DIR, 'debug_dftb_pentacene' if args.basis == 'mio-1-1' else 'debug_dftb_pentacene_3ob')
    out_dir = os.path.join(_THIS_DIR, f'afm_fitted_{args.basis.replace("-", "_")}')

    _, atomPos = load_xyz(_XYZ_PATH)
    run_afm_with_fitted_params(src_dir, out_dir, A_pauli, beta_pauli, atomPos)

if __name__ == "__main__":
    main()
