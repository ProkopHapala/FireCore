#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path:
    sys.path.append(ROOT)

from pyBall.OCL.GridFFRelaxedScan import GridFFRelaxedScan, ensure_dir


def parse_type_map(txt):
    if not txt:
        return {}
    out = {}
    for part in txt.split(','):
        if '=' not in part:
            continue
        k, v = part.split('=', 1)
        out[k.strip()] = v.strip()
    return out


def parse_plq(txt):
    vals = [float(v.strip()) for v in txt.split(',') if v.strip()]
    if len(vals) != 3:
        raise ValueError(f'PLQ must have 3 comma-separated values, got: {txt}')
    return tuple(vals)


def default_paths():
    return {
        'mol': os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', 'PTCDA.xyz'),
        'sub': os.path.join(ROOT, 'cpp', 'common_resources', 'Substrates', 'CaF2_6L_Ni3.xyz'),
    }


def main():
    dfl = default_paths()
    ap = argparse.ArgumentParser(description='Headless relaxed GridFF scan of PTCDA on CaF2 with one oxygen constrained by harmonic tether.')
    ap.add_argument('--mol', default=dfl['mol'])
    ap.add_argument('--sub', default=dfl['sub'])
    ap.add_argument('--gridff', default=None, help='Path to existing GridFF Bspline_PLQd*.npy (overrides --gridff_name/--grid_sigma)')
    ap.add_argument('--outdir', default='out_ptcda_caf2_relaxed_scan')
    ap.add_argument('--mol_types', default='C=C_R,O=O_2,H=H')
    ap.add_argument('--grid_p0', type=float, nargs=3, default=[0.0, 0.0, 0.0])
    ap.add_argument('--grid_step', type=float, nargs=3, default=[23.175/225.0, 20.070/200.0, 48.472/382.0], help='GridFF step (dx, dy, dz)')
    ap.add_argument('--grid_alpha', type=float, default=1.5, help='GridFF Morse alpha')
    ap.add_argument('--grid_sigma', type=float, default=0.0, help='GridFF Gaussian sigma (selects cached file if --gridff not specified)')
    ap.add_argument('--gridff_name', default='CaF2_6L_Ni3_rect_nx2_nz1_L2_top', help='Name of substrate cache folder under tests/tMMFF/data/')
    ap.add_argument('--dt', type=float, default=0.01, help='Relaxation dt')
    ap.add_argument('--damp', type=float, default=0.01, help='Relaxation damping')
    ap.add_argument('--Fconv', type=float, default=1e-4, help='Force convergence threshold')
    ap.add_argument('--nstep_max', type=int, default=8000, help='Max relaxation steps per scan point')
    ap.add_argument('--out_stride', type=int, default=10, help='Save stride during relaxation')
    ap.add_argument('--anchor_k', type=float, default=2000.0, help='Anchor spring stiffness')
    ap.add_argument('--anchor_z_above', type=float, default=6.0, help='Initial anchor height above substrate')
    ap.add_argument('--warm_start', action='store_true', help='Warm start each scan point from previous relaxed geometry')
    ap.add_argument('--plq_plot', type=float, nargs=3, default=[1.0, 1.0, 0.0], help='PLQ coefficients for XZ GridFF plot')
    ap.add_argument('--test_particle', type=str, default='H', choices=['H', 'O', 'C'], help='Test particle type for GridFF plotting (H, O, C)')
    ap.add_argument('--test_charge', type=float, default=0.0, help='Test particle charge override (e.g., +0.2 for H, -0.4 for O)')
    ap.add_argument('--plot_xy_z', type=float, default=3.0, help='Z height for XY top view plot [A]')
    ap.add_argument('--plot_vmax', type=float, default=0.1, help='Colorbar range for GridFF plots [eV] (symmetric +/-)')
    ap.add_argument('--nscan', type=int, default=101, help='Number of scan points')
    ap.add_argument('--x0', type=float, default=2.0, help='Scan start x')
    ap.add_argument('--x1', type=float, default=12.0, help='Scan end x (0.1A steps for nscan=101)')
    ap.add_argument('--y0', type=float, default=2.0, help='Scan y coordinate')
    ap.add_argument('--z0', type=float, default=None, help='Scan z coordinate (default from anchor)')
    ap.add_argument('--cold_start', type=int, default=0, help='If nonzero, do not warm-start between scan points')
    ap.add_argument('--verbose', action='store_true')
    args = ap.parse_args()

    outdir = os.path.abspath(args.outdir)
    ensure_dir(outdir)

    # Always reuse existing GridFF; never regenerate unless explicitly deleted
    if args.gridff is not None:
        expected_gridff = os.path.abspath(args.gridff)
    else:
        sigma_tag = f"{float(args.grid_sigma):.3f}".replace('.', 'p')
        expected_gridff = os.path.join(ROOT, 'tests', 'tMMFF', 'data', args.gridff_name, f'Bspline_PLQd_sigma_{sigma_tag}.npy')
    if not os.path.exists(expected_gridff):
        print(f'ERROR: GridFF cache not found at {expected_gridff}')
        print('Generate it once (and only once) by running:')
        print('  cd tests/tMMFF')
        print('  python3 run_test_GridFF_CaF2.py --src_xyz ../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx2_nz1_L2_top.xyz --save_name CaF2_6L_Ni3_rect_nx2_nz1_L2_top --job PLQ --Q 1.0')
        print('Then re-run this script; it will reuse the cached .npy and never regenerate.')
        return 1
    scanner = GridFFRelaxedScan(
        mol_path=args.mol,
        sub_xyz_path=args.sub,
        gridff_path=expected_gridff,
        out_dir=outdir,
        mol_type_map={'C': 'C_R', 'O': 'O_2', 'H': 'H'},
        grid_alpha=args.grid_alpha,
        grid_sigma=args.grid_sigma,
        grid_step=args.grid_step,
        dt=args.dt,
        damp=args.damp,
        Fconv=args.Fconv,
        nstep_max=args.nstep_max,
        out_stride=args.out_stride,
        anchor_k=args.anchor_k,
    )
    scanner.prepare(anchor_z_above=args.anchor_z_above, lateral_shift=(args.x0, args.y0), generate_grid=False)

    z0 = float(args.z0) if (args.z0 is not None) else float(scanner.get_anchor_position()[2])
    anchor_path = scanner.make_linear_anchor_path(args.x0, args.x1, args.y0, z0, args.nscan)

    xyz_path = os.path.join(outdir, 'ptcda_caf2_relaxed_scan.xyz')
    npz_path = os.path.join(outdir, 'ptcda_caf2_relaxed_scan.npz')
    # Derive PLQ from test particle type and charge
    def get_test_particle_plq(particle_type, charge_override):
        """Get PLQ coefficients for test particle based on atom type and charge"""
        # Actual REQ parameters from AtomTypes.dat
        atom_reqs = {
            'H': {'R': 1.4430, 'E': np.sqrt(0.00190802059), 'Q': 0.0},
            'O': {'R': 1.7500, 'E': np.sqrt(0.00260184625), 'Q': -0.4},
            'C': {'R': 1.9255, 'E': np.sqrt(0.00455323095), 'Q': 0.0}
        }
        
        req = atom_reqs.get(particle_type, atom_reqs['H'])
        if charge_override != 0.0:
            req['Q'] = charge_override
        
        # Proper Morse-derived PLQ coefficients
        # GridFF stores: V_P = Σ_i exp(-2α*(r-Ri)), V_L = Σ_i exp(-α*(r-Ri))  
        # Test particle provides: P = exp(2α*Rj), L = 2E * exp(α*Rj), Q = qj
        alpha = 1.5  # Default alphaMorse from GridFF generation
        P = np.exp(2 * alpha * req['R'])      # Pauli prefactor
        L = 2 * req['E'] * np.exp(alpha * req['R'])  # London prefactor  
        Q = req['Q']                           # Coulomb charge
        
        return (P, L, Q), req
    
    plq_plot, test_req = get_test_particle_plq(args.test_particle, args.test_charge)
    
    xz_png = os.path.join(outdir, 'ptcda_caf2_gridff_xz.png')
    xy_png = os.path.join(outdir, 'ptcda_caf2_gridff_xy.png')
    
    # Plot XZ slice with test particle parameters
    scanner.plot_gridff_xz(xz_png, iy=-1, plq=plq_plot, vmax=args.plot_vmax, 
                          test_req=test_req, test_particle=args.test_particle)
    
    # Plot XY slice at specified height
    scanner.plot_gridff_xy(xy_png, z=args.plot_xy_z, plq=plq_plot, vmax=args.plot_vmax,
                          test_req=test_req, test_particle=args.test_particle)
    opp_png = os.path.join(outdir, 'ptcda_caf2_opposite_oxygen_path.png')

    scan_result = scanner.run_path_scan(anchor_path, xyz_path=xyz_path, warm_start=(args.cold_start == 0), verbose=args.verbose)
    scanner.save_scan_npz(scan_result, npz_path)
    scanner.plot_opposite_trajectory(scan_result, opp_png)

    print('DONE')
    print(f'xyz      : {xyz_path}')
    print(f'npz      : {npz_path}')
    print(f'gridff xz: {xz_png}')
    print(f'gridff xy: {xy_png}')
    print(f'opposite : {opp_png}')


if __name__ == '__main__':
    main()
