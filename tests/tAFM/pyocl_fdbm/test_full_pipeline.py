#!/usr/bin/env python3
"""
Test full AFM pipeline from .xyz to AFM images using AFM_utils.py
"""

import sys, os, argparse
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL import AFM_utils as afm_utils
from pyBall.OCL import AFM as afm_mod
from pyBall.DFTB.DFTBplusParser import parse_basis_hsd_ang
from pyBall import dftb_utils as du
import numpy as np


'''
HOW TO USE IT:

python3 test_full_pipeline.py pentacene.xyz --output_dir YOUR_OUTPUT_FOLDER --step 0.1 --margin 2.0 --z_extra 2.0 --scan_range 3.0 --scan_step 0.1 --height_range 2.8 3.6 --height_step 0.1 --plot_steps
'''

def main():
    parser = argparse.ArgumentParser(description='Full AFM pipeline from .xyz to AFM images')
    parser.add_argument('xyz_file', type=str, help='Input .xyz file')
    parser.add_argument('--output_dir', type=str, default='test_pipeline_output', help='Output directory')
    parser.add_argument('--basis', type=str, default='mio-1-1', help='DFTB basis set')
    parser.add_argument('--co_tip_dir', type=str, default=None, help='CO tip directory (if not provided, computed on-the-fly)')
    parser.add_argument('--step', type=float, default=0.1, help='Density grid step (A)')
    parser.add_argument('--margin', type=float, default=4.0, help='Grid margin (A)')
    parser.add_argument('--z_extra', type=float, default=6.0, help='Extra z-space (A)')
    parser.add_argument('--scan_range', type=float, default=3.0, help='Scan range around molecule (A)')
    parser.add_argument('--scan_step', type=float, default=0.1, help='Scan step (A)')
    parser.add_argument('--height_range', type=float, nargs=2, default=[2.8, 3.6], help='Height range (A)')
    parser.add_argument('--height_step', type=float, default=0.1, help='Height step (A)')
    parser.add_argument('--pauli_A', type=float, default=None, help='Pauli A parameter (default: use PAULI_FITTED_DEFAULTS for basis)')
    parser.add_argument('--pauli_beta', type=float, default=None, help='Pauli beta parameter (default: use PAULI_FITTED_DEFAULTS for basis)')
    parser.add_argument('--vdw_C6', type=float, default=30.0, help='vdW C6 parameter')
    parser.add_argument('--relax_K', type=float, default=0.5, help='Relaxation K_LAT parameter')
    parser.add_argument('--plot_steps', action='store_true', default=True, help='Plot intermediate steps')
    parser.add_argument('--no_plot_steps', action='store_false', dest='plot_steps', help='Skip intermediate plots')
    parser.add_argument('--fit_pauli', action='store_true', help='Fit Pauli parameters using DFTB z-scan reference')
    parser.add_argument('--fit_target_indices', type=int, nargs='+', default=[0], help='Target atom indices for fitting (default: 0)')
    parser.add_argument('--fit_z_min', type=float, default=2.0, help='Fit range minimum z (Å)')
    parser.add_argument('--fit_z_max', type=float, default=3.0, help='Fit range maximum z (Å)')
    parser.add_argument('--fit_zscan_dir', type=str, default=None, help='Pre-computed DFTB z-scan directory')
    args = parser.parse_args()

    SK = args.basis
    slako_prefix = du.SK_PATHS[SK]

    # Make work_dir specific to molecule and basis
    mol_name = os.path.splitext(os.path.basename(args.xyz_file))[0]
    work_dir = os.path.join(args.output_dir, f'dftb_work_{mol_name}_{SK}')

    # Resolve Pauli params: CLI > PAULI_FITTED_DEFAULTS[basis] > error
    pauli_A    = args.pauli_A
    pauli_beta = args.pauli_beta
    if pauli_A is None or pauli_beta is None:
        defaults = afm_mod.PAULI_FITTED_DEFAULTS.get(SK)
        if defaults is None or defaults.get('A') is None:
            if not args.fit_pauli:
                raise ValueError(f"No fitted Pauli defaults for basis '{SK}'. Provide --pauli_A and --pauli_beta, or use --fit_pauli.")
        else:
            pauli_A    = pauli_A    if pauli_A    is not None else defaults['A']
            pauli_beta = pauli_beta if pauli_beta is not None else defaults['beta']
            print(f"Using PAULI_FITTED_DEFAULTS[{SK}]: A={pauli_A}, beta={pauli_beta}")
    
    # Prepare fit_pauli_params if fitting is requested
    fit_pauli_params = None
    if args.fit_pauli:
        if args.fit_zscan_dir is None:
            raise ValueError("--fit_zscan_dir is required when using --fit_pauli")
        fit_pauli_params = {
            'zscan_dir': args.fit_zscan_dir,
            'target_indices': args.fit_target_indices,
            'z_min': args.fit_z_min,
            'z_max': args.fit_z_max,
            'basis': SK
        }
        print(f"Will fit Pauli parameters using DFTB reference from {args.fit_zscan_dir}")
    
    results = afm_utils.run_afm_from_xyz(
        xyz_file=args.xyz_file,
        output_dir=args.output_dir,
        basis=SK,
        slako_prefix=slako_prefix,
        co_tip_dir=args.co_tip_dir,
        work_dir=work_dir,
        step=args.step, margin=args.margin, z_extra=args.z_extra,
        scan_range=args.scan_range, scan_step=args.scan_step,
        height_range=tuple(args.height_range), height_step=args.height_step,
        pauli_params={'A': pauli_A, 'beta': pauli_beta} if pauli_A is not None else None,
        fit_pauli=args.fit_pauli,
        fit_pauli_params=fit_pauli_params,
        vdw_params={'C6_CO': args.vdw_C6},
        relax_params={'K_LAT': args.relax_K},
        plot_steps=args.plot_steps
    )

    print(f"\ndf shape: {results['df'].shape}")
    print(f"df range: [{results['df'].min():.4f}, {results['df'].max():.4f}]")

    # Plot diagnostic panel with field components
    inter = results['intermediates']
    heights = np.arange(args.height_range[0], args.height_range[1], args.height_step)
    afm_utils.plot_diagnostic_panel(
        inter['E_pauli_field'], inter['E_ES_field'], inter['E_vdw'],
        inter['E_pauli_field'] + inter['E_ES_field'] + inter['E_vdw'],
        results['grid_spec']['origin'], args.step, heights,
        args.output_dir
    )
    afm_utils.plot_diagnostic_slices(
        inter['E_pauli_field'], inter['E_ES_field'], inter['E_vdw'],
        results['grid_spec']['origin'], args.step, heights,
        args.output_dir
    )

if __name__ == "__main__":
    main()
