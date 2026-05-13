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
from pyBall.globals import debug_plot_enabled
from pyBall.DFTB.TestUtils import CheckpointManager


'''
HOW TO USE IT:

python3 test_full_pipeline.py pentacene.xyz --output_dir YOUR_OUTPUT_FOLDER --step 0.1 --margin 2.0 --z_extra 2.0 --scan_range 3.0 --scan_step 0.1 --height_range 2.8 3.6 --height_step 0.1 --plot_steps


# Simple STM on TBTAP
python test_full_pipeline.py TBTAP.xyz --use_dense_projection --compute_stm

# Bond-resolved STM (AFM tip displacement affects STM)
python test_full_pipeline.py TBTAP.xyz --use_dense_projection --compute_stm --stm_bond_resolved

# Custom STM parameters
python test_full_pipeline.py TBTAP.xyz --use_dense_projection --compute_stm --stm_exp_beta 2.0 --stm_lumo_offsets 1 2


python test_full_pipeline.py TBTAP.xyz --basis 3ob-3-1 --use_dense_projection --max_shells 3 --compute_stm --stm_exp_beta 1.0 --stm_lumo_offsets 1


python test_full_pipeline.py /home/prokop/git/FireCore/cpp/common_resources/xyz/TBTAP.xyz --basis 3ob-3-1 --use_dense_projection --max_shells 3 --compute_stm --stm_bond_resolved --stm_exp_beta 1.0 --stm_lumo_offsets 1
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
    parser.add_argument('--plot_tip_disp', action='store_true', default=False, help='Plot tip displacement during relaxation')
    parser.add_argument('--use_dense_projection', action='store_true', default=False, help='Use dense matrix projection (supports d-orbitals, faster)')
    parser.add_argument('--max_shells', type=int, default=None, help='Max angular momentum shells (2=sp, 3=spd); auto-detected if not set')
    # STM arguments
    parser.add_argument('--compute_stm', action='store_true', default=False, help='Compute STM signal')
    parser.add_argument('--stm_lumo_offsets', type=int, nargs='+', default=[1, 2, 3], help='LUMO offsets from HOMO (default: 1 2 3)')
    parser.add_argument('--stm_use_exp_basis', action='store_true', default=True, help='Use exponential radial decay for STM')
    parser.add_argument('--stm_exp_beta', type=float, default=1.0, help='STM exponential decay constant (Å^-1)')
    parser.add_argument('--stm_exp_r0', type=float, default=3.0, help='STM reference distance (Å)')
    parser.add_argument('--stm_bond_resolved', action='store_true', default=False, help='Compute bond-resolved STM (STM at displaced tip positions)')
    # Parity checkpoint system (for optimization validation)
    parser.add_argument('--save-checkpoints', action='store_true', default=False, help='Save parity checkpoints to --checkpoint-dir')
    parser.add_argument('--verify-checkpoints', action='store_true', default=False, help='Verify parity against checkpoints in --checkpoint-dir')
    parser.add_argument('--checkpoint-dir', type=str, default=None, help='Directory for reference checkpoints')
    args = parser.parse_args()

    def _load_xyz_for_types(xyz_path: str):
        # Keep element mapping aligned with pyBall/OCL/AFM_utils.py:run_afm_from_xyz().
        ELEM_Z = {'H': 1, 'C': 6, 'N': 7, 'O': 8, 'P': 15, 'S': 16, 'Br': 35, 'I': 53}
        with open(xyz_path, "r") as f:
            lines = f.readlines()
        natoms = int(lines[0])
        pos = []
        types = []
        for line in lines[2:2 + natoms]:
            p = line.split()
            types.append(ELEM_Z.get(p[0], 6))
            pos.append([float(p[1]), float(p[2]), float(p[3])])
        return np.array(pos, dtype=np.float64), np.array(types, dtype=np.int32)

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
    
    # Prepare STM parameters
    stm_params = None
    if args.compute_stm:
        stm_params = {
            'compute': True,
            'lumo_offsets': args.stm_lumo_offsets,
            'use_exp_basis': args.stm_use_exp_basis,
            'exp_beta': args.stm_exp_beta,
            'exp_r0': args.stm_exp_r0,
            'bond_resolved': args.stm_bond_resolved
        }

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
        plot_steps=args.plot_steps,
        use_dense_projection=args.use_dense_projection,
        max_shells=args.max_shells,
        stm_params=stm_params
    )

    print(f"\ndf shape: {results['df'].shape}")
    print(f"df range: [{results['df'].min():.4f}, {results['df'].max():.4f}]")

    # Plot diagnostic panel with field components
    inter = results['intermediates']
    scan_xs = results['scan_xs']
    scan_ys = results['scan_ys']
    heights = results['heights']

    atomPos, atomTypes = _load_xyz_for_types(args.xyz_file)
    if args.plot_steps and debug_plot_enabled(2):
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
    
    # Plot tip displacement
    if args.plot_tip_disp and 'tip_disp' in inter and args.plot_steps and debug_plot_enabled(2):
        afm_utils.plot_tip_displacement(
            inter['tip_disp'],
            scan_xs,
            scan_ys,
            heights,
            args.output_dir
        )

    # Plot STM results
    if args.compute_stm and 'stm_grid' in inter and args.plot_steps and debug_plot_enabled(1):
        print(f"\nSTM range: [{inter['stm_grid'].min():.4e}, {inter['stm_grid'].max():.4e}]")
        afm_utils.plot_stm(
            inter['stm_grid'],
            scan_xs,
            scan_ys,
            heights,
            args.output_dir,
            prefix='stm_bond_resolved' if args.stm_bond_resolved else 'stm'
        )

    # Optional: save/verify parity checkpoints
    if args.save_checkpoints or args.verify_checkpoints:
        if not args.checkpoint_dir:
            raise ValueError("--checkpoint-dir is required when using --save-checkpoints/--verify-checkpoints")

        # Checkpoint loading requires intermediate .npy files to be saved
        from pyBall import globals
        if globals.DEBUG_SAVE_LEVEL < 2:
            raise ValueError(
                "Checkpoints require AFM_DEBUG_SAVE_LEVEL >= 2 to save intermediate .npy files. "
                "Set AFM_DEBUG_SAVE_LEVEL=2 in environment."
            )

        # Use the output_dir itself as base input for arrays saved during the run.
        # The library already writes the required auxiliary npy files into args.output_dir.
        checkpoint_mgr = CheckpointManager(args.checkpoint_dir, 'afm_pipeline_ref')

        grid_spec = results['grid_spec']
        base_config = {
            'xyz_file': args.xyz_file,
            'step': args.step,
            'origin': grid_spec['origin'].tolist(),
            'ngrid': grid_spec['ngrid'].tolist(),
            'pauli_params': {'A': pauli_A, 'beta': pauli_beta} if pauli_A is not None else None,
            'vdw_params': {'C6_CO': args.vdw_C6},
            'relax_params': {'K_LAT': args.relax_K},
        }

        # Load arrays produced by run_afm_pipeline (guarded by DEBUG_SAVE_LEVEL).
        # For checkpoint runs we expect DEBUG_SAVE_LEVEL >= 2.
        overlap_raw = np.load(os.path.join(args.output_dir, 'overlap_raw.npy'))
        # V_ES is an input from DFTB, not saved by AFM pipeline - get from intermediates
        V_ES = inter['V_ES']
        E_pauli_field = np.load(os.path.join(args.output_dir, 'E_Pauli_field.npy'))
        grads_pauli = np.load(os.path.join(args.output_dir, 'grads_E_Pauli.npy'))
        E_ES_field = np.load(os.path.join(args.output_dir, 'E_ES_field.npy'))
        grads_ES = np.load(os.path.join(args.output_dir, 'grads_E_ES.npy'))
        E_vdw = np.load(os.path.join(args.output_dir, 'E_vdw_field.npy'))
        grads_vdw = np.load(os.path.join(args.output_dir, 'grads_E_vdw.npy'))

        df = np.load(os.path.join(args.output_dir, 'df.npy'))
        tip_disp_dx = np.load(os.path.join(args.output_dir, 'tip_disp_dx.npy'))
        tip_disp_dy = np.load(os.path.join(args.output_dir, 'tip_disp_dy.npy'))

        checkpoints_to_handle = [
            ('pauli_overlap', {'overlap_raw': overlap_raw}, {'tip_rolled': True}),
            ('pauli_field', {'E_pauli_field': E_pauli_field, 'grads_pauli': grads_pauli}, {'step': args.step}),
            ('es_field', {'E_ES_field': E_ES_field, 'grads_ES': grads_ES, 'V_ES': V_ES}, {'step': args.step, 'tip_rolled': True}),
            ('dispersion_field', {
                'E_vdw': E_vdw,
                'grads_vdw': grads_vdw,
                'atomPos': atomPos,
                'atomTypes': atomTypes,
            }, {'C6_CO': args.vdw_C6, 'RA': 1.5}),
            ('relaxation', {
                'df': df,
                'tip_disp_dx': tip_disp_dx,
                'tip_disp_dy': tip_disp_dy,
            }, {'K_LAT': args.relax_K, 'N_RELAX': 50, 'step': args.step}),
        ]

        def _save_checkpoint(name, data):
            if not checkpoint_mgr.exists(name):
                checkpoint_mgr.save(name, data, config=base_config, overwrite=False)

        def _verify_checkpoint(name, data, tolerances):
            result = checkpoint_mgr.compare(
                name, data, config=base_config, tolerances=tolerances, strict_config=True, verbose=True
            )
            if not result['match']:
                raise AssertionError(f"Checkpoint verification failed for '{name}': {result['details']}")

        if args.save_checkpoints:
            for name, data, _conf in checkpoints_to_handle:
                _save_checkpoint(name, data)

        if args.verify_checkpoints:
            # Use generous tolerances for float32 CPU/GPU parity.
            tolerances = {
                'default': 1e-3,
            }
            for name, data, _conf in checkpoints_to_handle:
                _verify_checkpoint(name, data, tolerances=tolerances)

if __name__ == "__main__":
    main()
