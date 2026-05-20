#!/usr/bin/env python3
"""
Test Modular AFM/STM pipeline with disk-based caching.
"""

import sys, os, argparse
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

import numpy as np
import matplotlib.pyplot as plt
from pyBall.OCL import ModularPipeline as mp_mod
from pyBall import dftb_utils as du
from pyBall.OCL import AFM as afm_mod

def main():
    parser = argparse.ArgumentParser(description='Modular AFM/STM pipeline with disk-based caching')
    parser.add_argument('xyz_file', type=str, help='Input .xyz file')
    parser.add_argument('--output_dir', type=str, default='test_modular_output', help='Output directory')
    parser.add_argument('--basis', type=str, default='mio-1-1', help='DFTB basis set')
    parser.add_argument('--step', type=float, default=0.1, help='Density grid step (A)')
    parser.add_argument('--margin', type=float, default=4.0, help='Grid margin (A)')
    parser.add_argument('--z_extra', type=float, default=6.0, help='Extra z-space (A)')
    parser.add_argument('--scan_range', type=float, default=3.0, help='Scan range around molecule (A)')
    parser.add_argument('--scan_step', type=float, default=0.1, help='Scan step (A)')
    parser.add_argument('--height_range', type=float, nargs=2, default=[2.5, 3.3], help='Height range (A)')
    parser.add_argument('--height_step', type=float, default=0.1, help='Height step (A)')
    parser.add_argument('--pauli_A', type=float, default=None, help='Pauli A parameter')
    parser.add_argument('--pauli_beta', type=float, default=None, help='Pauli beta parameter')
    parser.add_argument('--vdw_C6', type=float, default=30.0, help='vdW C6 parameter')
    parser.add_argument('--relax_K', type=float, default=0.5, help='Relaxation K_LAT parameter')
    parser.add_argument('--ppm_mode', action='store_true', default=True, help='Use PPM mode')
    parser.add_argument('--no_ppm_mode', action='store_false', dest='ppm_mode', help='Disable PPM mode')
    parser.add_argument('--co_tip_dir', type=str, default=None, help='CO tip directory (if not provided, uses global cache or computes on-the-fly)')
    
    # Force recompute flags for modularity
    parser.add_argument('--force_scf', action='store_true', help='Force SCF (Stage 1) recompute')
    parser.add_argument('--force_project', action='store_true', help='Force projection (Stage 2) recompute')
    parser.add_argument('--force_potentials', action='store_true', help='Force potentials (Stage 3) recompute')
    parser.add_argument('--force_relax', action='store_true', help='Force relaxation (Stage 4) recompute')
    parser.add_argument('--force_all', action='store_true', help='Force all stages recompute')
    
    # STM options
    parser.add_argument('--compute_stm', action='store_true', default=False, help='Compute standard STM')
    parser.add_argument('--compute_br_stm', action='store_true', default=False, help='Compute Bond-Resolved STM')
    parser.add_argument('--stm_lumo_offsets', type=int, nargs='+', default=[1, 2, 3], help='LUMO offsets from HOMO')
    parser.add_argument('--stm_use_exp_basis', action='store_true', default=True, help='Use exponential basis for STM')
    parser.add_argument('--stm_exp_beta', type=float, default=1.0, help='STM exponential decay factor')
    parser.add_argument('--stm_exp_r0', type=float, default=3.0, help='STM reference distance r0')
    
    args = parser.parse_args()
    
    SK = args.basis
    slako_prefix = du.SK_PATHS[SK]
    
    # Resolve Pauli parameters
    pauli_A = args.pauli_A
    pauli_beta = args.pauli_beta
    if pauli_A is None or pauli_beta is None:
        defaults = afm_mod.PAULI_FITTED_DEFAULTS.get(SK)
        if defaults is not None:
            pauli_A = pauli_A if pauli_A is not None else defaults['A']
            pauli_beta = pauli_beta if pauli_beta is not None else defaults['beta']
        else:
            pauli_A = 787.22
            pauli_beta = 1.2371
            
    pauli_params = {'A': pauli_A, 'beta': pauli_beta}
    vdw_params = {'C6_CO': args.vdw_C6}
    relax_params = {'K_LAT': args.relax_K}
    
    # Initialize modular pipeline
    pipeline = mp_mod.ModularAFMPipeline(
        xyz_file=args.xyz_file,
        output_dir=args.output_dir,
        basis=SK,
        slako_prefix=slako_prefix,
        step=args.step,
        margin=args.margin,
        z_extra=args.z_extra,
        scan_range=args.scan_range,
        scan_step=args.scan_step,
        height_range=args.height_range,
        height_step=args.height_step,
        co_tip_dir=args.co_tip_dir
    )
    
    # Determine recomputation needs
    f_scf = args.force_all or args.force_scf
    f_proj = args.force_all or args.force_project
    f_pots = args.force_all or args.force_potentials
    f_relax = args.force_all or args.force_relax
    
    # Execute stages
    import time
    
    t0 = time.time()
    dm_dense, eigvecs, eigvals = pipeline.stage1_scf(force_recompute=f_scf)
    t1 = time.time()
    print(f"Stage 1 SCF time: {t1 - t0:.2f} s")
    
    rho_scf, rho_na, rho_diff = pipeline.stage2_project(dm_dense, force_recompute=f_proj)
    t2 = time.time()
    print(f"Stage 2 Grid Projection time: {t2 - t1:.2f} s")
    
    V_ES, E_pauli, E_ES, E_vdw, F_total = pipeline.stage3_potentials(
        rho_scf, rho_na, rho_diff, force_recompute=f_pots,
        pauli_params=pauli_params, vdw_params=vdw_params
    )
    t3 = time.time()
    print(f"Stage 3 FDBM Potentials time: {t3 - t2:.2f} s")
    
    df, tip_disp, FEs_relax = pipeline.stage4_relax(
        F_total, force_recompute=f_relax,
        relax_params=relax_params, ppm_mode=args.ppm_mode
    )
    t4 = time.time()
    print(f"Stage 4 Tip Relaxation time: {t4 - t3:.2f} s")
    
    # Debug: Print CO tip displacement statistics per height slice
    print("\n--- CO Tip Displacement Statistics ---")
    for iz, h in enumerate(pipeline.heights):
        dx = tip_disp['dx'][:, :, iz]
        dy = tip_disp['dy'][:, :, iz]
        disp_mag = np.sqrt(dx**2 + dy**2)
        print(f"  z={h:.2f} A: max|dx|={np.abs(dx).max():.4f} A, max|dy|={np.abs(dy).max():.4f} A, max|disp|={disp_mag.max():.4f} A, mean|disp|={disp_mag.mean():.4f} A")
    
    print("\n--- Plots generation ---")
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Plot AFM frequency shift for all heights
    nz = len(pipeline.heights)
    ncols = min(4, nz)
    nrows = int(np.ceil(nz / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3, nrows*2.5))
    if nz == 1:
        axes = np.array([axes])
    axes = axes.ravel()
    for iz, h in enumerate(pipeline.heights):
        im = axes[iz].imshow(df[:, :, iz].T, origin='lower', extent=[pipeline.scan_xs[0], pipeline.scan_xs[-1], pipeline.scan_ys[0], pipeline.scan_ys[-1]])
        axes[iz].set_title(f'z={h:.2f} A')
        axes[iz].set_xlabel('x (A)')
        axes[iz].set_ylabel('y (A)')
        plt.colorbar(im, ax=axes[iz], label='df (Hz)')
    for iz in range(nz, len(axes)):
        axes[iz].axis('off')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'afm_df_all_heights.png'), dpi=150)
    plt.close()
    print(f"  Saved AFM plots (all heights) to {os.path.join(args.output_dir, 'afm_df_all_heights.png')}")
    
    if args.compute_stm:
        t_stm0 = time.time()
        stm_map = pipeline.stage5_stm(
            eigvecs, eigvals, lumo_offsets=args.stm_lumo_offsets,
            use_exp_basis=args.stm_use_exp_basis, exp_beta=args.stm_exp_beta, exp_r0=args.stm_exp_r0
        )
        print(f"Stage 5 STM time: {time.time() - t_stm0:.2f} s")
        
        # Plot STM for all heights
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3, nrows*2.5))
        if nz == 1:
            axes = np.array([axes])
        axes = axes.ravel()
        for iz, h in enumerate(pipeline.heights):
            im = axes[iz].imshow(stm_map[:, :, iz].T, origin='lower', extent=[pipeline.scan_xs[0], pipeline.scan_xs[-1], pipeline.scan_ys[0], pipeline.scan_ys[-1]], cmap='inferno')
            axes[iz].set_title(f'z={h:.2f} A')
            axes[iz].set_xlabel('x (A)')
            axes[iz].set_ylabel('y (A)')
            plt.colorbar(im, ax=axes[iz], label='LDOS (a.u.)')
        for iz in range(nz, len(axes)):
            axes[iz].axis('off')
        plt.tight_layout()
        plt.savefig(os.path.join(args.output_dir, 'stm_all_heights.png'), dpi=150)
        plt.close()
        print(f"  Saved Standard STM plots (all heights) to {os.path.join(args.output_dir, 'stm_all_heights.png')}")
        
    if args.compute_br_stm:
        t_br0 = time.time()
        br_stm_map = pipeline.stage6_br_stm(
            eigvecs, eigvals, tip_disp, lumo_offsets=args.stm_lumo_offsets,
            use_exp_basis=args.stm_use_exp_basis, exp_beta=args.stm_exp_beta, exp_r0=args.stm_exp_r0
        )
        print(f"Stage 6 BR-STM time: {time.time() - t_br0:.2f} s")
        
        # Plot BR-STM for all heights
        fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*3, nrows*2.5))
        if nz == 1:
            axes = np.array([axes])
        axes = axes.ravel()
        for iz, h in enumerate(pipeline.heights):
            im = axes[iz].imshow(br_stm_map[:, :, iz].T, origin='lower', extent=[pipeline.scan_xs[0], pipeline.scan_xs[-1], pipeline.scan_ys[0], pipeline.scan_ys[-1]], cmap='inferno')
            axes[iz].set_title(f'z={h:.2f} A')
            axes[iz].set_xlabel('x (A)')
            axes[iz].set_ylabel('y (A)')
            plt.colorbar(im, ax=axes[iz], label='LDOS (a.u.)')
        for iz in range(nz, len(axes)):
            axes[iz].axis('off')
        plt.tight_layout()
        plt.savefig(os.path.join(args.output_dir, 'stm_bond_resolved_all_heights.png'), dpi=150)
        plt.close()
        print(f"  Saved Bond-Resolved STM plots (all heights) to {os.path.join(args.output_dir, 'stm_bond_resolved_all_heights.png')}")
        
    print(f"\n[ModularPipeline] Completed pipeline run. Total elapsed: {time.time() - t0:.2f} s")

if __name__ == '__main__':
    main()
