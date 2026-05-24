#!/usr/bin/env python3
"""
Test script for electrostatics comparison: GridFF vs Ewald2D vs Brute Force

Tests electrostatic potential from three methods on NaCl_1x1_L3.xyz:
1. GridFF (OpenCL B-spline sampling)
2. Ewald2D (2D Fourier representation)
3. Brute force (direct Coulomb sum over periodic replicas)

Usage:
  python test_electrostatics_comparison.py --xyz NaCl_1x1_L3.xyz --grid Bspline_PLQd.npy
  python test_electrostatics_comparison.py --xyz NaCl_1x1_L3.xyz --grid Bspline_PLQd.npy --N_rep 30 --coarse_spacing 0.5

See: doc/Topics/OnSurfaceAssembly/Ewald_2D.md for theory
"""

import sys
import os
import argparse
import numpy as np

sys.path.insert(0, os.path.expanduser("~/git/FireCore"))

from pyBall.AtomicSystem import AtomicSystem
from pyBall.Ewald2D import Ewald2D
from pyBall.OCL.Surface_utils import compare_electrostatics_methods

# OpenCL Ewald (optional)
try:
    from pyBall.OCL.SurfaceEwald import SurfaceEwaldCL
    HAS_OPENCL_EWALD = True
except Exception as e:
    HAS_OPENCL_EWALD = False
    print(f"Warning: OpenCL Ewald not available: {e}")


def main():
    parser = argparse.ArgumentParser( description="Compare electrostatics: GridFF vs Ewald2D vs Brute Force")
    parser.add_argument("--xyz", type=str, required=True,help="Path to substrate .xyz file (must have lattice vectors and charges)")
    parser.add_argument("--grid", type=str, required=True,  help="Path to GridFF Bspline_PLQd.npy file")
    parser.add_argument("--n_harm", type=int, default=4, help="Ewald2D G-vector half-width (default: 4)")
    parser.add_argument("--N_rep", type=int, default=20, help="Brute force periodic replicas (default: 20)")
    parser.add_argument("--coarse_spacing", type=float, default=1.0,  help="Point spacing (Å) for 1D brute force (default: 1.0)")
    parser.add_argument("--grid_res", type=int, default=200,  help="Resolution for 2D slices (default: 200)")
    parser.add_argument("--outdir", type=str, default="results_electrostatics",  help="Output directory (default: results_electrostatics)")
    parser.add_argument("--prefix", type=str, default=None, help="Output file prefix (default: derived from xyz filename)")
    parser.add_argument("--test_opencl", action="store_true", help="Also test OpenCL Ewald implementation")
    parser.add_argument("--test_cl_only", action="store_true", help="Only test OpenCL Ewald (skip GridFF comparison)")

    args = parser.parse_args()

    # Derive prefix from xyz filename if not provided
    if args.prefix is None:
        args.prefix = os.path.splitext(os.path.basename(args.xyz))[0]

    print(f"\n{'='*60}")
    print(f"Electrostatics Comparison Test")
    print(f"{'='*60}")
    print(f"Substrate: {args.xyz}")
    print(f"GridFF:    {args.grid}")
    print(f"Output:    {args.outdir}/{args.prefix}_*.png")
    print(f"{'='*60}\n")

    # ======================================================================
    # Step 1: Load substrate and initialize Ewald2D
    # ======================================================================
    print("Step 1: Loading substrate and initializing Ewald2D...")
    sys_at = AtomicSystem(fname=args.xyz)
    assert sys_at.lvec is not None, f"No lattice vectors in {args.xyz}"
    assert sys_at.qs is not None, f"No charges in {args.xyz}"

    ew = Ewald2D.from_AtomicSystem(sys_at, n_harm=args.n_harm)
    ew.print_info()

    # ======================================================================
    # Step 2: Verify GridFF file exists and load metadata
    # ======================================================================
    print("\nStep 2: Verifying GridFF file...")

    from pyBall.OCL.Surface_utils import load_gridff_array, load_gridff_metadata
    grid_data = load_gridff_array(args.grid)
    print(f"  Grid shape: {grid_data.shape}")

    # Try to load metadata from JSON
    meta = load_gridff_metadata(args.grid)
    if meta is not None:
        print(f"  Loaded metadata from JSON")
        if 'g0' in meta:
            print(f"  Grid origin (g0): {meta['g0']}")
        if 'dg' in meta:
            print(f"  Grid spacing (dg): {meta['dg']}")
    else:
        print(f"  No metadata JSON found")

    # ======================================================================
    # Step 3: Run comparison
    # ======================================================================
    print("\nStep 3: Running electrostatics comparison...")

    report = compare_electrostatics_methods(
        sys_at=sys_at,
        gridff_path=args.grid,
        sub_xyz_path=args.xyz,
        ew=ew,
        save_dir=args.outdir,
        prefix=args.prefix,
        N_rep=args.N_rep,
        coarse_spacing=args.coarse_spacing,
        grid_res=args.grid_res,
        z_min_for_error=2.0,
        verbose=True
    )

    # ======================================================================
    # Summary
    # ======================================================================
    print(f"\n{'='*60}")
    print(f"Summary of Results")
    print(f"{'='*60}")

    print("\n1D Line Scan Errors (vs Brute Force Reference):")
    for name in ['z_on_Na', 'z_on_Cl', 'z_midpoint']:
        if name in report['stats']:
            s = report['stats'][name]
            print(f"  {name}:")
            print(f"    Ewald2D: RMSE={s['ewald_vs_brute']['rmse']:.4e}, max={s['ewald_vs_brute']['max_err']:.4e}")
            print(f"    GridFF:  RMSE={s['gridff_vs_brute']['rmse']:.4e}, max={s['gridff_vs_brute']['max_err']:.4e}")

    print("\n2D Slice Errors (GridFF vs Ewald2D):")
    for key in report['stats']:
        if key.startswith('xy_z'):
            s = report['stats'][key]['gridff_vs_ewald']
            print(f"  {key}: RMSE={s['rmse']:.4e}, max={s['max_err']:.4e}")
    if 'xz' in report['stats']:
        s = report['stats']['xz']['gridff_vs_ewald']
        print(f"  xz_slice: RMSE={s['rmse']:.4e}, max={s['max_err']:.4e}")

    # ======================================================================
    # Step 4: Optional OpenCL Ewald test
    # ======================================================================
    if args.test_opencl or args.test_cl_only:
        print(f"\n{'='*60}")
        print(f"Step 4: OpenCL Ewald Test")
        print(f"{'='*60}")

        if not HAS_OPENCL_EWALD:
            print("ERROR: OpenCL Ewald not available. Install pyopencl.")
        else:
            try:
                # Initialize OpenCL Ewald
                print("\nInitializing OpenCL Ewald...")
                ew_cl = SurfaceEwaldCL()

                # Prepare system
                ion_data = np.column_stack([
                    sys_at.apos[:, 0],  # rx
                    sys_at.apos[:, 1],  # ry
                    sys_at.apos[:, 2],  # rz
                    sys_at.qs           # q
                ])
                a_vec = sys_at.lvec[0, :2]
                b_vec = sys_at.lvec[1, :2]

                ew_cl.prepare_system(ion_data, a_vec, b_vec, n_harm=args.n_harm)

                # Test 1: Vacuum evaluation
                print("\nOpenCL Vacuum Evaluation Test:")
                import matplotlib.pyplot as plt

                # Create test grid
                xv = np.linspace(0, 4.0, 50)
                yv = np.linspace(0, 4.0, 50)
                X_test, Y_test = np.meshgrid(xv, yv)
                z_test = 2.0

                # OpenCL evaluation
                phi_cl = ew_cl.eval_vacuum(X_test, Y_test, z_test)

                # Python reference
                phi_py = ew.phi_vacuum_xy(X_test, Y_test, z_test)

                # Compare
                diff = phi_cl - phi_py
                rmse = np.sqrt(np.mean(diff**2))
                max_err = np.max(np.abs(diff))

                print(f"  RMSE: {rmse:.6e} eV")
                print(f"  Max error: {max_err:.6e} eV")

                if rmse < 1e-5:
                    print("  ✓ PASS: OpenCL Ewald agrees with Python")
                else:
                    print("  ✗ FAIL: OpenCL Ewald differs from Python")

                # Test 2: Full evaluation (1D line)
                print("\nOpenCL Full Evaluation Test (1D line):")
                x0, y0 = 0.5, 0.5
                z_arr = np.linspace(-0.5, 5.0, 100)

                X_line = np.full((1, len(z_arr)), x0, dtype=np.float32)
                Y_line = np.full((1, len(z_arr)), y0, dtype=np.float32)
                Z_line = z_arr.reshape(1, -1).astype(np.float32)

                phi_cl_line = ew_cl.eval_full(X_line, Y_line, Z_line)[0, :]
                phi_py_line = ew.phi_full_1d(x0, y0, z_arr)

                diff_line = phi_cl_line - phi_py_line
                rmse_line = np.sqrt(np.mean(diff_line**2))
                max_err_line = np.max(np.abs(diff_line))

                print(f"  RMSE: {rmse_line:.6e} eV")
                print(f"  Max error: {max_err_line:.6e} eV")

                if rmse_line < 1e-5:
                    print("  ✓ PASS: OpenCL Full agrees with Python")
                else:
                    print("  ✗ FAIL: OpenCL Full differs from Python")

                # Generate comparison plots
                print("\nGenerating OpenCL comparison plots...")
                import matplotlib.pyplot as plt

                # Figure 5: XY slice comparison (Python vs OpenCL)
                fig5, axes5 = plt.subplots(2, 3, figsize=(15, 8))

                for idx, (z_h, phi_cl, phi_py) in enumerate(zip([z_test], [phi_cl], [phi_py])):
                    ax_py = axes5[0, idx]
                    ax_cl = axes5[1, idx]

                    vmin = min(phi_py.min(), phi_cl.min())
                    vmax = max(phi_py.max(), phi_cl.max())

                    im_py = ax_py.pcolormesh(X_test, Y_test, phi_py, shading='auto',
                                              cmap='RdBu_r', vmin=vmin, vmax=vmax)
                    ax_py.set_title(f'Python Ewald z={z_h:.2f}Å')
                    ax_py.set_xlabel('x (Å)')
                    ax_py.set_ylabel('y (Å)')
                    plt.colorbar(im_py, ax=ax_py, label='φ (eV)')

                    im_cl = ax_cl.pcolormesh(X_test, Y_test, phi_cl, shading='auto',
                                              cmap='RdBu_r', vmin=vmin, vmax=vmax)
                    ax_cl.set_title(f'OpenCL Ewald z={z_h:.2f}Å')
                    ax_cl.set_xlabel('x (Å)')
                    ax_cl.set_ylabel('y (Å)')
                    plt.colorbar(im_cl, ax=ax_cl, label='φ (eV)')

                # Hide unused subplots
                for idx in range(1, 3):
                    axes5[0, idx].axis('off')
                    axes5[1, idx].axis('off')

                fig5.suptitle(f'Python vs OpenCL Ewald (Vacuum)\nRMSE: {rmse:.2e} eV, Max: {max_err:.2e} eV')
                fig5.tight_layout()
                fig5_path = os.path.join(args.outdir, f'{args.prefix}_fig5_opencl_vacuum.png')
                fig5.savefig(fig5_path, dpi=150)
                print(f"  Saved {fig5_path}")

                # Figure 6: 1D line comparison (Python vs OpenCL)
                fig6, axes6 = plt.subplots(2, 1, figsize=(12, 8))

                ax1 = axes6[0]
                ax1.plot(z_arr, phi_py_line, 'b-', label='Python Ewald', linewidth=2)
                ax1.plot(z_arr, phi_cl_line, 'r--', label='OpenCL Ewald', linewidth=1.5)
                ax1.set_xlabel('z (Å)')
                ax1.set_ylabel('φ (eV)')
                ax1.set_title(f'1D Line Scan at ({x0:.2f}, {y0:.2f})')
                ax1.legend()
                ax1.grid(True, alpha=0.3)

                ax2 = axes6[1]
                ax2.plot(z_arr, phi_cl_line - phi_py_line, 'g-', label='Difference (OpenCL - Python)')
                ax2.axhline(y=0, color='k', linestyle='--', alpha=0.3)
                ax2.set_xlabel('z (Å)')
                ax2.set_ylabel('Δφ (eV)')
                ax2.set_title(f'Difference (RMSE: {rmse_line:.2e} eV, Max: {max_err_line:.2e} eV)')
                ax2.legend()
                ax2.grid(True, alpha=0.3)

                fig6.suptitle('Python vs OpenCL Ewald (Full Evaluation)')
                fig6.tight_layout()
                fig6_path = os.path.join(args.outdir, f'{args.prefix}_fig6_opencl_1d.png')
                fig6.savefig(fig6_path, dpi=150)
                print(f"  Saved {fig6_path}")

                # Store results in report
                report['opencl_ewald'] = {
                    'vacuum_rmse': float(rmse),
                    'vacuum_max_err': float(max_err),
                    'full_rmse': float(rmse_line),
                    'full_max_err': float(max_err_line),
                    'pass': rmse < 1e-5 and rmse_line < 1e-5,
                    'figures': {
                        'fig5_opencl_vacuum': fig5_path,
                        'fig6_opencl_1d': fig6_path
                    }
                }

            except Exception as e:
                print(f"ERROR in OpenCL Ewald test: {e}")
                import traceback
                traceback.print_exc()

    print(f"\nOutput files:")
    for fig_name, fig_path in report['figures'].items():
        print(f"  {fig_name}: {fig_path}")
    if 'opencl_ewald' in report and 'figures' in report['opencl_ewald']:
        for fig_name, fig_path in report['opencl_ewald']['figures'].items():
            print(f"  {fig_name}: {fig_path}")
    print(f"  Report: {os.path.join(args.outdir, args.prefix + '_report.json')}")

    print(f"\n{'='*60}")
    print(f"Test Complete")
    print(f"{'='*60}")

    return report


if __name__ == "__main__":
    main()
