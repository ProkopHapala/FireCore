#!/usr/bin/env python3
"""
Test script for GridFF alignment verification on NaCl systems.

Uses existing GridFF files in cpp/common_resources/:
- NaCl_1x1_L3/
- NaCl_8x8_L3_ClHole/

Generates visualization plots to verify grid-atom alignment.
"""

import os
import sys
import numpy as np

# Add FireCore to path
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))
if ROOT not in sys.path: sys.path.insert(0, ROOT)

from pyBall.OCL.Surface_utils import run_alignment_verification


def test_system(system_name, test_conventions=None, atom_req=(1.487, 0.0006808, 0.0, 0.0), alpha_morse=1.5, z_atom_range=2.0, verbose=True):
    """
    Test alignment for a given system.
    
    CRITICAL: atom_req E parameter MUST be sqrt(EvdW), NOT raw EvdW.
    If reading from ElementTypes.dat, you MUST sqrt the E value.
    This matches the GridFF generation convention in ocl_GridFF_new.py::make_atoms_arrays(bSqrtEvdw=True).
    
    Args:
        system_name: Name of system (e.g., 'NaCl_1x1_L3', 'NaCl_8x8_L3_ClHole')
        test_conventions: List of conventions to test (None = test all)
        atom_req: Tuple (R, sqrt(EvdW), Q, H) for test atom - E MUST be sqrt(EvdW)
        alpha_morse: Alpha parameter for REQ to PLQ conversion (must match GridFF generation)
        z_atom_range: Z-range below top atom to show in XY plots (default 2.0 A)
        verbose: Print progress
        
    Returns:
        report dict or None if error
    """
    print("\n" + "="*70)
    print(f"TEST: {system_name}")
    print("="*70)
    
    # Auto-detect paths
    grid_path = os.path.join(ROOT, 'cpp', 'common_resources', system_name, 'Bspline_PLQd.npy')
    sub_path  = os.path.join(ROOT, 'cpp', 'common_resources', 'xyz', f'{system_name}.xyz')
    out_dir   = os.path.join(ROOT, 'tests', 'tMMFF', 'output_alignment', system_name)
    
    if not os.path.exists(grid_path):
        print(f"ERROR: GridFF not found: {grid_path}")
        return None
    if not os.path.exists(sub_path):
        print(f"ERROR: Substrate not found: {sub_path}")
        return None
    
    print(f"GridFF:   {grid_path}")
    print(f"Substrate: {sub_path}")
    print(f"Output:   {out_dir}")
    
    # Run verification
    report = run_alignment_verification(grid_path, sub_path, out_dir,  test_conventions=test_conventions, 
                                       atom_req=atom_req, alpha_morse=alpha_morse, z_atom_range=z_atom_range, verbose=verbose)
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Test GridFF alignment on systems')
    parser.add_argument('system', nargs='?', default='NaCl_1x1_L3',  help='System name (e.g., NaCl_1x1_L3, NaCl_8x8_L3_ClHole)')
    parser.add_argument('--conventions', nargs='+', help='Specific conventions to test (default: all)')
    parser.add_argument('--z-atom-range', type=float, default=2.0, help='Z-range below top atom to show in XY plots (default 2.0 A)')
    parser.add_argument('--profile', action='store_true', help='Run with cProfile profiler')
    args = parser.parse_args()
    
    if args.profile:
        import cProfile
        import pstats
        from io import StringIO
        pr = cProfile.Profile()
        pr.enable()
        report = test_system(args.system, test_conventions=args.conventions, z_atom_range=args.z_atom_range, verbose=True)
        pr.disable()
        
        s = StringIO()
        ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
        ps.print_stats(30)
        print("\n" + "="*70)
        print("PROFILE (top 30 by cumulative time)")
        print("="*70)
        print(s.getvalue())
    else:
        report = test_system(args.system, test_conventions=args.conventions, z_atom_range=args.z_atom_range, verbose=True)
    
    if report:
        print("\n" + "="*70)
        print("SUMMARY")
        print("="*70)
        print(f"Best convention: {report['best_convention']}")
        print(f"g0: {report['best_g0']}")
        print(f"Aligned: {report['alignment_stats']['aligned']}")
        print(f"Max error: {report['alignment_stats']['max_error']:.4f} Å")
        print(f"Output: {report['summary_plot']}")
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
