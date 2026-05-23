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


def test_system(system_name, test_conventions=None, atom_req=(1.487, 0.0006808, 0.0, 0.0), alpha_morse=1.5, verbose=True):
    """
    Test alignment for a given system.
    
    Args:
        system_name: Name of system (e.g., 'NaCl_1x1_L3', 'NaCl_8x8_L3_ClHole')
        test_conventions: List of conventions to test (None = test all)
        atom_req: Tuple (R, E, Q, H) for test atom
        alpha_morse: Alpha parameter for REQ to PLQ conversion
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
                                       atom_req=atom_req, alpha_morse=alpha_morse, verbose=verbose)
    return report


def main():
    import argparse
    parser = argparse.ArgumentParser(description='Test GridFF alignment on systems')
    parser.add_argument('system', nargs='?', default='NaCl_1x1_L3',  help='System name (e.g., NaCl_1x1_L3, NaCl_8x8_L3_ClHole)')
    parser.add_argument('--conventions', nargs='+', help='Specific conventions to test (default: all)')
    args = parser.parse_args()
    report = test_system(args.system, test_conventions=args.conventions, verbose=True)
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
