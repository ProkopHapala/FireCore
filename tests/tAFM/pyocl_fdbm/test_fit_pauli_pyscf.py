#!/usr/bin/env python3
"""
Test script for pySCF Pauli parameter fitting.

Usage:
    # Fit with pre-computed FDBM grids and z-scan reference
    python3 test_fit_pauli_pyscf.py pentacene.xyz --basis sto-3g --target_indices 0,1,20,21 \
        --fdbm_dir ./fdbm_grids_pyscf_sto3g --zscan_dir ./zscan_pyscf_sto3g

    # Generate FDBM grids and z-scan reference on-the-fly, then fit
    python3 test_fit_pauli_pyscf.py pentacene.xyz --basis sto-3g --target_indices 0,1,20,21 \
        --generate_ref

    # Use DFTB method instead of RHF
    python3 test_fit_pauli_pyscf.py pentacene.xyz --basis sto-3g --method RKS --xc lda,vwn \
        --target_indices 0,1,20,21 --generate_ref
"""

import sys, os, argparse
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.realpath(os.path.join(_THIS_DIR, '..', '..', '..'))
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)

from pyBall.OCL import AFM_utils as afm_utils

def main():
    parser = argparse.ArgumentParser(description='Fit pySCF Pauli parameters')
    parser.add_argument('xyz_file',  type=str, help='Input XYZ file')
    parser.add_argument('--basis',   type=str, default='sto-3g',  help='pySCF basis set (sto-3g, 6-31g, etc.)')
    parser.add_argument('--method',  type=str, default='RHF', choices=['RHF', 'RKS'],  help='pySCF SCF method')
    parser.add_argument('--xc',      type=str, default=None,help='DFT XC functional for RKS (e.g., lda,vwn, pbe)')
    parser.add_argument('--target_indices', type=str, default='0', help='Comma-separated atom indices to fit (e.g., "0,1,20,21")')
    parser.add_argument('--fdbm_dir',   type=str, default=None, help='Pre-computed FDBM grid directory')
    parser.add_argument('--zscan_dir',  type=str, default=None, help='Pre-computed pySCF z-scan directory')
    parser.add_argument('--output_dir', type=str, default='fit_pauli_pyscf',help='Output directory for fitting results')
    parser.add_argument('--generate_ref', action='store_true', help='Generate pySCF z-scan reference if missing')
    parser.add_argument('--z_min',   type=float, default=2.0,  help='Minimum z for fit range [Å]')
    parser.add_argument('--z_max',   type=float, default=3.0,  help='Maximum z for fit range [Å]')
    parser.add_argument('--step',    type=float, default=0.15, help='Grid step for FDBM generation [Å]')
    parser.add_argument('--margin',  type=float, default=4.0,  help='Grid margin for FDBM generation [Å]')
    parser.add_argument('--z_extra', type=float, default=5.0,  help='Extra z-space for FDBM generation [Å]')
    parser.add_argument('--tip_xyz', type=str,   default='CO.xyz',  help='CO tip XYZ file')
    args = parser.parse_args()

    target_indices = [int(x.strip()) for x in args.target_indices.split(',')]
    
    print("="*70)
    print("pySCF Pauli Parameter Fitting")
    print("="*70)
    print(f"XYZ file: {args.xyz_file}")
    print(f"Basis: {args.basis}")
    print(f"Method: {args.method}")
    print(f"XC: {args.xc if args.xc else 'N/A'}")
    print(f"Target atoms: {target_indices}")
    print(f"Fit range: z=[{args.z_min}, {args.z_max}] Å")
    print(f"Generate reference: {args.generate_ref}")
    print("="*70)

    results = afm_utils.fit_pauli_parameters_pyscf(
        xyz_file=args.xyz_file,
        pyscf_basis=args.basis,
        pyscf_method=args.method,
        pyscf_xc=args.xc,
        target_indices=target_indices,
        fdbm_dir=args.fdbm_dir,
        zscan_dir=args.zscan_dir,
        output_dir=args.output_dir,
        z_min=args.z_min,
        z_max=args.z_max,
        generate_ref=args.generate_ref,
        step=args.step,
        margin=args.margin,
        z_extra=args.z_extra,
        tip_xyz=args.tip_xyz
    )
    
    print("\n" + "="*70)
    print("Fitting Results Summary")
    print("="*70)
    print(f"Basis: {results['basis']}")
    print(f"Method: {results['method']}")
    print(f"XC: {results['xc'] if results['xc'] else 'N/A'}")
    if results['A_mean'] is not None:
        print(f"A_pauli: {results['A_mean']:.2f} ± {results['A_std']:.2f}")
        print(f"beta:    {results['beta_mean']:.4f} ± {results['beta_std']:.4f}")
    print("="*70)

if __name__ == "__main__":
    main()
