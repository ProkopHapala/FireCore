#!/usr/bin/env python
"""
Example: Compute Hessian matrix from DFTB+ using native SecondDerivatives driver

This demonstrates Hessian extraction using DFTB+'s built-in SecondDerivatives driver.
The hessian.out file is read programmatically with numpy.
"""

import sys
sys.path.insert(0, '/home/prokophapala/git/FireCore')

import numpy as np
import argparse
import subprocess
import os
from pathlib import Path
from ase import Atoms
from ase.io import read

# Import helper functions from dftb_utils
from pyBall import dftb_utils as dftbu

def main():
    parser = argparse.ArgumentParser(
        description='Compute Hessian matrix from DFTB+ using SecondDerivatives driver'
    )
    parser.add_argument('molecule', nargs='?', help='XYZ file to load (default: built-in H2O)',default=None )
    parser.add_argument('--use-ase', action='store_true', help='Use ASE to load XYZ file (default: AtomicSystem)' )
    parser.add_argument('--sk-path', type=str,default='/home/prokophapala/git_SW/asi/tests/testcases/test_expdmhs.dftbp/', help='Path to Slater-Koster files (default: /home/prokophapala/git_SW/asi/tests/testcases/test_expdmhs.dftbp/)')
    parser.add_argument('--delta',type=float, default=1e-4, help='Finite difference step in atomic units (default: 1e-4)'  )
    parser.add_argument('--no-freq', action='store_true',help='Skip frequency calculation (only output Hessian)' )
    parser.add_argument('--workdir', type=str, default='hessian_calc', help='Working directory for DFTB+ calculation (default: hessian_calc)' )
    parser.add_argument('--modes-output', type=str, default='vibration_modes.xyz', help='Output file for vibration modes in Jmol XYZ format (default: vibration_modes.xyz)'  )
    parser.add_argument('--vector-scale',type=float,default=1.0,help='Scaling factor for vibration vectors (default: 1.0)')
    args = parser.parse_args()
    
    print("=== DFTB+ Hessian Extraction Example (SecondDerivatives Driver) ===\n")
    
    # Load molecule
    if args.molecule:
        print(f"Loading molecule from: {args.molecule}")
        atoms = dftbu.load_molecule(args.molecule, use_ase=args.use_ase)
    else:
        # Default: H2O molecule (proper geometry)
        print("Using built-in H2O molecule")
        atoms = Atoms('H2O', positions=[
            [0.0, 0.0, 0.0],  # O at origin
            [0.957, 0.0, 0.0],  # H1 (bond length 0.957 Å)
            [-0.240, 0.927, 0.0]  # H2 (bond length 0.957 Å, angle 104.5°)
        ])
    
    print(f"System: {len(atoms)} atoms")
    print(f"Formula: {atoms.get_chemical_formula()}")
    print()
    
    # Create working directory
    workdir = Path(args.workdir).absolute()
    workdir.mkdir(exist_ok=True, parents=True)
    os.chdir(workdir)
    
    # Write XYZ file
    from pyBall import atomicUtils as au
    au.saveXYZ(es=atoms.get_chemical_symbols(), xyzs=atoms.positions, fname="geo.xyz")
    
    # Write DFTB+ input with SecondDerivatives driver manually
    print("Writing DFTB+ input file...")
    
    dftbu.write_dftb_input_hessian(atoms.get_chemical_symbols(), gname="geo.xyz", fname='dftb_in.hsd', basis_path=args.sk_path, delta=args.delta)
    
    print("Running DFTB+ with SecondDerivatives driver...")
    print(f"  Delta: {args.delta} atomic units")
    print(f"  Working directory: {workdir.absolute()}")
    print()
    
    # Run DFTB+
    result = subprocess.run(['dftb+'], capture_output=True, text=True)
    
    if result.returncode != 0:
        print("DFTB+ failed!")
        print(result.stderr)
        return
    
    # Read Hessian
    if not os.path.exists('hessian.out'):
        print("Error: hessian.out not found!")
        print("DFTB+ output:")
        print(result.stdout)
        return
    
    print("Reading Hessian from hessian.out...")
    H_hartree_bohr = dftbu.read_hessian('hessian.out', n_atoms=len(atoms))
    
    # Convert to eV/Å²
    H = dftbu.hessian_hartree_bohr_to_eV_angstrom(H_hartree_bohr)
    
    print(f"Hessian shape: {H.shape}")
    print(f"Hessian units: eV/Å²")
    print()
    
    # Print Hessian matrix
    print("Hessian matrix (3N × 3N):")
    np.set_printoptions(precision=6, suppress=True, linewidth=120)
    print(H)
    print()
    
    # Convert to mass-weighted dynamical matrix
    masses = atoms.get_masses()
    print(f"Atomic masses: {masses}")
    
    D, im = dftbu.hessian_to_mass_weighted(H, masses)
    print(f"Mass-weighted Hessian shape: {D.shape}")
    print()
    
    # Convert to vibrational frequencies
    if not args.no_freq:
        frequencies, modes = dftbu.hessian_to_frequencies(H, masses)
        print(f"All vibrational frequencies (cm⁻¹):")
        for i, freq in enumerate(frequencies):
            print(f"  Mode {i:2d}: {freq:10.2f} cm⁻¹")
        print()
        
        # Save vibration modes in Jmol format (only modes above threshold)
        modes_file = workdir / args.modes_output
        print(f"Saving vibration modes to {modes_file}")
        dftbu.write_vibration_modes_jmol( str(modes_file),  atoms,  frequencies,  modes,  scale=args.vector_scale, min_freq=10.0 )
        print()
    
    # Return to original directory
    os.chdir('..')
    
    print("\n=== Example Complete ===")
    print(f"Results saved in: {workdir.absolute()}")

if __name__ == '__main__':
    main()
