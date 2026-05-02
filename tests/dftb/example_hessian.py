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
from pyBall.AtomicSystem import AtomicSystem

def load_molecule_xyz(filename, use_ase=True):
    """Load molecule from XYZ file using ASE or AtomicSystem."""
    if use_ase:
        atoms = read(filename)
        return atoms
    else:
        # Load using AtomicSystem
        mol = AtomicSystem(fname=filename)
        # Convert to ASE Atoms
        atoms = Atoms(
            symbols=mol.enames,
            positions=mol.apos,
            cell=mol.lvec if mol.lvec is not None else None,
            pbc=(mol.lvec is not None)
        )
        return atoms

def read_hessian(filename='hessian.out', n_atoms=None):
    """Read DFTB+ hessian.out file programmatically.
    
    DFTB+ writes the Hessian in a formatted text file with varying columns.
    We need to read all numbers and reshape them.
    """
    # Read all numbers from the file
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # Extract all floating point numbers
    numbers = []
    for line in lines:
        parts = line.split()
        for part in parts:
            try:
                numbers.append(float(part))
            except ValueError:
                pass
    
    data = np.array(numbers)
    
    if n_atoms is None:
        n_total = int(np.sqrt(len(data)))
        if n_total * n_total != len(data):
            raise ValueError(f"Cannot infer n_atoms from {len(data)} elements")
        n_atoms = n_total // 3
    
    expected = (3 * n_atoms) ** 2
    if len(data) != expected:
        raise ValueError(f"Expected {expected} elements for {n_atoms} atoms, got {len(data)}")
    
    hessian = data.reshape((3 * n_atoms, 3 * n_atoms))
    return hessian  # Hartree/Bohr²

def hessian_hartree_bohr_to_eV_angstrom(hessian):
    """Convert Hessian from Hartree/Bohr² to eV/Å²."""
    # 1 Hartree = 27.2114 eV
    # 1 Bohr = 0.529177 Å
    # So Hartree/Bohr² = 27.2114 / 0.529177² = 97.207 eV/Å²
    conversion = 27.2114 / (0.529177 ** 2)
    return hessian * conversion

def main():
    parser = argparse.ArgumentParser(
        description='Compute Hessian matrix from DFTB+ using SecondDerivatives driver'
    )
    parser.add_argument(
        'molecule',
        nargs='?',
        help='XYZ file to load (default: built-in H2O)',
        default=None
    )
    parser.add_argument(
        '--use-ase',
        action='store_true',
        help='Use ASE to load XYZ file (default: AtomicSystem)'
    )
    parser.add_argument(
        '--sk-path',
        type=str,
        default='/home/prokophapala/git_SW/asi/tests/testcases/test_expdmhs.dftbp/',
        help='Path to Slater-Koster files (default: /home/prokophapala/git_SW/asi/tests/testcases/test_expdmhs.dftbp/)'
    )
    parser.add_argument(
        '--delta',
        type=float,
        default=1e-4,
        help='Finite difference step in atomic units (default: 1e-4)'
    )
    parser.add_argument(
        '--no-freq',
        action='store_true',
        help='Skip frequency calculation (only output Hessian)'
    )
    parser.add_argument(
        '--workdir',
        type=str,
        default='hessian_calc',
        help='Working directory for DFTB+ calculation (default: hessian_calc)'
    )
    parser.add_argument(
        '--modes-output',
        type=str,
        default='vibration_modes.xyz',
        help='Output file for vibration modes in Jmol XYZ format (default: vibration_modes.xyz)'
    )
    parser.add_argument(
        '--vector-scale',
        type=float,
        default=1.0,
        help='Scaling factor for vibration vectors (default: 1.0)'
    )
    
    args = parser.parse_args()
    
    print("=== DFTB+ Hessian Extraction Example (SecondDerivatives Driver) ===\n")
    
    # Load molecule
    if args.molecule:
        print(f"Loading molecule from: {args.molecule}")
        atoms = load_molecule_xyz(args.molecule, use_ase=args.use_ase)
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
    
    # Determine angular momenta from elements
    elem_set = set(atoms.get_chemical_symbols())
    max_angular = {}
    for elem in elem_set:
        if elem in ['H']:
            max_angular[elem] = '"s"'
        elif elem in ['C', 'N', 'O']:
            max_angular[elem] = '"p"'
        elif elem in ['S', 'Si']:
            max_angular[elem] = '"d"'
        else:
            max_angular[elem] = '"p"'  # default
    
    # Create working directory
    workdir = Path(args.workdir).absolute()
    workdir.mkdir(exist_ok=True, parents=True)
    os.chdir(workdir)
    
    # Write XYZ file
    from pyBall import atomicUtils as au
    au.saveXYZ(es=atoms.get_chemical_symbols(), xyzs=atoms.positions, fname="geo.xyz")
    
    # Write DFTB+ input with SecondDerivatives driver manually
    print("Writing DFTB+ input file...")
    
    # Build angular momentum string
    max_angular_str = '\n'.join([f'        {elem} = {max_angular[elem]}' for elem in max_angular])
    
    hsd_content = f'''Geometry = xyzFormat {{
    <<< "geo.xyz"
}}

Driver = SecondDerivatives {{
    Delta = {args.delta}
    Atoms = 1:-1
}}

Hamiltonian = DFTB {{
    Scc = Yes
    SlaterKosterFiles = Type2FileNames {{
        Prefix = {args.sk_path}
        Separator = "-"
        Suffix = ".skf"
    }}
    MaxAngularMomentum {{
{max_angular_str}
    }}
    SCCTolerance = 1.000000e-07
}}
'''
    
    with open('dftb_in.hsd', 'w') as f:
        f.write(hsd_content)
    
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
    H_hartree_bohr = read_hessian('hessian.out', n_atoms=len(atoms))
    
    # Convert to eV/Å²
    H = hessian_hartree_bohr_to_eV_angstrom(H_hartree_bohr)
    
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
        dftbu.write_vibration_modes_jmol(
            str(modes_file), 
            atoms, 
            frequencies, 
            modes, 
            scale=args.vector_scale,
            min_freq=10.0  # Higher threshold to skip translations/rotations
        )
        print()
    
    # Return to original directory
    os.chdir('..')
    
    print("\n=== Example Complete ===")
    print(f"Results saved in: {workdir.absolute()}")

if __name__ == '__main__':
    main()
