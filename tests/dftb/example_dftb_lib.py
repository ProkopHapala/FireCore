#!/usr/bin/env python
"""
Lightweight example using dftb_lib.py for Hamiltonian/DM extraction
"""

import sys
sys.path.insert(0, '/home/prokophapala/git/FireCore')

from pyBall.dftb_lib import DftbPlusCalculator, calculate_with_matrices

print("=== DFTB+ Hamiltonian/DM Extraction Example ===\n")

# Method 1: Using the class (full control)
print("Method 1: Using DftbPlusCalculator class")
print("-" * 50)

calc = DftbPlusCalculator()
calc.initialize(input_file="/home/prokophapala/git_SW/dftbplus/test/input.dftb")
print(f"System: {calc.nr_atoms} atoms, basis size {calc.basis_size}")

calc.register_callbacks()
energy = calc.calculate()
print(f"Energy: {energy:.6f} Ha")

# Get matrices
H = calc.get_hamiltonian()
S = calc.get_overlap()
DM = calc.get_density_matrix()

print(f"\nHamiltonian shape: {H.shape}")
print(f"Overlap shape: {S.shape}")
print(f"Density matrix shape: {DM.shape}")

# Calculate properties
electron_count = calc.get_electron_count()
print(f"Electron count (Tr(S*DM)): {electron_count:.6f}")

eigenvalues = calc.get_eigenvalues()
print(f"Eigenvalues: {eigenvalues}")

calc.finalize()
print("\n" + "="*50 + "\n")

# Method 2: Using convenience function (one-liner)
print("Method 2: Using calculate_with_matrices() convenience function")
print("-" * 50)

result = calculate_with_matrices(
    input_file="/home/prokophapala/git_SW/dftbplus/test/input.dftb"
)

print(f"Energy: {result['energy']:.6f} Ha")
print(f"Electron count: {result['electron_count']:.6f}")
print(f"Hamiltonian shape: {result['hamiltonian'].shape}")

print("\n=== Example Complete ===")
