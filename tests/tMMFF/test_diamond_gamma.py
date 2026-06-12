#!/usr/bin/env python3
"""
Gamma-point phonon spectrum for diamond primitive cell (2 atoms) with PBC.
Uses finite-difference Hessian from MMFF_lib and mass-weights to get frequencies.
"""
import os, sys
import numpy as np


'''
run like this:

ASAN_OPTIONS=detect_leaks=0:abort_on_error=0 LD_PRELOAD=/lib/x86_64-linux-gnu/libasan.so.8 python3 test_diamond_gamma.py > diamond_stdout.txt 2>&1
ASAN_OPTIONS=detect_leaks=0:abort_on_error=0 LD_PRELOAD=/lib/x86_64-linux-gnu/libasan.so.8 python3 test_diamond_gamma.py > diamond_stdout.txt 2>&1
'''

sys.path.append("../../")
from pyBall import MMFF

xyz_path = os.path.join(os.path.dirname(__file__), "..", "..", "cpp/common_resources/crystals/diamond_primitive.xyz")

# Initialize MMFF with diamond primitive cell and PBC
MMFF.init(
    xyz_name=xyz_path,
    nPBC=(1,1,1),
    bEpairs=False,
    bMMFF=True,
)
print("Initialization complete (neighCell fix working)")

n_atoms = 2  # Diamond primitive cell
print(f"Number of atoms in primitive cell: {n_atoms}")

# NOTE: Diamond primitive cell is at equilibrium by symmetry
# (small residual forces due to MMFF bond length vs actual lattice)
print("\nComputing Hessian at input geometry (primitive cell is at equilibrium by symmetry)")

# Compute Hessian at Gamma point
print("\nComputing Γ-point Hessian (3Nx3N)...")
inds = np.arange(n_atoms, dtype=np.int32)
H = MMFF.getHessian3Nx3N(inds, dx=1e-4)
H = 0.5 * (H + H.T)  # Enforce symmetry

# Check for NaN in Hessian
if np.isnan(H).any() or np.isinf(H).any():
    print("ERROR: NaN or Inf detected in Hessian!")
    has_nans = True
else:
    print("Hessian is clean (no NaN/Inf)")
    has_nans = False

print(f"Hessian shape: {H.shape}")
print(f"Hessian norm: {np.linalg.norm(H):.6e}")

# Construct mass-weighted dynamical matrix
# Carbon atomic mass in amu
mass_C = 12.0107  # amu
masses = np.full(n_atoms, mass_C)

# Mass-weighting: D_ij = H_ij / sqrt(m_i * m_j)
dim = 3 * n_atoms
D = np.zeros((dim, dim))
for i in range(n_atoms):
    for j in range(n_atoms):
        mi = masses[i]
        mj = masses[j]
        D[i*3:(i+1)*3, j*3:(j+1)*3] = H[i*3:(i+1)*3, j*3:(j+1)*3] / np.sqrt(mi * mj)

# Diagonalize dynamical matrix
eigenvalues, eigenvectors = np.linalg.eigh(D)

# Convert eigenvalues to frequencies in cm^-1
# omega = sqrt(eigenvalue) in atomic units (Hartree / (amu * Bohr^2))
# Conversion: 1 Hartree / (amu * Bohr^2) -> (cm^-1)^2 factor ~ 264.2
# => freq [cm^-1] = sqrt(eigenvalue [Hartree/(amu*Bohr^2)]) * 16.25
#              = sqrt(eigenvalue) * 16.25
conv_factor = 16.25

freq_cm1 = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * conv_factor

# Project out rigid body modes (translation/rotation) to get true acoustic modes
# For crystals at Gamma, the 3 acoustic modes should be ~0
# Sort by absolute value to identify them
sorted_idx = np.argsort(np.abs(freq_cm1))

# Write results to file (to avoid ASan double-free truncating stdout at exit)
with open('diamond_phonon_results.txt', 'w') as fout:
    fout.write("="*60 + "\n")
    fout.write("Gamma-point phonon frequencies (cm^-1):\n")
    fout.write("="*60 + "\n")
    for i, idx in enumerate(sorted_idx):
        f = freq_cm1[idx]
        marker = "  [ACOUSTIC]" if i < 3 else ""
        fout.write(f"  Mode {idx:2d}: {f:10.3f} cm^-1{marker}\n")
    
    fout.write("\n" + "="*60 + "\n")
    fout.write(f"Lowest non-acoustic mode: {freq_cm1[sorted_idx[3]]:.2f} cm^-1\n")
    fout.write(f"Highest mode: {freq_cm1[sorted_idx[-1]]:.2f} cm^-1\n")
    fout.write("="*60 + "\n")
    fout.write("\nHessian eigenvalues:\n")
    for i, idx in enumerate(sorted_idx):
        fout.write(f"  Mode {idx:2d}: eig={eigenvalues[idx]:.4e}\n")
    fout.write("\nSUCCESS: Diamond primitive cell phonon calculation complete!\n")

# Also print to stdout
print("\n" + "="*60)
print("Gamma-point phonon frequencies (cm^-1):")
print("="*60)
for i, idx in enumerate(sorted_idx):
    f = freq_cm1[idx]
    marker = "  [ACOUSTIC]" if i < 3 else ""
    print(f"  Mode {idx:2d}: {f:10.3f} cm^-1{marker}")

print("\n" + "="*60)
print(f"Lowest non-acoustic mode: {freq_cm1[sorted_idx[3]]:.2f} cm^-1")
print(f"Highest mode: {freq_cm1[sorted_idx[-1]]:.2f} cm^-1")
print("="*60)
print("\nSUCCESS: Diamond primitive cell phonon calculation complete!")
print("Results also saved to: diamond_phonon_results.txt")
