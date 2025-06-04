#!/usr/bin/env python3

# Script to generate Coulomb potential scan (handles both 1D and 2D scans)
import sys
import os
import numpy as np
import shutil

sys.path.append("../../")
from pyBall import MMFF as mmff

# Import scan functions from run.py
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from run import scanPlot_uff, scanPlot2D_uff

# Default parameters
molecule = "data/xyz/old_mol_old_sub_PTCDA"
substrate = "data/xyz/Na_0.9_Cl_-0.9"
output_dir = "PTCDA_data"

# Parse command line arguments for paths
if len(sys.argv) > 1:
    molecule = sys.argv[1]
if len(sys.argv) > 2:
    substrate = sys.argv[2]
if len(sys.argv) > 3:
    output_dir = sys.argv[3]

# Create output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Extract molecule name for file naming
mol_name = os.path.basename(molecule)

# Copy the molecule and substrate files to output directory
mol_file = molecule
sub_file = substrate

# Check if the paths exist directly, if not try with .xyz extension
if not os.path.exists(mol_file) and not mol_file.endswith('.xyz'):
    mol_file = mol_file + '.xyz'
    
if not os.path.exists(sub_file) and not sub_file.endswith('.xyz'):
    sub_file = sub_file + '.xyz'

# Copy the molecule file if it exists
mol_basename = os.path.basename(mol_file)
if os.path.exists(mol_file):
    shutil.copy2(mol_file, os.path.join(output_dir, mol_basename))
    print(f"Copied molecule file to {os.path.join(output_dir, mol_basename)}")

# Copy the substrate file if it exists
sub_basename = os.path.basename(sub_file)
if os.path.exists(sub_file):
    shutil.copy2(sub_file, os.path.join(output_dir, sub_basename))
    print(f"Copied substrate file to {os.path.join(output_dir, sub_basename)}")

# Determine scan mode based on arguments
# If argument 5 is a digit, assume it's nscan2 for 2D scan
# Otherwise, assume it's for span in 1D scan
scan_mode = '2d' if len(sys.argv) > 5 and sys.argv[5].isdigit() else '1d'

# Initialize MMFF with molecule and substrate
print(f"Initializing MMFF for COULOMB potential {scan_mode} scan...")
mmff.init(xyz_name=molecule, surf_name=substrate, bUFF=True, bSimple=True)
mmff.getBuffs()

# Set only Coulomb components (zero out Pauli/London)
print("Setting Coulomb components only (zeroing out Pauli/London)...")
mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)

# Zero out Pauli (P) and London (L) components in PLQs
original_p = np.count_nonzero(mmff.PLQs[:,0])
original_l = np.count_nonzero(mmff.PLQs[:,1])
mmff.PLQs[:,0] = 0.0
mmff.PLQs[:,1] = 0.0

print(f"PLQs shape: {mmff.PLQs.shape}, non-zero values: P={np.count_nonzero(mmff.PLQs[:,0])}, L={np.count_nonzero(mmff.PLQs[:,1])}, Q={np.count_nonzero(mmff.PLQs[:,2])}")
print(f"Zeroed out {original_p} Pauli (P) and {original_l} London (L) components")

if scan_mode == '1d':
    # 1D scan parameters
    nscan = int(sys.argv[4]) if len(sys.argv) > 4 else 250
    
    # Extract span parameters
    span_parts = sys.argv[5].split(',') if len(sys.argv) > 5 and ',' in sys.argv[5] else ['2.6', '15.0']
    span = (float(span_parts[0]), float(span_parts[1]))
    
    # Extract direction vector
    dir_parts = sys.argv[6].split(',') if len(sys.argv) > 6 and ',' in sys.argv[6] else ['0.0', '0.0', '1.0']
    dir_vec = (float(dir_parts[0]), float(dir_parts[1]), float(dir_parts[2]))
    
    # Extract starting position
    p0_parts = sys.argv[7].split(',') if len(sys.argv) > 7 and ',' in sys.argv[7] else ['0.0', '0.0', '0.0']
    p0 = (float(p0_parts[0]), float(p0_parts[1]), float(p0_parts[2]))
    
    # Set output paths
    output_prefix = f"{output_dir}/{mol_name}_coul"
    
    # Run 1D scan
    print(f"Running 1D scan for Coulomb potential...")
    scanPlot_uff(
        nscan=nscan,
        span=span,
        dir=dir_vec,
        p0=p0,
        label=f"{mol_name} Coulomb Potential Scan",
        saveFig=f"{output_prefix}.png",
        saveData=f"{output_prefix}.dat"
    )
    print(f"Coulomb 1D scan completed. Data saved in {output_prefix}.dat")
    
else:  # 2D scan
    # 2D scan parameters
    nscan1 = int(sys.argv[4]) if len(sys.argv) > 4 else 50
    nscan2 = int(sys.argv[5]) if len(sys.argv) > 5 else 50
    
    # Extract span parameters
    span1_parts = sys.argv[6].split(',') if len(sys.argv) > 6 and ',' in sys.argv[6] else ['0.0', '4.0']
    span1 = (float(span1_parts[0]), float(span1_parts[1]))
    
    span2_parts = sys.argv[7].split(',') if len(sys.argv) > 7 and ',' in sys.argv[7] else ['0.0', '4.0']
    span2 = (float(span2_parts[0]), float(span2_parts[1]))
    
    # Extract direction vectors
    dir1_parts = sys.argv[8].split(',') if len(sys.argv) > 8 and ',' in sys.argv[8] else ['1.0', '0.0', '0.0']
    dir1 = (float(dir1_parts[0]), float(dir1_parts[1]), float(dir1_parts[2]))
    
    dir2_parts = sys.argv[9].split(',') if len(sys.argv) > 9 and ',' in sys.argv[9] else ['0.0', '1.0', '0.0']
    dir2 = (float(dir2_parts[0]), float(dir2_parts[1]), float(dir2_parts[2]))
    
    # Extract starting position
    p0_parts = sys.argv[10].split(',') if len(sys.argv) > 10 and ',' in sys.argv[10] else ['0.0', '0.0', '0.0']
    p0 = (float(p0_parts[0]), float(p0_parts[1]), float(p0_parts[2]))
    
    # Set output paths
    output_prefix = f"{output_dir}/{mol_name}_coul_2d"
    
    # Run 2D scan
    print(f"Running 2D scan for Coulomb potential...")
    scanPlot2D_uff(
        nscan1=nscan1,
        nscan2=nscan2,
        span1=span1,
        span2=span2,
        dir1=dir1,
        dir2=dir2,
        p0=p0,
        label=f"{mol_name} Coulomb Potential 2D Scan",
        saveFig=f"{output_prefix}.png",
        saveData=f"{output_prefix}.npz",
        saveDataTxt=f"{output_prefix}.dat"
    )
    print(f"Coulomb 2D scan completed. Data saved in {output_prefix}.npz")
