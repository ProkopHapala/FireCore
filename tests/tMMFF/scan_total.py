#!/usr/bin/env python3

# Script to generate total potential scan
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def scanPlot_uff(nscan=1000, span=(0.0,8.0), dir=(0.0,0.0,1.0), p0=(0.0,0.0,0.0), label="E", saveFig=None, saveData=None):
    ts = np.linspace(span[0], span[1], nscan, endpoint=False)
  
    poss = np.zeros((nscan, 3))
    poss[:, 0] = p0[0] + ts * dir[0]
    poss[:, 1] = p0[1] + ts * dir[1]
    poss[:, 2] = p0[2] + ts * dir[2]

    Es, Fs, Ps = mmff.scan_rigid_uff(poss, bF=True, bP=True)
    
    Es_end = Es[-1]

    if saveData is not None:
        np.savetxt(saveData, np.column_stack((ts, Es-Es_end)), header="ts\tEnergy", comments="# ")

    if saveFig is not None:
        plt.figure(figsize=(10, 6))
        plt.title(label)
        plt.plot(ts, Es-Es_end, '-', lw=1.5, label=label)
        plt.xlabel(f"Distance along ({dir[0]}_{dir[1]}_{dir[2]}) direction (\u00c5)")
        plt.ylabel(f"Energy (eV)")
        plt.grid(True, linestyle='--', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig(saveFig, dpi=300)
        plt.close()

# Default parameters
molecule = "data/xyz/old_mol_old_sub_PTCDA"
substrate = "data/xyz/Na_0.9_Cl_-0.9"
output_dir = "PTCDA_data"

# Parse command line arguments if provided
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

# Default parameters if not specified
nscan = 250  # Default value
span = (2.6, 15.0)  # Default value
dir_vec = (0.0, 0.0, 1.0)  # Default value (along z-axis)
p0 = (0.0, 0.0, 0.0)  # Default value (origin)

# Simple parameter extraction from command line
if len(sys.argv) > 4:
    nscan = int(sys.argv[4])

if len(sys.argv) > 5 and ',' in sys.argv[5]:
    span_parts = sys.argv[5].split(',')
    span = (float(span_parts[0]), float(span_parts[1]))

if len(sys.argv) > 6 and ',' in sys.argv[6]:
    dir_parts = sys.argv[6].split(',')
    if len(dir_parts) == 3:
        dir_vec = (float(dir_parts[0]), float(dir_parts[1]), float(dir_parts[2]))

if len(sys.argv) > 7 and ',' in sys.argv[7]:
    p0_parts = sys.argv[7].split(',')
    if len(p0_parts) == 3:
        p0 = (float(p0_parts[0]), float(p0_parts[1]), float(p0_parts[2]))

# Set scan parameters
scan_params = {"nscan": nscan, "span": span, "dir": dir_vec, "p0": p0, "label": f"{mol_name} total potential scan"}

# Initialize MMFF with molecule and substrate
print(f"Initializing MMFF for TOTAL potential scan with {molecule} molecule and {substrate} substrate...")
mmff.init(xyz_name=molecule, surf_name=substrate, bUFF=True, bSimple=True)
mmff.getBuffs()
mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
print(f"PLQs shape: {mmff.PLQs.shape}, non-zero values: P={np.count_nonzero(mmff.PLQs[:,0])}, L={np.count_nonzero(mmff.PLQs[:,1])}, Q={np.count_nonzero(mmff.PLQs[:,2])}")

# Generate data for total potential
print("Generating scan for total potential...")
scanPlot_uff(
    **scan_params,
    saveFig=f"{output_dir}/{mol_name}_total.png",
    saveData=f"{output_dir}/{mol_name}_total.dat"
)

print(f"Total scan completed. Data saved in {output_dir}/{mol_name}_total.dat")

for key, value in scan_params.items():
    print(f"  {key}: {value}")

# Initialize MMFF with molecule and substrate
print(f"Initializing MMFF for TOTAL potential scan with {molecule} molecule and {substrate} substrate...")
mmff.init(xyz_name=molecule, surf_name=substrate, bUFF=True, bSimple=True)
mmff.getBuffs()
mmff.setSwitches2(NonBonded=1, MMFF=1, SurfAtoms=1, GridFF=1)
print(f"PLQs shape: {mmff.PLQs.shape}, non-zero values: P={np.count_nonzero(mmff.PLQs[:,0])}, L={np.count_nonzero(mmff.PLQs[:,1])}, Q={np.count_nonzero(mmff.PLQs[:,2])}")

# Generate data for total potential
print("Generating scan for total potential...")
scanPlot_uff(**scan_params, saveFig=f"{output_dir}/{mol_name}_total.png", saveData=f"{output_dir}/{mol_name}_total.dat")

print(f"Total scan completed. Data saved in {output_dir}/{mol_name}_total.dat")
