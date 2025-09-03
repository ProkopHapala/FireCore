#!/usr/bin/env python3
"""
Test script for comparing CPU and GPU implementations of the Universal Force Field (UFF)
using a unified library interface.
"""

import sys
import os
import argparse
import numpy as np

# Add FireCore to the Python path
base_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(base_path)

from pyBall import MMFF_multi as uff

# Define the path to common resource files
data_dir = os.path.join(base_path, "cpp/common_resources")

def run_uff(use_gpu, components):
    """
    Run a single UFF evaluation step on either CPU or GPU.

    Args:
        use_gpu (bool): If True, run on GPU (OpenCL), otherwise run on CPU.
        components (dict): A dictionary of flags to enable/disable UFF components.

    Returns:
        tuple: A tuple containing the calculated energy (placeholder) and forces.
    """
    print(f"\n--- Running UFF on {'GPU' if use_gpu else 'CPU'} ---")
    
    # Set which components to evaluate for this run
    uff.setSwitchesUFF(
        DoBond=components.get('bonds', -1),
        DoAngle=components.get('angles', -1),
        DoDihedral=components.get('dihedrals', -1),
        DoInversion=components.get('inversions', -1),
        DoAssemble=1,
        SubtractBondNonBond=1, # This is usually needed for UFF
        ClampNonBonded=1
    )
    
    # Select CPU or GPU path via the iParalel flag
    # Convention: iParalel > 1 means GPU, otherwise CPU
    iParalel = 3 if use_gpu else 0
    
    # Run a single evaluation step
    uff.run(nstepMax=1, iParalel=iParalel)
    
    # Download results from GPU if necessary. The CPU path modifies host buffers directly.
    if use_gpu:
        uff.download(bForces=True)
    
    # Get access to the host-side buffers
    uff.getBuffs()
    
    # Extract results. For both CPU and GPU, the result is in the first system's buffer.
    # Note: Energy comparison is not fully implemented yet.
    energy = 0 
    forces = uff.gpu_aforces[0, :, :3].copy()
    
    return energy, forces

def compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces, tol=1e-5, component_name=""):
    """Compare energy and forces from CPU and GPU, and report differences."""
    print(f"\n--- {component_name} Comparison ---")
    print(f"CPU Energy: {cpu_energy:.6f} | GPU Energy: {gpu_energy:.6f} (NOTE: GPU energy not fully implemented for comparison)")
    
    if cpu_forces.shape != gpu_forces.shape:
        print(f"ERROR: Force array shapes mismatch! CPU: {cpu_forces.shape}, GPU: {gpu_forces.shape}")
        return False

    force_diff = np.abs(cpu_forces - gpu_forces)
    max_force_diff = np.max(force_diff) if force_diff.size > 0 else 0.0
    
    print(f"Max force component difference: {max_force_diff:.6f}")
    
    if max_force_diff > tol:
        print(f"Validation FAILED for {component_name}: Forces differ more than tolerance.")
        print("CPU Forces:\n", cpu_forces)
        print("GPU Forces:\n", gpu_forces)
        print("Difference:\n", force_diff)
        return False
    else:
        print(f"Validation PASSED for {component_name}: Forces are consistent.")
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='UFF CPU vs. GPU Validation Test')
    default_mol = os.path.join(data_dir, 'mol', 'formic_acid.mol2')
    parser.add_argument('-m', '--molecule', type=str, default=default_mol, help='Molecule file (.mol2, .xyz)')
    parser.add_argument('-t', '--tolerance', type=float, default=1e-5, help='Numerical tolerance for comparison')
    args = parser.parse_args()

    # --- Initialize the library once ---
    print("--- Initializing MMFF_multi library ---")
    uff.init(
        nSys_=1,
        xyz_name=args.molecule,
        sElementTypes=os.path.join(data_dir, "ElementTypes.dat"),
        sAtomTypes=os.path.join(data_dir, "AtomTypes.dat"),
        sBondTypes=os.path.join(data_dir, "BondTypes.dat"),
        sAngleTypes=os.path.join(data_dir, "AngleTypes.dat"),
        sDihedralTypes=os.path.join(data_dir, "DihedralTypes.dat"),
        bMMFF=True, # to use UFF MMFF should be True!
        bUFF=True
    )
    # NOTE: setSwitches() does not accept the 'NonBonded' keyword; non-bonded handling for UFF
    # is controlled via setSwitchesUFF (e.g., ClampNonBonded, SubtractBondNonBond). We already
    # set ClampNonBonded/SubtractBondNonBond inside run_uff(), so we disable this invalid call.  # TODO
    # uff.setSwitches(NonBonded=-1)

    all_components = ['bonds', 'angles', 'dihedrals', 'inversions']
    all_passed = True

    for component in all_components:
        print(f"\n================ Testing component: {component} ================")
        
        # Define which components to enable for this test run
        component_flags = {c: (1 if c == component else -1) for c in all_components}
        
        # Run CPU version
        cpu_energy, cpu_forces = run_uff(use_gpu=False, components=component_flags)
        
        # Run GPU version
        gpu_energy, gpu_forces = run_uff(use_gpu=True, components=component_flags)
        
        # Compare results
        if not compare_results(cpu_energy, cpu_forces, gpu_energy, gpu_forces, tol=args.tolerance, component_name=component):
            all_passed = False

    print("\n================= SUMMARY =================")
    if all_passed:
        print("All UFF components PASSED validation!")
    else:
        print("One or more UFF components FAILED validation.")
    print("=========================================")
