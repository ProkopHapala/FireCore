#!/usr/bin/env python3
"""
Compare scan results between different grid generation methods:
1. Paolo/Prokop (June 2025) - z0 = max(z)
2. Debug branch (Milan's z0 algorithm)
3. Indranil branch (Milan's z0 algorithm)

This script runs 1D z-scans using each grid dataset and compares the results.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

def run_scan_with_grid(grid_dir, molecule_xyz, output_prefix, nscan=200, z_range=(2.5, 15.0)):
    """
    Run a 1D z-scan using a specific grid dataset.
    
    Args:
        grid_dir: Directory containing the grid .npy files
        molecule_xyz: Path to molecule xyz file
        output_prefix: Prefix for output files
        nscan: Number of scan points
        z_range: (min_z, max_z) range for scan
    
    Returns:
        tuple: (z_values, energies) or None if failed
    """
    print(f"\n{'='*60}")
    print(f"Running scan with grid: {grid_dir}")
    print(f"{'='*60}")
    
    # Check if grid files exist
    required_files = ["Bspline_PLQd.npy"]
    for f in required_files:
        fpath = os.path.join(grid_dir, f)
        if not os.path.exists(fpath):
            print(f"ERROR: Required file not found: {fpath}")
            return None
    
    try:
        # Initialize MMFF with the grid
        mmff.init(
            xyz_name=molecule_xyz,
            surf_name=None,  # No surface atoms, using grid
            bMMFF=False,
            bUFF=True,
            gridFF_name=os.path.join(grid_dir, "Bspline_PLQd.npy")
        )
        
        # Set up scan positions along z-axis
        zs = np.linspace(z_range[0], z_range[1], nscan)
        poss = np.zeros((nscan, 3))
        poss[:, 2] = zs  # Scan along z
        
        # Run the scan
        Es, Fs, Ps = mmff.scan_rigid_uff(poss, bF=True, bP=True)
        
        # Normalize energies (subtract last value)
        Es_normalized = Es - Es[-1]
        
        # Save data
        output_file = f"{output_prefix}_zscan.dat"
        np.savetxt(output_file, np.column_stack((zs, Es_normalized)), 
                   header="z(A)\tEnergy(eV)", comments="# ")
        print(f"Saved scan data to: {output_file}")
        
        return zs, Es_normalized
        
    except Exception as e:
        print(f"ERROR during scan: {e}")
        import traceback
        traceback.print_exc()
        return None


def compare_scans(results, labels, output_file="scan_comparison.png"):
    """
    Create comparison plot of multiple scan results.
    
    Args:
        results: List of (z_values, energies) tuples
        labels: List of labels for each result
        output_file: Path to save the comparison plot
    """
    fig, axes = plt.subplots(2, 1, figsize=(12, 10))
    
    # Top plot: All scans overlaid
    ax1 = axes[0]
    colors = ['blue', 'red', 'green', 'orange', 'purple']
    
    for i, (result, label) in enumerate(zip(results, labels)):
        if result is not None:
            zs, Es = result
            ax1.plot(zs, Es, '-', linewidth=2, color=colors[i % len(colors)], label=label)
    
    ax1.set_xlabel("Z (Å)", fontsize=12)
    ax1.set_ylabel("Energy (eV)", fontsize=12)
    ax1.set_title("Z-Scan Comparison: Different Grid Generation Methods", fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Bottom plot: Differences from first result
    ax2 = axes[1]
    
    if results[0] is not None:
        ref_zs, ref_Es = results[0]
        
        for i, (result, label) in enumerate(zip(results[1:], labels[1:]), 1):
            if result is not None:
                zs, Es = result
                # Interpolate to match reference z values if needed
                if len(zs) == len(ref_zs) and np.allclose(zs, ref_zs):
                    diff = Es - ref_Es
                else:
                    diff = np.interp(ref_zs, zs, Es) - ref_Es
                
                ax2.plot(ref_zs, diff, '-', linewidth=2, color=colors[i % len(colors)], 
                        label=f"{label} - {labels[0]}")
                
                # Print statistics
                print(f"\n{label} vs {labels[0]}:")
                print(f"  Max diff: {np.abs(diff).max():.6e} eV")
                print(f"  Mean diff: {np.abs(diff).mean():.6e} eV")
                print(f"  RMS diff: {np.sqrt(np.mean(diff**2)):.6e} eV")
    
    ax2.set_xlabel("Z (Å)", fontsize=12)
    ax2.set_ylabel("Energy Difference (eV)", fontsize=12)
    ax2.set_title(f"Difference from {labels[0]}", fontsize=14)
    ax2.legend(fontsize=10)
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nComparison plot saved to: {output_file}")
    plt.close()


if __name__ == "__main__":
    # Define grid directories to compare
    data_base = "./data"
    
    grid_dirs = {
        "Paolo (Prokop z0)": os.path.join(data_base, "Paolo_Na_0.9_Cl_-0.9_Cl_hole_3"),
        "Indranil (Milan z0)": os.path.join(data_base, "Na_0.9_Cl_-0.9_Cl_hole_3_indranil"),
        "Debug (Milan z0)": os.path.join(data_base, "Na_0.9_Cl_-0.9_Cl_hole_3_debug"),
        "Prokop z0 (new gen)": os.path.join(data_base, "Na_0.9_Cl_-0.9_Cl_hole_3_prokop_z0"),
    }
    
    # Molecule to scan
    molecule_xyz = "./data/xyz/PTCDA.xyz"
    
    # Check which grid directories exist
    available_grids = {}
    for label, path in grid_dirs.items():
        if os.path.exists(path):
            bspline_file = os.path.join(path, "Bspline_PLQd.npy")
            if os.path.exists(bspline_file):
                available_grids[label] = path
                print(f"✓ Found: {label} -> {path}")
            else:
                print(f"✗ Missing Bspline_PLQd.npy in: {path}")
        else:
            print(f"✗ Directory not found: {path}")
    
    if len(available_grids) < 2:
        print("\nERROR: Need at least 2 grid datasets to compare")
        sys.exit(1)
    
    print(f"\nWill compare {len(available_grids)} grid datasets")
    
    # Run scans
    results = []
    labels = []
    
    for label, grid_dir in available_grids.items():
        output_prefix = os.path.join("scan_results", label.replace(" ", "_").replace("(", "").replace(")", ""))
        os.makedirs("scan_results", exist_ok=True)
        
        result = run_scan_with_grid(grid_dir, molecule_xyz, output_prefix)
        results.append(result)
        labels.append(label)
    
    # Compare results
    if any(r is not None for r in results):
        compare_scans(results, labels, "scan_results/scan_comparison.png")
    else:
        print("\nERROR: No successful scans to compare")
