#!/usr/bin/env python3
"""
Plot force optimization across all improvements using saved Fcon data
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
from matplotlib.cm import get_cmap
import matplotlib.colors as mcolors

def load_force_data(opt_dir):
    """Load all force data from best_XXXXX_Fcon.npy files"""
    pattern = os.path.join(opt_dir, "best_*_Fcon.npy")
    fcon_files = sorted(glob.glob(pattern))
    
    force_data = []
    attempt_numbers = []
    
    for fcon_file in fcon_files:
        # Extract attempt number from filename
        basename = os.path.basename(fcon_file)
        parts = basename.split('_')
        
        if parts[1] == 'final':
            attempt_num = 999999  # Put final at the end
        else:
            attempt_num = int(parts[1].split('_')[0])
        
        # Load force data
        Fcon = np.load(fcon_file)
        force_data.append(Fcon)
        attempt_numbers.append(attempt_num)
        
    return force_data, attempt_numbers

def plot_force_optimization(opt_dir, ia_anchor=28, out_file="force_optimization.png", f_safe=1.0):
    """Plot force magnitude evolution across improvements"""
    
    # Load force data
    force_data, attempt_numbers = load_force_data(opt_dir)
    
    if not force_data:
        print(f"No force data found in {opt_dir}")
        return
    
    print(f"Found {len(force_data)} force data files")
    print(f"Attempt range: {min(attempt_numbers)} - {max(attempt_numbers)}")
    
    # Create figure
    fig, ax = plt.subplots(figsize=(8,6))
    
    # Create colormap for rainbow effect using reversed jet
    cmap = plt.colormaps.get_cmap('jet_r')  # Reversed jet colormap
    n_improvements = len(force_data)
    colors = [cmap(i / (n_improvements - 1)) for i in range(n_improvements)]
    
    # Plot force magnitude for each improvement
    max_f_forces = []
    
    for i, (Fcon, attempt_num, color) in enumerate(zip(force_data, attempt_numbers, colors)):
        # Fcon is already force on anchor atom (shape: n_configs × 3)
        F_mag = np.linalg.norm(Fcon, axis=1)  # Force magnitude
        
        # Check if this is the final optimization (highest attempt number)
        if attempt_num == max(attempt_numbers):
            # Final curve: extra bold and black
            ax.plot(F_mag, color='black', alpha=1.0, linewidth=5, label='Final optimized')
        else:
            # Other curves: full opacity colors
            ax.plot(F_mag, color=color, alpha=1.0, linewidth=1.5, label=f'Attempt {attempt_num}' if i % 3 == 0 else "")
        
        # Store max force for this improvement
        max_f_forces.append(np.max(F_mag))
    
    # Add penalty threshold line
    ax.axhline(y=f_safe, color='red', linestyle='--', linewidth=2, alpha=0.7, label='Penalty threshold')
    
    # Formatting
    ax.set_xlabel('Configuration index along trajectory', fontsize=10)
    ax.set_ylabel('Force |F| on anchor (eV/Å)', fontsize=10)
    ax.set_title('Force Optimization Evolution Across Improvements', fontsize=10, fontweight='bold')
    
    # Set fixed limits 0-2 eV/Å
    ax.set_ylim(0, 2)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add legend
    ax.legend(loc='upper right', fontsize=10)
    
    plt.tight_layout()
    
    # Save plot
    out_path = os.path.join(opt_dir, out_file)
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    
    print(f"Saved force optimization plot to: {out_path}")
    
    # Print summary statistics
    print(f"\n=== FORCE OPTIMIZATION SUMMARY ===")
    print(f"Initial max force: {max_f_forces[0]:.3f} eV/Å")
    print(f"Final max force: {max_f_forces[-1]:.3f} eV/Å")
    print(f"Force reduction: {(max_f_forces[0] - max_f_forces[-1]):.3f} eV/Å")
    print(f"Improvement: {(1 - max_f_forces[-1]/max_f_forces[0])*100:.1f}%")
    print(f"Final force vs threshold: {max_f_forces[-1]:.3f} vs {f_safe:.1f} eV/Å")

if __name__ == "__main__":
    opt_dir = "/home/prokophapala/git/FireCore/tests/tMMFF/opt_3d_target"
    plot_force_optimization(opt_dir)
