#!/usr/bin/env python3
"""combined_plot.py

Creates a publication-quality 2x2 composite figure showing:
- Top row: Energy comparison plots (perfect vs defect systems)
- Bottom row: Molecular trajectory visualizations
"""

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from pathlib import Path
import matplotlib as mpl

# Import modified plotting functions (see notes below)
from energy_comparison_1d import plot_comparison
from visualize_top_layer_xy import plot_top_layer_projections

# Set global style parameters matching energy_comparison_1d.py
l_w = 2.5
f_s = 25
plt.style.use('default')
mpl.rcParams.update({
    'font.family': 'Times New Roman',
    'axes.linewidth': l_w,
    'xtick.major.width': l_w,
    'ytick.major.width': l_w,
    'xtick.major.size': 6,
    'ytick.major.size': 6,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'xtick.minor.width': 1.5,
    'ytick.minor.width': 1.5,
    'xtick.minor.size': 4,
    'ytick.minor.size': 4,
    'xtick.labelsize': f_s,
    'ytick.labelsize': f_s,
    'axes.labelsize': f_s,
    'legend.fontsize': f_s,
    'axes.grid.which': 'both'
})

def create_composite_figure():
    """Create and save the 2x2 composite figure with consistent styling."""
    # Create figure with exact dimensions for 8x6.5 inch subplots
    # (16.5x13.5 total figure size accounts for spacing between subplots)
    fig = plt.figure(figsize=(16.5, 13.5))
    
    # Use GridSpec with exact dimensions and spacing
    gs = GridSpec(2, 2, figure=fig, 
                 width_ratios=[8, 8], 
                 height_ratios=[6.5, 6.5],
                 wspace=0.4, hspace=0.4)
    
    # --- Panel A: Perfect System Energy ---
    ax1 = fig.add_subplot(gs[0, 0])
    plot_comparison(
        lammps_file=Path("/home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan/angle45/fixedatom27/nx20/total.dat"),
        firecore_file=Path("/home/indranil/git/FireCore/tests/tMMFF/relax_perfect_line/dir_1.0_1.0_0.0/cons_26/PTCDA_20x20_26_total.dat"),
        x_range=(5, 20),
        y_range=(0.2, 0.4),
        error_range=(-0.001, 0.005),
        ax=ax1
    )
    ax1.text(0.02, 0.98, "(a)", transform=ax1.transAxes, va='top', ha='left', fontsize=f_s)
    
    # --- Panel B: Defect System Energy ---
    ax2 = fig.add_subplot(gs[0, 1])
    plot_comparison(
        lammps_file=Path("/home/indranil/Documents/Project_1/Lammps/5-relaxed_linescan_defect/angle45/defect_aligned/fixedatom27/nx20/total.dat"),
        firecore_file=Path("/home/indranil/git/FireCore/tests/tMMFF/relax_defect_aligned_line/dir_1.0_1.0_0.0/cons_26/PTCDA_20x20_26_total.dat"),
        x_range=(5, 113),
        error_range=(-0.001, 0.005),
        ax=ax2
    )
    ax2.text(0.02, 0.98, "(b)", transform=ax2.transAxes, va='top', ha='left', fontsize=f_s)
    
    # --- Panel C: Perfect System Trajectory ---
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.set_xlim(-18, 90)
    ax3.set_ylim(-18, 90)
    plot_top_layer_projections(
        trajectory_file=Path("/home/indranil/git/FireCore/tests/tMMFF/relax_perfect_line/dir_1.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz"),
        fixed_atom_idx=26,
        opposite_atom_idx=29,
        ax=ax3,
    )
    ax3.text(0.02, 0.98, "(c)", transform=ax3.transAxes, va='top', ha='left', fontsize=f_s)
    
    # --- Panel D: Defect System Trajectory ---
    ax4 = fig.add_subplot(gs[1, 1])
    plot_top_layer_projections(
        trajectory_file=Path("/home/indranil/git/FireCore/tests/tMMFF/relax_defect_aligned_line/dir_1.0_1.0_0.0/cons_26/PTCDA_20x20_26_total_trajectory.xyz"),
        fixed_atom_idx=26,
        opposite_atom_idx=29,
        ax=ax4
    )
    ax4.text(0.02, 0.98, "(d)", transform=ax4.transAxes, va='top', ha='left', fontsize=f_s)
    
    # Final adjustments and save
    plt.tight_layout()
    fig.savefig("combined_results.png", dpi=600, bbox_inches='tight')
    print("Saved composite figure to combined_results.png")
    plt.close(fig)

if __name__ == "__main__":
    create_composite_figure()