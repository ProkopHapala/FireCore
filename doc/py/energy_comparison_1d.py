import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LightSource

# Re-use the generic visualiser that already exists in FireCore
from plot_utils import MoleculeTrajectoryVisualizer
from scipy.spatial import ConvexHull
import matplotlib as mpl
from matplotlib import font_manager
import os

# 2. Register all Times New Roman variants explicitly
font_files = [
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Italic.ttf',
    '/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman_Bold_Italic.ttf'
]

for font_file in font_files:
    if os.path.exists(font_file):
        font_manager.fontManager.addfont(font_file)
    else:
        print(f"Warning: Font file not found - {font_file}")

mpl.rcParams['font.family'] = 'Times New Roman'

l_w=2.5
f_s=25
# Adjusting Matplotlib's default settings
mpl.rcParams.update({
    'axes.linewidth': l_w,          # Thickness of the axix lines
    'xtick.major.width': l_w,       # Thickness of major ticks on x-axis
    'ytick.major.width': l_w,       # Thickness of major ticks on y-axis
    'xtick.minor.width': 1.5,       # Thickness of minor ticks on x-axis
    'ytick.minor.width': 1.5,       # Thickness of minor ticks on y-axis
    'xtick.major.size': 6,          # Length of major ticks on x-axis
    'ytick.major.size': 6,          # Length of major ticks on y-axis
    'xtick.minor.size': 4,          # Length of minor ticks on x-axis
    'ytick.minor.size': 4,          # Length of minor ticks on y-axis
    'xtick.labelsize': f_s,         # Font Size of x-axis tick labels
    'ytick.labelsize': f_s,         # Font Size of y-axis tick labels
    'axes.labelsize': f_s,          # Font Size of axis labels
    'legend.fontsize': f_s,         # Font Size of legend
    'axes.titlesize': f_s,          # Font Size of axis titles
    'figure.titlesize': f_s,        # Font Size of figure titles
    'xtick.direction': 'in',        # X-axis ticks point inward
    'ytick.direction': 'in',        # Y-axis ticks point inward
    'xtick.top': True,              # Show ticks on top axis
    'xtick.bottom': True,           # Show ticks on bottom axis (default)
    'ytick.left': True,             # Show ticks on left axis (default)
    'ytick.right': True,            # Show ticks on right axis
    'axes.grid.which': 'both',      # Grid lines at both major and minor ticks
    'xtick.major.size': 4,          # # Optional: control tick length for inward ticks Slightly longer major ticks
    'xtick.minor.size': 2,          # Minor ticks
    # 'ytick.major.size': 8,
    # 'ytick.minor.size': 4
    # 'axes.titlepad': 20,
    # 'axes.labelpad': 20,
    
    
})

def plot_comparison(lammps_file, firecore_file, title=None, out_file=None, x_range=None, y_range=None, error_range=None, ax=None):
    """
    Loads, processes, and plots energy data from LAMMPS and FireCore files.
    """
    lammps_data = np.loadtxt(lammps_file, skiprows=1)
    firecore_data = np.loadtxt(firecore_file, skiprows=2)

    # --- Shift both datasets so their first point is (0,0) ---
    x_l_rel = lammps_data[:, 0] - lammps_data[0, 0]
    y_l_rel = lammps_data[:, 1] - lammps_data[0, 1]

    x_f_rel = firecore_data[:, 0] - firecore_data[0, 0]
    y_f_rel = firecore_data[:, 1] - firecore_data[0, 1]

    # Interpolate shifted LAMMPS energies onto shifted FireCore x-grid for difference
    lammps_x_rel_sorted_indices = np.argsort(x_l_rel)
    lammps_x_rel_sorted = x_l_rel[lammps_x_rel_sorted_indices]
    lammps_y_rel_sorted = y_l_rel[lammps_x_rel_sorted_indices]

    y_l_interp_rel = np.interp(x_f_rel, lammps_x_rel_sorted, lammps_y_rel_sorted)
    difference_data_rel = y_f_rel - y_l_interp_rel

    # --- Plotting ---
    # fig, ax1 = plt.subplots(figsize=(8, 6.5))
    if ax is None:
        fig, ax1 = plt.subplots(figsize=(8, 6.5))
    else:
        fig = ax.figure
        ax1 = ax
    
    ax1.plot(x_l_rel, y_l_rel, linewidth=2, label='LAMMPS')
    ax1.plot(x_f_rel, y_f_rel, linestyle='--', linewidth=2, label='FireCore')
    
    ax1.set_xlabel(r'Along Scan direction ($\mathrm{\AA}$)', fontsize=f_s)
    ax1.set_ylabel('Energy (eV)', fontsize=f_s)
    ax1.set_title(title, fontsize=f_s)
    if y_range: ax1.set_ylim(y_range)
    if x_range: ax1.set_xlim(x_range)
    
    ax2 = ax1.twinx()
    ax2.plot(x_f_rel, difference_data_rel, linestyle='-', color='r', label='FireCore - LAMMPS')
    ax2.set_ylabel('Difference (eV)', fontsize=f_s, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    if error_range: ax2.set_ylim(error_range)
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, frameon=False,fontsize=f_s-5, loc='upper left',ncol=2,columnspacing=0.5,handletextpad=0.3)

    
    ax1.tick_params(axis='both', labelsize=f_s)
    ax1.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    if out_file is not None:  # Only save if output file specified
        plt.savefig(out_file, dpi=600)
    if ax is None:  # Only show/close if we created the figure
        plt.show()
        plt.close(fig)
    print(f"Saved plot to {out_file}" if out_file else "Plot created")

def main():
    parser = argparse.ArgumentParser(description="Compare FireCore and LAMMPS energy data.")
    parser.add_argument("--lammps", required=True, help="Path to LAMMPS data file.")
    parser.add_argument("--firecore", required=True, help="Path to FireCore data file.")
    parser.add_argument("--title", default=None, help="Plot title.")
    parser.add_argument("--out", required=True, help="Path to save the output PNG plot.")
    parser.add_argument("--error-range", type=float, nargs=2, help="Y-axis range for the error plot.")
    parser.add_argument("--x-range", type=float, nargs=2, help="X-axis range for the plot.")
    parser.add_argument("--y-range", type=float, nargs=2, help="Y-axis range for the plot.")
    args = parser.parse_args()

    plot_comparison(
        lammps_file=args.lammps,
        firecore_file=args.firecore,
        title=args.title,
        out_file=args.out,
        error_range=args.error_range,
        x_range=args.x_range,
        y_range=args.y_range,
    )

if __name__ == "__main__":
    main()



'''
python /home/indranil/git/FireCore/doc/py/energy_comparison_1d.py \
--lammps /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan/angle26.565051177078/fixedatom27/nx20/total.dat \
--firecore /home/indranil/git/FireCore/tests/tMMFF/relax_perfect_line/dir_2.0_1.0_0.0/cons_26/PTCDA_20x20_26_total.dat \
--title "Energy Comparison" --out /home/indranil/git/FireCore/tests/tMMFF/relax_perfect_line/dir_2.0_1.0_0.0/cons_26/trial.png


python /home/indranil/git/FireCore/doc/py/energy_comparison_1d.py \
--lammps /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan/angle45/fixedatom27/nx20/total.dat \
--firecore /home/indranil/git/FireCore/tests/tMMFF/relax_perfect_line/dir_1.0_1.0_0.0/cons_26/PTCDA_20x20_26_total.dat \
--title "Energy Comparison" --out trial.png --error-range -0.01 0.01 --x-range 0 10 --y-range -0.01 0.01



'''