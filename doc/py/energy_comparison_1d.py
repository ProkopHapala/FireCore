import argparse
import numpy as np
import matplotlib.pyplot as plt

def plot_comparison(lammps_file, firecore_file, title, out_file, x_range=None, y_range=None, error_range=None):
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
    fig, ax1 = plt.subplots(figsize=(8, 6.5))
    
    ax1.plot(x_l_rel, y_l_rel, linewidth=2, label='LAMMPS')
    ax1.plot(x_f_rel, y_f_rel, linestyle='--', linewidth=2, label='FireCore')
    
    ax1.set_xlabel(r'Relative Z ($\mathrm{\AA}$)', fontsize=14)
    ax1.set_ylabel('Relative Energy (eV)', fontsize=14)
    ax1.set_title(title, fontsize=16)
    if y_range: ax1.set_ylim(y_range)
    if x_range: ax1.set_xlim(x_range)
    
    ax2 = ax1.twinx()
    ax2.plot(x_f_rel, difference_data_rel, linestyle='-', color='r', label='FireCore - LAMMPS')
    ax2.set_ylabel('Difference (eV)', fontsize=14, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    if error_range: ax2.set_ylim(error_range)
    
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=12, loc='best')
    
    ax1.tick_params(axis='both', labelsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    plt.savefig(out_file, dpi=300)
    plt.show()
    plt.close(fig)
    print(f"Saved plot to {out_file}")

def main():
    parser = argparse.ArgumentParser(description="Compare FireCore and LAMMPS energy data.")
    parser.add_argument("--lammps", required=True, help="Path to LAMMPS data file.")
    parser.add_argument("--firecore", required=True, help="Path to FireCore data file.")
    parser.add_argument("--title", default="Energy Comparison", help="Plot title.")
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