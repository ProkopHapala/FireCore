from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable, Sequence

import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Re-use the generic visualiser that already exists in FireCore
# from plot_utils import MoleculeTrajectoryVisualizer
from scipy.spatial import ConvexHull
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




def find_local_extrema(y_values):
    """Find indices of local minima and maxima in a 1D array"""
    minima = []
    maxima = []
    for i in range(1, len(y_values)-1):
        if y_values[i] < y_values[i-1] and y_values[i] < y_values[i+1]:
            minima.append(i)
        elif y_values[i] > y_values[i-1] and y_values[i] > y_values[i+1]:
            maxima.append(i)
    return np.array(minima), np.array(maxima)

def plot_single_comparison_error(lammps_file, firecore_file, title, x_range=None, y_range=None, error_range=None):
    lammps_data = np.loadtxt(lammps_file)
    firecore_data = np.loadtxt(firecore_file)

    # Find all local extrema in FireCore data
    minima_indices, maxima_indices = find_local_extrema(firecore_data[:, 1])
    minima_points = firecore_data[minima_indices]
    maxima_points = firecore_data[maxima_indices]
    
    difference_data = firecore_data[:, 1] - lammps_data[:, 1]

    fig, ax1 = plt.subplots(figsize=(8.0, 6.5))
 
    # Plot main data
    ax1.plot(firecore_data[:, 0], firecore_data[:, 1], '-g', lw=1.5, label='FireCore') 
    ax1.plot(lammps_data[:, 0], lammps_data[:, 1], ':g', lw=3.0, label='LAMMPS')
    
    # Highlight and label minima (blue)
    minima_points = [[ 2.8,-0.94101623],[ 8.4,-0.21400479],[10.1,-0.46179296],[12.7,-0.09324949],[1.279999999999999893e+01 ,-4.507061617454444225e-01],[2.000000000000000000e+01-0.3 ,-5.971518845582224344e-05]]
    cmap = plt.cm.tab10
    point_colors = [cmap(i) for i in range(10)] 
    for i, (x, y) in enumerate(minima_points):
        ax1.scatter(x, y, color=point_colors[i], s=100, zorder=5)
        ax1.text(x, y, i+1, fontsize=f_s, ha='center', va='bottom',color=point_colors[i])

    # for x, y in minima_points:
    #     ax1.scatter(x, y, color='blue', s=100, zorder=5)
    #     # ax1.text(x, y, f'Min: ({x:.2f}, {y:.4f})', fontsize=10, ha='center', va='bottom',bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
    # # Highlight and label maxima (red)
    # for x, y in maxima_points:
    #     ax1.scatter(x, y, color='red', s=100, zorder=5)
    #     # ax1.text(x, y, f'Max: ({x:.2f}, {y:.4f})',fontsize=10, ha='center', va='top',bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
    # print(minima_points)
    # print(maxima_points)

    ax1.set_xlabel(r'Z (\AA)', fontsize=f_s)
    ax1.set_ylabel('Energy (eV)', fontsize=f_s)
    if title: ax1.set_title(title, fontsize=f_s)
    
    if y_range: ax1.set_ylim(y_range)
    if x_range: ax1.set_xlim(x_range)
    
    ax2 = ax1.twinx()
    ax2.plot(firecore_data[:, 0], difference_data*1000, linestyle='-',  
             markerfacecolor='none', color='r', lw=3.0, alpha=0.2, 
             label='Energy Difference')
    
    ax2.set_ylabel('Energy difference (meV)', fontsize=f_s, color='red')
    ax2.tick_params(axis='y', labelcolor='red')
    
    if error_range: ax2.set_ylim(error_range)
    
    # Add legends
    lines1, labels1 = ax1.get_legend_handles_labels()
    ax1.legend(lines1, labels1, fontsize=f_s, loc='best', frameon=False)
    
    # Formatting
    ax1.tick_params(axis='both', labelsize=f_s)
    ax1.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    
    return plt

lammps_file = '/home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan/total.dat'  
firecore_file = '/home/indranil/git/FireCore/tests/tMMFF/PTCDA_data_trial_1d_relax_z/old_mol_old_sub_PTCDA_total.dat' 

plt_c_na = plot_single_comparison_error(lammps_file, firecore_file, None, x_range=(0, 20),y_range=(-1, 0.2),error_range=(-1e-0, 2e-0)) 

plt.tight_layout()
# savepath='/home/indranil/Documents/Project_1/results/relaxed_1d_z.png'
# plt.savefig(savepath, format='png', dpi=600)


plt.show() 