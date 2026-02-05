#!/usr/bin/env python3
"""
Plot Thermodynamic Integration Results

This script reads the TI.dat file produced by the thermodynamic integration
calculation and plots dE/dlambda vs lambda, then integrates to get the
total free energy change.
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend to avoid Qt issues
import matplotlib.pyplot as plt
import argparse
from scipy import integrate
import re

# Constants
kB = 8.617333262145e-5  # eV/K

def read_ti_data(filename):
    """Read thermodynamic integration data from file"""
    try:
        data = np.loadtxt(filename, comments='#')
        if data.ndim == 1:
            data = data.reshape(1, -1)
        return data
    except Exception as e:
        print(f"Error reading file {filename}: {e}")
        sys.exit(1)

def detect_columns(filename):
    """
    Parse header to determine column mapping.
    Returns a dictionary mapping 'lambda', 'dE_dlambda', 'sigma_dE', 'cumulative_FE', 'cumulative_err', 'distance' to indices.
    """
    col_map = {'lambda': 0, 'dE_dlambda': 1, 'sigma_dE': 2, 'cumulative_FE': 3, 'cumulative_err': 4, 'distance': 5}
    header_found = False
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#') and 'lambda' in line and 'dE' in line:
                    parts = line.replace('#','').split()
                    # Map standard names
                    temp_map = {}
                    for i, p in enumerate(parts):
                        if p == 'lambda': temp_map['lambda'] = i
                        elif p == 'dE/dlambda': temp_map['dE_dlambda'] = i
                        elif p == 'sigma_dE': temp_map['sigma_dE'] = i
                        elif p == 'cumulative_FE': temp_map['cumulative_FE'] = i
                        elif p == 'cumulative_err': temp_map['cumulative_err'] = i
                        elif p == 'distance': temp_map['distance'] = i

                    if 'lambda' in temp_map and 'dE_dlambda' in temp_map:
                        col_map.update(temp_map)
                        header_found = True
                    break
    except Exception as e:
         print(f"Warning: Could not parse header, using default columns: {e}")

    return col_map, header_found

def plot_ti_results(filename, output_prefix=None, N_segments=None, T=300.0, b=1.198):
    """Plot TI results and compute free energy"""
    
    col_map, header_found = detect_columns(filename)
    if header_found:
        print(f"Detected columns: {col_map}")
    else:
        print("No header detected, using default column mapping.")

    data = read_ti_data(filename)
    ncols = data.shape[1]
    
    def get_col(name):
        idx = col_map.get(name)
        if idx is not None and idx < ncols:
            return data[:, idx]
        return None

    lambda_vals = get_col('lambda')
    dE_dlambda = get_col('dE_dlambda')
    
    if lambda_vals is None or dE_dlambda is None:
        print("Error: Could not extract lambda or dE/dlambda columns.")
        sys.exit(1)
        
    errors = get_col('sigma_dE')
    cumulative_FE_read = get_col('cumulative_FE')
    cumulative_err_read = get_col('cumulative_err')
    distances = get_col('distance')

    if not header_found:
        if ncols >= 3: errors = data[:, 2]
        else: errors = None
        if ncols >= 4: cumulative_FE_read = data[:, 3]
        else: cumulative_FE_read = None
        if ncols >= 5: distances = data[:, 4]
        else: distances = None
        cumulative_err_read = None

    if errors is None:
        errors = np.zeros_like(lambda_vals)

    sort_idx = np.argsort(lambda_vals)
    lambda_vals = lambda_vals[sort_idx]
    dE_dlambda = dE_dlambda[sort_idx]
    if errors is not None: errors = errors[sort_idx]
    if distances is not None: distances = distances[sort_idx]
    if cumulative_FE_read is not None: cumulative_FE_read = cumulative_FE_read[sort_idx]
    if cumulative_err_read is not None: cumulative_err_read = cumulative_err_read[sort_idx]

    Force_from_data, dE_ref_dlambda, F_ref, k_spring = None, None, None, None
    if distances is not None and len(distances) > 1:
        dR_dlambda = np.gradient(distances, lambda_vals)
        dR_dlambda_safe = np.copy(dR_dlambda)
        dR_dlambda_safe[np.abs(dR_dlambda_safe) < 1e-9] = 1e-9
        Force_from_data = dE_dlambda / dR_dlambda_safe
        
        if N_segments is not None:
            k_spring = kB * T / (N_segments * b**2)
            F_ref = k_spring * distances
            dE_ref_dlambda = F_ref * dR_dlambda

    cumulative_FE_calc = np.zeros_like(lambda_vals)
    cumulative_error_calc = np.zeros_like(lambda_vals)
    if errors is not None:
        for i in range(1, len(lambda_vals)):
            dx = lambda_vals[i] - lambda_vals[i-1]
            cumulative_FE_calc[i] = cumulative_FE_calc[i-1] + 0.5 * (dE_dlambda[i] + dE_dlambda[i-1]) * dx
            var_contribution = (0.5 * dx)**2 * (errors[i]**2 + errors[i-1]**2)
            cumulative_error_calc[i] = np.sqrt(cumulative_error_calc[i-1]**2 + var_contribution)

    cumulative_FE = cumulative_FE_read if cumulative_FE_read is not None else cumulative_FE_calc
    cumulative_error = cumulative_err_read if cumulative_err_read is not None else cumulative_error_calc
    
    delta_F_trapz = np.trapz(dE_dlambda, lambda_vals)
    delta_F_simps = integrate.simpson(y=dE_dlambda, x=lambda_vals) if len(lambda_vals) >= 3 else delta_F_trapz

    # Create figure
    num_plots = 4 if distances is not None else 2
    fig, axes = plt.subplots(num_plots, 1, figsize=(10, 5 * num_plots), sharex=True)
    
    # Plot 1: dE/dlambda
    ax1 = axes[0]
    ax1.plot(lambda_vals, dE_dlambda, 'o-', linewidth=2, markersize=6, label='dE/dλ (data)')
    if errors is not None:
        ax1.fill_between(lambda_vals, dE_dlambda - errors, dE_dlambda + errors, alpha=0.3)
    if dE_ref_dlambda is not None:
        ax1.plot(lambda_vals, dE_ref_dlambda, '--', color='purple', linewidth=2, label='dE/dλ (ref)')
    ax1.set_ylabel('dE/dλ [eV]')
    ax1.set_title('Thermodynamic Integration Analysis')
    ax1.grid(True, alpha=0.5)
    ax1.legend()

    textstr = f'ΔF (trapz) = {delta_F_trapz:.4f} eV\nΔF (Simpson) = {delta_F_simps:.4f} eV'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=12, verticalalignment='top', bbox=props)

    if distances is not None:
        # Plot 2: Force
        ax2 = axes[1]
        ax2.plot(lambda_vals, Force_from_data, 'o-', color='orange', markersize=6, label='Force (from data)')
        if F_ref is not None:
            ax2.plot(lambda_vals, F_ref, '--', color='purple', label='Force (ref spring)')
        ax2.set_ylabel('Force [eV/Å]')
        ax2.grid(True, alpha=0.5)
        ax2.legend()

        # Plot 3: Cumulative FE
        ax3 = axes[2]
        ax3.plot(lambda_vals, cumulative_FE, 'd-', color='green', markersize=6, label='Cumulative ΔF (data)')
        if cumulative_error is not None:
            ax3.fill_between(lambda_vals, cumulative_FE - cumulative_error, cumulative_FE + cumulative_error, alpha=0.3, color='green')
        if F_ref is not None:
            FE_ref = 0.5 * k_spring * distances**2
            FE_ref_shifted = FE_ref - FE_ref[0] + cumulative_FE[0]
            ax3.plot(lambda_vals, FE_ref_shifted, '--', color='purple', label='ΔF(λ) (ref spring)')
        ax3.set_ylabel('ΔF(λ) [eV]')
        ax3.grid(True, alpha=0.5)
        ax3.legend()

        # Plot 4: Distance
        ax4 = axes[3]
        ax4.plot(lambda_vals, distances, 's-', color='red', markersize=6, label='Si-Si Distance')
        ax4.set_xlabel('λ')
        ax4.set_ylabel('Distance [Å]')
        ax4.grid(True, alpha=0.5)
        ax4.legend()
    else:
        # Plot 2: Cumulative FE (no distance)
        ax2 = axes[1]
        ax2.plot(lambda_vals, cumulative_FE, 'd-', color='green', markersize=6, label='Cumulative ΔF')
        if cumulative_error is not None:
            ax2.fill_between(lambda_vals, cumulative_FE - cumulative_error, cumulative_FE + cumulative_error, alpha=0.3, color='green')
        ax2.set_xlabel('λ')
        ax2.set_ylabel('ΔF(λ) [eV]')
        ax2.grid(True, alpha=0.5)
        ax2.legend()
    
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.suptitle(f"TI Analysis for {filename}", fontsize=16)

    if output_prefix is None:
        output_prefix = filename.replace('.dat', '')
    
    output_png = f"{output_prefix}_TI_plot.png"
    plt.savefig(output_png, dpi=300)
    print(f"\nPlot saved to: {output_png}")
    
    print("\n" + "="*60)
    print("  Thermodynamic Integration Summary")
    print("="*60)
    print(f"  Number of λ windows: {len(lambda_vals)}")
    print(f"  Free energy change (trapezoidal): {delta_F_trapz:.6f} eV")
    print(f"  Free energy change (Simpson):     {delta_F_simps:.6f} eV")
    if distances is not None:
        print(f"  Distance range: {distances[0]:.3f} → {distances[-1]:.3f} Å")
    if N_segments is not None:
         print(f"  Entropic Spring Reference (N={N_segments}, T={T}K, b={b}A)")
    print("="*60)

    return delta_F_trapz, delta_F_simps

def main():
    parser = argparse.ArgumentParser(description="Plot Thermodynamic Integration Results")
    parser.add_argument("--input", type=str, default="entropic_spring_20_TI.dat", 
                       help="Input TI data file")
    parser.add_argument("--output", type=str, default=None, 
                       help="Output prefix for plots (default: derived from input)")
    parser.add_argument("--N", type=int, default=None,
                       help="Number of segments (for entropic spring reference). If None, tries to guess from filename.")
    parser.add_argument("--T", type=float, default=300.0,
                       help="Temperature in Kelvin (default: 300.0)")
    parser.add_argument("--b", type=float, default=1.198,
                       help="Segment length in Angstroms (default: 1.198)")
    
    args = parser.parse_args()
    
    # Try to guess N from filename if not provided
    if args.N is None:
        match = re.search(r'_(\d+)_', args.input)
        if match:
             args.N = int(match.group(1))
             print(f"Guessed N={args.N} from filename.")
    
    print(f"Reading TI data from: {args.input}")
    plot_ti_results(args.input, args.output, N_segments=args.N, T=args.T, b=args.b)

if __name__ == "__main__":
    main()
