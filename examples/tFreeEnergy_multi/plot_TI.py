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
    
    # Detect columns
    col_map, header_found = detect_columns(filename)
    if header_found:
        print(f"Detected columns: {col_map}")
    else:
        print("No header detected, using default column mapping.")

    # Read data
    data = read_ti_data(filename)
    
    # Check bounds
    ncols = data.shape[1]
    
    # Safe extraction helper
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
    cumulative_FE = get_col('cumulative_FE')
    cumulative_err_read = get_col('cumulative_err')
    distances = get_col('distance')

    # If using default mapping but data has fewer columns, we might have gotten None or wrong data.
    # The default mapping assumes 6 columns.
    # If we didn't find a header, we should apply fallback logic based on shape for backward compatibility.
    if not header_found:
        if ncols < 6:
            # Old format: lambda, dE/dlambda, sigma_dE, [cumulative_FE], [distance]
            # Reset based on observed shape
            if ncols >= 3: errors = data[:, 2]
            else: errors = None

            if ncols >= 4: cumulative_FE = data[:, 3]
            else: cumulative_FE = None

            if ncols >= 5: distances = data[:, 4]
            else: distances = None

            cumulative_err_read = None

    # Sort by lambda (in case data is not sorted)
    sort_idx = np.argsort(lambda_vals)
    lambda_vals = lambda_vals[sort_idx]
    dE_dlambda = dE_dlambda[sort_idx]
    if errors is not None:
        errors = errors[sort_idx]
    if distances is not None:
        distances = distances[sort_idx]
    if cumulative_FE is not None:
        cumulative_FE = cumulative_FE[sort_idx]
    if cumulative_err_read is not None:
        cumulative_err_read = cumulative_err_read[sort_idx]
    
    # Integrate to get free energy
    # Using trapezoidal rule
    delta_F_trapz = np.trapz(dE_dlambda, lambda_vals)
    
    # Using Simpson's rule if we have enough points
    if len(lambda_vals) >= 3:
        delta_F_simps = integrate.simpson(y=dE_dlambda, x=lambda_vals)
    else:
        delta_F_simps = delta_F_trapz
    
    # Calculate cumulative error if not present in file
    # Use the cumulative error from file if available, otherwise compute from local errors
    cumulative_error_calc = None
    if cumulative_err_read is not None:
        cumulative_error_calc = cumulative_err_read
    elif errors is not None:
        cumulative_error_calc = np.zeros_like(lambda_vals)
        for i in range(1, len(lambda_vals)):
            dx = lambda_vals[i] - lambda_vals[i-1]
            var_contribution = (0.5 * dx)**2 * (errors[i]**2 + errors[i-1]**2)
            cumulative_error_calc[i] = np.sqrt(cumulative_error_calc[i-1]**2 + var_contribution)
            
    # Create figure with subplots (3 rows)
    fig, axes = plt.subplots(3, 1, figsize=(10, 15))
    
    # Plot 1: dE/dlambda vs lambda
    ax1 = axes[0]
    ax1.plot(lambda_vals, dE_dlambda, 'o-', linewidth=2, markersize=8, label='dE/dλ')
    if errors is not None and np.any(errors > 0):
        ax1.fill_between(lambda_vals, dE_dlambda - errors, dE_dlambda + errors, 
                         alpha=0.3, label='Error estimate')
    
    # Entropic Spring Reference (Force)
    if distances is not None and N_segments is not None:
        # F = (3kT / Nb^2) * R
        k_spring = kB * T / (N_segments * b**2)
        F_ref = k_spring * distances
        ax1.plot(lambda_vals, F_ref, '--', color='purple', linewidth=2, 
                 label=f'Ref Force (N={N_segments}, T={T}K)')
                 
    ax1.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax1.set_xlabel('λ', fontsize=14)
    ax1.set_ylabel('dE/dλ [eV]', fontsize=14)
    ax1.set_title('Thermodynamic Integration: Force Derivative', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=12)
    
    # Add text box with integration results
    textstr = f'ΔF (trapz) = {delta_F_trapz:.4f} eV\nΔF (Simpson) = {delta_F_simps:.4f} eV'
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax1.text(0.05, 0.95, textstr, transform=ax1.transAxes, fontsize=12,
             verticalalignment='top', bbox=props)
             
    # Plot 2: Cumulative Free Energy vs lambda
    ax2 = axes[1]
    # If we read cumulative FE from file, use it, otherwise compute it
    if cumulative_FE is None:
        cumulative_FE_calc = np.zeros_like(lambda_vals)
        for i in range(1, len(lambda_vals)):
             cumulative_FE_calc[i] = cumulative_FE_calc[i-1] + 0.5 * (dE_dlambda[i] + dE_dlambda[i-1]) * (lambda_vals[i] - lambda_vals[i-1])
        cumulative_FE = cumulative_FE_calc
             
    ax2.plot(lambda_vals, cumulative_FE, 'd-', linewidth=2, markersize=8, 
             color='green', label='Cumulative ΔF')
             
    if cumulative_error_calc is not None:
         ax2.fill_between(lambda_vals, cumulative_FE - cumulative_error_calc, cumulative_FE + cumulative_error_calc,
                          alpha=0.3, color='green', label='Error estimate')

    # Entropic Spring Reference (Free Energy)
    if distances is not None and N_segments is not None:
        # FE = 0.5 * k * R^2 = 0.5 * (3kT / Nb^2) * R^2
        k_spring = kB * T / (N_segments * b**2)
        FE_ref = 0.5 * k_spring * distances**2
        
        # Shift to match start
        FE_ref_shifted = FE_ref - FE_ref[0] + cumulative_FE[0]
        
        ax2.plot(lambda_vals, FE_ref_shifted, '--', color='purple', linewidth=2,
                 label=f'Ref FE (shifted)')
                          
    ax2.set_xlabel('λ', fontsize=14)
    ax2.set_ylabel('ΔF(λ) [eV]', fontsize=14)
    ax2.set_title('Cumulative Free Energy', fontsize=16)
    ax2.grid(True, alpha=0.3)
    ax2.legend(fontsize=12)
    
    # Plot 3: Distance vs lambda (if available)
    ax3 = axes[2]
    if distances is not None:
        # Check if we should use distance or error analysis in the third plot.
        # Original code put distance here if available.
        ax3.plot(lambda_vals, distances, 's-', linewidth=2, markersize=8, 
                color='red', label='Distance between constrained atoms')
        ax3.set_xlabel('λ', fontsize=14)
        ax3.set_ylabel('Distance [Å]', fontsize=14)
        ax3.set_title('Constraint Distance vs λ', fontsize=16)
        ax3.grid(True, alpha=0.3)
        ax3.legend(fontsize=12)
    else:
        ax3.text(0.5, 0.5, 'No distance data available', 
                transform=ax3.transAxes, ha='center', va='center', fontsize=14)
        ax3.axis('off')
    
    plt.tight_layout()
    
    # Save figure
    if output_prefix is None:
        output_prefix = filename.replace('_TI.dat', '').replace('.dat', '')
    
    output_png = f"{output_prefix}_TI_plot.png"
    output_pdf = f"{output_prefix}_TI_plot.pdf"
    
    plt.savefig(output_png, dpi=300, bbox_inches='tight')
    plt.savefig(output_pdf, bbox_inches='tight')
    print(f"\nPlots saved to:")
    print(f"  {output_png}")
    print(f"  {output_pdf}")
    
    # Print summary
    print(f"\n{'=' * 60}")
    print(f"  Thermodynamic Integration Summary")
    print(f"{'=' * 60}")
    print(f"  Number of λ windows: {len(lambda_vals)}")
    print(f"  λ range: {lambda_vals[0]:.4f} → {lambda_vals[-1]:.4f}")
    print(f"  Free energy change (trapezoidal): {delta_F_trapz:.6f} eV")
    print(f"  Free energy change (Simpson):     {delta_F_simps:.6f} eV")
    if N_segments is not None:
         print(f"  Entropic Spring Reference (N={N_segments}, T={T}K, b={b}A)") 
    print(f"{'=' * 60}\n")
    
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
