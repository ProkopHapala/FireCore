#!/usr/bin/env python3

# Environment setup and compilation section
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import subprocess
import time
import shutil

# Regular imports
sys.path.append("../../")
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

# Import the all_scan module for scan functions
from all_scan import run_scan


def generate_scan(molecule, substrate, output_dir, scan_type='total', scan_params=None, skip_init=False):
    """
    Generate data for a specific potential component
    
    Args:
        molecule (str): Path to the molecule xyz file
        substrate (str): Path to the substrate xyz file
        output_dir (str): Directory to save output files
        scan_type (str): Type of scan to perform (one of the keys in SCAN_TYPES)
        scan_params (dict): Parameters for the scan
        skip_init (bool): If True, skip MMFF initialization (assumes it's already initialized)
    
    Returns:
        bool: True if scan completed successfully
    """
    # Validate scan_type is a string and one of the allowed types
    if not isinstance(scan_type, str):
        print(f"Error: scan_type must be a string, not {type(scan_type)}")
        return False
        
    # Only allow total, morse, and coulomb scan types
    if scan_type.lower() not in ['total', 'morse', 'coulomb']:
        print(f"Error: Invalid scan type '{scan_type}'. Only 'total', 'morse', and 'coulomb' are supported.")
        return False
        
    # Call the run_scan function from all_scan.py
    try:
        return run_scan(molecule, substrate, output_dir, scan_type, scan_params, skip_init=skip_init)
    except ValueError as e:
        print(f"Error: {e}")
        return False
    except Exception as e:
        print(f"Unexpected error during {scan_type} scan: {e}")
        return False


def plot_multiple_comparisons(data_pairs, titles, output_filename, x_range=None, y_ranges=None, error_ranges=None):
    """
    Create a figure with multiple comparison plots in a single row.
    """
    # Create the main figure
    n_plots = len(data_pairs)
    fig = plt.figure(figsize=(6*n_plots, 6))
    
    # Default ranges if not provided
    if y_ranges is None:
        y_ranges = [None] * n_plots
    if error_ranges is None:
        error_ranges = [None] * n_plots
    
    # Process each plot
    for i in range(n_plots):
        # Get data for this pair
        lammps_file, firecore_file = data_pairs[i]
        title = titles[i]
        y_range = y_ranges[i] if y_ranges else None
        error_range = error_ranges[i] if error_ranges else None
        
        # Create subplot
        ax1 = fig.add_subplot(1, n_plots, i+1)
        
        # Load data
        lammps_data = np.loadtxt(lammps_file)
        firecore_data = np.loadtxt(firecore_file)
        
        # Calculate difference
        difference_data = firecore_data[:, 1] - lammps_data[:, 1]
        
        # Plot on primary axis
        ax1.plot(lammps_data[:, 0], lammps_data[:, 1], linewidth=2, label='LAMMPS')
        ax1.plot(firecore_data[:, 0], firecore_data[:, 1], marker='o', markersize=8, 
                linestyle='', label='FireCore', markerfacecolor='none')
        
        # Set labels and title
        ax1.set_xlabel(r'Z ($\mathrm{\AA}$)', fontsize=12)
        if i == 0:  # Only add y-label to the first subplot
            ax1.set_ylabel('Energy (eV)', fontsize=12)
        ax1.set_title(title, fontsize=14)
        
        # Apply range limits
        if y_range:
            ax1.set_ylim(y_range)
        if x_range:
            ax1.set_xlim(x_range)
        
        # Secondary axis for difference data
        ax2 = ax1.twinx()
        ax2.plot(firecore_data[:, 0], difference_data, marker='*', markersize=5, 
                linestyle='', color='r', label='Difference')
        
        # Set secondary y-axis label only for the last subplot
        if i == n_plots - 1:
            ax2.set_ylabel('Difference (eV)', fontsize=12, color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        
        # Set error range
        if error_range:
            ax2.set_ylim(error_range)
        
        # Add legend
        lines1, labels1 = ax1.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax1.legend(lines1 + lines2, labels1 + labels2, fontsize=10, loc='best')
        
        # Grid
        ax1.grid(True, linestyle='--', alpha=0.7)
    
    # Adjust layout and save
    plt.tight_layout()
    if output_filename:
        plt.savefig(output_filename, dpi=600, bbox_inches='tight')
    
    return fig

def compare_with_lammps(molecule, lammps_dir, firecore_dir, output_dir):
    """
    Create comparison plots between FireCore and LAMMPS data
    
    Args:
        molecule (str): Molecule name for file naming
        lammps_dir (str): Directory containing LAMMPS data files
        firecore_dir (str): Directory containing FireCore data files
        output_dir (str): Directory to save comparison plots
    """
    # Extract molecule name for file naming
    mol_name = os.path.basename(molecule)
    
    # Define data pairs for Z-direction scan
    data_pairs = [
        (f'{lammps_dir}/morse.dat', f'{firecore_dir}/{mol_name}_morse.dat'),
        (f'{lammps_dir}/coul.dat', f'{firecore_dir}/{mol_name}_coul.dat'),
        (f'{lammps_dir}/total.dat', f'{firecore_dir}/{mol_name}_total.dat')
    ]

    # Plot titles
    titles = [
        "Morse Potential",
        "Coulomb Potential",
        "Total Potential"
    ]

    # Create and save the comparison figure
    fig = plot_multiple_comparisons(data_pairs, titles, 
                                   f"{output_dir}/{mol_name}_comparison.png")
    
    # Display the figure
    plt.show()
    
    return True

# Main execution - Examples
if __name__ == "__main__":
    import argparse
    
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Generate molecular scan data and compare with LAMMPS')

    parser.add_argument('--molecule', default='data/xyz/old_mol_old_sub_PTCDA', help='Path to molecule xyz file')
    parser.add_argument('--substrate', default='data/xyz/Na_0.9_Cl_-0.9', help='Path to substrate xyz file')
    # parser.add_argument('--molecule', default='data/xyz/old_PTCDA', help='Path to molecule xyz file')
    # parser.add_argument('--substrate', default='data/xyz/old_NaCl', help='Path to substrate xyz file')

    parser.add_argument('--output-dir', default='PTCDA_data_step_0.05', help='Output directory for scan data')
    parser.add_argument('--scan-types', nargs='+', 
                        choices=['total', 'morse', 'coulomb', 'all'], 
                        default=['all'], help='Types of scans to run')
    parser.add_argument('--compare', action='store_true', help='Compare with LAMMPS data after generating scans')
    # parser.add_argument('--lammps-dir', default='/home/indranil/Documents/Project_1/Lammps_backup15.05/3-rigid_PES/1-zscan_on_Na', help='Directory containing LAMMPS data files')
    parser.add_argument('--lammps-dir', default='/home/indranil/Documents/Project_1/Lammps/1-rigid_zscan_step_0.05', help='Directory containing LAMMPS data files')
    # parser.add_argument('--lammps-dir', default='/home/indranil/Documents/Project_1/Lammps/1-rigid_zscan_step_0.2', help='Directory containing LAMMPS data files')
    
    parser.add_argument('--nscan', type=int, default=250, help='Number of scan points')
    parser.add_argument('--span-min', type=float, default=2.6, help='Minimum scan distance')
    parser.add_argument('--span-max', type=float, default=15.1, help='Maximum scan distance')
    
    args = parser.parse_args()
    
    # Process the requested scan types
    if 'all' in args.scan_types:
        run_types = ['total', 'morse', 'coulomb']  # Only include these three types
    else:
        # Filter to only include valid scan types
        run_types = [t for t in args.scan_types if t in ['total', 'morse', 'coulomb']]
        
    # Ensure we have at least one valid scan type
    if not run_types:
        print("Error: No valid scan types specified. Must be one of: 'total', 'morse', 'coulomb', or 'all'")
        sys.exit(1)
    
    # Set scan parameters
    scan_params = {
        "nscan": args.nscan,
        "span": (args.span_min, args.span_max),
        "dir": (0.0, 0.0, 1.0),  # Default direction along z-axis
        "p0": (0, 0, 0),        # Default starting position
    }
    
    # Print banner
    print("=" * 80)
    print(f"Generating {', '.join(run_types)} scans for {args.molecule}")
    print(f"Substrate: {args.substrate}")
    print(f"Output directory: {args.output_dir}")
    print(f"Scan parameters: {scan_params}")
    print("=" * 80)
    
    # Generate the scans
    results = {}
    
    # If 'all' was specified and there are multiple scan types,
    # and this is not a recursive call, handle each scan type separately
    if 'all' in args.scan_types and len(run_types) > 1:
        # For each scan type, call the script again with just that scan type
        for scan_type in run_types:
            print(f"\n=====================================================================")
            print(f"RUNNING {scan_type.upper()} SCAN")
            print(f"====================================================================\n")
            
            # Run the scan for this type individually
            result = generate_scan(
                args.molecule, 
                args.substrate, 
                args.output_dir, 
                scan_params=scan_params,
                scan_type=scan_type  # This must be a string, not a list
            )
            
            results[scan_type] = result
    else:
        # Process a single scan type - this is either a direct call with a single type
        # or a recursive call from the 'all' handler above
        for scan_type in run_types:
            print(f"\n=====================================================================")
            print(f"RUNNING {scan_type.upper()} SCAN")
            print(f"====================================================================\n")
            
            # Run the scan for this type (passing a single string, not a list)
            result = generate_scan(
                args.molecule, 
                args.substrate, 
                args.output_dir, 
                scan_params=scan_params,
                scan_type=scan_type  # This must be a string, not a list
            )
            
            results[scan_type] = result
    
    print(f"\n=====================================================================")
    print(f"ALL REQUESTED SCANS COMPLETED")
    print(f"Results:") 
    for scan_type, success in results.items():
        status = "SUCCESS" if success else "FAILED"
        print(f"  {scan_type}: {status}")
    print(f"====================================================================\n")

    
    # Compare with LAMMPS data if requested
    if args.compare and results and all(results.values()):
        print("\nComparing with LAMMPS data...")
        mol_name = os.path.basename(args.molecule)
        compare_with_lammps(mol_name, args.lammps_dir, args.output_dir, args.output_dir)
        print(f"Comparison plots saved to {args.output_dir}")
