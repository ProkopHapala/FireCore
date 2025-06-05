#!/usr/bin/env python3

# Environment setup and compilation section
import sys
import numpy as np
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import subprocess
import time
import shutil

# Regular imports
sys.path.append("../../")
from all_scan import run_scan, scanPlot2D_uff, scanPlot_uff
from pyBall import atomicUtils as au
from pyBall import MMFF as mmff

# Import the all_scan module for scan functions


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
    return run_scan(molecule, substrate, output_dir, scan_type, scan_params, skip_init=skip_init)
    
    


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





def preprocess_lammps_data(input_file, output_file):
    """Remove blank lines from LAMMPS data file"""
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            if line.strip():  # Check if the line is not empty after stripping whitespace
                f_out.write(line)

def compare_2d_with_lammps(molecule, lammps_dir, firecore_dir, output_dir, scan_type):
    """Compare 2D scan results between FireCore and LAMMPS
    
    Args:
        molecule (str): Name of the molecule for file naming
        lammps_dir (str): Directory containing LAMMPS data
        firecore_dir (str): Directory containing FireCore data
        output_dir (str): Directory to save comparison plots
        scan_type (str): Type of scan (total, morse, coulomb)
    """
    # Map scan types to file names (both LAMMPS and FireCore use 'coul')
    file_name = {'morse': 'morse', 'coulomb': 'coul', 'total': 'total'}[scan_type]
    
    # Define file paths
    lammps_raw = os.path.join(lammps_dir, f'{file_name}.dat')
    lammps_processed = os.path.join(lammps_dir, f'{file_name}_no_blank.dat')
    firecore_file = os.path.join(firecore_dir, f'{molecule}_{file_name}_2d.dat')
    
    # Ensure files exist (this should already be checked by the caller)
    assert os.path.exists(lammps_raw), f"LAMMPS data file not found: {lammps_raw}"
    assert os.path.exists(firecore_file), f"FireCore data file not found: {firecore_file}"
    
    # Preprocess LAMMPS data to remove blank lines
    preprocess_lammps_data(lammps_raw, lammps_processed)
    
    # Load data
    lammps_data = np.loadtxt(lammps_processed)
    firecore_data = np.loadtxt(firecore_file, skiprows=4)  # Skip header rows
    
    # Extract coordinates and energies
    x_lammps, y_lammps, e_lammps = lammps_data.T
    x_firecore, y_firecore, e_firecore = firecore_data.T
    
    # Normalize FireCore energies to start at 0
    firecore_min = np.min(e_firecore)
    e_firecore -= firecore_min
    print(f"FireCore {scan_type} data range: {firecore_min:.6f} to {np.max(e_firecore):.6f} eV")
    
    # Determine grid dimensions
    nx = int(np.sqrt(len(x_firecore)))  # Assuming square grid
    ny = nx
    
    # Print coordinate ranges for debugging
    print(f"X range: {np.min(x_firecore):.12f} to {np.max(x_firecore):.12f}")
    print(f"Y range: {np.min(y_firecore):.12f} to {np.max(y_firecore):.12f}")
    
    # Reshape data into grids (transpose to match coordinate convention)
    e_firecore_grid = e_firecore.reshape(nx, ny).T
    e_lammps_grid = e_lammps.reshape(nx, ny).T
    
    # Store grid coordinates for plotting
    x_grid = x_firecore.reshape(nx, ny).T
    y_grid = y_firecore.reshape(nx, ny).T
    
    # Print detailed energy statistics
    print(f"Grid dimensions: {nx}x{ny}")
    print(f"FireCore {scan_type} energies:")
    print(f"  Range: {np.min(e_firecore_grid):.12f} to {np.max(e_firecore_grid):.12f} eV")
    print(f"  Mean: {np.mean(e_firecore_grid):.12f} eV")
    print(f"LAMMPS {scan_type} energies:")
    print(f"  Range: {np.min(e_lammps_grid):.12f} to {np.max(e_lammps_grid):.12f} eV")
    print(f"  Mean: {np.mean(e_lammps_grid):.12f} eV")
    
    # Calculate difference with full precision
    difference_grid = np.subtract(e_firecore_grid, e_lammps_grid, dtype=np.float64)
    print(f"Difference statistics:")
    print(f"  Range: {np.min(difference_grid):.12f} to {np.max(difference_grid):.12f} eV")
    print(f"  Mean: {np.mean(difference_grid):.12f} eV")
    print(f"  RMS: {np.sqrt(np.mean(np.square(difference_grid))):.12f} eV")
    
    # Set common scale for plots
    vmin = min(np.min(e_lammps_grid), np.min(e_firecore_grid))
    vmax = max(np.max(e_lammps_grid), np.max(e_firecore_grid))
    extent = (x_lammps.min(), x_lammps.max(), y_lammps.min(), y_lammps.max())
    
    # Create figure and subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Plot LAMMPS data
    im1 = axes[0].imshow(e_lammps_grid, extent=extent, origin='lower', 
                        cmap='viridis', vmin=vmin, vmax=vmax)
    axes[0].set_title('LAMMPS')
    fig.colorbar(im1, ax=axes[0], label='Energy')
    
    # Plot FireCore data
    im2 = axes[1].imshow(e_firecore_grid, extent=extent, origin='lower', 
                        cmap='viridis', vmin=vmin, vmax=vmax)
    axes[1].set_title('FireCore')
    fig.colorbar(im2, ax=axes[1], label='Energy')
    
    # Plot difference
    im3 = axes[2].imshow(difference_grid, extent=extent, origin='lower', 
                        cmap='coolwarm')
    axes[2].set_title('Difference (FireCore - LAMMPS)')
    fig.colorbar(im3, ax=axes[2], label='Energy Difference')
    
    # Set title and save
    fig.suptitle(f"{scan_type.capitalize()} potential at Z=3.3 Å", fontsize=22)
    plt.tight_layout()
    
    # Save plot
    output_file = os.path.join(output_dir, f'{molecule}_{scan_type}_2d_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

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
    parser.add_argument('--lammps-1d-dir', type=str, help='Directory containing LAMMPS 1D scan data for comparison')
    parser.add_argument('--lammps-2d-dir', type=str, help='Directory containing LAMMPS 2D scan data for comparison')
    
    parser.add_argument('--scan-mode', default='1d', help="Scan mode: '1d' for one-dimensional scan, '2d' for two-dimensional scan.")
    parser.add_argument('--scan-dir', default='0,0,1', help='Direction vector for 1d scan. Used if --scan-mode is 1d')
    parser.add_argument('--scan-dir1', default='1,0,0', help='Direction vector for 2d scan. Required if --scan-mode is 2d')
    parser.add_argument('--scan-dir2', default='0,1,0', help='Second direction vector for 2d scan. Required if --scan-mode is 2d')
    parser.add_argument('--scan-origin', default='0,0,0', help='Scan origin for 2d scan. Required if --scan-mode is 2d')
    
    parser.add_argument('--nscan', type=int, default=250, help='Number of scan points for 1d scan')
    parser.add_argument('--span-min', type=float, default=2.6, help='Minimum scan distance for 1d scan')
    parser.add_argument('--span-max', type=float, default=15.1, help='Maximum scan distance for 1d scan')

    parser.add_argument('--nscan-1', type=int, default=250, help='Number of scan points along x for 2d scan')
    parser.add_argument('--span-min-1', type=float, default=2.6, help='Minimum scan distance along x for 2d scan')
    parser.add_argument('--span-max-1', type=float, default=15.1, help='Maximum scan distance along x for 2d scan')
    parser.add_argument('--nscan-2', type=int, default=250, help='Number of scan points along y for 2d scan')
    parser.add_argument('--span-min-2', type=float, default=2.6, help='Minimum scan distance along y for 2d scan')
    parser.add_argument('--span-max-2', type=float, default=15.1, help='Maximum scan distance along y for 2d scan')
    
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
    
    if args.scan_mode == '2d':
        # Validate required 2d scan parameters
        if not args.scan_dir1 or not args.scan_dir2 or not args.scan_origin:
            raise ValueError('For 2d scan, --scan-dir1, --scan-dir2, and --scan-origin must be provided')
        
        # Parse direction vectors and origin from string to float arrays
        dir1 = np.array([float(x) for x in args.scan_dir1.split(',')])
        dir2 = np.array([float(x) for x in args.scan_dir2.split(',')])
        p0 = np.array([float(x) for x in args.scan_origin.split(',')])
        
        # Prepare scan parameters for 2D scan
        scan_params = {
            'nscan1': args.nscan_1,
            'nscan2': args.nscan_2,
            'span1': (args.span_min_1, args.span_max_1),
            'span2': (args.span_min_2, args.span_max_2),
            'p0': p0,
            'dir1': dir1,
            'dir2': dir2
        }
        print('Running 2D scan using scanPlot2D_uff...')
    elif args.scan_mode == '1d':
        # Parse direction vector and origin from string to float arrays
        dir1 = np.array([float(x) for x in args.scan_dir1.split(',')])
        p0 = np.array([float(x) for x in args.scan_origin.split(',')])
        
        # Prepare scan parameters for 1D scan
        scan_params = {
            'nscan': args.nscan_1,
            'span': (args.span_min_1, args.span_max_1),
            'dir': dir1,
            'p0': p0
        }
        print('Running 1D scan using scanPlot_uff...')
    else:
        raise ValueError(f'Unsupported scan mode: {args.scan_mode}')

    # Run scans for each requested type
    results = {}
    for scan_type in run_types:
        print(f"\n=====================================================================")
        print(f"RUNNING {scan_type.upper()} {'2D' if args.scan_mode == '2d' else '1D'} SCAN")
        try:
            success = run_scan(
                molecule=args.molecule,
                substrate=args.substrate,
                output_dir=args.output_dir,
                scan_type=scan_type,
                scan_params=scan_params,
                skip_init=False
            )
            results[scan_type] = 'SUCCESS' if success else 'FAILED'
        except Exception as e:
            print(f"Error running {scan_type} scan: {e}")
            results[scan_type] = 'ERROR'

    print("\n=====================================================================\nALL REQUESTED SCANS COMPLETED")
    print("Results:")
    for scan_type, result in results.items():
        print(f"  {scan_type}: {result}")
    print("====================================================================\n")

    # Compare with LAMMPS data if requested
    if args.compare and results and all(r == 'SUCCESS' for r in results.values()):
        print("\nComparing with LAMMPS data...")
        mol_name = os.path.basename(args.molecule)
        
        # Early return if required LAMMPS directory is not provided
        if args.scan_mode == '2d' and not args.lammps_2d_dir:
            print("Warning: No LAMMPS 2D data directory provided. Skipping comparison.")
        elif args.scan_mode == '1d' and not args.lammps_1d_dir:
            print("Warning: No LAMMPS 1D data directory provided. Skipping comparison.")
        else:
            # Proceed with comparison based on scan mode
            if args.scan_mode == '2d':
                # For 2D scans, generate separate plots for each component
                for scan_type in ['morse', 'coulomb', 'total']:
                    # Map scan types to file names (both LAMMPS and FireCore use 'coul')
                    file_name = {'morse': 'morse', 'coulomb': 'coul', 'total': 'total'}[scan_type]
                    
                    # Check file existence
                    lammps_file = os.path.join(args.lammps_2d_dir, f'{file_name}.dat')
                    firecore_file = os.path.join(args.output_dir, f'{mol_name}_{file_name}_2d.dat')
                    
                    if not os.path.exists(firecore_file):
                        print(f"  Skipping {scan_type} comparison: FireCore data file not found")
                        continue
                        
                    if not os.path.exists(lammps_file):
                        print(f"  Skipping {scan_type} comparison: LAMMPS data file not found")
                        continue
                    
                    # Generate comparison plot
                    compare_2d_with_lammps(
                        mol_name, 
                        args.lammps_2d_dir,
                        args.output_dir, 
                        args.output_dir,
                        scan_type
                    )
                    print(f"  Generated comparison plot for {scan_type} potential")
            else:
                # For 1D scans, use existing comparison function
                compare_with_lammps(mol_name, args.lammps_1d_dir, args.output_dir, args.output_dir)
            
            print(f"Comparison plots saved to {args.output_dir}")

    fig_path = 'fig_path'
    data_path = 'data_path'
