import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

sys.path.append("../../")

# Energy component mapping between different implementations
# Format: (component_name, homebrew_name, lammps_name, transform_function)
ENERGY_MAPPINGS = [
    # Basic direct mappings
    ('total', 'Etot', 'total', lambda x: x),
    ('kinetic', 'Eke', 'eke', lambda x: x),
    ('pauli', 'Epauli', 'epauli', lambda x: x),
    ('coulomb', 'Ecoul', 'ecoul', lambda x: x),
    ('residual', 'Erres', 'erres', lambda x: x),
    # Composite mappings can be added here if needed
    # Example: ('electronic', 'Eke', 'eke', lambda x: x * 0.5)
]

def extract_data_from_xyz(xyz_file):
    """Extract energy and parameter data from XYZ file with comments.
    Returns a dictionary with extracted parameters and energy components.
    """
    params = {}
    frames = []
    
    with open(xyz_file, 'r') as f:
        line_num = 0
        frame_num = 0
        
        while True:
            # Read number of atoms
            line = f.readline()
            if not line:
                break
            
            num_atoms = int(line.strip())
            
            # Read comment line with parameters
            comment = f.readline().strip()
            
            # Extract energy and parameter data
            if comment.startswith('#'):
                # Format: #param1 value1 param2 value2 ...
                parts = comment[1:].strip().split()
                for i in range(0, len(parts) - 1, 2):
                    key = parts[i]
                    val = float(parts[i + 1])
                    
                    if key not in params:
                        params[key] = []
                    
                    # Extend list if needed
                    while len(params[key]) <= frame_num:
                        params[key].append(np.nan)
                    
                    params[key][frame_num] = val
            elif 'Etot' in comment:
                # Format like: Etot(value) Eke(value) ...
                for part in comment.split():
                    if '(' in part and ')' in part:
                        key = part.split('(')[0]
                        val = float(part.split('(')[1].split(')')[0])
                        
                        if key not in params:
                            params[key] = []
                        
                        # Extend list if needed
                        while len(params[key]) <= frame_num:
                            params[key].append(np.nan)
                        
                        params[key][frame_num] = val
            
            # Skip atom positions
            for _ in range(num_atoms):
                f.readline()
            
            frames.append(frame_num)
            frame_num += 1
    
    # Convert lists to numpy arrays
    for key in params:
        params[key] = np.array(params[key])
    
    params['frames'] = np.array(frames)
    return params

def extract_data_from_npy(npy_file):
    """Extract data from numpy file.
    Expected format: columns are [ang, dist, total, eke, epauli, ecoul, erres]
    """
    data = np.load(npy_file)
    params = {
        'ang': data[:, 0],
        'dist': data[:, 1],
        'total': data[:, 2],
        'eke': data[:, 3],
        'epauli': data[:, 4],
        'ecoul': data[:, 5],
        'erres': data[:, 6],
        'frames': np.arange(len(data))
    }
    return params

def load_data(homebrew_source, lammps_source):
    """Load data from either XYZ files or NPY files."""
    # Try to load homebrew data
    if homebrew_source.endswith('.xyz'):
        homebrew_data = extract_data_from_xyz(homebrew_source)
    elif homebrew_source.endswith('.npy'):
        homebrew_data = extract_data_from_npy(homebrew_source)
    else:
        raise ValueError(f"Unsupported homebrew data source: {homebrew_source}")
    
    # Try to load LAMMPS data
    if lammps_source.endswith('.xyz'):
        lammps_data = extract_data_from_xyz(lammps_source)
    elif lammps_source.endswith('.npy'):
        lammps_data = extract_data_from_npy(lammps_source)
    else:
        raise ValueError(f"Unsupported LAMMPS data source: {lammps_source}")
    
    return homebrew_data, lammps_data

def create_comparison_grid(homebrew_data, lammps_data, param_x, param_y, energy_mapping):
    """Create a 2D grid comparison based on parameters and energy components."""
    # Extract mapping information
    component_name, homebrew_name, lammps_name, transform_fn = energy_mapping
    
    # Get parameters for grid creation
    x_values = homebrew_data.get(param_x, None)
    y_values = homebrew_data.get(param_y, None)
    
    if x_values is None or y_values is None:
        raise ValueError(f"Parameters {param_x} or {param_y} not found in data")
    
    # Get energy values
    homebrew_energy = homebrew_data.get(homebrew_name, None)
    lammps_energy = lammps_data.get(lammps_name, None)
    
    if homebrew_energy is None:
        raise ValueError(f"Energy component {homebrew_name} not found in homebrew data")
    if lammps_energy is None:
        raise ValueError(f"Energy component {lammps_name} not found in LAMMPS data")
    
    # Apply transformation if any
    homebrew_energy = transform_fn(homebrew_energy)
    lammps_energy = transform_fn(lammps_energy)
    
    # Create grid matrices
    unique_x = np.unique(x_values)
    unique_y = np.unique(y_values)
    homebrew_grid = np.full((len(unique_x), len(unique_y)), np.nan)
    lammps_grid = np.full((len(unique_x), len(unique_y)), np.nan)
    diff_grid = np.full((len(unique_x), len(unique_y)), np.nan)
    
    # Fill grids
    for i, x in enumerate(x_values):
        x_idx = np.where(unique_x == x)[0][0]
        y_idx = np.where(unique_y == y_values[i])[0][0]
        
        homebrew_grid[x_idx, y_idx] = homebrew_energy[i]
        lammps_grid[x_idx, y_idx] = lammps_energy[i]
        diff_grid[x_idx, y_idx] = homebrew_energy[i] - lammps_energy[i]
    
    return {
        'x': unique_x,
        'y': unique_y,
        'homebrew': homebrew_grid,
        'lammps': lammps_grid,
        'diff': diff_grid
    }

def plot_comparison(grid_data, param_x, param_y, component_name, output_dir=None):
    """Create comparison plots for a single energy component."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Determine appropriate color scale limits
    homebrew_min = np.nanmin(grid_data['homebrew'])
    homebrew_max = np.nanmax(grid_data['homebrew'])
    lammps_min = np.nanmin(grid_data['lammps'])
    lammps_max = np.nanmax(grid_data['lammps'])
    
    # Use the same scale for homebrew and LAMMPS
    energy_min = min(homebrew_min, lammps_min)
    energy_max = max(homebrew_max, lammps_max)
    
    # Compute absolute max for difference
    diff_abs_max = np.nanmax(np.abs(grid_data['diff']))
    
    # Plot homebrew data
    im0 = axes[0].imshow(grid_data['homebrew'], 
                         extent=[grid_data['y'].min(), grid_data['y'].max(), 
                                grid_data['x'].max(), grid_data['x'].min()],
                         aspect='auto', cmap='inferno', vmin=energy_min, vmax=energy_max)
    axes[0].set_xlabel(f'{param_y} (Å)')
    axes[0].set_ylabel(f'{param_x} (rad)')
    axes[0].set_title(f'Homebrew {component_name} Energy')
    plt.colorbar(im0, ax=axes[0], label='Energy (eV)')
    
    # Plot LAMMPS data
    im1 = axes[1].imshow(grid_data['lammps'], 
                         extent=[grid_data['y'].min(), grid_data['y'].max(), 
                                grid_data['x'].max(), grid_data['x'].min()],
                         aspect='auto', cmap='inferno', vmin=energy_min, vmax=energy_max)
    axes[1].set_xlabel(f'{param_y} (Å)')
    axes[1].set_ylabel(f'{param_x} (rad)')
    axes[1].set_title(f'LAMMPS {component_name} Energy')
    plt.colorbar(im1, ax=axes[1], label='Energy (eV)')
    
    # Plot difference
    im2 = axes[2].imshow(grid_data['diff'], 
                         extent=[grid_data['y'].min(), grid_data['y'].max(), 
                                grid_data['x'].max(), grid_data['x'].min()],
                         aspect='auto', cmap='coolwarm', 
                         vmin=-diff_abs_max, vmax=diff_abs_max)
    axes[2].set_xlabel(f'{param_y} (Å)')
    axes[2].set_ylabel(f'{param_x} (rad)')
    axes[2].set_title(f'Difference (Homebrew - LAMMPS)')
    plt.colorbar(im2, ax=axes[2], label='Energy Difference (eV)')
    
    plt.tight_layout()
    
    # Save figure if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, f'comparison_{component_name}.png')
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved comparison plot to {filename}")
    
    # Calculate and return statistics
    flat_diff = grid_data['diff'].flatten()
    flat_diff = flat_diff[~np.isnan(flat_diff)]
    
    stats = {
        'component': component_name,
        'mean_diff': np.mean(flat_diff),
        'abs_mean_diff': np.mean(np.abs(flat_diff)),
        'min_diff': np.min(flat_diff),
        'max_diff': np.max(flat_diff),
        'std_diff': np.std(flat_diff),
        'rms_diff': np.sqrt(np.mean(flat_diff**2))
    }
    
    return stats

def compare_energy_components(homebrew_data, lammps_data, param_x, param_y, output_dir=None):
    """Compare all energy components based on the defined mappings."""
    stats_summary = []
    
    for mapping in ENERGY_MAPPINGS:
        component_name = mapping[0]
        try:
            grid_data = create_comparison_grid(homebrew_data, lammps_data, param_x, param_y, mapping)
            stats = plot_comparison(grid_data, param_x, param_y, component_name, output_dir)
            stats_summary.append(stats)
        except ValueError as e:
            print(f"Warning: Could not compare {component_name} - {str(e)}")
    
    # Create and save a summary table
    if output_dir:
        summary_file = os.path.join(output_dir, 'comparison_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("Component\tMean Diff\tAbs Mean Diff\tMin Diff\tMax Diff\tStd Dev\tRMS Diff\n")
            for stats in stats_summary:
                f.write(f"{stats['component']}\t{stats['mean_diff']:.6f}\t{stats['abs_mean_diff']:.6f}\t"
                        f"{stats['min_diff']:.6f}\t{stats['max_diff']:.6f}\t{stats['std_diff']:.6f}\t"
                        f"{stats['rms_diff']:.6f}\n")
        print(f"Saved comparison summary to {summary_file}")
    
    return stats_summary

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare energy calculations between homebrew EFF and LAMMPS')
    parser.add_argument('--homebrew-source', required=True, 
                        help='Source file for homebrew EFF data (.xyz or .npy)')
    parser.add_argument('--lammps-source', required=True, 
                        help='Source file for LAMMPS EFF data (.xyz or .npy)')
    parser.add_argument('--output-dir', default='comparison_results', 
                        help='Directory for output files')
    parser.add_argument('--param-x', default='ang', 
                        help='Parameter to use for x-axis (default: ang)')
    parser.add_argument('--param-y', default='dist', 
                        help='Parameter to use for y-axis (default: dist)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data from sources
    print(f"Loading homebrew data from: {args.homebrew_source}")
    print(f"Loading LAMMPS data from: {args.lammps_source}")
    homebrew_data, lammps_data = load_data(args.homebrew_source, args.lammps_source)
    
    # Compare energy components
    print("\nComparing energy components...")
    stats_summary = compare_energy_components(
        homebrew_data, lammps_data, 
        args.param_x, args.param_y,
        args.output_dir
    )
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print("-" * 80)
    print(f"{'Component':<10} | {'Mean Diff':>10} | {'Abs Mean':>10} | {'Min Diff':>10} | {'Max Diff':>10} | {'Std Dev':>10} | {'RMS Diff':>10}")
    print("-" * 80)
    for stats in stats_summary:
        print(f"{stats['component']:<10} | {stats['mean_diff']:>10.6f} | {stats['abs_mean_diff']:>10.6f} | "
              f"{stats['min_diff']:>10.6f} | {stats['max_diff']:>10.6f} | {stats['std_diff']:>10.6f} | "
              f"{stats['rms_diff']:>10.6f}")
    
    plt.show()
