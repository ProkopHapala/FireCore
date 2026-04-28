import sys
import numpy as np
import os
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import re

# Energy component mapping between different implementations
# Format: (component_name, homebrew_name, lammps_name, transform_function)
ENERGY_MAPPINGS = [
    # Basic direct mappings
    ('total', 'Etot', 'Etot', lambda x: x),
    ('kinetic', 'Eke', 'Eke', lambda x: x),
    ('pauli', 'Epauli', 'Epauli', lambda x: x),
    ('coulomb', 'Ecoul', 'Ecoul', lambda x: x),
    ('residual', 'Erres', 'Erres', lambda x: x),
    # Composite mappings can be added here if needed
    # Example: ('electronic', 'Eke', 'Eke', lambda x: x * 0.5)
]

def extract_data_from_xyz(xyz_file):
    """
    Extract relevant parameters from XYZ comments.
    Supports LAMMPS (#key val ...) and homebrew key(val) formats.
    """
    import re
    import numpy as np
    keys = ['ang','dist','Etot','Eke','Epauli','Ecoul','Erres','T','ee','ea','aa']
    # compile patterns
    p_paren = {k: re.compile(rf"{k}\((-?\d+\.\d+)\)") for k in keys}
    p_space = {k: re.compile(rf"{k}\s+(-?\d+\.\d+)" ) for k in keys}
    rows = []
    with open(xyz_file) as f:
        for line in f:
            s = line.strip()
            if not (s.startswith('#') or '(' in s):
                continue
            txt = s.lstrip('#')
            row = {}
            for k in keys:
                m = p_paren[k].search(txt)
                if m:
                    row[k] = float(m.group(1))
                else:
                    m2 = p_space[k].search(txt)
                    row[k] = float(m2.group(1)) if m2 else np.nan
            rows.append(row)
    # assemble arrays
    raw = {k: np.array([r[k] for r in rows], dtype=float) for k in keys}
    return raw

def load_data(homebrew_source, lammps_source):
    """Load and map homebrew and LAMMPS datasets."""
    raw_h = extract_data_from_xyz(homebrew_source)
    homebrew_data = {
        'ang': raw_h['ang'], 'dist': raw_h['dist'],
        'Etot': raw_h['Etot'], 'Eke': raw_h['T'],
        'Epauli': raw_h['ee'], 'Ecoul': raw_h['ea'], 'Erres': raw_h['aa']
    }
    raw_l = extract_data_from_xyz(lammps_source)
    lammps_data = {
        'ang': raw_l['ang'], 'dist': raw_l['dist'],
        'Etot': raw_l['Etot'], 'Eke': raw_l['Eke'],
        'Epauli': raw_l['Epauli'], 'Ecoul': raw_l['Ecoul'], 'Erres': raw_l['Erres']
    }
    return homebrew_data, lammps_data

def create_comparison_grid(homebrew_data, lammps_data, param_x, param_y, energy_mapping):
    """Create a 2D grid comparison based on parameters and energy components."""
    # Extract mapping information
    component_name, homebrew_name, lammps_name, transform_fn = energy_mapping
    
    # Get parameters for grid creation
    homebrew_x_values = homebrew_data.get(param_x, None)
    homebrew_y_values = homebrew_data.get(param_y, None)
    lammps_x_values = lammps_data.get(param_x, None)
    lammps_y_values = lammps_data.get(param_y, None)
    
    if homebrew_x_values is None or homebrew_y_values is None:
        raise ValueError(f"Parameters {param_x} or {param_y} not found in homebrew data")
    if lammps_x_values is None or lammps_y_values is None:
        raise ValueError(f"Parameters {param_x} or {param_y} not found in lammps data")
    
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
    
    # Create a unified grid using all unique values from both datasets
    unique_x = np.unique(np.concatenate((homebrew_x_values, lammps_x_values)))
    unique_y = np.unique(np.concatenate((homebrew_y_values, lammps_y_values)))
    homebrew_grid = np.full((len(unique_x), len(unique_y)), np.nan)
    lammps_grid = np.full((len(unique_x), len(unique_y)), np.nan)
    diff_grid = np.full((len(unique_x), len(unique_y)), np.nan)
    
    # Fill homebrew grid
    for i in range(len(homebrew_energy)):
        if i < len(homebrew_x_values) and i < len(homebrew_y_values):
            x = homebrew_x_values[i]
            y = homebrew_y_values[i]
            
            # Find the nearest values instead of exact matches (more robust)
            if len(unique_x) > 0 and len(unique_y) > 0:
                x_idx = np.abs(unique_x - x).argmin()
                y_idx = np.abs(unique_y - y).argmin()
                homebrew_grid[x_idx, y_idx] = homebrew_energy[i]
    
    # Fill lammps grid
    for i in range(len(lammps_energy)):
        if i < len(lammps_x_values) and i < len(lammps_y_values):
            x = lammps_x_values[i]
            y = lammps_y_values[i]
            
            # Find the nearest values instead of exact matches (more robust)
            if len(unique_x) > 0 and len(unique_y) > 0:
                x_idx = np.abs(unique_x - x).argmin()
                y_idx = np.abs(unique_y - y).argmin()
                lammps_grid[x_idx, y_idx] = lammps_energy[i]
    
    # Calculate differences where both values exist
    for x_idx in range(len(unique_x)):
        for y_idx in range(len(unique_y)):
            if not np.isnan(homebrew_grid[x_idx, y_idx]) and not np.isnan(lammps_grid[x_idx, y_idx]):
                diff_grid[x_idx, y_idx] = homebrew_grid[x_idx, y_idx] - lammps_grid[x_idx, y_idx]
    
    return {
        'x_vals': unique_x,
        'y_vals': unique_y,
        'homebrew_vals': homebrew_grid.flatten(),
        'lammps_vals': lammps_grid.flatten(),
        'diff': diff_grid.flatten()
    }

def plot_comparison(grid_data, param_x, param_y, component_name, output_dir=None):
    """Create comparison plots for a single energy component."""
    # Extract unique values for each parameter
    x_vals = np.unique(grid_data['x_vals'])
    y_vals = np.unique(grid_data['y_vals'])
    nx, ny = len(x_vals), len(y_vals)
    
    # Reshape grids to 2D for plotting
    homebrew_grid = np.reshape(grid_data['homebrew_vals'], (nx, ny))
    lammps_grid = np.reshape(grid_data['lammps_vals'], (nx, ny))
    diff_grid = np.reshape(grid_data['diff'], (nx, ny))
    
    # We'll return these for plotting in the combined function
    return {
        'component_name': component_name,
        'homebrew_grid': homebrew_grid,
        'lammps_grid': lammps_grid,
        'diff_grid': diff_grid,
        'x_vals': x_vals,
        'y_vals': y_vals,
        'stats': {
            'component': component_name,
            'mean_diff': np.mean(diff_grid[~np.isnan(diff_grid)]),
            'abs_mean_diff': np.mean(np.abs(diff_grid[~np.isnan(diff_grid)])),
            'min_diff': np.min(diff_grid[~np.isnan(diff_grid)]),
            'max_diff': np.max(diff_grid[~np.isnan(diff_grid)]),
            'std_diff': np.std(diff_grid[~np.isnan(diff_grid)]),
            'rms_diff': np.sqrt(np.mean(diff_grid[~np.isnan(diff_grid)]**2))
        }
    }

def compare_energy_components(homebrew_data, lammps_data, param_x, param_y, output_dir=None):
    # Fallback to simple Etot time series if no valid ang/dist
    if (param_x not in homebrew_data or param_y not in homebrew_data \
        or homebrew_data[param_x].size == 0 or homebrew_data[param_y].size == 0 \
        or np.isnan(homebrew_data[param_x]).all() or np.isnan(homebrew_data[param_y]).all()):
         print("Using fallback Etot vs frame index plot")
         plt.figure()
         plt.plot(homebrew_data['Etot'], '-o', label='homebrew')
         plt.plot(lammps_data['Etot'], '-x', label='lammps')
         plt.xlabel('Frame index')
         plt.ylabel('Etot')
         plt.title('Total Energy Comparison')
         plt.legend()
         plt.show()
         return []
    
    valid_components = []
    for mapping in ENERGY_MAPPINGS:
        component_name = mapping[0]
        try:
            grid_data = create_comparison_grid(homebrew_data, lammps_data, param_x, param_y, mapping)
            result = plot_comparison(grid_data, param_x, param_y, component_name)
            valid_components.append(result)
        except ValueError as e:
            print(f"Warning: Could not compare {component_name} - {str(e)}")
    
    if not valid_components:
        print("No valid energy components to compare.")
        return []
        
    # Now create a single figure with 3 rows and N columns (one for each component)
    n_components = len(valid_components)
    fig, axs = plt.subplots(3, n_components, figsize=(5*n_components, 15))
    
    # If there's only one component, make sure axs is 2D
    if n_components == 1:
        axs = axs.reshape(3, 1)
    
    # Plot all components in the grid
    for i, result in enumerate(valid_components):
        component_name = result['component_name']
        homebrew_grid = result['homebrew_grid']
        lammps_grid = result['lammps_grid']
        diff_grid = result['diff_grid']
        x_vals = result['x_vals']
        y_vals = result['y_vals']
        
        # Row 1: Homebrew values
        im0 = axs[0, i].imshow(homebrew_grid, origin='lower', aspect='auto',
                           extent=[min(x_vals), max(x_vals), min(y_vals), max(y_vals)])
        axs[0, i].set_title(f'Homebrew {component_name}')
        axs[0, i].set_xlabel(param_x)
        if i == 0:  # Only add y label on the leftmost plot
            axs[0, i].set_ylabel(param_y)
        plt.colorbar(im0, ax=axs[0, i])
        
        # Row 2: LAMMPS values
        im1 = axs[1, i].imshow(lammps_grid, origin='lower', aspect='auto',
                           extent=[min(x_vals), max(x_vals), min(y_vals), max(y_vals)])
        axs[1, i].set_title(f'LAMMPS {component_name}')
        axs[1, i].set_xlabel(param_x)
        if i == 0:  # Only add y label on the leftmost plot
            axs[1, i].set_ylabel(param_y)
        plt.colorbar(im1, ax=axs[1, i])
        
        # Row 3: Difference
        im2 = axs[2, i].imshow(diff_grid, origin='lower', aspect='auto', cmap='seismic',
                           extent=[min(x_vals), max(x_vals), min(y_vals), max(y_vals)])
        axs[2, i].set_title(f'Difference {component_name}')
        axs[2, i].set_xlabel(param_x)
        if i == 0:  # Only add y label on the leftmost plot
            axs[2, i].set_ylabel(param_y)
        plt.colorbar(im2, ax=axs[2, i])
    
    plt.tight_layout()
    
    # Save figure if output directory is provided
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
        filename = os.path.join(output_dir, 'all_components_comparison.png')
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Saved comparison plot to {filename}")
        
        # Also save individual component plots
        for result in valid_components:
            component_name = result['component_name']
            fig_comp, axs_comp = plt.subplots(1, 3, figsize=(15, 5))
            
            # Plot homebrew values
            im0 = axs_comp[0].imshow(result['homebrew_grid'], origin='lower', aspect='auto', extent=[min(result['x_vals']), max(result['x_vals']),  min(result['y_vals']), max(result['y_vals'])])
            axs_comp[0].set_title(f'Homebrew {component_name}')
            axs_comp[0].set_xlabel(param_x)
            axs_comp[0].set_ylabel(param_y)
            plt.colorbar(im0, ax=axs_comp[0])
            
            # Plot LAMMPS values
            im1 = axs_comp[1].imshow(result['lammps_grid'], origin='lower', aspect='auto', extent=[min(result['x_vals']), max(result['x_vals']), min(result['y_vals']), max(result['y_vals'])])
            axs_comp[1].set_title(f'LAMMPS {component_name}')
            axs_comp[1].set_xlabel(param_x)
            axs_comp[1].set_ylabel(param_y)
            plt.colorbar(im1, ax=axs_comp[1])
            
            # Plot difference
            im2 = axs_comp[2].imshow(result['diff_grid'], origin='lower', aspect='auto', cmap='seismic', extent=[min(result['x_vals']), max(result['x_vals']), min(result['y_vals']), max(result['y_vals'])])
            axs_comp[2].set_title(f'Difference {component_name}')
            axs_comp[2].set_xlabel(param_x)
            axs_comp[2].set_ylabel(param_y)
            plt.colorbar(im2, ax=axs_comp[2])
            
            plt.tight_layout()
            component_filename = os.path.join(output_dir, f'comparison_{component_name}.png')
            plt.savefig(component_filename, dpi=300, bbox_inches='tight')
            print(f"Saved individual component plot to {component_filename}")
            plt.close(fig_comp)
        
        # Create and save a summary table
        summary_file = os.path.join(output_dir, 'comparison_summary.txt')
        with open(summary_file, 'w') as f:
            f.write("Component\tMean Diff\tAbs Mean Diff\tMin Diff\tMax Diff\tStd Dev\tRMS Diff\n")
            for stats in [result['stats'] for result in valid_components]:
                f.write(f"{stats['component']}\t{stats['mean_diff']:.6f}\t{stats['abs_mean_diff']:.6f}\t"
                        f"{stats['min_diff']:.6f}\t{stats['max_diff']:.6f}\t{stats['std_diff']:.6f}\t"
                        f"{stats['rms_diff']:.6f}\n")
        print(f"Saved comparison summary to {summary_file}")
    
    return [result['stats'] for result in valid_components]

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Compare energy calculations between homebrew EFF and LAMMPS')
    parser.add_argument('--homebrew-source', required=True,  help='Source file for homebrew EFF data (.xyz or .npy)')
    parser.add_argument('--lammps-source', required=True,  help='Source file for LAMMPS EFF data (.xyz or .npy)')
    parser.add_argument('--output-dir', default='comparison_results',    help='Directory for output files')
    parser.add_argument('--param-x', default='ang',   help='Parameter to use for x-axis (default: ang)')
    parser.add_argument('--param-y', default='dist',  help='Parameter to use for y-axis (default: dist)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data and plot total energy comparison
    print(f"Loading data from: {args.homebrew_source} and {args.lammps_source}")
    homebrew_data, lammps_data = load_data(args.homebrew_source, args.lammps_source)
    hb = homebrew_data.get('Etot', np.array([]))
    lp = lammps_data.get('Etot', np.array([]))
    if hb.size == 0 or lp.size == 0:
        print("Error: Etot data not found in one of the sources.")
        sys.exit(1)
    plt.figure()
    plt.plot(hb, '-o', label='homebrew')
    plt.plot(lp, '-x', label='lammps')
    plt.xlabel('Frame index')
    plt.ylabel('Etot')
    plt.title('Total Energy Comparison')
    plt.legend()
    output_png = os.path.join(args.output_dir, 'etot_comparison.png')
    plt.savefig(output_png)
    print(f"Saved total energy comparison plot to {output_png}")
    plt.show()
    sys.exit(0)
