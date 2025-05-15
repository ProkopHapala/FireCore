import sys
import numpy as np
import os
import shutil
import tempfile
import matplotlib.pyplot as plt
from pathlib import Path

sys.path.append("../../")
from pyBall import lammps_utils as lu
from pyBall.atomicUtils import load_xyz_movie

# EFF LAMMPS script template generator with relaxation options
def generate_lammps_script(data_file, output_file, relax_nuclei=False, relax_electrons=False, max_steps=4000):
    """Generate LAMMPS script with appropriate relaxation settings"""
    script = """
units           electron
newton          on
boundary        f f f
atom_style      electron

read_data       {data_file}

pair_style      eff/cut 100.0 ecp 1 O
pair_coeff      * *

compute         energies all pair eff/cut
variable        eke equal c_energies[1]
variable        epauli equal c_energies[2]
variable        ecoul equal c_energies[3]
variable        erres equal c_energies[4]

thermo          1
thermo_style    custom step etotal v_eke v_epauli v_ecoul v_erres
thermo_modify   format float %23.15g
"""
    
    # Create groups for nuclei and electrons
    script += """
group           nuclei type 1 2
group           electrons type 3
"""
    
    # Add fixes based on relaxation flags
    if not relax_nuclei:
        script += "fix             freeze_nuclei nuclei setforce 0.0 0.0 0.0\n"
    if not relax_electrons:
        script += "fix             freeze_electrons electrons setforce 0.0 0.0 0.0\n"
    
    # Skip minimization if both nuclei and electrons are fixed
    if relax_nuclei or relax_electrons:
        script += """
min_style       cg
compute         1 all property/atom spin eradius erforce
dump            2 all custom 1 mini.lammpstrj id type q c_1[1] c_1[2] x y z fx fy fz c_1[3]
minimize        0 1e-6 2000 {max_steps}
"""
    else:
        script += """
# No minimization needed - both nuclei and electrons are fixed
run             0
"""
    
    # Add final write_data
    script += "write_data      {output_file}\n"

    text = script.format(data_file=data_file, output_file=output_file, max_steps=max_steps)
    print(text)
    return text

# EFF-specific field mapping for energy extraction
EFF_FIELD_MAP = {
    'total': {'pattern': 'TotEng', 'index': 1},
    'eke': {'pattern': 'v_eke', 'index': 2},
    'epauli': {'pattern': 'v_epauli', 'index': 3},
    'ecoul': {'pattern': 'v_ecoul', 'index': 4},
    'erres': {'pattern': 'v_erres', 'index': 5}
}

def create_eff_input(elements, positions, charges, radii, **kwargs):
    """Callback function to create eFF input from XYZ data"""
    return lu.create_lammps_data_eff(elements, positions, charges, radii, **kwargs)

def process_eff_output(output_file, stdout):
    """Callback function to process eFF output"""
    # Read final geometry
    final_elements, final_positions = lu.read_lammps_data(output_file)
    
    # Return processed data
    return {
        'elements': final_elements,
        'positions': final_positions
    }

def plot_energy_landscape(Xs, Ys, Es, energy_component="total", Espan=None, save_path=None, reshape=None):
    """Plot energy landscape from arrays using imshow (simple and robust)."""
    # Ensure inputs are numpy arrays
    Xs = np.asarray(Xs)
    Ys = np.asarray(Ys)
    Es = np.asarray(Es)
    # Filter out NaN entries
    mask = ~(np.isnan(Xs) | np.isnan(Ys) | np.isnan(Es))
    Xs, Ys, Es = Xs[mask], Ys[mask], Es[mask]
    # Debug: print ranges
    if Xs.size == 0:
        print(f"No valid data to plot {energy_component}")
        return None
    print(f"DEBUG {energy_component}: Xs min {Xs.min():.3f}, max {Xs.max():.3f}, Ys min {Ys.min():.3f}, max {Ys.max():.3f}, Es min {Es.min():.3f}, max {Es.max():.3f}")
    # Determine grid axes and heatmap grid
    plt.figure(figsize=(10, 8))
    if reshape and Es.size == reshape[0] * reshape[1]:
        # reshape array directly
        grid = Es.reshape(reshape)
        xs = np.unique(Xs)
        ys = np.unique(Ys)
    else:
        xs = np.unique(Xs)
        ys = np.unique(Ys)
        grid = np.full((len(xs), len(ys)), np.nan)
        xi = {v: i for i, v in enumerate(xs)}
        yi = {v: i for i, v in enumerate(ys)}
        for x, y, e in zip(Xs, Ys, Es): grid[xi[x], yi[y]] = e
    # Check for degenerate axes and warn
    if len(xs) == 1:
        print(f"WARNING: Only one unique X value: {xs[0]:.3f}")
    if len(ys) == 1:
        print(f"WARNING: Only one unique Y value: {ys[0]:.3f}")
    # Determine extents, expanding degenerate axes
    x_min, x_max = xs.min(), xs.max()
    if x_min == x_max:
        delta = x_min * 0.05 if x_min != 0 else 1.0
        x_min -= delta; x_max += delta
        print(f"WARNING: X axis degenerate; expanded by {delta:.3f}")
    y_min, y_max = ys.min(), ys.max()
    if y_min == y_max:
        delta = y_min * 0.05 if y_min != 0 else 1.0
        y_min -= delta; y_max += delta
        print(f"WARNING: Y axis degenerate; expanded by {delta:.3f}")
    # Debug: inspect grid and extents
    valid = np.count_nonzero(~np.isnan(grid)); total = grid.size
    print(f"DEBUG {energy_component}: grid shape {grid.shape}, valid entries {valid}/{total}")
    print(f"DEBUG {energy_component}: extent x [{x_min:.3f}, {x_max:.3f}], y [{y_min:.3f}, {y_max:.3f}]")
    # Plot heatmap
    vmin = np.nanmin(grid)
    vmax = vmin + Espan if Espan is not None else None
    im = plt.imshow(grid, extent=[x_min, x_max, y_min, y_max], aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax)
    plt.colorbar(im, label=f'{energy_component} Energy')
    plt.xlabel('Distance (Å)'); plt.ylabel('Angle (rad)')
    plt.title(f'LAMMPS eFF {energy_component} Energy')
    
    if save_path:
        plt.savefig(save_path)
    # No need to return grid; return None to avoid unbound variable
    return None

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate 2D energy map with eFF LAMMPS')
    parser.add_argument('xyz_file', help='XYZ file containing frame data with parameters in comments')
    parser.add_argument('--output-dir', default='lammps_results', help='Directory for output files')
    parser.add_argument('--scratch-dir', default=None, help='Temporary directory for intermediate files')
    parser.add_argument('--box-size', type=float, default=500, help='Box size for simulation')
    parser.add_argument('--max-frames', type=int, default=None, help='Maximum number of frames to process')
    parser.add_argument('--lammps-exec', default='lmp_serial', help='LAMMPS executable name/path')
    parser.add_argument('--electron-type', type=int, default=None, help='Atom type to use for electrons')
    parser.add_argument('--relax-nuclei', action='store_true', help='Allow nuclei (atoms) to relax during minimization')
    parser.add_argument('--relax-electrons', action='store_true', help='Allow electrons to relax during minimization')
    parser.add_argument('--max-steps', type=int, default=4000, help='Maximum number of minimization steps')
    parser.add_argument('--keep-scratch', action='store_true', help='Keep scratch directory after completion')
    parser.add_argument('--energy-span', type=float, default=5.0, help='Energy span for plot in eV')
    parser.add_argument('--scan-file', default=None, help='Original scan XYZ with ang/dist comments to infer grid shape')
    
    args = parser.parse_args()
    
    # Auto-detect original scan file from OUTeFF movie
    if not args.scan_file and args.xyz_file.endswith('_OUTeFF.xyz'):
        candidate = args.xyz_file.replace('_OUTeFF', '')
        # check same dir
        if os.path.exists(candidate):
            args.scan_file = candidate
            print(f"INFO: Inferred scan file {candidate}")
        else:
            # check export/scan_data subdir
            base = os.path.basename(candidate)
            alt = os.path.join('export', 'scan_data', base)
            if os.path.exists(alt):
                args.scan_file = alt
                print(f"INFO: Inferred scan file {alt}")
    
    # Load frames and parse parameters
    print(f"Loading XYZ movie from: {args.xyz_file}")
    xyz_movie = load_xyz_movie(args.xyz_file)
    print(f"Loaded {len(xyz_movie)} frames")
    
    # Compute geometric parameters: H-O-H angle and O-H bond length
    angles, dists = [], []
    for idx, (es, apos, *_) in enumerate(xyz_movie):
        if '8' not in es:
            raise ValueError(f"Frame {idx} missing oxygen atom")
        iO = es.index('8')
        iH = [i for i, e in enumerate(es) if e == '1']
        if len(iH) < 2:
            raise ValueError(f"Frame {idx} has {len(iH)} hydrogens; expected at least 2")
        pO = apos[iO]
        pH1, pH2 = apos[iH[0]], apos[iH[1]]
        # O-H bond length
        d = np.linalg.norm(pH1 - pO)
        # H-O-H angle
        v1, v2 = pH1 - pO, pH2 - pO
        cosang = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cosang = np.clip(cosang, -1.0, 1.0)
        a = np.arccos(cosang)
        angles.append(a); dists.append(d)
    
    params = {'ang': np.array(angles), 'dist': np.array(dists)}
    nrec = len(angles)
    # Trim parameters to max_frames if set
    if args.max_frames is not None:
        for key in params:
            params[key] = params[key][:args.max_frames]
        nrec = args.max_frames
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Set up scratch directory
    if args.scratch_dir:
        scratch_dir = args.scratch_dir
        os.makedirs(scratch_dir, exist_ok=True)
    else:
        scratch_dir = tempfile.mkdtemp(prefix='lammps_scratch_')
    
    print(f"Using scratch directory: {scratch_dir}")
    
    try:
        # Process frames with LAMMPS using the modular function
        script_template = generate_lammps_script(
            "{data_file}", "{output_file}",
            relax_nuclei=args.relax_nuclei,
            relax_electrons=args.relax_electrons,
            max_steps=args.max_steps
        )
        results = lu.process_xyz_with_lammps(
            xyz_movie,
            create_eff_input,
            process_eff_output,
            script_template,
            field_map=EFF_FIELD_MAP,
            work_dir=scratch_dir,
            box_size=args.box_size,
            max_frames=args.max_frames,
            lammps_exec=args.lammps_exec,
            electron_type=args.electron_type
        )
        
        # Extract and organize energy data
        frames = []
        total_energies = []
        kinetic_energies = []
        pauli_energies = []
        coulomb_energies = []
        residual_energies = []
        
        for result in results:
            frames.append(result['frame'])
            energies = result['energies']
            if energies:
                total_energies.append(energies.get('total', 0))
                kinetic_energies.append(energies.get('eke', 0))
                pauli_energies.append(energies.get('epauli', 0))
                coulomb_energies.append(energies.get('ecoul', 0))
                residual_energies.append(energies.get('erres', 0))
            else:
                total_energies.append(np.nan)
                kinetic_energies.append(np.nan)
                pauli_energies.append(np.nan)
                coulomb_energies.append(np.nan)
                residual_energies.append(np.nan)
        
        # Output a summary
        print(f"\nSuccessfully processed {len(results)} frames")
        
        # Write summary of energies to file
        summary_file = os.path.join(args.output_dir, "lammps_energy_summary.txt")
        with open(summary_file, 'w') as f:
            # Include relaxation info in header
            relax_info = []
            if args.relax_nuclei: relax_info.append("nuclei")
            if args.relax_electrons: relax_info.append("electrons")
            relax_str = f" (relaxed: {', '.join(relax_info)})" if relax_info else " (fully fixed)"
            f.write(f"# LAMMPS eFF simulation{relax_str}\n")
            f.write("Frame Ang Dist Total_Energy E_ke E_pauli E_coul E_rres\n")
            for i, result in enumerate(results):
                energies = result['energies']
                if energies and 'ang' in params and 'dist' in params:
                    f.write(f"{result['frame']} {params['ang'][i]:.6f} {params['dist'][i]:.6f} "
                            f"{energies.get('total', 0):.6f} {energies.get('eke', 0):.6f} "
                            f"{energies.get('epauli', 0):.6f} {energies.get('ecoul', 0):.6f} "
                            f"{energies.get('erres', 0):.6f}\n")
        
        print(f"Energy summary written to {summary_file}")
        
        # Save data as NumPy array for later comparison
        energy_data = np.column_stack((params['ang'], params['dist'], total_energies,
                                      kinetic_energies, pauli_energies, coulomb_energies,
                                      residual_energies))
        np.save(os.path.join(args.output_dir, "lammps_energy_data.npy"), energy_data)
        
        # Create plots for each energy component
        energy_components = {
            'total': total_energies,
            'kinetic': kinetic_energies,
            'pauli': pauli_energies,
            'coulomb': coulomb_energies,
            'residual': residual_energies
        }
        
        if args.scan_file:
            # Parse original scan file to infer grid dims matching processed frames
            nres = len(results)
            scan_ang = []
            scan_dist = []
            with open(args.scan_file) as f:
                for line in f:
                    if line.startswith('#') and 'ang' in line and 'dist' in line:
                        parts = line.lstrip('#').split()
                        ia = parts.index('ang') + 1; id = parts.index('dist') + 1
                        scan_ang.append(float(parts[ia])); scan_dist.append(float(parts[id]))
                        if len(scan_ang) >= nres:
                            break
            ua = np.unique(scan_ang)
            ud = np.unique(scan_dist)
            dims = (len(ua), len(ud))
            # Reshape and plot for each component
            for comp, vals in energy_components.items():
                grid = np.array(vals[:len(scan_ang)]).reshape(dims)
                plt.figure(figsize=(10, 8))
                vmin = np.nanmin(grid)
                vmax = vmin + args.energy_span
                im = plt.imshow(grid, extent=[ud.min(), ud.max(), ua.min(), ua.max()], aspect='auto', cmap='inferno', vmin=vmin, vmax=vmax)
                plt.colorbar(im, label=f'{comp} Energy')
                plt.xlabel('Distance (Å)'); plt.ylabel('Angle (rad)')
                plt.title(f'LAMMPS eFF {comp} Energy')
                plot_file = os.path.join(args.output_dir, f"lammps_map2d_{comp}.png")
                plt.savefig(plot_file)
                print(f"Plot for {comp} energy saved to {plot_file}")
        else:
            for component_name, energy_values in energy_components.items():
                plot_file = os.path.join(args.output_dir, f"lammps_map2d_{component_name}.png")
                plot_energy_landscape(params['ang'], params['dist'], energy_values,
                                     component_name, Espan=args.energy_span,
                                     save_path=plot_file)
                print(f"Plot for {component_name} energy saved to {plot_file}")
        
        # Generate XYZ file with energy data
        xyz_file = os.path.join(args.output_dir, "lammps_final_geometries.xyz")
        with open(xyz_file, 'w') as f:
            for i, result in enumerate(results):
                if i >= len(params['ang']) or i >= len(params['dist']):
                    continue  # Skip if we don't have parameters for this frame
                    
                processed = result.get('processed_data', {})
                elements = processed.get('elements', [])
                positions = processed.get('positions', [])
                energy = result['energies'] if 'energies' in result else None
                
                if elements and positions and len(elements) == len(positions):
                    f.write(f"{len(elements)}\n")
                    # Add parameter data
                    if energy is not None:
                        # Include relaxation info in XYZ comments
                        relax_info = []
                        if args.relax_nuclei: relax_info.append("nuclei")
                        if args.relax_electrons: relax_info.append("electrons")
                        relax_str = f"relax {','.join(relax_info)} " if relax_info else "fixed "
                        
                        # Format the comment line with parameters and energy values
                        comment = (f"#ang {params['ang'][i]:.6f} dist {params['dist'][i]:.6f} {relax_str}"
                                  f"Etot {energy.get('total', 0):.6f} "
                                  f"Eke {energy.get('eke', 0):.6f} "
                                  f"Epauli {energy.get('epauli', 0):.6f} "
                                  f"Ecoul {energy.get('ecoul', 0):.6f} "
                                  f"Erres {energy.get('erres', 0):.6f}")
                        f.write(f"{comment}\n")
                    else:
                        f.write("\n")
                    
                    for elem, pos in zip(elements, positions):
                        f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
        
        print(f"Final geometries with parameters written to {xyz_file}")
        
        # Show plots at the end
        plt.show()
        
    finally:
        # Clean up scratch directory unless requested to keep it
        if not args.keep_scratch and (args.scratch_dir is None):
            shutil.rmtree(scratch_dir)
            print(f"Removed scratch directory: {scratch_dir}")
        elif args.keep_scratch:
            print(f"Kept scratch directory: {scratch_dir}")
