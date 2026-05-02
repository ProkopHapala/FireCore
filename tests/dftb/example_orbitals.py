#!/usr/bin/env python
"""
Example: Plot molecular orbitals from DFTB+ using waveplot

This demonstrates how to:
1. Run DFTB+ with eigenvector output
2. Use waveplot to generate cube files
3. Read and visualize molecular orbitals
"""

import sys
sys.path.insert(0, '/home/prokophapala/git/FireCore')

import argparse
import subprocess
import os
from pathlib import Path
from ase import Atoms
from ase.io import read

from pyBall import dftb_utils as dftbu
from pyBall.AtomicSystem import AtomicSystem

def plot_orbitals_2d(cube_files, atoms, args):
    """Plot orbitals as 2D projections using matplotlib imshow."""
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    # Parse orbital selection
    if args.orbital == 'all':
        selected_cubes = cube_files
    elif args.orbital == 'occupied':
        selected_cubes = [c for c in cube_files if '1-1-1' in str(c) or '1-1-2' in str(c) or '1-1-3' in str(c)]
    elif args.orbital == 'virtual':
        selected_cubes = [c for c in cube_files if '1-1-4' in str(c) or '1-1-5' in str(c) or '1-1-6' in str(c)]
    else:
        indices = [int(i) for i in args.orbital.split(',')]
        selected_cubes = [cube_files[i-1] for i in indices if i <= len(cube_files)]
    
    # Also add density and potential if requested
    if args.plot_density or args.electrostatic:
        density_cubes = [c for c in cube_files if 'abs2' in str(c) and 'diff' not in str(c)]
        selected_cubes.extend(density_cubes)
    
    if args.electrostatic:
        esp_cubes = [c for c in cube_files if 'esp' in str(c).lower()]
        selected_cubes.extend(esp_cubes)
            
    # Plot each orbital
    for cube_file in selected_cubes:
        try:
            # Read cube with grid info
            data, atoms_cube, origin, spacing = dftbu.read_cube_with_grid(cube_file)
            
            # Get grid dimensions
            nx, ny, nz = data.shape
            
            # Extract 2D slice based on plane
            if args.plane == 'xy':
                # Take middle z slice
                slice_data = data[:, :, nz // 2]
                xlabel, ylabel = 'x (Å)', 'y (Å)'
                extent = [origin[0], origin[0] + nx*spacing[0], 
                         origin[1], origin[1] + ny*spacing[1]]
                # Project atoms to 2D - use original atoms, not cube atoms
                atom_pos_2d = atoms.positions[:, :2]
            elif args.plane == 'xz':
                # Take middle y slice
                slice_data = data[:, ny // 2, :]
                xlabel, ylabel = 'x (Å)', 'z (Å)'
                extent = [origin[0], origin[0] + nx*spacing[0], 
                         origin[2], origin[2] + nz*spacing[2]]
                atom_pos_2d = atoms.positions[:, [0, 2]]
            else:  # yz
                # Take middle x slice
                slice_data = data[nx // 2, :, :]
                xlabel, ylabel = 'y (Å)', 'z (Å)'
                extent = [origin[1], origin[1] + ny*spacing[1], 
                         origin[2], origin[2] + nz*spacing[2]]
                atom_pos_2d = atoms.positions[:, [1, 2]]
            
            # Plot
            plt.figure(figsize=(8, 6))
            im = plt.imshow(slice_data.T, origin='lower', cmap='RdBu_r', aspect='auto', extent=extent)
            plt.colorbar(im, label='Orbital value')
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(f'{Path(cube_file).name} - {args.plane.upper()} projection')
            
            # Overlay atoms
            if args.plane == 'xy':
                atom_pos_2d = atoms.positions[:, :2]
            elif args.plane == 'xz':
                atom_pos_2d = atoms.positions[:, [0, 2]]
            else:
                atom_pos_2d = atoms.positions[:, [1, 2]]
            
            for i, (x, y) in enumerate(atom_pos_2d):
                symbol = atoms.get_chemical_symbols()[i]
                plt.scatter(x, y, c='red', s=20, marker='+', zorder=10)
                plt.text(x, y, f'{symbol}{i}', ha='center', va='center', 
                         color='white', fontsize=8, fontweight='bold',
                         bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5), zorder=11)
            
            # Save
            output_file = output_dir / f"{Path(cube_file).stem}_{args.plane}.png"
            plt.savefig(output_file, dpi=150, bbox_inches='tight')
            plt.close()
            
            print(f"  Saved: {output_file}")
            
        except Exception as e:
            print(f"  Error plotting {cube_file}: {e}")

def plot_density_at_points(cube_files, atoms, args):
    """Plot electron density on 2D grid with atoms overlay (like orbital projections)."""
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    
    # Get density cube file
    density_cube = None
    for cube in cube_files:
        cube_str = str(cube)
        if 'abs2' in cube_str and 'diff' not in cube_str:
            density_cube = cube
            break
    
    if not density_cube:
        print("  No density cube file found")
        return
    
    # Read density data with grid info
    data, atoms_cube, origin, spacing = dftbu.read_cube_with_grid(density_cube)
    
    # Get grid dimensions
    nx, ny, nz = data.shape
    
    # Extract 2D slice based on plane (same as orbital plotting)
    if args.plane == 'xy':
        # Take middle z slice
        slice_data = data[:, :, nz // 2]
        xlabel, ylabel = 'x (Å)', 'y (Å)'
        extent = [origin[0], origin[0] + nx*spacing[0], 
                 origin[1], origin[1] + ny*spacing[1]]
        # Project atoms to 2D
        atom_positions_2d = atoms.positions[:, :2]
    elif args.plane == 'xz':
        # Take middle y slice
        slice_data = data[:, ny // 2, :]
        xlabel, ylabel = 'x (Å)', 'z (Å)'
        extent = [origin[0], origin[0] + nx*spacing[0], 
                 origin[2], origin[2] + nz*spacing[2]]
        atom_positions_2d = atoms.positions[:, [0, 2]]
    else:  # yz
        # Take middle x slice
        slice_data = data[nx // 2, :, :]
        xlabel, ylabel = 'y (Å)', 'z (Å)'
        extent = [origin[1], origin[1] + ny*spacing[1], 
                 origin[2], origin[2] + nz*spacing[2]]
        atom_positions_2d = atoms.positions[:, [1, 2]]
    
    # Plot with imshow
    output_dir = Path(args.output_dir)
    output_dir.mkdir(exist_ok=True, parents=True)
    
    plt.figure(figsize=(8, 6))
    im = plt.imshow(slice_data.T, origin='lower', cmap='viridis', aspect='auto', extent=extent)
    plt.colorbar(im, label='Electron density')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(f'Electron density - {args.plane.upper()} projection')
    
    # Overlay atoms - small dot with symbol and number
    for i, (x, y) in enumerate(atom_positions_2d):
        symbol = atoms.get_chemical_symbols()[i]
        plt.scatter(x, y, c='red', s=20, marker='+', zorder=10)
        plt.text(x, y, f'{symbol}{i}', ha='center', va='center', 
                 color='white', fontsize=8, fontweight='bold',
                 bbox=dict(boxstyle='round,pad=0.2', facecolor='black', alpha=0.5), zorder=11)
    
    output_file = output_dir / f"density_{args.plane}.png"
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_file}")

def load_molecule_xyz(filename, use_ase=True):
    """Load molecule from XYZ file using ASE or AtomicSystem."""
    if use_ase:
        atoms = read(filename)
        return atoms
    else:
        mol = AtomicSystem(fname=filename)
        atoms = Atoms(
            symbols=mol.enames,
            positions=mol.apos,
            cell=mol.lvec if mol.lvec is not None else None,
            pbc=(mol.lvec is not None)
        )
        return atoms

def main():
    parser = argparse.ArgumentParser(
        description='Plot molecular orbitals from DFTB+ using waveplot'
    )
    parser.add_argument(
        '--molecule',
        type=str,
        help='XYZ file to load (short names: H2O, PTCDA; or full path to .xyz file)'
    )
    parser.add_argument(
        '--use-ase',
        action='store_true',
        help='Use ASE to load XYZ file (default: AtomicSystem)'
    )
    parser.add_argument(
        '--sk-path',
        type=str,
        default='/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/',
        help='Path to Slater-Koster files (default: /home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/)'
    )
    parser.add_argument(
        '--waveplot-exe',
        type=str,
        default='waveplot',
        help='Path to waveplot executable (default: waveplot)'
    )
    parser.add_argument(
        '--sk-wfc-path',
        type=str,
        default='/home/prokophapala/SIMULATIONS/dftbplus/recipes/defect/carbon2d-elec/data/wfc.mio-0-1.hsd',
        help='Path to wavefunction coefficient file (default: /home/prokophapala/SIMULATIONS/dftbplus/recipes/defect/carbon2d-elec/data/wfc.mio-0-1.hsd)'
    )
    parser.add_argument(
        '--plotted-levels',
        type=str,
        default='1:-1',
        help='Which orbitals to plot (default: 1:-1 for all)'
    )
    parser.add_argument(
        '--n-points',
        type=int,
        nargs=3,
        default=[50, 50, 50],
        help='Grid resolution (nx ny nz, default: 50 50 50)'
    )
    parser.add_argument(
        '--workdir',
        type=str,
        default='orbital_calc',
        help='Working directory (default: orbital_calc)'
    )
    parser.add_argument(
        '--skip-dftb',
        action='store_true',
        help='Skip DFTB+ calculation (use existing detailed.xml and eigenvec.bin)'
    )
    parser.add_argument(
        '--skip-waveplot',
        action='store_true',
        help='Skip waveplot (use if WFC files not available)'
    )
    parser.add_argument(
        '--plot-2d',
        action='store_true',
        help='Plot orbitals as 2D projections after waveplot'
    )
    parser.add_argument(
        '--plane',
        type=str,
        default='xy',
        choices=['xy', 'xz', 'yz'],
        help='Projection plane (default: xy)'
    )
    parser.add_argument(
        '--orbital',
        type=str,
        default='all',
        help='Which orbitals to plot: "all", "occupied", "virtual", or comma-separated indices'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default='orbital_plots',
        help='Output directory for PNG plots'
    )
    parser.add_argument(
        '--electrostatic',
        action='store_true',
        help='Calculate and plot electrostatic potential'
    )
    parser.add_argument(
        '--plot-density',
        action='store_true',
        help='Plot electron density (already generated by default)'
    )
    parser.add_argument(
        '--points',
        type=str,
        help='File with points (x y z per line) or comma-separated coordinates'
    )
    parser.add_argument(
        '--grid-points',
        type=int,
        nargs=2,
        default=[20, 20],
        help='Grid resolution for point sampling (nx ny, default: 20 20)'
    )
    
    args = parser.parse_args()
    
    print("=== DFTB+ Molecular Orbital Plotting Example ===\n")
    
    # Load molecule
    if args.molecule:
        # Check if it's a short name like "H2O" or "PTCDA"
        short_names = {
            'H2O': '/home/prokophapala/git/FireCore/cpp/common_resources/xyz/H2O.xyz',
            'PTCDA': '/home/prokophapala/git/FireCore/cpp/common_resources/xyz/PTCDA.xyz'
        }
        if args.molecule in short_names:
            molecule_path = short_names[args.molecule]
            print(f"Loading {args.molecule} from: {molecule_path}")
        elif args.molecule.endswith('.xyz'):
            molecule_path = args.molecule
            print(f"Loading molecule from: {molecule_path}")
        else:
            # Try to find in common_resources
            molecule_path = f'/home/prokophapala/git/FireCore/cpp/common_resources/xyz/{args.molecule}.xyz'
            if os.path.exists(molecule_path):
                print(f"Loading {args.molecule} from: {molecule_path}")
            else:
                molecule_path = args.molecule
                print(f"Loading molecule from: {molecule_path}")
        
        atoms = load_molecule_xyz(molecule_path, use_ase=args.use_ase)
    else:
        print("Using built-in H2O molecule")
        atoms = Atoms('H2O', positions=[
            [0.0, 0.0, 0.0],
            [0.957, 0.0, 0.0],
            [-0.240, 0.927, 0.0]
        ])
    
    print(f"System: {len(atoms)} atoms")
    print(f"Formula: {atoms.get_chemical_formula()}")
    print()
    
    # Create working directory
    workdir = Path(args.workdir).absolute()
    workdir.mkdir(exist_ok=True, parents=True)
    os.chdir(workdir)
    
    # Determine angular momenta from elements
    elem_set = set(atoms.get_chemical_symbols())
    max_angular = {}
    for elem in elem_set:
        if elem in ['H']:
            max_angular[elem] = '"s"'
        elif elem in ['C', 'N', 'O']:
            max_angular[elem] = '"p"'
        elif elem in ['S', 'Si']:
            max_angular[elem] = '"d"'
        else:
            max_angular[elem] = '"p"'
    
    # Write XYZ file
    from pyBall import atomicUtils as au
    au.saveXYZ(es=atoms.get_chemical_symbols(), xyzs=atoms.positions, fname="geo.xyz")
    
    if not args.skip_dftb:
        # Write DFTB+ input with eigenvector output
        print("Writing DFTB+ input file...")
        max_angular_str = '\n'.join([f'        {elem} = {max_angular[elem]}' for elem in max_angular])
        
        hsd_content = f'''Geometry = xyzFormat {{
    <<< "geo.xyz"
}}

Options {{
  WriteDetailedXml = Yes
}}

Analysis {{
  WriteEigenvectors = Yes
}}

Hamiltonian = DFTB {{
  Scc = Yes
  SlaterKosterFiles = Type2FileNames {{
    Prefix = "{args.sk_path}"
    Separator = "-"
    Suffix = ".skf"
  }}
  MaxAngularMomentum {{
{max_angular_str}
  }}
  SCCTolerance = 1.000000e-07
}}
'''
        
        with open('dftb_in.hsd', 'w') as f:
            f.write(hsd_content)
        
        # Run DFTB+
        print("Running DFTB+ calculation...")
        result = subprocess.run(['dftb+'], capture_output=True, text=True)
        
        if result.returncode != 0:
            print("DFTB+ failed!")
            print("STDERR:", result.stderr)
            print("STDOUT:", result.stdout)
            return
        
        print("DFTB+ calculation complete")
        print()
    
    # Check for required files
    if not os.path.exists('detailed.xml'):
        print("Error: detailed.xml not found!")
        print("Run DFTB+ with Options { WriteDetailedXml = Yes }")
        return
    
    if not os.path.exists('eigenvec.bin'):
        print("Error: eigenvec.bin not found!")
        print("Run DFTB+ with Analysis { WriteEigenvectors = Yes }")
        return
    
    # Run waveplot
    if not args.skip_waveplot:
        print("Running waveplot to generate cube files...")
        nx, ny, nz = args.n_points
        
        try:
            cubes = dftbu.run_waveplot(
                workdir='.',
                waveplot_exe=args.waveplot_exe,
                sk_wfc_path=args.sk_wfc_path,
                plotted_levels=args.plotted_levels,
                n_points=(nx, ny, nz),
                resolution=0.01,
                electrostatic_potential=args.electrostatic
            )
            print(f"Generated {len(cubes)} cube files:")
            for cube in cubes:
                print(f"  {cube}")
            print()
        except RuntimeError as e:
            print(f"Waveplot failed: {e}")
            print()
            print("Note: waveplot requires wavefunction coefficient files (wfc.*.hsd)")
            print("These are typically in the Slater-Koster directory or a separate wfc directory")
            print("Use --skip-waveplot to skip this step")
            return
    else:
        print("Skipping waveplot (--skip-waveplot specified)")
        cubes = []
        print()
    
    # Read and display cube file info
    if cubes:
        print("Reading cube files...")
        for cube_file in cubes[:3]:  # Show first 3
            try:
                data, atoms_cube = dftbu.read_cube(cube_file)
                print(f"  {cube_file}: grid shape {data.shape}")
            except Exception as e:
                print(f"  {cube_file}: error reading - {e}")
        print()
        
        # Plot 2D projections if requested
        if args.plot_2d:
            print("Plotting 2D orbital projections...")
            plot_orbitals_2d(cubes, atoms, args)
            print()
        
        # Plot density at specific points if requested
        if args.points:
            print("Evaluating density at specified points...")
            plot_density_at_points(cubes, atoms, args)
            print()
    
    # Return to original directory
    os.chdir('..')
    
    print("=== Example Complete ===")
    print(f"Results saved in: {workdir.absolute()}")
    print()
    print("Cube files can be visualized with:")
    print("  - VMD: vmd wp-*.cube")
    print("  - PyMOL: load wp-*.cube")
    print("  - ASE: from ase.io.cube import read_cube")

if __name__ == '__main__':
    main()
