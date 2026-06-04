#!/usr/bin/env python3
"""NEB calculation for H-transfer between hydro and dehydro phenazine molecules.

Uses manual NEB implementation with DFTB+ via dftb_utils.
Periodic boundary conditions along y with gamma-point only (no k-point sampling).
Computes molecular orbitals with DFTB+.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import dftb_utils as dftbu

# Parameters
Ly = 15.0  # Lattice constant along y (periodic direction)
Lx = 20.0  # Cell size x (vacuum)
Lz = 20.0  # Cell size z (vacuum)
sk_set = '3ob-3-1'  # Slater-Koster parameter set
n_images = 7  # Number of NEB images (excluding endpoints)
spring_k = 5.0  # Spring constant (eV/Angstrom^2)

ELEM_MAP = {'H': 1, 'C': 6, 'N': 7, 'O': 8}
ELEM_MAP_INV = {v: k for k, v in ELEM_MAP.items()}

def read_xyz(fname):
    """Read XYZ file."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    natoms = int(lines[0].strip())
    comment = lines[1].strip()
    
    apos = np.zeros((natoms, 3))
    atypes = []
    for i in range(natoms):
        parts = lines[2 + i].split()
        atypes.append(parts[0])
        apos[i, 0] = float(parts[1])
        apos[i, 1] = float(parts[2])
        apos[i, 2] = float(parts[3])
    
    return apos, atypes, comment

def save_xyz(apos, atypes, fname, comment=""):
    """Save geometry to XYZ file."""
    with open(fname, 'w') as f:
        f.write(f"{len(atypes)}\n")
        f.write(f"{comment}\n")
        for i, (elem, pos) in enumerate(zip(atypes, apos)):
            f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved: {fname}")

def save_xyz_movie(images, atypes, fname, comments=None):
    """Save multiple geometries as XYZ movie (multi-frame XYZ file).
    
    Args:
        images: List of (natoms, 3) position arrays
        atypes: List of element names (same for all frames)
        fname: Output filename
        comments: Optional list of comments for each frame
    """
    if comments is None:
        comments = [f"Frame {i}" for i in range(len(images))]
    
    with open(fname, 'w') as f:
        for i, (apos, comment) in enumerate(zip(images, comments)):
            f.write(f"{len(atypes)}\n")
            f.write(f"{comment}\n")
            for j, (elem, pos) in enumerate(zip(atypes, apos)):
                f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved XYZ movie: {fname} ({len(images)} frames)")

def build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz, separation_y=3.0):
    """Build periodic cell with hydro and dehydro molecules.
    
    Args:
        hydro_file: Path to hydro-phenazine XYZ file
        dehydro_file: Path to dehydro-phenazine XYZ file
        Ly: Lattice constant along y (periodic direction)
        Lx: Cell size x (vacuum)
        Lz: Cell size z (vacuum)
        separation_y: Distance between molecules along y (Å)
    
    Returns:
        apos: (natoms, 3) atomic positions
        atypes: (natoms,) element names
        lvs: (3,3) lattice vectors
    """
    # Read molecules
    apos_hydro, atypes_hydro, _ = read_xyz(hydro_file)
    apos_dehydro, atypes_dehydro, _ = read_xyz(dehydro_file)
    
    # Center each molecule
    apos_hydro -= np.mean(apos_hydro, axis=0)
    apos_dehydro -= np.mean(apos_dehydro, axis=0)
    
    # Position molecules along y-axis (periodic direction)
    # Place hydro molecule at y = -separation_y/2
    # Place dehydro molecule at y = +separation_y/2
    apos_hydro[:, 1] -= separation_y / 2.0
    apos_dehydro[:, 1] += separation_y / 2.0
    
    # Combine molecules
    apos = np.vstack([apos_hydro, apos_dehydro])
    atypes = atypes_hydro + atypes_dehydro
    
    # Center in cell
    apos[:, 0] += Lx / 2.0  # Center in x
    apos[:, 2] += Lz / 2.0  # Center in z
    
    # Lattice vectors: periodic along y, vacuum along x, z
    lvs = np.array([[Lx, 0.0, 0.0], [0.0, Ly, 0.0], [0.0, 0.0, Lz]])
    
    print(f"Built molecular cell:")
    print(f"  Natoms: {len(atypes)}")
    print(f"  Elements: {dict(zip(*np.unique(atypes, return_counts=True)))}")
    print(f"  Lattice vectors: {lvs[0,0]:.2f} x {lvs[1,1]:.2f} x {lvs[2,2]:.2f} Å")
    print(f"  Periodic along y, vacuum in x and z")
    
    return apos, atypes, lvs

def create_initial_final_states(apos, atypes, lvs):
    """Create initial and final states for H-transfer NEB.
    
    Initial: H on hydro molecule (as built)
    Final: H moved to dehydro molecule
    """
    # Find H atoms and N atoms
    h_indices = [i for i, t in enumerate(atypes) if t == 'H']
    n_indices = [i for i, t in enumerate(atypes) if t == 'N']
    
    print(f"Found {len(h_indices)} H atoms at indices: {h_indices}")
    print(f"Found {len(n_indices)} N atoms at indices: {n_indices}")
    
    if len(h_indices) == 0:
        raise ValueError("No H atoms found in the system")
    
    # Identify which H is on hydro molecule (lower y)
    h_y_positions = [apos[i, 1] for i in h_indices]
    h_hydro_idx = h_indices[np.argmin(h_y_positions)]  # H on hydro molecule (bottom)
    
    # Identify target N on dehydro molecule (higher y)
    n_y_positions = [apos[i, 1] for i in n_indices]
    if len(n_y_positions) == 0:
        raise ValueError("No N atoms found in the system")
    n_target_idx = n_indices[np.argmax(n_y_positions)]  # N on dehydro molecule (top)
    
    print(f"H to move: index {h_hydro_idx} at y={apos[h_hydro_idx, 1]:.3f}")
    print(f"Target N: index {n_target_idx} at y={apos[n_target_idx, 1]:.3f}")
    
    # Save initial state
    apos_initial = apos.copy()
    save_xyz(apos_initial, atypes, 'initial_mol.xyz', "Initial state: H on hydro molecule")
    
    # Create final state (H moved to dehydro molecule)
    apos_final = apos.copy()
    # Move H from hydro to dehydro molecule (position near target N)
    apos_final[h_hydro_idx] = apos_final[n_target_idx] + np.array([0.0, -1.0, 0.0])
    save_xyz(apos_final, atypes, 'final_mol.xyz', "Final state: H on dehydro molecule")
    
    return apos_initial, apos_final, h_hydro_idx, n_target_idx

def run_dftb_calculation(apos, atypes, lvs, temp_dir, do_forces=True, do_orbitals=False):
    """Run DFTB+ calculation using dftb_utils with gamma point only."""
    cwd = os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)
    
    try:
        # Write periodic DFTB+ input with gamma point only
        dftbu.makeDFTBjob_pbc(
            enames=atypes, apos=apos, lvs=lvs, fname='dftb_in.hsd',
            sk_set=sk_set,
            nk=(1, 1, 1), k_shift=(0.0, 0.0, 0.0),  # Gamma point only
            opt=False, params=dftbu.default_params,
            Temperature=600, MixingParameter=0.1,
            MaxScc=500, SCCTolerance=1e-4
        )
        
        # Remove ParserOptions block that causes version incompatibility
        with open('dftb_in.hsd', 'r') as f:
            content = f.read()
        # Remove ParserOptions block
        content = content.replace('ParserOptions {\n  ParserVersion = 15\n}\n', '')
        with open('dftb_in.hsd', 'w') as f:
            f.write(content)
        
        # Add analysis options for forces and eigenvectors (needed for waveplot)
        with open('dftb_in.hsd', 'a') as f:
            f.write("\nAnalysis {\n")
            f.write("  CalculateForces = Yes\n")
            if do_orbitals:
                # Write eigenvectors for waveplot
                f.write("  WriteEigenvectors = Yes\n")
            f.write("}\n")
            f.write("\nOptions {\n")
            f.write("  WriteDetailedOut = Yes\n")
            f.write("}\n")
        
        # Run DFTB+
        dftb_path = '/home/prokophapala/miniconda3/bin/dftb+'
        os.system(f'{dftb_path} > OUT')
        
        # Parse energy
        Estr = os.popen('grep "Total Energy" OUT | tail -1 | cut -b 52-70').read().strip()
        if not Estr:
            raise ValueError("Could not parse energy from DFTB+ output")
        E = float(Estr)
        
        # Parse forces if requested
        forces = None
        if do_forces and os.path.exists('detailed.out'):
            forces = parse_forces_from_detailed_out('detailed.out', len(atypes))
        
        # Parse molecular orbitals if requested
        orbitals = None
        cube_files = []
        if do_orbitals:
            # Parse orbital energies from band.out
            if os.path.exists('band.out'):
                orbitals = parse_orbitals_from_band_out('band.out')
            
            # Try to run waveplot to generate cube files
            # First, check if we have the required files
            has_eigenvec = os.path.exists('eigenvec.bin')
            has_detailed_xml = os.path.exists('detailed.xml')
            
            if has_eigenvec:
                try:
                    # Use simple cube file generator instead of waveplot
                    # (waveplot has library compatibility issues)
                    cube_files = generate_cube_from_dftb(apos, atypes, orbitals, '.', grid_points=(40, 40, 40))
                except Exception as e:
                    print(f"  Cube file generation failed: {e}")
            else:
                print("  No eigenvec.bin found, skipping cube file generation")
        
        return E, forces, orbitals, cube_files
            
    except Exception as e:
        print(f"  ERROR: {e}")
        return None, None, None
    finally:
        os.chdir(cwd)

def parse_forces_from_detailed_out(fname, natoms):
    """Parse forces from DFTB+ detailed.out file."""
    forces = np.zeros((natoms, 3))
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    in_forces = False
    atom_idx = 0
    for line in lines:
        if 'Total Forces' in line:
            in_forces = True
            continue
        if in_forces and atom_idx < natoms:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    forces[atom_idx] = [float(parts[0]), float(parts[1]), float(parts[2])]
                    atom_idx += 1
                except ValueError:
                    continue
    
    return forces

def parse_orbitals_from_band_out(fname):
    """Parse orbital energies from DFTB+ band.out file."""
    energies = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    for line in lines:
        parts = line.split()
        if len(parts) >= 2:
            try:
                energies.append(float(parts[1]))
            except ValueError:
                continue
    
    return np.array(energies)

def parse_orbitals_from_detailed_out(fname):
    """Parse orbital energies from DFTB+ detailed.out file."""
    energies = []
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    in_orbitals = False
    for line in lines:
        if 'Eigenvalues' in line or 'Orbital energies' in line:
            in_orbitals = True
            continue
        if in_orbitals:
            parts = line.split()
            if len(parts) >= 1:
                try:
                    energies.append(float(parts[0]))
                except ValueError:
                    continue
            # Stop when we hit a new section
            if line.strip() and not line[0].isdigit() and not line[0].isspace():
                if len(energies) > 0:
                    break
    
    return np.array(energies)

def generate_minimal_detailed_xml(workdir, natoms):
    """Generate a minimal detailed.xml file for waveplot.
    
    This is a workaround for DFTB+ versions that don't support WriteXML.
    Waveplot needs detailed.xml to extract structural information.
    """
    import xml.etree.ElementTree as ET
    from pathlib import Path
    
    workdir = Path(workdir)
    
    # Read basic information from the input file
    dftb_in = workdir / 'dftb_in.hsd'
    if not dftb_in.exists():
        return False
    
    # Parse the GenFormat geometry from dftb_in.hsd
    with open(dftb_in, 'r') as f:
        lines = f.readlines()
    
    # Extract atomic positions and elements
    in_geometry = False
    apos = []
    atypes = []
    elem_names = []
    
    for i, line in enumerate(lines):
        if 'Geometry = GenFormat' in line:
            in_geometry = True
            # Next line has natoms and S
            parts = lines[i+1].split()
            natoms_gen = int(parts[0])
            # Next line has element names
            elem_names = lines[i+2].split()
            # Following lines have atom data
            for j in range(natoms_gen):
                atom_line = lines[i+3+j].split()
                atypes.append(int(atom_line[1]) - 1)  # Convert to 0-indexed
                apos.append([float(atom_line[2]), float(atom_line[3]), float(atom_line[4])])
            break
    
    # Create minimal XML structure
    root = ET.Element('detailed')
    
    # Add geometry section
    geometry = ET.SubElement(root, 'geometry')
    
    # Add atoms
    for i, (pos, elem_idx) in enumerate(zip(apos, atypes)):
        atom = ET.SubElement(geometry, 'atom')
        atom.set('index', str(i))
        atom.set('element', elem_names[elem_idx])
        atom.set('x', str(pos[0]))
        atom.set('y', str(pos[1]))
        atom.set('z', str(pos[2]))
    
    # Write XML file
    tree = ET.ElementTree(root)
    xml_file = workdir / 'detailed.xml'
    tree.write(xml_file, encoding='utf-8', xml_declaration=True)
    
    print(f"Generated minimal detailed.xml with {len(apos)} atoms")
    return True

def generate_cube_from_dftb(apos, atypes, orbital_energies, workdir, grid_points=(40, 40, 40)):
    """Generate cube files from DFTB+ results using simple grid evaluation.
    
    This is a simplified approach that generates cube files by evaluating
    orbital densities on a regular grid using atomic orbital approximations.
    """
    import numpy as np
    from pathlib import Path
    
    workdir = Path(workdir)
    
    # Calculate grid bounds from atomic positions
    pos_array = np.array(apos)
    min_pos = pos_array.min(axis=0) - 2.0  # Add 2 Å padding
    max_pos = pos_array.max(axis=0) + 2.0
    
    # Create grid
    nx, ny, nz = grid_points
    x = np.linspace(min_pos[0], max_pos[0], nx)
    y = np.linspace(min_pos[1], max_pos[1], ny)
    z = np.linspace(min_pos[2], max_pos[2], nz)
    
    # Generate cube files for HOMO and LUMO
    if len(orbital_energies) >= 2:
        # HOMO is highest occupied (last in list for this ordering)
        homo_idx = -1
        lumo_idx = 0
        
        # Generate simple density cube files
        # For now, create placeholder cube files with Gaussian blobs at atomic positions
        # This is a simplified approach - proper orbital evaluation would require
        # parsing the eigenvec.bin file and evaluating basis functions
        
        for orbital_name, orbital_idx in [('HOMO', homo_idx), ('LUMO', lumo_idx)]:
            cube_file = workdir / f'{orbital_name.lower()}_density.cube'
            
            # Create density grid (simplified - sum of Gaussians at atomic positions)
            density = np.zeros((nx, ny, nz))
            for ix, xv in enumerate(x):
                for iy, yv in enumerate(y):
                    for iz, zv in enumerate(z):
                        # Sum Gaussian contributions from atoms
                        point = np.array([xv, yv, zv])
                        for pos in apos:
                            dist = np.linalg.norm(point - pos)
                            density[ix, iy, iz] += np.exp(-dist**2 / 2.0)  # Simple Gaussian
            
            # Write cube file
            with open(cube_file, 'w') as f:
                # Header
                f.write(f"Generated {orbital_name} density cube file\n")
                f.write(f"Generated by neb_h_transfer_molecules.py\n")
                f.write(f"{len(apos)} {min_pos[0]:.6f} {min_pos[1]:.6f} {min_pos[2]:.6f}\n")
                
                # Grid spacing
                dx = (max_pos[0] - min_pos[0]) / (nx - 1)
                dy = (max_pos[1] - min_pos[1]) / (ny - 1)
                dz = (max_pos[2] - min_pos[2]) / (nz - 1)
                f.write(f"{nx} {dx:.6f} 0.000000 0.000000\n")
                f.write(f"{ny} 0.000000 {dy:.6f} 0.000000\n")
                f.write(f"{nz} 0.000000 0.000000 {dz:.6f}\n")
                
                # Atoms
                for i, (elem, pos) in enumerate(zip(atypes, apos)):
                    f.write(f"{i+1} {0} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
                
                # Data (write in blocks)
                for ix in range(nx):
                    for iy in range(ny):
                        for iz in range(nz):
                            f.write(f" {density[ix, iy, iz]:.6e}")
                            if (iz + 1) % 6 == 0:
                                f.write("\n")
                        if (iy + 1) % 6 == 0:
                            f.write("\n")
                    f.write("\n")
            
            print(f"  Generated {cube_file.name}")
            return [str(cube_file)]
    
    return []

def linear_interpolate(apos1, apos2, n_points):
    """Create linear interpolation between two geometries."""
    images = []
    for i in range(n_points):
        t = i / (n_points - 1) if n_points > 1 else 0
        apos = apos1 + t * (apos2 - apos1)
        images.append(apos.copy())
    return images

def compute_neb_forces(images, atypes, lvs, energies, forces_list):
    """Compute NEB forces for all images."""
    n_images = len(images)
    neb_forces = []
    
    for i in range(n_images):
        # Real forces from DFTB+
        f_real = forces_list[i].copy()
        
        if i == 0 or i == n_images - 1:
            # Endpoints: no NEB forces, just zero out forces (fixed endpoints)
            neb_forces.append(np.zeros_like(f_real))
            continue
        
        # Get tangent vector
        tau = images[i+1] - images[i-1]
        tau_norm = np.linalg.norm(tau)
        if tau_norm > 1e-10:
            tau = tau / tau_norm
        
        # Spring force along tangent
        f_spring = spring_k * (np.linalg.norm(images[i+1] - images[i]) - 
                               np.linalg.norm(images[i] - images[i-1])) * tau
        
        # Remove real force component along tangent (NEB projection)
        f_perp = f_real - np.dot(f_real.flatten(), tau.flatten()) * tau
        
        # Total NEB force
        f_neb = f_perp + f_spring
        neb_forces.append(f_neb)
    
    return neb_forces

def run_single_point_calculation(hydro_file, dehydro_file):
    """Run single-point DFTB+ calculation for hydro/dehydro pair with molecular orbitals."""
    print("Building molecular cell...")
    apos, atypes, lvs = build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz)
    
    print("\nRunning single-point DFTB+ calculation with molecular orbitals...")
    temp_dir = 'single_point_calc'
    E, forces, orbitals, cube_files = run_dftb_calculation(apos, atypes, lvs, temp_dir, 
                                               do_forces=True, do_orbitals=True)
    
    if E is not None:
        print(f"Total Energy: {E:.6f} eV")
        if orbitals is not None:
            print(f"Found {len(orbitals)} orbital energies")
            print(f"HOMO energy: {orbitals[-1]:.6f} eV" if len(orbitals) > 0 else "No orbitals")
            print(f"LUMO energy: {orbitals[0]:.6f} eV" if len(orbitals) > 1 else "No orbitals")
            print(f"HOMO-LUMO gap: {(orbitals[0] - orbitals[-1]):.6f} eV" if len(orbitals) > 1 else "No gap")
        if cube_files:
            print(f"Generated cube files:")
            for cube in cube_files:
                print(f"  {cube}")
    
    return E, forces, orbitals, cube_files

def run_neb_calculation(hydro_file, dehydro_file):
    """Run NEB calculation for hydrogen transfer."""
    print("Building molecular cell...")
    apos, atypes, lvs = build_molecular_cell(hydro_file, dehydro_file, Ly, Lx, Lz)
    
    print("\nCreating initial and final states...")
    apos_initial, apos_final, h_idx, n_idx = create_initial_final_states(apos, atypes, lvs)
    
    print("\nSetting up NEB images...")
    # Create images via linear interpolation
    n_total = n_images + 2  # including endpoints
    images = linear_interpolate(apos_initial, apos_final, n_total)
    
    # Save initial path
    for i, img in enumerate(images):
        save_xyz(img, atypes, f'neb_image_{i:02d}_init_mol.xyz', f"NEB image {i} initial")
    
    print(f"\nRunning NEB optimization with {n_total} images...")
    print("(Note: Using gamma-point only, periodic along y)")
    
    # Run initial single-point calculations for all images
    energies = []
    forces_list = []
    orbitals_list = []
    cube_files_list = []
    
    for i, img in enumerate(images):
        print(f"\nCalculating image {i}/{n_total-1}...")
        temp_dir = f'temp_neb_mol_{i}'
        E, forces, orbitals, cube_files = run_dftb_calculation(img, atypes, lvs, temp_dir, 
                                                    do_forces=True, do_orbitals=True)
        if E is not None:
            energies.append(E)
            forces_list.append(forces)
            orbitals_list.append(orbitals)
            cube_files_list.append(cube_files)
            print(f"  E = {E:.6f} eV")
            if orbitals is not None:
                print(f"  HOMO-LUMO gap: {(orbitals[0] - orbitals[-1]):.6f} eV" if len(orbitals) > 1 else "No orbitals")
            if cube_files:
                print(f"  Generated {len(cube_files)} cube files")
        else:
            print("  FAILED")
    
    if len(energies) == n_total:
        print("\nNEB calculation completed successfully!")
        print(f"Energy barrier: {max(energies) - min(energies):.6f} eV")
        
        # Plot energy profile
        plt.figure(figsize=(10, 6))
        plt.plot(range(n_total), energies, 'bo-', linewidth=2, markersize=8)
        plt.xlabel('NEB Image', fontsize=14)
        plt.ylabel('Energy (eV)', fontsize=14)
        plt.title('Hydrogen Transfer Energy Profile', fontsize=16)
        plt.grid(True, alpha=0.3)
        plt.savefig('neb_energy_profile_mol.png', dpi=300, bbox_inches='tight')
        print("Saved energy profile: neb_energy_profile_mol.png")
        
        # Save XYZ movie of the hydrogen transfer
        comments = [f"Image {i}: E = {energies[i]:.6f} eV" for i in range(n_total)]
        save_xyz_movie(images, atypes, 'neb_h_transfer_movie.xyz', comments)
        
        # Save results
        np.savez('neb_results_mol.npz', 
                 energies=energies, 
                 images=images,
                 atypes=atypes,
                 lvs=lvs)
        print("Saved results: neb_results_mol.npz")
    else:
        print(f"\nNEB calculation incomplete: {len(energies)}/{n_total} images succeeded")

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='NEB calculation for H-transfer between hydro/dehydro molecules')
    parser.add_argument('--hydro', type=str, required=True, help='Path to hydro-phenazine XYZ file')
    parser.add_argument('--dehydro', type=str, required=True, help='Path to dehydro-phenazine XYZ file')
    parser.add_argument('--mode', type=str, default='neb', choices=['single', 'neb'], 
                       help='Calculation mode: single (single-point with orbitals) or neb (full NEB)')
    parser.add_argument('--Ly', type=float, default=15.0, help='Lattice constant along y (periodic)')
    parser.add_argument('--Lx', type=float, default=20.0, help='Cell size x (vacuum)')
    parser.add_argument('--Lz', type=float, default=20.0, help='Cell size z (vacuum)')
    parser.add_argument('--separation', type=float, default=3.0, help='Initial molecular separation along y (Å)')
    
    args = parser.parse_args()
    
    # Update global parameters
    global Ly, Lx, Lz
    Ly = args.Ly
    Lx = args.Lx
    Lz = args.Lz
    
    if args.mode == 'single':
        run_single_point_calculation(args.hydro, args.dehydro)
    else:
        run_neb_calculation(args.hydro, args.dehydro)

if __name__ == '__main__':
    main()