#!/usr/bin/env python3
"""NEB calculation for H-transfer between N-passivated and NH-passivated ribbons.

Uses manual NEB implementation with DFTB+ via dftb_utils (working approach).
Periodic boundary conditions along x with k-point sampling.
"""
import sys
import os
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("../../")
from pyBall import dftb_utils as dftbu

# Parameters
Lx = 2.4  # Lattice constant
nk_x = 16  # k-points along x
basis_path = '/home/prokophapala/SIMULATIONS/dftbplus/slakos/3ob-3-1/'
n_images = 7  # Number of NEB images (excluding endpoints)
spring_k = 5.0  # Spring constant (eV/Angstrom^2)

def read_gen_format(fname):
    """Read DFTB+ GenFormat file."""
    with open(fname, 'r') as f:
        lines = f.readlines()
    
    # Parse header
    natoms, pbc = lines[0].split()
    natoms = int(natoms)
    enames = lines[1].split()
    
    # Parse atoms
    apos = np.zeros((natoms, 3))
    atypes = []
    for i in range(natoms):
        parts = lines[2+i].split()
        idx = int(parts[1]) - 1  # Convert to 0-indexed
        apos[i, 0] = float(parts[2])
        apos[i, 1] = float(parts[3])
        apos[i, 2] = float(parts[4])
        atypes.append(enames[idx])
    
    # Parse lattice vectors
    origin = np.array([float(x) for x in lines[2+natoms].split()])
    lvs = np.zeros((3, 3))
    for i in range(3):
        lvs[i] = np.array([float(x) for x in lines[3+natoms+i].split()])
    
    return apos, atypes, lvs

def save_xyz(apos, atypes, fname, comment=""):
    """Save geometry to XYZ file."""
    with open(fname, 'w') as f:
        f.write(f"{len(atypes)}\n")
        f.write(f"{comment}\n")
        for i, (elem, pos) in enumerate(zip(atypes, apos)):
            f.write(f"{elem} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
    print(f"Saved: {fname}")

def create_initial_final_states():
    """Create initial and final states for H-transfer NEB.
    
    Initial: H on NH ribbon (as built)
    Final: H moved to N ribbon (forming N-H bond)
    """
    # Read initial geometry
    apos, atypes, lvs = read_gen_format('two_ribbons.gen')
    
    # Find H atoms and N atoms
    h_indices = [i for i, t in enumerate(atypes) if t == 'H']
    n_indices = [i for i, t in enumerate(atypes) if t == 'N']
    
    print(f"Found {len(h_indices)} H atoms at indices: {h_indices}")
    print(f"Found {len(n_indices)} N atoms at indices: {n_indices}")
    
    # Identify which H is on NH ribbon (higher y)
    h_y_positions = [apos[i, 1] for i in h_indices]
    h_nh_idx = h_indices[np.argmax(h_y_positions)]  # H on NH ribbon (top)
    
    # Identify target N on N ribbon (lower y)
    n_y_positions = [apos[i, 1] for i in n_indices]
    n_target_idx = n_indices[np.argmin(n_y_positions)]  # N on N ribbon (bottom)
    
    print(f"H to move: index {h_nh_idx} at y={apos[h_nh_idx, 1]:.3f}")
    print(f"Target N: index {n_target_idx} at y={apos[n_target_idx, 1]:.3f}")
    
    # Save initial state
    apos_initial = apos.copy()
    save_xyz(apos_initial, atypes, 'initial.xyz', "Initial state: H on NH ribbon")
    
    # Create final state (H moved to N ribbon)
    apos_final = apos.copy()
    # Move H from NH ribbon to N ribbon
    apos_final[h_nh_idx] = apos_final[n_target_idx] + np.array([0.0, 1.0, 0.0])
    save_xyz(apos_final, atypes, 'final.xyz', "Final state: H on N ribbon")
    
    return apos_initial, apos_final, atypes, lvs, h_nh_idx, n_target_idx

def run_dftb_calculation(apos, atypes, lvs, temp_dir, do_forces=True):
    """Run DFTB+ calculation using dftb_utils."""
    cwd = os.getcwd()
    os.makedirs(temp_dir, exist_ok=True)
    os.chdir(temp_dir)
    
    try:
        # Write periodic DFTB+ input
        dftbu.makeDFTBjob_pbc(
            enames=atypes, apos=apos, lvs=lvs, fname='dftb_in.hsd',
            basis_path=basis_path,
            nk=(nk_x, 1, 1), k_shift=(0.5, 0.0, 0.0),
            opt=False, params=dftbu.default_params,
            Temperature=600, MixingParameter=0.1,
            MaxScc=500, SCCTolerance=1e-4
        )
        
        # Add force calculation and analysis options
        with open('dftb_in.hsd', 'a') as f:
            f.write("\nAnalysis {\n")
            f.write("  CalculateForces = Yes\n")
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
        
        return E, forces
            
    except Exception as e:
        print(f"  ERROR: {e}")
        return None, None
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

def run_neb():
    """Run NEB calculation using dftb_utils."""
    print("Creating initial and final states...")
    apos_initial, apos_final, atypes, lvs, h_idx, n_idx = create_initial_final_states()
    
    print("\nSetting up NEB images...")
    # Create images via linear interpolation
    n_total = n_images + 2  # including endpoints
    images = linear_interpolate(apos_initial, apos_final, n_total)
    
    # Save initial path
    for i, img in enumerate(images):
        save_xyz(img, atypes, f'neb_image_{i:02d}_init.xyz', f"NEB image {i} initial")
    
    print(f"\nRunning NEB optimization with {n_total} images...")
    print("(Note: Using dftb_utils for DFTB+ calculations)")
    
    # Run initial single-point calculations for all images
    energies = []
    forces_list = []
    
    for i, img in enumerate(images):
        print(f"\nCalculating image {i}/{n_total-1}...")
        temp_dir = f'temp_neb_{i}'
        E, forces = run_dftb_calculation(img, atypes, lvs, temp_dir, do_forces=True)
        
        if E is None:
            print(f"  FAILED for image {i}")
            return
        
        energies.append(E)
        forces_list.append(forces)
        print(f"  E = {E:.4f} eV")
    
    # Simple optimization loop (no climbing image for now)
    max_steps = 50
    dt = 0.1  # Time step for steepest descent
    
    for step in range(max_steps):
        print(f"\n--- Optimization step {step+1}/{max_steps} ---")
        
        # Compute NEB forces
        neb_forces = compute_neb_forces(images, atypes, lvs, energies, forces_list)
        
        # Check convergence
        max_force = max([np.linalg.norm(f) for f in neb_forces[1:-1]])
        print(f"Max force: {max_force:.4f} eV/Ang")
        
        if max_force < 0.05:
            print("Converged!")
            break
        
        # Update images (steepest descent, excluding endpoints)
        for i in range(1, n_total - 1):
            images[i] -= dt * neb_forces[i]
        
        # Re-calculate energies and forces
        new_energies = []
        new_forces = []
        for i, img in enumerate(images):
            temp_dir = f'temp_neb_{i}_step{step}'
            E, forces = run_dftb_calculation(img, atypes, lvs, temp_dir, do_forces=True)
            new_energies.append(E)
            new_forces.append(forces)
        
        energies = new_energies
        forces_list = new_forces
        
        # Print energies
        for i, E in enumerate(energies):
            print(f"  Image {i}: E = {E:.4f} eV")
    
    # Save final images
    for i, img in enumerate(images):
        save_xyz(img, atypes, f'neb_image_{i:02d}_final.xyz', f"NEB image {i} final, E={energies[i]:.4f} eV")
    
    # Plot results
    plot_neb_results(energies)
    
    print("\nNEB calculation complete!")
    print("Output files: neb_image_*.xyz, neb_energy_profile.png")

def plot_neb_results(energies):
    """Plot NEB energy profile."""
    distances = np.arange(len(energies))
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(distances, energies, 'o-', markersize=8, linewidth=2)
    ax.set_xlabel('Image number')
    ax.set_ylabel('Energy (eV)')
    ax.set_title('NEB Energy Profile for H-Transfer')
    ax.grid(True, alpha=0.3)
    
    # Mark barrier height
    E_min = min(energies[0], energies[-1])
    E_max = max(energies)
    barrier = E_max - E_min
    ax.axhline(y=E_max, color='r', linestyle='--', alpha=0.5, label=f'Barrier: {barrier:.3f} eV')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig('neb_energy_profile.png', dpi=150)
    print("\nSaved NEB energy profile: neb_energy_profile.png")
    plt.close()

if __name__ == '__main__':
    print("=" * 60)
    print("NEB calculation for H-transfer between ribbons")
    print("=" * 60)
    print(f"Lx = {Lx} Å")
    print(f"k-points along x: {nk_x}")
    print(f"Number of images: {n_images}")
    print(f"Basis path: {basis_path}")
    print("=" * 60)
    
    run_neb()
