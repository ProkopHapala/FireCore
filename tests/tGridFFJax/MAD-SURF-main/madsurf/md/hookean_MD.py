import os
import subprocess
from ase import Atoms
from ase.io import read
from ase.io.aims import read_aims
from tblite.ase import TBLite
from ase.optimize import BFGS
from ase.md.verlet import VelocityVerlet
from ase.md.velocitydistribution import (MaxwellBoltzmannDistribution)
from ase.units import fs
from ase.constraints import Hookean, FixAtoms
import numpy as np
import argparse

def get_adsorbate_neighbors(atoms, metal_surface_size, molecule_size, threshold=2.0):
    metal_surface = atoms[:metal_surface_size]
    adsorbates = atoms[metal_surface_size:]
    closest_bond_indices_per_adsorbate = []
    global_index_offset = metal_surface_size

    for i in range(0, len(adsorbates), molecule_size):
        current_adsorbate = adsorbates[i:i + molecule_size]
        connectivity_matrix = current_adsorbate.get_all_distances(mic=True)
        closest_bond_indices = []

        for atom_idx in range(len(current_adsorbate)):
            global_atom_idx = atom_idx + i + global_index_offset
            neighbors_mask = np.where(connectivity_matrix[atom_idx] < threshold)[0]
            neighbors_mask = neighbors_mask[neighbors_mask != atom_idx]

            for neighbor_idx in neighbors_mask:
                global_neighbor_idx = neighbor_idx + i + global_index_offset
                distance = connectivity_matrix[atom_idx, neighbor_idx]
                atom_type = current_adsorbate[atom_idx].symbol
                neighbor_type = current_adsorbate[neighbor_idx].symbol
                pair = tuple(sorted([global_atom_idx, global_neighbor_idx]))
                closest_bond_indices.append((pair, atom_type, neighbor_type, distance))

        closest_bond_indices_per_adsorbate.append(closest_bond_indices)

    return closest_bond_indices_per_adsorbate


def apply_constraints(atoms, closest_bond_indices_per_adsorbate, metal_surface_size, fix_slab):
    constraints = []

    # Fix the metal slab if fix_slab is set
    if fix_slab:
        metal_indices = list(range(metal_surface_size))
        slab_constraint = FixAtoms(indices=metal_indices)
        constraints.append(slab_constraint)

    bond_types = {
        'C-C': 9.0,
        'C-H': 7.0,
        'C-O': 6.0,
        'O-H': 5.0,
        'C-N': 7.0,
        'O-N': 7.0
    }

    for bond_type, spring_constant in bond_types.items():
        atom1, atom2 = bond_type.split('-')
        for i, closest_bond_indices in enumerate(closest_bond_indices_per_adsorbate):
            for bond in closest_bond_indices:
                pair = bond[0]
                atom_type, neighbor_type = bond[1], bond[2]

                if (atom_type == atom1 and neighbor_type == atom2) or (atom_type == atom2 and neighbor_type == atom1):
                    idx1, idx2 = pair
                    length = atoms.get_distance(idx1, idx2) + 0.5
                    c = Hookean(a1=int(idx1), a2=int(idx2), rt=float(length), k=spring_constant)
                    constraints.append(c)

    atoms.set_constraint(constraints)


def main():
    parser = argparse.ArgumentParser(description="Molecular Dynamics Simulation with TBLite")
    parser.add_argument('input_file', type=str, help="Path to the input geometry file")
    parser.add_argument('--metal_size', type=int, default=192, help="Number of atoms in the metal surface")
    parser.add_argument('--molecule_size', type=int, default=35, help="Number of atoms in each adsorbate")
    parser.add_argument('--threshold', type=float, default=2.0, help="Distance threshold for neighbor detection")
    parser.add_argument('--temp', type=float, default=300, help="Temperature in Kelvin")
    parser.add_argument('--timestep', type=float, default=1.0, help="Timestep in femtoseconds for MD simulation")
    parser.add_argument('--steps', type=int, default=100, help="Number of steps for geometry optimization")
    parser.add_argument('--md_steps', type=int, default=5000, help="Number of steps for MD simulation")
    parser.add_argument('--fix_slab', action='store_true', help="Fix the metal surface atoms")

    args = parser.parse_args()

    # Create an ASE Atoms object with the initial geometry
    atoms = read_aims(args.input_file)
    metal_surface_size = args.metal_size
    molecule_size = args.molecule_size

    # Get closest bond indices for each adsorbate
    closest_bond_indices_per_adsorbate = get_adsorbate_neighbors(atoms, metal_surface_size, molecule_size, args.threshold)

    # Apply Hookean constraints to the bonds and fix slab if --fix_slab is set
    apply_constraints(atoms, closest_bond_indices_per_adsorbate, metal_surface_size, args.fix_slab)

    # Set up the TBLite calculator
    tblite_calculator = TBLite(method="GFN1-xTB")
    atoms.set_calculator(tblite_calculator)

    # Optimize geometry
    save_file = os.path.splitext(args.input_file)[0]
    optimizer = BFGS(atoms, trajectory=f"{save_file}_tblite_opt.traj")
    optimizer.run(fmax=0.005, steps=args.steps)

    # Set momenta corresponding to the specified temperature
    MaxwellBoltzmannDistribution(atoms, temperature_K=args.temp)

    # Set up the Molecular Dynamics run
    dyn = VelocityVerlet(atoms, timestep=args.timestep * fs,
                         trajectory=f"{save_file}_md_tblite.traj", logfile=f"{save_file}_md_tblite.log")
    dyn.run(args.md_steps)


if __name__ == "__main__":
    main()
