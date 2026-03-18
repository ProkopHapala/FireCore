import os
import shutil
import subprocess
from ase.io import read, write
from ase import Atoms
from ase.constraints import FixAtoms
from ase.calculators.aims import Aims
import numpy as np
from ase import units

# Specify the path to your configuration file (replace with your actual path)
config_file_path = 'au111_input.extxyz'

# Read the configurations from the extxyz file
configs = read(config_file_path, ':')

# Create a folder for aims_out_files
aims_out_folder = 'relaxation_files'
os.makedirs(aims_out_folder, exist_ok=True)

# DFT results file (append mode, so no need to create it explicitly)
dft_results_file = 'dft_relaxations.extxyz'

# Checkpoint file to track progress
checkpoint_file = 'dft_checkpoint.txt'
if os.path.exists(checkpoint_file):
    with open(checkpoint_file, 'r') as f:
        last_completed_index = int(f.read().strip()) + 1
else:
    last_completed_index = 0

# Define the maximum number of substrate atoms to fix
FixIndex = 50  # Example: fix the bottom 49 atoms of the substrate

# Loop over each configuration and perform DFT calculation
for i, config in enumerate(configs[last_completed_index:], start=last_completed_index):
    # Apply constraints if FixIndex is defined
    if FixIndex is not None:
        fixed_indices = list(range(FixIndex))  # Indices of atoms to fix
        fix_constraint = FixAtoms(indices=fixed_indices)
        config.set_constraint(fix_constraint)
        
    calc = Aims(
        aims_command="srun /scratch/project_2008059/FHI-aims/build_211214/aims.211214.x",
        species_dir="/scratch/project_2008059/FHI-aims/AIMS_SPECIES/defaults_2020/tight/",
        xc="pbe",
        charge=0,
        spin="none",
        k_grid=[1, 1, 1],
        # mixer="pulay",
        # n_max_pulay=8,
        relativistic=["atomic_zora", "scalar"],
        relax_geometry=["trm", 0.01],
        compute_forces=True,
        occupation_type=["gaussian", 0.01],
        #output=['dipole'],
        vdw_ts='vdw_params_kind=tssurf',
        #sc_accuracy_eev=1e-3,
        #sc_accuracy_rho=1e-5,
        #sc_accuracy_etot=1e-6,
        #sc_accuracy_forces=5e-4,
        sc_iter_limit=400,
        compensate_multipole_errors=True,
        load_balancing=True
    )

    # Set the calculator for the configuration
    config.set_calculator(calc)

    try:
        # Run DFT calculation
        energy = config.get_potential_energy()
        forces = config.get_forces()
        positions = config.get_positions()
        #dipole = config.get_dipole_moment()

        # Save the aims.out file with a unique name
        aims_out_file = os.path.join(aims_out_folder, f'aims_{i+1:03d}.out')
        shutil.copy('aims.out', aims_out_file)

        # Append the DFT results to the extxyz file
        with open(dft_results_file, 'a') as dft_file:
            write(dft_file, config, format='extxyz')

        # Update the checkpoint file
        with open(checkpoint_file, 'w') as f:
            f.write(str(i))

    except Exception as e:
        print(f"Calculation did not converge for configuration {i+1}: {e}")
        # Save the aims.out file with a unique name even if the calculation failed
        aims_out_file = os.path.join(aims_out_folder, f'aims_{i+1:03d}.out')
        if os.path.exists('aims.out'):
            shutil.copy('aims.out', aims_out_file)
        else:
            print(f"No aims.out file found for configuration {i+1}")


# filter structures based on their energy
# to remove the very high energy outliers
# from the initial data

# Load all the structures from the XYZ file
#structures = read('dft_results.xyz', index=':')

#print(len(structures))
# Extract the energy of the first (optimized) structure
#optimized_energy = structures[0].get_potential_energy()

# Filter structures based on the energy difference with the optimized structure
#filtered_structures = [structure for structure in structures if abs(structure.get_potential_energy() - optimized_energy) < 0.3]

# Save the filtered structures to a new XYZ file
#write('dft_results_filtered.xyz', filtered_structures)

#print(len(filtered_structures))
