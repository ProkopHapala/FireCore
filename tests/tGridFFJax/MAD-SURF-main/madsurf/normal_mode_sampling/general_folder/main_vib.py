import os
import numpy as np
import shutil
from ase import Atoms
from ase.optimize import BFGS
from ase.io import read, write
from ase.io.trajectory import Trajectory
from ase.calculators.aims import Aims
from ase.vibrations import Infrared
from ase.vibrations import Vibrations

# Read the input geometry file
file_CO = read('geometry_opt.in', format='aims')

aims_command="srun /scratch/project_2008059/FHI-aims/build_211214/aims.211214.x aims.211214.x "
species_dir="/scratch/project_2008059/FHI-aims/AIMS_SPECIES/defaults_2020/tight/"

# Define the AIMS calculator
calc = Aims(xc='pbe', relativistic='atomic_zora scalar', spin='none',
            #sc_accuracy_etot=1e-6, sc_accuracy_rho=1e-5,
            #sc_accuracy_eev=1e-3, sc_accuracy_forces=5e-4,
            charge=0, vdw_ts='vdw_params_kind=tssurf',
            sc_iter_limit=300,
            compute_forces=True,
            #k_grid = (12,12,1),
            #charge= -0.5,
            species_dir=species_dir,
            aims_command =aims_command)

# Assign the calculator to the initial atoms object
file_CO.set_calculator(calc)

ir = Vibrations(file_CO, delta=0.002)
ir.run()
ir.summary()
ir.write_jmol()
ir.write_mode()
ir.write_dos(out='vib-dos.dat', start=50, end=4000, npts=None, width=10, type='Gaussian', method='standard', direction='central')

# for modifying the trajectory files, it saves the information
# about normal modes and store them as velocities 
# velcoities is just given as a name and any other name can be used


# Create a directory to store the modified trajectories
output_directory = 'ase_traj'
os.makedirs(output_directory, exist_ok=True)

# read the geometry.in file
geometry_atoms = read('geometry_opt.in')
num_atoms = len(geometry_atoms)

def read_velocities_from_xyz(filename, mode_number, num_atoms):
    velocities = []

    with open(filename, 'r') as file:
        while True:
            line = file.readline()
            if not line:
                break

            # Check if this line contains the mode number
            if f'Mode #{mode_number}' in line:
                # Extract the frequency part from the line
                parts = line.split(',')
                if len(parts) > 1:
                    frequency_part = parts[1].strip()  # 'f = 0.0i cm^-1.' or 'f = 0.0 cm^-1.'
                    frequency_value = frequency_part.split('=')[1].strip().split()[0]  # '0.0i' or '0.0'

                    # Check if the frequency contains 'i' indicating an imaginary part
                    if 'i' not in frequency_value:
                        # Read the next num_atoms lines containing velocities for the specific mode
                        for _ in range(num_atoms):
                            data = file.readline().split()
                            velocities.append([float(data[4]), float(data[5]), float(data[6])])
                break

    return np.array(velocities)

for mode_number in range(6, 3 * num_atoms):
    traj_path = f'vib.{mode_number}.traj'

    # Check if the trajectory file exists
    if os.path.exists(traj_path):
        # Load the trajectory for the current mode
        traj_mode = Trajectory(traj_path)

        # Load velocities from vib.xyz for the current mode
        new_velocities_mode = read_velocities_from_xyz('vib.xyz', mode_number, num_atoms)

        # Only proceed if new_velocities_mode is not empty (i.e., frequency was real)
        if new_velocities_mode.size > 0:
            # Create an empty trajectory for writing in the new directory
            modified_traj_path = os.path.join(output_directory, f'modified_traj_mode{mode_number}.traj')
            modified_traj_mode = Trajectory(modified_traj_path, 'w')

            # Iterate through frames, set velocities, and write to the new trajectory
            for atoms in traj_mode:
                atoms.set_velocities(new_velocities_mode)
                modified_traj_mode.write(atoms)

            # Close the new trajectory for the current mode
            modified_traj_mode.close()
        else:
            print(f"Skipping mode {mode_number} due to imaginary frequency.")
    else:
        print(f"File {traj_path} does not exist, skipping mode {mode_number}")


# Convert geometry.in to geometry.xyz using ASE
geometry_atom = read('geometry_opt.in')
write("geometry.xyz", geometry_atom, format="xyz")

# Copy geometry.xyz to the ase_traj folder
shutil.copy("geometry.xyz", output_directory)

# Copy vib.xyz to the ase_traj folder
shutil.copy('vib.xyz', output_directory)
