import ase.io.aims
from ase import io
from ase.io.aims import read_aims_output
from ase.io.extxyz import write_extxyz
import os
import glob

# Set the folder containing the relaxation files
relaxation_folder = os.path.join(os.getcwd(), 'relaxation_files')

# Output file name for the extracted geometries
output_opt = "collected_opt_data.extxyz"
output_rel = "collected_relaxation_data.extxyz"

# Find all output files in the relaxation_files folder
output_files = sorted(glob.glob(os.path.join(relaxation_folder, 'aims_*.out')))

# Loop through each file in the folder
for output_path in output_files:
    # Print the file being processed
    print(f"Extracting geometries from {output_path}")
    
    # Determine the number of geometries
    num_geometries = 0
    with open(output_path, 'r') as f:
        for line in f:
            if "Geometry optimization" in line:
                num_geometries += 1
    
    # Save only the final optimized geometry in extxyz format
    output = read_aims_output(output_path, index=num_geometries-1)
    write_extxyz(output_opt, images=output, append=True)

    # Loop through geometries in output and save all of them in extxyz format
    for i in range (0, num_geometries):
        output = read_aims_output(output_path, index=i)
        write_extxyz(output_rel, images=output, append=True)


