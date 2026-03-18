#!/bin/bash
#SBATCH --account=project_2008059
#SBATCH -p test # test - for testing 1h ; medium - up to 20 nodes/36 hours ; large - 20-200 nodes/36 hours #
#SBATCH --time=1:00:00      # dd-hh:mm:ss
#SBATCH -J process_aims_aiss # name 
#SBATCH -o sbatch-%j.out # where the outputs & errors are written
#SBATCH -N 1                  # N nodes to run (N x 64 = n); max 192 ; max debug 48
#SBATCH --ntasks-per-node=1   # n processes to run (N x 64 = n); max 192 ; max debug 48

# Check modules.txt for updating following lines: # 
#module purge

# Base directory containing the calc_* directories
BASE_DIR="/scratch/project_2008059/mgonzalez/aims_calc_from_aiss"

# Loop through all directories with the prefix calc_*
for dir in $BASE_DIR/calc_*; do
    # Check if the directory exists and contains acq.out
    if [[ -d "$dir" && -f "$dir/acq.out" ]]; then
        echo "Processing directory: $dir"
        # Enter the directory
        cd "$dir" || { echo "Failed to enter $dir"; continue; }
        # Create a symbolic link to the Python script
        ln -sf "$BASE_DIR/aims_to_extxyz_dataset.py" aims_to_extxyz_dataset.py
        # Run the Python script to generate dataset.extxyz
        python aims_to_extxyz_dataset.py acq.out dataset.extxyz
        # Return to the original directory
        cd "$BASE_DIR" || { echo "Failed to return to $BASE_DIR"; exit 1; }
    else
        echo "Skipping $dir: Either not a directory or acq.out not found."
    fi
done
