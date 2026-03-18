#!/bin/bash

molecule_dir=$PWD  # Set the base directory to the current working directory
slurm_script_dir="./general_folder"  # Directory containing the SLURM script
slurm_script="run_all_NMS.sh"  # Name of the SLURM script
counter=1
total=$(find "$molecule_dir"/* -maxdepth 1 -type d | wc -l)

# Check if the SLURM script exists in the specified folder
if [[ ! -f "$slurm_script_dir/$slurm_script" ]]; then
    echo "Error: SLURM script $slurm_script not found in $slurm_script_dir"
    exit 1
fi

for molecule_subdir in "$molecule_dir"/*/; do
    echo "Processing directory $counter of $total: $molecule_subdir"

    # Submit the SLURM script from the specified folder
    cd "$molecule_subdir" || { echo "Failed to access $molecule_subdir"; exit 1; }
    sbatch "$slurm_script_dir/$slurm_script"
    echo "Submitted $slurm_script from $slurm_script_dir for $molecule_subdir"

    # Return to the original directory
    cd "$molecule_dir" || { echo "Failed to return to $molecule_dir"; exit 1; }

    counter=$((counter + 1))
done

echo "Job submission process completed."

