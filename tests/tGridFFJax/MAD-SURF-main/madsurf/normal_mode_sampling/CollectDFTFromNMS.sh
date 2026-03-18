#!/bin/bash

molecule_dir=$PWD  # Replace with the path to the molecule directory
counter=1
total=$(find "$molecule_dir"/* -maxdepth 1 -type d | wc -l)
collected_file="$molecule_dir/collected_dft_results.extxyz"

for molecule_subdir in "$molecule_dir"/*/; do
    echo "Processing directory $counter of $total: $molecule_subdir"

    # Path to the DFT results file
    results_file="$molecule_subdir/DFT_cal/selected_indices/dft_results.xyz"

    # Check if the file exists
    if [[ -f "$results_file" ]]; then
        echo "Extracting from $results_file"
        cat "$results_file" >> "$collected_file"
        echo "Collected results from $results_file"
    else
        echo "No DFT_results.extxyz file found in $molecule_subdir"
    fi

    counter=$((counter + 1))
done

echo "All NMS results collected."

