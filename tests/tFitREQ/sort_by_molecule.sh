#!/bin/bash
# Script to move .xyz files into molecule directories based on filename matching

# Use current directory, or accept optional argument
DIR="${1:-.}"
cd "$DIR" || exit 1

# Get list of molecule directories (skip non-directories)
for mol_dir in */; do
    # Remove trailing slash to get molecule name
    mol_name="${mol_dir%/}"
    
    echo "Processing molecule: $mol_name"
    
    # Find and move .xyz files containing the molecule name
    # Use find to avoid issues with glob expansion
    find "$DIR" -maxdepth 1 -type f -name "*${mol_name}*.xyz" | while read -r file; do
        # Get just the filename
        filename=$(basename "$file")
        echo "  Moving $filename -> $mol_dir"
        mv "$file" "$mol_dir/"
    done
done

echo "Done sorting files by molecule."
