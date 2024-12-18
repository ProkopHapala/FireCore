#!/bin/bash

path="./molecules"
molecules=("$path"/*.xyz)

for name in "${molecules[@]}"; do
    #name=$(basename "$mol")  # Extract filename without path
    echo "Processing: $name"
    python3 relax_molecule.py "$name" relaxed_mols  # Pass the molecule name to the Python script
done