#!/bin/bash

# Check if directories are passed as arguments
if [ $# -eq 0 ]; then
  echo "No directories supplied. Usage: ./script.sh locmin_*"
  exit 1
fi

# Output file
output_file="energies.out"

# Clear the output file if it exists
> $output_file

# Read substrate and adsorbate energies from their respective files
if [ -f "substrate_opt_energy.txt" ]; then
  substrate_energy=$(grep -i "TOTAL ENERGY (eV):" substrate_opt_energy.txt | awk -F ': *' '{print $2}')
else
  echo "File substrate_opt_energy.txt does not exist. Exiting."
  exit 1
fi
if [ -f "adsorbate_opt_energy.txt" ]; then
  adsorbate_energy=$(grep -i "TOTAL ENERGY (eV):" adsorbate_opt_energy.txt | awk -F ': *' '{print $2}')
else
  echo "File adsorbate_opt_energy.txt does not exist. Exiting."
  exit 1
fi
# Convert energies to float to ensure bc handles them correctly
substrate_energy=$(printf "%.15f" "$substrate_energy")
adsorbate_energy=$(printf "%.15f" "$adsorbate_energy")
# Loop through each directory provided as arguments
for dir in "$@"
do
  # Check if acq.out exists in the directory
  if [ -f "$dir/acq.out" ]; then
    # Extract energy from acq.out
    total_energy=$(grep -i "| Total energy of the DFT / Hartree-Fock s.c.f. calculation" "$dir/acq.out" | awk -F ': *' '{print $2}' | awk '{print
 $1}')
    total_energy=$(printf "%.15f" "$total_energy")
    
    # Calculate adsorption energy using bc for high precision
    adsorption_energy=$(echo "scale=15; $total_energy - $substrate_energy - $adsorbate_energy" | bc -l)
    
    # Append directory name and adsorption energy to energies.out
    echo "$dir    $adsorption_energy eV" >> $output_file
  else
    echo "File $dir/acq.out does not exist, skipping directory."
  fi
done
echo "Adsorption energies have been written to $output_file"