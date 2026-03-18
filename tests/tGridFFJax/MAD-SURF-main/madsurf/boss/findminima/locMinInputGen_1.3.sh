#!/bin/bash -l
#set -e
# Example of use: ./locMinInputGen_1.3.sh adsorbate_opt.in substrate_opt.in minimum_predictions 
# Remove duplicates and equivalent points:
python /scratch/project_2008059/mgonzalez/PUBLIC/scripts/python/uniqueSymmetry_fcc111_1.2.py "$3".dat
awk '{print $(NF), $0}' "$3"_unique.dat | sort -n | cut -d' ' -f2- > "$3"_unique_sorted.dat
# Generates DFT input files (for fhi-aims) from BOSS output (pp_local_min). Remember to filter out duplicates!
# Reads in local min variables from datafile line by line (command line input two "$2"), creates the corresponding aims geometry.in using the 
# adsorbate.xyz file (command line input one "$1") with the variables. 
counter=0
while read -r x y z Rx Ry Rz _ _; do
  echo $x $y $z $Rx $Ry $Rz > variables.dat
  ((counter++))
  output_dir="locmin_${counter}"
  input_file="geometry.in"
  mkdir -p "$output_dir"
  python /scratch/project_2008059/mgonzalez/PUBLIC/scripts/python/interface_builder_boss_1.4.py -ia "$1" -is "$2" -t "$x" "$y" -z "$z" -r "$Rx" "
$Ry" "$Rz" -ff aims -ss 4 4 3 -ca 2 -lc 3.632  
  mv "$input_file" "$output_dir"
  mv "variables.dat" "$output_dir"
  cp "control.in" "$output_dir"
done < "$3"_unique_sorted.dat
