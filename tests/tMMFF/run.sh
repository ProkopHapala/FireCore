#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
pwd
rm libMMFF_lib.so
make MMFF_lib
rm   libLattice2D_lib.so
make Lattice2D_lib
cd $wd

cd ../../cpp/Build/libs_SDL
rm libMolGUIlib.so
make MolGUIlib
cd $wd


# ---- Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS
export PYOPENCL_CTX=0

#rm *.bin

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
LD_PRELOAD=$LD_PRELOAD  $(g++ -print-file-name=libfftw3.so)
echo   $LD_PRELOAD
export LD_PRELOAD
# --- ignore memory leaks in ASAM
#export LSAN_OPTIONS=detect_leaks=0

# python3 run.py 
# # For 1D scan in a custom direction
# python3 generate_scans.py --scan-types total --scan-mode 1d --scan-dir1 "0,0,1" --output-dir PTCDA_data_trial_1d --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --scan-origin "0,0,0" --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types morse --scan-mode 1d --scan-dir1 "0,0,1" --output-dir PTCDA_data_trial_1d --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --scan-origin "0,0,0" --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types coulomb --compare --scan-mode 1d --scan-dir1 "0,0,1" --output-dir PTCDA_data_trial_1d --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --scan-origin "0,0,0" --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9

# # For 2D scan in a custom plane
# python3 generate_scans.py --scan-types total --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial_2d --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.1 --nscan-2 41 --span-min-2 0 --span-max-2 4.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types morse --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial_2d --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.1 --nscan-2 41 --span-min-2 0 --span-max-2 4.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types coulomb --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial_2d --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.1 --nscan-2 41 --span-min-2 0 --span-max-2 4.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9


# For 1D relaxed scan
# python3 generate_scans.py \
#   --scan-types total \
#   --scan-mode 1d \
#   --scan-dir1 "0.0,0.0,1.0" \
#   --scan-origin "0.0,0.0,0.0" \
#   --output-dir PTCDA_data_trial_1d_relax_z \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan \
#   --nscan 201 \
#   --span-min 1.3 \
#   --span-max 21.4 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26\

# python3 generate_scans.py \
#   --scan-types morse \
#   --scan-mode 1d \
#   --scan-dir1 "0.0,0.0,1.0" \
#   --scan-origin "0.0,0.0,0.0" \
#   --output-dir PTCDA_data_trial_1d_relax_z \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan \
#   --nscan 201 \
#   --span-min 1.3 \
#   --span-max 21.4 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26\

# python3 generate_scans.py \
#   --scan-types coulomb \
#   --scan-mode 1d \
#   --scan-dir1 "0.0,0.0,1.0" \
#   --scan-origin "0.0,0.0,0.0" \
#   --output-dir PTCDA_data_trial_1d_relax_z \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan \
#   --nscan 201 \
#   --span-min 1.3 \
#   --span-max 21.4 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 10000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26 \
#   --compare

## For 1D relaxed scan in angle
# python3 generate_scans.py \
#   --scan-types total \
#   --scan-mode 1d \
#   --scan-dir1 "0.894427190999916,0.447213595499958,0.0" \
#   --scan-origin "0.0,0.0,2.8" \
#   --output-dir trial_20_1d_relax_angle \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
#   --nscan (2*448) \
#   --span-min 0.0 \
#   --span-max (2*44.778640450004204) \
#   --molecule data/xyz/PTCDA_20x20_26 \
#   --substrate data/xyz/NaCl_perfect_20x20 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26
 

# python3 generate_scans.py \
#   --scan-types morse \
#   --scan-mode 1d \
#   --scan-dir1 "2.0,1.0,0.0" \
#   --scan-origin "0.0,0.0,2.8" \
#   --output-dir PTCDA_data_trial_1d_relax_line  \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
#   --nscan 357 \
#   --span-min 0.0 \
#   --span-max 35.7 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26\

# python3 generate_scans.py \
#   --scan-types coulomb \
#   --scan-mode 1d \
#   --scan-dir1 "2.0,1.0,0.0" \
#   --scan-origin "0.0,0.0,2.8" \
#   --output-dir PTCDA_data_trial_1d_relax_line  \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
#   --nscan 357 \
#   --span-min 0.0 \
#   --span-max 35.7 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 10000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26 \
#   --compare
# #  

# python3 generate_scans.py \
#   --scan-types total \
#   --scan-mode 1d \
#   --scan-dir1 "2.0,1.0,0.0" \
#   --scan-origin "0.0,0.0,2.8" \
#   --output-dir trial_1d_relax_line \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
#   --nscan 357 \
#   --span-min 0.0 \
#   --span-max 35.7 \
#   --molecule data/xyz/PTCDA_12x12_26 \
#   --substrate data/xyz/NaCl_45_defect_aligned_12x12 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26




### Multiple scans with loops for all configurations

SUB_SIZES=("20x20" "16x16" "12x12" "8x8")   #

for SUB_SIZE in "${SUB_SIZES[@]}"; do

SIZE_X=${SUB_SIZE%x*}         # first number (12, 20, 32 …)
LATTICE_UNIT=4               # Å per repeat of NaCl (4 × SIZE_X)
BASE_SUBS=("NaCl_perfect" \
           "NaCl_45_defect_aligned" "NaCl_45_defect_perpendicular" \
           "NaCl_210_defect_aligned" "NaCl_210_defect_perpendicular")

SUBSTRATES=()
for b in "${BASE_SUBS[@]}"; do
  SUBSTRATES+=("${b}_${SUB_SIZE}")
done

CONS_ATOMS=("26" "24")
DIRECTIONS=("1.0,1.0,0.0" "2.0,1.0,0.0")

# Step size for scan (in Angstroms)
STEP_SIZE=0.1

function get_lattice_x() {
  local size=${1##*_}        # “…_20x20” -> “20x20”
  local nx=${size%x*}        # “20”
  echo "$(bc <<< "$nx * $LATTICE_UNIT")"
}

# Base parameters that stay the same (without span and nscan which will be calculated)
BASE_PARAMS="--scan-types total \
  --scan-mode 1d \
  --scan-origin \"0.0,0.0,2.8\" \
  --relaxed \
  --niter-max 100000 \
  --dt 0.02 \
  --fconv 1e-3"

# Loop through constraint atoms and directions first
for cons_atom in "${CONS_ATOMS[@]}"; do
  # Set molecule based on cons_atom
  molecule="data/xyz/PTCDA_${SUB_SIZE}_${cons_atom}"
  
  for direction in "${DIRECTIONS[@]}"; do
    # Determine which defect substrates to use based on direction
    case "$direction" in
      "1.0,1.0,0.0") defect_tag="45"  ;;   # 110
      "2.0,1.0,0.0") defect_tag="210" ;;   # 210
      *) echo "Unknown direction $direction"; exit 1 ;;
    esac

    direction_substrates=(
      "NaCl_perfect_${SUB_SIZE}"
      "NaCl_${defect_tag}_defect_aligned_${SUB_SIZE}"
      "NaCl_${defect_tag}_defect_perpendicular_${SUB_SIZE}"
    )

    # Loop through the appropriate substrates for this direction
    for substrate in "${direction_substrates[@]}"; do
      # Skip if this substrate is not in the main SUBSTRATES array
      if [[ ! " ${SUBSTRATES[@]} " =~ " ${substrate} " ]]; then
        continue
      fi
      # Determine substrate type for directory structure
      if [[ "$substrate" == *"perfect"* ]]; then
        substrate_type="perfect"
      elif [[ "$substrate" == *"aligned"* ]]; then
        substrate_type="defect_aligned"
      elif [[ "$substrate" == *"perpendicular"* ]]; then
        substrate_type="defect_perpendicular"
      fi
      # Format direction for directory naming (replace commas and dots)
      dir_name=$(echo $direction | tr ',' '_')
      
      # Parse the direction vector components
      IFS=',' read -r x_comp y_comp z_comp <<< "$direction"
      
      # Calculate the angle between direction vector and x-axis (1,0,0) # First, normalize the direction vector
      magnitude=$(echo "sqrt($x_comp^2 + $y_comp^2 + $z_comp^2)" | bc -l)
      x_norm=$(echo "$x_comp / $magnitude" | bc -l)
      
      # Calculate cosine of the angle (dot product with (1,0,0) / magnitudes = x_norm)
      cosine_angle=$x_norm
      
      # Get the lattice constant for this substrate
      lattice_x=$(get_lattice_x "$substrate")
      
      # Calculate span-max: lattice_x / cos(angle) # Handle the case where cosine is very small to avoid division by zero
      if (( $(echo "$cosine_angle < 0.01" | bc -l) )); then
        cosine_angle=0.01
      fi
      
      span_max=$(echo "$lattice_x / $cosine_angle" | bc -l)
      # Round to 1 decimal place
      span_max=$(printf "%.1f" $span_max)
      
      # Calculate nscan to maintain step size of 0.1
      nscan=$(echo "$span_max / $STEP_SIZE" | bc)
      
      # Create output directory structure
      output_dir="relax_${substrate_type}_line/dir_${dir_name}/cons_${cons_atom}"
      
      echo "\n\n==== Running scan for: ===="
      echo "Substrate: $substrate"
      echo "Molecule: $molecule"
      echo "Direction: $direction"
      echo "Constraint atom: $cons_atom"
      echo "Output directory: $output_dir"
      echo "Calculated span_max: $span_max"
      echo "Calculated nscan: $nscan"
      
      log_file="timings.dat"
      # Run the command
      SECONDS=0
      python3 generate_scans.py \
        $BASE_PARAMS \
        --scan-dir1 "$direction" \
        --output-dir "$output_dir" \
        --lammps-1d-dir "/home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan" \
        --molecule "$molecule" \
        --substrate "data/xyz/$substrate" \
        --cons-atom "$cons_atom" \
        --span-min 0.0 \
        --span-max "$span_max" \
        --nscan "$nscan" \
        --scan-origin "0.0,0.0,2.8"
      elapsed=$SECONDS
      echo "⏱  run time: ${elapsed}s"
      # ----- append to log -----
      printf "%(%F %T)T,%s,%s,%s,%s,%s,%s,%s\n" -1 "$SUB_SIZE" "$substrate" "$direction" "$cons_atom" "$span_max" "$nscan" "$elapsed" >> "$log_file"
      # Optional: add a small delay between runs
      sleep 10
    done
  done
done
done




#python3 run_surf_lattice.py
#python3 run_propandiol.py
#python3 run_sample.py
#python3 run_Hbonds.py
#python3 run_sample_func.py
#python3 run_sample_Bsplines.py
#python3 run_sample_Hermite.py
#python3 run_test_ewald.py
#python3 run_test_GridFF.py
#python3 run_test_GridFF_ocl.py
#python3 run_test_GridFF_ocl_new.py
#python3 run_test_Multipole.py

# python3 run_sample_surf.py

#python3 run_sample_tricubic.py
#python3 run_sample_nonBond.py
#python3 run_lat_scan.py

#python3 run_collision_damp.py
#python3 run_collision_damp_scan.py
#python3 run_test_clear.py

#python3 run_opt_poly.py BB.HNH-hh.NHO-hp
#python3 run_opt_poly.py BB.HNH-hp.OHO-h_1,BB.HNH-hh.NHO-hp

#python3 run_test_ewald.py
