#!/bin/bash

ln -sf ../../cpp/common_resources data
ln -sf ../../cpp/common_resources common_resources 

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
LD_PRELOAD=$LD_PRELOAD  $(g++ -print-file-name=libfftw3.so)
echo   $LD_PRELOAD
export LD_PRELOAD

scan_types=("total" "morse" "coulomb")
for scan_type in "${scan_types[@]}"; do
SUB_SIZES=("20x20" "16x16" "12x12" "8x8")   #
#UB_SIZES=("12x12")   #

for SUB_SIZE in "${SUB_SIZES[@]}"; do

SIZE_X=${SUB_SIZE%x*}         # first number (12, 20, 32 …)
LATTICE_UNIT=4               # Å per repeat of NaCl (4 × SIZE_X)
BASE_SUBS=("NaCl_perfect")

SUBSTRATES=()
for b in "${BASE_SUBS[@]}"; do
  SUBSTRATES+=("${b}_${SUB_SIZE}")
done

CONS_ATOMS=("26")
DIRECTIONS=("1.0,1.0,0.0")

# Step size for scan (in Angstroms)
STEP_SIZE=0.1

function get_lattice_x() {
  local size=${1##*_}        # “…_20x20” -> “20x20”
  local nx=${size%x*}        # “20”
  echo "$(bc <<< "$nx * $LATTICE_UNIT")"
}

# Base parameters that stay the same (without span and nscan which will be calculated)
BASE_PARAMS="--scan-types $scan_type \
  --scan-mode 1d \
  --scan-origin \"0.0,0.0,2.8\" \
  --relaxed \
  --niter-max 5000 \
  --dt 0.02 \
  --fconv 1e-30"

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
      output_dir="perfor_${substrate_type}_line/dir_${dir_name}/cons_${cons_atom}"
      
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
        --lammps-1d-dir "/storage/praha1/home/indranil/work/4-relaxed_linescan" \
        --molecule "$molecule" \
        --substrate "data/xyz/$substrate" \
        --cons-atom "$cons_atom" \
        --span-min 0.0 \
        --span-max "$span_max" \
        --nscan "$nscan" \
        --scan-origin "0.0,0.0,2.8"
      elapsed=$SECONDS
      echo "run time: ${elapsed}s"
      # ----- append to log -----
      printf "%(%F %T)T,%s,%s,%s,%s,Time= %s\n" -1 "$SUB_SIZE" "$substrate" "$direction" "$cons_atom"   "$elapsed" >> "$log_file"  
      # Optional: add a small delay between runs
      # sleep 10
    done
  done
done
done
done

