#!/bin/bash

set -e  # Exit on error

# Default Mode
MODE="TI"  # Options: TI, JE, BOTH
FF="MMFF"   # Options: MMFF, UFF
K=10.0     # Default JE force constant
SURF_NAME="none"
HARD_ATOMS=""
SOFT_ATOMS=""
HARD_DIST=""
SOFT_DIST=""
NSYS=200
NLAMBDA=200
NMDSTEPS=20000000
NEQSTEPS=1000
DT=0.05
TDAMP=150
TEMP=300
XYZ_NAME="../../cpp/common_resources/polymers/gui_builder/output/generated_system.xyz"
CONSTRAINTS="constraints_combined_systems.txt"
OUT_BASE="generated_system"

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --mode) MODE="$2"; shift ;;
        --ff)   FF="$2"; shift ;;
        --k)    K="$2"; shift ;;
        --surf_name|--surface|--surf) SURF_NAME="$2"; shift ;;
        --nSys) NSYS="$2"; shift ;;
        --nLambda) NLAMBDA="$2"; shift ;;
        --nMDsteps) NMDSTEPS="$2"; shift ;;
        --nEQsteps) NEQSTEPS="$2"; shift ;;
        --dt) DT="$2"; shift ;;
        --t_damp) TDAMP="$2"; shift ;;
        --temperature|-T) TEMP="$2"; shift ;;
        --constraints) CONSTRAINTS="$2"; shift ;;
        --xyz_name|--xyz) XYZ_NAME="$2"; shift ;;
        --hard_atoms) HARD_ATOMS="--hard_atoms";;
        --soft_atoms) SOFT_ATOMS="--soft_atoms";;
        --hard_dist) HARD_DIST="--hard_dist";;
        --soft_dist) SOFT_DIST="--soft_dist";;
        *) ;; # Ignore unknown args
    esac
    shift
done

if [[ -z "$HARD_ATOMS" && -z "$SOFT_ATOMS" && -z "$HARD_DIST" && -z "$SOFT_DIST" ]]; then
    SOFT_ATOMS="--soft_atoms"
fi

# Ensure we are in the script directory
cd "$(dirname "$0")"
ln -sfn ../../cpp/common_resources data
ln -sfn ../../cpp/common_resources common_resources
wd=`pwd`

if [ ! -f "$XYZ_NAME" ]; then
    echo "ERROR: XYZ file not found: $XYZ_NAME"
    exit 1
fi
if [ ! -f "$CONSTRAINTS" ]; then
    echo "ERROR: Constraints file not found: $CONSTRAINTS"
    exit 1
fi

# Build the library
echo "Step 1: Building libMMFFmulti_lib.so..."
cd ../../cpp/Build/libs_OCL/
rm -f libMMFFmulti_lib.so
make MMFFmulti_lib
if [ $? -ne 0 ]; then
    echo "ERROR: Build failed!"
    exit 1
fi
cd $wd
echo "Build successful!"
echo ""

# Run the calculation
echo "Step 2: Running Free Energy Calculation for combined systems (Force field: $FF, Mode: $MODE)..."
echo "----------------------------------------"
echo "XYZ: $XYZ_NAME"
echo "Constraints: $CONSTRAINTS"
echo "Surface: $SURF_NAME"
echo "Constraint mode: ${HARD_ATOMS}${SOFT_ATOMS}${HARD_DIST}${SOFT_DIST}"
python3 run_ES.py \
    --mode $MODE \
    --ff $FF \
    --nSys $NSYS \
    --xyz_name "$XYZ_NAME" \
    --surf_name "$SURF_NAME" \
    --nLambda $NLAMBDA \
    --nMDsteps $NMDSTEPS \
    --nEQsteps $NEQSTEPS \
    --Fconv 1e-6 \
    --constraints "$CONSTRAINTS" \
    --K $K \
    --dt $DT \
    -T $TEMP \
    --t_damp $TDAMP \
    $HARD_ATOMS \
    $SOFT_ATOMS \
    $HARD_DIST \
    $SOFT_DIST

if [ $? -ne 0 ]; then
    echo "ERROR: Calculation failed!"
    exit 1
fi
if [ ! -s "${OUT_BASE}_free_energy.dat" ] || ! grep -q '^[[:space:]]*[-+0-9.]' "${OUT_BASE}_free_energy.dat"; then
    echo "ERROR: Calculation did not write data rows to ${OUT_BASE}_free_energy.dat"
    exit 1
fi
echo ""

# Plot the results
echo "Step 3: Plotting results..."
echo "----------------------------------------"
python3 plot_F_interactive.py --input ${OUT_BASE}_free_energy.dat
if [ $? -ne 0 ]; then
    echo "ERROR: Plotting failed!"
    exit 1
fi
echo ""

echo "=========================================="
echo "  combined_systems completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - ${OUT_BASE}_free_energy.dat (raw data)"
echo "  - ${OUT_BASE}_free_energy_F_interactive.html (interactive plot)"
echo ""
echo "To view the interactive plot, open ${OUT_BASE}_free_energy_F_interactive.html in a web browser"
echo ""
