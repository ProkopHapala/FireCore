#!/bin/bash

set -e  # Exit on error

# Default Mode
MODE="BOTH"  # Options: TI, JE, BOTH
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
# Get working directory
wd=`pwd`

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
echo "Step 2: Running Free Energy Calculation for DA (Force field: $FF, Mode: $MODE)..."
echo "----------------------------------------"
echo "Surface: $SURF_NAME"
python3 run_ES.py \
    --mode $MODE \
    --ff $FF \
    --nSys $NSYS \
    --xyz_name "../../cpp/common_resources/xyz/DA.xyz" \
    --surf_name "$SURF_NAME" \
    --nLambda $NLAMBDA \
    --nMDsteps $NMDSTEPS \
    --nEQsteps $NEQSTEPS \
    --Fconv 1e-6 \
    --constraints "constraints_DA.txt" \
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
echo ""

# Plot the results
echo "Step 3: Plotting results..."
echo "----------------------------------------"
python3 plot_F_interactive.py --input DA_free_energy.dat
if [ $? -ne 0 ]; then
    echo "ERROR: Plotting failed!"
    exit 1
fi
echo ""

echo "=========================================="
echo "  DA Completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - DA_free_energy.dat (raw data)"
echo "  - DA_free_energy_interactive.html (interactive plot)"
echo ""
echo "To view the interactive plot, open DA_free_energy_interactive.html in a web browser"
echo ""
