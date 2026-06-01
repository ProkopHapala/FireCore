#!/bin/bash

ln -sfn ../../cpp/common_resources data
ln -sfn ../../cpp/common_resources common_resources 

rm -f trajectory.xyz
set -e  # Exit on error

# Default comparison setup aligned with examples/tFreeEnergy/{CPU,GPU-debug}
MODE="TI"   # Options: TI, JE, BOTH
FF="MMFF"   # Options: MMFF, UFF
K=10.0
SURF_NAME="none"
HARD_ATOMS=""
SOFT_ATOMS=""
HARD_DIST=""
SOFT_DIST=""
NSYS=100
NLAMBDA=100
NMDSTEPS=10000000
NEQSTEPS=5000
DT=0.05
TDAMP=100
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
    SOFT_DIST="--soft_dist"
fi

N=30

OUT_DIR="results/entropic_spring_${N}"
mkdir -p "$OUT_DIR"

# Ensure we are in the script directory
cd "$(dirname "$0")"
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

# Run calculation
echo "Step 2: Running Free Energy Calculation (Force field: $FF, Mode: $MODE)..."
echo "----------------------------------------"
echo "Surface: $SURF_NAME"
python3 scripts/run/run_ES.py \
    --mode $MODE \
    --ff $FF \
    --nSys $NSYS \
    --xyz_name "data/entropic_spring_$N.xyz" \
    --surf_name "$SURF_NAME" \
    --nLambda $NLAMBDA \
    --nMDsteps $NMDSTEPS \
    --nEQsteps $NEQSTEPS \
    --Fconv 1e-6 \
    --constraints "configs/constraints_ES.txt" \
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
mv entropic_spring_${N}_free_energy.dat "$OUT_DIR/"
[ -f jarzynski_work.dat ] && mv jarzynski_work.dat "$OUT_DIR/"
echo ""

# Plot the results
echo "Step 3: Plotting results..."
echo "----------------------------------------"
python3 scripts/analysis/plot_F_interactive.py --input "$OUT_DIR/entropic_spring_${N}_free_energy.dat"
if [ $? -ne 0 ]; then
    echo "ERROR: Plotting failed!"
    exit 1
fi
echo ""

echo "=========================================="
echo "  Completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - $OUT_DIR/entropic_spring_${N}_free_energy.dat (raw data)"
echo "  - $OUT_DIR/entropic_spring_${N}_free_energy_F_interactive.html (interactive plot)"
echo ""
echo "To view the interactive plot, open $OUT_DIR/entropic_spring_${N}_free_energy_F_interactive.html in a web browser"
echo ""
