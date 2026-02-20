#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

rm trajectory.xyz
set -e  # Exit on error

# Default Mode
MODE="BOTH"  # Options: TI, JE, BOTH
JE_K=5.0     # Default JE force constant

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --mode) MODE="$2"; shift ;;
        --k)    JE_K="$2"; shift ;;
        *) ;; # Ignore unknown args
    esac
    shift
done

N=30

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
echo "Step 2: Running Free Energy Calculation (Mode: $MODE)..."
echo "----------------------------------------"
python3 run_ES.py \
    --mode $MODE \
    --nSys 10 \
    --xyz_name "../tMMFF/data/entropic_spring_$N.xyz" \
    --system_name "entropic_spring_$N" \
    --nLambda 100000 \
    --nMDsteps 10000000 \
    --nEQsteps 50000 \
    --Fconv 1e-6 \
    --constraints "constraints_ES.txt" \
    --JEforceconst $JE_K \
    --nPerVFs 1

if [ $? -ne 0 ]; then
    echo "ERROR: Calculation failed!"
    exit 1
fi
echo ""

# Plot the results
echo "Step 3: Plotting results..."
echo "----------------------------------------"
python3 plot_F_interactive.py --input entropic_spring_${N}_free_energy.dat
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
echo "  - entropic_spring_${N}_TI.dat (raw data)"
echo "  - entropic_spring_${N}_TI_interactive.html (interactive plot)"
echo ""
echo "To view the interactive plot, open entropic_spring_${N}_TI_interactive.html in a web browser"
echo ""
