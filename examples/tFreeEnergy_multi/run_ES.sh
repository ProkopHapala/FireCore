#!/bin/bash

set -e  # Exit on error

N=30

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

# Run the TI calculation
echo "Step 2: Running Thermodynamic Integration..."
echo "----------------------------------------"
python3 run_ES.py \
    --nSys 100 \
    --xyz_name "../tMMFF/data/entropic_spring_$N.xyz" \
    --system_name "entropic_spring_$N" \
    --nLambda 100 \
    --nMDsteps 1000000 \
    --nEQsteps 5000 \
    --Fconv 1e-6 \
    --constraints "constraints_ES.txt" # Values for the distances of the two end atoms to pull

if [ $? -ne 0 ]; then
    echo "ERROR: TI calculation failed!"
    exit 1
fi
echo ""

# Plot the results
echo "Step 3: Plotting results..."
echo "----------------------------------------"
python3 plot_TI_interactive.py --input entropic_spring_${N}_TI.dat
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
