#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources

set -e  # Exit on error

# Default Mode
MODE="BOTH"  # Options: TI, JE, BOTH
JE_K=20.0     # Default JE force constant

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --mode) MODE="$2"; shift ;;
        --k)    JE_K="$2"; shift ;;
        *) ;; # Ignore unknown args
    esac
    shift
done

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

# Run the calculation
echo "Step 2: Running Free Energy Calculation for DA (Mode: $MODE)..."
echo "----------------------------------------"
# python3 run_ES.py \
#     --mode $MODE \
#     --nSys 80 \
#     --xyz_name "../../cpp/common_resources/xyz/DA.xyz" \
#     --nLambda 80 \
#     --nMDsteps 4000000 \
#     --nEQsteps 50000 \
#     --Fconv 1e-6 \
#     --constraints "constraints_DA.txt" \
#     --JEforceconst $JE_K \
#     --dt 0.05 \
#     -T 300 \
#     --t_damp 150
python3 run_ES.py \
    --mode $MODE \
    --nSys 100 \
    --xyz_name "../../cpp/common_resources/xyz/DA.xyz" \
    --nLambda 100000 \
    --nMDsteps 40000000 \
    --nEQsteps 50000 \
    --Fconv 1e-6 \
    --constraints "constraints_DA.txt" \
    --JEforceconst $JE_K \
    --dt 0.05 \
    -T 300 \
    --t_damp 150

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
