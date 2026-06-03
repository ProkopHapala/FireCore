#!/bin/bash

set -e  # Exit on error

# Get working directory
wd=`pwd`

# Default variables
ANALYZE_PAIRS=""
HBOND_CUT=2.5
LOG_TEMPERATURE=""

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --analyze_pairs|--hbond_pairs) ANALYZE_PAIRS="$2"; shift ;;
        --hbond_cut) HBOND_CUT="$2"; shift ;;
        --log_temperature|--log_temp) LOG_TEMPERATURE="--log_temperature";;
        *) ;; # Ignore unknown args
    esac
    shift
done

OUT_DIR="results/nHexadecan"
mkdir -p "$OUT_DIR"

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
echo "Step 2: Running Thermodynamic Integration for nHexadecane..."
echo "----------------------------------------"
python3 scripts/run/run_ES.py \
    --nSys 100 \
    --xyz_name "../../cpp/common_resources/xyz/nHexadecan.xyz" \
    --nLambda 100 \
    --nMDsteps 2000000 \
    --nEQsteps 5000 \
    --Fconv 1e-6 \
    --constraints "configs/constraints_nHex.txt" \
    --analyze_pairs "$ANALYZE_PAIRS" \
    --hbond_cut $HBOND_CUT \
    $LOG_TEMPERATURE

if [ $? -ne 0 ]; then
    echo "ERROR: TI calculation failed!"
    exit 1
fi
mv nHexadecan_free_energy.dat "$OUT_DIR/"
[ -f jarzynski_work.dat ] && mv jarzynski_work.dat "$OUT_DIR/"
echo ""

# Plot the results
echo "Step 3: Plotting results..."
echo "----------------------------------------"
python3 scripts/analysis/plot_F_interactive.py --input "$OUT_DIR/nHexadecan_free_energy.dat"
if [ $? -ne 0 ]; then
    echo "ERROR: Plotting failed!"
    exit 1
fi
echo ""

echo "=========================================="
echo "  nHexadecane TI Completed successfully!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  - $OUT_DIR/nHexadecan_free_energy.dat (raw data)"
echo "  - $OUT_DIR/nHexadecan_free_energy_F_interactive.html (interactive plot)"
echo ""
echo "To view the interactive plot, open $OUT_DIR/nHexadecan_free_energy_F_interactive.html in a web browser"
echo ""
