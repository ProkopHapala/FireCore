#!/bin/bash

# Configuration
CONFIG="configs/default_free_energy_config.json"
ANALYZE_PAIRS=""
HBOND_CUT=""
LOG_TEMPERATURE=""

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --config) CONFIG="$2"; shift ;;
        --analyze_pairs|--hbond_pairs) ANALYZE_PAIRS="$2"; shift ;;
        --hbond_cut) HBOND_CUT="$2"; shift ;;
        --log_temperature|--log_temp) LOG_TEMPERATURE="--log_temperature";;
        *)
            # Backward compatibility: first non-option is treated as CONFIG
            if [[ -z "$1" || "$1" == --* ]]; then
                echo "Warning: Unknown argument $1"
            else
                CONFIG="$1"
            fi
            ;;
    esac
    shift
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure required symlinks exist
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/data"              2>/dev/null || true
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/common_resources"  2>/dev/null || true

# Step 1: Build the library
echo "Step 1: Building libMMFFmulti_lib.so..."
cd ../../cpp/Build/libs_OCL/
make MMFFmulti_lib
if [ $? -ne 0 ]; then
    echo "ERROR: Build failed!"
    exit 1
fi
cd "$SCRIPT_DIR"
echo "Build successful!"
echo ""

# Step 2: Run the sweep
echo "Step 2: Running Free Energy Calculation Sweep with config: $CONFIG"
echo "----------------------------------------"
RUNNER_ARGS=""
if [ -n "$ANALYZE_PAIRS" ]; then
    RUNNER_ARGS="$RUNNER_ARGS --analyze_pairs $ANALYZE_PAIRS"
fi
if [ -n "$HBOND_CUT" ]; then
    RUNNER_ARGS="$RUNNER_ARGS --hbond_cut $HBOND_CUT"
fi
if [ -n "$LOG_TEMPERATURE" ]; then
    RUNNER_ARGS="$RUNNER_ARGS --log_temperature"
fi

python3 scripts/run/free_energy_runner.py --config "$CONFIG" $RUNNER_ARGS
if [ $? -ne 0 ]; then
    echo "ERROR: Sweep failed!"
    exit 1
fi
echo ""

# Step 3: Run temperature sweep plotting
OUTPUT_ROOT=$(python3 -c "import json; print(json.load(open('$CONFIG')).get('output_root', 'results_sweep'))" 2>/dev/null || echo "results_sweep")
echo "Step 3: Plotting temperature dependence from $OUTPUT_ROOT..."
python3 scripts/analysis/plot_temperature_sweep.py --results_dir "$OUTPUT_ROOT"
echo ""

echo "=========================================="
echo "  Sweep completed successfully!"
echo "=========================================="
echo ""
