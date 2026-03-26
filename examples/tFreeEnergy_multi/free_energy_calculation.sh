#!/bin/bash

# Configuration
CONFIG=${1:-default_free_energy_config.json}
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
python3 run_general_sweep.py --config "$CONFIG"
if [ $? -ne 0 ]; then
    echo "ERROR: Sweep failed!"
    exit 1
fi
echo ""

echo "=========================================="
echo "  Sweep completed successfully!"
echo "=========================================="
echo ""
