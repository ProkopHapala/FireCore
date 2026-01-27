#!/bin/bash

# Free Energy Calculation Script for DA.mol2
# This script runs the Python script that computes free energy
# for the DA.mol2 molecule pair

echo "=========================================="
echo "Free Energy Calculation for DA.mol2"
echo "=========================================="
echo ""

# Check if Python script exists
if [ ! -f "run.py" ]; then
    echo "Error: run.py not found!"
    exit 1
fi

# Check if DA.mol2 exists
if [ ! -f "cpp/common_resources/DA.mol2" ]; then
    echo "Error: cpp/common_resources/DA.mol2 not found!"
    exit 1
fi

# Run the Python script
echo "Running free energy calculation..."
echo ""

python3 run.py

# Check exit status
if [ $? -eq 0 ]; then
    echo ""
    echo "=========================================="
    echo "Calculation completed successfully!"
    echo "=========================================="
    echo ""
    echo "Output files generated:"
    ls -lh loaded_DA_*.xyz 2>/dev/null || echo "  (No XYZ files found - check for errors)"
else
    echo ""
    echo "=========================================="
    echo "Calculation failed!"
    echo "=========================================="
    exit 1
fi

