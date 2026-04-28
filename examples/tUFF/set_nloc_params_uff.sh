#!/bin/bash

# Script to modify nloc parameters in OCL_UFF.h
# Usage: ./set_nloc_params_uff.sh [parameter_name] [value]
# Example: ./set_nloc_params_uff.sh nlocNBFF 64
# Example: ./set_nloc_params_uff.sh nlocGridFFbSpline 128
# Example: ./set_nloc_params_uff.sh all 64

OCL_UFF_PATH="../../cpp/common/OpenCL/OCL_UFF.h"

# Make a backup of the original file
cp "$OCL_UFF_PATH" "${OCL_UFF_PATH}.bak"

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 [parameter_name] [value]"
    echo "Available parameters:"
    echo "  nlocNBFF              - Non-bonded force field (default: 32)"
    echo "  nlocGridFFbSpline     - Grid force field B-spline (default: 32)"
    echo "  nlocSurf              - Surface parameters (default: 32)"
    echo "  nlocBonds             - UFF bonds (default: 32)"
    echo "  nlocAngles            - UFF angles (default: 32)"
    echo "  nlocDihedrals         - UFF dihedrals (default: 32)"
    echo "  nlocInversions        - UFF inversions (default: 32)"
    echo "  nlocAssemble          - Force assembly (default: 32)"
    echo "  nlocClearFapos        - Clear fapos (default: 32)"
    echo "  nlocClearFint         - Clear fint (default: 32)"
    echo "  nlocUpdateAtoms       - Update atoms (default: 32)"
    echo "  all                   - Change all parameters"
    exit 1
fi

PARAM=$1
VALUE=$2

update_param() {
    local param=$1
    local value=$2
    local function_name=$3
    
    echo "Updating $param to $value in function $function_name..."
    
    # Find the function and update the nloc parameter within it
    # We use a range-based sed to only modify within the specific function
    
    # First, find the line number where the function starts
    start_line=$(grep -n "OCLtask\* ${function_name}(" "$OCL_UFF_PATH" | head -1 | cut -d: -f1)
    
    if [ -z "$start_line" ]; then
        # Try alternative pattern for void functions
        start_line=$(grep -n "void ${function_name}(" "$OCL_UFF_PATH" | head -1 | cut -d: -f1)
    fi
    
    if [ -z "$start_line" ]; then
        echo "Error: Could not find function $function_name in $OCL_UFF_PATH"
        return 1
    fi
    
    # Find the closing brace of the function (approximate - find next function or end)
    end_line=$(tail -n +$((start_line + 1)) "$OCL_UFF_PATH" | grep -n "^    OCLtask\*\|^    void\|^};" | head -1 | cut -d: -f1)
    
    if [ -z "$end_line" ]; then
        end_line=50  # Default range if we can't find the end
    else
        end_line=$((start_line + end_line))
    fi
    
    echo "  Function range: lines $start_line to $end_line"
    
    # Update nloc within this range
    sed -i "${start_line},${end_line}s/int nloc\s*=\s*[0-9]\+;/int nloc = $value;/" "$OCL_UFF_PATH"
    
    # Verify the update
    if sed -n "${start_line},${end_line}p" "$OCL_UFF_PATH" | grep -q "int nloc = $value;"; then
        echo "  Successfully updated $param to $value"
        return 0
    else
        echo "  Warning: Could not verify update for $param"
        return 1
    fi
}

# Function to update all params
update_all() {
    local value=$1
    
    echo "Updating all nloc parameters to $value..."
    
    # Update all nloc parameters in the file
    sed -i "s/int nloc\s*=\s*[0-9]\+;/int nloc = $value;/g" "$OCL_UFF_PATH"
    
    echo "Updated all nloc parameters to $value"
}

# Main switch case to handle different parameters
case $PARAM in
    "nlocNBFF")
        update_param "nlocNBFF" "$VALUE" "setup_getNonBond"
        ;;
    "nlocGridFFbSpline")
        update_param "nlocGridFFbSpline" "$VALUE" "setup_getNonBond_GridFF_Bspline"
        ;;
    "nlocSurf")
        update_param "nlocSurf" "$VALUE" "setup_getSurfMorse"
        ;;
    "nlocBonds")
        update_param "nlocBonds" "$VALUE" "setup_kernels"
        # This will update the evalBonds section
        ;;
    "nlocAngles")
        update_param "nlocAngles" "$VALUE" "setup_kernels"
        # This will update the evalAngles section
        ;;
    "nlocDihedrals")
        update_param "nlocDihedrals" "$VALUE" "setup_kernels"
        # This will update the evalDihedrals section
        ;;
    "nlocInversions")
        update_param "nlocInversions" "$VALUE" "setup_kernels"
        # This will update the evalInversions section
        ;;
    "nlocAssemble")
        update_param "nlocAssemble" "$VALUE" "setup_kernels"
        # This will update the assembleForces section
        ;;
    "nlocClearFapos")
        update_param "nlocClearFapos" "$VALUE" "setup_kernels"
        # This will update the clear_fapos section
        ;;
    "nlocClearFint")
        update_param "nlocClearFint" "$VALUE" "setup_kernels"
        # This will update the clear_fint section
        ;;
    "nlocUpdateAtoms")
        update_param "nlocUpdateAtoms" "$VALUE" "setup_updateAtomsMMFFf4"
        ;;
    "all")
        update_all "$VALUE"
        ;;
    *)
        echo "Unknown parameter: $PARAM"
        echo "Use one of: nlocNBFF, nlocGridFFbSpline, nlocSurf, nlocBonds, nlocAngles, nlocDihedrals, nlocInversions, nlocAssemble, nlocClearFapos, nlocClearFint, nlocUpdateAtoms, all"
        exit 1
        ;;
esac

echo "Done. Original file backed up as ${OCL_UFF_PATH}.bak"

