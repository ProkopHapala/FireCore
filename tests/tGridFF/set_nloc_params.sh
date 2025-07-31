#!/bin/bash

# Script to modify nloc parameters in OCL_MM.h
# Usage: ./set_nloc_params.sh [parameter_name] [value]
# Example: ./set_nloc_params.sh nlocNBFF 64
# Example: ./set_nloc_params.sh nlocSurf 512
# Example: ./set_nloc_params.sh all 64

OCL_MM_PATH="/home/kocimil1/FireCore/cpp/common/OpenCL/OCL_MM.h"
RELAX_MULTI_PATH="/home/kocimil1/FireCore/cpp/common_resources/cl/relax_multi.cl"

# Make a backup of the original files
cp "$OCL_MM_PATH" "${OCL_MM_PATH}.bak"
cp "$RELAX_MULTI_PATH" "${RELAX_MULTI_PATH}.bak"

if [ "$#" -lt 2 ]; then
    echo "Usage: $0 [parameter_name] [value]"
    echo "Available parameters:"
    echo "  nlocNBFF         - Non-bonded force field (default: 32)"
    echo "  nlocGridFF       - Grid force field (default: 32)"
    echo "  nlocGridFFbSpline - Grid force field B-spline (default: 32)"
    echo "  nlocMMFF         - Molecular mechanics force field (default: 1)"
    echo "  nlocmove         - Move atoms (default: 1)"
    echo "  nlocSurf         - Surface parameters (default: 32)"
    echo "  all              - Change all parameters"
    exit 1
fi

PARAM=$1
VALUE=$2

update_param() {
    local param=$1
    local value=$2
    local pattern=$3
    local comment=$4
    
    # Special case for nlocmove which has a different format
    if [ "$param" = "nlocmove" ]; then
        sed -i "s/int nloc=[0-9]\+; \/\/nlocmove/int nloc=$value; \/\/nlocmove/" "$OCL_MM_PATH"
        if grep -q "nloc=$value.*nlocmove" "$OCL_MM_PATH"; then
            echo "Updated $param to $value in OCL_MM.h"
            return 0
        fi
    fi
    
    # Find the line with the pattern and update it in OCL_MM.h
    if [ -n "$comment" ]; then
        # Nejprve najděme příslušnou funkci, abychom se ujistili, že upravujeme správný parametr
        if [ "$param" = "nlocGridFF" ]; then
            # Hledáme konkrétně v setup_getNonBond_GridFF
            sed -i "/setup_getNonBond_GridFF(/,/task->local.x/s/int\s\+nloc\s\+=\s\+[0-9]\+;\s\+\/\/\s*$comment/int nloc = $value; \/\/$comment/" "$OCL_MM_PATH"
        elif [ "$param" = "nlocGridFFbSpline" ]; then
            # Hledáme konkrétně v setup_getNonBond_GridFF_Bspline
            sed -i "/setup_getNonBond_GridFF_Bspline(/,/task->local.x/s/int\s\+nloc\s\+=\s\+[0-9]\+;\s\+\/\/\s*$comment/int nloc = $value; \/\/$comment/" "$OCL_MM_PATH"
        else
            # Pro ostatní parametry používáme obecný vzor
            sed -i "s/int\s\+nloc\s\+=\s\+[0-9]\+;\s\+\/\/\s*$comment/int nloc = $value; \/\/$comment/" "$OCL_MM_PATH"
        fi
    else
        # Více flexibilní vzor pro řádky bez komentářů
        sed -i "s/int\s\+nloc\s\+=\s\+[0-9]\+;/int nloc = $value;/" "$OCL_MM_PATH"
    fi
    
    # Check if the update was successful in OCL_MM.h
    if grep -q "nloc.*=.*$value.*$comment" "$OCL_MM_PATH"; then
        echo "Updated $param to $value in OCL_MM.h"
    else
        echo "Failed to update $param in OCL_MM.h. Trying alternative method..."
        # Try a more direct approach by looking for the specific pattern
        case "$param" in
            "nlocSurf")
                # Specifically target the getSurfMorse function's nloc line
                sed -i "/getSurfMorse/,/task->local.x/s/int\s\+nloc\s\+=\s\+[0-9]\+;\s\+\/\/\s*nlocSurf/int nloc = $value; \/\/nlocSurf/" "$OCL_MM_PATH"
                if grep -q "nloc.*=.*$value.*nlocSurf" "$OCL_MM_PATH"; then
                    echo "Successfully updated $param to $value in OCL_MM.h with alternative method"
                else
                    echo "All update attempts failed for $param in OCL_MM.h"
                fi
                
                # Also update the LATOMS and LCLJS array sizes in relax_multi.cl
                update_relax_multi_arrays "$value" "getSurfMorse"
                ;;
            "nlocmove")
                # Try another pattern for nlocmove
                sed -i "/setup_updateAtomsMMFFf4/,/task->local.x/s/int nloc=[0-9]\+;/int nloc=$value;/" "$OCL_MM_PATH"
                if grep -q "nloc=$value" "$OCL_MM_PATH"; then
                    echo "Successfully updated $param to $value in OCL_MM.h with alternative method"
                else
                    echo "All update attempts failed for $param in OCL_MM.h"
                fi
                ;;
        esac
    fi
    
    # For nlocSurf we already updated the relax_multi.cl file in the case statement above
    # No need to do it again here
}

# Function to update LATOMS and LCLJS array sizes in a specific kernel function
update_relax_multi_arrays() {
    local value=$1
    local kernel_name=$2
    local updated=false
    
    echo "Searching for kernel '${kernel_name}'..."
    
    # Create a simpler approach - create a marker file for each parameter
    # Also reset the values for getNonBond_GridFF to 16
    sed -i "/^__kernel void getNonBond_GridFF(/,/__local float4 LATOMS\[[0-9]\+\];/s/__local float4 LATOMS\[[0-9]\+\];/__local float4 LATOMS[16];/" "$RELAX_MULTI_PATH"
    sed -i "/^__kernel void getNonBond_GridFF(/,/__local float4 LCLJS \[[0-9]\+\];/s/__local float4 LCLJS \[[0-9]\+\];/__local float4 LCLJS [16];/" "$RELAX_MULTI_PATH"
    echo "Reset getNonBond_GridFF arrays to 16"

    case "$kernel_name" in
        "getNonBond")
            # For getNonBond, we only want to update the exact getNonBond function,
            # not getNonBond_GridFF or getNonBond_GridFF_Bspline
            
            # First update - around line 1055
            if grep -A 20 "^__kernel void getNonBond(" "$RELAX_MULTI_PATH" | grep -q "__local float4 LATOMS\["; then
                sed -i "/^__kernel void getNonBond(/,/__local float4 LATOMS\[[0-9]\+\];/s/__local float4 LATOMS\[[0-9]\+\];/__local float4 LATOMS[$value];/" "$RELAX_MULTI_PATH"
                sed -i "/^__kernel void getNonBond(/,/__local float4 LCLJS \[[0-9]\+\];/s/__local float4 LCLJS \[[0-9]\+\];/__local float4 LCLJS [$value];/" "$RELAX_MULTI_PATH"
                updated=true
                echo "Updated first instance of getNonBond"
            fi
            
            # Second update - around line 1642
            if grep -A 20 "//__attribute.*__kernel void getNonBond(" "$RELAX_MULTI_PATH" | grep -q "__local float4 LATOMS\["; then
                line_num=$(grep -n "//__attribute.*__kernel void getNonBond(" "$RELAX_MULTI_PATH" | cut -d: -f1)
                if [ -n "$line_num" ]; then
                    start_line=$((line_num))
                    end_line=$((line_num + 50))
                    sed -i "${start_line},${end_line}s/__local float4 LATOMS\[[0-9]\+\];/__local float4 LATOMS[$value];/" "$RELAX_MULTI_PATH"
                    sed -i "${start_line},${end_line}s/__local float4 LCLJS \[[0-9]\+\];/__local float4 LCLJS [$value];/" "$RELAX_MULTI_PATH"
                    updated=true
                    echo "Updated second instance of getNonBond"
                fi
            fi
            ;;
            
        "getSurfMorse")
            # Only update getSurfMorse function
            if grep -A 20 "^__kernel void getSurfMorse(" "$RELAX_MULTI_PATH" | grep -q "__local float4 LATOMS\["; then
                sed -i "/^__kernel void getSurfMorse(/,/__local float4 LATOMS\[[0-9]\+\];/s/__local float4 LATOMS\[[0-9]\+\];/__local float4 LATOMS[$value];/" "$RELAX_MULTI_PATH"
                sed -i "/^__kernel void getSurfMorse(/,/__local float4 LCLJS \[[0-9]\+\];/s/__local float4 LCLJS \[[0-9]\+\];/__local float4 LCLJS [$value];/" "$RELAX_MULTI_PATH"
                updated=true
                echo "Updated getSurfMorse"
            fi
            ;;
            
        "getNonBond_GridFF_Bspline")
            # Only update getNonBond_GridFF_Bspline function
            if grep -A 20 "^__kernel void getNonBond_GridFF_Bspline(" "$RELAX_MULTI_PATH" | grep -q "__local float4 LATOMS\["; then
                sed -i "/^__kernel void getNonBond_GridFF_Bspline(/,/__local float4 LATOMS\[[0-9]\+\];/s/__local float4 LATOMS\[[0-9]\+\];/__local float4 LATOMS[$value];/" "$RELAX_MULTI_PATH"
                sed -i "/^__kernel void getNonBond_GridFF_Bspline(/,/__local float4 LCLJS \[[0-9]\+\];/s/__local float4 LCLJS \[[0-9]\+\];/__local float4 LCLJS [$value];/" "$RELAX_MULTI_PATH"
                updated=true
                echo "Updated getNonBond_GridFF_Bspline"
            fi
            ;;
            
        *)
            # Generic handling for other functions
            echo "Using generic handling for $kernel_name"
            if grep -A 20 "^__kernel void ${kernel_name}(" "$RELAX_MULTI_PATH" | grep -q "__local float4 LATOMS\["; then
                sed -i "/^__kernel void ${kernel_name}(/,/__local float4 LATOMS\[[0-9]\+\];/s/__local float4 LATOMS\[[0-9]\+\];/__local float4 LATOMS[$value];/" "$RELAX_MULTI_PATH"
                sed -i "/^__kernel void ${kernel_name}(/,/__local float4 LCLJS \[[0-9]\+\];/s/__local float4 LCLJS \[[0-9]\+\];/__local float4 LCLJS [$value];/" "$RELAX_MULTI_PATH"
                updated=true
                echo "Updated $kernel_name"
            fi
            ;;
    esac
    
    if [ "$updated" = true ]; then
        echo "Updated LATOMS and LCLJS array sizes to $value in ${kernel_name} function(s)"
        return 0
    else
        echo "No updates made to ${kernel_name} function(s)"
        return 1
    fi
}

# Function to update all params
update_all() {
    local value=$1
    
    # Update all nloc parameters to the same value
    sed -i "s/int nloc = [0-9]\+;/int nloc = $value;/g" "$OCL_MM_PATH"
    sed -i "s/int nloc=[0-9]\+;/int nloc=$value;/g" "$OCL_MM_PATH"
    sed -i "s/int nloc = [0-9]\+; \/\/nloc/int nloc = $value; \/\/nloc/g" "$OCL_MM_PATH"
    
    # Also update LATOMS and LCLJS array sizes in relax_multi.cl for each main kernel
    update_relax_multi_arrays "$value" "getSurfMorse"
    update_relax_multi_arrays "$value" "getNonBond_GridFF_Bspline"
    # Pridejte dalsi kernely, ktere potrebujete aktualizovat pri "all"
    
    echo "Updated all nloc parameters to $value"
}

# Main switch case to handle different parameters
case $PARAM in
    "nlocNBFF")
        update_param "nlocNBFF" "$VALUE" "setup_getNonBond" "nlocNBFF"
        update_relax_multi_arrays "$VALUE" "getNonBond"
        ;;
    "nlocGridFF")
        update_param "nlocGridFF" "$VALUE" "setup_getNonBond_GridFF" "nlocGridFF"
        
        # Aktualizace polí v getNonBond_GridFF
        echo "Updating arrays in getNonBond_GridFF..."
        update_relax_multi_arrays "$VALUE" "getNonBond_GridFF"
        ;;
    "nlocGridFFbSpline")
            update_param "nlocGridFFbSpline" "$VALUE" "setup_getNonBond_GridFF_Bspline" "nlocGridFFbSpline"
            
            # Manuální úprava polí v getNonBond_GridFF_Bspline
            echo "Manually updating arrays in getNonBond_GridFF_Bspline..."
            
            # Tento kernel je specifický svým formátem a nemá standardní deklaraci
            # Musíme přesně cílit na řádky 2331 a 2332
            
            # Přímá úprava souboru na daných řádcích (bez regexu)
            # Tento kernel je specifický svým formátem
            # Cílíme na řádky po deklaraci kernelu
            
            # Zjistíme číslo řádku, kde začíná deklarace kernelu getNonBond_GridFF_Bspline
            start_line=$(grep -n "getNonBond_GridFF_Bspline" "$RELAX_MULTI_PATH" | head -1 | cut -d: -f1)
            
            if [ -z "$start_line" ]; then
                echo "Error: Could not find kernel getNonBond_GridFF_Bspline in $RELAX_MULTI_PATH"
                return 1
            fi
            
            # Hledáme řádky s definicemi lokálních polí po tomto řádku
            latoms_line=$(grep -n "__local float4 LATOMS" "$RELAX_MULTI_PATH" | awk -v start=$start_line '$1 > start' | head -1 | cut -d: -f1)
            lcljs_line=$(grep -n "__local float4 LCLJS" "$RELAX_MULTI_PATH" | awk -v start=$start_line '$1 > start' | head -1 | cut -d: -f1)
            
            echo "Found GridFFbSpline kernel at line $start_line"
            echo "Found LATOMS at line $latoms_line and LCLJS at line $lcljs_line"
            
            if [ -n "$latoms_line" ]; then
                sed -i "${latoms_line}c\\    __local float4 LATOMS[$VALUE];         // local memory chumk of positions of atoms " "$RELAX_MULTI_PATH"
                echo "Updated LATOMS in getNonBond_GridFF_Bspline to size $VALUE"
            else
                echo "Could not find LATOMS definition"
            fi
            
            if [ -n "$lcljs_line" ]; then
                sed -i "${lcljs_line}c\\    __local float4 LCLJS [$VALUE];         // local memory chumk of atom parameters" "$RELAX_MULTI_PATH"
                echo "Updated LCLJS in getNonBond_GridFF_Bspline to size $VALUE"
            else
                echo "Could not find LCLJS definition"
            fi
            
            echo "Updated arrays in getNonBond_GridFF_Bspline to $VALUE"
            ;;
    "nlocMMFF")
        update_param "nlocMMFF" "$VALUE" "setup_getMMFFf4" "nlocMMFF"
        ;;
    "nlocmove")
        update_param "nlocmove" "$VALUE" "setup_updateAtomsMMFFf4" "nlocmove"
        ;;
    "nlocSurf")
        update_param "nlocSurf" "$VALUE" "getSurfMorse" "nlocSurf"
        update_relax_multi_arrays "$VALUE" "getSurfMorse"
        ;;
    "all")
        update_all "$VALUE"
        ;;
    *)
        echo "Unknown parameter: $PARAM"
        echo "Use one of: nlocNBFF, nlocGridFF, nlocGridFFbSpline, nlocMMFF, nlocmove, nlocSurf, all"
        exit 1
        ;;
esac

echo "Done. Original files backed up as ${OCL_MM_PATH}.bak and ${RELAX_MULTI_PATH}.bak"
