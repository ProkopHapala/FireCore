#!/bin/bash

# Create symbolic links if they don't exist
if [ ! -d data ]; then
    ln -s ../../cpp/common_resources data
fi
if [ ! -d common_resources ]; then
    ln -s ../../cpp/common_resources common_resources
fi

# Get working directory
wd=`pwd`

# Default values
bNonBonded=1
doSurfAtoms=1
bGridFF=6
bTex=0
bSaveToDatabase=-1
Fconv=1e-4
replicas_str="5000"
perframes_str="100"
perVF_str="100"
nlocMMFFs_str="32"
nlocmoves_str="32"
nlocNBFFs_str="--"
nlocSurfs_str="--"
nlocGridFFs_str="--"
nlocGridFFbSplines_str="--"
nPBC_str="(1,1,0)"
xyz_name="data/xyz/xylitol_WO_gridFF"
surf_name_template="data/xyz/surfaces_for_throughput/NaCl_%dx%d_Cl_hole"
Ns_str="16"
elapse_time=0.0

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bNonBonded) bNonBonded="$2"; shift ;;
        --doSurfAtoms) doSurfAtoms="$2"; shift ;;
        --bGridFF) bGridFF="$2"; shift ;;
        --bTex) bTex="$2"; shift ;;
        --bSaveToDatabase) bSaveToDatabase="$2"; shift ;;
        --Fconv) Fconv="$2"; shift ;;
        --replicas) replicas_str="$2"; shift ;;
        --perframes) perframes_str="$2"; shift ;;
        --perVF) perVF_str="$2"; shift ;;
        --nlocMMFFs) nlocMMFFs_str="$2"; shift ;;
        --nlocmoves) nlocmoves_str="$2"; shift ;;
        --nlocNBFFs) nlocNBFFs_str="$2"; shift ;;
        --nlocSurfs) nlocSurfs_str="$2"; shift ;;
        --nlocGridFFs) nlocGridFFs_str="$2"; shift ;;
        --nlocGridFFbSplines) nlocGridFFbSplines_str="$2"; shift ;;
        --nPBC) nPBC_str="$2"; shift ;;
        --xyz_name) xyz_name="$2"; shift ;;
        --surf_name_template) surf_name_template="$2"; shift ;;
        --Ns) Ns_str="$2"; shift ;;
        --elapse_time) elapse_time="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

lscpu

# Set up multithreading
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu=$ncpu"
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# Create results directory if needed
mkdir -p results/best
mkdir -p results/all

# Setup parameters
if (( bGridFF == 1 )); then
    bSpline=0
elif (( bGridFF == 5 || bGridFF == 6 )); then
    bSpline=1
fi

# Generate bit number from flags
bNonBonded_bit=$(( bNonBonded > 0 ? 1 : 0 ))
doSurfAtoms_bit=$(( doSurfAtoms > 0 ? 1 : 0 ))
bGridFF_bit=$(( bGridFF > 0 ? 1 : 0 ))
bSpline_bit=$(( bSpline > 0 ? 1 : 0 ))
bTex_bit=$(( bTex > 0 ? 1 : 0 ))
flags_bitnum=$(( (bNonBonded_bit << 4) | (doSurfAtoms_bit << 3) | (bGridFF_bit << 2) | (bSpline_bit << 1) | bTex_bit ))

# Convert string inputs to arrays
read -r -a replicas <<< "$replicas_str"
read -r -a perframes <<< "$perframes_str"
read -r -a perVF <<< "$perVF_str"
read -r -a nlocMMFFs <<< "$nlocMMFFs_str"
read -r -a nlocmoves <<< "$nlocmoves_str"
read -r -a nlocNBFFs <<< "$nlocNBFFs_str"
read -r -a nlocSurfs <<< "$nlocSurfs_str"
read -r -a nlocGridFFs <<< "$nlocGridFFs_str"
read -r -a nlocGridFFbSplines <<< "$nlocGridFFbSplines_str"
read -r -a nPBC <<< "$nPBC_str"
read -r -a Ns <<< "$Ns_str"

# Initialize best result tracking
best_value=0
best_name=""
best_line=""
best_nSys=0
best_perframe=0
best_perVF=0

# Create or clear the results file
echo "Parameter optimization results" > fastest_parameters_results.dat

for nPBC in "${nPBC[@]}"; do
    # Loop through all parameter combinations
    for nlocMMFF in "${nlocMMFFs[@]}"; do
        ./set_nloc_params.sh nlocMMFF $nlocMMFF
        for nlocmove in "${nlocmoves[@]}"; do
            ./set_nloc_params.sh nlocmove $nlocmove
            for nlocNBFF in "${nlocNBFFs[@]}"; do
                if [ "$nlocNBFF" != "--" ]; then ./set_nloc_params.sh nlocNBFF $nlocNBFF; fi
                for nlocSurf in "${nlocSurfs[@]}"; do
                    if [ "$nlocSurf" != "--" ]; then ./set_nloc_params.sh nlocSurf $nlocSurf; fi
                    for nlocGridFF in "${nlocGridFFs[@]}"; do
                        if [ "$nlocGridFF" != "--" ]; then ./set_nloc_params.sh nlocGridFF $nlocGridFF; fi
                        for nlocGridFFbSpline in "${nlocGridFFbSplines[@]}"; do
                            if [ "$nlocGridFFbSpline" != "--" ]; then ./set_nloc_params.sh nlocGridFFbSpline $nlocGridFFbSpline; fi
                            for nSys in "${replicas[@]}"; do
                                for perframe in "${perframes[@]}"; do
                                    for pvf in "${perVF[@]}"; do  
                                        # Build the library
                                        cd ../../cpp/Build/libs_OCL/
                                        pwd
                                        rm -f libMMFFmulti_lib.so
                                        make MMFFmulti_lib
                                        cd $wd 
                                        if (( pvf > perframe )); then continue; fi       
                                        # Run simulation with the current parameters
                                        if (( doSurfAtoms > 0 )); then
                                            echo "======================================================="
                                            echo "Testing with nSys=$nSys, perframe=$perframe, perVF=$pvf"
                                            
                                            # Remove previous results
                                            rm -f minima.dat gopt.xyz
                                            touch minima.dat
                                            
                                            for N in "${Ns[@]}"; do
                                                surf_name=$(printf "$surf_name_template" "$N" "$N")
                                                python3 run_throughput_MD.py \
                                                    --nSys "$nSys" \
                                                    --xyz_name "$xyz_name" \
                                                    --surf_name "$surf_name" \
                                                    --bNonBonded "$bNonBonded" \
                                                    --doSurfAtoms "$doSurfAtoms" \
                                                    --GridFF "$bGridFF" \
                                                    --Fconv "$Fconv" \
                                                    --perframe "$perframe" \
                                                    --perVF "$pvf" \
                                                    --gridnPBC "$nPBC" \
                                                    --elapse_time "$elapse_time"
                                                            
                                                # Generate name for the result file
                                                name="minima__$(printf '%04d' $(echo "obase=2;$flags_bitnum" | bc))_surf:_NaCl_${N}x${N}_nPBC_${nPBC}_nloc:_MMFF_${nlocMMFF}_move_${nlocmove}_NBFF_${nlocNBFF}_surf_${nlocSurf}_gridFF_${nlocGridFF}_gridFFbSpline_${nlocGridFFbSpline}___replica:_${nSys}_perframe:_${perframe}_perVF:_${pvf}"
                                                
                                                # Move minima file to results/all directory
                                                mv minima.dat results/all/${name}.dat
                                                
                                                # Extract the last line from the results file
                                                last_line=$(tail -n 1 results/all/${name}.dat)
                                                
                                                # Append results to the all results file
                                                echo "$name $last_line" >> results/all/results.dat
                                                
                                                # Extract value from column 6 (assuming space-separated columns and 0-based indexing)
                                                current_value=$(echo "$last_line" | awk '{print $6}')
                                                
                                                echo "Current value: $current_value, Best so far: $best_value"
                                                
                                                # Check if this is the best result so far
                                                if (( $(echo "$current_value > $best_value" | bc -l) )); then
                                                    best_value=$current_value
                                                    best_name=$name
                                                    best_line="$name $last_line"
                                                    best_nSys=$nSys
                                                    best_perframe=$perframe
                                                    best_perVF=$pvf
                                                    
                                                    echo "New best result found: $best_value with nSys=$nSys, perframe=$perframe, perVF=$pvf"
                                                fi
                                            done
                                        else
                                            echo "======================================================="
                                            echo "Testing with nSys=$nSys, perframe=$perframe, perVF=$pvf"
                                            
                                            # Remove previous results
                                            rm -f minima.dat gopt.xyz
                                            touch minima.dat
                                            N=0
                                            python3 run_throughput_MD.py \
                                                --nSys "$nSys" \
                                                --xyz_name "$xyz_name" \
                                                --bNonBonded "$bNonBonded" \
                                                --GridFF -1 \
                                                --Fconv "$Fconv" \
                                                --perframe "$perframe" \
                                                --perVF "$pvf" \
                                                --elapse_time "$elapse_time"

                                            # Generate name for the result file
                                            name="minima__$(printf '%04d' $(echo "obase=2;$flags_bitnum" | bc))_surf:_NaCl_${N}x${N}_nPBC_${nPBC}_nloc:_MMFF_${nlocMMFF}_move_${nlocmove}_NBFF_${nlocNBFF}_surf_${nlocSurf}_gridFF_${nlocGridFF}_gridFFbSpline_${nlocGridFFbSpline}___replica:_${nSys}_perframe:_${perframe}_perVF:_${pvf}"
                                            
                                            # Move minima file to results/all directory
                                            mv minima.dat results/all/${name}.dat
                                            
                                            # Extract the last line from the results file
                                            last_line=$(tail -n 1 results/all/${name}.dat)
                                            
                                            # Append results to the all results file
                                            echo "$name $last_line" >> results/all/results.dat
                                            
                                            # Extract value from column 6 (assuming space-separated columns and 0-based indexing)
                                            current_value=$(echo "$last_line" | awk '{print $6}')
                                            
                                            echo "Current value: $current_value, Best so far: $best_value"
                                            
                                            # Check if this is the best result so far
                                            if (( $(echo "$current_value > $best_value" | bc -l) )); then
                                                best_value=$current_value
                                                best_name=$name
                                                best_line="$name $last_line"
                                                best_nSys=$nSys
                                                best_perframe=$perframe
                                                best_perVF=$pvf
                                                
                                                echo "New best result found: $best_value with nSys=$nSys, perframe=$perframe, perVF=$pvf"
                                            fi
                                        fi
                                    done
                                done
                            done
                        done
                    done
                done
            done
        done
    done
done

cp results/all/${best_name}.dat results/best/
                
echo "$best_line" >> results/best/fastest_parameters_results.dat
echo "======================================================="
echo "Best result:"
echo "$best_line"
echo "Parameters: replica=$best_nSys, perframe=$best_perframe, perVF=$best_perVF"
echo "Value: $best_value"
echo "Results saved to fastest_parameters_results.dat"
