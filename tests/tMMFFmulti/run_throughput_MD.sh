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



lscpu

# Set up multithreading
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# Create results directory if needed
mkdir -p results/best
mkdir -p results/all

# Setup parameters
dovdW=1
doSurfAtoms=1
bGridFF=-6   # 1 for linear, 5or6 for bSpline=1 (works only 6)
if (( bGridFF == 1 )); then
    bSpline=0
elif (( bGridFF == 5 || bGridFF == 6 )); then
    bSpline=1
fi
bTex=0
bSaveToDatabase=-1

# Generate bit number from flags
dovdW_bit=$(( dovdW > 0 ? 1 : 0 ))
doSurfAtoms_bit=$(( doSurfAtoms > 0 ? 1 : 0 ))
bGridFF_bit=$(( bGridFF > 0 ? 1 : 0 ))
bSpline_bit=$(( bSpline > 0 ? 1 : 0 ))
bTex_bit=$(( bTex > 0 ? 1 : 0 ))
flags_bitnum=$(( (dovdW_bit << 4) | (doSurfAtoms_bit << 3) | (bGridFF_bit << 2) | (bSpline_bit << 1) | bTex_bit ))

# # Set XYZ file based on GridFF
# if (( bGridFF > 0 )); then
#     xyz_name="data/xyz/xylitol_for_gridFF"
# else
    xyz_name="data/xyz/xylitol_WO_gridFF"
# fi

# Set force convergence criteria
Fconv=1e-4

# Arrays of parameter values to test
replicas=(1000) # (1000 5000) # (1000 2000 3000 4000 5000)
perframes=(100) # (20 500) # (100 500)
perVF=(100) # (20 50) # (50 100)

# Arrays of local memory parameters to test
nlocMMFFs=(32)
nlocmoves=(32)
nlocNBFFs=("--") # (1 2 4 8 16 32 64 128)
nlocSurfs=("--") # (1 2 4 8 16 32 64 128)
nlocGridFFs=("--")
nlocGridFFbSplines=("--")

nPBC=("(1,1,0)") # ("(1,1,0)" "(2,2,0)" "(3,3,0)" "(4,4,0)" "(5,5,0)")


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
                                            
                                            Ns=(16) # (1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16)
                                            for N in "${Ns[@]}"; do
                                                surf_name="data/xyz/surfaces_for_throughput/NaCl_${N}x${N}_L3"
                                                python3 run_throughput_MD.py \
                                                    --nSys "$nSys" \
                                                    --xyz_name "$xyz_name" \
                                                    --surf_name "$surf_name" \
                                                    --dovdW "$dovdW" \
                                                    --doSurfAtoms "$doSurfAtoms" \
                                                    --GridFF "$bGridFF" \
                                                    --Fconv "$Fconv" \
                                                    --perframe "$perframe" \
                                                    --perVF "$pvf" \
                                                    --gridnPBC "$nPBC"
                                                            
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
                                                --dovdW "$dovdW" \
                                                --GridFF -1 \
                                                --Fconv "$Fconv" \
                                                --perframe "$perframe" \
                                                --perVF "$pvf"

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