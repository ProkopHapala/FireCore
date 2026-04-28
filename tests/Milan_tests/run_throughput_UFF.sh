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

# Setup parameters for UFF
dovdW=1 # DO NOT CHANGE (but has no effect)
doSurfAtoms=1 # negative for no surface, 1 for surface
bGridFF=6   # negative for no GridFF, 6 for bSpline (other values are not tested)
if (( bGridFF == 1 )); then
    bSpline=0
elif (( bGridFF == 5 || bGridFF == 6 )); then
    bSpline=1
fi
bTex=0
bSaveToDatabase=-1
# xyz_name="common_resources/xyz/xylitol_WO_gridFF"

# Generate bit number from flags DO NOT CHANGE
dovdW_bit=$(( dovdW > 0 ? 1 : 0 ))
doSurfAtoms_bit=$(( doSurfAtoms > 0 ? 1 : 0 ))
bGridFF_bit=$(( bGridFF > 0 ? 1 : 0 ))
bSpline_bit=$(( bSpline > 0 ? 1 : 0 ))
bTex_bit=$(( bTex > 0 ? 1 : 0 ))
flags_bitnum=$(( (dovdW_bit << 4) | (doSurfAtoms_bit << 3) | (bGridFF_bit << 2) | (bSpline_bit << 1) | bTex_bit ))

# Set force convergence criteria
Fconv=1e-4

# Arrays of parameter values to test (for article we used replicas=(5000), perframes=(100, 20), perVF=(100, 20))
replicas=(2) # (1000 5000) # (1000 2000 3000 4000 5000)
perframes=(100) # (20 500) # (100 500) #niter
perVF=(1) # (20 50) # (50 100) #nPerVFs

# Arrays of local memory parameters DO NOT CHANGE
nlocNBFFs=("--") # (1 2 4 8 16 32 64 128)
nlocSurfs=("--") # (1 2 4 8 16 32 64 128)
nlocGridFFbSplines=(32) # (1 2 4 8 16 32 64 128)

# Set nPBC for surface
nPBC=("(1,1,0)") # ("(1,1,0)" "(2,2,0)" "(3,3,0)" "(4,4,0)" "(5,5,0)")


# Initialize best result tracking
best_value=0
best_name=""
best_line=""
best_nSys=0
best_perframe=0
best_perVF=0

# Create or clear the results file
echo "UFF Parameter optimization results" > fastest_parameters_results.dat

for nPBC in "${nPBC[@]}"; do
    # Loop through all parameter combinations
    for nlocNBFF in "${nlocNBFFs[@]}"; do
        if [ "$nlocNBFF" != "--" ]; then ./set_nloc_params_uff.sh nlocNBFF $nlocNBFF; fi
        for nlocSurf in "${nlocSurfs[@]}"; do
            if [ "$nlocSurf" != "--" ]; then ./set_nloc_params_uff.sh nlocSurf $nlocSurf; fi
            for nlocGridFFbSpline in "${nlocGridFFbSplines[@]}"; do
                if [ "$nlocGridFFbSpline" != "--" ]; then ./set_nloc_params_uff.sh nlocGridFFbSpline $nlocGridFFbSpline; fi
                for nSys in "${replicas[@]}"; do
                    for perframe in "${perframes[@]}"; do
                        for pvf in "${perVF[@]}"; do
                            # Build the library
                            cd ../../cpp/Build/libs_OCL/
                            pwd
                            rm -f libMMFFmulti_lib.so
                            make MMFFmulti_lib
                            cd $wd

                            # Skip if perVF is greater than perframe
                            if (( pvf > perframe )); then continue; fi

                            if (( doSurfAtoms > 0 )); then
                                # Remove previous results
                                rm -f minima.dat gopt.xyz
                                touch minima.dat

                                Ns=(3) # (3 4 5 6 7 8 9 10 11 12 13 14 15 16)
                                for N in "${Ns[@]}"; do    
                                    if [ $bGridFF -lt 0 ]; then
                                        xyz_name="common_resources/xyz/molecules_for_throughput/xylitol_${N}x${N}"
                                    else
                                        xyz_name="common_resources/xyz/molecules_for_throughput/xylitol_${N}x${N}_grid"
                                    fi
                                    surf_name="common_resources/xyz/surfaces_for_throughput/NaCl_${N}x${N}_Cl_hole"

                                    # Run UFF simulation
                                    python3 run_throughput_UFF.py \
                                        --nSys "$nSys" \
                                        --xyz_name "$xyz_name" \
                                        --surf_name "$surf_name" \
                                        --dovdW "$dovdW" \
                                        --doSurfAtoms "$doSurfAtoms" \
                                        --bGridFF "$bGridFF" \
                                        --Fconv "$Fconv" \
                                        --perframe "$perframe" \
                                        --perVF "$pvf" \
                                        --gridnPBC "$nPBC"

                                    # Generate name for the result file
                                    name="minima__$(printf '%04d' $(echo "obase=2;$flags_bitnum" | bc))_surf:_NaCl_${N}x${N}_nPBC_${nPBC}_nloc:_NBFF_${nlocNBFF}_surf_${nlocSurf}_gridFFbSpline_${nlocGridFFbSpline}___replica:_${nSys}_perframe:_${perframe}_perVF:_${pvf}"

                                    # Move minima file to results/all directory
                                    mv minima.dat results/all/${name}.dat
                                    last_line=$(tail -n 1 results/all/${name}.dat)
                                    echo "$name $last_line" >> results/all/results.dat
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
                            else # doSurfAtoms=0
                                # Remove previous results
                                rm -f minima.dat gopt.xyz
                                touch minima.dat

                                # Define surface size as 0 for no surface
                                N=0
                                xyz_name="common_resources/xyz/molecules_for_throughput/xylitol_1x1"

                                # Run UFF simulation without surface
                                python3 run_throughput_UFF.py \
                                    --nSys "$nSys" \
                                    --xyz_name "$xyz_name" \
                                    --dovdW "$dovdW" \
                                    --bGridFF -1 \
                                    --Fconv "$Fconv" \
                                    --perframe "$perframe" \
                                    --perVF "$pvf"

                                # Generate name for the result file
                                name="minima__$(printf '%04d' $(echo "obase=2;$flags_bitnum" | bc))_surf:_NaCl_${N}x${N}_nPBC_${nPBC}_nloc:_NBFF_${nlocNBFF}_surf_${nlocSurf}_gridFFbSpline_${nlocGridFFbSpline}___replica:_${nSys}_perframe:_${perframe}_perVF:_${pvf}"
                                mv minima.dat results/all/${name}.dat
                                last_line=$(tail -n 1 results/all/${name}.dat)
                                echo "$name $last_line" >> results/all/results.dat
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

if [ -n "$best_name" ]; then
    cp results/all/${best_name}.dat results/best/
    echo "$best_line" >> results/best/fastest_parameters_results.dat
    echo "======================================================="
    echo "Best result:"
    echo "$best_line"
    echo "Parameters: replica=$best_nSys, perframe=$best_perframe, perVF=$best_perVF"
    echo "Time: $best_value"
    echo "Results saved to fastest_parameters_results.dat"
else
    echo "No valid results found"
fi