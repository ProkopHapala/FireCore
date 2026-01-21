#!/bin/bash

# This script runs throughput tests with preset configurations for MMFF and UFF.
# Each test is on a single line for simplicity.
# You can comment out the one you don't want to run.

echo "Running MMFF test"


# echo "Running UFF test"
bash run_throughput.sh --ff "uff" --bNonBonded 1 --doSurfAtoms 1 --bGridFF 6 --bTex 0 --bSaveToDatabase 1 --Fconv 1e-4 --replicas "2000" --perframes "100" --perVF "100" --nlocGridFFbSplines "32" --nPBC "(1,1,0)" --xyz_name_template "common_resources/xyz/molecules_for_throughput/xylitol_%dx%d_grid" --surf_name_template "common_resources/xyz/surfaces_for_throughput/NaCl_%dx%d_Cl_hole" --Ns "3" --elapse_time 10.0 --loops 50000 --dt 0.1


# Working tests
: <<'COMMENT'
# check if even GridFF works
bash run_throughput.sh --ff "mmff" --bNonBonded 1 --doSurfAtoms 1 --bGridFF 6 --bTex 0 --bSaveToDatabase 1 --Fconv 1e-4 --replicas "5000" --perframes "100" --perVF "100" --nlocMMFFs "32" --nlocmoves "32" --nlocNBFFs "--" --nlocSurfs "--" --nlocGridFFs "--" --nlocGridFFbSplines "--" --nPBC "(1,1,0)" --xyz_name_template "data/xyz/xylitol_WO_gridFF" --surf_name_template "data/xyz/surfaces_for_throughput/NaCl_%dx%d_Cl_hole" --Ns "16" --elapse_time 5.0 --loops 10000
# check if surfAtoms work (should be slower)
bash run_throughput.sh --ff "mmff" --bNonBonded 1 --doSurfAtoms 1 --bGridFF -6 --bTex 0 --bSaveToDatabase 1 --Fconv 1e-4 --replicas "500" --perframes "100" --perVF "100" --nlocMMFFs "32" --nlocmoves "32" --nlocNBFFs "--" --nlocSurfs "--" --nlocGridFFs "--" --nlocGridFFbSplines "--" --nPBC "(1,1,0)" --xyz_name_template "data/xyz/xylitol_WO_gridFF" --surf_name_template "data/xyz/surfaces_for_throughput/NaCl_%dx%d_Cl_hole" --Ns "3" --elapse_time 10.0 --loops 10000
COMMENT

echo "Tests completed"