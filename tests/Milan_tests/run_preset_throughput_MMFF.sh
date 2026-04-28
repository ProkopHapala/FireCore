#!/bin/bash

# This script calls the parametric version of the throughput test
# with a preset configuration.

bash run_throughput_MMFF.sh \
    --bNonBonded 1 \
    --doSurfAtoms 1 \
    --bGridFF 6 \
    --bTex 0 \
    --bSaveToDatabase -1 \
    --Fconv 1e-4 \
    --replicas "5000" \
    --perframes "100" \
    --perVF "100" \
    --nlocMMFFs "32" \
    --nlocmoves "32" \
    --nlocNBFFs "--" \
    --nlocSurfs "--" \
    --nlocGridFFs "--" \
    --nlocGridFFbSplines "--" \
    --nPBC "(1,1,0)" \
    --xyz_name "data/xyz/xylitol_WO_gridFF" \
    --surf_name_template "data/xyz/surfaces_for_throughput/NaCl_%dx%d_Cl_hole" \
    --Ns "16" \
    --elapse_time 5.0
