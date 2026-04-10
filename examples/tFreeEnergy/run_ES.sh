#!/bin/bash

set -e

cd "$(dirname "$0")"

if [ ! -e data ]; then
    ln -s ../../cpp/common_resources data
fi
if [ ! -e common_resources ]; then
    ln -s ../../cpp/common_resources common_resources
fi

N=30
DEBUG_TI=0
LAMDA1=0.5
LAMDA2=11.0
NBSTEP=100
NMDSTEPS=10000
NEQSTEPS=5000
DT=0.05
TDAMP=100
TEMP=300

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --debug_TI) DEBUG_TI=1 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
    shift
done

if [[ "$DEBUG_TI" -eq 1 ]]; then
    echo "Step 1: Building libMMFFmulti_lib.so..."
    cd ../../cpp/Build/libs_OCL
    make MMFFmulti_lib
    cd - >/dev/null
    echo "Build successful."
    echo

    echo "Step 2: Running GPU debug entropic-spring TI..."
    python3 run_ES_gpu_debug.py \
        --nSys 1 \
        --xyz_name "../tMMFF/data/entropic_spring_${N}.xyz" \
        --constr_name "../tMMFF/data/entropic_spring_${N}.cons" \
        --constraints "constraints_ES.txt" \
        --lamda1 ${LAMDA1} \
        --lamda2 ${LAMDA2} \
        --nbStep ${NBSTEP} \
        --nMDsteps ${NMDSTEPS} \
        --nEQsteps ${NEQSTEPS} \
        --dt ${DT} \
        --t_damp ${TDAMP} \
        -T ${TEMP}
    echo

    echo "Output files:"
    echo "  - results/entropic_spring_${N}_TI_gpu_debug_T${TEMP}K.dat"
    echo "  - results/entropic_spring_${N}_TI_gpu_debug_T${TEMP}K.png"
    exit 0
fi

echo "Step 1: Building libMMFF_lib.so..."
cd ../../cpp/Build/libs/Molecular
make MMFF_lib
cd - >/dev/null
echo "Build successful."
echo

echo "Step 2: Running CPU entropic-spring TI..."
python3 run_ES.py \
    --xyz_name "../tMMFF/data/entropic_spring_${N}.xyz" \
    --constr_name "../tMMFF/data/entropic_spring_${N}.cons" \
    --constraints "constraints_ES.txt" \
    --lamda1 ${LAMDA1} \
    --lamda2 ${LAMDA2} \
    --nbStep ${NBSTEP} \
    --nMDsteps ${NMDSTEPS} \
    --nEQsteps ${NEQSTEPS} \
    --dt ${DT} \
    --t_damp ${TDAMP} \
    -T ${TEMP}
echo

echo "Output files:"
echo "  - results/entropic_spring_${N}_TI_cpu_T${TEMP}K.dat"
echo "  - results/entropic_spring_${N}_TI_cpu_T${TEMP}K.png"
