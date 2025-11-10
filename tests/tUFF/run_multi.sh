#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources

wd=`pwd`
cd ../../cpp/Build/libs_OCL/
pwd
rm libMMFFmulti_lib.so
make MMFFmulti_lib
cd $wd

#cd ../../cpp/Build/libs_SDL
#rm libMolGUIlib.so
#make MolGUIlib
#cd $wd


# ---- Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

#rm *.bin

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD

lscpu

#python3 run.py

# === UFF CPU/GPU validation presets =========================================
python3 -u test_UFF_multi.py --preset bonded --non-bonded --grid-ff --use-scan-relaxed 0   --nconf 10000          2>&1 | tee OUT-UFF-multi-all
#python3 -u test_UFF_multi.py --preset bonded                                     2>&1 | tee OUT-UFF-multi-bonded
#python3 -u test_UFF_multi.py --preset bonded --non-bonded                        2>&1 | tee OUT-UFF-multi-bonded-nb
#python3 -u test_UFF_multi.py --preset bonded --grid-ff                           2>&1 | tee OUT-UFF-multi-bonded-grid
#python3 -u test_UFF_multi.py --preset none   --non-bonded                        2>&1 | tee OUT-UFF-multi-nonbonded

#python3 -u test_UFF_multi.py --preset grid-only --grid-ff                        2>&1 | tee OUT-UFF-multi-gridff
#python3 -u test_UFF_multi.py --preset none   --non-bonded --grid-ff              2>&1 | tee OUT-UFF-multi-nb-grid
#python3 -u test_UFF_multi.py --preset bonded --non-bonded                        2>&1 | tee OUT-UFF-multi-bonded-nb
