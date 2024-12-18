#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources

name=OCL_GridFF

wd=`pwd`
cd ../../cpp/Build/libs_OCL/
pwd
rm lib$name.so
make $name
cd $wd


# ---- Multiprocesing
#ncpu=`nproc`
#ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
#echo "compile using ncpu="$ncpu
#OMP_NUM_THREADS=$ncpu
#export OMP_NUM_THREADS

#rm *.bin

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD

python3 run.py

