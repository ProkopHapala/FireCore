#!/bin/bash

#echo "#=========== Compile C++"
wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
rm libMMFF_lib.so
make -j4 MMFF_lib
cd $wd

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD

stty cols 1000   # set terminal width


echo "#=========== RUN TEST "
#python3 run.py "./data_UFF/xyz/ethylene" | tee out
python3 -u test_UFF_ocl.py 2>1 | tee OUT-UFF-ocl