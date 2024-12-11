#!/bin/bash

wd=`pwd`

#LD_LIBRARY_PATH=/home/prokop/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
#LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/11/:$LD_LIBRARY_PATH


#echo "#=========== Compile C++"
cd ../../cpp/Build/libs/Molecular/
pwd
rm   libFitREQ_lib.so
make -j4 FitREQ_lib
cd $wd

# ------- asan (Memory Sanitizer)
#LD_PRELOAD=$(g++ -print-file-name=libasan.so)
#echo   $LD_PRELOAD
#export LD_PRELOAD


stty cols 200   # set terminal width

echo "#=========== RUN"
#> debug.xyz
#python3 -u opt_mini.py
python3 -u opt_mini.py 2> asan.log | tee OUT
python3 plot_DOF_trj.py
#python3 opt_mini.py