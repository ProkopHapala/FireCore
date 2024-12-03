#!/bin/bash

wd=`pwd`

echo "#=========== Compile C++"
cd ../../cpp/Build/libs/Molecular/
pwd
rm   libSchroedingerGreen2D_lib.so
make SchroedingerGreen2D_lib
cd $wd

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

echo "#=========== RUN"
python3 run.py