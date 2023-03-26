#!/bin/bash

wd=`pwd`

#LD_LIBRARY_PATH=/home/prokop/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/11/:$LD_LIBRARY_PATH


#echo "#=========== Compile C++"
#cd ../../cpp/Build_OCL/libs/Molecular/
#pwd
#rm libMMFF_lib.so
#make MMFF_lib
#cd $wd

#echo "#=========== Compile Fortran"
#cd ../../build
#pwd
#rm   libFireCore.so
#make libFireCore
#cd $wd

#/usr/lib/gcc/x86_64-linux-gnu/11/libgfortran.so

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

echo "#=========== RUN"
python3 run.py