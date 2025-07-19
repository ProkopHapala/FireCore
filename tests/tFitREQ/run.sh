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
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

#> FitREQ_debug.xyz
#python3 run.py
#python3 fit_manual.py
#python3 fit_manual_2.py
#python3 fit_manual_2d.py
#python3 fit_manual_OH.py

#python3 fit_manual_samp.py

stty cols 200   # set terminal width

echo "#=========== RUN"
#> debug.xyz
#python3 -u opt_mini.py
#python3 -u opt_mini.py 2> asan.log | tee OUT
#python3 -u opt_2D.py 
#python3 -u opt_check_derivs.py
#python3 -u opt_check_derivs.py 2> asan.log | tee OUT
#python3 opt_2D.py


echo "Current PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"

python3.12 opt_2D.py  2> asan.log | tee OUT
#python3 -u opt_2D.py 2> asan.log | tee OUT
#python3 -u opt_2D_2.py 2> asan.log | tee OUT
#python3 plot_DOF_trj.py #2> void | tee OUT
#python3 opt_mini.py