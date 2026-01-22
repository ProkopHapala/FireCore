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
# LD_PRELOAD=$(g++ -print-file-name=libasan.so)
# echo   $LD_PRELOAD
# export LD_PRELOAD


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

#python3 -u opt_2D.py    2>&1 | tee OUT-fit
#python3 -u opt_2D_new.py # 2>&1 | tee OUT-fit-new
python3 -u opt_2D_new.py --input /home/niko/work/HBOND/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m/H2O-A1_H2O-D1-y.xyz \
   --mode scan --scan_dofs 1 --scan_range 0.0 1.0 100 2>&1 | tee OUT-fit-multi
#python3 sample_damped_coulomb.py | tee OUT-sample

#python3 -u opt_2D.py 2> asan.log | tee OUT
#python3 -u opt_2D_2.py 2> asan.log | tee OUT
#python3 plot_DOF_trj.py #2> void | tee OUT
#python3 opt_mini.py