#!/bin/bash

wd=`pwd`

ln -s ../../cpp/common_resources data

#LD_LIBRARY_PATH=/home/prokop/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
#LD_LIBRARY_PATH=/usr/lib/gcc/x86_64-linux-gnu/11/:$LD_LIBRARY_PATH

#export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libmkl_def.so:/usr/lib/x86_64-linux-gnu/libmkl_avx2.so:/usr/lib/x86_64-linux-gnu/libmkl_core.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_lp64.so:/usr/lib/x86_64-linux-gnu/libmkl_intel_thread.so:/usr/lib/x86_64-linux-gnu/libiomp5.so

# echo "#=========== Compile Fortran"
# cd ../../build
# pwd
# rm   libFireCore.so
# make libFireCore
# cd $wd

#/usr/lib/gcc/x86_64-linux-gnu/11/libgfortran.so

echo "#=========== RUN"
#python3 relax_molecules.py
#python3 density_along_line.py
python3 export_HS_sparse.py | tee OUT