#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources

dir="../../cpp/Build/apps/MolecularEditor/"
name="FireCoreVisual"


#MKL_PATH=$PATH:/home/prokop/SW/intel/mkl/lib/intel64
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_PATH
#export $LD_LIBRARY_PATH
#echo $LD_LIBRARY_PATH

# ---- Compilation
wd=`pwd`
cd $dir
pwd
rm $name
make -j$ncpu $name   # 2>$wd/compile_err.log
cd $wd

#cd ../../build/
#make libFireCore
#cd $wd

ln -s $dir/$name .
ln -s ../../build/libFireCore.so .

#export $LD_LIBRARY_PATH:../../build/

rm answer.bas answer.xyz params.dat CHARGES *.out

./$name