#!/bin/bash

wd=`pwd`
cd cpp/Build
make FireCoreVisual
cd $wd

MKL_PATH=$PATH:/home/prokop/SW/intel/mkl/lib/intel64
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_PATH
export $LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

readelf -s  /home/prokop/git/FireCore/build/libFireCore.so | grep FUNC | grep init

./cpp/Build/apps/MolecularEditor/FireCoreVisual