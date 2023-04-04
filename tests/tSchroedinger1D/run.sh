#!/bin/bash

wd=`pwd`

echo "#=========== Compile C++"
cd ../../cpp/Build/libs/Molecular/
pwd
rm   libSchroedingerGreen1D_lib.so
make SchroedingerGreen1D_lib
cd $wd

echo "#=========== RUN"
python3 run.py