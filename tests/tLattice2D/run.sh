#!/bin/bash

wd=`pwd`

#echo "#=========== Compile C++"
echo $wd
cd ../../cpp/Build/libs/Molecular/
pwd
rm   libLattice2D_lib.so
make Lattice2D_lib
cd $wd

echo "#=========== RUN"
python3 run.py