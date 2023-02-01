#!/bin/bash

wd=`pwd`

ln -s ../../cpp/sketches_SDL/Molecular/data

echo "#=========== Compile C++"
cd ../../cpp/Build/libs/Molecular/
pwd
rm libeFF_lib.so
make eFF_lib
cd $wd

echo "#=========== RUN"
#python3 run_tests.py
python3 run_dynamics.py 

#echo "#=========== compare_components.py"
#python3  run_evalPieces.py | tee output.txt
#python3 compare_components.py