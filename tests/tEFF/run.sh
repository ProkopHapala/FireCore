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
python3 run.py