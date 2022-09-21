#!/bin/bash

wd=`pwd`
cd ../../cpp/Build_OCL/libs/Molecular/
pwd
rm libMMFF_lib.so
make MMFF_lib
cd $wd

python3 run.py