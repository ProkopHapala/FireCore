#!/bin/bash

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/

pwd
rm libMMFF_lib.so
make MMFF_lib
cd $wd

python3 run.py