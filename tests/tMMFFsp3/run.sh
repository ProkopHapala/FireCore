#!/bin/bash

name=MMFFsp3_lib
ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
pwd
rm lib$name.so
make -j4 $name
cd $wd

#python3 run.py
#python3 HBscan.py

python3 run_Eprofiles.py