#!/bin/bash

name=FF2D_lib

#echo "#=========== Compile C++"
wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
#pwd
rm lib$name.so
make $name
cd $wd

echo "#=========== RUN"
#python3 run.py
python3 ff2dGUI.py