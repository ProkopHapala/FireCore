#!/bin/bash

rm answer.bas answer.xyz params.dat CHARGES *.out

fireball="../../build/fireball.x"
fireball_py="python ../../pyBall/FireCore.py"
fireball_cpp="../../cpp/Build/apps/MolecularEditor/FireCoreVisual"



MKL_PATH=$PATH:/home/prokop/SW/intel/mkl/lib/intel64
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_PATH
export $LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

rm answer.bas answer.xyz params.dat CHARGES *.out

$fireball    | tee relaxation.out
$fireball_py | tee relaxation.py.out
$fireball_cpp | tee relaxation.cpp.out
