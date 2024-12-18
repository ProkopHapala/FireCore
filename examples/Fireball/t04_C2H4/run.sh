#!/bin/bash

rm answer.bas answer.xyz params.dat CHARGES *.out
fireball="../../build/fireball.x"
#PATH=$PATH:/home/prokop/SW/intel/mkl/lib/intel64
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prokop/SW/intel/mkl/lib/intel64
#LD_LIBRARY_PATH=/home/prokop/SW/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/prokop/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
export $LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH


rm answer.bas answer.xyz params.dat CHARGES *.out


$fireball | tee relaxation.out