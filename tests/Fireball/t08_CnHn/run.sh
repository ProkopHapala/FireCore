#!/bin/bash

fireball="../../build/fireball.x"

#PATH=$PATH:/home/prokop/SW/intel/mkl/lib/intel64
#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/prokop/SW/intel/mkl/lib/intel64
#LD_LIBRARY_PATH=/home/prokop/SW/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
LD_LIBRARY_PATH=/home/prokop/intel/mkl/lib/intel64:$LD_LIBRARY_PATH
export $LD_LIBRARY_PATH
echo $LD_LIBRARY_PATH

$fireball | tee relaxation.out

image answer 3 2 2