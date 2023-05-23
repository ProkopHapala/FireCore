#!/bin/bash

name=MMFFsp3_lib
ln -sf ../../cpp/common_resources data
ln -sf ../../cpp/common_resources common_resources 

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
pwd
rm lib$name.so
make -j4 $name
cd $wd

#python3 run.py
#python3 HBscan.py

python3 run_Eprofiles.py

#dir=/home/niko/work/CARBSIS/2-molecules/data_Paolo/uff_files/test_UFF/CHONH2

#mol='CHONH2'

#ds='-0.30 -0.25 -0.20 -0.15 -0.10 -0.05 0.00 0.05 0.10 0.15 0.20 0.25 0.30'

#rm -f x.xyz
#mkdir -p $mol
#for d in $ds
#do
#    echo $mol $d
#    ln -sf $dir/improper$d/run.xyz x.xyz
#    python3 eval_E_xyz.py x | grep MYOUTPUT 1> $mol/out.$d
#done

#rm x.xyz

exit
