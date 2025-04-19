#!/bin/bash

wd=`pwd`

ln -s ../../cpp/sketches_SDL/Molecular/data
ln -s ../../cpp/common_resources

echo "#=========== Compile C++"
cd ../../cpp/Build/libs/Molecular/
pwd
rm libeFF_lib.so
make eFF_lib
cd $wd

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

echo "#=========== RUN"

#python3 plot_EA.py 2> ERR | tee OUT
#python3 plot_EE.py 2> ERR | tee OUT

#python3 plot_EA.py 
#python3 plot_EE.py 
python3 -u run_scan_constr.py 2>ERR | tee OUT
#python3 run_energyToBondlength_Gabriel.py 2>ERR | tee OUT

#python3 -u run_process_xyz_1d.py 2>ERR | tee OUT
#python3 -u run_process_xyz.py 2>ERR | tee OUT
#python3 -u run_scan_Oe_ECP.py 2>ERR | tee OUT

#python3 run_tests.py 2> ERR
#python3 run_dynamics.py 

#echo "#=========== compare_components.py"
#python3  run_evalPieces.py | tee output.txt
#python3 compare_components.py