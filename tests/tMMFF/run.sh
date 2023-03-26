#!/bin/bash

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
pwd
rm libMMFF_lib.so
make MMFF_lib
rm   libLattice2D_lib.so
make Lattice2D_lib
cd $wd

cd ../../cpp/Build/libs_SDL
rm libMolGUIlib.so
make MolGUIlib
cd $wd

#rm *.bin

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

#python3 run.py
#python3 run_gui.py
#python3 run_surf_lattice.py
python3 run_propandiol.py
#python3 run_sample.py