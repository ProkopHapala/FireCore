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


# ---- Multiprocesing
ncpu=`nproc`
ncpu=$(($ncpu - 1))     # let one CPU free for user interaction
echo "compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

#rm *.bin

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

#python3 run.py
#python3 run_gui.py
#python3 run_surf_lattice.py
#python3 run_propandiol.py
#python3 run_sample.py
#python3 run_Hbonds.py
#python3 run_sample_func.py
#python3 run_sample_Bsplines.py
python3 run_sample_Hermite.py
#python3 run_sample_surf.py

#python3 run_sample_tricubic.py
#python3 run_sample_nonBond.py
#python3 run_lat_scan.py

#python3 run_collision_damp.py
#python3 run_collision_damp_scan.py
#python3 run_test_clear.py

#python3 run_opt_poly.py BB.HNH-hh.NHO-hp
#python3 run_opt_poly.py BB.HNH-hp.OHO-h_1,BB.HNH-hh.NHO-hp

