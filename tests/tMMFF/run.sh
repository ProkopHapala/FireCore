#!/bin/bash

ln -s ../../cpp/common_resources data 2>/dev/null
ln -s ../../cpp/common_resources common_resources 2>/dev/null

wd=`pwd`
cd ../../cpp/Build-opt/libs/Molecular/
pwd
rm libMMFF_lib.so
make MMFF_lib
rm   libLattice2D_lib.so
make Lattice2D_lib
cd $wd

cd ../../cpp/Build-opt/libs_SDL
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
LD_PRELOAD=$LD_PRELOAD:$(g++ -print-file-name=libfftw3.so)
echo   $LD_PRELOAD
export LD_PRELOAD
# --- ignore memory leaks in ASAM
export LSAN_OPTIONS=detect_leaks=0

#python3 run.py
#python3 run_hessian.py

#python3 run_gui.py
#python3 run_surf_lattice.py
#python3 run_propandiol.py
#python3 run_sample.py
#python3 run_Hbonds.py
#python3 run_sample_func.py
#python3 run_sample_Bsplines.py
#python3 run_sample_Hermite.py
#python3 run_test_ewald.py
#python3 run_test_GridFF.py
#python3 run_test_GridFF_ocl.py
#python3 run_test_Multipole.py

#python3 run_sample_surf.py

#python3 run_tipSpline_scan.py

#python3 run_tipSpline_scan.py --optimize 1 --nconf 100 --opt-attempts 1000 --opt-outdir opt_3d_target
#python3 run_tipSpline_scan.py --optimize 1 --nconf 100 --opt-attempts 100 --opt-outdir opt_3d_target

#echo "=== Running Vibration Spectra Test ==="
#python3 test_vibration_spectra.py

echo "=== Running Diamond Phonon Bands Test ==="
PY=python3
if [ -x /home/prokop/venvs/ML/bin/python3 ]; then PY=/home/prokop/venvs/ML/bin/python3; fi
if [ -n "$1" ]; then
  $PY "$@"
else
  $PY test_diamond_phonon_bands.py --unit THz --asr
fi

# Test without substrate
#echo "=== Testing without substrate ==="
#python3 run_tipSpline_scan_no_substrate.py

#python3 run_relax_surf.py

#python3 run_sample_tricubic.py
#python3 run_sample_nonBond.py
#python3 run_lat_scan.py

#python3 run_collision_damp.py
#python3 run_collision_damp_scan.py
#python3 run_test_clear.py

#python3 run_opt_poly.py BB.HNH-hh.NHO-hp
#python3 run_opt_poly.py BB.HNH-hp.OHO-h_1,BB.HNH-hh.NHO-hp

