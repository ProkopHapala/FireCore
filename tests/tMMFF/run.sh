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
export PYOPENCL_CTX=0

#rm *.bin

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
LD_PRELOAD=$LD_PRELOAD  $(g++ -print-file-name=libfftw3.so)
echo   $LD_PRELOAD
export LD_PRELOAD
# --- ignore memory leaks in ASAM
#export LSAN_OPTIONS=detect_leaks=0

# python3 run.py 
# python3 generate_scans.py --scan-types all --compare --nscan 63 --span-min 2.6 --span-max 15.2
# python3 generate_scans.py --scan-types all --compare --output-dir PTCDA_data_trial --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 215 --span-min 2.6 --span-max 15.2
# python3 generate_scans.py --scan-types all --compare --output-dir PTCDA_data_trial --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 215 --span-min 2.6 --span-max 15.1
# Run each scan type separately to avoid MMFF reinitialization issues
# python3 generate_scans.py --scan-types total --output-dir PTCDA_data_trial --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types morse --output-dir PTCDA_data_trial --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types coulomb --compare --output-dir PTCDA_data_trial --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# For 1D scan in a custom direction
# python3 generate_scans.py --scan-types total --scan-mode 1d --scan-dir1 "1,0,0" --output-dir PTCDA_data_trial_2d --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types morse --scan-mode 1d --scan-dir1 "1,0,0" --output-dir PTCDA_data_trial_2d --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types coulomb --compare --scan-mode 1d --scan-dir1 "1,0,0" --output-dir PTCDA_data_trial_2d --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-mode 1d --scan-dir "0,0,1" --scan-types morse --output-dir PTCDA_data_trial --lammps-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9

# For 2D scan in a custom plane
# python generate_scans.py --scan-mode xy --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,4.0" --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
python3 generate_scans.py --scan-types total --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.0732885303996796 --nscan-2 41 --span-min-2 0 --span-max-2 4.0732885303996796 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
python3 generate_scans.py --scan-types morse --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.0732885303996796 --nscan-2 41 --span-min-2 0 --span-max-2 4.0732885303996796 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
python3 generate_scans.py --scan-types coulomb --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.0732885303996796 --nscan-2 41 --span-min-2 0 --span-max-2 4.0732885303996796 --molecule data/xyz/new_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
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
#python3 run_test_GridFF_ocl_new.py
#python3 run_test_Multipole.py

# python3 run_sample_surf.py

#python3 run_sample_tricubic.py
#python3 run_sample_nonBond.py
#python3 run_lat_scan.py

#python3 run_collision_damp.py
#python3 run_collision_damp_scan.py
#python3 run_test_clear.py

#python3 run_opt_poly.py BB.HNH-hh.NHO-hp
#python3 run_opt_poly.py BB.HNH-hp.OHO-h_1,BB.HNH-hh.NHO-hp

#python3 run_test_ewald.py