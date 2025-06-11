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
# # For 1D scan in a custom direction
# python3 generate_scans.py --scan-types total --scan-mode 1d --scan-dir1 "0,0,1" --output-dir PTCDA_data_trial_1d --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --scan-origin "0,0,0" --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types morse --scan-mode 1d --scan-dir1 "0,0,1" --output-dir PTCDA_data_trial_1d --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --scan-origin "0,0,0" --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types coulomb --compare --scan-mode 1d --scan-dir1 "0,0,1" --output-dir PTCDA_data_trial_1d --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/1-rigid_zscan --scan-origin "0,0,0" --nscan 125 --span-min 2.6 --span-max 15.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9

# # For 2D scan in a custom plane
# python3 generate_scans.py --scan-types total --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial_2d --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.1 --nscan-2 41 --span-min-2 0 --span-max-2 4.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types morse --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial_2d --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.1 --nscan-2 41 --span-min-2 0 --span-max-2 4.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9
# python3 generate_scans.py --scan-types coulomb --compare --scan-mode 2d --scan-dir1 "1,0,0" --scan-dir2 "0,1,0" --scan-origin "0,0,3.3" --output-dir PTCDA_data_trial_2d --lammps-2d-dir /home/indranil/Documents/Project_1/Lammps/2-rigid_xyscan --nscan-1 41 --span-min-1 0 --span-max-1 4.1 --nscan-2 41 --span-min-2 0 --span-max-2 4.1 --molecule data/xyz/old_mol_old_sub_PTCDA --substrate data/xyz/Na_0.9_Cl_-0.9


# For 1D relaxed scan
# python3 generate_scans.py \
#   --scan-types total \
#   --scan-mode 1d \
#   --scan-dir1 "0.0,0.0,1.0" \
#   --scan-origin "0.0,0.0,0.0" \
#   --output-dir PTCDA_data_trial_1d_relax_z \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan \
#   --nscan 201 \
#   --span-min 1.3 \
#   --span-max 21.4 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26\

# python3 generate_scans.py \
#   --scan-types morse \
#   --scan-mode 1d \
#   --scan-dir1 "0.0,0.0,1.0" \
#   --scan-origin "0.0,0.0,0.0" \
#   --output-dir PTCDA_data_trial_1d_relax_z \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan \
#   --nscan 201 \
#   --span-min 1.3 \
#   --span-max 21.4 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26\

# python3 generate_scans.py \
#   --scan-types coulomb \
#   --scan-mode 1d \
#   --scan-dir1 "0.0,0.0,1.0" \
#   --scan-origin "0.0,0.0,0.0" \
#   --output-dir PTCDA_data_trial_1d_relax_z \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/3-relaxed_zscan \
#   --nscan 201 \
#   --span-min 1.3 \
#   --span-max 21.4 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 10000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26 \
#   --compare

## For 1D relaxed scan in angle
python3 generate_scans.py \
  --scan-types total \
  --scan-mode 1d \
  --scan-dir1 "0.894427190999916,0.447213595499958,0.0" \
  --scan-origin "0.0,0.0,2.8" \
  --output-dir PTCDA_data_trial_1d_relax_line \
  --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
  --nscan 357 \
  --span-min 0.0 \
  --span-max 35.7 \
  --molecule data/xyz/old_mol_old_sub_PTCDA \
  --substrate data/xyz/Na_0.9_Cl_-0.9 \
  --relaxed \
  --niter-max 100000 \
  --dt 0.02 \
  --fconv 1e-3 \
  --cons-atom 26
 

# python3 generate_scans.py \
#   --scan-types morse \
#   --scan-mode 1d \
#   --scan-dir1 "2.0,1.0,0.0" \
#   --scan-origin "0.0,0.0,2.8" \
#   --output-dir PTCDA_data_trial_1d_relax_line  \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
#   --nscan 357 \
#   --span-min 0.0 \
#   --span-max 35.7 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 100000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26\

# python3 generate_scans.py \
#   --scan-types coulomb \
#   --scan-mode 1d \
#   --scan-dir1 "2.0,1.0,0.0" \
#   --scan-origin "0.0,0.0,2.8" \
#   --output-dir PTCDA_data_trial_1d_relax_line  \
#   --lammps-1d-dir /home/indranil/Documents/Project_1/Lammps/4-relaxed_linescan \
#   --nscan 357 \
#   --span-min 0.0 \
#   --span-max 35.7 \
#   --molecule data/xyz/old_mol_old_sub_PTCDA \
#   --substrate data/xyz/Na_0.9_Cl_-0.9 \
#   --relaxed \
#   --niter-max 10000 \
#   --dt 0.02 \
#   --fconv 1e-3 \
#   --cons-atom 26 \
#   --compare
# #  
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
