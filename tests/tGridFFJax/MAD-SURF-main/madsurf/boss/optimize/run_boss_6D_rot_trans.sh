#!/bin/bash
#SBATCH --account=project_2008059
#SBATCH -p medium # test - for testing 1h ; medium - up to 20 nodes/36 hours ; large - 20-200 nodes/36 hours #
#SBATCH --time=30:00:00      # dd-hh:mm:ss
#SBATCH -J BOSS        # name 
#SBATCH -o sbatch-%j.out # where the outputs & errors are written
#SBATCH -N 1                  # N nodes to run (N x 64 = n); max 192 ; max debug 48
#SBATCH --ntasks-per-node=128                # n processes to run (N x 64 = n); max 192 ; max debug 48
# check modules.txt for updating following lines: # 
#module purge
module load gcc/11.2.0 openmpi/4.1.2 openblas/0.3.18-omp netlib-scalapack/2.1.0 
export OMP_NUM_THREADS=1
ulimit -s unlimited
boss o in_6D_rot_trans