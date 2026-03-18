#!/bin/bash
#SBATCH --account=project_2008059
#SBATCH -p medium # medium - for testing 1h ; medium - up to 20 nodes/36 hours ; large - 20-200 nodes/36 hours #
#SBATCH --time=03:00:00      # dd-hh:mm:ss
#SBATCH -J aimsTest-array        # name 
#SBATCH -o sbatch-%j.out # where the outputs & errors are written
#SBATCH -N 1                  # N nodes to run (N x 64 = n); max 192 ; max debug 48
#SBATCH --ntasks=128                # n processes to run (N x 64 = n); max 192 ; max debug 48

# check modules.txt for updating following lines: # module purge
module load gcc/11.2.0 openmpi/4.1.2 openblas/0.3.18-omp netlib-scalapack/2.1.0 
export OMP_NUM_THREADS=1
ulimit -s unlimited
export PATH="/scratch/project_2008059/dscribe/bin:$PATH"
#. /scratch/project_2008059/bossEnv_activate.sh

# Create necessary folders
mkdir -p mode_analysis
mkdir -p DFT_cal/selected_indices

# Copy the necessary files/scripts
cp /scratch/project_2008059/NMS/general_folder/NM_sampling.py .
cp /scratch/project_2008059/NMS/general_folder/main_vib.py mode_analysis/main_vib.py
cp /scratch/project_2008059/NMS/general_folder/DFT_cal_aims.py DFT_cal/selected_indices/DFT_cal_aims.py
cp geometry_opt.in mode_analysis/geometry_opt.in

echo "Performing vibrational calculations on system"
cd mode_analysis
python main_vib.py > output.txt
echo "done doing vibrational calcs on system"

echo "Running NMS"
cd ..
python NM_sampling.py
echo "done running NMS"

echo "Running DFT calculations" 
cd DFT_cal/selected_indices
python DFT_cal_aims.py > output_dft.txt
echo "done doing DFT calcs on system"

