#!/bin/bash
#SBATCH --account=project_2008059
#SBATCH -p test 	       # test - for testing 1h ; medium - up to 20 nodes/36 hours ; large - 20-200 nodes/36 hours #
#SBATCH --time=0-01:00:00      # dd-hh:mm:ss
#SBATCH -J BOSS_pp 	       # name 
#SBATCH -o sbatch-%j.out       # where the outputs & errors are written
#SBATCH -N 1                   # N nodes to run (N x 64 = n); max 192 ; max debug 48
#SBATCH --ntasks-per-node=128  # n processes to run (N x 64 = n); max 192 ; max debug 48
# #SBATCH --array=1,2,3,4,...   # uncomment and adjust if you want to run multiple instances of this, for instance for multiple systems. 

module load gcc/11.2.0 openmpi/4.1.2 openblas/0.3.18-omp netlib-scalapack/2.1.0 
export OMP_NUM_THREADS=1
ulimit -s unlimited

. /scratch/project_2008059/bossEnv_activate.sh

## Run script from folder containing BOSS run i.e. where you have the folders energy_calc, FHIaims_opt etc. ##

# Path for post-processing scripts
PP_PATH="/scratch/project_2008059/BOSS/pp_scripts"

mkdir -p results FHIaims_opt
cp boss.* results/
cd results

# Check if the .dat file already exists
if [ "$(find ../FHIaims_opt/locmin -maxdepth 1 -type f -name '*.dat' -print -quit)" ]; then
    echo "Local minima .dat file already exists. Skipping post-processing and creation of local minima files."
    cd ../FHIaims_opt
    mkdir -p locmin
    cd locmin
else
    # Running post_processing to get local minima. Using z variable as default
    bash $PP_PATH/pprepper_6D_unchanged_postprocessing_locmin.sh z
    echo "Mining local minima:."
    boss p boss_mod_z_locmin.rst boss.out 
    echo "Done mining local minima!"
    
    cd ../FHIaims_opt
    mkdir -p locmin
    cd locmin
    
    # Creating input files for the relaxation of local minima from the BOSS surrogate model
    echo "Creating input files for relaxation from BOSS predictions:"
    bash $PP_PATH/locMinInputGen_1.2_Mahti.sh ../../energy_calc/adsorbate_opt.in ../../energy_calc/substrate_opt.in ../../results/postprocessing/data_local_minima/*.dat
    echo "Input files created!"
fi

# Remove input duplicate structures
echo "Removing input duplicates:"
python $PP_PATH/removeDFTInputDuplicates_1.1_nequip_kabsch_PBC.py
echo "Input duplicates removed!"

# Remove input equivalent structures
echo "Removing input equivalents:"
python $PP_PATH/removeDFTInputEquivalents_1.1.py
echo "Input duplicates removed!"

# Relax the local minima
echo "Relaxing all predictions:"
bash $PP_PATH/DFT_locmin_relax_all_Mahti_timed.sh
echo "Done relaxing!"

## Removing duplicates/equivalent structures and analysing the data: 

##I have not installed povray on Mahti so this following part will not work currently

#echo "Removing duplicates/equivalents from relaxed data:"
#bash $PP_PATH/RunPostProcForDFT_conformers_with_povray_Mahti.sh

echo "Done with all post-processing!"


