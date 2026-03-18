#!/bin/bash
#SBATCH --account=XXX
#SBATCH --partition=medium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=128
#SBATCH --time=01:00:00
 
export PATH=""

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS

mkdir processed_data

python PATH_TO_MACE_ENVIRONMENT/mace/cli/preprocess_data.py \
--train_file="PATH_TO_TRAINING_SET/train_dataset_std.extxyz" \
--valid_fraction=0.05 \
--test_file="PATH_TO_TEST_SET/test_dataset_std.extxyz" \
--atomic_numbers="[1, 6, 7, 8, 16, 17, 29, 35, 47, 79]" \
--r_max=5.0 \
--h5_prefix="processed_data/" \
--E0s='{1:-13.6053084484, 6:-1029.1044589780, 7:-1485.3138688893, 8:-2043.2261372665, 16:-10875.6934658087, 17:-12577.6050159480, 29:-45242.4130555924, 35:-71401.3775787682, 47:-146385.3119095800, 79:-535650.8971237560}' \
--seed=123 \

