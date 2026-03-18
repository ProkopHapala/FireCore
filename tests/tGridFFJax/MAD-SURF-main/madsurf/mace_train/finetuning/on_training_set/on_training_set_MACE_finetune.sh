#!/bin/bash -l
### - HPC-specific parameters - ### 
#SBATCH --account=XXX
##SBATCH --partition=gpusmall 
##SBATCH --nodes=1
##SBATCH --ntasks-per-node=1
##SBATCH --cpus-per-task=32
##SBATCH --time=02:15:00
##SBATCH --gres=gpu:a100:1

#export PATH=""

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#echo $OMP_NUM_THREADS
### - ###

mace_run_train \
  --name MACE_finetuned \
  --train_file ./all_training_mol_rels_std_config_types.extxyz \
  --foundation_model ../mace-mpa-0-medium.model \
  --pt_train_file "mp" \
  --num_samples_pt 1800 \
  --filter_type_pt combinations \
  --subselect_pt fps \
  --weight_pt 1.0 \
  --foundation_filter_elements=True \
  --atomic_numbers="[1, 6, 7, 8, 16, 17, 29, 35, 47, 79]" \
  --E0s='{1:-13.6053084484, 6:-1029.1044589780, 7:-1485.3138688893, 8:-2043.2261372665, 16:-10875.6934658087, 17:-12577.6050159480, 29:-45242.4130555924, 35:-71401.3775787682, 47:-146385.3119095800, 79:-535650.8971237560}' \
  --multiheads_finetuning True \
  --force_mh_ft_lr False \
  --error_table='PerAtomRMSE' \
  --device=cuda \
  --enable_cueq=True \
  --forces_weight=10 \
  --energy_weight=1 \
  --stress_weight=0 \
  --max_num_epochs=40 \


