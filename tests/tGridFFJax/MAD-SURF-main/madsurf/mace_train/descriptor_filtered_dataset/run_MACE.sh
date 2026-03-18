#!/bin/bash -l
#SBATCH --account=XXX # Your project
#SBATCH --partition=YYY # your partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=36:00:00
#SBATCH --gres=gpu:a100:1

export PATH="PATH_TO_MACE_ENV:$PATH"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS

mace_run_train \
--name="MACE_model" \
--train_file="./train_set.extxyz" \
--valid_fraction=0.05 \
--test_file="./test_set.extxyz" \
--num_workers=32 \
--atomic_numbers="[1, 6, 7, 8, 16, 17, 29, 35, 47, 79]" \
--E0s='{1:-13.6053084484, 6:-1029.1044589780, 7:-1485.3138688893, 8:-2043.2261372665, 16:-10875.6934658087, 17:-12577.6050159480, 29:-45242.4130555924, 35:-71401.3775787682, 47:-146385.3119095800, 79:-535650.8971237560}' \
--model="MACE" \
--energy_key="REF_energy" \
--forces_key="REF_forces" \
--hidden_irreps="128x0e + 128x1o" \
--r_max=5.0 \
--forces_weight=10 \
--energy_weight=1 \
--batch_size=64 \
--valid_batch_size=64 \
--max_num_epochs=206 \
--swa \
--start_swa=160 \
--ema \
--ema_decay=0.99 \
--amsgrad \
--restart_latest \
--error_table='PerAtomRMSE' \
--device=cuda \
--enable_cueq=True \
--seed=123 \

