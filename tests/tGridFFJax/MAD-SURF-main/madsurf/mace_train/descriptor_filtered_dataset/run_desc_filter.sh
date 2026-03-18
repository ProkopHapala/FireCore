#!/bin/bash
#SBATCH --account=XXX # your project number
#SBATCH --partition=YY # Your partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=01:15:00
#SBATCH --gres=gpu:a100:1

export PATH="PATH_TO_MACE_ENV:$PATH"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS

echo "GPU is available/Torch version:"
python3 -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'
#python3 -c 'import faiss; print(faiss.__version__)'

# Parameters are described in the python script # 
python descriptor_filter_batch_faiss.py \
  --input PATH_TO_DATASET_TO_FILTER/full_train_test_std_config_types.extxyz \
  --model PATH_TO_MACE_MODEL_USED_TO_CALC_DESCRIPTORS/MACE_model.model \
  --output filtered_structures.extxyz \
  --threshold 0.01 \ 
  --device cuda \
  --nn_index 320 \
  --nn_k 5 \
  --plot_distances \
  --hist_output my_histogram.png

python save_extxyz_frames_train_test_split.py filtered_structures.extxyz

