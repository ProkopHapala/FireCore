#!/usr/bin/env python
# coding: utf-8

# In[10]:


import os


# In[21]:


def generate_submission_script(run_name, lambda_value, output_dir=None ):
    """
    Create a MACE training script on a lambda value for energy/force weighting.
    """
    

    # Modify loss_coeffs based on lambda_value
    if lambda_value == 0 or lambda_value == "only_energy":
        force_coeff = 0
        energy_coeff = 1
        
    elif lambda_value == "infinity" or lambda_value == "only_forces":
        force_coeff = 1
        energy_coeff = 0
    else:
        force_coeff = lambda_value
        energy_coeff = 1

    run_template = f"""#!/bin/bash -l
#SBATCH --account=XXX
#SBATCH --partition=YYY
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --time=24:15:00
#SBATCH --gres=gpu:a100:1

export PATH="PATH_TO_MACE_ENVIRONMENT/bin:$PATH"

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo $OMP_NUM_THREADS

mace_run_train \\
--name="MACE_model" \\
--train_file="train_small_dataset_std.extxyz" \\
--valid_fraction=0.05 \\
--test_file="test_small_dataset_std.extxyz" \\
--config_type_weights='{"Default":1.0}' \\
--E0s='{{1:-13.6053084484, 6:-1029.1044589780, 7:-1485.3138688893, 8:-2043.2261372665, 16:-10875.6934658087, 17:-12577.6050159480, 29:-45242.4130555924, 35:-71401.3775787682, 47:-146385.3119095800, 79:-535650.8971237560}}' \\
--model="MACE" \\
--energy_key="REF_energy" \\
--forces_key="REF_forces" \\
--hidden_irreps="128x0e + 128x1o" \\
--r_max=5.0 \\
--forces_weight={force_coeff} \\
--energy_weight={energy_coeff} \\
--batch_size=10 \\
--max_num_epochs=1500 \\
--swa \\
--start_swa=1200 \\
--ema \\
--ema_decay=0.99 \\
--amsgrad \\
--restart_latest \\
--device=cuda \\
--enable_cueq=True \\
--seed=123 \\
"""

    # Ensure output directory exists
    if output_dir is None:
        output_dir = '.'
   
    folder_path =  os.path.abspath(os.path.join(output_dir, run_name))
    os.makedirs(folder_path, exist_ok=True)

    
    submit_file_path = os.path.join(folder_path, "run_MACE.sh")
    # Write the config to file
    with open(submit_file_path, "w") as f:
        f.write(run_template)

    print(f"Saved submission file to: {submit_file_path}")


# In[22]:


# Dictionary with {run_name: lambda}
runs = {
    "lambda_1": 1,
    "lambda_10": 10,
    "lambda_1000": 1000,
    "only_energy": 0,
    "only_forces": "infinity"
}

current_path = '/adsorbgmlip/mace_train/small_dataset'
# dataset paths 
train_path = 'PATH_TO_TRAINING_SET/train_small_dataset_std.extxyz'
test_path = 'PATH_TO_TEST_SET/test_small_dataset_std.extxyz'

for run_name, lamb in runs.items():
    
    generate_submission_script(run_name, lamb, output_dir=current_path )
    folder_path = os.path.join(current_path, run_name,)
    try:
        os.symlink(train_path, os.path.join(folder_path, 'train_small_dataset_std.extxyz'))
        os.symlink(test_path, os.path.join(folder_path, 'test_small_dataset_std.extxyz'))
    except:
        print('Linksym in folder_path already done')
    


# In[1]:


#get_ipython().system('pwd')

