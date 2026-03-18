#!/bin/bash
#SBATCH --partition=XXX
#SBATCH --account=YYY
#SBATCH --time=01:30:00
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --ntasks-per-node=1

export OMP_NUM_THREADS=32
echo $OMP_NUM_THREADS

export MPICH_GPU_SUPPORT_ENABLED=1

## path to MACE env (bin) ##
export PATH=":$PATH"

echo "GPU is available/Torch version:"
python3 -c 'import torch; print(torch.cuda.is_available()); print(torch.__version__)'

echo "MACE Version:"
python3 -c "import mace; print(mace.__version__)"

echo "Running annealing script..."
srun python mace_annealer_OH_network.py "$1" # input structure as command line input
echo "Done with annealing!"
