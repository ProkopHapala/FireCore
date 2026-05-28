#!/bin/bash
# ORCA MPI run script with OpenMPI 4.1.8 compatibility fixes
# Usage: ./run_orca_mpi.sh <input.inp>

set -e

# 1. Purge path pollution from default system MPI
export PATH="/home/prokop/sw/openmpi-418/bin:$PATH"
export LD_LIBRARY_PATH="/home/prokop/sw/openmpi-418/lib:$LD_LIBRARY_PATH"

# 2. Add ORCA binaries to path
export PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$PATH"
export LD_LIBRARY_PATH="/home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418:$LD_LIBRARY_PATH"

# 3. Explicitly tell OpenMPI to use shared memory loops for a single desktop workstation
# This prevents it from looking for network cluster hardware it won't find
export OMPI_MCA_btl="vader,self"
export OMPI_MCA_pml="ob1"

# Input file
INPUT="${1:-adamantane_mpi.inp}"
INPUT_DIR=$(dirname "$INPUT")
INPUT_BASE=$(basename "$INPUT" .inp)

echo "========================================"
echo "ORCA MPI run script (OpenMPI 4.1.8)"
echo "========================================"
echo "Input: $INPUT"
echo "OMPI_MCA_btl: $OMPI_MCA_btl"
echo "OMPI_MCA_pml: $OMPI_MCA_pml"
echo "========================================"

cd "$INPUT_DIR"

# 4. Execute using the full absolute path of your compiled mpirun
# Using 12 cores with safe memory settings
/home/prokop/sw/openmpi-418/bin/mpirun --use-hwthread-cpus -np 12 /home/prokop/SW/orca_6_1_1_linux_x86-64_shared_openmpi418/orca "$INPUT_BASE.inp"

echo "========================================"
echo "ORCA MPI completed"
echo "========================================"
