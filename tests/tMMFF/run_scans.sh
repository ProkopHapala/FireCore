#!/bin/bash

# Script to run all three scan types in sequence
# Each scan runs in its own process to avoid memory management issues

ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources 

wd=`pwd`
cd ../../cpp/Build/libs/Molecular/
pwd
rm libMMFF_lib.so
make MMFF_lib
rm   libLattice2D_lib.so
make Lattice2D_lib
cd $wd

cd ../../cpp/Build/libs_SDL
rm libMolGUIlib.so
make MolGUIlib
cd $wd


# Create symlinks if they don't exist
if [ ! -e data ]; then
    ln -s ../../cpp/common_resources data
fi

if [ ! -e common_resources ]; then
    ln -s ../../cpp/common_resources common_resources
fi

# Set up multiprocessing
ncpu=$(nproc)
ncpu=$((ncpu - 1))     # leave one CPU free for user interaction
echo "Using ncpu=$ncpu"
export OMP_NUM_THREADS=$ncpu
export PYOPENCL_CTX=0

# Set up ASan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
LD_PRELOAD="$LD_PRELOAD $(g++ -print-file-name=libfftw3.so)"
echo "LD_PRELOAD=$LD_PRELOAD"
export LD_PRELOAD

# Default parameters
MOLECULE="data/xyz/old_mol_old_sub_PTCDA"
SUBSTRATE="data/xyz/Na_0.9_Cl_-0.9"
OUTPUT_DIR="PTCDA_data"

# Parse command line arguments
if [ $# -ge 1 ]; then
    MOLECULE=$1
fi

if [ $# -ge 2 ]; then
    SUBSTRATE=$2
fi

if [ $# -ge 3 ]; then
    OUTPUT_DIR=$3
fi

# Create output directory
mkdir -p "$OUTPUT_DIR"

# Display parameters
echo "====================================================================="
echo "Running scans with parameters:"
echo "Molecule:       $MOLECULE"
echo "Substrate:      $SUBSTRATE"
echo "Output Directory: $OUTPUT_DIR"
echo "====================================================================="

# Run total potential scan
echo "\n====================================================================="
echo "STEP 1/3: RUNNING TOTAL POTENTIAL SCAN"
echo "====================================================================="
python3 scan_total.py "$MOLECULE" "$SUBSTRATE" "$OUTPUT_DIR"

# Check if the previous scan completed successfully
# if [ $? -ne 0 ]; then
#     echo "ERROR: Total potential scan failed!"
#     exit 1
# fi

# Run London-Pauli (Morse) scan
echo "\n====================================================================="
echo "STEP 2/3: RUNNING LONDON-PAULI (MORSE) SCAN"
echo "====================================================================="
python3 scan_morse.py "$MOLECULE" "$SUBSTRATE" "$OUTPUT_DIR"

# Check if the previous scan completed successfully
# if [ $? -ne 0 ]; then
#     echo "ERROR: London-Pauli (Morse) scan failed!"
#     exit 1
# fi

# Run Coulomb scan
echo "\n====================================================================="
echo "STEP 3/3: RUNNING COULOMB SCAN"
echo "====================================================================="
python3 scan_coulomb.py "$MOLECULE" "$SUBSTRATE" "$OUTPUT_DIR"

# Check if the previous scan completed successfully
# if [ $? -ne 0 ]; then
#     echo "ERROR: Coulomb scan failed!"
#     exit 1
# fi

echo "\n====================================================================="
echo "ALL SCANS COMPLETED SUCCESSFULLY"
echo "Data files saved in $OUTPUT_DIR"
echo "====================================================================="
