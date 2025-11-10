#!/bin/bash

# This script runs an interactive GUI for UFF scanning.

# Set up symbolic links to common resources.
ln -s ../../cpp/common_resources data
ln -s ../../cpp/common_resources common_resources

# Get the current working directory.
wd=`pwd`

# Build the required shared library.
cd ../../cpp/Build/libs_OCL/
pwd
rm libMMFFmulti_lib.so
make MMFFmulti_lib
cd $wd

# Set the number of threads for OpenMP.
ncpu=`nproc`
ncpu=$(($ncpu - 1)) # Leave one CPU core free for user interaction.
echo "Compile using ncpu="$ncpu
OMP_NUM_THREADS=$ncpu
export OMP_NUM_THREADS

# Display CPU information.
lscpu

# Run the GUI application with surface name as argument
if [ -z "$1" ]; then
    echo "Usage: $0 <surface_name>"
    echo "Example: $0 NaCl_6x6_L3"
    exit 1
fi

SURF_NAME="$1"
echo "Starting GUI with surface: $SURF_NAME"
python3 -u scan_files/scan_gui.py --surf "xyz/surfaces_for_throughput/$SURF_NAME"

