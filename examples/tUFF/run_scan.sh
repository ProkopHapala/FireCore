#!/bin/bash

# This script runs a Python script that emulates scans of GridFF in real time.

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

# z-scan UFF
export Z_RANGE="0.0,10.0,0.12"
echo "Z_RANGE="$Z_RANGE
python3 -u scan_files/scan.py --preset grid-only --grid-ff --nconf 100 --scan_dim z --z_range "$Z_RANGE" --surf xyz/surfaces_for_throughput/NaCl_6x6_L3 --disable bonds --scan_pos "0.0,0.0" 2>&1 | tee OUT-UFF-z-scan

# xy-scan UFF
export XY_RANGE="0.0,10.0,1.0,0.0,10.0,1.0"
python3 -u scan_files/scan.py --preset grid-only --grid-ff --use-scan-relaxed 0 --nconf 100 --scan_dim xy --xy_range "$XY_RANGE" --surf xyz/surfaces_for_throughput/NaCl_6x6_L3 --disable bonds --scan_pos "2.0" 2>&1 | tee OUT-UFF-xy-scan


# z-scan MMFF
export Z_RANGE="0.0,10.0,0.12"
echo "Z_RANGE="$Z_RANGE
python3 -u scan_files/scan_MMFF.py --preset grid-only --grid-ff --nsys 50 --scan_dim z --z_range "$Z_RANGE" --surf xyz/surfaces_for_throughput/NaCl_6x6_L3 --disable bonds --scan_pos "0.0,0.0" 2>&1 | tee OUT-UFF-z-scan

# xy-scan MMFF
export XY_RANGE="0.0,10.0,1.0,0.0,10.0,1.0"
python3 -u scan_files/scan_MMFF.py --preset grid-only --grid-ff --nconf 100 --scan_dim xy --xy_range "$XY_RANGE" --surf xyz/surfaces_for_throughput/NaCl_6x6_L3 --disable bonds --scan_pos "2.0" 2>&1 | tee OUT-UFF-xy-scan
