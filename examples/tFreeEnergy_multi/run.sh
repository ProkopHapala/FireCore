#!/bin/bash

# Create symbolic links if they don't exist
if [ ! -d data ]; then
    ln -s ../../cpp/common_resources data
fi

# Get working directory
wd=`pwd`

# Build the library
echo "Building libMMFFmulti_lib.so..."
cd ../../cpp/Build/libs_OCL/
rm -f libMMFFmulti_lib.so
make MMFFmulti_lib
cd $wd

echo "Running run.py with parameters..."
python3 run.py \
    --nSys 10 \
    --xyz_name "data/DA.mol2" \
    --nCV 1.0 \
    --nLambda 10 \
    --nbStep 100 \
    --nMDsteps 100000 \
    --nEQsteps 10000 \
    --dt 0.5 \
    --tdamp 100.0 \
    --T 300.0
