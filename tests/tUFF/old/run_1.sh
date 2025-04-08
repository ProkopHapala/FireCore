#!/bin/bash

cd ../../cpp/Build/libs/Molecular/
make -j4 MMFF_lib
cd - 1> /dev/null

python3 run.py "./data_UFF/xyz/ethylene" | tee out

