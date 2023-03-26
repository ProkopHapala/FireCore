#!/bin/bash

#LD_PRELOAD=path/to/asan/runtime/lib
#LD_PRELOAD=/usr/lib/gcc/x86_64-linux-gnu/11/libasan.so
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

python3 run_asan.py
