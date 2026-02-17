#!/bin/bash
# run_fdbm.sh  –  FDBM AFM simulation of PTCDA using Fireball density
# Usage:
#   ./run_fdbm.sh             – default params
#   ./run_fdbm.sh --noplot    – skip plotting
#   ./run_fdbm.sh --nxy 80 80 – finer scan

cd "$(dirname "$0")"
ROOT="$(realpath ../../)"
export PYTHONPATH="$ROOT:$PYTHONPATH"
echo "=== run_fdbm.sh  ROOT=$ROOT ==="
python3 test_fdbm.py "$@" 2>&1 | tee run_fdbm.log
echo "=== run_fdbm.sh DONE ==="
