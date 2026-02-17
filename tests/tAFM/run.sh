#!/bin/bash
# run.sh  –  test AFM simulation of PTCDA
# Usage:
#   ./run.sh             – LJ + fixed charges
#   ./run.sh --morse     – Morse + fixed charges
#   ./run.sh --qeq       – LJ + QEq charges
#   ./run.sh --noplot    – skip plotting

cd "$(dirname "$0")"

ROOT="$(realpath ../../)"
export PYTHONPATH="$ROOT:$PYTHONPATH"

echo "=== tAFM run.sh  ROOT=$ROOT ==="
python3 test_ptcda.py "$@" 2>&1 | tee run.log
echo "=== run.sh DONE ==="
