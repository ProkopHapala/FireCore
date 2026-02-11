#!/usr/bin/env bash
set -euo pipefail

# UFF relaxation runner over reference molecules
# Usage (env overrides):
#   STEPS=500 DT=0.001 DAMP=0.05 DMAX=0.2 WRITE_STRIDE=5 FCONV=1e-3 VERBOSITY=1 QUIET=0 \
#   ./run_relax_all_uff.sh
# Optional: OUT_DIR to redirect outputs (defaults to CWD)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON=${PYTHON:-python3}
OUT_DIR=${OUT_DIR:-${PWD}}

# Tunables
STEPS=${STEPS:-2000}
DT=${DT:-0.01}
DAMP=${DAMP:-0.02}
DMAX=${DMAX:-0.3}
WRITE_STRIDE=${WRITE_STRIDE:-5}
FCONV=${FCONV:-1e-3}
FMAX_ABORT=${FMAX_ABORT:-1e3}
RMAX_ABORT=${RMAX_ABORT:-1e3}
VERBOSITY=${VERBOSITY:-1}
QUIET=${QUIET:-0}
SEED=${SEED:-0}

molecules=(
  "${SCRIPT_DIR}/../../cpp/common_resources/mol/methanol.mol2"
  "${SCRIPT_DIR}/../../cpp/common_resources/xyz/HCONH2.xyz"
  "${SCRIPT_DIR}/../../cpp/common_resources/xyz/uracil.xyz"
  "${SCRIPT_DIR}/../../cpp/common_resources/mol/xylitol.mol2"
  "${SCRIPT_DIR}/../../cpp/common_resources/xyz/guanine.xyz"
  "${SCRIPT_DIR}/../../cpp/common_resources/xyz/Si10_H.xyz"
)

log_summary="${OUT_DIR}/relax_summary_uff.txt"
: >"${log_summary}"

for mol in "${molecules[@]}"; do
  base=$(basename "${mol}")
  out_xyz="${OUT_DIR}/traj_${base}.UFF.xyz"
  echo "=== UFF ${base} ==="
  "${PYTHON}" -u "${SCRIPT_DIR}/relax_UFF_pyocl.py" \
    --molecule "${mol}" \
    --steps "${STEPS}" \
    --dt "${DT}" \
    --damp "${DAMP}" \
    --dmax "${DMAX}" \
    --seed "${SEED}" \
    --write-stride "${WRITE_STRIDE}" \
    --fconv "${FCONV}" \
    --fmax-abort "${FMAX_ABORT}" \
    --rmax-abort "${RMAX_ABORT}" \
    --verbosity "${VERBOSITY}" \
    --quiet "${QUIET}" \
    --out "${out_xyz}" | tee -a "${log_summary}"
done

echo "Summary written to ${log_summary}"
