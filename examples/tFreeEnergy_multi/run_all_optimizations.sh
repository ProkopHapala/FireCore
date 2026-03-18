#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure required symlinks exist (OpenCL kernels need common_resources/cl)
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/data"              2>/dev/null || true
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/common_resources"  2>/dev/null || true

# JE configs get --reset to clear old failed results from broken thermalization
TI_CONFIGS=(
  opt_nLambda_TI.json
  opt_nEQsteps_TI.json
  opt_nMDsteps_TI.json
  opt_dt_TI.json
  opt_t_damp_TI.json
)

JE_CONFIGS=(
  opt_nLambda_JE.json
  opt_nEQsteps_JE.json
  opt_nMDsteps_JE.json
  opt_dt_JE.json
  opt_t_damp_JE.json
  opt_K.json
  opt_nSys_JE.json
)

ALL_CONFIGS=("${TI_CONFIGS[@]}" "${JE_CONFIGS[@]}")
TOTAL=${#ALL_CONFIGS[@]}
FAILED=()

# Build a set of JE config names for quick lookup
declare -A IS_JE
for jc in "${JE_CONFIGS[@]}"; do
  IS_JE["$jc"]=1
done

for i in "${!ALL_CONFIGS[@]}"; do
  cfg="${ALL_CONFIGS[$i]}"
  n=$((i + 1))
  echo ""
  echo "============================================================"
  echo "  [$n/$TOTAL] Running: $cfg"
  echo "============================================================"
  EXTRA_ARGS=()
  if [[ -n "${IS_JE[$cfg]+x}" ]]; then
    EXTRA_ARGS+=(--reset)
    echo "  (JE config — resetting previous results)"
  fi
  if bash "$SCRIPT_DIR/run_ES_params_optimization.sh" "$SCRIPT_DIR/$cfg" "${EXTRA_ARGS[@]}"; then
    echo "  [$n/$TOTAL] $cfg finished OK"
  else
    echo "  [$n/$TOTAL] $cfg FAILED (continuing)"
    FAILED+=("$cfg")
  fi
done

echo ""
echo "============================================================"
echo "  All optimization sweeps complete."
echo "  Succeeded: $((TOTAL - ${#FAILED[@]}))/$TOTAL"
if [ ${#FAILED[@]} -gt 0 ]; then
  echo "  Failed:"
  for f in "${FAILED[@]}"; do
    echo "    - $f"
  done
fi
echo "============================================================"

