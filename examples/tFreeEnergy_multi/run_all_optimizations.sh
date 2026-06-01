#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure required symlinks exist (OpenCL kernels need common_resources/cl)
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/data"              2>/dev/null || true
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/common_resources"  2>/dev/null || true

ALL_CONFIGS=(
  configs/optimize_parameters/opt_all.json
)
TOTAL=${#ALL_CONFIGS[@]}
FAILED=()

for i in "${!ALL_CONFIGS[@]}"; do
  cfg="${ALL_CONFIGS[$i]}"
  n=$((i + 1))
  echo ""
  echo "============================================================"
  echo "  [$n/$TOTAL] Running: $cfg"
  echo "============================================================"
  EXTRA_ARGS=(--reset)
  if bash "$SCRIPT_DIR/run_ES_params_optimization.sh" "$cfg" "${EXTRA_ARGS[@]}"; then
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
