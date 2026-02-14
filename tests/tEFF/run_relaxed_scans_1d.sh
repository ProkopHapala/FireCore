#!/usr/bin/env bash
# Batch runner for relaxed scan parity (CPU ialg=3 vs GPU localMD) plus plotting.
# Requires working directory tests/tEFF and prebuilt libeFF_lib.so (use ./run.sh if needed).
set -euo pipefail
cd "$(dirname "$0")"

# name  xyz_path                                        dt     damping  steps
SCANS=(
  "scan_parity_h2_distscan   export/scan_data/distscan_H2__spins_fc.xyz   0.01 0.1 2000"
  "scan_parity_ch4_distscan  export/scan_data/distscan_CH4__spins_fc.xyz  0.01 0.1 2000"
  "scan_parity_h2o_distscan  export/scan_data/distscan_H2O__spins_fc.xyz  0.01 0.1 2000"
)

for entry in "${SCANS[@]}"; do
  read -r NAME XYZ DT DAMP STEPS <<<"$entry"
  OUTDIR="export/${NAME}"
  echo "=== Running scan ${NAME} ==="
  mkdir -p "$OUTDIR"
  python3 -u run_relax_parity_protocol.py \
    --scan-xyz "$XYZ" \
    --fix-atoms \
    --dt "$DT" \
    --damping "$DAMP" \
    --scan-steps "$STEPS" \
    --short-steps 0 \
    --long-steps 0 \
    --outdir "$OUTDIR"

  if [[ -f "$OUTDIR/Es5_cpu.npy" && -f "$OUTDIR/Es5_gpu.npy" ]]; then
    python3 -u plot_scan_parity.py \
      --cpu "$OUTDIR/Es5_cpu.npy" \
      --gpu "$OUTDIR/Es5_gpu.npy" \
      --out "$OUTDIR/scan_parity.png" \
      --title "${NAME} (dt=${DT}, steps=${STEPS})"
  else
    echo "[WARN] Missing Es5_cpu.npy or Es5_gpu.npy in ${OUTDIR}, skipping plot"
  fi
done

echo "All scans finished. Outputs in export/<name>/"
