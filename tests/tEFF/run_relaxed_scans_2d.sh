#!/usr/bin/env bash
# Batch relaxed 2D scans (frozen-core, spins/pairs) with CPU-vs-GPU parity and plotting.
# Runs relax (ialg=3 vs localMD) via run_relax_parity_protocol.py, then plots using
# eval_plot_cpu_gpu_maps.py consuming the precomputed Es5_* arrays.
# Run from tests/tEFF.
set -euo pipefail
cd "$(dirname "$0")"

# Parameters
DT=${DT:-0.01}
DAMP=${DAMP:-0.1}
STEPS=${STEPS:-2000}
NLOC=${NLOC:-32}
DEVICE=${DEVICE:-0}

# Discover 2D scan inputs (frozen-core spins/pairs)
inputs=()
for f in export/scan_data/*scan*__spins_fc.xyz export/scan_data/*scan*__pairs_fc.xyz; do
  [[ -f "$f" ]] || continue
  inputs+=("$f")
done

if [[ ${#inputs[@]} -eq 0 ]]; then
  echo "No scan inputs found matching spins/pairs fc patterns" >&2
  exit 1
fi

echo "Found ${#inputs[@]} scan inputs"

for xyz in "${inputs[@]}"; do
  base=$(basename "$xyz")
  tag=${base%.xyz}
  outdir="export/relaxed_${tag}"
  echo "=== Running relaxed scan ${tag} ==="
  mkdir -p "$outdir"

  if ! python3 -u run_relax_parity_protocol.py \
    --scan-xyz "$xyz" \
    --fix-atoms \
    --dt "$DT" \
    --damping "$DAMP" \
    --scan-steps "$STEPS" \
    --short-steps 0 \
    --long-steps 0 \
    --outdir "$outdir"; then
    echo "[WARN] relax failed for $xyz, skipping" >&2
    continue
  fi

  if [[ -f "$outdir/Es5_cpu.npy" && -f "$outdir/Es5_gpu.npy" ]]; then
    python3 -u eval_plot_cpu_gpu_maps.py \
      "$xyz" \
      --from-es5-cpu "$outdir/Es5_cpu.npy" \
      --from-es5-gpu "$outdir/Es5_gpu.npy" \
      --png \
      --noshow \
      --outdir "$outdir" \
      --nloc "$NLOC" \
      --device "$DEVICE"
  else
    echo "[WARN] Missing Es5_cpu.npy or Es5_gpu.npy in ${outdir}, skipping plot"
  fi

  echo "Done ${tag} -> ${outdir}"
done

echo "All relaxed 2D scans finished. Outputs in export/relaxed_*"
