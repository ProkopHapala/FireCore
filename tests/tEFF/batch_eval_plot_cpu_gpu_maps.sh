#!/usr/bin/env bash
set -euo pipefail

# Run from tests/tEFF (important for relative kernel/resource paths)

OUTDIR=${OUTDIR:-export/plots_cpu_gpu}
NOSHOW=${NOSHOW:-1}
PNG=${PNG:-1}
REBUILD=${REBUILD:-0}                   # set REBUILD=1 to force one rebuild before batch
CPP_BUILD_PATH=${CPP_BUILD_PATH:-$(realpath ../../cpp/Build-opt/libs)}

# If libeFF_lib.so was linked with AddressSanitizer, Python will abort unless libasan
# is preloaded before the shared library is loaded.
LIBASAN=$(gcc -print-file-name=libasan.so || true)
if [[ -n "${LIBASAN}" && -f "${LIBASAN}" ]]; then
    export LD_PRELOAD="${LIBASAN}${LD_PRELOAD:+:$LD_PRELOAD}"
    export ASAN_OPTIONS=${ASAN_OPTIONS:-detect_leaks=0}
fi

# Optional one-time rebuild of eFF_lib
if [[ "$REBUILD" == "1" ]]; then
    echo "== Rebuilding eFF_lib in $(dirname "$CPP_BUILD_PATH")"
    ( cd "$(dirname "$CPP_BUILD_PATH")" && make -j"$(nproc)" eFF_lib )
fi

# Disable per-run rebuilds in Python; point to chosen build
export CPP_RECOMPILE=0
export CPP_BUILD_PATH

mkdir -p "$OUTDIR"

args=("--outdir" "$OUTDIR")
if [[ "$NOSHOW" == "1" ]]; then args+=("--noshow"); fi
if [[ "$PNG" == "1" ]]; then args+=("--png"); fi

for f in export/scan_data/*.xyz; do
    # skip non-scan xyz if any
    base=$(basename "$f")
    if [[ "$base" != *scan* ]]; then
        continue
    fi
    if ! grep -q "na,ne" "$f"; then
        echo "SKIP (no na,ne header): $f"
        continue
    fi
    echo "==== $f"
    PYTHONPATH=../.. python3 -u eval_plot_cpu_gpu_maps.py "$f" "${args[@]}" || {
        echo "FAILED $f" >&2
    }
done

echo "DONE. Outputs in $OUTDIR"
