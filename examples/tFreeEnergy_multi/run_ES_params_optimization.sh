#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Ensure required symlinks exist (OpenCL kernels need common_resources/cl)
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/data"              2>/dev/null || true
ln -sf ../../cpp/common_resources "$SCRIPT_DIR/common_resources"  2>/dev/null || true

CONFIG_FILE="$SCRIPT_DIR/es_optimization_config.json"
RESET=0
EXTRA_ARGS=()

for arg in "$@"; do
    case "$arg" in
        --reset)
            RESET=1
            ;;
        *.json)
            CONFIG_FILE="$arg"
            ;;
        *)
            EXTRA_ARGS+=("$arg")
            ;;
    esac
done

if [[ ! -f "$CONFIG_FILE" ]]; then
    echo "ERROR: Config file not found: $CONFIG_FILE" >&2
    exit 1
fi

echo "Using config: $CONFIG_FILE"
OUTPUT_ROOT="$(python3 - "$CONFIG_FILE" <<'PY'
import json,sys
with open(sys.argv[1],"r",encoding="utf-8") as f:
    cfg=json.load(f)
print(cfg.get("output_root","bench_ES"))
PY
)"

if [[ "$RESET" -eq 1 ]]; then
    TARGET="$SCRIPT_DIR/$OUTPUT_ROOT"
    echo "Reset requested: removing $TARGET"
    rm -rf "$TARGET"
fi

python3 "$SCRIPT_DIR/run_ES_optimization.py" --config "$CONFIG_FILE" "${EXTRA_ARGS[@]}"
python3 "$SCRIPT_DIR/collect_ES_summary.py" --output-root "$OUTPUT_ROOT"
echo "Done. Output root: $SCRIPT_DIR/$OUTPUT_ROOT"
