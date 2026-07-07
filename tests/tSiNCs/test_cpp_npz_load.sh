#!/bin/bash
# Agent A — NPZ load verify for MolecularBrowser (headless by default; --gui opens browser)
set -euo pipefail

OPEN_GUI=0
for arg in "$@"; do
    case "$arg" in
        --gui|--open) OPEN_GUI=1 ;;
    esac
done

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIX_DIR="$REPO_ROOT/tests/tSiNCs/fixtures/npz_viewer"
BUILD_DIR="$REPO_ROOT/cpp/Build/apps/MolecularEditor"
BROWSER="$BUILD_DIR/MolecularBrowser"

echo "=== test_cpp_npz_load: bootstrap fixtures ==="
node "$FIX_DIR/bootstrap_fixtures.mjs"

echo "=== test_cpp_npz_load: build MolecularBrowser ==="
if [ ! -x "$BROWSER" ]; then
    mkdir -p "$REPO_ROOT/cpp/Build"
    (cd "$REPO_ROOT/cpp/Build" && cmake .. -DWITH_SDL=ON -DCMAKE_BUILD_TYPE=Release && make MolecularBrowser -j"$(nproc)")
fi
if [ ! -x "$BROWSER" ]; then
    echo "ERROR: MolecularBrowser not found at $BROWSER"
    exit 1
fi

if [ -x "$(command -v g++)" ]; then
    LD_PRELOAD="$(g++ -print-file-name=libasan.so 2>/dev/null || true)"
    if [ -n "$LD_PRELOAD" ] && [ -f "$LD_PRELOAD" ]; then export LD_PRELOAD; fi
fi

verify() {
    local f="$1"
    echo "--- verify-npz $f ---"
    "$BROWSER" --verify-npz "$f"
}

verify "$FIX_DIR/01_init.npz"
verify "$FIX_DIR/03_topology.npz"

echo "=== test_cpp_npz_load: PASS ==="

if [[ "$OPEN_GUI" -eq 1 ]]; then
    echo "=== launching MolecularBrowser GUI ==="
    exec "$(dirname "$0")/run_cpp_mol_browser.sh" "$FIX_DIR"
fi
