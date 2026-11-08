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

# Do not LD_PRELOAD libasan unless MolecularBrowser was linked with it (NVIDIA GL + ASan preload aborts).

verify() {
    local f="$1"
    echo "--- verify-npz $f ---"
    "$BROWSER" --verify-npz "$f"
}

verify "$FIX_DIR/01_init.npz"
verify "$FIX_DIR/03_topology.npz"

echo "=== test_cpp_npz_load: Python np.savez ZIP64 extra (force_zip64=True) ==="
PY_NPZ="$(mktemp --suffix=.npz)"
TRUNC_NPZ="$(mktemp --suffix=.npz)"
python3 - "$PY_NPZ" <<'PY'
import sys
import numpy as np
path = sys.argv[1]
pos = np.array([[0.0, 0.0, 0.0], [1.5, 0.0, 0.0], [0.0, 1.5, 0.0]], dtype=np.float64)
Z = np.array([6, 6, 6], dtype=np.int32)
bonds = np.array([[0, 1], [0, 2]], dtype=np.int32)
np.savez(path, pos=pos, Z=Z, natoms=np.int32(3), bonds_ij=bonds)
print(f"wrote ZIP64 np.savez {path}")
PY
verify "$PY_NPZ"

echo "=== test_cpp_npz_load: JS reader must parse the same ZIP64 npz ==="
node --input-type=module - "$REPO_ROOT" "$PY_NPZ" <<'JS'
import { pathToFileURL } from 'node:url';
import fs from 'node:fs';
const repo = process.argv[2];
const npzPath = process.argv[3];
const { readNpzFile } = await import(pathToFileURL(repo + '/web/common_js/npzIO.js').href);
const { arrays } = readNpzFile(fs, npzPath);
if (arrays.pos.length !== 9) throw new Error(`pos length ${arrays.pos.length}`);
if (arrays.Z.length !== 3) throw new Error(`Z length ${arrays.Z.length}`);
console.log('JS ZIP64 read OK keys=', Object.keys(arrays).join(','));
JS

ATLAS_NPZ="$REPO_ROOT/tests/tSiNCs/OUT_chem_atlas/atlas/L1_dft/cube_7ring/C/02_relaxed.npz"
if [[ -f "$ATLAS_NPZ" ]]; then
    echo "=== test_cpp_npz_load: atlas Python 02_relaxed.npz ==="
    verify "$ATLAS_NPZ"
fi

echo "=== test_cpp_npz_load: genuinely truncated npz must still fail loud ==="
python3 - "$PY_NPZ" "$TRUNC_NPZ" <<'PY'
import sys
src, dst = sys.argv[1], sys.argv[2]
data = open(src, 'rb').read()
open(dst, 'wb').write(data[:80])
print(f"wrote truncated {dst} bytes=80 of {len(data)}")
PY
set +e
"$BROWSER" --verify-npz "$TRUNC_NPZ"
trunc_rc=$?
set -e
rm -f "$PY_NPZ" "$TRUNC_NPZ"
if [[ "$trunc_rc" -eq 0 ]]; then
    echo "ERROR: truncated npz was accepted"
    exit 1
fi
echo "truncated npz failed loud (rc=$trunc_rc) as required"

echo "=== test_cpp_npz_load: PASS ==="

if [[ "$OPEN_GUI" -eq 1 ]]; then
    echo "=== launching MolecularBrowser GUI ==="
    exec "$(dirname "$0")/run_cpp_mol_browser.sh" "$FIX_DIR"
fi
