#!/bin/bash
# Launch C++ MolecularBrowser (SDL). From FireCore repo root:
#
#   ./tests/tSiNCs/run_cpp_mol_browser.sh
#   ./tests/tSiNCs/run_cpp_mol_browser.sh tests/tSiNCs/OUT_chem_atlas/atlas
#
# Keys: Esc = up (VIEW→folder, then parent). Ctrl+Q or Ctrl+D = quit.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
FIX="${ROOT}/tests/tSiNCs/fixtures"
NPZ_FIX="${FIX}/npz_viewer"
BUILD_DIR="${ROOT}/cpp/Build/apps/MolecularEditor"
BROWSER="${BUILD_DIR}/MolecularBrowser"
RES="${ROOT}/cpp/common_resources"

if [[ ! -f "${NPZ_FIX}/01_init.npz" ]]; then
  node "${NPZ_FIX}/bootstrap_fixtures.mjs"
fi

if [[ ! -x "${BROWSER}" ]]; then
  echo "Building MolecularBrowser..."
  mkdir -p "${ROOT}/cpp/Build"
  (cd "${ROOT}/cpp/Build" && cmake .. -DWITH_SDL=ON && make MolecularBrowser -j"$(nproc)")
fi

# Only preload ASan if this binary was actually linked with it.
# LD_PRELOAD libasan on a Release (WITH_ASAN=OFF) build + NVIDIA GL aborts with:
#   "Shadow memory range interleaves with an existing memory mapping"
if ldd "${BROWSER}" 2>/dev/null | grep -q libasan; then
  export LSAN_OPTIONS="${LSAN_OPTIONS:-detect_leaks=0}"
  echo "MolecularBrowser is ASan-linked; LSAN_OPTIONS=${LSAN_OPTIONS}"
fi

DIR="${1:-${FIX}}"
if [[ ! -d "${DIR}" ]]; then
  echo "ERROR: directory does not exist: ${DIR}  (cwd=$(pwd))" >&2
  exit 1
fi
DIR="$(cd "${DIR}" && pwd)"
echo "cmd: ${BROWSER} -res ${RES} -dir ${DIR}"
echo "keys: Esc=up   Ctrl+Q / Ctrl+D=quit"
exec "${BROWSER}" -res "${RES}" -dir "${DIR}"
