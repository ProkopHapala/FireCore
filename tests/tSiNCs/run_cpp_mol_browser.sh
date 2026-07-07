#!/bin/bash
# Launch C++ MolecularBrowser on NPZ viewer fixtures (or any -dir).
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

if [[ -x "$(command -v g++)" ]]; then
  ASAN="$(g++ -print-file-name=libasan.so 2>/dev/null || true)"
  if [[ -n "${ASAN}" && -f "${ASAN}" ]]; then export LD_PRELOAD="${ASAN}"; fi
fi

DIR="${1:-${FIX}}"
echo "MolecularBrowser: -dir '${DIR}'  (contains npz_viewer/, si_1nm_passivation/, ...)"
cd "${BUILD_DIR}"
exec ./MolecularBrowser -res "${RES}" -dir "${DIR}"
