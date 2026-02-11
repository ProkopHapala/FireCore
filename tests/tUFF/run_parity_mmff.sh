#!/bin/bash
# Comprehensive MMFF PyOpenCL vs C++ OpenCL parity sweep with rebuild
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Rebuild C++ OpenCL library (same flow as run_multi.sh)
(
  cd "${ROOT_DIR}/cpp/Build/libs_OCL" || exit 1
  rm -f libMMFFmulti_lib.so
  make MMFFmulti_lib
  mkdir -p ../Build/libs_OCL
  ln -sf ../../Build-opt/libs_OCL/libMMFFmulti_lib.so ../Build/libs_OCL/libMMFFmulti_lib.so
)

mols=(
  "${ROOT_DIR}/cpp/common_resources/mol/methanol.mol2"
  "${ROOT_DIR}/cpp/common_resources/xyz/HCONH2.xyz"
  "${ROOT_DIR}/cpp/common_resources/xyz/uracil.xyz"
  "${ROOT_DIR}/cpp/common_resources/mol/xylitol.mol2"
  "${ROOT_DIR}/cpp/common_resources/xyz/guanine.xyz"
  "${ROOT_DIR}/cpp/common_resources/xyz/Si10_H.xyz"
)
names=(
  "methanol.mol2"
  "HCONH2.xyz"
  "uracil.xyz"
  "xylitol.mol2"
  "guanine.xyz"
  "Si10_H.xyz"
)
# MMFF presets: bonds-only, bonds+angles/pi, bonds+angles+pi+nonbond
presets=(
  "--angles 0 --pisigma 0 --pipii 0 --nonbond 0 --label bonds"
  "--angles 1 --pisigma 1 --pipii 1 --nonbond 0 --label bonds,angles,pi"
  "--angles 1 --pisigma 1 --pipii 1 --nonbond 1 --label bonds,angles,pi+nonbond"
)

overall=0
cd "${SCRIPT_DIR}"
for idx in "${!mols[@]}"; do
  m="${mols[$idx]}"
  name="${names[$idx]}"
  for preset in "${presets[@]}"; do
    if python3 -u "${SCRIPT_DIR}/test_MMFF_cpp_vs_pyocl_parity.py" \
        --xyz "${m}" \
        --tolBuf 1e-6 \
        --tolF 5e-5 \
        --force-node-all 1 \
        --fast-exit 1 \
        ${preset} >/dev/null 2>&1; then
      status=PASS
    else
      status=FAIL
      overall=1
    fi
    label=${preset#*--label }
    echo "SUMMARY | ${name} | ${label} | ${status}"
  done
  echo "---------"
done

exit ${overall}
