#!/bin/bash
# Comprehensive UFF PyOpenCL vs C++ OpenCL parity sweep with rebuild
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# Rebuild C++ OpenCL library (mirror run_multi.sh flow)
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
comps=("bonds" "bonds,angles" "bonds,angles,dihedrals" "bonds,angles,dihedrals,inversions")

overall=0
cd "${SCRIPT_DIR}"
for idx in "${!mols[@]}"; do
  m="${mols[$idx]}"
  name="${names[$idx]}"
  for c in "${comps[@]}"; do
    if python3 -u "${SCRIPT_DIR}/test_UFF_cppocl_vs_pyocl_parity.py" \
        --molecule "${m}" \
        --components "${c}" \
        --tolerance 5e-5 \
        --tol-buf-f 1e-6 \
        --check-torsion-inputs 0 \
        --fast-exit 1 >/dev/null 2>&1; then
      status=PASS
    else
      status=FAIL
      overall=1
    fi
    echo "SUMMARY | ${name} | ${c} | ${status}"
  done
  echo "---------"
done

exit ${overall}
