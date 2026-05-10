#!/bin/bash
# Run 2D Ewald surface potential test
# Usage: ./run.sh [xyz_file] [n_harm]

N_HARM=${2:-3}

XYZ_NACL="../../tests/tQMMM_diacetylene/inputs/xyz/NaCl_1x1_L2.xyz"
XYZ_CAF2="../../cpp/common_resources/Substrates/generated_rect/CaF2_6L_Ni3_rect_nx1_nz1_fx0.666666_fy0.333333_L2_top.xyz"

XYZ=${1:-$XYZ_NACL}

echo "=== 2D Ewald Test: xyz=$XYZ n_harm=$N_HARM ==="

python3 test_ewald_2d.py --xyz $XYZ --n_harm $N_HARM --N_rep 30 --grid 200 --noPlot 2>&1 | tee run.log

echo ""
echo "=== Also running CaF2 ==="
python3 test_ewald_2d.py --xyz $XYZ_CAF2 --n_harm 4 --N_rep 20 --grid 150 --prefix CaF2 --noPlot 2>&1 | tee -a run.log
