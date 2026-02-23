#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# The script itself lives in tests/tFitREQ, so WORKDIR is already ROOT
WORKDIR="$ROOT"
OUT_DIR="$WORKDIR/grid_search_out"
XYZ="/home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz"
#XYZ="/home/niko/work/HBOND/REFERENCE/2-pairs_small_small/4-to_firecore/confs_wb97m/H2O-A1_H2O-D1-y.xyz"
DOF="dofSelection_MorseSR_H2O_epairOnly.dat"
NPZ="$OUT_DIR/ensemble_data.npz"
OUT_PREFIX="$OUT_DIR/ensemble"
# Endpoints for string scan (will be auto-detected if not set manually)
DOF_A=""
DOF_B=""
PATH_PREFIX="$OUT_DIR/path"

echo "[1/3] Running grid search -> ensemble and combined plots ..."
cd "$WORKDIR"
python grid_search.py --xyz "$XYZ" --dof_file "$DOF" --out_dir "$OUT_DIR"

echo "[2/3] Running DR (PCA/UMAP/PacMAP if available) on ensemble ..."
if [[ ! -f "$NPZ" ]]; then
  echo "ERROR: Ensemble NPZ not found: $NPZ"
  exit 1
fi
python plot_dr.py --trajectory_npz "$NPZ" --out_prefix "$OUT_PREFIX"

echo "[3/3] String scan + Hessian (NEB-like) ..."
# Auto-detect best and worst endpoints from grid_search if not set
if [[ -z "$DOF_A" || -z "$DOF_B" ]]; then
  echo "Auto-detecting endpoints from $NPZ ..."
  # Use a quick python one-liner to find min/max error indices
  IDX_A=$(python3 -c "import numpy as np; d=np.load('$NPZ'); print(np.argmin(d['error']))")
  IDX_B=$(python3 -c "import numpy as np; d=np.load('$NPZ'); print(np.argmax(d['error']))")
  
  # Find the corresponding run_*.dat files
  DOF_A=$(ls $OUT_DIR/run_${IDX_A}_*.dat 2>/dev/null | head -n 1 || true)
  DOF_B=$(ls $OUT_DIR/run_${IDX_B}_*.dat 2>/dev/null | head -n 1 || true)
  
  echo "Selected DOF_A (min error): $DOF_A"
  echo "Selected DOF_B (max error): $DOF_B"
fi

if [[ -n "$DOF_A" && -n "$DOF_B" && -f "$DOF_A" && -f "$DOF_B" ]]; then
  # RIGID SCAN
  #python string_scan.py --dofA "$DOF_A" --dofB "$DOF_B" --xyz "$XYZ" --n_points 11 --out_prefix "$PATH_PREFIX"
  # RELAXED SCAN & HESSIAN
  python string_scan.py --dofA "$DOF_A" --dofB "$DOF_B" --xyz "$XYZ" --n_points 11 --relax --hessian --out_prefix "$PATH_PREFIX"
  # DR on path trajectory (PCA/UMAP/PacMAP if available)
  if [[ -f "${PATH_PREFIX}_trajectory.npz" ]]; then
    python plot_dr.py --traj_file "${PATH_PREFIX}_trajectory.npz" --out_prefix "${PATH_PREFIX}_dr"
  else
    echo "WARNING: path trajectory not found at ${PATH_PREFIX}_trajectory.npz"
  fi
else
  echo "WARNING: DOF_A or DOF_B not found; string_scan skipped. Edit DOF_A/DOF_B in run_barrier_pipeline.sh"
fi

echo "Done. Check outputs in: $OUT_DIR"
echo "  - run_*_kM*_L*_wa*_hS*.png (combined 2D+1D per hyperparam)"
echo "  - ensemble_pca.png (and ensemble_umap/pacmap if libs installed)"
echo "  - path_string.png + path_hessian_* (if DOF_A/DOF_B existed)"