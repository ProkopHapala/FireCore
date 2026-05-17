#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="$ROOT"
OUT_DIR="$WORKDIR/grid_search_out"
XYZ="../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz"
DOF="dofSelection_MorseSR_H2O_epairOnly.dat"
NPZ="$OUT_DIR/ensemble_data.npz"
OUT_PREFIX="$OUT_DIR/ensemble"

# Add PN flag argument parsing
PN_FLAG=""
if [[ "${1:-}" == "--pn" ]]; then
    PN_FLAG="--pn"
    echo "Running in PN mode"
fi

# Endpoints for string scan (will be auto-detected if not set manually)
DOF_A=""
DOF_B=""
PATH_PREFIX="$OUT_DIR/path"

echo "[1/3] Running grid search -> ensemble and combined plots ..."
cd "$WORKDIR"
#python grid_search.py --xyz "$XYZ" --dof_file "$DOF" --out_dir "$OUT_DIR" $PN_FLAG || { echo "grid_search failed"; exit 1; }

echo "[2/3] Running DR (PCA/UMAP/PacMAP if available) on ensemble ..."
if [[ ! -f "$NPZ" ]]; then
  echo "ERROR: Ensemble NPZ not found: $NPZ"
  exit 1
fi
python plot_dr.py --trajectory_npz "$NPZ" --out_prefix "$OUT_PREFIX" $PN_FLAG || { echo "plot_dr failed"; exit 1; }

echo "[2.5/3] Clustering ensemble (UMAP/PCA) ..."
if [[ ! -f "$NPZ" ]]; then
  echo "ERROR: Ensemble NPZ not found: $NPZ"
  exit 1
fi
python cluster_ensemble.py --data "$NPZ" --out_prefix "$OUT_PREFIX" $PN_FLAG || { echo "cluster_ensemble failed"; exit 1; }

echo "[3/3] String scan + Hessian (NEB-like) ..."
if [[ -z "$DOF_A" || -z "$DOF_B" ]]; then
  echo "Auto-detecting endpoints from $NPZ ..."
  IDX_A=$(python3 -c "import numpy as np; d=np.load('$NPZ'); print(np.argmin(d['error']))")
  IDX_B=$(python3 -c "import numpy as np; d=np.load('$NPZ'); print(np.argmax(d['error']))")
  
  DOF_A=$(ls $OUT_DIR/run_${IDX_A}_*.dat 2>/dev/null | head -n 1 || true)
  DOF_B=$(ls $OUT_DIR/run_${IDX_B}_*.dat 2>/dev/null | head -n 1 || true)
  
  echo "Selected DOF_A (min error): $DOF_A"
  echo "Selected DOF_B (max error): $DOF_B"
fi

if [[ -n "$DOF_A" && -n "$DOF_B" && -f "$DOF_A" && -f "$DOF_B" ]]; then
  python string_scan.py --dofA "$DOF_A" --dofB "$DOF_B" --xyz "$XYZ" --n_points 11 --relax --hessian --out_prefix "$PATH_PREFIX" $PN_FLAG || { echo "string_scan failed"; exit 1; }
  
  if [[ -f "${PATH_PREFIX}_trajectory.npz" ]]; then
    python plot_dr.py --traj_file "${PATH_PREFIX}_trajectory.npz" --out_prefix "${PATH_PREFIX}_dr" $PN_FLAG || { echo "plot_dr failed"; exit 1; }
  else
    echo "ERROR: path trajectory not found at ${PATH_PREFIX}_trajectory.npz"
    exit 1
  fi
else
  echo "ERROR: DOF_A or DOF_B not found; string_scan cannot proceed."
  exit 1
fi

echo "Done. Check outputs in: $OUT_DIR"
