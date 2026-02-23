# FitREQ Barrier Analysis Tutorial

A practical guide for students diagnosing hydrogen-bond fitting runs in FitREQ. Focus: spotting competing minima (“basins”), visualizing barriers, and understanding which parameters drive conflicts.

## 1) Background (why barriers appear)
- The objective (fitting error vs DOFs) is non-convex; different hyperparameters (e.g., `kMorse`, `Lepairs`, `weight_alpha`, `hScale`) can land you in distinct basins.
- A “barrier” or “ridge” is revealed when interpolating between two solutions: energy spikes in the middle → competing physics.
- Diagnostics:
  - 2D maps (angle vs distance) and 1D minima (r_min, E_min)
  - Ensemble clustering (PCA/UMAP/PacMAP) on optimized DOFs
  - String/NEB-like scans in parameter space
  - Hessian eigenvectors/correlations at key points (soft modes show redundant/competing parameters)

## 2) Prerequisites
- Dataset: `tests/tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz` (double-well reference)
- DOF file: epair-only (no Q), e.g., [tests/tFitREQ/dofSelection_MorseSR_H2O_epairOnly.dat](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ/dofSelection_MorseSR_H2O_epairOnly.dat:0:0-0:0)
- Python deps: numpy, matplotlib; optional: umap-learn, pacmap, scikit-learn
- Ensure [pyBall/FitREQ.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/pyBall/FitREQ.py:0:0-0:0) uses the combined plotting with correct normalization (shift to min, shared vmin/vmax from ref).
- DOF metadata: each `.dat` may carry hyperparameters in a `#HYPER kMorse=... Lepairs=... weight_alpha=... hScale=...` line. Grid search now writes these automatically; string scans read and interpolate them.

## 3) Scripts overview (all in `tests/tFitREQ/`)
- [grid_search.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ/grid_search.py:0:0-0:0): Hyperparameter sweep; outputs combined 2D+1D figures per run and `ensemble_data.npz`.
- `compare_endpoints.py`: Side-by-side rigid/relaxed comparison of two DOF files (2D maps + 1D minima).
- `string_scan.py`: Interpolate between endpoints (DOFs + hyperparams), evaluate rigid/tethered paths, compute Hessians at α points, save `path_trajectory.npz`.
- [plot_dr.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ/plot_dr.py:0:0-0:0): Dimensionality reduction (PCA/UMAP/PacMAP) of trajectories/ensemble, colored by log(error).
- `barrier_utils.py`: Shared utilities (setup, evaluate, tethered relax, FD Hessian).

## 4) Quickstart: hyperparameter sweep (ensemble generation)
```bash
cd tests/tFitREQ
python grid_search.py \
  --xyz ../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz \
  --dof_file dofSelection_MorseSR_H2O_epairOnly.dat \
  --out_dir grid_search_out
```
- Sweeps `kMorse=[1.6,1.7,1.8]`, `Lepairs=[0.6,0.8,1.0,1.2]`, `weight_alpha=[4.0]`, `hScale=[0.0,0.5,1.0]` (edit in script).
- Outputs: one combined PNG per run in `grid_search_out/`, plus `ensemble_data.npz` and `run_*.dat` files with `#HYPER` metadata embedded.
- Check logs for min/max, NaN/Inf diagnostics; verbosity is set to 0.

## 5) Compare two endpoints
Pick two optimized DOF files (e.g., best/worst from grid):
```bash
python compare_endpoints.py \
  --dofA grid_search_out/run_47_kM1.8_L1.2_wa4.0.dat \
  --dofB grid_search_out/run_0_kM1.6_L0.6_wa1.0.dat \
  --xyz H2O-A1_H2O-D1-z.xyz \
  --out_dir compare_out
```
Outputs combined ref/model/diff + minima plots for rigid and relaxed variants.

## 6) Path (string/NEB-like) scan with Hessians
Interpolate DOFs **and** hyperparameters between two endpoints (reads `#HYPER` from the `.dat` files):
```bash
python string_scan.py \
  --dofA grid_search_out/run_MIN.dat \
  --dofB grid_search_out/run_MAX.dat \
  --xyz ../tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-y.xyz \
  --n_points 11 \
  --relax --hessian \
  --out_prefix grid_search_out/path
```
Outputs:
- `path_string.png` (energy barrier, **relaxed** parameter evolution with y-lim ±1, gradients)
- Hessian analyses at α=0,0.5,1:
  - `path_hessian_softmode_alpha_*.png` (soft-mode eigenvector composition → which params co-vary)
  - `path_hessian_corr_alpha_*.png` (correlation heatmaps from pseudo-inverse Hessian)
- Trajectory data: `path_trajectory.npz` (includes kMorse/Lepairs/weight_alpha/hScale paths)

Interpretation:
- Soft-mode bars: large components on two params → redundancy/competition.
- Correlation heatmaps: red = positively correlated, blue = anti-correlated; compare α=0 vs α=1 to see physics shift across basins.

## 7) Dimensionality reduction (ensemble or path)
Use ensemble (`ensemble_data.npz`) or path (`path_trajectory.npz`):
```bash
python plot_dr.py \
  --trajectory_npz grid_search_out/path_trajectory.npz \
  --out_prefix grid_search_out/dr
```
Outputs: `dr_pca.png` (and `dr_umap.png` / `dr_pacmap.png` if libs installed), showing clustering/bifurcation; colored by log(error), with initial/relaxed path overlays.

## 8) One-click pipeline
To run grid search, DR, auto-select farthest endpoints (min/max error), and do a relaxed + Hessian string scan with hyperparameter interpolation:
```bash
cd tests/tFitREQ
./run_barrier_pipeline.sh
```
Outputs land in `grid_search_out/`: grid run PNGs, ensemble DR, `path_string.png`, `path_hessian_*`, and `path_dr_*`.

## 9) Normalization and zero-line sanity
- All plots use ref-based vmin/vmax (shifted to ref min) shared across ref/model/diff.
- 2D panels include the angle-dependent zero line (from raw DFT grid) mapped to pixel indices.
- r_min is clamped to (1.5, 3.0 Å); E_min y-lims use ref min ± margins.

## 10) Troubleshooting
- If you see 1e+37 or flat model maps: check console min/max diagnostics; ensure epair-Q DOFs are removed; rerun with clean solver init; consider subprocess isolation if OpenMP state is dirty.
- Endpoints in `path_string.png` must match between linear/relaxed at α=0,1; relax only interior points.
- If parameters differ by orders of magnitude, apply manual scaling (*0.1 or *10) and note in legend (string scripts already plot raw values with y-lim ±1).

## 11) Suggested workflow
1. Run [grid_search.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ/grid_search.py:0:0-0:0) → inspect best/worst PNGs and `ensemble_data.npz`.
2. Run [plot_dr.py](cci:7://file:///home/prokophapala/git/FireCore-fitREQH/tests/tFitREQ/plot_dr.py:0:0-0:0) on `ensemble_data.npz` → pick separated basins.
3. Run `string_scan.py` between those basins (hyperparams auto-interpolated from `#HYPER`) → read `path_string.png`, Hessian soft-mode/correlation plots.
4. Optionally re-run `compare_endpoints.py` for a clean side-by-side reference/model/diff check.
5. Or just run `./run_barrier_pipeline.sh` for the full automated path (min/max endpoints).
