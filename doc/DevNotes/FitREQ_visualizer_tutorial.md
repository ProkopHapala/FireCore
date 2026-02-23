---
description: FitREQ/PN visualization quickstart and internals
---

# FitREQ/PN Visualization Tutorial

Audience: PhD users (run and interpret plots) and developers (extend or debug).

## What was added
- C++ API: `getSampleGeom` now returns geometry plus fragment split `n0` (includes epairs), `evalSamplePairs` returns per-pair components (Morse/vdW, Coulomb, Hcorr, SR) for a sample.
- Python ctypes (`pyBall/FitREQ_PN.py`): bindings updated to expose `getSampleGeom` (with `n0`) and pairwise evaluation; wrappers used by visualizers.
- Static visualizer (`tests/tFitREQ_PN/visualize_static.py`): non-interactive plots and diagnostics, label modes, interaction matrices, argparse for inputs/params.
- GUI visualizer (`tests/tFitREQ_PN/visualizer_gui.py`): interactive maps/slices/3D geometry with clickable matrix highlighting and label-mode switching; anchored colorbar, equal 3D aspect, smaller markers; argparse with sane defaults.

## Prerequisites
- Data files: `tests/tFitREQ_PN/data/AtomTypes.dat`, `ElementTypes.dat`, `dofSelection_run.dat`.
- Test XYZs: `tests/tFitREQ_PN/wb97m-split/H2O-A1_H2O-D1-z.xyz` (default), `H2O-D1_H2O-A1-y.xyz` (alt).
- Epairs are added on load (`bAddEpairs=True`); `Lepairs` default 0.8.

## Quick start (recommended)
From repo root:
```bash
cd tests/tFitREQ_PN

# Static stills + diagnostics (saves PNGs in cwd)
python3 visualize_static.py \
  --input wb97m-split/H2O-A1_H2O-D1-z.xyz \
  --lepairs 0.8 --kmorse 1.8 --label_mode chem_no_unders

# Interactive GUI
python3 visualizer_gui.py \
  --input wb97m-split/H2O-A1_H2O-D1-z.xyz \
  --lepairs 0.8 --kmorse 1.8 --label_mode chem_no_unders
```
Common CLI options (both scripts):
- `--input`: xyz path (abs or relative to script dir)
- `--dof_selection`: DOF file (default `data/dofSelection_run.dat`)
- `--atom_types`, `--element_types`: type tables
- `--label_mode`: `chem_no_unders` | `chem` | `numeric`
- `--lepairs`: epair distance (default 0.8)
- `--kmorse`: Morse curvature (default 1.8)

## What the static visualizer produces
- 2D maps: DFT, model, diff with physical ticks (angle, distance).
- 1D slices: distance and angle cuts showing baseline (Coul+vdW), +Hcorr, +SR, and DFT.
- 3D geometry snapshot: fragments split by `n0`, labeled atoms, epairs included, equal aspect.
- Interaction matrices: components per atom pair.
- Diagnostics: printed minima, component tables, etc.
- Saved PNGs: `plot_2d_maps.png`, `plot_1d_slices.png`, `plot_3d_geom.png`, `plot_mat_*.png`.

## How to drive the GUI
- Click any map (DFT/model/diff) to move the crosshair and update slices/3D/matrix.
- Click the interaction matrix to set the focus atom; corresponding 3D lines highlight.
- Radio buttons: choose interaction labeling on 3D lines (Distance, Morse, Coulomb, H-corr, SR).
- Label-mode radio: `chem_no_unders` / `chem` / `numeric`.
- Prev/Next buttons: cycle focus atom. Crosshair shows current (angle, distance) index.
- Colorbar is anchored; 3D view has equal x/y/z scale; atom markers are small dots.

## Data flow (developer view)
1) Setup: `setVerbosity` → `setModel(ivdW=4, iCoul=1, iHbond=3, Epairs=1, iEpairs=2, kMorse, Lepairs, bPN=True)` → `loadTypes` → `loadDOFSelection`.
2) `loadXYZ(..., bAddEpairs=True)` parses `n0`, adds/reorders epairs; `addAndReorderEpairs` updates `atoms->n0` (includes epairs).
3) Energies: `getEs` (total) and `getEs_components` (Coul, vdW/Morse, Hcorr, SR) per sample.
4) Geometry/per-pair: `getSampleGeom` returns positions, types, charges, host map, and `n0`; `evalSamplePairs` returns per-pair components for the current sample.
5) Grids: `parse_xyz_mapping` maps sequence → (idist, iang) so energies/components can be placed on 2D angle–distance arrays.
6) Visuals: plots use physical axes; label modes come from AtomTypes (`chem`, `chem_no_unders`, `numeric`).

## Notes on epairs and fragments
- `n0` from `getSampleGeom` is the fragment split including epairs. Fragment 1: indices `< n0`; Fragment 2: `>= n0`.
- Epairs are placed at distance `Lepairs` along host directions; interactions with epairs contribute to SR component.

## Typical troubleshooting
- Colorbar/matrix drift: fixed by anchored cbar in GUI; ensure current version.
- Fragment split wrong: confirm `bAddEpairs=True` and using `n0` returned by `getSampleGeom`.
- Missing labels: ensure AtomTypes path is correct; label_mode governs formatting.

## Extending
- Add more CLI knobs (ivdW/iCoul/iHbond/iEpairs) mirroring `run.py` if needed.
- For larger systems, compute per-pair data on-demand via `evalSamplePairs` to avoid big tensors.
- To change default test, pass `--input` or edit the default in the scripts.

## File pointers
- Core C++: `cpp/common/molecular/FitREQ_PN.h` (model, epairs, eval), `cpp/libs/Molecular/FitREQ_PN_lib.cpp` (C API exports).
- Python bridge: `pyBall/FitREQ_PN.py` (ctypes bindings, helpers).
- Visualizers: `tests/tFitREQ_PN/visualize_static.py`, `tests/tFitREQ_PN/visualizer_gui.py`.

Happy plotting!
