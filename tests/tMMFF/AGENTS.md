# MMFF & Surface Interaction Tests

## Purpose

Tests for MMFFsp3 force field, GridFF surface potentials, molecule-surface interactions, rigid-body dynamics, assembly, and AFM/STM simulation pipelines.

## Ownership

- MMFFsp3 CPU/GPU parity and phonon validation
- GridFF generation, evaluation, and diagnostics
- XYZ rigid surface scanning (Morse/LJ/Coulomb)
- Surface sampling and height-map extraction
- Molecular assembly (packing, clash detection, pose generation)
- Rigid-body dynamics on surfaces
- Path optimization with causal tethers
- Interactive GUI scripts (VisPy-based explorers)

## Local Contracts

- **Run from this directory** — scripts use relative paths to `common_resources/cl/*.cl`.
- **Use `run.sh` / `make.sh` scripts** when available; never invoke `make` directly.
- **GridFF precomputation:** B-spline grids (`Bspline_PLQd.npy`) are generated once and reused.
- **MMFF topology:** Use `MMFF.toMMFFsp3_loc()` for GPU-ready topology; set `nPBC` before calling.
- **Node vs cap atoms:** Only node atoms carry angular terms; caps feel recoil only.
- **Capping / reorder_nodes:** For large aromatics (e.g., DiTetracenoHelice), use `reorder_nodes_first=False` and `capping_atoms=set()` (force_node_all mode) to avoid NaN.
- **PBC shifts:** `nPBC=(1,1,0)` gives 9 shifts; zero shift is at index 4, not 0.

## Work Guidance

### Force Field Parity
- `test_diamond_phonon_bands.py` — phonon band structure validation
- `test_MMFF_ocl_parity.py` — CPU vs GPU parity
- `run_parity_mmff.sh` — automated parity suite

### GridFF
- `run_test_GridFF.py` / `run_test_GridFF_CaF2.py` — GridFF generation and evaluation
- `render_surface_iso.py` — headless surface height-map rendering
- `test_electrostatics_comparison.py` — compare electrostatics implementations

### Surface Scanning & Assembly
- `test_interaction_scan.py` — 1D z-scan (XYZ rigid)
- `test_assembly.py` — molecular assembly CLI driver
- `test_rigid_gridff_surface.py` / `test_rigid_gridff_ptcda_batch.py` — rigid-body relaxation

### Path Optimization
- `ManipulationPathOpt.py` — batch replica optimization with causal tethers

### GUI
- `run_gui.py` — launch interactive VisPy explorer
- `gui_fdbm_fit.py` — FDBM parameter fitting GUI

## Verification

- Run relevant `run.sh` or test script from this directory
- GridFF parity: compare GPU B-spline vs CPU reference for identical substrate/molecule
- Assembly: clash penalties should decrease monotonically with optimization
- Rigid-body: force norms should decay after perturbations

## Child DOX Index

- No child subtrees
