# pyTruss module map

This directory collects the Python-side building blocks for the truss / cloth solvers that were migrated from the separate pyTruss project.  The code is split between reusable back-ends (geometry, sparse linear algebra, CPU and GPU solvers) and a handful of command-line utilities/tests that exercise those components.

The sections below list every file, group them by responsibility, note how they relate to each other, and provide concrete command lines that we verified on this workspace.  All commands assume you start in `pyBall/pyTruss/` and have `PYTHONPATH` pointing at the repository root when running from elsewhere.

---

## Core data structures and utilities

| File | Summary | Entry point? | Notes |
| --- | --- | --- | --- |
| `truss.py` | Defines the `Truss` data structure plus helpers to procedurally build ropes, grids, wheels, graph colouring, etc. | Optional demo block (`python truss.py`) that renders with Matplotlib. | Imports from within this folder only; GUI run requires a Matplotlib backend. |
| `plot_utils.py` | Simple Matplotlib helper to draw trusses. | Module only. | Used by demos/tests. |
| `sparse.py` | Sparse linear algebra utilities (Jacobi/GS iterations, neighbour tables, colouring) operating on truss connectivity. | Contains a diagnostic main (`python sparse.py`) that compares dense vs sparse ops and uses Matplotlib for plots. | Provides the neighbour-building helpers used throughout. |
| `IterativeLinearSolvers.py` | Stand-alone iterative solver building blocks plus momentum/Chebyshev variants. | Module only. | Shared by tests and CPU solver. |


## Projective-dynamics CPU prototypes

| File | Summary | Entry point? | Notes |
| --- | --- | --- | --- |
| `projective_dynamics.py` | CPU reference implementation of projective dynamics. Demonstrates solver workflow, builds dense matrices. | Yes: `python projective_dynamics.py` (prints iteration log, attempts to plot). | Exports `make_pd_matrix`, `make_pd_rhs`, etc. **Does not** currently expose `build_neighbor_list`; other modules import that from `sparse.py`. |
| `projective_dynamics_iterative.py` | Alternative projective-dynamics loop with Jacobi/Chebyshev acceleration. | _Main currently fails_: tries to import `build_neighbor_list` from `projective_dynamics`. Needs refactor to import from `sparse`. | Keep if we plan to fix the import; otherwise archive. |
| `truss_solver.py` | CPU production solver class `TrussSolver` with multiple linear solvers (VBD, Jacobi/GS, Chebyshev, momentum), logging, etc. | Module only (instantiated by CLI scripts). | Wraps sparse utilities above; used by `run_vbd_cloth.py` and `run_solver_debug.py`. |


## OpenCL back-ends

| File | Summary | Entry point? | Notes |
| --- | --- | --- | --- |
| `truss_solver_ocl.py` | New GPU backend built on `pyBall/OCL/OpenCLBase.OpenCLBase`. Handles buffer allocation, kernel header parsing, shared utility path with other OCL modules. Exposes solver callbacks (`vbd_serial`, `jacobi_fly`, `jacobi_diff`). | Module only. | Imports the shared OpenCL infrastructure through `from OCL.OpenCLBase import OpenCLBase`. We added a `sys.path` shim so this resolves when running scripts from this folder. **Preferred path going forward**; consider deprecating the legacy file once its features are ported. |
| `truss_ocl.py` | Lower-level OpenCL experiments (buffer-level solvers, convergence comparisons). | Module contains several `test_*` functions; `python truss_ocl.py` runs `test_diff_vs_vbd()` by default. | Uses its own wrapper (`TrussOpenCLSolver`). Mostly for benchmarking; not referenced by the new solver. |
| `truss.cl` | Shared OpenCL kernel source for GPU solvers. | N/A | Consumed by both GPU back-ends. |


## Command-line utilities / entry scripts

These are intended for direct execution.  We verified each command below (using `MPLBACKEND=Agg` for scripts that try to plot).

| Script | Purpose | Verified command |
| --- | --- | --- |
| `run_vbd_cloth.py` | Original cloth runner comparing CPU and legacy GPU solvers. | `python run_vbd_cloth.py --nx 1 --ny 0 --nsteps 1 --niter 1 --cpu 1 --no-plot 1` |
| `run_vbd_cloth_new.py` | Updated cloth runner for the new `truss_solver_ocl` backend; supports toggling CPU/GPU. | `python run_vbd_cloth_new.py --nx 1 --ny 0 --nsteps 1 --niter 1 --cpu 1 --gpu 1 --no-plot 1` |
| `run_solver_debug.py` | CPU-only solver harness with perturbations, comparisons, plotting. | `python run_solver_debug.py --nx 1 --ny 0 --niter 1 --solver vbd --chain 1 --no-plot --perturb-sigma 0` |
| `test_graph_coloring.py` | Quick graph-colouring sanity check. | `python test_graph_coloring.py` |
| `test_Chebyshev_accel.py` | Standalone linear-system convergence demo with Jacobi vs Chebyshev acceleration. | `python test_Chebyshev_accel.py --size 5 --diag-dominance 1.1 --max-iters 5` |
| `projective_dynamics.py` | Dense PD reference with plotting. | `MPLBACKEND=Agg python projective_dynamics.py` |
| `truss.py` | Geometry demo producing plots of a wheel truss. | `MPLBACKEND=Agg python truss.py` |
| `sparse.py` | Sparse/dense comparison, plotting residual curves. | `MPLBACKEND=Agg python sparse.py` |

Scripts with current issues when launched directly:

- `projective_dynamics_iterative.py`: import error (`build_neighbor_list` expected in `projective_dynamics`).  Needs to import from `sparse` or re-export the helper.
- `test_Jacobi_Chebyshev_convergence.py`: same missing `build_grid_2d` re-export; requires a small fix to import directly from `truss` / helper file.

These failures are noted here so we can decide whether to patch or retire those demos.


## How to run from outside the directory

If you prefer calling these scripts from another working directory, set up the module path first:

```bash
export PYTHONPATH=/home/prokop/git/FireCore:$PYTHONPATH
cd /home/prokop/git/FireCore/pyBall/pyTruss
python run_vbd_cloth_new.py --help
```

For GPU runs ensure the OpenCL drivers are installed and the desired platform/device is visible; the new runner prints detected platforms during startup.  Matplotlib-based demos will require a GUI backend unless `MPLBACKEND=Agg` is specified.


## Next steps / cleanup ideas

2. Fix the import paths in `projective_dynamics_iterative.py` and `test_Jacobi_Chebyshev_convergence.py` if we want their demos working out-of-the-box.
3. Consider sharing plotting utilities / CLI options between the old and new cloth runners to avoid drift.

This README should serve as the quick-start guide for anyone exploring the migrated pyTruss components.
