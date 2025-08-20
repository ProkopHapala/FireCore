# tFitREQ run.sh — quick guide

This directory contains small tests/demos for REQ fitting. The `run.sh` helper script rebuilds the C++ shared library and runs the selected Python test.

## What it does
- Rebuilds `libFitREQ_lib.so` in `cpp/Build/libs/Molecular/`:
  - removes the old `libFitREQ_lib.so`
  - runs `make -j4 FitREQ_lib`
- Enables AddressSanitizer via `LD_PRELOAD=$(g++ -print-file-name=libasan.so)`
- Widens terminal for readable logs: `stty cols 200`
- Runs the chosen Python test (default: `test_export.py`) with unbuffered output, piping stdout to a file while keeping it on screen:
  - stderr -> `asan.log`
  - stdout -> `OUT-FitREQ` (via `tee`)

Script: `tests/tFitREQ/run.sh`

```bash
#!/bin/bash

wd=`pwd`

cd ../../cpp/Build/libs/Molecular/
pwd
rm   libFitREQ_lib.so
make -j4 FitREQ_lib
cd $wd

# ------- asan (Memory Sanitizer)
LD_PRELOAD=$(g++ -print-file-name=libasan.so)
echo   $LD_PRELOAD
export LD_PRELOAD

stty cols 200   # set terminal width

echo "#=========== RUN"
python3 -u test_export.py 2> asan.log | tee OUT-FitREQ
```

## Usage
Run from this directory:
```bash
bash run.sh
```
Artifacts:
- `OUT-FitREQ` — captured stdout
- `asan.log` — captured stderr (ASan reports, tracebacks)
- rebuilt `cpp/Build/libs/Molecular/libFitREQ_lib.so`

## Prerequisites
- Linux, Bash, `make`, `g++`
- Python 3
- Build deps for `FitREQ_lib` (see top-level build docs)
- AddressSanitizer available via GCC (`libasan.so`). If missing, disable ASan (see below).

Optional environment (commented in script):
- `LD_LIBRARY_PATH` for MKL or GCC runtime if needed

## Customization
- Switch test: uncomment one of the Python lines in `run.sh` (e.g., `opt_mini.py`, `opt_2D.py`, `opt_check_derivs.py`, `plot_DOF_trj.py`).
- Change parallelism: edit `make -j4`.
- Disable ASan: comment out the `LD_PRELOAD` lines or run `env -u LD_PRELOAD bash run.sh`.
- Non-interactive terminals: if `stty` fails, comment that line.

## Troubleshooting
- Build fails: re-run the build step manually in `cpp/Build/libs/Molecular/` and inspect errors.
- `libasan.so` not found: install GCC with ASan support or disable ASan as above.
- Permission/TTY issues with `stty`: comment the line.
- Hanging test: run the Python script directly with `-u` and consider adding timeouts; check `asan.log` for leaks/invalid memory.

## Related files
- C++: `cpp/common/molecular/FitREQ.h`
- Python: `pyBall/FitREQ.py`
- Notes: `doc/DevNotes/FitREQ.md`

---

# Parameter reference (concise)

All functions live in `pyBall/FitREQ.py` and wrap the C++ library.

* __setVerbosity(verbosity, idebug=0, PrintDOFs=0, PrintfDOFs=0, PrintBeforReg=0, PrintAfterReg=0)__
  Controls logging. In tests we often set `PrintDOFs=1`, `PrintfDOFs=1` to trace DOFs and forces.

* __setGlobalParams(kMorse=1.6, Lepairs=0.5)__
  Global constants: `kMorse` (Morse well stiffness), `Lepairs` (electrostatic-pair scale). Used before `setup()`.

* __setup(imodel, EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, AddRegError=0, Epairs=0, BroadcastFDOFs=0, UdateDOFbounds=0, EvalOnlyCorrections=0)__
  Selects model and toggles features.
  - imodel: model variant. Used values in tests: 1 (LJ/LJSR), 2 or 5 (Morse/MorseSR).
  - EvalJ/WriteJ: evaluate and/or print Jacobian-like contributions; keep `=1` during optimization for diagnostics.
  - Regularize/AddRegError: parameter regularization.
  - Epairs: include electrostatic pair terms.
  - CheckRepulsion/UdateDOFbounds/BroadcastFDOFs/EvalOnlyCorrections: expert toggles, rarely needed for basic runs.

* __loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat")__
  Loads element/atom-type tables.

* __loadDOFSelection(fname="dofSelection.dat")__
  Selects which type DOFs (R,E,Q,H,...) are fitted. See `dofSelection_*.dat` files here.

* __loadXYZ(fname, bAddEpairs=False, bOutXYZ=False, bEvalOnlyCorrections=False)__
  Loads reference xyz with per-structure energy/comments. `bAddEpairs=True` auto-adds `E_*` epair pseudo-types.

* __read_xyz_data(fname)__ → `(Erefs, x0s)`
  Parses comment lines (e.g., `# n0 ... Etot <val> x0 <val> ...`). Used for weighting/plots.

* __split_and_weight_curves(Erefs, x0s, n_before_min=..., weight_func=...)__
  Splits curves by monotonic `x0` segments and builds weights emphasizing minima regions.

* __setWeights(weights)__
  Assigns sample weights; length must match `nbatch` from `loadXYZ()`.

* __getBuffs()__
  Exposes internal buffers/metadata as numpy views: `nDOFs, ntype, nbatch, DOFs, fDOFs, vDOFs, fDOFbounds, typeREQs*, weights`.

* __setFilter(EmodelCut, EmodelCutStart=..., iWeightModel=..., List/Save/Print/DiscardOverRepulsive, WeightByEmodel=0)__
  Optional runtime filter against over-repulsive samples; handy for debugging bad regions.

* __getEs(Es=None, Fs=None, bOmp=False, bDOFtoTypes=True, bEs=True, bFs=False, xyz_name=None)__ → `(Eerr, Es, Fs)`
  Evaluates energies/forces; `bOmp` toggles OpenMP; `bDOFtoTypes=False` keeps DOFs decoupled from type table; `xyz_name` dumps xyz.

* __run(ialg=2, iparallel=1, nstep=100, Fmax=1e-8, dt=0.01, max_step=0.05, damping=0.0, bClamp=False)__ → `Err`
  Optimizer. In tests: `ialg=0` (GD), `1` (MD/dynamical descent), `2/3` (Barzilai–Borwein short/long per comments). Tune `dt`, `damping`, `max_step`, `bClamp`.

* __scanParam(iDOF, xs, ...), scanParam2D(iDOFx, iDOFy, xs, ys, ...)__
  1D/2D scans of energy and forces along chosen DOFs for diagnostics.

# Test workflow internals

Typical flow used by `opt_mini.py` and `opt_2D.py`:

1. __Verbosity/plots__: `fit.plt = plt`; `fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, ...)`.
2. __Types/DOFs__: `fit.loadTypes()`; `fit.loadDOFSelection("dofSelection_*.dat")`.
3. __Data__: `fit.loadXYZ(fname, bAddEpairs, bOutXYZ, bEvalOnlyCorrections)`; also `Erefs,x0s = fit.read_xyz_data(fname)`.
4. __Model__: `fit.setGlobalParams(kMorse=..., Lepairs=...)`; `fit.setup(imodel=..., EvalJ=1, WriteJ=1, Regularize=1, Epairs=bAddEpairs)`.
5. __Weights__: `weights, lens = fit.split_and_weight_curves(Erefs, x0s, ...)`; `fit.setWeights(weights)`.
6. __Buffers__: `fit.getBuffs()` exposes `fit.DOFs`, `fit.fDOFbounds`, etc., for inspection.
7. __Filtering (opt)__: `fit.setFilter(EmodelCut=..., ...)` to avoid over-repulsive outliers.
8. __Baseline__: `E,Es,Fs = fit.getEs(bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False)`; quick plot via `fit.plotEWs(...)`.
9. __Optimize__: `Err = fit.run(ialg=..., nstep=..., dt=..., damping=..., max_step=..., bClamp=...)`.
10. __Post__: re-evaluate `getEs`, plot panels, optionally `exportAllSamplesToXYZ()`.

# Script roles in this folder

* __`test_export.py`__ — basic loader/export sanity-check for types + xyz + optional epairs; prints summary.
* __`opt_mini.py`__ — minimal single-dataset optimization; demonstrates weights, filters, and GD/MD choices.
* __`opt_2D.py`__ — multi-system 2D H-bond maps; concatenates xyz across donors/acceptors, slices panels, runs a short optimization, and saves per-pair figures.
* __`opt_check_derivs.py`__ — derivative checks: scans DOFs, compares analytic vs numeric derivatives, plots forces panels.
* __`plot_DOF_trj.py`__ — parses logs and plots DOF and fDOF trajectories across optimization steps.

# Data/config files here

* __DOF selections__: `dofSelection_*.dat` control which type-components are active (R,E,Q,H,...). Use `fit.comment_non_matching_lines()` if your selection includes types absent in the current dataset; `fit.add_epair_types()` can extend with `E_*` pseudo-types when `bAddEpairs=True`.
* __Input xyz__: see `all.xyz`, `input_example*.xyz`, scans like `scan_OH.xyz`. Comments hold energy and geometry parameters consumed by `read_xyz_data()`.
* __Types tables__: default paths `data/ElementTypes.dat`, `data/AtomTypes.dat` under repo root.

# Tips and pitfalls

* __Imodel choices__: scripts use `imodel=1` (LJ/LJSR) and `imodel=2/5` (Morse variants). Match this with a consistent `dofSelection_*.dat` (e.g., LJ vs Morse SR files).
* __Epairs__: set `bAddEpairs=True` in `loadXYZ()` and `Epairs=1` in `setup()` to include epair terms; `setGlobalParams(Lepairs=...)` scales them.
* __Weights__: start with `split_and_weight_curves()`; tweak `n_before_min` and `exp_weight_func(alpha, a)` to emphasize minima.
* __Filtering__: `setFilter(EmodelCut=...)` can suppress exploding repulsive samples while tuning DOF bounds.
* __Units__: plotting helpers use `ev2kcal=23.060548` for kcal/mol; raw energies are in eV unless you convert.
* __OpenMP parity__: use `getEs(..., bOmp=False/True, bFs=True)` and compare arrays to verify CPU/OMP parity (see `test_getEs_openmp()` helper in `FitREQ.py`).

# Minimal example (outline)

```python
from pyBall import FitREQ as fit
fit.loadTypes()
fit.loadDOFSelection('dofSelection_LJSR2.dat')
nbatch = fit.loadXYZ('all.xyz', bAddEpairs=True, bOutXYZ=False)
Erefs,x0s = fit.read_xyz_data('all.xyz')
fit.setGlobalParams(kMorse=1.8, Lepairs=1.0)
fit.setup(imodel=1, EvalJ=1, WriteJ=1, Regularize=1, Epairs=1)
weights,_ = fit.split_and_weight_curves(Erefs, x0s, n_before_min=4)
fit.setWeights(weights); fit.getBuffs()
E,Es,Fs = fit.getEs(bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False)
Err = fit.run(ialg=1, nstep=100, dt=0.1, damping=0.1, max_step=-1, bClamp=True)
```

