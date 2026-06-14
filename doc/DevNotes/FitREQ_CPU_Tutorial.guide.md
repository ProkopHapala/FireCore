# FitREQ CPU Tutorial (student guide)

This tutorial explains how to use the CPU version of FitREQ for H-bond 2D scans using the script `tests/tFitREQ/opt_2D.py`. It focuses on what to change, why, and the physical/chemical meaning behind parameters.

References:
- CPU core: `cpp/common/molecular/FitREQ.h`
- Python API: `pyBall/FitREQ.py`
- User script: `tests/tFitREQ/opt_2D.py`
- Types: `cpp/common_resources/AtomTypes.dat`
- DOF selection (Morse): `tests/tFitREQ/dofSelection_MorseSR.dat`

---

## 1) Physical model (REQ + H)

Per atom type we fit 4 parameters `REQH = (R, sqrt(E), Q, H)`:
- R: vdW radius (Å) — repulsive core/contact distance.
- sqrt(E): sqrt of Lennard-Jones/Morse well depth (sqrt(eV)) — controls attraction strength.
- Q: base partial charge (e) — electrostatics magnitude/sign.
- H: H-bond correction amplitude — short-range, directional enhancement for donors/acceptors.

Global parameters (`fit.setGlobalParams`):
- kMorse: curvature of Morse well near minimum (stiffness proxy).
- Lepairs: scale of electron-pair terms (strength of lone-pair/σ-hole contribution).

Electron pairs (epairs): virtual sites representing lone pairs or σ-holes are added to emphasize H-bond directionality. Enable via `bAddEpairs=True` (Python) and `Epairs=1` (C++/setup).

---

## 2) Inputs and where you control things

- `AtomTypes.dat` (types and defaults)
  - Columns include: `RvdW` (used as R), `EvdW` (well depth, used as sqrt(E)=sqrt(EvdW)), `Qbase` (Q), `Hb` (H). See `FitREQ.h: initAllTypes()` where `typeREQs0 = {RvdW, sqrt(EvdW), Qbase, Hb}`.
  - Examples: `N_3, O_3, H_O, H_N` specialize chemistry (sp3 N, sp3 O, H on O, H on N...).

- `dofSelection_MorseSR.dat` (which DOFs are active and how)
  - Format per `FitREQ.h: loadDOFSelection()`:
    `typename comp Min Max xlo xhi Klo Khi K0 xstart invMass`
    - `typename`: e.g., `N_3`, `O_2`, `H_O`, or pseudo `E_O3`, `E_HO`.
    - `comp`: 0=R, 1=sqrt(E), 2=Q, 3=H.
    - `Min/Max`: hard bounds (parameter limits).
    - `xlo/xhi`: regularization targets (soft walls).
    - `Klo/Khi/K0`: regularization stiffness (penalty strength near/between targets).
    - `xstart`: initial value (default to type table if omitted).
    - `invMass`: inverse mass for MD-like optimizer (bigger → faster response).
  - Lines with `E_*` control pseudo-types for epairs (their Q/H or R as needed).
  - Lines with `N_*`, `O_*`, `H_*` control host atom types.

- Reference geometries and energies (`all.xyz`)
  - Built by `opt_2D.py` from donor/acceptor directories using `fit.combine_fragments()` and `fit.concatenate_xyz_files()`.
  - Comment lines hold energies and scan coordinates; parsed by `fit.read_xyz_data()` to get `Erefs, x0s` for weighting/plots.

---

## 3) What `opt_2D.py` does (step-by-step)

1. Select dataset
   - Edit `ref_path`, `donors`, `acceptors`. The script builds all donor–acceptor pairs and concatenates their `.xyz` scans into `all.xyz`.
   - `fit.extract_comments_and_types()` reads types from xyz; `fit.add_epair_types(type_names)` extends the type table with `E_*` pseudo-types if needed.

2. Choose model and verbosity
   - `bMorse=True` selects Morse-SR; else Lennard-Jones-SR.
   - `fit.setVerbosity(verbosity, PrintDOFs=1, PrintfDOFs=1, ...)` prints DOF values/forces each step.

3. Load types and DOFs
   - `fit.loadTypes()` reads `AtomTypes.dat` into `typeREQs0`.
   - `fit.loadDOFSelection(...)` activates specific `(type, component)` DOFs with bounds/regularization from `dofSelection_*.dat`.

4. Load data
   - `fit.loadXYZ('all.xyz', bAddEpairs, bOutXYZ)` loads all samples; with `bAddEpairs=True` lone pairs/σ-holes are added and re-ordered (see `FitREQ.h:addAndReorderEpairs`).
   - `Erefs, x0s = fit.read_xyz_data('all.xyz')` pulls target energies and scan coordinates for weighting.

5. Global physics and setup
   - `fit.setGlobalParams(kMorse=1.8, Lepairs=1.0)` sets global constants.
   - `fit.setup(imodel=5, EvalJ=1, WriteJ=1, Regularize=1, Epairs=bAddEpairs)` picks the CPU model variant and toggles options.
     - Typical values: LJ/LJSR → `imodel=1`; Morse/MorseSR → `imodel=2 or 5` (script uses 5).

6. Weighting of samples
   - `split_and_weight_curves(Erefs, x0s, n_before_min, weight_func)` increases weight near minima and near smooth parts of curves.
   - Examples in script:
     - Morse: `n_before_min=100` (wide lead-in) for smoother wells.
     - LJ: `n_before_min=2` (steeper near-contact).
   - `fit.setWeights(weights0)` applies weights.

7. Optional filtering
   - `fit.setFilter(EmodelCutStart, EmodelCut, ...)`: suppress or down-weight over-repulsive points (too small R).
   - Useful when initial parameters cause diverging repulsion.

8. Baseline evaluation
   - `E, Es, Fs = fit.getEs(bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False)` computes model energies before fitting.

9. Optimization
   - `fit.run(iparallel=0, ialg=1, nstep=10, Fmax=1e-8, dt, damping, max_step, bClamp)` updates DOFs.
     - `ialg=0`: gradient descent.
     - `ialg=1`: dynamical descent (heavy-ball), uses `dt`, `damping` and `invMass` from DOF selection.
     - `bClamp`: clamp steps to respect `Min/Max`.

10. Plotting
   - After optimization, the script recomputes energies and plots per-pair panels: `plot_Epanels_diff_separate(..., save_prefix="opt_2D_")` and saves `opt_2D.png`.

---

## 4) How to modify for your system

- Donor/acceptor sets: change `donors`, `acceptors` lists. Use names consistent with your directory structure under `ref_path`.
- Model family: set `bMorse=True/False` and pick matching DOF file (Morse vs LJ SR).
- DOFs to fit: edit `dofSelection_*.dat` to activate/deactivate types and components (R/E/Q/H), set sensible bounds and regularization.
- Epairs: keep `bAddEpairs=True` for H-bonds; set `Lepairs` to scale their effect (0–1 typical).
- Weighting: adjust `n_before_min` and `exp_weight_func(alpha, a)` to emphasize the well; check weights visually with `fit.plotEWs()`.
- Optimizer: start with `ialg=1` (damped dynamics). Tune `dt`, `damping`, `max_step`; increase `nstep` as needed. Use `ialg=0` for small clean updates.
- Filtering: raise `EmodelCut` if many early points are discarded; lower it to be stricter.

---

## 5) Chemical notes for H-bonds

- Donors: X–H groups (`H_O`, `H_N`) with H pointing to acceptor.
- Acceptors: lone pairs on O/N (`O_2`, `O_3`, `N_3`, etc.) represented explicitly via epairs `E_O*`, `E_N*`.
- Q controls electrostatic attraction/repulsion magnitude; H adds short-range directionality beyond isotropic vdW.
- R, sqrt(E) tune sterics and dispersion balance. For Morse, `kMorse` affects well curvature; for LJ-SR, short-range regularization softens r^(-12).

---

## 6) Key functions you will touch (Python)

- Types/DOFs: `fit.loadTypes()`, `fit.loadDOFSelection()`.
- Data: `fit.combine_fragments()`, `fit.concatenate_xyz_files()`, `fit.loadXYZ()`, `fit.read_xyz_data()`.
- Globals: `fit.setGlobalParams(kMorse, Lepairs)`, `fit.setup(imodel, ..., Epairs)`.
- Weight/filter: `fit.split_and_weight_curves()`, `fit.setWeights()`, `fit.setFilter()`.
- Eval/opt: `fit.getEs()`, `fit.run()`.
- Plots: `fit.slice_and_reshape()`, `fit.plot_Epanels_diff_separate()`.

---

## 7) Minimal example (Morse-SR)

```python
from pyBall import FitREQ as fit
import numpy as np
fit.loadTypes()
fit.loadDOFSelection('tests/tFitREQ/dofSelection_MorseSR.dat')
nb = fit.loadXYZ('all.xyz', bAddEpairs=True, bOutXYZ=False)
Erefs, x0s = fit.read_xyz_data('all.xyz')
fit.setGlobalParams(kMorse=1.8, Lepairs=1.0)
fit.setup(imodel=5, EvalJ=1, WriteJ=1, Regularize=1, Epairs=1)
ws,_ = fit.split_and_weight_curves(Erefs, x0s, n_before_min=100,
    weight_func=lambda E: fit.exp_weight_func(E, a=1.0, alpha=4.0))
fit.setWeights(ws)
fit.setFilter(EmodelCutStart=0.0, EmodelCut=0.5, PrintOverRepulsive=-1,
              DiscardOverRepulsive=-1, SaveOverRepulsive=-1, ListOverRepulsive=-1)
E,Es,Fs = fit.getEs(bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False)
Err = fit.run(iparallel=0, ialg=1, nstep=50, Fmax=1e-8, dt=0.5, damping=0.1,
              max_step=-1, bClamp=True)
```

---

## 8) Troubleshooting

- No type found when loading DOFs: ensure your `typename` exists (or add via `fit.add_epair_types(...)` before `loadDOFSelection`).
- Diverging repulsion: increase `R` start, enable filtering (`EmodelCut ~ 0.5–2.0`), clamp steps, or tighten `Min`.
- Flat optimization: increase weights near minima, raise `dt` slightly, or reduce `K*` regularization.
- Inconsistent panels: verify `marks, angle_data = fit.mark_molecule_blocks(comments)` and that directories are combined as expected.

## 9) Python API cheatsheet (exact signatures)

- __verbosity__: `fit.setVerbosity(verbosity=1, idebug=0, PrintDOFs=0, PrintfDOFs=0, PrintBeforReg=0, PrintAfterReg=0)`
- __globals__: `fit.setGlobalParams(kMorse=1.6, Lepairs=0.5)`
- __setup__: `fit.setup(imodel=1, EvalJ=0, WriteJ=0, CheckRepulsion=0, Regularize=0, AddRegError=0, Epairs=0, BroadcastFDOFs=0, UdateDOFbounds=0, EvalOnlyCorrections=0)`
- __filter__: `fit.setFilter(EmodelCut=1.0, EmodelCutStart=None, EmodelCutFactor=0.75, iWeightModel=1, ListOverRepulsive=0, SaveOverRepulsive=0, PrintOverRepulsive=0, DiscardOverRepulsive=0, WeightByEmodel=0)`
- __types__: `fit.loadTypes(fEtypes="data/ElementTypes.dat", fAtypes="data/AtomTypes.dat")`
- __DOFs__: `fit.loadDOFSelection(fname="dofSelection.dat") -> int`
- __weights__: `fit.setWeights(weights)` / `fit.loadWeights(fname="weights.dat") -> int`
- __data__: `fit.loadXYZ(fname, bAddEpairs=False, bOutXYZ=False, bEvalOnlyCorrections=False) -> nbatch`
- __evaluate__: `fit.getEs(Es=None, Fs=None, bOmp=False, bDOFtoTypes=True, bEs=True, bFs=False, xyz_name=None) -> (Eerr, Es, Fs)`
- __opt__: `fit.run(ialg=2, iparallel=1, nstep=100, Fmax=1e-8, dt=0.01, max_step=0.05, damping=0.0, bClamp=False) -> Err`
- __scans__: `fit.scanParam(iDOF, xs, Es=None, Fs=None, bEvalSamples=True) -> (Es, Fs)`
- __2D scans__: `fit.scanParam2D(iDOFx, iDOFy, xs, ys, Es=None, Fx=None, Fy=None, bEvalSamples=True) -> (Es, Fx, Fy)`
- __export xyz__: `fit.exportSampleToXYZ(i, fname)`; `fit.exportAllSamplesToXYZ(fname)`
- __type<->DOF__: `fit.setTypeToDOFs(i, REQ)`; `fit.getTypeFromDOFs(i, REQ=None) -> REQ`
- __helpers__:
  - `fit.read_xyz_data(fname) -> (Erefs, x0s)`
  - `fit.extract_comments_and_types(fname) -> (type_names, comment_lines)`; `fit.add_epair_types(type_names)`
  - `fit.combine_fragments(frags1, frags2, path=None, ext="") -> dirs`;
    `fit.concatenate_xyz_files(directories=None, base_path='./', fname="all.xyz", output_file="all.xyz", mode='w') -> marks`
  - weighting: `fit.split_and_weight_curves(Erefs, x0s, n_before_min=4, weight_func=None, EminMin=-0.02) -> (weights, lens)`; `fit.exp_weight_func(Erefs, a=1.0, alpha=3.0, Emin0=0.1)`
  - plotting: `fit.plotEWs(...)`, `fit.plot_Epanels_diff_separate(Emodels, Erefs, ref_dirs, ...)`

This guide is intentionally concise. For deeper internals, inspect `FitREQ.h` (model flags, DOF machinery) and `FitREQ.py` (wrappers, helpers).
