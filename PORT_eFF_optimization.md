# Port: eFF optimization improvements (prokop) -> Zverina branch

Working notes / plan for porting Prokop's eFF-optimization work into Tomas Zverina's
FireCore branch and testing it on Tomas's H2O test system.

**Decisions (user, 2026-11-24):**
- Port **everything**: optimizer machinery + physics bugfixes.
- Keep the **original implementation runnable** for A/B comparison
  -> build `libeFF_lib.so` BEFORE any edits and keep it as `libeFF_lib_origZ.so`.
- Work on branch `Zverina`, **do NOT commit** until user approves.
- `origin/prokop` has an unrelated (rewritten) git history => no cherry-pick;
  port is done file-by-file.

## What is being ported (prokop commits 075f0122..85cf4f32)

### Optimizer machinery
- `cpp/common/math/DynamicOpt.{h,cpp}` — wholesale copy (prokop version is a superset:
  already contains Tomas's `bfixmask`, `dr_limit` uniform step scaling, `f2sum` returns).
  Adds: `autoMass` trajectory-secant diagonal preconditioner, `massMetric` mass-consistent
  FIRE metric, `modeSupp` stiff-mode suppression, gated two-timescale momentum,
  `run_GH` (BB1 / Anderson line-search optimizers), `saveState/restoreState` checkpoints,
  `ForceFunction` now returns energy (double).
- `cpp/common/math/Lingebra.{h,cpp}` — wholesale; adds `leastSquareFit_QR`.
- `cpp/common/math/fastmath.h` — wholesale; `erfx_e6` constrained refit
  (dy(0)=0, removes |r|-cusp at coincident charges).

### Physics bugfixes in `cpp/common/molecular/eFF.h` (manual merge — Tomas heavily edited)
- `echarge` is positive occupancy: `_realloc0(echarge, ne, 1.0)` (was QE=-1);
  assign `fabs(c)` in `readSingleFGO` (+ closed-shell copy).
- Internal Coulomb self-repulsion of paired `e2` orbital in `evalEE()`
  (`qi*(qi-1)/2`, missing ~32.5 eV per e2 pair => energies shift, rho2 fits may need re-check!).
- `evalAE`: gate `bEvalAEPauli` on `aPar.y>1e-8` (0/0 NaN for point nuclei).
  (fsi/fsj swap already fixed by Tomas; keep his verbose-free version.)
- `evalAE_MultiGauss` (Tomas's function): same swap bug in `bCoreCoul` term => fix arg order.
- `evalECP` + `evalAE_ECP`: ECP Coulomb used `aPar.z` (a charge) as Gaussian width and
  repulsive `+(Z-Zeff)`; fix to width `aPar.y`, charge `-qj*(aPar.x-aPar.z)`.
  ECP Pauli size/radial forces sign-reversed; fix `fs += fe_`, `f.add_mul(dR,-fr_/rc)`.
  Bare-nucleus branch: guard `aPar.z>0` (H has no core => 0/0 NaN).

### Library driver `cpp/libs/Molecular/eFF_lib.cpp` (manual merge — keep Tomas's API)
- `init_buffers()`: `buffers[k]=v` instead of `insert` (stale pointers after realloc).
- `run()`: keep Tomas's `bool* bConv` signature, but port:
  convergence on RAW free-DOF force tested BEFORE the move; fresh FIRE restart per call
  (cleanVel + damping=damp_max + lastNeg=0 + am_init=false); stale-buffer rebind after
  ff.realloc; isfinite aborts (exit(1)); esize floor -> `ff.min_esize`; `*bConv` set on
  conv + NULL-safe local fallback (ctypes callers pass only 6 args!).
- New API: `setAutoMass`, `setMassMetric`, `setOptParam`, `run_GH`, `_eff_eval_E`,
  branch-lab step API (`optBind/optInitOpt/optEval/optMoveFIRE/optNDOFs/optBuffer/
  optStateSize/optSaveState/optRestoreState`).
- `processXYZ`: process last conf at EOF (`while(fgets() || last-conf-pending)`).

### OpenCL `cpp/common_resources/cl/eFF.cl`
- ~~Port only `erfx_e6` refit coefficients~~ => **superseded**: user confirmed Tomas's
  physics is CPU-only, so `eFF.cl` was **wholesale-overwritten with prokop's version**
  (466 -> 1214 lines: `localMD` gains `fout`/`fixmask` args, plus `localMD_energy`,
  `evalEnergy`, `eval_density_grid`, `fit_density_fire`). GPU physics now matches the
  ported CPU code (e2 internal Coulomb, qe=2 occupancy, erfx refit); Tomas's
  rho2/MultiGauss additions never existed on GPU anyway.
- `pyBall/OCL/eFF_ocl.py` replaced with prokop's driver (required — kernel signature
  changed). It is an API superset: same `EFF_OCL` class, `load_xyzs`/`relax_systems`
  keep compatible signatures.
- `pyBall/OCL/OpenCLBase.py`: added `build_options` kwarg to `load_program` (needed by
  new driver's `-DEFF_EVAL_ENERGY`/`-DDBG_*` build flags; other OCL modules unaffected).
- Verified: `python -m pyBall.OCL.eFF_ocl` smoke test compiles on RTX 3090 and runs
  `localMD` on 2 test systems.

### Python `pyBall/eFF.py`
- Add bindings: `setAutoMass`, `setMassMetric`, `setOptParam`, `run_GH`.
- Fix `lib.run` argtypes for Tomas's extra `bConv` param (pass a real pointer;
  current 6-arg call on a 7-arg C function is latent UB).
- Library selection for A/B: `EFF_LIB` env var (default `eFF_lib`,
  `eFF_lib_origZ` = frozen original).

### Tests (new files, no conflicts)
- `tests/tEFF/fire_study_lib.py`, `tests/tEFF/run_fire_study.py`,
  `tests/tEFF/test_erfx_refit.py` (regression test).

## NOT ported
- MMFF-side optimizer plumbing from commit 0c9fb61d (MolWorld_sp3.h, NBFF.h, UFF.h,
  MMFF_lib.cpp, pyBall/MMFF.py) — unrelated to the eFF H2O test; new DynamicOpt is
  backward-compatible so they keep compiling.
- `doc/Topics/Learning_Optimization_Heuristics/*`, `web/`, Draw3D.js, audit docs.

## Build & test
- Build dir: `cd cpp && mkdir -p Build && cd Build && cmake .. -DRELEASE=ON`
  (mirrors prokop repo; eFF_lib does not need SDL/OCL).
- `run.sh` in `tests/tEFF` then does `make eFF_lib` inside `cpp/Build/libs/Molecular`.
- Backup: after first build, `cp libeFF_lib.so libeFF_lib_origZ.so`.
- H2O test: Tomas's `tests/tEFF` scripts (processXYZ_e / minimum.py / H2O_*.xyz files).

## Original (Tomas) backup — three independent layers

1. **Git tag `zverina-orig`** → marks Tomas's last commit (`201caa6a`); `HEAD` of this
   branch IS his tree since the port is not committed yet.
2. **Git worktree `/home/prokop/git/FireCore-origZ`** → full clean checkout of his
   original tree, buildable standalone (`mkdir Build && cmake .. && make eFF_lib`).
3. **In-branch side-by-side `cpp/libs/Molecular/origZ/`** → the working backup used for
   A/B testing: original `eFF_lib.cpp` + `molecular/eFF.h` + a *minimal* set of
   original `math/` headers (12 files). All of Tomas's physics (rho2, MultiGauss Pauli,
   ECP/`paramsSi`/`eConst` params, XYZ format) is preserved there.

### Why origZ/math/ is a minimal set (and not a full copy)

The shadow mechanism relies on `#include "file"` same-dir lookup: any common header on the
path to a *differing* header would pull the ported version. Only files that transitively
include a differing header (`fastmath.h`, `DynamicOpt.h`) need shadowing:
`Vec3.h`, `Vec2.h`, `Mat3.h`, `quaternion.h`, `Forces.h`, `VecN.h`, `Vec3Utils.h`,
`raytrace.h` + the differing `DynamicOpt.{h,cpp}`, `fastmath.h` themselves.
Pure-leaf headers (`gonioApprox.h`, `integerOps.h`, `physics_constants.h`) and
`Lingebra.*` (not in the origZ include chain) resolve to the identical common copies.
Verified: `gcc -MM` dep check + rebuild reproduces E=-897.8811 @2k it on dimer_90,
bit-identical to the pre-trim build.

### Build + switch
- `make eFF_lib` → `libeFF_lib.so` (ported); `make eFF_lib_origZ` → `libeFF_lib_origZ.so`
  (original) in the same Build dir (`origZ/` includes prepend BEFORE global common dirs).
- **Runtime A/B switch**: `EFF_LIB=eFF_lib` (default) vs `EFF_LIB=eFF_lib_origZ`.
- Verified identical: rebuilt origZ lib reproduces -510.81013377 eV / 358 it on
  H2O_pairs_fc — bit-identical to the previously frozen binary.

## Physics & behavior reference
- **`PHYSICS_DIFFS.md`** — detailed three-way per-term physics comparison
  (origZ vs prokop vs ported): rho2 Pauli, MultiGauss, evalAA, ECP signs, e2 term, params.
- **`tests/tEFF/tutorial_eFF_H2O.md`** — didactic report of the H2O/H-bond relaxation
  experiments incl. the dimer A/B between the two libs (§7.5); reproduced by
  `tests/tEFF/tutorial_H2O_relax.py` (outputs → `tests/tEFF/artifacts/`, gitignored).

## Known caveats to report
- e2 internal repulsion changes relaxed energies => Tomas's rho2/KRSrho fits are
  against the old Hamiltonian and may need refitting (that's why we keep origZ lib).
- `min_esize` = 0.0 in Zverina (prokop used 0.1) — kept Tomas's value.
- Python `lib.run` 6-arg call on 7-arg C++ fn = latent bug in Tomas's branch (fixed here).

---

## Test results (2026-05, this session)

Build: `make eFF_lib` in `cpp/Build/libs/Molecular` — clean, no warnings on ported code.
Release build (`-Ofast -march=native -O3 -DNDEBUG`); run.sh's `LD_PRELOAD=libasan.so`
does not instrument (lib has no `-fsanitize` flag) — same for orig baseline.

Test protocol: `pyBall/eFF.py` loads `libeFF_lib.so` (ported) or
`libeFF_lib_origZ.so` (frozen original) via `EFF_LIB` env var. Driver:
`eff.processXYZ_e(xyz, nstepMax, dt, Fconv=1e-3, optAlg=2)` (FIRE) +
`eff.setAutoMass(1,100)` / `eff.run_GH(...)` for new-lib-only features.

### H2O_pairs_fc.xyz (frozen core, 4 paired e2 orbitals — easy case)

| lib | iters | Etot [eV] | result |
|-----|-------|-----------|--------|
| origZ | 358 | -510.81 | converged |
| ported | 365 | -401.95 | converged |

Energy shift +108.86 eV ≈ 4 orbitals × ~27 eV internal e2 self-repulsion now
correctly included (expected physics fix, see eFF_optimizer_speedup.md).
Iteration count similar — this case is not stiff.

### H2O_spins.xyz (all-electron, 9 e incl. sz=0.2 core e — stiff case)

| config | iters | Etot [eV] | |F| end | result |
|--------|-------|-----------|------|--------|
| origZ, dt=0.002 | 30000 | -1915.75 | nan | NOT converged; trajectory diverged (nan DOF, esize=[nan, 5.9e240, ...]) |
| ported FIRE, dt=0.002 | 30000 | -1785.56 | 0.024 | NOT converged but stable, still descending |
| ported FIRE+autoMass(1,100), dt=0.002 | 30000 | -1787.77 | 9.15 | NOT converged (destabilized late) |
| ported FIRE+autoMass(1,100), **dt=0.01** | **7735** | **-1785.55** | 0.001 | **CONVERGED** — same minimum as slow FIRE |
| ported FIRE+autoMass(1,100), dt=0.05 | 30000 | -1497.02 | 2.4e4 | diverged (dt too big) |

=> autoMass preconditioner + larger dt converges the case that diverged on the
original library; ~4x fewer evals than the still-running unconditioned FIRE.

### run_GH on H2O_spins (autoMass=1, mHist=5, eta=0.01, a_init=0.01, adapt=0)

| iMethod | evals | E end | |F| end | status |
|---------|-------|-------|------|--------|
| 1 BB1 | 1249 | -1799.53 | 72 | LINE-SEARCH-FAIL (770 rejects) — found deeper min than FIRE! |
| 2 Anderson | 1626 | -1789.96 | 735 | LINE-SEARCH-FAIL (1371 rejects) |
| 3 Anderson+precond | 399 | -1780.22 | 78 | LINE-SEARCH-FAIL (352 rejects, 16 resets) |

All three run mechanically correctly (evals, bounds, history, no crash) but bail
out on this stiff system at eta=0.01 — needs parameter tuning (smaller eta,
adapt=1). Note BB1 wandered into a *lower* basin (-1799.5 < -1785.5).

### Three-way verification (origZ vs ported vs prokop-repo lib)

| lib | H2O_pairs_fc | H2O_spins |
|-----|--------------|-----------|
| origZ (rebuilt from origZ/ sources) | -510.8101 (358 it) | -1915.75, NaN divergence @30k |
| ported | -401.9505 (365 it) | -1785.5546 CONVERGED 7735 it (autoMass,dt=.01) |
| prokop repo (`~/git/FireCore`, self-built) | **-343.31951053** — matches eFF_optimizer_speedup.md documented value | -1719.54 @30k |

- origZ rebuilt source ≡ frozen binary (bit-identical) => origZ reproduces previous results.
- ported ≠ prokop-repo numbers BY DESIGN (keeps Tomas's Hamiltonian: rho2, point-nucleus
  evalAE, 4-term evalAA, O R_eff=0.1 vs 0.1678); the +108.9 eV origZ→ported shift =
  4×~27.2 eV internal e2 self-Coulomb — matches the analytic expectation of the fix.

### EMolecules/H2O/15.mol2 — 15-water H-bond chain (processMol2, free nuclei, 615 DOF)

| lib | iters | Etot | |F| end | H-bonds |
|-----|-------|------|------|---------|
| origZ, dt=0.01 | 30000 | -6737.81 | 54.1 | 11 |
| ported, dt=0.01 | 30000 | -6737.79 | 384.5 | 11 |
| ported+autoMass(1,100), dt=0.01 | 30000 | -6742.89 | 0.457 | 10 |

- **H-bond bridges DO form** (what Tomas was testing): chain self-organizes into
  O-O ~3.0 Å segments with each molecule's second H pointing to the next O at
  H···O ~2.17 Å. Both libs produce the same structure — physics preserved.
- Chain breaks into 3 segments of 5 (12.37 Å gaps at positions 4-5, 9-10) in all runs.
- autoMass improves |F| ~800x (385 -> 0.46) and finds slightly deeper E.
- **WARNING**: both libs (incl. origZ) crash AFTER output with
  `free(): invalid next size` / `double free or corruption` — pre-existing heap
  corruption in the processMol2/builder teardown path, NOT introduced by the port.

## Open items / follow-ups
1. `min_esize=0.0` in Zverina eFF.h (prokop: 0.1) — kept Tomas's 0.0, so the
   esize floor in run()/run_GH is effectively inactive. If size collapse is
   ever observed, set min_esize>0.
2. ~~GPU parity~~ RESOLVED: `eFF.cl` + `eFF_ocl.py` wholesale-replaced by prokop's
   versions (Tomas's physics was CPU-only). GPU now has the same fixed physics.
3. ECP paths (`evalECP`/`evalAE_ECP`) got prokop's sign/width fixes — Tomas's
   ECP fits (paramsSi, eConst) were fitted against the buggy version and may
   need refitting. Same for rho2/KRSrho vs the +27 eV/orbital e2 fix.
4. setAutoMass called before opt is bound (ff.na=ff.ne=0) creates a 0-length
   am_mass0 snapshot — harmless for use via processXYZ_e (run() re-derives
   nBlocks), minor edge case inherited from prokop's code.
5. New pyBall bindings are guarded by `_bNewOpt` so `EFF_LIB=eFF_lib_origZ`
   still imports; calling them on the old lib raises a clear error.
