# MMFFsp3_loc scalar-energy evaluation: bugs found & fixed (2026-11-08)

Triggered by energy-only optimizer benchmarking (`tests/tMMFF/run_valleyopt.py`). MMFF energy was never used before — FIRE is force-based — so the scalar-energy path had silently rotted. This document records every bug found, the fix, and the remaining caveats so they are not rediscovered.

**Verification gate:** `tests/tMMFF/test_energy_consistency.py` — checks `numgrad(E) == -fDOFs` and `eval() == eval_no_omp() == run() outE`. Run:

```bash
cd tests/tMMFF
LD_PRELOAD=$(gcc -print-file-name=libasan.so) ASAN_OPTIONS=detect_leaks=0 \
LD_LIBRARY_PATH=../../cpp/Build/libs/Molecular:$PWD/.. python3 test_energy_consistency.py
# 3/3 PASS: non-PBC rel ~7e-9; PBC has documented numgrad caveat (see below)
```

---

## DOF layout reference (source of several bugs)

Two different MMFF force-field objects coexist in `MolWorld_sp3`:

| object | class | `nvecs` | `DOFs` layout | dims fields |
|---|---|---|---|---|
| `W.ff`  | `MMFFsp3`     | `natoms + npi`   | `[apos(natoms) | pipos(npi)]`   | `nDOFs,natoms,nnode,ncap,npi,nbonds,nvecs,ne,ie0` (consecutive ints) |
| `W.ffl` | `MMFFsp3_loc` | `natoms + nnode` | `[apos(natoms) | pipos(nnode)]` | `nDOFs,nnode,ncap,nvecs,ntors` (different field order!) |

- `ffl.pipos`/`ffl.fpipos` are `[nnode]` — **every node atom gets a pi vector**, whether or not it has real pi orbitals.
- For all-sp3 molecules (propandiol, xylitol) `Ksp=Kpp=0` on all nodes → pi DOFs are **dead flat directions** (verified: E invariant under pi scramble). `run_valleyopt.py` masks them out via `ArrayEnergy(free_mask)`; `setup()` prints `active_pi`.
- `mmff.DOFs` aliases `ffl.DOFs`; `mmff.apos` ≡ `DOFs[:natoms]`; `mmff.pipos` ≡ `DOFs[natoms:]`.
- Pi energies (`Eps`, `EppT`, `EppI`) depend only on pi **direction**; `ffl.normalizePis()` projects to unit sphere.

## Bugs found & fixed

### 1. `ndims` bound to `W.ff`, buffers point to `W.ffl` — HIDDEN STATE (the "E drift" mystery)
- `MMFF_lib.cpp::init_buffers`: `ibuffers["ndims"] = &W.ff.nDOFs` while `buffers["DOFs"]=W.ffl.DOFs`.
- Python shaped `DOFs`/`pipos`/`fDOFs`/`vDOFs` from `W.ff` dims (e.g. nHexadecan: saw 162 of the real 210 DOFs) → `pipos[4..19]` invisible → `run()` mutated state Python could not restore → **energy changed at "identical" positions** (drift localized to pi terms, `Epi` 12.08→9.88).
- Fix: explicit `static int ndims[9]` populated from `ffl` when `bMMFF && !bUFF`; `'npi'` slot reports `nnode` (ffl's pipos size). `Es` rebound `&W.ff.Etot` → `&W.ffl.Etot` (same 6-double layout, now live).
- Symptoms seen: nDOFs 162→210 (nHexadecan), 66→96 (xylitol), 39→54 (propandiol) — before masking dead pi.

### 2. `MolWorld_sp3::eval()` — wrong nonbonded path + missing accumulation
- Called `evalLJQs_ng4_PBC_omp()` unconditionally (PBC variant even for non-PBC), and that function in `NBFF.h` never did `E +=` → returned 0 while writing NB forces → `numgrad(E) ≠ -fDOFs` by construction.
- Fix: `eval()` now mirrors `eval_no_omp()` branching on `bPBC`/`bNonBondNeighs`/`bExclusion2` against `ffl` (or `ffu` for UFF).

### 3. `setNonBondStrategy` ordering — first-eval garbage
- `eval()` called it *after* `ffl.eval()`; the bonded NB-subtraction inside `eval_atom` depends on those flags → first `eval()` returned garbage (e.g. −16983 instead of 0.176) and "fixed itself" after any other eval call. Moved before `ffl.eval()`, matching `eval_no_omp` order.

### 4. Per-atom nonbonded energy double-counted
- `evalLJQs_*_atom*` loops all `j` from each `ia` → each pair's `eij` summed at both endpoints → `E_reported = E_b + 2·E_NB` while forces were correct 1×. `run()`'s `outE` was likewise inflated (e.g. propandiol E0 reported 15905 instead of 0.176; the "glassy landscape where FIRE oscillates forever" was substantially this artifact).
- Fix: return `E*0.5` (standard half-share partition) in `NBFF.h` per-atom functions (`evalLJQs_ng4_atom`, `evalLJQs_ng4_atom_omp`, `evalLJQs_ng4_PBC_atom_omp`, `evalLJQs_atom_omp`, `evalLJQs_PBC_atom_omp`, `evalLJQs_ex2_atom`, `evalLJQs_atom_avx`).
- Whole-pair routines (`evalLJQs_ng4`, `evalLJQs`, `evalLJQs_PBC`) already count once — NOT halved.

### 5. `eval_no_omp`/`run_no_omp` never normalized pi vectors — unbounded unphysical DOF
- `ffl.eval()` calls `normalizePis()` (MMFFsp3_loc.h:1051); the per-atom loop paths did not → energy scaled with `|pipos|` → scipy Powell legitimately drove `|pipos|≈5·10⁴` → `E=-54716`.
- Fix: `ffl.normalizePis()` added at top of `eval_no_omp()` and `run_no_omp()`.
- **GPU precedent (already solved years ago):** `common_resources/cl/relax_multi.cl::updateAtomsMMFFf4` line ~1392 does `pe.xyz = normalize(pe.xyz)` **every integration step**, and the (commented) FIRE variant at ~1238-1240 additionally projects force and velocity onto the tangent plane of the unit sphere. The CPU path lacked both. Still open: port the tangent-projection (proper rotational dynamics) into `move_atom_FIRE` so CPU motion matches GPU exactly.

### 6. `nPBC=(0,0,0)` vs `lvs` in xyz — wrong PBC activation
- `MolWorld_sp3::loadGeom` copied `builder.bPBC`; an xyz `lvs` line forced `bPBC=true` even for explicit `(0,0,0)` → entered PBC paths with `neighCell` still `-1` → `shifts[-1]` wild read (ASan heap-use-after-free at `MMFFsp3_loc::eval_atom`, freed via `MMFFBuilder::assignSpecialTypes` scope-exit).
- Fix: explicit `(0,0,0)` now overrides lattice-derived PBC.
- **Still latent:** `neighCell` is populated only via `toMMFFsp3_loc` (from bond `ipbc`) — `makeNeighCells` is called **only for UFF**. Under *real* PBC with MMFF, `neighCell` correctness relies on builder bond assignment; do NOT assume it is valid for arbitrary PBC setups. Also `ffl.bindShifts` binds `ffl.shifts`→`W.pbc_shifts` which `_realloc` may free/realloc → stale-binding risk; verify before relying on MMFF+PBC.

### 7. `evalLJQs_ng4_PBC_omp` missing `E +=` — fixed (see bug 2).

## Remaining caveats / open issues

- **PBC force clamping:** `bClampNonBonded` caps close-contact `fij` at `FmaxNonBonded` while `eij` keeps the true LJ → `numgrad(E) ≠ -fDOFs` on clashed pairs under PBC (test shows rel≈0.8 at E≈240; non-PBC `ng4` path does not clamp → consistent to ~7e-9 after relaxation). Deliberate FF feature, not a bug — but it means PBC energy/force parity can't hold in the clashed regime.
- **Repeated `mmff.init()` in one process** frees C++ buffers while Python numpy views still point at them → stale-pointer UAF under ASan. Run each configuration in a separate process (test_energy_consistency.py does this via subprocess).
- **ASan + matplotlib:** `LD_PRELOAD=libasan.so` crashes matplotlib's `ft2font` (`__cxa_throw` interceptor). `run_valleyopt.py` saves `.npz` then replots in a clean child process (`--plot` flag; module imports are lazy for this reason).
- **`Es` decomposition:** bound to `ffl.Etot..EppI`; `Etot` is bonded+NB via separate accumulators — currently only `Es[0]` (bonded part) is reliably nonzero; the full per-term decomposition is not populated by the per-atom path.
- **UFF path:** `ibuffers["ndims"]` UFF case still reads `W.ff` fields (unchanged, original behavior); verify `init_buffers_UFF` separately if UFF benchmarking is needed.

## Benchmark metrics used (for interpreting `run_valleyopt.py` tables)

- `E0` = energy at the initial (obabel/clashed) geometry.
- `Eref` = minimum **visited** energy over a long FIRE run (50k iters; prints iter count — early exit = converged to `Fconv`).
- Column `E_inc`/`dE%`: best energy each method found within `neval` evals; `dE% = 100·(E_best − Eref)/(E0 − Eref)` = fraction of the start→reference gap remaining (negative = deeper than reference).
- Milestone columns (`0.5 0.1 0.02 0.005`): first eval index reaching `Eref + L·(E0 − Eref)` — i.e. **evals to relative-accuracy level** (the convergence-speed metric).
- `FIRE*` = FIRE best-visited at the same eval count — force-based, not information-equivalent to energy-only methods.
- Per-method curves are saved to `<mol>_convergence.npz` and plotted to `<mol>_convergence.png` (gitignored artifact; regenerate by rerunning the benchmark or `--plot <name>`). The plot shows `(E_best(n) − Eref)/(E0 − Eref)` vs evals on log-y — the actual convergence-rate curves.

## Post-fix benchmark snapshot (non-PBC, neval=2000, step=0.05, seed=3, dead pi DOFs masked)

| molecule | free DOFs | E0 | FIRE ref (iters) | FIRE* 2k | powell | secant | ridge |
|---|---:|---:|---:|---:|---:|---:|---:|
| propandiol | 39 | 0.176 | 0.1075 (18669, conv) | 0.9% | 13.4% | 27.8% | 35.2% |
| xylitol | 66 | 0.927 | 0.3338 (21669, conv) | 4.4% | 24.0% | 30.7% | 37.6% |
| nHexadecan-2COOH | 168 | 11.499 | −1.227 (50000, best) | 0.7% | 12.0% | 25.2% | 43.3% |

(dE% remaining gap; lower better.) With the corrected smooth energy FIRE converges on small molecules — the previous "rugged landscape" conclusion was an energy bug artifact.
