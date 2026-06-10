# ASAN follow-up: double-free vs ctypes, and `eval_atom` OOB

Companion to `Debug_ASan.md`. Two separate bugs.

**Ownership / debugging reference:** [Memory_Ownership_and_Deallocation.md](../../dev_notes/MMFF/Memory_Ownership_and_Deallocation.md) · skill: `.cursor/skills/cpp-memory-ownership/SKILL.md`

---

## 1) Double-free — is it Python GC / `getBuffs`?

**Short answer:** Mostly **no**. Python's garbage collector does **not** free the C++ buffers. `getBuffs()` builds **non-owning** NumPy views:

```python
ptr = lib.getBuff(name)
return np.ctypeslib.as_array(ptr, shape=sh)   # view into C++ memory; no free on GC
```

`init_buffers()` only fills a `std::map` of raw pointers into the global `MolWorld_sp3 W` (`W.ffl.apos`, `W.nbmol.fapos`, …). Python holds aliases; C++ owns allocation.

### Where the double-free actually comes from (C++)

Global singleton:

```cpp
// MMFF_lib.cpp
MolWorld_sp3 W;
```

After `makeFFs()` → `initNBmol(&ffl)`:

```cpp
// NBFF.h — nbmol.apos aliases ffl.apos (no copy)
nbmol.bindOrRealloc(ff->natoms, ff->apos, ff->fapos, ...);
```

On teardown or `clearFFs()`:

```cpp
ffl.dealloc();    // frees DOFs → apos/fapos
nbmol.dealloc();  // _dealloc(apos); _dealloc(fapos);  AGAIN on same pointers
```

`NBFF::dealloc()` always calls `_dealloc(apos)` even when `apos` was **borrowed** via `_bindOrRealloc`. That is the teardown **double-free** (seen as `free(): double free detected in tcache 2` on process exit).

### Second `MMFF.init()` in one process

`MMFF_lib::init()` does **not** call `W.clear()` first. It always runs:

```cpp
W.initParams(...);   // params.init() appends — no clear
W.init();            // makeFFs() reallocates without full teardown
```

Second `initParams` → `ElementType[22](E) is duplicated => Exit()`. Then corrupted heap / another double-free.

**`getBuffs` stale views:** If you `init()` again without `clear()`, C++ may `realloc` while Python still holds old `apos`/`fapos` pointers → use-after-free on **read/write**, not Python GC freeing memory.

### Robust fixes (recommended order)

| Fix | What |
|-----|------|
| **A. `W.clear()` at start of `MMFF_lib::init()`** | Reset builder, `clearFFs()`, `params.clear()` before rebuild |
| **B. `NBFF::dealloc()` respect aliasing** | If `bindOrRealloc` borrowed pointers, set `apos=fapos=0` before dealloc, or track `bOwnsApos` flag; only free owned buffers |
| **C. `clearFFs()` order** | Null `nbmol.apos/fapos` (or call `nbmol.dealloc()` **before** `ffl.dealloc()`, with nbmol not freeing aliased ptrs) |
| **D. Python contract** | One geometry per process **or** subprocess per `init()`; after re-init call `getBuffs()` again; never cache ctypes arrays across `init()` |
| **E. Optional `MMFF.shutdown()`** | Expose `W.clear()` to Python for explicit reset without unloading `.so` |

Subprocess-per-job (your phonon harness) is a valid **workaround**, not the root fix.

---

## 2) ASAN `eval_atom` — root cause (confirmed)

### ASAN message

```
heap-buffer-overflow in MMFFsp3_loc::eval_atom
READ 0 bytes after 32-byte region
```

Opaque because the bad access is `shifts[ipbc]` with **`ipbc = -1`** (one element before the array).

### Confirmed with debug print (2026-06-10)

```
neighCell{ -1, -1, -1, -1 }
ERROR MMFFsp3_loc::eval_atom ia=0 ing=1 bond=0: ipbc=-1 npbc=1 bPBC=1
```

### Inconsistent flags (diamond primitive, `nPBC=(0,0,0)`)

| Layer | Value | Effect |
|-------|-------|--------|
| XYZ has `lvs` | `builder.bPBC = true` | `autoBondsPBC()`, bonds tagged `pbc(-1,0,0)` … |
| Python `nPBC=(0,0,0)` | `W.nPBC = 0` | `makePBCshifts` → **npbc=1** (single shift) |
| `toMMFFsp3_loc(..., nPBC)` | only fills `neighCell` if `nPBC.x>0` | **neighCell stays -1** |
| `makeMMFFs` | `makeNeighCells()` **commented out** | indices never repaired |
| `ff.bPBC = builder.bPBC` | **true** | `eval_atom` uses `shifts[ingC[i]]` |

Code path:

```cpp
// MMFFsp3_loc.h eval_atom
int ipbc = ingC[i];           // -1
h.f.add( shifts[ipbc] );       // shifts[-1]  → ASAN OOB
```

Opt “works” because reading `shifts[-1]` is UB but may not crash; ASAN catches it in `eval_check()` during `init()`.

### Why phonon supercell uses `nPBC=(0,0,0)`

Explicit 3×3×3 geometry — no image summation. Lattice in XYZ is for cell metadata / bond topology, **not** for shift-table PBC during FD. Runtime PBC must be off unless `W.nPBC>0` and `neighCell` is valid.

### Fix options (pick one coherent policy)

1. **Runtime PBC gate (minimal):** In `toMMFFsp3_loc` or `makeMMFFs`:
   `ff.bPBC = builder.bPBC && (nPBC.x>0 || nPBC.y>0 || nPBC.z>0);`
   For phonon cluster path, `eval_atom` uses direct `apos[ing]-apos[ia]` (no shifts).

2. **Fill `neighCell`:** Call `ffl.makeNeighCells(npbc, pbc_shifts)` when `builder.bPBC` (was disabled for primitive cells). Needs correct `nPBC` for shift table size.

3. **Map bond `ipbc` when building:** Extend `toMMFFsp3_loc` to map `Bond.ipbc` → shift index even for `W.nPBC=(0,0,0)` (only works if shift table includes those cells).

4. **Guard in `eval_atom`:** If `ipbc<0 || ipbc>=npbc`, skip shift or use `invLvec.nearestCell` branch — defensive, does not fix wrong physics.

**Recommended for phonon work:** (1) + eventually (2) when `hessian_pbc=True` with real `nPBC=(1,1,1)`.

### How to test / bisect

```bash
cd tests/tMMFF
bash run_asan_minimal.sh asan skip-full   # fails step 1 until fix
python3 test_asan_minimal.py --step 1 --no-fast-exit   # teardown double-free
```

Temporary C++ guard (already validated):

```cpp
if(ipbc<0 || ipbc>=npbc) printf("ERROR ... ia=%i ipbc=%i npbc=%i ...\n", ...);
```

After fix, expect Build-asan step 1–6 PASS in `test_asan_minimal.py`.

### Files to change (when implementing)

| File | Change |
|------|--------|
| `MMFF_lib.cpp` | `W.clear()` before `W.init()` |
| `NBFF.h` | `dealloc()` skip aliased `apos`/`fapos` |
| `MolWorld_sp3.h` `makeMMFFs` | `ff.bPBC` vs `W.nPBC`; optionally enable `makeNeighCells` |
| `MMFFBuilder.h` `toMMFFsp3_loc` | `neighCell` mapping when lattice but `nPBC=0` |
| `MMFFsp3_loc.h` `eval_atom` | bounds check on `ipbc` (fail loudly) |


---

after implementing new memory ownership described in `/doc/dev_notes/MMFF/Memory_Ownership_and_Deallocation.md`

## Test results

Rebuilt `libMMFF_lib.so` on **Build-opt** and ran basic tests from `tests/tMMFF` and `tests/tUFF`.

### MMFF (`tests/tMMFF`)

| Test | Result | Notes |
|------|--------|-------|
| `bash run_asan_minimal.sh opt` | **ALL PASS** (steps 1–7) | prim init, eval, 6×6 Hessian, 3×3×3 supercell, 162×162 Hessian, double `init` |
| H2O `init → clear → init` (opt) | **PASS** | Same energy both passes |
| H2O `init → clear → init` (ASAN) | **PASS** | Ownership/teardown OK on ASAN build |
| H2O + benzene `eval()` | **PASS** | `E(H2O)=-5679.70`, `E(benzene)=-10419.56` |
| ASAN diamond primitive (`test_asan_minimal --step 1`) | **FAIL** | Heap OOB — known `ipbc=-1` / PBC bug, not ownership |
| `run_test_clear.py` | **Not run** | Missing `data/hydropentacene_cross.xyz` |

### UFF (`tests/tUFF`)

| Test | Result | Notes |
|------|--------|-------|
| `run.py` water, ethylene, benzene | **exit 0** | No crash; `clearFFs()` on teardown |
| Batch: water, ethane, benzene, adenine, guanine | **5/5 exit 0** | NaNs in angle/nonbond forces (UFF + `nPBC=(1,1,1)` physics, separate from memory) |

### How to run

**Recommended MMFF smoke test (ownership + Hessian):**
```bash
cd tests/tMMFF
bash run_asan_minimal.sh opt          # full suite, no matplotlib
bash run_asan_minimal.sh asan skip-full   # ASAN without 162×162 (diamond still fails on step 1)
```

**UFF basic eval:**
```bash
cd tests/tUFF
unset LD_PRELOAD ASAN_OPTIONS LSAN_OPTIONS
python3 run.py ./data_UFF/xyz/water.xyz
python3 run.py ./data_UFF/xyz/benzene.xyz
```

**Note on `run.sh`:** It builds in `Build-opt` but always sets `LD_PRELOAD=libasan` and runs `test_diamond_phonon_bands.py` — mixed opt/ASAN. Prefer `run_asan_minimal.sh` for library/memory validation.

### Takeaways

1. **Memory ownership fix looks solid** on opt and ASAN for H2O recycle and the full `test_asan_minimal` opt suite.
2. **Diamond primitive + ASAN** still fails on the PBC path (documented separately).
3. **UFF tests run without crashes** but show NaN forces — worth a separate physics/debug pass, not a blocker for the ownership work.
4. `Build` symlink reset to **Build-opt**.

I can next run a broader UFF batch via `run_firecore.sh`, fix `run_test_clear.py` to use `H2O`, or add a small `run_basic.sh` that runs the opt + UFF smoke tests in one command.