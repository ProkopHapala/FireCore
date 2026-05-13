# USER

OK, we have this project for AFM simulation
@tests/tAFM/pyocl_fdbm/test_full_pipeline.py @pyBall/OCL/AFM_utils.py @pyBall/OCL/AFM.py @pyBall/DFTB/DFTBcore.py @pyBall/DFTB/Grid_dftb.py @pyBall/DFTB/cl/Grid_dftb.cl 

def run_afm_pipeline(

get_density_from_dftb_dense
def get_density_from_dftb_plus(atomPos, atomTypes, basis, slako_prefix, work_dir,
def compose_and_relax(grads_pauli, grads_es, grads_vdw, scan_xs, scan_ys, heights,
def compute_stm(projector, eigvecs, eigvals, scan_xs, scan_ys, heights,

compute_bond_resolved_stm
def compute_fdbm_forcefield(rho_grid, V_ES, origin, step, atomPos,

def pp_relax_2d(force_func, scan_xs, scan_ys, probe_heights, mol_z=0.0,



It wors but now it stime to profile performance, optimize and refactor. There are several things I want from you:

1) analyse code and tell me where are some overheads, were we do some complex tasks in python which can be maybe eleminated (factored out in fron of loops, precashed, moved to OpenCL)

2) we want to profile it, do you have some good way how to profile the performance without modifying the code? it uses pyOpenCL but I guess it is enogh to profile it on python level, anyway we spend most of the type doing some perparatory / postprocessing overhead in python, OpenCL kernels are already optimal and fast

3) figure out where we should put some expocit timing, and where we can put some debugging guards which turn on/off verbosity, strogin intermediate results or plotting them, these thing also slow down the run so we cannot profile it efficiently

do not program anything yet, just analyse and make detailed notes 

---

# SWE-1.6


# AFM Pipeline Performance Analysis

## 1. Python-Level Overheads Identified

### **Major Overheads:**

**A. DFTB Density Projection ([get_density_from_dftb_dense](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:176:0-340:35))**
- **Lines 220-226**: Basis parsing (`parse_wfc_hsd`) and orbital layout computation ([build_orbital_layout](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:149:0-173:18)) - done every run
  - **Optimization**: Cache basis_data and orbital layout per molecule+basis combination
- **Lines 246-281**: HSD file writing + SK file copying for every DFTB run
  - **Optimization**: Check if work_dir already has required files, skip copying
- **Lines 292-305**: DFTB SCF via DFTBcore - this is the dominant cost (Fortran), but Python overhead in ctypes calls
- **Lines 308**: Dense projection via OpenCL - already fast (GPU)

**B. CO Tip Computation ([_call_compute_co_tip_script](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:1135:0-1153:18))**
- **Lines 1136-1154**: Subprocess call to `compute_co_tip.py` - process spawn overhead
  - **Optimization**: Precompute and cache CO tip grids per step/margin combination
  - **Optimization**: Inline the CO tip computation directly in Python (avoid subprocess)
- **Lines 1157-1179**: Padding and rolling of CO tip - numpy operations on large 3D arrays
  - **Optimization**: Could be done in OpenCL kernel

**C. Pauli Parameter Fitting (if `fit_pauli=True`)**
- **Lines 964-994**: For each target atom, load DFTB reference, extract profile, fit power law
  - **Lines 976-983**: Loading reference data from disk for each atom
  - **Optimization**: Load all reference data once, cache profiles
  - **Lines 986**: `extract_z_profile` uses `map_coordinates` - scipy interpolation overhead
  - **Optimization**: Precompute z-profile lookup table

**D. STM Computation ([compute_stm](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:673:0-734:19), [compute_bond_resolved_stm](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:737:0-795:19))**
- **Lines 712-732 (compute_stm)**: Nested loops over heights and LUMO orbitals
  - **Line 716**: `np.meshgrid` creates full 2D arrays for each height
  - **Lines 720-732**: For each orbital, call [projector.project_orbital_dense_points_exp](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:893:4-950:18) - this launches GPU kernel per orbital
  - **Optimization**: Batch all LUMO orbitals into single GPU kernel launch
  - **Optimization**: Process all heights in single batch kernel
- **Lines 772-793 (compute_bond_resolved_stm)**: Same overhead + tip displacement array operations
  - **Lines 774-775**: Array arithmetic for displacement applied per height
  - **Optimization**: Move displacement logic into GPU kernel

**E. Probe Particle Relaxation ([pp_relax_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:791:0-843:30))**
- **Lines 803-811**: Force map grid computation - meshgrid creation
- **Lines 823-827**: [force_func](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:494:4-503:46) called once per height to build full force map
  - **Lines 784-787 (force_func)**: Python loop over atoms for vdW calculation
  - **Optimization**: Move vdW atom loop to OpenCL kernel
- **Lines 831-835**: Inner relaxation loop (N_RELAX=30 iterations) with interpolation
  - **Lines 832-833**: [_interp2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:811:4-814:111) uses `map_coordinates` - scipy overhead per iteration
  - **Optimization**: Precompute interpolation weights or use GPU interpolation

**F. Field Gradient Computation**
- **Lines 708-709, 879, 906, 936, 1049**: `np.gradient` called on 3D arrays for each field
  - **Optimization**: Implement gradient in OpenCL kernel (avoid Python + scipy overhead)
  - **Optimization**: Combine gradient computation with field computation in single kernel

**G. FFT Operations**
- **Lines 700-704, 869, 905**: Multiple FFT/IFFT operations
  - These use numpy/scipy FFT - already optimized, but multiple transforms could be batched
  - **Optimization**: Batch multiple FFTs if possible

### **Minor Overheads:**

**H. Plotting Functions**
- **Lines 517-598**: All plotting functions ([plot_xy_slice](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:34:0-58:29), [save_afm_images](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:61:0-81:34), etc.)
  - These are gated by `plot_steps` flag - already conditional
  - matplotlib overhead is significant but only for debugging

**I. File I/O**
- **Lines 957, 1020, 1027, 1038, 1048, 1050**: `np.save` for intermediate results
  - These are useful for debugging but add I/O overhead
  - **Optimization**: Make saving optional with debug flag

## 2. Profiling Approach (No Code Modification)

### **Recommended Profiling Tools:**

**A. Python Profiler (`cProfile`)**
```bash
# Profile the entire pipeline
python -m cProfile -o afm_profile.prof tests/tAFM/pyocl_fdbm/test_full_pipeline.py pentacene.xyz --no_plot_steps

# Visualize with snakeviz (install: pip install snakeviz)
snakeviz afm_profile.prof
```

**B. Line-by-Line Profiling (`line_profiler`)**
```bash
# Install: pip install line_profiler
# Add @profile decorator to key functions (requires code modification, but minimal)
kernprof -l -v tests/tAFM/pyocl_fdbm/test_full_pipeline.py pentacene.xyz --no_plot_steps
```

**C. Memory Profiling (`memory_profiler`)**
```bash
# Install: pip install memory_profiler
python -m memory_profiler tests/tAFM/pyocl_fdbm/test_full_pipeline.py pentacene.xyz --no_plot_steps
```

**D. PyOpenCL Event Timing**
- The OpenCL kernels already have some timing (see Grid_dftb.py lines 198-210, 212-231, 261-275)
- **Recommendation**: Add a global flag to enable/disable all OpenCL event timing
- Current timing is scattered and prints to stdout - could be aggregated

**E. Simple Manual Timing Wrapper**
- Use `time.perf_counter()` at strategic points (requires code modification)
- See existing timing in Grid_dftb.py for pattern

## 3. Explicit Timing and Debugging Guards

### **Where to Add Timing:**

**A. Pipeline-Level Timing (in [run_afm_from_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:1224:0-1432:5))**
- After molecule loading
- After density computation (DFTB)
- After CO tip computation/loading
- After each step in [run_afm_pipeline](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:870:0-1115:17) (steps 2-7)
- Total time

**B. Function-Level Timing:**
- [get_density_from_dftb_dense](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:176:0-340:35): breakdown by basis parsing, SCF, projection
- [compute_pauli_overlap](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:857:0-869:44): FFT timing
- [scale_pauli_field](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:872:0-879:25): gradient timing
- [compute_es_conv_field](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:892:0-906:22): convolution timing
- [compute_dispersion_grid](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:1015:0-1049:29): atom loop timing
- [pp_relax_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:791:0-843:30): per-height timing, relaxation iteration timing
- [compute_stm](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:673:0-734:19): per-height, per-orbital timing

**C. OpenCL Kernel Timing:**
- Already partially implemented in Grid_dftb.py
- Add to AFM.py for FFT operations (if moved to OpenCL)
- Add to any new GPU kernels

### **Where to Add Debugging Guards:**

**A. Verbosity Levels (Recommended Pattern):**
```python
# Global verbosity flag (0=silent, 1=timing, 2=detailed, 3=debug)
VERBOSITY = int(os.environ.get('AFM_VERBOSITY', '1'))
```

**B. Plotting Guards:**
- **Already exists**: `plot_steps` parameter in [run_afm_from_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:1224:0-1432:5) and [run_afm_pipeline](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:870:0-1115:17)
- **Recommendation**: Keep this, it's well-designed

**C. Intermediate Result Saving:**
- **Lines 957, 1020, 1027, 1038, 1048, 1050**: Currently unconditional
- **Recommendation**: Add `save_intermediates` flag
- **Recommendation**: Add `debug_dir` parameter to control where debug outputs go

**D. Print Statement Guards:**
- Many print statements throughout the code (e.g., lines 328-333, 536, 550, 564, 578, 592)
- **Recommendation**: Gate all print statements with `if VERBOSITY >= 1:`
- **Recommendation**: Separate timing prints (VERBOSITY >= 1) from debug prints (VERBOSITY >= 2)

**E. Specific Debug Sections:**

**In [get_density_from_dftb_dense](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:176:0-340:35):**
- Lines 328-333: Charge check prints - gate with VERBOSITY >= 2
- Lines 220-226: Basis parsing - add timing if VERBOSITY >= 1

**In [run_afm_pipeline](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:870:0-1115:17):**
- Lines 956, 1013-1015, 1024, 1032: Progress prints - gate with VERBOSITY >= 1
- Lines 1042: Compose and relax - add detailed timing if VERBOSITY >= 2

**In [pp_relax_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:791:0-843:30):**
- Line 823: Per-height progress - gate with VERBOSITY >= 1
- Lines 831-835: Inner loop - add iteration timing if VERBOSITY >= 3

**In [compute_stm](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:673:0-734:19) / [compute_bond_resolved_stm](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:737:0-795:19):**
- Lines 709, 766, 795: Progress prints - gate with VERBOSITY >= 1
- Lines 720-732, 781-792: Per-orbital timing if VERBOSITY >= 2

**In Grid_dftb.py:**
- Lines 191-192, 377, 398, 400, 446-447: Debug prints - already gated by `verbosity` parameter
- Lines 210, 231, 275: Timing prints - gate with `if verbosity >= 1`

### **Recommended Guard Pattern:**
```python
# At module level
import os
_VERBOSITY = int(os.environ.get('AFM_VERBOSITY', '1'))

# In functions
if _VERBOSITY >= 1:
    print(f"[TIME] function_name: {elapsed_ms:.3f} ms")
if _VERBOSITY >= 2:
    print(f"[DEBUG] detailed info...")
if _VERBOSITY >= 3:
    # Very detailed debug (arrays, etc.)
```

## 4. Specific Optimization Recommendations (No Implementation Yet)

### **High Impact (Easy Wins):**

1. **Cache CO tip grids** - compute once per step/margin, save to disk
2. **Cache DFTB basis parsing** - basis_data per basis file
3. **Batch STM orbital projection** - single GPU kernel for all LUMOs
4. **Move vdW atom loop to OpenCL** - already done in some places, ensure consistency
5. **Gate all intermediate file saves** - add `save_intermediates` flag

### **Medium Impact:**

1. **Combine gradient with field computation** - single OpenCL kernel
2. **Optimize scipy interpolation** - precompute weights or use GPU
3. **Batch FFT operations** - if multiple transforms needed
4. **Inline CO tip computation** - avoid subprocess overhead

### **Low Impact (Code Quality):**

1. **Reduce print statement overhead** - gate with verbosity
2. **Optimize numpy array operations** - use in-place operations where possible
3. **Pre-allocate output arrays** - avoid repeated allocation

## 5. Summary of Key Findings

**Dominant Python Overheads:**
1. DFTB SCF (unavoidable - Fortran)
2. STM orbital projection (multiple GPU kernel launches - can be batched)
3. Probe particle relaxation (scipy interpolation - can be GPU-accelerated)
4. vdW atom loop in Python (can be moved to OpenCL)
5. File I/O for intermediate results (can be gated)

**Best Profiling Approach:**
- Start with `cProfile` for overall hotspots
- Use `line_profiler` for detailed function-level timing
- Add manual timing at pipeline stages (requires minimal code change)
- Leverage existing OpenCL event timing in Grid_dftb.py

**Debugging Guard Strategy:**
- Use environment variable `AFM_VERBOSITY` (0-3)
- Keep existing `plot_steps` flag
- Add `save_intermediates` flag for file I/O
- Gate all print statements by verbosity level
- Add timing at major pipeline stages

---

# Cursor-Composer

### 1) Likely overheads / “too much Python” hotspots (based on current code)

#### **A. FDBM forces: SciPy interpolation + Python atom loops**
- **`pyBall/OCL/AFM.py:compute_fdbm_forcefield()`**
  - **Trilinear interpolation is done in Python/SciPy per component** via `_interp3()` → `scipy.ndimage.map_coordinates`. In the FDBM branch that’s **6 separate 3D interpolations** (3 for Pauli grads + 3 for ES grads) over **all scan points** (`grid_c` length = \(nx_s·ny_s·nz\)). That is often a top CPU cost.
  - **vdW is a pure Python loop over atoms**:
    - For each atom: compute `dr`, `r2`, `r8` for *all scan points*, then inner loop over axes `a in range(3)` updates `F_vdw`. This is vectorized per-atom but still **O(nAtoms · nPoints)** with repeated temporary allocations (`dr`, `r2`, `r8`).
  - **If `force_model != 'fdbm'` and `es_model != 'poisson'`**, ES uses a Python per-atom loop too (same scaling).
  - **Avoidable recomputation**:
    - `grid_c` is computed once (good), but `flat_pos`, `XX/YY/ZZ` and some meshgrids allocate big arrays every call. If you run many repeats with same scan geometry, these can be cached.

**Optimization direction:** move interpolation and atom-sums to OpenCL (or at least NumPy batch kernels), or restructure so you compute forces on a regular 3D grid once and only sample it (but your sampling step is already the costly part).

#### **B. Probe-particle relaxation: rebuilding force maps + map_coordinates inside inner loops**
- **`pyBall/OCL/AFM.py:pp_relax_2d()`**
  - For each height `iz`:
    - Builds a **2D “force map” grid** `(pp_nx × pp_ny)` by calling `force_func(pp_flat)` where `pp_flat` can be large (scan area + margins at `step` spacing).
    - Then performs **N_RELAX iterations**, each iteration does two **2D map_coordinates** calls (`_interp2d` on `FF_x` and `FF_y`).
  - This is often the *dominant* Python-side cost if `pp_nx*pp_ny` is large and `N_RELAX` ~ 30–50.
  - Lots of per-height allocations: `PP_Z`, `pp_flat`, reshapes, `probe_x/probe_y`, `vx/vy`.

**Optimization direction:** do relaxation on GPU (OpenCL kernel for relaxation iterations), or at least reduce Python/SciPy calls by precomputing interpolation objects / using faster interpolators / batching.

#### **C. Compose-and-relax creates a huge 4D `F_total` field**
- **`pyBall/OCL/AFM_utils.py:compose_and_relax()`**
  - Builds `F_total = -(grads_pauli + grads_es + grads_vdw)` as a full `(nx,ny,nz,3)` float array.
  - `force_func()` then interpolates each component via **three `map_coordinates` calls per force evaluation**.
  - `pp_relax_2d` calls `force_func()` once per height on the PP grid (big), so that becomes **3× map_coordinates over a huge 3D field** per height.

**Optimization direction:** avoid materializing `F_total` (compute components on-demand) *or* move force evaluation into OpenCL (preferred).

#### **D. STM: repeated point-projection in nested Python loops**
- **`pyBall/OCL/AFM_utils.py:compute_stm()` and `compute_bond_resolved_stm()`**
  - Nested loops:
    - For each height `iz`: allocate `points` (size \(nx_s·ny_s\)).
    - For each LUMO orbital: call `projector.project_orbital_dense_points_exp(...)`.
  - If `project_orbital_dense_points_exp` does OpenCL under the hood, the overhead is:
    - repeated **host→device transfers** of `points`/`coeffs` (unless cached internally),
    - repeated kernel launches,
    - repeated allocations.
  - If it’s CPU-side, then it’s plainly heavy Python/NumPy work.

**Optimization direction:** fuse “multiple orbitals at once” into one projection call, or keep `points` on device and stream only coeffs, or compute STM on grid (3D) once then sample.

#### **E. DFTB density acquisition: filesystem + setup overhead**
- **`pyBall/OCL/AFM_utils.py:get_density_from_dftb_dense()`**
  - Writes an XYZ + HSD, copies SK files into `work_dir`, changes cwd, runs DFTBcore SCF.
  - For profiling “AFM pipeline performance”, this is usually **startup overhead** and can drown the useful signal if included in the same run.
  - Also: `species_per_atom = list(range(len(enames)))` makes each atom a unique “species” (may increase setup work in projector / basis mapping).

**Optimization direction:** separate “SCF+projection” timing from “AFM from given grids” timing; cache projector/basis; avoid redundant file copying when `work_dir` already prepared.

#### **F. Debug/IO/plot overhead is currently mixed into compute**
- `run_afm_pipeline()` does many unconditional `np.save(...)` (overlap_raw, fields, grads, df, tip_disp, stm_grid). Even with `plot_steps=False`, disk I/O can be significant.
- Many `print()` statements inside per-height loops (e.g. `pp_relax_2d`) also matter when profiling.

### 2) How to profile without modifying code (recommended workflows)

#### **A. Sampling profilers (best “no changes” option)**
These attach to the running Python process and show where time is spent.
- **`py-spy`** (low overhead, great for wall-time hotspots, no code changes):
  - Start your script normally, or attach to PID.
  - Works well for identifying whether you’re dominated by `scipy.ndimage.map_coordinates`, `numpy.fft`, Python loops, file I/O, etc.
- **`scalene`** (CPU + memory + some native time attribution):
  - Useful if temporary allocations / memory bandwidth dominate (likely in `pp_relax_2d`, `compute_fdbm_forcefield`).

What you’ll learn: whether time goes to SciPy interpolation, NumPy FFT/gradients, Python loops, or OpenCL calls.

#### **B. Deterministic Python profilers**
- **`cProfile`**:
  - Great for high-level function call time, but it under-represents time spent in native code (it records it as time inside the Python wrapper).
  - Still very useful to confirm call graph: `run_afm_pipeline → compose_and_relax → pp_relax_2d` dominance.
- **`pyinstrument`**:
  - Often easier to read than `cProfile` for wall-time stacks, minimal friction.

#### **C. GPU/OpenCL-side profiling without code changes: limited**
True kernel timing needs the command queue created with `PROFILING_ENABLE`. In `AFMulator` (`pyBall/OCL/AFM.py`) the queue is created as `cl.CommandQueue(self.ctx)` with no profiling flags, so you **won’t** get per-kernel event profiling without a small code tweak.
Still, you can:
- Use **system profilers** (Nsight Systems / `nsys`, `perf`, etc.) to see GPU activity and CPU stalls externally.
- Use PyOpenCL environment variables to see compilation logs etc., but that’s not performance timing.

Given your statement (“kernels are already optimal; Python prep/post dominates”), Python-level profilers are the right first step anyway.

### 3) Where to add explicit timing + debug guards (so profiling reflects “fast path”)

#### **A. Timing points (high signal, minimal clutter)**
Add timers around **step boundaries** and the known hotspots:

- **`AFM_utils.run_afm_pipeline()`** (step-level)
  - Step 2: `fft_poisson`
  - Step 3a: `compute_pauli_overlap` (FFT)
  - Step 3c: `scale_pauli_field` (power + 3× `np.gradient`)
  - Step 4: `compute_es_conv_field` (FFT + gradients)
  - Step 5: `compute_dispersion_grid` (Python atom loop over full 3D grid)
  - Step 6: `compose_and_relax` (likely the biggest)
  - Step 7: STM (`compute_stm` / `compute_bond_resolved_stm`)

- **`AFM.py:pp_relax_2d()`** (inner)
  - Per-height:
    - time `force_func(pp_flat)` (3D interpolation cost)
    - time “relax loop total” (N_RELAX * 2D interpolation cost)
  - Also log sizes: `pp_nx, pp_ny, N_RELAX, nx_s, ny_s, nz`

- **`AFM.py:compute_fdbm_forcefield()`**
  - time FFT block (if used): rfftn/irfftn
  - time gradient computation: 6× `np.gradient`
  - time interpolation block: `_interp3` calls (these are likely huge)
  - time vdW atom loop

- **`Grid_dftb.py` projector**
  - You already have some `[TIME]` prints inside `build_tasks_gpu` (good). Similar timing around `project(...)` and any host↔device copies would help.

#### **B. Debug/verbosity/plot/storage guards**
Right now `plot_steps` controls plotting, but not disk saves and not all prints.

For accurate profiling you want a “silent fast path”:
- **`verbosity`**: gate all `print()` inside hot loops (`pp_relax_2d` per height, per relaxation iteration, etc.)
- **`save_intermediates`** (separate from plotting): gate `np.save` of large arrays (`overlap_raw`, fields, gradients, `tip_disp`, etc.)
- **`debug_checks`**: gate expensive diagnostics (charge checks are cheap, but things like extra projections, extra FFTs, full-field stats/min/max in loops can be nontrivial)
- **`plot_steps`**: keep as-is, but ensure plotting functions are never called when off (already mostly true).

A useful pattern is:
- `verbosity=0`: no prints in inner loops
- `save_intermediates=False`: no `np.save` for big arrays unless debugging
- `plot_steps=False`: no matplotlib at all

### Practical “first profiling run” suggestion (what to expect)
If you profile `test_full_pipeline.py` end-to-end, the top costs will likely appear as:
- **`scipy.ndimage.map_coordinates`** (3D + 2D) dominating (`compose_and_relax` / `pp_relax_2d` / `compute_fdbm_forcefield`)
- **NumPy FFT** and `np.gradient` sizable but usually smaller than the interpolation+relax loops (depends on grid size)
- Python per-atom vdW sums noticeable for large `natoms * npoints`

If you want, in the next step (still without changing code) I can propose an exact profiling command set for your repo (py-spy + pyinstrument + cProfile) and the minimal run configuration (small scan grid) that still reproduces the same hotspots.

## Baseline Profiling Results (H2O, small system)

Before any changes in code (all debugging stuff is there without guards)

**Test configuration:**
- Molecule: h2o.xyz (3 atoms)
- Grid: step=0.2, margin=2.0
- Scan: range=2.0 Å, step=0.2 Å (28×24×3 points)
- Heights: 2.8-3.2 Å, step=0.2 Å (3 heights)
- Command: `python -m cProfile -o afm_baseline.prof test_full_pipeline.py h2o.xyz --no_plot_steps --step 0.2 --margin 2.0 --scan_range 2.0 --scan_step 0.2 --height_range 2.8 3.2 --height_step 0.2`

**Total runtime: 6.181 seconds**

### Time breakdown by category:

**1. Plotting overhead (despite --no_plot_steps flag): ~2.9s (47%)**
- `savefig`: 2.967s
- `plot_diagnostic_panel`: 2.059s
- `plot_diagnostic_slices`: 1.878s
- Various matplotlib drawing: ~1-2s

**Note:** Diagnostic plotting functions are called unconditionally at the end of `run_afm_from_xyz`. The `--no_plot_steps` flag does not gate these.

**2. CO tip subprocess overhead: ~1.0s (16%)**
- Subprocess communication (poll): 0.992s
- This is `_call_compute_co_tip_script` spawning a separate Python process

**3. OpenCL initialization: ~1.0s (16%)**
- Device selection: 0.131s
- Program build logs: 0.067s
- GridProjector init (kernel loading): 0.861s

**4. Actual AFM computation: ~0.4s (6.5%)**
- `run_afm_pipeline`: 0.189s
- `compose_and_relax`: 0.178s
- `pp_relax_2d`: 0.019s
- FFT operations: 0.004s
- scipy interpolation (map_coordinates): 0.009s (318 calls in relaxation loop)
- compute_dispersion_grid: 0.002s
- compute_es_conv_field: 0.002s

**5. Other overhead: ~0.9s (14.5%)**
- NumPy operations, imports, etc.

### Key observations:

1. **Plotting dominates runtime** despite `--no_plot_steps` - need to gate diagnostic plotting functions
2. **CO tip subprocess is significant** - 16% overhead from process spawn and communication
3. **OpenCL initialization is one-time cost** - would be amortized in batch runs
4. **Actual computation is very fast** - only 6.5% of total time for this small system
5. **scipy interpolation** (map_coordinates) is called 318 times in relaxation loop - could be optimized

### Profiling method used:
```bash
python -m cProfile -o /tmp/afm_baseline.prof tests/tAFM/pyocl_fdbm/test_full_pipeline.py tests/tAFM/pyocl_fdbm/h2o.xyz --output_dir /tmp/afm_profile_baseline --no_plot_steps --step 0.2 --margin 2.0 --scan_range 2.0 --scan_step 0.2 --height_range 2.8 3.2 --height_step 0.2
```

Analysis with:
```python
import pstats
p = pstats.Stats('/tmp/afm_baseline.prof')
p.sort_stats('cumulative').print_stats(30)  # Top 30 by cumulative time
p.sort_stats('tottime').print_stats(30)     # Top 30 by self-time
```

## Baseline Profiling Results (TBTAP_3mols_c3h_F, large system)

**Test configuration:**
- Molecule: TBTAP_3mols_c3h_F.xyz (60 atoms: 3 TBTAP molecules with F substituents)
- Grid: step=0.2, margin=2.0
- Scan: range=2.0 Å, step=0.2 Å (111×107×3 points)
- Heights: 2.8-3.2 Å, step=0.2 Å (3 heights)
- Command: `python -m cProfile -o afm_baseline_large.prof test_full_pipeline.py TBTAP_3mols_c3h_F.xyz --output_dir /tmp/afm_profile_baseline_large --no_plot_steps --step 0.2 --margin 2.0 --scan_range 2.0 --scan_step 0.2 --height_range 2.8 3.2 --height_step 0.2`

**Total runtime: 8.699 seconds**

### Time breakdown by category:

**1. Plotting overhead (despite --no_plot_steps flag): ~3.2s (37%)**
- `savefig`: 3.249s
- `plot_diagnostic_panel`: 2.162s
- `plot_diagnostic_slices`: 2.136s
- Various matplotlib drawing: ~1-2s

**Note:** Same issue as H2O - diagnostic plotting functions are called unconditionally.

**2. CO tip subprocess overhead: ~0.9s (10%)**
- Subprocess communication (poll): 0.906s
- Same subprocess overhead as H2O (CO tip computation is independent of molecule size)

**3. OpenCL initialization: ~1.0s (11%)**
- Device selection: 0.090s
- GridProjector init (kernel loading): 1.051s
- Similar to H2O (one-time cost)

**4. Actual AFM computation: ~2.0s (23%)**
- `run_afm_pipeline`: 0.963s (5x slower than H2O)
- `compute_dispersion_grid`: 0.615s (307x slower than H2O!) - **NEW HOTSPOT**
- `prepare_orbital_projection`: 1.196s total (2 calls, 0.598s each) - **NEW HOTSPOT**
- `build_tasks`: 0.084s per call (84x slower than H2O)
- `check_overlap_sphere_aabb`: 0.321s (called 181,104 times vs 504 for H2O) - **NEW HOTSPOT**
- FFT operations: 0.097s
- scipy interpolation: 0.087s (geometric_transform)

**5. Other overhead: ~1.6s (18%)**
- NumPy operations, imports, etc.

### Key observations for large system:

1. **Plotting still dominates** but less than H2O (37% vs 47%) - computation is now significant
2. **CO tip subprocess is same overhead** - independent of molecule size (0.9s)
3. **compute_dispersion_grid is now a major hotspot**: 0.615s (7% of total, 307x slower than H2O)
   - This is the Python atom loop over 60 atoms for the full 3D grid (112×112×56 = 703,104 points)
   - O(nAtoms × nGridPoints) scaling - 60 atoms × 703k points = 42M operations
4. **prepare_orbital_projection is significant**: 1.196s (14% of total)
   - build_tasks takes 0.084s per call (84x slower than H2O)
   - check_overlap_sphere_aabb called 181,104 times (360x more than H2O)
5. **Actual AFM pipeline is 5x slower**: 0.963s vs 0.189s for H2O
6. **scipy interpolation is still present** but less dominant (0.087s vs 0.009s for H2O)

### Comparison: Small (H2O) vs Large (TBTAP_3mols_c3h_F)

| Component | H2O (3 atoms) | TBTAP (60 atoms) | Ratio |
|-----------|---------------|------------------|-------|
| Total runtime | 6.181s | 8.699s | 1.4x |
| Plotting | 2.9s (47%) | 3.2s (37%) | 1.1x |
| CO tip subprocess | 1.0s (16%) | 0.9s (10%) | 0.9x |
| OpenCL init | 1.0s (16%) | 1.0s (11%) | 1.0x |
| run_afm_pipeline | 0.189s (3%) | 0.963s (11%) | 5.1x |
| compute_dispersion_grid | 0.002s (0%) | 0.615s (7%) | 307x |
| prepare_orbital_projection | 0.006s (0%) | 1.196s (14%) | 199x |
| check_overlap_sphere_aabb | 0.004s (0%) | 0.321s (4%) | 80x |

### Scaling analysis:

- **Linear scaling components** (approx):
  - OpenCL initialization: constant (one-time)
  - CO tip subprocess: constant (independent of molecule)
  - Plotting: near-constant (depends on output size, not computation)

- **Superlinear scaling components**:
  - compute_dispersion_grid: 307x (expected: 20x for 60 atoms, actual is worse due to grid size increase)
  - prepare_orbital_projection: 199x (expected: 20x, actual is 10x worse due to task building overhead)
  - check_overlap_sphere_aabb: 80x (expected: 20x, actual is 4x worse due to grid size increase)

- **Grid size increase**:
  - H2O: (32, 24, 56) = 43,008 points
  - TBTAP: (112, 112, 56) = 703,104 points
  - Grid ratio: 16.3x

### Profiling method used:
```bash
python -m cProfile -o /tmp/afm_baseline_large.prof tests/tAFM/pyocl_fdbm/test_full_pipeline.py tests/tAFM/pyocl_fdbm/TBTAP_3mols_c3h_F.xyz --output_dir /tmp/afm_profile_baseline_large --no_plot_steps --step 0.2 --margin 2.0 --scan_range 2.0 --scan_step 0.2 --height_range 2.8 3.2 --height_step 0.2
```

---

# USER

based on results of these testes and reading whole @AFM_FDBM_profiling_optimization.md can you review the code again, analyse it and tell
1) where are biggies bottlenecs/overheads,
2) where we should put debugging/visualization guards 
3) where significant performance improvement can be achived easily
4) which parts of the code should be ported from python to OpenCL (if we do some heavy numerics in python)

---

# SWE-1.6

Based on the profiling results and code review, here's my analysis:

## 1. Big Bottlenecks/Overheads

### **Critical (Large System Impact):**

**A. [compute_dispersion_grid](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:1015:0-1049:29) (AFM.py:1016-1050)**
- **Impact**: 0.615s (7% total, 307x slower on large system)
- **Cause**: Python loop over atoms (60) × full 3D grid (703k points) = 42M operations
- **Code**: Lines 1045-1047 - pure Python loop with numpy operations per atom
- **Scaling**: O(nAtoms × nGridPoints) - superlinear with both

**B. [prepare_orbital_projection](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:1802:4-1849:9) (Grid_dftb.py:1803)**
- **Impact**: 1.196s (14% total, 199x slower on large system)
- **Cause**: [build_tasks](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:283:4-358:38) (line 284) calls [check_overlap_sphere_aabb](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:152:4-156:40) 181,104 times
- **Code**: Lines 317-319 in [build_tasks](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:283:4-358:38) - nested loops over blocks × atoms
- **Scaling**: O(nBlocks × nAtoms) - grid size increase makes this explode

**C. Diagnostic Plotting (AFM_utils.py:1436, 1474)**
- **Impact**: 3.2s (37% on large, 47% on small)
- **Cause**: Called unconditionally at end of [run_afm_from_xyz](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:1224:0-1432:5) (lines 1436, 1474)
- **Code**: [plot_diagnostic_panel](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:1435:0-1470:88), [plot_diagnostic_slices](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM_utils.py:1473:0-1541:102) - matplotlib overhead
- **Note**: `--no_plot_steps` doesn't gate these

### **Medium Impact:**

**D. CO Tip Subprocess (AFM_utils.py:1136)**
- **Impact**: 0.9s (10-16%)
- **Cause**: Process spawn + communication overhead
- **Code**: `subprocess.run` with `communicate()` - blocking I/O
- **Scaling**: Constant (independent of molecule size)

**E. `np.gradient` on 3D Fields**
- **Impact**: ~0.1-0.2s per field
- **Cause**: Called in [scale_pauli_field](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:872:0-879:25), [compute_es_conv_field](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:892:0-906:22), [compute_dispersion_grid](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:1015:0-1049:29)
- **Code**: Lines 879, 906, 1049 in AFM.py
- **Scaling**: O(nGridPoints) per field

## 2. Debugging/Visualization Guards

### **Where to Add Guards:**

**A. Plotting Gates (Priority 1):**
```python
# In run_afm_from_xyz (line 1225)
# Add parameter: diagnostic_plots=True
# Gate lines 1436, 1474:
if diagnostic_plots:
    plot_diagnostic_panel(...)
    plot_diagnostic_slices(...)
```

**B. Verbosity System (Priority 1):**
```python
# At module level in AFM_utils.py, AFM.py, Grid_dftb.py
import os
_VERBOSITY = int(os.environ.get('AFM_VERBOSITY', '1'))

# Gate all print statements:
if _VERBOSITY >= 1:
    print(f"[INFO] ...")
if _VERBOSITY >= 2:
    print(f"[DEBUG] ...")
```

**C. Intermediate Save Gates (Priority 2):**
```python
# In run_afm_pipeline (line 871)
# Add parameter: save_intermediates=False
# Gate np.save calls (lines 957, 1020, 1027, 1038, 1048, 1050):
if save_intermediates:
    np.save(...)
```

**D. Specific Hotspot Guards:**

**In [compute_dispersion_grid](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:1015:0-1049:29) (AFM.py:1016):**
- Gate print statement (line 1011) with VERBOSITY >= 2
- Add timing if VERBOSITY >= 1

**In [build_tasks](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:283:4-358:38) (Grid_dftb.py:284):**
- Gate block statistics prints (lines 305-327) with VERBOSITY >= 2
- Add timing per block if VERBOSITY >= 3

**In [pp_relax_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:791:0-843:30) (AFM.py:792):**
- Gate per-height print (line 823) with VERBOSITY >= 1
- Gate iteration timing with VERBOSITY >= 3

**In Grid_dftb.py:**
- Lines 191-192, 377, 398, 400, 446-447: Already gated by `verbosity` parameter
- Lines 210, 231, 275: Timing prints - gate with `if verbosity >= 1`

## 3. Easy Performance Wins

### **High Impact (Low Effort):**

**A. Gate Diagnostic Plotting**
- **Effort**: 2 lines (add if condition)
- **Impact**: 3.2s saved (37% on large system)
- **Location**: AFM_utils.py:1436, 1474

**B. Cache CO Tip Grids**
- **Effort**: Check if file exists before subprocess call
- **Impact**: 0.9s saved (10%)
- **Location**: AFM_utils.py:1136-1154
```python
co_tip_file = os.path.join(co_tip_dir, 'co_tip_grid.npz')
if os.path.exists(co_tip_file):
    data = np.load(co_tip_file)
    rho_tip_total = data['rho_total']
    rho_tip_delta = data['rho_delta']
else:
    # compute subprocess
```

**C. Use GPU Task Building**
- **Effort**: Set `use_gpu_tasks=True` in calls
- **Impact**: 0.3-0.5s saved (build_tasks faster on GPU)
- **Location**: Grid_dftb.py already has [build_tasks_gpu](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:158:4-281:38) (line 159)
- **Note**: Currently using CPU [build_tasks](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:283:4-358:38) by default

### **Medium Impact (Medium Effort):**

**D. Combine Gradient with Field Computation**
- **Effort**: Write OpenCL kernel that computes field + gradient in one pass
- **Impact**: 0.1-0.2s per field (avoid separate np.gradient)
- **Location**: AFM.py functions: [scale_pauli_field](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:872:0-879:25), [compute_es_conv_field](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:892:0-906:22), [compute_dispersion_grid](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:1015:0-1049:29)

**E. Precompute Interpolation Weights**
- **Effort**: In [pp_relax_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:791:0-843:30), precompute grid coordinates
- **Impact**: 0.05-0.1s (avoid repeated coordinate calculations)
- **Location**: AFM.py:812-814

## 4. Python → OpenCL Porting Priorities

### **Priority 1 (Critical):**

**A. [compute_dispersion_grid](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:1015:0-1049:29) → OpenCL**
- **Why**: 307x slowdown on large system, pure Python atom loop
- **Complexity**: Low - simple C6/r^6 formula per atom per grid point
- **Implementation**: 
  - Kernel: Each work item = grid point, loop over atoms
  - Or: Each work item = atom, add to all grid points (reduction)
- **Expected Speedup**: 10-50x (GPU parallelism)

**B. [build_tasks](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:283:4-358:38) → Use GPU Version**
- **Why**: 199x slowdown, 181k sphere-AABB checks
- **Complexity**: Very Low - [build_tasks_gpu](cci:1://file:///home/prokop/git/FireCore/pyBall/DFTB/Grid_dftb.py:158:4-281:38) already exists (line 159)
- **Implementation**: Just change default from CPU to GPU
- **Expected Speedup**: 5-20x (GPU already implemented)

### **Priority 2 (High Impact):**

**C. `np.gradient` → OpenCL**
- **Why**: Called on every field (Pauli, ES, vdW), 3D stencil
- **Complexity**: Low - central difference stencil
- **Implementation**: Single kernel computing 3 gradients simultaneously
- **Expected Speedup**: 2-5x (memory bandwidth bound)

**D. Field Gradient Computation in Kernels**
- **Why**: Avoid separate gradient pass
- **Complexity**: Medium - modify existing kernels
- **Implementation**: Compute gradient alongside field in FFT/convolution kernels
- **Expected Speedup**: 1.5-2x (reduce memory passes)

### **Priority 3 (Medium Impact):**

**E. [pp_relax_2d](cci:1://file:///home/prokop/git/FireCore/pyBall/OCL/AFM.py:791:0-843:30) Interpolation → OpenCL**
- **Why**: scipy interpolation called 318 times (H2O) to 1000+ times (large)
- **Complexity**: Medium - trilinear interpolation
- **Implementation**: GPU kernel for force field interpolation at probe positions
- **Expected Speedup**: 5-10x (batch processing)

**F. STM Orbital Projection Batching**
- **Why**: Multiple kernel launches per orbital
- **Complexity**: Medium - batch orbital coefficients
- **Implementation**: Single kernel processing all orbitals
- **Expected Speedup**: 2-5x (reduce kernel launch overhead)

## Summary Table

| Priority | Component | Location | Impact | Effort | Type |
|----------|-----------|----------|--------|--------|------|
| 1 | Diagnostic plotting | AFM_utils.py:1436,1474 | 37% (3.2s) | Very Low | Guard |
| 1 | compute_dispersion_grid | AFM.py:1016 | 7% (0.6s) | Low | OpenCL |
| 1 | build_tasks (use GPU) | Grid_dftb.py:284 | 4% (0.3s) | Very Low | Config |
| 2 | CO tip caching | AFM_utils.py:1136 | 10% (0.9s) | Low | Cache |
| 2 | np.gradient | AFM.py:879,906,1049 | 2-3% | Low | OpenCL |
| 3 | pp_relax_2d interp | AFM.py:812 | 1-2% | Medium | OpenCL |
| 3 | STM batching | AFM_utils.py:720 | Variable | Medium | OpenCL |

**Recommended Action Order:**
1. Add diagnostic plotting guard (immediate, 3.2s saved)
2. Enable GPU task building (immediate, 0.3s saved)
3. Cache CO tip grids (quick, 0.9s saved)
4. Port compute_dispersion_grid to OpenCL (medium effort, 0.6s saved)
5. Add verbosity system (cleanup for profiling)


# Cursort-Composer

## 1) Big bottlenecks / overheads (validated against current code + your baseline tests)

### A. Plotting + disk I/O dominates (especially for small/medium cases)
- **`tests/tAFM/pyocl_fdbm/test_full_pipeline.py`** calls:
  - `afm_utils.plot_diagnostic_panel(...)` and
  - `afm_utils.plot_diagnostic_slices(...)`
  **unconditionally at the end of the script** (no `args.plot_steps` guard).
- Even when you run with `--no_plot_steps`, this file still plots, so matplotlib `savefig()` overhead can easily become the majority of runtime (matches the baseline “~3s plotting”).
- Separately, in **`pyBall/OCL/AFM_utils.py:run_afm_pipeline()`**, a lot of `np.save(...)` for large arrays is done **unconditionally** (overlap/fields/grads/df/tip_disp/stm_grid). Even if plotting is off, this can be costly.

### B. CO tip computation subprocess + large array pad/roll
- **`pyBall/OCL/AFM_utils.py:_call_compute_co_tip_script()`**
  - `subprocess.run(...)` spawn + IPC shows up as a sizable constant overhead in your baseline (~0.9–1.0s).
- **`pyBall/OCL/AFM_utils.py:_pad_and_roll_co_tip()`**
  - Allocates a full `target_shape` 3D array and does multiple `np.roll` operations (memory-bandwidth heavy for large grids).
  - This is “one-time per run” but still significant.

### C. Dense-orbital preparation is CPU-heavy (`prepare_orbital_projection` → `build_tasks`)
- The baseline hotspot **`prepare_orbital_projection`** aligns strongly with:
  - **`pyBall/DFTB/Grid_dftb.py:GridProjector.prepare_orbital_projection()`**
    - It calls **CPU** `self.build_tasks(...)`, not `build_tasks_gpu(...)`.
  - CPU `build_tasks()` contains nested loops over grid blocks and atoms and calls:
    - **`check_overlap_sphere_aabb()`**
  repeatedly, i.e. exactly the pattern your baseline reports as “`check_overlap_sphere_aabb` called 181k+ times”.

This is a “huge win” because an OpenCL/GPU version (`build_tasks_gpu`) already exists in the same file.

### D. `compute_dispersion_grid` scales badly: atom-loop × full 3D grid
- **`pyBall/OCL/AFM.py:compute_dispersion_grid()`**
  - Builds `XX,YY,ZZ` via `np.meshgrid(...)` and then does:
    - `for ia in range(len(atomPos)):` and updates a full `E_vdw` grid each time.
  - Complexity is effectively **O(nAtoms * nGridPoints)** with large temporary arrays.
  - Your baseline showed it as the #1/new hotspot for large systems (307x slower vs H2O), which is exactly what this structure implies.

### E. Relaxation + interpolation overhead from SciPy `map_coordinates`
Even though your baseline suggests it’s not the *top* hotspot after the above are addressed, it is still a real overhead:
- **`pyBall/OCL/AFM_utils.py:compose_and_relax()`**
  - Builds `F_total = -(grads_pauli + grads_es + grads_vdw)` (large 4D allocation).
  - `force_func()` interpolates 3 components using SciPy `map_coordinates` (multiple calls).
- **`pyBall/OCL/AFM.py:pp_relax_2d()`**
  - Calls `force_func(pp_flat)` once per height (big map work),
  - then runs `N_RELAX` iterations where each iteration calls `_interp2d()` (SciPy interpolation again).

So: it’s “important” for correctness and sometimes can become dominant at larger scan sizes, but in your baseline it’s lower than plotting + CPU projection preparation + vdW.

---

## 2) Where debugging/visualization guards should go (to make profiling meaningful)

### Priority 1: Make *all* plotting conditional on a single flag
1. **Fix test-level plotting leak**
   - In **`tests/tAFM/pyocl_fdbm/test_full_pipeline.py`**, wrap:
     - `plot_diagnostic_panel(...)`
     - `plot_diagnostic_slices(...)`
     - STM plotting (`plot_stm`) and tip displacement plotting
     behind `if args.plot_steps:`.

2. **Ensure diagnostic plotting inside library code is also gated consistently**
   - `pyBall/OCL/AFM_utils.py:run_afm_from_xyz(...)` already gates many things with `plot_steps`, but the test currently doesn’t rely on that.

### Priority 2: Gate intermediate `np.save(...)` independently from plotting
- In **`pyBall/OCL/AFM_utils.py:run_afm_pipeline()`**, add (conceptually) a separate boolean like `save_intermediates` and gate `np.save(...)` for:
  - `V_ES.npy`, `overlap_raw.npy`
  - `E_Pauli_field.npy`, `grads_E_Pauli.npy`
  - `E_ES_field.npy`, `grads_E_ES.npy`
  - `E_vdw_field.npy`, `grads_E_vdw.npy`
  - `df.npy`, `tip_disp_dx.npy`, `tip_disp_dy.npy`, `stm_grid.npy`

This keeps profiling focused on compute rather than filesystem writes.

### Priority 3: Verbosity guards for inner loops
- Wherever you have frequent prints (notably relaxation per-height and per-iteration), gate them by a verbosity level so they cannot distort timing:
  - **`pyBall/OCL/AFM.py:pp_relax_2d()`** (per-height progress, and potentially iteration-level details if added later)
  - **`pyBall/OCL/AFM.py:compute_dispersion_grid()`** doesn’t print in the inner atom loop now, but gradient/FFT routines may print elsewhere—ensure those are verbosity-gated.

---

## 3) Significant performance improvements that look *easy* (highest ROI, lowest complexity)

### A. Stop doing plotting when `--no_plot_steps`
- The biggest “easy win” is the test-level unconditional plotting (matplotlib can dominate runtime entirely).
- This should make profiling results reflect compute reality.

### B. Gate intermediate saving (`np.save`) separately from plotting
- Prevent large-array writes during profiling runs.

### C. Use GPU task building for dense-orbital projection
- In **`pyBall/DFTB/Grid_dftb.py:prepare_orbital_projection()`**, switch from CPU `build_tasks()` to GPU `build_tasks_gpu()`.
- Rationale: your baseline shows `prepare_orbital_projection` and `check_overlap_sphere_aabb` explode in large systems, and a GPU implementation already exists right next to it.

### D. Cache CO tip outputs + avoid subprocess in repeated runs
- If `co_tip_dir` is provided and contains `co_rho_total.npy` / `co_rho_delta.npy`, that already bypasses subprocess.
- Otherwise, you can at least make “compute once per step/margin” behavior explicit and persistent on disk (your baseline shows subprocess overhead is ~constant).

### E. Optimize vdW grid construction to avoid `meshgrid` temp explosion
- In **`compute_dispersion_grid`**, the easiest micro-optimization is to avoid materializing full `XX,YY,ZZ` if possible (use coordinate vectors + broadcasting, or compute `r2` without allocating three full 3D arrays).
- Larger structural win is OpenCL (next section), but even NumPy-side changes can be meaningful.

---

## 4) What should be ported from Python → OpenCL (for heavy numerics)

I’d prioritize by **(a)** how much time they take and **(b)** how well they fit GPU parallelism.

### Priority 1 (strong candidates)
1. **`pyBall/OCL/AFM.py:compute_dispersion_grid()` → OpenCL**
   - Per-grid-point work; formula is simple (`C6/r^6` with damping/RA).
   - Kernel strategy: each work-item computes one grid point and loops over atoms (classic).
   - Also compute gradients in-kernel or as a second kernel.

2. **Interpolation inside relaxation (`compose_and_relax` / `pp_relax_2d`) → OpenCL**
   - SciPy `map_coordinates` is hard to beat with CPU in inner loops.
   - GPU version should implement trilinear (3D) and bilinear (2D) sampling over regular grids.
   - Likely implement:
     - “sample forces at probe positions” in one batched kernel,
     - and “sample 2D FF map at probe positions” in the relaxation iterations.

### Priority 2 (good, but more invasive)
3. **Replace `np.gradient` stencils for Pauli/ES/vdW fields**
   - Right now gradients are computed via `np.gradient` after FFT/convolution.
   - OpenCL gradient kernels can reduce memory traffic by:
     - computing energy field and gradients in one pass (best),
     - or at least keep everything device-side (still helpful).

4. **STM orbital projection batching**
   - STM already uses OpenCL projection kernels, but **`compute_stm` / `compute_bond_resolved_stm`** still launches projection repeatedly per height and per orbital.
   - A more GPU-friendly approach would batch multiple orbitals (or multiple heights) into fewer kernel launches.

### Priority 3 (lower ROI)
5. **CO tip pad+roll**
   - It’s memory heavy but not usually the top compute bottleneck once plotting is fixed.
   - Still could be moved to GPU if you’re already restructuring GPU workflows.


---

# Visual Code Flow with Bottlenecks & Debug Output

## Legend
- 🔴 **Critical bottleneck** (>5% total time on large system)
- 🟠 **Medium bottleneck** (1-5% total time)
- 🔵 **Print statements** (verbosity needed)
- 🟣 **Plotting/figures** (needs gating)
- 🟢 **File saves (.npy)** (needs gating)

## AFM Pipeline Flow (run_afm_from_xyz)

```mermaid
flowchart TD
    Start([run_afm_from_xyz<br/>AFM_utils.py:1225]) --> LoadXYZ[Load XYZ file<br/>au.load_xyz]

    LoadXYZ --> GetDensity{Dense projection?}
    GetDensity -->|Yes| GetDense[get_density_from_dftb_dense<br/>AFM_utils.py]
    GetDensity -->|No| GetSparse[get_density_from_dftb_plus<br/>AFM_utils.py]

    subgraph DFTB_Density ["DFTB Density Calculation"]
        GetDense --> DFTBRun[DFTB+ SCF]
        GetSparse --> DFTBRun
        DFTBRun --> ProjTotal[project_dftb_density<br/>Grid_dftb.py]
        ProjTotal --> ProjSetup[prepare_orbital_projection<br/>Grid_dftb.py:1803]
        ProjSetup --> BuildTasks[build_tasks<br/>Grid_dftb.py:284]
        BuildTasks --> BuildTasksBottleneck🔴[🔴 check_overlap_sphere_aabb<br/>0.3s, 181k calls]
        BuildTasksBottleneck🔴 --> ProjNeutral[project_neutral_density<br/>Grid_dftb.py]
        ProjNeutral --> ProjSave🟢[save: rho_total.npy<br/>rho_neutral.npy]
    end

    ProjSave🟢 --> GridSpecSave🟢[save: grid_spec.txt]

    GridSpecSave🟢 --> CheckCO{CO tip cached?}
    CheckCO -->|No| CoGrid[_compute_co_tip_grid<br/>AFM_utils.py:1119]
    CoGrid --> Subprocess[_call_compute_co_tip_script<br/>AFM_utils.py:1136]
    Subprocess --> SubprocessBottleneck🟠[🟠 subprocess.spawn<br/>0.9s]
    SubprocessBottleneck🟠 --> PadRoll[_pad_and_roll_co_tip<br/>AFM_utils.py:1157]
    PadRoll --> COSave🟢[save: co_rho_*.npy]
    CheckCO -->|Yes| LoadCO[load cached CO tip]

    LoadCO --> RunPipeline[run_afm_pipeline<br/>AFM_utils.py:871]

    subgraph AFM_Pipeline ["run_afm_pipeline - Core AFM Computation"]
        RunPipeline --> ES_Check{V_ES provided?}
        ES_Check -->|No| FFT_Poisson[afm.fft_poisson<br/>AFM.py]
        ES_Check -->|Yes| ES_Provided[use provided V_ES]
        FFT_Poisson --> ES_Save🟢[save: V_ES.npy]
        ES_Provided --> ES_Save🟢

        ES_Save🟢 --> PauliOverlap[afm.compute_pauli_overlap<br/>AFM.py:858]
        PauliOverlap --> PauliFFT[np.fft.fftn/ifftn]
        PauliFFT --> PauliSave�[save: overlap_raw.npy]

        PauliSave� --> PauliScale[afm.scale_pauli_field<br/>AFM.py:872]
        PauliScale --> PauliGradBottleneck🟠[🟠 np.gradient<br/>0.1-0.2s]
        PauliGradBottleneck🟠 --> PauliSave2🟢[save: E_Pauli_field.npy<br/>grads_E_Pauli.npy]

        PauliSave2🟢 --> ESConv[afm.compute_es_conv_field<br/>AFM.py:892]
        ESConv --> ESFFT[np.fft.fftn/ifftn]
        ESFFT --> ESGradBottleneck🟠[🟠 np.gradient<br/>0.1-0.2s]
        ESGradBottleneck🟠 --> ESSave🟢[save: E_ES_field.npy<br/>grads_E_ES.npy]

        ESSave🟢 --> DispGrid[afm.compute_dispersion_grid<br/>AFM.py:1016]
        DispGrid --> DispBottleneck🔴[🔴 Python atom loop<br/>0.6s, 42M ops]
        DispBottleneck🔴 --> DispGradBottleneck🟠[🟠 np.gradient<br/>0.1-0.2s]
        DispGradBottleneck🟠 --> DispSave🟢[save: E_vdw_field.npy<br/>grads_E_vdw.npy]

        DispSave🟢 --> Compose[compose_and_relax<br/>AFM_utils.py:467]
        Compose --> ForceFunc[define force_func<br/>uses map_coordinates]
        ForceFunc --> PP_Relax[afm.pp_relax_2d<br/>AFM.py:792]
        PP_Relax --> InterpBottleneck🟠[🟠 scipy map_coordinates<br/>0.05-0.1s]
        InterpBottleneck🟠 --> ComputeDF[afm.compute_df<br/>AFM.py]
        ComputeDF --> RelaxSave🟢[save: df.npy<br/>tip_disp_*.npy]

        RelaxSave🟢 --> STM_Check{STM requested?}
        STM_Check -->|Yes| ComputeSTM[compute_stm<br/>AFM_utils.py:684]
        ComputeSTM --> STM_Proj[projector.project_orbital_dense_points<br/>Grid_dftb.py]
        STM_Proj --> STM_Save🟢[save: stm_grid.npy]
        STM_Check -->|No| NoSTM[skip STM]
    end

    NoSTM --> DiagPanel[plot_diagnostic_panel<br/>AFM_utils.py:1436]
    STM_Save🟢 --> DiagPanel
    DiagPanel --> DiagPanelBottleneck🟠[🟠 matplotlib<br/>2.2s]
    DiagPanelBottleneck🟠 --> DiagSlices[plot_diagnostic_slices<br/>AFM_utils.py:1474]
    DiagSlices --> DiagSlicesBottleneck🟠[🟠 matplotlib<br/>2.1s]

    DiagSlicesBottleneck🟠 --> End([Done])

    style Start fill:#e1f5e1
    style End fill:#e1f5e1
    style BuildTasksBottleneck🔴 fill:#ffcdd2
    style SubprocessBottleneck🟠 fill:#ffe0b2
    style PauliGradBottleneck🟠 fill:#ffe0b2
    style ESGradBottleneck🟠 fill:#ffe0b2
    style DispBottleneck🔴 fill:#ffcdd2
    style DispGradBottleneck🟠 fill:#ffe0b2
    style InterpBottleneck🟠 fill:#ffe0b2
    style DiagPanelBottleneck🟠 fill:#ffe0b2
    style DiagSlicesBottleneck🟠 fill:#ffe0b2
    style ProjSave� fill:#c8e6c9
    style GridSpecSave� fill:#c8e6c9
    style COSave🟢 fill:#c8e6c9
    style ES_Save🟢 fill:#c8e6c9
    style PauliSave🟢 fill:#c8e6c9
    style PauliSave2🟢 fill:#c8e6c9
    style ESSave🟢 fill:#c8e6c9
    style DispSave🟢 fill:#c8e6c9
    style RelaxSave🟢 fill:#c8e6c9
    style STM_Save� fill:#c8e6c9
    style DiagPanel fill:#e1bee7
    style DiagSlices fill:#e1bee7
```

## Bottleneck Impact Summary (Large System: TBTAP, 60 atoms)

| Bottleneck | Time | % of Total | Location | Type |
|------------|------|------------|----------|------|
| 🔴 compute_dispersion_grid | 0.615s | 7% | AFM.py:1045-1047 | Python atom loop |
| 🔴 build_tasks (check_overlap) | 0.321s | 4% | Grid_dftb.py:317-319 | Python nested loops |
| 🟠 Diagnostic plotting (panel) | 2.162s | 25% | AFM_utils.py:1436 | Matplotlib |
| 🟠 Diagnostic plotting (slices) | 2.136s | 25% | AFM_utils.py:1474 | Matplotlib |
| 🟠 CO tip subprocess | 0.906s | 10% | AFM_utils.py:1136 | Process spawn |
| 🟠 np.gradient (×3 fields) | ~0.4s | 5% | AFM.py:879,906,1049 | NumPy |
| 🟠 scipy interpolation | ~0.1s | 1% | AFM.py:812 | SciPy |

**Total avoidable (plotting + subprocess + saves): ~4.5s (52%)**
**Total computation bottlenecks: ~1.4s (16%)**

## Debug Output Locations (Needs Verbosity Guards)

### Print Statements (🔵)
- `Grid_dftb.py:191-192, 377, 398, 400, 446-447` - Already gated by `verbosity`
- `Grid_dftb.py:210, 231, 275` - Timing prints (needs gating)
- `AFM.py:1011` - rho_grid shape (needs VERBOSITY >= 2)
- `AFM.py:823` - per-height progress (needs VERBOSITY >= 1)
- `AFM_utils.py:1136` - CO tip subprocess info (needs gating)

### Plotting/Figures (🟣)
- `AFM_utils.py:1220` - CO tip diagnostic plot (needs `plot_co_tip` flag)
- `AFM_utils.py:1436` - plot_diagnostic_panel (needs `diagnostic_plots` flag) **CRITICAL**
- `AFM_utils.py:1474` - plot_diagnostic_slices (needs `diagnostic_plots` flag) **CRITICAL**
- `AFM_utils.py:720-722` - STM plots (needs `plot_steps` flag)

### File Saves (🟢)
- `Grid_dftb.py:377` - grid_spec.txt (useful, keep)
- `AFM_utils.py:1154` - CO tip grids (already cached, good)
- `AFM.py:957, 1020, 1027, 1038, 1048, 1050` - Field/gradient saves (needs `save_intermediates` flag)
- `AFM_utils.py:1412` - df.npy, tip_disp.npy (needs `save_intermediates` flag)

## Optimization Impact Visualization

```
Before Optimization (Large System):
┌─────────────────────────────────────────────────┐
│ Total: 8.699s                                   │
├─────────────────────────────────────────────────┤
│ 🟣 Plotting:       3.2s  (37%)  ← GATE THIS      │
│ 🟠 Subprocess:     0.9s  (10%)  ← CACHE THIS     │
│ 🔴 Dispersion:     0.6s  (7%)   ← PORT TO GPU   │
│ 🔴 Task building:   0.3s  (4%)   ← USE GPU       │
│ 🟠 Gradients:      0.4s  (5%)   ← PORT TO GPU   │
│ 🟢 File saves:     0.3s  (3%)   ← GATE THIS     │
│ 🔵 Prints:         0.2s  (2%)   ← VERBOSITY     │
│ Actual compute:    2.8s  (32%)  ← KEEP          │
└─────────────────────────────────────────────────┘

After Quick Wins (Gates + Cache + GPU tasks):
┌─────────────────────────────────────────────────┐
│ Total: ~3.5s (60% faster)                       │
├─────────────────────────────────────────────────┤
│ Actual compute:    2.8s  (80%)                  │
│ 🟠 Gradients:      0.4s  (11%)  ← PORT TO GPU   │
│ 🔴 Dispersion:     0.3s  (9%)   ← PORT TO GPU   │
└─────────────────────────────────────────────────┘

After Full Optimization (All GPU):
┌─────────────────────────────────────────────────┐
│ Total: ~2.0s (77% faster)                       │
├─────────────────────────────────────────────────┤
│ Actual compute:    2.0s  (100%)                 │
└─────────────────────────────────────────────────┘
```

---

# Parity Checkpoint System for Optimization Validation

## Overview

A reusable parity checkpoint system has been implemented in `pyBall/DFTB/TestUtils.py` to save reference data at key pipeline checkpoints and verify parity after optimizations.

**Available functions:**
- `save_checkpoint()` - Save data and configuration to disk
- `load_checkpoint()` - Load checkpoint and verify configuration hash
- `compare_checkpoint()` - Compare new data against reference
- `CheckpointManager` - High-level manager for checkpoint operations

## Checkpoint Locations in AFM Pipeline

### 1. DFTB Density Calculation (AFM_utils.py:1319-1325)

**Location:** After `get_density_from_dftb_plus()` or `get_density_from_dftb_dense()` returns

**Checkpoint name:** `dftb_density`

**Data to save:**
```python
{
    'rho_scf': d['rho_scf'],
    'rho_na': d['rho_na'],
    'rho_diff': d['rho_diff'],
    'V_ES': d['V_ES'],
    'origin': d['origin'],
    'step': step,
    'ngrid': d['ngrid']
}
```

**Configuration:**
```python
{
    'xyz_file': xyz_file,
    'step': step,
    'margin': margin,
    'z_extra': z_extra,
    'basis': basis,
    'slako_prefix': slako_prefix,
    'use_dense_projection': use_dense_projection,
    'max_shells': max_shells
}
```

**Rationale:** Foundation - all subsequent computations depend on density. Critical for validating GPU task building changes.

---

### 2. Pauli Overlap (AFM_utils.py:956-958)

**Location:** After `afm.compute_pauli_overlap()` in `run_afm_pipeline()`

**Checkpoint name:** `pauli_overlap`

**Data to save:**
```python
{
    'overlap_raw': overlap_raw,
    'origin': origin,
    'step': step
}
```

**Configuration:**
```python
{
    'step': step,
    'tip_rolled': True
}
```

**Rationale:** Tests FFT convolution. Critical for validating any FFT operation changes.

---

### 3. Pauli Field Scaling (AFM_utils.py:1010-1020)

**Location:** After `afm.scale_pauli_field()` in `run_afm_pipeline()`

**Checkpoint name:** `pauli_field`

**Data to save:**
```python
{
    'E_pauli_field': E_pauli_field,
    'grads_pauli': grads_pauli,
    'A_pauli': A_pauli,
    'beta_pauli': beta_pauli
}
```

**Configuration:**
```python
{
    'A_pauli': A_pauli,
    'beta_pauli': beta_pauli,
    'step': step
}
```

**Rationale:** Tests power-law scaling and gradient computation. Validates `np.gradient` replacements.

---

### 4. Electrostatic Convolution (AFM_utils.py:1024-1028)

**Location:** After `afm.compute_es_conv_field()` in `run_afm_pipeline()`

**Checkpoint name:** `es_field`

**Data to save:**
```python
{
    'E_ES_field': E_ES_field,
    'grads_ES': grads_ES,
    'V_ES': V_ES
}
```

**Configuration:**
```python
{
    'step': step,
    'tip_rolled': True
}
```

**Rationale:** Tests ES convolution and gradient computation.

---

### 5. Dispersion Grid (AFM_utils.py:1032-1039)

**Location:** After `afm.compute_dispersion_grid()` in `run_afm_pipeline()`

**Checkpoint name:** `dispersion_field`

**Data to save:**
```python
{
    'E_vdw': E_vdw,
    'grads_vdw': grads_vdw,
    'atomPos': atomPos,
    'atomTypes': atomTypes,
    'origin': origin,
    'step': step,
    'ngrid': ngrid
}
```

**Configuration:**
```python
{
    'C6_CO': vdw_params['C6_CO'],
    'RA': 1.5,
    'step': step,
    'margin': margin
}
```

**Rationale:** **Critical** - this is the #1 bottleneck to be ported to OpenCL. Must validate exact parity before/after GPU porting.

---

### 6. Force Field Composition (AFM_utils.py:493-494)

**Location:** In `compose_and_relax()`, after `F_total = -(grads_pauli + grads_es + grads_vdw)`

**Checkpoint name:** `force_field`

**Data to save:**
```python
{
    'F_total': F_total,
    'grads_pauli': grads_pauli,
    'grads_es': grads_es,
    'grads_vdw': grads_vdw
}
```

**Configuration:**
```python
{
    'origin': origin,
    'step': step
}
```

**Rationale:** Tests force field composition. Validates gradient computations.

---

### 7. Probe Relaxation (AFM_utils.py:507-509)

**Location:** After `afm.pp_relax_2d()` returns in `compose_and_relax()`

**Checkpoint name:** `relaxation`

**Data to save:**
```python
{
    'df': df,
    'tip_disp': tip_disp,
    'FEs_relax': FEs_relax
}
```

**Configuration:**
```python
{
    'K_LAT': K_LAT,
    'N_RELAX': 50,
    'PP_MARGIN': 2.0,
    'step': step,
    'scan_xs': scan_xs.tolist(),
    'scan_ys': scan_ys.tolist(),
    'heights': heights.tolist()
}
```

**Rationale:** Tests full relaxation loop including scipy interpolation. Validates interpolation optimizations.

---

### 8. STM Grid (AFM_utils.py:1084-1086)

**Location:** After `compute_stm()` or `compute_bond_resolved_stm()` in `run_afm_pipeline()`

**Checkpoint name:** `stm_grid`

**Data to save:**
```python
{
    'stm_grid': stm_grid,
    'tip_disp': tip_disp
}
```

**Configuration:**
```python
{
    'lumo_offsets': lumo_offsets,
    'use_exp_basis': use_exp_basis,
    'exp_beta': exp_beta,
    'exp_r0': exp_r0,
    'bond_resolved': bond_resolved
}
```

**Rationale:** Tests STM orbital projection. Validates projection batching optimizations.

---

## Integration Example

```python
from pyBall.DFTB.TestUtils import CheckpointManager

# Initialize checkpoint manager
checkpoint_mgr = CheckpointManager(output_dir, 'afm_pipeline_ref')

# Base configuration for all checkpoints
base_config = {
    'step': step,
    'origin': origin.tolist(),
    'ngrid': ngrid.tolist(),
    'pauli_params': pauli_params,
    'vdw_params': vdw_params
}

# Save checkpoint (baseline run)
if not checkpoint_mgr.exists('dispersion_field'):
    checkpoint_mgr.save('dispersion_field', {
        'E_vdw': E_vdw,
        'grads_vdw': grads_vdw
    }, config=base_config)

# Verify checkpoint (after optimization)
if verify_mode:
    result = checkpoint_mgr.compare('dispersion_field', {
        'E_vdw': E_vdw,
        'grads_vdw': grads_vdw
    }, config=base_config, tolerances={'E_vdw': 1e-5, 'grads_vdw': 1e-5})
    assert result['match'], f"Dispersion parity failed: {result['details']}"
```

## Verification Workflow

**Step 1: Baseline (Before Optimization)**
```bash
python test_full_pipeline.py h2o.xyz --save-checkpoints --checkpoint-dir ./checkpoints/h2o_baseline
```

**Step 2: Optimize**
Make code changes (e.g., port `compute_dispersion_grid` to OpenCL).

**Step 3: Verify Parity**
```bash
python test_full_pipeline.py h2o.xyz --verify-checkpoints --checkpoint-dir ./checkpoints/h2o_baseline
```

**Step 4: Debug**
If parity fails, examine RMS/max error values and fix the bug.

## Tolerance Guidelines

**Important:** These tolerances account for:
- **CPU vs GPU precision differences** (float64 vs float32)
- **Operation ordering differences** (different algorithms can produce slightly different results)
- **Numerical accumulation errors** (different summation orders)

The hash is only used for configuration parameters (grid size, physics parameters, etc.) to detect when the computation setup changes. Data comparison uses numerical tolerances, not bitwise comparison.

| Data Type | CPU→CPU Tolerance | CPU↔GPU Tolerance | Rationale |
|-----------|------------------|-------------------|-----------|
| Densities (e/Å³) | 1e-6 | 1e-4 | Float32 precision ~7 decimal digits |
| Energies (eV) | 1e-5 | 1e-3 | Accumulation errors, float32 |
| Forces (eV/Å) | 1e-5 | 1e-3 | Gradient precision, float32 |
| Positions (Å) | 1e-8 | 1e-6 | Should be exact for same grid |
| Gradients | 1e-6 | 1e-4 | Numerical differentiation, float32 |
| STM signals | 1e-7 | 1e-5 | Squared amplitudes, more sensitive |

**Note:** When porting from CPU (float64) to GPU (float32), use the CPU↔GPU tolerances. For CPU-only optimizations (e.g., algorithm changes), CPU→CPU tolerances are sufficient.

---

# USER

OK, now plese implement the debug guards, make sure they can be easily controled globaly, there should be some debug_print level some dedicated debug_print function in module globals.py which is imported from everywhere, and maybe separately some debug_save and debug_plot level to control how many figures (.png) and .npy or ther axuliary data fliles we plot during the calculation which are not strictly necessary

This should be first level of optimization, safe without side effects. Then you should run it

before you start anything we should save some data-files which we can use for reference. I dded system for that in @pyBall/DFTB/TestUtils.py 
and there is documentation of it in @doc/Topics/AFM/AFM_FDBM_profiling_optimization.md , read it before you start and use it

---

# Cursor-Composer

Implemented global debug guards + auxiliary output throttling, using a new shared module:

- Added `pyBall/globals.py` with:
  - `debug_print(level, msg)`
  - `debug_save_enabled(level)` (controls `.npy` auxiliary saves)
  - `debug_plot_enabled(level)` (controls extra `.png` diagnostics)
  - Defaults are set to preserve current behavior unless you explicitly lower env vars:
    - `AFM_DEBUG_PRINT_LEVEL` (default `1`)
    - `AFM_DEBUG_SAVE_LEVEL` (default `2`)
    - `AFM_DEBUG_PLOT_LEVEL` (default `2`)

### What changed

1) **AFM pipeline auxiliary `.npy` saving is now guarded**
- File: `pyBall/OCL/AFM_utils.py`
- In `run_afm_pipeline()`, all intermediate `np.save(...)` calls (V_ES / overlap / fields / gradients / df / tip_disp / stm_grid) are now wrapped with `debug_save_enabled(2)`.

2) **High-frequency console output is now `debug_print`-guarded**
- File: `pyBall/OCL/AFM.py`
- In `pp_relax_2d()`, the per-height `iz=... building ...` and final `Fz_relax...` prints now use `debug_print(1, ...)`.

3) **Extra plotting is now controlled by global debug_plot**
- File: `tests/tAFM/pyocl_fdbm/test_full_pipeline.py`
- `plot_diagnostic_panel(...)` / `plot_diagnostic_slices(...)` are only called if:
  - `args.plot_steps` is true AND `debug_plot_enabled(2)` is true
- `plot_tip_displacement(...)` and `plot_stm(...)` are also gated (to avoid side-effect plots during profiling).

4) **Checkpoint CLI support (using `CheckpointManager`)**
- File: `tests/tAFM/pyocl_fdbm/test_full_pipeline.py`
- Added arguments:
  - `--save-checkpoints`
  - `--verify-checkpoints`
  - `--checkpoint-dir`
- When saving/verifying, it loads the arrays from the run output directory (e.g. `overlap_raw.npy`, `E_Pauli_field.npy`, `df.npy`, etc.) and stores/compares them via `pyBall/DFTB/TestUtils.py::CheckpointManager`.

### Sanity checks
- I ran `py_compile` on all modified files and also `ReadLints`; no Python syntax/lint issues.

### Running it
I **cannot run the full AFM pipeline in this sandbox** because OpenCL kernels fail to compile on the available CPU OpenCL device (`pocl`). On your GTX3090 environment it should work.

On your machine, use for example:
```bash
cd /home/prokop/git/FireCore/tests/tAFM/pyocl_fdbm
AFM_DEBUG_SAVE_LEVEL=2 AFM_DEBUG_PLOT_LEVEL=2 AFM_DEBUG_PRINT_LEVEL=1 \
python test_full_pipeline.py /home/prokop/git/FireCore/cpp/common_resources/xyz/TBTAP.xyz \
  --basis 3ob-3-1 --use_dense_projection --max_shells 3 --compute_stm --stm_bond_resolved \
  --stm_exp_beta 1.0 --stm_lumo_offsets 1
```

If you want *minimal side effects* for profiling, try:
```bash
AFM_DEBUG_SAVE_LEVEL=0 AFM_DEBUG_PLOT_LEVEL=0 AFM_DEBUG_PRINT_LEVEL=0 ...
```

If you want, tell me what exact checkpoint keys you want to match from the doc (the current implementation saves a practical subset based on what the pipeline already writes to disk).

---

# USER





Based on the profiling documentation, here are the optimization phases from easiest/best ROI to more elaborate:

## **Phase 1: Debug Guards & Easy Wins (DONE + Quick Fixes)**

**Status:** Debug guards implemented by Cursor-Composer ✅

**Remaining quick fixes:**
1. **Gate diagnostic plotting in test script** (AFM_utils.py:1436,1474)
   - Effort: 2 lines
   - Impact: 3.2s saved (37% on large system)
   - Currently called unconditionally despite `--no_plot_steps`

2. **Cache CO tip grids** (AFM_utils.py:1136-1154)
   - Effort: Check file existence before subprocess
   - Impact: 0.9s saved (10%)
   - Cache per step/margin combination

3. **Gate intermediate np.save calls** (AFM_utils.py:run_afm_pipeline)
   - Effort: Add `save_intermediates` flag
   - Impact: Reduces I/O overhead for profiling
   - Files: V_ES, overlap, fields, gradients, df, tip_disp, stm_grid

**Expected total gain:** ~4.1s (47% on large system)

---

## **Phase 2: Enable GPU Task Building (Very Low Effort)**

**Target:** `build_tasks` in Grid_dftb.py:284

**Change:** Use `build_tasks_gpu` instead of CPU `build_tasks`
- Effort: Change default parameter or add flag
- Impact: 0.3-0.5s saved (4% on large system, 199x bottleneck addressed)
- Speedup: 5-20x (GPU version already exists)
- Why: CPU version does 181k sphere-AABB checks, GPU version parallelizes this

**Expected gain:** ~0.4s (5% on large system)

---

## **Phase 3: Port compute_dispersion_grid to OpenCL (Low Effort, High ROI)**

**Target:** AFM.py:1016-1050

**Current bottleneck:** Python loop over atoms × 3D grid points
- 60 atoms × 703k points = 42M operations
- 0.615s on large system (7%, 307x slower than H2O)

**Implementation:**
- Kernel: Each work item = grid point, loop over atoms
- Or: Each work item = atom, add to all grid points (reduction)
- Formula: C6/r^6 per atom per grid point

**Expected gain:** 0.6s → 0.01-0.06s (10-50x speedup)

---

## **Phase 4: Gradient Computation Optimization (Low-Medium Effort)**

**Target:** np.gradient calls in AFM.py:879,906,1049

**Options:**
1. **Port np.gradient to OpenCL** (simpler)
   - Central difference stencil kernel
   - Speedup: 2-5x (memory bandwidth bound)

2. **Combine with field computation** (better)
   - Modify FFT/convolution kernels to compute gradient alongside field
   - Speedup: 1.5-2x (reduce memory passes)

**Expected gain:** 0.1-0.2s per field × 3 fields = 0.3-0.6s

---

## **Phase 5: Interpolation Optimization (Medium Effort)**

**Target:** scipy interpolation in AFM.py:812 (pp_relax_2d)

**Current:** 318 calls (H2O) to 1000+ calls (large) in relaxation loop

**Options:**
1. **Precompute interpolation weights** (quick fix)
   - Precompute grid coordinates
   - Speedup: 0.05-0.1s

2. **Port to OpenCL** (better)
   - GPU kernel for trilinear interpolation
   - Batch process all probe positions
   - Speedup: 5-10x

**Expected gain:** 0.05-0.1s (weights) or 0.2-0.5s (OpenCL)

---

## **Phase 6: STM Orbital Batching (Medium Effort)**

**Target:** AFM_utils.py:720-732 (compute_stm)

**Current:** Nested loops over heights × LUMO orbitals, separate GPU kernel per orbital

**Implementation:**
- Batch all LUMO orbitals into single GPU kernel
- Process all heights in single batch
- Reduce kernel launch overhead

**Expected gain:** 2-5x speedup (variable depending on #orbitals)

---

## **Phase 7: Inline CO Tip Computation (Medium-High Effort)**

**Target:** AFM_utils.py:1136-1154 (_call_compute_co_tip_script)

**Current:** Subprocess spawn + IPC overhead (0.9s constant)

**Implementation:**
- Inline CO tip computation directly in Python
- Avoid subprocess overhead
- Still need DFTB SCF but can reuse density projection

**Expected gain:** 0.9s saved (eliminates subprocess overhead)

---

## **Phase 8: Advanced Optimizations (High Effort)**

**Targets:**
1. **Full relaxation on GPU**
   - Port entire pp_relax_2d to OpenCL
   - Eliminate all Python/SciPy interpolation
   - Speedup: 10-50x

2. **Force field composition optimization**
   - Avoid materializing F_total (4D array)
   - Compute components on-demand
   - Reduce memory bandwidth

3. **FFT batching**
   - Batch multiple FFT operations if possible
   - Reduce kernel launch overhead

---

## **Summary Table**

| Phase | Component | Effort | Impact | Cumulative Gain |
|-------|-----------|--------|--------|-----------------|
| 1 | Plotting guards | Very Low | 3.2s (37%) | 3.2s |
| 1 | CO tip cache | Low | 0.9s (10%) | 4.1s |
| 1 | Gate np.save | Low | 0.1-0.2s | 4.2-4.3s |
| 2 | GPU task building | Very Low | 0.4s (5%) | 4.6-4.7s |
| 3 | dispersion_grid OpenCL | Low | 0.6s (7%) | 5.2-5.3s |
| 4 | Gradient OpenCL | Low-Med | 0.3-0.6s | 5.5-5.9s |
| 5 | Interpolation | Medium | 0.2-0.5s | 5.7-6.4s |
| 6 | STM batching | Medium | Variable | 5.7-6.4s+ |
| 7 | Inline CO tip | Med-High | 0.9s (10%) | 6.6-7.3s+ |
| 8 | Advanced | High | Variable | 6.6-7.3s+ |

**Baseline (large system):** 8.7s
**After Phase 1-3:** ~4.0s (54% reduction)
**After Phase 1-5:** ~3.0s (65% reduction)
**After all phases:** ~1.5-2.0s (77-82% reduction)

---

## **Recommended Implementation Order**

1. **Phase 1** (fix plotting leak, cache CO tip, gate saves) - immediate wins
2. **Phase 2** (enable GPU tasks) - trivial config change
3. **Phase 3** (dispersion_grid OpenCL) - biggest computation win
4. **Phase 4** (gradient optimization) - cleans up remaining Python bottlenecks
5. **Phase 5-8** (advanced) - for further optimization if needed

This order maximizes ROI while keeping each phase manageable and testable with the checkpoint system.

