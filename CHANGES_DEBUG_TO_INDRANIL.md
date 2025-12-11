# Comprehensive Changes: debug/grid_geration_and_test → Indranil Branch

## Document Information
- **Generated**: December 11, 2025
- **Source Branch**: `debug/grid_geration_and_test` (commit `b68be3ac`)
- **Target Branch**: `Indranil` (commit `b88cdd47`)
- **Purpose**: Document what was REMOVED or MODIFIED from debug branch to create clean Indranil branch

---

## Table of Contents
1. [Summary of Differences](#1-summary-of-differences)
2. [Debug Prints Removed](#2-debug-prints-removed)
3. [Profiler Files Excluded](#3-profiler-files-excluded)
4. [Plot/Data Files Excluded](#4-plotdata-files-excluded)
5. [File-by-File Detailed Changes](#5-file-by-file-detailed-changes)

---

## 1. Summary of Differences

### Philosophy
The Indranil branch contains **all functional code** from the debug branch, but with:
- Debug print statements removed
- Profiler files excluded
- Data-specific plot files excluded
- Code cleaned for production use

### Statistics
| Category | Debug Branch | Indranil Branch | Difference |
|----------|--------------|-----------------|------------|
| Total Code Files | ~170 | ~65 | -105 files |
| Debug Prints | ~200+ | ~20 | -180 prints |
| Profiler Files | ~40 | 0 | -40 files |
| Plot/Data Files | ~15 | 0 | -15 files |

---

## 2. Debug Prints Removed

### 2.1 `pyBall/OCL/GridFF.py`
**Total Debug Prints Removed**: 88

#### In `__init__()` method:
```python
# REMOVED from debug:
if DEBUG: print(" local_memory_per_workgroup() size=", local_size, " __local []  ", ...)
```

#### In `fit3D()` method:
```python
# REMOVED from debug:
print(f"GridFF::fit3D() nG={nG} nL={nL} nxyz={nxyz}")
print(f"GridFF::fit3D() gsh={gsh} lsh={lsh}")
# ... (multiple debug prints)
```

#### In `make_MorseFF()` method:
```python
# REMOVED from debug:
print(f"make_MorseFF() na={na} nG={nG}")
print(f"make_MorseFF() REQs shape: {REQs.shape}")
print(f"make_MorseFF() xyzq shape: {xyzq.shape}")
# ... (multiple debug prints)
```

#### In `makeCoulombEwald_slab()` method:
**NOTE**: Debug prints in this section were RETAINED per user request
```python
# RETAINED in Indranil (per user request):
print("\nDebug before slabPotential:")
print(f"V_Coul_buff exists with shape: {debug_before.shape}")
print(f"Range: {debug_before.min():.3f} to {debug_before.max():.3f}")
print(f"Verification buffer shape: {verify_vcoul.shape}")
print(f"Verification buffer range: [{verify_vcoul.min():.3f}, {verify_vcoul.max():.3f}]")
```

#### In `poisson()` method:
```python
# REMOVED from debug:
print(f"poisson() nG={nG} nL={nL}")
print(f"poisson() Vin shape: {Vin.shape}")
```

#### In `slabPotential()` method:
```python
# REMOVED from debug:
print(f"slabPotential() nz_slab={nz_slab}")
print(f"slabPotential() dipol={dipol}")
```

#### Profiler Function Removed:
```python
# REMOVED entirely from debug:
def profile_kernel(name, queue, kernel, global_size, local_size, *args, **kwargs):
    """Profile a kernel execution"""
    start_time = time.time()
    event = kernel(queue, global_size, local_size, *args, **kwargs)
    event.wait()
    end_time = time.time()
    duration = (end_time - start_time) * 1000
    if name not in kernel_times:
        kernel_times[name] = []
    kernel_times[name].append(duration)
    print(f"Kernel: {name}, Duration: {duration:.3f} ms")
    return event

def save_profiling_results(output_dir="firecore_profile_results"):
    # ... profiling save logic
```

#### Profile Kernel Calls Replaced:
```python
# BEFORE (debug):
profile_kernel("setMul", self.queue, self.prg.setMul, (nG,), (nL,), nxyz, Ref_buff, self.Gs_buff, np.float32(1.0))

# AFTER (Indranil):
self.prg.setMul(self.queue, (nG,), (nL,), nxyz, Ref_buff, self.Gs_buff, np.float32(1.0))
```

---

### 2.2 `pyBall/MMFF.py`
**Total Debug Prints Removed**: ~20

#### In `getBuffs()` function:
```python
# REMOVED from debug (Lines 924-943):
print("Debug: Library path:", lib._name)
print("Debug: Library path:", lib._name)
print("Debug: ndims pointer address:", hex(ctypes.cast(ndims_ptr, ctypes.c_void_p).value))
print("Debug: ndims array after conversion:", ndims)
print("Debug: Raw ndims array from C++:", ndims)
print("getBuffs(): nDOFs=%i nvecs=%i natoms=%i nnode=%i ncap=%i npi=%i nbonds=%i nvecs=%i ne=%i ie0=%i " % ...)
print("getBuffs(): natoms=%i nnode=%i ncap=%i npi=%i nbonds=%i nvecs=%i ne=%i ie0=%i " % ...)
print("getBuffs DONE")
```

#### In `init()` function:
```python
# REMOVED from debug (Line 1039):
print(f"Initializing with: xyz={xyz_name}, surf={surf_name}, smile={smile_name}")
```

---

### 2.3 `pyBall/tests/ocl_GridFF_new.py`
**Total Debug Prints Removed**: ~30

#### In `make_atoms_arrays()` function:
```python
# REMOVED from debug:
print("Debug: atoms loaded from", fname)
print("Debug: atypes =", atoms.atypes)
print("Debug: apos shape =", atoms.apos.shape)
print("Debug: qs =", atoms.qs)
print("Debug: REvdW =", REvdW)
print("Debug: REQs =", REQs)
print("Debug: xyzq =", xyzq)
```

#### In `test_gridFF_ocl()` function:
```python
# REMOVED from debug:
print("Debug: z_coords =", z_coords)
print("Debug: unique_zs =", unique_zs)
print("Debug: counts =", counts)
print("Debug: most_frequent_zs =", most_frequent_zs)
print("Debug: z0 =", z0)
print("Debug: grid.ns =", grid.ns)
print("Debug: grid.dg =", grid.dg)
print("Debug: g0 =", g0)
# ... (many more debug prints)
```

---

### 2.4 `pyBall/atomicUtils.py`
**Changes**: Minimal (only functional changes, no debug prints)

The only differences are the functional changes:
- Line 862: `bDict=False` → `bDict=True`
- Line 1019: `rec[0]` → `rec[1]`

---

## 3. Profiler Files Excluded

### 3.1 Files in `pyBall/OCL/` (EXCLUDED)
| File | Lines | Purpose |
|------|-------|---------|
| `GridFF_profiled.py` | ~1200 | Profiled version of GridFF |
| `profiler.py` | ~300 | Profiling utilities |

### 3.2 Files in `tests/tMMFF/` (EXCLUDED)
| File | Lines | Purpose |
|------|-------|---------|
| `grid_profiler.py` | ~200 | Grid profiling |
| `kernel_profiler.py` | ~150 | Kernel profiling |
| `profile_gridff.py` | ~180 | GridFF profiling |
| `run_test_GridFF_ocl_profile.py` | ~100 | Profiled OCL test |
| `run_test_GridFF_ocl_profiled.py` | ~120 | Profiled OCL test |
| `run_test_GridFF_ocl_profiled_direct.py` | ~90 | Direct profiled test |
| `run_with_enhanced_profiling.py` | ~150 | Enhanced profiling |
| `simple_profiler.py` | ~100 | Simple profiler |

### 3.3 Files in `tests/tMMFF/test_profiler/` (EXCLUDED - entire folder)
| File | Lines | Purpose |
|------|-------|---------|
| `basic_profiler.py` | 95 | Basic profiling |
| `detailed_profiler.py` | 120 | Detailed profiling |
| `direct_gridff_profiler.py` | 140 | Direct GridFF profiler |
| `direct_profiling.py` | 110 | Direct profiling |
| `grid_profiler.py` | 180 | Grid profiler |
| `gridff_monitor.py` | 200 | GridFF monitor |
| `gridff_path_profiler.py` | 160 | Path profiler |
| `gridff_profiler.py` | 220 | GridFF profiler |
| `inspect_pyopencl.py` | 80 | PyOpenCL inspection |
| `kernel_profiler.py` | 150 | Kernel profiler |
| `manual_profiling.py` | 130 | Manual profiling |
| `nvidia_opencl_profiler.py` | 170 | NVIDIA profiler |
| `opencl_env_profiler.py` | 140 | OpenCL env profiler |
| `opencl_kernel_profiler.py` | 160 | OpenCL kernel profiler |
| `opencl_profiler.py` | 190 | OpenCL profiler |
| `profile_direct.py` | 100 | Direct profiling |
| `profile_gpu_cpu.py` | 220 | GPU/CPU profiling |
| `profile_gpu_detailed.py` | 250 | Detailed GPU profiling |
| `profile_grid.py` | 180 | Grid profiling |
| `profile_gridff_cpu.py` | 259 | CPU GridFF profiling |
| `profile_gridff_detailed.py` | 300 | Detailed GridFF profiling |
| `profile_gridff_memory.py` | 376 | Memory profiling |
| `profile_gridff_ocl.py` | 145 | OCL GridFF profiling |
| `profile_opencl_kernels.py` | 170 | OpenCL kernel profiling |
| `profiled_test.py` | 69 | Profiled test |
| `profiler.py` | 531 | Main profiler |
| `run_gridff_with_paths.py` | 96 | GridFF with paths |
| `run_profiler.py` | 175 | Run profiler |
| `run_profiling.py` | 230 | Run profiling |
| `run_test_GridFF_ocl_profiled.py` | 61 | Profiled OCL test |
| `simple_profiler.py` | 264 | Simple profiler |
| `test_gpu_detection.py` | 153 | GPU detection test |

**Total Profiler Files Excluded**: ~40 files, ~6,000+ lines

---

## 4. Plot/Data Files Excluded

### 4.1 Debug-specific Files in `tests/tMMFF/` (EXCLUDED)
| File | Lines | Purpose |
|------|-------|---------|
| `debug_grid_shift.py` | ~150 | Debug grid shift |
| `debug_pipeline.py` | ~200 | Debug pipeline |
| `debug_slab_mapping.py` | ~180 | Debug slab mapping |
| `batch_plot.py` | ~120 | Batch plotting |
| `plot_numpy_object.py` | ~80 | Numpy object plotting |
| `trial_ploting.py` | 157 | Trial plotting |

### 4.2 Data-specific Plot Files in `cpp/common_resources/` (EXCLUDED)
| File | Purpose |
|------|---------|
| `NaCl_1x1_L1/plot_grid.py` | NaCl 1x1 L1 plotting |
| `NaCl_1x1_L3/plot_grid.py` | NaCl 1x1 L3 plotting |
| `NaCl_8x8_L3_ClHole/plot.py` | NaCl 8x8 ClHole plotting |
| `NaCl_8x8_L3_final/plot_grid.py` | NaCl 8x8 final plotting |
| `Na_0.9_Cl_-0.9_Cl_hole/plot.py` | Cl hole plotting |
| `Na_0.9_Cl_-0.9_Cl_hole_3/plot.py` | Cl hole 3 plotting |
| `Na_0.9_Cl_-0.9_Cl_hole_3/plot_grid.py` | Cl hole 3 grid plotting |
| `Na_0.9_Cl_-0.9_old_grod_step/plot_grid.py` | Old grid step plotting |

### 4.3 Throughput/Evaluation Files (EXCLUDED)
| File | Lines | Purpose |
|------|-------|---------|
| `tests/tGridFF/number_evaluations_vs_time/evaluation_vs_time.py` | 209 | Evaluation timing |
| `tests/tGridFF/plot_data.py` | ~500 | Data plotting |
| `tests/tGridFF/run_throughput_MD.py` | ~100 | Throughput MD |
| `tests/tGridFF/throughput_gui.py` | ~400 | Throughput GUI |
| `tests/tMMFFmulti/number_evaluations_vs_time/evaluation_vs_time.py` | 209 | Evaluation timing |
| `tests/tMMFFmulti/plot_data.py` | 1060 | Data plotting |
| `tests/tMMFFmulti/run_throughput_MD.py` | 47 | Throughput MD |
| `tests/tMMFFmulti/throughput_gui.py` | 880 | Throughput GUI |

---

## 5. File-by-File Detailed Changes

### 5.1 `pyBall/OCL/GridFF.py`

#### Lines 1-50: Imports and Global Variables
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| 3 | `import time` | (removed) | REMOVED |
| 4 | `import json` | (removed) | REMOVED |
| 6 | `import matplotlib.pyplot as plt` | (removed) | REMOVED |
| 8 | `import pyopencl.tools` | (removed) | REMOVED |
| 9 | `import matplotlib.pyplot as plt` | (removed) | REMOVED (duplicate) |
| 10-12 | `kernel_times = {}` | (removed) | REMOVED |

#### Lines 13-50: `profile_kernel()` function
| Status | Description |
|--------|-------------|
| REMOVED | Entire function (~40 lines) removed |

#### Lines 51-90: `save_profiling_results()` function
| Status | Description |
|--------|-------------|
| REMOVED | Entire function (~40 lines) removed |

#### Lines ~100-150: `__init__()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~105 | `DEBUG = False` | (removed) | REMOVED |
| ~110 | `if DEBUG: print(...)` | (removed) | REMOVED |
| ~115-130 | `try: ... except: ...` for kernel loading | Direct file open | SIMPLIFIED |

#### Lines ~400-460: `fit3D()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~405 | `print(f"GridFF::fit3D() nG={nG}...")` | (removed) | REMOVED |
| ~410 | `profile_kernel("setMul", ...)` | `self.prg.setMul(...)` | REPLACED |
| ~415 | `profile_kernel("setMul", ...)` | `self.prg.setMul(...)` | REPLACED |
| ~420 | `profile_kernel("BsplineConv3D", ...)` | `self.prg.BsplineConv3D(...)` | REPLACED |
| ~425 | `profile_kernel("BsplineConv3D", ...)` | `self.prg.BsplineConv3D(...)` | REPLACED |
| ~430 | `profile_kernel("move", ...)` | `self.prg.move(...)` | REPLACED |

#### Lines ~600-700: `projectAtoms_on_grid_d()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~605 | `print(f"projectAtoms_on_grid_d()...")` | (removed) | REMOVED |
| ~620 | `print(f"projectAtoms_on_grid_d() na={na}...")` | (removed) | REMOVED |

#### Lines ~800-900: `poisson()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~805 | `print(f"poisson() nG={nG}...")` | (removed) | REMOVED |
| ~820 | `profile_kernel("poissonW", ...)` | `self.prg.poissonW(...)` | REPLACED |

#### Lines ~900-990: `makeCoulombEwald_slab()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~949-984 | Debug verification prints | **RETAINED** | KEPT (per user request) |

#### Lines ~1000-1100: `make_MorseFF()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~1005 | `print(f"make_MorseFF() na={na}...")` | (removed) | REMOVED |
| ~1010 | `print(f"make_MorseFF() REQs shape...")` | (removed) | REMOVED |
| ~1015 | `print(f"make_MorseFF() xyzq shape...")` | (removed) | REMOVED |
| ~1050 | `profile_kernel("make_MorseFF_f4", ...)` | `self.prg.make_MorseFF_f4(...)` | REPLACED |

#### Lines ~1100-1120: `release_unused_buffs()` method
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~1105-1115 | `try: val.release() except: ...` | `if hasattr(val, 'release'): val.release()` | SIMPLIFIED |

---

### 5.2 `pyBall/MMFF.py`

#### Lines 922-972: `getBuffs()` function
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| 925 | `print("Debug: Library path:", lib._name)` | (removed) | REMOVED |
| 926 | `print("Debug: Library path:", lib._name)` | (removed) | REMOVED (duplicate) |
| 930 | `print("Debug: ndims pointer address:", ...)` | (removed) | REMOVED |
| 934 | `print("Debug: ndims array after conversion:", ndims)` | (removed) | REMOVED |
| 943 | `print("Debug: Raw ndims array from C++:", ndims)` | (removed) | REMOVED |
| 950 | `print("getBuffs(): nDOFs=%i nvecs=%i...")` | (removed) | REMOVED |
| 971 | `print("getBuffs(): natoms=%i nnode=%i...")` | (removed) | REMOVED |
| 972 | `print("getBuffs DONE")` | (removed) | REMOVED |

#### Lines 1035-1041: `init()` function
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| 1039 | `print(f"Initializing with: xyz=...")` | (removed) | REMOVED |

---

### 5.3 `pyBall/tests/ocl_GridFF_new.py`

#### Lines 197-237: `make_atoms_arrays()` function
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~200-210 | Multiple debug prints for atoms | (removed) | REMOVED |
| ~215-220 | Debug prints for REvdW, REQs | (removed) | REMOVED |

#### Lines 239-450: `test_gridFF_ocl()` function
| Line | Debug Branch | Indranil Branch | Change Type |
|------|--------------|-----------------|-------------|
| ~250-260 | Debug prints for z_coords | (removed) | REMOVED |
| ~270-280 | Debug prints for grid params | (removed) | REMOVED |
| ~350-352 | `print("Starting Coulomb...")` | **RETAINED** | KEPT |
| ~378-381 | `print("Starting Morse...")` | **RETAINED** | KEPT |

---

## 6. Summary: What Was Preserved vs Removed

### PRESERVED (Functional Code)
- All OpenCL kernel implementations
- All force field calculation methods
- All grid generation algorithms
- All scan functions (rigid, constrained, UFF)
- Essential status prints (e.g., "Starting Coulomb...", "Starting Morse...")
- Debug verification in `makeCoulombEwald_slab()` (per user request)

### REMOVED (Debug Artifacts)
- 88 debug print statements in `GridFF.py`
- 20 debug print statements in `MMFF.py`
- 30 debug print statements in `ocl_GridFF_new.py`
- `profile_kernel()` function and all profiling infrastructure
- `kernel_times` dictionary and profiling data collection
- `save_profiling_results()` function
- All `try-except` blocks replaced with `hasattr` checks where appropriate
- Duplicate imports (e.g., double `matplotlib.pyplot`)

### EXCLUDED (Separate Files)
- ~40 profiler files (~6,000+ lines)
- ~15 data-specific plot files
- ~8 debug-specific test files

---

## 7. Verification Checklist

To verify the Indranil branch has all functionality from debug:

### Core Functionality Tests
- [ ] Grid generation works (`generate_grid.py`)
- [ ] Coulomb potential calculation works
- [ ] Morse potential calculation works
- [ ] UFF scan works (`scan_rigid_uff`)
- [ ] Constrained scan works (`scan_constr`)
- [ ] Results match debug branch (within numerical precision)

### Expected Differences
- [ ] Less console output (debug prints removed)
- [ ] No profiling data generated
- [ ] Cleaner code structure

---

*Document generated for verification and audit purposes.*
*All line numbers are approximate and may vary slightly due to code modifications.*
