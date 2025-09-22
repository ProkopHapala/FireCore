
# UFF multi‑system (GPU/OpenCL) design notes

These notes capture the target design for running UFF over `nSystems` replicas in parallel on the GPU, and the action items to fix the current issues discovered by tests.

## Problem recap

- **Initial State:** The `scan()` test showed system 0 matching the CPU, while system 1 returned zero or garbage forces.
- **Root Causes Identified:**
  - Kernels were not launched with a 2D NDRange, so `get_global_id(1)` was invalid.
  - Buffer indexing was not multi-system aware, causing data overlaps.
  - Host-side packing of neighbor indices (`angNgs`, `dihNgs`) into `hneigh` was missing per-system offsets.
  - The `assembleForces_UFF` kernel used an incorrect base offset for the `a2f_indices` buffer.

## Relevant Files

- `/home/prokop/git/FireCore/cpp/common/molecular/UFF.h` - C++ implementation of UFF forcefield
- `/home/prokop/git/FireCore/cpp/common_resources/cl/UFF.cl` - OpenCL kernels for UFF forcefield
- `/home/prokop/git/FireCore/cpp/common/OpenCL/OCL_UFF.h` - OpenCL wrapper for UFF forcefield interfacing UFF.cl
- `/home/prokop/git/FireCore/cpp/common/molecular/MolWorld_sp3_multi.h` - Main application class for parallel simulations using MMFF and UFF implemented both on CPU and GPU. Interface both UFF.h and OCL_UFF.h. 
- `/home/prokop/git/FireCore/cpp/libs_OCL/MMFFmulti_lib.cpp` - extern "C" library interface for MolWorld_sp3_multi class, allowing it to be called from Python (avoiding C++ name mangling)
- `/home/prokop/git/FireCore/pyBall/MMFF_multi.py` - Python interface for MMFFmulti_lib exposing both MMFF and UFF multi-system forcefields (using both CPU and GPU implementations).
- `/home/prokop/git/FireCore/tests/tUFF/test_UFF_multi.py` - Test script for UFF multi-system forcefield ( compares CPU and GPU results, either single-point or relaxation).

## Target execution model

- Use a 2D NDRange for all kernels that operate per atom or per interaction: `global.x = nLocal * ceil(nItems / nLocal)`, `global.y = nSystems`.
- In kernels, obtain both indices:
  - `int iG   = get_global_id(0);`  // local item index within a system (atom or interaction id)
  - `int iSys = get_global_id(1);`  // system replica index in [0, nSystems)
- Compute per‑system offsets and index into the correct slice of each buffer:
  - Atoms:  `i0a = iSys * nAtoms;`
  - Bonds:  `i0b = iSys * nBonds;`
  - Angles: `i0ang = iSys * nAngles;`
  - Dihedrals: `i0dih = iSys * nDihedrals;`
  - Inversions: `i0inv = iSys * nInversions;`
  - H‑neigh: `i0h = iSys * (nAtoms*4);`
  - Fint:    `i0f = iSys * nf_per_system;` (exact CPU layout)
  - A2F:     `i0a2f_offs = iSys * nAtoms;` (for `a2f_offsets`/`counts`) and `i0a2f_inds = iSys * nf_per_system;` (for the `a2f_indices` buffer itself)
- All reads/writes use these base offsets:
  - `apos[i0a + ia]`, `fapos[i0a + ia]`, `neighs[i0a + ia]`, `neighCell[i0a + ia]`, `REQs[i0a + ia]`, etc.
  - `bonAtoms[i0b + ib]`, `bonParams[i0b + ib]`.
  - `angAtoms[i0ang + iang]`, `angNgs[i0ang + iang]`, ...
  - `dihAtoms[i0dih + id]`, `dihNgs[i0dih + id]`, `dihParams[i0dih + id]` (as float4), ...
  - `invAtoms[i0inv + ii]`, `invNgs[i0inv + ii]`, ...
  - `hneigh[i0h + ia*4 + slot]`.
  - `fint[i0f + ...]` as per CPU’s `i0` scheme for each interaction class.

## NDRange and OCL task configuration

- Ensure tasks are created as 2D when `nSystems > 1` (or always 2D for simplicity):
  - In `OCL_UFF::makeKernels()`, `newTask(name, program, /*nDim=*/2, ...)` for all UFF kernels.
  - In `OCL_UFF::setup_kernels()`, set `task->local = {nloc, 1, 1, 1}`, `task->global = {GX, nSystems, 1, 1}` where `GX` is rounded up to a multiple of `nloc`.
  - If `nSystems == 1`, 2D is also fine; `global.y = 1` is valid.

## Host/Device buffer layout and uploads

- Allocate GPU buffers sized for all systems (`nSystems * nPerSystem`).
- On the host, maintain packed, per‑system concatenated buffers:
  - Topology: `bonAtoms`, `angAtoms`, `dihAtoms`, `invAtoms`, `neighs`, `neighCell`, `neighBs`.
  - Parameters: `bonParams`, `angParams1`, `angParams2_w`, `dihParams`, `invParams`, `REQs`.
  - Aux: `a2f_offsets`, `a2f_counts`, `a2f_indices`, `pbc_shifts`, `lvecs`.
- Upload per‑system slices with matching host pointer and device offset (rule from memory):
  - Example: `upload(ibuff_bonAtoms, host_bonAtoms + i0b, nBonds, i0b)`.
  - Do this consistently for every buffer.

## Kernel specifics

- `evalBondsAndHNeigh_UFF` (atom‑centric):
  - Read `neighs[i0a + ia]`, `neighCell[i0a + ia]`, `apos[i0a + ia]`.
  - Use precomputed `neighBs[i0a + ia]` to fetch bond index `ib` (global within system), then fetch `bonParams[i0b + ib]`, `bonAtoms[i0b + ib]`.
  - Write `hneigh[i0h + ia*4 + slot]`.
  - Accumulate atomic forces into `fapos[i0a + ia]` and write any per‑bond contributions into `fint[i0f + (i0bon + ib*2) + {0,1}]` if needed.
- `evalAngles_UFF`, `evalDihedrals_UFF`, `evalInversions_UFF` (interaction‑centric):
  - Use `i0ang`, `i0dih`, `i0inv` to index atoms, params, and `hneigh`.
  - Write their force pieces into `fint[i0f + i0<component> + offset]`.
- `assembleForces_UFF` (atom‑centric):
  - Base offsets `i0a`, `i0a2f`, `i0f`.
  - For `ia` in `0..nAtoms`: read `off=a2f_offsets[i0a + ia]`, `cnt=a2f_counts[i0a + ia]`.
  - Accumulate over `cnt` entries at `a2f_indices[i0a2f_inds + off + k]` which index into `fint[i0f + idx]`.
  - Write into `fapos[i0a + ia]`.

## Debug instrumentation

- Keep the compact oneliners.
- Gate prints by `isys == IDBG_SYS`.
- For assemble kernel, also include `isys` and print `A2F` rows for the selected system only.

## CPU ↔ GPU parity in scan()

- The `scan()` GPU path should:
  - Pack/Upload topology/params for all systems once (bParams=true, blvec=true).
  - For each batch of configurations: pack positions for systems 0..nBatch-1 (bParams=false), Upload (bParams=false), Eval once, Download forces, Unpack per system to host UFF for copying out.
- Verify that `pack_uff_system()` populates A2F, neigh*, and params for each `isys` slice correctly before the first upload.

## Action items

1) Switch all UFF tasks to 2D NDRange in `OCL_UFF` (nDim=2) so `get_global_id(1)` is valid.
2) Add per‑system offsets to every UFF kernel and use them consistently for all buffers.
3) Ensure host `upload_uff()` uses matching host slice and device offset for every buffer.
4) Instrument `assembleForces_UFF` prints with `isys` and gate by `IDBG_SYS`.
5) Re‑run tests with `IDBG_SYS=1` and confirm non‑zero forces for system 1.
6) Validate A2F mapping tables per system (offsets/counts/indices) are non‑empty and consistent.

## Testing plan

- Unit: create a tiny 2‑atom bond system; duplicate into `nSystems=2`; run only bonds and verify both systems produce identical non‑zero forces; print both systems.
- Integration: `tests/tUFF/test_UFF_multi.py --nsys 2 --nconf 2 --use-scan 1` on formic acid and xylitol.
- Stress: increase `nSystems` to 8 and `nconf` to 8; spot‑check `IDBG_SYS` across different systems.

## Notes

- It’s acceptable to always set tasks as 2D; `global.y=1` works for single system.
- Do not rely on implicit per‑system strides inside kernels; compute offsets explicitly from `iSys`.
- Keep debug prints minimal and focused; use `IDBG_SYS` to isolate a problematic replica.

---

## Lab Book / Progress Log

### Session 1: Initial Debugging (Bonds & Angles)

**Initial State:**
- `tests/tUFF/test_UFF_multi.py` with `nSystems=2` showed system 0 matching CPU, but system 1 forces were zero or garbage.
- Kernels were not using 2D NDRange, and buffer indexing was not multi-system aware.

**Key Fixes & Insights:**
1.  **Enabled 2D NDRange:** Switched all UFF kernels in `OCL_UFF.cpp` to use `nDim=2` to get a valid `isys = get_global_id(1)`.
2.  **Per-System Host Packing:** Implemented `pack_uff_system()` in `MolWorld_sp3_multi.h` to create host-side concatenated buffers for all topology and parameters. This resolved the issue of systems overwriting each other's data.
3.  **Corrected `hneigh` Indexing:** The `angNgs`, `dihNgs`, and `invNgs` arrays contain indices into the `hneigh` buffer. On the host, these indices were corrected to include the per-system offset `i0h = isys * (nAtoms*4)` before being packed.
4.  **Fixed `assembleForces_UFF` Indexing:** The `a2f_indices` buffer was being indexed with a base of `isys * nAtoms`, which was incorrect. The correct base is `isys * nf_per_system`. This was the final fix needed to get the `['bonds', 'angles']` test to pass.
5.  **Added Verbosity Guards:** To reduce log spam, verbose CPU debug prints in `UFF.h` were gated behind `DBG_UFF > 3`. The `scan()` function in `MMFFmulti_lib.cpp` was updated to set this verbosity level only for the system matching `IDBG_SYS` on the GPU.

**End State (Session 1):**
- The test `['bonds', 'angles']` now **PASSES** with a tolerance of `1e-3`.
- Both system 0 and system 1 produce correct forces.

### Session 2: Dihedrals & Inversions

**Initial State:**
- `['bonds', 'angles']` are working correctly.
- Enabling `['dihedrals']` causes the test to fail with `max|ΔF| ≈ 1.8`.

**Key Fixes & Insights:**
1.  **Dihedral Parameter Type Mismatch:** The host (`MolWorld_sp3_multi.h`) packs dihedral parameters as `float4`, but the `evalDihedrals_UFF` kernel was reading them as `float3`. This caused a data stride mismatch for `isys > 0`. The kernel was corrected to accept `__global float4* dihParams` and read the `.xyz` components.

**Current State (End of Session):**
- The test with `['bonds', 'angles', 'dihedrals']` still fails, but the error has changed. The `max|ΔF|` is now smaller (`~0.24`), indicating the parameter type fix was a step in the right direction.
- The `evalDihedrals_UFF` kernel appears to have an internal logic error in how it calculates forces, as the debug prints show non-physical force distribution (e.g., `fi` and `fj` are identical in one case).
- The `['inversions']` component also causes a failure.
- Next step is to analyze the force calculation logic inside the `evalDihedrals_UFF` and `evalInversions_UFF` kernels and compare it line-by-line with the CPU implementation in `UFF.h`.

### Session 3: Dihedral + Inversion parity (float, multi-system)

**Goal:** Achieve exact CPU↔GPU parity for dihedrals and inversions across `nSystems` replicas, while keeping single precision (float) on GPU.

**Changes (Dihedrals):**
1.  **Switch to Prokop formulation (CPU parity):**
    - Use unit h-vectors from `hneigh`: `h12=q12.xyz (ji)`, `h32=q32.xyz (jk)`, `h43=q43.xyz (kl)`.
    - Plane normals: `n123=cross(h12,h32)`, `n234=cross(h43,h32)`.
    - `il2_123=1/|n123|^2`, `il2_234=1/|n234|^2`, `inv_n12=sqrt(il2_123*il2_234)`.
    - `cs={ dot(n123,n234)*inv_n12, -dot(n123,h43)*inv_n12 }`, `csn = cs^n` via complex multiply, `n=(int)par.z`.
    - Energy: `E = V*(1 + d*csn.x)`; force scalar: `f = -V*d*n*csn.y`.
    - End-atom forces: `fp1 = n123*(-f*il2_123*q12.w)`, `fp4 = n234*( f*il2_234*q43.w)`.
    - Central atoms (conservation of angular momentum around jk):
      - `c123 = dot(h32,h12)*(q32.w/q12.w)`, `c432 = dot(h32,h43)*(q32.w/q43.w)`.
      - `fp3 = (-c123)*fp1 + (-c432 - 1)*fp4`, `fp2 = (c123 - 1)*fp1 + (c432)*fp4`.
    - Debug phi uses the same normals and `clamp(...,-1.0f,1.0f)`.
2.  **Type fixes (GPU compilers are strict):**
    - All literals are float-suffixed: `1.0f`, `1e-30f`; `clamp` uses float overload.
    - Avoided mixing float vectors with double scalars.
3.  **Topology buffer typing:**
    - `dihNgs` changed to `__global int4*` to match host packing and avoid stride mismatch.
    - `dihParams` is `__global float4*`, read `.xyz`.
4.  **Multi-system indexing:**
    - Use `isys=get_global_id(1)` and per-system bases (`i0dih`, `i0h`, `i0f`, etc.) consistently.
    - Optional 1–4 NB subtraction remains off by default; if enabled, we will add `natoms` to kernel args and index REQs/apos with `i0a`.

**Changes (Inversions):**
1.  **Match `UFF.h::evalInversion_Prokop` exactly in float:**
    - Plane normal `n123=cross(q21.f, q31.f)`, `il123=1/|n123|`; `s=-dot(n123,q41.f)`, `c=sqrt(1-s*s+eps)`.
    - Params: `par={K,c0,c1,c2}`; `E=K*(c0 + c1*c + c2*cos(2w))` via complex multiply; `f = -K*(c1*s + 2*c2*sin(2w))/c`.
    - Forces: `fp4=fq41*n123 + s*fq41*q41.f`, `tq=s*fi123*n123 + fi123*q41.f`, `fp2=cross(q31,tq)*q21.e`, `fp3=cross(tq,q21)*q31.e`, `fp1=-(fp2+fp3+fp4)`.
2.  **Single precision throughout:** float math and literals; no doubles.

**Results:**
- Parity achieved for `['bonds','angles','dihedrals','inversions']` across `nSystems=10` (tested via `tests/tUFF/test_UFF_multi.py`).
- GPU debug gated by `IDBG_SYS` matches CPU prints for the selected system.

## Rules & Conventions (for future contributors)

- __Float-only on GPU__: Do everything in single precision in OpenCL. No `double` in kernels (slow on consumer GPUs, different math paths). Use float-suffixed literals (`1.0f`, `1e-14f`), float overloads (e.g., `clamp(x,-1.0f,1.0f)`).
- __Follow CPU Prokop math exactly__: Mirror `UFF.h` “Prokop” variants (`evalDihedral_Prokop`, `evalInversion_Prokop`) line-by-line. Use the same variable names/flow whenever possible to avoid drift.
- __Multi-system indexing__: Always compute per-system bases: `i0a,i0b,i0ang,i0dih,i0inv,i0h,i0f,i0a2f`. Index every buffer with its correct base.
- __Host ↔ Device slicing__: When uploading per-system slices, apply the same offset to the host pointer and the device buffer offset. Example: `upload(ibuff_dihAtoms, host + i0dih, nDihedrals, i0dih)`.
- __Topology buffer types__: Use `int4` for `*_Atoms`/`*_Ngs` that are 4-tuples; `float4` for params that are padded. Kernels should read `.xyz` or `.w` as needed.
- __Debugging__: Gate prints by `IDBG_SYS` and print only 1 system per kernel call. Keep oneliners to avoid interleaving.

## Future optimization opportunities

- __Local memory tiling for `hneigh`__: For interaction-centric kernels (angles/dihedrals/inversions), prefetch `hneigh` for the subset of atoms referenced by the current work-group into `__local` memory to reduce global reads.
- __Vectorized loads__: Use `float4`/`int4` aligned loads uniformly to reduce memory transactions (already done for params/atoms/ngs).
- __Occupancy & register pressure__: Tune work-group size (e.g., 64 vs 32) to balance occupancy against register usage in kernels with heavier math (dihedrals/inversions).
- __Kernel fusion (selective)__: Optionally fuse angle+dihedral+inversion into a single pass writing to `fint` when it improves cache locality; keep bonds separate due to different access pattern.
- __Math micro-opts__: Precompute `q*.w` reciprocals or reuse dot-products; consider fast-math flags if acceptable (ensure parity first).
- __Optional 1–4 NB subtraction__: If enabled, add `natoms` to kernel args and use `i0a` for `REQs/apos` indexing to maintain multi-system correctness.
- __Better debug controls__: Make `IDBG_SYS` and per-component `IDBG_*` runtime kernel args instead of macros to avoid rebuilds.

## Testing checklist (regression)

- `['bonds']` only, `nSystems={1,2,10}`.
- `['bonds','angles']`, parity within `1e-3`.
- `['dihedrals']` only, compare phi and force components per interaction.
- `['inversions']` only, compare `w`, `fi..fl` per interaction.
- Full set `['bonds','angles','dihedrals','inversions']` at `nSystems=10`, `nconf=10`.

## GPU Execution and Initialization Logic for UFF

This section details the initialization process and the runtime logic for evaluating forces on the GPU, especially concerning non-bonded and molecule-substrate interactions.

### Control Flags

The behavior of the GPU execution is controlled by several boolean flags in the `MolWorld_sp3_multi` class:

- `bUFF`: When `true`, the system uses the UFF force field and the `OCL_UFF` OpenCL handler.
- `bSurfAtoms`: When `true`, enables interactions with a substrate (surface).
- `bGridFF`: When `true` (and `bSurfAtoms` is also true), the substrate interaction is calculated by interpolating a pre-computed force field grid.
- `bNonBonded`: When `true`, enables pairwise Lennard-Jones and Coulomb interactions between atoms in the molecule.

### Initialization Flow

The initialization of substrate interactions is designed to be robust, even if the OpenCL context is not immediately available when the substrate is loaded.

1.  **`loadSurf()`**: When a surface is loaded, the `loadSurf` method is called. This triggers the `initGridFF` method.
2.  **`initGridFF()`**: This method is responsible for loading the grid parameters. Crucially, it checks if the required OpenCL context (`uff_ocl->context`) is available.
    - **Context Ready**: If the context exists, it immediately calls `surf2ocl_uff()` to prepare and upload the substrate geometry to the GPU.
    - **Context Not Ready**: If the context is not available (e.g., `loadSurf` is called before `init()`), it sets the `bGridFF_pending` flag to `true` and defers the GPU upload.
3.  **`completeGridFFInit()`**: This method is called later in the main `init()` function. If `bGridFF_pending` is true, it now calls `surf2ocl_uff()`, completing the deferred initialization.

### Substrate Interaction Modes

When `bSurfAtoms` is enabled, there are two mutually exclusive modes for calculating molecule-substrate interactions, controlled by the `bGridFF` flag:

1.  **GridFF Interpolation Mode (`bGridFF = true`)**
    - This mode relies on a pre-calculated force field grid stored on disk (e.g., `.lvs`, `.xsf` files).
    - The `initGridFF` function will attempt to load these files. **For UFF, on-the-fly grid generation is not supported.** If the grid files are not found, the program will print a fatal error and exit.
    - At runtime, the `run_uff_ocl` function enqueues the `task_NBFF_Grid_Bspline` kernel to calculate forces by interpolating the grid texture.

2.  **Direct Pairwise Mode (`bGridFF = false`)**
    - In this mode, substrate interactions are calculated directly (on-the-fly) between the molecule's atoms and the substrate's atoms.
    - The `setup_UFF_ocl` function prepares the `task_SurfAtoms` by calling `uff_ocl->getSurfMorse(...)`.
    - At runtime, `run_uff_ocl` enqueues the `task_SurfAtoms` kernel to perform the pairwise calculation.

### Runtime Logic in `run_uff_ocl`

The `run_uff_ocl` function has been updated to orchestrate the evaluation of all force components in the correct order:

1.  **Clear Forces**: `task_clear_fapos` is enqueued to zero out force accumulators.
2.  **Covalent Forces**: `uff_ocl->eval(false)` is called to enqueue the kernels for bonds, angles, dihedrals, and inversions.
3.  **Non-Covalent Forces**: Based on the control flags, the appropriate non-covalent kernel is enqueued:
    - If `bSurfAtoms` is true:
        - If `bGridFF` is true, `task_NBFF_Grid_Bspline` is used.
        - If `bGridFF` is false, `task_SurfAtoms` is used.
    - If `bSurfAtoms` is false but `bNonBonded` is true, the standard `task_NBFF` is enqueued for intramolecular non-bonded interactions.
4.  **Update Positions**: Finally, `task_updateAtoms` is enqueued to update particle positions and velocities based on the total calculated forces.
