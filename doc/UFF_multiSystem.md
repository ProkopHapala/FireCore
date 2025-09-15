---
description: UFF multi-system GPU design and debugging plan
---

# UFF multi‑system (GPU/OpenCL) design notes

These notes capture the target design for running UFF over `nSystems` replicas in parallel on the GPU, and the action items to fix the current issues discovered by tests.

## Problem recap

- **Initial State:** The `scan()` test showed system 0 matching the CPU, while system 1 returned zero or garbage forces.
- **Root Causes Identified:**
  - Kernels were not launched with a 2D NDRange, so `get_global_id(1)` was invalid.
  - Buffer indexing was not multi-system aware, causing data overlaps.
  - Host-side packing of neighbor indices (`angNgs`, `dihNgs`) into `hneigh` was missing per-system offsets.
  - The `assembleForces_UFF` kernel used an incorrect base offset for the `a2f_indices` buffer.

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
