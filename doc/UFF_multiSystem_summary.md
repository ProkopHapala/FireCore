
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

# UFF OpenCL Multi-System Implementation: Final Summary and Future Reference

This document summarizes the critical insights, design rules, and caveats discovered during the debugging and verification of the multi-system Universal Force Field (UFF) implementation on OpenCL. It serves as a permanent reference for future development, distilling the lessons learned from the extensive debugging log.

## 1. Core Architecture: The Multi-System Execution Model

The fundamental requirement for correct multi-system (multi-replica) GPU execution is the strict segregation of data per system.

### A. Kernel Launch Configuration (NDRange)

*   **Always Use 2D Ranges:** All kernels (covalent, non-bonded, utility) must be launched with a 2D NDRange when `nSystems >= 1`.
    *   `global.x`: Parallelization within a system (atoms or interactions). Rounded up to the nearest `local.x`.
    *   `global.y`: The system index (`nSystems`).
*   **Local Workgroup Size:**
    *   `local.x`: Typically 32 or 64.
    *   `local.y`: Must be `1`.

### B. Universal Kernel Indexing

Every kernel must retrieve both the intra-system index and the system ID:

```c
int iG   = get_global_id(0);  // Atom or interaction index within the system
int iSys = get_global_id(1);  // System replica index [0, nSystems)
```

### C. Per-System Base Offsets (The Golden Rule)

**Never access a buffer without adding the correct per-system offset.** These offsets must be calculated in every kernel:

*   **Atoms/Positions/Velocities:** `i0a = iSys * nAtoms;`
*   **Bonds:** `i0b = iSys * nBonds;`
*   **Angles:** `i0ang = iSys * nAngles;`
*   **Dihedrals:** `i0dih = iSys * nDihedrals;`
*   **Inversions:** `i0inv = iSys * nInversions;`
*   **H-Neigh (4 per atom):** `i0h = iSys * (nAtoms * 4);`
*   **Forces (fint):** `i0f = iSys * nf_per_system;` (Matches CPU layout)
*   **Atom-to-Force (A2F) Mapping:**
    *   Offsets/Counts: `i0a2f_offs = iSys * nAtoms;`
    *   Indices Buffer: `i0a2f_inds = iSys * nf_per_system;` (Crucial fix: Not `nAtoms`!)

**Example (Assemble Forces Kernel):**
To read the force index for atom `ia`:
`int off = a2f_offsets[i0a + ia];`
`int idx = a2f_indices[i0a2f_inds + off + k];`
`float4 f = fint[i0f + idx];`

### D. Host-Side Management

*   **Concatenated Buffers:** Host-side data (topology, parameters, positions) must be packed into single, large concatenated buffers (`System0 | System1 | ...`).
*   **Offset Packing (Crucial):** When packing neighbor indices (e.g., `angNgs`, `dihNgs` which point into `hneigh`), the host code must add the `i0h` offset *before* packing.
*   **Slice Uploads:** Uploads must correctly specify both the host pointer offset and the device buffer offset:
    `upload(device_buffer, host_ptr + i0_offset, count, i0_offset);`

---

## 2. Covalent Force Field Implementation (Bonds, Angles, Dihedrals, Inversions)

Parity with the CPU reference is paramount.

### A. Precision and Math

*   **Float Only on GPU:** Use single precision (`float`, `float4`) exclusively. `double` is slow and causes math path divergence.
*   **Explicit Literals:** Always use `f` suffixes (`1.0f`, `1e-30f`).
*   **Exact CPU Match:** GPU kernels must mirror the "Prokop" formulation in `UFF.h` line-by-line.

### B. Buffer Type Matching (Host/Kernel)

Mismatches cause stride errors in multi-system mode.

*   **Parameters (e.g., `dihParams`):** If packed as `float4` on the host for alignment, the kernel must accept `__global float4*` and use `.xyz`. A `float3` kernel argument will misalign for `iSys > 0`.
*   **Topology (e.g., `dihNgs`):** If packed as 4-tuples on the host, use `__global int4*` in the kernel.

---

## 3. Non-Bonded and GridFF Implementation (Context Management)

The interaction between UFF/MMFF contexts and GridFF (substrate forces) was the most complex debugging area.

### A. Context Independence (UFF vs. MMFF)

`OCL_UFF` and `OCL_MM` must be treated as separate, independent silos.

*   **No Sharing:** Do not copy buffer indices (`ibuff...`) or parameters (GridFF steps, origins) from the `ocl` (MMFF) context to the `uff_ocl` context.
*   **Initialization:** If UFF is active, `uff_ocl` must be initialized directly from the source data.

### B. GridFF: The Parity Requirement

The `OCL_UFF` implementation of GridFF must exactly mirror `OCL_MM`.

*   **Variable Parity:** `OCL_UFF.h` must contain all GridFF-related members present in `OCL_MM.h` (e.g., `grid_p0`, `grid_invStep`, `cl_dGrid`, `cl_diGrid`, `itex_BsplinePLQ`, etc.).
*   **`setGridShape` Method:** `OCL_UFF` must implement `setGridShape` identically to `OCL_MM` to correctly calculate inverse matrices and set the grid origin (`grid_p0`).

### C. Critical GridFF/Non-Bonded Fixes

If GridFF or Non-Bonded forces fail:

1.  **GridFF Buffer Overwrite (Session 12 Fix):** Ensure `MolWorld_sp3_multi.h` does *not* overwrite the UFF GridFF buffer index (`uff_ocl->ibuff_BsplinePLQ`) with the potentially uninitialized MMFF index (`ocl.ibuff_BsplinePLQ`).
2.  **GridFF Parameter Passing (Session 11 Fix):** In `MolWorld_sp3_multi.h`, configure UFF GridFF via `uff_ocl->setGridShape(gridFF.grid)`, not by copying from `ocl`.
3.  **Duplicate Interactions (Session 7 Fix):** The non-bonded kernel (`getNonBond`) is susceptible to logic errors that cause double-counting of atom pairs, often visible via debug prints showing the same pair `(i, j)` twice.

---

## 4. Debugging and Testing Protocol

### A. Effective Debugging

*   **Gate Prints:** Always gate kernel debug prints with `if (iSys == IDBG_SYS)`.
*   **CPU Parity Debugging:** To see CPU non-bonded debug prints, ensure `dbg_sys` (in `MolWorld_sp3_multi`) is set to the target system index (e.g., `0`), not `-1`.
*   **A2F Verification:** If forces are zero, verify `a2f_counts` and `a2f_indices` (with correct offsets) are populated.

### B. Regression Testing Checklist

To validate future changes, run tests that cover:

1.  **Multi-System Covalent:** `nSystems >= 2`, verify bonds, angles, dihedrals, and inversions independently and together.
2.  **Non-Bonded Parity:** `nSystems >= 2`, isolate `getNonBond` kernel.
3.  **GridFF Parity:** `nSystems >= 2`, isolate `getNonBond_GridFF_Bspline` kernel.

---

## 5. Summary of Key Takeaways

1.  **The Offset is Everything:** Multi-system failures are almost always due to incorrect or missing `iSys * N` offsets in kernels or host packing.
2.  **Type Strictness:** GPU compilers and memory layouts are unforgiving; host `float4` must meet kernel `float4`.
3.  **Context Isolation:** UFF and MMFF OpenCL contexts must not share resources or initialization paths.
4.  **Trust the `setGridShape`:** For GridFF, the `setGridShape` method in the OCL helper class is the single source of truth for coordinate transformation.

Of course. Here is a comprehensive summary that preserves the essential, permanent information from your debugging log, structured as a reference guide for future development. It omits the temporary states and repeated findings, focusing on the final, correct implementation and the lessons learned.