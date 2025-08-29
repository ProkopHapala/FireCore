# Status: Verifying eFF Force Calculations (CPU vs. GPU)

## 1. Overview

The primary objective is to validate that the pyOpenCL (GPU) implementation of the Electron Force Field (eFF) produces results identical to the reference C++ (CPU) implementation. This document tracks the plan, actions taken, and current findings.

## 2. Key Files

*   **CPU**: `cpp/common/molecular/eFF.h` (Logic), `cpp/libs/Molecular/eFF_lib.cpp` (C-Interface), `pyBall/eFF.py` (Python Wrapper)
*   **GPU**: `pyBall/OCL/eFF_ocl.py` (Python Host), `cpp/common_resources/cl/eFF.cl` (OpenCL Kernel)

## 3. Actions Taken

1.  **Test Script Created**: A script, `tests/tEFF/test_ocl_vs_cpu.py`, was created to programmatically compare the two implementations. It loads a test H₂ molecule, runs both CPU and GPU calculations, and computes the maximum absolute difference between the resulting force arrays.

2.  **GPU Force Output**: The OpenCL pipeline was modified to retrieve force arrays. This involved adding a `force_buff` to `pyBall/OCL/eFF_ocl.py` and a corresponding `__global float4* fout` argument to the `localMD` kernel, allowing the script to read forces directly from the GPU.

3.  **Pairwise Debug Prints**: To facilitate detailed debugging, `printf` statements were added to the core logic of both implementations, guarded by debug flags. The consistent format is: `CPU/GPU INTERACTION(i,j) TYPE: (fx,fy,fz) | fsi,fsj`.

4.  **Kernel Logic Correction**: The initial tests revealed a major bug in the `localMD` kernel's force accumulation.
    *   **The Problem**: The kernel functions `getCoulombGauss` and `getPauliGauss_New` return a `float4` containing `{Energy, dE/dr, dE/ds_i, dE/ds_j}`. The original kernel code was incorrectly adding this entire `float4` directly to the force vector accumulator (e.g., `forcei += fg`). This meant energy and force derivatives were being erroneously summed into the force components.
    *   **The Correction**: The accumulation logic was changed to correctly calculate the force vector from the radial derivative and add the components to the appropriate accumulators.
        *   **Before**: `forcei += fg;`
        *   **After**: `forcei.xyz += dR * fg.y; forcei.w += fg.z;` (where `fg.y` is `dE/dr` and `fg.z` is `fsi`).

5.  **Shell Redirection**: The `run.sh` script was modified from `2>ERR` to `2>&1` to merge the `stderr` and `stdout` streams, ensuring all diagnostic messages from the C++ library and OpenCL driver are captured in the output file.

## 4. Current Status & Findings

**Result: FAIL**

The test script successfully runs, but it continues to report a large difference between the CPU and GPU forces. The correction of the kernel logic significantly changed the GPU results, but they still do not match the CPU reference.

**Analysis & Theories for Discrepancy:**

The most significant remaining issue is the fundamental difference in how the two implementations handle pairwise force accumulation due to their serial vs. parallel nature.

1.  **Primary Theory: Asymmetric Force Application in GPU Kernel.**
    *   **CPU (Serial)**: The C++ code iterates through pairs `(i, j)`. When it calculates the force `f` between them, it correctly applies Newton's third law by adding `+f` to particle `i` and `-f` to particle `j` within the same loop iteration. This is possible because the operations are sequential.
    *   **GPU (Parallel)**: The `localMD` kernel uses a simple parallel model where each thread is responsible for one particle `i`. When thread `i` calculates its interaction with particle `j`, it can easily add the force to its own `forcei` accumulator. However, it **cannot safely write to the force accumulator for particle `j`**, as another thread is responsible for that memory location. Doing so would create a race condition.
    *   **Consequence**: The current GPU kernel only ever calculates one half of the interaction. The reaction force on particle `j` is never accumulated, leading to fundamentally incorrect total forces.

2.  **Secondary Theory: Incorrect Size-Force (`fsj`) Handling.**
    *   This is a direct result of the primary problem. The force on the size of particle `j` (`fsj`) is calculated by thread `i` but is never applied because thread `i` cannot write to particle `j`'s data. My correction `forcei.w += fg.z;` only accumulates the `fsi` component. This explains the large, persistent error in the 4th component (the size force) of the electron force vectors.

## 5. Next Steps

The immediate task is to resolve the asymmetric force application in the `localMD` kernel. A simple fix is not possible without redesigning the kernel's parallel strategy.

*   **Recommended Action**: Before attempting a complex kernel rewrite (e.g., using atomic operations or a two-pass approach), the **pairwise debug prints should be used to verify that the single-sided force contributions are correct**. We need to run a test where we manually inspect the `printf` output from both CPU and GPU for a single pair interaction (e.g., `CPU EE(0,1) Coul` vs. `GPU EE(0,1) Coul`). If these match, it confirms the physics equations are correct and the problem lies solely in the accumulation of the reaction force. This is a crucial step to isolate the problem and prevent debugging the wrong part of the code.

---

## 6. Lab Book: Latest Findings (2025-08-29)

**Observation:** The `run.sh` script was redirecting `stderr` to a separate `ERR` file, causing the main `OUT_ocl_vs_cpu` log to be incomplete. All diagnostic messages from the C++ library and OpenCL driver were being hidden.

**Action Taken:** Modified `run.sh` to redirect `stderr` to `stdout` (`2>&1`). This merges all output streams, providing a complete log for analysis in `OUT_ocl_vs_cpu`.

**Detailed Analysis of Kernel Correction:**

The `localMD` kernel in `eFF.cl` was indeed accumulating forces incorrectly. The fix involved changing how the return values from `getCoulombGauss` and `getPauliGauss_New` are used.

*   **The Problem Revisited**: These functions return a `float4` containing `{E, fr, fsi, fsj}`, which represent `{Energy, dE/dr, dE/ds_i, dE/ds_j}`. The original code was adding this `float4` directly to the `forcei` vector, mixing incompatible physical quantities.

*   **The Correction Implemented**: The code was changed to correctly compute the force vector from the radial derivative (`fr`) and accumulate the appropriate components.
    *   **Before**: `forcei += fg;`
    *   **After**: `forcei.xyz += dR * fg.y;` and `forcei.w += fg.z;` (for the force on the size of particle `i`).

**Current Hypothesis for Discrepancy:**

Even with the corrected force calculation, the test **still fails**. The analysis in Section 4 remains the most likely explanation: the GPU kernel's parallel design prevents it from applying the reaction force to the second particle in a pair. The C++ code does `force_j -= f` while the GPU code does not, and cannot in its current form. This asymmetric application of forces is the primary bug.

**Additional Note on ASAN:** The logs show `CompilerWarning`s when PyOpenCL is used in a process instrumented with ASAN (`LD_PRELOAD`). This is likely due to the OpenCL driver performing low-level memory operations that ASAN interprets as unsafe. While this is worth noting, it is considered a separate issue from the core numerical discrepancy.
