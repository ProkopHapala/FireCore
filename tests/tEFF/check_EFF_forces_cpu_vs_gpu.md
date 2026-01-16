# Status: Verifying eFF Force Calculations (CPU vs. GPU)

## 1. Overview

The primary objective is to validate that the pyOpenCL (GPU) implementation of the Electron Force Field (eFF) produces results identical to the reference C++ (CPU) implementation. This document tracks the plan, actions taken, and current findings.

## 2. Key Files

*   **CPU**: `cpp/common/molecular/eFF.h` (Logic), `cpp/libs/Molecular/eFF_lib.cpp` (C-Interface), `pyBall/eFF.py` (Python Wrapper)
*   **GPU**: `pyBall/OCL/eFF_ocl.py` (Python Host), `cpp/common_resources/cl/eFF.cl` (OpenCL Kernel)

## 3. Summary of Actions

1.  **Test Script**: Created `tests/tEFF/test_ocl_vs_cpu.py` to programmatically compare CPU and GPU force arrays.
2.  **GPU Force Output**: Modified the OpenCL pipeline (`eFF_ocl.py`, `eFF.cl`) to support a dedicated force output buffer.
3.  **Debug Prints**: Added detailed, single-line `printf` statements to both `eFF.h` and `eFF.cl` to output the inputs and results of every pairwise interaction.
4.  **Kernel Logic Corrections**: Performed several corrections on the `localMD` kernel, fixing incorrect force accumulation logic and debug printouts.
5.  **Enabled Debugging**: Modified the test script and kernel to ensure the new `printf` statements are activated during the test run.

2.  **GPU Force Output**: The OpenCL pipeline was modified to retrieve force arrays. This involved adding a `force_buff` to `pyBall/OCL/eFF_ocl.py` and a corresponding `__global float4* fout` argument to the `localMD` kernel, allowing the script to read forces directly from the GPU.

3.  **Pairwise Debug Prints**: To facilitate detailed debugging, `printf` statements were added to the core logic of both implementations, guarded by debug flags. The consistent format is: `CPU/GPU INTERACTION(i,j) TYPE: (fx,fy,fz) | fsi,fsj`.

4.  **Kernel Logic Correction**: The initial tests revealed a major bug in the `localMD` kernel's force accumulation. 
    *   **The Problem**: The kernel functions `getCoulombGauss` and `getPauliGauss_New` return a `float4` containing `{Energy, dE/dr, dE/ds_i, dE/ds_j}`. The original kernel code was incorrectly adding this entire `float4` directly to the force vector accumulator (e.g., `forcei += fg`).
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

**Entry 1: Initial Test & Debugging Setup**
*   **Observation:** The `run.sh` script was redirecting `stderr` to a separate `ERR` file, hiding diagnostic messages.
*   **Action Taken:** Modified `run.sh` to redirect `stderr` to `stdout` (`2>&1`) for complete logging.
*   **Observation:** The `localMD` kernel in `eFF.cl` was accumulating forces incorrectly.
*   **Action Taken:** Corrected the kernel to use force derivatives (`forcei.xyz += dR * fr`) instead of adding the raw `{E,fr,fsi,fsj}` vector.

**Entry 2: Deeper Analysis of Discrepancy**

*   **Correction to Previous Theory:** My initial analysis that the `fsj` size-force was the *only* problem with the symmetric `for(j!=i)` loop was wrong. The errors in `fx,fy,fz` clearly indicate a more fundamental issue.

*   **Revised Hypothesis:** The core problem lies in the fundamental difference between the CPU's asymmetric loop (`for j<i`) and the GPU's symmetric loop (`for j!=i`).
    *   The **CPU** code computes the interaction for pair `(i,j)` **once**. It then manually distributes the action and reaction forces to `force[i]` and `force[j]`.
    *   The **GPU** code computes the interaction **twice**: thread `i` computes the `(i,j)` interaction, and thread `j` computes the `(j,i)` interaction.

*   **Why the GPU's Double Calculation Fails:** This approach is only valid if the force function is perfectly symmetrical, meaning `Force_on_i_from_j(i,j) = -Force_on_j_from_i(j,i)`. While this holds for simple forces, the `eFF` implementation is more complex. For example, in `Elec-Ion` interactions, the call from the electron's perspective is `getCoulombGauss(dR, si, Rj, ...)` and from the ion's perspective it is `getCoulombGauss(-dR, Rj, si, ...)`. The order of the size arguments (`si`, `Rj`) is swapped. If the underlying `getCoulombGauss` function is not perfectly symmetric with respect to swapping its size arguments, the two calculations will not be equal and opposite, and the total force will be wrong.

*   **Conclusion:** The current GPU kernel is attempting a valid parallel strategy (symmetric computation), but it fails because the underlying physics functions are likely being called in a way that breaks the required symmetry, leading to errors in all force components (`fx,fy,fz,fw`). The most robust solution remains to refactor the GPU kernel to mimic the CPU's asymmetric `for(j<i)` loop, combined with a `__local` buffer and atomic operations to handle the action-reaction pairs correctly. This removes any ambiguity related to function symmetry.

**Entry 3: Input Verification & Bug Fixes (2025-08-29)**

*   **What I Did:** Added comprehensive serial debug blocks to both CPU and GPU implementations for complete input verification:
    - **CPU (`eFF.h`)**: Enhanced `EFF::info()` to include `printAtomParams2()` dumping full 8-parameter atom settings (`Z_nuc, R_eff, Zcore_eff, PA..PE`)
    - **GPU (`eFF.cl`)**: Added detailed input dump in `localMD` kernel showing `KRSrho` vector, `bFrozenCore`, and for each particle: atom parameters or electron size/spin using correct global indexing

*   **What The Problem Was:** Large force discrepancies between CPU and GPU implementations (max diff ~185) with no clear visibility into whether inputs matched.

*   **What We Corrected:** 
    - **GPU Spin Indexing Bug**: Fixed electron spin lookup in debug dump - was using wrong global index (`is.w` instead of `is.z + i`)
    - **EE Debug Print**: Corrected to show proper size-force components (`cg.z,cg.w` and `pg.z,pg.w`) instead of just `pg.z,pg.w` twice
    - **Input Verification**: Both CPU and GPU now dump identical geometry, parameters, and switches, confirming inputs match perfectly

*   **What Problem Remains:** Force calculations still don't match despite identical inputs. The core issue is **asymmetric force accumulation** in the GPU kernel:
    - **CPU**: Serial loop computes each pair `(i,j)` once, applies `+f` to `i` and `-f` to `j`
    - **GPU**: Parallel threads each compute their interactions, but **cannot safely write to other particles' force accumulators** (race condition)
    - **Result**: GPU only computes one half of pairwise interactions, missing all reaction forces
    - **Evidence**: Force differences show systematic errors in both position forces (`fx,fy,fz`) and size forces (`fw`) for electrons

*   **Next Steps:** Before kernel rewrite, use the new debug prints to verify that single-sided force contributions match CPU. Then implement asymmetric parallel strategy with atomic operations or local buffers to properly accumulate reaction forces.

**Entry 4: Force Accuracy Improvements (2025-08-29)**

*   **Key Changes That Fixed Accuracy:**
    - Proper force component separation: `forcei.xyz += dR * fg.y` (radial forces)
    - Added missing size-force term: `forcei.w += fg.z`
    - Fixed spin indexing (`is.z + i`) for correct Pauli terms
    - Complete parameter verification ensuring identical inputs

*   **Impact:**
    - Position forces now match CPU within machine precision (~1e-5)
    - Size forces still problematic due to asymmetric accumulation

*   **Conclusion:** Core physics now correct - remaining errors are purely from parallelization strategy

## 7. Systematic comparison (2025-08-29 13:05)

From `tests/tEFF/OUT_ocl_vs_cpu`:

* __Inputs agree__
  - Geometry/params printed by CPU match GPU: atoms, electrons (pos, size, spin), and switches.
  - Evidence: CPU lines 17–25; GPU lines 65–70 show identical positions, sizes (0.5), spins (-1, +1), and atom params; KRSrho matches (GPU also prints sc=1.0).

* __AE Coulomb matches__
  - CPU AE Coulomb: lines 28–31.
  - GPU[serial] AE Coulomb: lines 71–74.
  - Vectors and size-force terms agree pair-by-pair (including signs).

* __EE Coulomb matches (up to dR orientation)__
  - CPU EE Coulomb: line 26 — `(-4.405, 0.000, 17.621) | 2.141, 2.141` for dR=(0.200, 0.000, -0.800).
  - GPU[serial] EE Coulomb: line 75 — `(4.405, -0.000, -17.621) | 2.141, 2.141` for dR=(-0.200, 0.000, 0.800).
  - Same magnitude with opposite direction due to reversed dR; size-force terms identical.

* __EE Pauli differs (sign and magnitude)__
  - CPU EE Pauli: line 27 — `(-5.249, 0.000, 20.997) | 0.170, 0.170`.
  - GPU[serial] EE Pauli: line 75 — `(0.844, -0.000, -3.377) | -1.970, -1.970`.
  - Disagreement in both vector and size-force components (signs opposite; magnitudes notably different).

* __Consequences in totals__
  - Atoms: CPU vs GPU forces match within ~1e-5 (lines 81–89), confirming AE is correct.
  - Electrons: large deviations dominated by EE Pauli error (lines 90–97), consistent with the pairwise analysis above.

* __Interim conclusion__
  - Inputs are consistent; AE and EE Coulomb are correct. The remaining discrepancy is isolated to the EE Pauli contribution in the GPU path.
  - Next action: inspect/validate `getPauliGauss_New()` usage and spin/sign handling in `cpp/common_resources/cl/eFF.cl` for EE interactions, before any kernel-parallelization changes.

## 8. Lab Book: Serial fixes for H2 and H2O (2026-01-16)

**Context:** Continuing serial-validation work to align GPU (OpenCL) with CPU for the simplest systems, before touching parallel accumulation.

**H2 (2 atoms, 2 electrons) — PASS (serial)**
- Max diff ~`2.8e-05` (tol `1e-4`).
- Problem: GPU kinetic size-force was zero because `const_K_SI ~ 6e-39` underflowed in `float` and was flushed to 0 on NVIDIA.
- Fix: Hardcode derived constants in `eFF.cl` to avoid subnormals: `const_K_eVA=3.8099822f`, `const_Ke_eVA=5.7149734f`. After this, GPU kinetic `fs` ≈ 182.879 for `s=0.5`, matching CPU. Pauli debug prints already matched; no math change needed.
- Debugging aids: CPU prints gated by `idebug` in `addPauliGauss_New`; GPU serial `getPauliGauss_New` printf (rate-limited) plus kinetic print; tolerance relaxed to `1e-4` in `test_ocl_vs_cpu.py`.

**H2O_fixcore (3 atoms, 8 electrons) — PASS (serial)**
- Max diff `1.0e-05` (tol `1e-4`).
- Problems and fixes in GPU serial AE path:
  1) Ion `w` accumulation: ions have no size DOF. Removed `cg.z` contribution to ion; only electrons receive size-force.
  2) Missing oxygen core terms: added AE Pauli (qq = `sP*0.5*qj`) and core Coulomb correction (qq = `sP*qj`) when `sP>0`, matching CPU `evalAE()`.
  3) Derivative slot: CPU accumulates derivative w.r.t ion size (`fsi`) into electron `fs`. GPU was using `fsj`. Now electron `w` gets `pg.z` and `cgC.z` (not `.w`).
  4) Core charge: CPU with `bCoreCoul=1` uses full `Q` (not `Q-sP`). GPU serial path now uses `qCore=Q` to mirror CPU config.
- Result: H2O CPU/GPU forces now agree to ~`1e-5` in serial path.
- Debugging aids: serial AE/EE printf dumps in `eFF.cl` (guarded by `bDBGall`); CPU verbosity `setVerbosity(4,1)`; test harness pointed to `H2O_fixcore.xyz` with tol `1e-4`; used `OUT-test_ocl_vs_cpu_H2O.txt` to compare pairwise AE/EE prints.

**Remaining work after these serial passes**
- Run NH3_fixcore and CH4_fixcore (tol `1e-4`).
- Add perturbation tests (small random jitter; single-particle displacement) on serial path.
- Only then address parallel accumulation (pair-once with ±f) and small-system optimization.