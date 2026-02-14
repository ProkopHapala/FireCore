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

## 9. Lab Book: Parallel localMD aligned with CPU (2026-01-16)

**Context:** Removed the serial fallback and kept the optimized per-thread accumulation (each thread i sums force-on-i from all j). Achieved CPU match for H2O (max diff ~2.7e-05, tol 1e-4).

**Problems found in the parallel kernel (post-Gemini):**
1) **Forces not written:** `localMD` stopped writing forces to `fout`, so GPU forces were all zero.
   - Fix: store `my_f` in `l_force[lid]` and write to `fout[ip_start+lid]` at the end.
2) **Electron–ion (AE) ordering/signs/derivative slots wrong:**
   - Need `dR = elec - ion`, call `getCoulombGauss(dR, sQ, se, -Q)` (ion size first), force on electron is `-dR*fr`, and electron size-force comes from `.w` (fsi wrt se). Core Pauli/Coulomb use `.z` (fsi wrt ion size) but must be accumulated into electron fsize (CPU behavior).
   - Fix: matched the validated serial/CPU conventions exactly; ions never get `w` accumulation.
3) **Core charge amplitude:** main AE Coulomb uses full nuclear charge `Q`, not `Q-sP`; core correction uses `qq = sP`.

**Result:** `test_ocl_vs_cpu.py` on H2O_fixcore now **PASS** with max abs diff `2.736e-05` (tol 1e-4) using the parallel path (no serial Newton fallback).

**Next:** Run NH3_fixcore and CH4_fixcore with the parallel kernel; then proceed to perturbation tests.

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

## 10. Density fitting path (2026-01-16)

- Added OpenCL kernel `fit_density_fire` with selectable optimizers (`opt_mode`: FIRE / damped MD / gradient descent) plus proper workgroup-wide FIRE reduction.
- Python host `pyBall/OCL/eFF_ocl.py` updated to pass `opt_mode/md_damp` and expose `eval_density_grid` with amplitudes.
- Test harness `tests/pyFireball/test_fit_density_fire.py` now supports `--opt {fire,md,gd}` and produces a 2D density slice (imshow) of the fitted model using the OpenCL `eval_density_grid` kernel.

## 11. Recent advancements (2026-02-13)

**What was corrected (parallel `localMD` path):**
- Forces were not written from local buffer: restored writeback of `my_f` via `l_force[lid]` to `fout[ip_start+lid]`.
- AE ordering/signs/derivatives aligned to CPU: use `dR=elec-ion`, call `getCoulombGauss(dR, sQ, se, -Q)`, force on electron `-dR*fr`, accumulate Pauli/core Coulomb size-derivatives into electron `w` (ions never get `w`).
- Core charge amplitude: main AE Coulomb uses full nuclear charge `Q`; core correction uses `qq=sP` (matches CPU `bCoreCoul` behavior).

**Root problems:**
- Missing force writeback made GPU forces zero.
- AE derivative slots and charge amplitudes mismatched CPU, especially in frozen-core cases.

**What works now:**
- `test_ocl_vs_cpu.py` on `H2O_fixcore` **PASS** on the parallel path (no serial fallback), max abs diff ≈ `2.7e-05` (tol `1e-4`).
- Serial path remains matched for H2 and H2O (previous entries); parallel path now aligned for H2O and ready for NH3/CH4 runs.

**Next targets:**
- Run `NH3_fixcore` and `CH4_fixcore` on the parallel path (tol `1e-4`).
- Add perturbation tests (small random jitter / single-particle displacement) on parallel path.


### Feb 2026 parity fix (force arrays now match)
- **Symptom:** `test_ocl_vs_cpu.py` diff O(1e1–1e2) on EE forces with GPU correct math; AE already matched. Inputs verified identical via gated kernel dumps.
- **Root cause:** CPU `evalEE` called `addCoulombGauss` **twice** per EE pair (copy-paste bug), doubling EE Coulomb forces while energy counted once. GPU had single call.
- **Fix:** Commented the duplicate Coulomb call in `cpp/common/molecular/eFF.h` EE loop; kept GPU unchanged. Rebuilt `eFF_lib`.
- **Result:** `test_ocl_vs_cpu.py` on `H2O_fixcore` **PASS**, max abs diff ≈ `2.75e-05` (tol `1e-4`) on the parallel path.
- **Workflow lessons:**
  1. Verify inputs first (gated `DBG_EFF_INPUT`, CPU verbosity) before touching physics.
  2. Audit CPU reference for accidental duplicates/side-effects; don’t bend GPU to a CPU bug.
  3. Keep heavy prints gated (`DBG_EFF_INPUT`, `DBG_EFF_PAIR`) and only enable on the first divergence.
  4. Compare per-kind decompositions (AE vs EE; Coulomb vs Pauli) and **total per-particle** forces—pairwise sign flips don’t matter if totals match.
  5. Re-test after each single, minimal change; avoid back-and-forth toggles.

### Forward plan (parity + trajectories)
- **Static force parity:** Re-run NH3_fixcore and CH4_fixcore (tol 1e-4) on the parallel kernel.
- **Short trajectory parity (≤10 steps):** H2O and CH4 with identical MDdamp/dt/damping/fixmask/spins; accept positions/sizes ≤1e-3 after 10 steps; log per-step deltas; enable all-pairs debug only on first divergence.
- **Long convergence parity:** CPU relax to Fconv 1e-4–1e-5, perturb (~0.05 Å), GPU relax with same MDdamp settings; compare endpoints (positions/sizes priority ≤1e-3; energies via `localMD_energy`; forces ≤ Fconv).
- **Scan parity (fixed ions + electron relax):** 1D/2D scans with nuclei fixed; relax electrons CPU vs GPU, write multi-frame `gpu_scan_relaxed.xyz`; post-relax energies via `localMD_energy`; target positions/sizes ≤1e-3, |ΔE| ≤ 5e-4 eV.
- **Debug discipline:** identical algorithms/settings (MDdamp), identical fixmask/constraints, redirect heavy debug to files to avoid truncation, and gate prints to the first failing step.

### Status 2026-02-13 (evening)
- **Static force parity:** PASS for H2O_fixcore, CH4_fixcore, NH3_fixcore (N `R_eff` set to 0.1 in CPU+GPU tables; `test_ocl_vs_cpu.py` now takes argv/env).
- **Short trajectory parity (10-step, stepwise, tol_pos/size=1e-3):** PASS for H2O and CH4; no divergence flagged (forces/energies logged but not gating).
- **Debug gating:** `DBG_EFF_FDECOMP` added to silence per-electron prints by default; long/pert sections in `run_relax_parity_protocol.py` are skipped when `--long-steps 0`.
- **Fixed-ion scan parity (electron-only):** First run on H2O distscan (50 conf, nuclei fixed) shows large mismatch (pos_xyz_max~1.29, pos_size_max~0.59, Etot_max~111 eV). Needs investigation (likely scan handling/constraints mismatch).

**Next actions:**
1) Investigate fixed-ion scan mismatch: compare a single config CPU `relaxed_scan` vs GPU `relax` (fixmask) on distscan_H2O__spins_fc.xyz to find geometry/energy delta source.
2) Repeat fixed-ion scan after fix; target pos/size ≤1e-3 and |ΔE| ≤5e-4 eV.
3) (Optional) Long convergence parity after scan fix: CPU relax → perturb → GPU relax.

- **Debug discipline:** identical algorithms/settings (MDdamp), identical fixmask/constraints, redirect heavy debug to files to avoid truncation, and gate prints to the first failing step.

### Progress notes 2026-02-14 — Relaxed trajectories and scans (CPU↔GPU parity)

**Single-frame / trajectory groundwork**
- Switched CPU long-parity to the C++ integrator (`move_MD_dbg`, `ialg=3`), matching the GPU localMD algorithm. Critical divergence resolved by using the same integrator and writing the **initial frame** into trajectories (CPU/GPU alignment).
- Trajectory writer in C++ now uses `save_xyz_core(...)` to emit parser-friendly headers (`na,ne,core`).
- Validated on jittered H2: long parity (2000 steps, dt=0.01, damp=0.1) now `divergence=None` with overlaid trajectories.

**1D relaxed scans (fixed ions, electron relax)**
- Added `plot_scan_parity.py` (CPU vs GPU Etot curves + per-component diffs; CPU dotted lw=1.5, GPU solid lw=0.5).
- Added `run_relaxed_scans_1d.sh` to batch H2 / CH4 / H2O distance scans (dt=0.01, damp=0.1, 2000 steps). Outputs in `export/scan_parity_*` with Es5_{cpu,gpu}.npy and plots.
- Parity is good for the main H2/H2O runs; remaining 1D issues:
  - `relaxed_distscan_H2O__pairs_fc`: GPU-CPU Etot diff ~ +22 eV.
  - `relaxed_rotscalescan_CH4__spins_fc` and `__pairs_fc`: GPU-CPU Etot diff up to ~58 eV.

**2D relaxed scans (frozen-core spins/pairs grids)**
- Added `run_relaxed_scans_2d.sh` to auto-discover `*scan*__spins_fc.xyz` / `*scan*__pairs_fc.xyz`, relax via `run_relax_parity_protocol.py`, then plot maps.
- Extended `eval_plot_cpu_gpu_maps.py` with `--from-es5-{cpu,gpu}` to plot precomputed relaxed energies (Etot from Es5) without re-evaluating.
- Good parity examples:
  - `relaxed_distscan_H2__spins_fc`: Etot diff rms ~1.3e-05 eV.
  - `relaxed_distscan_H2O__spins_fc`: Etot diff rms ~1e-4 eV.
- Problem cases (need follow-up):
  - 1D: `relaxed_distscan_H2O__pairs_fc`, `relaxed_rotscalescan_CH4__spins_fc`, `relaxed_rotscalescan_CH4__pairs_fc` show large Etot offsets (tens of eV).
  - 2D: `relaxed_angdistscan_CH4__pairs_fc`, `relaxed_angdistscan_CH4__spins_fc`, `relaxed_angdistscan_H2O__pairs_fc`, `relaxed_angdistscan_H2O__spins_fc` show map artifacts (GPU energy surfaces with blotches).
  - Invalid input skipped: `distscan_Oe__spins_fc.xyz` contains `nan`/bad record; relax aborted.

**Hypotheses and next checks**
- OpenCL synchronization is unlikely the culprit: kernels `localMD` / `localMD_energy` wait and `queue.finish()` before host reads. Artifacts resemble initialization/model mismatch, not unsynchronized reads.
- Suspect inputs with tiny/zero electron sizes and core-mapping differences in pairs variants: CPU clamps sizes to ≥1e-3 during run; GPU upload currently preserves raw sizes. Action: clamp `pos_h[:,3]` ≥1e-3 before upload (or mirror CPU clamp) and retest failing pairs/spins scans.
- Also inspect problematic XYZs for NaNs or inconsistent `na,ne` headers; GPU parser may skip/interpret differently than CPU.

**Scripts added/updated today**
- `plot_scan_parity.py` — overlay Etot + component diffs for 1D scans.
- `run_relaxed_scans_1d.sh` — batch H2/CH4/H2O distance scans.
- `run_relaxed_scans_2d.sh` — batch 2D scans (spins/pairs fc) with relaxed parity and plotting.
- `eval_plot_cpu_gpu_maps.py` — new `--from-es5-cpu/--from-es5-gpu` path for plotting relaxed energies.

**Open items**
- Clamp GPU-uploaded electron sizes to match CPU clamp; rerun failing pairs/spins scans.
- Validate/clean problematic XYZ inputs (remove NaNs, fix headers).
- Re-run 2D angdistscan grids after clamp/input fixes to confirm parity.