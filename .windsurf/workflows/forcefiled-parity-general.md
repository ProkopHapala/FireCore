## Skill: scientific_debug_protocol
**Goal:** Achieve numerical parity between a Reference Implementation (e.g., C++/CPU) and a Target Implementation (e.g., OpenCL/CUDA/WebGPU) by following a strict scientific method.

**Rules of Engagement:**

1.  **Complexity Gradient:**
    *   Always start with the simplest possible system (e.g., isolated `H2O`, `NH2OH`, `formaldimine` ).
    *   Only move to aromatic/conjugated systems (`benzene`, `guanine`) once basic geometry is validated.
    *   molecular data are in `/cpp/common_resources/xyz` and `/cpp/common_resources/mol`

2.  **Phase 1: Input Topology & Parameter Check (The "Static" Check)**
    *   *axiom:* "Garbage In, Garbage Out."
    *   Before running *any* physics, compare and ensure the input data structures in memory (layout, sorting, indexing, values).
    *   **Topology:** Do both implementations have the same neighbor lists? Same bond connectivity?
    *   **Parameters:** Check stiffness ($k$), equilibrium lengths/angles/dihedrals ($l_0$), and partial charges ($Q$).
    *   **Indexing:** Watch for 0-based vs 1-based indexing errors and struct padding (float3 vs float4), sorting (node/capping atoms first?, workgroups, cluster).

3.  **Phase 2: Divide and Conquer (Component Isolation)**
    *   Compare total forces on atoms as the main criterum of match between reference and tested impementation to have overview it it works, but do not expect match from the start.
    *   Do not debug "The Forcefield" as a whole at once. Debug piece-by-piece one component at a time (bonds, angles, dihedrals, non-covalent etc.), use gating and switch on/off individual components"
    *   Disable **ALL** other forces. Run *only* Bonds. Then *only* Angles. Then *only* Non-Bonded (Van der Waals/Coulomb).
    *   *Crucial:* Verify exclusion masks (e.g., 1-2, 1-3, 1-4 exclusions). If one implementation excludes 1-3 interactions and the other doesn't, they will never match.

4.  **Phase 3: Synchronized Tracing (The "Dynamic" Check)**
    *   If Inputs match but Outputs differ, we must trace the execution line-by-line using debug prints in the code.
    *   **Instrumentation:** Inject `printf` (or `console.log`) in the inner-most loop of both implementations.
    *   **Strict Formatting:** The output string format must be *identical* (e.g., `"%d %f %f"`) to allow for automated `diff` checks.
    *   **Gating:** Gate prints by Atom ID and Step Number to avoid stdout flooding (e.g., `if(id==0 && step==0)`).

5.  **Phase 4: Single-Step Output Verification**
    *   Run exactly **ONE** iteration (integration step or force evaluation).
    *   Compare Forces ($F$) first, then Energies ($E$).

6.  **Switch Semantics (FireCore UFF/MMFF)**
    *   C++ switch APIs use **0 = keep**, **>0 = force-on**, **<0 = force-off**. Passing 0 does nothing. For UFF, use `setSwitchesUFF_NB(NonBonded=-1, NonBondNeighs=-1, SubtractAngleNonBond=-1)` to truly disable nonbonded/subtractions during bonded parity.
    *   Ensure bonds stay on when angles/dihedrals/inversions are tested (they build neighbor/hneigh tables).

7.  **Buffer-first discipline**
    *   Always compare buffers before forces: topology (neighs, a2f offsets/counts/indices), parameters (bon/ang/dih/inv params) with correct packing (e.g., UFF angParams CPU layout `[K,c0..c3]` vs kernel split).
    *   Use deterministic downloads or debug prints gated by atom/system to validate mapping.

8.  **Tolerance & precision expectations**
    *   GPU float32 vs CPU float64: expect few-e-6 for bonds/angles; up to ~5e-5 for dihedrals/inversions/trig-heavy parts. Set tolerances accordingly; don’t chase post-parity noise.

9.  **Run context and builds**
    *   Run from the intended test directory so relative kernel/data paths resolve (e.g., `tests/tUFF`).
    *   Build in `cpp/Build` (symlink to active variant such as Build-opt/asan/dbg) so loader paths stay stable; avoid ad-hoc extra symlinks.
    *   Use provided bash scripts (e.g., `tests/tUFF/run_1.sh`) to rebuild and run C++/OpenCL; do **not** call raw `make` in random dirs. Scripts clean/rebuild `libMMFF_lib.so` as needed and then run the test.