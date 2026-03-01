
Here is a detailed architectural design document and implementation plan for the **"Asymmetric Causal Band Softening & Dynamic n-section"** algorithm. 

This document is formatted as a set of direct instructions, considerations, and structural blueprints for an LLM coding assistant to implement the system into your PyOpenCL `FireCore` framework.

# System Architecture & Implementation Plan

## 1. Global Memory Map & 2D Grid Setup
The parallelization strategy maps $N_{pop}$ (number of candidate trajectories, e.g., 50) and $N_{steps}$ (points along the trajectory, e.g., 100) into a single 1D array of replicas: `N_rep = N_pop * N_steps`.

*   **OpenCL 2D Grid:** 
    *   `global_work_size = (atoms_padded_to_32, N_rep)`
    *   `local_work_size = (32, 1)`
*   **Replica Indexing:** Inside any kernel, `replica_id = get_global_id(1)`. 
    *   `pop_idx = replica_id / N_steps` (Which trajectory this belongs to)
    *   `step_idx = replica_id % N_steps` (Which time-slice this represents)
    *   *Crucial Logic:* If `step_idx == 0`, this replica is the causal anchor. It has no predecessor and must remain fixed to the initial molecular state (or only relax under local forces, with no causal band applied).

## 2. New OpenCL Kernels to Implement

### A. `init_trajectory_spline_kernel`
*   **Purpose:** Eliminates Python-side interpolation overhead. Generates both the AFM tip trajectory and the initial molecular coordinate guesses directly in GPU memory.
*   **Inputs:** Spline control points `(N_pop, num_cp, 3)`, Target/Initial molecule coordinates `(N_atoms, 3)`.
*   **Outputs:** `tip_positions (N_rep, 3)`, `atom_positions (N_rep, N_atoms, 3)`.
*   **Logic:**
    *   Calculate B-spline or NURBS interpolation for the tip using `step_idx / N_steps` as the parameter $t$.
    *   Perform linear interpolation for the molecule coordinates between `Initial` and `Target` states based on $t$. (The causal band will fix the unphysical parts of this guess later).

### B. `getTipMorse_kernel`
*   **Purpose:** Evaluates the interaction between the AFM tip and the specific "handle" atom of the molecule (e.g., the oxygen in PTCDA).
*   **Inputs:** `tip_positions`, `atom_positions`.
*   **Outputs:** Adds to the atom `forces` array. Updates a global `max_tip_force (N_rep)` array.
*   **Considerations:** 
    *   Identify the handle atom via an index passed as an argument. For all other atoms, this thread returns early.
    *   Use an `atomic_max` or a local float array to keep a running maximum of the tip force magnitude *across the entire MD loop*. This is critical for the bond-breaking penalty.

### C. `reduce_trajectory_gaps_kernel`
*   **Purpose:** Calculates the distance between adjacent replicas ($|R_i - R_{i-1}|$) entirely on the GPU to avoid downloading the massive `atom_positions` array to the CPU.
*   **Inputs:** `atom_positions`.
*   **Outputs:** `gap_sizes (N_rep)`.
*   **Logic:** 
    *   For each replica, compute the distance of its atoms to the atoms of `replica_id - 1`. 
    *   *(Physics Note for Coder):* A simple maximum per-atom distance or an RMSD across the 32 threads in the workgroup using local memory reduction is required. 
    *   If `step_idx == 0`, `gap_sizes = 0`.

## 3. Modifications to Existing Kernels

### A. `updateAtomsMMFFf4` (The Propagator & Causal Band)
*   **Purpose:** This kernel currently handles FIRE updates (leapfrog + velocity damping). It must be modified to include the **Asymmetric Causal Band Force**.
*   **New Arguments:** `float K_band`, `float L_allowed`. (These will be passed from Python, allowing $K_{band}$ to anneal to 0 over the MD loop).
*   **Logic to Add (Before Leapfrog):**
    *   If `step_idx > 0`:
        *   Read position of *this* atom from `replica_id - 1`: `R_prev = atom_positions[(replica_id - 1) * N_atoms + atom_id]`.
        *   Calculate distance $d = |\mathbf{R}_{current} - \mathbf{R}_{prev}|$.
        *   If $d > L_{allowed}$: apply strain-limiting force: $\mathbf{F}_{causal} = -K_{band} \times (d - L_{allowed}) \times \frac{\mathbf{R}_{current} - \mathbf{R}_{prev}}{d}$.
        *   Add $\mathbf{F}_{causal}$ to the total force.
    *   *(Physics Note for Coder):* Applying this tether *per-atom* rather than via center-of-mass is highly deliberate. It prevents the molecule from rotating unphysically out of the causal basin, avoids thread synchronization, and is computationally trivial.

## 4. CPU / Python Orchestration (The Boundary)

The Python side should be strictly reserved for high-level logic (Genetic Algorithm / CMA-ES) and orchestrating the GPU passes.

### Phase 1: Coarse Regularized Path (Python Logic)
1.  Optimizer generates 50 trial splines (control points).
2.  Upload control points. Call `init_trajectory_spline_kernel` with $N_{steps} = 20$ (coarse grid).
3.  **The Annealing MD Loop (Python):**
    ```python
    for md_step in range(TOTAL_MD_STEPS):
        # Calculate current K_band (linear decay to 0)
        current_K = K_max * (1.0 - (md_step / TOTAL_MD_STEPS))
        
        # 1. Internal molecular forces
        run_getMMFFf4_kernel()
        # 2. Surface & non-bonded
        run_getNonBond_GridFF_Bspline_kernel()
        # 3. Tip pulling force
        run_getTipMorse_kernel()
        # 4. Causal tethering & FIRE update
        run_updateAtomsMMFFf4_kernel(K_band=current_K, L_allow=L_allow)
    ```

### Phase 2: Gap Detection (CPU / GPU Boundary)
1.  Call `reduce_trajectory_gaps_kernel`.
2.  Download the lightweight `gap_sizes` array (size: $50 \times 20$ floats) to Python/NumPy.
3.  **NumPy Logic:** Find indices where `gap_sizes > L_allowed`. These are the unphysical "teleportation" tears (bifurcations).
4.  **Dynamic Allocation ($n$-section):** For each trajectory $P$, redistribute a budget of $N_{fine\_steps} = 100$ *only* into the parameter intervals $t \in [t_{tear-1}, t_{tear}]$. Create a new, dense array of spline parameters `t_dense`.

### Phase 3: Dense Resolution & Penalties
1.  Upload `t_dense`. Call `init_trajectory_spline_kernel` (using `t_dense` instead of uniform spacing).
2.  Run the exact same Annealing MD loop as Phase 1. (Because replicas are densely packed around the saddles, they will converge extremely fast).
3.  Download `max_tip_force` and `gap_sizes`.
4.  **Fitness Evaluation (Python/NumPy):**
    *   `Penalty_Target`: `Distance(R_final, Target_Coords)`
    *   `Penalty_Blockage`: `Sum(max(0, gap_sizes - Physical_Snap_Limit)^2)`. If a large gap remains despite the dense grid and $K_{band}=0$, the path is physically impossible (tip cannot overcome barrier).
    *   `Penalty_BondBreak`: `Sum(max(0, max_tip_force - F_limit)^2)`. Evaluated accurately exactly at the bifurcation limit.

## 5. Performance & Coding Warnings for the LLM

*   **Avoid Thread Divergence in MD:** The `if (step_idx > 0)` check inside `updateAtomsMMFFf4` is uniform across the entire workgroup (all 32 atoms of a replica share the same `step_idx`). This ensures zero thread divergence on the NVIDIA warp level.
*   **Memory Coalescing:** When reading `R_prev` in the causal band calculation, the memory access pattern is `(replica_id - 1) * N_atoms + atom_id`. Because `atom_id` maps to the fastest-moving dimension (workgroup local ID), this read will be perfectly coalesced on NVIDIA hardware. Do not change this memory layout.
*   **Recoil Forces:** Make sure the Asymmetric Causal Band force *does not* apply a recoil force to `replica_id - 1`. The entire point of the physics is that information only flows forward in time.
*   **FIRE Thermostat Sync:** Remember that each replica has its own FIRE thermostat state. A sudden snap across a barrier might cause a spike in velocity. The FIRE parameters (alpha, dt) must be carefully tuned or allowed to reset locally if the `gap_sizes` spike, to ensure the molecule slides down the new basin smoothly without numerical explosion.

---

Here’s a structured note set based on the requested documents and referenced code. (I followed the PathDiffusion discussion/design docs and the existing PyOpenCL MMFF/MD code.)

## 1) Algorithm concept: one-directional dissipative elastic-band softening + iterative refinement
- Goal: Find physically feasible AFM manipulation trajectories (tip spline) that drag a molecule over a surface, respecting hysteresis and bond limits.
- Key idea: Parallel relaxation of many replicas (pop × steps) with asymmetric causal tethering between consecutive replicas; gradually soften the tether to expose bifurcations. Then adaptively refine only where snapping occurs.
- One-directional causal tether (asymmetric spring):
  - Force on replica i+1 depends on (Ri+1 − Ri); no force back to Ri (arrow of time).
  - Use rest-length penalty: E = ½ K · (|Ri+1−Ri| − L_allowed)² with one-sided clamp (only if gap exceeds L_allowed). This prevents collapse of all replicas into one point and keeps parallel load balanced.
  - L_allowed ≈ |ΔR_tip| / sqrt(3N_atoms) (causality radius).
  - Anneal K_band (or K_asym) from large to 0 over FIRE/MD steps; when K → 0, remaining large gaps mark true bifurcations (snaps).
- Iterative refinement (n-section, not strict bisection):
  - Start coarse sampling along spline (e.g., 10–20 points).
  - After coarse pass, compute gaps |Ri+1−Ri|; intervals with gaps > threshold are “torn”.
  - Redistribute a fixed fine-point budget into torn intervals proportionally (dynamic n-section); smooth intervals may get few or no refinement points.
  - Optionally, for critical torn intervals, do a small number of n-section steps to localize the snap point and evaluate forces just before the jump (bond-breaking check).
- Fitness / penalties:
  - P_target: distance of final state to desired target.
  - P_block: gaps that remain > physical_snap_distance after refinement with K=0 (means tip cannot pull across barrier).
  - P_bond: max tip force before snap exceeding F_limit.
  - Optionally continuity penalty on coarse pass (teleportation) before refinement.
- Monte-Carlo/GA loop:
  - Optimize spline control points; reuse best-forward-relaxed paths as improved initial guesses for next generation (“self-healing”).
- Physical constraints to enforce:
  - Causality/hysteresis: only forward pull (Ri → Ri+1); no backward mixing.
  - No teleportation: large gaps imply unphysical jump → penalize.
  - Bond integrity: tip force before snap must stay below F_limit.
  - Tip stiffness vs. step: algorithm must be dt-independent by testing full interval displacement (coarse + refinement).
  - Rest-length tether prevents artificial collapse; annealing ensures realistic end states.

## 2) What to reuse from existing code
- pyBall/OCL/MMFF.py:
  - Topology build: MMFF.toMMFFsp3_loc(...) to generate bonded terms, neighs, REQs, etc.
  - Exclusions: _make_excl_1_2_3 for nonbonded kernel.
  - Node/cap handling and permutation helpers.
- pyBall/OCL/MolecularDynamics.py:
  - OpenCL program loader (relax_multi.cl), buffer allocation, work sizes.
  - Kernels available: getMMFFf4 (bonded), getNonBond_ex2, getNonBond_GridFF_Bspline, getSurfFlat, updateAtomsMMFFf4 (FIRE integrator), cleanForceMMFFf4, etc.
  - MDparams upload, pack_system, upload_all_systems, run_* kernel wrappers.
  - init_with_atoms for simple atom-only systems (useful for tiny tests).
- tests/tMMFF/test_ditetraceno_surface.py:
  - Example of building MMFF from mol2, packing into MD, setting surface params, running MD loop, downloading positions/forces, writing XYZ, sanity checks.
  - Pattern for downloading buffers (download function) and reporting force magnitudes.
- GridFF substrate:
  - Qgrid_ocl.npy at tests/tUFF/data/NaCl_1x1_L3 for nonbonded grid force; leverage existing getNonBond_GridFF_Bspline if needed.

## 3) What new to implement (ManipulationPathOpt.py)
- High-level orchestrator importing MMFF and MolecularDynamics.
- Data structures:
  - Population of splines (control points), step sampling (coarse + fine).
  - Tip trajectory arrays per replica (N_pop × N_steps × 3).
  - Buffers for atom positions per replica (N_pop × N_steps × natoms × 3) laid out as a flat replica dimension (as per design doc mapping).
- GPU-side changes (requires relax_multi.cl edits later; plan now):
  - Add asymmetric causal tether in updateAtomsMMFFf4: force += causal term if |ΔR| > L_allowed; only applied to replica_id>0; no recoil to previous replica.
  - Pass K_band (annealed) and L_allowed as kernel args.
  - Optional: track max tip force per replica (tip-handle Morse) in a small buffer.
  - Optional: gap reduction kernel to compute |Ri−Ri−1| on GPU to avoid downloading full coords.
- Python orchestration (ManipulationPathOpt.py):
  - Build MMFF for the molecule (could be from AtomicSystem, maybe small H2O or supplied mol2).
  - Initialize MD with per-replica packing (N_rep = N_pop × N_steps).
  - Generate tip spline samples (coarse/fine) and initial molecular guesses (spline interpolation between start/target).
  - Phase 1 (coarse, with annealing K_band): run fixed MD steps; download small summaries (gaps, max tip force maybe).
  - Gap detection → allocate fine sampling (n-section) intervals per trajectory.
  - Phase 2 (fine, K_band annealed to 0): rerun MD; compute fitness (P_target, P_block, P_bond).
  - Expose a function run_population(pop_ctrl_pts, ...) returning fitness array and diagnostics.
  - Keep debug prints for key variables per general-debug-guidelines.
- Test script test_ManipulationPathOpt.py:
  - Small molecule (e.g., H2O) placed above NaCl GridFF (load Qgrid_ocl.npy).
  - Simple tip path: move laterally across one saddle (one expected snap).
  - Run single-population, coarse+fine; log gaps, max tip force, final position.
  - Checks:
    - Gaps after fine pass near saddle should be small except at a single snap location.
    - Max tip force before snap < F_limit (test threshold).
    - Final position plausibly lands in the neighboring minimum (e.g., by a short sequential relaxed scan to verify the endpoint).
  - Output minimal npy/xyz for inspection (coarse/fine paths, gap arrays, max tip force).

## 4) Test/validation strategy
- Unit-like small scenario:
  - H2O (3 atoms) in vacuum first: verify causal tether behaves (no collapse, small gaps with interpolation; when K_band→0, path matches sequential relaxation).
  - H2O above NaCl GridFF: tip pulls O over a surface hollow; ensure a single snap is detected.
- Metrics to check:
  - Continuity: gaps ≤ L_allowed except at intended snap; after fine sampling and K=0, remaining large gaps → mark unfeasible.
  - Bond safety: max_tip_force just before snap < F_limit.
  - Target accuracy: final RMSD to desired site.
  - dt/K_tip independence: varying coarse dt should not change feasibility when full interval is considered (use coarse+refinement).
  - Optional regression: compare parallel causal-anneal result vs. explicit sequential forward scan on the same tip path (should match endpoints and snap location).

## 5) Module skeleton proposal (ManipulationPathOpt.py)
- Public entrypoints:
  - build_system(mol_or_path, gridff=None, tip_handle_idx=0, …) → (mmff, md, dims, buffers)
  - sample_tip_spline(ctrl_pts, n_steps) → tip_positions[N_steps,3]
  - init_replica_states(mmff, tip_positions, interp='linear') → apos_guess[N_rep, natoms, 3]
  - run_relax(md, tip_positions, apos_guess, K_band_sched, L_allowed, n_steps_md, track_max_tip=True) → apos_rel, max_tip_force (maybe gaps GPU)
  - compute_gaps(apos_rel, natoms) → gaps[N_pop, N_steps-1]
  - refine_intervals(gaps, tip_positions, budget_fine) → tip_positions_fine, mapping
  - evaluate_fitness(...)
- Keep debug prints for K_band, max_tip_force, largest gaps, and any NaN detection.

## 6) Next steps (for implementation phase)
- Extend relax_multi.cl (or a new variant) with causal tether args (K_band, L_allowed) in updateAtomsMMFFf4 and tip-force tracker.
- Add gap reduction kernel (or initial CPU version if download cost is acceptable for small tests).
- Implement ManipulationPathOpt.py orchestration and test_ManipulationPathOpt.py.
- Create minimal test assets: small mol2 for H2O, reuse NaCl Qgrid_ocl.npy.

These notes should complement doc/Topics/PathDiffusion/PathDiffusion_desing.md and align with the discussion’s final asymmetric rest-length tether plus adaptive refinement strategy.