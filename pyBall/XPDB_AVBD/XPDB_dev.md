## USER

OK, this is good, so I pushed it to git. now we need to do something completely different perhaps in new stest script. We have new kernels in the end of @XPDB.cl and also extensive dicussion at the end of @XPDB_AVBD.md about two topics / posible optimizations
1) Gauss-Seidel with coloring using local memory for efficient implementation. 
2) modification to use boxes and lists of particles inisde the box, and particles outside the box which interact with the box (which are closer to the box than margin 2*Rmax)

please read both @XPDB.cl and @XPDB_AVBD.md and sumarize usefull insights which hel us for implementation of these two features,
Doing so list all kernels and also their interface whcih we already added in @XPDB.py , what is what, how they differ
Also mention all the new datastructures (GPU buffers) which which need to be implemented
for the graph coloring - do we have some functionaity, perhaps in python to do that?
I realized that Gauss-Seidel for collision would be difficulet because we cannot precompute the colors if we do not know the neighborlist in the first hand, right?
We would still need to compute the neighborlist, to do graph coloring but then not using it.

What I forgot actually that the box-based tiled desing (ommiting neighbor list for colisions) can be used also for jacobi not only for Gauss-Seidle so we should also do it

please make comprehesive discussion of these issues

---

Here’s a consolidated summary of the two optimization tracks plus what’s already in code:

## 1) Gauss–Seidel with coloring (local-memory)
**Kernels in XPDB.cl**
- `solve_gauss_seidel_local` (under “GS_Solver.cl” section)
  - Inputs: `curr_pos (RW)`, `pred_pos (R)`, `params (R)`, bond CSR (`bond_start, bond_count, bond_neighbors, bond_lengths, bond_stiffness`), `atom_colors`, `num_colors`, `inner_iterations`, `dt`, `omega`.
  - Uses `WORKGROUP_SIZE` threads, loads a chunk of positions into `__local float4 l_pos`, loops over colors sequentially, updates only active-color atoms, writes back to local per color, with barriers between colors. Neighbor fetch via `get_neighbor_pos` chooses local vs global.
  - Needs atom coloring precomputed and bond CSR in GPU buffers.

- `solve_block_gs` (under “BlockGS.cl” section)
  - Config: `WG_SIZE`, `NUM_COLORS`, `BLOCK_CAPACITY = WG_SIZE * NUM_COLORS`.
  - Inputs: `curr_pos (RW)`, `pred_pos (R)`, `params (R)`, bond CSR (`bond_start, bond_count, bond_neighbors, bond_lengths, bond_stiffness`), `num_atoms`, `inner_iters`, `dt`, `omega`.
  - Workgroup loads a block of `BLOCK_CAPACITY` atoms into `__local l_pos`, iterates colors inside the block, uses `get_pos_cached` (local vs global) for neighbors, writes back cooperatively to global.
  - Also needs bond CSR; no coloring array, but fixed NUM_COLORS layout (block partition).

**Python interfaces already present in XPDB.py**
- [init_GS](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:462:4-479:73), [upload_topology_GS](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:481:4-507:85), [solve_GS](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:493:4-502:9) (for `solve_gauss_seidel_local`) – but these are the older GS bindings; they set up bond buffers and `cl_colors`.
- [init_GSblock](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:524:4-543:81), [upload_topology_GSblock](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:561:4-588:133), [solve_GS_block](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:590:4-607:9) (for `solve_block_gs`) – manage padded sizes and bond CSR; coloring is implicit by block layout (NUM_COLORS, WG_SIZE).
- New shared utilities usable by future code: [ensure_neighbors](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:311:4-313:54), [reset_iteration_buffers](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:315:4-325:69), [jacobi_iteration](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:311:4-372:33), [run_jacobi_iterations](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB.py:374:4-440:19).

**Data structures needed**
- Bond CSR: `bond_start, bond_count, bond_neighbors, bond_lengths, bond_stiffness` (already used throughout).
- Coloring: `atom_colors` buffer (for `solve_gauss_seidel_local`), or block layout (implicit) for `solve_block_gs`.
- For block GS: padded sizes to multiples of `BLOCK_CAPACITY`.

**Graph coloring**
- We don’t currently have a Python coloring utility in this file. You’d need a graph-coloring step over the bond graph (for bonds) before calling `solve_gauss_seidel_local`. For collisions, pre-coloring is hard because the neighbor list is dynamic; you’d still need the neighbor list to color, so it defeats the purpose for per-frame collisions.

## 2) Box-based tiled design (ghosts, no neighbor list)
**Kernels in XPDB.cl (ghost/tiled section)**
- `build_local_topology`:
  - Inputs: `curr_pos`, `bboxes_min`, `bboxes_max`, `global_bond_indices` (int4 per atom, -1 when absent), `margin_sq`, `bbox_margin`, `num_atoms`, `num_groups`.
  - Outputs: `ghost_indices_flat` (per-cluster list, size `MAX_GHOSTS`), `ghost_counts`, `local_bond_indices` (remapped bonds into local/ghost indices).
  - Uses `MAX_GHOSTS`, `WG_SIZE`, builds per-cluster ghost list by bbox overlap + distance check (< margin_sq), then remaps bonds to local/ghost indices.

- (Likely following kernels not fully shown but implied) would consume `ghost_indices_flat`, `ghost_counts`, `local_bond_indices` to do tiled Jacobi/GS without a full global neighbor list. Current snippet only shows builder.

**New GPU buffers implied**
- `bboxes_min`, `bboxes_max` (already exist in XPDB for Verlet update).
- `global_bond_indices` (int4 per atom) for quick remap, separate from CSR.
- `ghost_indices_flat` [clusters * MAX_GHOSTS], `ghost_counts` [clusters].
- `local_bond_indices` (int4 per atom) mapped to local/ghost indexing inside a cluster kernel.

**Use for Jacobi too**
- The tiled design (bbox + ghosts) can be used for Jacobi as well; you’d run per-cluster kernels that fetch neighbors from local cache plus ghosts instead of a global neighbor list.

**Collision + coloring concern**
- For collisions, precomputing colors is impractical because neighbor lists change each rebuild; coloring needs the neighbor graph, so for collision-GS it defeats the goal of avoiding neighbor lists. Bonds are static, so coloring is feasible there.

## Actionable insights / differences
- `solve_gauss_seidel_local` vs `solve_block_gs`: both do color-sequential updates in local memory; the former uses explicit `atom_colors`, the latter uses fixed layout (NUM_COLORS, BLOCK_CAPACITY) and a cached block.
- Ghost/tile path removes the need for a global neighbor list; instead uses bboxes + ghost lists and remapped local bonds.
- To support these in Python, additional buffers and setup routines are needed:
  - For coloring GS: a coloring routine (Python) to produce `atom_colors` for bond graph; upload buffer; launch `solve_gauss_seidel_local` with `num_colors`.
  - For block GS: padding and block-level bond CSR are already in XPDB.py (GSblock).
  - For ghost/tile Jacobi/GS: upload `global_bond_indices`, allocate ghost buffers, run `build_local_topology`, then a tiled solver kernel (not yet wired) that consumes remapped local indices.

## Next steps (suggested)
1) Add Python helpers to:
   - Compute/ upload `atom_colors` for bonds (simple greedy coloring) to use `solve_gauss_seidel_local`.
   - Upload `global_bond_indices` int4 per atom for ghost builder.
   - Allocate ghost buffers (`ghost_indices_flat`, `ghost_counts`, `local_bond_indices`) and run `build_local_topology`.
2) Add/verify a tiled Jacobi (or GS) kernel that uses the ghost lists (if not yet present) and expose a Python wrapper.
3) For Jacobi using tiles: mirror the block-GS approach but Jacobi-style updates inside the block using cached/ghost positions.
4) For collisions: accept that coloring is impractical without a neighbor list; use tiled Jacobi or stick with standard Jacobi plus neighbor list.

---

## USER

can you read what you found in the discusion about optimal sizes WORKGROUPSIZE and size of local memory / size of cluster (number of atoms in the box) to optimize usage of threads on nvidia GPU? read @XPDB_AVBD.md once more to find it

---

Here is the breakdown of the optimal layout for an NVIDIA GPU (considering the 48KB Shared Memory limit).

### 1. The Constraints
*   **Warp Size:** 32 threads. Your Workgroup Size ($M$) **must** be a multiple of 32 (e.g., 32, 64, 128, 256).
*   **Max Threads per SM:** Typically 1024 or 1536 (depending on architecture).
*   **Shared Memory (L1) per SM:** 48KB is the safe baseline.
*   **Goal:** We want enough "Active Warps" to hide memory latency. A good target is to have at least **50% to 75%** of the max threads active.

### 2. The Math (Memory per Thread)
In the Block-GS approach, one thread manages $K$ atoms (one for each color phase).
If you store `float4` (16 bytes) per atom in Local Memory:
$$ \text{Bytes Per Thread} = K \times 16 $$

Let's calculate the max threads we can fit into 48KB for different Color counts ($K$):

| Colors ($K$) | Bytes/Thread | Max Threads in 48KB | Occupancy (approx) | Verdict |
| :--- | :--- | :--- | :--- | :--- |
| **2** (Red-Black) | 32 | 1500 | 100% | Excellent |
| **4** | 64 | 750 | 75% | **The Sweet Spot** |
| **8** | 128 | 375 | 37% | Poor |
| **16** | 256 | 187 | 18% | Bad (SM mostly idle) |

**Conclusion on Colors:**
You should aim for **$K=4$**.
*   4 colors are sufficient for almost any molecular topology (chains require 2, 3D crystals require 2, messy tetrahedrons require 4).
*   8 colors kill your occupancy. If you absolutely need 8, try to optimize your graph or use "Structure of Arrays" to save space (store only X,Y,Z floats, drop the W/Padding).

### 3. The Optimal Workgroup Size ($M$)
Now that we target $K=4$ (64 bytes/thread), we can fit **750 threads** per Compute Unit.
We need to divide these 750 threads into Workgroups ($M$).

*   **If $M=32$ (1 Warp):**
    *   Shared Mem per WG: $32 \times 64 = 2$ KB.
    *   WGs per SM: $48 / 2 = 24$.
    *   *Issue:* Hardware limit on "Max Workgroups per SM" is usually 16 or 32. 32 is risky. Overhead of launching 24 tiny groups is high.
*   **If $M=64$ (2 Warps):**
    *   Shared Mem per WG: 4 KB.
    *   WGs per SM: 12.
    *   Total Threads: $12 \times 64 = 768$.
    *   *Verdict:* Very Good.
*   **If $M=128$ (4 Warps):**
    *   Shared Mem per WG: 8 KB.
    *   WGs per SM: 6.
    *   Total Threads: $6 \times 128 = 768$.
    *   *Verdict:* **Best.** 128 is the standard recommendation for NVIDIA. It balances register usage and scheduling perfectly.
*   **If $M=256$ (8 Warps):**
    *   Shared Mem per WG: 16 KB.
    *   WGs per SM: 3.
    *   Total Threads: $3 \times 256 = 768$.
    *   *Verdict:* Good, but less flexible. If a kernel uses a few too many registers, the GPU might drop from 3 WGs to 2, losing 33% performance instantly.

### 4. The Optimized Layout Strategy

**Configuration:**
*   **Workgroup Size ($M$):** 128
*   **Colors ($K$):** 4
*   **Total Atoms per Block:** $128 \times 4 = 512$
*   **Shared Memory:** ~8 KB per Block (Leaves plenty of room for L1 Cache).

---

## 2026-01 Updates (implemented & tested)

### Tiled Jacobi (ghost-based) now wired and instrumented
- `build_local_topology` + `solve_cluster_jacobi` path is active via `XPDB.init_tiled`, `upload_global_bonds`, `build_tiled_topology`, `solve_tiled_jacobi`.
- New per-atom bond buffers (float4) uploaded: `cl_global_bond_L4`, `cl_global_bond_K4` so tiled solver uses correct rest lengths/stiffness (CSR no longer misaligned).
- Debug outputs from kernel:
  - Per-atom bond force `cl_debug_force_bond`, collision force `cl_debug_force_coll`.
  - Collision neighbor dump: indices (`cl_debug_neigh_indices`, 64 slots/atom) and positions (`cl_debug_neigh_pos`).
- Test harness `test_Tiled_Jacobi.py`:
  - Options: `--debug_viz` (force vectors), `--force_scale`, `--print_mapping` (per-cluster internal/ghost maps + local bond targets), `--label_atoms`, `--stretch_tol` (post-run overstretch report), adjustable dt/omega/momentum, etc.
  - Post-run overstretched bond report prints global ids, distances vs rest, cluster id, and local target (int/ghost mapping) to flag bad constraints.
  - Mapping verifier confirms local/ghost remap; labels on animation for tracking atoms.

### Key lesson on stability
- The Jacobi RHS uses `alpha = mass/dt^2 * (pred_pos - p)` each inner iter. If `pred_pos` stays at initial positions and `dt` is small, bonds lose to that pull. Using larger `dt` or seeding `pred_pos := curr_pos` (or scaling down the alpha term) is essential when testing bonds-only or weak collisions.

### Current default params (as of tests)
- Default `dt` in tiled test raised to 1.0; `bond_len` optionally 0.4 when experimenting.
- `group_size` default 64; `max_ghosts` 128; collision defaults still available.

### What is working
- Ghost/topology build and reindex verification passes.
- Tiled Jacobi solver runs with force visualization and neighbor dumps; overstretch report highlights remaining geometry issues for tuning dt/omega/momentum/bond_k.

### Remaining follow-ups
- Add a `pred_scale` or auto copy `curr_pos -> pred_pos` before tiled solves to remove unintended anchor forces when desired.
- Optional: expose debug neighbor/force downloads in XPDB.py helpers for non-visual runs.


---

## USER

aha, I see the solve_cluster_jacobi
have this
inner_iters
right? but can this work? I can imagine that if clusters are small (the whole molecule fit into one workgroup) we can run Jacobi solver completely in clocal memory, but when we have larger molecule (which does not fit into one workgroup) we can run only single iteration because we need to synchronize (read updated positions of atoms from other workgroups) and there is no outher sycnrhonization mechanism in OpenCL between workgroups than over global memory and kernel pass, right?

We need to clearly distinguis that.

Even for small molecules running local loop would be problematic because although the bonds can be local, the collisions are non-local (molecule which fit to one workgroup may have bonds only to atoms form the same group, but it colide with atoms of different molecule updated by different group), which is main blocker from having pure local memory solver. right?

We should explain this clearly
There are few solutions
1) make purely local solver for small molecules which does just local interactions, where the collisions are not present.
2) make local solver which treat interaction with other molecules as weak "external" forces which are solved by predictor (like vdW and electronstatics) not by the linear solver (corrector) like hard-spahere colisions.
3) do loop in global memroy calling kernel every iteration and keeping momentum in global buffers. That can run for molecule of any size.

The gobal kernel should be default as it is general, the other should be more experimental atempst to optimize performance further for specialized cases. 

Am I right? Make clear layout of these architectorual alternatives  (traeoff) optimized for different situations.

---

## GPT-5.1-Codex

You’re exactly right: OpenCL provides no cross–work-group barrier, so any data that must become visible to other groups (bond constraints spanning groups, collisions, momentum history) has to flow through global memory between kernel dispatches. Here’s the layout of viable strategies and their trade-offs:

### 1. Purely local solver (single work-group fits entire system)
- **How it works**: `inner_iters > 1` loop inside kernel; positions stay in local memory.
- **Requirements**:
  - All bonded neighbors reside in the same work-group.
  - Collisions are either disabled or only involve atoms inside that work-group.
- **Pros**: Maximum performance (no global traffic mid-iteration).
- **Cons**: Only works for tiny molecules or contrived test cases; breaks as soon as collisions or bonds span multiple groups because other groups never see the updates.

### 2. Hybrid local solver with “soft” external forces
- **Idea**: Keep pure local iterations for intra-group bonds while treating inter-group interactions (collisions, vdW, electrostatics) as predictor forces applied outside the Jacobi corrector.
- **Pros**: Keeps fast local behavior for stiff constraints, offloads far-field interactions to less stringent updates.
- **Cons**: Still requires a consistent predictor/corrector split; inter-group hard-sphere collisions (what we currently model) don’t fit this approximation, so it’s only applicable if we relax collision modeling.

### 3. Global loop (current default)  
- **How it works**: Run one outer iteration per kernel launch, write `curr_pos` (and, after fixing the bug, `prev_pos`) to global memory, then relaunch for the next iteration.
- **Pros**:
  - Works for any molecule size.
  - Naturally synchronizes all work-groups via the host between iterations.
  - Allows correct heavy-ball momentum once we store `prev_pos` persistently.
  - Lets the host insert diagnostics, adjust `pred_pos`, early-stop, etc.
- **Cons**: Each iteration pays a kernel launch and global memory traffic, so it’s slower than an all-local loop on small systems.

Given those realities:
- The **global kernel-per-iteration loop** must be the baseline path; it’s the only correct approach when molecules span multiple work-groups or when collisions couple groups.
- The **local inner-loop optimization** is a specialized fast path we can keep for cases where the entire system truly fits in one group and collisions are disabled.
- Diagnostic kernels (like `bond_residuals_fixed_global`) remain separate passes so they don’t pollute the hot solver.

Next concrete step is to fix [solve_cluster_jacobi](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.py:157:4-167:9) so it writes/reads `prev_pos` in global memory. Then we can drop the helper kernel and benchmarks will stay on the main path while still supporting per-iteration control.


---

## USER

OK, so now we should make clear plan how to modify @XPDB_new.cl @XPDB_new.py do reflect thise tradeoffs specialization. We do not want to have milion of legacy kernels, and we do not want to rewrite kernels from scratch. We want to do minimal changes to existing kernels so that 
1) one is global single iteration with persistent momentum doing colisions
2) second is local multiple iterations within single workgrou, but without collisions in between workgroups (i.e. molecules does not seem each other, they are like independent systems only eventually connected by external forces in predictor). For that kernel we need to write clear warrning about the limitations  
3) third one is local multiple iterations for molecules smaller than workgroup but without collisions. We should also write warning about the limitations.

make a clear plan

---

## GPT-5.1-Codex

Here’s a minimal-change plan to satisfy the three solver variants while avoiding duplicate kernels:

1. **Global single-iteration solver (default)**
   - Rename/repurpose [solve_cluster_jacobi](cci:1://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.py:197:4-207:9) so it always performs exactly one outer iteration per launch.
   - Add a mandatory `prev_pos` buffer argument; kernel loads previous positions from global memory and writes updated ones back, so heavy-ball momentum persists across launches.
   - Host ([XPDB_new.py](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.py:0:0-0:0)) loops for `iterations` times, optionally stopping early once diagnostics fall below tolerance. Collisions remain fully supported because every iteration commits to global memory.

2. **Local multi-iteration solver without cross-cluster collisions (specialized)**
   - Keep the current `inner_iters` loop inside the kernel but gate it behind a new entry point, e.g., `solve_cluster_jacobi_local`.
   - This kernel assumes:
     - All bonded neighbors of atoms in the work-group also live in that group (already required).
     - Collisions with atoms outside the group are disabled; we either skip collision handling or only draw from local memory.
   - Add a loud warning/comment plus a Python wrapper docstring explaining that this kernel is only valid when molecules do not interact via collisions (e.g., independent clusters, or interactions handled by predictor forces).

3. **Local micro-molecule solver (entire molecule per group, no collisions)**
   - Provide another entry point, e.g., `solve_cluster_jacobi_local_small`, which is essentially Variant 2 but also assumes there are *no* ghosts at all and `num_atoms <= GROUP_SIZE`.
   - Strip collision code entirely to reduce register pressure; document explicitly that this is only for standalone small molecules (e.g., H₂O, NH₃) and breaks if ghosts/collisions are needed.

**Implementation notes**
- Keep shared math (bond loops, alpha computation) identical by factoring helper inline functions/macros inside the OpenCL file; each kernel calls the same helper to avoid divergence.
- In Python, expose three methods: `solve_cluster_jacobi_global(...)`, `solve_cluster_jacobi_cluster_local(...)`, `solve_cluster_jacobi_small_local(...)`, each performing input validation (e.g., check `num_atoms <= group_size` for the small variant, enforce `k_coll = 0` when a kernel doesn’t handle collisions).
- Diagnostics (residual kernel) stay independent; benchmarks call whichever solver variant then run the residual kernel afterward when needed.

Let me know if you want me to start applying these changes or adjust the split before modifying [XPDB_new.cl](cci:7://file:///home/prokop/git/FireCore/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0) / `.py`.




---


### Status after refactor (Jan 2026)

- `solve_cluster_jacobi_step` is now the canonical **global** kernel (one outer iteration per launch) and keeps heavy-ball momentum by loading/storing `prev_pos` in global memory.
- The original multi-iteration kernel has been repurposed as `solve_cluster_jacobi_local(...)`, which caches momentum in local memory but **syncs** the start/end state through global `prev_pos`. This variant is documented as “isolated clusters only”.
- Added `solve_cluster_jacobi_nocoll(...)`, a collision-free local variant for very small molecules (entire molecule <= work-group). Warned clearly in the Python wrapper + docstring.
- Python wrapper (`XPDB_new.py`) now exposes three explicit entry points: `solve_cluster_jacobi` (global loop), `solve_cluster_jacobi_local`, and `solve_cluster_jacobi_local_nocoll`. CLI scripts (`test_XPDB_new_relax_benchmark.py`, `test_TiledJacobi_molecules.py`) select among them via `--solver {global,local,local_nocoll}`.
- Interactive viewer improvements:
  - Dragging a node now multiplies its mass by `--pin_mass_mul` (default `1e8`) while mouse-down, giving a stiff “pin” without extra kernels. Mass is restored on mouse-up.
  - Bonds are drawn in multiple layers (primary/angle/pi/pi-align/epair/epair-pair) with distinct colors + line widths so pi/epair scaffolding is easy to inspect.
- Benchmarks run using the same default global solver, so residual plots are apples-to-apples with any future GPU/JS implementation.

No new files were introduced for this iteration; only existing modules/scripts were updated.

### Jacobi + heavy-ball math (matches C++ ProjectiveDynamics_d)

For each atom `i`, the Jacobi corrector solves the PD linear system

\[ (\alpha_i + \sum_j k_{ij})\,\mathbf{p}_i^{\text{corr}} = \alpha_i\,\mathbf{p}_i^{\text{pred}} + \sum_j k_{ij}\bigl(\mathbf{p}_j^{\text{curr}} + \mathbf{d}_{ij}\bigr) + \sum_c k^{\text{coll}}_{ic}\,\mathbf{c}_{ic} \]

where

- \( \alpha_i = m_i / \Delta t^2 \) is the inertial diagonal from the projective term \(m_i (\mathbf{p}_i - \mathbf{p}_i^{\text{pred}})/\Delta t^2\).
- \(k_{ij}\) are bond stiffnesses, \(\mathbf{d}_{ij}\) encodes the projected edge direction (rest-length correction), and collision terms are treated as virtual springs \(k^{\text{coll}}_{ic}\).

Solving the scalar denominator explicitly gives the update used in all kernels:

\[ \mathbf{p}_i^{\text{corr}} = \frac{\alpha_i\,\mathbf{p}_i^{\text{pred}} + \sum_j k_{ij}\bigl(\mathbf{p}_j^{\text{curr}} + \mathbf{d}_{ij}\bigr) + \sum_c k^{\text{coll}}_{ic}\,\mathbf{c}_{ic}}{\alpha_i + \sum_j k_{ij} + \sum_c k^{\text{coll}}_{ic}}. \]

Heavy-ball momentum is applied after the pure Jacobi correction exactly like `ProjectiveDynamics_d::updateIterativeMomentum`:

\[ \mathbf{p}_i^{n+1} = \mathbf{p}_i^{\text{corr}} + \omega\bigl(\mathbf{p}_i^{\text{corr}} - \mathbf{p}_i^{\text{prev}}\bigr), \qquad \mathbf{p}_i^{\text{prev}} \leftarrow \mathbf{p}_i^{n+1}. \]

The displacement \(\Delta\mathbf{p}_i = \mathbf{p}_i^{\text{corr}} - \mathbf{p}_i^{\text{prev}}\) is persisted per atom in `momentum`/`prev_pos` buffers so that the next Jacobi launch can reuse the exact same recurrence as the C++ reference.

### Kernel specializations (trade-offs)

1. **Global step (`solve_cluster_jacobi_step`)** – one Jacobi iteration per kernel launch; all work-groups read/write `curr_pos`, `prev_pos`, and `momentum` in global memory before the host re-launches. Supports full collisions and arbitrary molecule sizes; this is the default solver used by benchmarks.

2. **Local multi-iter with collisions isolated (`solve_cluster_jacobi_local`)** – runs `inner_iters` loops entirely inside one work-group, only synchronizing with other groups at entry/exit. Valid when each work-group owns an isolated cluster (no cross-group collisions); we document this limitation prominently.

3. **Local multi-iter, no collisions (`solve_cluster_jacobi_nocoll`)** – stripped-down version for molecules that fit entirely inside one work-group and do not collide with anything else (e.g., micro benchmarks). Drops collision code to minimize register pressure; host enforces `k_coll = 0`.

The Python wrapper guards these modes, emitting warnings if a caller requests a specialized kernel while still enabling collisions or providing molecules larger than the group size. Benchmarks keep using Mode 1 to match the C++ PD reference, while Mode 2/3 remain optional fast paths for carefully controlled experiments.


### Explicit force-based rigid solver notes (30.Jan 2026)

**Goals.** Provide a faithful force/inertia fallback path for rigid XPBD that (a) keeps linear/angular momentum symmetric via action–reaction recoils, (b) exposes damping to let oscillations decay after perturbations, and (c) offers diagnostics + visualization so we can debug port geometry and force directions on tiny systems (H₂O, CH₂NH).

> Working name: **RAsp3ff** (Rigid-atom sp³ force-field), since every rigid atom exports up to four sp³-style ports. We will keep using this label for both the position-based XPBD solver and the explicit force backup to make logs/tests consistent.

**Implementation summary.**

1. **Kernel interface cleanup**
   - `clear_rigid_forces` now zeroes `force`, `fneigh`, and the cached lever arms `pneigh`.
   - `gather_port_forces` takes per-atom **local** port vectors (`port_local`, `port_n`) and simply loops over the four neighbor slots. For each slot we rotate the stored vector by the atom quaternion, compute the spring force `K (p_j - (p_i + r_arm))`, add it to the atom and store the recoil (`-F`) plus the lever arm so the integrator can reconstruct torque.
   - `integrate_rigid_explicit` converts the cached recoils into linear force + torque (`tau += r × F`) and updates velocity/quaternion with user-provided damping.

2. **Python wrapper (`XPDB_new.py`)**
   - Added `cl_rport_local` / `cl_rport_n` buffers and `upload_rigid_ports_local()` so the test harness can feed per-atom port tables without relying on type heuristics.
   - Explicit force step now accepts `damp` argument from the caller rather than hard-coding 1.0.

3. **Test harness (`test_rigid_XPBD_molecules.py`)**
   - Added perturbation helpers so we can randomly displace positions (Å) and quaternions (rad) before running either solver.
   - Added CLI flags for damping (`--damp_force`), perturbation scales, seed, optional trajectory dump (`--dump_xyz`, writes `{base}_{mol}.xyz` with ports as pseudo-atoms), and optional Matplotlib visualization (`--viz_force`).
   - Diagnostics now print `|F|` every ~5 % of the force iterations in log scale, along with linear/angular momentum drift and bond residuals for both projective and force pipelines.

**Status / results.** With `dt_force=1e-3`, `iters_force=2000`, perturb_pos=perturb_rot=0.1, and `damp_force=0.98`:

| Molecule | `|F|` start | `|F|` end | `max |d-L0|` | `|P|` / `|L|` final |
|----------|-------------|-----------|--------------|---------------------|
| H₂O      | ~5.8e1      | ~2.2e-4   | 6e-7         | 3e-11 / 2e-11       |
| CH₂=NH   | ~9.9e1      | ~9.4e-2   | 3e-4         | 4e-10 / 9e-10       |

The damping drives oscillations down quickly, so both molecules return to their rest bond lengths even after the random perturbation. Visualization shows ports and forces aligned with expectations.

**Known issues / next steps.**

- **XPBD (projective) solver:** with undistorted initial geometry it looked stable, but once we add perturbations we see angular momentum drift—likely because we still rely on the legacy port tables (per-type) instead of the per-atom local ports now used by the force solver. Need to port the same RAsp3ff port data over to the XPBD kernels and re-test from distorted starts.
- Ports are currently built in the test script from the rest geometry. If we want the force solver to be agnostic to initial orientation without depending on damping, we should store the lever arms in each atom’s body frame (i.e., rotate `rest_vector` by `q⁻¹` before upload).
- Need to double-check that explicit-force momentum remains bounded for more complex molecules once we add collision constraints; might need per-axis damping or constraint projection.
- Visualization is Python-side only; we may want a lightweight CLI flag inside `XPDB_new.py` to dump ports/forces directly from OpenCL for regression tests.

Overall the rigid-force path now shares the same design vocabulary as the XPBD solver: simple slot mapping, symmetric recoils, explicit damping, and tooling to inspect what each port is doing.


### Status after refactor (Jan 2026)

- **Rigid-force solver**: Implemented a faithful force/inertia fallback path for rigid XPBD that keeps linear/angular momentum symmetric via action–reaction recoils, exposes damping to let oscillations decay after perturbations, and offers diagnostics + visualization so we can debug port geometry and force directions on tiny systems (H₂O, CH₂NH).
- **Kernel interface cleanup**: Cleaned up kernel interfaces for `clear_rigid_forces`, `gather_port_forces`, and `integrate_rigid_explicit` to improve readability and maintainability.
- **Python wrapper updates**: Updated the Python wrapper (`XPDB_new.py`) to accept `damp` argument from the caller and added `cl_rport_local` / `cl_rport_n` buffers and `upload_rigid_ports_local()` to feed per-atom port tables without relying on type heuristics.
- **Test harness updates**: Updated the test harness (`test_rigid_XPBD_molecules.py`) to add perturbation helpers, CLI flags for damping and perturbation scales, and optional trajectory dump and Matplotlib visualization.
- **Benchmarks**: Ran benchmarks with `dt_force=1e-3`, `iters_force=2000`, perturb_pos=perturb_rot=0.1, and `damp_force=0.98` and observed that the damping drives oscillations down quickly, returning both molecules to their rest bond lengths even after the random perturbation.

**Next steps**

- Store lever arms in each atom’s body frame to make the force solver agnostic to initial orientation without depending on damping.
- Double-check that explicit-force momentum remains bounded for more complex molecules once we add collision constraints.
- Add a lightweight CLI flag inside `XPDB_new.py` to dump ports/forces directly from OpenCL for regression tests.
