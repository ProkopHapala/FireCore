## USER

look we have this forcefield based on rigid rotation of atom ports ("As Rigid as possible" ARAP method) combined with collision accelerated by clusters and bounding boxes in framework of positions based dynamics with jacobi solver in pyOpenCL

 we more or less debugged 2D version here
@XPBD_2D.cl @XPBD_2D.py @test_XPBD_2D.py 
in paticular methods

'pbd_relax', 'pbd_cluster_relax', 'pbd_cluster_fused', works well take inspiration, 'xpbd_md', 'xpbd_relax', are wrong avoid them.

see especially kernels
build_local_topology_2d
compute_collision_cluster
compute_corrections_2d
apply_corrections_2d

and eventually also 
compute_cluster_fused_2d
compute_collision_cluster_scattered
build_local_topology_2d_scattered

in @XPBD_2D.cl , this is our reference but be critical

Now we want to write very similar solver in 3D in
@XPDB_new.cl @XPDB_new.py 
and test it in
@test_rigid_XPBD_molecules.py 

this means we also need to write function to prepare topology (neibhgor) and place the molecules into groups and eventually permut indexes. See @XPBD_2D_utils.py and replicate or import it in @XPTB_utils.py 

notice we already have kernel solve_cluster_jacobi, but that is something different, it replace angles by fake bons, instead we want to use rigid ports, in 3D we rotate them by quarnions rather then by complex number, but the rest of mechanis should be the same. Especially using dpos_neigh and ghost_indices_flat

for how to do rotation using quaternions you may see kernels
integrate_rigid_explicit
integrate_and_project
apply_corrections
 @XPDB_new.cl 
we have many solvers there, and I'm not sure which works and which not, so be critical/skeptical

finally you have to implement and run test similar to those in @test_XPBD_2D.py with similar CLI, especailly I like 
    parser.add_argument('--method',     default='pbd_relax', choices=['force', 'pbd_relax', 'pbd_md', 'pbd_cluster_relax', 'pbd_cluster_relax_ports', 'pbd_cluster_fused', 'xpbd_md', 'xpbd_relax',])

I also like the mouse-picking in @test_XPBD_2D.py and in @XPBD_2D_utils.py 
mouse picking in 3D is more difficulat, so we should have option to choose 2D or 3D view in @test_rigid_XPBD_molecules.py , where 2D has mouse picking, notice we already use something like that in @test_TiledJacobi_molecules.py (but that is for fake angle-derived bonds instead of ports)


---

### Plan
I’ll first map the current 2D “known-good” implementation and the existing 3D AVBD code, then implement a **3D PBD Jacobi solver using rigid ports + quaternion rotation** mirroring the 2D kernels (`build_local_topology`, `compute_collision_cluster`, `compute_corrections`, `apply_corrections`). After that I’ll add the missing **topology/group preparation utilities** in `XPTB_utils.py` (reusing/porting from `XPBD_2D_utils.py`), and finally update `test_rigid_XPBD_molecules.py` to provide a CLI similar to the 2D test, including a 2D view mode with mouse picking.

I’m going to start by locating/reading the relevant files in the *worktree* and identifying which kernels/paths are already correct vs. obsolete.

### Findings so far (what’s already in place)
- **2D reference (good)**: `compute_corrections_2d` + `apply_corrections_2d` implement the rigid-port Jacobi step with:
  - per-port world lever arm (`r_arm`)
  - symmetric recoil into `dpos_neigh[ i*4+k ]`
  - the crucial **`*0.5`** factor because each bond is processed twice (i→j and j→i)
  - special handling for caps (`i>=nnode`) in apply (`corr*=2`) because caps don’t get double-counted.
- **3D already has a similar “ports PBD” path**: in [XPDB_new.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0) there is `compute_corrections` + `apply_corrections` and Python wrapper `rigid_ports_pbd_step()`. This is very close to what we want for the *port constraints* side (and already uses quaternion rotation).
- **What’s missing for your request** is the **cluster collision accelerated path** in 3D analogous to the 2D `pbd_cluster_relax` / `pbd_cluster_fused`, plus the **topology/group/permute utilities** for that path, plus updating `test_rigid_XPBD_molecules.py` with the 2D-like CLI and 2D/3D view modes.

### Next step
I’m going to inspect the current 3D collision cluster infrastructure in [XPDB_new.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0) + `XPDB_new.py` (AABB update, ghost discovery, collision solve), and identify what we can reuse vs. what we must implement to match the 2D “known-good” cluster pipeline (especially `ghost_indices_flat` + `dpos_neigh` usage).

## What I found in 3D (`XPDB_new.cl/.py`) and what it means
You already have **two largely separate solvers** living in [XPDB_new.cl](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.cl:0:0-0:0):

1. **Cluster collision + “fake bonds/angles” Jacobi solver (XPDB-style)**
   - Kernels: [update_bboxes](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.py:839:4-849:9), [build_local_topology](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.py:851:4-866:9), [solve_cluster_jacobi](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.py:868:4-879:13), [solve_cluster_jacobi_step](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.py:914:4-924:9)
   - Data model: fixed-size `bond_indices_*` arrays with `N_MAX_BONDED=16`
   - Ghost infra: `ghost_indices_flat` + `ghost_counts` already implemented (contiguous groups of size `GROUP_SIZE`).
   - This is **not** the rigid-port ARAP solver you want, but the **collision clustering + ghost discovery mechanism is directly reusable**.

2. **Rigid-port quaternion solver (ARAP-like)**
   - Kernels: `compute_corrections`, `apply_corrections`, plus explicit-force helpers like `integrate_rigid_explicit`.
   - Python wrapper: [rigid_ports_pbd_step()](cci:1://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPDB_new.py:593:4-628:13) is already the *correct shape* for the port constraint Jacobi step:
     - clears `delta_neigh`, `pos_delta`, `omega_delta`
     - runs `compute_corrections(nnode, ...)`
     - runs `apply_corrections(natoms, ...)`
   - This currently runs **without any cluster collision acceleration**.

So the missing piece is basically: **port-Jacobi (rigid solver) + cluster collisions (ghost-based)**, combined the same way your 2D `pbd_cluster_relax` / `pbd_cluster_fused` do.

## Critical mismatch with 2D reference (needs an explicit decision)
In 2D, the combined path works by:
- collision solver writes into `dpos_node[natoms]`
- port solver writes node + recoil into `dpos_node/dpos_neigh`
- apply gathers recoils via `bkSlots` and updates both nodes and caps, with the special cap scaling

In 3D rigid code, your current `apply_corrections` only gathers:
- `dpos_node` (nodes only)
- `dpos_neigh` via `bkSlots` (all atoms)
and applies rotation via `drot_node`.

**But there is no place yet where “collision corrections” live** for rigid mode. We need to add that in one of two ways:

- **Option A (mirror 2D split-path exactly)**:
  - add a rigid collision kernel that writes **per-atom** `dpos_coll[natoms]` (float4)
  - then call existing rigid `compute_corrections(..., accumulate=1)`-style (we’d add an `accumulate` flag) so ports add on top
  - then apply once
- **Option B (fused step)**:
  - implement a fused cluster kernel like `compute_cluster_fused_2d` but in 3D with quaternion lever arms
  - this is faster but more complex/higher-risk.

Given your “be critical/skeptical” request, I strongly suggest **Option A first** (split collision + ports), because it isolates bugs.

## Open question I need you to answer before I implement grouping/permutation utilities
You said: *“place the molecules into groups and eventually permute indexes … replicate/import `XPBD_2D_utils.py` into [XPTB_utils.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPTB_utils.py:0:0-0:0)”*.

In 3D cluster collision code (already implemented), groups are currently **contiguous blocks**:
- group `g` contains atoms `[g*GROUP_SIZE .. g*GROUP_SIZE+GROUP_SIZE)`

That works well only if you **permute** atoms so each cluster corresponds to (roughly) a molecule or spatial chunk.

### Please choose one (so I implement the right utility):
1. **Molecule-contiguous grouping (recommended for your test harness)**  
   - permute atoms so each molecule is packed contiguously
   - then clusters are just contiguous slices
2. **Scattered grouping (like 2D `*_scattered` kernels)**  
   - keep a `group_indices_flat` list and build ghost lists based on that
   - requires adding `update_bboxes_scattered_3d` + `build_local_topology_scattered_3d` equivalents in 3D

If you say “(1)”, implementation is smaller: we reuse existing `update_bboxes/build_local_topology` *as-is*.

## Status
- **Done**: extracted the 2D “known-good” algorithmic pattern and verified the 3D code already contains usable parts (quat port solve + ghost-based cluster infra).
- **In progress**: confirming exact reuse points and deciding on grouping/permutation path (your input needed).

## Next steps (once you answer the grouping question)
- **Implement split collision + rigid-port Jacobi in 3D**:
  - add rigid-mode collision kernel that uses `ghost_indices_flat` (same as cluster solver)
  - add accumulation plumbing so port corrections add on top of collision corrections
- **Port/replicate topology/group preparation into [XPTB_utils.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/XPTB_utils.py:0:0-0:0)**:
  - likely a `setup_from_bonds_3d(...)` similar to 2D’s `setup_from_bonds`
  - plus permutation helpers (nodes-first, molecule-packed, etc.)
- **Update `test_rigid_XPBD_molecules.py` CLI**:
  - add `--method` choices mirroring 2D (and avoid xpbd_* ones)
  - add `--view {2d,3d}`; 2D includes mouse picking (reuse from [test_TiledJacobi_molecules.py](cci:7://file:///home/prokop/.windsurf/worktrees/FireCore/FireCore-095810d1/pyBall/XPDB_AVBD/test_TiledJacobi_molecules.py:0:0-0:0) style)

## Question (blocking)
Which grouping strategy do you want for the 3D rigid collision clusters?

- **A)** permute atoms so groups are contiguous (simpler; reuse existing cluster kernels)  
- **B)** support scattered groups with `group_indices_flat` (more flexible; more kernels)

Once you pick A/B, I’ll start implementing the actual kernels + Python wiring.