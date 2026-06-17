# USER

I want you to explore 
/home/prokophapala/git/FireCore/pyBall/XPDB_AVBD

wit the goal of making prop[er documentation and codemap 

you should 
1) describe role and purpose of each test script and each module
2) describe role and purpose of each kernel in @RRsp3.cl 

we have there also verions of files which are "old" and "new" so we somehow need to understand the history, chronology, when they were created and modified, which are deprecated, which are experimental etc. 
What are the most important to review

we should also split the discussion of kernells into few sub problems
1) basics : definition of ports and their rotation, torque, angular velocity with quaternions
1b) this is related to how evaluate local roatation frame, we can also do it in massless-way analytically to solve optimal rotation, that needs some intertia tensor and diagonalization of 3x3 matrix or quaternion fitting, right, we should describe that.
2) dealing with recoil - on GPU we do not wan async write (atomic write) we want pure gather approach, that is why we store recoils (3rd newton law) into axuliary slots from which they are assembled later to avoid collisions
3) we have version with evaluation of forces (impulses) and version which is "position based" using projective dynamics, position based dynamics, and Jacobi momentum accelerated solver,m right? 
4) then we have non-covalent interactions, if I remember properly we have some extended "halo" which we store to local memory for each workgroupto evaluate this effiicently. this involve some reindexing
5) we have also way to exclude bonded neihbors and maybe also 2nd neighbors from non-covalent interactions (like 1-2 and 1-3 exclusion)

plese revie the files and try to write into chat some broad overview of these topics about the core

---

# CODEMAP
https://windsurf.com/codemaps/fa0ddc4f-c1eb-477c-8a09-ed60d232e522-fe86ab10a43f3d18

---


# SWE-1.6

# XPDB_AVBD CodeMap Documentation

## Directory Overview

The [RigidAtomFF](../../RigidAtomFF) directory contains implementations of **Position-Based Dynamics (PBD)** and **Extended Position-Based Dynamics (XPBD)** for molecular simulations, with focus on:
- Cluster-sorted rigid body dynamics with "ports" (ARAP-style constraints)
- Collision detection with halo/ghost atom approach
- Multiple rotation solver variants (Newton substeps, shape matching, eigen methods)
- Jacobi iterative solvers with momentum acceleration

## Directory Structure

```
pyBall/RigidAtomFF/
├── RRsp3/           # Active implementation (cluster-sorted rigid ports + collisions)
├── XPDB_new/        # Experimental implementation (force-based + position-based)
├── XPDB_legacy/     # Deprecated simple PBD implementation
├── XPBD_2D/         # 2D Position-Based Dynamics implementation
└── shared/          # Shared utilities and documentation
```

## File Chronology and Status

### **Current/Active Files (Primary Focus - RRsp3/)**
- **[RRsp3.cl](RRsp3.cl)** (1311 lines) - Consolidated cluster-sorted rigid ports + collisions kernel
- **[RRsp3.py](RRsp3.py)** (624 lines) - Python harness for RRsp3.cl

### **Experimental/Larger Versions (XPDB_new/)**
- **[XPDB_new.cl](../XPDB_new/XPDB_new.cl)** (3584 lines) - Larger experimental version with both force-based and position-based dynamics
- **[XPDB_new.py](../XPDB_new/XPDB_new.py)** (1303 lines) - Python harness for XPDB_new.cl

### **Older/Deprecated Versions (XPDB_legacy/)**
- **[XPDB.cl](../XPDB_legacy/XPDB.cl)** (1130 lines) - Older simpler PBD implementation
- **[XPDB.py](../XPDB_legacy/XPDB.py)** (734 lines) - Older Python harness

### **2D Implementation (XPBD_2D/)**
- **[XPBD_2D.cl](../XPBD_2D/XPBD_2D.cl)** - 2D PBD kernels
- **[XPBD_2D.py](../XPBD_2D/XPBD_2D.py)** - 2D PBD harness
- **[XPBD_2D_utils.py](../XPBD_2D/XPBD_2D_utils.py)** - 2D utilities

### **Shared Utilities (shared/)**
- **[XPTB_utils.py](../shared/XPTB_utils.py)** (546 lines) - Utilities for molecule packing, trajectory I/O
- **[generate_viz.py](../shared/generate_viz.py)** - Visualization generation for scan outputs

### **Documentation Files (RRsp3/)**
- [RRsp3.progress.md](RRsp3.progress.md) - Progress tracking and session summaries
- [RRsp3_momentum_design.md](RRsp3_momentum_design.md) - Momentum design documentation
- Multiple `.chat.md` files recording development discussions

### **Documentation Files (shared/)**
- [Analytic_Procrustes_doc.md](../shared/Analytic_Procrustes_doc.md) - Procrustes problem documentation
- Multiple `.chat.md` files recording historical discussions

## Test Scripts - Role and Purpose

### **RRsp3-Specific Tests (RRsp3/)**
- **[test_RRsp3_convergence.py](test_RRsp3_convergence.py)** - Convergence testing with parameter sweeps (dt, relaxation, stiffness), generates error plots and trajectories
- **[test_RRsp3_debug.py](test_RRsp3_debug.py)** - Debug testing parsing kernel debug logs for collision/exclusion behavior
- **[test_RRsp3_momentum.py](test_RRsp3_momentum.py)** - Momentum conservation testing (linear and angular invariants)
- **[test_RRsp3_smoke.py](test_RRsp3_smoke.py)** - Basic smoke tests
- **[test_RRsp3_vispy.py](test_RRsp3_vispy.py)** - Interactive GUI with Vispy for pick/drag, run/stop relaxation

### **XPDB_new Tests (XPDB_new/)**
- **[test_TiledJacobi_molecules.py](../XPDB_new/test_TiledJacobi_molecules.py)** - Tiled Jacobi solver testing on molecular systems
- **[test_rigid_XPBD_molecules.py](../XPDB_new/test_rigid_XPBD_molecules.py)** - Rigid XPBD testing on molecules (H2O, CH2NH, etc.)
- **[test_XPDB_new_dump_headless.py](../XPDB_new/test_XPDB_new_dump_headless.py)** - Headless XPDB_new testing
- **[test_XPDB_new_relax_benchmark.py](../XPDB_new/test_XPDB_new_relax_benchmark.py)** - XPDB relaxation benchmarking
- **[test_procrustes_ch4.py](../XPDB_new/test_procrustes_ch4.py)** - Procrustes rotation testing on CH4 molecule
- **[ffdebug_xpdb_rigid_viz.py](../XPDB_new/ffdebug_xpdb_rigid_viz.py)** - Debug visualization for rigid XPDB
- **[ffdebug_xpdb_tiled_viz.py](../XPDB_new/ffdebug_xpdb_tiled_viz.py)** - Debug visualization for tiled XPDB
- **[test_viz_rigid_XPBD.py](../XPDB_new/test_viz_rigid_XPBD.py)** - Rigid XPDB visualization testing

### **XPDB_legacy Tests (XPDB_legacy/)**
- **[test_XPDB.py](../XPDB_legacy/test_XPDB.py)** - Basic XPDB testing with chain systems
- **[test_Jacobi_Convergence.py](../XPDB_legacy/test_Jacobi_Convergence.py)** - Jacobi convergence analysis with momentum beta sweeps
- **[test_Tiled_Jacobi.py](../XPDB_legacy/test_Tiled_Jacobi.py)** - Tiled Jacobi testing
- **[test_TiledJacobi_new.py](../XPDB_legacy/test_TiledJacobi_new.py)** - New tiled Jacobi variant testing
- **[test_TiledJacobi_dump_headless.py](../XPDB_legacy/test_TiledJacobi_dump_headless.py)** - Headless tiled Jacobi testing
- **[test_XPDB_relax_benchmark.py](../XPDB_legacy/test_XPDB_relax_benchmark.py)** - XPDB relaxation benchmarking

### **XPBD_2D Tests (XPBD_2D/)**
- **[test_XPBD_2D.py](../XPBD_2D/test_XPBD_2D.py)** - Main 2D PBD testing
- **[test_viz_XPBD_2D.py](../XPBD_2D/test_viz_XPBD_2D.py)** - 2D visualization testing
- **[ffdebug_xpbd2d_viz.py](../XPBD_2D/ffdebug_xpbd2d_viz.py)** - Debug visualization for 2D XPBD

### **Shared Utilities (shared/)**
- **[generate_viz.py](../shared/generate_viz.py)** - Visualization generation for scan outputs (implementation-agnostic)

## Modules - Role and Purpose

### **[RRsp3.py](RRsp3.py) - Primary Harness**
Functions:
- `build_neighs_bk_from_bonds()` - Build neighbor lists and back-slots from bonds
- `make_bk_slots_clustered()` - Create back-slot mapping for clustered layout
- `make_exclusions_1st_2nd()` - Build 1st and 2nd neighbor exclusion lists
- `RRsp3` class - Main OpenCL context manager:
  - Buffer allocation (pos, quat, radius, neighs, exclusions, bboxes, ghosts)
  - State upload/download (pos, quat, fixmask, radius)
  - Topology upload (neighs, exclusions, ports, stiffness, bkSlots)
  - Solver execution (step_cluster, compute_cluster_deltas, apply_cluster_corrections)
  - Momentum management (reset_momentum)
  - Diagnostics (download bboxes, ghosts, deltas)

### **[XPDB_new.py](../XPDB_new/XPDB_new.py) - Experimental Harness**
Functions:
- Same topology helpers as RRsp3.py
- `XPDB_new` class - Larger experimental harness:
  - Both force-based and position-based dynamics paths
  - Rigid body buffers (pos, vel, quat, omega)
  - Port-based rigid solver
  - Procrustes/shape matching methods
  - Clustered and global node buffers
  - Projective dynamics and force-based step methods

### **[XPDB.py](../XPDB_legacy/XPDB.py) - Older Harness**
Functions:
- `XPDB` class - Simpler PBD implementation:
  - Basic position-based dynamics
  - Bond constraints with CSR format
  - Collision neighbor lists
  - Jacobi solver with Chebyshev acceleration
  - Gauss-Seidel variants (with coloring and block approaches)

### **[XPTB_utils.py](../shared/XPTB_utils.py) - Utilities**
Functions:
- `invert_permutation()` - Permutation inversion
- `apply_permutation_to_bonds()` - Remap bonds after permutation
- `pack_molecules_contiguous()` - Pack molecules into cluster-sorted layout
- `load_xyz()` - Load XYZ files
- `masses_from_elems()` - Get atomic masses from elements
- `perturb_state()` - Apply random perturbations
- `write_xyz_with_ports()` - Write XYZ with port information
- `write_pdb_trajectory()` - Write PDB trajectories
- `plot_state_with_ports()` - Visualize state with ports

## RRsp3.cl Kernels - Detailed Breakdown

### **1. Basics: Ports, Rotation, Torque, Angular Velocity with Quaternions**

**Helper Functions (lines 51-176):**
- `quat_rotate(q, v)` - Rotate vector v by quaternion q
- `quat_from_axis_angle(axis, angle)` - Create quaternion from axis-angle
- `quat_mul(a, b)` - Quaternion multiplication
- `quat_conj(q)` - Quaternion conjugate
- `quat_delta_rotvec(q_new, q_old)` - Convert quaternion difference to rotation vector
- `quat_normalize(q)` - Normalize quaternion
- `apply_delta_rot(q, dtheta)` - Apply rotation vector to quaternion

**Physical Model:**
- Each atom has position `pos` (float4: x,y,z,inv_mass) and orientation `quat` (float4)
- Ports are local vectors `r_local` attached to atoms, rotated by atom's quaternion
- Port tip position: `tip = pos + quat_rotate(quat, r_local)`
- Torque computed as cross product: `tau = r_arm × force`
- Angular velocity related to rotation vector: `omega ≈ dtheta / dt`

### **1b. Local Rotation Frame Evaluation (Massless Analytic Solution)**

**Three rotation solver variants:**

**`compute_ports_cluster_rigid_shapematch` (lines 894-1058):**
- Builds weighted covariance matrix: `A = Σ w * (x_j - x_i) * r_local^T`
- Performs polar decomposition via Newton-Schulz iteration (4 iterations)
- Extracts optimal rotation matrix, converts to quaternion
- Then standard XPBD linear pass with that optimal rotation
- **Pros:** Near-exact best-fit rotation, very stable
- **Cons:** Heaviest register/ALU footprint, massless rotation (no angular momentum conservation)

**`compute_ports_cluster_rigid_eigen` (lines 1060-1214):**
- Centers targets by weighted centroid
- Accumulates covariance matrix (6 unique terms)
- Applies implicit K-matrix via power iteration in quaternion space (4 iterations)
- Directly returns unit quaternion (no mat→quat conversion)
- Then standard XPBD linear pass
- **Pros:** Lighter than shapematch, direct quaternion eigensolve, stable
- **Cons:** Still massless rotation, only own inertia in linear denom

**`compute_ports_cluster_rigid_substep_optimized` (lines 726-892):**
- Two-phase approach:
  1. Rotation-only Newton substeps in registers (ignoring linear masses/recoil)
  2. Full XPBD linear pass using updated quaternion
- Uses Hessian accumulation and 3x3 solve for rotation update
- **Pros:** Cheap, register-only inner loop, aggressive rotational tightening
- **Cons:** Angular momentum not conserved, quality depends on substep count

**Mathematical Background (from [Analytic_Procrustes_doc.md](../shared/Analytic_Procrustes_doc.md)):**
- 2D case: Optimal rotation angle θ = atan2(B, A) where:
  - A = Σ k_j (n_j · p_j) (weighted dot products)
  - B = Σ k_j (p_j × n_j) (weighted cross products = torque)
- 3D case: Requires eigenvalue problem on 4×4 K-matrix or SVD/polar decomposition
- The "massless" approach solves for optimal geometry without considering inertia

### **2. Recoil Handling (Gather Approach, Auxiliary Slots)**

**Problem:** In GPU Jacobi scheme, both thread i and thread j process the same bond. To avoid atomic writes, we use a gather approach.

**Solution in `compute_ports_cluster_rigid` (lines 472-605):**
- Each thread writes its recoil to `dpos_neigh[slot]` where `slot = inode * 4 + k`
- `bkSlots` array maps each atom to which recoil slots it should gather from
- In `apply_corrections_rigid_ports` (lines 1216-1311):
  - Each atom gathers recoils: `dx_port += dpos_neigh[slot]` for each slot in bkSlots[i]
  - This avoids atomic operations and maintains momentum conservation

**Double-Counting Logic:**
- When both i and j are nodes, both threads process the bond
- Linear impulse halved: `impulse_lin *= 0.5` for node-node bonds
- Angular impulse kept full (torque is local, not shared)
- Linear reconstruction: i gets 0.5*P (direct) + 0.5*P (recoil from j) = P

### **3. Force-Based vs Position-Based Dynamics**

**Position-Based (XPBD) - Current Implementation:**
- Uses constraint projection instead of forces
- Impulse magnitude: `impulse = dist / w_total`
- Weight denominator: `w_total = w_i + w_j + w_ang + alpha`
  - `w_i = inv_mass_i`, `w_j = inv_mass_j`
  - `w_ang = dot(cross(r_arm, n), cross(r_arm, n)) * inv_inertia`
  - `alpha = 1 / (K * dt^2)` (compliance term)
- Position update: `dx = n * impulse * w_i`
- Rotation update: `dtheta = cross(r_arm, n * impulse) * inv_inertia`

**Force-Based (in XPDB_new.cl, not in RRsp3.cl):**
- Explicit force/impulse integration
- Verlet-style integration step
- Then constraint projection
- Used in `rigid_force_explicit_step` in XPDB_new.py

**Jacobi Momentum Acceleration:**
- Heavy-ball method: `move = dx * relaxation + d_mom * beta`
- Momentum buffer: `dpos_mom`, `dquat_mom`
- Momentum updated after each step: `d_mom = move`

### **4. Non-Covalent Interactions with Halo Approach**

**`build_local_topology_rigid` (lines 272-376):**
- **Broad phase:** Check bounding box overlap between clusters
- **Narrow phase:** For overlapping clusters, find atoms within margin
- **Ghost atoms:** External atoms within margin are stored in `ghost_indices_flat`
- **Local indexing:** Global indices remapped to local (0-63 for cluster, 64+ for ghosts)
- **Exclusion mapping:** Global exclusion lists remapped to local indices

**`compute_collision_cluster_rigid` (lines 378-470):**
- Loads cluster + ghost positions into local memory (LDS)
- Iterates over all local atoms (cluster + ghosts)
- Checks collision: `if (dist < r_i + r_j)` apply repulsive impulse
- Uses local exclusion lists to skip bonded neighbors
- Writes collision deltas to `dpos_coll`

**Reindexing:**
- `map_global_to_local()` function (lines 261-270) maps global atom index to local index
- If atom is in same cluster: return local index (0-63)
- If atom is ghost: search ghost list, return 64+ghost_index
- This enables pure local memory access in collision/port kernels

### **5. Bonded Neighbor Exclusion (1-2, 1-3)**

**`make_exclusions_1st_2nd()` in RRsp3.py (lines 58-99):**
- **excl1:** Up to 4 first neighbors (direct bonds from neighs)
- **excl2:** Up to 4 second neighbors (neighbors of neighbors, excluding self and first neighbors)

**Usage in collision kernel:**
- `excluded8()` helper (lines 127-131) checks if j is in exclusion lists
- Collision kernel loads `excl1_local` and `excl2_local` (remapped to local indices)
- Skips collision check if neighbor is in exclusion list
- This prevents bonded atoms from colliding with each other

**Exclusion mapping in topology kernel:**
- Global exclusion lists remapped to local indices along with neighbor lists
- Enables fast register-based exclusion check in collision loop

## Most Important Files to Review

### **For Core Algorithm Understanding:**
1. **[RRsp3.cl](RRsp3.cl)** - Primary kernel implementation (all variants)
2. **[RRsp3.py](RRsp3.py)** - Primary harness and data flow
3. **[RRsp3.progress.md](RRsp3.progress.md)** - Development history and design decisions

### **For Mathematical Background:**
1. **[Analytic_Procrustes_doc.md](../shared/Analytic_Procrustes_doc.md)** - Rotation optimization theory
2. **[RRsp3_momentum_design.md](RRsp3_momentum_design.md)** - Momentum conservation design

### **For Testing and Validation:**
1. **[test_RRsp3_convergence.py](test_RRsp3_convergence.py)** - Convergence testing framework
2. **[test_RRsp3_momentum.py](test_RRsp3_momentum.py)** - Momentum conservation validation
3. **[test_RRsp3_debug.py](test_RRsp3_debug.py)** - Kernel behavior verification

### **For Experimental Features:**
1. **[XPDB_new.cl](../XPDB_new/XPDB_new.cl)** - Force-based dynamics and additional variants
2. **[XPDB_new.py](../XPDB_new/XPDB_new.py)** - Experimental harness


## Key Design Decisions (from progress logs)

1. **Cluster-Sorted Layout:** Atoms packed into groups of 64, with nodes first, then caps, then padding
2. **Local Memory Strategy:** Ghost atoms loaded into LDS for collision/port evaluation
3. **Recoil Buffers:** Auxiliary slots avoid atomic operations for momentum conservation
4. **Multiple Rotation Solvers:** Trade-offs between accuracy, speed, and momentum conservation
5. **Jacobi with Momentum:** Heavy-ball acceleration for faster convergence
6. **Exclusion Mapping:** Global→local index remapping for fast exclusion checks

This documentation provides a comprehensive overview of the RigidAtomFF system, its components, and the mathematical/algorithmic foundations of the rigid body dynamics implementation.