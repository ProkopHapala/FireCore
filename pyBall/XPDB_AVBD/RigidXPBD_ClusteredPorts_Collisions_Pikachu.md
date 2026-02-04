---
description: Rigid clustered PBD ports+collisions (3D) overview
---

# Rigid clustered ports + collisions (3D) — framework notes

This document explains the **clustered rigid solver** exposed in:

- `pyBall/XPDB_AVBD/XPDB_new.cl` (OpenCL kernels)
- `pyBall/XPDB_AVBD/XPDB_new.py` (host-side orchestration)
- `pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py` (CLI test harness + visualization)

The target audience is computational physicists/chemists and developers.

## 1) What is being solved

We represent a molecule as atoms with:

- position `x_i`
- inverse mass `w_i = 1/m_i` (stored as `pos[i].w` in device buffers)
- rigid orientation for **node atoms** as quaternions `q_i` (stored as `quat[i]`)

We split atoms into:

- **node atoms**: have rigid orientation DOFs and ports
- **cap atoms**: only have position DOFs (no orientation updates)

Ports define constraints that connect a node to its neighbors via **attachment points** ("ports") that rotate with the node.

### Rigid-port distance constraints

For a node `i` and neighbor `j`, we define a local attachment vector `r_local(i,k)` in the node’s body frame.

World-space port tip:

- `r_arm = R(q_i) * r_local`
- `tip = x_i + r_arm`

Constraint (one per neighbor/port):

- `C = | x_j - tip | = 0`

In the current implementation we treat this as a **zero-length spring-like constraint** with stiffness `K` (see solver weighting below). This keeps the neighbor `j` at the port tip.

### Collision constraints (cluster-local)

For any pair of non-bonded atoms `(i,j)` within a cluster (and discovered ghosts), we enforce non-penetration using a soft PBD-style correction if:

- `|x_i - x_j| < r_i + r_j`

with user-controlled collision stiffness `k_coll`.

## 2) Why clustering is required (parallelization)

The OpenCL implementation is designed around **workgroups of fixed size**:

- `GROUP_SIZE = 64`

A **cluster** is exactly one workgroup-sized contiguous block of atoms.

This enables:

- efficient local-memory caching of positions/radii
- per-cluster AABB computation
- cluster-local collision loops

In 3D we currently expose the clustered method as:

- `cluster_relax_ports`

from `XPDB_new.py`:

- `rigid_cluster_relax_ports_step(...)`

## 3) Data layout / packing invariants

### 3.1) Group packing

The clustered kernels assume:

- atoms are stored contiguously by cluster
- total number of atoms is a multiple of 64

The test harness provides this via:

- `pack_molecules_contiguous(..., group_size=64, pad_to_group=True)`

Padding atoms are marked by element `'X'`.

### 3.2) Nodes-first ordering (per cluster)

Inside each cluster:

- the first `nnode_per_group` atoms are **node atoms**
- remaining atoms are caps/padding

This invariant is required by kernels:

- node corrections are indexed by `inode = group_id*nnode_per_group + lid` where `lid < nnode_per_group`

### 3.3) Neighbor buffers

Global neighbor list:

- `neighs_global[i] : int4` (up to 4 neighbors)

Cluster-local neighbor list:

- `neighs_local[i] : int4`

where each entry is mapped to:

- `0..GROUP_SIZE-1` for local atoms
- `GROUP_SIZE..GROUP_SIZE+nghosts-1` for ghosts
- `-1` for not present

The mapping is built by `build_local_topology_rigid`.

### 3.4) Port buffers

For clustered solver ports are stored **node-only** in contiguous inode-space:

- `port_local_cl[inode, k]` (float4, xyz used)
- `Kflat_cl[inode, k]` (float)

Host helper:

- `upload_rigid_cluster_ports_from_atoms(port_local_atoms, bKs_atoms, nnode_per_group=...)`

extracts data for the first `nnode_per_group` atoms in each cluster.

### 3.5) Recoil slot mapping (`bkSlots`)

When a port constraint applies a correction to a neighbor atom `j`, that neighbor needs to “receive” it.

To avoid atomics, each neighbor correction is written into a unique slot. The slot index for atom `j` is found from `bkSlots[j,k]` which points into the packed `dpos_neigh` array.

For clustered solver we must build `bkSlots` consistent with inode-space size:

- slots are in `[0 .. ngroups*nnode_per_group*4)`

Host helper:

- `upload_rigid_bk_slots_clustered(neighs, nnode_per_group=...)`

## 4) Momentum conservation: what is conserved and how

This solver is Jacobi-like and uses symmetric pair corrections with explicit mass weighting.

### 4.1) Linear momentum in PBD context

Pure position projection does not explicitly integrate velocities, so **linear momentum is not directly advanced** as a velocity variable.

However, the correction rules are built from equal-and-opposite impulses (and mass-weighted splits) such that **center-of-mass drift is minimized**.

For collisions the pair correction uses:

- `dl = (rsum - dist) / (w_i + w_j)`
- each pair is seen twice → `dl *= 0.5`
- `Δx_i += n * dl * w_i`
- `Δx_j -= n * dl * w_j` (in effect)

For port constraints the correction includes angular compliance through an effective weight term:

- `w_ang = | r_arm × n |^2 * invI`

so that translations and rotations share the correction.

### 4.2) Angular momentum / rigid rotation update

Node orientation is updated by accumulating a small-angle rotation vector `Δθ`:

- `Δθ += (r_arm × P) * invI`

and then applying quaternion increment `dq(axis, angle)`.

Caps do not update orientation.

### 4.3) Double counting compensation for caps

In `apply_corrections_rigid_ports`:

- node atoms participate in node kernels
- caps do not

To match the 2D scheme, caps apply:

- `dx *= 2` for caps (so that total pair contribution matches the “seen twice” convention)

This is a critical implementation detail.

## 5) Algorithm (host-side orchestration)

The method `XPDB_new.rigid_cluster_relax_ports_step(...)` runs:

Outer loop (`outer_iters`):

1. **AABB per cluster** (`update_bboxes_rigid`)
2. **Ghost discovery + neighbor remap** (`build_local_topology_rigid`)

Inner loop (`inner_iters`):

3. **Collisions** (`compute_collision_cluster_rigid`) → writes `dpos_coll`
4. **Rigid ports** (`compute_ports_cluster_rigid`) → writes/accumulates
   - `dpos_node`, `drot_node`, `dpos_neigh`
5. **Apply corrections** (`apply_corrections_rigid_ports`)

Parameters:

- `k_coll`: collision stiffness (0 disables collisions)
- `relaxation`: PBD relaxation factor (0..1)
- `bbox_margin`: enlarges AABB overlap region for ghost discovery

## 6) CLI usage (3D test harness)

The test harness is:

- `python3 pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py ...`

Relevant flags:

- `--method cluster_relax_ports`
- `--copies N` and `--copy_shift S` (multi-molecule packing)
- `--coll_radius R` (if `R>0`, fixed radius; if `R<=0`, use `ELEMENT_DICT[e][index_Rvdw]`)
- `--k_coll K` collision stiffness
- `--outer_iters`, `--inner_iters`
- `--view none|2d|3d`, `--viz_every`, `--viz_force`

Example (no visualization, quick smoke test):

- `python3 pyBall/XPDB_AVBD/test_rigid_XPBD_molecules.py --method cluster_relax_ports --copies 2 --iters 10 --view none --viz_force 0 --noshow 1`

## 7) Current limitations / TODO

- The “projective” method in the CLI is distance constraints only (no ports).
- XPBD ports variants exist in `XPDB_new.py` but were removed from the CLI menu because they were not producing meaningful results in this harness.
- Padding atoms are inserted for packing; for stable collisions it is recommended to move them away or give them zero radius and/or huge mass (harness currently uses huge mass + zero radius by default).

