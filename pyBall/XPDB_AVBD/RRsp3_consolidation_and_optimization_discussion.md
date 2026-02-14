# USER

I have there few OpenCL forcefields, mostly I'm interested in those using "ports" and "as rigid as possible" strategy, where node atoms has stored vectros in local system where neighbor atoms should be, we have also capping atoms which are kinda like second-rate citizens which just carry on, now what is importaint we have parallelization strategy of using clusters with common bounding box, and these are processed by single workgroup. For that reason and for reason of having node atoms and capping atoms we need to distinguish what is the bst strategy of packng the atoms, performance wise, from poijt of view of flexibility but also preventing to make errors. Please review our code and list all different kernels, which goes together (e.g. one kernel produce attayd fneigh which is then consumed by other kernel), and what are all the paralleization strategies using cluster/boxes/workgroups and sorting.

I guess most efficient strategy is sort atoms within one workgroup, so we have fist node atoms, and then capping atoms within each workgroup.

Other option could be consider all atoms as node atoms.

Other strategy could be sort node atoms fist globally, then capping globally, but then each cluster need to have explicit list of atoms which are in which means one more indexing indirection, which is a bit less favoudrable performance wise, however, since we already do local-reindexing in make_local_topology kernel perhaps it is not such a bit problem. Do you actually understand local_topology kernel? It is absolutely crusial. Can you read attached files .cl with the kernels source?

---

# Gemini-3-pro

This is a sophisticated implementation of **Cluster-based / Tiled Position Based Dynamics (PBD)**. You are essentially implementing a "Ghost Fluid" or "Halo exchange" method entirely on the GPU within Local Memory (LDS/Shared Memory).

Here is the breakdown of your kernel ecosystem, the parallelization strategy, and a critical look at atom packing strategies.

### 1. Kernel Workflow & Producer-Consumer Relationships

Your code has two distinct "flavors" of solvers: **Global Explicit** (Force-based) and **Clustered XPBD** (Constraint-based). Based on your interest in "ports" and "clusters", the **Clustered XPBD** pipeline is the relevant one.

#### The Clustered XPBD Pipeline
This pipeline runs one Workgroup per Cluster.

1.  **`update_bboxes_rigid`** (or `_2d` variants)
    *   **Input:** `curr_pos`, `radius`.
    *   **Output:** `bboxes_min`, `bboxes_max`.
    *   **Role:** Reduces atom positions to find the AABB of each cluster. Used for broad-phase culling.

2.  **`build_local_topology_rigid`** (The "Brain")
    *   **Input:** `curr_pos`, `bboxes`, `neighs_global`.
    *   **Output:** `ghost_indices_flat`, `ghost_counts`, `neighs_local`.
    *   **Role:**
        1.  Checks if my cluster overlaps neighbors.
        2.  Finds specific "Ghost" atoms from neighbors that are within interaction range (`margin`).
        3.  **Crucially:** Remaps `neighs_global` (global indices) to `neighs_local` (indices into the Local Memory cache: `0..GROUP_SIZE` are internal, `GROUP_SIZE..+Ghosts` are external).

3.  **`compute_collision_cluster_rigid`**
    *   **Input:** `pos`, `neighs_local`, `ghosts`.
    *   **Output:** `dpos_coll` (Correction vector for collisions).
    *   **Role:** Loads cluster + ghosts into Local Memory. Solves particle-particle collisions. Writes a position delta.

4.  **`compute_ports_cluster_rigid`**
    *   **Input:** `pos`, `quat`, `neighs_local`, `port_local`, `stiffness`.
    *   **Output:**
        *   `dpos_node`: Linear correction for the node itself.
        *   `drot_node`: Angular correction (torque) for the node itself.
        *   `dpos_neigh`: **Recoil** correction to apply to the neighbor (stored in a slot specific to the bond).
    *   **Role:** Solves the Port/ARAP constraints. This is the Jacobi "Gather" phase.

5.  **`apply_corrections_rigid_ports`**
    *   **Input:** `dpos_node`, `drot_node`, `dpos_neigh`, `dpos_coll`, `bkSlots` (Back-slots).
    *   **Output:** `pos`, `quat` (Updated state).
    *   **Role:** Sums up all deltas:
        *   Collision delta.
        *   Node constraint delta.
        *   **Recoil deltas:** It looks at `bkSlots` to find where *neighbors* wrote corrections intended for *this* atom.
    *   Updates position and orientation.

### 2. Parallelization Strategy Analysis

Your strategy is **Spatial Clustering / Tiling**.

*   **Granularity:** 1 Workgroup = 1 Cluster (fixed size, e.g., 64 atoms).
*   **Memory Model:**
    *   **Load Phase:** Threads cooperatively load their own atoms (0..63) and then the Ghost atoms (64..64+N) into `__local` memory.
    *   **Compute Phase:** All physics (collisions and bonds) happens using `__local` indices. This avoids random access to global memory, which is the bottleneck in graph-based physics.
    *   **Write Phase:** Threads write results back to global memory.
*   **Topology:** The `neighs_local` array effectively creates a "mini-mesh" that exists only inside the Workgroup's Local Memory.

**Why `build_local_topology` is crucial:**
Without this kernel, your solver inside `compute_ports_cluster_rigid` would have to check: *"Is neighbor J inside my workgroup? If yes, read Local Mem. If no, read Global Mem."* This divergence ruins performance.
By pre-calculating `neighs_local`, the solver simply does `l_pos[neighs_local[k]]`. It doesn't care if the data came from the cluster or a ghost; it's all just data in Local Memory.

### 3. Atom Packing & Sorting Strategies

You asked about the best strategy for packing Node atoms vs. Capping atoms.

#### Option A: Sort within Workgroup (Nodes First, then Caps)
**Recommendation: Strong Yes.**

*   **Layout:** `[ Node_0, ... Node_k, Cap_0, ... Cap_m ]` inside the cluster.
*   **Why it works:**
    *   **Divergence:** In your solver kernels, you have logic like `if (i < nnode) update_quat`. If all nodes are packed at the start of the group (Local ID 0 to `nnode_per_group`), the GPU executes the rotation logic for the first warp(s) and skips it entirely for the capping warps.
    *   **Memory Coalescing:** When loading `quats`, only the first $K$ threads issue loads. The others don't.
    *   **Complexity:** You need to maintain this sort order. If an atom changes type (rare in MD) or moves clusters, you must re-sort.

#### Option B: All atoms are Nodes (Treat Caps as Nodes with Identity Inertia)
**Recommendation: Only for debugging or if ratio is low.**

*   **Pros:** Code simplicity. No branching `if (i < nnode)`.
*   **Cons:**
    *   **Memory Bloat:** Capping atoms usually don't need quaternions. You are storing/loading 4x floats (`quat`) that are useless.
    *   **Compute Waste:** You perform quaternion multiplication and normalization on atoms that shouldn't rotate.
    *   If you have 1 Node per 3 Hydrogens (Caps), you are doing 3x unnecessary rotation math.

#### Option C: Global Sorting vs. Indirection (Your "Scatter" kernel approach)
**Recommendation: Avoid Indirection (`group_indices_flat`) for the Solver.**

*   **Current Scatter approach (`_scattered` kernels in 2D file):**
    *   You keep atoms fixed in global arrays.
    *   Clusters are just lists of indices: `Cluster 1 = [Atom 5, Atom 900, Atom 2]`.
    *   **Performance penalty:** Inside the kernel, `pos[group_indices[i]]` is a **gather** operation. It breaks memory coalescing.
*   **Preferred Approach (Data Movement):**
    *   Periodically (every N steps), physically move atoms in the global `pos`/`quat` arrays so that atoms in Cluster 1 are at indices 0..63, Cluster 2 at 64..127.
    *   This makes `group_indices` implicit (`gid = group_id * 64 + lid`).
    *   This maximizes memory bandwidth.

### 4. Recommendations for Your Code

#### 1. Optimization for `build_local_topology`
You correctly identified this as critical. Currently, it does an $O(N^2)$ check between bounding boxes of clusters.
*   **Current:** Loop `other_g` from 0 to `num_groups`.
*   **Improvement:** Use a **Grid** or **BVH** for the clusters themselves. Only check clusters that are spatially adjacent.
*   **Optimization:** Since you already calculate `bboxes`, you can run a quick "Broadphase" kernel that outputs a list of `potential_interacting_clusters` for each cluster. `build_local_topology` then only iterates that list.

#### 2. Handling "Back Slots" (`bkSlots` / `bkNeighs`)
In `apply_corrections`, you sum up `dpos_neigh` from neighbors.
*   **Issue:** How do you know *which* slot of the neighbor points back to *me*?
*   **Solution:** You currently use `bkSlots`. Ensure this is updated whenever topology changes.
*   **Alternative:** In `compute_ports`, you can use **Atomics** in global memory to add to the neighbor's `pos`.
    *   *Current:* `dpos_neigh[my_slot] = recoil`. Consumer gathers.
    *   *Atomic:* `atomic_add_float4( &dpos_node[neighbor_idx], recoil )`.
    *   *Trade-off:* Atomics are fast on modern GPUs if collisions are low, but the Gather approach (your current one) is deterministic and generally safer for stable constraints. Stick with your current Gather approach.

#### 3. Specific Packing Recommendation
Go with **Sorted Coalesced Clusters**.

1.  **Global Array Layout:**
    `[ Cluster 0 Data | Cluster 1 Data | ... ]`
2.  **Inside Cluster 0 Data:**
    `[ Node 0, Node 1... Node K | Cap 0, Cap 1... Cap M | Padding... ]`
3.  **Kernel Parameters:** Pass `int nnode_per_group` to the kernels.
4.  **Solver Loop:**
    ```c
    // Load data...
    // ...
    // Solve
    if (lid < nnode_per_group) {
        // I am a Node: Solve rotation + translation
    } else if (lid < count_valid_atoms) {
        // I am a Cap: Solve translation only
        // (Constraint logic simplifies because I have no ports, just a connection to a parent)
    }
    ```

#### 4. The 2D vs 3D Discrepancy
Your 2D code (`XPBD_2D.cl`) uses a slightly different strategy (`_scattered`) with indirection arrays (`group_indices_flat`).
*   **Recommendation:** Align the 2D code to the 3D strategy (`XPDB_new.cl`). The 3D strategy (implicit indexing `grp * SIZE + lid`) is significantly faster on hardware because the memory unit can fetch contiguous blocks of 1024 bytes (position vectors) in a single transaction. The scattered approach forces the memory unit to fetch distinct cache lines for every atom.

### Summary
Your `XPDB_new.cl` logic is solid. It implements a Tiled Jacobi solver. To maximize performance:
1.  Physically sort atoms in global memory so clusters are contiguous.
2.  Inside the cluster block, put Nodes first, Caps second.
3.  Keep the `build_local_topology` kernel; it is the bridge that allows you to use fast Local Memory for irregular graph physics.

---

# USER

OK, I think you understand the algorithm very well. You also understand the issues. I think we have too many altarnatives and need to focus and consolidate now. So 1) focus only on 3D 2) foxus only on PDB 3) focus only on version with are cluster/box accelerated and include collisions

There are few issues we need to solve
1) sorting. You say that cluster-sorted layout is the best. I also think so. However, I'm not so sure if so much. Lets analyse it more deeply. The think is that once I create the local reindexing (which I must do anyway because of ghost atoms from other cells) we will read this from local memory anyway. 

the cost of l_pos[neighs_local[k]] is perhaps quite small if l_pos is in local memory and neighs_local[k] is in local memory or even in registers (the 1st neighbors for each atom (max 4 ints) we can store in registers of given thread (thread=atom i) in single int4, what you think?

2) to reduce cost of reading atoms to local memory we should prehaps merge the bonding and collision kernells into one, I mean compute_collision_cluster_rigid and compute_collision_cluster_rigid

3) there is some question of collision exclusion - bonded atoms should not colide. But maybe also secons-neighbors atoms should not colide. We should create special exculision list for that and use efficient algorithm for this exclusion, which rely on sorting. Not sure hot the algorithm reliying on sorting is compactible with our cluster-soring ? please see code attached at the end

4) what could be a bit problem is how to prepare the system - intitial remaping. E.g. imagine we have multiple molecules. we should create clusters from them. For testing at leas we should consider molecules are smaller than 64 atoms, so they fit into single cluster. Then preparation of initial sorted buffers should be simple right?

Few notes which makes it easy and we should keep it in minde:

1) the cost of local topology may be substantial (seach in other clusters other bounding boxes) but it is less performance critical because we do not run it every jacobi update. We typically 

2) the splitting of atoms into clusters is stable over the whole simulation (based on bonding topology, e.g. benzene runs or other chemical groups form cluster), we set clustes at the start on host. There is not exchange or movement of atoms between clusters (only ghost atoms are reasignemment as the system moves). It is importaint that the groups does not occupy always all 64 slots, there may be some padding. This paadin is removed in local mapping (the in-cluster atoms are followed by out-cluster=ghost atoms in the local memmory - we do not want to have voids in local memory as it is precious resource and ofthen size of  group local memory is main bottleneck)

---

# Gemini-3-pro

This is the right moment to consolidate. We will focus purely on **3D**, **XPBD**, and the **Cluster-Accelerated Fused Solver**.

Here is the definitive architecture to solve the issues of sorting, kernel merging, and exclusions efficiently.

### 1. Data Layout & Sorting Strategy
You are correct: **Cluster-Sorted** is the only way to go for high performance.

**Global Memory Layout:**
Physical arrays (`pos`, `quat`, `vel`) are sorted by Cluster ID.
```text
Index: [ 0 .. 63 ] | [ 64 .. 127 ] | [ 128 .. 191 ] ...
Data:  [ Cluster 0 ] | [ Cluster 1 ] | [ Cluster 2 ] ...
```

**Internal Cluster Layout (The "Split" Sort):**
Inside every 64-atom block, we strictly order Node atoms first, then Capping atoms.
```text
[ Node_0, Node_1 ... Node_k | Cap_0, Cap_1 ... Cap_m | (Padding) ]
```
*   **Why:** This allows the kernel to branch efficiently. `if (lid < nNodes)` executes ARAP logic. The rest execute simpler logic. Warps don't diverge as much.

**Preparation (Host Side):**
*   **Molecules < 64 atoms:** Trivial. Assign Molecule ID = Cluster ID.
*   **Molecules > 64 atoms:** Use a graph partitioner (like METIS) or a simple spatial greedy splitter to chop the molecule into chunks of 64.
*   **Padding:** You *must* pad arrays to multiples of 64. The extra atoms can be marked "inactive" (mass = 0 or infinite), so they don't move.

---

### 2. The Exclusion Strategy (Solved via Registers)
You raised a critical point: *How to exclude bonded atoms from colliding without a slow global exclusion list?*

**Solution: Register-Based Exclusion Checks.**
Since we are using "Ports" (Topology), every atom *already knows* its bonded neighbors. We typically have a max of ~4 bonds per atom.
We load these bond indices into **Private Registers** at the start of the kernel. When checking collisions in Local Memory, we compare indices against these registers.

**No global lookups. No sorted exclusion lists.**

---

### 3. The "Holy Grail" Fused Kernel
We will merge **Load**, **Topology**, **Collision**, and **Bond Constraints** into one kernel. This keeps data in L1/LDS cache.

**Producer:** `search_ghosts` (Replaces `build_local_topology`)
**Consumer:** `solve_cluster_fused`

#### A. The Producer: `search_ghosts`
*   **Input:** `bboxes` (of clusters).
*   **Output:** `ghost_indices` (list of global indices of atoms from *neighboring* clusters that are close enough to interact).
*   **Logic:** Simple AABB overlap check. If Cluster A overlaps Cluster B, add atoms from B to A's ghost list.

#### B. The Consumer: `solve_cluster_fused` (The Code)
Here is the consolidated kernel logic. This handles packing, exclusions, and solving.

```c
#define GROUP_SIZE 64
#define MAX_GHOSTS 32
#define LDS_SIZE (GROUP_SIZE + MAX_GHOSTS)

// Helper: Check if index 'j' is in the exclusion list stored in registers
inline bool is_excluded(int j, int4 bonds) {
    return (j == bonds.x) || (j == bonds.y) || (j == bonds.z) || (j == bonds.w);
}

__kernel void solve_cluster_fused(
    // 1. Data Buffers (Sorted by Cluster)
    __global float4* pos,
    __global float4* quat,
    __global float4* prev_pos, // For XPBD velocity
    
    // 2. Topology (Indices are Local: 0..63)
    __global int4*   bonds,    // Bonded neighbors (up to 4)
    __global float*  params,   // Stiffness, etc.
    
    // 3. Ghost Data (Produced by previous kernel)
    __global int*    ghost_indices, 
    __global int*    ghost_counts,
    
    // 4. Constants
    const int n_nodes_per_cluster // e.g. how many are ARAP nodes vs Caps
) {
    // --- 1. SETUP LOCAL MEMORY ---
    __local float4 l_pos[LDS_SIZE];  
    __local float4 l_quat[LDS_SIZE]; // Only needed for first GROUP_SIZE usually
    
    int gid = get_global_id(0);
    int lid = get_local_id(0); // 0..63
    int grp = get_group_id(0);

    // --- 2. COOPERATIVE LOADING ---
    // A. Load my own cluster (0..63)
    float4 my_pos  = pos[gid];
    float4 my_quat = quat[gid];
    l_pos[lid]  = my_pos;
    l_quat[lid] = my_quat;

    // B. Load Ghosts (64..96)
    // Threads cooperatively read ghost_indices array, then fetch from global pos
    int n_ghosts = ghost_counts[grp];
    if (lid < n_ghosts) {
        int ghost_global_idx = ghost_indices[grp * MAX_GHOSTS + lid];
        l_pos[GROUP_SIZE + lid] = pos[ghost_global_idx]; 
        // We usually don't need ghost quats unless ghost is also an ARAP node we constrain against
    }
    
    // C. Load Topology (Exclusions) into REGISTERS
    // Bonds are stored as LOCAL indices. 
    // e.g., if bonded to atom 5 in cluster, val is 5. 
    // If bonded to ghost, val is 64+. 
    // If no bond, val is -1.
    int4 my_bonds = bonds[gid]; 

    barrier(CLK_LOCAL_MEM_FENCE);

    // --- 3. XPBD SOLVER LOOP ---
    // We can iterate multiple sub-steps here without leaving the kernel
    for (int iter=0; iter<NUM_SUBSTEPS; iter++) {
        
        float4 dpos = (float4)(0.0f);
        float4 drot = (float4)(0.0f); // Torque/Rotation delta

        // --- PHASE A: COLLISIONS ---
        // Iterate over ALL atoms in Local Memory (Own + Ghosts)
        int total_atoms = GROUP_SIZE + n_ghosts;
        
        for (int j=0; j<total_atoms; j++) {
            if (j == lid) continue; // Don't collide with self
            
            // FAST EXCLUSION CHECK (Registers)
            if (is_excluded(j, my_bonds)) continue; 
            
            // Check Collision
            float4 other_pos = l_pos[j];
            float dist_sq = distance_squared(my_pos, other_pos);
            float rad_sum = MY_RADIUS + OTHER_RADIUS; // Simplify for example
            
            if (dist_sq < rad_sum*rad_sum) {
                // Compute collision impulse
                // Add to dpos
            }
        }

        // --- PHASE B: CONSTRAINTS (ARAP / PORTS) ---
        // Branch based on packing: Nodes vs Caps
        
        if (lid < n_nodes_per_cluster) {
            // I am a NODE atom. I have ports.
            // Solve ARAP rotational constraints using my_bonds
            // Access neighbors via l_pos[my_bonds.x], etc.
            // Calculate dpos and drot
        } 
        else {
            // I am a CAP atom.
            // I strictly follow my parent node.
            // Solve distance constraint or rigid attachment
        }

        // --- PHASE C: UPDATE & SYNC ---
        // Apply deltas
        my_pos += dpos;
        my_quat = normalize(my_quat + drot); // Simplified quaternion update

        // Update Local Memory for next iteration/neighbor visibility
        barrier(CLK_LOCAL_MEM_FENCE);
        l_pos[lid] = my_pos;
        l_quat[lid] = my_quat; // If neighbors need my orientation
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- 4. WRITE BACK ---
    pos[gid] = my_pos;
    quat[gid] = my_quat;
}
```

### 4. Why this fixes your problems

1.  **Sorting:** By putting Nodes 0..K and Caps K..63, the `if (lid < n_nodes)` creates a clean split. The first half of the warp does complex math, the second half does simple math.
2.  **Kernel Merging:** We removed the overhead of writing `dpos` to global memory and reading it back.
3.  **Exclusions:**
    *   Old way: Read `excl` array from global memory (Cache miss).
    *   New way: `my_bonds` is in **Registers**. Comparing `j == my_bonds.x` is essentially free.
    *   Second neighbors (1-3 exclusions): You can pack these into another `int4` register if needed, or simply assume that ARAP constraints are stiff enough that 1-3 collisions are impossible to trigger (often true in rigid clusters).
4.  **Local Topology:** The `bonds` array stores indices relative to the local cache.
    *   `0..63`: Inside cluster.
    *   `64+`: Ghost slot.
    *   This makes `l_pos[index]` safe and fast.

### 5. Summary of Parallelization Strategy

1.  **Broadphase (Kernel 1):** Compute Bounding Boxes of Clusters.
2.  **Ghost Search (Kernel 2):** Check Bounding Box overlaps. Fill `ghost_indices` buffer.
3.  **Physics (Kernel 3):** One workgroup per cluster.
    *   Load Cluster + Ghosts -> LDS.
    *   Load Topology -> Registers.
    *   Solve -> Write Global.

This is the most efficient mapping for an OpenCL GPU implementation of your specific forcefield requirements.

---

# USER

now critical consideration if we should use 
solve_cluster_fused
or split compute_collision_cluster_rigid and compute_collision_cluster_rigid

it is not quite obvious to me both have some pros and cons.

1) the advantage of fused apporahc is that we avoid the overhead of enquiing more kernels (if taks size is small cost of enquing can be substantial part of total budget), also big difference it that we pre-load the atoms into local memory only once tand then reuse them. It is typicall that collision and bonding neighbors has large overlap between the atoms within the cluster, so we better reuse data loaded from global memory. 

2) if we split the collision kernel and bonding kernel we can use more efficinet layout of local memory and maybe also registers. Only the bonding interactions need the explicit port geometry and rotation quaternions, and export the neighbors recoils, on the other hand bonding geometry does not need the second neighbor exclusions, does not need the ghost atoms etc. Actually this kernel can run efficiently without clusters. Perhaps optimal usage of the kernel is achieved if we simply sort all atoms node first and then take workgroup by 32-threads (nVidia wrap), But this would mean either reading negibhbors from global memory, or having separate mapping (commin neighbor) to preload them. It may be more efficient to make sparate local mapping common neighbors for bonding interaction and ghost_atoms for collisions optimized for each layout

---

# Gemini-3-pro

You have hit on the classic GPU optimization trade-off: **Occupancy (Registers) vs. Memory Bandwidth**.

Let's evaluate this based on the specific math of ARAP/XPBD and the hardware characteristics of modern GPUs (AMD RDNA / NVIDIA Ampere+).

### The Verdict: Use the FUSED Kernel
**Why?** In Forcefield/MD codes, **Memory Bandwidth is the bottleneck, not Compute.**

If you split the kernels, you are forced to:
1.  Read `pos`, `quat` (Bond Kernel) $\rightarrow$ Write `dpos_bond`.
2.  Read `pos` (Collision Kernel) $\rightarrow$ Write `dpos_coll`.
3.  Read `dpos_bond`, `dpos_coll` (Integrate Kernel) $\rightarrow$ Update `pos`.

That is **3x Global Memory Reads and 3x Writes** per time-step.
A Fused kernel does **1x Read and 1x Write**.

Furthermore, there is a synergy you missed in point #3: **Exclusions**.
*   The **Bonding** logic loads the indices of bonded neighbors (to calculate constraints).
*   The **Collision** logic needs those exact same indices (to *exclude* them).
*   By fusing, you load the topology into registers *once* and use it for both physics steps.

---

### How to make the Fused Kernel Efficient (Addressing your concerns)

You are worried about register pressure and the complexity of local memory mapping. Here is the strategy to solve that.

#### 1. The Unified Local Indexing Strategy
You shouldn't have "separate mappings" for bonding and collision. That kills performance. You need **One Local Map to Rule Them All**.

*   **Concept:**
    *   Atoms `0..63`: The Atoms in this Workgroup (Cluster).
    *   Atoms `64..N`: The Ghost Atoms (Neighbors from other clusters).

*   **The Trick:**
    *   For **Collision**: You iterate `0..N` (brute force or spatial hash in LDS).
    *   For **Bonding**: You need to know which index represents your bonded neighbor.
    *   **Pre-computation:** The `build_local_topology` kernel must convert global bond indices into *Local* indices.
        *   If Bond-Neighbor is in my cluster $\rightarrow$ Index `0..63`.
        *   If Bond-Neighbor is outside $\rightarrow$ It *must* be in the Ghost list (or physics is broken). Find its index in the Ghost list (e.g., `65`).
        *   Store this converted index in `bonds_local`.

#### 2. Register Optimization (The "State Machine" Approach)
To avoid running out of registers (SGPR/VGPR), do not try to do everything simultaneously. Use the fact that `dpos` is an accumulator.

**Kernel Pseudocode:**
```c
__kernel void solve_fused(...) {
    // 1. Load Data to Local Memory (LDS)
    // Coalesced load of Pos, Quat.
    // Load Ghosts.
    // Load Converted Topology (int4 bonds_local) into REGISTERS.
    barrier(CLK_LOCAL_MEM_FENCE);

    float4 acc_dpos = (float4)(0.0f);
    float4 acc_drot = (float4)(0.0f);

    // --- BLOCK A: BONDING (ARAP) ---
    {
        // Compiler can optimize registers here because collision vars aren't alive yet
        if (is_node) {
           // Loop over my 4 bonds
           // int neighbor_idx = bonds_local.x; // e.g., 65 (Ghost)
           // float4 p_neigh = LDS[neighbor_idx];
           // ... Compute ARAP ...
           // acc_dpos += result;
           // acc_drot += result;
        }
    }
    
    // --- BLOCK B: COLLISIONS ---
    {
        // Registers used for ARAP math (quaternions) can be released/reused here
        // We only hold 'acc_dpos' and 'pos'
        
        // Loop j from 0 to n_ghosts
        // if (j == bonds_local.x) continue; // EXCLUSION using register data!
        // ... Compute Collision ...
        // acc_dpos += result;
    }

    // 3. XPBD Update
    // Apply acc_dpos to pos
    // Write Global
}
```

### Why "Split" is actually slower for your case

You suggested splitting because Bonding doesn't need Ghosts/Clusters in the same way.
*   **False assumption:** If an atom is bonded to a neighbor, that neighbor is physically close. Therefore, that neighbor is *already* in the Ghost list required for collisions.
*   If you split, the Bonding kernel would typically need to do random access to Global Memory (`pos[bond_idx]`) because you wouldn't have the pre-loaded Cluster+Ghost cache. Random global access is 10x slower than LDS access.

### Recommendation: The "Cluster-Accelerated Fused Solver"

1.  **Host/Topology Kernel:**
    *   Sort atoms by cluster.
    *   Translate `bonds_global_id` to `bonds_local_id` (Indices into the `64+Ghost` array).
2.  **Solver Kernel (One Workgroup per Cluster):**
    *   Load Cluster + Ghosts -> LDS.
    *   **XPBD Loop (e.g., 2-5 substeps inside kernel):**
        *   Accumulate Bond Constraints (using LDS indices).
        *   Accumulate Collision Constraints (using LDS iteration + Register exclusion).
        *   Update `l_pos` (Local Memory Position).
        *   `barrier()`
    *   Write Final Result -> Global Memory.

**This enables "Sub-stepping"**:
The biggest advantage of Fused is that you can run 5 iterations of the solver *entirely in L1 cache / Shared Memory* without touching Global RAM. This makes the simulation incredibly stiff and stable (rigid-as-possible) effectively for free. You cannot do this with Split kernels.

---

# USER

I'm not so sure I can run this entirely in local cash (single kernel) this I can do only if I can do completely "gather" operations which means
1) the worgroups are independnet (no interactions with updates of atom positions from other clusters)
2) as you can see in my original kenrls I was using explicit recouls fneigh or dpos_neigh, which complicates the situation, but I found it necessary to properly ensure conservation of linear and angular momentum (maybe it can be done without it but we were struggling a lot and in the end resorted to this solution storing explicit recoils and writing them assembling them to avoid atomic writes - which we must avoid in OpenCL), please review the whole algorithm in my original implementtion.

=> I think we cannot do local-memory sub-stepping inside one kernel call, as we need global synchronization of updates between groups, and the assemble recoil (dpos_neigh) step

Now about  tradefoff - I agree that reading and writing global memory is substantial bottleneck. But having unbalanced ernels and branch-divergence is often even higher. And I just see that fused kenrel (where some atoms are node and some capping) is just to heterogenous, and this may lead to problamting blaancing and divergences. If I separate the two I can make differnt workgroups for the bonding kernel (only for node atoms) and for collisions (for all atoms), I can use the local memory more efficiently (I do not need collision ghost atoms for bonding kernel, I do not need port-direction vectors for collisions). You say that the neighbors indexes are same for both, but are not. For bonding kernel I need only the 1st neighbors, while for collision I want to exclude the second neighbors, that is substantial difference (4 1st neighs, 16 2nd neighs).

In the end I think it is not at all obvious what is faster. Also splited kernels allow easier debugging by pieces. 

For that reason I'm thinking we should try both. 

The think is however when I sort node first globally (to satisfy bonding kernel) I would need explicit indexes of atoms in groups, which may be inefficient for other puroposes (collision kernel, the local topology builder ...)

so perhaps what is best is to keep cluster-sorting, and that each cluster have initial index (offset) and number f node and capping atoms.
The bonding kernel will be still split into worgroups based on these clustrs, but the workgroups can be smaller (as they operate only on node  atoms). Dowes it make sense?

Think about this deeply and give me deteild analysis of both solutions fused, and split (bonding, collisions), after this analysis we may sketch the actuall kernels.

---

# Gemini-3-pro

You have successfully identified the critical bottleneck: **The Recoil / Global Synchronization requirement.**

You are absolutely correct. If your physics model requires an explicit `dpos_neigh` (recoil) write that must be visible to other clusters to conserve momentum, **you cannot do Local Memory sub-stepping.** The moment you need to "communicate" a force to a ghost atom (which lives in another Workgroup), you must hit Global Memory and finish the kernel.

Since we lose the main benefit of Fusing (sub-stepping in L1 cache), the **Split Kernel strategy becomes superior.** It allows us to optimize register usage and execution masks (divergence) for the distinct physics types.

Here is the detailed analysis and the concrete Kernel Design for the **Cluster-Sorted Split Approach**.

---

### Analysis: Why Split is Better for Your Case

| Feature | **Bonding Kernel (ARAP/Ports)** | **Collision Kernel (Sphere-Sphere)** | **Conflict in Fused Kernel** |
| :--- | :--- | :--- | :--- |
| **Active Atoms** | **Nodes Only** (0..K) | **All Atoms** (0..64) | **Divergence:** 50% of threads (Caps) are idle during heavy math. |
| **Math Intensity** | **High** (Quaternions, Cross products, Mat3 mult) | **Low** (Distance checks, simple deltas) | **Register Pressure:** Fused kernel needs registers for *both*, reducing occupancy. |
| **Neighbor Access** | **Specific** (Explicit IDs: `bond_idx_1, bond_idx_2`) | **Spatial** (All atoms within radius $R$) | **Data Waste:** Fused kernel loads spatial ghosts that Bonding doesn't need. |
| **Output** | `dpos_node`, `drot`, **`dpos_recoil`** (Scatter) | `dpos_coll` (Local Accumulation) | **Memory Patterns:** Bonding scatters writes; Collision is local. |
| **Exclusions** | None (calculates forces *for* them) | **Complex** (Must skip 1st & 2nd neighbors) | **Logic Complexity:** Managing exclusion lists inside a giant loop is messy. |

### The Winning Strategy: Cluster-Sorted Split

We maintain your **Cluster-Sorted Layout** in Global Memory:
`[ Node_0...Node_k | Cap_0...Cap_m ]` (padded to 64).

We split the physics into three kernels:
1.  **`solve_bonds_rigid`**: Runs on **Nodes only**. Reads Global. Writes Recoil.
2.  **`solve_collisions_rigid`**: Runs on **Clusters**. Uses LDS. Reads Ghosts. Writes Collision.
3.  **`apply_updates_rigid`**: Consumes everything. Integrates.

---

### Kernel 1: Bonding (The High-Precision Kernel)

**Optimization:**
*   **Grid Size:** `GlobalSize = Total_Nodes`. (Not Total Atoms).
*   **Workgroup:** Standard (e.g., 64 or 32).
*   **Memory:** We **do not** use Local Memory (LDS) for positions here.
    *   *Why?* Bond neighbors are specific (sparse). Loading a whole spatial tile (dense) is wasteful if we only need 4 specific atoms. Modern GPU L2 caches handle these "nearby" global reads very well since atoms are sorted.

```c
__kernel void solve_bonds_rigid(
    // Input
    __global float4* pos,
    __global float4* quat,
    __global int4*   bonds,       // Indices of neighbors
    __global float*  stiffness,
    // Output
    __global float4* dpos_node,   // Delta for me
    __global float4* drot_node,   // Torque for me
    __global float4* dpos_recoil  // Delta for neighbor (SCATTER WRITE)
    // No Ghost/Cluster info needed! Pure topology.
) {
    int gid = get_global_id(0);
    // No "if (is_cap)" check needed. We only launched threads for Nodes.

    float4 p_me = pos[gid];
    float4 q_me = quat[gid];
    int4 my_bonds = bonds[gid];

    // Accumulators
    float4 dp_sum = (float4)(0.0f);
    float4 dr_sum = (float4)(0.0f);

    // Loop over 4 potential bonds (unrolled)
    int neighbors[4] = {my_bonds.x, my_bonds.y, my_bonds.z, my_bonds.w};

    for(int i=0; i<4; i++) {
        int nid = neighbors[i];
        if (nid == -1) continue;

        // READ GLOBAL (Cache hit likely due to sorting)
        float4 p_neigh = pos[nid];
        float4 q_neigh = quat[nid]; // If needed

        // --- ARAP / PORT MATH HERE ---
        // Compute constraint C(p_me, p_neigh)
        // ...
        
        // WRITE RECOIL (The "Scatter")
        // We write to a specific slot reserved for this bond
        // e.g. dpos_recoil[ bond_id_offset + i ] = recoil_force;
        // OR better: use atomic_add if precision allows, 
        // OR use the 'bkSlots' logic you mentioned later in Apply.
        
        // For now, let's assume we output the node's contribution
        dp_sum += calculated_delta_pos;
        dr_sum += calculated_delta_rot;
    }

    dpos_node[gid] = dp_sum;
    drot_node[gid] = dr_sum;
}
```

---

### Kernel 2: Collisions (The Bandwidth Kernel)

**Optimization:**
*   **Grid Size:** `GlobalSize = Total_Clusters * 64`.
*   **Memory:** Heavy use of **Local Memory (LDS)**.
*   **Exclusions:** We load a compressed exclusion bitmask or list into **Registers**.

```c
__kernel void solve_collisions_rigid(
    // Input
    __global float4* pos,
    __global int*    ghost_indices, // From make_local_topology
    __global int*    ghost_counts,
    __global int4*   exclusions,    // 1st & 2nd neighbors to ignore
    // Output
    __global float4* dpos_coll      // Collision delta
) {
    // 1. Setup Local Memory
    __local float4 l_pos[64 + 32]; // Cluster + Ghosts
    
    int gid = get_global_id(0);
    int lid = get_local_id(0);     // 0..63
    int grp = get_group_id(0);

    // 2. Cooperative Load (Cluster + Ghosts)
    // ... (Same logic as before) ...
    barrier(CLK_LOCAL_MEM_FENCE);

    // 3. Load Exclusions to REGISTERS
    // We need to ignore bond neighbors (1st) and angles (2nd)
    int4 my_excl = exclusions[gid]; 
    // If you have more than 4, use two int4s. Registers are fast.

    float4 p_me = l_pos[lid];
    float4 dp = (float4)(0.0f);
    int n_ghosts = ghost_counts[grp];

    // 4. Brute Force against Local Cache
    // Iterate all potential colliders in LDS
    int total = 64 + n_ghosts;
    
    for (int j=0; j<total; j++) {
        if (j == lid) continue;

        // --- FAST EXCLUSION CHECK ---
        // Note: 'exclusions' stores GLOBAL indices. 
        // But 'j' is a LOCAL index. 
        // We need to map j -> global_id to compare? 
        // ERROR TRAP: This is the hard part of splitting.
        
        // BETTER STRATEGY: 
        // In 'make_local_topology', convert global exclusion IDs 
        // to LOCAL indices (if they are in the cluster/ghost list).
        // If an excluded atom is NOT in the ghost list, it's far away anyway!
        // So: 'exclusions' input should be mapped to LOCAL indices [-1 if not present].
        
        if (j == my_excl.x || j == my_excl.y || j == my_excl.z || j == my_excl.w) continue;

        float4 p_other = l_pos[j];
        
        // Sphere-Sphere Check
        float dist2 = dot(p_me.xyz - p_other.xyz, p_me.xyz - p_other.xyz);
        if (dist2 < (R_sum * R_sum)) {
            // Solve Collision
            dp += correction_vector;
        }
    }

    dpos_coll[gid] = dp;
}
```

---

### Kernel 3: Apply & Integrate

This kernel gathers the results. This is where your `bkSlots` (Back Slots) come in to handle the Recoil.

```c
__kernel void apply_updates_rigid(
    __global float4* pos,
    __global float4* quat,
    __global float4* dpos_node,     // From Kernel 1
    __global float4* drot_node,     // From Kernel 1
    __global float4* dpos_recoil,   // From Kernel 1 (The buffer neighbors wrote to)
    __global float4* dpos_coll,     // From Kernel 2
    __global int4*   bkSlots,       // "Who writes to me?" indices
    float dt
) {
    int gid = get_global_id(0);

    // 1. Gather Recoil
    // My neighbors calculated forces *on me* and stored them in dpos_recoil.
    // bkSlots tells me where to look.
    float4 recoil_sum = (float4)(0.0f);
    
    // Unroll lookups
    if (bkSlots[gid].x != -1) recoil_sum += dpos_recoil[bkSlots[gid].x];
    if (bkSlots[gid].y != -1) recoil_sum += dpos_recoil[bkSlots[gid].y];
    // ... etc

    // 2. Sum Everything
    float4 total_dpos = dpos_node[gid] + recoil_sum + dpos_coll[gid];
    float4 total_drot = drot_node[gid]; // Recoil might generate torque too, handled similarly

    // 3. XPBD Update / Euler Integration
    pos[gid] += total_dpos;
    quat[gid] = integrate_quat(quat[gid], total_drot);
}
```

---

### Critical Implementation Detail: The Exclusion Mapping

In the collision kernel, I noted a tricky part: **Exclusion Indices.**
*   The `exclusions` array usually stores **Global IDs** (e.g., "Atom 500 should not hit Atom 501").
*   But inside `solve_collisions_rigid`, we loop over **Local Indices** (0..96).
*   Comparing `Global_ID` vs `Local_ID` requires looking up the global ID of the local atom `j` every time. `if (local_to_global[j] == excluded_global_id)`.
*   **Performance Hit:** This adds an indirect memory read inside the inner collision loop. Bad.

**Solution:** Update your `make_local_topology` kernel.
When you build the `ghost_indices` list, you should also produce a **`local_exclusions`** buffer.
1.  Read `global_exclusion_id` for Atom I.
2.  Is `global_exclusion_id` inside my Cluster? If yes, store its `local_id` (0..63).
3.  Is it inside my Ghost List? If yes, store its ghost `local_id` (64+).
4.  Is it nowhere? Store `-1` (It's too far away to collide anyway, so we don't need to exclude it).

This allows the Collision Kernel to check `if (j == my_local_exclusion.x)` purely in registers.

### Summary of Recommendation

1.  **Use Split Kernels.** It solves your Recoil synchronization problem and handles the heterogeneity of Node vs. Cap atoms.
2.  **Kernel 1 (Bonding):** Use Global Memory reads. Run only on Nodes.
3.  **Kernel 2 (Collision):** Use Local Memory (LDS). Run on Cluster Groups. Use **Mapped Local Exclusions** (Computed during topology step).
4.  **Sorting:** Global Cluster Sort `[Nodes | Caps]` is essential. Keep it.

This approach is modular, easier to debug, and likely faster because it maximizes register occupancy for the complex bonding math and memory bandwidth for the collision checks.