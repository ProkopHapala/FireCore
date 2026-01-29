// XPDB_new.cl - Stripped down Position Based Dynamics with Tiled Jacobi
// Only essential kernels: bounding boxes, ghost building, and single Jacobi iteration

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define GROUP_SIZE     64      // Workgroup size = atoms per cluster
#define MAX_NEIGH_COLL 64      // Max collision neighbors for debug output
#define N_MAX_BONDED   16      // Max bonded neighbors per atom (fixed-size)
#define MAX_GHOSTS     128     // Max external atoms per cluster

// ------------------------------------------------------------------
// HELPER: Bounding Box Intersection
// ------------------------------------------------------------------
bool bboxes_overlap(float4 minA, float4 maxA, float4 minB, float4 maxB, float margin) {
    if (maxA.x + margin < minB.x || minA.x > maxB.x + margin) return false;
    if (maxA.y + margin < minB.y || minA.y > maxB.y + margin) return false;
    if (maxA.z + margin < minB.z || minA.z > maxB.z + margin) return false;
    return true;
}

// ------------------------------------------------------------------
// KERNEL 1: Update Bounding Boxes (Reduction)
// ------------------------------------------------------------------
// Assumes 1 WorkGroup = 1 Atom Group (Cluster)
// Local memory size must be at least GROUP_SIZE * sizeof(float4)
__kernel void update_bboxes(
    __global const float4* curr_pos,
    __global const float4* params, // .x = radius
    __global float4* bboxes_min,
    __global float4* bboxes_max,
    __local float4* local_min,
    __local float4* local_max,
    int num_atoms
) {
    int lid = get_local_id(0);
    int gid = get_global_id(0);
    int group_id = get_group_id(0);

    // 1. Load into local memory
    float4 p = (gid < num_atoms) ? curr_pos[gid] : (float4)(0.0f);
    float r = (gid < num_atoms) ? params[gid].x : 0.0f;
    
    // Init min/max with atom extent (pos +/- radius)
    // If out of bounds, set to infinity/-infinity to be ignored
    if (gid < num_atoms) {
        local_min[lid] = (float4)(p.x - r, p.y - r, p.z - r, 0.0f);
        local_max[lid] = (float4)(p.x + r, p.y + r, p.z + r, 0.0f);
    } else {
        local_min[lid] = (float4)(1e10f);
        local_max[lid] = (float4)(-1e10f);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // 2. Parallel Reduction (Tree)
    // GROUP_SIZE must be power of 2 for this simple implementation
    for (int stride = GROUP_SIZE / 2; stride > 0; stride >>= 1) {
        if (lid < stride) {
            local_min[lid] = min(local_min[lid], local_min[lid + stride]);
            local_max[lid] = max(local_max[lid], local_max[lid + stride]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // 3. Write result (Thread 0 of group)
    if (lid == 0) {
        bboxes_min[group_id] = local_min[0];
        bboxes_max[group_id] = local_max[0];
    }
}

// ------------------------------------------------------------------
// KERNEL 2: Build Local Topology (Ghost Discovery + Bond Re-indexing)
// ------------------------------------------------------------------
// Work distribution: 1 WorkGroup = 1 Cluster
__kernel void build_local_topology(
    // Input
    __global const float4* curr_pos,
    __global const float4* bboxes_min,
    __global const float4* bboxes_max,
    __global const int* bond_indices_global,   // [num_atoms * N_MAX_BONDED], -1 for no bond
    
    // Output
    __global int* ghost_indices_flat,          // [Cluster * MAX_GHOSTS + k]
    __global int* ghost_counts,                // [Cluster]
    __global int* bond_indices_local,         // [num_atoms * N_MAX_BONDED]
    
    int num_atoms,
    int num_groups,
    float margin_sq, // (2*Rmax)^2 used for BBox check
    float bbox_margin
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    
    // Local storage for the ghost list being built
    // We build it in Local Mem first for speed, then flush to global
    __local int l_ghost_list[MAX_GHOSTS];
    __local int l_ghost_counter;
    __local float4 l_my_bbox_min;
    __local float4 l_my_bbox_max;

    if (lid == 0) {
        l_ghost_counter = 0;
        l_my_bbox_min = bboxes_min[grp];
        l_my_bbox_max = bboxes_max[grp];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    float4 my_min = l_my_bbox_min;
    float4 my_max = l_my_bbox_max;

    // ------------------------------------------------
    // STEP 1: Find Ghosts (Cooperative)
    // ------------------------------------------------
    // Iterate over all other groups
    for (int other_g = 0; other_g < num_groups; other_g++) {
        if (other_g == grp) continue;

        // Check BBox Overlap
        bool overlap = false;
        float4 o_min = bboxes_min[other_g];
        float4 o_max = bboxes_max[other_g];
        
        if (my_max.x + bbox_margin >= o_min.x && my_min.x <= o_max.x + bbox_margin &&
            my_max.y + bbox_margin >= o_min.y && my_min.y <= o_max.y + bbox_margin &&
            my_max.z + bbox_margin >= o_min.z && my_min.z <= o_max.z + bbox_margin) {
            overlap = true;
        }

        // If overlap, all threads load atoms from 'other_g' and check distance
        if (overlap) {
            int global_idx = other_g * GROUP_SIZE + lid;
            
            if (global_idx < num_atoms) {
                float4 p = curr_pos[global_idx];
                
                // Distance to My BBox (Clamped Distance)
                float dx = max(0.0f, max(my_min.x - p.x, p.x - my_max.x));
                float dy = max(0.0f, max(my_min.y - p.y, p.y - my_max.y));
                float dz = max(0.0f, max(my_min.z - p.z, p.z - my_max.z));
                float dist_sq = dx*dx + dy*dy + dz*dz;
                
                // If close enough to potentially collide
                if (dist_sq < margin_sq) {
                    int slot = atomic_inc(&l_ghost_counter);
                    if (slot < MAX_GHOSTS) {
                        l_ghost_list[slot] = global_idx;
                    }
                }
            }
        }
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Cap counter
    int total_ghosts = min(l_ghost_counter, MAX_GHOSTS);

    // ------------------------------------------------
    // STEP 2: Flush Ghosts to Global
    // ------------------------------------------------
    int base_offset = grp * MAX_GHOSTS;
    for (int i = lid; i < total_ghosts; i += GROUP_SIZE) {
        ghost_indices_flat[base_offset + i] = l_ghost_list[i];
    }
    if (lid == 0) ghost_counts[grp] = total_ghosts;

    // ------------------------------------------------
    // STEP 3: Re-Index Bonds (The "Map")
    // ------------------------------------------------
    // Each thread handles its own atom's bonds
    int my_global_id = grp * GROUP_SIZE + lid;
    
    if (my_global_id < num_atoms) {
        // Clear local bond indices
        for (int slot = 0; slot < N_MAX_BONDED; slot++) {
            bond_indices_local[my_global_id * N_MAX_BONDED + slot] = -1;
        }

        // Process each bond slot
        for (int slot = 0; slot < N_MAX_BONDED; slot++) {
            int target = bond_indices_global[my_global_id * N_MAX_BONDED + slot];
            if (target == -1) break;
            
            int target_grp = target / GROUP_SIZE;
            if (target_grp == grp) {
                // Internal: Simple math
                bond_indices_local[my_global_id * N_MAX_BONDED + slot] = target % GROUP_SIZE;
            } else {
                // External: Search Ghost List
                // Linear search is OK because list is small and in Shared Mem
                for (int k = 0; k < total_ghosts; k++) {
                    if (l_ghost_list[k] == target) {
                        bond_indices_local[my_global_id * N_MAX_BONDED + slot] = GROUP_SIZE + k; // Offset by Internal Size
                        break;
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------------
// KERNEL 3: Tiled Jacobi Solver (Single Iteration)
// ------------------------------------------------------------------
__kernel void solve_cluster_jacobi(
    __global float4*       curr_pos,          // RW (in-place Jacobi iterate)
    __global float4*       prev_pos,          // RW (stores previous iterate, updated each launch)
    __global float4*       momentum,          // RW (stores heavy-ball displacement)
    __global const float4* params,            // R: .x = radius, .w = mass
    
    // Topology (Fixed-size bond buffers)
    __global const int*   bond_indices_local,   // [num_atoms * N_MAX_BONDED]
    __global const float* bond_lengths,         // [num_atoms * N_MAX_BONDED]
    __global const float* bond_stiffness,       // [num_atoms * N_MAX_BONDED]
    
    // Ghost Data
    __global const int* ghost_indices_flat,
    __global const int* ghost_counts,

    int num_atoms,
    int inner_iters,
    float dt,
    float k_coll,
    float omega
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float4 l_pos_new[GROUP_SIZE];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];

    // Load Internal
    if (my_global_id < num_atoms) {
        l_pos[lid] = curr_pos[my_global_id];
    } else {
        l_pos[lid] = (float4)(0.0f);
    }

    // Load Constants
    float my_mass  = (my_global_id < num_atoms) ? params[my_global_id].w : 1.0f;
    float my_rad   = (my_global_id < num_atoms) ? params[my_global_id].x : 0.0f;
    float alpha    = my_mass / (dt*dt);

    if (my_global_id < num_atoms) {
        l_rad[lid] = my_rad;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int iter = 0; iter < inner_iters; iter++) {
        // Load Ghosts
        int g_count = ghost_counts[grp];
        int g_offset = grp * MAX_GHOSTS;
        for (int k = lid; k < g_count; k += GROUP_SIZE) {
            int g_idx = ghost_indices_flat[g_offset + k];
            l_pos[GROUP_SIZE + k] = curr_pos[g_idx];
            l_rad[GROUP_SIZE + k] = params[g_idx].x;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        if (my_global_id < num_atoms) {
            float3 p = l_pos[lid].xyz;
            float3 rhs = p * alpha;
            float k_sum = alpha;

            // Bonds
            for (int slot = 0; slot < N_MAX_BONDED; slot++) {
                int idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
                if (idx == -1) break;
                
                float3 n_pos = l_pos[idx].xyz;
                float3 diff = p - n_pos;
                float dist = length(diff);
                if (dist > 1e-6f) {
                    float L = bond_lengths[my_global_id * N_MAX_BONDED + slot];
                    float st = bond_stiffness[my_global_id * N_MAX_BONDED + slot];
                    float coeff = st * (L / dist);
                    rhs += diff * coeff;
                    rhs += n_pos * st;
                    k_sum += st;
                }
            }

            // Collisions
            int check_count = GROUP_SIZE + g_count;
            for (int j = 0; j < check_count; j++) {
                if (j == lid) continue;
                
                // Skip if this is a bonded neighbor
                bool is_bonded = false;
                for (int slot = 0; slot < N_MAX_BONDED; slot++) {
                    int idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
                    if (idx == -1) break;
                    if (idx == j) {
                        is_bonded = true;
                        break;
                    }
                }
                if (is_bonded) continue;
                
                float4 other = l_pos[j];
                float3 diff = p - other.xyz;
                float dist_sq = dot(diff, diff);
                float r_sum = my_rad + l_rad[j];
                if (dist_sq < r_sum * r_sum && dist_sq > 1e-12f) {
                    float dist = sqrt(dist_sq);
                    float coeff = k_coll * (r_sum / dist);
                    rhs += diff * coeff;
                    rhs += other.xyz * k_coll;
                    k_sum += k_coll;
                }
            }

            float3 p_new = rhs / k_sum;
            l_pos_new[lid] = (float4)(p_new, 0.0f);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        
        if (my_global_id < num_atoms) {
            l_pos[lid] = l_pos_new[lid];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (my_global_id < num_atoms) {
        float3 p_corr = l_pos[lid].xyz;
        float3 p_prev = prev_pos[my_global_id].xyz;
        float3 dp_prev = momentum[my_global_id].xyz;
        float3 p_new = p_corr + dp_prev * omega;
        float3 dp_new = p_new - p_prev;
        momentum[my_global_id] = (float4)(dp_new, 0.0f);
        prev_pos[my_global_id] = (float4)(p_new, 0.0f);
        curr_pos[my_global_id] = (float4)(p_new, 0.0f);
    }
}

// ------------------------------------------------------------------
// KERNEL 3b: Local multi-iteration Jacobi WITHOUT collisions
// ------------------------------------------------------------------
// WARNING: intended only for isolated small systems where num_atoms <= GROUP_SIZE
// and where collision constraints are disabled/irrelevant. This kernel performs
// multiple inner iterations purely on bonds.
__kernel void solve_cluster_jacobi_nocoll(
    __global float4*       curr_pos,          // RW
    __global float4*       prev_pos,          // RW
    __global float4*       momentum,          // RW
    __global const float4* params,            // R: .w = mass

    __global const int*   bond_indices_local, // [num_atoms * N_MAX_BONDED]
    __global const float* bond_lengths,       // [num_atoms * N_MAX_BONDED]
    __global const float* bond_stiffness,     // [num_atoms * N_MAX_BONDED]

    int num_atoms,
    int inner_iters,
    float dt,
    float omega
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE];
    __local float4 l_pos_new[GROUP_SIZE];

    if (my_global_id < num_atoms) {
        l_pos[lid] = curr_pos[my_global_id];
    } else {
        l_pos[lid] = (float4)(0.0f);
    }

    float my_mass  = (my_global_id < num_atoms) ? params[my_global_id].w : 1.0f;
    float alpha    = my_mass / (dt*dt);

    barrier(CLK_LOCAL_MEM_FENCE);

    for (int iter = 0; iter < inner_iters; iter++) {
        if (my_global_id < num_atoms) {
            float3 p = l_pos[lid].xyz;
            float3 rhs = alpha * p;
            float k_sum = alpha;

            for (int slot = 0; slot < N_MAX_BONDED; slot++) {
                int idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
                if (idx == -1) break;
                float3 n_pos = l_pos[idx].xyz;
                float3 diff = p - n_pos;
                float dist = length(diff);
                if (dist > 1e-6f) {
                    float L = bond_lengths[my_global_id * N_MAX_BONDED + slot];
                    float st = bond_stiffness[my_global_id * N_MAX_BONDED + slot];
                    float coeff = st * (L / dist);
                    rhs += diff * coeff;
                    rhs += n_pos * st;
                    k_sum += st;
                }
            }

            float3 p_new = rhs / k_sum;
            l_pos_new[lid] = (float4)(p_new, 0.0f);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        if (my_global_id < num_atoms) {
            l_pos[lid] = l_pos_new[lid];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (my_global_id < num_atoms) {
        curr_pos[my_global_id] = l_pos[lid];
    }
}

// ------------------------------------------------------------------
// KERNEL 4: Bond residuals for fixed-slot global topology
// ------------------------------------------------------------------
// Writes absolute residual |dist - L0| for each (atom,slot) into out_abs_residual.
// Layout: out_abs_residual[ i*N_MAX_BONDED + slot ]
__kernel void bond_residuals_fixed_global(
    __global const float4* curr_pos,
    __global const int*    bond_indices_global,
    __global const float*  bond_lengths,
    __global float*        out_abs_residual,
    int num_atoms
){
    int i = get_global_id(0);
    if (i >= num_atoms) return;
    float3 p = curr_pos[i].xyz;
    int base = i * N_MAX_BONDED;
    for (int slot = 0; slot < N_MAX_BONDED; slot++) {
        int j = bond_indices_global[base + slot];
        if (j == -1) break;
        float3 q = curr_pos[j].xyz;
        float L0 = bond_lengths[base + slot];
        float r = length(p - q) - L0;
        out_abs_residual[base + slot] = fabs(r);
    }
    // Fill remaining slots with 0 to keep reductions clean
    for (int slot = 0; slot < N_MAX_BONDED; slot++) {
        if (bond_indices_global[base + slot] == -1) {
            out_abs_residual[base + slot] = 0.0f;
        }
    }
}

// ------------------------------------------------------------------
// KERNEL 5: Single-step tiled Jacobi with persistent prev_pos
// ------------------------------------------------------------------
// This is equivalent to one inner iteration of solve_cluster_jacobi, but with momentum
// based on a global prev_pos buffer so that momentum works across kernel launches.
__kernel void solve_cluster_jacobi_step(
    __global float4*       curr_pos,
    __global float4*       prev_pos,
    __global float4*       momentum,
    __global const float4* params,

    __global const int*   bond_indices_local, // [num_atoms * N_MAX_BONDED]
    __global const float* bond_lengths,       // [num_atoms * N_MAX_BONDED]
    __global const float* bond_stiffness,     // [num_atoms * N_MAX_BONDED]

    __global const int* ghost_indices_flat,
    __global const int* ghost_counts,

    int num_atoms,
    float dt,
    float k_coll,
    float omega
){
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];

    float my_mass  = (my_global_id < num_atoms) ? params[my_global_id].w : 1.0f;
    float my_rad   = (my_global_id < num_atoms) ? params[my_global_id].x : 0.0f;
    float alpha    = my_mass / (dt*dt);

    if (my_global_id < num_atoms) {
        l_pos[lid] = curr_pos[my_global_id];
        l_rad[lid] = my_rad;
    } else {
        l_pos[lid] = (float4)(0.0f);
        l_rad[lid] = 0.0f;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int g_idx = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = curr_pos[g_idx];
        l_rad[GROUP_SIZE + k] = params[g_idx].x;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        float3 p = l_pos[lid].xyz;
        float3 rhs = p * alpha;
        float k_sum = alpha;

        // Bonds
        for (int slot = 0; slot < N_MAX_BONDED; slot++) {
            int idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
            if (idx == -1) break;
            float3 n_pos = l_pos[idx].xyz;
            float3 diff = p - n_pos;
            float dist = length(diff);
            if (dist > 1e-6f) {
                float L = bond_lengths[my_global_id * N_MAX_BONDED + slot];
                float st = bond_stiffness[my_global_id * N_MAX_BONDED + slot];
                float coeff = st * (L / dist);
                rhs += diff * coeff;
                rhs += n_pos * st;
                k_sum += st;
            }
        }

        // Collisions
        if (k_coll != 0.0f) {
            int check_count = GROUP_SIZE + g_count;
            for (int j = 0; j < check_count; j++) {
                if (j == lid) continue;
                bool is_bonded = false;
                for (int slot = 0; slot < N_MAX_BONDED; slot++) {
                    int idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
                    if (idx == -1) break;
                    if (idx == j) { is_bonded = true; break; }
                }
                if (is_bonded) continue;
                float3 diff = p - l_pos[j].xyz;
                float dist_sq = dot(diff, diff);
                float r_sum = my_rad + l_rad[j];
                if (dist_sq < r_sum * r_sum && dist_sq > 1e-12f) {
                    float dist = sqrt(dist_sq);
                    float coeff = k_coll * (r_sum / dist);
                    rhs += diff * coeff;
                    rhs += l_pos[j].xyz * k_coll;
                    k_sum += k_coll;
                }
            }
        }

        float3 p_corr = rhs / k_sum;
        float3 p_prev = prev_pos[my_global_id].xyz;
        float3 dp_prev = momentum[my_global_id].xyz;
        float3 p_new = p_corr + dp_prev * omega;
        float3 dp_new = p_new - p_prev;
        momentum[my_global_id] = (float4)(dp_new, 0.0f);
        prev_pos[my_global_id] = (float4)(p_new, 0.0f);
        curr_pos[my_global_id] = (float4)(p_new, 0.0f);
    }
}
