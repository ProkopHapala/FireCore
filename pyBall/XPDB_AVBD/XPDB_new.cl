// XPDB_new.cl - Stripped down Position Based Dynamics with Tiled Jacobi
// Only essential kernels: bounding boxes, ghost building, and single Jacobi iteration

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define GROUP_SIZE     64      // Workgroup size = atoms per cluster
#define MAX_NEIGH_COLL 64      // Max collision neighbors for debug output
#define N_MAX_BONDED   16      // Max bonded neighbors per atom (fixed-size)
#define MAX_GHOSTS     128     // Max external atoms per cluster

#define N_PORT_TYPES 4

// ------------------------------------------------------------------
// HELPER: Bounding Box Intersection
// ------------------------------------------------------------------
bool bboxes_overlap(float4 minA, float4 maxA, float4 minB, float4 maxB, float margin) {
    if (maxA.x + margin < minB.x || minA.x > maxB.x + margin) return false;
    if (maxA.y + margin < minB.y || minA.y > maxB.y + margin) return false;
    if (maxA.z + margin < minB.z || minA.z > maxB.z + margin) return false;
    return true;
}

// inline void porttab_load(__local float4* ldirs, __local uchar* lns, __global const float4* gdirs, __global const uchar* gns);
// inline int  porttab_n  (__local const uchar* lns, int typ);
// inline float3 porttab_dir(__local const float4* ldirs, int typ, int k);


inline void porttab_load(__local float4* ldirs, __local uchar* lns, __global const float4* gdirs, __global const uchar* gns){
    int lid = get_local_id(0);
    int lsz = get_local_size(0);
    for (int i = lid; i < N_PORT_TYPES*4; i += lsz) { ldirs[i] = gdirs[i]; }
    for (int i = lid; i < N_PORT_TYPES;   i += lsz) { lns[i]   = gns[i];   }
    barrier(CLK_LOCAL_MEM_FENCE);
}

inline int    porttab_n  (__local const uchar* lns, int typ){ return (typ>=0 && typ<N_PORT_TYPES) ? (int)lns[typ] : 0; }
inline float3 porttab_dir(__local const float4* ldirs, int typ, int k){ return ldirs[typ*4 + k].xyz; }


inline float4 quat_from_vec(float3 w);
inline float4 quat_mul(float4 a, float4 b);

// Safe component fetch from __global float4 without address-space cast (Intel OpenCL complains otherwise)
inline float fetch4(__global const float4* arr, int idx, int k){
    float4 v = arr[idx];
    // TODO: replace cascading if with vector load macros if we need speed; kept explicit for debuggability
    if(k==0) return v.x;
    if(k==1) return v.y;
    if(k==2) return v.z;
    return v.w;
}


// Rotate vector v by quaternion q
inline float3 quat_rotate(float4 q, float3 v) {
    float3 t = 2.0f * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}

// Convert small rotation vector to quaternion delta
inline float4 quat_from_vec(float3 w) {
    float angle = length(w);
    if (angle < 1e-8f) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    float s = sin(angle * 0.5f) / angle;
    return (float4)(w * s, cos(angle * 0.5f));
}

// Quaternion from axis-angle
inline float4 quat_from_axis_angle(float3 axis, float angle) {
    float a = length(axis);
    if (a < 1e-8f || fabs(angle) < 1e-8f) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    float3 n = axis / a;
    float s = sin(angle * 0.5f);
    return (float4)(n * s, cos(angle * 0.5f));
}

// Quaternion multiplication
inline float4 quat_mul(float4 a, float4 b) {
    return (float4)(
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    );
}

inline int port_n(int typ){
    if(typ==1) return 2;   // sp1
    if(typ==2) return 3;   // sp2
    if(typ==3) return 4;   // sp3
    return 0;
}

__kernel void clear_rigid_forces(
    const int natoms,
    __global float4* force
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    force[i]  = (float4)(0.0f);
}

__kernel void clear_rigid_node_buffers(
    const int nnode,
    __global float4* fneigh,
    __global float4* pneigh
) {
    int i = get_global_id(0);
    int n = nnode * 4;
    if (i >= n) return;
    fneigh[i] = (float4)(0.0f);
    pneigh[i] = (float4)(0.0f);
}


__kernel void gather_port_forces(
    const int nnode,
    __global const float4* pos,
    __global const float4* quat,
    __global int4*   neighs,
    __global float4* bKs,
    __global const float4* port_local,
    __global const uchar*  port_n,
    __global float4* force,
    __global float4* fneigh,
    __global float4* pneigh
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float3 p_i = pos[i].xyz;
    float4 q_i = quat[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;

    int npi = (int)port_n[i];

    float3 fi = (float3)(0.0f);

    int i4 = i * 4;
    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float Kij = fetch4(bKs, i, k);
        float3 p_j = pos[j].xyz;

        float3 r_arm = (float3)(0.0f);
        if (k < npi) r_arm = quat_rotate(q_i, port_local[i4 + k].xyz);

        float3 tip = p_i + r_arm;
        float3 f   = (p_j - tip) * Kij;
        pneigh[i4 + k] = (float4)(r_arm, 0.0f);

        fi += f;
        fneigh[i4 + k] = (float4)(-f, 0.0f);
    }

    force[i]  = (float4)(fi, 0.0f);
}


__kernel void integrate_rigid_explicit(
    const int natoms,
    const int nnode,
    __global float4* pos,
    __global float4* vel,
    __global float4* quat,
    __global float4* omega,
    __global int4*   bkSlots,
    __global float4* force,
    __global float4* fneigh,
    __global float4* pneigh,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float invMass = pos[i].w;
    float mass    = (invMass > 1e-12f) ? (1.0f / invMass) : 1e12f;
    float invI    = 1.0f / (0.4f * mass + 1e-12f);

    float3 f = force[i].xyz;

    int4 bk = bkSlots[i];
    int* islots = (int*)&bk;
    for (int k = 0; k < 4; k++) {
        int islot = islots[k];
        if (islot >= 0) {
            f += fneigh[islot].xyz;
        }
    }

    float3 v = vel[i].xyz;
    v = v * damp + (f * invMass) * dt;
    float3 p = pos[i].xyz + v * dt;
    vel[i] = (float4)(v, 0.0f);
    pos[i] = (float4)(p, invMass);

    if (i < nnode) {
        float4 q = quat[i];
        float3 tau = (float3)(0.0f);
        int i4 = i * 4;
        for (int k = 0; k < 4; k++) {
            float3 r_arm = pneigh[i4 + k].xyz;
            float3 Fk    = -fneigh[i4 + k].xyz;
            tau += cross(r_arm, Fk);
        }

        float3 w = omega[i].xyz;
        w = w * damp + (tau * invI) * dt;
        omega[i] = (float4)(w, 0.0f);
        float4 dq = quat_from_vec(w * dt);
        quat[i] = normalize(quat_mul(dq, q));
    }
}

// forward declarations (used by apply_deltas_rigid before helper block below)
// inline float4 quat_from_vec(float3 w);
// inline float4 quat_mul(float4 a, float4 b);


__kernel void apply_deltas_rigid(
    const int natoms,
    const int nnode,
    __global float4* pos,
    __global float4* quat,
    __global float4* delta_p,
    __global float4* delta_rot,
    const float relax_p,
    const float relax_q
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float4 p4 = pos[i];
    p4.xyz += delta_p[i].xyz * relax_p;
    pos[i] = p4;

    if (i < nnode) {
        float4 q = quat[i];
        float3 drot = delta_rot[i].xyz * relax_q;
        float4 dq = quat_from_vec(drot);
        quat[i] = normalize(quat_mul(dq, q));
    }
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



// ==================
// ==================
//   Rigid Atom Rotation
// ==================
// ==================




// =============================================================
// MATH HELPERS
// =============================================================



inline float3 port_dir(int typ, int k){
    // Unit directions in LOCAL frame
    if(typ==1){ // sp1
        if(k==0) return (float3)( 1.0f, 0.0f, 0.0f);
        if(k==1) return (float3)(-1.0f, 0.0f, 0.0f);
        return (float3)(0.0f);
    }
    if(typ==2){ // sp2
        if(k==0) return (float3)( 1.0f,  0.0f, 0.0f);
        if(k==1) return (float3)(-0.5f,  0.8660254f, 0.0f);
        if(k==2) return (float3)(-0.5f, -0.8660254f, 0.0f);
        if(k==3) return (float3)(0.0f, 0.0f, 1.0f);  // optional pi
        return (float3)(0.0f);
    }
    if(typ==3){ // sp3 tetrahedron
        // 4 corners of tetrahedron: (±1,±1,±1)/sqrt(3)
        const float s = 0.57735026919f;
        if(k==0) return (float3)( s, s, s);
        if(k==1) return (float3)( s,-s,-s);
        if(k==2) return (float3)(-s, s,-s);
        if(k==3) return (float3)(-s,-s, s);
        return (float3)(0.0f);
    }
    return (float3)(0.0f);
}

// =============================================================
// THE KERNEL (Ping-Pong)
// =============================================================



// Rotate vector v by quaternion q
// inline float3 quat_rotate(float4 q, float3 v) {
//     float3 t = 2.0f * cross(q.xyz, v);
//     return v + q.w * t + cross(q.xyz, t);
// }

// Convert small rotation vector to quaternion delta
// inline float4 quat_from_vec(float3 w) {
//     float angle = length(w);
//     if (angle < 1e-8f) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
//     float s = sin(angle * 0.5f) / angle;
//     return (float4)(w * s, cos(angle * 0.5f));
// }

// // Quaternion multiplication
// inline float4 quat_mul(float4 a, float4 b) {
//     return (float4)(
//         a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
//         a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
//         a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
//         a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
//     );
// }

// Helper: XPBD Compliance (alpha) = 1 / (stiffness * dt^2)
inline float get_alpha(float stiffness, float dt) {
    return 1.0f / (stiffness * dt * dt + 1e-12f);
}

__kernel void xpbd_predict_and_external_forces(
    const int natoms,
    const int nnode,
    __global float4* apos,      // Positions [atoms] + [pi-vectors]
    __global float4* avel,      // Velocities [linear] + [angular]
    __global float4* aforce,    // External forces (Electrostatics, Pi-Pi, VdW)
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms + nnode) return;

    bool is_pi = i >= natoms;
    
    float4 p = apos[i];
    float4 v = avel[i];
    float4 f = aforce[i]; // Calculated by your existing evalPiAling, VdW, etc.
    
    if (!is_pi) {
        // --- Linear Prediction ---
        // v = v + f * dt / m (assuming mass=1 for simplicty, add invMass array if needed)
        v.xyz += f.xyz * dt; 
        
        // p_pred = p + v * dt
        p.xyz += v.xyz * dt;
        p.w = 0.0f; // Clear w for accumulation later if needed
    } else {
        // --- Angular Prediction (Vector-based) ---
        // Simple Euler integration for the normal vector n
        // v.xyz here stores angular velocity omega
        v.xyz += f.xyz * dt; // Torque -> Omega
        
        // Rotate n by omega * dt
        // Using your taylor expansion or simple cross product approx for small steps
        float3 w = v.xyz * dt;
        float angle = length(w);
        if (angle > 1e-8f) {
             float3 axis = w / angle;
             // Rodrigues rotation or simplified
             // For PBD stability, simple addition + normalize often suffices for prediction
             float3 p_cross_w = cross(p.xyz, w);
             p.xyz += p_cross_w; 
        }
        p.xyz = normalize(p.xyz);
    }

    // Write predicted positions/orientations back
    apos[i] = p;
    avel[i] = v;
}

__kernel void xpbd_solve_constraints_jacobi(
    const int4 nDOFs,            // x=natoms, y=nnode
    __global float4* apos,       // Predicted positions/vectors
    __global float4* delta_p,    // Output: Position corrections (accumulator)
    __global int*    delta_n,    // Output: Denominator (constraint count) for averaging
    __global int4*   neighs,
    __global float4* bLs,        // Bond Lengths
    __global float4* bKs,        // Bond Stiffness
    __global float4* Ksp,        // Pi-Sigma Stiffness
    __global float4* apars,      // Angle Params
    const float dt
) {
    int i = get_global_id(0); // Atom Index
    int natoms = nDOFs.x;
    int nnode = nDOFs.y;

    if (i >= natoms) return; // Only process atoms, not pi-vectors directly here

    float3 pi = apos[i].xyz;     // My Position
    float3 dp = (float3)(0.0f);  // Accumulator for my position correction
    int count = 0;               // How many constraints affected me

    // 1. Load my Pi-Vector (if I am a node)
    bool i_is_node = i < nnode;
    float3 ni = (float3)(0.0f);
    if (i_is_node) ni = apos[natoms + i].xyz;

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float4 bond_lens = bLs[i];
    float4 bond_stiff = bKs[i];
    float4 k_pi_sigma = Ksp[i];

    // ITERATE NEIGHBORS
    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float3 pj = apos[j].xyz;
        float3 r_ij = pj - pi;
        float d_ij = length(r_ij);
        float3 dir = r_ij / d_ij;

        // --- A. Bond Constraint (Distance) ---
        {
            float rest_len = ((float*)&bond_lens)[k];
            float stiffness = ((float*)&bond_stiff)[k];
            float alpha = get_alpha(stiffness, dt);
            
            // C = d_ij - rest
            float C = d_ij - rest_len;
            
            // XPBD Lambda (assuming unit masses for simplicity)
            // w_i = 1, w_j = 1. denominator = w_i + w_j + alpha
            float lambda = -C / (2.0f + alpha);
            
            // Apply to me (pi)
            dp -= dir * lambda; 
            count++;
        }

        // --- B. Pi-Sigma Constraint (Orthogonalization) ---
        // This is the implementation of your "Vector per atom" idea
        
        // 1. My Pi-Vector vs This Bond (Aligns ME to plane of neighbor)
        if (i_is_node) {
            float stiffness = ((float*)&k_pi_sigma)[k];
            float alpha = get_alpha(stiffness, dt);

            // Constraint: ni . r_ij = 0
            float C = dot(ni, r_ij);
            
            // Gradients
            // grad_pi = -ni
            // grad_ni = r_ij
            
            // Inverse Inertia for Normal Vector (w_n)
            // Treat vector tip as particle with mass ~1? 
            // Better: w_n = 1.0 / (moment_of_inertia). Let's assume w_n = 1.0
            float w_p = 1.0f; 
            float w_n = 1.0f; 
            
            // Denominator: w_p * |grad_p|^2 + w_n * |grad_n|^2 + alpha
            // |grad_p|^2 = |-ni|^2 = 1
            // |grad_n|^2 = |r_ij|^2 = d_ij*d_ij
            float denom = w_p * 1.0f + w_n * (d_ij * d_ij) + alpha;
            
            float lambda = -C / denom;

            // Apply to me (Position update)
            // delta_pi = w_p * grad_pi * lambda = 1.0 * (-ni) * lambda
            dp += -ni * lambda;
            
            // Apply to my Pi-Vector (Orientation update)
            // We need to write this to a separate buffer because multiple neighbors affect my Normal
            // To avoid atomics on floats, we can recompute this in a separate "Update Pi" kernel 
            // OR use atomic_add (available in OpenCL 2.0 or via extension)
            // For now, let's just focus on position p.
            count++;
        }

        // 2. Neighbor's Pi-Vector vs This Bond (Aligns ME to plane of neighbor)
        if (j < nnode) {
            float3 nj = apos[natoms + j].xyz;
            // Stiffness? Need to look up J's stiffness. 
            // Approximation: Use my stiffness for symmetry or look up global array.
            // Let's assume symmetric stiffness.
            float stiffness = ((float*)&k_pi_sigma)[k]; 
            float alpha = get_alpha(stiffness, dt);

            // Constraint: nj . (pi - pj) = 0  => nj . (-r_ij) = 0
            float C = dot(nj, -r_ij); // = - dot(nj, r_ij)
            
            // Gradient w.r.t pi (Me):
            // grad_pi = nj
            
            float denom = 1.0f + 1.0f * (d_ij * d_ij) + alpha; // Assuming neighbors w_p = 1
            float lambda = -C / denom;
            
            // Apply to me
            dp += nj * lambda;
            count++;
        }
    }
    
    // Store Result
    delta_p[i] = (float4)(dp, 0.0f);
    delta_n[i] = count;
}


// Separate Kernel to update Pi-Vectors (Normals)
// This avoids race conditions. Each Pi-Vector looks at its own neighbors.
__kernel void xpbd_solve_pi_orientation(
    const int nnode,
    const int natoms,
    __global float4* apos,
    __global float4* delta_pi_vec, // Output: Delta for pi vector
    __global int4*   neighs,
    __global float4* Ksp,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float3 ni = apos[natoms + i].xyz;
    float3 pi = apos[i].xyz;
    float3 d_ni = (float3)(0.0f);

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float4 k_pi_sigma = Ksp[i];

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float3 pj = apos[j].xyz;
        float3 r_ij = pj - pi; // Vector from me to neighbor
        float d2 = dot(r_ij, r_ij);

        float stiffness = ((float*)&k_pi_sigma)[k];
        float alpha = get_alpha(stiffness, dt);

        // Constraint: ni . r_ij = 0
        float C = dot(ni, r_ij);

        // We only care about updating 'ni' here
        // Gradient w.r.t ni is r_ij
        
        // Denominator (same as in position kernel)
        float w_p = 1.0f;
        float w_n = 1.0f;
        float denom = w_p + w_n * d2 + alpha;

        float lambda = -C / denom;

        // delta_ni = w_n * grad_ni * lambda
        d_ni += r_ij * lambda;
    }

    delta_pi_vec[i] = (float4)(d_ni, 0.0f);
}

__kernel void xpbd_finalize_step(
    const int natoms,
    const int nnode,
    __global float4* apos,
    __global float4* avel,
    __global float4* delta_p,
    __global int*    delta_n,
    __global float4* delta_pi_vec,
    __global float4* apos_old, // Previous positions for velocity update
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms + nnode) return;

    if (i < natoms) {
        // --- Update Atom Position ---
        int count = delta_n[i];
        if (count > 0) {
            float3 dp = delta_p[i].xyz / (float)count; // Average the deltas (Jacobi)
            // SOR (Successive Over-Relaxation) factor could be applied here: dp *= 1.2f;
            apos[i].xyz += dp;
        }

        // --- Update Velocity (v = (p - p_old) / dt) ---
        float3 p_new = apos[i].xyz;
        float3 p_old = apos_old[i].xyz;
        avel[i].xyz = (p_new - p_old) / dt;
        
        // Update history
        apos_old[i] = apos[i];
        
    } else {
        // --- Update Pi-Vector ---
        int pi_idx = i - natoms;
        float3 dni = delta_pi_vec[pi_idx].xyz; // Usually averaged by neighbor count (4)
        // Since we didn't count neighbors in the Pi kernel, assume standard valence or just add
        // For relaxation, averaging is safer:
        dni /= 3.0f; // Approx valence
        
        float3 n_new = apos[i].xyz + dni;
        n_new = normalize(n_new); // CONSTANT LENGTH CONSTRAINT IS HARD CODED HERE
        
        // Update Angular Velocity (approx)
        // w ~ (n_new - n_old) / dt ... rough approx for vector dynamics
        float3 n_old = apos_old[i].xyz;
        avel[i].xyz = cross(n_old, n_new) / dt; // Approximate angular velocity
        
        apos[i].xyz = n_new;
        apos_old[i] = apos[i];
    }
}

// =============================================================
// HELPER FUNCTIONS (Inline)
// =============================================================

// Rotate vector v by quaternion q
// NOTE: quat_rotate is defined earlier in this file; keep single definition to avoid OpenCL redefinition.
#if 0
inline float3 quat_rotate(float4 q, float3 v) {
    float3 t = 2.0f * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}
#endif

// Inverse rotate: rotate v by conjugate of q
inline float3 quat_inv_rotate(float4 q, float3 v) {
    float4 q_inv = (float4)(-q.xyz, q.w); 
    return quat_rotate(q_inv, v);
}

// Convert a rotation vector (angle * axis) to a quaternion change
// Used for applying angular corrections
// NOTE: quat_from_vec/quat_mul are defined earlier in this file; keep single definition to avoid OpenCL redefinition.
#if 0
inline float4 quat_from_vec(float3 w) {
    float angle = length(w);
    if (angle < 1e-8f) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    float s = sin(angle * 0.5f) / angle;
    return (float4)(w * s, cos(angle * 0.5f));
}

inline float4 quat_mul(float4 a, float4 b) {
    return (float4)(
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    );
}
#endif

// =============================================================
// THE KERNEL
// =============================================================

__kernel void rigid_atom_xpbd_update(
    const int natoms,
    const int nnode, // Atoms with hybridization (sp2)
    __global float4* apos,      // Position (xyz) + invMass (w)
    __global float4* aquat,     // Quaternion (xyzw) - REPLACES PI VECTOR
    __global float4* avel,      // Linear Velocity
    __global float4* aomega,    // Angular Velocity
    __global int4*   neighs,    // Neighbor Indices
    __global float4* bLs,       // Bond Lengths (rest)
    __global float4* bKs,       // Bond Stiffness
    __global float4* bond_dirs, // **NEW**: Local direction of bonds (e.g., (1,0,0), (-0.5, 0.866,0)...)
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // 1. Load State
    float4 p4 = apos[i];
    float3 p = p4.xyz;
    float invMass = p4.w; // Store invMass in w for efficiency
    
    // Rigid body state
    float4 q = (float4)(0,0,0,1); 
    if (i < nnode) q = aquat[i]; // Only nodes have orientation

    float3 v = avel[i].xyz;
    float3 omega = aomega[i].xyz;

    // 2. Prediction Step (Integrate Velocities)
    // External forces (Gravity, etc) could be added here
    if (invMass > 0.0f) {
        p += v * dt;
        
        // Angular prediction (apply omega to quat)
        if (i < nnode) {
            float4 dq = quat_from_vec(omega * dt);
            q = normalize(quat_mul(dq, q));
        }
    }

    // 3. CONSTRAINT SOLVE (Jacobi Gather)
    // We accumulate corrections and apply them at the end of the kernel
    // No barriers needed because we read 'prev_p' (from neighbors) if we did strict PBD,
    // but for XPBD we can often use current predicted 'p' if we accept slight asymmetry 
    // or (better) we read the neighbors' predicted positions but DON'T write to them.
    
    float3 dp = (float3)(0.0f); // Position correction accumulator
    float3 domega = (float3)(0.0f); // Orientation correction accumulator (axis-angle)
    float w_sum = 0.0f; // To average position updates
    float w_rot_sum = 0.0f; // To average rotation updates

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    
    // If I am a node, I have up to 4 neighbors
    // We define local bond slots. 
    // Slot 0: (1,0,0), Slot 1: (-0.5, 0.866, 0), Slot 2: (-0.5, -0.866, 0) (Planar)
    // Slot 3: (0,0,1) (Pi-bond or capping)
    
    if (i < nnode) {
        float4 local_dirs_packed = bond_dirs[i]; // Let's assume this stores angles or we hardcode sp2
        // Hardcoded sp2 for example clarity:
        float3 local_vectors[3];
        local_vectors[0] = (float3)(1.0f, 0.0f, 0.0f);
        local_vectors[1] = (float3)(-0.5f, 0.866f, 0.0f);
        local_vectors[2] = (float3)(-0.5f, -0.866f, 0.0f);

        for (int k = 0; k < 3; k++) { // Loop over Sigma bonds
            int j = neighbors[k];
            if (j < 0) continue;

            // Read Neighbor Pos (Read-Only access to global memory)
            float3 pj = apos[j].xyz + avel[j].xyz * dt; // Use predicted pos
            
            float3 r_ij = pj - p; // Vector to neighbor
            float dist = length(r_ij);
            float3 dir = r_ij / dist;
            
            // --- A. Bond Constraint (Distance) ---
            // Standard PBD
            // float rest_len = ((float*)&bLs[i])[k]; // DEBUG: address-space cast was rejected by Intel
            float rest_len = fetch4(bLs, i, k);
            float compliance = 0.0001f; // inverse stiffness
            float C_dist = dist - rest_len;
            float lambda_dist = -C_dist / (invMass + apos[j].w + compliance);
            
            dp += -dir * lambda_dist * invMass; // Move ME towards optimal distance
            w_sum += 1.0f;

            // --- B. Angular/Planar Constraint (Orientation) ---
            // "My local vector[k], when rotated by Q, should align with direction to neighbor"
            
            float3 target_dir = dir; // The bond direction in world space
            float3 current_dir = quat_rotate(q, local_vectors[k]);
            
            // Cross product gives the rotation axis needed to align them
            float3 rot_axis = cross(current_dir, target_dir);
            float sin_angle = length(rot_axis);
            
            // Stiffness for angular part
            float alpha_rot = 0.1f; // Depends on your units
            
            // The correction:
            // We want to rotate 'current_dir' to match 'target_dir'.
            // Gradient w.r.t quaternion is related to the cross product.
            // Simplified: apply torque to align vectors.
            
            // Weighting: Inverse Inertia. Assume unit inertia for sphere approximation.
            float w_q = 1.0f; 
            
            // Correction magnitude
            // C = dot(current_dir, target_dir) - 1.0 (Cos constraint) or just angle
            // Simple approach: apply correction proportional to cross product
            float3 d_angle = rot_axis * (1.0f / (w_q + alpha_rot)); 
            
            domega += d_angle;
            w_rot_sum += 1.0f;
            
            // Note: In a full solver, this orientation change effectively applies a 
            // tangential force to the position too, but for relaxation, separating 
            // linear (A) and angular (B) is stable enough.
        }
    }

    // 4. Apply Corrections & Update State
    
    // Normalize updates (Jacobi Averaging)
    if (w_sum > 0) p += dp / w_sum; // Or just dp if stiffness is tuned for it
    if (w_rot_sum > 0 && i < nnode) {
        float3 axis_angle = domega / w_rot_sum;
        float4 dq = quat_from_vec(axis_angle);
        q = normalize(quat_mul(dq, q));
    }

    // 5. Finalize Velocity (XPBD/PBD Standard)
    avel[i].xyz = (p - p4.xyz) / dt;
    apos[i].xyz = p;
    
    if (i < nnode) {
        // Compute angular velocity from quaternion difference
        // w = 2 * dq * inv(q) / dt roughly
        float4 q_old = aquat[i];
        float4 q_diff = quat_mul(q, (float4)(-q_old.xyz, q_old.w)); // q * q_old^-1
        // Extract vector part roughly corresponds to omega * dt / 2
        aomega[i].xyz = q_diff.xyz * (2.0f / dt);
        if (q_diff.w < 0) aomega[i].xyz *= -1.0f; // Shortest path check
        
        aquat[i] = q;
    }
}


__kernel void rigid_solver_pingpong(
    const int natoms,
    const int nnode,
    // INPUT BUFFERS (Read-Only "Old State")
    __global float4* pos_in,    // xyz, w=invMass
    __global float4* quat_in,   // xyzw
    __global float4* vel_in,    // xyz (Linear Velocity)
    __global float4* omega_in,  // xyz (Angular Velocity)
    
    // OUTPUT BUFFERS (Write-Only "New State")
    __global float4* pos_out,
    __global float4* quat_out,
    __global float4* vel_out,
    __global float4* omega_out,

    // PARAMETERS
    __global int4*   neighs,    // Neighbors
    __global float4* bLs,       // Rest lengths
    __global float4* bKs,       // Stiffness
    __global float4* port_dirs, // Predefined local bond directions (The "Ports")
    
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // 1. LOAD STATE (Read from IN buffers)
    float4 p_curr = pos_in[i];
    float4 q_curr = (i < nnode) ? quat_in[i] : (float4)(0,0,0,1);
    float3 v      = vel_in[i].xyz;
    float3 omega  = omega_in[i].xyz;
    float invMass = p_curr.w;

    // 2. PREDICTION (Leap-Frog Stage 1)
    // Apply gravity or external forces here if needed to v/omega
    float3 p_pred = p_curr.xyz + v * dt;
    float4 q_pred = q_curr;
    
    if (i < nnode) {
        // Integrate Rotation: q_new = q_old + 0.5 * w * q_old * dt
        // Or simply multiply by rotation quaternion
        float4 dq = quat_from_vec(omega * dt);
        q_pred = normalize(quat_mul(dq, q_curr));
    }

    // 3. SOLVE CONSTRAINTS (Jacobi Gather)
    // We compute corrections based on the PREDICTED state of myself
    // and the PREDICTED state of neighbors (approximated by their current state + vel*dt)
    
    float3 dp = (float3)(0.0f);      // Position Delta Accumulator
    float3 d_phi = (float3)(0.0f);   // Rotation Delta Accumulator (Axis-Angle)
    float w_p_sum = 0.0f;            // Weight accumulator for position
    float w_q_sum = 0.0f;            // Weight accumulator for rotation

    if (i < nnode) {
        int4 ng = neighs[i];
        int* neighbors = (int*)&ng;
        
        // Load my local ports (Hardcoded sp2 example or loaded from memory)
        // Let's say port_dirs contains 3 vectors packed sequentially in global memory or constant memory
        // For efficiency, let's assume standard sp2 angles: 0, 120, 240 deg in XY plane
        float3 ports[3];
        ports[0] = (float3)(1.0f, 0.0f, 0.0f);
        ports[1] = (float3)(-0.5f, 0.866f, 0.0f);
        ports[2] = (float3)(-0.5f, -0.866f, 0.0f);

        for (int k = 0; k < 3; k++) {
            int j = neighbors[k];
            if (j < 0) continue;

            // READ NEIGHBOR (From IN buffer - Read Only!)
            float4 p_neigh_ref = pos_in[j];
            float3 v_neigh_ref = vel_in[j].xyz;
            
            // Predict neighbor position (Jacobi guess)
            float3 pj_pred = p_neigh_ref.xyz + v_neigh_ref * dt;

            // --- A. LINEAR CONSTRAINT (Bond Length) ---
            float3 r_ij = pj_pred - p_pred;
            float dist = length(r_ij);
            float3 dir = r_ij / dist;
            
            // float rest_len = ((float*)&bLs[i])[k]; // DEBUG
            // float stiffness = ((float*)&bKs[i])[k]; // DEBUG
            float rest_len = fetch4(bLs, i, k);
            float stiffness = fetch4(bKs, i, k);
            
            // Compliance alpha = 1 / (k * dt^2)
            float alpha = 1.0f / (stiffness * dt * dt + 1e-6f);
            
            float C = dist - rest_len;
            float w_j = p_neigh_ref.w;
            
            // XPBD Correction
            float lambda = -C / (invMass + w_j + alpha);
            dp += -dir * lambda * invMass; // Move ME
            w_p_sum += 1.0f;

            // --- B. ANGULAR CONSTRAINT (Port Alignment) ---
            // Rotate my local port into world space
            float3 port_world = quat_rotate(q_pred, ports[k]);
            
            // The bond direction is 'dir'. 
            // We want 'port_world' to align with 'dir'.
            // Cross product gives rotation axis & magnitude (sin theta)
            float3 rot_axis = cross(port_world, dir);
            
            // Stiffness for angle (can be different from bond)
            float alpha_rot = 0.0f; // Rigid alignment
            
            // Angular correction (simplified)
            // Ideally we use proper inertia tensor, but sphere approx (w=1) is fine here
            // w_rot + alpha
            float lambda_rot = 1.0f / (1.0f + alpha_rot); 
            d_phi += rot_axis * lambda_rot;
            w_q_sum += 1.0f;
        }
    }

    // 4. APPLY CORRECTIONS (Update Step)
    
    // Apply averaged position correction
    if (w_p_sum > 0) p_pred += dp; // / w_p_sum (if you want strict averaging, otherwise it's stiffer)
    
    // Apply averaged rotation correction
    if (w_q_sum > 0) {
        // d_phi is a sum of small rotation vectors.
        // We can apply it directly to the quaternion.
        float4 dq_correct = quat_from_vec(d_phi); // / w_q_sum
        q_pred = normalize(quat_mul(dq_correct, q_pred));
    }

    // 5. UPDATE VELOCITIES & WRITE OUTPUT
    // v_new = (p_new - p_old) / dt
    
    float3 v_new = (p_pred - p_curr.xyz) / dt;
    
    // omega_new approx (diff between quaternions)
    float4 q_diff = quat_mul(q_pred, (float4)(-q_curr.xyz, q_curr.w)); // q_new * q_old^-1
    float3 omega_new = q_diff.xyz * 2.0f / dt;
    if (q_diff.w < 0) omega_new *= -1.0f; // prevent shortest path flip

    // WRITE TO OUT BUFFERS (No Race Conditions!)
    pos_out[i]   = (float4)(p_pred, invMass);
    vel_out[i]   = (float4)(v_new, 0.0f);
    
    if (i < nnode) {
        quat_out[i]  = q_pred;
        omega_out[i] = (float4)(omega_new, 0.0f);
    }
}


// =============================================================
// UNIFIED SOLVER KERNEL
// =============================================================
__kernel void solve_rigid_symmetry(
    const int natoms,
    const int nnode,
    // READ-ONLY OLD STATE
    __global float4* pos_in,   // xyz, w=invMass
    __global float4* quat_in,  // xyzw
    __global float4* vel_in,   // xyz (linear vel)
    __global float4* omega_in, // xyz (angular vel)
    // WRITE-ONLY NEW STATE
    __global float4* pos_out,
    __global float4* quat_out,
    __global float4* vel_out,
    __global float4* omega_out,
    // PARAMS
    __global int4*   neighs,
    __global float4* bLs,       // Rest Lengths
    __global float4* bKs,       // Stiffness
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // 1. PREDICTED STATE OF SELF
    float4 p4_i = pos_in[i];
    float3 p_i  = p4_i.xyz;
    float  w_i  = p4_i.w;  // Inverse Mass
    float4 q_i  = (i < nnode) ? quat_in[i] : (float4)(0,0,0,1);
    float3 v_i  = vel_in[i].xyz;
    float3 om_i = omega_in[i].xyz;

    // Step 1: Apply Inertia (Prediction)
    p_i += v_i * dt;
    if (i < nnode) {
        float4 dq = quat_from_vec(om_i * dt);
        q_i = normalize(quat_mul(dq, q_i));
    }

    // 2. ACCUMULATORS
    float3 sum_dp = (float3)(0.0f);    // Sum of position corrections
    float3 sum_dq_vec = (float3)(0.0f);// Sum of rotation corrections (torque-like)
    float  count = 0.0f;               // Number of constraints handled

    // 3. DEFINE PORTS (Local bond directions for sp2)
    // In a real kernel, load these from constant memory or 'apars'
    float3 local_ports[3];
    local_ports[0] = (float3)( 1.0f,  0.0f, 0.0f);
    local_ports[1] = (float3)(-0.5f,  0.866f, 0.0f);
    local_ports[2] = (float3)(-0.5f, -0.866f, 0.0f);

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;

    // 4. SOLVE CONSTRAINTS
    for (int k = 0; k < 4; k++) { // 4 max neighbors
        int j = neighbors[k];
        if (j < 0) break;

        // READ NEIGHBOR (Prediction)
        float4 p4_j = pos_in[j];
        float3 p_j  = p4_j.xyz + vel_in[j].xyz * dt; // Predict neighbor pos
        float  w_j  = p4_j.w;

        // --- CONSTRAINT TYPE 1: I am a Node, J is a neighbor attached to Port[k] ---
        // This enforces directionality AND distance
        if (i < nnode && k < 3) {
            float3 r_port_local = local_ports[k] * 1.4f; // 1.4A arm length approx, or separate param
            float3 r_arm = quat_rotate(q_i, r_port_local);
            float3 tip_pos = p_i + r_arm; // Position of the bond "port" in world space

            // The Constraint: The Tip of my Port should touch the Neighbor Center
            // (Or simpler: Bond direction aligns with port. Let's do Distance constraint to Tip)
            
            float3 diff = p_j - tip_pos;
            float dist = length(diff);
            float3 n = diff / dist; // Direction from Me to Him

            // XPBD Logic
            // C = dist - 0 (They should overlap? Or use rest length?)
            // If r_arm is the full bond length, dist should be 0.
            // If r_arm is unit vector, dist should be RestLength - 1.0.
            
            // Let's assume r_arm is unit vector (Direction only)
            // C = dist - (RestLength - 1.0) is tricky.
            
            // BETTER: "Rigid Rod" constraint.
            // My Port vector (rotated) should point to J.
            // And the distance should be RestLength.
            
            // Part A: Linear Fix (Standard PBD)
            // Move p_i towards p_j to satisfy distance
            // float rest_len = ((float*)&bLs[i])[k];
            // float stiffness = ((float*)&bKs[i])[k];
            float rest_len = fetch4(bLs, i, k);
            float stiffness = fetch4(bKs, i, k);
            float alpha = 1.0f / (stiffness * dt * dt + 1e-6f);
            
            float3 r_ij = p_j - p_i;
            float d_ij = length(r_ij);
            float3 dir_ij = r_ij / d_ij;
            
            float C_lin = d_ij - rest_len;
            float lambda_lin = -C_lin / (w_i + w_j + alpha);
            
            // Apply Linear Recoil to Self
            sum_dp += dir_ij * lambda_lin * w_i; 

            // Part B: Angular Fix (Torque)
            // Rotate q_i so 'r_arm' aligns with 'dir_ij'
            // My Port: r_arm (unit). Target: dir_ij (unit).
            float3 curr_dir = quat_rotate(q_i, local_ports[k]); // Unit vector
            float3 axis = cross(curr_dir, dir_ij);
            float sin_theta = length(axis);
            
            // Angular Stiffness (Planarization force)
            float alpha_rot = alpha * 0.1f; // Usually stiffer or softer?
            
            // Inverse Inertia (simplified sphere)
            float w_rot = 10.0f; // Approx 1/I
            
            // Lambda for rotation
            // We want to rotate by 'theta'
            // C_ang = theta ~ sin_theta
            // Generalized: d_omega = C / (w_rot + alpha)
            // Correction vector (axis-angle)
            if (sin_theta > 1e-5f) {
                float3 d_om = axis * (1.0f / (w_rot + alpha_rot));
                sum_dq_vec += d_om * w_rot; 
            }
        }
        
        // --- CONSTRAINT TYPE 2: J is a Node, I am attached to J's port ---
        // (Symmetry check: Does J treat me as a port attachment?)
        // In your system, if I am a Hydrogen, I just have a distance constraint.
        // If I am a Node, and J is a Node, we likely have a bond that controls BOTH angles.
        
        // FOR SIMPLICITY: 
        // We only processed "My Angles" above. 
        // We rely on J processing "His Angles" in his thread.
        // We PROCESSED distance above. J will also process distance. 
        // This effectively applies the distance constraint twice (once in my thread, once in his).
        // This is fine! It just means effective stiffness is 2x. We can halve k in params.
        
        count += 1.0f;
    }

    // 5. FINALIZE UPDATE (Jacobi Average)
    if (count > 0.0f) {
        // Average the linear suggestions
        p_i += sum_dp / count; 
        
        // Average the angular suggestions
        if (i < nnode) {
           float3 final_rot_vec = sum_dq_vec / count;
           float4 dq = quat_from_vec(final_rot_vec);
           q_i = normalize(quat_mul(dq, q_i));
        }
    }

    // 6. UPDATE VELOCITIES & WRITE
    float3 v_new = (p_i - pos_in[i].xyz) / dt;
    
    // Angular vel update
    float4 q_old = (i < nnode) ? quat_in[i] : (float4)(0,0,0,1);
    float4 q_diff = quat_mul(q_i, (float4)(-q_old.xyz, q_old.w));
    float3 om_new = q_diff.xyz * (2.0f / dt);
    if (q_diff.w < 0) om_new *= -1.0f;

    pos_out[i]   = (float4)(p_i, w_i);
    vel_out[i]   = (float4)(v_new, 0.0f);
    if (i < nnode) {
        quat_out[i]  = q_i;
        omega_out[i] = (float4)(om_new, 0.0f);
    }
}

__kernel void solve_rigid_bk_symmetry(
    const int natoms,
    const int nnode,
    __global float4* pos_in,
    __global float4* quat_in,
    __global float4* vel_in,
    __global float4* omega_in,
    __global float4* pos_out,
    __global float4* quat_out,
    __global float4* vel_out,
    __global float4* omega_out,
    __global int4*   neighs,
    __global int4*   bkNeighs,
    __global float4* bLs,
    __global float4* bKs,
    const float dt,
    const float k_rot
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float4 p4_i = pos_in[i];
    float3 p_i  = p4_i.xyz;
    float  w_i  = p4_i.w;
    float4 q_i  = (i < nnode) ? quat_in[i] : (float4)(0,0,0,1);

    // Pure relaxation (no inertial prediction) for Jacobi iterations
    float3 p_pred = p_i;
    float4 q_pred = q_i;

    float3 dp = (float3)(0.0f);
    float3 dphi = (float3)(0.0f);
    float wp = 0.0f;
    float wq = 0.0f;

    float3 ports[3];
    ports[0] = (float3)( 1.0f,  0.0f, 0.0f);
    ports[1] = (float3)(-0.5f,  0.866f, 0.0f);
    ports[2] = (float3)(-0.5f, -0.866f, 0.0f);

    int4 ng  = neighs[i];
    int4 bkg = bkNeighs[i];
    int* neighbors = (int*)&ng;
    int* backSlots = (int*)&bkg;

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float4 p4_j = pos_in[j];
        float3 p_j  = p4_j.xyz;
        float  w_j  = p4_j.w;
        float3 pj_pred = p_j;

        float3 r_ij = pj_pred - p_pred;
        float dist = length(r_ij);
        if (dist < 1e-8f) continue;
        float3 dir = r_ij / dist;

        float rest_len = fetch4(bLs, i, k);
        float stiff    = fetch4(bKs, i, k);
        float alpha    = 1.0f / (stiff + 1e-12f);
        float C        = dist - rest_len;
        float lambda   = -C / (w_i + w_j + alpha);
        dp += dir * lambda * w_i;
        wp += 1.0f;

        if (i < nnode && k < 3) {
            float3 port_world = quat_rotate(q_pred, ports[k]);
            float3 axis = cross(port_world, dir);
            dphi += axis;
            wq += 1.0f;

            // Tangential position correction (node-side) so that r_ij aligns with port_world
            // This is the translational recoil partner for the neighbor-side correction below.
            float proj = dot(r_ij, port_world);
            float3 t = r_ij - port_world * proj;
            float t2 = dot(t, t);
            if (t2 > 1e-16f) {
                float alpha_rot = 1.0f / (k_rot + 1e-12f);
                float inv = 1.0f / (w_i + w_j + alpha_rot);
                dp += t * (w_i * inv);
                wp += 1.0f;
            }
        }

        int bk = backSlots[k];
        if (bk >= 0 && bk < 3 && j < nnode) {
            float4 q_j  = quat_in[j];
            float3 port_j_world = quat_rotate(q_j, ports[bk]);

            // Tangential position correction (neighbor-side): move self so that r_ji aligns with node's port
            float3 r_ji = p_pred - pj_pred;
            float proj = dot(r_ji, port_j_world);
            float3 t = r_ji - port_j_world * proj;
            float t2 = dot(t, t);
            if (t2 > 1e-16f) {
                float alpha_rot = 1.0f / (k_rot + 1e-12f);
                float inv = 1.0f / (w_i + w_j + alpha_rot);
                dp += (-t) * (w_i * inv);
                wp += 1.0f;
            }
        }
    }

    const float relax_p = 0.2f;
    const float relax_q = 0.1f;

    // NOTE: Do NOT normalize by wp/wq here. Per-atom normalization breaks symmetry between i and j
    // (different neighbor counts), which causes linear momentum drift.
    if (wp > 0) p_pred += dp * relax_p;
    if (wq > 0 && i < nnode) {
        float4 dq_corr = quat_from_vec(dphi * relax_q);
        q_pred = normalize(quat_mul(dq_corr, q_pred));
    }

    float3 v_new = (p_pred - p_i) / dt;

    float3 om_new = (float3)(0.0f);
    if (i < nnode) {
        float4 q_diff = quat_mul(q_pred, (float4)(-q_i.xyz, q_i.w));
        om_new = q_diff.xyz * (2.0f / dt);
        if (q_diff.w < 0) om_new *= -1.0f;
    }

    pos_out[i]   = (float4)(p_pred, w_i);
    vel_out[i]   = (float4)(v_new, 0.0f);
    if (i < nnode) {
        quat_out[i]  = q_pred;
        omega_out[i] = (float4)(om_new, 0.0f);
    }
}

__kernel void integrate_and_project(
    const int natoms,
    const int nnode,
    __global float4* pos,       // xyz, w=invMass
    __global float4* quat,      // xyzw (Nodes only)
    __global float4* vel,       // xyz (Linear)
    __global float4* omega,     // xyz (Angular)
    __global float4* global_ports, // OUTPUT: World space ports
    __global uchar*  atom_types,   // INPUT: port type per atom (0 none, 1 sp1, 2 sp2, 3 sp3)
    __global const float4* port_dirs,
    __global const uchar*  port_ns,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // 1. Integration (Predictor)
    float4 p4 = pos[i];
    float3 p = p4.xyz;
    float3 v = vel[i].xyz;
    
    // Apply Gravity/External Forces here if needed
    v += (float3)(0.0f, 0.0f, 0.0f) * dt; // Placeholder
    p += v * dt; // Prediction

    // Write back predicted position
    pos[i] = (float4)(p, p4.w);
    vel[i] = (float4)(v, 0.0f);

    __local float4 l_port_dirs[N_PORT_TYPES*4];
    __local uchar  l_port_ns  [N_PORT_TYPES];
    porttab_load(l_port_dirs, l_port_ns, port_dirs, port_ns);

    // 2. Rotation & Projection (Nodes Only)
    if (i < nnode) {
        float4 q = quat[i];
        float3 w = omega[i].xyz;
        
        // Integrate Rotation
        // q_new = q + 0.5 * w * q * dt
        float4 dq = (float4)(w * dt * 0.5f, 0.0f);
        // Quaternion mul logic (simplified for readability)
        float4 q_update = (float4)(
            dq.w*q.x + dq.x*q.w + dq.y*q.z - dq.z*q.y,
            dq.w*q.y - dq.x*q.z + dq.y*q.w + dq.z*q.x,
            dq.w*q.z + dq.x*q.y - dq.y*q.x + dq.z*q.w,
            dq.w*q.w - dq.x*q.x - dq.y*q.y - dq.z*q.z
        );
        q += q_update;
        q = normalize(q); // Important!
        quat[i] = q;

        // PROJECT PORTS
        // We calculate where the ports ARE right now based on predicted P and Q
        // and write them to global memory for the solver to see.
        
        // Standard SP2/SP3 geometry can be hardcoded or loaded
        // Assuming 'port_local_defs' holds 4 float4s per node type, or generic
        
        int typ = (int)atom_types[i];
        for (int k=0; k<4; k++) {
            float3 local_p = porttab_dir(l_port_dirs, typ, k);
            float3 t = 2.0f * cross(q.xyz, local_p);
            float3 rot_p = local_p + q.w * t + cross(q.xyz, t);
            global_ports[i*4 + k] = (float4)(p + rot_p, 0.0f);
        }
    }
}


// =============================================================
// KERNEL 2: SOLVE CONSTRAINTS (GATHER)
// =============================================================
__kernel void solve_ports_xpbd(
    const int natoms,
    const int nnode,
    __global float4* pos,        // Predicted P
    __global float4* quat,       // Predicted Q
    __global float4* global_ports, // The "Visual" Ports
    __global int4*   neighs,     // Neighbor Indices
    __global int4*   bkNeighs,   // Back-Indices (Which port on neighbor connects to me?)
    __global float4* bLs,        // Rest lengths
    __global float4* bKs,        // Stiffness
    __global float4* delta_p,    // OUTPUT: Position correction accumulator
    __global float4* delta_rot,  // OUTPUT: Rotation vector accumulator
    __global uchar*  atom_types, // INPUT: port type per atom
    __global const float4* port_dirs,
    __global const uchar*  port_ns,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float3 p_i = pos[i].xyz;
    float w_i  = pos[i].w;
    
    float3 sum_dp = (float3)(0.0f);
    float3 sum_drot = (float3)(0.0f); // Torque/Rotation axis
    float count = 0.0f;

    int4 ng = neighs[i];
    int4 bk = bkNeighs[i]; // Needed!
    int* neighbors = (int*)&ng;
    int* back_indices = (int*)&bk;

    __local float4 l_port_dirs[N_PORT_TYPES*4];
    __local uchar  l_port_ns  [N_PORT_TYPES];
    porttab_load(l_port_dirs, l_port_ns, port_dirs, port_ns);

    for (int k=0; k<4; k++) {
        int j = neighbors[k];
        if (j < 0) break;
        
        int port_idx_on_neighbor = back_indices[k]; // Which port of J connects to I?
        
        float3 target_pos;
        float3 current_pos;
        bool is_angular_constraint = false;
        bool uses_port_tip = false;
        
        // --- LOGIC SPLIT ---
        
        // SCENARIO 1: I am a Node, and I have a Port [k] for this neighbor.
        // My Port [k] should be at Neighbor's Position.
        int typ_i = (i < nnode) ? (int)atom_types[i] : 0;
        int npi = porttab_n(l_port_ns, typ_i);

        if (i < nnode && k < npi) {
            // My "Hand" is the port
            current_pos = global_ports[i*4 + k].xyz; 
            // The "Wall" I am holding is the neighbor
            target_pos = pos[j].xyz; 
            is_angular_constraint = true;
            uses_port_tip = true;
        }
        else {
            // I am just a particle (or this is a non-directional bond).
            // Use centers.
            current_pos = p_i;
            target_pos = pos[j].xyz;
        }
        
        // SCENARIO 2: The Neighbor is a Node.
        // I should be at the location of His Port.
        // (If both are nodes, this adds a second constraint, which is fine/symmetric)
        int typ_j = (j < nnode) ? (int)atom_types[j] : 0;
        int npj = porttab_n(l_port_ns, typ_j);
        if (j < nnode && port_idx_on_neighbor >= 0 && port_idx_on_neighbor < npj) {
            // Override Target: Instead of his center, I want his Port location.
            // We READ the global buffer.
            target_pos = global_ports[j*4 + port_idx_on_neighbor].xyz;
            uses_port_tip = true;
        }

        // --- SOLVE SPRING/CONSTRAINT ---
        
        float3 diff = target_pos - current_pos;
        float dist = length(diff);
        
        // XPBD Params
        // float rest_len = uses_port_tip ? 0.0f : ((float*)&bLs[i])[k]; // DEBUG
        // float stiffness = ((float*)&bKs[i])[k]; // DEBUG
        float rest_len = uses_port_tip ? 0.0f : fetch4(bLs, i, k);
        float stiffness = fetch4(bKs, i, k);
        float alpha = 1.0f / (stiffness * dt * dt + 1e-6f);

        // If I am a node using a port (Scenario 1), current_pos is offset from my center.
        // This generates TORQUE.
        
        if (is_angular_constraint) {
            // Vector from My Center to My Port
            float3 r = current_pos - p_i; 
            
            // Constraint C = dist - rest
            float C = dist - rest_len;
            
            if (dist > 1e-8f) {
               float3 n = diff / dist; // Direction I need to move/rotate
               
               // XPBD Generalized Inverse Mass
               // w_generalized = w_i + w_j + |cross(r, n)|^2 * w_rot
               // (Simplified for diagonal inertia w_rot)
               float w_rot = 1.0f; // Inverse Inertia (debug-stable)
               float3 r_cross_n = cross(r, n);
               float w_ang = dot(r_cross_n, r_cross_n) * w_rot;
               
               // Lambda
               // We assume J also moves (w_j), but we only apply to I here.
               // For symmetry, we include w_j in denominator.
               float w_j = pos[j].w;
               float lambda = -C / (w_i + w_j + w_ang + alpha);
               
               // Apply to Self
               float3 impulse = n * lambda;
               
               // Linear Update
               sum_dp += impulse * w_i;
               
               // Angular Update: output small rotation vector (axis-angle) for apply_deltas_rigid
               // Approx: dphi ~= dt * I^-1 * (r x impulse)
               sum_drot += cross(r, impulse) * (w_rot * dt);
               
               count += 1.0f;
            }
        } 
        else {
            // Simple Linear Pull (Scenario 2 or plain bond)
            float C = dist - rest_len;
            if (dist > 1e-8f) {
                float3 n = diff / dist;
                float w_j = pos[j].w;
                float lambda = -C / (w_i + w_j + alpha);
                
                sum_dp += n * lambda * w_i;
                count += 1.0f;
            }
        }
    }

    // Write Accumulators
    if (count > 0.0f) {
        // Average the updates (Jacobi style)
        delta_p[i] = (float4)(sum_dp / count, 0.0f);
        delta_rot[i] = (float4)(sum_drot / count, 0.0f);
    } else {
        delta_p[i] = (float4)(0.0f);
        delta_rot[i] = (float4)(0.0f);
    }
}

__kernel void project_ports(
    const int nnode,
    __global float4* pos_pred,      // Predicted position (p + v*dt)
    __global float4* quat_pred,     // Predicted orientation
    __global float4* global_ports,  // OUTPUT: World-space port positions
    __global float4* bLs,           // Rest lengths per slot (arm length)
    __global uchar*  atom_types,    // port type per atom
    __global const float4* port_dirs,
    __global const uchar*  port_ns
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    __local float4 l_port_dirs[N_PORT_TYPES*4];
    __local uchar  l_port_ns  [N_PORT_TYPES];
    porttab_load(l_port_dirs, l_port_ns, port_dirs, port_ns);

    float3 p = pos_pred[i].xyz;
    float4 q = quat_pred[i];

    int typ = (int)atom_types[i];
    // Project all 4 potential ports
    for (int k = 0; k < 4; k++) {
        float3 local_p = porttab_dir(l_port_dirs, typ, k);
        // float L = ((float*)&bLs[i])[k]; // DEBUG
        float L = fetch4(bLs, i, k);

        // 2. Rotate to World Space
        // rot = q * local_p * q_inv
        float3 t = 2.0f * cross(q.xyz, local_p);
        float3 uW = local_p + q.w * t + cross(q.xyz, t);
        float3 r_arm = uW * L;

        // 3. Write "Ideal" Port Position
        // This is where the Neighbor SHOULD be.
        global_ports[i*4 + k] = (float4)(p + r_arm, 0.0f);
    }
}

__kernel void jacobi_solve_rigid(
    const int natoms,
    const int nnode,
    // PREDICTED STATE (Input)
    __global float4* pos_pred,    // xyz, w=invMass
    __global float4* quat_pred,   // xyzw
    // PROJECTED PORTS (Input)
    __global float4* global_ports,
    // OUTPUT STATE
    __global float4* pos_new,
    __global float4* quat_new,
    // TOPOLOGY
    __global int4*   neighs,
    __global int4*   bkNeighs,    // Needed to find "Case 2"
    __global float4* bKs,         // Stiffness
    __global float4* bLs,         // Rest lengths (used as arm length)
    __global uchar*  atom_types,  // port type per atom
    __global const float4* port_dirs,
    __global const uchar*  port_ns,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    __local float4 l_port_dirs[N_PORT_TYPES*4];
    __local uchar  l_port_ns  [N_PORT_TYPES];
    porttab_load(l_port_dirs, l_port_ns, port_dirs, port_ns);

    // --- 1. INITIALIZE WEIGHTED SUMS ---
    // Alpha = M / dt^2 = 1 / (invMass * dt^2)
    float invMass = pos_pred[i].w;
    float alpha = (invMass > 1e-6f) ? (1.0f / (invMass * dt * dt)) : 1e6f;
    
    // Linear Accumulators
    float3 sum_pos_numerator = pos_pred[i].xyz * alpha;
    float  sum_pos_weight    = alpha;

    // Angular Accumulators (Torque-like)
    // For rotation, Projective Dynamics often linearizes, but we can use 
    // a weighted average of "Target Orientations" or "Torque Vectors".
    // Simple approach: Accumulate torque vector to apply to Q later.
    float3 sum_torque = (float3)(0.0f); 
    float  sum_rot_weight = alpha; // Rotational inertia approx

    float4 q_i = (i < nnode) ? quat_pred[i] : (float4)(0,0,0,1);
    float3 p_i = pos_pred[i].xyz;

    int4 ng = neighs[i];
    int4 bk = bkNeighs[i];
    int* neighbors = (int*)&ng;
    int* back_indices = (int*)&bk;

    // --- 2. GATHER CONSTRAINTS ---
    
    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;
        
        // float stiffness = ((float*)&bKs[i])[k]; // K_ij // DEBUG
        // float Lij       = ((float*)&bLs[i])[k]; // bond length used as arm length // DEBUG
        float stiffness = fetch4(bKs, i, k);
        float Lij       = fetch4(bLs, i, k);
        
        float3 p_j = pos_pred[j].xyz;

        // ============================================================
        // CASE 0: Plain distance constraint (always active)
        // ============================================================
        // target for i: p_j - dir * Lij
        float3 r0 = p_j - p_i;
        float d0 = length(r0);
        if (d0 > 1e-8f) {
            float3 dir0 = r0 / d0;
            float3 target_i = p_j - dir0 * Lij;
            sum_pos_numerator += target_i * stiffness;
            sum_pos_weight    += stiffness;
        }

        // ============================================================
        // CASE A: I am a Node, J is connected to My Port[k]
        // ============================================================
        int typ_i = (i < nnode) ? (int)atom_types[i] : 0;
        int npi = porttab_n(l_port_ns, typ_i);
        if (i < nnode && k < npi) {
            // The constraint says: "My Port[k] should be at J's position"
            // To satisfy this, I must move and rotate.
            
            float3 p_j = pos_pred[j].xyz; // Target Position
            
            // 1. Where is my port relative to me?
            // Recompute locally or read from global_ports[i*4+k] - p_i
            // Better to recompute to avoid precision drift from global buffer
            // Use rotated unit port dir scaled by bond length as arm
            float3 local_port = porttab_dir(l_port_dirs, typ_i, k);
            float3 t2 = 2.0f * cross(q_i.xyz, local_port);
            float3 uW = local_port + q_i.w * t2 + cross(q_i.xyz, t2);
            float3 r_arm = uW * Lij;
            float3 port_world_pos = p_i + r_arm;

            // 2. Linear Contribution to Me
            // Ideally, p_i should be at (p_j - r_arm)
            float3 target_pos_for_me = p_j - r_arm;
            
            sum_pos_numerator += target_pos_for_me * stiffness;
            sum_pos_weight    += stiffness;

            // 3. Angular Contribution to Me
            // Ideally, r_arm should point towards (p_j - p_i).
            // This creates a "Target Rotation".
            // Torque = r_arm x Force_direction. 
            // Force_direction ~ (p_j - port_world_pos).
            float3 diff = p_j - port_world_pos; 
            
            // Add torque: r x (K * displacement)
            // This fits the "Momentum" view: Force at arm tip.
            sum_torque += cross(r_arm, diff) * stiffness;
            // Weight? Rotational weight depends on radius^2.
            sum_rot_weight += stiffness * dot(r_arm, r_arm); 
        }

        // ============================================================
        // CASE B: Neighbor J is a Node, I am connected to His Port
        // ============================================================
        // Check if J is a node and if I am attached to a directional port
        int port_idx_on_j = back_indices[k];
        int typ_j = (j < nnode) ? (int)atom_types[j] : 0;
        int npj = porttab_n(l_port_ns, typ_j);
        if (j < nnode && port_idx_on_j >= 0 && port_idx_on_j < npj) {
            // The constraint says: "I should be at J's Port Position"
            
            // Read J's Port from Global Buffer (Crucial!)
            float3 target_pos = global_ports[j*4 + port_idx_on_j].xyz;
            
            // This is a pure linear pull on Me towards that red sphere
            sum_pos_numerator += target_pos * stiffness;
            sum_pos_weight    += stiffness;
            
            // No torque on Me (unless I am also a node trying to align, 
            // covered by Case A of the reciprocal bond)
        }
        
        // If neither is a node (H-H bond?), standard distance logic applies
        // ... (omitted for brevity, similar to Case B but target is p_j)
    }

    // --- 3. FINALIZE (JACOBI UPDATE) ---
    
    // Linear Update
    float3 p_final = sum_pos_numerator / sum_pos_weight;
    
    // Angular Update
    // q_new = Integrate(q_old, angular_velocity_from_torque)
    // Effective "displacement" in angle = sum_torque / sum_rot_weight
    float3 delta_rot = sum_torque / sum_rot_weight;
    
    // Apply small rotation
    float4 dq = (float4)(delta_rot * 0.5f, 0.0f); 
    // Simplified quat add (valid for small angles)
    float4 q_final = q_i + (float4)(
            dq.w*q_i.x + dq.x*q_i.w + dq.y*q_i.z - dq.z*q_i.y,
            dq.w*q_i.y - dq.x*q_i.z + dq.y*q_i.w + dq.z*q_i.x,
            dq.w*q_i.z + dq.x*q_i.y - dq.y*q_i.x + dq.z*q_i.w,
            dq.w*q_i.w - dq.x*q_i.x - dq.y*q_i.y - dq.z*q_i.z
    );
    q_final = normalize(q_final);

    // Write Output
    pos_new[i]  = (float4)(p_final, invMass);
    if (i < nnode) quat_new[i] = q_final;
}

__kernel void apply_deltas(
    const int natoms,
    __global float4* pos_pred,    // xyz, w=invMass
    __global float4* quat_pred,   // xyzw
    __global float4* delta_p,     // Linear impulse
    __global float4* delta_rot,   // Angular impulse
    __global float4* pos_new,     // New position
    __global float4* quat_new,    // New orientation
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float3 p_i = pos_pred[i].xyz;
    float4 q_i = quat_pred[i];

    float3 dp = delta_p[i].xyz;
    float4 dq = delta_rot[i];

    // Linear Update
    float3 p_final = p_i + dp;

    // Angular Update
    // q_new = Integrate(q_old, angular_velocity_from_torque)
    // Effective "displacement" in angle = sum_torque / sum_rot_weight
    float4 q_final = q_i + (float4)(
            dq.w*q_i.x + dq.x*q_i.w + dq.y*q_i.z - dq.z*q_i.y,
            dq.w*q_i.y - dq.x*q_i.z + dq.y*q_i.w + dq.z*q_i.x,
            dq.w*q_i.z + dq.x*q_i.y - dq.y*q_i.x + dq.z*q_i.w,
            dq.w*q_i.w - dq.x*q_i.x - dq.y*q_i.y - dq.z*q_i.z
    );
    q_final = normalize(q_final);

    // Write Output
    pos_new[i]  = (float4)(p_final, pos_pred[i].w);
    quat_new[i] = q_final;
}

// ===============================
//  Solution from Kimi 2.5
// ===============================


__kernel void compute_corrections(
    const int nnode,
    __global const float4* pos,        // w = invMass
    __global const float4* quat,
    __global const int4*   neighs,
    __global const float4* port_local, // Local offsets (xyz, 0)
    __global const float*  stiffness,  // K weight
    __global float4* dpos_node,        // Accumulated COM correction for node i
    __global float4* drot_node,        // Accumulated rotation vector for node i
    __global float4* dpos_neigh        // Correction to write to neighbor's slot (recoils)
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float3 xi = pos[i].xyz;
    float4 qi = quat[i];
    float invMi = pos[i].w;
    float mi = (invMi > 1e-12f) ? 1.0f/invMi : 1e12f;
    float3 invIi_vec = (float3)(2.5f * invMi); // Simplified: 1/(0.4*m) for sphere

    float3 dx_accum = (float3)(0.0f);
    float3 dtheta_accum = (float3)(0.0f);
    int i4 = i * 4;

    for (int k = 0; k < 4; k++) {
        int j = neighs[i][k];
        if (j < 0) break;

        float3 aj = port_local[i4 + k].xyz;
        float3 ri = quat_rotate(qi, aj); // Arm in world space
        float3 pi = xi + ri;             // Port position
        
        float3 xj = pos[j].xyz;
        float invMj = pos[j].w;
        
        float3 diff = pi - xj;  // tip - neighbor
        float dist2 = dot(diff, diff);
        if (dist2 < 1e-12f) continue;
        float dist = sqrt(dist2);
        float3 n = diff / dist;
        float K = stiffness[i4 + k]; // weight
        
        // Generalized inverse masses
        float W_lin_i = invMi;
        float W_lin_j = invMj;
        
        float3 rxn = cross(ri, n);
        // W_rot = (r x n)^T * I^-1 * (r x n). Using diagonal approx:
        float W_rot_i = dot(rxn*rxn, invIi_vec); 
        
        float W_total = W_lin_i + W_lin_j + W_rot_i + 1e-12f;
        float delta_lambda = -dist / W_total;
        
        // --- Correction for node i (accumulated) ---
        // Translational part
        float3 dx_i = delta_lambda * n * W_lin_i;
        dx_accum += dx_i;
        
        // Rotational part: dtheta = I^-1 * (r x n) * scale
        float3 dtheta_i = (invIi_vec * rxn) * (W_rot_i / W_total) * delta_lambda * dist;
        dtheta_accum += dtheta_i;
        
        // --- Correction for neighbor j (recoil) ---
        // This is stored in the slot so j can pick it up
        float3 dx_j = -delta_lambda * n * W_lin_j;
        dpos_neigh[i4 + k] = (float4)(dx_j, 0.0f);
    }
    
    dpos_node[i] = (float4)(dx_accum, 0.0f);
    drot_node[i] = (float4)(dtheta_accum, 0.0f);
}

__kernel void apply_corrections(
    const int natoms,
    const int nnode,
    __global float4* pos,
    __global float4* quat,
    __global const int4* bkSlots,      // Back-mapping: which slots affect me?
    __global const float4* dpos_node,  // From kernel 1 (nodes only)
    __global const float4* drot_node,  // From kernel 1 (nodes only)
    __global const float4* dpos_neigh, // From kernel 1 (all bonds)
    const float relaxation             // 0.2 - 0.8 for Jacobi stability
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float3 dx = (float3)(0.0f);
    if (i < nnode) {
        dx = dpos_node[i].xyz;
    }
    
    // Gather recoils from all neighbors that wrote to my slots
    int4 bk = bkSlots[i];
    for (int k = 0; k < 4; k++) {
        int slot = bk[k];
        if (slot >= 0) {
            dx += dpos_neigh[slot].xyz;
        }
    }
    
    // Apply position update (Jacobi step with under-relaxation)
    float3 xi = pos[i].xyz;
    xi += dx * relaxation;
    pos[i].xyz = xi;
    
    // Apply rotation update (only nodes)
    if (i < nnode) {
        float3 dtheta = drot_node[i].xyz * relaxation;
        
        float angle = length(dtheta);
        if (angle > 1e-8f) {
            float3 axis = dtheta / angle;
            float4 dq = (float4)(axis * sin(angle * 0.5f), cos(angle * 0.5f));
            float4 qi = quat[i];
            // Apply rotation: q_new = dq * q_old (local rotation applied to current)
            quat[i] = normalize(quat_mul(dq, qi));
        }
    }
}

// ===============================
//  Solution from Claude 4.5 Sonet 
// ===============================

__kernel void compute_and_apply_corrections(
    const int nnode,
    __global const float4* pos,          // w = invMass
    __global const float4* quat,
    __global const int4*   neighs,
    __global const float4* bKs,
    __global const float4* port_local,
    __global float4*       delta_neigh,  // size nnode*4
    __global float4*       pos_delta,    // size nnode
    __global float4*       omega_delta,  // size nnode
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float4 pi4 = pos[i];
    float invMi = pi4.w;
    float mi = (invMi > 1e-12f) ? (1.0f / invMi) : 1e12f;
    float3 p_i = pi4.xyz;
    float4 q_i = quat[i];

    float I_i = 0.4f * mi;
    (void)dt; // placeholder to match signature; dt may be used for compliance later

    float3 dp_i_total = (float3)(0.0f);
    float3 dw_i_total = (float3)(0.0f);

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float K_ij = fetch4(bKs, i, k);
        float4 pj4 = pos[j];
        float invMj = pj4.w;
        float mj = (invMj > 1e-12f) ? (1.0f / invMj) : 1e12f;

        float3 r_k = quat_rotate(q_i, port_local[i * 4 + k].xyz);
        float3 port_world = p_i + r_k;

        float3 v_k = pj4.xyz - port_world;
        float C = length(v_k);
        if (C < 1e-8f) continue;
        float3 n = v_k / C;

        float3 dr_dq = cross(r_k, n);
        float w_i_trans = invMi;
        float w_i_rot = (I_i > 1e-12f) ? dot(dr_dq, dr_dq) / I_i : 0.0f;
        float w_j = invMj;

        float w_total = w_i_trans + w_i_rot + w_j + 1.0f / (K_ij + 1e-12f);
        float lambda = -C / w_total;

        float3 dp_i = lambda * n * w_i_trans;
        float3 dw_i = (I_i > 1e-12f) ? lambda * cross(r_k, n) / I_i : (float3)(0.0f);
        float3 dp_j = -lambda * n * w_j;

        dp_i_total += dp_i;
        dw_i_total += dw_i;

        delta_neigh[i * 4 + k] = (float4)(dp_j, 0.0f);
    }

    pos_delta[i] = (float4)(dp_i_total, 0.0f);
    omega_delta[i] = (float4)(dw_i_total, 0.0f);
}

__kernel void gather_and_apply(
    const int natoms,
    const int nnode,
    __global float4*       pos,
    __global float4*       quat,
    __global const int4*   bkSlots,
    __global const float4* delta_neigh,
    __global const float4* pos_delta,
    __global const float4* omega_delta
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float4 pi4 = pos[i];
    float3 dp = pos_delta[i].xyz;

    int4 bk = bkSlots[i];
    int* islots = (int*)&bk;
    for (int k = 0; k < 4; k++) {
        int islot = islots[k];
        if (islot >= 0) {
            dp += delta_neigh[islot].xyz;
        }
    }

    pos[i] = (float4)(pi4.xyz + dp, pi4.w);

    if (i < nnode) {
        float3 dw = omega_delta[i].xyz;
        float angle = length(dw);
        if (angle > 1e-8f) {
            float3 axis = dw / angle;
            float4 dq = quat_from_axis_angle(axis, angle);
            quat[i] = normalize(quat_mul(dq, quat[i]));
        }
    }
}

// ================================================================
//  Solution from Claude 4.5 Sonet XPBD with λ accumulation
// ==============================================================


__kernel void reset_lambda(
    const int nnode,
    __global float* lambda
) {
    int i = get_global_id(0);
    if (i >= nnode * 4) return;
    lambda[i] = 0.0f;
}

__kernel void compute_and_apply_corrections_xpbd(
    const int nnode,
    __global const float4* pos,          // w = invMass
    __global const float4* quat,
    __global const int4*   neighs,
    __global const float4* bKs,          // now interpreted as stiffness K (compliance α = 1/K)
    __global const float4* port_local,
    __global float*        lambda,       // size nnode*4, accumulated per constraint
    __global float4*       delta_neigh,  // size nnode*4
    __global float4*       pos_delta,    // size nnode
    __global float4*       omega_delta,  // size nnode
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float4 pi4 = pos[i];
    float invMi = pi4.w;
    float mi = (invMi > 1e-12f) ? (1.0f / invMi) : 1e12f;
    float3 p_i = pi4.xyz;
    float4 q_i = quat[i];

    float I_i = 0.4f * mi;
    float dt2 = dt * dt;

    float3 dp_i_total = (float3)(0.0f);
    float3 dw_i_total = (float3)(0.0f);

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float K_ij = fetch4(bKs, i, k);
        if (K_ij < 1e-12f) continue;  // skip if no stiffness
        
        float4 pj4 = pos[j];
        float invMj = pj4.w;
        float mj = (invMj > 1e-12f) ? (1.0f / invMj) : 1e12f;

        // Port position in world frame
        float3 r_k = quat_rotate(q_i, port_local[i * 4 + k].xyz);
        float3 port_world = p_i + r_k;

        // Constraint violation (consistent with force kernels): diff = tip - xj
        float3 v_k = port_world - pj4.xyz;
        float C = length(v_k);
        if (C < 1e-8f) continue;
        float3 n = v_k / C;

        // Compute generalized mass (accounts for translation + rotation coupling)
        float3 dr_dq = cross(r_k, n);
        float w_i_trans = invMi;
        float w_i_rot = (I_i > 1e-12f) ? dot(dr_dq, dr_dq) / I_i : 0.0f;
        float w_j = invMj;

        // XPBD: add compliance term α/dt² where α = 1/K
        float alpha = 1.0f / K_ij;  // compliance
        float alpha_tilde = alpha / dt2;
        
        float w_total = w_i_trans + w_i_rot + w_j + alpha_tilde;

        int constraint_idx = i * 4 + k;
        float lambda_prev = lambda[constraint_idx];
        float delta_lambda = (-C - alpha_tilde * lambda_prev) / w_total;
        lambda[constraint_idx] = lambda_prev + delta_lambda;

        // Apply corrections using Δλ (not total λ!)
        float3 dp_i = delta_lambda * n * w_i_trans;
        float3 dw_i = (I_i > 1e-12f) ? delta_lambda * cross(r_k, n) / I_i : (float3)(0.0f);
        float3 dp_j = -delta_lambda * n * w_j;

        dp_i_total += dp_i;
        dw_i_total += dw_i;

        delta_neigh[i * 4 + k] = (float4)(dp_j, 0.0f);
    }

    pos_delta[i] = (float4)(dp_i_total, 0.0f);
    omega_delta[i] = (float4)(dw_i_total, 0.0f);
}

__kernel void gather_and_apply_xpbd(
    const int natoms,
    const int nnode,
    __global float4*       pos,
    __global float4*       quat,
    __global const int4*   bkSlots,
    __global const float4* delta_neigh,
    __global const float4* pos_delta,
    __global const float4* omega_delta
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float4 pi4 = pos[i];
    float3 dp = (float3)(0.0f);
    
    // Get delta for this atom (0 if i >= nnode, since pos_delta only allocated for nnode)
    if (i < nnode) {
        dp = pos_delta[i].xyz;
    }

    // Gather recoil corrections from neighbors
    int4 bk = bkSlots[i];
    int* islots = (int*)&bk;
    for (int k = 0; k < 4; k++) {
        int islot = islots[k];
        if (islot >= 0) {
            dp += delta_neigh[islot].xyz;
        }
    }

    pos[i] = (float4)(pi4.xyz + dp, pi4.w);

    if (i < nnode) {
        float3 dw = omega_delta[i].xyz;
        float angle = length(dw);
        if (angle > 1e-8f) {
            float3 axis = dw / angle;
            float4 dq = quat_from_axis_angle(axis, angle);
            quat[i] = normalize(quat_mul(dq, quat[i]));
        }
    }
}


__kernel void compute_and_apply_corrections_xpbd_vector(
    const int nnode,
    __global const float4* pos,
    __global const float4* quat,
    __global const int4*   neighs,
    __global const float4* bKs,
    __global const float4* port_local,
    __global float*        lambda,
    __global float4*       delta_neigh,
    __global float4*       pos_delta,
    __global float4*       omega_delta,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float4 pi4 = pos[i];
    float invMi = pi4.w;
    float mi = (invMi > 1e-12f) ? (1.0f / invMi) : 1e12f;
    float3 p_i = pi4.xyz;
    float4 q_i = quat[i];

    float I_i = 0.4f * mi;
    float dt2 = dt * dt;

    float3 dp_i_total = (float3)(0.0f);
    float3 dw_i_total = (float3)(0.0f);

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) break;

        float K_ij = fetch4(bKs, i, k);
        if (K_ij < 1e-12f) continue;
        
        float4 pj4 = pos[j];
        float invMj = pj4.w;
        float mj = (invMj > 1e-12f) ? (1.0f / invMj) : 1e12f;

        // Port position in world frame
        float3 r_k = quat_rotate(q_i, port_local[i * 4 + k].xyz);
        float3 port_world = p_i + r_k;

        // Constraint violation (consistent with force kernels): diff = tip - xj
        float3 diff = port_world - pj4.xyz;
        float diff2 = dot(diff, diff);
        if (diff2 < 1e-16f) continue;

        // For vector form, we need to compute effective mass along diff direction
        // This is more complex for rotations, so we'll use a scalar constraint
        float C = sqrt(diff2);
        float3 n = diff / C;

        float3 dr_dq = cross(r_k, n);
        float w_i_trans = invMi;
        float w_i_rot = (I_i > 1e-12f) ? dot(dr_dq, dr_dq) / I_i : 0.0f;
        float w_j = invMj;

        float alpha = 1.0f / K_ij;
        float alpha_tilde = alpha / dt2;
        
        int constraint_idx = i * 4 + k;
        float lambda_prev = lambda[constraint_idx];
        
        float denom = w_i_trans + w_i_rot + w_j + alpha_tilde + 1e-12f;
        float delta_lambda = (-C - alpha_tilde * lambda_prev) / denom;
        
        lambda[constraint_idx] += delta_lambda;

        // Apply using vector diff (equivalent to delta_lambda * n since n = diff/C)
        float3 dp_i = delta_lambda * n * w_i_trans;
        float3 dw_i = (I_i > 1e-12f) ? delta_lambda * cross(r_k, n) / I_i : (float3)(0.0f);
        float3 dp_j = -delta_lambda * n * w_j;

        dp_i_total += dp_i;
        dw_i_total += dw_i;

        delta_neigh[i * 4 + k] = (float4)(dp_j, 0.0f);
    }

    pos_delta[i] = (float4)(dp_i_total, 0.0f);
    omega_delta[i] = (float4)(dw_i_total, 0.0f);
}