// XPBD_2D.cl - 2D Position Based Dynamics with Complex Number Rotation
// In 2D, rotation is represented by complex numbers z = (cosθ, sinθ)
// Rotating a vector v by z: z * v = (z.x*v.x - z.y*v.y, z.y*v.x + z.x*v.y)

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define MAX_DEGREE     4      // Max ports per node
 
 // Cluster collision solver configuration (mirrors XPDB_new.cl)
 #define GROUP_SIZE     64      // Workgroup size = atoms per cluster
 #define MAX_GHOSTS     128     // Max external atoms per cluster

// ------------------------------------------------------------------
// 2D COMPLEX NUMBER HELPERS (rotation representation)
// ------------------------------------------------------------------

// Rotate vector v by complex rotation z = (cosθ, sinθ)
inline float2 cmplx_rotate(float2 z, float2 v) {
    return (float2)(
        z.x * v.x - z.y * v.y,  // real part
        z.y * v.x + z.x * v.y   // imag part
    );
}


__kernel void compute_cluster_fused_2d(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const float*  radius,
    __global const float2* mass,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const int*    cluster_counts,
    __global const float2* port_local,
    __global const float*  stiffness_flat,
    __global float2*       dpos_node,
    __global float*        dtheta_node,
    __global float2*       dpos_neigh,
    const int num_atoms,
    const float dt,
    const float k_coll,
    const int group_stride
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * group_stride + lid;

    __local float2 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float2 l_rot[GROUP_SIZE + MAX_GHOSTS];
    __local float2 l_mass[GROUP_SIZE + MAX_GHOSTS];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];
    __local float2 l_pos_new[GROUP_SIZE];

    if (my_global_id < num_atoms) {
        l_pos[lid]  = pos[my_global_id];
        l_rot[lid]  = rot[my_global_id];
        l_mass[lid] = mass[my_global_id];
        l_rad[lid]  = radius[my_global_id];
    } else {
        l_pos[lid]  = (float2)(0.0f, 0.0f);
        l_rot[lid]  = (float2)(1.0f, 0.0f);
        l_mass[lid] = (float2)(1.0f, 1.0f);
        l_rad[lid]  = 0.0f;
    }

    int g_count  = ghost_counts[grp];
    int c_count  = cluster_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int base = grp * MAX_GHOSTS;
        int gid = ghost_indices_flat[base + k];
        int lj = c_count + k;
        l_pos[lj]  = pos[gid];
        l_rot[lj]  = rot[gid];
        l_mass[lj] = mass[gid];
        l_rad[lj]  = radius[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        float my_mass = l_mass[lid].x;
        float alpha = my_mass / (dt*dt + 1e-20f);
        float my_rad = l_rad[lid];

        float2 p = l_pos[lid];
        float2 rhs = p * alpha;
        float  k_sum = alpha;

        int4 ng = neighs_local[my_global_id];
        int* neigh = (int*)&ng;
        int n_ext = c_count + g_count;
        for (int j = 0; j < n_ext; j++) {
            if (j == lid) continue;
            if (j >= c_count && j < group_stride) continue;
            bool is_bonded = false;
            for (int kk = 0; kk < MAX_DEGREE; kk++) {
                int nb = neigh[kk];
                if (nb < 0) continue;
                if (nb == j) { is_bonded = true; break; }
            }
            if (is_bonded) continue;
            float2 q = l_pos[j];
            float2 d = p - q;
            float dist_sq = dot(d, d);
            float r_sum = my_rad + l_rad[j];
            float r2 = r_sum * r_sum;
            if (dist_sq < r2 && dist_sq > 1e-16f) {
                float dist = sqrt(dist_sq);
                float coeff = k_coll * (r_sum / dist);
                rhs += d * coeff;
                rhs += q * k_coll;
                k_sum += k_coll;
            }
        }
        l_pos_new[lid] = rhs / (k_sum + 1e-20f);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        float2 p0 = l_pos[lid];
        float2 p1 = l_pos_new[lid];
        float2 dpos = p1 - p0;
        dpos_node[my_global_id] = dpos;
    }

    if (my_global_id < nnode) {
        float2 p_i = l_pos_new[lid];
        float2 z_i = l_rot[lid];

        float M_i = l_mass[lid].x;
        float I_i = l_mass[lid].y;
        float w_i     = (M_i > 1e-9f) ? 1.0f / M_i : 0.0f;
        float w_rot_i = (I_i > 1e-9f) ? 1.0f / I_i : 0.0f;
        float dt2 = dt * dt + 1e-16f;

        float2 sum_dpos = (float2)(0.0f, 0.0f);
        float  sum_dtheta = 0.0f;
        int4 ng = neighs_local[my_global_id];
        int* neighbors = (int*)&ng;
        int i4 = my_global_id * MAX_DEGREE;

        for (int k = 0; k < MAX_DEGREE; k++) {
            int idx = i4 + k;
            dpos_neigh[idx] = (float2)(0.0f, 0.0f);
            int jloc = neighbors[k];
            if (jloc < 0) continue;
            if (jloc >= c_count && jloc < group_stride) continue;
            float K = stiffness_flat[idx];
            if (K <= 0.0f) continue;

            float2 r_local = port_local[idx];
            float2 r_arm   = cmplx_rotate(z_i, r_local);
            float2 tip     = p_i + r_arm;
            float2 p_j     = l_pos[jloc];

            float2 diff = p_j - tip;
            float dist  = length(diff);
            if (dist < 1e-9f) continue;
            float2 n = diff / dist;

            float alpha = 1.0f / (K * dt2);
            float M_j = l_mass[jloc].x;
            float w_j = (M_j > 1e-9f) ? 1.0f / M_j : 0.0f;

            float rn        = r_arm.x * n.y - r_arm.y * n.x;
            float w_angular = rn * rn * w_rot_i;
            float w_total   = w_i + w_j + w_angular + alpha;

            float impulse_mag = dist / w_total;
            impulse_mag *= 0.5f;

            float2 P = n * impulse_mag;
            sum_dpos += P * w_i;
            sum_dtheta += (rn * impulse_mag) * w_rot_i;
            dpos_neigh[idx] = -P * w_j;
        }

        dpos_node[my_global_id] += sum_dpos;
        dtheta_node[my_global_id] = sum_dtheta;
    }
}

// Complex multiplication: z1 * z2
inline float2 cmplx_mul(float2 z1, float2 z2) {
    return (float2)(
        z1.x * z2.x - z1.y * z2.y,
        z1.y * z2.x + z1.x * z2.y
    );
}

// Get complex rotation from angle (radians)
inline float2 cmplx_from_angle(float angle) {
    return (float2)(cos(angle), sin(angle));
}

// Get angle from complex rotation
inline float cmplx_angle(float2 z) {
    return atan2(z.y, z.x);
}

inline float2 cmplx_normalize(float2 z) {
    float r2 = z.x*z.x + z.y*z.y;
    if (r2 < 1e-20f) return (float2)(1.0f, 0.0f);
    float invr = rsqrt(r2);
    return z * invr;
}

// ------------------------------------------------------------------
// CLUSTER COLLISION HELPERS
// ------------------------------------------------------------------

inline float4 bbox_init_min_2d(float2 p, float r){ return (float4)(p.x - r, p.y - r, 0.0f, 0.0f); }
inline float4 bbox_init_max_2d(float2 p, float r){ return (float4)(p.x + r, p.y + r, 0.0f, 0.0f); }

// ------------------------------------------------------------------
// KERNEL: Initialize heavy-ball buffers
// ------------------------------------------------------------------
__kernel void init_hb_pos_2d(
    const int natoms,
    __global const float2* pos,
    __global float2* hb_pos
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    hb_pos[i] = pos[i];
}

__kernel void init_hb_rot_2d(
    const int nnode,
    __global const float2* rot,
    __global float2* hb_rot
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    hb_rot[i] = cmplx_normalize(rot[i]);
}

// ------------------------------------------------------------------
// KERNEL: Clear forces and node buffers
// ------------------------------------------------------------------
__kernel void clear_2d_forces(
    const int natoms,
    __global float2* force
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    force[i] = (float2)(0.0f, 0.0f);
}

__kernel void clear_2d_node_buffers(
    const int nnode,
    __global float2* fneigh,   // size nnode*4
    __global float2* pneigh    // size nnode*4 (lever arms)
) {
    int i = get_global_id(0);
    int n = nnode * MAX_DEGREE;
    if (i >= n) return;
    fneigh[i] = (float2)(0.0f, 0.0f);
    pneigh[i] = (float2)(0.0f, 0.0f);
}

// ==================================================================================
// ==================================================================================
//    As Rigid As possible (ARAP) https://igl.ethz.ch/projects/ARAP/arap_web.pdf
// ==================================================================================
// ================================================================

// ------------------------------------------------------------------
// KERNEL: Explicit Rigid Shape Matching (2D)
// Solves optimal Translation (t) and Rotation (z) analytically.
// ------------------------------------------------------------------
__kernel void solve_shape_matching_2d(
    const int nnode,
    __global const float2* pos,          // Current world positions
    __global const int4*   neighs,       // Neighbor indices
    __global const float2* port_local,   // Local offsets (unrotated)
    __global const float*  stiffness,    // Weights
    __global float2*       target_pos,   // OUTPUT: Optimal Center Position
    __global float2*       target_rot    // OUTPUT: Optimal Rotation (Complex Unit)
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    // 1. Calculate Weighted Centroids
    // ------------------------------------
    float2 c_local = (float2)(0.0f); // Center of mass of ports
    float2 c_world = (float2)(0.0f); // Center of mass of neighbors
    float sum_w = 0.0f;

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    int i4 = i * 4; // Assuming max 4 neighbors

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) continue;

        float w = stiffness[i4 + k];
        if (w <= 1e-9f) continue;

        c_local += port_local[i4 + k] * w;
        c_world += pos[j] * w;
        sum_w   += w;
    }

    if (sum_w < 1e-9f) return; // Degenerate
    float inv_sum_w = 1.0f / sum_w;
    c_local *= inv_sum_w;
    c_world *= inv_sum_w;

    // 2. Calculate Covariance (Optimal Rotation)
    // ------------------------------------
    // We want to maximize: Real( Sum( w * n_centered * conj(p_centered) ) )
    // effectively calculating the "Average Rotation" required.
    
    float2 cov_sum = (float2)(0.0f);

    for (int k = 0; k < 4; k++) {
        int j = neighbors[k];
        if (j < 0) continue;
        float w = stiffness[i4 + k];
        
        // Centered vectors
        float2 p = port_local[i4 + k] - c_local; // Local arm
        float2 n = pos[j] - c_world;             // World arm
        
        // Complex Multiply: n * conj(p)
        // (nx + i*ny) * (px - i*py) = (nx*px + ny*py) + i(ny*px - nx*py)
        float dot_prod = dot(n, p);              // Real part (Alignment)
        float cross_prod = n.y * p.x - n.x * p.y; // Imag part (Torque)
        
        cov_sum.x += w * dot_prod;
        cov_sum.y += w * cross_prod;
    }

    // 3. Extract Rotation & Translation
    // ------------------------------------
    float2 optimal_rot;
    float norm = length(cov_sum);
    
    if (norm > 1e-12f) {
        optimal_rot = cov_sum / norm;
    } else {
        optimal_rot = (float2)(1.0f, 0.0f); // Identity fallback
    }

    // optimal_pos = c_world - R * c_local
    // Complex multiplication for rotation
    float2 rot_c_local;
    rot_c_local.x = optimal_rot.x * c_local.x - optimal_rot.y * c_local.y;
    rot_c_local.y = optimal_rot.x * c_local.y + optimal_rot.y * c_local.x;

    float2 optimal_pos = c_world - rot_c_local;

    // Write output
    target_pos[i] = optimal_pos;
    target_rot[i] = optimal_rot;
}

// ================================================================
// ================================================================
//   Explicit Force-based Molecular Dynamics (MD)
// ================================================================
// ================================================================

// ------------------------------------------------------------------
// KERNEL: Gather port forces (2D version)
// Each node gathers forces from its ports connected to neighbors
// ------------------------------------------------------------------
__kernel void gather_port_forces_2d(
    const int nnode,
    __global const float2* pos,           // w component stores invMass in .x (reuse)
    __global const float2* rot,           // complex rotation per node
    __global const int4*   neighs,        // neighbor indices (up to 4)
    __global const float4* stiffness,     // K per port (only .x used)
    __global const float2* port_local,    // local port offsets [nnode*4]
    __global const uchar*  port_n,        // number of active ports per node
    __global float2*       force,         // output force per atom
    __global float2*       fneigh,        // recoil force per port slot
    __global float2*       pneigh         // lever arm per port slot
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    int npi = (int)port_n[i];

    float2 fi = (float2)(0.0f, 0.0f);
    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        if (j < 0) continue;

        float Kij = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (Kij <= 0.0f) continue;

        float2 p_j = pos[j];

        // Rotate local port offset to global frame
        float2 r_arm = (float2)(0.0f, 0.0f);
        if (k < npi) {
            r_arm = cmplx_rotate(z_i, port_local[i4 + k]);
        }

        float2 tip = p_i + r_arm;
        float2 f = (p_j - tip) * Kij;

        pneigh[i4 + k] = r_arm;
        fi += f;
        fneigh[i4 + k] = -f;
    }

    force[i] = fi;
}

// ------------------------------------------------------------------
// KERNEL: Explicit integration with force gather (2D)
// ------------------------------------------------------------------
__kernel void integrate_2d_explicit(
    const int natoms,
    const int nnode,
    __global float2*       pos,        // .x = invMass (stored separately for capping atoms)
    __global float2*       vel,
    __global float2*       rot,        // complex rotation
    __global float*        omega,      // angular velocity (scalar in 2D)
    __global const int4*   bkSlots,    // back-slot indices for recoil gather
    __global float2*       force,
    __global float2*       fneigh,     // recoil forces from nodes
    __global float2*       pneigh,     // lever arms
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // For simplicity: uniform mass. invMass passed via separate buffer or constant
    float invMass = 1.0f;  // Could be passed as additional buffer
    float mass = 1.0f;
    float invI = 1.0f / (0.4f * mass);  // Simplified inertia for 2D

    float2 f = force[i];

    // Gather recoil forces from connected nodes
    int4 bk = bkSlots[i];
    if (bk.x >= 0) f += fneigh[bk.x];
    if (bk.y >= 0) f += fneigh[bk.y];
    if (bk.z >= 0) f += fneigh[bk.z];
    if (bk.w >= 0) f += fneigh[bk.w];

    // Linear integration
    float2 v = vel[i];
    v = v * damp + f * invMass * dt;
    pos[i] = pos[i] + v * dt;
    vel[i] = v;

    // Rotational integration (only for nodes)
    if (i < nnode) {
        float2 z = rot[i];

        // Compute torque from recoil forces
        float tau = 0.0f;
        int i4 = i * MAX_DEGREE;
        for (int k = 0; k < MAX_DEGREE; k++) {
            float2 r_arm = pneigh[i4 + k];
            float2 Fk = -fneigh[i4 + k];
            // 2D cross product: tau = r.x * F.y - r.y * F.x
            tau += r_arm.x * Fk.y - r_arm.y * Fk.x;
        }

        // Integrate angular velocity
        float w = omega[i];
        w = w * damp + tau * invI * dt;
        omega[i] = w;

        // Update rotation: integrate angle then convert to complex
        float angle = w * dt;
        float2 dz = cmplx_from_angle(angle);
        rot[i] = cmplx_mul(dz, z);
    }
}

// ------------------------------------------------------------------
// KERNEL: Compute real velocities from position changes (for MD mode)
// v = (x_new - x_old) / dt
// ------------------------------------------------------------------
__kernel void compute_velocities_from_positions(
    const int natoms,
    __global const float2* pos_new,
    __global const float2* pos_old,
    __global float2* vel,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    
    float2 dx = pos_new[i] - pos_old[i];
    float inv_dt = 1.0f / (dt + 1e-20f);
    vel[i] = dx * inv_dt * damp;
}

// ------------------------------------------------------------------
// KERNEL: Compute angular velocities from rotation changes (for MD mode)
// omega = (angle_new - angle_old) / dt
// ------------------------------------------------------------------
__kernel void compute_angular_velocities_from_rotations(
    const int nnode,
    __global const float2* rot_new,
    __global const float2* rot_old,
    __global float* omega,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    
    float2 z_new = rot_new[i];
    float2 z_old = rot_old[i];
    
    // Compute relative rotation: z_rel = conj(z_old) * z_new
    // In 2D: conj(z) = (z.x, -z.y)
    float2 z_old_conj = (float2)(z_old.x, -z_old.y);
    float2 z_rel = cmplx_mul(z_old_conj, z_new);
    
    // Extract angle from z_rel = (cos(dtheta), sin(dtheta))
    float dtheta = atan2(z_rel.y, z_rel.x);
    
    float inv_dt = 1.0f / (dt + 1e-20f);
    omega[i] = dtheta * inv_dt * damp;
}


// ================================================================================
// ================================================================================
//      Collision using Position Based Dynamics (PBD) / Projective Dynamics (PD)
// ================================================================
// ================================================================


// ------------------------------------------------------------------
// CLUSTER COLLISION KERNELS (2D)
// ------------------------------------------------------------------

__kernel void update_bboxes_2d(
    __global const float2* curr_pos,
    __global const float*  radius,
    __global float4*       bboxes_min,
    __global float4*       bboxes_max,
    __local  float4*       local_min,
    __local  float4*       local_max,
    const int num_atoms
) {
    int lid = get_local_id(0);
    int gid = get_global_id(0);
    int grp = get_group_id(0);

    float2 p = (gid < num_atoms) ? curr_pos[gid] : (float2)(0.0f, 0.0f);
    float  r = (gid < num_atoms) ? radius[gid]   : 0.0f;

    if (gid < num_atoms) {
        local_min[lid] = bbox_init_min_2d(p, r);
        local_max[lid] = bbox_init_max_2d(p, r);
    } else {
        local_min[lid] = (float4)(1e10f);
        local_max[lid] = (float4)(-1e10f);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int stride = GROUP_SIZE / 2; stride > 0; stride >>= 1) {
        if (lid < stride) {
            local_min[lid] = min(local_min[lid], local_min[lid + stride]);
            local_max[lid] = max(local_max[lid], local_max[lid + stride]);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if (lid == 0) {
        bboxes_min[grp] = local_min[0];
        bboxes_max[grp] = local_max[0];
    }
}

__kernel void build_local_topology_2d(
    __global const float2* curr_pos,
    __global const float4* bboxes_min,
    __global const float4* bboxes_max,
    __global const int4*   neighs_global,        // natoms int4
    __global int*          ghost_indices_flat,   // [ngroup*MAX_GHOSTS]
    __global int*          ghost_counts,         // [ngroup]
    __global const int*    cluster_counts,       // [ngroup] number of real atoms in each cluster
    __global int4*         neighs_local,         // natoms int4 (local indices 0..GROUP_SIZE+nghost)
    const int num_atoms,
    const int num_groups,
    const float margin_sq,
    const float bbox_margin,
    const int group_stride
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);

    __local int   l_ghost_list[MAX_GHOSTS];
    __local int   l_ghost_counter;
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

    for (int other_g = 0; other_g < num_groups; other_g++) {
        if (other_g == grp) continue;

        float4 o_min = bboxes_min[other_g];
        float4 o_max = bboxes_max[other_g];

        bool overlap = false;
        if (my_max.x + bbox_margin >= o_min.x && my_min.x <= o_max.x + bbox_margin &&
            my_max.y + bbox_margin >= o_min.y && my_min.y <= o_max.y + bbox_margin) {
            overlap = true;
        }

        if (overlap) {
            int global_idx = other_g * group_stride + lid;
            if (global_idx < num_atoms) {
                float2 p = curr_pos[global_idx];
                float dx = max(0.0f, max(my_min.x - p.x, p.x - my_max.x));
                float dy = max(0.0f, max(my_min.y - p.y, p.y - my_max.y));
                float dist_sq = dx*dx + dy*dy;
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

    int total_ghosts = min(l_ghost_counter, MAX_GHOSTS);
    int base_offset = grp * MAX_GHOSTS;
    for (int i = lid; i < total_ghosts; i += GROUP_SIZE) {
        ghost_indices_flat[base_offset + i] = l_ghost_list[i];
    }
    if (lid == 0) ghost_counts[grp] = total_ghosts;

    int my_global_id = grp * group_stride + lid;
    if (my_global_id < num_atoms) {
        int4 ng = neighs_global[my_global_id];
        int* neighbors = (int*)&ng;
        int4 nl = (int4)(-1,-1,-1,-1);
        int* nloc = (int*)&nl;

        int c_count = cluster_counts[grp];
        for (int k = 0; k < MAX_DEGREE; k++) {
            int target = neighbors[k];
            if (target < 0) continue;
            int target_grp = target / group_stride;
            if (target_grp == grp) {
                nloc[k] = target % group_stride;
            } else {
                int loc = -1;
                for (int g = 0; g < total_ghosts; g++) {
                    if (l_ghost_list[g] == target) { loc = c_count + g; break; }
                }
                nloc[k] = loc;
            }
        }
        neighs_local[my_global_id] = nl;
    }
}

 // DEPRECATED __kernel void solve_cluster_jacobi_2d(


__kernel void compute_collision_cluster(
    __global const float2* curr_pos,           // R
    __global const float*  radius,             // R
    __global const float2* mass,               // R: mass[i].x = M
    __global const int4*   neighs_local,       // natoms int4 local indices
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const int*    cluster_counts,
    __global float2*       dpos_node,          // W: collision correction per atom (natoms-sized)
    const int num_atoms,
    const float dt,
    const float k_coll,
    const int group_stride
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * group_stride + lid;

    __local float2 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float2 l_pos_new[GROUP_SIZE];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];

    if (my_global_id < num_atoms) {
        l_pos[lid] = curr_pos[my_global_id];
        l_rad[lid] = radius[my_global_id];
    } else {
        l_pos[lid] = (float2)(0.0f, 0.0f);
        l_rad[lid] = 0.0f;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    float my_mass = (my_global_id < num_atoms) ? mass[my_global_id].x : 1.0f;
    float alpha = my_mass / (dt*dt + 1e-20f);
    float my_rad = (my_global_id < num_atoms) ? l_rad[lid] : 0.0f;

    int g_count  = ghost_counts[grp];
    int c_count  = cluster_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int g_idx = ghost_indices_flat[g_offset + k];
        int lj = c_count + k;
        l_pos[lj] = curr_pos[g_idx];
        l_rad[lj] = radius[g_idx];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        float2 p = l_pos[lid];
        float2 rhs = p * alpha;
        float  k_sum = alpha;

        // Collisions against internal + ghosts, skip if bonded via neighs_local
        // NOTE: Bonded neighbors are excluded from collision detection.
        // TODO: May need 2nd neighbor exclusion list (up to 16 entries per atom)
        //       for more complex force field interactions.
        int4 ng = neighs_local[my_global_id];
        int* neigh = (int*)&ng;

        int n_ext = c_count + g_count;
        for (int j = 0; j < n_ext; j++) {
            if (j == lid) continue;
            if (j >= c_count && j < GROUP_SIZE) continue; // skip padding lanes beyond real atoms
            bool is_bonded = false;
            for (int kk = 0; kk < MAX_DEGREE; kk++) {
                int nb = neigh[kk];
                if (nb < 0) continue;
                if (nb == j) { is_bonded = true; break; }
            }
            if (is_bonded) continue;

            float2 q = l_pos[j];
            float2 d = p - q;
            float dist_sq = dot(d, d);
            float r_sum = my_rad + l_rad[j];
            float r2 = r_sum * r_sum;
            if (dist_sq < r2 && dist_sq > 1e-16f) {
                float dist = sqrt(dist_sq);
                float coeff = k_coll * (r_sum / dist);
                rhs += d * coeff;
                rhs += q * k_coll;
                k_sum += k_coll;
            }
        }

        l_pos_new[lid] = rhs / (k_sum + 1e-20f);
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        l_pos[lid] = l_pos_new[lid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        float2 p0   = curr_pos[my_global_id];
        float2 p1   = l_pos[lid];
        dpos_node[my_global_id] = p1 - p0;
    }
}

// ================================================================
// ================================================================
//      POSITION BASED DYNAMICS (PBD) / Projective Dynamics (PD)
// ================================================================
// ================================================================


// ------------------------------------------------------------------
// KERNEL: Initialize momentum buffers to zero
// ------------------------------------------------------------------
__kernel void init_mom_pos_2d(
    const int natoms,
    __global float2* mom_pos
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    mom_pos[i] = (float2)(0.0f, 0.0f);
}

__kernel void init_mom_rot_2d(
    const int nnode,
    __global float2* mom_rot
) {
    int i = get_global_id(0);
    if (i >= nnode) return;
    mom_rot[i] = (float2)(0.0f, 0.0f);
}

/*
### The Fundamental Problem

You sensed that standard Jacobi-PD (your original code) breaks momentum conservation. Here is why:

In Projective Dynamics, the global step is a linear system: $(M/dt^2 + L) \Delta q = f_{rhs}$.
If you solve this system **exactly** (e.g., Cholesky, CG), momentum is conserved.
However, on the GPU, we often use **Jacobi iteration** (diagonal approximation).

*   **Jacobi PD:** Node $i$ calculates its next position assuming neighbors stay fixed. Node $j$ calculates its position assuming node $i$ stays fixed.
*   **The math:** Node $i$ moves by $\approx \frac{F_{ij}}{M_i/dt^2 + K}$. Node $j$ moves by $\approx \frac{-F_{ij}}{M_j/dt^2 + K}$.
*   **The Error:** Unless $M_i = M_j$, the center of mass moves! $\Delta P \propto F (\frac{M_i}{M_i + \dots} - \frac{M_j}{M_j + \dots}) \neq 0$.

**Conclusion:** You cannot use standard node-based Jacobi iteration (row-based solve) if you require strict momentum conservation per iteration. You must use a **constraint-based (edge-based)** approach.

### Is the solution "Projective Dynamics"?

**Yes, mathematically.**
XPBD (the code I provided previously) is exactly the analytic solution to the PD minimization problem for a single constraint.
*   PD minimizes Energy $E = \frac{1}{2} \frac{M}{dt^2} ||\Delta x||^2 + \frac{1}{2} K ||C(x)||^2$.
*   If you take the derivative of this energy with respect to $x$ and set it to zero for a single constraint, you get **exactly** the XPBD update formula.

So, to keep the "PD spirit" (using $M/dt^2$ and $K$ explicitly) while fixing momentum, we just need to ensure we calculate the **Coupled Impulse** using the correct inertial terms.

### The Solution: Coupled "Block-Jacobi" Solver

Instead of decoupling translation and rotation (which you correctly identified as a big issue), we solve the local system for the constraint exactly. This handles the lever arm $r$ correctly, coupling translation and rotation instantly.

Here is the 2D Kernel. It uses your specific PD inputs ($M, I, K, dt$) but solves the constraint pairwise to guarantee $\Delta P = 0$.

### Why this works (and why it is better than Force-based)

1.  **Implicit Scaling (`alpha`):** The term `alpha = 1.0f / (K * dt2)` provides the infinite stiffness stability of PD/Implicit integration. If `K` is huge, `alpha` is small, and the constraint is solved rigidly. If `K` is small, the constraint is springy.
2.  **Coupling (`w_angular`):** By including `rn * rn * w_rot_i` in the denominator (`w_total`), we acknowledge that it is harder to satisfy a constraint if it requires spinning a heavy object. This prevents the "vibration" often seen when rotation and translation are solved separately.
3.  **Conservation:**
    *   Node $i$ adds $P$.
    *   Neighbor $j$ adds $-P$.
    *   Linear Momentum change: $\sum \Delta p = M_i(P/M_i) + M_j(-P/M_j) = P - P = 0$.

*/

// ------------------------------------------------------------------
// KERNEL: Projective Dynamics / XPBD with Momentum Conservation (2D)
// ------------------------------------------------------------------
__kernel void compute_corrections_2d(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float2* port_local,
    __global const float*  stiffness_flat,
    __global const float2* mass,            // mass[i].x = M, mass[i].y = I
    __global float2*       dpos_node,       // Output: Delta P for Node i
    __global float*        dtheta_node,     // Output: Delta Theta for Node i
    __global float2*       dpos_neigh,      // Output: Delta P for Neighbor j (Recoil)
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    // 1. Get Inertial Properties
    // In PD, the diagonal term is M / dt^2.
    // In the update formula, we need inverse mass: w = 1/M.
    float M_i = mass[i].x;
    float I_i = mass[i].y;
    
    // Inverse masses
    float w_i     = (M_i > 1e-9f) ? 1.0f / M_i : 0.0f;
    float w_rot_i = (I_i > 1e-9f) ? 1.0f / I_i : 0.0f;
    
    float dt2 = dt * dt + 1e-16f;

    // Accumulators for Node i
    float2 sum_dpos   = (float2)(0.0f, 0.0f);
    float  sum_dtheta = 0.0f;

    int4 ng        = neighs[i];
    int* neighbors = (int*)&ng;
    int i4         = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int idx = i4 + k;
        int j = neighbors[k];
        
        // Reset recoil slot
        dpos_neigh[idx] = (float2)(0.0f, 0.0f);

        if (j < 0) continue;

        float K = stiffness_flat[idx];
        if (K <= 0.0f) continue;

        // --- Geometry & Constraint Evaluation ---
        // Port location in World Space
        float2 r_local = port_local[idx];
        float2 r_arm   = cmplx_rotate(z_i, r_local); 
        float2 tip     = p_i + r_arm;
        float2 p_j     = pos[j]; // Neighbor position

        // Constraint vector: We want Tip and Neighbor to coincide.
        // Vector from Tip -> Neighbor
        float2 diff = p_j - tip; 
        float dist  = length(diff);

        if (dist < 1e-9f) continue;

        float2 n = diff / dist; // Direction of the constraint force

        // --- The "PD" Momentum Conserving Solve ---

        // 1. Compliance (The PD Regularizer)
        // PD Energy: 0.5 * K * x^2. 
        // This corresponds to compliance alpha = 1 / (K * dt^2)
        // If K is infinite (rigid), alpha is 0.
        float alpha = 1.0f / (K * dt2);

        // 2. Neighbor Mass (Assuming Point Mass for neighbor j for simplicity)
        // If j is also a node, you technically need its rotation too, but 
        // standard implementation often treats the "other end" as a point 
        // in the scatter step to avoid complex locking.
        float M_j = mass[j].x;
        float w_j = (M_j > 1e-9f) ? 1.0f / M_j : 0.0f;

        // 3. Generalized Inverse Mass (Coupling Translation & Rotation)
        // This solves the "Fundamental Problem" of coupling.
        // How much does the system move along 'n' if we apply unit impulse?
        // Linear part: w_i + w_j
        // Angular part: (r x n)^2 / I
        
        // 2D Cross product magnitude (r_arm cross n)
        float rn        = r_arm.x * n.y - r_arm.y * n.x; 
        float w_angular = rn * rn * w_rot_i;

        // Total effective inverse mass of this constraint
        float w_total   = w_i + w_j + w_angular + alpha;

        // 4. Calculate Impulse Magnitude (Lagrange Multiplier)
        // We want to correct the gap 'dist'.
        // lambda = -C / w_total
        float impulse_mag = dist / w_total;

        // 5. Symmetric Distribution (Factor of 0.5)
        // Since we process bond i-j AND bond j-i (in separate threads),
        // we essentially calculate this force twice. We must halve it 
        // to avoid double stiffness.
        impulse_mag *= 0.5f;

        // 6. Apply Impulses (Momentum Conserved!)
        // Impulse Vector P
        float2 P = n * impulse_mag;

        // Correct Node i (Move towards neighbor)
        // Linear: dx = P / M_i
        sum_dpos += P * w_i;
        
        // Angular: dtheta = (r x P) / I_i
        // r x P = r x (n * mag) = (r x n) * mag = rn * impulse_mag
        sum_dtheta += (rn * impulse_mag) * w_rot_i;

        // Correct Neighbor j (Recoil - Move towards Node)
        // Linear: dx = -P / M_j
        dpos_neigh[idx] = -P * w_j;
    }

    // Store accumulated corrections for Node i
    dpos_node[i]   = sum_dpos;
    dtheta_node[i] = sum_dtheta;
}


// ------------------------------------------------------------------
// KERNEL: Apply corrections with relaxation and proper momentum mixing (2D)
// ------------------------------------------------------------------
__kernel void apply_corrections_2d(
    const int natoms,
    const int nnode,
    __global float2*       pos,
    __global float2*       rot,
    __global const int4*   bkSlots,
    __global const float2* dpos_node,
    __global const float*  dtheta_node,
    __global const float2* dpos_neigh,
    const float relaxation,
    __global float2*       hb_pos,
    __global float2*       hb_rot,
    __global float2*       mom_pos,
    __global float2*       mom_rot,
    const float bmix_pos,
    const float bmix_rot
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    // Apply position correction
    float2 corr = (float2)(0.0f, 0.0f);

    // Gather corrections from neighbors
    int4 bk = bkSlots[i];
    if (bk.x >= 0) corr += dpos_neigh[bk.x];
    if (bk.y >= 0) corr += dpos_neigh[bk.y];
    if (bk.z >= 0) corr += dpos_neigh[bk.z];
    if (bk.w >= 0) corr += dpos_neigh[bk.w];

    // note - capping atoms (i>=nnode) are double-weighted becaous they are not double counted in compute_corrections_2d
    if (i>=nnode){ corr*=2.0f;}
    corr += dpos_node[i];

    // Relaxation step
    float2 p_corr = pos[i] + corr * relaxation;

    // Heavy-ball momentum mixing for position (C++ style)
    if (bmix_pos > 1e-6f) {
        float2 d_prev = mom_pos[i];        // momentum direction from previous iter
        float2 p_prev = hb_pos[i];         // previous iterate
        float2 p_next = p_corr + d_prev * bmix_pos;
        mom_pos[i] = p_next - p_prev;      // store new momentum direction
        hb_pos[i] = p_next;
        pos[i] = p_next;
    } else {
        pos[i] = p_corr;
        if (bmix_pos > -0.5f) {            // Only update hb if not explicitly disabled
            hb_pos[i] = p_corr;
        }
    }

    // Apply rotation correction (nodes only)
    if (i < nnode) {
        float dtheta = dtheta_node[i] * relaxation;
        float2 dz = cmplx_from_angle(dtheta);
        float2 r_corr = cmplx_mul(dz, rot[i]);
        
        // Heavy-ball momentum mixing for rotation
        if (bmix_rot > 1e-6f) {
            float2 d_prev = mom_rot[i];
            float2 r_prev = hb_rot[i];
            // For rotation, we apply momentum in the tangent space
            // r_next = normalize(r_corr + bmix * momentum direction)
            float2 r_next = r_corr + d_prev * bmix_rot;
            mom_rot[i] = r_next - r_prev;
            hb_rot[i] = r_next;
            rot[i] = cmplx_normalize(r_next);
        } else {
            rot[i] = cmplx_normalize(r_corr);
            if (bmix_rot > -0.5f) {
                hb_rot[i] = cmplx_normalize(r_corr);
            }
        }
    }
}

// ================================================================
// ================================================================
//    XPBD with Lagrange Multipliers ( XPBD-LM ) - Does not work
// ================================================================
// ================================================================

// ------------------------------------------------------------------
// KERNEL: Reset lambda for XPBD
// ------------------------------------------------------------------
__kernel void reset_lambda_2d(
    const int n,
    __global float* lambda
) {
    int i = get_global_id(0);
    if (i >= n) return;
    lambda[i] = 0.0f;
}

// ------------------------------------------------------------------
// KERNEL: XPBD constraint solve with lambda accumulation (2D)
// ------------------------------------------------------------------
__kernel void compute_xpbd_corrections_2d(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float4* stiffness,     // per-node stiffness
    __global const float2* port_local,
    __global const float2* mass,          // mass.x=M, mass.y=I
    __global float*        lambda,        // accumulated multipliers [nnode*4]
    __global float2*       dpos_neigh,
    __global float2*       dpos_node,
    __global float*        dtheta_node,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        if (j < 0) continue;

        float K = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (K <= 0.0f) continue;

        int idx = i4 + k;
        float lambda_prev = lambda[idx];

        float2 r_local = port_local[idx];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;

        float2 p_j = pos[j];

        float2 diff = tip - p_j;
        float C2 = dot(diff, diff);
        if (C2 < 1e-16f) continue;

        float C = sqrt(C2);
        float2 n = diff / C;

        // Generalized mass (translation + rotation coupling) + XPBD compliance
        // Mass properties from buffers: mass[i].x = M, mass[i].y = I
        float M_i = mass[i].x;
        float I_i = mass[i].y;
        float M_j = mass[j].x;
        float dt2 = dt * dt;

        float r_cross_n = r_world.x * n.y - r_world.y * n.x;

        // Reduced-mass style weights for proper momentum conservation
        float w_i_trans = (M_i > 1e-9f) ? 1.0f / M_i : 0.0f;
        float w_i_rot   = (I_i > 1e-12f) ? (r_cross_n * r_cross_n) / I_i : 0.0f;
        float w_j_trans = (M_j > 1e-9f) ? 1.0f / M_j : 0.0f;

        float alpha = 1.0f / (K + 1e-12f);
        float alpha_tilde = alpha / (dt2 + 1e-20f);

        float w_total = w_i_trans + w_i_rot + w_j_trans + alpha_tilde;
        float delta_lambda = (-C - alpha_tilde * lambda_prev) / (w_total + 1e-12f);
        lambda[idx] = lambda_prev + delta_lambda;

        float2 dp_i = delta_lambda * n * w_i_trans;
        float2 dp_j = -delta_lambda * n * w_j_trans;
        float dtheta_k = (I_i > 1e-12f) ? (delta_lambda * r_cross_n / I_i) : 0.0f;

        dpos += dp_i;
        dtheta += dtheta_k;
        dpos_neigh[idx] = dp_j;
    }

    dpos_node[i] = dpos;
    dtheta_node[i] = dtheta;
}

// ------------------------------------------------------------------
// KERNEL: XPBD Corrections with DEBUG PRINTS (diagnostic version)
// ------------------------------------------------------------------
__kernel void compute_xpbd_corrections_2d_debug(
    const int nnode,
    __global const float2* pos,
    __global const float2* rot,
    __global const int4*   neighs,
    __global const float4* stiffness,     // per-node stiffness
    __global const float2* port_local,
    __global const float2* mass,          // mass.x=M, mass.y=I
    __global float*        lambda,        // accumulated multipliers [nnode*4]
    __global float2*       dpos_neigh,
    __global float2*       dpos_node,
    __global float*        dtheta_node,
    const float dt,
    const int debug_step,                 // current step for gating output
    const int max_debug_steps,            // max steps to output (e.g., 5)
    // Diagnostic output buffers (only written for debug_step < max_debug_steps)
    __global float*        dbg_data,      // [nnode*4*5] = C, lambda, dtheta, K, alpha per constraint
    __global float2*       dbg_dpos_i,    // [nnode*4] position correction on node i
    __global float2*       dbg_dpos_j,    // [nnode*4] position correction on neighbor j
    __global float2*       dbg_r_world,   // [nnode*4] rotated port vector
    __global float2*       dbg_n          // [nnode*4] constraint normal
) {
    int i = get_global_id(0);
    if (i >= nnode) return;

    float2 p_i = pos[i];
    float2 z_i = rot[i];

    int4 ng = neighs[i];
    int* neighbors = (int*)&ng;
    float2 dpos = (float2)(0.0f, 0.0f);
    float dtheta = 0.0f;

    int i4 = i * MAX_DEGREE;
    int base_offset = i4;  // For diagnostic arrays

    for (int k = 0; k < MAX_DEGREE; k++) {
        int j = neighbors[k];
        int idx = i4 + k;

        // Initialize diagnostics to zero (for inactive constraints)
        if (debug_step < max_debug_steps) {
            dbg_data[base_offset + k] = 0.0f;           // C
            dbg_data[base_offset + k + nnode*4] = 0.0f; // lambda
            dbg_data[base_offset + k + nnode*8] = 0.0f; // dtheta
            dbg_data[base_offset + k + nnode*12] = 0.0f; // K
            dbg_data[base_offset + k + nnode*16] = 0.0f; // alpha
            dbg_dpos_i[idx] = (float2)(0.0f, 0.0f);
            dbg_dpos_j[idx] = (float2)(0.0f, 0.0f);
            dbg_r_world[idx] = (float2)(0.0f, 0.0f);
            dbg_n[idx] = (float2)(0.0f, 0.0f);
        }

        if (j < 0) continue;

        float K = (k == 0) ? stiffness[i].x : (k == 1) ? stiffness[i].y : (k == 2) ? stiffness[i].z : stiffness[i].w;
        if (K <= 0.0f) continue;

        float lambda_prev = lambda[idx];

        float2 r_local = port_local[idx];
        float2 r_world = cmplx_rotate(z_i, r_local);
        float2 tip = p_i + r_world;

        float2 p_j = pos[j];

        float2 diff = tip - p_j;
        float C2 = dot(diff, diff);
        if (C2 < 1e-16f) continue;

        float C = sqrt(C2);
        float2 n = diff / C;

        // Generalized mass (translation + rotation coupling) + XPBD compliance
        // Mass properties from buffers: mass[i].x = M, mass[i].y = I
        float M_i = mass[i].x;
        float I_i = mass[i].y;
        float M_j = mass[j].x;
        float dt2 = dt * dt;

        float r_cross_n = r_world.x * n.y - r_world.y * n.x;

        // Reduced-mass style weights for proper momentum conservation
        float w_i_trans = (M_i > 1e-9f) ? 1.0f / M_i : 0.0f;
        float w_i_rot   = (I_i > 1e-12f) ? (r_cross_n * r_cross_n) / I_i : 0.0f;
        float w_j_trans = (M_j > 1e-9f) ? 1.0f / M_j : 0.0f;

        float alpha = 1.0f / (K + 1e-12f);
        float alpha_tilde = alpha / (dt2 + 1e-20f);

        float w_total = w_i_trans + w_i_rot + w_j_trans + alpha_tilde;
        float delta_lambda = (-C - alpha_tilde * lambda_prev) / (w_total + 1e-12f);
        lambda[idx] = lambda_prev + delta_lambda;

        float2 dp_i = delta_lambda * n * w_i_trans;
        float2 dp_j = -delta_lambda * n * w_j_trans;
        float dtheta_k = (I_i > 1e-12f) ? (delta_lambda * r_cross_n / I_i) : 0.0f;

        // Store diagnostics for active constraints
        if (debug_step < max_debug_steps) {
            dbg_data[base_offset + k] = C;
            dbg_data[base_offset + k + nnode*4] = lambda[idx];
            dbg_data[base_offset + k + nnode*8] = dtheta_k;
            dbg_data[base_offset + k + nnode*12] = K;
            dbg_data[base_offset + k + nnode*16] = alpha;
            dbg_dpos_i[idx] = dp_i;
            dbg_dpos_j[idx] = dp_j;
            dbg_r_world[idx] = r_world;
            dbg_n[idx] = n;
        }

        dpos += dp_i;
        dtheta += dtheta_k;
        dpos_neigh[idx] = dp_j;
    }

    dpos_node[i] = dpos;
    dtheta_node[i] = dtheta;
}

// ------------------------------------------------------------------
// KERNEL: Gather and apply XPBD corrections (2D) with momentum buffers
// ------------------------------------------------------------------
__kernel void gather_and_apply_xpbd_2d(
    const int natoms,
    const int nnode,
    __global float2*       pos,
    __global float2*       rot,
    __global const int4*   bkSlots,
    __global const float2* dpos_neigh,
    __global const float2* dpos_node,
    __global const float*  dtheta_node,
    __global float2*       hb_pos,
    __global float2*       hb_rot,
    __global float2*       mom_pos,
    __global float2*       mom_rot,
    const float bmix_pos,
    const float bmix_rot
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float2 corr = (float2)(0.0f, 0.0f);

    if (i < nnode) {
        corr += dpos_node[i];
    }

    int4 bk = bkSlots[i];
    if (bk.x >= 0) corr += dpos_neigh[bk.x];
    if (bk.y >= 0) corr += dpos_neigh[bk.y];
    if (bk.z >= 0) corr += dpos_neigh[bk.z];
    if (bk.w >= 0) corr += dpos_neigh[bk.w];

    // For XPBD, no relaxation factor - corrections are already scaled
    float2 p_corr = pos[i] + corr;

    // Heavy-ball momentum mixing for position
    if (bmix_pos > 1e-6f) {
        float2 d_prev = mom_pos[i];
        float2 p_prev = hb_pos[i];
        float2 p_next = p_corr + d_prev * bmix_pos;
        mom_pos[i] = p_next - p_prev;
        hb_pos[i] = p_next;
        pos[i] = p_next;
    } else {
        pos[i] = p_corr;
        if (bmix_pos > -0.5f) {
            hb_pos[i] = p_corr;
        }
    }

    if (i < nnode) {
        float dtheta = dtheta_node[i];
        float2 dz = cmplx_from_angle(dtheta);
        float2 r_corr = cmplx_mul(dz, rot[i]);
        
        // Heavy-ball momentum mixing for rotation
        if (bmix_rot > 1e-6f) {
            float2 d_prev = mom_rot[i];
            float2 r_prev = hb_rot[i];
            float2 r_next = r_corr + d_prev * bmix_rot;
            mom_rot[i] = r_next - r_prev;
            hb_rot[i] = r_next;
            rot[i] = cmplx_normalize(r_next);
        } else {
            rot[i] = cmplx_normalize(r_corr);
            if (bmix_rot > -0.5f) {
                hb_rot[i] = cmplx_normalize(r_corr);
            }
        }
    }
}
