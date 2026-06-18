// RRsp3.cl - Cluster-sorted rigid ports + collisions (PBD/Jacobi) with recoils
//
// Kernel dependency graph (state buffers only; topology / index params omitted)
//
// Common frame (all methods):
//   update_bboxes_rigid(pos,radius) -> bboxes
//   build_local_topology_rigid(pos,neighs,excl,bboxes) -> neighs_local, excl_local
//   compute_collision_cluster_rigid(pos,radius,excl_local) -> dpos_coll
//   apply_corrections_rigid_ports(pos,quat,dpos_node,drot_node,dpos_neigh,dpos_coll) -> pos,quat[,tips]
//
// Massfull rotation (physical rotational inertia):
//   current:
//     compute_ports_cluster_rigid(pos,quat,port_local) -> dpos_node, drot_node, dpos_neigh
//   orig:
//     compute_ports_cluster_rigid_orig(pos,quat,port_local) -> dpos_node, drot_node, dpos_neigh
//
// Massless rotation (geometric alignment, zero inertia):
//   substep_optimized:
//     compute_ports_cluster_rigid_substep_optimized(pos,quat,port_local) -> dpos_node, drot_node, dpos_neigh
//   shapematch:
//     compute_ports_cluster_rigid_shapematch(pos,quat,port_local) -> dpos_node, drot_node, dpos_neigh
//   eigen (two-pass):
//     compute_tips(pos,quat,port_local) -> tips
//     compute_optimal_rotation_eigen(pos,quat,port_local) -> drot_node, quat_opt(dquat_mom)
//     compute_ports_cluster_rigid_eigen_tips(pos,tips) -> dpos_node, dpos_neigh


// Only enable prints if this flag is set during compilation
#ifdef ENABLE_DEBUG_PRINTS

// Verbosity: 0=none, 1=errors, 2=summary per WG, 3=per-atom
#ifndef DEBUG_VERBOSITY
#define DEBUG_VERBOSITY 1
#endif

// Component bitmask: 1=collision, 2=port, 4=topology, 8=correction
#ifndef DEBUG_COMPONENTS
#define DEBUG_COMPONENTS 0xFFFF
#endif

// Workgroup targeting: -1 = all workgroups
#ifndef DEBUG_TARGET_WG
#define DEBUG_TARGET_WG -1
#endif

// Atom range targeting
#ifndef DEBUG_GID_START
#define DEBUG_GID_START 0
#endif
#ifndef DEBUG_GID_END
#define DEBUG_GID_END   8
#endif

#define DBG_COLL  1
#define DBG_PORT  2
#define DBG_TOPO  4
#define DBG_CORR  8

#define DBG_WG_OK   (DEBUG_TARGET_WG < 0 || get_group_id(0) == DEBUG_TARGET_WG)
#define DBG_ATOM_OK(gid) (DEBUG_GID_START <= gid && gid < DEBUG_GID_END)
#define DBG_COMP_OK(comp) (DEBUG_COMPONENTS & comp)

#define LOG_TOPOLOGY_SUMMARY(grp, n_ghosts)  \
    if (DEBUG_VERBOSITY >= 2 && DBG_WG_OK && DBG_COMP_OK(DBG_TOPO)) \
        printf("TOPO_SUM: WG=%d grp=%d n_ghosts=%d\n", get_group_id(0), grp, n_ghosts);

#define LOG_TOPOLOGY(gid, lid, n_ghosts)  \
    if (DEBUG_VERBOSITY >= 3 && DBG_WG_OK && DBG_ATOM_OK(gid) && DBG_COMP_OK(DBG_TOPO)) \
        printf("TOPOLOGY: GID=%d LID=%d n_ghosts=%d WG=%d\n", gid, lid, n_ghosts, get_group_id(0));

#define LOG_MAPPING(gid, local_idx, mapped_global_idx, type)  \
    if (DEBUG_VERBOSITY >= 3 && DBG_WG_OK && DBG_ATOM_OK(gid) && DBG_COMP_OK(DBG_TOPO)) \
        printf("MAP: GID=%d L_IDX=%d G_IDX=%d TYPE=%s WG=%d\n", gid, local_idx, mapped_global_idx, type, get_group_id(0));

#define LOG_COLLISION_CHECK(gid, my_global, other_local, other_global, dist, action) \
    if (DEBUG_VERBOSITY >= 3 && DBG_WG_OK && DBG_ATOM_OK(gid) && DBG_COMP_OK(DBG_COLL)) \
        printf("COLL: MeG=%d OtherL=%d OtherG=%d Dist=%.4f Action=%s WG=%d\n", my_global, other_local, other_global, dist, action, get_group_id(0));

#define LOG_PORT(gid, inode, k, dist, impulse) \
    if (DEBUG_VERBOSITY >= 3 && DBG_WG_OK && DBG_ATOM_OK(gid) && DBG_COMP_OK(DBG_PORT)) \
        printf("PORT: GID=%d Inode=%d k=%d Dist=%.6f Impulse=%.6f WG=%d\n", gid, inode, k, dist, impulse, get_group_id(0));

#define LOG_CORRECTION(gid, dx, dtheta) \
    if (DEBUG_VERBOSITY >= 3 && DBG_WG_OK && DBG_ATOM_OK(gid) && DBG_COMP_OK(DBG_CORR)) \
        printf("CORR: GID=%d dx=(%.6f,%.6f,%.6f) dtheta=(%.6f,%.6f,%.6f) WG=%d\n", gid, dx.x, dx.y, dx.z, dtheta.x, dtheta.y, dtheta.z, get_group_id(0));

#else
// Empty macros for production
#define LOG_TOPOLOGY(gid, lid, n_ghosts)
#define LOG_MAPPING(gid, local_idx, mapped_global_idx, type)
#define LOG_COLLISION_CHECK(gid, my_global, other_local, other_global, dist, action)
#define LOG_PORT(gid, inode, k, dist, impulse)
#define LOG_CORRECTION(gid, dx, dtheta)
#define LOG_TOPOLOGY_SUMMARY(grp, n_ghosts)
#endif

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#ifndef GROUP_SIZE
#define GROUP_SIZE     64
#endif
#ifndef MAX_GHOSTS
#define MAX_GHOSTS     128
#endif

#ifndef ENABLE_COLL
#define ENABLE_COLL 1
#endif

#ifndef ENABLE_PORT
#define ENABLE_PORT 1
#endif

// =========================================================
// =========================================================
//  Helper, Mmatrix, Quaternion, AABB 
// =========================================================
// =========================================================



bool bboxes_overlap(float4 minA, float4 maxA, float4 minB, float4 maxB, float margin) {
    if (maxA.x + margin < minB.x || minA.x > maxB.x + margin) return false;
    if (maxA.y + margin < minB.y || minA.y > maxB.y + margin) return false;
    if (maxA.z + margin < minB.z || minA.z > maxB.z + margin) return false;
    return true;
}

inline float3 quat_rotate(float4 q, float3 v) {
    float3 t = 2.0f * cross(q.xyz, v);
    return v + q.w * t + cross(q.xyz, t);
}
inline float3 q_rot(float4 q, float3 v) { return quat_rotate(q, v); }

inline float4 quat_from_axis_angle(float3 axis, float angle) {
    float a = length(axis);
    if (a < 1e-8f || fabs(angle) < 1e-8f) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    float3 n = axis / a;
    float s = sin(angle * 0.5f);
    return (float4)(n * s, cos(angle * 0.5f));
}

inline float4 quat_mul(float4 a, float4 b) {
    return (float4)(
        a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
        a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
        a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w,
        a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z
    );
}

inline float4 quat_conj(float4 q){ return (float4)(-q.x, -q.y, -q.z, q.w); }

inline float3 quat_delta_rotvec(float4 q_new, float4 q_old){
    float4 dq = quat_mul(q_new, quat_conj(q_old));
    float sign = (dq.w < 0.0f) ? -1.0f : 1.0f;
    return 2.0f * dq.xyz * sign;
}

inline float4 quat_normalize(float4 q){
    float n2 = dot(q,q);
    if(n2<1e-16f) return (float4)(0.0f,0.0f,0.0f,1.0f);
    return q * rsqrt(n2);
}

inline float3 solve_3x3(float3 r0, float3 r1, float3 r2, float3 b) {
    // Solve A x = b for A with rows r0,r1,r2 using adjugate.
    float3 c0 = cross(r1, r2);
    float3 c1 = cross(r2, r0);
    float3 c2 = cross(r0, r1);
    float det = dot(r0, c0);
    if (fabs(det) < 1e-20f) return (float3)(0.0f);
    float invDet = 1.0f / det;
    // inv(A) = [c0 c1 c2]^T / det
    return (float3)( dot(b, (float3)(c0.x,c1.x,c2.x)), dot(b, (float3)(c0.y,c1.y,c2.y)), dot(b, (float3)(c0.z,c1.z,c2.z)) ) * invDet;
}

inline float4 apply_delta_rot(float4 q, float3 dtheta) {
    float angle = length(dtheta);
    if (angle < 1e-8f) return q;
    float3 axis = dtheta / angle;
    float4 dq = quat_from_axis_angle(axis, angle);
    return quat_normalize(quat_mul(dq, q));
}

inline void add_hessian_rr(float3* H0, float3* H1, float3* H2, float w, float3 r){
    float r2 = dot(r, r);
    (*H0).x += w * (r2 - r.x*r.x);
    (*H0).y += w * (   - r.x*r.y);
    (*H0).z += w * (   - r.x*r.z);
    (*H1).x += w * (   - r.y*r.x);
    (*H1).y += w * (r2 - r.y*r.y);
    (*H1).z += w * (   - r.y*r.z);
    (*H2).x += w * (   - r.z*r.x);
    (*H2).y += w * (   - r.z*r.y);
    (*H2).z += w * (r2 - r.z*r.z);
}

inline void mat3_outer_add(float* A, float w, float3 a, float3 b){
    A[0] += w * a.x * b.x;  A[1] += w * a.x * b.y;  A[2] += w * a.x * b.z;
    A[3] += w * a.y * b.x;  A[4] += w * a.y * b.y;  A[5] += w * a.y * b.z;
    A[6] += w * a.z * b.x;  A[7] += w * a.z * b.y;  A[8] += w * a.z * b.z;
}

inline int excluded8(int j, int4 a, int4 b){
    if( (j==a.x) || (j==a.y) || (j==a.z) || (j==a.w) ) return 1;
    if( (j==b.x) || (j==b.y) || (j==b.z) || (j==b.w) ) return 1;
    return 0;
}

// Multiply 3x3 matrices (Column Major or Row Major consistent)
inline void mat3_mul(float* A, float* B, float* Out) {
    // Unrolled for performance
    for (int r = 0; r < 3; r++) {
        for (int c = 0; c < 3; c++) {
            Out[r*3 + c] = A[r*3 + 0] * B[0*3 + c] + 
                           A[r*3 + 1] * B[1*3 + c] + 
                           A[r*3 + 2] * B[2*3 + c];
        }
    }
}

// ------------------------------------------------------------------
// DYNAMICS (leapfrog/PBD-style)
// ------------------------------------------------------------------

__kernel void predict_dynamics(
    const int natoms,
    const int nnode_per_group,
    __global float4* pos,
    __global float4* quat,
    __global const int* fixmask,
    __global float4* vel,
    __global float4* omega,
    __global float4* pos_prev,
    __global float4* quat_prev,
    const float dt
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float4 p4 = pos[i];
    float invMi = p4.w;
    pos_prev[i] = p4;
    quat_prev[i] = quat[i];

    if (invMi <= 1e-12f) return; // padding only

    int msk = (fixmask != 0) ? fixmask[i] : 0;

    float3 v = vel[i].xyz;
    float3 p = p4.xyz + v * dt;

    if (msk & 1) { p.x = p4.x; v.x = 0.0f; }
    if (msk & 2) { p.y = p4.y; v.y = 0.0f; }
    if (msk & 4) { p.z = p4.z; v.z = 0.0f; }
    if (msk & 8) { p.z = 0.0f; v.z = 0.0f; }

    pos[i] = (float4)(p, invMi);
    vel[i] = (float4)(v, 0.0f);

    int lid = i & (GROUP_SIZE - 1);
    if (lid < nnode_per_group) {
        if (msk & (1|2|4)) {
            omega[i] = (float4)(0.0f);
            return;
        }
        float3 w = omega[i].xyz;
        float3 dtheta = w * dt;
        quat[i] = apply_delta_rot(quat[i], dtheta);
    } else {
        omega[i] = (float4)(0.0f);
    }
}

__kernel void update_velocities_dynamics(
    const int natoms,
    const int nnode_per_group,
    __global const float4* pos,
    __global const float4* quat,
    __global const int* fixmask,
    __global float4* vel,
    __global float4* omega,
    __global const float4* pos_prev,
    __global const float4* quat_prev,
    const float dt,
    const float damp
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float invMi = pos[i].w;
    if (invMi <= 1e-12f) return; // padding only

    int msk = (fixmask != 0) ? fixmask[i] : 0;

    float inv_dt = 1.0f / (dt + 1e-16f);
    float3 v = (pos[i].xyz - pos_prev[i].xyz) * inv_dt;

    if (msk & 1) v.x = 0.0f;
    if (msk & 2) v.y = 0.0f;
    if (msk & 4) v.z = 0.0f;
    if (msk & 8) v.z = 0.0f;
    v *= damp;
    vel[i] = (float4)(v, 0.0f);

    int lid = i & (GROUP_SIZE - 1);
    if (lid < nnode_per_group) {
        if (msk & (1|2|4)) {
            omega[i] = (float4)(0.0f);
            return;
        }
        float3 dtheta = quat_delta_rotvec(quat[i], quat_prev[i]);
        float3 w = dtheta * inv_dt;
        w *= damp;
        omega[i] = (float4)(w, 0.0f);
    } else {
        omega[i] = (float4)(0.0f);
    }
}

// Convert Rotation Matrix to Quaternion
// Robust method handling the trace singularities
inline float4 mat3_to_quat(float* R) {
    float4 q;
    float tr = R[0] + R[4] + R[8];
    if (tr > 0.0f) {
        float S = sqrt(tr + 1.0f) * 2.0f; // S=4*qw 
        q.w = 0.25f * S;
        q.x = (R[5] - R[7]) / S;
        q.y = (R[6] - R[2]) / S;
        q.z = (R[1] - R[3]) / S;
    } else if ((R[0] > R[4]) && (R[0] > R[8])) {
        float S = sqrt(1.0f + R[0] - R[4] - R[8]) * 2.0f; 
        q.w = (R[5] - R[7]) / S;
        q.x = 0.25f * S;
        q.y = (R[1] + R[3]) / S;
        q.z = (R[6] + R[2]) / S;
    } else if (R[4] > R[8]) {
        float S = sqrt(1.0f + R[4] - R[0] - R[8]) * 2.0f; 
        q.w = (R[6] - R[2]) / S;
        q.x = (R[1] + R[3]) / S;
        q.y = 0.25f * S;
        q.z = (R[5] + R[7]) / S;
    } else {
        float S = sqrt(1.0f + R[8] - R[0] - R[4]) * 2.0f; 
        q.w = (R[1] - R[3]) / S;
        q.x = (R[6] + R[2]) / S;
        q.y = (R[5] + R[7]) / S;
        q.z = 0.25f * S;
    }
    return normalize(q);
}

// This calculates q_out = K * q_in without storing the 4x4 matrix K.
// It uses the 9 elements of the Covariance Matrix B (Sxx, Sxy...) directly.
// K = [ Trace(B)   z^T          ]
//     [ z          B+B^T-Tr(B)I ]
inline float4 apply_K_matrix(
    float Sxx, float Sxy, float Sxz,
    float Syx, float Syy, float Syz,
    float Szx, float Szy, float Szz,
    float4 q
) {
    float4 q_out;
    
    // 1. Calculate helper terms
    float trace = Sxx + Syy + Szz;
    float3 z = (float3)(Syz - Szy, Szx - Sxz, Sxy - Syx);
    
    // 2. Row 0 (Scalar w component)
    // K_00 * w + K_01 * x + K_02 * y + K_03 * z
    q_out.w = trace * q.w + dot(z, q.xyz);
    
    // 3. Vector components (x, y, z)
    // q_vec_out = z * w + (B + B^T - tr*I) * v
    
    float3 v = q.xyz;
    
    // Diagonal terms of the block (B + B^T - tr*I)
    // Sxx + Sxx - (Sxx + Syy + Szz) = Sxx - Syy - Szz
    float3 diag = (float3)(Sxx - Syy - Szz, 
                           Syy - Sxx - Szz, 
                           Szz - Sxx - Syy);
                           
    q_out.x = z.x * q.w + diag.x * v.x + (Sxy + Syx) * v.y + (Sxz + Szx) * v.z;
    q_out.y = z.y * q.w + (Sxy + Syx) * v.x + diag.y * v.y + (Syz + Szy) * v.z;
    q_out.z = z.z * q.w + (Sxz + Szx) * v.x + (Syz + Szy) * v.y + diag.z * v.z;

    return q_out;
}

// =========================================================
// =========================================================
//  Kernels
// =========================================================
// =========================================================

/// Common frame. Axis-aligned bounding boxes from atom positions + radii.
/// Consumes: pos, radius  |  Produces: bboxes_min, bboxes_max
__kernel void update_bboxes_rigid(
    __global const float4* curr_pos,
    __global const float*  radius,
    __global float4*       bboxes_min,
    __global float4*       bboxes_max,
    __local float4*        local_min,
    __local float4*        local_max,
    const int num_atoms
) {
    int lid = get_local_id(0);
    int gid = get_global_id(0);
    int group_id = get_group_id(0);

    float4 p = (gid < num_atoms) ? curr_pos[gid] : (float4)(0.0f);
    float  r = (gid < num_atoms) ? radius[gid]   : 0.0f;
    float invM = (gid < num_atoms) ? curr_pos[gid].w : 0.0f;
    if ((gid < num_atoms) && (invM > 1e-12f) && (r > 0.0f)) {
        local_min[lid] = (float4)(p.x - r, p.y - r, p.z - r, 0.0f);
        local_max[lid] = (float4)(p.x + r, p.y + r, p.z + r, 0.0f);
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
        bboxes_min[group_id] = local_min[0];
        bboxes_max[group_id] = local_max[0];
    }
}

inline int map_global_to_local(int t, int grp, int total_ghosts, __local const int* l_ghost_list){
    if (t < 0) return -1;
    int tgrp = t / GROUP_SIZE;
    if (tgrp == grp) return t % GROUP_SIZE;
    int found = -1;
    for (int g = 0; g < total_ghosts; g++) {
        if (l_ghost_list[g] == t) { found = GROUP_SIZE + g; break; }
    }
    return found;
}

/// Common frame. Cluster-sorted local neighbor/exclusion lists from global topology + AABBs.
/// Consumes: pos, neighs, excl, bboxes  |  Produces: neighs_local, excl_local
__kernel void build_local_topology_rigid(
    __global const float4* curr_pos,
    __global const float4* bboxes_min,
    __global const float4* bboxes_max,
    __global const int4*   neighs_global,
    __global const int4*   excl1_global,
    __global const int4*   excl2_global,
    __global int*          ghost_indices_flat,
    __global int*          ghost_counts,
    __global int4*         neighs_local,
    __global int4*         excl1_local,
    __global int4*         excl2_local,
    const int num_atoms,
    const int num_groups,
    const float margin_sq,
    const float bbox_margin
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
            my_max.y + bbox_margin >= o_min.y && my_min.y <= o_max.y + bbox_margin &&
            my_max.z + bbox_margin >= o_min.z && my_min.z <= o_max.z + bbox_margin) {
            overlap = true;
        }
        if (overlap) {
            int global_idx = other_g * GROUP_SIZE + lid;
            if (global_idx < num_atoms) {
                float4 p = curr_pos[global_idx];
                // Robust skip for padding/fixed atoms (invM<=0). Avoid relying on NaN behavior.
                if (p.w > 1e-12f) {
                    float dx = max(0.0f, max(my_min.x - p.x, p.x - my_max.x));
                    float dy = max(0.0f, max(my_min.y - p.y, p.y - my_max.y));
                    float dz = max(0.0f, max(my_min.z - p.z, p.z - my_max.z));
                    float dist_sq = dx*dx + dy*dy + dz*dz;
                    if (dist_sq < margin_sq) {
                        int slot = atomic_inc(&l_ghost_counter);
                        if (slot < MAX_GHOSTS) {
                            l_ghost_list[slot] = global_idx;
                        }
                    }
                }
            }
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    int total_ghosts = min(l_ghost_counter, MAX_GHOSTS);
    if(lid==0){ LOG_TOPOLOGY_SUMMARY(grp, total_ghosts); LOG_TOPOLOGY(grp * GROUP_SIZE + 0, lid, total_ghosts); }
    int base_offset = grp * MAX_GHOSTS;
    for (int i = lid; i < total_ghosts; i += GROUP_SIZE) {
        ghost_indices_flat[base_offset + i] = l_ghost_list[i];
        if(lid==0){ LOG_MAPPING(grp * GROUP_SIZE + 0, GROUP_SIZE + i, l_ghost_list[i], "GHOST"); }
    }
    if (lid == 0) ghost_counts[grp] = total_ghosts;
    barrier(CLK_LOCAL_MEM_FENCE);

    int my_global_id = grp * GROUP_SIZE + lid;
    if (my_global_id < num_atoms) {
        int4 ng = neighs_global[my_global_id];
        neighs_local[my_global_id] = (int4)(
            map_global_to_local(ng.x, grp, total_ghosts, l_ghost_list),
            map_global_to_local(ng.y, grp, total_ghosts, l_ghost_list),
            map_global_to_local(ng.z, grp, total_ghosts, l_ghost_list),
            map_global_to_local(ng.w, grp, total_ghosts, l_ghost_list)
        );

        int4 e1 = excl1_global[my_global_id];
        int4 e2 = excl2_global[my_global_id];
        excl1_local[my_global_id] = (int4)(
            map_global_to_local(e1.x, grp, total_ghosts, l_ghost_list),
            map_global_to_local(e1.y, grp, total_ghosts, l_ghost_list),
            map_global_to_local(e1.z, grp, total_ghosts, l_ghost_list),
            map_global_to_local(e1.w, grp, total_ghosts, l_ghost_list)
        );
        excl2_local[my_global_id] = (int4)(
            map_global_to_local(e2.x, grp, total_ghosts, l_ghost_list),
            map_global_to_local(e2.y, grp, total_ghosts, l_ghost_list),
            map_global_to_local(e2.z, grp, total_ghosts, l_ghost_list),
            map_global_to_local(e2.w, grp, total_ghosts, l_ghost_list)
        );
    }
}

/// Common frame. Jacobi positional corrections for pairwise sphere collisions.
/// Consumes: pos, radius, excl_local  |  Produces: dpos_coll
__kernel void compute_collision_cluster_rigid(
    __global const float4* pos,
    __global const float*  radius,
    __global const int4*   excl1_local,
    __global const int4*   excl2_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global float4*       dpos_coll,
    const int num_atoms,
    const float k_coll
) {
    if(ENABLE_COLL==0){
        int gid0 = get_global_id(0);
        if(gid0 < num_atoms) dpos_coll[gid0] = (float4)(0.0f);
        return;
    }
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];

    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    float invMi = (my_global_id < num_atoms) ? pi4.w : 0.0f;
    float ri = (my_global_id < num_atoms) ? radius[my_global_id] : 0.0f;
    // Never return before barriers; mask invalid lanes instead.
    if (invMi <= 1e-12f || ri <= 0.0f) {
        pi4 = (float4)(0.0f);
        invMi = 0.0f;
        ri = 0.0f;
    }
    l_pos[lid] = pi4;
    l_rad[lid] = ri;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
        l_rad[GROUP_SIZE + k] = radius[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id < num_atoms) {
        if (invMi <= 1e-12f || ri <= 0.0f) { dpos_coll[my_global_id] = (float4)(0.0f); return; }

        float3 p = l_pos[lid].xyz;
        float w_i = invMi;

        float3 sum = (float3)(0.0f);
        int4 e1 = excl1_local[my_global_id];
        int4 e2 = excl2_local[my_global_id];

        int n_ext = GROUP_SIZE + g_count;
        for (int j = 0; j < n_ext; j++) {
            if (j == lid) continue;
            int other_global = (j < GROUP_SIZE) ? (grp * GROUP_SIZE + j) : ghost_indices_flat[g_offset + (j - GROUP_SIZE)];
            float4 pj4_dbg = l_pos[j];
            float3 d_dbg = p - pj4_dbg.xyz;
            float dist_dbg = sqrt(dot(d_dbg,d_dbg));

            if (excluded8(j, e1, e2)) {
                LOG_COLLISION_CHECK(my_global_id, my_global_id, j, other_global, dist_dbg, "SKIP_EXCL");
                continue;
            }

            float4 pj4 = l_pos[j];
            float invMj = pj4.w;
            float rj = l_rad[j];
            if (invMj <= 1e-12f || rj <= 0.0f) continue;
            float3 q = pj4.xyz;
            float3 d = p - q;
            float d2 = dot(d, d);
            float rsum = ri + rj;
            float r2 = rsum * rsum;
            if (d2 < r2 && d2 > 1e-16f) {
                float dist = sqrt(d2);
                float3 n = d / dist;
                float w_j = invMj;
                float w_tot = w_i + w_j + 1e-12f;
                float dl = (rsum - dist) / w_tot;
                dl *= 0.5f;
                sum += n * (dl * w_i);
                LOG_COLLISION_CHECK(my_global_id, my_global_id, j, other_global, dist, "COLLIDE");
            } else {
                LOG_COLLISION_CHECK(my_global_id, my_global_id, j, other_global, dist_dbg, "TOO_FAR");
            }
        }
        dpos_coll[my_global_id] = (float4)(sum, 0.0f);
    }
}

/// Massfull (current). Full rigid-body XPBD: tip->atom distance with physical inertia.
/// Distributes impulse into linear (dpos_node) and angular (drot_node) recoil.
/// Consumes: pos, quat, port_local  |  Produces: dpos_node, drot_node, dpos_neigh
__kernel void compute_ports_cluster_rigid(
    __global const float4* pos,
    __global const float4* quat,
    __global const float*  radius,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float4* port_local,
    __global const float*  stiffness_flat,
    __global float4*       dpos_node,
    __global float4*       drot_node,
    __global float4*       dpos_neigh,
    const int num_atoms,
    const int nnode_per_group,
    const float dt,
    const int accumulate_dpos,
    const float rot_mass_scale
) {
    if(ENABLE_PORT==0){
        return;
    }
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];

    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    l_rad[lid] = (my_global_id < num_atoms) ? radius[my_global_id] : 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
        l_rad[GROUP_SIZE + k] = radius[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return;
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    int inode = grp * nnode_per_group + lid;

    float3 xi = pi4.xyz;
    float4 qi = quat[my_global_id];

    float mi = 1.0f / invMi;
    float invI = rot_mass_scale / (0.4f * mi + 1e-12f);
    float dt2 = dt * dt + 1e-16f;

    float3 sum_dpos = (float3)(0.0f);
    float3 sum_dtheta = (float3)(0.0f);
    if (accumulate_dpos) {
        sum_dpos = dpos_node[inode].xyz;
    }

    int4 ng = neighs_local[my_global_id];
    int* neighbors = (int*)&ng;
    int i4 = inode * 4;

    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        dpos_neigh[idx] = (float4)(0.0f);
        int jloc = neighbors[k];
        if (jloc < 0) continue;
        float K = stiffness_flat[idx];
        if (K <= 0.0f) continue;

        float invMj;
        float3 xj;
        if (jloc < (GROUP_SIZE + g_count)) {
            xj = l_pos[jloc].xyz;
            invMj = l_pos[jloc].w;
        } else {
            continue;
        }
        if (invMj <= 1e-12f) continue;

        float3 r_local = port_local[idx].xyz;
        float3 r_arm = quat_rotate(qi, r_local);
        float3 tip = xi + r_arm;

        float3 diff = xj - tip;
        float dist2 = dot(diff, diff);
        if (dist2 < 1e-16f) continue;
        float dist = sqrt(dist2);
        float3 n = diff / dist;

        float w_i = invMi;
        float w_j = invMj;
        float3 rxn = cross(r_arm, n);
        float w_ang = dot(rxn, rxn) * invI;
        float alpha = 1.0f / (K * dt2);

        // Determine neighbor type early to adjust weighting
        int j_global = (jloc < GROUP_SIZE) ? (grp * GROUP_SIZE + jloc) : ghost_indices_flat[g_offset + (jloc - GROUP_SIZE)];
        int j_lid = j_global & (GROUP_SIZE - 1);
        int j_isnode = (j_lid < nnode_per_group);

        float w_total = w_i + w_j + w_ang + alpha + 1e-12f;
        // If neighbor is also a node, it also rotates. Assume symmetric inertia (w_ang_j ~ w_ang_i) 
        // to avoid over-stiffening the constraint.
        if (j_isnode) { w_total += w_ang; }

        float impulse_mag = dist / w_total;

        // If neighbor is also a node, the constraint is evaluated by both node threads.
        // Therefore BOTH linear impulse and torque contribution must be halved to avoid double-counting.
        if (j_isnode) { impulse_mag *= 0.5f; }

        float3 P = n * impulse_mag;

        sum_dpos   += P * w_i;
        sum_dtheta += cross(r_arm, P) * invI;
        dpos_neigh[idx] = (float4)(-P * w_j, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos, 0.0f);
    drot_node[inode] = (float4)(sum_dtheta, 0.0f);
}


// compute_ports_cluster_rigid before changes implemented by Gemini 
/// Massfull (orig). Same physics as 'current' but without rot_mass_scale tuning.
/// Consumes: pos, quat, port_local  |  Produces: dpos_node, drot_node, dpos_neigh
__kernel void compute_ports_cluster_rigid_orig(
    __global const float4* pos,
    __global const float4* quat,
    __global const float*  radius,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float4* port_local,
    __global const float*  stiffness_flat,
    __global float4*       dpos_node,
    __global float4*       drot_node,
    __global float4*       dpos_neigh,
    const int num_atoms,
    const int nnode_per_group,
    const float dt,
    const int accumulate_dpos,
    __global const float4* quat_opt,
    const int skip_rotation
) {
    if(ENABLE_PORT==0){
        return;
    }
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    __local float  l_rad[GROUP_SIZE + MAX_GHOSTS];

    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    l_rad[lid] = (my_global_id < num_atoms) ? radius[my_global_id] : 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
        l_rad[GROUP_SIZE + k] = radius[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return;
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    int inode = grp * nnode_per_group + lid;

    float3 xi = pi4.xyz;
    float4 qi = (quat_opt != (__global const float4*)0) ? quat_opt[my_global_id] : quat[my_global_id];

    float mi = 1.0f / invMi;
    float invI = 1.0f / (0.4f * mi + 1e-12f);
    float dt2 = dt * dt + 1e-16f;

    float3 sum_dpos = (float3)(0.0f);
    float3 sum_dtheta = (float3)(0.0f);
    if (accumulate_dpos) {
        sum_dpos = dpos_node[inode].xyz;
    }

    int4 ng = neighs_local[my_global_id];
    int* neighbors = (int*)&ng;
    int i4 = inode * 4;

    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        dpos_neigh[idx] = (float4)(0.0f);
        int jloc = neighbors[k];
        if (jloc < 0) continue;
        float K = stiffness_flat[idx];
        if (K <= 0.0f) continue;

        float invMj;
        float3 xj;
        if (jloc < (GROUP_SIZE + g_count)) {
            xj = l_pos[jloc].xyz;
            invMj = l_pos[jloc].w;
        } else {
            continue;
        }
        if (invMj <= 1e-12f) continue;

        float3 r_local = port_local[idx].xyz;
        float3 r_arm = quat_rotate(qi, r_local);
        float3 tip = xi + r_arm;

        float3 diff = xj - tip;
        float dist2 = dot(diff, diff);
        if (dist2 < 1e-16f) continue;
        float dist = sqrt(dist2);
        float3 n = diff / dist;

        float w_i = invMi;
        float w_j = invMj;
        float3 rxn = cross(r_arm, n);
        float w_ang = dot(rxn, rxn) * invI;
        float alpha = 1.0f / (K * dt2);
        float w_total = w_i + w_j + w_ang + alpha + 1e-12f;
        float impulse_mag = dist / w_total;
        // Only halve if the same constraint will also be evaluated by the neighbor thread.
        // For node-cap bonds, only node threads run => no double counting.
        int j_global = (jloc < GROUP_SIZE) ? (grp * GROUP_SIZE + jloc) : ghost_indices_flat[g_offset + (jloc - GROUP_SIZE)];
        int j_lid = j_global & (GROUP_SIZE - 1);
        int j_isnode = (j_lid < nnode_per_group);
        if (j_isnode) { impulse_mag *= 0.5f; }

        float3 P = n * impulse_mag;
        sum_dpos += P * w_i;
        if (!skip_rotation) {
            sum_dtheta += cross(r_arm, P) * invI;
        }
        dpos_neigh[idx] = (float4)(-P * w_j, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos, 0.0f);
    if (!skip_rotation) {
        drot_node[inode] = (float4)(sum_dtheta, 0.0f);
    }
}

/// Massless. Iterative Newton-Raphson in omega-space for port alignment, then
/// central linear recoil along the center-center line.
/// Consumes: pos, quat, port_local  |  Produces: dpos_node, drot_node, dpos_neigh
__kernel void compute_ports_cluster_rigid_substep_optimized(
    __global const float4* pos,
    __global const float4* quat,
    __global const float*  radius,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float4* port_local,
    __global const float*  stiffness_flat,
    __global float4*       dpos_node,
    __global float4*       drot_node,
    __global float4*       dpos_neigh,
    const int num_atoms,
    const int nnode_per_group,
    const float dt,
    const int accumulate_dpos,
    const int n_rot_substeps,
    const float rot_eps,
    const float theta_max
) {
    if(ENABLE_PORT==0) return;

    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    // --- 1. LOCAL MEMORY LOAD ---
    // (We reuse the existing logic, which is efficient)
    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    
    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- 2. EARLY EXITS ---
    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return; 
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    // --- 3. PRELOAD DATA TO REGISTERS (Essential for Speed) ---
    int inode = grp * nnode_per_group + lid;
    int4 ng = neighs_local[my_global_id];
    int neighbors[4] = {ng.x, ng.y, ng.z, ng.w};
    
    // Cache for Ports and Stiffness to avoid Global Reads in loop
    float3 r_local_cache[4];
    float  K_cache[4];
    int    valid_bond[4]; // 0 or 1
    
    int i4 = inode * 4;
    
    for(int k=0; k<4; k++) {
        int idx = i4 + k;
        int jloc = neighbors[k];
        
        if (jloc >= 0 && stiffness_flat[idx] > 0.0f) {
             float K = stiffness_flat[idx];
             r_local_cache[k] = port_local[idx].xyz;
             K_cache[k]       = K;
             valid_bond[k]    = 1;
        } else {
            valid_bond[k] = 0;
        }
    }

    // Initialize State Registers
    float3 xi = pi4.xyz;
    float4 qi = quat[my_global_id];
    float4 qi_orig = qi;

    float mi = 1.0f / invMi;
    float invI = 1.0f / (0.4f * mi + 1e-12f); // Approx sphere inertia
    float dt2 = dt * dt + 1e-16f;

    float3 sum_dpos = (float3)(0.0f);
    if (accumulate_dpos) sum_dpos = dpos_node[inode].xyz;

    for (int sub=0; sub < n_rot_substeps; sub++) {
        // Newton step in rotation-vector space: H * omega = torque
        float3 torque = (float3)(0.0f);
        float3 H0 = (float3)(0.0f);
        float3 H1 = (float3)(0.0f);
        float3 H2 = (float3)(0.0f);

        for (int k = 0; k < 4; k++) {
            if (!valid_bond[k]) continue;
            int jloc = neighbors[k];
            float3 xj = l_pos[jloc].xyz;

            float3 r_arm = q_rot(qi, r_local_cache[k]);
            float3 diff  = xj - (xi + r_arm);

            // Weight: using stiffness only (no mass coupling)
            float w = K_cache[k];

            torque += cross(r_arm, diff) * w;

            add_hessian_rr(&H0, &H1, &H2, w, r_arm);
        }

        float damp = (rot_eps > 0.0f) ? rot_eps : 0.0f;
        H0.x += (1e-8f + damp);
        H1.y += (1e-8f + damp);
        H2.z += (1e-8f + damp);

        float3 omega = solve_3x3(H0, H1, H2, torque);

        if (theta_max > 0.0f) {
            float a = length(omega);
            if (a > theta_max) omega *= (theta_max / (a + 1e-16f));
        }

        qi = apply_delta_rot(qi, omega);
    }

    float3 sum_dpos_final = (float3)(0.0f);
    if (accumulate_dpos) sum_dpos_final = dpos_node[inode].xyz;

    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        dpos_neigh[idx] = (float4)(0.0f); // Reset
        if (!valid_bond[k]) continue;

        int jloc = neighbors[k];
        float3 xj = l_pos[jloc].xyz;
        float invMj = l_pos[jloc].w;

        // Recalculate with optimized qi
        float3 r_arm = q_rot(qi, r_local_cache[k]);
        float3 diff  = xj - (xi + r_arm);

        float3 rij = xj - xi;
        float r2 = dot(rij, rij);
        if (r2 < 1e-16f) continue;
        float3 n = rij * rsqrt(r2);
        float d = dot(diff, n);
        if (fabs(d) < 1e-12f) continue;

        float3 rxn = cross(r_arm, n);
        float w_ang = dot(rxn, rxn) * invI;
        float alpha = 1.0f / (K_cache[k] * dt2);

        // Full XPBD denominator
        float w_total = invMi + invMj + w_ang + alpha + 1e-12f;
        float impulse = d / w_total;
        
        // Double counting prevention (same as original logic)
        int j_global = (jloc < GROUP_SIZE) ? (grp * GROUP_SIZE + jloc) : ghost_indices_flat[g_offset + (jloc - GROUP_SIZE)];
        int j_lid = j_global & (GROUP_SIZE - 1);
        if (j_lid < nnode_per_group) { impulse *= 0.5f; }

        float3 P = n * impulse;
        sum_dpos_final += P * invMi;
        dpos_neigh[idx] = (float4)(-P * invMj, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos_final, 0.0f);
    
    drot_node[inode] = (float4)(quat_delta_rotvec(qi, qi_orig), 0.0f);
}

/// Massless. Kabsch/Polar-decomposition orientation solve from covariance matrix,
/// then central linear recoil along the center-center line.
/// Consumes: pos, quat, port_local  |  Produces: dpos_node, drot_node, dpos_neigh
__kernel void compute_ports_cluster_rigid_shapematch(
    __global const float4* pos,
    __global const float4* quat,
    __global const float*  radius,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float4* port_local,
    __global const float*  stiffness_flat,
    __global float4*       dpos_node,
    __global float4*       drot_node,
    __global float4*       dpos_neigh,
    const int num_atoms,
    const int nnode_per_group,
    const float dt,
    const int accumulate_dpos
) {
    if(ENABLE_PORT==0) return;

    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    // --- 1. LOCAL MEMORY SETUP (Standard) ---
    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- 2. EARLY EXITS & PRELOAD ---
    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return; 
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    // Load Indices
    int inode = grp * nnode_per_group + lid;
    int4 ng = neighs_local[my_global_id];
    int neighbors[4] = {ng.x, ng.y, ng.z, ng.w};
    int i4 = inode * 4;

    // Preload Data to Registers
    float3 r_local_cache[4];
    float  K_cache[4];
    int    valid_bond[4];
    
    for(int k=0; k<4; k++) {
        int idx = i4 + k;
        int jloc = neighbors[k];
        
        if (jloc >= 0 && stiffness_flat[idx] > 0.0f) {
             float K = stiffness_flat[idx];
             r_local_cache[k] = port_local[idx].xyz;
             K_cache[k]       = K;
             valid_bond[k]    = 1;
        } else {
             valid_bond[k] = 0;
        }
    }

    float3 xi = pi4.xyz;
    float4 qi_orig = quat[my_global_id];

    // --- 3. EXPLICIT SHAPE MATCHING (Fixed Pivot) ---
    // We want Rotation R that minimizes sum( w * || R*r_local - (x_j - x_i) ||^2 )
    // This is solved by Polar Decomposition of Covariance Matrix A.
    // A = Sum( w * (x_j - x_i) * r_local^T )

    float A[9] = {0.0f}; // Covariance Matrix
    
    for(int k=0; k<4; k++) {
        if(!valid_bond[k]) continue;
        
        int jloc = neighbors[k];
        float3 xj = l_pos[jloc].xyz;
         float3 vec_target = xj - xi; // Target vector relative to node center
         float3 vec_source = r_local_cache[k];
         float w = K_cache[k]; // Use stiffness as weight

        mat3_outer_add(A, w, vec_target, vec_source);
     }

    float nF = 0.0f;
    for(int i=0; i<9; i++) nF += A[i]*A[i];
    float scale = (nF > 1e-12f) ? native_rsqrt(nF) * 1.7f : 1.0f; // 1.7 approx sqrt(3)

    float R[9];
    for(int i=0; i<9; i++) R[i] = A[i] * scale;
    
    float Rt[9], M[9], T[9];
    // 4 Iterations is usually sufficient for single precision
    for(int iter=0; iter<4; iter++) {
        // Transpose R -> Rt
        Rt[0]=R[0]; Rt[1]=R[3]; Rt[2]=R[6];
        Rt[3]=R[1]; Rt[4]=R[4]; Rt[5]=R[7];
        Rt[6]=R[2]; Rt[7]=R[5]; Rt[8]=R[8];
        
        // M = Rt * R
        mat3_mul(Rt, R, M);
        
        // T = 3I - M
        for(int i=0; i<9; i++) T[i] = -M[i];
        T[0] += 3.0f; T[4] += 3.0f; T[8] += 3.0f;
        
        // R_new = 0.5 * R * T
        mat3_mul(R, T, M); // Use M as temp
        for(int i=0; i<9; i++) R[i] = 0.5f * M[i];
    }
    
    float4 qi_opt = mat3_to_quat(R);

    // --- 5. LINEAR RECOIL (Using Exact Rotation) ---
    float3 sum_dpos = (float3)(0.0f);
    if (accumulate_dpos) sum_dpos = dpos_node[inode].xyz;
    
    float invI = 1.0f / (0.4f * (1.0f/invMi) + 1e-12f);
    float dt2 = dt*dt;

    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        dpos_neigh[idx] = (float4)(0.0f);
        if (!valid_bond[k]) continue;

        int jloc = neighbors[k];
        float3 xj = l_pos[jloc].xyz;
        float invMj = l_pos[jloc].w;

        // Use qi_opt
        float3 r_arm = q_rot(qi_opt, r_local_cache[k]);

        float3 rij = xj - xi;
        float r2 = dot(rij, rij);
        if (r2 < 1e-16f) continue;
        float3 n = rij * rsqrt(r2);
        float3 diff = xj - (xi + r_arm);
        float d = dot(diff, n);
        if (fabs(d) < 1e-12f) continue;

        // Same XPBD weighting as substep kernel for fair comparison
        float3 rxn = cross(r_arm, n);
        float w_ang = dot(rxn, rxn) * invI;
        float alpha = 1.0f / (K_cache[k] * dt2);

        float w_total = invMi + invMj + w_ang + alpha + 1e-12f;
        float impulse = d / w_total;
        
        int j_global = (jloc < GROUP_SIZE) ? (grp * GROUP_SIZE + jloc) : ghost_indices_flat[g_offset + (jloc - GROUP_SIZE)];
        int j_lid = j_global & (GROUP_SIZE - 1);
        if (j_lid < nnode_per_group) { impulse *= 0.5f; }

        float3 P = n * impulse;
        sum_dpos += P * invMi;
        dpos_neigh[idx] = (float4)(-P * invMj, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos, 0.0f);

    drot_node[inode] = (float4)(quat_delta_rotvec(qi_opt, qi_orig), 0.0f);
}

/// Massless Pass-1. Davenport q-method: 4x4 eigenproblem from weighted tip-neighbor vectors.
/// Writes optimal quaternion delta to drot_node and quat_opt into temp buffer (dquat_mom).
/// Consumes: pos, quat, port_local  |  Produces: drot_node, quat_opt(dquat_mom)
__kernel void compute_optimal_rotation_eigen(
    __global const float4* pos,
    __global const float4* quat,
    __global const float*  radius,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float4* port_local,
    __global const float*  stiffness_flat,
    __global float4*       drot_node,
    __global float4*       quat_opt,
    const int num_atoms,
    const int nnode_per_group,
    const float dt
) {
    if(ENABLE_PORT==0) return;

    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return;
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    int inode = grp * nnode_per_group + lid;
    int4 ng = neighs_local[my_global_id];
    int neighbors[4] = {ng.x, ng.y, ng.z, ng.w};
    int i4 = inode * 4;

    float3 r_local_cache[4];
    float  K_cache[4];
    int    valid_bond[4];

    for(int k=0; k<4; k++) {
        int idx = i4 + k;
        int jloc = neighbors[k];
        if (jloc >= 0 && stiffness_flat[idx] > 0.0f) {
             float K = stiffness_flat[idx];
             r_local_cache[k] = port_local[idx].xyz;
             K_cache[k]       = K;
             valid_bond[k]    = 1;
        } else {
             valid_bond[k] = 0;
        }
    }

    float3 centroid_neigh = (float3)(0.0f);
    float sum_w = 0.0f;
    for(int k=0; k<4; k++) {
        if(!valid_bond[k]) continue;
        int jloc = neighbors[k];
        float3 xj = l_pos[jloc].xyz;
        float w = K_cache[k];
        centroid_neigh += xj * w;
        sum_w += w;
    }
    if (sum_w > 1e-9f) centroid_neigh *= (1.0f / sum_w);
    else centroid_neigh = pi4.xyz;

    float S[9] = {0.0f};
    for(int k=0; k<4; k++) {
        if(!valid_bond[k]) continue;
        float w = K_cache[k];
        float3 p = r_local_cache[k];
        float3 n = l_pos[neighbors[k]].xyz - centroid_neigh;
        mat3_outer_add(S, w, n, p);
    }

    float4 q = quat[my_global_id];
    for(int iter=0; iter<4; iter++) {
        q = apply_K_matrix(S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], q);
        q = normalize(q);
    }

    float4 qi_opt = q;
    float4 qi_orig = quat[my_global_id];

    drot_node[inode] = (float4)(quat_delta_rotvec(qi_opt, qi_orig), 0.0f);
    quat_opt[my_global_id] = qi_opt;
}

/// Legacy single-pass eigen (deprecated). Superseded by compute_ports_cluster_rigid_eigen_tips.
__kernel void compute_ports_cluster_rigid_eigen(
    __global const float4* pos,
    __global const float4* quat,
    __global const float*  radius,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float4* port_local,
    __global const float*  stiffness_flat,
    __global float4*       dpos_node,
    __global float4*       drot_node,
    __global float4*       dpos_neigh,
    const int num_atoms,
    const int nnode_per_group,
    const float dt,
    const int accumulate_dpos
) {
    if(ENABLE_PORT==0) return;

    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    // --- 1. LOCAL MEMORY SETUP ---
    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- 2. EARLY EXITS & PRELOAD ---
    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return; 
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    int inode = grp * nnode_per_group + lid;
    int4 ng = neighs_local[my_global_id];
    int neighbors[4] = {ng.x, ng.y, ng.z, ng.w};
    int i4 = inode * 4;

    float3 r_local_cache[4];
    float  K_cache[4];
    int    valid_bond[4];
    
    for(int k=0; k<4; k++) {
        int idx = i4 + k;
        int jloc = neighbors[k];
        
        if (jloc >= 0 && stiffness_flat[idx] > 0.0f) {
             float K = stiffness_flat[idx];
             r_local_cache[k] = port_local[idx].xyz;
             K_cache[k]       = K;
             valid_bond[k]    = 1;
        } else {
             valid_bond[k] = 0;
        }
    }

    float3 centroid_neigh = (float3)(0.0f);
    float sum_w = 0.0f;
    
    for(int k=0; k<4; k++) {
        if(!valid_bond[k]) continue;
        
        int jloc = neighbors[k];
        float3 xj = l_pos[jloc].xyz;
         float w = K_cache[k];
         centroid_neigh += xj * w;
         sum_w += w;
    }
    
    if (sum_w > 1e-9f) centroid_neigh *= (1.0f / sum_w);
    else centroid_neigh = pi4.xyz; // Fallback

    // --- 4. ACCUMULATE COVARIANCE MATRIX B ---
    // B = Sum( w_i * (n_i - centroid) * p_i^T )
    // Note: p_i (ports) are already local (centered at 0,0,0).
     float S[9] = {0.0f};

    for(int k=0; k<4; k++) {
        if(!valid_bond[k]) continue;
        float w = K_cache[k];
        float3 p = r_local_cache[k];
        float3 n = l_pos[neighbors[k]].xyz - centroid_neigh; // Centered target
        mat3_outer_add(S, w, n, p);
     }

    float4 q = quat[my_global_id];
    
    // 4 Iterations is robust for standard physics
    // Each step: q = K * q; q = normalize(q);
    for(int iter=0; iter<4; iter++) {
        q = apply_K_matrix(S[0], S[1], S[2], S[3], S[4], S[5], S[6], S[7], S[8], q);
        q = normalize(q);
    }
    
    float4 qi_opt = q;
    float4 qi_orig = quat[my_global_id];
    float3 xi = pi4.xyz;

    // --- 6. LINEAR RECOIL ---
    // Same linear logic as before, using the optimal orientation

    float3 sum_dpos = (float3)(0.0f);
    if (accumulate_dpos) sum_dpos = dpos_node[inode].xyz;
    
    float invI = 1.0f / (0.4f * (1.0f/invMi) + 1e-12f);
    float dt2 = dt*dt;

    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        dpos_neigh[idx] = (float4)(0.0f);
        if (!valid_bond[k]) continue;

        int jloc = neighbors[k];
        float3 xj = l_pos[jloc].xyz;
        float invMj = l_pos[jloc].w;

        float3 r_arm = q_rot(qi_opt, r_local_cache[k]);
        
        float3 diff = xj - (xi + r_arm);
        float dist2 = dot(diff, diff);
        if (dist2 < 1e-16f) continue;
        float dist = sqrt(dist2);
        float3 n = diff / dist;

        float3 rxn = cross(r_arm, n);
        float w_ang = dot(rxn, rxn) * invI;
        float alpha = 1.0f / (K_cache[k] * dt2);
        
        float w_total = invMi + invMj + w_ang + alpha + 1e-12f;
        float impulse = dist / w_total;
        
        // Double-count check
        int j_global = (jloc < GROUP_SIZE) ? (grp * GROUP_SIZE + jloc) : ghost_indices_flat[g_offset + (jloc - GROUP_SIZE)];
        int j_lid = j_global & (GROUP_SIZE - 1);
        if (j_lid < nnode_per_group) { impulse *= 0.5f; }

        float3 P = n * impulse;
        sum_dpos += P * invMi;
        dpos_neigh[idx] = (float4)(-P * invMj, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos, 0.0f);

    drot_node[inode] = (float4)(quat_delta_rotvec(qi_opt, qi_orig), 0.0f);
}

/// Massless Pass-2. Reads precomputed tip positions, computes symmetric linear recoil
/// with impulses constrained to the center-center line.
/// Consumes: pos, tips  |  Produces: dpos_node, dpos_neigh
__kernel void compute_ports_cluster_rigid_eigen_tips(
    __global const float4* pos,
    __global const int4*   neighs_local,
    __global const int*    ghost_indices_flat,
    __global const int*    ghost_counts,
    __global const float*  stiffness_flat,
    __global const float4* tips,
    __global float4*       dpos_node,
    __global float4*       dpos_neigh,
    const int num_atoms,
    const int nnode_per_group,
    const float dt
) {
    if(ENABLE_PORT==0) return;

    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * GROUP_SIZE + lid;

    __local float4 l_pos[GROUP_SIZE + MAX_GHOSTS];
    float4 pi4 = (my_global_id < num_atoms) ? pos[my_global_id] : (float4)(0.0f);
    l_pos[lid] = pi4;
    barrier(CLK_LOCAL_MEM_FENCE);

    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    for (int k = lid; k < g_count; k += GROUP_SIZE) {
        int gid = ghost_indices_flat[g_offset + k];
        l_pos[GROUP_SIZE + k] = pos[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (my_global_id >= num_atoms) return;
    if (lid >= nnode_per_group) return;
    float invMi = pi4.w;
    if (invMi <= 1e-12f) return;

    int inode = grp * nnode_per_group + lid;
    int4 ng = neighs_local[my_global_id];
    int* neighbors = (int*)&ng;
    int i4 = inode * 4;

    float3 sum_dpos = (float3)(0.0f);
    float dt2 = dt * dt + 1e-16f;

    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        dpos_neigh[idx] = (float4)(0.0f);
        int jloc = neighbors[k];
        if (jloc < 0) continue;
        float K = stiffness_flat[idx];
        if (K <= 0.0f) continue;

        float invMj;
        float3 xj;
        if (jloc < (GROUP_SIZE + g_count)) {
            xj = l_pos[jloc].xyz;
            invMj = l_pos[jloc].w;
        } else {
            continue;
        }
        if (invMj <= 1e-12f) continue;

        float3 my_tip = tips[idx].xyz;
        int j_global = (jloc < GROUP_SIZE) ? (grp * GROUP_SIZE + jloc) : ghost_indices_flat[g_offset + (jloc - GROUP_SIZE)];
        int j_lid = j_global & (GROUP_SIZE - 1);
        int j_isnode = (j_lid < nnode_per_group);

        float3 rij = xj - pi4.xyz;
        float r2 = dot(rij, rij);
        if (r2 < 1e-16f) continue;
        float3 n = rij * rsqrt(r2);
        float d = dot(xj - my_tip, n);
        if (fabs(d) < 1e-12f) continue;

        float w_i = invMi;
        float w_j = invMj;
        float alpha = 1.0f / (K * dt2);
        float w_total = w_i + w_j + alpha + 1e-12f;
        float impulse_mag = d / w_total;
        if (j_isnode) { impulse_mag *= 0.5f; }

        float3 P = n * impulse_mag;
        sum_dpos += P * w_i;
        dpos_neigh[idx] = (float4)(-P * w_j, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos, 0.0f);
}

/// Helper. Rotates each node's local ports into world-space tip positions.
/// Consumes: pos, quat, port_local  |  Produces: tips
__kernel void compute_tips(
    const int natoms,
    const int nnode_per_group,
    __global const float4* pos,
    __global const float4* quat,
    __global const float4* port_local,
    __global float4* tips
) {
    int i = get_global_id(0);
    if (i >= natoms) return;
    float invMi = pos[i].w;
    if (invMi <= 1e-12f) return;
    int lid = i & (GROUP_SIZE - 1);
    int grp = i / GROUP_SIZE;
    int isnode = (lid < nnode_per_group);
    if (!isnode) return;
    int inode = grp * nnode_per_group + lid;
    int i4 = inode * 4;
    float3 xi = pos[i].xyz;
    float4 qi = quat[i];
    for (int k = 0; k < 4; k++) {
        int idx = i4 + k;
        float3 r_local = port_local[idx].xyz;
        float3 r_arm = quat_rotate(qi, r_local);
        tips[idx] = (float4)(xi + r_arm, 0.0f);
    }
}

/// Integrator. Applies Jacobi/XPBD corrections to pos/quat, optionally with
/// Heavy-Ball momentum. For massless modes, also writes updated tips.
/// Consumes: pos, quat, dpos_node, drot_node, dpos_neigh, dpos_coll
/// Produces: pos, quat, dpos_mom, dquat_mom[, tips]
__kernel void apply_corrections_rigid_ports(
    const int natoms,
    const int nnode_per_group,
    __global float4* pos,
    __global float4* quat,
    __global const int* fixmask,
    __global const int4* bkSlots,
    __global const float4* dpos_node,
    __global const float4* drot_node,
    __global const float4* dpos_neigh,
    __global const float4* dpos_coll,
    __global float4* dpos_mom,
    __global float4* dquat_mom,
    const float relaxation,
    const float beta,
    const int massless_rot,
    __global const float4* port_local,
    __global float4* tips
) {
    int i = get_global_id(0);
    if (i >= natoms) return;

    float invMi = pos[i].w;
    if (invMi <= 1e-12f) return; // padding/fixed atoms are not updated

    int lid = i & (GROUP_SIZE - 1);
    int grp = i / GROUP_SIZE;
    int inode = grp * nnode_per_group + lid;
    int isnode = (lid < nnode_per_group);

    float3 dx_coll = dpos_coll[i].xyz;
    float3 dx_port = (float3)(0.0f);
    if (isnode) {
        dx_port += dpos_node[inode].xyz;
    }

    int4 bk = bkSlots[i];
    int* pbk = (int*)&bk;
    for (int k = 0; k < 4; k++) {
        int slot = pbk[k];
        if (slot >= 0) {
            dx_port += dpos_neigh[slot].xyz;
        }
    }

    float3 dx = dx_coll + dx_port;

    int msk = (fixmask != 0) ? fixmask[i] : 0;
    if (msk & 1) dx.x = 0.0f;
    if (msk & 2) dx.y = 0.0f;
    if (msk & 4) dx.z = 0.0f;

    float3 d_mom = dpos_mom[i].xyz;
    float3 move = dx * relaxation + d_mom * beta;
    
    if (msk & 8) { // constrain Z to 0 plane if requested (e.g. 2D mode)
        float3 p_old = pos[i].xyz;
        if (fabs(p_old.z) > 1e-6f) { // snap to plane if drifted
             move.z -= p_old.z;
        } else {
             move.z = 0.0f;
        }
    }
    
    pos[i].xyz += move;
    dpos_mom[i] = (float4)(move, 0.0f);

    if (isnode) {
        if (msk & (1|2|4)) return; // pinned node: do not rotate
        
        float3 dtheta = drot_node[inode].xyz * relaxation;
        float angle = length(dtheta);
        
        float4 q_old = quat[i];
        float4 q_jacobi = q_old;
        
        if (angle > 1e-8f) {
            float3 axis = dtheta / angle;
            float4 dq = quat_from_axis_angle(axis, angle);
            q_jacobi = quat_normalize(quat_mul(dq, q_old));
        }
        
        float4 q_new;
        if (massless_rot) {
            q_new = q_jacobi;
            dquat_mom[i] = (float4)(0.0f);
        } else {
            float4 dq_mom = dquat_mom[i];
            q_new = q_jacobi + dq_mom * beta;
            q_new = quat_normalize(q_new);
        }
        
        if (dot(q_new, q_old) < 0.0f) q_new = -q_new;
        
        quat[i] = q_new;
        if (!massless_rot) dquat_mom[i] = q_new - q_old;

        if (massless_rot && tips != (__global float4*)0 && port_local != (__global const float4*)0) {
            int i4 = inode * 4;
            float3 xi = pos[i].xyz;
            for (int k = 0; k < 4; k++) {
                int idx = i4 + k;
                float3 r_local = port_local[idx].xyz;
                float3 r_arm = quat_rotate(q_new, r_local);
                tips[idx] = (float4)(xi + r_arm, 0.0f);
            }
        }
    }
}