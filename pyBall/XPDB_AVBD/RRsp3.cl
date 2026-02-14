// RRsp3.cl - Cluster-sorted rigid ports + collisions (PBD/Jacobi) with recoils


// Only enable prints if this flag is set during compilation
#ifdef ENABLE_DEBUG_PRINTS

// Filter: Only print for specific atoms to save buffer space
#define DEBUG_GID_START 0
#define DEBUG_GID_END   8 

#define LOG_TOPOLOGY(gid, lid, n_ghosts)  if (gid >= DEBUG_GID_START && gid < DEBUG_GID_END) printf("TOPOLOGY: GID=%d LID=%d n_ghosts=%d\n", gid, lid, n_ghosts);
#define LOG_MAPPING(gid, local_idx, mapped_global_idx, type)  if (gid >= DEBUG_GID_START && gid < DEBUG_GID_END) printf("MAP: GID=%d L_IDX=%d G_IDX=%d TYPE=%s\n", gid, local_idx, mapped_global_idx, type);
#define LOG_COLLISION_CHECK(gid, my_global, other_local, other_global, dist, action) if (gid >= DEBUG_GID_START && gid < DEBUG_GID_END) printf("COLL: MeG=%d OtherL=%d OtherG=%d Dist=%.4f Action=%s\n", my_global, other_local, other_global, dist, action);

#else
// Empty macros for production
#define LOG_TOPOLOGY(gid, lid, n_ghosts)
#define LOG_MAPPING(gid, local_idx, mapped_global_idx, type)
#define LOG_COLLISION_CHECK(gid, my_global, other_local, other_global, dist, action)
#endif

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define GROUP_SIZE     64
#define MAX_GHOSTS     128

#ifndef ENABLE_COLL
#define ENABLE_COLL 1
#endif

#ifndef ENABLE_PORT
#define ENABLE_PORT 1
#endif

// ------------------------------------------------------------------
// HELPER: Bounding Box Intersection
// ------------------------------------------------------------------
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

inline float4 quat_normalize(float4 q){
    float n2 = dot(q,q);
    if(n2<1e-16f) return (float4)(0.0f,0.0f,0.0f,1.0f);
    return q * rsqrt(n2);
}

inline int excluded8(int j, int4 a, int4 b){
    if( (j==a.x) || (j==a.y) || (j==a.z) || (j==a.w) ) return 1;
    if( (j==b.x) || (j==b.y) || (j==b.z) || (j==b.w) ) return 1;
    return 0;
}

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
    if(lid==0){ LOG_TOPOLOGY(grp * GROUP_SIZE + 0, lid, total_ghosts); }
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
    const int accumulate_dpos
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
        sum_dtheta += cross(r_arm, P) * invI;
        dpos_neigh[idx] = (float4)(-P * w_j, 0.0f);
    }

    dpos_node[inode] = (float4)(sum_dpos, 0.0f);
    drot_node[inode] = (float4)(sum_dtheta, 0.0f);
}

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
    const float relaxation
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

    // Port constraints include explicit recoils (dpos_neigh), so we do not apply any extra scaling for caps.
    float3 dx = dx_coll + dx_port;

    int msk = (fixmask != 0) ? fixmask[i] : 0;
    if (msk & 1) dx.x = 0.0f;
    if (msk & 2) dx.y = 0.0f;
    if (msk & 4) dx.z = 0.0f;

    float3 xi = pos[i].xyz;
    xi += dx * relaxation;
    if (msk & 8) xi.z = 0.0f;
    pos[i].xyz = xi;

    if (isnode) {
        if (msk & (1|2|4)) return; // pinned node: do not rotate
        float3 dtheta = drot_node[inode].xyz * relaxation;
        float angle = length(dtheta);
        if (angle > 1e-8f) {
            float3 axis = dtheta / angle;
            float4 dq = quat_from_axis_angle(axis, angle);
            quat[i] = quat_normalize(quat_mul(dq, quat[i]));
        }
    }
}
