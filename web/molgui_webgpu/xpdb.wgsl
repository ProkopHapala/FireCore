// CONSTANTS (Must match JS configuration)
const GROUP_SIZE: u32 = 64;
const MAX_NEIGH_COLL: u32 = 64;
const N_MAX_BONDED: u32 = 16;
const MAX_GHOSTS: u32 = 128;

struct Params {
    num_atoms: u32,
    num_groups: u32,
    inner_iters: u32,
    _pad1: u32, // Alignment padding
    dt: f32,
    k_coll: f32,
    omega: f32,
    momentum_beta: f32,
    coll_margin_sq: f32,
    bbox_margin: f32,
    bbox_margin_sq: f32, // unused but keeps alignment
    _pad2: f32,
};

// Bindings
@group(0) @binding(0) var<uniform> u_params: Params;
@group(0) @binding(1) var<storage, read_write> curr_pos: array<vec4f>;
@group(0) @binding(2) var<storage, read_write> pred_pos: array<vec4f>;
@group(0) @binding(3) var<storage, read> atom_params: array<vec4f>; // .x=radius, .w=mass

@group(0) @binding(4) var<storage, read_write> bboxes: array<vec4f>; // packed [min,max] per group (2 vec4 per group)
@group(0) @binding(5) var<storage, read_write> ghost_packed: array<i32>; // per group: [indices... count] length=(MAX_GHOSTS+1)

// Topology (Fixed size arrays flattened)
@group(0) @binding(6) var<storage, read> bond_indices_global: array<i32>;
@group(0) @binding(7) var<storage, read_write> bond_indices_local: array<i32>;
@group(0) @binding(8) var<storage, read> bond_len_stiff: array<vec2f>; // packed (L, K)

// Shared Memory
var<workgroup> local_min: array<vec4f, GROUP_SIZE>;
var<workgroup> local_max: array<vec4f, GROUP_SIZE>;
var<workgroup> l_ghost_list: array<i32, MAX_GHOSTS>;
var<workgroup> l_ghost_counter: atomic<i32>;
var<workgroup> l_my_bbox_min: vec4f;
var<workgroup> l_my_bbox_max: vec4f;

// Solver Shared Memory
// Need to combine size: GROUP_SIZE + MAX_GHOSTS. WGSL requires const size for arrays.
// 64 + 128 = 192
var<workgroup> l_pos: array<vec4f, 192>; 
var<workgroup> l_pos_prev: array<vec4f, GROUP_SIZE>;
var<workgroup> l_pos_new: array<vec4f, GROUP_SIZE>;
var<workgroup> l_rad: array<f32, 192>;

// ------------------------------------------------------------------
// KERNEL 1: Update Bounding Boxes
// ------------------------------------------------------------------
@compute @workgroup_size(64, 1, 1) // GROUP_SIZE
fn update_bboxes(
    @builtin(local_invocation_id) lid: vec3u,
    @builtin(global_invocation_id) gid: vec3u,
    @builtin(workgroup_id) grp: vec3u
) {
    let l_idx = lid.x;
    let g_idx = gid.x;
    let group_id = grp.x;

    // 1. Load into local memory
    if (g_idx < u_params.num_atoms) {
        let p = curr_pos[g_idx];
        let r = atom_params[g_idx].x;
        local_min[l_idx] = vec4f(p.x - r, p.y - r, p.z - r, 0.0);
        local_max[l_idx] = vec4f(p.x + r, p.y + r, p.z + r, 0.0);
    } else {
        local_min[l_idx] = vec4f(1e10);
        local_max[l_idx] = vec4f(-1e10);
    }
    workgroupBarrier();

    // 2. Parallel Reduction
    for (var stride: u32 = GROUP_SIZE / 2; stride > 0; stride >>= 1) {
        if (l_idx < stride) {
            local_min[l_idx] = min(local_min[l_idx], local_min[l_idx + stride]);
            local_max[l_idx] = max(local_max[l_idx], local_max[l_idx + stride]);
        }
        workgroupBarrier();
    }

    // 3. Write result
    if (l_idx == 0) {
        // each group has two vec4 slots: 2*group -> min, 2*group+1 -> max
        let base = group_id * 2u;
        bboxes[base] = local_min[0];
        bboxes[base + 1u] = local_max[0];
    }
}

// ------------------------------------------------------------------
// KERNEL 2: Build Local Topology
// ------------------------------------------------------------------
@compute @workgroup_size(64, 1, 1)
fn build_local_topology(
    @builtin(local_invocation_id) lid: vec3u,
    @builtin(global_invocation_id) gid: vec3u,
    @builtin(workgroup_id) grp: vec3u
) {
    let l_idx = lid.x;
    let group_id = grp.x;

    if (l_idx == 0) {
        atomicStore(&l_ghost_counter, 0);
        let base = group_id * 2u;
        l_my_bbox_min = bboxes[base];
        l_my_bbox_max = bboxes[base + 1u];
    }
    workgroupBarrier();

    let my_min = l_my_bbox_min;
    let my_max = l_my_bbox_max;
    let margin = u_params.bbox_margin;

    // Step 1: Find Ghosts
    for (var other_g: u32 = 0; other_g < u_params.num_groups; other_g++) {
        if (other_g == group_id) { continue; }

        // BBox Check (packed: 2 vec4 per group)
        let o_base = other_g * 2u;
        let o_min = bboxes[o_base];
        let o_max = bboxes[o_base + 1u];
        
        var overlap = false;
        if (my_max.x + margin >= o_min.x && my_min.x <= o_max.x + margin &&
            my_max.y + margin >= o_min.y && my_min.y <= o_max.y + margin &&
            my_max.z + margin >= o_min.z && my_min.z <= o_max.z + margin) {
            overlap = true;
        }

        if (overlap) {
            let global_target_idx = other_g * GROUP_SIZE + l_idx;
            if (global_target_idx < u_params.num_atoms) {
                let p = curr_pos[global_target_idx];

                // Dist to My BBox
                let dx = max(0.0, max(my_min.x - p.x, p.x - my_max.x));
                let dy = max(0.0, max(my_min.y - p.y, p.y - my_max.y));
                let dz = max(0.0, max(my_min.z - p.z, p.z - my_max.z));
                let dist_sq = dx*dx + dy*dy + dz*dz;

                if (dist_sq < u_params.coll_margin_sq) {
                    let slot = atomicAdd(&l_ghost_counter, 1);
                    if (slot < i32(MAX_GHOSTS)) {
                        l_ghost_list[slot] = i32(global_target_idx);
                    }
                }
            }
        }
    }
    workgroupBarrier();

    // Cap ghost counter
    let total_ghosts_i = min(atomicLoad(&l_ghost_counter), i32(MAX_GHOSTS));
    let total_ghosts: u32 = u32(max(total_ghosts_i, 0));

    // Step 2: Flush ghosts to global memory (packed: [indices... count])
    let base_offset: u32 = group_id * (u32(MAX_GHOSTS) + 1u);
    if (l_idx < total_ghosts) {
        ghost_packed[i32(base_offset + l_idx)] = l_ghost_list[l_idx];
    }
    if (l_idx == 0u) {
        ghost_packed[i32(base_offset + u32(MAX_GHOSTS))] = i32(total_ghosts);
    }

    // Step 3: Map Bonds
    let my_global_id = gid.x;
    if (my_global_id < u_params.num_atoms) {
        // Init with -1
        for (var slot: u32 = 0; slot < N_MAX_BONDED; slot++) {
            bond_indices_local[my_global_id * N_MAX_BONDED + slot] = -1;
        }

        for (var slot: u32 = 0; slot < N_MAX_BONDED; slot++) {
            let t_idx = bond_indices_global[my_global_id * N_MAX_BONDED + slot];
            if (t_idx < 0) { break; }

            let target_grp = u32(t_idx) / GROUP_SIZE;
            
            if (target_grp == group_id) {
                // Internal
                bond_indices_local[my_global_id * N_MAX_BONDED + slot] = t_idx % i32(GROUP_SIZE);
            } else {
                // External - Linear Search in Shared Mem
                for (var k: u32 = 0u; k < total_ghosts; k++) {
                    if (l_ghost_list[k] == t_idx) {
                        // Offset by Internal Size (GROUP_SIZE)
                        bond_indices_local[my_global_id * N_MAX_BONDED + slot] = i32(GROUP_SIZE) + i32(k); 
                        break;
                    }
                }
            }
        }
    }
}

// ------------------------------------------------------------------
// KERNEL 3: Jacobi Solver
// ------------------------------------------------------------------
@compute @workgroup_size(64, 1, 1)
fn solve_cluster_jacobi(
    @builtin(local_invocation_id) lid: vec3u,
    @builtin(global_invocation_id) gid: vec3u,
    @builtin(workgroup_id) grp: vec3u
) {
    let l_idx = lid.x;
    let group_id = grp.x;
    let my_global_id = gid.x;

    // Load Internal
    if (my_global_id < u_params.num_atoms) {
        l_pos[l_idx] = curr_pos[my_global_id];
        l_pos_prev[l_idx] = l_pos[l_idx];
    } else {
        l_pos[l_idx] = vec4f(0.0);
        l_pos_prev[l_idx] = vec4f(0.0);
    }

    let my_pred = select(vec4f(0.0), pred_pos[my_global_id], my_global_id < u_params.num_atoms);
    let params_vec = select(vec4f(0.0), atom_params[my_global_id], my_global_id < u_params.num_atoms);
    let my_rad = params_vec.x;
    let my_mass = params_vec.w;
    
    // alpha = m / dt^2
    let alpha = my_mass / (u_params.dt * u_params.dt);

    if (my_global_id < u_params.num_atoms) {
        l_rad[l_idx] = my_rad;
    } else {
        // IMPORTANT: avoid uninitialized l_rad for inactive lanes (partial groups).
        // Otherwise collision loop may see garbage radii and create huge phantom repulsions.
        l_rad[l_idx] = 0.0;
    }
    workgroupBarrier();

    for (var iter: u32 = 0; iter < u_params.inner_iters; iter++) {
        // Load Ghosts (Cooperative loading)
        let g_offset = group_id * (MAX_GHOSTS + 1u);
        let g_count = u32(ghost_packed[i32(g_offset + MAX_GHOSTS)]);
        
        for (var k = l_idx; k < u32(g_count); k += GROUP_SIZE) {
            let g_idx = ghost_packed[g_offset + k];
            l_pos[GROUP_SIZE + k] = curr_pos[g_idx];
            l_rad[GROUP_SIZE + k] = atom_params[g_idx].x;
        }
        workgroupBarrier();

        if (my_global_id < u_params.num_atoms) {
            let p = l_pos[l_idx].xyz;
            var force = alpha * (my_pred.xyz - p);
            var k_sum = alpha;

            // Bonds
            for (var slot: u32 = 0; slot < N_MAX_BONDED; slot++) {
                let idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
                if (idx == -1) { break; }

                let n_pos = l_pos[idx].xyz;
                let diff = p - n_pos;
                let dist = length(diff);
                if (dist > 1e-6) {
                    let lk = bond_len_stiff[my_global_id * N_MAX_BONDED + slot];
                    let L = lk.x;
                    let st = lk.y;
                    let C = dist - L;
                    let f_vec = -st * C * (diff / dist);
                    force += f_vec;
                    k_sum += st;
                }
            }

            // Collisions
            let check_count = GROUP_SIZE + u32(g_count);
            for (var j: u32 = 0; j < check_count; j++) {
                if (j == l_idx) { continue; }

                // Skip Bonded
                var is_bonded = false;
                for (var slot: u32 = 0; slot < N_MAX_BONDED; slot++) {
                    let idx = bond_indices_local[my_global_id * N_MAX_BONDED + slot];
                    if (idx == -1) { break; }
                    if (idx == i32(j)) { is_bonded = true; break; }
                }
                if (is_bonded) { continue; }

                let other = l_pos[j];
                let diff = p - other.xyz;
                let dist_sq = dot(diff, diff);
                let r_sum = my_rad + l_rad[j];

                if (dist_sq < r_sum * r_sum && dist_sq > 1e-12) {
                    let dist = sqrt(dist_sq);
                    let overlap = dist - r_sum;
                    let f_vec = -u_params.k_coll * overlap * (diff / dist);
                    force += f_vec;
                    k_sum += u_params.k_coll;
                }
            }

            var p_new = p + (force / k_sum) * u_params.omega;
            if (u_params.momentum_beta != 0.0) {
                p_new += (p - l_pos_prev[l_idx].xyz) * u_params.momentum_beta;
            }
            l_pos_new[l_idx] = vec4f(p_new, 0.0);
        }
        workgroupBarrier();

        if (my_global_id < u_params.num_atoms) {
            l_pos_prev[l_idx] = l_pos[l_idx];
            l_pos[l_idx] = l_pos_new[l_idx];
        }
        workgroupBarrier();
    }

    if (my_global_id < u_params.num_atoms) {
        curr_pos[my_global_id] = l_pos[l_idx];
    }
}