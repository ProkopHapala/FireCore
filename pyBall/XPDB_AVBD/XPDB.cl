// XPDB.cl

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define GROUP_SIZE 64      // Must match local_work_size in Python
#define MAX_NEIGHBORS 64   // Hard limit for collision list storage

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
// KERNEL 2: Update Verlet List (Broadphase + Narrowphase)
// ------------------------------------------------------------------
__kernel void update_verlet_list(
    __global const float4* curr_pos,
    __global const float4* params,       // .x = radius
    __global const float4* bboxes_min,
    __global const float4* bboxes_max,
    __global int* coll_neighbors,        // Flat array: [atom0_n0, atom0_n1..., atom1_n0...]
    __global int* coll_count,            // Count per atom
    int num_atoms,
    int num_groups,
    float margin_sq,    // (2*Rmax + margin)^2
    float bbox_margin   // margin for BBox check
) {
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float4 p_i = curr_pos[i];
    float r_i = params[i].x;
    int my_group = i / GROUP_SIZE; // Implicit grouping

    // Cache my group's bbox to avoid re-reading for optimization? 
    // Actually, we check OTHER groups, so we read bboxes_min[g]
    float4 my_box_min = bboxes_min[my_group];
    float4 my_box_max = bboxes_max[my_group];

    int count = 0;
    int base_offset = i * MAX_NEIGHBORS;

    // Iterate over ALL groups to find candidates
    // (This is N_atoms * N_groups complexity - much better than N^2)
    for (int g = 0; g < num_groups; g++) {
        
        // Skip self-group? 
        // No, we need to check self-collisions within the group too, 
        // unless we rely on bonds. But let's check everything except self-index.
        
        // 1. Broadphase: Box vs Box
        if (!bboxes_overlap(my_box_min, my_box_max, bboxes_min[g], bboxes_max[g], bbox_margin)) {
            continue;
        }

        // 2. Narrowphase: Atom vs Atoms in Group g
        int start_j = g * GROUP_SIZE;
        int end_j = min(start_j + GROUP_SIZE, num_atoms);

        for (int j = start_j; j < end_j; j++) {
            if (i == j) continue; // Skip self

            float4 p_j = curr_pos[j];
            float r_j = params[j].x;

            float dx = p_i.x - p_j.x;
            float dy = p_i.y - p_j.y;
            float dz = p_i.z - p_j.z;
            float dist_sq = dx*dx + dy*dy + dz*dz;

            // Check against Verlet Margin (R_sum + margin)^2
            // Note: margin_sq passed in should be roughly (R_max*2 + physical_margin)^2
            // Or we compute precise threshold:
            float threshold = (r_i + r_j) + 1.0f; // Arbitrary margin e.g. 1.0
            // Using the user provided squared margin for simplicity:
            
            if (dist_sq < margin_sq) {
                if (count < MAX_NEIGHBORS) {
                    coll_neighbors[base_offset + count] = j;
                    count++;
                }
            }
        }
    }
    coll_count[i] = count;
}

// ------------------------------------------------------------------
// KERNEL 3: Predict (Momentum)
// ------------------------------------------------------------------
__kernel void predict(
    __global const float4* curr_pos,
    __global float4* velocity,
    __global const float4* params,   // .x = radius, .w = mass
    __global float4* pred_pos,
    int num_atoms,
    float dt,
    float gx, float gy, float gz,
    float bxmin, float bymin, float bzmin,
    float bxmax, float bymax, float bzmax,
    float box_k,
    float bottom_y,
    float bottom_c,
    float bottom_x0,
    float bottom_a
){
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float4 p = curr_pos[i];
    float4 v4 = velocity[i];
    float3 v = v4.xyz;

    float mass = params[i].w;
    float3 force = mass * (float3)(gx, gy, gz);

    float r = params[i].x; // radius
    float xmin = bxmin + r;
    float xmax = bxmax - r;
    float ymin = bymin + r;
    float ymax = bymax - r;
    float zmin = bzmin + r;
    float zmax = bzmax - r;

    // Parabolic bottom force: y_floor = bottom_y + bottom_a*(x-bottom_x0)^2
    if (bottom_a != 0.0f) {
        float y_floor = bottom_y + bottom_a * (p.x - bottom_x0) * (p.x - bottom_x0);
        if (p.y < y_floor) {
            float d = y_floor - p.y;
            float dy_dx = 2.0f * bottom_a * (p.x - bottom_x0);
            float floor_k = (bottom_c > 0.0f) ? bottom_c : box_k;
            force.y += floor_k * d;
            force.x -= floor_k * d * dy_dx;
        }
    }

    // Flat wall springs (sides/z only; no flat floor/top when parabola active)
    if (p.x < xmin) force.x += box_k * (xmin - p.x);
    else if (p.x > xmax) force.x -= box_k * (p.x - xmax);
    if (bottom_c <= 0.0f && bottom_a == 0.0f) {
        if (p.y < ymin) force.y += box_k * (ymin - p.y);
        else if (p.y > ymax) force.y -= box_k * (p.y - ymax);
    } // otherwise parabola handles floor; no top wall
    if (p.z < zmin) force.z += box_k * (zmin - p.z);
    else if (p.z > zmax) force.z -= box_k * (p.z - zmax);

    float inv_m = (mass > 0.0f) ? (1.0f / mass) : 0.0f;
    v          += force * inv_m * dt;
    float3 y    = p.xyz + v * dt;

    velocity[i] = (float4)(v, v4.w);
    pred_pos[i] = (float4)(y, p.w);
}

__kernel void predict_position_based(
    __global const float4* curr_pos,
    __global const float4* velocity,
    __global const float4* params,   // .x = radius
    __global float4* pred_pos,
    int num_atoms,
    float dt,
    float gx, float gy, float gz,
    float bxmin, float bymin, float bzmin,
    float bxmax, float bymax, float bzmax,
    float box_k,
    float bottom_y,
    float bottom_c,
    float bottom_x0,
    float bottom_a
){
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float4 p = curr_pos[i];
    float4 v = velocity[i];

    float3 g = (float3)(gx, gy, gz);
    float3 y = p.xyz + v.xyz * dt + g * (dt * dt); // inertial + gravity

    // Simple harmonic walls: if outside, pull back inside with spring displacement
    float r = params[i].x; // radius
    float xmin = bxmin + r;
    float xmax = bxmax - r;
    float ymin = bymin + r;
    float ymax = bymax - r;
    float zmin = bzmin + r;
    float zmax = bzmax - r;

    // Parabolic bottom push (optional) floor: y_floor = bottom_y + bottom_a*(x-bottom_x0)^2
    if (bottom_c > 0.0f) {
        float y_floor = bottom_y + bottom_a * (y.x - bottom_x0) * (y.x - bottom_x0);
        if (y.y < y_floor) {
            float d = y_floor - y.y;
            float dy_dx = 2.0f * bottom_a * (y.x - bottom_x0);
            // Vertical correction
            y.y += bottom_c * d * d * dt * dt;
            // Lateral correction along gradient of parabola
            y.x -= bottom_c * d * dy_dx * dt * dt;
        }
    }

    // Flat walls
    if (y.x < xmin) y.x += box_k * (xmin - y.x) * dt * dt;
    else if (y.x > xmax) y.x -= box_k * (y.x - xmax) * dt * dt;
    if (y.y < ymin) y.y += box_k * (ymin - y.y) * dt * dt;
    else if (y.y > ymax) y.y -= box_k * (y.y - ymax) * dt * dt;
    if (y.z < zmin) y.z += box_k * (zmin - y.z) * dt * dt;
    else if (y.z > zmax) y.z -= box_k * (y.z - zmax) * dt * dt;

    pred_pos[i] = (float4)(y, p.w);
}

// ------------------------------------------------------------------
// KERNEL 4: Unified Solver (Pure Gather, Stateless)
// ------------------------------------------------------------------
__kernel void solve_simple_jacobi(
    __global const float4* iter_pos_in,   // Read: Position from Prev Iter
    __global const float4* pred_pos,      // Read: Inertial Target
    __global const float4* params,        // Read: .w = mass, .x = radius
    __global float4* iter_pos_out,        // Write: Next Position
    
    // Bonds
    __global const int* bond_start,
    __global const int* bond_count,
    __global const int* bond_neighbors,
    __global const float* bond_lengths,
    __global const float* bond_stiffness,
    
    // Collisions
    __global const int* coll_neighbors,
    __global const int* coll_count,
    
    int num_atoms,
    float dt,
    float k_collision,
    float omega, // Relaxation factor
    float pd_scale // multiply mass/dt^2 diagonal (projective inertia term)
) {
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float4 p_i = iter_pos_in[i];
    float4 tgt = pred_pos[i];
    float mass = params[i].w;
    float r_i = params[i].x;

    // 1. Inertial Potential (Momentum) with projective diagonal
    // F = alpha * (target - x), alpha = pd_scale * mass/dt^2
    float alpha = pd_scale * mass / (dt * dt);
    float3 total_force = alpha * (tgt.xyz - p_i.xyz);
    float total_stiffness = alpha + 1e-8f; // avoid divide by zero

    // 2. Bonds
    int b_start = bond_start[i];
    int b_cnt = bond_count[i];
    
    for (int k = 0; k < b_cnt; k++) {
        int idx = b_start + k;
        int j = bond_neighbors[idx];
        float L = bond_lengths[idx];
        float k_bond = bond_stiffness[idx];
        
        float4 p_j = iter_pos_in[j];
        float3 diff = p_i.xyz - p_j.xyz;
        float dist = length(diff);
        
        if (dist > 1e-6f) {
            float3 n = diff / dist;
            float C = dist - L;
            
            // F = -k * C
            total_force -= k_bond * C * n;
            total_stiffness += k_bond;
        }
    }

    // 3. Collisions
    int c_cnt = coll_count[i];
    int c_base = i * MAX_NEIGHBORS;
    
    for (int k = 0; k < c_cnt; k++) {
        int j = coll_neighbors[c_base + k];
        float4 p_j = iter_pos_in[j];
        float r_j = params[j].x;
        
        float3 diff = p_i.xyz - p_j.xyz;
        float dist_sq = dot(diff, diff);
        float r_sum = r_i + r_j;
        
        if (dist_sq < r_sum * r_sum) {
            float dist = sqrt(dist_sq);
            if (dist > 1e-6f) {
                float3 n = diff / dist;
                float overlap = dist - r_sum; // Negative
                
                // Repulsion
                total_force -= k_collision * overlap * n;
                total_stiffness += k_collision;
            }
        }
    }

    // 4. Update
    float3 dx = total_force / total_stiffness;
    float3 final_pos = p_i.xyz + dx * omega;
    
    iter_pos_out[i] = (float4)(final_pos, p_i.w);
}

// ------------------------------------------------------------------
// KERNEL 5: Integrate
// ------------------------------------------------------------------
__kernel void integrate(
    __global const float4* start_pos, // Position before solver loop
    __global const float4* final_pos, // Result of solver loop
    __global float4* velocity,
    __global float4* out_pos_main,    // Write back to main buffer
    int num_atoms,
    float dt
) {
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float3 p_s = start_pos[i].xyz;
    float3 p_f = final_pos[i].xyz;
    
    // v = (x_new - x_old) / dt
    float3 v = (p_f - p_s) / dt;
    
    // Damping (Drag)
    v *= 0.995f; 
    
    // Update State
    velocity[i] = (float4)(v, 0.0f);
    out_pos_main[i] = final_pos[i];
}

// ------------------------------------------------------------------
// KERNEL 6: Diagnostics (per-atom min/max violation)
// ------------------------------------------------------------------
__kernel void diag_constraints(
    __global const float4* pos,
    __global const float4* params,        // .x radius
    // Bonds
    __global const int* bond_start,
    __global const int* bond_count,
    __global const int* bond_neighbors,
    __global const float* bond_lengths,
    // Collisions
    __global const int* coll_neighbors,
    __global const int* coll_count,
    // Output
    __global float2* out_bond, // (min,max) violation per atom
    __global float2* out_coll,
    int num_atoms
) {
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float4 p_i = pos[i];
    float rad_i = params[i].x;

    float bond_min = 1e20f;
    float bond_max = -1e20f;
    float coll_min = 1e20f;
    float coll_max = -1e20f;

    // Bonds: violation = dist - L
    int b_start = bond_start[i];
    int b_cnt = bond_count[i];
    for (int k = 0; k < b_cnt; k++) {
        int j = bond_neighbors[b_start + k];
        float L = bond_lengths[b_start + k];
        float4 p_j = pos[j];
        float3 diff = p_i.xyz - p_j.xyz;
        float dist = length(diff);
        float v = dist - L;
        bond_min = fmin(bond_min, v);
        bond_max = fmax(bond_max, v);
    }

    // Collisions: penetration = dist - (r_i + r_j)
    int c_cnt = coll_count[i];
    int c_base = i * MAX_NEIGHBORS;
    for (int k = 0; k < c_cnt; k++) {
        int j = coll_neighbors[c_base + k];
        float4 p_j = pos[j];
        float rad_j = params[j].x;
        float3 diff = p_i.xyz - p_j.xyz;
        float dist = length(diff);
        float v = dist - (rad_i + rad_j);
        coll_min = fmin(coll_min, v);
        coll_max = fmax(coll_max, v);
    }

    if (b_cnt == 0) { bond_min = 0.0f; bond_max = 0.0f; }
    if (c_cnt == 0) { coll_min = 0.0f; coll_max = 0.0f; }

    out_bond[i] = (float2)(bond_min, bond_max);
    out_coll[i] = (float2)(coll_min, coll_max);
}

// ------------------------------------------------------------------
// KERNEL 7: Clamp Z to zero (for 2D debug)
// ------------------------------------------------------------------
__kernel void clamp_z0(
    __global float4* pos,
    __global float4* vel,
    int num_atoms
){
    int i = get_global_id(0);
    if (i >= num_atoms) return;
    float4 p = pos[i];
    float4 v = vel[i];
    p.z = 0.0f;
    v.z = 0.0f;
    pos[i] = p;
    vel[i] = v;
}