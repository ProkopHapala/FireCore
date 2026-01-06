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
    __global const float4* prev_pos_in,   // Read: Position from Iter-1 (for momentum)
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
    float omega,   // Relaxation factor
    float pd_scale, // multiply mass/dt^2 diagonal (projective inertia term)
    float momentum_beta // heavy-ball momentum coefficient
) {
    int i = get_global_id(0);
    if (i >= num_atoms) return;

    float4 p_i = iter_pos_in[i];
    float4 p_prev = prev_pos_in[i];
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
    // Heavy-ball: add momentum based on last displacement (p_i - p_prev)
    if (momentum_beta != 0.0f) {
        final_pos += (p_i.xyz - p_prev.xyz) * momentum_beta;
    }
    
    iter_pos_out[i] = (float4)(final_pos, p_i.w);
}

// ------------------------------------------------------------------
// KERNEL 4b: Chebyshev Mix (out = omega*(tmp - prev2) + prev2)
// ------------------------------------------------------------------
__kernel void cheby_mix(
    __global const float4* prev2,
    __global const float4* tmp,
    __global float4* out,
    float omega
) {
    int i = get_global_id(0);
    float4 p2 = prev2[i];
    float4 t = tmp[i];
    out[i] = p2 + (t - p2) * omega;
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









// GS_Solver.cl

#define WORKGROUP_SIZE 64

// Helper to get position from Local or Global based on index
inline float4 get_neighbor_pos(
    int neighbor_idx, 
    int group_start_idx, 
    __local float4* local_pos, 
    __global const float4* global_pos
) {
    // Check if neighbor is inside this workgroup
    if (neighbor_idx >= group_start_idx && neighbor_idx < group_start_idx + WORKGROUP_SIZE) {
        // Fast path: Read from Local Memory
        return local_pos[neighbor_idx - group_start_idx];
    } else {
        // Slow path: Read from Global Memory (Halo/Ghost region)
        // Note: This data remains static during the local iterations (boundary condition)
        return global_pos[neighbor_idx];
    }
}

__kernel void solve_gauss_seidel_local(
    // Global Data
    __global float4* curr_pos,          // In/Out: current positions
    __global const float4* pred_pos,    // Inertial target positions (y)
    __global const float4* params,      // .w = mass, .x = radius
    
    // Topology (Bonds)
    __global const int* bond_start,     // Start index in neighbor list
    __global const int* bond_count,     // Number of neighbors
    __global const int* bond_neighbors, // The neighbor indices
    __global const float* bond_lengths,
    __global const float* bond_stiffness,

    // Coloring Data
    __global const int* atom_colors,    // Color of each atom (0..num_colors-1)
    int num_colors,
    int inner_iterations,               // How many GS loops to run inside kernel
    float dt,
    float omega                         // Relaxation factor (SOR)
) {
    // -----------------------------------------------------------
    // 1. SETUP & LOAD TO LOCAL MEMORY
    // -----------------------------------------------------------
    int gid = get_global_id(0);
    int lid = get_local_id(0);
    int group_id = get_group_id(0);
    int group_start = group_id * WORKGROUP_SIZE;

    // Local Memory Buffer
    __local float4 l_pos[WORKGROUP_SIZE];

    // Cooperative Load
    float4 my_pos = curr_pos[gid];
    l_pos[lid] = my_pos;
    
    // Load static data to registers
    float4 my_pred = pred_pos[gid];
    float my_mass = params[gid].w;
    float inertia_alpha = my_mass / (dt * dt);
    
    // Topology (Read once)
    int my_bond_start = bond_start[gid];
    int my_bond_count = bond_count[gid];
    int my_color = atom_colors[gid];

    // Sync to ensure everyone is loaded
    barrier(CLK_LOCAL_MEM_FENCE);

    // -----------------------------------------------------------
    // 2. GAUSS-SEIDEL LOOP (Inside Kernel)
    // -----------------------------------------------------------
    for (int iter = 0; iter < inner_iterations; iter++) {
        
        // Loop through colors sequentially
        for (int c = 0; c < num_colors; c++) {
            
            // Only threads with the active color work
            if (my_color == c) {
                
                float3 total_force = (float3)(0.0f);
                float total_stiffness = 0.0f;

                // A. Inertial Constraint (Target)
                // F = alpha * (y - x)
                // Note: We use l_pos[lid] (current state)
                float3 curr_xyz = l_pos[lid].xyz;
                total_force += inertia_alpha * (my_pred.xyz - curr_xyz);
                total_stiffness += inertia_alpha;

                // B. Bond Constraints
                for (int k = 0; k < my_bond_count; k++) {
                    int idx = my_bond_start + k;
                    int neighbor_gid = bond_neighbors[idx];
                    float L = bond_lengths[idx];
                    float k_bond = bond_stiffness[idx];

                    // Access Neighbor (Local or Global)
                    float4 n_pos = get_neighbor_pos(neighbor_gid, group_start, l_pos, curr_pos);
                    
                    float3 diff = curr_xyz - n_pos.xyz;
                    float dist = length(diff);

                    if (dist > 1e-8f) {
                        float3 n = diff / dist;
                        float C = dist - L;
                        
                        // Force
                        total_force -= k_bond * C * n;
                        total_stiffness += k_bond;
                    }
                }

                // C. Solve & Update Local Memory
                // dx = Force / Stiffness
                float3 dx = total_force / total_stiffness;
                
                // Write back to LOCAL memory immediately
                // This makes the update visible to the next color in the next Barrier
                l_pos[lid].xyz += dx * omega;
            }
            
            // CRITICAL BARRIER
            // Wait for Color 'c' to finish writing to l_pos before Color 'c+1' starts reading
            barrier(CLK_LOCAL_MEM_FENCE); 
        }
    }

    // -----------------------------------------------------------
    // 3. WRITE BACK TO GLOBAL
    // -----------------------------------------------------------
    curr_pos[gid] = l_pos[lid];
}


// BlockGS.cl

// ------------------------------------------------------------------
// CONFIGURATION
// ------------------------------------------------------------------
#define WG_SIZE 64        // M: Threads per WorkGroup
#define NUM_COLORS 8      // K: Number of Colors/Partitions
#define BLOCK_CAPACITY (WG_SIZE * NUM_COLORS) // Total atoms in Local Mem (512)

// Helper: Smart Fetch from Cache or Global
inline float4 get_pos_cached(
    int neighbor_gid,
    int block_gid_start,
    __local float4* l_pos,
    __global const float4* g_pos
) {
    // Is neighbor inside our cached block?
    // Optimization: (unsigned) cast trick handles negative check implicitly if needed
    int local_offset = neighbor_gid - block_gid_start;
    
    if (local_offset >= 0 && local_offset < BLOCK_CAPACITY) {
        return l_pos[local_offset]; // HIT: Fast Local Memory
    } else {
        return g_pos[neighbor_gid]; // MISS: Slow Global Memory (Halo)
    }
}

__kernel void solve_block_gs(
    __global float4* curr_pos,        // RW: Current Positions
    __global const float4* pred_pos,  // R: Inertial Target
    __global const float4* params,    // R: .w = mass
    
    // Topology
    __global const int* bond_start,
    __global const int* bond_count,
    __global const int* bond_neighbors,
    __global const float* bond_lengths,
    __global const float* bond_stiffness,
    
    int num_atoms,
    int inner_iters,
    float dt,
    float omega
) {
    // -------------------------------------------------------
    // 1. COOPERATIVE LOAD (Global -> Local)
    // -------------------------------------------------------
    // The WorkGroup loads a chunk of (WG_SIZE * NUM_COLORS) atoms
    __local float4 l_pos[BLOCK_CAPACITY];
    
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int block_start_gid = grp * BLOCK_CAPACITY; // Start of this block in global mem

    // Each thread loads K atoms (one per color slot)
    // Coalesced access: Thread 0 reads 0, 64, 128...
    for (int k = 0; k < NUM_COLORS; k++) {
        int local_idx = lid + k * WG_SIZE;
        int global_idx = block_start_gid + local_idx;
        
        if (global_idx < num_atoms) {
            l_pos[local_idx] = curr_pos[global_idx];
        } else {
            l_pos[local_idx] = (float4)(0.0f); // Padding
        }
    }
    
    // Also preload mass/pred if you have enough local memory (Recommended!)
    // For now, we read those from global inside loop to save register pressure, 
    // assuming they are cached L2.

    barrier(CLK_LOCAL_MEM_FENCE);

    // -------------------------------------------------------
    // 2. INNER RELAXATION LOOP
    // -------------------------------------------------------
    for (int iter = 0; iter < inner_iters; iter++) {
        
        // Iterate over colors sequentially
        for (int c = 0; c < NUM_COLORS; c++) {
            
            // Thread 'lid' is responsible for the atom at color 'c'
            int my_local_idx = lid + c * WG_SIZE;
            int my_global_idx = block_start_gid + my_local_idx;
            
            if (my_global_idx < num_atoms) {
                
                // Read my own data from LOCAL
                float4 p_i = l_pos[my_local_idx];
                
                // Read constant data (L2 cache hit usually)
                float4 tgt = pred_pos[my_global_idx];
                float mass = params[my_global_idx].w;
                float inertia = mass / (dt*dt);
                
                float3 total_force = inertia * (tgt.xyz - p_i.xyz);
                float total_k = inertia;
                
                // Process Bonds
                int start = bond_start[my_global_idx];
                int count = bond_count[my_global_idx];
                
                for (int n = 0; n < count; n++) {
                    int neighbor_gid = bond_neighbors[start + n];
                    float L = bond_lengths[start + n];
                    float k_bond = bond_stiffness[start + n];
                    
                    // SMART FETCH: Local or Global?
                    float4 p_j = get_pos_cached(neighbor_gid, block_start_gid, l_pos, curr_pos);
                    
                    float3 diff = p_i.xyz - p_j.xyz;
                    float dist = length(diff);
                    
                    if (dist > 1e-8f) {
                        float3 dir = diff / dist;
                        float C = dist - L;
                        total_force -= k_bond * C * dir;
                        total_k += k_bond;
                    }
                }
                
                // Update Position
                float3 dx = total_force / total_k;
                p_i.xyz += dx * omega;
                
                // WRITE BACK TO LOCAL IMMEDIATELY
                // So the next color sees the updated position
                l_pos[my_local_idx] = p_i;
            }
            
            // Wait for color 'c' to finish before color 'c+1' starts
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

    // -------------------------------------------------------
    // 3. COOPERATIVE STORE (Local -> Global)
    // -------------------------------------------------------
    for (int k = 0; k < NUM_COLORS; k++) {
        int local_idx = lid + k * WG_SIZE;
        int global_idx = block_start_gid + local_idx;
        
        if (global_idx < num_atoms) {
            curr_pos[global_idx] = l_pos[local_idx];
        }
    }
}


#define WG_SIZE 64
#define MAX_GHOSTS 128  // Max external atoms per cluster
#define TOTAL_CAP (WG_SIZE + MAX_GHOSTS)

// ------------------------------------------------------------------
// KERNEL 1: GHOST BUILDER & RE-INDEXER
// ------------------------------------------------------------------
// Work distribution: 1 WorkGroup = 1 Cluster
__kernel void build_local_topology(
    // Input
    __global const float4* curr_pos,
    __global const float4* bboxes_min,
    __global const float4* bboxes_max,
    __global const int4* global_bond_indices, // -1 for no bond
    
    // Output
    __global int* ghost_indices_flat,     // [Cluster * MAX_GHOSTS + k]
    __global int* ghost_counts,           // [Cluster]
    __global int4* local_bond_indices,    // The translated map (0..TOTAL_CAP)
    
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
    // To parallelize, threads stride through the group list
    for (int other_g = 0; other_g < num_groups; other_g++) {
        if (other_g == grp) continue;

        // 1. Check BBox Overlap (Only Thread 0 needs to check? 
        // Or we check cooperatively? Let's have Thread 0 check to save bandwidth)
        // Optimization: This serialization of BBox check is fine because num_groups is small.
        
        bool overlap = false;
        // Simple AABB check
        float4 o_min = bboxes_min[other_g];
        float4 o_max = bboxes_max[other_g];
        
        if (my_max.x + bbox_margin >= o_min.x && my_min.x <= o_max.x + bbox_margin &&
            my_max.y + bbox_margin >= o_min.y && my_min.y <= o_max.y + bbox_margin &&
            my_max.z + bbox_margin >= o_min.z && my_min.z <= o_max.z + bbox_margin) {
            overlap = true;
        }

        // 2. If overlap, all threads load atoms from 'other_g' and check distance
        if (overlap) {
            // "Vote" to check this group
            // Every thread loads 1 atom from the other group
            // Assuming WG_SIZE == atoms per group
            int global_idx = other_g * WG_SIZE + lid;
            
            if (global_idx < num_atoms) {
                float4 p = curr_pos[global_idx];
                
                // Distance to My BBox (Clamped Distance)
                float dx = max(0.0f, max(my_min.x - p.x, p.x - my_max.x));
                float dy = max(0.0f, max(my_min.y - p.y, p.y - my_max.y));
                float dz = max(0.0f, max(my_min.z - p.z, p.z - my_max.z));
                float dist_sq = dx*dx + dy*dy + dz*dz;
                
                // If close enough to potentially collide
                if (dist_sq < margin_sq) {
                    // Add to Ghost List using Atomic
                    int slot = atomic_inc(&l_ghost_counter);
                    if (slot < MAX_GHOSTS) {
                        l_ghost_list[slot] = global_idx;
                    }
                }
            }
        }
        // No barrier needed here strictly, as atomic_inc handles safety
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);
    
    // Cap counter
    int total_ghosts = min(l_ghost_counter, MAX_GHOSTS);

    // ------------------------------------------------
    // STEP 2: Flush Ghosts to Global
    // ------------------------------------------------
    // Threads cooperatively copy local list to global
    // (So the solver can read them later)
    int base_offset = grp * MAX_GHOSTS;
    for (int i = lid; i < total_ghosts; i += WG_SIZE) {
        ghost_indices_flat[base_offset + i] = l_ghost_list[i];
    }
    if (lid == 0) ghost_counts[grp] = total_ghosts;

    // ------------------------------------------------
    // STEP 3: Re-Index Bonds (The "Map")
    // ------------------------------------------------
    // Each thread handles its own atom's bonds
    int my_global_id = grp * WG_SIZE + lid;
    
    if (my_global_id < num_atoms) {
        int4 g_bonds = global_bond_indices[my_global_id];
        int4 l_bonds = (int4)(-1);

        // Helper macro to process one bond slot
        #define TRANSLATE(comp) \
        if (g_bonds.comp != -1) { \
            int target = g_bonds.comp; \
            int target_grp = target / WG_SIZE; \
            if (target_grp == grp) { \
                /* Internal: Simple math */ \
                l_bonds.comp = target % WG_SIZE; \
            } else { \
                /* External: Search Ghost List */ \
                /* Linear search is OK because list is small and in Shared Mem */ \
                for (int k = 0; k < total_ghosts; k++) { \
                    if (l_ghost_list[k] == target) { \
                        l_bonds.comp = WG_SIZE + k; /* Offset by Internal Size */ \
                        break; \
                    } \
                } \
            } \
        }

        TRANSLATE(x);
        TRANSLATE(y);
        TRANSLATE(z);
        TRANSLATE(w);
        
        // Write map to global
        local_bond_indices[my_global_id] = l_bonds;
    }
}


// ------------------------------------------------------------------
// KERNEL 2: THE SOLVER (Purely Local)
// ------------------------------------------------------------------
__kernel void solve_cluster_local(
    __global float4* curr_pos,          // RW
    __global const float4* pred_pos,    // R
    __global const float4* params,      // R
    
    // Topology (Pre-cached)
    __global const int4* local_bond_indices, // 0..TOTAL_CAP
    __global const float4* bond_lengths,
    __global const float4* bond_stiffness,
    
    // Ghost Data
    __global const int* ghost_indices_flat,
    __global const int* ghost_counts,

    // Coloring
    __global const int* atom_colors,
    
    int num_atoms,
    int inner_iters,
    float dt,
    float k_coll,
    float omega
) {
    int lid = get_local_id(0);
    int grp = get_group_id(0);
    int my_global_id = grp * WG_SIZE + lid;

    // -------------------------------------------
    // 1. UNIFIED LOCAL MEMORY LOAD
    // -------------------------------------------
    // Indices [0..63] are Internal
    // Indices [64..TOTAL_CAP] are Ghosts
    __local float4 l_pos[TOTAL_CAP];

    // A. Load Internal
    if (my_global_id < num_atoms) {
        l_pos[lid] = curr_pos[my_global_id];
    } else {
        l_pos[lid] = (float4)(0.0f); // Padding
    }

    // B. Load Ghosts
    int g_count = ghost_counts[grp];
    int g_offset = grp * MAX_GHOSTS;
    
    // Cooperative load
    for (int k = lid; k < g_count; k += WG_SIZE) {
        int g_idx = ghost_indices_flat[g_offset + k];
        l_pos[WG_SIZE + k] = curr_pos[g_idx];
    }

    // Load Constants to Registers
    float4 my_pred = (my_global_id < num_atoms) ? pred_pos[my_global_id] : (float4)(0.0f);
    float my_mass = (my_global_id < num_atoms) ? params[my_global_id].w : 1.0f;
    float my_rad = (my_global_id < num_atoms) ? params[my_global_id].x : 0.0f;
    
    int4 my_bonds = (my_global_id < num_atoms) ? local_bond_indices[my_global_id] : (int4)(-1);
    float4 my_L = (my_global_id < num_atoms) ? bond_lengths[my_global_id] : (float4)(0.0f);
    float4 my_K = (my_global_id < num_atoms) ? bond_stiffness[my_global_id] : (float4)(0.0f);
    
    int my_color = (my_global_id < num_atoms) ? atom_colors[my_global_id] : -1;
    float alpha = my_mass / (dt*dt);

    barrier(CLK_LOCAL_MEM_FENCE);

    // -------------------------------------------
    // 2. SOLVER LOOP (No Global Reads!)
    // -------------------------------------------
    // We iterate over colors to ensure stability
    
    for (int iter = 0; iter < inner_iters; iter++) {
        
        // Assume 2 Colors (Red/Black) usually sufficient
        for (int c = 0; c < 2; c++) {
            
            if (my_color == c && my_global_id < num_atoms) {
                
                float3 p = l_pos[lid].xyz;
                float3 force = alpha * (my_pred.xyz - p);
                float k_sum = alpha;

                // --- A. Bonds (Using Local Indices) ---
                #define APPLY_BOND(comp) \
                if (my_bonds.comp != -1) { \
                    int idx = my_bonds.comp; \
                    float3 n_pos = l_pos[idx].xyz; \
                    float3 diff = p - n_pos; \
                    float dist = length(diff); \
                    if (dist > 1e-6f) { \
                        float C = dist - my_L.comp; \
                        float st = my_K.comp; \
                        force -= st * C * (diff/dist); \
                        k_sum += st; \
                    } \
                }

                APPLY_BOND(x);
                APPLY_BOND(y);
                APPLY_BOND(z);
                APPLY_BOND(w);

                // --- B. Collisions (All-to-All in Local Mem) ---
                // Check against everything: Internal + Ghosts
                int check_count = WG_SIZE + g_count;
                
                for (int j = 0; j < check_count; j++) {
                    if (j == lid) continue; // Skip self
                    
                    // Skip if bonded? 
                    // Optimization: We can check against my_bonds registers
                    if (j == my_bonds.x || j == my_bonds.y || j == my_bonds.z || j == my_bonds.w) continue;

                    float4 other = l_pos[j];
                    // Optimization: If ghost is padded/invalid, other.w might be 0, assume mass/rad > 0
                    // Let's assume params are stored in .x (radius) .w (mass) of pos vector for simplicity here
                    // Or we assume radius is uniform? 
                    // Let's assume other.x is position, and we lack radius info for ghost?
                    // **CRITICAL**: Ghost load needs to load radius too. 
                    // We assume curr_pos.w contains mass? No, typically radius is in params.
                    // For ghosts, we loaded curr_pos. We need ghost params.
                    // *Fix*: We can pack Radius into curr_pos.w for the Solver phase.
                    
                    float3 diff = p - other.xyz;
                    float dist_sq = dot(diff, diff);
                    
                    // Assume radius stored in .w of pos for this kernel context
                    float r_sum = my_rad + 0.5f; // Hack: need to load ghost radius properly
                    
                    if (dist_sq < r_sum * r_sum) {
                        float dist = sqrt(dist_sq);
                        float overlap = dist - r_sum;
                        force -= k_coll * overlap * (diff/dist);
                        k_sum += k_coll;
                    }
                }

                // Update
                p += (force / k_sum) * omega;
                l_pos[lid] = (float4)(p, 0.0f); // Preserve w?
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
    }

    // -------------------------------------------
    // 3. WRITE BACK
    // -------------------------------------------
    if (my_global_id < num_atoms) {
        curr_pos[my_global_id] = l_pos[lid];
    }
}