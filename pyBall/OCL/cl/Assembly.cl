typedef struct {
    float4 r0;
    float4 r1;
    float4 r2;
    float4 shift;
} Transform;

// Helper function to apply the rigid body transformation to a 3D atom
inline float4 apply_transform(float4 atom, Transform t) {
    float4 res;
    // 3x3 Matrix multiplication + Shift
    res.x = t.r0.x * atom.x + t.r0.y * atom.y + t.r0.z * atom.z + t.shift.x;
    res.y = t.r1.x * atom.x + t.r1.y * atom.y + t.r1.z * atom.z + t.shift.y;
    res.z = t.r2.x * atom.x + t.r2.y * atom.y + t.r2.z * atom.z + t.shift.z;
    res.w = atom.w; // Preserve the Van der Waals radius
    return res;
}


// 1 Thread = 1 Atom in the final exported system
__kernel void emit_configuration_xyz(
    __global const float4* base_atoms,    // The original molecule (natoms)
    int natoms,                           // Number of atoms in one molecule
    __global const Transform* transforms, // Array of 54 transforms for THIS specific config
    int nmols,                            // Total molecules (e.g., 54)
    __global float4* out_atoms            // Output buffer, size: (natoms * nmols)
) {
    // Get the global 1D thread ID
    int global_id = get_global_id(0);
    
    // Ensure we don't read/write out of bounds
    if (global_id >= (natoms * nmols)) return;

    // Determine WHICH molecule and WHICH atom this thread is responsible for
    int mol_idx  = global_id / natoms;
    int atom_idx = global_id % natoms;

    // Read the base atom and the corresponding transformation
    float4 base_atom = base_atoms[atom_idx];
    Transform t = transforms[mol_idx];

    // Transform and write straight to global memory
    out_atoms[global_id] = apply_transform(base_atom, t);
}

__kernel void evaluate_packing_3d(
    __constant float4* base_atoms,        // - Base molecule coordinates
    int natoms,                           // Number of atoms per molecule
    __global const Transform* transforms, // - Precomputed transforms
    int nmols,                            // Total molecules (e.g., 54)
    float max_clash_penalty,              // Threshold for early exit
    __local float4* local_replica,        // Dynamically sized local memory:
    __local float* local_scores,          // Dynamically sized local memory:
    __local float* local_min_dist,        // Dynamically sized local memory:
    __global float* results,              // - Output array (clash sums)
    __global float* results_min           // - Output array (min distances)
) {
    // 1. Identify Workgroup (Configuration) and Thread
    int conf_id = get_group_id(0);
    int lid     = get_local_id(0);
    int wg_size = get_local_size(0);

    // Early exit flag shared across the workgroup
    __local int abort_flag;
    if (lid == 0) abort_flag = 0;
    barrier(CLK_LOCAL_MEM_FENCE);

    // Locate the start of the transforms for THIS configuration
    int conf_transform_offset = conf_id * nmols;

    // Load Transform for Molecule 0 (The reference molecule)
    Transform t0 = transforms[conf_transform_offset + 0];

    float my_score = 0.0f;
    float my_min_dist2 = 1e20f;

    // 2. OUTER LOOP: Tiles of Molecule 0
    for (int t_m0 = 0; t_m0 < natoms; t_m0 += wg_size) {
        
        int i = t_m0 + lid;
        bool valid_i = (i < natoms);
        float4 my_atom;
        
        // Transform the reference atom and keep it in ultra-fast private registers
        if (valid_i) {
            my_atom = apply_transform(base_atoms[i], t0);
        }

        // 3. MIDDLE LOOP: Iterate over all replicas (skip 0, which is the reference itself)
        for (int k = 1; k < nmols; k++) {
            
            Transform tk = transforms[conf_transform_offset + k];

            // 4. INNER LOOP: Tiles of Replica K
            for (int t_rep = 0; t_rep < natoms; t_rep += wg_size) {
                
                int j = t_rep + lid;
                
                // Cooperatively load transformed Replica K into local memory
                if (j < natoms) {
                    local_replica[lid] = apply_transform(base_atoms[j], tk);
                }
                barrier(CLK_LOCAL_MEM_FENCE); // Wait for all threads to finish loading tile

                // Compute pairwise interactions against the cached tile
                if (valid_i) {
                    int num_j = min(wg_size, natoms - t_rep);
                    
                    for (int lj = 0; lj < num_j; lj++) {
                        float4 partner = local_replica[lj];
                        
                        // Compute squared distance (saves square root operation if no clash)
                        float dx = my_atom.x - partner.x;
                        float dy = my_atom.y - partner.y;
                        float dz = my_atom.z - partner.z;
                        float dist_sq = dx*dx + dy*dy + dz*dz;
                        
                        if (dist_sq < my_min_dist2) {
                            my_min_dist2 = dist_sq;
                        }
                        
                        float r_sum = my_atom.w + partner.w;
                        float r_sum_sq = r_sum * r_sum;

                        // Check for collision
                        if (dist_sq < r_sum_sq) {
                            float dist = sqrt(dist_sq);
                            float overlap = r_sum - dist;
                            my_score += overlap * overlap; // Or your specific penalty function
                        }
                    }
                }
                barrier(CLK_LOCAL_MEM_FENCE); // Wait for compute to finish before loading next tile
            } // end tile replica
            
            // Check for Workgroup Early Exit
            // If just ONE thread sees a massive clash, the whole configuration is dead.
            if (my_score > max_clash_penalty) {
                abort_flag = 1; 
            }
            barrier(CLK_LOCAL_MEM_FENCE);
            if (abort_flag) break;

        } // end replica
        
        if (abort_flag) break;
        
    } // end m0 tile

    // 5. PARALLEL REDUCTION: Sum up scores from all threads in the workgroup
    local_scores[lid] = my_score;
    local_min_dist[lid] = my_min_dist2;
    barrier(CLK_LOCAL_MEM_FENCE);

    for (int stride = wg_size / 2; stride > 0; stride >>= 1) {
        if (lid < stride) {
            local_scores[lid] += local_scores[lid + stride];
            if (local_min_dist[lid + stride] < local_min_dist[lid]) {
                local_min_dist[lid] = local_min_dist[lid + stride];
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // 6. Write Final Result
    if (lid == 0) {
        if (abort_flag) {
            results[conf_id] = 999999.0f; // Write a massive penalty for aborted configs
            results_min[conf_id] = 0.0f;
        } else {
            results[conf_id] = local_scores[0]; // Write the total accumulated score
            results_min[conf_id] = sqrt(local_min_dist[0]);
        }
    }
}