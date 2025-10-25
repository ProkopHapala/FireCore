// Assume the helper functions above are in the same file.



// OpenCL Kernel Code (.cl file)

// The only allowed structure for a 3x3 matrix.
// The __attribute__((packed)) ensures it has no padding.
typedef struct __attribute__ ((packed)){
    float4 a; // Row 0
    float4 b; // Row 1
    float4 c; // Row 2
} cl_Mat3;

// --- Helper Functions ---

// Multiplies two quaternions (q2 * q1). No changes needed.
inline float4 quat_mult(float4 q1, float4 q2) {
    return (float4)(
        q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
        q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
        q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
        q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
    );
}

// Rotates a vector by a matrix, using the cl_Mat3 structure.
inline float4 rotate_vec_by_matrix(const float4 v, __local const cl_Mat3* R) {
    return (float4)(
        dot(R->a.xyz, v.xyz),
        dot(R->b.xyz, v.xyz),
        dot(R->c.xyz, v.xyz),
        0.0f
    );
}

// if (lid == 0) R.a = (float4)(1.0f - (yy + zz), xy - wz,           xz + wy,          0.0f);
// if (lid == 1) R.b = (float4)(xy + wz,          1.0f - (xx + zz), yz - wx,          0.0f);
// if (lid == 2) R.c = (float4)(xz - wy,          yz + wx,          1.0f - (xx + yy), 0.0f);
inline float3 quat_to_a( float4 q ){  return(float3)(1.0f - (q.y*q.y + q.z*q.z), q.x*q.y - q.z*q.w,           q.x*q.z + q.y*q.w         );}
inline float3 quat_to_b( float4 q ){  return(float3)(q.x*q.y + q.z*q.w,          1.0f - (q.x*q.x + q.z*q.z), q.y*q.z - q.x*q.w          );}
inline float3 quat_to_c( float4 q ){  return(float3)(q.x*q.z - q.y*q.w,          q.y*q.z + q.x*q.w,          1.0f - (q.x*q.x + q.y*q.y) );}




// Assume the struct and helpers above are in the same file.

#define WORKGROUP_SIZE     32
#define MAX_ATOMS_PER_BODY 128
#define ATOMS_PER_THREAD   4

__kernel
void rigid_body_dynamics_kernel(
    // --- Inputs ---
    __global const float4*  pos_in,              // Center of mass positions
    __global const float4*  quats_in,            // Orientations (quaternions)
    __global const float4*  lin_mom_in,          // Linear momenta
    __global const float4*  ang_mom_in,          // Angular momenta
    __global const float*   mass,                 // Mass of each body
    __global const float*   inv_mass,             // Inverse mass of each body
    __global const cl_Mat3* I_body_inv,         // Body-space inverse inertia tensor
    __global const float4*  all_atom_pos_body,   // Body-space positions of atoms

    // --- Outputs ---
    __global float4* pos_out,
    __global float4* quats_out,
    __global float4* lin_mom_out,
    __global float4* ang_mom_out,
    __global float4* all_atom_pos_world,        // Final world-space positions of atoms

    // --- Parameters ---
    const int num_atoms,
    const int num_integration_steps,
    const float dt
) {
    // --- 1. Initialization and Data Loading ---

    const int gid   = get_group_id(0);
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);

    // Local memory for the state of the single rigid body this workgroup processes.
    __local float4 local_pos;
    __local float4 local_quat;
    __local float4 local_lin_mom;
    __local float4 local_ang_mom;
    __local float  local_inv_mass;
    
    __local cl_Mat3 R; // World-space rotation matrix
    __local float4 local_torque[WORKGROUP_SIZE];
    __local float4 local_force [WORKGROUP_SIZE];

    __private float4 private_pos[ATOMS_PER_THREAD];

    // Thread 0 loads the rigid body's state into fast local memory.
    if (lid == 0) {
        local_pos     = pos_in    [gid];
        local_quat    = quats_in  [gid];
        local_lin_mom = lin_mom_in[gid];
        local_ang_mom = ang_mom_in[gid];
        local_inv_mass = inv_mass [gid];
    }

    // All threads cooperate to load the body-space atomic positions.
    //for (int i = lid; i < num_atoms; i += lsize) {local_atom_pos_body[i] = all_atom_pos_body[gid * MAX_ATOMS_PER_BODY + i];}

    for (int i=0; i<ATOMS_PER_THREAD; i++) {
        int atom_idx = lid+i*lsize;
        if(atom_idx < num_atoms){ private_pos[i] = all_atom_pos_body[ gid*MAX_ATOMS_PER_BODY + atom_idx]; }
    }

    // --- 2. Main Integration Loop ---
    for (int step = 0; step < num_integration_steps; ++step) {
        // --- initialized the rotation matrix
        if (lid == 0) R.a = (float4){ quat_to_a(local_quat), 0.f };
        if (lid == 1) R.b = (float4){ quat_to_b(local_quat), 0.f };
        if (lid == 2) R.c = (float4){ quat_to_c(local_quat), 0.f };
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- Step C: Update atom positions and calculate forces/torques ---
        float4 total_torque = (float4)(0.0f);
        float4 total_force  = (float4)(0.0f); // If you were also calculating this

        // --- evaluate position dependent forces on atoms in world coordinates
        for (int i=0; i<ATOMS_PER_THREAD; i+=lsize) {
            int atom_idx = lid+i*lsize;
            if(atom_idx >= num_atoms){ break; }
            // 1. Update world position of this atom relative to center of mass
            float4 r_world = rotate_vec_by_matrix(private_pos[i], &R);
            float4 p_world = local_pos + r_world;
            // 2. ==========================================================
            //    == YOUR FORCE CALCULATION LOGIC GOES HERE ==
            //    Calculate the force 'f' acting on the atom at p_world.
            // ==========================================================
            float4 f = (float4)(0.1f * r_world.y, 0.0f, 0.0f, 0.0f); // Placeholder spinning force
            total_torque += cross(r_world, f);
            total_force  += f;
        }

        // --- Step D: Reduce torque and force across the workgroup ---
        local_torque[lid] = total_torque;
        local_force [lid] = total_force;
        const int stride = WORKGROUP_SIZE/4;
        barrier(CLK_LOCAL_MEM_FENCE);
        const int lid_ = lid & (stride-1);
        if ( lid_==0 ){ // threads divisible by stride
            float4 tq = local_torque[lid+1];
            for(int i=2; i<stride; i++){ tq+=local_torque[lid+i]; }
            local_torque[lid]+=tq;
        }else if ( lid_==1 ) {
            float4 f = local_force[lid+1];
            for(int i=2; i<stride; i++){ f+=local_force[lid+i]; }
            local_force[lid]+=f;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // --- Step E: Update momenta (final part of leapfrog) ---
        if       ( lid == 0 ){
            float4 final_force = local_force[0]  +local_force[stride]+local_force[stride*2]+local_force[stride*3]; 
            local_lin_mom += final_force * dt;
            local_pos += local_lin_mom * local_inv_mass * dt;
        }else if ( lid ==1  ){
            float4 final_torque = local_torque[0]+local_torque[stride]+local_torque[stride*2]+local_torque[stride*3];
            local_ang_mom  += final_torque * dt;
            float4 omega    = local_ang_mom; // Placeholder! You need ω = I⁻¹L
            float4 q_deriv  = quat_mult(local_quat, (float4)(omega.x, omega.y, omega.z, 0.0f));
            local_quat     += 0.5f * dt * q_deriv;
            local_quat      = normalize(local_quat);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- store output atom positions (if needed) ---
    for (int i=0; i<ATOMS_PER_THREAD; i+=lsize) {
        int atom_idx = lid+i*lsize;
        if(atom_idx >= num_atoms){ break; }
        float4 r_world = rotate_vec_by_matrix(private_pos[i], &R);
        float4 p_world = local_pos + r_world;
        all_atom_pos_world[gid * MAX_ATOMS_PER_BODY + atom_idx] = p_world;
    }
    // --- store output body state ---
    if (lid == 0) {
        pos_out    [gid] = local_pos;
        quats_out  [gid] = local_quat;
        lin_mom_out[gid] = local_lin_mom;
        ang_mom_out[gid] = local_ang_mom;
    }
}