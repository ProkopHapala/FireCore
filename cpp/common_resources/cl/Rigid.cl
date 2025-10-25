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
    __global const int*      mols,        // [nworkgroups] {i0} starting index of first atom in each molecule in the bbuffer
    __global       float4*   poss,        // [nworkgroups] {x,y,z,q} Center of mass positions, q used for evaluation of forces
    __global       float4*   qrots,       // [nworkgroups] Orientations (quaternions)
    __global       float4*   vposs,       // [nworkgroups] Linear  momenta
    __global       float4*   vrots,       // [nworkgroups] Angular momenta
    __global const cl_Mat3*  I_body_inv,  // [nworkgroups] Body-space inverse inertia tensor
    __global const float4*   apos_body,   // [natomsTot] Body-space positions of atoms
    __global       float4*   apos_world,  // [natomsTot] World  Final world-space positions of atoms
    // --- Parameters ---
    const int   natoms, // number of atoms in the system
    const int   niter,  // number of integration steps
    const float dt      // time step
) {
    // --- 1. Initialization and Data Loading ---

    const int gid   = get_group_id(0);
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);

    // Local memory for the state of the single rigid body this workgroup processes.
    __local float4 pos;
    __local float4 qrot;
    __local float4 vpos;
    __local float4 vrot;
    __local float  inv_mass;
    __local cl_Mat3 R; // World-space rotation matrix
    __local cl_Mat3 Iinv_body; // Body-space inverse inertia tensor
    __local float4 Ltorq [WORKGROUP_SIZE];
    __local float4 Lforce[WORKGROUP_SIZE];

    __private float4 private_pos[ATOMS_PER_THREAD];

    // Thread 0 loads the rigid body's state into fast local memory.
    if (lid == 0) {
        pos      = poss   [gid];
        qrot     = qrots  [gid];
        vpos     = vposs  [gid];
        vrot     = vrots  [gid];
        inv_mass = 1.0f;
        Iinv_body.a = I_body_inv[gid].a;
        Iinv_body.b = I_body_inv[gid].b;
        Iinv_body.c = I_body_inv[gid].c;
        //printf("GPU gid=%d pos(%6e,%6e,%6e,%6e) qrot(%6e,%6e,%6e,%6e)|%6e| vpos(%6e,%6e,%6e,%6e) vrot(%6e,%6e,%6e,%6e)\n", 
        //    gid, pos.x, pos.y, pos.z, pos.w, qrot.x, qrot.y, qrot.z, qrot.w, length(qrot), vpos.x, vpos.y, vpos.z, vpos.w, vrot.x, vrot.y, vrot.z, vrot.w);
        printf("GPU gid=%d pos(%6e,%6e,%6e,%6e) qrot(%6e,%6e,%6e,%6e)|%6e| \n", gid, pos.x, pos.y, pos.z, pos.w, qrot.x, qrot.y, qrot.z, qrot.w, length(qrot));
    }

    // All threads cooperate to load the body-space atomic positions.
    //for (int i = lid; i < num_atoms; i += lsize) {local_atom_pos_body[i] = all_atom_pos_body[gid * MAX_ATOMS_PER_BODY + i];}

    const int i0 = mols[gid];
    const int na = mols[gid+1]-i0;
    for (int i=0; i<ATOMS_PER_THREAD; i++) {
        int atom_idx = lid+i*lsize;
        if(atom_idx < na){ private_pos[i] = apos_body[ i0 + atom_idx]; }
    }

    // --- 2. Main Integration Loop ---
    for (int step = 0; step < niter; ++step) {
        // --- initialized the rotation matrix
        if (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
        if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
        if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- Step C: Update atom positions and calculate forces/torques ---
        float4 total_torque = (float4)(0.0f);
        float4 total_force  = (float4)(0.0f); // If you were also calculating this

        // --- evaluate position dependent forces on atoms in world coordinates
        for (int i=0; i<ATOMS_PER_THREAD; i++) {
            int atom_idx = lid+i*lsize;
            if(atom_idx >= na){ break; }
            // 1. Update world position of this atom relative to center of mass
            float4 r_world = rotate_vec_by_matrix(private_pos[i], &R);
            float4 p_world = pos + r_world;
            // 2. ==========================================================
            //    == YOUR FORCE CALCULATION LOGIC GOES HERE ==
            //    Calculate the force 'f' acting on the atom at p_world.
            // ==========================================================
            float4 f = (float4)(0.0f, 0.0f, 0.0f, 0.0f); // Placeholder spinning force
            total_torque += cross(r_world, f);
            total_force  += f;
        }

        // --- Step D: Reduce torque and force across the workgroup ---
        Ltorq [lid] = total_torque;
        Lforce[lid] = total_force;
        const int stride = WORKGROUP_SIZE/4;
        barrier(CLK_LOCAL_MEM_FENCE);
        const int lid_ = lid & (stride-1);
        if ( lid_==0 ){ // threads divisible by stride
            float4 tq = Ltorq[lid+1];
            for(int i=2; i<stride; i++){ tq+=Ltorq[lid+i]; }
            Ltorq[lid]+=tq;
        }else if ( lid_==1 ) {
            float4 f = Lforce[lid+1];
            for(int i=2; i<stride; i++){ f+=Lforce[lid+i]; }
            Lforce[lid]+=f;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        // --- Step E: Update momenta (final part of leapfrog) ---
        if       ( lid == 0 ){
            float4 f   = Lforce[0]  +Lforce[stride]+Lforce[stride*2]+Lforce[stride*3]; 
            vpos.xyz  += f.xyz * dt;
            pos.xyz   += vpos.xyz * inv_mass * dt;
        }else if ( lid ==1  ){
            float4 tq      = Ltorq[0]+Ltorq[stride]+Ltorq[stride*2]+Ltorq[stride*3];
            //vrot.xyz      += tq.xyz   * dt;
            float3 omega   = vrot.xyz * dt*0.1; // Placeholder! You need ω = I⁻¹L
            // float3 omega3 = (float3)(
            //     dot(Iinv_body.a.xyz, vrot.xyz),
            //     dot(Iinv_body.b.xyz, vrot.xyz),
            //     dot(Iinv_body.c.xyz, vrot.xyz)
            // );
            //float4 dq     = quat_mult(qrot, (float4)(omega.x, omega.y, omega.z, 0.0f));
            //qrot          += 0.5f * dt * dq;
            qrot    = quat_mult(qrot, (float4)(omega.x, omega.y, omega.z, 0.0f));
            qrot    = normalize(qrot);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- store output atom positions (if needed) ---
    for (int i=0; i<ATOMS_PER_THREAD; i++) {
        int atom_idx = lid+i*lsize;
        if(atom_idx >= na){ break; }
        float4 r_world = rotate_vec_by_matrix(private_pos[i], &R);
        float4 p_world = pos + r_world;
        apos_world[i0 + atom_idx] = p_world;
    }
    // --- store output body state ---
    if (lid == 0) {
        poss   [gid] = pos;
        qrots  [gid] = qrot;
        vposs  [gid] = vpos;
        vrots  [gid] = vrot;
    }
}