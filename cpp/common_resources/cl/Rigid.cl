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

float2 sinc_div_r2_taylor(float r2){
    // series up to r^6 terms (i.e. up to r2^3)
    // s = sin(r)/r      = 1   - r^2/6  + r^4/120 - r^6/5040
    // c = (1-cos r)/r^2 = 1/2 - r^2/24 + r^4/720 - r^6/40320
    const float s = 1.0f + r2 * ( (-1.0f/6.0f)  + r2 * ( (1.0f/120.0f) + r2 * (-1.0f/5040.0f  ) ) );
    const float c = 0.5f + r2 * ( (-1.0f/24.0f) + r2 * ( (1.0f/720.0f) + r2 * (-1.0f/40320.0f ) ) );
    return (float2){s, c};
}

inline float2 quat_factors_taylor(float r2){
    // Series up to r^6 terms (i.e., up to r2^3)
    // s = sin(r/2)/r = 1/2 - r2/48 + r2^2/3840 - r2^3/645120
    // c = cos(r/2)   = 1   - r2/8  + r2^2/384  - r2^3/46080
    const float s = 0.5f + r2 * ((-1.0f/48.0f)  + r2 * ((1.0f/3840.0f) + r2 * (-1.0f/645120.0f)));
    const float c = 1.0f + r2 * ((-1.0f/8.0f)   + r2 * ((1.0f/384.0f)  + r2 * (-1.0f/46080.0f)));
    return (float2)(s, c);
}

float4 qrot_omega_taylor( float4 qrot, float3 omega){
    const float r2 = dot(omega,omega);
    const float2 sc = quat_factors_taylor(r2);
    return quat_mult(qrot, (float4)(omega * sc.x, sc.y) );
}

float4 make_qrot_taylor(  float3 omega){
    const float r2 = dot(omega,omega);
    const float2 sc = quat_factors_taylor(r2);
    return (float4)(omega * sc.x, sc.y);
}

inline float4 make_qrot(float3 omega){
    const float angle = length(omega);
    if(angle < 1e-12f) return (float4)(0.0f, 0.0f, 0.0f, 1.0f);
    const float3 axis = omega / angle;
    const float s = sin(0.5f * angle);
    float c = cos(0.5f * angle);
    return (float4)(axis * s, c);
}

float4 qrot_omega( float4 qrot, float3 omega){
    const float4 dq  = make_qrot(omega);
    return quat_mult(qrot, dq);
}



// Rotates a vector by a matrix, using the cl_Mat3 structure.
inline float3 rotate_vec_by_matrix(const float3 v, __local const cl_Mat3* R) {
    return (float3)(
        dot(R->a.xyz, v.xyz),
        dot(R->b.xyz, v.xyz),
        dot(R->c.xyz, v.xyz)
    );
}

inline float3 rotate_vec_by_matrix_T(const float3 v, __local const cl_Mat3* R) {  return R->a.xyz*v.x + R->b.xyz*v.y + R->c.xyz*v.z;}

inline float3 quat_to_a( float4 q ){  return(float3)(1.0f-2.0f*(q.y*q.y + q.z*q.z),      2.0f*(q.x*q.y - q.z*q.w),      2.0f*(q.x*q.z + q.y*q.w));}
inline float3 quat_to_b( float4 q ){  return(float3)(     2.0f*(q.x*q.y + q.z*q.w), 1.0f-2.0f*(q.x*q.x + q.z*q.z),      2.0f*(q.y*q.z - q.x*q.w));}
inline float3 quat_to_c( float4 q ){  return(float3)(     2.0f*(q.x*q.z - q.y*q.w),      2.0f*(q.y*q.z + q.x*q.w), 1.0f-2.0f*(q.x*q.x + q.y*q.y));}


// Assume the struct and helpers above are in the same file.

#define WORKGROUP_SIZE     32
#define MAX_ATOMS_PER_BODY 128
#define ATOMS_PER_THREAD   4

__kernel
void rigid_body_dynamics_kernel(
    // --- Inputs ---
    __global const int*      mols,        // [nworkgroups] {i0} starting index of first atom in each molecule in the bbuffer
    __global       float4*   poss,        // [nworkgroups] {x,y,z,m} Center of mass positions, m is mass of the whole molecule
    __global       float4*   qrots,       // [nworkgroups] Orientations (quaternions)
    __global       float4*   vposs,       // [nworkgroups] Linear  momenta
    __global       float4*   vrots,       // [nworkgroups] Angular momenta
    __global const cl_Mat3*  I_body_inv,  // [nworkgroups] Body-space inverse inertia tensor
    __global const float4*   apos_body,   // [natomsTot] {x,y,z,q} Body-space positions of atoms, anch charge of each atom 
    __global       float4*   apos_world,  // [natomsTot] World  Final world-space positions of atoms
    __global const float4*   anchors,     // [natomsTot] {x,y,z,K} per-atom anchors (K<0 => free)
    // --- Parameters ---
    const int   natoms,   // number of atoms in the system
    const int   niter,    // number of integration steps
    const float dt,       // time step
    const float3  Efield  // Homogeneous electric field
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

    //__private float4 private_pos[ATOMS_PER_THREAD];

    const int ia0 = mols[gid];
    const int na  = mols[gid+1]-ia0;
    // ---- This optimization we need to benchmark later, it is not obvious if the private memory registers lost are worth it
    // for (int i=0; i<ATOMS_PER_THREAD; i++) {
    //     int atom_idx = lid+i*lsize;
    //     if(atom_idx < na){ private_pos[i] = apos_body[ i0 + atom_idx]; }
    // }

    // Thread 0 loads the rigid body's state into fast local memory.
    if (lid == 0) {
        pos      = poss   [gid];
        qrot     = qrots  [gid];  qrot=normalize(qrot);
        vpos     = vposs  [gid];
        vrot     = vrots  [gid];
        inv_mass = 1.0f;
        Iinv_body.a = I_body_inv[gid].a;
        Iinv_body.b = I_body_inv[gid].b;
        Iinv_body.c = I_body_inv[gid].c;
        //printf("GPU gid=%d pos(%6e,%6e,%6e,%6e) qrot(%6e,%6e,%6e,%6e)|%6e| vpos(%6e,%6e,%6e,%6e) vrot(%6e,%6e,%6e,%6e)\n", 
        //    gid, pos.x, pos.y, pos.z, pos.w, qrot.x, qrot.y, qrot.z, qrot.w, length(qrot), vpos.x, vpos.y, vpos.z, vpos.w, vrot.x, vrot.y, vrot.z, vrot.w);
        //printf("GPU gid=%d pos(%6e,%6e,%6e,%6e) qrot(%6e,%6e,%6e,%6e)|%6e| \n", gid, pos.x, pos.y, pos.z, pos.w, qrot.x, qrot.y, qrot.z, qrot.w, length(qrot));
        //for(int i=0; i<na; i++){ const float4 a=anchors[ia0+i];printf("GPU atom %d anchor(%6e,%6e,%6e,%6e)\n",i,a.x,a.y,a.z,a.w); }
        //for(int i=0; i<na; i++){ const float4 a=apos_body[ia0+i];printf("GPU atom %d xyz|q(%6e,%6e,%6e|%6e)\n",i,a.x,a.y,a.z,a.w); }
        //printf("GPU Efield (%6e,%6e,%6e)\n", Efield.x, Efield.y, Efield.z);
        //for(int i=0; i<na; i++){ const float4 a=apos_body[ia0+i]; const float4 b=anchors[ia0+i]; printf("GPU atom %d xyz|q(%6e,%6e,%6e|%6e) anchor(%6e,%6e,%6e|%6e)\n",i,a.x,a.y,a.z,a.w,b.x,b.y,b.z,b.w); }
    }

    // --- 2. Main Integration Loop ---
    for (int step = 0; step < niter; ++step) {
        // --- initialized the rotation matrix
        if      (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
        else if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
        else if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- Step C: Update atom positions and calculate forces/torques ---
        float4 total_torque = (float4)(0.0f);
        float4 total_force  = (float4)(0.0f); // If you were also calculating this

        // --- evaluate position dependent forces on atoms in world coordinates
        for (int i=0; i<ATOMS_PER_THREAD; i++) {
            int atom_idx = lid+i*lsize;
            if(atom_idx >= na){ break; }
            // 1. Update world position of this atom relative to center of mass
            //float4 r_world = rotate_vec_by_matrix(private_pos[i], &R);
            //float4 r_world = rotate_vec_by_matrix_T(private_pos[i], &R);        // Using private memory? We test it later
            float4 p_body  = apos_body[ia0+atom_idx];
            float3 p_world = pos.xyz + rotate_vec_by_matrix(p_body.xyz, &R); 

            // ------ Force evaluation ------
            float4 f = (float4)(0.0f, 0.0f, 0.0f, 0.0f); // Placeholder spinning force
            {// constant electronstatic force
                f.xyz += Efield.xyz*p_body.w; // E*q   // p_body.w is the charge of the atom
            }
            float4 anchor   = anchors[ia0+atom_idx];
            if(anchor.w > 0.0f){ // anchors
                float3 d  = p_world.xyz - anchor.xyz;
                float3 fa = d * -anchor.w;
                f.xyz    += fa;
                printf("ia %3i l_anchor=%8e fa(%8e,%8e,%8e)  anchor(%8e,%8e,%8e|%8e)\n", atom_idx, length(d), fa.x,fa.y,fa.z, anchor.x,anchor.y,anchor.z,anchor.w);
            }

            // ---- Project force to rigid body
            float3 tq = cross(p_world.xyz, f.xyz);
            //if(dot(f.xyz,f.xyz) > 1.e-16 ){ printf("ia %3i f(%8e,%8e,%8e) tq(%8e,%8e,%8e) l_anchor=%8e anchor(%8e,%8e,%8e|%8e)\n", atom_idx, f.x,f.y,f.z, tq.x,tq.y,tq.z, anchor.x,anchor.y,anchor.z,anchor.w); }
            total_torque.xyz += tq;
            total_force .xyz += f.xyz;
        }

        // --- Reduce torque and force across the workgroup ---
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
            const int id=lid-1;
            float4 f = Lforce[id];
            for(int i=2; i<stride; i++){ f+=Lforce[id+i]; }
            Lforce[id]+=f;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- Update Ridig body pose and momentum using Leapfrog ---
        if       ( lid == 0 ){
            float4 f   = Lforce[0]+Lforce[stride]+Lforce[stride*2]+Lforce[stride*3]; 
            vpos.xyz*=0.9f;
            vpos.xyz  += f.xyz               * dt;
            pos.xyz   += vpos.xyz * inv_mass * dt;
        //}else if ( lid ==1  ){
            float4 tq  = Ltorq[0]+Ltorq[stride]+Ltorq[stride*2]+Ltorq[stride*3];
            //printf( "force(%8e,%8e,%8e) tq(%8e,%8e,%8e)\n", f.x,f.y,f.z, tq.x,tq.y,tq.z );
            vrot.xyz*=0.9f;
            float      Iinv = 1.f/25.f;
            vrot.xyz   += tq.xyz  *dt*Iinv; // for now we want to keep vrot constant
            float3 ot   = vrot.xyz*dt;      // Placeholder! You need ω = I⁻¹L
            //float4 dq_= make_qrot       ( ot ); 
            float4 dq   = make_qrot_taylor( ot ); 
            //printf( "w*dt(%8e,%8e,%8e) dq_sincos(%8e,%8e,%8e,%8e) dq_taylor(%8e,%8e,%8e,%8e) \n", ot.x,ot.y,ot.z, dq_.x,dq_.y,dq_.z, dq_.w, dq.x,dq.y,dq.z,dq.w );
            qrot        = quat_mult(qrot, dq );
            //qrot       = qrot_omega_taylor(qrot, vrot.xyz*dt*0.1f);
            qrot      = normalize(qrot);
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if      (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
    else if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
    else if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
    barrier(CLK_LOCAL_MEM_FENCE);
    //if(lid==0){ printf("GPU gid=%d R.a(%6e,%6e,%6e)|%6e| R.b(%6e,%6e,%6e)|%6e| R.c(%6e,%6e,%6e)|%6e|\n",   gid, R.a.x,R.a.y,R.a.z,length(R.a),  R.b.x,R.b.y,R.b.z,length(R.b),   R.c.x,R.c.y,R.c.z,length(R.c)); }
    // --- store output atom positions (if needed) ---
    for (int i=0; i<ATOMS_PER_THREAD; i++) {
        int atom_idx = lid+i*lsize;
        if(atom_idx >= na){ break; }
        //float4 r_world = rotate_vec_by_matrix(private_pos[i], &R);
        //float4 r_world = rotate_vec_by_matrix_T(private_pos[i], &R);
        int ia = ia0+atom_idx;
        float4 p_body  = apos_body[ia];
        float3 p_world = pos.xyz + rotate_vec_by_matrix(p_body.xyz, &R); 
        apos_world[ia] = (float4){p_world, 0.f};
    }
    // --- store output body state ---
    if (lid == 0) {
        poss   [gid] = pos;
        qrots  [gid] = qrot;
        vposs  [gid] = vpos;
        vrots  [gid] = vrot;
    }
}