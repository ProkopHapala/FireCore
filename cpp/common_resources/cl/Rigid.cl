typedef struct __attribute__ ((packed)){
    float4 a; // Row 0
    float4 b; // Row 1
    float4 c; // Row 2
} cl_Mat3;

#ifndef RIGID_DBG
#define RIGID_DBG 0
#endif

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

inline float3 mat3_dot(const cl_Mat3 M, const float3 v){
    return (float3)( dot(M.a.xyz, v), dot(M.b.xyz, v), dot(M.c.xyz, v) );
}

inline float3 mat3_dot_T(const cl_Mat3 M, const float3 v){
    return M.a.xyz*v.x + M.b.xyz*v.y + M.c.xyz*v.z;
}

inline int modulo(const int i, const int m) {
    int result = i % m;
    if (result < 0) result += m;
    return result;
}

inline int4 make_inds_pbc(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int4)(0, 1,   2,   3  ); }
        case 1: { return (int4)(0, 1,   2,   3-n); }
        case 2: { return (int4)(0, 1,   2-n, 3-n); }
        case 3: { return (int4)(0, 1-n, 2-n, 3-n); }
    }
    return (int4)(-100, -100, -100, -100);
}

inline int4 choose_inds_pbc(const int i, const int n, const __local int4* iqs) {
    if (i >= (n-3)) {
        const int ii = i + 4 - n;
        return iqs[ii];
    }
    return (int4)(0, +1, +2, +3);
}

inline float4 basis(float u) {
    const float inv6 = 1.0f / 6.0f;
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        inv6 * t * t * t,
        inv6 * (3.0f * u2 * (u - 2.0f) + 4.0f),
        inv6 * (3.0f * u * (1.0f + u - u2) + 1.0f),
        inv6 * u2 * u
    );
}

inline float4 dbasis(float u) {
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        -0.5f * t * t,
        0.5f * (3.0f * u2 - 4.0f * u),
        0.5f * (-3.0f * u2 + 2.0f * u + 1.0f),
        0.5f * u2
    );
}

inline float2 fe1Dcomb(__global const float4* E, const float4 C, const float4 p, const float4 d) {
    const float4 cs = (float4)(dot(C, E[0]), dot(C, E[1]), dot(C, E[2]), dot(C, E[3]));
    return (float2)(dot(p, cs), dot(d, cs));
}

inline float3 fe2d_comb(int nz, __global const float4* E, int4 di, const float4 C, const float4 pz, const float4 dz, const float4 by, const float4 dy) {
    const float2 fe0 = fe1Dcomb(E + di.x, C, pz, dz);
    const float2 fe1 = fe1Dcomb(E + di.y, C, pz, dz);
    const float2 fe2 = fe1Dcomb(E + di.z, C, pz, dz);
    const float2 fe3 = fe1Dcomb(E + di.w, C, pz, dz);
    return (float3)(
        fe0.x * dy.x + fe1.x * dy.y + fe2.x * dy.z + fe3.x * dy.w,
        fe0.y * by.x + fe1.y * by.y + fe2.y * by.z + fe3.y * by.w,
        fe0.x * by.x + fe1.x * by.y + fe2.x * by.z + fe3.x * by.w
    );
}

inline float4 fe3d_pbc_comb(const float3 u, const int3 n, __global const float4* Es, const float4 PLQH, __local const int4* xqis, __local int4* yqis) {
    int ix = (int)u.x;
    int iy = (int)u.y;
    int iz = (int)u.z;
    if (u.x < 0) ix--;
    if (u.y < 0) iy--;
    const float tx = u.x - ix;
    const float ty = u.y - iy;
    const float tz = u.z - iz;
    if ((iz < 1) || (iz >= n.z - 2)) return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    ix = modulo(ix-1, n.x);
    iy = modulo(iy-1, n.y);
    const int nyz = n.z * n.y;
    int4 qx = choose_inds_pbc(ix, n.x, xqis);
    const int4 qy = choose_inds_pbc(iy, n.y, yqis) * n.z;
    const float4 bz = basis(tz);
    const float4 dz = dbasis(tz);
    const float4 by = basis(ty);
    const float4 dy = dbasis(ty);
    const int i0 = (iz - 1) + n.z * (iy + n.y * ix);
    qx *= nyz;
    float3 E1 = fe2d_comb(n.z, Es + (i0 + qx.x), qy, PLQH, bz, dz, by, dy);
    float3 E2 = fe2d_comb(n.z, Es + (i0 + qx.y), qy, PLQH, bz, dz, by, dy);
    float3 E3 = fe2d_comb(n.z, Es + (i0 + qx.z), qy, PLQH, bz, dz, by, dy);
    float3 E4 = fe2d_comb(n.z, Es + (i0 + qx.w), qy, PLQH, bz, dz, by, dy);
    const float4 bx = basis(tx);
    const float4 dx = dbasis(tx);
    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))
    );
}

#define WORKGROUP_SIZE     32
#define MAX_ATOMS_PER_BODY 128
#define ATOMS_PER_THREAD   4

__kernel
void rigid_body_dynamics_kernel(
    __global const int*      mols,
    __global       float4*   poss,
    __global       float4*   qrots,
    __global       float4*   vposs,
    __global       float4*   vrots,
    __global const cl_Mat3*  I_body_inv,
    __global const float4*   apos_body,
    __global       float4*   apos_world,
    __global const float4*   anchors,
    const int   natoms,
    const int   niter,
    const float dt,
    const float3  Efield
) {
    const int gid   = get_group_id(0);
    const int lid   = get_local_id(0);
    const int lsize = get_local_size(0);
    __local float4 pos;
    __local float4 qrot;
    __local float4 vpos;
    __local float4 vrot;
    __local float  inv_mass;
    __local cl_Mat3 R;
    __local cl_Mat3 Iinv_body;
    __local float4 Ltorq [WORKGROUP_SIZE];
    __local float4 Lforce[WORKGROUP_SIZE];
    const int ia0 = mols[gid];
    const int na  = mols[gid+1]-ia0;
    if (lid == 0) {
        pos      = poss   [gid];
        qrot     = qrots  [gid];  qrot=normalize(qrot);
        vpos     = vposs  [gid];
        vrot     = vrots  [gid];
        inv_mass = (pos.w > 1e-8f) ? (1.0f / pos.w) : 1.0f;
        Iinv_body.a = I_body_inv[gid].a;
        Iinv_body.b = I_body_inv[gid].b;
        Iinv_body.c = I_body_inv[gid].c;
    }
    for (int step = 0; step < niter; ++step) {
        if      (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
        else if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
        else if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
        barrier(CLK_LOCAL_MEM_FENCE);
        float4 total_torque = (float4)(0.0f);
        float4 total_force  = (float4)(0.0f);
        for (int i=0; i<ATOMS_PER_THREAD; i++) {
            int atom_idx = lid+i*lsize;
            if(atom_idx >= na){ break; }
            float4 p_body  = apos_body[ia0+atom_idx];
            float3 r_world = rotate_vec_by_matrix(p_body.xyz, &R);
            float3 p_world = pos.xyz + r_world;
            float4 f = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
            f.xyz += Efield.xyz*p_body.w;
            float4 anchor   = anchors[ia0+atom_idx];
            if(anchor.w > 0.0f){
                float3 d  = p_world.xyz - anchor.xyz;
                float3 fa = d * -anchor.w;
                f.xyz    += fa;
            }
            float3 tq = cross(r_world, f.xyz);
            total_torque.xyz += tq;
            total_force .xyz += f.xyz;
        }
        Ltorq [lid] = total_torque;
        Lforce[lid] = total_force;
        const int stride = WORKGROUP_SIZE/4;
        barrier(CLK_LOCAL_MEM_FENCE);
        const int lid_ = lid & (stride-1);
        if ( lid_==0 ){
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
        if ( lid == 0 ){
            float3 f = (Lforce[0]+Lforce[stride]+Lforce[stride*2]+Lforce[stride*3]).xyz;
            float3 tq_world = (Ltorq[0]+Ltorq[stride]+Ltorq[stride*2]+Ltorq[stride*3]).xyz;
            float3 tq_body  = mat3_dot_T(R, tq_world);
            float3 alpha_body  = mat3_dot(Iinv_body, tq_body);
            float3 alpha_world = mat3_dot(R, alpha_body);
            vpos.xyz *= 0.90f;
            vrot.xyz *= 0.90f;
            vpos.xyz += f * (dt * inv_mass);
            vrot.xyz += alpha_world * dt;
            pos.xyz  += vpos.xyz * dt;
            float4 dq = make_qrot_taylor(vrot.xyz * dt);
            qrot = normalize(quat_mult(qrot, dq));
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if      (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
    else if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
    else if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i=0; i<ATOMS_PER_THREAD; i++) {
        int atom_idx = lid+i*lsize;
        if(atom_idx >= na){ break; }
        int ia = ia0+atom_idx;
        float4 p_body  = apos_body[ia];
        float3 p_world = pos.xyz + rotate_vec_by_matrix(p_body.xyz, &R); 
        apos_world[ia] = (float4){p_world, 0.f};
    }
    if (lid == 0) {
        poss   [gid] = pos;
        qrots  [gid] = qrot;
        vposs  [gid] = vpos;
        vrots  [gid] = vrot;
    }
}

__kernel
void rigid_body_gridff_kernel(
    __global const int*      mols,
    __global       float4*   poss,
    __global       float4*   qrots,
    __global       float4*   vposs,
    __global       float4*   vrots,
    __global const cl_Mat3*  I_body_inv,
    __global const float4*   apos_body,
    __global       float4*   apos_world,
    __global const float4*   atom_PLQ,
    __global const float4*   BsplinePLQ,
    __global       float4*   atom_force,
    __global       float4*   body_force,
    __global       float4*   body_torque,
    __global const float4*   anchors,
    const int4               grid_ns,
    const float4             grid_invStep,
    const float4             grid_p0,
    const float              dt,
    const float4             md_params,
    const int                niter
) {
    const int gid   = get_group_id(0);
    const int lid   = get_local_id(0);
    __local float4 pos;
    __local float4 qrot;
    __local float4 vpos;
    __local float4 vrot;
    __local float  inv_mass;
    __local cl_Mat3 R;
    __local cl_Mat3 Iinv_body;
    __local float4 Ltorq [WORKGROUP_SIZE];
    __local float4 Lforce[WORKGROUP_SIZE];
    __local int4 xqs[4];
    __local int4 yqs[4];
    const int ia0 = mols[gid];
    const int na  = mols[gid+1]-ia0;
    if      (lid < 4){ xqs[lid]   = make_inds_pbc(grid_ns.x, lid); }
    else if (lid < 8){ yqs[lid-4] = make_inds_pbc(grid_ns.y, lid-4); }
    if (lid == 0) {
        pos      = poss[gid];
        qrot     = normalize(qrots[gid]);
        vpos     = vposs[gid];
        vrot     = vrots[gid];
        inv_mass = (pos.w > 1e-8f) ? (1.0f / pos.w) : 1.0f;
        Iinv_body.a = I_body_inv[gid].a;
        Iinv_body.b = I_body_inv[gid].b;
        Iinv_body.c = I_body_inv[gid].c;
    }
    barrier(CLK_LOCAL_MEM_FENCE);
    const float3 inv_dg = grid_invStep.xyz;
    for (int step = 0; step < niter; ++step) {
        if      (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
        else if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
        else if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
        barrier(CLK_LOCAL_MEM_FENCE);
        float4 total_torque = (float4)(0.0f);
        float4 total_force  = (float4)(0.0f);
        for (int i=0; i<ATOMS_PER_THREAD; i++) {
            const int atom_idx = lid + i*WORKGROUP_SIZE;
            if(atom_idx >= na) break;
            const int ia = ia0 + atom_idx;
            const float4 p_body = apos_body[ia];
            const float3 r_world = rotate_vec_by_matrix(p_body.xyz, &R);
            const float3 p_world = pos.xyz + r_world;
            const float4 fe = fe3d_pbc_comb((p_world - grid_p0.xyz) * inv_dg, grid_ns.xyz, BsplinePLQ, atom_PLQ[ia], xqs, yqs);
            float3 f = fe.xyz * (-inv_dg);
            float4 anchor = anchors[ia];
            if(anchor.w > 0.0f){
                float3 d = p_world - anchor.xyz;
                f += d * -anchor.w;
            }
            total_force.xyz += f;
            total_torque.xyz += cross(r_world, f);
            apos_world[ia] = (float4)(p_world, fe.w);
            atom_force[ia] = (float4)(f, fe.w);
#if RIGID_DBG
            if((gid==0)&&(step==0)&&(atom_idx<4)){
                printf("RIGID_DBG atom %i p(%g,%g,%g) PLQ(%g,%g,%g,%g) fe(%g,%g,%g,%g)\n", atom_idx, p_world.x,p_world.y,p_world.z, atom_PLQ[ia].x,atom_PLQ[ia].y,atom_PLQ[ia].z,atom_PLQ[ia].w, f.x,f.y,f.z,fe.w);
            }
#endif
        }
        Ltorq[lid] = total_torque;
        Lforce[lid] = total_force;
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int stride = WORKGROUP_SIZE >> 1; stride > 0; stride >>= 1) {
            if (lid < stride) {
                Ltorq[lid]  += Ltorq [lid + stride];
                Lforce[lid] += Lforce[lid + stride];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        if (lid == 0) {
            const float lin_damp = md_params.x;
            const float ang_damp = md_params.y;
            const float force_scale = md_params.z;
            const float torque_scale = md_params.w;
            const float3 f = Lforce[0].xyz * force_scale;
            const float3 tq_world = Ltorq[0].xyz * torque_scale;
            body_force [gid] = (float4)(f, 0.0f);
            body_torque[gid] = (float4)(tq_world, 0.0f);
#if RIGID_DBG
            if(gid==0){
                printf("RIGID_DBG body f(%g,%g,%g) tq(%g,%g,%g) pos(%g,%g,%g)\n", f.x,f.y,f.z, tq_world.x,tq_world.y,tq_world.z, pos.x,pos.y,pos.z);
            }
#endif
            const float3 tq_body  = mat3_dot_T(R, tq_world);
            const float3 alpha_body  = mat3_dot(Iinv_body, tq_body);
            const float3 alpha_world = mat3_dot(R, alpha_body);
            vpos.xyz *= lin_damp;
            vrot.xyz *= ang_damp;
            vpos.xyz += f * (dt * inv_mass);
            vrot.xyz += alpha_world * dt;
            pos.xyz  += vpos.xyz * dt;
            qrot = normalize(quat_mult(qrot, make_qrot_taylor(vrot.xyz * dt)));
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if      (lid == 0) R.a = (float4){ quat_to_a(qrot), 0.f };
    else if (lid == 1) R.b = (float4){ quat_to_b(qrot), 0.f };
    else if (lid == 2) R.c = (float4){ quat_to_c(qrot), 0.f };
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i=0; i<ATOMS_PER_THREAD; i++) {
        const int atom_idx = lid + i*WORKGROUP_SIZE;
        if(atom_idx >= na) break;
        const int ia = ia0 + atom_idx;
        const float4 p_body = apos_body[ia];
        const float3 p_world = pos.xyz + rotate_vec_by_matrix(p_body.xyz, &R);
        apos_world[ia].xyz = p_world;
    }
    if (lid == 0) {
        poss [gid] = pos;
        qrots[gid] = qrot;
        vposs[gid] = vpos;
        vrots[gid] = vrot;
    }
}