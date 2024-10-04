
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

typedef struct {
    float4 a;
    float4 b;
    float4 c;
} Mat3;

inline int modulo(int i, int m) {
    int result = i % m;
    return (result < 0) ? (result + m) : result;
}

inline int4 make_inds_pbc(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int4)(0, 1,   2,   3  ); }
        case 1: { return (int4)(0, 1,   2,   3-n); }
        case 2: { return (int4)(0, 1,   2-n, 3-n); }
        case 3: { return (int4)(0, 1-n, 2-n, 3-n); }
    }
    return (int4)(-100, -100, -100, -100);
    // iqs[0] = (int4)(0, 1,   2,   3  );
    // iqs[1] = (int4)(0, 1,   2,   3-n);
    // iqs[2] = (int4)(0, 1,   2-n, 3-n);
    // iqs[3] = (int4)(0, 1-n, 2-n, 3-n);
}

inline int4 choose_inds_pbc(const int i, const int n, const int4* iqs) {
    if (i >= (n-3)) {
        const int ii = i + 4 - n;
        return iqs[ii];
    }
    return (int4)(0, +1, +2, +3);
}

inline int4 choose_inds_pbc_3( const int i, const int n, const int4* iqs ){
    if(i>=(n-3)){ 
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        const int4 d = iqs[ii];
        return (int4){ i+d.x, i+d.y, i+d.z, i+d.w }; 
    }
    return (int4){ i, i+1, i+2, i+3 };
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

    if ((iz < 1) || (iz >= n.z - 2)) {
        return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }

    ix = modulo(ix-1, n.x);
    iy = modulo(iy-1, n.y);

    const int nyz = n.z * n.y;
    // int4 qx = xqis[ix%4] * nyz;
    // int4 qy = yqis[iy%4] * n.z;

    int4 qx = choose_inds_pbc( ix, n.x, xqis );
    //const int4 qx = choose_inds_pbc( ix, n.x, xqis )*nyz;
    const int4 qy = choose_inds_pbc( iy, n.y, yqis )*n.z;

    const float4 bz = basis(tz);
    const float4 dz = dbasis(tz);
    const float4 by = basis(ty);
    const float4 dy = dbasis(ty);
    
    const int i0 = (iz - 1) + n.z * (iy + n.y * ix);

    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) n(%i,%i,%i) \n", u.x,u.y,u.z, ix,iy,iz, n.x,n.y,n.z );
    printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) qx(%i,%i,%i,%i) nyz=%i\n", u.x,u.y,u.z, ix,iy,iz, qx.x,qx.y,qx.z,qx.w, nyz );
    qx*=nyz;
    
    //return (float4){ 0.0f, 0.0f, 0.0f, dot(PLQH, Es[ i0 ])  };

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

__kernel void sample3D_comb(
    const float4 g0,
    const float4 dg,
    const int4 ng,
    __global const float4* Eg,
    const int n,
    __global const float4* ps,
    __global float4* fes,
    const float4 C
    //__global int4* xqs,
    //__global int4* yqs
) {
    const int iG = get_global_id(0);
    const int iL = get_local_id(0);
    if (iG >= n) return;

    __local int4 xqs[4];
    __local int4 yqs[4];
    if      (iL<4){             xqs[iL]=make_inds_pbc(ng.x,iL); }
    else if (iL<8){ int i=iL-4; yqs[i ]=make_inds_pbc(ng.y,i ); };
    barrier(CLK_LOCAL_MEM_FENCE);

    const float3 inv_dg = 1.0f / dg.xyz;


    if( iG==0 ){
        printf("GPU sample3D_comb() ng(%i,%i,%i) g0(%g,%g,%g) dg(%g,%g,%g) C(%g,%g,%g) \n", ng.x,ng.y,ng.z,   g0.x,g0.y,g0.z,   dg.x,dg.y,dg.z,   C.x,C.y,C.z );
        printf("GPU xqs[0](%i,%i,%i,%i) xqs[1](%i,%i,%i,%i) xqs[2](%i,%i,%i,%i) xqs[3](%i,%i,%i,%i)\n", xqs[0].x, xqs[0].y, xqs[0].z, xqs[0].w,   xqs[1].x, xqs[1].y, xqs[1].z, xqs[1].w,   xqs[2].x, xqs[2].y, xqs[2].z, xqs[2].w,  xqs[3].x, xqs[3].y, xqs[3].z, xqs[3].w   );
        //for(int i=0; i<ng; i++){  printf("Gs[%i]=%f\n", i, Gs[i]); }
        for(int i=0; i<n; i++){
            
            float3 p = ps[i].xyz;
            //printf( "ps[%3i] ( %8.4f, %8.4f, %8.4f,) \n", i, p.x,p.y,p.z );

            float3 u = (p - g0.xyz) * inv_dg;

            // int ix = (int)u.x; 
            // int iy = (int)u.y;
            // int iz = (int)u.z;
            // int ixyz = iz + ng.z*( iy + ng.y*ix);
            // float4 Es = Eg[ixyz];
            // //printf( "Eg[%3i,%3i,%3i]=(%g,%g,%g,%g) \n", ix,iy,iz, Es.x,Es.y,Es.z,Es.w );
            // float E = dot(Es,C);
            // float4 fe  = (float4){E,E,E,E};

            float4 fe = fe3d_pbc_comb(u, ng.xyz, Eg, C, xqs, yqs);

            fe.xyz *= -inv_dg;

            //printf( "GPU sample3D_comb()[%i] fe(%g,%g,%g | %g) \n",i, fe.x, fe.y, fe.z, fe.w );
            fes[i] = fe;
        }

    }

    /*
    float3 inv_dg = 1.0f / dg.xyz;
    float3 p = ps[iG].xyz;
    float3 u = (p - g0.xyz) * inv_dg;
    float4 fe = fe3d_pbc_comb(u, ng.xyz, Eg, C, xqs, yqs);
    fe.xyz *= inv_dg;
    fes[iG] = fe;
    */
}


inline float2 fe1d_pbc_macro(float x, int n, __global const float* Es, __local const int4* xqis ){
    int i = (int)x;
    if (x < 0) i--;
    float t = x - i;
    i = modulo(i - 1, n);
    int4 q = choose_inds_pbc_3(i, n, xqis);
    //printf( "fe1d_pbc_macro(x=%8.4f) %3i/%3i qi(%i,%i,%i,%i)     q0(%i,%i,%i,%i) q1(%i,%i,%i,%i) q2(%i,%i,%i,%i) q3(%i,%i,%i,%i) \n", x, i, n,  q.x,q.y,q.z,q.w,   xqis[0].x,xqis[0].y,xqis[0].z,xqis[0].w,   xqis[1].x,xqis[1].y,xqis[1].z,xqis[1].w,  xqis[2].x,xqis[2].y,xqis[2].z,xqis[2].w,  xqis[3].x,xqis[3].y,xqis[3].z,xqis[3].w     );
    float4 b = basis(t);
    float4 d = dbasis(t);
    float4 cs = (float4)(Es[q.x], Es[q.y], Es[q.z], Es[q.w]);
    
    return (float2)(dot(b, cs), dot(d, cs));
}

__kernel void sample1D_pbc(
    const float g0,
    const float dg,
    const int ng,
    __global const float* Gs,
    const int n,
    __global const float* ps,
    __global float2* fes
    //__global int4* xqs
) {
    const int iG = get_global_id(0);
    if (iG >= n) return;

    
    __local int4 xqs[4];
    const int iL = get_local_id(0);
    if      (iL<4){ xqs[iL]=make_inds_pbc(ng,iL); }
    barrier(CLK_LOCAL_MEM_FENCE);

    // if( (iG==0) ){
    //     printf("xqs[0](%i,%i,%i,%i)\n xqs[1](%i,%i,%i,%i)\n xqs[2](%i,%i,%i,%i)\n xqs[3](%i,%i,%i,%i)\n", xqs[0].x, xqs[0].y, xqs[0].z, xqs[0].w,   xqs[1].x, xqs[1].y, xqs[1].z, xqs[1].w,   xqs[2].x, xqs[2].y, xqs[2].z, xqs[2].w,  xqs[3].x, xqs[3].y, xqs[3].z, xqs[3].w   );
    //     for(int i=0; i<ng; i++){  printf("Gs[%i]=%f\n", i, Gs[i]); }
    // }

    // local memory barrire
    //int4 xqis[4]; make_inds_pbc(ng, xqis);   // this should be pre-calculated globaly

    float inv_dg = 1.0f / dg;
    float p = ps[iG];
    float2 fe = fe1d_pbc_macro(  (p - g0) * inv_dg, ng, Gs, xqs);
    fe.y *= inv_dg;
    fes[iG] = fe;
}


// =============================================
// ===================== Fitting
// =============================================

float conv3x3_pbc( __global const float* Gs, const float3 B, const int iiz, const int3 ix, const int3 iy ){
    return  Gs[ix.x+iy.x+iiz]*B.z + Gs[ix.y+iy.x+iiz]*B.y + Gs[ix.z+iy.x+iiz]*B.z  +
            Gs[ix.x+iy.y+iiz]*B.y + Gs[ix.y+iy.y+iiz]*B.x + Gs[ix.z+iy.y+iiz]*B.y  +
            Gs[ix.x+iy.z+iiz]*B.z + Gs[ix.y+iy.z+iiz]*B.y + Gs[ix.z+iy.z+iiz]*B.z  ;
}

float conv_3x3_tex( sampler_t samp, __read_only image3d_t tex, float3 B, int4 coord ){
    return
      read_imagef(tex, samp, coord + (int4)(-1,-1,0,0) ).x * B.z
    + read_imagef(tex, samp, coord + (int4)( 0,-1,0,0) ).x * B.y
    + read_imagef(tex, samp, coord + (int4)( 1,-1,0,0) ).x * B.z

    + read_imagef(tex, samp, coord + (int4)(-1, 0,0,0) ).x * B.y
    + read_imagef(tex, samp, coord                     ).x * B.x
    + read_imagef(tex, samp, coord + (int4)( 1, 0,0,0) ).x * B.y

    + read_imagef(tex, samp, coord + (int4)(-1, 1,0,0) ).x * B.z
    + read_imagef(tex, samp, coord + (int4)( 0, 1,0,0) ).x * B.y
    + read_imagef(tex, samp, coord + (int4)( 1, 1,0,0) ).x * B.z;

}

__kernel void BsplineConv3D(
    const int4 ns,
    __global const float* Gs,
    __global const float* G0,
    __global       float* out
) {
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    if (ix >= ns.x || iy >= ns.y || iz >= ns.z) return;


    const float  B0 = 2.0/3.0;
    const float  B1 = 1.0/6.0;
    const float3 Bs = (float3){B0*B0, B0*B1, B1*B1 };
    
    const int3 ixs =  (int3){ modulo(ix-1, ns.x)*ns.x,  modulo(ix, ns.x)*ns.x,   modulo(ix+1, ns.x)  };
    const int3 iys = ((int3){ modulo(iy-1, ns.y)*ns.x,  modulo(iy, ns.y)*ns.x,   modulo(iy+1, ns.y)  })*ns.x;

    float val=0;
    const int                 iiz = modulo(iz, ns.z)*ns.y*ns.x;  val += conv3x3_pbc( Gs, Bs, iiz, ixs, iys ) * B0;
    if(iz>0     ){ const int  biz = modulo(iz, ns.z)*ns.y*ns.x;  val += conv3x3_pbc( Gs, Bs, biz, ixs, iys ) * B1; }
    if(iz<ns.z-1){ const int  diz = modulo(iz, ns.z)*ns.y*ns.x;  val += conv3x3_pbc( Gs, Bs, diz, ixs, iys ) * B1; }

    const int i = iiz + iys.y + ixs.y;

    if (G0 != NULL) { val-=G0[i]; }

    out[i] =  val;
}

// =============================================
// ===================== Fitting  Texture Version
// =============================================

__constant sampler_t samp_pbc = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

__kernel void BsplineConv3D_tex(
    const int4 ns,
    __read_only image3d_t Gs,
    __global const float* G0,
    __global       float* out    
) {

    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    
    if (ix >= ns.x || iy >= ns.y || iz >= ns.z) return;

    const float  B0 = 2.0/3.0;
    const float  B1 = 1.0/6.0;
    const float3 Bs = (float3){B0*B0, B0*B1, B1*B1 };

    int4 coord = (int4){ix, iy, iz, 0};

    float          val  = conv_3x3_tex( samp_pbc, Gs, Bs, coord                  ) * B0;
    if(iz>0     ){ val += conv_3x3_tex( samp_pbc, Gs, Bs, coord-(int4){0,0,0,-1} ) * B1; } 
    if(iz<ns.z-1){ val += conv_3x3_tex( samp_pbc, Gs, Bs, coord-(int4){0,0,0, 1} ) * B1; }

    const int i = ix + ns.x * ( iy* + iz*ns.y );
    
    if (G0 != NULL) { val-=G0[i]; }
    out[i] =  val;

    //out[i] =  E_fit - G0[i];

}

__kernel void move(
    const int  ntot,
    __global float* p,
    __global float* v,
    __global float* f,  
    const float4 par
) {

    const int i = get_global_id(0);
    if (i > ntot ) return;

    // leap frog
    float vi =  v[i];
    float pi =  p[i];
    float fi  = f[i];

    vi *=    par.z;
    vi += fi*par.x;
    pi += vi*par.y;

    v[i]=vi;
    p[i]=pi;
}



