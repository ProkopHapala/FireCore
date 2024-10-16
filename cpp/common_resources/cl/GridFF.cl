
//#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

typedef struct __attribute__ ((packed)){
    float4 a;
    float4 b;
    float4 c;
} cl_Mat3;

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f       // [ eV*Ang/e^2 ]
#define const_kB        8.617333262145e-5f   // [ eV/K ]


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
    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) qx(%i,%i,%i,%i) nyz=%i\n", u.x,u.y,u.z, ix,iy,iz, qx.x,qx.y,qx.z,qx.w, nyz );
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
    const float3 inv_dg = 1.0f / dg.xyz;
    barrier(CLK_LOCAL_MEM_FENCE);

    // if( iG==0 ){
    //     printf("GPU sample3D_comb() ng(%i,%i,%i) g0(%g,%g,%g) dg(%g,%g,%g) C(%g,%g,%g) \n", ng.x,ng.y,ng.z,   g0.x,g0.y,g0.z,   dg.x,dg.y,dg.z,   C.x,C.y,C.z );
    //     printf("GPU xqs[0](%i,%i,%i,%i) xqs[1](%i,%i,%i,%i) xqs[2](%i,%i,%i,%i) xqs[3](%i,%i,%i,%i)\n", xqs[0].x, xqs[0].y, xqs[0].z, xqs[0].w,   xqs[1].x, xqs[1].y, xqs[1].z, xqs[1].w,   xqs[2].x, xqs[2].y, xqs[2].z, xqs[2].w,  xqs[3].x, xqs[3].y, xqs[3].z, xqs[3].w   );
    //     //for(int i=0; i<ng; i++){  printf("Gs[%i]=%f\n", i, Gs[i]); }
    //     for(int i=0; i<n; i++){
    //         float3 p = ps[i].xyz;
    //         //printf( "ps[%3i] ( %8.4f, %8.4f, %8.4f,) \n", i, p.x,p.y,p.z );
    //         float3 u = (p - g0.xyz) * inv_dg;
    //         // int ix = (int)u.x; 
    //         // int iy = (int)u.y;
    //         // int iz = (int)u.z;
    //         // int ixyz = iz + ng.z*( iy + ng.y*ix);
    //         // float4 Es = Eg[ixyz];
    //         // //printf( "Eg[%3i,%3i,%3i]=(%g,%g,%g,%g) \n", ix,iy,iz, Es.x,Es.y,Es.z,Es.w );
    //         // float E = dot(Es,C);
    //         // float4 fe  = (float4){E,E,E,E};
    //         float4 fe = fe3d_pbc_comb(u, ng.xyz, Eg, C, xqs, yqs);
    //         fe.xyz *= -inv_dg;
    //         //printf( "GPU sample3D_comb()[%i] fe(%g,%g,%g | %g) \n",i, fe.x, fe.y, fe.z, fe.w );
    //         fes[i] = fe;
    //     }
    // }

    float3 p = ps[iG].xyz;
    float3 u = (p - g0.xyz) * inv_dg;
    float4 fe = fe3d_pbc_comb(u, ng.xyz, Eg, C, xqs, yqs);
    fe.xyz *= -inv_dg;
    fes[iG] = fe;
    
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
    __global       float* out,
    const float2 coefs
) {
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);

    //if( (ix==0)&&(iy==0)&&(iz==0) ){ printf("GPU BsplineConv3D() ns{%i,%i,%i,%i}\n", ns.x,ns.y,ns.z,ns.w); }
    if( (ix>=ns.x) || (iy>=ns.y) || (iz>=ns.z) ) return;

    const float  B0 = 2.0/3.0;
    const float  B1 = 1.0/6.0;
    const float3 Bs = (float3){B0*B0, B0*B1, B1*B1 };
    
    const int3 ixs =  (int3){ modulo(ix-1,ns.x),  ix,   modulo(ix+1,ns.x)  };
    const int3 iys = ((int3){ modulo(iy-1,ns.y),  iy,   modulo(iy+1,ns.y)  })*ns.x;

    const int nxy = ns.x*ns.y;

    float val=0;
    const int iiz =iz*nxy;  val += conv3x3_pbc( Gs, Bs, iiz                    , ixs, iys ) * B0;
    if(iz>0     ){          val += conv3x3_pbc( Gs, Bs, modulo(iz-1, ns.z)*nxy , ixs, iys ) * B1; }
    if(iz<ns.z-1){          val += conv3x3_pbc( Gs, Bs, modulo(iz+1, ns.z)*nxy , ixs, iys ) * B1; }
    
    const int i = iiz + iys.y + ixs.y;
    val*=coefs.x;
    if (G0 != NULL) { val+=G0[i]*coefs.y; }
    out[i] =  val;

    // const int i = ix + ns.x*( iy + ns.y*iz);
    // // out[i] =  Gs[i];
    // // out[i] =  G0[i];
    // out[i] =  G0[i] - Gs[i];


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
    
    //if( (ix==0)&&(iy==0)&&(iz==0) ){ printf("GPU BsplineConv3D_tex() ns{%i,%i,%i,%i}\n", ns.x,ns.y,ns.z,ns.w); }
    if( (ix>=ns.x) || (iy>=ns.y) || (iz>=ns.z) ) return;

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

}

__kernel void move(
    const int  ntot,
    __global float* p,
    __global float* v,
    __global float* f,  
    const float4 MDpar
) {

    const int i = get_global_id(0);
    //if( i==0 ){ printf("GPU move() ntot=%i MDpar{%g,%g,%g,%g}\n", ntot,  MDpar.x, MDpar.y, MDpar.z,MDpar.w); }
    if (i > ntot ) return;

    // leap frog
    float vi =  v[i];
    float pi =  p[i];
    float fi  = f[i];

    vi *=    MDpar.z;
    vi += fi*MDpar.x;
    pi += vi*MDpar.y;

    v[i]=vi;
    p[i]=pi;
}

__kernel void setMul(
    const int  ntot,
    __global float* v,
    __global float* out,  
    float c
) {
    const int i = get_global_id(0);
    //if( i==0 ){ printf("GPU move() ntot=%i MDpar{%g,%g,%g,%g}\n", ntot,  MDpar.x, MDpar.y, MDpar.z,MDpar.w); }
    if (i > ntot ) return;
    out[i] = v[i]*c;
}


__attribute__((reqd_work_group_size(32,1,1)))
__kernel void make_MorseFF(
    const int nAtoms,                // 1
    __global float4*  atoms,         // 2
    __global float4*  REQs,          // 3
    __global float* FE_Paul,         // 4
    __global float* FE_Lond,         // 5
    //__global * FE_Coul,
    const int4     nPBC,             // 6
    const int4     nGrid,            // 7
    //const cl_Mat3  lvec,           
    const float4  lvec_a,            // 8
    const float4  lvec_b,            // 9
    const float4  lvec_c,            // 10
    const float4  grid_p0,           // 11
    const float4  GFFParams          // 12
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  =  iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  =  iG/nab; 

    const float  alphaMorse = GFFParams.y;
    const float  R2damp     = GFFParams.x*GFFParams.x;
    const float3 dGrid_a = lvec_a.xyz*(1.f/(float)nGrid.x);
    const float3 dGrid_b = lvec_b.xyz*(1.f/(float)nGrid.y);
    const float3 dGrid_c = lvec_c.xyz*(1.f/(float)nGrid.z); 
    const float3 shift_b = lvec_b.xyz + lvec_a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec_c.xyz + lvec_b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz) 
    
    //if( (ia==0)&&(ib==0) ){  printf(  "GPU ic %i nGrid(%i,%i,%i)\n", ic, nGrid.x,nGrid.y,nGrid.z );}

    const int nMax = nab*nGrid.z;
    if(iG>=nMax) return;

    const float3 pos    = grid_p0.xyz  + dGrid_a.xyz*ia      + dGrid_b.xyz*ib      + dGrid_c.xyz*ic       // grid point within cell
                                       +  lvec_a.xyz*-nPBC.x + lvec_b.xyz*-nPBC.y + lvec_c.xyz*-nPBC.z;  // most negative PBC-cell

    //const float3  shift0 = lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
    float fe_Paul = 0.0f;
    float fe_Lond = 0.0f;
    //float4 fe_Coul = float4Zero;
    for (int j0=0; j0<nAtoms; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = REQs [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( ja<nAtoms ){ 
                const float4 REQK =       LCLJS [jl];
                float3       dp   = pos - LATOMS[jl].xyz;
            
                //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                //dp+=lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;

                //float3 shift=shift0; 
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){

                            //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc[%i,%i,%i] dp(%g,%g,%g)\n", ix,iy,iz, dp.x,dp.y,dp.z );   
                            float  r2  = dot(dp,dp);
                            float  r   = sqrt(r2+1e-32 );
                            // ---- Electrostatic
                            //float ir2  = 1.f/(r2+R2damp); 
                            //float   E  = COULOMB_CONST*REQK.z*sqrt(ir2);
                            //fe_Coul   += (float4)(dp*(E*ir2), E );
                            // ---- Morse ( Pauli + Dispersion )
                            float    e = exp( -alphaMorse*(r-REQK.x) );
                            float   eM = REQK.y*e;
                            fe_Paul += eM * e;
                            fe_Lond += eM * -2.0f;

                            // if((iG==0)&&(j==0)){
                            //     //float3 sh = dp - pos + LCLJS[j].xyz + lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
                            //     float3 sh = shift;
                            //     printf( "GPU(%2i,%2i,%2i) sh(%7.3f,%7.3f,%7.3f)\n", ix,iy,iz, sh.x,sh.y,sh.z  );
                            // }
                            //ipbc++; 
                            
                            dp   +=lvec_a.xyz;
                            //shift+=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                        //shift+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp   +=shift_c;
                    //shift+=shift_c;
                    //dp+=lvec.b.xyz*(nPBC.y*-2.f-1.f);
                    //dp+=lvec.c.xyz;
                }

            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE_Paul[iG] = fe_Paul;
    FE_Lond[iG] = fe_Lond;
    //FE_Coul[iG] = fe_Coul;
    //int4 coord = (int4){ia,ib,ic,0};
    //write_imagef( FE_Paul, coord, (float4){pos,(float)iG} );
    //write_imagef( FE_Paul, coord, fe_Paul );
    //write_imagef( FE_Lond, coord, fe_Lond );
    //write_imagef( FE_Coul, coord, fe_Coul );
}


__attribute__((reqd_work_group_size(32,1,1)))
__kernel void make_MorseFF_f4(
    const int nAtoms,                // 1
    __global float4*  atoms,         // 2
    __global float4*  REQs,          // 3
    __global float4* FE_Paul,
    __global float4* FE_Lond,
    // __global float4* FE_Coul,
    const int4     nPBC,             // 7
    const int4     nGrid,            // 8
    const cl_Mat3  lvec,             // 9
    const float4   grid_p0,          // 10
    const float4   GFFParams         // 11
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  =  iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  =  iG/nab; 

    const float  alphaMorse = GFFParams.y;
    const float  R2damp     = GFFParams.x*GFFParams.x;
    const float3 dGrid_a = lvec.a.xyz*(1.f/(float)nGrid.x);
    const float3 dGrid_b = lvec.b.xyz*(1.f/(float)nGrid.y);
    const float3 dGrid_c = lvec.c.xyz*(1.f/(float)nGrid.z); 
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz) 
    
    const int nMax = nab*nGrid.z;
    if(iG>=nMax) return;

    const float3 pos    = grid_p0.xyz  + dGrid_a.xyz*ia      + dGrid_b.xyz*ib      + dGrid_c.xyz*ic       // grid point within cell
                                       +  lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell

    //const float3  shift0 = lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
    float4 fe_Paul = float4Zero;
    float4 fe_Lond = float4Zero;
    //float4 fe_Coul = float4Zero;
    for (int j0=0; j0<nAtoms; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = REQs [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( ja<nAtoms ){ 
                const float4 REQK =       LCLJS [jl];
                float3       dp   = pos - LATOMS[jl].xyz;
            
                //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                //dp+=lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;

                //float3 shift=shift0; 
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){

                            //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc[%i,%i,%i] dp(%g,%g,%g)\n", ix,iy,iz, dp.x,dp.y,dp.z );   
                            float  r2  = dot(dp,dp);
                            float  r   = sqrt(r2+1e-32 );
                            // ---- Electrostatic
                            //float ir2  = 1.f/(r2+R2damp); 
                            //float   E  = COULOMB_CONST*REQK.z*sqrt(ir2);
                            //fe_Coul   += (float4)(dp*(E*ir2), E );
                            // ---- Morse ( Pauli + Dispersion )
                            float    e = exp( -alphaMorse*(r-REQK.x) );
                            float   eM = REQK.y*e;
                            float   de = 2.f*alphaMorse*eM/r;
                            float4  fe = (float4)( dp*de, eM );
                            fe_Paul += fe * e;
                            fe_Lond += fe * (float4)( -1.0f,-1.0f,-1.0f, -2.0f );

                            // if((iG==0)&&(j==0)){
                            //     //float3 sh = dp - pos + LCLJS[j].xyz + lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
                            //     float3 sh = shift;
                            //     printf( "GPU(%2i,%2i,%2i) sh(%7.3f,%7.3f,%7.3f)\n", ix,iy,iz, sh.x,sh.y,sh.z  );
                            // }
                            //ipbc++; 
                            
                            dp   +=lvec.a.xyz;
                            //shift+=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                        //shift+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp   +=shift_c;
                    //shift+=shift_c;
                    //dp+=lvec.b.xyz*(nPBC.y*-2.f-1.f);
                    //dp+=lvec.c.xyz;
                }

            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    FE_Paul[iG] = fe_Paul;
    FE_Lond[iG] = fe_Lond;
    //FE_Coul[iG] = fe_Coul;

    //int4 coord = (int4){ia,ib,ic,0};
    //write_imagef( FE_Paul, coord, (float4){pos,(float)iG} );
    //write_imagef( FE_Paul, coord, fe_Paul );
    //write_imagef( FE_Lond, coord, fe_Lond );
    //write_imagef( FE_Coul, coord, fe_Coul );
}


int pbc_ifw(int i, int n){ i++; return (i<n )?  i :  i-n; };
int pbc_ibk(int i, int n){ i--; return (i>=0)?  i :  i+n; };


// float4 Bspline_basis(const float u) {
//     const float inv6 = 1.0f / 6.0f;
//     const float u2 = u * u;
//     const float t = 1.0f - u;
//     return (float4)(
//         inv6 * t * t * t,
//         inv6 * (3.0f * u2 * (u - 2.0f) + 4.0f),
//         inv6 * (3.0f * u * (1.0f + u - u2) + 1.0f),
//         inv6 * u2 * u
//     );
// }
// float4 Bspline_dbasis(const float u) {
//     const float u2 = u * u;
//     const float t = 1.0f - u;
//     return (float4)(
//         -0.5f * t * t,
//         0.5f * (3.0f * u2 - 4.0f * u),
//         0.5f * (-3.0f * u2 + 2.0f * u + 1.0f),
//         0.5f * u2
//     );
// }

void Bspline_basis(const float u, float * ws) {
    const float inv6 = 1.0f / 6.0f;
    const float u2 = u * u;
    const float t = 1.0f - u;
    //return (float4)(
    ws[0]=    inv6 * t * t * t;
    ws[1]=    inv6 * (3.0f * u2 * (u - 2.0f) + 4.0f);
    ws[2]=    inv6 * (3.0f * u * (1.0f + u - u2) + 1.0f);
    ws[3]=    inv6 * u2 * u;
    //);
}

void Bspline_dbasis(const float u, float * ws) {
    const float u2 = u * u;
    const float t = 1.0f - u;
    //return (float4)(
    ws[0]=    -0.5f * t * t;
    ws[1]=     0.5f * ( 3.0f * u2 - 4.0f * u);
    ws[2]=     0.5f * (-3.0f * u2 + 2.0f * u + 1.0f);
    ws[3]=     0.5f * u2;
    //);
}


void Bspline_basis5(const float t, float * ws){
    const float inv6 = 1.f/6.f;
    const float t2 = t*t;
    const float t3 = t2*t;
    const float t4 = t2*t2;
    const float t5 = t3*t2;
    //return (float8){                                                  
    ws[0]=  -0.008333333333333333*t5  +0.041666666666666666*t4  -0.08333333333333333*t3 +0.08333333333333333*t2  -0.041666666666666666*t   +0.008333333333333333;
    ws[1]=   0.041666666666666666*t5  -0.166666666666666666*t4  +0.16666666666666666*t3 +0.16666666666666666*t2  -0.416666666666666666*t   +0.216666666666666666;        
    ws[2]=  -0.083333333333333333*t5  +0.250000000000000000*t4                          -0.50000000000000000*t2                            +0.550000000000000000;  
    ws[3]=   0.083333333333333333*t5  -0.166666666666666666*t4  -0.16666666666666666*t3 +0.16666666666666666*t2  +0.416666666666666666*t   +0.216666666666666666;
    ws[4]=  -0.041666666666666666*t5  +0.041666666666666666*t4  +0.08333333333333333*t3 +0.08333333333333333*t2  +0.041666666666666666*t   +0.008333333333333333; 
    ws[5]=   0.008333333333333333*t5;
    //     0.f,0.f,
    //};
}


void Bspline_dbasis5(const float t, float * ws){
    const float inv6 = 1.f/6.f;
    const float t2 = t*t;
    const float t3 = t2*t;
    const float t4 = t2*t2;
    //return (float8){           
    ws[0]=    -0.0416666666666667*t4	+0.166666666666667*t3	-0.25*t2   +0.166666666666667*t	-0.041666666666666666;	
    ws[1]=     0.2083333333333333*t4	-0.666666666666667*t3	+0.50*t2   +0.333333333333333*t -0.416666666666666666;	
    ws[2]=    -0.4166666666666667*t4	+1.000000000000000*t3	           -1.000000000000000*t                      ;	
    ws[3]=     0.4166666666666667*t4	-0.666666666666667*t3	-0.50*t2   +0.333333333333333*t	+0.416666666666666666;	
    ws[4]=    -0.2083333333333333*t4	+0.166666666666667*t3	+0.25*t2   +0.166666666666667*t	+0.041666666666666666;	
    ws[5]=     0.0416666666666667*t4;
    //};
}


__kernel void project_atom_on_grid_cubic_pbc(
    const int na,                   // 1 number of atoms
    __global const float4* atoms,   // 2 Atom positions and charges
    __global       float*  Qgrid,   // 3 Output grid
    const int4 ng,                  // 4 grid size
    const float3 g0,                // 5 grid orgin
    const float3 dg                 // 6 grid dimensions
) {
    int iG = get_global_id(0);
    const int iL = get_local_id(0);
    if (iG >= na) return;

    __local int4 xqs[4];
    __local int4 yqs[4];
    __local int4 zqs[4];
    if      (iL<4 ){             xqs[iL]=make_inds_pbc(ng.x,iL); }
    else if (iL<8 ){ int i=iL-4; yqs[i ]=make_inds_pbc(ng.y,i ); }
    else if (iL<12){ int i=iL-8; yqs[i ]=make_inds_pbc(ng.y,i ); };
    barrier(CLK_LOCAL_MEM_FENCE);


    // Load atom position and charge
    float4 atom = atoms[iG];
    //float3 pos  = (float3)(atom_data.x, atom_data.y, atom_data.z);
    //float charge = atom_data.w;

    // Convert to grid coordinates
    float3      g = (atom.xyz - g0) / dg;
    int3       gi = (int3  ){(int)g.x, (int)g.y, (int)g.z};
    float3 t      = (float3){     g.x - gi.x, g.y - gi.y, g.z - gi.z};

    // Compute weights for cubic B-spline interpolation
    float wx[4], wy[4], wz[4];
    Bspline_basis(t.x, wx);
    Bspline_basis(t.y, wy);
    Bspline_basis(t.z, wz);

    const int nxy = ng.x * ng.y;
    // Pre-calculate periodic boundary condition indices for each dimension
    gi.x=modulo(gi.x-1,ng.x); const int4 xq = choose_inds_pbc_3(gi.x, ng.x, xqs );  const int* xq_ = (int*)&xq;
    gi.y=modulo(gi.y-1,ng.y); const int4 yq = choose_inds_pbc_3(gi.y, ng.y, yqs );  const int* yq_ = (int*)&xq;
    gi.z=modulo(gi.z-1,ng.z); const int4 zq = choose_inds_pbc_3(gi.z, ng.z, zqs );  const int* zq_ = (int*)&xq;

    //float4 Bspline_dbasis();

    for (int dz = 0; dz < 4; dz++) {
        const int gz  = zq_[dz];
        const int iiz = gz * nxy;
        for (int dy = 0; dy < 4; dy++) {
            const int gy = yq_[dy];
            const int iiy = iiz + gy * ng.x;
            const double qbyz = atom.w * wy[dy] * wz[dz];
            for (int dx = 0; dx < 4; dx++) {
                const int gx = xq_[dx];
                const int ig = gx + iiy;
                double qi = qbyz * wx[dx];
                Qgrid[ig] += qi;
            }
        }
    }

}

inline void make_inds_pbc_5(const int n, const int iG, __local int inds[6]) {
    switch (iG) {
        case 0:  inds[0]=0;    inds[1]=1;    inds[2]=2;    inds[3]=3;    inds[4]=4;    inds[5]=5;    break;
        case 1:  inds[0]=0;    inds[1]=1;    inds[2]=2;    inds[3]=3;    inds[4]=4;    inds[5]=5-n;  break;
        case 2:  inds[0]=0;    inds[1]=1;    inds[2]=2;    inds[3]=3;    inds[4]=4-n;  inds[5]=5-n;  break;
        case 3:  inds[0]=0;    inds[1]=1;    inds[2]=2;    inds[3]=3-n;  inds[4]=4-n;  inds[5]=5-n;  break;
        case 4:  inds[0]=0;    inds[1]=1;    inds[2]=2-n;  inds[3]=3-n;  inds[4]=4-n;  inds[5]=5-n;  break;
        case 5:  inds[0]=0;    inds[1]=1-n;  inds[2]=2-n;  inds[3]=3-n;  inds[4]=4-n;  inds[5]=5-n;  break;
        default: inds[0]=-100; inds[1]=-100; inds[2]=-100; inds[3]=-100; inds[4]=-100; inds[5]=-100; break;
    }
}

inline void choose_inds_pbc_5(const int i, const int n, __local const int iqs[6][6], int out[6]) {
    if (i >= (n - 5)) {
        const int ii  = i+6-n;
        const int* qi = iqs[ii];
              out[0]=i+qi[0];    out[1]=i+qi[1];    out[2]=i+qi[2];    out[3]=i+qi[3];    out[4]=i+qi[4];    out[5]=i+qi[5];
    } else {  out[0]=i;          out[1]=i+1;        out[2]=i+2;        out[3]=i+3;        out[4]=i+4;        out[5]=i+5;      }
}


__kernel void project_atoms_on_grid_quintic_pbc(
    const int na,                   // 1 number of atoms
    __global const float4* atoms,   // 2 Atom positions and charges
    __global       float2* Qgrid,   // 3 Output grid (complex, in order to be compatible with poisson)
    const int4   ng,                // 4 Grid size
    const float4 g0,                // 5 Grid origin
    const float4 dg                 // 6 Grid dimensions
) {
    int iG = get_global_id(0);
    const int iL = get_local_id(0);
    
    // Declare and initialize shared memory for periodic boundary condition indices
    __local int xqs[6][6];
    __local int yqs[6][6];
    __local int zqs[6][6];
    if      (iL<6 ) { const int i=iL;    make_inds_pbc_5(ng.x,i,xqs[i]); }
    else if (iL<12) { const int i=iL-6;  make_inds_pbc_5(ng.y,i,yqs[i]); }
    else if (iL<18) { const int i=iL-12; make_inds_pbc_5(ng.z,i,zqs[i]); }
    barrier(CLK_LOCAL_MEM_FENCE);
    if (iG >= na) return;

    // if( iG==0 ){
    //     printf("GPU project_atoms_on_grid_quintic_pbc() ng(%i,%i,%i) g0(%g,%g,%g) dg(%g,%g,%g) \n", ng.x,ng.y,ng.z,   g0.x,g0.y,g0.z,   dg.x,dg.y,dg.z );
    //     for(int i=0; i<6; i++){ int* q=xqs[i]; printf("GPU xqs[0](%4i,%4i,%4i,%4i,%4i,%4i) \n", q[0],  q[1], q[2], q[3], q[4], q[5] ); }
    //     for(int i=0; i<6; i++){ int* q=yqs[i]; printf("GPU yqs[0](%4i,%4i,%4i,%4i,%4i,%4i) \n", q[0],  q[1], q[2], q[3], q[4], q[5] ); }
    //     for(int i=0; i<6; i++){ int* q=zqs[i]; printf("GPU zqs[0](%4i,%4i,%4i,%4i,%4i,%4i) \n", q[0],  q[1], q[2], q[3], q[4], q[5] ); }
    //     for(int ia=0; ia<na; ia++){ printf("GPU atom[%i](%8.4f,%8.4f,%8.4f |%8.4f) \n", ia, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w ); }
    // }

    // Load atom position and charge
    float4 atom = atoms[iG];
    float3 g    = (atom.xyz - g0.xyz) / dg.xyz;
    int3   gi   = (int3  ){(int)g.x,(int)g.y,(int)g.z};
    float3 t    = (float3){g.x-gi.x, g.y-gi.y, g.z-gi.z};

    // Compute weights for quintic B-spline interpolation
    float bx[6], by[6], bz[6];
    Bspline_basis5(t.x, bx);
    Bspline_basis5(t.y, by);
    Bspline_basis5(t.z, bz);

    const int nxy = ng.x * ng.y;
    
    int xq[6];
    int yq[6];
    int zq[6];
    // Pre-calculate periodic boundary condition indices for each dimension
    gi.x = modulo( gi.x-2, ng.x ); choose_inds_pbc_5(gi.x,ng.x, xqs, xq );
    gi.y = modulo( gi.y-2, ng.y ); choose_inds_pbc_5(gi.y,ng.y, yqs, yq );
    gi.z = modulo( gi.z-2, ng.z ); choose_inds_pbc_5(gi.z,ng.z, zqs, zq );

    for (int dz = 0; dz < 6; dz++) {
        const int gz  = zq[dz];
        const int iiz = gz * nxy;
        for (int dy = 0; dy < 6; dy++) {
            const int gy  = yq[dy];
            const int iiy = iiz + gy * ng.x;
            const float qbyz = atom.w * by[dy] * bz[dz];
            for (int dx = 0; dx < 6; dx++) {
                const int gx = xq[dx];
                const int ig = gx + iiy;
                float qi = qbyz * bx[dx];
                //Qgrid[ig].x += qi;
                Qgrid[ig] = (float2){qi,0.0f};
            }
        }
    }
    //const int ig = gi.z*nxy + gi.y*ng.x + gi.x;
    //Qgrid[ig] = (float2){gi.y*1.0f,0.0f};
}

__kernel void poissonW(
    const int4   ns,         // (nx,ny,nz,nxyz)
    __global float2* rho_k,  // input array  rho(k) - fourier coefficients (complex)
    __global float2* V_k,    // output array V(k)   - fourier coefficients (complex)
    const float4 coefs       // (0,0,0, 4*pi*eps0*dV)
){    
    const int iG = get_global_id (0);
    if(iG>=ns.w) return;
    const int nab = ns.x*ns.y;
    const int ix  =  iG%ns.x; 
    const int iy  = (iG%nab)/ns.x;
    const int iz  =  iG/nab; 
    float4 k = (float4){ ix/(0.5f*ns.x), iy/(0.5f*ns.y), iz/(0.5f*ns.z), 0};
    k = 1.0f-fabs(k-1.0f); 
    float  f = coefs.w/dot( k, k );    // dCell.w = 4*pi*eps0*dV - rescaling constant
    if(iG==0)f=0;
    if(iG<ns.w){ 
        V_k[iG] = rho_k[iG]*f;
    }
};


__kernel void laplace_real_pbc( 
    int4 ng,
    __global float* Vin, 
    __global float* Vout, 
    __global float* vV, 
    float cSOR, 
    float cV
){
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    if( (ix>=ng.x) || (iy>=ng.y) || (iz>=ng.z) ) return;

    int nxy = ng.x * ng.y;
    const float fac = 1.0f/6.0f;
    
    const int iiz =          iz       *nxy;
    const int ifz =  pbc_ifw(iz, ng.z)*nxy;
    const int ibz =  pbc_ibk(iz, ng.z)*nxy;
    
    const int iiy =          iy       *ng.x;
    const int ify =  pbc_ifw(iy, ng.y)*ng.x;
    const int iby =  pbc_ibk(iy, ng.y)*ng.x;
    const int ifx =  pbc_ifw(ix, ng.x);
    const int ibx =  pbc_ibk(ix, ng.x);

    float vi = 
    Vin[ ibx + iiy + iiz ] + Vin[ ifx + iiy + iiz ] + 
    Vin[ ix  + iby + iiz ] + Vin[ ix  + ify + iiz ] + 
    Vin[ ix  + iiy + ibz ] + Vin[ ix  + iiy + ifz ];

    vi *= fac;
    const int i = ix + iiy + iiz;
    const float vo = Vin[ i ];
    vi += (vi-vo)*cSOR; 
    
    if(vV != 0){   // inertia
        float v = vi - vo;                 // velocity ( change between new and old potential )
        v       = v*cV + vV[i]*(1.0f-cV);  // inertia ( mixing of new and old change )
        vV[i]   = v;                       // store updated velocity ( change )
        vi      = v + vo;                  // new potantial corrected by intertia
    }

    Vout[i] = vi;

    // double v = V_[i]-V[i];
    // if(iter>0){ v = v*cV + vV[i]*(1-cV); }
    // vV[i] = v; 
    // V_[i] = V[i] + v;

}

__kernel void slabPotential( 
    int4 ns,
    __global float2* Vin,   // real and imag part ( result of Particle Mesh Ewald poissin solever using clFFT
    __global float*  Vout,
    float4 params,          // (dz, Vol, dVcor, Vcor0)
    int4 ng
){
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    if( (ix>=ng.x) || (iy>=ng.y) || (iz>=ns.w) ) return;

    float dz    = params.x;
    float dVcor = params.z;
    float Vcor0 = params.w;

    float Vcor_z = Vcor0 + dVcor * (iz*dz);
    const int i = ix + ng.x*(iy + ng.y*iz);
    Vout[i] = Vin[i].x + Vcor_z;
}

