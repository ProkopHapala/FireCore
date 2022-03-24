
#define R2SAFE          1e-4f
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

//#define N_RELAX_STEP_MAX  64
#define N_RELAX_STEP_MAX  16
#define F2CONV  1e-8

#include "splines.cl"

inline float3 rotMat ( float3 v,  float3 a, float3 b, float3 c  ){ return (float3)(dot(v,a),dot(v,b),dot(v,c)); }
inline float3 rotMatT( float3 v,  float3 a, float3 b, float3 c  ){ return a*v.x + b*v.y + c*v.z; }

float4 getCoulomb( float4 atom, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     float   ir2 = 1.0f/( dot(dp,dp) +  R2SAFE );
     float   ir  = sqrt(ir2);
     float   E   = atom.w*sqrt(ir2);
     return (float4)(dp*(E*ir2), E );
}

float4 getLJ( float3 apos, float2 cLJ, float3 pos ){
     float3  dp  =  pos - apos;
     float   ir2 = 1.0f/( dot(dp,dp) + R2SAFE );
     float   ir6 = ir2*ir2*ir2;
     float   E   =  (    cLJ.y*ir6 -   cLJ.x )*ir6;
     float3  F   = (( 12.0f*cLJ.y*ir6 - 6.0f*cLJ.x )*ir6*ir2)*dp;
     return (float4)(F, E);
}

float4 getMorse( float3 dp, float3 REA ){
    //float3  dp  =  pos - apos;
    float   r     = sqrt( dot(dp,dp) + R2SAFE );
    float   expar = exp( REA.z*(r-REA.x) );
    float   E     = REA.y*expar*( expar - 2 );
    float   fr    = REA.y*expar*( expar - 1 )*2*REA.z;
    return (float4)(dp*(fr/r), E);
}

float8 getLJC( float4 atom, float2 cLJ, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     float   ir2 = 1.0/( dot(dp,dp) +  R2SAFE );
     float   ir6 = ir2*ir2*ir2;
     float   ELJ =  (    cLJ.y*ir6 -   cLJ.x )*ir6;
     float3  FLJ = (( 12.0f*cLJ.y*ir6 - 6.0f*cLJ.x )*ir6*ir2)*dp;
     float   ir  = sqrt(ir2);
     float   Eel = atom.w*sqrt(ir2);
     return (float8)(FLJ, ELJ, dp*(Eel*ir2), Eel );
}

float getLorenz( float4 atom, float4 coefs, float3 pos ){
     float3  dp  =  pos - atom.xyz;
     return coefs.x/( dot(dp,dp) +  coefs.y*coefs.y );
     //return 1.0/( dot(dp,dp) +  0.000 );
}


__kernel void evalCoulomb(
    int nAtoms, 
    __global float4* atoms, 
    __global float4*  poss,
    __global float4*    FE
){
    __local float4 LATOMS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        if(i>=nAtoms) break;
        //if(iL==0) printf("%i (%f,%f,%f)  %f \n", i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].w );
        
        LATOMS[iL] = atoms[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            fe += getCoulomb( LATOMS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;
    //FE[iG] = poss[iG];
}

__kernel void evalLJ(
    const int nAtoms, 
    __global float4* atoms,
    __global float2*  cLJs,
    __global float4*  poss,
    __global float4*    FE
){
    __local float4 LATOMS[32];
    __local float2 LCLJS  [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        if(i>=nAtoms) break;
        //if(iL==0) printf("%i (%f,%f,%f)  %f \n", i, atoms[i].x, atoms[i].y, atoms[i].z, atoms[i].w );
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            fe += getLJ( LATOMS[j].xyz, LCLJS [j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;
    //FE[iG] = poss[iG];
}

__kernel void evalLJC_Q(
    const int nAtoms, 
    __global float4* atoms,
    __global float2*  cLJs,
    __global float4*  poss,
    __global float4*    FE,
    float Qmix
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    //if(iG==0){ printf("evalLJC: nAtoms: %i \n", nAtoms ); }
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    //fe.hi  = fe.hi*COULOMB_CONST;
    Qmix *= COULOMB_CONST;
    FE[iG] = fe.lo + Qmix * fe.hi;
    //FE[iG] = poss[iG];
}

__kernel void evalLJC_Q_noPos(
    const int nAtoms, 
    __global float4* atoms,
    __global float2*  cLJs,
    __global float4*    FE,
    int4 nGrid,
    float4 grid_p0,
    float4 grid_dA,
    float4 grid_dB,
    float4 grid_dC,
    float Qmix
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    //float3 pos = poss[iG].xyz;
    //float3 pos = grid_p0 + grid_dA*get_global_id(0) + grid_dA*get_global_id(1)  + grid_dA*get_global_id (2);      // there would be more problematic local_id optimization

    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;

    if(iG>nMax) return;

    float3 pos    = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;

    // if(iG==430662){ 
    //     printf("evalLJC_Q_noPos: nAtom %i nGrid(%i,%i,%i,%i)   nL %i  nG %i nG_ %i  \n", nAtoms, nGrid.x, nGrid.y, nGrid.z, nGrid.w, nL, get_global_size(0), nGrid.x*nGrid.y*nGrid.z ); 
    //     printf("evalLJC_Q_noPos: grid.p0 (%g,%g,%g) \n", grid_p0.x, grid_p0.y, grid_p0.z ); 
    //     printf("evalLJC_Q_noPos: grid.p0 (%g,%g,%g) \n", grid_dA.x, grid_dA.y, grid_dA.z ); 
    //     printf("evalLJC_Q_noPos: grid.p0 (%g,%g,%g) \n", grid_dB.x, grid_dB.y, grid_dB.z ); 
    //     printf("evalLJC_Q_noPos: grid.p0 (%g,%g,%g) \n", grid_dC.x, grid_dC.y, grid_dC.z ); 
    //     // 150*150*19 + 150*21 + 12 = 430662
    //     printf("evalLJC_Q_noPos: iG %i ia,ib,ic(%i,%i,%i)   pos(%g,%g,%g) \n", iG, ia, ib, ic,   pos.x,pos.y,pos.z  ); 
    // }

    float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    //fe.hi  = fe.hi*COULOMB_CONST;

    //if(iG==430662){      printf("evalLJC_Q_noPos: iG %i ia,ib,ic(%i,%i,%i)   pos(%g,%g,%g) Qmix %g fe.lo(%g,%g,%g)  fe.hi(%g,%g,%g) \n", iG, ia, ib, ic,   pos.x,pos.y,pos.z,   Qmix,   fe.lo.x,fe.lo.y,fe.lo.z,   fe.hi.x,fe.hi.y,fe.hi.z  ); }

    Qmix *= COULOMB_CONST;
    FE[iG] = fe.lo + Qmix * fe.hi;
    //FE[iG] = poss[iG];
}

__kernel void evalLJC_QZs_noPos(
    const int nAtoms, 
    __global float4* atoms,
    __global float2*  cLJs,
    __global float4*    FE,
    int4 nGrid,
    float4 grid_p0,
    float4 grid_dA,
    float4 grid_dB,
    float4 grid_dC,
    float4 Qs,
    float4 QZs
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;

    if(iG>nMax) return;
    //if(iG==0) printf( " Qs (%g,%g,%g,%g) QZs (%g,%g,%g,%g) \n", Qs.x,Qs.y,Qs.z,Qs.w,   QZs.x,QZs.y,QZs.z,QZs.w   );

    float3 pos    = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;

    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    
    Qs *= COULOMB_CONST;

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                //fe += getLJC( LATOMS[j], LCLJS[j], pos );
                float4 xyzq = LATOMS[j];
                fe += getLJ     ( xyzq.xyz, LCLJS[j], pos );
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.x) ) * Qs.x;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.y) ) * Qs.y;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.z) ) * Qs.z;
                fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.w) ) * Qs.w;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    //if ( (ia==75)&&(ib==75) ) { printf(" iz %i fe %g,%g,%g,%g \n", ic, fe.x, fe.y, fe.z, fe.w ); }

    FE[iG] = fe;
}




__kernel void evalLJC(
    int nAtoms, 
    __global float4*   atoms,
    __global float2*    cLJs,
    __global float4*    poss,
    __global float8*    FE
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    //if(iG==0){ printf("evalLJC: nAtoms: %i \n", nAtoms ); }
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    fe.hi  = fe.hi*COULOMB_CONST;
    FE[iG] = fe;
}

__kernel void evalLJC_QZs(
    int nAtoms, 
    __global float4*   atoms,
    __global float2*    cLJs,
    __global float4*    poss,
    __global float8*    FE,
    float4 Qs,
    float4 QZs
){
    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);
    //if(iG==0){ printf("evalLJC: nAtoms: %i \n", nAtoms ); }
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = cLJs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            //if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
            if( (j+i0)<nAtoms ){
                float4 xyzq = LATOMS[j];
                fe.lo += getLJ     ( xyzq.xyz, LCLJS[j], pos );
                fe.hi += getCoulomb( xyzq, pos+(float3)(0,0,QZs.x) ) * Qs.x;
                fe.hi += getCoulomb( xyzq, pos+(float3)(0,0,QZs.y) ) * Qs.y;
                fe.hi += getCoulomb( xyzq, pos+(float3)(0,0,QZs.z) ) * Qs.z;
                fe.hi += getCoulomb( xyzq, pos+(float3)(0,0,QZs.w) ) * Qs.w;
            } 
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    fe.hi  = fe.hi*COULOMB_CONST;
    FE[iG] = fe;
}




float3 tipForce( float3 dpos, float4 stiffness, float4 dpos0 ){
    float r = sqrt( dot( dpos,dpos) );
    return  (dpos-dpos0.xyz) * stiffness.xyz        // harmonic 3D
         + dpos * ( stiffness.w * (r-dpos0.w)/r );  // radial
}


__kernel void relaxStrokesDirect(
    int nAtoms, 
    __global float4*    atoms,
    __global float2*    cLJs,
    __global float4*    points,
    __global float4*    FEs,
    float4 dTip,
    float4 stiffness,
    float4 dpos0,
    float4 relax_params,
    float Q,
    int nz
){

    __local float4 LATOMS[32];
    __local float2 LCLJS [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);


    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0.xyz; 
    
    const float dt   = relax_params.x;
    const float damp = relax_params.y;
    //printf( " %li (%f,%f,%f)  \n",  get_global_id(0), tipPos.x, tipPos.y, tipPos.z);

    for(int iz=0; iz<nz; iz++){
        float4 fe_;
        float3 vel   = 0.0f;

        
        for(int i=0; i<N_RELAX_STEP_MAX; i++){
            //fe        = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );

            // Get Atomic Forces
            float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);

            for (int i0=0; i0<nAtoms; i0+= nL ){
                int i = i0 + iL;
                //if(i>=nAtoms) break;  // wrong !!!!
                LATOMS[iL] = atoms[i];
                LCLJS [iL] = cLJs[i];
                barrier(CLK_LOCAL_MEM_FENCE);
                for (int j=0; j<nL; j++){
                    if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
            fe_ = fe.lo + fe.hi * Q;


            float3 f  = fe_.xyz;
            f        += tipForce( pos-tipPos, stiffness, dpos0 );
            vel      *=       damp;
            vel      += f   * dt;
            pos.xyz  += vel * dt;
            if(dot(f,f)<F2CONV) break;

        }

        /*
            // Get Atomic Forces
            float8 fe  = (float8) (0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f);

            for (int i0=0; i0<nAtoms; i0+= nL ){
                int i = i0 + iL;
                //if(i>=nAtoms) break;  // wrong !!!!
                LATOMS[iL] = atoms[i];
                LCLJS [iL] = cLJs[i];
                barrier(CLK_LOCAL_MEM_FENCE);
                for (int j=0; j<nL; j++){
                    if( (j+i0)<nAtoms ) fe += getLJC( LATOMS[j], LCLJS[j], pos );
                }
                barrier(CLK_LOCAL_MEM_FENCE);
            }
        fe_ = fe.lo + fe.hi * -0.1f;
        */

        FEs[get_global_id(0)*nz + iz] = fe_;
        //FEs[get_global_id(0)*nz + iz].xyz = pos;
        tipPos += dTip.xyz;
        pos    += dTip.xyz;

    }
    


}




__kernel void evalMorse(
    const int nAtoms, 
    __global float4*   atoms,
    __global float4*   REAs,
    __global float4*   poss,
    __global float4*   FE
){
    __local float4 lATOMs[32];
    __local float4 lREAs [32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    float3 pos = poss[iG].xyz;
    float4 fe  = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break; // wrong !!!!
        lATOMs[iL] = atoms[i];
        lREAs [iL] = REAs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ) fe += getMorse( pos - lATOMs[j].xyz, lREAs[j].xyz );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // http://www.informit.com/articles/article.aspx?p=1732873&seqNum=3
    FE[iG] = fe;
}


__kernel void evalLorenz(
    int nAtoms, 
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     FE
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float fe = 0.0f;

    //if(iG==0){ 
    //    printf( "nAtoms %i \n", nAtoms );
    //}
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){

            //if(iG==64*128 + 64) printf( "iatom %i atom (%g,%g,%g) pos(%g,%g,%g) \n", i0+j, LATOMS[j].x, LATOMS[j].y, LATOMS[j].z, pos.x,pos.y,pos.z );
            //if( (j+i0)<nAtoms ) fe += getLorenz( LATOMS[j], coefs[j], pos );
            if( (j+i0)<nAtoms ) fe += getLorenz( LATOMS[j], LCOEFS[j], pos );
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;
    //FE[iG] = pos.z;
}



__kernel void evalDisk(
    int nAtoms, 
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     FE,
    float dzmax,
    float dzmax_s,
    float offset,
    float Rpp,
    float4 rotA,
    float4 rotB,
    float4 rotC
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float fe = 0.0f;
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf( " xyzq (%g,%g,%g,%g) coef (%g,%g,%g,%g) \n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w );  } }
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp  =  pos - LATOMS[j].xyz;
                float3 abc = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                //float   R  = coefs[j].z;
                float   R  = LCOEFS[j].z;
                float   r2 = dot(abc.xy,abc.xy);
                float dxy2 = r2/( (R*R) );
                float  z   = abc.z - offset;
                float  z_s = abc.z - LCOEFS[j].w - Rpp;
                if( ( dxy2 < 1.0f ) && ( z < dzmax ) && ( z_s < dzmax_s ) ){
                    //fe += 1-dxy2;
                    //fe += ( 1-(abc.z/dzmax) ) * sqrt( 1- dxy2 );
                    fe += ( 1.0f-(z/dzmax) ) * ( 1.0f- sqrt(dxy2) );
                    //fe += coefs[j].w;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;
}

__kernel void evalSpheres(
    int nAtoms, 
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     FE,
    float Rpp,
    float zmin,
    float4 rotA,
    float4 rotB,
    float4 rotC
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;

    //float Rpp  =  1.0;
    //float Rpp  =  0.0;
    //float Rpp  = -0.7;

    float mask = 1.0;
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf( " xyzq (%g,%g,%g,%g) coef (%g,%g,%g,%g) \n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w );  } }

    float ztop = zmin;
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                //float  Rvdw  = coefs[j].w + Rpp;
                float  Rvdw  = LCOEFS[j].w + Rpp;
                //float  Rvdw  = coefs[j].w;
                //float  Rvdw  = coefs[j].z;
                float r2xy   =  dot(abc.xy,abc.xy);
                float  z     = -abc.z + sqrt( Rvdw*Rvdw - r2xy );
                if(z>ztop){
                    ztop=z;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = ztop;
}


__kernel void evalSphereCaps(
    int nAtoms, 
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     FE,
    float Rpp,
    float zmin,
    float tgMax,
    float tgWidth,
    float4 rotA,
    float4 rotB,
    float4 rotC
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;

    //float Rpp  =  1.0;
    //float Rpp  =  0.0;
    //float Rpp  = -0.7;

    //float mask = 1.0;
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf( " xyzq (%g,%g,%g,%g) coef (%g,%g,%g,%g) \n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w );  } }

    float ztop  = zmin;
    //float ztop  = 0;
    float tgtop = -1.0;
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                //float  Rvdw  = coefs[j].w + Rpp;
                float  Rvdw  = LCOEFS[j].w + Rpp;
                float r2xy   =  dot(abc.xy,abc.xy);

                float dz     =  sqrt( Rvdw*Rvdw - r2xy );
                float  z     = -abc.z + dz;

                //float tgWidth = 0.1f;
                float tg      = sqrt(r2xy)/Rvdw;
                z  = zmin + (z-zmin)*(1.0f-smoothstep( tgMax-tgWidth, tgMax, tg ));
                ztop = fmax( z, ztop );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = ztop;
    
}


__kernel void evalDisk_occlusion(
    int nAtoms, 
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     FE,
    float Rpp,
    float zmin,
    float zmargin,
    float dzmax,
    float4 rotA,
    float4 rotB,
    float4 rotC
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;

    //float dzmax = 0.2;
    //float Rpp  =  1.0;
    //float Rpp  =  0.0;
    //float Rpp  = -0.7;

    float mask = 1.0;
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf( " xyzq (%g,%g,%g,%g) coef (%g,%g,%g,%g) \n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w );  } }

    float ztop = zmin;
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                //float  Rvdw  = coefs[j].w + Rpp;
                float  Rvdw  = LCOEFS[j].w + Rpp;
                float  r2xy   =  dot(abc.xy,abc.xy);
                float z = 2.0*Rvdw - r2xy/(Rvdw*Rvdw) - abc.z; // osculation parabola to atom
                if(z>ztop){
                    ztop=z;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //FE[iG] = ztop;

    float fe   = 0.0f;
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                //float  Rvdw  = coefs[j].w + Rpp;
                float  Rvdw  = LCOEFS[j].w + Rpp;
                float  r2xy   =  dot(abc.xy,abc.xy);
                float z   = 2.0f*Rvdw - r2xy/(Rvdw*Rvdw) - abc.z; // osculation parabola to atom
                float cdz = 1.0 - ((ztop-z)/zmargin);
                //fe += cdz;
                //fe = fmax(ztop,z);
                //fe = fmax(cdz,fe);
                if( cdz>0.0f ){
                    //float   R     = coefs[j].z * 2.0;
                    float   R     = LCOEFS[j].z * 2.0;
                    float dxy2 = r2xy/( R*R );
                    if( ( dxy2 < 1.0f ) && ( abc.z < dzmax ) ){
                        //fe += 1-dxy2;
                        //fe += ( 1-(abc.z/dzmax) ) * sqrt( 1- dxy2 );
                        fe += ( 1.0f-(abc.z/dzmax) ) * ( 1.0f- sqrt(dxy2) ) * cdz;
                        //fe += coefs[j].w;
                    }
                }
            }
        }

        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = fe;


}



__kernel void evalMultiMapSpheres(
    int nAtoms,
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     MultMap,
    float Rpp,
    float zmin,
    float4 rotA,
    float4 rotB,
    float4 rotC,
    int bOccl,
    int   nChan,
    float Rmin,
    float Rstep
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    float ztops[20];

    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;

    float mask = 1.0;
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf( " xyzq (%g,%g,%g,%g) coef (%g,%g,%g,%g) \n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w );  } }

    float ztop = zmin;
    for (int i=0; i<nChan; i++){
        ztops[i]=zmin;
    }

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                float  Rvdw  = LCOEFS[j].w + Rpp;
                float r2xy   =  dot(abc.xy,abc.xy);
                float  z     = -abc.z + sqrt( Rvdw*Rvdw - r2xy );

                int ityp = (int)( (LCOEFS[j].w - Rmin)/Rstep );
                if(ityp>=nChan)ityp=nChan-1;
                if(z>ztop){
                    ztop=z;
                }
                if(z>ztops[ityp]){
                    ztops[ityp  ]=z;
                    //ztops[ityp+4]=z;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    for (int i=0; i<nChan; i++){
        //ztops[i]=zmin;
        if(bOccl>0){ if( ztop-ztops[i] < 0.02 ){ MultMap[iG*nChan+i] = ztops[i]; } else { MultMap[iG*nChan+i] = zmin; }; }
        else       { MultMap[iG*nChan+i] = ztops[i]; }
    }
}

__kernel void evalSpheresType(
    int nAtoms,
    int nType,
    __global float4*    atoms,
    __global int*       itypes,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     MultMap,
    float Rpp,
    float zmin,
    float4 rotA,
    float4 rotB,
    float4 rotC,
    int bOccl
){
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    float ztops[20];

    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    float3 pos = poss[iG].xyz;
    float mask = 1.0;
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf( " xyzq (%g,%g,%g,%g) coef (%g,%g,%g,%g) \n", atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   coefs[i].x,coefs[i].y,coefs[i].z,coefs[i].w );  } }
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf("itypes[%i] = %i \n" , i, itypes[i] ); } }
    
    float ztop = zmin;
    for (int i=0; i<nType; i++){
        ztops[i]=zmin;
    }

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                float  Rvdw  = LCOEFS[j].w + Rpp;
                float r2xy   =  dot(abc.xy,abc.xy);
                float  z     = -abc.z + sqrt( Rvdw*Rvdw - r2xy );

                int ityp = itypes[j+i0];
                //if( iG==0 ){ printf( " %i : %i ityp %i \n", j+i0, itypes[j+i0], ityp ); }
                //ityp = 0;
                if(z>ztop){
                    ztop=z;
                }
                if(ityp>=0){
                    if(z>ztops[ityp]){
                        ztops[ityp]=z;
                    }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    for (int i=0; i<nType; i++){
        if(bOccl>0){ if( ztop-ztops[i] < 0.02 ){ MultMap[iG*nType+i] = ztops[i]; } else { MultMap[iG*nType+i] = zmin; }; }
        else       { MultMap[iG*nType+i] = ztops[i]; }
    }
}

__kernel void evalBondEllipses(
    int nBonds,
    __global float8*    bondPoints,
    __global float4*    poss,
    __global float*     out,
    __constant float*   Rfunc,     // definition of radial function
    float drStep,
    float Rmax,
    float elipticity,
    float zmin,
    float4 rotA,
    float4 rotB,
    float4 rotC
){
    #define  nFuncMax    128
    __local float  RFUNC [nFuncMax];
    __local float8 LATOMS[32];

    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);

    //float3 pos = poss[iG].xyz;
    float3 pos = rotMat ( poss[iG].xyz,  rotA.xyz, rotB.xyz, rotC.xyz );
    float invStep = 1/drStep;
    //if( iG==0 ){ for(int i=0; i<nBonds; i++){ printf( " a1(%g,%g,%g,%g) a2(%g,%g,%g,%g) \n", bondPoints[i].lo.x,bondPoints[i].lo.y,bondPoints[i].lo.z,bondPoints[i].lo.w,   bondPoints[i].hi.x,bondPoints[i].hi.y,bondPoints[i].hi.z,bondPoints[i].hi.w );  } }
    //if( iG==0 ){ for(int i=0; i<nAtoms; i++){ printf("itypes[%i] = %i \n" , i, itypes[i] ); } }

    // local copy of Rfunc
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i0=0; i0<nFuncMax; i0+= nL ){ int i = i0+iL; RFUNC[i] = Rfunc[i]; }
    barrier(CLK_LOCAL_MEM_FENCE);

    float R2max = Rmax*Rmax;
    float ftot=0;
    //nBonds = 1;
    for (int i0=0; i0<nBonds; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = bondPoints[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nBonds ){ 
                float3 p1 = rotMat ( LATOMS[j].lo.xyz,  rotA.xyz, rotB.xyz, rotC.xyz );
                float3 p2 = rotMat ( LATOMS[j].hi.xyz,  rotA.xyz, rotB.xyz, rotC.xyz );
                float2 hbond   = p2.xy - p1.xy;
                float  rbond  = length(hbond);
                hbond /= rbond;
                float2 d  = pos.xy - p1.xy;
                float  c  = dot(d,hbond);
                float2 dT = d   - hbond*c;
                //float dl = ( fabs(c) + fabs(rbond-c) - rbond )*0.5;
                //float r2 = dot(dT,dT) + dl*dl;
                //float r2 = dot(dT,dT)* max(0.0, 4.0*c*(rbond-c)/(rbond*rbond)   );
                float c_ = (c-0.5*rbond)*elipticity;
                float r2 = dot(dT,dT) + c_*c_ ;
                float r  = sqrt(r2); 
                //float z       = reciprocalInterp( p1.z, p2.z, r2xy1, r2xy2  );
                float z =  mix ( p1.z, p2.z, c/rbond ) - pos.z;
                if( (r < Rmax) && (z>zmin) ){
                    float fi   = spline_Hermite_y_arr( r*invStep, RFUNC );
                    float zfunc = max( (zmin-z)/zmin, 0.0f );
                    fi*=zfunc;
                    //if( iG==0 ){ printf( " i0 %i j %i   fi %g   r2xy %g  beta %g   dp1(%g,%g,%g), dp1(%g,%g,%g) \n", i0, j,  fi, r2xy,    beta,   dp1.x,dp1.y,dp1.z, dp2.x,dp2.y,dp2.z ); }
                    //if( iG==0 ){ printf( " i0 %i j %i   z,zmin(%g,%g) zfunc %g    fi %g  fi*zfunc %g  rxy %g   sx %g  \n", i0, j,  z, zmin,  zfunc,  fi, zfunc*fi, r, r*invStep ); }
                    //ftot += fi  ;
                    ftot = fmax( fi, ftot );
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    out[iG] = ftot;
}


__kernel void evalAtomRfunc(
    int nAtoms, 
    __global float4*    atoms,
    __global float4*    coefs,
    __global float4*    poss,
    __global float*     FE,
    __constant float*   Rfunc,     // definition of radial function
    float drStep,
    float Rmax,
    float zmin,
    float4 rotA,
    float4 rotB,
    float4 rotC
){
    #define  nFuncMax    128
    __local float  RFUNC [nFuncMax];
    __local float4 LATOMS[32];
    __local float4 LCOEFS[32];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    float3 pos = poss[iG].xyz;
    float invStep = 1/drStep;

    // local copy of Rfunc
    barrier(CLK_LOCAL_MEM_FENCE);
    for (int i0=0; i0<nFuncMax; i0+= nL ){ int i = i0+iL; RFUNC[i] = Rfunc[i]; }
    barrier(CLK_LOCAL_MEM_FENCE);

    float ftot = 0;
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float3 dp    = pos - LATOMS[j].xyz;
                float3 abc   = (float3)( dot(dp,rotA.xyz), dot(dp,rotB.xyz), dot(dp,rotC.xyz) );
                float  RvdW  = LCOEFS[j].w;
                float r2xy   =  dot(abc.xy,abc.xy);
                float r      = sqrt(r2xy);
                float  z     = -dp.z + RvdW - 1.4;  // subtract hydrogen ?
                if( (r < Rmax) && (z>zmin) ){
                    //float fi = exp( beta*rxy );
                    float fi   = spline_Hermite_y_arr( r*invStep, RFUNC );
                    float zfunc = max( (zmin-z)/zmin, 0.0f );
                    fi*=zfunc;
                    //if( iG==0 ){ printf( " i0 %i j %i   fi %g   r2xy %g  beta %g   dp1(%g,%g,%g), dp1(%g,%g,%g) \n", i0, j,  fi, r2xy,    beta,   dp1.x,dp1.y,dp1.z, dp2.x,dp2.y,dp2.z ); }
                    //if( iG==0 ){ printf( " i0 %i j %i   z,zmin(%g,%g) zfunc %g    fi %g  fi*zfunc %g  r %g   sx %g  \n", i0, j,  z, zmin,  zfunc,  fi, zfunc*fi, r, r*invStep ); }
                    //ftot += fi  ;
                    ftot = fmax( fi, ftot );
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    FE[iG] = ftot;
}
