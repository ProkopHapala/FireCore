#define R2SAFE          1e-4f
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

#define N_LOCAL  32
#define NORB_MAX 64

//__constant sampler_t sampler_wrf =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_wrf =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_wrf =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;
__constant sampler_t sampler_wrf =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
//__constant sampler_t sampler_wrf =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;
//__constant sampler_t sampler_wrf =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

__constant sampler_t sampler_2 =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;

__constant float pref_s = 0.28209479177; // sqrt(1.0f/(4.0f*M_PI));
__constant float pref_p = 0.4886025119;  //sqrt(3.0f/(4.0f*M_PI));
__constant float pref_d = 1.09254843059; //sqrt(15.0f/(4.0f*M_PI));
__constant float wf_tiles_per_angstroem = 20.0f;

int wrap( int i, int n ){ if(i<0)i+=n; if(i>=n)i-=n; return i; }

/*
float4 read_imagef_trilin( __read_only image3d_t imgIn, float4 coord ){
    float4 d = (float4)(0.00666666666,0.00666666666,0.00666666666,1.0); 
    float4 icoord;
    float4 fc     =  fract( coord/d, &icoord );
    icoord*=d;
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    // NOTE AMD-GPU seems to not accept CLK_NORMALIZED_COORDS_FALSE
    //return read_imagef( imgIn, sampler_2, icoord );
    //return read_imagef( imgIn, sampler_1, coord );
    return  
     (( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,0.0,0.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,d.y,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,d.y,0.0,0.0) ) * fc.x )*fc.y )*mc.z
    +(( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,d.z,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,0.0,d.z,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,d.y,d.z,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,d.y,d.z,0.0) ) * fc.x )*fc.y )*fc.z;
}; 
*/

/*
float4 read_imagef_lin( __read_only image3d_t imgIn, float x, float y ){
    float4 d = (float4)(0.00666666666,0.00666666666,0.00666666666,1.0); 
    float4 icoord;
    float4 fc     =  fract( coord/d, &icoord );
    icoord*=d;
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    return  
     (( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(,0.0,0.0,0.0) ) * fc.x )
}; 

float4 interpFE_prec( float2 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef_trilin( imgIn, coord ); 
    // read_imagef( imgIn, sampler_1, coord );
    //return coord;
}

float4 interpFE( float2 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef( imgIn, sampler_1, coord );
    //return coord;
}
*/

float2 lerp_basis( float x, float slot, __read_only image2d_t imgIn ){
    //float d = 0.01;
    float icoord;
    float fc     =  fract( x, &icoord );
    //printf( "iG %i x %g | iu %f du %f \n", iG, x, icoord, fc );
    return read_imagef( imgIn, sampler_wrf, (float2){ icoord  , slot } ).xy*(1.f-fc)
        +  read_imagef( imgIn, sampler_wrf, (float2){ icoord+1, slot } ).xy*     fc ;
    //return read_imagef( imgIn, sampler_wrf, (float2){ x  , slot } )*(1.f-fc)
    //    +  read_imagef( imgIn, sampler_wrf, (float2){ x+1, slot } )*    fc ;
    //return coord;
}

__kernel void mul(
    const int N,
    __global float* A,
    __global float* B,
    __global float* out
){
    const size_t i = get_global_id(0);
    //if(i==0)printf( "DEBUG_GPU mul N %i \n", N );
    if(i<N){ 
        out[i] = A[i] * B[i]; 
        //out[i] = sin( i*0.1 ); 
        //if((i/100)<100)printf( "DEBUG_GPU mul[%i] A=%g B=%g out=%g \n", i, A[i], B[i], out[i] );
    }
};

/*
__kernel void roll(
    __global float2*  InBuff,    //4
    __global float2*  OutBuff,   //4
    int4   shift,                //5
    int4   nGrid                //6
){
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;
    int ja = ia + shift.x; 
    int jb = ib + shift.y;
    int jc = ic + shift.z;
    int jG = ia + nGrid.x*( ib + nGrid.y*ic);
    OutBuff[iG] = InBuff[jG];
}
*/

__kernel void roll(
    __global float2*  InBuff,    //1
    __global float2*  OutBuff,   //2
    int4   shift,                //3
    int4   nGrid                 //4
){
    const int ia  = get_global_id (0);
    const int ib  = get_global_id (1);
    const int ic  = get_global_id (2);
    const int na  = get_global_size (0);
    const int nb  = get_global_size (1);
    const int nc  = get_global_size (2);
    int iG = ia + nGrid.x*( ib + nGrid.y*ic);
    if(iG==0){ printf( "GPU roll() size(%i,%i,%i) nGrid(%i,%i,%i) shift(%i,%i,%i)\n", na,nb,nc,   nGrid.x,nGrid.y,nGrid.z, shift.x,shift.y,shift.z );}
    //OutBuff[iG] = InBuff[iG];
    if( (ia<nGrid.x) && (ib<nGrid.y) && (ic<nGrid.z) ){
        int ja = wrap( ia + shift.x, nGrid.x );
        int jb = wrap( ib + shift.y, nGrid.y );
        int jc = wrap( ic + shift.z, nGrid.z );
        int iG = ia + nGrid.x*( ib + nGrid.y*ic);
        int jG = ja + nGrid.x*( jb + nGrid.y*jc);
        OutBuff[iG] = InBuff[jG];
    }
}

__kernel void lincomb(
    const int N,
    __global float* A,
    __global float* B,
    __global float* out,
    const float2 coefs
){
    const size_t i = get_global_id(0);
    if(i<N){ 
        out[i] = A[i]*coefs.x + B[i]*coefs.y; 
    }
}

// https://github.com/ProkopHapala/ProbeParticleModel/blob/OpenCL_py3/pyProbeParticle/fieldFFT.py
//
//  E = np.real(np.fft.ifftn(convFFT * (dd[0]*dd[1]*dd[2]) / (detLmatInv) ) )
//

__kernel void poissonW(
    const int   N,
    __global float2* A,
    __global float2* out,
    const float4 dCell
){
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    //const int iw = get_global_id(3);
    const int nx = get_global_size(0);
    const int ny = get_global_size(1);
    const int nz = get_global_size(2);
    //const int nw = get_global_size(3);
    //if( (ix==0)&&(iy==0)&&(iz==0) ){ printf( "GPU poissonW size(%i,%i,%i) n %i dCell(%g,%g,%g)\n", nx,ny,nz, N, dCell.x,dCell.y,dCell.z  ); };
    int i = ix + nx*( iy + ny*iz );
    //float4 k = (float4){ dCell.x*ix, dCell.y*iy, dCell.z*iz, 0};
    float4 k = (float4){ ix/(0.5f*nx), iy/(0.5f*ny), iz/(0.5f*nz), 0};
    k = 1.0f-fabs(k-1.0f);
    float  f = 1/dot( k, k ); 
    if(i==0)f=0;
    if(i<N){ 
        out[i] = A[i]*f;
        //out[i] = f;
        //out[i] = k.x;
    }
    //if( (ix==(nx/2))&&(iy==(ny/2)) ){ printf( "GPU iz,i[%i,%i|%i] k %g A[i] %g out %g \n", iz, i, N, f, A[i].x, out[i].x ); };
};

#define ixyz(ix,iy,iz)  ix + nx*( iy + ny*iz )

__kernel void gradient(
    const int4       off,    
    __global float2* A,
    __global float4* out,
    const float4     inv_d
){
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    //const int iw = get_global_id(3);
    const int nx = get_global_size(0);
    const int ny = get_global_size(1);
    const int nz = get_global_size(2);
    //__local float4 LATOMS[64];    // ToDo - LATER use local memory ?
    int i = ixyz(ix,iy,iz);
    int jx = ix+off.x;
    int jy = iy+off.y;
    int jz = iz+off.z;
    //int j = ixyz(ix,iy,iz);
    out[ i ] = (float4){
        (A[ ixyz(ix+1,iy  ,iz  ) ].x - A[ixyz(ix-1,iy  ,iz  )].x)*inv_d.x,
        (A[ ixyz(ix  ,iy+1,iz  ) ].x - A[ixyz(ix  ,iy-1,iz  )].x)*inv_d.y,
        (A[ ixyz(ix  ,iy  ,iz+1) ].x - A[ixyz(ix  ,iy  ,iz-1)].x)*inv_d.z,
        A[ ixyz(ix,iy,iz) ].x
    }; 
};


// Grid projection

float sp3( float3 dp, float4 c, float beta ){
    float   r2  = ( dot(dp,dp) +  R2SAFE );
    float   r   = sqrt(r2);
    float   ir  = 1/r;
    float   fr  = exp( -beta*r );
    float4  v = (float4) ( dp*ir, 1 );
    //v*=fr;
    return  dot( v ,c )*fr;
}

float4 sp3_tex( float3 dp, float slot, __read_only image2d_t imgIn ){
    float   r2  = ( dot(dp,dp) +  R2SAFE );
    float   r   = sqrt(r2);
    float   ir  = 1/r;
    float4  fr  = lerp_basis( r, slot, imgIn ).yyyx;
    return (float4) ( dp*ir*pref_p, pref_s )*fr;
}

__kernel void projectAtomsToGrid(
    const int nAtoms,            //1
    __global float4*  atoms,     //2
    __global float4*  coefs,     //3
    __global float2*  outGrid,   //4
    int4   nGrid,                //5
    float4 grid_p0,              //6
    float4 grid_dA,              //7
    float4 grid_dB,              //8
    float4 grid_dC               //9
){
    __local float4 LATOMS[N_LOCAL];
    __local float4 LCOEFS[N_LOCAL];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;
    
    if(iG>nMax) return;
    //if(iG==0) for(int ia=0; ia<nAtoms;ia++)printf( "DEBUG_GPU atoms[%i](%g,%g,%g,%g) coefs[0](%g,%g,%g,%g) \n", ia, atoms[ia].x,atoms[ia].y,atoms[ia].z,atoms[ia].w,  coefs[ia].x,coefs[ia].y,coefs[ia].z,coefs[ia].w  );
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;
    outGrid[iG].x = sin(pos.x+pos.y+pos.z);
    //outGrid[iG].x = ia*grid_dA.x;
    //outGrid[iG].x = 1;
    //outGrid[iG] = nAtoms;
    //return;
    float2 wf   = (float2) (0.0f,0.0f);
    //wf.x += sp3( pos-(float3){0.0,0.0,0.0}, (float4){0.0,0.0,0.0,1.0}, 1.0 );
    //wf.x += sp3( pos-atoms[0].xyz, (float4){0.0,0.0,0.0,1.0}, 1.0 );
    //wf.x = atoms[0].z;
    
    //for (int i=0; i<nAtoms; i++ ){
    //    wf.x += sp3( pos-atoms[i].xyz, coefs[i], 1.0 );
    //}
    
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        // TODO : we can optimize it here - load just the atom which are close to center of the block !!!!!
        //        BUT it will be better/faster to do this on CPU ... atoms-buffer should be already segment to blocks each considering to local work-group (3D tile)
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float4 xyzq  = LATOMS[j];
                float4 cs    = LCOEFS[j];
                //fe += getCoulomb( xyzq, pos+(float3)(0,0,QZs.x) ) * Qs.x;
                wf.x += sp3( pos-xyzq.xyz, cs, xyzq.w );
                //wf.x += sp3( pos-atoms[i].xyz, coefs[i], 1.0 );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    //if ( (ia==75)&&(ib==75) ) { printf(" iz %i fe %g,%g,%g,%g \n", ic, fe.x, fe.y, fe.z, fe.w ); }
    outGrid[iG] = wf;
}

__kernel void projectAtomsToGrid_texture(
    const int nAtoms,            //1
    __global float4*  atoms,     //2
    __global float4*  coefs,     //3
    __global float2*  outGrid,   //4
    __read_only image2d_t imgIn, //5 
    int4   nGrid,                //6
    float4 grid_p0,              //7
    float4 grid_dA,              //8
    float4 grid_dB,              //9
    float4 grid_dC               //10
){
    __local float4 LATOMS[N_LOCAL];
    __local float4 LCOEFS[N_LOCAL];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;
    
    if(iG>nMax) return;
    //if(iG==0) for(int ia=0; ia<nAtoms;ia++)printf( "DEBUG_GPU atoms[%i](%g,%g,%g,%g) coefs[0](%g,%g,%g,%g) \n", ia, atoms[ia].x,atoms[ia].y,atoms[ia].z,atoms[ia].w,  coefs[ia].x,coefs[ia].y,coefs[ia].z,coefs[ia].w  );
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;

    float2 wf   = (float2) (0.0f,0.0f);
    //wf.x = read_imagef( imgIn, sampler_wrf, pos.xy*10.f ).x;
    //wf.x = read_imagef( imgIn, sampler_wrf, (float2){ pos.x*10.0f, 0.5 } ).x;

    //if(iG==0){    for(int i=0; i<100; i++){ wf = read_imagef( imgIn, sampler_wrf, (float2){ i*10.f, 0.5 } ).xy; printf( "wf [%i] %g %g \n", i, wf.x, wf.y ); };}

    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        // TODO : we can optimize it here - load just the atom which are close to center of the block !!!!!
        //        BUT it will be better/faster to do this on CPU ... atoms-buffer should be already segment to blocks each considering to local work-group (3D tile)
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float4 xyzq  = LATOMS[j];
                float4 cs    = LCOEFS[j];
                //wf.x += sp3( pos-xyzq.xyz, cs, xyzq.w );
                //wf.x += sp3_tex( (pos-xyzq.xyz)*10.f, cs, xyzq.w, imgIn );
                //wf.x = read_imagef( imgIn, sampler_wrf, pos.xy*8 ).x;
                //wf.x +=        sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem, cs, xyzq.w, imgIn );
                wf.x += dot( cs, sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem,     xyzq.w, imgIn ) );
                //if((ia==16)&&(ib==16)){ printf("DEBUG_GPU pos(%g,%g,%g) xyzq(%g,%g,%g,%g) cs(%g,%g,%g,%g) wf(%g,%g) \n", pos.x,pos.y,pos.z,  xyzq.x,xyzq.y,xyzq.z,xyzq.w, cs.x,cs.y,cs.z,cs.w, wf.x, wf.y ); }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //if((ia==16)&&(ib==16)){ printf("DEBUG_GPU pos(%g,%g,%g) wf(%g,%g) \n", pos.x,pos.y,pos.z, wf.x, wf.y ); }
    outGrid[iG] = wf;
}




__kernel void projectOrbDenToGrid_texture(
    const int nAtoms,            //1
    const int iorb0,            //2
    const int iorb1,            //3
    __global float4*  atoms,     //4
    __global float4*  coefs,     //5
    __global float2*  outGrid,   //6
    __read_only image2d_t imgIn, //7 
    int4   nGrid,                //8
    float4 grid_p0,              //9
    float4 grid_dA,              //10
    float4 grid_dB,              //11
    float4 grid_dC,              //12
    float2 acumCoef
){
    __local float4 LATOMS[N_LOCAL];
    __local float4 LCOEFS[N_LOCAL];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;

    if(iG==0){  printf("!!!!!!!!!!!!!!!! projectOrbDenToGrid_texture acumCoef %g,%g nAtoms %i iorb(%i,%i)  nMax %i \n", acumCoef.x, acumCoef.y, nAtoms, iorb0, iorb1, nMax ); }
    if(iG>nMax) return;

    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;

    float dens = 0.0;
    
    // ToDo : Later we have to change the order of the loops
    for(int iorb=iorb0; iorb<iorb1; iorb++){
        //if(iG==0){ printf( "GPU iorb %i \n", iorb ); }
        int icoef0 = iorb*nAtoms;
        float2 wf   = (float2) (0.0f,0.0f);
        for (int i0=0; i0<nAtoms; i0+=nL ){
            int i = i0 + iL;
            LATOMS[iL] = atoms[i];
            LCOEFS[iL] = coefs[i+icoef0];
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int j=0; j<nL; j++){
                if( (j+i0)<nAtoms ){ 
                    float4 xyzq  = LATOMS[j];
                    float4 cs    = LCOEFS[j];
                    //wf.x +=        sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem, cs, xyzq.w, imgIn );
                    wf.x += dot( cs, sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem,     xyzq.w, imgIn ) );
                }
            } // j
            barrier(CLK_LOCAL_MEM_FENCE);
        } // i0
        dens += wf.x*wf.x;
    } // iorb
    //if(iG==0){ printf( "GPU loop DONE ! \n" ); }
    //outGrid[iG] = (float2){dens,0.0f};
    if(fabs(acumCoef.x)<1e-8){
        outGrid[iG] = (float2){dens,0.0f}*acumCoef.y;
    }else{
        outGrid[iG] = outGrid[iG]*acumCoef.x + ((float2){dens,0.0f})*acumCoef.y;
    }
    //if(iG==0){ printf( "GPU all DONE ! \n" ); }
    //if(iG==0){ printf("projectOrbDenToGrid_texture END \n"); }
}

__kernel void projectAtomDenToGrid_texture(
    const int nAtoms,            //1
    __global float4*  atoms,     //4
    __global float4*  coefs,     //5
    __global float2*  outGrid,   //6
    __read_only image2d_t imgIn, //7 
    int4   nGrid,                //8
    float4 grid_p0,              //9
    float4 grid_dA,              //10
    float4 grid_dB,              //11
    float4 grid_dC,              //12
    float2 acumCoef
){
    __local float4 LATOMS[N_LOCAL];
    __local float4 LCOEFS[N_LOCAL];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
   
    const int nab  =  nGrid.x*nGrid.y;
    const int ia   =  iG%nGrid.x; 
    const int ib   = (iG%nab)/nGrid.x;
    const int ic   =  iG/nab; 
    const int nMax =  nab*nGrid.z;
    //if(iG==0){  printf("GPU_DEBUG projectAtomDenToGrid_texture() acumCoef %g,%g nAtoms %i nMax %i \n", acumCoef.x, acumCoef.y, nAtoms, nMax ); }
    if(iG>nMax) return;
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;
    float dens = 0.0;
    for (int i0=0; i0<nAtoms; i0+=nL ){
        int i = i0 + iL;
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float4 xyzq  = LATOMS[j];
                float4 cs    = LCOEFS[j];
                float4 wfs = sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem,     xyzq.w, imgIn );
                wfs*=wfs;  
                dens += dot( cs, wfs );
            }
        } // j
        barrier(CLK_LOCAL_MEM_FENCE);
    } // i0
    //if(iG==0){ printf( "GPU loop DONE ! \n" ); }
    //outGrid[iG] = (float2){dens,0.0f};
    if(fabs(acumCoef.x)<1e-8){
        outGrid[iG] = (float2){dens,0.0f};
    }else{
        outGrid[iG] = outGrid[iG]*acumCoef.x + ((float2){dens,0.0f})*acumCoef.y;
    }
    //if(iG==0){ printf( "GPU all DONE ! \n" ); }
    //if(iG==0){ printf("projectOrbDenToGrid_texture END \n"); }
}


__kernel void projectOrbDenToGrid_texture_2(
    const int nAtoms,            //1
    const int nOrb,              //2
    __global float4*  atoms,     //4
    __global float4*  Cijs,      //5
    __global float2*  outGrid,   //6
    __read_only image2d_t imgIn, //7 
    int4   nGrid,                //8
    float4 grid_p0,              //9
    float4 grid_dA,              //10
    float4 grid_dB,              //11
    float4 grid_dC               //12
){
    // NOTE: 
    //  * each workgroup works on the same voxel
    //  * 
    //__local float4 LATOMS[N_LOCAL];
    //__local float4 LCOEFS[N_LOCAL];
    __local float4 LBAS[N_LOCAL];
    __local float  LDEN[N_LOCAL];
    //__local float4 LMOS[NORB_MAX];
    const int iL = get_local_id(0);
    if(iL>nAtoms) return;

    // === evaluate grid position
    const int iG = get_group_id(0);
    //const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    float3 pos    = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;

    // === evaluate basisfunctions of each atom
    float4 xyzq  = atoms[iL];
    float4 bas_i = sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem, xyzq.w, imgIn );
    LBAS[iL]     = bas_i;    
    barrier(CLK_LOCAL_MEM_FENCE);
    // === sum density matrix
    float dens_i = 0;
    int k = iL*nAtoms*4;
    for(int j=0; j<nAtoms; j++){
        float4 bas_j = LBAS[iL];
        dens_i += bas_i.x * dot( bas_j, Cijs[k  ] )   // TODO : Cijs should be also pre-loaded to local memory
               +  bas_i.y * dot( bas_j, Cijs[k+1] )
               +  bas_i.z * dot( bas_j, Cijs[k+2] )
               +  bas_i.w * dot( bas_j, Cijs[k+3] );
        k+=4;
    }
    LDEN[iL]=dens_i;
    barrier(CLK_LOCAL_MEM_FENCE);
    // === REDUCE from local memory 
    float dens=0;
    if(iL==0){
        for(int i=0; i<nAtoms; i++) dens+=LDEN[i];
        outGrid[iG] = (float2){dens,0.0f};
    }
}



__kernel void projectWfAtPoints_tex(
    const int nAtoms,            //1
    __global float4*  atoms,     //2
    __global float4*  coefs,     //3
    int   nPos,                  //4
    __global float4*  poss,      //5
    __global float2*  out,       //6
    __read_only image2d_t imgIn  //7
){
    __local float4 LATOMS[N_LOCAL];
    __local float4 LCOEFS[N_LOCAL];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    
    if(iG>nPos) return;
    float3 pos  = poss[iG].xyz;
    float2 wf   = (float2) (0.0f,0.0f);
    //printf( "projectWfAtPoints_tex %i (%g, %g,%g) \n", iG, poss[iG].x, poss[iG].y, poss[iG].z  );
    //if(iG==0){ for(int i=0; i<nAtoms; i++){ printf( "atom[%i] atom(%g,%g,%g,%g) coefs(%g,%g,%g,%g)\n", i, atoms[i].x, atoms[i].y, atoms[i].z,atoms[i].w,  coefs[i].x, coefs[i].y, coefs[i].z,coefs[i].w ); }  }
    //if(iG==0){ for(int i=0; i<30; i++){ float x=i*0.1; float2 xy=lerp_basis( x*wf_tiles_per_angstroem, 1.1, imgIn ); printf( "i %i x %g (%g,%g) \n",i, x, xy.x, xy.y ); } }
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        // TODO : we can optimize it here - load just the atom which are close to center of the block !!!!!
        //        BUT it will be better/faster to do this on CPU ... atoms-buffer should be already segment to blocks each considering to local work-group (3D tile)
        LATOMS[iL] = atoms[i];
        LCOEFS[iL] = coefs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                float4 xyzq  = LATOMS[j];
                float4 cs    = LCOEFS[j];
                //wf.x += sp3( pos-xyzq.xyz, cs, xyzq.w );
                //wf.x += sp3_tex( (pos-xyzq.xyz)*10.f, cs, xyzq.w, imgIn );
                //wf.x = read_imagef( imgIn, sampler_wrf, pos.xy*8 ).x;
                float3 dp = (pos-xyzq.xyz)*wf_tiles_per_angstroem; ///0.529177210903f);
                //printf( "iG %i x %g r %g \n", iG, pos.x, length(dp) );
                //wf.x +=         sp3_tex( dp, cs, xyzq.w, imgIn );
                wf.x += dot( cs,  sp3_tex( dp,     xyzq.w, imgIn ) );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    out[iG] = wf;
}