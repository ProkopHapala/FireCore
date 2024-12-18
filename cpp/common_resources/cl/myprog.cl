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
        //out[i] = A[i];
        //out[i] = B[i];
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

__kernel void makeForceField(
    __global   float2*  InBuff,   //1
    __global   float4*  OutBuff,  //2
    //__write_only image3d_t imgOut,  // 2
    float4  mask,                   //3
    int4   nGrid                    //4
){
    const int ia  = get_global_id (0);
    const int ib  = get_global_id (1);
    const int ic  = get_global_id (2);
    float4 FE=(float4){0.0f,0.0f,0.0f,0.0f};
    int nab = nGrid.x*nGrid.y;
    if( (ia<nGrid.x) && (ib<nGrid.y) && (ic<nGrid.z) ){
        int i0; 
        int i1 = nGrid.x*( ib + nGrid.y*ic);
        FE.x =  mask.y * InBuff[ i0 + wrap( ia+2, nGrid.x) ].x
             +  mask.x * InBuff[ i0 + wrap( ia+1, nGrid.x) ].x
             -  mask.x * InBuff[ i0 + wrap( ia-1, nGrid.x) ].x
             -  mask.y * InBuff[ i0 + wrap( ia-2, nGrid.x) ].x;
        int i2 = ia + nab*ic;
        FE.y =  mask.y * InBuff[ i0 + wrap( ib+2, nGrid.y)*nGrid.x ].x
             +  mask.x * InBuff[ i0 + wrap( ib+1, nGrid.y)*nGrid.x ].x
             -  mask.x * InBuff[ i0 + wrap( ib-1, nGrid.y)*nGrid.x ].x
             -  mask.y * InBuff[ i0 + wrap( ib-2, nGrid.y)*nGrid.x ].x;
        int i3 = ia + nGrid.x*ib;
        FE.x =  mask.y * InBuff[ i0 + wrap( ic+2, nGrid.z)*nab ].x
             +  mask.x * InBuff[ i0 + wrap( ic+1, nGrid.z)*nab ].x
             -  mask.x * InBuff[ i0 + wrap( ic-1, nGrid.z)*nab ].x
             -  mask.y * InBuff[ i0 + wrap( ic-2, nGrid.z)*nab ].x;
        int iG = nGrid.x*( ib + nGrid.y*ic);
        OutBuff[iG] = FE;
        //write_imagef( imgOut, (int4){ia,ib,ic,0}, FE );
    }
}

#define ixyz(ix,iy,iz)  (ix) + nx*( (iy) + ny*(iz) )

__kernel void gradient(
    const int4       nGrid,    
    __global float2* A,
    __global float4* out,
    const float4     inv_d
){
    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    const int nx = nGrid.x;
    const int ny = nGrid.y;
    const int nz = nGrid.z;
    //__local float4 LATOMS[64];    // ToDo - LATER use local memory ?
    int i = ixyz(ix,iy,iz);
    //if(i==0){ printf( "GPU gradient() nxyz(%i,%i,%i) \n", nx,ny,nz ); }
    out[ i ] = (float4){
        (A[ ixyz(ix+1,iy  ,iz  ) ].x - A[ixyz(ix-1,iy  ,iz  )].x)*inv_d.x,
        (A[ ixyz(ix  ,iy+1,iz  ) ].x - A[ixyz(ix  ,iy-1,iz  )].x)*inv_d.y,
        (A[ ixyz(ix  ,iy  ,iz+1) ].x - A[ixyz(ix  ,iy  ,iz-1)].x)*inv_d.z,
        A[ i ].x
    }; 
    //if( (ix==60)&&(iy==32)&&(iz==15) ){ printf( "GPU out(%g,%g,%g,%g) iy(%i,%i) iz(%i,%i)\n", out[i].x,out[i].y,out[i].z,out[i].w,  ixyz(ix,iy+1,iz),ixyz(ix,iy-1,iz),   ixyz(ix,iy,iz+1),ixyz(ix,iy,iz-1)   ); }
};

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
    const float4 dCell     // 
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
    k = 1.0f-fabs(k-1.0f);  // 
    float  f = dCell.w/dot( k, k );    // dCell.w = 4*pi*eps0*dV - rescaling constant
    if(i==0)f=0;
    if(i<N){ 
        out[i] = A[i]*f;
        //out[i] = f;
        //out[i] = k.x;
    }
    //if( (ix==(nx/2))&&(iy==(ny/2)) ){ printf( "GPU iz,i[%i,%i|%i] k %g A[i] %g out %g \n", iz, i, N, f, A[i].x, out[i].x ); };
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
    
    if(iG>=nMax) return;
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
    
    if(iG>=nMax) return;
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




/**
 * @brief Kernel function to project electron density from sum of molecular orbitals to a grid using basis functions store in texture memory. 
 *
* @param[in] nAtoms   : The number of atoms
* @param[in] iorb0    : index of the first orbital to project
* @param[in] iorb1    : index of the last orbital to project
* @param[in] atoms    : array of atom positions and slot (x,y,z,islot), 'islot' is used to select the basis function from the texture 
* @param[in] coefs    : array of molecular orbital coefficients (px,py,pz,s), slot 's' is used to select the basis function from the texture
* @param[out] outGrid : output grid to store the projected density
* @param[in] imgIn    : input image containing the basis functions
* @param[in] nGrid    : size of the grid (nx,ny,nz,0)
* @param[in] grid_p0  : origin of the grid
* @param[in] grid_dA  : grid step along axis A
* @param[in] grid_dB  : grid step along axis B
* @param[in] grid_dC  : grid step along axis C
* @param[in] acumCoef : coefficients of the linear combination of the density (c0,c1), c0 is the coefficient of the old density stored in the grid, c1 is the coefficient of the new density to be added to the grid.
*/

__kernel void projectOrbDenToGrid_texture(
    const int nAtoms,            //1
    const int iorb0,             //2
    const int iorb1,             //3
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

    if(iG==0){  printf("GPU: projectOrbDenToGrid_texture acumCoef %g,%g nAtoms %i iorb(%i,%i) nL %i nMax %i  nGrid(%i|%i,%i,%i)\n", acumCoef.x, acumCoef.y, nAtoms, iorb0, iorb1, nL, nMax, nGrid.x*nGrid.y*nGrid.z, nGrid.x,nGrid.y,nGrid.z ); }

    
    if(iG>=nMax) return;
    
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

                    //float3 dp = pos-xyzq.xyz;
                    //float wfi = dot( cs, sp3_tex( dp*wf_tiles_per_angstroem,     xyzq.w, imgIn ) );
                    //if( (iorb==iorb0) && (i0==0) && ( dot(dp,dp)<0.1f ) ){ printf( "R=%f wfi=%f \n", dot(dp,dp), wfi ); };
                    //if( isnan(wfi) || (wfi>1e+3) || (wfi<-1e+3) ){ printf( "GPU wfi=%g orb=%i,atom=%i ig(%i,%i,%i)\n",wfi, iorb,i,  ia,ib,ic ); };
                    //wf.x += wfi; // accumulate atomic orbital contribution to molecular orbital
                }
            } // j
            barrier(CLK_LOCAL_MEM_FENCE);
        } // i0
        dens += wf.x*wf.x;   // accumulate the square of molecular orbital to the density
    } // iorb
    //if(iG==0){ printf( "GPU loop DONE ! \n" ); }
    //outGrid[iG] = (float2){dens,0.0f};
    if(fabs(acumCoef.x)<1e-8){                         // if c0==0 then we can just overwrite the grid (no density0 is used)
        outGrid[iG] = (float2){dens,0.0f}*acumCoef.y;    
    }else{                                             // if c0!=0 then we have to add the density to the grid rho_tot = c0*rho_0 + c1*rho_scf
        outGrid[iG] = outGrid[iG]*acumCoef.x + ((float2){dens,0.0f})*acumCoef.y;
    }
    //outGrid[iG] = (float2){1.2f,2.3f};
    //if(iG==0){ printf( "GPU all DONE ! \n" ); }
    //if(iG==0){ printf("projectOrbDenToGrid_texture END \n"); }    
}

__kernel void projectDenmatToGrid_simp(
    const int nAtoms,            //1
    __global int*     sel,       //2
    __global float4*  atoms,     //3
    __global float16* denmat,    //4
    __global float2*  outGrid,   //5
    __read_only image2d_t imgIn, //6 
    int4   nGrid,                //7
    float4 grid_p0,              //8
    float4 grid_dA,              //9
    float4 grid_dB,              //10
    float4 grid_dC,              //11
    float2 acumCoef              //12
){
    const int iG = get_global_id (0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;
    if(iG==0){  printf("GPU: projectDenmatToGrid_texture_simp acumCoef %g,%g nAtoms %i nMax %i  nGrid(%i|%i,%i,%i)\n", acumCoef.x, acumCoef.y, nAtoms,nMax, nGrid.x*nGrid.y*nGrid.z, nGrid.x,nGrid.y,nGrid.z ); }    
    if(iG>=nMax) return;
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;    
    // ToDo : Later we have to change the order of the loops
    float dens   = 0.0f;
    for (int i=0; i<nAtoms; i++ ){
        const int    ia    = sel[i]; 
        const float4 atomi = atoms [ia];
        const float4 wfi = sp3_tex( (pos-atomi.xyz)*wf_tiles_per_angstroem, atomi.w, imgIn );
        for (int j=0; j<nAtoms; j++ ){
            const int ja = sel[j];
            const float4  atomj  = atoms [ja];
            const float16 coefs  = denmat[ia*nAtoms+ja];
            const float4 wfj = sp3_tex( (pos-atomj.xyz)*wf_tiles_per_angstroem, atomj.w, imgIn );
            dens += wfi.x*dot( coefs.lo.lo, wfj )
                 +  wfi.y*dot( coefs.lo.hi, wfj )
                 +  wfi.z*dot( coefs.hi.lo, wfj )
                 +  wfi.w*dot( coefs.hi.hi, wfj );
        } // j
    } // i
    if(fabs(acumCoef.x)<1e-8){                         // if c0==0 then we can just overwrite the grid (no density0 is used)
        outGrid[iG] = (float2){dens,0.0f}*acumCoef.y;    
    }else{                                             // if c0!=0 then we have to add the density to the grid rho_tot = c0*rho_0 + c1*rho_scf
        outGrid[iG] = outGrid[iG]*acumCoef.x + ((float2){dens,0.0f})*acumCoef.y;
    }    
}

__kernel void projectDenmatToGrid(
    const int nAtoms,            //1
    __global int*     sel,       //2
    __global float4*  atoms,     //3
    __global float16* denmat,    //4
    __global float2*  outGrid,   //5
    __read_only image2d_t imgIn, //6 
    int4   nGrid,                //7
    float4 grid_p0,              //8
    float4 grid_dA,              //9
    float4 grid_dB,              //10
    float4 grid_dC,              //11
    float2 acumCoef              //12
){
    __local float4  LATOMS[N_LOCAL];
    __local float16 LCOEFS[N_LOCAL];
    const int iG = get_global_id (0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x; 
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab; 
    const int nMax = nab*nGrid.z;
    if(iG==0){  printf("GPU: projectDenmatToGrid_texture acumCoef %g,%g nAtoms %i iorb(%i,%i) nL %i nMax %i  nGrid(%i|%i,%i,%i)\n", acumCoef.x, acumCoef.y, nAtoms, nL, nMax, nGrid.x*nGrid.y*nGrid.z, nGrid.x,nGrid.y,nGrid.z ); }    
    if(iG>=nMax) return;
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;    
    // ToDo : Later we have to change the order of the loops
    float dens   = 0.0f;
    for (int i=0; i<nAtoms; i++ ){
        const int    ia    = sel[i]; 
        const float4 atomi = atoms [ia];
        const float4 wfi = sp3_tex( (pos-atomi.xyz)*wf_tiles_per_angstroem, atomi.w, imgIn );
        for (int j0=0; j0<nAtoms; j0+=nL ){
            const int j  = j0 + iL;
            const int ja = sel[j];
            LATOMS[iL] = atoms [ja];
            LCOEFS[iL] = denmat[ia*nAtoms+ja];
            barrier(CLK_LOCAL_MEM_FENCE);
            for (int j=0; j<nL; j++){
                if( (j+j0)<nAtoms ){ 
                    const float4  atomj  = LATOMS[j];
                    const float16 coefs  = LCOEFS[j];
                    const float4 wfj = sp3_tex( (pos-atomj.xyz)*wf_tiles_per_angstroem, atomj.w, imgIn );
                    dens += wfi.x*dot( coefs.lo.lo, wfj )
                         +  wfi.y*dot( coefs.lo.hi, wfj )
                         +  wfi.z*dot( coefs.hi.lo, wfj )
                         +  wfi.w*dot( coefs.hi.hi, wfj );
                }
            } // j
            barrier(CLK_LOCAL_MEM_FENCE);
        } // i0
    } // ia
    if(fabs(acumCoef.x)<1e-8){                         // if c0==0 then we can just overwrite the grid (no density0 is used)
        outGrid[iG] = (float2){dens,0.0f}*acumCoef.y;    
    }else{                                             // if c0!=0 then we have to add the density to the grid rho_tot = c0*rho_0 + c1*rho_scf
        outGrid[iG] = outGrid[iG]*acumCoef.x + ((float2){dens,0.0f})*acumCoef.y;
    }    
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
    if(iG>=nMax) return;
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
    //if(iG==0){ printf("projectAtomDenToGrid_texture END \n"); }
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



// ========================================================================= 
// ========================================================================= 
// ===================== Project Point Charges on Grid  ====================
// ========================================================================= 
// ========================================================================= 

// void make_inds_pbc(const int n, int4* indices) {
//     indices[0] = (int4){0, 1, 2, 3};
//     indices[1] = (int4){0, 1, 2, 3 - n};
//     indices[2] = (int4){0, 1, 2 - n, 3 - n};
//     indices[3] = (int4){0, 1 - n, 2 - n, 3 - n};
// }

inline int modulo(int i, int m) {
    int result = i % m;
    return (result < 0) ? (result + m) : result;
}


int4 make_inds_pbc(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int4)(0, 1,   2,   3  ); }
        case 1: { return (int4)(0, 1,   2,   3-n); }
        case 2: { return (int4)(0, 1,   2-n, 3-n); }
        case 3: { return (int4)(0, 1-n, 2-n, 3-n); }
    }
    return (int4)(-100, -100, -100, -100);
}

int8 make_inds_pbc_5(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int8){0, 1,   2,   3,   4,   5  ,-100,-100 }; }
        case 1: { return (int8){0, 1,   2,   3,   4,   5-n,-100,-100 }; }
        case 2: { return (int8){0, 1,   2,   3,   4-n, 5-n,-100,-100 }; }
        case 3: { return (int8){0, 1,   2,   3-n, 4-n, 5-n,-100,-100 }; }
        case 4: { return (int8){0, 1,   2-n, 3-n, 4-n, 5-n,-100,-100 }; }
        case 5: { return (int8){0, 1-n, 2-n, 3-n, 4-n, 5-n,-100,-100 }; }
    }
    return (int8)(-100,-100,-100,-100, -100,-100,-100,-100);
}

inline int4 choose_inds_pbc(const int i, const int n, const int4* iqs) {
    if (i >= (n-3)) {
        const int ii = i + 4 - n;
        return iqs[ii];
    }
    return (int4)(0, +1, +2, +3);
}

// inline int4 choose_inds_pbc_3( const int i, const int n, const int4* iqs ){
//     if(i>=(n-3)){ 
//         const int ii = i+4-n;
//         //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
//         const int4 d = iqs[ii];
//         return (int4){ i+d.x, i+d.y, i+d.z, i+d.w }; 
//     }
//     return (int4){ i, i+1, i+2, i+3 };
// }

inline int4 choose_inds_pbc_3( const int i, const int n, const int4* iqs ){
    if(i>=(n-3)){ 
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        const int4 d = iqs[ii];
        return (int4){ i+d.x, i+d.y, i+d.z, i+d.w }; 
    }
    return (int4){ i, i+1, i+2, i+3 };
}



int pbc_ifw(int i, int n){ i++; return (i<n )?  i :  i-n; };
int pbc_ibk(int i, int n){ i--; return (i>=0)?  i :  i+n; };

__kernel void laplace_real_pbc( 
    int4 ng,
    __global float* Vin, 
    __global float* Vout, 
    float cSOR
){

    const int ix = get_global_id(0);
    const int iy = get_global_id(1);
    const int iz = get_global_id(2);
    if( (ix>=ng.x) || (iy>=ng.y) || (iz>=ng.z) ) return;

    int nxy = ng.x * ng.y;
    const float fac = 1.0f/6.0f;
    //float err2 = 0.0; 
    
    const int iiz =          iz       *nxy;
    const int ifz =  pbc_ifw(iz, ng.z)*nxy;
    const int ibz =  pbc_ibk(iz, ng.z)*nxy;
    
    const int iiy =          iy       *ng.x;
    const int ify =  pbc_ifw(iy, ng.y)*ng.x;
    const int iby =  pbc_ibk(iy, ng.y)*ng.x;
    const int ifx =  pbc_ifw(ix, ng.x);
    const int ibx =  pbc_ibk(ix, ng.x);

    double vi = 
    Vin[ ibx + iiy + iiz ] + Vin[ ifx + iiy + iiz ] + 
    Vin[ ix  + iby + iiz ] + Vin[ ix  + ify + iiz ] + 
    Vin[ ix  + iiy + ibz ] + Vin[ ix  + iiy + ifz ];

    vi*=fac;
    const int i = ix + iiy + iiz;
    const float vo = Vin[ i ];
    vi += (vi-vo)*cSOR; 
    const float dv = vi - vo;
    
    //err2 += dv*dv;
    Vout[i] = vi;
}




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
    const int num_atoms,
    __global const float4* atoms,   // Atom positions and charges
    __global       float*  Qgrid,   // Output grid
    const int4 ng,                  // grid size
    const float3 g0,       // grid orgin
    const float3 dg    // grid dimensions
) {
    int iG = get_global_id(0);
    const int iL = get_local_id(0);
    if (iG >= num_atoms) return;

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
