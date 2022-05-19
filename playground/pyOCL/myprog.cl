#define R2SAFE          1e-4f
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

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
    if(i<N){ 
        //printf( "DEBUG_GPU mul[%i] A=%g B=%g \n", i, A[i], B[i] );
        out[i] = A[i] * B[i]; 
        //out[i] = sin( i*0.1 ); 
    }
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


float sp3_tex( float3 dp, float4 c, float slot, __read_only image2d_t imgIn ){
    float   r2  = ( dot(dp,dp) +  R2SAFE );
    float   r   = sqrt(r2);
    float   ir  = 1/r;
    float4  fr  = lerp_basis( r, slot, imgIn ).xyyy;
    //float4  fr  = read_imagef( imgIn, sampler_wrf, (float2){ r, slot } );
    float4  v = (float4) ( dp*ir*pref_p, pref_s )*fr;
    //float4  v = (float4) ( dp*ir, 0.0f )*fr;
    return  dot( v ,c );
    //return r;
}



#define N_LOCAL  32

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
    if(iG==0) for(int ia=0; ia<nAtoms;ia++)printf( "DEBUG_GPU atoms[%i](%g,%g,%g,%g) coefs[0](%g,%g,%g,%g) \n", ia, atoms[ia].x,atoms[ia].y,atoms[ia].z,atoms[ia].w,  coefs[ia].x,coefs[ia].y,coefs[ia].z,coefs[ia].w  );
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
                wf.x += sp3_tex( (pos-xyzq.xyz)*wf_tiles_per_angstroem, cs, xyzq.w, imgIn );
                //if((ia==16)&&(ib==16)){ printf("DEBUG_GPU pos(%g,%g,%g) xyzq(%g,%g,%g,%g) cs(%g,%g,%g,%g) wf(%g,%g) \n", pos.x,pos.y,pos.z,  xyzq.x,xyzq.y,xyzq.z,xyzq.w, cs.x,cs.y,cs.z,cs.w, wf.x, wf.y ); }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //if((ia==16)&&(ib==16)){ printf("DEBUG_GPU pos(%g,%g,%g) wf(%g,%g) \n", pos.x,pos.y,pos.z, wf.x, wf.y ); }
    outGrid[iG] = wf;
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
                wf.x += sp3_tex( dp, cs, xyzq.w, imgIn );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    out[iG] = wf;
}