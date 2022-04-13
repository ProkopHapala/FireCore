#define R2SAFE          1e-4f
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

__constant sampler_t sampler_1 =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_2 =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;


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

float4 lerp_basis( float x, float slot, __read_only image2d_t imgIn ){
    float d = 0.01;
    float icoord;
    float fc     =  fract( x/d, &icoord );
    return read_imagef( imgIn, sampler_1, (float2){ x  , slot } )*(1.f-fc)
        +  read_imagef( imgIn, sampler_1, (float2){ x+d, slot } )*    fc ;
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
        printf( "DEBUG_GPU mul[%i] A=%g B=%g \n", i, A[i], B[i] );
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
    float4  fr  = lerp_basis( r, slot, imgIn );
    float4  v = (float4) ( dp*ir, 1 )*fr;
    return  dot( v ,c );
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
    //if(iG==0) for(int ia=0; ia<nAtoms;ia++)printf( "DEBUG_GPU atoms[%i](%g,%g,%g,%g) coefs[0](%g,%g,%g,%g) \n", ia, atoms[ia].x,atoms[ia].y,atoms[ia].z,atoms[ia].w,  coefs[ia].x,coefs[ia].y,coefs[ia].z,coefs[ia].w  );
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;
    outGrid[iG].x = sin(pos.x+pos.y+pos.z);
    float2 wf   = (float2) (0.0f,0.0f);
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
                wf.x += sp3_tex( (pos-xyzq.xyz)*0.3f, cs, 0.1, imgIn );

                //wf.x += sp3_tex( pos, cs, 0.1, imgIn );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    outGrid[iG] = wf;
}