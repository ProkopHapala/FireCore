#define R2SAFE          1e-4f
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

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