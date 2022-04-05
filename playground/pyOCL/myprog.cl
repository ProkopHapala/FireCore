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
    float4  v = (float4) ( 1, dp*ir );
    //v*=fr;
    return  dot( v ,c )*fr;
}


#define N_LOCAL  32

__kernel void projectAtomsToGrid(
    const int nAtoms, 
    __global float4* atoms,
    __global float4* coefs,
    //__global float4*    FE,
    __global float2*     outGrid,
    int4 nGrid,
    float4 grid_p0,
    float4 grid_dA,
    float4 grid_dB,
    float4 grid_dC
    //float  Rcut2,
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
    //if(iG==0) printf( " Qs (%g,%g,%g,%g) QZs (%g,%g,%g,%g) \n", Qs.x,Qs.y,Qs.z,Qs.w,   QZs.x,QZs.y,QZs.z,QZs.w   );
    float3 pos  = grid_p0.xyz + grid_dA.xyz*ia + grid_dB.xyz*ib  + grid_dC.xyz*ic;
    float2 wf   = (float2) (0.0f,0.0f);
    for (int i0=0; i0<nAtoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
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
                wf.x += sp3( xyzq.xyz, cs, xyzq.w );
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //if ( (ia==75)&&(ib==75) ) { printf(" iz %i fe %g,%g,%g,%g \n", ic, fe.x, fe.y, fe.z, fe.w ); }
    outGrid[iG] = wf;
}