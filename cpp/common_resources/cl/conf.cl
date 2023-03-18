
#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

__kernel void conf_force(
    const int4 ns,                  // 1
    // Dynamical
    __global int2*    conf_pairs,   // 2  // pairs of configurations (offsets)
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global int4*    atom2group,   // 3   // every atom belong to max 4 groups
    // Parameters
    __global int2*    group_ns,    // 4  range of configuration (i0,n)
    __global int*     group_ia,    // 5  atom indexes belonging to the group
    __global float4*  group_cs,    // 6  projection coefficients for the group
){
    __local float4 ALs[32]; // 
    __local float4 GFs[32]; // distance for this group - should be buffer
    
    const int iG = get_group_id  (0);  // pairs of configurations   (independnet system, should not overlap !!!! )
    const int iL = get_local_id  (0);  // atoms in group

    int2  i0s      = conf_pairs[iG];

    int ia = iL;
    const int4 a2g = atom2group[ia];
    const int* igs = (int*)&a2g;

    //==== (( 1 ))   set(AL2s)  ... first we store partial distances over atoms in each group  AL2s ( each atom and each group has separate bucket)       loop[ nAtom, nGroupPerAtom ]  

    for(int i=0; i<4; i++){     // loop over groups of this atom
        int ig = igs[i];        // index of group
        if(ig<0) continue;

        float4 L2;  // partial square distance for this atom 
        
        float3 p1 = atoms[ ia + i0s.x ].xyz;
        float3 p2 = atoms[ ia + i0s.y ].xyz;
        // ----- distnace
        float3 dp =  p1 - p2;
        L2.x  = dot( dp,dp );
        // ----- Orientation
        
        ALs[ia] = L2;
    }

    //==== (( 2 )) GFs=F(Sum( AL2s)) ...  then we sum ALs over atoms and define GFs for each groups   ( loop[ nGroup, nAtomOfGroup ] )

    // ------------- reduction here
    barrier(CLK_LOCAL_MEM_FENCE);

    // ------------- end reduction
    barrier(CLK_LOCAL_MEM_FENCE);
    
    //==== (( 3 )) apply(GFs)       ...  then we apply GFs back to the atoms   ( loop[ nAtom, nGroupPerAtom ] )   

    float4 fe = float4Zero;

    for(int i=0; i<4; i++){     // loop over groups
        int ig = igs[i];        // index of group
        if(ig<0) continue;
        float4 GF = GFs[ig];
        // ---- atom distance force
        float fr  = GF.x;    // group force
        fe.xyz   += dp * fr;
        // --- direction cosine force
        float fu  = GF.y;

        float fv  = GF.z;
    }

    forces[ ia + i0s.x ] -= fe;
    forces[ ia + i0s.y ] += fe;

}