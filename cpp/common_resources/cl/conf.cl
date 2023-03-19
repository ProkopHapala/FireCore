
#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}



/*
This is OpenCL version of   /common/molecular/Confs.h


F_i = d_( L - L0 ) / d_xi

L^2 = Sum_ig{   Kc*| cg1_ig - cg2_ig |^2  +    Ku* ( 1-cos(u1_ig,u2_ig) )  + Kv*( 1 - cos(u1_ig,u2_ig) }


cos(u1,u2) =  <u1|u2>/(|u1||u2|)


cg1 = Sum_i{ p_i * wc_i  } / Sum_i{ wc_i }

d{c1_ig}/d{p_i} =  wc_i / Sum_i{ wc_i  }

u1 = Sum_i{ p_i * wu_i  } - cg1 * Sum_i{wu_i}

d{u1}/d{x1i} =    wu_i -   wc_i  Sum_i{wu_i}/Sum_i{wc_i}   =   Wu_i     .... just constant (renormalized)
d{u1}/d{x1i} =    wu_i -   wc_i   ... if vec{wu} and vec{wc} are normalized 

d{ cos(u1,u2) }/dx1i = d{ <u1|u2>/(|u1|*|u2|) }/dx1i =     1/(|u1||u2|)  *  d{ <u1|u2> }/dx1i     +   <u1|u2>/|u2|  *  d{ 1/|u1| }/dx1i 
d{ <u1|u2> }/dx1i = x2i
d{ 1/|u1|  }/dx1i = -1/|u1|^3 * x1i
d{ cos(u1,u2) }/dx1_i =  x2_i /(|u1||u2|)  -  x1_i / <u1|u2>/( |u1|^3 * |u2| )

d{ cos(u1,u2) }/dx1_i =  x2_i * ilu1*ilu2  -  x1_i * cos(u1,u2) * ilu1*ilu1 =  Wui * ilu1*( x2_i*ilu2  -  x1_i*ilu1* cos(u1,u2) )   

=> need to store 
ilu1 = 1/|u1| 
ilu2 = 1/|u2|
ilv1 = 1s/|v1|
ilv2 = 1/|v2| 


d< Sum wi * (p1i - p2i) | Sum wi * (p1i - p2i) >/dpi
= wi * Sum wi * (p1i - p2i) >



*/


__kernel void conf_force_BAK(
    const int4 ns,                  // 1
    // Dynamical
    __global int2*    conf_pairs,   // 2  // pairs of configurations (offsets)
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global int*     atom2group,   // 3   // every atom belong to max 1 group ( -1 if none)
    // Parameters
    __global int2*    group_ns,    // 4  range of configuration (i0,n)
    __global int*     group_ia,    // 5  atom indexes belonging to the group
    __global float4*  group_cs,    // 6  projection coefficients for the group
){
    __local float4 ALs[32]; // 
    __local float4 GFs[32]; // distance for this group - should be buffer
    
    const int iG = get_group_id  (0);  // pairs of configurations   (independnet system, should not overlap !!!! )
    const int iL = get_local_id  (0);  // atoms in group

    int2  gi0s      = conf_pairs[iS];
    const int2  i0n = group_ns  [iL];

    //==== (( 1 ))  we need to sum atom in group => iL  

    float3 cg1 = float3Zero;
    float3 u1  = float3Zero;
    float3 v1  = float3Zero;
    float3 cg2 = float3Zero;
    float3 u2  = float3Zero;
    float3 v2  = float3Zero;
    for(int i=0; i<i0n.x; i++){     // loop over groups of this atom
        int ia = group_ia[i];  
        // system 1
        float3 p1 = atoms[ ia + gi0s.x ].xyz;
        cg1 += p1 * W.z; // W.z sum to 1
        u1  += p1 * W.x; // W.x sum to 1
        v1  += p1 * W.y; // W.y sum to 1
        
        // system 2
        float3 p2 = atoms[ ia + gi0s.y ].xyz;
        cg1 += p1 * W.z; // W.z sum to 1
        u1  += p1 * W.x; // W.x sum to 1
        v1  += p1 * W.y; // W.y sum to 1
        
    }
    float3 dcg = cg1-cg2;
    GC[iL]  = dcg;
    float ilu1 = 1.f/length(u1);
    float ilv1 = 1.f/length(v1);
    float ilu2 = 1.f/length(u2);
    float ilv2 = 1.f/length(v2);
    float cu12 = dot(u1,u2)*ilu1*ilu2;
    float cv12 = dot(v1,v2)*ilv1*ilv2;
    GU[iL] = (float3){ ilu1,ilu2,cu12 };
    GV[iL] = (float3){ ilv1,ilv2,cv12 };
    GL2[iL] =  K.z*dot(dcg,dcg)  +   K.x*(1-cos(u1,u2)) + K.y*(1-cos(v1,v2));
    
    //==== (( 2 ))  Sum total distance over all groups 

    barrier(CLK_LOCAL_MEM_FENCE);
    float l2 = 0;
    for(int ig=0; ig<ng; ig++ ){    // sum length over all the  groups
        l2 += GL2[ig];
    }
    float l   = sqrt(l2);
    float f_l = K*( l - l0 )/l; 
    // ------------- end reduction
    //barrier(CLK_LOCAL_MEM_FENCE);
    
    //==== (( 3 ))  apply force ... loop over atoms belonging to any of the groups

    float4 fe = float4Zero;

    int ig = atom2group[ia];

    float3 gu  = GU[ig];
    float3 gv  = GV[ig];
    float3 dcg = GV[ig];

    float4 fe1.xyz = (Wui * W.x* gu.x*( p2 *gu.y  -  p1*gu.x*gu.z )    // u-cos force 
                   +  Wvi * W.y* gv.x*( p2 *gv.y  -  p1*gv.x*gv.z )    // v-cos force
                   +  dcg*W.z) * f_l;                                  // cog force
    forces[ ia + i0s.x ] -= fe1;

    float4 fe2.xyz = (Wui * W.x* gu.y*( p1 *gu.x  -  p2*gu.y*gu.z )    // u-cos force 
                   +  Wvi * W.y* gv.y*( p1 *gv.x  -  p2*gv.y*gv.z )    // v-cos force
                   -  dcg*W.z) * f_l;                                  // cog force1
    forces[ ia + i0s.y ] += fe2;

}
























__kernel void conf_force_BAK(
    const int4 ns,                  // 1
    // Dynamical
    __global int2*    conf_pairs,   // 2  // pairs of configurations (offsets)
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global int*     atom2group,   // 3   // every atom belong to max 1 group ( -1 if none)
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