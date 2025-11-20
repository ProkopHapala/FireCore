// Linearized Force-Field (LFF) OpenCL kernel
// Jacobi iterations with Chebyshev/heavy-ball mixing performed entirely on local data.

#define COULOMB_CONST   14.3996448915f       // [ eV*Ang/e^2 ]
#define _R2damp         1.e-8f

// evaluate non-covalent interaction force and energy for Lennard-Jones (Q) and Coulomb interactions of charges (Q) and hydrogen bond correction (pseudo-charges H), damping R2damp is used to avoid singularity at r=0
inline float4 getLJQH( float3 dp, float4 REQ, float R2damp ){
    // ---- Electrostatic (damped Coulomb potential)
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 +  R2damp);              // inverse distance squared and damped
    float   Ec    =  COULOMB_CONST*REQ.z*sqrt( ir2_ ); // Ec = Q1*Q2/sqrt(r^2+R2damp)
    // --- Lennard-Jones and Hydrogen bond correction
    float  ir2 = 1.f/r2;          // inverse distance squared
    float  u2  = REQ.x*REQ.x*ir2; // u2 = (R0/r)^2
    float  u6  = u2*u2*u2;        // u6 = (R0/r)^6
    float vdW  = u6*REQ.y;        // vdW = E0*(R0/r)^6
    float E    =       (u6-2.f)*vdW     + Ec  ;     // E = E0*(R0/r)^6 - E0*(R0/r)^12 + Q1*Q2/sqrt(r^2+R2damp)
    float fr   = -12.f*(u6-1.f)*vdW*ir2 - Ec*ir2_;  // fr = -12*E0*( (R0/r)^8/r + 12*E0*(R0/r)^14) - Q1*Q2/(r^2+R2damp)^1.5
    return  (float4){ dp*fr, E };
}


//#define LFF_MAX_VERTS   512
//#define LFF_MAX_VERTS   512

#define LFF_WG_SIZE    32
#define MAX_NEIGHBORS   8

__kernel void lff_jacobi(
    __global const int*    mols,       // 1  [nsys, nMol] {i0} starting index of first atom in each molecule in the bbuffer
    __global       float4* apos,       // 2  [nSys * nMol * natoms] positions  (xyz | charge )
    __global       float4* avel,       // 3  [nSys * nMol * natoms] velocities (xyz | mass   )
    __global const int*    neighs,     // 4  [nSys * nMol * natoms * LFF_MAX_NEIGHBORS] neighbor indices
    __global const float2* KLs,        // 5  [nSys * nMol * natoms * LFF_MAX_NEIGHBORS] {K,l0} bond params per neighbor, K is stiffness, l0 is rest length
    __global const int*    fixed_mask, // 6  optional 0/1 mask of fixed vertices
    const float3           Efield,     // 7  global external electric field
    const float            dt,         // 8  global time step
    const int              nOuter,     // 9  number of outer iterations
    const int              nInner,     // 10 number of inner iterations
    const float            bMix        // 11 momentum mixing parameter
){
    // NOTE: molecules are processed per workgroup in dimension 0; systems/replicas iterate over dimension 1
    const int imol = get_group_id(0);
    const int isys = get_group_id(1);
    const int lid  = get_local_id(0);
    const int lsz  = get_local_size(0);

    const int nMolPerSys  = get_num_groups(0);
    const int atomsPerSys = mols[nMolPerSys];
    const int sys_atom0   = isys * atomsPerSys;

    const int ia0    = sys_atom0 + mols[imol    ];
    const int ia1    = sys_atom0 + mols[imol + 1];
    const int natoms = ia1 - ia0;

    if (lid >= natoms) return;

    const int idx        = ia0 + lid;          // vertex index in global buffers
    const int neigh_base = idx * MAX_NEIGHBORS;

    __local float4 lpos[LFF_WG_SIZE];  // xyz: y, w: diag mass term
    float3 pi = apos[idx].xyz;
    float3 vi = avel[idx].xyz;
    float  mi = avel[idx].w;
    float  Qi = apos[idx].w;

    lpos[lid] = (float4){pi,mi};
    const float  inv_mass = 1.0f/mi;
    const float  inv_dt   = 1.0f/dt;
    const float  inv_dt2  = inv_dt*inv_dt;
    const float  Ii       = mi*inv_dt2;
    const int is_fixed = fixed_mask[idx] != 0;
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- debug if neighbor params are correct
    if(lid == 0){
        printf("GPU: lff_projective_jacobi() dt=%f bMix=%f nOuter=%i nInner=%i Efield=(%f,%f,%f)\n", dt, bMix, nOuter, nInner, Efield.x, Efield.y, Efield.z);
        for(int i=0; i<natoms; i++){
            int ing0 = i * MAX_NEIGHBORS;
            printf("GPU: atom %3i (%6.4f,%6.4f,%6.4f|%6.4f) neigs:", i, apos[i].x, apos[i].y, apos[i].z, apos[i].w) ;
            for (int jj = 0; jj < MAX_NEIGHBORS; jj++) {
                int j     = neighs[ing0 + jj];
                float2 kl = KLs   [ing0 + jj];
                if (j < 0) break;
                printf(" (%3i,%6.4f,%6.4f)", j, kl.y, kl.x);
            }
            printf("\n");
        }
    }

    // --- private copy of neighbor data
    int    ng_idx[MAX_NEIGHBORS];
    float2 ng_KLs[MAX_NEIGHBORS];
    int nneigh = 0;
    
    for (int jj = 0; jj < MAX_NEIGHBORS; ++jj) {
        const int j = neighs[neigh_base + jj];
        ng_idx[jj]  = j; // make sure j is local (in-molecule, in-group index)
        ng_KLs[jj]  = KLs [neigh_base + jj];
        if (j < 0) break;
        ++nneigh;
    }

    // --- Projective Dynamics MD-loop
    for (int outer = 0; outer < nOuter; ++outer) {

        if(is_fixed) continue;

        // --- external force
        float3 fex    = Efield * Qi;
        // TODO: non-covalent forces can be accumulated between outer iterations.


        // DEBUG - turn of dynamics
        //fex    = (float3)(0.0f, 0.0f, 0.0f);
        //vi     = (float3)(0.0f, 0.0f, 0.0f);

        // --- Leapfrog predictor
        float3 pi_old = pi;
        vi   += fex *  (dt * inv_mass);
        pi   += vi *dt;
        lpos[lid] = (float4){pi,mi};

        // ---- Jacobi-fly Linear Constraint Solver loop
        float3 mom_vec  = (float3)(0.0f, 0.0f, 0.0f);        
        for (int iter = 0; iter < nInner; ++iter) {
            if(lid==0)printf( "GPU iter %i\n", iter );
            barrier(CLK_LOCAL_MEM_FENCE);
            //const float  mi = lpos[lid].w;
            // inertial contribution: M_i/dt^2 * p_i
            //const float  Ii = mi*inv_dt2;
            float3  bi = pi*Ii;
            float  Aii = Ii;
            // ---- Loop over constraints (neighbors of given vertex)
            for (int jj=0; jj<nneigh; ++jj){
                int j      = ng_idx[jj];
                if (j < 0) break;
                // load parameters (float)
                float2  kl = ng_KLs[jj];
                float3 pj  = lpos[j - ia0].xyz; // use local index within molecule
                float3 dij = pi - pj;   
                float len  = length(dij);
                if(lid==0)printf( "GPU bond (%3i,%3i) l: %6.4f l0: %6.4f K: %6.4f\n", lid, j, len, kl.y, kl.x );
                float inv_len = len > 1e-8f ? 1.0f/len : 0.0f;
                float3 rest_pos = pj + dij * (kl.y * inv_len);
                bi  += rest_pos * kl.x;   // accumulate projected target
                Aii += kl.x;              // accumulate diagonal weight
            }
            const float   invA   = 1.0f / Aii;
            const float3  pi_new = bi * invA;

            // momentum acceleration
            // const float3  pi_    = pi_new + mom_vec * bMix;
            // mom_vec              = pi_ - pi;
            // pi                   = pi_;
            pi = pi_new;
            
            barrier(CLK_LOCAL_MEM_FENCE);
            lpos[lid]            = (float4)(pi, mi);
        }

        vi = (pi - pi_old)*inv_dt;
        
    }

    apos[idx] = (float4)(pi, Qi);
    avel[idx] = (float4)(vi, mi);
}




__kernel void lff_nb_jacobi(
    __global const int*    mols,       // 1  [nsys, nMol] {i0} starting index of first atom in each molecule in the bbuffer
    __global       float4* apos,       // 2  [nSys * nMol * natoms] positions  (xyz | charge )
    __global       float4* avel,       // 3  [nSys * nMol * natoms] velocities (xyz | mass   )
    __global       float4* REQHs,      // 4  [nSys * nMol * natoms] {RvdW,EvdW,Q,H}   non-covalent interaction parameters (radius | energy)
    __global const int*    neighs,     // 5  [nSys * nMol * natoms * LFF_MAX_NEIGHBORS] neighbor indices
    __global const float2* KLs,        // 6  [nSys * nMol * natoms * LFF_MAX_NEIGHBORS] {K,l0} bond params per neighbor, K is stiffness, l0 is rest length
    __global const int*    fixed_mask, // 7  optional 0/1 mask of fixed vertices
    const float3           Efield,     // 8  global external electric field
    const float            dt,         // 9  global time step
    const int              nOuter,     // 10  number of outer iterations
    const int              nInner,     // 11 number of inner iterations
    const float            bMix        // 12 momentum mixing parameter
){

    // NOTE: we have muliple molecules in the same system each molecule processed by one workgroup, the molecules interact just non-covalently, 
    //    we also paralelize over the system replicas but there is no interaction between replicas
    const int imol = get_group_id(0);
    const int isys = get_group_id(1);
    const int lid  = get_local_id(0);
    const int lsz  = get_local_size(0);

    const int nMolPerSys  = get_num_groups(0);
    const int atomsPerSys = mols[nMolPerSys];
    const int sys_atom0   = isys * atomsPerSys;

    const int ia0    = sys_atom0 + mols[imol    ];
    const int ia1    = sys_atom0 + mols[imol + 1];
    const int natoms = ia1 - ia0;

    if (lid >= natoms) return;

    const int idx        = ia0 + lid;          // vertex index in global buffers
    const int neigh_base = idx * MAX_NEIGHBORS;

    __local float4 lpos[LFF_WG_SIZE];  // xyz: y, w: diag mass term
    __local float4 lREQ[LFF_WG_SIZE];  // scratch for streamed REQ parameters
    const float4 REQKi    = REQHs[idx];

    float3 pi = apos[idx].xyz;
    float3 vi = avel[idx].xyz;
    float  mi = avel[idx].w;
    float  Qi = apos[idx].w;

    lpos[lid] = (float4){pi,mi};
    const float  inv_mass = 1.0f/mi;
    const float  inv_dt   = 1.0f/dt;
    const float  inv_dt2  = inv_dt*inv_dt;
    const float  Ii       = mi*inv_dt2;
    const int is_fixed = fixed_mask[idx] != 0;
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- debug if neighbor params are correct
    if(lid == 0){
        printf("GPU: lff_projective_jacobi() dt=%f bMix=%f nOuter=%i nInner=%i Efield=(%f,%f,%f)\n", dt, bMix, nOuter, nInner, Efield.x, Efield.y, Efield.z);
        for(int i=0; i<natoms; i++){
            const int ig = ia0 + i;
            int ing0 = ig * MAX_NEIGHBORS;
            printf("GPU: atom %3i (%6.4f,%6.4f,%6.4f|%6.4f) neigs:", ig, apos[ig].x, apos[ig].y, apos[ig].z, apos[ig].w) ;
            for (int jj = 0; jj < MAX_NEIGHBORS; jj++) {
                int j     = neighs[ing0 + jj];
                float2 kl = KLs   [ing0 + jj];
                if (j < 0) break;
                printf(" (%3i,%6.4f,%6.4f)", j, kl.y, kl.x);
            }
            printf("\n");
        }
    }

    // --- private copy of neighbor data
    int    ng_idx[MAX_NEIGHBORS];
    float2 ng_KLs[MAX_NEIGHBORS];
    int nneigh = 0;
    
    for (int jj = 0; jj < MAX_NEIGHBORS; ++jj) {
        const int j = neighs[neigh_base + jj];
        ng_idx[jj]  = j; // make sure j is local (in-molecule, in-group index)
        ng_KLs[jj]  = KLs [neigh_base + jj];
        if (j < 0) break;
        ++nneigh;
    }

    // --- Projective Dynamics MD-loop
    for (int outer = 0; outer < nOuter; ++outer) {
        // --- external force
        float4 fe    = (float4){Efield * Qi, 0.0f};
        // TODO: non-covalent forces can be accumulated between outer iterations.

        // --- Non-Covalent interactions - we are re-using local memory for temporary storage
        for (int j0=0; j0<atomsPerSys; j0+=LFF_WG_SIZE){     // loop over all atoms in the system, by chunks of size of local memory
            const int chunkStart = sys_atom0 + j0;
            const int jLoad = chunkStart + lid;            // index of atom in global buffers for this chunk
            if(jLoad < sys_atom0 + atomsPerSys){                  // guard against reading beyond current system block
                lpos[lid] = apos [jLoad];  // read atom position to local memory 
                lREQ[lid] = REQHs[jLoad];  // read atom parameters to local memory
            }else{ // just for debug - can be remove later
                lpos[lid] = (float4)(0.0f,0.0f,0.0f,0.0f);
                lREQ[lid] = (float4)(0.0f,0.0f,0.0f,0.0f);
            }
            barrier(CLK_LOCAL_MEM_FENCE);        // wait until all atoms are read to local memory
            const int chunkCount = (j0 + LFF_WG_SIZE <= atomsPerSys) ? LFF_WG_SIZE : (atomsPerSys - j0);
            for (int jl=0; jl<chunkCount; jl++){         // loop over all atoms in local memory (like 32 atoms)
                const int ja = chunkStart + jl;              // index of atom in global memory
                if( ja == idx ) continue; // if atom is not the same as current atom
                if((ja>=ia0)&&(ja<ia1)){  // exclusion of non-covalent interactions with bonded atoms (i.e. within the same molecule)
                    // NOTE: we do not need to care about exclusions of angles because angles are already represented by bonds
                    bool bonded=false;
                    for(int jj=0;jj<nneigh;jj++){ if(ng_idx[jj]==ja){bonded=true;break;} }
                    if(bonded) continue;
                }
                const float4 aj = lpos[jl];  // read atom position   from local memory
                float4 REQK     = lREQ[jl];  // read atom parameters from local memory
                float3 dp       = aj.xyz - pi; // vector between atoms
                REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy
                float4 fij = getLJQH( dp, REQK, _R2damp ); 
                fe += fij;
            }
            //barrier(CLK_LOCAL_MEM_FENCE);
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        
        
        lpos[lid]=(float4){pi,mi};  // recover lpos
        barrier(CLK_LOCAL_MEM_FENCE);

        // --- Leapfrog predictor
        float3 pi_old = pi;
        vi   += fe.xyz *  (dt * inv_mass);
        pi   += vi *dt;
        lpos[lid] = (float4){pi,mi};

        // ---- Jacobi-fly Linear Constraint Solver loop
        float3 mom_vec  = (float3)(0.0f, 0.0f, 0.0f);        
        for (int iter = 0; iter < nInner; ++iter) {
            if(lid==0)printf( "GPU iter %i\n", iter );
            barrier(CLK_LOCAL_MEM_FENCE);
            //const float  mi = lpos[lid].w;
            // inertial contribution: M_i/dt^2 * p_i
            //const float  Ii = mi*inv_dt2;
            float3  bi = pi*Ii;
            float  Aii = Ii;
            // ---- Loop over constraints (neighbors of given vertex)
            for (int jj=0; jj<nneigh; ++jj){
                int j      = ng_idx[jj];
                if (j < 0) break;
                // load parameters (float)
                float2  kl = ng_KLs[jj];
                float3 pj  = lpos[j].xyz; // make sure j is local (in-molecule, in-group index)
                float3 dij = pi - pj;   
                float len  = length(dij);
                if(lid==0)printf( "GPU bond (%3i,%3i) l: %6.4f l0: %6.4f K: %6.4f\n", lid, j, len, kl.y, kl.x );
                float inv_len = len > 1e-8f ? 1.0f/len : 0.0f;
                float3 rest_pos = pj + dij * (kl.y * inv_len);
                bi  += rest_pos * kl.x;   // accumulate projected target
                Aii += kl.x;              // accumulate diagonal weight
            }
            const float   invA   = 1.0f / Aii;
            const float3  pi_new = bi * invA;

            // momentum acceleration
            // const float3  pi_    = pi_new + mom_vec * bMix;
            // mom_vec              = pi_ - pi;
            // pi                   = pi_;
            pi = pi_new;
            
            barrier(CLK_LOCAL_MEM_FENCE);
            lpos[lid]            = (float4)(pi, mi);
        }

        vi = (pi - pi_old)*inv_dt;
        
    }

    apos[idx] = (float4)(pi, Qi);
    avel[idx] = (float4)(vi, mi);
}