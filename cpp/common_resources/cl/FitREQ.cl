#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f  // [eV*Ang/e^2]


__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond(
    const int4 ns,                  // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    const float4      GFFParams     // 9 // Grid-Force-Field parameters
){
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LCLJS [32];   // local buffer for atom parameters
    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory

    const int natoms=ns.x;  // number of atoms
    const int nnode =ns.y;  // number of node atoms
    //const int nAtomCeil =ns.w;
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)

    //const int i0n = iS*nnode; 
    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int i0v = iS*nvec;    // index of first atom in vectors array
    //const int ian = iG + i0n;
    const int iaa = iG + i0a; // index of atom in atoms array
    const int iav = iG + i0v; // index of atom in vectors array
    
    const float4 REQKi = REQKs    [iaa];  // non-bonded parameters
    const float3 posi  = atoms    [iav].xyz; // position of atom
    const float  R2damp = GFFParams.x*GFFParams.x; // squared damping radius
    float4 fe          = float4Zero;  // force on atom

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){      // loop over all atoms in the system, by chunks of size of local memory
        const int i=j0+iL;              // index of atom in local memory
        if(i<natoms){                   // j0*nL may be larger than natoms, so we need to check if we are not reading from invalid address
            LATOMS[iL] = atoms [i+i0v]; // read atom position to local memory 
            LCLJS [iL] = REQKs [i+i0a]; // read atom parameters to local memory
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( (ja!=iG) && (ja<natoms) ){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float4 aj = LATOMS[jl];    // read atom position   from local memory
                float4 REQK     = LCLJS [jl];    // read atom parameters from local memory
                float3 dp       = aj.xyz - posi; // vector between atoms
                REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy
                fe += getLJQH( dp, REQK, R2damp ); 
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    if(iG<natoms){

        forces[iav] = fe;           // If we do    run it as first forcefield, we can just store force (non need to clean it before in that case)

    }
}

__attribute__((reqd_work_group_size(32,1,1)))
__kernel void evalSampleDerivatives(
    const int4 ns,                  // (n1,n2,0,0) - number of atoms in fragments
    // Dynamical
    __global float4*  atoms,        // positions of atoms
    __global float4*  dEdREQs,      // output derivatives
    __global int*     hosts,        // index of host atoms (it electron pair), -1 if it is not electron pair
    // Parameters  
    __global float4*  REQHs,        // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int*     types,        // atom types
    const float4      GFFParams     // Grid-Force-Field parameters
){
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LREQKS[32];   // local buffer for REQ parameters

    const int iG = get_global_id  (0); // index of atom in fragment 1
    const int nG = get_global_size(0); // total number of atoms in fragment 1
    const int iL = get_local_id   (0); // local thread index
    const int nL = get_local_size (0); // workgroup size (32)

    const int n1 = ns.x;  // number of atoms in fragment 1
    const int n2 = ns.y;  // number of atoms in fragment 2
    
    // Load data for atom i from fragment 1
    float4 REQi = float4Zero;
    float3 posi = float3Zero;
    int    ti   = 0;
    float4 fREQi = float4Zero;
    
    if(iG >= n1) return;
    posi = atoms[iG].xyz;
    ti   = types[iG];
    REQi = REQHs[ti];
    int ih = hosts[iG];
    if( ih >= 0 ){
        int th = types[ih];
        REQi.z = REQHs[th].z - REQi.z;
    }


    // Process fragment 2 atoms in chunks of workgroup size
    for(int j0=0; j0<n2; j0+=nL){
        const int j = j0+iL;
        
        // Load chunk of fragment 2 atoms to local memory
        if(j < n2){
            const int jj = n1+j;
            LATOMS[iL] = atoms  [jj];
            int    ti   = types [jj];
            float4 REQj =  REQHs[ti];
            int ih = hosts[jj];
            if( ih >= 0 ){
                int th = types[ih];
                REQj.z = REQHs[th].z - REQj.z;
            }
            LREQKS[iL]  = REQj;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // Process atoms in local memory
        for(int jl=0; jl<nL; jl++){
            const int j = j0+jl;
            if(j < n2){
                // Load atom j data from local memory
                const float4 REQj  = LREQKS[jl];
                const float3 posj  = LATOMS[jl].xyz;


                // Compute LJQH2 derivatives
                float3 dij = posj - posi;
                float  r   = length(dij);
                float  ir  = 1.f/r;

                // Mixing rules
                float R0 = REQi.x + REQj.x;
                float E0 = REQi.y * REQj.y;
                float Q  = REQi.z * REQj.z;
                float H2 = REQi.w * REQj.w;
                float sH = (H2 < 0.f) ? 1.f : 0.f;
                H2 *= sH;

                // Electrostatic
                float dE_dQ  = ir * COULOMB_CONST;
                float Eel    = Q * dE_dQ;

                // Lennard-Jones
                float u   = R0*ir;
                float u3  = u*u*u;
                float u6  = u3*u3;
                float u6p = (1.f + H2)*u6;

                float dE_dE0 = u6 * (u6p - 2.f);
                float dE_dR0 = 12.f * (E0/R0) * u6 * (u6p - 1.f);
                float dE_dH2 = -E0 * u6 * u6;
                float ELJ    = E0 * dE_dE0;

                // Accumulate derivatives for atom i
                fREQi.x += dE_dR0;
                fREQi.y += dE_dE0 * 0.5f * REQj.y;
                fREQi.z += dE_dQ  *        REQj.z;
                fREQi.w += dE_dH2 *        REQj.w * sH;
                //atomic_add_float4(&dEdREQs[n1+j], fREQj);
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // Write accumulated derivatives for atom i
    // if(iG < n1){
    //     atomic_add_float4(&dEdREQs[iG], fREQi);
    // }


}

// Helper function for atomic float4 addition
inline void atomic_add_float4(__global float4* addr, float4 val){
    atomic_add(&addr->x, val.x);
    atomic_add(&addr->y, val.y); 
    atomic_add(&addr->z, val.z);
    atomic_add(&addr->w, val.w);
}
