#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f  // [eV*Ang/e^2]

__attribute__((reqd_work_group_size(32,1,1)))
__kernel void evalSampleDerivatives(
    //const int4 ns,                  
    __global int4*    ranges,    // 1: [nsample]  (i0,ni, j0,nj) star and number of atoms in fragments 1,2 
    __global float4*  tREQHs,    // 2: [ntypes]   non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int*     atypes,    // 3: [natomTot] atom types
    //__global int*     hosts,   // [natomTot] index of host atoms (it electron pair), -1 if it is not electron pair
    __global int2*    ieps,      // 4: [natomTot] {iep1,iep2} index of electron pair type ( used to subtract  charge)
    __global float4*  atoms,     // 5: [natomTot] positions and REPS charge of each atom (x,y,z,Q) 
    __global float4*  dEdREQs,   // 6: [natomTot] output derivatives of type REQH parameters
    __global float2*  ErefW      // 7: [nsamp]   {E,W} reference energies and weights for each sample (molecule,system)
){
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LREQKS[32];   // local buffer for REQ parameters
    __local float  LdE;          // local variable to store energy difference

    const int iS = get_group_id   (0); // index of system
    const int iG = get_global_id  (0); // index of atom in fragment 1
    const int iL = get_local_id   (0); // local thread index
    const int nL = get_local_size (0); // workgroup size (32)

    const int4 nsi = ranges[iS];
    const int i0   = nsi.x;  // first     atom  in fragment 1
    const int ni   = nsi.z;  // number of atoms in fragment 1 
    const int j0   = nsi.y;  // first     atom  in fragment 2
    const int nj   = nsi.w;  // number of atoms in fragment 2
    
    if( iG-i0 >= ni) return;
    
    const int    ti    = atypes[iG];
    const float4 atomi = atoms [iG];
          float4 REQi  = tREQHs[ti];
    REQi.z = atomi.w;
    const int2   iep   = ieps  [iG];
    if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; } // subtract charge of electron pair 1
    if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; } // subtract charge of electron pair 2
    float4 fREQi = float4Zero;

    float Ei=0;
    // Process fragment 2 atoms in chunks of workgroup size
    for(int j0=0; j0<nj; j0+=nL){
        const int j = j0+iL;
        
        // Load chunk of fragment 2 atoms to local memory
        if(j<nj){
            const int jj       = j0+j;
            const int    tj    = atypes[jj];
            const float4 atomj = atoms [jj];
                  float4 REQj  = tREQHs[tj];
            REQj.z = atomj.w;
            const int2   jep  = ieps[jj];
            if( jep.x >= 0 ){ REQj.z -= tREQHs[jep.x].z; } // subtract charge of electron pair 1
            if( jep.y >= 0 ){ REQj.z -= tREQHs[jep.y].z; } // subtract charge of electron pair 2
            LATOMS[iL]  = atomj;
            LREQKS[iL]  = REQj;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        // Process atoms in local memory
        for(int jl=0; jl<nL; jl++){
            const int j = j0+jl;
            if(j < nj){
                // Load atom j data from local memory
                const float4 REQj  = LREQKS[jl];
                const float4 atomj = LATOMS[jl];

                // Compute LJQH2 derivatives
                float3 dij = atomj.xyz - atomi.xyz;
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
                float ELJ    =  E0 * dE_dE0;

                // Accumulate derivatives for atom i
                fREQi.x += dE_dR0;
                fREQi.y += dE_dE0 * 0.5f * REQj.y;
                fREQi.z += dE_dQ  *        REQj.z;
                fREQi.w += dE_dH2 *        REQj.w * sH;
                //atomic_add_float4(&dEdREQs[n1+j], fREQj);

                Ei += ELJ + Eel;
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    barrier(CLK_LOCAL_MEM_FENCE); // all local threads must finish loop above so we reuse  LATOMS

    LATOMS[iL].x = Ei;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL == 0) { 
        float Emol = 0.0f;
        for(int i=0; i<ni; i++){ Emol += LATOMS[i].x; }  // sum total energy from atomic
        float2 EW = ErefW[iS];
        LdE = (Emol - EW.x)*EW.y;
    } 
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL < ni) {
        dEdREQs[i0 + iL] = fREQi * LdE;
    }
}


#define NLOC_assembleDOFderivatives  128

__attribute__((reqd_work_group_size(NLOC_assembleDOFderivatives,1,1)))
// __kernel void assembleDOFderivatives(
//     int nDOFs,         
//     __global float*   fDOFs,      // [nDOFs]    derivatives of REQH parameters
//     __global int2*    DOFnis,     // [nDOFs]    (i0,ni) star and number of atoms in fragments 1,2    
//     __global int*     DOFtoAtom,  // [nInds]    list of atom indexes relevant for each DOF (non-uniform ranges indexed by DOFnis)
//     __global float4*  DOFcofefs,  // [nInds]    factors for update of each DOF from the atom dEdREQH parameters   fDOFi = dot( DOFcofefs[i], dEdREQH[i] )
//     __global float4*  dEdREQs     // [nAtomTot] derivatives of REQH parameters
// ){
//     __local float LfDOFi[NLOC_assembleDOFderivatives];

//     // this is reduction kernel, we want to use local memory to 
//     const int iG   = get_global_id  (0); // index of atom in fragment 1
//     const int iDOF = get_group_id   (0); // index of system
//     const int iL   = get_local_id   (0); // local thread index
//     const int nL   = get_local_size (0); // workgroup size (NLOC_assembleDOFderivatives)

//     const int2 nsi = DOFnis[iDOF];

//     // partial reduction to local memory - each 
//     float fDOFj = 0.0f;
//     for(int i0=0; i0<nsi.y; i0+=NLOC_assembleDOFderivatives){
//         int i = iL + i0;
//         if(i < nsi.y){
//             int j = DOFtoAtom[i + nsi.x];
//             fDOFj += dot( DOFcofefs[j], dEdREQs[j] );
//         }
//     } 
//     LfDOFi[iL] = fDOFj;
//     barrier(CLK_LOCAL_MEM_FENCE);

//     // final reduction from local memory
//     if(iL == 0){ 
//         float fDOFi = 0.0f;
//         for(int jl=0; jl<nL; jl++){  fDOFi += LfDOFi[jl];} 
//         fDOFs[iDOF] = fDOFi;
//     }

// }

__attribute__((reqd_work_group_size(NLOC_assembleDOFderivatives,1,1)))
__kernel void assembleAndRegularize(
     int nDOFs,                       // 1: number of DOFs
     __global float*   fDOFs,         // 2: [nDOFs]    derivatives of REQH parameters
     __global int2*    DOFnis,        // 3: [nDOFs]    (i0,ni) star and number of atoms in fragments 1,2    
     __global int*     DOFtoAtom,     // 4: [nInds]    list of atom indexes relevant for each DOF (non-uniform ranges indexed by DOFnis)
     __global float4*  DOFcofefs,     // 5: [nInds]    factors for update of each DOF from the atom dEdREQH parameters   fDOFi = dot( DOFcofefs[i], dEdREQH[i] )
     __global float4*  dEdREQs,       // 6: [nAtomTot] derivatives of REQH parameters
    __global const float*   DOFs,     // 7: [nDOFs]    current values of DOFs (need only for regularization)
    __global const float8*  regParams // 8: [nDOFs]   {min,max, xlo,xhi, Klo,Khi, K0,x0}
){
    __local float LfDOFi[NLOC_assembleDOFderivatives];

    const int iDOF = get_group_id(0);
    const int iL   = get_local_id(0);
    const int nL   = get_local_size(0);

    // --- Part 1: Assemble Physical Forces (same as before) ---
    const int2 nsi = DOFnis[iDOF];
    float fl = 0.0f;
    for(int i = iL; i < nsi.y; i += nL){
        int j   = i + nsi.x;
        int ia  = DOFtoAtom[j];
        fl     += dot(DOFcofefs[j], dEdREQs[ia]);
    }
    LfDOFi[iL] = fl;
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- Part 2: High-Performance Parallel Reduction ---
    for (int offset = nL / 2; offset > 0; offset /= 2) {
        if (iL < offset) { LfDOFi[iL] += LfDOFi[iL + offset]; }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- Part 2: Final Reduction & Regularization by Thread 0 ---
    if(iL == 0){
        //float fDOF = 0.0f; for(int i = 0; i < nL; i++){ fDOF += LfDOFi[i]; } // serial reduction - probably slow, but backup
        float fDOF = LfDOFi[0];
        // --- Regularization
        float  x = DOFs     [iDOF];
        float8 p = regParams[iDOF];          // p.s0=min, p.s1=max, p.s2=xlo, p.s3=xhi, p.s4=Klo, p.s5=Khi, p.s6=K0, p.s7=x0
        x = clamp(x, p.s0, p.s1);            // Clamp to hard limits (simple projection, more advanced methods exist)
        fDOF             -= p.s6*(x-p.s7);   // Spring force toward initial value K0*(x-x0)
        if(x<p.s2){ fDOF -= p.s4*(x-p.s2); } // Lower soft wall
        if(x>p.s3){ fDOF -= p.s5*(x-p.s3); } // Upper soft wall
        // --- Store
        fDOFs[iDOF] = fDOF;
    }
}