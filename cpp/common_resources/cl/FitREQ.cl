#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f  // [eV*Ang/e^2]

#define iDBG 0

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
    __local float  LdE;          // shared energy diff across workgroup          // local variable to store energy difference

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

    if( iG == iDBG ){
        printf("GPU: OCL evalSampleDerivatives() iG %i i0,ni %3i,%3i j0,nj %3i,%3i ErefW %g %g\n", iG, i0,ni, j0,nj, ErefW[iS].x, ErefW[iS].y);
        for(int i=0; i<ni; i++){ int ia=i0+i; int it=atypes[ia]; printf("GPU: atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
        for(int i=0; i<nj; i++){ int ia=j0+i; int it=atypes[ia]; printf("GPU: atom j %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    }
    
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
    for(int off=0; off<nj; off+=nL){
        const int local_j = off + iL;
        
        // Load chunk of fragment 2 atoms to local memory
        if(local_j<nj){
            const int jj       = j0 + local_j;
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
            const int local_j2 = off + jl;
            const int jj2 = j0 + local_j2;
            if(local_j2 < nj){
                // Load atom j data from local memory
                const float4 REQj  = LREQKS[jl];
                const float4 atomj = LATOMS[jl];

                //if( iG == iDBG ){  printf("GPU: j %2i  pi( %16.8f %16.8f %16.8f ) pj( %16.8f %16.8f %16.8f )\n", jl, atomi.x, atomi.y, atomi.z, atomj.x, atomj.y, atomj.z); }
                //printf("GPU: iG %2i iS %2i iL %2i j %2i pi( %16.8f %16.8f %16.8f ) pj( %16.8f %16.8f %16.8f )\n", iG, iS, iL, jj2, atomi.x, atomi.y, atomi.z, atomi.w, atomj.x, atomj.y, atomj.z, atomj.w); 

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

                

                //if( iG == iDBG ){  printf("GPU: j %2i ELJ %16.8e Eel %16.8e r %16.8e    \n", jl, ELJ, Eel, r); }
                printf("GPU: iG %2i iS %2i iL %2i j %2i ELJ %16.8e Eel %16.8e r %16.8e  R0 %16.8e E0 %16.8e Q %16.8e H2 %16.8e     dE_dR0 %16.8e dE_dE0 %16.8e dE_dQ %16.8e dE_dH2 %16.8e \n", iG, iS, iL, jl, ELJ, Eel, r,  R0, E0, Q, H2, dE_dR0, dE_dE0, dE_dQ, dE_dH2); 
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    printf("GPU: iG %2i iS %2i iL %2i Ei %16.8e fREQi( %16.8e %16.8e %16.8e %16.8e )\n", iG, iS, iL, Ei, fREQi.x, fREQi.y, fREQi.z, fREQi.w);
    barrier(CLK_LOCAL_MEM_FENCE); // all local threads must finish loop above so we reuse  LATOMS

    LATOMS[iL].x = Ei;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL == 0) { 
        float Emol = 0.0f;
        for(int i=0; i<ni; i++){ Emol += LATOMS[i].x; }  // sum total energy from atomic
        float2 EW = ErefW[iS];
        LdE = (Emol - EW.x)*EW.y;  // thread 0 writes shared value
        if( iG == iDBG ){  printf("GPU: Emol %16.8e Eref %16.8e LdE %16.8e\n", Emol, EW.x, LdE); }
    } 
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL < ni) {
        //if( iG == iDBG ){  printf("GPU: fREQi %16.8e LdE %16.8e\n", fREQi, LdE); }
        printf("GPU: iG %2i iS %2i iL %2i dE %16.8e fREQi( %16.8e %16.8e %16.8e %16.8e )\n", iG, iS, iL, LdE, fREQi.x, fREQi.y, fREQi.z, fREQi.w);
        dEdREQs[i0 + iL] = fREQi * LdE;
    }
}

// -----------------------------------------------------------------------------
// Template variant of evalSampleDerivatives.
// The pairwise model code is injected at the marker line:
//     //<<<MODEL_PAIR_ACCUMULATION
// from Python via OpenCLBase.preprocess_opencl_source(macros={"MODEL_PAIR_ACCUMULATION": code}).
// The injected code runs inside the i–j pair loop and is expected to update:
//   - fREQi (float4): per-atom-i derivatives (dE/dR_i, dE/dE_i, dE/dQ_i, dE/dH_i)
//   - Ei    (float ): accumulated pair energy for atom i
// using the following in-scope variables:
//   atomi, atomj (float4), REQi, REQj (float4), dij (float3), r (float), ir (float)
// -----------------------------------------------------------------------------
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void evalSampleDerivatives_template(
    __global int4*    ranges,
    __global float4*  tREQHs,
    __global int*     atypes,
    __global int2*    ieps,
    __global float4*  atoms,
    __global float4*  dEdREQs,
    __global float2*  ErefW,
    float4   globParams  // {alpha,?,?,?} global parameters (min,max, xlo,xhi, Klo,Khi, K0,x0)
){
    __local float4 LATOMS[32];
    __local float4 LREQKS[32];
    __local float  LdE;

    const int iS = get_group_id   (0);
    const int iG = get_global_id  (0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int4 nsi = ranges[iS];
    const int i0   = nsi.x;
    const int ni   = nsi.z;
    const int j0   = nsi.y;
    const int nj   = nsi.w;

    if( iG - i0 >= ni ) return;

    const int    ti    = atypes[iG];
    const float4 atomi = atoms [iG];
          float4 REQi  = tREQHs[ti];
    REQi.z = atomi.w;
    const int2   iep   = ieps  [iG];
    if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
    if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }

    float4 fREQi = float4Zero;
    float  Ei    = 0.0f;

    for(int off=0; off<nj; off+=nL){
        const int local_j = off + iL;
        if(local_j<nj){
            const int jj       = j0 + local_j;
            const int    tj    = atypes[jj];
            const float4 atomj = atoms [jj];
                  float4 REQj  = tREQHs[tj];
            REQj.z = atomj.w;
            const int2   jep  = ieps[jj];
            if( jep.x >= 0 ){ REQj.z -= tREQHs[jep.x].z; }
            if( jep.y >= 0 ){ REQj.z -= tREQHs[jep.y].z; }
            LATOMS[iL]  = atomj;
            LREQKS[iL]  = REQj;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        for(int jl=0; jl<nL; jl++){
            const int local_j2 = off + jl;
            if(local_j2 < nj){
                const float4 atomj = LATOMS[jl];
                const float4 REQj  = LREQKS[jl];

                float3 dij  = atomj.xyz - atomi.xyz;
                float r    = length(dij);
                float inv_r = 1.f/fmax(r, R2SAFE);
                float R0 = REQi.x + REQj.x;
                float E0 = REQi.y * REQj.y;
                float Q  = REQi.z * REQj.z;
                float H  = REQi.w * REQj.w;
                // Mixing rules and hydrogen-bond gating
                float sH = (H < 0.f) ? 1.f : 0.f; // only apply H2 when negative
                H *= sH;

                //<<<MODEL_PAIR_ACCUMULATION
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    LATOMS[iL].x = Ei;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL == 0) {
        float Emol = 0.0f;
        for(int i=0; i<ni; i++){ Emol += LATOMS[i].x; }
        float2 EW = ErefW[iS];
        LdE = (Emol - EW.x)*EW.y;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL < ni) {
        dEdREQs[i0 + iL] = fREQi * LdE;
    }
}

// -----------------------------------------------------------------------------
// Energy-only template variant.
// Inject model-specific pair energy at the marker:
//     //<<<MODEL_PAIR_ENERGY
// The injected code should only accumulate into Ei using in-scope variables:
//   atomi, atomj (float4), REQi, REQj (float4), dij (float3), r (float), ir (float)
// and should NOT reference derivatives.
// Produces one scalar energy per sample (work-group) in Emols[iS].
// -----------------------------------------------------------------------------
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void evalSampleEnergy_template(
    __global int4*    ranges,    // (i0, j0, ni, nj) but stored as (x=i0, y=j0, z=ni, w=nj)
    __global float4*  tREQHs,    // [ntypes] REQH parameters
    __global int*     atypes,    // [nAtomTot]
    __global int2*    ieps,      // [nAtomTot]
    __global float4*  atoms,     // [nAtomTot]
    __global float*   Emols,      // [nSamples] output molecular energies
    float4   globParams  // {alpha,?,?,?} global parameters (min,max, xlo,xhi, Klo,Khi, K0,x0)
){
    __local float4 LATOMS[32];
    __local float4 LREQKS[32];

    const int iS = get_group_id   (0);
    const int iG = get_global_id  (0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int4 nsi = ranges[iS];
    const int i0   = nsi.x;
    const int ni   = nsi.z;
    const int j0   = nsi.y;
    const int nj   = nsi.w;

    // Each work-item handles one atom i within fragment 1 of the current sample.
    // Use local index within the work-group to cover [i0, i0+ni). Avoid early returns to keep barriers safe.
    const int    i     = i0 + iL;
    const int    active = (iL < ni);
    int    ti    = 0;
    float4 atomi = (float4)(0.0f,0.0f,0.0f,0.0f);
    float4 REQi  = (float4)(0.0f,0.0f,0.0f,0.0f);
    if(active){
        ti    = atypes[i];
        atomi = atoms [i];
        REQi  = tREQHs[ti];
        REQi.z = atomi.w;
        const int2   iep   = ieps  [i];
        if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
        if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }
    }

    float  Ei    = 0.0f;

    float alpha = globParams.x;

    // --- Debug: print config for a chosen sample and few lanes ---
    if((iS==iDBG) && (iL==0)){ 
        printf("GPU: evalSampleEnergy_template() nG %7i nL %2i nS %6i | i0=%d ni=%d j0=%d nj=%d | i=%d ti=%d\n", get_global_size(0), get_local_size(0), get_num_groups(0), i0, ni, j0, nj, i, ti); 
        for(int i=0; i<ni; i++){
            int ia=i0+i;
            int ti=atypes[ia];
            float4 atomi = atoms [ia];
            float4 REQi  = tREQHs[ti];
            REQi.z = atomi.w;
            const int2   iep   = ieps  [ia];
            if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
            if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }
            printf("GPU: frag1 atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, ti, atomi.x, atomi.y, atomi.z, atomi.w, REQi.x, REQi.y, REQi.z, REQi.w);
        }
        for(int i=0; i<nj; i++){
            int ia=j0+i;
            int ti=atypes[ia];
            float4 atomi = atoms [ia];
            float4 REQi  = tREQHs[ti];
            REQi.z = atomi.w;
            const int2   iep   = ieps  [ia];
            if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
            if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }
            printf("GPU: frag2 atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, ti, atomi.x, atomi.y, atomi.z, atomi.w, REQi.x, REQi.y, REQi.z, REQi.w);
        }
    }

    for(int off=0; off<nj; off+=nL){
        const int local_j = off + iL;
        if(local_j<nj){
            const int jj       = j0 + local_j;
            const int    tj    = atypes[jj];
            const float4 atomj = atoms [jj];
                  float4 REQj  = tREQHs[tj];
            REQj.z = atomj.w;
            const int2   jep  = ieps[jj];
            if( jep.x >= 0 ){ REQj.z -= tREQHs[jep.x].z; }
            if( jep.y >= 0 ){ REQj.z -= tREQHs[jep.y].z; }
            LATOMS[iL]  = atomj;
            LREQKS[iL]  = REQj;
        }
        barrier(CLK_LOCAL_MEM_FENCE);

        for(int jl=0; jl<nL; jl++){
            const int local_j2 = off + jl;
            if(local_j2 < nj){
                const float4 atomj = LATOMS[jl];
                const float4 REQj  = LREQKS[jl];

                float3 dij  = atomj.xyz - atomi.xyz;
                float r     = length(dij);
                float inv_r = 1.f/fmax(r, R2SAFE);

                float R0 = REQi.x + REQj.x;
                float E0 = REQi.y * REQj.y;
                float Q  = REQi.z * REQj.z;
                float H  = REQi.w * REQj.w;
                // Mixing rules and hydrogen-bond gating
                float sH = (H < 0.f) ? 1.f : 0.f; // only apply H2 when negative
                H *= sH;

                if((iS==iDBG) && (iL==0)){ 
                    printf("GPU: evalSampleEnergy_template() ia,ja (%3i,%3i) R0 %16.8f E0 %16.8f Q %16.8f H %16.8f\n", i, jl+j0, R0, E0, Q, H);
                }

                if(active){
                    //<<<MODEL_PAIR_ENERGY
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    LATOMS[iL].x = active ? Ei : 0.0f;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL == 0) {
        float Emol = 0.0f;
        for(int i=0; i<ni; i++){ Emol += LATOMS[i].x; }
        Emols[iS] = Emol;
        //if(iS == iDBG){ printf("GPU: iS=%d Emol=%g\n", iS, Emol); }
    }
}

// -----------------------------------------------------------------------------
// Rest of the code remains the same
// -----------------------------------------------------------------------------

#define NLOC_assembleDOFderivatives  128

__attribute__((reqd_work_group_size(NLOC_assembleDOFderivatives,1,1)))
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