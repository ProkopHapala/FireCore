#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f  // [eV*Ang/e^2]

#define iDBG 0
//<<<HBOND_GATE_DEFINE
#ifndef HBOND_GATE
#define HBOND_GATE 1
#endif
#if HBOND_GATE
#define APPLY_H_GATE(Hval) ((Hval) < 0.f ? (Hval) : 0.f)
#define S_H(Hval)          ((Hval) < 0.f ? 1.f     : 0.f)
#else
#define APPLY_H_GATE(Hval) (Hval)
#define S_H(Hval)          (1.f)
#endif

//<<<CHARGE_SOURCE_DEFINE
// Runtime charge-source switch: useTypeQ kernel arg selects between atom.w and type.z
// Old compile-time option kept for reference (disabled):
// #ifndef USE_CHARGE_FROM_TYPE
// #define USE_CHARGE_FROM_TYPE 0
// #endif
// #if USE_CHARGE_FROM_TYPE
// #define ASSIGN_Q_FROM_SOURCE(REQ, atom, REQtype) (REQ).z = (REQtype).z
// #else
// #define ASSIGN_Q_FROM_SOURCE(REQ, atom, REQtype) (REQ).z = (atom).w
// #endif
#define ASSIGN_Q_FROM_SOURCE(REQ, atom, REQtype, useTypeQ) (REQ).z = ((useTypeQ) ? (REQtype).z : (atom).w)

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
    __global float2*  ErefW,     // 7: [nsamp]   {E,W} reference energies and weights for each sample (molecule,system)
    int useTypeQ                 // 8: runtime switch: 0=per-atom atoms.w, 1=type-based REQH.z
){
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LREQKS[32];   // local buffer for REQ parameters
    __local int    LTYPES[32];   // local buffer for type indices (tj) to avoid global reads in inner loop
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

    // if( iG == iDBG ){
    //     printf("-----------------\nGPU: OCL evalSampleDerivatives() iG %i i0,ni %3i,%3i j0,nj %3i,%3i ErefW %g %g\n", iG, i0,ni, j0,nj, ErefW[iS].x, ErefW[iS].y);
    //     for(int i=0; i<ni; i++){ int ia=i0+i; int it=atypes[ia]; printf("GPU: atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    //     for(int i=0; i<nj; i++){ int ia=j0+i; int it=atypes[ia]; printf("GPU: atom j %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    // }
    
    const int    ti    = atypes[iG];
    const float4 atomi = atoms [iG];
          float4 REQi  = tREQHs[ti];
    ASSIGN_Q_FROM_SOURCE(REQi, atomi, REQi, useTypeQ);
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
            ASSIGN_Q_FROM_SOURCE(REQj, atomj, REQj, useTypeQ);
            const int2   jep  = ieps[jj];
            if( jep.x >= 0 ){ REQj.z -= tREQHs[jep.x].z; } // subtract charge of electron pair 1
            if( jep.y >= 0 ){ REQj.z -= tREQHs[jep.y].z; } // subtract charge of electron pair 2
            LATOMS[iL]  = atomj;
            LREQKS[iL]  = REQj;
            LTYPES[iL]  = tj;
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
                // If atom types are identical across fragments, scaling by 2 accounts for both sides sharing the same type
                {
                    int tj = LTYPES[jl];
                    float eScale = (ti == tj) ? 2.f : 1.f;
                    fREQi.y += eScale * dE_dE0 * REQj.y;
                }
                fREQi.z += dE_dQ  *        REQj.z;
                fREQi.w += dE_dH2 *        REQj.w * sH;
                //atomic_add_float4(&dEdREQs[n1+j], fREQj);

                Ei += ELJ + Eel;

                //if( iG == iDBG ){  printf("GPU: j %2i ELJ %16.8e Eel %16.8e r %16.8e    \n", jl, ELJ, Eel, r); }
                //if( iG == iDBG ){  printf("GPU: iG %2i iS %2i iL %2i j %2i ELJ %16.8e Eel %16.8e r %16.8e  R0 %16.8e E0 %16.8e Q %16.8e H2 %16.8e     dE_dR0 %16.8e dE_dE0 %16.8e dE_dQ %16.8e dE_dH2 %16.8e \n", iG, iS, iL, jl, ELJ, Eel, r,  R0, E0, Q, H2, dE_dR0, dE_dE0, dE_dQ, dE_dH2);}
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    //if( iG == iDBG ){ printf("GPU: iG %2i iS %2i iL %2i Ei %16.8e fREQi( %16.8e %16.8e %16.8e %16.8e )\n", iG, iS, iL, Ei, fREQi.x, fREQi.y, fREQi.z, fREQi.w);}
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
        //if( iG == iDBG ){      printf("GPU: iG %2i iS %2i iL %2i dE %16.8e fREQi( %16.8e %16.8e %16.8e %16.8e )\n", iG, iS, iL, LdE, fREQi.x, fREQi.y, fREQi.z, fREQi.w); }
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
    __global int4*    ranges,  // 1: [nsample] (i0, j0, ni, nj) but stored as (x=i0, y=j0, z=ni, w=nj)
    __global float4*  tREQHs,  // 2: [ntypes] REQH parameters
    __global int*     atypes,  // 3: [nAtomTot] atom types
    __global int2*    ieps,    // 4: [nAtomTot] {iep1,iep2} index of electron pair type ( used to subtract  charge)
    __global float4*  atoms,   // 5: [nAtomTot] positions and REPS charge of each atom (x,y,z,Q) 
    __global float4*  dEdREQs, // 6: [nAtomTot] derivatives of REQH parameters
    __global float2*  ErefW,   // 7: [nsamp] {E,W} reference energies and weights for each sample (molecule,system)
    __global float*   Emols,   // 8: [nsample] per-sample energies
    __global float*   Jmols,   // 9: [nSamples] per-sample objective contributions: 0.5*(Emol - Eref)*LdE
    int      useTypeQ,         //10: runtime switch: 0=per-atom atoms.w, 1=type-based REQH.z
    float4   globParams        //11: {alpha,?,?,?} global parameters (min,max, xlo,xhi, Klo,Khi, K0,x0)
){
    __local float4 LATOMS[32];
    __local float4 LREQKS[32];
    __local int    LTYPES[32];
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

    // if((iS==iDBG) && (iL==0)){  
    //     printf("-----------------\nGPU: evalSampleDerivatives_template() nG %7i nL %2i nS %6i | i0=%d ni=%d j0=%d nj=%d \n", (int)get_global_size(0), (int)get_local_size(0), (int)get_num_groups(0), i0, ni, j0, nj);
    //     //printf("GPU: evalSampleDerivatives_template() i0=%d ni=%d j0=%d nj=%d \n", i0, ni, j0, nj);
    //     for(int i=0; i<ni; i++){ int ia=i0+i; int it=atypes[ia]; printf("GPU: atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    //     for(int i=0; i<nj; i++){ int ia=j0+i; int it=atypes[ia]; printf("GPU: atom j %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    // }

    if( iG >= ni ) return;

    const int    ia    = i0 + iG; 
    const int    ti    = atypes[ia];
    const float4 atomi = atoms [ia];
          float4 REQi  = tREQHs[ti];
    
    // DEBUG: Print tREQHs values for first atom
    //if(iG == 0) { printf("GPU eval: tREQHs[0] = (%.8e, %.8e, %.8e, %.8e)\n",  REQi.x, REQi.y, REQi.z, REQi.w); }
    
    ASSIGN_Q_FROM_SOURCE(REQi, atomi, REQi, useTypeQ);
    const int2   iep   = ieps  [iG];
    if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
    if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }

    float4 fREQi = float4Zero;
    float  Ei    = 0.0f;
    int    cH    = 0;              // count how many j-pairs pass H-gate

    for(int off=0; off<nj; off+=nL){
        const int local_j = off + iL;
        if(local_j<nj){
            const int jj       = j0 + local_j;
            const int    tj    = atypes[jj];
            const float4 atomj = atoms [jj];
                  float4 REQj  = tREQHs[tj];
            ASSIGN_Q_FROM_SOURCE(REQj, atomj, REQj, useTypeQ);
            const int2   jep  = ieps[jj];
            if( jep.x >= 0 ){ REQj.z -= tREQHs[jep.x].z; }
            if( jep.y >= 0 ){ REQj.z -= tREQHs[jep.y].z; }
            LATOMS[iL]  = atomj;
            LREQKS[iL]  = REQj;
            LTYPES[iL]  = tj;
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
                // Mixing rules and hydrogen-bond gating (compile-time switch)
                float sH = S_H(H);
                H = APPLY_H_GATE(H);
                cH += (sH > 0.f);

                float4 fij; float Eij;
                //<<<MODEL_PAIR_ACCUMULATION

                Ei   +=Eij;
                fREQi+=fij;        
                
                // DEBUG: Print energy after pair
                //if(iS == 0 && iL == 0 && jl == 0) { printf("GPU: Energy after pair: Ei=%.8e, dEi=%.8e\\n", Ei, Ei - Ei0); }
                // --- Debug: print per-pair contributions (delta Ei and delta fREQi)
                if(iS==iDBG){
                    int tj = LTYPES[jl];
                    // printf("GPU: [%3i,%3i] tij(%2i,%2i) r %+10.6f REQH(%+10.6f,%+12.6e,%+12.6e,%+12.6e) | dE %+12.6e | dREQH( %+12.6e, %+12.6e, %+12.6e, %+12.6e)\n",  ia, j0+local_j2,  ti,tj, r, R0, E0, Q, H, Eij, fij.x, fij.y, fij.z, fij.w );
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    barrier(CLK_LOCAL_MEM_FENCE);
    printf("GPU: iS %3i atom %3i ti %3i Ei %12.6e fREQi( %+12.6e, %+12.6e, %+12.6e, %+12.6e )\n", iS, iL, ti, Ei, fREQi.x, fREQi.y, fREQi.z, fREQi.w ); 
    LATOMS[iL].x = Ei;
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL == 0) {
        float Emol = 0.0f;
        for(int i=0; i<ni; i++){ Emol += LATOMS[i].x; }
        float2 EW  = ErefW[iS];
        float dE   = (Emol - EW.x);
        LdE        = dE * EW.y;                // W*(Emol - Eref)
        float Jmol = 0.5f * dE * LdE;          // 0.5 * W * (Emol - Eref)^2
        Emols[iS]  = Emol;
        Jmols[iS]  = Jmol;
        if( (iS==iDBG) ){ printf("GPU: iS %2i Emol %16.8e Eref %16.8e dE %16.8e LdE %16.8e Jmol %16.8e\n", iS, Emol, EW.x, dE, LdE, Jmol);}
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    if (iL < ni) {
        //if( iG == iDBG ){  printf("GPU: fREQi %16.8e LdE %16.8e\n", fREQi, LdE); }
        //if( iG == iDBG ){ printf("GPU: iG %2i iS %2i iL %2i dE %16.8e fREQi( %16.8e %16.8e %16.8e %16.8e )  cH %d HBOND_GATE %d\n", iG, iS, iL, LdE, fREQi.x, fREQi.y, fREQi.z, fREQi.w, cH, HBOND_GATE); }
        //printf("GPU: iG %2i iS %2i iL %2i dE %16.8e fREQi( %16.8e %16.8e %16.8e %16.8e )\n", iG, iS, iL, LdE, fREQi.x, fREQi.y, fREQi.z, fREQi.w);
        dEdREQs[ia] = fREQi * LdE;
        // Debug: dump a few resulting dEdREQs entries for the debug sample
        //if( (iS==iDBG) && (iL < 4) ){float4 v = dEdREQs[i0 + iL]; printf("GPU: iS %2i lane %2i write dEdREQs[%3d] = ( %16.8e %16.8e %16.8e %16.8e )\n", iS, iL, i0+iL, v.x, v.y, v.z, v.w);}
    }
}


// -----------------------------------------------------------------------------
// Standalone SERIAL variant for deterministic debugging.
// Executes entirely on lane 0 of each workgroup to avoid barrier hazards.
// Performs two passes:
//   1) Compute Emol over all i–j pairs
//   2) Recompute pairs, accumulate per-atom derivatives, and write dEdREQs = fREQi * LdE
// Also writes per-sample Jmols[iS].
// Model-specific pair accumulation is injected at:
//     //<<<MODEL_PAIR_ACCUMULATION
// using the same in-scope variables as the templated kernel.
// -----------------------------------------------------------------------------
__attribute__((reqd_work_group_size(1,1,1)))
__kernel void evalSampleDerivatives_template_serial(
    __global const int4*   ranges,    // 1: [nsample] (i0, j0, ni, nj)
    __global const float4* tREQHs,    // 2: [ntype] (R0, E0, Q, H)
    __global const int*    atypes,    // 3: [nAtomTot] atom types
    __global const int2*   ieps,      // 4: [nAtomTot] {iep1,iep2} index of electron pair type ( used to subtract  charge)
    __global const float4* atoms,     // 5: [nAtomTot] positions and REPS charge of each atom (x,y,z,Q) 
    __global float4*       dEdREQs,   // 6: [nAtomTot] derivatives of REQH parameters
    __global const float2* ErefW,     // 7: [nsample] (Eref, W) per-sample reference energy and weight
    __global float*        Emols,     // 8: [nsample] per-sample energies
    __global float*        Jmols,     // 9: [nsample] per-sample objective contributions: 0.5*(Emol - Eref)*LdE
    const int              useTypeQ,  //10: runtime switch: 0=per-atom atoms.w, 1=type-based REQH.z
    const float4           globParams //11: {alpha,?,?,?} global parameters (min,max, xlo,xhi, Klo,Khi, K0,x0)
){
    const int iS = get_global_id(0);
    
    const int4 nsi = ranges[iS];
    const int i0 = nsi.x;
    const int ni = nsi.z;
    const int j0 = nsi.y;
    const int nj = nsi.w;

    printf("-----------------\nGPU: evalSampleDerivatives_template_serial() useTypeQ %i iS %i/3i | i0=%i ni=%i j0=%i nj=%i \n", useTypeQ, iS, get_global_size(0), i0, ni, j0, nj);
    for(int i=0; i<ni; i++){ int ia=i0+i; int it=atypes[ia]; printf("GPU: atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    for(int i=0; i<nj; i++){ int ia=j0+i; int it=atypes[ia]; printf("GPU: atom j %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, atoms[ia].x, atoms[ia].y, atoms[ia].z, atoms[ia].w, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }


    float Emol = 0.0f;
    
    int i0s = j0;
    int j0s = i0;
    int nis = nj;
    int njs = ni;
    float EmolS = 0.0f;
    for(int ii=0; ii<nis; ii++){
        int ia = i0s + ii;
        int ti = atypes[ia];
        float4 atomi = atoms [ia];
        float4 REQi  = tREQHs[ti];
        ASSIGN_Q_FROM_SOURCE(REQi, atomi, REQi, useTypeQ);
        int2  iep = ieps[ia];
        if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
        if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }

        float4 fREQi = float4Zero;
        float  Ei    = 0.0f;
        int    cH    = 0;

        for(int jj=0; jj<njs; jj++){
            int ja = j0s + jj;
            int tj = atypes[ja];
            float4 atomj = atoms [ja];
            float4 REQj  = tREQHs[tj];
            ASSIGN_Q_FROM_SOURCE(REQj, atomj, REQj, useTypeQ);
            int2  jep = ieps[ja];
            if( jep.x >= 0 ){ REQj.z -= tREQHs[jep.x].z; }
            if( jep.y >= 0 ){ REQj.z -= tREQHs[jep.y].z; }

            float3 dij = atomj.xyz - atomi.xyz;
            float  r   = length(dij);
            float  inv_r  = 1.f/fmax(r, R2SAFE);


            float R0 = REQi.x + REQj.x;
            float E0 = REQi.y * REQj.y;
            float Q  = REQi.z * REQj.z;
            float H  = REQi.w * REQj.w;
            float sH = S_H(H);
            H = APPLY_H_GATE(H);

            float  Eij = 0.f;
            float4 fij = float4Zero;

                //<<<MODEL_PAIR_ACCUMULATION
                
                Ei+=Eij;
                fREQi+=fij;        
                
                //printf("GPU: [%3i,%3i] tij(%2i,%2i) r %+10.6f REQH(%+10.6f,%+12.6e,%+12.6e,%+12.6e) | dE %+12.6e | dREQH( %+12.6e, %+12.6e, %+12.6e, %+12.6e)\n",  ia, ja,  ti,tj, r, R0, E0, Q, H, Eij, fij.x, fij.y, fij.z, fij.w );
                
        }
        Emol += Ei;
        dEdREQs[ia] = fREQi;
        printf("GPU: atom i %3i ti %3i Ei %12.6e fREQi( %+12.6e, %+12.6e, %+12.6e, %+12.6e )\n", ia, ti, Ei, fREQi.x, fREQi.y, fREQi.z, fREQi.w );        
    }
    
    float2 EW = ErefW[iS];
    float dE  = (Emol - EW.x);
    float LdE = dE * EW.y;
    float Jmol = 0.5f * dE * LdE;
    Emols[iS] = Emol;
    Jmols[iS] = Jmol;
    printf("GPU: evalSampleDerivatives_template_serial() iS %2i Emol %16.8e Eref %16.8e dE %16.8e LdE %16.8e Jmol %16.8e\n",  iS, Emol, EW.x, dE, LdE, Jmol);

    for(int ii=0; ii<nis; ii++){
        int ia = i0s + ii;
        dEdREQs[ia] *= LdE;
        int ti = atypes[ia];
        printf("GPU: atom i %3i ti %3i i0 %i dEdREQs( %+12.6e, %+12.6e, %+12.6e, %+12.6e ) \n", ia, ti, i0s, dEdREQs[ia].x, dEdREQs[ia].y, dEdREQs[ia].z, dEdREQs[ia].w);
    }
    
    //printf("GPU: evalSampleDerivatives_template_serial() iS %2i Emol %16.8e Eref %16.8e dE %16.8e LdE %16.8e Jmol %16.8e\n",  iS, Emol, EW.x, dE, LdE, Jmol);
    //for(int i=0; i<ni; i++){ int ia=i0+i; int it=atypes[ia]; printf("GPU: atom i %3i it %3i REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
    //for(int i=0; i<nj; i++){ int ia=j0+i; int it=atypes[ia]; printf("GPU: atom j %3i it %3i REQH %16.8f %16.8f %16.8f %16.8f \n", ia, it, tREQHs[it].x, tREQHs[it].y, tREQHs[it].z, tREQHs[it].w); }
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
    int      useTypeQ,           // runtime switch: 0=per-atom atoms.w, 1=type-based REQH.z
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
        
        // DEBUG: Print tREQHs values for first active work item
        if(iL == 0 && iS == 0) { printf("GPU energy_eval: tREQHs[%d] = (%.8e, %.8e, %.8e, %.8e)\n",   ti, REQi.x, REQi.y, REQi.z, REQi.w); }
        
        ASSIGN_Q_FROM_SOURCE(REQi, atomi, REQi, useTypeQ);
        const int2   iep   = ieps  [i];
        if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
        if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }
    }

    float  Ei    = 0.0f;

    float alpha = globParams.x;

    // --- Debug: print config for a chosen sample and few lanes ---
    //if((iS==iDBG) && (iL==0)){ 
    //     printf("-----------------\n GPU: evalSampleEnergy_template() nG %7i nL %2i nS %6i | i0=%d ni=%d j0=%d nj=%d | i=%d ti=%d\n", get_global_size(0), get_local_size(0), get_num_groups(0), i0, ni, j0, nj, i, ti); 
    //     for(int i=0; i<ni; i++){
    //         int ia=i0+i;
    //         int ti=atypes[ia];
    //         float4 atomi = atoms [ia];
    //         float4 REQi  = tREQHs[ti];
    //         REQi.z = atomi.w;
    //         const int2   iep   = ieps  [ia];
    //         if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
    //         if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }
    //         printf("GPU: frag1 atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, ti, atomi.x, atomi.y, atomi.z, atomi.w, REQi.x, REQi.y, REQi.z, REQi.w);
    //     }
    //     for(int i=0; i<nj; i++){
    //         int ia=j0+i;
    //         int ti=atypes[ia];
    //         float4 atomi = atoms [ia];
    //         float4 REQi  = tREQHs[ti];
    //         REQi.z = atomi.w;
    //         const int2   iep   = ieps  [ia];
    //         if( iep.x >= 0 ){ REQi.z -= tREQHs[iep.x].z; }
    //         if( iep.y >= 0 ){ REQi.z -= tREQHs[iep.y].z; }
    //         printf("GPU: frag2 atom i %3i it %3i pos %16.8f %16.8f %16.8f %16.8f  REQH %16.8f %16.8f %16.8f %16.8f \n", ia, ti, atomi.x, atomi.y, atomi.z, atomi.w, REQi.x, REQi.y, REQi.z, REQi.w);
    //     }
    //}

    for(int off=0; off<nj; off+=nL){
        const int local_j = off + iL;
        if(local_j<nj){
            const int jj       = j0 + local_j;
            const int    tj    = atypes[jj];
            const float4 atomj = atoms [jj];
                  float4 REQj  = tREQHs[tj];
            ASSIGN_Q_FROM_SOURCE(REQj, atomj, REQj, useTypeQ);
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
                // Mixing rules and hydrogen-bond gating (compile-time switch)
                float sH = S_H(H);
                H = APPLY_H_GATE(H);

                //if((iS==iDBG) && (iL==0)){  printf("GPU: evalSampleEnergy_template() ia,ja (%3i,%3i) R0 %16.8f E0 %16.8f Q %16.8f H %16.8f\n", i, jl+j0, R0, E0, Q, H);}

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
__kernel void assembleAndRegularize(
     int nDOFs,                       // 1: number of DOFs
     __global float*   fDOFs,         // 2: [nDOFs]    derivatives of REQH parameters
     __global int2*    DOFnis,        // 3: [nDOFs]    (i0,ni) star and number of eleemnts in block of DOFtoAtom corresponding to given DOF
     __global int*     DOFtoAtom,     // 4: [nInds]    list of atom indexes relevant for each DOF (non-uniform ranges indexed by DOFnis)
     __global float4*  DOFcofefs,     // 5: [nInds]    factors for update of each DOF from the atom dEdREQH parameters   fDOFi = dot( DOFcofefs[i], dEdREQH[i] )
     __global float4*  dEdREQs,       // 6: [nAtomTot] derivatives of REQH parameters
    __global const float*   DOFs,     // 7: [nDOFs]    current values of DOFs (need only for regularization)
    __global const float8*  regParams, // 8: [nDOFs]   {min,max, xlo,xhi, Klo,Khi, K0,x0}
    __global float4*        tREQHs,    // 9: [ntypes]  REQH parameters (to update with fitted values)
    __global int2*          DOFtoTypeComp // 10: [nDOFs] (atom_type, component) for each DOF
){
    __local float LfDOFi[NLOC_assembleDOFderivatives];


    const int iDOF = get_group_id(0);
    const int iL   = get_local_id(0);
    const int nL   = get_local_size(0);

    if( (iDOF==iDBG) && (iL==0) ){ 
        printf("-----------------\nGPU assembleAndRegularize().1 iDOF %2i / nDOFs %2i iL %2i nL %2i\n", iDOF, nDOFs, iL, nL); 
    //     for(int i=0; i<nDOFs; i++){
    //         printf("GPU assembleAndRegularize().2 DOFnis %2i %2i\n", DOFnis[i].x, DOFnis[i].y);
    //         int nsi = DOFnis[i].y; if(nsi>5){ nsi=5; }
    //         for(int jj=0; jj<nsi; jj++){
    //             int j = DOFnis[i].x + jj;
    //             printf("GPU assembleAndRegularize().3 DOF %2i j %2i DOFtoAtom %2i (%16.8f, %16.8f, %16.8f, %16.8f)\n", i,j, DOFtoAtom[j], DOFcofefs[j].x, DOFcofefs[j].y, DOFcofefs[j].z, DOFcofefs[j].w);
    //         }
    //     }
    }

    // --- Part 1: Assemble Physical Forces (same as before) ---
    const int2 nsi = DOFnis[iDOF];
    float fl = 0.0f;
    for(int i = iL; i < nsi.y; i += nL){
        int j   = i + nsi.x;
        int ia  = DOFtoAtom[j];
        fl     += dot(DOFcofefs[j], dEdREQs[ia]);
        //if( (iDOF==iDBG) && (iL==0) ){ printf("GPU assembleAndRegularize().4 iDOF %2i iL %2i ia %2i fl %16.8f dEdREQs(%16.8f, %16.8f, %16.8f, %16.8f) DOFcofefs(%16.8f, %16.8f, %16.8f, %16.8f)\n", iDOF, iL, ia, fl, dEdREQs[ia].x,dEdREQs[ia].y,dEdREQs[ia].z,dEdREQs[ia].w, DOFcofefs[j].x, DOFcofefs[j].y, DOFcofefs[j].z, DOFcofefs[j].w); }
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
        
        // --- Update tREQHs with fitted parameter values ---
        int2 type_comp = DOFtoTypeComp[iDOF];
        int atom_type = type_comp.x;
        int component = type_comp.y;
        float xDOF = DOFs[iDOF];    
        // Update the corresponding component in tREQHs
        if     (component == 0) tREQHs[atom_type].x = xDOF; // R
        else if(component == 1) tREQHs[atom_type].y = xDOF; // E (sqrt)
        else if(component == 2) tREQHs[atom_type].z = xDOF; // Q
        else if(component == 3) tREQHs[atom_type].w = xDOF; // H

        // DEBUG: Print before update
        if(iDOF == 0) {
            printf("GPU assembleAndRegularize: Updating tREQHs[%d].%c from %.8e to %.8e (DOF value: %.8e)\n", 
                   atom_type, "REQH"[component], 
                   (component==0 ? tREQHs[atom_type].x : 
                    component==1 ? tREQHs[atom_type].y :
                    component==2 ? tREQHs[atom_type].z : tREQHs[atom_type].w),
                   xDOF, DOFs[iDOF]);
        }
    }
}