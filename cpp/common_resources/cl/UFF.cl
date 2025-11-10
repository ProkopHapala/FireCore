#pragma OPENCL EXTENSION cl_khr_fp64 : enable // Enable if float precision is needed


// ======================================================
// Debug Controls (compile-time macros)
// ======================================================
// Enable concise debug prints without changing C++ host interface.
#define DBG_UFF    0         // 0/1 master switch
#define IDBG_ATOM  (-5)    // atom index to trace
#define IDBG_BOND  (0)    // bond index to trace (global bond id), -1 disables
#define IDBG_ANGLE (0)    // angle index to trace
#define IDBG_DIH   (0)    // dihedral index to trace
#define IDBG_INV   (0)    // inversion index to trace
#define IDBG_SYS   (0)    // system index to trace
//#define IDBG_SYS   (8)    // system index to trace

// ======================================================
// Utility Functions (Translate/Provide these accurately)
// ======================================================
// Assuming getLJQH, mixREQ_*, clampForce, float3Zero are defined as before
// and match the C++ behavior precisely, including precision (float/float).
// Using 'float' here for consistency with previous examples, adjust if using float.

typedef struct __attribute__ ((packed)){
    float4 a;
    float4 b;
    float4 c;
} cl_Mat3;

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f       // [ eV*Ang/e^2 ]
#define const_kB        8.617333262145e-5f   // [ eV/K ]


// Simple Arithmetic Mixing Rule (Placeholder for _mixREQ)
// REQ = {R0, E0, Q}
inline float4 mixREQ_arithmetic( float4 REQi, float4 REQj ){
    //float R0 = 0.5f * (REQi.x + REQj.x); // Arithmetic mean for R0 (sigma)
    //float E0 = sqrt(REQi.y * REQj.y); // Geometric mean for E0 (epsilon) // THIS IS WRONG, WE DO NOT DO THIS FOR PERFORMANCE REASONS
    float R0 = REQi.x + REQj.x;          // This is correct
    float E0 = REQi.y * REQj.y;          // This is correct we assume Ei,Ej are already pre-sqrt like Ei=sqrt(Eii)
    float Q  = REQi.z * REQj.z;          // Product of charges for Coulomb term in getLJQH
    float H  = REQi.w * REQj.w; if(H>0){H=0;}
    return (float4)(R0, E0, Q, H);
}

// Force clamping (Placeholder)
inline float3 clampForce( float3 f, float Fmax2 ){
    float f2 = dot(f,f);
    if(f2 > Fmax2){
        f *= sqrt(Fmax2/f2);
    }
    return f;
}

//inline float3 float3Zero(){ return (float3)(0.0f, 0.0f, 0.0f); }

// Complex number division for unit complex numbers (backward rotation)
inline float2 udiv_cmplx( float2 a, float2 b ){ return (float2)( a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y ); }

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

inline float4 getMorseQH( float3 dp,  float4 REQH, float K, float R2damp ){
    float r2    = dot(dp,dp);
    float ir2_  = 1/(r2+R2damp);
    float r     = sqrt( r2   );
    float ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    float e     = exp ( K*(r-REQH.x));
    //double e2    = e*e;
    //double fMors =  E0*  2*K*( e2 -   e ); // Morse
    //double EMors =  E0*      ( e2 - 2*e );
    float   Ae  = REQH.y*e;
    float fMors = Ae*  2*K*(e - 1); // Morse
    float EMors = Ae*      (e - 2);
    float Eel   = COULOMB_CONST*REQH.z*ir_;
    float fr    = fMors/r - Eel*ir2_ ;
    return  (float4){ dp*fr, EMors+Eel };
}

// inline float4 getLJQH( float3 dp, float4 REQij, float R2damp ){ /* ... As defined before ... */ }
// inline float4 mixREQ_arithmetic( float4 REQi, float4 REQj ){ /* ... As defined before ... */ } // Use the correct UFF mixing rule
// inline float3 clampForce( float3 f, float Fmax2 ){ /* ... As defined before ... */ }
// inline float3 float3Zero(){ return (float3)(0.0f, 0.0f, 0.0f); }
// inline float2 udiv_cmplx( float2 a, float2 b ){ return (float2)( a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y ); }

//#define GPU_PREFIX "GPU"

// ======================================================
// Kernels (Following C++ Structure Closely)
// ======================================================

// --- Clear kernels: zero force/energy buffers ---
__kernel void clear_fapos_UFF(const int n, __global float4* fapos){    
    int i = get_global_id(0);
    if(i>=n) return;
    fapos[i] = (float4)(0.f,0.f,0.f,0.f);
}

__kernel void clear_fint_UFF(const int n, __global float4* fint){
    int i = get_global_id(0);
    if(i>=n) return;
    fint[i] = (float4)(0.f,0.f,0.f,0.f);
}

// --- Kernel 1: Evaluate Bonds and Calculate H-Neigh vectors (Atom-Centric) ---
// Corresponds to the loop calling C++ evalAtomBonds
__kernel void evalBondsAndHNeigh_UFF(
    const int        natoms,
    const int        npbc,
    const int        i0bon,       // Offset for bond forces in fint array
    const int        bSubtractBondNonBond, // Use int 0/1 for bool
    const float      Rdamp,
    const float      FmaxNonBonded,
    // --- Input Arrays ---
    __global float4* apos,        // Atom positions {x,y,z,w}
    __global float4* fapos,       // [natoms] {fx,fy,fz,0} Accumulated bond forces
    __global int4*   neighs,      // [natoms] {ix,iy,iz,iw} Neighbor atom indices
    __global int4*   neighCell,   // [natoms] {icx,icy,icz,icw} Neighbor cell indices
    __global float4* pbc_shifts,  // [npbc] PBC shift vectors
    __global int4*   neighBs,     // [natoms] {ibx,iby,ibz,ibw} Precomputed bond indices for neighbors
    __global float2* bonParams,   // [nbonds] { K, l0 }
    __global float4* REQs,        // [natoms] {R0, E0, Q} (Passed if bSubtractBondNonBond=1)
    __global int2*   bonAtoms,    // [nbonds] {i, j}
    // --- Output Arrays ---
    __global float4* hneigh,      // Output: [natoms*4] bond vectors {nx,ny,nz, 1/L}
    __global float4* fint        // Output: [nf] Stores force pieces {fx,fy,fz, E_contrib}
    //__global float*  Eb_contrib   // Output: [natoms] Optional per-atom energy contribution buffer
) {
    int ia    = get_global_id(0);
    int isys = get_global_id(1);
    int i0a = isys * natoms;          // atoms per system (base index)
    int i0h = isys * (natoms*4);      // hneigh per system (base index)
    // Header print by a single work-item to avoid async interleaving (selected system only)
    // if ((DBG_UFF!=0) && ia==IDBG_ATOM && (isys==IDBG_SYS)){
    //     printf("GPU evalBondsAndHNeigh_UFF() natoms=%d npbc=%d i0bon=%d Rdamp=% .6e Fmax=% .6e bSubtractBondNonBond=%d iDBG=%d isys=%d\n", natoms, npbc, i0bon, Rdamp, FmaxNonBonded, bSubtractBondNonBond, IDBG_ATOM, isys);
    //     // Print first 64 bond parameter rows on one line per bond (safe upper bound without host arg)
    //     printf("GPU BOND-TABLE  ib   ia   ja           K           l0\n");

    //     // this is loop over bonds, but we should loop over bonds of atoms instead
    //     // for(int ib=0; ib<N; ++ib){
    //     //     int2 a   = bonAtoms[ib];
    //     //     float2 p = bonParams[ib];
    //     //     if(a.x==-1) break; // End of bonds
    //     //     printf("[%s][BOND-TABLE] %4d %4d %4d % .9e % .9e\n", GPU_PREFIX, ib, a.x, a.y, p.x, p.y);
    //     // }

    //     for (int ia=0; ia<natoms; ++ia){
    //         int4 ng = neighs[i0a + ia];
    //         int4 ngC= neighCell[i0a + ia];
    //         printf("GPU ATOM %3d : ng={%3d,%3d,%3d,%3d} ngC={%3d,%3d,%3d,%3d} pos{%8.4f,%8.4f,%8.4f} isys=%d ", ia, ng.x, ng.y, ng.z, ng.w, ngC.x, ngC.y, ngC.z, ngC.w, apos[i0a+ia].x, apos[i0a+ia].y, apos[i0a+ia].z, isys);
    //         for(int in=0; in<4; ++in){
    //             int ing = ng[in];
    //             if(ing<0) break;
    //             // bond params
    //             float2 bp = bonParams[ing];
    //             printf(" k,l[%i](%8.4f,%8.4f)", in, bp.x, bp.y);
    //         }
    //         printf("\n");
    //     }
    //     printf("GPU evalBondsAndHNeigh_UFF().eval \n");
    // }
    if (ia >= natoms) return;

    float E_b = 0.0f; // Accumulate bond energy for this atom
    float3 pa = apos[i0a + ia].xyz;
    int4   ng = neighs[i0a + ia];
    int4   ngC= neighCell[i0a + ia];
    int4   nB = neighBs[i0a + ia];   // Precomputed bond indices for neighbors (per-system)
    int*   ings = (int*)&ng;
    int*   ingC = (int*)&ngC;
    int*   inbs = (int*)&nB;   // Bond indices ibx, iby, ibz, ibw

    const float4 REQi = (0) ? REQs[i0a + ia] : (float4)(0.0f); // Load only if needed
    const float R2damp = Rdamp * Rdamp;
    const float Fmax2  = FmaxNonBonded * FmaxNonBonded;

    int i4_h = (i0a + ia) * 4; // Base index for hneigh output for this atom
    float3 ftot = (float3){0.0f,0.0f,0.0f};

    for (int in = 0; in < 4; ++in) { // Loop over neighbor slots 0..3
        int ing = ings[in]; // Neighbor atom index
        if (ing < 0) {
            hneigh[i4_h + in] = (float4)(0.0f, 0.0f, 0.0f, 1e+6f); // Mark as invalid/unused hneigh
            continue;
        }

        // --- Bond vectors (Calculate and Store HNeigh) ---
        int hneigh_idx = i4_h + in;
        float3 png = apos[i0a + ing].xyz;
        float3 dp = png - pa;

        // Apply PBC shift
        int ipbc_idx = ingC[in];
        if (ipbc_idx >= 0 && ipbc_idx < npbc) {
             dp += pbc_shifts[ipbc_idx].xyz;
        }

        float l2 = dot(dp, dp);
        if (l2 < 1e-16f) l2 = 1e-16f;
        float l = sqrt(l2);
        float inv_l = 1.0f / l;
        float3 h = dp * inv_l; // Normalized bond vector `h`
        hneigh[hneigh_idx] = (float4)(h, inv_l); // Store {nx,ny,nz, 1/L}

        // --- Bond Energy & Force ---
        int ib = inbs[in];    // Precomputed bond index
        if (ib < 0) continue; // Should not happen for valid neighbor

        float2 param = bonParams[ib]; // { K, l0 } (ib is already per-system via neighBs packing)
        float k  = param.x;
        float l0 = param.y;
        float dl = l - l0;
        float fr = 2.0f * k * dl; // Force magnitude dE/dl
        float3 f_bond = h * fr; // Magnitude-direction along (ia->ing)

        float E = k * dl * dl; // Harmonic energy E = k*dl^2

        // CPU (UFF.h::evalAtomBonds): f set along dp/l and added to fapos[ia]
        // Therefore, force ON atom ia is +f_bond, and on neighbor ing is -f_bond
        float3 fi =  f_bond; // Force component ON atom `ia` from this bond
        float3 fj = -f_bond; // Force component ON atom `ing` from this bond

        float Enb = 0.0f;
        // --- Subtract Non-Bonded Interaction ---
        if (0*bSubtractBondNonBond != 0) {
            float4 REQj = REQs[i0a + ing];
            float4 REQij = mixREQ_arithmetic(REQi, REQj); // Use correct mixing rule
            //printf("GPU: REQi=(%g,%g,%g) REQj=(%g,%g,%g) REQij=(%g,%g,%g)\n", REQi.x,REQi.y,REQi.z,  REQj.x,REQj.y,REQj.z,  REQij.x,REQij.y,REQij.z );
            // dp vector already calculated
            float4 fnb = getLJQH(dp, REQij, R2damp); // Returns {force_on_j, Energy}
            Enb = fnb.w;
            float3 fnb_clamped = clampForce(fnb.xyz, Fmax2);

            // CPU does: f (on ia) -= fnb
            // Using our convention fi=+f_bond (on ia), fj=-f_bond (on ing):
            // fi =  f_bond - fnb_clamped
            // fj = -f_bond + fnb_clamped
            fi =  f_bond - fnb_clamped;
            fj = -f_bond + fnb_clamped;
        }

        // Per-DOF concise, aligned debug print: only for selected bond index
        // if ((DBG_UFF!=0) && (isys==IDBG_SYS)){
        //     int ib_dbg = inbs[in];
        //     if (ib_dbg == IDBG_BOND){
        //         int2 aij = bonAtoms[ib_dbg];
        //         //printf("GPU BOND %3i : ia=%3i ja=%3i  k=% .4e l0=% .4e  l=% .4e dl=% .4e  fr=% .4e  Enb=% .4e  fi=(% .4e % .4e % .4e)  fj=(% .4e % .4e % .4e)  E=% .4e isys=%d\n",
        //         //       ib_dbg, aij.x, aij.y, k, l0, l, dl, fr, Enb, fi.x,fi.y,fi.z, fj.x,fj.y,fj.z, (E-Enb), isys);
        //     }
        // }

        ftot += fi;
        float E_contrib = (E - Enb) * 0.5f; // Energy contribution per atom
        E_b += (E - Enb); // Accumulate total bond energy for atom ia

        // this is not done in UFF.cl assembleAtomForce, so I don't see why should we do it here. Keep it simple !!!
        // // Store bond recoil forces into fint
        // {
        //     int idx0 = i0bon + ib * 2;
        //     int idx1 = idx0 + 1;
        //     int2 atoms_ib = bonAtoms[ib];
        //     if (atoms_ib.x == ia) {
        //         fint[idx0] = (float4)(fi.x, fi.y, fi.z, E_contrib);
        //         fint[idx1] = (float4)(fj.x, fj.y, fj.z, E_contrib);
        //     } else {
        //         fint[idx1] = (float4)(fi.x, fi.y, fi.z, E_contrib);
        //         fint[idx0] = (float4)(fj.x, fj.y, fj.z, E_contrib);
        //     }
        // }
    } // End loop over neighbors


    // After loop:
    //fapos[i0a + ia] += (float4)(ftot, E_b); // Accumulate bond force directly into this system slice

    fapos[i0a + ia] = (float4)(ftot, E_b); // If we call this kernel first, we can set it, and ommit call to clean-forces, but be carefull when calling other contribution (Non-Bonded, Angles, Dihedrals, Assemble etc.) to call it AFTER this kernel.

    //if (Eb_contrib) Eb_contrib[ia] = E_b;
    //if (Eb_contrib) Eb_contrib[ia] = E_b * 0.5f; // Energy contribution is per-bond, associate half with ia? Check C++ E calculation. C++ adds E per bond.
                                                 // Let's store total bond energy involving ia.


}


// --- Kernel 2: Evaluate Angles (Interaction-Centric, Inlined) ---
__kernel void evalAngles_UFF(
    const int        nangles,
    const int        i0ang,       // Offset for angle forces in fint array
    const int        bSubtractAngleNonBond, // 0/1
    const float      Rdamp,
    const float      FmaxNonBonded,
    // --- Input Arrays ---
    __global int4*   angAtoms,    // [nangles] {ia, ja, ka, 0} where ja is central
    __global int2*   angNgs,      // [nangles] Precomputed {hneigh_idx_ji, hneigh_idx_kj}
    __global float4* angParams1,  // [nangles] {c0, c1, c2, c3}
    __global float*  angParams2_w,// [nangles] {K}
    __global float4* hneigh,      // Input: [natoms*4] Precomputed h-vectors
    // --- Optional NB Subtraction Inputs (pass only if bSubtractAngleNonBond=1) ---
    __global float4* REQs,
    __global float4* apos,
    __global float4* pbc_shifts,
    __global int4*   neighs,
    __global int4*   neighCell,
    const int        npbc,
    // --- Output Arrays ---
    __global float4* fint,        // Output: [nf] Stores force pieces {fx,fy,fz, E_contrib}
    __global float*  Ea_contrib,  // Output: [nangles] Optional per-angle energy buffer
    const int        nf_per_system
){
    int iang = get_global_id(0);
    int isys = get_global_id(1);
    int i0A = isys * nangles;       // per-system base for angle arrays
    int i0f = isys * nf_per_system; // per-system base for fint
    // if ((DBG_UFF!=0) && iang==IDBG_ANGLE && (isys==IDBG_SYS)){
    //     printf("GPU evalAngles_UFF() nangles=%3i i0ang=%3i Rdamp=% .4e Fmax=% .4e bSubtractAngleNonBond=%d iDBG=%d isys=%d\n", nangles, i0ang, Rdamp, FmaxNonBonded, bSubtractAngleNonBond, IDBG_ANGLE, isys);
    //     printf("GPU ANG-TABLE  id   ia   ja   ka            K          c0          c1          c2          c3\n");
    //     int N = (nangles<64)?nangles:64;
    //     for(int i=0;i<N;i++){
    //         int4 a = angAtoms[i0A + i]; int ia0=a.x, ja0=a.y, ka0=a.z;
    //         float4 cs=angParams1[i0A + i]; float K=angParams2_w[i0A + i];
    //         printf("GPU ANG %3i : ia=%3i ja=%3i ka=%3i  K=% .4e c0=% .4e c1=% .4e c2=% .4e c3=% .4e\n", i, ia0,ja0,ka0,K,cs.x,cs.y,cs.z,cs.w);
    //     }
    //     printf("evalAngles_UFF().eval\n");
    // }
    if (iang >= nangles) return;

    // --- Get Data ---
    int4 a = angAtoms[i0A + iang];
    int ia = a.x; // Atom i
    int ja = a.y; // Atom j (central)
    int ka = a.z; // Atom k

    int2 ngs = angNgs[i0A + iang];         // Precomputed hneigh indices {ji, kj}
    float4 qij = hneigh[ngs.x]; // indices already include per-system base from host
    float4 qkj = hneigh[ngs.y];

    // Non-bonded energy placeholder for debug printing
    float Enb = 0.0f;

    float4 par1 = angParams1[i0A + iang];   // {c0,c1,c2,c3}
    float  K     = angParams2_w[i0A + iang]; // K

    // --- Angle calculation exactly as CPU (UFF.h::evalAngle_Prokop) ---
    float3 fpi, fpj, fpk; // Forces on i, j, k
    float E = 0.0f;
    {
        // CPU computes c via h = qij + qkj
        float3 h = qij.xyz + qkj.xyz;
        float c = 0.5f * ( dot(h,h) - 2.0f );   // cos(theta)
        c = clamp(c, -1.0f, 1.0f);
        float s = sqrt(fmax(0.0f, 1.0f - c*c) + 1e-14f);
        float inv_s = (s > 1e-12f) ? (1.0f / s) : 0.0f;

        // E and f via Fourier series using complex multiplication, matching CPU
        float c0 = par1.x;
        float c1 = par1.y;
        float c2 = par1.z;
        float c3 = par1.w;

        // Represent cos/sin as a complex number: (x,y) = (cos, sin)
        float2 cs  = (float2)(c, s);
        float2 csn = cs; // will hold cos(n*theta), sin(n*theta)

        // Start with coefficients
        float Eloc = c0;
        float fmag = c1;

        // n=2 term
        // csn *= cs  => (cos2, sin2)
        csn = (float2)( csn.x*cs.x - csn.y*cs.y,  csn.x*cs.y + csn.y*cs.x );
        Eloc += c2 * csn.x;
        fmag += c2 * csn.y * inv_s * 2.0f;

        // n=3 term
        csn = (float2)( csn.x*cs.x - csn.y*cs.y,  csn.x*cs.y + csn.y*cs.x );
        Eloc += c3 * csn.x;
        fmag += c3 * csn.y * inv_s * 3.0f;

        // Scale by K
        E = K * Eloc;
        fmag *= K;

        // Assemble forces exactly as CPU
        float fi = fmag * qij.w;
        float fk = fmag * qkj.w;
        float fic = fi * c;
        float fkc = fk * c;

        // fpi = fic*qij - fi*qkj
        fpi = fic * qij.xyz - fi * qkj.xyz;
        // fpk = -fk*qij + fkc*qkj
        fpk = -fk * qij.xyz + fkc * qkj.xyz;
        // fpj = (fk-fic)*qij + (fi-fkc)*qkj
        fpj = (fk - fic) * qij.xyz + (fi - fkc) * qkj.xyz;
    }

    // if ((DBG_UFF!=0) && iang==IDBG_ANGLE && (isys==IDBG_SYS)){
    //     int4 a_dbg = angAtoms[iang];
    //     int ia0=a_dbg.x, ja0=a_dbg.y, ka0=a_dbg.z;
    //     float theta = acos(clamp(dot(qij.xyz,qkj.xyz),-1.0f,1.0f));
    //     printf("GPU ANG %3d : ia=%3d ja=%3d ka=%3d  K=% .4e c0=% .4e c1=% .4e c2=% .4e c3=% .4e  ang=% .4e  Enb=% .4e  fi=(% .4e % .4e % .4e)  fj=(% .4e % .4e % .4e)  fk=(% .4e % .4e % .4e)  E=% .4e isys=%d\n",
    //            iang, ia0,ja0,ka0,
    //            K, par1.x,par1.y,par1.z,par1.w,
    //            theta,Enb,
    //            fpi.x,fpi.y,fpi.z,
    //            fpj.x,fpj.y,fpj.z,
    //            fpk.x,fpk.y,fpk.z,
    //            (E-Enb), isys);
    // }
    // --- Subtract 1-3 Non-Bonded Interaction ---
    if (0*bSubtractAngleNonBond != 0) {
        // Needs REQs, apos, pbc_shifts, neighs, neighCell passed
        float4 REQi = REQs[ia];
        float4 REQk = REQs[ka];
        float4 REQik = mixREQ_arithmetic(REQi, REQk); // Use correct mixing rule

        // Use C++ Prokop reconstruction: dp = ik = ji - jk
        // Vec3d dp; dp.set_lincomb( (1./qij.w), qij.f, (-1./qkj.w), qkj.f );
        // Note: C++ qij.f is vector ji. So this is (1/lij)*ji - (1/lkj)*kj = r_ji - r_jk = r_ji + r_jk = r_ik ? Needs verification.
        // Let's assume this formula calculates the ik vector correctly *without* PBC initially.
        float3 dp_ik = (qij.xyz / qij.w) - (qkj.xyz / qkj.w); // Vector ik = ji - ki ? No, kj. ji - kj = ik. Seems correct.

        // Add PBC shift: shift(ik) = shift(jk) - shift(ji)
        // Need neighbor slots for i and k relative to central atom j (ja)
        // This still requires lookup if we don't pass slots/use apos.
        // Let's assume slots are looked up here (suboptimal but matches C++ possibility).
        int4 ngC_ja = neighCell[ja];
        int4 nbor_ja= neighs[ja];
        int in_slot = -1, kn_slot = -1;
        if(nbor_ja.x == ia) in_slot=0; else if(nbor_ja.y == ia) in_slot=1; else if(nbor_ja.z == ia) in_slot=2; else if(nbor_ja.w == ia) in_slot=3;
        if(nbor_ja.x == ka) kn_slot=0; else if(nbor_ja.y == ka) kn_slot=1; else if(nbor_ja.z == ka) kn_slot=2; else if(nbor_ja.w == ka) kn_slot=3;

        float3 shift_ji = (float3){0.0f,0.0f,0.0f};
        float3 shift_jk = (float3){0.0f,0.0f,0.0f};
        if(in_slot != -1) { int ipbc_ji = ((int*)&ngC_ja)[in_slot]; if(ipbc_ji>=0 && ipbc_ji<npbc) shift_ji = pbc_shifts[ipbc_ji].xyz; }
        if(kn_slot != -1) { int ipbc_jk = ((int*)&ngC_ja)[kn_slot]; if(ipbc_jk>=0 && ipbc_jk<npbc) shift_jk = pbc_shifts[ipbc_jk].xyz; }
        dp_ik += (shift_jk - shift_ji); // Apply combined shift

        float R2damp = Rdamp * Rdamp;
        float4 fnb = getLJQH(dp_ik, REQik, R2damp);
        Enb = fnb.w;

        float Fmax2 = FmaxNonBonded * FmaxNonBonded;
        float3 fnb_clamped = clampForce(fnb.xyz, Fmax2);

        // fnb is force ON k FROM i. Subtract this from angle forces.
        fpi += fnb_clamped; // Add force ON i FROM k (i.e., -fnb)
        fpk -= fnb_clamped; // Add force ON k FROM i (i.e., +fnb)
    }

    // --- Store forces into `fint` array ---
    float E_contrib = (E - Enb) / 3.0f;
    int idx_base = i0f + i0ang + iang * 3; // per-system fint offset
    fint[idx_base + 0] = (float4)(fpi, E_contrib); // Store force on atom i
    fint[idx_base + 1] = (float4)(fpj, E_contrib); // Store force on atom j
    fint[idx_base + 2] = (float4)(fpk, E_contrib); // Store force on atom k

    if (Ea_contrib) Ea_contrib[iang] = E - Enb; // Store per-angle energy
}


// --- Kernel 3: Evaluate Dihedrals (Interaction-Centric, Inlined) ---
__kernel void evalDihedrals_UFF(
    const int        ndihedrals,
    const int        i0dih,       // Offset for dihedral forces in fint array
    const float      SubNBTorsionFactor, // Use 0.0f to disable
    const float      Rdamp,
    const float      FmaxNonBonded,
    // --- Input Arrays ---
    __global int*    dihAtoms,    // [ndihedrals*4] {ia, ja, ka, la}
    __global int4*   dihNgs,      // [ndihedrals] Precomputed {hneigh_idx_ji, hneigh_idx_kj, hneigh_idx_lk, 0}
    __global float4* dihParams,   // [ndihedrals] { V, d=cos(n*phi0), n, w(ignored) }
    __global float4* hneigh,      // Input: [natoms*4] Precomputed h-vectors
    // --- Optional NB Subtraction Inputs (pass only if SubNBTorsionFactor > 0) ---
    __global float4* REQs,
    __global float4* apos,
    __global float4* pbc_shifts,
    __global int4*   neighs,
    __global int4*   neighCell,
    const int        npbc,
    // --- Output Arrays ---
    __global float4* fint,        // Output: [nf] Stores force pieces {fx,fy,fz, E_contrib}
    __global float*  Ed_contrib,  // Output: [ndihedrals] Optional per-dihedral energy buffer
    const int        nf_per_system
) {
    int id = get_global_id(0);
    int isys = get_global_id(1);
    int i0D = isys * ndihedrals;    // per-system base for dihedral arrays
    int i0f = isys * nf_per_system; // per-system base for fint
    // if ((DBG_UFF!=0) && id==IDBG_DIH && (isys==IDBG_SYS)){
    //     printf("GPU evalDihedrals_UFF() ndihedrals=%d i0dih=%d Rdamp=% .6e Fmax=% .6e SubNBTorsionFactor=% .6e iDBG=%d isys=%d\n", ndihedrals, i0dih, Rdamp, FmaxNonBonded, SubNBTorsionFactor, IDBG_DIH, isys);
    //     printf("GPU DIH-TABLE  ia   ja   ka   la            V           d           n\n");
    //     int N=(ndihedrals<64)?ndihedrals:64;
    //     for(int i=0;i<N;i++){ int ia=dihAtoms[(i0D+i)*4+0],ja=dihAtoms[(i0D+i)*4+1],ka=dihAtoms[(i0D+i)*4+2],la=dihAtoms[(i0D+i)*4+3]; float4 p=dihParams[i0D+i];
    //         printf("GPU DIH %3d : %3d %3d %3d %3d % .4e % .4e % .3f\n", i, ia,ja,ka,la, p.x,p.y,p.z);
    //     }
    //     printf("GPU evalDihedrals_UFF().eval\n");
    // }
    if (id >= ndihedrals) return;

    // --- Get Data --- (per-system)
    int i4a = (i0D + id) * 4;
    int ia = dihAtoms[i4a + 0]; // Atom i
    int ja = dihAtoms[i4a + 1]; // Atom j
    int ka = dihAtoms[i4a + 2]; // Atom k
    int la = dihAtoms[i4a + 3]; // Atom l
    int4 ngs = dihNgs[i0D + id];
    float3 par = dihParams[i0D + id].xyz;

    // Non-bonded energy placeholder for debug printing (used before subtraction block)
    float Enb = 0.0f;

    // Get precomputed hneigh indices {ji, kj, lk}
    float4 q12 = hneigh[ngs.x]; // ji {h_ji, 1/l_ji}
    float4 q32 = hneigh[ngs.y]; // kj {h_kj, 1/l_kj}
    float4 q43 = hneigh[ngs.z]; // lk {h_lk, 1/l_lk}

    //float3 par = dihParams[id]; // { V, d, n }

    // --- Inlined Dihedral Calculation (exact UFF.h::evalDihedral_Prokop) ---
    float3 fi, fj, fk, fl; // Forces on i, j, k, l
    float E = 0.0f;
    {
        // Compute in float precision to match CPU reference
        float3 h12 = convert_float3(q12.xyz);
        float3 h32 = convert_float3(q32.xyz);
        float3 h43 = convert_float3(q43.xyz);
        float3 n123 = cross(h12, h32);
        float3 n234 = cross(h43, h32);
        float n123_2 = dot(n123,n123);
        float n234_2 = dot(n234,n234);
        if (n123_2 < 1e-30f || n234_2 < 1e-30f){ fi=fj=fk=fl=(float3)(0.0f); E=0.0f; }
        else{
            float il2_123 = 1.0f / n123_2;
            float il2_234 = 1.0f / n234_2;
            float inv_n12 = sqrt(il2_123 * il2_234);
            float2 cs = (float2)( dot(n123,n234)*inv_n12,  -dot(n123,h43)*inv_n12 );
            float2 csn = cs;
            int   nint = (int)(par.z);
            for(int i=1;i<nint;i++){ csn = (float2)( csn.x*cs.x - csn.y*cs.y,  csn.x*cs.y + csn.y*cs.x ); }
            float Ed = (float)par.x * ( 1.0f + (float)par.y * csn.x );
            E = (float)Ed;
            float f = -(float)par.x * (float)par.y * (float)par.z * csn.y;
            // Forces on end atoms
            float3 fp1_loc_d = n123 * (-f * il2_123 * (float)q12.w);
            float3 fp4_loc_d = n234 * ( f * il2_234 * (float)q43.w);
            // Recoil on axis atoms
            float c123 = dot(h32,h12) * ((float)q32.w / (float)q12.w);
            float c432 = dot(h32,h43) * ((float)q32.w / (float)q43.w);
            float3 fp3_loc_d = fp1_loc_d * (-c123) + fp4_loc_d * (-c432 - 1.0f);
            float3 fp2_loc_d = fp1_loc_d * ( c123 - 1.0f) + fp4_loc_d * ( c432      );

            float3 fp1_loc = convert_float3(fp1_loc_d);
            float3 fp2_loc = convert_float3(fp2_loc_d);
            float3 fp3_loc = convert_float3(fp3_loc_d);
            float3 fp4_loc = convert_float3(fp4_loc_d);

            fi = fp1_loc; fj = fp2_loc; fk = fp3_loc; fl = fp4_loc;
        }
    }

    // if ((DBG_UFF!=0) && id==IDBG_DIH && (isys==IDBG_SYS)){
    //     int ia0=dihAtoms[(i0D+id)*4+0], ja0=dihAtoms[(i0D+id)*4+1], ka0=dihAtoms[(i0D+id)*4+2], la0=dihAtoms[(i0D+id)*4+3];
    //     // Recompute phi exactly as CPU Prokop path
    //     float3 h12d = convert_float3(q12.xyz); float3 h32d = convert_float3(q32.xyz); float3 h43d = convert_float3(q43.xyz);
    //     float3 n123d = cross(h12d, h32d);
    //     float3 n234d = cross(h43d, h32d);
    //     float n1_2 = dot(n123d,n123d); float n2_2 = dot(n234d,n234d);
    //     float cphi_d = 1.0f;
    //     if(n1_2>1e-30f && n2_2>1e-30f){ cphi_d = dot(n123d,n234d)/sqrt(n1_2*n2_2); cphi_d = clamp(cphi_d,-1.0f,1.0f); }
    //     float cphi = (float)cphi_d;
    //     float phi = acos(cphi);
    //     float3 par = dihParams[i0D + id].xyz;
    //     printf("GPU DIH %4d : ia=%4d ja=%4d ka=%4d la=%4d  V=% .4e d=% .4e n=% .3f  phi=% .4e  Enb=% .4e  fi=(% .4e % .4e % .4e)  fj=(% .4e % .4e % .4e)  fk=(% .4e % .4e % .4e)  fl=(% .4e % .4e % .4e)  E=% .4e isys=%d\n",
    //            id, ia0,ja0,ka0,la0,
    //            par.x,par.y,par.z,
    //            phi,Enb,
    //            fi.x,fi.y,fi.z,
    //            fj.x,fj.y,fj.z,
    //            fk.x,fk.y,fk.z,
    //            fl.x,fl.y,fl.z,
    //            (float)(E-Enb), isys);
    // }
    // --- Subtract 1-4 Non-Bonded Interaction ---
    if (SubNBTorsionFactor > 1e-6f) {
        // Needs REQs, apos, pbc_shifts, neighs, neighCell passed
        float4 REQi = REQs[ia];
        float4 REQl = REQs[la];
        float4 REQil = mixREQ_arithmetic(REQi, REQl); // Use correct mixing rule

        // Use C++ Prokop reconstruction: dp = il = -ji + jk - lk
        // Vec3d dp; dp.set_lincomb( (1./q12.w), (-1./q43.w), (-1./q32.w), q12.f, q43.f, q32.f );
        // Let's re-evaluate: (1/lji)*ji + (-1/llk)*lk + (-1/lkj)*kj
        // = r_ji - r_lk - r_kj = r_ji + r_kl + r_jk = r_ik + r_kl = r_il. Seems correct.
        float3 dp_il = (q12.xyz / q12.w) - (q43.xyz / q43.w) - (q32.xyz / q32.w);

        // Add PBC shift: shift(il) = shift(ij) + shift(jk) + shift(kl) = -shift(ji) + shift(jk) - shift(lk)
        // Need neighbor slots to find cell indices for ji, kj, lk bonds. Requires lookups.
        // Look up slots for (i relative to j), (k relative to j), (l relative to k)
        int4 ngC_ja = neighCell[ja]; int4 nbor_ja = neighs[ja];
        int4 ngC_ka = neighCell[ka]; int4 nbor_ka = neighs[ka];
        int in_slot_j = -1, kn_slot_j = -1, ln_slot_k = -1;
        // Find slots... (lookup logic omitted for brevity - assume we get them)
        // float3 shift_ji = ...; float3 shift_jk = ...; float3 shift_lk = ...;
        // dp_il += (-shift_ji + shift_jk - shift_lk); // Apply combined shift

        float R2damp = Rdamp * Rdamp;
        float4 fnb = getLJQH(dp_il, REQil, R2damp);
        Enb = fnb.w;

        float Fmax2 = FmaxNonBonded * FmaxNonBonded;
        float3 fnb_clamped = clampForce(fnb.xyz, Fmax2) * SubNBTorsionFactor;

        // fnb is force ON l FROM i. Subtract this from torsion forces.
        fi += fnb_clamped; // Add force ON i FROM l
        fl -= fnb_clamped; // Add force ON l FROM i
    }

    // --- Store forces into `fint` array ---
    float E_contrib = (E - Enb) * 0.25f;
    int idx_base = i0f + i0dih + id * 4; // per-system fint offset
    fint[idx_base + 0] = (float4)(fi, E_contrib); // Store force on atom i
    fint[idx_base + 1] = (float4)(fj, E_contrib); // Store force on atom j
    fint[idx_base + 2] = (float4)(fk, E_contrib); // Store force on atom k
    fint[idx_base + 3] = (float4)(fl, E_contrib); // Store force on atom l

    if (Ed_contrib) Ed_contrib[id] = E - Enb; // Store per-dihedral energy
}


// --- Kernel 4: Evaluate Inversions (Interaction-Centric, Inlined) ---
__kernel void evalInversions_UFF(
    const int        ninversions,
    const int        i0inv,       // Offset for inversion forces in fint array
    // --- Input Arrays ---
    __global int*    invAtoms,    // [ninversions*4] {ia, ja, ka, la} where ia is central
    __global int*    invNgs,      // [ninversions*3] Precomputed {hneigh_idx_ji, hneigh_idx_ki, hneigh_idx_li}
    __global float4* invParams,   // [ninversions] { K, c0, c1, c2 } -> Actually Quat4d in C++? {K, c0, c1, c2} assume float4
    __global float4* hneigh,      // Input: [natoms*4] Precomputed h-vectors
    // --- Output Arrays ---
    __global float4* fint,        // Output: [nf] Stores force pieces {fx,fy,fz, E_contrib}
    __global float*  Ei_contrib,  // Output: [ninversions] Optional per-inversion energy buffer
    const int        nf_per_system
) {
    int ii = get_global_id(0);
    int isys = get_global_id(1);
    int i0I = isys * ninversions;    // per-system base for inversion arrays
    int i0f = isys * nf_per_system;  // per-system base for fint
    if ((DBG_UFF!=0) && ii==IDBG_INV && (isys==IDBG_SYS)){
        printf("GPU INV  ninversions=%d i0inv=%d iDBG=%d isys=%d\n", ninversions, i0inv, IDBG_INV, isys);
        printf("GPU INV-TABLE  ia   ja   ka   la            K          c0          c1          c2\n");
        int N=(ninversions<64)?ninversions:64;
        for(int i=0;i<N;i++){ int ia=invAtoms[(i0I+i)*4+0],ja=invAtoms[(i0I+i)*4+1],ka=invAtoms[(i0I+i)*4+2],la=invAtoms[(i0I+i)*4+3]; float4 p=invParams[i0I+i];
            printf("GPU INV %3i : ia=%3i ja=%3i ka=%3i la=%3i  K=% .4e c0=% .4e c1=% .4e c2=% .4e\n", i, ia,ja,ka,la,p.x,p.y,p.z,p.w);
        }
        printf("GPU evalInversions_UFF().eval: \n");
    }
    if (ii >= ninversions) return;

    // --- Get Data ---
    int i4a = (i0I + ii) * 4;
    // Atoms ia, ja, ka, la (ia is central)
    int ia = invAtoms[i4a + 0]; // Atom i (central)
    int ja = invAtoms[i4a + 1]; // Atom j
    int ka = invAtoms[i4a + 2]; // Atom k
    int la = invAtoms[i4a + 3]; // Atom l

    // Get precomputed hneigh indices {ji, ki, li} relative to central atom ia
    int3 ngs = ((__global int3*)invNgs)[i0I + ii]; // Read as int3 (already offset by i0I)
    float4 q21 = hneigh[ngs.x]; // ji {h_ji, 1/l_ji}
    float4 q31 = hneigh[ngs.y]; // ki {h_ki, 1/l_ki}
    float4 q41 = hneigh[ngs.z]; // li {h_li, 1/l_li}

    float4 par = invParams[i0I + ii]; // { K, c0, c1, c2 }

    // --- Inlined Inversion Calculation (Prokop Style) ---
    float3 fp1, fp2, fp3, fp4; // Forces on i(1), j(2), k(3), l(4)
    float E = 0.0f;
    { // Start of inlined evalInversionUFF block
        // --- normal to plane jik (atoms 2,3,1)
        float3 n123 = cross(-q21.xyz, -q31.xyz); // cross(ij, ik)
        float n123_mag2 = dot(n123, n123);
        float il123 = 0.0f; // 1/|n123|
        if (n123_mag2 > 1e-16f) {
             il123 = rsqrt(n123_mag2);
             n123 *= il123; // normalize
        } else {
             // Handle degenerate case if needed
             n123 = (float3){0.0f,0.0f,0.0f};
        }

        // --- energy and force based on angle between n123 and il (-q41.xyz)
        // C++ uses: s = -n123.dot(q41.f); c = sqrt(1-s*s); E = K*(par.y + par.z*c + par.w*cs2.x)
        // Assuming C++ c = cos(w), s = sin(w) where w is the inversion angle.
        float s_w = -dot(n123, q41.xyz); // sin(w) : angle between normal(ijk) and vector il
        s_w = clamp(s_w, -1.0f, 1.0f);
        float c_w = sqrt(1.0f - s_w*s_w); // cos(w) (Should be >= 0 based on geometry?)

        float K  = par.x;
        float c0 = par.y;
        float c1 = par.z;
        float c2 = par.w;

        // E = K*(c0 + c1*cos(w) + c2*cos(2w))
        float c_2w = 2.0f*c_w*c_w - 1.0f; // cos(2w)
        E = K * ( c0 + c1*c_w + c2*c_2w );

        // Force calculation based on C++ Prokop code:
        // dE/dw = -K * ( c1*sin(w) + 2*c2*sin(2w) )
        float s_2w = 2.0f*s_w*c_w; // sin(2w)
        float dEdw = -K * ( c1*s_w + 2.0f*c2*s_2w );

        // const float f = -par.x * ( par.z * s + 2.0 * par.w * cs2.y ) / c; -> dEdw / c_w ?
        float f_term = (c_w > 1e-7f) ? dEdw / c_w : 0.0f; // Need check for c_w near zero

        float fq41 = f_term * q41.w; // f_term / l_li
        float fi123 = f_term * il123; // f_term / |n123|

        // Forces based on C++ evalInversion_Prokop
        float3 tq = s_w*fi123*n123 + fi123*q41.xyz; // Check sign of q41.xyz term? C++ uses +q41.f
        fp4 = fq41*n123 + s_w*fq41*q41.xyz;       // Check sign of q41.xyz term? C++ uses +q41.f
        fp2 = cross( q31.xyz, tq) * q21.w;        // CPU: fp2 = cross(q31, tq) * q21.e
        fp3 = cross( tq, q21.xyz) * q31.w;        // CPU: fp3 = cross(tq, q21) * q31.e
        fp1 = -(fp2 + fp3 + fp4);                 // Force on central atom i(1)
    } // End of inlined evalInversionUFF block

    if ((DBG_UFF!=0) && ii==IDBG_INV && (isys==IDBG_SYS)){
        int ia0=invAtoms[(i0I+ii)*4+0], ja0=invAtoms[(i0I+ii)*4+1], ka0=invAtoms[(i0I+ii)*4+2], la0=invAtoms[(i0I+ii)*4+3];
        float3 n123_dbg = cross(-q21.xyz, -q31.xyz);
        float n1 = length(n123_dbg);
        float s_w = (n1>1e-12f)? -dot(n123_dbg*(1.0f/n1), q41.xyz) : 0.0f; s_w=clamp(s_w,-1.0f,1.0f);
        float w = asin(s_w);
        printf("GPU INV %4d : ia=%4d ja=%4d ka=%4d la=%4d  K=% .9e c0=% .6e c1=% .6e c2=% .6e  w=% .9e  fi=(% .9e % .9e % .9e)  fj=(% .9e % .9e % .9e)  fk=(% .9e % .9e % .9e)  fl=(% .9e % .9e % .9e)  E=% .9e\n",
               ii, ia0,ja0,ka0,la0,
               par.x,par.y,par.z,par.w,
               w,
               fp1.x,fp1.y,fp1.z,
               fp2.x,fp2.y,fp2.z,
               fp3.x,fp3.y,fp3.z,
               fp4.x,fp4.y,fp4.z,
               E);
    }

    // --- Store forces into `fint` array ---
    float E_contrib = E * 0.25f;
    int idx_base = i0f + i0inv + ii * 4; // per-system fint offset
    fint[idx_base + 0] = (float4)(fp1, E_contrib); // Store force on atom i (central)
    fint[idx_base + 1] = (float4)(fp2, E_contrib); // Store force on atom j
    fint[idx_base + 2] = (float4)(fp3, E_contrib); // Store force on atom k
    fint[idx_base + 3] = (float4)(fp4, E_contrib); // Store force on atom l

    if (Ei_contrib) Ei_contrib[ii] = E; // Store per-inversion energy
}


// --- Kernel 5: Assemble Forces (Atom-Centric) ---
// Assumes C++ `a2f`-like map is provided.
// NOW needs to handle the directly accumulated bond forces as well.
__kernel void assembleForces_UFF(
    const int natoms,
    // --- Input force pieces ---
    __global float4* fint,         // Input: [nf] angle/dihedral/inversion force pieces {fx,fy,fz, E_contrib}
    // --- Assembly map (like C++ a2f) for fint pieces ---
    __global int*    a2f_offsets, // [natoms] Start index in a2f_indices
    __global int*    a2f_counts,  // [natoms] Number of fint entries for this atom
    __global int*    a2f_indices, // [total_fint_refs] Index into fint array for ang/dih/inv forces
    // --- Bond forces (if calculated directly) ---
    __global float4* fapos,       // Input buffer containing accumulated bond forces {fx,fy,fz,0} from Kernel 1
                                  // OR pass the temporary `fbond_on_atom` buffer if used.
                                  // Let's assume `fapos` was used by Kernel 1.
    // --- Output forces and energy ---
    // Output forces are written back to `fapos`
    const int        bClearForce,  // 1 to set fapos=sum(bonds+fint), 0 to do fapos+=sum(fint)
    const int        nf_per_system
                                   // If Kernel 1 already added bonds, use bClearForce=0 for fint part.
) {
    int ia   = get_global_id(0);
    int isys = get_global_id(1);
    // Per-system base offsets
    int i0a      = isys * natoms;          // per-system base for atom arrays (offsets/counts)
    int i0a2f    = isys * nf_per_system;   // per-system base for a2f_indices (length = nf_per_system)
    int i0f      = isys * nf_per_system;   // per-system base for fint
    // NOTE: a2f_indices values are per-system local indices into fint; add i0f below if needed
    // Debug dump by a single thread to avoid interleaved prints (selected system only)
    if ( (DBG_UFF!=0) &&  (ia == 0) && (isys==IDBG_SYS) ){
        printf("GPU assembleForces_UFF() ia=%3d natoms=%3d isys=%d\n", ia, natoms, isys);
        printf("GPU A2F TABLE natoms=%3d\n", natoms);
        for (int ia0 = 0; ia0 < natoms; ++ia0) {
            int off = a2f_offsets[i0a + ia0];
            int cnt = a2f_counts [i0a + ia0];
            float4 f = fapos[i0a + ia0];
            printf("GPU A2F ia=%3d f=( % .6e % .6e % .6e) off=%6d cnt=%4d idxs:", ia0, f.x, f.y, f.z, off, cnt);
            for (int k = 0; k < cnt; ++k) { int j = a2f_indices[i0a2f + off + k]; printf(" %d", j); }
            printf("\n");
        }
        printf("GPU FINT BY ATOM natoms=%3d\n", natoms);
        for (int ia0 = 0; ia0 < natoms; ++ia0) {
            int off = a2f_offsets[i0a + ia0];
            int cnt = a2f_counts [i0a + ia0];
            printf("GPU FINT ia=%3d:", ia0);
            for (int k = 0; k < cnt; ++k) {
                int j = a2f_indices[i0a2f + off + k];
                if(j>=0){ float4 v = fint[i0f + j]; printf(" [%3d]( % .6e % .6e % .6e)", j, v.x, v.y, v.z); }
                else{ printf(" [%3d]( skip )", j); }
            }
            printf("\n");
        }
    }
    if (ia >= natoms) return;

    float3 f_local = (float3){0.0f,0.0f,0.0f}; // Accumulate forces from fint for this atom
    float  E_local = 0.0f;         // Accumulate energy from fint for this atom

    int i0 = a2f_offsets[i0a + ia];
    int n  = a2f_counts [i0a + ia];
    int i1 = i0 + n;

    // Loop through all fint force pieces contributing to this atom (angles, dihedrals, inversions)
    for (int i = i0; i < i1; ++i) {
        int j = a2f_indices[i0a2f + i]; // per-system fint index (local)
        if(j<0) continue; // guard against invalid indices
        float4 force_piece = fint[i0f + j];
        f_local += force_piece.xyz;
        E_local += force_piece.w; // Accumulate energy contribution
    }

    // --- Combine with bond forces and write final force ---
    if (bClearForce != 0) {
        // Read bond force accumulated by Kernel 1 (assuming it's in fapos)
        float4 bond_force = fapos[i0a + ia];
        // Final force = bond_force + fint_forces
        fapos[i0a + ia] = (float4)(bond_force.xyz + f_local, E_local); // Store total E in .w ? Or keep bond energy separate? C++ doesn't store E in fapos.
                                                                 // Let's store total E contrib in .w here.
    } else {
        // Accumulate angle/dih/inv forces onto existing forces (which should include bonds from Kernel 1)
        fapos[i0a + ia] += (float4)(f_local, E_local);
    }
}



// Assemble recoil forces from neighbors and  update atoms positions and velocities
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void updateAtomsMMFFf4(
    const int4        n,            // 1 // (natoms,nnode) dimensions of the system
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  avel,         // 3 // velocities of atoms
    __global float4*  aforce,       // 4 // forces on atoms
    __global float4*  cvfs,         // 5 // damping coefficients for velocity and force
    __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
    __global float4*  TDrives       // 11 // Thermal driving (T,gamma_damp,seed,?)
    //__global cl_Mat3* bboxes,     // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
){
    const int natoms=n.x;           // number of atoms
    const int iG = get_global_id  (0); // index of atom

    if(iG>=natoms) return;

    const int iS = get_global_id  (1); // index of system
    const int nG = get_global_size(0); // number of atoms
    const int nS = get_global_size(1); // number of systems

    //const int ian = iG + iS*nnode;
    const int iaa = iG + iS*natoms;  // index of atom in atoms array

    const float4 MDpars  = MDparams[iS]; // (dt,damp,Flimit)
    const float4 TDrive = TDrives[iS];

    if((DBG_UFF>0) && (iS==IDBG_SYS)&&(iG==IDBG_ATOM)){
        // printf("GPU updateAtomsMMFFf4[isys=%i]: MDpars(dt=%g,Flim=%g,vel_damp_factor=%g,?= %g) \n", iS, MDpars.x,MDpars.y,MDpars.z,MDpars.w);
        // for(int is=0; is<nS; is++){
        //     //printf( "GPU::TDrives[%i](%g,%g,%g,%g)\n", i, TDrives[i].x,TDrives[i].y,TDrives[i].z,TDrives[i].w );
        //     //printf( "GPU::bboxes[%i](%g,%g,%g)(%g,%g,%g)(%g,%g,%g)\n", is, bboxes[is].a.x,bboxes[is].a.y,bboxes[is].a.z,   bboxes[is].b.x,bboxes[is].b.y,bboxes[is].b.z,   bboxes[is].c.x,bboxes[is].c.y,bboxes[is].c.z );
        //     for(int ia=0; ia<natoms; ia++){
        //         int ic = ia+is*natoms;
        //         if(constr[ia+is*natoms].w>0) printf( "GPU:sys[%i]atom[%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g|%g)\n", is, ia, constr[ic].x,constr[ic].y,constr[ic].z,constr[ic].w,   constrK[ic].x,constrK[ic].y,constrK[ic].z,constrK[ic].w  );
        //     }
        // }
    }

    const int iS_DBG = 5; // debug system
    //const int iG_DBG = 0;
    const int iG_DBG = 1; // debug atom

    //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "updateAtomsMMFFf4() natoms=%i nnode=%i nvec=%i nG %i iS %i/%i  dt=%g damp=%g Flimit=%g \n", natoms,nnode, nvec, iS, nG, nS, MDpars.x, MDpars.y, MDpars.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     int i0a = iS*natoms;
    //     for(int i=0; i<natoms; i++){
    //         printf( "GPU:constr[%i](%7.3f,%7.3f,%7.3f |K= %7.3f) \n", i, constr[i0a+i].x,constr[i0a+i].y,constr[i0a+i].z,  constr[i0a+i].w   );
    //     }
    // }

    if(iG>=(natoms)) return; // make sure we are not out of bounds of current system

    //aforce[iav] = float4Zero;

    float4 fe      = aforce[iaa]; // force on atom or pi-orbital

    // ------ Gather Forces from back-neighbors
    // ---- Limit Forces
    float Flimit = MDpars.y;
    float fr2 = dot(fe.xyz,fe.xyz);
    if( fr2 > (Flimit*Flimit) ){  fe.xyz *= (Flimit/sqrt(fr2)); }

    // =============== FORCE DONE
    aforce[iaa] = fe;             // store force before limit
    //aforce[iav] = float4Zero;   // clean force   : This can be done in the first forcefield run (best is NBFF)

    // =============== DYNAMICS

    float4 ve = avel[iaa]; // velocity of atom or pi-orbital
    float4 pe = apos[iaa]; // position of atom or pi-orbital

    // -------- Fixed Atoms and Bounding Box
    // if(iG<natoms){                  // only atoms have constraints, not pi-orbitals
    //     // ------- bboxes
    //     // const cl_Mat3 B = bboxes[iS];
    //     // if(B.c.z>0.0f){ if(pe.z<B.a.z){ fe.z+=(B.a.z-pe.z)*B.c.z; }else if(pe.z>B.b.z){ fe.z+=(B.b.z-pe.z)*B.c.z; }; }
    //     // ------- constrains
    //     float4 cons = constr[ iaa ]; // constraints (x,y,z,K)
    //     if( cons.w>0.f ){            // if stiffness is positive, we have constraint
    //         float4 cK = constrK[ iaa ];
    //         cK = max( cK, (float4){0.0f,0.0f,0.0f,0.0f} );
    //         const float3 fc = (cons.xyz - pe.xyz)*cK.xyz;
    //         fe.xyz += fc; // add constraint force
    //         if(iS==0){printf( "GPU::constr[ia=%i|iS=%i] (%g,%g,%g|K=%g) fc(%g,%g,%g) cK(%g,%g,%g)\n", iG, iS, cons.x,cons.y,cons.z,cons.w, fc.x,fc.y,fc.z , cK.x, cK.y, cK.z ); }
    //     }
    // }

    // Thermal driving  - Langevin thermostat, see C++ MMFFsp3_loc::move_atom_Langevin()
    if( TDrive.y > 0.0f ){ // if gamma>0
        fe.xyz    += ve.xyz * -TDrive.y ;  // damping,  check the untis  ... cdamp/dt = gamma
        //const float3 rnd = (float3){ hashf_wang(ve.x+TDrive.w,-1.0,1.0),hashf_wang(ve.y+TDrive.w,-1.0,1.0),hashf_wang(ve.z+TDrive.w,-1.0,1.0)};
        __private float3 ix;
        // + (float3){TDrive.w,TDrive.w,TDrive.w}
        //const float4 rnd = fract( (ve*541547.1547987f + TDrive.wwww), &ix )*2.f - (float4){1.0,1.0,1.0,1.0};  // changes every frame
        const float3 rvec = (float3){  // random vector depending on the index
            (((iG+136  + (int)(1000.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f,
            (((iG+778  + (int)(1013.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f,
            (((iG+4578 + (int)( 998.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f
        };
        //const float3 rnd = fract( ( rvec + TDrive.www)*12.4565f, &ix )*2.f - (float3){1.0,1.0,1.0};
        const float3 rnd = sin( ( rvec + TDrive.www )*124.4565f );
        //if(iS==3){  printf( "atom[%i] seed=%g rvec(%g,%g,%g) rnd(%g,%g,%g) \n", iG, TDrive.w, rvec.x,rvec.y,rvec.z, rnd.x,rnd.y,rnd.z ); }
        fe.xyz    += rnd.xyz * sqrt( 2*const_kB*TDrive.x*TDrive.y/MDpars.x );
    }

    float4 cvf = (float4){ dot(fe.xyz,fe.xyz),dot(ve.xyz,ve.xyz),dot(fe.xyz,ve.xyz), 0.0f };    // accumulate |f|^2 , |v|^2  and  <f|v>  to calculate damping coefficients for FIRE algorithm outside of this kernel
    cvfs[iaa] += cvf;
    //if(!bDrive){ ve.xyz *= MDpars.z; } // friction, velocity damping

    if ((DBG_UFF>0) && (iS==IDBG_SYS) && (iG==IDBG_ATOM)){
        printf( "GPU isys=%i iG=%i pos(%10.4e,%10.4e,%10.4e) vel(%10.4e,%10.4e,%10.4e) fe(%10.4e,%10.4e,%10.4e) cvf(%10.4e,%10.4e,%10.4e) \n", iS, iG, pe.x,pe.y,pe.z, ve.x,ve.y,ve.z, fe.x,fe.y,fe.z, cvf.x,cvf.y,cvf.z );
    }

    ve.xyz *=        MDpars.z;      // friction, velocity damping
    ve.xyz += fe.xyz*MDpars.x;      // acceleration
    pe.xyz += ve.xyz*MDpars.x;      // move
    //ve     *= 0.99f;              // friction, velocity damping
    //ve.xyz += fe.xyz*0.1f;        // acceleration
    //pe.xyz += ve.xyz*0.1f;        // move
    pe.w=0;ve.w=0;    // This seems to be needed, not sure why ?????
    avel[iaa] = ve;   // store velocity
    apos[iaa] = pe;   // store position

}



// ======================================================================
//                           getNonBond()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond(
    const int4 ns,                  // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
    __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
    __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams     // 9 // Grid-Force-Field parameters
){

    // we use local memory to store atomic position and parameters to speed up calculation, the size of local buffers should be equal to local workgroup size
    //__local float4 LATOMS[2];
    //__local float4 LCLJS [2];
    //__local float4 LATOMS[4];
    //__local float4 LCLJS [4];
    //__local float4 LATOMS[8];
    //__local float4 LCLJS [8];
    //__local float4 LATOMS[16];
    //__local float4 LCLJS [16];
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LCLJS [32];   // local buffer for atom parameters
    //__local float4 LATOMS[64];
    //__local float4 LCLJS [64];
    //__local float4 LATOMS[128];
    //__local float4 LCLJS [128];

    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory

    const int natoms=ns.x;  // number of atoms
    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int iaa = iG + i0a; // index of atom in atoms array

    // if((iG==IDBG_ATOM)&&(iS==IDBG_SYS)){  printf( "GPU::getNonBond() natoms(%i) nS,nG,nL(%i,%i,%i) \n", natoms, nS,nG,nL ); }
    // if((iG==IDBG_ATOM)&&(iS==IDBG_SYS)){
    //    printf( "GPU::getNonBond() natoms(%i) nS,nG,nL(%i,%i,%i) \n", natoms, nS,nG,nL );
    //    for(int is=0; is<nS; is++){
    //         cl_Mat3 lvec = lvecs[is];
    //         printf( "GPU[is:%i]   lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  );
    //         for(int ia=0; ia<natoms; ia++){
    //             int i = ia+is*natoms;
    //             printf( "GPU atom [is:%i,ia:%i] p(% .3f,% .3f,% .3f,% .3f) REQ(% .3f,% .3f,% .3f,% .3f) neighs{%3i,%3i,%3i,%3i} \n", is,ia, atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z,REQKs[i].w, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w );
    //         }
    //    }
    // }


    //if( iL==0 ){ for(int i=0; i<nL; i++){  LATOMS[i]=(float4){ 10000.0, (float)i,(float)iG,(float)iS }; LCLJS[i]=(float4){ 20000.0, (float)i,(float)iG,(float)iS }; } }


    // NOTE: if(iG>=natoms) we are reading from invalid adress => last few processors produce crap, but that is not a problem
    //       importaint is that we do not write this crap to invalid address, so we put   if(iG<natoms){forces[iav]+=fe;} at the end
    //       we may also put these if(iG<natoms){ .. } around more things, but that will unnecessarily slow down other processors
    //       we need these processors with (iG>=natoms) to read remaining atoms to the local memory.

    //if(iG<natoms){
    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;  // PBC is used if any of the PBC dimensions is >0
    //const bool bPBC=false;

    const int4   ng    = neighs   [iaa];  // neighbors indices
    const int4   ngC   = neighCell[iaa];  // neighbors cell indices
    const float4 REQKi = REQKs    [iaa];  // non-bonded parameters
    const float3 posi  = atoms    [iaa].xyz; // position of atom
    const float  R2damp = GFFParams.x*GFFParams.x; // squared damping radius
    float4 fe          = float4Zero;  // force on atom

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors for this system

    //if(iG==0){ printf("GPU[iS=%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f} \n", iS, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );  }

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;   // shift of PBC image 0
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                      // shift of PBC image in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                      // shift of PBC image in the outer loop
    //}

    /*
    if((iG==iG_DBG)&&(iS==iS_DBG)){
        printf( "GPU::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) bPBC=%i nPBC(%i,%i,%i)\n", natoms,nnode,nvec, nS,nG,nL, bPBC, nPBC.x,nPBC.y,nPBC.z );
        for(int i=0; i<natoms; i++){
            printf( "GPU a[%i] ", i);
            printf( "p{%6.3f,%6.3f,%6.3f} ", atoms[i0v+i].x,atoms[i0v+i].y,atoms[i0v+i].z  );
            printf( "ng{%i,%i,%i,%i} ", neighs[i0a+i].x,neighs[i0a+i].y,neighs[i0a+i].z,neighs[i0a+i].w );
            printf( "ngC{%i,%i,%i,%i} ", neighCell[i0a+i].x,neighCell[i0a+i].y,neighCell[i0a+i].z,neighCell[i0a+i].w );
            printf( "\n");
        }
    }
    */

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){      // loop over all atoms in the system, by chunks of size of local memory
        const int i=j0+iL;              // index of atom in local memory
        if(i<natoms){                   // j0*nL may be larger than natoms, so we need to check if we are not reading from invalid address
            LATOMS[iL] = atoms [i+i0a]; // read atom position to local memory
            LCLJS [iL] = REQKs [i+i0a]; // read atom parameters to local memory
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( (ja!=iG) && (ja<natoms) ){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float4 aj = LATOMS[jl];    // read atom position   from local memory
                float4 REQK = mixREQ_arithmetic( REQKi, LCLJS [jl] );   // read atom parameters from local memory
                float3 dp       = aj.xyz - posi; // vector between atoms
                //if((iG==44)&&(iS==0))printf( "[i=%i,ja=%i/%i,j0=%i,jl=%i/%i][iG/nG/na %i/%i/%i] aj(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,ja,nG,j0,jl,nL,   iG,nG,natoms,   aj.x,aj.y,aj.z,aj.w,  REQK.x,REQK.y,REQK.z,REQK.w  );

                // Debug print for pairwise interaction - print all interactions for atom 0
                //if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){ printf("GPU pair [%i,%i] dp(% .6e,% .6e,% .6e) r % .6e REQi(% .6e,% .6e,% .6e) REQij(% .6e,% .6e,% .6e) \n", iG, ja, dp.x, dp.y, dp.z, length(dp), REQKi.x,REQKi.y,REQKi.z,  REQK.x,REQK.y,REQK.z); }

                //REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                //REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy
                //REQK = mixREQ_arithmetic( REQKi, REQK );

                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
                //if( (j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                //if( (ji==1)&&(iG==0) )printf( "2 non-bond[%i,%i] bBonded %i\n",iG,ji,bBonded );

                if(bPBC){         // ===== if PBC is used, we need to loop over all PBC images of the atom
                    int ipbc=0;   // index of PBC image
                    dp += shift0; // shift to PBC image 0
                    // Fixed PBC size
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "GPU ipbc %i(%i,%i) shift(%7.3g,%7.3g,%7.3g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z );
                            // Without these IF conditions if(bBonded) time of evaluation reduced from 61 [ms] to 51[ms]
                           if( !( bBonded &&(                     // if atoms are bonded, we do not want to calculate non-covalent interaction between them
                                   ((ja==ng.x)&&(ipbc==ngC.x))    // check if this PBC image is not the same as one of the bonded atoms
                                   ||((ja==ng.y)&&(ipbc==ngC.y))  // i.e. if ja is neighbor of iG, and ipbc is its neighbor cell index then we skip this interaction
                                   ||((ja==ng.z)&&(ipbc==ngC.z))
                                   ||((ja==ng.w)&&(ipbc==ngC.w))
                           ))){
                                //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "ipbc %i(%i,%i) shift(%g,%g,%g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z );
                                float4 fij = getLJQH( dp, REQK, R2damp );  // calculate non-bonded force between atoms using LJQH potential
                                //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); }

                                // Debug print for force calculation - print all interactions for atom 0
                                // if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){ printf("GPU fij(% .6e,% .6e,% .6e) E % .6e REQij(% .6e,% .6e,% .6e) bBonded %i bPBC %i\n", fij.x, fij.y, fij.z, fij.w, REQK.x, REQK.y, REQK.z, bBonded, bPBC);  }

                                fe += fij;
                            }
                            ipbc++;
                            dp    += lvec.a.xyz;
                        }
                        dp    += shift_a;
                    }
                }else                                                  // ===== if PBC is not used, it is much simpler
                if( !bBonded ){
                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate non-bonded force between atoms using LJQH potential

                    // Debug print for force calculation (non-PBC case) - print all interactions for atom 0
                    // if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){
                    //     printf("GPU fij [i:%3i,j:%3i|isys:%i] dp(% .6e,% .6e,% .6e) r % .6e REQi(% .6e,% .6e,% .6e) REQij(% .6e,% .6e,% .6e) fij(% .6e,% .6e,% .6e) E % .6e bBonded %i bPBC %i\n", iG, ja, iS, dp.x, dp.y, dp.z, length(dp),  REQKi.x, REQKi.y, REQKi.z, REQK.x, REQK.y, REQK.z,  fij.x, fij.y, fij.z, fij.w, bBonded, bPBC );
                    // }

                    fe += fij;
                }  // if atoms are not bonded, we calculate non-bonded interaction between them
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(iG<natoms){
        //if(iS==0){ printf( "GPU::getNonBond(iG=%i) fe(%g,%g,%g,%g)\n", iG, fe.x,fe.y,fe.z,fe.w ); }
        forces[iaa] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //forces[iav] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //forces[iav] = fe*(-1.f);
    }
}





inline int modulo(const int i, const int m) {
    int result = i % m;
    if (result < 0) {
        result += m;
    }
    return result;
}

inline int4 make_inds_pbc(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int4)(0, 1,   2,   3  ); }
        case 1: { return (int4)(0, 1,   2,   3-n); }
        case 2: { return (int4)(0, 1,   2-n, 3-n); }
        case 3: { return (int4)(0, 1-n, 2-n, 3-n); }
    }
    return (int4)(-100, -100, -100, -100);
    // iqs[0] = (int4)(0, 1,   2,   3  );
    // iqs[1] = (int4)(0, 1,   2,   3-n);
    // iqs[2] = (int4)(0, 1,   2-n, 3-n);
    // iqs[3] = (int4)(0, 1-n, 2-n, 3-n);
}

inline int4 choose_inds_pbc(const int i, const int n, const int4* iqs) {
    if (i >= (n-3)) {
        const int ii = i + 4 - n;
        return iqs[ii];
    }
    return (int4)(0, +1, +2, +3);
}

inline int4 choose_inds_pbc_3( const int i, const int n, const int4* iqs ){
    if(i>=(n-3)){
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        const int4 d = iqs[ii];
        return (int4){ i+d.x, i+d.y, i+d.z, i+d.w };
    }
    return (int4){ i, i+1, i+2, i+3 };
}


inline float4 basis(float u) {
    const float inv6 = 1.0f / 6.0f;
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        inv6 * t * t * t,
        inv6 * (3.0f * u2 * (u - 2.0f) + 4.0f),
        inv6 * (3.0f * u * (1.0f + u - u2) + 1.0f),
        inv6 * u2 * u
    );
}

inline float4 dbasis(float u) {
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        -0.5f * t * t,
        0.5f * (3.0f * u2 - 4.0f * u),
        0.5f * (-3.0f * u2 + 2.0f * u + 1.0f),
        0.5f * u2
    );
}

inline float2 fe1Dcomb(__global const float4* E, const float4 C, const float4 p, const float4 d) {
    const float4 cs = (float4)(dot(C, E[0]), dot(C, E[1]), dot(C, E[2]), dot(C, E[3]));
    return (float2)(dot(p, cs), dot(d, cs));
}

inline float3 fe2d_comb(int nz, __global const float4* E, int4 di, const float4 C, const float4 pz, const float4 dz, const float4 by, const float4 dy) {
    const float2 fe0 = fe1Dcomb(E + di.x, C, pz, dz);
    const float2 fe1 = fe1Dcomb(E + di.y, C, pz, dz);
    const float2 fe2 = fe1Dcomb(E + di.z, C, pz, dz);
    const float2 fe3 = fe1Dcomb(E + di.w, C, pz, dz);

    return (float3)(
        fe0.x * dy.x + fe1.x * dy.y + fe2.x * dy.z + fe3.x * dy.w,
        fe0.y * by.x + fe1.y * by.y + fe2.y * by.z + fe3.y * by.w,
        fe0.x * by.x + fe1.x * by.y + fe2.x * by.z + fe3.x * by.w
    );
}

inline float4 fe3d_pbc_comb(const float3 u, const int3 n, __global const float4* Es, const float4 PLQH, __local const int4* xqis, __local int4* yqis) {
    int ix = (int)u.x;
    int iy = (int)u.y;
    int iz = (int)u.z;
    if (u.x < 0) ix--;
    if (u.y < 0) iy--;
    const float tx = u.x - ix;
    const float ty = u.y - iy;
    const float tz = u.z - iz;

    if ((iz < 1) || (iz >= n.z - 2)) {
        return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }

    ix = modulo(ix-1, n.x);
    iy = modulo(iy-1, n.y);

    const int nyz = n.z * n.y;
    // int4 qx = xqis[ix%4] * nyz;
    // int4 qy = yqis[iy%4] * n.z;

    int4 qx = choose_inds_pbc( ix, n.x, xqis );
    //const int4 qx = choose_inds_pbc( ix, n.x, xqis )*nyz;
    const int4 qy = choose_inds_pbc( iy, n.y, yqis )*n.z;

    const float4 bz = basis(tz);
    const float4 dz = dbasis(tz);
    const float4 by = basis(ty);
    const float4 dy = dbasis(ty);

    const int i0 = (iz - 1) + n.z * (iy + n.y * ix);

    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) n(%i,%i,%i) \n", u.x,u.y,u.z, ix,iy,iz, n.x,n.y,n.z );
    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) qx(%i,%i,%i,%i) nyz=%i\n", u.x,u.y,u.z, ix,iy,iz, qx.x,qx.y,qx.z,qx.w, nyz );
    qx*=nyz;

    //return (float4){ 0.0f, 0.0f, 0.0f, dot(PLQH, Es[ i0 ])  };

    float3 E1 = fe2d_comb(n.z, Es + (i0 + qx.x), qy, PLQH, bz, dz, by, dy);
    float3 E2 = fe2d_comb(n.z, Es + (i0 + qx.y), qy, PLQH, bz, dz, by, dy);
    float3 E3 = fe2d_comb(n.z, Es + (i0 + qx.z), qy, PLQH, bz, dz, by, dy);
    float3 E4 = fe2d_comb(n.z, Es + (i0 + qx.w), qy, PLQH, bz, dz, by, dy);

    const float4 bx = basis(tx);
    const float4 dx = dbasis(tx);

    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))
    );
}

inline int gridIndex(int ix, int iy, int iz, int3 nxyz) {
    return iz + nxyz.z * (iy + nxyz.y * ix);
}

// ======================================================================
//                           getNonBond_GridFF_Bspline()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system, and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by texture interpolation
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
// NOTE: this modified version of getNonBond() just by including Grid-Force-Field (GFF) forces
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF_Bspline(
    const int4 ns,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW,EvdW,Q,H)
    // GridFF
    __global float4*  BsplinePLQ,   // 10 // Grid-Force-Field (GFF) for Pauli repulsion
    const int4     grid_ns,         // 11 // origin of the grid
    const float4   grid_invStep,    // 12 // origin of the grid
    const float4   grid_p0          // 13 // origin of the grid
){
    __local float4 LATOMS[32];         // local memory chumk of positions of atoms
    __local float4 LCLJS [32];         // local memory chumk of atom parameters
    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms=ns.x;         // number of atoms in the system

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system

    const float4 REQKi = REQKs    [iaa];           // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H) of the atom
    const float3 posi  = atoms    [iaa].xyz;       // position of the atom
    float4 fe          = float4Zero;              // forces on the atom


    // =================== Non-Bonded interaction ( molecule-molecule )

    if(ns.w>=0)
    { // insulate nbff

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors of the system

    // if((iG==IDBG_ATOM)&&(iS==IDBG_SYS)){  }
    // if((iG==IDBG_ATOM)&&(iS==IDBG_SYS)){   }
    // if((iG==iG_DBG)&&(iS==iS_DBG)){  }
    // if((iG==iG_DBG)&&(iS==iS_DBG)) 
    if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){
        printf( "GPU::getNonBond_GridFF_Bspline() natoms(%i) nS,nG,nL(%i,%i,%i) \n", natoms, nS,nG,nL );
        printf( "GPU::getNonBond_GridFF_Bspline() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
        printf( "GPU::getNonBond_GridFF_Bspline() grid_p0(%12.4f,%12.4f,%12.4f) grid_invStep(%12.4f,%12.4f,%12.4f) \n", grid_p0.x,grid_p0.y,grid_p0.z, grid_invStep.x, grid_invStep.y,grid_invStep.z );
        printf( "GPU::getNonBond_GridFF_Bspline() grid_ns(%4i,%4i,%4i|%6i) GFFParams(%12.4f,%12.4f,%12.4f,%12.4f) \n", grid_ns.x,grid_ns.y,grid_ns.z,grid_ns.w, GFFParams.x,GFFParams.y,GFFParams.z,GFFParams.w );
        float4 gp;
        gp=BsplinePLQ[gridIndex( 0, 0, 0, grid_ns.xyz)]; printf( "GPU:BsplinePLQ[0,0,0] (%12.4f,%12.4f,%12.4f,%12.4f) \n", gp.x,gp.y,gp.z,gp.w );
        gp=BsplinePLQ[gridIndex( 0,10,10, grid_ns.xyz)]; printf( "GPU:BsplinePLQ[1,1,1] (%12.4f,%12.4f,%12.4f,%12.4f) \n", gp.x,gp.y,gp.z,gp.w );
        gp=BsplinePLQ[gridIndex(10,10,10, grid_ns.xyz)]; printf( "GPU:BsplinePLQ[2,2,2] (%12.4f,%12.4f,%12.4f,%12.4f) \n", gp.x,gp.y,gp.z,gp.w ); 
        // printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL );
        // for(int i=0; i<nS*nG; i++){
        //     int ia = i%nS;
        //     int is = i/nS;
        //     if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
        //     //printf( "GPU[%i,%i] \n", is,ia,  );
        // }
        for(int ia=0; ia<natoms; ia++){
            float4 pos = atoms[ia+i0a];
            float4 req = REQKs[ia+i0a];
             printf( "GPU[%3i,%3i] pos(%12.4f,%12.4f,%12.4f) REQ(%12.4f,%12.4f,%12.4f,%12.4f) \n", iS,ia, pos.x,pos.y,pos.z, req.x,req.y,req.z,req.w );
        }
    }

    //if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0; // Periodic boundary conditions if any of nPBC.x,nPBC.y,nPBC.z is non-zero
    const int4   ng    = neighs   [iaa];           // indexes of neighbors of the atom
    const int4   ngC   = neighCell[iaa];           // indexes of cells of neighbors of the atom

    const float  R2damp = GFFParams.x*GFFParams.x; // damping radius for Lenard-Jones potential

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;  // shift of the first PBC image
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                  // shift of lattice vector in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                  // shift of lattice vector in the outer loop

    // ========= Atom-to-Atom interaction ( N-body problem )     - we do it by chunks of nL atoms in order to reuse data and reduce number of global memory reads
    for (int j0=0; j0<natoms; j0+= nL ){ // loop over atoms in the system by chunks of nL atoms which fit into local memory
        const int i = j0 + iL;           // global index of atom in the system
        LATOMS[iL] = atoms [i+i0a];      // load positions  of atoms into local memory
        LCLJS [iL] = REQKs [i+i0a];      // load parameters of atoms into local memory
        barrier(CLK_LOCAL_MEM_FENCE);    // wait until all atoms are loaded into local memory
        for (int jl=0; jl<nL; jl++){     // loop over atoms in the local memory chunk
            const int ja=jl+j0;          // global index of atom in the system
            if( (ja!=iG) && (ja<natoms) ){ // atom should not interact with himself, and should be in the system ( j0*nL+iL may be out of range of natoms )
                const float4 aj   = LATOMS[jl]; // position of the atom
                float4       REQK = LCLJS [jl]; // parameters of the atom
                float3 dp   = aj.xyz - posi;    // vector between atoms
                REQK = mixREQ_arithmetic( REQKi, REQK ); // Correct mixing
                
                // Debug print for pairwise interaction - print all interactions for atom 0
                //if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){ printf("GPU pair [%i,%i] dp(% .6e,% .6e,% .6e) r % .6e REQi(% .6e,% .6e,% .6e) REQij(% .6e,% .6e,% .6e) \n", iG, ja, dp.x, dp.y, dp.z, length(dp), REQKi.x,REQKi.y,REQKi.z,  REQK.x,REQK.y,REQK.z); }
                
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));   // atom is bonded if it is one of the neighbors
                if(bPBC){       // ==== with periodic boundary conditions we need to consider all PBC images of the atom
                    int ipbc=0; // index of PBC image
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                    dp -= shift0;  // shift to the first PBC image
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(  // if bonded in any of PBC images, then we have to check both index of atom and index of PBC image to decide if we should skip this interaction
                                          ((ja==ng.x)&&(ipbc==ngC.x)) // check 1. neighbor and its PBC cell
                                        ||((ja==ng.y)&&(ipbc==ngC.y)) // check 2. neighbor and its PBC cell
                                        ||((ja==ng.z)&&(ipbc==ngC.z)) // ...
                                        ||((ja==ng.w)&&(ipbc==ngC.w)) // ...
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); }
                                    
                                    // Debug print for force calculation - print all interactions for atom 0
                                    // if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){ printf("GPU fij(% .6e,% .6e,% .6e) E % .6e REQij(% .6e,% .6e,% .6e) bBonded %i bPBC %i\n", fij.x, fij.y, fij.z, fij.w, REQK.x, REQK.y, REQK.z, bBonded, bPBC);  }
                                    
                                    fe += fij; // accumulate forces
                                }
                                ipbc++;         // increment index of PBC image
                                dp+=lvec.a.xyz; // shift to the next PBC image
                            }
                            dp+=shift_a;        // shift to the next PBC image
                        }
                        dp+=shift_b;            // shift to the next PBC image
                    }
                }else{ //  ==== without periodic boundary it is much simpler, not need to care about PBC images
                    if(bBonded) continue;  // Bonded ?
                    float4 fij = getLJQH( dp, REQK, R2damp ); // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                    
                    // Debug print for force calculation (non-PBC case) - print all interactions for atom 0
                    //if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){
                    //    printf("GPU fij [i:%3i,j:%3i|isys:%i] dp(% .6e,% .6e,% .6e) r % .6e REQi(% .6e,% .6e,% .6e) REQij(% .6e,% .6e,% .6e) fij(% .6e,% .6e,% .6e) E % .6e bBonded %i bPBC %i\n", iG, ja, iS, dp.x, dp.y, dp.z, length(dp),  REQKi.x, REQKi.y, REQKi.z, REQK.x, REQK.y, REQK.z,  fij.x, fij.y, fij.z, fij.w, bBonded, bPBC );
                    //}
                    
                    fe += fij;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // wait until all atoms are processed, ToDo: not sure if it is needed here ?
    }

    } // insulate nbff

    if(iG>=natoms) return; // natoms <= nG, because nG must be multiple of nL (loccal kernel size). We cannot put this check at the beginning of the kernel, because it will break reading of atoms to local memory

    __local int4 xqs[4];
    __local int4 yqs[4];
    { // insulate gridff
        if      (iL<4){             xqs[iL]=make_inds_pbc(grid_ns.x,iL); }
        else if (iL<8){ int i=iL-4; yqs[i ]=make_inds_pbc(grid_ns.y,i ); };
        //const float3 inv_dg = 1.0f / grid_d.xyz;
        barrier(CLK_LOCAL_MEM_FENCE);

        const float ej = exp( GFFParams.y * REQKi.x ); // exp(-alphaMorse*RvdW) pre-factor for factorized Morse potential
        const float4 PLQH = (float4){
            ej*ej*REQKi.y,                   // prefactor London dispersion (attractive part of Morse potential)
            ej*   REQKi.y,                   // prefactor Pauli repulsion   (repulsive part of Morse potential)
            REQKi.z,
            0.0f
        };
        //const float3 p = ps[iG].xyz;
        const float3 u = (posi - grid_p0.xyz) * grid_invStep.xyz;

        float4 fg = fe3d_pbc_comb(u, grid_ns.xyz, BsplinePLQ, PLQH, xqs, yqs);

        // GridFF-specific debug prints for test atom
    //    if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){
    //         printf("GPU GridFF [i:%3i|isys:%i] posi(% .6e,% .6e,% .6e) u(% .6e,% .6e,% .6e) PLQH(% .6e,% .6e,% .6e,% .6e) fg_raw(% .6e,% .6e,% .6e,% .6e)\n",
    //                iG, iS, posi.x, posi.y, posi.z, u.x, u.y, u.z, PLQH.x, PLQH.y, PLQH.z, PLQH.w, fg.x, fg.y, fg.z, fg.w);
    //    }

        // for(int ia=0; ia<natoms; ia++){
        //     float4 pos = atoms[ia+i0a];
        //     float4 req = REQKs[ia+i0a];
        //     printf( "GPU[%3i,%3i] pos(%12.4f,%12.4f,%12.4f) fg(%12.4f,%12.4f,%12.4f,%12.4f)\n", iS,ia, pos.x,pos.y,pos.z, fg.x,fg.y,fg.z,fg.w );
        // }

        //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() fg(%g,%g,%g|%g) u(%g,%g,%g) posi(%g,%g,%g) grid_invStep(%g,%g,%g)\n", fg.x,fg.y,fg.z,fg.w,  u.x,u.y,u.z, posi.x,posi.y,posi.z, grid_invStep.x, grid_invStep.y, grid_invStep.z  ); }

        fg.xyz *= -grid_invStep.xyz;
        
        // Debug print for final GridFF force after scaling
        // if((DBG_UFF!=0) && (iG==IDBG_ATOM) && (iS==IDBG_SYS)){
        //     printf("GPU GridFF scaled [i:%3i|isys:%i] fg_scaled(% .6e,% .6e,% .6e,% .6e) grid_invStep(% .6e,% .6e,% .6e)\n",
        //            iG, iS, fg.x, fg.y, fg.z, fg.w, grid_invStep.x, grid_invStep.y, grid_invStep.z);
        // }
        
        fe += fg;
        //fes[iG] = fe;
    }  // insulate gridff
    forces[iaa] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
    //forces[iav] += fe;     // If we don't run it as first forcefield, we need to add forces to the forces calculated by previous forcefields
    //forces[iav] = fe*(-1.f);


}




// ======================================================================
//                           getSurfMorse()
// ======================================================================
// This is brute-force alternative to getNonBond_GridFF_Bspline - it describe interactin of molecule with substrate by pairwise interactins with multiple replicas

__kernel void getSurfMorse(
    const int4 ns,                // 1
    __global float4*  atoms,      // 2
    __global float4*  REQs,       // 3
    __global float4*  forces,     // 4
    __global float4*  atoms_s,    // 5
    __global float4*  REQ_s,      // 6
    const int4     nPBC,          // 7
    const cl_Mat3  lvec,          // 8
    const float4   pos0,          // 9
    const float4   GFFParams      // 10
){

    __local float4 LATOMS[32];
    __local float4 LCLJS [32];

    const int nAtoms  = ns.x;

    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms  = ns.x;         // number of atoms in the system
    const int nnode   = ns.y;         // number of nodes in the system
    const int nvec    = natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system
    const int na_surf = ns.z;         //

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system

    if( (iG==IDBG_ATOM) && (iS==IDBG_SYS) ){
        printf("GPU::getSurfMorse() nglob(%i,%i) nloc(%i) ns(%i,%i,%i) \n", nG,nS, nL, ns.x,ns.y,ns.z  );
        printf("GPU::lvec.a=(%g,%g,%g) lvec.b=(%g,%g,%g) lvec.c=(%g,%g,%g)\n", lvec.a.x, lvec.a.y, lvec.a.z, lvec.b.x, lvec.b.y, lvec.b.z, lvec.c.x, lvec.c.y, lvec.c.z);
        printf("GPU::pos0=(%g,%g,%g) \n", pos0.x, pos0.y, pos0.z);
        printf("GPU::GFFParams=(%g,%g,%g,%g) \n", GFFParams.x, GFFParams.y, GFFParams.z, GFFParams.w);
        printf("GPU::nPBC=(%i,%i,%i) \n", nPBC.x, nPBC.y, nPBC.z);
        //printf("GPU::getSurfMorse() nglob(%i,%i) nloc(%i) ns(%i,%i,%i) nPBC(%i,%i,%i)\n", nG,nS, nL, ns.x,ns.y,ns.z,  nPBC.x,nPBC.y,nPBC.z  );
        //for(int i=0; i<ns.x; i++){ printf( "forces[%i] (%g,%g,%g) \n", i, forces[i].x,forces[i].y,forces[i].z );   }
        for(int i=0; i<na_surf; i++){ printf( "GPU.atoms_s[%i](%g,%g,%g) REQs[%i](%g,%g,%g,%g) \n", i, atoms_s[i].x,atoms_s[i].y,atoms_s[i].z, i, REQ_s[i].x,REQ_s[i].y,REQ_s[i].z,REQ_s[i].w );   }
    }
    float4 fe   = (float4){0.0f,0.0f,0.0f,0.0f};

    // ========== BEGIN test loading atoms to local memory

    if((iG>nL) || (iS>0)) return;

    // for (int j0=0; j0<na_surf; j0+= nL ){
    //     const int i = j0 + iL;
    //     LATOMS[iL] = atoms_s[i];
    //     LCLJS [iL] = REQ_s  [i];
    //     barrier(CLK_LOCAL_MEM_FENCE);
    //     if(iG==0){
    //         for(int j=0; j<nL; j++){
    //             int ji = j0 + j;
    //             printf("j: %i i: %i , LATOMS = (%16.6f,%16.6f,%16.6f), atoms_s = (%16.6f,%16.6f,%16.6f)\n", j, ji, LATOMS[j].x, LATOMS[j].y, LATOMS[j].z, atoms_s[ji].x, atoms_s[ji].y, atoms_s[ji].z);
    //         }
    //         //printf("iL = %i, LATOMS = (%16.6f,%16.6f,%16.6f), atoms_s = (%16.6f,%16.6f,%16.6f)\n", iL, LATOMS[iL].x, LATOMS[iL].y, LATOMS[iL].z, atoms_s[i].x, atoms_s[i].y, atoms_s[i].z);
    //     }
    
    // }


    // ========== END test loading atoms to local memory


    const float  K          = -GFFParams.y;
    const float  R2damp     =  1e-30;
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz)

    const float3 pos  = atoms[iav].xyz - pos0.xyz +  lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell
    const float4 REQi = REQs [iaa];
    // if((iaa<22)&&(iS==0)) printf("REQi[%i].z = %f\n", iaa, REQi.z);
    // ToDo: perhaps it is efficient to share surface along isys direction ( all system operate with the same surface atoms )
    // for (int j0 = 0; j0 < na_surf; j0 ++) {
    //     LATOMS[j0   ] = atoms_s[j0];
    //     LCLJS[j0] = REQ_s[j0];
    // }
    for (int j0=0; j0<na_surf; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms_s[i];
        LCLJS [iL] = REQ_s  [i];
    //        if(iS==0 && iL > 21)printf("iL = %i, LATOMS = (%f,%f,%f), atoms_s = (%f,%f,%f)\n", iL, LATOMS[iL].x, LATOMS[iL].y, LATOMS[iL].z, atoms_s[i].x, atoms_s[i].y, atoms_s[i].z);
        barrier(CLK_LOCAL_MEM_FENCE);
        if(iG>=nAtoms) continue;
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( ja<na_surf ){
                float4 REQH =       LCLJS [jl];
                float3 dp   = pos - LATOMS[jl].xyz;
                REQH.x   += REQi.x;
                REQH.yzw *= REQi.yzw;
                //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                //dp+=lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
                //float3 shift=shift0;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                            const float4 fej = getMorseQH( dp,  REQH, K, R2damp );
                            fe -= fej;
                            if( (iG==IDBG_ATOM) && (iS==IDBG_SYS) && (iz==0)&&(iy==0)&&(ix==0)){
                            //if( (iG==0) && (iS==0) && (jl==0) ){
                            // if( (iG==0) && (jl==0)   && (iz==0)&&(iy==0)&&(ix==0)  ){
                            //    printf( "GPU[%3i|%3i/%3i] dp(%10.6f,%10.6f,%10.6f) LATOMS[%i](%10.6f,%10.6f,%10.6f) pos(%10.6f,%10.6f,%10.6f) pos0(%10.6f,%10.6f,%10.6f)  \n", ja,0,0,   dp.x,dp.y,dp.z,  jl,  LATOMS[jl].x,LATOMS[jl].y,LATOMS[jl].z,   pos.x,pos.y,pos.z,  pos0.x,pos0.y,pos0.z );
                                printf( "GPU[%3i(%3i,%3i,%3i)] K,R2damp(%10.6f,%10.6f) l=%10.6f dp(%10.6f,%10.6f,%10.6f) REQij(%10.6f,%10.6f,%10.6f,%10.6f) fij(%10.6f,%10.6f,%10.6f) \n", ja, ix,iy,iz,  K,R2damp, length(dp), dp.x,dp.y,dp.z,  REQH.x, REQH.y, REQH.z, REQH.w,  fej.x,fej.y,fej.z );
                            }
                            dp   +=lvec.a.xyz;
                            //shift+=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                        //shift+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp   +=shift_c;
                    //shift+=shift_c;
                    //dp+=lvec.b.xyz*(nPBC.y*-2.f-1.f);
                    //dp+=lvec.c.xyz;
                }
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    if( (iG==IDBG_ATOM) && (iS==IDBG_SYS) ){
         printf( "GPU[iG=%i,iS=%i] fe(%10.6f,%10.6f,%10.6f)\n", iG,iS, fe.x,fe.y,fe.z );
    }

    // if( (iG==0) ){
    //       printf( "GPU[iG=%i,iS=%i] fe(%10.6f,%10.6f,%10.6f)\n", iG,iS, fe.x,fe.y,fe.z );
    // }
    forces[iav] += fe;

}
