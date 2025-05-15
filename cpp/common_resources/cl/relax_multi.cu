#ifndef MMFF_KERNELS_CUH
#define MMFF_KERNELS_CUH


#define iGdbg 0
#define iSdbg 0

#define NNEIGH 4

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <vector_types.h>
#include <math.h> // For sqrtf, expf, sinf, cosf, fractf, fmaxf, fminf
#include <stdio.h> // For device-side printf
#include "relax_multi.cuh"

// cu_Mat3 struct is defined in relax_multi.cuh
// Use standard CUDA vector literals or initializer lists
#define  float4Zero  make_float4(0.f,0.f,0.f,0.f)
#define  float3Zero  make_float3(0.f,0.f,0.f)
#define  float2Zero  make_float2(0.f,0.f)

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f       // [ eV*Ang/e^2 ]
#define const_kB        8.617333262145e-5f   // [ eV/K ]

// --- CUDA Vector Math Helpers ---
// Define swizzle macros
#define XYZ(v) make_float3((v).x,(v).y,(v).z)
#define YZ(v) make_float2((v).y,(v).z)
#define XY(v) make_float2((v).x,(v).y)
#define XZ(v) make_float2((v).x,(v).z)

// Dot product
__device__ inline float dot(float3 a, float3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
__device__ inline float length(float3 v) { return sqrtf(dot(v, v)); }
__device__ inline float3 normalize(float3 v) { float invLen = 1/length(v); return make_float3(v.x * invLen, v.y * invLen, v.z * invLen); }

// Scalar multiplication
__device__ inline float3  operator* (float3  v, float s  ){ return make_float3(v.x * s, v.y * s, v.z * s); }
__device__ inline float3  operator* (float   s, float3 v ){ return make_float3(v.x * s, v.y * s, v.z * s); }
__device__ inline float3& operator*=(float3& v, float s  ){ v.x *= s; v.y *= s; v.z *= s; return v;}
__device__ inline float3  operator/ (float3  v, float s  ){ float invS = 1.0f/s; return make_float3(v.x * invS, v.y * invS, v.z * invS); }
__device__ inline float3& operator/=(float3& v, float s  ){ float invS = 1.0f/s; v.x *= invS; v.y *= invS; v.z *= invS; return v; }

// Vector addition/subtraction
__device__ inline float3  operator+ (float3  a,       float3  b) { return make_float3(a.x + b.x, a.y + b.y, a.z + b.z); }
__device__ inline float3  operator- (float3  a,       float3  b) { return make_float3(a.x - b.x, a.y - b.y, a.z - b.z); }
__device__ inline float3  operator* (float3  a,       float3  b) { return make_float3(a.x * b.x, a.y * b.y, a.z * b.z); }

__device__ inline float3& operator+=(float3& a, const float3& b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
__device__ inline float3& operator-=(float3& a, const float3& b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
__device__ inline float4& operator+=(float4& a, const float4& b) { a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w; return a; }

// --- End Helper Functions ---

// Use __device__ for functions callable from kernels
__device__ inline float2 udiv_cmplx( float2 a, float2 b ){ return make_float2(a.x*b.x + a.y*b.y, a.y*b.x - a.x*b.y);  }
__device__ inline float3 rotMat    ( float3 v, float3 a, float3 b, float3 c ){  return make_float3(dot(v,a), dot(v,b), dot(v,c));  }
__device__ inline float3 rotMatT   ( float3 v, float3 a, float3 b, float3 c ){  return a*v.x + b*v.y + c*v.z;  }

__device__ inline float evalAngCos( float4 hr1, float4 hr2, float K, float c0, float3* f1, float3* f2 ){
    float3 h1 = make_float3(hr1.x,hr1.y,hr1.z);
    float3 h2 = make_float3(hr2.x,hr2.y,hr2.z);
    float  c = dot(h1,h2);
    float3 hf1 = h2 - h1*c;
    float3 hf2 = h1 - h2*c;
    float c_   = c-c0;
    float E    = K*c_*c_;
    float fang = -K*c_*2;
    hf1 *= fang*hr1.w;
    hf2 *= fang*hr2.w;
    *f1=hf1;
    *f2=hf2;
    return E;
}

__device__ inline float evalAngleCosHalf( float4 hr1, float4 hr2, float2 cs0, float k, float3* f1, float3* f2 ){
    float3 h  = make_float3(hr1.x,hr1.y,hr1.z) + make_float3(hr2.x,hr2.y,hr2.z);
    float  c2 = dot(h,h)*0.25f;
    float  s2 = 1.f-c2 + 1e-7f;
    float2 cso = make_float2(sqrtf(c2), sqrtf(s2));
    float2 cs = udiv_cmplx(cs0, cso);
    float  E  = k*(1 - cs.x);
    float  fr = -k*cs.y;
    c2 *= -2.f;
    fr /= 4.f*cso.x*cso.y;
    float fr1 = fr*hr1.w;
    float fr2 = fr*hr2.w;
    *f1 = h*fr1 + make_float3(hr1.x,hr1.y,hr1.z)*(fr1*c2);
    *f2 = h*fr2 + make_float3(hr2.x,hr2.y,hr2.z)*(fr2*c2);
    return E;
}

__device__ inline float evalPiAling( float3 h1, float3 h2, float K, float3* f1, float3* f2 ){
    float c = dot(h1,h2);
    float3 hf1 = h2 - h1*c;
    float3 hf2 = h1 - h2*c;
    bool sign = c<0; if(sign) c=-c;
    float E = -K*c;
    float fang = K;
    if(sign) fang=-fang;
    hf1 *= fang;
    hf2 *= fang;
    *f1=hf1;
    *f2=hf2;
    return E;
}

__device__ inline float evalBond( float3 h, float dl, float k, float3* f ){
    float fr = dl*k;
    *f = h * fr;
    return fr*dl*0.5f;
}

__device__ inline float4 getLJQH( float3 dp, float4 REQ, float R2damp ){
    float r2 = dot(dp,dp);
    float ir2_ = 1.f/(r2 + R2damp);
    float Ec = COULOMB_CONST*REQ.z*sqrtf(ir2_);
    float ir2 = 1.f/r2;
    float u2 = REQ.x*REQ.x*ir2;
    float u6 = u2*u2*u2;
    float vdW = u6*REQ.y;
    float E = (u6-2.f)*vdW + Ec;
    float fr = -12.f*(u6-1.f)*vdW*ir2 - Ec*ir2_;
    return make_float4(dp.x*fr, dp.y*fr, dp.z*fr, E);
}

// --- Placeholder definitions for missing functions ---
// IMPORTANT: Replace these with your actual CUDA implementations!
// __device__ inline float evalPiAling( float3 hpi, float3 hpi_ng, float kpp, float3* f1, float3* f2 ){ *f1=float3Zero; *f2=float3Zero; return 0.0f; }
// __device__ inline float evalAngCos( float4 hr1, float4 hr2, float K, float C0, float3* f1, float3* f2 ){ *f1=float3Zero; *f2=float3Zero; return 0.0f; }
// __device__ inline float evalAngleCosHalf( float4 hr1, float4 hr2, float2 cs0, float K, float3* f1, float3* f2 ){ *f1=float3Zero; *f2=float3Zero; return 0.0f; }
// __device__ inline float4 getLJQH( float3 dp, float4 REQij, float R2damp ){ return float4Zero; }
// __device__ inline float evalBond(float3 h_xyz, float l, float bK, float3* f1) {
//     float3 f1_xyz = h_xyz * bK * (l - 1.0f); // f = h_xyz * bK * (l - 1.0f)
//     *f1 = f1_xyz;
//     return 0.5f * bK * (l - 1.0f) * (l - 1.0f); // E = 0.5f * bK * (l - 1.0f) * (l - 1.0f)
// }
// --- End Placeholder definitions ---
// ======================================================================
// ======================================================================
//                           MMFF kernells
// ======================================================================
// ======================================================================


// ======================================================================
//                     cleanForceMMFFf4()
// ======================================================================
// Clean aforce on atoms/pi-orbitals and neighbors to prepare for next forcefield evaluation
// Based on analysis, the OpenCL version likely had incorrect fneigh indexing for clearing.
// This CUDA version clears fapos/aforce for ALL vectors (atoms + pi) and fneigh for ALL node recoil aforce (sigma + pi).
// Assuming launched over grid(nvec, nS)

__global__ void cleanForceMMFFf4(
    int4        n,           // 2 // (natoms,nnode,?,?)
    float4*  aforce,      // 5 // aforce on atoms and pi-orbitals
    float4*  fneigh       // 6 // recoil aforce on neighbors (sigma + pi)
){
    const int natoms = n.x;
    const int nnode  = n.y;
    const int nvec   = natoms+nnode;
    const int iG     = blockIdx.x * blockDim.x + threadIdx.x; // index of vector (atom or pi)
    const int iS     = blockIdx.y * blockDim.y + threadIdx.y; // index of system

    const int iav = iG + iS*nvec; // global index in aforce


    if( (iG==0) &&(iS==0) ){ printf("CUDA cleanForceMMFFf4(): natoms=%i, nnode=%i \n", natoms, nnode );}

    // Clear force on this vector (atom or pi-orbital)
    // This should be done for ALL vectors (0..nvec-1)
    if(iG < nvec) {
       aforce[iav] = float4Zero;
       // Optional: aforce[iav] = make_float4(__int_as_float(iG), __int_as_float(iS), __int_as_float(iav), 0.0f); // Debugging
    }

    // Clear recoil aforce in fneigh.
    // fneigh stores recoil aforce for node atoms ONLY.
    // This means only threads with iG < nnode need to clear fneigh entries.
    // The indices to clear are those written by getMMFFf4.
    // fneigh layout derived from getMMFFf4: [ iS * nnode*8 + iG*8 + sigma/pi*4 + neigh_idx ]
    if(iG < nnode){
        // Base index for sigma aforce of node iG in system iS: iS * nnode*8 + iG*8
        const int fneigh_node_sigma_base_idx = iS * nnode * 8 + iG * 8;
        // Base index for pi aforce of node iG in system iS: iS * nnode*8 + iG*8 + 4
        const int fneigh_node_pi_base_idx    = fneigh_node_sigma_base_idx + 4;

        // Clear sigma recoil aforce for this node's 4 neighbors
        fneigh[fneigh_node_sigma_base_idx + 0] = float4Zero;
        fneigh[fneigh_node_sigma_base_idx + 1] = float4Zero;
        fneigh[fneigh_node_sigma_base_idx + 2] = float4Zero;
        fneigh[fneigh_node_sigma_base_idx + 3] = float4Zero;

        // Clear pi recoil aforce for this node's 4 neighbors
        fneigh[fneigh_node_pi_base_idx + 0] = float4Zero;
        fneigh[fneigh_node_pi_base_idx + 1] = float4Zero;
        fneigh[fneigh_node_pi_base_idx + 2] = float4Zero;
        fneigh[fneigh_node_pi_base_idx + 3] = float4Zero;
    }
}


// ======================================================================
//                          getMMFFf4()
// ======================================================================

// Map OpenCL (get_global_id(0), get_global_id(1)) to CUDA (blockIdx.x*blockDim.x + threadIdx.x, blockIdx.y*blockDim.y + threadIdx.y)
// Assuming 2D launch grid (x=atoms, y=systems)

__global__ void getMMFFf4(
    int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    float4*  fapos,        // 3  [natoms]     aforce on    atoms (just node atoms are evaluated)
    float4*  fneigh,       // 4  [nnode*4*2]  recoil aforce on neighbors (and pi-orbitals)
    // parameters
    int4*    neighs,       // 5  [nnode]  neighboring atoms
    int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge
    float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    cu_Mat3* lvecs,        // 12 lattice vectors         for each system
    cu_Mat3* ilvecs,       // 13 inverse lattice vectors for each system // Unused in original kernel
    float4*  pbc_shifts,   // 14 PBC shifts precalculated from lvecs and cell indices
    const int npbc,        // 15 number of unique PBC shifts (must match size of pbc_shifts per system)
    const int bSubtractVdW // 16 Flag to enable/disable subtracting bonded non-bonded interactions
){

    const int iG = blockIdx.x * blockDim.x + threadIdx.x; // index of atom
    const int iS = blockIdx.y * blockDim.y + threadIdx.y; // index of system

    const int nAtoms = nDOFs.x;  // number of atoms in the system
    const int nnode  = nDOFs.y;  // number of nodes in the system
    const int nvec   = nAtoms+nnode; // total vectors (atoms + pi-orbitals)

    if(iG>=nnode) return; // Only process node atoms for this kernel

    const int i0a   = iS*nAtoms;         // index of first atom      in the system (in atoms array)
    const int i0n   = iS*nnode;          // index of first node atom in the system (in apars, bLs, bKs, Ksp, Kpp arrays)
    const int i0v   = iS*nvec;           // index of first vector    in the system (in apos, avel, aforce arrays)

    const int iaa = iG + i0a;  // index of current atom (in atoms, neighs, neighCell, REQKs)
    const int ian = iG + i0n;  // index of current node atom (in apars, bLs, bKs, Ksp, Kpp)
    const int iav = iG + i0v;  // index of current vector (in apos, fapos)

    //if(iG==iGdbg){ printf("CUDA getMMFFf4(): natoms=%i, nnode=%i nvec=%i iS=%i npbc=%i\n", nAtoms, nnode, nvec, iS, npbc ); }
    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("CUDA getMMFFf4(): natoms=%i, nnode=%i nvec=%i, npbc=%i\n", nAtoms, nnode, nvec, npbc ); }
    // if((iG==iGdbg)&&(iS==iSdbg)){
    //     for(int ia=0; ia<nnode; ia++){
    //         int4   ng=neighs[iaa+ia];
    //         float4 pi=apos[iav+ia];
    //         float4 bk=bLs[ian];
    //         float4 bK=bKs[ian];
    //         float4 Ks=Ksp[ian];
    //         float4 Kp=Kpp[ian];
    //         float4 apar=apars[ian];
    //         //printf("OCL getMMFFf4(): ia %3i: pos=(%10.5f,%10.5f,%10.5f) ngs=(%3i,%3i,%3i,%3i) bLs=(%10.5f,%10.5f,%10.5f,%10.5f) bKs=(%10.5f,%10.5f,%10.5f,%10.5f) Ks=(%10.5f,%10.5f,%10.5f,%10.5f) Kp=(%10.5f,%10.5f,%10.5f,%10.5f) apar=(%10.5f,%10.5f,%10.5f,%10.5f)\n", 
    //         //   ia, pi.x, pi.y, pi.z, ng.x, ng.y, ng.z, ng.w, bk.x, bk.y, bk.z, bk.w, bK.x, bK.y, bK.z, bK.w, Ks.x, Ks.y, Ks.z, Ks.w, Kp.x, Kp.y, Kp.z, Kp.w, apar.x, apar.y, apar.z, apar.w );
    //         printf("CUDA getMMFFf4(): ia %3i: pos=(%10.5f,%10.5f,%10.5f) ngs=(%3i,%3i,%3i,%3i) bLs=(%10.5f,%10.5f,%10.5f,%10.5f) bKs=(%10.5f,%10.5f,%10.5f,%10.5f) apar=(%10.5f,%10.5f,%10.5f,%10.5f)\n", 
    //             ia, pi.x, pi.y, pi.z, ng.x, ng.y, ng.z, ng.w, bk.x, bk.y, bk.z, bk.w, bK.x, bK.y, bK.z, bK.w, apar.x, apar.y, apar.z, apar.w );
    //     } 
    //     for(int ipbc=0; ipbc<npbc; ipbc++){
    //         float4 shift = pbc_shifts[ipbc];
    //         printf("CUDA getMMFFf4(): pbc %i: shift=(%10.5f,%10.5f,%10.5f,%10.5f)\n", ipbc, shift.x, shift.y, shift.z, shift.w); 
    //     }  
    // }

    // ---- Dynamical
    float4  hs [NNEIGH];         // direction vectors of bonds (h) and inverse bond lengths (h.w)
    float3  fbs[NNEIGH];         // force on neighbor sigma    (fbs[i] is sigma recoil force on i-th neighbor)
    float3  fps[NNEIGH];         // force on neighbor pi       (fps[i] is pi    recoil force on i-th neighbor)
    float3  fa  = float3Zero;    // force on center atom positon

    float E=0;                   // Total Energy of this atom
    // ---- Params
    const int4   ng    = neighs[iaa];    // neighboring atoms (indices relative to system start 0..natoms-1)
    const float3 pa    = XYZ(apos[iav]);  // position of current atom
    const float4  par  = apars[ian];     // (xy=s0_ss,z=ssK,w=piC0 ) forcefield parameters for current atom
    const int*   ings  = (const int*  )&ng; // neighboring atoms, we cast it to int[] to be index it in for loop

    for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }   // clear recoil aforce on neighbors

    float3 f1,f2;         // working aforce

    { // ========= BONDS - here we evaluate pairwise interactions of node atoms with its 4 neighbors

        float3  fpi = float3Zero;                // force on pi-orbital
        const int4   ngC = neighCell[iaa];       // neighboring atom cell index (index into pbc_shifts array *per system*)
        const float3 hpi = XYZ(apos[iav+nAtoms]); // direction of pi-orbital (index iav+nAtoms is the pi-orbital associated with this node atom)
        const float4 vbL = bLs[ian];             // bond lengths
        const float4 vbK = bKs[ian];             // bond stiffness
        const float4 vKs = Ksp[ian];             // stiffness of sigma-pi othogonalization
        const float4 vKp = Kpp[ian];             // stiffness of pi-pi    alignment

        const int*   ingC  = (const int*  )&ngC;   // neighboring atom cell index
        const float* bL    = (const float*)&vbL;   // bond lengths
        const float* bK    = (const float*)&vbK;   // bond stiffness
        const float* Kspi  = (const float*)&vKs;   // stiffness of sigma-pi othogonalization
        const float* Kppi  = (const float*)&vKp;   // stiffness of pi-pi    alignment

        const int ipbc0 = iS*npbc;  // index of first PBC shift for current system

        for(int i=0; i<NNEIGH; i++){  // loop over 4 neighbors
            float4 h;                 // direction vector of bond
            const int ing  = ings[i]; // index of i-th neighbor atom (relative to system start 0..natoms-1)
            if(ing<0) break; // if neighbor index is negative, this neighbor slot is empty

            const int ingv = ing+i0v; // index of i-th neighbor vector (in apos, avel, aforce)
            // Note: Original code used inga = ing+i0a, but REQKs, neighs, etc. are indexed by iaa, not iav. It should be fine.
            //const int inga = ing+i0a; // index of i-th neighbor atom (in atoms, neighs, neighCell, REQKs) - unused in this loop body

            // --- Compute bond direction vector and inverse bond length
            float3 h_xyz = XYZ(apos[ingv]) - pa;  // direction vector of bond
            { // shift bond to the proper PBC cell
                int ic  = ingC[i];                  // index of i-th neighbor cell (0..npbc-1)
                if(ic >= 0 && ic < npbc) { // Check boundary for safety
                    h_xyz  += XYZ(pbc_shifts[ipbc0+ic]); // shift bond to the proper PBC cell
                }
            }
            float  l = length(h_xyz);  // compute bond length
            h.w      = 1.f/l;           // store inverse bond length (use 1.f)
            h_xyz   *= h.w;            // normalize bond direction vector
            h.x = h_xyz.x; h.y = h_xyz.y; h.z = h_xyz.z; // Store back into float4
            hs[i]    = h;              // store bond direction vector (xyz) and inverse bond length (w)

            // --- Evaluate bond-length stretching energy and aforce
            if(iG<ing){
                float elb = evalBond(h_xyz, l-bL[i], bK[i], &f1);
                //if((iG==iGdbg)&&(iS==iSdbg)){ printf("CUDA getMMFFf4(): bond-length:    iG %3i ing %3i elb %10.5f f(%10.5f,%10.5f,%10.5f) l %10.5f bL %10.5f bK %10.5f \n", iG, ing, elb, f1.x, f1.y, f1.z, l, bL[i], bK[i] ); }
                fbs[i] -= f1;
                fa     += f1;
                E      +=elb;

                // pi-pi alignment interaction - only evaluated if both atoms are nodes
                float kpp = Kppi[i];
                if( (ing<nnode) && (kpp>1.e-6f) ){   // Only node atoms have pi-pi alignemnt interaction (ing < nnode)
                    float epp = evalPiAling( hpi, XYZ(apos[ingv+nAtoms]), kpp,  &f1, &f2 ); // apos[ingv+nAtoms] is the pi-orbital of neighbor ing
                    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("CUDA getMMFFf4(): cos(pi,pi):     iG %3i ing %3i epp %10.5f f(%10.5f,%10.5f,%10.5f) c %10.5f kpp %10.5f \n", iG, ing, epp, f1.x, f1.y, f1.z, dot(hpi,XYZ(apos[ingv+nAtoms])), kpp ); }
                    E     +=epp;
                    fpi   +=f1;
                    fps[i]+=f2; // Uses defined operator+=
                }
            }

            // pi-sigma orthogonalization interaction
            float ksp = Kspi[i];
            if(ksp>1.e-6f){
                float esp = evalAngCos( make_float4(hpi.x, hpi.y, hpi.z, 1.f), h, ksp, par.w, &f1, &f2 ); // hpi is direction, hr1.w=1.f; h is bond vector, hr2.w=1./l;
                //if((iG==iGdbg)&&(iS==iSdbg)){ printf("CUDA getMMFFf4(): cos(pi,sigma):      iG %3i ing %3i esp %10.5f f1(%10.5f,%10.5f,%10.5f) c %10.5f ksp %10.5f par.w %10.5f \n", iG, ing, esp, f1.x, f1.y, f1.z, dot(hpi,XYZ(h)), ksp, par.w ); }
                E     +=esp;
                fpi   +=f1;
                fa    -=f2; // Uses defined operator-=
                fbs[i]+=f2; // Uses defined operator+=
            }
        }

        // --- Store Pi-aforce
        // fneigh layout: [system][sigma/pi][node][neighbor]
        // Index for pi force on neighbor i of node iG in system iS: iS * nnode*8 + nnode*4 + iG*4 + i
        // Base index for pi aforce of node iG in system iS: iS * nnode*8 + nnode*4 + iG*4
        const int pi_fneigh_base_idx = iS * nnode * 8 + nnode * 4 + iG * 4;

        for(int i=0; i<NNEIGH; i++){
            if(ings[i]<0) break; // Only store if neighbor exists
            fneigh[pi_fneigh_base_idx + i] = make_float4(fps[i].x, fps[i].y, fps[i].z, 0.f); // store recoil pi-force on i-th neighbor
        }
        fapos[iav+nAtoms]  = make_float4(fpi.x, fpi.y, fpi.z, 0.f);  // store pi-force on pi-orbital of current atom (iav+nAtoms points to pi-orbital)

    }

    { //  ============== Angles   - here we evaluate angular interactions between pair of sigma-bonds of node atoms with its 4 neighbors

        
        //const float   ssC0 = par.x; // cos(ang0/2)
        //const float   ssS0 = par.y; // sin(ang0/2)
        //const float   ssK  = par.z; // sigma-sigma stiffness
        //const float   piC0 = par.w; // pi-sigma cos(ang0)

        for(int i=0; i<NNEIGH; i++){ // loop over first bond
            int ing = ings[i];
            if(ing<0) break; // if there is no i-th neighbor we break the loop
            const float4 hi = hs[i];
            //const int ingv = ing+i0v; // unused
            const int inga = ing+i0a; // index for REQKs

            for(int j=i+1; j<NNEIGH; j++){ // loop over second bond
                int jng  = ings[j];
                if(jng<0) break; // if there is no j-th neighbor we break the loop
                //const int jngv = jng+i0v; // unused
                const int jnga = jng+i0a; // index for REQKs
                const float4 hj = hs[j];

                float ea = evalAngleCosHalf( hi, hj, make_float2(par.x, par.y), par.z, &f1, &f2 ); // evaluate angular force and energy using cos(angle/2) formulation // NOTE: evalAngleCosHalf must be defined
                fa    -= (f1+f2); // Uses defined operator-= // total angular force on center atom is -(f1+f2)
                E     +=ea;

                //if((iG==iGdbg)&&(iS==iSdbg)){ printf("CUDA getMMFFf4(): angle():        iG %3i ing %3i jng %3i ea=%10.5f f1(%10.5f,%10.5f,%10.5f) f2(%10.5f,%10.5f,%10.5f) cos %10.5f apar(%10.5f,%10.5f,%10.5f) \n", iG, ing, jng, ea, f1.x, f1.y, f1.z, f2.x, f2.y, f2.z, dot(XYZ(hi),XYZ(hj)), par.x, par.y, par.z ); }

                if(bSubtractVdW){ // Remove non-bonded interactions from atoms that are bonded to common neighbor
                    float4 REQi=REQKs[inga];   // non-bonding parameters of i-th neighbor
                    float4 REQj=REQKs[jnga];   // non-bonding parameters of j-th neighbor
                    // combine non-bonding parameters of i-th and j-th neighbors using mixing rules
                    float4 REQij = float4Zero; // Initialize
                    REQij.x  = REQi.x  + REQj.x;
                    REQij.y = REQi.y * REQj.y;
                    REQij.z = REQi.z * REQj.z;
                    float3 hi_xyz = XYZ(hi);
                    float3 hj_xyz = XYZ(hj);
                    float3 dp = hi_xyz/hi.w - hj_xyz/hj.w; // vector from neighbor j to neighbor i
                    float4 fij = getLJQH( dp, REQij, 1.0f ); // fij is force on i due to j // NOTE: getLJQH must be defined
                    float3 fij_xyz = XYZ(fij);
                    f1 -=  fij_xyz; // subtract force on i // Uses defined operator-=
                    f2 +=  fij_xyz; // add force on j (f_on_j = -f_on_i = -fij) // Uses defined operator+=
                }

                fbs[i]+= f1; // Add modified angular force on neighbor i // Uses defined operator+=
                fbs[j]+= f2; // Add modified angular force on neighbor j // Uses defined operator+=
            }
        }
    }

    // ========= Save results - store aforce on atoms and recoil on its neighbors
    const int sigma_fneigh_base_idx = iS * nnode * 8 + iG * 8;

    for(int i=0; i<NNEIGH; i++){
        if(ings[i]<0) break; // Only store if neighbor exists
        fneigh[sigma_fneigh_base_idx + i] = make_float4(fbs[i].x, fbs[i].y, fbs[i].z, 0.f);
    }

    atomicAdd(&fapos[iav].x, fa.x);
    atomicAdd(&fapos[iav].y, fa.y);
    atomicAdd(&fapos[iav].z, fa.z);

}

// ======================================================================
//                           getNonBond()
// ======================================================================
// Calculate non-bonded aforce on atoms (including both node atoms and capping atoms), considering periodic boundary conditions
// Assuming launched over grid(natoms, nS)
__global__ void getNonBond(
    int4 ns,                  // 1 // (natoms,nnode,?,nAtomCeil_for_local_buffer) dimensions of the system and local buffer size info
    // Dynamical
    float4*  apos,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms]) - Note: pi-orbitals are NOT included in non-bonded.
    float4*  aforce,       // 3 // aforce on atoms (accumulated)
    // Parameters
    float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond?) per atom
    int4*    neighs,       // 5 // neighbors indices (0..natoms-1) per atom
    int4*    neighCell,    // 6 // neighbors cell indices (0..npbc-1) per atom
    cu_Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z) - Note: not used in loop bounds in original
    const float4      GFFParams     // 9 // Grid-Force-Field parameters (R2damp, Rcut, etc.)
){
    __shared__ float4 LATOMS[32];  // Shared memory for positions
    __shared__ float4 LCLJS [32];  // Shared memory for parameters
    const int iG = blockIdx.x * blockDim.x + threadIdx.x; // index of atom
    const int iS = blockIdx.y * blockDim.y + threadIdx.y; // index of system
    // total number of threads along y direction
    //const int nS = gridDim.y * blockDim.y;
    const int natoms=ns.x;               // number of atoms PER SYSTEM
    const int nnode =ns.y;               // number of node atoms PER SYSTEM
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)

    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int i0v = iS*nvec;    // index of first atom in vectors array
    //const int ian = iG + i0n;
    const int iaa = iG + i0a; // index of atom in atoms array
    const int iav = iG + i0v; // index of atom in vectors array

    const bool   bPBC   = (nPBC.x+nPBC.y+nPBC.z)>0;
    const float  R2damp = GFFParams.x*GFFParams.x; // R2damp = GFFParams.x^2

    float4 fe = float4Zero; // force on atom (accumulated locally)

    // Get data for the current atom if it's valid
    int4   ng    = { -1, -1, -1, -1 }; // Initialize neighbors to -1
    int4   ngC   = { -1, -1, -1, -1 }; // Initialize neighbor cells to -1
    float4 REQKi = float4Zero;
    float3 posi  = float3Zero;
    float4 posi4 = float4Zero; // Store original float4 for shared memory
    if (iG < natoms) {
        ng    = neighs   [iaa];
        ngC   = neighCell[iaa];
        REQKi = REQKs    [iaa];
        posi  = XYZ(apos [iaa]); // Use XYZ macro
    }

    const cu_Mat3 lvec = lvecs[iS]; // lattice vectors for this system

    //if((iG==iGdbg)&&(iS==iSdbg)){ printf("CUDA getNonBond(): natoms=%i, nnode=%i nSys=(%i,%i)  nPBC=(%i,%i,%i)\n", natoms, nnode, blockDim.y, gridDim.y, nPBC.x, nPBC.y, nPBC.z); }
    // if((iG==iGdbg)&&(iS==iSdbg)){ 
    //     printf("CUDA getNonBond(): lvec.a=(%g,%g,%g) lvec.b=(%g,%g,%g) lvec.c=(%g,%g,%g)\n", lvec.a.x, lvec.a.y, lvec.a.z, lvec.b.x, lvec.b.y, lvec.b.z, lvec.c.x, lvec.c.y, lvec.c.z);
    //     printf("CUDA getNonBond(): GFFParams=(%g,%g,%g,%g) \n", GFFParams.x, GFFParams.y, GFFParams.z, GFFParams.w);
    //     for(int i=0; i<natoms; i++){
    //         float4 pi = apos[i];
    //         int4   ng = neighs[i];
    //         int4   ngC = neighCell[i];
    //         float4 REQKi = REQKs[i];
    //         printf("CUDA getNonBond(): atom %i: ng=(%i,%i,%i,%i), ngC=(%i,%i,%i,%i), REQKi=(%10.5f,%10.5f,%10.5f|%10.5f), posi=(%10.5f,%10.5f,%10.5f,%10.5f)\n", i, ng.x, ng.y, ng.z, ng.w, ngC.x, ngC.y, ngC.z, ngC.w, REQKi.x, REQKi.y, REQKi.z, REQKi.w, pi.x, pi.y, pi.z, pi.w);
    //     }   
    // }
    // ========= Atom-to-Atom interaction in chunks
    for (int j0=0; j0<natoms; j0+=blockDim.x){ 
        int i = j0 + threadIdx.x + iS*natoms; 
        if (j0 + threadIdx.x < natoms) { 
            LATOMS[threadIdx.x] = apos  [i+i0v];
            LCLJS [threadIdx.x] = REQKs [i+i0a];
        }
        __syncthreads(); // Wait for all threads in the block to load data

        // Compute aforce between atom i (this thread) and atoms j in the shared memory chunk
        if (iG < natoms) { // Only compute for valid atom i
            for (int jl=0; jl<blockDim.x; jl++){    // loop over all atoms in local memory (like 32 atoms)
                const int ja = j0+jl; // index of atom j within the current system (0..natoms-1)
                if (ja >= natoms) continue; // Ensure atom j is valid

                const float4 aj4 = LATOMS[jl];    // read atom position from local memory
                const float3 aj  = XYZ(aj4);       // Use XYZ macro
                float4 REQK      = LCLJS [jl];    // read atom parameters from local memory

                // Mix parameters
                REQK.x  += REQKi.x;   // mixing rules for vdW Radius (Rij = Ri + Rj)
                REQK.y *= REQKi.y; // E0_ij = E0_i * E0_j
                REQK.z *= REQKi.z; // Q_ij = Q_i * Q_j

                const bool bSameAtom = (ja == iG);
                // Check if atom j is a direct neighbor of atom i
                const bool bBonded = ((ja == ng.x)||(ja == ng.y)||(ja == ng.z)||(ja == ng.w));


                if(!bSameAtom){ // No interaction with itself

                    float3 dp = aj - posi; // vector from atom i to atom j in the origin cell

                    if(bPBC){ // If PBC is enabled
                        //int ipbc = 0; // image index 0..8 (for 3x3 grid)
                        for(int iy=0; iy<3; iy++){ // Loop over Y images (-1, 0, 1) -> Map iy to -1,0,1
                            int shift_iy = iy - 1;
                            for(int ix=0; ix<3; ix++){ // Loop over X images (-1, 0, 1) -> Map ix to -1,0,1
                                int shift_ix = ix - 1;

                                float3 lvec_a_xyz = XYZ(lvec.a);
                                float3 lvec_b_xyz = XYZ(lvec.b);
                                float3 shift_vec = shift_ix * lvec_a_xyz + shift_iy * lvec_b_xyz; // No Z shift based on original loop
                                int ipbc_ = (shift_iy+1)*3 + (shift_ix+1);
                                if( !( bBonded && (                    // if atoms are bonded, we do not want to calculate non-covalent interaction between them
                                            ((ja==ng.x)&&(ipbc_==ngC.x)) || // check if this image corresponds to a bonded neighbor's cell
                                            ((ja==ng.y)&&(ipbc_==ngC.y)) ||
                                            ((ja==ng.z)&&(ipbc_==ngC.z)) ||
                                            ((ja==ng.w)&&(ipbc_==ngC.w)) ))
                                ){
                                    float3 dp_ = dp + shift_vec;
                                    float4 fij = getLJQH( dp, REQK, R2damp ); // NOTE: getLJQH must be defined
                                    fe += fij; // Uses defined operator+=
                                }
                            }
                        }
                    } else { // If PBC is not used
                        if( !bBonded ){ // if atoms are not bonded, calculate interaction (only in origin image)
                            float4 fij = getLJQH( dp, REQK, R2damp ); // NOTE: getLJQH must be defined
                            //if((iG==iGdbg)&&(iS==iSdbg)){   printf("CUDA getNonBond(): ia,ja %3i %3i aj(%10.5f,%10.5f,%10.5f) dp( %10.5f | %10.5f,%10.5f,%10.5f)  fij( %10.5f,%10.5f,%10.5f|%10.5f)\n", iG, ja, aj.x,aj.y,aj.z, length(dp), dp.x, dp.y, dp.z, fij.x, fij.y, fij.z, fij.w); }
                            fe += fij;
                        }
                    }
                }
            }
        }
        __syncthreads(); // Wait for all threads in the block to finish processing the chunk
    }

    
    if(iG<natoms){
        //if(iS==0){ printf( "GPU::getNonBond(iG=%i) fe(%g,%g,%g,%g)\n", iG, fe.x,fe.y,fe.z,fe.w ); }
        aforce[iav] = fe;           // If we do    run it as first forcefield, we can just store force (non need to clean it before in that case)
        //aforce[iav] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //aforce[iav] = fe*(-1.f);
    }
}


// ======================================================================
//                     updateAtomsMMFFf4()
// ======================================================================

__device__ unsigned int hash_wang(unsigned int bits) {
    bits = (bits ^ 61) ^ (bits >> 16);
    bits *= 9;
    bits = bits ^ (bits >> 4);
    bits *= 0x27d4eb2d;
    bits = bits ^ (bits >> 15);
    return bits;
}

__device__ float hashf_wang( float val, float xmin, float xmax) {
    unsigned int bits = __float_as_int(val);
    // Multiply by 2^-31 (approx 4.6566129e-10) to get a float in [0, 1)
    float u = (float)hash_wang( bits ) * 2.3283064e-10f; // Use 2^-32 instead of 2^-31
    return u *(xmax-xmin)+ xmin;
}


// Assemble recoil aforce from neighbors and  update atoms positions and velocities
// Assuming 2D launch grid (x=vectors, y=systems)
__global__ void updateAtomsMMFFf4(
    int4        n,            // 1 // (natoms,nnode,nsys,nMaxSysNeighs) dimensions and counts
    float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    float4*  avel,         // 3 // velocities of atoms
    float4*  aforce,       // 4 // aforce on atoms (cleared before getMMFFf4 and getNonBond)
    float4*  cvf,          // 5 // accumulated |f|^2, |v|^2, <f|v> per atom/pi-orbital
    float4*  fneigh,       // 6 // recoil aforce on neighbors (and pi-orbitals)
    int4*    bkNeighs,     // 7 // back neighbors indices (for recoil aforce) - Global indices into fneigh array
    float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit,seed_inc)
    float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    cu_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz) per system
    int*     sysneighs,    // 13 // for each system contains array int[nMaxSysNeighs] of nearby other systems
    float4*  sysbonds      // 14 // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
){
    const int natoms=n.x;           // number of atoms
    const int nnode =n.y;           // number of node atoms
    //const int nsys  =n.z;         // number of systems - unused directly, handled by gridDim.y
    const int nMaxSysNeighs = n.w;  // max number of inter-system interactions; if <0 switch inter system interactions off
    const int nvec  = natoms+nnode; // number of vectors (atoms+node atoms)
    const int iG = blockIdx.x * blockDim.x + threadIdx.x; // index of vector (atom or pi-orbital)

    if(iG>=nvec) return; // make sure we are not out of bounds of current system's vectors

    const int iS = blockIdx.y * blockDim.y + threadIdx.y; // index of system

    const int iaa = iG + iS*natoms;  // index of atom in constr, constrK arrays (only first natoms elements per system are relevant)
    const int iav = iG + iS*nvec;    // index of current vector (atom or pi) in apos, avel, aforce, cvf arrays

    const float4 MDpars  = MDparams[iS]; // (dt,damp,Flimit,seed_inc)
    const float4 TDrive = TDrives[iS]; // (T,gamma_damp,seed_base,?)

    float4 fe      = aforce[iav]; // force on atom or pi-orbital (includes aforce calculated by previous kernels)
    const bool bPi = iG>=natoms;  // is it pi-orbital ?

    // ------ Gather Forces from back-neighbors
    // bkNeighs contains GLOBAL indices into fneigh array
    // Example: bkNeighs[iav].x gives the index in fneigh where the recoil force from neighbor x of atom iav is stored.
    // This index must be correct regardless of which system atom iav belongs to.
    int4 ngs = bkNeighs[ iav ];

    float4 f_neigh_sum = float4Zero;
    // sum all recoil aforce from back neighbors
    if(ngs.x>=0){ f_neigh_sum += fneigh[ngs.x]; } // Uses defined operator+=
    if(ngs.y>=0){ f_neigh_sum += fneigh[ngs.y]; } // Uses defined operator+=
    if(ngs.z>=0){ f_neigh_sum += fneigh[ngs.z]; } // Uses defined operator+=
    if(ngs.w>=0){ f_neigh_sum += fneigh[ngs.w]; } // Uses defined operator+=
    fe += f_neigh_sum; // Add accumulated neighbor aforce

    // ---- Limit Forces (optional, can cause drift)
    float Flimit = MDpars.z; // Flimit is stored in MDparams.z
    float3 fe_xyz = XYZ(fe); // Use XYZ macro
    float fr2 = dot(fe_xyz, fe_xyz);
    if( fr2 > (Flimit*Flimit) ){
        fe_xyz *= (Flimit/sqrtf(fr2)); // Use sqrtf and defined operator*=
        fe.x = fe_xyz.x; fe.y = fe_xyz.y; fe.z = fe_xyz.z; // Store back
    }

    // =============== FORCE DONE

    // Note: The original OpenCL stores the limited force back into aforce, then reads it again for dynamics.
    // This is redundant if aforce was only used for accumulation. Let's just use the local fe for dynamics.
    // aforce[iav] = fe; // Optional: store limited force if needed elsewhere

    // =============== DYNAMICS

    float4 ve = avel[iav]; // velocity of atom or pi-orbital
    float4 pe = apos[iav]; // position of atom or pi-orbital

    float3 pe_xyz = XYZ(pe); // Use XYZ macro
    // fe_xyz is already defined and potentially limited

    // -------- Fixed Atoms and Bounding Box
    if(!bPi){ // only atoms (iG < natoms) have constraints and bboxes, not pi-orbitals
        // ------- bboxes
        const cu_Mat3 B = bboxes[iS];
        // Only apply if stiffness is positive (OpenCL had this check for z, let's apply to all)
        if(B.c.x>0.0f){ if(pe.x<B.a.x){ fe_xyz.x+=(B.a.x-pe.x)*B.c.x; }else if(pe.x>B.b.x){ fe_xyz.x+=(B.b.x-pe.x)*B.c.x; }; }
        if(B.c.y>0.0f){ if(pe.y<B.a.y){ fe_xyz.y+=(B.a.y-pe.y)*B.c.y; }else if(pe.y>B.b.y){ fe_xyz.y+=(B.b.y-pe.y)*B.c.y; }; }
        if(B.c.z>0.0f){ if(pe.z<B.a.z){ fe_xyz.z+=(B.a.z-pe.z)*B.c.z; }else if(pe.z>B.b.z){ fe_xyz.z+=(B.b.z-pe.z)*B.c.z; }; }

        // ------- constrains
        float4 cons = constr[ iaa ]; // constraints (x,y,z,K)
        // Use constr.w as a master switch, and constrK.xyz for stiffness per dimension
        if( cons.w > 0.f ){
            float4 cK = constrK[ iaa ];
            // Apply fmaxf component-wise
            cK.x = fmaxf( cK.x, 0.0f );
            cK.y = fmaxf( cK.y, 0.0f );
            cK.z = fmaxf( cK.z, 0.0f );
            // cK.w = fmaxf( cK.w, 0.0f ); // If needed
            const float3 fc = (XYZ(cons) - pe_xyz) * XYZ(cK); // Element-wise multiply
            fe_xyz += fc; // add constraint force // Uses defined operator+=
        }
    }

    // -------- Inter system interactions
    // This applies to ALL vectors (atoms and pi-orbitals) if nMaxSysNeighs > 0
    if( nMaxSysNeighs > 0 && iG < natoms){ // only atoms (not pi-orbitals) participate in inter-system interactions
        for(int i=0; i<nMaxSysNeighs; i++){
            const int j_bond_idx = iS*nMaxSysNeighs + i; // index in sysneighs and sysbonds
            const int    jS      = sysneighs[j_bond_idx]; // index of neighbor system
            if (jS < 0) continue; // Safety check

            const float4 bj = sysbonds [j_bond_idx]; // {Lmin,Lmax,Kpres,Ktens}
            const float4 pj = apos[jS*nvec + iG]; // position of corresponding atom in neighbor system
            float3 d        = XYZ(pj) - pe_xyz; // vector from current atom to neighbor system's atom

            float  l = length( d );
            float3 inter_f = float3Zero;
            if      (l<bj.x && bj.z > 1.e-6f){ // Lmin, Kpress
                inter_f = d * ( (l-bj.x)*bj.z / l ); // f = d/l * (l-Lmin)*Kpress // Uses defined operators
            }else if(l>bj.y && bj.w > 1.e-6f){ // Lmax, Ktens
                inter_f = d * ( (bj.y-l)*bj.w / l ); // f = d/l * (Lmax-l)*Ktens // Uses defined operators
            }
            fe_xyz += inter_f; // add inter-system force // Uses defined operator+=
        }
    }

    // Accumulate |f|^2, |v|^2, <f|v> for FIRE or other algorithms
    // Use atomic operations as multiple threads *could* theoretically write to the same cvf entry
    // if the grid/block setup was different, but for 1 thread per vector, atomics are not strictly needed,
    // but they don't hurt and make it safe if launch config changes.
    // Let's just use direct assignment as it's 1:1 thread:vector.
    float3 ve_xyz = XYZ(ve);
    cvf[iav].x = dot(fe_xyz, fe_xyz); // |f|^2
    cvf[iav].y = dot(ve_xyz, ve_xyz); // |v|^2
    cvf[iav].z = dot(fe_xyz, ve_xyz); // <f|v>
    cvf[iav].w = 0.0f; // Not used? Clear it.

    const bool bDrive = TDrive.y > 1.e-6f; // Check if damping is positive

    // ------ Move (Leap-Frog)
    float dt = MDpars.x; // dt is in MDparams.x
    // MDpars.z is Flimit, MDpars.y is damp in OpenCL example? Let's assume damp is MDpars.y
    float damp = MDpars.y; // damp is in MDpars.y (friction coeff)

    if(bPi){ // if pi-orbital, we need to make sure that it has unit length
        // Subtract force component that changes length (force along the vector pe)
        fe_xyz -= pe_xyz * dot( pe_xyz, fe_xyz ); // Uses defined operators
        // Subtract velocity component that changes length (velocity along the vector pe)
        ve_xyz -= pe_xyz * dot( pe_xyz, ve_xyz ); // Uses defined operators
    }else{ // For atoms (iG < natoms)
        // Thermal driving - Langevin thermostat
        if( bDrive ){
            // Damping force: -gamma * v
            fe_xyz    -= ve_xyz * TDrive.y; // TDrive.y is gamma_damp // Uses defined operators
            unsigned int seed_int = __float_as_int(TDrive.w); // Seed from TDrive.w
            float r1 = hashf_wang(__int_as_float(hash_wang(iG*136 + seed_int)), -1.0f, 1.0f);
            float r2 = hashf_wang(__int_as_float(hash_wang(iG*778 + seed_int)), -1.0f, 1.0f);
            float r3 = hashf_wang(__int_as_float(hash_wang(iG*4578 + seed_int)), -1.0f, 1.0f);
            float3 rand_vec = make_float3(r1, r2, r3);

            float random_mag = sqrtf( 2.0f * const_kB * TDrive.x * TDrive.y / dt ); // TDrive.x is T, TDrive.y is gamma, dt is MDpars.x
            fe_xyz += rand_vec * random_mag; // Uses defined operators
        }
    }

    ve_xyz *= damp; // MDparams.y is damp // Uses defined operators
    ve_xyz += fe_xyz * dt; // acceleration * dt // Uses defined operators
    pe_xyz += ve_xyz * dt; // velocity * dt // Uses defined operators

    if(bPi){ // if pi-orbital, maintain unit length after update
        pe_xyz = normalize(pe_xyz); // Use defined normalize
    }

    // Store updated xyz parts back into float4
    pe.x = pe_xyz.x; pe.y = pe_xyz.y; pe.z = pe_xyz.z;
    ve.x = ve_xyz.x; ve.y = ve_xyz.y; ve.z = ve_xyz.z;

    // Clear w components (often used for energy/mass/charge, cleared here for safety)
    pe.w = 0.f; // Use 0.f
    ve.w = 0.f; // Use 0.f

    // Store updated velocity and position
    avel[iav] = ve;
    apos[iav] = pe;
}


// ======================================================================
//                     printOnGPU()
// ======================================================================

__global__ void printOnGPU(
    int4        n,            // 1 // (natoms,nnode,isys_to_print,?)
    int4        mask,         // 2 // (print_atom_aforce, print_pi_aforce, print_fneigh, ?)
    float4*  apos,         // 3
    float4*  avel,         // 4
    float4*  aforce,       // 5
    float4*  fneigh,       // 6
    int4*    bkNeighs,     // 7
    float4*  constr        // 8
){
    // Assuming launched with grid(1, nS), block(1,1) or similar, and n.z is the target system index
    // If launched with a larger grid, only the first thread per system should print to avoid spam.
    const int thread_id_in_block = threadIdx.x + threadIdx.y * blockDim.x;
    if (thread_id_in_block != 0) return; // Only first thread of block prints

    const int isys  = n.z; // System index to print
    const int natoms= n.x;
    const int nnode = n.y;
    //const int nS = gridDim.y; // Total number of systems (if launched one block per system)
    const int nS_ = n.w; // Assuming n.w is passed as total systems for fneigh indexing

    const int i0a = isys * natoms; // base index for atoms data (constr, REQKs, neighs, etc.)
    const int i0v = isys * (natoms+nnode); // base index for vector data (apos, avel, aforce, cvf)

    printf( "CUDA::printOnGPU(isys=%i) natoms=%i nnode=%i nS_=%i\n", isys,  natoms, nnode, nS_ );
    if(mask.x){ // Print atom aforce/positions/constraints
        printf("--- Atoms (%i -- %i) ---\n", i0v, i0v+natoms-1);
        for(int i=0; i<natoms; i++){
            int ia = i + i0a; // index for constr (indexed by natoms)
            int iv = i + i0v; // index for vector data (indexed by nvec)
            printf( "CUDA[%i|isys=%i] i_vec=%i, i_atom=%i :: ", i, isys, iv, ia );
            //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[iv].x, bkNeighs[iv].y, bkNeighs[iv].z, bkNeighs[iv].w );
            printf( "aforce{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[iv].x, aforce[iv].y, aforce[iv].z, aforce[iv].w );
            //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[iv].x, avel[iv].y, avel[iv].z, avel[iv].w );
            printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[iv].x, apos[iv].y, apos[iv].z, apos[iv].w );
            if (i < natoms) { // constr is only for atoms
               printf(  "constr{%6.3f,%6.3f,%6.3f,%6.3f} ", constr[ia].x, constr[ia].y, constr[ia].z, constr[ia].w );
            }
            printf( "\n" );
        }
    }
    if(mask.y){ // Print pi aforce/positions
        printf("--- Pi Orbitals (%i -- %i) ---\n", i0v+natoms, i0v+(natoms+nnode)-1);
        for(int i=0; i<nnode; i++){ // Pi orbitals are associated with node atoms
            int ipi = i + natoms + i0v; // index for pi-orbital vector data
            printf( "CUDA[%i|isys=%i] i_vec=%i :: ", i, isys, ipi );
            printf(  "aforce_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[ipi].x, aforce[ipi].y, aforce[ipi].z, aforce[ipi].w );
            //printf(  "avel_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[ipi].x, avel[ipi].y, avel[ipi].z, avel[ipi].w );
            printf(   "apos_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[ipi].x, apos[ipi].y, apos[ipi].z, apos[ipi].w );
            printf( "\n" );
        }
    }
    if(mask.z){ // Print fneigh (recoil aforce)
        printf("--- Recoil Forces (fneigh) for nodes (%i) ---\n", nnode);
        // Based on getMMFFf4 layout: fneigh[ iS * nnode*8 + iG*8 + sigma/pi*4 + neigh_idx ]
        for(int i=0; i<nnode; i++){ // loop over nodes
            // Get base index for this node's recoil aforce in this system
            const int ingf = isys * nnode * 8 + i * 8;
            for(int j=0; j<4; j++){ // loop over neighbors
                int isigma = ingf + j;
                int ipi    = ingf + 4 + j;
                printf( "CUDA[node%i,neigh%i|isys=%i] :: ", i, j, isys );
                printf( "fneigh_sigma{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[isigma].x, fneigh[isigma].y, fneigh[isigma].z, fneigh[isigma].w );
                printf( "fneigh_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[ipi].x, fneigh[ipi].y, fneigh[ipi].z, fneigh[ipi].w );
                printf( "\n" );
            }
        }
    }
}



#endif // MMFF_KERNELS_CUH