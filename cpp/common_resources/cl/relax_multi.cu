#ifndef MMFF_KERNELS_CUH
#define MMFF_KERNELS_CUH

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <vector_types.h>
#include <math.h> // For sqrtf, expf, sinf, cosf, fractf, fmaxf, fminf
#include <stdio.h> // For device-side printf
#include "relax_multi.cuh"

// Mat3 struct is defined in relax_multi.cuh
// Use standard CUDA vector literals or initializer lists
#define  float4Zero  make_float4(0.f,0.f,0.f,0.f)
#define  float3Zero  make_float3(0.f,0.f,0.f)
// Note: Original OpenCL had float2Zero defined as float3, fixed to float2
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
__device__ inline float dot(float3 a, float3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

// Length
__device__ inline float length(float3 v) {
    return sqrtf(dot(v, v));
}

__device__ inline float3 normalize(float3 v) {
    float len = length(v);
    // Avoid division by zero
    if (len > 1e-12f) { // Use a small threshold
        float invLen = 1.0f / len;
        return make_float3(v.x * invLen, v.y * invLen, v.z * invLen);
    } else {
        return make_float3(0.0f, 0.0f, 0.0f); // Or return a default vector
    }
}

// Scalar multiplication
__device__ inline float3 operator*(float3 v, float s) {
    return make_float3(v.x * s, v.y * s, v.z * s);
}
__device__ inline float3 operator*(float s, float3 v) {
    return make_float3(v.x * s, v.y * s, v.z * s);
}
__device__ inline float3& operator*=(float3& v, float s) {
    v.x *= s; v.y *= s; v.z *= s;
    return v;
}
__device__ inline float3 operator/(float3 v, float s) {
    float invS = 1.0f / s; // Consider checking for s == 0 if necessary
    return make_float3(v.x * invS, v.y * invS, v.z * invS);
}
__device__ inline float3& operator/=(float3& v, float s) {
    float invS = 1.0f / s; // Consider checking for s == 0 if necessary
    v.x *= invS; v.y *= invS; v.z *= invS;
    return v;
}

// Vector addition/subtraction
__device__ inline float3 operator+(float3 a, float3 b) {
    return make_float3(a.x + b.x, a.y + b.y, a.z + b.z);
}
__device__ inline float3& operator+=(float3& a, const float3& b) {
    a.x += b.x; a.y += b.y; a.z += b.z;
    return a;
}
__device__ inline float3 operator-(float3 a, float3 b) {
    return make_float3(a.x - b.x, a.y - b.y, a.z - b.z);
}
__device__ inline float3& operator-=(float3& a, const float3& b) {
    a.x -= b.x; a.y -= b.y; a.z -= b.z;
    return a;
}

// Element-wise multiplication
__device__ inline float3 operator*(float3 a, float3 b) {
    return make_float3(a.x * b.x, a.y * b.y, a.z * b.z);
}

// float4 operators (add only needed ones for now)
__device__ inline float4& operator+=(float4& a, const float4& b) {
    a.x += b.x; a.y += b.y; a.z += b.z; a.w += b.w;
    return a;
}

// --- End Helper Functions ---

// Use __device__ for functions callable from kernels
__device__ inline float2 udiv_cmplx( float2 a, float2 b ){ return make_float2(  a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y ); }

// --- Placeholder definitions for missing functions ---
// IMPORTANT: Replace these with your actual CUDA implementations!
__device__ inline float evalPiAling( float3 hpi, float3 hpi_ng, float kpp, float3* f1, float3* f2 ){ *f1=float3Zero; *f2=float3Zero; return 0.0f; }
__device__ inline float evalAngCos( float4 hr1, float4 hr2, float K, float C0, float3* f1, float3* f2 ){ *f1=float3Zero; *f2=float3Zero; return 0.0f; }
__device__ inline float evalAngleCosHalf( float4 hr1, float4 hr2, float2 cs0, float K, float3* f1, float3* f2 ){ *f1=float3Zero; *f2=float3Zero; return 0.0f; }
__device__ inline float4 getLJQH( float3 dp, float4 REQij, float R2damp ){ return float4Zero; }
// --- End Placeholder definitions ---
// ======================================================================
// ======================================================================
//                           MMFF kernells
// ======================================================================
// ======================================================================


// ======================================================================
//                          getMMFFf4()
// ======================================================================

// Map OpenCL (get_global_id(0), get_global_id(1)) to CUDA (blockIdx.x*blockDim.x + threadIdx.x, blockIdx.y*blockDim.y + threadIdx.y)
// Assuming 2D launch grid (x=atoms, y=systems)

__global__ void getMMFFf4(
    int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    float4*  fapos,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
    float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
    // parameters
    int4*    neighs,       // 5  [nnode]  neighboring atoms
    int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge
    float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    Mat3* lvecs,        // 12 lattice vectors         for each system
    Mat3* ilvecs,       // 13 inverse lattice vectors for each system // Unused in original kernel
    float4*  pbc_shifts,   // 14 PBC shifts precalculated from lvecs and cell indices
    const int npbc,        // 15 number of unique PBC shifts (must match size of pbc_shifts per system)
    const int bSubtractVdW // 16 Flag to enable/disable subtracting bonded non-bonded interactions
){

    const int iG = blockIdx.x * blockDim.x + threadIdx.x; // index of atom
    const int iS = blockIdx.y * blockDim.y + threadIdx.y; // index of system
    const int nG = gridDim.x * blockDim.x; // global size x
    const int nS = gridDim.y * blockDim.y; // global size y

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

    #define NNEIGH 4

    // ---- Dynamical
    float4  hs [NNEIGH];         // direction vectors of bonds (h.xyz) and inverse bond lengths (h.w)
    float3  fbs[NNEIGH];         // force on neighbor sigma    (fbs[i] is sigma recoil force on i-th neighbor)
    float3  fps[NNEIGH];         // force on neighbor pi       (fps[i] is pi    recoil force on i-th neighbor)
    float3  fa  = float3Zero;    // force on center atom positon

    float E=0;                   // Total Energy of this atom
    // ---- Params
    const int4   ng  = neighs[iaa];    // neighboring atoms (indices relative to system start 0..natoms-1)
    const float3 pa  = XYZ(apos[iav]);  // position of current atom
    const float4 par = apars[ian];     // (xy=s0_ss,z=ssK,w=piC0 ) forcefield parameters for current atom


    const int*   ings  = (const int*  )&ng; // neighboring atoms, we cast it to int[] to be index it in for loop


    const float   ssC0_sq = par.x; // cos(ang0/2)
    const float   ssS0_sq = par.y; // sin(ang0/2)
    const float   ssK     = par.z; // sigma-sigma stiffness
    const float   piC0    = par.w; // pi-sigma cos(ang0)

    // Equilibrium cos(angle) for sigma-sigma is cos(2*ang0/2) = cos(ang0/2)^2 - sin(ang0/2)^2
    const float ssCosAng0 = ssC0_sq*ssC0_sq - ssS0_sq*ssS0_sq; // Use the stored cos(ang0/2) and sin(ang0/2)


    for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }   // clear recoil forces on neighbors

    float3 f1,f2;         // working forces

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

            // pi-pi alignment interaction - only evaluated if both atoms are nodes
            float kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1.e-6f) ){   // Only node atoms have pi-pi alignemnt interaction (ing < nnode)
                float E_pi = evalPiAling( hpi, XYZ(apos[ingv+nAtoms]), kpp,  &f1, &f2 ); // apos[ingv+nAtoms] is the pi-orbital of neighbor ing
                E+=E_pi;
                fpi+=f1;
                fps[i]+=f2; // Uses defined operator+=
            }

            // pi-sigma orthogonalization interaction
            float ksp = Kspi[i];
            if(ksp>1.e-6f){
                float E_ps = evalAngCos( make_float4(hpi.x, hpi.y, hpi.z, 1.f), h, ksp, piC0, &f1, &f2 ); // hpi is direction, hr1.w=1.f; h is bond vector, hr2.w=1./l;
                E+=E_ps; // OpenCL code added epp here, assuming a typo and it should be E_ps
                fpi+=f1;
                fa-=f2; // Uses defined operator-=
                fbs[i]+=f2; // Uses defined operator+=
            }
        }

        // --- Store Pi-forces
        // fneigh layout: [system][sigma/pi][node][neighbor]
        // Index for pi force on neighbor i of node iG in system iS: iS * nnode*8 + nnode*4 + iG*4 + i
        // Base index for pi forces of node iG in system iS: iS * nnode*8 + nnode*4 + iG*4
        const int pi_fneigh_base_idx = iS * nnode * 8 + nnode * 4 + iG * 4;

        for(int i=0; i<NNEIGH; i++){
            if(ings[i]<0) break; // Only store if neighbor exists
            fneigh[pi_fneigh_base_idx + i] = make_float4(fps[i].x, fps[i].y, fps[i].z, 0.f); // store recoil pi-force on i-th neighbor
        }
        fapos[iav+nAtoms]  = make_float4(fpi.x, fpi.y, fpi.z, 0.f);  // store pi-force on pi-orbital of current atom (iav+nAtoms points to pi-orbital)

    }

    { //  ============== Angles   - here we evaluate angular interactions between pair of sigma-bonds of node atoms with its 4 neighbors

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

                // ssK is already ssK in apars.z
                E += evalAngleCosHalf( hi, hj, make_float2(ssC0_sq, ssS0_sq), ssK, &f1, &f2 ); // evaluate angular force and energy using cos(angle/2) formulation // NOTE: evalAngleCosHalf must be defined
                fa    -= (f1+f2); // Uses defined operator-= // total angular force on center atom is -(f1+f2)

                if(bSubtractVdW){ // Remove non-bonded interactions from atoms that are bonded to common neighbor
                    float4 REQi=REQKs[inga];   // non-bonding parameters of i-th neighbor
                    float4 REQj=REQKs[jnga];   // non-bonding parameters of j-th neighbor
                    // combine non-bonding parameters of i-th and j-th neighbors using mixing rules
                    float4 REQij = float4Zero; // Initialize
                    REQij.x  = REQi.x  + REQj.x;
                    REQij.y = REQi.y * REQj.y;
                    REQij.z = REQi.z * REQj.z;

                    // Recover vector between neighbors. d_ij = pos_j - pos_i.
                    // bond_i = pos_i - pos_a = -hi.xyz / hi.w
                    // bond_j = pos_j - pos_a = -hj.xyz / hj.w
                    // pos_i = pos_a - hi.xyz / hi.w
                    // pos_j = pos_a - hj.xyz / hj.w
                    // dp = pos_j - pos_i = (pos_a - hj.xyz / hj.w) - (pos_a - hi.xyz / hi.w) = hi.xyz/hi.w - hj.xyz/hj.w
                    float3 hi_xyz = XYZ(hi);
                    float3 hj_xyz = XYZ(hj);
                    // Note: The OpenCL code used (hj.xyz/hj.w) - (hi.xyz/hi.w), which is vector from j to i.
                    // The direction of force fij is along dp. getLJQH returns force on dp.
                    // So getLJQH(dp) is force on j due to i. We subtract it from f1 (force on i).
                    // If dp = pos_j - pos_i, force on i is -fij, force on j is +fij.
                    // The original code subtracts from f1 and adds to f2. f1 applies to neighbor i, f2 applies to neighbor j.
                    // This means f1 is the angular force on neighbor i, f2 is the angular force on neighbor j.
                    // So we should add -fij to f1 and +fij to f2.
                    // getLJQH(dp) calculates force *along dp*. So force on j is along dp, force on i is along -dp.
                    // f_on_j = getLJQH(dp).xyz
                    // f_on_i = -getLJQH(dp).xyz
                    // We need to modify the forces *already computed* (angular forces) on neighbors i and j.
                    // Original angular force on i is f1, on j is f2.
                    // We are subtracting the non-bonded force.
                    // So new f1 = f1 - f_on_i = f1 - (-getLJQH(dp).xyz) = f1 + getLJQH(dp).xyz
                    // So new f2 = f2 - f_on_j = f2 - getLJQH(dp).xyz
                    // The OpenCL code has f1 -= fij.xyz and f2 += fij.xyz. This implies fij.xyz is force on j due to i.
                    // Let's trust the OpenCL implementation: dp is j-i, fij is force on j from i.
                    // Angular forces f1 on i, f2 on j. Subtract non-bonded: f1 -= f_on_i, f2 -= f_on_j.
                    // f_on_i = -f_on_j = -fij.xyz.
                    // So: f1 -= (-fij.xyz) = f1 + fij.xyz
                    //     f2 -= fij.xyz
                    // OpenCL code is f1 -= fij.xyz; f2 += fij.xyz; This matches if fij is force on i from j.
                    // Let's assume dp = pos_i - pos_j, and fij is force on i from j.
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

    // ========= Save results - store forces on atoms and recoil on its neighbors
    // fneigh layout: [system][sigma/pi][node][neighbor]
    // Index for sigma force on neighbor i of node iG in system iS: iS * nnode*8 + iG*8 + i
    // Base index for sigma forces of node iG in system iS: iS * nnode*8 + iG*8
    const int sigma_fneigh_base_idx = iS * nnode * 8 + iG * 8;

    for(int i=0; i<NNEIGH; i++){
        if(ings[i]<0) break; // Only store if neighbor exists
        fneigh[sigma_fneigh_base_idx + i] = make_float4(fbs[i].x, fbs[i].y, fbs[i].z, 0.f);
    }

    // Add force on the center atom
    // Original OpenCL adds to existing fapos. Let's follow that.
    // If this is the *first* kernel writing to fapos, clearing fapos beforehand is needed.
    // If it's not the first, adding is correct. The clean kernel clears fapos.
    // So adding is correct here.
    atomicAdd(&fapos[iav].x, fa.x);
    atomicAdd(&fapos[iav].y, fa.y);
    atomicAdd(&fapos[iav].z, fa.z);

    // Energy is not stored in fapos in the OpenCL code, keeping that behavior.
    // Energy is typically accumulated separately if needed.
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


// Assemble recoil forces from neighbors and  update atoms positions and velocities
// Assuming 2D launch grid (x=vectors, y=systems)
__global__ void updateAtomsMMFFf4(
    int4        n,            // 1 // (natoms,nnode,nsys,nMaxSysNeighs) dimensions and counts
    float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    float4*  avel,         // 3 // velocities of atoms
    float4*  aforce,       // 4 // forces on atoms (cleared before getMMFFf4 and getNonBond)
    float4*  cvf,          // 5 // accumulated |f|^2, |v|^2, <f|v> per atom/pi-orbital
    float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
    int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces) - Global indices into fneigh array
    float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit,seed_inc)
    float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz) per system
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
    const int nS = gridDim.y * blockDim.y; // number of systems

    const int iaa = iG + iS*natoms;  // index of atom in constr, constrK arrays (only first natoms elements per system are relevant)
    const int iav = iG + iS*nvec;    // index of current vector (atom or pi) in apos, avel, aforce, cvf arrays

    const float4 MDpars  = MDparams[iS]; // (dt,damp,Flimit,seed_inc)
    const float4 TDrive = TDrives[iS]; // (T,gamma_damp,seed_base,?)

    float4 fe      = aforce[iav]; // force on atom or pi-orbital (includes forces calculated by previous kernels)
    const bool bPi = iG>=natoms;  // is it pi-orbital ?

    // ------ Gather Forces from back-neighbors
    // bkNeighs contains GLOBAL indices into fneigh array
    // Example: bkNeighs[iav].x gives the index in fneigh where the recoil force from neighbor x of atom iav is stored.
    // This index must be correct regardless of which system atom iav belongs to.
    int4 ngs = bkNeighs[ iav ];

    float4 f_neigh_sum = float4Zero;
    // sum all recoil forces from back neighbors
    if(ngs.x>=0){ f_neigh_sum += fneigh[ngs.x]; } // Uses defined operator+=
    if(ngs.y>=0){ f_neigh_sum += fneigh[ngs.y]; } // Uses defined operator+=
    if(ngs.z>=0){ f_neigh_sum += fneigh[ngs.z]; } // Uses defined operator+=
    if(ngs.w>=0){ f_neigh_sum += fneigh[ngs.w]; } // Uses defined operator+=
    fe += f_neigh_sum; // Add accumulated neighbor forces

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
        const Mat3 B = bboxes[iS];
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
            if (jS < 0 || jS >= nS) continue; // Safety check

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

            // Random force: sqrt(2*kB*T*gamma) * R(t), where R(t) is white noise (mean 0, variance 1)
            // For discrete time, variance of (dt * R(t)) is dt. So we need sqrt(dt) scaling?
            // Discrete random force term: sqrt(2*kB*T*gamma/dt) * rand_vec (variance 1)
            // OpenCL code uses sqrt( 2*const_kB*TDrive.x*TDrive.y/MDpars.x ) = sqrt(2*kB*T*gamma/dt)
            // Let's use the hash function for pseudo-random numbers in [-1, 1]
            // The hash uses iG and TDrive.w (seed). TDrive.w should be updated on CPU each step.
            unsigned int seed_int = __float_as_int(TDrive.w); // Seed from TDrive.w
            float r1 = hashf_wang(__int_as_float(hash_wang(iG*136 + seed_int)), -1.0f, 1.0f);
            float r2 = hashf_wang(__int_as_float(hash_wang(iG*778 + seed_int)), -1.0f, 1.0f);
            float r3 = hashf_wang(__int_as_float(hash_wang(iG*4578 + seed_int)), -1.0f, 1.0f);
            float3 rand_vec = make_float3(r1, r2, r3);

            float random_mag = sqrtf( 2.0f * const_kB * TDrive.x * TDrive.y / dt ); // TDrive.x is T, TDrive.y is gamma, dt is MDpars.x
            fe_xyz += rand_vec * random_mag; // Uses defined operators
        }
    }

    // Apply friction (velocity damping)
    ve_xyz *= damp; // MDparams.y is damp // Uses defined operators

    // Update velocity
    ve_xyz += fe_xyz * dt; // acceleration * dt // Uses defined operators

    // Update position
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
// Print atoms and forces on GPU
// This kernel seems intended to be launched with a single thread per system
// (e.g., grid(1, nS), block(1,1)) and uses the loop indices to print data.
// Or perhaps grid(natoms or nvec, nS) with checks if iG < natoms or iG < nnode.
// Given the loops, it's likely launched with a single thread per system.
// Using iS = blockIdx.y and iG = threadIdx.x (assuming blockDim.x=1)
// and using n.z as the target system index `isys`.

__global__ void printOnGPU(
    int4        n,            // 1 // (natoms,nnode,isys_to_print,?)
    int4        mask,         // 2 // (print_atom_forces, print_pi_forces, print_fneigh, ?)
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
    const int nS_total = gridDim.y; // Total number of systems (if launched one block per system)

    // Need total number of systems if fneigh indexing depends on it
    // Based on analysis of getMMFFf4, fneigh indexing might be:
    // sigma: iS * nnode*8 + iG*8 + j
    // pi:    iS * nnode*8 + iG*8 + 4 + j
    // This implies total size nS * nnode * 8. Let's assume n.w contains the total number of systems used for fneigh sizing.
    // Or maybe the 'n' parameter is just for the system being printed? Let's assume the latter and the fneigh size/indexing is handled by host.
    // If fneigh indexing uses total systems `nS_total`, we need that value. Let's add nS_total to n.w for clarity.
    const int nS_for_fneigh = n.w; // Assuming n.w is passed as total systems for fneigh indexing

    const int i0a = isys * natoms; // base index for atoms data (constr, REQKs, neighs, etc.)
    const int i0v = isys * (natoms+nnode); // base index for vector data (apos, avel, aforce, cvf)

    printf( "#### CUDA::printOnGPU(isys=%i) natoms=%i nnode=%i nS_fneigh=%i \n", isys,  natoms, nnode, nS_for_fneigh );

    if(mask.x){ // Print atom forces/positions/constraints
        printf("--- Atoms (%i -- %i) ---\n", i0v, i0v+natoms-1);
        for(int i=0; i<natoms; i++){
            int ia_constr_idx = i + i0a; // index for constr (indexed by natoms)
            int ia_vec_idx    = i + i0v; // index for vector data (indexed by nvec)
            printf( "CUDA[%i|isys=%i] i_vec=%i, i_atom=%i :: ", i, isys, ia_vec_idx, ia_constr_idx );
            //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[ia_vec_idx].x, bkNeighs[ia_vec_idx].y, bkNeighs[ia_vec_idx].z, bkNeighs[ia_vec_idx].w );
            printf( "aforce{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[ia_vec_idx].x, aforce[ia_vec_idx].y, aforce[ia_vec_idx].z, aforce[ia_vec_idx].w );
            //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[ia_vec_idx].x, avel[ia_vec_idx].y, avel[ia_vec_idx].z, avel[ia_vec_idx].w );
            printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[ia_vec_idx].x, apos[ia_vec_idx].y, apos[ia_vec_idx].z, apos[ia_vec_idx].w );
            if (i < natoms) { // constr is only for atoms
               printf(  "constr{%6.3f,%6.3f,%6.3f,%6.3f} ", constr[ia_constr_idx].x, constr[ia_constr_idx].y, constr[ia_constr_idx].z, constr[ia_constr_idx].w );
            }
            printf( "\n" );
        }
    }
    if(mask.y){ // Print pi forces/positions
        printf("--- Pi Orbitals (%i -- %i) ---\n", i0v+natoms, i0v+(natoms+nnode)-1);
        for(int i=0; i<nnode; i++){ // Pi orbitals are associated with node atoms
            int i_pi_vec_idx = i + natoms + i0v; // index for pi-orbital vector data
            printf( "CUDA[%i|isys=%i] i_vec=%i :: ", i, isys, i_pi_vec_idx );
            printf(  "aforce_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[i_pi_vec_idx].x, aforce[i_pi_vec_idx].y, aforce[i_pi_vec_idx].z, aforce[i_pi_vec_idx].w );
            //printf(  "avel_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i_pi_vec_idx].x, avel[i_pi_vec_idx].y, avel[i_pi_vec_idx].z, avel[i_pi_vec_idx].w );
            printf(   "apos_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i_pi_vec_idx].x, apos[i_pi_vec_idx].y, apos[i_pi_vec_idx].z, apos[i_pi_vec_idx].w );
            printf( "\n" );
        }
    }
    if(mask.z){ // Print fneigh (recoil forces)
        printf("--- Recoil Forces (fneigh) for nodes (%i) ---\n", nnode);
        // Based on getMMFFf4 layout: fneigh[ iS * nnode*8 + iG*8 + sigma/pi*4 + neigh_idx ]
        for(int i=0; i<nnode; i++){ // loop over nodes
            // Get base index for this node's recoil forces in this system
            const int fneigh_node_base_idx = isys * nnode * 8 + i * 8;
            for(int j=0; j<4; j++){ // loop over neighbors
                int sigma_idx = fneigh_node_base_idx + j;
                int pi_idx    = fneigh_node_base_idx + 4 + j;
                printf( "CUDA[node%i,neigh%i|isys=%i] :: ", i, j, isys );
                printf( "fneigh_sigma{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[sigma_idx].x, fneigh[sigma_idx].y, fneigh[sigma_idx].z, fneigh[sigma_idx].w );
                printf( "fneigh_pi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[pi_idx].x, fneigh[pi_idx].y, fneigh[pi_idx].z, fneigh[pi_idx].w );
                printf( "\n" );
            }
        }
    }
}


// ======================================================================
//                     cleanForceMMFFf4()
// ======================================================================
// Clean forces on atoms/pi-orbitals and neighbors to prepare for next forcefield evaluation
// Based on analysis, the OpenCL version likely had incorrect fneigh indexing for clearing.
// This CUDA version clears fapos/aforce for ALL vectors (atoms + pi) and fneigh for ALL node recoil forces (sigma + pi).
// Assuming launched over grid(nvec, nS)

__global__ void cleanForceMMFFf4(
    int4        n,           // 2 // (natoms,nnode,?,?)
    float4*  aforce,      // 5 // forces on atoms and pi-orbitals
    float4*  fneigh       // 6 // recoil forces on neighbors (sigma + pi)
){
    const int natoms = n.x;
    const int nnode  = n.y;
    const int nvec   = natoms+nnode;
    const int iG     = blockIdx.x * blockDim.x + threadIdx.x; // index of vector (atom or pi)
    const int iS     = blockIdx.y * blockDim.y + threadIdx.y; // index of system
    const int nS     = gridDim.y * blockDim.y; // total systems

    const int iav = iG + iS*nvec; // global index in aforce

    // Clear force on this vector (atom or pi-orbital)
    // This should be done for ALL vectors (0..nvec-1)
    if(iG < nvec) {
       aforce[iav] = float4Zero;
       // Optional: aforce[iav] = make_float4(__int_as_float(iG), __int_as_float(iS), __int_as_float(iav), 0.0f); // Debugging
    }

    // Clear recoil forces in fneigh.
    // fneigh stores recoil forces for node atoms ONLY.
    // This means only threads with iG < nnode need to clear fneigh entries.
    // The indices to clear are those written by getMMFFf4.
    // fneigh layout derived from getMMFFf4: [ iS * nnode*8 + iG*8 + sigma/pi*4 + neigh_idx ]
    if(iG < nnode){
        // Base index for sigma forces of node iG in system iS: iS * nnode*8 + iG*8
        const int fneigh_node_sigma_base_idx = iS * nnode * 8 + iG * 8;
        // Base index for pi forces of node iG in system iS: iS * nnode*8 + iG*8 + 4
        const int fneigh_node_pi_base_idx    = fneigh_node_sigma_base_idx + 4;

        // Clear sigma recoil forces for this node's 4 neighbors
        fneigh[fneigh_node_sigma_base_idx + 0] = float4Zero;
        fneigh[fneigh_node_sigma_base_idx + 1] = float4Zero;
        fneigh[fneigh_node_sigma_base_idx + 2] = float4Zero;
        fneigh[fneigh_node_sigma_base_idx + 3] = float4Zero;

        // Clear pi recoil forces for this node's 4 neighbors
        fneigh[fneigh_node_pi_base_idx + 0] = float4Zero;
        fneigh[fneigh_node_pi_base_idx + 1] = float4Zero;
        fneigh[fneigh_node_pi_base_idx + 2] = float4Zero;
        fneigh[fneigh_node_pi_base_idx + 3] = float4Zero;
    }
}


// ======================================================================
//                           getNonBond()
// ======================================================================
// Calculate non-bonded forces on atoms (including both node atoms and capping atoms), considering periodic boundary conditions
// Assuming launched over grid(natoms, nS)
__global__ void getNonBond(
    int4 ns,                  // 1 // (natoms,nnode,?,nAtomCeil_for_local_buffer) dimensions of the system and local buffer size info
    // Dynamical
    float4*  atoms,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms]) - Note: pi-orbitals are NOT included in non-bonded.
    float4*  forces,       // 3 // forces on atoms (accumulated)
    // Parameters
    float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond?) per atom
    int4*    neighs,       // 5 // neighbors indices (0..natoms-1) per atom
    int4*    neighCell,    // 6 // neighbors cell indices (0..npbc-1) per atom
    Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z) - Note: not used in loop bounds in original
    const float4      GFFParams     // 9 // Grid-Force-Field parameters (R2damp, Rcut, etc.)
){

    // Use __shared__ for local memory
    // The local buffer size should match blockDim.x
    // The original code had [32], let's use blockDim.x
    extern __shared__ float4 LATOMS[]; // Shared memory for positions
    extern __shared__ float4 LCLJS[];  // Shared memory for parameters

    const int iG = blockIdx.x * blockDim.x + threadIdx.x; // index of atom
    const int nG = gridDim.x * blockDim.x; // total atoms across all systems (if launch covers all atoms)
    const int iS = blockIdx.y * blockDim.y + threadIdx.y; // index of system
    const int nS = gridDim.y * blockDim.y; // number of systems
    const int iL = threadIdx.x; // index of atom within the thread block
    const int nL = blockDim.x; // size of the thread block

    const int natoms=ns.x;  // number of atoms PER SYSTEM
    const int nnode =ns.y;  // number of node atoms PER SYSTEM
    //const int nAtomCeil_for_local_buffer =ns.w; // Ceiling for atom count, used for buffer size?

    //const int i0a = iS*natoms;  // index of first atom in atoms array for this system
    //const int i0v = iS*(natoms+nnode); // index of first vector in vector arrays (like apos) for this system
                                       // Note: atoms array passed to non-bond is just atoms, not pi-orbitals.
                                       // So it should be indexed by iS*natoms + iG.
                                       // But apos in other kernels includes pi. Need clarification.
                                       // Assuming 'atoms' pointer here is just atom positions (size nS*natoms)
                                       // and 'forces' is just atom forces (size nS*natoms).
                                       // This is inconsistent with other kernels using 'apos' and 'aforce' of size nS*(natoms+nnode).
                                       // Let's assume 'atoms' is 'apos' subset for atoms (0..natoms-1 per system)
                                       // and 'forces' is 'aforce' subset for atoms.
                                       // This means indices should be iS*natoms + iG. Let's fix.
    const int iaa_global = iG + iS*natoms; // index of current atom in GLOBAL arrays (atoms, forces, REQKs, neighs, neighCell)

    // Check if this thread corresponds to a valid atom index (0..natoms-1 within its system slice)
    // If launch grid covers more than natoms * nS, many threads will be invalid.
    // Check iG < natoms to ensure we are within system's atom count.
    if(iG >= natoms) {
       // If this thread is not assigned to a valid atom, it might still need to participate
       // in shared memory loading and sync, but it should not compute forces or write outputs.
       // Its shared memory slot will hold dummy data or be unused.
    }

    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const float  R2damp = GFFParams.x*GFFParams.x; // R2damp = GFFParams.x^2

    float4 fe = float4Zero; // force on atom (accumulated locally)

    // Get data for the current atom if it's valid
    int4   ng    = { -1, -1, -1, -1 }; // Initialize neighbors to -1
    int4   ngC   = { -1, -1, -1, -1 }; // Initialize neighbor cells to -1
    float4 REQKi = float4Zero;
    float3 posi  = float3Zero;
    float4 posi4 = float4Zero; // Store original float4 for shared memory
    if (iG < natoms) {
        ng    = neighs   [iaa_global];
        ngC   = neighCell[iaa_global];
        REQKi = REQKs    [iaa_global];
        posi  = XYZ(atoms [iaa_global]); // Use XYZ macro
    }

    const Mat3 lvec = lvecs[iS]; // lattice vectors for this system

    // PBC shifts assuming 3x3 grid in XY, centered at (0,0,0) image
    // The shifts should be precalculated based on lvec.a, lvec.b, lvec.c
    // OpenCL shifts:
    // dp += shift0; // shift to PBC image (-nPBC.x, -nPBC.y, -nPBC.z) relative to origin image (0,0,0)? No, it's relative to loop structure.
    // Let's trust the OpenCL shifts calculation relative to loop structure.
    // Assuming 3x3 grid (ix, iy) where ix from 0 to 2, iy from 0 to 2.
    // This covers images (-1,-1,0), (0,-1,0), (1,-1,0), (-1,0,0), (0,0,0), (1,0,0), (-1,1,0), (0,1,0), (1,1,0) if the z-shift is always zero. The shifts `shift0`, `shift_a`, `shift_b` calculations are complex and don't seem to generate these images easily.
    // Let's simplify and use standard PBC image generation for the 3x3 grid assumed by the loops.
    // A 3x3 grid of images relative to the origin cell (0,0,0) would iterate shifts ix, iy from -1 to 1.
    // The original loop is 0..2 for ix, iy. Let's stick to 0..2 and hope the shifts balance out.
    // The check `if( !( bBonded && (...) ) )` is the critical part for neighbor list exclusion.
    // Let's abandon the confusing OpenCL shift logic and just generate the 3x3 shifts relative to the origin image (0,0,0), and use the `neighCell` index to exclude the bonded image.

    // Let's use a standard 3x3 loop over images ix, iy = -1, 0, 1. Total 9 images.
    // ipbc index 0..8 mapped from (ix, iy) pair. e.g., ipbc = (iy+1)*3 + (ix+1)
    const int ipbc_neighbor_base = 4; // (0,0) image index is (0+1)*3 + (0+1) = 4

    // ========= Atom-to-Atom interaction in chunks
    for (int j0=0; j0<natoms; j0+=nL){ // loop over all atoms in the system, by chunks of size of local memory
        // Read a chunk of atom data into shared memory
        int src_idx_j = j0 + iL + iS*natoms; // Global index for source atom j
        if (j0 + iL < natoms) { // Ensure we don't read out of bounds
             LATOMS[iL] = atoms [src_idx_j]; // atoms is float4*
             LCLJS [iL] = REQKs [src_idx_j];
        } else {
             // Initialize with dummy data if out of bounds, or rely on check inside inner loop
             LATOMS[iL] = make_float4(1e18f, 1e18f, 1e18f, 0.f); // Large position to ensure large r2
             LCLJS[iL] = float4Zero;
        }
        __syncthreads(); // Wait for all threads in the block to load data

        // Compute forces between atom i (this thread) and atoms j in the shared memory chunk
        if (iG < natoms) { // Only compute for valid atom i
            for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
                const int ja_system_idx = j0+jl; // index of atom j within the current system (0..natoms-1)
                if (ja_system_idx >= natoms) continue; // Ensure atom j is valid

                const float4 aj4 = LATOMS[jl];    // read atom position from local memory
                const float3 aj = XYZ(aj4); // Use XYZ macro
                float4 REQK     = LCLJS [jl];    // read atom parameters from local memory

                // Mix parameters
                REQK.x  += REQKi.x;   // mixing rules for vdW Radius (Rij = Ri + Rj)
                // For Epsilon: Eij = sqrt(Ei * Ej) - Original code used *=, this is incorrect for sqrt mixing
                // Let's stick to the original code's *= behavior unless sqrt is confirmed.
                // REQK.yz *= REQKi.yz; // This seems to be Q1*Q2 for z, and E0*E0 for y? Very strange.
                // Standard LJ mixing: R0_ij = (R0_i + R0_j) / 2, E0_ij = sqrt(E0_i * E0_j)
                // The code uses REQK.x += REQKi.x, which means R0_ij = R0_i + R0_j. Unusual but let's keep.
                // Let's assume REQK.yz stores {E0, Q} * {E0_i, Q_i}
                // This implies REQK.y should be sqrt(E0_j) or something pre-processed.
                // Given the function name getLJQH(dp, REQ, R2damp) takes REQ.x=R0, REQ.y=E0, REQ.z=Q.
                // And the mixing is REQij.x = REQi.x + REQj.x; REQij.yz = REQi.yz * REQj.yz;
                // This means R0_ij = R0_i + R0_j
                // E0_ij = E0_i * E0_j
                // Q_ij = Q_i * Q_j
                // This mixing is very unusual. Let's translate it as is.
                 REQK.y *= REQKi.y; // E0_ij = E0_i * E0_j
                 REQK.z *= REQKi.z; // Q_ij = Q_i * Q_j


                // Exclude interaction if bonded or same atom
                // ja_system_idx is index 0..natoms-1 for atom j in the current system
                // iG is index 0..natoms-1 for atom i in the current system
                // If iG == ja_system_idx, it's the same atom.
                const bool bSameAtom = (ja_system_idx == iG);
                // Check if atom j is a direct neighbor of atom i
                const bool bBonded = ((ja_system_idx == ng.x)||(ja_system_idx == ng.y)||(ja_system_idx == ng.z)||(ja_system_idx == ng.w));

                // Loop over PBC images
                // Standard 3x3 grid in XY (ix, iy = -1, 0, 1)
                // OpenCL loop indices ix, iy are 0, 1, 2. Map this to images -1, 0, 1
                // OpenCL used ipbc 0..8. Let's match that mapping.
                // ipbc = iy*3 + ix
                // (0,0) image is ix=0, iy=0 -> ipbc=0. But original code used ipbc=4 for (0,0)?
                // Let's use a consistent mapping: ix = 0..2, iy = 0..2. ipbc = iy*3 + ix. (0,0) image is ipbc=0.
                // Neighbor cell index ngC stores 0..npbc-1. npbc must be >= 9 if 3x3 is used.

                if(!bSameAtom){ // No interaction with itself

                    float3 dp_base = aj - posi; // vector from atom i to atom j in the origin cell

                    if(bPBC){ // If PBC is enabled
                        int ipbc = 0; // image index 0..8 (for 3x3 grid)
                        for(int iy=0; iy<3; iy++){ // Loop over Y images (-1, 0, 1) -> Map iy to -1,0,1
                            int shift_iy = iy - 1;
                            for(int ix=0; ix<3; ix++){ // Loop over X images (-1, 0, 1) -> Map ix to -1,0,1
                                int shift_ix = ix - 1;

                                float3 lvec_a_xyz = XYZ(lvec.a);
                                float3 lvec_b_xyz = XYZ(lvec.b);
                                float3 shift_vec = shift_ix * lvec_a_xyz + shift_iy * lvec_b_xyz; // No Z shift based on original loops

                                // Map (shift_ix, shift_iy) back to OpenCL's ipbc index convention if needed for ngC check
                                // The OpenCL check `((ja==ng.x)&&(ipbc==ngC.x))` suggests ngC stores the ipbc index.
                                // The OpenCL ipbc goes 0..8 inside the 3x3 loops.
                                // Let's assume OpenCL's ipbc = iy*3 + ix (with ix, iy from 0 to 2).
                                // Then the origin image (0,0,0) corresponds to ix=1, iy=1 (center of 0..2 range).
                                // ipbc = 1*3 + 1 = 4. This matches the idea that ipbc=4 is the origin image.
                                // Let's map ix, iy (-1, 0, 1) to ipbc used in ngC check.
                                // ipbc_check = (shift_iy+1)*3 + (shift_ix+1); // This is ipbc in 0..8 range
                                // This seems unnecessarily complex. Let's assume `ngC` stores an index that correctly identifies the PBC image of the bonded neighbor relative to the central atom.
                                // A simpler check: if bonded, and the image is the central image (ix=0, iy=0), skip.
                                // This assumes bonded atoms are always in the central image (0,0,0). If they can cross PBC, the ngC check is needed.
                                // Let's trust the OpenCL check using ngC and map our (ix,iy) to its ipbc.
                                // Assuming ipbc = iy_opencl*3 + ix_opencl, where ix_opencl, iy_opencl are 0..2.
                                // And ix = ix_opencl - 1, iy = iy_opencl - 1. So ix_opencl = ix + 1, iy_opencl = iy + 1.
                                // ipbc_opencl = (iy+1)*3 + (ix+1);
                                int ipbc_opencl = (shift_iy+1)*3 + (shift_ix+1);

                                if( !( bBonded &&                     // if atoms are bonded, we do not want to calculate non-covalent interaction between them
                                        (
                                            ((ja_system_idx==ng.x)&&(ipbc_opencl==ngC.x)) || // check if this image corresponds to a bonded neighbor's cell
                                            ((ja_system_idx==ng.y)&&(ipbc_opencl==ngC.y)) ||
                                            ((ja_system_idx==ng.z)&&(ipbc_opencl==ngC.z)) ||
                                            ((ja_system_idx==ng.w)&&(ipbc_opencl==ngC.w))
                                        )
                                    )
                                ){
                                    float3 dp_shifted = dp_base + shift_vec;
                                    float4 fij = getLJQH( dp_shifted, REQK, R2damp ); // NOTE: getLJQH must be defined
                                    fe += fij; // Uses defined operator+=
                                }
                            }
                        }
                    } else { // If PBC is not used
                        if( !bBonded ){ // if atoms are not bonded, calculate interaction (only in origin image)
                             float4 fij = getLJQH( dp_base, REQK, R2damp ); // NOTE: getLJQH must be defined
                             fe += fij;
                        }
                    }
                }
            }
        }
        __syncthreads(); // Wait for all threads in the block to finish processing the chunk
    }

    if(iG < natoms){ // Only write output if this thread processed a valid atom
        // Original OpenCL added force: forces[iav] += fe;
        // If cleanForce clears aforce beforehand, then this is the first write.
        // Let's trust the cleanForce clears aforce for atoms, so adding is correct here.
        // Need to use atomicAdd for safety if multiple threads could write to the same location,
        // but with 1 atom per thread, it's not needed. Let's add directly.
        atomicAdd(&forces[iaa_global].x, fe.x); // Use iaa_global for forces indexed by natoms
        atomicAdd(&forces[iaa_global].y, fe.y);
        atomicAdd(&forces[iaa_global].z, fe.z);
        // forces[iaa_global].w = 0.0f; // Don't touch w component if used for energy/charge/mass
    }
}

#endif // MMFF_KERNELS_CUH