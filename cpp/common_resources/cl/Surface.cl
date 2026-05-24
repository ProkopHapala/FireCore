// Surface.cl
// Kernels for surface electrostatics and folded-basis potential evaluation
// Copied from relax_multi.cl for use in pyOpenCL surface potential calculations
//
// ============================================================
//  DESIGN NOTES: OpenCL Ewald2D Implementation Plan
// ============================================================
//
// Goal: Efficient GPU evaluation of 2D Ewald electrostatic potential
// Reference: pyBall/Ewald2D.py (Python implementation)
//
// The Ewald2D potential is computed in two parts:
//
// Part 1: Charge Density Projection (compute C_G coefficients)
// -----------------------------------------------------------
// Input:
//   - Ion positions: rx, ry, rz (N_ions)
//   - Ion charges: q (N_ions)
//   - Lattice vectors: a_vec, b_vec (2D)
//   - Reciprocal vectors: b1, b2 (computed from a,b)
//   - G-vectors: Gx, Gy, Gn (N_G, with |h|,|k| <= n_harm, G≠0)
//
// Output:
//   - Complex coefficients C_G (N_G) for vacuum evaluation
//   - Complex per-ion weights w[g,i] (N_G, N_ions) for full evaluation
//
// Formula (complex form, matches fixed Python implementation):
//   C_G = (2π/(A|G|)) * Σ_i q_i * exp(|G| z_i) * exp(-i G·ρ_i)
//   w[g,i] = (2π/(A|G|)) * q_i * exp(-i G·ρ_i)
//
// Kernel design:
//   - Kernel: compute_ewald_coefficients
//   - Work items: 1D, one per G-vector (N_G work items)
//   - Each work item computes C_G by summing over all ions
//   - Use parallel reduction if N_ions is large (for now, simple loop)
//   - Store results as float2 (real, imag) in global memory
//
// Data structures:
//   - G_vectors: float4 array (Gx, Gy, Gn, 0) - N_G elements
//   - ion_data: float4 array (rx, ry, rz, q) - N_ions elements
//   - C_G_out: float2 array (real, imag) - N_G elements
//   - w_out: float2 array (N_G * N_ions) - flattened [g*N_ions + i]
//
// Part 2: Potential Evaluation
// -----------------------------
// Two modes:
//
// 2a) Vacuum evaluation (phi_vacuum_xy):
//    - Input: C_G coefficients, evaluation points (X, Y, z)
//    - Formula: φ = Re( Σ_G C_G * exp(i G·ρ) * exp(-|G| z) )
//    - Use for: XY slices at z > max(z_i)
//    - Kernel: eval_potential_vacuum
//    - Work items: 2D grid over evaluation points (Nx, Ny)
//    - Each work item sums over all G-vectors
//
// 2b) Full evaluation (phi_full_2d, phi_full_1d):
//    - Input: w[g,i] per-ion weights, evaluation points (X, Y, Z)
//    - Formula: φ = -(2π/A) Σ_i q_i |z - z_i| + Re( Σ_G Σ_i w[g,i] * exp(i G·ρ) * exp(-|G||z-z_i|) )
//    - Use for: XZ slices, 1D line scans, any z
//    - Kernel: eval_potential_full
//    - Work items: 2D grid over evaluation points (Nx, Nz for XZ slice)
//    - Each work item sums over G-vectors and ions
//
// Kernel design considerations:
//   - For small N_G (typical n_harm=3-5, N_G~48), simple loop is fine
//   - For large N_G, consider using local memory for G-vectors
//   - Use native_exp, native_cos, native_sin for performance
//   - COULOMB_CONST = 14.3996448915 (eV·Å/e²) must be applied at end
//
// Part 3: Brute Force Reference (for validation)
// ---------------------------------------------
// Kernel: eval_potential_brute
// - Direct Coulomb sum over periodic images
// - Formula: φ = Σ_{n,m} Σ_i q_i / |r - R_{nm} - r_i|
// - Use N_rep circular shells of PBC images
// - Slow but exact (within N_rep)
// - For validation only, not production
//
// ============================================================
//  Implementation Steps
// ============================================================
//
// Step 1: Compute reciprocal lattice and G-vectors (Python side)
//   - Use pyBall/Ewald2D.make_reciprocal_2d()
//   - Use pyBall/Ewald2D.generate_G_vectors()
//   - Upload to GPU as float4 array
//
// Step 2: Compute coefficients (GPU kernel)
//   - Kernel: compute_ewald_coefficients
//   - Input: ion_data (rx,ry,rz,q), G_vectors (Gx,Gy,Gn), area
//   - Output: C_G (float2), w (float2 flattened)
//   - Compare with Python: pyBall/Ewald2D.compute_C_G(), compute_w_per_ion()
//
// Step 3: Vacuum potential evaluation (GPU kernel)
//   - Kernel: eval_potential_vacuum
//   - Input: C_G, G_vectors, evaluation grid (X,Y,z)
//   - Output: potential (float)
//   - Compare with Python: pyBall/Ewald2D.phi_vacuum_xy()
//
// Step 4: Full potential evaluation (GPU kernel)
//   - Kernel: eval_potential_full
//   - Input: w, G_vectors, ion_data, evaluation grid (X,Y,Z)
//   - Output: potential (float)
//   - Compare with Python: pyBall/Ewald2D.phi_full_2d(), phi_full_1d()
//
// Step 5: Brute force reference (GPU kernel)
//   - Kernel: eval_potential_brute
//   - Input: ion_data, lattice vectors, evaluation points, N_rep
//   - Output: potential (float)
//   - Compare with Python: pyBall/Ewald2D.phi_brute_1d()
//
// Step 6: Python wrapper
//   - Create pyBall/OCL/SurfaceEwald.py
//   - Class SurfaceEwaldCL:
//     - __init__: load CL program, build kernels
//     - compute_coefficients(): call compute_ewald_coefficients
//     - eval_vacuum(): call eval_potential_vacuum
//     - eval_full(): call eval_potential_full
//     - eval_brute(): call eval_potential_brute
//
// Step 7: Validation test
//   - Test script: tests/tMMFF/test_ewald_cl.py
//   - Compare GPU results with Python reference (pyBall/Ewald2D)
//   - Test cases:
//     - NaCl surface (same as GridFF comparison)
//     - Random ion configurations
//     - Check RMSE < 1e-6 eV
//
// ============================================================
//  Notes on Data Layout and Memory Access
// ============================================================
//
// Ion data layout (float4):
//   - .x = rx
//   - .y = ry
//   - .z = rz
//   - .w = q
//
// G-vector layout (float4):
//   - .x = Gx
//   - .y = Gy
//   - .z = Gn (|G|)
//   - .w = unused
//
// Complex numbers stored as float2:
//   - .x = real part
//   - .y = imag part
//
// Evaluation grid:
//   - For XY slice: 2D grid, flatten to 1D for kernel
//   - For XZ slice: 2D grid, flatten to 1D for kernel
//   - Use work-item ID to index into flattened arrays
//
// ============================================================

// ============================================================
//  Constants and Helper Functions
// ============================================================

#define COULOMB_CONST 14.3996448915f  // eV·Å/e²

// Type definitions (from relax_multi.cl)
typedef struct cl_Mat3 {
    float4 a;
    float4 b;
    float4 c;
} cl_Mat3;

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}

inline float4 getMorsePLQH( float3 dp, float4 REQH, float4 PLQH, float K, float R2damp ){
    float r2    = dot(dp,dp);
    float ir2_  = 1/(r2+R2damp);
    float r     = sqrt( r2   );
    float ir_   = sqrt( ir2_ );
    float e     = exp ( K*(r-REQH.x));
    float Ee    = REQH.y*e;
    float EP    = Ee*e      * PLQH.x;
    float EL    = -2.0f*Ee  * PLQH.y;
    float EQ    = COULOMB_CONST*REQH.z*ir_ * PLQH.z;
    float frP   = (2.0f*K*EP)/r;
    float frL   = (-2.0f*K*Ee*PLQH.y)/r;
    float frQ   = -EQ*ir2_;
    return (float4){ dp*(frP+frL+frQ), EP+EL+EQ };
}

// evaluate damped Coulomb potential and force
inline float4 getCoulomb( float3 dp, float R2damp ){
    // ---- Electrostatic
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 + R2damp);
    float   E    = COULOMB_CONST*sqrt( ir2_ );
    return  (float4){ dp*-E*ir2_, E };
}

inline float macro_phi_rect_dipole(float3 p, float4 Pz, float4 AB) {
    float Ax = AB.x;
    float Bx = AB.y;
    float x = p.x;
    float y = p.y;
    float z = p.z;
    float sumOmega = 0.0f;
    float sumLogY  = 0.0f;
    float sumLogX  = 0.0f;
    float xs[2] = {-Ax, Ax};
    float ys[2] = {-Bx, Bx};
    for (int ix=0; ix<2; ix++) {
        for (int iy=0; iy<2; iy++) {
            float X = x - xs[ix];
            float Y = y - ys[iy];
            float R = sqrt(X*X + Y*Y + z*z);
            float s = ((ix==0)?-1.0f:1.0f) * ((iy==0)?-1.0f:1.0f);
            sumOmega += s * atan2( X*Y, z * R + 1e-12f );
            sumLogY  += s * log( Y + R + 1e-12f );
            sumLogX  += s * log( X + R + 1e-12f );
        }
    }
    return (Pz.z * sumOmega) - (Pz.x * sumLogY) - (Pz.y * sumLogX);
}

inline float rect_sheet_F(float X, float Y, float Z){
    float R = sqrt(X*X + Y*Y + Z*Z);
    return X*log(Y + R + 1e-12f) + Y*log(X + R + 1e-12f) - Z*atan2(X*Y, Z*R + 1e-12f);
}

inline float macro_phi_rect_charge(float3 p, float4 AB){
    float Ax = AB.x;
    float By = AB.y;
    float x0 = p.x + Ax;
    float x1 = p.x - Ax;
    float y0 = p.y + By;
    float y1 = p.y - By;
    return rect_sheet_F(x0,y0,p.z) - rect_sheet_F(x1,y0,p.z) - rect_sheet_F(x0,y1,p.z) + rect_sheet_F(x1,y1,p.z);
}

inline float4 getMacroRectLayers( float3 pos, float q, float4 bounds, float4 L0, float4 L1, float4 L2, float4 S0, float4 Q0, float4 Q1, float4 Q2, int nlayer ){
    float Ax = 0.5f*(bounds.y - bounds.x);
    float By = 0.5f*(bounds.w - bounds.z);
    float cx = 0.5f*(bounds.y + bounds.x);
    float cy = 0.5f*(bounds.w + bounds.z);
    float3 p = pos - (float3)(cx,cy,0.0f);
    float phi = 0.0f;
    float4 ls[3] = {L0,L1,L2};
    float sigmas[3] = {S0.x,S0.y,S0.z};
    float4 qs[3] = {Q0,Q1,Q2};
    for(int i=0; i<nlayer; i++){
        float4 Li = ls[i];
        float3 pp = (float3)(p.x,p.y,p.z-Li.w);
        float4 AB = (float4)(Ax,By,0.0f,0.0f);
        phi += sigmas[i] * macro_phi_rect_charge( pp, AB );
        // dipole contribution
        float4 Pz = (float4)(qs[i].x, qs[i].y, qs[i].z, 0.0f);
        phi += q * macro_phi_rect_dipole( pp, Pz, AB );
    }
    // potential gradient (force) - TODO: implement gradient
    return (float4){0.0f, 0.0f, 0.0f, phi};
}

inline float folded_eval_basis(float u, float v, float z, float4 prm){
    float bx = cos( (2.0f*M_PI_F) * prm.x * u );
    float by = cos( (2.0f*M_PI_F) * prm.y * v );
    float dz = fmax(0.0f, z - prm.w);
    float bz = exp( -prm.z * dz );
    return bx * by * bz;
}

inline float3 folded_eval_grad(float u, float v, float z, float4 prm, float4 invLvec2d){
    float phix = (2.0f*M_PI_F) * prm.x;
    float phiy = (2.0f*M_PI_F) * prm.y;
    float cu = cos(phix*u);
    float su = sin(phix*u);
    float cv = cos(phiy*v);
    float sv = sin(phiy*v);
    float dz = fmax(0.0f, z - prm.w);
    float bz = exp(-prm.z * dz);
    float dEdu = -phix * su * cv * bz;
    float dEdv = -phiy * cu * sv * bz;
    float dEdz = (z > prm.w) ? (-prm.z * cu * cv * bz) : 0.0f;
    float dudx = invLvec2d.x;
    float dudy = invLvec2d.z;
    float dvdx = invLvec2d.y;
    float dvdy = invLvec2d.w;
    return (float3)( dEdu*dudx + dEdv*dvdx, dEdu*dudy + dEdv*dvdy, dEdz );
}

// limit force magnitude to fmax
float3 limnitForce( float3 f, float fmax ){
    float fr2 = dot(f,f);                         // force magnitude squared
    if( fr2>(fmax*fmax) ){ f*=(fmax/sqrt(fr2)); } // if force magnitude is larger than fmax we scale it down to fmax
    return f;
}

float4 getR4repulsion( float3 d, float R, float Rcut, float A ){
    // we use R4blob(r) = A * (1-r^2)^2
    // such that at distance r=R we have force f = fmax
    // f = -dR4blob/dr = 4*A*r*(1-r^2) = fmax
    // A = fmax/(4*R*(1-R^2))
    float R2    = R*R;
    float R2cut = Rcut*Rcut;
    float r2 = dot(d,d);
    if( r2>R2cut ){
        return (float4){0.0f,0.0f,0.0f,0.0f};
    }else if( r2>R2 ){
        float mr2 = R2cut-r2;
        float fr = A*mr2;
        return (float4){ d*(-4*fr), fr*mr2 };
    }else{
        float mr2 = R2cut-R2;
        float fr = A*mr2;
        return (float4){ d*(-4*fr), fr*mr2 };
    }
}

inline int4 make_inds_pbc(const int n, const int iG) {
    // Generate PBC index patterns for B-spline interpolation
    // Returns 4 indices: (i0, i1, i2, i3) for 4-point B-spline
    // Handles wrapping at boundaries
    int4 inds;
    int i = iG % n;
    inds.x = (i - 1 + n) % n;
    inds.y = i;
    inds.z = (i + 1) % n;
    inds.w = (i + 2) % n;
    return inds;
}

// ============================================================
//  Brute Force Surface Interaction (getSurfMorse)
// ============================================================
// This is brute-force alternative to GridFF - describes interaction
// of molecule with substrate by pairwise interactions with multiple replicas

__kernel void getSurfMorse(
    const int4 ns,                // 1
    __global float4*  atoms,      // 2
    __global float4*  REQs,       // 3
    __global float4*  forces,     // 4
    __global float4*  atoms_s,    // 5
    __global float4*  REQ_s,      // 6
    __global float4*  surf_mpos,  // 7  (xmin,xmax,ymin,ymax)
    __global float4*  surf_mdip,  // 8  (mx,my,mz,0)
    __global float4*  surf_mQa,   // 9  Q row a
    __global float4*  surf_mQb,   // 10 Q row b
    __global float4*  surf_mQc,   // 11 (sigma0,sigma1,sigma2,Qtot)
    __global float4*  surf_qQa,   // 12 layer quadrupole (Qxx,Qxy,Qyy,z0)
    __global float4*  surf_qQb,   // 13 layer quadrupole (Qxx,Qxy,Qyy,z1)
    __global float4*  surf_qQc,   // 14 layer quadrupole (Qxx,Qxy,Qyy,z2)
    const int4     nPBC,          // 15
    const cl_Mat3  lvec,          // 16
    const float4   pos0,          // 17
    const float4   GFFParams,     // 18
    const float4   PLQH           // 19   (Pauli, London, Coulomb, HBond)
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

    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system

    float4 fe   = (float4){0.0f,0.0f,0.0f,0.0f};

    if(iG>=nAtoms) return;

    const float  K          = -GFFParams.y;
    const float  R2damp     =  GFFParams.x*GFFParams.x;
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz)
    const int bMacro      = (int)(GFFParams.z>0.5f);

    const float3 pos  = atoms[iav].xyz - pos0.xyz +  lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell
    const float4 REQi = REQs [iaa];

    for (int j0=0; j0<na_surf; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms_s[i];
        LCLJS [iL] = REQ_s  [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( ja<na_surf ){
                float4 REQH =       LCLJS [jl];
                float3 dp   = pos - LATOMS[jl].xyz;
                REQH.x   += REQi.x;
                REQH.yzw *= REQi.yzw;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                            float4 fej = getMorsePLQH( dp, REQH, PLQH, K, R2damp );
                            fe -= fej;
                            dp   +=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                    }
                    dp   +=shift_c;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if( bMacro && (fabs(PLQH.z) > 1e-12f) && (fabs(REQi.z) > 1e-12f) ){
        int nlayer = (int)(GFFParams.w + 0.5f);
        float4 fm = getMacroRectLayers( atoms[iav].xyz, REQi.z, surf_mpos[iS], surf_mdip[iS], surf_mQa[iS], surf_mQb[iS], surf_mQc[iS], surf_qQa[iS], surf_qQb[iS], surf_qQc[iS], nlayer );
        fe.xyz += fm.xyz;
        fe.w   += fm.w;
    }

    forces[iav] += fe;
}

// ============================================================
//  Folded Basis Evaluation (getSurfFolded)
// ============================================================

__kernel void getSurfFolded(
    const int4 ns,                     // 1
    __global float4*  atoms,           // 2
    __global float4*  REQs,            // 3
    __global float4*  forces,          // 4
    __global float*   folded_coeffs,   // 5  [ntypeMax*nbasisMax]
    __global float4*  folded_kxyz,     // 6  [nbasisMax]
    __global int*     folded_atom_type,// 7  [natoms]
    const int4        folded_meta,     // 8  (nbasis, ntypes, 0, 0)
    const float4      folded_lvec2d    // 9  (ax,bx,ay,by)
){
    __local float4 LBASIS[64];
    __local float  LCOEFFS[8*64];

    const int iG = get_global_id(0);
    const int iS = get_global_id(1);
    const int iL = get_local_id(0);
    const int nL = get_local_size(0);

    const int natoms = ns.x;
    const int nnode  = ns.y;
    const int nvec   = natoms + nnode;
    const int i0a    = iS*natoms;
    const int i0v    = iS*nvec;
    const int iaa    = iG + i0a;
    const int iav    = iG + i0v;
    if(iG>=natoms) return;

    const int nbasis = folded_meta.x;
    const int ntypes = folded_meta.y;
    if(nbasis<=0) return;
    if(nbasis>64){ return; }
    if(ntypes>8 ){ return; }

    for(int j=iL; j<nbasis; j+=nL){
        LBASIS[j] = folded_kxyz[j];
    }
    for(int j=iL; j<nbasis*ntypes; j+=nL){
        LCOEFFS[j] = folded_coeffs[j];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    float ax = folded_lvec2d.x;
    float bx = folded_lvec2d.y;
    float ay = folded_lvec2d.z;
    float by = folded_lvec2d.w;
    float det = ax*by - bx*ay;
    if(fabs(det) < 1e-12f) return;
    float4 invLvec2d = (float4)( by/det, -bx/det, -ay/det, ax/det );

    float3 pos = atoms[iav].xyz;
    float u = invLvec2d.x*pos.x + invLvec2d.y*pos.y;
    float v = invLvec2d.z*pos.x + invLvec2d.w*pos.y;
    u = u - floor(u);
    v = v - floor(v);
    int ityp = folded_atom_type[iG];
    if(ityp < 0 || ityp >= ntypes) return;

    float E = 0.0f;
    float3 F = (float3)(0.0f,0.0f,0.0f);
    int ioff = ityp*nbasis;
    for(int ib=0; ib<nbasis; ib++){
        float c = LCOEFFS[ioff + ib];
        float4 prm = LBASIS[ib];
        float  b = folded_eval_basis(u, v, pos.z, prm);
        float3 g = folded_eval_grad (u, v, pos.z, prm, invLvec2d);
        E += c * b;
        F -= c * g;
    }
    forces[iav] += (float4)(F.x, F.y, F.z, -E);
}

// ============================================================
//  Folded Basis Workgroup-Optimized (getSurfFolded_workgroup)
// ============================================================

#define MAX_ATOMS 64
#define MAX_XY 4
#define MAX_Z  8

__kernel void getSurfFolded_workgroup(
    const int4 ns,                     // (natoms, nnode, 0, 0)
    __global float4*  atoms,           
    __global float4*  REQs,            
    __global float4*  forces,          
    __global float*   folded_coeffs,   
    __global float4*  folded_kxyz,     // [Nxy params, Nz params]
    __global int*     folded_atom_type,
    const int4        folded_meta,     // (N_xy, N_z, ntypes, 0) 
    const float4      folded_lvec2d    
){
    const int iG = get_global_id(0);
    const int iS = get_global_id(1);
    const int iL = get_local_id(0);    // Thread ID (0 to 63) maps to Atom index within batch
    const int nL = get_local_size(0);  // 64

    const int natoms = ns.x;
    const int Nxy = folded_meta.x; 
    const int Nz  = folded_meta.y;
    const int ntypes = folded_meta.z;
    const int nbasis_total = Nxy * Nxy * Nz;

    // ==================================================================
    // 1. ALLOCATE __LOCAL MEMORY FOR EXPLICIT PRECALCULATION STORAGE
    // ==================================================================
    // Coefficients and parameters
    __local float  LCOEFFS[MAX_XY * MAX_XY * MAX_Z * 8]; 
    __local float4 LPARAMS_XY[MAX_XY]; 
    __local float4 LPARAMS_Z[MAX_Z];

    // Evaluated 1D Basis Arrays [Atom_Index][Basis_Index]
    __local float L_BX [MAX_ATOMS][MAX_XY];
    __local float L_dBX[MAX_ATOMS][MAX_XY];
    __local float L_BY [MAX_ATOMS][MAX_XY];
    __local float L_dBY[MAX_ATOMS][MAX_XY];
    __local float L_BZ [MAX_ATOMS][MAX_Z];
    __local float L_dBZ[MAX_ATOMS][MAX_Z];

    // Cooperative parameter loading
    for(int j = iL; j < Nxy; j += nL) LPARAMS_XY[j] = folded_kxyz[j];
    for(int j = iL; j < Nz;  j += nL) LPARAMS_Z[j]  = folded_kxyz[Nxy + j];
    for(int j = iL; j < nbasis_total * ntypes; j += nL) LCOEFFS[j] = folded_coeffs[j];

    barrier(CLK_LOCAL_MEM_FENCE);

    int active = (iG < natoms);
    int ityp = active ? folded_atom_type[iG] : -1;
    active = active && (ityp >= 0) && (ityp < ntypes);

    // Geometry transforms
    float det = folded_lvec2d.x * folded_lvec2d.w - folded_lvec2d.y * folded_lvec2d.z;
    float4 invLvec = (float4)(folded_lvec2d.w/det, -folded_lvec2d.y/det, -folded_lvec2d.z/det, folded_lvec2d.x/det);

    int iav = iG + iS * (natoms + ns.y);
    float3 pos = (float3)(0.0f, 0.0f, 0.0f);
    if(active){ pos = atoms[iav].xyz; }
    
    float u = invLvec.x * pos.x + invLvec.y * pos.y;
    float v = invLvec.z * pos.x + invLvec.w * pos.y;
    u -= floor(u);
    v -= floor(v);

    // ==================================================================
    // 2. PARALLEL PRECALCULATION -> SAVE TO LOCAL MEMORY
    // Every thread calculates its own atom's basis and explicitly saves 
    // it to its dedicated row in the Local Memory array.
    // ==================================================================
    for(int i = 0; i < Nxy; i++){
        float k = LPARAMS_XY[i].x; 
        float phi = 2.0f * M_PI_F * k;
        
        float phix_u = phi * u;
        L_BX[iL][i]  = active ? native_cos(phix_u) : 0.0f;
        L_dBX[iL][i] = active ? (-phi * native_sin(phix_u)) : 0.0f;
        
        float phiy_v = phi * v;
        L_BY[iL][i]  = active ? native_cos(phiy_v) : 0.0f;
        L_dBY[iL][i] = active ? (-phi * native_sin(phiy_v)) : 0.0f;
    }

    for(int i = 0; i < Nz; i++){
        float kz = LPARAMS_Z[i].z;
        float z0 = LPARAMS_Z[i].w;
        float dz = fmax(0.0f, pos.z - z0);
        float bz = active ? native_exp(-kz * dz) : 0.0f;
        L_BZ[iL][i]  = bz;
        L_dBZ[iL][i] = active && (pos.z > z0) ? (-kz * bz) : 0.0f;
    }

    barrier(CLK_LOCAL_MEM_FENCE);

    // ==================================================================
    // 3. THE TRIPLE LOOP
    // Thread streams its precalculated 1D factors from Local Memory,
    // avoiding the risk of register spilling entirely.
    // ==================================================================
    float E_tot = 0.0f;
    float dEdu_tot = 0.0f;
    float dEdv_tot = 0.0f;
    float dEdz_tot = 0.0f;

    int ic = active ? (ityp * nbasis_total) : 0; // Pointer to coefficients

    for(int iz = 0; iz < Nz; iz++){
        float bz  = L_BZ[iL][iz];
        float dbz = L_dBZ[iL][iz];

        for(int iy = 0; iy < Nxy; iy++){
            float by  = L_BY[iL][iy];
            float dby = L_dBY[iL][iy];
            
            // Outer loop multipliers
            float bz_by  = bz * by;
            float dbz_by = dbz * by;
            float bz_dby = bz * dby;

            for(int ix = 0; ix < Nxy; ix++){
                float bx  = L_BX[iL][ix];
                float dbx = L_dBX[iL][ix];

                float c = LCOEFFS[ic++]; 

                // Dynamic 3D Basis Construction
                E_tot    += c * (bx * bz_by);
                dEdu_tot += c * (dbx * bz_by);
                dEdv_tot += c * (bx * bz_dby);
                dEdz_tot += c * (bx * dbz_by);
            }
        }
    }

    // Map gradients back to forces
    float3 F_tot;
    F_tot.x = -(dEdu_tot * invLvec.x + dEdv_tot * invLvec.z);
    F_tot.y = -(dEdu_tot * invLvec.y + dEdv_tot * invLvec.w);
    F_tot.z = -dEdz_tot;

    if(active){ forces[iav] += (float4)(F_tot.x, F_tot.y, F_tot.z, -E_tot); }
}

// ============================================================
//  Folded Basis Harmonics (getSurfFolded_harmonics)
// ============================================================

__kernel void getSurfFolded_harmonics(
    const int4 ns,                     
    __global float4*  atoms,           
    __global float4*  REQs,            
    __global float4*  forces,          
    __global float*   folded_coeffs,   
    __global float4*  folded_kxyz,     // Now stores 1D params: [Nx params, Ny params, Nz params]
    __global int*     folded_atom_type,
    const int4        folded_meta,     // (Nx, Ny, Nz, ntypes)
    const float4      folded_lvec2d    
){    
    // Local memory for coefficients and 1D parameters
    __local float  LCOEFFS[MAX_XY * MAX_XY * MAX_Z * 8];
    __local float4 LBASIS[(2 * MAX_XY) + MAX_Z];

    const int iG = get_global_id(0);
    const int iS = get_global_id(1);
    const int iL = get_local_id(0);
    const int nL = get_local_size(0);
    const int natoms = ns.x;
    
    if(iG >= natoms) return;

    // Tensor product dimensions
    const int Nx = folded_meta.x;
    const int Ny = folded_meta.y;
    const int Nz = folded_meta.z;
    const int ntypes = folded_meta.w;
    const int nbasis_total = Nx * Ny * Nz;
    const int nparams_1d = Nx + Ny + Nz;

    // TODO: Complete harmonics kernel implementation
}

// ============================================================
//  OpenCL Ewald2D Kernels (GPU-accelerated surface electrostatics)
// ============================================================
// Reference: pyBall/Ewald2D.py (Python implementation)
// Key optimization: Use complex multiplication to compute e^{iG·ρ}
//
// For G = h*b1 + k*b2:
//   e^{iG·ρ} = e^{i(h*b1 + k*b2)·ρ} = e^{ih*b1·ρ} * e^{ik*b2·ρ}
//
// Precompute z1_b1 = e^{i*b1·ρ}, z1_b2 = e^{i*b2·ρ}
// Then:
//   e^{ih*b1·ρ} = z1_b1^h (by repeated multiplication)
//   e^{ik*b2·ρ} = z1_b2^k
//
// This reduces N_G cos/sin evaluations to just 2 per point!

// Helper: Complex multiply
inline float2 cmul(float2 a, float2 b) {
    return (float2)(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

// ------------------------------------------------------------------
// Kernel 1: Compute C_G coefficients (vacuum) and w[g,i] (full)
// ------------------------------------------------------------------
// Each work item computes coefficients for one G-vector
// Work size: N_G (number of G-vectors)
__kernel void compute_ewald_coefficients(
    __global const float4* ion_data,
    __global const float4* G_data,
    __global const float2* b_vectors,
    const float area,
    const int N_ions,
    const int N_G,
    __global float2* C_G_out,
    __global float2* w_out
){
    const int ig = get_global_id(0);
    if(ig >= N_G) return;

    float4 G = G_data[ig];
    int h = (int)G.x;
    int k = (int)G.y;
    float Gn = G.z;

    float2 b1 = b_vectors[0];
    float2 b2 = b_vectors[1];

    float Gx = h * b1.x + k * b2.x;
    float Gy = h * b1.y + k * b2.y;

    float prefactor = (2.0f * M_PI_F) / (area * Gn);

    float2 C_G = (float2)(0.0f, 0.0f);

    for(int i = 0; i < N_ions; i++){
        float4 ion = ion_data[i];
        float rx = ion.x;
        float ry = ion.y;
        float rz = ion.z;
        float q = ion.w;

        float Gdotr = Gx * rx + Gy * ry;
        float cos_gr = cos(Gdotr);
        float sin_gr = sin(Gdotr);
        float2 phase = (float2)(cos_gr, -sin_gr);

        float decay_ion = exp(Gn * rz);
        float2 contrib = (float2)(q * decay_ion * phase.x, q * decay_ion * phase.y);
        C_G += contrib;

        if(w_out != NULL){
            float2 w_gi = (float2)(q * phase.x * prefactor, q * phase.y * prefactor);
            w_out[ig * N_ions + i] = w_gi;
        }
    }

    C_G_out[ig] = (float2)(C_G.x * prefactor, C_G.y * prefactor);
}

// ------------------------------------------------------------------
// Kernel 2: Vacuum potential evaluation
// ------------------------------------------------------------------
__kernel void eval_potential_vacuum(
    __global const float4* eval_points,
    __global const float2* C_G,
    __global const float4* G_data,
    __global const float2* b_vectors,
    const int N_points,
    const int N_G,
    const int n_harm,
    __global float* phi_out
){
    const int ip = get_global_id(0);
    if(ip >= N_points) return;

    float4 p = eval_points[ip];
    float x = p.x;
    float y = p.y;
    float z = p.z;

    float2 b1 = b_vectors[0];
    float2 b2 = b_vectors[1];

    float b1dotr = b1.x * x + b1.y * y;
    float b2dotr = b2.x * x + b2.y * y;
    float2 z1_b1 = (float2)(cos(b1dotr), sin(b1dotr));
    float2 z1_b2 = (float2)(cos(b2dotr), sin(b2dotr));

    float phi = 0.0f;

    for(int ig = 0; ig < N_G; ig++){
        float4 G = G_data[ig];
        int h = (int)G.x;
        int k = (int)G.y;
        float Gn = G.z;

        float2 zh_b1 = (float2)(1.0f, 0.0f);
        int h_abs = abs(h);
        for(int i = 0; i < h_abs; i++){
            zh_b1 = cmul(zh_b1, z1_b1);
        }
        if(h < 0) zh_b1.y = -zh_b1.y;

        float2 zk_b2 = (float2)(1.0f, 0.0f);
        int k_abs = abs(k);
        for(int i = 0; i < k_abs; i++){
            zk_b2 = cmul(zk_b2, z1_b2);
        }
        if(k < 0) zk_b2.y = -zk_b2.y;

        float2 phase = cmul(zh_b1, zk_b2);
        float decay = exp(-Gn * z);
        float2 C = C_G[ig];
        float2 contrib = cmul(C, phase);

        phi += contrib.x * decay;
    }

    phi_out[ip] = phi * COULOMB_CONST;
}

// ------------------------------------------------------------------
// Kernel 3: Full potential evaluation (any z)
// ------------------------------------------------------------------
__kernel void eval_potential_full(
    __global const float4* eval_points,
    __global const float2* w,
    __global const float4* ion_data,
    __global const float4* G_data,
    __global const float2* b_vectors,
    const float area,
    const int N_points,
    const int N_ions,
    const int N_G,
    __global float* phi_out
){
    const int ip = get_global_id(0);
    if(ip >= N_points) return;

    float4 p = eval_points[ip];
    float x = p.x;
    float y = p.y;
    float z = p.z;

    float2 b1 = b_vectors[0];
    float2 b2 = b_vectors[1];

    float b1dotr = b1.x * x + b1.y * y;
    float b2dotr = b2.x * x + b2.y * y;
    float2 z1_b1 = (float2)(cos(b1dotr), sin(b1dotr));
    float2 z1_b2 = (float2)(cos(b2dotr), sin(b2dotr));

    float phi0 = 0.0f;
    for(int i = 0; i < N_ions; i++){
        float4 ion = ion_data[i];
        float q = ion.w;
        float rz = ion.z;
        phi0 -= q * fabs(z - rz);
    }
    phi0 *= (2.0f * M_PI_F / area);

    float phi_G = 0.0f;

    for(int ig = 0; ig < N_G; ig++){
        float4 G = G_data[ig];
        int h = (int)G.x;
        int k = (int)G.y;
        float Gn = G.z;

        float2 zh_b1 = (float2)(1.0f, 0.0f);
        int h_abs = abs(h);
        for(int i = 0; i < h_abs; i++){
            zh_b1 = cmul(zh_b1, z1_b1);
        }
        if(h < 0) zh_b1.y = -zh_b1.y;

        float2 zk_b2 = (float2)(1.0f, 0.0f);
        int k_abs = abs(k);
        for(int i = 0; i < k_abs; i++){
            zk_b2 = cmul(zk_b2, z1_b2);
        }
        if(k < 0) zk_b2.y = -zk_b2.y;

        float2 phase = cmul(zh_b1, zk_b2);

        for(int i = 0; i < N_ions; i++){
            float4 ion = ion_data[i];
            float rz = ion.z;
            float decay = exp(-Gn * fabs(z - rz));
            float2 w_gi = w[ig * N_ions + i];
            float2 contrib = cmul(w_gi, phase);
            phi_G += contrib.x * decay;
        }
    }

    phi_out[ip] = (phi0 + phi_G) * COULOMB_CONST;
}

// ------------------------------------------------------------------
// Kernel 4: Brute force Coulomb sum (reference/validation)
// ------------------------------------------------------------------
__kernel void eval_potential_brute(
    __global const float4* eval_points,
    __global const float4* ion_data,
    __global const float2* a_vec,
    __global const float2* b_vec,
    const int N_points,
    const int N_ions,
    const int N_rep,
    __global float* phi_out
){
    const int ip = get_global_id(0);
    if(ip >= N_points) return;

    float4 p = eval_points[ip];
    float3 r = (float3)(p.x, p.y, p.z);

    float2 a = a_vec[0];
    float2 b = b_vec[0];

    float phi = 0.0f;

    for(int n = -N_rep; n <= N_rep; n++){
        for(int m = -N_rep; m <= N_rep; m++){
            if(n*n + m*m > N_rep*N_rep) continue;

            float3 R = (float3)(n*a.x + m*b.x, n*a.y + m*b.y, 0.0f);

            for(int i = 0; i < N_ions; i++){
                float4 ion = ion_data[i];
                float3 ri = (float3)(ion.x, ion.y, ion.z);
                float q = ion.w;

                float3 dr = r - (ri + R);
                float r_mag = sqrt(dr.x*dr.x + dr.y*dr.y + dr.z*dr.z);

                if(r_mag > 1e-12f){
                    phi += q / r_mag;
                }
            }
        }
    }

    phi_out[ip] = phi * COULOMB_CONST;
}
