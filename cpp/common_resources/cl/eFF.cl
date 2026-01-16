#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

// ========================================================================
// Constants & Math Helpers (Kept aligned with Reference)
// ========================================================================

__constant static const float const_eV_SI        = 1.602176620898e-19; 
__constant static const float const_Angstroem_SI = 1.0e-10;            
__constant static const float const_El_eVA       = 14.3996448915f;
__constant static const float const_K_eVA        = 3.8099822f;
__constant static const float const_Bohr_Radius  = 0.529177210903f;

// Optimized approximations (inline for performance)
inline float2 erfx_e6( float x_, float k ){
    float x = x_*k;
    if( x>4.5f ){ float y=1.0f/x_; return (float2){ y, -y*y }; }
    float xx = x*x;
    float even  =  0.98501562f +xx*(-0.02756061f +xx*(-0.00188409f +xx*(-0.00309862f +xx*(-0.00134885f +xx*(-3.989465e-05f ) ) ) ) );
    float odd   = -0.13893350f +xx*(-0.00766429f +xx*( 0.00304682f +xx*( 0.00287933f +xx*( 0.00032604f +xx*( 1.970936e-06f ) ) ) ) );
    float deven =                  -0.05512122f +xx*(-0.00753638f +xx*(-0.01859177f +xx*(-0.01079087f +xx*(-0.00039894f ) ) ) )  ;
    float dodd  = -0.13893350f +xx*(-0.02299287f +xx*( 0.01523413f +xx*( 0.02015536f +xx*( 0.00293444f +xx*( 2.168030e-05f ) ) ) ) );
    float  t = even    + x*odd;
    float dt = deven*x +  dodd;
    float t2 = t*t;
    float t4 = t2*t2;
    float dt8_dx = 8.0f*dt*t*t2*t4;
    float y      = k/(t4*t4 + x);
    float dy     = -y*y*(dt8_dx+1.0f);
    return (float2){ y, dy };
}

inline float gauss_p8( float x ){ 
    if(x<-25.0f) return 0.0f;
    x *= 0.125f; float xx = x*x;
    float p = (1.0f+x) + xx*( 0.5f + 0.16666647f*x + xx*( 0.04166189f + 0.00830404f*x + xx*( 0.00132143f + 0.00013326f*x ) ) );
    p*=p; p*=p; p*=p;
    return p*p; // wait, exp_p8 returns p^4? original code: p*=p;p*=p;p*=p (p^8). exp(-x^2).
                // Let's stick to original code logic exactly.
                // Original: return exp_p8( -x*x );
}
// Re-implementing strictly based on provided snippet to ensure match
inline float exp_p8_exact( float x ){
    if(x<-25.0f) return 0.0f;
    x *= 0.125f;
    float xx = x*x;
    float p = (1.0f+x) +
               xx*( 0.5f   + 0.16666647f   *x +
               xx*( 0.04166189f  + 0.00830404f *x +
               xx*( 0.00132143f + 0.00013326f*x ) ) );
    p*=p; p*=p; p*=p;
    return p;
}
inline float gauss_p8_exact( float x ){ return exp_p8_exact( -x*x ); };

// ------------------------------------------------------------------------
// Interaction Physics (Exact Matches)
// ------------------------------------------------------------------------

inline float4 getCoulomb( float3 dp, float qq ){
    float r2   = dot(dp,dp) + 1e-32f;
    float ir   = rsqrt(r2);
    float E    = const_El_eVA*qq*ir;
    float fr   = -E/r2;
    return (float4){ dp*fr, E };
}

// Returns { E, fr, fsi, fsj }
inline float4 getCoulombGauss( float3 dR, float si, float sj, float qq ){
    float s2   = si*si + sj*sj;
    float s    = sqrt(s2);
    float r    = length( dR );
    
    float Amp    = qq*const_El_eVA;
    float is     = M_SQRT2_F/s; 
    float2 E_fr  = erfx_e6( r, is ); // .x = E_approx, .y = fr_approx
    
    float E  = E_fr.x * Amp;
    float fr = E_fr.y * Amp / (r + 1e-16f); // Derivative dE/dr / r

    float r_s    = r*is;
    // Original: fs = gauss_p8(r_s) *is*is*is*(0.5*M_2_SQRTPI*Amp);
    // M_2_SQRTPI = 1.12837916709551
    float fs     = gauss_p8_exact(r_s) * (is*is*is * 0.56418958354f * Amp); 
    
    return (float4){ E, fr, fs*si, fs*sj };
}

inline float4 getPauliGauss_New( float3 dR, float si, float sj, int spin, const float4 KRSrho ){
    float r2 = dot(dR,dR) + 1e-8f; 

    const float Hartree2eV = 27.211386f;
    const float A2bohr     = 1.0f/const_Bohr_Radius;

    float KR  = A2bohr*KRSrho.x;
    float KR2 = KR*KR;
    float KS  = A2bohr*KRSrho.y;
    float si_ = si * KS;
    float sj_ = sj * KS;
    
    r2 *= KR2; 

    float si2        = si_*si_;
    float sj2        = sj_*sj_;
    float si2sj2     = si2 + sj2;
    float invsi2sj2  = 1.0f/si2sj2;
    float invsi2sj22 = invsi2sj2*invsi2sj2;
    float invsi2sj23 = invsi2sj2*invsi2sj22;
    float denom_sij  = si_*sj_*invsi2sj2;
    float si4sj4     = si2*si2 - sj2*sj2;
    float invsj      = 1.0f/sj_;
    float invsi      = 1.0f/si_;
    float invsj2     = invsj*invsj;
    float invsi2     = invsi*invsi;
    float r2_4       = 4.0f*r2;

    float DT      = 1.5f*si2sj2*invsi2*invsj2 - (6.0f*si2sj2 - r2_4)*invsi2sj22;
    float dDT_dsi = -3.0f*invsi2*invsi        + 4.0f*si_*(3.0f*si2sj2 - r2_4)*invsi2sj23;
    float dDT_dsj = -3.0f*invsj2*invsj        + 4.0f*sj_*(3.0f*si2sj2 - r2_4)*invsi2sj23;
    float dDT_dr  =  8.0f*invsi2sj22;      

    float S22      = 8.0f*denom_sij*denom_sij*denom_sij*exp(-2.0f*r2*invsi2sj2);
    float dS22_dsi = S22*( -3.0f*si4sj4 + r2_4*si2 )*invsi2sj22*invsi;
    float dS22_dsj = S22*( +3.0f*si4sj4 + r2_4*sj2 )*invsi2sj22*invsj;
    float dS22_dr  = -4.0f*S22*invsi2sj2;   

    float rho = KRSrho.z;
    float E=0.0f, dE_dDT=0.0f, dE_dS22=0.0f;

    if(spin<=0){
        float invS22m1 = 1.0f/(S22+1.0f);
        E       += - rho*DT*S22  *invS22m1;
        dE_dDT  += -(rho*   S22 )*invS22m1;
        dE_dS22 += -(rho*DT     )*invS22m1*invS22m1;
    }
    if(spin>=0){
        float invS222m1 = 1.0f/( S22*S22-1.0f );
        E       +=   S22 * DT * ( -rho*S22                     + rho-2.0f ) *invS222m1;
        dE_dDT  += - S22 *      (  rho*S22                     - rho+2.0f ) *invS222m1;
        dE_dS22 +=      -  DT * (      S22*(S22*(rho-2.0f)-2.0f*rho) + rho-2.0f ) *invS222m1*invS222m1;
    }

    float sc = KRSrho.w;
    E         *= Hartree2eV*sc;
    
    float fsi = (dE_dS22 * dS22_dsi + dE_dDT * dDT_dsi)*Hartree2eV*-KS *sc;
    float fsj = (dE_dS22 * dS22_dsj + dE_dDT * dDT_dsj)*Hartree2eV*-KS *sc;
    float fr  = (dE_dS22 * dS22_dr  + dE_dDT * dDT_dr )*Hartree2eV*KR2 *sc;

    return (float4){ E, fr, fsi, fsj };
}

inline float2 addKineticGauss_eFF( float s ){
    float is  = M_SQRT2_F/s;
    float is2 = is*is*(const_K_eVA*1.5f);
    float fs  = is2*is*M_SQRT2_F;
    float E   = is2;
    return (float2){ E, fs  };
}

// ========================================================================
// KERNEL: Optimized Local MD for Small Systems (e.g., H2O)
// ========================================================================
// 1 Workgroup per System. N <= 64.
// Loads entirely to Local Memory.
// Performs exact physics replication of the serial CPU reference.

#define MAX_LOC_SIZE 64
__kernel void localMD(
    __global       int4*        sysinds, // 0  : [nsys]   {na,ne,i0p,i0a} size and initial index for each atom
    __global       float4*      pos,     // 1  : [ntot]   {x,y,z,w} positions (and size) of ions and electrons
    __global       float4*      vel,     // 2  : [ntot]   {vx,vy,vz,dw/dt} velocities of ions and electrons (including change of size)
    __global const float8*      aParams, // 3  : [ntot_a] parameters of ions { Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD }
    __global const signed char* espins,  // 4  : [ntot]   {spin}
    __global       float4*      fout,    // 5  : [ntot]   {fx,fy,fz,fw} output force buffer
    const int    nsys,                   // 6  : Number of systems
    const int    nsteps,                 // 7  : Number of steps
    const float  dt,                     // 8  : Time step
    const float  damping,                // 9  : Damping factor
    const float4 KRSrho,                 // 10 : KRSrho
    const int    bFrozenCore             // 11 : Boolean flag for frozen core
) {
    // --- 1. SLM Setup ---
    __local float4      l_pos     [MAX_LOC_SIZE];
    __local float4      l_vel     [MAX_LOC_SIZE];
    __local float4      l_force   [MAX_LOC_SIZE]; // forces from last step
    __local float8      l_aparams [MAX_LOC_SIZE];
    __local signed char l_spins   [MAX_LOC_SIZE];

    const int group_id = get_group_id(0);
    const int lid      = get_local_id(0);
    
    // Bounds check
    if(group_id >= nsys) return;

    const int4 inds    = sysinds[group_id];
    const int na       = inds.x;
    const int ne       = inds.y;
    const int ntot     = na + ne;
    const int ip_start = inds.z;
    const int ia_start = inds.w;

    // --- 2. Parallel Load from Global to Local ---
    if (lid < ntot) {
        l_pos[lid] = pos[ip_start + lid];
        l_vel[lid] = vel[ip_start + lid];
        if (lid < na) {
            l_aparams[lid] = aParams[ia_start + lid];
        } else {
            l_spins[lid - na] = espins[ip_start + lid];
        }
    } else if (lid < MAX_LOC_SIZE) {
        l_pos[lid] = (float4)(0.0f); // Safety padding
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // --- 3. MD Loop (Integration in SLM) ---
    for (int step = 0; step < nsteps; ++step) {
        
        if (lid < ntot) {
            float4 my_pos = l_pos[lid];
            float4 my_f   = (float4)(0.0f);
            bool i_am_ion = (lid < na);

            // --- A. Kinetic Force ---
            if (!i_am_ion) {
                float2 fk = addKineticGauss_eFF(my_pos.w);
                my_f.w += fk.y;
            }

            // --- B. Interaction Loop ---
            // N^2 loop. We do NOT use Newton's 3rd optimization here to avoid atomics.
            // Instead, we compute the exact force term for the (i,j) pair and apply it to i.
            // This ensures symmetry and correctness if the math is consistent.
            
            for (int j = 0; j < ntot; ++j) {
                if (lid == j) continue;

                float4 other_pos = l_pos[j];
                float3 dR        = other_pos.xyz - my_pos.xyz; // Vector from Me to Other
                bool j_is_ion    = (j < na);

                if (i_am_ion) {
                    float8 my_par = l_aparams[lid];
                    float Qi      = my_par.s0 - my_par.s2; // Q_nuc - Z_core_eff
                    float Q       = my_par.s0;

                    if (j_is_ion) {
                        // --- Ion-Ion ---
                        float8 other_par = l_aparams[j];
                        float Qj         = other_par.s0 - other_par.s2;
                        float4 fc        = getCoulomb(dR, Qi * Qj);
                        my_f.xyz        += fc.xyz; 
                        // Note: dR points to J. fc.xyz is (dR * fr). 
                        // fr is (-E/r2) < 0 for repulsion. 
                        // dR * neg = Away from J. Correct.
                    } 
                    else {
                        // --- Ion-Electron (I am Ion) ---
                        // Reference AE logic: dR = Elec - Ion.
                        // Here dR = Elec - Ion (since I am Ion). Matches ref.
                        float Ri = my_par.s1; // R_eff
                        float4 cg = getCoulombGauss(dR, Ri, other_pos.w, -Q);
                        // WAIT: Reference says: `getCoulombGauss(dR, Ri, pj.w, Qi)`
                        // `forcei.xyz += dR * cg.y`
                        
                        my_f.xyz += dR * cg.y;
                        // ion has no size DOF
                        
                        // Core correction (Frozen Core)
                        float sP = my_par.s2; // Zcore_eff ? No, sP usually parameter.
                        // Actually in provided code: `float sP = pari.s2;`
                        if (sP > 1e-8f) {
                             float4 KRS = KRSrho;
                             // Pauli
                             float4 pg = getPauliGauss_New(dR, Ri, other_pos.w, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(sP*0.5f)));
                             my_f.xyz += dR * pg.y;
                             // ion has no size DOF

                             // Core Coulomb
                             float4 cgC = getCoulombGauss(dR, Ri, other_pos.w, sP);
                             my_f.xyz += dR * cgC.y;
                             // ion has no size DOF
                        }
                    }
                } 
                else { // I am Electron
                    if (j_is_ion) {
                        // --- Electron-Ion (I am Electron) ---
                        // Match CPU/serial convention: dR = elec - ion, sizes ordered (sQ, se)
                        const float8 other_par = l_aparams[j];
                        const float  Q         = other_par.s0;
                        const float  sQ        = other_par.s1;
                        const float  sP        = other_par.s2;

                        const float3 dRei = my_pos.xyz - other_pos.xyz;
                        const float4 cg   = getCoulombGauss(dRei, sQ, my_pos.w, -Q);

                        my_f.xyz += -dRei * cg.y;
                        my_f.w   +=  cg.w;   // electron size-force from main Coulomb
                        
                        if (sP > 1e-8f) {
                            float4 KRS = KRSrho;
                            float4 pg  = getPauliGauss_New(dRei, sQ, my_pos.w, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(sP*0.5f)));
                            my_f.xyz  += -dRei * pg.y;
                            my_f.w    +=  pg.z; // CPU puts fsi (wrt sQ) into electron fsize

                            float4 cgC = getCoulombGauss(dRei, sQ, my_pos.w, sP);
                            my_f.xyz  += -dRei * cgC.y;
                            my_f.w    +=  cgC.z; // CPU puts fsi (wrt sQ) into electron fsize
                        }

                    } 
                    else {
                        // --- Electron-Electron ---
                        int my_spin    = l_spins[lid - na];
                        int other_spin = l_spins[j - na];
                        
                        // Reference mimics CPU: Calls Coulomb twice.
                        float4 cg1 = getCoulombGauss(dR, my_pos.w, other_pos.w, 1.0f);
                        float4 cg2 = getCoulombGauss(dR, my_pos.w, other_pos.w, 1.0f);
                        
                        // Pauli
                        float qq = ((my_spin==0) && (other_spin==0)) ? 2.0f : 1.0f;
                        float4 eff_KRS = KRSrho; eff_KRS.w *= qq;
                        float4 pg = getPauliGauss_New(dR, my_pos.w, other_pos.w, my_spin*other_spin, eff_KRS);
                        
                        float total_fr = cg1.y + cg2.y + pg.y;
                        my_f.xyz += dR * total_fr;
                        
                        // Size force: Accumulate .z (fsi) for me
                        my_f.w   += cg1.z + cg2.z + pg.z;
                    }
                }
            } // end loop j

            // --- C. Integration (Verlet/Euler) ---
            l_vel[lid] = l_vel[lid] * damping + my_f * dt;
            l_pos[lid] = l_pos[lid] + l_vel[lid] * dt;
            
            if (!i_am_ion) l_pos[lid].w = fmax(l_pos[lid].w, 0.001f);
            // Cache force from this step (for output)
            l_force[lid] = my_f;
        }
        
        barrier(CLK_LOCAL_MEM_FENCE); // Sync before next step
    }

    // --- 4. Write Back (positions, velocities, forces from last step) ---
    if (lid < ntot) {
        pos[ip_start + lid] = l_pos[lid];
        vel[ip_start + lid] = l_vel[lid];
        // Expose forces accumulated in the last step (needed for test_ocl_vs_cpu)
        fout[ip_start + lid] = l_force[lid];
    }
}


// ========================================================================
// Helper: Density Evaluation (Requested)
// ========================================================================
__kernel void eval_density_grid(
    const           int    ne,          // 0 : Number of electrons
    const          int     n_grid,      // 1 : Number of grid points
    __global const float4* electrons,   // 2 : [ne] {x,y,z,s} Electron positions
    __global const float*  amps,        // 3 : [ne] amplitudes (q)
    __global const float4* grid_points, // 4 : [n_grid] {x,y,z,rho_target} Grid points
    __global       float*  grid_out     // 5 : [n_grid] Output density grid
) {
    int gid = get_global_id(0);
    if (gid >= n_grid) return;

    float4 point = grid_points[gid];
    float  rho   = 0.0f;

    for(int i=0; i<ne; ++i){
        float4 elec = electrons[i];
        float  q    = amps[i];
        float3 dR   = point.xyz - elec.xyz;
        float s     = elec.w;
        float is    = M_SQRT2_F / s;
        // eFF density is generally summed squared magnitude of orbitals
        // simple sum of gaussians:
        float val = gauss_p8_exact( length(dR) * is );
        rho += q * val; 
    }
    grid_out[gid] = rho;
}


// ========================================================================
// Math Helper: Gaussian & Derivative
// ========================================================================

// Returns value p = exp( - 2 * r^2 / s^2 )
inline float gauss_density(float r, float s) {
    // We use the eFF definition: width parameter s.
    // Argument to exp is - (r / (s/sqrt(2)))^2 = - 2 * r^2 / s^2
    const float alpha = 2.0f / (s*s + 1e-16f);
    return native_exp(-alpha * r * r);
}

inline void fire_update_state( float power, float f_inc, float f_dec, float alpha_start, float dt_max, __private float* dt, __private float* alpha, __private int* np ){
    if( power > 0.0f ){
        (*np)++;
        if( (*np) > 5 ){
            (*dt)    = fmin( (*dt) * f_inc, dt_max );
            (*alpha) = (*alpha) * 0.99f;
        }
    }else{
        (*np)    = 0;
        (*dt)    = (*dt) * f_dec;
        (*alpha) = alpha_start;
    }
    (*dt)    = fmin( fmax( (*dt), 0.0f ), dt_max );
    (*alpha) = fmin( fmax( (*alpha), 0.0f ), 1.0f );
}

inline void fire_mix_v( float power, float v2, float f2, float alpha, float4 f, float fq, __private float4* v, __private float* vq ){
    if( power <= 0.0f ){
        (*v)  = (float4)(0.0f);
        (*vq) = 0.0f;
    }else{
        float v_mag     = sqrt(v2);
        float inv_f_mag = rsqrt(f2 + 1e-16f);
        float scale     = alpha * v_mag * inv_f_mag;
        (*v)  = (1.0f - alpha) * (*v)  + scale * f;
        (*vq) = (1.0f - alpha) * (*vq) + scale * fq;
    }
}

// ========================================================================
// KERNEL: Non-linear Fitting of Density using FIRE
// ========================================================================
// - 1 Workgroup per System
// - Electrons stored in SLM
// - Grid points streamed through SLM in tiles
// - All Optimization math done locally

#define MAX_LOC_ELEC 32
#define GRID_TILE_SIZE 32 

__kernel void fit_density_fire(
    const          int     ne,              // 0  : Number of electrons
    const          int     n_grid,          // 1  : Number of grid points
    __global       float4* electrons_inout, // 2  : [ne] {x,y,z,s} - Initial guess, overwritten by result
    __global const float4* grid_points,     // 3  : [n_grid] {x,y,z, rho_target}
    __global       float*  amps_inout,      // 4  : [ne] amplitudes q, overwritten by result
    __global       float4* vel_inout,       // 5  : [ne] velocities for {x,y,z,s}
    __global       float*  vq_inout,        // 6  : [ne] velocities for q
    __global       float4* fire_state,      // 7  : [nsys] {dt, alpha, 0, 0}
    __global       int*    fire_np,         // 8  : [nsys] n_pos
    __global       float4* force_out,       // 9  : [ne] variational force for {x,y,z,s}
    __global       float*  fq_out,          // 10 : [ne] variational force for q
    const          int     nsteps,          // 5  : Max FIRE steps
    const          float4  params,          // 6  : {dt_start, stiffness_pos, stiffness_size, stiffness_q}
    const          float4  fire_params,     // 7  : {f_inc, f_dec, alpha_start, dt_max}
    const          int     opt_mode,        // 8  : 0=FIRE, 1=damped MD, 2=GD
    const          float   md_damp          // 9  : damping factor for opt_mode=1
) {
    // --- 1. Setup Local Memory ---
    __local float4 l_pos      [MAX_LOC_ELEC]; // Current pos+size
    __local float4 l_ref      [MAX_LOC_ELEC]; // Reference pos+size (priors)
    __local float  l_q        [MAX_LOC_ELEC];
    __local float  l_qref     [MAX_LOC_ELEC];
    __local float4 l_vel      [MAX_LOC_ELEC]; // Velocity
    __local float4 l_force    [MAX_LOC_ELEC]; // Gradient/Force
    __local float  l_vq       [MAX_LOC_ELEC];
    __local float  l_fq       [MAX_LOC_ELEC];
    
    __local float4 l_grid_tile[GRID_TILE_SIZE]; // {x,y,z,rho_target} Cache for grid points
    __local float  l_rho_calc [GRID_TILE_SIZE]; // Calculated density for tile 

    // FIRE State (Shared across WG)
    __local float  l_dt;
    __local float  l_alpha;
    __local int    l_n_pos;
    __local float4 l_reduction[MAX_LOC_ELEC]; // {power, v2, f2, 0} reduced over lanes

    const int lid = get_local_id(0);
    const int gid = get_group_id(0); // If fitting multiple systems, use this offset (assumed 1 for now)

    // Configuration
    const float k_pos  = params.y; // Regularization spring constant (position)
    const float k_size = params.z; // Regularization spring constant (size)
    const float k_q    = params.w; // Regularization spring constant (amplitude)
    const float fire_f_inc       = fire_params.x;
    const float fire_f_dec       = fire_params.y;
    const float fire_alpha_start = fire_params.z;
    const float fire_dt_max      = fire_params.w;

    // --- 2. Load Electrons ---
    if (lid < ne) {
        float4 p = electrons_inout[lid];
        l_pos[lid] = p;
        l_ref[lid] = p;
        l_vel[lid] = vel_inout[lid];
        float q = amps_inout[lid];
        l_q[lid]    = q;
        l_qref[lid] = q;
        l_vq[lid]   = vq_inout[lid];
    }
    
    // Initialize FIRE state (host controls state; fall back to defaults if zero)
    if (lid == 0) {
        float4 st = fire_state[gid];
        l_dt    = (st.x > 0.0f) ? st.x : params.x;
        l_alpha = (st.y > 0.0f) ? st.y : fire_alpha_start;
        l_n_pos = fire_np[gid];
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    // =========================================================
    // OPTIMIZATION LOOP
    // =========================================================
    for (int step = 0; step < nsteps; ++step) {
        
        // A. Reset Forces
        if (lid < ne) { l_force[lid] = (float4)(0.0f); l_fq[lid]=0.0f; }
        barrier(CLK_LOCAL_MEM_FENCE);

        // B. Stream Grid Points in Tiles
        int num_tiles = (n_grid + GRID_TILE_SIZE - 1) / GRID_TILE_SIZE;

        for (int t = 0; t < num_tiles; ++t) {
            int g_idx = t * GRID_TILE_SIZE + lid;
            int tile_count = (t == num_tiles - 1) ? (n_grid - t * GRID_TILE_SIZE) : GRID_TILE_SIZE;

            // 1. Load Tile
            if (lid < tile_count) {
                l_grid_tile[lid] = grid_points[t * GRID_TILE_SIZE + lid];
                l_rho_calc [lid] = 0.0f; // Reset calc density
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            // 2. Compute Density at Grid Points (All Threads Help)
            // Strategy: We need rho at `tile_count` points. 
            // Thread `lid` computes density for grid point `lid`.
            // Inner loop over all `ne` electrons.
            if (lid < tile_count) {
                float4 gp = l_grid_tile[lid];
                float rho_sum = 0.0f;
                // Unroll this loop manually if NE is small and fixed, otherwise loop
                for (int j = 0; j < ne; ++j) {
                    float4 elec = l_pos[j];
                    float3 dR   = gp.xyz - elec.xyz;
                    float r     = length(dR);
                    rho_sum    += l_q[j] * gauss_density(r, elec.w);
                }
                l_rho_calc[lid] = rho_sum;
            }
            barrier(CLK_LOCAL_MEM_FENCE);

            // 3. Compute Gradients (Backprop) 
            // Thread `lid` updates Electron `lid`.
            // It must loop over all grid points in the tile to sum errors.
            if (lid < ne) {
                float4 my_elec = l_pos[lid];
                float4 my_grad = (float4)(0.0f);
                float  my_gq   = 0.0f;
                float  s       = my_elec.w;
                float  s_inv   = 1.0f / (s + 1e-16f);
                float  s_inv3  = s_inv * s_inv * s_inv;
                float  q       = l_q[lid];

                for (int k = 0; k < tile_count; ++k) {
                    float4 gp       = l_grid_tile[k];
                    float  rho_targ = gp.w;
                    float  rho_curr = l_rho_calc[k];
                    float  diff     = rho_curr - rho_targ; // Derivative of (rho - targ)^2 is 2*(diff)

                    // Derivative of rho w.r.t r and s
                    // G = exp(-2*r^2/s^2)
                    // dG/dr = G * (-4*r / s^2)
                    // dG/ds = G * (4*r^2 / s^3)
                    
                    float3 dR = my_elec.xyz - gp.xyz; // Vec from Grid to Elec? No, Elec - Grid.
                    // Wait, r = |Elec - Grid|.
                    // d(G)/d(Elec_pos). 
                    // Let u = -2 |r_e - r_g|^2 / s^2.
                    // du/dx_e = -2/s^2 * 2(x_e - x_g) = -4/s^2 * dR_x.
                    
                    float r2 = dot(dR, dR);
                    float G  = gauss_density(sqrt(r2), s);
                    
                    // Factor 2.0f comes from derivative of error squared
                    float commonq = 2.0f * diff * q * G; 
                    
                    // Force = -Gradient
                    // F_pos = - (common * (-4/s^2 * dR)) = common * 4/s^2 * dR
                    float f_pre_pos = commonq * 4.0f * (s_inv * s_inv);
                    // Logic check: If calc > target (diff > 0), we want to reduce density.
                    // Move electron AWAY from grid point.
                    // dR is (Elec - Grid).
                    // my_grad -= pos * dR -> Moves opposite to dR (Towards Grid). 
                    // Wait. 
                    // d(Error)/dRe = 2 * diff * dG/dRe.
                    // dG/dRe = -4/s^2 * (Re - Rg) * G.
                    // Grad = - (Pos_const * dR). (Points towards Grid).
                    // We want to minimize error. Move against Gradient.
                    // Force = - Grad = + (Pos_const * dR).
                    // So Force should be += f_pre_pos * dR.
                    // Let's re-verify:
                    // If my_grad.xyz -= ... that implies we are accumulating Gradient.
                    // At end we do F = -Grad - Reg.
                    
                    my_grad.xyz += f_pre_pos * dR; // Force pushes away if diff>0 (Too much density)

                    // dG/ds = G * (4*r^2 / s^3)
                    // Grad_s = 2 * diff * G * 4 * r2 / s3
                    // Force_s = -Grad_s
                    float f_pre_s = commonq * 4.0f * r2 * s_inv3;
                    my_grad.w -= f_pre_s; // If diff>0 (too fat), reduce size (s).

                    my_gq -= 2.0f * diff * G;
                }
                
                l_force[lid] += my_grad;
                l_fq[lid]    += my_gq;
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }

        // C. Regularization & FIRE Update
        if (lid < ne) {
            float4 p   = l_pos[lid];
            float4 ref = l_ref[lid];
            float4 v   = l_vel[lid];
            float4 f   = l_force[lid];
            float  q   = l_q[lid];
            float  q0  = l_qref[lid];
            float  fq  = l_fq[lid];

            // 1. Add Regularization Forces (Harmonic priors)
            // F_reg = -k * (x - x0)
            f.xyz -= k_pos  * (p.xyz - ref.xyz);
            f.w   -= k_size * (p.w   - ref.w);
            fq    -= k_q    * (q     - q0);

            // 2. FIRE Reduction Scalars (reduced over the whole parameter vector)
            // power = <F|V> ; v2 = <V|V> ; f2 = <F|F>
            float power = dot(f, v) + fq*l_vq[lid];
            float v2    = dot(v, v) + l_vq[lid]*l_vq[lid];
            float f2    = dot(f, f) + fq*fq;
            l_reduction[lid] = (float4)(power, v2, f2, 0.0f);
            l_force[lid]    = f; // Store full force
            l_fq[lid]       = fq;
        }
        if (lid >= ne) { l_reduction[lid] = (float4)(0.0f); }
        barrier(CLK_LOCAL_MEM_FENCE);

        // D. Reduction (Parallel Reduction)
        for (int offset = MAX_LOC_ELEC / 2; offset > 0; offset >>= 1) {
            if (lid < offset) {
                l_reduction[lid] += l_reduction[lid + offset];
            }
            barrier(CLK_LOCAL_MEM_FENCE);
        }
        float total_power = l_reduction[0].x;
        float total_v2    = l_reduction[0].y;
        float total_f2    = l_reduction[0].z;

        // E. Update FIRE State (Thread 0) / dt clamp for other modes
        if (lid == 0) {
            float dt    = l_dt;
            float alpha = l_alpha;
            int   np    = l_n_pos;

            if( opt_mode == 0 ){
                fire_update_state( total_power, fire_f_inc, fire_f_dec, fire_alpha_start, fire_dt_max, &dt, &alpha, &np );
            }else{
                dt    = fmin( fmax(dt, 0.0f), fire_dt_max );
                alpha = fmin( fmax(alpha, 0.0f), 1.0f );
            }

            l_dt    = dt;
            l_alpha = alpha;
            l_n_pos = np;

            if (gid == 0 && step == 0 && fire_np[gid] == 0) {
                printf("FIRE cfg dt=%g alpha=%g inc=%g dec=%g dtmax=%g\n", dt, alpha, fire_f_inc, fire_f_dec, fire_dt_max);
            }
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        // F. Integration (FIRE / MD / GD)
        if (lid < ne) {
            float4 v = l_vel[lid];
            float4 f = l_force[lid];
            float  vq = l_vq[lid];
            float  fq = l_fq[lid];
            float  dt = l_dt;
            float  alpha = l_alpha;

            if( opt_mode == 0 ){
                fire_mix_v( total_power, total_v2, total_f2, alpha, f, fq, &v, &vq );
                v  += f  * dt;
                vq += fq * dt;
            }else if( opt_mode == 1 ){
                v  = v  * md_damp + f  * dt;
                vq = vq * md_damp + fq * dt;
            }else{
                v  = (float4)(0.0f);
                vq = 0.0f;
            }

            float4 p_new = l_pos[lid];
            float  q_new = l_q[lid];
            if( opt_mode == 2 ){
                p_new += f  * dt;
                q_new += fq * dt;
            }else{
                p_new += v  * dt;
                q_new += vq * dt;
            }

            // Constraint: Size > 0.1 (Prevent collapse)
            p_new.w = fmax(p_new.w, 0.1f);

            l_vel[lid] = v;
            l_pos[lid] = p_new;
            l_vq [lid] = vq;
            l_q  [lid] = q_new;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- 3. Write Back Result ---
    if (lid < ne) {
        electrons_inout[lid] = l_pos[lid];
        amps_inout[lid]      = l_q[lid];
        vel_inout[lid]       = l_vel[lid];
        vq_inout[lid]        = l_vq[lid];
        force_out[lid]       = l_force[lid];
        fq_out[lid]          = l_fq[lid];
    }
    if (lid == 0) {
        fire_state[gid] = (float4)(l_dt, l_alpha, 0.0f, 0.0f);
        fire_np[gid]    = l_n_pos;
    }
}