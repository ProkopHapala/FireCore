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
        float denom = (S22*S22 - 1.0f);
        if(fabs(denom) < 1e-12f){ denom = (denom>=0.0f) ? 1e-12f : -1e-12f; }
        float invS222m1 = 1.0f/denom;
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
#ifndef DBG_EFF_PAIR
#define DBG_EFF_PAIR 0
#endif
#ifndef IDBG_SYS
#define IDBG_SYS 0
#endif
 #ifndef IDBG_STEP
 #define IDBG_STEP 0
 #endif
 #ifndef IDBG_I
 #define IDBG_I 0
 #endif
 #ifndef IDBG_J
 #define IDBG_J 1
 #endif
 #ifndef DBG_EFF_ALLPAIRS
 #define DBG_EFF_ALLPAIRS 0
 #endif
#ifndef DBG_EFF_INPUT
#define DBG_EFF_INPUT 0
#endif
#ifndef DBG_EFF_FDECOMP
#define DBG_EFF_FDECOMP 0
#endif
__kernel void localMD(
    __global       int4*        sysinds, // 0  : [nsys]   {na,ne,i0p,i0a} size and initial index for each atom
    __global       float4*      pos,     // 1  : [ntot]   {x,y,z,w} positions (and size) of ions and electrons
    __global       float4*      vel,     // 2  : [ntot]   {vx,vy,vz,dw/dt} velocities of ions and electrons (including change of size)
    __global const float8*      aParams, // 3  : [ntot_a] parameters of ions { Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD }
    __global const signed char* espins,  // 4  : [ntot]   {spin}
    __global       float4*      fout,    // 5  : [ntot]   {fx,fy,fz,fw} output force buffer
    __global const  uchar*      fixmask, // 6  : [ntot]   bitmask of fixed DOFs: 1=x 2=y 4=z 8=w(size)
    const int    nsys,                   // 7  : Number of systems
    const int    nsteps,                 // 8  : Number of steps
    const float  dt,                     // 9  : Time step
    const float  damping,                // 10 : Damping factor
    const float4 KRSrho,                 // 11 : KRSrho
    const int    bFrozenCore             // 12 : Boolean flag for frozen core
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

    const int bDbg = (DBG_EFF_PAIR && (group_id==IDBG_SYS));

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
            float4 dbg_fAE = (float4)(0.0f);
            float4 dbg_fEE = (float4)(0.0f);
            float4 dbg_fK  = (float4)(0.0f);
            bool i_am_ion = (lid < na);
            const uchar msk = fixmask ? fixmask[ip_start + lid] : (uchar)0;

            if ((get_global_id(0)==0) && (lid==0) && (step==0)) {
                printf("SWITCH localMD: nsys=%d group=%d na=%d ne=%d ntot=%d bFrozenCore=%d KRSrho=(%.8g %.8g %.8g %.8g)\n",
                    nsys, group_id, na, ne, ntot, bFrozenCore,
                    (double)KRSrho.x, (double)KRSrho.y, (double)KRSrho.z, (double)KRSrho.w );
#if DBG_EFF_INPUT
                // DEBUG: Dump all particle inputs (compile-time gated)
                for (int dd=0; dd<ntot; dd++){
                    float4 pp = l_pos[dd];
                    if (dd < na){
                        float8 ap = l_aparams[dd];
                        printf("GPU_INPUT ion[%d] pos(%12.6f %12.6f %12.6f) aP(Q=%g sQ=%g sP=%g s3=%g s4=%g s5=%g s6=%g s7=%g)\n",
                            dd, (double)pp.x,(double)pp.y,(double)pp.z,
                            (double)ap.s0,(double)ap.s1,(double)ap.s2,(double)ap.s3,(double)ap.s4,(double)ap.s5,(double)ap.s6,(double)ap.s7);
                    } else {
                        int sp = l_spins[dd - na];
                        printf("GPU_INPUT elec[%d] pos(%12.6f %12.6f %12.6f) size=%12.8f spin=%d\n",
                            dd, (double)pp.x,(double)pp.y,(double)pp.z, (double)pp.w, sp);
                    }
                }
#endif
            }

            // --- A. Kinetic Force ---
            if (!i_am_ion) {
                float2 fk = addKineticGauss_eFF(my_pos.w);
                my_f.w += fk.y;
                dbg_fK.w += fk.y;
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
                    float Qi      = my_par.s0 - ( (bFrozenCore!=0) ? my_par.s2 : 0.0f ); // Q_nuc - Z_core_eff (only in frozen-core mode)
                    float Q       = my_par.s0;

                    if (j_is_ion) {
                        // --- Ion-Ion ---
                        float8 other_par = l_aparams[j];
                        float Qj         = other_par.s0 - ( (bFrozenCore!=0) ? other_par.s2 : 0.0f );
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
                        const int   se = l_spins[j - na];
                        const float qe = 1.0f;
                        float4 cg = getCoulombGauss(dR, Ri, other_pos.w, -Q*qe);
                        my_f.xyz += dR * cg.y;
                        // ion has no size DOF
                        
                        // Core correction (Frozen Core)
                        float sP = (bFrozenCore!=0) ? my_par.s2 : 0.0f;
                        if (sP > 1e-8f) {
                             float4 KRS = KRSrho;
                            // Pauli
                            float4 pg = getPauliGauss_New(dR, Ri, other_pos.w, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(qe*sP*0.5f)));
                            my_f.xyz += dR * pg.y;
                            // ion has no size DOF

                             // Core Coulomb
                            float4 cgC = getCoulombGauss(dR, Ri, other_pos.w, qe*sP);
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
                        const float  sP        = (bFrozenCore!=0) ? other_par.s2 : 0.0f;

                        const float3 dRei = my_pos.xyz - other_pos.xyz;
                        const int    se   = l_spins[lid - na];
                        const float  qe   = 1.0f;
                        const float4 cg   = getCoulombGauss(dRei, sQ, my_pos.w, -Q*qe);

                        my_f.xyz += -dRei * cg.y;
                        my_f.w   += cg.w;   // size-force for electron (fsj)
                        dbg_fAE.xyz += -dRei * cg.y;
                        dbg_fAE.w   += cg.w;
                        
                        if (sP > 1e-8f) {
                            float4 KRS = KRSrho;
                            float4 pg  = getPauliGauss_New(dRei, sQ, my_pos.w, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(qe*sP*0.5f)));
                            my_f.xyz  += -dRei * pg.y;
                            my_f.w    += pg.z;
                            dbg_fAE.xyz += -dRei * pg.y;
                            dbg_fAE.w   += pg.z;

                            float4 cgC = getCoulombGauss(dRei, sQ, my_pos.w, sP*qe);
                            my_f.xyz  += -dRei * cgC.y;
                            my_f.w    += cgC.z;
                            dbg_fAE.xyz += -dRei * cgC.y;
                            dbg_fAE.w   += cgC.z;
                        }

                        if (bDbg && (step==IDBG_STEP) && ( (DBG_EFF_ALLPAIRS && (lid==0)) || ((!DBG_EFF_ALLPAIRS) && (lid==IDBG_I) && (j==IDBG_J)) ) ) {
                            const float r = sqrt(dot(dRei,dRei));
                            float4 pg_dbg = (float4)(0.0f);
                            float4 cgC_dbg = (float4)(0.0f);
                            if (sP > 1e-8f) {
                                float4 KRS = KRSrho;
                                pg_dbg  = getPauliGauss_New(dRei, sQ, my_pos.w, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(qe*sP*0.5f)));
                                cgC_dbg = getCoulombGauss(dRei, sQ, my_pos.w, sP*qe);
                            }
                            const float3 fC  = -dRei * cg.y;
                            const float3 fP  = -dRei * pg_dbg.y;
                            const float3 fCC = -dRei * cgC_dbg.y;
                            printf("DBG_PAIR step=%d kind=AE(iElec) i=%d j=%d r=%.8g sQ=%.8g se=%.8g Qion=%.8g sP=%.8g | Ec=%.10g Ep=%.10g EcC=%.10g | fC(%.6g %.6g %.6g) fsC=%.6g | fP(%.6g %.6g %.6g) fsP=%.6g | fCC(%.6g %.6g %.6g) fsCC=%.6g\n",
                                step, lid, j, (double)r, (double)sQ, (double)my_pos.w, (double)Q, (double)sP,
                                (double)cg.x, (double)pg_dbg.x, (double)cgC_dbg.x,
                                (double)fC.x, (double)fC.y, (double)fC.z, (double)cg.w,
                                (double)fP.x, (double)fP.y, (double)fP.z, (double)pg_dbg.z,
                                (double)fCC.x,(double)fCC.y,(double)fCC.z,(double)cgC_dbg.z );
                        }
                    }
 
                    else {
                        // --- Electron-Electron ---
                        int my_spin    = l_spins[lid - na];
                        int other_spin = l_spins[j - na];
                        const float qij = 1.0f;
                        
                        // Per-particle force accumulation: this work-item computes force on i.
                        // Therefore each pair contribution must be added once (CPU i<j loop adds to both i and j).
                        float4 cg1 = getCoulombGauss(dR, my_pos.w, other_pos.w, qij);
                        
                        // Pauli
                        float qq = ((my_spin==0) && (other_spin==0)) ? 2.0f : 1.0f;
                        float4 eff_KRS = KRSrho; eff_KRS.w *= qq;
                        float4 pg = getPauliGauss_New(dR, my_pos.w, other_pos.w, my_spin*other_spin, eff_KRS);
                        
                        float total_fr = cg1.y + pg.y;
                        my_f.xyz +=  dR * total_fr;
                        
                        // Size force: Accumulate .z (fsi) for me
                        my_f.w   += cg1.z + pg.z;
                        dbg_fEE.xyz += dR * total_fr;
                        dbg_fEE.w   += cg1.z + pg.z;

                        if (bDbg && (step==IDBG_STEP) && ( (DBG_EFF_ALLPAIRS && (lid==0)) || ((!DBG_EFF_ALLPAIRS) && (lid==IDBG_I) && (j==IDBG_J)) ) ) {
                            const float r = sqrt(dot(dR,dR));
                            printf("DBG_PAIR step=%d kind=EE i=%d j=%d r=%.8g si=%.8g sj=%.8g spinij=%d qq=%.8g | Ec1=%.10g Ec2=%.10g Ep=%.10g fr=(%.10g %.10g %.10g)\n",
                                step, lid, j, (double)r, (double)my_pos.w, (double)other_pos.w, (int)(my_spin*other_spin), (double)qq,
                                (double)cg1.x, 0.0, (double)pg.x, (double)dR.x, (double)dR.y, (double)dR.z );
                        }
                    }
                }
            } // end loop j

            if (bDbg && (step==IDBG_STEP) && (lid==IDBG_I)) {
                printf("DBG_SUM step=%d i=%d kind=SUM | fTot(%.8g %.8g %.8g %.8g) fAE(%.8g %.8g %.8g %.8g) fEE(%.8g %.8g %.8g %.8g) fK(%.8g %.8g %.8g %.8g)\n",
                    step, lid,
                    (double)my_f.x,(double)my_f.y,(double)my_f.z,(double)my_f.w,
                    (double)dbg_fAE.x,(double)dbg_fAE.y,(double)dbg_fAE.z,(double)dbg_fAE.w,
                    (double)dbg_fEE.x,(double)dbg_fEE.y,(double)dbg_fEE.z,(double)dbg_fEE.w,
                    (double)dbg_fK.x,(double)dbg_fK.y,(double)dbg_fK.z,(double)dbg_fK.w
                );
            }
            // DEBUG: print per-electron force decomposition for ALL electrons at step 0 (compile-time gated)
#if DBG_EFF_FDECOMP
            if ((step==0) && (!i_am_ion) && (group_id==0)) {
                printf("GPU_FDECOMP elec[%d] fTot(%12.6f %12.6f %12.6f %12.6f) fAE(%12.6f %12.6f %12.6f %12.6f) fEE(%12.6f %12.6f %12.6f %12.6f) fK(%12.6f %12.6f %12.6f %12.6f)\n",
                    lid,
                    (double)my_f.x,(double)my_f.y,(double)my_f.z,(double)my_f.w,
                    (double)dbg_fAE.x,(double)dbg_fAE.y,(double)dbg_fAE.z,(double)dbg_fAE.w,
                    (double)dbg_fEE.x,(double)dbg_fEE.y,(double)dbg_fEE.z,(double)dbg_fEE.w,
                    (double)dbg_fK.x,(double)dbg_fK.y,(double)dbg_fK.z,(double)dbg_fK.w);
            }
#endif

            // --- C. Integration (Verlet/Euler) ---
            float4 vnew = l_vel[lid] * damping + my_f * dt;
            float4 pnew = l_pos[lid] + vnew * dt;

            // Apply fixmask (simple, per-DOF)
            if (msk) {
                if (msk & 1) { vnew.x = 0.0f; pnew.x = l_pos[lid].x; my_f.x = 0.0f; }
                if (msk & 2) { vnew.y = 0.0f; pnew.y = l_pos[lid].y; my_f.y = 0.0f; }
                if (msk & 4) { vnew.z = 0.0f; pnew.z = l_pos[lid].z; my_f.z = 0.0f; }
                if (msk & 8) { vnew.w = 0.0f; pnew.w = l_pos[lid].w; my_f.w = 0.0f; }
            }

            l_vel[lid] = vnew;
            l_pos[lid] = pnew;
            
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


#ifdef EFF_EVAL_ENERGY
// ========================================================================
// KERNEL: localMD_energy (no-relax energy eval using localMD data path)
// ========================================================================
// Conservative: does NOT modify localMD. This is a separate kernel with identical
// force path signature, plus an energy output buffer.
//
// Intended use: nsteps=1, dt=0, damping=0. We compute energies from positions only.
__kernel void localMD_energy(
    __global       int4*        sysinds, // 0  : [nsys]   {na,ne,i0p,i0a}
    __global       float4*      pos,     // 1  : [ntot]   {x,y,z,w}
    __global       float4*      vel,     // 2  : [ntot]   {vx,vy,vz,dw/dt}
    __global const float8*      aParams, // 3  : [ntot_a] parameters of ions { Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD }
    __global const signed char* espins,  // 4  : [ntot]   {spin}
    __global       float4*      fout,    // 5  : [ntot]   {fx,fy,fz,fw} output force buffer (kept for compatibility)
    __global       float8*      Es,      // 6  : [nsys]   {Etot,Ek,Eee,Eae,Eaa,0,0,0}
    __global const   char*      coreMap,
    __global const  uchar*      fixmask,
    const int    nsys,                   // 7
    const int    nsteps,                 // 8  (unused; kept for interface compatibility)
    const float  dt,                     // 9  (unused)
    const float  damping,                // 10 (unused)
    const float4 KRSrho,                 // 11
    const int    bFrozenCore,            // 12
    const int    bOffloadCore
){
    __local float4      l_pos     [MAX_LOC_SIZE];
    __local float8      l_aparams [MAX_LOC_SIZE];
    __local signed char l_spins   [MAX_LOC_SIZE];

    __local float l_Ek   [MAX_LOC_SIZE];
    __local float l_Eee  [MAX_LOC_SIZE];
    __local float l_Eae  [MAX_LOC_SIZE];
    __local float l_Eaa  [MAX_LOC_SIZE];
    __local float l_EeeP [MAX_LOC_SIZE];
    __local float l_EaeP [MAX_LOC_SIZE];

    const int group_id = get_group_id(0);
    const int lid      = get_local_id(0);
    if(group_id >= nsys) return;

    const int4 inds    = sysinds[group_id];
    const int na       = inds.x;
    const int ne       = inds.y;
    const int ntot     = na + ne;
    const int ip_start = inds.z;
    const int ia_start = inds.w;

    if (lid < ntot) {
        l_pos[lid] = pos[ip_start + lid];
        if (lid < na) {
            l_aparams[lid] = aParams[ia_start + lid];
        } else {
            l_spins[lid - na] = espins[ip_start + lid];
        }
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    double Ek      = 0.0;
    float EeeCoul = 0.0f;
    float EeePaul = 0.0f;
    float EaeCoul = 0.0f;
    float EaePaul = 0.0f;
    float Eaa     = 0.0f;

    if(lid < ntot){
        const float4 pi    = l_pos[lid];
        const int    i_ion = (lid < na);

        // Kinetic (electrons only)  (C++ scales by electron charge magnitude)
        if(!i_ion){
            const int si = l_spins[lid - na];
            const float qe = (si==0) ? 2.0f : 1.0f;
            float2 fk = addKineticGauss_eFF(pi.w);
            if( !(bOffloadCore && (coreMap[ip_start + lid] >= 0)) ){
                Ek += (double)(fk.x * qe);
            }
        }

        for(int j=lid+1; j<ntot; j++){
            const float4 pj    = l_pos[j];
            const int    j_ion = (j < na);
            const float3 dR    = pj.xyz - pi.xyz;

            if(i_ion && j_ion){
                const float8 pari = l_aparams[lid];
                const float8 parj = l_aparams[j];
                const float Qi = pari.s0 - ( (bFrozenCore!=0) ? pari.s2 : 0.0f );
                const float Qj = parj.s0 - ( (bFrozenCore!=0) ? parj.s2 : 0.0f );
                float4 fc = getCoulomb(dR, Qi*Qj);
                Eaa += fc.w;
            }else if(i_ion != j_ion){
                // electron-ion (order: dR_ei = epos - ipos)
                const int   ia = i_ion ? lid : j;
                const int   ie = i_ion ? j   : lid;
                const float4 pa = l_pos[ia];
                const float4 pe = l_pos[ie];
                const float3 dRei = pe.xyz - pa.xyz;

                const signed char se = l_spins[ie - na];
                const float qe = (se==0) ? 2.0f : 1.0f;

                const float8 par = l_aparams[ia];
                const float  Q   = par.s0;
                const float  sQ  = par.s1;
                const float  sP  = (bFrozenCore!=0) ? par.s2 : 0.0f;

                float4 cg = getCoulombGauss(dRei, sQ, pe.w, -Q*qe);
                EaeCoul += cg.x;
                if(sP > 1e-8f){
                    float4 KRS = KRSrho;
                    // C++: addPauliGauss_New(..., KRSrho, qj*sP*0.5)
                    float4 pg  = getPauliGauss_New(dRei, sQ, pe.w, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(qe*sP*0.5f)));
                    EaePaul += pg.x;
                    // C++ counts core-correction Coulomb into Eee (not Eae)
                    float4 cgC = getCoulombGauss(dRei, sQ, pe.w, sP*qe);
                    EeeCoul += cgC.x;
                }
            }else{
                // electron-electron
                const int si = l_spins[lid - na];
                const int sj = l_spins[j   - na];
                const float qi = (si==0) ? 2.0f : 1.0f;
                const float qj = (sj==0) ? 2.0f : 1.0f;
                // Match CPU energy bookkeeping: Coulomb energy counted once
                float4 cg = getCoulombGauss(dR, pi.w, pj.w, qi*qj);
                float qq = ((si==0) && (sj==0)) ? 2.0f : 1.0f;
                float4 eff_KRS = KRSrho; eff_KRS.w *= qq;
                float4 pg = getPauliGauss_New(dR, pi.w, pj.w, si*sj, eff_KRS);
                EeeCoul += cg.x;
                EeePaul += pg.x;
            }
        }
    }

    // reduction
    if(lid < MAX_LOC_SIZE){
        l_Ek  [lid] = (lid < ntot) ? Ek      : 0.0f;
        l_Eee [lid] = (lid < ntot) ? EeeCoul : 0.0f;
        l_Eae [lid] = (lid < ntot) ? EaeCoul : 0.0f;
        l_Eaa [lid] = (lid < ntot) ? Eaa     : 0.0f;
        l_EeeP[lid] = (lid < ntot) ? EeePaul : 0.0f;
        l_EaeP[lid] = (lid < ntot) ? EaePaul : 0.0f;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    for(int stride=MAX_LOC_SIZE/2; stride>0; stride>>=1){
        if(lid < stride){
            l_Ek [lid] += l_Ek [lid+stride];
            l_Eee[lid] += l_Eee[lid+stride];
            l_Eae[lid] += l_Eae[lid+stride];
            l_Eaa[lid] += l_Eaa[lid+stride];
            l_EeeP[lid] += l_EeeP[lid+stride];
            l_EaeP[lid] += l_EaeP[lid+stride];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(lid == 0){
        double Etot = l_Ek[0] + l_Eee[0] + l_EeeP[0] + l_Eae[0] + l_EaeP[0] + l_Eaa[0];
        Es[group_id] = (float8)(Etot, l_Ek[0], l_Eee[0], l_Eae[0], l_Eaa[0], l_EeeP[0], l_EaeP[0], 0.0f);
    }
}
#endif


// ========================================================================
// KERNEL: Energy evaluation (no relaxation)
// ========================================================================
// One workgroup per system. Only lid==0 does the work (serial per system).
// Outputs: {Etot,Ek,Eee,Eae,Eaa} for each system.
__kernel void evalEnergy(
    __global const int4*        sysinds,   // 0  : [nsys] {na,ne,i0p,i0a}
    __global const float4*      pos,       // 1  : [ntot] {x,y,z,w}
    __global const float8*      aParams,   // 2  : [ntot_a]
    __global const signed char* espins,    // 3  : [ntot]
    __global       float8*      Es,        // 4  : [nsys] {Etot,Ek,Eee,Eae,Eaa,0,0,0}
    __global const   char*      coreMap,
    const int                 nsys,       // 5
    const float4              KRSrho,     // 6
    const int                 bFrozenCore,// 7
    const int                 bOffloadCore
){
    const int group_id = get_group_id(0);
    const int lid      = get_local_id(0);
    if(group_id >= nsys) return;
    if(lid != 0) return;

    const int4 inds    = sysinds[group_id];
    const int na       = inds.x;
    const int ne       = inds.y;
    const int ntot     = na + ne;
    const int ip_start = inds.z;
    const int ia_start = inds.w;

    const int bDbg = (DBG_EFF_PAIR && (group_id==IDBG_SYS));

    double Ek      = 0.0;
    float EeeCoul = 0.0f; // matches C++ ff.Eee (includes core-correction Coulomb)
    float EeePaul = 0.0f; // matches C++ ff.EeePaul
    float EaeCoul = 0.0f; // matches C++ ff.Eae
    float EaePaul = 0.0f; // matches C++ ff.EaePaul
    float Eaa     = 0.0f;

    // Kinetic (electrons only)  (C++ scales by electron charge magnitude)
    for(int ie=0; ie<ne; ie++){
        const int ip = ip_start + na + ie;
        const float4 pi = pos[ip];
        const float si = pi.w;
        const char spin_i = espins[ip];
        const float qe = (spin_i==0) ? 2.0f : 1.0f;
        if( !(bOffloadCore && (coreMap[ip] >= 0)) ){
            Ek += (double)( addKineticGauss_eFF(si).x * qe );
        }
    }

    // AA
    for(int i=0; i<na; i++){
        float4 pi   = pos[ip_start + i];
        float8 pari = aParams[ia_start + i];
        float Qi    = pari.s0 - ( (bFrozenCore!=0) ? pari.s2 : 0.0f );
        for(int j=i+1; j<na; j++){
            float4 pj   = pos[ip_start + j];
            float8 parj = aParams[ia_start + j];
            float Qj    = parj.s0 - ( (bFrozenCore!=0) ? parj.s2 : 0.0f );
            float3 dR   = pj.xyz - pi.xyz;
            float4 fc = getCoulomb( dR, Qi*Qj );
            Eaa += (double)fc.w;
            if(bDbg){ printf("GPU AA(%d,%d) r=%g E=%g\n", i, j, sqrt((double)(dot(dR,dR))), (double)fc.w ); }
        }
    }

    // AE + (optional frozen-core correction terms)
    for(int ia=0; ia<na; ia++){
        const int ia_glb = ia_start + ia;
        const float4 ai = pos[ip_start + ia];
        const float8 pari = aParams[ia_glb];
        const float3 ri = ai.xyz;
        const float  Q  = pari.s0;
        const float  sQ = pari.s1;
        const float  sP = pari.s2;
        const float  Qi = Q - sP;
        for(int ie=0; ie<ne; ie++){
            const int ip = ip_start + na + ie;
            const float4 pi = pos[ip];
            const float3 re = pi.xyz;
            const float  se = pi.w;
            const char spin_e = espins[ip];
            const float qe = (spin_e==0) ? 2.0f : 1.0f;
            const float3 dR = re - ri;
            const char cme = coreMap[ip];
            const int skip_core_self = (bOffloadCore && (cme == (char)ia));
            if(!skip_core_self){
                float4 cg = getCoulombGauss( dR, sQ, se, -Q*qe );
                EaeCoul += (double)cg.x;

                float4 pg = (float4)(0.0f);
                float4 cgC = (float4)(0.0f);
                if(sP > 1e-8f){
                    float4 KRS = KRSrho;
                    pg  = getPauliGauss_New( dR, sQ, se, 0, (float4)(KRS.x, KRS.y, KRS.z, KRS.w*(sP*0.5f*qe)) );
                    cgC = getCoulombGauss( dR, sQ, se, sP*qe );
                    EaePaul += (double)pg.x;
                    EeeCoul += (double)cgC.x;
                }

                if(bDbg){
                    printf("GPU AE(a=%d,e=%d) r=%g se=%g sQ=%g qe=%g Ec=%g Ep=%g EcoreC=%g\n", ia, ie, sqrt((double)(dot(dR,dR))), (double)se, (double)sQ, (double)qe, (double)cg.x, (double)pg.x, (double)cgC.x );
                }
            }else{
                if(bDbg){ printf("GPU AE_SKIP_CORE_SELF(a=%d,e=%d)\n", ia, ie ); }
            }
        }
    }

    // EE
    for(int ie=0; ie<ne; ie++){
        const int ipi = ip_start + na + ie;
        const float4 pi = pos[ipi];
        const float3 ri = pi.xyz;
        const float  si = pi.w;
        const char spin_i = espins[ipi];
        const float qe_i = (spin_i==0) ? 2.0f : 1.0f;
        for(int je=0; je<ie; je++){
            const int ipj = ip_start + na + je;
            if(bOffloadCore){
                const char cmi = coreMap[ipi];
                const char cmj = coreMap[ipj];
                if( (cmi>=0) && (cmi==cmj) ){
                    if(bDbg){ printf("GPU EE_SKIP_CORE_SELF(%d,%d) ia=%d\n", ie, je, (int)cmi ); }
                    continue;
                }
            }
            const float4 pj = pos[ipj];
            const float3 dR = pj.xyz - ri;
            const float  sj = pj.w;
            const char spin_j = espins[ipj];
            const float qe_j = (spin_j==0) ? 2.0f : 1.0f;
            const float qij = qe_i*qe_j;
            if(1){
                float4 cg = getCoulombGauss( dR, si, sj, qij );
                EeeCoul += (double)cg.x;
            }
            if(1){
                const int spinij = ((int)spin_i)*((int)spin_j);
                const float qq = ((spin_i==0) && (spin_j==0)) ? 2.0f : 1.0f;
                float4 eff_KRS = KRSrho; eff_KRS.w *= qq;
                float4 pg = getPauliGauss_New( dR, si, sj, spinij, eff_KRS );
                EeePaul += (double)pg.x;
                if(bDbg){
                    printf("GPU EE(%d,%d) r=%g si=%g sj=%g spin(%d,%d) q(%g,%g) Ec=%g Ep=%g\n", ie, je, sqrt((double)(dot(dR,dR))), (double)si, (double)sj, (int)spin_i, (int)spin_j, (double)qe_i, (double)qe_j, (double)getCoulombGauss(dR,si,sj,qij).x, (double)pg.x );
                }
            }
        }
    }

    double Etot = Ek + EeeCoul + EeePaul + EaeCoul + EaePaul + Eaa;
    // Match C++ processXYZ_e outputs (first 5): {Etot,Ek,Eee,Eae,Eaa}
    Es[group_id] = (float8)(Etot, Ek, EeeCoul, EaeCoul, Eaa, EeePaul, EaePaul, 0.0f);
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