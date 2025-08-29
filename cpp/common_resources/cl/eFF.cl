#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

// TODO : this should go elsewhere (physical constants or something)
__constant static const float const_hbar_SI      = 1.054571817e-34;    ///< [J.s]  #6.582119569e-16 # [eV/s]static const float const_Me_SI        = 9.10938356e-31;     ///< [kg]
__constant static const float const_Me_SI        = 9.10938356e-31;     ///< [kg]
__constant static const float const_Matom_SI     = 1.6605402e-27;      ///< [kg]
__constant static const float const_e_SI         = 1.602176620898e-19; ///< [Coulomb]
__constant static const float const_eps0_SI      = 8.854187812813e-12; ///< [F.m = Coulomb/(Volt*m)]
__constant static const float const_eV_SI        = 1.602176620898e-19; ///< [J]
__constant static const float const_Angstroem_SI = 1.0e-10;            ///< [m]
__constant static const float const_K_SI   =  const_hbar_SI*const_hbar_SI/(2*const_Me_SI);   // this is correct see schroedinger equation : https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Preliminaries
__constant static const float const_El_SI  =  const_e_SI*const_e_SI/(4.*M_PI*const_eps0_SI);
__constant static const float const_Ry_SI  = 0.5 * const_El_SI*const_El_SI/const_K_SI;
__constant static const float const_Ry_eV  = 13.6056925944;
__constant static const float const_El_eVA = const_El_SI/( const_eV_SI*const_Angstroem_SI                    );
__constant static const float const_K_eVA  = const_K_SI /( const_eV_SI*const_Angstroem_SI*const_Angstroem_SI );
__constant static const float const_Ke_eVA = const_K_eVA*1.5;
__constant static const float au_Me           = 1822.88973073;
__constant static const float eV_MeAfs        = 17.5882001106;   // [ Me*A^2/fs^2]
__constant static const float const_Coulomb_eVA = 14.3996448915;
__constant static const float const_Bohr_Radius = 0.529177210903;

float2 erfx_e6( float x_, float k ){
    // approximation of erf(k*x)/x and its derivative with maximum error ~ 1e-6
    float x =x_*k;
    if( x>4.5 ){ float y=1/x_; float dy=-y*y; return (float2){ y, dy }; }
    float xx = x*x;
    float even  =  0.9850156202961753  +xx*(-0.02756061032579559  +xx*(-0.00188409579491924  +xx*(-0.003098629936170076 +xx*(-0.001348858853909826  +xx*(-3.98946569988845e-05 ) ) ) ) );
    float odd   = -0.13893350387140332 +xx*(-0.007664292475021448 +xx*( 0.003046826535877866 +xx*( 0.002879338499080343 +xx*( 0.0003260490382458129 +xx*( 1.97093650414204e-06 ) ) ) ) );
    float deven =                           -0.05512122065159118  +xx*(-0.00753638317967696  +xx*(-0.01859177961702045  +xx*(-0.01079087083127861   +xx*(-0.000398946569988845 ) ) ) )  ;
    float dodd  = -0.1389335038714033  +xx*(-0.02299287742506434  +xx*( 0.01523413267938933  +xx*( 0.0201553694935624   +xx*( 0.002934441344212316  +xx*(2.168030154556244e-05 ) ) ) ) );
    float  t = even    + x*odd;
    float dt = deven*x +  dodd;
    float t2 = t*t;
    float t4 = t2*t2;
    float dt8_dx = 8*dt*t*t2*t4;
    float y      = k/(t4*t4 + x);
    float dy     = -y*y*(dt8_dx+1);
    // ToDo : We will need rather (dy/x) for Gauss:Coulomb() => we need to fit it onece more
    return (float2){ y, dy };
}

float exp_p8( float x ){
    if(x<-25) return 0;
    x *= 0.125;
    float xx = x*x;
    float p = (1+x) +
               xx*( 0.5000000000000000   + 0.1666664718006032   *x +
               xx*( 0.04166189077950237  + 0.008304046626191663 *x +
               xx*( 0.001321435070258156 + 0.0001332637951696261*x ) ) );
    p*=p; p*=p; p*=p;
    return p;
}

float gauss_p8( float x ){ return exp_p8( -x*x ); };

// float3 CoulombGauss( float r, float s, float qq ){
//     float Amp = qq*const_El_eVA;
//     float is  = M_SQRT2/s;  // Original from paper (eq.2c)        http://aip.scitation.org/doi/10.1063/1.3272671
//     float2 E_fr  = erfx_e6( r, is );
//     float E      = E_fr.x;
//     float fr     = E_fr.y;
//     float r_s = r*is;
//     float fs  = gauss_p8(r_s) *is*is*is*(0.5*M_2_SQRTPI*Amp);  // How is it possible that "is" was added ?
//     //fs  = is*is*is*(0.5*M_2_SQRTPI*Amp);                 //   1/is^3   because it is multiplied by si and sj later to get (si/(si^2+sj^2)^(3/2) )
//     E *= Amp;
//     fr*= Amp*(1/(r+1e-16));
//     return (float3){ E, fr, fs };
// }


float4 getCoulombGauss( float3 dR, float si, float sj, float qq ){
    float s2   = si*si + sj*sj;
    float s    = sqrt(s2);
    float r    = length( dR );
    //float r    = sqrt( dR.norm2() + 1e-8 );

    // --- copied from CoulombGauss()
    float Amp = qq*const_El_eVA;
    float is  = M_SQRT2/s;  // Original from paper (eq.2c)        http://aip.scitation.org/doi/10.1063/1.3272671
    float2 E_fr  = erfx_e6( r, is );
    float E      = E_fr.x;
    float fr     = E_fr.y;
    float r_s    = r*is;
    float fs     = gauss_p8(r_s) *is*is*is*(0.5*M_2_SQRTPI*Amp);  // How is it possible that "is" was added ?
    E *= Amp;
    fr*= Amp*(1/(r+1e-16));

    return (float4){ E, fr, fs*si, fs*sj };
}

float4 getCoulomb( float3 dp, float qq ){
    float ir2  = 1/( dot(dp,dp) + 1e-32 );
    float ir   = sqrt(ir2);
    float E    = const_Coulomb_eVA*qq*ir;
    float fr   = -E*ir2;
    return (float4){ dp*fr, E };
}

float4 getPauliGauss_New( float3 dR, float si, float sj, int spin, const float4 KRSrho ){
    float r2         = dot(dR,dR) + 1e-8;  // for r=0 there are numercial instabilities

    const float Hartree2eV = 27.211386245988;
    const float A2bohr     = 1/const_Bohr_Radius;

    const float KR=A2bohr*KRSrho.x;
    const float KR2=KR*KR;
    const float KS =A2bohr*KRSrho.y;
    si*=KS; sj*=KS; r2*=KR2;

    float si2        = si*si;
    float sj2        = sj*sj;
    float si2sj2     = si2 + sj2;
    float invsi2sj2  = 1/si2sj2;
    float invsi2sj22 = invsi2sj2*invsi2sj2;
    float invsi2sj23 = invsi2sj2*invsi2sj22;
    float denom_sij  = si*sj*invsi2sj2;
    float si4sj4     = si2*si2 - sj2*sj2;
    float invsj      = 1/sj;
    float invsi      = 1/si;
    float invsj2     = invsj*invsj;
    float invsi2     = invsi*invsi;

    float r2_4   =  4*r2;

    // ------- Kinetic Energy Difference
    float DT      = 1.5*si2sj2*invsi2*invsj2 -      (6*si2sj2 - r2_4)*invsi2sj22;
    float dDT_dsi =  -3*invsi2*invsi         + 4*si*(3*si2sj2 - r2_4)*invsi2sj23;
    float dDT_dsj =  -3*invsj2*invsj         + 4*sj*(3*si2sj2 - r2_4)*invsi2sj23;
    float dDT_dr  =   8*invsi2sj22;      // missing 'r' it is in |dR|

    // ------- Overlap  ..... actually S22 = 2*S**2
    float S22      = 8*denom_sij*denom_sij*denom_sij*exp(-2*r2*invsi2sj2);
    float dS22_dsi = S22*( -3*si4sj4 + r2_4*si2 )*invsi2sj22*invsi;
    float dS22_dsj = S22*( +3*si4sj4 + r2_4*sj2 )*invsi2sj22*invsj;
    float dS22_dr  = -4*S22*invsi2sj2;   // missing 'r' it is in |dR|

    float rho = KRSrho.z;

    float E=0, dE_dDT=0, dE_dS22=0;
    if(spin<=0){
        float invS22m1 = 1/(S22+1);
        E       += - rho*DT*S22  *invS22m1;
        dE_dDT  += -(rho*   S22 )*invS22m1;
        dE_dS22 += -(rho*DT     )*invS22m1*invS22m1;
    }
    if(spin>=0){
        float invS222m1 = 1/( S22*S22-1 );
        E       +=   S22 * DT * ( -rho*S22                     + rho-2 ) *invS222m1;
        dE_dDT  += - S22 *      (  rho*S22                     - rho+2 ) *invS222m1;
        dE_dS22 +=      -  DT * (      S22*(S22*(rho-2)-2*rho) + rho-2 ) *invS222m1*invS222m1;
    }

    float sc = KRSrho.w;

    E         *= Hartree2eV*sc;
    float fsi = (dE_dS22 * dS22_dsi + dE_dDT * dDT_dsi)*Hartree2eV*-KS *sc;
    float fsj = (dE_dS22 * dS22_dsj + dE_dDT * dDT_dsj)*Hartree2eV*-KS *sc;
    float fr  = (dE_dS22 * dS22_dr  + dE_dDT * dDT_dr )*Hartree2eV*KR2 *sc;

    //f = dR * fr;

    //printf( "r %g si %g sj %g DT %g S22 %g E %g anti(%i) \n", sqrt(r2), si,sj, DT,S22, E, anti );
    //return E;
    return (float4){ E, fr, fsi, fsj };
}

float2 addKineticGauss_eFF( float s ){
    // This is currently used in eFF
    float is  = M_SQRT2/s;
    float is2 = is*is*(const_K_eVA*1.5);
    float fs  = is2*is*M_SQRT2;
    float E   = is2;
    //printf( "addKineticGauss s %g is %g is2 %g const_Ke_eVA %20.10f const_K_eVA %20.10f \n", s, is, is2, const_Ke_eVA, const_K_eVA );
    return (float2){ E, fs  };
}

// ==============================================================
//      VERSION 4: High-Performance Kernel with Loop Splitting
#define idDBG    0 // debug thread index
#define nLocal   32
#define nIonMax  8
// ==============================================================

#define bDBGall true


__kernel void localMD(
    __global       int4*        sysinds, // [nsys]   {na,ne,i0p,i0a} size and initial index for each atom
    __global       float4*      pos,     // [ntot]   {x,y,z,w} positions (and size) of ions and electrons
    __global       float4*      vel,     // [ntot]   {vx,vy,vz,dw/dt} velocities of ions and electrons (including change of size)
    __global const float8*      aParams, // [ntot_a] parameters of ions { Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD }
    __global const signed char* espins,  // [ntot]   {spin}
    __global       float4*      fout,    // [ntot]   {fx,fy,fz,fw} output force buffer
    const int    nsys,
    const int    nsteps,
    const float  dt,
    const float  damping,
    const float4 KRSrho,
    const int    bFrozenCore
) {
    const int isys  = get_group_id(0);
    const int4 inds = sysinds[isys];
    const int lid   = get_local_id(0);
    const int ip    = inds.z + lid;
    const int na    = inds.x;
    const int ne    = inds.y;
    const int ntot  = inds.x + inds.y;

    if (get_global_id(0) == idDBG){
        int nSys = get_num_groups(0);
        int nL   = get_local_size(0);
        printf("GPU: nSys=%d nL=%d na=%d ne=%d ntot=%d   nsteps=%d dt=%f damping=%f  | KRSrho(%.6f,%.6f,%.6f,%.6f) bFrozenCore=%d\n", nSys, nL, na, ne, ntot, nsteps, dt, damping, KRSrho.x, KRSrho.y, KRSrho.z, KRSrho.w, bFrozenCore);
        for(int isys=0; isys<nSys; isys++){
            int4 is = sysinds[isys];
            int nt = is.x + is.y;
            printf("System %2i na %2i ne %2i ntot %2i  | i0 %4i  i0a %4i \n", isys, is.x, is.y, nt, is.z, is.w );
            for(int i=0; i<nt; i++){
                float4 pi = pos[i];
                printf("sys %2i p# %2i p(%12.8f,%12.8f,%12.8f,%12.8f)", isys, i, pi.x, pi.y, pi.z, pi.w );
                if(i<na){
                    float8 api = aParams[i];
                    printf(" Atom params(%12.8f,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f,%12.8f)\n", api.s0, api.s1, api.s2, api.s3, api.s4, api.s5, api.s6, api.s7 );
                }else{
                    int sp = espins[ is.z + i ];
                    printf(" Electron sz %12.8f spin %d\n", pi.w, sp);
                }
            }
        }
    }

    if (lid >= ntot) return;

    // --- Local Memory State ---
    __local float4      l_pos     [nLocal];
    __local float8      l_aparams [nIonMax];
    __local signed char l_spins   [nLocal];

    // --- Load initial state from Global to Local memory ---
    float4 posi = pos[ip];
    float4 veli = vel[ip];

    l_pos[lid] = posi;
    if (lid < inds.x){ l_aparams[lid       ] = aParams[inds.w+lid];   }
    else             { l_spins  [lid-inds.x] = espins[ip]; }

    barrier(CLK_LOCAL_MEM_FENCE);

    // --- VALIDATION PATH: serial pairwise accumulation matching CPU (Newton's 3rd law)
    // Run only on lid==0. All other work-items return now to avoid executing the parallel path.
    if (lid != 0) { return; }
    {
        float4 F[nLocal];
        for (int i=0; i<ntot; ++i) { F[i] = (float4)(0.0f); }

        // Kinetic size forces for electrons
        for (int i=na; i<ntot; ++i) {
            float2 fk = addKineticGauss_eFF(l_pos[i].w);
            F[i].w += fk.y;
        }

        // Ion-Ion (AA)
        for (int i=0; i<na; ++i) {
            const float8 pari = l_aparams[i];
            const float  Qi   = pari.s0 - pari.s2;
            for (int j=i+1; j<na; ++j) {
                const float8 parj = l_aparams[j];
                const float  Qj   = parj.s0 - parj.s2;
                const float3 dR   = l_pos[j].xyz - l_pos[i].xyz;
                float4 c = getCoulomb(dR, Qi*Qj);
                float3 f = c.xyz;
                F[i].xyz += f;
                F[j].xyz -= f;
            }
        }

        // Ion-Electron (AE)
        for (int i=0; i<na; ++i) {
            const float8 pari = l_aparams[i];
            const float  Qi   = pari.s0 - pari.s2;
            const float  Ri   = pari.s1;
            for (int j=na; j<ntot; ++j) {
                const float4 ej = l_pos[j];
                const float3 dR = ej.xyz - l_pos[i].xyz;
                float4 cg = getCoulombGauss(dR, Ri, ej.w, -Qi);
                float3 f = dR * cg.y;
                if(bDBGall){ printf("GPU[serial] AE(%i,%i) dR(%.3f,%.3f,%.3f) s(%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f) | %.3f,%.3f\n", i, j-na, dR.x, dR.y, dR.z, Ri, ej.w, f.x, f.y, f.z, 0.0, cg.w); }
                F[i].xyz += f;
                F[j].xyz -= f;
                F[i].w   += cg.z;   // size force on ion
                F[j].w   += cg.w;   // size force on electron
                if (bFrozenCore) {
                    // Pauli with frozen-core (if used) would go here
                }
            }
        }

        // Electron-Electron (EE): Coulomb + Pauli
        for (int i=na; i<ntot; ++i) {
            const float4 ei   = l_pos[i];
            const int    si   = l_spins[i-na];
            for (int j=i+1; j<ntot; ++j) {
                const float4 ej = l_pos[j];
                const int    sj = l_spins[j-na];
                const float3 dR = ej.xyz - ei.xyz;
                float4 cg = getCoulombGauss  (dR, ei.w, ej.w, 1.0f);
                // CPU scales Pauli by qq = 2 if both spins are 0 (paired), else 1.0
                float  qq = ((si==0) && (sj==0)) ? 2.0f : 1.0f;
                float4 KRS = KRSrho; KRS.w *= qq;
                float4 pg = getPauliGauss_New(dR, ei.w, ej.w, si*sj, KRS);
                float  fr = cg.y + pg.y;
                float3 f  = dR * fr;
                if(bDBGall){ printf("GPU[serial] EE(%i,%i) dR(%.3f,%.3f,%.3f) s(%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f) | %.3f,%.3f | Paul:(%.3f,%.3f,%.3f) | %.3f,%.3f\n",
                    i-na, j-na,
                    dR.x, dR.y, dR.z, ei.w, ej.w,
                    dR.x*cg.y, dR.y*cg.y, dR.z*cg.y, cg.z, cg.w,
                    dR.x*pg.y, dR.y*pg.y, dR.z*pg.y, pg.z, pg.w ); }
                F[i].xyz += f;
                F[j].xyz -= f;
                F[i].w   += cg.z + pg.z;   // size force on i
                F[j].w   += cg.w + pg.w;   // size force on j
            }
        }
        for (int i=0; i<ntot; ++i) { fout[inds.z + i] = F[i]; }
        return;
    }

    for (int i_step = 0; i_step < nsteps; ++i_step) {

        float4 forcei = (float4)(0.0f);
        const bool i_am_ion = (lid < na);

        // --- PART 1: FORCE CALCULATION ---

        if (i_am_ion) {  // =============   I am an ION thread
            const float8 pari    = l_aparams[lid];
            const float  Qi      = pari.s0 - pari.s2;
            const float  Ri      = pari.s1;

            //if (get_global_id(0) == idDBG){   printf("Ion %3i Qi %12.8e Ri %12.8e pos(%12.8e,%12.8e,%12.8e,%12.8e)\n", lid, Qi, Ri, posi.x,posi.y,posi.z,posi.w );}
            if(bDBGall){ printf("Ion %3i Qi %12.8e Ri %12.8e pos(%12.8e,%12.8e,%12.8e,%12.8e)\n", lid, Qi, Ri, posi.x,posi.y,posi.z,posi.w );}
            // --- Ion-Ion Interactions ---
            for (int j = 0; j < na; ++j) {
                if (lid == j) continue;
                const float3 dR            = l_pos[j].xyz - posi.xyz;
                const float8 parj          = l_aparams[j];
                const float  Qj            = parj.s0 - parj.s2;
                //const float  Rj            = parj.s1;
                float4 fij = getCoulomb(dR, Qi * Qj);
                //if (get_global_id(0) == idDBG){ printf("GPU AA(%i,%i) dR(%.3f,%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f)\n", lid, j, dR.x, dR.y, dR.z, fij.x,fij.y,fij.z); }
                if(bDBGall){ printf("GPU AA(%i,%i) dR(%.3f,%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f)\n", lid, j, dR.x, dR.y, dR.z, fij.x,fij.y,fij.z); }
                forcei.xyz += fij.xyz;
            }
            // --- Ion-Electron Interactions ---
            for (int j = na; j < ntot; ++j) {
                const float4 pj = l_pos[j];
                const float3 dR = pj.xyz - posi.xyz;
                float4 fg = getCoulombGauss (dR, Ri, pj.w, Qi );
                //if (get_global_id(0) == idDBG){ printf("GPU AE(%i,%i) dR(%.3f,%.3f,%.3f) s(%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f) | %.3f,%.3f\n", lid, j-na, dR.x, dR.y, dR.z, Ri, pj.w, dR.x*fg.y, dR.y*fg.y, dR.z*fg.y, 0.0, fg.w); }
                if(bDBGall){ printf("GPU AE(%i,%i) dR(%.3f,%.3f,%.3f) s(%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f) | %.3f,%.3f\n", lid, j-na, dR.x, dR.y, dR.z, Ri, pj.w, dR.x*fg.y, dR.y*fg.y, dR.z*fg.y, 0.0, fg.w); }
                forcei.xyz += dR * fg.y;
                forcei.w   += fg.w;
                if( bFrozenCore ){
                    //float4 fp = getPauliGauss_New(dR, Ri, pj.w, 0, KRSrho);
                }
            }
        } else {  // =========== I am an ELECTRON thread
            const int spini = l_spins[lid - na];

            // --- Kinetic Force ---
            float2 fk = addKineticGauss_eFF(posi.w);
            forcei.w += fk.y;

            // --- Electron-Ion Interactions ---
            for (int j = 0; j < na; ++j) {
                const float8 parj = l_aparams[j];
                const float3 dR   = l_pos[j].xyz - posi.xyz;
                const float Qj    = parj.s0 - parj.s2;
                const float Rj    = parj.s1;
                float4 fg = getCoulombGauss  (dR, posi.w, Rj, -1.0f * Qj);
                //if (get_global_id(0) == idDBG){ printf("GPU AE(%i,%i) dR(%.3f,%.3f,%.3f) s(%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f) | %.3f,%.3f\n", j, lid-na, dR.x, dR.y, dR.z, posi.w, Rj, dR.x*fg.y, dR.y*fg.y, dR.z*fg.y, fg.z, 0.0); }
                if(bDBGall){ printf("GPU AE(%i,%i) dR(%.3f,%.3f,%.3f) s(%.3f,%.3f) -> Coul:(%.3f,%.3f,%.3f) | %.3f,%.3f\n", j, lid-na, dR.x, dR.y, dR.z, posi.w, Rj, dR.x*fg.y, dR.y*fg.y, dR.z*fg.y, fg.z, 0.0); }
                forcei.xyz -= dR * fg.y;
                forcei.w   += fg.z;
                if( bFrozenCore ){
                    //float4 fp = getPauliGauss_New(dR, posi.w, Rj, 0, KRSrho);
                }
            }

            // --- Electron-Electron Interactions ---
            for (int j = na; j < ntot; ++j) {
                if (lid == j) continue;
                const float4 pj = l_pos[j];
                const float3 dR           = pj.xyz - posi.xyz;
                const signed char spinj   = l_spins[j - na];
                float4 fg = getCoulombGauss  (dR, posi.w, pj.w, 1.0f);
                float4 fp = getPauliGauss_New(dR, posi.w, pj.w, spinj * spini, KRSrho);
                //if (get_global_id(0) == idDBG){
                if(bDBGall){
                    printf("GPU EE(%i,%i) dR(%g,%g,%g) s(%g,%g)   Coul: (%g,%g,%g) | %g,%g  | Paul: (%g,%g,%g) | %g,%g\n",    lid-na, j-na, dR.x, dR.y, dR.z, posi.w, pj.w, dR.x*fg.y, dR.y*fg.y, dR.z*fg.y, fg.z, dR.x*fp.y, dR.y*fp.y, dR.z*fp.y, fp.z); 
                }
                forcei.xyz += dR * (fg.y + fp.y);
                forcei.w   += fg.z + fp.z; // f_si
                // TODO: How to handle fsj ? The C++ version adds it to the other particle.
            }
        }

        barrier(CLK_LOCAL_MEM_FENCE);

        // --- PART 2: INTEGRATION STEP ---

        //if (get_global_id(0) == idDBG)
        { printf("update iter %4i il %4i pos(%12.8e,%12.8e,%12.8e,%12.8e) vel(%12.8e,%12.8e,%12.8e,%12.8e) force(%12.8e,%12.8e,%12.8e,%12.8e)\n", i_step, lid, posi.x,posi.y,posi.z,posi.w, veli.x,veli.y,veli.z,veli.w, forcei.x,forcei.y,forcei.z,forcei.w );}

        fout[ip] = forcei;

        veli *= damping;
        veli += forcei * dt;
        posi += veli   * dt;
        if (!i_am_ion) { posi.w = fmax(posi.w, 0.001f); } // Electron size cannot be zero or negative

        // Update local memory for the next step
        //l_pos[lid] = pos;
        //l_vel[lid] = vel;

        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- Write final state back to Global Memory ---
    pos[ip] = posi;
    vel[ip] = veli;
}




// =======================================
// ========   Eval Ions           ========
// =======================================

__kernel void eval_ions(
    const int na, const int ne,
    __global const float4* apos,
    __global const float4* apar,
    __global const float4* epos,
    __global const int   * espin,
    __global       float4* aforce,
    const float4 KRSrho)
{
    const int gid = get_global_id(0);
    if (gid >= na) return;

    const float4 ai = apos[gid];
    const float4 pi = apar[gid];
    float4 fi = (float4)(0.0f);

    __local float4 LAPOS[32];
    __local float4 LAPAR[32];

    // ion–ion interactions
    for (int base=0; base<na; base+=32) {
        int lid = get_local_id(0);
        int j = base + lid;
        if (j < na) { LAPOS[lid] = apos[j]; LAPAR[lid] = apar[j]; }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int l=0; l<32 && base+l<na; ++l) {
            if (base+l == gid) continue;
            float4 aj = LAPOS[l];
            float4 pj = LAPAR[l];
            float3 dp = aj.xyz - ai.xyz;
            float4 c  = getCoulomb(dp, pj.w * pi.w);
            fi.xyz += c.xyz;
            fi.w   += c.w;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // ion–electron interactions
    __local float4 LEPOS[32];
    __local int    LSPIN[32];
    for (int base=0; base<ne; base+=32) {
        int lid = get_local_id(0);
        int j = base + lid;
        if (j < ne) { LEPOS[lid] = epos[j]; LSPIN[lid] = espin[j]; }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int l=0; l<32 && base+l<ne; ++l) {
            float4 ej = LEPOS[l];
            int sj = LSPIN[l];
            float3 dp = ej.xyz - ai.xyz;
            float4 c  = getCoulombGauss(dp, ej.w, pi.w, 1.0f);
            float4 p  = getPauliGauss_New(dp, ej.w, pi.w, sj * 0, KRSrho);
            fi.xyz += c.xyz;
            fi.w   += c.w;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    aforce[gid] = fi;
}

// =======================================
// ========   Eval Electrons      ========
// =======================================

__kernel void eval_electrons(
        const int na, const int ne,
        __global const float4* apos,
        __global const float4* apar,
        __global const float4* epos,
        __global const int   * espin,
        __global       float4* eforce,
        const float4 KRSrho)
{
    const int gid = get_global_id(0);          // electron index
    if (gid >= ne) return;

    const float4 ei   = epos [gid];
    const int    si   = espin[gid];
    float4 fi        = (float4)(0.0f);

    __local float4 LEPOS[32];
    __local int    LSPIN[32];

    // -------- electron–electron ----------
    for (int base=0; base<ne; base+=32) {
        const int lid = get_local_id(0);
        const int j   = base + lid;
        if (j < ne) { LEPOS[lid] = epos[j];  LSPIN[lid]=espin[j]; }
        barrier(CLK_LOCAL_MEM_FENCE);

        for (int l=0; l<32 && base+l<ne; ++l) {
            if (base+l == gid) continue;
            float4 ej = LEPOS[l];
            int    sj = LSPIN[l];
            float3 dp = ej.xyz - ei.xyz;
            // coulomb / pauli
            float4 c  = getCoulombGauss(dp, ei.w, ej.w, 1.0f);
            float  qq = ((si==0) && (sj==0)) ? 2.0f : 1.0f; // CPU scales Pauli by qq when both spins are 0
            float4 KRS = KRSrho; KRS.w *= qq;
            float4 p  = getPauliGauss_New(dp, ei.w, ej.w, si*sj, KRS);
            fi.xyz += dp*(c.y + p.y);
            fi.w   += c.z + p.z;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // -------- electron–ion ----------
    __local float4 LAPOS[32];
    __local float4 LAPAR[32];
    for (int base=0; base<na; base+=32) {
        const int lid = get_local_id(0);
        const int j   = base + lid;
        if (j < na) { LAPOS[lid]=apos[j];  LAPAR[lid]=apar[j]; }
        barrier(CLK_LOCAL_MEM_FENCE);

        for (int l=0; l<32 && base+l<na; ++l) {
            float4 aj = LAPOS[l];
            float4 pj = LAPAR[l];
            float3 dp = aj.xyz - ei.xyz;
            float4 c  = getCoulombGauss(dp, ei.w, pj.y, pj.x);
            fi.xyz += dp*(c.y+c.y);
            fi.w   += c.z;
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    eforce[gid] = fi;
}
