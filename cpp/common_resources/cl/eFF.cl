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

// =======================================
// ========   Eval Electrons      ========
// =======================================

__kernel void eval_electrons(
    int na,                    // 1 
    int ne,                    // 2
    __global float4*  apos,    // 3      
    __global float4*  aforce,  // 4  
    __global float4*  aParams, // 5 
    __global float4*  epos,    // 6     
    __global float4*  eforce,  // 7
    __global int*     espin,   // 8
    float4 KRSrho              // 9
){
    __local float4 LATOM[32];   // local buffer for atom positions
    __local int    LSPIN[32];   // local buffer for atom parameters
    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory
    const int i0e = iS*ne;    // index of first atom in atoms array
    const int ie  = iG + i0e; // index of atom in atoms array
    
    const float4 ai    = epos [ie];
    const int    spini = espin[ie];
    float4 fi          = (float4){0.0f,0.0f,0.0f,0.0f};
    float  Ei          = 0.0f; 

    if(iG>=ne) return;

    // ========= Kinetic energy
    float2 efk = addKineticGauss_eFF( ai.w );
    fi.w += efk.y;
    Ei   += efk.x;

    // ========= Electron-Electron  ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){      
        const int i=j0+iL;             
        if(i<ne){                       
            LATOM[iL] = epos [i+i0e]; 
            LSPIN[iL] = espin[i+i0e]; 
        }
        barrier(CLK_LOCAL_MEM_FENCE);  
        for (int jl=0; jl<nL; jl++){    
            const int ja=j0+jl;        
            if( (ja!=iG) && (ja<ne) ){  
                const float4 aj     = LATOM[jl];    
                const int    spinj  = LSPIN[jl];   
                const float3 dp     = aj.xyz - ai.xyz; // vector between atoms
                const float4 coul   = getCoulombGauss  ( dp, ai.w, aj.w, 1.0 );
                const float4 paul   = getPauliGauss_New( dp, ai.w, aj.w, spini*spinj, KRSrho );
                Ei      +=    (coul.x + paul.x);
                fi.xyz  += dp*(coul.y+coul.y);
                fi.w    +=     coul.z+paul.z;
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }

    // ========= Electron-Core  ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    for (int j0=0; j0<nG; j0+=nL){      
        const int i=j0+iL;             
        if(i<ne){                       
            LATOM[iL] = epos   [i+i0e]; 
            LAPAR[iL] = aParams[i+i0e]; 
        }
        barrier(CLK_LOCAL_MEM_FENCE);   
        for (int jl=0; jl<nL; jl++){   
            const int ja=j0+jl;         
            if( (ja!=iG) && (ja<ne) ){  
                const float4 aj  = LATOM[jl];    
                const float4 Pj  = LAPAR[jl];    // { x=Q,y=sQ,z=sP,w=cP }
                const float3 dp  = aj.xyz - ei.xyz; 

                const float4 coul  = getCoulombGauss  ( dp, ei.w, Pj.y, Pj.x );
                //const float4 paul   = getPauliGauss_New( dp, ai.w, aj.w, spini*spinj, KRSrho );
                const float4 paul  = getPauliGauss_New( dp, ei.w, Pj.z, 0, KRSrho );    

                Ei      +=    (coul.x + paul.x);
                fi.xyz  += dp*(coul.y+coul.y);
                fi.w    +=     coul.z+paul.z;
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    eforce[ie] += fi;

    eforce[ie] += fi;
}

// =======================================
// ========   Eval Ions           ========
// =======================================

__kernel void eval_ions(
    int na,                    // 1
    int ne,                    // 2
    __global float4*  apos,    // 3       
    __global float4*  aforce,  // 4    
    __global float4*  aParams, // 5  
    __global float4*  epos,    // 6      
    __global float4*  eforce,  // 7 
     float4 KRSrho            // 8
){
    __local float4 LATOM[32];   // local buffer for atom positions
    __local float4 LAPAR[32];   // local buffer for atom parameters
    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory
    const int i0e = iS*ne;    // index of first atom in atoms array
    const int ie  = iG + i0e; // index of atom in atoms array
    
    const float4 ei    = epos   [ie];
    float4 fi          = (float4){0.0f,0.0f,0.0f,0.0f};
    float  Ei          = 0.0f; 

    if(iG>=ne) return;

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){      
        const int i=j0+iL;              
        if(i<na){                      
            LATOM[iL] = apos [i+i0a];  
        }
        barrier(CLK_LOCAL_MEM_FENCE);   
        for (int jl=0; jl<nL; jl++){   
            const int ja=j0+jl;        
            if( (ja!=iG) && (ja<na) ){  
                const float4 aj     = LATOM[jl];    
                const float3 dp     = aj.xyz - ai.xyz; 
                fe += getCoulomb( dp, aj.w*ai.w );
                { // ToDo - core electron interaction

                }
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){      
        const int i=j0+iL;             
        if(i<na){                       
            LATOM[iL] = apos [i+i0a];  
        }
        barrier(CLK_LOCAL_MEM_FENCE);  
        for (int jl=0; jl<nL; jl++){    
            const int ja=j0+jl;        
            if( (ja!=iG) && (ja<na) ){   
                const float4 aj     = LATOM[jl];    
                const float3 dp     = aj.xyz - ai.xyz; 
                fe += getCoulomb( dp, aj.w*ai.w );
                { // ToDo - core electron interaction

                }
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }



}


/*
__kernel void eFF_evalAE(
    int na, int ne,
    __global float4*  apos,        
    __global float4*  aforce,      
    __global float4*  epos,         
    __global float4*  eforce,       
    __global int   *  espin         
){

}

*/

