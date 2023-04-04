
import numpy as np

const_Coulomb =  14.3996448915   # [eV/A]

'''
inline float getLJQ( const Vec3f& dp, Vec3f& f, const Vec3f& REQ, const float R2damp ){
    const float   r2   = dp.norm2();
    // ---- Coulomb
    const float   ir2_ = 1/( r2 + R2damp  );
    float E = COULOMB_CONST*REQ.z*sqrt( ir2_ );
    float F = E*ir2_ ;
    // --- LJ 
    const float  ir2 = 1/r2;
    const float  u2  = REQ.x*REQ.x*ir2;
    const float  u6  = u2*u2*u2;
    const float vdW  = u6*REQ.y;
    E    =      (u6-2.)*vdW     ;
    F    =  12.*(u6-1.)*vdW*ir2 ;
    f.set_mul( dp, -F );
    return E;
}

inline double getMorseQH( const Vec3d& dp, Vec3d& f, const Quat4d& REQH, const double K, double R2damp, const double Cr2_cut, const double Mr2_cut ){
    // Morse Cutoff:
    // E_cut                  = e0*exp( -K * r_cut )
    // ln(E_cut/e0)           = -K * r_cut
    // ln(e0_max/E_cut)/K     = r_cut
    // (ln(e0_max/E_cut)/K)^2 = r2_cut
    //
    // Coulomb Cutoff:
    // E_cut = COULOMB_CONST*Q_max/r_cut
    // r2_cut = (E_cut/(COULOMB_CONST*Q_max))^2
    const double r2    = dp.norm2();
    double E,F;
    // --- Coulomb
    if( r2 > Cr2_cut ){
        const double ir2_  = 1/( r2 + R2damp );
        E = COULOMB_CONST*REQH.z*sqrt( ir2_ );
        F = E*ir2_ ;
    }else{ E=0; F=0; };
    // --- Morse
    if( r2 > Mr2_cut ){
        const double  r  = sqrt( r2   );
        const double  e  = exp( -K*(r-REQH.x) );
        const double  Ae = REQH.y*e;
        const double  He = REQH.w*e; // H-bond correction
        E +=  Ae*(e - 2)   + He;
        F += (Ae*(e - 1)*2 + He)*K/r;
    }
    f.set_mul( dp, F );
    return E;
}

'''

def Coulomb( r, qq=1 ):
    E = const_Coulomb*qq/r
    F = E/r
    return E,F

def LenardJones( r, E0, R0=3.5 ):
    u6   = (R0/r)**6
    e6   =  E0*u6
    E    =      (u6-2.)*e6     
    F    =  12.*(u6-1.)*e6/(r*r)
    return E,F

def Vdw( r, E0, R0=3.5 ):
    u6   =  (R0/r)**6
    E    = -E0*u6     
    F    =  6.*E/(r*r)
    return E,F