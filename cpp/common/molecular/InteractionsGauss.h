
#ifndef InteractionsGauss_h
#define InteractionsGauss_h
/// @file InteractionsGauss.h @brief Implements functions for evaluation of interaction between gaussian functions such as product, overlap, kinetic energy and electrostatic interaction
/// @ingroup Electron_Forcefield

#include "math.h"
#include "fastmath.h"

/*!  
eFF : Electron Force Field
---------------------------

[1] http://aip.scitation.org/doi/10.1063/1.3272671
The dynamics of highly excited electronic systems: Applications of the electron force field
Julius T. Su, William A. Goddard

[2] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
Non-adiabatic dynamics modeling framework for materials in extreme conditions
Hai Xiao, Andrés Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard,

### NOTES:

1) It seems that decrease of kinetic energy by sharing electron between two atoms is dominant contribution which makes formation of bonds fabourable
    * H2 Molecule perhaps cannot be stable without this contribution ( i.e. with fixed radius of electron blobs )

### Params
 \verbatim
            a             b             c           d           e           s-core
Al     0.486000       1.049000      0.207000                              1.660000
Si     0.320852       2.283269      0.814857                              1.691398
C     22.721015       0.728733      1.103199     17.695345    6.693621    0.621427
O     25.080199       0.331574      1.276183     12.910142    3.189333    0.167813

Erf approximation:
 Gaussian:    F = (x2-1)**2 / sqrtPi
 Erf          E = x*(1 + x2 * ( -0.66666666666 + 0.2*x2 ) ) * (2/(16.0/15.0))
\endverbatim
*/

/*
const double const_hbar_SI      = 1.054571817e-34;    ///< [J.s]  #6.582119569e-16 # [eV/s]
const double const_Me_SI        = 9.10938356e-31;     ///< [kg]
const double const_e_SI         = 1.602176620898e-19; ///< [Coulomb]
const double const_eps0_SI      = 8.854187812813e-12; ///< [F.m = Coulomb/(Volt*m)]
const double const_eV_SI        = 1.602176620898e-19; ///< [J]
const double const_Angstroem_SI = 1.0e-10;

const double const_K_SI     =  const_hbar_SI*const_hbar_SI/const_Me_SI;
const double const_El_SI    =  const_e_SI*const_e_SI/(4.*M_PI*const_eps0_SI);
const double const_Ry_SI    = 0.5 * const_El_SI*const_El_SI/const_K_SI;

const double const_Ry_eV  = 13.6056925944;
const double const_El_eVA = const_El_SI/( const_e_SI*const_Angstroem_SI );
const double const_K_eVA  = (const_El_eVA*const_El_eVA)/(2*const_Ry_eV);
const double const_Ke_eVA = const_K_eVA*1.5;
*/

#include "physics_constants.h"


// ================================================================
// ========= Kinetic
// ================================================================

inline double addKineticGauss( double s, double& fs ){
    // WARRNING: NOT SURE THIS IS USED NOW
    double is  = 1/s;
    double is2 = is*is*(const_K_eVA*1.5);
    fs += 2.*is2*is;
    //printf( "addKineticGauss s %g is %g is2 %g const_Ke_eVA %20.10f const_K_eVA %20.10f \n", s, is, is2, const_Ke_eVA, const_K_eVA );
    return is2;
}

inline double addKineticGauss_eFF( double s, double& fs ){
    // This is currently used in eFF
    double is  = M_SQRT2/s;
    double is2 = is*is*(const_K_eVA*1.5);
    fs += is2*is*M_SQRT2;
    //printf( "addKineticGauss s %g is %g is2 %g const_Ke_eVA %20.10f const_K_eVA %20.10f \n", s, is, is2, const_Ke_eVA, const_K_eVA );
    return is2;
}

// ================================================================
// ========= Coulomb
// ================================================================

/*
inline double Coulomb( double r, double s, double& fr, double& fs ){
    constexpr const double const_F2 = M_2_SQRTPI * M_SQRT2;
    double ir   = 1./r;
    double is   = 1./s;
    double r_s  = r*is;
    double r_2s = M_SQRT1_2 * r_s; // This is for charge-density blobs (assuming si,sj comes from charge denisty)
    //double r_2s = r_s;
    //double r_2s = M_SQRT2   * r_s; // This is for wavefunction blobs (assuming si,sj comes from wavefunction)
    double e1   = ir * const_El_eVA;
    double e2   = erf(  r_2s      );            // ToDo : this should be possible to compute together !!!
    double g    = exp( -r_2s*r_2s ) * const_F2;
    double f1   = -e1*ir;
    double f2   = g*is*0.5;
    double e1f2 = e1*f2;
    fr          = (f1*e2 + e1f2)*ir      ;
    //printf( "r %g fr %g = (f1 %g * e2 %g )+(e1 %g *f2 %g) r_2s %g r %g s %g\n", r, fr, f1, e2, e1, f2, r_2s, r, s );
    fs          =          e1f2 *r_s * is;
    return e1 * e2;
}
*/
inline double CoulombGauss( double r, double s, double& fr, double& fs, double qq ){
    // ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    /*
    //constexpr const double const_F2 = -2.*sqrt(2./M_PI);
    constexpr const double const_F2 = M_2_SQRTPI * M_SQRT2;
    //double ir   = 1/r; //(r+1.e-8);
    //double is   = 1/s; //(s+1.e-8);
    double ir   = 1/(r+1.e-8);
    double is   = 1/(s+1.e-8);
    double r_s  = r*is;
    //double r_2s = M_SQRT2 * r_s; // Original from paper (eq.2c)        http://aip.scitation.org/doi/10.1063/1.3272671
    double r_2s = r_s;             // modified according to derivation   following eq.11 http://doi.wiley.com/10.1002/wcms.78
    double e1   = qq*ir * const_El_eVA;
    double e2   = erf(  r_2s      );
    double g    = exp( -r_2s*r_2s ) * const_F2;
    double f1   = -e1*ir;
    double f2   = g*is;
    double e1f2 = e1*f2;
    fr = (f1*e2 + e1f2)*ir;
    fs =          e1f2 *r_s * is;
    */

    double Amp = qq*const_El_eVA;
    //double Amp = 1;

    //double is  = M_SQRT1_2/s;
    double is  = M_SQRT2/s;  // Original from paper (eq.2c)        http://aip.scitation.org/doi/10.1063/1.3272671
    double E   = erfx_e6( r, is, fr ); // This is for charge-density blobs (assuming si,sj comes from charge denisty)
    //double E   = erf( r*is )/r;
    double r_s = r*is;
    fs  = gauss_p8(r_s) *is*is*is*(0.5*M_2_SQRTPI*Amp);  // How is it possible that "is" was added ?
    //fs  = is*is*is*(0.5*M_2_SQRTPI*Amp);                 //   1/is^3   because it is multiplied by si and sj later to get (si/(si^2+sj^2)^(3/2) )
    E *= Amp;
    //fr*=-Amp;
    fr*= Amp*(1/(r+1e-16));

    //printf( "CoulombGauss(r %g,s %g)-> E %g fr %g fs %g qq %g \n", r, s, E, fr, fs, qq );
    //exit(0);
    return E;
}

/*
inline double CoulombGauss_FixSize( double r, double s, double& fr, double qq ){
    // ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    //constexpr const double const_F2 = -2.*sqrt(2./M_PI);
    constexpr const double const_F2 = M_2_SQRTPI * M_SQRT2;
    double ir   = 1./r; //(r+1.e-8);
    double is   = 1./s; //(s+1.e-8);
    double r_s  = r*is;
    double r_2s = M_SQRT2 * r_s;
    double e1   = qq*ir * const_El_eVA;
    double e2   = erf(  r_2s );
    double g    = exp( -r_2s*r_2s ) * const_F2;
    double f1   = -e1*ir;
    double f2   = g*is;
    fr = (f1*e2 + e1*f2)*ir;
    return e1 * e2;
}
*/

// version of CoulombGauss optimized for fixed size of electron clouds
inline double CoulombGauss_FixSize( double r, double beta, double& fr, double qq ){
    // ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    //constexpr const double const_F2 = -2.*sqrt(2./M_PI);
    constexpr const double const_F2 = M_2_SQRTPI * M_SQRT2;
    double ir = 1./r;
    double r_b  = beta * r;
    double e1   = qq*ir * const_El_eVA;
    double e2   = erf(  r_b  );
    double f1   = -e1*ir;
    double f2   = exp( -r_b*r_b ) * const_F2 * beta;
    fr = (f1*e2 + e1*f2)*ir;
    return e1 * e2;
}

inline double addCoulombGauss( const Vec3d& dR, double s, Vec3d& f, double& fsi, double qq ){
    double r    = dR.norm();
    //double r    = sqrt( dR.norm2() + 1e-8 );
    double fr,fs;
    double E = CoulombGauss( r, s, fr, fs, qq );
    //printf( "addCoulombGauss E %g s %g r %g fr %g fs %g \n", E, s, r, fs, fr );
    fsi += fs*s;
    f.add_mul( dR, fr );
    return E;
}

inline double addCoulombGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, double qq ){
    double s2   = si*si + sj*sj;
    double s    = sqrt(s2);
    double r    = dR.norm();
    //double r    = sqrt( dR.norm2() + 1e-8 );
    double fs,fr;
    double E = CoulombGauss( r, s, fr, fs, qq );
    //printf( "addCoulombGauss: fs %g s[i,j](%g,%g) fs[i,j](%g,%g) \n", fs, si,sj, fs*si, fs*sj );
    fsi += fs*si;
    fsj += fs*sj;
    f.add_mul( dR, fr );
    return E;

}

// ================================================================
// =========   Pauli Repulsion - New - Dens Overlap
// ================================================================

//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double DensOverlapGauss_S( double r2, double amp, double si, double sj, double& dSr, double& dSsi, double& dSsj,
    double si2, double sj2, double is2, double is4
){
    // eq. 12 in (Xiao, H., et. al. Mechanics of Materials, 90, 243–252 (2015). https://doi.org/10.1016/j.mechmat.2015.02.008 )
    // E = (2/(si/sj+si/sj))^3 * exp( -2*r^2/(si^2+sj^2) )
    double a    = 2.*(si*sj)*is2;
    double a2   = a*a;
    double Aa2  = amp*a2;
    double e1   = Aa2*a;              // (2/(si/sj+si/sj))^3
    double e2   = exp( -2*r2*is2 );   // exp( -2*r^2/(si^2+sj^2) )

    double E    = e1*e2;
    double f1e2 = 6.*e2*Aa2*(si2-sj2)*is4;
    double e1f2 = 4.*E     * r2      *is4;

    //if(i_DEBUG>0) printf( "E %g e1 %g e2 %g(%g) r %g s%g(%g,%g) \n", E, e1, e2,-2*r2*is2, sqrt(r2), is2, si, sj );

    dSsi = -(e1f2*si - f1e2*sj);
    dSsj = -(e1f2*sj + f1e2*si);
    dSr  = E   *(-4.*is2);

    return E;
}

inline double DensOverlapGauss_Snorm( double r2, double amp, double si, double sj, double& dSr, double& dSsi, double& dSsj,
    double si2, double sj2, double is2, double is4
){
    /// eq. 12 in (Xiao, H., et. al. Mechanics of Materials, 90, 243–252 (2015). https://doi.org/10.1016/j.mechmat.2015.02.008 )
    /// E = (2/(si/sj+si/sj))^3 * exp( -2*r^2/(si^2+sj^2) )

    //amp*= const_K_eVA;

    double a    = 2.*(si*sj)*is2;
    double a2   = a*a;
    double e1   = a2*a;              // (2/(si/sj+si/sj))^3
    double e2   = exp( -2*r2*is2 );  // exp( -2*r^2/(si^2+sj^2) )

    /// prefactr derived from T :  e0    = (     si^4 +   sj^4 - 1.25*(si^2*sj^2) )/((si^2*sj^2)*(si^2+sj^2))
    /// and its derivative        de0/da = ( 4.5*si^4 - 2*sj^4 - 4.0 *(si^2*sj^2) )/((si^3     )*(si^2+sj^2)**2)
    double sisj   = si*sj;
    double sisj2  = sisj*sisj;
    //double isisj2 = 1/sisj2;
    double si4    = si2*si2;
    double sj4    = sj2*sj2;
    double pre    = amp * 3.3; // What is this 3.3 constant ?
    double e0     = pre * (     si4 +  sj4 - 1.25*sisj2 )*is2/sisj2;
    //double e0     = amp;
    double e0si   = pre * ( 4.5*si4 -2*sj4 - 4*sisj2 )*is4/(si2*si);
    double e0sj   = pre * ( 4.5*sj4 -2*si4 - 4*sisj2 )*is4/(sj2*sj);

    double e1e2 = e1*e2;
    double E    = e1e2*e0;
    //if(i_DEBUG>0) printf( "E %g e1 %g e2 %g(%g)  e0 %g r %g s%g(%g,%g) \n", E, e1, e2,-2*r2*is2, e0, sqrt(r2), is2, si, sj );
    double f1e2 = 6.*e2*a2*(si2-sj2)*is4;
    double e1f2 = 4.*e1e2 * r2      *is4;

    dSsi = -((e1f2*si - f1e2*sj)*e0 + e1e2*e0si);
    dSsj = -((e1f2*sj + f1e2*si)*e0 + e1e2*e0sj);
    dSr  = E   *(-4.*is2);

    return E;
}

//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double DensOverlapGauss_P( double r2, double amp, double si, double sj, double& dSr, double& dSsi, double& dSsj,
    double si2, double sj2, double is2, double is4
){
    /// eq. 12 in (Xiao, H., et. al. Mechanics of Materials, 90, 243–252 (2015). https://doi.org/10.1016/j.mechmat.2015.02.008 )
    /// E = (2/(si/sj+si/sj))^5 * (r12-s2/sqrt2)^2 * exp( -2*(r^2-sj/sqrt2)/(si^2+sj^2) )
    double isi  = 1/si;
    double a    = 2.*(si*sj)*is2;
    double a2   = a *a;
    double a4   = a2*a2;
    double e1   = amp*a4*a;   //  (2/(si/sj+si/sj))^5
    double e2   = exp( -2*r2*is2 );       //  exp( -2*(r^2-sj/sqrt2)/(si^2+sj^2) )
    double e3   = r2*isi*isi;             //  (r12-s2/sqrt2)^2

    double de1 = 5*(si2-sj2)*is2;
    double de2 = 4*r2*is4;
    double E   = e1*e2*e3;
    dSsi       = E*( de2*si - (2 + de1)*isi );
    dSsj       = E*( de2*sj +    + de1 /sj  );
    dSr        = E*( 2/r2     -4*is2        );

    return E;

    /*
    double e1si =  5*(-si2 + sj2)/(si*(si2 + sj2));
    double e1sj =  5*( si2 - sj2)/(sj*(si2 + sj2));
    double e2si =  4*r2*si/sq(si2 + sj2);
    double e2sj =  4*r2*sj/sq(si2 + sj2);
    double e2r  = -4/(si2 + sj2);
    double e3si = -2/si;
    double e3r  =  2/r2;
    double E  = e1*e2*e3;
    dSsi     = E*( e2si + e3si + e1si );
    dSsj     = E*( e2sj +      + e1sj );
    dSr      = E*( e2r  + e3r         );
    */
}

inline double addDensOverlapGauss_S( const Vec3d& dR, double si, double sj, double amp, Vec3d& f, double& fsi, double& fsj ){
    double r2 = dR.norm2();
    double si2  = si*si;
    double sj2  = sj*sj;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;


    //double amp_ = amp * const_K_eVA * 2.2*( 1.5*( (si2 + sj2)/(si2*sj2) ) - 4.9/( si2 + sj2 ) );
    //double amp_ = amp * const_K_eVA * ( 3.3*s2/(si2*sj2) - 10.78*is2 );
    //double amp_ = amp * const_K_eVA * ( 3.3*s2*s2 - 10.78*si2*sj2 )/(si2*sj2*is2)  );


    double fr,fi,fj;
    double E = DensOverlapGauss_S( r2,amp,si,sj,    fr,fi,fj, si2, sj2, is2, is4 );
    //double E = DensOverlapGauss_Snorm( r2,amp,si,sj,    fr,fi,fj, si2, sj2, is2, is4 );
    //printf(  " dR.x %g amp %g si %g sj %g -> %g \n ", dR.x, amp, si, sj, E );
    fsi += fi;
    fsj += fj;
    f.add_mul( dR, fr );
    return E;
}

inline double addDensOverlapGauss_P( const Vec3d& dR, double si, double sj, double amp, Vec3d& f, double& fsi, double& fsj ){
    double r2 = dR.norm2();
    double si2  = si*si;
    double sj2  = sj*sj;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    double fr,fi,fj;
    double E = DensOverlapGauss_P( r2,amp,si,sj,    fr,fi,fj, si2, sj2, is2, is4 );
    fsi += fi;
    fsj += fj;
    f.add_mul( dR, fr );
    return E;
}

/*
inline double PauliGauss_S( double r, double s, double* params, double& fr, double& fs ){
    double a = params[0];

    double w  = 2.*(si*sj)*is2;

    double b = params[1];
    double c = params[2];
    double r2 = r*r;
    double s2 = s*s;
    double invs = 1/( c + s2 );
    double beta = -b*invs;
    double E =  a * exp( beta*r2 );
    fr = E  * beta * 2;   // dont forget multiply by vec_r
    fs = fr * r2*invs2;   // dont forget muliply by s
    return E;
}

inline double PauliGauss_P( double r, double s, double& fr, double& fs ){
    double a = params[0];
    double b = params[1];
    double c = params[2];
    double d = params[3];
    double e = params[4];
    //double r2 = r*r;
    double s2 = s*s;

    double r_   = (r-c*s);
    double r_2  = r_*r_;
    double invs = 1/( c + s2 );
    double beta = -d*invs;

    double e1 = 2*(b*b+s*s)/(s*b);

    double e2 =  exp( -beta*r_2 );

    double E =  a * e1  * r_2 *  e2;
    return

}
*/

// ======================================================
// =========   Pauli Repulsion - Older
// ======================================================

//inline double getDeltaTGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getDeltaTGauss( double r2, double si, double sj,  double& dTr, double& dTsi, double& dTsj,
    double isi2, double isj2, double s2, double is2, double is4
  ){

    // see eq.4 in http://aip.scitation.org/doi/10.1063/1.3272671
    // and eq.2 in http://doi.wiley.com/10.1002/jcc.21637

    double B   =  4.*( 3.*s2 - 4.*r2 )*is4*is2;
    dTsi = const_K_eVA*( -3.*isi2*isi2 + B )*si;
    dTsj = const_K_eVA*( -3.*isj2*isj2 + B )*sj;
    //f.add_mul( dR, const_K_eVA * (8.*is4) );
    dTr  = const_K_eVA*(8.*is4);
    //double T = const_K_eVA * ( 1.5*s2*isi2*isj2 -2.*( 3.*s2 - 2.*r2 )*is4 );  // (t11+t22-t12/s12) eq.4 in http://aip.scitation.org/doi/10.1063/1.3272671
    double T = const_K_eVA * ( -2.*( 3.*s2 - 2.*r2 )*is4 );                     // WARRNING : this is just debugging - ommit t11,t22 from (  )

    // This is kinetic energy change normalized by overlap
    // see e.g.   eq.4   in http://aip.scitation.org/doi/10.1063/1.3272671
    //       or   eq.3   in https://link.aps.org/doi/10.1103/PhysRevLett.99.185003
    // see also   eq.3,4 in https://iopscience.iop.org/article/10.1088/1367-2630/14/8/083023

    //printf( "getDeltaTGauss T %g r %g si %g sj %g ", T, sqrt(r2), si, sj );
    //printf( "getDeltaTGauss: e1 %g \n", 1.5*s2*isi2*isj2 );
    //printf( "getDeltaTGauss: e2 %g \n", -2.*( 3.*s2 - 2.*r2 )*is4 );
    //printf( "getDeltaTGauss: T  %g \n", const_K_eVA * ( 1.5*s2*isi2*isj2 -2.*( 3.*s2 - 2.*r2 )*is4 ) );
    //printf( "getDeltaTGauss: B  %g \n", B    );
    //printf( "getDeltaTGauss: dTsi %g \n", dTsi );
    //printf( "getDeltaTGauss: dTsj %g \n", dTsj );
    //printf( "getDeltaTGauss: dTr  %g \n", dTr  );
    //printf( "getDeltaTGauss: fr   %g \n", dTr*sqrt(r2) );

    return T;
}


//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getOverlapSGauss( double r2, double si, double sj, double& dSr, double& dSsi, double& dSsj,
    double si2, double sj2, double is2, double is4
){

    double a2   = 2.*(si*sj)*is2;
    double a    = sqrt(a2);
    double e1   = a2*a;
    double e2   = exp( -r2*is2 );

    double f1   = 3.*a  * (si2-sj2)*is4;
    double f2   = 2.*e2 * r2*is4;

    dSsi = e1*f2*si - e2*f1*sj;
    dSsj = e1*f2*sj + e2*f1*si;
    //f.add_mul( dR, e1*e2*(-2.*is2) );
    dSr  = e1*e2*(-2.*is2);

    //printf( "getOverlapSGauss r %g s(%g,%g) S %g ", sqrt(r2), si, sj, e1*e2 );
    return e1 * e2;
}

inline double getDeltaTGauss( double r2, double si, double sj,  double& dTr, double& dTsi, double& dTsj ){
    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    //printf( "getDeltaTGauss: r2, si2, sj2, s2, s4 %g %g %g %g %g \n", r2, si2, sj2, s2, is4 );
    return getDeltaTGauss  ( r2, si, sj, dTr, dTsi, dTsj, 1./si2,1./sj2,s2,is2,is4 );
}

//inline double getOverlapSGauss( const Vec3d& dR, Vec3d& f, double si, double sj, double& fsi, double& fsj ){
inline double getOverlapSGauss( double r2, double si, double sj, double& dSr, double& dSsi, double& dSsj ){
    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;
    //printf( "getDeltaSGauss: r2, si2, sj2, s2, s4 %g %g %g %g %g \n", r2, si2, sj2, s2, is4 );

    return getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj, si2,sj2,is2,is4 );
}

inline double PauliSGauss_anti( double S, double& fS, double rho ){
    double S2    = S*S;
    double D     = 1./(1.+S2);
    double SDrho = rho*S*D;
    fS  = 2.*D*SDrho;
    return   S*SDrho;
    //printf( "EpAnti %g a  %g \n", S*SDrho, rho );
    //printf( "rho %g \n", rho);
    //return rho;
}

inline double PauliSGauss_syn( double S, double& fS, double rho ){
    double S2   = S*S;
    double D    = 1./(1.+S2);
    double Dm   = 1./(1.-S2);
    //double S2 = S*S;
    double rhom = 1-rho;

    double SDm    = S*Dm;
    double SDrhom = S*D*rhom;

    fS   = 2.*( SDm*Dm  + SDrhom*D );
    return  S*( SDm     + SDrhom   );

}

inline double addPauliGauss( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, bool anti, const Vec3d& KRSrho ){

    double r2 = dR.norm2();
    
    //WARRNING : switched off  KRSrho-scalling 
    //r2 *= KRSrho.x*KRSrho.x;
    //si *= KRSrho.y;
    //sj *= KRSrho.y;

    //si*=M_SQRT2; sj*=M_SQRT2;  // DEBUG - with this it is OK

    double r = sqrt(r2 + 1e-16);

    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;

    double dTr,dTsi,dTsj;
    double dSr,dSsi,dSsj;
    double T = getDeltaTGauss  ( r2, si, sj, dTr, dTsi, dTsj, 1./si2, 1./sj2, s2,is2,is4 );
    double S = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj,    si2,    sj2    ,is2,is4 );

    

    // Eq.3 in http://aip.scitation.org/doi/10.1063/1.3272671
    double eS,fS;
    if(anti){ eS = PauliSGauss_anti( S, fS, KRSrho.z ); }
    else    { eS = PauliSGauss_syn ( S, fS, KRSrho.z ); }
    double TfS = T*fS;
    fsi +=          (dTsi*eS + TfS*dSsi)*KRSrho.y;
    fsj +=          (dTsj*eS + TfS*dSsj)*KRSrho.y;
    f.add_mul( dR, -(dTr *eS + TfS*dSr )*KRSrho.x*KRSrho.x ); // second *KRSrho.x because dR is not multiplied

    //printf( "addPauliGauss r2, si, sj %g %g %g \n", r2, si, sj );
    //printf( "addPauliGauss T,dTr,dTsi,dTsj %g %g %g %g \n", T,dTr,dTsi,dTsj );
    //printf( "addPauliGauss S,dSr,dSsi,dSs  %g %g %g %g \n", S,dSr,dSsi,dSsj );
    //printf( "addPauliGauss S, eS, fS  %g %g %g \n", S, eS, fS  );
    //printf( "addPauliGauss fr1, fr2 fr %g %g %g %g \n", dTr*eS, TfS*dSr, (dTr *eS + TfS*dSr )*KRSrho.x, dR.z );
    //printf( "---------------------- \n" );

    //printf( " E %g T %g eS %g S %g ", T*eS, T, eS, S );
    //printf( "addPauliGauss() E %g T %g S %g | r %g anti %i \n", T*eS, T, S, sqrt(r2), (int)anti  );

    double E = T * eS;

    printf( "Epaul_ij %g T %g S %g r %g si %g sj %g \n",  E, T, S, r, si, sj  );

    return E;

    //return T * eS;
    //return eS;
    //return S;
    //return eS*-3;
    //return T;
    //return T;
}

//inline double addPauliGauss_New( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, bool anti, const Vec3d& KRSrho ){
inline double addPauliGauss_New( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, int spin, const Vec3d& KRSrho, double sc=1.0 ){
    double r2         = dR.norm2() + 1e-8;  // for r=0 there are numercial instabilities

    constexpr const double Hartree2eV = 27.211386245988;
    constexpr const double A2bohr     = 1/0.5291772105638411;

    /*
    r2*=(A2bohr*A2bohr);
    si*=A2bohr;
    sj*=A2bohr;
    double KR2=KRSrho.x*KRSrho.x;
    r2*=KR2;  si*=KRSrho.y; sj*=KRSrho.y;
    */

    double KR2=A2bohr*KRSrho.x; KR2*=KR2;
    const double KS =A2bohr*KRSrho.y;
    si*=KS; sj*=KS; r2*=KR2;

    double si2        = si*si;
    double sj2        = sj*sj;
    double si2sj2     = si2 + sj2;
    double invsi2sj2  = 1/si2sj2;
    double invsi2sj22 = invsi2sj2*invsi2sj2;
    double invsi2sj23 = invsi2sj2*invsi2sj22;
    double denom_sij  = si*sj*invsi2sj2;
    double si4sj4     = si2*si2 - sj2*sj2; 
    double invsj      = 1/sj;
    double invsi      = 1/si; 
    double invsj2     = invsj*invsj;
    double invsi2     = invsi*invsi; 

    //double r     = sqrt(r2 + 1e-16);
    //double expr  = exp(-r2*invsi2sj2);
    //double expr2 = expr*expr;
    //double expr2  = exp(-2*r2*invsi2sj2);

    double r2_4   =  4*r2;

    // ------- Kinetic Energy Difference
    double DT      = 1.5*si2sj2*invsi2*invsj2 -      (6*si2sj2 - r2_4)*invsi2sj22;
    double dDT_dsi =  -3*invsi2*invsi         + 4*si*(3*si2sj2 - r2_4)*invsi2sj23;   // TESTED with eff.py .check_DerivsPauli()
    double dDT_dsj =  -3*invsj2*invsj         + 4*sj*(3*si2sj2 - r2_4)*invsi2sj23;   // TESTED with eff.py .check_DerivsPauli() 
    double dDT_dr  =   8*invsi2sj22;      // missing 'r' it is in |dR|               // TESTED with eff.py .check_DerivsPauli()

    // ------- Overlap  ..... actually S22 = 2*S**2
    double S22      = 8*denom_sij*denom_sij*denom_sij*exp(-2*r2*invsi2sj2);         
    double dS22_dsi = S22*( -3*si4sj4 + r2_4*si2 )*invsi2sj22*invsi;                // TESTED with eff.py .check_DerivsPauli()
    double dS22_dsj = S22*( +3*si4sj4 + r2_4*sj2 )*invsi2sj22*invsj;                // TESTED with eff.py .check_DerivsPauli()
    double dS22_dr  = -4*S22*invsi2sj2;   // missing 'r' it is in |dR|              // TESTED with eff.py .check_DerivsPauli()

    //fsi=dDT_dsi;    fsj=dDT_dsj;  f=dR*dDT_dr;   return DT;  
    //fsi=dS22_dsi;   fsj=dS22_dsj; f=dR*dS22_dr;  return S22;

    double rho = KRSrho.z;

    double E=0, dE_dDT=0, dE_dS22=0;
    if(spin<=0){
        double invS22m1 = 1/(S22+1);
        E       += - rho*DT*S22  *invS22m1;
        dE_dDT  += -(rho*   S22 )*invS22m1;
        dE_dS22 += -(rho*DT     )*invS22m1*invS22m1;
    }
    if(spin>=0){
        double invS222m1 = 1/( S22*S22-1 );
        //printf( "DEBUG S22*DT %g invS222m1 %g (-rho*S22+rho-2) %g \n", S22*DT, invS222m1, (-rho*S22+rho-2) );
        E       +=   S22 * DT * ( -rho*S22                     + rho-2 ) *invS222m1;
        dE_dDT  += - S22 *      (  rho*S22                     - rho+2 ) *invS222m1;
        dE_dS22 +=      -  DT * (      S22*(S22*(rho-2)-2*rho) + rho-2 ) *invS222m1*invS222m1;
    }

    // double E, dE_dDT, dE_dS22;
    // if( anti ){
    //     double invS22m1 = 1/(S22+1);
    //     E       = - rho*DT*S22  *invS22m1;
    //     dE_dDT  = -(rho*   S22 )*invS22m1;
    //     dE_dS22 = -(rho*DT     )*invS22m1*invS22m1;
    // }else{
    //     double invS222m1 = 1/( S22*S22-1 );
    //     //printf( "DEBUG S22*DT %g invS222m1 %g (-rho*S22+rho-2) %g \n", S22*DT, invS222m1, (-rho*S22+rho-2) );
    //     E       =   S22 * DT * ( -rho*S22                     + rho-2 ) *invS222m1;
    //     dE_dDT  = - S22 *      (  rho*S22                     - rho+2 ) *invS222m1;
    //     dE_dS22 =      -  DT * (      S22*(S22*(rho-2)-2*rho) + rho-2 ) *invS222m1*invS222m1;
    // }

    //fsi       += dE_dS22 * dS22_dsi + dE_dDT * dDT_dsi;
    //fsj       += dE_dS22 * dS22_dsj + dE_dDT * dDT_dsj;
    //double fr  = dE_dS22 * dS22_dr  + dE_dDT * dDT_dr;

    // E         *= Hartree2eV;
    // fsi       += (dE_dS22 * dS22_dsi + dE_dDT * dDT_dsi)*Hartree2eV*A2bohr*KRSrho.y;
    // fsj       += (dE_dS22 * dS22_dsj + dE_dDT * dDT_dsj)*Hartree2eV*A2bohr*KRSrho.y;
    // double fr  = (dE_dS22 * dS22_dr  + dE_dDT * dDT_dr )*Hartree2eV*A2bohr*A2bohr*KR2;

    E         *= Hartree2eV*sc;
    fsi       += (dE_dS22 * dS22_dsi + dE_dDT * dDT_dsi)*Hartree2eV*-KS *sc;
    fsj       += (dE_dS22 * dS22_dsj + dE_dDT * dDT_dsj)*Hartree2eV*-KS *sc;
    double fr  = (dE_dS22 * dS22_dr  + dE_dDT * dDT_dr )*Hartree2eV*KR2*sc;

    f.add_mul( dR, fr  );

    //printf( "r %g si %g sj %g DT %g S22 %g E %g anti(%i) \n", sqrt(r2), si,sj, DT,S22, E, anti );
    return E;
    
}

inline double addPauliGaussVB( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj ){

    //if (anti) return 0;

    double r2 = dR.norm2();
    double r = sqrt(r2 + 1e-16);

    double si2  = si*si;
    double sj2  = sj*sj;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;

    double dTr,dTsi,dTsj;
    double dSr,dSsi,dSsj;
    double T = getDeltaTGauss  ( r2, si, sj, dTr, dTsi, dTsj, 1./si2, 1./sj2, s2,is2,is4 );
    double S = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj,    si2,    sj2    ,is2,is4 );

    // Delta Ek = S^2/(1-S^2) * DeltaT   Eq.2 in  http://doi.wiley.com/10.1002/jcc.21637
    //double S2  = S*S;
    double D   = 1./(1.-S*S);
    double SD  = S*D;
    double fS  = 2.*( D*SD );
    double eS  =      S*SD  ;

    double TfS = T*fS;
    fsi +=          (dTsi*eS + TfS*dSsi);
    fsj +=          (dTsj*eS + TfS*dSsj);
    f.add_mul( dR, -(dTr *eS + TfS*dSr ) ); // second *KRSrho.x because dR is not multiplied

    //printf( " E %g T %g eS %g S %g ", T*eS, T, eS, S );

    return T * eS;

}

/*
inline double addPauliGaussS( const Vec3d& dR, double si, double sj, Vec3d& f, double& fsi, double& fsj, bool anti, const Vec3d& KRSrho ){

    double r2 = dR.norm2();
    r2 *= KRSrho.x*KRSrho.x;
    si *= KRSrho.y;
    sj *= KRSrho.y;
    double r = sqrt(r2+1e-16);

    double si2  = si*si;
    double sj2  = sj*sj;
    //double isi2 = 1./si2;
    //double isj2 = 1./sj2;
    double s2   = si2 + sj2;
    double is2  = 1./s2;
    double is4  = is2 * is2;

    double dTr,dTsi,dTsj;
    double dSr,dSsi,dSsj;
    double T = getDeltaTGauss  ( r2, si, sj, dTr, dTsi, dTsj, 1./si2, 1./sj2,s2,is2,is4 );
    double S = getOverlapSGauss( r2, si, sj, dSr, dSsi, dSsj,    si2,    sj2,is2,is4 );

    double eS,fS;
    if(anti){
        eS = PauliSGauss_anti( S, fS, KRSrho.z );
    }else{
        eS = PauliSGauss_anti( S, fS, KRSrho.z );
    }

    double TfS = T*fS;
    fsi +=         (dTsi*eS + TfS*dSsi)*KRSrho.y;
    fsj +=         (dTsj*eS + TfS*dSsj)*KRSrho.y;
    f.add_mul( dR, (dTr *eS + TfS*dSr )*KRSrho.x*KRSrho.x ); // second *KRSrho.x because dR is not multiplied

    //printf( "addPauliGauss r2, si, sj %g %g %g \n", r2, si, sj );
    //printf( "addPauliGauss T,dTr,dTsi,dTsj %g %g %g %g \n", T,dTr,dTsi,dTsj );
    //printf( "addPauliGauss S,dSr,dSsi,dSs  %g %g %g %g \n", S,dSr,dSsi,dSsj );
    //printf( "addPauliGauss S, eS, fS  %g %g %g \n", S, eS, fS  );
    //printf( "addPauliGauss fr1, fr2 fr %g %g %g %g \n", dTr*eS, TfS*dSr, (dTr *eS + TfS*dSr )*KRSrho.x, dR.z );
    //printf( "---------------------- \n" );

    return T * eS;

}
*/



inline void PauliCoreElec_Orig(double rc, double re2, double *epauli, double *frc, double *fre2, double PAULI_CORE_A, double PAULI_CORE_B, double PAULI_CORE_C){
    double E, dEdrc, dEdre2, rcsq, ssq;
    rcsq = rc * rc;
    ssq = re2 * re2;
    // (eq.5) in http://dx.doi.org/10.1016/j.mechmat.2015.02.008
    // A*exp( -B*r^2/(c+s^2)  )
    E      =      PAULI_CORE_A *                             exp(-PAULI_CORE_B * rcsq / (ssq + PAULI_CORE_C));
    dEdrc  = -2 * PAULI_CORE_A * PAULI_CORE_B * rc *         exp(-PAULI_CORE_B * rcsq / (ssq + PAULI_CORE_C)) /  (PAULI_CORE_C + ssq);
    dEdre2 = 2 *  PAULI_CORE_A * PAULI_CORE_B * re2 * rcsq * exp(-PAULI_CORE_B * rcsq / (ssq + PAULI_CORE_C)) / ((PAULI_CORE_C + ssq)*(PAULI_CORE_C + ssq));
    *epauli += E;
    *frc -= dEdrc;
    *fre2 -= dEdre2;
}

inline void PauliCoreElec(double r, double re2, double& epauli, double& frc, double& fre2, double A, double B, double C){
    double rcsq   = r * r;
    double ssq    = re2 * re2;
    double invssC = 1/(ssq + C);
    double expr   =  exp(-B * rcsq *invssC );
    double ABr    = 2 * A * B * r;
    // (eq.5) in http://dx.doi.org/10.1016/j.mechmat.2015.02.008
    // A*exp( -B*r^2/(c+s^2)  )
    double E       =  A   *           expr;
    double dE_drc  = -ABr *           expr * invssC;
    double dE_dre2 =  ABr * r * re2 * expr * invssC*invssC;
    epauli += E;
    frc    -= dE_drc;
    fre2   -= dE_dre2;
}

inline void PauliCorePElec_Orig(double rc, double re2, double *epauli, double *frc, double *fre2, double PAULI_CORE_P_A, double PAULI_CORE_P_B, double PAULI_CORE_P_C, double PAULI_CORE_P_D, double PAULI_CORE_P_E) {
    double E, dEdrc, dEdre2;

    //(eq.5) in http://dx.doi.org/10.1016/j.mechmat.2015.02.008
    // A * (2/(B/s + s/B) )^5 * (r-C*s)^2 *  exp( -D*(r-C*s)^2/(E+s^2) )          
    E     = PAULI_CORE_P_A * pow((2.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B)), 5.0) *   pow((rc - PAULI_CORE_P_C * re2), 2.0) * exp(-PAULI_CORE_P_D * pow((rc - PAULI_CORE_P_C * re2), 2.0) / (PAULI_CORE_P_E + re2 * re2));
    dEdrc = PAULI_CORE_P_A * pow((2.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B)), 5.0) *      (rc - PAULI_CORE_P_C * re2) * 2.0  * exp(-PAULI_CORE_P_D * pow((rc - PAULI_CORE_P_C * re2), 2.0) / (PAULI_CORE_P_E + re2 * re2))
            + E *                                                                                   ((rc - PAULI_CORE_P_C * re2)* 2.0 *      -PAULI_CORE_P_D                                         / (PAULI_CORE_P_E + re2 * re2));

    dEdre2 = E * (-5.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B) * (-PAULI_CORE_P_B / (re2 * re2) + 1.0 / PAULI_CORE_P_B))
            + PAULI_CORE_P_A *  pow((2.0 / (PAULI_CORE_P_B / re2 + re2 / PAULI_CORE_P_B)), 5.0) *2.0 * (rc - PAULI_CORE_P_C * re2) * (-PAULI_CORE_P_C) * exp(-PAULI_CORE_P_D * pow((rc - PAULI_CORE_P_C * re2), 2.0) / (PAULI_CORE_P_E + re2 * re2)) 
            + E * (2.0 * PAULI_CORE_P_D * (rc - PAULI_CORE_P_C * re2) * (PAULI_CORE_P_C * PAULI_CORE_P_E + rc * re2) / pow((PAULI_CORE_P_E + re2 * re2), 2.0));

    *epauli += E;
    *frc -= dEdrc;
    *fre2 -= dEdre2;
}

inline void PauliCorePElec(double rc, double re2, double& epauli, double& frc, double& fre2, double A, double B, double C, double D, double E_ ) {
    double invBre  =  2.0 / (B / re2 + re2 / B);
    double invBre2 = invBre*invBre;
    double invBre5 = invBre2*invBre2*invBre; 
    double rrC      = rc - C * re2;
    double rrC2     = rrC*rrC;
    double invErr   = 1/(E_ + re2 * re2);
    double expArg   = -D * rrC2 * invErr;
    double expr     = exp( expArg );

    //(eq.5) in http://dx.doi.org/10.1016/j.mechmat.2015.02.008
    // A * (2/(B/s + s/B) )^5 * (r-C*s)^2 *  exp( -D*(r-C*s)^2/(E+s^2) )          
    //   double e       = A * invBre5 * rrC2 * expr;
    //   double dE_drc  = A * invBre5 * 2.0 * rrC * expr   +    e * D * -2.0 * rrC * invErr ;
    //   double dE_dre2 = e * (-2.5 * invBre *     ( 1.0/B - B/(re2*re2) )                                     // e * (-2.5 * invBre *     ((re2/B) - (B/re2)))   / re2
    //                  + A * invBre5 * -2.0 * rrC * C * expr
    //                  + e * (2.0 * D * rrC * (C * E + rc * re2) * invErr* invErr;

    double erc     = A * invBre5 * rrC * expr;
    double e       = erc*rrC;
    double dE_drc  = 2.0 * erc * ( 1   + expArg  );
    double dE_dre2 = erc * ( -2.0*C + rrC * ( -2.5 * invBre  * ( 1.0/B - B/(re2*re2) )                                     // e * (-2.5 * invBre *     ((re2/B) - (B/re2)))   / re2   
                   +                  rrC * (  2.0 * D *  (C * E_ + rc * re2) * invErr* invErr ) ) );

    epauli += e;
    frc    -= dE_drc;
    fre2   -= dE_dre2;
}


///  @}

#endif



