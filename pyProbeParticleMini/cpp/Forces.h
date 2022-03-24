
#ifndef Forces_h
#define Forces_h

#include "fastmath_light.h"
#include "Vec2.h"
#include "Vec3.h"

#define COULOMB_CONST  14.3996448915f

//const double kcoulomb   = 14.3996448915; 
//const double R2SAFE     = 1.0e-8;

#define RSAFE   1.0e-4f
#define R2SAFE  1.0e-8f
#define F2MAX   10.0f

void sum(int n, Vec3d* ps, Vec3d& psum){ for(int i=0;i<n;i++){ psum.add(ps[i]); } };

void sumTroq(int n, Vec3d* fs, Vec3d* ps, const Vec3d& cog, const Vec3d& fav, Vec3d& torq){
    for(int i=0;i<n;i++){  torq.add_cross(ps[i]-cog,fs[i]-fav);  }
    //for(int i=0;i<n;i++){  torq.add_cross(ps[i],fs[i]);  }
}

void checkForceInvariatns( int n, Vec3d* fs, Vec3d* ps, Vec3d& cog, Vec3d& fsum, Vec3d& torq ){
    cog =Vec3dZero;
    fsum=Vec3dZero;
    torq=Vec3dZero;
    double dw = 1./n;
    sum(n, ps, cog ); cog.mul(dw);
    sum(n, fs, fsum); //cog.mul(dw);
    sumTroq(n, fs, ps, cog, fsum*dw, torq );
}

inline double boxForce1D(double x, double xmin, double xmax, double k){
    double f=0;
    if(k<0) return 0;
    if(x>xmax){ f+=k*(xmax-x); }
    if(x<xmin){ f+=k*(xmin-x); }
    return f;
}

inline void boxForce(const Vec3d p, Vec3d& f,const Vec3d& pmin, const Vec3d& pmax, const Vec3d& k){
    f.x+=boxForce1D( p.x, pmin.x, pmax.x, k.x);
    f.y+=boxForce1D( p.y, pmin.y, pmax.y, k.y);
    f.z+=boxForce1D( p.z, pmin.z, pmax.z, k.z);
}

inline double evalCos2(const Vec3d& hi, const Vec3d& hj, Vec3d& fi, Vec3d& fj, double k, double c0){
    double c    = hi.dot(hj) - c0;
    double dfc  =  k*-2*c;
    fi.add_mul(hj,dfc);
    fj.add_mul(hi,dfc);
    return k*c*c;
}

inline double evalCos2_o(const Vec3d& hi, const Vec3d& hj, Vec3d& fi, Vec3d& fj, double k, double c0){
    double c    = hi.dot(hj) - c0;
    double dfc  =  k*-2*c;
    double dfcc = -c*dfc;
    fi.add_lincomb( dfc,hj, dfcc,hi );
    fj.add_lincomb( dfc,hi, dfcc,hj );
    return k*c*c;
}

inline double evalCosHalf(const Vec3d& hi, const Vec3d& hj, Vec3d& fi, Vec3d& fj, double k, Vec2d cs ){
    Vec3d h; h.set_add( hi, hj );
    double c2 = h.norm2()*0.25;               // cos(a/2) = |ha+hb|
    double s2 = 1-c2;
    double c = sqrt(c2);
    double s = sqrt(s2);
    cs.udiv_cmplx({c,s});
    double E         =  k*( 1 - cs.x );  // just for debug ?
    double fr        = -k*(     cs.y );
    fr /= 2*c*s;  // 1/sin(2a)
    c2 *=-2*fr;
    Vec3d fa,fb;
    fi.set_lincomb( fr,h,  c2,hi ); 
    fj.set_lincomb( fr,h,  c2,hj );
    return E;
}


// ================= BEGIN:  From ProbeParticle.cpp

// radial spring constrain - force length of vector |dR| to be l0
inline Vec3d forceRSpring( const Vec3d& dR, double k, double l0 ){
    double l = sqrt( dR.norm2() );
    Vec3d f; f.set_mul( dR, k*( l - l0 )/l );
    return f;
}

inline Vec3d forceSpringRotated( const Vec3d& dR, const Vec3d& Fw, const Vec3d& Up, const Vec3d& R0, const Vec3d& K ){
    // dR - vector between actual PPpos and anchor point (in global coords)
    // Fw - forward diraction of anchor coordinate system (previous bond direction; e.g. Tip->C for C->O) (in global coords)
    // Up - Up vector --,,-- ; e.g. x axis (1,0,0), defines rotation of your tip (in global coords)
    // R0 - equlibirum position of PP (in local coords)
    // K  - stiffness (ka,kb,kc) along local coords
    // return force (in global coords)
    Mat3d rot; Vec3d dR_,f_,f;
    rot.fromDirUp( Fw*(1/Fw.norm()), Up );  // build orthonormal rotation matrix
    rot.dot_to  ( dR, dR_   );              // transform dR to rotated coordinate system
    f_ .set_mul ( dR_-R0, K );              // spring force (in rotated system)
    // here you can easily put also other forces - e.g. Torsion etc. 
    rot.dot_to_T( dR_, f );                 // transform force back to world system
    return f;
}

// Lenard-Jones force between two atoms a,b separated by vector dR = Ra - Rb
inline double addAtomLJ( const Vec3d& dR, Vec3d& fout, double c6, double c12 ){
    double ir2  = 1.0/ ( dR.norm2( ) + R2SAFE ); 
    double ir6  = ir2*ir2*ir2;
    double E6   = c6  * ir6;
    double E12  = c12 * ir6*ir6;
    //return dR * ( ( 6*ir6*c6 -12*ir12*c12 ) * ir2  );
    fout.add_mul( dR , ( 6*E6 -12*E12 ) * ir2 );
    //fout.add_mul( dR , -12*E12 * ir2 );
    //fout.add_mul( dR , 6*E6 * ir2 );
    //printf(" (%g,%g,%g)  (%g,%g)  %g \n", dR.x,dR.y,dR.z, c6, c12,  E12 - E6);
    //printf(" (%g,%g,%g)  %f %f  (%g,%g,%g) \n", dR.x,dR.y,dR.z, c6, c12,  fout.x,fout.y,fout.z);
    return E12 - E6;
}

// Lenard-Jones force between two atoms a,b separated by vector dR = Ra - Rb
inline double addAtomVdW( const Vec3d& dR, Vec3d& fout, double c6 ){
    double r2 = dR.norm2(); r2*=r2; r2*=r2;
    //fout.add_mul( dR , 6*c6 /( r2 + 1.0 ) );
    //fout.add_mul( dR , 6*c6 /( r2 + 60*c6 ) );
    fout.add_mul( dR , 6*c6 /( r2 + 180*c6 ) );
    return 0;
}

// Morse force between two atoms a,b separated by vector dR = Ra - Rb
inline double addAtomMorse( const Vec3d& dR, Vec3d& fout, double r0, double eps, double alpha ){
    double r     = sqrt( dR.norm2() + R2SAFE );
    double expar = exp( alpha*(r-r0));
    double E     = eps*( expar*expar - 2*expar );
    double fr    = eps*2*alpha*( expar*expar - expar );
    fout.add_mul( dR, fr/r );
    return E;
}

// coulomb force between two atoms a,b separated by vector dR = R1 - R2, with constant kqq should be set to kqq = - k_coulomb * Qa * Qb 
inline double addAtomCoulomb( const Vec3d& dR, Vec3d& fout, double kqq ){
    double ir2   = 1.0/( dR.norm2() + R2SAFE );
    double ir    = sqrt(ir2); 
    double E     = ir * kqq;
    fout.add_mul( dR , E * ir2 );
    //printf("(%g,%g,%g) %g %g (%g,%g,%g)", dR.x,dR.y,dR.z, kqq, ir, fout.x,fout.y,fout.z );
    return E;
}

// ================= END: From ProbeParticle.cpp



inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq ){
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*qq*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceMorse( const Vec3d& dp, Vec3d& f, double r0, double eps, double beta ){
    const double R2ELEC = 1.0;
    double r     = sqrt( dp.norm2()+R2SAFE );
    double expar = exp ( beta*(r-r0) );
    double fr    = eps*2*beta*( expar*expar - expar );
    f.add_mul( dp, fr/r );
}

inline void addAtomicForceMorseQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq, double alpha ){
    const double R2ELEC = 1.0;
    double r     = sqrt( dp.norm2()+R2SAFE );
    double expar = exp( alpha*(r-r0));
    double fr    = eps*2*alpha*( expar*expar - expar ) + COULOMB_CONST*qq/( r*r + R2ELEC );
    f.add_mul( dp, fr/r );
}

inline void addAtomicForceQ( const Vec3d& dp, Vec3d& f, double qq ){
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double fr   = ( ir*qq*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceLJ( const Vec3d& dp, Vec3d& f, double r0, double eps ){;
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps )*ir2;
    f.add_mul( dp, fr );
}

inline void addAtomicForceExp( const Vec3d& dp, Vec3d& f, double r0, double eps, double alpha ){
    double r    = sqrt(dp.norm2() + R2SAFE );
    double E    = eps*exp( alpha*(r-r0) );
    double fr   = alpha*E/r;
    f.add_mul( dp, fr );
}

inline Vec3d REQ2PLQ( const Vec3d& REQ, double alpha ){
    double eps   = REQ.y;
    double expar = exp(-alpha*REQ.x);
    double CP    =    eps*expar*expar;
    double CL    = -2*eps*expar;
    return (Vec3d){ CP, CL, REQ.z };
}

inline Vec3d REnergyQ2PLQ( const Vec3d& REQ, double alpha ){
    return REQ2PLQ( {REQ.x, sqrt(REQ.y), REQ.z}, alpha );
}

inline Vec3d getForceSpringPlane( const Vec3d& p, const Vec3d& normal, double c0, double k ){
    double cdot = normal.dot(p) - c0;
    return normal * (cdot * k);
}

inline Vec3d getForceHamakerPlane( const Vec3d& p, const Vec3d& normal, double c0, double e0, double r0 ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    double cdot = normal.dot(p) - c0;
    double ir   = r0/cdot;
    double ir3  = ir*ir*ir;
    double f    = e0*(ir/r0)*ir3*(ir3-1);
    return normal * f;
}

inline Vec3d getForceSpringRay( const Vec3d& p, const Vec3d& hray, const Vec3d& ray0, double k ){
    Vec3d dp; dp.set_sub( p, ray0 );
    double cdot = hray.dot(dp);
    dp.add_mul(hray,-cdot);
    return dp*k;
}

#endif
