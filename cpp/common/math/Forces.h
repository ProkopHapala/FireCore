
#ifndef Forces_h
#define Forces_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

#include "physics_constants.h"

#define RSAFE   1.0e-4f
#define R2SAFE  1.0e-8f
#define F2MAX   10.0f




inline bool clampForce( Vec3d& f, const double f2max ){
    const double f2   = f.norm2();
    const bool bClamp = f2>f2max;
    if( bClamp )[[unlikely]]{
        f.mul( sqrt(f2max/f2) );
    }
    return bClamp;
}


// ================= Trashold functions

double smoothstep_up(double x_, double xmin, double xmax) {
    if      (x_<xmin){ return 0; }
    else if (x_>xmax){ return 1; }
    double x = (x_-xmin)/(xmax-xmin);
    return x*x*(3-2*x);
}

double smoothstep_down(double x_, double xmin, double xmax) {
    if      (x_<xmin){ return 1; }
    else if (x_>xmax){ return 0; }
    double x = (x_-xmin)/(xmax-xmin);
    return 1-x*x*(3-2*x);
}

double R4blob(double r2) { r2=1-r2; return r2*r2; }   // simplest and fastest cutoff function which depends only on r2 (i.e. eliminate need for sqrt)

double R8func(double r2, double R, double Rnod ){
    //This functions should is C1-continuous smoothstep function which is use only r2 (i.e. eliminate need for sqrt), it can be used in 3 ways:
    //  1) smoothstep from 0.0 to 1.0 at interval [Rnod,R] if Rnod<R
    //  2) smoothstep from 1.0 to 0.0 at interval [R,Rnod] if Rnod>R
    //  3) smooth bumb at interval [R1,R2] with peak at R, where 2R^2 = R1^2 + R2^2
    double R2   = R*R;              // 1 mul
    double R2n  = Rnod*Rnod;        // 1 mul
    double y1   =       R2 - r2;    // 1 add
    double y2   = R2n + R2 - y1*y1; // 2 add, 1 mul 
    y2*= R2/(R2+R2n); // rescale to have maximum at y=1   // 1 add, 1 div, 1 mul
    return y2*y2;                   // 1 mul .... in total cost 5 mul, 3 add, 1 div
}

double R8down(double r2, double R, double Rnod ){
    //This functions should is C1-continuous smoothstep function which is use only r2 (i.e. eliminate need for sqrt), it can be used in 3 ways:
    //  1) smoothstep from 0.0 to 1.0 at interval [Rnod,R] if Rnod<R
    double R2   = R*R;             
    double R2n  = Rnod*Rnod;       
    if     ( r2<R2  ) return 1;
    else if( r2>R2n ) return 0;
    double y1   =       R2 - r2;    
    double y2   = R2n + R2 - y1*y1; 
    y2*= R2/(R2+R2n); 
    return y2*y2;   
}

double finiteLorenz( double r2, double w2, double R2cut ){
    if( r2>R2cut ) return 0;
    double fcut = (R2cut-r2);
    return fcut*fcut/(R2cut*R2cut*(r2+w2));
}

double repulsion_R4( Vec3d d, Vec3d& f, double R, double Rcut, double A ){
    // we use R4blob(r) = A * (1-r^2)^2
    // such that at distance r=R we have force f = fmax
    // f = -dR4blob/dr = 4*A*r*(1-r^2) = fmax
    // A = fmax/(4*R*(1-R^2))
    double R2    = R*R;
    double R2cut = Rcut*Rcut;
    double r2 = d.norm2();
    if( r2>R2cut ){ 
        return 0;
        // f = Vec3dZero;
    }else if( r2>R2 ){ 
        double mr2 = R2cut-r2;
        double fr = A*mr2;
        f.add_mul( d, -4*fr );
        return fr*mr2;
    }else{
        double mr2 = R2cut-R2;
        double fr  = A*mr2;
        double r    = sqrt(r2);
        double fmax = 4*R*fr;
        f.add_mul( d, -fmax/r );
        return fmax*(R-r) + fr*mr2;
    }
}

// ================ Zero Torque ( for torsion angles )

// inline void zeroTorque( const Vec3d& a, const Vec3d& b, const Vec3d& fa, const Vec3d& fb, const Mat3d rot, double l, Vec3d& f ){
//     // This function find force f which counteracts torque from fa and fb ( at the end points A,B of bonds a and b, a = A-C, b = B-D )
//     // We consider torsion angle between two bonds a and b between atoms A,B,C,D  (A,B are end points of a, C,D are on the axis, C is at origin adjecent to A and D is adjecent to B in distance l=|C-D| from C)
//     // we consider action of force in plane defined by orientation matrix rot (rot.a,rot.b,rot.c), rot.a is normal to plane, rot.c is axis of torsion, rot.b is up vector such that b is in plane of rot.b and rot.c  
//     // for force component fI in the    plane we have: fI*l = faI*aI                (because fbII=0 by choice of plane)
//     // for force component fN normal to plane we have: fN*l = faN*aI + fbN*(l+bI)
//     double aI  = rot.c.dot(a);  // will be probably negative
//     double bI  = rot.c.dot(b);
//     double faI = rot.b.dot(fa);
//     double faN = rot.a.dot(fa);
//     double fbN = rot.a.dot(fb);
//     double il  = -1./l;
//     double fI  = ( faI*aI              )*il;
//     double fN  = ( faN*aI + fbN*(l+bI) )*il;
//     f = rot.a*fN + rot.b*fI;
// }

// inline void torsion( double fang, Vec3d a, Vec3d b, Quat4d ax,   Vec3d& fa, Vec3d& fb, Vec3d& fc, Vec3d& fd ){
//     Mat3d rota,rotb;
//     rota.fromDirUp( ax.f, a );
//     rotb.fromDirUp( ax.f, b );
//     double ila = 1./rota.b.dot(a);
//     double ilb = 1./rotb.b.dot(b);
//     Vec3d   fa = rota.a*(fang*ila);
//     Vec3d   fb = rotb.a*(fang*ilb);
//     zeroTorque( a, b, fa, fb, rota, ax.w, fd );
//     zeroTorque( b, a, fb, fa, rota, ax.w, fc );
// }

// ================ Angular Forces (MMFF) 

inline double evalBond( const Vec3d& h, double dl, double k, Vec3d& f ){
    double fr = dl*k;
    f.set_mul ( h, fr );
    //f1.add( h );
    //f2.sub( h );
    return fr*dl*0.5;
}

inline double springbound( double x, double l, double k, double& f ){
    double E;
    if(x<0){
        f =-x*k;
        E = 0.5*x*x*k;
    }else if( x>l ){
        x-=l;
        f=-x*k;
        E = 0.5*x*x*k;
    }else{
        f=0;
        E=0;
    }
    return E;
}

inline double spring( double l, Vec2d ls, Vec2d ks, double flim, double& f ){
    double E=0;
    f       =0;
    if  (l>ls.x){
        double dl=l-ls.x;
        f=dl*ks.x;
        //printf( "l(%g)>ls.x(%g) \n", l, ls.x, f  );
        if(f>flim){
            f=flim;
            double dlim = flim/ks.x;
            E = (0.5*dlim*dlim)*ks.x + (dl-dlim)*flim;
        }else{
            E = 0.5*dl*dl*ks.x; 
        }
    }else if(l<ls.y){
        double dl=l-ls.y;
        f=dl*ks.y;
        //printf( "l(%g)<ls.x(%g) \n", l, ls.y, f  );
        if(f<-flim){
            f=-flim;
            double dlim = -flim/ks.y;
            E = (0.5*dlim*dlim)*ks.y - (dl-dlim)*flim;
        }else{
            E = 0.5*dl*dl*ks.y; 
        }
    } 
    return E;
}

inline double evalAngleCosHalf( const Vec3d& h1, const Vec3d& h2, double ir1, double ir2, const Vec2d& cs0, double k, Vec3d& f1, Vec3d& f2 ){
    //printf( " ir1 %g ir2 %g \n", ir1, ir2 );
    // This is much better angular function than evalAngleCos() with just a little higher computational cost ( 2x sqrt )
    Vec3d h; h.set_add( h1, h2 );
    double c2 = h.norm2()*0.25;   // cos(a/2) = |ha+hb|
    double s2 = 1-c2 + 1e-14;      // s2 must be positive number !!!
    double c  = sqrt(c2);
    double s  = sqrt(s2);
    Vec2d cs  = cs0;
    cs.udiv_cmplx({c,s});
    //Vec2d cs{c,s};
    //cs.mul_cmplx(cs0);
    double E         =  k*( 1 - cs.x );  // just for debug ?
    double fr        = -k*(     cs.y );
    c2 *=-2;
    fr /= 4*c*s;   //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    double fr1    = fr*ir1;
    double fr2    = fr*ir2;
    double _nan = E+fr+fr1+fr2;
    //ckeckNaN( 1,1, &_nan, [&]{ printf("evalAngleCosHalf{ E %g fr %g fr1 %g fr2 %g c %g s %g c2 %g s2 %g }", E, fr, fr1, fr2, c, s, c2, s2 ); } );
    f1.set_lincomb( fr1, h,  fr1*c2, h1 );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    f2.set_lincomb( fr2, h,  fr2*c2, h2 );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );
    //printf( "evalAngleCosHalf fr=%g cs2(%g,%g) cso(%g,%g) cs(%g,%g) cs0(%g,%g) \n", fr, c2,s2,  c,s,  cs.x,cs.y, cs0.x,cs0.y );

    return E;
}

inline double evalAngleCos( const Vec3d& h1, const Vec3d& h2, double ir1, double ir2, double K, double c0, Vec3d& f1, Vec3d& f2 ){
    double c = h1.dot(h2);
    f1.set_add_mul( h2,h1,-c );
    f2.set_add_mul( h1,h2,-c );
    double c_   = c-c0;
    double E    =  K*c_*c_;
    //printf( "evalAngleCos() cos=%g E=%g \n", c, E );
    double fang = -K*c_*2;
    f1.mul( fang*ir1 );
    f2.mul( fang*ir2 );
    return E;
}

inline double evalPiAling( const Vec3d& h1, const Vec3d& h2, double ir1, double ir2, double K, Vec3d& f1, Vec3d& f2 ){  // interaction between two pi-bonds
    double c = h1.dot(h2);
    f1.set_add_mul( h2,h1,-c );
    f2.set_add_mul( h1,h2,-c );
    bool sign = c<0; if(sign) c=-c;
    double E    = K*(1-c);
    double fang =  K;
    if(sign)fang=-fang;
    f1.mul( fang );
    f2.mul( fang );
    return E;
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


// ================ Non-Covalent Forces (Lenard-Jones, Morse, etc.)

inline void combineREQ(const Quat4d& a, const Quat4d& b, Quat4d& out){
    out.x=a.x+b.x; // radius
    out.y=a.y*b.y; // epsilon
    out.z=a.z*b.z; // q*q
    out.w=a.w*b.w; // Hbond
}

#define _mixREQ(A,B)    Quat4d{ A.x+B.x, A.y*B.y, A.z*B.z, A.w*B.w }

// evaluate energy and force using Lennard-Jones and Coulomb potential, single precision
inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*qq*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

// evaluate energy and force using Lennard-Jones and Coulomb potential, double precision
inline double addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, const Quat4d& REQ ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + 1e-4 );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*REQ.x*REQ.x;
    double ir6  = ir2_*ir2_*ir2_;
    //double fr   = ( ( 1 - ir6 )*ir6*12*REQ.b + ir*REQ.c*-COULOMB_CONST )*ir2;
    double Eel  = ir*REQ.z*COULOMB_CONST;
    double vdW  = ir6*REQ.y;
    double fr   = ( ( 1 - ir6 )*12*vdW - Eel )*ir2;
    //printf( " (%g,%g,%g) r %g fr %g \n", dp.x,dp.y,dp.z, 1/ir, fr );
    //printf( " r %g fr %g vdW %g Eel %g \n", 1/ir, fr, vdW, Eel  );
    f.add_mul( dp, fr );
    return  ( ir6 - 2 )*vdW + Eel;
}

// evaluate energy and force using Lennard-Jones and Coulomb potential, single precision
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
    E    +=      (u6-2.)*vdW     ;
    F    +=  12.*(u6-1.)*vdW*ir2 ;
    f.set_mul( dp, -F );
    return E;
}

// evaluate energy and force using Lennard-Jones and Coulomb potential, double precision
inline double getLJQ( const Vec3d& dp, Vec3d& f, const Quat4d& REQ, const double R2damp ){
    const double   r2   = dp.norm2();
    // ---- Coulomb
    const double  ir2_ = 1/( r2 + R2damp  );
    double E = COULOMB_CONST*REQ.z*sqrt( ir2_ );
    double F = E*ir2_ ;
    // --- LJ 
    const double  ir2 = 1/r2;
    const double  u2  = REQ.x*REQ.x*ir2;
    const double  u6  = u2*u2*u2;
    const double vdW  = u6*REQ.y;
    E    +=      (u6-2.)*vdW     ;
    F    +=  12.*(u6-1.)*vdW*ir2 ;
    f.set_mul( dp, -F );
    return E;
}

// evaluate energy and force using Lennard-Jones and Coulomb potential and Hydrogen bond pseudo-charges
inline double getLJQH( const Vec3d& dp, Vec3d& f, const Quat4d& REQH, const double R2damp ){
    const double  r2  = dp.norm2();
    double E,F;
    // ---- Electrostatic
    const double ir2_ = 1/( r2 + R2damp  );
    E =  COULOMB_CONST*REQH.z*sqrt( ir2_ );
    F =  E*ir2_ ;
    // --- LJ 
    const double  ir2 = 1/r2;
    const double  u2  = REQH.x*REQH.x*ir2;
    const double  u6  = u2*u2*u2;
    const double vdW  = u6*REQH.y;
    const double   H  = u6*u6* ((REQH.w<0) ? REQH.w*REQH.y : 0.0);  // H-bond correction
    E   +=  (u6-2.)*vdW + H             ;
    F   += ((u6-1.)*vdW + H )*ir2*12 ;
    f.set_mul( dp, -F );
    return E;
}

// evaluate energy and force using Morse and Coulomb potential and Hydrogen bond pseudo-charges
inline double getMorseQH( const Vec3d& dp, Vec3d& f, const Quat4d& REQH, const double K, double R2damp ){
    const double r2    = dp.norm2();
    double E,F;
    // --- Coulomb
    const double ir2_  = 1/( r2 + R2damp );
    E = COULOMB_CONST*REQH.z*sqrt( ir2_ );
    F = E*-ir2_ ;
    // --- Morse
    const double  r  = sqrt( r2   );
    const double  e  = exp( -K*(r-REQH.x) );
    const double  Ae = REQH.y*e;
    const double  He  = e * ( (REQH.w<0) ? REQH.w*REQH.y : 0.0 );  // H-bond correction
    E +=  Ae*(e - 2)   + He;
    F += (Ae*(e - 1)*2 + He)*-K/r;
    f.set_mul( dp, F );
    return E;
}

// evaluate energy and force using Lennard-Jones and Coulomb potential and Hydrogen bond pseudo-charges, with cutoffs
inline double getLJQH_cut( const Vec3d& dp, Vec3d& f, const Quat4d& REQH, const double R2damp, const double Cr2_cut, const double LJr2_cut ){
    // E_cut = e_max * (R/r_cut)^6
    // u2_cut = (e_max/E_cut)^(1/6)
    // u     = R/r
    // if( u2 < u2_cut ){ E=0 }
    const double  r2  = dp.norm2();
    double E,F;
    // ---- Electrostatic
    if( r2 > Cr2_cut ){
        const double ir2_ = 1/( r2 + R2damp  );
        E =  COULOMB_CONST*REQH.z*sqrt( ir2_ );
        F =  E*ir2_ ;
    }else{ E=0; F=0; };
    // --- LJ 
    if( r2 > LJr2_cut  ){
        const double  ir2 = 1/r2;
        const double  u2  = REQH.x*REQH.x*ir2;
        const double  u6  = u2*u2*u2;
        const double vdW  = u6*REQH.y;
        const double   H  = u6*( (REQH.w<0) ? REQH.w*REQH.y : 0.0 );  // H-bond correction
        E   +=       (u6-2.)*vdW            ;
        F   +=  (12.*(u6-1.)*vdW + H*6.)*ir2;
    }
    f.set_mul( dp, -F );
    return E;
}

// evaluate energy and force using Morse and Coulomb potential and Hydrogen bond pseudo-charges, with cutoffs
inline double getMorseQH_cut( const Vec3d& dp, Vec3d& f, const Quat4d& REQH, const double K, double R2damp, const double Cr2_cut, const double Mr2_cut ){
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

// evaluate energy and force using Morse 
inline void addAtomicForceMorse( const Vec3d& dp, Vec3d& f, double r0, double eps, double beta ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    const double R2ELEC = 1.0;
    double r     = sqrt( dp.norm2()+R2SAFE );
    double expar = exp ( beta*(r-r0) );
    //double E     = eps*( expar*expar - 2*expar );
    double fr    = eps*2*beta*( expar*expar - expar );
    //printf( " %g -> %g | (%g,%g,%g) %g\n" , r, fr,  r0, eps,  q, alpha );
    //printf( " r %g expar %g fr %g kqq %g a %g eps %g \n" , r, expar, fr, COULOMB_CONST*qq, alpha, eps );
    f.add_mul( dp, fr/r );
}

// evaluate energy and force using Morse and Coulomb potential
inline double addAtomicForceMorseQ( const Vec3d& dp, Vec3d& f, double r0, double E0, double qq, double K=-1., double R2damp=1. ){
    double r2    = dp.norm2();
    double ir2_  = 1/(r2+R2damp);
    double r     = sqrt( r2   + 1e-32 );
    double ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    double e     = exp( K*(r-r0));
    //double e2    = e*e;
    //double fMors =  E0*  2*K*( e2 -   e ); // Morse
    //double EMors =  E0*      ( e2 - 2*e );
    float   Ae  =  E0*e;
    double fMors =  Ae*  2*K*(e - 1); // Morse
    double EMors =  Ae*      (e - 2);
    //if(idebug>0) printf("CPU expar %g ", e  );
    //printf(  "addAtomicForceMorseQ() k %g r %g e %g E0 %g E %g \n", K, r, e/exp( K*(1.487)), E0, EMors );
    //printf(  "addAtomicForceMorseQ() k %g r %g e %g E0 %g E %g \n", K, r, e, E0, EMors );
    //fr          += COULOMB_CONST*qq/( r*r + R2ELEC );   // Comlomb cheal_damp : Problem - it would reqire asinh() to get energy
    double Eel   = COULOMB_CONST*qq*ir_;
    f.add_mul( dp, fMors/r - Eel*ir2_ );
    //printf( "r %g E0 %g E %g  e2 %g -2*e %g  \n ", r, E0, EMors, e2, -2*e );
    return EMors + Eel;
}

// evaluate energy and force using Harmonic potential
inline double addAtomicForceQ_R2( const Vec3d& dp, Vec3d& f, double qq, double K=-1., double R2damp=1. ){
    double r2    = dp.norm2();
    double ir2_  = 1/(r2+R2damp);
    double ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    double Eel   = COULOMB_CONST*qq*ir_;
    f.add_mul( dp, -Eel*ir2_ );
    return Eel;
}

// evaluate energy and force using Coulomb potential
inline double addAtomicForceQ( const Vec3d& dp, Vec3d& f, double qq ){
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double E    = COULOMB_CONST*qq*ir;
    double fr   = -E*ir2;
    f.add_mul( dp, fr );
    return E;
}

// evaluate energy and force using Lennard-Jones potential
inline void addAtomicForceLJ( const Vec3d& dp, Vec3d& f, double r0, double eps ){
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps )*ir2;
    f.add_mul( dp, fr );
}

// evaluate energy and force using Exponential potential
inline void addAtomicForceExp( const Vec3d& dp, Vec3d& f, double r0, double eps, double alpha ){
    double r    = sqrt(dp.norm2() + R2SAFE );
    double E    = eps*exp( alpha*(r-r0) );
    double fr   = alpha*E/r;
    f.add_mul( dp, fr );
    //f.add_mul( dp, 1/(dp.norm2()+R2SAFE) ); // WARRNING DEBUG !!!!
}

// transform from non-covalent parameters from (R_vdw,E_vdw,Q) to (Pauli,London,Q), single precision
inline Quat4f REQ2PLQ( const Quat4d& REQ, double K ){
    float e  = (float) exp(K*REQ.x);
    float cL = (float) e*REQ.y;
    float cP = (float) e*cL;
    float cH = (float) e*e*REQ.w;
    return Quat4f{ cP, cL, REQ.z, cH };
}

// transform from non-covalent parameters from (R_vdw,E_vdw,Q) to (Pauli,London,Q), double precision
inline Quat4d REQ2PLQ_d( const Quat4d& REQ, double K ){
    double e  = exp(K*REQ.x);
    double cL = e*REQ.y;
    double cP = e*cL;
    double cH = e*e*REQ.w;
    return Quat4d{ cP, cL, REQ.z, cH };
}

// transform from non-covalent parameters from (Pauli,London,Q) to (R_vdw,E_vdw,Q),  double precision
inline Quat4f REnergyQ2PLQ( const Quat4d& REQ, double alpha ){
    return REQ2PLQ( Quat4d{REQ.x, sqrt(REQ.y), REQ.z, REQ.w}, alpha );
}
// inline Quat4d REnergyQ2PLQ( const Quat4d& REQ, double alpha ){
//     return REQ2PLQ( Quat4d{REQ.x, sqrt(REQ.y), REQ.z, REQ.w}, alpha );
// }

// ================= Force Bounding Box, Plane etc.

// evaluate force from bounding box forces in 1D
inline double boxForce1D(double x, double xmin, double xmax, double k){
    double f=0;
    if(k<0) return 0;
    if(x>xmax){ f+=k*(xmax-x); }
    if(x<xmin){ f+=k*(xmin-x); }
    return f;
}

// evaluate force from bounding box forces in 3D
inline void boxForce(const Vec3d& p, Vec3d& f,const Vec3d& pmin, const Vec3d& pmax, const Vec3d& k){
    f.x+=boxForce1D( p.x, pmin.x, pmax.x, k.x);
    f.y+=boxForce1D( p.y, pmin.y, pmax.y, k.y);
    f.z+=boxForce1D( p.z, pmin.z, pmax.z, k.z);
}

// evaluate spring force in given direction (normal) from plane and point c0
inline Vec3d getForceSpringPlane( const Vec3d& p, const Vec3d& normal, double c0, double k ){
    double cdot = normal.dot(p) - c0;
    return normal * (cdot * k);
}

// ================= Force from Surface & Plan

// evaluate force from plane using Hamaker potential
inline Vec3d getForceHamakerPlane( const Vec3d& p, const Vec3d& normal, double z0, double amp, double R ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    //printf(  " normal %g %g %g \n", normal.x, normal.y, normal.z );
    double cdot = normal.dot(p) - z0;
    double ir   = R/cdot;
    double ir3  = ir*ir*ir;
    double f    = amp*(ir/R)*ir3*(ir3-1);
    //printf( "%g %g %g %g %g %g %g \n", f, cdot, ir, ir3, e0, c0, r0  );
    return normal * f;
}

// evaluate force from plane using Morse potential
inline Vec3d getForceMorsePlane( const Vec3d& p, const Vec3d& normal, double amp, double R, double beta ){
    // https://en.wikipedia.org/wiki/Lennard-Jones_potential
    //printf(  " normal %g %g %g \n", normal.x, normal.y, normal.z );
    double r       = normal.dot(p) - R;
    double expar   = exp( beta*r );
    double m_expar = 1-expar;
    double E       =  amp*m_expar*m_expar;
    double f       = -amp*m_expar*  expar * beta;
    //printf( "%g %g %g %g %g %g %g \n", f, cdot, ir, ir3, e0, c0, r0  );
    return normal * f;
}

// ================= Simple Radial Forces & Pulling

// evaluate spring force on particle p, from line (ray) defined by point ray0 and hray (normalized direction vector)
inline Vec3d getForceSpringRay( const Vec3d& p, const Vec3d& hray, const Vec3d& ray0, double k ){
    Vec3d dp; dp.set_sub( p, ray0 );
    double cdot = hray.dot(dp);
    dp.add_mul(hray,-cdot);
    return dp*k;
}

// evaluate spring force on particl
inline double addForceR2( const Vec3d& dp, Vec3d& f, double R2, double K ){
    double r2 = dp.norm2();
    if(r2<R2){
        double u  = R2-r2;
        f.add_mul( dp, 4*K*u );
        return K*u*u;
    }
    return 0;
}

inline double addForceR2inv( const Vec3d& dp, Vec3d& f, double R2, double K, double w2 ){
    //  E =   K*(1-R2/r2)^2
    //  f = 4*K*(1-R2/r2)*(R2/(r2*r2))*x
    //  k = 4*    5R4/r6 - 3R2/r4
    //  k(r=R) 4*(5      - 3 )/r2 = 8/R2
    //  => K_ = K*R2/8
    double r2 = dp.norm2();
    if(r2<R2){
        double R2_ = R2+w2;
        double K_  = R2_*0.125;
        double ir2 = 1/(r2+w2);
        double u   = R2_*ir2-1;
        f.add_mul( dp, K_*u*ir2 );
        return K_*u*u;
    }
    return 0;
}

inline double addForceR2mix( const Vec3d& dp, Vec3d& f, double R2, double K, double w2 ){
    //  E =   K*(R2-r2)(1-R2/r2)
    //  f = 2*K*(1-(R2/r2)^2)*x
    //  k = 2*K*(3*R4/r4-1)
    //  k(r=R) 4*K
    double r2 = dp.norm2();
    if(r2<R2){
        double K_  = K*0.5;
        double R2_ = R2+w2;
        double u   = R2_/(r2+w2);
        f.add_mul( dp, K_*(1-u*u) );
        return K_*u*(R2-r2)*0.5;
    }
    return 0;
}

// =============== Force Checking

inline void sum(int n, Vec3d* ps, Vec3d& psum){ for(int i=0;i<n;i++){ psum.add(ps[i]); } };

inline void sumTroq(int n, Vec3d* fs, Vec3d* ps, const Vec3d& cog, const Vec3d& fav, Vec3d& torq){
    for(int i=0;i<n;i++){  torq.add_cross(ps[i]-cog,fs[i]-fav);  }
    //for(int i=0;i<n;i++){  torq.add_cross(ps[i],fs[i]);  }
}

inline void checkForceInvariatns( int n, Vec3d* fs, Vec3d* ps, Vec3d& cog, Vec3d& fsum, Vec3d& torq ){
    cog =Vec3dZero;
    fsum=Vec3dZero;
    torq=Vec3dZero;
    double dw = 1./n;
    sum(n, ps, cog ); cog.mul(dw);
    sum(n, fs, fsum); //cog.mul(dw);
    sumTroq(n, fs, ps, cog, fsum*dw, torq );
}

#endif
