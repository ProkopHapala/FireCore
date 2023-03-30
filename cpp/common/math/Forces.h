
#ifndef Forces_h
#define Forces_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

#define COULOMB_CONST         14.3996448915

#define RSAFE   1.0e-4f
#define R2SAFE  1.0e-8f
#define F2MAX   10.0f

// ================ Angular Forces (MMFF) 

inline double evalBond( const Vec3d& h, double dl, double k, Vec3d& f ){
    double fr = dl*k;
    f.set_mul ( h, fr );
    //f1.add( h );
    //f2.sub( h );
    return fr*dl*0.5;
}

inline double spring( double l, Vec2d ls, Vec2d ks, double flim, double& f ){
    double E=0;
    f       =0;
    if  (l>ls.x){
        double dl=l-ls.x;
        f=dl*ks.x;
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
    double c2 = h.norm2()*0.25;               // cos(a/2) = |ha+hb|
    double s2 = 1-c2;
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
    f1.set_lincomb( fr1, h,  fr1*c2, h1 );  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    f2.set_lincomb( fr2, h,  fr2*c2, h2 );  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );
    //printf( "evalAngleCosHalf fr=%g cs2(%g,%g) cso(%g,%g) cs(%g,%g) cs0(%g,%g) \n", fr, c2,s2,  c,s,  cs.x,cs.y, cs0.x,cs0.y );

    return E;
}
/*
inline double evalAngleCosHalf( const float4 hr1, const float4 hr2, const float2 cs0, double k, __private float3* f1, __private float3* f2 ){
    // This is much better angular function than evalAngleCos() with just a little higher computational cost ( 2x sqrt )
    float3 h  = hr1.xyz + hr2.xyz;
    float  c2 = dot(h,h)*0.25f;              // cos(a/2) = |ha+hb|
    float  s2 = 1.f-c2;
    float2 cso = (float2){ sqrt(c2), sqrt(s2) };
    float2 cs = udiv_cmplx( cs0, cso );
    float  E         =  k*( 1.f - cs.x );  // just for debug ?
    float  fr        = -k*(       cs.y );
    c2 *= -2.f;
    fr /=  4.f*cso.x*cso.y;   //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    float  fr1    = fr*hr1.w;
    float  fr2    = fr*hr2.w;
    *f1 =  h*fr1  + hr1.xyz*(fr1*c2);  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    *f2 =  h*fr2  + hr2.xyz*(fr2*c2);  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );
    return E;
}
*/
inline double evalAngleCos( const Vec3d& h1, const Vec3d& h2, double ir1, double ir2, double K, double c0, Vec3d& f1, Vec3d& f2 ){
    double c = h1.dot(h2);
    //f1 = h2 - h1*c;
    //f2 = h1 - h2*c;
    f1.set_add_mul( h2,h1,-c );
    f2.set_add_mul( h1,h2,-c );
    double c_   = c-c0;
    double E    =  K*c_*c_;
    double fang = -K*c_*2;
    f1.mul( fang*ir1 );
    f2.mul( fang*ir2 );
    return E;
}

inline double evalPiAling( const Vec3d& h1, const Vec3d& h2, double ir1, double ir2, double K, Vec3d& f1, Vec3d& f2 ){  // interaction between two pi-bonds
    double c = h1.dot(h2);
    //f1 = h2 - h1*c;
    //f2 = h1 - h2*c;
    f1.set_add_mul( h2,h1,-c );
    f2.set_add_mul( h1,h2,-c );
    bool sign = c<0; if(sign) c=-c;
    double E    = -K*c;
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

inline void combineREQ(const Vec3d& a, const Vec3d& b, Vec3d& out){
    out.a=a.a+b.a; // radius
    out.b=a.b*b.b; // epsilon
    out.c=a.c*b.c; // q*q
}

inline void addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, double r0, double eps, double qq ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*r0*r0;
    double ir6  = ir2_*ir2_*ir2_;
    double fr   = ( ( 1 - ir6 )*ir6*12*eps + ir*qq*-COULOMB_CONST )*ir2;
    f.add_mul( dp, fr );
}

inline double addAtomicForceLJQ( const Vec3d& dp, Vec3d& f, const Vec3d& REQ ){
    //Vec3f dp; dp.set_sub( p2, p1 );
    double ir2  = 1/( dp.norm2() + 1e-4 );
    double ir   = sqrt(ir2);
    double ir2_ = ir2*REQ.a*REQ.a;
    double ir6  = ir2_*ir2_*ir2_;
    //double fr   = ( ( 1 - ir6 )*ir6*12*REQ.b + ir*REQ.c*-COULOMB_CONST )*ir2;
    double Eel  = ir*REQ.c*COULOMB_CONST;
    double vdW  = ir6*REQ.b;
    double fr   = ( ( 1 - ir6 )*12*vdW - Eel )*ir2;
    //printf( " (%g,%g,%g) r %g fr %g \n", dp.x,dp.y,dp.z, 1/ir, fr );
    f.add_mul( dp, fr );
    return  ( ir6 - 2 )*vdW + Eel;
}

inline float getLJQ( Vec3f dp, Vec3f REQ, float R2damp, Vec3f& f ){
    // ---- Electrostatic
    float   r2   = dp.norm2();
    float   ir2_ = 1/( r2 + R2damp  );
    float   Ec   = COULOMB_CONST*REQ.z*sqrt( ir2_ );
    // --- LJ 
    float  ir2 = 1/r2;
    float  u2  = REQ.x*REQ.x*ir2;
    float  u6  = u2*u2*u2;
    float vdW  = u6*REQ.y;
    float E    =      (u6-2.)*vdW     + Ec      ;
    float fr   = -12.*(u6-1.)*vdW*ir2 - Ec*ir2_ ;
    f.set_mul( dp, fr );
    return E;
}

inline double getLJQ( Vec3d dp, Vec3d REQ, double R2damp, Vec3d& f ){
    // ---- Electrostatic
    double   r2   = dp.norm2();
    double   ir2_ = 1/( r2 + R2damp  );
    double   Ec   = COULOMB_CONST*REQ.z*sqrt( ir2_ );
    // --- LJ 
    double  ir2 = 1/r2;
    double  u2  = REQ.x*REQ.x*ir2;
    double  u6  = u2*u2*u2;
    double vdW  = u6*REQ.y;
    double E    =      (u6-2.)*vdW     + Ec      ;
    double fr   = -12.*(u6-1.)*vdW*ir2 - Ec*ir2_ ;
    f.set_mul( dp, fr );
    return E;
}

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

inline double addAtomicForceMorseQ( const Vec3d& dp, Vec3d& f, double r0, double E0, double qq, double K=-1., double R2damp=1. ){
    double r2    = dp.norm2();
    double ir2_  = 1/(r2+R2damp);
    double r     = sqrt( r2   );
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
    //fr          += COULOMB_CONST*qq/( r*r + R2ELEC );   // Comlomb cheal_damp : Problem - it would reqire asinh() to get energy
    double Eel   = COULOMB_CONST*qq*ir_;
    f.add_mul( dp, fMors/r - Eel*ir2_ );
    //printf( "r %g E0 %g E %g  e2 %g -2*e %g  \n ", r, E0, EMors, e2, -2*e );
    return EMors + Eel;
}

inline double addAtomicForceQ_R2( const Vec3d& dp, Vec3d& f, double qq, double K=-1., double R2damp=1. ){
    double r2    = dp.norm2();
    double ir2_  = 1/(r2+R2damp);
    double ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    double Eel   = COULOMB_CONST*qq*ir_;
    f.add_mul( dp, -Eel*ir2_ );
    return Eel;
}

inline double addAtomicForceQ( const Vec3d& dp, Vec3d& f, double qq ){
    double ir2  = 1/( dp.norm2() + R2SAFE );
    double ir   = sqrt(ir2);
    double E    = COULOMB_CONST*qq*ir;
    double fr   = -E*ir2;
    f.add_mul( dp, fr );
    return E;
}

inline void addAtomicForceLJ( const Vec3d& dp, Vec3d& f, double r0, double eps ){
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
    //f.add_mul( dp, 1/(dp.norm2()+R2SAFE) ); // WARRNING DEBUG !!!!
}

inline Vec3d REQ2PLQ( const Vec3d& REQ, double K ){
    double e   = exp(K*REQ.x);
    double cL  = e*REQ.y;
    double cP  = e*cL;
    return Vec3d{ cP, cL, REQ.z };
}

inline Vec3d REnergyQ2PLQ( const Vec3d& REQ, double alpha ){
    return REQ2PLQ( {REQ.x, sqrt(REQ.y), REQ.z}, alpha );
}

// ================= Force Bounding Box, Plane etc.

inline double boxForce1D(double x, double xmin, double xmax, double k){
    double f=0;
    if(k<0) return 0;
    if(x>xmax){ f+=k*(xmax-x); }
    if(x<xmin){ f+=k*(xmin-x); }
    return f;
}

inline void boxForce(const Vec3d& p, Vec3d& f,const Vec3d& pmin, const Vec3d& pmax, const Vec3d& k){
    f.x+=boxForce1D( p.x, pmin.x, pmax.x, k.x);
    f.y+=boxForce1D( p.y, pmin.y, pmax.y, k.y);
    f.z+=boxForce1D( p.z, pmin.z, pmax.z, k.z);
}

inline Vec3d getForceSpringPlane( const Vec3d& p, const Vec3d& normal, double c0, double k ){
    double cdot = normal.dot(p) - c0;
    return normal * (cdot * k);
}

// ================= Force from Surface & Plan

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

inline Vec3d getForceSpringRay( const Vec3d& p, const Vec3d& hray, const Vec3d& ray0, double k ){
    Vec3d dp; dp.set_sub( p, ray0 );
    double cdot = hray.dot(dp);
    dp.add_mul(hray,-cdot);
    return dp*k;
}

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
