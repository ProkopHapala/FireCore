
#ifndef  Multipoles_h
#define  Multipoles_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"
#include "constants.h"

// https://en.wikipedia.org/wiki/Multipole_expansion
// https://en.wikipedia.org/wiki/Axial_multipole_moments
// https://en.wikipedia.org/wiki/Spherical_multipole_moments
// http://physics.stackexchange.com/questions/106564/clarification-of-multipole-expansion-for-a-point-charge


// fast multipole method
// https://github.com/davidson16807/fast-multipole-method/blob/master/fast-multipole-method.js
// https://github.com/barbagroup/pyfmm/tree/master/pyFMM

/*

From Maxima:

# ------- Monopole
Em_0: 1/r

Fm_x = Epx  :  x/r^3 = xh/r^2
Fm_y = Epy  :  y/r^3 = yh/r^2
Fm_z = Epz  :  z/r^3 = zh/r^2

# ------- Dipole
Epx = Fm_x  :  x/r^3 = xh/r^2
Epy = Fm_y  :  y/r^3 = yh/r^2
Epz = Fm_z  :  z/r^3 = zh/r^2

Fpx_x = Eqxx =        :  ( 2x^2 -( z^2 + y^2 ) )/r^5   =  (3x^2 - r^2)/r^5
Fpx_y = Fpy_x = Eqxy  :  ( 3xy )/r^5
Fpx_z = Fpz_x = Eqxz  :  ( 3xz )/r^5
Fpy_y = Eqyy =        :  ( 2y^2 -( z^2 + x^2 ) )/r^5   =  (3y^2 - r^2)/r^5
Fpy_x = Fpx_y = Eqxy  :  ( 3xy )/r^5
Fpy_z = Fpz_y = Eqyz  :  ( 3yz )/r^5
Fpz_z = Eqzz =        :  ( 2z^2 -( y^2 + x^2 ) )/r^5   =  (3z^2 - r^2)/r^5
Fpz_x = Fpx_z = Eqxz  :  ( 3xz )/r^5
Fpz_y = Fpy_z = Eqyz  :  ( 3yz )/r^5

# ------- Quadrupole

Eqxx = Fpx_x          :  ( 2x^2 -( z^2 + y^2 ) )/r^5   =  (3x^2 - r^2)/r^5
Eqyx = Fpx_y = Fpy_x  :  ( 3xy )/r^5

(3z^2 + 3y^2 - 2*x^2 ) = 3(z^2 + y^2 + x^2 ) - 5*r^2
( z^2 - 4*y^2 + x^2 )  = r^2 - 5*y^2 


Fqxx_x :      3x(3z^2 + 3y^2 - 2*x^2 )/r^7    =  3x(3r^2 - 5*x^2 )/r^7
Fqxy_y :      3x( z^2 - 4*y^2 + x^2 )         =  3x( r^2 - 5*y^2 )/r^7
Fqxy_z :      3y( z^2 - 4*y^2 + x^2 )         =  15xyz/r^7

*/


namespace Multiplole{

__attribute__((hot)) 
inline double Coulomb( const Vec3d& d, double qq, Vec3d* f=0 ){
    qq *= COULOMB_CONST;
    const double r2  = d.norm2();
    const double ir2 = 1/r2;
    const double ir  = 1/sqrt(r2);
    if( f ) *f = d*qq*ir2*ir;
    return qq*ir;
}

__attribute__((hot)) 
double EF_brute( const Vec3d& p, double Q, Vec3d* f,  int n, Vec3d * ps, double * Qs ){
    double E = 0;
    for( int i=0; i<n; i++){
        E += Coulomb( ps[i]-p, Q*Qs[i], f );
    }
    return E;
}

__attribute__((hot)) 
inline void project( const Vec3d& d, const double Q, const int order, double * cs ){
    cs[0] += Q;
    if(order<1) return;
    cs[1] += d.x*Q;
    cs[2] += d.y*Q;
    cs[3] += d.z*Q;
    if(order<2) return;
    cs[4] += d.x*d.x*Q; // xx
    cs[5] += d.y*d.y*Q; // yy
    cs[6] += d.z*d.z*Q; // zz
    cs[7] += d.y*d.z*Q; // yz
    cs[8] += d.x*d.z*Q; // xz
    cs[9] += d.x*d.y*Q; // xy
}


__attribute__((hot)) 
Vec3d center( int n, const Vec3d * ps, const double * Qs ){
    Vec3d c = Vec3dZero;
    double Q = 0;
    for( int i=0; i<n; i++){
        const double qi = Qs[i];
        c.add_mul( ps[i], qi );
        Q += qi;
    }
    c.mul( 1/Q );
    return c;
}

__attribute__((hot)) 
void project( const Vec3d* p0_, int n, const Vec3d * ps, const double * Qs, int order, double * cs, bool bClear=true ){
    Vec3d p0;
    if(p0_){ p0=*p0_; }else{ p0=center(n, ps, Qs); }
    if(bClear) for( int i=0; i<10; i++ ) cs[i]=0;
    for( int i=0; i<n; i++){
        project( ps[i]-center, Qs[i], order, cs );
    }
}

double Emultipole( const Vec3d& d, int order, double * cs ){
    //double r   = dR.norm();
    //double ir  = 1 / r;
    //double ir2 = ir*ir;
    double ir2 = 1/d.norm2();
    double E   = cs[0];
    if( order>0 ) E += ir2    *( cs[1]*d.x + cs[2]*d.y + cs[3]*d.z );
    if( order>1 ) E += ir2*ir2*((cs[4]*d.x + cs[5]*d.y)*d.x +
                                (cs[6]*d.y + cs[7]*d.z)*d.y +
                                (cs[8]*d.z + cs[9]*d.x)*d.z );
    return sqrt(ir2)*E;
}

/*

Fpx_x = Eqxx =        :  ( 2x^2 -( z^2 + y^2 ) )/r^5   =  (3x^2 - r^2)/r^5
Fpx_y = Fpy_x = Eqxy  :  ( 3xy )/r^5
Fpx_z = Fpz_x = Eqxz  :  ( 3xz )/r^5

Fqxx_x :      3x(3z^2 + 3y^2 - 2*x^2 )/r^7    =  3x(3r^2 - 5*x^2 )/r^7
Fqxx_z :      3z(z^2  + y^2  - 4*x^2 )        =  3z( r^2 - 5*x^2 )/r^7
Fqxy_y :      3x( z^2 - 4*y^2 + x^2 )         =  3x( r^2 - 5*y^2 )/r^7
Fqxy_z :      3y( z^2 - 4*y^2 + x^2 )         =  15xyz/r^7


*/
double EFmultipole( const Vec3d& d, Vec3d& f, double * cs, int order=2 ){
    //double r   = dR.norm();
    //double ir  = 1 / r;
    //double ir2 = ir*ir;
    const double r2  =  d.norm2();
    const double ir2 = 1/r2;
    const double ir  = sqrt(ir2); 
    // ---------- Monopole
    const double  Em = cs[0]*ir;
    double E  = Em;
    f         = d*(Em*ir2);
    if( order<1 ) return E;

    // ---------- Dipole
    const Vec3d p   = *(Vec3d*)(cs+1);
    const double pd =  p.dot(d);
    E += ir2 * pd;
    const double scp = ir2*ir2*ir;
    const Vec3d Fp{
        //   cc               
        // Question: are the signs correct? It is strange it can be reduced to dot-product
        //(( p.x*d.x   +    p.y*d.y + p.z*d.z   )*3-r2)*d.x*scp,
        //(( p.y*d.y   +    p.x*d.x + p.z*d.z   )*3-r2)*d.y*scp,
        //(( p.z*d.z   +    p.x*d.x + p.y*d.y   )*3-r2)*d.z*scp
        ( pd*3-r2)*d.x*scp,
        ( pd*3-r2)*d.y*scp,
        ( pd*3-r2)*d.z*scp
    };
    f.add(Fp);
    if( order<2 ) return E*ir;

    // ---------- Quadrupole
    //E += ir2*ir2*ir2();
    Vec6d q  = *(Vec6d*)(cs+4);
    Vec6d dd{ d.x*d.x, d.y*d.y, d.z*d.z, d.y*d.z, d.x*d.z, d.x*d.y };
    const Vec3d Fp{
        //   cc_c           ab_b  ab_a               
        (( q.xx*dd.xx*5   +    (q.yy)*dd.xy + q.zz*d.z  )*3 - 9*r2 )*d.x*scp,
        (( q.yy*dd.yy*5   +    *dd.xz + cs[3]*d.z       )*3 - 9*r2 )*d.y*scp,
        (( q.zz*dd.zz*5   +    *dd.x + cs[2]*d.y        )*3 - 9*r2 )*d.z*scp
    };
    f.add(Fp);
    return E;
}

}

#endif

