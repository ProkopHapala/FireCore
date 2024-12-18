
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
Vec3d project( const Vec3d* p0_, int n, const Vec3d * ps, const double * Qs, int order, double * cs, bool bClear=true ){
    Vec3d p0;
    if(p0_){ p0=*p0_; }else{ p0=center(n, ps, Qs); }
    if(bClear) for( int i=0; i<10; i++ ) cs[i]=0;
    for( int i=0; i<n; i++){
        project( ps[i]-p0, Qs[i], order, cs );
        //printf( "project[ia=%3i] p(%+10.5e,%+10.5e,%+10.5e) Q=%+10.5e \n", i, ps[i].x,ps[i].y,ps[i].z,  Qs[i] );
        //printf( "project[ia=%3i] Q=%+10.5e p(%+10.5e,%+10.5e,%+10.5e) Qxx,yy,zz(%+10.5e,%+10.5e,%+10.5e)|yz,xz,xy(%+10.5e,%+10.5e,%+10.5e)\n", i, cs[0], cs[1],cs[2],cs[3],  cs[4],cs[5],cs[6], cs[7],cs[8],cs[9] );
    }
    return p0;
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

Fp_x   :	-(   Px*z^2  + Px*y^2 -2*Px*x^2  -3*Py*x*y  -3*Pz*x*z  )/r^5 =   ( Px*( r^2 - 3*x^2 ) - 3*(Py*x*y + Pz*x*z) )/r^5
Fp_y   :	-(   Py*z^2  + Py*x^2 -2*Py*y^2  -3*Pz*y*z  -3*Px*x*y  )/r^5 =   ( Py*( r^2 - 3*y^2 ) - 3*(Px*x*y + Pz*y*z) )/r^5
Fp_z   :	-(   Pz*y^2  + Pz*x^2 -2*Pz*z^2  -3*Py*y*z  -3*Px*x*z  )/r^5 =   ( Pz*( r^2 - 3*z^2 ) - 3*(Px*x*z + Py*y*z) )/r^5


Fq_x  : -3*(Q_xz*z^3+Q_xy*y*z^2-4*Q_zz*x*z^2+Q_yy*x*z^2+3*Q_xx*x*z^2+Q_xz*y^2*z-5*Q_yz*x*y*z-4*Q_xz*x^2*z+Q_xy*y^3+Q_zz*x*y^2-4*Q_yy*x*y^2+3*Q_xx*x*y^2-4*Q_xy*x^2*y+Q_zz*x^3+Q_yy*x^3-2*Q_xx*x^3)
Fq_y  : -3*(Q_yz*z^3-4*Q_zz*y*z^2+3*Q_yy*y*z^2+Q_xx*y*z^2+Q_xy*x*z^2-4*Q_yz*y^2*z-5*Q_xz*x*y*z+Q_yz*x^2*z+Q_zz*y^3-2*Q_yy*y^3+Q_xx*y^3-4*Q_xy*x*y^2+Q_zz*x^2*y+3*Q_yy*x^2*y-4*Q_xx*x^2*y+Q_xy*x^3)
Fq_z  :  3*(Q_yy*z^3+Q_xx*z^3-4*Q_zz*y*z^2-4*Q_yz*y*z^2-4*Q_yy*y^2*z+Q_xx*y^2*z-5*Q_xz*x*y*z-5*Q_xy*x*y*z+Q_yy*x^2*z-4*Q_xx*x^2*z+Q_zz*y^3+Q_yz*y^3+Q_zz*x^2*y+Q_yz*x^2*y)


Fq_x  : -3*(  -2*Q_xx*x^3  +3*Q_xx*x*z^2    +3*Q_xx*x*y^2            +Q_yy*x^3 +Q_yy*x*z^2 -4*Q_yy*x*y^2             +Q_zz*x^3  +Q_zz*x*y^2  -4*Q_zz*x*z^2                     +Q_xy*y^3  +Q_xy*y*z^2  -4*Q_xy*x^2*y                    +Q_xz*z^3   +Q_xz*y^2*z  -4*Q_xz*x^2*z       -5*Q_yz*x*y*z        )
Fq_x  : -3*(   Q_xx*(-2*x^3 + 3*x*z^2 + 3*x*y^2)                     +Q_yy*(x^3 + x*z^2 - 4*x*y^2)                   +Q_zz*(x^3  +x*y^2  -4*x*z^2)                             +Q_xy*(y^3  +y*z^2  -4*x^2*y)                            +Q_xz*(z^3   +z*y^2  -4*z*x^2)               -5*Q_yz*(x*y*z)      )
Fq_x  : -3*(   Q_xx*(-5*x^2 + 3*r^2 )*x                              +Q_yy*(r^2 - 5*y^2)*x                           +Q_zz*(r^2  -5*z^2)*x                                     +Q_xy*(r^2  -5*x^2)*y                                    +Q_xz*(r^2  -5*x^2)*z                        -5*Q_yz*(x*y*z)      )
Fq_x  : -3*(   r^2*(  (3*Q_xx + Q_yy + Q_zz)*x   + (Q_xy*y+Q_xz*z)  )    -5*x*(Q_xx*x^2 + Q_yy*y^2 +  Q_zz*z^2)                                                                -5*x^2*(   Q_xy*y + Q_xz*z )                                                                          -5*Q_yz*(x*y*z)      )


Fq_y  : -3*(  -2*Q_yy*y^3  +3*Q_yy*x^2*y   +3*Q_yy*y*z^2            +Q_xx*y^3       +Q_xx*y*z^2   -4*Q_xx*x^2*y      +Q_zz*x^2*y    +Q_zz*y^3     -4*Q_zz*y*z^2            +Q_xy*x*z^2     +Q_xy*x^3   -4*Q_xy*x*y^2                   +Q_yz*x^2*z      +Q_yz*z^3     -4*Q_yz*y^2*z                                          -5*Q_xz*x*y*z          )
Fq_y  : -3*(  Q_yy*( -2*y^3  +3*( x^2*y   +y*z^2)                   +Q_xx*( y^3     +y*z^2   -4*x^2*y)                +Q_zz*( x^2*y       +y^3     -4y*z^2      )           +Q_xy*( x*z^2   +x^3        -4*x*y^2  )                     +Q_yz*( x^2*z      + z^3       -4* y^2*z     )                                        -5*Q_xz*(x*y*z)        )       
Fq_y  : -3*(  Q_yy*( -3*y^2 + r^2  )*y                              +Q_xx*( r^2   -5*x^2)*y                           +Q_zz*( r^2                  -5z^2        )*y         +Q_xy*( r^2                 -5*y^2    )*x                   +Q_yz*( r^2                    -5* y^2       )*z                                      -5*Q_xz*(x*y*z)        )     


Fq_x  : -3*(   r^2*( (3*Q_xx + Q_yy + Q_zz)*x  + (Q_xy*y+Q_xz*z)  )    -5*x*(Q_xx*x^2 + Q_yy*y^2 + Q_zz*z^2)   -5*x^2*( Q_xy*y + Q_xz*z )     -5*Q_yz*(x*y*z)  )
Fq_y  : -3*(   r^2*( (3*Q_yy + Q_xx + Q_zz)*y  + (Q_xy*x+Q_yz*z)  )    -5*y*(Q_xx*x^2 + Q_yy*y^2 + Q_zz*z^2)   -5*y^2*( Q_xy*x + Q_yz*z )     -5*Q_xz*(x*y*z)  )
Fq_z  : -3*(   r^2*( (3*Q_zz + Q_xx + Q_yy)*z  + (Q_xz*x+Q_yz*y)  )    -5*z*(Q_xx*x^2 + Q_yy*y^2 + Q_zz*z^2)   -5*z^2*( Q_xz*x + Q_yz*y )     -5*Q_xy*(x*y*z)  )


Qd2 = (Q_xx*x^2 + Q_yy*y^2 + Q_zz*z^2)
Qsum

Fq_x  : -3*(   r^2*( (3*Q_xx + Q_yy + Q_zz)*x  + (Q_xy*y+Q_xz*z)  )    -5*( x*Qd2    +x^2*( Q_xy*y + Q_xz*z )     +Q_yz*(x*y*z) ) )
Fq_y  : -3*(   r^2*( (3*Q_yy + Q_xx + Q_zz)*y  + (Q_xy*x+Q_yz*z)  )    -5*( y*Qd2    +y^2*( Q_xy*x + Q_yz*z )     +Q_xz*(x*y*z) ) )
Fq_z  : -3*(   r^2*( (3*Q_zz + Q_xx + Q_yy)*z  + (Q_xz*x+Q_yz*y)  )    -5*( z*Qd2    +z^2*( Q_xz*x + Q_yz*y )     +Q_xy*(x*y*z) ) )


Fq_x  : -3*(   r^2*(3*Q_xx + Q_yy + Q_zz)*x      -5*x*Qd2 + (r^2-5*x^2)*( Q_xy*y + Q_xz*z )     +Q_yz*(x*y*z) ) )
Fq_y  : -3*(   r^2*(3*Q_yy + Q_xx + Q_zz)*y      -5*y*Qd2 + (r^2-5*y^2)*( Q_xy*x + Q_yz*z )     +Q_xz*(x*y*z) ) )
Fq_z  : -3*(   r^2*(3*Q_zz + Q_xx + Q_yy)*z      -5*z*Qd2 + (r^2-5*z^2)*( Q_xz*x + Q_yz*y )     +Q_xy*(x*y*z) ) )



 r^2*(Q_xy*y+Q_xz*z)    -5*x^2*( Q_xy*y + Q_xz*z )  =   (r^2-5*x^2)*( Q_xy*y + Q_xz*z )       =   Q_xy*y*(  r^2 -5*x^2  ) + Q_xz*z*( r^2 -5*x^2  )


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
    // const Vec3d Fp{
    //     //   cc               
    //     // Question: are the signs correct? It is strange it can be reduced to dot-product
    //     //(p.x*(r2-3*dd.xx) - 3*(p.y*dd.xy + p.z*dd.xz) )*scp,
    //     //(p.y*(r2-3*dd.yy) - 3*(p.x*dd.xy + p.z*dd.yz) )*scp,
    //     //(p.z*(r2-3*dd.zz) - 3*(p.x*dd.xz + p.y*dd.yz) )*scp,

    //     //( p.x*r2 - 3*( p.x*d.x + p.y*d.y + p.z*d.z)*d.x )*scp,
    //     //( p.y*r2 - 3*( p.y*d.y + p.x*d.x + p.z*d.z)*d.y )*scp,
    //     //( p.z*r2 - 3*( p.z*d.z + p.x*d.x + p.y*d.y)*d.z )*scp,

    //     ( p.x*r2 - 3*pd*d.x )*scp,
    //     ( p.y*r2 - 3*pd*d.y )*scp,
    //     ( p.z*r2 - 3*pd*d.z )*scp,
    //     //( pd*3-r2)*d.x*scp,
    //     //( pd*3-r2)*d.y*scp,
    //     //( pd*3-r2)*d.z*scp
    // };
    //f.add(Fp);

    f.add_mul( d, (pd*3-r2)*scp );
    if( order<2 ) return E*ir;

    // ---------- Quadrupole
    //E += ir2*ir2*ir2();
    const Vec6d q  = *(Vec6d*)(cs+4);
    //const Vec6d dd{ d.x*d.x, d.y*d.y, d.z*d.z, d.y*d.z, d.x*d.z, d.x*d.y };
    const Vec3d dd{ d.x*d.x, d.y*d.y, d.z*d.z };
    const double Qd2 = -5*( q.xx*dd.x + q.yy*dd.y + q.zz*dd.z );
    const double xyz = 5*d.x*d.y*d.z;
    const double Qsum = q.xx + q.yy + q.zz;
    const double scq = scp *3* ir2;
    const Vec3d Fp{
        //   cc_c           ab_b  ab_a               
        //(   r2*d.x*(2*q.xx + Qsum)  +d.x*Qd2 + (r2-5*dd.xx)*( q.xy*d.y + q.xz*d.z )     + q_yz*xyz ) * scq,
        //(   r2*d.y*(2*q.yy + Qsum)  +d.y*Qd2 + (r2-5*dd.yy)*( q.xy*d.x + q.yz*d.z )     + q_xz*xyz ) * scq,
        //(   r2*d.z*(2*q.zz + Qsum)  +d.z*Qd2 + (r2-5*dd.zz)*( q.xz*d.x + q.yz*d.y )     + q_xy*xyz ) * scq
        (  d.x*( r2*(2*q.xx + Qsum) + Qd2) + (r2-5*dd.x)*( q.xy*d.y + q.xz*d.z )     + q.yz*xyz ) * scq,
        (  d.y*( r2*(2*q.yy + Qsum) + Qd2) + (r2-5*dd.y)*( q.xy*d.x + q.yz*d.z )     + q.xz*xyz ) * scq,
        (  d.z*( r2*(2*q.zz + Qsum) + Qd2) + (r2-5*dd.z)*( q.xz*d.x + q.yz*d.y )     + q.xy*xyz ) * scq
    };
    f.add(Fp);
    return E;
}

}

#endif

