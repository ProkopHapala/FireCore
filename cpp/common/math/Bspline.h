
#ifndef  Bspline_h
#define  Bspline_h

#include "quaternion.h"
//#include "CG.h"
#include "globals.h"


inline bool checkIndexRange( int imin, int imax, const char* name, int i, bool bExit=true, bool bPrint=true ){ 
    if( (i<imin) || (i>imax) ){ 
        if(bPrint)printf( "ERROR: %s(%i) out of range(%i,%i)\n", name, i, imin, imax ); 
        if(bExit )exit(0); 
        return true; 
    }
    return false;
}

inline bool checkIndexRange2( int imin, int imax, const char* name1, int i1, const char* name2, int i2, bool bExit=true, bool bPrint=true ){ 
    int i=i1+i2;
    if( (i<imin) || (i>imax) ){ 
        if(bPrint)printf( "ERROR: i(%i)=%s(%i)+%s(%i) out of range(%i,%i)\n", i, name1, i1, name2, i2, imin, imax ); 
        if(bExit )exit(0); 
        return true; 
    }
    return false;
}

inline bool outrange(int i, int imin, int imax ){
    if     (i<=imin){ return true; }
    else if(i>=imax){ return true; }
    return false;
}

inline int biwrap( int i, int n ){ return (i<=0    )? n-1 : -1; }
inline int diwrap( int i, int n ){ return (i>=(n-1))? 1-n :  1; }

/// Cubic B-spline 
/// See:  https://www.studocu.com/en-us/document/university-of-michigan/computer-graphics/lecture-notes-lecture-37/839972
/// https://web.eecs.umich.edu/~sugih/courses/eecs487/lectures/37-B-splines+NURBS.pdf
// https://dokumen.tips/documents/eecs-487-interactive-computer-graphics-sugihcourseseecs487-eecs-487-interactive.html?page=1


/*

Basis:

(1/6) *  u^3
(1/6) * -3u^3 + 3u^2 + 3u + 1   = 3*u*  ( 1 + u - u^2 )   +  1
(1/6) * +3u^3 - 6u^2 + 4        = 3*u^2*( u - 2       )   +  4
(1/6) * -u^3 + 3*u^2 - 3*u - 1  = (1-u)^3 



-3(1-t)^3 + 3(1-t)^2 + 3(1-t) + 1   = -3( -t^3 + 3*t^2 - 3*t - 1 ) + 3*( 1 - 2*t + t^2 ) + 3 - 3*t +1
= +3t^3     



dBasis:

(1/2)   u^2
(1/2)  -3u^2 + 2u + 1 
(1/2)   3u^2  - 4u
(1/2)  -(1-u)^2

ddBasis:

a1 = u
a2 =  -3*u + 1
a3 =   3*u - 2
a4 = 1-u

Points:

  Basis: [ 0.0, 1./6.,   2./3.,   1./6. ]
 dBasis: [ 0.0, 0.5,     0.0,    -0.5   ]
ddBasis: [ 0.0, 1.0,    -2.0,     1.0   ]
*/


/*

Quintic bspline

Here's the explicit quintic B-spline basis function with real number coefficients instead of fractions:

\[
B_i(x) = 

   0.00833 * x^5                                                              if  x < 1
  -0.04167 * x^5 + 0.25 * x^4 - 0.5  * x^3 +  0.50 * x^2 -  0.25 * x + 0.05   if  x < 2 
   0.08333 * x^5 - 0.75 * x^4 + 2.5  * x^3 -  4.25 * x^2 +  3.75 * x - 1.25   if  x < 3 
  -0.08333 * x^5 + 1.00 * x^4 - 4.25 * x^3 +  9.50 * x^2 - 11.25 * x + 5.25   if  x < 4 
   0.04167 * x^5 - 0.75 * x^4 + 4.25 * x^3 - 12.75 * x^2 + 18.75 * x - 10.5   if  x < 5
  -0.00833 * x^5 + 0.25 * x^4 - 3.0  * x^3 + 20.00 * x^2 - 62.50 * x + 75.0   if  x < 6 


        t5/120
-  (( 5*t5   -5*t4  -10*t3 -10*t2   -5*t  -1  )/120)
   (  5*t5  -10*t4  -10*t3 +10*t2  +25*t  +13 )/60
- ((  5*t5  -15*t4         +30*t2         -33 )/60)
  (   5*t5  -20*t4  +20*t3 +20*t2  -50*t  +26 )/120
 -((    t5   -5*t4  +10*t3 -10*t2   +5*t   -1 )/120)

       t5                                     /120
 (( -5*t5   +5*t4  +10*t3 +10*t2   +5*t  +1  )/120)
  (  5*t5  -10*t4  -10*t3 +10*t2  +25*t  +13 )/60
 (( -5*t5  +15*t4         -30*t2         +33 )/60)
  (  5*t5  -20*t4  +20*t3 +20*t2  -50*t  +26 )/120
 (( -1*t5   +5*t4  -10*t3 +10*t2   -5*t   +1 )/120)



       t5                                     /120
 (( -5*t5   +5*t4  +10*t3 +10*t2   +5*t  +1  )/120)
  ( 10*t5  -20*t4  -20*t3 +20*t2  +50*t  +26 )/120
 ((-10*t5  +30*t4  -60*t2                +66 )/120)
  (  5*t5  -20*t4  +20*t3 +20*t2  -50*t  +26 )/120
 (( -1*t5   +5*t4  -10*t3 +10*t2   -5*t   +1 )/120)

   0.008333333333333333*t5                                     
  -0.041666666666666666*t5  +0.041666666666666666*t4  +0.08333333333333333*t3 +0.08333333333333333*t2  +0.041666666666666666*t   +0.008333333333333333
   0.083333333333333333*t5  -0.166666666666666666*t4  -0.16666666666666666*t3 +0.16666666666666666*t2  +0.416666666666666666*t   +0.216666666666666666
  -0.083333333333333333*t5  +0.250000000000000000*t4                          -0.50000000000000000*t2                            +0.550000000000000000
   0.041666666666666666*t5  -0.166666666666666666*t4  +0.16666666666666666*t3 +0.16666666666666666*t2  -0.416666666666666666*t   +0.216666666666666666
  -0.008333333333333333*t5  +0.041666666666666666*t4  -0.08333333333333333*t3 +0.08333333333333333*t2  -0.041666666666666666*t   +0.008333333333333333


These are the same quintic B-spline basis functions, but now expressed with real numbers for the coefficients.

*/



//#include <math.h>
//#include <cstdlib>
//#include <stdio.h>

// ============ optimized

namespace Bspline{

// ==============================================
// ===================   Basis   ================
// ==============================================

constexpr static const double CPs[]{ 0.0,  1.0/6.0,   2.0/3.0,  1.0/6.0 };
constexpr static const double DPs[]{ 0.0,  0.5,       0.0,     -0.5     };
constexpr static const double APs[]{ 0.0,  1.0,      -2.0,      1.0     };


template<typename T>
Vec6T<T> basis5(T t){
    constexpr T inv6 = 1./6.;
    const T t2 = t*t;
    const T t3 = t2*t;
    const T t4 = t2*t2;
    const T t5 = t3*t2;
    return Vec6T<T>{
         0.008333333333333333*t5,                                     
        -0.041666666666666666*t5  +0.041666666666666666*t4  +0.08333333333333333*t3 +0.08333333333333333*t2  +0.041666666666666666*t   +0.008333333333333333, 
         0.083333333333333333*t5  -0.166666666666666666*t4  -0.16666666666666666*t3 +0.16666666666666666*t2  +0.416666666666666666*t   +0.216666666666666666,
        -0.083333333333333333*t5  +0.250000000000000000*t4                          -0.50000000000000000*t2                            +0.550000000000000000,
         0.041666666666666666*t5  -0.166666666666666666*t4  +0.16666666666666666*t3 +0.16666666666666666*t2  -0.416666666666666666*t   +0.216666666666666666,
        -0.008333333333333333*t5  +0.041666666666666666*t4  -0.08333333333333333*t3 +0.08333333333333333*t2  -0.041666666666666666*t   +0.008333333333333333
    };
}

template<typename T>
inline Quat4T<T> basis(T u){
    constexpr T inv6 = 1./6.;
    const T u2 = u*u;
    const T t  = 1-u;
    return Quat4T<T>{
        inv6*t*t*t,
        inv6*( 3*u2*( u - 2      )  +4 ),
        inv6*( 3*u *( 1 + u - u2 )  +1 ),
        inv6*u2*u,
    };
}

template<typename T>
inline Quat4T<T> dbasis(T u){
    const T u2 = u*u;
    const T t  = 1-u;
    return Quat4T<T>{
       -0.5*t*t,
        0.5*(  3*u2 - 4*u     ),
        0.5*( -3*u2 + 2*u + 1 ),
        0.5*u2,
    };
}

template<typename T>
inline Quat4T<T> ddbasis(T u){
    return Quat4T<T>{
        1-u,
        3*u-2,
        -3*u+1,
        u,
    };
}

// ======================================================  
// ===================   Interpolation   ================
// ======================================================

__attribute__((hot)) 
void sample1D( const double g0, const double dg, const int ng, const double* Gs, const int n, const double* ps, Vec2d* fes ){
    const double inv_dg = 1/dg; 
    for(int i=0; i<n; i++ ){
        const double x  = (ps[i] - g0)*inv_dg;  
        const int    ix = (int)x;
        if( ((ix<1)||(ix>=(ng-2))) )[[unlikely]]{ 
            continue;
            //printf( "ERROR: Bspline::sample1D() ixyz(%i,%i) out of range 0 .. (%i,%i) u(%g,%g)\n", ix,iy, n.x,n.y, u.x,u.y );   exit(0); 
        }
        const double tx = x-ix; 
        const Quat4d p  =  basis(tx);
        const Quat4d d  = dbasis(tx);
        const Quat4d gs = *(Quat4d*)(Gs+ix-1);
        fes[i]=Vec2d{
            p.dot( gs ),
            d.dot( gs )*inv_dg,
        };
    }
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d( const double tx, const double ty, const Quat4i i, const double* Es ){
    alignas(32) Quat4d e,fx;
    {
        const Quat4d bx =  basis( tx );
        const Quat4d dx = dbasis( tx );
        {
            alignas(32) const Quat4d p = *(Quat4d*)(Es+i.x); // read 4 doubles from global memory at a time ( 4*8 = 32 bytes = 256 bits ) ideal for SIMD AVX2
            e.x  = bx.dot(p);                                // not sure how dot() is SIMD optimized => maybe we should flip the order of x and y strides ?
            fx.x = dx.dot(p);
        }
        {
            alignas(32) const Quat4d p = *(Quat4d*)(Es+i.y); 
            e.y  = bx.dot(p);
            fx.y = dx.dot(p);
        }
        {
            alignas(32) const Quat4d p = *(Quat4d*)(Es+i.z); 
            e.z  = bx.dot(p);
            fx.z = dx.dot(p);
        }
        {
            alignas(32) const Quat4d p = *(Quat4d*)(Es+i.w); 
            e.w  = bx.dot(p);
            fx.w = dx.dot(p);
        }
    }
    alignas(32) const Quat4d by =  basis( ty );
    alignas(32) const Quat4d dy = dbasis( ty );
    return Vec3d{
        by.dot(fx), // Fx
        dy.dot(e ), // Fy
        by.dot(e )  // E
    };
}

//__attribute__((always_inline))
__attribute__((pure))
__attribute__((hot)) 
inline Vec2d fe1Dcomb3( const Vec3d* E, const Vec3d& C, const Quat4d& p, const Quat4d& d ){
    alignas(32) const Quat4d cs{ C.dot(E[0]), C.dot(E[1]), C.dot(E[2]), C.dot(E[3]) };
    return Vec2d{ p.dot( cs ), d.dot( cs ) };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d_comb3( int nz, const Vec3d* E, Quat4i di, const Vec3d& C, const Quat4d& pz, const Quat4d& dz, const Quat4d& by, const Quat4d& dy ){
    //const Quat4d* FEx = FE + ( i.x*n.y  + i.y )*n.z;
    //alignas(32) const Vec2d fe0 = fe1Dcomb3( E+di.x*nz, C, pz, dz );
    //alignas(32) const Vec2d fe1 = fe1Dcomb3( E+di.y*nz, C, pz, dz );
    //alignas(32) const Vec2d fe2 = fe1Dcomb3( E+di.z*nz, C, pz, dz );
    //alignas(32) const Vec2d fe3 = fe1Dcomb3( E+di.w*nz, C, pz, dz );
    // ---- more efficient version expacts di is already multiplied by nz
    alignas(32) const Vec2d fe0 = fe1Dcomb3( E+di.x, C, pz, dz );
    alignas(32) const Vec2d fe1 = fe1Dcomb3( E+di.y, C, pz, dz );
    alignas(32) const Vec2d fe2 = fe1Dcomb3( E+di.z, C, pz, dz );
    alignas(32) const Vec2d fe3 = fe1Dcomb3( E+di.w, C, pz, dz );
    return Vec3d{
        fe0.x*dy.x  +  fe1.x*dy.y  +  fe2.x*dy.z  +  fe3.x*dy.w,  // Fy
        fe0.y*by.x  +  fe1.y*by.y  +  fe2.y*by.z  +  fe3.y*by.w,  // Fz
        fe0.x*by.x  +  fe1.x*by.y  +  fe2.x*by.z  +  fe3.x*by.w   // E
    };
}


__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d( const Vec2d u, const Vec2i n, const double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ;
    if( ((ix<1)||(ix>=n.x-2)) || ((iy<1)||(iy>=n.y-2))  )[[unlikely]]{   
        return Vec3dZero;
        //printf( "ERROR: Bspline::fe2d() ixyz(%i,%i) out of range 0 .. (%i,%i) u(%g,%g)\n", ix,iy, n.x,n.y, u.x,u.y );   exit(0); 
    }
    //int i0 = ix + n.x*iy;
    int i0 = (ix-1) + n.x*(iy-1); 
    return fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
} 

__attribute__((hot)) 
void sample2D( const Vec2d g0, const Vec2d dg, const Vec2i ng, const double* Eg, const int n, const Vec2d* ps, Vec3d* fes ){
    printf( "Bspline::sample2D() ng[%i,%i] dg(%g,%g) g0(%g,%g)\n",   ng.x,ng.y, dg.x,dg.y,   g0.x,g0.y   );
    printf( "Bspline::sample2D() ps[0](%g,%g) ps[%i](%g,%g)\n",   ps[0].x,ps[0].y,  n-1,  ps[n-1].x,ps[n-1].y );
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        Vec3d fe = fe2d( (ps[i]-g0)*inv_dg, ng, Eg );        // sample3D(n=10000) time=2009.44[kTick] 200.944[tick/point]
        fe.x*=inv_dg.x;
        fe.y*=inv_dg.y;
        fes[i] = fe;
    }
}

__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d( const Vec3d u, const Vec3i n, const double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    if( 
        ((ix<1)||(ix>=n.x-2)) ||
        ((iy<1)||(iy>=n.y-2)) ||
        ((iz<1)||(iz>=n.z-2))        
    )[[unlikely]]{ 
         return Quat4dZero;
        //printf( "ERROR: Bspline::fe3d() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) u(%g,%g,%g)\n", ix,iy,iz, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); 
    }
    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;

    //printf( "ixyz(%i,%i,%i)  (%g,%g,%g)\n", ix,iy,iz,  u.x,u.y,u.z );

    int i0 = (ix-1) + n.x*( (iy-1) + n.y*(iz-1) );  
    
    const Vec3d Exy1 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
    int i1 = i0+nxy  ;                  const Vec3d Exy2 = fe2d(tx,ty, {i1,i1+n.x,i1+n.x*2,i1+3*n.x}, Es );
    int i2 = i0+nxy*2;                  const Vec3d Exy3 = fe2d(tx,ty, {i2,i2+n.x,i2+n.x*2,i2+3*n.x}, Es );
    int i3 = i0+nxy*3;                  const Vec3d Exy4 = fe2d(tx,ty, {i3,i3+n.x,i3+n.x*2,i3+3*n.x}, Es );

    const Quat4d bz =  basis( tz );
    const Quat4d dz = dbasis( tz );
    return Quat4d{
        bz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        bz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        dz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        bz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
    };
    //return Quat4d{0.0,0.0,0.0,Es[i0]};
} 

// static const Quat4i cubic_qis[4] = {
//     {0,1  ,2  ,3  },
//     {0,1  ,2  ,3-n},
//     {0,1  ,2-n,3-n},
//     {0,1-n,2-n,3-n}
// }; 

void make_inds_pbc( const int n, Quat4i* iqs ){
    iqs[0]={0,1  ,2  ,3  };
    iqs[1]={0,1  ,2  ,3-n};
    iqs[2]={0,1  ,2-n,3-n};
    iqs[3]={0,1-n,2-n,3-n};
}

inline Quat4i choose_inds_pbc( const int i, const int n, const Quat4i* iqs ){
    if(i>=(n-3))[[unlikely]]{ return iqs[i+4-n]; }
    return Quat4i{ 0, 1, 2, 3 };
}

inline int modulo( const int i, const int m ){
    int result = i % m;
    if (result < 0) {
        result += m;
    }
    return result;
}

__attribute__((pure))
__attribute__((hot)) 
Vec3d fe2d_pbc_comb3( const Vec2d u, const Vec2i n, const Vec3d* Es, const Vec3d PLQ, const Quat4i* yqs ){
	int          ix = (int)u.x  ,  iy = (int)u.y  ;
    if(u.y<0) iy--;
    const double tx = u.x - ix  ,  ty = u.y - iy  ;

    //printf( "Bspline::fe2d_pbc_comb3() ixy(%i,%i) u(%g,%g) n(%i,%i)\n", ix,iy, u.x,u.y,  n.x,n.y  );

    if(  ((ix<1)||(ix>=n.x-2))  )[[unlikely]]{  return Vec3dZero; }

    iy=modulo(iy-1,n.y);
    const Quat4i qy = choose_inds_pbc( iy, n.y, yqs );
    // if(ix==1){
    // printf( "fe2d_pbc_comb3() iy=%i iq=%i iy+qy{%i,%i,%i,%i}  qy{%i,%i,%i,%i}\n", iy,  (iy>=(n.y-3))?(iy+4-n.y):(-1),  iy+qy.x,iy+qy.y,iy+qy.z,iy+qy.w,  qy.x,qy.y,qy.z,qy.w  );
    // }

    const Quat4d bx =  basis( tx );
    const Quat4d dx = dbasis( tx );
    const Quat4d by =  basis( ty );
    const Quat4d dy = dbasis( ty );

    int i = (ix-1) + n.x*iy;  
    return fe2d_comb3( n.x, Es+i,  qy*n.x, PLQ, bx, dx, by, dy );


    //Vec3d Ei = Es[i]; 
    //printf( "fe2d_pbc_comb3() ixy[%i,%i] Etot=%g \n", ix,iy, PLQ.dot(Ei)  );
    //return Vec3d{  0.0, 0.0,  PLQ.dot(Ei) };
    //return Quat4d{0.0,0.0,0.0,Es[i0]};
} 


__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d_pbc_comb3( const Vec3d u, const Vec3i n, const Vec3d* Es, const Vec3d PLQ, const Quat4i* xqis, const Quat4i* yqis ){
    // We assume there are boundary added to simplify the index calculations
	int          ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    if(u.x<0) ix--;
    if(u.y<0) iy--;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;

    // printf( "Bspline::fe3d_pbc_comb3() ixyz(%3i,%3i,%3i)/n(%3i,%3i,%3i)  u(%g,%g,%g) \n", ix,iy,iz, n.x,n.y,n.z, u.x,u.y,u.z  ); 

    // ---- boundary conditions
    if(  ((iz<1)||(iz>=n.z-2))  )[[unlikely]]{  return Quat4dZero; }
    //if(  ((ix<1)||(ix>=n.x-2))  )[[unlikely]]{  return Quat4dZero; }
    //if(  ((iy<1)||(iy>=n.y-2))  )[[unlikely]]{  return Quat4dZero; }

    ix=modulo(ix-1,n.x);
    iy=modulo(iy-1,n.y);

    //printf( "Bspline::fe3d_pbc_comb3() ixyz(%3i,%3i,%3i)/n(%3i,%3i,%3i)  u(%g,%g,%g) \n", ix,iy,iz, n.x,n.y,n.z, u.x,u.y,u.z  ); 

    const int nyz = n.z*n.y;
    const Quat4i qx = choose_inds_pbc( ix, n.x, xqis )*nyz;
    const Quat4i qy = choose_inds_pbc( iy, n.y, yqis )*n.z;

    const Quat4d bz =  basis( tz );
    const Quat4d dz = dbasis( tz );
    const Quat4d by =  basis( ty );
    const Quat4d dy = dbasis( ty );
    //int i0 = (iz-2) + n.z*( iy + n.y*ix);  
    int i0 = (iz-1) + n.z*( iy + n.y*ix); 
    const Vec3d E1 = fe2d_comb3( n.z, Es+(i0+qx.x ), qy, PLQ, bz, dz, by, dy );
    const Vec3d E2 = fe2d_comb3( n.z, Es+(i0+qx.y ), qy, PLQ, bz, dz, by, dy );
    const Vec3d E3 = fe2d_comb3( n.z, Es+(i0+qx.z ), qy, PLQ, bz, dz, by, dy );;
    const Vec3d E4 = fe2d_comb3( n.z, Es+(i0+qx.w ), qy, PLQ, bz, dz, by, dy );
    const Quat4d bx =  basis( tx );
    const Quat4d dx = dbasis( tx );
    return Quat4d{
        dx.dot( {E1.z, E2.z, E3.z, E4.z} ), // Fx
        bx.dot( {E1.x, E2.x, E3.x, E4.x} ), // Fy
        bx.dot( {E1.y, E2.y, E3.y, E4.y} ), // Fz
        bx.dot( {E1.z, E2.z, E3.z, E4.z} ), // E
    };
} 

__attribute__((hot)) 
void sample2D_comb3( const Vec2d g0, const Vec2d dg, const Vec2i ng, const Vec3d* Eg, const int n, const Vec2d* ps, Vec3d* fes, Vec3d C ){
    printf( "Bspline::sample2D_comb3() ng[%i,%i] dg(%g,%g) g0(%g,%g) C(%g,%g,%g) n=%i \n",   ng.x,ng.y,   dg.x,dg.y,   g0.x,g0.y,   C.x,C.y,C.z, n  );
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    Quat4i yqs[4];
    make_inds_pbc( ng.y, yqs );
    for(int i=0; i<n; i++ ){
        //printf( "sample2D_comb3()[%i] \n", i );
        Vec3d fe = fe2d_pbc_comb3( (ps[i]-g0)*inv_dg, ng, Eg, C, yqs ); 
        fe.x *= inv_dg.x;
        fe.y *= inv_dg.y;
        fes[i] = fe;
    }
}

__attribute__((hot)) 
void sample3D_comb3( const Vec3d g0, const Vec3d dg, const Vec3i ng, const Vec3d* Eg, const int n, const Vec3d* ps, Quat4d* fes, Vec3d C ){
    printf( "Bspline::sample3D_comb3() ng[%i,%i,%i] dg(%g,%g,%g) g0(%g,%g,%g) C(%g,%g,%g) n=%i \n",   ng.x,ng.y,ng.z,   dg.x,dg.y,dg.z,   g0.x,g0.y,g0.z,   C.x,C.y,C.z, n  );
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    Quat4i xqs[4]; make_inds_pbc( ng.x, xqs );
    Quat4i yqs[4]; make_inds_pbc( ng.y, yqs );
    for(int i=0; i<n; i++ ){
        //printf( "sample3D_comb3()[%i] \n", i );
        Quat4d fe = fe3d_pbc_comb3( (ps[i]-g0)*inv_dg, ng, Eg, C, xqs, yqs ); 
        fe.f.mul( inv_dg );
        fes[i] = fe;
    }
}


__attribute__((hot)) 
void sample3D( const Vec3d g0, const Vec3d dg, const Vec3i ng, const double* Eg, const int n, const Vec3d* ps, Quat4d* fes ){
    printf( "Bspline::sample3D() ng[%i,%i,%i] dg(%g,%g,%g) g0(%g,%g,%g)\n",   ng.x,ng.y,ng.z, dg.x,dg.y,dg.z,   g0.x,g0.y,g0.z   );
    printf( "Bspline::sample3D() ps[0](%g,%g,%g) ps[%i](%g,%g,%g)\n",   ps[0].x,ps[0].y,ps[0].z,  n-1,  ps[n-1].x,ps[n-1].y,ps[n-1].z );
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        //printf( "Bspline::sample3D()[%i] p(%g,%g,%g) g0(%g,%g,%g) dg(%g,%g,%g)\n", i,  ps[i].x,ps[i].y,ps[i].z,    g0.x,g0.y,g0.z,    dg.x,dg.y,dg.z  );
        Quat4d fe = fe3d( (ps[i]-g0)*inv_dg, ng, Eg );        // sample3D(n=10000) time=2009.44[kTick] 200.944[tick/point]
        //Quat4d fe = fe3d_v2( (ps[i]-g0)*inv_dg, ng, Eg );   // sample3D(n=10000) time=2175.84[kTick] 217.584[tick/point]
        fe.f.mul(inv_dg);
        fes[i] = fe;
    }
}






// ================================================  
// ===================   Fitting   ================
// ================================================  

__attribute__((hot)) 
Vec3d move( double dt, int n, double* gs, double* fs, double* vs ){
    //constexpr const int m = 2;
    // --- eval velocity-to-force projection 
    double vf = 0.0;
    double ff = 0.0;
    double vv = 0.0;
    for(int i=0; i<n; i++){
        ff += fs[i]*fs[i];
        vf += vs[i]*fs[i];
        vv += vs[i]*vs[i];
    }
    //printf( "|F[%i]|=%g \n",itr,  sqrt(ff) );
    //printf( "p=%i v=%g f=%g \n",  Gs[1], vs[1], fs[1] );
    //if(ff<F2max){ break; }
    if(vf<0){ for(int i=0; i<n; i++){ vs[i]=0; }; }
    // --- update the points
    for(int i=0; i<n; i++){
        vs[i] += fs[i]*dt;
        gs[i] += vs[i]*dt;
        //Gs[i] += fs[i]*dt;  // GD
    }
    return Vec3d{vf,ff,vv};
}


__attribute__((hot)) 
Vec3d move_GD( double dt, int n, double* gs, double* fs ){
    double ff = 0.0;
    for(int i=0; i<n; i++){
        double f = fs[i];
        gs[i]   += f*dt;
        ff      += f*f;
    }
    return Vec3d{0.0,ff,0.0};
}


__attribute__((hot)) 
int fit1D_old( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1 ){
    if(verbosity>1)printf("Bspline::fit1D_old() \n");
    const double F2max = Ftol*Ftol;
    double* ps = new double[n];
    double* fs = new double[n];
    double* vs = new double[n];
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    int itr=0;
    for(int i=0; i<n; i++){ vs[i]=0.0; } // clear velocity
    //if( ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) && ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) ){ printf("ERROR in fit1D()=> exit\n"); exit(0); };
    bool bErr = false;
    double err2sum=0;
    for(itr=0; itr<nmaxiter; itr++){
        err2sum=0.0;
        // --- evaluate current spline
        for(int i=0; i<n; i++){
        //for(int i=1; i<n-1; i++){
            double p = Gs[i]*B0;
            if(i>0  )[[likely]]{ p += B1*Gs[i-1]; }
            if(i<n-1)[[likely]]{ p += B1*Gs[i+1]; }

            //p += B1*Gs[i-1];
            //p += B1*Gs[i+1];

            //bErr|=checkNumRange( i, T val, T min, T max, const char* pre, bool bPrint=true, bool bExit=false ){
            //printf( "p[%i]=%g Gs=%g\n", i, p, Gs[i] );
            //bErr|=  checkNumRange( i, p, -1.e+6, +1.e+6, "p=B*Gs", true, true );
            ps[i] = p;   //Ws[i] = p;
            fs[i] = 0;
        }
        // --- evaluate variatiaonal derivatives
        //for(int i=0; i<n; i++){
        for(int i=1; i<n-1; i++){
            double dp  = Es[i] - ps[i];
            if( Ws ){ dp *= Ws[i]; } // weighting 
            err2sum += dp*dp;
            fs[i] += dp*B0;
            //bErr|=  checkNumRange( i, dp, -1.e+6, +1.e+6, "dp=Es-ps", true, true );
            //if(i>0  )[[likely]]{ fs[i-1] += B1*dp; }
            //if(i<n-1)[[likely]]{ fs[i+1] += B1*dp; }
            fs[i-1] += B1*dp;
            fs[i+1] += B1*dp;
        }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "Bspline::fit1D_old() iter=%i err=%g \n", itr, sqrt(err2sum) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
double getVariations1D( const int n, const double* Gs, const double* Es, const double* Ws, double* fs, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    // --- evaluate current spline
    for(int i=0; i<n; i++){ fs[i]=0; ps[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int i=1; i<n-1; i++){
        double val = Gs[i-1]*B1 + Gs[i]*B0 + Gs[i+1]*B1;
        double err = Es[i] - val;
        if(Ws){ err*=Ws[i]; }
        err2sum += err*err;
        ps[i] = err;
        //Ws[j] = err;
    }
    // --- distribute variational derivatives of approximation error
    for(int i=0; i<n; i++){
        if( (i>0)&&(i<(n-1)) ) [[likely]] {
            double val = ps[i-1]*B1 + ps[i]*B0 + ps[i+1]*B1;
            fs[i] = val;
        }else{
            fs[i] = 0;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
double getVariations1D_half( const int n, const double* Gs, const double* Es, const double* Ws, double* fs, double* ers ){
    constexpr double B00 =2.0/3.0;   // 0.66666666666666666666666666666667
    constexpr double B05 =23.0/48.0; // 0.47916666666666666666666666666667
    constexpr double B10 =1.0/6.0;   // 0.16666666666666666666666666666667
    constexpr double B15 =1.0/48.0;  // 0.020833333333333333333333333333333
    
    int n2 = 2*n;
    // --- evaluate current spline
    for(int i=0; i<n;  i++){ fs[i]=0;   }
    for(int i=0; i<n2; i++){ ers[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int i=2; i<n-2; i++){
        double e0 = Gs[i-1]*B10 + Gs[i]*B00 + Gs[i+1]*B10;
        double e1 = Gs[i-1]*B15 + Gs[i]*B05 + Gs[i+1]*B05 + Gs[i+2]*B15;
        //double e1 = Gs[i-2]*B15 + Gs[i-1]*B05 + Gs[i]*B05 + Gs[i+1]*B15;

        int i2=i*2;
        double d0 = Es[i2  ] - e0;   ers[i2  ] = d0;
        double d1 = Es[i2+1] - e1;   ers[i2+1] = d1;
        
        err2sum += e0*e0 + e1*e1;
        
        //Ws[j] = err;
    }
    // --- distribute variational derivatives of approximation error
    for(int i=0; i<n; i++){
        if( (i>0)&&(i<(n-1)) ) [[likely]] {
            int i2=i*2;
            //double f = ers[i2-3]*B15 + ers[i2-2]*B10 + ers[i2-1]*B05 + ers[i2]*B00 + ers[i2+1]*B05 + ers[i2+2]*B10 + ers[i2+3]*B15;

            double f = ers[i2-2]*B10 + ers[i2]*B00 + ers[i2+2]*B10;

            f += ers[i2-1]*B15 + ers[i2+1]*B15;
            //if(i2>=3){ f+=ers[i2-3]*B15; }
            //if(i2<n ){ f+=ers[i2+3]*B15; }

            fs[i]    = f;
        }else{
            fs[i] = 0;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
int fit1D( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bHalf=false ){
    if(verbosity>1)printf("Bspline::fit1D() !!!!!!\n");
    const double F2max = Ftol*Ftol;
    int n2=n; if(bHalf){ n2=2*n; }
    double* ps = new double[n2];
    double* fs = new double[n];
    double* vs = new double[n];
    int itr=0;
    for(int i=0; i<n; i++){ vs[i]=0.0; } // clear velocity
    //if( ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) && ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) ){ printf("ERROR in fit1D()=> exit\n"); exit(0); };
    bool bErr = false;
    double err2sum=0;
    for(itr=0; itr<nmaxiter; itr++){        
        if(bHalf){ err2sum = getVariations1D_half( n, Gs, Es, Ws, fs, ps ); }
        else     { err2sum = getVariations1D(      n, Gs, Es, Ws, fs, ps ); }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "Bspline::fit1D() iter=%i err=%g \n", itr, sqrt(err2sum) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
int fit1D_EF( const double dg, const int n, double* Gs,  Vec2d* fes, Vec2d* Ws, double Ftol, int nmaxiter=1000, double dt=0.1 ){
    if(verbosity>1)printf("Bspline::fit1D_EF() \n");
    const double inv_dg = 1/dg;
    const double F2max = Ftol*Ftol;
    Vec2d*  ps = new Vec2d[n];
    double* fs = new double[n];
    double* vs = new double[n];
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double D1=0.5;
    int itr=0;
    for(int i=0; i<n; i++){ 
        vs[i] = 0; 
        ps[i] = Vec2dZero;   
        //Ws[i] = Vec2dZero;
    };
    double err=0;
    for(itr=0; itr<nmaxiter; itr++){
        err=0.0;
        // --- evaluate current spline
        for(int i=0; i<n; i++){ fs[i]=0; };
        for(int i=1; i<n-1; i++){
            Vec2d p = Vec2d{ Gs[i]*B0, 0.0};
            // if(i>0  )[[likely]]{ p.add_mul( Vec2d{B1,-D1*inv_dg}, Gs[i-1] ); }
            // if(i<n-1)[[likely]]{ p.add_mul( Vec2d{B1,+D1*inv_dg}, Gs[i+1] ); }
            p.add_mul( Vec2d{B1,-D1*inv_dg}, Gs[i-1] );
            p.add_mul( Vec2d{B1,+D1*inv_dg}, Gs[i+1] ); 
            ps[i] = p;    
            //Ws[i] = p;
            //Ws[i].x = p.y;
        }
        // --- evaluate variatiaonal derivatives
        for(int i=2; i<n-2; i++){
            Vec2d fei = fes[i];
            //fei.y*=-1;
            //const Vec2d w   = Ws[i];
            Vec2d dp  = fei - ps[i];
            if(Ws){ dp.mul( Ws[i] ); }
            err+= dp.norm2();
            fs[i] += dp.x*B0;
            //if(i>0  )[[likely]]{ fs[i-1] += B1*dp.x*0 + -D1*dp.y; }
            //if(i<n-1)[[likely]]{ fs[i+1] += B1*dp.x*0 + +D1*dp.y; }
            //fs[i-1] += B1*dp.x + -D1*dp.y*0;
            //fs[i+1] += B1*dp.x + +D1*dp.y*0;
            fs[i  ] += B0*dp.x;
            fs[i-1] += B1*dp.x + -D1*dp.y;
            fs[i+1] += B1*dp.x + +D1*dp.y;
        }
        //for(int i=0; i<n; i++){  Ws[i].y = fs[i];  } // Debug
        // --- move
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "Bspline::fit1D_EF() iter=%i err=%g \n", itr, err );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}


__attribute__((pure)) 
__attribute__((hot)) 
inline double assemleBound2D( const double B00, const double B01, const double B11, const double* Gs, const int i, const int nx, const bool ylo, const bool yhi){
    const bool xlo = i > 0;
    const bool xhi = i < nx-1; 
    const int  i0  = i-nx;
    const int  i1  = i+nx;
    double val=0;
    if( xhi && xlo && ylo && yhi) [[likely]] { 
        val = Gs[i0-1]*B11 + Gs[i0]*B01 + Gs[i0+1]*B11
            + Gs[i -1]*B01 + Gs[i ]*B00 + Gs[i +1]*B01
            + Gs[i1-1]*B11 + Gs[i1]*B01 + Gs[i1+1]*B11;
    }else{
                       val =Gs[i ]*B00;
        if( ylo ){     val+=Gs[i0]*B01; };
        if( yhi ){     val+=Gs[i1]*B01; };
        if( xlo ){ 
                       val+=Gs[i -1]*B01; 
            if( ylo ){ val+=Gs[i0-1]*B11; };
            if( yhi ){ val+=Gs[i1-1]*B11; };
        };
        if( xhi ){
                       val+=Gs[i +1]*B01;  
            if( ylo ){ val+=Gs[i0+1]*B11; };
            if( yhi ){ val+=Gs[i1+1]*B11; };
        };
    }
    return val;
}

__attribute__((pure)) 
__attribute__((hot)) 
inline double assemleBound2D_pbc( const double B00, const double B01, const double B11, const double* Gs, int i, int ibx, int idx, int iby, int idy ){
    return   Gs[i+ibx+iby]*B11 + Gs[i+iby]*B01 + Gs[i+idx+iby]*B11
           + Gs[i+ibx    ]*B01 + Gs[i    ]*B00 + Gs[i+idx    ]*B01
           + Gs[i+ibx+idy]*B11 + Gs[i+idy]*B01 + Gs[i+idx+idy]*B11;
}


__attribute__((hot)) 
double getVariations2D( const Vec2i ns, double* Gs,  const double* Es, double* Ws, double* fs, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B00=B0*B0;
    constexpr double B01=B0*B1;
    constexpr double B11=B1*B1;
    //const int nxy  = ns.x*ns.y;
    //for(int i=0; i<nxy; i++){ fs[i]=0; ps[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int iy=0; iy<ns.y; iy++){
        int iiy = iy*ns.x;
        const bool ylo  = iy > 0;
        const bool yhi  = iy < ns.y-1;
        for(int ix=0; ix<ns.x; ix++){
            double val = assemleBound2D( B00,B01,B11, Gs+iiy, ix, ns.x, ylo, yhi );
            const int i = ix + iiy;
            double err = Es[i] - val;
            //if((ix==6)&&(iy==10)){  printf("getVariations2D_mod()[%i,%i] E=%g val=%g err=%g \n", ix,iy, Es[i], val, err ); }
            //if(Ws){ err*=Ws[i]; }
            err2sum += err*err;
            ps[i] = err;
            Ws[i] = err;
        }
    }
    // --- distribute variational derivatives of approximation error
    for(int iy=0; iy<ns.y; iy++){
        int iiy = iy*ns.x;
        const bool ylo  = iy > 0;
        const bool yhi  = iy < ns.y-1;
        for(int ix=0; ix<ns.x; ix++){
            //printf("c2 ix,iy: %3i %3i \n", ix,iy );
            double f= assemleBound2D( B00,B01,B11, ps+iiy, ix, ns.x, ylo, yhi );
            const int i = ix + iiy;
            //Ws[i] = f;
            fs[i] = f;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
double getVariations2D_pbc( const Vec2i ns, double* Gs,  const double* Es, double* Ws, double* fs, double* ps ){
    //printf("getVariations2D_pbc()\n");
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B00=B0*B0;
    constexpr double B01=B0*B1;
    constexpr double B11=B1*B1;
    //const int nxy  = ns.x*ns.y;
    //for(int i=0; i<nxy; i++){ fs[i]=0; ps[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    int nxy=ns.x*ns.y;
    for(int iy=0; iy<ns.y; iy++){
        //printf("getVariations2D_pbc() iy=%i \n", iy );
        const int iiy = iy*ns.x;
        const int iby = biwrap(iy,ns.y)*ns.x;
        const int idy = diwrap(iy,ns.y)*ns.x;
        //checkIndexRange2( 0, (ns.y-1)*ns.x, "iiy", iiy, "iby", iby ); 
        //checkIndexRange2( 0, (ns.y-1)*ns.x, "iiy", iiy, "idy", idy ); 
        for(int ix=0; ix<ns.x; ix++){
            const int i   = ix + iiy;
            const int ibx = biwrap(ix,ns.x);
            const int idx = diwrap(ix,ns.x);
            //printf("getVariations2D_pbc() [%2i,%2i] %4i  xbd(%2i,%2i)  ybd(%4i,%4i)\n", iy, ix,  i, ibx+ix,idx+ix,   iby+iy*ns.x,idy+iy*ns.x,     i+ibx+iby,  i+idx+idy, nxy );
            //printf("getVariations2D_pbc() [%2i,%2i] %4i  xbd(%2i,%2i)  ybd(%4i,%4i)   {%4i,%4i,%4i, %4i,%4i,%4i, %4i,%4i,%4i}   nxy=%i \n", iy, ix,  i, ibx,idx, iby,idy,     i+iby+ibx,i+iby,i+iby+idx,     i+ibx,i,i+idx,  i+idy+ibx,i+idy,i+idy+idx, nxy );
            //checkIndexRange2( 0, ns.x-1, "ix", ix, "ibx", ibx );
            //checkIndexRange2( 0, ns.x-1, "ix", ix, "idx", idx );
            //checkIndexRange( 0, nxy-1, "i", i );
            double val = assemleBound2D_pbc( B00,B01,B11,Gs,i,ibx,idx,iby,idy );
            double err = Es[i] - val;
            //if((ix==6)&&(iy==10)){  printf("getVariations2D_mod()[%i,%i] E=%g val=%g err=%g \n", ix,iy, Es[i], val, err ); }
            //if(Ws){ err*=Ws[i]; }
            err2sum += err*err;
            ps[i] = err;
            Ws[i] = err;
        }
    }
    // --- distribute variational derivatives of approximation error
    for(int iy=0; iy<ns.y; iy++){
        const int iiy = iy*ns.x;
        const int iby = biwrap(iy,ns.y)*ns.x;
        const int idy = diwrap(iy,ns.y)*ns.x;
        for(int ix=0; ix<ns.x; ix++){
            const int i   = ix + iiy;
            const int ibx = biwrap(ix,ns.x);
            const int idx = diwrap(ix,ns.x);;
            double f= assemleBound2D_pbc( B00,B01,B11,ps,i,ibx,idx,iby,idy );
            //Ws[i] = f;
            fs[i] = f;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
double error_iz( int iz, const Vec3i ns, const double* Gs, const double* Es, const double* Ws, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B000=B0*B0*B0;
    constexpr double B001=B0*B0*B1;
    constexpr double B011=B0*B1*B1;
    constexpr double B111=B1*B1*B1;
    const int nxy  = ns.x*ns.y;
    double err2sum=0.0;
    const int  iiz  = iz*nxy;
    const bool zlo  = iz > 0;
    const bool zhi  = iz < ns.z-1;
    for(int iy=0; iy<ns.y; iy++){
        const int  iiy = iy*ns.x;
        const int  iyz = iiz+iiy;
        const bool ylo = iy > 0;
        const bool yhi = iy < ns.y-1;
        for(int ix=0; ix<ns.x; ix++){
            
            double  val  = assemleBound2D( B000,B001,B011, Gs+iyz    , ix, ns.x, ylo, yhi ); 
            if(zlo) val += assemleBound2D( B001,B011,B111, Gs+iyz-nxy, ix, ns.x, ylo, yhi ); 
            if(zhi) val += assemleBound2D( B001,B011,B111, Gs+iyz+nxy, ix, ns.x, ylo, yhi ); 
            const int i = ix + iyz;
            double err = Es[i] - val;
            //if((ix==6)&&(iy==10)&&(iz==10)){  printf("getVariations3D_mod()[%i,%i,%i] E=%g val=%g err=%g \n", ix,iy,iz, Es[i], val, err ); }
            //Ws[i] = err;
            if(Ws){ err*=Ws[i]; }
            err2sum += err*err;
            ps[i]    = err;
            //Ws[i] = err;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
void force_iz( int iz, const Vec3i ns, const double* ps, double* fs ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B000=B0*B0*B0;
    constexpr double B001=B0*B0*B1;
    constexpr double B011=B0*B1*B1;
    constexpr double B111=B1*B1*B1;
    const int nxy  = ns.x*ns.y;
    //double err2sum=0.0;
    const int  iiz  = iz*nxy;
    const bool zlo  = iz > 0;
    const bool zhi  = iz < ns.z-1;
    for(int iy=0; iy<ns.y; iy++){
        const int iiy  = iy*ns.x;
        const int iyz  = iiz+iiy;
        const bool ylo = iy > 0;
        const bool yhi = iy < ns.y-1;
        for(int ix=0; ix<ns.x; ix++){
            double  val  = assemleBound2D( B000,B001,B011, ps+iyz    , ix, ns.x, ylo, yhi ); 
            if(zlo) val += assemleBound2D( B001,B011,B111, ps+iyz-nxy, ix, ns.x, ylo, yhi ); 
            if(zhi) val += assemleBound2D( B001,B011,B111, ps+iyz+nxy, ix, ns.x, ylo, yhi );     
            const int i = ix + iyz;
            //val*=-1;
            fs[i] = val;
        }
    }
    //return err2sum;
}

__attribute__((hot)) 
double error_iz_pbc( int iz, const Vec3i ns, const double* Gs, const double* Es, double* Ws, double* ps ){
    //printf("error_iz_pbc(iz=%i) \n");
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B000=B0*B0*B0;
    constexpr double B001=B0*B0*B1;
    constexpr double B011=B0*B1*B1;
    constexpr double B111=B1*B1*B1;
    const int nxy  = ns.x*ns.y;
    double err2sum=0.0;
    const int  iiz  = iz*nxy;
    const bool zlo  = iz > 0;
    const bool zhi  = iz < ns.z-1;
    for(int iy=0; iy<ns.y; iy++){
        const int iiy = iy*ns.x + iiz;
        const int iby = biwrap(iy,ns.y)*ns.x;
        const int idy = diwrap(iy,ns.y)*ns.x;
        for(int ix=0; ix<ns.x; ix++){
            
            const int i   = ix + iiy;
            const int ibx = biwrap(ix,ns.x);
            const int idx = diwrap(ix,ns.x);

            // if( outrange(i+iby+ibx,0,nxy)
            // ||  outrange(i+iby    ,0,nxy)
            // ||  outrange(i+iby+idx,0,nxy)
            // ||  outrange(i    +ibx,0,nxy)
            // ||  outrange(i        ,0,nxy)
            // ||  outrange(i    +idx,0,nxy)
            // ||  outrange(i+idy+ibx,0,nxy)
            // ||  outrange(i+idy    ,0,nxy)
            // ||  outrange(i+idy+idx,0,nxy) 
            // )[[unlikely]]{ printf("getVariations2D_pbc() [%2i,%2i] %4i  xbd(%2i,%2i)  ybd(%4i,%4i)   {%4i,%4i,%4i, %4i,%4i,%4i, %4i,%4i,%4i}   nxy=%i \n", iy, ix,  i, ibx,idx, iby,idy,     i+iby+ibx,i+iby,i+iby+idx,     i+ibx,i,i+idx,  i+idy+ibx,i+idy,i+idy+idx, nxy );  }
            //printf("getVariations2D_pbc() [%2i,%2i] %4i  xbd(%2i,%2i)  ybd(%4i,%4i)   {%4i,%4i,%4i, %4i,%4i,%4i, %4i,%4i,%4i}   nxy=%i \n", iy, ix,  i, ibx,idx, iby,idy,     i+iby+ibx,i+iby,i+iby+idx,     i+ibx,i,i+idx,  i+idy+ibx,i+idy,i+idy+idx, nxy );

            double  val  = assemleBound2D_pbc( B000,B001,B011, Gs    , i, ibx,idx,  iby,idy ); 
            if(zlo) val += assemleBound2D_pbc( B001,B011,B111, Gs-nxy, i, ibx,idx,  iby,idy ); 
            if(zhi) val += assemleBound2D_pbc( B001,B011,B111, Gs+nxy, i, ibx,idx,  iby,idy );

            double err = Es[i] - val;
            //if((ix==6)&&(iy==10)&&(iz==10)){  printf("getVariations3D_mod()[%i,%i,%i] E=%g val=%g err=%g \n", ix,iy,iz, Es[i], val, err ); }
            //Ws[i] = err;
            //if(Ws){ err*=Ws[i]; }
            err2sum += err*err;
            ps[i]    = err;
            //Ws[i] = err;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
void force_iz_pbc( int iz, const Vec3i ns, const double* ps, double* fs ){
    //printf("force_iz_pbc(iz=%i) \n");
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B000=B0*B0*B0;
    constexpr double B001=B0*B0*B1;
    constexpr double B011=B0*B1*B1;
    constexpr double B111=B1*B1*B1;
    const int nxy  = ns.x*ns.y;
    //double err2sum=0.0;
    const int  iiz  = iz*nxy;
    const bool zlo  = iz > 0;
    const bool zhi  = iz < ns.z-1;
    for(int iy=0; iy<ns.y; iy++){ 
        const int iiy = iy*ns.x + iiz;
        const int iby = biwrap(iy,ns.y)*ns.x;
        const int idy = diwrap(iy,ns.y)*ns.x;
        for(int ix=0; ix<ns.x; ix++){

            const int i   = ix + iiy;
            const int ibx = biwrap(ix,ns.x);
            const int idx = diwrap(ix,ns.x);

            double  val  = assemleBound2D_pbc( B000,B001,B011, ps    , i, ibx,idx,  iby,idy ); 
            if(zlo) val += assemleBound2D_pbc( B001,B011,B111, ps-nxy, i, ibx,idx,  iby,idy ); 
            if(zhi) val += assemleBound2D_pbc( B001,B011,B111, ps+nxy, i, ibx,idx,  iby,idy );    
            
            //val*=-1;
            fs[i] = val;
        }
    }
    //return err2sum;
}

__attribute__((hot)) 
double getVariations3D( const Vec3i ns, const double* Gs, const double* Es, double* Ws, double* fs, double* ps, bool bPBC ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B000=B0*B0*B0;
    constexpr double B001=B0*B0*B1;
    constexpr double B011=B0*B1*B1;
    constexpr double B111=B1*B1*B1;
    //double sum_B = B000 + B001*6 + B011*12 + B111*8;
    //printf( "getVariations3D_mod() sum_B %g \n", sum_B );
    const int nxy  = ns.x*ns.y;
    const int nxyz = nxy*ns.z;
    // --- evaluate current spline
    //for(int i=0; i<nxyz; i++){ fs[i]=0; ps[i]=0;  }
    double err2sum = 0;
    // --- evaluate current spline (in order to evelauet approximation error)
    for(int iz=0; iz<ns.z; iz++){ 
        //printf("getVariations3D_mod2().error[iz=%i] \n", iz);
        if(bPBC){ err2sum += error_iz_pbc( iz, ns, Gs, Es, Ws, ps ); }
        else    { err2sum += error_iz    ( iz, ns, Gs, Es, Ws, ps ); }            
    }
    for(int iz=0; iz<ns.z; iz++){ 
        //printf("getVariations3D_mod2().force[iz=%i] \n", iz);
        if(bPBC){ force_iz_pbc( iz, ns, ps, fs ); }
        else    { force_iz    ( iz, ns, ps, fs );     }
    }
    //printf("getVariations3D_mod() err %g ns(%i,%i,%i) \n", err2sum, ns.x, ns.y, ns.z );
    return err2sum;
}

__attribute__((hot)) 
int fit2D( const Vec2i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false ){
    //if(verbosity>1)
    printf( "Bspline::fit2D() ns(%i,%i) Ftol=%g dt=%g nmaxiter=%i bPBC=%i \n", ns.x,ns.y, Ftol, dt, nmaxiter, bPBC );
    const double F2max = Ftol*Ftol;
    const int nxy  = ns.x*ns.y;
    double* ps = new double[nxy];
    double* fs = new double[nxy];
    double* vs = new double[nxy];
    int itr=0;
    //while(false){
    double err;    
    Vec3d  cfv;
    for(int i=0; i<nxy; i++){ vs[i]=0; };
    for(itr=0; itr<nmaxiter; itr++){
        //for(int i=0; i<nxy; i++){ fs[i]=0; };  // Not necessary - force is cleared inside getVariations2D
        if(bPBC){ err = getVariations2D_pbc( ns, Gs, Es, Ws, fs, ps ); }
        else    { err = getVariations2D    ( ns, Gs, Es, Ws, fs, ps ); }
        cfv = move(dt,nxy,Gs,fs,vs);
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g Error=%g \n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z), sqrt(err) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "|F[%i]|=%g Error=%g \n",itr,sqrt(cfv.y), sqrt(err) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
int fit3D( const Vec3i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false, bool bInitGE=false ){
    long t0 = getCPUticks();
    //if(verbosity>1)
    printf( "Bspline::fit3D() ns(%i,%i,%i)  Ftol=%g dt=%g nmaxiter=%i bPBC=%i \n", ns.x,ns.y,ns.z, Ftol, dt, nmaxiter, bPBC  );
    const double F2max = Ftol*Ftol;
    const int nxy  = ns.x*ns.y;
    const int nxyz = ns.x*ns.y*ns.z;
    double* ps = new double[nxyz];
    double* fs = new double[nxyz];
    double* vs = new double[nxyz];
    int itr=0;
    //while(false){
    double err=0; 
    Vec3d  cfv;
    if(bInitGE){ for(int i=0; i<nxyz; i++){ Gs[i]=Es[i]; }; };
    for(int i=0; i<nxyz; i++){ vs[i]=0; };
    //dt = 0.3;
    for(itr=0; itr<nmaxiter; itr++){
        //for(int i=0; i<nxyz; i++){ fs[i]=0; };
        err = getVariations3D( ns, Gs, Es, Ws, fs, ps, bPBC );
        //err = getVariations3D_omp( ns, Gs, Es, Ws, fs, ps );
        cfv = move(dt,nxyz,Gs,fs,vs);
        //cfv = move_GD( dt, nxyz, Gs, fs );
        //if(verbosity>2)
        printf( "|F[%i]|=%g Error=%g \n",itr,sqrt(cfv.y), sqrt(err) );
        ///printf( "|F[%i]|=%g cos(f,v)=%g Error=%g \n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z), sqrt(err) );
        if(cfv.y<F2max){ break; };
    }
    //if(verbosity>1)
    printf( "Bspline::fit3D() DONE |F[%i]|=%g Error=%g \n",itr,sqrt(cfv.y), sqrt(err) );    
    double t = getCPUticks()-t0; printf( "Bspline::fit3D() niter=%i nxyz=%i time %g[GTicks] %g[tick/(nxyz*iter)] \n", itr, t*1e-9, t/(nxyz*itr) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}




__attribute__((hot)) 
int fit3D_omp( const Vec3i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false, bool bInitGE=false ){
    //if(verbosity>1)
    printf( "Bspline::fit3D_omp() ns(%i,%i,%i) bPBC=%i dt=%g Ftol=%g nmaxiter=%i bInitGE=%i \n", ns.x,ns.y,ns.z, bPBC, dt, Ftol, nmaxiter, bInitGE );
    const int nxy  = ns.x*ns.y;
    const int nxyz = nxy*ns.z;
    const double F2max = Ftol*Ftol;
    double* ps = new double[nxyz];
    double* fs = new double[nxyz];
    double* vs = new double[nxyz];

    // omp shared
    double vf = 0.0;
    double ff = 0.0;
    double vv = 0.0;
    double err2sum = 0;
    int    itr =0;

    if(bInitGE){ for(int i=0; i<nxyz; i++){ Gs[i]=Es[i]; }; };
    for(int i=0; i<nxyz; i++){ vs[i]=0; };
    //dt = 0.3;

    //int ix=0,iy=0,iz=0;
    int niterdone = 0;
    long t0 = getCPUticks();
    //#pragma omp parallel shared(Gs,fs,ps,vs,itr)
    #pragma omp parallel shared(itr,niterdone,nmaxiter,nxyz,vv,ff,vf,err2sum)
    {
    //for(itr=0; itr<nmaxiter; itr++){
    while(itr<nmaxiter){
        
        #pragma omp single
        { err2sum=0.0; } 
        
        #pragma omp for reduction(+:err2sum)
        for(int iz=0; iz<ns.z; iz++){ 
            if(bPBC){ err2sum += error_iz_pbc( iz, ns, Gs, Es, Ws, ps ); }
            else    { err2sum += error_iz    ( iz, ns, Gs, Es, Ws, ps ); }            
        }
        #pragma omp for
        for(int iz=0; iz<ns.z; iz++){ 
            if(bPBC){ force_iz_pbc( iz, ns, ps, fs ); }
            else    { force_iz    ( iz, ns, ps, fs );     }
        }
        #pragma omp single
        { vf=0; ff=0; vv=0; }
        #pragma omp for reduction(+:vf,ff,vv)
        for(int i=0; i<nxyz; i++){
            ff += fs[i]*fs[i];
            vf += vs[i]*fs[i];
            vv += vs[i]*vs[i];
        }
        #pragma omp single
        { 
        if(vf<0){ for(int i=0; i<nxyz; i++){ vs[i]=0; }; }
        }
        #pragma omp for
        for(int i=0; i<nxyz; i++){
            vs[i] += fs[i]*dt;
            Gs[i] += vs[i]*dt;
        }
        #pragma omp single
        {
            //cfv = move_GD( dt, nxyz, Gs, fs );
            //if(verbosity>2)
            printf( "|F[%i]|=%g Error=%g \n",itr,sqrt(ff), sqrt(err2sum) );
            ///printf( "|F[%i]|=%g cos(f,v)=%g Error=%g \n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z), sqrt(err) );
            if(ff<F2max){ 
                printf( "Bspline::fit3D_omp() DONE |F[%i]|=%g Error=%g \n",itr,sqrt(ff), sqrt(err2sum) );   
                niterdone=itr;
                itr=nmaxiter+1; 
            };
            itr++;
        }
    }
    }
    //if(verbosity>1)
     
    double t = getCPUticks()-t0; printf( "Bspline::fit3D_omp() niter=%i nxyz=%i time %g[GTicks] %g[tick/(nxyz*iter)] \n", niterdone, t*1e-9, t/(nxyz*niterdone) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return niterdone;
}

} // namespace SplineBcub{

#endif



