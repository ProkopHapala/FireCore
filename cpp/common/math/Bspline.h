
#ifndef  Bspline_h
#define  Bspline_h

#include "Vec3.h"
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


void make_inds_pbc( const int n, Quat4i* iqs ){
    iqs[0]={0,1  ,2  ,3  };
    iqs[1]={0,1  ,2  ,3-n};
    iqs[2]={0,1  ,2-n,3-n};
    iqs[3]={0,1-n,2-n,3-n};
}

inline Quat4i choose_inds_pbc( const int i, const int n, const Quat4i* iqs ){
    if(i>=(n-3))[[unlikely]]{ 
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        return iqs[ii]; 
    }
    return Quat4i{ 0, +1, +2, +3 };
}

inline Quat4i choose_inds_pbc_3( const int i, const int n, const Quat4i* iqs ){
    if(i>=(n-3))[[unlikely]]{ 
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        const Quat4i& d = iqs[ii];
        return Quat4i{ i+d.x, i+d.y, i+d.z, i+d.w }   ; 
    }
    return Quat4i{ i, i+1, i+2, i+3 };
}

void make_inds_pbc_5(const int n, Vec6i* iqs) {
    iqs[0] = {0, 1,   2,   3,   4,   5  };
    iqs[1] = {0, 1,   2,   3,   4,   5-n};
    iqs[2] = {0, 1,   2,   3,   4-n, 5-n};
    iqs[3] = {0, 1,   2,   3-n, 4-n, 5-n};
    iqs[4] = {0, 1,   2-n, 3-n, 4-n, 5-n};
    iqs[5] = {0, 1-n, 2-n, 3-n, 4-n, 5-n};
}

inline Vec6i choose_inds_pbc_5(const int i, const int n, const Vec6i* iqs) {
    if (i >= (n - 5)) [[unlikely]] {
        const int ii = i+6-n;
        //return iqs[i+6-n];
        const Vec6i& d = iqs[ii];
        return Vec6i{ i+d.a, i+d.b, i+d.c, i+d.d, i+d.e, i+d.f }   ; 
    }
    return Vec6i{ i, i+1, i+2, i+3, i+4, i+5 };
}

inline int modulo( const int i, const int m ){
    int result = i % m;
    if (result < 0) {
        result += m;
    }
    return result;
}


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
constexpr Vec6T<T> basis5(T t){
    constexpr T inv6 = 1./6.;
    const T t2 = t*t;
    const T t3 = t2*t;
    const T t4 = t2*t2;
    const T t5 = t3*t2;
    return Vec6T<T>{
        //  0.008333333333333333*t5,                                     
        // -0.041666666666666666*t5  +0.041666666666666666*t4  +0.08333333333333333*t3 +0.08333333333333333*t2  +0.041666666666666666*t   +0.008333333333333333, 
        //  0.083333333333333333*t5  -0.166666666666666666*t4  -0.16666666666666666*t3 +0.16666666666666666*t2  +0.416666666666666666*t   +0.216666666666666666,
        // -0.083333333333333333*t5  +0.250000000000000000*t4                          -0.50000000000000000*t2                            +0.550000000000000000,
        //  0.041666666666666666*t5  -0.166666666666666666*t4  +0.16666666666666666*t3 +0.16666666666666666*t2  -0.416666666666666666*t   +0.216666666666666666,
                                                  
        -0.008333333333333333*t5  +0.041666666666666666*t4  -0.08333333333333333*t3 +0.08333333333333333*t2  -0.041666666666666666*t   +0.008333333333333333,
         0.041666666666666666*t5  -0.166666666666666666*t4  +0.16666666666666666*t3 +0.16666666666666666*t2  -0.416666666666666666*t   +0.216666666666666666,        
        -0.083333333333333333*t5  +0.250000000000000000*t4                          -0.50000000000000000*t2                            +0.550000000000000000,  
         0.083333333333333333*t5  -0.166666666666666666*t4  -0.16666666666666666*t3 +0.16666666666666666*t2  +0.416666666666666666*t   +0.216666666666666666,
        -0.041666666666666666*t5  +0.041666666666666666*t4  +0.08333333333333333*t3 +0.08333333333333333*t2  +0.041666666666666666*t   +0.008333333333333333, 
         0.008333333333333333*t5
    };
}

template<typename T>
constexpr inline Vec6T<T> dbasis5(T t){
    constexpr T inv6 = 1./6.;
    const T t2 = t*t;
    const T t3 = t2*t;
    const T t4 = t2*t2;
    return Vec6T<T>{
        // -0.008333333333333333*5*t4  +0.041666666666666666*4*t3  -0.08333333333333333*3*t2 +0.08333333333333333*2*t  -0.041666666666666666 ,
        //  0.041666666666666666*5*t4  -0.166666666666666666*4*t3  +0.16666666666666666*3*t2 +0.16666666666666666*2*t  -0.416666666666666666 ,
        // -0.083333333333333333*5*t4  +0.250000000000000000*4*t3                            -0.50000000000000000*2*t                        ,
        //  0.083333333333333333*5*t4  -0.166666666666666666*4*t3  -0.16666666666666666*3*t2 +0.16666666666666666*2*t  +0.416666666666666666 ,
        // -0.041666666666666666*5*t4  +0.041666666666666666*4*t3  +0.08333333333333333*3*t2 +0.08333333333333333*2*t  +0.041666666666666666 ,   
        //  0.008333333333333333*5*t4               
        -0.0416666666666667*t4	+0.166666666666667*t3	-0.25*t2   +0.166666666666667*t	-0.041666666666666666,	
         0.2083333333333333*t4	-0.666666666666667*t3	+0.50*t2   +0.333333333333333*t -0.416666666666666666,	
        -0.4166666666666667*t4	+1.000000000000000*t3	           -1.000000000000000*t                      ,	
         0.4166666666666667*t4	-0.666666666666667*t3	-0.50*t2   +0.333333333333333*t	+0.416666666666666666,	
        -0.2083333333333333*t4	+0.166666666666667*t3	+0.25*t2   +0.166666666666667*t	+0.041666666666666666,	
         0.0416666666666667*t4
    };
}

template<typename T>
constexpr inline Quat4T<T> basis(T u){
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
constexpr inline Quat4T<T> dbasis(T u){
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
constexpr inline Quat4T<T> ddbasis(T u){
    return Quat4T<T>{
        1-u,
        3*u-2,
        -3*u+1,
        u,
    };
}



inline double project1D_cubic( double w, double x, double x0, double dx, int n, double* ys, Quat4i* xqs ){
    // Convert atomic position to grid position
    double gp = (x-x0)/dx;
    int ix = (int) gp;
    if(gp<0) ix--; // because (int) does not work as floor() for negative numbers
    const double tx = gp - ix;
    const Quat4d bx = Bspline::basis(tx);
    ix=modulo(ix-1,n); 
    const Quat4i xq = choose_inds_pbc_3(ix, n, xqs );
    //printf( "project1D_cubic() ix(%3i) tx(%8.4f) x(%8.4f) x0(%8.4f) dx(%8.4f) \n", ix, tx, x,  x0,dx );
    double sum = 0;
    for (int dx = 0; dx < 4; dx++) {
        const int gx = xq.array[dx];
        double    yi = bx.array[dx] * w;
        //printf( "project1D_cubic()[%i] gx(%3i) yi(%10.5f) \n", dx, gx, yi );
        ys[gx]   += yi;
        sum += yi;
    }
    return sum;
}

inline double project1D_quintic( double w, double x, double x0, double dx, int n, double* ys, Vec6i* xqs ){
    // Convert atomic position to grid position
    double gp = (x-x0)/dx;
    int ix = (int) gp;
    if(gp<0) ix--; // because (int) does not work as floor() for negative numbers
    const double tx = gp - ix;
    const Vec6d bx = Bspline::basis5(tx);
    ix=modulo(ix-2,n); 
    const Vec6i xq = choose_inds_pbc_5(ix, n, xqs );
    double sum = 0;
    for (int dx = 0; dx < 6; dx++) {
        const int gx = xq.array[dx];
        double yi = bx.array[dx] * w;
        ys[gx]   += yi;
        sum += yi;
    }
    return sum;
}

void project2D_cubic( double w, const Vec2d pi, const Vec2d g0, const Vec2d inv_dg, const Vec2i n, double* ys, Quat4i* xqs, Quat4i* yqs ){
    Vec2d gp = (pi-g0)*inv_dg;
    int ix = (int) gp.x;     if(gp.x<0) ix--;
    int iy = (int) gp.y;     if(gp.y<0) iy--;
    const double tx = gp.x - ix;
    const double ty = gp.y - iy;
    const Quat4d bx = Bspline::basis(tx);
    const Quat4d by = Bspline::basis(ty);
    const int nxy = n.x * n.y;
    ix=modulo(ix-1,n.x); const Quat4i xq = choose_inds_pbc_3(ix, n.x, xqs );  // Assuming you pre-calculate xqs, yqs, zqs
    iy=modulo(iy-1,n.y); const Quat4i yq = choose_inds_pbc_3(iy, n.y, yqs );

    printf( "project2D_cubic() ixy(%i,%i) xqs(%i,%i,%i,%i) yqs(%i,%i,%i,%i) nxyz(%i,%i) \n", ix,iy,  xq.x,xq.y,xq.z,xq.w,   yq.x,yq.y,yq.z,yq.w,   n.x,n.y );

    for (int dy = 0; dy < 4; dy++) {
        const int gy  = yq.array[dy];
        const int iiy = gy * n.x;
        const double qbyz = w * by.array[dy];
        for (int dx = 0; dx < 4; dx++) {
            const int gx = xq.array[dx];
            const int ig = gx + iiy;
            ys[ig] += qbyz * bx.array[dx];
        }
    }
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

__attribute__((hot)) 
void sample1D_o5( const double g0, const double dg, const int ng, const double* Gs, const int n, const double* ps, Vec2d* fes ){
    const double inv_dg = 1/dg; 
    for(int i=0; i<n; i++ ){
        const double x  = (ps[i] - g0)*inv_dg;  
        const int    ix = (int)x;
        if( ((ix<2)||(ix>=(ng-3))) )[[unlikely]]{ 
            continue;
            //printf( "ERROR: Bspline::sample1D() ixyz(%i,%i) out of range 0 .. (%i,%i) u(%g,%g)\n", ix,iy, n.x,n.y, u.x,u.y );   exit(0); 
        }
        const double tx = x-ix; 
        const Vec6d p  =  basis5(tx);
        const Vec6d d  = dbasis5(tx);
        const Vec6d gs = *(Vec6d*)(Gs+ix-2);
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
inline Vec3d fe2d_pbc( int nz, const double* E, Quat4i di, const Quat4d& pz, const Quat4d& dz, const Quat4d& by, const Quat4d& dy ){
    alignas(32) Vec2d fe0,fe1,fe2,fe3;
    { alignas(32) const Quat4d cs = *(Quat4d*)(E+di.x); fe0.x=pz.dot( cs ); fe0.y=dz.dot( cs ); } // z-stride 0
    { alignas(32) const Quat4d cs = *(Quat4d*)(E+di.y); fe1.x=pz.dot( cs ); fe1.y=dz.dot( cs ); } // z-stride 1
    { alignas(32) const Quat4d cs = *(Quat4d*)(E+di.z); fe2.x=pz.dot( cs ); fe2.y=dz.dot( cs ); } // z-stride 2
    { alignas(32) const Quat4d cs = *(Quat4d*)(E+di.w); fe3.x=pz.dot( cs ); fe3.y=dz.dot( cs ); } // z-stride 3
    return Vec3d{
        fe0.x*dy.x  +  fe1.x*dy.y  +  fe2.x*dy.z  +  fe3.x*dy.w,  // Fy
        fe0.y*by.x  +  fe1.y*by.y  +  fe2.y*by.z  +  fe3.y*by.w,  // Fz
        fe0.x*by.x  +  fe1.x*by.y  +  fe2.x*by.z  +  fe3.x*by.w   // E
    };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec2d fe1d_pbc_macro( double x, int n, const double* Es, const Quat4i* xqis ){
    int          i = (int)x;
    if(x<0) i--;
    const double t = x - i;
    i=modulo(i-1,n);
    const Quat4i q = choose_inds_pbc( i, n, xqis );
    const Quat4d b =  basis( t );
    const Quat4d d = dbasis( t );
    alignas(32) const Quat4d cs = {Es[q.x],Es[q.y],Es[q.z],Es[q.w]};
    return Vec2d{
        b.dot( cs ),
        d.dot( cs )
    };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec2d fe1d_pbc_imacro( double i, int n, const double* Es, const Quat4i* xqis, Quat4d b, Quat4d d ){
    //int          i = (int)x;  if(x<0) i--;
    //const double t = x - i;
    i=modulo(i-1,n);
    const Quat4i q = choose_inds_pbc( i, n, xqis );
    //const Quat4d b =  basis( t );
    //const Quat4d d = dbasis( t );
    alignas(32) const Quat4d cs = {Es[q.x],Es[q.y],Es[q.z],Es[q.w]};
    return Vec2d{
        b.dot( cs ),
        d.dot( cs )
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


__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d_pbc( const Vec3d u, const Vec3i n, const double* Es, const Quat4i* xqis, const Quat4i* yqis ){
	int          ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    if(u.x<0) ix--;
    if(u.y<0) iy--;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;

    // ---- boundary conditions
    //if(  ((iz<1)||(iz>=n.z-2))  )[[unlikely]]{  
    if(  ((iz<2)||(iz>=n.z-3))  )[[unlikely]]{ 
        //printf( "Bspline::fe3d_pbc_comb3() iz=%i n.z=%i  ixy(%i,%i) \n", iz, n.z, ix, iy );    
        return Quat4dZero; 
    }

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
    //int i0 = (iz-1) + n.z*( iy + n.y*ix); 
    int i0 = iz + n.z*( iy + n.y*ix);  
    const Vec3d E1 = fe2d_pbc( n.z, Es+(i0+qx.x ), qy, bz, dz, by, dy );
    const Vec3d E2 = fe2d_pbc( n.z, Es+(i0+qx.y ), qy, bz, dz, by, dy );
    const Vec3d E3 = fe2d_pbc( n.z, Es+(i0+qx.z ), qy, bz, dz, by, dy );;
    const Vec3d E4 = fe2d_pbc( n.z, Es+(i0+qx.w ), qy, bz, dz, by, dy );
    const Quat4d bx =  basis( tx );
    const Quat4d dx = dbasis( tx );
    return Quat4d{
        dx.dot( {E1.z, E2.z, E3.z, E4.z} ), // Fx
        bx.dot( {E1.x, E2.x, E3.x, E4.x} ), // Fy
        bx.dot( {E1.y, E2.y, E3.y, E4.y} ), // Fz
        bx.dot( {E1.z, E2.z, E3.z, E4.z} ), // E
    };
} 



__attribute__((pure))
__attribute__((hot)) 
Vec3d fe2d_pbc_macro( const Vec2d u, const Vec2i n, const double* Es, const Quat4i* xqis, const Quat4i* yqis ){
	int          ix = (int)u.x  ,  iy = (int)u.y;
    if(u.x<0) ix--;
    if(u.y<0) iy--;
    const double tx = u.x - ix  ,  ty = u.y - iy;

    ix=modulo(ix-1,n.x);
    iy=modulo(iy-1,n.y);

    const Quat4i qx = choose_inds_pbc( ix, n.x, xqis );
    const Quat4i qy = choose_inds_pbc( iy, n.y, yqis )*n.x;

    const Quat4d bx =  basis( tx );
    const Quat4d dx = dbasis( tx );
    const Quat4d by =  basis( ty );
    const Quat4d dy = dbasis( ty );
    int i0 = iy*n.x + ix;  
    return fe2d_pbc( n.x, Es+i0, qy, bx, dx, by, dy );
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
    //if(  ((iz<1)||(iz>=n.z-2))  )[[unlikely]]{  
    if(  ((iz<2)||(iz>=n.z-3))  )[[unlikely]]{ 
        //printf( "Bspline::fe3d_pbc_comb3() iz=%i n.z=%i  ixy(%i,%i) \n", iz, n.z, ix, iy );    
        return Quat4dZero; 
    }
    //printf( "Bspline::fe3d_pbc_comb3() ixyz(%5i,%5i,%5i) @Es=%li \n", ix,iy,iz, (long)Es );  

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
    //int i0 = (iz-1) + n.z*( iy + n.y*ix); 
    int i0 = iz + n.z*( iy + n.y*ix);  
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

} // namespace SplineBcub{

#endif



