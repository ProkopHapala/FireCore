#ifndef InterpolateTricubic_h
#define InterpolateTricubic_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "spline_hermite.h"

#include <immintrin.h>

inline int fold_cubic( int i, int n ){
    if( i==0 )[[unlikely]]{  return n-1; } else if (i>n)[[unlikely]]{ return i-n; } else [[likely]] {  return i-1; };
}


namespace Spline_Hermite{

template<typename T>
inline Quat4T<T> basis_val( T x  ){
    //         x3     x2   x      1  
    //  p-1   -0.5   1.0  -0.5  -0.0
    //  p+0    1.5  -2.5   0.0   1.0  
    //  p+1   -1.5   2.0   0.5   0.0 
    //  p+2    0.5  -0.5   0.0   0.0 
    const double x2 = x*x;
    const double x3 = x*x2;
    return Quat4T<T>{
	   -0.5*x3 +  1.0*x2 +  -0.5*x,  
        1.5*x3 + -2.5*x2 +   0.0*x + 1.0, 
	   -1.5*x3 +  2.0*x2 +   0.5*x,  
        0.5*x3 + -0.5*x2 +   0.0*x
    };
}

template <class T>
inline Quat4T<T> dbasis_val( T x ){
    //  p-1  -1.5   2.0  -0.5   
    //  p+0   4.5  -5.0   0.0   
    //  p+1  -4.5   4.0   0.5   
    //  p+2   1.5  -1.0   0.0  
    const double x2 = x*x;
    return Quat4T<T>{
       -1.5*x2 +  2.0*x + -0.5,  
        4.5*x2 + -5.0*x,  
       -4.5*x2 +  4.0*x +  0.5,  
        1.5*x2 + -1.0*x
    };
}

template <class T>
inline Quat4T<T> dbasis_val_mod( T x ){
//       x3   x2   x   1
//  ----------------------
//  y0   2   -3        1
//  y1  -2   +3
// dy0   1   -2    1
// dy1   1   -1
//       x3   x2   x   1
//  ----------------------
//  y0   6   -6        0
//  y1  -6   +6
// dy0   3   -4    1
// dy1   3   -2

    // T K    =  3*x*(x - 1);
	// c0        =  2*K        ;   //    6*x2 - 6*x
	// c1        = -2*K        ;   //   -6*x2 + 6*x
	// d0        =    K - x + 1;   //    3*x2 - 4*x + 1
	// d1        =    K + x    ;   //    3*x2 - 2*x
    const double f  = 0.5;

    // ToDo: Derivative should consider convexity/concavity - For example assume it is parabola passing the 3 points

    const double x2 = x*x;
    //const double x3 = x2*x;
    const double dy0 = 3*x2 - 4*x + 1;
    const double dy1 = 3*x2 - 2*x;
    return Quat4T<T>{
        0                -f*dy0   ,    // -0.5*dy0  
        6*x2 + -6*x      -f*dy1   ,    // -0.5*dy1
       -6*x2 +  6*x      +f*dy0   ,    // +0.5*dy0
        0                +f*dy1   ,    // +0.5*dy0
    };

    // const double x2 = x*x;
    // return Quat4T<T>{
    //    (-1.5*x2 +  2.0*x + -0.5 )*1.0,  
    //    ( 4.5*x2 + -5.0*x)*1.1,  
    //    (-4.5*x2 +  4.0*x +  0.5)*1.1,  
    //     (1.5*x2 + -1.0*x)*1.0
    // };
}

template <class T>
inline Quat4T<T> ddbasis_val( T x ){
    //  p-1   -3.0   2.  
    //  p+0   +9.0  -5.   
    //  p+1   -9.0   4.  
    //  p+2    3.0  -1.  
    return Quat4T<T>{
       -3.0*x +  2.0,  
        9.0*x + -5.0,  
       -9.0*x +  4.0,  
        3.0*x + -1.0
    };

}

template <class T>
inline Quat4T<T> basis( T x ){
	T x2   = x*x;
	T K    =  x2*(x - 1);
    return Quat4T<T>{
	  2*K - x2 + 1,   //  c0 =  2*x3 - 3*x2 + 1
	 -2*K + x2    ,   //  c1 = -2*x3 + 3*x2
	    K - x2 + x,   //  d0 =    x3 - 2*x2 + x
	    K         ,   //  d1 =    x3 -   x2
    };
}

template <class T>
inline Quat4T<T> dbasis( T x ){
	T K    =  3*x*(x - 1);
    return Quat4T<T>{
	  2*K        ,   //  c0 =  6*x2 - 6*x
	 -2*K        ,   //  c1 = -6*x2 + 6*x
	    K - x + 1,   //  d0 =  3*x2 - 4*x + 1
	    K + x    ,   //  d1 =  3*x2 - 2*x
    };
}

template <class T>
inline Quat4T<T> ddbasis( T x ){
	T x6   =  6*x;
    return Quat4T<T>{
    //  x3     x2    x  1
	    x6 + x6 -  6,  //  c0 =  12*x - 6
	    6  - x6 - x6,  //  c1 = -12*x + 6
	    x6 -  4,      //  d0 =   6*x - 4
	    x6 -  2,      //  d1 =   6*x - 2
    };
}


template<typename T>
inline Quat4T<T> basis_val_old( T x  ){
	const T x2 = x*x;
	const T K  =  x2*(x - 1);
	const T d0 =    K - x2 + x;       //      x3 - 2*x2 + x
	const T d1 =    K         ;       //      x3 -   x2   
    return Quat4T<T>{
	                d0*-0.5, //  p-1 =      -0.5*d0
     2*K - x2 + 1 + d1*-0.5, //  p+0 = c0 + -0.5*d1
	-2*K + x2     + d0* 0.5,  //  p+1 = c1 + +0.5*d0
                    d1* 0.5   //  p+2 =      +0.5*d1
    };
}

template <class T>
inline Quat4T<T> dbasis_val_old( T x ){
	const T K    =  3*x*(x - 1);
    const T d0   =    K - x + 1;   //    3*x2 - 4*x + 1
	const T d1   =    K + x    ;   //    3*x2 - 2*x
    return Quat4T<T>{
              d0*-0.5, //  p-1 =      -0.5*d0
	  2*K   + d1*-0.5, //  p+0 = c0 + -0.5*d1
	 -2*K   + d0* 0.5, //  p+1 = c1 + +0.5*d0
	          d1* 0.5  //  p+2 =      +0.5*d1
    };
}

template <class T>
inline Quat4T<T> ddbasis_val_old( T x ){
	const T x6  =  6*x;
    const T d0  =  x6 -  4;        //     6*x - 4
	const T d1  =  x6 -  2;        //     6*x - 2
    return Quat4T<T>{
	                  d0*-0.5,  //  p-1 =      -0.5*d0
       x6 + x6 -  6 + d1*-0.5,  //  p+0 = c0 + -0.5*d1
	    6 - x6 - x6 + d0* 0.5,  //  p+1 = c1 + +0.5*d0
	                  d1* 0.5   //  p+2 =      +0.5*d1
    };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec2d fe2d_deriv( const double tx, const double ty, const Vec2i i, const Quat4d* FEs, const Quat4d* dFs, Quat4d& out ){
    alignas(32) Quat4d fe1,fe2;
    //alignas(32) Quat4d dfe;
    double dxy1,dxy2,dzy1,dzy2;  //   dxy1 = dE/dxdy(1),  
    double dxyz1,dxyz2;
    {
        const Quat4d bx =  basis( tx ); // 4 Hermite basis functions along x-axis
        const Quat4d dx = dbasis( tx ); // 4 derivatives of Hermite basis functions along x-axis
        { // x-line y1
            alignas(32) const Quat4d p1 = FEs[i.x  ];  // (dx,dy,dz,E)       Energy and forces at point_1
            alignas(32) const Quat4d d1 = dFs[i.x  ];  // (dyz,dxy,dxz,dxyz) 2nd derivatives at point_1
            alignas(32) const Quat4d p2 = FEs[i.x+1];  // (dx,dy,dz,E)       Energy and forces at point_1
            alignas(32) const Quat4d d2 = dFs[i.x+1];  // (dyz,dxy,dxz,dxyz) 2nd derivatives at point_2
            
            fe1.e  = bx.x*p1.e +  bx.y*p2.e + bx.z*p1.x + bx.w*p2.x;  // E       (x,0)
            fe1.x  = dx.x*p1.e +  dx.y*p2.e + dx.z*p1.x + dx.w*p2.x;  // (dE/dx) (x,0)
            fe1.y  = bx.x*p1.y +  bx.y*p2.y + bx.z*d1.y + bx.w*d2.y;  // (dE/dy) (x,0)
            fe1.z  = bx.x*p1.z +  bx.y*p2.z + bx.z*d1.z + bx.w*d2.z;  // (dE/dz) (x,0)

            dxy1   = dx.x*p1.y +  dx.y*p2.y + dx.z*d1.y  + dx.w*d2.y; // (d2E/dxdy)  (x,0)
            dzy1   = bx.x*d1.x +  bx.y*d2.x + bx.z*d1.e  + bx.w*d2.e; // (d2E/dzdy)  (x,0)
            dxyz1  = dx.x*d1.x +  dx.y*d2.x + dx.z*d1.e  + dx.w*d2.e; // (d3E/dxdydz)(x,0)
        }
        { // x-line y2
            alignas(32) const Quat4d p1 = FEs[i.y  ]; 
            alignas(32) const Quat4d d1 = dFs[i.y  ];  // (dxy,dxz,dyz,dxyz)
            alignas(32) const Quat4d p2 = FEs[i.y+1]; 
            alignas(32) const Quat4d d2 = dFs[i.y+1]; 

            fe2.e  = bx.x*p1.e +  bx.y*p2.e + bx.z*p1.x  + bx.w*p2.x;  // E       (x,1)
            fe2.x  = dx.x*p1.e +  dx.y*p2.e + dx.z*p1.x  + dx.w*p2.x;  // (dE/dx) (x,1)
            fe2.y  = bx.x*p1.y +  bx.y*p2.y + bx.z*d1.y + bx.w*d2.y;   // (dE/dy) (x,1)
            fe2.z  = bx.x*p1.z +  bx.y*p2.z + bx.z*d1.z + bx.w*d2.z;   // (dE/dz) (x,1)

            dxy2  = dx.x*p1.y +  dx.y*p2.y + dx.z*d1.y  + dx.w*d2.y;   // (d2E/dxdy)   (x,1)
            dzy2  = bx.x*d1.x +  bx.y*d2.x + bx.z*d1.e  + bx.w*d2.e;   // (d2E/dxdy)   (x,1)
            dxyz2 = dx.x*d1.x +  dx.y*d2.x + dx.z*d1.e  + dx.w*d2.e;   // (d3E/dxdydz) (x,0)
        }
    }
    alignas(32) const Quat4d by =  basis( ty );
    alignas(32) const Quat4d dy = dbasis( ty );
    out.w = by.x*fe1.e +  by.y*fe2.e + by.z*fe1.y + by.w*fe2.y; // E       (x,y)
    out.y = dy.x*fe1.e +  dy.y*fe2.e + dy.z*fe1.y + dy.w*fe2.y; // (dE/dy) (x,y)
    out.x = by.x*fe1.x +  by.y*fe2.x + by.z*dxy1  + by.w*dxy2;  // (dE/dx) (x,y)
    out.z = by.x*fe1.z +  by.y*fe2.z + by.z*dzy1  + by.w*dzy2;  // (dE/dz) (x,y)
    return Vec2d{
        by.x*dzy1  + by.y*dzy2  + by.z*dxyz1 + by.w*dxyz2,      // (d2E/dxdz) (x,y)
        dy.x*fe1.z + dy.y*fe2.z + dy.z*dzy1  + dy.w*dzy2,       // (d2E/dydz) (x,y)
    };
}

__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d_deriv( const Vec3d u, const Vec3i n, const Quat4d* FEs, const Quat4d* dFs  ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    if( 
        ((ix<0)||(ix>=n.x-3)) ||
        ((iy<0)||(iy>=n.y-3)) ||
        ((iz<0)||(iz>=n.z-3))        
    )[[unlikely]]{ printf( "ERROR: Spline_Hermite::fe3d() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }
    const int nxy = n.x*n.y;
    Quat4d fe1,fe2;
    const int i0  = ix + n.x*(iy + n.y*iz ); const Vec2d dfe1 = fe2d_deriv(tx,ty, {i0,i0+n.x}, FEs, dFs, fe1 );
    const int i1  = i0+nxy;                  const Vec2d dfe2 = fe2d_deriv(tx,ty, {i1,i1+n.x}, FEs, dFs, fe2 );
    const Quat4d bz =  basis( tz );
    const Quat4d dz = dbasis( tz );
    return Quat4d{
        bz.dot( {fe1.x, fe2.x, dfe1.x, dfe2.x} ), // Fx
        bz.dot( {fe1.y, fe2.y, dfe1.y, dfe2.y} ), // Fy
        dz.dot( {fe1.e, fe2.e,  fe1.z,  fe2.z} ), // Fz
        bz.dot( {fe1.e, fe2.e,  fe1.z,  fe2.z} ), // E
    };
} 

__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d_v4( const double tx, const double ty, const Quat4i i, const Quat4d PLQH, const Quat4d* V ){
    alignas(32) Quat4d e,fx;
    {
        const Quat4d bx =  basis_val( tx );
        const Quat4d dx = dbasis_val( tx );
        {
            alignas(32) Quat4d p;
            const Quat4d* Vi = V+i.x;
            p.x = PLQH.dot( *(Vi  ) );
            p.y = PLQH.dot( *(Vi+1) );
            p.z = PLQH.dot( *(Vi+2) );
            p.w = PLQH.dot( *(Vi+3) );
            e.x  = bx.dot(p); 
            fx.x = dx.dot(p);
        }
        {
            alignas(32) Quat4d p;
            const Quat4d* Vi = V+i.y;
            p.x = PLQH.dot( *(Vi  ) );
            p.y = PLQH.dot( *(Vi+1) );
            p.z = PLQH.dot( *(Vi+2) );
            p.w = PLQH.dot( *(Vi+3) );
            e.y  = bx.dot(p);
            fx.y = dx.dot(p);
        }
        {
            alignas(32) Quat4d p;
            const Quat4d* Vi = V+i.z;
            p.x = PLQH.dot( *(Vi  ) );
            p.y = PLQH.dot( *(Vi+1) );
            p.z = PLQH.dot( *(Vi+2) );
            p.w = PLQH.dot( *(Vi+3) );
            e.z  = bx.dot(p);
            fx.z = dx.dot(p);
        }
        {
            alignas(32) Quat4d p;
            const Quat4d* Vi = V+i.w;
            p.x = PLQH.dot( *(Vi  ) );
            p.y = PLQH.dot( *(Vi+1) );
            p.z = PLQH.dot( *(Vi+2) );
            p.w = PLQH.dot( *(Vi+3) );
            e.w  = bx.dot(p);
            fx.w = dx.dot(p);
        }
    }
    alignas(32) const Quat4d by =  basis_val( ty );
    alignas(32) const Quat4d dy = dbasis_val( ty );
    return Vec3d{
        by.dot(fx), // Fx
        dy.dot(e ), // Fy
        by.dot(e )  // E
    };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d( const double tx, const double ty, const Quat4i i, const double* Es ){
    alignas(32) Quat4d e,fx;
    {
        const Quat4d bx =  basis_val( tx );
        const Quat4d dx = dbasis_val( tx );
        {
            alignas(32) const Quat4d p = *(Quat4d*)(Es+i.x); // read 4 doubles from global memory at a time ( 4*8 = 32 bytes = 256 bits ) ideal for SIMD AVX2
            e.x  = bx.dot(p);   // not sure how dot() is SIMD optimized => maybe we should flip the order of x and y strides ?
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
    alignas(32) const Quat4d by =  basis_val( ty );
    alignas(32) const Quat4d dy = dbasis_val( ty );
    return Vec3d{
        by.dot(fx), // Fx
        dy.dot(e ), // Fy
        by.dot(e )  // E
    };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d_v2( const double tx, const double ty, const Quat4i i, const double* Es ){
    Quat4d e,fy;
    {
        alignas(32) const Quat4d by =  basis_val( ty );
        alignas(32) const Quat4d dy = dbasis_val( ty );

        alignas(32) const Quat4d p1 = *(Quat4d*)(Es+i.x);
        alignas(32) const Quat4d p2 = *(Quat4d*)(Es+i.y);
        alignas(32) const Quat4d p3 = *(Quat4d*)(Es+i.z); 
        alignas(32) const Quat4d p4 = *(Quat4d*)(Es+i.w); 

        e.x = by.x*p1.x + by.y*p2.x + by.z*p3.x + by.w*p4.x;
        e.y = by.x*p1.y + by.y*p2.y + by.z*p3.y + by.w*p4.y;
        e.z = by.x*p1.z + by.y*p2.z + by.z*p3.z + by.w*p4.z;
        e.w = by.x*p1.w + by.y*p2.w + by.z*p3.w + by.w*p4.w;

        fy.x = dy.x*p1.x + dy.y*p2.x + dy.z*p3.x + dy.w*p4.x;
        fy.y = dy.x*p1.y + dy.y*p2.y + dy.z*p3.y + dy.w*p4.y;
        fy.z = dy.x*p1.z + dy.y*p2.z + dy.z*p3.z + dy.w*p4.z;
        fy.w = dy.x*p1.w + dy.y*p2.w + dy.z*p3.w + dy.w*p4.w;

    }
    alignas(32) const Quat4d bx =  basis_val( tx );
    alignas(32) const Quat4d dx = dbasis_val( tx );
    return Vec3d{
        dx.dot(e ), // Fx
        bx.dot(fy), // Fy
        bx.dot(e )  // E
    };
}


__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d_v3( const Quat4d bx, const Quat4d dx, const Quat4d by, const Quat4d dy, const Quat4i i, const double* Es ){
    Quat4d e,fy;
    {
        alignas(32) const Quat4d p1 = *(Quat4d*)(Es+i.x);
        alignas(32) const Quat4d p2 = *(Quat4d*)(Es+i.y);
        alignas(32) const Quat4d p3 = *(Quat4d*)(Es+i.z); 
        alignas(32) const Quat4d p4 = *(Quat4d*)(Es+i.w); 

        e.x = by.x*p1.x + by.y*p2.x + by.z*p3.x + by.w*p4.x;
        e.y = by.x*p1.y + by.y*p2.y + by.z*p3.y + by.w*p4.y;
        e.z = by.x*p1.z + by.y*p2.z + by.z*p3.z + by.w*p4.z;
        e.w = by.x*p1.w + by.y*p2.w + by.z*p3.w + by.w*p4.w;

        fy.x = dy.x*p1.x + dy.y*p2.x + dy.z*p3.x + dy.w*p4.x;
        fy.y = dy.x*p1.y + dy.y*p2.y + dy.z*p3.y + dy.w*p4.y;
        fy.z = dy.x*p1.z + dy.y*p2.z + dy.z*p3.z + dy.w*p4.z;
        fy.w = dy.x*p1.w + dy.y*p2.w + dy.z*p3.w + dy.w*p4.w;

    }
    return Vec3d{
        dx.dot(e ), // Fx
        bx.dot(fy), // Fy
        bx.dot(e )  // E
    };
}



__attribute__((pure))
__attribute__((hot)) 
Vec3f fe2f( const float tx, const float ty, const Quat4i i, const float* Es ){
    alignas(32) Quat4f e,fx;
    {
        const Quat4f bx =  basis_val( tx );
        const Quat4f dx = dbasis_val( tx );
        {
            alignas(32) const Quat4f p = *(Quat4f*)(Es+i.x); // read 4 doubles from global memory at a time ( 4*8 = 32 bytes = 256 bits ) ideal for SIMD AVX2
            e.x  = bx.dot(p);   // not sure how dot() is SIMD optimized => maybe we should flip the order of x and y strides ?
            fx.x = dx.dot(p);
        }
        {
            alignas(32) const Quat4f p = *(Quat4f*)(Es+i.y); 
            e.y  = bx.dot(p);
            fx.y = dx.dot(p);
        }
        {
            alignas(32) const Quat4f p = *(Quat4f*)(Es+i.z); 
            e.z  = bx.dot(p);
            fx.z = dx.dot(p);
        }
        {
            alignas(32) const Quat4f p = *(Quat4f*)(Es+i.w); 
            e.w  = bx.dot(p);
            fx.w = dx.dot(p);
        }
    }
    alignas(32) const Quat4f by =  basis_val( ty );
    alignas(32) const Quat4f dy = dbasis_val( ty );
    return Vec3f{
        by.dot(fx), // Fx
        dy.dot(e ), // Fy
        by.dot(e )  // E
    };
}

__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d_v2( const Vec3d u, const Vec3i n, const double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    if( 
        ((ix<0)||(ix>=n.x-3)) ||
        ((iy<0)||(iy>=n.y-3)) ||
        ((iz<0)||(iz>=n.z-3))        
    )[[unlikely]]{ printf( "ERROR: Spline_Hermite::fe3d_v2() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }

    alignas(32) const Quat4d bx =  basis_val( tx );
    alignas(32) const Quat4d dx = dbasis_val( tx );
    alignas(32) const Quat4d by =  basis_val( ty );
    alignas(32) const Quat4d dy = dbasis_val( ty );
    alignas(32) const Quat4d bz =  basis_val( tz );
    alignas(32) const Quat4d dz = dbasis_val( tz );

    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    int i0 = ix + n.x*( iy + n.y*iz );  const Vec3d Exy1 = fe2d_v3(bx,dx,by,dy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es ); // i0 += nxy;
    int i1 = i0+nxy  ;                  const Vec3d Exy2 = fe2d_v3(bx,dx,by,dy, {i1,i1+n.x,i1+n.x*2,i1+3*n.x}, Es );  //i0 += nxy;
    int i2 = i0+nxy*2;                  const Vec3d Exy3 = fe2d_v3(bx,dx,by,dy, {i2,i2+n.x,i2+n.x*2,i2+3*n.x}, Es );  //i0 += nxy;
    int i3 = i0+nxy*3;                  const Vec3d Exy4 = fe2d_v3(bx,dx,by,dy, {i3,i3+n.x,i3+n.x*2,i3+3*n.x}, Es );
    return Quat4d{
        bz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        bz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        dz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        bz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
        // Exy1.x*bz.x + Exy2.x*bz.y + Exy3.x*bz.z + Exy4.x*bz.w , // Fx
        // Exy1.y*bz.x + Exy2.y*bz.y + Exy3.y*bz.z + Exy4.y*bz.w , // Fy
        // Exy1.z*dz.x + Exy2.z*dz.x + Exy3.z*dz.x + Exy4.z*dz.x , // Fz
        // Exy1.z*bz.x + Exy2.z*bz.x + Exy3.z*bz.x + Exy4.z*bz.x   // E
    };
} 

__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d( const Vec3d u, const Vec3i n, const double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    if( 
        ((ix<0)||(ix>=n.x-3)) ||
        ((iy<0)||(iy>=n.y-3)) ||
        ((iz<0)||(iz>=n.z-3))        
    )[[unlikely]]{ printf( "ERROR: Spline_Hermite::fe3d() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }

    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    //int i0 = ix + n.x*( iy + n.y*iz ); 
    // const Vec3d Exy1 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    // const Vec3d Exy2 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    // const Vec3d Exy3 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    // const Vec3d Exy4 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );

    int i0 = ix + n.x*( iy + n.y*iz );  const Vec3d Exy1 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
    int i1 = i0+nxy  ;                  const Vec3d Exy2 = fe2d(tx,ty, {i1,i1+n.x,i1+n.x*2,i1+3*n.x}, Es );
    int i2 = i0+nxy*2;                  const Vec3d Exy3 = fe2d(tx,ty, {i2,i2+n.x,i2+n.x*2,i2+3*n.x}, Es );
    int i3 = i0+nxy*3;                  const Vec3d Exy4 = fe2d(tx,ty, {i3,i3+n.x,i3+n.x*2,i3+3*n.x}, Es );

    const Quat4d bz =  basis_val( tz );
    const Quat4d dz = dbasis_val( tz );
    return Quat4d{
        bz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        bz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        dz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        bz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
    };
} 

__attribute__((pure))
__attribute__((hot)) 
Quat4d fe3d_v4( const Quat4d PLQH, const Vec3d u, const Vec3i n, const Quat4d* V ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    // if( 
    //     ((ix<0)||(ix>=n.x-3)) ||
    //     ((iy<0)||(iy>=n.y-3)) ||
    //     ((iz<0)||(iz>=n.z-3))        
    // )[[unlikely]]{ printf( "ERROR: Spline_Hermite::fe3d_v4() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }
    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    int i0 = ix + n.x*( iy + n.y*iz );  const Vec3d Exy1 = fe2d_v4(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, PLQH, V );
    int i1 = i0+nxy  ;                  const Vec3d Exy2 = fe2d_v4(tx,ty, {i1,i1+n.x,i1+n.x*2,i1+3*n.x}, PLQH, V );
    int i2 = i0+nxy*2;                  const Vec3d Exy3 = fe2d_v4(tx,ty, {i2,i2+n.x,i2+n.x*2,i2+3*n.x}, PLQH, V );
    int i3 = i0+nxy*3;                  const Vec3d Exy4 = fe2d_v4(tx,ty, {i3,i3+n.x,i3+n.x*2,i3+3*n.x}, PLQH, V );
    const Quat4d bz =  basis_val( tz );
    const Quat4d dz = dbasis_val( tz );
    return Quat4d{
        bz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        bz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        dz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        bz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
    };
    //return Quat4d{ 0.0,0.0,0.0, PLQH.dot( V[ix + n.x*( iy + n.y*iz )]  )   };
} 

__attribute__((pure))
__attribute__((hot)) 
Quat4f fe3f( const Vec3f u, const Vec3i n, const float* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const float tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    const float mx = 1-tx      ,  my = 1-ty      ,  mz = 1-tz      ;

    if( 
        ((ix<0)||(ix>=n.x-3)) ||
        ((iy<0)||(iy>=n.y-3)) ||
        ((iz<0)||(iz>=n.z-3))        
    )[[unlikely]]{ printf( "ERROR: Spline_Hermite::fe3f() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }

    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    int i0 = ix + n.x*( iy + n.y*iz ); 
    const Vec3f Exy1 = fe2f(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3f Exy2 = fe2f(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3f Exy3 = fe2f(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3f Exy4 = fe2f(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
    const Quat4f bz =  basis_val( tz );
    const Quat4f dz = dbasis_val( tz );
    return Quat4f{
        bz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        bz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        dz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        bz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
    };
}


inline Vec2d fe1Dcomb2( const int ix, const Quat4d* FE, const Vec2d& C, const Quat4d& p, const Quat4d& d ){
    //const Quat4d p =  basis(x);
    //const Quat4d d = dbasis(x);
    const Quat4d a = FE[ix  ];
    const Quat4d b = FE[ix+1];
    const Quat4d cs{ a.x*C.x+a.z*C.y, b.x*C.x+b.z*C.y, a.y*C.x+a.w*C.y, b.y*C.x+b.w*C.y };
    return Vec2d{ p.dot( cs ), d.dot( cs ) };
}


__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d_comb2( int iz, int nz, const Quat4d* FE, const Vec2d& C, const Quat4d& pz, const Quat4d& dz, const Quat4d& by, const Quat4d& dy ){
    //const Quat4d* FEx = FE + ( i.x*n.y  + i.y )*n.z;
    alignas(32) const Vec2d fe0 = fe1Dcomb2( iz, FE     , C, pz, dz );
    alignas(32) const Vec2d fe1 = fe1Dcomb2( iz, FE+nz  , C, pz, dz );
    alignas(32) const Vec2d fe2 = fe1Dcomb2( iz, FE+nz*2, C, pz, dz );
    alignas(32) const Vec2d fe3 = fe1Dcomb2( iz, FE+nz*3, C, pz, dz );
    return Vec3d{
        fe0.x*dy.x  +  fe1.x*dy.y  +  fe2.x*dy.z  +  fe3.x*dy.w,  // Fy
        fe0.y*by.x  +  fe1.y*by.y  +  fe2.y*by.z  +  fe3.y*by.w,  // Fz
        fe0.x*by.x  +  fe1.x*by.y  +  fe2.x*by.z  +  fe3.x*by.w   // E
    };
}


__attribute__((pure))
__attribute__((hot)) 
inline Quat4d fe3d_comb2( const double tx, const double ty, const double tz, const Quat4i i, const Quat4i n, const Quat4d* FE, const Vec2d& C ){
    alignas(32) const Quat4d pz  =  basis(tz);
    alignas(32) const Quat4d dz  = dbasis(tz);
    alignas(32) const Quat4d by  =  basis_val( ty );
    alignas(32) const Quat4d dy  = dbasis_val( ty );
    alignas(32) const Quat4d bx  =  basis_val( tx );
    alignas(32) const Quat4d dx  = dbasis_val( tx );

    const int nyz = n.y*n.z;
    const Quat4d* FEx = FE + (i.z + ( (i.y-1) + (i.x-1)*n.y )*n.z);
    alignas(32) const Vec3d fe0 = fe2d_comb2( i.z, n.z, FEx      , C, pz, dz, by, dy );
    alignas(32) const Vec3d fe1 = fe2d_comb2( i.z, n.z, FEx+nyz  , C, pz, dz, by, dy );
    alignas(32) const Vec3d fe2 = fe2d_comb2( i.z, n.z, FEx+nyz*2, C, pz, dz, by, dy );
    alignas(32) const Vec3d fe3 = fe2d_comb2( i.z, n.z, FEx+nyz*3, C, pz, dz, by, dy );

    return Quat4d{
        fe0.z*bx.x  +  fe1.z*bx.y  +  fe2.z*bx.z  +  fe3.z*bx.w, // E
        fe0.x*dx.x  +  fe1.x*dx.y  +  fe2.x*dx.z  +  fe3.x*dx.w, // Fx
        fe0.x*bx.x  +  fe1.x*bx.y  +  fe2.x*bx.z  +  fe3.x*bx.w, // Fy
        fe0.y*bx.x  +  fe1.y*bx.y  +  fe2.y*bx.z  +  fe3.y*bx.w  // Fz  
    };
}

//__attribute__((always_inline))
__attribute__((pure))
__attribute__((hot)) 
inline Vec2d fe1Dcomb3( const Vec2d* FE, const Vec3d& C, const Quat4d& p, const Quat4d& d ){
    //const Quat4d p =  basis(x);
    //const Quat4d d = dbasis(x);
    alignas(32) const Vec2d a1 = FE[0];
    alignas(32) const Vec2d b1 = FE[1];
    alignas(32) const Vec2d c1 = FE[2];
    alignas(32) const Vec2d a2 = FE[3];
    alignas(32) const Vec2d b2 = FE[4];
    alignas(32) const Vec2d c2 = FE[5];
    alignas(32) const Quat4d cs{ 
        a1.x*C.x + b1.x*C.y + c1.x*C.z,   // y1 
        a2.x*C.x + b2.x*C.y + c2.x*C.z,   // y2
        a1.y*C.x + b1.y*C.y + c1.y*C.z,   // dy1
        a2.y*C.x + b2.y*C.y + c2.y*C.z};  // dy2
    return Vec2d{ p.dot( cs ), d.dot( cs ) };
    //printf( "i %i ", ix );
    //return Vec2d{ a1.x*C.x + b1.x*C.y + c1.x*C.z, 0.0 };
}

__attribute__((pure))
__attribute__((hot)) 
inline Vec3d fe2d_comb3( int nz, const Vec2d* FE, const Vec3d& C, const Quat4d& pz, const Quat4d& dz, const Quat4d& by, const Quat4d& dy ){
    //const Quat4d* FEx = FE + ( i.x*n.y  + i.y )*n.z;
    alignas(32) const Vec2d fe0 = fe1Dcomb3( FE     , C, pz, dz );
    alignas(32) const Vec2d fe1 = fe1Dcomb3( FE+nz*3, C, pz, dz );
    alignas(32) const Vec2d fe2 = fe1Dcomb3( FE+nz*6, C, pz, dz );
    alignas(32) const Vec2d fe3 = fe1Dcomb3( FE+nz*9, C, pz, dz );
    return Vec3d{
        fe0.x*dy.x  +  fe1.x*dy.y  +  fe2.x*dy.z  +  fe3.x*dy.w,  // Fy
        fe0.y*by.x  +  fe1.y*by.y  +  fe2.y*by.z  +  fe3.y*by.w,  // Fz
        fe0.x*by.x  +  fe1.x*by.y  +  fe2.x*by.z  +  fe3.x*by.w   // E
    };
}


__attribute__((pure))
__attribute__((hot)) 
inline Quat4d fe3d_comb3( const Vec3d& t, const Vec3i i, const Vec3i n, const Vec2d* FE, const Vec3d& C ){
    alignas(32) const Quat4d pz  =  basis(t.z);
    alignas(32) const Quat4d dz  = dbasis(t.z);
    alignas(32) const Quat4d by  =  basis_val( t.y );
    //alignas(32) const Quat4d dy  = dbasis_val( t.y );
    alignas(32) const Quat4d dy  = dbasis_val_mod( t.y );
    alignas(32) const Quat4d bx  =  basis_val( t.x );
    //alignas(32) const Quat4d dx  = dbasis_val( t.x );
    alignas(32) const Quat4d dx  = dbasis_val_mod( t.x );

    const int nyz = n.y*n.z*3;
    //const Vec2d* FEx = FE + (i.z + ( (i.y-1) + (i.x-1)*n.y )*n.z)*3;
    const Vec2d* FEx = FE + (i.z + ( i.y + i.x*n.y )*n.z)*3;

    // constant approx
    //const Vec3d fe0 = fe2d_comb3( n.z, FEx      , C, pz, dz, by, dy );
    //return Quat4d{ 0.0, fe0.x, fe0.y, fe0.z };

    // linear approx
    // const Vec3d fe0 = fe2d_comb3( n.z, FEx    , C, pz, dz, by, dy );
    // const Vec3d fe1 = fe2d_comb3( n.z, FEx+nyz, C, pz, dz, by, dy );
    // return Quat4d{ 
    //     0.0, 
    //     fe0.x*(1-t.x)+fe1.x*(t.x), 
    //     fe0.y*(1-t.x)+fe1.y*(t.x),
    //     fe0.z*(1-t.x)+fe1.z*(t.x),
    // };

    alignas(32) const Vec3d fe0 = fe2d_comb3( n.z, FEx      , C, pz, dz, by, dy );
    alignas(32) const Vec3d fe1 = fe2d_comb3( n.z, FEx+nyz  , C, pz, dz, by, dy );
    alignas(32) const Vec3d fe2 = fe2d_comb3( n.z, FEx+nyz*2, C, pz, dz, by, dy );
    alignas(32) const Vec3d fe3 = fe2d_comb3( n.z, FEx+nyz*3, C, pz, dz, by, dy );
    return Quat4d{
        fe0.z*dx.x  +  fe1.z*dx.y  +  fe2.z*dx.z  +  fe3.z*dx.w, // Fx
        fe0.x*bx.x  +  fe1.x*bx.y  +  fe2.x*bx.z  +  fe3.x*bx.w, // Fy
        fe0.y*bx.x  +  fe1.y*bx.y  +  fe2.y*bx.z  +  fe3.y*bx.w, // Fz  
        fe0.z*bx.x  +  fe1.z*bx.y  +  fe2.z*bx.z  +  fe3.z*bx.w  // E
    };
    
}


//#endif WITH_AVX

__attribute__((hot)) 
void sample1D( const double g0, const double dg, const int ng, const double* Eg, const int n, const double* xs, Vec2d* FEs ){
    const double inv_dg = 1/dg;
    for(int i=0; i<n; i++ ){
        const double t = ((xs[i]-g0)*inv_dg);
        const int it   = (int)t;
        //printf( "sample_SplineHermite() x[%i]=%g t=%g it=%i | g0=%g dg=%g \n", i, xs[i], t, it, g0, dg );
        if( (it<0)||(it>=ng-3) )[[unlikely]]{ printf( "ERROR: Spline_Hermite::sample1D it(%i) out of range (0,%i) | xs[%i]=%g -> t=%g \n", it, ng, i, xs[i], t ); exit(0); }
        double dt = t-it;
        const Quat4d bs = basis_val ( dt );
        const Quat4d ds = dbasis_val( dt );
        const Quat4d v = *(Quat4d*)(Eg+it-1);
        FEs[i].x = bs.dot(v);
        FEs[i].y = ds.dot(v)*inv_dg;
    }
}


__attribute__((hot)) 
void sample2D( const Vec2d g0, const Vec2d dg, const Vec2i ng, const double* Eg, const int n, const Vec2d* ps, Vec3d* fes ){
    printf( "sample2D() g0=(%g,%g) dg=(%g,%g) ng=(%i,%i) n=%i \n", g0.x,g0.y, dg.x,dg.y, ng.x,ng.y, n );
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        const Vec2d t  = (ps[i] - g0)*inv_dg; 
        const int ix = ((int)t.x);
        const int iy = ((int)t.y);
        if( ((ix<1)||(ix>=ng.x-2)) || ((iy<1)||(iy>=ng.y-2)) )[[unlikely]]{ 
            //printf( "ERROR: Spline_Hermite::sample2D() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); 
            fes[i]=Vec3dZero;
            continue;
        }
        //const int i0 = ix + ng.x*iy;
        //const int i0 = (ix+1) + ng.x*(iy+1);
        const int i0 = (ix-1) + ng.x*(iy-1);
        
        //Vec3d fe = fe2d( t.x-ix,t.y-iy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );        // sample2D(n=10000) time=527.478[kTick] 52.7478[tick/point]
        Vec3d fe = fe2d_v2( t.x-ix,t.y-iy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );   // sample2D(n=10000) time=553.47[kTick] 55.347[tick/point]
        //fe.x*=inv_dg.x;
        //fe.y*=inv_dg.x;

        //Vec3d fe{  0.0,0.0, Eg[ iy*ng.x + ix ] }; // Nearest interpolation
        //Vec3d fe{  0.0,0.0, Eg[ ix*ng.y + iy ] }; // Nearest interpolation
        //Vec3d fe{  0.0,0.0, iy  }; // Nearest interpolation
        //Vec3d fe{  0.0,0.0, ix  }; // Nearest interpolation

        fes[i]=fe;
        //printf( "sample2D()[%i] ps(%g,%g) E=%g Fxy(%g,%g)\n", i, ps[i].x,ps[i].y,  fes[i].z,fes[i].x,fes[i].y );
    }
}

__attribute__((hot)) 
void sample3D( const Vec3d g0, const Vec3d dg, const Vec3i ng, const double* Eg, const int n, const Vec3d* ps, Quat4d* fes ){
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        Quat4d fe = fe3d( (ps[i]-g0)*inv_dg, ng, Eg );        // sample3D(n=10000) time=2009.44[kTick] 200.944[tick/point]
        //Quat4d fe = fe3d_v2( (ps[i]-g0)*inv_dg, ng, Eg );   // sample3D(n=10000) time=2175.84[kTick] 217.584[tick/point]
        fe.f.mul(inv_dg);
        fes[i] = fe;
    }
}

__attribute__((hot)) 
void sample2D_deriv( const Vec2d g0, const Vec2d dg, const Vec2i ng, const Quat4d* Eg, const Quat4d* dEg, const int n, const Vec2d* ps, Quat4d* fes ){
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        const Vec2d t  = (ps[i] - g0)*inv_dg; 
        const int ix = (int)t.x;
        const int iy = (int)t.y;
        //if( ((ix<0)||(ix>=ng.x-3)) || ((iy<0)||(iy>=ng.y-3)) )[[unlikely]]{ printf( "ERROR: Spline_Hermite::sample2D() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); }
        const int i0 = ix + ng.x*iy;
        Quat4d fe;
        fe2d_deriv( t.x-ix,t.y-iy, {i0,i0+ng.x}, Eg, dEg, fe );        // sample2D(n=10000) time=527.478[kTick] 52.7478[tick/point]
        //Vec3d fe = fe2d_v2( t.x-ix,t.y-iy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );   // sample2D(n=10000) time=553.47[kTick] 55.347[tick/point]
        fe.x*=inv_dg.x;
        fe.y*=inv_dg.x;
        fes[i]=fe;
        //printf( "sample2D()[%i] ps(%g,%g) E=%g Fxy(%g,%g)\n", i, ps[i].x,ps[i].y,  fes[i].z,fes[i].x,fes[i].y );
    }
}

__attribute__((hot)) 
void sample3D_deriv( const Vec3d g0, const Vec3d dg, const Vec3i ng, const Quat4d* Eg, const Quat4d* dEg, const int n, const Vec3d* ps, Quat4d* fes ){
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        Quat4d fe = fe3d_deriv( (ps[i]-g0)*inv_dg, ng, Eg, dEg );        // sample3D(n=10000) time=2009.44[kTick] 200.944[tick/point]
        //Quat4d fe = fe3d_v2( (ps[i]-g0)*inv_dg, ng, Eg );   // sample3D(n=10000) time=2175.84[kTick] 217.584[tick/point]
        fe.f.mul(inv_dg);
        fes[i] = fe;
    }
}

__attribute__((hot)) 
void sample3D( const Vec3f g0, const Vec3f dg, const Vec3i ng, const float* Eg, const int n, const Vec3f* ps, Quat4f* fes ){
    Vec3f inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        Quat4f fe = fe3f( (ps[i]-g0)*inv_dg, ng, Eg );
        fe.f.mul(inv_dg);
        fes[i] = fe;
    }
}

//__attribute__((hot)) 
void sample1D_deriv( const double g0, const double dg, const int ng, const Vec2d* FE, const int n, const double* ps, Vec2d* fes ){
    const double inv_dg = 1/dg; 
    for(int i=0; i<n; i++ ){
        const double x = (ps[i] - g0)*inv_dg;  
        const int    ix = (int)x;
        const double tx = x-ix; 
        const Quat4d p =  basis(tx);
        const Quat4d d = dbasis(tx);
        const Vec2d fe1 = FE[ix  ];
        const Vec2d fe2 = FE[ix+1];
        const Quat4d cs{ fe1.x, fe2.x, fe1.y*dg, fe2.y*dg };
        //const Quat4d cs{ fe1.x, fe2.x, -fe1.y, -fe2.y };
        //const Quat4d cs{ fe1.x, fe2.x, 0, 0 };
        fes[i]=Vec2d{
            p.dot( cs ),
            d.dot( cs )*inv_dg,
            //fe1.x,
            //fe1.y,
        };
    }
}






__attribute__((hot)) 
void sample3D_deriv_comb3( const Vec3d g0, const Vec3d dg, const Vec3i ng, const Vec2d* FEg, const int n, const Vec3d* ps, Quat4d* fes, Vec3d C ){
    //printf( "sample3D_deriv_comb3() g0=(%g,%g,%g) dg=(%g,%g,%g) ng=(%i,%i,%i) n=%i C(%g,%g,%g)\n", g0.x,g0.y,g0.z, dg.x,dg.y,dg.z, ng.x,ng.y,ng.z, n, C.x,C.y,C.z );
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    Vec3d inv_dg2=inv_dg; inv_dg2.mul(-1);
    for(int i=0; i<n; i++ ){
        const Vec3d t  = (ps[i] - g0)*inv_dg; 
        const int ix = (int)t.x;
        const int iy = (int)t.y;
        const int iz = (int)t.z;
        // if( ((ix<1)||(ix>=ng.x-2)) || ((iy<0)||(iy>=ng.y-1)) || ((iz<0)||(iz>=ng.z-1)) )[[unlikely]]{ 
        //     //printf( "ERROR: Spline_Hermite::sample2D() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); 
        //     //printf( "Spline_Hermite::sample2D() ixyz(%i,%i,%i) out of range ng(%i,%i,%i) ... p[%i](%g,%g,%g)-> t(%g,%g,%g)\n", ix,iy,iz, ng.x,ng.y,ng.z, i, ps[i].x,ps[i].y,ps[i].z, t.x,t.y,t.z );
        //     //fes[i] = Quat4dZero;
        //     fes[i] = Quat4d{NAN,NAN,NAN,NAN};
        //     continue;
        // }

        Quat4d fe = fe3d_comb3( Vec3d{t.x-ix-1,t.y-iy-1,t.z-iz}, Vec3i{ix,iy,iz}, ng, FEg, C );
        fe.f.mul(inv_dg2);
        fes[i]=fe;

        // int ixy = ng.z*(iy + ix*ng.y)*3;
        // int ii = (iz + ng.z*(iy + ix*ng.y))*3;
        // alignas(32) const Quat4d pz  =  basis(t.z-iz);
        // alignas(32) const Quat4d dz  = dbasis(t.z-iz);
        // Vec2d efz = fe1Dcomb3( FEg+ii, C, pz, dz );
        // fes[i]=Quat4d{0.0,0.0,0.0, efz.x };

        //int ii = (iz + ng.z*(iy + ix*ng.y))*3;
        //printf( "iz=%i ii=%i \n", iz, ii );
        //fes[i]=Quat4d{0.0,0.0,  0.0, FEg[ii].x*C.x + FEg[ii+1].x*C.y + FEg[ii+2].x*C.z };
        //printf( "sample2D()[%i] ps(%g,%g) E=%g Fxy(%g,%g)\n", i, ps[i].x,ps[i].y,  fes[i].z,fes[i].x,fes[i].y );
    }
}


void sample1D_deriv_comb2( const double g0, const double dg, const int ng, const Quat4d* FE, const int n, const double* ps, Vec2d* fes, Vec2d C  ){
    const double inv_dg = 1/dg; 
    for(int i=0; i<n; i++ ){
        const double x = (ps[i] - g0)*inv_dg;  
        const int    ix = (int)x;
        const double tx =  x-ix; 
        const Quat4d p  = Spline_Hermite::basis(tx);
        const Quat4d d  = Spline_Hermite::dbasis(tx);
        Vec2d fe = Spline_Hermite::fe1Dcomb2( ix, FE, C, p, d );
        fe.y *= inv_dg;
        fes[i] = fe;
    }
}

__attribute__((hot)) 
void sample2D_deriv_comb( const Vec2d g0, const Vec2d dg, const Vec2i ng, const Quat4d* FEg, const int n, const Vec2d* ps, Vec3d* fes, Vec2d C ){
    printf( "sample2D_deriv_comb() g0=(%g,%g) dg=(%g,%g) ng=(%i,%i) n=%i C(%g,%g)\n", g0.x,g0.y, dg.x,dg.y, ng.x,ng.y, n, C.x,C.y );
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        const Vec2d t  = (ps[i] - g0)*inv_dg; 
        const int ix = (int)t.x;
        const int iy = (int)t.y;
        if( ((ix<1)||(ix>=ng.x-2)) || ((iy<0)||(iy>=ng.y-1)) )[[unlikely]]{ 
            //printf( "ERROR: Spline_Hermite::sample2D() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); 
            fes[i] = Vec3d{0.0,0.0,0.0};
            continue;
        }
        
        const double tx = t.x-ix;
        const double ty = t.y-iy;
        alignas(32) const Quat4d py  = Spline_Hermite::basis(ty);
        alignas(32) const Quat4d dy  = Spline_Hermite::dbasis(ty);
        alignas(32) const Quat4d bx  = Spline_Hermite::basis_val( tx );
        alignas(32) const Quat4d dx  = Spline_Hermite::dbasis_val( tx );
        Vec3d fe = Spline_Hermite::fe2d_comb2( iy, ng.y, FEg+(ix-1)*ng.y, C, py, dy, bx, dx );
        //Vec3d fe = Spline_Hermite::fe2d_comb2( iy-1, ng.y, FEg+(ix-1)*ng.y, C, py, dy, bx, dx );
        

        //ec3d fe{0.0,0.0, ix };
        //Vec3d fe{0.0,0.0, iy };
        //Vec3d fe{0.0,0.0, FEg[ ix*ng.y+iy ].x };
        //Vec3d fe{0.0,0.0, FEg[ iy*ng.x+ix ].x };

        //fe2d_deriv( t.x-ix,t.y-iy, {i0,i0+ng.x}, Eg, dEg, fe );        // sample2D(n=10000) time=527.478[kTick] 52.7478[tick/point]
        //Vec3d fe = fe2d_v2( t.x-ix,t.y-iy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );   // sample2D(n=10000) time=553.47[kTick] 55.347[tick/point]
        fe.x*=inv_dg.x;
        fe.y*=inv_dg.y;
        fes[i]=fe;
        //printf( "sample2D()[%i] ps(%g,%g) E=%g Fxy(%g,%g)\n", i, ps[i].x,ps[i].y,  fes[i].z,fes[i].x,fes[i].y );
    }
}

// =============== AVX2 version

__attribute__((hot)) 
inline void evalHermiteBasis_avx2( const Vec3d t, __m256d& bx, __m256d& by, __m256d& bz, __m256d& dx, __m256d& dy, __m256d& dz ){
    const __m256d mx = _mm256_set1_pd( t.x );
    const __m256d my = _mm256_set1_pd( t.y );
    const __m256d mz = _mm256_set1_pd( t.z );
    const __m256d c3 = _mm256_set_pd(  0.5, -1.5,  +1.5,  -0.5 );
    const __m256d c2 = _mm256_set_pd( -0.5,  2.0,  -2.5,   1.0 );
    const __m256d c1 = _mm256_set_pd(  0.0,  0.5,   0.0,  -0.5 );
    const __m256d c0 = _mm256_set_pd(  0.0,  0.0,   1.0,   0.0 );
    // poly = c1 + x*( c2 + x*( c3 + x*c4)))  -- Horner's scheme
    bx = _mm256_fmadd_pd( _mm256_fmadd_pd( _mm256_fmadd_pd( c3, mx, c2 ), mx,  c1 ), mx, c0  );
    by = _mm256_fmadd_pd( _mm256_fmadd_pd( _mm256_fmadd_pd( c3, my, c2 ), my,  c1 ), my, c0  );
    bz = _mm256_fmadd_pd( _mm256_fmadd_pd( _mm256_fmadd_pd( c3, mz, c2 ), mz,  c1 ), mz, c0  );
    const __m256d d2 = _mm256_set_pd(  1.5,  -4.5,    4.5,   -1.5  );
    const __m256d d1 = _mm256_set_pd( -1.0,   4.0,   -5.0,    2.0  );
    const __m256d d0 = _mm256_set_pd(  0.0,   0.5,    0.0,   -0.5  );
    dx = _mm256_fmadd_pd( _mm256_fmadd_pd( d2, mx, d1 ),  mx,  d0 );
    dy = _mm256_fmadd_pd( _mm256_fmadd_pd( d2, my, d1 ),  my,  d0 );
    dz = _mm256_fmadd_pd( _mm256_fmadd_pd( d2, mz, d1 ),  mz,  d0 );
}

__attribute__((hot)) 
inline void evalHermiteBasis_avx2( const int n, const double* ts, __m256d* bs, __m256d* ds, __m256d* dds=0 ){
    // -------- Value 
    //         x3     x2     x   
    //  p-1   -0.5   1.0  -0.5  -0.0
    //  p+0    1.5  -2.5   0.0   1.0  
    //  p+1   -1.5   2.0   0.5   0.0 
    //  p+2    0.5  -0.5   0.0   0.0  
    if(  bs )[[likely]]{
        // const __m256d c3 = _mm256_set_pd( -0.5, +1.5, -1.5,  0.5 );
        // const __m256d c2 = _mm256_set_pd(  1.0, -2.5,  2.0, -0.5 );
        // const __m256d c1 = _mm256_set_pd( -0.5,  0.0,  0.5,  0.0 );
        // const __m256d c4 = _mm256_set_pd(  0.0,  1.0,  0.0,  0.0 );
        const __m256d c3 = _mm256_set_pd(  0.5, -1.5,  +1.5,  -0.5 );
        const __m256d c2 = _mm256_set_pd( -0.5,  2.0,  -2.5,   1.0 );
        const __m256d c1 = _mm256_set_pd(  0.0,  0.5,   0.0,  -0.5 );
        const __m256d c0 = _mm256_set_pd(  0.0,  0.0,   1.0,   0.0 );
        // poly = c1 + x*( c2 + x*( c3 + x*c4)))  -- Horner's scheme
        for( int i=0; i<n; i++ ){
            const __m256d mx = _mm256_set1_pd( ts[i] );
            bs[i] = _mm256_fmadd_pd( _mm256_fmadd_pd( _mm256_fmadd_pd( c3, mx, c2 ), mx,  c1 ), mx, c0  );
        }
    }
    // -------- Derivative
    //         x2   x    1
    //  p-1  -1.5   2.0  -0.5   
    //  p+0   4.5  -5.0   0.0   
    //  p+1  -4.5   4.0   0.5   
    //  p+2   1.5  -1.0   0.0  
    if(  ds )[[likely]]{
        //const __m256d d2 = _mm256_set_pd( -1.5,  4.5, -4.5,  1.5 );
        //const __m256d d1 = _mm256_set_pd(  2.0, -5.0,  4.0, -1.0 );
        //const __m256d d0 = _mm256_set_pd( -0.5,  0.0,  0.5,  0.0 );
        const __m256d d2 = _mm256_set_pd(  1.5,  -4.5,    4.5,   -1.5  );
        const __m256d d1 = _mm256_set_pd( -1.0,   4.0,   -5.0,    2.0  );
        const __m256d d0 = _mm256_set_pd(  0.0,   0.5,    0.0,   -0.5  );
        for( int i=0; i<n; i++ ){
            const __m256d mx = _mm256_set1_pd( ts[i] );
            ds[i] = _mm256_fmadd_pd( _mm256_fmadd_pd( d2, mx, d1 ),  mx,  d0 );
            //ds[i] = _mm256_mul_pd( _mm256_mul_pd( d2, mx ),  mx );
        }
    }
    // -------- Acceleration
    //          x    1
    //  p-1   -3.0   2.  
    //  p+0   +9.0  -5.   
    //  p+1   -9.0   4.  
    //  p+2    3.0  -1.  
    if( dds )[[unlikely]]{
        //const __m256d d1 = _mm256_set_pd( -3.0,  9.0, -9.0,  3.0 );
        //const __m256d d0 = _mm256_set_pd(  2.0, -5.0,  4.0, -1.0 );
        const __m256d d1 = _mm256_set_pd(  3.0, -9.0,  9.0, -3.0 );
        const __m256d d0 = _mm256_set_pd( -1.0,  4.0, -5.0,  2.0 );
        for( int i=0; i<n; i++ ){
            const __m256d mx = _mm256_set1_pd( ts[i] );
            dds[i] = _mm256_fmadd_pd( d1, mx, d0 );
        }
    }
}

__attribute__((hot)) 
Vec3d fe2d_avx( const __m256d mbx, const __m256d mdx, const Quat4d by, const Quat4d dy, const Quat4i i, const double* Eg ){
    // --- read data from global memory ( 4 doubles at a time along x-axis )
    // --- aligned load
    //const __m256d e1 = _mm256_load_pd( Eg+i.x );
    //const __m256d e2 = _mm256_load_pd( Eg+i.y );
    //const __m256d e3 = _mm256_load_pd( Eg+i.z );
    //const __m256d e4 = _mm256_load_pd( Eg+i.w );
    // --- un-aligned load
    const __m256d e1 = _mm256_loadu_pd( Eg+i.x );
    const __m256d e2 = _mm256_loadu_pd( Eg+i.y );
    const __m256d e3 = _mm256_loadu_pd( Eg+i.z );
    const __m256d e4 = _mm256_loadu_pd( Eg+i.w );
    // --- interpolate along y-axis ( for x-lines in parallel )
    const __m256d mE = 
        _mm256_fmadd_pd(e4, _mm256_set1_pd(by.w), 
        _mm256_fmadd_pd(e3, _mm256_set1_pd(by.z), 
        _mm256_fmadd_pd(e2, _mm256_set1_pd(by.y), 
        _mm256_mul_pd  (e1, _mm256_set1_pd(by.x)))));
    const __m256d mFy =
        _mm256_fmadd_pd(e4, _mm256_set1_pd(dy.w), 
        _mm256_fmadd_pd(e3, _mm256_set1_pd(dy.z), 
        _mm256_fmadd_pd(e2, _mm256_set1_pd(dy.y), 
        _mm256_mul_pd  (e1, _mm256_set1_pd(dy.x)))));
    alignas(32) Quat4d E ; _mm256_store_pd( (double*)&E , _mm256_mul_pd( mbx, mE  ) );
    alignas(32) Quat4d Fx; _mm256_store_pd( (double*)&Fx, _mm256_mul_pd( mdx, mE  ) );
    alignas(32) Quat4d Fy; _mm256_store_pd( (double*)&Fy, _mm256_mul_pd( mbx, mFy ) );
    return Vec3d{
        Fx.x+Fx.y+Fx.z+Fx.w, // Fx
        Fy.x+Fy.y+Fy.z+Fy.w, // Fy
        E .x+E .y+E .z+E .w  // E
    };
}

__attribute__((hot)) 
Quat4d fe3d_avx( const Vec3d u, const Vec3i n, const double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int   ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const Vec3d dt{ u.x - ix  ,  u.y - iy  ,  u.z - iz  };
    // if( 
    //     ((ix<0)||(ix>=n.x-3)) ||
    //     ((iy<0)||(iy>=n.y-3)) ||
    //     ((iz<0)||(iz>=n.z-3))        
    // )[[unlikely]]{ printf( "ERROR: Spline_Hermite::interpolateTricubic() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }
    __m256d bx,by,bz,dx,dy,dz;
    evalHermiteBasis_avx2( dt, bx,by,bz,dx,dy,dz );
    alignas(32) Quat4d qby; _mm256_store_pd( (double*)&qby, by );
    alignas(32) Quat4d qdy; _mm256_store_pd( (double*)&qdy, dy );
    // __m256d bs[3];
    // __m256d ds[3];
    // evalHermiteBasis_avx2( 3, (double*)&dt, bs, ds  );
    //alignas(32) Quat4d qby; _mm256_store_pd( (double*)&qby, bs[1] );
    //alignas(32) Quat4d qdy; _mm256_store_pd( (double*)&qdy, ds[1] );
    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    int i0 = ix + n.x*( iy + n.y*iz ); 
    const Vec3d Exy1 = fe2d_avx( bx,dx, qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy2 = fe2d_avx( bx,dx, qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy3 = fe2d_avx( bx,dx, qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy4 = fe2d_avx( bx,dx, qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
    alignas(32) Quat4d qbz; _mm256_store_pd( (double*)&qbz, bz );
    alignas(32) Quat4d qdz; _mm256_store_pd( (double*)&qdz, dz );
    return Quat4d{
        qbz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        qbz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        qdz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        qbz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
    };
} 

__attribute__((hot)) 
void sample2D_avx( const Vec2d g0, const Vec2d dg, const Vec2i ng, const double* Eg, const int n, const Vec2d* ps, Vec3d* fes ){
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    bool err = false;
    for(int i=0; i<n; i++ ){
        const Vec2d t  = (ps[i] - g0)*inv_dg; 
        const int ix = (int)t.x;
        const int iy = (int)t.y;
        if( ((ix<0)||(ix>=ng.x-3)) || ((iy<0)||(iy>=ng.y-3)) )[[unlikely]]{ printf( "ERROR: Spline_Hermite::interpolateTricubic() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); }
        const int i0 = ix + ng.x*iy;
        const Vec3d dt{ t.x-ix, t.y-iy, 0.0 };
        // __m256d bs[3];
        // __m256d ds[3];
        // //__m256d as[3];
        // evalHermiteBasis_avx2( 2, (double*)&dt, bs, ds );
        // alignas(32) Quat4d qby; _mm256_store_pd( (double*)&qby, bs[1] );
        // alignas(32) Quat4d qdy; _mm256_store_pd( (double*)&qdy, ds[1] );

        __m256d bx,by,bz,dx,dy,dz;
        evalHermiteBasis_avx2( dt, bx,by,bz,dx,dy,dz );
        alignas(32) Quat4d qby; _mm256_store_pd( (double*)&qby, by );
        alignas(32) Quat4d qdy; _mm256_store_pd( (double*)&qdy, dy );

        // {  // Check basis
        //     alignas(32) Quat4d qbx; _mm256_store_pd( (double*)&qbx, bs[0] );
        //     alignas(32) Quat4d qdx; _mm256_store_pd( (double*)&qdx, ds[0] );

        //     // alignas(32) Quat4d qax; _mm256_store_pd( (double*)&qax, as[0] );
        //     // const Quat4d ax =  ddbasis_val( dt.x );
        //     // if( (qax-ax).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qax(%g,%g,%g,%g) != ax(%g,%g,%g,%g)\n",i, qax.x,qax.y,qax.z,qax.w, ax.x,ax.y,ax.z,ax.w ); err=true; }

        //     const Quat4d bx =  basis_val( dt.x );
        //     const Quat4d dx = dbasis_val( dt.x );
        //     const Quat4d by =  basis_val( dt.y );
        //     const Quat4d dy = dbasis_val( dt.y );

        //     const Quat4d bx_ =  basis_val_2( dt.x );
        //     const Quat4d by_ =  basis_val_2( dt.y );
        //     const Quat4d dx_ = dbasis_val_2( dt.x );
        //     const Quat4d dy_ = dbasis_val_2( dt.y );

        //     if( (bx-bx_).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] bx_(%g,%g,%g,%g) != bx(%g,%g,%g,%g)\n",i, bx_.x,bx_.y,bx_.z,bx_.w, bx.x,bx.y,bx.z,bx.w ); err=true; }
        //     if( (by-by_).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] by_(%g,%g,%g,%g) != by(%g,%g,%g,%g)\n",i, by_.x,by_.y,by_.z,by_.w, by.x,by.y,by.z,by.w ); err=true; }
        //     if( (dx-dx_).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] dx_(%g,%g,%g,%g) != dx(%g,%g,%g,%g)\n",i, dx_.x,dx_.y,dx_.z,dx_.w, dx.x,dx.y,dx.z,dx.w ); err=true; }
        //     if( (dy-dy_).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] dy_(%g,%g,%g,%g) != dy(%g,%g,%g,%g)\n",i, dy_.x,dy_.y,dy_.z,dy_.w, dy.x,dy.y,dy.z,dy.w ); err=true; }

        //     if( (bx-qbx).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qbx(%g,%g,%g,%g) != bx(%g,%g,%g,%g)\n",i, qbx.x,qbx.y,qbx.z,qbx.w, bx.x,bx.y,bx.z,bx.w ); err=true; }
        //     if( (dx-qdx).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qdx(%g,%g,%g,%g) != dx(%g,%g,%g,%g)\n",i, qdx.x,qdx.y,qdx.z,qdx.w, dx.x,dx.y,dx.z,dx.w ); err=true; }
        //     if( (by-qby).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qby(%g,%g,%g,%g) != by(%g,%g,%g,%g)\n",i, qby.x,qby.y,qby.z,qby.w, by.x,by.y,by.z,by.w ); err=true; }
        //     if( (dy-qdy).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qdy(%g,%g,%g,%g) != dy(%g,%g,%g,%g)\n",i, qdy.x,qdy.y,qdy.z,qdy.w, dy.x,dy.y,dy.z,dy.w ); err=true; } 

        //     if( (bx_-qbx).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qbx(%g,%g,%g,%g) != bx_(%g,%g,%g,%g)\n",i, qbx.x,qbx.y,qbx.z,qbx.w, bx_.x,bx_.y,bx_.z,bx_.w ); err=true; }
        //     if( (dx_-qdx).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qdx(%g,%g,%g,%g) != dx_(%g,%g,%g,%g)\n",i, qdx.x,qdx.y,qdx.z,qdx.w, dx_.x,dx_.y,dx_.z,dx_.w ); err=true; }
        //     if( (by_-qby).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qby(%g,%g,%g,%g) != by_(%g,%g,%g,%g)\n",i, qby.x,qby.y,qby.z,qby.w, by_.x,by_.y,by_.z,by_.w ); err=true; }
        //     if( (dy_-qdy).norm2() > 1e-14 ){ printf( "ERROR: sample2D_avx()[%i] qdy(%g,%g,%g,%g) != dy_(%g,%g,%g,%g)\n",i, qdy.x,qdy.y,qdy.z,qdy.w, dy_.x,dy_.y,dy_.z,dy_.w ); err=true; } 

        //     if(err){ printf( "ERROR: sample2D_avx()[%i] p(%g,%g) t(%g,%g) \n", i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); }
        // }

        //Vec3d fe = fe2d_avx( bs[0],ds[0], qby,qdy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );
        Vec3d fe = fe2d_avx( bx,dx, qby,qdy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );
        fe.x*=inv_dg.x;
        fe.y*=inv_dg.x;
        fes[i]=fe;
        
        //printf( "sample2D()[%i] ps(%g,%g) E=%g Fxy(%g,%g)\n", i, ps[i].x,ps[i].y,  fes[i].z,fes[i].x,fes[i].y );
    }
}

__attribute__((hot)) 
void sample3D_avx( const Vec3d g0, const Vec3d dg, const Vec3i ng, const double* Eg, const int n, const Vec3d* ps, Quat4d* fes ){
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        Quat4d fe = fe3d_avx( (ps[i]-g0)*inv_dg, ng, Eg );
        fe.f.mul(inv_dg);
        fes[i] = fe;
    }
}



}; // namespace Spline_Hermite

#endif
