#ifndef InterpolateTricubic_h
#define InterpolateTricubic_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "spline_hermite.h"

#include <immintrin.h>

namespace Spline_Hermite{

template<typename T>
inline Quat4T<T> basis_val_2( T x  ){
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
inline Quat4T<T> dbasis_val_2( T x ){
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
inline Quat4T<T> ddbasis_val_2( T x ){
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


template<typename T>
inline Quat4T<T> basis_val( T x  ){
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
inline Quat4T<T> dbasis_val( T x ){
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
inline Quat4T<T> ddbasis_val( T x ){
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

__attribute__ ((pure))
__attribute__((hot)) 
Vec3d fe2d( const double tx, const double ty, const Quat4i i, const double* Es ){
    Quat4d e,fx;
    {
        const Quat4d bx =  basis_val( tx );
        const Quat4d dx = dbasis_val( tx );
        {
            const Quat4d p = *(Quat4d*)(Es+i.x); // read 4 doubles from global memory at a time ( 4*8 = 32 bytes = 256 bits ) ideal for SIMD AVX2
            e.x  = bx.dot(p);   // not sure how dot() is SIMD optimized => maybe we should flip the order of x and y strides ?
            fx.x = dx.dot(p);
            //fy.y = bx.dot(p);
        }
        {
            const Quat4d p = *(Quat4d*)(Es+i.y); 
            e.y  = bx.dot(p);
            fx.y = dx.dot(p);
            //fy.y = bx.dot(p);
        }
        {
            const Quat4d p = *(Quat4d*)(Es+i.z); 
            e.z  = bx.dot(p);
            fx.z = dx.dot(p);
            //fy.z = bx.dot(p);
        }
        {
            const Quat4d p = *(Quat4d*)(Es+i.w); 
            e.w  = bx.dot(p);
            fx.w = dx.dot(p);
            //fy.w = bx.dot(p);
        }
    }
    const Quat4d by =  basis_val( ty );
    const Quat4d dy = dbasis_val( ty );
    return Vec3d{
        by.dot(fx), // Fx
        dy.dot(e ), // Fy
        by.dot(e )  // E
    };
}

__attribute__ ((pure))
__attribute__((hot)) 
Quat4d fe3d( const Vec3d u, const Vec3i n, const double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    const double mx = 1-tx      ,  my = 1-ty      ,  mz = 1-tz      ;

    if( 
        ((ix<0)||(ix>=n.x-3)) ||
        ((iy<0)||(iy>=n.y-3)) ||
        ((iz<0)||(iz>=n.z-3))        
    )[[unlikely]]{ printf( "ERROR: Spline_Hermite::interpolateTricubic() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }

    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    int i0 = ix + n.x*( iy + n.y*iz ); 
    const Vec3d Exy1 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy2 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy3 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy4 = fe2d(tx,ty, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
    const Quat4d bz =  basis_val( tz );
    const Quat4d dz = dbasis_val( tz );
    return Quat4d{
        bz.dot( {Exy1.x, Exy2.x, Exy3.x, Exy4.x} ), // Fx
        bz.dot( {Exy1.y, Exy2.y, Exy3.y, Exy4.y} ), // Fy
        dz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // Fz
        bz.dot( {Exy1.z, Exy2.z, Exy3.z, Exy4.z} ), // E
    };
} 


//#endif WITH_AVX

__attribute__((hot)) 
void sample1D( const double g0, const double dg, const int ng, const double* Eg, const int n, const double* xs, double* Es, double* Fs ){
    const double inv_dg = 1/dg;
    for(int i=0; i<n; i++ ){
        const double t = ((xs[i]-g0)*inv_dg);
        const int it   = (int)t;
        //printf( "sample_SplineHermite() x[%i]=%g t=%g it=%i | g0=%g dg=%g \n", i, xs[i], t, it, g0, dg );
        if( (it<0)||(it>=ng-3) )[[unlikely]]{ printf( "ERROR: sample_SplineHermite it(%i) out of range (0,%i) | xs[%i]=%g -> t=%g \n", it, ng, i, xs[i], t ); exit(0); }
        double dt = t-it;
        const Quat4d bs = basis_val ( dt );
        const Quat4d ds = dbasis_val( dt );
        const Quat4d v = *(Quat4d*)(Eg+it);
        Es[i] = bs.dot(v);
        Fs[i] = ds.dot(v)*inv_dg;
    }
}

__attribute__((hot)) 
void sample2D( const Vec2d g0, const Vec2d dg, const Vec2i ng, const double* Eg, const int n, const Vec2d* ps, Vec3d* fes ){
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        const Vec2d t  = (ps[i] - g0)*inv_dg; 
        const int ix = (int)t.x;
        const int iy = (int)t.y;
        if( ((ix<0)||(ix>=ng.x-3)) || ((iy<0)||(iy>=ng.y-3)) )[[unlikely]]{ printf( "ERROR: Spline_Hermite::interpolateTricubic() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); }
        const int i0 = ix + ng.x*iy;
        Vec3d fe = fe2d( t.x-ix,t.y-iy, {i0,i0+ng.x,i0+ng.x*2,i0+ng.x*3}, Eg );
        fe.x*=inv_dg.x;
        fe.y*=inv_dg.x;
        fes[i]=fe;
        //printf( "sample2D()[%i] ps(%g,%g) E=%g Fxy(%g,%g)\n", i, ps[i].x,ps[i].y,  fes[i].z,fes[i].x,fes[i].y );
    }
}

__attribute__((hot)) 
void sample3D( const Vec3d g0, const Vec3d dg, const Vec3i ng, const double* Eg, const int n, const Vec3d* ps, Quat4d* fes ){
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        fes[i] = fe3d( (ps[i]-g0)*inv_dg, ng, Eg );
    }
}

// =============== AVX2 version

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
    if( dds ){
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
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const Vec3d t{ u.x - ix  ,  u.y - iy  ,  u.z - iz  };

    if( 
        ((ix<0)||(ix>=n.x-3)) ||
        ((iy<0)||(iy>=n.y-3)) ||
        ((iz<0)||(iz>=n.z-3))        
    )[[unlikely]]{ printf( "ERROR: Spline_Hermite::interpolateTricubic() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) t(%g,%g,%g)\n", ix,iz,iy, n.x,n.y,n.z, u.x,u.y,u.z ); exit(0); }
    //__m256d bx,by,bz,dx,dy,dz;
    __m256d bs[3];
    __m256d ds[3];

    evalHermiteBasis_avx2( 3, (double*)&t, bs, ds  );
    alignas(32) Quat4d qby; _mm256_store_pd( (double*)&qby, bs[1] );
    alignas(32) Quat4d qdy; _mm256_store_pd( (double*)&qdy, bs[1] );
    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;
    int i0 = ix + n.x*( iy + n.y*iz ); 
    const Vec3d Exy1 = fe2d_avx( bs[0],ds[0], qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy2 = fe2d_avx( bs[0],ds[0], qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy3 = fe2d_avx( bs[0],ds[0], qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );  i0 += nxy;
    const Vec3d Exy4 = fe2d_avx( bs[0],ds[0], qby,qdy, {i0,i0+n.x,i0+n.x*2,i0+3*n.x}, Es );
    alignas(32) Quat4d qbz; _mm256_store_pd( (double*)&qbz, bs[2] );
    alignas(32) Quat4d qdz; _mm256_store_pd( (double*)&qdz, ds[2] );
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



}; // namespace Spline_Hermite

#endif
