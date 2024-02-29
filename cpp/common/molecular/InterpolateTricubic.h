#ifndef InterpolateTricubic_h
#define InterpolateTricubic_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "spline_hermite.h"

namespace Spline_Hermite{

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
	 -2*K   + d0* 0.5,  //  p+1 = c1 + +0.5*d0
	          d1* 0.5   //  p+2 =      +0.5*d1
    };
}

template <class T>
inline Quat4T<T> ddbasis_val( T x ){
	const T x6  =  6*x;
    const T d0  =  x6 -  4;        //     6*x - 4
	const T d1  =  x6 -  2;        //     6*x - 2
    return Quat4T<T>{
	                  d0*-0.5, //  p-1 =      -0.5*d0
       x6 + x6 -  6 + d1*-0.5, //  p+0 = c0 + -0.5*d1
	    6 - x6 - x6 + d0* 0.5,  //  p+1 = c1 + +0.5*d0
	                  d1* 0.5   //  p+2 =      +0.5*d1
    };
}

__attribute__((hot)) 
Vec3d fe2d( double tx, double ty, Quat4i i, double* Es ){
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
            e.x  = bx.dot(p);
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

__attribute__((hot)) 
Quat4d fe3d( Vec3d u, Vec3i n, double* Es ){
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

//__attribute__((hot)) 
void sample1D( double g0, double dg, int ng, double* Eg, int n, double* xs, double* Es, double* Fs ){
    const double inv_dg = 1/dg;
    for(int i=0; i<n; i++ ){
        const double t = (xs[i]-g0)*inv_dg;
        const int it   = (int)t;
        if( (it<0)||(it>=ng-3) )[[unlikely]]{ printf( "ERROR: sample_SplineHermite it(%i) out of range (0,%i) | xs[%i]=%g -> t=%g \n", it, ng, i, xs[i], t ); exit(0); }
        double dt = t-it;
        const Quat4d bs = basis_val ( dt );
        const Quat4d ds = dbasis_val( dt );
        const Quat4d v = *(Quat4d*)(Eg+it);
        Es[i] = bs.dot(v);
        Fs[i] = ds.dot(v);
    }
}

//__attribute__((hot)) 
void sample2D( Vec2d g0, Vec2d dg, Vec2i ng, double* Eg, int n, Vec2d* ps, Vec3d* fes ){
    Vec2d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        const Vec2d t  = (ps[i] - g0)*inv_dg; 
        const int ix = (int)t.x;
        const int iy = (int)t.y;
        if( ((ix<0)||(ix>=ng.x-3)) || ((iy<0)||(iy>=ng.y-3)) )[[unlikely]]{ printf( "ERROR: Spline_Hermite::interpolateTricubic() ixyz(%i,%i) out of range 0 .. (%i,%i) p[%i](%g,%g)-> t(%g,%g)\n", ix,iy, ng.x,ng.y, i, ps[i].x,ps[i].y, t.x,t.y ); exit(0); }
        const int i0 = ix + ng.x*iy;
        fes[i] = fe2d( t.x,t.y, {i0,i0+ng.x,i0+ng.x*2,i0+3*ng.x}, Eg );
    }
}

//__attribute__((hot)) 
void sample3D( Vec3d g0, Vec3d dg, Vec3i ng, double* Eg, int n, Vec3d* ps, Quat4d* fes ){
    Vec3d inv_dg; inv_dg.set_inv(dg); 
    for(int i=0; i<n; i++ ){
        fes[i] = fe3d( (ps[i]-g0)*inv_dg, ng, Eg );
    }
}

}; // namespace Spline_Hermite

#endif
