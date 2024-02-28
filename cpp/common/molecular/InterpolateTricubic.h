#ifndef InterpolateTricubic_h
#define InterpolateTricubic_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "spline_hermite.h"


template<typename T>
inline Quat4T<T>& basis_val( T x  ){
	const T x2 = x*x;
	const T K  =  x2*(x - 1);
	const T d0 =    K - x2 + x;       //      x3 - 2*x2 + x
	const T d1 =    K         ;       //      x3 -   x2   
    return Quat4T<T>{
	                d0*-0.5, //  p-1 =      -0.5*d0
     2*K - x2 + 1 + d1*-0.5, //  p+0 = c0 + -0.5*d1
	-2*K + x2     + d0* 0.5,  //  p+1 = c1 + +0.5*d0
                    d1* 0.5   //  p+2 =      +0.5*d1
    }
}

template <class T>
inline Quat4T<T>& dbasis_val( T x ){
	const T K    =  3*x*(x - 1);
    const T d0   =    K - x + 1;   //    3*x2 - 4*x + 1
	const T d1   =    K + x    ;   //    3*x2 - 2*x
    return Quat4T<T>{
              d0*-0.5, //  p-1 =      -0.5*d0
	  2*K   + d1*-0.5, //  p+0 = c0 + -0.5*d1
	 -2*K   + d0* 0.5,  //  p+1 = c1 + +0.5*d0
	          d1* 0.5   //  p+2 =      +0.5*d1
    }
}

template <class T>
inline Quat4T<T>& ddbasis_val( T x ){
	const T x6  =  6*x;
    const T d0  =  x6 -  4;        //     6*x - 4
	const T d1  =  x6 -  2;        //     6*x - 2
    return Quat4T<T>{
	                  d0*-0.5, //  p-1 =      -0.5*d0
       x6 + x6 -  6 + d1*-0.5, //  p+0 = c0 + -0.5*d1
	    6 - x6 - x6 + d0* 0.5,  //  p+1 = c1 + +0.5*d0
	                  d1* 0.5   //  p+2 =      +0.5*d1
    }
}

__attribute__((hot)) 
Vec3d interpolateTricubic( double tx, double ty, Quat4i i, double* Es ){
    Quat4d e,fx;
    {
        const Quat4d bx =  basis_val( tx );
        const Quat4d dx = dbasis_val( tx );
        {
            const Quat4d p = *(Quat4d*)(Es+i.x); 
            e.x  = bx.dot(p);
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
    //E.x  = by.dot(e );
    //Fy.x = dy.dot(e );
    //Fx.x = by.dot(fx);
    return Vec3d{
        by.dot(e ), // E
        dy.dot(e ), // Fy
        by.dot(fx)  // Fx
    };
}



__attribute__((hot)) 
Quat4d interpolateTricubic( Vec3d u, Vec3i n, double* Es ){
    // We assume there are boundary added to simplify the index calculations
	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    const double mx = 1-tx      ,  my = 1-ty      ,  mz = 1-tz      ;

    //Quat4d E,Fx,Fy;

    int i0 = ix + n.x*( iy + n.y*iz );

    Vec3d Exy1 = interpolateTricubic(tx,ty, {i0,i0+n.y,i0+n.y*2,i0+3*n.y}, Es );
    Vec3d Exy2 = interpolateTricubic(tx,ty, {i0,i0+n.y,i0+n.y*2,i0+3*n.y}, Es );
    Vec3d Exy3 = interpolateTricubic(tx,ty, {i0,i0+n.y,i0+n.y*2,i0+3*n.y}, Es );
    Vec3d Exy4 = interpolateTricubic(tx,ty, {i0,i0+n.y,i0+n.y*2,i0+3*n.y}, Es );

    const Quat4d bz =  basis_val( tz );
    const Quat4d dz = dbasis_val( tz );

    // e.w  = 
    // fx.w = dz.dot(p);

    // return {
    //     bz.dot({ Exy1.x  });

    // }

    return Quat4dZero;

} 

#endif
