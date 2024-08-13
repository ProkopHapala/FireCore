
#ifndef  NURBS_h
#define  NURBS_h

#include "quaternion.h"
#include "CG.h"
#include "globals.h"
#include "Bspline.h"

/*

NURBS : Non-uninorm Besier spline

C(t) = sum_i{  w_i B_i(t) y_i  } / sum_i{  w_i B_i(t) }


*/

namespace NURBS{

// ======================================================  
// ===================   Interpolation   ================
// ======================================================

inline Vec2d fe1d( double t, const Quat4d& ys, const Quat4d& ws ){
    const Quat4d p  =  Bspline::basis (t)*ws;
    const Quat4d d  =  Bspline::dbasis(t)*ws;
    double db    = d.dot( ys );
    double b     = p.dot( ys );
    double f     = p.x+p.y+p.z+p.w;
    double F     = 1/f;
    double df    = d.x+d.y+d.z+d.w; 
    return Vec2d{ b*F, (db*f - b*df)*F*F };
}

__attribute__((hot)) 
void sample1D( const double g0, const double dg, const int ng, const double* Gs, const double* Ws, const int n, const double* ps, Vec2d* fes ){
    //printf("NURBS::sample1D() \n");
    const double inv_dg = 1/dg; 
    for(int i=0; i<n; i++ ){
        const double x  = (ps[i] - g0)*inv_dg;  
        const int    ix = (int)x;
        if( ((ix<1)||(ix>=(ng-2))) )[[unlikely]]{ 
            continue;
            //printf( "ERROR: Bspline::sample1D() ixyz(%i,%i) out of range 0 .. (%i,%i) u(%g,%g)\n", ix,iy, n.x,n.y, u.x,u.y );   exit(0); 
        }
        const double tx = x-ix; 
        Vec2d fe = fe1d( tx, *(Quat4d*)(Gs+ix-1), *(Quat4d*)(Ws+ix-1) );
        fe.y *= inv_dg;
        fes[i]=fe;
    }
}

} // namespace SplineBcub{

#endif



