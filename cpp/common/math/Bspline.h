
#ifndef  Bspline_h
#define  Bspline_h

#include "quaternion.h"

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

//#include <math.h>
//#include <cstdlib>
//#include <stdio.h>

// ============ optimized

namespace Bspline{

constexpr static const double CPs[]{ 0.0,  1.0/6.0,   2.0/3.0,  1.0/6.0 };
constexpr static const double DPs[]{ 0.0,  0.5,       0.0,     -0.5     };
constexpr static const double APs[]{ 0.0,  1.0,      -2.0,      1.0     };

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

__attribute__((hot)) 
void sample1D( const double g0, const double dg, const int ng, const double* Gs, const int n, const double* ps, Vec2d* fes ){
    const double inv_dg = 1/dg; 
    for(int i=0; i<n; i++ ){
        const double x = (ps[i] - g0)*inv_dg;  
        const int    ix = (int)x;
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
int fit1D( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1 ){
    const double F2max = Ftol*Ftol;
    double* ps = new double[n];
    double* fs = new double[n];
    double* vs = new double[n];
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    int itr=0;
    for(itr=0; itr<nmaxiter; itr++){
        // --- evaluate current spline
        for(int i=0; i<n; i++){
            double p = Gs[i]*B0;
            if(i>0  )[[likely]]{ p += B1*Gs[i-1]; }
            if(i<n-1)[[likely]]{ p += B1*Gs[i+1]; }
            ps[i] = p;   Ws[i] = p;
            fs[i] = 0;
        }
        // --- evaluate variatiaonal derivatives
        for(int i=0; i<n; i++){
            const double dp  = Es[i] - ps[i];
            fs[i] += dp*B0;
            if(i>0  )[[likely]]{ fs[i-1] += B1*dp; }
            if(i<n-1)[[likely]]{ fs[i+1] += B1*dp; }
        }
        // --- move
        double vf = 0.0;
        double ff = 0.0;
        for(int i=0; i<n; i++){
            ff += fs[i]*fs[i];
            vf += vs[i]*fs[i];
        }
        printf( "|F[%i]|=%g \n",itr,  sqrt(ff) );
        //printf( "p=%i v=%g f=%g \n",  Gs[1], vs[1], fs[1] );
        if(ff<F2max){ break; }
        if(vf<0){ for(int i=0; i<n; i++){ vs[i]=0; }; }
        for(int i=0; i<n; i++){
            vs[i] += fs[i]*dt;
            Gs[i] += vs[i]*dt;
            //Gs[i] += fs[i]*dt;  // GD
        }
        
    }
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
int fit1D_EF( const double dg, const int n, double* Gs,  Vec2d* fes, Vec2d* Ws, double Ftol, int nmaxiter=1000, double dt=0.1 ){
    const double inv_dg = 1/dg;
    const double F2max = Ftol*Ftol;
    Vec2d*  pfs = new Vec2d[n];
    double* fs  = new double[n];
    double* vs  = new double[n];
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double D1=0.5;
    int itr=0;

    for(int i=0; i<n; i++){  Ws[i]= Vec2dZero; pfs[i] = Vec2dZero;   };

    for(itr=0; itr<nmaxiter; itr++){
        // --- evaluate current spline
        for(int i=0; i<n-1; i++){
            Vec2d p = Vec2d{ Gs[i]*B0, 0.0};
            // if(i>0  )[[likely]]{ p.add_mul( Vec2d{B1,-D1*inv_dg}, Gs[i-1] ); }
            // if(i<n-1)[[likely]]{ p.add_mul( Vec2d{B1,+D1*inv_dg}, Gs[i+1] ); }
            p.add_mul( Vec2d{B1,-D1*inv_dg}, Gs[i-1] );
            p.add_mul( Vec2d{B1,+D1*inv_dg}, Gs[i+1] ); 
            pfs[i] = p;    //Ws[i] = p;
            fs [i] = 0;
            Ws[i].x = p.y;
        }
        // --- evaluate variatiaonal derivatives
        for(int i=2; i<n-2; i++){
            Vec2d fei = fes[i];
            //fei.y*=-1;
            const Vec2d dp  = fei - pfs[i];
            fs[i] += dp.x*B0*0;
            //if(i>0  )[[likely]]{ fs[i-1] += B1*dp.x*0 + -D1*dp.y; }
            //if(i<n-1)[[likely]]{ fs[i+1] += B1*dp.x*0 + +D1*dp.y; }
            fs[i-1] += B1*dp.x*0 + -D1*dp.y;
            fs[i+1] += B1*dp.x*0 + +D1*dp.y;
        }
        for(int i=0; i<n; i++){  Ws[i].y = fs[i];  } // Debug

        // --- move
        double vf = 0.0;
        double ff = 0.0;
        for(int i=0; i<n; i++){
            ff += fs[i]*fs[i];
            vf += vs[i]*fs[i];
        }
        printf( "|F[%i]|=%g \n",itr,  ff );
        //printf( "p=%i v=%g f=%g \n",  Gs[1], vs[1], fs[1] );
        if(ff<F2max){ break; }
        //if(vf<0){ for(int i=0; i<n; i++){ vs[i]=Vec2dZero; }; }
        for(int i=2; i<n-2; i++){
            //vs[i] += fs[i]*dt;
            //Gs[i] += vs[i]*dt;

            Gs[i] += fs[i]*dt;  // GD
        }
        
    }
    delete [] pfs;
    delete [] fs;
    delete [] vs;
    return itr;
}

} // namespace SplineBcub{

#endif



