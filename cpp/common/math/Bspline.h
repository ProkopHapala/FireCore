
#ifndef  Bspline_h
#define  Bspline_h

#include "quaternion.h"
#include "CG.h"



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

__attribute__((pure))
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
    )[[unlikely]]{ 
        printf( "ERROR: Spline_Hermite::fe3d() ixyz(%i,%i,%i) out of range 0 .. (%i,%i,%i) u(%g,%g,%g)\n", ix,iy,iz, n.x,n.y,n.z, u.x,u.y,u.z ); 
        //printf( "DETAILS:",   u.x ); 
        exit(0); 
    }
    //Quat4d E,Fx,Fy;
    const int nxy = n.x*n.y;

    //printf( "ixyz(%i,%i,%i)  (%g,%g,%g)\n", ix,iy,iz,  u.x,u.y,u.z );

    int i0 = ix + n.x*( iy + n.y*iz );  
    
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


Vec3d move( double dt, int n, double* gs, double* fs, double* vs ){
        constexpr const int m = 2;
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
int fit1D( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1 ){
    printf("Bspline::fit1D() \n");
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
    for(itr=0; itr<nmaxiter; itr++){
        // --- evaluate current spline
        for(int i=0; i<n; i++){
            double p = Gs[i]*B0;
            if(i>0  )[[likely]]{ p += B1*Gs[i-1]; }
            if(i<n-1)[[likely]]{ p += B1*Gs[i+1]; }
            //bErr|=checkNumRange( i, T val, T min, T max, const char* pre, bool bPrint=true, bool bExit=false ){
            //printf( "p[%i]=%g Gs=%g\n", i, p, Gs[i] );
            //bErr|=  checkNumRange( i, p, -1.e+6, +1.e+6, "p=B*Gs", true, true );
            ps[i] = p;   //Ws[i] = p;
            fs[i] = 0;
        }
        // --- evaluate variatiaonal derivatives
        for(int i=0; i<n; i++){
            double dp  = Es[i] - ps[i];
            if( Ws ){ dp *= Ws[i]; } // weighting 
            fs[i] += dp*B0;
            //bErr|=  checkNumRange( i, dp, -1.e+6, +1.e+6, "dp=Es-ps", true, true );
            if(i>0  )[[likely]]{ fs[i-1] += B1*dp; }
            if(i<n-1)[[likely]]{ fs[i+1] += B1*dp; }
        }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        printf( "|F[%i]|=%g \n",itr,sqrt(cfv.y) );
        if(cfv.y<F2max){ break; };
    }
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
int fit1D_EF( const double dg, const int n, double* Gs,  Vec2d* fes, Vec2d* Ws, double Ftol, int nmaxiter=1000, double dt=0.1 ){
    printf("Bspline::fit1D_EF() \n");
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
    for(itr=0; itr<nmaxiter; itr++){
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
        printf( "|F[%i]|=%g \n",itr,sqrt(cfv.y) );
        if(cfv.y<F2max){ break; };
    }
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}


void getVariations3D( const Vec3i ns, double* Gs,  double* Es, double* Ws, double* fs, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    // constexpr double B00=B0*B0;
    // constexpr double B01=B0*B1;
    // constexpr double B11=B1*B1;
    constexpr double B000=B0*B0*B0;
    constexpr double B001=B0*B0*B1;
    constexpr double B011=B0*B1*B1;
    constexpr double B111=B1*B1*B1;
    const int nxy  = ns.x*ns.y;
    const int nxyz = ns.x*ns.y*ns.z;

    // --- evaluate current spline
    for(int i=0; i<nxyz; i++){ fs[i]=0; ps[i]=0;  }

    // --- evaluate current spline (in order to evelauet approximation error)
    for(int iz=1; iz<ns.z-1; iz++){
        int iiz = iz*nxy;
        for(int iy=1; iy<ns.y-1; iy++){
            int iiy = iy*ns.x;
            for(int ix=1; ix<ns.x-1; ix++){
                double val=0; 
                int        i = ix + iiy + iiz;
                const int j = i;
                int i0 = i-ns.x;
                int i1 = i+ns.x;
                
                val += 
                    + Gs[i0-1]*B011 + Gs[i0]*B001 + Gs[i0+1]*B011
                    + Gs[i -1]*B001 + Gs[i ]*B000 + Gs[i +1]*B001
                    + Gs[i1-1]*B011 + Gs[i1]*B001 + Gs[i1+1]*B011;

                i   -= nxy;    i0 = i-ns.x;  i1 = i+ns.x;
                val +=    
                    + Gs[i0-1]*B111 + Gs[i0]*B011 + Gs[i0+1]*B111
                    + Gs[i -1]*B011 + Gs[i ]*B001 + Gs[i +1]*B011
                    + Gs[i1-1]*B111 + Gs[i1]*B011 + Gs[i1+1]*B111;  

                i   += 2*nxy;  i0 = i-ns.x;  i1 = i+ns.x;
                val +=
                    + Gs[i0-1]*B111 + Gs[i0]*B011 + Gs[i0+1]*B111
                    + Gs[i -1]*B011 + Gs[i ]*B001 + Gs[i +1]*B011
                    + Gs[i1-1]*B111 + Gs[i1]*B011 + Gs[i1+1]*B111;        
            
                //val +=  Gs[i ];
                double err = Es[j] - val;
                //double err = val;
                //double d =  Gs[i ];
                //printf( "[%i,%i,%i] %g %g %g \n", ix,iy,iz, d, val, Es[i] );
                ps[j] = err;
                Ws[j] = err;

                // if(itr==0){
                //     printf( "Gs[%i,%i,%i] G %g E %g p %g \n",    ix,iy,iz,    Gs[i], Es[i], ps[i] );
                // }
            }
        }
    }

    // --- distribute variational derivatives of approximation error
    for(int iz=0; iz<ns.z; iz++){
        int iiz = iz*nxy;
        for(int iy=0; iy<ns.y; iy++){
            int iiy = iy*ns.x;
            for(int ix=0; ix<ns.x; ix++){
                double val=0; 
                int i  = ix + iiy + iiz;
                const int j = i;
                int i0 = i-ns.x;
                int i1 = i+ns.x;


                if(  
                    (ix>0)&&(ix<(ns.x-1)) &&
                    (iy>0)&&(ix<(ns.y-1)) &&
                    (iz>0)&&(ix<(ns.z-1))           
                ){

                val += 
                    + ps[i0-1]*B011 + ps[i0]*B001 + ps[i0+1]*B011
                    + ps[i -1]*B001 + ps[i ]*B000 + ps[i +1]*B001
                    + ps[i1-1]*B011 + ps[i1]*B001 + ps[i1+1]*B011;

                i   -= nxy;    i0 = i-ns.x;  i1 = i+ns.x;
                val +=    
                    + ps[i0-1]*B111 + ps[i0]*B011 + ps[i0+1]*B111
                    + ps[i -1]*B011 + ps[i ]*B001 + ps[i +1]*B011
                    + ps[i1-1]*B111 + ps[i1]*B011 + ps[i1+1]*B111; 
                i   += 2*nxy;  i0 = i-ns.x;  i1 = i+ns.x;
                val +=
                    + ps[i0-1]*B111 + ps[i0]*B011 + ps[i0+1]*B111
                    + ps[i -1]*B011 + ps[i ]*B001 + ps[i +1]*B011
                    + ps[i1-1]*B111 + ps[i1]*B011 + ps[i1+1]*B111; 

                //val = ps[i];
                //printf( "[%i,%i,%i] %g \n", ix,iy,iz, val );
                fs[j] = val;
                // if(itr==0){
                //     printf( "Gs[%i,%i,%i] G %g p %g f %g\n", ix,iy,iz,  Gs[i], ps[i], fs[i]   );
                // }

                }else{
                    fs[j] = 0;
                    Gs[j] = 0;
                }

            }
        }
    }

}


__attribute__((hot)) 
int fit3D( const Vec3i ns, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1 ){
    printf( "Bspline::fit3D() ns(%i,%i,%i) \n", ns.x,ns.y,ns.z  );
    const double F2max = Ftol*Ftol;
    const int nxy  = ns.x*ns.y;
    const int nxyz = ns.x*ns.y*ns.z;
    double* ps = new double[nxyz];
    double* fs = new double[nxyz];
    double* vs = new double[nxyz];
    int itr=0;
    //while(false){
    for(itr=0; itr<nmaxiter; itr++){
        getVariations3D( ns, Gs, Es, Ws, fs, ps );
        // --- move
        double vf = 0.0;
        double ff = 0.0;
        for(int i=0; i<nxyz; i++){
            ff += fs[i]*fs[i];
            vf += vs[i]*fs[i];
        }
        printf( "|F[%i]|=%g  <v|f>=%g \n",itr,  sqrt(ff), vf );
        //printf( "p=%i v=%g f=%g \n",  Gs[1], vs[1], fs[1] );
        if(ff<F2max){ break; }
        if(vf<0){ 
            //printf( "<v|f>[%i]=%g < 0  => v=0 \n", itr, vf );
            for(int i=0; i<nxyz; i++){ vs[i]=0; }; 
        }
        for(int i=0; i<nxyz; i++){
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
int fit3D_CG( const Vec3i ns, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1 ){
    printf( "Bspline::fit3D() ns(%i,%i,%i) \n", ns.x,ns.y,ns.z  );
    const double F2max = Ftol*Ftol;
    const int nxy  = ns.x*ns.y;
    const int nxyz = ns.x*ns.y*ns.z;
    double* ps = new double[nxyz];
    double* fs = new double[nxyz];
    double* vs = new double[nxyz];
    int itr=0;
    //while(false){
    for(itr=0; itr<nmaxiter; itr++){
        getVariations3D( ns, Gs, Es, Ws, fs, ps );
        // --- move
        double vf = 0.0;
        double ff = 0.0;
        for(int i=0; i<nxyz; i++){
            ff += fs[i]*fs[i];
            vf += vs[i]*fs[i];
        }
        printf( "|F[%i]|=%g  <v|f>=%g \n",itr,  sqrt(ff), vf );
        //printf( "p=%i v=%g f=%g \n",  Gs[1], vs[1], fs[1] );
        if(ff<F2max){ break; }
        if(vf<0){ 
            //printf( "<v|f>[%i]=%g < 0  => v=0 \n", itr, vf );
            for(int i=0; i<nxyz; i++){ vs[i]=0; }; 
        }
        for(int i=0; i<nxyz; i++){
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


} // namespace SplineBcub{

#endif



