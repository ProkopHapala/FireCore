
#ifndef  Bspline_fit_3D_h
#define  Bspline_fit_3D_h

#include "Bspline_fit_2D.h"
#include "Vec3.h"
#include "quaternion.h"
//#include "CG.h"
#include "globals.h"

#include "Bspline.h"
#include "Bspline_fit.h"
#include "testUtils.h"

namespace Bspline{

    __attribute__((hot)) 
inline double error_iz( int iz, const Vec3i ns, const double* Gs, const double* Es, const double* Ws, double* ps, const double Esafe=-1.0 ){
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
            const double Ei = Es[i];
            double err = Ei - val;
            //if((ix==6)&&(iy==10)&&(iz==10)){  printf("getVariations3D_mod()[%i,%i,%i] E=%g val=%g err=%g \n", ix,iy,iz, Es[i], val, err ); }
            //Ws[i] = err;
            if(Esafe>0){ err /= ( fabs(Ei) + Esafe ); }
            if(Ws){ err*=Ws[i]; }
            err2sum += err*err;
            ps[i]    = err;
            //Ws[i] = err;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
inline void force_iz( int iz, const Vec3i ns, const double* ps, double* fs ){
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

/// Calculates the point-wise approximation error `ps` of the cubic B-spline for given spline coefficients `Gs` for a single z-slice, taking into account periodic boundary conditions.
///
/// @param iz The index of the z-slice to process.
/// @param ns The dimensions of the 3D grid {nx, ny, nz}
/// @param Gs The spline coefficients on the 3D grid.
/// @param Es The target values on the 3D grid.
/// @param Ws The weights for the target values (optional).
/// @param ps The output pre-computed approximation errors for the 3D grid.
/// @param Esafe A small value added to the denominator to avoid division by zero (optional).
/// @return The total sum of square error for the z-slice.
__attribute__((hot)) 
inline double error_iz_pbc( int iz, const Vec3i ns, const double* Gs, const double* Es, double* Ws, double* ps, const double Esafe=-1.0 ){
    //printf("error_iz_pbc(iz=%i) \n");
    constexpr double B0=2.0/3.0; // value of cubic B-spline basis at x= 0.0
    constexpr double B1=1.0/6.0; // value of cubic B-spline basis at x=-1.0 or 1.0
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

            double  val  = assembleBound2D_pbc( B000,B001,B011, Gs    , i, ibx,idx,  iby,idy ); 
            if(zlo) val += assembleBound2D_pbc( B001,B011,B111, Gs-nxy, i, ibx,idx,  iby,idy ); 
            if(zhi) val += assembleBound2D_pbc( B001,B011,B111, Gs+nxy, i, ibx,idx,  iby,idy );

            const double Ei = Es[i];
            double err = Ei - val;

            if(Esafe>0){ err /= ( fabs(Ei) + Esafe ); }

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

/**
 * @brief Calculates the variational derivatives of the total mean-square error of the cubic B-spline fit with respect to the spline coefficients on the grid `Gs`, for a single z-slice using the point-wise approximation errors stored in `ps`
 *
 * @param iz The index of the z-slice to process.
 * @param ns The dimensions of the 3D grid {nx, ny, nz}
 * @param ps The pre-computed approximation errors for the 3D grid.
 * @param fs The output variational derivatives (3D array). i.e. derivative of MSE with respect to Gs 
 */
__attribute__((hot)) 
inline void force_iz_pbc( int iz, const Vec3i ns, const double* ps, double* fs ){
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

            double  val  = assembleBound2D_pbc( B000,B001,B011, ps    , i, ibx,idx,  iby,idy ); 
            if(zlo) val += assembleBound2D_pbc( B001,B011,B111, ps-nxy, i, ibx,idx,  iby,idy ); 
            if(zhi) val += assembleBound2D_pbc( B001,B011,B111, ps+nxy, i, ibx,idx,  iby,idy );    

            //if(Esafe>0){ val /= ( fabs(Ei) + Esafe ); }
            //val*=-1;
            fs[i] = val;
        }
    }
    //return err2sum;
}



// ========== Regularization for fitting

/**
 * Computes the variational derivative of Force difference F(pi)-F(pi+d) represented by B-spline with respect to B-spline basis function centered at grid pint point pi.
 *
 * @param d The displacement vector from the grid point pi.
 * @return The vector of variational derivatives of force components dFx,dFy,dFz 
 */
inline Vec3d varG_ForceMidpoint( const Vec3d d ){
    double dB  = dbasis(0.0).y;
    double  B  =  basis(0.0).y;
    double dBB = B*B*dB;
    Vec2d Bx  = valDval( d.x );
    Vec2d By  = valDval( d.y );
    Vec2d Bz  = valDval( d.z );
    return Vec3d{
       dBB  -  Bx.y * By.x * Bz.x  , 
       dBB  -  Bx.x * By.y * Bz.x  ,
       dBB  -  Bx.x * By.x * Bz.y 
    };
}

/**
 * Computes the variational derivative of Force difference F(pi)-F(pi+di*0.5) represented by B-spline with respect to B-spline basis function centered at grid point point pi. 
 *
 * @param di The discretized displacement vector from the grid point pi. The true displacement is d=di*0.5.
 * @return The vector of variational derivatives of force components dFx,dFy,dFz
 */
inline Vec3d varG_ForceMidpoint( const Vec3i di ){
    double dB = dBs[0];
    double  B =  Bs[0];
    double dBB = B*B*dB;
    int sx = (di.x>>31); int ix = (di.x^sx)+(sx<0); 
    int sy = (di.y>>31); int iy = (di.y^sy)+(sy<0);
    int sz = (di.z>>31); int iz = (di.z^sz)+(sz<0);
    //int s = sx*sy*sz;
    Vec3d B_ {  Bs[di.x],     Bs[di.y],     Bs[di.z]     };
    Vec3d dB_{ dBs[di.x]*sx, dBs[di.y]*sy, dBs[di.z]*sz  };
    return Vec3d{
       dBB  -  dB_.x * B_.y *B_.z  , 
       dBB  -   B_.x *dB_.y *B_.z  ,
       dBB  -   B_.x * B_.y *dB_.z 
    };
}

/**
 * Computes the variational derivative of the sum of force differences sum_i{ F(p) - F(p + d_i*0.5) } represented by a B-spline with respect to the B-spline basis function centered at the grid point p.
 *
 * @param p    The position of the grid-point on which the relevant basis function G_i is centerd.
 * @param ndi  The number of displacement vectors to sum
 * @param dis  The array of displacement vectors.
 * @param Gs   The B-spline coefficients on the 3D grid ( which respect to which we calculate the derivative)
 * @param n    The dimensions of the 3D grid {nx, ny, nz}.
 * @param xqis The array of 4 index offsets to efficieltly evaluate periodic boundary conditions along x-axis.
 * @param yqis The array of 4 index offsets to efficieltly evaluate periodic boundary conditions along y-axis.
 * @return The vector of variational derivatives of the force components dFx, dFy, dFz.
 */
inline double regularizationForceMidpoint( Vec3d p, int ndi, Vec3i* dis, const double* Gs, Vec3i n, const Quat4i* xqis, const Quat4i* yqis ){
    const Vec3d f_ref = fe3d_pbc( p, n, Gs, xqis,yqis ).f;      //   Force at point colocated with the grid point  G_i
    //Vec3d     reg_var = Vec3dZero;
    double reg_var=0;
    for(int i=0; i<ndi; i++){
        const Vec3i& di = dis[i];
        const Vec3d  vG = varG_ForceMidpoint( di );  // variational deriv of basis function B_i  at point  p
        const Vec3d  p_ = p + Vec3d{di.x*0.5, di.y*0.5, di.z*0.5};
        const Vec3d  dF = f_ref - fe3d_pbc( p_, n, Gs, xqis,yqis ).f;
        //reg_var += vG*dF;
        reg_var += 2.0 * vG.dot(dF);
    }
    return reg_var;
}

/**
 * @brief Computes the variational derivatives of the force differences represented by a B-spline with respect to the B-spline basis function centered at the grid point summed over number of points offset with respect to the grid point.
 *
 * @param iz The current z-index of the grid.
 * @param ns The dimensions of the 3D grid {nx, ny, nz}.
 * @param ps The B-spline coefficients on the 3D grid.
 * @param fs The output variational derivatives of the force components.
 * @param ndi The number of displacement vectors to sum.
 * @param dis The array of displacement vectors.
 * @param Gs The B-spline coefficients on the 3D grid (with respect to which we calculate the derivative).
 * @param xqis The array of 4 index offsets to efficiently evaluate periodic boundary conditions along the x-axis.
 * @param yqis The array of 4 index offsets to efficiently evaluate periodic boundary conditions along the y-axis.
 */
__attribute__((hot)) 
inline void force_iz_pbc_regForce( int iz, const Vec3i ns, const double* ps, double* fs,      int ndi, Vec3i* dis, double* Gs, const Quat4i* xqis, const Quat4i* yqis ){
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

            double  val  = assembleBound2D_pbc( B000,B001,B011, ps    , i, ibx,idx,  iby,idy ); 
            if(zlo) val += assembleBound2D_pbc( B001,B011,B111, ps-nxy, i, ibx,idx,  iby,idy ); 
            if(zhi) val += assembleBound2D_pbc( B001,B011,B111, ps+nxy, i, ibx,idx,  iby,idy );    

            if( (iz>3) && (iz<ns.z-4) ){
                val    += regularizationForceMidpoint( Vec3d{(double)ix,(double)iy,(double)iz}, ndi, dis, Gs, ns, xqis, yqis );
            }

            //if(Esafe>0){ val /= ( fabs(Ei) + Esafe ); }
            //val*=-1;
            fs[i] = val;
        }
    }
    //return err2sum;
}

/**
 * @brief Calculates the variational derivatives of the total mean-square error of the cubic B-spline fit of the reference (Es) with respect to the spline coefficients on the grid (Gs).
 *
 * @details Mathematical derivation:
 * Total mean-square error MSE = sum_i{ (E_i - y_i)^2 }/ntot = sum_i{ p_i^2 }/ntot  
 * spline approximation: y_i = y[ix,iy,iz] = sum_(jx,jy,jz){  B(dx[jx]) * B(dx[jx]) * B(dx[jx]) * G[ix+jx,iy+jy,iz+jz]
 * p_i = E_i - y_i = G[ix,iy,iz]   -  sum_(jx,jy,jz){  B(dx[jx]) * B(dx[jx]) * B(dx[jx]) * G[ix+jx,iy+jy,iz+jz]
 * d(MSE)/dG_i =  d(MSE)/dG[ix,iy,iz] =   sum_(jx,jy,jz) {  2 * p[ix+jx,iy+jy,iz+jz] * B(dx[jx]) * B(dx[jx]) * B(dx[jx]) }
 * 
 * @param ns dimensions of the grid {nx,ny,nz}.
 * @param Gs spline coefficients (3D array)
 * @param Es reference values (3D array)
 * @param Ws weights for the reference (3D array) (currently not used)
 * @param fs output variational derivatives (3D array)
 * @param ps output approximation errors. (3D array)
 * @param bPBC Whether to use periodic boundary conditions.
 * @param Esafe A safety value to avoid division by zero. If Esafe>0 we use relative error metric err/=(fabs(Ei)+Esafe) instead of absolute error to calculate RMS.
 * @return The sum of the squared approximation errors.
 */
static double getVariations3D( const Vec3i ns, const double* Gs, const double* Es, double* Ws, double* fs, double* ps, bool bPBC, const double Esafe=-1.0 ){
    /*
    We calculate the variational derivatives of the error of the Bspline fit of the refferentce (Es) with respect to spline coefficiet on the grid (Gs)
    */
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
        if(bPBC){ err2sum += error_iz_pbc( iz, ns, Gs, Es, Ws, ps, Esafe ); }
        else    { err2sum += error_iz    ( iz, ns, Gs, Es, Ws, ps, Esafe ); }            
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
inline int fit3D( const Vec3i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false, bool bInitGE=false, const double Esafe=-1.0 ){
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
        err = getVariations3D( ns, Gs, Es, Ws, fs, ps, bPBC, Esafe );
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
inline int fit3D_omp( const Vec3i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.3, double damp=0.0, bool bPBC=false, bool bInitGE=false, const double Esafe=-1.0 ){
    //if(verbosity>1)
    //printf( "Bspline::fit3D_omp() ns(%i,%i,%i) bPBC=%i dt=%g damp=%g Ftol=%g nmaxiter=%i bInitGE=%i \n", ns.x,ns.y,ns.z, bPBC, dt, damp, Ftol, nmaxiter, bInitGE );
    const int nxy  = ns.x*ns.y;
    const int nxyz = nxy*ns.z;
    const double F2max = Ftol*Ftol;
    double* ps = new double[nxyz];
    double* fs = new double[nxyz];
    double* vs = new double[nxyz];

    double cdamp = 1-damp;

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
            if(bPBC){ err2sum += error_iz_pbc( iz, ns, Gs, Es, Ws, ps, Esafe ); }
            else    { err2sum += error_iz    ( iz, ns, Gs, Es, Ws, ps, Esafe ); }            
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
            double vi=vs[i];
            vi    *= cdamp;
            vi    += fs[i]*dt;
            Gs[i] += vi*dt;
            vs[i] = vi;
        }
        #pragma omp single
        {
            //cfv = move_GD( dt, nxyz, Gs, fs );
            //if(verbosity>2)
            //printf( "|F[%i]|=%7.2e Error=%7.2e \n",itr,sqrt(ff), sqrt(err2sum) );
            ///printf( "|F[%i]|=%g cos(f,v)=%g Error=%g \n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z), sqrt(err) );
            niterdone=itr;
            if(ff<F2max){ 
                //printf( "Bspline::fit3D_omp() DONE |F[%i]|=%7.2e Error=%7.2e \n",itr,sqrt(ff), sqrt(err2sum) );   
                itr=nmaxiter+1; 
            };
            itr++;
        }
    }
    }
    //if(verbosity>1)
    double t = getCPUticks()-t0; printf( "Bspline::fit3D_omp() RMS=%7.2e |F|=%7.2e niter(%4i) nxyz=(%8i) time %g[GTicks] %g[tick/(nxyz*iter)] \n", sqrt(err2sum), sqrt(ff), niterdone,nxyz, t*1e-9, t/(nxyz*niterdone) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return niterdone;
}

} // namespace SplineBcub{

#endif



