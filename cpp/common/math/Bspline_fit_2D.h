
#ifndef  Bspline_fit_2D_h
#define  Bspline_fit_2D_h

#include "Vec3.h"
#include "quaternion.h"
//#include "CG.h"
#include "globals.h"

#include "Bspline.h"
#include "Bspline_fit.h"

namespace Bspline{

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

/**
 * @brief  a 2D B-spline value at the given index, with periodic boundary conditions.
 *
 * @details function computes the value of a 2D B-spline at the given index `i`, using the  provided B-spline coefficients `Gs` and the B-spline basis function coefficients `B00`, `B01`, and `B11`. The function assumes periodic boundary conditions in both the x and y dimensions.
 *
 * @param B00 value of the B-spline basis function B0 * B0 
 * @param B01 value of the B-spline basis function B0 * B1
 * @param B11 value of the B-spline basis function B1 * B1
 * @param Gs  array of B-spline coefficients
 * @param i   index of the B-spline value to assemble (already folded from ix,iy (iz)
 * @param ibx x-offset for the left  neighbor in the x-dimension
 * @param idx x-offset for the right neighbor in the x-dimension
 * @param iby y-offset for the lower neighbor in the y-dimension
 * @param idy y-offset for the upper neighbor in the y-dimension
 * @return The assembled B-spline value at the given index.
 */
__attribute__((pure)) 
__attribute__((hot)) 
inline double assembleBound2D_pbc( const double B00, const double B01, const double B11, const double* Gs, int i, int ibx, int idx, int iby, int idy ){
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
            double val = assembleBound2D_pbc( B00,B01,B11,Gs,i,ibx,idx,iby,idy );
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
            double f= assembleBound2D_pbc( B00,B01,B11,ps,i,ibx,idx,iby,idy );
            //Ws[i] = f;
            fs[i] = f;
        }
    }
    return err2sum;
}

// ========== Regularization for fitting


/**
 * Computes the variational derivative of Force difference F(pi)-F(pi+d) represented by B-spline with respect to B-spline basis function centered at grid pint point pi.
 *
 * @param d The displacement vector from the grid point pi.
 * @return The vector of variational derivatives of force components dFx,dFy,dFz 
 */
inline Vec2d varG_2D_ForceMidpoint( const Vec2d d ){
    double dB  = dbasis(0.0).y;
    double  B  =  basis(0.0).y;
    double dBB = B*B*dB;
    Vec2d Bx  = valDval( d.x );
    Vec2d By  = valDval( d.y );
    return Vec2d{
       dBB  -  Bx.y * By.x, 
       dBB  -  Bx.x * By.y,
    };
}

/**
 * Computes the variational derivative of Force difference F(pi)-F(pi+di*0.5) represented by B-spline with respect to B-spline basis function centered at grid point point pi. 
 *
 * @param di The discretized displacement vector from the grid point pi. The true displacement is d=di*0.5.
 * @return The vector of variational derivatives of force components dFx,dFy,dFz
 */
inline Vec2d varG_2D_ForceMidpoint( const Vec2i di ){
    double dB = dBs[0];
    double  B =  Bs[0];
    double dBB = B*dB;
    int sx = (di.x>>31); int ix = (di.x^sx)+(sx<0); 
    int sy = (di.y>>31); int iy = (di.y^sy)+(sy<0);
    //int s = sx*sy*sz;
    Vec2d B_ {  Bs[di.x],     Bs[di.y],   };
    Vec2d dB_{ dBs[di.x]*sx, dBs[di.y]*sy };
    return Vec2d{
       dBB  -  dB_.x * B_.y , 
       dBB  -   B_.x *dB_.y ,
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
inline double regularizationForceMidpoint_2D( Vec2d p, int ndi, Vec2i* dis, const double* Gs, Vec2i n, const Quat4i* xqis, const Quat4i* yqis ){
    const Vec2d f_ref = fe2d_pbc_macro( p, n, Gs, xqis,yqis ).xy();      //   Force at point colocated with the grid point  G_i
    //Vec3d     reg_var = Vec3dZero;
    double reg_var=0;
    for(int i=0; i<ndi; i++){
        const Vec2i& di = dis[i];
        const Vec2d  vG = varG_2D_ForceMidpoint( di );  // variational deriv of basis function B_i  at point  p
        const Vec2d  p_ = p + Vec2d{di.x*0.5, di.y*0.5};
        const Vec2d  dF = f_ref - fe2d_pbc_macro( p_, n, Gs, xqis,yqis ).xy();
        //reg_var += vG*dF;
        reg_var += 2.0 * vG.dot(dF);
    }
    return reg_var;
}

__attribute__((hot)) 
double getVariations2D_regForce( const Vec2i ns, double* Gs,  const double* Es, double* Ws, double* fs, double* ps,    int ndi, Vec3i* dis, const Quat4i* xqis, const Quat4i* yqis ){
    //printf("getVariations2D_pbc()\n");
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    constexpr double B00=B0*B0;
    constexpr double B01=B0*B1;
    constexpr double B11=B1*B1;
    double err2sum=0;
    int nxy=ns.x*ns.y;
    for(int iy=0; iy<ns.y; iy++){
        //printf("getVariations2D_pbc() iy=%i \n", iy );
        const int iiy = iy*ns.x;
        const int iby = biwrap(iy,ns.y)*ns.x;
        const int idy = diwrap(iy,ns.y)*ns.x;
        for(int ix=0; ix<ns.x; ix++){
            const int i   = ix + iiy;
            const int ibx = biwrap(ix,ns.x);
            const int idx = diwrap(ix,ns.x);
            double val = assembleBound2D_pbc( B00,B01,B11,Gs,i,ibx,idx,iby,idy );
            double err = Es[i] - val;
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
            double f= assembleBound2D_pbc( B00,B01,B11,ps,i,ibx,idx,iby,idy );
            //Ws[i] = f;
            fs[i] = f;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
int fit2D( const Vec2i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false, bool bRegForce=false ){
    //if(verbosity>1)
    printf( "Bspline::fit2D() ns(%i,%i) Ftol=%g dt=%g nmaxiter=%i bPBC=%i \n", ns.x,ns.y, Ftol, dt, nmaxiter, bPBC );
    const double F2max = Ftol*Ftol;
    const int nxy  = ns.x*ns.y;
    double* ps = new double[nxy];
    double* fs = new double[nxy];
    double* vs = new double[nxy];

    if( bRegForce ){

    }


    int itr=0;
    //while(false){
    double err;    
    Vec3d  cfv;
    for(int i=0; i<nxy; i++){ vs[i]=0; };
    for(itr=0; itr<nmaxiter; itr++){
        //for(int i=0; i<nxy; i++){ fs[i]=0; };  // Not necessary - force is cleared inside getVariations2D
        // if(bRegForce){
        //     err = getVariations2D_regForce( ns, Gs, Es, Ws, fs, ps,    ndi, dis, xqis, yqis );
        // }else
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

} // namespace SplineBcub{

#endif



