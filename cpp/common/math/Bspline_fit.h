
#ifndef  Bspline_fit_h
#define  Bspline_fit_h

#include "Vec3.h"
#include "quaternion.h"
//#include "CG.h"
#include "globals.h"

#include "Bspline.h"

namespace Bspline{

__attribute__((hot)) 
Vec3d move( double dt, int n, double* gs, double* fs, double* vs ){
    //constexpr const int m = 2;
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
Vec3d move_GD( double dt, int n, double* gs, double* fs ){
    double ff = 0.0;
    for(int i=0; i<n; i++){
        double f = fs[i];
        gs[i]   += f*dt;
        ff      += f*f;
    }
    return Vec3d{0.0,ff,0.0};
}


__attribute__((hot)) 
int fit1D_old( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1 ){
    if(verbosity>1)printf("Bspline::fit1D_old() \n");
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
    double err2sum=0;
    for(itr=0; itr<nmaxiter; itr++){
        err2sum=0.0;
        // --- evaluate current spline
        for(int i=0; i<n; i++){
        //for(int i=1; i<n-1; i++){
            double p = Gs[i]*B0;
            if(i>0  )[[likely]]{ p += B1*Gs[i-1]; }
            if(i<n-1)[[likely]]{ p += B1*Gs[i+1]; }

            //p += B1*Gs[i-1];
            //p += B1*Gs[i+1];

            //bErr|=checkNumRange( i, T val, T min, T max, const char* pre, bool bPrint=true, bool bExit=false ){
            //printf( "p[%i]=%g Gs=%g\n", i, p, Gs[i] );
            //bErr|=  checkNumRange( i, p, -1.e+6, +1.e+6, "p=B*Gs", true, true );
            ps[i] = p;   //Ws[i] = p;
            fs[i] = 0;
        }
        // --- evaluate variatiaonal derivatives
        //for(int i=0; i<n; i++){
        for(int i=1; i<n-1; i++){
            double dp  = Es[i] - ps[i];
            if( Ws ){ dp *= Ws[i]; } // weighting 
            err2sum += dp*dp;
            fs[i] += dp*B0;
            //bErr|=  checkNumRange( i, dp, -1.e+6, +1.e+6, "dp=Es-ps", true, true );
            //if(i>0  )[[likely]]{ fs[i-1] += B1*dp; }
            //if(i<n-1)[[likely]]{ fs[i+1] += B1*dp; }
            fs[i-1] += B1*dp;
            fs[i+1] += B1*dp;
        }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "Bspline::fit1D_old() iter=%i err=%g \n", itr, sqrt(err2sum) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
double getVariations1D( const int n, const double* Gs, const double* Es, const double* Ws, double* fs, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    // --- evaluate current spline
    for(int i=0; i<n; i++){ fs[i]=0; ps[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int i=1; i<n-1; i++){
        double val = Gs[i-1]*B1 + Gs[i]*B0 + Gs[i+1]*B1;
        double err = Es[i] - val;
        if(Ws){ err*=Ws[i]; }
        err2sum += err*err;
        ps[i] = err;
        //Ws[j] = err;
    }
    // --- distribute variational derivatives of approximation error
    for(int i=0; i<n; i++){
        if( (i>0)&&(i<(n-1)) ) [[likely]] {
            double val = ps[i-1]*B1 + ps[i]*B0 + ps[i+1]*B1;
            fs[i] = val;
        }else{
            fs[i] = 0;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
double getVariations1D_half( const int n, const double* Gs, const double* Es, const double* Ws, double* fs, double* ers ){
    constexpr double B00 =2.0/3.0;   // 0.66666666666666666666666666666667
    constexpr double B05 =23.0/48.0; // 0.47916666666666666666666666666667
    constexpr double B10 =1.0/6.0;   // 0.16666666666666666666666666666667
    constexpr double B15 =1.0/48.0;  // 0.020833333333333333333333333333333
    
    int n2 = 2*n;
    // --- evaluate current spline
    for(int i=0; i<n;  i++){ fs[i]=0;   }
    for(int i=0; i<n2; i++){ ers[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int i=2; i<n-2; i++){
        double e0 = Gs[i-1]*B10 + Gs[i]*B00 + Gs[i+1]*B10;
        double e1 = Gs[i-1]*B15 + Gs[i]*B05 + Gs[i+1]*B05 + Gs[i+2]*B15;
        //double e1 = Gs[i-2]*B15 + Gs[i-1]*B05 + Gs[i]*B05 + Gs[i+1]*B15;

        int i2=i*2;
        double d0 = Es[i2  ] - e0;   ers[i2  ] = d0;
        double d1 = Es[i2+1] - e1;   ers[i2+1] = d1;
        
        err2sum += e0*e0 + e1*e1;
        
        //Ws[j] = err;
    }
    // --- distribute variational derivatives of approximation error
    for(int i=0; i<n; i++){
        if( (i>0)&&(i<(n-1)) ) [[likely]] {
            int i2=i*2;
            //double f = ers[i2-3]*B15 + ers[i2-2]*B10 + ers[i2-1]*B05 + ers[i2]*B00 + ers[i2+1]*B05 + ers[i2+2]*B10 + ers[i2+3]*B15;

            double f = ers[i2-2]*B10 + ers[i2]*B00 + ers[i2+2]*B10;

            f += ers[i2-1]*B15 + ers[i2+1]*B15;
            //if(i2>=3){ f+=ers[i2-3]*B15; }
            //if(i2<n ){ f+=ers[i2+3]*B15; }

            fs[i]    = f;
        }else{
            fs[i] = 0;
        }
    }
    return err2sum;
}

__attribute__((hot)) 
int fit1D( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bHalf=false ){
    if(verbosity>1)printf("Bspline::fit1D() !!!!!!\n");
    const double F2max = Ftol*Ftol;
    int n2=n; if(bHalf){ n2=2*n; }
    double* ps = new double[n2];
    double* fs = new double[n];
    double* vs = new double[n];
    int itr=0;
    for(int i=0; i<n; i++){ vs[i]=0.0; } // clear velocity
    //if( ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) && ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) ){ printf("ERROR in fit1D()=> exit\n"); exit(0); };
    bool bErr = false;
    double err2sum=0;
    for(itr=0; itr<nmaxiter; itr++){        
        if(bHalf){ err2sum = getVariations1D_half( n, Gs, Es, Ws, fs, ps ); }
        else     { err2sum = getVariations1D(      n, Gs, Es, Ws, fs, ps ); }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "Bspline::fit1D() iter=%i err=%g \n", itr, sqrt(err2sum) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}

__attribute__((hot)) 
int fit1D_EF( const double dg, const int n, double* Gs,  Vec2d* fes, Vec2d* Ws, double Ftol, int nmaxiter=1000, double dt=0.1 ){
    if(verbosity>1)printf("Bspline::fit1D_EF() \n");
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
    double err=0;
    for(itr=0; itr<nmaxiter; itr++){
        err=0.0;
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
            err+= dp.norm2();
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
        if(verbosity>2)printf( "|F[%i]|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    if(verbosity>1)printf( "Bspline::fit1D_EF() iter=%i err=%g \n", itr, err );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}


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

__attribute__((hot)) 
double error_iz( int iz, const Vec3i ns, const double* Gs, const double* Es, const double* Ws, double* ps, const double Esafe=-1.0 ){
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
void force_iz( int iz, const Vec3i ns, const double* ps, double* fs ){
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
double error_iz_pbc( int iz, const Vec3i ns, const double* Gs, const double* Es, double* Ws, double* ps, const double Esafe=-1.0 ){
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
void force_iz_pbc( int iz, const Vec3i ns, const double* ps, double* fs ){
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



//                               [0]=0.0          [1]=0.5          [2]=1.0          [3]=1.5
constexpr const double  Bs[4]{  basis(0.0).y,  basis(0.5).y,   basis(0.0).x,   basis(0.5).x };  // values B-spline centered at 0 at points {0.0, 0.5, 1.0, 1.5}
constexpr const double dBs[4]{ dbasis(0.0).y, dbasis(0.5).y,  dbasis(0.0).x,  dbasis(0.5).x };  // 1st derivatives B-spline centered at 0 at points {0.0, 0.5, 1.0, 1.5}
constexpr const double ddBs[4]{ddbasis(0.0).y,ddbasis(0.5).y, ddbasis(0.0).x, ddbasis(0.5).x };  // 2nd derivatives B-spline centered at 0 at points {0.0, 0.5, 1.0, 1.5}

/**
 * Computes the value and derivative of a B-spline basis function at a given point.
 *
 * @param x The input point at which to evaluate the B-spline basis function.
 * @return A Vec2d containing the value and derivative of the B-spline basis function at the given point.
 */
inline Vec2d valDval(double x){
    double s;
    if(x>0){ s=1.; }else{ s=-1.; x=-x; };
    //                       value           derivative
    if (x<1){ return Vec2d{ basis(x).y,   s*dbasis(x).y   }; }
    else    { return Vec2d{ basis(x-1).x, s*dbasis(x-1).x }; }
}

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



__attribute__((hot)) 
void force_iz_pbc_refForce( int iz, const Vec3i ns, const double* ps, double* fs,      int ndi, Vec3i* dis, double* Gs, const Quat4i* xqis, const Quat4i* yqis ){
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
double getVariations3D( const Vec3i ns, const double* Gs, const double* Es, double* Ws, double* fs, double* ps, bool bPBC, const double Esafe=-1.0 ){
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
int fit2D( const Vec2i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false ){
    //if(verbosity>1)
    printf( "Bspline::fit2D() ns(%i,%i) Ftol=%g dt=%g nmaxiter=%i bPBC=%i \n", ns.x,ns.y, Ftol, dt, nmaxiter, bPBC );
    const double F2max = Ftol*Ftol;
    const int nxy  = ns.x*ns.y;
    double* ps = new double[nxy];
    double* fs = new double[nxy];
    double* vs = new double[nxy];
    int itr=0;
    //while(false){
    double err;    
    Vec3d  cfv;
    for(int i=0; i<nxy; i++){ vs[i]=0; };
    for(itr=0; itr<nmaxiter; itr++){
        //for(int i=0; i<nxy; i++){ fs[i]=0; };  // Not necessary - force is cleared inside getVariations2D
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

__attribute__((hot)) 
int fit3D( const Vec3i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false, bool bInitGE=false, const double Esafe=-1.0 ){
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
int fit3D_omp( const Vec3i ns, double* Gs, const double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bPBC=false, bool bInitGE=false, const double Esafe=-1.0 ){
    //if(verbosity>1)
    //printf( "Bspline::fit3D_omp() ns(%i,%i,%i) bPBC=%i dt=%g Ftol=%g nmaxiter=%i bInitGE=%i \n", ns.x,ns.y,ns.z, bPBC, dt, Ftol, nmaxiter, bInitGE );
    const int nxy  = ns.x*ns.y;
    const int nxyz = nxy*ns.z;
    const double F2max = Ftol*Ftol;
    double* ps = new double[nxyz];
    double* fs = new double[nxyz];
    double* vs = new double[nxyz];

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
            vs[i] += fs[i]*dt;
            Gs[i] += vs[i]*dt;
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



