
#ifndef  Bspline_fit_h
#define  Bspline_fit_h

#include "Vec3.h"
#include "quaternion.h"
//#include "CG.h"
#include "globals.h"

#include "Bspline.h"

namespace Bspline{

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
    //---------              value           derivative
    if (x<1){ return Vec2d{ basis(x).y,   s*dbasis(x).y   }; }
    else    { return Vec2d{ basis(x-1).x, s*dbasis(x-1).x }; }
}

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
double getVariations1D_pbc( const int n, const double* Gs, const double* Es, const double* Ws, double* fs, double* ps ){
    constexpr double B0=2.0/3.0;
    constexpr double B1=1.0/6.0;
    // --- evaluate current spline
    for(int i=0; i<n; i++){ fs[i]=0; ps[i]=0;  }
    // --- evaluate current spline (in order to evelauet approximation error)
    double err2sum=0;
    for(int i=0; i<n; i++){
            const int ib = biwrap(i,n);
            const int id = diwrap(i,n);
            double val = Gs[i+ib]*B1 + Gs[i]*B0 + Gs[i+id]*B1;
            //double val = Gs[i-1]*B1 + Gs[i]*B0 + Gs[i+1]*B1;
        double err = Es[i] - val;
        if(Ws){ err*=Ws[i]; }
        err2sum += err*err;
        ps[i] = err;
        //Ws[j] = err;
    }
    // --- distribute variational derivatives of approximation error
    for(int i=0; i<n; i++){
        const int ib = biwrap(i,n);
        const int id = diwrap(i,n);
        double val = ps[i+ib]*B1 + ps[i]*B0 + ps[i+id]*B1;
        //double val = ps[i-1]*B1 + ps[i]*B0 + ps[i+1]*B1;
        fs[i] = val;
    }
    return err2sum;
}


inline double regularizationForceMidpoint_1D( int n, const double* Gs, double* fs, double Kreg, const Quat4i* xqis ){
    // E_reg= sum_i Kreg*( f(x_i) - f(x_i+0.5) )^2 + Kreg*( f(x_i) - f(x_i-0.5) )^2  
    // dE_reg/dG_i = 2*Kreg*( f(x_i) - f(x_i+0.5) )*df(x_i)/dG_i + 2*Kreg*( f(x_i) - f(x_i-0.5) )*df(x_i)/dG_i
    double dB   =  dBs[0]; // derivative of B-spline basis function at point  p
    double dBfw =  dBs[1]; // derivative of B-spline basis function at point  p+0.5
    double dBbk = -dBs[1]; // derivative of B-spline basis function at point  p-0.5
    double f2reg=0;
    for(int i=0; i<n; i++){
        const double x = i;
        const double f_ref = fe1d_pbc_macro( x, n, Gs, xqis ).y;      //   Force at point colocated with the grid point  G_i
        double freg = 0;
        {  // backward point
            const double  vG = dB - dBbk;                                     // difference of B-spline basis functions at point  p and p-0.5
            const double  dF = f_ref - fe1d_pbc_macro( x-0.5, n, Gs, xqis).y; // differece of spline derivative (force) at point  p and p-0.5
            freg += 2.0 * vG * dF;
        }
        {   // forward point
            const double  vG = dB - dBfw;                                      // difference of B-spline basis functions at point  p and p+0.5
            const double  dF = f_ref - fe1d_pbc_macro( x+0.5, n, Gs, xqis).y; // differece of spline derivative (force) at point  p and p+0.5
            freg += 2.0 * vG * dF;
        }
        freg *= -Kreg;
        fs[i] += freg;
        f2reg += freg*freg;
    }
    return f2reg;
}


__attribute__((hot)) 
int fit1D( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, double Kreg=-1.0, bool bPBC=false ){
    //if(verbosity>1)
    printf("Bspline::fit1D() bPBC=%i Kreg=%g Ftol=%g nmaxiter=%i dt=%g \n", bPBC, Kreg, Ftol, nmaxiter, dt );
    const double F2max = Ftol*Ftol;
    double* ps = new double[n];
    double* fs = new double[n];
    double* vs = new double[n];
    bool bRegForce = Kreg>0;
    Quat4i xqis[4];
    if(bRegForce){ make_inds_pbc( n, xqis ); }
    int itr=0;
    for(int i=0; i<n; i++){ vs[i]=0.0; } // clear velocity
    //if( ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) && ckeckRange(n, 1, Gs, -1.e+6, +1.e+6, "Gs", true ) ){ printf("ERROR in fit1D()=> exit\n"); exit(0); };
    bool bErr = false;
    double err2sum=0;
    double f2reg=0;
    for(itr=0; itr<nmaxiter; itr++){        
        //if(bPBC){ err2sum = getVariations1D_pbc( n, Gs, Es, Ws, fs, ps ); }
        //else    { err2sum = getVariations1D    ( n, Gs, Es, Ws, fs, ps ); }
        if(bPBC){ err2sum = getVariations1D_pbc( n, Gs, Es, 0, fs, ps ); }
        else    { err2sum = getVariations1D    ( n, Gs, Es, 0, fs, ps ); }
        if( bRegForce ){
            //f2reg = regularizationForceMidpoint_1D( n, Gs, fs, Kreg, xqis );
            for(int i=0; i<n; i++){ Ws[i]=0.0; };
            f2reg = regularizationForceMidpoint_1D( n, Gs, Ws, Kreg, xqis );
            //for(int i=0; i<n; i++){ fs[i]+=Ws[i]; };
        }
        Vec3d cfv = move(dt,n,Gs,fs, vs );
        if(verbosity>2)printf( "fit1D |F[%i]|=%g |Err|=%g |Freg|=%g cos(f,v)=%g\n",itr,sqrt(cfv.y), sqrt(err2sum), sqrt(f2reg), cfv.x/sqrt(cfv.y*cfv.z) );
        if(cfv.y<F2max){ break; };
    }
    //if(verbosity>1)
    printf( "Bspline::fit1D() iter=%i err=%g \n", itr, sqrt(err2sum) );
    delete [] ps;
    delete [] fs;
    delete [] vs;
    return itr;
}


__attribute__((hot)) 
int fit1D_old2( const int n, double* Gs,  double* Es, double* Ws, double Ftol, int nmaxiter=100, double dt=0.1, bool bHalf=false ){
    if(verbosity>1)printf("Bspline::fit1D_old2() !!!!!!\n");
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
    if(verbosity>1)printf( "Bspline::fit1D_old2() iter=%i err=%g \n", itr, sqrt(err2sum) );
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

} // namespace SplineBcub{

#endif



