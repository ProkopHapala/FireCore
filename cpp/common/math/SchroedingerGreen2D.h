
#ifndef  SchroedingerGreen2D_h
#define  SchroedingerGreen2D_h

#include "macroUtils.h"
#include "VecN.h"
//#include "fastmath.h"

#include "Lingebra.h"

const double hbar_eVfs =  0.6582119569;  // [eV*fs]


//const double hbar = 6.582119569e-16;   // [eV*s]
//const double Me   = 9.1093837e-31;     // [kg]  


// hbar = 6.582119569e-16  [eV*s] = 6.582119569e-16 = 0.6582119569 eV*fs
// h    = 4.1356676966e-15 [eV*s] = 4.1356676966 eV*fs 


// E = h*f = h*(1/t)
// t = h/E  =>    eV = 4.135667696e-15 [s] = 4.135667696 [fs]    // fs = 1e-15 s

// eV = 1.602176634e-19 [J]
// eV = 1.602176634e-19 [kg*m^2/s^2] = 1.602176634e-19 * ( (1e-15)^2/((9.1093837e-31)*(1e-10)^2)   * [Me*A^2/fs^2]
// eV = 17.5882001106   ((Me*A^2)/(fs^2))   


// J  = kg*m^2/s^2


// (J^2 * s^2)/(kg * m^2 ) =  J * (kg*m^2/s^2) * (s^2) / (kg * m^2) = J * 1 

//   ( hbar**2)/(2*Me) =   0.6582119569^2    (eV^2 * fs^2) / ( 2*Me * A^2 ) =  0.6582119569^2 * 0.5 * (eV^2) * (fs^2/(Me*A^2)) =  0.6582119569^2 * 0.5 *    (eV^2)/(eV/17.5882001106)  =      3.80998211619   eV       

//   


/*

E0 ... target energy
E  ... current Energy
L  ... is laplace operator

minimize{ (E-E0)^2 } by yi
(d/dyi) * (E-E0)^2    =     2*(E-E0) * dE/dyi 

E = <Y|H|Y> / <Y|Y>
Q = <Y|Y>

dE/dyi = ( Q* (d<Y|H|Y>/dyi) - <Y|H|Y>*(dQ/dyi) )/ Q^2

H = (L+V)
d<Y|H|Y>/dyi =   d<Y|L+V|Y>/dyi 


L[ix  ,iy]*y[ix  ,iy] = y[ix  ,iy]* ( - 4*y[ix  ,iy]  + y[ix+1,iy] + y[ix-1,iy] + y[ix  ,iy+1] + y[ix  ,iy-1] )
L[ix+1,iy]*y[ix+1,iy] = y[ix+1,iy]* ( - 4*y[ix+1,iy]  + y[ix+2,iy] + y[ix  ,iy] + y[ix+1,iy+1] + y[ix+1,iy-1] )

dL[ix,iy]/dy[ix  ,iy] =  -8*y[ix,iy] +     (y[ix+1,iy] + y[ix-1,iy] + y[ix,iy+1] + y[ix,iy-1])
dL[ix,iy]/dy[ix+1,iy] =  +y[ix,iy]

dL[ix+1,iy  ]/dy[ix,iy] =  +y[ix+1,iy  ]
dL[ix-1,iy  ]/dy[ix,iy] =  +y[ix-1,iy  ]
dL[ix  ,iy+1]/dy[ix,iy] =  +y[ix  ,iy+1]
dL[ix  ,iy-1]/dy[ix,iy] =  +y[ix  ,iy-1]

dL/dy[ix,iy] = -8*y[ix,iy] +    2*(y[ix+1,iy] + y[ix-1,iy] + y[ix,iy+1] + y[ix,iy-1])

*/

class SchroedingerGreen2D : public LinSolver{ public:
    int iter=0;
    int nx=0,ny=0, ntot=0;
    double*  V  =0; // potential
    double*  psi=0; // psi
    double* fpsi=0; // derivatives of
    double* vpsi=0; // velocity of psi change
    double* source=0;
    double E=0,Q=0,F2sum=0,E0=0.0,dstep=1.,cL=1.0;

    bool   bNromForce=false;
    double KnormForce=0.0;

//void init(int nx_, int ny_, double* V_=0, double* source_=0, double* psi_=0, double* fpsi_=0){
double init(int nx_, int ny_, double dstep_=0.1, double m_Me=1.0 ){
    nx=nx_; ny=ny_; ntot=nx*ny;
    iter=0;
    double* source_ = b;
    //fpsi_   = b;
    //_bindOrRealloc( ntot,    V_,    V );
    //_bindOrRealloc( ntot,  psi_,  psi );
    //_bindOrRealloc( ntot, fpsi_, fpsi );
    //_bindOrRealloc( ntot, source_, source );
    _realloc( V,      ntot );
    _realloc( source, ntot );
    _realloc( psi,    ntot ); 
    _realloc( fpsi,   ntot );                           // ToDo: Needed just for gradient descent
    _realloc( vpsi,   ntot ); VecN::set(ntot,0.0,vpsi); // ToDo: Needed just for MDdamp
    LinSolver::setLinearProblem(ntot, psi, source );
    return setStep( dstep_, m_Me );
}

double setStep( double dstep_, double m_Me=1.0 ){
    // hbar = 6.582119569e-16  [eV*s] = 6.582119569e-16 = 0.6582119569 eV*fs
    // h    = 4.1356676966e-15 [eV*s] = 4.1356676966 eV*fs 
    // eV   = 17.5882001106    [(Me*A^2)/(fs^2)]  
    //   (hbar**2)/(2*Me) =   0.6582119569^2    (eV^2 * fs^2) / ( 2*Me * A^2 ) =  0.6582119569^2 * 0.5 * (eV^2) * (fs^2/(Me*A^2)) =  0.6582119569^2 * 0.5 *    (eV^2)/(eV/17.5882001106)  =      3.80998211619   eV    
    dstep = dstep_;
    //cL    = -3.80998211619/(m_Me*dstep*dstep);
    cL    = -3.8099821121133326/(m_Me*dstep*dstep);
    printf( "setStep(): dstep %g m_Me %g cL %g \n", dstep, m_Me, cL );
    return cL;
}

double sumE(){ // energy (hamiltonian)
    double E  =0;
    double Q_ =0;
    //double cL = -1/(dstep*dstep);
    double Lsum_DEBUG = 0;
    double Vsum_DEBUG = 0;
    for(int iy=0; iy<ny; iy++){
        int i0y = iy*nx;
        for(int ix=0; ix<nx; ix++){
            int i = i0y+ix;
            double yi = psi[i];
            double vi = V  [i];
            //---- Laplacian
            // https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
            double li = -4*yi;
            if((ix+1)<nx){ li+=psi[i+1 ]; }
            if((ix-1)>=0){ li+=psi[i-1 ]; }
            if((iy+1)<ny){ li+=psi[i+nx]; }
            if((iy-1)>=0){ li+=psi[i-nx]; }
            li*=cL;
            // ---- Total energy
            Lsum_DEBUG+=  li*yi;     // DEBUG
            Vsum_DEBUG+=  vi*yi*yi;  // DEBUG
            double ei  = yi *( vi*yi + li );
            E  += ei;
            Q_ += yi*yi;
        }
    }
    //printf( "sumE()  Q %g  E %g   E/Q  %g \n", Q_, E,  E/Q_, dstep );
    printf( "iter[%i] sumE()  Q %g  E %g   E/Q  %g Lsum %g Vsum %g  cL %g dstep %g \n", iter, Q_, E,  E/Q_, Lsum_DEBUG/Q_, Vsum_DEBUG/Q_, cL, dstep );
    E*=1/Q_;
    return E;
}

void sumF( double factor, bool bQderiv=true ){ // derivatives of energy
    //double cL = -1/(dstep*dstep);
    for(int iy=0; iy<ny; iy++){
        int i0y = iy*nx;
        for(int ix=0; ix<nx; ix++){
            int i = i0y+ix;
            double yi = psi[i];
            double vi = V  [i];
            //---- Laplacian
            // https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
            double li = -4*yi;
            if((ix+1)<nx){ li+=psi[i+1 ]; }
            if((ix-1)>=0){ li+=psi[i-1 ]; }
            if((iy+1)<ny){ li+=psi[i+nx]; }
            if((iy-1)>=0){ li+=psi[i-nx]; }
            li*=cL;
            // ---- Total energy
            double f = (li + yi*vi)*2;           // d<Y|L-V|Y>/dyi
            //if( (li*li)>1e-4 ) printf( "[%i,%i] dY %g Y %g \n", iy,ix, f, yi );
            /*
            if(bQderiv){
                double ei = yi *( vi*yi + li );  // <Y|L-V|Y>
                f         = f*Q - 2*yi*ei;       // dE/dyi = ( Q* (d<Y|H|Y>/dyi) - Y*(dQ/dyi) )/ Q^2
            }
            if(source){
               f += source[i]; // ToDo: should it be here ?
            }
            */
            fpsi[i] = f*factor;
            //fpsi[i] = 0;
        }
    }
}

virtual void dotFunc( int n, double * psi, double * Ax ) override {
    // Operator is b = A*y = (H-e0*I)*Y   =    ( L + V -e0*I )*Y

    //bNromForce=true;
    //KnormForce=0.01;
    double KQ =0;
    //if(bNromForce){
    //    double Q  = VecN::norm2(ntot,psi);
    //    printf("Q[%i] %g\n", istep, Q);
    //    KQ = (Q - 1)*-KnormForce;
    //}

    for(int iy=0; iy<ny; iy++){
        int i0y = iy*nx;
        for(int ix=0; ix<nx; ix++){
            int i = i0y+ix;
            double yi = psi[i];
            double vi = V  [i];
            //---- Laplacian
            // https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
            double li = -4*yi;
            if((ix+1)<nx){ li+=psi[i+1 ]; }
            if((ix-1)>=0){ li+=psi[i-1 ]; }
            if((iy+1)<ny){ li+=psi[i+nx]; }
            if((iy-1)>=0){ li+=psi[i-nx]; }
            li*=cL;
            // ---- Total energy
            //if(verbosity>2){
            //    printf( "[%i,%i] li %g vi*yi %g yi*E0 %g \n", ix,iy,  li, yi*vi, -yi*E0  );
            //}
            double f  = (li + yi*(vi - E0))*-1;
            f        += KQ*yi;
            Ax[i]     = f;      //  (H-e0*I)*Y   = ( L + V -e0*I )*Y
            //fpsi[i] = 0;
        }
    }
};

double moveGD(double dt){
    Q=0;
    F2sum=0;
    for(int i=0; i<ntot; i++){ 
        double y=psi[i];
        double f=fpsi[i];
        y     += f*dt;
        Q     += y*y;
        F2sum += f*f;
        psi[i] = y;
    }
    Q    /=ntot;
    F2sum/=ntot;
    printf( "iter[%i] Q %g F2sum %g \n", iter, Q, F2sum );
    return F2sum;
}

double moveMDdamp(double dt, double damping){
    double cD = 1- damping;
    Q=0;
    F2sum=0;
    for(int i=0; i<ntot; i++){ 
        double y=psi[i];
        double v=vpsi[i];
        double f=fpsi[i];
        v     *= cD;
        v     += f*dt;
        y     += v*dt;
        Q     += y*y;
        F2sum += f*f;
        vpsi[i] = v;
        psi [i] = y;
    }
    Q    /=ntot;
    F2sum/=ntot;
    printf( "iter[%i] Q %g F2sum %g \n", iter, Q, F2sum );
    return F2sum;
}

double step( double E0, double dt ){
    E = sumE();
    double dE = (E-E0)*-1;
    sumF( dE, true );
    //moveGD(dt);
    moveMDdamp(dt,0.1);
    VecN::mul(ntot, 1./sqrt(Q), psi,psi );
    iter++;
    return F2sum;
}


};

#endif

