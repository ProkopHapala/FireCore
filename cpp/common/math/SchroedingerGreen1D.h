
#ifndef  SchroedingerGreen1D_h
#define  SchroedingerGreen1D_h

#include "macroUtils.h"
#include "VecN.h"

const double hbar_eVfs =  0.6582119569;  // [eV*fs]

class SchroedingerGreen1D { public:
    int iter=0;
    int n=0;
    double*  V  =0; // potential
    double*  psi=0; // psi
    double* fpsi=0; // derivatives of
    double* vpsi=0; // velocity of psi change
    double* Apsi=0;
    double* source=0;
    double E=0,Q=0,F2sum=0,E0=0.0,dstep=1.,cL=1.0;

    bool   bNromForce=false;
    double KnormForce=0.0;

double init(int n_, double dstep_=0.1, double m_Me=1.0 ){
    n=n_;
    iter=0;
    _realloc( V,      n );
    _realloc( source, n );
    _realloc( psi,    n ); 
    _realloc( Apsi,   n ); 
    _realloc( fpsi,   n );                            // ToDo: Needed just for gradient descent
    _realloc( vpsi,   n ); VecN::set(n,0.0,vpsi);     // ToDo: Needed just for MDdamp
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

double applyGreen( double* Yin, double* Yout, double factor=1.0 ){ // derivatives of energy
    double Qo = 0;
    for(int i=0; i<n; i++){
        double yi  = psi[i];
        double vi  = V  [i];
        //---- Laplacian
        // https://en.wikipedia.org/wiki/Discrete_Laplace_operator#Finite_differences
        double li = -2*yi;
        if((i+1)<n ){ li+=psi[i+1 ]; }
        if((i-1)>=0){ li+=psi[i-1 ]; }
        li*=cL;
        double y_ = (li + (vi-E0)*yi)*factor;
        Qo += y_*y_;
        Yout[i] = y_;
    }
    return Qo;
}

double moveGD(double dt){
    Q=0;
    F2sum=0;
    for(int i=0; i<n; i++){ 
        double y=psi[i];
        double f=fpsi[i];
        y     += f*dt;
        Q     += y*y;
        F2sum += f*f;
        psi[i] = y;
    }
    Q    /=n;
    F2sum/=n;
    //printf( "iter[%i] Q %g F2sum %g \n", iter, Q, F2sum );
    return F2sum;
}

double moveMDdamp(double dt, double damping){
    double cD = 1- damping;
    Q=0;
    F2sum=0;
    for(int i=0; i<n; i++){ 
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
    Q    /=n;
    F2sum/=n;
    //printf( "iter[%i] Q %g F2sum %g \n", iter, Q, F2sum );
    return F2sum;
}

double stepResonance( double E0_, double dt ){
    /*
        The reason why this does not seem to work is that LLy and dy<Ly|Ly> have different prefactors 1/d^4 vs 1/d^3 therefore if we combine it with (V+e) it may need rescaling
    */
    E0 = E0_;
    double Q1 = applyGreen(  psi, Apsi,  1 );
    double Q2 = applyGreen( Apsi, fpsi, -1 );
    //moveMDdamp(dt,0.1);
    bool bSubstractNormChange=true;
    if(bSubstractNormChange){
        double c_ab = VecN::dot(n, psi, fpsi );
        double Qpsi = VecN::dot(n, psi,  psi );
        VecN::fma( n, fpsi, psi, -c_ab/sqrt(Qpsi), fpsi );
    }
    
    //moveGD(dt);
    moveMDdamp( dt, 0.1 );

    printf( "[%i] Q1 %g Q2 %g Q %g |F| %g \n", iter, Q1, Q2, Q, sqrt(F2sum) );
    E=Q1;
    VecN::mul(n, 1./sqrt(Q), psi,psi );
    iter++;
    return F2sum;
}


};

#endif

