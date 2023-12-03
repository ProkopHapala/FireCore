

#ifndef  QuadratureCoefsOptimizer_h
#define  QuadratureCoefsOptimizer_h

// Use quadpy as a reference https://github.com/sigma-py/quadpy


#include "Vec3.h"
#include "DynamicOpt.h"



/**
 * @class QuadratureCoefsOptimizer
 * @brief A class for optimizing quadrature coefficients.
 * 
 * Algorithm:
 * -# we choose distribution of quadrature points over certain 3d domain ( e.g. tetrahedrons, pyramids, prism, quadrilaterals )
 * -# we generate large number of test functions ( e.g. polynomials P(x), or rational functions P(x)/Q(y) )
 * -# we calculate reference values of the integrals of test functions over the domain using high precision quadrature with large number of points, and store them in a vector
 * -# we run the optimization which adjusts the quadrature coefficients to minimize the error of the integrals of test functions over the domain
 */
class QuadratureCoefsOptimizer : public DynamicOpt{ public:


static double polynominal3d( Vec3d p, Vec3i np, double * params ){
    double P = 0;
    int ip = 0;
    double zn = 1;
    //#pragma omp parallel for simd collapse(2) reduction(+:P)
    for(int iz=0; iz<np.z; iz++){
        double yn = 1;
        for(int iy=0; iy<np.y; iy++){
            double xn = 1;
            for(int ix=0; ix<np.x; ix++){
                double w = params[ip];
                P += w*xn*yn*zn;
                ip++;
                xn*=p.y;
            }
            yn*=p.y;   
        }
        zn*=p.y;    
    }
    return P;
}



    // --- from DynamicOpt.h ---
    //ForceFunction getForce = 0;

    // --- Sizes ---
    int nqps     = 0; // number of quadrature points
    int nqps_ref = 0; // number of reference quadrature points
    int ntrain   = 0; // number of training functions
    int nparams  = 0; // total number of parameters of the training function
    int nparPtot = 0; // total number of parameters for the numerator if the training function
    int nparQtot = 0; // total number of parameters for the denominator if the training function
    Vec3i nparP; // (nx,ny,nz) number of parameters along each axix for numerator of the training function
    Vec3i nparQ; // (nx,ny,nz) number of parameters along each axix for denominator of the training function

    // --- training functions ---
    bool  bRationalTrainingFunction = true;
    double* train_I  = 0; // [ntrain] integral of training function
    double* train_cs = 0; // [ntrain*nparams] derivative of training function

    // --- quadrature points ---
    Vec3d*  qps  = 0; // [nqps*3]  quadrature points
    double* qys  = 0; // [nqps]   function values at quadrature points
    double* qws  = 0; // [nqps]   weights of quadrature points

    // --- reference quadrature points ---
    Vec3d*  qps_ref = 0; // [nqps_ref*3]  quadrature points
    double* qys_ref = 0; // [nqps]   function values at quadrature points
    double* qws_ref = 0; // [nqps]   weights of quadrature points

    // =============== functions
    
    void reallocQuadOpt( int nqps_, int ntrain_, Vec3i nparP_, Vec3i nparQ_ ){
        ntrain = ntrain_;
        nqps   = nqps_;
        nparP  = nparP_;
        nparQ  = nparQ_;
        nparPtot = nparP.x*nparP.y*nparP.z;
        nparQtot = nparQ.x*nparQ.y*nparQ.z;
        nparams  = nparPtot;
        _realloc(train_cs, ntrain*nparams);
        _realloc(train_I,  ntrain);
        _realloc(qps,       nqps);
    }





    /**
     * Calculates the value of a quadrature using the given quadrature points and parameters of the test function.
     * 
     * @param n The number of points in the quadrature.
     * @param ps An array of quadrature points.
     * @param params An array of parameters of the test function.
     * @return The calculated value of the quadrature.
     */
    double evaluateQuadrature( int n, Vec3d* ps, double * params, double* ws=0, double* ys=0 ){
        double I = 0;
        for(int i=0; i<n; i++){ // loop over quadrature points
            double y = polynominal3d( ps[i], nparP, params );
            if(bRationalTrainingFunction){
                y/= polynominal3d( ps[i], nparQ, params+nparPtot );
            }
            if(ws) y*=ws[i] = y;   // multiply by quadrature weights
            if(ys) ys[i]    = y;   // store function values at quadrature points
            I              += y;   // accumulate integral 
        }
        return I;
    }



    /**
     * Evaluates the reference values for the given number of points
     */
    void evaluateReferenceValues(){
        for(int i=0; i<ntrain; i++){
            train_I[i] = evaluateQuadrature( nqps_ref, qps_ref, train_cs+i*nparams, qys_ref );
        }
    }

    double evaluateVariationalDerivs( int n, double * ws, double * dws ){
        double totalError = 0;
        for(int i=0; i<ntrain; i++){
            double I = evaluateQuadrature( nqps, qps, train_cs+i*nparams, qys );
            double err = train_I[i] - I;  // error of the integral
            // I = sum_i( w_i * y_i )
            //   dI/dw_i = y_i
            // E = (I-I0)^2  = ( sum_i(w_i*y_i) - I0 )^2 
            //  dE/dw_i = 2*( sum_i(w_i*y_i) - I0 )*y_i = 2*(I-I0)*y_i
            //if(dws) dws[i] = y; // derivative of the integral over the quadrature weights
            for(int j=0; j<nqps; j++){
                dws[j] += err*qys[j];
            }
            totalError += err*err;
        }
        return totalError;
    }

    void allocQuadraturePoints( int nqps_, bool bRef ){
        if(bRef){ _realloc(qps_ref, nqps_); } else { _realloc(qps,     nqps_); }
        if(bRef){ _realloc(qws_ref, nqps_); } else { _realloc(qws,     nqps_); }
        if(bRef){ _realloc(qys_ref, nqps_); } else { _realloc(qys,     nqps_); }
    }

    void setReferenceQuadraturePoints( int nqps_, Vec3d* qps_, bool bAlloc=0, double* qws_=0, double* qys_=0 ){
        nqps_ref = nqps_;
        qps_ref  = qps_;
        if(qws_==0) { if(bAlloc)_realloc(qws_ref,nqps_ref); for(int i=0; i<nqps_ref; i++){ qws_ref[i]=1; } }
        else        {qws_ref = qws_; }
        if(qys_==0) { if(bAlloc)_realloc(qys_ref,nqps_ref); }
        else        {qys_ref = qys_; }
    }

    void setQuadraturePoints( int nqps_, Vec3d* qps_,  bool bAlloc=0,  double* qws_=0, double* qys_=0 ){
        nqps = nqps_;
        qps  = qps_;
        if(qws_==0) { if(bAlloc)_realloc(qws,nqps); for(int i=0; i<nqps; i++){ qws_ref[i]=1; } }
        else        {qws_ref = qws_; }
        if(qys_==0) { if(bAlloc)_realloc(qys,nqps); }
        else        {qys_ref = qys_; }
    }
    
    double distributPointsTetrahedron_open( int n, bool bRef, bool bAlloc ){
        double step = 1.0/n;
        int np = (n+1)*(n+2)*(n+3)/6;
        if(bAlloc) allocQuadraturePoints( np, bRef );
        double* ws = bRef ? qws_ref : qws;
        Vec3d*  ps = bRef ? qps_ref : qps;
        int ip=0;
        double wsum = 0;
        double w0 = 1/np;
        for(int ix=0; ix<n; ix++){
            Vec3d p;
            p.x = ix*step + step*0.5;
            for(int iy=0; iy<n-ix; iy++){
                p.y = iy*step + step*0.5;
                for(int iz=0; iz<n-ix-iy; iz++){
                    p.z = iz*step + step*0.5;
                    ps[ip] = p;
                    ws[ip] = w0;
                    ip++;
                }
            }
        }
        if(fabs(wsum-1)>1e-8){ printf("ERROR in QuadratureCoefsOptimizer::distributPointsTetrahedron_open wsum=%g\n", wsum ); exit(0); }
        return wsum;
    }

};

#endif

