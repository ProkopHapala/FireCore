

#ifndef  QuadratureCoefsOptimizer_h
#define  QuadratureCoefsOptimizer_h



// Use quadpy as a reference https://github.com/sigma-py/quadpy


/*
A simple way to generate quadrature points for terahedron is to start from FCC lattice.
The problem is that poits on the boundary of the tetrahedron should distribute its contribution to the neighboring solids.
The contribution to each solid is proportional to the solid angle occupied by the each solid.
Therefore in orde to caculate the weights we need to calculate the solid angles.

* Face: weight of points on face is 1/2  ( because the face splits to space to two half-spaces )
* Edge: for the edge we need to calculate the diheral angle between the two faces
* Vertex: for the vertex we need to calculate the solid angle from the three edges, this can be done using the triple product or polar sine fuction 

*/

// https://en.wikipedia.org/wiki/Triple_product
// https://en.wikipedia.org/wiki/Polar_sine



#include "Vec3.h"
#include "Mat3.h"
#include "DynamicOpt.h"

#include "integration.h"

/**
 * Calculates the value of a 3-dimensional polynomial at a given point.
 * 
 * @param p The point at which to evaluate the polynomial.
 * @param np The number of terms in each dimension of the polynomial.
 * @param params The coefficients of the polynomial.
 * @return The value of the polynomial at the given point.
 */
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

/**
 * Calculates the solid angle of a tetrahedron defined by three vectors.
 *
 * The solid angle is calculated using the formula described in the Wikipedia article:
 * https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
 *
 * @param a The first vector defining the tetrahedron.
 * @param b The second vector defining the tetrahedron.
 * @param c The third vector defining the tetrahedron.
 * @param bNormalize Flag indicating whether to normalize the vectors before calculation.
 *                   Default is false.
 * @return The solid angle of the tetrahedron in radians.
 */
double solidAngle_tetrahedron( const Vec3d& a, const Vec3d& b, const Vec3d& c, const bool bNormalize=false ){
    // https://en.wikipedia.org/wiki/Solid_angle#Tetrahedron
    double prod = a.triple_product( b, c );
    double denom;
    if(bNormalize){ 
        double ra,rb,rc;
        ra = a.norm(); rb = b.norm(); rc = c.norm();   
        denom = ra*rb*rc + a.dot(b)*rc + b.dot(c)*ra + c.dot(a)*rb;  
    }else{ 
        denom = 1 + a.dot(b) + b.dot(c) + c.dot(a);
    }
    prod/=denom;
    return 2*atan2( prod, denom );
};

double dihedral_angle( Vec3d edge, Vec3d a, Vec3d b, const bool bNormalize=true, const bool bNormalizeEdge=true, const bool bOrthogonalize=true ){
    // https://en.wikipedia.org/wiki/Dihedral_angle
    if(bNormalizeEdge) edge.normalize();
    if(bOrthogonalize){
        a.add_mul( edge, -edge.dot(a) );
        b.add_mul( edge, -edge.dot(b) );
    }
    double ca = a.dot(b);
    if(bNormalize){ ca/= a.norm()*b.norm(); }
    return acos( ca );
}

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

    // --- from DynamicOpt.h ---
    //ForceFunction getForce = 0;

    // --- Sizes ---
    int ntrain   = 0; // number of training functions
    int nqps     = 0; // number of quadrature points
    int nqps_ref = 0; // number of reference quadrature points
    int nparams  = 0; // total number of parameters of the training function
    int nparPtot = 0; // total number of parameters for the numerator if the training function
    int nparQtot = 0; // total number of parameters for the denominator if the training function
    Vec3i nparP{0,0,0}; // (nx,ny,nz) number of parameters along each axix for numerator of the training function
    Vec3i nparQ{0,0,0}; // (nx,ny,nz) number of parameters along each axix for denominator of the training function

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

    void allocQuadraturePoints( int n, bool bRef ){
        printf( "allocQuadraturePoints n=%i bRef=%i \n", n, bRef );
        if(bRef){ 
            nqps_ref = n;
            _realloc(qps_ref, n);
            _realloc(qws_ref, n);
            _realloc(qys_ref, n);
        } else { 
            nqps = n;
            _realloc(qps, n);
            _realloc(qws, n);
            _realloc(qys, n);
        }
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
        int np = (n)*(n+1)*(n+2)/6;
        if(bAlloc) allocQuadraturePoints( np, bRef );
        double* ws = bRef ? qws_ref : qws;
        Vec3d*  ps = bRef ? qps_ref : qps;
        int ip=0;
        double wsum = 0;
        double w0 = 1./np;
        for(int ix=0; ix<n; ix++){
            Vec3d p;
            p.x = ix*step + step*0.5;
            for(int iy=0; iy<n-ix; iy++){
                p.y = iy*step + step*0.5;
                for(int iz=0; iz<n-ix-iy; iz++){
                    p.z = iz*step + step*0.5;
                    printf( "distributPointsTetrahedron_open[%i] p(%i,%i,%i) (%10.8f,%10.8f,%10.8f) \n", ip, ix,iy,iz, p.x,p.y,p.z );
                    ps[ip] = p;
                    ws[ip] = w0;
                    wsum  += w0;
                    ip++;
                }
            }
        }
        if(   ip!=np        ){ printf("ERROR in QuadratureCoefsOptimizer::distributPointsTetrahedron_open ip(%i)!=np(%i)\n", ip, np ); exit(0); }
        if(fabs(wsum-1)>1e-8){ printf("ERROR in QuadratureCoefsOptimizer::distributPointsTetrahedron_open wsum=%g\n", wsum ); exit(0); }
        return wsum;
    }

    double distributPointsTetrahedron_Shunn( int ishun, bool bRef, bool bAlloc ){
        const Quat4d* ps=0;
        const double* ws=0;
        int np=0;
        switch (ishun){
            case 1: np=4;  ps = (const Quat4d*) TetrahedronIntegral::ps_4;  ws = TetrahedronIntegral::ws_4;  break;
            case 2: np=10; ps = (const Quat4d*) TetrahedronIntegral::ps_10; ws = TetrahedronIntegral::ws_10; break;
            case 3: np=20; ps = (const Quat4d*) TetrahedronIntegral::ps_20; ws = TetrahedronIntegral::ws_20; break;
            case 4: np=35; ps = (const Quat4d*) TetrahedronIntegral::ps_35; ws = TetrahedronIntegral::ws_35; break;
        }
        if(bAlloc) allocQuadraturePoints( np, bRef );
        double wsum = 0;
        if(bRef){ for(int i=0; i<np; i++){ qps_ref[i] = ps[i].f; qws_ref[i] = ws[i]; wsum+=ws[i]; }; }
        else    { for(int i=0; i<np; i++){ qps[i]     = ps[i].f; qws[i]     = ws[i]; wsum+=ws[i]; }; }
        return wsum;
    }


    double distributPointsTetrahedron_close( int n, bool bRef, bool bAlloc, Vec3d* vs=0 ){
        printf( "distributPointsTetrahedron_close n=%i bRef=%i bAlloc=%i \n", n, bRef, bAlloc );
        double vws[4];
        double ews[6];
        if(vs){ // calculate weights of vertices and edges
            double wS = 1./(4*M_PI); // angle of full sphere
            double wC = 1./(2*M_PI); // angle of full circle
            // weights of vertices
            vws[0] = solidAngle_tetrahedron( vs[1]-vs[0], vs[2]-vs[0], vs[3]-vs[0] )*wS;
            vws[1] = solidAngle_tetrahedron( vs[0]-vs[1], vs[2]-vs[1], vs[3]-vs[1] )*wS;
            vws[2] = solidAngle_tetrahedron( vs[0]-vs[2], vs[1]-vs[2], vs[3]-vs[2] )*wS;
            vws[3] = solidAngle_tetrahedron( vs[0]-vs[3], vs[1]-vs[3], vs[2]-vs[3] )*wS;
            // weights of edges
            ews[0] = dihedral_angle( vs[2]-vs[0], vs[1]-vs[0], vs[3]-vs[0] )*wC;   // w01 
            ews[1] = dihedral_angle( vs[1]-vs[0], vs[2]-vs[0], vs[3]-vs[0] )*wC;   // w02 
            ews[2] = dihedral_angle( vs[1]-vs[0], vs[3]-vs[0], vs[2]-vs[0] )*wC;   // w03 
            ews[3] = dihedral_angle( vs[3]-vs[1], vs[2]-vs[1], vs[0]-vs[1] )*wC;   // w12 
            ews[4] = dihedral_angle( vs[2]-vs[1], vs[3]-vs[1], vs[0]-vs[1] )*wC;   // w13
            ews[5] = dihedral_angle( vs[1]-vs[2], vs[3]-vs[2], vs[0]-vs[2] )*wC;   // w23
        }else{
            double wv = 0.043869914;  // acos(23./27)      / ( 4*M_PI )
            double we = 0.195913276;  // 2*atan(sqrt(2)/2) / ( 2*M_PI )
            for(int i=0; i<4; i++){ vws[i]=wv; }
            for(int i=0; i<6; i++){ ews[i]=we; }
        }
        double step = 1.0/n;
        int np = (n+1)*(n+2)*(n+3)/6;
        if(bAlloc) allocQuadraturePoints( np, bRef );
        double* ws = bRef ? qws_ref : qws;
        Vec3d*  ps = bRef ? qps_ref : qps;
        { // DEBUG
            for(int i=0; i<np; i++){  ws[i]=0; ps[i].set(NAN); } 
        }
        int ip=0;
        double wsum = 0;
        double w0 = 1./np;
        // --- vertices ---
        for(int i=0; i<4; i++){
            double w = vws[i]*w0;
            ps[ip].set(0.); if(i<3) ps[ip].array[i]=1;
            ws[ip] = w;
            wsum  += w;
            ip++;
        }
        // --- faces ---
        for(int i=1; i<n; i++){
            double ui = i*step;
            for(int j=1; j<n-i; j++){
                double uj = j*step;
                ps[ip].set(ui,uj,1-ui-uj); ws[ip]=w0*0.5; wsum+=ws[ip]; ip++;
                ps[ip].set(ui,uj,0); ws[ip]=w0*0.5; wsum+=ws[ip]; ip++;
                ps[ip].set(ui,0,uj); ws[ip]=w0*0.5; wsum+=ws[ip]; ip++;
                ps[ip].set(0,ui,uj); ws[ip]=w0*0.5; wsum+=ws[ip]; ip++;
            }
        }
        // --- edges ---  // ToDo: currently it is not correct
        for(int i=1; i<n; i++){
            double ui = i*step;
            ps[ip].set(0,ui,1-ui); ws[ip]=w0*ews[0]; wsum+=ws[ip]; ip++;
            ps[ip].set(ui,0,1-ui); ws[ip]=w0*ews[1]; wsum+=ws[ip]; ip++;
            ps[ip].set(ui,1-ui,0); ws[ip]=w0*ews[2]; wsum+=ws[ip]; ip++;
            ps[ip].set(ui,0,0); ws[ip]=w0*ews[3]; wsum+=ws[ip]; ip++;
            ps[ip].set(0,ui,0); ws[ip]=w0*ews[4]; wsum+=ws[ip]; ip++;
            ps[ip].set(0,0,ui); ws[ip]=w0*ews[5]; wsum+=ws[ip]; ip++;
        }
        // --- interior points ---
        for(int ix=1; ix<n; ix++){
            Vec3d p;
            p.x = ix*step;
            for(int iy=1; iy<n-ix; iy++){
                p.y = iy*step;
                for(int iz=1; iz<n-ix-iy; iz++){
                    p.z = iz*step;
                    ps[ip] = p;
                    double w = w0; 
                    ws[ip] = w;
                    wsum  += w;
                    ip++;
                }
            }
        }
        //if(   ip!=np        ){ printf("ERROR in QuadratureCoefsOptimizer::distributPointsTetrahedron_open ip(%i)!=np(%i)\n", ip, np ); exit(0); }
        //if(fabs(wsum-1)>1e-8){ printf("ERROR in QuadratureCoefsOptimizer::distributPointsTetrahedron_open wsum=%g\n", wsum ); exit(0); }
        return wsum;
    }




};

#endif

