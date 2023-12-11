

constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

int verbosity = 1;
int idebug    = 0;
//double tick2second=1e-9;

#include "testUtils.h"
#include "QuadratureCoefsOptimizer.h"

QuadratureCoefsOptimizer W;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
#include "libUtils.h"

extern "C"{

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}

void init_buffers(){
    ibuffers.insert( { "ndims",  &W.ntrain } );
    buffers.insert( { "train_cs",  W.train_cs } );
    buffers.insert( { "qps",     (double*)W.qps      } );
    buffers.insert( { "qws",     W.qws      } );
    buffers.insert( { "qys",     W.qys      } );
    buffers.insert( { "qps_ref", (double*)W.qps_ref  } );
    buffers.insert( { "qws_ref", W.qws_ref  } );
    buffers.insert( { "qys_ref", W.qys_ref  } );
    //ibuffers.insert( { "bond2atom",      (int*   )W.bond2atom      } );
}

void setParamsSize( int* cs ){
    W.nparP = { cs[0], cs[1], cs[2] };
    W.nparQ = { cs[3], cs[4], cs[5] };
}

void setQuadraturePoints( int nqps_, double* qps_, bool bRef, bool bAlloc, double* qws_, double* qys_ ){
    if(bRef){ W.setReferenceQuadraturePoints( nqps_, (Vec3d*)qps_, bAlloc, qws_, qys_ ); }
    else    { W.setQuadraturePoints         ( nqps_, (Vec3d*)qps_, bAlloc, qws_, qys_ ); }
};

double evaluateQuadrature( int n, double* ps, double* params, double* ws, double* ys ){
    return W.evaluateQuadrature( n, (Vec3d*)ps, params, ws, ys );
}

double distributPointsTetrahedron( int n, bool bRef, bool bAlloc, int imode ){
    switch (imode){
        case 0: return W.distributPointsTetrahedron_open ( n, bRef, bAlloc ); break;
        case 1: return W.distributPointsTetrahedron_close( n, bRef, bAlloc ); break;
        case 2: return W.distributPointsTetrahedron_Shunn( n, bRef, bAlloc ); break;
    }
    return 0;
}



} // extern "C"
