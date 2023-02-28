

//constexpr int ntmpstr=2048;
//char tmpstr[ntmpstr];

//int verbosity = 1;
//int idebug    = 0;
//double tick2second=1e-9;

#include "testUtils.h"
#include "FitREQ.h"

FitREQ fit;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
//#include "libUtils.h"

extern "C"{

void init_types(int ntyp, int* typeMask, double* typREQs ){
    fit.init_types( ntyp, (Vec3i*)typeMask, (Vec3d*)typREQs );
}

void setSystem( int isys, int na, int* types, double* ps, bool bCopy=false ){
    fit.setSystem( isys, na, types, (Vec3d*)ps, bCopy );
}

void setRigidSamples( int n, double* Es_, double* poses_, bool bCopy, bool bAlloc ){
    fit.setRigidSamples( n, Es_, (Mat3d*)poses_, bCopy );
}

double run( int nstep, double ErrMax, double dt, bool bRigid ){
    double Err=0;
    for(int i=0; i<nstep; i++){
        Err=0;
        fit.DOFsToTypes();
        if(bRigid){ Err= fit.evalDerivsRigid(); }
        else      { Err= fit.evalDerivsRigid(); }
        if( Err<ErrMax ){ break; }
        fit.move_GD(  dt );
    }
    return Err;
}

} // extern "C"
