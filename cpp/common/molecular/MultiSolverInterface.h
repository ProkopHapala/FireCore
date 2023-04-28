
#ifndef MultiSolverInterface_h
#define MultiSolverInterface_h

#include "Vec3.h"
#include "Mat3.h"
//#include "quaternion.h"

class SolverInterface{ public:
    virtual double solve  ( int nmax, double Tol )=0;
    virtual double getGeom( Vec3d* ps, Mat3d *lvec ) =0;    // bPrepared=true is used when whole population is downloaded before by downloadPop()
    virtual void   setGeom( Vec3d* ps, Mat3d *lvec ) =0;    // bPrepared=true is used when whole population is uploaded   before by uploadPop()
};

class MultiSolverInterface{ public:
    virtual int    paralel_size ( ) =0;                            // number of items which can run un paralel
    virtual double solve_multi   ( int nmax, double Tol ) =0;
    virtual double getGeom      ( int isys, Vec3d* ps, Mat3d *lvec, bool bPrepared ) =0;    // bPrepared=true is used when whole population is downloaded before by downloadPop()
    virtual void   setGeom      ( int isys, Vec3d* ps, Mat3d *lvec, bool bPrepared ) =0;    // bPrepared=true is used when whole population is uploaded   before by uploadPop()
    virtual void   downloadPop  () =0;
    virtual void   uploadPop    () =0;
};

#endif
