/*
Non-Bonded Force-Field
 should be easily plugged with any molecular dynamics, by either sharing pointer to same data buffer, or by copying data
*/

#ifndef RigidBodyFF_h
#define RigidBodyFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"
#include "NBFF.h"

class RigidBodyFF{
    int n=0; int naloc;
    // --- Rigid Body Pose
    Vec3d    *ps    =0;
    Vec3d    *frots =0;
    Vec3d    *vrots =0;
    Quat4d   *qrots =0;
    // --- Atomic Coords
    NBsystem *mols  =0;
    Vec3d    **pos0s =0;

    void projectAtoms(){
        Mat3d mrot;
        for(int i=0; i<n; i++){ 
            qrots[i].toMatrix(mrot);
            mols [i].fromRigid(pos0s[i],ps[i], mrot); 
        }
    }

    void evalTorqs(){
        for(int i=0; i<n; i++){
            Vec3d tq=Vec3dZero; 
            mols [i].torq( ps[i], tq ); 
            frots[i]=tq;
        }
    }

    void rotateQuaternions( double dt ){
        for(int i=0; i<n; i++){ qrots[i].dRot_exact ( dt, vrots[i] ); };
    }

    void realloc( int n_, Vec3d* ps_, Quat4d* qrots_, Vec3d* frots_, Vec3d* vrots_, Vec3d** pos0s_=0, NBsystem* mols_=0 ){
        naloc=n_;
        _bindOrRealloc(n,ps_,ps    );
        _bindOrRealloc(n,frots_,frots );
        _bindOrRealloc(n,vrots_,vrots );
        _bindOrRealloc(n,qrots_,qrots );
        _bindOrRealloc(n,pos0s_,pos0s );
        _bindOrRealloc(n,mols_, mols  );
    };

};

#endif

