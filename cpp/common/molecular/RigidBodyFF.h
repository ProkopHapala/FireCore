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
//#include "Vec3Utils.h"
//#include "Forces.h"
#include "NBFF.h"

class RigidBodyFF{ public:
    int n=0; //int naloc;
    // --- Rigid Body Pose
    Vec3d    *ps    =0;
    Vec3d    *frots =0;
    Vec3d    *vrots =0;
    Quat4d   *qrots =0;
    // --- Atomic Coords
    //int      *ns     =0;
    //Vec3d    **apos  =0;
    NBsystem *mols  =0;
    Vec3d    **apos0 =0;

    void projectAtoms(){
        Mat3d mrot;
        for(int i=0; i<n; i++){ 
            qrots[i].toMatrix(mrot);
            mols [i].fromRigid(apos0[i],ps[i], mrot); 
        }
    }

    void evalTorqs(){
        for(int i=0; i<n; i++){
            Vec3d tq=Vec3dZero; 
            mols [i].torq( ps[i], tq ); 
            //tq.add( torq( ns[i], const Vec3T<T>* ps, const Vec3T<T>* vs );
            frots[i]=tq;
        }
    }

    void rotateQuaternions( double dt ){
        for(int i=0; i<n; i++){ qrots[i].dRot_exact ( dt, vrots[i] ); };
    }

    //void realloc( int n_, Vec3d* ps_, Quat4d* qrots_, Vec3d* frots_, Vec3d* vrots_, Vec3d** apos0_=0, Vec3d** apos_=0 ){
    void realloc( int n_, Vec3d* ps_, Quat4d* qrots_, Vec3d* frots_, Vec3d* vrots_, Vec3d** apos0_=0, NBsystem* mols_=0 ){
        n=n_;
        _bindOrRealloc(n,ps_,ps    );
        _bindOrRealloc(n,frots_,frots );
        _bindOrRealloc(n,vrots_,vrots );
        _bindOrRealloc(n,qrots_,qrots );
        _bindOrRealloc(n,mols_ ,mols  );
        _bindOrRealloc(n,apos0_,apos0 ); // for(int i=0; i<i++)
    };

    void makePos0s(int na=-1){
        if(na<0){ na=0; for(int i=0; i<n; i++){ na+=mols[i].n; } }
        //printf( "RigidBodyFF::makePos0s() na=%i,n=%i \n", na, n );
        Vec3d* buff=new Vec3d[na];
        na=0;        
        for(int i=0; i<n; i++){
            Vec3d* ps0=buff+na;
            int ni=mols[i].n;
            Vec3d* ps=mols[i].ps;
            for(int j=0; j<ni; j++){ ps0[j]=ps[j]; }
            //for(int j=0; j<ni; j++){ ps0[j]=ps[j]; printf("# makePos0s[%i,%i] (%g,%g,%g) \n", i, j, ps0[j].x,ps0[j].y,ps0[j].z ); }
            apos0[i]=ps0;
            na+=ni;
        }
    }

};

#endif

