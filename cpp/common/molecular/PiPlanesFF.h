
#ifndef PiPlanes_h
#define PiPlanes_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

class PiPlane{  public:
    Vec3d dir; // this would make problem with optimizer 
    Vec3d pos;
    double K;
    double R;
    Vec3d fpos;
    Vec3d fdir; 
    int  natom;
    int* atoms=0;
};

class PiPlanes{ public:
    int natoms;
    Vec3d* apos=0;
    Vec3d* aforce=0;
    std::vector<PiPlane> planes;

    double evalPlane(int ip){
        double E = 0; 
        PiPlanes& pl = planes[ip];
        for(int i=0;i<pl.natoms;i++){
            int ia = pl.atoms[i];
            Vec3d dp; dp.set_sub(  pl.pos, apos[ia] );
            double cp = pl.dir.dot( dp ); 
            E        += K*cp*cp;
            Vec3d f; f.set_mul( pl.dir,  K*cp );
            aforce[ia].add( f );
            pl.fpos   .sub( f );
            pl.fdir   .add_mul( dp, K*cp );
            // if radial damp 
            // double r2 = dp.
            // er = 1-r2;
            // E *= er*er;
        }
        return E;
    }

}; // MMFF

#endif
