

#include "globals.h"

#include "testUtils.h"
#include "LatticeMatch2D.h"

LatticeMatch2D LM;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
//#include "libUtils.h"

extern "C"{

    void walk2D( double* lat0, double* lat1, double Rmax, double dRmax ){
        LM.lat0[0] = *   (Vec2d*)lat0;
        LM.lat0[1] = * (((Vec2d*)lat0) + 1);
        LM.lat1[0] = *   (Vec2d*)lat1;
        LM.lat1[1] = * (((Vec2d*)lat1) + 1);
        LM.walk2D( Rmax, dRmax );           
    }

    int match( double* lat0, double* lat1, double Rmax, double dRmax, double dAngMax ){
        LM.lat0[0] = *   (Vec2d*)lat0;
        LM.lat0[1] = * (((Vec2d*)lat0) + 1);
        LM.lat1[0] = *   (Vec2d*)lat1;
        LM.lat1[1] = * (((Vec2d*)lat1) + 1);
        LM.walk2D( Rmax, dRmax );           
        LM.angleToRange();
        LM.sort();
        return LM.matchAngles( dAngMax);
    }

    int getVecMatch( int nmax, int* inds, int* ns, double* errs, bool bSort, bool bV, bool bCleans, bool bPrint ){
        int n;
        printf( "getVecMatch nu %i nv %i \n", LM.match_u.size(), LM.match_v.size() );
        DEBUG
        if(bV){ if(bCleans)LM.cleansRedudat( LM.match_v ); LM.printVecs( LM.match_v, 'v' );  n=LM.exportVecMatch( LM.match_v, nmax, (int2*)inds, ns, errs, 0, bSort, bPrint ); }
        else  { if(bCleans)LM.cleansRedudat( LM.match_u ); LM.printVecs( LM.match_u, 'u' );  n=LM.exportVecMatch( LM.match_u, nmax, (int2*)inds, ns, errs, 0, bSort, bPrint );  }
        return n;
    }

    void getMatches( int* inds, int* ns, double* errs, bool bSort, double* Ks ){ LM.exportMatch( (int4*)inds, (int2*)ns, (double4*)errs, 0, bSort, *(double4*)Ks ); }

} // extern "C"
