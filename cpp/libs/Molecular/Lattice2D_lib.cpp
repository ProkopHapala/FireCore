

//constexpr int ntmpstr=2048;
//char tmpstr[ntmpstr];

//int verbosity = 1;
//int idebug    = 0;
//double tick2second=1e-9;

#include "LatticeMatch2D.h"

LatticeMatch2D LM;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
//#include "libUtils.h"

extern "C"{

    int match( double* lat0, double* lat1, double Rmax, double dRmax, double dAngMax ){
        LM.lat0[0] = *   (Vec2d*)lat0;
        LM.lat0[1] = * (((Vec2d*)lat0) + 1);
        LM.lat1[0] = *   (Vec2d*)lat1;
        LM.lat1[1] = * (((Vec2d*)lat1) + 1);
        LM.walk2D( Rmax, dRmax );
        LM.angleToRange();
        LM.sort();
        //for(int i=0; i<LM.match_u.size(); i++){ printf( "match_u[%i] ang %g n,d(%i,%g) \n", i, LM.match_u[i].alpha, LM.match_u[i].n, LM.match_u[i].d ); }
        //for(int i=0; i<LM.match_v.size(); i++){ printf( "match_v[%i] ang %g n,d(%i,%g) \n", i, LM.match_v[i].alpha, LM.match_v[i].n, LM.match_v[i].d ); }
        return LM.matchAngles( dAngMax);
    }

    void getMatches( int* inds, double* errs, bool bSort, double* Ks ){ LM.exportMatch( (int4*)inds, (double4*)errs, 0, bSort, *(double4*)Ks ); }

} // extern "C"
