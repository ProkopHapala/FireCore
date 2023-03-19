
#ifndef constrains_h
#define constrains_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"
#include "Forces.h"


struct DistConstr{
    Vec2i ias;
    Vec2d  ls;
    Vec2d  ks;
    double flim;
    bool active;

    DistConstr()=default;
    DistConstr( Vec2i ias_, Vec2d ls_, Vec2d ks_, double flim_=1e+300 ):ias(ias_),ls(ls_),ks(ks_),flim(flim_),active(true){ };

    inline double apply( Vec3d* ps, Vec3d* fs )const{
        Vec3d d   = ps[ias.b] -ps[ias.a];
        double l  = d.norm(); 
        double f,E;
        E = spring( l, ls, ks, flim, f );
        /*
        if    (l>ls.x){
            double dl=l-ls.x;
            f=dl*ks.x;
            if(f>flim){
                f=flim;
                double dlim = flim/ks.x;
                E = (0.5*dlim*dlim)*ks.x + (dl-dlim)*flim;
            }else{
                E = 0.5*dl*dl*ks.x; 
            }
        }else if(l<ls.y){
            double dl=l-ls.y;
            f=dl*ks.y;
            if(f<-flim){
                f=-flim;
                double dlim = -flim/ks.y;
                E = (0.5*dlim*dlim)*ks.y - (dl-dlim)*flim;
            }else{
                E = 0.5*dl*dl*ks.y; 
            }
        };
        */
        d.mul(f/l);
        fs[ias.b].sub(d);
        fs[ias.a].add(d);
        return E;
    }

};

struct AngleConstr{
    Vec3i  ias;
    Vec2d  cs0;
    double k;
    bool active;

    AngleConstr()=default;
    AngleConstr( Vec3i ias_, Vec2d cs0_, double k_ ):ias(ias_),cs0(cs0_),k(k_),active(true){ };

    inline double  apply( Vec3d* ps, Vec3d* fs )const{
        Vec3d f1,f2,h1,h2;
        h1.set_sub(ps[ias.x],ps[ias.y]); double ir1 = 1./h1.normalize();
        h2.set_sub(ps[ias.x],ps[ias.z]); double ir2 = 1./h2.normalize();
        double E =  evalAngleCosHalf( h1, h2, ir1,ir2, cs0, k, f1,f2 );
        fs[ias.x].sub(f1); fs[ias.y].add(f1);
        fs[ias.x].sub(f2); fs[ias.z].add(f2);
        return E;
    }

};

class Constrains{ public:
    std::vector<DistConstr>  bonds;
    std::vector<AngleConstr> angles;

    void apply( Vec3d* ps, Vec3d* fs  ){  
        for( const DistConstr&  c : bonds  ){ c.apply(ps,fs); }
        for( const AngleConstr& c : angles ){ c.apply(ps,fs); }
    }
    
};

#endif
