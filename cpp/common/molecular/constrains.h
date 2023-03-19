
#ifndef constrains_h
#define constrains_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"


struct DistConstr{

    Vec2i ias;
    Vec2d  ls;
    Vec2d  ks;
    double flim;

    DistConstr()=default;
    DistConstr( Vec2i ias_, Vec2d ls_, Vec2d ks_, double flim_=1e+300 ):ias(ias_),ls(ls_),ks(ks_),flim(flim_){};

    double apply( Vec3d* ps, Vec3d* fs )const{
        Vec3d d   = ps[ias.b] -ps[ias.a];
        double l  = d.norm(); 
        double f,E;
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
        d.mul(f/l);
        fs[ias.b].sub(d);
        fs[ias.a].add(d);
        return E;
    }

};

struct AngleConstr{

    Vec3i  ias;
    double l0;
    double kTens;
    double kPress;
    double fmax;

    double  apply( Vec3d* ps, Vec3d* fs )const{
        return 0;
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
