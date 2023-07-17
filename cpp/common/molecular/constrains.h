
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
    Vec3d  shift;   // like PBC_shift
    double flim;
    bool active;


    DistConstr()=default;
    DistConstr( Vec2i ias_, Vec2d ls_, Vec2d ks_, double flim_=1e+300, Vec3d shift_=Vec3dZero ):ias(ias_),ls(ls_),ks(ks_),flim(flim_),shift(shift_),active(true){ };

    inline double apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec =0 )const{
        Vec3d sh;
        if(lvec){ lvec->dot_to_T( shift, sh ); }else{ sh=shift; }
        Vec3d d   = ps[ias.b] -ps[ias.a] + shift;
        double l  = d.norm(); 
        double f,E;
        E = spring( l, ls, ks, flim, f );
        d.mul(f/l);
        fs[ias.b].sub(d);
        fs[ias.a].add(d);
        return E;
    }

    void print(){ printf( "ias(%i,%i) ls(%f,%f) ks(%lf,%lf) shift(%lf,%lf,%lf) flim=%lf \n",   ias.a,ias.b,  ls.a,ls.b,   ks.a,ks.b,    shift.a,shift.b,shift.c,    flim ); };

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

    void apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec ){  
        for( const DistConstr&  c : bonds  ){ c.apply(ps,fs, lvec ); }
        for( const AngleConstr& c : angles ){ c.apply(ps,fs); }
    }

    int loadBonds( const char* fname, int _0=1 ){
        FILE* pFile = fopen( fname, "r" );
        if(pFile==0){ printf("ERROR in Constrains::loadBonds(%s) - No Such File \n", fname ); return -1; }
        else{
            char buff[1024];            
            int i=0;
            for(i=0; i<10000; i++){
                char* line = fgets( buff, 1024, pFile );
                if(line==NULL)  break;
                if(line[0]=='#')continue;
                DistConstr cons; cons.active=true;
                int nret = sscanf( line, "%i %i   %lf %lf   %lf %lf   %lf %lf %lf  %lf ",   &cons.ias.a,&cons.ias.b,  &cons.ls.a,&cons.ls.b,   &cons.ks.a,&cons.ks.b,    &cons.shift.a,&cons.shift.b,&cons.shift.c,    &cons.flim   );
                cons.ias.a-=_0;
                cons.ias.b-=_0;
                if(nret<10){ printf("WARRNING : Constrains::loadBonds[%i] nret(%i)<10 line=%s", nret, line ); }
                cons.print();
                bonds.push_back( cons );
            }
            return i;
            fclose( pFile );
        }
    }
    
};

#endif
