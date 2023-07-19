
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

    inline double apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec =0, Mat3d* dlvec =0 )const{
        Vec3d sh;
        if(lvec){ lvec->dot_to_T( shift, sh ); }else{ sh=shift; }
        Vec3d d   = ps[ias.b] -ps[ias.a] + sh;
        double l  = d.norm(); 
        double f,E;
        E = spring( l, ls, ks, flim, f );
        //f = ks.x*l; E= 0.5*ks.x*l*l;
        d.mul(f/l);
        fs[ias.b].sub(d);
        fs[ias.a].add(d);

        if( dlvec ){  // Forces on Lattice Vector
            //  dE/dsh = (dE/dr) * ( dr/dsh ) = f * ( dr/dsh ) = f * d|a-b+sh|/dsh = f/|a-b+sh| * sh = sh*(f/l) 
            sh.mul( f/l );
            dlvec->a.add_mul( sh, shift.a );
            dlvec->b.add_mul( sh, shift.b );
            dlvec->c.add_mul( sh, shift.c );
        }

        //fs[ias.b].add(d);
        //fs[ias.a].sub(d);
        //printf( "DistConstr:apply(%i,%i) l %g E %g f %g | ls(%g,%g) ks(%g,%g) lvec %li |sh| %g \n", ias.b, ias.a, l, E,f, ls.x,ls.y, ks.x,ks.y, (long)lvec, sh.norm() );
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

    double apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec=0, Mat3d* dlvec=0 ){
        double E=0;  
        for( const DistConstr&  c : bonds  ){ E+= c.apply(ps,fs, lvec, dlvec ); }
        for( const AngleConstr& c : angles ){ E+= c.apply(ps,fs); }
        return E;
    }

    int loadBonds( const char* fname, int* atom_permut=0, int _0=1 ){
        //printf("Constrains::loadBonds() \n" );
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
                if( atom_permut ){
                    printf( "permut bond (%i->%i)-(%i->%i) \n", cons.ias.a, atom_permut[cons.ias.a], cons.ias.b, atom_permut[cons.ias.b] ); 
                    cons.ias.a=atom_permut[cons.ias.a];
                    cons.ias.b=atom_permut[cons.ias.b];
                }
                if(nret<10){ printf("WARRNING : Constrains::loadBonds[%i] nret(%i)<10 line=%s", nret, line ); }
                cons.print();
                bonds.push_back( cons );
            }
            return i;
            fclose( pFile );
        }
    }

    void clear( bool bShring=false){
        bonds.clear();
        angles.clear();
        if(bShring){
            bonds.shrink_to_fit();
            angles.shrink_to_fit();
        }
    }
    
};

#endif
