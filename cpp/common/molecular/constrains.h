
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
        //printf( "DistConstr:apply(%i,%i) l %g E %g f %g | ls(%g,%g) ks(%g,%g) flim %g lvec %li |sh| %g \n", ias.b, ias.a, l, E,f, ls.x,ls.y, ks.x,ks.y, flim, (long)lvec, sh.norm() );
        return E;
    }

    void print(){ printf( "bond_constr ias(%i,%i) ls(%f,%f) ks(%lf,%lf) shift(%lf,%lf,%lf) flim=%lf \n",   ias.a,ias.b,  ls.a,ls.b,   ks.a,ks.b,    shift.a,shift.b,shift.c,    flim ); };

};


struct AngleConstr{
    Vec3i  ias;     // (ia,ib,ic), ia is the center, ib and ic are the ends
    Vec2d  cs0;     // unitary complex number (cos,sin) of the equilibrium angle
    Vec3i  acell;   // indexes of PBC cell shift of bond a
    Vec3i  bcell;   // indexes of PBC cell shift of bond b
    double k;
    double flim;
    bool active;

    AngleConstr()=default;
    AngleConstr( Vec3i ias_, Vec2d cs0_, double k_ ):ias(ias_),cs0(cs0_),k(k_),active(true){ };

    inline double  apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec =0, Mat3d* dlvec =0 )const{
        Vec3d ashift,bshift;
        if(lvec){
            ashift=lvec->a*acell.a + lvec->b*acell.b + lvec->c*acell.c;
            bshift=lvec->a*bcell.a + lvec->b*bcell.b + lvec->c*bcell.c;
        }else{ ashift=Vec3dZero; bshift=Vec3dZero;}
        Vec3d f1,f2,h1,h2;
        h1.set_sub(ps[ias.x],ps[ias.y]+ashift); double ir1 = 1./h1.normalize();
        h2.set_sub(ps[ias.x],ps[ias.z]+bshift); double ir2 = 1./h2.normalize();
        double E =  evalAngleCosHalf( h1, h2, ir1,ir2, cs0, k, f1,f2 );
        //printf(  "AngleConstr::apply(%i,%i,%i) E= %g   cos(a)= %g \n", ias.x, ias.y, ias.z, E,  f1.dot(f2)/sqrt(f1.norm2()*f2.norm2()) );
        //fs[ias.x].sub(f1); fs[ias.y].add(f1);
        //fs[ias.x].sub(f2); fs[ias.z].add(f2);
        fs[ias.y].sub(f1); fs[ias.x].add(f1);
        fs[ias.z].sub(f2); fs[ias.x].add(f2);
        return E;
    }

    void print(){ printf( "angle_constr ias(%i,%i,%i) cs0(%f,%f) k(%lf) flim=%lf acell(%i,%i,%i) bcell(%i,%i,%i)\n",   ias.a,ias.b,ias.b,   cs0.x,cs0.y,      k,       flim, acell.a, acell.b, acell.c, bcell.a, bcell.b, bcell.c); };

};

class Constrains{ public:
    std::vector<DistConstr>  bonds;
    std::vector<AngleConstr> angles;

    double apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec=0, Mat3d* dlvec=0 ){
        double E=0;  
        int i=0;
        for( const DistConstr&  c : bonds  ){ E+= c.apply(ps,fs, lvec, dlvec ); }
        for( const AngleConstr& c : angles ){ E+= c.apply(ps,fs, lvec, dlvec ); }
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
                printf( "Constrains::loadBonds[i=%i](%s)\n", i, line );
                if(line==NULL)  break;
                if     (line[0]=='#'){continue;}
                else if(line[0]=='b'){

                    DistConstr cons; cons.active=true;
                    int nret = sscanf( line, "b %i %i   %lf %lf   %lf %lf   %lf    %lf %lf %lf",   &cons.ias.a,&cons.ias.b,  &cons.ls.a,&cons.ls.b,  &cons.ks.a,&cons.ks.b,  &cons.flim, &cons.shift.a,&cons.shift.b,&cons.shift.c );
                    cons.ias.a-=_0;
                    cons.ias.b-=_0;
                    if( atom_permut ){
                        printf( "permut bond (%i->%i)-(%i->%i) \n", cons.ias.a, atom_permut[cons.ias.a], cons.ias.b, atom_permut[cons.ias.b] ); 
                        cons.ias.a=atom_permut[cons.ias.a];
                        cons.ias.b=atom_permut[cons.ias.b];
                    }
                    if(nret<10){cons.shift=Vec3dZero;}
                    if(nret<7 ){ printf("WARRNING : Constrains::loadBonds[%i] bond nret(%i)<10 line=%s", i, nret, line ); }
                    cons.print();
                    bonds.push_back( cons );

                }else if(line[0]=='g'){
            
                    double ang;
                    AngleConstr cons; cons.active=true;
                    int nret = sscanf( line, "g %i %i %i   %lf  %lf %lf   %i %i %i  %i %i %i",    &cons.ias.a,&cons.ias.b,&cons.ias.c,   &ang,   &cons.k,   &cons.flim,   &cons.acell.a,&cons.acell.b,&cons.acell.c,   &cons.bcell.a,&cons.bcell.b,&cons.bcell.c );
                    cons.cs0.fromAngle( (ang/180.0)*M_PI*0.5 );
                    cons.ias.a-=_0;
                    cons.ias.b-=_0;
                    cons.ias.c-=_0;
                    if( atom_permut ){
                        printf( "permut angle (%i->%i)-(%i->%i)-(%i->%i) \n", cons.ias.b, atom_permut[cons.ias.b],   cons.ias.a, atom_permut[cons.ias.a],    cons.ias.c, atom_permut[cons.ias.c] ); 
                        cons.ias.a=atom_permut[cons.ias.a];
                        cons.ias.b=atom_permut[cons.ias.b];
                        cons.ias.c=atom_permut[cons.ias.c];
                    }
                    if(nret<12){cons.bcell=Vec3iZero;}
                    if(nret<9 ){cons.acell=Vec3iZero;}
                    if(nret<6){ printf("WARRNING : Constrains::loadBonds[%i] angle nret(%i)<6 line=%s", i, nret, line ); }
                    cons  .print();
                    angles.push_back( cons );
                    
                }
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
