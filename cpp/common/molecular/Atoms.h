
#ifndef Atoms_h
#define Atoms_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Mat3.h"

// ToDo later: make Atoms Child of  Points
class Points{ public:
    int     n=0;
    Vec3d * ps=0;
};

class Atoms{ public:
    int     natoms =0;
    int   * atypes =0;
    Vec3d * apos   =0;
    // --- for global optimization
    Mat3d * lvec   =0;  // ToDo: should this be pointer or full array ?
    double Energy  =0;
    long   id      =0;

    void realloc ( int n, bool bAtypes=true ){ natoms=n;  _realloc(apos,natoms); if(bAtypes)_realloc(atypes,natoms); }
    void allocNew( int n, bool bAtypes=true ){ natoms=n;  _alloc(apos,natoms);   if(bAtypes)_alloc(atypes,natoms);   }
    void dealloc (        bool bAtypes=true ){            _dealloc(apos);        if(bAtypes)_dealloc(atypes);        }
    void bind    ( int n, int* atypes_, Vec3d* apos_ ){ natoms=n; atypes=atypes_; apos=apos_; }

    //void bindOrRealloc(){}
    void copyOf(const Atoms& p){
        if(natoms!=p.natoms)realloc(p.natoms);
        if(lvec  !=p.lvec  ){ lvec=new Mat3d; *lvec=*(p.lvec); }
        memcpy( atypes, p.atypes, sizeof(int)  *natoms );
        memcpy( apos,   p.apos,   sizeof(Vec3d)*natoms );
        
    }

    Atoms() = default;
    Atoms(FILE* fin){ if(atomsFromXYZ(fin)<0){ natoms=-1; }; };
    Atoms(int n,bool bLvec=false, bool bAtypes=true ){ realloc( n, bAtypes ); if(bLvec){ lvec=new Mat3d; *lvec=Mat3dIdentity; };    };
    Atoms(const Atoms& As, bool bCopy=true){ if(bCopy){ copyOf(As); }else{  bind(As.natoms,As.atypes,As.apos); } };


    void fromRigid( Vec3d* ps0, const Vec3d& p0, const Mat3d& rot ){ for(int i=0; i<natoms; i++){ rot.dot_to_T( ps0[i], apos[i] ); apos[i].add(p0);         } }
    void shift    ( Vec3d d                                       ){ for(int i=0; i<natoms; i++){ apos[i].add(d); } }

    void print()const{ printf("Atoms::print() natom=%i\n", natoms); for(int i=0; i<natoms; i++){ printf( "[%i] atype %i apos(%6.3f,%6.3f,%6.3f)\n", i, atypes[i], apos[i].x, apos[i].y, apos[i].z ); } }


    bool cellFromString( char* s, Mat3d& lvec )const{
        char c[3]; Mat3d M;
        int n = sscanf( s, "%c%c%c %lf %lf %lf   %lf %lf %lf   %lf %lf %lf", c,c+1,c+2, &(M.a.x),&(M.a.y),&(M.a.z),   &(M.b.x),&(M.b.y),&(M.b.z),   &(M.c.x),&(M.c.y),&(M.c.z) );
        if( (n==12) && (c[0]=='l')&&(c[1]=='v')&&(c[2]=='s') ){ lvec=M; return true; }
        return false;
    }

    int atomsFromXYZ(FILE* file, bool bRealloc=true ){
        const int nbuf=1024;
        char buff[nbuf];
        char* line=0;
        int na=0;
        line = fgets(  buff, nbuf, file );  if(!line){ return -1; }
        int n = sscanf( line, "%i\n", &na );
        //printf( "atomsFromXYZ() na=%i \n", na );
        if( n != 1 ){ return -1; }
        if((na<=0)||(na>10000)){ printf( "ERROR in Atoms::atomsFromXYZ() natoms(%i) is invalid => Exit() \n", natoms ); return -1; };
        if(bRealloc)realloc( na );
        Mat3d M; char c[3];
        line = fgets(  buff, nbuf, file );
        n = sscanf( line, "%c%c%c %lf %lf %lf   %lf %lf %lf   %lf %lf %lf", c,c+1,c+2, &(M.a.x),&(M.a.y),&(M.a.z),   &(M.b.x),&(M.b.y),&(M.b.z),   &(M.c.x),&(M.c.y),&(M.c.z) );
        if( (n==12) && (c[0]=='l')&&(c[1]=='v')&&(c[2]=='s') ){  lvec=new Mat3d(); *lvec=M; }else{  printf("WARRNING: Atoms::atomsFromXYZ() lvec not read \n" ); }
        int typ; Vec3d p;
        //if(lvec)printf( "lvec %g %g %g  %g %g %g  %g %g %g  \n", lvec->a.x,lvec->a.y,lvec->a.z,  lvec->b.x,lvec->b.y,lvec->b.z,   lvec->c.x,lvec->c.y,lvec->c.z  );
        //printf( "atomsFromXYZ() natoms=%i \n", natoms );
        for(int i=0; i<natoms; i++){
            line = fgets(  buff, nbuf, file );
            //printf( "`%s`\n", line );
            n = sscanf( line, "%i %lf %lf %lf\n", &typ, &p.x, &p.y, &p.z );
            //printf( "atomsFromXYZ[%i] %i   %g %g %g   nread=%i\n", i, typ, p.x, p.y, p.z, n );
            atypes[i]=typ;
            apos  [i]=p;
        }
        return 0;
    }

    void atomsToXYZ(FILE* file, bool bN=false, bool bComment=false, Vec3i nPBC=Vec3i{1,1,1}, const char* comment="", bool bEnergy=true ){
        int npbc=nPBC.totprod();
        //printf( "atomsToXYZ() atypes=%li   natoms=%i npbc=%i natoms*npbc=%i \n", (long)atypes, natoms, npbc, natoms*npbc );
        if(bN      )fprintf( file, "%i\n", natoms*npbc );
        if(bComment){
           if(lvec){ fprintf( file, "lvs %g %g %g  %g %g %g  %g %g %g ", lvec->a.x,lvec->a.y,lvec->a.z,  lvec->b.x,lvec->b.y,lvec->b.z,   lvec->c.x,lvec->c.y,lvec->c.z ); }
           if(lvec){ fprintf( file, "E %g id %li ", Energy,id  ); }
           fprintf(file, "%s\n", comment );
        }
        Vec3d shift=Vec3dZero;
        for(int iz=0;iz<nPBC.z;iz++){for(int iy=0;iy<nPBC.y;iy++){for(int ix=0;ix<nPBC.x;ix++){  if(lvec)shift= lvec->c*iz + lvec->b*iy + lvec->a*ix;
            for(int i=0; i<natoms; i++){
                //printf( "atomsToXYZ[%i]\n", i );
                Vec3d p =  apos[i]+shift;
                fprintf( file, "%i %20.10f %20.10f %20.10f\n", atypes[i], p.x, p.y, p.z );
            }
        }}};
    }

    void toNewLattice( const Mat3d& lvec_new, Atoms* source=0 ){
        Vec3d* apos_  =apos;
        int*   atypes_=atypes;
        //printf("Atoms::toNewLattice() lvec\n");     printMat(*lvec);
        //printf("Atoms::toNewLattice() lvec_new\n"); printMat(lvec_new);
        if( source ){
            apos_   = source->apos;
            atypes_ = source->atypes;
            *lvec   = *(source->lvec);
        }
        Mat3d invLvec;
        lvec->invert_T_to( invLvec );
        for(int i=0;i<natoms;i++){
            Vec3d u; invLvec.dot_to( apos_[i], u );
            lvec_new.dot_to_T( u, apos[i] );
            atypes[i] = atypes_[i];
        }
        *(lvec) = lvec_new;
    }

};

#endif
