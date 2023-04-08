
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
    Mat3d * lvec   =0;
    double Energy  =0;
    long   id      =0;

    void realloc ( int n                             ){ natoms=n; _realloc(atypes,natoms); _realloc(apos,natoms); }
    void allocNew( int n                             ){ natoms=n;   _alloc(atypes,natoms);   _alloc(apos,natoms); }
    void dealloc (                                   ){             _dealloc(atypes);        _dealloc(apos);      }
    void bind    ( int n, int* atypes_, Vec3d* apos_ ){ natoms=n; atypes=atypes_; apos=apos_; }
    void copyOf(const Atoms& p){
        if(natoms!=p.natoms)realloc(p.natoms);
        memcpy( atypes, p.atypes, sizeof(int)  *natoms );
        memcpy( apos,   p.apos,   sizeof(Vec3d)*natoms );
    }

    Atoms() = default;
    Atoms(int n,bool bLvec=false          ){ realloc( n ); if(bLvec){ lvec=new Mat3d; *lvec=Mat3dIdentity; };    };
    Atoms(const Atoms& As, bool bCopy=true){ if(bCopy){ copyOf(As); }else{  bind(As.natoms,As.atypes,As.apos); } };


    void fromRigid( Vec3d* ps0, const Vec3d& p0, const Mat3d& rot ){ for(int i=0; i<natoms; i++){ rot.dot_to_T( ps0[i], apos[i] ); apos[i].add(p0);         } }
    void shift    ( Vec3d d                                       ){ for(int i=0; i<natoms; i++){ apos[i].add(d); } }

    void print()const{ printf("Atoms::print() natom=%i\n", natoms); for(int i=0; i<natoms; i++){ printf( "[%i] atype %i apos(%6.3f,%6.3f,%6.3f)\n", i, atypes[i], apos[i].x, apos[i].y, apos[i].z ); } }

    void atomsToXYZ(FILE* fout, bool bN, bool bComment){
        if(bN      )fprintf( fout, "%i\n", natoms );
        if(bComment){
           if(lvec){ fprintf( fout, "lvs %g %g %g  %g %g %g  %g %g %g  E %g id %li\n", lvec->a.x,lvec->a.y,lvec->a.z,  lvec->b.x,lvec->b.y,lvec->b.z,   lvec->c.x,lvec->c.y,lvec->c.z,  Energy,id  ); }
           else    { fprintf( fout, "E %g id %li\n", Energy,id  ); }
        }
        for(int i=0; i<natoms; i++){
            fprintf( fout, "%i %20.10f %20.10f %20.10f\n", atypes[i], apos[i].x,apos[i].y,apos[i].z );
        }
    }

    void toNewLattic( const Mat3d& lvec_new, Atoms* source=0 ){
        Vec3d* apos_  =apos;
        int*   atypes_=atypes;
        if( source ){
            apos_   = source->apos;
            atypes_ = source->atypes;
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
