
#ifndef Atoms_h
#define Atoms_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"

// ToDo later: make Atoms Child of  Points
class Points{ public:
    int     n=0;
    Vec3d * ps=0;
};

class Atoms{ public:
    int     natoms =0;
    int   * atypes =0;
    Vec3d * apos   =0;

    void realloc ( int n                             ){ natoms=n; _realloc(atypes,natoms); _realloc(apos,natoms); }
    void allocNew( int n                             ){ natoms=n;   _alloc(atypes,natoms);   _alloc(apos,natoms); }
    void dealloc (                                   ){             _dealloc(atypes);        _dealloc(apos);      }
    void bind    ( int n, int* atypes_, Vec3d* apos_ ){ natoms=n; atypes=atypes_; apos=apos_; }
    void copyOf(const Atoms& p){
        if(natoms!=p.natoms)realloc(p.natoms);
        memcpy( atypes, p.atypes, sizeof(int)  *natoms );
        memcpy( apos,   p.apos,   sizeof(Vec3d)*natoms );
    }

    Atoms(){};
    Atoms(int n){  realloc( n ); };
    Atoms(const Atoms& As, bool bCopy=true){ if(bCopy){ copyOf(As); }else{  bind(As.natoms,As.atypes,As.apos); } };


    void fromRigid( Vec3d* ps0, const Vec3d& p0, const Mat3d& rot ){ for(int i=0; i<natoms; i++){ rot.dot_to_T( ps0[i], apos[i] ); apos[i].add(p0);         } }
    void shift    ( Vec3d d                                       ){ for(int i=0; i<natoms; i++){ apos[i].add(d); } }

    void atomsToXYZ(FILE* fout){
        for(int i=0; i<natoms; i++){
            fprintf( fout, "%i %20.10f %20.10f %20.10f\n", atypes[i], apos[i].x,apos[i].y,apos[i].z );
        }
    }

};

#endif
