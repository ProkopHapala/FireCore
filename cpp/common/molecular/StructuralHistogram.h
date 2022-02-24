
#ifndef StructuralHistogram_h
#define StructuralHistogram_h

#include<math.h>
#include <Vec3.h>
#include "VecN.h"
#include "CG.h"

// Structural Histrogram should be used form
//  1)  Analysis of molecular geometry
//  2)  Usefull statistics for Machine Learning
//  3)  structure-hasing to estimate if two structires are similar - invariant to translation, rotation, and permutation of atoms
//          if histograms differs considerably, it is clear that the structures differs considerably
//          if histograms are the same, it is still possible that the structures are different => needs more carful examination by exact comparison method (e.g. grid-accelerated collision) 


inline void acumHistLinEqui( double x, double xmin, double xmax, int n, double* hist){
    if( (x<xmin)||(x>xmax) ) return;
    double xi = n*(x-xmin)/(xmin-xmax);
    int    ix = (int)xi; xi-=xi;
    //int ix; xi=modf(xi,ix);
    hist[ix  ]+=1-xi;
    hist[ix+1]+=  xi;
}

inline void acumHistLinNoeq( double x, double xmin, double xmax, int n, double* hist, double* xs ){
    for(int i=0; i<n; i++){
        double xi=xs[i];
        if(x>xi){
            double f   =(x-xi)/(xs[i+1]-xi);
            hist[i ]+=1-f;
            hist[i+1]+=  f;
        }
    }
}

inline int typeComb(int i, int j, int n){
    if(i>j) _swap(i,j);
    return j + (i*(i+1))/2;
}

class StructuralHistogram_h{ public:
    
    int ntypes;
    int ntypecomb; // ntypes*(ntypes+1)/2
    int nhist;
    Vec2d*   bondHistBounds=0;
    double** bondHistType=0; // 
    double** bondSamps=0;    // position of bond samples (bins boundaries) 
    int      excludeType=0;  // exclude type - e.g. Capping Hydrogen ?
    // perhaps most reasonable is to take
    //   1) histogram of bonds between structural atoms ( e.g. C,O,N ... ) 
    //   2) histogram of bonds between hydrogen and structural atom
    //  NOTE : but we do not make histograms only for bonded atoms, we make histrograms also for distant atoms

    int evalHist_bond( int natoms, Vec3d* apos, int* atypes ){
        int itypecomb=0;
        for(int i=0; i<natoms; i++){
            int ityp=atypes[i];
            if(ityp==excludeType) continue;
            const Vec3d& pi=apos[i];
            int typeoff=(ityp*(ityp+1))/2; 
            for(int j=i; j<natoms; j++){
                int jtyp=atypes[j];
                if(ityp==excludeType) continue;
                int itypecomb = typeoff + jtyp;
                //int itypecomb = typeComb( ityp,jtyp );
                Vec3d  dp; dp.set_sub( apos[j], pi );
                double r2 = dp.norm2();
                const Vec2d& xlim =  bondHistBounds[itypecomb];
                if( (r2<(xlim.a*xlim.a) )||(r2>(xlim.b*xlim.b)) ) continue;
                double x = sqrt(x);  // we can ommit this if bondSamps[] are in square
                acumHistLinEqui( x, xlim.a, xlim.b, nhist, bondHistType[itypecomb] );
                itypecomb++;
            } 
        }
    }

}; // DirStiff

bool structureOverlap( double Rmax, int natoms, Vec3d* aposA, int* atypesA, Vec3d* aposB, int* atypesB ){
    double R2max=Rmax*Rmax;
    for(int i=0; i<natoms; i++){
        int          ityp=atypesA[i];
        const Vec3d&   pi=aposA  [i];
        double R2min = 1e+300;
        for(int j=0; j<natoms; j++){
            if( atypesB[j]==ityp ){
                Vec3d dp; dp.set_sub(aposA[j],pi);
                double r2 = dp.norm2();
                R2min = fmax( R2min, r2 );
            }
        }
        if( R2min>R2max ) return false;
    }
    return true;
}

#endif
