
#ifndef MolecularGraph_h
#define MolecularGraph_h

#include <math.h>
#include "fastmath.h"
#include "Vec3.h"


class MolecularGraph{ public:
    int natom=0;
    int nbond=0;
    //Vec3d* apos=0;
    Vec2i* bond2atom=0;
    int*   ngNs=0;
    int*   ngIs=0;
    int*   atom2bond=0;
    int*   atom2neigh=0;
    int*   amask=0;
    int*   bmask=0;

    int    if0=0,if1=0,if2=0;
    int*   front_=0;
    //int*   front2=0;
    int    endAtom=0;
    bool   endFound=0;

    void bindOrRealloc(int natom_, int nbond_, Vec2i* bond2atom_=0 ){
        natom=natom_;
        nbond=nbond_;
        if(bond2atom_){ bond2atom=bond2atom_; }else{ realloc(bond2atom,nbond); }
        _realloc( ngNs,   natom );
        _realloc( ngIs,   natom );
        _realloc( amask,  natom );
        _realloc( bmask,  nbond );
        _realloc( front_, natom );
        //_realloc( front2, natom );
        _realloc( atom2bond, nbond*2 );
        _realloc( atom2neigh, nbond*2 );
    }

    void makeNeighbors(){
        // NOTE : this algorithm is very similar to atoms->cells
        // count number of neighbors
        for(int i=0; i<natom;i++){ ngNs[i]=0; };
        //printf( "nbond %i \n", nbond );
        for(int i=0; i<nbond;i++){
            const Vec2i& b=bond2atom[i];
            ngNs[ b.a ]++;
            ngNs[ b.b ]++;
        }
        // set initial position of neighbors
        int n=0;
        //ngIs[0]=n;
        for(int i=0; i<natom; i++){
            //printf( "atom[%i] N %i i0 %i \n", i, ngNs[i], n );
            ngIs[i]=n;
            n+=ngNs[i];
            ngNs[i]=0;
        }
        //printf( "nTOT : %i \n", n );
        // bond2atom -> atom2bond
        for(int i=0; i<nbond; i++){
            const Vec2i& b=bond2atom[i];
            int j;
            j=ngIs[b.a]+ngNs[b.a]; atom2bond[j]=i; atom2neigh[j]=b.b; ngNs[b.a]++;
            j=ngIs[b.b]+ngNs[b.b]; atom2bond[j]=i; atom2neigh[j]=b.a; ngNs[b.b]++;
        }
    }

    inline int getOtherAtom( int ib, int ia  ){
        const Vec2i& b=bond2atom[ib];
        return (b.a==ia) ? b.b : b.a;
    }
    void cleanMask(){ for(int i=0; i<natom; i++){ amask[i]=-1; }; }

    inline void addToFront(int ja, int color){
        amask[ja]=color;
        front_[if2]=ja;
        if2++;
        if(endAtom)endFound=true;
    }

    void addNeigh2front( int ia, int color, int bg){
        int  ng =ngNs[ia];
        int* ngs=atom2neigh+ngIs[ia];
        for(int i=0; i<ng; i++){
            int ja = ngs[i];
            if( amask[ja] == bg ){  addToFront(ja, color); }
        }
    }
    inline void clearFront(){ if0=0;if1=0;if2=0; }
    inline void swapFronts(){
        if0=if1;if1=if2;
        //_swap(front1,front2);
        //_swap(nfront1,nfront2);
    }

    void fillSubGraph( int ia, int color ){
        int if1=1;
        front_[if0]=ia;
        int bg = amask[ia];
        do{
            for(int i=if0; i<if1; i++){
                addNeigh2front( front_[i], color, bg);
            }
            swapFronts();
        }while(if1>if0);
    }

    void splitByBond(int ib, int color){
        const Vec2i& b=bond2atom[ib];
        amask[ b.b ]=color;
        endAtom=b.b;
        fillSubGraph( b.a, color );
    }

    void maskCaps( int color ){
        for(int i=0; i<natom; i++){
            if(ngNs[i]<2){
                amask[ i     ]=color;
                bmask[ngIs[i]]=color;
            }
        }
    }

    void selectBondNeighs(){

    }

    void printNeighs(){
        for(int i=0; i<natom; i++){
            int  ng =ngNs[i];
            int i0=ngIs[i];
            printf( "atom[%i] nneigh %i \n", i, ng );
            for(int j=0; j<ng; j++){
                printf( "   ia %i ib %i \n", atom2neigh[i0+j], atom2bond[i0+j] );
            }
        }
    }

};

#endif


