/*
Non-Bonded Force-Field
 should be easily plugged with any molecular dynamics, by either sharing pointer to same data buffer, or by copying data
*/

#ifndef NBFF_old_h
#define NBFF_old_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

#include "Forces.h"
#include "NBFF.h"

bool checkPairsSorted( int n, Vec2i* pairs ){
    int ia=-1,ja=-1;
    for(int i=0;i<n; i++){
        const Vec2i& b = pairs[i];
        //printf( "pair[%i] %i,%i | %i %i  | %i %i %i ", i, b.i, b.j,   ia,ja ,   b.i>=b.j,  b.i<ia, b.j<=ja );
        if(b.i>=b.j){ return false; }
        if(b.i<ia)  { return false; }
        else if (b.i>ia){ia=b.i; ja=-1; };
        if(b.j<=ja){ return false; }
        ja=b.j;
    }
    return true;
}

// =========== SR repulsion functions

//   good SR repulsion function should have spike in the center ... this is the case of Lorenz
//   (R2-r2)^2
//   (R2-r2)*((R2/r2)-1)^2
//   ((R2/r2)-1)^2
//   (((R2+w2)/(r2+w2))-1)^2


class NBFF_old : public NBFF{ public:
// non bonded forcefield
    //int n        = 0;
    //Quat4d* REQs  = 0;
    //Vec3d*  apos  = 0;
    //Vec3d*  fapos = 0;

    int nmask   = 0;

    Vec2i* pairMask = 0; // should be ordered ?

    double sr_w      = 0.1;
    double sr_dRmax  = 1.0; // maximum allowed atom movement before neighborlist update
    double sr_K      = 1.0; // collision stiffness
    double sr_Rscale = 0.5; // rescale atoms for short range
    Vec3d* psBack = 0;   // backup positon for neighbor list
    std::vector<Vec2i> neighList; // possibly too slow ?

void realloc(int n_, int nmask_){
    natoms=n_;
    nmask=nmask_;
    _realloc(REQs ,natoms);
    _realloc(apos ,natoms);
    _realloc(fapos,natoms);
    _realloc(pairMask,nmask);
}

void bindOrRealloc(int n_, int nmask_, Vec3d* apos_, Vec3d* fapos_, Quat4d* REQs_, Vec2i* pairMask_ ){
    natoms=n_;
    nmask=nmask_;
    _bindOrRealloc(natoms, apos_  ,apos  );
    _bindOrRealloc(natoms, fapos_ ,fapos  );
    _bindOrRealloc(natoms, REQs_  ,REQs);
    _bindOrRealloc(nmask,pairMask_,pairMask);
}

void setREQs(int i0,int i1, const Quat4d& REQ){
    for(int i=i0;i<i1;i++){ REQs[i]=REQ; }
}

void printAtomParams(){
    printf( " # NBFF::printAtomParams() \n" );
    for(int i=0;i<natoms;i++){
        printf(  "atom[%i] R %g E %g Q %g \n", i, REQs[i].x, REQs[i].y, REQs[i].z );
    }
}

void cleanForce(){
    for(int i=0; i<natoms; i++){ fapos[i].set(0.); }
}

double evalLJQ_sortedMask( const Vec3d& shift=Vec3dZero ){
    //printf(  "evalLJQ_sortedMask \n" );
    int im=0;
    const int N=natoms;
    double E=0;
    for(int i=0; i<N; i++){
        Vec3d fi = Vec3dZero;
        Vec3d pi = apos[i];
        pi.add( shift );
        const Quat4d& REQi = REQs[i];
        for(int j=i+1; j<N; j++){    // atom-atom
            // --- mask some atom pairs (e.g. those which are bonded), as long as atoms are sorted we can do it efficiently
            if( (im<nmask)&&(i==pairMask[im].i)&&(j==pairMask[im].j) ){
                im++; continue;
            }
            Vec3d fij = Vec3dZero;
            Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
            double E = addAtomicForceLJQ( apos[j]-pi, fij, REQij );

            /*
            if(E>100.0){
                //printf( "%i %i  %g \n", i, j, E );
                opengl1renderer.color3f(0.0,1.0,0.0); Draw3D::drawLine( apos[j],pi);
                //printf("%i %i %g \n", i, j, fij.norm());
                opengl1renderer.color3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fij*-1000,apos[j]);
                opengl1renderer.color3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fij*1000 ,apos[i]);
            }
            */

            fapos[j].sub(fij);
            fi   .add(fij);

        }
        fapos[i].add(fi);
    }
    //exit(0);
    return E;

}

double evalLJQ_pbc( Mat3d lvec, Vec3i npbc ){
    double E=0;
    for(int ix=-npbc.x;ix<=npbc.x;ix++){
        for(int iy=-npbc.y;iy<=npbc.y;iy++){
            for(int iz=-npbc.z;iz<=npbc.z;iz++){
                //E+=evalLJQ_shifted( lvec.a*ix + lvec.b*iy + lvec.c*iz );
                E+=evalLJQ_sortedMask( lvec.a*ix + lvec.b*iy + lvec.c*iz );
            }
        }
    }
    return E;
}

// ============= Short Range Force using neighbor list

bool checkListValid(){
    double R2=sq(sr_dRmax);
    const int N=natoms;
    for(int i=0; i<N; i++){
        Vec3d d; d.set_sub(apos[i],psBack[i]);
        if(d.norm2()>R2){ return false; }
    }
    return true;
}

void makeSRList(){
    //Vec2i* m=pairMask;
    int im=0;
    const int N=natoms;
    //printf( "N %i \n" );
    neighList.clear();
    for(int i=0; i<N; i++){
        const Vec3d& pi = apos[i];
        psBack[i]=pi;
        double Ri = REQs[i].x;
        const Quat4d& REQi = REQs[i];
        for(int j=i+1; j<N; j++){    // atom-atom
            if( (im<nmask)&&(i==pairMask[im].i)&&(j==pairMask[im].j) ){
                im++; continue;
            }
            Vec3d d   = apos[j]-pi;
            double r2 = d.norm2();
            double Rij = (Ri + REQs[j].x)*sr_Rscale + sr_dRmax;
            if( r2<Rij*Rij ){
                neighList.push_back({i,j});  //  TODO: is this fast enough?
            }
        }
    }
}

double evalSRlist(){
    double E=0;
    double w2 = sq(sr_w);
    for(const Vec2i& ij : neighList ){
        Vec3d fij = Vec3dZero;
        double Rij = (REQs[ij.i].x+REQs[ij.i].y)* sr_Rscale;
        //addForceR2   ( apos[ij.j]-apos[ij.i], fij, Rij*Rij, sr_K     );
        //addForceR2inv( apos[ij.j]-apos[ij.i], fij, Rij*Rij, sr_K, w2 );
        E += addForceR2inv( apos[ij.j]-apos[ij.i], fij, Rij*Rij, sr_K, w2 );
        fapos[ij.i].sub(fij);
        fapos[ij.j].add(fij);
    }
    return E;
}

};

#endif

