
#ifndef DirectionStiffness_h
#define DirectionStiffness_h

#include <Vec3.h>
//#include "VecN.h"

class DirStiff{ public:
// Rank Cartesian directions by Stiffness due to bond-lengh-constrains

    static constexpr int nneigh_max=16;
    int* neighs=0;
    int* nneighs=0;
    double* kneighs = 0;

    int  natoms=0, nbonds=0;
    Vec2i  * bond2atom = 0;
    //double * bond_0    = NULL;  // [A]
    double * bond_k    = 0;  // [eV/A] ?
    Vec3d  * apos      = 0;   // atomic position
    Vec3d  * bondDir   = 0;
    Vec3d  * dirStiff  = 0;

    bool   * amask = 0;

    void evalDirs(){
        // this function simply calcuate projected stiffness of cartesian direction of atomic movement due to bonds
        for(int ib=0; ib<nbonds; ib++){
            const Vec2i& b = bond2atom[ib];
            Vec3d vk; vk.set_mul( bondDir[ib], bond_k[ib] );
            vk.abs();
            dirStiff[b.a].add(vk);
            dirStiff[b.b].add(vk);
        }
    }

    void sparseLinearBondSmoother( Vec3d* dpos, Vec3d* dpos_, double trash ){
        // this shoud spread pull-vector between neighboring atoms depending on bond stiffness
        //  - in principle it is approximation of running force field relaxation, but without re-calculating stiffness matrix and bond directions
        for(int i=0; i<natoms; i++){ dpos_[i]=dpos[i]; }
        for(int ib=0; ib<nbonds; ib++){
            const Vec2i& b = bond2atom[ib];
            if( !(amask[b.a]||amask[b.b])  ) continue;
            Vec3d vk; vk.set_mul( bondDir[ib], bond_k[ib] );
            Vec3d fa; fa.set_div( vk, dirStiff[b.a] );
            Vec3d fb; fb.set_div( vk, dirStiff[b.b] );
            dpos_[b.a].add_mul();
            dpos_[b.b].add_mul();
        }
        return;
    }


    void sparseLinearBondSmoother( Vec3d* dpos, Vec3d* dpos_, double trash ){
        // this shoud spread pull-vector between neighboring atoms depending on bond stiffness
        //  - in principle it is approximation of running force field relaxation, but without re-calculating stiffness matrix and bond directions
        for(int ia=0; ia<natoms; ia++){
            for(int ing=0; ing<nneighs[ia]; ing++){
                int    ja=neighs[ing];
                double kb=kneighs[ing]; 

            /*
            const Vec2i& b = bond2atom[ib];
            if( !(amask[b.a]||amask[b.b])  ) continue;
            Vec3d vk; vk.set_mul( bondDir[ib], bond_k[ib] );
            Vec3d fa; fa.set_div( vk, dirStiff[b.a] );
            Vec3d fb; fb.set_div( vk, dirStiff[b.b] );
            dpos_[b.a].add_mul();
            dpos_[b.b].add_mul();
            */
        }
        return;
    }

    evalStiffnessSubSpace( int npick, int* iatoms  ){
        // calculate small stiffness matrix for several selected atoms, considering just bonds (no angles etc.)
        // Kij = E/(dxi*dxj)  
    }

}; // DirStiff

#endif
