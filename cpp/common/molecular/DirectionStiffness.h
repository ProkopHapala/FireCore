
#ifndef DirectionStiffness_h
#define DirectionStiffness_h

#include "macroUtils.h"
#include <Vec3.h>
#include "VecN.h"
#include "CG.h"

#include "MolecularGraph.h"

class Deformer{ public:
    int natoms;
    Vec3d* apos=0;
    Vec3d* aforce=0;
    DotFunc evalForce=0;
    MolecularGraph* graph=0;

    int     npick=0; // number of picked atoms;
    int*    picks=0; // picked atoms
    double* Ks=0;
    Vec3d*  pulls=0; // pull vectors



    int    nstep=1;   // number of relaxation steps
    double L    =1.0; // pull lenght
    double dt   =0.02; // relaxation speed

    void bind( int natoms_, Vec3d* apos_, Vec3d* aforce_ ){ natoms=natoms_; apos=apos_; aforce=aforce_;  }

    void initPicks(int npick_){
        npick=npick_;
        _realloc(picks,npick);
        _realloc(pulls,npick);
        _realloc(Ks,npick);
    }

    void genPicks(){
        for(int i=0; i<npick; i++){ picks[i]=rand()%npick; }
    }
    void genKs(double K){
        for(int i=0; i<npick; i++){ Ks[i]=randf(-K,K); }
    }

    void genPulls(){
        for(int i=0; i<npick; i++){ pulls[i].fromRandomSphereSample(); }
    }

    void addRotationForces( const Vec3d& center, const Vec3d& hdir, double K, double R ){
        double R2 = R*R;
        double iR2 = 1/R2;
        for(int i=0; i<npick; i++){
            Vec3d d = apos[i] - center;
            Vec3d f; f.set_cross( hdir, d );
            double r2 = f.norm2();
            if( (r2>R2) || (r2<1e-8) )continue;
            f.mul( (1-r2*iR2)*K );
            apos[i].add( f );
        }
    }

    void deform_L(){
        double dl = L/nstep;
        for(int istep=0; istep<nstep; istep++){
            for(int i=0; i<npick; i++){ apos[i].add_mul( pulls[i], dl); }     // pull the atoms
            evalForce(natoms*3, (const double*)apos, (double*)aforce );                              // eval relaxation forces
            for(int i=0; i<npick;  i++){ aforce[i].set(0.0); }                  //
            for(int i=0; i<natoms; i++){ apos  [i].add_mul( aforce[i], dt ); }  //
        }
    }

    void deform_F(){
        for(int istep=0; istep<nstep; istep++){
            evalForce(natoms*3, (const double*)apos, (double*)aforce );           // eval relaxation forces
            for(int i=0; i<npick;  i++){ aforce[i].add_mul(pulls[picks[i]], L );    }  //
            for(int i=0; i<natoms; i++){ apos  [i].add_mul( aforce[i], dt ); }  //
        }
    }

    void deform_Rot(){
        for(int istep=0; istep<nstep; istep++){
            evalForce(natoms*3, (const double*)apos, (double*)aforce );           // eval relaxation forces
            for(int i=0; i<npick; i++){
                addRotationForces( apos[ picks[i]], pulls[i], 1.0, L );
            }
            for(int i=0; i<natoms; i++){ apos  [i].add_mul( aforce[i], dt ); }  //
        }
    }

    double rotateNeighs( int ia, const Vec3d& hdir, double K ){
        const Vec3d& p0=apos[ia];
        int* ngs=graph->atom2neigh + graph->ngIs[ia];
        int  nng=graph->ngNs[ia];
        //moment;
        for(int i=0; i<nng; i++){
            int ja=ngs[i];
            Vec3d d; d.set_sub  (apos[ja],p0);
            Vec3d f; f.set_cross( hdir,   d );
            aforce[ja].add_mul( f, K );
            //printf( "ia,ja %i,%i f(%g,%g,%g) K %g \n", ia, ja, f.x,f.y,f.z, K  );
        }
        return 0;
    }

    void bondRot( int ib, double K ){
        const Vec2i& b = graph->bond2atom[ib];
        Vec3d hdir=apos[b.a]-apos[b.b];
        hdir.normalize();
        //printf( " === bond[%i](%i,%i) hdir(%g,%g,%g) \n", ib, b.a, b.b, hdir.x,hdir.y,hdir.z );
        rotateNeighs( b.a, hdir,  K );
        rotateNeighs( b.b, hdir, -K );
        //aforce[b.a].add( 0.,0.,+10.);
        //aforce[b.b].add( 0.,0.,-10.);
    }

    void deform_BondRot( double K ){
        for(int istep=0; istep<nstep; istep++){
            evalForce(natoms*3, (const double*)apos, (double*)aforce );           // eval relaxation forces
            for(int i=0; i<npick; i++){
                bondRot( picks[i], Ks[i]*K );
            }
            for(int i=0; i<natoms; i++){ apos[i].add_mul( aforce[i], dt ); }  //
        }
    }

};


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
        // this should spread pull-vector between neighboring atoms depending on bond stiffness
        //  - in principle it is approximation of running force field relaxation, but without re-calculating stiffness matrix and bond directions
        for(int i=0; i<natoms; i++){ dpos_[i]=dpos[i]; }
        for(int ib=0; ib<nbonds; ib++){
            const Vec2i& b = bond2atom[ib];
            if( !(amask[b.a]||amask[b.b])  ) continue;
            Vec3d vk; vk.set_mul( bondDir[ib], bond_k[ib] );
            Vec3d fa; fa.set_div( vk, dirStiff[b.a] );
            Vec3d fb; fb.set_div( vk, dirStiff[b.b] );
            dpos_[b.a].add(fa);
            dpos_[b.b].add(fb);
        }
        return;
    }

    void evalLinearForce( Vec3d* dpos, Vec3d* force ){  // plug this as    CG.dotA(n,x,r)
        //for(int i=0; i<natoms; i++){ amask[i] =( dpos.norm2()<R2min ); }
        for(int ib=0; ib<nbonds; ib++){
            const Vec2i& b = bond2atom[ib];
            // To Do we can set some masks and trasholds to speed things up
            //if( !(amask[b.a]||amask[b.b])  ) continue;
            Vec3d fv;
            fv.set_sub( dpos[b.a], dpos[b.b]   );         // d = di-dj
            // if( fv.norm2()<F2min ) contine;
            const Vec3d& hdir = bondDir[ib];
            fv.set_mul( hdir, hdir.dot(fv)*bond_k[ib] );  // f = k<d|h>h
            force[b.a].add(fv);
            force[b.b].sub(fv);
            //amask[b.a]=true; amask[b.b]=true;
        }
    }

    /*
    void sparseLinearBondSmoother( Vec3d* dpos, Vec3d* dpos_, double trash ){
        // this shoud spread pull-vector between neighboring atoms depending on bond stiffness
        //  - in principle it is approximation of running force field relaxation, but without re-calculating stiffness matrix and bond directions
        for(int ia=0; ia<natoms; ia++){
            for(int ing=0; ing<nneighs[ia]; ing++){
                int    ja=neighs[ing];
                double kb=kneighs[ing];
            }
        }
        return;
    }
    */

    void evalStiffnessSubSpace( int npick, int* iatoms  ){
        // calculate small stiffness matrix for several selected atoms, considering just bonds (no angles etc.)
        // Kij = E/(dxi*dxj)
    }

}; // DirStiff

#endif
