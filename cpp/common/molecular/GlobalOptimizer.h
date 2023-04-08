
#ifndef GlobalOptimizer_h
#define GlobalOptimizer_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Atoms.h"
#include "MultiSolverInterface.h"

class GlobalOptimizer{ public:

    MultiSolverInterface* msolver =0;
    SolverInterface*      solver  =0;

    int npop=0;
    Atoms**             population=0;
    //std::vector<Atoms*> history;

    Vec3d dpos_min{-0.5,-0.5,-0.5};
    Vec3d dpos_max{+0.5,+0.5,+0.5};

    Mat3d dLvec_min{ {-1.0,-1.0,-0.0},  {-0.0,-0.1,-0.0},  {-0.5,-0.5,-0.5}  };
    Mat3d dLvec_max{ {+1.0,+1.0,+0.0},  {+0.0,+0.1,+0.0},  {+0.5,+0.5,+0.5}  };

    // ================= Functions

    void reallocPop(int npop_, int nvec ){
        npop=npop_;
        _realloc( population, ntot );
        for(int i=0; i<npop; i++){
            population[i] = new Atoms(nvec,true);
        }
    }

    void solve(int ipop ){
        Atoms atoms = population[ipop];
        solver->setGeom( atoms->apos,  atoms->lvec );
        solver->solve();
        solver->getGeom( atoms->apos,            0 );
    }

    void download_multi(int npara){
        msolver->downloadPop();
        for(int i=0; i<npara; i++){
            //solver->getGeom( i, population[i]->apos, population[i]->lvec, true );
            msolver->getGeom( i, population[i]->apos,                   0, true );
        }
    }

    void upload_multi(int npara){
        for(int i=0; i<npara; i++){
            msolver->setGeom( i, population[i]->apos, population[i]->lvec, true );
        }
        msolver->uploadPop();
    }

    void download_multi(int npara){
        for(int i=0; i<npara; i++){
            //solver->getGeom( i, population[i]->apos, population[i]->lvec, true );
            msolver->getGeom( i, population[i]->apos,                   0, true );
        }
        msolver->uploadPop();
    }

/*
    void toNewLattic( int ipop, const Mat3d& lvec ){
        Atoms* atoms = population[ipop];
        Mat3d invLvec;
        atom->lvec->lvec.invert_T_to( invLvec );
        Vec3d* apos = atoms->apos;
        for(int i=0;i<atoms->natoms;i++){
            Vec3d u; invLvec.dot_to( apos[i], u );
            lvec.dot_to_T( u, apos[i] )
        }
        *(atoms->lvec) = lvec;
    }
*/

    void lattice_scan_a( Vec3i ns, Mat3d lvec, Vec3d da1, Vec3d da2 ){

        FILE* fout = fopen("lattice_scan_a.xyz","w");
        if(!fout){ printf("ERROR in GlobalOptimizer::lattice_scan_a() fopen failed \n"); exit(0); }

        int npara = solver->paralel_size();
        if( npop!==npara )reallocPop( npara );
        int ipop=0;
        //int ipara = 0;
        // ---- first we scan 1st DOF using serial sover (CPU) 
        for(int i=0; i<ns.a; i++){
            if(ia>0){
                lvec.a.add(da1);
                population[i]->toNewLattic( lvec, population[i-1] );
            }
            solve      ( ia );
            population[i]->atomsToXYZ (fout,true,true);
        }
        // ---- Then we scan 2nd DOF using parallel solver (GPU), we start form  
        for(int j=0;j<ns.b;j++){ // loop over paralle steps
            for(int i=0;i<npop;i++){
                lvec = *(population[ia]->lvec);
                lvec.a.add( da2 );
                population[i]->toNewLattic( lvec );
            }
            upload();
            msolver->sove_multi();
            download();
            for(int ia=0;ia<ns.a;ia++){
                population[i]->atomsToXYZ(fout,true,true);
            }
        }
    }
    fclose(fout);
};

#endif
