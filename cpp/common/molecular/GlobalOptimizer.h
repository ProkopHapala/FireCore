
#ifndef GlobalOptimizer_h
#define GlobalOptimizer_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Atoms.h"
#include "MultiSolverInterface.h"

class GlobalOptimizer{ public:

    int* atypes=0;
    MultiSolverInterface* msolver =0;
    SolverInterface*      solver  =0;

    int npop=0;
    Atoms**             population=0;
    //std::vector<Atoms*> history;

    Vec3d dpos_min{-0.5,-0.5,-0.5};
    Vec3d dpos_max{+0.5,+0.5,+0.5};

    Mat3d dLvec_min{ -1.0,-1.0,-0.0,   -0.0,-0.1,-0.0,   -0.5,-0.5,-0.5  };
    Mat3d dLvec_max{ +1.0,+1.0,+0.0,   +0.0,+0.1,+0.0,   +0.5,+0.5,+0.5  };

    // ================= Functions

    void reallocPop(int npop_, int nvec, bool bAtypes ){
        npop=npop_;
        if(bAtypes)_realloc(atypes,nvec);
        _realloc( population, npop );
        for(int i=0; i<npop; i++){
            population[i] = new Atoms(nvec,true, false);
            population[i]->atypes=atypes;
        }
    }

    void solve(int ipop ){
        Atoms* atoms = population[ipop];
        solver->setGeom( atoms->apos,  atoms->lvec );   
        solver->solve( 10000, 1.0e-6 );                 
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

    void lattice_scan_1d( int n, Mat3d lvec0, Mat3d dlvec, int ipop0=0, const char* outfname=0 ){
        FILE* fout=0;
        if(outfname){ fout = fopen(outfname,"w"); if(!fout){ printf("ERROR in GlobalOptimizer::lattice_scan_a() cannot open %s \n", outfname ); exit(0); } }
        Mat3d lvec=lvec0;
        Atoms* atoms0 = population[0];
        solver->getGeom( atoms0->apos, atoms0->lvec );
        printf("atoms0->lvec\n"    ); printMat( *(atoms0->lvec) );
        for(int i=0; i<n; i++){
            if(i>0){    
                lvec.add(dlvec);
                population[i]->copyOf( *population[i-1] ); *(population[i]->lvec)=lvec;
                //population[i]->toNewLattic( lvec, population[i-1] );
            }
            printf( "### lattice_scan_1d(%i) lvec{{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}} \n", i, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,  lvec.c.x,lvec.c.y,lvec.c.z  );
            solve( i ); DEBUG
            if(fout){ population[ipop0+i]->atomsToXYZ(fout,true,true, {2,2,1}); } DEBUG
        }
        if(fout){fclose(fout);}
    }

/*
    void lattice_scan_2d( Vec3i ns, Mat3d lvec0, Mat3d dlvec1, Mat3d dlvec2 ){

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
        fclose(fout);
    }
*/

};

#endif
