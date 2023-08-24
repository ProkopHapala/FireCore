
#ifndef GlobalOptimizer_h
#define GlobalOptimizer_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Atoms.h"
#include "MultiSolverInterface.h"

class GlobalOptimizer{ public:
    //double tolerance = 1.0e-6;
    double tolerance = 0.01;
    int    nmaxiter  = 10000; 
    int    initMode  =0;

    int* atypes=0;
    MultiSolverInterface* msolver =0;
    SolverInterface*      solver  =0;

    //int npop=0;
    //Atoms**             population=0;
    //std::vector<Atoms*> history;

    std::vector<Atoms*> population;

    Vec3d dpos_min{-0.5,-0.5,-0.5};
    Vec3d dpos_max{+0.5,+0.5,+0.5};

    Mat3d dLvec_min{ -1.0,-1.0,-0.0,   -0.0,-0.1,-0.0,   -0.5,-0.5,-0.5  };
    Mat3d dLvec_max{ +1.0,+1.0,+0.0,   +0.0,+0.1,+0.0,   +0.5,+0.5,+0.5  };

    // ================= Functions

    void reallocPop(int npop, int nvec, bool bAtypes ){
        //npop=npop_;
        if(bAtypes)_realloc(atypes,nvec);
        //_realloc( population, npop );
        population.resize(npop);
        for(int i=0; i<npop; i++){
            population[i] = new Atoms(nvec,true, false);
            population[i]->atypes=atypes;
        }
    }

    void setGeom(int i){ Atoms* atoms = population[i];               solver->setGeom( atoms->apos,  atoms->lvec ); };
    void getGeom(int i){ Atoms* atoms = population[i]; atoms->Energy=solver->getGeom( atoms->apos,  atoms->lvec ); };

    void solve(int ipop ){
        Atoms* atoms = population[ipop];
        solver->setGeom( atoms->apos,  atoms->lvec );   
        solver->solve  ( nmaxiter, tolerance       );                 
        atoms->Energy = solver->getGeom( atoms->apos, 0 );   
        atoms->id=ipop;
    }

    void download_multi(int n, int i0){
        msolver->downloadPop();
        for(int i=0; i<n; i++){
            population[i+i0]->Energy = msolver->getGeom( i, population[i+i0]->apos, 0, true );
        }
    }

    void upload_multi(int n, int i0, bool bGeom, bool blvec){
        for(int i=0; i<n; i++){
            Vec3d* apos=0;
            Mat3d* lvec=0;
            if(bGeom)apos=population[i+i0]->apos; 
            if(blvec)lvec=population[i+i0]->lvec;
            msolver->setGeom( i, apos, lvec, true );
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

    void loadPopXYZ( const char* fname ){
        FILE* file = fopen(fname,"r");    if(!file){ printf("ERROR in GlobalOptimizer::lattice_scan_a() cannot open \n", fname ); exit(0); }
        while(true){
            //printf( "GlobalOptimizer::loadPopXYZ[%i]\n", population.size() );
            Atoms* atoms = new Atoms(file);
            if( atoms->natoms <=0 ){ break; };
            population.push_back( atoms );
        }
        fclose(file);
    }

    void popToXYZ( const char* fname, int imin=0, int imax=-1, Vec3i nPBC=Vec3i{1,1,1}){
        FILE* file=0;
        file = fopen(fname,"w");    if(!file){ printf("ERROR in GlobalOptimizer::lattice_scan_a() cannot open %s \n", fname ); exit(0); }
        if(imax<0){ imax=population.size(); }
        for(int i=imin; i<imax; i++){
            if(file){ population[i]->atomsToXYZ(file,true,true,nPBC); }
        }
        fclose(file);
    }

    void lattice_scan_1d( int n, Mat3d lvec0, Mat3d dlvec, const char* outfname=0, int ipop0=0, int istep=1 ){
        FILE* fout=0;
        if(outfname){ fout = fopen(outfname,"w"); if(!fout){ printf("ERROR in GlobalOptimizer::lattice_scan_1d() cannot open %s \n", outfname ); exit(0); } }
        Mat3d lvec=lvec0;
        Atoms* atoms0 = population[ipop0];
        atoms0->Energy = solver->getGeom( atoms0->apos, atoms0->lvec );
        //printf("atoms0->lvec\n"    ); printMat( *(atoms0->lvec) );
        for(int i=0; i<n; i++){
            int ipop=i*istep+ipop0;
            long t0=getCPUticks();
            int io=ipop0;
            if(i>0){ lvec.add(dlvec); io=ipop-istep; }
            switch(initMode){
                case 0:{ population[ipop]->copyOf(            *population[io] ); *(population[ipop]->lvec)=lvec; }break;
                case 1:{ population[ipop]->toNewLattice( lvec, population[io] );                                 }break;
            }
            //printf( "lattice_scan_1d[%i] ipop %i lvec\n",  i,ipop ); printMat( lvec );
            printf( "lattice_scan_1d[%i] ipop %i lvec(%6.2f,%6.2f,%6.2f) \n",  i,ipop, lvec.a.x,lvec.a.y,lvec.a.z ); //printMat( lvec );
            //printf( "### lattice_scan_1d(%i) lvec{{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}} \n", i, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,  lvec.c.x,lvec.c.y,lvec.c.z  );
            solve( ipop );
            //printf( "time lattice_scan_1d[%i] %g[Mtick]\n", i, (getCPUticks()-t0)*1e-6 );
            if(fout){ population[ipop]->atomsToXYZ(fout,true,true, {2,2,1}); }
        }
        if(fout){fclose(fout);}
    }

    void lattice_scan_2d_multi( int n, Mat3d dlvec, int initMode=0, const char* outfname=0, int imin=0, int imax=-1 ){
        char comment[256];
        FILE* fout=0;
        if(outfname){ fout = fopen(outfname,"w"); if(!fout){ printf("ERROR in GlobalOptimizer::lattice_scan_2d_multi() cannot open %s \n", outfname ); exit(0); } }
        if(imax<0){ imax=population.size(); };
        int nmult = imax-imin;
        int npara = msolver->paralel_size(); if( nmult!=npara ){ printf("ERROR in GlobalOptimizer::lattice_scan_1d_multi(): (imax-imin)=(%i) != solver.paralel_size(%i) => Exit() \n", nmult, npara ); exit(0); }
        //std::vector<Mat3d> lvec0s;/if(bReturn0){ for(int i=imin;i<imax;i++){ lvec0s[i]= *(population[i]->lvec); } } 
        for(int j=0;j<n;j++){ // main loop ove lvec steps
            if(j>0){
                for(int ipop=imin;ipop<imax;ipop++){ // change cell of all replicas
                    Mat3d new_lvec = *(population[ipop]->lvec)  +  dlvec;   // ToDo: We may want to use lvec-trajectroy different for each of the replicas
                    switch(initMode){
                        case 0:{ *(population[ipop]->lvec) =       new_lvec;   }break;
                        case 1:{   population[ipop]->toNewLattice( new_lvec ); }break;
                    }
                }
            }
            upload_multi  (nmult,imin, initMode!=0, true ); 
            msolver->solve_multi( nmaxiter, tolerance );
            download_multi(nmult,imin);

            //atomsToXYZ(FILE* file, bool bN=false, bool bComment=false, Vec3i nPBC=Vec3i{1,1,1}, const char* comment="", bool bEnergy=true ){
            if(fout) for(int ipop=imin;ipop<imax;ipop++){  sprintf(comment, "step %i replica %i ",  j, ipop  );   population[ipop]->atomsToXYZ(fout,true,true, {2,2,1}, comment, true );  }
        }
        fclose(fout);
        //if(bReturn0){ for(int i=imin;i<imax;i++){ *(population[i]->lvec) = lvec0s[i]; } }
    }


};

#endif
