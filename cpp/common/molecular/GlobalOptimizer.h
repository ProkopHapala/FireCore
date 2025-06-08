
#ifndef GlobalOptimizer_h
#define GlobalOptimizer_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Atoms.h"
#include "MultiSolverInterface.h"
#include "constrains.h"
#include <numeric>
#include "distribution_generator.h"
#include "MolecularDatabase.h"

class GlobalOptimizer{ public:

    // bool bExploring = false;
    // int  istep   =0;

    // double vel_kick = 1.0;
    // double pos_kick = 0.25;

    // double bThermalSampling = 0.0;   // if >0 then we do thermal sampling
    // double T_target         = 300.0;   // target temperature for thermal sampling (if bThermalSampling>0)
    // double T_current        = 0.0;   // current temperature for thermal sampling (if bThermalSampling>0)
    // double gamma_damp       = 0.01;  // damping factor for thermal sampling (if bThermalSampling>0)

 

    // void apply_kick( int n, Vec3d* pos=0, Vec3d* vel=0 ){
    //     //printf( "GOpt::apply_kick() n=%i |pos|=%g |vel=%g|\n", n, pos_kick, vel_kick );
    //     Vec3d vcog = Vec3dZero;
    //     for(int i=0; i<n; i++){
    //         if(pos) pos[i].add( randf(-pos_kick,pos_kick), randf(-pos_kick,pos_kick), randf(-pos_kick,pos_kick) );
    //         if(vel){
    //             Vec3d dv = Vec3d{ randf(-vel_kick,vel_kick), randf(-vel_kick,vel_kick), randf(-vel_kick,vel_kick) };
    //             vel[i].add( dv );
    //             vcog  .add( dv );
    //         }
    //     }
    //     vcog.mul( 1.0/n );
    //     if(vel)for(int i=0; i<n; i++){
    //         vel[i].sub( vcog );
    //     }
    // }

    // void clear(){
    //     bExploring = false;
    //     istep=0;
    //     constrs.clear();
    // }

    // void copy(const GOpt& go){
    //     //printf( "GOpt::copy() \n" );
    //     nExplore = go.nExplore;
    //     nRelax   = go.nRelax;
    //     vel_kick = go.vel_kick;
    //     pos_kick = go.pos_kick;
    //     T_target = go.T_target;
    //     gamma_damp = go.gamma_damp;
    //     constrs.copy( go.constrs );
    //     printf( "GOpt::copy() nExplore=%i nRelax=%i vel_kick=%g pos_kick=%g T_target=%g gamma_damp=%g \n", nExplore, nRelax, vel_kick, pos_kick, T_target, gamma_damp );
    // }


    // void print(){ printf( "GOpt::print() nExplore=%i nRelax=%i vel_kick=%g pos_kick=%g T_target=%g gamma_damp=%g \n", nExplore, nRelax, vel_kick, pos_kick, T_target, gamma_damp ); }
















    //double tolerance = 1.0e-6;
    double tolerance = 0.01;
    int    nmaxiter  = 10000; 
    int    initMode  =0;

    int* atypes=0;
    MultiSolverInterface* msolver =0; // multi solver
    SolverInterface*      solver  =0; // single solver

    //int npop=0;
    //Atoms**             population=0;
    //std::vector<Atoms*> history;

    std::vector<Atoms*> population;

    Vec3d dpos_min{-0.5,-0.5,-0.5};
    Vec3d dpos_max{+0.5,+0.5,+0.5};

    Mat3d dLvec_min{ -1.0,-1.0,-0.0,   -0.0,-0.1,-0.0,   -0.5,-0.5,-0.5  };
    Mat3d dLvec_max{ +1.0,+1.0,+0.0,   +0.0,+0.1,+0.0,   +0.5,+0.5,+0.5  };

    // ================= Functions

    /**
     * Reallocates the population with the specified number of individuals and vectors.
     * 
     * @param npop The number of individuals in the population.
     * @param nvec The number of vectors in each individual.
     * @param bAtypes A boolean flag indicating whether to reallocate the atom types.
     */
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

    /**
     * Downloads the population from the parallel solver and updates the energy of each individual.
     * 
     * @param n The number of individuals to update.
     * @param i0 The starting index of the individuals to update.
     */
    void download_multi(int n, int i0){
        msolver->downloadPop();
        for(int i=0; i<n; i++){
            population[i+i0]->Energy = msolver->getGeom( i, population[i+i0]->apos, 0, true );
        }
    }

    /**
     * Uploads multiple molecular configurations to the parallel solver.
     * 
     * @param n The number of configurations to upload.
     * @param i0 The starting index of the configurations in the population array.
     * @param bGeom Flag indicating whether to upload atomic positions.
     * @param blvec Flag indicating whether to upload lattice vectors.
     */
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

    /**
     * @brief Loads the XYZ file containing atomic coordinates for the population.
     * 
     * @param fname The name of the file to load.
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

    /**
     * Writes the atomic coordinates of the population to an XYZ file.
     * 
     * @param fname The name of the output file.
     * @param imin The index of the first individual in the population to write (default: 0).
     * @param imax The index of the last individual in the population to write (default: -1, which means all individuals).
     * @param nPBC The number of periodic boundary conditions in each direction (default: {1, 1, 1}).
     */
    void popToXYZ( const char* fname, int imin=0, int imax=-1, Vec3i nPBC=Vec3i{1,1,1}){
        FILE* file=0;
        file = fopen(fname,"w");    if(!file){ printf("ERROR in GlobalOptimizer::lattice_scan_a() cannot open %s \n", fname ); exit(0); }
        if(imax<0){ imax=population.size(); }
        for(int i=imin; i<imax; i++){
            if(file){ population[i]->atomsToXYZ(file,true,true,nPBC); }
        }
        fclose(file);
    }

    /**
     * Performs a 1-dimensional lattice scan.
     *
     * @param n The number of iterations for the scan.
     * @param lvec0 The initial lattice vector.
     * @param dlvec The change in lattice vector for each iteration.
     * @param outfname The output file name (optional).
     * @param ipop0 The initial population index (optional).
     * @param istep The step size for population index (optional).
     */
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

    /**
     * Performs a 2D lattice scan for multiple replicas.
     *
     * @param n The number of steps in the lattice scan.
     * @param dlvec The displacement vector for each step.
     * @param initMode The initialization mode for changing the cell of replicas. Default is 0.
     * @param outfname The output file name. Default is nullptr.
     * @param imin The minimum index of replicas to consider. Default is 0.
     * @param imax The maximum index of replicas to consider. Default is -1, which means all replicas.
     */
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



    // void startExploring(){
    //     if(verbosity) printf( "GOpt::startExploring()\n" );
    //     bExploring = true;
    //     //istep=0;
    //     constrs.update_drives();
    //     constrs.update_drives();
    // }

    // bool update(){
    //     //istep++;
    //     if(bExploring){
    //         // if(istep>=nExplore){ 
    //         //     if(verbosity) printf( "GOpt::update() stop exploring istep(%i)>nExplore(%i) \n" );
    //         //     bExploring=false; istep=0; return true; 
    //         // }
    //     }else{
    //        // if(istep>=nRelax  ){ bExploring=true; istep=0; return true; }
    //     }
    //     //bExploring       = go.bExploring; 
    //     //bThermalSampling = bExploring;
    //     //bConstrains      = bExploring;
    //     //if(bExploring) bConverged = false;
    //     return false;
    // }

































    int DoF         = 0;
    int Findex      = 0;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<int> boundaryRules; // 0 - periodic, 1 - mirroring, 2 - free

    double Fstar    = -__DBL_MAX__;
    int maxeval     = 1e5;
    int  nExplore   = 0;
    int  nExploring = 0;
    int  nRelax     = 0;


    bool bExploring = false;

    std::vector<double> x_opt;
    double F_opt    = __DBL_MAX__;

    int neval       = 0;
    int bShow       = 0;  // 0 - no output, 1 - only better, 2 - all
    int bSave       = 0;  // 0 - no output, 1 - only better, 2 - all
    int bDatabase   = 0;  // 0 - no output, 1 - only better, 2 - all

    std::vector<double> currentStructure;
    Atoms* currentAtoms;
    MMFFparams* params;

    MolecularDatabase* database = 0;
    Constrains constrs;
    DistributionGenerator* d = 0;

bool containsNaN(const std::vector<double>& vec) {
    for (const auto& val : vec) {
        if (std::isnan(val)) {
            return true;
        }
    }
    return false;
}


    void init_heur(Atoms* _atoms, MMFFparams* _params, int _Findex, std::vector<double>* _a=0, std::vector<double>* _b=0, std::vector<int>* _boundaryRules=0,
     double _Fstar=-__DBL_MAX__, int _maxeval=0, int _nRelax=500, int _nExploring=__INT_MAX__, int _bShow=1, int _bSave=1, int _bDatabase=0)
    {
        currentAtoms    = _atoms;
        params          = _params;
        Findex          = _Findex;

        if(_a){a = (*_a);}
        if(_b){b = (*_b);}
        if(_boundaryRules){boundaryRules = (*_boundaryRules);}
        generateborder();

        Fstar           = _Fstar;
        maxeval         = _maxeval;
        nRelax          = _nRelax;
        nExploring      = _nExploring;
        nExplore        = _nExploring;
        printf("nExplore %i\n", nExplore);

        
        bShow           = _bShow;
        bSave           = _bSave;
        bDatabase       = _bDatabase;

        DoF             = a.size();
        for(int i=0; i<DoF; i++) x_opt.push_back(NAN);
        F_opt           = __DBL_MAX__;
        bExploring      = true;
        neval           = 0;
        
        structureToVector(_atoms);
        srand(time(NULL));
    }

    void randomBrutal(int index_mut=-1, std::vector<double>* par_mut=0, std::vector<double>* par_alg=0){
            printf("randomBrutal\n");
            //init_heur(_atoms, _params, _Findex, _a, _b, _boundaryRules, _Fstar, _maxeval, _bShow);
            while(bExploring){
                mutate(&currentStructure, index_mut, par_mut);
                double F_new = evaluate(&currentStructure);
            }
        }

// needs par_alg repeater (nb of steps before giving up on a mutation)
    void randomDescent(int index_mut=-1, std::vector<double>* par_mut=0, std::vector<double>* par_alg=0){
            printf("randomDescent\n");
            if(!par_alg){
                par_alg = new std::vector<double>();
                (*par_alg).push_back(1);
            }
            int repeater = (int)(*par_alg)[0];
            std::vector<double> x_new;
            double F_new, F;
            while(bExploring){
                mutate(&currentStructure, -1, nullptr);
                F = evaluate(&currentStructure);
                while (bExploring)
                {
                    int i = 0;
                    for(; i<repeater; i++){
                        mutate(&x_new, index_mut, par_mut);
                        F_new = evaluate(&x_new);
                        if(F_new < F){
                            F = F_new;
                            currentStructure = x_new;
                            break;
                        }
                    }
                    if(i==repeater){
                        break;
                    }
                }
            }
            results();
        }

    std::vector<double>* gradSPSA(std::vector<double>* x, double ck){
        std::vector<double>* grad = new std::vector<double>();
        std::vector<double> delta, x_plus, x_minus;
        double F_plus, F_minus, slope;
        x_plus = (*x);
        x_minus = (*x);
        for(int i=0; i<DoF; i++){
            if(randf(0, 1) > 0.5){
                delta.push_back(1.0);
            }else{
                delta.push_back(-1.0);}
            x_plus[i] += ck*delta[i];
            x_minus[i] -= ck*delta[i];
        }
        F_plus = evaluate(&x_plus);
        F_minus = evaluate(&x_minus);
        slope = (F_plus-F_minus)/(2*ck);
        for(int i=0; i<DoF; i++){
            (*grad).push_back(slope*delta[i]);
        }
        return grad;
    }
    void SPSA(int index_mut=-1, std::vector<double>* par_mut=0, std::vector<double>* par_alg=0){
            printf("SPSA\n");
            printf("DoF %i\n", DoF);
            if(!par_alg){
                par_alg = new std::vector<double>();
                (*par_alg).push_back(0.602);
                (*par_alg).push_back(0.101);
                (*par_alg).push_back(maxeval);
            }
            double alpha = (*par_alg)[0];
            double gamma = (*par_alg)[1];
            double N = (*par_alg)[2];
            std::vector<double> grad;
            mutate(&currentStructure, 2, par_mut);
            double F_new, F;
            double A = N*0.1;
            double c = 1e-2;
            grad = (*gradSPSA(&currentStructure, c));
            double magnitude_g0 = std::reduce(grad.begin(), grad.end())/grad.size();
            magnitude_g0 = abs(magnitude_g0);
            double a = 2*pow(A+1, alpha)/magnitude_g0;
            double ak, ck;
            //printf("beggining of iteration\n");
            for (int k = 0; k < N; k++){
                ak = a/pow(k+A+1, alpha);
                ck = c/pow(k+1, gamma);
if (containsNaN(currentStructure)) {
    std::cout << "Vector contains NaN values. Before grad\n";exit(1);
}                 
                grad = (*gradSPSA(&currentStructure, ck));
if (containsNaN(grad)) {
    std::cout << "Vector contains NaN values. After grad\n";exit(1);
}
                for(int i=0; i<DoF; i++){
                    currentStructure[i] -= ak*grad[i];
                }
if (containsNaN(currentStructure)) {
    std::cout << "Vector contains NaN values. Before checkBorders\n";exit(1);
} 
                checkBorders(&currentStructure);
if (containsNaN(currentStructure)) {
    std::cout << "Vector contains NaN values. After checkBorders\n";exit(1);
} 
            }

    }

    // void simpleLevyFlights(int index_mut=3, std::vector<double>* par_mut=0, std::vector<double>* par_alg=0){
    //         printf("simpleLevyFlights\n");
    //         if(!par_alg){
    //             par_alg = new std::vector<double>();
    //             (*par_alg).push_back(1e-5);
    //             (*par_alg).push_back(100);
    //             (*par_alg).push_back(0.1);
    //             (*par_alg).push_back(0.1);
    //         }
    //         double alpha = (*par_alg)[0];
    //         double T0    = (*par_alg)[1];
    //         double n0    = (*par_alg)[2];
    //         int waiting  = (*par_alg)[3];

    //         mutate(&currentStructure, 0, nullptr);
    //         double F, F_trial;
    //         F = evaluate(&currentStructure);
    //         int k = 0;
    //         double Tmut;
    //         int counter = 0;
    //         std::vector<double> x_trial;
    //         printf("beggining of iteration\n");
    //         while(bExploring){
    //             Tmut = T0/pow(1+k/n0, 1/alpha);
    //             k++;
    //             counter = waiting;
    //             while(counter>0 && bExploring){
    //                 x_trial = currentStructure;
    //                 mutate(&x_trial, 3, nullptr);
    //                 F_trial = evaluate(&x_trial);
    //                 if(F_trial < F){
    //                     currentStructure = x_trial;
    //                     F = F_trial;
    //                     counter = waiting;
    //                 }
    //                 else{
    //                     counter--;
    //                 }
    //             }
    //         }
    //         results();
    //     }



    double evaluate(std::vector<double>* x){
        //printf("evaluate\n");
        if(neval>=maxeval || !bExploring){
            bExploring = false;
            return __DBL_MAX__;
        }

        
        if(nExplore<=0){        // maybe this should be done differently, but for this exploration it works
            Findex = 1;
        }
        double F = calculateEnergy(x);
        if(nExplore<=0){
            Findex = 0;
            nExplore = nExploring;
        }
        nExplore--;
        informUser(F, x, F<F_opt);
        if (F < F_opt)
        {
            F_opt = F;
            x_opt = (*x);
            if (F <= Fstar)
            {
                bExploring = false;
            }
        }       
        return F;
    }

    void informUser(double F, std::vector<double>* x, bool better){
        if(bShow == 2){
            printf("neval %8i F %10.3f", neval, F);
            // for (int i = 0; i < DoF; i++)
            // {
            //     printf(" %10.3f", (*x)[i]);
            // }
            printf("\n");
        }
        else if(bShow == 1 && better){
            printf("neval %8i F %10.3f", neval, F);
            // for (int i = 0; i < DoF; i++)
            // {
            //     printf(" %10.3f", (*x)[i]);
            // }
            printf("\n");
        }
        if(bSave == 2){
            saveXYZ(F, x);
        }
        else if(bSave == 1 && better){
            saveXYZ(F, x);
        }
        if(bDatabase == 2){
            addSnapshot();
        }
        else if(bDatabase == 1 && better){
            addSnapshot();
        }
    }


    void mutate(std::vector<double>* x, int index_mut, std::vector<double>* par_mut=0){
        //printf("mutate\n");
        if ((*x).size() == 0)
        {
            for (int i = 0; i < DoF; i++)
            {
                (*x).push_back(randf(a[i], b[i]));
            }
        }
        std::vector<double> x_new;
        double l = 0;
        double rndf = randf(0, 1);
        int rnd = 1;
        switch (index_mut)
        {
        case 0: // L1 mutation
            for (int i = 0; i < DoF; i++)
            {
                x_new.push_back(randf((*x)[i]-(*par_mut)[0], (*x)[i]+(*par_mut)[0]));
                l += fabs((*x)[i]-x_new[i]);
            }
            
            for (int i = 0; i < DoF; i++)
            {
                (*x)[i] += (x_new[i]-(*x)[i])/l*rndf*(*par_mut)[0];
            }
            break;
        case 1: // Hamming mutation
            rnd = rand()%DoF;
            (*x)[rnd] = randf(a[rnd], b[rnd]);
            break;
        case 2: // Gaussian mutation
            if(!d || d->distributionType != 1){
                if(d) delete d;
                d = new DistributionGenerator();
                d->init(1); // Gaussian distribution
            }
            for (int i = 0; i < DoF; i++)
            {
                rnd = d->generateRandom();
                (*x)[i] += (*par_mut)[0]*rnd;
            }
            break;
        case 3: // Student mutation
            if(!d || d->distributionType != 3){
                if(d) delete d;
                d = new DistributionGenerator();
                d->init(3); // Student's t-distribution with 1 Dof
            }
            l = d->generateRandom();

            d->randomPointOnNdimSphere(x_new, DoF);
            for(int i=0; i<DoF; i++){
                (*x)[i] += (*par_mut)[0]*x_new[i]*l;
            }
            break; 
        default:
                for (int i = 0; i < DoF; i++)
                {
                    (*x)[i] = randf(a[i], b[i]);
                }                
            break;
        }
        
        checkBorders(x);
    }

    void results(){
        printf("Results\n");
        printf("F_opt %f\n", F_opt);
        printf("x_opt\n");
        for(int i=0; i<DoF; i++){
            printf("%f ", x_opt[i]);
        }
        printf("\n");
        vectorToStructure(&x_opt);
    }


    void generateborder(){
        switch (Findex)
        {
        case -1:
            if(!a.size()){for(int i=0; i<3; i++){a.push_back(-10);}}
            if(!b.size()){for(int i=0; i<3; i++){b.push_back( 10);}}
            break;
        case 0 || 1:
            DoF=currentAtoms->natoms*3;
            if(!a.size()){for(int i=0; i<DoF; i++){a.push_back(-DoF);}} // natoms*3A (all atoms in one line times 3A=upper bound for bond length)
            if(!b.size()){for(int i=0; i<DoF; i++){b.push_back(DoF);}}
            break;
        default:
            break;
        }
        if(!boundaryRules.size()){
            for(int i=0; i<a.size(); i++){
                boundaryRules.push_back(0);
            }
        }
    }

    void checkBorders(std::vector<double>* x){
        for(int i=0; i<DoF; i++){
            double range = b[i]-a[i];
            double xi;
            switch (boundaryRules[i])
            {
            case 0: //periodic boundaries
                if((*x)[i] > b[i]){
                    (*x)[i] = a[i]+std::fmod((*x)[i]-a[i], range);
                }
                else if((*x)[i] < a[i]){
                    (*x)[i] = b[i]-std::fmod(b[i]-(*x)[i], range);
                }
                break;
            case 1: //mirroring boundaries
                    if((*x)[i] > b[i]){
                        xi =std::fmod((*x)[i]-a[i], 2*range);
                        (*x)[i] = a[i]+std::min(xi, 2*range-xi);
                    }
                    else if((*x)[i] < a[i]){
                        xi =std::fmod(b[i]-(*x)[i], 2*range);
                        (*x)[i] = b[i]-std::min(xi, 2*range-xi);
                    }
                break;
            default://free boundaries
                break;
            }
        }
    }

    void structureToVector(Atoms* atoms=0){ // missing lvec DoF
        bool notAtoms = false;
        Vec3d* apos = new Vec3d[currentStructure.size()/3];
        Mat3d lvec=Mat3dIdentity;
        if(!atoms){
            notAtoms = true;
            solver->getGeom(apos, &lvec);
            for (int i = 0; i < DoF ; i+=3)
            {
                currentStructure[i] = apos[i/3].x;
                currentStructure[i + 1] = apos[i/3].y;
                currentStructure[i + 2] = apos[i/3].z;
            }
        }
        else{
            for(int i=0; i<atoms->natoms; i++){
                apos[i].x = atoms->apos[i].x;
                apos[i].y = atoms->apos[i].y;
                apos[i].z = atoms->apos[i].z;
            }
            for (int i = 0; i < atoms->natoms; i++)
            {
                currentStructure.push_back(apos[i].x);
                currentStructure.push_back(apos[i].y);
                currentStructure.push_back(apos[i].z);
            }
        }

        delete[] apos;
    }

    void vectorToStructure(std::vector<double>* x=0){ // missing lvec DoF
        if(!x){
            x = &currentStructure;
        }
        Vec3d* apos = new Vec3d[x->size()/3];
        Mat3d lvec=Mat3dIdentity;
        for(int i=0; i<x->size(); i+=3){
            apos[i/3].x = (*x)[i];
            apos[i/3].y = (*x)[i+1];
            apos[i/3].z = (*x)[i+2];
        }
        solver->setGeom(apos, &lvec);
        delete[] apos;
    }

    double function(std::vector<double>* x){ // Findex = -1
        double Energy=0;
        for(int i=0; i < (*x).size(); i++){
            Energy += (*x)[i]*(*x)[i];
        }
        return Energy;
    }

    double calculateEnergy(std::vector<double>* x){
        double E = __DBL_MAX__;
        switch (Findex)
        {
        case -1: // debugging "energy" only a 3D function with one minimum
            return function(x);
            break;
        case 0: // evaluate the energy of the current structure
            vectorToStructure(x);
            neval += 1;
            return solver->solve(1, tolerance);
            break;
        case 1: // minimize the energy of the current structure
            vectorToStructure(x);
            E = solver->solve(nRelax, tolerance);
            structureToVector();
            neval += nRelax;
            return E;
            break;
        default:
            break;
        }
        return __DBL_MAX__;
    }



















    void saveXYZ(double F, std::vector<double>* x){
            char str_tmp[1024];
            bool bPBC = false;
            char comment[256];
            sprintf(comment, "neval %8i F %10.3f", neval, F);
            if(bPBC){ sprintf( str_tmp, "lvs %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f %s", currentAtoms->lvec->a.x, currentAtoms->lvec->a.y, currentAtoms->lvec->a.z, currentAtoms->lvec->b.x, currentAtoms->lvec->b.y, currentAtoms->lvec->b.z, currentAtoms->lvec->c.x, currentAtoms->lvec->c.y, currentAtoms->lvec->c.z, comment ); }
            else    { sprintf( str_tmp, "%s", comment ); }
            Vec3i nPBC = {1,1,1};
            if(!currentAtoms->lvec) {currentAtoms->lvec = new Mat3d();currentAtoms->lvec->set(Mat3dIdentity);}
            params->saveXYZ("gopt_trajectory.xyz", currentAtoms->natoms, currentAtoms->atypes, currentAtoms->apos, str_tmp, 0, "a", true, nPBC, *(currentAtoms->lvec) );
    }

    bool addSnapshot(bool ifNew = false, char* fname = 0)
    {
        if(!database){
            database = new MolecularDatabase();
            database->setDescriptors();
        }
        Atoms* atoms = new Atoms(DoF/3); 
        atoms->Energy=solver->getGeom( atoms->apos,  atoms->lvec );
        if (fname)
        {
            //loadNBmol(fname);
            //buildMolecule_xyz(fname);
        }
        if (ifNew)
        {
            int ID = database->addIfNewDescriptor(atoms);
            if (ID != -1)
            {
                printf("Same as %d\n", ID);
                return false;
            }
        }
        else
        {
            database->addMember(atoms);
        }
        return true;
    }

    void printDatabase()
    {
        if (database)
        {
            database->print();
        }
    }


};

#endif
