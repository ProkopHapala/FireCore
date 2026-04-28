
#ifndef OptRandomWalk_h
#define OptRandomWalk_h

#include "globals.h"
//#include <cstddef>
#include <math.h>
#include "macroUtils.h"
#include "VecN.h"

//typedef void (*ForceFunction)( int n, double * xs, double * dfs );

typedef double (*EnergyFunction)( int n, double * Xs );

class OptRandomWalk{ public:

    // ====  Variables

    int n=0;
    double * X     = 0;
    double * Xbest = 0; // currently best solution
    double * Xdir  = 0; // decent solution at reasonable distance from Xbest (not too far, not too close) which approximates soft-mode direction
    double * scales = 0;

    // limits 
    double * Xmin = 0;
    double * XMax = 0;

    int nTryDir = 0;

    double cDir = 1.2;
    double cBack = 0.2;

    //double biasDir;
    //double biasSpread;

    EnergyFunction getEnergy = 0;

    double stepSize=0.1;
    double Ebest,E;
    double Edir,Rdir; // dir bias, like second best solution

    int  iMutAlg = 0;
    bool bLimit  = true;


    // ====  Functions

    void realloc( int n_, double* X_ = 0 ){
        printf( "OptRandomWalk::realloc    *X %li %li \n", X_, X_+(n_-1) );
        n=n_;
        if(X_){ X=X_; }else{ _realloc(X,n); };
        _realloc(Xbest ,n);
        _realloc(Xdir  ,n);
        _realloc(scales,n);
        for(int i=0; i<n; i++){
            //Xbest[i]=
            Xbest [i]=X[i];
            scales[i]=1;
        }
    }

    void dealloc( double bX=true ){
        _dealloc(X     );
        _dealloc(Xbest );
        _dealloc(Xdir  );
        _dealloc(scales);
    }

    void mutate( double step ){
        double rndDir = randf();
        double cfw    = 0.5 + 1.0/nTryDir;
        double cbak   = 0.5;
        double c      = cfw+cbak;
        switch ( iMutAlg ){
            case -1:{ int i = rand()%n;        X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step; break; }        // random step along single random axis
            case  0:{ for( int i=0; i<n; i++){ X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step; }; break; }     // random step along all axes
            case  1:{ for( int i=0; i<n; i++){ X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*cDir-cBack); }; break; } // random step along soft-mode direction
            default:{
                printf( "OptRandomWalk::mutate() iMutAlg %i not implemented \n", iMutAlg );
                exit(0);
            } 
        }
        //X[i] = X[i] + (randf()-0.5);
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step;
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*rndDir;
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*rndDir;
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*1.2-0.2);
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*3.0-1.0);
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*0.75-0.25);
        //X[i] = Xbest[i] + (randf()-0.5)*scales[i]*step  +  (Xdir[i]-Xbest[i])*(rndDir*c-cbak);
    }

    void limit(){ for(int i=0; i<n; i++){ _clamp(X[i],Xmin[i],XMax[i]);} }

    void run( int nstep, bool bStart=true ){
        if(bStart){ start(); }
        for(int i=0; i<nstep;i++){
            mutate( stepSize );
            if(bLimit)limit();
            E = getEnergy( n, X );
            if(E<Ebest){
                if(verbosity>0)printf( "Energy improved %g -> %g \n", Ebest, E );
                VecN::set(n,X,Xbest);
                //Rdir = VecN::err2( n, X, Xbest );
                Edir = 1e+300;
                Ebest=E;
            }else{
                double R = VecN::err2( n, X   , Xbest );
                Rdir     = VecN::err2( n, Xdir, Xbest );
                double curv    = (E   -Ebest)/(R*R);
                double curvDir = (Edir-Ebest)/(Rdir*Rdir);
                if( curv<curvDir ){
                    //printf( "dir update: (%g|dE %g R %g ) -> ( %g|dE %g R %g ) \n", curvDir,(Edir-Ebest), Rdir,    curv,(E   -Ebest), R    );
                    VecN::set(n,X,Xdir);
                    Rdir=R;
                    Edir=E;
                    nTryDir=0;
                }else{
                    nTryDir++;
                }
            }
        }
    }

    void start(){
        Ebest = getEnergy( n, X );
    }

};

#endif
