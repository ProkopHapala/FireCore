
#ifndef DynamicOpt_h
#define DynamicOpt_h

#include <cstdio>// DEBUG
//#include <cstddef>
#include <math.h>
#include "macroUtils.h"

#include <globals.h>

typedef void (*ForceFunction)( int n, double * xs, double * dfs );

class OptLog { public:
    int n=0;
    double* cos  = 0;
    double* f    = 0;
    double* v    = 0;
    double* dt   = 0;
    double* damp = 0;

    void set(int i, double ci, double fi, double vi, double dti, double Di ){
        if(n==0)return;
        i = i%n;
        if(cos )cos [i] = ci;
        if(v   )v   [i] = vi;
        if(f   )f   [i] = fi;
        if(dt  )dt  [i] = dti;
        if(damp)damp[i] = Di;
    }

};


class DynamicOpt{ public:
    // variables
    //int verbosity = 0;
    int n=0;
    double * pos       = 0;
    double * vel       = 0;
    double * force     = 0;
    double * invMasses = 0;

    bool  bfixmask=false;
    bool* fixmask = 0;

    double * avs =0;
    double * avs2=0;


    // parameters
    double dt           = 0.1;
    double damping      = 0.2;

    //double cvf_min = -0.5;  // minimum cosine for velocity damping in  move_FIRE_smooth()
    //double cvf_max = +0.5;  // maximum cosine for velocity damping in  move_FIRE_smooth()

    double cvf_min = -0.1;  // minimum cosine for velocity damping in  move_FIRE_smooth()
    double cvf_max = +0.1;  // maximum cosine for velocity damping in  move_FIRE_smooth()
    double cv_kill =  0.0;

    //double cvf_min = -0.30;  // minimum cosine for velocity damping in  move_FIRE_smooth()
    //double cvf_max =  0.05;  // maximum cosine for velocity damping in  move_FIRE_smooth()

    //double cvf_min = -0.0;  // minimum cosine for velocity damping in  move_FIRE_smooth()
    //double cvf_max =  0.1;  // maximum cosine for velocity damping in  move_FIRE_smooth()

    //double cvf_min = -0.01;  // minimum cosine for velocity damping in  move_FIRE_smooth()
    //double cvf_max = +0.01;  // maximum cosine for velocity damping in  move_FIRE_smooth()

    double fTminmax     = 0.1;
    double f_limit      = 10.0;
    double v_limit      = 10.0;
    //double l_limit      = 0.1;
    double dr_limit     = 0.2;
    //double fscale_safe  = 1;
    double scale_dt  = 1;
    double ff_safety    = 1e-32;

    // FIRE
    //int    minLastNeg   = 5;  // Original paper
    int    minLastNeg   = 3;
    double finc         = 1.1;
    double fdec         = 0.5;
    //double falpha       = 0.98; // Original paper
    double falpha       = 0.8;
    double kickStart    = 1.0;

    double dt_max       = dt;
    double dt_min       = dt*fTminmax;
    double damp_max     = damping;

    int    lastNeg      = 0;

    // other
    int method    = 2;
    int stepsDone = 0;
    double t      = 0.0;

    double ff=0,vv=0,vf=0;

    double f_len      = 0;
    double v_len      = 0;
    double cos_vf     = 0;
    double renorm_vf  = 0;
    double cv,cf;


    ForceFunction getForce = 0;

    // ==== function declarations

    double move_LeapFrog( double dt_loc );
    //void   move_LeapFrog_vlimit();
    double move_GD      ( double dt_loc );
    //double move_GD_safe ( double dt_loc );
    double move_MD( double dt_loc, double damp);
    double move_VSpread(double dt_loc,double damp, double beta );
    //double move_MD_safe ( double dt_loc );
    inline double move_MDquench(){ return move_MD(dt,damping);};
    double move_FIRE();
    double move_FIRE_bak();
    double move_FIRE_smooth();
    double damp_func     ( double c, double& cv );
    double damp_func_FIRE( double c, double& cv );
    void   FIRE_update_params();
    double optStep();
    bool   optimize( double convF, int nMaxSteps );

    double getFmaxAbs( );
    double getFsqSum( );

    // ==== inline functions

    inline void project_vf(){
        ff=0,vv=0,vf=0;
        for(int i=0; i<n; i++){
            if(bfixmask){ if(fixmask[i])continue; }
            //double fi = force[i];
            double fi = force[i]*invMasses[i];
            double vi = vel[i];
            ff += fi*fi;
            vv += vi*vi;
            vf += vi*fi;
            //printf( "move_FIRE %i f %g v %g p %g \n", i, force[i], vel[i], pos[i] );
        } 
    }

    inline void eval_cos_vf(){
        project_vf();
        f_len      = sqrt(ff);
        v_len      = sqrt(vv);
        cos_vf     = vf    / ( f_len*v_len + ff_safety );
        renorm_vf  = v_len / ( f_len       + ff_safety );
    }

    inline void setInvMass(double invM){  if(invMasses==0){ _realloc(invMasses,n);}  for(int i=0;i<n;i++){invMasses[i]=invM;} };

    inline void bindArrays( int n_, double * pos_, double * vel_, double * force_, double * invMasses_ ){
        n = n_; pos=pos_;  vel=vel_; force=force_; invMasses=invMasses_;
    }

    inline void bindOrAlloc( int n_, double * pos_, double * vel_, double * force_, double * invMasses_ ){
        n = n_;
        if(pos_  ==0){ _realloc(pos  ,n); }else{ pos   = pos_;   };
        if(vel_  ==0){ _realloc(vel  ,n); }else{ vel   = vel_;   };
        if(force_==0){ _realloc(force,n); }else{ force = force_; };
        if(invMasses_==0) { _realloc(invMasses,n); setInvMass(1.0); }else{ invMasses=invMasses_; }
        //if(invMasses_==0) setInvMass(1.0);

        _realloc(avs  ,n);
        _realloc(avs2 ,n);
        cleanAV   ( );
    }

    inline void realloc( int n_ ){
        n = n_;
        _realloc(pos    ,n);
        _realloc(vel    ,n);
        _realloc(force  ,n);
        _realloc(invMasses,n);

        _realloc(avs  ,n);
        _realloc(avs2 ,n);
    }

    inline void dealloc( ){
        _dealloc(pos);
        _dealloc(vel);
        _dealloc(force);
        _dealloc(invMasses);
    }

    inline void cleanForce( ){  for(int i=0; i<n; i++){ force[i]=0; } }
    inline void cleanVel  ( ){  for(int i=0; i<n; i++){ vel  [i]=0; } }
    inline void cleanAV   ( ){  for(int i=0; i<n; i++){ avs[i]=0; avs2[i]=0; } }

    // f ~ sqrt(k/m)    dpos ~ f*dt*dt
    inline double limit_dt_x2 (double xx,double xmax){ double sc=1.0; if( xx > (xmax*xmax) ){ sc= fmin( sc, sqrt(xmax/sqrt(xx)) ); }; return sc;       }
    //inline double limit_dt_x2 (double xx,double xmax){ double sc=1.0; if( xx > (xmax*xmax) ){ sc= fmin( sc, xmax/sqrt(xx) ); }; return sc;       }
    inline double limit_dt_vf2(double ff, double vv ){ scale_dt=fmin(limit_dt_x2(ff,f_limit),limit_dt_x2(vv,v_limit));         return scale_dt; }

    inline void setTimeSteps(double dt_  ){ dt = dt_max = dt_; dt_min=dt_max*fTminmax; }
    inline void setDamping  (double damp_){ damping = damp_max = damp_; }

    inline void initOpt( double dt_, double damp_=0.1 ){
        setTimeSteps(dt_ );
        setDamping  (damp_);
        cleanForce( );
        cleanVel  ( );
    }

    inline double spreadDamp( int i, double v, double beta ){
        double mb = 1-beta;
        double v2   = v*v;
        double av   = avs [i]*mb + v   * beta;
        double av2  = avs2[i]*mb + v*v * beta;
        double cdamp = 1;
        //if(i==1)printf(  "spreadDamp[%i] %g ratio %g av(%g,%g) v(%g,%g) \n", i, cdamp,   av*av/(av2+1e-32), av,av2, v,v2  );
        //if(i==1)printf(  "%i   %15.10f    %15.10f %15.10f     %15.10f %15.10f \n", i, av*av/(av2+1e-32), av,av2, v,v2  );
        if(i==1)printf(  "%i   %15.10f    %15.10f %15.10f     %15.10f \n", i, av*av/(av2+1e-32), av*av,av2, v  );
        avs [i]     = av;
        avs2[i]     = av2;
        return cdamp;
    }

};

#endif
