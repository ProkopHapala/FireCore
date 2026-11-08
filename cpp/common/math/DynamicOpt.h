
#ifndef DynamicOpt_h
#define DynamicOpt_h

#include <cstdio>// DEBUG
//#include <cstddef>
#include <math.h>
#include "macroUtils.h"

#include <globals.h>

typedef double (*ForceFunction )( int n, double * xs, double * dfs );  // fills dfs with forces; returns energy (may be ignored by methods that only need forces)

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
    bool bOwnPos=true, bOwnVel=true, bOwnForce=true, bOwnInvMasses=true; // false when bound to ffl.DOFs/fDOFs in setOptimizer

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

    // --- opt-in relaxation preconditioner: trajectory-secant diagonal masses (see eFF_optimizer_speedup.md §7.2)
    // Separate stiffness per particle xyz-block and per scalar (width) DOF. Fictitious relaxation masses, NOT physical MD masses.
    int     autoMass    = 0;      // 0=off; 1=block-norm secant EMA; 2=asymmetric log-stiffness (fast-up/slow-down memory)
    int     am_nBlocks  = 0;      // leading 3-component blocks (atoms + electron positions)
    int     am_nScalars = 0;      // trailing scalar DOFs (electron sizes)
    double  am_kref     = 100.0;  // invM = kref/(kref+k_b); k_b starts at 0 => unit masses
    double  am_ema      = 0.1;    // EMA weight of each new stiffness sample
    double  am_emaUp    = 0.2;    // autoMass=2: log-stiffness learning rate when observation is STIFFER  (remembers stiff modes fast)
    double  am_emaDn    = 0.01;   // autoMass=2: log-stiffness learning rate when observation is SOFTER (slowly recovers mobility)
    double* am_k        = 0;      // [am_nBlocks+am_nScalars] smoothed stiffness estimate
    double* am_xold     = 0;      // [n] previous EVALUATED positions
    double* am_fold     = 0;      // [n] previous evaluated forces
    double* am_mass0    = 0;      // [n] base invMasses saved before first autoMass activation; restored when autoMass is disabled (§11.2B: otherwise mode-0 candidates inherit adaptive masses)
    bool    am_init     = false;
    int     massMetric  = 0;      // 0=legacy metric (vf=v.(invM.f), vv=v.v); 1=mass-consistent (P=v.f, V2=v.(v/invM), G2=f.(invM.f)) — §7.3

    // --- opt-in gradient-history optimizers (line-search, need energy => getEnergy must be set) — audit §8
    // lsMethod: 0=off, 1=BB1 spectral GD, 2=Anderson acceleration of x->x-eta*g (multisecant, memory mHist<=mMax)
    //           3=Anderson preconditioned by current invMasses: map x->x-eta*D*g, D=diag(invMasses)
    static const int mMax = 8;                       // absolute maximum Anderson history (typically m<=5)
    int     lsMethod  = 0;
    int     mHist     = 5;
    double  lsEta     = 1e-3;                        // Anderson base step scale
    double  lsAinit   = 1e-3;                        // BB1 first step length
    double  lsArmijo  = 1e-4;                        // Armijo slope coefficient
    int     lsEwin    = 5;                           // nonmonotone window (max of last N accepted energies)
    // online adaptation (all opt-in, zero extra force evals — signals are already on the trajectory):
    int     lsAdapt   = 0;                           // bit1 eta<-EMA(BB1 s.s/s.y), bit2 m_eff shrinks on resets, bit4 Armijo warm-start (a may grow up to lsAmax), bit8 momentum fallback on reset, bit16 drop stale secants on sharp turns
    double  lsAmax    = 8.0;                         // cap on extrapolated Armijo multiplier when lsAdapt&4
    double  lsMomBeta = 0.5;                         // fallback-momentum EMA decay (lsAdapt&8): pm <- beta*pm + (1-beta)*s_accepted
    double  lsMomGain = 1.0;                         // gain of momentum term in the reset fallback direction (lsAdapt&8)
    double  lsTurnCos = 0.5;                         // sharp-turn cosine: accepted-step tangent turns harder than this -> drop oldest half of secant history (lsAdapt&16)
    double  lsDescEps = 1e-3;                        // descent margin for bounded correction (lsAdapt&32): theta makes g.d <= -lsDescEps*q; larger = stricter descent, smaller theta
    // rank-1 stiff-mode suppression (move_FIRE): force component along smoothed transverse-ΔF direction is damped by (1-msMu)
    int     modeSupp   = 0;                          // 0=off, 1=suppress force along smoothed stiff-mode direction u
    double  msMu       = 0.2;                        // residual force fraction along u (1=off, 0=fully suppressed)
    double  msEma      = 0.3;                        // EMA for u direction
    double* msU        = 0;                          // [n] smoothed unit stiff-mode direction (lazy alloc)
    double* msFe       = 0;                          // [n] mode-suppressed effective force
    double* msFold     = 0;                          // [n] force snapshot for ΔF
    double* msXold     = 0;                          // [n] pos  snapshot for Δx
    bool    msInit     = false;
    // gated two-timescale momentum (§15.6): slow-EMA velocity register, boost v only while cos(F,v)>gate_cos
    double  gate_cos   = 0.7;                        // alignment gate threshold
    double  gate_g     = 0.0;                        // boost gain (0=off)
    double  vs_ema     = 0.9;                        // EMA decay for vSlow register
    double* vSlow      = 0;                          // [n] slow-EMA velocity register (lazy alloc when gate_g!=0)
    bool    vSlowInit  = false;
    double  lsEmaCoef = 0.3;                         // EMA coefficient for eta_ad <- BB1 step scale
    double  lsEtaLo   = 0.02;                        // eta_ad bounds, relative to lsEta
    double  lsEtaHi   = 200.0;
    double* xmin      = 0;                           // optional per-DOF lower bounds (caller-owned, e.g. esize floor); 0=none
    double* lsXc=0;   // [n]   accepted state          (allocated once in allocLS, only when lsMethod>0)
    double* lsGc=0;   // [n]   gradient at xc (-force)
    double* lsXp=0;   // [n]   previous accepted state
    double* lsGp=0;   // [n]   previous gradient
    double* lsD =0;   // [n]   search direction
    double* lsPm=0;   // [n]   fallback momentum: EMA of accepted steps (lsAdapt&8)
    double* lsV =0;   // [n]   work vector
    double* lsS =0;   // [n*mMax] accepted step history (Anderson)
    double* lsY =0;   // [n*mMax] gradient-diff history
    double* lsQ =0;   // [n*mMax] QR basis workspace
    double* lsR =0;   // [mMax*mMax] triangular factor
    double* lsB =0;   // [mMax] rhs/gamma
    int*    lsPiv=0;  // [mMax] pivot map
    double  lsEwin_[5];                              // energy ring buffer

    ForceFunction getForce = 0;   // fills force; returns energy (unused by FIRE/GD, required by run_GH line search)

    // ==== function declarations

    double move_LeapFrog( double dt_loc, const double* farr=0 );
    //void   move_LeapFrog_vlimit();
    double move_GD      ( double dt_loc );
    //double move_GD_safe ( double dt_loc );
    double move_MD( double dt_loc, double damp);
    double move_MD_dbg( double dt_loc, double damping_factor );
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
    void   updateAutoMasses();
    void   allocLS();
    int    run_GH( int nstepMax, double Fconv, double* outE, double* outF );

    double getFmaxAbs( );
    double getFsqSum( );

    // Fat checkpoint (§13.4/§13.8: branch experiments need identical complete states). Serializes all scalars
    // and the CONTENTS of every state buffer into a flat double blob; restoreState copies contents back into the
    // currently-bound arrays (pos/vel/force/invMasses are never re-bound — they may alias force-field buffers).
    int  stateSize ( );                                    // number of doubles required for the blob
    int  saveState ( double* buf );                        // returns stateSize, or -1 on mismatch
    void restoreState( const double* buf );                // fails loudly on magic/n mismatch

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

    // Fresh-state contract for autoMass (§11.2B): save the base invMasses on first activation, restore on
    // disable — otherwise a mode-0 run inherits the previous candidate's adaptive masses (A->B->A contamination).
    inline void saveMassesOnce(){ if(!am_mass0){ _realloc(am_mass0,n); for(int i=0;i<n;i++)am_mass0[i]=invMasses[i]; } }
    inline void restoreMasses(){ if(am_mass0){ for(int i=0;i<n;i++)invMasses[i]=am_mass0[i]; delete[] am_mass0; am_mass0=0; } }

    inline void bindArrays( int n_, double * pos_, double * vel_, double * force_, double * invMasses_ ){
        n = n_; pos=pos_;  vel=vel_; force=force_; invMasses=invMasses_;
        bOwnPos=false; bOwnVel=false; bOwnForce=false; bOwnInvMasses=false;
    }

    inline void bindOrAlloc( int n_, double * pos_, double * vel_, double * force_, double * invMasses_ ){
        n = n_;
        if(pos_  ==0){ _realloc(pos  ,n); bOwnPos=true;   }else{ pos   = pos_;        bOwnPos=false;   };
        if(vel_  ==0){ _realloc(vel  ,n); bOwnVel=true;   }else{ vel   = vel_;        bOwnVel=false;   };
        if(force_==0){ _realloc(force,n); bOwnForce=true; }else{ force = force_;      bOwnForce=false; };
        if(invMasses_==0) { _realloc(invMasses,n); bOwnInvMasses=true; setInvMass(1.0); }else{ invMasses=invMasses_; bOwnInvMasses=false; }
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
        bOwnPos=true; bOwnVel=true; bOwnForce=true; bOwnInvMasses=true;

        _realloc(avs  ,n);
        _realloc(avs2 ,n);
    }

    inline void dealloc( ){ // never free borrowed pos/force (alias ffl.DOFs); avs/avs2 always owned here
        if(bOwnPos      ){ _dealloc(pos);       }else{ pos=0;       }
        if(bOwnVel      ){ _dealloc(vel);       }else{ vel=0;       }
        if(bOwnForce    ){ _dealloc(force);     }else{ force=0;     }
        if(bOwnInvMasses){ _dealloc(invMasses); }else{ invMasses=0; }
        _dealloc(avs);
        _dealloc(avs2);
        _dealloc(am_k); _dealloc(am_xold); _dealloc(am_fold); _dealloc(am_mass0); am_init=false;
        _dealloc(lsXc); _dealloc(lsGc); _dealloc(lsXp); _dealloc(lsGp); _dealloc(lsD); _dealloc(lsV); _dealloc(lsPm);
        _dealloc(msU); _dealloc(msFe); _dealloc(msFold); _dealloc(msXold);
        _dealloc(vSlow); vSlowInit=false;
        _dealloc(lsS);  _dealloc(lsY);  _dealloc(lsQ);  _dealloc(lsR);  _dealloc(lsB); _dealloc(lsPiv);
        n=0; bOwnPos=true; bOwnVel=true; bOwnForce=true; bOwnInvMasses=true;
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
        vSlowInit=false;
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
