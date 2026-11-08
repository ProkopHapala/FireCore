
#include "DynamicOpt.h" // THE HEADER

//#include <cstdio>// DEBUG
//#include "VecN.h"

#include "fastmath.h"
#include "Lingebra.h"   // leastSquareFit_QR (Anderson direction)
#include <string.h>     // memcpy
#include <cmath>        // isfinite
#include <cfloat>       // DBL_MAX (rejected-trial sentinel in outF, §11.2A)

// ===============  MoveSteps

double DynamicOpt::move_LeapFrog(double dt_loc, const double* farr){
    //double dt_ = dt*fscale_safe;
    //printf("DynamicOpt::move_LeapFrog() dt %g \n", dt_loc);
    if(!farr)farr=force;
    double f2sum=0;
    for ( int i=0; i<n; i++ ){
        if(bfixmask){ if(fixmask[i])continue; }
        //printf( "%3i %c v %12.5e f %12.5e p %12.5e iM %12.5e \n", i, (fixmask[i])?'.':' ', vel[i],force[i],pos[i],invMasses[i]  );
        double f = farr[i];
        double v = vel[i];
        f2sum+=f*f;
        v += f*invMasses[i]*dt_loc;
        pos[i]  += v*dt_loc;
        vel[i]   = v;
        //printf("DynamicOpt::move_LeapFrog[%i] invM %g \n", i, invMasses[i] );
    }
    stepsDone++;
    t += dt_loc;
    return f2sum;
}

/*
void DynamicOpt::move_LeapFrog_vlimit(){
    double dtv = dt*fscale_safe;
    double vmax = 0.0d;
    for ( int i=0; i<n; i++ ){
        double v = vel[i] + invMasses[i]*dtv*force[i];
        vmax = fmax( fabs(v), vmax );
        vel[i]=v;
    }
    double dtp = dt;
    if( vmax>v_limit ) dtp=v_limit/vmax;
    //printf("vmax %g dtp  %g dtv %g\n", vmax);
    for ( int i=0; i<n; i++ ){
        pos[i] += dtp*vel[i];
    }
    stepsDone++;
    t += dt;
}
*/

double DynamicOpt::move_GD(double dt_loc){
    //double dt_ = dt*fscale_safe;
    double f2sum=0;
    for ( int i=0; i<n; i++ ){
        if(bfixmask){ if(fixmask[i])continue; }
        double f = force[i];
        double v = vel[i];
        f2sum+=f*f;
        v += f*invMasses[i]*dt_loc;
        pos[i]  += v*dt_loc;
        vel[i]   = v;
    }
    stepsDone++;
    t += dt_loc;
    return f2sum;
}

/*
double DynamicOpt::move_GD_safe(double dt_loc){
    double fmax = VecN::absmax(n,force);
    scale_dt = fmin( 1, f_limit/fmax );
    dt_loc*=scale_dt;
    move_GD(dt_loc);
    stepsDone++;
    t += dt_loc;
    return fmax;
}
*/

double DynamicOpt::move_MD(double dt_loc,double damp){
    double cdamp = 1 - damp; if(cdamp<0)cdamp=0;
    //printf( "DynamicOpt::move_MD() cdamp %g \n", cdamp );
    double f2sum=0;
    for ( int i=0; i<n; i++ ){
        if(bfixmask){ if(fixmask[i])continue; }
        double f = force[i];
        double v = vel[i];
        f2sum+=f*f;
        v*= cdamp;
        v+= invMasses[i]*(f*dt_loc);
        pos[i] += v*dt_loc;
        vel[i]  = v;
    }
    stepsDone++;
    t += dt;
    return f2sum;
}

double DynamicOpt::move_MD_dbg(double dt_loc,double damping_factor){
    // Mirror OpenCL localMD integrator: v = v*damping + f*dt ; p = p + v*dt
    // Here we keep invMasses scaling consistent with other integrators in this class.
    double f2sum=0;
    for ( int i=0; i<n; i++ ){
        if(bfixmask){ if(fixmask[i])continue; }
        double fi = force[i]*invMasses[i];
        double vi = vel[i]*damping_factor + fi*dt_loc;
        pos[i] += vi*dt_loc;
        vel[i]  = vi;
        f2sum  += fi*fi;
    }
    stepsDone++;
    t += dt_loc;
    return f2sum;
}



double DynamicOpt::move_VSpread(double dt_loc,double damp, double beta ){
    //   damping should be higher when     Sum_i{ |v_i|^2 } >> |Sum{v_i}|^2 
    double f2sum=0;
    for ( int i=0; i<n; i++ ){
        if(bfixmask){ if(fixmask[i])continue; }
        double f = force[i];
        double v = vel  [i];
        f2sum+=f*f;
        v*= 1. - damp*spreadDamp( i, v, beta ); 
        v+= invMasses[i]*(f*dt_loc);
        pos[i] += v*dt_loc;
        vel[i]  = v;
    }
    stepsDone++;
    t += dt;
    return f2sum;
}


/*
double DynamicOpt::move_MD_safe(double dt_loc){
    double fmax = VecN::absmax(n,force);
    scale_dt = fmin(1,fmin( v_limit/VecN::absmax(n,vel), f_limit/fmax ));
    dt_loc*=scale_dt;
    move_MD(dt_loc);
    return fmax;
}
*/

double DynamicOpt::move_FIRE_bak(){
    // ToDo: according to FIRE implementation in LAMMPS we should update velocity by force first !!!!
    // see:   https://github.com/lammps/lammps/blob/730e5d2e64106f3e5357fd739b44c7eec19c7d2a/src/min_fire.cpp#L393
    project_vf();
    //printf( "DEBUG move_FIRE %i |f| %g |v| %g <v|f> %g n %i\n", sqrt(ff), sqrt(vv), vf, n );
	if( vf < 0.0 ){
	//if( (vf<0.0)||bOverLimit){
		//dt       = dt * fdec;
		dt       = fmax( dt * fdec, dt_min );
	  	damping  = damp_max;
		lastNeg  = 0;
		cleanVel  ( );
		//for(int i=0; i<n; i++){ vel[i] = kickStart*dt*force[i]; }
		//for(int i=0; i<n; i++){ vel[i] = force[i] * 0.5*sqrt(vv/(ff+ff_safety)); }
		//for(int i=0; i<n; i++){ vel[i] = dmax*force[i]*sqrt(1/ff)/dt_var; }
	}else{
		double cf  =      damping * sqrt(vv/(ff+ff_safety));
		//double cf     =     damping * sqrt(vv/ff);
		double cv     = 1 - damping;
		for(int i=0; i<n; i++){
            if(bfixmask){ if(fixmask[i])continue; }
            //vel[i]    = cv * vel[i]  + cf * force[i];
			vel[i]    = cv * vel[i]  + cf * force[i]*invMasses[i];
		}
		if( lastNeg > minLastNeg ){
			dt        = fmin( dt * finc, dt_max );
			damping   = damping  * falpha;
		}
		lastNeg++;
	}
	//printf( "DEBUG 5.5.3 \n" );

	double dt_=dt;

    //if( ff>(f_limit*f_limit )){
    //    dt_*=sqrt(f_limit/sqrt(ff));
    //   //printf( "force too large: %g => limit dt: %g \n", f, dt_ );
    //};

    if( ff>(f_limit*f_limit) ){
        double f = sqrt(ff);
        /*
        if( ff>(100*f_limit*f_limit) ){
            cleanVel();
            move_GD( dr_limit/f ); // do GD step of length == l_limit
            return ff;
        }
        */
        //printf( "force too large: %g => limit dt: %g \n", f, dt_ );
        dt_*=sqrt( f_limit/f );
    };

    /*
    // dr_limit
    // dr = (v+f*dt)*dt
    {
        double v=sqrt(vv);
        double f=sqrt(ff+ff_safety);
        double x1,x2;
        if( ((v+f*dt)*dt) > dr_limit*4.0 ){
            quadratic_roots( f+ff_safety, v, -dr_limit, x1, x2 );

            printf( "v %g f %g  dt %g dr %g  | roots %g %g \n", v, f, dt, ((v+f*dt)*dt), x1,x2 );
            dt_=x2;
        }
    }
    */
    //printf( "DynamicOpt::move_FIRE dt %g \n", dt);

    //if(verbosity>1) printf( "dt %g damp %g n+ %i | cfv %g |f| %g |v| %g \n", dt,damping, lastNeg, vf/sqrt(vv*ff), sqrt(ff), sqrt(vv) );
    if(verbosity>1) printf( "dt %7.5f damp %7.5f n+ %4i | cfv %7.5f |f| %12.5e |v| %12.5e \n", dt,damping, lastNeg, vf/sqrt(vv*ff), sqrt(ff), sqrt(vv) );
    
    move_LeapFrog( dt_ );
	//move_LeapFrog();
	//move_LeapFrog_vlimit();  // this does not seem to help

	//printf( " %i f v vf  %f %f %f   dt damp  %f %f \n",  stepsDone,   sqrt(ff), sqrt(vv), vf/sqrt(vv*ff),   dt_var, damp_var  );
	//stepsDone++;
	return ff;
}

double DynamicOpt::damp_func( double c, double& cv ){
    double cf;
    if      (c < cvf_min){
        cv = 0.;
        cf = 0.;
    }else if(c > cvf_max){
        cv = 1-damping;
        cf =   damping;
    }else{  // cos(v,f) from [ cvf_min .. cvf_max ]
        double f = (c-cvf_min)/(cvf_max-cvf_min);
        cv = (1.-damping)*f;
        cf =     damping *f;
    }
    return cf;
}

double DynamicOpt::damp_func_FIRE( double c, double& cv ){ // Original FIRE damping
    double cf;
    if( c < 0.0 ){
        cv = cv_kill;
        cf = 0.;
    }else{  
        cv = 1-damping;
        cf = damping;
        //cv = 1.;
        //cf = 0.;
    }
    return cf;
}

double DynamicOpt::move_FIRE_smooth(){
    if(autoMass) updateAutoMasses();
    eval_cos_vf();
	double cv;
    double cf  = renorm_vf * damp_func( cos_vf, cv );
    for(int i=0; i<n; i++){ 
        if(bfixmask){ if(fixmask[i])continue; }
        vel[i]  = vel[i]*cv  + force[i]*invMasses[i]*cf;  
    }
	if( cos_vf < 0.0 ){
        dt       = fmax( dt * fdec, dt_min );
		lastNeg  = 0;
        damping  = damp_max;
	}else{
		if( lastNeg > minLastNeg ){
			dt        = fmin( dt * finc, dt_max );
			damping   = damping  * falpha;
		}
		lastNeg++;
	}
	double dt_=dt;
    if( ff>(f_limit*f_limit) ){
        double f = sqrt(ff);
        dt_*=sqrt( f_limit/f );
    };
    if(verbosity>1) printf( "dt %7.5f damp %7.5f n+ %4i | cfv %7.5f |f| %12.5e |v| %12.5e \n", dt,damping, lastNeg, vf/sqrt(vv*ff), sqrt(ff), sqrt(vv) );
    move_LeapFrog( dt_ );
	return ff;
}

void DynamicOpt::FIRE_update_params(){
    f_len      = sqrt(ff);
    v_len      = sqrt(vv);
    cos_vf     = vf    / ( f_len*v_len + ff_safety );
    renorm_vf  = v_len / ( f_len       + ff_safety );
    cf  = renorm_vf * damp_func_FIRE( cos_vf, cv );
    if( cos_vf < 0.0 ){
        dt       = fmax( dt * fdec, dt_min );
		lastNeg  = 0;
        damping  = damp_max;
	}else{
		if( lastNeg > minLastNeg ){
			dt        = fmin( dt * finc, dt_max );
			damping   = damping  * falpha;
		}
		lastNeg++;
	}
}

// Trajectory-secant diagonal preconditioner (relaxation-only fictitious masses, opt-in).
// One scalar stiffness per leading 3-component position block (atoms + electron centers)
// and one per trailing scalar DOF (electron sizes) — xyz and size are preconditioned separately.
// k_b = |dF_b|/|dx_b| between two consecutively EVALUATED states, EMA-smoothed and clamped;
// invM_b = kref/(kref+k_b), k_b init 0 => unit masses. No extra force evaluations.
void DynamicOpt::updateAutoMasses(){
    int nb = am_nBlocks + am_nScalars;
    if( am_nBlocks*3 + am_nScalars != n ){ printf("ERROR DynamicOpt::updateAutoMasses: block layout %i*3 + %i != n=%i\n", am_nBlocks, am_nScalars, n); exit(1); }
    if(!am_init){
        _realloc(am_k,nb); _realloc(am_xold,n); _realloc(am_fold,n);
        for(int i=0;i<nb;i++)am_k[i]=0;
        for(int i=0;i<n ;i++){ am_xold[i]=pos[i]; am_fold[i]=force[i]; }
        // k_b=0 => formula gives unit inverse masses; write them explicitly so a fresh am_init reset
        // does not silently inherit adaptive masses left over from a previous run() (§8.5 review note)
        for(int i=0;i<n ;i++){ if(bfixmask && fixmask[i])continue; invMasses[i]=am_kref/(am_kref+0.0); }
        am_init=true; return;
    }
    for(int b=0;b<nb;b++){
        int i0 = (b<am_nBlocks) ? b*3 : am_nBlocks*3 + (b-am_nBlocks);
        int ni = (b<am_nBlocks) ? 3   : 1;
        double dx2=0, df2=0;
        for(int k=0;k<ni;k++){
            int i=i0+k; if(bfixmask && fixmask[i])continue;   // free subspace only (partially constrained particles ok)
            double dx = pos[i]-am_xold[i];
            double df = force[i]-am_fold[i];
            dx2 += dx*dx; df2 += df*df;
        }
        if( dx2 > 1e-12 ){                                    // |dx_b| > 1e-6 => measurable step
            double knew = fmin( fmax( sqrt(df2/dx2), 1e-2 ), 1e7 );
            if( autoMass==2 ){                                // asymmetric log-stiffness memory: stiffening is learned fast (eta_up),
                double zold = log( am_k[b] + 1e-12 );         // softening is forgotten slowly (eta_dn) => remebers transient stiff modes
                double dz   = log( knew ) - zold;             // without being permanently immobilized by one bad observation
                am_k[b] = exp( zold + (dz>0?am_emaUp:am_emaDn)*dz );
            }else{
                am_k[b] = (1.0-am_ema)*am_k[b] + am_ema*knew;
            }
        }
        double invM = am_kref/(am_kref + am_k[b]);
        for(int k=0;k<ni;k++){ int i=i0+k; if(bfixmask && fixmask[i])continue; invMasses[i]=invM; }
    }
    for(int i=0;i<n;i++){ am_xold[i]=pos[i]; am_fold[i]=force[i]; }
}

double DynamicOpt::move_FIRE(){
    if(autoMass) updateAutoMasses();
    const double* fe = force;
    if(modeSupp){   // rank-1 stiff-mode suppression: u = EMA of normalized transverse ΔF (q = ΔF - proj_s(ΔF)); fe = f - (1-mu)(f·u)u
        if(!msInit){ _realloc(msU,n); _realloc(msFe,n); _realloc(msFold,n); _realloc(msXold,n);
            for(int i=0;i<n;i++){ msU[i]=0; msFold[i]=force[i]; msXold[i]=pos[i]; } msInit=true; }
        else{
            double ss=0, dqs=0, df2=0;
            for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i])continue; double s_=pos[i]-msXold[i], df_=force[i]-msFold[i]; ss+=s_*s_; dqs+=df_*s_; df2+=df_*df_; }
            if( ss>1e-24 && df2>1e-24 ){
                double q2=0; for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i]){msFe[i]=0;continue;} double q_= (force[i]-msFold[i]) - (dqs/ss)*(pos[i]-msXold[i]); msFe[i]=q_; q2+=q_*q_; }
                if( q2 > 1e-12*df2 ){ double iq=1.0/sqrt(q2), w=1.0-msEma, u2=0;
                    double qu=0; for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i])continue; qu+=msU[i]*msFe[i]; }
                    double sgn = (qu<0) ? -1.0 : 1.0;          // §11.6: u and -u are the same projector u*u^T; align sign so opposite observations don't cancel in the EMA
                    for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i]){msU[i]=0;continue;} msU[i]=msEma*msU[i]+w*sgn*msFe[i]*iq; u2+=msU[i]*msU[i]; }
                    double iu=1.0/sqrt(u2+1e-300); for(int i=0;i<n;i++)msU[i]*=iu;
                }
            }
            double fu=0; for(int i=0;i<n;i++)fu+=force[i]*msU[i];
            for(int i=0;i<n;i++)msFe[i]=force[i]-(1.0-msMu)*fu*msU[i];
            fe=msFe;
        }
        for(int i=0;i<n;i++){ msFold[i]=force[i]; msXold[i]=pos[i]; }
    }
    if(massMetric){
        // Mass-consistent metric (§7.3): P=v.f  V2=sum(v*v/invM)  G2=sum(f*invM*f).
        // With invM==1 this is identical to the legacy eval_cos_vf() reduction.
        double P=0,V2=0,G2=0;
        for(int i=0;i<n;i++){
            if(bfixmask && fixmask[i])continue;
            double fi=force[i], vi=vel[i], mi=invMasses[i];
            P += vi*fi; V2 += vi*vi/mi; G2 += fi*fi*mi;
        }
        vf=P; vv=V2; ff=G2;
        v_len = sqrt(V2); f_len = sqrt(G2);
        cos_vf    = vf    / ( f_len*v_len + ff_safety );
        renorm_vf = v_len / ( f_len       + ff_safety );
    }else{
        eval_cos_vf();
    }
	double cv;
    double cf  = renorm_vf * damp_func_FIRE( cos_vf, cv );
    for(int i=0; i<n; i++){
        //printf( "DynamicOpt::move_FIRE() [%i] fix %i \n", i, fixmask[i] );
        if(bfixmask){ if(fixmask[i])continue; }
        vel[i]  = vel[i]*cv  + fe[i]*invMasses[i]*cf;
    }
	if( cos_vf < 0.0 ){
        dt       = fmax( dt * fdec, dt_min );
		lastNeg  = 0;
        damping  = damp_max;
	}else{
		if( lastNeg > minLastNeg ){
			dt        = fmin( dt * finc, dt_max );
			damping   = damping  * falpha;
		}
		lastNeg++;
	}
	double dt_=dt;
    if( ff>(f_limit*f_limit) ){
        double f = sqrt(ff);
        dt_*=sqrt( f_limit/f );
    };
    // dr_limit: cap the step so that the largest per-DOF displacement |dr_i|=|(v_i+f_i*invM_i*dt)*dt| <= dr_limit.
    // The whole velocity is scaled by one common factor => the step direction is preserved (dr'_i = s*dr_i).
    {
        double drmax=0;
        for(int i=0;i<n;i++){ if(bfixmask && fixmask[i])continue; double dr=(vel[i]+fe[i]*invMasses[i]*dt_)*dt_; drmax=fmax(drmax,fabs(dr)); }
        if( drmax > dr_limit ){
            double s = dr_limit/drmax;
            for(int i=0;i<n;i++){ if(bfixmask && fixmask[i])continue; vel[i] = s*(vel[i]+fe[i]*invMasses[i]*dt_) - fe[i]*invMasses[i]*dt_; }
        }
    }
    if(verbosity>1) printf( "dt %7.5f damp %7.5f n+ %4i | cfv %7.5f |f| %12.5e |v| %12.5e \n", dt,damping, lastNeg, vf/sqrt(vv*ff), sqrt(ff), sqrt(vv) );
    move_LeapFrog( dt_, fe );
	return ff;
}

void DynamicOpt::allocLS(){
    // preallocate all gradient-history buffers once (mMax is a fixed small constant <=8)
    _realloc(lsXc,n); _realloc(lsGc,n); _realloc(lsXp,n); _realloc(lsGp,n); _realloc(lsD,n); _realloc(lsV,n); _realloc(lsPm,n);
    _realloc(lsS ,n*mMax); _realloc(lsY,n*mMax); _realloc(lsQ,n*mMax);
    _realloc(lsR ,mMax*mMax); _realloc(lsB,mMax); _realloc(lsPiv,mMax);
}

// Gradient-history optimizers (audit §8): lsMethod 1 = BB1 spectral GD, 2 = Anderson multisecant.
// Needs getForce callback returning energy (fills force). Raw free-DOF force convergence; step cap
// dr_limit on largest per-DOF displacement; optional lower bounds xmin (projected/shortened BEFORE eval);
// nonmonotone Armijo vs max of last lsEwin accepted energies; every trial evaluation counted.
// Returns number of energy/force evaluations.
int DynamicOpt::run_GH( int nstepMax, double Fconv, double* outE, double* outF ){
    if(!getForce){ printf("ERROR DynamicOpt::run_GH: getForce callback not set\n"); exit(1); }
    if( lsMethod<1 || lsMethod>3 ){ printf("ERROR DynamicOpt::run_GH: lsMethod=%i (1=BB1, 2=Anderson, 3=Anderson+diagPrecond)\n", lsMethod); exit(1); }
    if( lsMethod>=2 && (mHist<1||mHist>mMax) ){ printf("ERROR DynamicOpt::run_GH: mHist=%i > mMax=%i\n", mHist, mMax); exit(1); }
    if(!autoMass) restoreMasses();                 // fresh-state contract: don't inherit a previous candidate's adaptive invMasses (§11.2B)
    if(!lsXc)allocLS();
    double *xc=lsXc,*gc=lsGc,*xp=lsXp,*gp=lsGp,*d=lsD,*v=lsV,*S=lsS,*Y=lsY,*pm=lsPm;
    for(int i=0;i<n;i++)pm[i]=0;                 // fallback momentum starts cold each call
    double *Ewin=lsEwin_; int nEwin=0, iEwin=0, nhist=0;
    const double F2conv=Fconv*Fconv;
    auto free_F2=[&](){ double s=0; for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i])continue; double f=force[i]; s+=f*f; } return s; };
    auto set_g =[&](double*g){ for(int i=0;i<n;i++){ g[i]=(bfixmask&&fixmask[i])?0.0:-force[i]; } };
    auto eval_E=[&](){ return getForce(n,pos,force); };   // nonfinite allowed: handled per-context (fatal at anchor, Armijo-reject during trials)
    double Ecur=eval_E(), F2=free_F2();
    if(!std::isfinite(Ecur)){ printf("ERROR run_GH: non-finite E=%g at eval 0\n",Ecur); exit(1); }
    if(!std::isfinite(F2))  { printf("ERROR run_GH: non-finite force at eval 0\n"); exit(1); }
    int it=0;
    if(outE)outE[it]=Ecur; if(outF)outF[it]=sqrt(F2); it++;
    Ewin[nEwin++]=Ecur; iEwin=nEwin%lsEwin;
    for(int i=0;i<n;i++){ xc[i]=pos[i]; xp[i]=xc[i]; }
    set_g(gc); for(int i=0;i<n;i++)gp[i]=gc[i];
    bool bConv=(F2<F2conv), bFail=false;
    int nResets=0, nRejects=0, nBounded=0, nLowTheta=0, kcLast=0;
    double eta_ad=lsEta, a_last=1.0;                 // adaptive step scale + warm-start Armijo memory
    double rfracLast=0, thetaLast=1;                 // diagnostics: LS residual fraction |r|/|g|, bounded-correction factor
    int    m_eff=mHist, nClean=0;                    // effective Anderson rank (adaptive, lsAdapt&2)
    bool   bWasReset=false;
    while( !bConv && it<nstepMax ){
        double a;
        bWasReset=false;
        if(lsMethod==1){                            // ---- BB1 spectral step: a = (s.s)/(s.y)
            double ss=0,sy=0;
            for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i])continue; double s_=xc[i]-xp[i], y_=gc[i]-gp[i]; ss+=s_*s_; sy+=s_*y_; }
            a = (sy>1e-20) ? ss/sy : lsAinit;
            a = fmin( fmax(a,1e-6), 0.1 );
            for(int i=0;i<n;i++)d[i]=-gc[i];        // d = -g : descent
        }else{                                      // ---- Anderson: d = d0 + theta*c; d0=-eta*D*g, c=-(S-eta*D*Y)*gamma, r=g-Y*gamma (D=I or diag(invMasses) for lsMethod=3)
            bool bPC = (lsMethod==3);
            int kc = Lingebra::leastSquareFit_QR( n, nhist, Y, gc, lsB, lsQ, lsR, lsPiv, lsV, 1e-8 );
            kcLast=kc;
            double *c=lsV;                           // vwork free after the solve -> reuse as correction scratch
            double gg=0, qq=0, gcc=0, rr=0, gd=0;
            for(int i=0;i<n;i++){
                double eta_i = bPC ? eta_ad*invMasses[i] : eta_ad;
                double cy=0, cs=0; for(int k=0;k<kc;k++){ int j=lsPiv[k]; cy+=Y[j*n+i]*lsB[j]; cs+=S[j*n+i]*lsB[j]; }
                double r_=gc[i]-cy; rr+=r_*r_;
                c[i] = eta_i*cy - cs;                // correction c = -(S - eta*D*Y)*gamma ; theta=1 reproduces the original direction
                d[i] = -eta_i*gc[i] + c[i];
                qq  += gc[i]*eta_i*gc[i];            // q = g.(eta*D).g > 0 on the free subspace
                gcc += gc[i]*c[i];
                gd  += gc[i]*d[i]; gg += gc[i]*gc[i];
            }
            rfracLast = sqrt(rr/(gg+1e-300));
            double dd=0; for(int i=0;i<n;i++)dd+=d[i]*d[i]; dd=sqrt(dd); gg=sqrt(gg);
            if( gd >= -0.001*gg*dd ){                // proposed direction is not descent
                bWasReset=true; if(lsAdapt&2){ if(m_eff>1)m_eff--; nClean=0; }
                if( (lsAdapt&32) && qq>0 ){          // §11.4: bound the secant correction (theta) instead of erasing history
                    double theta = (gcc>0) ? fmin( 1.0, (1.0-lsDescEps)*qq/gcc ) : 1.0;   // ensures g.d <= -lsDescEps*q
                    thetaLast=theta; if(theta<1.0)nBounded++;
                    for(int i=0;i<n;i++){ double eta_i = bPC ? eta_ad*invMasses[i] : eta_ad; d[i]=-eta_i*gc[i] + theta*c[i] + ((lsAdapt&8)?lsMomGain*pm[i]:0.0); }
                    double gd2=0,dd2=0; for(int i=0;i<n;i++){ gd2+=gc[i]*d[i]; dd2+=d[i]*d[i]; }                 // descent guard incl. momentum part
                    if( gd2 >= -0.001*gg*sqrt(dd2) ){ theta=0; for(int i=0;i<n;i++){ double eta_i = bPC ? eta_ad*invMasses[i] : eta_ad; d[i]=-eta_i*gc[i]; } }
                    int drop = (theta<0.1) ? nhist/2 : ((nLowTheta++ >= 3) ? 1 : 0);                                 // catastrophic direction: drop half; chronic bounding: drip-forget oldest column
                    if(drop>0){ for(int jj=0;jj<nhist-drop;jj++){ memcpy(S+jj*n,S+(jj+drop)*n,n*sizeof(double)); memcpy(Y+jj*n,Y+(jj+drop)*n,n*sizeof(double)); }
                        nhist-=drop; nResets++; if(theta<0.1)nLowTheta=0; }
                }else{                               // legacy: non-descent -> clear ALL history; fallback = -eta*g (+ momentum pm when lsAdapt&8, guarded)
                    nhist=0; nResets++;
                    for(int i=0;i<n;i++){ double eta_i = bPC ? eta_ad*invMasses[i] : eta_ad; d[i]=-eta_i*gc[i] + ((lsAdapt&8)?lsMomGain*pm[i]:0.0); }
                    if(lsAdapt&8){ double gd2=0,dd2=0; for(int i=0;i<n;i++){ gd2+=gc[i]*d[i]; dd2+=d[i]*d[i]; }
                        if( gd2 >= -0.001*gg*sqrt(dd2) ) for(int i=0;i<n;i++){ double eta_i = bPC ? eta_ad*invMasses[i] : eta_ad; d[i]=-eta_i*gc[i]; } }
                }
            }
            a = (lsAdapt&4) ? fmin( 1.3*a_last, lsAmax ) : 1.0;   // warm start: extrapolate along a good-but-short direction
        }
        // ---- caps: largest |d_i*a| <= dr_limit; lower bounds xmin (project at-bound comps, shorten otherwise)
        double dmax=0; for(int i=0;i<n;i++){ double ad=fabs(d[i]); if(ad>dmax)dmax=ad; }
        if(dmax>0) a=fmin(a, dr_limit/dmax);
        if(xmin){ for(int i=0;i<n;i++){ if(d[i]<0){ if(xc[i]<=xmin[i])d[i]=0; else a=fmin(a,(xc[i]-xmin[i])/(-d[i])); } } }
        // ---- nonmonotone Armijo line search (all trials counted)
        double Emax=Ewin[0]; for(int k=1;k<nEwin;k++)Emax=fmax(Emax,Ewin[k]);
        double gd=0; for(int i=0;i<n;i++)gd+=gc[i]*d[i];
        bool ok=false; double Et=0;
        while(a>1e-14 && it<nstepMax){                 // budget caps evals inside the line search too (was: it overshot nstepMax)
            for(int i=0;i<n;i++)pos[i]=xc[i]+a*d[i];
            Et=eval_E(); it++;
            double Ft=free_F2();
            bool acc = std::isfinite(Et) && std::isfinite(Ft) && ( Et <= Emax + lsArmijo*a*gd );   // require finite ENERGY and FORCES; nonfinite = outside domain -> reject
            if(outE&&it<=nstepMax)outE[it-1]=Et;
            if(outF&&it<=nstepMax)outF[it-1]= acc ? sqrt(Ft) : DBL_MAX;                          // §11.2A: rejected trials must never count as convergence
            if(acc){ ok=true; break; }
            a*=0.5; nRejects++;
        }
        if(!ok){ for(int i=0;i<n;i++)pos[i]=xc[i]; if(it>=nstepMax)break; bFail=true; break; }   // restore last accepted state; budget-out != line-search fail
        // ---- accept: online adaptation (free signals), then update state/history/energy window
        double F2t=free_F2();
        if(lsAdapt&1){                                  // eta <- EMA(BB1 s.s/s.y) — 1-D optimal scale along -g, zero extra evals
            double ss=0,sy=0;
            for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i])continue; double s_=pos[i]-xc[i], y_=-force[i]-gc[i]; ss+=s_*s_; sy+=s_*y_; }
            if(sy>0){ eta_ad += lsEmaCoef*( ss/sy - eta_ad ); }
            eta_ad = fmin( fmax( eta_ad, lsEta*lsEtaLo ), lsEta*lsEtaHi );
        }
        if( (lsAdapt&2) && !bWasReset ){ nClean++; if(nClean>=8 && m_eff<mHist){ m_eff++; nClean=0; } }   // reset already shrank m_eff+nClean; regrow on clean streak
        if( lsAdapt&(8|16) ){
            double ss=0,sp=0,pp=0;
            for(int i=0;i<n;i++){ if(bfixmask&&fixmask[i])continue; double st=pos[i]-xc[i], pv=xc[i]-xp[i]; ss+=st*st; sp+=st*pv; pp+=pv*pv; }
            if(lsAdapt&8){ double w=1.0-lsMomBeta; for(int i=0;i<n;i++) pm[i]=lsMomBeta*pm[i]+w*(pos[i]-xc[i]); }   // EMA of accepted steps = valley-following momentum
            if( (lsAdapt&16) && nhist>1 && pp>0 && sp < lsTurnCos*sqrt(ss*pp) ){                                     // sharp turn => old secants stale: drop oldest half before they produce non-descent directions
                int drop=nhist/2;
                for(int jj=0;jj<nhist-drop;jj++){ memcpy(S+jj*n,S+(jj+drop)*n,n*sizeof(double)); memcpy(Y+jj*n,Y+(jj+drop)*n,n*sizeof(double)); }
                nhist-=drop;
            }
        }
        a_last = a;
        if(lsMethod>=2){ int j; if(nhist==m_eff){ j=m_eff-1; for(int jj=0;jj<m_eff-1;jj++){ memcpy(S+jj*n,S+(jj+1)*n,n*sizeof(double)); memcpy(Y+jj*n,Y+(jj+1)*n,n*sizeof(double)); } } else j=nhist++;
            for(int i=0;i<n;i++){ bool fx=(bfixmask&&fixmask[i]); S[j*n+i]=fx?0.0:a*d[i]; Y[j*n+i]=fx?0.0:(-force[i]-gc[i]); } }   // §11.2D: projected histories — masked reaction force must not enter the secant fit
        for(int i=0;i<n;i++){ xp[i]=xc[i]; gp[i]=gc[i]; xc[i]=pos[i]; }
        if(autoMass) updateAutoMasses();                 // feeds invMasses -> preconditioned Anderson (lsMethod=3); same secant-stiffness machinery as move_FIRE
        Ecur=Et; Ewin[iEwin]=Ecur; iEwin=(iEwin+1)%lsEwin; if(nEwin<lsEwin)nEwin++;
        set_g(gc);
        F2=F2t;
        if(!std::isfinite(F2)){ printf("ERROR run_GH: non-finite force\n"); exit(1); }
        bConv=(F2<F2conv);
    }
    if(verbosity>0) printf("run_GH(lsMethod=%i m=%i eta=%g adapt=%i) %s evals=%i E=%.8f |F|=%.3e resets=%i bounded=%i rejects=%i eta_ad=%g m_eff=%i kc=%i rfrac=%.3f theta=%.2f\n", lsMethod, mHist, lsEta, lsAdapt, bConv?"CONVERGED":(bFail?"LINE-SEARCH-FAIL":"BUDGET-OUT"), it, Ecur, sqrt(F2), nResets, nBounded, nRejects, eta_ad, m_eff, kcLast, rfracLast, thetaLast );
    return it;
}

double DynamicOpt::optStep(){
    //cleanForce( );
    getForce( n, pos, force );
    switch( method ){
        //case 0: move_LeapFrog(dt);
        case 0: move_GD(dt);     break;
        case 1: move_MDquench(); break;
        case 2: move_FIRE();     break;
    }
    return getFmaxAbs( );
}

bool DynamicOpt::optimize( double convF, int nMaxSteps ){
    for( int i=0; i<nMaxSteps; i++ ){
        double f = optStep();
        if( f < convF ) return true;
    }
    return false;
}

// =============== common rutines

double DynamicOpt::getFmaxAbs( ){
    double fmax = 0;
    for(int i=0; i<n; i++){
        if(bfixmask){ if(fixmask[i])continue; }
        double fi = fabs( force[i] );
        fmax=(fi>fmax)?fi:fmax;
    }
    return fmax;
}

double DynamicOpt::getFsqSum( ){
    double ff = 0;
    for(int i=0; i<n; i++){
        if(bfixmask){ if(fixmask[i])continue; }
        double fi = force[i];
        ff += fi*fi;
    }
    return ff;
}
