
#include "DynamicOpt.h" // THE HEADER

//#include <cstdio>// DEBUG
//#include "VecN.h"

#include "fastmath.h"

// ===============  MoveSteps

void DynamicOpt::move_LeapFrog(double dt_loc){
    //double dt_ = dt*fscale_safe;
    for ( int i=0; i<n; i++ ){
        //printf( "i %i v %g f %g p %g iM %g \n", i, vel[i],force[i],pos[i],invMasses[i]  );
        double v = vel[i];
        v += force[i]*invMasses[i]*dt_loc;
        pos[i]  += v*dt_loc;
        vel[i]   = v;
        //printf("DynamicOpt::move_LeapFrog[%i] invM %g \n", i, invMasses[i] );
    }
    stepsDone++;
    t += dt_loc;
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

void DynamicOpt::move_GD(double dt_loc){
    //double dt_ = dt*fscale_safe;
    for ( int i=0; i<n; i++ ){
        pos[i] += force[i]*dt_loc;
    }
    stepsDone++;
    t += dt_loc;
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



double DynamicOpt::move_VSpread(double dt_loc,double damp, double beta ){
    //   damping should be higher when     Sum_i{ |v_i|^2 } >> |Sum{v_i}|^2 
    double f2sum=0;
    for ( int i=0; i<n; i++ ){
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
    eval_cos_vf();
	double cv;
    double cf  = renorm_vf * damp_func( cos_vf, cv );
    for(int i=0; i<n; i++){ vel[i]  = vel[i]*cv  + force[i]*invMasses[i]*cf;  }
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

double DynamicOpt::move_FIRE(){
    eval_cos_vf();
	double cv;
    double cf  = renorm_vf * damp_func_FIRE( cos_vf, cv );
    for(int i=0; i<n; i++){ vel[i]  = vel[i]*cv  + force[i]*invMasses[i]*cf;  }
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
        double fi = fabs( force[i] );
        fmax=(fi>fmax)?fi:fmax;
    }
    return fmax;
}

double DynamicOpt::getFsqSum( ){
    double ff = 0;
    for(int i=0; i<n; i++){
        double fi = force[i];
        ff += fi*fi;
    }
    return ff;
}
