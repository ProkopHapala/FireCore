
#ifndef ForceField_h
#define ForceField_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Atoms.h"

#include "Forces.h"
#include <functional>


// check if "vals" are within limits "vmin","vmax"
bool checkLimits( int n, int m, const double* vals, const double* vmin, const double* vmax, const char* message, bool bPrint=true ){
    //for(int j=0; j<m; j++){ printf( "checkLimits[%i] [%g,%g]\n", j, vmin[j], vmax[j] ); }
    bool b=false;
    for(int i=0; i<n; i++){
        const double* vali = vals+i*m;
        for(int j=0; j<m; j++){
            double v = vali[j];
            if( v<vmin[j] || v>vmax[j] || isnan(v) ){  
                b=true; 
                if(bPrint){
                    printf( "%s[%i/%i,%i/%i] %g out of limits [%g,%g] \n", message, i,n, j,m,   v, vmin[j], vmax[j]  );
                }
            }
        }
    }
    return b;
}

// damping functions for velocity and forces depending on cos(angle)=<v|f>  (v,f) - vectors of velocity and force
inline double cos_damp_lin( double c, double& cv, double D, double cmin, double cmax  ){
    double cf;
    if      (c < cmin){
        cv = 0.;
        cf = 0.;
    }else if(c > cmax){
        cv = 1-D;
        cf =   D;
    }else{  // cos(v,f) from [ cmin .. cmax ]
        double u = (c-cmin)/(cmax-cmin);
        cv = (1.-D)*u;
        cf =     D *u;
    }
    return cf;
}

class CollisionDamping{ public:
    bool    bBond = false; // if true we use collision damping
    bool    bAng  = false;
    bool    bNonB = false;  // if true we use collision damping for non-bonded interactions

    int     nstep  = 10;    // how many steps it takes to decay velocity to to 1/e of the initial value

    double  medium = 0.0;   // cdamp       = 1 -(damping_medium     /ndampstep     )
    double  bond   = 0.0;   // col_damp    =     collisionDamping   /(dt*ndampstep )
    double  ang    = 0.0;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    double  nonB   = 0.0;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    
    double  dRcut1  = -0.2;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]
    double  dRcut2  =  0.3;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]

    double cdampB    = 0.0;  //  collisionDamping   /(dt*ndampstep );
    double cdampAng  = 0.0;
    double cdampNB   = 0.0;  //  collisionDamping_NB/(dt*ndampstep );

    double cos_vf_acc = 0.5;
    int nstep_acc     = 0;
    int nstep_acc_min = 20;

    inline bool   canAccelerate(){ return nstep_acc>nstep_acc_min; }
    inline double tryAccel     (){ if(    nstep_acc>nstep_acc_min  ){ return 1-medium;  }else{ return 1.0; } }

    inline double update( double dt ){
        if( bBond ){ cdampB   = bond /(dt*nstep ); }else{ cdampB=0;   }
        if( bBond ){ cdampAng = ang  /(dt*nstep ); }else{ cdampAng=0; }
        if( bNonB ){ cdampNB  = nonB /(dt*nstep ); }else{ cdampNB=0;  }
        double  cdamp = 1-(medium/nstep );  if(cdamp<0)cdamp=0;
        //printf( "update_collisionDamping(dt=%g,nstep=%i|medium=%g,bond=%g,ang=%g,nonB=%g) cdamp=%g cdampAng=%g cdampNB=%g \n", dt, nstep, medium,bond,ang,nonB, cdamp, cdampB, cdampAng, cdampNB  );
        return cdamp;
    }

    void set( int nstep_, double medium_, double bond_, double ang_, double nonB_, double dRcut1_, double dRcut2_ ){
        nstep   = nstep_;
        //medium  = fmax( medium_,0 );
        medium  = medium_; // we want also acceleration
        bond    = fmax( bond_  ,0 );
        ang     = fmax( ang_   ,0 );
        nonB    = fmax( nonB_  ,0 );
        dRcut1  = dRcut1_;
        dRcut2  = dRcut2_;
        bBond = bond>0;
        bAng  = ang >0;
        bNonB = nonB>0;
    }

    void setup_accel(int nstep_acc_min_, double cos_vf_acc_ ){
        nstep_acc_min = nstep_acc_min_;
        cos_vf_acc    = cos_vf_acc_;
        printf( "CollisionDamping::setup_accel nstep_acc_min %i cos_vf_acc %g \n", nstep_acc_min, cos_vf_acc );
    }

    inline double update_acceleration( Vec3d cvf ){
        double cos_fv = cvf.x/sqrt(cvf.z*cvf.y+1e-32);
        if(cvf.x<0){ 
            nstep_acc = 0;
            //printf( "MolWorld_sp3::run_no_omp(itr=%i).cleanVelocity() \n", itr, ffl.cvf.x/sqrt(ffl.cvf.z*ffl.cvf.y+1e-32), sqrt(ffl.cvf.y), sqrt(ffl.cvf.z) ); 
        }else {
            if( cos_fv>cos_vf_acc ){
                nstep_acc++;
                //printf( "MolWorld_sp3::run_no_omp(itr=%i) Acceleration: cdamp=%g(%g) nstepOK=%i cos(f,v)\n", itr, cdamp, ffl.colDamp.nstepOK, ffl.cvf.x/sqrt(ffl.cvf.z*ffl.cvf.y+1e-32) );
            }else{
                //double vf_damp = (0.5-cos_fv)*0.0;
                //double renorm_vf  = sqrt( ffl.cvf.y / ( ffl.cvf.z  + 1e-32 ) );
                //for(int i=0; i<ffl.nvecs; i++){ ffl.vapos[i] = ffl.vapos[i]*(1-vf_damp) + ffl.fapos[i]*( renorm_vf * vf_damp ); }
                nstep_acc=0;
            }
        }
        return cos_fv;
    }

};

class ForceField: public Atoms{ public:

    // ---  inherited from Atoms
    //int     natoms =0; // from Atoms
    //int    *atypes =0; // from Atoms
    //Vec3d  *apos   =0; // from Atoms

    Vec3d  *fapos __attribute__((aligned(64))) = 0; // forces on atomic positions
    Vec3d * vapos __attribute__((aligned(64))) = 0; // [natoms] velocities of atoms

    Vec3d cvf=Vec3dZero;   // <f|f>, <v|v>, <v|f> 

    Mat3d   lvec __attribute__((aligned(64)));  // lattice vectors
    Vec3i   nPBC;  // number of periodic images in each direction 
    bool    bPBC=false; // periodic boundary conditions ?
    int    npbc   =0;  // total number of periodic images
    Vec3d* shifts __attribute__((aligned(64))) =0;  // array of bond vectors shifts in periodic boundary conditions

    bool  bNonBonded            = true;  // if true we calculate non-bonded interactions
    bool  bNonBondNeighs        = false;
    bool  bSubtractBondNonBond  = true;  // if true we subtract bond energy from non-bonded energy
    bool  bSubtractAngleNonBond = true;  // if true we subtract angle energy from non-bonded energy
    bool  bClampNonBonded       = true;  // if true we clamp non-bonded energy to zero
    double  FmaxNonBonded       = 10.0;  // if bClampNonBonded>0 then we clamp non-bonded forces to this value
    std::function<double(int,const Vec3d,Vec3d&)> atomForceFunc = nullptr;
    double time_per_iter = 0;

    CollisionDamping colDamp;

    // ==================== Functions

    inline void setNonBondStrategy( int imode=0 ){
        if      ( imode>0 ){ bNonBondNeighs = true;  }
        else if ( imode<0 ){ bNonBondNeighs = false; }
        if( !bNonBonded ){
            bSubtractAngleNonBond = false;
            bSubtractBondNonBond  = false;
        }else if(bNonBondNeighs){
            bSubtractAngleNonBond = true;
            bSubtractBondNonBond  = false;
            bClampNonBonded       = false;
        }else{ 
            bSubtractAngleNonBond = true;
            bSubtractBondNonBond  = true;
            bClampNonBonded       = true;
        }
        //printf( "ForceField::setNonBondStrategy() imode=%i bNonBonded=%i bNonBondNeighs=%i bSubtractBondNonBond=%i bSubtractAngleNonBond=%i bClampNonBonded=%i\n", imode, bNonBonded, bNonBondNeighs, bSubtractBondNonBond, bSubtractAngleNonBond, bClampNonBonded );
    }

    // update atom positions using molecular dynamics (damped leap-frog)
    __attribute__((hot))  
    inline Vec3d move_atom_MD( int i, const double dt, const double Flim, const double cdamp=0.9 ){
        Vec3d f = fapos[i];
        Vec3d v = vapos[i];
        Vec3d p = apos [i];
        const Vec3d  cvf{ v.dot(f), v.norm2(), f.norm2() };
        if(cvf.z>(Flim*Flim)){ f.mul(Flim/sqrt(cvf.z)); };
        v.mul    ( cdamp );
        v.add_mul( f, dt );                         //b|=ckeckNaN_d( 1,3, (double*)&v, "v.2" );
        p.add_mul( v, dt );
        apos [i] = p;
        vapos[i] = v;
        //fapos[i] = Vec3dZero;
        return cvf;
    }

    __attribute__((hot))  
    inline Vec3d move_atom_Langevin( int i, const float dt, const double Flim,  const double gamma_damp=0.1, double T=300 ){
        Vec3d f = fapos[i];
        Vec3d v = vapos[i];
        Vec3d p = apos [i];
        const Vec3d  cvf{ v.dot(f), v.norm2(), f.norm2() };
        // ----- Langevin
        f.add_mul( v, -gamma_damp );  // check the untis  ... cdamp/dt = gamma
        // --- generate randomg force from normal distribution (i.e. gaussian white noise)
        // -- this is too costly
        //Vec3d rnd = {-6.,-6.,-6.};
        //for(int i=0; i<12; i++){ rnd.add( randf(), randf(), randf() ); }    // ToDo: optimize this
        // -- keep uniform distribution for now
        Vec3d rnd = {randf(-1.0,1.0),randf(-1.0,1.0),randf(-1.0,1.0)};
        //if(i==0){ printf( "dt=%g[arb.]  dt=%g[fs]\n", dt, dt*10.180505710774743  ); }
        f.add_mul( rnd, sqrt( 2*const_kB*T*gamma_damp/dt ) );
        v.add_mul( f, dt );                         
        p.add_mul( v, dt );
        apos [i] = p;
        vapos[i] = v;
        return cvf;
    }

    // update atom positions using FIRE (Fast Inertial Relaxation Engine)
    __attribute__((hot))  
    inline double move_atom_FIRE( int i, const double dt, const double Flim, const double cv, const double cf ){
        Vec3d  f   = fapos[i];
        Vec3d  v   = vapos[i];
        Vec3d  p   = apos [i];
        double fr2 = f.norm2();
         if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
        v.mul(             cv );
        v.add_mul( f, dt + cf );
        p.add_mul( v, dt );
        apos [i] = p;
        vapos[i] = v;
        //fapos[i] = Vec3dZero;
        return fr2;
    }

    // update atom positions using modified FIRE algorithm
    __attribute__((hot))  
    inline double move_atom_kvaziFIRE( int i, const double dt, const double Flim ){
        Vec3d  f   = fapos[i];
        Vec3d        v   = vapos[i];
        Vec3d        p   = apos [i];
        double fr2 = f.norm2();
        if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
        double vv  = v.norm2();
        double ff  = f.norm2();
        double vf  = v.dot(f);
        double c   = vf/sqrt( vv*ff + 1.e-16    );
        double v_f =    sqrt( vv/( ff + 1.e-8)  );
        double cv;
        double cf = cos_damp_lin( c, cv, 0.01, -0.7,0.0 );
        v.mul(                 cv );
        v.add_mul( f, dt + v_f*cf );
        p.add_mul( v, dt );
        apos [i] = p;
        vapos[i] = v;
        //fapos[i] = Vec3dZero;
        return fr2;
    }

    // update atom positions using gradient descent
    __attribute__((hot))  
    inline double move_atom_GD(int i, double dt, double Flim){
        Vec3d  f   = fapos[i];
        Vec3d  p = apos [i];
        double fr2 = f.norm2();
        if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
        p.add_mul( f, dt );
        apos [i] = p;
        //fapos[i] = Vec3dZero;
        return fr2;
    }

    // update atom positions using gradient descent
    __attribute__((hot))  
    double move_GD(float dt, double Flim=100.0 ){
        double F2sum=0;
        for(int i=0; i<natoms; i++){
            F2sum += move_atom_GD(i, dt, Flim);
        }
        return F2sum;
    }

    // update atom positions using gradient descent
    __attribute__((hot))  
    double move_MD(float dt, const double Flim=100.0, const double cdamp=0.9 ){
        cvf = Vec3dZero;
        for(int i=0; i<natoms; i++){
            cvf.add( move_atom_MD( i, dt,Flim,cdamp ) );
        }
        return cvf.z;
    }

    void cleanForce(){ 
        Energy=0;
        for(int i=0; i<natoms; i++){ fapos[i]=Vec3dZero;  } 
    }
    void cleanVelocity(){ 
        for(int i=0; i<natoms; i++){ vapos[i]=Vec3dZero;  } 
    }

};

#endif

