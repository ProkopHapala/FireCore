
#ifndef  ProjectiveDynamics_d_h
#define  ProjectiveDynamics_d_h

/// @brief Projective Dynamics is fast end robust solver of dynamics of deformable objects consisitng of mass-points connected by very stiff springs. It is position-based solver based on implicit integration avoiding numerical instability due to large timesteps by wich suffer explicit force-based dynamical integration schemes. 
/// see: * [Projective Dynamics: Fusing Constraint Projections for Fast Simulation](https://www.projectivedynamics.org/Projective_Dynamics/index.html)
//       * [Parallel iterative solvers for real-time elastic deformations](https://mfratarcangeli.github.io/publication/sa2018course/)
//       * [Position Based Dynamics](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)
//       * https://github.com/ProkopHapala/FireCore/wiki/Constrain-Solvers
//       

#include "datatypes.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  
#include "VecN.h"

//#include "Buckets.h"
//#include "raytrace.h"
//#include "geom3D.h"
//#include "Interfaces.h"

#include "arrayAlgs.h"

#include "SparseMatrix.h"
#include "SparseMatrix2.h"


double checkDist(int n, const Vec3d* vec, const Vec3d* ref, int verb=1, double tol=1e-12 );
void print_vector(int n, double * a, int pitch, int j0, int j1 );

// ========== Helper classes ==========

struct EdgeVertBond{ 
    Vec3i verts; 
    double c; 
    double K; 
    Vec3d f=Vec3dZero; 
};

struct SmartMixer{
    float b_end   = 0.75f;
    float b_start = 0.55f;
    int   istart = 3;
    int   iend   = 10; 

    // float get_bmix(int itr){
    //     if      ( itr < istart ) { return 0; }
    //     else if ( itr > iend   ) { return b_end; }
    //     else                     { return b_start + (b_end - b_start ) * (itr - istart) / (iend - istart); }
    // }

    inline float get_bmix(int itr){
        if ( itr < istart ) { return 0; }
        else                { return b_end; }
    }
};

// ========== main class ==========

class ProjectiveDynamics_d { public:
    double time=0;
    double Cdrag   = 0.0;
    //Vec3d Gravity = {0,-9.81,0};
    // Rotating frame
    Vec3d pos0{0.,0.,0.};
    Vec3d ax{0.0,0.0,1.0};
    //double omega = 0.05;
    //Quat4d accel{ 0.0,-9.81 , 0.0 , 0.0 };    // acceleration
    Quat4d accel{ 0.0, 0.0 , 0.0 , 0.0  };    // acceleration
    Quat4d rot0 { 0.0, 0.0 , 0.0 , 0.0  };    // center of rotation
    Quat4d omega{ 0.0, 0.0 , 0.0 , 0.05 };    // angular velocity

    //double dt      = 2e-3; //double kGlobal = 1e+6;
    //double dt      = 2e-3;    double kGlobal = 1e+7;
    double dt      = 0.5e-3;  double kGlobal = 1e+8;
    //double damping = 1e-4;
    double damping  = 0.05;
    int nSolverIters = 10;

    //std::vector<int> damped_bonds;
    //std::unordered_set<int> damped_points;


    int nPoint=0, nNeighMax=0, nNeighTot=0;
    // cpu buffers
    Quat4d* points=0;    // position and mass
    Quat4d* forces=0;    // force and energy
    Quat4d* vel   =0;    // velocity
    Quat4d* vel0 = 0;
    //Quat4d* impuls=0;  // accumulated impulse from corrector
    Quat4d* bvec  =0;    // right hand side of linear system Ap=b ( it has meaning of internal force due to projective dybamics matrix A = K + D + I ) 
    Vec3d * bvec0 =0;
    
    // Binding to external linear solver of system Ap=b ( where A is matrix of projective dynamics )
    Quat4f* extern_b = 0;
    Quat4f* extern_x = 0;
    void (*extern_solve)();   // function pointer extern_solve

    double* kFix=0;   // force constant for fixed points
    double  kLinRegularize = 1.0;

    Quat4d* params=0;    // neighbor parameters (l0,kP,kT,damp)
    int*    neighs=0;    // neighbor indices
    Vec2i*  neighBs=0;   // neighbor bond indices
    int*    neighB2s=0;  // neighbor indices

    // Cholesky / Projective Dynamics 
    double* PDmat=0;
    double* LDLT_L=0;
    double* LDLT_D=0; 
    int*    neighsLDLT=0;
    int     nNeighMaxLDLT=0;
    SparseMatrix2<double>  Lsparse; 
    SparseMatrix2<double>  LsparseT; 
    SparseMatrix<double>   PDsparse;
    
    SmartMixer mixer;

    // choice of linear solver method
    int linSolveMethod = 2; // 0=CG, 1=CGsparse, 2=Cholesky
    enum class LinSolveMethod{ Cholesky=2, CholeskySparse=3, Jacobi=4, GaussSeidel=5, JacobiMomentum=6, GSMomentum=7, JacobiFlyMomentum=8, GSFlyMomentum=9, Force=10, JacobiDiff=11, MomentumDiff=12, ExternDiff=13 };
    bool   bApplyResudualForce = true;
    double residualForceFactor = 1.0;

    Vec3d*  ps_cor      =0; // new Vec3d[nPoint];
    Vec3d*  ps_pred     =0; // new Vec3d[nPoint]; 
    Vec3d*  linsolve_b  =0; // new Vec3d[nPoint];
    Vec3d*  linsolve_yy =0; // new Vec3d[nPoint];
    Vec3d*  ps_0        =0;

    int     nBonds  =0;    // number of bonds
    Quat4d* bparams =0;    // bond parameters (l0,kPress,kTens,damp)
    Vec2i*  bonds   =0;    // indices of bonded points (i,j)
    double* strain =0;    // strain
    Vec2d*  maxStrain=0;

    // ====== Invairiants

    double mass = 0;
    Vec3d cog  = Vec3dZero;
    Vec3d vcog = Vec3dZero;
    Mat3d I    = Mat3dZero;
    Vec3d L    = Vec3dZero;
    Vec3d torq = Vec3dZero;

    // callback function pointer what to do in between iterations
    void (*user_update)(double dt);

    // ============ inline  Functions

    inline void set_time_step( float dt_ ){
        dt      = dt_;
        //inv_dt2 = 1.0f / (dt * dt);
    }

    inline Vec3d getPointForce( int i ){
        return (vel[i].f*Cdrag) + (accel.f*points[i].w);
    }

    inline void reallocFixed(){ _realloc0( kFix, nPoint, 0.0 ); }
    inline void cleanForce(){ for (int i=0; i<nPoint; i++){ forces[i]=Quat4dZero; } };
    inline void cleanVel  ( Quat4d v0=Quat4dZero ){ for (int i=0; i<nPoint; i++){ vel[i]=v0; } };

    inline double getLinearBondStiffness( int ib ){
        // (l0,kPress,kTens,damp)
        //return  bparams[ib].y;  //k = kPress;
        return  bparams[ib].z;  //k = kTens;  
    }

    // ============ Functions
 
void recalloc(int nPoint_, int nNeighMax_, int nBonds_=0);
void realloc_LinearSystem( bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true );
double norm_butFixed(Vec3d* ps );
double norm_butFixed(Quat4d* ps );
void make_PD_Matrix(double* A, double dt );
void make_PDmat_sparse(SparseMatrix<double>& A, double dt, bool bRealloc );
void rhs_ProjectiveDynamics(const Vec3d* pnew, Vec3d* b);
Vec3d rhs_ProjectiveDynamics_i(int i, const Vec3d* pnew);
void rhs_ProjectiveDynamics_(const Vec3d* pnew, Vec3d* b);
void updatePD_RHS(const Vec3d* pnew, Quat4d* bvec );
void dotPD(const Vec3d* p, Vec3d* f );
void updatePD_dRHS(const Vec3d* pnew, Quat4d* bvec );
void updateJacobi_lin(Vec3d* ps_in, Vec3d* ps_out, Quat4d* bvec, Vec3d* rs=0 );
void updateGaussSeidel_lin(Vec3d* ps, Quat4d* bvec );
void evalTrussForce(Vec3d* ps_in, Vec3d* force );
void updateJacobi_fly(Vec3d* ps_in, Vec3d* ps_out );
void updateGaussSeidel_fly(Vec3d* ps );
void updateIterativeMomentum(Vec3d* psa, Vec3d* psb );
void updateIterativeJacobi(Vec3d* psa, Vec3d* psb );
void updateIterativeJacobiDiff(Vec3d* psa, Vec3d* psb );
void updateIterativeMomentumDiff(Vec3d* psa, Vec3d* psb );
void updateIterativeExternDiff(Vec3d* psa, Vec3d* psb );
void prepare_LinearSystem(bool bRealloc=true, bool bCholesky=true, int nNeighMaxLDLT_=32, bool bDens=true );
void run_LinSolve(int niter);
void updateInveriants(bool bPrint=false);



double evalBondTension();
void printNeighs(int i);
void printAllNeighs();
double getFmax();
void setFixPoints(int n, int* fixPoints, double Kfix=1e12, bool bRealloc=true );

//void setOpt(double dt_, double damp_ );


};   // ProjectiveDynamics_d

#endif
