
#ifndef MMFFsp3_loc_h
#define MMFFsp3_loc_h

#define ANG_HALF_COS  1

#include <omp.h>

#include "fastmath.h"   // fast math operations
#include "Vec2.h"       // 2D vector
#include "Vec3.h"       // 3D vector
#include "quaternion.h" // quaternions

#include "constants.h"   
#include "Forces.h"     // various physical interactions
#include "SMat3.h"             // Symmetric Matrix
#include "molecular_utils.h"   // various molecular utilities

#include "NBFF.h" // Non-Bonded Force Field

//#include <cstdio>

// ========================
// ====   MMFFsp3_loc  ====
// ========================
/*
    This is a local version of MMFFsp3, i.e. we store all parameters per atom, and we compute energy and forces for each atom separately. 
    The recoil forces on neighbors are stored in a temporary array fneigh. Which is latter assembled into the global force array fapos.
    This allows for efficient parallelization, since we avoid synchronization of the global force array fapos or using atomic operations.
    The drawback is that we need to allocate additional memory for fneigh and possibly lower the performance due to cache misses.
*/

//class MMFFsp3_loc: public NBFF { public:

class MMFFsp3_loc : public NBFF { public:
    static constexpr const int nneigh_max = 4; // maximum number of neighbors

    // === inherited from NBFF
    // int natoms=0;        // [natoms] // from Atoms
    //int   * atypes  =0; // [natom]  atom types
    //Vec3d * apos  =0;     // [natoms] // from Atoms
    //Vec3d * vapos = 0;    // [natom]  velocities of atoms
    //Vec3d * fapos  =0;    // [natoms] // from NBFF
    //Quat4i* neighs =0;    // [natoms] // from NBFF
    //Quat4i* neighCell=0;  // [natoms] // from NBFF
    //Quat4d* REQs =0;      // [nnode]  // from NBFF
    //bool    bPBC=false;       // from NBFF
    //Vec3i   nPBC;         // from NBFF 
    //Mat3d   lvec;         // from NBFF
    //double  Rdamp  = 1.0; // from NBFF

    // dimensions of the system
    int  nDOFs=0,nnode=0,ncap=0,nvecs=0,ntors=0;
    double Etot,Eb,Ea, Eps,EppT,EppI;  // total energy, bond energy, angle energy, pi-sigma energy, pi-pi torsion energy, pi-pi interaction energy

    double *  DOFs = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    bool doBonds  =true; // compute bonds
    bool doNeighs =true; // compute neighbors
    bool doPiPiI  =true; // compute pi-pi interaction parallel
    bool doPiPiT  =true; // compute pi-pi torsion     perpendicular
    bool doPiSigma=true; // compute pi-sigma interaction
    bool doAngles =true; // compute angles
    //bool doEpi    =true; // compute pi-electron interaction

    // bool    bCollisionDamping        = false; // if true we use collision damping
    // bool    bCollisionDampingAng     = false;
    // bool    bCollisionDampingNonBond = false;  // if true we use collision damping for non-bonded interactions
    // double  damping_medium           = 1.0;   // cdamp       = 1 -(damping_medium     /ndampstep     )
    // double  collisionDamping         = 0.0;   // col_damp    =     collisionDamping   /(dt*ndampstep )
    // double  collisionDamping_ang     = 0.0;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    // double  collisionDamping_NB      = 0.0;   // col_damp_NB =     collisionDamping_NB/(dt*ndampstep )
    // int     ndampstep                = 10;    // how many steps it takes to decay velocity to to 1/e of the initial value
    // double  col_damp_dRcut1          = -0.2;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]
    // double  col_damp_dRcut2          =  0.3;   // non-covalent collision damping interaction goes between 1.0 to 0.0 on interval  [ Rvdw , Rvdw+col_damp_dRcut ]
    // double col_damp      = 0.0;  //  collisionDamping   /(dt*ndampstep );
    // double col_damp_NB   = 0.0;  //  collisionDamping_NB/(dt*ndampstep );
    // double col_dampAng   = 0.0;  //  collisionDamping   /(dt*ndampstep );

    CollisionDamping colDamp;

    Vec3d cvf=Vec3dZero; // <f|f>, <v|v>, <v|f> 

    bool bEachAngle = false; // if true we compute angle energy for each angle separately, otherwise we use common parameters for all angles
    bool bTorsion   = false; // if true we compute torsion energy
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 }; // default parameters for neighbors, c0 is cosine of equilibrium angle, Kss is bond stiffness, Ksp is pi-sigma stiffness, c0_e is cos of equilibrium angle for pi-electron interaction

    // Dynamical Varaibles;
    //Vec3d *   apos=0;   // [natom]
    //Vec3d *  fapos=0;   // [natom]
    Vec3d *  pipos=0;   // [nnode]
    Vec3d * fpipos=0;   // [nnode]

    // Aux Dynamil
    Vec3d * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Vec3d * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    //Quat4i*  neighs =0;   // [natoms] // from NBFF
    Quat4i*  bkneighs=0;   // [natoms]  inverse neighbors
    
    Quat4d*  apars __attribute__((aligned(64))) =0;  // [nnode] angle parameters
    Quat4d*  bLs   __attribute__((aligned(64))) =0;  // [nnode] bond lengths
    Quat4d*  bKs   __attribute__((aligned(64))) =0;  // [nnode] bond stiffness
    Quat4d*  Ksp   __attribute__((aligned(64))) =0;  // [nnode] stiffness of pi-alignment
    Quat4d*  Kpp   __attribute__((aligned(64))) =0;  // [nnode] stiffness of pi-planarization

    Vec3d*   angles=0; // [nnode*6]  angles between bonds

    Quat4i*  tors2atom  __attribute__((aligned(64))) =0; // [ntors]  torsion atoms
    Quat4d*  torsParams __attribute__((aligned(64))) =0; // [ntors]  torsion parameters

    Quat4d* constr  __attribute__((aligned(64))) = 0; // [natom]  constraints
    Quat4d* constrK __attribute__((aligned(64))) = 0; // [natom]  constraints
    //Vec3d * vapos  __attribute__((aligned(64))) = 0; // [natom]  velocities of atoms

    Mat3d   invLvec; // inverse lattice vectors

    bool    bAngleCosHalf         = true;   // if true we use evalAngleCosHalf() instead of evalAngleCos() to compute anglular energy
    // these are defined in ForceFiled.h
    //bool    bSubtractAngleNonBond = false;  // if true we subtract angle energy from non-bonded energy
    //bool    bSubtractBondNonBond  = false;  // if true we subtract bond energy from non-bonded energy

    //int itr_DBG=0;

// =========================== Functions

// reallcoate MMFFsp3_loc
void realloc( int nnode_, int ncap_, int ntors_=0 ){
    nnode=nnode_; ncap=ncap_; ntors=ntors_;
    natoms= nnode + ncap; 
    nvecs = natoms+nnode;  // each atom as also pi-orientiation (like up-vector)
    nDOFs = nvecs*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    int ipi0=natoms;
    
    _realloc0(  DOFs    , nDOFs , (double)NAN );
    _realloc0( fDOFs    , nDOFs , (double)NAN );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc0( fneigh  , nnode*4, Vec3dNAN );
    _realloc0( fneighpi, nnode*4, Vec3dNAN );
    // ----- Params [natom]
    _realloc0( atypes    , natoms, -1 );
    _realloc0( neighs    , natoms, Quat4iMinusOnes );
    _realloc0( neighCell , natoms, Quat4iMinusOnes );
    _realloc0( bkneighs  , natoms, Quat4iMinusOnes);
    _realloc0( apars     , nnode, Quat4dNAN );
    _realloc0( bLs       , nnode, Quat4dNAN );
    _realloc0( bKs       , nnode, Quat4dNAN );
    _realloc0( Ksp       , nnode, Quat4dNAN );
    _realloc0( Kpp       , nnode, Quat4dNAN );

    // Additional:
    // Angles
    _realloc0( angles, nnode*6, Vec3dNAN );   // 6=4*3/2
    // Torsions
    _realloc0( tors2atom,  ntors, Quat4iZero );
    _realloc0( torsParams, ntors, Quat4dNAN  ); 

    _realloc0( constr    , natoms, Quat4dOnes*-1. );
    _realloc0( constrK   , natoms, Quat4dOnes*-1. );
}

// clone from another MMFFsp3_loc
void clone( MMFFsp3_loc& from, bool bRealloc, bool bREQsDeep=true ){
    realloc( from.nnode, from.ncap  );
    lvec   =from.lvec;
    invLvec=from.invLvec;
    for(int i=0; i<nDOFs; i++){
        DOFs[i]=from. DOFs[i];
     //fDOFs[i]=from.fDOFs[i];
    }
    for(int i=0; i<natoms; i++){
        atypes   [i]=from.atypes   [i];
        neighs   [i]=from.neighs   [i];
        neighCell[i]=from.neighCell[i];
        bkneighs [i]=from.bkneighs [i];
        constr   [i]=from.constr   [i];
    }
    for(int i=0; i<nnode; i++){
        apars[i]=from.apars[i];
        bLs  [i]=from.bLs  [i];    
        bKs  [i]=from.bKs  [i];    
        Ksp  [i]=from.Ksp  [i]; 
        Kpp  [i]=from.Kpp  [i];    
    }
    if(from.REQs){
        if(bREQsDeep){ _realloc(REQs, natoms ); for(int i=0; i<natoms; i++){  REQs[i]=from.REQs[i]; } }
        else         {          REQs=from.REQs;                                                       }
    }
}

// deallcoate MMFFsp3_loc
void dealloc(){
    _dealloc(DOFs );
    _dealloc(fDOFs);
    apos   = 0;
    fapos  = 0;
    pipos  = 0;
    fpipos = 0;

    //vDOFs  = 0;   // to-do (not sure if we should deallocate it here ?)
    vapos  = 0;   // to-do (not sure if we should deallocate it here ?)
    //vpipos = 0;   // to-do (not sure if we should deallocate it here ?)
    _dealloc(atypes);
    _dealloc(neighs);
    _dealloc(neighCell);
    _dealloc(bkneighs);
    _dealloc(apars);
    _dealloc(bLs);
    _dealloc(bKs);
    _dealloc(Ksp);
    _dealloc(Kpp);
    _dealloc(angles);

    _dealloc(tors2atom  );
    _dealloc(torsParams );
    nnode=0; ncap=0; ntors=0; natoms=0; nvecs =0; nDOFs =0; int ipi0=natoms;
}

// set lattice vectors
void setLvec(const Mat3d& lvec_){ lvec=lvec_; lvec.invert_T_to( invLvec ); }

// find optimal time-step for dy FIRE optimization algorithm 
double optimalTimeStep(double m=1.0){
    double Kmax = 1.0;
    for(int i=0; i<nnode; i++){ 
        Kmax=fmax(Kmax, bKs[i].x ); 
        Kmax=fmax(Kmax, bKs[i].y ); 
        Kmax=fmax(Kmax, bKs[i].z ); 
        Kmax=fmax(Kmax, bKs[i].w ); 
    }
    return M_PI*2.0*sqrt(m/Kmax)/10.0;  // dt=T/10;   T = 2*pi/omega = 2*pi*sqrt(m/k)
}

// ============== Evaluation

// evaluate energy and forces for single atom (ia) depending on its neighbors
__attribute__((hot))   
double eval_atom(const int ia){
    //printf( "MMFFsp3_loc::eval_atom(%i) bSubtractBondNonBond=%i \n", ia, bSubtractBondNonBond );
    double E=0;
    const double Fmax2  = FmaxNonBonded*FmaxNonBonded;
    const double R2damp = Rdamp*Rdamp;
    const Vec3d pa  = apos [ia]; 
    const Vec3d hpi = pipos[ia]; 

    //printf( "apos[%i](%g,%g,%g)\n",ia,apos[ia].x,apos[ia].y,apos[ia].z );
    //return E;

    //Vec3d& fa  = fapos [ia]; 
    //Vec3d& fpi = fpipos[ia];
    Vec3d fa   = Vec3dZero;
    Vec3d fpi  = Vec3dZero; 
    
    //--- array aliases
    const int*    ings = neighs   [ia].array; // neighbors
    const int*    ingC = neighCell[ia].array; // neighbors cell index
    const double* bK   = bKs      [ia].array; // bond stiffness
    const double* bL   = bLs      [ia].array; // bond length
    const double* Kspi = Ksp      [ia].array; // pi-sigma stiffness
    const double* Kppi = Kpp      [ia].array; // pi-pi stiffness
    Vec3d* fbs  = fneigh   +ia*4;             // forces on bonds
    Vec3d* fps  = fneighpi +ia*4;             // forces on pi vectors

    const bool bColDampB   = colDamp.bBond && vapos; // if true we use collision damping
    const bool bColDampAng = colDamp.bAng  && vapos; // if true we use collision damping for non-bonded interactions

    // // --- settings
    // double  ssC0 = apars[ia].x;
    // double  ssK  = apars[ia].y;
    // double  piC0 = apars[ia].z;
    // //bool    bPi  = ings[3]<0;   we distinguish this by Ksp, otherwise it would be difficult for electron pairs e.g. (-O-C=)

    const Quat4d& apar  = apars[ia]; // [c0, Kss, Ksp, c0_e] c0 is cos of equilibrium angle, Kss is bond stiffness, Ksp is pi-sigma stiffness, c0_e is cos of equilibrium angle for pi-electron interaction
    const double  piC0 = apar.w;     // cos of equilibrium angle for pi-electron interaction

    //printf( "ang0 %g cs0(%g,%g)\n", atan2(cs0_ss.y,cs0_ss.x)*180/M_PI, cs0_ss.x,cs0_ss.x );

    //--- Aux Variables 
    Quat4d  hs[4]; // bond vectors (normalized in .xyz ) and their inverse length in .w
    Vec3d   f1,f2; // working forces
    
    //bool bErr=0;
    //const int ia_DBG = 0;
    //if(ia==ia_DBG)printf( "ffl[%i] neighs(%i,%i,%i,%i) \n", ia, ings[0],ings[1],ings[2],ings[3] );
    //printf("lvec %i %i \n", ia, ia_DBG ); // printMat(lvec);

    //if( ia==5 ){ printf( "ffls[%2i] atom[%2i] ng(%3i,%3i,%3i,%3i) ngC(%3i,%3i,%3i,%3i) shifts=%li bPBC=%i\n", id, ia,   ings[0],ings[1],ings[2],ings[3],   ingC[0],ingC[1],ingC[2],ingC[3], (long)shifts, bPBC ); };

    for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break in the loop

    // double Eb=0;
    // double Ea=0;
    // double EppI=0;
    // double Eps=0;

    // --------- Bonds Step
    for(int i=0; i<4; i++){ // loop over bonds
        int ing = ings[i];
        //printf( "bond[%i|%i=%i]\n", ia,i,ing );
        //fbs[i]=Vec3dOne; fps[i]=Vec3dOne;
        //fbs[i]=Vec3dZero; fps[i]=Vec3dZero; // NOTE: wee need to initialize it before, because of the break
        if(ing<0) break;

        //printf("ia %i ing %i \n", ia, ing ); 
        Vec3d  pi = apos[ing];
        Quat4d h; 
        h.f.set_sub( pi, pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        
        //Vec3d h_bak = h.f;    
        //shifts=0; 
        // Periodic Boundary Conditions
        if(bPBC) [[likely]]  {   
            if(shifts) [[likely]] {   // if we have bond shifts vectors we use them
                int ipbc = ingC[i]; 
                //Vec3d sh = shifts[ipbc]; //apbc[i]  = pi + sh;
                h.f.add( shifts[ipbc] );
                //if( (ia==0) ){ printf("ffls[%i] atom[%i,%i=%i] ipbc %i shifts(%g,%g,%g)\n", id, ia,i,ing, ipbc, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z); };
                //if( (ipbc!=4) ){ printf("ffls[%i] atom[%i,%i=%i] ipbc %i shifts(%g,%g,%g)\n", id, ia,i,ing, ipbc, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z); };
            }else{  // if we don't have bond shifts vectors we use lattice vectors
                Vec3i g  = invLvec.nearestCell( h.f );
                // if(ia==ia_DBG){
                //     Vec3d u; invLvec.dot_to(h.f,u);
                //     printf( "CPU:bond[%i,%i] u(%6.3f,%6.3f,%6.3f) shi(%6.3f,%6.3f,%6.3f) \n", ia, ing, u.x,u.y,u.z,   (float)g.x,(float)g.y,(float)g.z );
                // }
                Vec3d sh = lvec.a*g.x + lvec.b*g.y + lvec.c*g.z;
                h.f.add( sh );
                //apbc[i] = pi + sh;
            }
        }

        //wrapBondVec( h.f );
        //printf( "h[%i,%i] r_old %g r_new %g \n", ia, ing, h_bak.norm(), h.f.norm() );
       
        // initial bond vectors
        const double l = h.f.normalize(); 
        h.e    = 1/l;
        hs [i] = h;

        //bErr|=ckeckNaN( 1,4, (double*)(hs+i), [&]{ printf("atom[%i]hs[%i]",ia,i); } );
        // bond length force
        //continue; 

        //if(ia==ia_DBG) printf( "ffl:h[%i|%i=%i] l %g h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", ia,i,ing, l, h.x,h.y,h.z, apos[ing].x,apos[ing].y,apos[ing].z, pa.x,pa.y,pa.z );

        if(ia<ing){     // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            if(doBonds)[[likely]]{  

                if(bColDampB){ // 
                    //printf( "MMFFsp3_loc::eval_atom() bCollisionDamping=%i col_damp=%g \n", bCollisionDamping, col_damp ); exit(0);
                    // double invL = 1./l;
                    // double dv   = d.dot( vel[b.y].f - vel[b.x].f )*invL;
                    // double mcog = pj.w + pi.w;
                    // double reduced_mass = pj.w*pi.w/mcog;
                    // double imp  = collision_damping * reduced_mass * dv;
                    // for masses = 1.0 we have reduced_mass = 1*1/(1+1) = 0.5

                    double dv   = h.f.dot( vapos[ing] - vapos[ia] );
                    double fcol = colDamp.cdampB * 0.5 * dv;   // col_damp ~ collision_damping/(dt*ndampstep);     f = m*a = m*dv/dt

                    f1.set_mul( h.f, fcol );
                    fbs[i].sub(f1);  fa.add(f1); 

                    // double w    = damp_rate * 0.5 * smoothstep_down(sqrt(r2),R, R+dRcut ); // ToDo : we can optimize this by using some other cutoff function which depends only on r2 (no sqrt)
                    // //double w    = damp_rate * 0.5 * R8down         (r2,      R, R+dRcut );
                    // double fcol = w * d.dot( vapos[j]-vi );                                // collisionDamping ~ 1/(dt*ndampstep);     f = m*a = m*dv/dt
                    // Vec3d fij; fij.set_mul( d, fcol/r2 ); //  vII = d*d.fot(v)/|d|^2 
                    // fx+=fij.x;
                    // fy+=fij.y;
                    // fz+=fij.z;
                }

                // bond length force
                //printf( "bond[%i,%i] l=%g bL=%g bK=%g f=%g \n", ia, ing, l, bL[i], bK[i], f1 );
                E+= evalBond( h.f, l-bL[i], bK[i], f1 ); 
                // { 
                //     f1=Vec3dZero; // DEBUG
                // } 

                if(bSubtractBondNonBond) [[likely]] { // subtract non-bonded interactions between atoms which have common neighbor
                    Vec3d fij=Vec3dZero;
                    //Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                    Quat4d REQij = _mixREQ(REQs[ia],REQs[ing]);  // combine van der Waals parameters for the pair of atoms
                    Vec3d dp = h.f*l;
                    E -= getLJQH( dp, fij, REQij, R2damp ); // subtract non-bonded interactions 
                    if(bClampNonBonded)[[likely]] { clampForce( fij, Fmax2 ); }
                    //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                    //bErr|=ckeckNaN( 1,3, (double*)&fij, [&]{ printf("atom[%i]fLJ2[%i,%i]",ia,i,j); } );
                    f1.sub(fij);
                    //printf( "ffl:SubtractBondNonBond[%i|%i] r=%g fij(%g,%g,%g) REQ(%g,%g,%g)\n", ia,ing, dp.norm(), fij.x,fij.y,fij.z,  REQij.x,REQij.y,REQij.z);
                }
                
                fbs[i].sub(f1);  fa.add(f1); 

                //double Ebi = evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);
                //E +=Ebi; 
                //Eb+=Ebi; 

                //printf( "ffl:bond[%i|%i=%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,i,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
                //if(ia==ia_DBG)printf( "ffl:bond[%i|%i=%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,i,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
                //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fbond[%i]",ia,i); } );
            }

            double kpp = Kppi[i];
            if( (doPiPiI) && (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                // pi-pi interaction (make them parallel)
                E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                
                //double EppIi = evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                //E +=EppIi; 
                //EppI+=EppIi;
                
                //if(verbosity>0)printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g), EppI=%g\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z, EppIi  );
                //printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g), EppI=%g, Epp=%g\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z, EppIi, EppI  );
                //if(ia==ia_DBG)printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
                //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fpp1[%i]",ia,i); } );
                //bErr|=ckeckNaN( 1,3, (double*)&f2, [&]{ printf("atom[%i]fpp2[%i]",ia,i); } );
            }
            // ToDo: triple bonds ?
            
        } 
        
        // pi-sigma 
        //if(bPi){    
        double ksp = Kspi[i];
        if( doPiSigma && (ksp>1e-6) ){  
            // pi-sigma interaction (make them orthogonal)
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)

            //double Epsi = evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)
            //E +=Epsi; 
            //Eps+=Epsi;

            //printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z  );
            //if(kpp<-1e-6){  E += evalPiAling( hpi, pipos[ing], 1., 1., -kpp, f1, f2 );   fpi.add(f1);  fps[i].add(f2);  }    //   align pi-electron pair (e.g. in Nitrogen)
            //if(ia==ia_DBG)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z  );
            //if(ia==ia_DBG)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fsp1[%i]",ia,i); } );
            //bErr|=ckeckNaN( 1,3, (double*)&f2, [&]{ printf("atom[%i]fsp2[%i]",ia,i); } );
        }
        //}
        
    }

    //printf( "MMFF_atom[%i] cs(%6.3f,%6.3f) ang=%g [deg]\n", ia, cs0_ss.x, cs0_ss.y, atan2(cs0_ss.y,cs0_ss.x)*180./M_PI );
    
    
    // ======= Angle Step : we compute forces due to angles between bonds, and also due to angles between bonds and pi-vectors

    // {
    //     doAngles = false; // DEBUG
    // }

    if(doAngles){

    double  ssK,ssC0;
    Vec2d   cs0_ss;
    Vec3d*  angles_i;
    if(bEachAngle)[[unlikely]] {  // if we compute angle energy for each angle separately, we need to store angle parameters for each angle
        angles_i = angles+(ia*6);
    }else{          // otherwise we use common parameters for all angles
        ssK    = apar.z;
        cs0_ss = Vec2d{apar.x,apar.y};
        ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
    }

    int iang=0;
    for(int i=0; i<3; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4d& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4d& hj = hs[j];    
            if(bEachAngle) [[unlikely]] {  //
                // 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
                //  0    1    2    3    4    5 
                cs0_ss = angles_i[iang].xy();  
                ssK    = angles_i[iang].z;
                //printf( "types{%i,%i,%i} ssK %g cs0(%g,%g) \n", atypes[ing], atypes[ia], atypes[jng], ssK, cs0_ss.x, cs0_ss.y  );
                iang++; 
            };
            //bAngleCosHalf = false;
            //double Eai;
            //printf( "bAngleCosHalf= %i\n", bAngleCosHalf);
            if( bAngleCosHalf )[[likely]] { 
                E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
                //Eai = evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            }else{ 
                E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
                //Eai = evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            }

            //E +=Eai; 
            //Ea+=Eai;
            
            //printf( "ffl:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //if(ia==ia_DBG)printf( "ffl:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fss1[%i,%i]",ia,i,j); } );
            //bErr|=ckeckNaN( 1,3, (double*)&f2, [&]{ printf("atom[%i]fss2[%i,%i]",ia,i,j); } );
            fa    .sub( f1+f2  );  // apply force on the central atom

            if(bColDampAng) [[unlikely]] {
                Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                double dv   = dp.dot( vapos[ing] - vapos[ia] );
                double fcol = colDamp.cdampAng * 0.5 * dv;   // col_damp ~ collision_damping/(dt*ndampstep);     f = m*a = m*dv/dt
                dp.mul( fcol/dp.norm2() );
                f1.sub(dp);
                f2.add(dp);
            }

            if(bSubtractAngleNonBond) [[likely]] { // subtract non-bonded interactions between atoms which have common neighbor
                Vec3d fij=Vec3dZero;
                //Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);  // combine van der Waals parameters for the pair of atoms
                Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f ); // compute vector between neighbors i,j
                //Vec3d dp   = hj.f*(1./hj.w) - hi.f*(1./hi.w);
                //Vec3d dp   = apbc[j] - apbc[i];
                E -= getLJQH( dp, fij, REQij, R2damp ); // subtract non-bonded interactions 
                if(bClampNonBonded)[[likely]] { clampForce( fij, Fmax2 ); }
                //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                //bErr|=ckeckNaN( 1,3, (double*)&fij, [&]{ printf("atom[%i]fLJ2[%i,%i]",ia,i,j); } );
                f1.sub(fij); // 
                f2.add(fij);
            }
            // apply forces on neighbors
            fbs[i].add( f1     ); 
            fbs[j].add( f2     );
            //if(ia==ia_DBG)printf( "ffl:ANG[%i|%i,%i] fa(%g,%g,%g) fbs[%i](%g,%g,%g) fbs[%i](%g,%g,%g)\n", ia,ing,jng, fa.x,fa.y,fa.z, i,fbs[i].x,fbs[i].y,fbs[i].z,   j,fbs[j].x,fbs[j].y,fbs[j].z  );
            // ToDo: subtract non-covalent interactions
        }
    }
    }    
    //if(bErr){ printf("ERROR in ffl.eval_atom[%i] => Exit() \n", ia ); exit(0); }
    
    /*
    double Kfix = constr[ia].w;
    if(Kfix>0){  
        //printf( "applyConstrain(i=%i,K=%g)\n", ia, Kfix );
        Vec3d d = constr[ia].f-pa;
        d.mul( Kfix );
        double fr2 = d.norm2();
        double F2max = 1.0;
        if( fr2>F2max ){
            d.mul( sqrt(F2max/fr2) );
        }
        fa.add( d );
        E += d.norm()*Kfix*0.5;
    }
    */
    
    //fapos [ia].add(fa ); 
    //fpipos[ia].add(fpi);
    fapos [ia]=fa; 
    fpipos[ia]=fpi;

    //printf( "MYOUTPUT atom %i Ebond= %g Eangle= %g Edihed= %g Eimpr = %g Etot=%g\n", ia+1, Eb, Ea, EppI, Eps, E );
    //fprintf( file, "atom %i Ebond= %g Eangle= %g Edihed= %g Eimpr= %g Etot=%g \n", ia+1, Eb, Ea, EppI, Eps, E );
    

    return E;
}



// evaluate energy and forces for single atom (ia) depending on its neighbors
template<bool _bPBC, bool _doBonds, bool _doPiPiI, bool _doPiSigma, bool _doAngles, bool _bSubtractBondNonBond, bool _bClampNonBonded, bool _bEachAngle, bool _bAngleCosHalf, bool _bSubtractAngleNonBond>
__attribute__((hot))   double eval_atom_t(const int ia){
    double E=0;
    const double Fmax2  = FmaxNonBonded*FmaxNonBonded;
    const double R2damp = Rdamp*Rdamp;
    const Vec3d pa  = apos [ia]; 
    const Vec3d hpi = pipos[ia]; 
    Vec3d fa   = Vec3dZero;
    Vec3d fpi  = Vec3dZero; 
    //--- array aliases
    const int*    ings = neighs   [ia].array; // neighbors
    const int*    ingC = neighCell[ia].array; // neighbors cell index
    const double* bK   = bKs      [ia].array; // bond stiffness
    const double* bL   = bLs      [ia].array; // bond length
    const double* Kspi = Ksp      [ia].array; // pi-sigma stiffness
    const double* Kppi = Kpp      [ia].array; // pi-pi stiffness
    Vec3d* fbs  = fneigh   +ia*4;             // forces on bonds
    Vec3d* fps  = fneighpi +ia*4;             // forces on pi vectors
    const bool bColDampB   = colDamp.bBond && vapos; // if true we use collision damping
    const bool bColDampAng = colDamp.bAng  && vapos; // if true we use collision damping for non-bonded interactions
    const Quat4d& apar  = apars[ia]; // [c0, Kss, Ksp, c0_e] c0 is cos of equilibrium angle, Kss is bond stiffness, Ksp is pi-sigma stiffness, c0_e is cos of equilibrium angle for pi-electron interaction
    const double  piC0 = apar.w;     // cos of equilibrium angle for pi-electron interaction
    //--- Aux Variables 
    Quat4d  hs[4]; // bond vectors (normalized in .xyz ) and their inverse length in .w
    Vec3d   f1,f2; // working forces
    for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break in the loop
    // --------- Bonds Step
    for(int i=0; i<4; i++){ // loop over bonds
        int ing = ings[i];
        if(ing<0) break;
        //printf("ia %i ing %i \n", ia, ing ); 
        Vec3d  pi = apos[ing];
        Quat4d h; 
        h.f.set_sub( pi, pa );
        // Periodic Boundary Conditions
        if constexpr(_bPBC){   
            int ipbc = ingC[i]; 
            h.f.add( shifts[ipbc] );
        }
        const double l = h.f.normalize(); 
        h.e    = 1/l;
        hs [i] = h;
        if(ia<ing){     // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            if constexpr(_doBonds){  
                E+= evalBond( h.f, l-bL[i], bK[i], f1 ); 
                if constexpr(_bSubtractBondNonBond) { // subtract non-bonded interactions between atoms which have common neighbor
                    Vec3d fij=Vec3dZero;
                    Quat4d REQij = _mixREQ(REQs[ia],REQs[ing]);  // combine van der Waals parameters for the pair of atoms
                    Vec3d dp = h.f*l;
                    E -= getLJQH( dp, fij, REQij, R2damp ); // subtract non-bonded interactions 
                    if constexpr(_bClampNonBonded){ clampForce( fij, Fmax2 ); }
                    f1.sub(fij);
                }
                fbs[i].sub(f1);  fa.add(f1); 
            }
            double kpp = Kppi[i];
            if constexpr(_doPiPiI) if( (ing<nnode) && (kpp>1e-6) ){
                E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2); 
            }            
        }  
        double ksp = Kspi[i];
        if constexpr(_doPiSigma) if (ksp>1e-6){  
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2); 
        }        
    }
    if(doAngles){
    double  ssK,ssC0;
    Vec2d   cs0_ss;
    Vec3d*  angles_i;
    if constexpr(_bEachAngle) { 
        angles_i = angles+(ia*6);
    }else{ 
        ssK    = apar.z;
        cs0_ss = Vec2d{apar.x,apar.y};
        ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;
    }
    int iang=0;
    for(int i=0; i<3; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4d& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4d& hj = hs[j];    
            if constexpr(_bEachAngle) {  //
                cs0_ss = angles_i[iang].xy();  
                ssK    = angles_i[iang].z;
                iang++; 
            };
            if constexpr(_bAngleCosHalf ){ 
                E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            }else{ 
                E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            }
            fa    .sub( f1+f2  );  // apply force on the central atom
            if constexpr(_bSubtractAngleNonBond){ 
                Vec3d fij=Vec3dZero;
                Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]); 
                Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                E -= getLJQH( dp, fij, REQij, R2damp );
                if constexpr(_bClampNonBonded){ clampForce( fij, Fmax2 ); }
                f1.sub(fij); // 
                f2.add(fij);
            }
            fbs[i].add( f1     ); 
            fbs[j].add( f2     );
        }
    }
    }    

    fapos [ia]=fa; 
    fpipos[ia]=fpi;
    return E;
}




    __attribute__((hot))  
    double eval_atom_opt(const int ia){

        double E=0;
        const Vec3d pa   = apos [ia]; 
        Vec3d       fa   = Vec3dZero;
        const int*  ings = neighs[ia].array;

        Vec3d* fbs  = fneigh   +ia*4;

        const Quat4d& apar  = apars[ia];
        const double  piC0  = apar.w;

        //--- Aux Variables 
        Quat4d  hs[4];
        Vec3d   f1,f2;

        // ========================= Bonds
        { // Bonds

            const Vec3d hpi = pipos[ia]; 
            Vec3d fpi       = Vec3dZero; 
            const int* ingC = neighCell[ia].array;

            const double* bK   = bKs      [ia].array;
            const double* bL   = bLs      [ia].array;
            const double* Kspi = Ksp      [ia].array;
            const double* Kppi = Kpp      [ia].array;

            Vec3d* fps  = fneighpi +ia*4;
            for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break

            for(int i=0; i<4; i++){
                const int ing = ings[i];
                if(ing<0) break;

                const Vec3d  pi = apos[ing];
                Quat4d h; 
                h.f.set_sub( pi, pa );
                
                const int ipbc = ingC[i]; 
                //Vec3d sh = shifts[ipbc]; //apbc[i]  = pi + sh;
                h.f.add( shifts[ipbc] );
                const double l = h.f.normalize();
                h.e    = 1/l;
                hs [i] = h;

                if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
                    E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);
                    const double kpp = Kppi[i];
                    if( (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                        E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                    }
                } 
                const double ksp = Kspi[i];
                if( (ksp>1e-6) ){  
                    E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)
                }            
            }
            fpipos[ia]=fpi;

        } // Bonds


        // ========================= Angles
        const double R2damp=Rdamp*Rdamp;
        if(doAngles){
            const double ssK    = apar.z;
            const Vec2d  cs0_ss = Vec2d{apar.x,apar.y};
            const double ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

            int iang=0;
            for(int i=0; i<3; i++){
                int ing = ings[i];
                if(ing<0) break;
                const Quat4d& hi = hs[i];
                for(int j=i+1; j<4; j++){
                    int jng  = ings[j];
                    if(jng<0) break;
                    const Quat4d& hj = hs[j];    
                    E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
                    fa    .sub( f1+f2  );

                    /*
                    // ----- Error is HERE
                    if(bSubtractAngleNonBond){
                        Vec3d fij=Vec3dZero;
                        //Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                        Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);
                        Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                        //Vec3d dp   = hj.f*(1./hj.w) - hi.f*(1./hi.w);
                        //Vec3d dp   = apbc[j] - apbc[i];
                        E -= getLJQH( dp, fij, REQij, R2damp );
                        //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                        //bErr|=ckeckNaN( 1,3, (double*)&fij, [&]{ printf("atom[%i]fLJ2[%i,%i]",ia,i,j); } );
                        f1.sub(fij);
                        f2.add(fij);
                    }
                    */
                    fbs[i].add( f1     );
                    fbs[j].add( f2     );
                }
            }
        }    // doAngles
        fapos [ia]=fa; 
        return E;   
    }

double eval_atoms( bool bDebug=false, bool bPrint=false ){
    //FILE *file = fopen("out","w");
    //Etot,
    //Eb=0;Ea=0;Eps=0;EppT=0;EppI=0;
    double E=0;
    if  (bDebug){ for(int ia=0; ia<nnode; ia++){  E+=eval_atom_debug(ia,bPrint); }}
    else        { for(int ia=0; ia<nnode; ia++){  E+=eval_atom(ia);              }}
    //printf( "MYOUTPUT Ebond= %g Eangle= %g Edihed= %g Eimpr = %g Etot=%g\n", Eb, Ea, EppI, Eps, E );
    //fprintf( file, "Ebond= %g Eangle= %g Edihed= %g Eimpr= %g Etot=%g \n", Eb, Ea, EppI, Eps, E );
    //fclose(file);
    return E;
}

double evalKineticEnergy(){
    // In SI units:
    double mass=1; 
    double Ek=0;
    double vabsmean = 0;
    double msum = 0;
    for(int ia=0; ia<natoms; ia++){ 
        vabsmean = fabs( vapos[ia].norm()/10.180505710774743 ); // to Angstrom/fs
        double mi = (1.66053906660e-27 * mass);
        msum  += mi;
        Ek    += 0.5* mi* vapos[ia].norm2()*sq( 1e-10/1.0180505710774743e-14 );
    }
    double Ek_eV = Ek/1.602176634e-19;
    double v2mean = sqrt(2*Ek/msum);
    //printf( "vabsmean %g[m/s](%g[au]) v2mean %g[m/s](%g[au]) Ek=%g[J](Ek=%g[eV]) msum=%g[kg] \n", vabsmean*1e+5,vabsmean,    v2mean,v2mean/1e+5,   Ek, Ek_eV, msum );
    return Ek_eV;

    /*
    double Ek=0;
    double mass=1; // ToDo: get some masses later
    for(int ia=0; ia<natoms; ia++){ 
        Ek+=0.5*mass*vapos[ia].norm2();
    }
    Ek * = const_fs_timeu*const_fs_timeu; // convert to eV
    double v2mean = sqrt(Ek/natoms);
    printf( "v2mean %g [A/fs] \n", v2mean );
    // units conversion 
    //  velocity is in Angstrom/fs   1 Angstrom = 1.0e-10 m, 1 fs = 1.0e-15 s =>   1 Angstrom/fs = 1.0e+5 m/s 
    //  mass is in amu      1 amu = 1.66053906660E-27 kg
    //  1/2 m v^2 = 1/2 * 1.66053906660E-27 * (1.0e+5)^2 = 1.66053906660E-27 * 1.0e+10 = 1.66053906660E-17 J
    // eV = 1.602176634E-19 J
    // 1.66053906660e-17 / 1.602176634e-19 = 103.642695
    // 1.602176634e-19 / 1.66053906660e-17 = 0.00964853321
    Ek*=0.00964853321;
    return Ek;
    */
}

// move cappping atom to equilibrium distance from the node atoms
void addjustAtomCapLenghs(int ia){
    const Vec3d pa  = apos [ia]; 
    const int*    ings = neighs   [ia].array;
    //const int*    ingC = neighCell[ia].array;
    const double* bL   = bLs      [ia].array;
    for(int i=0; i<4; i++){
        int ing = ings[i];
        if(ing<nnode)continue;
        Vec3d  pi = apos[ing];
        Vec3d h; 
        h.set_sub( pi, pa );
        h.mul(bL[ing]/h.norm());
        apos[ing].set_add( pa, h ); 
    } 
}
// move cappping atoms to equilibrium distance from the node atoms
void addjustCapLenghs(){
    for(int ia=0; ia<nnode; ia++){ 
        addjustAtomCapLenghs(ia);
    }
}

// initialize directions of pi orbitals
void initPi( Vec3d* pbc_shifts, double Kmin=0.0001, double r2min=1e-4, bool bCheck=true ){
    //printf( "MMFFsp3_loc::initPi()\n" );

    // first set pi-orbitals on atoms determined by orthogonality of sigma bonds
    for(int ia=0; ia<nnode; ia++){ 
        if(vapos)vapos[natoms+ia]=Vec3dZero;
        const int*    ngs = neighs   [ia].array;
        const int*    ngC = neighCell[ia].array;
        const double* ks  = Ksp      [ia].array;
        Vec3d u,v,p;
        int nfound=0; 
        int j=0;
        while(j<4){ if(ks[j]>Kmin){ u=apos[ngs[j]]+pbc_shifts[ngC[j]]; nfound++; j++;break; }; j++; }
        while(j<4){ if(ks[j]>Kmin){ v=apos[ngs[j]]+pbc_shifts[ngC[j]]; nfound++; j++;break; }; j++; }
        while(j<4){ if(ks[j]>Kmin){ p=apos[ngs[j]]+pbc_shifts[ngC[j]]; nfound++; j++;break; }; j++; }
        //printf("...ia=%i nfound=%i\n",ia,nfound);
        if      ( nfound>=2 ){
            if( nfound==2 ){  p=apos[ia]; }
            //printf( "::[%i] u(%6.3f,%6.3f,%6.3f) v(%6.3f,%6.3f,%6.3f) p(%6.3f,%6.3f,%6.3f) nf=%1i \n", ia, u.x,u.y,u.z,  v.x,v.y,v.z, p.x,p.y,p.z, nfound );
            u.sub(p); v.sub(p);
            Vec3d pi; pi.set_cross( u,v );
            double r2 = pi.norm2();
            //printf( "::pi[%i](%6.3f,%6.3f,%6.3f) u(%6.3f,%6.3f,%6.3f) v(%6.3f,%6.3f,%6.3f) nf=%1i r2 %g \n", ia, u.x,u.y,u.z,  v.x,v.y,v.z, pi.x,pi.y,pi.z, nfound, sqrt(r2) );
            if( r2>r2min ){
                pi.mul(1/sqrt(r2));
                pipos[ia]=pi;
                //printf( "pi[%i]=pi(%5.3f,%5.3f,%5.3f) nfound %1i r2 %g \n", ia, pi.x,pi.y,pi.z, nfound, sqrt(r2) );
                continue;
            }
        }
        //printf( "pi[%i]=Vec3dZero nf=%1i \n", ia, nfound );
        //pipos[ia]=Vec3dZ; // pi cannot be define
        pipos[ia]=Vec3dZero; // pi cannot be define
    }

    // then set remaining pi-orbitals to be parallel to the pi-orbital bonded by strongest pi-pi interaction (Kpp)
    for(int ia=0; ia<nnode; ia++){ 
        const int*    ngs = neighs   [ia].array;
        const double* ks  = Kpp      [ia].array;
        int imax=-1;
        double kmax=0.0;
        double r2 = pipos[ia].norm2();
        if( r2>0.1 )continue;
        for(int i=0; i<4; i++){
            double k = ks[i];
            if(k>kmax){ 
                if(bCheck){
                    int ing=ngs[i];
                    if     ( ing<0     ){ printf("ERROR in MMFFsp3_loc::initPi() atom[%i] Kpp=%g but neihgbor[%i] is undefined(%i) => Exit() \n", ia, k,i, ing    ); exit(0); }    
                    else if( ing>nnode ){ printf("ERROR in MMFFsp3_loc::initPi() atom[%i] Kpp=%g but neihgbor[%i] is cap(nnode=%i) => Exit() \n", ia, k,i, nnode  ); exit(0); }
                }
                imax=i; kmax=k; 
            };
        }
        if( kmax>Kmin ){ 
            pipos[ia]=pipos[ ngs[imax] ]; 
            //printf("pi[%i] set by neigh[%i==%i] with Kpp=%g\n", ia, imax, ngs[imax], kmax ); 
        }else{
            pipos[ia]=Vec3dZ; // pi cannot be define
        }
        
    }
    //for(int ia=0; ia<nnode; ia++){  pipos[ia]=Vec3dZ ; } // Debug
}

// initialize directions of pi orbitals iteratively
__attribute__((hot))  
int relax_pi( int niter, double dt, double Fconv, double Flim=1000.0 ){
    double F2conv = Fconv*Fconv;
    double E=0,F2=0;
    int    itr=0;
    for(itr=0; itr<niter; itr++){
        {E=0;F2=0;}
        cleanForce();
        normalizePis();
        {Eb=0;Ea=0;Eps=0;EppT=0;EppI=0;}
        for(int ia=0; ia<nnode; ia++){ 
            E += eval_atom(ia);
        }
        for(int ia=0; ia<nnode; ia++){
            assemble_atom( ia );
        }
        for(int i=natoms; i<nvecs; i++){
            F2 +=move_atom_MD( i, dt, Flim, 0.9 ).z;
        }
        if(F2<F2conv){ return itr; }
    }
    return itr;
}

// normalize pi orbitals
__attribute__((hot))  
void normalizePis(){ 
    for(int i=0; i<nnode; i++){ pipos[i].normalize(); } 
}

// constrain atom to fixed position
void constrainAtom( int ia, double Kfix=1.0 ){
    printf( "constrainAtom(i=%i,K=%g)\n", ia, Kfix );
    constr[ia].f=apos[ia];
    constr[ia].w=Kfix;
};

// clear forces on all atoms and other DOFs
void cleanForce(){ 
    Etot=0;
    //for(int i=0; i<natoms; i++){ fapos [i].set(0.0);  } 
    //for(int i=0; i<nnode;  i++){ fpipos[i].set(0.0);  } 
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0;  } 
    // NOTE: We do not need clean fneigh,fneighpi because they are set in eval_atoms 
}

void cleanVelocity(){  for(int i=0; i<nvecs; i++){ vapos[i]=Vec3dZero; }  }

// add neighbor recoil forces to an atom (ia)
__attribute__((hot))  
void assemble_atom(int ia){
    Vec3d fa=Vec3dZero,fp=Vec3dZero;
    const int* ings = bkneighs[ia].array;
    bool bpi = ia<nnode; // if this is node atom, it has pi-orbital
    for(int i=0; i<4; i++){
        int j = ings[i];
        if(j<0) break;
        //if(j>=(nnode*4)){ printf("ERROR bkngs[%i|%i] %i>=4*nnode(%i)\n", ia, i, j, nnode*4 ); exit(0); }
        fa.add(fneigh  [j]);
        if(bpi){
            //printf( "assemble[%i,%i|%i] pi(%g,%g,%g) fp(%g,%g,%g) fpng(%g,%g,%g) \n", ia,i,j, pipos[ia].x,pipos[ia].y,pipos[ia].z, fpipos[ia].x,fpipos[ia].y,fpipos[ia].z, fneighpi[j].x,fneighpi[j].y,fneighpi[j].z );
            fp.add(fneighpi[j]);
            //ckeckNaN( 1,3, (double*)(fneighpi+j), [&]{ printf("assemble.fneighpi[%i]",j); } );
        }
    }
    fapos [ia].add( fa ); 
    
    if(bpi){ // if this is node atom, apply recoil to pi-orbital as well
        //if( pipos[ia].norm2()>1.5 ){ printf("pipos[%i](%g,%g,%g) not normalized !!! (assemble.iteration=%i) => Exit() \n", ia, pipos[ia].x,pipos[ia].y,pipos[ia].z, itr_DBG ); exit(0); };
        fpipos[ia].add( fp );
        fpipos[ia].makeOrthoU( pipos[ia] );  // subtract force component which change pi-vector size
        //ckeckNaN( 1,3, (double*)(fpipos+ia), [&]{ printf("assemble.makeOrthoU[%i] c %g pipos(%g,%g,%g) fpipos(%g,%g,%g)",ia,  pi.dot(fi),   pi.x,pi.y,pi.z, fi.x,fi.y,fi.z  ); } );
    }
}

// add neighbor recoil forces to all atoms
__attribute__((hot))  
void asseble_forces(){
    for(int ia=0; ia<natoms; ia++){
        assemble_atom(ia);
    }
}

// Full evaluation of MMFFsp3 intramolecular force-field
__attribute__((hot))  
double eval( bool bClean=true, bool bCheck=true ){
    //if(bClean){ cleanAll(); }
    //printf( "print_apos() BEFORE\n" );print_apos();
    if(bClean)cleanForce();
    normalizePis();
    //printf( "print_apos() AFTER \n" ); print_apos();
    Etot += eval_atoms();
    //if(idebug){printf("CPU BEFORE assemble() \n"); printDebug();} 
    asseble_forces();
    //if(bTorsion) EppI +=eval_torsions();
    //Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}

// debug evaluation of MMFFsp3 intramolecular force-field
double eval_check(){
    if(verbosity>0){
        printf(" ============ check MMFFsp3_loc START\n " );
        printSizes();
    }
    //print_pipos();
    cleanForce();
    if(vapos)cleanVelocity();
    eval();
    checkNans();
    if(verbosity>0)printf(" ============ check MMFFsp3_loc DONE\n " );
    return Etot;
}





// Loop which iteratively evaluate MMFFsp3 intramolecular force-field and move atoms to minimize energy
// ToDo: OpenMP paraelization atempt
__attribute__((hot))  
int run( int niter, double dt, double Fconv, double Flim, double damping=0.1 ){
    double F2conv = Fconv*Fconv;
    double E=0,ff=0,vv=0,vf=0;

    double cdamp = colDamp.update( dt );

    //printf( "MMFFsp3_loc::run(bCollisionDamping=%i) niter %i dt %g Fconv %g Flim %g damping %g collisionDamping %g \n", bCollisionDamping, niter, dt, Fconv, Flim, damping, collisionDamping );
    printf( "MMFFsp3_loc::run(niter=%i,bCol(B=%i,A=%i,NB=%i)) dt %g damp(cM=%g,cB=%g,cA=%g,cNB=%g)\n", niter, colDamp.bBond, colDamp.bAng, colDamp.bNonB, dt, 1-cdamp, colDamp.cdampB*dt, colDamp.cdampAng*dt, colDamp.cdampNB*dt );

    int    itr;
    //if(itr_DBG==0)print_pipos();
    //bool bErr=0;
    for(itr=0; itr<niter; itr++){
        E=0;
        // ------ eval MMFF
        for(int ia=0; ia<natoms; ia++){ 
            {             fapos[ia       ] = Vec3dZero; } // atom pos force
            if(ia<nnode){ fapos[ia+natoms] = Vec3dZero; } // atom pi  force
            if(ia<nnode)E += eval_atom(ia);
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+ia), [&]{ printf("eval.MMFF[%i]",ia); } );
            //E += evalLJQs_ng4_PBC_atom( ia ); 

            if(bPBC){ E+=evalLJQs_ng4_PBC_atom_omp( ia ); }
            else    { 
                E+=evalLJQs_ng4_atom_omp    ( ia ); 
                if( colDamp.bNonB ){ evalCollisionDamp_atom_omp( ia, colDamp.nonB, colDamp.dRcut1, colDamp.dRcut2 ); }
            } 
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+ia), [&]{ printf("eval.NBFF[%i]",ia); } );
        }
        // ---- assemble (we need to wait when all atoms are evaluated)
        for(int ia=0; ia<natoms; ia++){
            assemble_atom( ia );
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+ia), [&]{ printf("assemble[%i]",ia); } );
        }
        // ------ move
        cvf = Vec3dZero;
        for(int i=0; i<nvecs; i++){
            //F2 += move_atom_GD( i, dt, Flim );
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+i), [&]{ printf("move[%i]",i); } );
            cvf.add( move_atom_MD( i, dt, Flim, cdamp ) );
            //move_atom_MD( i, 0.05, 1000.0, 0.9 );
            //F2 += move_atom_kvaziFIRE( i, dt, Flim );
        }
        if(cvf.z<F2conv)break;
        if(cvf.x<0){ cleanVelocity(); };
        //itr_DBG++;
    }
    return itr;
}

// Loop which iteratively evaluate MMFFsp3 intramolecular force-field and move atoms to minimize energy, wih OpenMP parallelization
__attribute__((hot))  
int run_omp( int niter, double dt, double Fconv, double Flim, double damping=0.1 ){
    double F2conv = Fconv*Fconv;
    double E=0,ff=0,vv=0,vf=0;
    double cdamp = 1-damping; if(cdamp<0)cdamp=0;
    int    itr=0;
    #pragma omp parallel shared(E,ff,vv,vf) private(itr)
    for(itr=0; itr<niter; itr++){
        // This {} should be done just by one of the processors
        #pragma omp single
        {E=0;ff=0;vv=0;vf=0;}
        // ------ eval MMFF
        #pragma omp for reduction(+:E)
        for(int ia=0; ia<natoms; ia++){ 
            if(verbosity>3)[[unlikely]]{ printf( "atom[%i]@cpu[%i/%i]\n", ia, omp_get_thread_num(), omp_get_num_threads()  ); }
            {             fapos[ia       ] = Vec3dZero; } // atom pos force
            if(ia<nnode){ fapos[ia+natoms] = Vec3dZero; } // atom pi  force
            if(ia<nnode)E += eval_atom(ia);
            //E += evalLJQs_ng4_PBC_atom( ia ); 
            //E += evalLJQs_ng4_PBC_atom_omp( ia ); 
            if(bPBC){ E+=evalLJQs_ng4_PBC_atom_omp( ia ); }
            else    { E+=evalLJQs_ng4_atom_omp    ( ia ); } 
        }
        // ---- assemble (we need to wait when all atoms are evaluated)
        #pragma omp for
        for(int ia=0; ia<natoms; ia++){
            assemble_atom( ia );
        }
        // ------ move
        #pragma omp for reduction(+:ff,vv,vf)
        for(int i=0; i<nvecs; i++){
            const Vec3d cvf_ = move_atom_MD( i, dt, Flim, cdamp );
            //const Vec3d cvf_ += move_atom_MD( i, 0.05, 1000.0, 0.9 );
            ff += cvf_.x; vv += cvf_.y; vf += cvf_.z;
        }
        #pragma omp single
        {
            cvf.x=ff; cvf.y=vv; cvf.z=vf;
            if(cvf.x<0){ cleanVelocity(); };
            //if(cvf.z<F2conv)break;
            if(verbosity>2)[[unlikely]]{printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(ff), omp_get_num_threads() );}
        }
    }
    return itr;
}

// flip pi-orbitals to be in the same half-space as the given reference vector (ax) 
void flipPis( Vec3d ax ){
    for(int i=0; i<nnode; i++){
        double c = pipos[i].dot(ax);
        if( c<0 ){ pipos[i].mul(-1); } 
    }
}

// update atom positions using gradient descent
__attribute__((hot))  
inline double move_atom_GD(int i, float dt, double Flim){
    Vec3d  f   = fapos[i];
    Vec3d  p = apos [i];
    double fr2 = f.norm2();
    if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    if(bPi){ f.add_mul( p, -p.dot(f) ); }
    p.add_mul( f, dt );
    if(bPi)p.normalize();
    apos [i] = p;
    //fapos[i] = Vec3dZero;
    return fr2;
}
// update atom positions using molecular dynamics (damped leap-frog)
__attribute__((hot))  
inline Vec3d move_atom_Langevin( int i, const float dt, const double Flim,  const double gamma_damp=0.1, double T=300 ){
    Vec3d f = fapos[i];
    Vec3d v = vapos[i];
    Vec3d p = apos [i];
    const bool bPi = i>=natoms;
    if(bPi)f.add_mul( p, -p.dot(f) ); 
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
    if(bPi)v.add_mul( p, -p.dot(v) );                     
    p.add_mul( v, dt );
    if(bPi)p.normalize();
    
    apos [i] = p;
    vapos[i] = v;

    return cvf;
}
Vec3d move_Langevin( const float dt, const double Flim,  const double gamma_damp=0.1, double T=300 ){
    Vec3d cvf = Vec3dZero;
    for(int i=0; i<nvecs; i++){
        cvf.add( move_atom_Langevin( i, dt, Flim, gamma_damp, T ) );
    }
    return cvf;
};


// update atom positions using molecular dynamics (damped leap-frog)
__attribute__((hot))  
inline Vec3d move_atom_MD( int i, const float dt, const double Flim, const double cdamp=0.9 ){
    Vec3d f = fapos[i];
    Vec3d v = vapos[i];
    Vec3d p = apos [i];
    //const double fr2 = f.norm2();
    //if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    // bool b=false;
    // b|=ckeckNaN_d( 1,3, (double*)&p, "p.1" );
    // b|=ckeckNaN_d( 1,3, (double*)&v, "v.1" );
    // b|=ckeckNaN_d( 1,3, (double*)&f, "f.1" );
    if(bPi)f.add_mul( p, -p.dot(f) );           //b|=ckeckNaN_d( 1,3, (double*)&f, "f.2" );
    const Vec3d  cvf{ v.dot(f), v.norm2(), f.norm2() };
    v.mul    ( cdamp );
    v.add_mul( f, dt );                         //b|=ckeckNaN_d( 1,3, (double*)&v, "v.2" );
    if(bPi)v.add_mul( p, -p.dot(v) );           //b|=ckeckNaN_d( 1,3, (double*)&v, "v.3" );
    p.add_mul( v, dt );
    if(bPi)p.normalize();
    // if( bPi &&( p.norm2()>1.5 )){ printf("pipos[%i/%i](%g,%g,%g) not normalized !!! (move.iteration=%i) => Exit() \n", i, natoms, p.x,p.y,p.z, itr_DBG ); exit(0); };
    // b|=ckeckNaN_d( 1,3, (double*)&p, "p.4" );
    // b|=ckeckNaN_d( 1,3, (double*)&v, "v.4" );
    // b|=ckeckNaN_d( 1,3, (double*)&f, "f.4" );
    // if(b){ printf("ERROR NaNs in move_atom_MD[%i] => Exit() \n", i ); exit(0);}
    apos [i] = p;
    vapos[i] = v;
    //fapos[i] = Vec3dZero;
    return cvf;
}

// update atom positions using FIRE (Fast Inertial Relaxation Engine)
__attribute__((hot))  
inline double move_atom_FIRE( int i, float dt, double Flim, double cv, double cf ){
    Vec3d  f   = fapos[i];
    Vec3d  v   = vapos[i];
    Vec3d  p   = apos [i];
    double fr2 = f.norm2();
    //if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    if(bPi)f.add_mul( p, -p.dot(f) ); 
    v.mul(             cv );
    v.add_mul( f, dt + cf );
    if(bPi) v.add_mul( p, -p.dot(v) );
    p.add_mul( v, dt );
    if(bPi)  p.normalize();
    apos [i] = p;
    vapos[i] = v;
    //fapos[i] = Vec3dZero;
    return fr2;
}

// update atom positions using modified FIRE algorithm
__attribute__((hot))  
inline double move_atom_kvaziFIRE( int i, float dt, double Flim ){
    Vec3d  f   = fapos[i];
    Vec3d        v   = vapos[i];
    Vec3d        p   = apos [i];
    double fr2 = f.norm2();
    if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    if(bPi)f.add_mul( p, -p.dot(f) ); 
    double vv  = v.norm2();
    double ff  = f.norm2();
    double vf  = v.dot(f);
    double c   = vf/sqrt( vv*ff + 1.e-16    );
    double v_f =    sqrt( vv/( ff + 1.e-8)  );
    double cv;
    double cf = cos_damp_lin( c, cv, 0.01, -0.7,0.0 );
    v.mul(                 cv );
    v.add_mul( f, dt + v_f*cf );
    if(bPi) v.add_mul( p, -p.dot(v) );
    p.add_mul( v, dt );
    if(bPi)  p.normalize();
    apos [i] = p;
    vapos[i] = v;
    //fapos[i] = Vec3dZero;
    return fr2;
}

// update atom positions using gradient descent
__attribute__((hot))  
double move_GD(float dt, double Flim=100.0 ){
    double F2sum=0;
    for(int i=0; i<nvecs; i++){
        F2sum += move_atom_GD(i, dt, Flim);
    }
    return F2sum;
}

Vec3d shiftBack(bool bPBC=false){
    Vec3d cog=Vec3dZero;
    for(int i=0; i<natoms; i++){ cog.add( apos[i] ); }
    cog.mul( 1.0/natoms );
    //Vec3d cog=average( ffl.natoms, ffl.apos );  
    if( bPBC ){
        Vec3i g  = invLvec.nearestCell( cog );
        cog.add_lincomb( g.x,g.y,g.z, lvec.a, lvec.b, lvec.c );
    }
    for(int i=0; i<natoms; i++){ apos[i].sub( cog ); }
    return cog;
}

// make list of back-neighbors
void makeBackNeighs( bool bCapNeighs=true ){
    for(int i=0; i<natoms; i++){ bkneighs[i]=Quat4i{-1,-1,-1,-1}; }; // clear back-neighbors to -1
    for(int ia=0; ia<nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = neighs[ia].array[j];
            if( ja<0 )continue;
            //NOTE: We deliberately ignore back-neighbors from caping atoms 
            bool ret = addFirstEmpty( bkneighs[ja].array, 4, ia*4+j, -1 ); // add neighbor to back-neighbor-list on the first empty place
            if(!ret){ printf("ERROR in MMFFsp3_loc::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
    //for(int i=0; i<natoms; i++){printf( "bkneigh[%i] (%i,%i,%i,%i) \n", i, bkneighs[i].x, bkneighs[i].y, bkneighs[i].z, bkneighs[i].w );}
    //checkBkNeighCPU();
    if(bCapNeighs){   // set neighbors for capping atoms
        for(int ia=nnode; ia<natoms; ia++){ 
            if( bkneighs[ia].x<0 ){
                printf( "ERROR makeBackNeighs() capping atom[%i] has not back-neighbor (bkneighs[%i]==%i) => exit()\n", ia, ia, bkneighs[ia].x ); exit(0);
                //continue;
            }
            neighs[ia]=Quat4i{-1,-1,-1,-1};  neighs[ia].x = bkneighs[ia].x/4;  
            //printf("makeBackNeighs().bCapNeighs [ia=%i] ng.x=%i bkng.x=%i \n", ia, neighs[ia].x, bkneighs[ia].x );
        }
    }
}

// make list of neighbors cell index (in periodic boundary conditions), by going through all periodic images
void makeNeighCells( const Vec3i nPBC_ ){ 
    nPBC=nPBC_;
    //printf( "makeNeighCells() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    for(int ia=0; ia<natoms; ia++){
        Quat4i ngC = Quat4i{-1,-1,-1,-1};
        for(int j=0; j<4; j++){
            const int ja = neighs[ia].array[j];
            if( ja<0 )continue;
            const Vec3d d0 = apos[ja] - apos[ia];
            int ipbc =  0;
            int imin = -1;
            double r2min = 1e+300;
            // go through all periodic images and find nearest distance
            for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){   
                Vec3d d = d0 + (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz); 
                double r2 = d.norm2();
                //printf( "[%i,%i][%i,%i,%i] %g   (%g,%g,%g)\n", ia,ja, ix,iy,iz, sqrt(r2), d.x,d.y,d.z );
                if(r2<r2min){   // find nearest distance
                    r2min = r2;
                    imin  = ipbc;
                }
                ipbc++; 
            }}}
            ngC.array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
        }
        //printf("\n", ngC.x,ngC.y,ngC.z,ngC.w);
        neighCell[ia]=ngC; // set neighbor cell index
    }
    //printNeighs();
}

// make list of neighbors cell index (in periodic boundary conditions) using precomputed pbc_shifts
void makeNeighCells( int npbc, Vec3d* pbc_shifts ){ 
    for(int ia=0; ia<natoms; ia++){
        for(int j=0; j<4; j++){
            //printf("ngcell[%i,j=%i] \n", ia, j);
            int ja = neighs[ia].array[j];
            //printf("ngcell[%i,ja=%i] \n", ia, ja);
            if( ja<0 )continue;
            const Vec3d d = apos[ja] - apos[ia];

            // ------- Brute Force method
            int imin=-1;
            float r2min = 1.e+300;
            for( int ipbc=0; ipbc<npbc; ipbc++ ){   
                Vec3d shift = pbc_shifts[ipbc]; 
                shift.add(d);
                float r2 = shift.norm2();
                if(r2<r2min){   // find nearest distance
                    r2min=r2;
                    imin=ipbc;
                }
            }

            /*
            // -------- Fast method
            { //Debug
                Vec3d u;
                invLvec.dot_to( d, u );
                int ix = 1-(int)(u.x+1.5);
                int iy = 1-(int)(u.y+1.5);
                int iz = 1-(int)(u.z+1.5);
                Vec3d dpbc = lvec.a*ix + lvec.b*iy + lvec.c*iz;
                printf( "NeighCell[%i,%i] ipbc %i pbc_shift(%6.3f,%6.3f,%6.3f) dpbc(%6.3f,%6.3f,%6.3f) iabc(%2i,%2i,%2i) u(%6.3f,%6.3f,%6.3f)\n", ia, ja, imin, pbc_shifts[imin].x,pbc_shifts[imin].y,pbc_shifts[imin].z,  dpbc.x,dpbc.y,dpbc.z, ix,iy,iz, u.x,u.y,u.z  );
            }
            */

            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
            neighCell[ia].array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i ---- \n", ia, ja, imin);
        }
    }
}

// ================== print functions  

void printSizes     (      ){ printf( "MMFFf4::printSizes(): nDOFs(%i) natoms(%i) nnode(%i) ncap(%i) nvecs(%i) npbc(%i)\n", nDOFs,natoms,nnode,ncap,nvecs,npbc ); };
void printAtomParams(int ia){ printf("atom[%i] t(%i) ngs{%3i,%3i,%3i,%3i} par(%5.3f,%5.3f,%5.3f,%5.3f)  bL(%5.3f,%5.3f,%5.3f,%5.3f) bK(%6.3f,%6.3f,%6.3f,%6.3f)  Ksp(%5.3f,%5.3f,%5.3f,%5.3f) Kpp(%5.3f,%5.3f,%5.3f,%5.3f) \n", ia, atypes[ia], neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,    apars[ia].x,apars[ia].y,apars[ia].z,apars[ia].w,    bLs[ia].x,bLs[ia].y,bLs[ia].z,bLs[ia].w,   bKs[ia].x,bKs[ia].y,bKs[ia].z,bKs[ia].w,     Ksp[ia].x,Ksp[ia].y,Ksp[ia].z,Ksp[ia].w,   Kpp[ia].x,Kpp[ia].y,Kpp[ia].z,Kpp[ia].w  ); };
void printNeighs    (int ia){ printf("atom[%i] neigh{%3i,%3i,%3i,%3i} neighCell{%3i,%3i,%3i,%3i} \n", ia, neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,   neighCell[ia].x,neighCell[ia].y,neighCell[ia].z,neighCell[ia].w ); };
void printBKneighs  (int ia){ printf("atom[%i] bkngs{%3i,%3i,%3i,%3i} \n", ia, bkneighs[ia].x,bkneighs[ia].y,bkneighs[ia].z,bkneighs[ia].w ); };
void printAtomParams(      ){ printf("MMFFsp3_loc::printAtomParams()\n" ); for(int i=0; i<nnode;  i++){ printAtomParams(i); }; };
void printNeighs    (      ){ printf("MMFFsp3_loc::printNeighs()\n"     ); for(int i=0; i<natoms; i++){ printNeighs    (i);     }; };
void printBKneighs  (      ){ printf("MMFFsp3_loc::printBKneighs()\n"   ); for(int i=0; i<natoms; i++){ printBKneighs  (i);   }; };
void print_pipos    (      ){ printf("MMFFsp3_loc::print_pipos()\n"     ); for(int i=0; i<nnode;  i++){ printf( "pipos[%i](%g,%g,%g) r=%g\n", i, pipos[i].x,pipos[i].y,pipos[i].z, pipos[i].norm() ); } }
void print_apos     (      ){ printf("MMFFsp3_loc::print_apos()\n"      ); for(int i=0; i<natoms; i++){ printf( "apos [%i](%g,%g,%g)\n",      i, apos[i].x ,apos[i].y ,apos[i].z                   ); } }
void print_pbc_shifts(     ){ printf("MMFFsp3_loc::print_pbc_shifts()\n"); for(int i=0; i<npbc;   i++){ printf( "pbc_shifts[%i](%g,%g,%g)\n", i, shifts[i].x,shifts[i].y,shifts[i].z                   ); } }

void print_constrains(){ printf("MMFFsp3_loc::print_constrains()\n"   );    for(int i=0; i<natoms;   i++){ if(constr[i].w>0)printf( "constr[%3i](%8.5f,%8.5f,%8.5f|%g)K(%8.5f,%8.5f,%8.5f|%8.5f)\n", i, constr[i].x,constr[i].y,constr[i].z,constr[i].w,  constrK[i].x,constrK[i].y,constrK[i].z,constrK[i].w   ); }  }

void printAngles(int ia){
    Vec3d* angles_i = angles+(ia*6);
    int* ings       = neighs[ia].array; 
    int iang=0;
    for(int i=0; i<3; i++){
        int ing = ings[i];
        if(ing<0) break;
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            if(bEachAngle){
                Vec2d cs0_ss = angles_i[iang].xy();  
                double ssK    = angles_i[iang].z;
                printf( "atom[%i|%i]types{%i,%i,%i} ssK %g cs0(%g,%g) \n", ia,iang,  atypes[ing], atypes[ia], atypes[jng], ssK, cs0_ss.x, cs0_ss.y  );
                iang++; 
            }
        }
    }
}
void printAngles(      ){ printf("MMFFsp3_loc::printAngles()\n"); for(int i=0; i<nnode; i++){ printAngles(i); } }

void printTorsions(      ){ printf("MMFFsp3_loc::printTorsions()\n"); for(int i=0; i<ntors; i++){ printf( "torsion[%i]{%i,%i,%i,%i} {%i,%i,%i,%i}\n", i, tors2atom[i].x,tors2atom[i].y,tors2atom[i].z,tors2atom[i].w,   torsParams[i].x,torsParams[i].y,torsParams[i].z,torsParams[i].w  ); } }

void printAtomsConstrains( bool bWithOff=false ){ printf("MMFFsp3_loc::printAtomsConstrains()\n"); for(int i=0; i<natoms; i++){ if(bWithOff || (constr[i].w>0.0f) )printf( "consrt[%i](%g,%g,%g|K=%g)\n", i, constr[i].x,constr[i].y,constr[i].z,constr[i].w ); } }

void printDebug(  bool bNg=true, bool bPi=true, bool bA=true ){
    printf( "MMFFsp3_loc::printDebug()\n" );
    if(bA)for(int i=0; i<natoms; i++){
        printf( "CPU[%i] ", i );
        //printf( "bkngs{%2i,%2i,%2i,%2i,%2i} ",         bkNeighs[i].x, bkNeighs[i].y, bkNeighs[i].z, bkNeighs[i].w );
        printf( "fapos{%6.3f,%6.3f,%6.3f} ", fapos[i].x, fapos[i].y, fapos[i].z );
        //printf(  "avel{%6.3f,%6.3f,%6.3f} ", avel[i].x, avel[i].y, avel[i].z);
        printf(  "apos{%6.3f,%6.3f,%6.3f} ", apos[i].x, apos[i].y, apos[i].z );
        printf( "\n" );
    }
    if(bPi)for(int i=0; i<nnode; i++){
        int i1=i+natoms;
        printf( "CPU[%i] ", i1 );
        printf(  "fpipos{%6.3f,%6.3f,%6.3f} ", fapos[i1].x, fapos[i1].y, fapos[i1].z );
        //printf(  "vpipos{%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z );
        printf(   "pipos{%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z );
        printf( "\n" );
    }
    if(bNg)for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
        int i1=i*4+j;
        //int i2=(i+natoms)*4+j;
        printf( "CPU[%i,%i] ", i, j );
        printf( "fneigh  {%6.3f,%6.3f,%6.3f} ", fneigh  [i1].x, fneigh  [i1].y, fneigh  [i1].z );
        printf( "fneighpi{%6.3f,%6.3f,%6.3f} ", fneighpi[i1].x, fneighpi[i1].y, fneighpi[i1].z );
        printf( "\n" );
    }}
}

// check if there are NaNs in the arrays (and exit with error-message if there are)
bool checkNans( bool bExit=true, bool bNg=true, bool bPi=true, bool bA=true ){
    bool ret = false;
    int ia = whereNaN_d( natoms,  3, (double*)fapos,   "fapos"   ); if(ia>=0){ eval_atom_debug(ia,true);  printf("ERROR NaNs detected in fapos[%i]\n", ia ); exit(0); }
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*) apos,   "apos"  );
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*)fapos,  "fapos"  );
    if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*) pipos,  "pipos" );
    if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*)fpipos, "fpipos" );
    if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneigh,  "fneigh"   );
    if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneighpi,"fneighpi" );
    if(bExit&&ret){ printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); 
        printSizes();
        printAtomParams();
        printNeighs();
        print_pbc_shifts();
        printDebug( false, false );
        eval_atoms(true,true);
        printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); 
        exit(0); 
    };
    return ret;
}

// rotate selected atoms including their caps
void rotateNodes(int n, int* sel, Vec3d p0, Vec3d ax, double phi ){
    //printf( "MMFFsp3_loc::rotateNodes() nsel=%i phi=%g ax(%g,%g,%g) p0(%g,%g,%g) \n", n,  phi, ax.x,ax.y,ax.z,   p0.x,p0.y,p0.z );
    ax.normalize();
    double ca=cos(phi);
    double sa=sin(phi);
    for(int i=0;i<n; i++){
        int ia = sel[i];
        //printf( "MMFFsp3_loc::rotateNodes() atom[%i](%g,%g,%g) cs(%g,%g) \n", ia, apos[ia].x,apos[ia].y,apos[ia].z, ca,sa );
        apos [ia].rotate_csa( ca, sa, ax, p0 );
        if(ia>=nnode)continue;
        pipos[ia].rotate_csa( ca, sa, ax     );
        int* ngs=neighs[ia].array; 
        for(int j=0;j<4;j++){
            int ja = ngs[j];
            if(ja>=0){ if(ja>nnode) apos[ ja  ].rotate_csa( ca, sa, ax, p0 ); }
        }
    }
}

// redistribute charges from atom to capping electron pair
void chargeToEpairs( Quat4d* REQs, int* atypes, double cQ=-0.2, int etyp=-1 ){
    printf( "MMFFsp3_loc::chargeToEpairs() cQ %g etyp %i \n", cQ, etyp );
    // ToDo: we should operate on element-type not atom-type basis ( ielem = params.atypes[ityp].element ),   this is because we can have multiple electron-pair-types (depending on the atom-type to which they are attached)
    for( int ia=0; ia<nnode; ia++ ){
        int* ngs=neighs[ia].array; 
        for( int j=0; j<4; j++ ){
            int ja = ngs[j]; 
            if(ja<0) continue;
            if( atypes[ja]==etyp ){ REQs[ja].z+=cQ; REQs[ia].z-=cQ; };
        }
    }
}
void chargeToEpairs( double cQ=-0.2, int etyp=-1 ){ chargeToEpairs( REQs, atypes, cQ, etyp ); }

// ========== Measure functions

// measure cos of angle between two pi-orbitals
inline double measureCosPiPi(int ia, int ib, bool bRenorm=true){
    double c = pipos[ia].dot(pipos[ib]);
    if(bRenorm){ c/=sqrt( pipos[ia].norm2()* pipos[ib].norm2() ); }
    return c;
}
inline double measureAnglePiPi(int ia, int ib, bool bRenorm=true ){ return acos( measureCosPiPi(ia, ib, bRenorm ) ); }
inline double measureCosSigmaPi(int ipi, int ia, int ib){
    Vec3d b = apos[ib]-apos[ia];
    return pipos[ipi].dot(b)/sqrt( pipos[ipi].norm2()*b.norm2() );
}
inline double measureAngleSigmaPi(int ipi, int ia, int ib){ return acos( measureCosSigmaPi(ipi, ia, ib ) ); }




double eval_atom_debug(const int ia, bool bPrint=true){
    if(bPrint)printf( "#### MMFFsp3_loc::eval_atom_debug(%i)\n", ia );
    double E=0;
    const Vec3d pa  = apos [ia]; 
    const Vec3d hpi = pipos[ia]; 

 
    Vec3d fa   = Vec3dZero;
    Vec3d fpi  = Vec3dZero; 
    
    //--- array aliases
    const int*    ings = neighs   [ia].array;
    const int*    ingC = neighCell[ia].array;
    const double* bK   = bKs      [ia].array;
    const double* bL   = bLs      [ia].array;
    const double* Kspi = Ksp      [ia].array;
    const double* Kppi = Kpp      [ia].array;
    Vec3d* fbs  = fneigh   +ia*4;
    Vec3d* fps  = fneighpi +ia*4;

    const Quat4d& apar  = apars[ia];
    const double  piC0 = apar.w;

    const bool bColDampB   = colDamp.bBond && vapos; // if true we use collision damping
    const bool bColDampAng = colDamp.bAng  && vapos; // if true we use collision damping for non-bonded interactions

    //--- Aux Variables 
    Quat4d  hs[4];
    Vec3d   f1,f2;

    for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break

    double Eb=0;
    double Ea=0;
    double EppI=0;
    double Eps=0;

    // --------- Bonds Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        if(ing<0) break;

        //printf("ia %i ing %i \n", ia, ing ); 
        Vec3d  pi = apos[ing];
        Quat4d h; // {x,y,z|w} = {f.x,f.y,f.z|e}
        h.f.set_sub( pi, pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        
        if(bPBC){   
            if(shifts){
                int ipbc = ingC[i]; 
                //Vec3d sh = shifts[ipbc]; //apbc[i]  = pi + sh;
                h.f.add( shifts[ipbc] );
                //if( (ia==0) ){ printf("ffls[%i] atom[%i,%i=%i] ipbc %i shifts(%g,%g,%g)\n", id, ia,i,ing, ipbc, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z); };
                //if( (ipbc!=4) ){ printf("ffls[%i] atom[%i,%i=%i] ipbc %i shifts(%g,%g,%g)\n", id, ia,i,ing, ipbc, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z); };
            }else{
                Vec3i g  = invLvec.nearestCell( h.f );
                // if(ia==ia_DBG){
                //     Vec3d u; invLvec.dot_to(h.f,u);
                //     printf( "CPU:bond[%i,%i] u(%6.3f,%6.3f,%6.3f) shi(%6.3f,%6.3f,%6.3f) \n", ia, ing, u.x,u.y,u.z,   (float)g.x,(float)g.y,(float)g.z );
                // }
                Vec3d sh = lvec.a*g.x + lvec.b*g.y + lvec.c*g.z;
                h.f.add( sh );
                //apbc[i] = pi + sh;
            }
        }

        //printf( "h[%i,%i] r_old %g r_new %g \n", ia, ing, h_bak.norm(), h.f.norm() );
        double l = h.f.normalize();

        h.e    = 1/l;
        hs [i] = h;



        //if(ia==ia_DBG) printf( "ffl:h[%i|%i=%i] l %g h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", ia,i,ing, l, h.x,h.y,h.z, apos[ing].x,apos[ing].y,apos[ing].z, pa.x,pa.y,pa.z );

        if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            if(doBonds){

                if(bColDampB){ // 
                    double dv   = h.f.dot( vapos[ing] - vapos[ia] );
                    double fcol = colDamp.cdampB * 0.5 * dv;   // col_damp ~ collision_damping/(dt*ndampstep);     f = m*a = m*dv/dt
                    f1.set_mul( h.f, fcol );
                    if(bPrint)printf( "ffl:bondCol[%i|%i=%i] colDamp=%g dv=%g fcol=%g h(%g,%g,%g) f(%g,%g,%g) E %g\n", ia,i,ing, colDamp.cdampB, dv, fcol, h.x,h.y,h.z, f1.x,f1.y,f1.z);
                    fbs[i].sub(f1);  fa.add(f1); 
                }

                double Ebi = evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);
                E +=Ebi; 
                Eb+=Ebi; 
                if(bPrint)printf( "ffl:bond[%i|%i=%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) E %g\n", ia,i,ing, bK[i],bL[i], l, h.x,h.y,h.z, f1.x,f1.y,f1.z,  Ebi);
            }
            double kpp = Kppi[i];
            if( (doPiPiI) && (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                double EppIi = evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                E   +=EppIi; 
                EppI+=EppIi;
                if(bPrint)printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g), EppI=%g\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z, EppIi  );
            }
            // ToDo: triple bonds ?
            
        } 
        
        // pi-sigma 
        //if(bPi){    
        double ksp = Kspi[i];
        if( doPiSigma && (ksp>1e-6) ){  
            double Epsi = evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)
            E  +=Epsi; 
            Eps+=Epsi;
            if(bPrint)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g) E %g \n", ia,ing, ksp,piC0, hpi.dot(h.f), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z, Epsi );
        }
        //}
        
    }

    //printf( "MMFF_atom[%i] cs(%6.3f,%6.3f) ang=%g [deg]\n", ia, cs0_ss.x, cs0_ss.y, atan2(cs0_ss.y,cs0_ss.x)*180./M_PI );
    // --------- Angle Step
    const double R2damp=Rdamp*Rdamp;
    if(doAngles){

    double  ssK,ssC0;
    Vec2d   cs0_ss;
    Vec3d*  angles_i;
    if(bEachAngle){
        angles_i = angles+(ia*6);
    }else{
        ssK    = apar.z;
        cs0_ss = Vec2d{apar.x,apar.y};
        ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
    }

    int iang=0;
    for(int i=0; i<3; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4d& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4d& hj = hs[j];    
            if(bEachAngle){
                // 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
                //  0    1    2    3    4    5 
                cs0_ss = angles_i[iang].xy();  
                ssK    = angles_i[iang].z;
                //printf( "types{%i,%i,%i} ssK %g cs0(%g,%g) \n", atypes[ing], atypes[ia], atypes[jng], ssK, cs0_ss.x, cs0_ss.y  );
                iang++; 
            };
            //bAngleCosHalf = false;
            //double Eai;
            //printf( "bAngleCosHalf= %i\n", bAngleCosHalf);
            double Eai;
            if( bAngleCosHalf ){
                Eai = evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            }else{ 
                Eai = evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            }
            E +=Eai; 
            Ea+=Eai;
            if(bPrint)printf( "ffl:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g) E %g\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z, f2.x,f2.y,f2.z, Eai );
            fa    .sub( f1+f2  );


            if(bColDampAng){
                Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                double dv   = dp.dot( vapos[ing] - vapos[ia] );
                double fcol = colDamp.cdampAng * 0.5 * dv;   // col_damp ~ collision_damping/(dt*ndampstep);     f = m*a = m*dv/dt
                dp.mul( fcol/dp.norm2() );
                f1.sub(dp);
                f2.add(dp);
                if(bPrint)printf( "ffl:colDampAng[%i|%i,%i] cDampAng=%g dv=%g fcol=%g fij(%g,%g,%g)\n", ia,ing,jng, colDamp.cdampAng, dv, fcol, dp.x,dp.y,dp.z );
            }

            // ----- Error is HERE
            if(bSubtractAngleNonBond){
                Vec3d fij=Vec3dZero;
                //Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);
                Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                //Vec3d dp   = hj.f*(1./hj.w) - hi.f*(1./hi.w);
                //Vec3d dp   = apbc[j] - apbc[i];
                double Evdwi = getLJQH( dp, fij, REQij, R2damp );
                //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                //bErr|=ckeckNaN( 1,3, (double*)&fij, [&]{ printf("atom[%i]fLJ2[%i,%i]",ia,i,j); } );
                if(bPrint)printf( "ffl:vdw[%i|%i,%i] REQH{%g,%g,%g,%g} r %g fij(%g,%g,%g) E %g\n", ia,ing,jng, REQij.x,REQij.y,REQij.z,REQij.w,  dp.norm(), fij.x,fij.y,fij.z, Evdwi );
                f1.sub(fij);
                f2.add(fij);
                E-=Evdwi;
            }
            //if(bPrint)printf( "ffl:ang-vdw[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            fbs[i].add( f1     );
            fbs[j].add( f2     );
            //if(ia==ia_DBG)printf( "ffl:ANG[%i|%i,%i] fa(%g,%g,%g) fbs[%i](%g,%g,%g) fbs[%i](%g,%g,%g)\n", ia,ing,jng, fa.x,fa.y,fa.z, i,fbs[i].x,fbs[i].y,fbs[i].z,   j,fbs[j].x,fbs[j].y,fbs[j].z  );
            // ToDo: subtract non-covalent interactions
        }
    }
    }    
    //if(bErr){ printf("ERROR in ffl.eval_atom[%i] => Exit() \n", ia ); exit(0); }
    
    /*
    double Kfix = constr[ia].w;
    if(Kfix>0){  
        //printf( "applyConstrain(i=%i,K=%g)\n", ia, Kfix );
        Vec3d d = constr[ia].f-pa;
        d.mul( Kfix );
        double fr2 = d.norm2();
        double F2max = 1.0;
        if( fr2>F2max ){
            d.mul( sqrt(F2max/fr2) );
        }
        fa.add( d );
        E += d.norm()*Kfix*0.5;
    }
    */
    
    //fapos [ia].add(fa ); 
    //fpipos[ia].add(fpi);
    fapos [ia]=fa; 
    fpipos[ia]=fpi;

    //printf( "MYOUTPUT atom %i Ebond= %g Eangle= %g Edihed= %g Eimpr = %g Etot=%g\n", ia+1, Eb, Ea, EppI, Eps, E );
    //fprintf( file, "atom %i Ebond= %g Eangle= %g Edihed= %g Eimpr= %g Etot=%g \n", ia+1, Eb, Ea, EppI, Eps, E );
    

    return E;
}


};


#endif
