#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

//#include "testUtils.h"
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"

//#include "raytrace.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "NBFF.h"

#include "Molecule.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"

//#include "QEq.h"
//#include "NBSRFF.h"
//#include "IO_utils.h"

Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;
int i_DEBUG=0;

// ============ Global Variables

Molecule    mol;
MMFFparams  params;
MM::Builder builder;

MMFFsp3      ff;
NBFF        nff;
DynamicOpt  opt;


bool doBonded = false;

bool bNonBonded  = false;
bool bOptimizer  = true; 
bool bPBC        = false;
bool bCheckInvariants = true;
Vec3d cog,vcog,fcog,tqcog;

double Etot=0;
double  maxVcog = 1e-9;
double  maxFcog = 1e-9;
double  maxTg   = 1e-1;

//============================

#include "libUtils.h"

extern "C"{

void init_buffers(){
    /*
    double * DOFs  = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces
    double Kpi = 0.1;
    bool doPi=1;
    int  ipi0=0;
    int   * atype=0;
    Vec3d * apos=0;
    Vec3d * fapos=0;
    Vec3d * pipos=0;
    Vec3d * fpipos=0;
    //Vec3d * cappos=0;
    //Vec3d * fcappos=0;
    // --- Parameters
    Vec2i  * bond2atom = 0;
    double * bond_l0   = 0;  // [A]
    double * bond_k    = 0;  // [eV/A] ?
    Vec3d  * pbcShifts = 0;  // [A]
    int*    aneighs = 0;  // [natom*nneigh_max] neigbors of atom
    double* Kneighs = 0;  // [natom*nneigh_max] angular stiffness for given neighbor
    */
    ibuffers.insert( { "ndims",    &ff.nDOFs } );
    buffers.insert( { "Es",        &ff.Etot  } );
    buffers.insert( { "DOFs",    ff.DOFs } );
    buffers.insert( { "fDOFs",   ff.fDOFs } );
    buffers.insert( { "apos",   (double*)ff.apos   } );
    buffers.insert( { "fapos",  (double*)ff.fapos } );
    buffers.insert( { "pipos",  (double*)ff.pipos   } );
    buffers.insert( { "fpipos", (double*)ff.fpipos } );
    buffers.insert( { "bond_l0",   (double*)ff.bond_l0   } );
    buffers.insert( { "bond_k",    (double*)ff.bond_k    } );
    buffers.insert( { "pbcShifts", (double*)ff.pbcShifts } );
    buffers.insert( { "Kneighs",   (double*)ff.Kneighs   } );
    ibuffers.insert( { "bond2atom",    (int*)ff.bond2atom  } );
    ibuffers.insert( { "aneighs",      (int*)ff.aneighs  } );
}

void init_params(char* fatomtypes, char* fbondtypes){
    //params.loadAtomTypes( "AtomTypes.dat" );
    //params.loadBondTypes( "BondTypes.dat" );
    params.loadAtomTypes( fatomtypes );
    params.loadBondTypes( fbondtypes );
    builder.bindParams(&params);
}

void init_nonbond(){
    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.fapos, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );
    if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
        printf( "ERROR: nff.pairMask is not sorted => exit \n" );
        exit(0);
    };
}

void buildFF( bool bNonBonded_, bool bOptimizer_ ){
    bOptimizer=bOptimizer_;
    bNonBonded=bNonBonded_;
    builder.autoBonds();
    builder.autoAngles( 10.0, 10.0 );
    builder.sortConfAtomsFirst();
    builder.makeAllConfsSP();
    builder.assignAllBondParams();
    builder.toMMFFsp3( ff );
    if(bNonBonded)init_nonbond();
    if(bOptimizer){
        opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
        opt.cleanVel();
    }
    init_buffers();
}

int loadmol(char* fname_mol ){
    int imol = builder.loadMolType( fname_mol, "molecule" );
    builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
    return imol;
}

void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ ){
    init_params("data/AtomTypes.dat", "data/BondTypes.dat");
    loadmol( fname_mol );
    buildFF( bNonBonded_, bOptimizer_ );
}

void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( ff.doAngles , doAngles  );
    _setbool( ff.doPiPiT  , doPiPiT   );
    _setbool( ff.doPiSigma, doPiSigma );
    _setbool( ff.doPiPiI  , doPiPiI   );
    _setbool( doBonded  , doBonded_  );
    _setbool( bPBC        , PBC );
    _setbool( bPBC        , CheckInvariants  );
    #undef _setbool
}

bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){
    cog   = average( ff.natoms, ff.apos  );
    vcog  = sum    ( ff.natoms, (Vec3d*)opt.vel  );
    fcog  = sum    ( ff.natoms, ff.fapos );
    tqcog = torq   ( ff.natoms, ff.apos, ff.fapos, cog );
    //tqcog.add( ff.evalPiTorq() );
    return ( vcog.norm()>maxVcog ) || ( fcog.norm()>maxFcog ) || ( tqcog.norm() );
}

double eval(){
    double E=0;
    E+=ff.eval();
    if(bNonBonded)E+= nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
    return E;
}

bool relax( int niter, double Ftol ){
    Etot=0.0;
    double f2tol=Ftol*Ftol;
    bool bConverged=false; 
    for(int itr=0; itr<niter; itr++){
        Etot=eval();
        if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
        double f2 = opt.move_FIRE();
        if(f2<f2tol){ bConverged=true; break; }
    }
    return bConverged;
}

// ========= Manipulation with the molecule

void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      ){ rotate( n, selection, ff.apos, *(Vec3d*)p0, *(Vec3d*)ax, phi ); };
void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ rotate( n, selection, ff.apos, ff.apos[ia0], (ff.apos[iax1]-ff.apos[iax0]).normalized(), phi ); };

} // extern "C"
