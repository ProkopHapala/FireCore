#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

//#include "testUtils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

//#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFsp3.h"
#include "NBFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

//#include "NBSRFF.h"
#include "IO_utils.h"


Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;
int i_DEBUG=0;

#include "InteractionsGauss.h"
#include "eFF.h"

// ============ Global Variables

Molecule    mol;
MMFFparams  params;
MM::Builder builder;

//NBFF       nff;
//DynamicOpt  opt;

MMFFsp3 ff;

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

void init( const char* fname_mol ){
    //ff.realloc( int nnode_, int nbonds_, int npi_, int ncap_, bool bNeighs=true );
    params.loadAtomTypes( "AtomTypes.dat" );
    params.loadBondTypes( "BondTypes.dat" );
    readMatrix( "polymer-2.lvs", 3, 3, (double*)&builder.lvec );
    builder.insertFlexibleMolecule( builder.loadMolType( fname_mol, "molecule" ), {0,0,0}, Mat3dIdentity, -1 );

    //if(verbosity>1)builder.printAtomConfs();
    //builder.export_atypes(atypes);
    builder.autoBondsPBC();             //if(verbosity>1)builder.printBonds ();  // exit(0);
    builder.autoAngles( 10.0, 10.0 );   //if(verbosity>1)builder.printAngles();
    builder.sortConfAtomsFirst();
    builder.makeAllConfsSP();
    //if(verbosity>1)builder.printAtomConfs();
    builder.toMMFFsp3( ff );

    init_buffers();
}

} // extern "C"
