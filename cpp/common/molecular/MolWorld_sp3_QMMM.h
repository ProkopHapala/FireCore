
#ifndef MolWorld_sp3_QMMM_h
#define MolWorld_sp3_QMMM_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <omp.h>

#include "IO_utils.h"

//#include "testUtils.h"
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"

#include "MMFFparams.h"

//#include "raytrace.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "MMFFsp3_loc.h"
#include "MMFFf4.h"

#include "NBFF.h"
#include "GridFF.h"
#include "RigidBodyFF.h"
#include "QEq.h"
#include "constrains.h"
#include "molecular_utils.h"

#include "Molecule.h"
#include "MMFFBuilder.h"
#include "SMILESparser.h"
#include "DynamicOpt.h"

#include "MultiSolverInterface.h"
#include "GlobalOptimizer.h"

#include "datatypes_utils.h"

#include "MolWorld_sp3.h"
#include "MarchingCubes.h"
#include "FireCoreAPI.h"

class MolWorld_sp3_QMMM : public MolWorld_sp3 { public:
    
    FireCore::Lib  fireCore;
    FireCore::QMMM qmmm;
    bool bPrepared_mm = false;
    bool bPrepared_qm = false;

    //void   InitQMMM();
    //virtual int projectOrbital(int iMO, double*& ewfaux ) override;  
    //virtual int projectDensity( double*& ewfaux ) override;

    // ==== Functions Definition

virtual void init( bool bGrid ){
    MolWorld_sp3::init( bGrid );
    if( file_exist("Fdata/info.dat") ){ 
        InitQMMM(); bPrepared_qm=true; 
        printf("QM preparation DONE \n");
    }

}


void InitQMMM(){
    printf( "TestAppFireCoreVisual::InitQMMM()\n" );
    // ----- QMMM setup
    if( file_exist("QM.dat") ){ 
        qmmm.init_file( "QM.dat" );
    }else{
        qmmm.init(ffl.natoms);
        qmmm.params=&params;
        //_vec2arr(qmmm.imms,  {4,5,10,11,    8,23} ); // map atom indexes from QM-system to MM-system 
        //_vec2arr(qmmm.isCap, {0,0, 0, 0,    1, 1} ); // which atoms are caps (not to be optimized)
        for(int i=0; i<ffl.natoms; i++){ qmmm.imms[i]=i; qmmm.isCap[i]=0; }
    }
    //ff.reallocMask();                           
    //ff.bondMasked[0] = true; printf( "bondMasked[0] %i,%i \n", ff.bond2atom[0].a,ff.bond2atom[0].b );
    //qmmm.maskMMFF(ff);
    qmmm.setAtypes( ffl.atypes );
    qmmm.load_apos( ffl.apos   );

    // ----- FireCore setup

    //fireCore.loadLib( "/home/prokop/git/FireCore/build/libFireCore.so" );
    //fireCore.loadLib( "/home/prokophapala/git/FireCore/build/libFireCore.so" );
    fireCore.loadLib( "./libFireCore.so" );
    qmmm.nmax_scf = 100;
    //qmmm.nmax_scf = 1;

    bool bInitFromDir = false;
    //bool bInitFromDir = true;

    if(bInitFromDir){
        //fireCore.setVerbosity(1,1);
        fireCore.initdir( );
        //fireCore.setVerbosity(2,1);
    }else{
        //fireCore.setVerbosity(1,1);
        fireCore.preinit( );
        //fireCore.setVerbosity(2,1);
        fireCore.set_lvs( (double*)&(builder.lvec) );
        fireCore.init( qmmm.nqm, qmmm.atypeZ, (double*)qmmm.apos );
    }

    double tmp[3]{0.,0.,0.};
    fireCore.setupGrid( 100.0, 0, tmp, (int*)&MOgrid.n, (double*)&MOgrid.dCell );
    MOgrid.updateCell_2();
    printf("MOgrid.printCell()\n");MOgrid.printCell();
    qmmm.bindFireCoreLib( fireCore );
    printf( "TestAppFireCoreVisual::InitQMMM() DONE !!! \n" );
}

virtual int projectOrbital(int iMO, double*& ewfaux ) override{  
    printf( "TestAppFireCoreVisual::projectOrbital() \n" );
    printf( "projectDensity() ffl.natoms %i qmmm.nqm %i \n", ffl.natoms, qmmm.nqm );
    qmmm.evalQM( ffl.apos, ffl.fapos );
    int ntot = MOgrid.n.x*MOgrid.n.y*MOgrid.n.z;
    ewfaux = new double[ ntot ];
    fireCore.getGridMO( iMO, ewfaux );
    return ntot;

    /*
    ogl_MO  = glGenLists(1);
    Vec3d p=Vec3d{0.4,2.5,0.0};
    glNewList(ogl_MO, GL_COMPILE);
    glTranslatef( p.x, p.y, p.z );
    int ntris=0;  
    glColor3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( MOgrid,  iso, ewfaux, isoSurfRenderType);
    glColor3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( MOgrid, -iso, ewfaux, isoSurfRenderType);
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawTriclinicBox(builder.lvec.transposed(), Vec3dZero, Vec3dOne );
    glTranslatef( -p.x, -p.y, -p.z );
    glEndList();
    delete [] ewfaux;
    */
}

virtual  int projectDensity( double*& ewfaux )override{
    printf( "TestAppFireCoreVisual::projectDensity() \n" );
    printf( "projectDensity() ffl.natoms %i qmmm.nqm %i \n", ffl.natoms, qmmm.nqm );
    qmmm.evalQM( ffl.apos, ffl.fapos );
    int ntot = MOgrid.n.x*MOgrid.n.y*MOgrid.n.z;
    ewfaux = new double[ ntot ];
    fireCore.getGridDens( 0, 0, ewfaux );
    return ntot;

    /*
    ogl_MO  = glGenLists(1);
    glNewList(ogl_MO, GL_COMPILE);
    int ntris = Draw3D::MarchingCubesCross( MOgrid, iso, ewfaux, isoSurfRenderType  );
    //printf( "renderOrbital() ntris %i \n", ntris );
    glEndList();
    delete [] ewfaux;
    */
}

virtual int getHOMO() override { return qmmm.nelec/2; };

};

#endif
