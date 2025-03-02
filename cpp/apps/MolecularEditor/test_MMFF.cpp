
#include <globals.h>
//int verbosity = 0;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "raytrace.h"
#include "Molecule.h"
#include "MMFF.h"
#include "MMFFBuilder.h"

#include "IO_utils.h"

#include "geom3D.h"

#include "DynamicOpt.h"

// ==== Global Variables

MMFFparams  params;
MMFF        world;
MM::Builder builder;
DynamicOpt  opt;

//int     fontTex;
//int     ogl_sph;

char str[256];

void   initRigidSubstrate();
void   bakeMMFF();
void   prepareOpt();
double relaxNsteps( int nsteps, double F2conf );

// ==== Function Implementation

void initRigidSubstrate(){
    // ---- Rigid Substrate
    printf( "params.atypNames:\n" );
    for(auto kv : params.atomTypeDict) { printf(" %s %i \n", kv.first.c_str(), kv.second ); }
    world.gridFF.grid.n    = Vec3i{60,60,100};
    world.gridFF.grid.pos0 = Vec3d{0.0,0.0,0.0};
    world.gridFF.loadCell ( "inputs/cel.lvs" );
    world.gridFF.grid.printCell();
    //world.gridFF.loadXYZ  ( "inputs/NaCl_wo4.xyz", params );
    params.loadXYZ( "inputs/NaCl_wo4.xyz", world.gridFF.natoms, &world.gridFF.apos, &world.gridFF.REQs, &world.gridFF.atypes );
    world.translate( {0.0,0.0,4.5} );

    Quat4d testREQ;
    Quat4f testPLQ;
    testREQ = Quat4d{ 1.487, sqrt(0.0006808), 0., 0.}; // H
    testPLQ = REQ2PLQ( testREQ, -1.6 );//
    world.genPLQ();
    world.gridFF.allocateFFs();
    bool recalcFF = false;
    if( recalcFF ){
        world.gridFF.makeGridFF();
        if(world.gridFF.FFelec) saveBin( "data/FFelec.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFelec );
        if(world.gridFF.FFPaul) saveBin( "data/FFPaul.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFPaul );
        if(world.gridFF.FFLond) saveBin( "data/FFLond.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFLond );
    }else{
        if(world.gridFF.FFelec) loadBin( "data/FFelec.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFelec );
        if(world.gridFF.FFPaul) loadBin( "data/FFPaul.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFPaul );
        if(world.gridFF.FFLond) loadBin( "data/FFLond.bin", world.gridFF.grid.getNtot()*sizeof(Vec3d), (char*)world.gridFF.FFLond );
    }
    int iatom = 11;
    printf( "testREQ   (%g,%g,%g) -> PLQ (%g,%g,%g) \n",        testREQ.x, testREQ.y, testREQ.z, testPLQ.x, testPLQ.y, testPLQ.z   );
    printf( "REQs[%i] (%g,%g,%g) -> PLQ (%g,%g,%g) \n", iatom, world.REQ[iatom].x, world.REQ[iatom].y, world.REQ[iatom].z, world.PLQ[iatom].x, world.PLQ[iatom].y, world.PLQ[iatom].z );
    Quat4f * FFtot = new Quat4f[world.gridFF.grid.getNtot()];
    world.gridFF.evalCombindGridFF( testREQ, FFtot );
    world.gridFF.grid.saveXSF<float>( "FFtot_z.xsf", (float*)FFtot, 4,3, world.gridFF.natoms, world.gridFF.atypes,  world.gridFF.apos );

    //isoOgl = glGenLists(1);
    //glNewList(isoOgl, GL_COMPILE);
    //    renderSubstrate_( world.gridFF.grid, FFtot, world.gridFF.FFelec, 0.01, true );
    //    Draw3D::drawAxis(1.0);
    //glEndList();
    //cam.pos.z = +5.0;
}

/*
void initParams( char* fname_atomTypes, char* fname_bondTypes ){
    builder.params = &params;
    if(fname_atomTypes) params.loadAtomTypes( fname_atomTypes );
    if(fname_bondTypes) params.loadBondTypes( fname_bondTypes );
    printf( "params.atypNames.size() %i \n", params.atypNames.size() );
    //printf("initParams done! \n");
}
*/

//int loadMolType   ( char* fname ){ return builder.loadMolType(fname ); };
//int insertMolecule( int itype, double* pos, double* rot, bool rigid ){ return builder.insertMolecule( itype, *(Vec3d*)pos, *(Mat3d*)rot, rigid ); };

void bakeMMFF(){
    builder.toMMFF( &world, &params );
    world.genPLQ();
    world.printAtomInfo(); //exit(0);
    //world.allocFragment( nFrag );
    //opt.bindArrays( 8*world.nFrag, (double*)world.poses, new double[8*world.nFrag], (double*)world.poseFs );
}

void prepareOpt(){
    //opt.bindArrays( 8*world.nFrag, world.poses, world.poseVs, world.poseFs );
    //printf("DEBUG a.0\n");
    world.allocateDyn();   
    world.initDyn();       
    opt.bindArrays( world.nDyn, world.dynPos, world.dynVel, world.dynForce, world.dynInvMass );
    opt.setInvMass( 1.0 ); 
    opt.cleanVel  ( );     
    //exit(0);
    //printf("POSE_pos   : \n"); printPoses( world.nFrag, world.poses  );
    //printf("POSE_Force : \n"); printPoses( world.nFrag, world.poseFs );
}

double relaxNsteps( int nsteps, double F2conf ){
    double F2=1e+300;
    for(int itr=0; itr<nsteps; itr++){
        //printf( "===== relaxNsteps itr %i \n", itr );
        world.cleanAtomForce();
        world.frags2atoms();
        if( world.gridFF.FFPaul ) world.eval_FFgrid();
        world.eval_MorseQ_On2_fragAware();

        world.cleanPoseTemps();
        world.aforce2frags();

        world.toDym(true);
        F2 = opt.move_FIRE();
        //printf( "F2 %g dt %g \n", F2, opt.dt );
        if(F2<F2conf) break;
        world.checkPoseUnitary();
        world.fromDym();
        printf( ">> itr %i F2 %g dt %g qrot (%g,%g,%g,%g) int %li \n", itr, F2, opt.dt, world.poses[4], world.poses[5], world.poses[6], world.poses[7], world.gridFF.FFPaul );
        //printf( ">> itr %i F2 %g dt %g poses (%g,%g,%g,%g, %g,%g,%g,%g) \n", itr, F2, world.poses[0], world.poses[1], world.poses[2], world.poses[3], world.poses[4], world.poses[5], world.poses[6], world.poses[7] );
    }
    return F2;
}

int main(){

    // ======= common Potentials etc.
    //builder.params = &params;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );

    Mat3d rot0; rot0.setOne();
    int itype=-1;

    printf( "// =========== System 1 \n" );
    builder.clear();
    itype = builder.loadMolTypeXYZ( "inputs/water_T5_ax.xyz", &params );
    builder.insertMolecule( itype, Vec3d{5.78, 6.7, 12.24}, rot0, true );
    bakeMMFF();
    prepareOpt();
    print("DEBUG prepareOpt() -> relaxNsteps \n");
    relaxNsteps( 3, 0.0 );

    // =========== System 2
    printf( "// =========== System 2 \n" );
    builder.clear();
    itype = builder.loadMolTypeXYZ( "inputs/Campher.xyz", &params );
    builder.insertMolecule( itype, Vec3d{5.78, 6.7, 12.24}, rot0, true );
    bakeMMFF();
    prepareOpt();
    relaxNsteps( 3, 0.0 );


    for(Molecule* m : builder.molTypes ){ m->dealloc(); delete m; };


    printf( "ALL DONE !\n" );
    exit(0);

}

