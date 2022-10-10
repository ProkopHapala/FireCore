
/*
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
#include "molecular_utils.h"

#include "Molecule.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "SMILESparser.h"
#include "DynamicOpt.h"

//#include "QEq.h"
//#include "NBSRFF.h"
//#include "IO_utils.h"

Vec3d DEBUG_dQdp;
int DEBUG_iter     = 0;
int DEBUG_log_iter = 0;
int i_DEBUG=0;
*/

constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

#include "MolWorld_sp3.h"

// ============ Global Variables

MolWorld_sp3 W;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
#include "libUtils.h"

extern "C"{

void init_buffers(){
    ibuffers.insert( { "ndims",    &W.ff.nDOFs } );
    buffers.insert( { "Es",        &W.ff.Etot  } );
    buffers.insert( { "DOFs",    W.ff.DOFs } );
    buffers.insert( { "fDOFs",   W.ff.fDOFs } );
    buffers.insert( { "apos",   (double*)W.ff.apos   } );
    buffers.insert( { "fapos",  (double*)W.ff.fapos } );
    buffers.insert( { "pipos",  (double*)W.ff.pipos   } );
    buffers.insert( { "fpipos", (double*)W.ff.fpipos } );
    buffers.insert( { "bond_l0",   (double*)W.ff.bond_l0   } );
    buffers.insert( { "bond_k",    (double*)W.ff.bond_k    } );
    buffers.insert( { "pbcShifts", (double*)W.ff.pbcShifts } );
    buffers.insert( { "Kneighs",   (double*)W.ff.Kneighs   } );
    ibuffers.insert( { "bond2atom",    (int*)W.ff.bond2atom  } );
    ibuffers.insert( { "aneighs",      (int*)W.ff.aneighs  } );
    ibuffers.insert( { "selection", W.manipulation_sel  } );
}

int loadmol(char* fname_mol ){ return W.loadmol(fname_mol ); }

void init(){
    W.tmpstr=tmpstr;
    //sprintf("%s/AtomTypes.dat");
    W.params.init("data/AtomTypes.dat", "data/BondTypes.dat", "data/AngleTypes.dat" );
}

void initWithMolFile(char* fname_mol, bool bNonBonded_, bool bOptimizer_ ){
    W.tmpstr=tmpstr;
    W.params.init("data/AtomTypes.dat", "data/BondTypes.dat", "data/AngleTypes.dat" );
    loadmol( fname_mol );
    W.buildFF( bNonBonded_, bOptimizer_ );
}

void insertSMILES( char* s ){  W.insertSMILES(s); };

void initWithSMILES(char* s, bool bPrint, bool bCap, bool bNonBonded_, bool bOptimizer_ ){
    initWithSMILES(s,bPrint,bCap,bNonBonded_,bOptimizer_ );
    init_buffers();
}


void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.ff.doAngles , doAngles  );
    _setbool( W.ff.doPiPiT  , doPiPiT   );
    _setbool( W.ff.doPiSigma, doPiSigma );
    _setbool( W.ff.doPiPiI  , doPiPiI   );
    _setbool( W.doBonded    , doBonded_  );
    _setbool( W.bPBC        , PBC );
    _setbool( W.bCheckInvariants, CheckInvariants  );
    #undef _setbool
}

bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){ return W.checkInvariants( maxVcog, maxFcog, maxTg ); }
//void open_xyzFile (const char* fname){ xyz_file=fopen( fname,"w" ); };
//void close_xyzFile(){fclose(xyz_file)};

int toXYZ(const char* comment="#comment"){ return W.toXYZ(comment); }

double eval(){ return W.eval(); };

bool relax( int niter, double Ftol, bool bWriteTrj ){
    return W.relax( niter, Ftol, bWriteTrj );
}

// ========= Manipulation with the molecule

void shift_atoms_ax ( int n, int* selection, double* d                               ){ W.shift_atoms ( n, selection,*(Vec3d*)d); };
void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      ){ W.rotate_atoms( n, selection,*(Vec3d*)p0, *(Vec3d*)ax, phi ); };
void shift_atoms    ( int n, int* selection, int ia0, int ia1, double l              ){ W.shift_atoms ( n, selection, ia0, ia1, l ); };
void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ W.rotate_atoms( n, selection, ia0, iax0, iax1, phi ); };

int splitAtBond( int ib, int* selection ){ return W.splitAtBond( ib, selection ); }

void sampleNonBond(int n, double* rs, double* Es, double* fs, int kind, double*REQi_,double*REQj_, double K, double Rdamp ){
    Vec3d REQi = *(Vec3d*)REQi_;
    Vec3d REQj = *(Vec3d*)REQj_;
    Vec3d REQij; combineREQ( REQi, REQj, REQij );
    REQij.y = sqrt(REQij.y);
    Vec3d pi=Vec3dZero;
    Vec3d pj=Vec3dZero;
    double R2damp=Rdamp*Rdamp;
    for(int i=0; i<n; i++){
        double E;
        Vec3d  f=Vec3dZero;
        pj.x=rs[i];
        switch(kind){
            case 1: E=addAtomicForceMorseQ( pj-pi, f, REQij.x, REQij.y, REQij.z, K, R2damp );      break;  // Morse
            case 2: E=addAtomicForceLJQ   ( pj-pi, f, REQij );                                              break;  // Lenard Jones
            case 3: double fr; E=erfx_e6( pj.x, K, fr ); f.x=fr; break;  // gauss damped electrostatics
        }
        //printf( "i %i r %g E %g f %g \n", i, pj.x, E, f.x );
        fs[i]=f.x;
        Es[i]=E;
    }
}

void sampleSurf(char* name, int n, double* rs, double* Es, double* fs, int kind, int atyp, double Q, double K, double Rdamp, double* pos0_, bool bInit, bool bSave){
    if(bInit){
        W.ff.realloc( 1,0,0,0, true );
        W.ff.apos [0] = *(Vec3d*)pos0_;
        W.ff.atype[0] = atyp;
        W.loadSurf( name, bSave );
        W.nbmol.REQs[0].z = Q;
        if(bSave){
            Quat4f* FFtot = new Quat4f[W.gridFF.grid.getNtot()];
            W.gridFF.evalCombindGridFF ( W.nbmol.REQs[0], FFtot );
            //void saveXSF( const char * fname, T* FF, int pitch, int offset, int natoms=0, int* atypes=0, Vec3d* apos=0  )const {
            //W.gridFF.grid.saveXSF<float>( "FFtot_E.xsf", (float*)FFtot, 4, 3, W.gridFF.natoms, W.gridFF.atypes, W.gridFF.apos );
            W.gridFF.grid.saveXSF<float>( "FFtot_E.xsf", (float*)FFtot, 4, 3 );
            printf("save DONE \n");
            delete [] FFtot;
        }
    }
    Vec3d PLQ = REQ2PLQ( W.nbmol.REQs[0], K );
    printf( "PLQ(%g,%g,%g) \n", PLQ.x, PLQ.y, PLQ.z );
    double R2damp=Rdamp*Rdamp;
    for(int i=0; i<n; i++){
        Quat4f fe=Quat4fZero;
        W.nbmol.ps[0].z=rs[i];
        W.ff.cleanAtomForce();
        switch(kind){
            case  1: fe.e=   W.nbmol.evalMorse(W.surf, false, K); fe.f=(Vec3f)W.nbmol.fs[0]; break; 
            case 11:         W.gridFF.addForce_surf( W.nbmol.ps[0], PLQ, fe );  break; // printf( "sampleSurf[%i] r %g fe(%g,%g,%g|%g)\n", i, rs[i], fe.x,fe.y,fe.z, fe.e ); 
            case 12:         W.gridFF.addForce     ( W.nbmol.ps[0], PLQ, fe );  break; // printf( "sampleSurf[%i] r %g fe(%g,%g,%g|%g)\n", i, rs[i], fe.x,fe.y,fe.z, fe.e ); 
            //case 3: double fr; E=erfx_e6( pj.x, K, fr ); f.x=fr; break;  // gauss damped electrostatics
        }
        fs[i]=fe.z;
        Es[i]=fe.e;
    }
}

void scanTranslation_ax( int n, int* selection, double* vec, int nstep, double* Es, bool bWriteTrj ){
    if(selection==0){ selection=W.manipulation_sel; n=W.manipulation_nsel; }
    W.scanTranslation_ax( n, selection, *(Vec3d*)vec, nstep, Es, bWriteTrj );
}
void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj ){ 
    W.scanTranslation( n, selection, ia0, ia1, l, nstep, Es, bWriteTrj ); 
}

void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, bool bWriteTrj ){
    if(p0==0) p0=(double*)&W.manipulation_p0;
    if(ax==0) ax=(double*)&W.manipulation_ax;
    if(selection==0){selection=W.manipulation_sel; n=W.manipulation_nsel; }
    W.scanRotation_ax( n, selection, *(Vec3d*)p0, *(Vec3d*)ax, phi, nstep, Es, bWriteTrj );
}
void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj ){ 
    W.scanRotation( n, selection,ia0, iax0, iax1, phi, nstep, Es, bWriteTrj );
}

} // extern "C"
