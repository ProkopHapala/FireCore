
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

int verbosity = 1;
int idebug    = 0;
double tick2second=1e-9;

#include "MolWorld_sp3_simple.h"

// ============ Global Variables

MolWorld_sp3_simple W;

//============================

//const std::vector<std::string>* atomTypeNames = 0;
#include "libUtils.h"

extern "C"{

void setVerbosity( int verbosity_, int idebug_ ){
    verbosity = verbosity_;
    idebug    = idebug_;
}

void init_buffers(){
    buffers .insert( { "apos",   (double*)W.nbmol.apos  } );
    buffers .insert( { "fapos",  (double*)W.nbmol.fapos } );
    if(W.bMMFF){
        buffers .insert( { "DOFs",      W.ffl.DOFs  } );
        buffers .insert( { "fDOFs",     W.ffl.fDOFs } );
        buffers .insert( { "pipos",  (double*)W.ffl.pipos   } );
        buffers .insert( { "fpipos", (double*)W.ffl.fpipos } );
        ibuffers.insert( { "neighs",   (int*)W.ffl.neighs  } );
    }
    ibuffers.insert( { "ndims",    &W.ffl.nDOFs } );
    buffers .insert( { "Es",       &W.ffl.Etot  } );
    //ibuffers.insert( { "selection", W.manipulation_sel  } );
    //for( auto c : buffers ){ printf("buff>>%s<<\n", c.first.c_str() ); }
}

void* init( char* xyz_name, char* smile_name, int* nPBC, char* sAtomTypes, char* sBondTypes, char* sAngleTypes ){
	//printf( "DEBUG nPBC(%i,%i,%i)\n", nPBC[0],nPBC[1],nPBC[2] );
    W.smile_name = smile_name;
	W.xyz_name   = xyz_name;
    W.nPBC       = *(Vec3i*)nPBC;
    W.tmpstr=tmpstr;
    W.params.init( sAtomTypes, sBondTypes, sAngleTypes );
	W.builder.bindParams(&W.params);
    W.init();
    //init_buffers();
    return &W;
}

void initParams       ( const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes ){ W.tmpstr=tmpstr; W.initParams(sAtomTypes,sBondTypes,sAngleTypes ); }
int  buildMolecule_xyz( const char* xyz_name, bool bEpairs, double fAutoCharges ){  W.builder.bDummyEpair=bEpairs; W.fAutoCharges=fAutoCharges;  return W.buildMolecule_xyz( xyz_name );  }
void makeFFs          ( ){ W.makeFFs(); }
void clear            ( ){ W.clear();   }

void setTrjName( char* trj_fname_, int savePerNsteps_ ){ W.trj_fname=trj_fname_; W.savePerNsteps=savePerNsteps_; if(verbosity>0)printf( "setTrjName(%s)\n", W.trj_fname ); }

void insertSMILES( char* s ){  W.insertSMILES(s); };

void initWithSMILES(char* s, bool bPrint, bool bCap, bool bNonBonded_, bool bOptimizer_ ){
    initWithSMILES(s,bPrint,bCap,bNonBonded_,bOptimizer_ );
    init_buffers();
}

const char* getType( int ia, bool fromFF){
    int it;
    if(fromFF){ if(ia>=W.nbmol.natoms)                return "?"; it = W.nbmol.atypes[ia];       }
    else      { if(ia>=W.builder.atoms.size())   return "?"; it = W.builder.atoms[ia].type; }
    if( (it<0) || (it>=W.params.atypes.size()) ) return "?";
    return W.params.atypes[it].name;
}

void setSwitches( int doAngles, int doPiPiT, int  doPiSigma, int doPiPiI, int doBonded_, int PBC, int CheckInvariants ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.ffl.doAngles , doAngles  );
    _setbool( W.ffl.doPiPiT  , doPiPiT   );
    _setbool( W.ffl.doPiSigma, doPiSigma );
    _setbool( W.ffl.doPiPiI  , doPiPiI   );
    _setbool( W.doBonded    , doBonded_  );
    _setbool( W.bPBC        , PBC );
    _setbool( W.bCheckInvariants, CheckInvariants  );
    #undef _setbool
}

void setOptLog( int n, double* cos, double* f, double* v, double* dt, double* damp ){
    W.opt_log.n    = n;
    W.opt_log.cos  = cos;
    W.opt_log.f    = f;
    W.opt_log.v    = v;
    W.opt_log.dt   = dt;
    W.opt_log.damp = damp;
}

bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){ return W.checkInvariants( maxVcog, maxFcog, maxTg ); }
//void open_xyzFile (const char* fname){ xyz_file=fopen( fname,"w" ); };
//void close_xyzFile(){fclose(xyz_file)};


int saveXYZ( const char* fname, const char* comment, int imod){ 
    int ret=-1;
    switch (imod){
        case 0: { ret=W.builder.save2xyz(fname,comment ); } break;
        case 1: { ret=W        .saveXYZ (fname,comment ); } break;
    }
    return ret;  
}


int toXYZ(const char* comment="#comment"){ return W.toXYZ(comment);  }
double eval(){ return W.eval(); };
int run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF ){ return W.run(nstepMax,dt,Fconv,ialg,outE,outF);  }

// ========= Manipulation with the molecule

void shift_atoms_ax ( int n, int* selection, double* d                               ){ W.shift_atoms ( n, selection,*(Vec3d*)d); };
void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      ){ W.rotate_atoms( n, selection,*(Vec3d*)p0, *(Vec3d*)ax, phi ); };
void shift_atoms    ( int n, int* selection, int ia0, int ia1, double l              ){ W.shift_atoms ( n, selection, ia0, ia1, l ); };
void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ W.rotate_atoms( n, selection, ia0, iax0, iax1, phi ); };

//int splitAtBond( int ib, int* selection ){ return W.splitAtBond( ib, selection ); }

void sampleNonBond(int n, double* rs, double* Es, double* fs, int kind, double*REQi_,double*REQj_, double K, double Rdamp ){
    Quat4d REQi = *(Quat4d*)REQi_;
    Quat4d REQj = *(Quat4d*)REQj_;
    Quat4d REQij; combineREQ( REQi, REQj, REQij );
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

int selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ, bool bOnlyFirstNeigh, int* atoms_ ){
    W.builder.selectBondsBetweenTypes( imin, imax, it1, it2, byZ, bOnlyFirstNeigh );
    Vec2i* atoms = (Vec2i*)atoms_;
    int i=0;
    for( int ib : W.builder.selection ){
        Vec2i b = W.builder.bonds[ib].atoms;
        if(b.b==it1){ b=b.swaped(); };
        //printf( "selectBondsBetweenTypes[%i] %i,%i \n", i, b.a,b.b );
        atoms[i]=b;
        i++;
    }
    return i;
}

int getFrament( int ifrag, int* bounds_, double* pose ){
    const MM::Fragment& frag = W.builder.frags[ifrag];
    if(bounds_){
        Vec2i* bounds=(Vec2i*)bounds_;
        bounds[0] = frag.atomRange;
        bounds[1] = frag.confRange;
        bounds[2] = frag.bondRange;
    }
    if( pose ){
        const double* ds = &(frag.pos.x);
        for(int i=0; i<7; i++){ pose[i] = ds[i]; }
    }
    return frag.imolType;
    //angRange;
    //dihRange;
}

void orient( const char* fname, int fw1,int fw2,  int up1,int up2,  int i0,  int imin, int imax ){
    FILE* fout = fopen(fname,"w");    
    Vec3d fw = W.builder.atoms[fw2].pos - W.builder.atoms[fw1].pos;  fw.normalize();
    Vec3d up = W.builder.atoms[up2].pos - W.builder.atoms[up1].pos;  up.normalize();
    W.builder.orient_atoms( fw, up, W.builder.atoms[i0].pos,  Vec3dZero,   imin,imax );
    W.builder.write2xyz(fout, "#scanHBond[0]" );
    fclose(fout);
}

void scanHBond( const char* fname, int n, double d,  int ifrag1, int ifrag2, int i1a,int i1b, int i2a,int i2b ){
    FILE* fout = fopen(fname,"w");
    const MM::Fragment& frag1 = W.builder.frags[ifrag1];
    const MM::Fragment& frag2 = W.builder.frags[ifrag2];
    
    Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    Vec3d fw2 = W.builder.atoms[i2b].pos - W.builder.atoms[i2a].pos;  fw2.normalize();

    W.builder.orient_atoms( fw1, Vec3dZ, W.builder.atoms[i1a].pos,  Vec3dZero,             frag1.atomRange.a, frag1.atomRange.b );

    W.builder.orient_atoms( fw2, Vec3dZ, W.builder.atoms[i2a].pos,  Vec3d{ 0.0,0.0,5.0 },  frag2.atomRange.a, frag2.atomRange.b );
    //W.builder.move();

    W.builder.write2xyz(fout, "#scanHBond[0]" );

    fclose(fout);
}


void printTypes     ( ){ W.params.printAtomTypes(); }
void printAtomConfs ( bool bOmmitCap, bool bNeighs ){ W.builder.printAtomConfs(bOmmitCap,bNeighs); }
void printAtomTypes ( ){ W.builder.printAtomTypes ( ); }
void printBonds     ( ){ W.builder.printBonds     ( ); }
void printBondParams( ){ W.builder.printBondParams( ); }

} // extern "C"
