
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

#include "libMMFF.h"
#include "libUtils.h"

extern "C"{

// void setVerbosity( int verbosity_, int idebug_ ){
//     verbosity = verbosity_;
//     idebug    = idebug_;
// }

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
    //printf( "MMFFsp3_lib::init_buffers() nDOFs=%i nnode=%i ncap=%i nvecs=%i \n", W.ffl.nDOFs, W.ffl.nnode, W.ffl.ncap, W.ffl.nvecs );
    //printf( "MMFFsp3_lib::init_buffers() nDOFs=%i nnode=%i ncap=%i nvecs=%i \n", (&W.ffl.nDOFs)[0], (&W.ffl.nDOFs)[1], (&W.ffl.nDOFs)[2], (&W.ffl.nDOFs)[3] );
    //ibuffers.insert( { "selection", W.manipulation_sel  } );
    //for( auto c : buffers ){ printf("buff>>%s<<\n", c.first.c_str() ); }
}

void* init( const char* xyz_name, const char* smile_name, int* nPBC,  const char* sElementTypes, const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes ){
	//printf( "DEBUG nPBC(%i,%i,%i)\n", nPBC[0],nPBC[1],nPBC[2] );
    W.smile_name = smile_name;
	W.xyz_name   = xyz_name;
    W.nPBC       = *(Vec3i*)nPBC;
    W.tmpstr=tmpstr;
    W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes );
	W.builder.bindParams(&W.params);
    W.init();
    //init_buffers();
    return &W;
}

// void initParams       ( const char* sElementTypes, const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes ){ W.tmpstr=tmpstr; W.initParams(sElementTypes,sAtomTypes,sBondTypes,sAngleTypes ); }
// int  buildMolecule_xyz( const char* xyz_name, bool bEpairs, double fAutoCharges ){  W.builder.bDummyEpair=bEpairs; W.fAutoCharges=fAutoCharges;  return W.buildMolecule_xyz( xyz_name );  }
// void makeFFs          ( ){ W.makeFFs(); }
// void clear            ( ){ W.clear();   }
// void setTrjName( const char* trj_fname_, int savePerNsteps_ ){ W.trj_fname=trj_fname_; W.savePerNsteps=savePerNsteps_; if(verbosity>0)printf( "setTrjName(%s)\n", W.trj_fname ); }
// void insertSMILES( const char* s ){  W.insertSMILES(s); };
// void initWithSMILES(const char* s, bool bPrint, bool bCap, bool bNonBonded_, bool bOptimizer_ ){
//     initWithSMILES(s,bPrint,bCap,bNonBonded_,bOptimizer_ );
//     init_buffers();
// }

// const char* getTypeName( int ia, bool fromFF){
//     int it;
//     if(fromFF){ if(ia>=W.nbmol.natoms)           return "?"; it = W.nbmol.atypes[ia];       }
//     else      { if(ia>=W.builder.atoms.size())   return "?"; it = W.builder.atoms[ia].type; }
//     if( (it<0) || (it>=W.params.atypes.size()) ) return "?";
//     return W.params.atypes[it].name;
// }

// void setSwitches( int CheckInvariants, int PBC, int NonBonded, int MMFF, int Angles, int PiSigma, int PiPiI ){
//     #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
//     _setbool( W.bCheckInvariants, CheckInvariants  );
//     _setbool( W.bPBC         , PBC       );
//     _setbool( W.bNonBonded   , NonBonded );
//     _setbool( W.bMMFF        , MMFF      );
//     _setbool( W.ffl.doAngles , Angles    );
//     _setbool( W.ffl.doPiSigma, PiSigma   );
//     _setbool( W.ffl.doPiPiI  , PiPiI     );
//     W.ffl.bSubtractAngleNonBond = W.bNonBonded;
//     #undef _setbool
// }

// void setOptLog( int n, double* cos, double* f, double* v, double* dt, double* damp ){
//     W.opt_log.n    = n;
//     W.opt_log.cos  = cos;
//     W.opt_log.f    = f;
//     W.opt_log.v    = v;
//     W.opt_log.dt   = dt;
//     W.opt_log.damp = damp;
// }

// bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){ return W.checkInvariants( maxVcog, maxFcog, maxTg ); }
//void open_xyzFile (const char* fname){ xyz_file=fopen( fname,"w" ); };
//void close_xyzFile(){fclose(xyz_file)};

// int toXYZ(const char* comment="#comment"){ return W.toXYZ(comment);  }
// double eval(){ return W.eval(); };
// int run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF ){ return W.run(nstepMax,dt,Fconv,ialg,outE,outF);  }

// ========= Manipulation with the molecule

// void shift_atoms_ax ( int n, int* selection, double* d                               ){ W.shift_atoms ( n, selection,*(Vec3d*)d); };
// void rotate_atoms_ax( int n, int* selection, double* p0, double* ax, double phi      ){ W.rotate_atoms( n, selection,*(Vec3d*)p0, *(Vec3d*)ax, phi ); };
// void shift_atoms    ( int n, int* selection, int ia0, int ia1, double l              ){ W.shift_atoms ( n, selection, ia0, ia1, l ); };
// void rotate_atoms   ( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ W.rotate_atoms( n, selection, ia0, iax0, iax1, phi ); };

// double measureBond        (int ia, int ib         ){ return W.ffl.measureBondLegth(ia,ib); }
// double measureAngle       (int ic, int ia, int ib ){ return W.ffl.measureAngle (ic,ia,ib); }
// double measureAnglePiPi   (int ia, int ib         ){ return W.ffl.measureAnglePiPi(ia, ib, true ); }
// double measureAngleSigmaPi(int ipi, int ia, int ib){ return W.ffl.measureAngleSigmaPi(ipi, ia, ib ); }

//int splitAtBond( int ib, int* selection ){ return W.splitAtBond( ib, selection ); }

// void sample_DistConstr( double lmin, double lmax, double kmin, double kmax, double flim , int n, double* xs, double* Es, double* Fs ){
//     DistConstr C( {0,1}, {lmax,lmin}, {kmax,kmin}, flim  );
//     Vec3d ps[2]{{.0,.0,.0},{.0,.0,.0}};
//     Vec3d fs[2];
//     for(int i=0; i<n; i++ ){
//         ps[1]={xs[i],0.0,0.0};
//         fs[0]=Vec3dZero;
//         fs[1]=Vec3dZero;
//         Es[i] = C.apply( ps, fs );
//         Fs[i] = fs[0].x;
//     }
// }

// void sample_evalPiAling( double k, double ang0, double r1, double r2, int n, double* angles, double* Es, double* Fs ){
//     Vec3d h1={1,0,0};
//     Vec3d f1,f2;
//     double c0 = cos(ang0);
//     for(int i=0; i<n; i++ ){
//         double a = angles[i];
//         Vec3d h2={cos(a),sin(a),0.0};
//         Es[i] = evalPiAling( h1, h2, 1./r1, 1./r2, k, f1, f2 );  
//         Fs[i] = f1.y;
//     }
// }

// void sample_evalAngleCos( double k, double ang0, double r1, double r2, int n, double* angles, double* Es, double* Fs ){
//     Vec3d h1={1,0,0};
//     Vec3d f1,f2;
//     double c0 = cos(ang0);
//     for(int i=0; i<n; i++ ){
//         double a = angles[i];
//         Vec3d h2={cos(a),sin(a),0.0};
//         Es[i] = evalAngleCos( h1, h2, 1./r1, 1./r2, k, c0, f1, f2 );
//         Fs[i] = f1.y;
//     }
// }

// void sample_evalAngleCosHalf( double k, double ang0, double r1, double r2, int n, double* angles, double* Es, double* Fs ){
//     Vec3d h1={1,0,0};
//     Vec3d f1,f2;
//     Vec2d cs0; cs0.fromAngle( ang0/2. );
//     for(int i=0; i<n; i++ ){
//         double a = angles[i];
//         Vec3d h2={cos(a),sin(a),0.0};
//         Es[i] = evalAngleCosHalf( h1, h2, 1./r1, 1./r2, cs0, k, f1, f2 );
//         Fs[i] = f1.y;
//     }
// }

// void sampleNonBond(int n, double* rs, double* Es, double* fs, int kind, double*REQi_,double*REQj_, double K, double Rdamp ){
//     Quat4d REQi = *(Quat4d*)REQi_;
//     Quat4d REQj = *(Quat4d*)REQj_;
//     Quat4d REQij; combineREQ( REQi, REQj, REQij );
//     REQij.y = sqrt(REQij.y);
//     Vec3d pi=Vec3dZero;
//     Vec3d pj=Vec3dZero;
//     double R2damp=Rdamp*Rdamp;
//     for(int i=0; i<n; i++){
//         double E;
//         Vec3d  f=Vec3dZero;
//         pj.x=rs[i];
//         switch(kind){
//             case 1: E=addAtomicForceMorseQ( pj-pi, f, REQij.x, REQij.y, REQij.z, K, R2damp );      break;  // Morse
//             case 2: E=addAtomicForceLJQ   ( pj-pi, f, REQij );                                              break;  // Lenard Jones
//             case 3: double fr; E=erfx_e6( pj.x, K, fr ); f.x=fr; break;  // gauss damped electrostatics
//         }
//         //printf( "i %i r %g E %g f %g \n", i, pj.x, E, f.x );
//         fs[i]=f.x;
//         Es[i]=E;
//     }
// }

// void scanTranslation_ax( int n, int* selection, double* vec, int nstep, double* Es, const char* trjName, bool bAddjustCaps ){
//     if(selection==0){ selection=W.manipulation_sel; n=W.manipulation_nsel; }
//     W.scanTranslation_ax( n, selection, *(Vec3d*)vec, nstep, Es, trjName, bAddjustCaps );
// }
// void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, const char* trjName, bool bAddjustCaps ){ 
//     W.scanTranslation( n, selection, ia0, ia1, l, nstep, Es, trjName, bAddjustCaps ); 
// }

// void scanRotation_ax( int n, int* selection, double* p0, double* ax, double phi, int nstep, double* Es, const char* trjName ){
//     if(p0==0) p0=(double*)&W.manipulation_p0;
//     if(ax==0) ax=(double*)&W.manipulation_ax;
//     if(selection==0){selection=W.manipulation_sel; n=W.manipulation_nsel; }
//     W.scanRotation_ax( n, selection, *(Vec3d*)p0, *(Vec3d*)ax, phi, nstep, Es, trjName );
// }
// void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, const char* trjName ){ 
//     W.scanRotation( n, selection,ia0, iax0, iax1, phi, nstep, Es, trjName );
// }

// void scanAngleToAxis_ax( int n, int* selection, double r, double R, double* p0, double* ax, int nstep, double* angs, double* Es, const char* trjName ){
//     W.scanAngleToAxis_ax( n, selection, r, R, *(Vec3d*)p0, *(Vec3d*)ax, nstep, angs, Es, trjName );
// }

// int selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ, bool bOnlyFirstNeigh, int* atoms_ ){
//     W.builder.selectBondsBetweenTypes( imin, imax, it1, it2, byZ, bOnlyFirstNeigh );
//     Vec2i* atoms = (Vec2i*)atoms_;
//     int i=0;
//     for( int ib : W.builder.selection ){
//         Vec2i b = W.builder.bonds[ib].atoms;
//         //int t1 = W.builder.atoms[b.a].type;
//         int t2 = W.builder.atoms[b.b].type;
//         //if( byZ ){ t1=W.params.atypes[t1].iZ; t2=W.params.atypes[t2].iZ; };
//         if( byZ ){ t2=W.params.atypes[t2].iZ; };
//         //printf( "b(%i,%i) types(%i,%i) its(%i,%i)\n",  b.a,b.b,  t1,t2, it1,it2 );
//         //if(t2==it1){ b=b.swaped(); printf( "bond swap b(%i,%i)\n", b.a, b.b ); };
//         if(t2==it1){ b=b.swaped(); };
//         //printf( "selectBondsBetweenTypes[%i] %i,%i \n", i, b.a,b.b );
//         atoms[i]=b;
//         i++;
//     }
//     return i;
// }

// int getFrament( int ifrag, int* bounds_, double* pose ){
//     const MM::Fragment& frag = W.builder.frags[ifrag];
//     if(bounds_){
//         Vec2i* bounds=(Vec2i*)bounds_;
//         bounds[0] = frag.atomRange;
//         bounds[1] = frag.confRange;
//         bounds[2] = frag.bondRange;
//     }
//     if( pose ){
//         const double* ds = &(frag.pos.x);
//         for(int i=0; i<7; i++){ pose[i] = ds[i]; }
//     }
//     return frag.imolType;
//     //angRange;
//     //dihRange;
// }

// void findMainAxes( double* rot, int ifrag, int imin, int imax, int* permut_, bool bRot){
//     Vec3i permut{2,1,0};
//     if(permut_){ permut=*((Vec3i*)permut_); };
//     if(ifrag>=0){ const MM::Fragment& frag = W.builder.frags[ifrag]; imin=frag.atomRange.a; imax=frag.atomRange.b; }
//     if(rot){ *((Mat3d*)rot)= W.builder.findMainAxes(imin,imax,true,bRot,permut); }
//     else   {                 W.builder.findMainAxes(imin,imax,true,bRot,permut); }    
// }

// void findSymmetry( int* found, int ifrag, int imin=0,int imax=-1, double tol=0.1 ){
//     if(ifrag>=0){ const MM::Fragment& frag = W.builder.frags[ifrag]; imin=frag.atomRange.a; imax=frag.atomRange.b; }
//     W.builder.findSymmetry( found, imin,imax, tol );
// }

// void orient( const char* fname, int fw1,int fw2,  int up1,int up2,  int i0,  int imin, int imax ){
//     FILE* fout = fopen(fname,"w");    
//     Vec3d fw = W.builder.atoms[fw2].pos - W.builder.atoms[fw1].pos;  fw.normalize();
//     Vec3d up = W.builder.atoms[up2].pos - W.builder.atoms[up1].pos;  up.normalize();
//     W.builder.orient_atoms( fw, up, W.builder.atoms[i0].pos,  Vec3dZero,   imin,imax );
//     W.builder.write2xyz(fout, "#scanHBond[0]" );
//     fclose(fout);
// }

void scanHBond( const char* fname, int n, double dl, double l0, int ifrag1, int ifrag2, int i1a,int i1b, int i2a,int i2b, bool isDonor1, bool isDonor2, double* ups_ ){
    FILE* fout = fopen(fname,"w");
    const MM::Fragment& frag1 = W.builder.frags[ifrag1];
    const MM::Fragment& frag2 = W.builder.frags[ifrag2];
    

    //Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    //Vec3d fw2 = W.builder.atoms[i2b].pos - W.builder.atoms[i2a].pos;  fw2.normalize();    
    // W.builder.orient_atoms( fw1, up1, W.builder.atoms[i1a].pos,  Vec3dZero,            frag1.atomRange.a, frag1.atomRange.b );
    // W.builder.orient_atoms( fw2, up2, W.builder.atoms[i2a].pos,  Vec3d{ 0.0,0.0,l0 },  frag2.atomRange.a, frag2.atomRange.b );

    // Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    // Vec3d fw2 = W.builder.atoms[i2a].pos - W.builder.atoms[i2b].pos;  fw2.normalize();  // reverse direction    
    // W.builder.orient_atoms( fw1, up1, W.builder.atoms[i1b].pos,  Vec3dZero,            frag1.atomRange.a, frag1.atomRange.b ); // donor    => center is capping hydrogen [i1b]
    // W.builder.orient_atoms( fw2, up2, W.builder.atoms[i2a].pos,  Vec3d{ 0.0,0.0,l0 },  frag2.atomRange.a, frag2.atomRange.b ); // acceptor => center is node atom        [i2a]
    //printf( "isDonor %i %i \n", isDonor1, isDonor2 );
    Vec3d p1,p2;
    if(isDonor1){ p1=W.builder.atoms[i1b].pos; }else{ p1=W.builder.atoms[i1a].pos; }
    if(isDonor2){ p2=W.builder.atoms[i2b].pos; }else{ p2=W.builder.atoms[i2a].pos; }
    Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    Vec3d fw2 = W.builder.atoms[i2b].pos - W.builder.atoms[i2a].pos;  fw2.normalize();   fw2.mul(-1.); // reverse direction   

    Vec3d up1,up2;
    if(ups_){ Vec3d* ups=(Vec3d*)ups_; up1=ups[0]; up2=ups[1]; up1.normalize(); up2.normalize(); }
    else{  // auto-up-vector (make sure it is not colinear with fw)
        if( fabs(fabs(fw1.z)-1.0) < 0.1 ){ up1=Vec3dY; }else{ up1=Vec3dZ; } 
        if( fabs(fabs(fw2.z)-1.0) < 0.1 ){ up2=Vec3dY; }else{ up2=Vec3dZ; }
    }
    printf( "fw1(%g,%g,%g) |fw1.z-1|=%g up1(%g,%g,%g) p1(%g,%g,%g)\n", fw1.x,fw1.y,fw1.z, fabs(fabs(fw1.z)-1.0),   up1.x,up1.y,up1.z,    p1.x,p1.y,p1.z );
    printf( "fw2(%g,%g,%g) |fw2.z-1|=%g up2(%g,%g,%g) p2(%g,%g,%g)\n", fw2.x,fw2.y,fw2.z, fabs(fabs(fw2.z)-1.0),   up2.x,up2.y,up2.z,    p2.x,p2.y,p2.z );

    W.builder.orient_atoms( fw1, up1, p1,  Vec3d{ 0.0,0.0,0  }, frag1.atomRange.a, frag1.atomRange.b ); // donor    => center is capping hydrogen [i1b]
    W.builder.orient_atoms( fw2, up2, p2,  Vec3d{ 0.0,0.0,l0 }, frag2.atomRange.a, frag2.atomRange.b ); // acceptor => center is node atom        [i2a]

    Vec3d d=Vec3d{ 0.0,0.0, dl };
    for(int i=0; i<n+1; i++){
        sprintf(tmpstr, "#scanHBond[%i] l=%7.3f[A]", i, l0+dl*i );
        W.builder.write2xyz(fout, tmpstr );
        W.builder.move_atoms( d, frag2.atomRange.a, frag2.atomRange.b );
    }
    fclose(fout);
}


// void printTypes     ( ){ W.params.printAtomTypes(true); }
// void printAtomConfs ( bool bOmmitCap, bool bNeighs ){ W.builder.printAtomConfs(bOmmitCap,bNeighs); }
// void printAtomTypes ( ){ W.builder.printAtomTypes ( ); }
// void printBonds     ( ){ W.builder.printBonds     ( ); }
// void printBondParams( ){ W.builder.printBondParams( ); }
// void printAtomParams( ){ W.ffl.printAtomParams( );     }
// void printSwitches  ( ){ W.printSwitches( );           }

} // extern "C"
