
#ifndef MolWorld_sp3_h
#define MolWorld_sp3_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "IO_utils.h"

//#include "testUtils.h"
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"

#include "MMFFparams.h"
static MMFFparams* params_glob;

//#include "raytrace.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "MMFFsp3_loc.h"
#include "MMFFf4.h"

#include "NBFF.h"
#include "GridFF.h"
#include "RigidBodyFF.h"
#include "QEq.h"
#include "molecular_utils.h"

#include "Molecule.h"
#include "MMFFBuilder.h"
#include "SMILESparser.h"
#include "DynamicOpt.h"

#include "datatypes_utils.h"

class MolWorld_sp3_simple{ public:
    //const char* data_dir     = "common_resources";
    const char* xyz_name     = "input";
    const char* substitute_name = 0;    int isubs;
    const char* smile_name   = 0;
    const char* trj_fname    = "trj.xyz";
    int savePerNsteps = 1;
    double fAutoCharges=-1;

    OptLog opt_log;

	// Building
	MMFFparams   params;
	MM::Builder  builder;
	SMILESparser smiles;

	// Force-Fields & Dynamics
    MMFFsp3_loc  ffl;
	NBFF     nbmol;
    //QEq          qeq;
	DynamicOpt   opt;

    Vec3i nPBC{1,1,0};
    int    npbc       = 0;
    Vec3d* pbc_shifts = 0;
    double RdampCoul = 1.0;


	// state
	bool    bConverged = false;
	double  Etot=0;
	double  maxVcog = 1e-9;
	double  maxFcog = 1e-9;
	double  maxTg   = 1e-1;
	double  Kmorse = -1.0;

	// force-eval-switchefs
    int  imethod=0;
	bool doBonded    = false;
	bool bNonBonded  = true;
    bool bMMFF       = true;
	bool bOptimizer  = true; 
	bool bPBC        = false;
	bool bCheckInvariants = true;
	Vec3d cog,vcog,fcog,tqcog;
    int nloop=0;
    //bool bChargeUpdated=false;

	// Selecteion & Manipulasion
	std::vector<int> selection;
	Vec3d manipulation_p0=Vec3dZero; 
	Vec3d manipulation_ax=Vec3dZ;
	int*  manipulation_sel=0;
	int   manipulation_nsel=0;

    int ipicked    = -1; // picket atom 
    int ibpicked   = -1; // picket bond
    int iangPicked = -1; // picket angle
    Vec3d* picked_lvec = 0;
    Vec3d pick_hray, pick_ray0;

	// IO temp & aux
	FILE* xyz_file=0;
	char* tmpstr;


// =================== Functions


// ======================================
// ========= Initialization
// ======================================

int makePBCshifts( Vec3i nPBC, const Mat3d& lvec ){
    npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
    pbc_shifts = new Vec3d[npbc];
    int ipbc=0;
    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
        pbc_shifts[ipbc] = (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz);   
        //printf( "shifts[%3i=%2i,%2i,%2i] (%7.3f,%7.3f,%7.3f)\n",  ipbc, ix,iy,iz, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z );
        ipbc++; 
    }}}
    return npbc;
}

void initNBmol( int na, Vec3d* apos, Vec3d* fapos, int* atypes, bool bCleanCharge=true ){
    if(verbosity>0)printf( "MolWorld_sp3::initNBmol() na %i \n", na  );
	nbmol.bindOrRealloc( na, apos, fapos, 0, atypes );    
    //nbmol.bindOrRealloc( na, apos, fapos, 0, 0 );   
    //builder.export_atypes( nbmol.atypes );     
	builder.export_REQs( nbmol.REQs  );    
    if(bCleanCharge)for(int i=builder.atoms.size(); i<na; i++){ nbmol.REQs[i].z=0; }  // Make sure that atoms not present in Builder has well-defined chanrge                       
    params.assignREs( na, nbmol.atypes, nbmol.REQs, true, false  );
    if(verbosity>1)nbmol.print();                              
}

int loadGeom( const char* name ){ // TODO : overlaps with buildFF()
    if(verbosity>0)printf("MolWorld_sp3::loadGeom(%s)\n",  name );
    // ------ Load geometry
    sprintf(tmpstr, "%s.xyz", name ); 
    //printf("tmpstr=`%s`\n", tmpstr);
    int imol  = builder.loadMolType( tmpstr, name );
    int iret  = builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
    int ifrag = builder.frags.size()-1;
    if(iret<0){ printf("!!! exit(0) in MolWorld_sp3::loadGeom(%s)\n", name); exit(0); }
    builder.addCappingTypesByIz(1);
    builder.tryAddConfsToAtoms( 0, -1 );
    builder.cleanPis();
    if(verbosity>2)builder.printAtomConfs(false);
    //builder.export_atypes(atypes);
    // ------- Load lattice vectros
    sprintf(tmpstr, "%s.lvs", name );
    //builder.printAtomConfs(true);
    if( file_exist(tmpstr) ){
        builder.bPBC=true;
        readMatrix( tmpstr, 3, 3, (double*)&builder.lvec );
    }
    bPBC=builder.bPBC;
    if( bPBC ){
        builder.autoBondsPBC();   if(verbosity>2)builder.printBonds ();
    }else{
        builder.autoBonds();      if(verbosity>2)builder.printBonds ();
    }
    return ifrag;
}

void insertSMILES(const char* s){
    smiles.builder=&builder;
    smiles.parseString( 10000, s );
}

void setOptimizer( int n, double* ps, double* fs ){
    //opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
    opt.bindOrAlloc( n, ps, 0, fs, 0 );
    double dtopt=ffl.optimalTimeStep(); 
    if(verbosity>0)printf("MolWorld_sp3::setOptimizer(): optimnal time step = %g \n", dtopt);
    opt.initOpt( dtopt );
    opt.cleanVel();
    //opt.verbosity=2;
}
void setOptimizer(){ setOptimizer( ffl.nDOFs, ffl.DOFs, ffl.fDOFs ); };

void initParams( const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes ){
    params.init( sAtomTypes, sBondTypes, sAngleTypes );
    builder.bindParams(&params);
    params_glob = &params;
    builder.capAtomEpair.type = params.getAtomType("E");
    //params.printAtomTypeDict();
    //params.printAtomTypes();
    //params.printBond();
}

int buildMolecule_xyz( const char* xyz_name ){
    int ifrag = loadGeom( xyz_name );
    if( fAutoCharges>0 )builder.chargeByNeighbors( true, fAutoCharges, 10, 0.5 );
    //if(substitute_name) substituteMolecule( substitute_name, isubs, Vec3dZ );
    if( builder.checkNeighsRepeat( true ) ){ printf( "ERROR: some atoms has repating neighbors => exit() \n"); exit(0); };
    builder.autoAllConfEPi( );        
    builder.printAtomConfs(false, true );
    builder.assignAllBondParams();    //if(verbosity>1)
    builder.finishFragment(ifrag);    
    //printf( "buildMolecule_xyz: nMulPBC(%i,%i,%i) \n",nMulPBC.x,nMulPBC.y,nMulPBC.z  );
    //if( nMulPBC    .totprod()>1 ){ PBC_multiply    ( nMulPBC, ifrag ); };
    //if( bCellBySurf             ){ changeCellBySurf( bySurf_lat[0], bySurf_lat[1], bySurf_ia0, bySurf_c0 ); };
    //printf("builder.lvec\n");builder.lvec.print();
    return ifrag;
}

void makeMMFF(){
    //builder.printAtoms();          
    //if( builder.checkBondsOrdered( false, true ) ) { printf("ERROR Bonds are not ordered => exit"); exit(0); };
    if( builder.checkBondsInNeighs(true) ) { 
        printf("ERROR some bonds are not in atom neighbors => exit"); 
        exit(0); 
    };
    DEBUG
    builder.sortConfAtomsFirst();
    //builder.printAtomConfs(false,true);
    builder.checkBondsOrdered( true, false );
    DEBUG
    builder.assignTypes( );
    DEBUG
    builder.toMMFFsp3_loc( ffl, &params ); //ffl.printAtomParams(); ffl.printBKneighs(); 
    DEBUG
    ffl.setLvec(       builder.lvec);
    nPBC=Vec3i{0,0,0};
    DEBUG
    ffl.makeNeighCells( nPBC );  
    //builder.printBonds();
    //printf("!!!!! builder.toMMFFsp3() DONE \n");
    DEBUG
    {   //printf(" ============ check MMFFsp3_loc START\n " );
        //printf("### ffl.apos:\n");  printVecs( ffl.natoms, ffl.apos  );
        //printf("### ffl.pipos:\n"); printVecs( ffl.nnode , ffl.pipos );
        idebug=1;
        ffl.eval();
        idebug=0;
        //printf("### ffl.fneigh  :\n"); printVecs( ffl.nnode*4, ffl.fneigh   );
        //printf("### ffl.fneighpi:\n"); printVecs( ffl.nnode*4, ffl.fneighpi );
        //printf("### ffl.fapos:\n");   printVecs( ffl.natoms, ffl.fapos  );
        //printf("### ffl.fpipos:\n");  printVecs( ffl.nnode,  ffl.fpipos );
        if( ckeckNaN_d( ffl.natoms, 3, (double*)ffl.fapos,  "ffl.apos"  ) || ckeckNaN_d( ffl.nnode, 3, (double*)ffl.fpipos,  "ffl.fpipos"  ) ) { printf("ERROR: NaNs produced in MMFFsp3_loc.eval() => exit() \n"); exit(0); };
        //printf(" ============ check MMFFsp3_loc DONE\n " );
    }
    DEBUG
}

void makeFFs(){
    makeMMFF();
    initNBmol( ffl.natoms, ffl.apos, ffl.fapos, ffl.atypes ); 
    ffl.bSubtractAngleNonBond=true;
    if(ffl.bSubtractAngleNonBond){
        ffl.REQs = nbmol.REQs;
        //ff.REQs=nbmol.REQs;
    }
    bool bChargeToEpair=true;
    //bool bChargeToEpair=false;                     
    if(bChargeToEpair){
        int etyp=-1; etyp=params.atomTypeDict["E"];
        ffl.chargeToEpairs( nbmol.REQs, nbmol.atypes, -0.2, etyp );  
    }
    //nbmol.evalPLQs(gridFF.alphaMorse);
    if(bOptimizer){ 
        setOptimizer( ffl.nDOFs, ffl.DOFs, ffl.fDOFs );
    }
    _realloc( manipulation_sel, nbmol.natoms );  
}


virtual void init_empty(){

}

virtual void init( bool bGrid=false ){
    if( params.atypes.size() == 0 ){
        initParams( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
    }
    if( nbmol.natoms>0 ){ clear(); } // re-initialization
    builder.verbosity=verbosity;
    if(verbosity>0){
        printf("\n#### MolWorld_sp3::init()\n");
        if(smile_name   )printf("smile_name  (%s)\n", smile_name );
        if(xyz_name     )printf("xyz_name    (%s)\n", xyz_name );
        //if(data_dir     )printf("data_dir    (%s)\n", data_dir );
        //if(surf_name    )printf("surf_name   (%s)\n", surf_name );
        //if(substitute_name)printf("substitute_name  (%s)\n", substitute_name );
    }
    if ( smile_name ){               
        insertSMILES( smile_name );    
        builder.addAllCapTopo();       
        builder.randomizeAtomPos(1.0); 
        bMMFF=true;
    }else if ( xyz_name ){
        buildMolecule_xyz( xyz_name );
    }
    if(bMMFF){ makeFFs(); }
    if(verbosity>0) printf( "... MolWorld_sp3::init() DONE \n");
}

virtual void clear(){
    printf("MolWorld_sp3_simple.clear() \n");
    builder.clear();
    ffl.dealloc();
    // --- nbmol
    nbmol.neighs=0;   // NOTE : if we set pointer to zero it does not try to deallocate it !!!
    nbmol.apos=0;  
    nbmol.fapos=0;  
    nbmol.atypes=0;
    nbmol.dealloc();
    // --- opt
    opt.pos = 0;
    opt.force = 0;
    opt.dealloc();
}

// ======================================
// ========= Evaluation & LOOPS
// ======================================

double eval(){
    double E=0;
    if(bMMFF){ 
        E += ffl.eval();  
    }else{ VecN::set( nbmol.natoms*3, 0.0, (double*)nbmol.fapos );  }
    if(bNonBonded){
        //if  (bPBC){ E += nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ffl.lvec, {1,1,0} ); }
        if  (bPBC){ E += nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, RdampCoul ); } 
        else      { E += nbmol.evalLJQs_ng4    ( ffl.neighs );                                                }
    }
    return E;
}

void pullAtom( int ia, float K=-2.0 ){ 
    Vec3d f = getForceSpringRay( nbmol.apos[ia], pick_hray, pick_ray0, K ); nbmol.fapos[ia].add( f );
}

virtual void MDloop( int nIter, double Ftol = 1e-6 ){
    for(int itr=0; itr<nIter; itr++){
        double E = eval();
        //ckeckNaN_d( nbmol.natoms, 3, (double*)nbmol.fapos, "nbmol.fapos" );
        if(ipicked>=0){ pullAtom( ipicked );  }; // printf( "pullAtom(%i) E=%g\n", ipicked, E ); };
        //opt.move_GD(0.001);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        double f2=opt.move_FIRE();   
        //printf( "[%i] E= %g [eV] |F|= %g [eV/A]\n", nloop, E, sqrt(f2) );
        //double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
        nloop++;
    }
}

//int run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0, double* outE, double* outF ){ 
virtual int run( int nstepMax, double dt=-1, double Fconv=1e-6, int ialg=2, double* outE=0, double* outF=0 ){ 
    double F2conv=Fconv*Fconv;
    double F2 = 1.0;
    double Etot;
    int itr=0;
    //if( (ialg!=0)&(!opt_initialized) ){ printf("ERROR ialg(%i)>0 but optimizer not initialized => call initOpt() first !"); exit(0); };
    if(dt>0){ opt.setTimeSteps(dt); }
    if(ialg>0){ opt.cleanVel( ); }
    for(itr=0; itr<nstepMax; itr++ ){
        //ff.clearForce();
        Etot = eval();
        switch(ialg){
            case  0: ffl.move_GD      (opt.dt);      break;
            case -1: opt.move_LeapFrog(opt.dt);      break;
            case  1: F2 = opt.move_MD (opt.dt,opt.damping); break;
            case  2: F2 = opt.move_FIRE();          break;
            case  3: F2 = opt.move_FIRE_smooth();   break;
        }
        opt_log.set(itr, opt.cos_vf, opt.f_len, opt.v_len, opt.dt, opt.damping );
        if(outE){ outE[itr]=Etot; }
        if(outF){ outF[itr]=F2;   }
        if(verbosity>0){ printf("[%i] Etot %g[eV] |F| %g [eV/A] \n", itr, Etot, sqrt(F2) ); };
        if(F2<F2conv){
            if(verbosity>0){ printf("Converged in %i iteration Etot %g[eV] |F| %g[eV/A] <(Fconv=%g) \n", itr, Etot, sqrt(F2), Fconv ); };
            break;
        }
        if( (trj_fname) && (itr%savePerNsteps==0) ){
            sprintf(tmpstr,"# %i E %g |F| %g", itr, Etot, sqrt(F2) );
            saveXYZ( trj_fname, tmpstr, false, "a" );
        }
    }
    //printShortestBondLengths();
    return itr;
}

// ======================================
// ========= Cehcking & Debugging
// ======================================


bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){
    cog   = average( ffl.natoms, ffl.apos  );
    vcog  = sum    ( ffl.natoms, (Vec3d*)opt.vel  );
    fcog  = sum    ( ffl.natoms, ffl.fapos );
    tqcog = torq   ( ffl.natoms, ffl.apos, ffl.fapos, cog );
    //tqcog.add( ff.evalPiTorq() );
    return ( vcog.norm()>maxVcog ) || ( fcog.norm()>maxFcog ) || ( tqcog.norm() );
}

int toXYZ(const char* comment="#comment", bool bNodeOnly=false){
    if(xyz_file==0){ printf("ERROR no xyz file is open \n"); return -1; }
    int n=ffl.natoms; if(bNodeOnly){ n=ffl.nnode; }
    params.writeXYZ( xyz_file, n, nbmol.atypes, nbmol.apos, comment );
    return 0;
}

int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false, const char* mode="w" ){ return params.saveXYZ( fname, (bNodeOnly ? ffl.nnode : ffl.natoms) , nbmol.atypes, nbmol.apos, comment, nbmol.REQs, mode ); }

// ==============================================================
// ========= Edit molecule (add/remove atoms, change topology)
// ==============================================================

int substituteMolecule( const char* fname,  int ib, Vec3d up, int ipivot=0, bool bSwapBond=false, const Vec3i* axSwap=0 ){
    builder.printAtomConfs(false);
    builder.printBonds();
    printf( " ===================== Substitute molecule START !!! \n");
    Molecule* mol = new Molecule(); mol->init_xyz( fname, &params, true );
    //Vec3i axSwap={2,1,0};
    //Vec3i axSwap={2,0,1}; // THIS ONE
    //Vec3i axSwap={0,1,2};
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false );
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &axSwap, &debug_rot );
    //int ja = builder.substituteMolecule( mol, Vec3dZ, ib, ipivot, false, 0, &debug_rot );
    int ja = builder.substituteMolecule( mol, Vec3dZ, ib, ipivot, false, 0 );
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &(Vec3i{2,1,0}), &debug_rot );
    builder.addCappingTypesByIz(1);
    builder.tryAddConfsToAtoms( 0, -1 );    
    builder.sortConfAtomsFirst();              
    builder.tryAddBondsToConfs( );      
    builder.finishFragment();       
    //builder.printAtomConfs(false);
    //builder.printBonds();
    //builder.printBondParams();
    printf( "====================== Substitute molecule DONE !!! \n");
    return ja;
}

// ======================================
// ========= Manipulation & Scan
// ======================================

void shift_atoms ( int n, int* selection, Vec3d d                                  ){ move  ( n, selection, nbmol.apos, d           ); };
void rotate_atoms( int n, int* selection, Vec3d p0, Vec3d ax, double phi          ){ rotate( n, selection, nbmol.apos, p0, ax, phi ); };
void shift_atoms ( int n, int* selection, int ia0, int ia1, double l              ){ Vec3d d=(nbmol.apos[ia1]-nbmol.apos[ia0]).normalized()*l; move( n, selection, nbmol.apos, d); };
void rotate_atoms( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ rotate( n, selection, nbmol.apos, nbmol.apos[ia0], (nbmol.apos[iax1]-nbmol.apos[iax0]).normalized(), phi ); };

/*
int splitAtBond( int ib, int* selection ){
    bool bGlob=(selection==0); 
    if(bGlob){ selection=manipulation_sel; }
    int n = MM::splitByBond( ib, ffl.nbonds, ffl.bond2atom, ffl.apos, selection, manipulation_ax, manipulation_p0 );
    if(bGlob){ manipulation_nsel=n; }
    return n;
}
*/

void selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){
    Vec3d Tp0,Tp1,Tp;
    //Mat3d rot = (Mat3d)cam.rot;
    rot.dot_to(p0,Tp0);
    rot.dot_to(p1,Tp1);
    _order(Tp0.x,Tp1.x);
    _order(Tp0.y,Tp1.y);
    Tp0.z=-1e+300;
    Tp1.z=+1e+300;
    selection.clear();
    for(int i=0; i<nbmol.natoms; i++ ){
        rot.dot_to(nbmol.apos[i],Tp);
        if( Tp.isBetween(Tp0,Tp1) ){
            selection.push_back( i );
        }
    }
}

void scanTranslation_ax( int n, int* selection, Vec3d d, int nstep, double* Es, bool bWriteTrj ){
    //if(selection==0){ selection=manipulation_sel; n=manipulation_nsel; }
    //Vec3d d=(*(Vec3d*)(vec)); 
	d.mul(1./nstep);
    if(bWriteTrj){ xyz_file=fopen( "scan_trans_trj.xyz","w" ); }
    for(int i=0; i<nstep; i++){
        if(bWriteTrj){ toXYZ(); };
        double E = eval();
        if(Es)Es[i]=E;
        move( n, selection, nbmol.apos, d);
    }
    if(bWriteTrj){ fclose(xyz_file); }
}
void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj ){ Vec3d d=(nbmol.apos[ia1]-nbmol.apos[ia0]).normalized()*l; scanTranslation_ax(n,selection, d, nstep, Es, bWriteTrj ); };


void scanRotation_ax( int n, int* selection, Vec3d p0, Vec3d ax, double phi, int nstep, double* Es, bool bWriteTrj ){
    //if(p0==0) p0=(double*)&manipulation_p0;
    //if(ax==0) ax=(double*)&manipulation_ax;
    //if(selection==0){selection=manipulation_sel; n=manipulation_nsel; }
    double dphi=phi/nstep;
    if(bWriteTrj){ xyz_file=fopen( "scan_rot_trj.xyz","w" ); }
    for(int i=0; i<nstep; i++){
        double E = eval();
        Vec3d tq = torq( n, nbmol.apos, nbmol.fapos, p0, selection );
        if(bWriteTrj){  sprintf(tmpstr,"# rotScan[%i] E=%g tq=(%g,%g,%g)", i, E, tq.x,tq.y,tq.z );  toXYZ(tmpstr); };
        if(Es)Es[i]=E;
        ffl.rotateNodes(n, selection, p0, ax, dphi );
    }
    if(bWriteTrj){ fclose(xyz_file); }
}
void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj ){ Vec3d ax=(nbmol.apos[iax1]-nbmol.apos[iax0]).normalized(); scanRotation_ax(n,selection, nbmol.apos[ia0], ax, phi, nstep, Es, bWriteTrj ); };

};

#endif
