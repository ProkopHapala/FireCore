
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

//#include "raytrace.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "NBFF.h"
#include "GridFF.h"
#include "QEq.h"
#include "molecular_utils.h"

#include "Molecule.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "SMILESparser.h"
#include "DynamicOpt.h"

class MolWorld_sp3{ public:

	// Building
	MMFFparams   params;
	MM::Builder  builder;
	SMILESparser smiles;

	// Force-Fields & Dynamics
	MMFFsp3    ff;
	NBFF       nff;
	NBsystem   surf, nbmol;
	GridFF     gridFF;
    QEq        qeq;
	DynamicOpt  opt;

	// state
	bool bConverged = false;
	double Etot=0;
	double  maxVcog = 1e-9;
	double  maxFcog = 1e-9;
	double  maxTg   = 1e-1;
	double  Kmorse = -1.0;

	// force-eval-switchefs
	bool doBonded         = false;
	bool bNonBonded       = false;
	bool bSurfAtoms       = false;
	bool bPlaneSurfForce  = false;
	bool bOptimizer  = true; 
	bool bPBC        = false;
	bool bCheckInvariants = true;
	Vec3d cog,vcog,fcog,tqcog;

	// Selecteion & Manipulasion
	std::vector<int> selection;
	Vec3d manipulation_p0=Vec3dZero; 
	Vec3d manipulation_ax=Vec3dZ;
	int*  manipulation_sel=0;
	int   manipulation_nsel=0;

	// IO temp & aux
	FILE* xyz_file=0;
	char* tmpstr;


// =================== Functions

void init_nonbond(){
    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.fapos, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );
    if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
        printf( "ERROR: nff.pairMask is not sorted => exit \n" );
        exit(0);
    };
}

void autoCharges(){
    printf("autoCharges() \n");
    qeq.realloc( ff.natoms );
    params.assignQEq ( ff.natoms, ff.atype, qeq.affins, qeq.hards );
    int iconstr = params.getAtomType("E");
    printf("constrain type %i \n", iconstr );
    qeq.constrainTypes( ff.atype, iconstr );
    qeq.relaxChargeMD( ff.apos, 1000, 1e-2, 0.1, 0.1 );
    copy( qeq.n, 1, 0, (double*)qeq.qs, 3, 2, (double*)nbmol.REQs );
}

void buildFF( bool bNonBonded_, bool bOptimizer_ ){
    bOptimizer=bOptimizer_;
    bNonBonded=bNonBonded_;
    builder.autoBonds();
    //builder.autoAngles( 10.0, 10.0 );
    builder.sortConfAtomsFirst();
    builder.makeAllConfsSP();
    builder.assignAllBondParams();
    builder.toMMFFsp3( ff );
    if(bNonBonded)init_nonbond();
    if(bOptimizer){ setOptimizer(); }
    _realloc( manipulation_sel, ff.natoms );
    //init_buffers();
}

bool loadSurf(const char* fname){
	int ret = params.loadXYZ( fname, surf.n, &surf.ps, &surf.REQs );
	if(ret<0)return false;
	nbmol.bindOrRealloc( ff.natoms, ff.apos,  ff.fapos, 0 );
	params.assignREs   ( ff.natoms, ff.atype, nbmol.REQs, true, true  );
    //surf .print();
    //autoCharges(); 
    //nbmol.print();
	bSurfAtoms=true;
	return true;
}

void loadGeom(){ // TODO : overlaps with buildFF()
    // ---- Load & Build Molecular Structure
    readMatrix( "cel.lvs", 3, 3, (double*)&builder.lvec );
    builder.insertFlexibleMolecule(  builder.loadMolType( "mm.xyz", "polymer1" ), {0,0,0}, Mat3dIdentity, -1 );
    //builder.lvec.a.x *= 2.3;
    builder.printAtomConfs();
    //builder.export_atypes(atypes);
    builder.verbosity = true;
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoAngles( 10.0, 10.0 );   builder.printAngles();
    builder.toMMFFsp3( ff, &params );
    //builder.saveMol( "builder_output.mol" );
    // ----- Non-bonded interactions setup 
    if(bNonBonded){
		nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.fapos, 0, ff.bond2atom );
    	builder.export_REQs( nff.REQs );
        if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
            printf( "ERROR: nff.pairMask is not sorted => exit \n" );
            exit(0);
        };
    }else{
        printf( "WARRNING : we ignore non-bonded interactions !!!! \n" );
    }
}

int loadmol(const char* fname_mol ){
    int imol = builder.loadMolType( fname_mol, "molecule" );
    builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
    return imol;
}

void insertSMILES(const char* s){
    smiles.builder=&builder;
    smiles.parseString( 10000, s );
}

void setOptimizer(){
    opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
    double dtopt=ff.optimalTimeStep(); 
    printf("setOptimizer(): optimnal time step = %g \n", dtopt);
    opt.initOpt( dtopt );
    opt.cleanVel();
    //opt.verbosity=2;
}

void initWithSMILES(const char* s, bool bPrint=false, bool bCap=true, bool bNonBonded_=false, bool bOptimizer_=true ){
    params.init("common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
    //params.printAtomTypeDict();
    params.printAtomTypes();
    //params.printBond();
	builder.bindParams(&params);
    insertSMILES( s );
    if(bCap)builder.addAllCapTopo();
    //builder.autoAngles( 10.0, 10.0 );
    builder.randomizeAtomPos(1.0);
    // if(bPrint){
    //     printf("=============\n");
    //     printf("%s\n", s);
    //     builder.printAtoms();
    //     builder.printBonds();
    //     builder.printAtomConfs(true);
    //     //builder.printAngles();
    // }
    builder.toMMFFsp3( ff );
    if(bPrint){
        printf("=============\n");
        printf("%s\n", s);
        ff.printBonds();
        ff.printNeighs();
    }
    if(bNonBonded)init_nonbond();
    if(bOptimizer){ setOptimizer(); }
    _realloc( manipulation_sel, ff.natoms );
	printf( "MolWorld_sp3::initWithSMILES() DONE\n" );
}

void ini_in_dir(){
    params.init( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
    builder.bindParams(&params);
    int nheavy = 0;  // ---- Load Atomic Type Parameters
    if( file_exist("cel.lvs") ){ 
        loadGeom(); 
        makeGridFF();
        // ----- Optimizer setup
        //opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.fapos, 0 );
        setOptimizer();
        //double E = ff.eval(true);
    }else{
		printf("WARNING: cel.lvs not found => molecular system not initialized in [MolWorld_sp3::ini_in_dir()] \n" );
	}
}

bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){
    cog   = average( ff.natoms, ff.apos  );
    vcog  = sum    ( ff.natoms, (Vec3d*)opt.vel  );
    fcog  = sum    ( ff.natoms, ff.fapos );
    tqcog = torq   ( ff.natoms, ff.apos, ff.fapos, cog );
    //tqcog.add( ff.evalPiTorq() );
    return ( vcog.norm()>maxVcog ) || ( fcog.norm()>maxFcog ) || ( tqcog.norm() );
}

//void open_xyzFile (const char* fname){ xyz_file=fopen( fname,"w" ); };
//void close_xyzFile(){fclose(xyz_file)};

int toXYZ(const char* comment="#comment"){
    if(xyz_file==0){ printf("ERROR no xyz file is open \n"); return -1; }
    //int* atypes = builder.molTypes[0]->atomType;
    writeXYZ( xyz_file, ff.natoms, ff.atype, ff.apos, params.atomTypeNames, comment );
    return 0;
}

double eval(){
    double E=0;
    E+=ff.eval();
    if(bNonBonded)E+= nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
    return E;
}

bool relax( int niter, double Ftol = 1e-6, bool bWriteTrj=false ){
    Etot=0.0;
    double f2tol=Ftol*Ftol;
    bConverged=false; 
    if(bWriteTrj){ xyz_file=fopen( "relax_trj.xyz","w" ); }
    for(int itr=0; itr<niter; itr++){
        Etot=eval();                                                  
        if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
        double f2 = opt.move_FIRE();
        //if(bWriteTrj){ toXYZ(); ;printf("DEBUB[%i] 4 \n", itr); };
        if(bWriteTrj){  sprintf(tmpstr,"# relax[%i] E=%g f2=%g", itr, Etot, sqrt(f2) );  toXYZ(tmpstr); };
        printf( "relax[%i] |F| %g (Ftol=%g) \n", itr, sqrt(f2), Ftol );
        if(f2<f2tol){ bConverged=true; break; }
    }
    if(bWriteTrj){ fclose(xyz_file); }
    return bConverged;
}

void MDloop( int nIter, double Ftol = 1e-6 ){
    //ff.doPiPiI  =false;
    //ff.doPiPiT  =false;
    //ff.doPiSigma=false;
    //ff.doAngles =false;
    ff.cleanAll();
    for(int itr=0; itr<nIter; itr++){
        double E=0;
		E += ff.eval();
        
		if(bNonBonded){ E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} ); }
		if(bSurfAtoms){ E+=nbmol.evalMorse(surf, false, Kmorse); };
		if( bPlaneSurfForce )for(int i=0; i<ff.natoms; i++){ ff.fapos[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) ); }
        
        if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
        // if(ipicked>=0){
        //     float K = -2.0;
        //     Vec3d f = getForceSpringRay( W->ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, K );
        //     ff.aforce[ipicked].add( f );
        // };
        //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
        //opt.move_GD(0.001);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        opt.move_FIRE();
        double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
    }
}

void makeGridFF() {
    gridFF.loadXYZ  ( "inputs/NaCl_sym.xyz", params );
    gridFF.grid.n    = (Vec3i){60,60,100};
    gridFF.grid.pos0 = (Vec3d){0.0,0.0,0.0};
    gridFF.loadCell ( "inputs/cel.lvs" );
    gridFF.grid.printCell();
    gridFF.allocateFFs();
    gridFF.tryLoad( "data/FFelec.bin", "data/FFPauli.bin", "data/FFLondon.bin" );
    gridFF.shift = (Vec3d){0.0,0.0,-8.0};
}

//inline int pickParticle( const Vec3d& ray0, const Vec3d& hRay, double R, int n, Vec3d * ps, bool* ignore=0 ){
//int pickParticle( Vec3d ray0, Vec3d hray, double R=0.5 ){ return pickParticle( ray0, hray, R, ff.natoms, ff.apos ); }

// ========= Manipulation with the molecule

void shift_atoms ( int n, int* selection, Vec3d d                          ){ move  ( n, selection, ff.apos, d           ); };
void rotate_atoms( int n, int* selection, Vec3d p0, Vec3d ax, double phi   ){ rotate( n, selection, ff.apos, p0, ax, phi ); };
void shift_atoms ( int n, int* selection, int ia0, int ia1, double l              ){ Vec3d d=(ff.apos[ia1]-ff.apos[ia0]).normalized()*l; move( n, selection, ff.apos, d); };
void rotate_atoms( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ rotate( n, selection, ff.apos, ff.apos[ia0], (ff.apos[iax1]-ff.apos[iax0]).normalized(), phi ); };

int splitAtBond( int ib, int* selection ){
    bool bGlob=(selection==0); 
    if(bGlob){ selection=manipulation_sel; }
    int n = MM::splitByBond( ib, ff.nbonds, ff.bond2atom, ff.apos, selection, manipulation_ax, manipulation_p0 );
    if(bGlob){ manipulation_nsel=n; }
    return n;
}

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
    for(int i=0; i<ff.natoms; i++ ){
        rot.dot_to(ff.apos[i],Tp);
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
        move( n, selection, ff.apos, d);
    }
    if(bWriteTrj){ fclose(xyz_file); }
}
void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, bool bWriteTrj ){ Vec3d d=(ff.apos[ia1]-ff.apos[ia0]).normalized()*l; scanTranslation_ax(n,selection, d, nstep, Es, bWriteTrj ); };


void scanRotation_ax( int n, int* selection, Vec3d p0, Vec3d ax, double phi, int nstep, double* Es, bool bWriteTrj ){
    //if(p0==0) p0=(double*)&manipulation_p0;
    //if(ax==0) ax=(double*)&manipulation_ax;
    //if(selection==0){selection=manipulation_sel; n=manipulation_nsel; }
    double dphi=phi/nstep;
    if(bWriteTrj){ xyz_file=fopen( "scan_rot_trj.xyz","w" ); }
    for(int i=0; i<nstep; i++){
        double E = eval();
        Vec3d tq = torq( n, ff.apos, ff.fapos, p0, selection );
        if(bWriteTrj){  sprintf(tmpstr,"# rotScan[%i] E=%g tq=(%g,%g,%g)", i, E, tq.x,tq.y,tq.z );  toXYZ(tmpstr); };
        if(Es)Es[i]=E;
        //rotate( n, selection, ff.apos, *(Vec3d*)p0, *(Vec3d*)ax, dphi );
        ff.rotateNodes(n, selection, p0, ax, dphi );
    }
    if(bWriteTrj){ fclose(xyz_file); }
}
void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, bool bWriteTrj ){ Vec3d ax=(ff.apos[iax1]-ff.apos[iax0]).normalized(); scanRotation_ax(n,selection, ff.apos[ia0], ax, phi, nstep, Es, bWriteTrj ); };


};

#endif
