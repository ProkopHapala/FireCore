
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
MMFFparams* params_glob;

//#include "raytrace.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "NBFF.h"
#include "GridFF.h"
#include "RigidBodyFF.h"
#include "QEq.h"
#include "molecular_utils.h"

#include "Molecule.h"
#include "MMFFBuilder.h"
#include "SMILESparser.h"
#include "DynamicOpt.h"

class MolWorld_sp3{ public:
    const char* data_dir     = "common_resources";
    const char* xyz_name     = "input";
    const char* surf_name    = "surf";
    const char* substitute_name = 0;    int isubs;
    //const char* lvs_name     ="input.lvs";
    //const char* surflvs_name ="surf.lvs";
    const char* smile_name   = 0;
    Vec3i nMulPBC  = Vec3iZero; 

    bool bCellBySurf=false;
    int   bySurf_ia0 =0;
    Vec2d bySurf_c0=Vec2dZero;
    Vec2d bySurf_lat[2];

    Mat3d new_lvec=Mat3dIdentity;

    Mat3d debug_rot; // This was used for debuging molecular orientation

	// Building
	MMFFparams   params;
	MM::Builder  builder;
	SMILESparser smiles;

	// Force-Fields & Dynamics
	MMFFsp3      ff;
	NBFF         nff;
	NBsystem     surf, nbmol;
	GridFF       gridFF;
    RigidBodyFF  rbff;
    QEq          qeq;
	DynamicOpt   opt;
    DynamicOpt   optRB;  // rigid body optimizer

    bool bOcl=false; // used only in Ocl version

    double gridStep = 0.1; 
    //double gridStep = 0.2; 
    //Vec3i nPBC{0,0,0};   // just debug
    Vec3i nPBC{1,1,0};

	// state
	bool bConverged = false;
	double Etot=0;
	double  maxVcog = 1e-9;
	double  maxFcog = 1e-9;
	double  maxTg   = 1e-1;
	double  Kmorse = -1.0;

	// force-eval-switchefs
    int  imethod=0;
	bool doBonded         = false;
	bool bNonBonded       = false;
	bool bSurfAtoms       = false;
    bool bGridFF          = false;
	bool bPlaneSurfForce  = false;
    bool bMMFF            = true;
    bool bRigid           = false;
	bool bOptimizer  = true; 
	bool bPBC        = false;
	bool bCheckInvariants = true;
	Vec3d cog,vcog,fcog,tqcog;
    int nloop=0;
    bool bChargeUpdated=false;

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

virtual void swith_method(){ bGridFF=!bGridFF; };
virtual char* info_str   ( char* str=0 ){ if(str==0)str=tmpstr; sprintf(str,"bGridFF %i \n", bGridFF ); return str; }

void init_nonbond(){
    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.fapos, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );
    if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
        printf( "ERROR: nff.pairMask is not sorted => exit \n" );
        exit(0);
    };
}

void autoCharges(){
    if(verbosity>0)printf("MolWorld_sp3::autoCharges() \n");
    qeq.realloc( ff.natoms );
    params.assignQEq ( ff.natoms, ff.atype, qeq.affins, qeq.hards );
    int iconstr = params.getAtomType("E");    //printf("constrain type %i \n", iconstr );
    qeq.constrainTypes( ff.atype, iconstr );
    qeq.relaxChargeMD ( ff.apos, 1000, 1e-2, 0.1, 0.1 );
    copy( qeq.n, 1, 0, (double*)qeq.qs, 3, 2, (double*)nbmol.REQs );
    bChargeUpdated=true;
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


static void autoNPBC( const Mat3d& cell, Vec3i& nPBC, double Lmin=30.0 ){
    if(nPBC.x!=0){ nPBC.x=(int)Lmin/cell.a.norm(); }
    if(nPBC.y!=0){ nPBC.y=(int)Lmin/cell.b.norm(); }
    if(nPBC.z!=0){ nPBC.z=(int)Lmin/cell.c.norm(); }
    printf("DEBUG autoNPBC()->(%i,%i,%i) \n", nPBC.x, nPBC.y, nPBC.z );
}


virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true ){
    if(verbosity>0)printf("MolWorld_sp3::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    sprintf(tmpstr, "%s.lvs", name );
    if( file_exist(tmpstr) ){  gridFF.grid.loadCell( tmpstr, gridStep );  gridFF.bCellSet=true; }
    if( !gridFF.bCellSet ){
        bGridFF=false; 
        printf( "WARRNING!!! GridFF not initialized because %s not found\n", tmpstr );
        return;
    }
    if(bGrid){
        gridFF.grid.center_cell( cel0 );
        bGridFF=true;
        gridFF.bindSystem(surf.n, surf.atypes, surf.ps, surf.REQs );
        if( isnan(z0) ){  z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
        gridFF.grid.pos0.z=z0;
        //gridFF.grid.pos0.z=-5;  // DEBUG
        if(verbosity>1)gridFF.grid.printCell();
        gridFF.allocateFFs();
        //gridFF.tryLoad( "FFelec.bin", "FFPauli.bin", "FFLondon.bin", false, {1,1,0}, bSaveDebugXSFs );
        if(bAutoNPBC){  autoNPBC( gridFF.grid.cell, nPBC, 20.0 ); }
        //nPBC = (Vec3i){0,0,0};
        //nPBC = (Vec3i){1,1,0};
        //nPBC = (Vec3i){10,10,0};
        bSaveDebugXSFs=true;
        gridFF.tryLoad( "FFelec.bin", "FFPauli.bin", "FFLondon.bin", false, nPBC, bSaveDebugXSFs, true );
        bGridFF   =true; 
        //bSurfAtoms=false;
    }
}

void initNBmol(){
    if(verbosity>0)printf( "MolWorld_sp3::initNBmol() ff.natoms %i \n", ff.natoms  );
	nbmol  .bindOrRealloc( ff.natoms, ff.apos,  ff.fapos, 0 );
    nbmol.atypes = ff.atype;              
	builder.export_REQs  ( nbmol.REQs   );   
    for(int i=builder.atoms.size(); i<ff.natoms; i++){ nbmol.REQs[i].z=0; }  // Make sure that atoms not present in Builder has well-defined chanrge                              
    params .assignREs    ( ff.natoms, ff.atype, nbmol.REQs, true, false  ); 
    nbmol  .makePLQs     ( gridFF.alpha );    
    if(verbosity>1)nbmol.print();                              
}

void loadNBmol( const char* name){
    if(verbosity>0)printf( "MolWorld_sp3::loadNBmol() \n" );
	sprintf(tmpstr, "%s.xyz", name );
    params.loadXYZ( tmpstr, nbmol.n, &nbmol.ps, &nbmol.REQs, &nbmol.atypes );
    _realloc(nbmol.fs,nbmol.n);
    nbmol  .makePLQs     ( gridFF.alpha );    
    if(verbosity>1)nbmol.print();                              
}


bool loadSurf(const char* name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0} ){
    sprintf(tmpstr, "%s.xyz", name );
	int ret = params.loadXYZ( tmpstr, surf.n, &surf.ps, &surf.REQs, &surf.atypes, 0, &gridFF.grid.cell );
    if(ret>0){ 
        gridFF.grid.updateCell(gridStep); gridFF.bCellSet=true;
        //gridFF.grid.printCell(); 
    };
    //updateCell(step);
    if(verbosity>0)printf("MolWorld_sp3::loadSurf(%s) 1 natoms %i apos %li atyps %li \n", name, surf.n, (long)surf.ps, (long)surf.atypes  );
    //surf.print();
	if(ret<0)return false; 
	bSurfAtoms=true;
    initGridFF( name,bGrid,bSaveDebugXSFs,z0,cel0 );
	return true;
}

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
    int ja = builder.substituteMolecule( mol, Vec3dZ, ib, ipivot, false, 0, &debug_rot );
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &(Vec3i{2,1,0}), &debug_rot );
    builder.tryAddConfsToAtoms( 0, -1, 1 );    
    builder.sortConfAtomsFirst();              
    builder.tryAddBondsToConfs( );             
    //builder.printAtomConfs(false);
    //builder.printBonds();
    //builder.printBondParams();
    printf( "====================== Substitute molecule DONE !!! \n");
    return ja;
}

int loadGeom( const char* name ){ // TODO : overlaps with buildFF()
    if(verbosity>0)printf("MolWorld_sp3::loadGeom(%s)\n",  name );
    // ------ Load geometry
    sprintf(tmpstr, "%s.xyz", name ); printf("tmpstr=`%s`\n", tmpstr);
    int imol  = builder.loadMolType( tmpstr, name );
    int iret  = builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
    int ifrag = builder.frags.size()-1;
    if(iret<0){ printf("!!! exit(0) in MolWorld_sp3::loadGeom(%s)\n", name); exit(0); }
    builder.tryAddConfsToAtoms( 0, -1, 1 );
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
    if(verbosity>0)printf("MolWorld_sp3::setOptimizer(): optimnal time step = %g \n", dtopt);
    opt.initOpt( dtopt );
    opt.cleanVel();
    //opt.verbosity=2;
}

void initRigid(){
    int nrb = builder.frags.size();
    //printf("# --- initRigid() nrb=%i \n", nrb);
    int n0rot=nrb*3;
    optRB.bindOrAlloc( n0rot + nrb*4, 0, 0, 0, 0);
    rbff.realloc( nrb, (Vec3d*)optRB.pos, (Quat4d*)(optRB.pos+n0rot), (Vec3d*)(optRB.force+n0rot), (Vec3d*)(optRB.vel+n0rot), 0, 0 );
    int natom=0;
    //printf("# --- initRigid() rbff.n=%i \n", rbff.n );
    for(int i=0; i<nrb; i++){
        const MM::Fragment& frag = builder.frags[i]; // problem - some atoms are not in builder - e.g. Epair
        int i0 = frag.atomRange.x;
        int ni =  frag.atomRange.y - i0;
        //printf("# initRigid[%i] i0 %i ni %i \n", i, i0, ni );
        nbmol.ps + i0;
        rbff.mols[i].bindOrRealloc(ni, nbmol.ps+i0, nbmol.fs+i0, nbmol.REQs+i0 );
        natom+=ni;
    }
    rbff.makePos0s();
    //printf("# --- initRigid() END \n");
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
    builder.toMMFFsp3( ff );
    if(bPrint){   
        printf("=============\n"); printf("%s\n", s);
        ff.printBonds();
        ff.printNeighs();
    }
    if(bNonBonded)init_nonbond();
    if(bOptimizer){ setOptimizer(); }
    _realloc( manipulation_sel, ff.natoms );
	printf( "... MolWorld_sp3::initWithSMILES() DONE\n" );
}

void ini_in_dir(){
    params.init( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
    builder.bindParams(&params);
    int nheavy = 0;  // ---- Load Atomic Type Parameters
    if( file_exist("cel.lvs") ){ 
        loadGeom( "mm" ); 
        if(bGridFF)makeGridFF();
        // ----- Optimizer setup
        //opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.fapos, 0 );
        setOptimizer();
        //double E = ff.eval(true);
    }else{
		printf("WARNING: cel.lvs not found => molecular system not initialized in [MolWorld_sp3::ini_in_dir()] \n" );
	}
}

void PBC_multiply( Vec3i& nMulPBC_, int ifrag ){
    //printf("surface  lattice:\n"); gridFF .grid.cell.print();
    //printf("molecule lattice:\n"); builder.lvec.print();
    builder.multFragPBC( ifrag, nMulPBC_, builder.lvec );
    //printf("molecule lattice:\n"); builder.lvec.print();
    //builder.printAtoms();
    //new_lvec.ax=builder.lvec.a.norm(); new_lvec.by=builder.lvec.b.norm(); new_lvec.cz=builder.lvec.c.norm();
    builder.correctPBCbonds( ifrag, builder.frags.size() ); // correct bonds for newly added fragments
    builder.checkBondsInNeighs(true); // DEBUG
    builder.sortConfAtomsFirst(); 
    //printf("molecule lattice:\n"); builder.lvec.print();
    //builder.printAtomConfs();
    //builder.printBonds();    
}

void changeCellBySurf( Vec2d a, Vec2d b, int ia0=-1, Vec2d c0=Vec2dZero ){
    //printf( "changeCellBySurf() a(%g,%g) b(%g,%g) \n", a.x,a.y, b.x,b.y  );
    double la0=builder.lvec.a.norm();
    double lb0=builder.lvec.b.norm();
    Mat3d lvs;
    lvs.a=gridFF.grid.cell.a*a.a + gridFF.grid.cell.b*a.b;
    lvs.b=gridFF.grid.cell.a*b.a + gridFF.grid.cell.b*b.b; 
    lvs.c=builder.lvec.c;
    builder.changeCell( lvs );
    //Vec3d pmin,pmax; builder.bbox(pmin,pmax); printf( "BBOX pmin(%g,%g,%g) pmax(%g,%g,%g)\n", pmin.x,pmin.y,pmin.z,  pmax.x,pmax.y,pmax.z ); builder.move_atoms(pmin*-1);
    //builder.move_atoms( builder.atoms[0].pos*-1.);
    if(ia0>=0){
        Vec3d shift =  builder.atoms[ia0].pos*-1 + gridFF .grid.cell.a*c0.a + gridFF .grid.cell.b*c0.b;
        builder.move_atoms( shift );
    }
    printf( "changeCellBySurf() DONE, |a,b|=%g,%g (old |a,b|=%g,%g) ", builder.lvec.a.norm(), builder.lvec.b.norm(), la0, lb0 );
    //builder.printAtomConfs();
    //builder.printBonds();    
}

virtual void init( bool bGrid ){
    params.init("common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
	builder.bindParams(&params);
    //params.printAtomTypeDict();
    //params.printAtomTypes();
    //params.printBond();

    params_glob = &params;

    builder.verbosity=verbosity;
    if(verbosity>0){
        printf("\n#### MolWorld_sp3::init()\n");
        if(smile_name   )printf("smile_name  (%s)\n", smile_name );
        if(data_dir     )printf("data_dir    (%s)\n", data_dir );
        if(xyz_name     )printf("xyz_name    (%s)\n", xyz_name );
        if(surf_name    )printf("surf_name   (%s)\n", surf_name );
        if(substitute_name)printf("substitute_name  (%s)\n", substitute_name );
        //if(lvs_name     )printf("lvs_name    (%s)\n", lvs_name );
        //if(surflvs_name )printf("surflvs_name(%s)\n", surflvs_name );
        printf("bMMFF %i bRigid %i \n", bMMFF, bRigid );
    }
    if(surf_name )loadSurf( surf_name, bGrid, idebug>0 );
    if ( smile_name ){               
        insertSMILES( smile_name );    
        builder.addAllCapTopo();       
        builder.randomizeAtomPos(1.0); 
        bMMFF=true;
    }else if ( xyz_name ){
        if( bMMFF ){ 
            int ifrag = loadGeom( xyz_name );

            //printf( "DEBUG substituteMolecule(%i,%s) \n", isubs, substitute_name );
            if(substitute_name) substituteMolecule( substitute_name, isubs, Vec3dZ );
            //int substituteMolecule( const char fname,  int ib, Vec3d up, int ipivot=0, bool bSwapBond=false, const Vec3i* axSwap=0 ){

            builder.printAtomConfs();
            if( builder.checkNeighsRepeat( true ) ){ printf( "ERROR: some atoms has repating neighbors => exit() \n"); exit(0); };

            //builder.printAtomConfs(true);
            builder.autoAllConfEPi( );          //builder.printAtomConfs(true);
            //builder.makeAllConfsSP(true);     if(verbosity>1)builder.printAtomConfs(true);
            //builder.printAtomConfs(true);
            builder.assignAllBondParams();    //if(verbosity>1)
            //builder.printBonds    ();
            //builder.autoAngles( 10.0, 10.0 );   builder.printAngles();
            //builder.toMMFFsp3( ff, &params );
            //builder.saveMol( "builder_output.mol" );
            //if(bNonBonded){ init_nonbond(); }else{ printf( "WARRNING : we ignore non-bonded interactions !!!! \n" ); }
            builder.finishFragment(ifrag);

            if( nMulPBC    .totprod()>1 ){ PBC_multiply    ( nMulPBC, ifrag ); };
            if( bCellBySurf             ){ changeCellBySurf( bySurf_lat[0], bySurf_lat[1], bySurf_ia0, bySurf_c0 ); };
        }else{
            loadNBmol( xyz_name ); 
            if(bRigid)initRigid();
        }
    }
    if(bMMFF){       
        //builder.printAtoms();          
        //if( builder.checkBondsOrdered( false, true ) ) { printf("ERROR Bonds are not ordered => exit"); exit(0); };
        if( builder.checkBondsInNeighs(true) ) { 
            printf("ERROR some bonds are not in atom neighbors => exit"); 
            exit(0); 
        };
        builder.checkBondsOrdered( true, false );
        builder.toMMFFsp3( ff, &params );
        //builder.printBonds();
        //printf("!!!!! builder.toMMFFsp3() DONE \n");
        if(verbosity>0){
            ff.printSizes();
            ff.printAtoms();
            ff.printNeighs();
            ff.printBonds();
            ff.printAtomPis();
        }
        if( ff.checkBonds( 1.5, true ) ){ printf("ERROR Bonds are corupted => exit"); exit(0); };
        { // check MMFF
            ff.eval();
            if(ff.checkNaNs()){ printf("ERROR: NaNs produced in MMFFsp3.eval() => exit() \n"); exit(0); };
        }
        initNBmol();
        int etyp=-1; etyp=params.atomTypeDict["E"];
        ff.chargeToEpairs( nbmol.REQs, -0.2, etyp );  
        nbmol.evalPLQs(gridFF.alpha);

        if(bOptimizer){ setOptimizer(); }                         
        _realloc( manipulation_sel, ff.natoms );  
    }
    if(verbosity>0) printf( "... MolWorld_sp3::init() DONE \n");
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

int toXYZ(const char* comment="#comment", bool bNodeOnly=false){
    if(xyz_file==0){ printf("ERROR no xyz file is open \n"); return -1; }
    int n=ff.natoms; if(bNodeOnly){ n=ff.nnode; }
    params.writeXYZ( xyz_file, n, ff.atype, ff.apos, comment );
    return 0;
}

//int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false){ return params.saveXYZ( fname, (bNodeOnly ? ff.nnode : ff.natoms) , ff.atype, ff.apos, comment, nbmol.REQs ); }
int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false){ return params.saveXYZ( fname, (bNodeOnly ? ff.nnode : ff.natoms) , nbmol.atypes, nbmol.ps, comment, nbmol.REQs ); }

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

double eval(){
    double E=0;
    if(bMMFF){ E += ff.eval();  }else{ VecN::set( nbmol.n*3, 0.0, (double*)nbmol.fs );  }
    if(bSurfAtoms){ 
        if   (bGridFF){ E+= gridFF.eval(nbmol.n, nbmol.ps, nbmol.PLQs, nbmol.fs ); }
        //else        { E+= nbmol .evalMorse   (surf, false,                   gridFF.alpha, gridFF.Rdamp );  }
        else          { E+= nbmol .evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp );  }
    }
    //printf( "eval() bSurfAtoms %i bGridFF %i \n", bSurfAtoms, bGridFF );
    //for(int i=0; i<nbmol.n; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fs[i].x,nbmol.fs[i].y,nbmol.fs[i].z ); }
    return E;
}

void pullAtom( int ia, float K=-2.0 ){ 
    Vec3d f = getForceSpringRay( ff.apos[ia], pick_hray, pick_ray0, K );
    ff.fapos[ia].add( f );
}

virtual void MDloop( int nIter, double Ftol = 1e-6 ){
    //ff.doPiPiI  =false;
    //ff.doPiPiT  =false;
    //ff.doPiSigma=false;
    //ff.doAngles =false;
    ff.cleanAll();
    for(int itr=0; itr<nIter; itr++){
        double E = eval();
        ckeckNaN_d( nbmol.n, 3, (double*)nbmol.fs, "nbmol.fs" );
		//if( bPlaneSurfForce )for(int i=0; i<ff.natoms; i++){ ff.fapos[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) ); }
        //printf( "apos(%g,%g,%g) f(%g,%g,%g)\n", ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.fapos[0].x,ff.fapos[0].y,ff.fapos[0].z );
        //if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
        if(ipicked>=0){ pullAtom( ipicked ); };
        //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
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
    bChargeUpdated=false;
}

void makeGridFF( bool bSaveDebugXSFs=false, Vec3i nPBC={1,1,0} ) {
    gridFF.bindSystem(surf.n, 0, surf.ps, surf.REQs );
    //if(verbosity>1)
    gridFF.grid.printCell();
    gridFF.allocateFFs();
    //double x0= ( gridFF.grid.cell.a.x + gridFF.grid.cell.b.x )*-0.5;
    //double y0= ( gridFF.grid.cell.a.y + gridFF.grid.cell.b.y )*-0.5;
    //gridFF.grid.pos0 = (Vec3d){ x0,y0,-8.0};
    //gridFF.shift   = (Vec3d){0.0,0.0,-8.0};
    //gridFF.tryLoad( "data/FFelec.bin", "data/FFPauli.bin", "data/FFLondon.bin", true, {0,0,0} );
    //gridFF.tryLoad( "FFelec.bin", "FFPauli.bin", "FFLondon.bin", false, {1,1,1} );
    gridFF.tryLoad( "FFelec.bin", "FFPauli.bin", "FFLondon.bin", false, nPBC, bSaveDebugXSFs );
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
