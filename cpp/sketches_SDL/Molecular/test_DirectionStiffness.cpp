
#include <globals.h>
//int verbosity = 0;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include <SDL2/SDL.h>

#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFmini.h"
#include "MMFFsp3.h"
#include "NBFF_old.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

#include "DirectionStiffness.h"
#include "Graph.h"
#include "MolecularGraph.h"



#include "Draw3D_Molecular.h"  // it needs to know MMFFparams

//#include "NBSRFF.h"
#include "IO_utils.h"

#include "AppSDL2OGL_3D.h"


Vec3d*    apos0=0;
DirStiff* p_Dstiff=0;
MMFFmini* p_mmff=0;

void CG_DotFunc( int n, const double * x, double * Ax ){
    VecN::add( n,  (double*)apos0,  x,  (double*)p_mmff->apos );
    p_mmff->eval_bonds();
    VecN::set( n, (double*)p_mmff->aforce, Ax );
};

void relax1(int nmax, double Fmax, double dt){
    //double Fmax2=Fmax*Fmax;
    for(int i=0; i<nmax; i++){
        p_mmff->eval_bonds();
        double f = p_mmff->getFmax();
        if(f<Fmax) break;
        for(int i=0; i<p_mmff->natoms; i++){ p_mmff->apos[i].add_mul( p_mmff->aforce[i], dt ); };
    }
}

double defeormer_evalForce(int n, const double * x, double * Ax){
    double E=0;
    p_mmff->cleanAtomForce();
    E+=p_mmff->eval_bonds();
    E+=p_mmff->eval_angles();
    return E;
}

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppDirectionStiffness : public AppSDL2OGL_3D { public:
	Molecule    mol;
	MMFFparams  params;
    MMFFmini    ff;
    NBFF_old    nff;
    MM::Builder builder;
    DynamicOpt  opt;

    Graph* graph2=0;
    MolecularGraph graph;
    Deformer deformer;
    int ndeform=0;
    int deform_Nsteps=50;
    int deform_K     =10.0;



    int* atypes = 0;

    bool bNonBonded = true;

    std::vector<Molecule*> molecules;

    std::vector<int> selection;
    bool bDragging = false;
    Vec3d  ray0_start;
    Vec3d  rotation_center = Vec3dZero;
    Vec3d  rotation_axis   = Vec3dZ;
    double rotation_step   = 15.0 * (M_PI/180.0);

    Vec3d lvec_a0;
    int icell = 0;
    int frameCountPrev=0;

    //std::unordered_map<std::string,int> atomTypeDict;

    //Mat3d lvec;

    bool bConverged = false;
    bool bRunRelax  = false;

    int     fontTex;
    GLMesh ogl_sph = Draw3D::makeSphereOgl( 5, 1.0 );
    int     ogl_mol=0;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1, iangPicked = -1;
    int perFrame =  50;


    double drndv =  10.0;
    double drndp =  0.5;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppDirectionStiffness( int& id, int WIDTH_, int HEIGHT_ );

    void MDloop(  );

	int makeMoleculeInline();
	int makeMoleculeInlineBuilder( bool bPBC );
	int loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int loadMoleculeXYZ( const char* fname, bool bAutoH );
	int loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );

	void drawSystem( );

	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

};

TestAppDirectionStiffness::TestAppDirectionStiffness( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ------ using molecule from mol-file does not seem to work now - There is some problem with AtomConfs
    // >> This algorithm assumes all atoms with conf precede atoms without confs in the array
    // >>   ERROR in builder.sortBonds() => exit

    int nheavy = 0;
    //nheavy = loadMoleculeXYZ( "common_resources/polymer.xyz", false );
    //nheavy = loadMoleculeXYZ    ( "common_resources/polymer-2.xyz", "common_resources/polymer-2.lvs", false );
    //nheavy = loadMoleculeXYZ    ( "common_resources/polymer-2-monomer.xyz", "common_resources/polymer-2.lvs", false );


    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);

    readMatrix( "common_resources/polymer-2.lvs", 3, 3, (double*)&builder.lvec );
    molecules.push_back( new Molecule() ); molecules[0]->atomTypeDict = builder.atomTypeDict; molecules[0]->load_xyz("common_resources/polymer-2.xyz", true);
    //molecules.push_back( new Molecule() ); molecules[1]->atomTypeDict = builder.atomTypeDict; molecules[1]->load_xyz("common_resources/polymer-2-monomer.xyz", true);
    builder.insertFlexibleMolecule(  molecules[0], {0,0,0}       , Mat3dIdentity, -1 );
    //builder.insertFlexibleMolecule(  molecules[1], builder.lvec.a*1.2, Mat3dIdentity, -1 );

    builder.lvec.a.x *= 2.3;

    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    builder.export_atypes(atypes);
    //builder.verbosity = 5;
    //builder.autoBonds ();             builder.printBonds ();
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     builder.printAngles();
    builder.toMMFFmini( ff, &params );
    builder.saveMol( "data/polymer.mol" );

    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.aforce, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );

    //bNonBonded = false;
    if(bNonBonded){
        if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
            printf( "ERROR: nff.pairMask is not sorted => exit \n" );
            exit(0);
        };
    }else{
        printf( "WARRNING : we ignore non-bonded interactions !!!! \n" );
    }

    graph.bindOrRealloc( ff.natoms, ff.nbonds, ff.bond2atom );
    graph.makeNeighbors();
    graph.printNeighs();
    //exit(0);

    graph2 = new Graph( builder.confs.size());
    graph2->fromBonds( ff.nbonds, ff.bond2atom );
    graph2->bridge();
    //exit(0);

    deformer.bind( ff.natoms, ff.apos, ff.aforce );
    deformer.initPicks( 2 );
    deformer.picks[0]=17; // select bond #17
    deformer.picks[1]=42; // select bond #17
    deformer.genKs( 1.0 );
    //deformer.genPicks();
    //deformer.genPulls();
    deformer.graph=&graph;
    deformer.evalForce = defeormer_evalForce;
    p_mmff = &ff;

    opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.aforce, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );
    opt.initOpt( 0.1, 0.1 );

    // ======== Test before we run
    nff.printAtomParams();
    ff.printAtomPos();
    ff.printBondParams();
    ff.printAngleParams();
    ff.printTorsionParams();
    double E = ff.eval(true);
    printf( "iter0 E = %g \n", E );
    printf("TestAppDirectionStiffness.init() DONE \n");
    //exit(0);


    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //opengl1renderer.lightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//opengl1renderer.lightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    //selection.insert( selection.end(), {12, 16, 14, 6, 2, 3,   20,18,31,25,26} );
    //selection.insert( selection.end(), {13,29,30} );
    //splitGraphs( ff.nbonds, ff.bond2atom, 12, 22 );
    std::unordered_set<int> innodes;
    //innodes.insert(12);
    innodes.insert(11);
    MM::splitGraphs( ff.nbonds, ff.bond2atom, 22, innodes );
    for( int i : innodes ){ selection.push_back(i); }

}


void TestAppDirectionStiffness::MDloop(  ){

    double Ftol = 1e-6;
    perFrame = 25;
    builder.lvec.a    = lvec_a0 + Vec3d{-1.0,0.0,0.0};

    for(int itr=0; itr<perFrame; itr++){

        deform_Nsteps=50;
        deform_K     =10.0;

        //ndeform=1;
        if(ndeform>0){
            //deformer.deform_F();
            //deformer.deform_Rot();
            deformer.deform_BondRot( deform_K );
            ndeform--;
        }else{

            double E=0;
            ff.cleanAtomForce();
            E += ff.eval(false);
            if(bNonBonded){
                //E += nff.evalLJQ_sortedMask();   // This is fast but does not work in PBC
                E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
            }
            if(ipicked>=0){
                Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rotMat().c, ray0, -1.0 );
                //printf( "f (%g,%g,%g)\n", f.x, f.y, f.z );
                ff.aforce[ipicked].add( f );
            };

            float K = -0.01;
            for(int i=0; i<ff.natoms; i++){
                ff.aforce[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) );

            }

            //float d = 0.05;
            //for(int i=0; i<ff.natoms; i++){ ff.aforce[i].add({randf(-d,d),randf(-d,d),randf(-d,d)});  };

            //double f2; opt.move_MD( 0.1, 0.005 );
            double f2 = opt.move_FIRE( );

            if(f2<sq(Ftol)){
                bConverged=true;
            }
        }
    }
}


void TestAppDirectionStiffness::draw(){
    //opengl1renderer.clearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    opengl1renderer.clearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	opengl1renderer.enable(GL_DEPTH_TEST);

    if(frameCount==1){

        cam.qrot.pitch( M_PI );
        ff.printAtomPos();
        ff.printBondParams();
        ff.printAngleParams();
        ff.printTorsionParams();

        lvec_a0 = builder.lvec.a;
        printf( "lvec_a0  (%g %g,%g) \n", lvec_a0.x, lvec_a0.y, lvec_a0.z );
    }

	if( ogl_mol ){
        opengl1renderer.callList( ogl_mol );
        return;
        //exit(0);
    }

    ray0 = (Vec3d)(cam.rotMat().a*mouse_begin_x + cam.rotMat().b*mouse_begin_y);
    Draw3D::drawPointCross( renderer, ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( renderer, ff.apos[ipicked], ray0, {0, 0, 0});



    bool makeScreenshot = false;
    if(bRunRelax){
        MDloop();
    }


    drawSystem();
    for(int i=0; i<deformer.npick; i++){  int ia=deformer.picks[i]; Draw3D::drawVecInPos( renderer, deformer.aforce[ia]*15.0,deformer.apos[ia], {0, 1, 0}); }

    opengl1renderer.disable(GL_DEPTH_TEST);
    for(const Vec2i& b : graph2->found ){
        Draw3D::drawLine( renderer, ff.apos[b.a], ff.apos[b.b], {0, 1, 1} );
    }

};

void TestAppDirectionStiffness::drawSystem( ){
    opengl1renderer.color3f(1.0f,0.0f,0.0f); Draw3D::vecsInPos( ff.natoms, ff.aforce,  ff.apos, 10.0 );
    //opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::drawLines ( ff.nbonds, (int*)ff.bond2atom, ff.apos );
    //opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC ( ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec ); // DEBUG
    Draw3D::bondsPBC ( renderer, ff.nbonds, ff.bond2atom, ff.apos, ff.pbcShifts, {0, 0, 0} ); // DEBUG
    //opengl1renderer.color3f(0.5f,0.0f,0.0f); Draw3D::atomLabels( ff.natoms, ff.apos, fontTex                     );                     //DEBUG
    //opengl1renderer.color3f(0.0f,0.0f,1.0f); Draw3D::bondLabels( ff.nbonds, ff.bond2atom, ff.apos, fontTex, 0.02 );                     //DEBUG
    //opengl1renderer.color3f(0.0f,0.0f,1.0f); Draw3D::atomPropertyLabel( ff.natoms, (double*)nff.REQs, ff.apos, 3,2, fontTex, 0.02, "%4.2f\0" );

    //opengl1renderer.color3f(1.0f,0.0f,0.0f); Draw3D::vecsInPoss( ff.natoms, ff.aforce, ff.apos, 300.0              );
    //Draw3D::atomsREQ  ( ff.natoms, ff.apos,   nff.REQs, ogl_sph, 1.0, 0.25, 1.0 );
    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 1.0, 1.0 );       //DEBUG
    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.5, 1.0 );       //DEBUG
    Draw3D::atoms( renderer, ff.natoms, ff.apos, atypes, params, &ogl_sph, 1.0, 0.25, 1.0 );       //DEBUG
}

int TestAppDirectionStiffness::loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH ){
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);
    int nheavy = builder.load_xyz( fname, bAutoH, true );
    readMatrix( fnameLvs, 3, 3, (double*)&builder.lvec );

    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    //exit(0);
    builder.export_atypes(atypes);

    //builder.verbosity = 5;
    //builder.autoBonds ();             builder.printBonds ();
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     builder.printAngles();

    //bNonBonded = false;
    //exit(0);
    builder.toMMFFmini( ff, &params );

    builder.saveMol( "data/polymer.mol" );

    return nheavy;
}


int TestAppDirectionStiffness::loadMoleculeMol( const char* fname, bool bAutoH, bool bLoadTypes ){

    /// should look in   test_SoftMolecularDynamics.cpp

    if(bLoadTypes){
        printf( "bLoadTypes==True : load atom and bond types from file \n" );
        params.loadAtomTypes( "common_resources/AtomTypes.dat" );
        params.loadBondTypes( "common_resources/BondTypes.dat");
        //builder.params = &params;
    }else{
        printf( "bLoadTypes==False : make default Atom names dictionary \n" );
        params.initDefaultAtomTypeDict( );
    }
    //mol.atomTypeDict  = &params.atomTypeDict;
    //mol.atomTypeNames = &params.atomTypeNames;
    mol.bindParams(&params);
    mol.loadMol( fname );
    //mol.loadMol_const( "common_resources/propylacid.mol");
    //mol.loadMol_const( "/home/prokop/Dropbox/TEMP/ERC2021/Molecules/chain--frag-4---N-----.mol" );
    //exit(0);
    //int iH = 1;
    int iH = params.atomTypeDict["H"];
    int nh     = mol.countAtomType(iH); printf( "nh %i\n", nh );
    int nheavy = mol.natoms - nh;
    if(bAutoH){
        printf( "bAutoH==True : MMFFBuilder will re-calculate hydrogens, pi-orbitals and free electron pairs \n" );
        builder.insertFlexibleMolecule_ignorH( &mol, Vec3dZero, Mat3dIdentity, iH );
        //builder.setConfs( 0, 0, 0, mol.natoms-nh );
        //for(int i=0;i<(mol.natoms-nh);i++){ builder.makeSPConf(i,0,0); }
        //for(int i=0;i<mol.natoms;i++)     { builder.makeSPConf(i,0,0); }
    }else{
        printf( "bAutoH==False : Angles assigned by simple algorithm Molecule::autoAngles \n" );
        //mol.bondsOfAtoms();
        params.assignREs( mol.natoms, mol.atomType, mol.REQs );
        mol.autoAngles(true);
        Vec3d cog = mol.getCOG_av();
        mol.addToPos( cog*-1.0 );
        builder.insertMolecule(&mol, Vec3dZero, Mat3dIdentity, false );
        builder.toMMFFmini( ff, &params );
    }

    //builder.sortAtomsOfBonds();
    builder.addCappingTypesByIz(1);
    builder.tryAddConfsToAtoms(0, nh);
    builder.tryAddBondsToConfs();
    //for(int i=0; i<nh; i++){ builder.addConfToAtom(i); }
    //builder.tryAddBondsToConfs();

    //mol.printAtomInfo();
    //mol.printAtom2Bond();
    //mol.printAngleInfo();
    builder.printAtoms();
    //builder.printBonds();
    //builder.printAngles();
    //builder.printConfs();

    //bNonBonded = false;      // ToDo : WARRNING : this is just hack, because builder.sortBonds() does not seem to work, we have to switch off Non-Bonding interactions
    builder.trySortBonds();
    //builder.sortBonds();
    builder.printBonds();
    builder.printAngles();
    builder.printConfs();
    builder.toMMFFmini( ff, &params );

    //Draw3D::shapeInPoss( ogl_sph, ff.natoms, ff.apos, 0 );

    return nheavy;
}

int TestAppDirectionStiffness::makeMoleculeInlineBuilder( bool bPBC ){
    //const int natom=4,nbond=3,nang=2,ntors=1;
    //const int natom=4,nbond=3,nang=0,ntors=0;
    const int natom=4,nbond=4;
    Vec3d apos0[] = {
        {-2.0,0.0,0.0},  // 0
        {-1.0,2.0,0.0},  // 1
        {+1.0,2.0,0.0},  // 2
        {+2.0,0.0,1.0},  // 3
    };
    Vec2i bond2atom[] = {
        {0,1},  // 0
        {1,2},  // 1
        {2,3},  // 2
        {3,0},  // 3  // PBC
    };
    // ============ Build molecule
    MM::Atom brushAtom        { -1, -1, -1 , Vec3dZero, MM::Atom::defaultREQ };
    MM::Bond brushBond        { -1, {-1,-1}, 1.5,  25.0 };
    builder.capBond = MM::Bond{ -1, {-1,-1}, 1.07, 15.0 };

    builder.insertAtoms( natom, brushAtom,  apos0    );
    builder.insertBonds( nbond, brushBond, bond2atom );
    builder.setConfs   ( 0, 0 );
    builder.autoAngles ( 2.5, 1.25 );

    // instert aditional dihedral
    MM::Dihedral brushDihedral{ -1,   Vec3i{-1,-1,-1},    0,3, 0.5 };  println(brushDihedral);
    builder.insertDihedralByAtom( {0,1,2,3}, brushDihedral );
    builder.trySortBonds();

    builder.toMMFFmini( ff, &params );

    builder.lvec.a = Vec3d{  5.0,0.0,0.0 };
    builder.lvec.b = Vec3d{  0.0,5.0,0.0 };
    builder.lvec.c = Vec3d{  0.0,0.0,5.0 };
    if(bPBC){    // --- Periodic Boundary Conditions
        ff.initPBC();                 // as far as good, pbc-shifts are curenlty zero, so no change
        ff.pbcShifts[1] = builder.lvec.a*-1.; // make bond 3 from nighboring cell
        ff.printBondParams();
    }

    return natom;
}

int TestAppDirectionStiffness::makeMoleculeInline(){
    printf( "----- makeMoleculeInline \n" );
    const int natom=5+2,nbond=4+3,nang=6, ntors=0;
    Vec3d apos0[] = {
        { 0.5, 0.5, 0.5},  // 0
        {-1.0,+1.0,+1.0},  // 1
        {+1.0,-1.0,+1.0},  // 2
        {-1.0,-1.0,-1.0},  // 3
        {+1.0,+1.0,-1.0},  // 4
        {-1.0,-1.0,-2.0},  // 5
        {+1.0,+1.0,-2.0}   // 6
    };
    Vec2i bong2atom[] = {
        {0,1},  // 0
        {0,2},  // 1
        {0,3},  // 2
        {0,4},  // 3
        {5,6},  // 4
        {3,5},  // 5
        {4,6}   // 6
    };
    Vec2i ang2bond[] = {
        {0,1},
        {0,2},
        {0,3},
        {1,2},
        {1,3},
        {2,3}
    };
    double a0s[] ={
        2.0,
        2.0,
        2.0,
        2.0,
        2.0,
        2.0
    };
    Vec3i tors2bond[] = { };
    double l0    = 1.5;
    double Kbond = 10.0;
    double Kang  = 3.0;
    double Ktors = 1.0;
    int tors_n = 3;
    ff.realloc(natom,nbond,nang,ntors);
    printf( "----- Atoms \n" );
    for(int i=0; i<ff.natoms; i++){
        ff.apos[i] = apos0[i];
    }
    printf( "----- Bonds \n" );
    for(int i=0; i<ff.nbonds; i++){
        ff.bond2atom[i] = bong2atom[i];
        ff.bond_k [i] = Kbond;
        ff.bond_l0[i] = l0;
    }
    printf( "----- Angles \n" );
    for(int i=0; i<ff.nang; i++){
        ff.ang2bond[i] = ang2bond[i];
        double a0 = -a0s[i]/2.0; // NOTE: we use half-angle
        ff.ang_cs0[i] = { cos(a0), sin(a0) };
        ff.ang_k  [i] = Kang;
    }
    printf( "----- Dihedrals \n" );
    for(int i=0; i<ff.ntors; i++){
        ff.tors2bond[i] = tors2bond[i];
        ff.tors_k   [i] = Ktors;
        ff.tors_n   [i] = tors_n;
    }
    ff.angles_bond2atom  ();
    ff.torsions_bond2atom();
    return natom;
};

//void TestAppDirectionStiffness::makeAtoms(){}
//template<typename T> std::function<T(const T&,const T&         )> F2;


void TestAppDirectionStiffness::saveScreenshot( int i, const char* fname ){
    //if(makeScreenshot){
        char str[64];
        sprintf( str, fname, i );               // DEBUG
        printf( "save to %s \n", str );
        unsigned int *screenPixels = new unsigned int[WIDTH*HEIGHT*4];  //DEBUG
        opengl1renderer.flush();                                                      //DEBUG
        opengl1renderer.finish();                                                     //DEBUG
        //opengl1renderer.readPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_INT, screenPixels);
        opengl1renderer.readPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, screenPixels);   //DEBUG
        //SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff );   //DEBUG
        SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000 );   //DEBUG
        SDL_SaveBMP(bitmap, str);    //DEBUG
        SDL_FreeSurface(bitmap);
        delete[] screenPixels;
}

void TestAppDirectionStiffness::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                //case SDLK_LEFTBRACKET:  ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;
                //case SDLK_RIGHTBRACKET: ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;

                //case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                //case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                //case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                //case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                //case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                //case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                //case SDLK_w: world.apos[1].mul( 1.1 ); break;
                //case SDLK_s: world.apos[1].mul( 0.9 ); break;

                //case SDLK_KP_7: builder.lvec.a.mul(   1.01); break;
                //case SDLK_KP_4: builder.lvec.a.mul( 1/1.01); break;

                case SDLK_KP_7: builder.lvec.a.x+=xstep; break;
                case SDLK_KP_4: builder.lvec.a.x-=xstep; break;

                case SDLK_KP_8: builder.lvec.a.y+=xstep; break;
                case SDLK_KP_5: builder.lvec.a.y-=xstep; break;

                case SDLK_KP_9: builder.lvec.a.z+=xstep; break;
                case SDLK_KP_6: builder.lvec.a.z-=xstep; break;

                case SDLK_RETURN:{ ndeform=deform_Nsteps; deformer.genKs(1); }break;

                case SDLK_f:
                    //selectShorterSegment( (Vec3d)(cam.rotMat().a*mouse_begin_x + cam.rotMat().b*mouse_begin_y + cam.rotMat().c*-1000.0), (Vec3d)cam.rotMat().c );
                    //selectShorterSegment( ray0, (Vec3d)cam.rotMat().c );
                    //selection.erase();
                    //for(int i:builder.selection){ selection.insert(i); };
                    break;

                case SDLK_LEFTBRACKET:
                    rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, +rotation_step );
                    break;
                case SDLK_RIGHTBRACKET:
                    rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, -rotation_step );
                    break;

                case SDLK_SPACE: bRunRelax=!bRunRelax; break;

                case SDLK_g: iangPicked=(iangPicked+1)%ff.nang;
                    printf( "ang[%i] cs(%g,%g) k %g (%i,%i,%i)\n", iangPicked, ff.ang_cs0[iangPicked].x, ff.ang_cs0[iangPicked].y, ff.ang_k[iangPicked],
                        ff.ang2atom[iangPicked].a,ff.ang2atom[iangPicked].b,ff.ang2atom[iangPicked].c );
                    break;

            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    /*
                    ipicked = pickParticle( ray0, (Vec3d)cam.rotMat().c, 0.5, ff.natoms, ff.apos );
                    selection.clear();
                    if(ipicked>=0){ selection.push_back(ipicked); };
                    printf( "picked atom %i \n", ipicked );
                    */
                    ray0_start = ray0;
                    bDragging = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = ff.pickBond( ray0, (Vec3d)cam.rotMat().c , 0.5 );
                    //printf("ibpicked %i \n", ibpicked);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = -1;
                    //ray0_start
                    if( ray0.dist2(ray0_start)<0.1 ){
                        ipicked = pickParticle( ray0, (Vec3d)cam.rotMat().c, 0.5, ff.natoms, ff.apos );
                        selection.clear();
                        if(ipicked>=0){ selection.push_back(ipicked); };
                        printf( "picked atom %i \n", ipicked );
                    }else{
                    //    selectRect( ray0_start, ray0 );
                    }
                    bDragging=false;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void TestAppDirectionStiffness::drawHUD(){
    opengl1renderer.disable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppDirectionStiffness * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppDirectionStiffness( junk , 1000, 800 );
	//thisApp = new TestAppDirectionStiffness( junk , 800, 600 );
	//thisApp = new TestAppDirectionStiffness( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















