
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFsp3.h"
#include "NBFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams

//#include "NBSRFF.h"
#include "IO_utils.h"

#include "AppSDL2OGL_3D.h"

bool bPrint = true;


// ================= Free Functions ==============

void drawBonds( const MM::Builder& builder ){
    //drawSystem( false, true, false );
    glBegin( GL_LINES );
    for(int ib=0; ib<builder.bonds.size(); ib++ ){
        const MM::Bond& b = builder.bonds[ib]; 
        //printf( "bond[%i] (%i,%i) \n)", ib, b.atoms.a, b.atoms.b );
        //Draw3D::drawLine( builder.atoms[b.atoms.a].pos, builder.atoms[b.atoms.b].pos );
        Draw3D::vertex(builder.atoms[b.atoms.a].pos);
        Draw3D::vertex(builder.atoms[b.atoms.b].pos);
    }
    glEnd();
}

void drawNeighs( const MM::Builder& builder ){
    glBegin( GL_LINES );
    //printf( "DEBUG drawNeighs() \n" );
    for(int ia=0; ia<builder.atoms.size(); ia++ ){
        const MM::Atom& a = builder.atoms[ia];
        if( a.iconf<0 ) continue;
        //printf( "a.iconf %i \n", a.iconf );
        Vec3d pa = builder.atoms[ia].pos;
        MM::AtomConf c = builder.confs[a.iconf];
        for(int j=0; j<N_NEIGH_MAX; j++ ){
            int ib = c.neighs[j];
            //printf( "neigh[%i] = %i \n", j, ib );
            if( ib>=0 ){
                int ja = builder.bonds[ib].getNeighborAtom(ia);
                //Draw3D::drawLine( pa, builder.atoms[c.neighs[j]].pos  );
                Draw3D::vertex(pa);
                Draw3D::vertex(builder.atoms[ja].pos);
            }
        }
    }
    glEnd();
}

void drawBonds( const MMFFsp3& ff, double Fsc=0.0 ){
    //drawSystem( false, true, false );
    glBegin( GL_LINES );
    for(int ib=0; ib<ff.nbonds; ib++ ){
        const Vec2i& b = ff.bond2atom[ib]; 
        Draw3D::vertex(ff.apos[b.i]);
        Draw3D::vertex(ff.apos[b.j]);
    }
    glEnd();
}

void drawNeighs( const MMFFsp3& ff, double Fsc=0.0 ){
    //drawSystem( false, true, false );
    for(int ia=0; ia<ff.nnode; ia++ ){
        //printf( "atom[%i]\n", ia );
        int* ngs = ff.aneighs + ia*ff.nneigh_max;
        for(int j=0; j<ff.nneigh_max; j++ ){
            //printf( "atom[%i]neigh[%i]=%i \n", ia, j, ngs[j] );
            if(ngs[j]>=0){
                glColor3f(0.,0.,0.); Draw3D::drawLine( ff.apos[ia], ff.apos[ngs[j]] );
                if(Fsc>0.0){ glColor3f(1.,0.,0.); Draw3D::drawVecInPos( ff.fapos[ia]*Fsc, ff.apos[ia] ); }
            }else{
                int ipi = -ngs[j]-1;
                glColor3f(0.,0.5,0.); Draw3D::drawVecInPos( ff.pipos[ipi], ff.apos[ia] );
                if(Fsc>0.0){ glColor3f(1.,0.5,0.); Draw3D::drawVecInPos( ff.fpipos[ipi]*Fsc, ff.apos[ia]+ff.pipos[ipi] ); }
            }
        }
    }
}

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppMMFFsp3 : public AppSDL2OGL_3D { public:

    int verbosity = 0;

	Molecule    mol;
	MMFFparams  params;
    MMFFsp3    ff;
    NBFF       nff;
    MM::Builder builder;
    DynamicOpt  opt;

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
    int     ogl_sph=0;
    int     ogl_mol=0;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1, iangPicked = -1;
    int perFrame =  50;


    double drndv =  10.0;
    double drndp =  0.5;

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMMFFsp3( int& id, int WIDTH_, int HEIGHT_ );

	int makeMoleculeInline();
	int makeMoleculeInlineBuilder( bool bPBC );
	int loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int loadMoleculeXYZ( const char* fname, bool bAutoH );
	int loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );

	//void drawSystem( );
    void drawSystem( bool bAtoms=true, bool bBonds=true, bool bForces=false );

	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

};

void TestAppMMFFsp3::selectRect( const Vec3d& p0, const Vec3d& p1 ){
    Vec3d Tp0,Tp1,Tp;
    Mat3d rot = (Mat3d)cam.rot;
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

void  TestAppMMFFsp3::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
    //int ib         = builder.rayBonds( ro, rd, 0.3 );
    int ib = pickBondCenter( ff.nbonds, ff.bond2atom, ff.apos, ro, rd, 0.5 );
    printf( "picked bond  %i \n", ib );
    const Vec2i& b = builder.bonds[ib].atoms;
    rotation_axis = (ff.apos[b.b]-ff.apos[b.a]).normalized();
    std::unordered_set<int> innodes1; innodes1.insert( b.a );
    std::unordered_set<int> innodes2; innodes2.insert( b.b );
    MM::splitGraphs( ff.nbonds, ff.bond2atom, ib, innodes1 );
    MM::splitGraphs( ff.nbonds, ff.bond2atom, ib, innodes2 );
    std::unordered_set<int>* sel;
    if( innodes1.size()<innodes2.size()  ){ sel=&innodes1;  rotation_center = ff.apos[b.a]; }else{ sel=&innodes2; rotation_center = ff.apos[b.b]; rotation_axis.mul(-1); }
    selection.clear();
    for( int i:*sel){ selection.push_back(i); };
}


int TestAppMMFFsp3::loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH ){
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);
    int nheavy = builder.load_xyz( fname, bAutoH, true );
    readMatrix( fnameLvs, 3, 3, (double*)&builder.lvec );

    //builder.printAtoms ();
    //builder.printConfs ();
    if(verbosity>2)builder.printAtomConfs();
    //exit(0);
    builder.export_atypes(atypes);

    builder.verbosity = verbosity;
    //builder.autoBonds ();             if(verbosity>2)builder.printBonds ();
    builder.autoBondsPBC();             if(verbosity>2)builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             if(verbosity>2)builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     if(verbosity>2)builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     if(verbosity>2)builder.printAngles();

    //bNonBonded = false;
    //exit(0);
    builder.toMMFFsp3( ff );

    builder.saveMol( "data/polymer.mol" );

    return nheavy;
}


int TestAppMMFFsp3::loadMoleculeMol( const char* fname, bool bAutoH, bool bLoadTypes ){
    if(bLoadTypes){
        printf( "bLoadTypes==True : load atom and bond types from file \n" );
        params.loadAtomTypes( "common_resources/AtomTypes.dat" );
        params.loadBondTypes( "common_resources/BondTypes.dat");
        //builder.params = &params;
    }else{
        printf( "bLoadTypes==False : make default Atom names dictionary \n" );
        params.initDefaultAtomTypeDict( );
    }
    mol.bindParams(&params);
    mol.loadMol( fname );
    int iH = params.atomTypeDict["H"];
    int nh     = mol.countAtomType(iH); printf( "nh %i\n", nh );
    int nheavy = mol.natoms - nh;
    if(bAutoH){
        printf( "bAutoH==True : MMFFBuilder will re-calculate hydrogens, pi-orbitals and free electron pairs \n" );
        builder.insertFlexibleMolecule_ignorH( &mol, Vec3dZero, Mat3dIdentity, iH );
    }else{
        printf( "bAutoH==False : Angles assigned by simple algorithm Molecule::autoAngles \n" );
        //mol.bondsOfAtoms();
        params.assignREs( mol.natoms, mol.atomType, mol.REQs );
        mol.autoAngles(true);
        Vec3d cog = mol.getCOG_av();
        mol.addToPos( cog*-1.0 );
        builder.insertMolecule(&mol, Vec3dZero, Mat3dIdentity, false );
        builder.toMMFFsp3( ff );
    }
    //builder.sortAtomsOfBonds();
    builder.tryAddConfsToAtoms(0, nh);
    builder.tryAddBondsToConfs();
    //for(int i=0; i<nh; i++){ builder.addConfToAtom(i); }
    //builder.tryAddBondsToConfs();
    //if(verbosity>2)mol.printAtomInfo();
    //if(verbosity>2)mol.printAtom2Bond();
    //if(verbosity>2)mol.printAngleInfo();
    if(verbosity>2)builder.printAtoms();
    //if(verbosity>2)builder.printBonds();
    //if(verbosity>2)builder.printAngles();
    //if(verbosity>2)builder.printConfs();
    //bNonBonded = false;      // ToDo : WARRNING : this is just hack, because builder.sortBonds() does not seem to work, we have to switch off Non-Bonding interactions
    builder.trySortBonds();
    //builder.sortBonds();
    if(verbosity>2)builder.printBonds();
    if(verbosity>2)builder.printAngles();
    if(verbosity>2)builder.printConfs();
    builder.toMMFFsp3( ff );
    //Draw3D::shapeInPoss( ogl_sph, ff.natoms, ff.apos, 0 );
    return nheavy;
}

int TestAppMMFFsp3::makeMoleculeInlineBuilder( bool bPBC ){
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
    MM::Dihedral brushDihedral{ -1,   {-1,-1,-1},    3, 0.5 };  println(brushDihedral);
    builder.insertDihedralByAtom( {0,1,2,3}, brushDihedral );
    builder.trySortBonds();

    builder.toMMFFsp3( ff );

    builder.lvec.a = (Vec3d){  5.0,0.0,0.0 };
    builder.lvec.b = (Vec3d){  0.0,5.0,0.0 };
    builder.lvec.c = (Vec3d){  0.0,0.0,5.0 };
    //if(bPBC){    // --- Periodic Boundary Conditions
    //    ff.initPBC();                 // as far as good, pbc-shifts are curenlty zero, so no change
    //    ff.pbcShifts[1] = builder.lvec.a*-1.; // make bond 3 from nighboring cell
    //    ff.printBondParams();
    //}

    return natom;
}

//void TestAppMMFFsp3::makeAtoms(){}
//template<typename T> std::function<T(const T&,const T&         )> F2;

TestAppMMFFsp3::TestAppMMFFsp3( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ------ using molecule from mol-file does not seem to work now - There is some problem with AtomConfs
    // >> This algorithm assumes all atoms with conf precede atoms without confs in the array
    // >>   ERROR in builder.sortBonds() => exit

    printf("DEBUG 1\n");
    int nheavy = 0;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);

    printf("DEBUG 2\n");
    readMatrix( "common_resources/polymer-2.lvs", 3, 3, (double*)&builder.lvec );
    molecules.push_back( new Molecule() ); molecules[0]->atomTypeDict = builder.atomTypeDict; molecules[0]->load_xyz("common_resources/polymer-2.xyz", verbosity );
    molecules.push_back( new Molecule() ); molecules[1]->atomTypeDict = builder.atomTypeDict; molecules[1]->load_xyz("common_resources/polymer-2-monomer.xyz", verbosity );
    builder.insertFlexibleMolecule(  molecules[0], {0,0,0}       , Mat3dIdentity, -1 );
    builder.insertFlexibleMolecule(  molecules[1], builder.lvec.a*1.2, Mat3dIdentity, -1 );

    printf("DEBUG 3\n");
    builder.lvec.a.x *= 2.3;

    //if(verbosity>1)builder.printAtoms ();
    //if(verbosity>1)builder.printConfs ();
    if(verbosity>1)builder.printAtomConfs();
    builder.export_atypes(atypes);
    builder.verbosity = verbosity;
    //builder.autoBonds ();             if(verbosity>1)builder.printBonds ();
    builder.autoBondsPBC();             if(verbosity>1)builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});   if(verbosity>1)builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     if(verbosity>1)builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     if(verbosity>1)builder.printAngles();

    builder.sortConfAtomsFirst();
    builder.makeAllConfsSP();
    builder.printAtomConfs();
    builder.toMMFFsp3( ff );
    builder.saveMol( "data/polymer.mol" );

    ff.printBonds();

    printf("DEBUG 4\n");
    bNonBonded = false;
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

    printf("DEBUG 5\n");

    opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    printf("DEBUG 6\n");
    // ======== Test before we run
    if(verbosity>1)nff.printAtomParams();
    ff.printNeighs();
    
    printf("DEBUG 7 \n");
    ckeckNaN_d( ff.natoms, ff.nneigh_max, ff.Kneighs, "ff.Kneighs" );
    ckeckNaN_d( ff.nbonds,             1, ff.bond_k,  "ff.bond_k"  );

    double E = ff.eval(true);
    printf( "iter0 E = %g \n", E );
    printf("TestAppMMFFsp3.init() DONE \n");
    //exit(0);

    printf("DEBUG 8\n");
    //Draw3D::makeSphereOgl( ogl_sph, 3, 1.0 );
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );

    printf("DEBUG 9\n");
    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    //selection.insert( selection.end(), {12, 16, 14, 6, 2, 3,   20,18,31,25,26} );
    //selection.insert( selection.end(), {13,29,30} );
    //splitGraphs( ff.nbonds, ff.bond2atom, 12, 22 );
    //std::unordered_set<int> innodes;
    //innodes.insert(12);
    //innodes.insert(11);
    //MM::splitGraphs( ff.nbonds, ff.bond2atom, 22, innodes );
    //for( int i : innodes ){ selection.push_back(i); }

    printf("DEBUG 10\n");
}

void TestAppMMFFsp3::drawSystem( bool bAtoms, bool bBonds, bool bForces ){
    if(bBonds){
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLines ( ff.nbonds, (int*)ff.bond2atom, ff.apos );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec ); // DEBUG
        glColor3f(0.0f,0.0f,0.0f);Draw3D::bonds( ff.nbonds, ff.bond2atom, ff.apos ); // DEBUG
        //glColor3f(0.0f,0.0f,1.0f); Draw3D::bondLabels( ff.nbonds, ff.bond2atom, ff.apos, fontTex, 0.02 );                     //DEBUG
    }
    if(bForces){ glColor3f(1.0f,0.0f,0.0f); Draw3D::vecsInPoss( ff.natoms, ff.fapos, ff.apos, 300.0              ); }
    if(bAtoms){
        //glColor3f(0.0f,0.0f,1.0f); Draw3D::atomPropertyLabel( ff.natoms, (double*)nff.REQs, ff.apos, 3,2, fontTex, 0.02, "%4.2f\0" );
        glColor3f(0.5f,0.0f,0.0f); Draw3D::atomLabels( ff.natoms, ff.apos, fontTex, 0.01                     );                     //DEBUG
        //Draw3D::atomsREQ  ( ff.natoms, ff.apos,   nff.REQs, ogl_sph, 1.0, 0.25, 1.0 );
        //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 1.0, 1.0 );       //DEBUG
        //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.5, 1.0 );       //DEBUG
        //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.25, 1.0 );       //DEBUG
    }
}

void TestAppMMFFsp3::draw(){
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( "====== Frame # %i \n", frameCount );
    //cam.qrot.rotate(M_PI*0.01,Vec3fX);

    //glColor3f(0.,0.,0.); drawBonds( builder );
    //glColor3f(0.,0.,0.); drawNeighs( builder );
    //drawNeighs( ff );
    //builder.toMMFFsp3( ff, false );

    if(frameCount==1){
        qCamera.pitch( M_PI );
        //ff.printAtomPos();
        //ff.printBondParams();
        //ff.printAngleParams();
        //ff.printTorsionParams();
        lvec_a0 = builder.lvec.a;
        printf( "lvec_a0  (%g %g,%g) \n", lvec_a0.x, lvec_a0.y, lvec_a0.z );
    }

    //Draw3D::drawAxis(  10. );

	if( ogl_mol ){
        glCallList( ogl_mol );
        return;
        //exit(0);
    }

	//ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);

    double Ftol = 1e-6;
    //double Ftol = 1e-2;

    //ff.apos[0].set(-2.0,0.0,0.0);
    perFrame = 1;
    //perFrame = 100;
    //perFrame = 20;
    //bRunRelax = false;
    bool makeScreenshot = false;

    if(bRunRelax){
        builder.lvec.a    = lvec_a0 + Vec3d{-1.0,0.0,0.0};
        for(int itr=0; itr<perFrame; itr++){
            double E=0;

            // --- Eval Forces
            ff.cleanAtomForce();
            E += ff.eval(true);
            /*
            if( isnan( E) ){ 
                printf("ERROR : E=ff.eval() is NaN  \n"); 
                ff.checkNaNs();
                exit(0); 
            }
            */
            if(bNonBonded){ E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} ); }
            if(ipicked>=0){ Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 ); ff.fapos[ipicked].add( f ); };
            for(int i=0; i<ff.natoms; i++){ ff.fapos[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) ); }
            
            // --- Move
            double f2;
            //opt.move_MD( 0.001, 0.005 );
            opt.move_GD( 0.01 );
            //printf( "E %g |F| %g |Ftol %g \n", E, sqrt(f2), Ftol );
            if(f2<sq(Ftol)){
                bConverged=true;
            }

        }
    }

    //drawSystem();
    //glColor3f(0.,0.,0.); drawBonds( ff );
    drawSystem(true,true,false);
    drawNeighs( ff, 0.0 );

    //glColor3f(0.,0.,0.); drawBonds( builder );
    //glColor3f(0.,0.,0.); drawNeighs( builder );
    //drawNeighs( ff );

    /*
    if(bDragging)Draw3D::drawTriclinicBox(cam.rot.transposed(), (Vec3f)ray0_start, (Vec3f)ray0 );
    //Draw3D::drawTriclinicBox(builder.lvec, Vec3dZero, Vec3dOne );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawTriclinicBox(builder.lvec.transposed(), Vec3dZero, Vec3dOne );
    //glColor3f(0.6f,0.6f,0.6f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    //glColor3f(0.95f,0.95f,0.95f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    if(builder.bPBC){
        //printf( "draw PBC \n" );
        //Draw3D::drawPBC( (Vec3i){1,1,0}, builder.lvec, [&](){drawSystem();} );
        Draw3D::drawPBC( (Vec3i){2,2,0}, builder.lvec, [&](){drawSystem();} );
    }else{
        drawSystem();
    }
    for(int i=0; i<selection.size(); i++){
        int ia = selection[i];
        glColor3f( 0.f,1.f,0.f );
        Draw3D::drawSphereOctLines( 8, 0.3, ff.apos[ia] );
    }
    if(makeScreenshot){
        saveScreenshot( icell );
        icell++;
    }
    */
};


void TestAppMMFFsp3::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_KP_7: builder.lvec.a.x+=xstep; break;
                //case SDLK_KP_4: builder.lvec.a.x-=xstep; break;
                //case SDLK_KP_8: builder.lvec.a.y+=xstep; break;
                //case SDLK_KP_5: builder.lvec.a.y-=xstep; break;
                //case SDLK_KP_9: builder.lvec.a.z+=xstep; break;
                //case SDLK_KP_6: builder.lvec.a.z-=xstep; break;

                case SDLK_KP_7: ff.iDEBUG_pick++; if(ff.iDEBUG_pick>=ff.iDEBUG_n)ff.iDEBUG_pick=0; break;
                case SDLK_KP_4: ff.iDEBUG_pick--; if(ff.iDEBUG_pick<0           )ff.iDEBUG_pick=ff.iDEBUG_n-1; break;

                case SDLK_f:
                    //selectShorterSegment( (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y + cam.rot.c*-1000.0), (Vec3d)cam.rot.c );
                    selectShorterSegment( ray0, (Vec3d)cam.rot.c );
                    //selection.erase();
                    //for(int i:builder.selection){ selection.insert(i); };
                    break;

                case SDLK_LEFTBRACKET:
                    Vec3d::rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, +rotation_step );
                    break;
                case SDLK_RIGHTBRACKET:
                    Vec3d::rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, -rotation_step );
                    break;
                case SDLK_SPACE: bRunRelax=!bRunRelax; break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                    ff.iDEBUG_pick=ipicked;
                    printf( "ipicked %i \n", ipicked );
                    /*
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                    selection.clear();
                    if(ipicked>=0){ selection.push_back(ipicked); };
                    printf( "picked atom %i \n", ipicked );
                    */
                    ray0_start = ray0;
                    bDragging = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = ff.pickBond( ray0, (Vec3d)cam.rot.c , 0.5 );
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
                        ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                        selection.clear();
                        if(ipicked>=0){ selection.push_back(ipicked); };
                        printf( "picked atom %i \n", ipicked );
                    }else{
                        selectRect( ray0_start, ray0 );
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

void TestAppMMFFsp3::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppMMFFsp3 * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	//thisApp = new TestAppMMFFsp3( junk , 800, 600 );
    thisApp = new TestAppMMFFsp3( junk , 1600, 1000 );
	//thisApp = new TestAppMMFFsp3( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















