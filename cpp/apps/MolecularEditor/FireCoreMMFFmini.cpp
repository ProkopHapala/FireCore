
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"
#include "IO_utils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
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
#include "NBFF.h"
#include "GridFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams

#include "FireCoreAPI.h"

//#include "NBSRFF.h"


#include "MarchingCubes.h"
#include "MolecularDraw.h"

#include "AppSDL2OGL_3D.h"

int idebug=0;


// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppMMFFmini : public AppSDL2OGL_3D { public:
	Molecule    mol;
	MMFFparams  params;
    MMFFmini    ff;
    NBFF       nff;
    GridFF     gridFF;
    MM::Builder builder;
    DynamicOpt  opt;
    FireCore::Lib  fireCore;
    FireCore::QMMM qmmm;

    GridShape MOgrid;

    int* atypes = 0;
    int* atypeZ = 0;

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

    bool bDoMM=true,bDoQM=true;
    bool bConverged = false;
    bool bRunRelax  = false;

    // Visualization params
    int which_MO  = 7; 
    double mm_Rsc =  0.25;
    double mm_Rsub = 1.0;
    bool   mm_bAtoms = false;
    bool makeScreenshot = false;
    bool renderType = 1;

    int  fontTex;
    int  ogl_sph=0;
    int  ogl_mol=0;
    int  ogl_isosurf=0;
    int  ogl_MO = 0;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1, iangPicked = -1;
    Vec3d* picked_lvec = 0;
    int perFrame =  50;

    double drndv =  10.0;
    double drndp =  0.5;

    Vec3d testREQ,testPLQ;



	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

    void MDloop();

	TestAppMMFFmini( int& id, int WIDTH_, int HEIGHT_ );

	//int makeMoleculeInline();
	//int makeMoleculeInlineBuilder( bool bPBC );
	int  loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int loadMoleculeXYZ( const char* fname, bool bAutoH );
	int  loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );
    void tryLoadGridFF();
    void makeGridFF   (bool recalcFF=false, bool bRenderGridFF=true);

	void drawSystem( Vec3d ixyz );
    void drawSystemQMMM();
    void renderOrbital(int i, double iso=0.1);
    void renderDensity(       double iso=0.1);

	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

};

//=================================================
//                   INIT()
//=================================================

TestAppMMFFmini::TestAppMMFFmini( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    int nheavy = 0;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);

    readMatrix( "common_resources/polymer-2.lvs", 3, 3, (double*)&builder.lvec );
    molecules.push_back( new Molecule() ); molecules[0]->atomTypeDict = builder.atomTypeDict; molecules[0]->load_xyz("common_resources/polymer-2.xyz", true);
    //molecules.push_back( new Molecule() ); molecules[1]->atomTypeDict = builder.atomTypeDict; molecules[1]->load_xyz("common_resources/polymer-2-monomer.xyz", true);
    builder.insertFlexibleMolecule(  molecules[0], {0,0,0}           , Mat3dIdentity, -1 );
    //builder.insertFlexibleMolecule(  molecules[1], builder.lvec.a*1.2, Mat3dIdentity, -1 );

    builder.lvec.a.x *= 2.3;

    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    builder.export_atypes(atypes);
    builder.verbosity = true;
    //builder.autoBonds ();             builder.printBonds ();
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     builder.printAngles();
    builder.toMMFFmini( ff, &params );
    builder.saveMol( "data/polymer.mol" );

    //builder.lvec.a.x *= 2.0;

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


    DEBUG
    makeGridFF();
    DEBUG

    qmmm.init(6);
    qmmm.params=&params;
    _vec2arr(qmmm.imms,  {4,5,10,11,    8,23} );
    _vec2arr(qmmm.isCap, {0,0, 0, 0,    1, 1} );
    ff.reallocMask();
    //ff.bondMasked[0] = true; printf( "bondMasked[0] %i,%i \n", ff.bond2atom[0].a,ff.bond2atom[0].b );
    qmmm.maskMMFF(ff);
    qmmm.setAtypes( atypes);
    qmmm.load_apos(ff.apos);

    
    fireCore.loadLib( "/home/prokop/git/FireCore/build/libFireCore.so" );
    fireCore.preinit( );
    fireCore.set_lvs( (double*)&(builder.lvec) );
    fireCore.init( qmmm.nqm, qmmm.atypeZ, (double*)qmmm.apos );
    double tmp[3]{0.,0.,0.};
    fireCore.setupGrid( 100.0, 0, tmp, (int*)&MOgrid.n, (double*)&MOgrid.dCell );
    MOgrid.printCell();
    //exit(0);

    qmmm.bindFireCoreLib( fireCore );

    //_list2array(int,qmmm.nqm,#{4,5,10,11,  8,23},qmmm.imms);
    //int imms_[] = ; _forN(i,qmmm.nqm) qmmm.imms[i]=imms_;

    opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.aforce, 0 );
    opt.initOpt( 0.05, 0.9 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    // ======== Test before we run
    /*
    nff.printAtomParams();
    ff.printAtomPos();
    ff.printBondParams();
    ff.printAngleParams();
    ff.printTorsionParams();
    */
    double E = ff.eval(true);
    printf( "iter0 ff.eval() E = %g \n", E );
    //for(int i=0; i<ff.natoms; i++){ printf( "ff.aforce[%i] %g %g %g \n", i, ff.aforce[i].x, ff.aforce[i].y, ff.aforce[i].z ); }
    printf("TestAppMMFFmini.init() DONE \n");
    //exit(0);

    //Draw3D::makeSphereOgl( ogl_sph, 3, 1.0 );
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );

    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    std::unordered_set<int> innodes;
    //innodes.insert(12);
    innodes.insert(11);
    MM::splitGraphs( ff.nbonds, ff.bond2atom, 22, innodes );
    for( int i : innodes ){ selection.push_back(i); }

    picked_lvec = &builder.lvec.a;

    params.printAtomTypeDict();
    printf( " # ==== SETUP DONE ==== \n" );
}


void TestAppMMFFmini::renderOrbital(int iMO, double iso ){
    if(ogl_MO){ glDeleteLists(ogl_MO,1); }
    qmmm.evalQM( ff.apos, ff.aforce );
    int ntot = MOgrid.n.x*MOgrid.n.y*MOgrid.n.z;
    double* ewfaux = new double[ ntot ];
    fireCore.getGridMO( iMO, ewfaux );
    //for(int i=0; i<ntot; i++) ewfaux[i]=0;
    //ewfaux[ntot/2 + MOgrid.n.x*(MOgrid.n.y/2) + MOgrid.n.x/2 ] = 1.0;
    //MOgrid.pos0 = MOgrid.dCell.a*(-MOgrid.n.a/3) + MOgrid.dCell.b*(-MOgrid.n.b/2) + MOgrid.dCell.c*(-MOgrid.n.c/2);
    //MOgrid.pos0 = Vec3dZero; ewfaux[10*MOgrid.n.x*MOgrid.n.y + MOgrid.n.x*10 + 10 ] = 1.0;
    ogl_MO  = glGenLists(1);
    Vec3d p=(Vec3d){0.4,2.5,0.0};
    glNewList(ogl_MO, GL_COMPILE);
    glTranslatef( p.x, p.y, p.z );
    int ntris=0;  
    glColor3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( MOgrid,  iso, ewfaux, renderType  );
    glColor3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( MOgrid, -iso, ewfaux, renderType  );
    //printf( "renderOrbital() ntris %i \n", ntris );
    //int ntris = Draw3D::Grid2Points( MOgrid, 0.1, ewfaux );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawTriclinicBox(builder.lvec.transposed(), Vec3dZero, Vec3dOne );
    glTranslatef( -p.x, -p.y, -p.z );
    glEndList();
    delete [] ewfaux;
}

void TestAppMMFFmini::renderDensity(double iso){
    if(ogl_MO){ glDeleteLists(ogl_MO,1); }
    qmmm.evalQM( ff.apos, ff.aforce );
    int ntot = MOgrid.n.x*MOgrid.n.y*MOgrid.n.z;
    double* ewfaux = new double[ ntot ];
    fireCore.getGridDens( 0, 0, ewfaux );
    ogl_MO  = glGenLists(1);
    glNewList(ogl_MO, GL_COMPILE);
    int ntris = Draw3D::MarchingCubesCross( MOgrid, iso, ewfaux, renderType  );
    //printf( "renderOrbital() ntris %i \n", ntris );
    glEndList();
    delete [] ewfaux;
}




//=================================================
//                   MDloop()
//=================================================
void TestAppMMFFmini::MDloop(){
    //builder.lvec.a    = lvec_a0 + Vec3d{-1.0,0.0,0.0};
    double Ftol = 1e-6;
    //double Ftol = 1e-2;
    //ff.apos[0].set(-2.0,0.0,0.0);
    perFrame = 1;
    //perFrame = 100;
    //perFrame = 20;
    //bRunRelax = false;
    for(int itr=0; itr<perFrame; itr++){
        double E=0;
        ff.cleanAtomForce();
        if(bDoQM){
            qmmm.evalQM      ( ff.apos, ff.aforce );
            qmmm.applyCharges( nff.REQs, true );
        }
        if(bDoMM){
            //qmmm.evalQM      ( ff.apos, ff.aforce );
            //qmmm.applyCharges( nff.REQs );
            E += ff.eval(false);
            if(bNonBonded){
                //E += nff.evalLJQ_sortedMask();   // This is fast but does not work in PBC
                E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
            }
            for(int i=0; i<ff.natoms; i++){
                //ff.aforce[i].add( getForceHamakerPlane( ff.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );
                ff.aforce[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) );
                //ff.aforce[i].z += ff.apos[i].z * K;
                //printf( "%g %g %g\n",  world.aforce[i].x, world.aforce[i].y, world.aforce[i].z );
            }
        }
        if(ipicked>=0){
            float K = -2.0;
            Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, K );
            ff.aforce[ipicked].add( f );
        };
        ff.aforce[  10 ].set(0.0); // This is Hack to stop molecule from moving
        //opt.move_LeapFrog(0.01);
        opt.move_MDquench();
        //opt.move_GD(0.001);
        double f2=1;
        //double f2 = opt.move_FIRE();
        // =========== Molecular Dynamics with random velicity
        //float d = 0.05;
        //for(int i=0; i<ff.natoms; i++){ ff.aforce[i].add({randf(-d,d),randf(-d,d),randf(-d,d)});  };
        //double f2; opt.move_MD( 0.1, 0.005 );
        //printf( "E %g |F| %g |Ftol %g \n", E, sqrt(f2), Ftol );
        if(f2<sq(Ftol)){
            bConverged=true;
        }
    }
}

//=================================================
//                   DRAW()
//=================================================

void TestAppMMFFmini::draw(){
    glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    //glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);
    //glDepthMask(false);
    //glLineWidth(lw);
    //glDrawArrays(GL_LINE_STRIP, offset, count);

    //printf( "====== Frame # %i \n", frameCount );
    //cam.qrot.rotate(M_PI*0.01,Vec3fX);
    //Draw3D::drawAxis(  10. );
    //if( ogl_mol ){ glCallList( ogl_mol ); return; }
    //printf( "builder.lvec: " ); builder.lvec.print();

    if(frameCount==1){
        qCamera.pitch( M_PI );
        //ff.printAtomPos();
        //ff.printBondParams();
        //ff.printAngleParams();
        //ff.printTorsionParams();
        //lvec_a0 = builder.lvec.a;
        //printf( "lvec_a0  (%g %g,%g) \n", lvec_a0.x, lvec_a0.y, lvec_a0.z );
    }

    //bDoQM=1; bDoMM=0;
    if(bRunRelax){ MDloop(); }

    // --- Mouse Interaction / Visualization
    //ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
	ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
    //ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*(HEIGHT-mouse_begin_y));
    Draw3D::drawPointCross( ray0, 0.1 );        // Mouse Cursor 
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0); // Mouse Dragging Visualization
    Vec3d ray0_ = ray0;            ray0_.y=-ray0_.y;
    Vec3d ray0_start_=ray0_start;  ray0_start_.y=-ray0_start_.y;
    if(bDragging)Draw3D::drawTriclinicBoxT(cam.rot, (Vec3f)ray0_start_, (Vec3f)ray0_ );   // Mouse Selection Box

    // ---  Isosurface Rendering ( Molecular Orbital, Density ) 
    //Draw3D::drawTriclinicBox(builder.lvec, Vec3dZero, Vec3dOne );
    if(ogl_MO){ 
        glPushMatrix();
        Vec3d c = builder.lvec.a*-0.5 + builder.lvec.b*-0.5 + builder.lvec.c*-0.5;
        glTranslatef( c.x, c.y, c.z );
          // Molecular Orbital ?
            glColor3f(1.0,1.0,1.0); 
            glCallList(ogl_MO); 
            //return; 
        glPopMatrix();
    }
    //glColor3f(0.6f,0.6f,0.6f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    //glColor3f(0.95f,0.95f,0.95f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    if(ogl_isosurf)viewSubstrate( 2, 2, ogl_isosurf, gridFF.grid.cell.a, gridFF.grid.cell.b, gridFF.shift );

    //printf( "bDoQM %i bDoMM %i \n", bDoQM, bDoMM );
    if(bDoQM)drawSystemQMMM();
    if(bDoMM)if(builder.bPBC){ Draw3D::drawPBC( (Vec3i){2,2,0}, builder.lvec, [&](Vec3d ixyz){drawSystem(ixyz);} ); } else { drawSystem({0,0,0}); }
    for(int i=0; i<selection.size(); i++){ int ia = selection[i];
        glColor3f( 0.f,1.f,0.f ); Draw3D::drawSphereOctLines( 8, 0.3, ff.apos[ia] );     }
    if(iangPicked>=0){
        glColor3f(0.,1.,0.);      Draw3D::angle( ff.ang2atom[iangPicked], ff.ang_cs0[iangPicked], ff.apos, fontTex );
    }
    //if(makeScreenshot){ saveScreenshot( icell ); icell++; }
};

void TestAppMMFFmini::selectRect( const Vec3d& p0, const Vec3d& p1 ){
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

void  TestAppMMFFmini::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
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

void TestAppMMFFmini::makeGridFF( bool recalcFF, bool bRenderGridFF ) {
    //world.substrate.grid.n    = (Vec3i){120,120,200};
    gridFF.loadXYZ  ( "inputs/NaCl_sym.xyz", params );
    gridFF.grid.n    = (Vec3i){60,60,100};
    //world.substrate.grid.n    = (Vec3i){12,12,20};
    gridFF.grid.pos0 = (Vec3d){0.0,0.0,0.0};
    gridFF.loadCell ( "inputs/cel.lvs" );
    //world.gridFF.loadCell ( "inputs/cel_2.lvs" );
    gridFF.grid.printCell();
    //testREQ = (Vec3d){ 2.181, 0.0243442, 0.0}; // Xe
    //genPLQ();
    gridFF.allocateFFs();
    //world.gridFF.evalGridFFs( {0,0,0} );
    //world.gridFF.evalGridFFs( {1,1,1} );
    //world.gridFF.evalGridFFs(int natoms, Vec3d * apos, Vec3d * REQs );
    //bool recalcFF = true;
    //if( recalcFF ){
    gridFF.tryLoad( "data/FFelec.bin", "data/FFPauli.bin", "data/FFLondon.bin" );
    gridFF.shift = (Vec3d){0.0,0.0,-8.0};
    if(bRenderGridFF){
        int iatom = 11;
        testREQ = (Vec3d){ 1.487, 0.0006808, 0.0}; // H
        testPLQ = REQ2PLQ( testREQ, -1.6 );
        //printf( "testREQ   (%g,%g,%g) -> PLQ (%g,%g,%g) \n",        testREQ.x, testREQ.y, testREQ.z, testPLQ.x, testPLQ.y, testPLQ.z   );
        //printf( "aREQs[%i] (%g,%g,%g) -> PLQ (%g,%g,%g) \n", iatom, aREQ[iatom].x, aREQ[iatom].y, aREQ[iatom].z, world.aPLQ[iatom].x, world.aPLQ[iatom].y, world.aPLQ[iatom].z );
        Vec3d * FFtot = new Vec3d[ gridFF.grid.getNtot() ];
        //world.gridFF.evalCombindGridFF_CheckInterp( (Vec3d){ 2.181, 0.0243442, 0.0}, FFtot );
        //saveXSF( "FFtot_z_CheckInterp.xsf", world.gridFF.grid, FFtot, 2, world.gridFF.natoms, world.gridFF.apos, world.gridFF.atypes );
        gridFF.evalCombindGridFF            ( testREQ, FFtot );
        if(idebug>1) saveXSF( "FFtot_z.xsf",  gridFF.grid, FFtot, 2, gridFF.natoms, gridFF.apos, gridFF.atypes );
        ogl_isosurf = glGenLists(1);
        glNewList(ogl_isosurf, GL_COMPILE);
        glShadeModel( GL_SMOOTH );
        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);
        //getIsovalPoints_a( world.gridFF.grid, 0.1, FFtot, iso_points );
        //renderSubstrate( iso_points.size(), &iso_points[0], GL_POINTS );
        //renderSubstrate_( world.gridFF.grid, FFtot, 0.1, true );
        //renderSubstrate_( world.gridFF.grid, FFtot, 0.01, true );
        renderSubstrate_( gridFF.grid, FFtot, gridFF.FFelec, 0.01, true, 0.1);
        Draw3D::drawAxis(1.0);
        glEndList();
        delete [] FFtot;
    }
}


int TestAppMMFFmini::loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH ){
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);
    int nheavy = builder.load_xyz( fname, bAutoH, true );
    readMatrix( fnameLvs, 3, 3, (double*)&builder.lvec );

    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    //exit(0);
    builder.export_atypes(atypes); // NOTE : these are not proton numbers !!!!

    builder.verbosity = 5;
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

int TestAppMMFFmini::loadMoleculeMol( const char* fname, bool bAutoH, bool bLoadTypes ){

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

    /*
    ogl_mol = glGenLists(1);
    glNewList( ogl_mol, GL_COMPILE );
        Draw3D::drawLines( mol.nbonds, (int*)mol.bond2atom, mol.pos );
    glEndList();
    */

    return nheavy;
}

void TestAppMMFFmini::drawSystem( Vec3d ixyz ){
    bool bOrig = (ixyz.x==0)&&(ixyz.y==0)&&(ixyz.z==0);
    //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLines ( ff.nbonds, (int*)ff.bond2atom, ff.apos );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec ); // DEBUG
    if(bOrig&&mm_bAtoms){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels( ff.natoms, ff.apos, fontTex                     ); }                     //DEBUG
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::bondLabels( ff.nbonds, ff.bond2atom, ff.apos, fontTex, 0.02 );                     //DEBUG
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::atomPropertyLabel( ff.natoms, (double*)nff.REQs, ff.apos, 3,2, fontTex, 0.02, "%4.2f\0" );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::vecsInPoss( ff.natoms, ff.aforce, ff.apos, 300.0              );
    //Draw3D::atomsREQ  ( ff.natoms, ff.apos,   nff.REQs, ogl_sph, 1.0, 0.25, 1.0 );
    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 1.0, 1.0 );       //DEBUG
    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.5, 1.0 );       //DEBUG
    Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, mm_Rsc, mm_Rsub );       //DEBUG
}

void TestAppMMFFmini::drawSystemQMMM(){
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    double Rsub = 1.0; 
    double Rsc  = 0.5;
    for(int i=0; i<qmmm.nqm; i++){
        int im = qmmm.imms[i];
        //const AtomType& atyp = params.atypes[atypes[i]];
        int ityp = qmmm.isCap[i]? 0 : 1;
        const AtomType& atyp = params.atypes[ ityp ];
        Draw::setRGB( atyp.color );
        Draw3D::drawShape( ogl_sph, ff.apos[im], Mat3dIdentity*((atyp.RvdW-Rsub)*Rsc) );
    }
    glColor3f(0.5f,0.0f,0.0f); 
    Draw3D::atomPropertyLabel( qmmm.nqm, qmmm.charges, qmmm.apos, 1,0, fontTex );

}

void TestAppMMFFmini::saveScreenshot( int i, const char* fname ){
    //if(makeScreenshot){
        char str[64];
        sprintf( str, fname, i );               // DEBUG
        printf( "save to %s \n", str );
        unsigned int *screenPixels = new unsigned int[WIDTH*HEIGHT*4];  //DEBUG
        glFlush();                                                      //DEBUG
        glFinish();                                                     //DEBUG
        //glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_INT, screenPixels);
        glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, screenPixels);   //DEBUG
        //SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff );   //DEBUG
        SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000 );   //DEBUG
        SDL_SaveBMP(bitmap, str);    //DEBUG
        SDL_FreeSurface(bitmap);
        delete[] screenPixels;
}

void TestAppMMFFmini::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    switch( event.type ){
        case SDL_KEYDOWN :
            //printf( "key: %c \n", event.key.keysym.sym );
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

                case SDLK_KP_1: picked_lvec = &builder.lvec.a; break;
                case SDLK_KP_2: picked_lvec = &builder.lvec.b; break;
                case SDLK_KP_3: picked_lvec = &builder.lvec.c; break;

                case SDLK_KP_7: picked_lvec->x+=xstep; break;
                case SDLK_KP_4: picked_lvec->x-=xstep; break;

                case SDLK_KP_8: picked_lvec->y+=xstep; break;
                case SDLK_KP_5: picked_lvec->y-=xstep; break;

                case SDLK_KP_9: picked_lvec->z+=xstep; break;
                case SDLK_KP_6: picked_lvec->z-=xstep; break;

                case SDLK_COMMA:  which_MO--; printf("which_MO %i \n", which_MO ); break;
                case SDLK_PERIOD: which_MO++; printf("which_MO %i \n", which_MO ); break;
                //case SDLK_LESS:    which_MO--; printf("which_MO %i \n"); break;
                //case SDLK_GREATER: which_MO++; printf("which_MO %i \n"); break;

                case SDLK_m: renderOrbital( which_MO ); break;
                case SDLK_r: renderDensity(          ); break;
                case SDLK_c: saveScreenshot( frameCount ); break;

                case SDLK_f:
                    //selectShorterSegment( (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y + cam.rot.c*-1000.0), (Vec3d)cam.rot.c );
                    selectShorterSegment( ray0, (Vec3d)cam.rot.c );
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

void TestAppMMFFmini::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppMMFFmini * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMMFFmini( junk , 800, 600 );
	//thisApp = new TestAppMMFFmini( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















