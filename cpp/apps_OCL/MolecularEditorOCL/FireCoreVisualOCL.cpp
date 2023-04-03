int verbosity = 0;
int idebug=0;

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
#include "quaternion.h"
#include "Vec3Utils.h"

#include "raytrace.h"
#include "Forces.h"

#include "MMFFparams.h"
static MMFFparams* params_glob=0;

#include "Molecule.h"
#include "MMFFmini.h"
#include "NBFF_old.h"
#include "GridFF.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

#include "OCL_DFT.h"
#include "OCL_PP.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams
#include "FireCoreAPI.h"
//#include "NBSRFF.h"

#include "MarchingCubes.h"
#include "MolecularDraw.h"
#include "GUI.h"
#include "EditorGizmo.h"
#include "SimplexRuler.h"
#include "AppSDL2OGL_3D.h"



// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

//void command_example(double x, void* caller);

class TestAppFireCoreVisual : public AppSDL2OGL_3D { public:
    // ---- Simulations objects
	Molecule    mol;
	MMFFparams  params;
    MMFFmini    ff;
    NBFF_old    nff;
    GridFF      gridFF;
    OCL_PP      ocl;
    MM::Builder builder;
    FireCore::Lib  fireCore;
    FireCore::QMMM qmmm;
    DynamicOpt  opt;
    GridShape MOgrid;

    int* atypes = 0;
    int* atypeZ = 0;

    bool bNonBonded = true;

    // ---- Interaction / Editation variables
    std::vector<int> selection;
    
    Vec3d  rotation_center = Vec3dZero;
    Vec3d  rotation_axis   = Vec3dZ;
    double rotation_step   = 15.0 * (M_PI/180.0);

    bool   bDragging = false;
    Vec3d  ray0_start;
    Vec3d  ray0;
    int ipicked    = -1; // picket atom 
    int ibpicked   = -1; // picket bond
    int iangPicked = -1; // picket angle
    Vec3d* picked_lvec = 0;
    int perFrame =  1;
    Quat4f qCamera0;

    bool bDoMM=true,bDoQM=true;
    bool bConverged = false;
    bool bRunRelax  = false;
    double cameraMoveSpeed = 1.0;
    bool useGizmo=true;
    bool bDrawHexGrid=true;
    bool bHexDrawing=false; 

    bool bPrepared_mm = false;
    bool bPrepared_qm = false;

    GUI gui;
    EditorGizmo  gizmo;
    SimplexRuler ruler; // Helps paiting organic molecules

    // ---- Visualization params
    int which_MO  = 7; 
    double mm_Rsc =  0.25;
    double mm_Rsub = 1.0;
    bool   mm_bAtoms = false;
    bool   isoSurfRenderType = 1;
    Quat4d testREQ,testPLQ;

    // ---- Graphics objects
    int  fontTex,fontTex3D;
    int  ogl_sph=0;
    int  ogl_mol=0;
    int  ogl_isosurf=0;
    int  ogl_MO = 0;

    char str[256];

    // --------- Functions 

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

    void MDloop();

	TestAppFireCoreVisual( int& id, int WIDTH_, int HEIGHT_ );

	//int  loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int  loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );
    void tryLoadGridFF();
    void makeGridFF   (bool recalcFF=false, bool bRenderGridFF=true);

	void drawSystem( Vec3d ixyz );
    void drawSystemQMMM();
    void renderOrbital(int i, double iso=0.1);
    void renderDensity(       double iso=0.1);

	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

    void InitQMMM();
    void loadGeom();
    void initGUI();
    void drawingHex(double z0);
};

//=================================================
//                   INIT()
//=================================================

void TestAppFireCoreVisual::InitQMMM(){
    // ----- QMMM setup
    qmmm.init(6);
    qmmm.params=&params;
    params_glob = &params;
    _vec2arr(qmmm.imms,  {4,5,10,11,    8,23} );
    _vec2arr(qmmm.isCap, {0,0, 0, 0,    1, 1} );
    ff.reallocMask();
    //ff.bondMasked[0] = true; printf( "bondMasked[0] %i,%i \n", ff.bond2atom[0].a,ff.bond2atom[0].b );
    qmmm.maskMMFF(ff);
    qmmm.setAtypes( atypes);
    qmmm.load_apos(ff.apos);
    // ----- FireCore setup
    //fireCore.loadLib( "/home/prokop/git/FireCore/build/libFireCore.so" );
    fireCore.loadLib( "/home/prokophapala/git/FireCore/build/libFireCore.so" );
    fireCore.preinit( );
    fireCore.set_lvs( (double*)&(builder.lvec) );
    fireCore.init( qmmm.nqm, qmmm.atypeZ, (double*)qmmm.apos );
    double tmp[3]{0.,0.,0.};
    fireCore.setupGrid( 100.0, 0, tmp, (int*)&MOgrid.n, (double*)&MOgrid.dCell );
    MOgrid.printCell();
    qmmm.bindFireCoreLib( fireCore );
}

void TestAppFireCoreVisual::loadGeom(){
    // ---- Load & Build Molecular Structure
    readMatrix( "cel.lvs", 3, 3, (double*)&builder.lvec );
    builder.insertFlexibleMolecule(  builder.loadMolType( "mm.xyz", "polymer1" ), {0,0,0}, Mat3dIdentity, -1 );
    //builder.lvec.a.x *= 2.3;
    builder.printAtomConfs();
    builder.export_atypes(atypes);
    builder.verbosity = true;
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    builder.autoAngles( 10.0, 10.0 );     builder.printAngles();
    builder.toMMFFmini( ff, &params );
    builder.saveMol( "builder_output.mol" );

    // ----- Non-bonded interactions setup 
    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.aforce, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );
    if(bNonBonded){
        if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
            printf( "ERROR: nff.pairMask is not sorted => exit \n" );
            exit(0);
        };
    }else{
        printf( "WARRNING : we ignore non-bonded interactions !!!! \n" );
    }
}

void TestAppFireCoreVisual::initGUI(){
    GUI_stepper ylay;
    ylay.step(6);
    Table* tab1 = new Table( 9, sizeof(builder.lvec.a), (char*)&builder.lvec );
    tab1->addColum( &(builder.lvec.a.x), 1, DataType::Double    );
    tab1->addColum( &(builder.lvec.a.y), 1, DataType::Double    );
    tab1->addColum( &(builder.lvec.a.z), 1, DataType::Double    );
    
    ((TableView*)gui.addPanel( new TableView( tab1, "lattice", 5, ylay.x0,  0, 0, 3, 3 ) ))->input = new GUITextInput();

    ylay.step(2); 
    ((GUIPanel*)gui.addPanel( new GUIPanel( "Zoom: ", 5,ylay.x0,5+100,ylay.x1, true, true ) ) )
        ->setRange(5.0,50.0)
        ->setValue(zoom)
        //->command = [&](GUIAbstractPanel* p){ zoom = ((GUIPanel *)p)->value; return 0; };
        ->setCommand( [&](GUIAbstractPanel* p){ zoom = ((GUIPanel *)p)->value; return 0; } );

    ylay.step(2); 
    ((DropDownList*)gui.addPanel( new DropDownList("Pick Mode:",5,ylay.x0,5+100, 3 ) ) )
        ->addItem("pick_atoms")
        ->addItem("pick_bonds")
        ->addItem("Item_angles");

    ylay.step(6); 
    ((DropDownList*)gui.addPanel( new DropDownList("View Side",5,ylay.x0,5+100, 3 ) ) )
        ->addItem("Top")
        ->addItem("Botton")
        ->addItem("Front")
        ->addItem("Back")
        ->addItem("Left")
        ->addItem("Right")
        ->setCommand( [&](GUIAbstractPanel* me_){ 
            DropDownList& me = *(DropDownList*)me_;
            printf( "old qCamera(%g,%g,%g,%g) -> %s \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w, me.labels[me.iSelected].c_str()  );
            switch(me.iSelected){
                case 0: qCamera=Quat4fTop;    break;
                case 1: qCamera=Quat4fBotton; break;
                case 2: qCamera=Quat4fFront;  break;
                case 3: qCamera=Quat4fBack;   break;
                case 4: qCamera=Quat4fLeft;   break;
                case 5: qCamera=Quat4fRight;  break;
            }
            printf( "->new qCamera(%g,%g,%g,%g) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
            qCamera.toMatrix(cam.rot);
            printf( "cam: aspect %g zoom %g \n", cam.aspect, cam.zoom);
            printMat((Mat3d)cam.rot);
            }
        );
}


TestAppFireCoreVisual::TestAppFireCoreVisual( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ---- Load Atomic Type Parameters
    int nheavy = 0;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);

    if( file_exist("cel.lvs")        ){ 
        loadGeom(); 
        makeGridFF();
        // ----- Optimizer setup
        opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.aforce, 0 );
        opt.initOpt( 0.3, 0.95 );
        opt.cleanVel( );
        double E = ff.eval(true);
        printf( "iter0 ff.eval() E = %g \n", E );
        bPrepared_mm=true; 
        printf("... MM preparation DONE \n");
    }
    if( file_exist("Fdata/info.dat") ){ 
        InitQMMM(); bPrepared_qm=true; 
        printf("... QM preparation DONE \n");
    }

    ocl.initPP( "common_resources/cl" );

    picked_lvec = &builder.lvec.a;

    // ---- Graphics setup
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );
    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    // ========== Gizmo
    cam.persp = false;
    gizmo.cam = &cam;
    gizmo.bindPoints(ff.natoms, ff.apos      );
    gizmo.bindEdges (ff.nbonds, ff.bond2atom );
    gizmo.pointSize = 0.5;
    //gizmo.iDebug    = 2;

    ruler.setStep( 1.5 * sqrt(3) );

    initGUI();

}

void TestAppFireCoreVisual::renderOrbital(int iMO, double iso ){
    if(ogl_MO){ glDeleteLists(ogl_MO,1); }
    qmmm.evalQM( ff.apos, ff.aforce );
    int ntot = MOgrid.n.x*MOgrid.n.y*MOgrid.n.z;
    double* ewfaux = new double[ ntot ];
    fireCore.getGridMO( iMO, ewfaux );
    ogl_MO  = glGenLists(1);
    Vec3d p=Vec3d{0.4,2.5,0.0};
    glNewList(ogl_MO, GL_COMPILE);
    glTranslatef( p.x, p.y, p.z );
    int ntris=0;  
    glColor3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( MOgrid,  iso, ewfaux, isoSurfRenderType  );
    glColor3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( MOgrid, -iso, ewfaux, isoSurfRenderType  );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawTriclinicBox(builder.lvec.transposed(), Vec3dZero, Vec3dOne );
    glTranslatef( -p.x, -p.y, -p.z );
    glEndList();
    delete [] ewfaux;
}

void TestAppFireCoreVisual::renderDensity(double iso){
    if(ogl_MO){ glDeleteLists(ogl_MO,1); }
    qmmm.evalQM( ff.apos, ff.aforce );
    int ntot = MOgrid.n.x*MOgrid.n.y*MOgrid.n.z;
    double* ewfaux = new double[ ntot ];
    fireCore.getGridDens( 0, 0, ewfaux );
    ogl_MO  = glGenLists(1);
    glNewList(ogl_MO, GL_COMPILE);
    int ntris = Draw3D::MarchingCubesCross( MOgrid, iso, ewfaux, isoSurfRenderType  );
    //printf( "renderOrbital() ntris %i \n", ntris );
    glEndList();
    delete [] ewfaux;
}

//=================================================
//                   MDloop()
//=================================================

void TestAppFireCoreVisual::MDloop(){
    double Ftol = 1e-6;
    for(int itr=0; itr<perFrame; itr++){
        double E=0;
        ff.cleanAtomForce();
        if(bDoQM){
            qmmm.evalQM      ( ff.apos, ff.aforce );
            qmmm.applyCharges( nff.REQs, true );
        }
        if(bDoMM){
            E += ff.eval(false);
            if(bNonBonded){
                E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
            }
            for(int i=0; i<ff.natoms; i++){
                ff.aforce[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) );
            }
        }
        if(ipicked>=0){
            float K = -2.0;
            Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, K );
            ff.aforce[ipicked].add( f );
        };
        ff.aforce[  10 ].set(0.0); // This is Hack to stop molecule from moving
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

//=================================================
//                   DRAW()
//=================================================

void TestAppFireCoreVisual::draw(){
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);

    if(frameCount==1){ qCamera.pitch( M_PI );  qCamera0=qCamera; }
    if(bRunRelax){ MDloop(); }

    // --- Mouse Interaction / Visualization
	ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
    Draw3D::drawPointCross( ray0, 0.1 );        // Mouse Cursor 
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0); // Mouse Dragging Visualization
    Vec3d ray0_ = ray0;            ray0_.y=-ray0_.y;
    Vec3d ray0_start_=ray0_start;  ray0_start_.y=-ray0_start_.y;
    if(bDragging)Draw3D::drawTriclinicBoxT(cam.rot, (Vec3f)ray0_start_, (Vec3f)ray0_ );   // Mouse Selection Box

    if(ogl_MO){ 
        glPushMatrix();
        Vec3d c = builder.lvec.a*-0.5 + builder.lvec.b*-0.5 + builder.lvec.c*-0.5;
        glTranslatef( c.x, c.y, c.z );
            glColor3f(1.0,1.0,1.0); 
            glCallList(ogl_MO); 
        glPopMatrix();
    }
    if(ogl_isosurf)viewSubstrate( 2, 2, ogl_isosurf, gridFF.grid.cell.a, gridFF.grid.cell.b, gridFF.shift );
    if(bDoQM)drawSystemQMMM();
    if(bDoMM)if(builder.bPBC){ Draw3D::drawPBC( (Vec3i){2,2,0}, builder.lvec, [&](Vec3d ixyz){drawSystem(ixyz);} ); } else { drawSystem({0,0,0}); }
    for(int i=0; i<selection.size(); i++){ int ia = selection[i];
        glColor3f( 0.f,1.f,0.f ); Draw3D::drawSphereOctLines( 8, 0.3, ff.apos[ia] );     }
    if(iangPicked>=0){
        glColor3f(0.,1.,0.);      Draw3D::angle( ff.ang2atom[iangPicked], ff.ang_cs0[iangPicked], ff.apos, fontTex3D );
    }

    if(useGizmo){
        gizmo.draw();
    }
    if(bHexDrawing)drawingHex(5.0);
};

void TestAppFireCoreVisual::drawHUD(){
    glDisable ( GL_LIGHTING );
    gui.draw();
}

void TestAppFireCoreVisual::drawingHex(double z0){
    Vec2i ip; Vec2d dp;
    Vec3d p3 = rayPlane_hit( ray0, (Vec3d)cam.rot.c, {0.0,0.0,1.0}, {0.0,0.0,z0} );
    Vec2d p{p3.x,p3.y};
    double off=1000.0;
    bool s = ruler.simplexIndex( p+(Vec2d){off,off}, ip, dp );
    //ruler.nodePoint( ip, p );    glColor3f(1.,1.,1.); Draw3D::drawPointCross(  {p.x,p.y, 5.0}, 0.5 );
    if(s){glColor3f(1.,0.2,1.);}else{glColor3f(0.2,1.0,1.);}
    ruler.tilePoint( ip, s, p ); Draw3D::drawPointCross(  {p.x-off,p.y-off, z0}, 0.2 );
    
    bool bLine=true;
    if(bDrawHexGrid){
        if(bLine){glBegin(GL_LINES);}else{glBegin(GL_POINTS);}
        ruler.simplexIndex( (Vec2d){off,off}, ip, dp );
        double sc = ruler.step/sqrt(3.0);
        for(int ix=0;ix<10;ix++ ){
            for(int iy=0;iy<10;iy++ ){
                Vec2i ip_{ip.x+ix,ip.y+iy};
                ruler.tilePoint( ip_, true,  p ); 
                p.sub(off,off);
                if(bLine){
                    glColor3f(1.0,0.2,1.0); 
                    Vec2d p2;
                    Draw3D::vertex(Vec3f{p.x,p.y,z0}); p2=p+ruler.lvecs[0]*sc; Draw3D::vertex(Vec3f{p2.x,p2.y,z0});
                    Draw3D::vertex(Vec3f{p.x,p.y,z0}); p2=p+ruler.lvecs[1]*sc; Draw3D::vertex(Vec3f{p2.x,p2.y,z0});
                    Draw3D::vertex(Vec3f{p.x,p.y,z0}); p2=p+ruler.lvecs[2]*sc; Draw3D::vertex(Vec3f{p2.x,p2.y,z0});
                }else{
                    glColor3f(1.0,0.2,1.0); Draw3D::vertex(Vec3f{p.x,p.y,z0}); ruler.tilePoint( ip_, false, p );  p.add(off,off);
                    glColor3f(0.2,1.0,1.0); Draw3D::vertex(Vec3f{p.x,p.y,z0});
                }
            }
        }
        glEnd();
    }
}

void TestAppFireCoreVisual::selectRect( const Vec3d& p0, const Vec3d& p1 ){
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

void  TestAppFireCoreVisual::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
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

void TestAppFireCoreVisual::makeGridFF( bool recalcFF, bool bRenderGridFF ) {
    //gridFF.loadXYZ  ( "inputs/NaCl_sym.xyz", params );
    int ret = params.loadXYZ( "inputs/NaCl_sym.xyz", gridFF.natoms, &gridFF.apos, &gridFF.aREQs, &gridFF.atypes );
    gridFF.grid.n    = Vec3i{60,60,100};
    gridFF.grid.pos0 = Vec3d{0.0,0.0,0.0};
    gridFF.loadCell ( "inputs/cel.lvs" );
    gridFF.grid.printCell();
    gridFF.allocateFFs();
    gridFF.tryLoad( "data/FFelec.bin", "data/FFPauli.bin", "data/FFLondon.bin" );
    gridFF.shift = Vec3d{0.0,0.0,-8.0};
    if(bRenderGridFF){
        int iatom = 11;
        testREQ = Quat4d{ 1.487, 0.0006808, 0., 0.}; // H
        testPLQ = REQ2PLQ( testREQ, -1.6 );
        Quat4f * FFtot = new Quat4f[ gridFF.grid.getNtot() ];
        gridFF.evalCombindGridFF            ( testREQ, FFtot );
        //if(idebug>1) saveXSF( "FFtot_z.xsf",  gridFF.grid, FFtot, 2, gridFF.natoms, gridFF.apos, gridFF.atypes );
        if(idebug>1)  gridFF.grid.saveXSF( "FFtot_E.xsf", (float*)FFtot, 4,3, gridFF.natoms, gridFF.atypes, gridFF.apos );
        ogl_isosurf = glGenLists(1);
        glNewList(ogl_isosurf, GL_COMPILE);
        glShadeModel( GL_SMOOTH );
        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);
        renderSubstrate_( gridFF.grid, FFtot, gridFF.FFelec, 0.01, true, 0.1);
        Draw3D::drawAxis(1.0);
        glEndList();
        delete [] FFtot;
    }
}

void TestAppFireCoreVisual::drawSystem( Vec3d ixyz ){
    bool bOrig = (ixyz.x==0)&&(ixyz.y==0)&&(ixyz.z==0);
    //glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec ); // DEBUG
    glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, ff.pbcShifts ); // DEBUG
    if(bOrig&&mm_bAtoms){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels( ff.natoms, ff.apos, fontTex3D                     ); }                     //DEBUG
    Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, mm_Rsc, mm_Rsub );       //DEBUG
}

void TestAppFireCoreVisual::drawSystemQMMM(){
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    double Rsub = 1.0; 
    double Rsc  = 0.5;
    for(int i=0; i<qmmm.nqm; i++){
        int im = qmmm.imms[i];
        int ityp = qmmm.isCap[i]? 0 : 1;
        const AtomType& atyp = params.atypes[ ityp ];
        Draw::setRGB( atyp.color );
        Draw3D::drawShape( ogl_sph, ff.apos[im], Mat3dIdentity*((atyp.RvdW-Rsub)*Rsc) );
    }
    glColor3f(0.5f,0.0f,0.0f); 
    Draw3D::atomPropertyLabel( qmmm.nqm, qmmm.charges, qmmm.apos, 1,0, fontTex3D );

}

void TestAppFireCoreVisual::saveScreenshot( int i, const char* fname ){
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

void TestAppFireCoreVisual::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    gui.onEvent( mouseX, mouseY, event );
    Vec2f pix = ((Vec2f){ 2*mouseX/float(HEIGHT) - ASPECT_RATIO,
                          2*mouseY/float(HEIGHT) - 1              });// *(1/zoom);
    if(useGizmo)gizmo.onEvent( pix, event );
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN :
            //printf( "key: %c \n", event.key.keysym.sym );
            if(gui.bKeyEvents) switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                //case SDLK_LEFTBRACKET:  ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;
                //case SDLK_RIGHTBRACKET: ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;

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

                case SDLK_KP_0: qCamera = qCamera0; break;

                case SDLK_COMMA:  which_MO--; printf("which_MO %i \n", which_MO ); break;
                case SDLK_PERIOD: which_MO++; printf("which_MO %i \n", which_MO ); break;
                //case SDLK_LESS:    which_MO--; printf("which_MO %i \n"); break;
                //case SDLK_GREATER: which_MO++; printf("which_MO %i \n"); break;

                case SDLK_m: renderOrbital( which_MO ); break;
                case SDLK_r: renderDensity(          ); break;
                case SDLK_c: saveScreenshot( frameCount ); break;

                case SDLK_g: useGizmo=!useGizmo; break;
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

                case SDLK_d: {
                    printf( "DEBUG Camera Matrix\n");
                    printf( "DEBUG qCamera(%g,%g,%g,%g) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
                    qCamera.toMatrix(cam.rot);
                    printf( "DEBUG cam aspect %g zoom %g \n", cam.aspect, cam.zoom);
                    printMat((Mat3d)cam.rot);
                } break;

                //case SDLK_g: iangPicked=(iangPicked+1)%ff.nang;
                //    printf( "ang[%i] cs(%g,%g) k %g (%i,%i,%i)\n", iangPicked, ff.ang_cs0[iangPicked].x, ff.ang_cs0[iangPicked].y, ff.ang_k[iangPicked],
                //        ff.ang2atom[iangPicked].a,ff.ang2atom[iangPicked].b,ff.ang2atom[iangPicked].c );
                //    break;

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
        case SDL_WINDOWEVENT:{
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE:{  quit(); }break;
            }   } break;
    } // switch( event.type ){
    //AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void  TestAppFireCoreVisual::keyStateHandling( const Uint8 *keys ){
    double dstep=0.025;
    //if( keys[ SDL_SCANCODE_X ] ){ cam.pos.z +=0.1; }
    //if( keys[ SDL_SCANCODE_Z ] ){ cam.pos.z -=0.1; }
    //if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
	//if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
    //if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rot.a, -cameraMoveSpeed ); }
	//if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rot.a,  cameraMoveSpeed ); }
    //if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rot.b,  cameraMoveSpeed ); }
	//if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rot.b, -cameraMoveSpeed ); }
    //if( keys[ SDL_SCANCODE_Q ] ){ cam.pos.add_mul( cam.rot.c, -cameraMoveSpeed ); }
	//if( keys[ SDL_SCANCODE_E ] ){ cam.pos.add_mul( cam.rot.c,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_LEFT  ] ){ cam.pos.add_mul( cam.rot.a, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ cam.pos.add_mul( cam.rot.a,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_UP    ] ){ cam.pos.add_mul( cam.rot.b,  cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ cam.pos.add_mul( cam.rot.b, -cameraMoveSpeed ); }
    //AppSDL2OGL_3D::keyStateHandling( keys );
};

//void command_example(double x, void* caller){
//    TestAppFireCoreVisual* app = (TestAppFireCoreVisual*)caller;
//    app->zoom = x;
//};

// ===================== MAIN

TestAppFireCoreVisual * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;

    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);
	thisApp = new TestAppFireCoreVisual( junk, DM.w-100, DM.h-100 );
	//thisApp = new TestAppFireCoreVisual( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















