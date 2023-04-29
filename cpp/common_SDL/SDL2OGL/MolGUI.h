
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
//#include "Solids.h"

#include "MolWorld_sp3.h"

#include "raytrace.h"
#include "Draw3D_Molecular.h"  // it needs to know MMFFparams
#include "MolecularDraw.h"
#include "MarchingCubes.h"
#include "GUI.h"
#include "EditorGizmo.h"
#include "SimplexRuler.h"
#include "AppSDL2OGL_3D.h"

#include <chrono>



// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

//void command_example(double x, void* caller);

class MolGUI : public AppSDL2OGL_3D { public:

    Vec3d  rotation_center = Vec3dZero;
    Vec3d  rotation_axis   = Vec3dZ;
    double rotation_step   = 15.0 * (M_PI/180.0);

    bool   bDragging = false;
    Vec3d  ray0_start;
    Vec3d  ray0;
    //int ipicked    = -1; // picket atom 
    //int ibpicked   = -1; // picket bond
    //int iangPicked = -1; // picket angle
    //Vec3d* picked_lvec = 0;
    int perFrame =  1;
    Quat4f qCamera0;

    bool bDoMM=true;//,bDoQM=true;
    bool bConverged = false;
    bool bRunRelax  = false;
    double cameraMoveSpeed = 1.0;
    //bool useGizmo=true;
    bool useGizmo=false;
    bool bDrawHexGrid=true;
    bool bHexDrawing=false; 

    bool bWriteOptimizerState = true;
    //bool bPrepared_mm = false;
    //bool bPrepared_qm = false;

    double z0_scan = 0.0;

    MolWorld_sp3* W=0;

    GUI gui;
    GUIPanel* Qpanel;
    EditorGizmo  gizmo;
    SimplexRuler ruler; // Helps paiting organic molecules

    // ---- Visualization params
    int iSystemCur = 0;
    int which_MO  = 7; 
    //double ForceViewScale = 1.0;
    //double mm_Rsc         = 0.25;
    //double mm_Rsub        = 1.0;

    double ForceViewScale = 100.0;
    double mm_Rsc         = 0.05;
    double mm_Rsub        = 0.0;

    bool   mm_bAtoms        = true;
    bool   bViewMolCharges  = false;
    bool   bViewAtomLabels  = true;
    bool   bViewAtomTypes  = false;
    bool   bViewBondLabels  = false;
    bool   bViewAtomSpheres = true;
    bool   bViewAtomForces  = true;
    bool   bViewPis         = true;
    bool   bViewSubstrate   = true;
    bool   isoSurfRenderType = 1;
    bool bDebug_scanSurfFF = false;
    Quat4d testREQ;
    Quat4f testPLQ;


    Mat3d dlvec{ 0.1,0.0,0.0,   0.0,0.0,0.0,  0.0,0.0,0.0 };
    Vec2f mouse_pix;

    // ----- Visualization Arrays - allows to switch between forcefields, and make it forcefield independnet
    int    natoms=0,nnode=0,nbonds=0;
    int*   atypes;
    Vec2i* bond2atom=0; 
    Vec3d* pbcShifts=0; 
    Vec3d* apos     =0;
    Vec3d* fapos    =0;
    Vec3d* pipos    =0;
    Vec3d* fpipos   =0;
    Quat4d* REQs     =0;

    Quat4i* neighs    = 0;
    Quat4i* neighCell = 0;

    int nsys=0,nvec=0;
    int    * M_neighs    =0;
    int    * M_neighCell =0;
    Quat4f * M_apos      =0;

    // ---- Graphics objects
    int  fontTex,fontTex3D;

    int  ogl_esp=0;
    int  ogl_sph=0;
    int  ogl_mol=0;
    int  ogl_isosurf=0;
    int  ogl_MO = 0;

    static const int nmaxstr=2048;
    char str[nmaxstr];

    std::vector<Quat4f> debug_ps;
    std::vector<Quat4f> debug_fs;
    std::vector<Quat4f> debug_REQs;

    std::vector<Vec2i> bondsToShow;
    Vec3d * bondsToShow_shifts = 0; 

    enum class Gui_Mode { base, edit, scan };
    Gui_Mode  gui_mode = Gui_Mode::base;
    //int gui_mode = Gui_Mode::edit;

    // ======================= Functions 

	virtual void draw   () override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( ) override;
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys )       override;

    void eventMode_default( const SDL_Event& event );
    void eventMode_scan   ( const SDL_Event& event );
    void eventMode_edit   ( const SDL_Event& event );
    void mouse_default( const SDL_Event& event );

	MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_ );
    void init();

	//int  loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int  loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );
    void tryLoadGridFF();
    //void makeGridFF   (bool recalcFF=false, bool bRenderGridFF=true);
    void renderGridFF( double isoVal=0.001, int isoSurfRenderType=0, double colorScale = 50. );
    void renderESP( Quat4d REQ=Quat4d{ 1.487, 0.02609214441, 1., 0.} );

    void bindMolecule(int natoms_, int nnode_, int nbonds_, int* atypes_,Vec3d* apos_,Vec3d* fapos_,Quat4d* REQs_, Vec3d* pipos_, Vec3d* fpipos_, Vec2i* bond2atom_, Vec3d* pbcShifts_);
	void drawSystem    ( Vec3i ixyz=Vec3iZero );
    void drawPi0s( float sc );
    void  showAtomGrid( char* s, int ia, bool bDraw=true );
    Vec3d showNonBond ( char* s, Vec2i b, bool bDraw=true );
    void  showBonds();
    void printMSystem( int isys, int perAtom, int na, int nvec, bool bNg=true, bool bNgC=true, bool bPos=true );
    //void flipPis( Vec3d ax );
    //void drawSystemQMMM();
    //void renderOrbital(int i, double iso=0.1);
    //void renderDensity(       double iso=0.1);
	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );
    void debug_scanSurfFF( int n, Vec3d p0, Vec3d p1, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.2, 0.}, double sc=NAN );

    //void InitQMMM();
    void initGUI();
    void initWiggets();
    void drawingHex(double z0);
};

//=================================================
//                   INIT()
//=================================================


void MolGUI::initWiggets(){

    GUI_stepper ylay;
    ylay.step(2);
    Qpanel = new GUIPanel( "Q_pick: ", 5,ylay.x0,5+100,ylay.x1, true, true ); 
    Qpanel->setRange(-1.0,1.0)
          ->setValue(0.0)
        //->command = [&](GUIAbstractPanel* p){ zoom = ((GUIPanel *)p)->value; return 0; };
          ->setCommand( [&](GUIAbstractPanel* p){ W->nbmol.REQs[W->ipicked].z = ((GUIPanel *)p)->value; return 0; } );
    (GUIPanel*)gui.addPanel( Qpanel );

    ylay.step(6);
    Table* tab1 = new Table( 9, sizeof(W->builder.lvec.a), (char*)&W->builder.lvec );
    tab1->addColum( &(W->builder.lvec.a.x), 1, DataType::Double    );
    tab1->addColum( &(W->builder.lvec.a.y), 1, DataType::Double    );
    tab1->addColum( &(W->builder.lvec.a.z), 1, DataType::Double    );
    
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
                case 0: qCamera=qTop;    break;
                case 1: qCamera=qBottom; break;
                case 2: qCamera=qFront;  break;
                case 3: qCamera=qBack;   break;
                case 4: qCamera=qLeft;   break;
                case 5: qCamera=qRight;  break;
            }
            printf( "->new qCamera(%g,%g,%g,%g) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
            qCamera.toMatrix(cam.rot);
            printf( "cam: aspect %g zoom %g \n", cam.aspect, cam.zoom);
            printMat((Mat3d)cam.rot);
            }
        );
}

MolGUI::MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    long T0=getCPUticks();
    float nseconds = 0.1;
    SDL_Delay( (int)(1000*0.1) );
    tick2second = nseconds/(getCPUticks()-T0);
    printf( "CPU speed calibration: tick=%g [s] ( %g GHz)\n", tick2second, 1.0e-9/tick2second );

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    if(W_==0){ W = new MolWorld_sp3(); }else{ W=W_; }
    W->tmpstr=str;

}

void MolGUI::initGUI(){
    // ---- Graphics setup
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );
    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );
    // ---- Gizmo
    cam.persp = false;
    gizmo.cam = &cam;
    //gizmo.bindPoints(W->ff.natoms, W->ff.apos      );
    //gizmo.bindEdges (W->ff.nbonds, W->ff.bond2atom );
    gizmo.bindPoints( natoms, apos      );
    gizmo.bindEdges ( nbonds, bond2atom );
    gizmo.pointSize = 0.5;
    //gizmo.iDebug    = 2;
    ruler.setStep( 1.5 * sqrt(3) );
    initWiggets();
}

void MolGUI::init(){
    if(verbosity>0)printf("MolGUI::init() \n");
    W->init( true );
    //MolGUI::bindMolecule( W->ff.natoms, W->ff.nbonds,W->ff.atypes,W->ff.bond2atom,Vec3d* fapos_,Quat4d* REQs_,Vec2i*  bond2atom_, Vec3d* pbcShifts_ );
    //MolGUI::bindMolecule( W->nbmol.natoms, W->ff.nbonds, W->nbmol.atypes, W->nbmol.apos, W->nbmol.fapos, W->nbmol.REQs,                         0,0, W->ff.bond2atom, W->ff.pbcShifts );
    MolGUI::bindMolecule( W->nbmol.natoms, W->ffl.nnode, W->ff.nbonds, W->nbmol.atypes, W->nbmol.apos, W->nbmol.fapos, W->nbmol.REQs, W->ffl.pipos, W->ffl.fpipos, W->ff.bond2atom, W->ff.pbcShifts );
    neighs    = W->ffl.neighs;
    neighCell = W->ffl.neighCell;
    initGUI();
    if(verbosity>0)printf("... MolGUI::init() DONE\n");
}


//=================================================
//                   DRAW()
//=================================================

void MolGUI::draw(){
    //printf( "MolGUI::draw() 1 \n" );
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);

    if( (ogl_isosurf==0) && W->bGridFF ){ renderGridFF(); }
    //if( ogl_esp==0 ){ renderESP(); }

    if(frameCount==1){ qCamera.pitch( M_PI );  qCamera0=qCamera; }

    //debug_scanSurfFF( 100, {0.,0.,z0_scan}, {0.0,3.0,z0_scan}, 10.0 );

    W->pick_hray = (Vec3d)cam.rot.c;
    W->pick_ray0 = ray0;

    if(bRunRelax){ W->MDloop(perFrame); }
    //if(bRunRelax){ W->relax( perFrame ); }

    // --- Mouse Interaction / Visualization
	ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
    Draw3D::drawPointCross( ray0, 0.1 );        // Mouse Cursor 
    //if(W->ipicked>=0) Draw3D::drawLine( W->ff.apos[W->ipicked], ray0); // Mouse Dragging Visualization
    if(W->ipicked>=0) Draw3D::drawLine( apos[W->ipicked], ray0); // Mouse Dragging Visualization
    Vec3d ray0_ = ray0;            ray0_.y=-ray0_.y;
    Vec3d ray0_start_=ray0_start;  ray0_start_.y=-ray0_start_.y;
    if(bDragging)Draw3D::drawTriclinicBoxT(cam.rot, (Vec3f)ray0_start_, (Vec3f)ray0_ );   // Mouse Selection Box

    if(ogl_MO){ 
        glPushMatrix();
        Vec3d c = W->builder.lvec.a*-0.5 + W->builder.lvec.b*-0.5 + W->builder.lvec.c*-0.5;
        glTranslatef( c.x, c.y, c.z );
            glColor3f(1.0,1.0,1.0); 
            glCallList(ogl_MO); 
        glPopMatrix();
    }

    //printf( "bViewSubstrate %i ogl_isosurf %i W->bGridFF %i \n", bViewSubstrate, ogl_isosurf, W->bGridFF );

    if( bViewSubstrate && W->bSurfAtoms ) Draw3D::atomsREQ( W->surf.natoms, W->surf.apos, W->surf.REQs, ogl_sph, 1., 0.1, 0., true, W->gridFF.shift0 );
    //if( bViewSubstrate && W->bSurfAtoms ) Draw3D::atomsREQ( W->surf.natoms, W->surf.apos, W->surf.REQs, ogl_sph, 1., 1., 0. );
    //if( bViewSubstrate                  ){ glColor3f(0.,0.,1.); Draw3D::drawTriclinicBoxT( W->gridFF.grid.cell, Vec3d{0.0, 0.0, 0.0}, Vec3d{1.0, 1.0, 1.0} ); }
    if( bViewSubstrate                  ){ glColor3f(0.,0.,1.); Draw3D::drawTriclinicBoxT( W->gridFF.grid.cell, Vec3d{-0.5, -0.5, 0.0}, Vec3d{0.5, 0.5, 1.0} ); }
    if( bViewSubstrate && ogl_isosurf   ) viewSubstrate( 3, 3, ogl_isosurf, W->gridFF.grid.cell.a, W->gridFF.grid.cell.b, W->gridFF.shift0 + W->gridFF.grid.pos0 );

    if( ogl_esp ){ glCallList(ogl_esp);  }

    //Draw3D::drawMatInPos( W->debug_rot, W->ff.apos[0] ); // DEBUG  

    //if(bDoQM)drawSystemQMMM();

    if(bDoMM){
        if(W->builder.bPBC){ 
            //Draw3D::drawPBC( (Vec3i){2,2,0}, W->builder.lvec, [&](Vec3i ixyz){drawSystem(ixyz);} ); 
            //printf( "draw() W->npbc=%i \n", W->npbc );
            Draw3D::drawShifts( W->npbc, W->pbc_shifts, 4, [&](Vec3i ixyz){drawSystem(ixyz);} ); 
            glColor3f(0.,0.5,0.5); Draw3D::drawTriclinicBoxT( W->builder.lvec, Vec3d{0.,0.,0.}, Vec3d{1.,1.,1.} );
        }else{ drawSystem(); }

        //drawSystem(); // DEBUG
        //Draw3D::drawNeighs( W->ff, -1.0 );    
        //Draw3D::drawVectorArray( W->ff.natoms, W->ff.apos, W->ff.fapos, 10000.0, 100.0 );
    }

    glColor3f(0.0f,0.5f,0.0f); showBonds();


    if(W->ipicked>-1){ 
        Vec3d p = W->ffl.apos[W->ipicked];
        Vec3d f = getForceSpringRay( W->ffl.apos[W->ipicked], W->pick_hray, W->pick_ray0, W->Kpick ); 
        glColor3d(1.f,0.f,0.f); Draw3D::drawVecInPos( f*-ForceViewScale, p );
    }

    if(bDebug_scanSurfFF){ // --- GridFF debug_scanSurfFF()
        Vec3d p0=W->gridFF.grid.pos0; 
        if(W->ipicked>-1){ p0.z = W->ffl.apos[W->ipicked].z; }{
            p0.z = W->ffl.apos[0].z;
        } 
        //printf( "p0.x %g \n", p0.z );
        int nx=W->gridFF.grid.n.x*3*5;
        int ny=W->gridFF.grid.n.y*3*5;
        Vec3d a = W->gridFF.grid.cell.a;
        Vec3d b = W->gridFF.grid.cell.b;
        //debug_scanSurfFF( nx, p0-a,       p0+a*2.       );
        //debug_scanSurfFF( ny, p0-b,       p0+b*2.       );
        debug_scanSurfFF( nx, p0-a+b*0.5, p0+a*2.+b*0.5 );
        debug_scanSurfFF( ny, p0-b+a*0.5, p0+b*2.+a*0.5 );
    }

    for(int i=0; i<W->selection.size(); i++){ 
        int ia = W->selection[i];
        glColor3f( 0.f,1.f,0.f ); Draw3D::drawSphereOctLines( 8, 0.3, W->nbmol.apos[ia] );     
    }

    // --- Drawing Population of geometies overlay
    if(frameCount>=1){ 
        if(frameCount==1){
            nsys=W->getMultiSystemPointers( M_neighs, M_neighCell, M_apos, nvec); 
            //for(int i=0;i<nsys; i++){ printMSystem(i, 4, natoms, nvec ); }
        }
        //printf( "nsys %i \n", nsys );
        if(nsys>0){ for(int isys=0; isys<nsys; isys++){ 
            //printf( "#=========== isys= %i \n", isys );
            Draw3D::neighs_multi(natoms,4,M_neighs,M_neighCell,M_apos, W->pbc_shifts, isys, nvec ); 
        } } 
    }
    

    //if(iangPicked>=0){
    //    glColor3f(0.,1.,0.);      Draw3D::angle( W->ff.ang2atom[iangPicked], W->ff.ang_cs0[iangPicked], W->ff.apos, fontTex3D );
    //}
    if(useGizmo){ gizmo.draw(); }
    if(bHexDrawing)drawingHex(5.0);
};

void MolGUI::printMSystem( int isys, int perAtom, int na, int nvec, bool bNg, bool bNgC, bool bPos ){
    int i0a = (isys * na  )*perAtom;
    int i0v = (isys * nvec); //*perAtom;
    printf( "### MolGUI::MSystem[%i] \n", isys );
    //if(blvec ){printf("lvec\n");  printMat(lvecs [isys]); }
    //if(bilvec){printf("ilvec\n"); printMat(ilvecs[isys]); }
    for(int i=0; i<na; i++){
        int ii=i*perAtom;
        printf("[%i]", i );
        if(bNg )printf("ng (%3i,%3i,%3i,%3i)",  M_neighs   [i0a+ii+0], M_neighs   [i0a+ii+1], M_neighs   [i0a+ii+2], M_neighs   [i0a+ii+3] );
        if(bNgC)printf("ngC(%2i,%2i,%2i,%2i)",  M_neighCell[i0a+ii+0], M_neighCell[i0a+ii+1], M_neighCell[i0a+ii+2], M_neighCell[i0a+ii+3] );
        if(bPos)printf("pa(%6.3f,%6.3f,%6.3f)", M_apos     [i0v+i].x , M_apos     [i0v+i].y , M_apos     [i0v+i].z );
        printf("\n" );
    }
}

void MolGUI::showBonds(  ){
    Vec3d* ps = W->ffl.apos;
    glBegin(GL_LINES);
    for(int i=0; i<bondsToShow.size(); i++ ){
        Vec2i b  = bondsToShow       [i];
        Vec3d d  = bondsToShow_shifts[i];
        Vec3d pi = ps[b.i];
        Vec3d pj = ps[b.j];
        printf( "b[%i,%i] pi(%7.3f,%7.3f,%7.3f) pj(%7.3f,%7.3f,%7.3f) d(%7.3f,%7.3f,%7.3f) \n", pi.x,pi.y,pi.z,    pj.x,pj.y,pj.z,  d.x,d.y,d.z );
        Draw3D::vertex(pi-d);  Draw3D::vertex(pj  );
        Draw3D::vertex(pi  );  Draw3D::vertex(pj+d);
    }
    glEnd();
}

void MolGUI::showAtomGrid( char* s, int ia, bool bDraw ){
    if( (ia>=W->ffl.natoms) ){ printf( "ERROR showAtomGrid(%i) out of atom range [0 .. %i] \n", ia, W->ffl.natoms ); }
    const GridShape& grid = W->gridFF.grid;

    Vec3d  pi  = W->ffl.apos[ia];
    Quat4f PLQ = W->ffl.PLQs[ia];

    Quat4f fp = W->gridFF.getForce( pi, Quat4f{ PLQ.x , 0.    , 0.    , 0. }, true );
    Quat4f fl = W->gridFF.getForce( pi, Quat4f{ 0.0   , PLQ.y , 0.    , 0. }, true );
    Quat4f fq = W->gridFF.getForce( pi, Quat4f{ 0.0   , 0.    , PLQ.z , 0. }, true );
    Quat4f fe = fp+fl+fq;

    s += sprintf(s, "GridFF(ia=%i) Etot %15.10f Epaul %15.10f; EvdW %15.10f Eel %15.10f p(%7.3f,%7.3f,%7.3f) \n", ia,   fe.e, fp.e, fl.e, fq.e,  pi.x,pi.y,pi.z  );
    //if( fabs(fe.e-fe_ref.e)>1e-8 ){ s += sprintf(s, "ERROR: getLJQH(%15.10f) Differs !!! \n", fe_ref.e ); }
    if( fe.e>0 ){ glColor3f(0.7f,0.f,0.f); }else{ glColor3f(0.f,0.f,1.f); }
    Draw::drawText( str, fontTex, fontSizeDef, {150,20} );
    glTranslatef( 0.0,fontSizeDef*2,0.0 );
    
}

Vec3d MolGUI::showNonBond( char* s, Vec2i b, bool bDraw ){
    int na = W->ffl.natoms ;
    if( (b.i>=na)||(b.j>=na) ){ printf( "ERROR showNonBond(%i,%i) out of atom range [0 .. %i] \n", b.i,b.j, W->ffl.natoms ); }
    Quat4d* REQs = W->ffl.REQs;
    Vec3d * ps   = W->ffl.apos;
    Mat3d& lvec  = W->ffl.lvec;
    Quat4d REQH  = _mixREQ( REQs[b.i], REQs[b.j] );
    Vec3d pi     = ps[b.i];
    Vec3d pj     = ps[b.j];
    Vec3d d      = pj-pi;
    // --- PBC
    
    Vec3i g      = W->ffl.invLvec.nearestCell( d );
    Vec3d shift  = lvec.a*g.x + lvec.b*g.y + lvec.c*g.z;
    d.add( shift );

    // --- Reference
    double R2damp = W->gridFF.Rdamp;
    Vec3d f;
    double Eref  = getLJQH( d, f, REQH, R2damp );
    // --- Decomposed
    double  r2  = d.norm2();
    // ---- Electrostatic
    double ir2_ = 1/( r2 + R2damp  );
    double Eel =  COULOMB_CONST*REQH.z*sqrt( ir2_ );
    double Fel =  Eel*ir2_ ;
    // --- LJ 
    double  ir2   = 1/r2;
    double  u2    = REQH.x*REQH.x*ir2;
    double  u6    = u2*u2*u2;
    double vdW    = u6*REQH.y;
    double EvdW   = -2.  *vdW;
    double Epaul  = u6   *vdW;
    double EH     = u6*u6* ((REQH.w<0) ? REQH.w : 0.0);  // H-bond correction
    double Etot   = Epaul + EvdW + EH + Eel;
    //F      +=  (12.*(u6-1.)*vdW + H*6.)*ir2;
    s += sprintf(s, "PLQH[%i,%i] r=%6.3f Etot %15.10f Epaul %15.10f; EvdW %15.10f EH %15.10f Eel %15.10f\n", b.i,b.j, sqrt(r2), Etot, Epaul, EvdW, EH, Eel );
    if( fabs(Etot-Eref)>1e-8 ){ s += sprintf(s, "ERROR: getLJQH(%15.10f) Differs !!! \n", Eref ); }
    if( Etot>0 ){ glColor3f(0.7f,0.f,0.f); }else{ glColor3f(0.f,0.f,1.f); }
    Draw::drawText( str, fontTex, fontSizeDef, {150,20} );
    glTranslatef( 0.0,fontSizeDef*2,0.0 );
    return shift;
}

void MolGUI::drawHUD(){
    glDisable ( GL_LIGHTING );
    gui.draw();

    if(W->bCheckInvariants){
        glTranslatef( 10.0,HEIGHT-20.0,0.0 );
        glColor3f(0.5,0.0,0.3);
        //char* s=str;
        //printf( "(%i|%i,%i,%i) cog(%g,%g,%g) vcog(%g,%g,%g) fcog(%g,%g,%g) torq (%g,%g,%g)\n", ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0,  cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, fcog.x,fcog.y,fcog.z, tq.x,tq.y,tq.z );
        //printf( "neval Ang %i nevalPiSigma %i PiPiT %i PiPiI %i v_av %g \n", ff.nevalAngles, ff.nevalPiSigma, ff.nevalPiPiT, ff.nevalPiPiI, v_av );
        // s += sprintf(s, "iSystemCur %i\n",  W->iSystemCur );
        // if(W->bMMFF) s += sprintf(s, "eval:Ang,ps,ppT,ppI(%i|%i,%i,%i)\n",  W->ff.nevalAngles>0, W->ff.nevalPiSigma>0, W->ff.nevalPiPiT>0, W->ff.nevalPiPiI>0 );
        // s += sprintf(s, "cog (%g,%g,%g)\n", W->cog .x,W->cog .y,W->cog .z);
        // s += sprintf(s, "vcog(%15.5e,%15.5e,%15.5e)\n", W->vcog.x,W->vcog.y,W->vcog.z);
        // s += sprintf(s, "fcog(%15.5e,%15.5e,%15.5e)\n", W->fcog.x,W->fcog.y,W->fcog.z);
        // s += sprintf(s, "torq(%15.5e,%15.5e,%15.5e)\n", W->tqcog.x,W->tqcog.y,W->tqcog.z);
        W->getStatusString( str, nmaxstr );
        Draw::drawText( str, fontTex, fontSizeDef, {100,20} );
    }

    if(bWriteOptimizerState){
        glTranslatef( 0.0,fontSizeDef*-5*2,0.0 );
        glColor3f(0.0,0.5,0.0);
        char* s=str;
        //dt 0.132482 damp 3.12175e-17 n+ 164 | cfv 0.501563 |f| 3.58409e-10 |v| 6.23391e-09
        double v=sqrt(W->opt.vv);
        double f=sqrt(W->opt.ff);
        s += sprintf(s,"dt %7.5f damp %7.5f n+ %4i | cfv %7.5f |f| %12.5e |v| %12.5e \n", W->opt.dt, W->opt.damping, W->opt.lastNeg, W->opt.vf/(v*f), f, v );
        Draw::drawText( str, fontTex, fontSizeDef, {100,20} );
        glTranslatef( 0.0,fontSizeDef*-5*2,0.0 );
        Draw::drawText( W->info_str(str), fontTex, fontSizeDef, {100,20} );

    }

    /*
    glTranslatef( 0.0,fontSizeDef*-2*2,0.0 );
    if( !bondsToShow_shifts ){
        bondsToShow.push_back( {18,36} );
        bondsToShow.push_back( {35,7}  );
        _alloc( bondsToShow_shifts, bondsToShow.size() );
        // ---- Setup H-bond corrections
        Quat4d* REQs = W->ffl.REQs;
        REQs[18].w = -0.7 * REQs[18].y;
        REQs[36].w = +0.7 * REQs[36].y;
        REQs[35].w = +0.7 * REQs[35].y;
        REQs[7 ].w = -0.7 * REQs[7 ].y;
    }
    DEBUG
    for(int i=0; i<bondsToShow.size(); i++ ){
        bondsToShow_shifts[i] =  MolGUI::showNonBond( str, bondsToShow[i] );
    }
    */

    mouse_pix = ((Vec2f){ 2*mouseX/float(HEIGHT) - ASPECT_RATIO,
                          2*mouseY/float(HEIGHT) - 1      });// *(1/zoom);
}

void MolGUI::drawingHex(double z0){
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

void MolGUI::selectRect( const Vec3d& p0, const Vec3d& p1 ){ W->selectRect( p0, p1, (Mat3d)cam.rot ); }

void  MolGUI::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
    int ib = pickBondCenter( W->ff.nbonds, W->ff.bond2atom, W->ff.apos, ro, rd, 0.5 );
    printf( "picked bond  %i \n", ib );
    W->selection.reserve(W->ff.natoms);
    W->splitAtBond( ib, &(W->selection[0]) );
}

void MolGUI::renderGridFF( double isoVal, int isoSurfRenderType, double colorSclae ){
    if(verbosity>0) printf( "MolGUI::renderGridFF()\n" );
    //int iatom = 11;
    testREQ = Quat4d{ 1.487, sqrt(0.0006808), 0., 0.}; // H
    testPLQ = REQ2PLQ( testREQ, W->gridFF.alphaMorse );
    Quat4f * FFtot = new Quat4f[ W->gridFF.grid.getNtot() ];
    W->gridFF.evalCombindGridFF ( testREQ, FFtot );
    //W->gridFF.grid.saveXSF( "E_renderGridFF.xsf",  (float*)FFtot, 4, 3, W->gridFF.natoms, W->gridFF.atypes, W->gridFF.apos );
    //if(idebug>0) W->gridFF.grid.saveXSF( "FFtot_z.xsf",  (float*)FFtot, 4, 2, W->gridFF.natoms, W->gridFF.atypes, W->gridFF.apos );
    ogl_isosurf = glGenLists(1);
    glNewList(ogl_isosurf, GL_COMPILE);
    glShadeModel( GL_SMOOTH );
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    //bool sign=false;
    bool sign=true;
    int nvert = renderSubstrate_( W->gridFF.grid, FFtot, W->gridFF.FFelec, +isoVal, sign, colorSclae );   printf("DEBUG renderGridFF() renderSubstrate() -> nvert= %i ", nvert );
    // ---- This seems still not work properly
    //int ntris=0;
    //glColor3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( W->gridFF.grid,  isoVal, (double*)FFtot, isoSurfRenderType,  3,2 );
    //glColor3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( W->gridFF.grid, -isoVal, (double*)FFtot, isoSurfRenderType,  3,2 );
    //glColor3f(0.,0.,1.); Draw3D::drawTriclinicBox( W->gridFF.grid.cell, Vec3d{0.0, 0.0, 0.0}, Vec3d{1.0, 1.0, 1.0} );
    //Draw3D::drawAxis(1.0);
    glEndList();
    delete [] FFtot;
    if(verbosity>0) printf( "... MolGUI::renderGridFF() DONE\n" );
}

void MolGUI::renderESP( Quat4d REQ){
    printf( "DEBUG MolGUI::renderESP() %li \n", ogl_esp ); //exit(0);
    ogl_esp = glGenLists(1);
    glNewList(ogl_esp, GL_COMPILE);
    glShadeModel( GL_SMOOTH );
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    const NBFF& nbmol = W->nbmol;
    int nvert = Draw3D::drawESP( nbmol.natoms, nbmol.apos, nbmol.REQs, REQ );
    glEndList();
};

void MolGUI::drawPi0s( float sc=1.0 ){
    const MMFFsp3& ff = W->ff;
    glBegin(GL_LINES);
    for(int ia=0; ia<ff.nnode; ia++){
        int* ngs = ff.neighs + ia*ff.nneigh_max;
        for(int j=0; j<ff.nneigh_max; j++){
            int ing = ngs[j]; 
            if(ing<0){
                int ipi = -ing-2;
                Vec3f p=(Vec3f)ff.apos[ia];      glVertex3f(p.x,p.y,p.z);
                p.add_mul( ff.pi0s[ipi].f, sc);  glVertex3f(p.x,p.y,p.z);
                //printf("drawPi0s[%i,%i|%i] (%g,%g,%g)\n", ia, j, ipi, ff.pi0s[ipi].f.z, ff.pi0s[ipi].f.y, ff.pi0s[ipi].f.z );
            }
        }
    }
    glEnd();
}

/*
void MolGUI::flipPis( Vec3d ax ){
    for(int i=0; i<nnode; i++){
        double c = pipos[i].dot(ax);
        if( c<0 ){ pipos[i].mul(-1); } 
    }
}
*/

void MolGUI::bindMolecule( int natoms_, int nnode_, int nbonds_, int* atypes_,Vec3d* apos_,Vec3d* fapos_,Quat4d* REQs_, Vec3d* pipos_, Vec3d* fpipos_, Vec2i* bond2atom_, Vec3d* pbcShifts_ ){
    natoms=natoms_; nnode=nnode_; nbonds=nbonds_;
    if(atypes_)atypes=atypes_;
    if(apos_  )apos  =apos_;
    if(fapos_ )fapos =fapos_;
    if(REQs_  )REQs  =REQs_;
    if(pipos_ )pipos =pipos_;
    if(fpipos_)fpipos=fpipos_;
    if(bond2atom_)bond2atom=bond2atom_;
    if(pbcShifts_)pbcShifts=pbcShifts_;
}

void MolGUI::drawSystem( Vec3i ixyz ){
    //float textSize=0.007;
    float textSize=1.0;
    glEnable(GL_DEPTH_TEST);
    bool bOrig = (ixyz.x==0)&&(ixyz.y==0)&&(ixyz.z==0);
    //printf( "bOrig %i ixyz(%i,%i,%i)\n", bOrig, ixyz.x,ixyz.y,ixyz.z );
    //printf( "DEBUG MolGUI::drawSystem() bViewMolCharges %i W->nbmol.REQs %li\n", bViewMolCharges, W->nbmol.REQs );
    //printf("DEBUG MolGUI::drawSystem()  bOrig %i W->bMMFF %i mm_bAtoms %i bViewAtomSpheres %i bViewAtomForces %i bViewMolCharges %i \n", bOrig, W->bMMFF, mm_bAtoms, bViewAtomSpheres, bViewAtomForces, bViewMolCharges  );
    if( bond2atom ){
        //if(W->builder.bPBC){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC    ( nbonds, bond2atom, apos, &W->builder.bondPBC[0], W->builder.lvec ); } 
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::bonds       ( nbonds,bond2atom,apos );  // DEBUG
        /*
        if(W->builder.bPBC){ glColor3f(0.0f,0.0f,0.0f); if(pbcShifts)Draw3D::bondsPBC          ( nbonds, bond2atom, apos,  pbcShifts                          );  
                             glColor3f(0.0f,0.0f,1.0f); if(pbcShifts)Draw3D::pbcBondNeighLabels( nbonds, bond2atom, apos,  pbcShifts, fontTex3D,        textSize );
        }else              { glColor3f(0.0f,0.0f,0.0f); Draw3D::bonds             ( nbonds, bond2atom, apos                                            );                                          
                             glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsLengths      ( nbonds, bond2atom, apos, fontTex );                                
        }
        */
        if(bOrig){
            if(bViewPis &&  fpipos ){ glColor3f(0.0f,1.0f,1.0f); Draw3D::drawVectorArray( nnode, apos, pipos, 1.0, 100.0 );          }
            //if(mm_bAtoms         ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels       ( natoms, apos, fontTex3D               ); }                    
            if( bViewBondLabels    ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondLabels( nbonds, bond2atom, apos, fontTex3D,  textSize  ); }
            //void bondLabels( int n, const Vec2i* b2a, const Vec3d* apos, int fontTex, float sz=0.02 );
            //Draw3D::atoms          ( natoms, apos, atypes, W->params, ogl_sph, 1.0, mm_Rsc, mm_Rsub );
            //Draw3D::drawVectorArray( natoms, apos, fapos, 100.0, 10000.0 );   
        }

    }
    if( neighs ){  glColor3f(0.0f,0.0f,0.0f);   Draw3D::neighs(  natoms, 4, (int*)neighs, (int*)neighCell, apos, W->pbc_shifts );   }
    //W->nbmol.print();
    if(bViewAtomSpheres&&mm_bAtoms                  ){                            Draw3D::atoms            ( natoms, apos, atypes, W->params, ogl_sph, 1.0, mm_Rsc, mm_Rsub ); }
    if(bOrig){
        //if(bViewAtomP0s     &&  fapos           ){ glColor3f(0.0f,1.0f,1.0f); Draw3D::drawVectorArray  ( natoms, apos, fapos, ForceViewScale, 10000.0 );  }
        if(bViewAtomForces    &&  fapos           ){ glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVectorArray  ( natoms, apos, fapos, ForceViewScale, 10000.0 );  }
        if(mm_bAtoms&&bViewAtomLabels             ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels       ( natoms, apos,                                    fontTex3D, textSize );  }
        if(mm_bAtoms&&bViewAtomTypes              ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomTypes        ( natoms, apos, atypes, &(params_glob->atypes[0]), fontTex3D, textSize );  }
        if(bViewMolCharges && (W->nbmol.REQs!=0)  ){ glColor3f(0.0,0.0,0.0);    Draw3D::atomPropertyLabel( natoms,  (double*)REQs,  apos, 4, 2,             fontTex3D, textSize ); }
        //if(W->ff.pi0s                           ){ glColor3f(0.0f,1.0f,1.0f); drawPi0s(1.0); }
    }

}

void MolGUI::saveScreenshot( int i, const char* fname ){
    char str[64];
    sprintf( str, fname, i );               
    printf( "save to %s \n", str );
    unsigned int *screenPixels = new unsigned int[WIDTH*HEIGHT*4];  
    glFlush();                                                      
    glFinish();                                                     
    //glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_INT, screenPixels);
    glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, screenPixels);  
    //SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff );   
    SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000 );  
    SDL_SaveBMP(bitmap, str);   
    SDL_FreeSurface(bitmap);
    delete[] screenPixels;
}

void MolGUI::debug_scanSurfFF( int n, Vec3d p0, Vec3d p1, Quat4d REQ, double sc ){
    if(isnan(sc)){ sc=ForceViewScale; }
    //printf("========= MolGUI::scanFF\n");
    //sc=10.0;
    //p0=Vec3d{0.0,0.0,z0_scan}; p1=Vec3d{5.0,5.0,z0_scan}; // Horizontal scan
    //p0=Vec3d{0.0,0.0,z0_scan}; p1=Vec3d{0.0,5.0,z0_scan}; // Horizontal scan
    //p0=Vec3d{0.0,z0_scan,0.0}; p1=Vec3d{0.0,z0_scan,10.0,}; // Vertical scan
    Vec3d dp=p1-p0; dp.mul(1./n);
    //Quat4f PLQ = REQ2PLQ( REQ, W->gridFF.alphaMorse );
    //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    debug_ps  .resize(n);
    debug_fs  .resize(n);
    debug_REQs.resize(n);
    for(int i=0; i<n; i++){ 
        debug_REQs[i]  =(Quat4f)REQ;
        debug_ps  [i].f=(Vec3f )(p0+dp*i);
    }
    W->scanSurfFF( n, &debug_ps[0], &debug_REQs[0], &debug_fs[0] );
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        //Vec3d  p  = p0 + dp*i;
        //Quat4f fe = Quat4fZero; W->gridFF.addForce_surf( p, PLQ, fe );
        //Quat4f fe = W->gridFF.getForce( p, PLQ, true );
        //printf( "debug_scanSurfFF[%i] p(%6.3f,%6.3f,%6.3f) fe(%g,%g,%g|%g)\n", i,  p.x,p.y,p.z, fe.x,fe.y,fe.z,fe.w );
        Vec3f f = debug_fs[i].f;
        Vec3f p = debug_ps[i].f;

        Draw3D::vertex( p ); Draw3D::vertex(p + (Vec3f)dp   );
        Draw3D::vertex( p ); Draw3D::vertex(p + f*sc        );

        //double E=0;
        //if   (W->bGridFF){ E+= W->gridFF.eval        (1, W->nbmol.apos, W->nbmol.PLQs,  W->nbmol.fapos );                        }
        //else           { E+= W->nbmol .evalMorse   (W->surf, false,                        W->gridFF.alphaMorse, W->gridFF.Rdamp ); }
        //else             { E+= W->nbmol .evalMorsePBC(W->surf, W->gridFF.grid.cell, W->nPBC, W->gridFF.alphaMorse, W->gridFF.Rdamp ); }
        //Draw3D::vertex( W->nbmol.apos[0] ); Draw3D::vertex( W->nbmol.apos[0]+W->nbmol.fapos[0]*sc );       // Force Vectro
        //Draw3D::vertex( W->nbmol.apos[0] ); Draw3D::vertex( W->nbmol.apos[0]+(Vec3d{0.0,0.0,E})*sc  );  // Energy -> z
        //Draw3D::vertex( W->nbmol.apos[0] ); Draw3D::vertex( W->nbmol.apos[0]+(Vec3d{0.0,E,0.0})*sc  );  // Energy -> x
    }
    glEnd();
}

void MolGUI::mouse_default( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT: { ray0_start = ray0;  bDragging = true; }break;
                case SDL_BUTTON_RIGHT:{ }break;
            } break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if( ray0.dist2(ray0_start)<0.1 ){
                        int ipick = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, W->nbmol.natoms, W->nbmol.apos );
                        if( ipick == W->ipicked ){ W->ipicked=-1; }else{ W->ipicked = ipick; }; 
                        W->selection.clear();
                        if(W->ipicked>=0){ 
                            W->selection.push_back(W->ipicked); 
                            Qpanel->value = W->nbmol.REQs[W->ipicked].z;
                            Qpanel->redraw=true;
                        };
                        printf( "picked atom %i \n", W->ipicked );
                    }else{
                        selectRect( ray0_start, ray0 );
                    }
                    bDragging=false;
                    break;
                case SDL_BUTTON_RIGHT:{ W->ipicked=-1; } break;
            }break;
    }
}

void MolGUI::eventMode_edit( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN : { if(gui.bKeyEvents)switch( event.key.keysym.sym ){
                case SDLK_t:{
                    affineTransform( W->ff.natoms, W->ff.apos, W->ff.apos, W->builder.lvec, W->new_lvec );
                    W->builder.updatePBC( W->ff.pbcShifts, &(W->new_lvec) );
                    _swap( W->builder.lvec, W->new_lvec );
                }break;
                case SDLK_i:
                    //selectShorterSegment( (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y + cam.rot.c*-1000.0), (Vec3d)cam.rot.c );
                    selectShorterSegment( ray0, (Vec3d)cam.rot.c );
                    //selection.erase();
                    //for(int i:builder.selection){ selection.insert(i); };
                    break;
            }}break;
        case SDL_WINDOWEVENT:{switch (event.window.event) {case SDL_WINDOWEVENT_CLOSE:{ quit(); }break;} } break;
    }
}

void MolGUI::eventMode_scan( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN : { if(gui.bKeyEvents)switch( event.key.keysym.sym ){
            //case SDLK_KP_1: picked_lvec = &W->builder.lvec.a; break;
            //case SDLK_KP_2: picked_lvec = &W->builder.lvec.b; break;
            //case SDLK_KP_3: picked_lvec = &W->builder.lvec.c; break;
            //case SDLK_KP_7: picked_lvec->x+=xstep; break;
            //case SDLK_KP_4: picked_lvec->x-=xstep; break;
            //case SDLK_KP_8: picked_lvec->y+=xstep; break;
            //case SDLK_KP_5: picked_lvec->y-=xstep; break;
            //case SDLK_KP_9: picked_lvec->z+=xstep; break;
            //case SDLK_KP_6: picked_lvec->z-=xstep; break;

            //case SDLK_COMMA:  W->change_lvec( W->ffl.lvec+dlvec    ); break;
            //case SDLK_PERIOD: W->change_lvec( W->ffl.lvec+dlvec*-1 ); break;

            case SDLK_COMMA:  W->add_to_lvec( dlvec    ); break;
            case SDLK_PERIOD: W->add_to_lvec( dlvec*-1 ); break;

            case SDLK_a: apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
            case SDLK_d: apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

            case SDLK_w: apos[1].mul( 1.1 ); break;
            case SDLK_s: apos[1].mul( 0.9 ); break;

            case SDLK_KP_PLUS:  z0_scan = z0_scan+0.1; break;
            case SDLK_KP_MINUS: z0_scan = z0_scan-0.1; break;

            case SDLK_SPACE: bRunRelax=!bRunRelax; break;

        }  }break;
        case SDL_MOUSEBUTTONDOWN:
        case SDL_MOUSEBUTTONUP:
            mouse_default( event );
            break;
        case SDL_WINDOWEVENT:{switch (event.window.event) {case SDL_WINDOWEVENT_CLOSE:{ quit(); }break;} } break;
    }
}

void MolGUI::eventMode_default( const SDL_Event& event ){
    if(useGizmo)gizmo.onEvent( mouse_pix, event );
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN : if(gui.bKeyEvents) switch( event.key.keysym.sym ){
                case SDLK_KP_0: qCamera = qCamera0; break;

                //case SDLK_COMMA:  which_MO--; printf("which_MO %i \n", which_MO ); break;
                //case SDLK_PERIOD: which_MO++; printf("which_MO %i \n", which_MO ); break;

                case SDLK_COMMA:  W->add_to_lvec( dlvec    ); break;
                case SDLK_PERIOD: W->add_to_lvec( dlvec*-1 ); break;

                //case SDLK_LESS:    which_MO--; printf("which_MO %i \n"); break;
                //case SDLK_GREATER: which_MO++; printf("which_MO %i \n"); break;

                //case SDLK_m: renderOrbital( which_MO ); break;
                //case SDLK_r: renderDensity(          ); break;
                case SDLK_s: W->saveXYZ( "out.xyz", "#comment", false ); break;
                case SDLK_p: saveScreenshot( frameCount ); break;
                
                //case SDLK_o: W->optimizeLattice_1d( 10,40, Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  } ); break;
                case SDLK_o:{
                        long T0 = getCPUticks();
                        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                        W->optimizeLattice_1d( 20,20, Mat3d{   0.0,0.5,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  } ); 
                        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                        double time_s     = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() * 1e-9;
                        double time_GTick = (getCPUticks()-T0)*1e-9;
                        printf( "Time{W->optimizeLattice_1d(20,20)} %g[s] %g[GTick]  %g[GTick/s] \n", time_s,time_GTick, time_GTick/time_s  );
                    }break;
                case SDLK_u:{ 
                        long T0 = getCPUticks();
                        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                        W->upload_pop        ( "population.xyz" ); 
                        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                        double time_s     = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() * 1e-9;
                        double time_GTick = (getCPUticks()-T0)*1e-9;
                        printf( "Time{W->upload_pop(population.xyz)} %g[s] %g[GTick]  %g[GTick/s] \n", time_s,time_GTick, time_GTick/time_s  );
                    }break;
                //case SDLK_u: W->upload_pop        ( "population.xyz" ); break;
                //case SDLK_o: W->optimizeLattice_1d( 0,2, Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  } ); break;
                //case SDLK_LEFTBRACKET:  {iSystemCur++; int nsys=W->gopt.population.size(); if(iSystemCur>=nsys)iSystemCur=0;  W->gopt.setGeom( iSystemCur ); } break;
                //case SDLK_RIGHTBRACKET: {iSystemCur--; int nsys=W->gopt.population.size(); if(iSystemCur<0)iSystemCur=nsys-1; W->gopt.setGeom( iSystemCur ); } break;

                case SDLK_LEFTBRACKET:  W->prevSystemReplica(); break;
                case SDLK_RIGHTBRACKET: W->nextSystemReplica(); break;

                //case SDLK_g: useGizmo=!useGizmo; break;
                //case SDLK_g: W->bGridFF=!W->bGridFF; break;
                //case SDLK_g: W->swith_gridFF(); break;
                //case SDLK_c: W->autoCharges(); break;

                case SDLK_g: W->bGridFF=!W->bGridFF; break;
                case SDLK_c: W->bOcl=!W->bOcl;       break;
                case SDLK_m: W->swith_method();      break;
                case SDLK_h: W->ff4.bAngleCosHalf = W->ffl.bAngleCosHalf = !W->ffl.bAngleCosHalf; break;
                case SDLK_k: bDebug_scanSurfFF ^=1; break;
                //case SDLK_q: W->autoCharges(); break;
                case SDLK_a: bViewAtomSpheres ^= 1; break;
                case SDLK_l: bViewAtomLabels  ^= 1; break;
                case SDLK_t: bViewAtomTypes   ^= 1; break;
                case SDLK_b: bViewBondLabels  ^= 1; break;  
                case SDLK_q: bViewMolCharges  ^= 1; break;
                case SDLK_f: bViewAtomForces  ^= 1; break;
                case SDLK_w: bViewSubstrate   ^= 1; break;
                //case SDLK_LEFTBRACKET: rotate( W->selection.size(), &W->selection[0], W->ff.apos, rotation_center, rotation_axis, +rotation_step ); break;
                //case SDLK_RIGHTBRACKET: rotate( W->selection.size(), &W->selection[0], W->ff.apos, rotation_center, rotation_axis, -rotation_step );  break;
                case SDLK_SPACE: bRunRelax=!bRunRelax;  if(bRunRelax)W->setConstrains();  break;
                // case SDLK_d: {
                //     printf( "DEBUG Camera Matrix\n");
                //     printf( "DEBUG qCamera(%g,%g,%g,%g) \n", qCamera.x,qCamera.y,qCamera.z,qCamera.w );
                //     qCamera.toMatrix(cam.rot);
                //     printf( "DEBUG cam aspect %g zoom %g \n", cam.aspect, cam.zoom);
                //     printMat((Mat3d)cam.rot);
                // } break;
                //case SDLK_g: iangPicked=(iangPicked+1)%ff.nang;
                //    printf( "ang[%i] cs(%g,%g) k %g (%i,%i,%i)\n", iangPicked, ff.ang_cs0[iangPicked].x, ff.ang_cs0[iangPicked].y, ff.ang_k[iangPicked],
                //        ff.ang2atom[iangPicked].a,ff.ang2atom[iangPicked].b,ff.ang2atom[iangPicked].c );
                //    break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
        case SDL_MOUSEBUTTONUP: mouse_default( event ); break;
        case SDL_WINDOWEVENT:{switch (event.window.event) {case SDL_WINDOWEVENT_CLOSE:{ quit(); }break;} } break;
    } // switch( event.type ){
}

void MolGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    gui.onEvent( mouseX, mouseY, event );
    switch (gui_mode){
        case  Gui_Mode::edit: eventMode_edit(event); break;
        case  Gui_Mode::scan: eventMode_scan(event); break;
        case  Gui_Mode::base: 
        default:              eventMode_default(event); break;   
    }
    //AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void MolGUI::keyStateHandling( const Uint8 *keys ){
    double dstep=0.025;
    switch (gui_mode){
        case Gui_Mode::edit: 
        case Gui_Mode::scan: 
        case Gui_Mode::base:
        default:{
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
            if( keys[ SDL_SCANCODE_KP_4 ] ){ W->nbmol.shift( {-0.1,0.,0.} ); }
            if( keys[ SDL_SCANCODE_KP_6 ] ){ W->nbmol.shift( {+0.1,0.,0.} ); }
            if( keys[ SDL_SCANCODE_KP_8 ] ){ W->nbmol.shift( {0.,+0.1,0.} ); }
            if( keys[ SDL_SCANCODE_KP_2 ] ){ W->nbmol.shift( {0.,-0.1,0.} ); }
            if( keys[ SDL_SCANCODE_KP_7 ] ){ W->nbmol.shift( {0.,0.,+0.1} ); }
            if( keys[ SDL_SCANCODE_KP_9 ] ){ W->nbmol.shift( {0.,0.,-0.1} ); }
            if( keys[ SDL_SCANCODE_LEFT  ] ){ cam.pos.add_mul( cam.rot.a, -cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_RIGHT ] ){ cam.pos.add_mul( cam.rot.a,  cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_UP    ] ){ cam.pos.add_mul( cam.rot.b,  cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_DOWN  ] ){ cam.pos.add_mul( cam.rot.b, -cameraMoveSpeed ); }
            //AppSDL2OGL_3D::keyStateHandling( keys );
        } break;   
    }
};
