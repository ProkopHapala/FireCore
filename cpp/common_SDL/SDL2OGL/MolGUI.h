
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
    bool   bViewBondLabels  = false;
    bool   bViewAtomSpheres = true;
    bool   bViewAtomForces  = false;
    bool   bViewSubstrate   = true;
    bool   isoSurfRenderType = 1;
    Vec3d testREQ,testPLQ;

    // ---- Graphics objects
    int  fontTex,fontTex3D;


    int  ogl_esp=0;
    int  ogl_sph=0;
    int  ogl_mol=0;
    int  ogl_isosurf=0;
    int  ogl_MO = 0;

    char str[2048];

    // --------- Functions 

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

	MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_ );
    void init();

	//int  loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int  loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );
    void tryLoadGridFF();
    //void makeGridFF   (bool recalcFF=false, bool bRenderGridFF=true);
    void renderGridFF( double isoVal=0.01, int isoSurfRenderType=0, double colorSclae = 30.0 );
    void renderESP( Vec3d REQ=Vec3d{1.487,0.0006808,1.0} );

	void drawSystem( Vec3i ixyz=Vec3iZero );
    void drawPi0s( float sc );
    //void drawSystemQMMM();
    //void renderOrbital(int i, double iso=0.1);
    //void renderDensity(       double iso=0.1);
	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );
    void debug_scanSurfFF( int n, Vec3d p0, Vec3d p1, double sc );

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
    SDL_Delay(100);
    tick2second =  0.1/(getCPUticks()-T0);
    printf( ">>tick2second = %g ( %g GHz)\n", tick2second, 1.0e-9/tick2second );

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
    gizmo.bindPoints(W->ff.natoms, W->ff.apos      );
    gizmo.bindEdges (W->ff.nbonds, W->ff.bond2atom );
    gizmo.pointSize = 0.5;
    //gizmo.iDebug    = 2;
    ruler.setStep( 1.5 * sqrt(3) );
    initWiggets();
}

void MolGUI::init(){
    if(verbosity>0)printf("MolGUI::init() \n");
    W->init( true );
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
    if( ogl_esp==0 ){ renderESP(); }

    if(frameCount==1){ qCamera.pitch( M_PI );  qCamera0=qCamera; }

    //debug_scanSurfFF( 100, {0.,0.,z0_scan}, {0.0,3.0,z0_scan}, 10.0 );

    W->pick_hray = (Vec3d)cam.rot.c;
    W->pick_ray0 = ray0;

    if(bRunRelax){ W->MDloop(perFrame); }
    //if(bRunRelax){ W->relax( perFrame ); }

    // --- Mouse Interaction / Visualization
	ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
    Draw3D::drawPointCross( ray0, 0.1 );        // Mouse Cursor 
    if(W->ipicked>=0) Draw3D::drawLine( W->ff.apos[W->ipicked], ray0); // Mouse Dragging Visualization
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

    if( bViewSubstrate && W->bSurfAtoms ) Draw3D::atomsREQ( W->surf.n, W->surf.ps, W->surf.REQs, ogl_sph, 1., 0.1, 0., true );
    //if( bViewSubstrate && W->bSurfAtoms ) Draw3D::atomsREQ( W->surf.n, W->surf.ps, W->surf.REQs, ogl_sph, 1., 1., 0. );
    //if( bViewSubstrate                  ){ glColor3f(0.,0.,1.); Draw3D::drawTriclinicBoxT( W->gridFF.grid.cell, (Vec3d){0.0, 0.0, 0.0}, (Vec3d){1.0, 1.0, 1.0} ); }
    if( bViewSubstrate                  ){ glColor3f(0.,0.,1.); Draw3D::drawTriclinicBoxT( W->gridFF.grid.cell, (Vec3d){-0.5, -0.5, 0.0}, (Vec3d){0.5, 0.5, 1.0} ); }
    if( bViewSubstrate && ogl_isosurf   ) viewSubstrate( 3, 3, ogl_isosurf, W->gridFF.grid.cell.a, W->gridFF.grid.cell.b, W->gridFF.shift + W->gridFF.grid.pos0 );

    if( ogl_esp ){ glCallList(ogl_esp);  }

    Draw3D::drawMatInPos( W->debug_rot, W->ff.apos[4] ); // DEBUG  

    //if(bDoQM)drawSystemQMMM();

    if(bDoMM){
        //if(W->builder.bPBC){ Draw3D::drawPBC( (Vec3i){2,2,0}, W->builder.lvec, [&](Vec3d ixyz){drawSystem(ixyz);} ); } 
        //else               { drawSystem(); }

        //if( W->builder.bPBC ){ glColor3f(0.,0.5,0.5); Draw3D::drawTriclinicBox( W->builder.lvec, (Vec3d){-0.5, -0.5, 0.0}, (Vec3d){0.5, 0.5, 1.0} ); }
        if( W->builder.bPBC ){ glColor3f(0.,0.5,0.5); Draw3D::drawTriclinicBoxT( W->builder.lvec, (Vec3d){ 0.0, 0.0, 0.0}, (Vec3d){1.0, 1.0, 1.0} ); }

        drawSystem(); // DEBUG
        //Draw3D::drawNeighs( W->ff, -1.0 );    
        //Draw3D::drawVectorArray( W->ff.natoms, W->ff.apos, W->ff.fapos, 10000.0, 100.0 );
    }

    for(int i=0; i<W->selection.size(); i++){ 
        int ia = W->selection[i];
        glColor3f( 0.f,1.f,0.f ); Draw3D::drawSphereOctLines( 8, 0.3, W->nbmol.ps[ia] );     
    }
    //if(iangPicked>=0){
    //    glColor3f(0.,1.,0.);      Draw3D::angle( W->ff.ang2atom[iangPicked], W->ff.ang_cs0[iangPicked], W->ff.apos, fontTex3D );
    //}
    if(useGizmo){ gizmo.draw(); }
    if(bHexDrawing)drawingHex(5.0);
};

void MolGUI::drawHUD(){

    glDisable ( GL_LIGHTING );
    gui.draw();

    if(W->bCheckInvariants){
        glTranslatef( 10.0,HEIGHT-20.0,0.0 );
        glColor3f(0.5,0.0,0.3);
        char* s=str;
        //printf( "(%i|%i,%i,%i) cog(%g,%g,%g) vcog(%g,%g,%g) fcog(%g,%g,%g) torq (%g,%g,%g)\n", ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0,  cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, fcog.x,fcog.y,fcog.z, tq.x,tq.y,tq.z );
        //printf( "neval Ang %i nevalPiSigma %i PiPiT %i PiPiI %i v_av %g \n", ff.nevalAngles, ff.nevalPiSigma, ff.nevalPiPiT, ff.nevalPiPiI, v_av );
        if(W->bMMFF) s += sprintf(s, "eval:Ang,ps,ppT,ppI(%i|%i,%i,%i)\n",  W->ff.nevalAngles>0, W->ff.nevalPiSigma>0, W->ff.nevalPiPiT>0, W->ff.nevalPiPiI>0 );
        s += sprintf(s, "cog (%g,%g,%g)\n", W->cog .x,W->cog .y,W->cog .z);
        s += sprintf(s, "vcog(%15.5e,%15.5e,%15.5e)\n", W->vcog.x,W->vcog.y,W->vcog.z);
        s += sprintf(s, "fcog(%15.5e,%15.5e,%15.5e)\n", W->fcog.x,W->fcog.y,W->fcog.z);
        s += sprintf(s, "torq(%15.5e,%15.5e,%15.5e)\n", W->tqcog.x,W->tqcog.y,W->tqcog.z);
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
                    Draw3D::vertex((Vec3f){p.x,p.y,z0}); p2=p+ruler.lvecs[0]*sc; Draw3D::vertex((Vec3f){p2.x,p2.y,z0});
                    Draw3D::vertex((Vec3f){p.x,p.y,z0}); p2=p+ruler.lvecs[1]*sc; Draw3D::vertex((Vec3f){p2.x,p2.y,z0});
                    Draw3D::vertex((Vec3f){p.x,p.y,z0}); p2=p+ruler.lvecs[2]*sc; Draw3D::vertex((Vec3f){p2.x,p2.y,z0});
                }else{
                    glColor3f(1.0,0.2,1.0); Draw3D::vertex((Vec3f){p.x,p.y,z0}); ruler.tilePoint( ip_, false, p );  p.add(off,off);
                    glColor3f(0.2,1.0,1.0); Draw3D::vertex((Vec3f){p.x,p.y,z0});
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
    testREQ = (Vec3d){ 1.487, 0.0006808, 0.0}; // H
    testPLQ = REQ2PLQ( testREQ, -1.6 );
    Quat4f * FFtot = new Quat4f[ W->gridFF.grid.getNtot() ];
    W->gridFF.evalCombindGridFF ( testREQ, FFtot );
    if(idebug>0) W->gridFF.grid.saveXSF( "FFtot_z.xsf",  (float*)FFtot, 4, 2, W->gridFF.natoms, W->gridFF.atypes, W->gridFF.apos );
    ogl_isosurf = glGenLists(1);
    glNewList(ogl_isosurf, GL_COMPILE);
    glShadeModel( GL_SMOOTH );
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    int nvert = renderSubstrate_( W->gridFF.grid, FFtot, W->gridFF.FFelec, isoVal, true, colorSclae );   printf("DEBUG renderGridFF() renderSubstrate() -> nvert= %i ", nvert );
    // ---- This seems still not work properly
    //int ntris=0;
    //glColor3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( W->gridFF.grid,  isoVal, (double*)FFtot, isoSurfRenderType,  3,2 );
    //glColor3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( W->gridFF.grid, -isoVal, (double*)FFtot, isoSurfRenderType,  3,2 );
    //glColor3f(0.,0.,1.); Draw3D::drawTriclinicBox( W->gridFF.grid.cell, (Vec3d){0.0, 0.0, 0.0}, (Vec3d){1.0, 1.0, 1.0} );
    //Draw3D::drawAxis(1.0);
    glEndList();
    delete [] FFtot;
    if(verbosity>0) printf( "... MolGUI::renderGridFF() DONE\n" );
}


void MolGUI::renderESP( Vec3d REQ){
    printf( "DEBUG MolGUI::renderESP() \n" ); //exit(0);
    glNewList(ogl_esp, GL_COMPILE);
    glShadeModel( GL_SMOOTH );
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    const NBsystem& nbmol = W->nbmol;
    int nvert = Draw3D::drawESP( nbmol.n, nbmol.ps, nbmol.REQs, REQ );
    glEndList();
};

void MolGUI::drawPi0s( float sc=1.0 ){
    const MMFFsp3& ff = W->ff;
    glBegin(GL_LINES);
    for(int ia=0; ia<ff.nnode; ia++){
        int* ngs = ff.aneighs + ia*ff.nneigh_max;
        for(int j=0; j<ff.nneigh_max; j++){
            int ing = ngs[j]; 
            if(ing<0){
                int ipi = -ing-1;
                Vec3f p=(Vec3f)ff.apos[ia];      glVertex3f(p.x,p.y,p.z);
                p.add_mul( ff.pi0s[ipi].f, sc);  glVertex3f(p.x,p.y,p.z);
                //printf("drawPi0s[%i,%i|%i] (%g,%g,%g)\n", ia, j, ipi, ff.pi0s[ipi].f.z, ff.pi0s[ipi].f.y, ff.pi0s[ipi].f.z );
            }
        }
    }
    glEnd();
}

void MolGUI::drawSystem( Vec3i ixyz ){
    glEnable(GL_DEPTH_TEST);
    bool bOrig = (ixyz.x==0)&&(ixyz.y==0)&&(ixyz.z==0);
    //printf( "DEBUG MolGUI::drawSystem() bViewMolCharges %i W->nbmol.REQs %li\n", bViewMolCharges, W->nbmol.REQs );
    //printf("DEBUG MolGUI::drawSystem()  bOrig %i W->bMMFF %i mm_bAtoms %i bViewAtomSpheres %i bViewAtomForces %i bViewMolCharges %i \n", bOrig, W->bMMFF, mm_bAtoms, bViewAtomSpheres, bViewAtomForces, bViewMolCharges  );
    if(W->bMMFF){
        //if(W->builder.bPBC){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC    ( W->ff.nbonds, W->ff.bond2atom, W->ff.apos, &W->builder.bondPBC[0], W->builder.lvec ); } 
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::bonds       ( W->ff.nbonds, W->ff.bond2atom, W->ff.apos );  // DEBUG
        if(W->builder.bPBC){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC          ( W->ff.nbonds, W->ff.bond2atom, W->ff.apos,  W->ff.pbcShifts );  
                             glColor3f(0.0f,0.0f,1.0f); Draw3D::pbcBondNeighLabels( W->ff.nbonds, W->ff.bond2atom, W->ff.apos,  W->ff.pbcShifts, fontTex3D,        0.007 );
        }else               { glColor3f(0.0f,0.0f,0.0f); Draw3D::bonds       ( W->ff.nbonds, W->ff.bond2atom, W->ff.apos );                                          
                             glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsLengths( W->ff.nbonds, W->ff.bond2atom, W->ff.apos, fontTex );                                
        }
        //Draw3D::atoms          ( W->ff.natoms, W->ff.apos, W->ff.atype, W->params, ogl_sph, 1.0, mm_Rsc, mm_Rsub );   
        //Draw3D::drawVectorArray( W->ff.natoms, W->ff.apos, W->ff.fapos, 100.0, 10000.0 );   
        //if(bOrig&&mm_bAtoms){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels       ( W->ff.natoms, W->ff.apos, fontTex3D                     ); }                    
        //if(bViewMolCharges && (W->nbmol.REQs!=0) ){ glColor3f(0.0,0.0,0.0);    Draw3D::atomPropertyLabel( W->ff.natoms,  (double*)W->nbmol.REQs, W->ff.apos, 3, 2, fontTex3D, 0.01 ); }   
        //void bondLabels( int n, const Vec2i* b2a, const Vec3d* apos, int fontTex, float sz=0.02 ){
        if(bOrig &&  bViewBondLabels     ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondLabels( W->ff.nbonds, W->ff.bond2atom, W->ff.apos, fontTex3D,        0.007              );       }

    }
    //W->nbmol.print();
    if(bViewAtomSpheres&&mm_bAtoms           ){                            Draw3D::atoms            ( W->nbmol.n, W->nbmol.ps, W->nbmol.atypes, W->params, ogl_sph, 1.0, mm_Rsc, mm_Rsub ); }
    //if(bViewAtomP0s                        ){ glColor3f(0.0f,1.0f,1.0f); Draw3D::drawVectorArray  ( W->nbmol.n, W->nbmol.ps, W->nbmol.fs, ForceViewScale, 10000.0 );  }
    if(bViewAtomForces                       ){ glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVectorArray  ( W->nbmol.n, W->nbmol.ps, W->nbmol.fs, ForceViewScale, 10000.0 );  }
    if(bOrig&&mm_bAtoms&&bViewAtomLabels     ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels       ( W->nbmol.n, W->nbmol.ps, fontTex3D,        0.007              );       }
    if(bViewMolCharges && (W->nbmol.REQs!=0) ){ glColor3f(0.0,0.0,0.0);    Draw3D::atomPropertyLabel( W->nbmol.n,  (double*)W->nbmol.REQs,  W->nbmol.ps, 3, 2, fontTex3D, 0.01 ); }
    if(W->ff.pi0s                            ){ glColor3f(0.0f,1.0f,1.0f); drawPi0s(1.0); }

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

void MolGUI::debug_scanSurfFF( int n, Vec3d p0, Vec3d p1, double sc ){
    //printf("========= MolGUI::scanFF\n");
    sc=10.0;
    //p0=(Vec3d){0.0,0.0,z0_scan}; p1=(Vec3d){5.0,5.0,z0_scan}; // Horizontal scan
    p0=(Vec3d){0.0,0.0,z0_scan}; p1=(Vec3d){0.0,5.0,z0_scan}; // Horizontal scan
    //p0=(Vec3d){0.0,z0_scan,0.0}; p1=(Vec3d){0.0,z0_scan,10.0,}; // Vertical scan
    Vec3d dp=p1-p0; dp.mul(1./n);
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        W->nbmol.ps[0]=p0+dp*i;
        W->nbmol.fs[0]=Vec3dZero;
        double E=0;
        if   (W->bGridFF){ E+= W->gridFF.eval     (1, W->nbmol.ps, W->nbmol.PLQs,  W->nbmol.fs );      }
        //else           { E+= W->nbmol .evalMorse   (W->surf, false,                        W->gridFF.alpha, W->gridFF.Rdamp ); }
        else             { E+= W->nbmol .evalMorsePBC(W->surf, W->gridFF.grid.cell, W->nPBC, W->gridFF.alpha, W->gridFF.Rdamp ); }
        Draw3D::vertex( W->nbmol.ps[0] ); Draw3D::vertex( W->nbmol.ps[0]+W->nbmol.fs[0]*sc );       // Force Vectro
        //Draw3D::vertex( W->nbmol.ps[0] ); Draw3D::vertex( W->nbmol.ps[0]+((Vec3d){0.0,0.0,E})*sc  );  // Energy -> z
        //Draw3D::vertex( W->nbmol.ps[0] ); Draw3D::vertex( W->nbmol.ps[0]+((Vec3d){0.0,E,0.0})*sc  );  // Energy -> x
    }
    glEnd();
}

void MolGUI::eventHandling ( const SDL_Event& event  ){
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

                /*
                case SDLK_KP_1: picked_lvec = &W->builder.lvec.a; break;
                case SDLK_KP_2: picked_lvec = &W->builder.lvec.b; break;
                case SDLK_KP_3: picked_lvec = &W->builder.lvec.c; break;

                case SDLK_KP_7: picked_lvec->x+=xstep; break;
                case SDLK_KP_4: picked_lvec->x-=xstep; break;

                case SDLK_KP_8: picked_lvec->y+=xstep; break;
                case SDLK_KP_5: picked_lvec->y-=xstep; break;

                case SDLK_KP_9: picked_lvec->z+=xstep; break;
                case SDLK_KP_6: picked_lvec->z-=xstep; break;
                */


                case SDLK_KP_PLUS:  z0_scan = z0_scan+0.1; break;
                case SDLK_KP_MINUS: z0_scan = z0_scan-0.1; break;
                

                case SDLK_KP_0: qCamera = qCamera0; break;

                case SDLK_COMMA:  which_MO--; printf("which_MO %i \n", which_MO ); break;
                case SDLK_PERIOD: which_MO++; printf("which_MO %i \n", which_MO ); break;
                //case SDLK_LESS:    which_MO--; printf("which_MO %i \n"); break;
                //case SDLK_GREATER: which_MO++; printf("which_MO %i \n"); break;

                //case SDLK_m: renderOrbital( which_MO ); break;
                //case SDLK_r: renderDensity(          ); break;
                case SDLK_s: W->saveXYZ( "out.xyz", "#comment", false ); break;
                case SDLK_p: saveScreenshot( frameCount ); break;
                
                //case SDLK_g: useGizmo=!useGizmo; break;
                //case SDLK_g: W->bGridFF=!W->bGridFF; break;
                //case SDLK_g: W->swith_gridFF(); break;
                //case SDLK_c: W->autoCharges(); break;

                case SDLK_g: W->bGridFF=!W->bGridFF; break;
                case SDLK_c: W->bOcl=!W->bOcl; break;
                case SDLK_m: W->swith_method(); break;

                case SDLK_a: bViewAtomSpheres=! bViewAtomSpheres; break;
                case SDLK_l: bViewAtomLabels =! bViewAtomLabels; break;
                case SDLK_b: bViewBondLabels =! bViewBondLabels; break;
                //case SDLK_q: W->autoCharges(); break;
                case SDLK_q: bViewMolCharges =! bViewMolCharges;  break;
                case SDLK_f: bViewAtomForces =! bViewAtomForces;  break;
                case SDLK_w: bViewSubstrate  =! bViewSubstrate;   break;


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

                case SDLK_LEFTBRACKET:
                    rotate( W->selection.size(), &W->selection[0], W->ff.apos, rotation_center, rotation_axis, +rotation_step );
                    break;
                case SDLK_RIGHTBRACKET:
                    rotate( W->selection.size(), &W->selection[0], W->ff.apos, rotation_center, rotation_axis, -rotation_step );
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
                        int ipick = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, W->ff.natoms, W->ff.apos );
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
                case SDL_BUTTON_RIGHT:
                    //ibpicked = -1;
                    break;
            }
            break;
        case SDL_WINDOWEVENT:
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE:
                    //SDL_Log("Window %d closed", event->window.windowID);
                    printf( "window[%i] SDL_WINDOWEVENT_CLOSE \n", id );
                    delete this;
                    printf( "window[%i] delete *this DONE \n", id );
                    return;
                    break;
            } break;
    };
    //AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void MolGUI::keyStateHandling( const Uint8 *keys ){
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
};
