
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
#include "MarchingCubes.h"
#include "MolecularDraw.h"
#include "GUI.h"
#include "EditorGizmo.h"
#include "SimplexRuler.h"
#include "AppSDL2OGL_3D.h"

int idebug=0;

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
    int ipicked    = -1; // picket atom 
    int ibpicked   = -1; // picket bond
    int iangPicked = -1; // picket angle
    Vec3d* picked_lvec = 0;
    int perFrame =  1;
    Quat4f qCamera0;

    bool bDoMM=true;//,bDoQM=true;
    bool bConverged = false;
    bool bRunRelax  = false;
    double cameraMoveSpeed = 1.0;
    bool useGizmo=true;
    bool bDrawHexGrid=true;
    bool bHexDrawing=false; 

    bool bPrepared_mm = false;
    bool bPrepared_qm = false;

    MolWorld_sp3* W=0;

    GUI gui;
    EditorGizmo  gizmo;
    SimplexRuler ruler; // Helps paiting organic molecules

    // ---- Visualization params
    int which_MO  = 7; 
    double mm_Rsc =  0.25;
    double mm_Rsub = 1.0;
    bool   mm_bAtoms = false;
    bool   isoSurfRenderType = 1;
    Vec3d testREQ,testPLQ;

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

	MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_, const char* ssmile );

	//int  loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int  loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );
    void tryLoadGridFF();
    //void makeGridFF   (bool recalcFF=false, bool bRenderGridFF=true);
    void renderGridFF();

	void drawSystem( Vec3d ixyz );
    //void drawSystemQMMM();
    //void renderOrbital(int i, double iso=0.1);
    //void renderDensity(       double iso=0.1);
	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

    //void InitQMMM();
    void initGUI();
    void drawingHex(double z0);
};

//=================================================
//                   INIT()
//=================================================


void MolGUI::initGUI(){

    GUI_stepper ylay;
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

MolGUI::MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_, const char* ssmile ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    if(W_==0){
        W = new MolWorld_sp3();
        if(ssmile){
            printf("MolGUI() init_smile(%s)\n", ssmile);
            W->initWithSMILES(ssmile);
        }else{
            W->ini_in_dir();
        }
        W->eval();
        W->ff.checkNaNs();
    }

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
    initGUI();
}


//=================================================
//                   DRAW()
//=================================================

void MolGUI::draw(){
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);

    if(frameCount==1){ qCamera.pitch( M_PI );  qCamera0=qCamera; }
    if(bRunRelax){ W->MDloop(perFrame); }

    // --- Mouse Interaction / Visualization
	ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
    Draw3D::drawPointCross( ray0, 0.1 );        // Mouse Cursor 
    if(ipicked>=0) Draw3D::drawLine( W->ff.apos[ipicked], ray0); // Mouse Dragging Visualization
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

    if(ogl_isosurf)viewSubstrate( 2, 2, ogl_isosurf, W->gridFF.grid.cell.a, W->gridFF.grid.cell.b, W->gridFF.shift );
    //if(bDoQM)drawSystemQMMM();
    if(bDoMM)if(W->builder.bPBC){ Draw3D::drawPBC( (Vec3i){2,2,0}, W->builder.lvec, [&](Vec3d ixyz){drawSystem(ixyz);} ); } else { drawSystem({0,0,0}); }
    for(int i=0; i<W->selection.size(); i++){ int ia = W->selection[i];
        glColor3f( 0.f,1.f,0.f ); Draw3D::drawSphereOctLines( 8, 0.3, W->ff.apos[ia] );     }
    //if(iangPicked>=0){
    //    glColor3f(0.,1.,0.);      Draw3D::angle( W->ff.ang2atom[iangPicked], W->ff.ang_cs0[iangPicked], W->ff.apos, fontTex3D );
    //}

    if(W->bSurfAtoms)Draw3D::atomsREQ( W->surf.n, W->surf.ps, W->surf.REQs, ogl_sph, 1., 1., 0. );

    if(useGizmo){
        gizmo.draw();
    }
    if(bHexDrawing)drawingHex(5.0);
};

void MolGUI::drawHUD(){
    glDisable ( GL_LIGHTING );
    gui.draw();
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

void MolGUI::renderGridFF(){
    int iatom = 11;
    testREQ = (Vec3d){ 1.487, 0.0006808, 0.0}; // H
    testPLQ = REQ2PLQ( testREQ, -1.6 );
    Vec3d * FFtot = new Vec3d[ W->gridFF.grid.getNtot() ];
    W->gridFF.evalCombindGridFF            ( testREQ, FFtot );
    if(idebug>1) saveXSF( "FFtot_z.xsf",  W->gridFF.grid, FFtot, 2, W->gridFF.natoms, W->gridFF.apos, W->gridFF.atypes );
    ogl_isosurf = glGenLists(1);
    glNewList(ogl_isosurf, GL_COMPILE);
    glShadeModel( GL_SMOOTH );
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    renderSubstrate_( W->gridFF.grid, FFtot, W->gridFF.FFelec, 0.01, true, 0.1);
    Draw3D::drawAxis(1.0);
    glEndList();
    delete [] FFtot;
}

void MolGUI::drawSystem( Vec3d ixyz ){
    bool bOrig = (ixyz.x==0)&&(ixyz.y==0)&&(ixyz.z==0);
    if(W->builder.bPBC){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC( W->ff.nbonds, W->ff.bond2atom, W->ff.apos, &W->builder.bondPBC[0], W->builder.lvec ); } // DEBUG
    else               { glColor3f(0.0f,0.0f,0.0f); Draw3D::bonds   ( W->ff.nbonds, W->ff.bond2atom, W->ff.apos );                                          }
    if(bOrig&&mm_bAtoms){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels( W->ff.natoms, W->ff.apos, fontTex3D                     ); }                     //DEBUG
    Draw3D::atoms( W->ff.natoms, W->ff.apos, W->ff.atype, W->params, ogl_sph, 1.0, mm_Rsc, mm_Rsub );       //DEBUG
}

void MolGUI::saveScreenshot( int i, const char* fname ){
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

                case SDLK_KP_1: picked_lvec = &W->builder.lvec.a; break;
                case SDLK_KP_2: picked_lvec = &W->builder.lvec.b; break;
                case SDLK_KP_3: picked_lvec = &W->builder.lvec.c; break;

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

                //case SDLK_m: renderOrbital( which_MO ); break;
                //case SDLK_r: renderDensity(          ); break;
                case SDLK_c: saveScreenshot( frameCount ); break;

                case SDLK_g: useGizmo=!useGizmo; break;
                case SDLK_f:
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
                        ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, W->ff.natoms, W->ff.apos );
                        W->selection.clear();
                        if(ipicked>=0){ W->selection.push_back(ipicked); };
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
        case SDL_WINDOWEVENT:
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE:
                    //SDL_Log("Window %d closed", event->window.windowID);
                    printf( "window[%i] SDL_WINDOWEVENT_CLOSE \n", id );
                    delete this;
                    printf( "window[%i] delete this done \n", id );
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
    if( keys[ SDL_SCANCODE_LEFT  ] ){ cam.pos.add_mul( cam.rot.a, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ cam.pos.add_mul( cam.rot.a,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_UP    ] ){ cam.pos.add_mul( cam.rot.b,  cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ cam.pos.add_mul( cam.rot.b, -cameraMoveSpeed ); }
    //AppSDL2OGL_3D::keyStateHandling( keys );
};
