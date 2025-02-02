
/*

NOTE:

WORKING COMMITS - ( Relax CH4 and C4H12 )
commit d8b1156f6ef8f283785dae1fee71105d591a3280    2021-Apr-26  testing CLCFGO.h vs eFF.h for Hydrogen atom with electron radius 0.5A…

NOT WORKING COMMITS
commit dae1d0b16d3b892f0c3d982c5f277df38bfb4179    2021-Jun-01    CLCFGO : option to out-project force components which breaks normaliz…
commit 94a94e956acad8e3d23a54acbd0f715fe0d1f827    2021-May-05    CLCFGO : tested H2 and H2O molecule with respect to eFF; total energy…


*/




#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include  <functional>


#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "testUtils.h"
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "SDL_utils.h"
#include "GUI.h"
#include "Plot2D.h"
#include "Console.h"

//#include "MMFF.h"

#define R2SAFE  1.0e-8f

int i_DEBUG = 0;

int DEBUG_i = 0;
int DEBUG_j = 0;

Vec3d* DEBUG_fe_ae =0;
Vec3d* DEBUG_fa_ae =0;
Vec3d* DEBUG_fe_ee =0;
Vec3d* DEBUG_fa_aa =0;

#include "DynamicOpt.h"
#include "InteractionsGauss.h"
#include "eFF.h"

//#include "InteractionsGauss_old.h"
//#include "eFF_old.h"

//#include "e2FF.h" // old currently not working

#include "Forces.h"

#include "eFF_plots.h"

#include "argparse.h"
LambdaDict funcs;
#ifdef WITH_LUA
//#include "LuaUtils.h"
//#include "LuaClass.h"
#include "LuaHelpers.h"
lua_State  * theLua=0;
#endif // WITH_LUA

void cleanDebugForce(int ne, int na){
    for(int i=0; i<ne; i++){
        DEBUG_fe_ae[i]=Vec3dZero;
        DEBUG_fe_ee[i]=Vec3dZero;
    }
    for(int i=0; i<na; i++){
        DEBUG_fa_ae[i]=Vec3dZero;
        DEBUG_fa_aa[i]=Vec3dZero;
    }
};

void applyCartesianBoxForce( const Vec3d& pmin, const Vec3d& pmax,const Vec3d& k, int n, const Vec3d* ps, Vec3d* fs ){
    for(int i=0;i<n; i++){ boxForce( ps[i], fs[i], pmin, pmax, k ); }
   //for(int i=0;i<n; i++){     if(i==0){ boxForce( ps[i], fs[i], pmin, pmax, k ); printf( " atom[%i] p(%g,%g,%g) f(%g,%g,%g) \n", i, ps[i].x, ps[i].y, ps[i].z,   fs[i].x, fs[i].y, fs[i].z ); }  }
}

void applyParabolicPotential( const Vec3d& p0, const Vec3d& k, int n, const Vec3d* ps, Vec3d* fs ){
    for(int i=0;i<n; i++){ Vec3d dp=ps[i]-p0; fs[i].add( dp*dp*k ); }
}

void draw2Dfunc( const EFF& ff, Vec2i ns, Vec2d spanx, Vec2d spany, double z_cut, Vec3i axes, bool bAtom=true, bool bElectron=true, bool bCoulomb=true, int spin=0.0, double s=0.0, double Q=1.0 ){
    double dx = (spanx.y-spany.x)/ns.x;
    double dy = (spany.y-spany.x)/ns.y;
    Vec3d p0=Vec3dZero;
    p0.array[axes.x] = spanx.x;
    p0.array[axes.y] = spany.x;
    double sign_Q = (Q>0)?1.0:-1.0;
    double absQ   = fabs(Q);
    Q*=bCoulomb; sign_Q*=bCoulomb;
    //printf( "draw2Dfunc Q %g sign_Q %g absQ %g \n", Q, sign_Q, absQ );
    for(int iy=0; iy<=ns.y; iy++){
        Vec3d p = p0;
        p.array[axes.y] = p0.y + iy*dy;
        glBegin(GL_LINE_STRIP);
        for(int ix=0; ix<=ns.x; ix++){
            p.array[axes.x] = p0.x + ix*dx;
            p.array[axes.z] = z_cut;
            double E = 0;
            if(bAtom    ){ E += ff.atomsPotAtPoint   ( p, s,      Q            ); }
            if(bElectron){ E += ff.electronPotAtPoint( p, s,-sign_Q, spin )*absQ; }
            p.array[axes.z] = E;
            glVertex3f( p.x, p.y, p.z );
        }
        glEnd();
    }
    glBegin(GL_LINE_LOOP);
    //p0.array[axes.z]=0;
    p0.array[axes.y]=spany.x; p0.array[axes.x]=spanx.x; glVertex3f( p0.x, p0.y, p0.z );
    p0.array[axes.y]=spany.x; p0.array[axes.x]=spanx.y; glVertex3f( p0.x, p0.y, p0.z );
    p0.array[axes.y]=spany.y; p0.array[axes.x]=spanx.y; glVertex3f( p0.x, p0.y, p0.z );
    p0.array[axes.y]=spany.y; p0.array[axes.x]=spanx.x; glVertex3f( p0.x, p0.y, p0.z );
    glEnd();
}

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:
    //enum PICK_MODE{ PICK_NONE, ATOM, ELECTRON };
    enum PICK_MODE{ ATOM, ELECTRON, N_PICK_MODE };
    bool bRun = false;

    //E2FF ff2;
    EFF  ff;
    DynamicOpt opt;

    Vec3d ray0;
    //int perFrame = 1;
    int perFrame = 10;
    //int perFrame = 20;
    //int perFrame = 50;

    Vec2i field_ns;
    Vec2d Erange;
    double E0,Espread;
    Vec3d  * field_ps=0;
    double * field_Es=0;
    std::function<void   (const Vec3d& p, Vec3d& f)>  FFfunc;
    std::function<double (const Vec3d& p)          >  Efunc ;

    bool bViewForce   = false;
    bool bDrawPointLabels = true;
    //bool bDrawPlots   = true;
    bool bDrawPlots   = false;
    bool bDrawPlotXY  = true;
    bool bDrawPlotXZ  = true;
    bool bDrawPlotYZ  = true;
    bool bDrawObjects = true;
    bool bMapElectron = false;
    int ipicked  = -1;
    int PickMode = ELECTRON;

    double xy_height = 0.0;
    double yz_height = 0.0;
    double xz_height = 0.0;
    double scale_V   = 1.0;
    double electron_size = 1.0;
    int    electron_spin = 1;
    int    electron_Q    = 1;

    // DEBUG STUFF
    GLint ogl_fs = 0;
    GLint oglSph = 0;

    // ToDO: maybe move to globals.h
    int    fontTex=-1;
    int  fontTex3D=-1;

    bool bConsole=false;
    GUI gui;
    Console console;
    Plot2D plot1;
    //DropDownList* panel_Frags=0;
    //GUIPanel*     panel_iMO  =0;
    CheckBoxList* panel_Plot=0;

    // ========== Functions

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

    void initEFFsystem( const char * fname, bool bTestEval=true, bool bDebug=false );
    void init2DMap( int n, double dx );
    void initWiggets();

};

void TestAppRARFF::initWiggets(){
    printf( "TestAppRARFF::initWiggets() \n" );
    // MultiPanel(const std::string& caption, int xmin, int ymin, int xmax, int dy, int nsubs)
    //GUI_stepper ylay;
    GUI_stepper gx(0,4);
    GUIPanel*   p=0;
    MultiPanel* mp=0;

    // mp = new MultiPanel( "Plot", gx.x0, 10, gx.x1, 0, 3, false, true, false, false );   gui.addPanel( mp );
    // p=mp->subs[0]; p->caption = "view"; p->command=[&](GUIAbstractPanel* p){ bDrawPlots=!bDrawPlots; };
    // p=mp->subs[1]; p->caption = "grid"; p->command=[&](GUIAbstractPanel* p){ plot1.bGrid=!plot1.bGrid;                             plot1.redraw=true; };
    // p=mp->subs[2]; p->caption = "axes"; p->command=[&](GUIAbstractPanel* p){ plot1.bAxes=!plot1.bAxes; plot1.bTicks=!plot1.bTicks; plot1.redraw=true; };

    // void initCheckBoxList( int xmin_, int ymin_, int xmax_, int dy=fontSizeDef*2 );
    //CheckBox* b;
    CheckBoxList* chk = new CheckBoxList( gx.x0, 10, gx.x1 );
    gui.addPanel( chk );
    panel_Plot = chk; 
    chk->caption = "Plot"; chk->bgColor = 0xFFE0E0E0;
    chk->addBox( "view", &bDrawPlots  );
    chk->addBox( "grid", &plot1.bGrid );
    chk->addBox( "axes", &plot1.bAxes );

    gx.x1+=5; gx.step( 15 );
    //mp = new MultiPanel( "PlotXY", gx.x0, 10, gx.x1, 0, 5, true, true, false, true, true );   gui.addPanel( mp );
    // p=mp->subs[0]; p->caption = "view";    p->command=[&](GUIAbstractPanel* p){ bDrawPlotXY=!bDrawPlotXY;            }; p->isSlider=false; p->viewVal=false;
    // p=mp->subs[1]; p->caption = "z_cut: "; p->command=[&](GUIAbstractPanel* p){ xy_height    =((GUIPanel*)p)->value; }; p->setRange(-10.0,10.0); p->setValue(0.0);
    // p=mp->subs[1]; p->caption = "scale: "; p->command=[&](GUIAbstractPanel* p){ scale_V      =pow(10.0,((GUIPanel*)p)->value); }; p->setRange(-2.0,2.0); p->setValue(0.0);
    // p=mp->subs[2]; p->caption = "size : "; p->command=[&](GUIAbstractPanel* p){ electron_size=((GUIPanel*)p)->value; }; p->setRange( 0.0 ,2.0);  p->setValue(1.0);
    // p=mp->subs[3]; p->caption = "spin : "; p->command=[&](GUIAbstractPanel* p){ electron_spin=((GUIPanel*)p)->value; }; p->setRange(-1.0 ,1.0);  p->setValue(1.0); p->isInt=true;
    // p=mp->subs[4]; p->caption = "Q    : "; p->command=[&](GUIAbstractPanel* p){ electron_Q   =((GUIPanel*)p)->value; }; p->setRange( 0.0 ,1.0);  p->setValue(1.0); p->isInt=true;
    
    mp = new MultiPanel( "PlotXY", gx.x0, 10, gx.x1, 0, -5, true, true, false, true, true );   gui.addPanel( mp );
    // GUIPanel* addPanel( const std::string& caption, Vec3d vals{min,max,val}, bool isSlider, bool isButton, bool isInt, bool viewVal, bool bCmdOnSlider );
    mp->addPanel( "view",    { 0.0 ,1.0 ,0.0}, 0,1,0,0,0 )->command=[&](GUIAbstractPanel* p){ bDrawPlotXY=!bDrawPlotXY;            };
    mp->addPanel( "z_cut: ", {-10.0,10.0,0.0}, 1,1,0,1,1 )->command=[&](GUIAbstractPanel* p){ xy_height    =((GUIPanel*)p)->value; };
    mp->addPanel( "scale: ", {-2.0 ,2.0 ,0.0}, 1,1,0,1,1 )->command=[&](GUIAbstractPanel* p){ scale_V      =pow(10.0,((GUIPanel*)p)->value); };
    mp->addPanel( "size : ", { 0.0 ,2.0 ,1.0}, 1,1,0,1,1 )->command=[&](GUIAbstractPanel* p){ electron_size=((GUIPanel*)p)->value; };
    mp->addPanel( "spin : ", {-1.0 ,1.0 ,1.0}, 1,1,1,1,1 )->command=[&](GUIAbstractPanel* p){ electron_spin=((GUIPanel*)p)->value; };
    mp->addPanel( "Q    : ", { 0.0 ,1.0 ,1.0}, 1,1,1,1,1 )->command=[&](GUIAbstractPanel* p){ electron_Q   =((GUIPanel*)p)->value; };
    
    printf( "TestAppRARFF::initWiggets() END \n" );
}

void TestAppRARFF::initEFFsystem( const char * fname, bool bTestEval, bool bDebug ){

    ff.loadFromFile_fgo( fname );

    if(bDebug){
        DEBUG_fe_ae = new Vec3d[ff.ne];
        DEBUG_fa_ae = new Vec3d[ff.na];
        DEBUG_fe_ee = new Vec3d[ff.ne];
        DEBUG_fa_aa = new Vec3d[ff.na];
    }
    //setGeom(ff);
    //double sz = 0.2;
    //for(int i=0; i<ff.na; i++){ ff.apos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    //for(int i=0; i<ff.ne; i++){ ff.epos[i].add( randf(-sz,sz),randf(-sz,sz),randf(-sz,sz) );  }
    //ff.autoAbWs( default_aAbWs, default_eAbWs );
    //VecN::set(ff.ne,4.0,ff.esize);

    // ==== Test Eval

    //makePlots( plot1, ff );
    //ff.loadFromFile_xyz( fname  );
    //init2DMap( 100, 0.1 );

    // ==== Bind Optimizer
    opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, 0, ff.fDOFs, 0 );
    opt.cleanVel( );
    //opt.initOpt( 0.01, 0.2 );
    opt.initOpt( 0.0015, 0.01 );
    opt.f_limit = 1000.0;

    //ff.iPauliModel = 0; // dens overlap
    ff.iPauliModel = 1; // addPauliGauss   from the article using  KRSrho
    //ff.iPauliModel = 2; // addPauliGaussVB valence bons
    ff.info();

    // {
    // Vec3d p = ff.epos[0];
    // double dx = 0.1;
    // double sz = 0.5;
    // double sign_Q = 0.0;
    // double absQ   = 1.0;
    // int    spin   = -1;
    // for(int ix=0; ix<10; ix++){ 
    //     double E = ff.electronPotAtPoint( p, sz, sign_Q, spin )*absQ;
    //     //printf( "E[%i] x=%g  E=%g absQ=%g \n", ix, ix*dx, E, absQ );
    //     p.x += 0.1;
    // }
    // exit(0);
    // }

    if(bTestEval){
        double E = ff.eval();
        //printf( "E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", E, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
        ff.printEnergies();
        ff.printEnergies();
        //printf( " test_eFF exits ... \n" ); exit(0);
    }

    initWiggets();
};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_, " test_eFF " ) {
    printf("EFFapp.cpp TestAppRARFF start \n");
    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    plot1.fontTex=fontTex;

    //checkDerivs( ff.KRSrho );   // exit(0);
    //makePlots( plot1, ff );     // exit(0);
    //makePlots2( plot1 );        // exit(0);
    //checkDerivs2();             // exit(0);

    //ff.bEvalAECoulomb = 0;
    //ff.bEvalAEPauli   = 0;
    //ff.bEvalCoulomb   = 0;
    //ff.bEvalPauli     = 0;
    //ff.bEvalKinetic   = 0;
    //ff.bEvalAA        = 0;

    oglSph=Draw::list(oglSph);
    //Draw3D::drawSphere_oct(3,1.0d,(Vec3d){0.,0.,0.});
    Draw3D::drawSphere_oct(5,1.0,(Vec3d){0.,0.,0.});
    glEndList();

    plot1.init();
    plot1.fontTex = fontTex;
    plot1.add( new DataLine2D( 200, -10.0, 0.1, 0xFF0000FF, "Vatom" ) );
    plot1.bAxes  = true;
    plot1.bTicks = true;
    plot1.update();
    plot1.render();
    plot1.view();

    console.init( 256, window );
    console.callback = [&](const char* s){ printf( "console.callback(%s)\n", s ); return 0; };
    console.fontTex = fontTex;
    printf("EFFapp.cpp TestAppRARFF end \n");

#ifdef WITH_LUA
    console.callback = [&](const char* str)->bool{
       lua_State* L=theLua;
        printf( "console.lua_callback: %s\n", str );
        if (luaL_dostring(L, str) != LUA_OK)[[unlikely]]{
            // If there's an error, print it
            //fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
            printf( "Error: %s\n", lua_tostring(L, -1) );
            lua_pop(L, 1);  // Remove error message from the stack
            return false;
        }
        return true;
    };
#endif // WITH_LUA
    
}

void TestAppRARFF::draw(){
    printf( " ==== frame %i \n", frameCount );
    // TADY SE VYKRESULJI JEDNOTLIVY FRAMES

    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    //return;

    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    double vminOK = 1e-6;
    double vmaxOK = 1e+3;

    //perFrame=10; // ToDo : why it does not work properly for  perFrame>1 ?
    //perFrame = 1;
    double sum = 0;
    if(bRun){
        double F2 = 1.0;
        double Etot;
        for(int itr=0;itr<perFrame;itr++){
            //printf( " ==== frame %i i_DEBUG  %i \n", frameCount, i_DEBUG );

            ff.clearForce();
            //ff.clearForce_noAlias();
            cleanDebugForce( ff.ne, ff.na);
            //VecN::sum( ff.ne*3, ff.eforce );

            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,50.0}, ff.na, ff.apos, ff.aforce );
            //applyCartesianBoxForce( {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0,0,5.0},  ff.ne, ff.epos, ff.eforce );

            //applyParabolicPotential( {0.0,0.0,0.0}, {0.1,0.1,0.1}, ff.na, ff.apos, ff.aforce );
            //applyParabolicPotential( {0.0,0.0,0.0}, {0.1,0.1,0.1}, ff.ne, ff.epos, ff.eforce );

            Etot = ff.eval();
            
            //ff.apos[0].set(.0);
            //checkFinite( ff, vminOK, vmaxOK );

            //printf( "fa1(%g,%g,%g) fe1(%g,%g,%g)\n", fa1.x,fa1.x,fa1.x,   fe1.x,fe1.x,fe1.x );

            if(ipicked>=0){                                // Drag atoms by Mouse
                if     (PickMode==ATOM    ){ Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 ); ff.aforce[ipicked].add( f ); }
                else if(PickMode==ELECTRON){ Vec3d f = getForceSpringRay( ff.epos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 ); ff.eforce[ipicked].add( f ); }
            }

            //VecN::set( ff.na*3, 0.0, (double*)ff.aforce );  // FIX ATOMS
            //VecN::set( ff.ne  , 0.0, ff.fsize  );           // FIX ELECTRON SIZE
            //VecN::set( ff.ne*3, 0.0, (double*)ff.eforce );  // FIX ELECTRON POS
            //if(bRun)ff.move_GD(0.001 );
            //ff.move_GD( 0.0001 );
            //if(bRun) ff.move_GD_noAlias( 0.0001 );
            F2 = opt.move_FIRE();

            //checkFinite( ff, vminOK, vmaxOK );

            //printf( "frame[%i] E %g pa[0](%g,%g,%g) pe[0](%g,%g,%g) \n", frameCount, E,   ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z );
            //printf( "frame[%i] E %g pe[0](%g,%g,%g) s %g fe[0](%g,%g,%g) fs %g \n", frameCount, E,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z,  ff.esize[0],   ff.eforce[0].x,ff.eforce[0].y,ff.eforce[0].z, ff.fsize[0] );

            //printf( "frame[%i] E %g pa[0](%g,%g,%g) pe[0](%g,%g,%g) s %g \n", frameCount, E, ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.epos[0].x,ff.epos[0].y,ff.epos[0].z,  ff.esize[0] );
            //printf( "frame[%i] E %g lHH %g lH1e1 %g se1 %g \n", frameCount, E, (ff.apos[0]-ff.apos[1]).norm(),   (ff.apos[0]-ff.epos[0]).norm(), ff.esize[0] );

            //printf( "frame[%i,%i] E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", frameCount, itr, Etot, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
            //printf( "frame[%i,%i] " );  ff.printEnergies();
            //printf( "E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", E, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
            //printf( "=== %i %i frame[%i][%i] |F| %g \n", ff.na, ff.ne, frameCount, itr, sqrt(F2) );
        }
        if( F2 < 1e-6 ){
            //printf( "Finished: E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", Etot, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
            //printf( "Finished:"); ff.printEnergies();
            //ff.info();
            //printDistFormAtom( ff.na, ff.apos, 0 );
            //bRun=false;
            //ff.save_xyz( "data/eff_relaxed.xyz" );
        }
    }

    glColor3f(0.0,0.5,0.0);
    Draw3D::drawPointCross( ray0, 0.1 );
    //if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);
    if(ipicked>=0){                                // Drag atoms by Mouse
        //glColor3f(0.0,0.5,0.0);
        if     (PickMode==ATOM    ){
            printf( "ipicked %i pi(%g,%g,%g) ray0(%g,%g,%g) \n", ipicked, ff.apos[ipicked].x, ff.apos[ipicked].y, ff.apos[ipicked].z, ray0.x, ray0.y, ray0.z );
            Draw3D::drawLine( ff.apos[ipicked], ray0); 
            }
        else if(PickMode==ELECTRON){ Draw3D::drawLine( ff.epos[ipicked], ray0); }
    }

    drawEFF( ff, oglSph, 1.0*bViewForce, 0.1, 0.1, 1.5 );
    if(bDrawPointLabels){
        glColor3f(1.0,0.0,0.0); Draw3D::pointLabels( ff.na, ff.apos, fontTex, 0.05 );
        //glColor3f(0.0,0.5,0.5); Draw3D::pointLabels( ff.ne, ff.epos, fontTex, 0.05 );
    }

    if(bDrawPlots){
        plotAtomsPot( ff, plot1.lines[0], (Vec3d){0.0,0.0,0.0}, (Vec3d){1.0,0.0,0.0}, -0.2, 0.1 );
        // plot1.bGrid=false;
        // plot1.bAxes=false;
        // plot1.bTicks=false;
        //plot1.update();
        //plot1.render();
        plot1.tryRender();
        plot1.view();
    }

    if( bDrawPlotXY ){
        //draw2Dfunc( ff, {20,20}, {-10.,10.0}, {-10.,10.0}, xy_height, {0,1,2}, true, true, true, 1.0, 0.5, 1.0 );
        draw2Dfunc( ff, {100,100}, {-10.,10.0}, {-10.,10.0}, xy_height, {0,1,2}, true, true, fabs(electron_Q)>1e-6, electron_spin, electron_size, scale_V );
    }

    //printf( "e[0] r %g s %g \n", ff.epos[0].norm(), ff.esize[0] );

    //printf( "r07 r %g s %g \n", ff.epos[0].norm(), ff.esize[0] );

    // --- Constrain in Z
    //double Kz = 1.0;
    //for(int i=0; i<ff.na; i++){  ff.aforce[i].z += -Kz*ff.apos[i].z;  };
    //for(int i=0; i<ff.ne; i++){  ff.eforce[i].z += -Kz*ff.epos[i].z;  };
    //printf( "na %i ne %i \n", ff.na, ff.ne );

    //Vec3d d = ff.apos[0]-ff.apos[1];

    glCallList(ogl_fs);
    //Draw3D::drawColorScale( 20, {0.0,0.0,0.0}, Vec3dY, Vec3dX, Draw::colors_rainbow, Draw::ncolors );
    //printf( "apos (%g,%g,%g) \n", ff.apos[0].x, ff.apos[0].y, ff.apos[0].z );


    char strtmp[256];
    double Qsz = 0.05;
    double fsc = 1.0;

    //for(int i=0; i<ff.ne; i+=2){
    //    Draw3D::drawLine(ff.epos[i],ff.epos[i+1] );
    //}

    //exit(0);

    //ff.aforce[0].set(0.);
    //ff.aforce[1].set(0.);
    //if(bRun) ff.move_GD( 0.01 );

    //glDisable(GL_DEPTH_TEST);
    //plot1.view();

};


void TestAppRARFF::drawHUD(){
    gui.draw();

    glPushMatrix();
    glTranslatef( 10.0,HEIGHT-20.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    //Draw::drawText( "AHOJ ", fontTex, fontSizeDef, {100,20} );
    //int nstr=2048;
	//char str[nstr];
	char* s=tmpstr;
	s+=ff.Eterms2str(s);
	ff.orbs2str(s);
    Draw::drawText( tmpstr, fontTex, fontSizeDef, {100,20} );
    glPopMatrix();

    if(bConsole) console.draw();

}

/*
void TestAppRARFF::mouseHandling( ){
    gui.onEvent( mouseX, mouseY, event );

    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {

    }
    AppSDL2OGL_3D::mouseHandling( );
};
*/

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    Vec3d pa0;
    GUIAbstractPanel* clicked_panel = gui.onEvent( mouseX, mouseY, event );
    if( clicked_panel ){ 
        int ichanged =  clicked_panel->clearChanged();
        if( clicked_panel==panel_Plot ){ if( (ichanged==1)||(ichanged==2)  ) plot1.redraw=true; }
    };
    switch( event.type ){
        case SDL_KEYDOWN :
            if (bConsole){ bConsole=console.keyDown( event.key.keysym.sym ); }
            else switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_BACKQUOTE:{ bConsole = !bConsole;}break;   // ` SDLK_ for key '`' 
                case SDLK_p: PickMode=(PickMode+1)%N_PICK_MODE; break;
                case SDLK_f: bViewForce=!bViewForce; break;
                case SDLK_q: bDrawPointLabels=!bDrawPointLabels; break;

                case SDLK_i: ff.info(); break;
                case SDLK_LEFTBRACKET :  Espread *= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_RIGHTBRACKET:  Espread /= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_e: bMapElectron=!bMapElectron; break;

                // plot // plot/grid options
                case SDLK_g:  
                //plot1.bGrid=false;
                // plot1.bAxes=false;
                // plot1.bTicks=false;


                case SDLK_m:{
                    pa0 = ff.apos[ipicked];
                    sampleScalarField( Efunc, field_ns, {-5.0,-5.0,+0.1}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, field_ps, field_Es, Erange );
                    E0 = field_Es[0];
                    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );
                    ff.apos[ipicked]= pa0;
                    }break;
                case SDLK_SPACE: bRun = !bRun; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_MOUSEBUTTONDOWN:{
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:{
                    if     ( PickMode==PICK_MODE::ATOM     ){ ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.na, ff.apos, 0 ); if(ipicked>=0)printf("picked atom %i\n", ipicked); }
                    else if( PickMode==PICK_MODE::ELECTRON ){ ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.ne, ff.epos, 0 ); if(ipicked>=0)printf("picked elec %i\n", ipicked); }
                    }break;
                case SDL_BUTTON_RIGHT:{
                    /*
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    printf( "remove atom %i \n", ipicked );
                    ff.ignoreAtoms[ ipicked ] = true;
                    */
                    }break;
                }
            } break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            } break;
    }
    AppSDL2OGL::eventHandling( event );
}

void TestAppRARFF::keyStateHandling( const Uint8 *keys ){
    if(bConsole){ return; }
    //if( keys[ SDL_SCANCODE_R ] ){}
	AppSDL2OGL_3D::keyStateHandling( keys );
};

void TestAppRARFF::init2DMap( int n, double dx ){

    int ipicked = 1;

    FFfunc = [&](const Vec3d& p, Vec3d& f)->void{
        if(bMapElectron){ ff.epos[ipicked] = p; }else{ ff.apos[ipicked] = p; }
         // force on one electron
        ff.clearForce();
        ff.eval();
        if(bMapElectron){ f = ff.eforce[ipicked]; }else{ f = ff.aforce[ipicked]; }
    };

    Efunc = [&](const Vec3d& p)->double{
        //ff.apos[ipicked] = p; // force on one electron
        if(bMapElectron){ ff.epos[ipicked] = p; }else{ ff.apos[ipicked] = p; }
        ff.clearForce();
        double E = ff.eval();
        return E;
    };

    field_ns = {n,n};

    Vec3d pa0;// = ff.apos[ipicked];
    if(bMapElectron){ pa0 = ff.epos[ipicked]; }else{ pa0 = ff.apos[ipicked]; }

    sampleScalarField( Efunc, field_ns, {-(dx*n)/2,-(dx*n)/2,+0.1}, {dx,0.0,0.0}, {0.0,dx,0.0}, field_ps, field_Es, Erange );
    if(bMapElectron){ ff.epos[ipicked] = pa0; }else{ ff.apos[ipicked] = pa0; }

    E0 = field_Es[0];
    printf( "val_range: %g %g %g \n", Erange.x, Erange.y, field_Es[0] );
    Espread = 3.0;

    ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );

    //return ogl_fs;
}

// ===================== MAIN

TestAppRARFF* app=0;

#ifdef WITH_LUA
int l_insertQuickCommand(lua_State *L){
    const char* s = Lua::getString(L,1);
    printf( "l_insertQuickCommand `%s`\n", s );
    app->console.quick_tab.table.push_back( s );
    app->console.quick_tab.sort();
    printf( "l_insertQuickCommand table.size() %i \n", app->console.quick_tab.table.size() );
    for(int i=0; i<app->console.quick_tab.table.size(); i++){
        printf( "l_insertQuickCommand[%i] `%s`\n", i, app->console.quick_tab.table[i].c_str() );
    }
    return 1; // number of return values to Lua environment
}

int initMyLua(){
    printf( "initMyLua()\n" );
    theLua         = luaL_newstate();
    lua_State  * L = theLua;
    luaL_openlibs(L);
    // lua_register(L, "fix",   l_fixAtom      );
    // lua_register(L, "natom", l_getAtomCount );
    // lua_register(L, "apos",  l_getAtomPos   );
    // lua_register(L, "run",   l_toggleStop   );
    lua_register(L, "command", l_insertQuickCommand  );
    printf( "initMyLua() DONE\n" );
    return 1;
}
#endif // WITH_LUA

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);
    app = new TestAppRARFF( junk , DM.w-100, DM.h-100 );
	//app = new TestAppRARFF( junk , 800, 600 );

    funcs["-f"]={1,[&](const char** ss){ app->initEFFsystem( ss[0], true, true ); }}; // molecule as .fgo
    
    //app->initEFFsystem( const char * fname, true, false );

#ifdef WITH_LUA
    initMyLua();
    funcs["-lua"]={1,[&](const char** ss){ if( Lua::dofile(theLua,ss[0]) )[[unlikely]]{ printf( "ERROR in funcs[-lua] no such file %s => exit()\n", ss[0] ); exit(0); }; }};
#endif // WITH_LUA
    //char str[40];  sprintf(str,  );
	//SDL_SetWindowTitle( app->child_windows[1]->window, " test_eFF " );
    process_args( argc, argv, funcs );
	app->loop( 1000000 );



	return 0;
}







    // ===== SETUP GEOM
    //char* fname = "data/H_eFF.xyz";
    //char* fname = "data/e2_eFF_singlet.xyz";
    //char* fname = "data/e2_eFF_triplet.xyz";
    //char* fname = "data/H2_eFF.xyz";
    //char* fname = "data/He_eFF_singlet.xyz";
    //char* fname = "data/He_eFF_triplet.xyz";
    //char* fname = "data/H2O_eFF.xyz";
    //char* fname = "data/H2_eFF_spin.xyz";
    //char* fname = "data/Ce1_eFF.xyz";
    //char* fname = "data/Ce2_eFF.xyz";
    //char* fname = "data/Ce4_eFF.xyz";
    //char* fname = "data/CH3_eFF_spin.xyz";
    //char* fname = "data/CH4_eFF_flat_spin.xyz";
    //char* fname = "data/CH4_eFF_spin.xyz";
    //char* fname = "data/C2_eFF_spin.xyz";
    //char* fname = "data/C2H4_eFF_spin.xyz";
    //char* fname = "data/C2H4_eFF_spin_.xyz";
    //char* fname = "data/C2H6_eFF_spin.xyz";
    //char* fname = "data/C2H6_eFF_spin_.xyz";
    //ff.loadFromFile_xyz( "data/C2H4_eFF_spin.xyz" );
    //ff.loadFromFile_xyz( fname );

    //ff.loadFromFile_fgo( "data/e2_1g_2o_singlet.fgo" );
    //ff.loadFromFile_fgo( "data/e2_1g_2o_triplet.fgo );
    //ff.loadFromFile_fgo( "data/H_1g_1o.fgo" );
    //ff.loadFromFile_fgo( "data/He_singlet.fgo" );
    //ff.loadFromFile_fgo( "data/He_triplet.fgo" );
    //ff.loadFromFile_fgo( "data/H2_1g_2o.fgo" );
    //ff.loadFromFile_fgo( "data/H2.fgo" );
    //ff.loadFromFile_fgo( "data/C_1g.fgo" );
    //ff.loadFromFile_fgo( "data/C_2g_o1.fgo" );
    //ff.loadFromFile_fgo( "data/N2.fgo" );
    //ff.loadFromFile_fgo( "data/O2.fgo" );
    //ff.loadFromFile_fgo( "data/O2_half.fgo" );
    //ff.loadFromFile_fgo( "data/H2O_1g_8o.fgo" );

    //ff.loadFromFile_fgo( "data/C_e4_1g.fgo" );
    //ff.loadFromFile_fgo( "data/CH4.fgo" );
    //ff.loadFromFile_fgo( "data/NH3.fgo" );
    //ff.loadFromFile_fgo( "data/H2O.fgo" );
    //ff.loadFromFile_fgo( "data/C2H4.fgo" );    // Atoms Fly Away
    //ff.loadFromFile_fgo( "data/C2H2.fgo" );










