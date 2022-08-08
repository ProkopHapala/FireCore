
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

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
//#include "MMFFmini.h"
//#include "NBFF.h"
//#include "GridFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
//#include "DynamicOpt.h"
//#include "QEq.h"

#include "browser.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams
//#include "MolecularDraw.h"
#include "GUI.h"
#include "EditorGizmo.h"
//#include "FireCoreAPI.h"
//#include "NBSRFF.h"
//#include "MarchingCubes.h"
//#include "SimplexRuler.h"
#include "AppSDL2OGL_3D.h"

int idebug=0;

Quat4f qFront {-M_SQRT1_2,      0.0f,     0.0f,M_SQRT1_2 };
Quat4f qBack  {      0.0f, M_SQRT1_2,M_SQRT1_2,      0.0f };
Quat4f qTop   { 1.0f, 0.0f, 0.0f, 0.0f};
Quat4f qBottom{ 0.0f, 0.0f, 0.0f, 1.0f};
Quat4f qLeft  {-0.5f, 0.5f, 0.5f, 0.5f};
Quat4f qRight { 0.5f, 0.5f, 0.5f,-0.5f};

inline bool file_exist(const char* fname) { if (FILE *file = fopen(fname, "r")) { fclose(file); return true; } else { return false; } }

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

//void command_example(double x, void* caller);

class TestAppMolecularBrowser : public AppSDL2OGL_3D, public Browser { public:
    // ---- Simulations objects
	//Molecule    mol;
	MMFFparams  params;
    MM::Builder builder;

    std::vector<Molecule*> molecules;

    int* atypes = 0;
    int* atypeZ = 0;

    bool bNonBonded = true;

    bool   bDragging = false;
    //Vec3d  ray0_start;
    //Vec3d  ray0;
    //int ipicked    = -1; // picket atom 
    Vec3d* picked_lvec = 0;
    int perFrame =  1;
    Quat4f qCamera0;


    //GUI gui;
    //EditorGizmo  gizmo;

    // ---- Graphics objects
    int  fontTex,fontTex3D;
    int  ogl_sph=0;
    int  ogl_mol=0;

    char str[256];

    // --------- Functions 

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

    void MDloop();

	TestAppMolecularBrowser( int& id, int WIDTH_, int HEIGHT_ );

    void readMoleculess();

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );
    void loadGeom();
    void initGUI();
    void drawingHex(double z0);
};

//=================================================
//                   INIT()
//=================================================

void TestAppMolecularBrowser::initGUI(){

}

void TestAppMolecularBrowser::readMoleculess(){
    for(int i=0; i<molecules.size(); i++){
        delete molecules[i];
    }
    molecules.clear();
    for(int i=0; i<fileNames.size(); i++){
        printf("==================\n", i );
        printf("readMoleculess[%i]\n", i );
            Molecule* mol = new Molecule();
            std::string path = work_dir+"/"+fileNames[i]; //printf("DEBUG 0 \n");
            mol->bindParams( &params );
            int ret = mol->loadByExt( path, 2 );             //printf("DEBUG 1 \n");
            if(ret>=0){
                mol->assignREQs( params );         //printf("DEBUG 2 \n");
                mol->findBonds_brute( 0.9, true ); //printf("DEBUG 3 \n");
                molecules.push_back(mol);
                printf( "loaded mol[%i] '%s' natom %i nbond %i \n", i, fileNames[i].c_str(), mol->natoms, mol->natoms  );
            }
        //}catch (...){
            //printf( "cannot load '%s' \n", fileNames[%i].c_str() ); 
            //current_exception_diagnostic_information() << std::endl;
        //}
    }
}

TestAppMolecularBrowser::TestAppMolecularBrowser( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ),Browser::Browser( "./common_resources" ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ---- Load Atomic Type Parameters
    int nheavy = 0;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);

    //extensions.insert( {"lvs",1} );
    extensions.insert( {"xyz",1} );
    extensions.insert( {"mol",1} );

    readDir( work_dir );
    for(int i=0; i<fileNames  .size(); i++ ){ printf("file[%i] '%s'\n",i,fileNames  [i].c_str()); }
    for(int i=0; i<subDirNames.size(); i++ ){ printf("dir [%i] '%s'\n",i,subDirNames[i].c_str()); }

    readMoleculess();


    picked_lvec = &builder.lvec.a;

    // ---- Graphics setup
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );
    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    initGUI();

}

//=================================================
//                   DRAW()
//=================================================

void TestAppMolecularBrowser::draw(){
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    //glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING );
    glEnable(GL_DEPTH_TEST);

    if(frameCount==1){ qCamera.pitch( M_PI );  qCamera0=qCamera; }
    

};

void TestAppMolecularBrowser::drawHUD(){
    glDisable ( GL_LIGHTING );
    //gui.draw();
}

void TestAppMolecularBrowser::saveScreenshot( int i, const char* fname ){
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

void TestAppMolecularBrowser::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    //gui.onEvent( mouseX, mouseY, event );
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
            } break;
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    break;
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    break;
                case SDL_BUTTON_RIGHT:
                    break;
            }
            break;
        case SDL_WINDOWEVENT:
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE:
                    quit();
                    break;
            } break;
    };
    AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void  TestAppMolecularBrowser::keyStateHandling( const Uint8 *keys ){
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
//    TestAppMolecularBrowser* app = (TestAppMolecularBrowser*)caller;
//    app->zoom = x;
//};


/*
void TestAppMolecularBrowser::loadGeom(){
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
*/

// ===================== MAIN

TestAppMolecularBrowser * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;

    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);
	thisApp = new TestAppMolecularBrowser( junk, DM.w-100, DM.h-100 );
	//thisApp = new TestAppMolecularBrowser( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















