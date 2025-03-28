

//int verbosity = 0;
#include <globals.h>
//int idebug=0;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"
#include "IO_utils.h"

#include <SDL2/SDL.h>
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



// RenderToTexture:   https://nehe.gamedev.net/tutorial/radial_blur__rendering_to_a_texture/18004/


void RenderToTexture( GLuint tx ){
    opengl1renderer.viewport(0,0,128,128);                    // Set Our Viewport (Match Texture Size)
    opengl1renderer.clearColor(0.0f, 0.0f, 0.5f, 0.5);                // Set The Clear Color To Medium Blue
    opengl1renderer.clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     // Clear The Screen And Depth Buffer
    // RENDER 
    opengl1renderer.bindTexture(GL_TEXTURE_2D, tx );           // Bind To The Blur Texture
    opengl1renderer.copyTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 0, 0, 128, 128, 0);  // Copy Our ViewPort To The Blur Texture (From 0,0 To 128,128... No Border)
    //opengl1renderer.clearColor(0.0f, 0.0f, 0.5f, 0.5);                // Set The Clear Color To Medium Blue
    //opengl1renderer.clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     // Clear The Screen And Depth Buffer 
    //opengl1renderer.viewport(0 , 0,640 ,480);                 // Set Viewport (0,0 to 640x480)
}

void makeTexture( GLuint& tx, int sz ){
    opengl1renderer.genTextures(1, &tx);                       // Create 1 Texture
    opengl1renderer.bindTexture(GL_TEXTURE_2D, tx);            // Bind The Texture
    opengl1renderer.texImage2D(GL_TEXTURE_2D, 0, 4, sz, sz, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0 );           // Build Texture Using Information In data
    opengl1renderer.texParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
    opengl1renderer.texParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
}

void renderMolecule(int na, int nb, const Vec3d* atoms, const Vec2i* b2a ){
    opengl1renderer.begin(GL_LINES);
    for(int i=0; i<nb; i++){
        const Vec2i& b = b2a[i];
        Draw3D::vertex( atoms[b.a] );
        Draw3D::vertex( atoms[b.b] );
    }
    opengl1renderer.end();
}

void drawThumbnail( int itex, Vec2d p0, Vec2d p1, float sz ){
    opengl1renderer.enable     ( GL_TEXTURE_2D );
    opengl1renderer.bindTexture( GL_TEXTURE_2D, itex );
    opengl1renderer.enable(GL_BLEND);
    opengl1renderer.enable(GL_ALPHA_TEST);
    opengl1renderer.blendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    opengl1renderer.begin(GL_QUADS);
    opengl1renderer.texCoord2f( 0, 0 ); opengl1renderer.vertex3f( p0.x, p1.y, 0.0f );
    opengl1renderer.texCoord2f( 1, 0 ); opengl1renderer.vertex3f( p1.x, p1.y, 0.0f );
    opengl1renderer.texCoord2f( 1, 1 ); opengl1renderer.vertex3f( p1.x, p0.y, 0.0f );
    opengl1renderer.texCoord2f( 0, 1 ); opengl1renderer.vertex3f( p0.x, p0.y, 0.0f );
    opengl1renderer.end();
    opengl1renderer.disable  ( GL_TEXTURE_2D );
};

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
    std::vector<GLuint>    thumbnails;

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

    int texture_size = 256;
    int thumb_size   = texture_size;


    //GUI gui;
    //EditorGizmo  gizmo;

    // ---- Graphics objects
    int  fontTex,fontTex3D;
    int  ogl_mol=0;

    char str[256];


    int nrow = 0;
    int ncol = 0;

    // --------- Functions 

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	virtual void keyStateHandling( const Uint8 *keys );

    void MDloop();

	TestAppMolecularBrowser( int& id, int WIDTH_, int HEIGHT_ );

    void readMoleculess( bool bOrientFlat=true );
    void renderThumbnails( int i0, int n, float zoom_, bool bNew=true );

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

void TestAppMolecularBrowser::readMoleculess( bool bOrientFlat ){
    Mat3d rot;
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
            int ret = mol->loadByExt( path, 0 );             //printf("DEBUG 1 \n");
            if(ret>=0){
                mol->addToPos( mol->getCOG_minmax()*-1.0 );
                if(bOrientFlat){
                    mol->FindRotation( rot );
                    printf( "Rot: " ); printMat(rot);
                    //mol->orient( {0,0,0}, rot.c, rot.b );
                    mol->orient( {0,0,0}, rot.a, rot.b );
                }
                //mol->addToPos( mol->getCOG_av()*-1.0 );
                //mol->addToPos( Vec3d{-5.0,0.0,-5.0} );
                mol->assignREQs( params );         //printf("DEBUG 2 \n");
                if(mol->nbonds==0)mol->findBonds_brute( 0.5, true ); //printf("DEBUG 3 \n");
                molecules.push_back(mol);
                printf( "loaded mol[%i] '%s' natom %i nbond %i \n", i, fileNames[i].c_str(), mol->natoms, mol->natoms  );
            }
        //}catch (...){
            //printf( "cannot load '%s' \n", fileNames[%i].c_str() ); 
            //current_exception_diagnostic_information() << std::endl;
        //}
    }
}

void TestAppMolecularBrowser::renderThumbnails( int i0, int n, float zoom_, bool bNew){
    GLuint tx;
    //opengl1renderer.viewport(0,0,WIDTH,HEIGHT);                    // Set Our Viewport (Match Texture Size)
    opengl1renderer.viewport(0,0,texture_size,texture_size);
    _swap(zoom_,zoom);
    float ASPECT_RATIO_=ASPECT_RATIO; ASPECT_RATIO = 1;
    camera();
    for(int i=0;i<n;i++){
        if(bNew){
            makeTexture( tx, texture_size );
            thumbnails.push_back(tx);
        }
        opengl1renderer.clearColor(0.9f, 0.9f, 0.9f, 0.0f);                // Set The Clear Color To Medium Blue
        opengl1renderer.clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     // Clear The Screen And Depth Buffer

        Molecule& mol =  *molecules[i0+i];
        //renderMolecule( mol.natoms, mol.nbonds, mol.pos, mol.bond2atom );
        //Draw3D::atomsREQ( mol.natoms, mol.pos, mol.REQs, ogl_sph, qsc=1, );
        Draw3D::atoms( mol.natoms, mol.pos, mol.atomType, params, 1.0, 0.5, 1.0 );
        opengl1renderer.color3f(0.0f,0.0f,0.0f);
        Draw3D::bonds( mol.nbonds, mol.bond2atom, mol.pos);
        Draw3D::drawText( fileNames[i0+i].c_str(), Vec3d{-9.0,+9.0,5.0}, fontTex, 0.07, 0);


        opengl1renderer.color3f(1.0f,1.0f,1.0f);
        opengl1renderer.bindTexture(GL_TEXTURE_2D, tx );           // Bind To The Blur Texture
        //opengl1renderer.copyTexImage2D(GL_TEXTURE_2D, 0, GL_LUMINANCE, 0, 0, texture_size, texture_size, 0);  // Copy Our ViewPort To The Blur Texture (From 0,0 To 128,128... No Border)
        opengl1renderer.copyTexImage2D  (GL_TEXTURE_2D, 0, GL_RGBA,      0, 0, texture_size, texture_size, 0);

    }
    _swap(zoom_,zoom);
    _swap(ASPECT_RATIO_,ASPECT_RATIO);
    opengl1renderer.clearColor(0.5f, 0.5f, 0.5f, 0.5);                // Set The Clear Color To Medium Blue
    opengl1renderer.clear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);     // Clear The Screen And Depth Buffer 
    opengl1renderer.viewport(0 , 0, WIDTH, HEIGHT );                 // Set Viewport (0,0 to 640x480)
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

    ncol=(int)( WIDTH /thumb_size );
    nrow=(int)( HEIGHT/thumb_size );

    /*
    GLuint tx;
    makeTexture( tx, texture_size );
    thumbnails.push_back(tx);
    */

    picked_lvec = &builder.lvec.a;

    // ---- Graphics setup
    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //opengl1renderer.lightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//opengl1renderer.lightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    initGUI();

}

//=================================================
//                   DRAW()
//=================================================

void TestAppMolecularBrowser::draw(){
    //opengl1renderer.clearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    opengl1renderer.clearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    //opengl1renderer.enable(GL_LINE_SMOOTH);
    opengl1renderer.enable(GL_BLEND);
    opengl1renderer.enable(GL_LIGHTING );
    opengl1renderer.enable(GL_DEPTH_TEST);

    if(frameCount==1){ 
        cam.qrot.pitch( M_PI );  qCamera0=cam.qrot; 
    }
    if(frameCount==2){
        renderThumbnails( 0, molecules.size(), 10.0, true );
    }
    //drawThumbnail( thumbnails[0], {-5.0,-5.0}, {5.0,5.0}, texture_size );
};

void TestAppMolecularBrowser::drawHUD(){
    opengl1renderer.disable ( GL_LIGHTING );
    //gui.draw();

    float dx=(thumb_size+1);
    int imol=0;
    for(int iy=0; iy<nrow; iy++){
        for(int ix=0; ix<ncol; ix++){
            if(imol<thumbnails.size())drawThumbnail( thumbnails[imol], {ix*dx,iy*dx+(dx-1)}, (Vec2d){ix*dx+(dx-1),iy*dx}, texture_size );
            imol++;
        }
    }
}

void TestAppMolecularBrowser::saveScreenshot( int i, const char* fname ){
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
        case SDL_WINDOWEVENT:{
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE:{ quit(); } break;
            } break;
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
    //if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rotMat().a, -cameraMoveSpeed ); }
	//if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rotMat().a,  cameraMoveSpeed ); }
    //if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rotMat().b,  cameraMoveSpeed ); }
	//if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rotMat().b, -cameraMoveSpeed ); }
    //if( keys[ SDL_SCANCODE_Q ] ){ cam.pos.add_mul( cam.rotMat().c, -cameraMoveSpeed ); }
	//if( keys[ SDL_SCANCODE_E ] ){ cam.pos.add_mul( cam.rotMat().c,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_LEFT  ] ){ cam.pos.add_mul( cam.rotMat().a, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ cam.pos.add_mul( cam.rotMat().a,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_UP    ] ){ cam.pos.add_mul( cam.rotMat().b,  cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ cam.pos.add_mul( cam.rotMat().b, -cameraMoveSpeed ); }
    //AppSDL2OGL_3D::keyStateHandling( keys );
};

//void command_example(double x, void* caller){
//    TestAppMolecularBrowser* app = (TestAppMolecularBrowser*)caller;
//    app->zoom = x;
//};

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
















