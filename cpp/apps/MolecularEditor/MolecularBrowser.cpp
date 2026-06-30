/// @file MolecularBrowser.cpp
/// @brief Entry point and orchestrator for the molecular file browser application.
///
/// MolecularBrowser is an ACDsee-like browser for molecular files (.xyz, .mol, .mol2).
/// It switches between two modes:
/// - BROWSE: thumbnail grid view (delegates to BrowserView)
/// - VIEW:   interactive 3D molecule view (delegates to MolView)
/// Press Enter to toggle between modes. Click thumbnails to view molecules.
/// Click subdirectory buttons or press Backspace to navigate directories.
///
/// Role in repo: Thin app/entry point that owns the SDL window, GL context,
/// shared resources (MMFFparams, fonts, sphere display list), and delegates
/// rendering/input to BrowserView and MolView components.

#include <globals.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <string>
#include <libgen.h>
#include <unistd.h>

#include "testUtils.h"
#include "IO_utils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "Draw2D.h"
#include "Draw.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFparams.h"

#include "Draw3D_Molecular.h"
#include "GUI.h"
#include "AppSDL2OGL_3D.h"

#include "argparse.h"

#include "BrowserView.h"
#include "MolView.h"

// ---- Path configuration (can be overridden via CLI args -res, -dir, -ini)
std::string g_res_dir  = "common_resources";
std::string g_work_dir = ".";

void loadIniFile( const char* fname ){
    FILE* f = fopen(fname, "r");
    if(!f){ printf("ERROR: cannot open ini file '%s'\n", fname); return; }
    char line[512];
    while(fgets(line, sizeof(line), f)){
        if(line[0]=='#' || line[0]=='\n' || line[0]=='\0' || line[0]==' ') continue;
        char key[256], val[256];
        if(sscanf(line, "%255[^=]=%255[^\n]", key, val)==2){
            char* p=key+strlen(key)-1; while(p>=key && (*p==' '||*p=='\t')) *p--='\0';
            char* v=val; while(*v==' '||*v=='\t') v++;
            if     (strcmp(key,"res_dir")==0)  g_res_dir  = v;
            else if(strcmp(key,"work_dir")==0) g_work_dir = v;
            else printf("WARNING: unknown ini key '%s'\n", key);
        }
    }
    fclose(f);
    printf("Ini loaded: res_dir='%s' work_dir='%s'\n", g_res_dir.c_str(), g_work_dir.c_str());
}

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppMolecularBrowser : public AppSDL2OGL_3D { public:
    // ---- Shared resources
    MMFFparams  params;
    int  fontTex, fontTex3D;
    int  ogl_sph = 0;

    // ---- Mode
    enum Mode { BROWSE, VIEW };
    Mode mode = BROWSE;

    // ---- Components
    BrowserView browser;
    MolView     molView;

    Quat4f qCamera0;

    // --------- Functions

    virtual void draw   ();
    virtual void drawHUD();
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppMolecularBrowser( int& id, int WIDTH_, int HEIGHT_ );
    ~TestAppMolecularBrowser(){}
};

//=================================================
//                   INIT()
//=================================================

TestAppMolecularBrowser::TestAppMolecularBrowser( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ), browser( g_work_dir ) {

    fontTex   = makeTextureHard( (char*)(g_res_dir+"/dejvu_sans_mono_RGBA_pix.bmp").c_str() ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( (char*)(g_res_dir+"/dejvu_sans_mono_RGBA_inv.bmp").c_str() );

    // ---- Load Atomic Type Parameters
    params.loadElementTypes( (g_res_dir+"/ElementTypes.dat").c_str() );
    params.loadAtomTypes   ( (g_res_dir+"/AtomTypes.dat").c_str() );
    params.loadBondTypes   ( (g_res_dir+"/BondTypes.dat").c_str() );

    // ---- Graphics setup
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );
    float l_specular[]{ 0.0f, 0.0f, 0.0f, 1.0f };
    glLightfv( GL_LIGHT0, GL_SPECULAR, l_specular );

    // ---- Init components
    browser.init( &params, fontTex, ogl_sph );
    molView .init( &params, fontTex, fontTex3D, ogl_sph );

    // ---- Initial directory read
    browser.readDir( browser.work_dir );
    browser.readMoleculess();
}

//=================================================
//                   DRAW()
//=================================================

void TestAppMolecularBrowser::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_BLEND);
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);

    if(frameCount==1){ qCamera.pitch( M_PI ); qCamera0=qCamera; }

    if(mode == BROWSE){
        if(browser.bNeedRerender){
            browser.bNeedRerender = false;
            browser.renderThumbnails( this );
        }
    }else{ // VIEW
        molView.draw3D( this );
    }
}

void TestAppMolecularBrowser::drawHUD(){
    if(mode == BROWSE){
        browser.drawHUD( this );
    }else{ // VIEW
        molView.drawHUD( this );
    }
}

void TestAppMolecularBrowser::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN:
            switch( event.key.keysym.sym ){
                case SDLK_ESCAPE:   quit(); break;
                case SDLK_RETURN:
                case SDLK_KP_ENTER:
                    if(mode == BROWSE){
                        if(browser.selectedThumb >= 0 && browser.selectedThumb < (int)browser.molecules.size()){
                            molView.setMolecule( browser.molecules[browser.selectedThumb], browser.fileNames[browser.selectedThumb], browser.selectedThumb );
                        }
                        mode = VIEW;
                        molView.bNeedFit = true;
                        printf("Switched to VIEW mode\n");
                    }else{
                        mode = BROWSE;
                        printf("Switched to BROWSE mode\n");
                    }
                    break;
                case SDLK_BACKSPACE:
                    if(mode == BROWSE) browser.navigateToDir("..");
                    break;
                // ---- Numpad view presets (VIEW mode) ----
                case SDLK_KP_7: if(mode == VIEW){ qCamera = Quat4fBack;   } break;  // back
                case SDLK_KP_9: if(mode == VIEW){ qCamera = Quat4fFront;  } break;  // front
                case SDLK_KP_8: if(mode == VIEW){ qCamera = Quat4fTop;    } break;  // top
                case SDLK_KP_2: if(mode == VIEW){ qCamera = Quat4fBotton; } break;  // bottom
                case SDLK_KP_4: if(mode == VIEW){ qCamera = Quat4fLeft;   } break;  // left
                case SDLK_KP_6: if(mode == VIEW){ qCamera = Quat4fRight;  } break;  // right
                case SDLK_KP_5: if(mode == VIEW){ qCamera = Quat4fIdentity;} break; // reset
                // ---- Label toggles (VIEW mode) ----
                case SDLK_l: if(mode == VIEW){ molView.bViewAtomLabels  ^= 1; } break;
                case SDLK_t: if(mode == VIEW){ molView.bViewAtomTypes   ^= 1; } break;
                case SDLK_b: if(mode == VIEW){ molView.bViewBondLabels  ^= 1; } break;
                case SDLK_i: if(mode == VIEW){ molView.bViewBondLenghts ^= 1; } break;
                case SDLK_LEFT:
                    if(mode == BROWSE) browser.moveCursor(-1, 0);
                    break;
                case SDLK_RIGHT:
                    if(mode == BROWSE) browser.moveCursor( 1, 0);
                    break;
                case SDLK_UP:
                    if(mode == BROWSE) browser.moveCursor( 0,-1);
                    break;
                case SDLK_DOWN:
                    if(mode == BROWSE) browser.moveCursor( 0, 1);
                    break;
            } break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if(mode == BROWSE){
                        int clicked = browser.handleClick( event.button.x, event.button.y, HEIGHT );
                        if(clicked >= 0 && clicked < (int)browser.molecules.size()){
                            molView.setMolecule( browser.molecules[clicked], browser.fileNames[clicked], clicked );
                            mode = VIEW;
                            printf("Opened mol[%i] '%s' in VIEW mode\n", clicked, browser.fileNames[clicked].c_str());
                        }
                    }else if(mode == VIEW){
                        molView.handleMouse( mouseX, mouseY, event );
                    }
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
}

void TestAppMolecularBrowser::keyStateHandling( const Uint8 *keys ){
    if(mode == VIEW){ AppSDL2OGL_3D::keyStateHandling( keys ); }
};

// ===================== MAIN

TestAppMolecularBrowser * thisApp;

int main(int argc, char *argv[]){
    SDL_Init(SDL_INIT_VIDEO);
    SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    int junk;

    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);

    // ---- Parse CLI arguments
    LambdaDict funcs;
    funcs["-res"]={1,[&](const char** ss){ g_res_dir  = ss[0]; printf("ARG -res  '%s'\n", g_res_dir.c_str() );  }};
    funcs["-dir"]={1,[&](const char** ss){ g_work_dir = ss[0]; printf("ARG -dir  '%s'\n", g_work_dir.c_str()); }};
    funcs["-ini"]={1,[&](const char** ss){ printf("ARG -ini  '%s'\n", ss[0]); loadIniFile(ss[0]); }};
    process_args( argc, argv, funcs, false );

    printf("MolecularBrowser: res_dir='%s' work_dir='%s'\n", g_res_dir.c_str(), g_work_dir.c_str());

    thisApp = new TestAppMolecularBrowser( junk, DM.w-100, DM.h-100 );
    thisApp->loop( 1000000 );
    return 0;
}
















