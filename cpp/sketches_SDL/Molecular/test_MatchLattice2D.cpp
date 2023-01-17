#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include  <functional>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"



#include "LatticeMatch2D.h"




class TestAppLatticeMatch2D: public AppSDL2OGL_3D { public:

    LatticeMatch2D LM;
    int      fontTex;
    int ogl1=-1;

    

    // ---- Functions

    virtual void draw   ();
    //virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppLatticeMatch2D( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppLatticeMatch2D::TestAppLatticeMatch2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_, " test_eFF " ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    //plot1.fontTex=fontTex;


    LM.lat0[0].set(  1.3,  0.0 );
    LM.lat0[1].set( -0.5,  0.6 );
    LM.lat1[0].set(  0.3,  0.2 );
    LM.lat1[1].set( -0.13, 0.9 );

    ogl1 = glGenLists (1);

    double Rmax = 10.0;

    glNewList (ogl1, GL_COMPILE);
    Vec2d p;
    p=LM.lat0[0]; glColor3f(1.0,0.0,0.0); Draw3D::drawLine(  Vec3dZero,  {p.x,p.y,0.0}  );
    p=LM.lat0[1]; glColor3f(0.0,0.0,1.0); Draw3D::drawLine(  Vec3dZero,  {p.x,p.y,0.0}  );
    p=LM.lat1[0]; glColor3f(0.0,0.5,0.0); Draw3D::drawLine(  Vec3dZero,  {p.x,p.y,0.0}  );
    p=LM.lat1[1]; glColor3f(0.0,0.5,1.0); Draw3D::drawLine(  Vec3dZero,  {p.x,p.y,0.0}  );
    Draw3D::drawCircleAxis(100,Vec3dZero, Vec3dX, Vec3dZ, Rmax  );
    //LM.walk2D( Rmax, 0.05 );
    LM.walk2D( Rmax, 0.02 );
    glEndList();

}

void TestAppLatticeMatch2D::draw(){
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    glCallList(ogl1);
};

/*
void TestAppLatticeMatch2D::drawHUD(){
	//glTranslatef( 100.0, 250.0, 0.0 );
	//glScalef    ( 100.0, 100.0, 1.0 );
	//plot1.view();
}
*/


void TestAppLatticeMatch2D::eventHandling ( const SDL_Event& event  ){

/*
    Vec3d pa0;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;

                case SDLK_i: ff.info(); break;
                case SDLK_LEFTBRACKET :  Espread *= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_RIGHTBRACKET:  Espread /= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                case SDLK_e: bMapElectron=!bMapElectron; break;
                case SDLK_f:{
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
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    */
    AppSDL2OGL::eventHandling( event );
}


// ===================== MAIN

TestAppLatticeMatch2D* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppLatticeMatch2D( junk , 800, 600 );
	thisApp->loop( 1000000 );

    //char str[40];  sprintf(str,  );
	//SDL_SetWindowTitle( thisApp->child_windows[1]->window, " test_eFF " );

	return 0;
}
















