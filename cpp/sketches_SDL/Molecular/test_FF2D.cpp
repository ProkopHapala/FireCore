
/*

NOTE:

WORKING COMMITS - ( Relax CH4 and C4H12 )
commit d8b1156f6ef8f283785dae1fee71105d591a3280    2021-Apr-26  testing CLCFGO.h vs eFF.h for Hydrogen atom with electron radius 0.5A…

NOT WORKING COMMITS
commit dae1d0b16d3b892f0c3d982c5f277df38bfb4179    2021-Jun-01    CLCFGO : option to out-project force components which breaks normaliz…
commit 94a94e956acad8e3d23a54acbd0f715fe0d1f827    2021-May-05    CLCFGO : tested H2 and H2O molecule with respect to eFF; total energy…


*/


int verbosity = 0;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include  <functional>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw2D.h"
//#include "Draw3D.h"
//#include "Solids.h"

//#include "fastmath.h"
//#include "Vec3.h"
//#include "Mat3.h"
//#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
//#include "Plot2D.h"
//#include "DynamicOpt.h"

#include "FF2D.h"


class TestAppFF2D: public AppSDL2OGL_3D { public:

    bool bRun = false;
    FF2D  ff;
    //DynamicOpt opt;

    int perFrame = 1;


    int ipicked  = 0;

    //GLint oglSph = 0;

    //Plot2D plot1;

    int      fontTex;

    // ---- Functions

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppFF2D( int& id, int WIDTH_, int HEIGHT_ );
    void init2DMap( int n, double dx );

};

TestAppFF2D::TestAppFF2D( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_, " test_eFF " ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    //plot1.fontTex=fontTex;
    
    //ff.insertString( "6 0 2 3 4; 6 0 1; 6 0 1; 6 0 1;" );  // neo-butan

    ff.insertString( 
    R"(
     6 1 6 2 7; 
     6 1 1 3; 
     6 1 2 4 8; 
     6 1 3 5;
     6 1 4 6 9; 
     6 1 5 1;
     6 0 1;
     6 0 3;
     6 0 5;
     )" );  // Benzene

    ff.try_realloc();
    ff.cleanVelocity();  
    srand(454);    
    ff.setPosRandom(5.0);    
    //return ff.atoms.size();  

}

void TestAppFF2D::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    perFrame = 1;
    if(bRun){
        double F2 = 1.0;
        double Etot;
        for(int itr=0;itr<perFrame;itr++){
            F2= ff.step(0.1,0.1);
        }
        if( F2 < 1e-4 ){
            //printf( "Finished: E %g | Ek %g Eee %g EeePaul %g Eaa %g Eae %g EaePaul %g \n", Etot, ff.Ek, ff.Eee, ff.EeePaul, ff.Eaa, ff.Eae, ff.EaePaul );
            bRun=false;
        }
    }

    for(int ia=0; ia<ff.atoms.size(); ia++){
        const Atom2D& A = ff.atoms[ia];
        Draw2D::drawPointCross((Vec2f)A.pos,0.1);
        for(int j=0; j<A.nb; j++ ){
            int ja=A.neighs[j];
            Draw2D::drawLine( (Vec2f)A.pos, (Vec2f)ff.atoms[ja].pos );
        }
    }

};

void TestAppFF2D::drawHUD(){
	//glTranslatef( 100.0, 250.0, 0.0 );
	//glScalef    ( 100.0, 100.0, 1.0 );
	//plot1.view();
    glTranslatef( 10.0,HEIGHT-20.0,0.0 );
	glColor3f(0.5,0.0,0.3);
}

void TestAppFF2D::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    Vec3d pa0;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_SPACE: bRun = !bRun; break;
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
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppFF2D* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppFF2D( junk , 800, 600 );
	thisApp->loop( 1000000 );

    //char str[40];  sprintf(str,  );
	//SDL_SetWindowTitle( thisApp->child_windows[1]->window, " test_eFF " );

	return 0;
}
















