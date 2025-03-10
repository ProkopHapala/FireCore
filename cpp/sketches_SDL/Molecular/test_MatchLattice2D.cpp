#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include  <functional>

#include <SDL2/SDL.h>

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
    int ipick=-1;

    

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

/*
    LM.lat0[0].set(  1.3,  0.0 );
    LM.lat0[1].set( -0.5,  0.6 );
    //LM.lat1[0].set(  0.3,  0.5 );
    //LM.lat1[1].set(  0.5, -0.3 );
    //LM.lat1[0].set(  0.5, -0.3 );
    //LM.lat1[1].set(  0.3,  0.5 );
*/  

    ogl1 = opengl1renderer.genLists (1);

    double Rmax = 10.0;
    double dmax = 0.1;
    double dang = 0.1;

    opengl1renderer.newList (ogl1, GL_COMPILE);
    Vec2d p;
    // p=LM.lat0[0]*5; opengl1renderer.color3f(1.0,0.0,0.0); Draw3D::drawLine(  Vec3dZero,  {p.x,p.y,0.0}  );
    // p=LM.lat0[1]*5; opengl1renderer.color3f(0.0,0.0,1.0); Draw3D::drawLine(  Vec3dZero,  {p.x,p.y,0.0}  );
    p=LM.lat1[0]*5; Draw3D::drawLine( renderer, Vec3dZero,  {p.x,p.y,0.0}, {0, 0.5, 0} );
    p=LM.lat1[1]*5; Draw3D::drawLine( renderer, Vec3dZero,  {p.x,p.y,0.0}, {0, 0.5, 1} );
    opengl1renderer.color3f(0.0,0.0,0.0); Draw3D::drawCircleAxis(100,Vec3dZero, Vec3dX, Vec3dZ, Rmax  );
    //LM.walk2D( Rmax, 0.05 );
    LM.walk2D( Rmax, dmax );
    opengl1renderer.endList();

    LM.angleToRange();
    LM.sort();
    
    //for(int i=0; i<LM.match_u.size(); i++){ printf( "match_u[%i] ang %g n,d(%i,%g) \n", i, LM.match_u[i].alpha, LM.match_u[i].n, LM.match_u[i].d ); }
    //for(int i=0; i<LM.match_v.size(); i++){ printf( "match_v[%i] ang %g n,d(%i,%g) \n", i, LM.match_v[i].alpha, LM.match_v[i].n, LM.match_v[i].d ); }

    int n = LM.matchAngles(dang);

    printf( "matchAngles: nfound = %i\n", n );
    for(int i=0; i<n; i++){
        Vec2i m = LM.matches[i];
        const Latmiss& Lu = LM.match_u[m.i];
        const Latmiss& Lv = LM.match_v[m.j];
        printf( "[%i][%i %i %i %i]\n", i, Lu.ia, Lu.ib, Lv.ia, Lv.ib );
    }
}

void TestAppLatticeMatch2D::draw(){
    opengl1renderer.clearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    opengl1renderer.enable(GL_DEPTH_TEST);
    opengl1renderer.disable(GL_LIGHTING);

    opengl1renderer.callList(ogl1);

    if(ipick>=0){
        
        Vec2d u,v;
        
        LM.makeMatch( LM.matches[ipick],u,v,true );
        Draw3D::drawLine( renderer, Vec3dZero,  {u.x,u.y,0.0}, {1.0,0.0,0.0} );
        Draw3D::drawLine( renderer, Vec3dZero,  {v.x,v.y,0.0}, {0.0,0.0,1.0} );

        LM.makeMatch( LM.matches[ipick],u,v,false );
        Draw3D::drawLine( renderer, Vec3dZero,  {u.x,u.y,0.0}, {1.0,0.7,0.0} );
        Draw3D::drawLine( renderer, Vec3dZero,  {v.x,v.y,0.0}, {0.0,0.7,1.0} );
        
        //v = LM.reproduce_grid( LM.match_u[ipick]                      );   opengl1renderer.color3f(0.0,0.0,1.0);  Draw3D::drawLine( Vec3dZero,  {v.x,v.y,0.0} );
        //v = LM.reproduce_vec ( LM.match_u[ipick], LM.lat1[0], LM.angU );   opengl1renderer.color3f(1.0,0.7,0.0);  Draw3D::drawLine( Vec3dZero,  {v.x,v.y,0.0} );
        //v = LM.reproduce_grid( LM.match_v[ipick]                      );   opengl1renderer.color3f(0.0,0.0,1.0);  Draw3D::drawLine( Vec3dZero,  {v.x,v.y,0.0} );
        //v = LM.reproduce_vec ( LM.match_v[ipick], LM.lat1[1], LM.angV+LM.angUV );   opengl1renderer.color3f(1.0,0.7,0.0);  Draw3D::drawLine( Vec3dZero,  {v.x,v.y,0.0} );
    }
};

/*
void TestAppLatticeMatch2D::drawHUD(){
	//opengl1renderer.translatef( 100.0, 250.0, 0.0 );
	//opengl1renderer.scalef    ( 100.0, 100.0, 1.0 );
	//plot1.view();
}
*/


void TestAppLatticeMatch2D::eventHandling ( const SDL_Event& event  ){

    int npick = LM.matches.size();
    //int npick = LM.match_u.size();
    //int npick = LM.match_v.size();
    Vec3d pa0;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_LEFTBRACKET :  ipick--; if(ipick<0     )ipick=npick-1; break;
                case SDLK_RIGHTBRACKET:  ipick++; if(ipick>=npick)ipick=0;       break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rotMat().c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rotMat().c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
        */
            break;
    };
    
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
















