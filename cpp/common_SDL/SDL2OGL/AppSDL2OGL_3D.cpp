
#include <SDL2/SDL.h>


#include "Draw3D.h"

#include "AppSDL2OGL_3D.h" // THE HEADER


void AppSDL2OGL_3D::camera(){
    cam.zoom   = zoom;
    cam.aspect = ASPECT_RATIO;
    //Cam::ortho( cam, true );
    //Cam::perspective( cam );
    if (perspective){ Cam::perspective( cam ); }
    else            { Cam::ortho( cam, true ); }

}

void AppSDL2OGL_3D::draw(){
    opengl1renderer.clearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	opengl1renderer.enable    ( GL_LIGHTING );
	opengl1renderer.shadeModel( GL_FLAT     );

	Draw3D::drawBox       ( -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 0.8f, 0.8f, 0.8f );

	opengl1renderer.shadeModel( GL_SMOOTH     );
	Draw3D::drawSphere( Vec3f{3.0,3.0,3.0}, 1 );

	opengl1renderer.disable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

};

void AppSDL2OGL_3D::eventHandling ( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_o:  perspective   = !perspective; break;
                case SDLK_p:  first_person  = !first_person ;   break;
            }
            break;
        /*
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //mouseSpinning = true;
					//SDL_GetMouseState( &spinning_start.x, &spinning_start.y );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
					//mouseSpinning = false;
					//int mx,my;
					//SDL_GetMouseState( &mx, &my );
					//qCamera.fromTrackballQ( spinning_start.x, spinning_start.y, mx, my );
                    break;
            }
            break;
        */
    };
    AppSDL2OGL::eventHandling( event );
}

void AppSDL2OGL_3D::keyStateHandling( const Uint8 *keys ){

    if( keys[ SDL_SCANCODE_LEFT  ] ){ cam.qrot.dyaw  (  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ cam.qrot.dyaw  ( -keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ cam.qrot.dpitch(  keyRotSpeed ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ cam.qrot.dpitch( -keyRotSpeed ); }

	if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rotMat().a, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rotMat().a,  cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rotMat().b,  cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rotMat().b, -cameraMoveSpeed ); }
    if( keys[ SDL_SCANCODE_Q ] ){ cam.pos.add_mul( cam.rotMat().c, -cameraMoveSpeed ); }
	if( keys[ SDL_SCANCODE_E ] ){ cam.pos.add_mul( cam.rotMat().c,  cameraMoveSpeed ); }

/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.yaw  (  keyRotSpeed ); qCamera.normalize(); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.yaw  ( -keyRotSpeed ); qCamera.normalize(); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.pitch(  keyRotSpeed ); qCamera.normalize(); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.pitch( -keyRotSpeed ); qCamera.normalize(); }
 */

/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
*/
/*
    if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.yaw  (  0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.yaw  ( -0.01 ); printf( "yaw   qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.pitch(  0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.pitch( -0.01 ); printf( "pitch qCamera (%3.3f,%3.3f,%3.3f,%3.3f) \n", qCamera.x, qCamera.y, qCamera.z, qCamera.w ); }
*/

};


void AppSDL2OGL_3D::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        cam.qrot.qmul_T( q );
    }
    //qCamera.qmul( q );
}

void AppSDL2OGL_3D::drawCrosshair( float sz ){
    opengl1renderer.begin( GL_LINES );
    float whalf = WIDTH *0.5;
    float hhalf = HEIGHT*0.5;
    opengl1renderer.vertex3f( whalf-10,hhalf, 0 ); opengl1renderer.vertex3f( whalf+10,hhalf, 0 );
    opengl1renderer.vertex3f( whalf,hhalf-10, 0 ); opengl1renderer.vertex3f( whalf,hhalf+10, 0 );
    opengl1renderer.end();
}

void AppSDL2OGL_3D::drawMuseSelectionBox(){
    //opengl1renderer.lineWidth(3.0);
    //opengl1renderer.color3f(1.0,0.5,0.0); Draw3D::drawPointCross( ray0_start, 0.1 );    
    //opengl1renderer.lineWidth(1.0);
    Vec3f ray0_;        cam.rotMat().dot_to( (Vec3f)ray0, ray0_);
    Vec3f ray0_start_;  cam.rotMat().dot_to( (Vec3f)ray0_start, ray0_start_);
    opengl1renderer.color3f(1.0,0.5,0.0); Draw3D::drawTriclinicBoxT(cam.rotMat(), (Vec3f)ray0_start_, (Vec3f)ray0_ );   // Mouse Selection Box
    //opengl1renderer.color3f(0.0,0.5,1.0); Draw3D::drawTriclinicBox (cam.rotMat(), (Vec3f)ray0_start_, (Vec3f)ray0_ ); // Mouse Selection Box
}

AppSDL2OGL_3D::AppSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_, const char* name ) : AppSDL2OGL( id, WIDTH_, HEIGHT_, name ) {
	cam.qrot.setOne();
	cam.pos.set(0.0);
	GLbyte* s;
	// http://stackoverflow.com/questions/40444046/c-how-to-detect-graphics-card-model
	printf( "GL_VENDOR  : %s \n", opengl1renderer.getString(GL_VENDOR)  );
	printf( "GL_VERSION : %s \n", opengl1renderer.getString(GL_VERSION) );
}

