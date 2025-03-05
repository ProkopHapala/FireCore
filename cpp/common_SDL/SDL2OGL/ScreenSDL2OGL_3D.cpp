
#include <SDL2/SDL.h>


#include "ScreenSDL2OGL_3D.h"
#include "quaternion.h"

void ScreenSDL2OGL_3D::camera(){
    cam.qrot = Quat4fIdentity;
    cam.zoom   = zoom;
    cam.aspect = ASPECT_RATIO;
    //Cam::ortho( cam, true );
    //Cam::perspective( cam );
    if (perspective){ Cam::perspective( cam ); }
    else            { Cam::ortho( cam, true ); }

}

void ScreenSDL2OGL_3D::draw   (){
    opengl1renderer.clearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	opengl1renderer.enable    ( GL_LIGHTING );
	opengl1renderer.shadeModel( GL_FLAT     );

	Draw3D::drawBox       ( -1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f, 0.8f, 0.8f, 0.8f );

	opengl1renderer.shadeModel( GL_SMOOTH     );
	Draw3D::drawSphere_oct( 5, 1.0f, Vec3f{3.0,3.0,3.0} );

	opengl1renderer.disable ( GL_LIGHTING );
	Draw3D::drawAxis ( 3.0f );

};

void ScreenSDL2OGL_3D::eventHandling ( const SDL_Event& event  ){
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
    ScreenSDL2OGL::eventHandling( event );
}

void ScreenSDL2OGL_3D::keyStateHandling( const Uint8 *keys ){

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

void ScreenSDL2OGL_3D::mouseHandling( ){
    SDL_GetMouseState( &mouseX, &mouseY );   mouseY=HEIGHT-mouseY;
    mouse_begin_x = (2*mouseX-WIDTH )*zoom/HEIGHT;
    mouse_begin_y = (2*mouseY-HEIGHT)*zoom/HEIGHT;
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    //printf( " %i %i \n", mx,my );
    if ( buttons & SDL_BUTTON(SDL_BUTTON_RIGHT)) {
        Quat4f q; q.fromTrackball( 0, 0, -mx*mouseRotSpeed, my*mouseRotSpeed );
        cam.qrot.qmul_T( q );
        printf("creenSDL2OGL_3D::mouseHandling( )\n");
    }
    //qCamera.qmul( q );
}

ScreenSDL2OGL_3D::ScreenSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) {
	cam.qrot = Quat4fIdentity;
	cam.pos.set(0.0);
}
