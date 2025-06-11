
#include "ScreenSDL2OGL.h" // THE HEADER
#include "Renderer.h"
#include <SDL2/SDL_video.h>
#include <iostream>
#include <chrono>

//#include "testUtils.h"

// ============== per frame

void setLightingNormal(){
	opengl1renderer.enable     ( GL_DEPTH_TEST       );
	opengl1renderer.polygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

void setLightingRGB(){
	opengl1renderer.enable     ( GL_DEPTH_TEST       );
	opengl1renderer.polygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

std::chrono::duration<long, std::ratio<1, 1000>> durationTotal;
int durationCount = 0;

void ScreenSDL2OGL::update( ){
	auto start = std::chrono::high_resolution_clock::now();

	if( GL_LOCK ){ printf("ScreenSDL2OGL::update GL_LOCK\n"); return; }
	GL_LOCK = true;
	//printf( " window[%i] SDL_GL_MakeCurrent \n", id );
    SDL_GL_MakeCurrent(window, glctx);
	SDL_GL_GetDrawableSize(window, &WIDTH, &HEIGHT);
	GLES::screen_size = { WIDTH, HEIGHT };
	glViewport(0, 0, GLES::screen_size.x, GLES::screen_size.y);
	glScissor (0, 0, GLES::screen_size.x, GLES::screen_size.y);
	if (GLES::active_camera) GLES::active_camera->update();
	camera();

	glClearColor(0.2f, 0.5f, 0.8f, 1.0f);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	draw();
	cameraHUD();
	drawHUD();
	
	//opengl1renderer.popMatrix();
	//opengl1renderer.flush();
	//SDL_RenderPresent( );
	frameCount++;
    SDL_GL_SwapWindow(window);
    //printf( " window[%i] SDL_GL_SwapWindow \n", id );
    GL_LOCK = false;

	auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::milliseconds>(stop - start);

	durationTotal += duration;
	durationCount++;

	if (durationCount == 100){
		auto ms = durationTotal.count()/(float)durationCount;
		std::cout << "update() took " << ms << " milliseconds (last "<<durationCount<<" avg) -> fps = "<<1000.0/ms<<"\n";
		durationCount = 0;
		durationTotal = std::chrono::duration<long, std::ratio<1, 1000>>::zero();
	}
};

void ScreenSDL2OGL::draw(){}; // virtual function, meant to be overriden

void ScreenSDL2OGL::drawHUD(){ };

//void ScreenSDL2OGL::inputHanding(){};

void ScreenSDL2OGL::eventHandling( const SDL_Event& event ){
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = true;  break;
                case SDL_BUTTON_RIGHT: RMB = true;  break;
            };  break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:  LMB = false; break;
                case SDL_BUTTON_RIGHT: RMB = false; break;
            }; break;
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_ESCAPE:   quit(); break;
                //case SDLK_SPACE:    STOP = !STOP; printf( STOP ? " STOPED\n" : " UNSTOPED\n"); break;
                case SDLK_KP_MINUS: zoom*=VIEW_ZOOM_STEP; break;
                case SDLK_KP_PLUS:  zoom/=VIEW_ZOOM_STEP; break;
            } break;
        case SDL_WINDOWEVENT:
            switch (event.window.event) {
                case SDL_WINDOWEVENT_CLOSE:{
                    printf( "window[%i] SDL_WINDOWEVENT_CLOSE \n", id );
                    delete this;
				} break;
            } break;
        //case SDL_QUIT: quit(); break;
    };
};

void ScreenSDL2OGL::keyStateHandling( const Uint8 *keys ){
    //printf( "camStep %g \n", camStep );
    if( keys[ SDL_SCANCODE_LEFT  ] ){ camX0 -= camStep; }
	if( keys[ SDL_SCANCODE_RIGHT ] ){ camX0 += camStep; }
	if( keys[ SDL_SCANCODE_UP    ] ){ camY0 += camStep; }
	if( keys[ SDL_SCANCODE_DOWN  ] ){ camY0 -= camStep; }
};

void ScreenSDL2OGL::mouseHandling( ){
    mouseState = SDL_GetMouseState   ( &mouseX, &mouseY ); //mouseY=HEIGHT-mouseY; // this is done in mouseUp()
    defaultMouseHandling( mouseX, mouseY );
}


void ScreenSDL2OGL::camera(){
    opengl1renderer.matrixMode( GL_PROJECTION );
    opengl1renderer.loadIdentity();
	opengl1renderer.ortho ( -zoom*ASPECT_RATIO, zoom*ASPECT_RATIO, -zoom, zoom, -VIEW_DEPTH, +VIEW_DEPTH );
	opengl1renderer.translatef( -camX0, -camY0, 0.0f );
	opengl1renderer.matrixMode (GL_MODELVIEW);
}

void ScreenSDL2OGL::cameraHUD(){
    opengl1renderer.matrixMode( GL_PROJECTION );
    opengl1renderer.loadIdentity();
	opengl1renderer.ortho ( 0, WIDTH, 0, HEIGHT, -VIEW_DEPTH, +VIEW_DEPTH );
	opengl1renderer.matrixMode (GL_MODELVIEW);
	opengl1renderer.loadIdentity();
}

//void ScreenSDL2OGL::updateMousePos ( int x, int y ){
//    mouse_begin_x = mouseRight( x );
//    mouse_begin_y = mouseUp   ( y );
//}

void ScreenSDL2OGL::defaultMouseHandling( const int& mouseX, const int& mouseY ){
	mouse_begin_x = mouseRight( mouseX ) + camX0;
	mouse_begin_y = mouseUp   ( mouseY ) + camY0;
    fWIDTH  = zoom*ASPECT_RATIO;
	fHEIGHT = zoom;
	camXmin = camX0 - fWIDTH; camYmin = camY0 - fHEIGHT;
	camXmax = camX0 + fWIDTH; camYmax = camY0 + fHEIGHT;
};

// ============== initialization



void ScreenSDL2OGL::setDefaults(){
	VIEW_DEPTH   = VIEW_DEPTH_DEFAULT;
	ASPECT_RATIO = WIDTH/(float)HEIGHT;
	zoom         = VIEW_ZOOM_DEFAULT;
	//printf(" %f %f %f \n", zoom, ASPECT_RATIO, VIEW_DEPTH  );
	mouse_begin_x  = 0;
	mouse_begin_y  = 0;
}

#define DEBUG_ printf( "DEBUG LINE %i %s %s \n", __LINE__, __FUNCTION__, __FILE__ );

void ScreenSDL2OGL::init( int& id_, int WIDTH_, int HEIGHT_, const char* name ){
	WIDTH  = WIDTH_;
	HEIGHT = HEIGHT_;
	setDefaults();
	// modified according to : http://forums.libsdl.org/viewtopic.php?p=40286
	window = SDL_CreateWindow( "Some_Window", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL | SDL_WINDOW_RESIZABLE);
	glctx  = SDL_GL_CreateContext(window);

	if (glctx == NULL) {
        printf("Failed to create GL context: %s\n", SDL_GetError());

        exit(-1);
    }

	
    //SDL_CreateWindowAndRenderer(WIDTH, HEIGHT, SDL_WINDOW_OPENGL, &window, &renderer);
    if(name==0){
        id = SDL_GetWindowID(window); printf( " win id %i \n", id );
        id_=id;
        char str[40];  sprintf(str, " Window id = %d", id );
        SDL_SetWindowTitle( window, str );
    }else{
        SDL_SetWindowTitle( window, name );
    }

	//setupRenderer();
	//setupOpenGLglobals();
	setLightingNormal();
	//printf( " ASPECT_RATIO %f \n", ASPECT_RATIO );
}

ScreenSDL2OGL::ScreenSDL2OGL( int& id, int WIDTH_, int HEIGHT_, const char* name ){
	init( id, WIDTH_, HEIGHT_, name );
}

ScreenSDL2OGL::~ScreenSDL2OGL(){
    //printf(" ScreenSDL2OGL::~ScreenSDL2OGL() \n");
    SDL_DestroyWindow(window);
    if(parent) parent->removeChild(this);
    //printf(" ScreenSDL2OGL::~ScreenSDL2OGL() DONE \n");
}

void ScreenSDL2OGL::removeChild(ScreenSDL2OGL* child){};
