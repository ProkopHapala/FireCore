
#include "ScreenSDL2OGL.h" // THE HEADER
#include "Renderer.h"
#include <SDL2/SDL_video.h>

//#include "testUtils.h"

// ============== per frame

void setLightingNormal(){
//float white    [] = { 1.0f, 1.0f,  1.0f,  1.0f };
	float ambient  []{ 1.0f, 1.0f, 1.0f, 1.0f };
	float diffuse  []{ 2.0f, 1.0f, 1.0f, 1.0f };
	float specular []{ 1.0f, 1.0f, 1.0f, 1.0f };
	float shininess[]{ 128, 1.0f }; // exponent for specular
	float emission []{ 0.0f, 0.0f, 0.0f, 1.0f }; // as light source

	opengl1renderer.enable     ( GL_COLOR_MATERIAL    );
    opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient   );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse   );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular  );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission  );

    float lightPos   []{ 1.0f, -1.0f, 1.0f, 0.0f  };
    opengl1renderer.lightfv    ( GL_LIGHT0, GL_POSITION,  lightPos );
    //float l_ambient  []{ 0.1f*0, 0.15f*0, 0.2f*0,  1.0f };
    float l_ambient  []{ 0.15f, 0.20f, 0.25f,  1.0f };
    //float l_ambient  []{ 0.2f, 0.2f, 0.2f,  1.0f };
	float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 1.0f, 1.0f,  1.0f,  1.0f };
    opengl1renderer.lightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

	opengl1renderer.lightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER, 1 );

	//float ambient_lm  []{ 0.15f, 0.2f, 0.25f, 1.0f };
	float ambient_lm  []{ 0.2f, 0.2f, 0.2f, 1.0f };
	opengl1renderer.lightModelfv( GL_LIGHT_MODEL_AMBIENT, ambient_lm );

	opengl1renderer.enable     ( GL_LIGHTING         );
	opengl1renderer.enable     ( GL_LIGHT0           );
	opengl1renderer.enable     ( GL_NORMALIZE        );

	opengl1renderer.enable     ( GL_DEPTH_TEST       );
	opengl1renderer.hint       ( GL_LINE_SMOOTH_HINT, GL_NICEST );
	opengl1renderer.shadeModel ( GL_SMOOTH           );
	opengl1renderer.polygonMode( GL_FRONT_AND_BACK, GL_FILL );
    //opengl1renderer.lightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );
}

void setLightingRGB(){
//float white    [] = { 1.0f, 1.0f,  1.0f,  1.0f };
	float ambient  []{ 1.0f, 1.0f, 1.0f, 1.0f };
	float diffuse  []{ 2.0f, 1.0f, 1.0f, 1.0f };
	float specular []{ 1.0f, 1.0f, 1.0f, 1.0f };
	float shininess[]{ 128, 1.0f }; // exponent for specular
	float emission []{ 0.0f, 0.0f, 0.0f, 1.0f }; // as light source

	opengl1renderer.enable     ( GL_COLOR_MATERIAL    );
    opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient   );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse   );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular  );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_SHININESS, shininess );
	opengl1renderer.materialfv( GL_FRONT_AND_BACK, GL_EMISSION,  emission  );

    float lightPos_r   []{ -1.0f,  0.0f,  0.0f, 0.0f  };
    float lightPos_g   []{  0.0f, -1.0f,  0.0f, 0.0f  };
    float lightPos_b   []{  0.0f,  0.0f, -1.0f, 0.0f  };

	opengl1renderer.enable     ( GL_LIGHTING         );

	opengl1renderer.enable     ( GL_NORMALIZE        );
	opengl1renderer.enable     ( GL_LIGHT0           );
    opengl1renderer.enable     ( GL_LIGHT1           );
    opengl1renderer.enable     ( GL_LIGHT2           );

    opengl1renderer.lightfv    ( GL_LIGHT0, GL_POSITION,  lightPos_r );
    opengl1renderer.lightfv    ( GL_LIGHT1, GL_POSITION,  lightPos_g );
    opengl1renderer.lightfv    ( GL_LIGHT2, GL_POSITION,  lightPos_b );

    float l_ambient  []{ 0.1f, 0.1f, 0.1f,  1.0f };
	float l_specular []{ 1.0f, 1.0f,  1.0f,  1.0f };

	float r_diffuse  []{ 0.8f, 0.2f, 0.2f,  1.0f };
	float g_diffuse  []{ 0.2f, 0.8f, 0.2f,  1.0f };
	float b_diffuse  []{ 0.2f, 0.2f, 0.8f,  1.0f };
    opengl1renderer.lightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_DIFFUSE,   r_diffuse  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    opengl1renderer.lightfv    ( GL_LIGHT1, GL_AMBIENT,   l_ambient  );
	opengl1renderer.lightfv    ( GL_LIGHT1, GL_DIFFUSE,   g_diffuse );
	opengl1renderer.lightfv    ( GL_LIGHT1, GL_SPECULAR,  l_specular );

    opengl1renderer.lightfv    ( GL_LIGHT2, GL_AMBIENT,   l_ambient  );
	opengl1renderer.lightfv    ( GL_LIGHT2, GL_DIFFUSE,   b_diffuse  );
	opengl1renderer.lightfv    ( GL_LIGHT2, GL_SPECULAR,  l_specular );

	opengl1renderer.lightModeli( GL_LIGHT_MODEL_LOCAL_VIEWER, 1 );

	//float ambient_lm  []{ 0.15f, 0.2f, 0.25f, 1.0f };
	float ambient_lm  []{ 0.2f, 0.2f, 0.2f, 1.0f };
	opengl1renderer.lightModelfv( GL_LIGHT_MODEL_AMBIENT, ambient_lm );

	opengl1renderer.enable     ( GL_DEPTH_TEST       );
	opengl1renderer.hint       ( GL_LINE_SMOOTH_HINT, GL_NICEST );
	opengl1renderer.shadeModel ( GL_SMOOTH           );
	opengl1renderer.polygonMode( GL_FRONT_AND_BACK, GL_FILL );
    //opengl1renderer.lightModeli( GL_LIGHT_MODEL_TWO_SIDE, 1 );
}

void ScreenSDL2OGL::update( ){
	if( GL_LOCK ){ printf("ScreenSDL2OGL::update GL_LOCK\n"); return; }
	GL_LOCK = true;
	//printf( " window[%i] SDL_GL_MakeCurrent \n", id );
    SDL_GL_MakeCurrent(window, glctx);
	camera();

	glClearColor(0.2f, 0.5f, 0.8f, 1.0f); // TODO: replace with renderer->clearColor
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	draw();
	cameraHUD();
	drawHUD();
	
	//opengl1renderer.popMatrix();
	//opengl1renderer.flush();
	//SDL_RenderPresent(renderer);
	frameCount++;
    SDL_GL_SwapWindow(window);
    //printf( " window[%i] SDL_GL_SwapWindow \n", id );
    GL_LOCK = false;
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
	window = SDL_CreateWindow( "Some_Window", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, WIDTH, HEIGHT, SDL_WINDOW_OPENGL);
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




	renderer = new Renderer();
	opengl1renderer.bind_renderer(renderer);

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
