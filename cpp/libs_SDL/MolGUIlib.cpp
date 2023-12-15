
#ifndef  GLView_cpp
#define  GLView_cpp

constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

#include <globals.h>
//int verbosity = 1;
int idebug    = 0;
double tick2second=1e-9;

#include "MolGUI.h"

MolGUI* gui=0;

extern "C"{

void init( void* W ){
    SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);

    gui = new MolGUI( junk, DM.w-100, DM.h-100, (MolWorld_sp3*)W );
	gui->initGUI();
	//gui->loop( nframes );
}

void run(int nframes){
    printf( "!!!!!!!!!!!!!!!! gui->W: %li \n", (long)(gui->W) );
    gui->loop( nframes );
}

} // extern "C"{

#endif
