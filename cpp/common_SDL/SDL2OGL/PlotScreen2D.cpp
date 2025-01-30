
#ifdef __EMSCRIPTEN__
#include <SDL.h>
#else
#include <SDL2/SDL.h>
#endif
#ifdef __EMSCRIPTEN__
#include <SDL_opengl.h>
#else
#include <SDL2/SDL_opengl.h>
#endif

#include "PlotScreen2D.h"

PlotScreen2D::PlotScreen2D( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( id, WIDTH_, HEIGHT_ ) { }
