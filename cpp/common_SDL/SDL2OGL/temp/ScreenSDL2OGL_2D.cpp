
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

#include "ScreenSDL2OGL_2D.h" // THE HEADER

ScreenSDL2OGL_2D::ScreenSDL2OGL_2D( int& id, int WIDTH_, int HEIGHT_ ) : ScreenSDL2OGL( int& id, int WIDTH_, int HEIGHT_ ) { }

