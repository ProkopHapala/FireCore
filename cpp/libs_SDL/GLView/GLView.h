/// @file GLView.h
/// @brief C API header for GLView — lightweight standalone OpenGL 3D viewer window.
///
/// @see GLView.hpp for the C++ class definition and full description.

#ifndef  GLView_h
#define  GLView_h

#include <stdbool.h>

extern "C" {
typedef void (*ProcedurePointer)();

void init( int w, int h );
bool draw();
bool pre_draw();
bool post_draw();
void run_Nframes(int nframes);
void set_draw_function( ProcedurePointer draw_func );
void getMousePos(double* x, double* y);

} // extern "C" {

#endif
