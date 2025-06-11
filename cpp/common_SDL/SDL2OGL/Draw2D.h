
#ifndef  Draw2D_h
#define  Draw2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "geom2D.h"
#include "Renderer.h"

#define COL2VEC(c) ((Vec3f){(c>>16 & 0xff)/255.0f, (c>>8 & 0xff)/255.0f, (c & 0xff)/255.0f} )

namespace Draw2D{

extern float z_layer; // this can be used for drawing in multiple overlaping layers
                      // should be initialized like this http://stackoverflow.com/questions/19136059/namespace-global-variable-losing-value-c


void drawPointCross ( const Vec2f& vec, float d);
void drawLine       ( const Vec2f& p1, const Vec2f& p2   );
void drawRectangle  ( float p1x, float p1y, float p2x, float p2y, Vec3f color, bool filled = true);
void drawRectangle  ( const Vec2f& p1, const Vec2f& p2, Vec3f color, bool filled=true);
void drawCircle     ( const Vec2f& center, float radius, bool filled );

void plot       ( int n, float dx,    double * ys );
void plot       ( int n, double * xs, double * ys );
void plot_dots  ( int n, double * xs, double * ys );
void plot_cross ( int n, double * xs, double * ys, double sz );
void plot_X     ( int n, double * xs, double * ys, double sz );
void plot_O     ( int n, double * xs, double * ys, double sz );

void drawGrid   ( int n, double * ticks, double lmin, double lmax, bool XorY );

void drawText ( const char * str, int nchar, Vec2d pos, float angle, int fontTex, float textSize );
void drawText ( const char * str,            Vec2d pos, Vec2d sz,    int fontTex, float textSize );


// ==== inline functions
inline void toGLMat( const Vec2d& pos, const Vec2d& rot, float* glMat ){
	glMat[0 ] = +rot.x; glMat[1 ] = +rot.y; glMat[2 ] = 0;  glMat[3 ] = 0;
	glMat[4 ] = -rot.y; glMat[5 ] = +rot.x; glMat[6 ] = 0;  glMat[7 ] = 0;
	glMat[8 ] = 0;      glMat[9 ] = 0;      glMat[10] = 1;  glMat[11] = 0;
	glMat[12] = pos.x;  glMat[13] = pos.y;  glMat[14] = 0;  glMat[15] = 1;
};

inline void toGLMatCam( const Vec2d& pos, const Vec2d& rot, float* glMat ){
	glMat[0 ] = +rot.x; glMat[1 ] = -rot.y; glMat[2 ] = 0;  glMat[3 ] = 0;
	glMat[4 ] = +rot.y; glMat[5 ] = +rot.x; glMat[6 ] = 0;  glMat[7 ] = 0;
	glMat[8 ] = 0;      glMat[9 ] = 0;      glMat[10] = 1;  glMat[11] = 0;
	glMat[12] = -pos.x; glMat[13] = -pos.y; glMat[14] = 0;  glMat[15] = 1;
};

}; // namespace Draw2D

#endif

