
#ifndef  cameraOGL_h
#define  cameraOGL_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "raytrace.h"
#include "Camera.h"
#include "Draw3D.h"

namespace Cam{

inline void ortho( const Camera& cam, bool zsym ){
    opengl1renderer.matrixMode( GL_PROJECTION );
    opengl1renderer.loadIdentity();
    //float zmin = cam.zmin; if(zsym) zmin=-cam.zmax;
	opengl1renderer.ortho( -cam.zoom()*cam.aspect(), cam.zoom()*cam.aspect(), -cam.zoom(), cam.zoom(), cam.zmin, cam.zmax );
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, cam.rotMat(), glMat );
	opengl1renderer.multMatrixf( glMat );

	opengl1renderer.matrixMode ( GL_MODELVIEW );
	opengl1renderer.loadIdentity();
	opengl1renderer.translatef(-cam.pos().x,-cam.pos().y,-cam.pos().z);
}

inline void perspective( const Camera& cam ){
    opengl1renderer.matrixMode( GL_PROJECTION );
    opengl1renderer.loadIdentity();
    //opengl1renderer.frustum( -ASPECT_RATIO, ASPECT_RATIO, -1, 1, camDist/zoom, VIEW_DEPTH );
    opengl1renderer.frustum( -cam.aspect()*cam.zoom(), cam.aspect()*cam.zoom(), -cam.zoom(), cam.zoom(), cam.zmin, cam.zmax );
    //opengl1renderer.frustum( -cam.zoom*cam.aspect, cam.zoom*cam.aspect, -cam.zoom, cam.zoom, cam.zmin, cam.zmax );
	float glMat[16];
	Draw3D::toGLMatCam( { 0.0f, 0.0f, 0.0f}, cam.rotMat(), glMat );
	opengl1renderer.multMatrixf( glMat );

	opengl1renderer.matrixMode ( GL_MODELVIEW );
	opengl1renderer.loadIdentity();
	opengl1renderer.translatef(-cam.pos().x,-cam.pos().y,-cam.pos().z);
    //opengl1renderer.translatef ( -camPos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
    //opengl1renderer.translatef ( -cam.pos.x+camMat.cx*camDist, -camPos.y+camMat.cy*camDist, -camPos.z+camMat.cz*camDist );
}

} // namespace Cam

#endif

