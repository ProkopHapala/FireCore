
#ifndef  AppSDL2OGL_3D_h
#define  AppSDL2OGL_3D_h

#include "Vec2.h"
#include "quaternion.h"

#include "Camera.h"
#include "cameraOGL.h"

#include "AppSDL2OGL.h"
//#include "Draw3D.h"

//#include "Camera.h"
//#include "cameraOGL.h"

class AppSDL2OGL_3D : public AppSDL2OGL{ public:
	bool mouseSpinning = false;

	float  mouseRotSpeed   = 0.001;
	float  keyRotSpeed     = 0.01;
    float  cameraMoveSpeed = 0.2f;

	Camera cam;

    bool   bDragging = false;
    Vec3f  ray0_start;
    Vec3f  ray0;
    int perFrame =  100;

	float camDist = 50.0;
	Vec2i spinning_start;

	bool perspective  = false;
	bool first_person = false;

// ============ function declarations

	//virtual void quit(       );
	//virtual void loop( int n );
	//virtual void inputHanding();
	virtual void keyStateHandling( const Uint8 *keys       ) override;
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void mouseHandling   (                         ) override;

	virtual void draw     () override;
	virtual void camera   () override;


	inline Vec3f mouseRay0(){ return cam.rotMat().a*mouse_begin_x + cam.rotMat().b*mouse_begin_y + cam.pos(); }
    inline Vec3f updateRay0(){ ray0 = mouseRay0(); return ray0; }
    //ray0 = (Vec3d)(  cam.rotMat().a*mouse_begin_x  +  cam.rotMat().b*mouse_begin_y  +  cam.pos );
    inline void mouseStartSelectionBox(){ ray0_start = ray0;  bDragging = true; }

	void drawCrosshair( float sz );
    void drawMuseSelectionBox();

	AppSDL2OGL_3D( int& id, int WIDTH_, int HEIGHT_, const char* name=0 );

};

#endif
