/// @file AppSDL2OGL.h
/// @brief Base class for SDL2/OpenGL applications with event loop and input handling.
///
/// AppSDL2OGL extends ScreenSDL2OGL with SDL event processing, key state polling,
/// child window management, and a main loop (wait/loop/quit). Provides virtual
/// methods for draw, drawHUD, eventHandling, keyStateHandling that subclasses override.
/// Used as base class by AppSDL2OGL_3D and directly by simple 2D OpenGL apps.

#ifndef  AppSDL2OGL_h
#define  AppSDL2OGL_h

#include <vector>
#include "ScreenSDL2OGL.h"



class AppSDL2OGL : public ScreenSDL2OGL{ public:
	bool LMB=false,RMB=false;
	int  upTime=0,delay=20,timeSlice=5; //,frameCount=0;
	bool loopEnd    = false, STOP = false;
	//float camStep   = VIEW_MOVE_STEP;
    const Uint8 *keys;

    std::vector<ScreenSDL2OGL*> child_windows;

// ============ function declarations

    void wait(float ms);
	virtual void quit(       );
	void         wait(int ms);
	virtual void loop( int n );
	virtual void inputHanding();
	virtual void eventHandling   ( const SDL_Event& event               );
	//virtual void keyStateHandling( const Uint8 *keys                    );
	//virtual void mouseHandling   ( );
	//void defaultMouseHandling    ( const int& mouseX, const int& mouseY );
    virtual void removeChild(ScreenSDL2OGL* child);

	AppSDL2OGL( int& id, int WIDTH_, int HEIGHT_, const char* name=0 );

};

#endif
