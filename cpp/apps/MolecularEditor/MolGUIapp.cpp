
#include <globals.h>

#include "MolGUI.h"
#include "argparse.h"

MolWorld_sp3* W   =0;
MolGUI*       app =0;

// There are 3 stages of initialization command line arguments can be processed (  1. After MolGUI initialization and before initMol, 2. after initMol, ) 
LambdaDict funcs; // functions to be called before initMol
LambdaDict funcs2; // functions to be called after initMol

int   prelat_nstep=0;
int   prelat_nItrMax=0;
Mat3d prelat_dlvec;

#ifdef WITH_LUA
//#include "LuaUtils.h"
//#include "LuaClass.h"
//#include "LuaHelpers.h"
//lua_State  * theLua=0;
#include "MolGUIapp_Lua.h"
#endif // WITH_LUA

#ifdef __EMSCRIPTEN__
void main_loop(){
    if (app->GL_LOCK) return;
    app->loop(1);
}
#endif

int main(int argc, char *argv[]){
#ifdef __EMSCRIPTEN__
    argc = 5;
    char* new_argv[5];
    argv = new_argv;
    argv[0] = "./MolGUIapp";
    argv[1] = "-x";
    argv[2] = "common_resources/xyz/PTCDA";
    argv[3] = "-g";
    argv[4] = "common_resources/xyz/NaCl_1x1_L3";
#endif

	SDL_Init(SDL_INIT_VIDEO);

    // === GLES 3 setup ===
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_PROFILE_MASK, SDL_GL_CONTEXT_PROFILE_ES);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 0);
    //SDL_GL_SetAttribute(SDL_GL_RED_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE, 8);
    //SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE, 16);
    //SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);

	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 8);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);

#ifdef DEBUG_ALLOCATOR
    debugAllocator_init();
#endif

    // --------- Initialize MolGUI and MolWorld_sp3

    W   = new MolWorld_sp3();
    app = new MolGUI( junk, DM.w-100, DM.h-100, W );

    //printf( "WE ARE HERE %s \n", _CODE_LOCATION ); 
    //printf( "WE ARE HERE %s \n", __FUNCTION__ ); 
    //exit(0);

    #include "MolGUIapp_argv.h"

    funcs["-prelat"]={2,[&](const char** ss){ 
        Mat3d* m=&prelat_dlvec;
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );
        sscanf(ss[0],"%i,%i", &prelat_nstep, &prelat_nItrMax );
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z );
        printf( "prelat_dlvec(%i,%i) latscan_dlvec \n", prelat_nstep, prelat_nItrMax  ); printMat(prelat_dlvec);
    } }; // test

    // ========== Before initMol
#ifdef WITH_LUA
    theLua = initMyLua();
    funcs["-lua0"]={1,[&](const char** ss){ if( Lua::dofile(theLua,ss[0]) ){ printf( "ERROR in funcs[-lua0] dofile(%s) => exit()\n", ss[0] ); exit(0); }; }};
#endif // WITH_LUA

	// --------- using argparse & LabdaDict;

    process_args( argc, argv, funcs, false );
    W->init();
    app->bindMolWorld( W );

    // ========== After initMol

#ifdef WITH_LUA
    app->console.callback = [&](const char* str)->bool{
       lua_State* L=theLua;
        printf( "console.lua_callback: %s\n", str );
        if (luaL_dostring(L, str) != LUA_OK) {
            // If there's an error, print it
            //fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
            printf( "Error: %s\n", lua_tostring(L, -1) );
            lua_pop(L, 1);  // Remove error message from the stack
            return false;
        }
        return true;
    };
    funcs2["-lua"]={1,[&](const char** ss){ if( Lua::dofile(theLua,ss[0]) ){ printf( "ERROR in funcs[-lua] dofile(%s) => exit()\n", ss[0] ); exit(0); }; }};
#endif // WITH_LUA
    process_args( argc, argv, funcs2, false );
    
    if(prelat_nstep>0)W->change_lvec_relax( prelat_nstep, prelat_nItrMax, 1e-3, prelat_dlvec );
    W->pre_loop();
    #ifdef __EMSCRIPTEN__
        emscripten_set_main_loop(main_loop, 0, 1);
    #else
        app->loop( 1000000 );
    #endif

    W->clear( true, true );
#ifdef DEBUG_ALLOCATOR
    debugAllocator_print( );
#endif    
	return 0;
}
