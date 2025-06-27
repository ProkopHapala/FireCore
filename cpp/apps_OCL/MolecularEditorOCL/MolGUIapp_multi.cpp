#include <globals.h>

#include "MolGUI.h"
#include "MolWorld_sp3_multi.h"
#include "argparse.h"


MolWorld_sp3_multi* W  =0;
MolGUI*            app =0;

LambdaDict funcs; // functions to be called before initMol
LambdaDict funcs2; // functions to be called after initMol

#ifdef WITH_LUA
    //#include "MolGUIapp_Lua.h"
    #include "MolGUIapp_multi_Lua.h"
#endif // WITH_LUA

//int iParalel=-100; 

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);

	W    = new MolWorld_sp3_multi();
    app  = new MolGUI( junk, DM.w-100, DM.h-100, W );

    #include "MolGUIapp_argv.h"
    DEBUG
    funcs["-m"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &W->nSystems ); }}; // number of systems
    DEBUG
    funcs["-tex"]={1,[&](const char** ss){ int i; sscanf( ss[0], "%i", &i ); printf("ARG -tex=%i \n", i); W->ocl.bUseTexture=i; printf("W->ocl.bUseTexture=%i \n", W->ocl.bUseTexture); }}; // use texture for gridFF, boolean
    DEBUG

#ifdef WITH_LUA
    theLua = initMyLua();
    initMyLua_multi( theLua );
    funcs["-lua0"]={1,[&](const char** ss){ if( Lua::dofile(theLua,ss[0]) ){ printf( "ERROR in funcs[-lua0] dofile(%s) => exit()\n", ss[0] ); exit(0); }; }};
#endif // WITH_LUA
    
    DEBUG
	process_args( argc, argv, funcs );
    DEBUG   
    //W->go.constrs.printSizes();
    W->init();
    DEBUG
    //W->go.constrs.printSizes();
    app->bindMolWorld( W );
    //W->go.constrs.printSizes();

    // // ========== After initMol
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
    //funcs["-iParalel"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &W->iParalel); printf( "ARG -iParalel %i \n", W->iParalel ); }};   
    process_args( argc, argv, funcs2, false );

    // --- Apply after-initialization settings and hacks 
    //if(iParalel>-100){ W->iParalel=iParalel; printf( "#### W->iParalel set to %i \n", W->iParalel ); };


    W->pre_loop();
	app->loop( 1000000 );
    W->database->print();
	return 0;
}
















