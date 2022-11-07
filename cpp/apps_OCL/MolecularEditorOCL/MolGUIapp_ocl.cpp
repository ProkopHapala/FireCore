
int verbosity = 1;
int idebug    = 0;

double tick2second=1e-9;

#include "MolGUI.h"
#include "MolWorld_sp3_ocl.h"

#include "argparse.h"


//MMFFsp3 W;
MolGUI* app=0;
LambdaDict funcs;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);
    //char* ssmile = 0; if(argc>1){ ssmile=argv[1]; } // initialize with smiles ?
	//app = new MolGUI( junk, DM.w-100, DM.h-100, NULL, ssmile );
	//if(argc>2){ app->W->loadSurf( argv[2] ); }      // load surface
	MolWorld_sp3_ocl* W = new MolWorld_sp3_ocl();
	// --------- using argparse & LabdaDict;
	app = new MolGUI( junk, DM.w-100, DM.h-100, W );
	funcs["-s"]={1,[&](const char** ss){ app->W->smile_name=ss[0]; }}; // molecule as SMILEs
	funcs["-x"]={1,[&](const char** ss){ app->W->xyz_name  =ss[0]; }}; // molecule as .xyz
	funcs["-g"]={1,[&](const char** ss){ app->W->surf_name =ss[0]; }}; // substrate as .xyz
	funcs["-r"]={0,[&](const char** ss){ app->W->bMMFF=false;      }}; // rigid
	process_args( argc, argv, funcs );
	app->init();
	app->loop( 1000000 );
	return 0;
}
















