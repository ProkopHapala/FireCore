
int verbosity = 0;

#include "MolGUI.h"
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
	
	// --------- using argparse & LabdaDict;
	app = new MolGUI( junk, DM.w-100, DM.h-100, NULL );
	funcs["-s"]=[&](const char* s){ app->W->smile_name=s; }; // molecule
	funcs["-x"]=[&](const char* s){ app->W->xyz_name=s;   }; // substrate
	funcs["-g"]=[&](const char* s){ app->W->surf_name=s;  }; // substrate
	//funcs["-l"]=[&](const char* s){ app->W->lvs_name=s;   }; // substrate
	process_args( argc, argv, funcs );

	app->init();
	app->loop( 1000000 );
	return 0;
}
















