
int verbosity = 1;
int idebug    = 0;
double tick2second=1e-9;

#include "MolGUI.h"
#include "MolWorld_sp3_multi.h"
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

	MolWorld_sp3_multi* W = new MolWorld_sp3_multi();
	// --------- using argparse & LabdaDict;
	app = new MolGUI( junk, DM.w-100, DM.h-100, W );
	funcs["-m"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &W->nSystems ); printf(  "!!!!!!!!!!!!!!! DEBUG MolGUIapp_multi: W->nSystems %i \n", W->nSystems ); }}; // number of systems
	funcs["-s"]={1,[&](const char** ss){ W->smile_name=ss[0]; }}; // molecule as SMILEs
	funcs["-x"]={1,[&](const char** ss){ W->xyz_name  =ss[0]; }}; // molecule as .xyz
	funcs["-g"]={1,[&](const char** ss){ W->surf_name =ss[0]; }}; // substrate as .xyz
	funcs["-r"]={0,[&](const char** ss){ W->bMMFF=false;      }}; // rigid
	funcs["-n"]={0,[&](const char** ss){ W->nMulPBC.x=(ss[0][0]-'0'); W->nMulPBC.y=(ss[0][1]-'0'); W->nMulPBC.z=(ss[0][2]-'0');      }}; // rigid
	funcs["-q"]={0,[&](const char** ss){ sscanf( ss[0], "%lf", &(W->fAutoCharges) ); }}; // AutoCharge
	funcs["-t"]={0,[&](const char** ss){ sscanf( ss[0], "%i", &(app->W->itest) ); }}; // test
    funcs["-c"]={0,[&](const char** ss){ int iconstr; sscanf( ss[0], "%i", &iconstr ); app->W->constrain_list.push_back(iconstr); }}; // test
	process_args( argc, argv, funcs );
	app->init();
	app->loop( 1000000 );
	return 0;
}
















