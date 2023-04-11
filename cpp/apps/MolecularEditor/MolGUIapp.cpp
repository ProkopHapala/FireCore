
int verbosity = 1;
int idebug    = 0;
double tick2second=1e-9;

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

	// --------- using argparse & LabdaDict;
	app = new MolGUI( junk, DM.w-100, DM.h-100, NULL );
	funcs["-s"]={1,[&](const char** ss){ app->W->smile_name=ss[0]; }}; // molecule as SMILEs
	funcs["-x"]={1,[&](const char** ss){ app->W->xyz_name  =ss[0]; }}; // molecule as .xyz
	funcs["-g"]={1,[&](const char** ss){ app->W->surf_name =ss[0]; }}; // substrate as .xyz
	funcs["-r"]={0,[&](const char** ss){ app->W->bMMFF=false;      }}; // rigid
	funcs["-n"]={0,[&](const char** ss){  app->W->nMulPBC.x=(ss[0][0]-'0'); app->W->nMulPBC.y=(ss[0][1]-'0'); app->W->nMulPBC.z=(ss[0][2]-'0'); }}; // PBC multiplication of molecule
	funcs["-ng"]={0,[&](const char** ss){ app->W->bCellBySurf=true; sscanf(ss[0],"%lf,%lf,%lf,%lf", &app->W->bySurf_lat[0].x,&app->W->bySurf_lat[0].y,  &app->W->bySurf_lat[1].x,&app->W->bySurf_lat[1].y ); }}; // change molecule cell by surface multiple
	funcs["-subs"]={1,[&](const char** ss){ app->W->substitute_name=new char[256]; sscanf(ss[0],"%i,%s", &app->W->isubs, app->W->substitute_name ); }}; // substitute group on molecule
	funcs["-q"]={0,[&](const char** ss){ sscanf( ss[0], "%lf", &(app->W->fAutoCharges) ); }}; // AutoCharge
	funcs["-t"]={0,[&](const char** ss){ sscanf( ss[0], "%i", &(app->W->itest) ); }}; // test

	process_args( argc, argv, funcs );
	app->init();
	app->loop( 1000000 );
	return 0;
}
















