
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
	funcs["-n"]={1,[&](const char** ss){  app->W->nMulPBC.x=(ss[0][0]-'0'); app->W->nMulPBC.y=(ss[0][1]-'0'); app->W->nMulPBC.z=(ss[0][2]-'0'); }}; // PBC multiplication of molecule
	funcs["-ng"]={1,[&](const char** ss){ app->W->bCellBySurf=true; sscanf(ss[0],"%lf,%lf,%lf,%lf", &app->W->bySurf_lat[0].x,&app->W->bySurf_lat[0].y,  &app->W->bySurf_lat[1].x,&app->W->bySurf_lat[1].y ); }}; // change molecule cell by surface multiple
	funcs["-subs"]={1,[&](const char** ss){ app->W->substitute_name=new char[256]; sscanf(ss[0],"%i,%s", &app->W->isubs, app->W->substitute_name ); }}; // substitute group on molecule
	funcs["-q"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &(app->W->fAutoCharges) ); }}; // AutoCharge
	funcs["-t"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &(app->W->itest) ); }}; // test
    funcs["-c"]={1,[&](const char** ss){ int iconstr; sscanf( ss[0], "%i", &iconstr ); app->W->constrain_list.push_back(iconstr); }}; // 
    //funcs["-b"]={1,[&](const char** ss){ app->W->bConstrains=true; app->W->constrs.loadBonds( ss[0], &app->W->builder.atom_permut[0] ); }}; // constrains must be loaded after initialization of geometry
    funcs["-b"]={1,[&](const char** ss){ app->W->bConstrains=true; app->W->constr_name=ss[0]; }}; // test
    funcs["-dlvec"]={1,[&](const char** ss){ Mat3d* m=new Mat3d(); app->W->dlvec=m; printf( "ss[0] `%s`\n", ss[0] ); sscanf(ss[0],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); printf( "app->W->dlvec set to " ); printMat(*(app->W->dlvec)); } }; // test
    funcs["-latscan"]={2,[&](const char** ss){ 
        app->W->bLatScan=true;
        Mat3d* m=new Mat3d(); app->W->latscan_dlvec=m; 
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] ); DEBUG
        sscanf(ss[0],"%i,%i", &app->W->latscan_n.x, &app->W->latscan_n.y ); DEBUG
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); DEBUG
        printf( "app->W->latscan_n(%i,%i) latscan_dlvec ", app->W->latscan_n.x, app->W->latscan_n.y ); printMat(*(app->W->latscan_dlvec)); DEBUG
    } }; // test

    funcs["-e"]={0,[&](const char** ss){ app->W->bEpairs=true; }}; // add explicit electron pair
    funcs["-EachAngle"]={0,[&](const char** ss){ app->W->ffl.bEachAngle=true;                          }};
    funcs["-torsions"]={0,[&](const char** ss){ app->W->ffl.bTorsion=true; app->W->ffl.doPiPiI=false;  }};
    



	process_args( argc, argv, funcs );
	app->init();
	app->loop( 1000000 );
	return 0;
}
















