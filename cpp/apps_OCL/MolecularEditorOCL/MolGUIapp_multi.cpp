#include <globals.h>
//int verbosity = 1;
int idebug    = 0;
double tick2second=1e-9;

#include "MolGUI.h"
#include "MolWorld_sp3_multi.h"
#include "argparse.h"

//MMFFsp3 W;
MolGUI* app=0;
LambdaDict funcs;

int iParalel=-100; 


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
	funcs["-m"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &W->nSystems ); }}; // number of systems
	funcs["-s"]={1,[&](const char** ss){ W->smile_name=ss[0]; }}; // molecule as SMILEs
	funcs["-x"]={1,[&](const char** ss){ W->xyz_name  =ss[0]; }}; // molecule as .xyz
	funcs["-g"]={1,[&](const char** ss){ W->surf_name =ss[0]; }}; // substrate as .xyz
	funcs["-r"]={0,[&](const char** ss){ W->bMMFF=false;      }}; // rigid
	funcs["-n"]={1,[&](const char** ss){ W->nMulPBC.x=(ss[0][0]-'0'); W->nMulPBC.y=(ss[0][1]-'0'); W->nMulPBC.z=(ss[0][2]-'0');      }}; // multiply lattice cell
	funcs["-q"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &(W->fAutoCharges) ); }}; // AutoCharge
	funcs["-t"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &(W->itest) ); }}; // test
    funcs["-c"]={1,[&](const char** ss){ int iconstr; sscanf( ss[0], "%i", &iconstr ); W->constrain_list.push_back(iconstr); }}; // constrain atom
    funcs["-e"]={0,[&](const char** ss){ W->bEpairs=true; }}; // add explicit electron pairs
	funcs["-ManipulAnim"]={0,[&](const char** ss){ W->bAnimManipulation=true; }}; // AnimateManipulation
    //funcs["-EachAngle"]={0,[&](const char** ss){ W->ffl.bEachAngle=true;                          }};
    //funcs["-torsions"]={0,[&](const char** ss){ W->ffl.bTorsion=true; W->ffl.doPiPiI=false;  }};

    funcs["-iParalel"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &iParalel        );                                printf( "#### -iParalel %i \n", iParalel ); }};         // paralelization model
    funcs["-perframe"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &W->iterPerFrame ); app->perFrame=W->iterPerFrame; printf( "#### -perframe %i \n", W->iterPerFrame ); }};  // interations per frame

    funcs["-b"  ]={1,[&](const char** ss){ W->bConstrains=true; W->constr_name=ss[0]; }};  // constrain bond lenghts
    funcs["-pop"]={1,[&](const char** ss){ W->uploadPopName=ss[0]; }};                     // upload population
    funcs["-dlvec"]={1,[&](const char** ss){ Mat3d* m=new Mat3d(); W->dlvec=m; printf( "ss[0] `%s`\n", ss[0] ); sscanf(ss[0],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); printf( "W->dlvec set to " ); printMat(*(W->dlvec)); } }; // modify lattice vector
    
    funcs["-substr_iso"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &app->subs_iso ); }};
    funcs["-iMO"       ]={1,[&](const char** ss){ sscanf( ss[0], "%i",  &app->which_MO ); }};

    funcs["-latscan"]={2,[&](const char** ss){  // scan lattice vector
        W->bLatScan=true;
        Mat3d* m=new Mat3d(); W->latscan_dlvec=m; 
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );
        sscanf(ss[0],"%i,%i", &W->latscan_n.x, &W->latscan_n.y );
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z );
        printf( "W->latscan_n(%i,%i) latscan_dlvec ", W->latscan_n.x, W->latscan_n.y ); printMat(*(W->latscan_dlvec));
    } }; // test

	process_args( argc, argv, funcs );
	app->init();

    // --- Apply after-initialization settings and hacks 
    if(iParalel>-100){ W->iParalel=iParalel; printf( "#### W->iParalel set to %i \n", W->iParalel ); };

	app->loop( 1000000 );
	return 0;
}
















