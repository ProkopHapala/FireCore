
#include <globals.h>
//int verbosity = 1;
//int idebug    = 0;
double tick2second=1e-9;

#include "MolGUI.h"
#include "argparse.h"


//MMFFsp3 W;
MolGUI* app=0;
LambdaDict funcs;


int   prelat_nstep=0;
int   prelat_nItrMax=0;
Mat3d prelat_dlvec;

#ifdef WITH_LUA
//#include "LuaUtils.h"
//#include "LuaClass.h"
#include "LuaHelpers.h"

lua_State  * theLua=0;

int l_fixAtom(lua_State *L){
    // LuaCall: fixAtom( ia, true )
    int ia = Lua::getInt(L,1);
    int b  = Lua::getInt(L,2);
    printf( "l_fixAtom(ia=%i,b=%i)\n", ia, b ); 
    //if     (b>0){ app->W->atomFixed[ia] = true;  }
    //else if(b<0){ app->W->atomFixed[ia] = false; }
    return 0; // number of return values to Lua environment
}

int initMyLua(){
    theLua         = luaL_newstate();
    lua_State  * L = theLua;
    luaL_openlibs(L);
    lua_register(L, "fixAtom", l_fixAtom  );
    return 1;
}

#endif // WITH_LUA

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);

	// --------- using argparse & LabdaDict;
	app = new MolGUI( junk, DM.w-100, DM.h-100, NULL );
    MolWorld_sp3* W = app->W;

    funcs["-col_damp"]={6,[&](const char** ss){
        //printf( "ss[0](%s)\n", ss[0] ); printf( "ss[1](%s)\n", ss[1] );printf( "ss[2](%s)\n", ss[2] );printf( "ss[3](%s)\n", ss[3] );printf( "ss[4](%s)\n", ss[4] );exit(0);
        int n; double cB, cNB, cAng=0, cm,dR1,dR2;
        sscanf( ss[0], "%i" , &n   ) ; sscanf( ss[1], "%lf", &cB  ); sscanf( ss[2], "%lf", &cNB ); sscanf( ss[3], "%lf", &cm  ); sscanf( ss[3], "%lf", &dR1  ); sscanf( ss[3], "%lf", &dR2  ); 
        printf( "W->ffl.ndampstep %i collisionDamping %g collisionDamping_NB %g damping_medium %g R1,2(%g,%g)\n", n, cB, cNB, cm, dR1, dR2 );
        W->ffl.colDamp.set( n, cm, cB, cAng, cNB, dR1, dR2 );
        // W->ffl.ndampstep                = n;
        // W->ffl.damping_medium           = cm;
        // W->ffl.collisionDamping         = fmax( 0.0, cB  );
        // W->ffl.collisionDamping_NB      = fmax( 0.0, cNB );
        // W->ffl.bCollisionDamping        = cB  > 0.0;
        // W->ffl.bCollisionDampingNonBond = cNB > 0.0;
        // W->ffl.col_damp_dRcut1          = dR1;
        // W->ffl.col_damp_dRcut2          = dR2;
        printf( "W->ffl.colDamp(n=%i,bond=%g,nonB=%g,medium=%g,R12(%g,%g)\n", W->ffl.colDamp.nstep, W->ffl.colDamp.bond, W->ffl.colDamp.nonB, W->ffl.colDamp.medium, W->ffl.colDamp.dRcut1, W->ffl.colDamp.dRcut2 );
    }};// collision damping parameters

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
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );
        sscanf(ss[0],"%i,%i", &app->W->latscan_n.x, &app->W->latscan_n.y );
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); 
        printf( "app->W->latscan_n(%i,%i) latscan_dlvec ", app->W->latscan_n.x, app->W->latscan_n.y ); printMat(*(app->W->latscan_dlvec)); 
    } }; // test

    funcs["-prelat"]={2,[&](const char** ss){ 
        Mat3d* m=&prelat_dlvec;  DEBUG
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );  DEBUG
        sscanf(ss[0],"%i,%i", &prelat_nstep, &prelat_nItrMax );  DEBUG
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z );  DEBUG
        printf( "prelat_dlvec(%i,%i) latscan_dlvec ", prelat_nstep, prelat_nItrMax  ); printMat(prelat_dlvec);  DEBUG
    } }; // test

    funcs["-uff"]={0,[&](const char** ss){ W->bUFF=true; }}; // AutoCharge

    funcs["-e"]={0,[&](const char** ss){ app->W->bEpairs=true; }}; // add explicit electron pair
    funcs["-EachAngle"]={0,[&](const char** ss){ app->W->ffl.bEachAngle=true;                          }};
    funcs["-torsions"]={0,[&](const char** ss){ app->W->ffl.bTorsion=true; app->W->ffl.doPiPiI=false;  }};
    
    funcs["-substr_iso"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &app->subs_iso ); }};
    funcs["-iMO"       ]={1,[&](const char** ss){ sscanf( ss[0], "%i",  &app->which_MO ); }};
    
    funcs["-perframe"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &W->iterPerFrame ); app->perFrame=W->iterPerFrame; printf( "#### -perframe %i \n", W->iterPerFrame ); }};  // interations per frame

	process_args( argc, argv, funcs );
	app->init();

#ifdef WITH_LUA
    initMyLua();
    app->console.callback = [&](const char* str)->bool{
       lua_State* L=theLua;
       // printf( "console.callback: %s\n", cmd );
        if (luaL_dostring(L, str) != LUA_OK) {
            // If there's an error, print it
            fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
            lua_pop(L, 1);  // Remove error message from the stack
            return false;
        }
        return true;
    };
#endif // WITH_LUA

    if(prelat_nstep>0)app->W->change_lvec_relax( prelat_nstep, prelat_nItrMax, 1e-3, prelat_dlvec );
	app->loop( 1000000 );
	return 0;
}
















