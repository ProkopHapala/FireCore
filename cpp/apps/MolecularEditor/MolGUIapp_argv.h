#ifndef MolGUIapp_argv_h
#define MolGUIapp_argv_h

    funcs["-col_damp"]={6,[&](const char** ss){
        //printf( "ss[0](%s)\n", ss[0] ); printf( "ss[1](%s)\n", ss[1] );printf( "ss[2](%s)\n", ss[2] );printf( "ss[3](%s)\n", ss[3] );printf( "ss[4](%s)\n", ss[4] );exit(0);
        int n; double cB, cNB, cAng=0, cm,dR1,dR2;
        sscanf( ss[0], "%i" , &n   ) ; sscanf( ss[1], "%lf", &cB  ); sscanf( ss[2], "%lf", &cNB ); sscanf( ss[3], "%lf", &cm  ); sscanf( ss[3], "%lf", &dR1  ); sscanf( ss[3], "%lf", &dR2  ); 
        printf( "W->ffl.ndampstep %i collisionDamping %g collisionDamping_NB %g damping_medium %g R1,2(%g,%g)\n", n, cB, cNB, cm, dR1, dR2 );
        W->ffl.colDamp.set( n, cm, cB, cAng, cNB, dR1, dR2 );
        printf( "W->ffl.colDamp(n=%i,bond=%g,nonB=%g,medium=%g,R12(%g,%g)\n", W->ffl.colDamp.nstep, W->ffl.colDamp.bond, W->ffl.colDamp.nonB, W->ffl.colDamp.medium, W->ffl.colDamp.dRcut1, W->ffl.colDamp.dRcut2 );
    }};// collision damping parameters

	funcs["-s"]={1,[&](const char** ss){ W->smile_name=ss[0]; }}; // molecule as SMILEs
	funcs["-x"]={1,[&](const char** ss){ W->xyz_name  =ss[0]; }}; // molecule as .xyz
	funcs["-g"]={1,[&](const char** ss){ W->surf_name =ss[0]; }}; // substrate as .xyz
	funcs["-r"]={0,[&](const char** ss){ W->bMMFF=false;      }}; // rigid
	funcs["-n"]={1,[&](const char** ss){  W->nMulPBC.x=(ss[0][0]-'0'); W->nMulPBC.y=(ss[0][1]-'0'); W->nMulPBC.z=(ss[0][2]-'0'); }}; // PBC multiplication of molecule
	funcs["-ng"]={1,[&](const char** ss){ W->bCellBySurf=true; sscanf(ss[0],"%lf,%lf,%lf,%lf", &W->bySurf_lat[0].x,&W->bySurf_lat[0].y,  &W->bySurf_lat[1].x,&W->bySurf_lat[1].y ); }}; // change molecule cell by surface multiple
	funcs["-subs"]={1,[&](const char** ss){ W->substitute_name=new char[256]; sscanf(ss[0],"%i,%s", &W->isubs, W->substitute_name ); }}; // substitute group on molecule
	funcs["-q"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &(W->fAutoCharges) ); }}; // AutoCharge
	funcs["-t"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &(W->itest) ); }}; // test
    funcs["-c"]={1,[&](const char** ss){ int iconstr; sscanf( ss[0], "%i", &iconstr ); W->constrain_list.push_back(iconstr); }}; // 
    //funcs["-b"]={1,[&](const char** ss){ W->bConstrains=true; W->constrs.loadBonds( ss[0], &W->builder.atom_permut[0] ); }}; // constrains must be loaded after initialization of geometry
    funcs["-b"]={1,[&](const char** ss){ W->bConstrains=true; W->constr_name=ss[0]; }}; // test
    funcs["-dlvec"]={1,[&](const char** ss){ Mat3d* m=new Mat3d(); W->dlvec=m; printf( "ss[0] `%s`\n", ss[0] ); sscanf(ss[0],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); printf( "W->dlvec set to " ); printMat(*(W->dlvec)); } }; // test
    funcs["-latscan"]={2,[&](const char** ss){ 
        W->bLatScan=true;
        Mat3d* m=new Mat3d(); W->latscan_dlvec=m; 
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );
        sscanf(ss[0],"%i,%i", &W->latscan_n.x, &W->latscan_n.y );
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); 
        printf( "W->latscan_n(%i,%i) latscan_dlvec ", W->latscan_n.x, W->latscan_n.y ); printMat(*(W->latscan_dlvec)); 
    } }; // test

    funcs["-prelat"]={2,[&](const char** ss){ 
        Mat3d* m=&prelat_dlvec;
        printf( "ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );
        sscanf(ss[0],"%i,%i", &prelat_nstep, &prelat_nItrMax );
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z );
        printf( "prelat_dlvec(%i,%i) latscan_dlvec ", prelat_nstep, prelat_nItrMax  ); printMat(prelat_dlvec);
    } }; // test

    // set verbosity
    funcs["-verb"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &verbosity ); }};

    funcs["-uff"]={0,[&](const char** ss){ W->bUFF=true; }}; // AutoCharge

    funcs["-e"]={0,[&](const char** ss){ W->bEpairs=true; }}; // add explicit electron pair
    funcs["-EachAngle"]={0,[&](const char** ss){ W->ffl.bEachAngle=true;                          }};
    funcs["-torsions"]={0,[&](const char** ss){ W->ffl.bTorsion=true; W->ffl.doPiPiI=false;  }};
    
    funcs["-substr_iso"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &app->subs_iso ); }};
    funcs["-iMO"       ]={1,[&](const char** ss){ sscanf( ss[0], "%i",  &app->which_MO ); }};
    funcs["-perframe"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &W->iterPerFrame ); app->perFrame=W->iterPerFrame; printf( "#### -perframe %i \n", W->iterPerFrame ); }};  // interations per frame

#endif // MolGUIapp_Lua



