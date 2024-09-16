#ifndef MolGUIapp_argv_h
#define MolGUIapp_argv_h

    funcs["-col_damp"]={6,[&](const char** ss){
        //printf( "ss[0](%s)\n", ss[0] ); printf( "ss[1](%s)\n", ss[1] );printf( "ss[2](%s)\n", ss[2] );printf( "ss[3](%s)\n", ss[3] );printf( "ss[4](%s)\n", ss[4] );exit(0);
        int n; double cB, cNB, cAng=0, cm,dR1,dR2;
        sscanf( ss[0], "%i" , &n   ) ; sscanf( ss[1], "%lf", &cB  ); sscanf( ss[2], "%lf", &cNB ); sscanf( ss[3], "%lf", &cm  ); sscanf( ss[3], "%lf", &dR1  ); sscanf( ss[3], "%lf", &dR2  ); 
        printf( "ARG W->ffl.ndampstep %i collisionDamping %g collisionDamping_NB %g damping_medium %g R1,2(%g,%g)\n", n, cB, cNB, cm, dR1, dR2 );
        W->ffl.colDamp.set( n, cm, cB, cAng, cNB, dR1, dR2 );
        printf( "ARG W->ffl.colDamp(n=%i,bond=%g,nonB=%g,medium=%g,R12(%g,%g)\n", W->ffl.colDamp.nstep, W->ffl.colDamp.bond, W->ffl.colDamp.nonB, W->ffl.colDamp.medium, W->ffl.colDamp.dRcut1, W->ffl.colDamp.dRcut2 );
    }};// collision damping parameters

    funcs["-T"]={2,[&](const char** ss){ sscanf( ss[0], "%lf", &W->go.T_target ); sscanf( ss[1], "%lf",&W->go.gamma_damp ); W->go.bExploring=true; W->bToCOG=true; }}; // run at non-zero temperature ( i.e. using Langevin dynamics termostat )

    funcs["-gopt"]={2,[&](const char** ss){ sscanf( ss[0], "%i,%i", &W->go.nExplore, &W->go.nRelax ); sscanf( ss[1], "%lf,%lf", &W->go.pos_kick, &W->go.vel_kick ); W->bGopt=true; W->bToCOG=true;  }}; // global optimization
    funcs["-drive"]={1,[&](const char** ss){ W->go.constrs.loadBonds(ss[0]); }}; // test

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
    funcs["-dlvec"]={1,[&](const char** ss){ Mat3d* m=new Mat3d(); W->dlvec=m; printf( "ARG ss[0] `%s`\n", ss[0] ); sscanf(ss[0],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); printf( "ARG W->dlvec set to " ); printMat(*(W->dlvec)); } }; // test
    funcs["-latscan"]={2,[&](const char** ss){ 
        W->bLatScan=true;
        Mat3d* m=new Mat3d(); W->latscan_dlvec=m; 
        printf( "ARG ss[0] `%s` ss[1] `%s`\n", ss[0], ss[1] );
        sscanf(ss[0],"%i,%i", &W->latscan_n.x, &W->latscan_n.y );
        sscanf(ss[1],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m->a.x,&m->a.y,&m->a.z,  &m->b.x,&m->b.y,&m->b.z,  &m->c.x,&m->c.y,&m->c.z ); 
        printf( "ARG W->latscan_n(%i,%i) latscan_dlvec ", W->latscan_n.x, W->latscan_n.y ); printMat(*(W->latscan_dlvec)); 
    } }; // test

    // set verbosity
    funcs["-verb"]={1,[&](const char** ss){ sscanf( ss[0], "%i", &verbosity ); }};

    funcs["-uff"]={0,[&](const char** ss){ W->bUFF=true; }}; // AutoCharge

    funcs["-e"]={0,[&](const char** ss){ W->bEpairs=true; }}; // add explicit electron pair
    funcs["-EachAngle"]={0,[&](const char** ss){ W->ffl.bEachAngle=true;                          }};
    funcs["-torsions"]={0,[&](const char** ss){ W->ffl.bTorsion=true; W->ffl.doPiPiI=false;  }};
    
    funcs["-substr_iso"]={1,[&](const char** ss){ sscanf( ss[0], "%lf", &app->subs_iso ); }};
    funcs["-iMO"       ]={1,[&](const char** ss){ sscanf( ss[0], "%i",  &app->which_MO ); }};
    funcs["-perframe"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &W->iterPerFrame ); app->perFrame=W->iterPerFrame; printf( "ARG -perframe %i \n", W->iterPerFrame ); }};  // interations per frame

    funcs["-Ftol"]={1,[&](const char** ss){ sscanf(ss[0],"%lf", &W->Ftol_default ); printf( "ARG -Ftol Ftol_default=%g\n", W->Ftol_default ); }};  
    funcs["-seed"]={1,[&](const char** ss){ int seed; sscanf(ss[0],"%i", &seed ); srand(seed); irand=0; }}; 
    funcs["-stuck"]={1,[&](const char** ss){ sscanf(ss[0],"%i,%i,%lf", &W->nStuckMax, &W->nStuckTrj, &W->RStuck );  W->bCheckStuck=true;  }};  
    funcs["-zspring"]={1,[&](const char** ss){ sscanf(ss[0],"%lf,%lf,%lf", &W->ConstrZ_xmin, &W->ConstrZ_l, &W->ConstrZ_k ); W->ConstrZ_l-=W->ConstrZ_xmin; W->bConstrZ=true;   printf( "ARG: ConstrZ_xmin %g ConstrZ_l %g ConstrZ_k %g\n", W->ConstrZ_xmin, W->ConstrZ_l, W->ConstrZ_k ); }};  
    funcs["-iParalel"]={1,[&](const char** ss){ sscanf(ss[0],"%i", &W->iParalel); printf( "ARG -iParalel %i \n", W->iParalel ); }};         // paralelization model
    funcs["-dt"]={1,[&](const char** ss){ sscanf(ss[0],"%lf", &W->dt_default); printf( "ARG -dt dt_default=%g \n", W->dt_default ); }};
    funcs["-NBneigh"]={0,[&](const char** ss){ W->bNonBondNeighs=true; }};
    funcs["-noNB"]={0,[&](const char** ss){ W->bNonBonded=false; }};
    funcs["-tricubic"]={0,[&](const char** ss){ W->bTricubic=true; }};
    funcs["-bbox"]={1,[&](const char** ss){ Mat3d m;  sscanf(ss[0],"%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",  &m.ax,&m.ay,&m.az,  &m.bx,&m.by,&m.bz,  &m.cx,&m.cy,&m.cz );  W->bbox=m; }};

    funcs["-nogridff"]={0,[&](const char** ss){ W->bGridFF=false; }}; // AutoCharge
    funcs["-group"]={3,[&](const char** ss){  }};
    
    funcs2["-c"]={1,[&](const char** ss){ 
        //printf( "MolGUIapp_argv.h :: [-c] `%s`\n", ss[0] );
        int ic=-1; Quat4d k=Quat4d{W->Kfix,W->Kfix,W->Kfix,0};  
        int nret = sscanf( ss[0], "%i,%lf,%lf,%lf", &ic, &k.x,&k.y,&k.z ); 
        if(nret==2){ k.y=k.x; k.z=k.x; }; 
        if( (ic<0)||(ic>W->ffl.natoms) ){  printf( "ERROR in MolGUIapp_argv.h `-c` ic(%i) out of range 0 .. natoms(%i) => exit() \n", ic, W->ffl.natoms ); exit(0); };
        W->ffl.constr [ic].w=1; 
        W->ffl.constrK[ic]=k;
        }}; //
    funcs2["-nogridff"]={0,[&](const char** ss){ W->bGridFF=false; }}; // AutoCharge
    funcs2["-group"]={3,[&](const char** ss){  
        printf( "MolGUIapp_argv.h :: [-group] \n" );
        W->bGroups = true;
        //int ig;               sscanf(ss[0], "%i", &ig);
        int i0;Vec2i ifw,iup; sscanf(ss[0], "%i,%i,%i,%i,%i", &ifw.x,&ifw.y,&iup.x,&iup.y, &i0);  printf( "MolGUIapp_argv.h::GROUPS ifw(%i,%i),iup(%i,%i) i0=%i \n", ifw.x,ifw.y,  iup.x,iup.y, i0 );
        int a2g[16]; int z; int   ng=sscanf(ss[1], "%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i,%i", &z,a2g+1,a2g+2,a2g+3,a2g+4,a2g+5,a2g+6,a2g+7,a2g+8,a2g+9,a2g+10,a2g+11,a2g+12,a2g+13,a2g+14,a2g+15);
        Vec2f fDrive;         sscanf(ss[2], "%lf,%lf", &fDrive.x,&fDrive.y );
        //printf( "--------- group[%i] n=%i\n", ig, ng );
        a2g[0]=z;
        if(W->atom2group.size()==0){ W->atom2group.resize( W->ffl.natoms ); for(int i=0;i<W->ffl.natoms;i++){ W->atom2group[i]=-1; } }

        // if( W->groups.a2g == 0 ){
        //     W->groups.bindAtoms(W->ffl.apos, W->ffl.fapos);
        //     printf( "MolGUIapp_argv.h::GROUPS W->groups.realloc(%i)\n", W->ffl.natoms );
        //     W->groups.realloc(W->ffl.natoms);

        // }
        int ig = W->groups.addGroup( ng, a2g );
        printf( "MolGUIapp_argv.h::GROUPS W->groups.addGroup(%i) \n", ig );
        W->groups.defPoseByAtoms( ig, ifw, iup, i0 );
        //W->groups.groups[ig].torq = Vec3d{1.0,0.0,0.0};   
        W->groups.groups[ig].torq = Vec3d{0.0,0.0,0.0};   
        W->groups.groups[ig].fDrive = fDrive;   
        //W->groups.evalRot(ig);
        for(int i=0; i<ng;i++){ W->atom2group[a2g[i]]=ig;}   
        //exit(0);
    }}; // AutoCharge
    


#endif // MolGUIapp_Lua



