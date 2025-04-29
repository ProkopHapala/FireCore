
#include  "globals.h"

#include "testUtils.h"
#include "MolWorld_sp3.h"

// ============ Global Variables

MolWorld_sp3 W;



//============================

#include "libMMFF.h"
#include "libUtils.h"


GridShape grid;
double* grid_data=0;

extern "C"{

void get_gridFF_info( int* int_data, double* float_data ){
    int_data[0]=W.gridFF.natoms;
    // int_data[1]=W.gridFF.natoms.size();
    // Add debug print here
    // printf("DEBUG: W.gridFF.shift0 = (%g, %g, %g)\n", W.gridFF.shift0.x, W.gridFF.shift0.y, W.gridFF.shift0.z);
    // printf("DEBUG: W.gridFF.grid.pos0 = (%g, %g, %g)\n", W.gridFF.grid.cell.ax, W.gridFF.grid.cell.ay, W.gridFF.grid.cell.az);
    // printf("DEBUG: W.gridFF.grid.new_voxels = (%g, %g, %g)\n", W.gridFF.grid.new_voxels.x, W.gridFF.grid.new_voxels.y, W.gridFF.grid.new_voxels.z);

    *((Vec3d*)(float_data  ))=W.gridFF.shift0;
    *((Vec3d*)(float_data+3))=W.gridFF.grid.pos0;
    *((Mat3d*)(float_data + 6))=W.gridFF.grid.cell;
    *((Mat3d*)(float_data + 15)) =W.gridFF.grid.dCell;
    

}
void get_atom_positions(double* apos_substrate, double* apos_molecule) {
    // Substrate atom positions
    if (apos_substrate != nullptr) {
        for (int i = 0; i < W.surf.natoms; ++i) {
            ((Vec3d*)apos_substrate)[i] = W.surf.apos[i];
        }
    }

    // Molecule atom positions
    if (apos_molecule != nullptr) {
        for (int i = 0; i < W.nbmol.natoms; ++i) {
            ((Vec3d*)apos_molecule)[i] = W.nbmol.apos[i];
        }
    }

    
}

int get_molecule_natoms() {
    return W.nbmol.natoms;
}

void init_buffers(){
    //printf( "init_buffers() \n" );

    // buffers .insert( { "gridff_apos",   (double*)W.gridFF.apos    } );
    // // buffers .insert( { "gridff_apos_",   (double*)W.gridFF.apos_  } );
    //  buffers.insert(std::make_pair(std::string("gridff_apos_"), (double*)W.gridFF.apos_.data()));
    // buffers .insert( { "gridff_Bspline_PLQ",   (double*)W.gridFF.Bspline_PLQ  } );


    buffers .insert( { "apos",   (double*)W.nbmol.apos  } );
    buffers .insert( { "fapos",  (double*)W.nbmol.fapos } );
    buffers .insert( { "REQs",   (double*)W.nbmol.REQs  } );
    if (W.bUFF){
        buffers .insert( { "DOFs",   (double*)   W.ffu.apos  } );
        buffers .insert( { "fDOFs",  (double*)   W.ffu.fapos } );
        buffers .insert( { "vDOFs",  (double*)   W.ffu.vapos   } );
        //buffers .insert( { "REQs",   (double*)W.ffl.REQs  } );
        buffers .insert( { "PLQs",   (double*)W.nbmol.PLQd  } );
        //buffers .insert( { "pipos",  (double*)W.ffl.pipos   } );
        //buffers .insert( { "fpipos", (double*)W.ffl.fpipos } );
        //ibuffers.insert( { "neighs",      (int*)W.ffl.neighs  } );

        ibuffers.insert( { "ndims",    &W.ffu.natoms } );
        buffers .insert( { "Es",       &W.ffu.Etot   } );

    }else if(W.bMMFF){
        ibuffers.insert( { "ndims",    &W.ff.nDOFs } );
        buffers .insert( { "Es",       &W.ff.Etot  } );

        buffers .insert( { "DOFs",      W.ffl.DOFs  } );
        buffers .insert( { "fDOFs",     W.ffl.fDOFs } );
        buffers .insert( { "vDOFs",     W.opt.vel  } );
        //buffers .insert( { "REQs",   (double*)W.ffl.REQs  } );
        buffers .insert( { "PLQs",   (double*)W.ffl.PLQd  } );
        if(!W.bUFF){
            buffers .insert( { "pipos",  (double*)W.ffl.pipos   } );
            buffers .insert( { "fpipos", (double*)W.ffl.fpipos }  );
            ibuffers.insert( { "neighs", (int*)W.ffl.neighs  }    );
        }
    } else{
        W.ff.natoms=W.nbmol.natoms;
    }
    printf( "before_MMFF_lib.cpp::init_buffers() ndims{nDOFs=%i,natoms=%i,nnode=%i,ncap=%i}\n", W.ffl.nDOFs, W.ffl.natoms, W.ffl.nnode, W.ffl.ncap );
    printf( "MMFF_lib.cpp::init_buffers() ndims{nDOFs=%i,natoms=%i,nnode=%i,ncap=%i,npi=%i,nbonds=%i,nvecs=%i,ne=%i,ie0=%i}\n", W.ff.nDOFs, W.ff.natoms, W.ff.nnode, W.ff.ncap, W.ff.npi, W.ff.nbonds, W.ff.nvecs, W.ff.ne, W.ff.ie0 );

    ibuffers.insert( { "selection", W.manipulation_sel  } );
    bbuffers.insert( { "ffflags", &W.doBonded  } );
    //printBuffNames();
}

// int loadmol(char* fname_mol ){ return W.loadmol(fname_mol ); }
//lib.init( cstr(xyz_name), cstr(surf_name), cstr(smile_name),      bMMFF,      bEpairs,      bUFF,      b141,      bSimple,      bConj,      bCumulene,      nPBC,        gridStep, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )
void* init( char* xyz_name, char* surf_name, char* smile_name, bool bMMFF, bool bEpairs, bool bUFF, bool b141, bool bSimple, bool bConj, bool bCumulene, int* nPBC, double gridStep, char* sElementTypes, char* sAtomTypes, char* sBondTypes, char* sAngleTypes, char* sDihedralTypes ){
	W.smile_name = smile_name;
	W.xyz_name   = xyz_name;
	W.surf_name  = surf_name;
	W.bMMFF      = bMMFF;
    W.bEpairs    = bEpairs;
    W.gridStep   = gridStep;
    W.nPBC       = *(Vec3i*)nPBC;
    W.bUFF       = bUFF; 
    W.b141       = b141;
    W.bSimple    = bSimple;
    W.bConj      = bConj;
    W.bCumulene  = bCumulene;
    // read and store parameters from tables
    // TBD pass bUFF to MMFFparams::init so that if true, no need to read bonds, angles nor dihedrals...
    //W.params.verbosity = verbosity;
    //W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
    // bring names of atom types into builder (H is capping atom, E is electron pair)
	//W.builder.bindParams(&W.params);
    W.initParams( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
    bool bGrid = gridStep>0;
    // initialize the main
    //W.init( bGrid, bUFF );
    W.bGridFF=bGrid;
    W.bUFF   =bUFF;
    W.init();
    init_buffers();
    return &W;
}

void makeGridFF( const char* name, int* ffshape, int mode, double z0, double* cel0, bool bSymmetrize, bool bAutoNPBC, bool bFit, bool bRefine ){
    bool bCheckEval=false;
    bool bUseEwald =true;
    printf("MMFF_lib::makeGridFF() bAutoNPBC=%i bCheckEval=%i bUseEwald=%i bFit=%i bRefine=%i \n", bAutoNPBC, bCheckEval, bUseEwald, bFit, bRefine );
    char fname[256];
    sprintf(fname, "%s.xyz", name );
    int ret = W.params.loadXYZ( fname, W.surf.natoms, &W.surf.apos, &W.surf.REQs, &W.surf.atypes, 0, &W.gridFF.grid.cell );
    if     ( ret<0 ){ getcwd(tmpstr,1024); printf("ERROR in MMFF_lib::makeGridFF() file(%s) not found in path(%s)=> Exit() \n", fname, tmpstr ); exit(0); }
    if     ( ret==0){                      printf("ERROR in MMFF_lib::makeGridFF() no lattice vectors in (%s) => Exit() \n",    fname ); exit(0); }
    else if( ret>0 ){ W.gridFF.grid.updateCell(W.gridStep); W.gridFF.bCellSet=true;  }
    //gridFF.grid.printCell(); 
    //if(verbosity>0)printf("MolWorld_sp3::loadSurf(%s) 1 natoms %i apos %li atyps %li \n", name, surf.natoms, (long)surf.apos, (long)surf.atypes  );
    //surf.print();
    //W.surf.print_nonbonded();
    W.gridFF.mode=(GridFFmod)mode;
    W.bSurfAtoms=true;
    //printf("MMFF_lib::makeGridFF() bAutoNPBC=%i bCheckEval=%i bUseEwald=%i bFit=%i bRefine=%i \n", bAutoNPBC, bCheckEval, bUseEwald, bFit, bRefine );
    W.initGridFF( name, z0, *(Vec3d*)cel0, bSymmetrize, bAutoNPBC, bCheckEval, bUseEwald, bFit, bRefine );
    ffshape[0]=W.gridFF.grid.n.x;
    ffshape[1]=W.gridFF.grid.n.y;
    ffshape[2]=W.gridFF.grid.n.z;
    ffshape[3]=W.gridFF.perVoxel;
    //return ff_ptr;
}

double* getArrayPointer( const char* name, int* shape  ){
    if(golbal_array_dict.find(name)!=golbal_array_dict.end()){
        NDArray arr = golbal_array_dict[name];
        (*(Quat4i*)shape) = arr.dims;  
        if(arr.data==0){ printf("ERROR in MMFF_lib::getArrayPointer() golbal_array_dict[%s].data==NULL \n", name ); }
        return arr.data;
    }else{
        printf("ERROR in MMFF_lib::getArrayPointer() golbal_array_dict[%s] not found \n", name );
        //exit(0);
    }
    return 0;
}

int setupEwaldGrid( double* pos0, double* dCell, int* ns, bool bPrint ){
    W.gewald.n     = *(Vec3i*)ns;
    W.gewald.pos0  = *(Vec3d*)pos0;
    W.gewald.dCell = *(Mat3d*)dCell;
    W.gewald.updateCell_2();
    if(bPrint){W.gewald.printCell();}
    return W.gewald.n.totprod();
}

void projectAtomsEwaldGrid( int na, double* apos, double* qs, double* dens, int order ){
    W.gewald.projectAtoms( na, (Vec3d*)apos, qs, dens, order );
}


#ifdef WITH_FFTW

void EwaldGridSolveLaplace( double* dens, int nz_slab, double* Vout, bool bPrepare, bool bDestroy, int flags, bool bOMP, int nBlur, double cSOR, double cV ){
    W.gewald.solve_laplace_macro( dens, nz_slab, Vout, bPrepare, bDestroy, flags, bOMP, nBlur, cSOR, cV );
}

void EwaldGridSolveLaplaceDebug( double* dens, double* Vout, double* densw, double* kerw, double* VwKer ){
    int ntot = W.gewald.n.totprod();
    
    W.gewald.prepare_laplace( );
    array2fftc( ntot, dens, W.gewald.V );
    fftw_execute(W.gewald.fft_plan);   fftc2array( ntot, W.gewald.Vw,  densw  );
    
    for(int i=0; i<ntot; i++){  W.gewald.V[i][0]=1.0; W.gewald.V[i][1]=1.0; }
    W.gewald.laplace_reciprocal_kernel( W.gewald.V  );  fftc2array( ntot, W.gewald.V,  kerw  );
    W.gewald.laplace_reciprocal_kernel( W.gewald.Vw );  fftc2array( ntot, W.gewald.Vw, VwKer );
    fftw_execute(W.gewald.ifft_plan);                   fftc2array( ntot, W.gewald.V,  Vout  );

    W.gewald.destroy_laplace( );
}

#endif // WITH_FFTW

void evalGridFFAtPoints( int n, double* ps, double* FFout, double* PLQH, bool bSplit, int* nPBC ){
    //long t0 = getCPUticks();
    if(bSplit){ W.gridFF.evalAtPoints_Split( n, (Vec3d*)ps, (Quat4d*)FFout, *(Quat4d*)PLQH, (Vec3i*)nPBC ); }
    else      { 
        //W.gridFF.evalAtPoints      ( n, (Vec3d*)ps, (Quat4d*)FFout, *(Quat4d*)PLQH, (Vec3i*)nPBC ); 
        W.gridFF.evalAtPoints_REQ  ( n, (Vec3d*)ps, (Quat4d*)FFout, *(Quat4d*)PLQH, (Vec3i*)nPBC ); 
    }
    //double T = getCPUticks()-t0; printf( "evalGridFFAtPoints(n=%i,bSplit=%i) DONE in %g[MTicks] %g[kTick/point] \n", n, bSplit, T*1e-6, (T*1e-3)/n  );
}

int    run( int nstepMax, double dt, double Fconv, int ialg, double damping, double* outE, double* outF, double* outV, double* outVF, bool omp ){
    //printf( "bOpenMP = %i \n", omp );
    //W.rum_omp_ocl( nstepMax, dt, Fconv, 1000.0, 1000 ); 
    // run_omp( int niter_max, double dt, double Fconv=1e-6, double Flim=1000, double timeLimit=0.02, double* outE=0, double* outF=0 ){
    if(omp){ return W.run_omp   (nstepMax,dt,Fconv,   10.0, -1.0, outE, outF, outV, outVF ); }
    else   { return W.run_no_omp(nstepMax,dt,Fconv, 1000.0,  damping, outE, outF, outV, outVF ); }
    //else   { return W.run       (nstepMax,dt,Fconv,ialg,       outE, outF, outV, outVF ); }
}

void  scan( int nconf, double* poss, double* rots, double* Es, double* aforces, double* aposs, bool omp, bool bRelax, int niter_max, double dt, double Fconv, double Flim ){
    // Add debug print here
    // printf("DEBUG: scan() function using shift0 = (%g, %g, %g)\n", W.gridFF.shift0.x, W.gridFF.shift0.y, W.gridFF.shift0.z);
    // --- DEBUG PRINT (Added previously) ---
    printf("DEBUG>> MMFF_lib: AFTER setSwitches2: W.ffl.PLQd[0].z = %.6f\n", W.ffl.PLQd[0].z); 
    // --- END DEBUG ---

    if(bRelax){
        if(omp){ printf("ERROR: scan_relaxed() not implemented witht OMP\n"); exit(0); } 
        else   { W.scan_relaxed( nconf, (Vec3d*)poss, (Mat3d*)rots, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp, niter_max, dt, Fconv, Flim );  }
    }else{
        if(omp){ printf("ERROR: scan_rigid() not implemented witht OMP\n"); exit(0); } 
        else   { W.scan_rigid( nconf, (Vec3d*)poss, (Mat3d*)rots, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp ); }
    }
}

// void scan_relaxed_constr( int nconf, int ncontr, int *icontrs, Quat4d* contrs, double* Es, Vec3d* aforces, Vec3d* aposs, bool bHardConstr, bool omp, int niter_max, double dt, double Fconv=1e-6, double Flim=1000 ){

void  scan_constr( int nconf, int ncontr, int *icontrs, double* contrs, double* Es, double* aforces, double* aposs, bool bHardConstr, bool bOmp, int niter_max, double dt, double Fconv, double Flim ){
    W.scan_constr( nconf, ncontr, icontrs, (Quat4d*)contrs, Es, (Vec3d*)aforces, (Vec3d*)aposs, bHardConstr, bOmp, niter_max, dt, Fconv, Flim);
}


/*void setSwitches_testUFF( int bBonds, int bAngles, int bDihedrals, int bInversions, int bNonBonded, int bNonBondNeighs, int bSubtractBondNonBond, int bSubtractAngleNonBond, int bClampNonBonded ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.ffu.bBonds                , bBonds                );
    _setbool( W.ffu.bAngles               , bAngles               );
    _setbool( W.ffu.bDihedrals            , bDihedrals            );
    _setbool( W.ffu.bInversions           , bInversions           );
    _setbool( W.ffu.bNonBonded            , bNonBonded            );
    _setbool( W.ffu.bNonBondNeighs        , bNonBondNeighs        );
    _setbool( W.ffu.bSubtractBondNonBond  , bSubtractBondNonBond  );
    _setbool( W.ffu.bSubtractAngleNonBond , bSubtractAngleNonBond );
    _setbool( W.ffu.bClampNonBonded       , bClampNonBonded       );
    printf( "setSwitches_testUFF() bBonds=%i bAngles=%i bDihedrals=%i bInversions=%i | bNonBonded=%i bNonBondNeighs=%i bSubtractBondNonBond=%i bSubtractAngleNonBond=%i bClampNonBonded=%i \n", W.ffu.bBonds, W.ffu.bAngles, W.ffu.bDihedrals, W.ffu.bInversions, W.ffu.bNonBonded, W.ffu.bNonBondNeighs, W.ffu.bSubtractBondNonBond, W.ffu.bSubtractAngleNonBond, W.ffu.bClampNonBonded );
    #undef _setbool
}*/

void test_UFF( int test ){
    W.ffu.bBonds                = 0; // consider bonds
    W.ffu.bAngles               = 0; // consider angles
    W.ffu.bDihedrals            = 0; // consider dihedrals
    W.ffu.bInversions           = 0; // consider impropers
    W.ffu.bNonBonded            = 0; // consider intra-molecular LJ (general)
    W.ffu.bNonBondNeighs        = 0; // consider intra-molecular LJ and does not compute 1-2 non-bonded interactions (ng4)
    W.ffu.bSubtractBondNonBond  = 0; // subtract 1-2 non-bonded interactions when computing bonds
    W.ffu.bSubtractAngleNonBond = 0; // subtract 1-3 non-bonded interactions when computing angles
    W.ffu.bClampNonBonded       = 0; // clamp forces before subtracting 
    W.ffu.bSubstrate            = 0; // consider non-bonded interaction with the substrate
    if     (test== 1){ printf("MMFF_lib.cpp::test_UFF(1): bonds\n");      
        W.ffu.bBonds = 1; }
    else if(test== 2){ printf("MMFF_lib.cpp::test_UFF(2): angles\n");     
        W.ffu.bAngles = 1; }
    else if(test== 3){ printf("MMFF_lib.cpp::test_UFF(3): dihedrals\n");  
        W.ffu.bDihedrals = 1; }
    else if(test== 4){ printf("MMFF_lib.cpp::test_UFF(4): inversions\n"); 
        W.ffu.bInversions = 1; }
    else if(test== 5){ printf("MMFF_lib.cpp::test_UFF(5): lj\n"); 
        W.ffu.bNonBonded = 1; }
    else if(test== 6){ printf("MMFF_lib.cpp::test_UFF(6): bonds+lj_excluding_1-2_subtract_noclamp\n"); 
        W.ffu.bBonds = 1; 
        W.ffu.bNonBonded = 1; 
        W.ffu.bSubtractBondNonBond = 1; }
    else if(test== 7){ printf("MMFF_lib.cpp::test_UFF(7): bonds+lj_excluding_1-2_subtract_clamp\n"); 
        W.ffu.bBonds = 1; 
        W.ffu.bNonBonded = 1; 
        W.ffu.bSubtractBondNonBond = 1; 
        W.ffu.bClampNonBonded = 1; }
    else if(test== 8){ printf("MMFF_lib.cpp::test_UFF(8): bonds+lj_excluding_1-2_ng4\n"); 
        W.ffu.bBonds = 1; 
        W.ffu.bNonBonded = 1;
        W.ffu.bNonBondNeighs = 1; }
    else if(test==9){ printf("MMFF_lib.cpp::test_UFF(9): all_intra_subtract_noclamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bSubtractBondNonBond = 1;
        W.ffu.bSubtractAngleNonBond = 1; }
    else if(test==10){ printf("MMFF_lib.cpp::test_UFF(10): all_intra_subtract_clamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bSubtractBondNonBond = 1;
        W.ffu.bSubtractAngleNonBond = 1;
        W.ffu.bClampNonBonded = 1; }
    else if(test==11){ printf("MMFF_lib.cpp::test_UFF(11): all_intra_ng4_noclamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bNonBondNeighs = 1; 
        W.ffu.bSubtractAngleNonBond = 1; }
    else if(test==12){ printf("MMFF_lib.cpp::test_UFF(12): all_intra_ng4_clamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bNonBondNeighs = 1; 
        W.ffu.bSubtractAngleNonBond = 1;
        W.ffu.bClampNonBonded = 1; }
    else if(test==13){ printf("MMFF_lib.cpp::test_UFF(13): substrate\n"); 
        W.ffu.bSubstrate = 1; }
    else if(test==14){ printf("MMFF_lib.cpp::test_UFF(14): all_subtract_noclamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bSubtractBondNonBond = 1;
        W.ffu.bSubtractAngleNonBond = 1;
        W.ffu.bSubstrate = 1; }
    else if(test==15){ printf("MMFF_lib.cpp::test_UFF(15): all_subtract_clamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bSubtractBondNonBond = 1;
        W.ffu.bSubtractAngleNonBond = 1;
        W.ffu.bClampNonBonded = 1;
        W.ffu.bSubstrate = 1; }
    else if(test==16){ printf("MMFF_lib.cpp::test_UFF(16): all_ng4_noclamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bNonBondNeighs = 1; 
        W.ffu.bSubtractAngleNonBond = 1;
        W.ffu.bSubstrate = 1; }
    else if(test==17){ printf("MMFF_lib.cpp::test_UFF(17): all_ng4_clamp\n"); 
        W.ffu.bBonds = 1;
        W.ffu.bAngles = 1;
        W.ffu.bDihedrals = 1;
        W.ffu.bInversions = 1;
        W.ffu.bNonBonded = 1;
        W.ffu.bNonBondNeighs = 1; 
        W.ffu.bSubtractAngleNonBond = 1;
        W.ffu.bClampNonBonded = 1;
        W.ffu.bSubstrate = 1; }
    W.ffu.test_UFF();
}


void setSwitches2( int CheckInvariants, int PBC, int NonBonded, int NonBondNeighs,  int SurfAtoms, int GridFF, int MMFF, int Angles, int PiSigma, int PiPiI ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.bCheckInvariants, CheckInvariants  );
    _setbool( W.bPBC           , PBC       );
    
    _setbool( W.bNonBonded     , NonBonded );
    _setbool( W.bNonBondNeighs , NonBondNeighs );
    
    _setbool( W.bSurfAtoms   , SurfAtoms );
    _setbool( W.bGridFF      , GridFF    );

    _setbool( W.bMMFF        , MMFF      );
    _setbool( W.ffl.doAngles , Angles    );
    _setbool( W.ffl.doPiSigma, PiSigma   );
    _setbool( W.ffl.doPiPiI  , PiPiI     );

    printf( "setSwitches2() W.bCheckInvariants==%i bPBC=%i | bNonBonded=%i bNonBondNeighs=%i | bSurfAtoms=%i bGridFF=%i | bMMFF=%i doAngles=%i doPiSigma=%i doPiPiI=%i \n", W.bCheckInvariants, W.bPBC,  W.bNonBonded, W.bNonBondNeighs, W.bSurfAtoms, W.bGridFF, W.bMMFF, W.ffl.doAngles, W.ffl.doPiSigma, W.ffl.doPiPiI );
    // --- DEBUG PRINT ---
    printf("DEBUG>> MMFF_lib: AFTER setSwitches2: W.ffl.PLQd[0].z = %.6f\n", W.ffl.PLQd[0].z);
    // --- END DEBUG ---
    //W.ffl.bSubtractAngleNonBond = W.bNonBonded;
    #undef _setbool
}



int substituteMolecule( const char* fname, int ib, double* up, int ipivot, bool bSwapBond ){
    return W.substituteMolecule( fname, ib, *(Vec3d*)up, ipivot, bSwapBond );
}

void set_opt( 
        double dt_max,  double dt_min, double damp_max, 
        double finc,    double fdec,   double falpha, int minLastNeg,
        double cvf_min, double cvf_max
    ){

    W.opt.dt_max  = dt_max;
    W.opt.dt_min  = dt_min;
    W.opt.dt      = dt_max;

    W.opt.damp_max   = damp_max;
    W.opt.damping    = damp_max;

    W.opt.cvf_min    = cvf_min;
    W.opt.cvf_max    = cvf_max;

    W.opt.minLastNeg =  minLastNeg;
    W.opt.finc       =  finc;
    W.opt.fdec       =  fdec;
    W.opt.falpha     =  falpha;
    
    //W.opt.f_limit  =  f_limit  ;
    //W.opt.v_limit  =  v_limit  ;
    //W.opt.dr_limit =  dr_limit ;

}

double* setupGrid( int* ns, double* cell, bool bAlloc ){
    if( ns     ){ grid.n    = *((Vec3i*)ns);   }
    if( cell   ){ grid.cell = *((Mat3d*)cell); }
    if( bAlloc ) _realloc( grid_data, grid.getNtot() );
     return grid_data;
}

int loadBin_d( const char* fname, double* data ){
    return loadBin( fname, grid.getNtot() * sizeof(double), (char*)data );
}

int loadBin_f( const char* fname, float* data ){
    return loadBin( fname, grid.getNtot() * sizeof(float), (char*)data );
}


int saveBin_d( const char* fname, double* data ){
    return saveBin( fname, grid.getNtot() * sizeof(double), (char*)data );
}

double* loadXSF( const char* fname, int* ns, double* cell ){
    grid_data = grid.loadXSF<double>( fname, 0 );
    if(ns  ){ *((Vec3i*)ns)   = grid.n;    }
    if(cell){ *((Mat3d*)cell) = grid.cell; }
    return grid_data;
}

void saveXSF( const char* fname, const double* data, int* ns, double* cell ){
    if(ns  ){ grid.n    = *((Vec3i*)ns);    }
    if(cell){ grid.cell = *((Mat3d*)cell); }
    if(data==0){  data=grid_data; }
    grid.saveXSF<double>( fname, data );
}

void sampleSurf(char* name, int n, double* rs, double* Es, double* fs, int kind, int atyp, double Q, double K, double RQ, double* pos0_, bool bSave){
    if(name){
        W.ff.realloc( 1,0,0,0, true );
        W.ff.apos [0] = *(Vec3d*)pos0_;
        W.ff.atype[0] = atyp;
        bool bGrid=(kind>10);
        if( kind==10 ) W.gridFF.iDebugEvalR=1;
        W.gridFF.alphaMorse = K;
        W.gridFF.Rdamp = RQ;
        W.loadSurf( name, bGrid, bSave );
        W.nbmol.REQs[0].z = Q;
        if(bSave){
            Quat4f* FFtot = new Quat4f[W.gridFF.grid.getNtot()];
            W.gridFF.evalCombindGridFF ( W.nbmol.REQs[0], FFtot );
            W.gridFF.grid.saveXSF<float>( "FFtot_E.xsf", (float*)FFtot, 4, 3, W.surf.natoms, W.surf.atypes, W.surf.apos );
            delete [] FFtot;
        }
    }
    Quat4d REQ=W.nbmol.REQs[0];
    Quat4f PLQ = REQ2PLQ( REQ, K );
    Quat4d PLQd= REQ2PLQ_d( REQ, K );
    printf( "DEBUG sampleSurf REQ(%g,%g,%g) \n", REQ.x, REQ.y, REQ.z );
    printf( "DEBUG sampleSurf PLQ(%g,%g,%g) \n", PLQ.x, PLQ.y, PLQ.z );
    //exit(0);
    double R2Q=RQ*RQ;
    for(int i=0; i<n; i++){
        Quat4f fe=Quat4fZero;
        //Quat4d fed=Quat4dZero;
        W.nbmol.apos[0].z=rs[i];
        W.ff.cleanAtomForce();
        switch(kind){
            case  0: fe.e=   W.nbmol.evalR         (W.surf ); break; 
            case  1: fe.e=   W.nbmol.evalMorse     (W.surf, false,                           K,RQ  ); fe.f=(Vec3f)W.nbmol.fapos[0]; break; 
            //case  5: fe.e=   W.nbmol.evalMorsePLQ  (W.surf, PLQ, W.gridFF.grid.cell, {1,1,0},K,R2Q ); fe.f=(Vec3f)W.nbmol.fapos[0]; break; 
            case 10:         W.gridFF.addForce_surf(W.nbmol.apos[0], {1.,0.,0.}, fe );  break;
            case 11:         W.gridFF.addForce_surf(W.nbmol.apos[0], PLQ, fe );  break;
            case 12:         W.gridFF.addForce     (W.nbmol.apos[0], PLQ, fe );  break;

            //case 13:         W.gridFF.addForce_surf(W.nbmol.apos[0], {1.,0.,0.}, fe );  break;

            case 13: fe = (Quat4f) W.gridFF.getForce_HHermit( W.nbmol.apos[0], PLQd );   break;
            case 14: fe = (Quat4f) W.gridFF.getForce_Bspline( W.nbmol.apos[0], PLQd );    break;  

        }
        fs[i]=fe.z;
        Es[i]=fe.e;
    }
}


void sampleSurf_new( int n, double* ps_, double* FEout_, int mode, double* PLQH_, double K, double RQ ){
    Vec3d*  ps   =((Vec3d*)ps_);
    Quat4d* FEout=(Quat4d*)FEout_;
    Quat4d  PLQd = *(Quat4d*) PLQH_;
    Quat4f  PLQ  =  (Quat4f) PLQd;
    double R2Q=RQ*RQ;

    printf( "sampleSurf_new() n=%i mode=%i PLQH(%g,%g,%g,%g) K=%g RQ=%g \n", n, mode, PLQd.x, PLQd.y, PLQd.z, PLQd.w, K, RQ  );
    //W.gridFF.grid.printCell();
    //printf( "sampleSurf_new() gff.shift0(%g,%g,%g) gff.pos0(%g,%g,%g)\n", W.gridFF.shift0.x, W.gridFF.shift0.y, W.gridFF.shift0.z, W.gridFF.grid.pos0.x, W.gridFF.grid.pos0.y, W.gridFF.grid.pos0.z );
    //printf( "sampleSurf_new() gridN(%i,%i,%i) \n", W.gridFF.gridN.x, W.gridFF.gridN.y, W.gridFF.gridN.z );

    //PLQd=Quat4d{1.0,0.0,0.0,0.0};
    // {
    //     Vec3i ng = W.gridFF.grid.n;
    //     Vec3d g0 = W.gridFF.grid.pos0;
    //     Vec3d dg = Vec3d{ W.gridFF.grid.dCell.xx, W.gridFF.grid.dCell.yy, W.gridFF.grid.dCell.zz };
    //     Quat4d C = PLQd;
    //     Quat4i* xqs = W.gridFF.cubic_xqis;
    //     printf("CPU sampleSurf_new() ng(%i,%i,%i) g0(%g,%g,%g) dg(%g,%g,%g) C(%g,%g,%g) \n", ng.x,ng.y,ng.z,   g0.x,g0.y,g0.z,   dg.x,dg.y,dg.z,   C.x,C.y,C.z );
    //     printf("CPU sampleSurf_new() xqs[0](%i,%i,%i,%i) xqs[1](%i,%i,%i,%i) xqs[2](%i,%i,%i,%i) xqs[3](%i,%i,%i,%i)\n", xqs[0].x, xqs[0].y, xqs[0].z, xqs[0].w,   xqs[1].x, xqs[1].y, xqs[1].z, xqs[1].w,   xqs[2].x, xqs[2].y, xqs[2].z, xqs[2].w,  xqs[3].x, xqs[3].y, xqs[3].z, xqs[3].w   );
    // }

    { // debug
        for(int i=0; i<4; i++ ){ 
            Quat4i xq = W.gridFF.cubic_xqis[i];
            Quat4i yq = W.gridFF.cubic_yqis[i];
            printf("sampleSurf_new() gridFF qs[%] xs(%3i,%3i,%3i,%3i) ys(%3i,%3i,%3i,%3i) \n", i, xq.x, xq.y, xq.z, xq.w,   yq.x, yq.y, yq.z, yq.w ); 
        }
    }
    
    //long t0 = getCPUticks();
    for(int i=0; i<n; i++){
        Quat4f fef=Quat4fZero;
        Quat4d fed=Quat4dZero;
        Vec3d pi = ps[i];
        switch(mode){
            case 1:   fef = W.gridFF.getForce( pi, PLQ );    fed=(Quat4d)fef; break;
            case 2:   fed = W.gridFF.getForce_d( pi, PLQd );       break;
            case 4:   fed = W.gridFF.getForce_HHermit( pi, PLQd ); break;
            case 6:   fed = W.gridFF.getForce_Bspline( pi, PLQd ); break;  
        }
        FEout[i]= fed;
    }
    //double t = (getCPUticks()-t0); printf( "sampleSurf_new(mode=%i,n=%i) time=%g[kTick] %g[tick/point]\n", mode, n, t*(1.e-3), t/n );
}


int findHbonds( double Rcut, double Hcut, double angMax ){
    W.Hbonds.clear();
    W.findHbonds_PBC( Rcut, Hcut, angMax*deg2rad );
    return W.Hbonds.size();
}

int sampleHbond( int ib, int n, double* rs, double* Es, double* fs, int kind, double maskQ, double maskH, double K, double Rdamp, double dcomp, char* s ){
    int nb = W.Hbonds.size();
    if( (ib<0) || (ib>nb)){  return nb;}
    Vec3i b = W.Hbonds[ib];
    Quat4d REQi = W.ffl.REQs[b.x];
    Quat4d REQj = W.ffl.REQs[b.y];
    Quat4d REQij; combineREQ( REQi, REQj, REQij );
    //printf( "@s=%li\n", (long)s );
    //printf( "type_name[%i]=%s\n", b.x, W.params.atypes[ W.ffl.atypes[b.x]].name );
    sprintf( s, "%s[%i]-%s[%i] (%4.2f,%5.4f,%4.2f,%4.2f) (%4.2f,%5.4f,%4.2f,%4.2f)", W.params.atypes[ W.ffl.atypes[b.x]].name, b.x, W.params.atypes[ W.ffl.atypes[b.y]].name, b.y, REQi.x,REQi.y,REQi.z,REQi.w,  REQj.x,REQj.y,REQj.z,REQj.w  );
    REQij.z*=maskQ; // Mask Electrostatics
    REQij.w*=maskH; // Mask HBond
    Vec3d pi=Vec3dZero;
    Vec3d pj=Vec3dZero;
    double R2damp=Rdamp*Rdamp;
    Vec3d d{dcomp,0.0,0.0};
    for(int i=0; i<n; i++){
        double E;
        Vec3d  f=Vec3dZero;
        pj.x=rs[i];
        switch(kind){
            case 1:{ E = getLJQH( pj-pi, f, REQij, R2damp ); } break;
            case 2:{ 
                Vec3d fi=Vec3dZero;
                E  = getLJQH( pj-pi    , fi, REQij                       , R2damp );  f.add(fi);
                E += getLJQH( pj-pi+d  , fi, Quat4d{1.0,0.0,-REQij.z,0.0}, R2damp );  f.add(fi);
                E += getLJQH( pj-pi+d  , fi, Quat4d{1.0,0.0,-REQij.z,0.0}, R2damp );  f.add(fi);
                E += getLJQH( pj-pi+d+d, fi, Quat4d{1.0,0.0, REQij.z,0.0}, R2damp );  f.add(fi);
                } break;
            //case 1: E=addAtomicForceMorseQ( pj-pi, f, REQij.x, REQij.y, REQij.z, K, R2damp ); break;  // Morse
            //case 2: E=addAtomicForceLJQ   ( pj-pi, f, REQij );                                break;  // Lenard Jones
            //case 3: double fr; E=erfx_e6( pj.x, K, fr ); f.x=fr; break;  // gauss damped electrostatics
            //case 4: E=repulsion_R4( pj-pi, f, REQij.x-Rdamp, REQij.x, K );
        }
        //printf( "i %i r %g E %g f %g \n", i, pj.x, E, f.x );
        fs[i]=f.x;
        Es[i]=E;
    }
    return nb;
}

void print_debugs( bool bParams, bool bNeighs, bool bShifts ){
    W.ffl.printSizes();
    if( bParams ) W.ffl.printAtomParams();
    if( bNeighs ) W.ffl.printNeighs();
    if( bShifts ) W.ffl.print_pbc_shifts();
}

//void sampleSurf_vecs(char* name, int n, double* poss_, double* FEs_, int kind, int ityp, double RvdW, double EvdW, double Q, double K, double RQ, double* pos0_, int npbc, bool bSave){
void sampleSurf_vecs(int n, double* poss_, double* FEs_, int kind, int ityp, double RvdW, double EvdW, double Q, double K, double RQ, int npbc, bool bSave){    
    //printf( "MMFF_lib::sampleSurf_vecs() kind=%i n=%i dCell(%g,%g,%g)\n", kind, n, W.gridFF.grid.dCell.xx,W.gridFF.grid.dCell.yy,W.gridFF.grid.dCell.zz );
    Vec3i nPBC{npbc,npbc,0};
    Vec3d* poss =(Vec3d*)poss_;
    //Vec3d* fs =(Vec3d*)fs_;
    Quat4d* FEs = (Quat4d*)FEs_;
    Quat4d  test_REQ{ RvdW, sqrt(EvdW), Q }; 
    if( ityp>0 ){
        AtomType atyp = W.params.atypes[ityp];
        test_REQ.x = atyp.RvdW;        // UFF natural bond radius
        test_REQ.y = sqrt(atyp.EvdW);  // LJ distance parameter
        //test_REQ.z = atyp.Qbase;
        test_REQ.y = atyp.Hb;          // LJ energy parameter
    }
    // if(name){
    //     // W.ff.realloc( 1,0,0,0, true );
    //     // W.ff.apos [0] = *(Vec3d*)pos0_;
    //     // W.ff.atype[0] = atyp;
    //     // bool bGrid=(kind>=10);
    //     // if( kind==10 ) W.gridFF.iDebugEvalR=1;
    //     // W.gridFF.alphaMorse = K;
    //     // W.gridFF.Rdamp = RQ;
    //     // W.loadSurf( name, bGrid, bSave );
    //     // W.nbmol.REQs[0].z = Q;
    // }
    if(bSave){
        Quat4f* FFtot = new Quat4f[W.gridFF.grid.getNtot()];
        W.gridFF.evalCombindGridFF ( test_REQ, FFtot );
        W.gridFF.grid.saveXSF<float>( "FFtot_E.xsf", (float*)FFtot, 4, 3, W.surf.natoms, W.surf.atypes, W.surf.apos );
        printf( "saveXSF() DONE \n" );
        delete [] FFtot;
    }
    Quat4f PLQ   = REQ2PLQ  ( test_REQ, K );
    Quat4d PLQ_d = REQ2PLQ_d( test_REQ, K );
    // Quat4f PLQ   {0.0,0.0,1.0,0.0};
    // Quat4d PLQ_d {0.0,0.0,1.0,0.0};
    printf( "MMFF_lib::sampleSurf_vecs() kind=%3i n=%6i dCell(%g,%g,%g) PLQ(%g,%g,%g,%g) test_REQ(%g,%g,%g,%g) K=%g alphaMorse=%g \n", kind, n, W.gridFF.grid.dCell.xx,W.gridFF.grid.dCell.yy,W.gridFF.grid.dCell.zz, PLQ.x,PLQ.y,PLQ.z,PLQ.w,  test_REQ.x,test_REQ.y,test_REQ.z,test_REQ.w, K, W.gridFF.alphaMorse );
    double R2Q=RQ*RQ;
    Quat4d bak_REQ;
    Quat4f bak_PLQ;
    Vec3d  bak_pos;
    int    bak_n;
    bool bModeNBmol = (kind==0)||(kind==1)||(kind==2)||(kind==3);
    if( bModeNBmol ){
        bak_n  =W.nbmol.natoms;  W.nbmol.natoms =1;
        bak_REQ=W.nbmol.REQs[0]; W.nbmol.REQs[0]=test_REQ;
        bak_PLQ=W.nbmol.PLQs[0]; W.nbmol.PLQs[0]=PLQ;
        bak_pos=W.nbmol.apos[0];
    }
    //W.gridFF.alphaMorse = 1.6;
    //printf( "!!!!!!!! MMFF_lib::sampleSurf_vecs() K=%g alphaMorse=%g \n", K, W.gridFF.alphaMorse  );
    if( fabs(K-W.gridFF.alphaMorse) > 1e-6 ){ printf("ERROR in sampleSurf_vecs K(%20.10f) does not match gridFF.alphaMorse(%20.10f) => exit()\n", K, W.gridFF.alphaMorse );  exit(0); }
    for(int i=0; i<n; i++){
        //printf( "sampleSurf_vecs()[%i]\n", i  );
        Quat4f fe  =Quat4fZero;
        Quat4d fe_d=Quat4dZero;
        Vec3d pos = poss[i];
        if(bModeNBmol){
            W.nbmol.apos[0]=pos;
            W.ffl.cleanForce();
        }
        //printf( "[%i] (%g,%g,%g)\n", i, W.nbmol.apos[0].x,W.nbmol.apos[0].y,W.nbmol.apos[0].z );
        switch(kind){
            case  0: fe_d.e= W.nbmol.evalR           (W.surf );                                                              FEs[i]=fe_d; break; 
            case  1: fe_d.e= W.nbmol.evalMorse       (W.surf, false, K,RQ  );                      fe_d.f=W.nbmol.fapos[0];  FEs[i]=fe_d; break; 
            case  2: fe_d.e= W.nbmol.evalMorsePBC    (W.surf, W.gridFF.grid.cell, nPBC, K, RQ  );  fe_d.f=W.nbmol.fapos[0];  FEs[i]=fe_d; break; 
            case  3: fe.e  = W.nbmol.evalMorsePLQ    (W.surf, W.gridFF.grid.cell, nPBC, K, R2Q );  fe_d.f=W.nbmol.fapos[0];  FEs[i]=fe_d; break; 
            // TODO: we should calculate interaction of test atom with  test_REQ
            // see gridFF.bindSystem(surf.natoms, surf.atypes, surf.apos, surf.REQs ) in MolWorld_sp3::initGridFF()
            //double evalMorsePLQ( NBFF& B, Mat3d& cell, Vec3i nPBC, double K=-1.0, double RQ=1.0 ){
            // nbmol .evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );
            
            // evalMorsePBC(  Vec3d pi, Quat4d REQi, Vec3d& fi, int natoms, Vec3d * apos, Quat4d * REQs ){
            case  9:         fe_d.e = W.gridFF.evalMorsePBC_sym( pos, test_REQ, fe_d.f );  FEs[i]=fe_d;  break;
            case 10:         W.gridFF.addForce_surf     (pos, {1.,0.,0.}, fe );   FEs[i]=(Quat4d)fe;  break;
            case 11:         W.gridFF.addForce_surf     (pos, PLQ, fe );          FEs[i]=(Quat4d)fe;  break;
            //case 12:        W.gridFF.addForce         (pos, PLQ, fe );         FEs[i]=(Quat4d)fe;  break;
            case 12: fe     = W.gridFF.getForce         (pos, PLQ        );      FEs[i]=(Quat4d)fe;  break;
            case 15: fe     = W.gridFF.getForce         (pos, {1.,0.,0.} );      FEs[i]=(Quat4d)fe;  break;
            case 16: fe     = W.gridFF.getForce         (pos, {0.,1.,0.} );      FEs[i]=(Quat4d)fe;  break;
            case 18: FEs[i] = W.gridFF.evalMorsePBC_PLQ_sym( pos, PLQ_d ); break;
            case 19: FEs[i] = W.gridFF.evalMorsePBC_PLQ_sym( pos, {1.,0.,0.,0.} ); break;
            case 20: FEs[i] = W.gridFF.evalMorsePBC_PLQ_sym( pos, {0.,1.,0.,0.} ); break;
            case  8: fe     = W.gridFF.getForce         (pos+Vec3d{W.gridFF.grid.dCell.xx*1.0,W.gridFF.grid.dCell.yy*1.0,W.gridFF.grid.dCell.zz*1.0}, PLQ     );         FEs[i]=(Quat4d)fe;  break;


            case 13: fe_d = W.gridFF.getForce_d       (pos, PLQ_d   );         FEs[i]=fe_d;        break;
            case 14: fe_d = W.gridFF.getForce_Tricubic(pos, PLQ_d   );         FEs[i]=fe_d;        break;
        }
        //fs[i]=(Vec3d)(fe.f);
        //Es[i]=fe.e;
        //FEs[i] = (Quat4d)fe; 
    }
    if( bModeNBmol ){
        W.nbmol.natoms=bak_n;    
        W.nbmol.REQs[0]=bak_REQ; 
        W.nbmol.PLQs[0]=bak_PLQ;
        W.nbmol.apos[0]=bak_pos; 
    }
}

void change_lvec( double* lvec, bool bAdd, bool  ){
    if(bAdd){ W.add_to_lvec( *(Mat3d*)lvec ); }
    else    { W.change_lvec( *(Mat3d*)lvec ); }
}

void optimizeLattice_1d( double* dlvec, int n1, int n2, int initMode, double tol ){
    printf("MMFF_lib::optimizeLattice_1d(n1=%i,n2=%i,initMode=%i,tol=%g) \n", n1, n2, initMode, tol );
    W.gopt.tolerance=tol;
    W.gopt.initMode =initMode; 
    W.optimizeLattice_1d( n1, n2, *(Mat3d*)dlvec );
}



void addSnapshot(bool ifNew = false, char* fname = 0){
    W.addSnapshot(ifNew, fname);
}

void printDatabase(){
    W.printDatabase();
}

void computeDistance(int i, int j, double* dist){
    *dist = W.computeDistance(i,j);
}

} // extern "C"
