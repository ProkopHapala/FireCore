
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

void init_buffers_UFF(){
    printf( "init_buffers_UFF() \n" );
    // Common buffers
    buffers["apos"]  = (double*)W.nbmol.apos;
    buffers["fapos"] = (double*)W.nbmol.fapos;
    buffers["REQs"]  = (double*)W.nbmol.REQs;
    //buffers.insert( { "PLQs",   (double*)W.nbmol.PLQs  } );
    // UFF-specific buffers
    if(W.bUFF){ // UFF-specific buffers
        buffers["hneigh"]    = (double*)W.ffu.hneigh;
        buffers["fint"]      = (double*)W.ffu.fint;
        buffers["bonParams"] = (double*)W.ffu.bonParams;
        buffers["angParams"] = (double*)W.ffu.angParams;
        buffers["dihParams"] = (double*)W.ffu.dihParams;
        buffers["invParams"] = (double*)W.ffu.invParams;

        ibuffers["neighs"]   = (int*)W.ffu.neighs;
        ibuffers["neighBs"]  = (int*)W.ffu.neighBs;
        ibuffers["bonAtoms"] = (int*)W.ffu.bonAtoms;
        ibuffers["angAtoms"] = (int*)W.ffu.angAtoms;
        ibuffers["dihAtoms"] = (int*)W.ffu.dihAtoms;
        ibuffers["invAtoms"] = (int*)W.ffu.invAtoms;

        // neighbor indices for angles, dihedrals, inversions
        ibuffers["angNgs"] = (int*)W.ffu.angNgs;
        ibuffers["dihNgs"] = (int*)W.ffu.dihNgs;
        ibuffers["invNgs"] = (int*)W.ffu.invNgs;
    }
    // UFF-specific dimensions
    if(W.bUFF){
        ibuffers["ndims"] = &W.ffu._natoms;
        buffers ["Es"]    = &W.ffu.Etot;
    }
    ibuffers["selection"] = W.manipulation_sel;
    bbuffers["ffflags"]   = &W.doBonded;
    // int _natoms, nbonds, nangles, ndihedrals, ninversions, nf; // 5
    // int i0dih,i0inv,i0ang,i0bon;                               // 4
    // double Etot, Eb, Ea, Ed, Ei;                               // 5
    printf( "MMFF_lib.cpp::init_buffers_UFF() ndims{natoms=%i, nbonds=%i, nangles=%i, ndihedrals=%i, ninversions=%i, nf=%i, i0dih=%i,i0inv=%i,i0ang=%i,i0bon=%i, }\n", W.ffu._natoms, W.ffu.nbonds, W.ffu.nangles, W.ffu.ndihedrals, W.ffu.ninversions, W.ffu.nf, W.ffu.i0dih, W.ffu.i0inv, W.ffu.i0ang, W.ffu.i0bon );
    printf( "MMFF_lib.cpp::init_buffers_UFF() Es{ Etot=%f, Eb=%f, Ea=%f, Ed=%f, Ei=%f, }\n", W.ffu.Etot, W.ffu.Eb, W.ffu.Ea, W.ffu.Ed, W.ffu.Ei );
}

void init_buffers(){
    printf( "init_buffers() \n" );
    buffers["apos"]  = (double*)W.nbmol.apos;
    buffers["fapos"] = (double*)W.nbmol.fapos;
    buffers["REQs"]  = (double*)W.nbmol.REQs;
    if(W.bMMFF){
        buffers["DOFs"]  = W.ffl.DOFs;
        buffers["fDOFs"] = W.ffl.fDOFs;
        buffers["vDOFs"] = W.opt.vel;
        //buffers .insert( { "REQs",   (double*)W.ffl.REQs  } );
        buffers["PLQs"] = (double*)W.ffl.PLQd;
        if(!W.bUFF){
            buffers["pipos"]  = (double*)W.ffl.pipos;
            buffers["fpipos"] = (double*)W.ffl.fpipos;
            ibuffers["neighs"] = (int*)W.ffl.neighs;
        } // else{ // UFF-specific}
    }else{
        W.ff.natoms=W.nbmol.natoms;
    }
    printf( "MMFF_lib.cpp::init_buffers() ndims{nDOFs=%i,natoms=%i,nnode=%i,ncap=%i,npi=%i,nbonds=%i,nvecs=%i,ne=%i,ie0=%i}\n", W.ff.nDOFs, W.ff.natoms, W.ff.nnode, W.ff.ncap, W.ff.npi, W.ff.nbonds, W.ff.nvecs, W.ff.ne, W.ff.ie0 );
    ibuffers["ndims"]    = &W.ff.nDOFs;
    buffers ["Es"]       = &W.ff.Etot;
    ibuffers["selection"] = W.manipulation_sel;
    bbuffers["ffflags"]   = &W.doBonded;
    //printBuffNames();
}


void print_debugs( bool bParams, bool bNeighs, bool bShifts, bool bAtoms ){
    printf("print_debugs() W.bUFF=%i, bParams=%i, bNeighs=%i, bShifts=%i, bAtoms=%i \n", W.bUFF, bParams, bNeighs, bShifts, bAtoms);
    if(W.bUFF){
        W.ffu.printSizes();
        if(bParams) W.ffu.printAllParams(true, true, true, true, true);
        if(bAtoms ) W.ffu.print();
    } else {
        W.ffl.printSizes();
        if( bParams ) W.ffl.printAtomParams();
        if( bNeighs ) W.ffl.printNeighs();
        if( bShifts ) W.ffl.print_pbc_shifts();
        if( bAtoms  ) W.ffl.print();
    }
}

void print_setup(){
    if(W.bUFF){
        W.ffu.printSimulationSetup();
    }else{
        printf("MMFF_lib::print_setup() bUFF=false\n");
       //W.ffl.printSimulationSetup();
    }
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

    // unbuffered printf()
    setbuf(stdout, NULL);
    setbuf(stderr, NULL);

    W.initParams( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
    bool bGrid = gridStep>0;
    // initialize the main
    //W.init( bGrid, bUFF );
    W.bGridFF=bGrid;
    W.bUFF   =bUFF;
    int ret = W.init();
    if(ret!=0){
        printf("MMFF_lib::init() failed ret=%i -> return nullptr\n", ret);
        return nullptr;
    }
    //init_buffers();
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

void setSurfFlatPlane( double* pos0, double* normal ){
    W.setSurfFlatPlane( *(Vec3d*)pos0, *(Vec3d*)normal );
}

void setSurfFlatParams( int mode, double* REQ, double K ){
    W.setSurfFlatParams( mode, *(Quat4d*)REQ, K );
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

    if(W.bUFF){
        if(omp){ return W.ffu.run_omp(nstepMax,dt,Fconv,   10.0, -1.0,     outE, outF, outV, outVF ); }
        else   { return W.ffu.run    (nstepMax,dt,Fconv, 1000.0,  damping, outE, outF, outV, outVF ); }
    }else{
        if(omp){ return W.run_omp   (nstepMax,dt,Fconv,   10.0, -1.0,     outE, outF, outV, outVF ); }
        else   { return W.run_no_omp(nstepMax,dt,Fconv, 1000.0,  damping, outE, outF, outV, outVF ); }
    }
    //else   { return W.run       (nstepMax,dt,Fconv,ialg,       outE, outF, outV, outVF ); }
}

void  scan( int nconf, double* poss, double* rots, double* Es, double* aforces, double* aposs, bool omp, bool bRelax, int niter_max, double dt, double Fconv, double Flim ){
    if(bRelax){
        if(omp){ printf("ERROR: scan_relaxed() not implemented witht OMP\n"); exit(0); } 
        else   { W.scan_relaxed( nconf, (Vec3d*)poss, (Mat3d*)rots, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp, niter_max, dt, Fconv, Flim );  }
    }else{
        if(omp){ printf("ERROR: scan_rigid() not implemented witht OMP\n"); exit(0); } 
        else   { W.scan_rigid( nconf, (Vec3d*)poss, (Mat3d*)rots, Es, (Vec3d*)aforces, (Vec3d*)aposs, omp ); }
    }
}

void scan_manipulation( int nconf, double* ts, const char* spline_fname, int iAnchor, double Kanchor, const char* trjName, int* nPBC, double* Es, double* aforces, double* aposs, double* fconstr, int niter_max, double dt, double Fconv, double Flim ){
    printf("MMFF_lib::scan_manipulation(nconf=%i) spline=%s iAnchor=%i Kanchor=%g trjName=%s nPBC(%i,%i,%i)\n",
        nconf, spline_fname?spline_fname:"NULL", iAnchor, Kanchor, trjName?trjName:"NULL", nPBC?nPBC[0]:-1, nPBC?nPBC[1]:-1, nPBC?nPBC[2]:-1
    );
    if(spline_fname==0){ printf("ERROR in MMFF_lib::scan_manipulation() spline_fname==NULL => exit()\n"); exit(0); }
    W.tipSpline.loadDat( spline_fname );
    W.tipSpline.setAnchor( iAnchor );
    W.tipSpline.Kanchor = Kanchor;
    W.tipSpline.bActive = true;
    Vec3i npbc_ = nPBC ? *(Vec3i*)nPBC : Vec3i{1,1,1};
    W.scan_manipulation( nconf, ts, Es, (Vec3d*)aforces, (Vec3d*)aposs, (Vec3d*)fconstr, trjName, npbc_, niter_max, dt, Fconv, Flim );
}

// In MMFF_lib.cpp before extern "C" closing
void scan_atoms_rigid(int nscan, int nsel, int* inds, double* scan_pos, double* out_forces, double* out_Es, bool bRelative ){
    int na = W.nbmol.natoms;
    std::vector<Vec3d> orig(na);
    for(int i=0; i<na; i++){ orig[i] = W.nbmol.apos[i]; }
    for(int i=0; i<nscan; i++){
        Vec3d* ps = (Vec3d*)(scan_pos+i*nsel*3);
        for(int j=0; j<nsel; j++){ 
            int ia = inds[j];  
            Vec3d p = ps[j]; 
            if(bRelative){ p.add(orig[ia]); }
            W.nbmol.apos[ia] = p; 
        }
        if(out_Es){ out_Es[i] = W.eval_no_omp(); }
        if(out_forces){ for(int j=0; j<nsel; j++){ int ia = inds[j];  out_forces[i*nsel*3 + j*3] = W.nbmol.fapos[ia].x; } }
    }
    for(int j=0; j<nsel; j++){ int ia = inds[j];  W.nbmol.apos[ia] = orig[ia]; }
}

// Hessian: independent 3×3 blocks for selected atoms
void getHessian3x3( int n, int* inds, double* Hess_, double dx, bool bDiag ){
    int na=W.nbmol.natoms;
    printf("getHessian3x3(n=%i) na=%i dx=%g    bMMFF=%i bNonBonded=%i bSurfAtoms=%i bGridFF=%i bPBC=%i bNonBondNeighs=%i \n", n, na, dx, W.bMMFF, W.bNonBonded, W.bSurfAtoms, W.bGridFF, W.bPBC, W.bNonBondNeighs);
    // save original positions
    DEBUG
    std::vector<Vec3d> orig(n);
    for(int i=0;i<n;i++){ int ia=inds[i]; orig[i]=W.nbmol.apos[ia]; }
    double denom = 1.0/(2*dx);
    DEBUG
    W.saveXYZ("getHessian3x3.xyz");
    DEBUG
    Mat3d H, U;
    Vec3d Ks;
    for(int i=0;i<n;i++){
        printf( "getHessian3x3() i=%i \n", i );
        int ia=inds[i];
        printf( "getHessian3x3() i=%i ia=%i  \n", i, ia );
        Vec3d  p0=orig[i];
        Vec3d& p=W.nbmol.apos[ia];
        //double* H = out_hessians+i*9;
        //Vec3d* Hess = ((Vec3d*)out_hessians)+(i*3); 
        printf( "getHessian3x3() ia=%i p=(%g,%g,%g)\n", ia, p.x,p.y,p.z );
        for(int k=0;k<3;k++){
            p.array[k]=p0.array[k]+dx; W.eval_no_omp(); Vec3d df=W.nbmol.fapos[ia];
            p.array[k]=p0.array[k]-dx; W.eval_no_omp(); df.sub(  W.nbmol.fapos[ia] );
            p.array[k]=p0.array[k];                     df.mul(denom);
            printf( "getHessian3x3() ia=%i k=%i p=(%g,%g,%g) df=(%g,%g,%g)\n", ia, k, p.x,p.y,p.z, df.x,df.y,df.z );
            //Hess[k]=df;
            H.vecs[k]=df;
            //for(int l=0;l<3;l++) H[l*3+k]=(fp.array[l]-fm.array[l])*denom;
        }
        if(bDiag){
            H.eigenvals(Ks);
            printf( "getHessian3x3() ia=%i Ks=(%g,%g,%g)\n", ia, Ks.x,Ks.y,Ks.z );
            H.eigenvec(Ks.x,U.a);
            H.eigenvec(Ks.y,U.b);
            H.eigenvec(Ks.z,U.c);
            printf( "getHessian3x3() ia=%i u1(%g,%g,%g) u2(%g,%g,%g) u3(%g,%g,%g)\n", ia, U.a.x,U.a.y,U.a.z,   U.b.x,U.b.y,U.b.z,   U.c.x,U.c.y,U.c.z );
            *((Mat3d*)(Hess_ + i*12   ))=U;
            *((Vec3d*)(Hess_ + i*12 +9))=Ks;
        }else{
            *((Mat3d*)(Hess_ + i*12))=H;
        }
    }
    DEBUG
    // restore original positions
    for(int i=0;i<n;i++){ int ia=inds[i]; W.nbmol.apos[ia]=orig[i]; }
}

void getHessian3Nx3N(int n,int* inds,double* out_hessian,double dx){
    // Phi(ia,ja,R) = -dF_ja(0)/du_ia via central FD.  Bloch D(k) applied downstream in Python.
    // Periodic crystal: init with lvs xyz + nPBC>0 so bonded/H-bond neighbors include cell images.
    printf("getHessian3Nx3N(n=%i) dx=%g bPBC=%i nPBC=(%i,%i,%i) bMMFF=%i\n",
        n, dx, (int)W.bPBC, W.nPBC.x, W.nPBC.y, W.nPBC.z, (int)W.bMMFF);
    if(!W.bPBC){
        printf("WARNING: getHessian3Nx3N bPBC=0 — cluster Hessian; for bulk phonons use lvs xyz + nPBC>0\n");
    }
    std::vector<Vec3d> orig(n);
    int dim = n * 3;
    for(int i=0;i<dim*dim;i++){ out_hessian[i]=0.0; }
    for(int i=0;i<n;i++){ int ia=inds[i]; orig[i]=W.nbmol.apos[ia]; }
    for(int p=0;p<n;p++){
        int ip=inds[p];
        for(int k=0;k<3;k++){
            double v=orig[p].array[k];
            W.nbmol.apos[ip].array[k]=v+dx; W.eval_no_omp();
            for(int o=0; o<n; o++){
                int io = inds[o];
                for(int l=0; l<3; l++){
                    out_hessian[(o*3+l)*dim + (p*3+k)] = -W.nbmol.fapos[io].array[l];
                }
            }
            W.nbmol.apos[ip].array[k]=v-dx; W.eval_no_omp();
            for(int o=0; o<n; o++){
                int io = inds[o];
                for(int l=0; l<3; l++){
                    out_hessian[(o*3+l)*dim + (p*3+k)] = (out_hessian[(o*3+l)*dim + (p*3+k)] + W.nbmol.fapos[io].array[l]) / (2*dx);
                }
            }
            W.nbmol.apos[ip].array[k]=v;
        }
    }
    for(int i=0;i<n;i++){ int ia=inds[i]; W.nbmol.apos[ia]=orig[i]; }
    for(int i=0;i<dim;i++){
        for(int j=i+1;j<dim;j++){
            double s = 0.5*(out_hessian[i*dim+j]+out_hessian[j*dim+i]);
            out_hessian[i*dim+j]=s;
            out_hessian[j*dim+i]=s;
        }
    }
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

    // Keep subtraction/clamp consistent with nonbond evaluation.
    // NOTE: switch semantics are 0=keep, >0=force true, <0=force false.
    if(NonBonded!=0){
        bool b = (NonBonded>0);
        W.ffl.bSubtractBondNonBond  = b;
        W.ffl.bSubtractAngleNonBond = b;
        W.ffl.bClampNonBonded       = b;
    }

    printf( "setSwitches2() W.bCheckInvariants==%i bPBC=%i | bNonBonded=%i bNonBondNeighs=%i | bSurfAtoms=%i bGridFF=%i | bMMFF=%i doAngles=%i doPiSigma=%i doPiPiI=%i | subNB(bond=%i,angle=%i,clamp=%i) \n",
        W.bCheckInvariants, W.bPBC,  W.bNonBonded, W.bNonBondNeighs, W.bSurfAtoms, W.bGridFF, W.bMMFF, W.ffl.doAngles, W.ffl.doPiSigma, W.ffl.doPiPiI,
        W.ffl.bSubtractBondNonBond, W.ffl.bSubtractAngleNonBond, W.ffl.bClampNonBonded
    );

    #undef _setbool
}

void setSwitchesUFF( int DoBond, int DoAngle, int DoDihedral, int DoInversion, int DoAssemble, int SubtractBondNonBond, int ClampNonBonded ){
    // bool bDoBond=true, bDoAngle=true, bDoDihedral=true, bDoInversion=true; bSubtractBondNonBond, bClampNonBonded
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.ffu.bDoBond,              DoBond );
    _setbool( W.ffu.bDoAngle,             DoAngle );
    _setbool( W.ffu.bDoDihedral,          DoDihedral );
    _setbool( W.ffu.bDoInversion,         DoInversion );
    _setbool( W.ffu.bDoAssemble,          DoAssemble );
    _setbool( W.ffu.bSubtractBondNonBond, SubtractBondNonBond );
    _setbool( W.ffu.bClampNonBonded,      ClampNonBonded );
    printf( "setSwitchesUFF() bDoBond=%i bDoAngle=%i bDoDihedral=%i bDoInversion=%i bDoAssemble=%i bSubtractBondNonBond=%i bClampNonBonded=%i \n", W.ffu.bDoBond, W.ffu.bDoAngle, W.ffu.bDoDihedral, W.ffu.bDoInversion, W.ffu.bDoAssemble, W.ffu.bSubtractBondNonBond, W.ffu.bClampNonBonded );
    #undef _setbool
}

void setSwitchesUFF_NB( int NonBonded, int NonBondNeighs, int SubtractAngleNonBond ){
    #define _setbool(b,i) { if(i>0){b=true;}else if(i<0){b=false;} }
    _setbool( W.ffu.bNonBonded           , NonBonded           );
    _setbool( W.ffu.bNonBondNeighs       , NonBondNeighs       );
    _setbool( W.ffu.bSubtractAngleNonBond, SubtractAngleNonBond);
    printf( "setSwitchesUFF_NB() bNonBonded=%i bNonBondNeighs=%i bSubtractAngleNonBond=%i \n", W.ffu.bNonBonded, W.ffu.bNonBondNeighs, W.ffu.bSubtractAngleNonBond );
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
    W.gridFF.grid.printCell();
    printf( "sampleSurf_new() gff.shift0(%g,%g,%g) gff.pos0(%g,%g,%g)\n", W.gridFF.shift0.x, W.gridFF.shift0.y, W.gridFF.shift0.z, W.gridFF.grid.pos0.x, W.gridFF.grid.pos0.y, W.gridFF.grid.pos0.z );
    printf( "sampleSurf_new() gridN(%i,%i,%i) \n", W.gridFF.grid.n.x, W.gridFF.grid.n.y, W.gridFF.grid.n.z );
    
    //PLQd=Quat4d{1.0,0.0,0.0,0.0};
    {
        Vec3i ng = W.gridFF.grid.n;
        Vec3d g0 = W.gridFF.grid.pos0;
        Vec3d dg = Vec3d{ W.gridFF.grid.dCell.xx, W.gridFF.grid.dCell.yy, W.gridFF.grid.dCell.zz };
        Quat4d C = PLQd;
        Quat4i* xqs = W.gridFF.cubic_xqis;
        printf("CPU sampleSurf_new() ng(%i,%i,%i) g0(%g,%g,%g) dg(%g,%g,%g) C(%g,%g,%g,%g) \n", ng.x,ng.y,ng.z,   g0.x,g0.y,g0.z,   dg.x,dg.y,dg.z,   C.x,C.y,C.z,C.w );
        printf("CPU sampleSurf_new() xqs[0](%i,%i,%i,%i) xqs[1](%i,%i,%i,%i) xqs[2](%i,%i,%i,%i) xqs[3](%i,%i,%i,%i)\n", xqs[0].x, xqs[0].y, xqs[0].z, xqs[0].w,   xqs[1].x, xqs[1].y, xqs[1].z, xqs[1].w,   xqs[2].x, xqs[2].y, xqs[2].z, xqs[2].w,  xqs[3].x, xqs[3].y, xqs[3].z, xqs[3].w   );
    }

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
    const bool bKindFlat = (kind==30)||(kind==31);
    if( (!bKindFlat) && ( fabs(K-W.gridFF.alphaMorse) > 1e-6 ) ){ printf("ERROR in sampleSurf_vecs K(%20.10f) does not match gridFF.alphaMorse(%20.10f) => exit()\n", K, W.gridFF.alphaMorse );  exit(0); }
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

            case 30: W.setSurfFlatParams( MolWorld_sp3::SurfFlat_HamakerLJ93, W.surfFlat_REQ, W.surfFlat_K ); fe_d.e = W.evalSurfFlat( pos, test_REQ, fe_d.f ); FEs[i]=fe_d; break;
            case 31: W.setSurfFlatParams( MolWorld_sp3::SurfFlat_Morse      , W.surfFlat_REQ, K            ); fe_d.e = W.evalSurfFlat( pos, test_REQ, fe_d.f ); FEs[i]=fe_d; break;
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

// =========================== Hessian Basis Matrix Functions (for parameter fitting)
// Consolidated API: just 2 functions - getHessianContext() and setParams()

/// @brief Get all information needed for Hessian parameter fitting in one call
/// 
/// Returns: basis matrices, atom indices, and current parameters for bonds and angles.
/// The linear model is: H = sum_i k_i * B_i where B_i are geometric basis matrices.
///
/// @param[out] n_bonds       Number of bonds (0 if not UFF)
/// @param[out] n_angles      Number of angles (0 if not UFF)
/// @param[out] n_atoms       Number of atoms
/// @param[out] bond_bases    Bond basis matrices, size nbonds*36 (4 blocks: Bii,Bij,Bji,Bjj each 3x3)
/// @param[out] bond_atoms    Bond atom indices [i,j] per bond, size nbonds*2
/// @param[out] bond_params   Current bond stiffness k values, size nbonds
/// @param[out] angle_bases   Angle basis matrices, size nangles*81 (9x9 block per angle)
/// @param[out] angle_atoms   Angle atom indices [i,j,k] per angle, size nangles*3  
/// @param[out] angle_params  Current angle stiffness k values, size nangles
///
/// Arrays can be NULL if that data is not needed. All outputs are only written if W.bUFF.
void getHessianContext(int* n_bonds, int* n_angles, int* n_atoms,
                       double* bond_bases, int* bond_atoms, double* bond_params,
                       double* angle_bases, int* angle_atoms, double* angle_params){
    if(!W.bUFF){
        if(n_bonds) *n_bonds = 0;
        if(n_angles) *n_angles = 0;
        if(n_atoms) *n_atoms = W.nbmol.natoms;
        return;
    }
    
    int nb = W.ffu.nbonds;
    int na = W.ffu.nangles;
    int nat = W.ffu.natoms;
    
    if(n_bonds) *n_bonds = nb;
    if(n_angles) *n_angles = na;
    if(n_atoms) *n_atoms = nat;
    
    // Bond basis matrices (geometric, without k factor)
    if(bond_bases){
        W.ffu.getAllBondHessianBases(bond_bases);
    }
    
    // Bond atom indices [i,j]
    if(bond_atoms){
        for(int ib=0; ib<nb; ib++){
            bond_atoms[ib*2+0] = W.ffu.bonAtoms[ib].x;
            bond_atoms[ib*2+1] = W.ffu.bonAtoms[ib].y;
        }
    }
    
    // Current bond parameters (k values)
    if(bond_params){
        for(int ib=0; ib<nb; ib++){
            bond_params[ib] = W.ffu.bonParams[ib].x;
        }
    }
    
    // Angle basis matrices (geometric, without k factor)
    if(angle_bases){
        W.ffu.getAllAngleHessianBases(angle_bases);
    }
    
    // Angle atom indices [i,j,k]
    if(angle_atoms){
        for(int ia=0; ia<na; ia++){
            angle_atoms[ia*3+0] = W.ffu.angAtoms[ia].x;
            angle_atoms[ia*3+1] = W.ffu.angAtoms[ia].y;
            angle_atoms[ia*3+2] = W.ffu.angAtoms[ia].z;
        }
    }
    
    // Current angle parameters (k values)
    if(angle_params){
        for(int ia=0; ia<na; ia++){
            angle_params[ia] = W.ffu.angParams[ia].k;
        }
    }
}

/// @brief Set force field stiffness parameters after fitting
/// 
/// @param n_bonds     Number of bonds (must match system, or 0 to skip)
/// @param bond_params New bond stiffness k values [k_0, k_1, ...], size n_bonds, or NULL
/// @param n_angles    Number of angles (must match system, or 0 to skip)  
/// @param angle_params New angle stiffness k values [k_0, k_1, ...], size n_angles, or NULL
///
/// Only updates parameters if W.bUFF and counts match.
void setParams(int n_bonds, const double* bond_params, int n_angles, const double* angle_params){
    if(!W.bUFF) return;
    
    // Update bond parameters
    if(bond_params && n_bonds > 0){
        if(n_bonds != W.ffu.nbonds){
            printf("WARNING setParams: n_bonds=%d != system nbonds=%d\n", n_bonds, W.ffu.nbonds);
            return;
        }
        for(int ib=0; ib<n_bonds; ib++){
            W.ffu.bonParams[ib].x = bond_params[ib];
        }
    }
    
    // Update angle parameters
    if(angle_params && n_angles > 0){
        if(n_angles != W.ffu.nangles){
            printf("WARNING setParams: n_angles=%d != system nangles=%d\n", n_angles, W.ffu.nangles);
            return;
        }
        for(int ia=0; ia<n_angles; ia++){
            W.ffu.angParams[ia].k = angle_params[ia];
        }
    }
}

/// @brief Assemble full Hessian from parameter values (for validation)
/// @param n_params total number of parameters (nbonds + nangles)
/// @param params array of parameter values [k_bonds..., k_angles...]
/// @param out_hessian output 3N x 3N Hessian matrix (row-major, size 3*natoms)
void assembleHessianFromParams(int n_params, const double* params, double* out_hessian){
    if(!W.bUFF) return;
    int dim = W.ffu.natoms * 3;
    // Initialize to zero
    for(int i=0; i<dim*dim; i++) out_hessian[i] = 0.0;
    
    // Add bond contributions
    for(int ib=0; ib<W.ffu.nbonds; ib++){
        W.ffu.assembleBondHessianBasis(ib, params[ib], out_hessian);
    }
    
    // Add angle contributions (similar pattern)
    int offset = W.ffu.nbonds;
    for(int ia=0; ia<W.ffu.nangles; ia++){
        double B[81], g[9];
        W.ffu.getAngleHessianBasis(ia, g, B);
        double k = params[offset + ia];
        const Vec3i& ijk = W.ffu.angAtoms[ia];
        int atoms[3] = {ijk.x, ijk.y, ijk.z};
        
        // Scale and add 9x9 block to global Hessian
        for(int r=0; r<9; r++){
            for(int c=0; c<9; c++){
                double val = B[r*9+c] * k;
                int atom_r = atoms[r/3];
                int atom_c = atoms[c/3];
                int idx_r = atom_r*3 + (r%3);
                int idx_c = atom_c*3 + (c%3);
                out_hessian[idx_r*dim + idx_c] += val;
            }
        }
    }
}

} // extern "C"
