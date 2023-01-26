﻿
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
#include "OCL_PP.h"


void     pack    (int n, Vec3d* fs, Quat4f* qs, float K ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];  qs[i].e=K; } }
void   unpack    (int n, Vec3d* fs, Quat4f* qs          ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f;           } }
double unpack_add(int n, Vec3d* fs, Quat4f* qs          ){ double E=0; for(int i=0; i<n; i++){ fs[i].add( (Vec3d)qs[i].f ); E+=qs[i].e; }; return E; }

void pack(int n, Vec3d* fs, double* es, Quat4f* qs ){
    for(int i=0; i<n; i++){ 
        float e; if(es){e=es[i];}else{e=0;}
        //float e; if(es){e=es[i];}else{e=0;}
        qs[i].f=(Vec3f)fs[i];
        qs[i].e=e; 
    } 
}
Quat4f* pack(int n, Vec3d* fs, double* es=0){
    Quat4f* qs = new Quat4f[n];
    pack( n, fs, es, qs );
    return qs;
}

void unpack(int n, Vec3d* fs, double* es, Quat4f* qs ){
    for(int i=0; i<n; i++){ 
        const Quat4f& q= qs[i];
        fs[i]      =(Vec3d)q.f;
        if(es)es[i]=       q.e; 
    } 
}

//tempate<typename T> bool addFirstEmpty( T* arr, n, T what, T empty=-1 ){
bool addFirstEmpty( int* arr, int n, int what, int empty=-1 ){
    for(int i=0; i<n; i++){
        if(arr[i]==empty){ arr[i]=what; return true; }
    }
    return false;
};

// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_ocl : public MolWorld_sp3 { public:
    OCL_PP     ocl;

    Quat4f * q_ps=0;
    Quat4f * q_fs=0; 
    Quat4f * q_vs=0; 
    int* bkneighs=0;

    bool bGPU_MMFF = true;


// ======== Functions

void surf2ocl(Vec3i nPBC, bool bSaveDebug=false){
    int ncell = (nPBC.x*2+1) * (nPBC.y*2+1) * (nPBC.z*2+1); 
    int n = surf.n*ncell;
    printf( "surf2ocl() na(%i) = ncell(%i) * natom(%i)\n", n, ncell, surf.n );
    long T0=getCPUticks();
    Quat4f* atoms = new Quat4f[n];
    Quat4f* coefs = new Quat4f[n];
    double R2damp = sq(gridFF.Rdamp);
    int ii=0;
    for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
        Vec3d  p0 = gridFF.grid.cell.a*ia + gridFF.grid.cell.b*ib + gridFF.grid.cell.c*ic;
        for(int i=0; i<surf.n; i++){
            //printf("ii %i i %i iabc(%i,%i,%i)\n", ii, i, ia,ib,ic);
            atoms[ii].f=(Vec3f)(surf.ps[i]+p0);
            atoms[ii].e=R2damp; 
            coefs[ii].f=(Vec3f)surf.REQs[i];
            coefs[ii].e=gridFF.alpha;
            ii++;
        }
    }}}
    printf( ">>time(surf_to_GPU) %g \n", (getCPUticks()-T0)*tick2second );
    long T1=getCPUticks();
    ocl.makeGridFF( n, (float4*)atoms, (float4*)coefs );
    printf( ">>time(ocl.makeGridFF() %g \n", (getCPUticks()-T1)*tick2second );
    ocl.download( ocl.itex_FE_Paul, gridFF.FFPauli );  
    ocl.download( ocl.itex_FE_Lond, gridFF.FFLondon );
    ocl.download( ocl.itex_FE_Coul, gridFF.FFelec );
    printf( ">>time(ocl.makeGridFF(); ocl.download() %g \n", (getCPUticks()-T1)*tick2second );
    delete [] atoms;
    delete [] coefs;
    if(bSaveDebug){
        gridFF.grid.saveXSF( "ocl_E_Paul.xsf", (float*)gridFF.FFPauli,  4, 3 );
        gridFF.grid.saveXSF( "ocl_E_Lond.xsf", (float*)gridFF.FFLondon, 4, 3 );
        gridFF.grid.saveXSF( "ocl_E_Coul.xsf", (float*)gridFF.FFelec,   4, 3 );
        // ---- Save combined forcefield
        Quat4f * FFtot = new Quat4f[gridFF.grid.getNtot() ];
        Vec3d testREQ = (Vec3d){ 1.487, 0.0006808, 0.0}; // H
        gridFF.evalCombindGridFF ( testREQ, FFtot );
        gridFF.grid.saveXSF( "ocl_E_H.xsf",  (float*)FFtot, 4, 3, gridFF.natoms, gridFF.atypes, gridFF.apos );
        delete [] FFtot;
    }
    ocl.buffers[ocl.ibuff_atoms].release();
    ocl.buffers[ocl.ibuff_coefs].release();
}

void init_ocl(){
    // WARRNING - DEBUG
    for(int i=0; i<nbmol.n; i++){ nbmol.ps->addRandomCube(0.1); }; // DEBUG - displace atoms to test relaxation forces

    printf( "init_ocl() \n" );
    long T0 = getCPUticks();
    ocl.init();
    ocl.initPP( "common_resources/cl" );
    // ---- Init Surface force-field grid
    printf( ">>time(ocl.init();ocl.initPP(): %g [s] \n", (getCPUticks()-T0)*tick2second  );
    long T1 = getCPUticks();
    //surf2ocl( nPBC, true );
    surf2ocl( nPBC, false );
    printf( ">>time(surf2ocl): %g [s] \n", (getCPUticks()-T1)*tick2second  );
    printf( "... init_ocl() END\n" );
}

void REQs2ocl(){
    Quat4f* q_cs= new Quat4f[nbmol.n];
    pack  ( nbmol.n, nbmol.REQs, q_cs, gridFF.alpha );
    ocl.upload( ocl.ibuff_coefs, q_cs );    
    delete [] q_cs;
};

int ithNeigh(int ia, int ja){
    int* aneighs = (int*)nbmol.neighs;
    int* ngs     = aneighs + ia*4;
    for(int i=0; i<4; i++){
        if(ngs[i]==ja) return i; 
    }
    return -1;
}

void MMFFparams2ocl(){
    //int n=nbmol.n;
    int n=ff.natoms;
    float*  blks= new float[n*8];
    Quat4f* a0ks= new Quat4f[ff.nnode];
    for(int i=0; i<n*8; i++){ blks[i]=1000.0; }
    for(int i=0; i<ff.nbonds; i++){
        Vec2i b = ff.bond2atom[i];
        double K  = ff.bond_k [i];
        double l0 = ff.bond_l0[i];
        int ii = ithNeigh( b.i, b.j );
        int jj = ithNeigh( b.j, b.i );
        printf( "bond[%i] (%i,%i)-(%i,%i) l0 %g K %g \n", i, b.i,ii,  b.j,jj,   l0, K );
        if(ii>=0){ blks[b.i*8+ii] = l0; blks[b.i*8+ii+4] = K; }else{ printf("ii<0 \n"); exit(0); };
        if(jj>=0){ blks[b.j*8+jj] = l0; blks[b.j*8+jj+4] = K; }else{ printf("jj<0 \n"); exit(0); };  
    }
    for(int i=0; i<ff.nnode; i++){
        a0ks[i]=(Quat4f)ff.NeighParams[i];
    }
    ocl.upload( ocl.ibuff_bondLK, blks );  
    ocl.upload( ocl.ibuff_ang0K,  a0ks );  
    delete [] blks;
    delete [] a0ks;
}

void mol2ocl(){
    int n   = nbmol.n;
    int npi = 0;
    if(bGPU_MMFF){ npi=ff.npi; }
    int n2 = n+npi; 
    // ---- Prepare for sampling atoms
    ocl.initAtomsForces( n, npi, ff.nnode, bGPU_MMFF );
    printf( "mol2ocl() n,pi(%i,%i) buffs atoms, coefs, aforces, neighs: %i %i %i %i \n", n,npi, ocl.ibuff_atoms, ocl.ibuff_coefs, ocl.ibuff_aforces, ocl.ibuff_neighs );
    q_ps= new Quat4f[n2];
    q_fs= new Quat4f[n2];
    // ---- nbmol coefs
    pack      ( n2, nbmol.ps, q_ps, sq(gridFF.Rdamp) );
    ocl.upload( ocl.ibuff_atoms, q_ps );  
    REQs2ocl();
    makeOCLNeighs( );
    if(bGPU_MMFF){
        Quat4f* q_vs= new Quat4f[n2];
        for(int i=0; i<n2; i++){ q_vs[i].set(0.); }
        ocl.upload( ocl.ibuff_avel,   q_vs );
        makeBackNeighs( nbmol.neighs );
        ocl.upload( ocl.ibuff_bkNeighs, bkneighs );
        MMFFparams2ocl();
        makePi0s();
    }
}

void makePi0s(){
    if(ff.pi0s==0){ ff.pi0s=new Quat4f[ff.npi]; };
    ff.evalPi0s();
    ocl.upload( ocl.ibuff_pi0s, ff.pi0s ); 
}

void  makeOCLNeighs( ){
    printf("makeOCLNeighs() n %i nnode %i \n", nbmol.n, ff.nnode );
    int natom=ff.natoms;
    int n=ff.nvecs;
    int* aneighs=new int[n*4];
    for(int i=0; i<n; i++){
        int i4=i*4;
        aneighs[i4+0]=-1; 
        aneighs[i4+1]=-1; 
        aneighs[i4+2]=-1; 
        aneighs[i4+3]=-1;
    }
    if(bMMFF){
        for(int i=0; i<ff.nnode; i++){
            int i4=i*4;
            for(int j=0; j<4; j++){
                int ngi=ff.aneighs[i4+j];
                if(ngi<0){ ngi= natom-ngi-1; }  // pi-bonds
                aneighs[i4+j]=ngi;
                if( (ngi>=ff.nnode) &&  (ngi<ff.natoms) ){   // back-neighbor
                    aneighs[ngi*4+0]=i;
                }
            }
        }
    }
    ff.makePiNeighs( aneighs+ff.natoms*4 );
    //if(verbosity>1) 
    for(int i=0; i<n; i++){ printf( "neighs[%i] (%i,%i,%i,%i) atyp %i R %g \n", i, aneighs[i*4+0],aneighs[i*4+1], aneighs[i*4+2],aneighs[i*4+3], nbmol.atypes[i], nbmol.REQs[i].x ); }
    ocl.upload( ocl.ibuff_neighs, aneighs );
    nbmol.neighs = (Quat4i*)aneighs;
    //delete [] aneighs;
    //return aneighs;
    //exit(0);
}

void makeBackNeighs( Quat4i* aneighs ){
    //int nbk = n+npi;
    bkneighs=new int[ff.nvecs*4];
    for(int i=0; i<ff.nvecs*4; i++){ bkneighs[i]=-1; };
    for(int ia=0; ia<ff.nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = aneighs[ia].array[j];
            if( ja<0 )continue;
            bool ret = addFirstEmpty( bkneighs+ja*4, 4, ia*4+j, -1 );
            if(!ret){ printf("ERROR in MolWorld_sp3_ocl::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
    for(int i=0; i<ff.nvecs; i++){
        printf( "bkneigh[%i] (%i,%i,%i,%i) \n", i, bkneighs[i*4+0], bkneighs[i*4+1], bkneighs[i*4+2], bkneighs[i*4+3] ); 
    }
    checkBkNeighCPU();
}

void checkBkNeighCPU(){
    int n=ff.nvecs;
    Quat4f* neighForce = new Quat4f[n*4];
    // ------ Write
    for(int i=0; i<n*4; i++){ neighForce[i].w=-1; };
    for(int i=0; i<n; i++){
        Quat4i ngs = nbmol.neighs[i];
        neighForce[i*4+0] = Quat4f{0,0,0,ngs.x+0.2f};
        neighForce[i*4+1] = Quat4f{0,0,0,ngs.y+0.2f};
        neighForce[i*4+2] = Quat4f{0,0,0,ngs.z+0.2f};
        neighForce[i*4+3] = Quat4f{0,0,0,ngs.w+0.2f};
    }
    // ------- Read
    bool match=true;
    Quat4i* bkngs = (Quat4i*)bkneighs;
    for(int i=0; i<n; i++){
        //float fi = i*1.f;
        Quat4i ng = bkngs[i];
        if(ng.x>0) match &= ((int)neighForce[ng.x].w == i); 
        if(ng.y>0) match &= ((int)neighForce[ng.y].w == i); 
        if(ng.z>0) match &= ((int)neighForce[ng.z].w == i); 
        if(ng.w>0) match &= ((int)neighForce[ng.w].w == i);
        //printf( "atom[%i] bkneighs(%i,%i,%i,%i) forMe(%i,%i,%i,%i) \n", i, ng.x,ng.y,ng.z,ng.w,  (int)(neighForce[ng.x].w-0.1),(int)(neighForce[ng.y].w-0.1),(int)(neighForce[ng.z].w-0.1),(int)(neighForce[ng.w].w-0.1) );
    }
    //printf( "!!! %i\n", (int)(0-0.1) );
    if(!match){ printf("ERROR in checkBkNeighCPU \n"); exit(0); }
    //exit(0);
}


virtual void init( bool bGrid ) override  {
    MolWorld_sp3::init(bGrid);
    mol2ocl();
    /*
    gridFF.grid.printCell();
    ocl.setNs(3, gridFF.grid.n.array );
    v2f4( gridFF.grid.pos0,ocl.pos0); 
    ocl.setGridShape( gridFF.grid.dCell );
    init_ocl();
    printf( "... MolWorld::init() DONE \n");
    */
    bGridFF=false;
    bOcl   =false;
}

double eval_gridFF_ocl( int n, Vec3d* ps,             Vec3d* fs ){ 
    //printf("eval_gridFF_ocl() \n");
    pack  ( n, ps, q_ps, sq(gridFF.Rdamp) );
    ocl.getNonBondForce_GridFF( n, (float4*)q_ps, 0, (float4*)q_fs, 0 );
    double E=unpack_add( n, fs, q_fs );
    //unpack( n, fs, q_fs );
    //double E=0;
    //for(int i=0; i<n; i++){ E+=q_fs[i].e; }; // sum energy
    return E;
}

double eval_MMFF_ocl( int niter, int n, Vec3d* ps, Vec3d* fs ){ 
    //printf( " ======= eval_MMFF_ocl() \n" );
    //pack  ( n, ps, q_ps, sq(gridFF.Rdamp) );
    OCLtask* task_getF = ocl.getTask("getMMFFsp3");
    OCLtask* task_move = ocl.getTask("gatherForceAndMove");
    OCLtask* task_pi0s = ocl.getTask("updatePiPos0");
    OCLtask* task_pipi = ocl.getTask("evalPiPi");
    //OCLtask* task_gff  = ocl.setup_getNonBondForce_GridFF( 0, n);
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    ocl.setup_gatherForceAndMove( ff.nvecs,  ff.natoms, task_move );
    ocl.setup_getMMFFsp3        ( ff.natoms, ff.nnode,  task_getF );
    ocl.setup_updatePiPos0      ( ff.natoms, ff.npi,    task_pi0s );
    ocl.setup_evalPiPi          ( ff.natoms, ff.npi,    task_pipi );
    //pack  ( n, ps, q_ps, sq(gridFF.Rdamp) );
    //ocl.upload( ocl.ibuff_atoms,  (float4*)q_ps, n ); // Note - these are other atoms than used for makeGridFF()
    //ocl.upload( ocl.ibuff_coefs,   coefs,  na);
    for(int i=0; i<niter; i++){
        //task_gff->enque_raw();
        task_getF->enque_raw();
        //task_pi0s->enque_raw();
        task_pipi->enque_raw();
        task_move->enque_raw();
    }
    printf( "ocl.download(n=%i) \n", n );
    ocl.download( ocl.ibuff_aforces, q_fs,    n );
    ocl.download( ocl.ibuff_atoms,   q_ps,    n );
    ocl.download( ocl.ibuff_pi0s,    ff.pi0s, ff.npi );
    ocl.finishRaw();
    printf( "*pi0s=%li\n", ff.pi0s ); for(int i=0; i<ff.npi; i++){printf( "CPU pi0s[%i](%g,%g,%g,%g)\n", i, ff.pi0s[i].x,ff.pi0s[i].y,ff.pi0s[i].z,ff.pi0s[i].w );}
    unpack             ( n, ps, q_ps );
    double E=unpack_add( n, fs, q_fs );
    //unpack( n, fs, q_fs );
    //double E=0;
    //for(int i=0; i<n; i++){  printf( "atom[%i] pos(%g,%g,%g) force(%g,%g,%g)\n",  i, ps[i].x,ps[i].y,ps[i].z,    fs[i].x,fs[i].y,fs[i].z ); }; // sum energy
    return E;
}

void eval(){
    //SDL_Delay(500);
    //SDL_Delay(100);
    //printf("#======= MDloop[%i] \n", nloop );
    double E=0;
    //bGPU_MMFF=false;

    //ff.doBonds  = false;  
    //ff.doNeighs =false;  
    //ff.doAngles =false;
    //ff.doPiSigma=false;
    ff.doPiPiI  =false;
    ff.doPiPiT  =false;

    //bSurfAtoms= false;
    if(bGPU_MMFF){
        printf( " ### GPU \n" );
        //E += ff.eval();
        //for(int i=0; i<ff.natoms; i++){ printf("CPU atom[%i] f(%g,%g,%g) \n", i, ff.fapos[i].x,ff.fapos[i].y,ff.fapos[i].z ); };
        //eval_MMFF_ocl( 1, nbmol.n, nbmol.ps, nbmol.fs );
        eval_MMFF_ocl( 1, ff.natoms+ff.npi, ff.apos, ff.fapos );
        //for(int i=0; i<ff.natoms; i++){ printf("OCL atom[%i] f(%g,%g,%g) \n", i, ff.fapos[i].x,ff.fapos[i].y,ff.fapos[i].z ); };
        //exit(0);
    }else{
        printf( " ### CPU \n" );
        if(bMMFF){ E += ff.eval();  } 
        else     { VecN::set( nbmol.n*3, 0.0, (double*)nbmol.fs );  }
        //if(bNonBonded){ E+= nff   .evalLJQ_pbc( builder.lvec, {1,1,1} ); }
        //bGridFF=true;
        //bOcl   =true;
        if(bSurfAtoms){ 
            if  (bGridFF){ 
                if(bOcl){ E+= eval_gridFF_ocl(nbmol.n, nbmol.ps,             nbmol.fs ); } 
                else    { E+= gridFF.eval    (nbmol.n, nbmol.ps, nbmol.PLQs, nbmol.fs ); 
                        E+= nbmol.evalNeighs();
                    }
            }else { 
                E+= nbmol.evalNeighs();   // Non-bonded interactions between atoms within molecule
              //E+= nbmol.evalMorse   (surf, false,                   gridFF.alpha, gridFF.Rdamp );
                E+= nbmol.evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp );
              //E+= nbmol.evalMorsePLQ( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp ); 
            }
        }
    }
    //for(int i=0; i<ff.natoms; i++){ printf("atom[%i] f(%g,%g,%g) \n", i, ff.fapos[i].x ,ff.fapos[i].y ,ff.fapos [i].z ); };
    //for(int i=0; i<ff.npi   ; i++){ printf("pvec[%i] f(%g,%g,%g) \n", i, ff.fpipos[i].x,ff.fpipos[i].y,ff.fpipos[i].z ); };
}

virtual void MDloop( int nIter, double Ftol = 1e-6 ) override {
    //printf( "MolWorld_sp3_ocl::MDloop(%i) bGridFF %i bOcl %i bMMFF %i \n", nIter, bGridFF, bOcl, bMMFF );
    //bMMFF=false;
    if(bMMFF)ff.cleanAll();
    if(bChargeUpdated){  REQs2ocl(); }
    for(int itr=0; itr<nIter; itr++){
        eval();
        //for(int i=0; i<nbmol.n; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fs[i].x,nbmol.fs[i].y,nbmol.fs[i].z ); }
        ckeckNaN_d( nbmol.n, 3, (double*)nbmol.fs, "nbmol.fs" );
        if(ipicked>=0){
             float K = -2.0;
             Vec3d f = getForceSpringRay( ff.apos[ipicked], pick_hray, pick_ray0, K );
             ff.fapos[ipicked].add( f );
        };
        if( !bGPU_MMFF){ // update atomic positions
            //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
            //opt.move_GD(0.001);
            //opt.move_LeapFrog(0.01);
            //opt.move_MDquench();
            opt.move_FIRE();
        }
        double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
        nloop++;
    }
    bChargeUpdated=false;
}


virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0} )override{
    printf( "MolWorld_sp3_ocl::initGridFF() \n");
    if(verbosity>0)printf("MolWorld_sp3::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    sprintf(tmpstr, "%s.lvs", name );
    if( file_exist(tmpstr) ){ 
        gridFF.grid.loadCell( tmpstr, gridStep );
        if(bGrid){
            gridFF.grid.center_cell( cel0 );
            bGridFF=true;
            gridFF.bindSystem(surf.n, surf.atypes, surf.ps, surf.REQs );
            if( isnan(z0) ){  z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
            gridFF.grid.pos0.z=z0;
            if(verbosity>1)gridFF.grid.printCell();
            gridFF.allocateFFs();
            //gridFF.tryLoad( "FFelec.bin", "FFPauli.bin", "FFLondon.bin", false, {1,1,0}, bSaveDebugXSFs );
            //gridFF.tryLoad( "FFelec.bin", "FFPauli.bin", "FFLondon.bin", false, nPBC, bSaveDebugXSFs );

            long T0 = getCPUticks();
            {// OpenCL-accelerated   GridFF initialization
              gridFF.grid.printCell();
              ocl.setNs(3, gridFF.grid.n.array );
              v2f4( gridFF.grid.pos0,ocl.pos0); 
              ocl.setGridShape( gridFF.grid.dCell );
              init_ocl();
            }
            printf( ">>time(init_ocl;GridFF_ocl): %g [s] \n", (getCPUticks()-T0)*tick2second  );
            bGridFF   =true; 
            //bSurfAtoms=false;
        }
    }else{ 
        bGridFF=false; 
        printf( "WARRNING!!! GridFF not initialized because %s not found\n", tmpstr );
    }
}

virtual void swith_method()override{ 
    bGPU_MMFF=!bGPU_MMFF;    bOcl=bGPU_MMFF;
    /*
    imethod=(imethod+1)%2; 
    switch (imethod){
        case 0: bGridFF=0; bOcl=0; break;
        case 1: bGridFF=1; bOcl=1; break;
    }
    */
}

virtual char* info_str   ( char* str=0 ){ if(str==0)str=tmpstr; sprintf(str,"bGridFF %i bOcl %i \n", bGridFF, bOcl ); return str; }

}; // class MolWorld_sp3_ocl

#endif