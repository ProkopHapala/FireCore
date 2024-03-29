﻿
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
#include "OCL_PP.h"
#include "datatypes_utils.h"

// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_ocl : public MolWorld_sp3 { public:
    OCL_PP     ocl;
    MMFFf4     ff4;

    Quat4f * q_ps=0;
    Quat4f * q_fs=0; 
    Quat4f * q_vs=0; 
    int* bkneighs=0;

    bool bGPU_MMFF    = true;
    //bool bGPU_MMFFsp3 = true;
    bool bGPU_MMFFsp3 = false;
    bool bGPU_MMFFf4  = true;
    //bool bGPU_MMFFf4  = false;

    OCLtask* task_cleanF=0;
    OCLtask* task_NBFF=0; //ocl.getTask("getMMFFf4");
    OCLtask* task_MMFF=0; //ocl.getTask("getMMFFf4");
    OCLtask* task_move=0;    //ocl.getTask("gatherForceAndMove");
    OCLtask* task_pi0s=0;    //ocl.getTask("updatePiPos0");
    OCLtask* task_pipi=0;    //ocl.getTask("evalPiPi");

// ======== Functions

void surf2ocl(Vec3i nPBC, bool bSaveDebug=false){
    int err=0;
    int ncell = (nPBC.x*2+1) * (nPBC.y*2+1) * (nPBC.z*2+1); 
    int n = surf.natoms*ncell;
    printf( "surf2ocl() na(%i) = ncell(%i) * natom(%i)\n", n, ncell, surf.natoms );
    long T0=getCPUticks();
    Quat4f* atoms = new Quat4f[n];
    Quat4f* coefs = new Quat4f[n];
    double R2damp = sq(gridFF.Rdamp);
    int ii=0;
    for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
        Vec3d  p0 = gridFF.grid.cell.a*ia + gridFF.grid.cell.b*ib + gridFF.grid.cell.c*ic;
        for(int i=0; i<surf.natoms; i++){
            //printf("ii %i i %i iabc(%i,%i,%i)\n", ii, i, ia,ib,ic);
            atoms[ii].f=(Vec3f)(surf.apos[i]+p0);
            atoms[ii].e=R2damp; 
            coefs[ii]  =(Quat4f)surf.REQs[i];
            //coefs[ii].e=gridFF.alphaMorse;
            ii++;
        }
    }}}
    printf( ">>time(surf_to_GPU) %g \n", (getCPUticks()-T0)*tick2second );
    long T1=getCPUticks();
    ocl.makeGridFF( n, (float4*)atoms, (float4*)coefs );
    err |=  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.finish");
    printf( ">>time(ocl.makeGridFF() %g \n", (getCPUticks()-T1)*tick2second );
    delete [] atoms;
    delete [] coefs;
    if(bSaveDebug){
        ocl.download( ocl.itex_FE_Paul, gridFF.FFPaul );  
        ocl.download( ocl.itex_FE_Lond, gridFF.FFLond );
        ocl.download( ocl.itex_FE_Coul, gridFF.FFelec );
        err |=  ocl.finishRaw();    OCL_checkError(err, "surf2ocl.download.finish");
        printf( ">>time(ocl.makeGridFF(); ocl.download() %g \n", (getCPUticks()-T1)*tick2second );
        saveGridXsfDebug();
    }
    ocl.buffers[ocl.ibuff_atoms].release();
    ocl.buffers[ocl.ibuff_coefs].release();
}

/*
void init_ocl(){
    //  ToDo : This is probably wrong
    printf( "MolWorld_sp3_ocl::init_ocl()\n" );
    //for(int i=0; i<nbmol.natoms; i++){ nbmol.apos->addRandomCube(0.1); }; // Debug - displace atoms to test relaxation forces
    long T0 = getCPUticks();
    //ocl.initPP( "common_resources/cl" );    // Debug
    // ---- Init Surface force-field grid
    //surf2ocl( nPBC, true );
    surf2ocl( nPBC, false );                // Debug
    printf( ">>time(surf2ocl): %g [s] \n", (getCPUticks()-T1)*tick2second  );
    printf( "MolWorld_sp3_ocl::init_ocl() END\n" );
}
*/

void REQs2ocl(){
    Quat4f* q_cs= new Quat4f[nbmol.natoms];
    pack  ( nbmol.natoms, nbmol.REQs, q_cs );
    ocl.upload( ocl.ibuff_coefs, q_cs );    
    delete [] q_cs;
};

int ithNeigh(int ia, int ja){
    int* neighs = (int*)nbmol.neighs;
    int* ngs     = neighs + ia*4;
    for(int i=0; i<4; i++){
        if(ngs[i]==ja) return i; 
    }
    return -1;
}

void MMFFsp3_to_ocl(){
    //int n=nbmol.natoms;
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
    int n   = nbmol.natoms;
    int npi = 0;
    if(bGPU_MMFFsp3){ npi=ff.npi;    }
    if(bGPU_MMFFf4 ){ npi=ff4.nnode; }
    int n2 = n+npi; 
    // ---- Prepare for sampling atoms
    ocl.initAtomsForces( n, npi, ff.nnode, bGPU_MMFFsp3, bGPU_MMFFf4 );
    printf( "mol2ocl() n,pi(%i,%i) buffs atoms, coefs, aforces, neighs: %i %i %i %i \n", n,npi, ocl.ibuff_atoms, ocl.ibuff_coefs, ocl.ibuff_aforces, ocl.ibuff_neighs );
    q_ps= new Quat4f[n2];
    q_fs= new Quat4f[n2];
    // ---- nbmol coefs
    pack      ( n2, nbmol.apos,  q_ps, sq(gridFF.Rdamp) );
    ocl.upload( ocl.ibuff_atoms, q_ps );  
    REQs2ocl();
    if(bGPU_MMFFsp3 || bGPU_MMFFf4 ){
        Quat4f* q_vs = new Quat4f[n2];
        for(int i=0; i<n2; i++){ q_vs[i].set(0.); }
        ocl.upload( ocl.ibuff_avel,   q_vs );
        if( bGPU_MMFFsp3 ){
            makeOCLNeighs( );
            makeBackNeighs( nbmol.neighs );
            ocl.upload( ocl.ibuff_bkNeighs, bkneighs );
            MMFFsp3_to_ocl();
            makePi0s();
        }
        if( bGPU_MMFFf4 ){

            //makeOCLNeighs( );
            //makeBackNeighs( nbmol.neighs );
            //ocl.upload( ocl.ibuff_bkNeighs, bkneighs );

            ocl.upload( ocl.ibuff_neighs   , ff4.neighs    );
            ocl.upload( ocl.ibuff_neighCell, ff4.neighCell );
            ocl.upload( ocl.ibuff_bkNeighs , ff4.bkneighs   );

            ocl.upload( ocl.ibuff_MMpars  , ff4.apars );  
            ocl.upload( ocl.ibuff_BLs     , ff4.bLs   ); 
            ocl.upload( ocl.ibuff_BKs     , ff4.bKs   ); 
            ocl.upload( ocl.ibuff_Ksp     , ff4.Ksp   ); 
            ocl.upload( ocl.ibuff_Kpp     , ff4.Kpp   ); 
/*
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,         // 2    [natoms]
    __global float4*  fapos,        // 3    [natoms]     
    __global float4*  fneigh,       // 4    [nnode*4]
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 8  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
    const cl_Mat3 lvec,             // 12
    const cl_Mat3 invLvec           // 13
*/
        }
    }
}

void makePi0s(){
    if(ff.pi0s==0){ ff.pi0s=new Quat4f[ff.npi]; };
    ff.evalPi0s();
    ocl.upload( ocl.ibuff_pi0s, ff.pi0s ); 
}

void  makeOCLNeighs( ){
    printf("makeOCLNeighs() n %i nnode %i \n", nbmol.natoms, ff.nnode );
    int natom=ff.natoms;
    int n=ff.nvecs;
    int* neighs=new int[n*4];
    for(int i=0; i<n; i++){
        int i4=i*4;
        neighs[i4+0]=-1; 
        neighs[i4+1]=-1; 
        neighs[i4+2]=-1; 
        neighs[i4+3]=-1;
    }
    if(bMMFF){
        for(int i=0; i<ff.nnode; i++){
            int i4=i*4;
            for(int j=0; j<4; j++){
                int ngi=ff.neighs[i4+j];
                if(ngi<0){ ngi= natom-ngi-1; }  // pi-bonds
                neighs[i4+j]=ngi;
                if( (ngi>=ff.nnode) &&  (ngi<ff.natoms) ){   // back-neighbor
                    neighs[ngi*4+0]=i;
                }
            }
        }
    }
    ff.makePiNeighs( neighs+ff.natoms*4 );
    //if(verbosity>1) 
    for(int i=0; i<n; i++){ printf( "neighs[%i] (%i,%i,%i,%i) atyp %i R %g \n", i, neighs[i*4+0],neighs[i*4+1], neighs[i*4+2],neighs[i*4+3], nbmol.atypes[i], nbmol.REQs[i].x ); }
    ocl.upload( ocl.ibuff_neighs, neighs );
    nbmol.neighs = (Quat4i*)neighs;
    //delete [] neighs;
    //return neighs;
    //exit(0);
}

void makeBackNeighs( Quat4i* neighs ){
    //int nbk = n+npi;
    bkneighs=new int[ff.nvecs*4];
    for(int i=0; i<ff.nvecs*4; i++){ bkneighs[i]=-1; };
    for(int ia=0; ia<ff.nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = neighs[ia].array[j];
            if( ja<0 )continue;
            bool ret = addFirstEmpty( bkneighs+ja*4, 4, ia*4+j, -1 );
            if(!ret){ printf("ERROR in MolWorld_sp3_ocl::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
    //for(int i=0; i<ff.nvecs; i++){ printf( "bkneigh[%i] (%i,%i,%i,%i) \n", i, bkneighs[i*4+0], bkneighs[i*4+1], bkneighs[i*4+2], bkneighs[i*4+3] ); }
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

virtual void init() override  {
    ocl.init();
    ocl.makeKrenels_PP("common_resources/cl" );
    MolWorld_sp3::init();
    mol2ocl();
    bGridFF=false;
    bOcl   =false;

    // initialization of ff4 is here because parrent MolWorld_sp3 does not contain ff4 anymore 
    builder.toMMFFf4( ff4, true, bEpairs );  //ff4.printAtomParams(); ff4.printBKneighs(); 
    ff4.flipPis( Vec3fOne );
    ff4.setLvec((Mat3f)builder.lvec);
    ff4.makeNeighCells  ( nPBC );
}

// ======================== SETUP OCL KERNEL functions

void setup_MMFFsp3_ocl( ){
    //printf( " ======= eval_MMFF_ocl() \n" );
    //pack  ( n, ps, q_ps, sq(gridFF.Rdamp) );
    //OCLtask* task_gff  = ocl.setup_getNonBondForce_GridFF( 0, n);
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    //printf( "CPU lvec    (%g,%g,%g)(%g,%g,%g)(%g,%g,%g) \n", ff.lvecT.a.x,ff.lvecT.a.y,ff.lvecT.a.z,    ff.lvecT.b.x,ff.lvecT.b.y,ff.lvecT.b.z,   ff.lvecT.c.x,ff.lvecT.c.y,ff.lvecT.c.z  );
    //printf( "CPU lvec    (%g,%g,%g)(%g,%g,%g)(%g,%g,%g) \n", ff.lvec.a.x,ff.lvec.a.y,ff.lvec.a.z,    ff.lvec.b.x,ff.lvec.b.y,ff.lvec.b.z,   ff.lvec.c.x,ff.lvec.c.y,ff.lvec.c.z  );
    //printf( "CPU cl_lvec (%g,%g,%g)(%g,%g,%g)(%g,%g,%g) \n", ocl.cl_lvec.a.s[0],ocl.cl_lvec.a.s[1],ocl.cl_lvec.a.s[2],    ocl.cl_lvec.b.s[0],ocl.cl_lvec.b.s[1],ocl.cl_lvec.b.s[2],   ocl.cl_lvec.c.s[0],ocl.cl_lvec.c.s[1],ocl.cl_lvec.c.s[2]  );
    Mat3_to_cl( ff.lvec   , ocl.cl_lvec    );
    Mat3_to_cl( ff.invLvec, ocl.cl_invLvec );
    /*
    task_MMFF = ocl.getTask("getMMFFsp3");
    task_move = ocl.getTask("gatherForceAndMove");
    task_pi0s = ocl.getTask("updatePiPos0");
    task_pipi = ocl.getTask("evalPiPi");
    ocl.setup_gatherForceAndMove( ff.nvecs,  ff.natoms,      task_move );
    ocl.setup_getMMFFsp3        ( ff.natoms, ff.nnode, bPBC, task_MMFF );     
    ocl.setup_updatePiPos0      ( ff.natoms, ff.npi,         task_pi0s );
    ocl.setup_evalPiPi          ( ff.natoms, ff.npi,         task_pipi );
    */
    task_move = ocl.setup_gatherForceAndMove( ff.nvecs,  ff.natoms      );
    task_MMFF = ocl.setup_getMMFFsp3        ( ff.natoms, ff.nnode, bPBC );     
    task_pi0s = ocl.setup_updatePiPos0      ( ff.natoms, ff.npi         );
    task_pipi = ocl.setup_evalPiPi          ( ff.natoms, ff.npi         );
    //pack  ( n, ps, q_ps, sq(gridFF.Rdamp) );
    //ocl.upload( ocl.ibuff_atoms,  (float4*)q_ps, n ); // Note - these are other atoms than used for makeGridFF()
    //ocl.upload( ocl.ibuff_coefs,   coefs,  na);
} 

void setup_MMFFf4_ocl(){
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    Mat3_to_cl( ff.lvec   , ocl.cl_lvec    );
    Mat3_to_cl( ff.invLvec, ocl.cl_invLvec );
    //printf( "na %i nnode %i \n", ff4.natoms, ff4.nnode );
    /*
    task_MMFF    = ocl.getTask("getMMFFf4");
    task_move    = ocl.getTask("updateAtomsMMFFf4");
    task_cleanF  = ocl.getTask("cleanForceMMFFf4");
    ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode,       task_move   );   
    ocl.setup_getMMFFf4        ( ff4.natoms, ff4.nnode, bPBC, task_MMFF   );
    ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode,       task_cleanF );
    */
    task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode       );   
    task_MMFF   = ocl.setup_getMMFFf4        ( ff4.natoms, ff4.nnode, bPBC );
    task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, nPBC );
    task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode       );
}

void setup_NBFF_ocl(){
    printf("======= setup_NBFF_ocl \n");
    ocl.nDOFs.x=ff.natoms;
    ocl.nDOFs.y=ff.nnode;
    Mat3_to_cl( ff.lvec   , ocl.cl_lvec    );
    Mat3_to_cl( ff.invLvec, ocl.cl_invLvec );
    //task_NBFF    = ocl.getTask("getNonBond");
    //task_move    = ocl.getTask("updateAtomsMMFFf4");
    //task_cleanF  = ocl.getTask("cleanForceMMFFf4");
    task_cleanF = ocl.setup_cleanForceMMFFf4 ( ff4.natoms, ff4.nnode );
    task_move   = ocl.setup_updateAtomsMMFFf4( ff4.natoms, ff4.nnode ); 
    task_NBFF   = ocl.setup_getNonBond       ( ff4.natoms, ff4.nnode, gridFF.nPBC );
}

// ======================== EVAL OCL KERNEL functions

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

double eval_MMFFsp3_ocl( int niter, int n, Vec3d* ps, Vec3d* fs ){ 
    if( task_MMFF==0 )setup_MMFFsp3_ocl();
    for(int i=0; i<niter; i++){
        //task_gff->enque_raw();
        task_MMFF->enque_raw();
        //task_pi0s->enque_raw();
        task_pipi->enque_raw();
        task_move->enque_raw();
    }
    //printf( "ocl.download(n=%i) \n", n );
    ocl.download( ocl.ibuff_aforces, q_fs,    n );
    ocl.download( ocl.ibuff_atoms,   q_ps,    n );
    ocl.download( ocl.ibuff_pi0s,    ff.pi0s, ff.npi );
    ocl.finishRaw();
    //printf( "*pi0s=%li\n", ff.pi0s ); for(int i=0; i<ff.npi; i++){printf( "CPU pi0s[%i](%g,%g,%g,%g)\n", i, ff.pi0s[i].x,ff.pi0s[i].y,ff.pi0s[i].z,ff.pi0s[i].w );}
    unpack             ( n, ps, q_ps );
    double E=unpack_add( n, fs, q_fs );
    //unpack( n, fs, q_fs );
    //double E=0;
    //for(int i=0; i<n; i++){  printf( "atom[%i] pos(%g,%g,%g) force(%g,%g,%g)\n",  i, ps[i].x,ps[i].y,ps[i].z,    fs[i].x,fs[i].y,fs[i].z ); }; // sum energy
    return E;
}


double eval_MMFFf4_ocl( int niter ){ 
    //printf( " ======= eval_MMFFf4() \n" );
    /*
    // ========= To Compare CPU and GPU forces
    ff4.eval();
    printf("CPU AFTER assemble() \n"); ff4.printDebug(  false,false );
    unpack( ff4.natoms, ffl. apos, ff4. apos  );
    unpack( ff4.natoms, ffl.fapos, ff4.fapos  );
    // ---- Check Invariatns
    fcog  = sum ( ffl.natoms, ffl.fapos   );
    tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
    if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4 |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }else{ printf("eval_MMFFf4 |fcog| OK\n"); }
    //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4 |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }  // NOTE: torq is non-zero because pi-orbs have inertia
    */
    //printf( " ======= eval_MMFFf4_ocl() \n" );
    if( task_MMFF==0 )setup_MMFFf4_ocl();
    for(int i=0; i<niter; i++){
        task_cleanF->enque_raw();  // Debug: this should be solved inside  task_move->enque_raw();
        task_MMFF  ->enque_raw();
        //task_NBFF  ->enque_raw();
        /*
        { // Debug
        ocl.download( ocl.ibuff_aforces,    ff4.fapos , ff4.nvecs );
        ocl.download( ocl.ibuff_atoms,      ff4.apos  , ff4.nvecs );
        ocl.download( ocl.ibuff_neighForce, ff4.fneigh, ff4.nvecs );
        ocl.finishRaw();   
        printf("GPU BEFORE assemble() \n"); ff4.printDdebug( false,false );
        }
        */
        task_move->enque_raw(); 
    }
    //printf( "ocl.download(n=%i) \n", n );
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs );
    //for(int i=0; i<ff4.natoms; i++){  printf("CPU[%i] p(%g,%g,%g) f(%g,%g,%g) \n", i, ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z ); }
    ocl.finishRaw();                              //Debug
    //printf("GPU AFTER assemble() \n"); ff4.printDebug( false,false );
    unpack( ff4.natoms, ffl.  apos, ff4.  apos );
    unpack( ff4.natoms, ffl. fapos, ff4. fapos );
    unpack( ff4.nnode,  ffl. pipos, ff4. pipos );
    unpack( ff4.nnode,  ffl.fpipos, ff4.fpipos );

    // ---- Check Invariatns
    fcog  = sum ( ffl.natoms, ffl.fapos   );
    tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
    if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia

    return 0;
}


double eval_NBFF_ocl( int niter ){ 
    //printf( " ======= eval_NBFF_ocl() \n" );
    if( task_NBFF==0 ){
        /*
        for(int i=0; i<ff4.natoms; i++){
            ff4.neighs[i] = Quat4i{-1,-1,-1,-1};
            ffl.neighs[i] = Quat4i{-1,-1,-1,-1};
        }
        ocl.upload( ocl.ibuff_neighs,  ff4.neighs, ff4.natoms );
        */
        setup_NBFF_ocl();
    }
    //pack      ( ff4.natoms,        ffl.apos,  ff4.apos );
    //ocl.upload( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs );
    for(int i=0; i<niter; i++){

        { // Debug
            ffl  .cleanForce();
            // nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, ffl.lvec, ffl.nPBC, gridFF.Rdamp );
            nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp );
            //ffl.printDebug( false, false );
            fcog  = sum ( ffl.natoms, ffl.fapos   );
            tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
            if(  fcog.norm2()>1e-8 ){ printf("WARRNING: ffl.cleanForce() |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
            //opt.move_FIRE();
            //ffl.move_GD(0.01); return 0;
        }
        
        task_cleanF->enque_raw(); // Debug: this should be solved inside  task_move->enque_raw();
        task_NBFF  ->enque_raw();
        task_move  ->enque_raw(); //Debug
        
    }
    
    //printf( "ocl.download(n=%i) \n", n );
    ocl.download( ocl.ibuff_aforces, ff4.fapos, ff4.nvecs );
    ocl.download( ocl.ibuff_atoms,   ff4.apos , ff4.nvecs );
    //for(int i=0; i<ff4.natoms; i++){  printf("CPU[%i] p(%g,%g,%g) f(%g,%g,%g) \n", i, ff4.apos[i].x,ff4.apos[i].y,ff4.apos[i].z,  ff4.fapos[i].x,ff4.fapos[i].y,ff4.fapos[i].z ); }
    ocl.finishRaw();                              //Debug
    //ff4.move_GD( 0.01 );
        
    // ---- Compare to ffl
    bool ret=false;
    printf("### Compare ffl.fapos,  ff4.fapos   \n"); ret |= compareVecs( ff4.natoms, ffl.fapos,  ff4.fapos,  1e-4, true );
    if(ret){ printf("ERROR: GPU task_NBFF.eval() != ffl.nbmol.evalLJQs_ng4_PBC() => exit() \n"); exit(0); }else{ printf("CHECKED: GPU task_NBFF.eval() == ffl.nbmol.evalLJQs_ng4_PBC() \n"); }

    //printf("GPU AFTER assemble() \n"); ff4.printDebug( false,false );
    unpack( ff4.natoms, ffl.  apos, ff4.  apos );
    unpack( ff4.natoms, ffl. fapos, ff4. fapos );
    //opt.move_FIRE();
    
    
    //ff4.printDebug( false, false );
    // ============== CHECKS
    // ---- Check Invariatns
    fcog  = sum ( ffl.natoms, ffl.fapos   );
    tqcog = torq( ffl.natoms, ffl.apos, ffl.fapos );
    if(  fcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |fcog| =%g; fcog=(%g,%g,%g)\n", fcog.norm(),  fcog.x, fcog.y, fcog.z ); exit(0); }
    //if( tqcog.norm2()>1e-8 ){ printf("WARRNING: eval_MMFFf4_ocl |torq| =%g; torq=(%g,%g,%g)\n", tqcog.norm(),tqcog.x,tqcog.y,tqcog.z ); exit(0); }   // NOTE: torq is non-zero because pi-orbs have inertia
    
    //exit(0);
    return 0;
}


// ======================== OTHER

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
        //printf( " ### GPU \n" );
        //E += ff.eval();
        //for(int i=0; i<ff.natoms; i++){ printf("CPU atom[%i] f(%g,%g,%g) \n", i, ff.fapos[i].x,ff.fapos[i].y,ff.fapos[i].z ); };
        //eval_MMFF_ocl( 1, nbmol.natoms, nbmol.apos, nbmol.fapps );
        //eval_MMFFsp3_ocl( 1, ff.natoms+ff.npi, ff.apos, ff.fapos );
        eval_MMFFf4_ocl( 1 );
        //eval_NBFF_ocl  ( 1 );
        //for(int i=0; i<ff.natoms; i++){ printf("OCL atom[%i] f(%g,%g,%g) \n", i, ff.fapos[i].x,ff.fapos[i].y,ff.fapos[i].z ); };
        //exit(0);
    }else{
        //printf( " ### CPU \n" );
        if(bMMFF){ E += ff.eval();  } 
        else     { VecN::set( nbmol.natoms*3, 0.0, (double*)nbmol.fapos );  }
        //if(bNonBonded){ E+= nff   .evalLJQ_pbc( builder.lvec, {1,1,1} ); }
        //bGridFF=true;
        //bOcl   =true;

        if(bSurfAtoms){ 
            if  (bGridFF){ 
                if(bOcl){ E+= eval_gridFF_ocl(nbmol.natoms, nbmol.apos,             nbmol.fapos ); } 
                else    { E+= gridFF.eval    (nbmol.natoms, nbmol.apos, nbmol.PLQs, nbmol.fapos ); 
                        E += nbmol.evalLJQs_ng4( ffl.neighs );
                    }
            }else { 
                E += nbmol.evalLJQs_ng4( ffl.neighs );  // Non-bonded interactions between atoms within molecule
              //E+= nbmol.evalMorse   (surf, false,                   gridFF.alphaMorse, gridFF.Rdamp );
                E+= nbmol.evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );
              //E+= nbmol.evalMorsePLQ( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp ); 
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
        ckeckNaN_d( nbmol.natoms, 3, (double*)nbmol.fapos, "nbmol.fs" );
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

virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true )override{
    printf( "MolWorld_sp3_ocl::initGridFF() \n");
    if(verbosity>0)printf("MolWorld_sp3::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    sprintf(tmpstr, "%s.lvs", name );
    // if( file_exist(tmpstr) ){  gridFF.grid.loadCell( tmpstr, gridStep );  gridFF.bCellSet=true; }
    // if( !gridFF.bCellSet ){
    //     bGridFF=false; 
    //     printf( "WARRNING!!! GridFF not initialized because %s not found\n", tmpstr );
    //     return;
    // }
    gridFF.grid.center_cell( cel0 );
    bGridFF=true;
    gridFF.bindSystem(surf.natoms, surf.atypes, surf.apos, surf.REQs );
    if( isnan(z0) ){  z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
    gridFF.grid.pos0.z=z0;
    if(verbosity>1)gridFF.grid.printCell();
    gridFF.allocateFFs();
    //gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false, {1,1,0}, bSaveDebugXSFs );
    //gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false, nPBC, bSaveDebugXSFs );
    if(bAutoNPBC){  autoNPBC( gridFF.grid.cell, nPBC, 30.0 ); }
    
    long T0 = getCPUticks();
    {// OpenCL-accelerated   GridFF initialization
        gridFF.grid.printCell();
        ocl.setNs(3, gridFF.grid.n.array );
        v2f4( gridFF.grid.pos0,ocl.pos0); 
        ocl.setGridShape( gridFF.grid.dCell );
        //init_ocl();
        surf2ocl( nPBC );

        if(bSaveDebugXSFs)saveGridXsfDebug();

        //surf2ocl( nPBC, false );
    }
    printf( ">>time(init_ocl;GridFF_ocl): %g [s] \n", (getCPUticks()-T0)*tick2second  );
    bGridFF   =true; 
    //bSurfAtoms=false;
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
