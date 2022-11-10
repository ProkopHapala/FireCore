
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
#include "OCL_PP.h"


void pack  (int n, Vec3d* fs, Quat4f* qs, float K ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];  qs[i].e=K; } }
void unpack(int n, Vec3d* fs, Quat4f* qs          ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f;                 } }
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
}

void init_ocl(){
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

void mol2ocl(){
    ocl.buffers[ocl.ibuff_atoms].release();
    ocl.buffers[ocl.ibuff_coefs].release();
    // ---- Prepare for sampling atoms
    ocl.initAtomsForces( nbmol.n, bGPU_MMFF );
    printf( "buffs atoms, coefs, aforces, neighs: %i %i %i %i \n", ocl.ibuff_atoms, ocl.ibuff_coefs, ocl.ibuff_aforces, ocl.ibuff_neighs );
    q_ps= new Quat4f[nbmol.n];
    q_fs= new Quat4f[nbmol.n];
    // ---- nbmol coefs
    REQs2ocl();
    makeOCLNeighs( );
    if(bGPU_MMFF){
      makeBackNeighs( nbmol.n, nbmol.neighs );
      
    }
}

void  makeOCLNeighs( ){
    printf("makeOCLNeighs() n %i nnode %i \n", nbmol.n, ff.nnode );
    int n=nbmol.n;
    int* aneighs=new int[n*4];
    for(int i=0; i<n; i++){
      int i4=i*4;
      aneighs[i4+0]=-1; aneighs[i4+1]=-1; aneighs[i4+2]=-1; aneighs[i4+3]=-1;
    }
    if(bMMFF){
      for(int i=0; i<ff.nnode; i++){
        int i4=i*4;
        for(int j=0; j<4; j++){
          int ngi=ff.aneighs[i4+j];
          aneighs[i4+j]=ngi;
          if(ngi>=ff.nnode){   // back-neighbor
            aneighs[ngi*4+0]=i;
          }
        }
      }
    }
    if(verbosity>1) for(int i=0; i<n; i++){ printf( "neighs[%i] (%i,%i,%i,%i) atyp %i R %g \n", i, aneighs[i*4+0],aneighs[i*4+1], aneighs[i*4+2],aneighs[i*4+3], nbmol.atypes[i], nbmol.REQs[i].x ); }
    ocl.upload( ocl.ibuff_neighs, aneighs );
    nbmol.neighs = (Quat4i*)aneighs;

    
    //delete [] aneighs;
    //return aneighs;
    //exit(0);
}

void makeBackNeighs( int n, Quat4i* aneighs ){
    bkneighs=new int[n*4];
    for(int i=0; i<n*4; i++){ bkneighs[i]=-1; };
    for(int ia=0; ia<n; ia++){
      for(int j=0; j<4; j++){   // 4 neighbors
        int ja = aneighs[ia].array[j];
        if( (ja<0)||(ja>=n) )continue;
        bool ret = addFirstEmpty( bkneighs+ja*4, 4, ia, -1 );
        if(!ret){ printf("ERROR in MolWorld_sp3_ocl::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
      };
    }
    for(int i=0; i<n; i++){
      printf( "bkneigh[%i] (%i,%i,%i,%i) \n", i, bkneighs[i*4+0], bkneighs[i*4+1], bkneighs[i*4+2], bkneighs[i*4+3] ); 
    }
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

virtual void MDloop( int nIter, double Ftol = 1e-6 ) override {
    //printf( "MolWorld_sp3_ocl::MDloop(%i) bGridFF %i bOcl %i bMMFF %i \n", nIter, bGridFF, bOcl, bMMFF );
    bool bMMFF_ = bMMFF;
    //bMMFF_=false;
    if(bMMFF_)ff.cleanAll();
    if(bChargeUpdated){  REQs2ocl(); }
    for(int itr=0; itr<nIter; itr++){
        
        //printf("#======= MDloop[%i] \n", nloop );
        double E=0;
        if(bMMFF_){ E += ff.eval();  } 
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
              E+= nbmol.evalNeighs();
            //E+= nbmol.evalMorse   (surf, false,                   gridFF.alpha, gridFF.Rdamp );
              E+= nbmol.evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp );
            //E+= nbmol.evalMorsePLQ( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp ); 
            }
        }
        //for(int i=0; i<nbmol.n; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fs[i].x,nbmol.fs[i].y,nbmol.fs[i].z ); }
        ckeckNaN_d( nbmol.n, 3, (double*)nbmol.fs, "nbmol.fs" );
        if(ipicked>=0){
             float K = -2.0;
             Vec3d f = getForceSpringRay( ff.apos[ipicked], pick_hray, pick_ray0, K );
             ff.fapos[ipicked].add( f );
        };
        
        //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
        //opt.move_GD(0.001);
        //opt.move_LeapFrog(0.01);
        //opt.move_MDquench();
        opt.move_FIRE();
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
  imethod=(imethod+1)%2; 
  switch (imethod){
    case 0: bGridFF=0; bOcl=0; break;
    case 1: bGridFF=1; bOcl=1; break;
  }
}

virtual char* info_str   ( char* str=0 ){ if(str==0)str=tmpstr; sprintf(str,"bGridFF %i bOcl %i \n", bGridFF, bOcl ); return str; }

}; // class MolWorld_sp3_ocl

#endif
