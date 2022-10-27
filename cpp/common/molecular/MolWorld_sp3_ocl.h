
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
#include "OCL_PP.h"


void pack  (int n, Vec3d* fs, Quat4f* qs, float K ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];  qs[i].e=K; } }
void unpack(int n, Vec3d* fs, Quat4f* qs          ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f;                 } }

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


// ======================================
// class:        MolWorld_sp3_ocl
// ======================================

class MolWorld_sp3_ocl : public MolWorld_sp3 { public:
  OCL_PP     ocl;
  bool bOcl=false;

  Quat4f * q_ps=0;
  Quat4f * q_fs=0; 


// ======== Functions

void surf2ocl(Vec3i nPBC, bool bSaveDebug=false){
  int ncell = (nPBC.x*2+1) * (nPBC.y*2+1) * (nPBC.z*2+1); 
  int n = surf.n*ncell;
  printf( "surf2ocl() na(%i) = ncell(%i) * natom(%i)\n", n, ncell, surf.n );
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

  ocl.makeGridFF( n, (float4*)atoms, (float4*)coefs );
  delete [] atoms;
  delete [] coefs;
  if(bSaveDebug){
    ocl.download( ocl.itex_FE_Paul, gridFF.FFPauli );  
    ocl.download( ocl.itex_FE_Lond, gridFF.FFLondon );
    ocl.download( ocl.itex_FE_Coul, gridFF.FFelec );
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
    ocl.init();
    ocl.initPP( "common_resources/cl" );
    // ---- Init Surface force-field grid
    surf2ocl( nPBC, true );
    ocl.buffers[ocl.ibuff_atoms].release();
    ocl.buffers[ocl.ibuff_coefs].release();
    // ---- Prepare for sampling atoms
    ocl.initAtomsForces( nbmol.n );
    printf( "buffs atoms, coefs, aforces: %i %i %i \n", ocl.ibuff_atoms, ocl.ibuff_coefs, ocl.ibuff_aforces );
    q_ps= new Quat4f[nbmol.n];
    q_fs= new Quat4f[nbmol.n];
    // ---- nbmol coefs
    Quat4f* q_cs= new Quat4f[nbmol.n];
    pack  ( nbmol.n, nbmol.REQs, q_cs, gridFF.alpha );
    ocl.upload( ocl.ibuff_coefs, q_cs );
    delete [] q_cs;
    printf( "... init_ocl() END\n" );
}

virtual void init( bool bGrid ) override  {
    MolWorld_sp3::init(bGrid);
    ocl.setNs(3, gridFF.grid.n.array );
    v2f4( gridFF.grid.pos0,ocl.pos0); 
    ocl.setGridShape( gridFF.grid.dCell );
    init_ocl();
    printf( "... MolWorld::init() DONE \n");
}

double eval_gridFF_ocl( int n, Vec3d* ps,             Vec3d* fs ){ 
    printf("eval_gridFF_ocl() \n");
    pack  ( n, ps, q_ps, sq(gridFF.Rdamp) );  DEBUG
    ocl.getNonBondForce_GridFF( n, (float4*)q_ps, 0, (float4*)q_fs );  DEBUG
    unpack( n, fs, q_fs );  DEBUG
    double E=0;
    for(int i=0; i<n; i++){ E+=q_fs[i].e; }; // sum energy
    return E;
}

virtual void MDloop( int nIter, double Ftol = 1e-6 ) override {
    printf( "MolWorld_sp3_ocl::MDloop(%i) bGridFF %i bOcl %i \n", nIter, bGridFF, bOcl );
    ff.cleanAll();
    for(int itr=0; itr<nIter; itr++){
        //printf("#======= MDloop[%i] \n", nloop );
        double E=0;
        ff.cleanAll();
		    //E += ff.eval();
		    //if(bNonBonded){ E+= nff   .evalLJQ_pbc( builder.lvec, {1,1,1} ); }
        //bGridFF=true;
        //bOcl   =true;
        if(bSurfAtoms){ 
            if  (bGridFF){ 
              if(bOcl){ E+= eval_gridFF_ocl(nbmol.n, nbmol.ps,             nbmol.fs ); } 
              else    { E+= gridFF.eval    (nbmol.n, nbmol.ps, nbmol.PLQs, nbmol.fs ); }
            }
            //else    { E+= nbmol .evalMorse   (surf, false,                   gridFF.alpha, gridFF.Rdamp );  }
            else      { E+= nbmol .evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alpha, gridFF.Rdamp );  }
        }
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
        //opt.move_FIRE();
        double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
        nloop++;
    }
}


virtual void swith_gridFF()override   { bGridFF=!bGridFF; if(bGridFF){ bOcl=!bOcl; }; };
virtual char* info_str   ( char* str=0 ){ if(str==0)str=tmpstr; sprintf(str,"bGridFF %i bOcl %i \n", bGridFF, bOcl ); return str; }

}; // class MolWorld_sp3_ocl

#endif
