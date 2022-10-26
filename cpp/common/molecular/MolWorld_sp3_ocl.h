
#ifndef MolWorld_sp3_ocl_h
#define MolWorld_sp3_ocl_h

#include "MolWorld_sp3.h"
#include "OCL_DFT.h"
#include "OCL_PP.h"


void pack  (int n, Vec3d* fs, Quat4f* qs ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];   } }
void unpack(int n, Vec3d* fs, Quat4f* qs ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f; } }

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
void init_ocl(){
    printf( "init_ocl() \n" );
    ocl.initPP( "common_resources/cl" );

    Quat4f* atoms = pack( surf.n, surf.ps   );
    Quat4f* coefs = pack( surf.n, surf.REQs );
    ocl.makeGridFF( surf.n, (float4*)atoms, (float4*)coefs );
    delete [] atoms;
    delete [] coefs; 
    //upload(ibuff_atoms,atoms);
    //upload(ibuff_coefs,coefs);

    printf( "... init_ocl() END\n" );
}

void init( bool bGrid ){
    MolWorld_sp3::init(bGrid);
    init_ocl();
    printf( "... MolWorld::init() DONE \n");
}

double eval_gridFF_ocl( int n, Vec3d* ps,             Vec3d* fs ){ 
    pack  ( n, ps, q_ps );
    ocl.getNonBondForce_GridFF( n, (float4*)q_ps, 0, (float4*)q_fs );
    unpack( n, fs, q_fs );
    double E=0;
    for(int i=0; i<n; i++){ E+=q_fs[i].e; } // sum energy
    return E;
}

void MDloop( int nIter, double Ftol = 1e-6 ){
    ff.cleanAll();
    for(int itr=0; itr<nIter; itr++){
        //printf("#======= MDloop[%i] \n", nloop );
        double E=0;
		    E += ff.eval();
		//if(bNonBonded){ E+= nff   .evalLJQ_pbc( builder.lvec, {1,1,1} ); }
        bGridFF=false;
        if(bSurfAtoms){ 
            if   (bGridFF){ 
              if(bOcl){ E+= gridFF.eval    (nbmol.n, nbmol.ps, nbmol.PLQs, nbmol.fs ); } 
              else    { E+= eval_gridFF_ocl(nbmol.n, nbmol.ps,             nbmol.fs ); }
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
        opt.move_FIRE();
        double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
        nloop++;
    }
}

}; // class MolWorld_sp3_ocl

#endif
