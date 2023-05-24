
constexpr int ntmpstr=2048;
char tmpstr[ntmpstr];

int verbosity = 1;
int idebug    = 0;
double tick2second=1e-9;

#include "testUtils.h"
#include "MolWorld_sp3_simple.h"

// ============ Global Variables

MolWorld_sp3_simple W;

//============================


#include "libMMFF.h"
#include "libUtils.h"

extern "C"{

// void setVerbosity( int verbosity_, int idebug_ ){
//     verbosity = verbosity_;
//     idebug    = idebug_;
// }

void init_buffers(){
    buffers .insert( { "apos",   (double*)W.nbmol.apos  } );
    buffers .insert( { "fapos",  (double*)W.nbmol.fapos } );
    if(W.bMMFF){
        buffers .insert( { "DOFs",      W.ffl.DOFs  } );
        buffers .insert( { "fDOFs",     W.ffl.fDOFs } );
        buffers .insert( { "pipos",  (double*)W.ffl.pipos   } );
        buffers .insert( { "fpipos", (double*)W.ffl.fpipos } );
        ibuffers.insert( { "neighs",   (int*)W.ffl.neighs  } );
    }
    ibuffers.insert( { "ndims",    &W.ffl.nDOFs } );
    buffers .insert( { "Es",       &W.ffl.Etot  } );
    //printf( "MMFFsp3_lib::init_buffers() nDOFs=%i nnode=%i ncap=%i nvecs=%i \n", W.ffl.nDOFs, W.ffl.nnode, W.ffl.ncap, W.ffl.nvecs );
    //printf( "MMFFsp3_lib::init_buffers() nDOFs=%i nnode=%i ncap=%i nvecs=%i \n", (&W.ffl.nDOFs)[0], (&W.ffl.nDOFs)[1], (&W.ffl.nDOFs)[2], (&W.ffl.nDOFs)[3] );
    //ibuffers.insert( { "selection", W.manipulation_sel  } );
    //for( auto c : buffers ){ printf("buff>>%s<<\n", c.first.c_str() ); }
}

void* init( const char* xyz_name, const char* smile_name, int* nPBC, const char* sElementTypes, const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes, const char* sDihedralTypes ){
	//printf( "DEBUG nPBC(%i,%i,%i)\n", nPBC[0],nPBC[1],nPBC[2] );
    W.smile_name = smile_name;
	W.xyz_name   = xyz_name;
    W.nPBC       = *(Vec3i*)nPBC;
    W.tmpstr=tmpstr;
    W.params.init( sElementTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
	W.builder.bindParams(&W.params);
    W.init();
    //init_buffers();
    return &W;
}
  
// int  buildMolecule_xyz( const char* xyz_name, bool bEpairs, double fAutoCharges, bool bAutoTypes, bool bRelaxPi ){  
//     W.builder.bDummyEpair=bEpairs; W.builder.bAutoTypes=bAutoTypes; W.bRelaxPi=bRelaxPi; W.fAutoCharges=fAutoCharges;  
//     //printf( "buildMolecule_xyz W.builder.bDummyEpair=%i bEpairs=%i \n", W.builder.bDummyEpair, bEpairs );
//     return W.buildMolecule_xyz( xyz_name ); 
// }
  
int    run( int nstepMax, double dt, double Fconv, int ialg, double* outE, double* outF ){
    //W.rum_omp_ocl( nstepMax, dt, Fconv, 1000.0, 1000 ); 
    return W.run(nstepMax,dt,Fconv,ialg,outE,outF);  
}

void scanHBond( const char* fname, int n, double dl, double l0, int ifrag1, int ifrag2, int i1a,int i1b, int i2a,int i2b, bool isDonor1, bool isDonor2, double* ups_ ){
    FILE* fout = fopen(fname,"w");
    const MM::Fragment& frag1 = W.builder.frags[ifrag1];
    const MM::Fragment& frag2 = W.builder.frags[ifrag2];
    

    //Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    //Vec3d fw2 = W.builder.atoms[i2b].pos - W.builder.atoms[i2a].pos;  fw2.normalize();    
    // W.builder.orient_atoms( fw1, up1, W.builder.atoms[i1a].pos,  Vec3dZero,            frag1.atomRange.a, frag1.atomRange.b );
    // W.builder.orient_atoms( fw2, up2, W.builder.atoms[i2a].pos,  Vec3d{ 0.0,0.0,l0 },  frag2.atomRange.a, frag2.atomRange.b );

    // Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    // Vec3d fw2 = W.builder.atoms[i2a].pos - W.builder.atoms[i2b].pos;  fw2.normalize();  // reverse direction    
    // W.builder.orient_atoms( fw1, up1, W.builder.atoms[i1b].pos,  Vec3dZero,            frag1.atomRange.a, frag1.atomRange.b ); // donor    => center is capping hydrogen [i1b]
    // W.builder.orient_atoms( fw2, up2, W.builder.atoms[i2a].pos,  Vec3d{ 0.0,0.0,l0 },  frag2.atomRange.a, frag2.atomRange.b ); // acceptor => center is node atom        [i2a]
    //printf( "isDonor %i %i \n", isDonor1, isDonor2 );
    Vec3d p1,p2;
    if(isDonor1){ p1=W.builder.atoms[i1b].pos; }else{ p1=W.builder.atoms[i1a].pos; }
    if(isDonor2){ p2=W.builder.atoms[i2b].pos; }else{ p2=W.builder.atoms[i2a].pos; }
    Vec3d fw1 = W.builder.atoms[i1b].pos - W.builder.atoms[i1a].pos;  fw1.normalize();
    Vec3d fw2 = W.builder.atoms[i2b].pos - W.builder.atoms[i2a].pos;  fw2.normalize();   fw2.mul(-1.); // reverse direction   

    Vec3d up1,up2;
    if(ups_){ Vec3d* ups=(Vec3d*)ups_; up1=ups[0]; up2=ups[1]; up1.normalize(); up2.normalize(); }
    else{  // auto-up-vector (make sure it is not colinear with fw)
        if( fabs(fabs(fw1.z)-1.0) < 0.1 ){ up1=Vec3dY; }else{ up1=Vec3dZ; } 
        if( fabs(fabs(fw2.z)-1.0) < 0.1 ){ up2=Vec3dY; }else{ up2=Vec3dZ; }
    }
    printf( "fw1(%g,%g,%g) |fw1.z-1|=%g up1(%g,%g,%g) p1(%g,%g,%g)\n", fw1.x,fw1.y,fw1.z, fabs(fabs(fw1.z)-1.0),   up1.x,up1.y,up1.z,    p1.x,p1.y,p1.z );
    printf( "fw2(%g,%g,%g) |fw2.z-1|=%g up2(%g,%g,%g) p2(%g,%g,%g)\n", fw2.x,fw2.y,fw2.z, fabs(fabs(fw2.z)-1.0),   up2.x,up2.y,up2.z,    p2.x,p2.y,p2.z );

    W.builder.orient_atoms( fw1, up1, p1,  Vec3d{ 0.0,0.0,0  }, frag1.atomRange.a, frag1.atomRange.b ); // donor    => center is capping hydrogen [i1b]
    W.builder.orient_atoms( fw2, up2, p2,  Vec3d{ 0.0,0.0,l0 }, frag2.atomRange.a, frag2.atomRange.b ); // acceptor => center is node atom        [i2a]

    Vec3d d=Vec3d{ 0.0,0.0, dl };
    for(int i=0; i<n+1; i++){
        sprintf(tmpstr, "#scanHBond[%i] l=%7.3f[A]", i, l0+dl*i );
        W.builder.write2xyz(fout, tmpstr );
        W.builder.move_atoms( d, frag2.atomRange.a, frag2.atomRange.b );
    }
    fclose(fout);
}

} // extern "C"
