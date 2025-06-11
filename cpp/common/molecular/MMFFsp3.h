
#ifndef MMFFsp3_h
#define MMFFsp3_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Forces.h"
#include "molecular_utils.h"

// ======================
// ====   MMFFsp3
// ======================

class MMFFsp3{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0,nnode=0,ncap=0,npi=0,nbonds=0,nvecs=0,ne=0,ie0=0;
    bool bPBC=false;
    double Etot,Eb,Ea, Eps,EppT,EppI;

    double * DOFs  = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 };
    double Kpipi = 0.25;
    int  ipi0=0;
    int   * atype=0;
    Vec3d * apos=0;
    Vec3d * fapos=0;
    Vec3d * pipos=0;
    Vec3d * fpipos=0;
    //Vec3d * cappos=0;
    //Vec3d * fcappos=0;
    Quat4f* pi0s=0;

    // --- Parameters
    Vec2i  * bond2atom = 0;
    double * bond_l0   = 0;  // [A]
    double * bond_k    = 0;  // [eV/A] ?
    double * bond_kPi  = 0;  // [eV/A] ?
    Vec3d  * pbcShifts = 0;  // [A]

    bool    bSubtractAngleNonBond=false;
    Quat4d*  REQs=0;  // this is used only when bSubtractAngleNonBond==true

    bool    bPBCbyLvec  =false;
    Mat3d   invLvec, lvec;

    int*    neighs = 0;  // [natom*nneigh_max] neigbors of atom
    int*    abonds  = 0;  // [natom*nneigh_max] bonds    of atom
    //double* Kneighs = 0;  // [natom*nneigh_max] angular stiffness for given neighbor
    Quat4d* NeighParams=0;

    // NOTE - For the moment we do not use temporary bonds - to keep things simple
    // --- Axuliary variables
    //double * lbond  = 0;   // bond lengths
    //Vec3d  * hbond  = 0;   // normalized bond unitary vectors

    bool doBonds  =true;
    bool doNeighs =true;
    bool doPiPiI  =true;
    bool doPiPiT  =true;
    bool doPiSigma=true;
    bool doAngles =true;
    bool doEpi    =true; 
    int nevalPiPiI  =0;
    int nevalPiPiT  =0;
    int nevalPiSigma=0;
    int nevalAngles =0;


    bool bDEBUG_plot=0;
    int iDEBUG_n=0,iDEBUG_pick=0;

void realloc( int nnode_, int nbonds_, int npi_, int ncap_, bool bNeighs=true ){
    nnode=nnode_; ncap=ncap_; nbonds=nbonds_; npi=npi_; 
    natoms=nnode_+ncap_; 
    nvecs = natoms + npi;
    nDOFs = nvecs*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    _realloc( DOFs     , nDOFs );
    _realloc( fDOFs    , nDOFs );
    //printf( "MMFFsp3::realloc() 1 \n" );
    ipi0=natoms;
    _realloc( atype, natoms );
    apos   = (Vec3d*)DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    //printf( "MMFFsp3::realloc() 2 \n" );
    _realloc( bond2atom , nbonds );
    _realloc( bond_l0   , nbonds );
    _realloc( bond_k    , nbonds );
    _realloc( bond_kPi  , nbonds );    for(int i=0; i<nbonds; i++){ bond_kPi[i]=0; }
    //_realloc( lbond     , nbonds );
    //_realloc( hbond     , nbonds );
    //printf( "MMFFsp3::realloc() 3 \n" );
    if(bNeighs){
        int nn = nnode*nneigh_max;
        _realloc( neighs, nn );
        _realloc( abonds,  nn );
        //_realloc( Kneighs, nn );
        _realloc( NeighParams, nnode );
        for(int i=0; i<nnode; i++){ NeighParams[i]=default_NeighParams; }
        //_realloc( Kpis       , nnode );
    }
    //printf( "MMFFsp3::realloc() DONE \n" );
}

void dealloc(  ){
    nnode=0; ncap=0; nbonds=0; npi=0; natoms=0; nvecs=0; nDOFs=0;
    _dealloc( DOFs      );
    _dealloc( fDOFs     );
    _realloc( atype, natoms );
    apos   = 0 ;
    fapos  = 0;
    pipos  = 0;
    fpipos = 0;
    _dealloc( bond2atom );
    _dealloc( bond_l0    );
    _dealloc( bond_k     );
    _dealloc( bond_kPi   );

    _dealloc( neighs );
    _dealloc( abonds );
    //_dealloc( Kneighs, nn );
    _dealloc( NeighParams );
    //_dealloc( Kpis      );
    
    //printf( "MMFFsp3::realloc() DONE \n" );
}

void initPBC(){
    _realloc(  pbcShifts, nbonds );
    for(int i=0; i<nbonds; i++){ pbcShifts[i]=Vec3dZero; }
}

void cleanEnergy   (){ Etot=0;Eb=0;Ea=0;Eps=0;EppT=0;EppI=0; };
void cleanAtomForce(){ for(int i=0; i<natoms; i++){ fapos [i].set(0.0); } }
void cleanPiForce  (){ for(int i=0; i<npi;    i++){ fpipos[i].set(0.0); } }
void normalizePi   (){ for(int i=0; i<npi;    i++){ pipos [i].normalize(); } }

void cleanAll(){
    cleanEnergy();
    cleanAtomForce();
    cleanPiForce();
}

// ============== Evaluation

void setLvec(const Mat3d& lvec_){
    lvec=lvec_; lvec.invert_T_to( invLvec );
    //lvecT.setT(lvec_);
    //lvecT.invert_to( invLvec );
}

int i_DEBUG = 0;

inline double evalSigmaSigma_dist( int ia, int ing, int jng, double K ){
    nevalAngles++;
    //K = -1.0;
    //  E = -K*|pi-pj|
    // ToDo: This may be made more efficient if weMMFFmini store hbonds
    Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();  // j->i
    Vec3d hi = apos[ing] - apos[ia];    double ri  = hi.normalize();
    Vec3d hj = apos[jng] - apos[ia];    double rj  = hj.normalize();
    Vec3d fi  = d* K;
    Vec3d fj  = d*-K;
    Vec3d fiI = hi* hi.dot(fi);   
    Vec3d fjI = hj* hj.dot(fj);   
    //printf(  "ang[%i|%i,%i] c1 %g c2 %g \n", ia, ing, jng, hi.dot(fi), hj.dot(fj)  );
    //opengl1renderer.color3f(1.,0.,0.);
    //Draw3D::drawVecInPos(  fi, apos[ing] );
    //Draw3D::drawVecInPos(  fj, apos[jng] );
    //Draw3D::drawVecInPos(  fjT+fiT, apos[ia] );
    //  force repel (ing,jng)
    //  but it is constrained so that is does not affect bond (ia,ing), (ia,jng)
    //  
    fapos[ia] .add(fiI);
    fapos[ia] .add(fjI);
    fi.sub(fiI); fapos[ing].add(fi);   // Why not  fiT ?
    fj.sub(fjI); fapos[jng].add(fj);   // Why not  fjT ?

    //fapos[ia] .sub(fiT);
    //fapos[ia] .sub(fjT);
    //fapos[ing].add(fiT);
    //fapos[jng].add(fjT);
    return K*rij;
}

inline double evalSigmaSigma_cos(  int ia, int ing, int jng, double K, double c0 ){
    nevalAngles++;
    //Vec3d h1,h2;
    //Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();
    Vec3d h1 = apos[ing] - apos[ia];    double ir1 = 1/h1.normalize();
    Vec3d h2 = apos[jng] - apos[ia];    double ir2 = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double c_ = c-c0;
    double E    =  K*c_*c_;
    double fang = -K*c_*2;
    hf1.mul( fang*ir1 );
    hf2.mul( fang*ir2 );
    //if(bDEBUG_plot){
    // //if(ia==7){
    //     //printf( "evalAngle_cos [%i|%i,%i,%i] c(%g,%g) \n", iDEBUG_pick, ia, ing, jng, hf1.dot(h1), hf1.dot(h1)  );
    //     //printf( "evalAngle_cos [%i|%i,%i,%i] c %g c_ %g fang_ %g E %g \n", iDEBUG_pick, ia, ing, jng, c, c_, fang, E );
    //     //printf( "evalAngle_cos bDEBUG_plot [%i|%i,%i] iDEBUG_pick %i \n", ia, ing, jng, iDEBUG_pick );
    //     opengl1renderer.color3f(1.,0.,0.);
    //     Draw3D::drawVecInPos(  hf1, apos[ing] );
    //     Draw3D::drawVecInPos(  hf2, apos[jng] );
    //     Draw3D::drawVecInPos(  (hf1+hf2)*-1.0, apos[ia] );
    // }
    fapos[ing].add( hf1     );
    fapos[jng].add( hf2     );
    //fapos[ia ].sub( hf1+hf2 );MMFFmini
    Vec3d fei=hf1+hf2; fapos[ia ].sub( fei );
    
    //if(ia==0) printf( "CPU atom[%i|%i,%i] c %g h1(%g,%g,%g) h2(%g,%g,%g) | hf1(%g,%g,%g) hf2(%g,%g,%g) \n", ia, ing,jng, c, h1.x,h1.y,h1.z,  h2.x,h2.y,h2.z,   hf1.x,hf1.y,hf1.z,   hf2.x,hf2.y,hf2.z );
    //if(ia==0) printf( "CPU atom[%i|%i,%i] c %g h1(%g,%g,%g) h2(%g,%g,%g) \n", ia, ing,jng, c, h1.x,h1.y,h1.z,  h2.x,h2.y,h2.z );
    //if(ia==0) printf( "CPU atom[%i|%i,%i] c %g c_ %g E %g fe(%g,%g,%g) fei(%g,%g,%g) \n", ia, ing,jng, c, c_, E, fei.x,fei.y,fei.z,   fapos[ia].x,fapos[ia].y,fapos[ia].z );

    
    //if( ing==0 ) printf( "ngFi[%i,%i,%i](%g,%g,%g) ir %g fang %g c %g c0 %g K %g \n", ia,ing,jng, hf1.x,hf2.y,hf2.z,  ir1, fang, c, c0, K );
    //if( ing==0 ) printf( "neighForce_i[0|%i|%i](%g,%g,%g) invl %g \n", ia, ia*4, hf1.x, hf2.y, hf2.z,  ir1 );
    //if( jng==0 ) printf( "neighForce_j[0|%i|%i](%g,%g,%g) invl %g \n", ia, ia*4, hf2.x, hf2.y, hf2.z,  ir2 );


    return E;
}

inline double evalSigmaPi( int ia, int ing, int ipi, double K ){
    nevalPiSigma++;
    //Vec3d h1,h2;
    //Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();
    Vec3d h1 = apos[ing] - apos[ia];    double ir1  = 1/h1.normalize();
    Vec3d h2 = pipos[ipi];              //double ir2  = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double E = K*c*c;
    double fang = -K*c*2;
    hf1.mul( fang*ir1 );
    hf2.mul( fang     );
    fapos [ing].add( hf1 );
    fpipos[ipi].add( hf2 );
    fapos [ia ].sub( hf1 );
    return E;
}

inline double evalPiPi_T( int ia, int ipi, int jpi, double K ){
    nevalPiPiT++;
    //Vec3d h1,h2;
    //Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();
    Vec3d h1 = pipos[ipi];   //double ir1  = 1/h1.normalize();
    Vec3d h2 = pipos[jpi];   //double ir2  = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double E = K*c*c;
    double fang = -K*c*2;
    hf1.mul( fang );
    hf2.mul( fang );
    fpipos[ipi].add( hf1 );
    fpipos[jpi].add( hf2 );
    //fapos [ia ].sub( hf1+hf2 );  // NOTE! We should not apply recoil here, because pi-vectors are already attached to fulcrum
    return E;
}

inline void orthogonalizePiPi( int ia, int ipi, int jpi ){
    nevalPiPiT++;
    Vec3d hi = pipos[ipi];   // double ir1  = 1/h1.normalize();
    Vec3d hj = pipos[jpi];   //double ir2  = 1/h2.normalize();
    Vec3d ha = hi + hj; ha.normalize(); // ToDo : this can be done more efficiently
    Vec3d hb = hi - hj; hb.normalize();
    pipos[ipi] = ha*0.70710678118 + hb* 0.70710678118;
    pipos[jpi] = ha*0.70710678118 + hb*-0.70710678118;    
}


inline double evalPiPi_I( const Vec2i& at, int ipi, int jpi, double K ){  // interaction between two pi-bonds
    nevalPiPiI++;
    // E = -K*<pi|pj>^2
    //Vec3d d = apos[at.a]-apos[at.b];
    Vec3d h1 = pipos[ipi];   // double ir1  = 1/h1.normalize();
    Vec3d h2 = pipos[jpi];   //double ir2  = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    bool sign = c<0; if(sign) c=-c;
    double E    = -K*c;
    double fang =  K;
    if(sign)fang=-fang;
    hf1.mul( fang );
    hf2.mul( fang );    
    fpipos[ipi].add( hf1 );
    fpipos[jpi].add( hf2 );
    //fapos[ia].sub( hf1+hf2 );
    return E;
}


double eval_bond(int ib){
    //printf( "eval_bond[%i/%i]\n", ib, nbonds );
    Vec2i at = bond2atom[ib];
    //printf( "eval_bond[%i/%i] at(%i,%i)\n", ib, nbonds, at.i, at.j );
    Vec3d f; f.set_sub( apos[at.y], apos[at.x] );
    if(pbcShifts)f.add( pbcShifts[ib] );
    //printf( "bond[%i|%i,%i] (%g,%g,%g) (%g,%g,%g) \n", ib,iat.a,iat.b, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
    
    // int iG=6;
    // if( (at.i==iG)||(at.j==6) ){
    //     printf( "CPU bond[%i,%i] dp(%g,%g,%g) l,k(%g,%g) \n",  at.i,at.j,  f.x,f.y,f.z, bond_l0[ib], bond_k[ib] );
    // }

    double l = f.normalize();
    //printf( "bond[%i|%i,%i] (%g,%g,%g) %g \n", ib,iat.a,iat.b, f.x, f.y, f.z, l );
    //lbond [ib] = l;
    //hbond [ib] = f;
    //if(bondMasked)if(bondMasked[ib])return 0;
    const double k = bond_k[ib];
    double dl = (l-bond_l0[ib]);
    double fr = dl*k*2;
    f.mul( fr );
    //if( isnan( f ) ){ printf("ERROR : eval_bond(%i); f(%g,%g,%g) is NaN fr %g k %g l0 %g \n", ib, f.x,f.y,f.z, fr, k, bond_l0[ib] ); }
    //opengl1renderer.color3f(1.,0.,0.);
    //Draw3D::drawVecInPos(  f, apos[iat.x] );
    //Draw3D::drawVecInPos(  f, apos[iat.y] );
    fapos[at.x].add( f );
    fapos[at.y].sub( f );
    double E = k*dl*dl;
    Eb+=E;

    /*
    if( (at.i==iDEBUG_pick)||(at.j==iDEBUG_pick) ){ 
        int i0=at.i*nneigh_max;
        int j0=at.j*nneigh_max;
        //printf( "eval_bond at(%i,%i) ai[%i,%i,%i,%i] aj[%i,%i,%i,%i] \n", at.i, at.j,  neighs[i0+0],neighs[i0+1],neighs[i0+2],neighs[i0+3],   neighs[j0+0],neighs[j0+1],neighs[j0+2],neighs[j0+3] );
    }
    */
    // --- Pi Bonds

    if(doPiPiI){ // interaction between pi-bonds of given atom
        double Epi = 0;
        if( (at.i<nnode) && (at.j<nnode) ){
            int i0=at.i*nneigh_max;
            int j0=at.j*nneigh_max;
            for(int i=2;i<nneigh_max;i++){
                int ipi   = neighs[i0+i];
                //double Ki = Kneighs[i0+i] * Kpipi;
                if(ipi>=0) continue;
                for(int j=2;j<nneigh_max;j++){
                    int jpi=neighs[j0+j];
                    if(jpi>=0) continue;
                    //if( (at.i==iDEBUG_pick)||(at.j==iDEBUG_pick) ){  printf(" i,j[%i,%i]->ipi,jpi[%i,%i,] \n",  i,j, ipi,jpi ); };
                    //double Kij = Ki*Kneighs[j0+j];
                    double Kij = Kpipi;
                    //printf( "bond[%i,] evalPiPi_I(%i,%i) Kij %g \n", ib, i,j, Kij );
                    Epi += evalPiPi_I( at, -ipi-2,-jpi-2, Kij );
                }
            }
        }
        EppI+=Epi;
        E   +=Epi;
    }

    return E;
}


double eval_bond_neigh(int ib, Vec3d h, double l){
    double k  = bond_k[ib];
    double dl = (l-bond_l0[ib]);
    double fr = dl*k*2;
    h.mul( fr );
    Vec2i B = bond2atom[ib];
    fapos[B.x].add( h );
    fapos[B.y].sub( h );
    double E = k*dl*dl;
    Eb+=E;

    double Kij   = bond_kPi[ib];
    const bool   bPiPi = Kij >  0.001;
    const bool   bPiEp = (Kij < -0.001)&&doEpi;
    if(bPiEp)Kij=-Kij;

    if(doPiPiI && (bPiPi||bPiEp) ){ // interaction between pi-bonds of given atom
        double Epi = 0;
        if( (B.i<nnode) && (B.j<nnode) ){
            const int i0=B.i*nneigh_max;
            const int j0=B.j*nneigh_max;
            for(int i=2;i<nneigh_max;i++){
                const int ing   = neighs[i0+i];
                const bool  bi = ing<0;
                Vec3d hi; double il1;
                if(bi){ hi=pipos[-ing-2];            il1=1;              }
                else  { hi=apos [ ing  ]-apos[B.i];  il1=hi.normalize(); }
                for(int j=2;j<nneigh_max;j++){
                    const int jng=neighs[j0+j];
                    const bool  bj = jng<0;
                    
                    Vec3d hj; double il2;
                    if(bj){ hj=pipos[-jng-2];            il2=1;              }
                    else  { hj=apos [ jng  ]-apos[B.j];  il2=hi.normalize(); }
                    Vec3d f1,f2;
                    Epi += evalPiAling( hi, hj, il1, il2, Kij, f1, f2 );  
                    if(bi&&bj){ // Pi-Pi
                        //Epi += evalPiPi_I( at, -ipi-2,-jpi-2, Kij );
                        fpipos[-ing-2].add( f1 );
                        fpipos[-jng-2].add( f2 );
                    }
                    else if( bPiEp ){
                        if      ( bi && (jng>=ie0) ){   
                            fpipos[-ing-2].add( f1 );
                            fapos [ jng  ].add( f2 );
                            fapos [ B.j  ].sub( f2 );
                        }else if( bj && (ing>=ie0) ){   
                            fpipos[-jng-2].add( f2 );
                            fapos [ ing  ].add( f1 );
                            fapos [ B.i  ].sub( f1 );
                        }
                    }
                }
            }
        }
        EppI+=Epi;
        E   +=Epi;
    }
    return E;
}

double eval_neighs(int ia){
    //printf( "eval_neighs(%i) \n", ia ); printNeigh(ia);
    double E=0;
    int ioff = ia*nneigh_max;
    //bool bDEBUG=(ia==7);
    //if(bDEBUG)printf( "#-atom[%i] [%i,%i,%i,%i] \n", ia, neighs[ioff+0],neighs[ioff+1],neighs[ioff+2],neighs[ioff+3] );
    const bool doPiPiT_   = doPiPiT; 
    const bool doAngles_  = doAngles;
    const bool doPiSigma_ = doPiSigma;
    Quat4d params = NeighParams[ia]; 
    double c0  = params.x;
    double Kss = params.y;
    double Ksp = params.z;
    for(int i=0; i<nneigh_max; i++){
        int    ing = neighs[ioff+i];
        //double Ki  = Kneighs[ioff+i];
        bool bipi=ing<0;
        //if( ing<0 ){ break; // empty
        //}else if( ing<ipi0 ){ // signa-bonds
        for(int j=i+1; j<nneigh_max; j++){
            int jng  = neighs[ioff+j];
            bool bjpi=jng<0;
            //double K = Kneighs[ioff+j]*Ki;
            //printf( "eval_neighs[%i,%i] ing,jng(%i,%i) \n", i,j, ing,jng );
            if( bipi ){
                if(bjpi){ if(doPiPiT_){
                    double e=evalPiPi_T   ( ia, -ing-2, -jng-2, Kpipi ); EppT+=e; E+=e;
                    //orthogonalizePiPi( ia, -ing-2, -jng-2 );
                }   }
                else    { if(doPiSigma_){ double e=evalSigmaPi( ia, jng, -ing-2, Ksp ); Eps+=e; E+=e; } }
            }else{
                if(bjpi){ if(doPiSigma_){ double e=evalSigmaPi( ia, ing, -jng-2, Ksp ); Eps+=e; E+=e; } } 
                else    {  if(doAngles_){
                    //E+= evalSigmaSigma_dist( ia, ing, jng, K );    Eps+=e; E+=e; 
                    double e= evalSigmaSigma_cos( ia, ing, jng, Kss, c0 );    Ea+=e; E+=e; 
                }  }
            }
        }
    }
    return E;
}

void wrapBondVec( Vec3d& d ){
    Vec3d u;
    invLvec.dot_to( d, u );
    //if(hcell.a>0.5){ h.sub( lvec.a ); }else{ h.add( lvec.a ); };
    //printf( "before wrap u=(%g,%g,%g) d=(%g,%g,%g) \n", u.x,u.y,u.z, d.x,d.y,d.z );
    u.a=u.a+(1-(int)(u.a+1.5));
    u.b=u.b+(1-(int)(u.b+1.5));
    u.c=u.c+(1-(int)(u.c+1.5));
    lvec.dot_to_T( u, d );
    //printf( "after  wrap u=(%g,%g,%g) d=(%g,%g,%g) \n", u.x,u.y,u.z, d.x,d.y,d.z );
    //lvecT.dot_to( u, d );
}

double eval_neighs_new(int ia){
    //printf( "MMFFsp3::eval_neighs_new(%i) \n", ia );
    double E=0;
    int ioff = ia*nneigh_max;
    const bool doPiPiT_   = doPiPiT; 
    const bool doAngles_  = doAngles;
    const bool doPiSigma_ = doPiSigma;
    Quat4d params = NeighParams[ia]; 
    double c0  = params.x;
    double Kss = params.y;
    double Ksp = params.z;

    const Vec3d pa = apos[ia]; 

    Vec3d    hs[4];
    double  ils[4];
    
    // ToDo: there could be computed also the bonds (similarly as done in OpenCL)
    //if(idebug>0)printf( "DEBUG eval_neighs_new() [%i] neighs[%i,%i,%i,%i] abonds[%i,%i,%i,%i] \n", ia,  neighs[ioff+0],neighs[ioff+1],neighs[ioff+2],neighs[ioff+3],     abonds[ ioff+0],abonds[ ioff+1],abonds[ ioff+2],abonds[ ioff+3] );
    for(int i=0; i<nneigh_max; i++){
        int  ing = neighs[ioff+i];
        if(ing==-1)continue;
        if(ing<0){
            ils[i] = 1;
            hs [i] = pipos[-ing-2];;
        }else{
            int  ib  = abonds[ ioff+i];
            Vec3d h; h.set_sub( apos[ing], pa );
            double c=1;
            
            if(bPBCbyLvec){  
                wrapBondVec( h );
                //opengl1renderer.color3f(1.,0.,1.); Draw3D::drawVecInPos( h, pa );
            }else if(pbcShifts){
                if( bond2atom[ib].a!=ia ){ c=-1; }; // bond should be inverted
                h.add_mul( pbcShifts[ib], c );
            }
            
            double l = h.normalize();
            ils[i] = 1/l;
            hs [i] = h;
            //printf( "eval_h ia(%i)[%i] h(%g,%g,%g|l=%g)\n", ia, i, h.x,h.y,h.z,  ils[i] );
            if(ia<=ing){
                //doPiPiI = false; // WARRNING BIG DEBUG
                double  e= eval_bond_neigh(ib, h*c, l);
                //if( isnan(e) ){ printf( "DEBUG atom[%i] neigh[%i|ib=%i,ja=%i] Eb=%g h(%g,%g,%g) l=%g \n", ia, i, ib, ing, e, h.x,h.y,h.z, l );}
                Eb+=e;
                //if(idebug>0)printf( "DEBUG bond[%i|%i,%i] l=%g \n", ib, ia, ing, l );
            }
        }
        // ToDd: Compute bonds here
    }
    Vec3d f1,f2;
    for(int i=0; i<nneigh_max; i++){
        int    ing = neighs[ioff+i];
        if(ing==-1)continue;
        bool bipi=ing<0;
        for(int j=i+1; j<nneigh_max; j++){
            int jng  = neighs[ioff+j];
            if(jng==-1)continue;
            bool bjpi=jng<0;
            //printf( "DEBUG angle[%i,%i] bpi(%i,%i)\n", i, j, bipi, bjpi );
            if( bipi ){
                if(bjpi){ 
                    if(doPiPiT_){ // Pi and Pi should be orthogonal
                        double e = evalAngleCos( hs[i], hs[j], ils[i], ils[j], Kpipi, 0, f1, f2 );    
                        fpipos[-ing-2].add( f1     );
                        fpipos[-jng-2].add( f2     );
                        //double e=evalPiPi_T   ( ia, -ing-2, -jng-2, Kpipi ); EppT+=e; E+=e;
                        //printf( "eval_PPT ia(%i)[%i,%i] f1(%g,%g,%g)f2(%g,%g,%g)\n", ia, i,j, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z );
                        EppT+=e; E+=e;
                    }   
                }else{ 
                    if(doPiSigma_){  // Pi and Sigma should be orthogonal
                        double e = evalAngleCos( hs[i], hs[j], 1, ils[j], Ksp, 0, f1, f2 ); 
                        fpipos[-ing-2].add( f1 );
                        fapos [ jng  ].add( f2 );
                        fapos [ ia   ].sub( f2 );
                        //double e=evalSigmaPi( ia, jng, -ing-2, Ksp ); Eps+=e; E+=e; } 
                        //printf( "eval_PS ia(%i)[%i,%i] f1(%g,%g,%g)f2(%g,%g,%g)\n", ia, i,j, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z );
                        Eps+=e; E+=e;
                    }
                }
            }else{
                if(bjpi){ 
                    if(doPiSigma_){ // Pi and Sigma should be orthogonal
                        double e = evalAngleCos( hs[i], hs[j], ils[i], 1, Ksp, 0, f1, f2 );  
                        fapos [ ing  ].add( f1 );
                        fpipos[-jng-2].add( f2 );
                        fapos [ ia   ].sub( f1 );
                        //double e=evalSigmaPi( ia, ing, -jng-2, Ksp ); Eps+=e; E+=e; 
                        //printf( "eval_SP ia(%i)[%i,%i] f1(%g,%g,%g)f2(%g,%g,%g)\n", ia, i,j, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z );
                        Eps+=e; E+=e; 
                    } 
                } else {  
                    if(doAngles_){ // Angle between two sigma bonds
                        double c0_;
                        if( (ing>=ie0)!=(jng>=ie0) ){  c0_=params.w; }else{ c0_=c0; } // special angle for electron pairs
                        if( (ing>=ie0)&&(jng>=ie0) ){  c0_=-1;       }                      // special angle for electron pairs
                        //if(idebug)printf( "atom[%i]ss[%i,%i] c0,K(%g,%g) \n", ia,ing,jng, c0_, Kss  );
                        //printf( "eval_SS ia(%i)[%i,%i] f1(%g,%g,%g)f2(%g,%g,%g)\n", ia, i,j, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z );
                        double e = evalAngleCos( hs[i], hs[j], ils[i], ils[j], Kss, c0_, f1, f2 );  
                        if(bSubtractAngleNonBond){
                            Vec3d fij=Vec3dZero;
                            Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                            e += addAtomicForceLJQ( apos[ing]-apos[jng], fij, REQij );
                            f1.add(fij);
                            f2.sub(fij);
                        }
                        fapos[ing].add( f1     );
                        fapos[jng].add( f2     );
                        fapos[ia ].sub( f1+f2 );
                        //double e= evalSigmaSigma_cos( ia, ing, jng, Kss, c0 );    Ea+=e; E+=e; 
                        Ea+=e; E+=e;
                    }  
                }
            }
        }
    }
    //printf( "DEBUG atom[%i] E=%g \n", ia, E );
    //printf( "MMFFsp3::eval_neighs_new(%i) DONE \n", ia );
    return E;
}

//inline double evalAngle( const Vec3d& h1, const Vec3d& h2, double ir1, double ir2, double K, double c0, Vec3d& f1, Vec3d& f2 );
   

double eval_bonds(){
    double E=0;
    for(int ib=0; ib<nbonds; ib++){ E+=eval_bond(ib); }
    return E;
}

double eval_neighs(){
    double E=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_neighs(ia); }
    return E;
}

double eval_neighs_new(){
    double E=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_neighs_new(ia); }
    return E;
}

double eval( bool bClean=true, bool bCheck=true ){
    //printf( "DEBUG MMFFsp3.eval() 1 \n" );
    //printf( "MMFFsp3::eval( nnode %i natoms %i nvecs %i ) \n", nnode, natoms, nvecs );
    if(bClean){ cleanAll(); }
    normalizePi(); 
    if(bCheck)ckeckNaN_d(npi, 3, (double*)pipos, "pipos" );
    //if(doBonds )eval_bonds();   if( isnan( Eb) ){ printf("ERROR : Eb = eval_bonds();  is NaN  \n"); checkNaNs(); exit(0); }
    //if(doNeighs)eval_neighs();  if( isnan( Ea) ){ printf("ERROR : Ea = eval_neighs(); is NaN  \n"); checkNaNs(); exit(0); }
    eval_neighs_new();
    Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}

double eval_check(){
    if(verbosity>0){
        printf(" ============ check MMFFsp3 START\n " );
        printSizes();
    }
    eval();
    checkNans();
    if(verbosity>0)printf(" ============ check MMFFsp3 DONE\n " );
    return Etot;
} 

void evalPi0s(){
    for(int ia=0; ia<nnode; ia++){
        int ioff = ia*nneigh_max;
        int* ngs = neighs+ioff;
        int ipi  = ngs[3];
        int nbond=0;
        if(ipi<0){
            Vec3d pi0=Vec3dZero;
            for(int j=0; j<3; j++){
                int ja  = ngs[j];
                if( (ja>=0) && (ja<nnode) ){
                    int jpi = neighs[ja*nneigh_max+3];
                    if(jpi<0){
                        pi0.add( pipos[-jpi-2] );
                        nbond++;
                    }
                }
            }
            //pi0.mul(1./nbond);
            pi0.normalize();
            pi0s[-ipi-2].f=(Vec3f)pi0;
        }
    }
}

void chargeToEpairs( Quat4d* REQs, double cQ=-0.2, int etyp=-1 ){
    //printf( "chargeToEpairs() \n" );
    for( int ib=0; ib<nbonds; ib++ ){
        Vec2i b = bond2atom[ib];
        if( atype[b.i]==etyp ){ REQs[b.i].z+=cQ; REQs[b.j].z-=cQ; }
        if( atype[b.j]==etyp ){ REQs[b.j].z+=cQ; REQs[b.i].z-=cQ;  }
    }
}

void makePiNeighs( int* pi_neighs ){
    printf( "makePiNeighs() \n" );
    for(int ia=0; ia<nnode; ia++){
        int* ngs = neighs+(ia*nneigh_max);
        int ipi  = ngs[3];
        if(ipi<0){
            int ioff = (-ipi-2)*4;
            int nbond=0;
            printf( "ngs[%i](%i,%i,%i) ->", ia, ngs[0],ngs[1],ngs[2] );
            for(int j=0; j<3; j++){
                int ja  = ngs[j];
                if( (ja>=0) && (ja<nnode) ){
                    int jpi = neighs[ja*nneigh_max+3];
                    printf( " [%i](%i,%i,%i,%i) ", ja, neighs[ja*nneigh_max+0],neighs[ja*nneigh_max+1],neighs[ja*nneigh_max+2],neighs[ja*nneigh_max+3] ); 
                    if(jpi<0){
                        pi_neighs[ioff+nbond]=natoms-jpi-2;
                        nbond++;
                    }
                }
            }
            printf( "nbond=%i \n", nbond); 
        }
    }
    for(int i=0; i<npi; i++){ printf("pi[%i] (%i,%i,%i,%i) \n", i, pi_neighs[i*4+0],pi_neighs[i*4+1],pi_neighs[i*4+2],pi_neighs[i*4+3] ); };
    //exit(0);
}


//double getFmax(){ double Fmax=0; for(int i=0; i<natoms;i++){ _max( Fmax, aforce[i].norm2() ); }; return Fmax; };
//void moveGD(double dt){ for(int i=0; i<natoms;i++){ apos[i].add_mul( aforce[i],dt); } };

inline bool insertNeigh( int i, int j ){
    int off = i*nneigh_max;
    int*   ngs   = neighs + off;
    //double* Kngs = Kneighs + off;
    int ii;
    for(ii=0; ii<nneigh_max; ii++){
        if(ngs[ii]<0){ ngs[ii]=j; break; }
    }
    return ii<nneigh_max;
}

void bond2neighs(){
    int nn=natoms*nneigh_max;
    for(int i=0; i<nn; i++){ neighs[i]=-1; };
    for(int i=0; i<nbonds; i++){
        Vec2i b = bond2atom[i];
        if(b.i<natoms)insertNeigh( b.i, b.j );
        if(b.j<natoms)insertNeigh( b.j, b.i );
    }
    int ipi=0;
    for(int i=0; i<nn; i++){ if(neighs[i]<0){ neighs[i]=ipi; ipi++; }; };
}

Vec3d evalPiTorq(){
    Vec3d tq; tq.set(0.0);
    for(int ia=0; ia<nnode; ia++){
        int ioff = ia*nneigh_max;
        for(int i=0; i<nneigh_max; i++){
            int ing = neighs[ioff+i];
            if(ing<0){
                int ipi=-ing-2;
                tq.add_cross( apos[ia]+pipos[ipi], fpipos[ipi] );
                //tq.add_cross( apos[ia], fpipos[ipi] );
            };
        }
    }
    return tq;
}

int selectCaps(int n, int* sel, int* selout){
    // this function finds cap-atoms (including e-pairs) attached to node atoms
    int nout=0;
    for(int i=0;i<n; i++){
        int ia = sel[i];
        if(ia>=nnode)continue;
        int* ngs=neighs + ia*nneigh_max; 
        for(int j=0;j<nneigh_max;j++){
            int ja = ngs[j];
            if(ja>=0){ nout=ja; nout++; } // ignore pi-bonds
        }
    }
    return nout;
}

void rotateNodes(int n, int* sel, Vec3d p0, Vec3d ax, double phi ){
    ax.normalize();
    double ca=cos(phi);
    double sa=sin(phi);
    for(int i=0;i<n; i++){
        int ia = sel[i];
        if(ia>=nnode)continue;
        apos[ia].rotate_csa( ca, sa, ax, p0 );
        int* ngs=neighs + ia*nneigh_max; 
        for(int j=0;j<nneigh_max;j++){
            int ja = ngs[j];
            if(ja>=0){ if(ja>nnode) apos[ ja  ].rotate_csa( ca, sa, ax, p0 ); } // cap atoms
            else     {             pipos[-ja-1].rotate_csa( ca, sa, ax     ); } // pi-bonds
        }
    }
}

bool checkNans( bool bExit=true, bool bPi=true, bool bA=true ){
    bool ret = false;
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*)  apos,   "apos" );
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*) fapos,  "fapos" );
    if(bPi) ret |= ckeckNaN_d(npi,     3, (double*) pipos,  "pipos" );
    if(bPi) ret |= ckeckNaN_d(npi,     3, (double*)fpipos, "fpipos" );
    if(bExit&&ret){ printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); exit(0); };
    return ret;
}

void printSizes(){ printf( "MMFFsp3::printSizes(): nDOFs(%i) natoms(%i) nnode(%i) ncap(%i) npi(%i) nbonds(%i) nvecs(%i) \n", nDOFs,natoms,nnode,ncap,npi,nbonds,nvecs ); };

void printAtom (int i){ printf( "Atom[%i] pos(%g,%g,%g) \n", i, apos[i].x,apos[i].y,apos[i].z ); }
void printParam(int i){ printf( "Atom[%i] par(%g,%g,%g,%g) \n", i, NeighParams[i].x, NeighParams[i].y, NeighParams[i].z, NeighParams[i].w ); }
void printBond (int i){ printf( "bond[%i|%i,%i] l0 %g k %g kPi %g \n", i, bond2atom[i].i,bond2atom[i].j, bond_l0[i], bond_k[i], bond_kPi[i] ); }

void printAtoms     (){     printf( "MMFFsp3::printAtoms() : \n" ); for(int i=0;i<natoms;i++){ printAtom(i); }}
void printAtomParams(){     printf( "MMFFsp3::printParams() : \n" ); for(int i=0;i<nnode;i++){ printParam(i); }}
void printBonds     (){     printf( "MMFFsp3::printBonds() : \n" ); for(int i=0;i<nbonds;i++){ printBond(i); }}


void printAtomPis(){ 
    printf( "MMFFsp3::printAtomPis() : \n" );
    for(int ia=0; ia<nnode; ia++ ){
        int* ngs = neighs + ia*nneigh_max;
        for(int j=0; j<nneigh_max; j++){
            int ing = ngs[j];
            if(ing<0){ int ipi=-ing-2; printf("pi[%i]atom[%i]ng[%i] pos(%g,%g,%g) force(%g,%g,%g)\n", ipi, ia, j, pipos[ipi].x,pipos[ipi].y,pipos[ipi].z,   fpipos[ipi].x,fpipos[ipi].y,fpipos[ipi].z );  };
        }
    }

}

bool checkBonds( double factor=1.5, bool bPrintPBC=false ){
    bool bErr=false;
    for(int ib=0;ib<nbonds;ib++){
        Vec2i b = bond2atom[ib];
        Vec3d d = apos[b.j] - apos[b.i]; 
        if( pbcShifts ){ d.add(pbcShifts[ib]); }
        double r2 = d.norm();
        double R  = bond_l0[ib]*factor;  
        if( bPrintPBC ){ if(pbcShifts[ib].norm2()>0.1) printf( "PBC bond[%i(%i,%i)] pbcShifts(%g,%g,%g) \n", ib, b.i,b.j, pbcShifts[ib].x,pbcShifts[ib].y,pbcShifts[ib].z ); };
        if( r2>(R*R)  ){ printf( "ERROR bond[%i(%i,%i)] r=%g > Lmax=%g pbcShifts(%g,%g,%g) \n", ib, b.i,b.j, sqrt(r2), R, pbcShifts[ib].x,pbcShifts[ib].y,pbcShifts[ib].z ); bErr=true; }
    }
    return bErr;
}

void printNeigh(int ia){
    printf( "atom[%i] neighs{", ia );
    for(int j=0;j<nneigh_max;j++){
        int ij=ia*nneigh_max + j;
        printf( "%i,", neighs[ij] );
    }
    /*
    printf( "] K(" );
    for(int j=0;j<nneigh_max;j++){
        int ij=ia*nneigh_max + j;
        printf( "%g,", Kneighs[ij] );
    }
    */
    printf( "}\n");
}
void printNeighs(){
    printf( "MMFFsp3::printNeighs(): nnode=%i \n", nnode );
    for(int i=0;i<nnode;i++){ printNeigh(i); }
}

void printPis(){
    printf( "MMFFsp3::printPis(): npi=%i \n", npi );
    for(int i=0;i<npi;i++){ printf( "pipos[%i] (%g,%g,%g)\n", i, pipos[i].x,pipos[i].x,pipos[i].x ); }
}


double optimalTimeStep(double m=1.0){
    double Kmax = 1.0;
    for(int i=0; i<nbonds; i++){ Kmax=fmax(Kmax, bond_k[i] ); }
    return M_PI*2.0*sqrt(m/Kmax)/15.0;  // dt=T/10;   T = 2*pi/omega = 2*pi*sqrt(m/k)
}


};


#endif
