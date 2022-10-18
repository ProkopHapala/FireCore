
#ifndef MMFFsp3_h
#define MMFFsp3_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "molecular_utils.h"

// ======================
// ====   MMFFsp3
// ======================


bool isnan(Vec3d& v){ return (isnan(v.x)||isnan(v.y)||isnan(v.z)); }

void ckeckNaN_d(int n, int m, double* xs, const char* pre ){
    for(int i=0; i<n;i++){
        bool b=false;
        for(int j=0; j<m;j++){
            int ij=i*m+j;
            b|=isnan( xs[ij] );
        }
        if(b){
            printf("%s[%i](", pre, i );
            for(int j=0; j<m;j++){
                int ij=i*m+j;
                printf("%g,", xs[ij] );
            }
            printf(")\n");
        }
    }
}


class MMFFsp3{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0,nnode=0,ncap=0,npi=0,nbonds=0;
    bool bPBC=false;
    double Etot,Eb,Ea, Eps,EppT,EppI;

    double * DOFs  = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    double Kpipi = 0.25;
    int  ipi0=0;
    int   * atype=0;
    Vec3d * apos=0;
    Vec3d * fapos=0;
    Vec3d * pipos=0;
    Vec3d * fpipos=0;
    //Vec3d * cappos=0;
    //Vec3d * fcappos=0;

    // --- Parameters
    Vec2i  * bond2atom = 0;
    double * bond_l0   = 0;  // [A]
    double * bond_k    = 0;  // [eV/A] ?
    Vec3d  * pbcShifts = 0;  // [A]

    int*    aneighs = 0;  // [natom*nneigh_max] neigbors of atom
    double* Kneighs = 0;  // [natom*nneigh_max] angular stiffness for given neighbor

    // NOTE - For the moment we do not use temporary bonds - to keep things simple
    // --- Axuliary variables
    //double * lbond  = 0;   // bond lengths
    //Vec3d  * hbond  = 0;   // normalized bond unitary vectors

    bool doPiPiI  =true;
    bool doPiPiT  =true;
    bool doPiSigma=true;
    bool doAngles =true;
    int nevalPiPiI  =0;
    int nevalPiPiT  =0;
    int nevalPiSigma=0;
    int nevalAngles =0;

    bool bDEBUG_plot=0;
    int iDEBUG_n=0,iDEBUG_pick=0;

void realloc( int nnode_, int nbonds_, int npi_, int ncap_, bool bNeighs=true ){
    nnode=nnode_; ncap=ncap_; nbonds=nbonds_; npi=npi_; 
    natoms=nnode_+ncap_; 
    nDOFs = (natoms + npi)*3;
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
    //_realloc( lbond     , nbonds );
    //_realloc( hbond     , nbonds );
    //printf( "MMFFsp3::realloc() 3 \n" );
    if(bNeighs){
        int nn = nnode*nneigh_max;
        _realloc( aneighs, nn );
        _realloc( Kneighs, nn );
    }
    //printf( "MMFFsp3::realloc() DONE \n" );
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


int i_DEBUG = 0;

inline double evalSigmaSigma_dist( int ia, int ing, int jng, double K ){
    nevalAngles++;
    //K = -1.0;
    //  E = -K*|pi-pj|
    // ToDo: This may be made more efficient if we store hbonds
    Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();  // j->i
    Vec3d hi = apos[ing] - apos[ia];    double ri  = hi.normalize();
    Vec3d hj = apos[jng] - apos[ia];    double rj  = hj.normalize();
    Vec3d fi  = d* K;
    Vec3d fj  = d*-K;
    Vec3d fiI = hi* hi.dot(fi);   
    Vec3d fjI = hj* hj.dot(fj);   
    //printf(  "ang[%i|%i,%i] c1 %g c2 %g \n", ia, ing, jng, hi.dot(fi), hj.dot(fj)  );
    //glColor3f(1.,0.,0.);
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

inline double evalSigmaSigma_cos(  int ia, int ing, int jng, double K ){
    nevalAngles++;
    //Vec3d h1,h2;
    //Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();
    Vec3d h1 = apos[ing] - apos[ia];    double ir1 = 1/h1.normalize();
    Vec3d h2 = apos[jng] - apos[ia];    double ir2 = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double c_ = c+1.0;
    double E    =  K*c_*c_;
    double fang = -K*c_*2;
    hf1.mul( fang*ir1 );
    hf2.mul( fang*ir2 );
    //if(bDEBUG_plot){
    // //if(ia==7){
    //     //printf( "evalAngle_cos [%i|%i,%i,%i] c(%g,%g) \n", iDEBUG_pick, ia, ing, jng, hf1.dot(h1), hf1.dot(h1)  );
    //     //printf( "evalAngle_cos [%i|%i,%i,%i] c %g c_ %g fang_ %g E %g \n", iDEBUG_pick, ia, ing, jng, c, c_, fang, E );
    //     //printf( "evalAngle_cos bDEBUG_plot [%i|%i,%i] iDEBUG_pick %i \n", ia, ing, jng, iDEBUG_pick );
    //     glColor3f(1.,0.,0.);
    //     Draw3D::drawVecInPos(  hf1, apos[ing] );
    //     Draw3D::drawVecInPos(  hf2, apos[jng] );
    //     Draw3D::drawVecInPos(  (hf1+hf2)*-1.0, apos[ia] );
    // }
    fapos[ing].add( hf1     );
    fapos[jng].add( hf2     );
    fapos[ia ].sub( hf1+hf2 );
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

    // if(bDEBUG_plot)
    // if(ipi==iDEBUG_pick)
    // {
    //     float fsc=100.0;
    //      //printf( "evalAngle_cos [%i|%i,%i,%i] c(%g,%g) \n", iDEBUG_pick, ia, ing, jng, hf1.dot(h1), hf1.dot(h1)  );
    //      //printf( "evalAngle_cos [%i|%i,%i,%i] c %g c_ %g fang_ %g E %g \n", iDEBUG_pick, ia, ing, jng, c, c_, fang, E );
    //      //printf( "evalAngle_cos bDEBUG_plot [%i|%i,%i] iDEBUG_pick %i \n", ia, ing, jng, iDEBUG_pick );
    //      glColor3f(1.,0.,0.);
    //      Draw3D::drawVecInPos(  hf1*fsc, apos [ing] );
    //      Draw3D::drawVecInPos(  hf2*fsc, pipos[ipi]+apos[ia] );
    //      Draw3D::drawVecInPos(  (hf1+hf2)*-fsc, apos[ia] );
    // }

    fapos [ing].add( hf1 );
    fpipos[ipi].add( hf2 );
    //fapos [ia ].sub( hf1+hf2 ); // NOTE: we should not apply recoil of pi-orbital here because it is in coordinates relative to  fulcrum
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
    //if( (ia==iDEBUG_pick) ){
    //    printf( "evalPiPi_I ia(%i) ipi,jpi(%i,%i) fang %g c %g \n", ia,  ipi,jpi, fang, c );
    //    glColor3f(1.,0.,0.);
    //    Draw3D::drawVecInPos(  hf1*100.0, apos[ia]+pipos[ipi] );
    //    Draw3D::drawVecInPos(  hf2*100.0, apos[ia]+pipos[jpi] );
    //}
    fpipos[ipi].add( hf1 );
    fpipos[jpi].add( hf2 );
    //fapos [ia ].sub( hf1+hf2 );  // NOTE! We should not apply recoil here, because pi-vectors are already attached to fulcrum
    return E;
}

inline void orthogonalizePiPi( int ia, int ipi, int jpi ){
    nevalPiPiT++;
    Vec3d hi = pipos[ipi];   // double ir1  = 1/h1.normalize();
    Vec3d hj = pipos[jpi];   //double ir2  = 1/h2.normalize();
    //double c = h1.dot(h2);
    //double c = hi.dot( hj );
    Vec3d ha = hi + hj; ha.normalize(); // ToDo : this can be done more efficiently
    Vec3d hb = hi - hj; hb.normalize();
    //double ca = sqrt( 0.5*(1+c) );
    //double sa = sqrt( 0.5*(1-c) );
    //glColor3f(0.,0.,1.);
    //Draw3D::drawVecInPos(  ha, apos[ia]);
    //Draw3D::drawVecInPos(  hb, apos[ia] );
    pipos[ipi] = ha*0.70710678118 + hb* 0.70710678118;
    pipos[jpi] = ha*0.70710678118 + hb*-0.70710678118;
    //glColor3f(0.,0.7,1.);
    //Draw3D::drawVecInPos(  pipos[ipi], apos[ia] );
    //Draw3D::drawVecInPos(  pipos[jpi], apos[ia] );
    
};


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
    //double c2  = c*c;
    //double c2_ = c2-1; 
    //double E    =  K*c2_;
    //double fang =  K*c*4;
    //double E    =  K*c2_*c2_;
    //double fang = -K*c2_*c*4;
    hf1.mul( fang );
    hf2.mul( fang );
    //if( isnan( hf1 )||isnan( hf1 ) ){ printf("ERROR : evalPiPi_I(%i,%i|%i,%i); hf1(%g,%g,%g) hf2(%g,%g,%g) is NaN c %g k %g \n", at.i,at.j, ipi,jpi,  hf1.x,hf1.y,hf1.z,  hf2.x,hf2.y,hf2.z,  c, K ); }
    
    //if( (at.i==iDEBUG_pick)||(at.j==iDEBUG_pick) ){
    //    //printf( "evalPiPi_I at(%i,%i) ipi,jpi(%i,%i) fang %g c %g \n", at.i, at.j,  ipi,jpi, fang, c );
    //    glColor3f(1.,0.,0.);
    //    Draw3D::drawVecInPos(  hf1*10.0, apos[at.i]+pipos[ipi] );
    //    Draw3D::drawVecInPos(  hf2*10.0, apos[at.j]+pipos[jpi] );
    //}
    
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
    //if(pbcShifts)f.add( pbcShifts[ib] );
    //printf( "bond[%i|%i,%i] (%g,%g,%g) (%g,%g,%g) \n", ib,iat.a,iat.b, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
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
    //glColor3f(1.,0.,0.);
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
        //printf( "eval_bond at(%i,%i) ai[%i,%i,%i,%i] aj[%i,%i,%i,%i] \n", at.i, at.j,  aneighs[i0+0],aneighs[i0+1],aneighs[i0+2],aneighs[i0+3],   aneighs[j0+0],aneighs[j0+1],aneighs[j0+2],aneighs[j0+3] );
    }
    */
    // --- Pi Bonds

    if(doPiPiI){ // interaction between pi-bonds of given atom
        double Epi = 0;
        if( (at.i<nnode) && (at.j<nnode) ){
            int i0=at.i*nneigh_max;
            int j0=at.j*nneigh_max;
            for(int i=2;i<nneigh_max;i++){
                int ipi   = aneighs[i0+i];
                double Ki = Kneighs[i0+i] * Kpipi;
                if(ipi>=0) continue;
                for(int j=2;j<nneigh_max;j++){
                    int jpi=aneighs[j0+j];
                    if(jpi>=0) continue;
                    //if( (at.i==iDEBUG_pick)||(at.j==iDEBUG_pick) ){  printf(" i,j[%i,%i]->ipi,jpi[%i,%i,] \n",  i,j, ipi,jpi ); };
                    double Kij = Ki*Kneighs[j0+j];
                    //printf( "bond[%i,] evalPiPi_I(%i,%i) Kij %g \n", ib, i,j, Kij );
                    Epi += evalPiPi_I( at, -ipi-1,-jpi-1, Kij );
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
    //if(bDEBUG)printf( "#-atom[%i] [%i,%i,%i,%i] \n", ia, aneighs[ioff+0],aneighs[ioff+1],aneighs[ioff+2],aneighs[ioff+3] );
    const bool doPiPiT_   = doPiPiT; 
    const bool doAngles_  = doAngles;
    const bool doPiSigma_ = doPiSigma;
    for(int i=0; i<nneigh_max; i++){
        int    ing = aneighs[ioff+i];
        double Ki  = Kneighs[ioff+i];
        bool bipi=ing<0;
        //if( ing<0 ){ break; // empty
        //}else if( ing<ipi0 ){ // signa-bonds
        for(int j=i+1; j<nneigh_max; j++){
            int jng  = aneighs[ioff+j];
            bool bjpi=jng<0;
            double K = Kneighs[ioff+j]*Ki;
            //printf( "eval_neighs[%i,%i] ing,jng(%i,%i) \n", i,j, ing,jng );
            if( bipi ){
                if(bjpi){ if(doPiPiT_){
                    double e=evalPiPi_T   ( ia, -ing-1, -jng-1, K ); EppT+=e; E+=e;
                    //orthogonalizePiPi( ia, -ing-1, -jng-1 );
                }   }
                else    { if(doPiSigma_){ double e=evalSigmaPi( ia, jng, -ing-1, K ); Eps+=e; E+=e; } }
            }else{
                if(bjpi){ if(doPiSigma_){ double e=evalSigmaPi( ia, ing, -jng-1, K ); Eps+=e; E+=e; } } 
                else    {  if(doAngles_){
                    //E+= evalSigmaSigma_dist( ia, ing, jng, K );    Eps+=e; E+=e; 
                    double e= evalSigmaSigma_cos( ia, ing, jng, K );    Ea+=e; E+=e; 
                }  }
            }
        }
    }
    return E;
}

double eval_bonds(){
    double E=0;
    for(int ib=0; ib<nbonds; ib++){ E+=eval_bond(ib); }
    return E;
}

double eval_neighs(){
    double E=0;
    iDEBUG_n=0;
    for(int ia=0; ia<nnode; ia++){ E+=eval_neighs(ia); }
    return E;
}

double eval( bool bClean=true ){
    //printf( "DEBUG MMFFsp3.eval() 1 \n" );
    if(bClean){ cleanAll(); }
    normalizePi(); 
    ckeckNaN_d(npi, 3, (double*)pipos, "pipos" );
    eval_bonds();   if( isnan( Eb) ){ printf("ERROR : Eb = eval_bonds();  is NaN  \n"); checkNaNs(); exit(0); }
    eval_neighs();  if( isnan( Ea) ){ printf("ERROR : Ea = eval_neighs(); is NaN  \n"); checkNaNs(); exit(0); }
    Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
};

//double getFmax(){ double Fmax=0; for(int i=0; i<natoms;i++){ _max( Fmax, aforce[i].norm2() ); }; return Fmax; };
//void moveGD(double dt){ for(int i=0; i<natoms;i++){ apos[i].add_mul( aforce[i],dt); } };

inline bool insertNeigh( int i, int j ){
    int off = i*nneigh_max;
    int*   ngs   = aneighs + off;
    //double* Kngs = Kneighs + off;
    int ii;
    for(ii=0; ii<nneigh_max; ii++){
        if(ngs[ii]<0){ ngs[ii]=j; break; }
    }
    return ii<nneigh_max;
}

void bond2neighs(){
    int nn=natoms*nneigh_max;
    for(int i=0; i<nn; i++){ aneighs[i]=-1; };
    for(int i=0; i<nbonds; i++){
        Vec2i b = bond2atom[i];
        if(b.i<natoms)insertNeigh( b.i, b.j );
        if(b.j<natoms)insertNeigh( b.j, b.i );
    }
    int ipi=0;
    for(int i=0; i<nn; i++){ if(aneighs[i]<0){ aneighs[i]=ipi; ipi++; }; };
}

Vec3d evalPiTorq(){
    Vec3d tq; tq.set(0.0);
    for(int ia=0; ia<nnode; ia++){
        int ioff = ia*nneigh_max;
        for(int i=0; i<nneigh_max; i++){
            int ing = aneighs[ioff+i];
            if(ing<0){
                int ipi=-ing-1;
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
        int* ngs=aneighs + ia*nneigh_max; 
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
        int* ngs=aneighs + ia*nneigh_max; 
        for(int j=0;j<nneigh_max;j++){
            int ja = ngs[j];
            if(ja>=0){ if(ja>nnode) apos[ ja  ].rotate_csa( ca, sa, ax, p0 ); } // cap atoms
            else     {             pipos[-ja-1].rotate_csa( ca, sa, ax     ); } // pi-bonds
        }
    }
}

void checkNaNs(){
    //printf( "checkNaNs\n" );
    ckeckNaN_d(natoms, 3, (double*)apos,  "apos"  );
    ckeckNaN_d(natoms, 3, (double*)fapos, "fapos" );
    ckeckNaN_d(npi, 3, (double*)pipos,    "pipos" );
    ckeckNaN_d(npi, 3, (double*)fpipos,   "fpipos");
    //for(int i=0; i<natoms; i++){
    //    if( isnan(fapos[i].x) || isnan(fapos[i].x) || isnan(fapos[i].x) ){ printf( "fatoms[%i] is NaN (%g,%g,%g)\n", i, fapos[i].x,fapos[i].y,fapos[i].z ); };
    //}
    //for(int i=0; i<npi; i++){
    //    if( isnan(fpipos[i].x) || isnan(fpipos[i].x) || isnan(fpipos[i].x) ){ printf( "fpipos[%i] is NaN (%g,%g,%g)\n", i, fpipos[i].x,fpipos[i].y,fpipos[i].z ); };
    //}
}

void printBond(int i){ printf( "bond[%i|%i,%i] l0 %g k %g \n", i, bond2atom[i].i,bond2atom[i].j, bond_l0[i], bond_k[i] ); }
void printBonds(){
    printf( "MMFFsp3::printBonds() : \n" );
    for(int i=0;i<nbonds;i++){ printBond(i); }
}

void printNeigh(int ia){
    printf( "atom[%i] neighs[", ia );
    for(int j=0;j<nneigh_max;j++){
        int ij=ia*nneigh_max + j;
        printf( "%i,", aneighs[ij] );
    }
    printf( "] K(", ia );
    for(int j=0;j<nneigh_max;j++){
        int ij=ia*nneigh_max + j;
        printf( "%g,", Kneighs[ij] );
    }
    printf( ")\n");
}
void printNeighs(){
    printf( "MMFFsp3::printNeighs() : \n" );
    for(int i=0;i<nnode;i++){ printNeigh(i); }
}

double optimalTimeStep(double m=1.0){
    double Kmax = 1.0;
    for(int i=0; i<nbonds; i++){ Kmax=fmax(Kmax, bond_k[i] ); }
    return M_PI*2.0*sqrt(m/Kmax)/10.0;  // dt=T/10;   T = 2*pi/omega = 2*pi*sqrt(m/k)
}


};


#endif
