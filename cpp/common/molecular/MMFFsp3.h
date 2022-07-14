
#ifndef MMFFsp3_h
#define MMFFsp3_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

// ======================
// ====   MMFFsp3
// ======================

class MMFFsp3{ public:
    static constexpr const int nneigh_max = 4;
    int  nDOFs=0,natoms=0,nnode=0,ncap=0,nbonds=0,npi=0;
    bool bPBC=false;
    double Eb,Ea,Ep;

    double * DOFs  = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    bool doPi=0;
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

    bool bDEBUG_plot=0;
    int iDEBUG_n=0,iDEBUG_pick=0;

void realloc( int nnode_, int nbonds_, int npi_, int ncap_, bool bNeighs=true ){
    nnode=nnode_; ncap=ncap_; nbonds=nbonds_; npi=npi_; 
    natoms=nnode_+ncap_; 
    nDOFs = (natoms + npi)*3;
    printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    _realloc( DOFs     , nDOFs );
    _realloc( fDOFs    , nDOFs );
    ipi0=natoms;
    _realloc( atype, natoms );
    apos   = (Vec3d*)DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;

    _realloc( bond2atom , nbonds );
    _realloc( bond_l0   , nbonds );
    _realloc( bond_k    , nbonds );
    //_realloc( lbond     , nbonds );
    //_realloc( hbond     , nbonds );

    if(bNeighs){
        int nn = nnode*nneigh_max;
        _realloc( aneighs, nn );
        _realloc( Kneighs, nn );
    }

}

void cleanAtomForce(){ for(int i=0; i<natoms; i++){ fapos [i].set(0.0); } }
void cleanPiForce  (){ for(int i=0; i<npi;    i++){ fpipos[i].set(0.0); } }
void normalizePi   (){ for(int i=0; i<npi;    i++){ pipos[i].normalize(); } }

// ============== Evaluation


int i_DEBUG = 0;


inline double evalPi( const Vec2i& at, int ipi, int jpi, double K ){  // interaction between two pi-bonds
    // E = -K*<pi|pj>^2
    //Vec3d d = apos[at.a]-apos[at.b];
    double c = pipos[ipi].dot(  pipos[jpi] );
    return K * c*c;
}

inline double evalAngle_dist( int ia, int ing, int jng, double K ){
    //K = -1.0;
    //  E = -K*|pi-pj|
    // ToDo: This may be made more efficient if we store hbonds
    Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();  // j->i
    Vec3d hi = apos[ing] - apos[ia];    double ri  = hi.normalize();
    Vec3d hj = apos[jng] - apos[ia];    double rj  = hj.normalize();
    Vec3d fi  = d* K;
    Vec3d fj  = d*-K;
    
    //Vec3d fiT = hi* hi.dot(fi);   fi.sub(fiT);
    //Vec3d fjT = hj* hj.dot(fj);   fj.sub(fjT);

    //printf(  "ang[%i|%i,%i] c1 %g c2 %g \n", ia, ing, jng, hi.dot(fi), hj.dot(fj)  );
    
    glColor3f(1.,0.,0.);
    Draw3D::drawVecInPos(  fi, apos[ing] );
    Draw3D::drawVecInPos(  fj, apos[jng] );
    //Draw3D::drawVecInPos(  fjT+fiT, apos[ia] );
    
    //fapos[ia] .add(fiT);
    //fapos[ia] .add(fjT);
    fapos[ing].add(fi );
    fapos[jng].add(fj );
    return K*rij;
}

inline double evalAngle_cos(  int ia, int ing, int jng, double K ){
    //Vec3d h1,h2;
    //Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();
    Vec3d h1 = apos[ing] - apos[ia];    double ir1  = 1/h1.normalize();
    Vec3d h2 = apos[jng] - apos[ia];    double ir2  = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double c_ = c+1.0;
    double E = K*c_*c_;
    double fang = -K*c_*2;
    hf1.mul( fang*ir1 );
    hf2.mul( fang*ir2 );
    
    //if(bDEBUG_plot){
    if(ia==7){
        //printf( "evalAngle_cos [%i|%i,%i,%i] c(%g,%g) \n", iDEBUG_pick, ia, ing, jng, hf1.dot(h1), hf1.dot(h1)  );
        //printf( "evalAngle_cos [%i|%i,%i,%i] c %g c_ %g fang_ %g E %g \n", iDEBUG_pick, ia, ing, jng, c, c_, fang, E );
        //printf( "evalAngle_cos bDEBUG_plot [%i|%i,%i] iDEBUG_pick %i \n", ia, ing, jng, iDEBUG_pick );
        glColor3f(1.,0.,0.);
        Draw3D::drawVecInPos(  hf1, apos[ing] );
        Draw3D::drawVecInPos(  hf2, apos[jng] );
        Draw3D::drawVecInPos(  (hf1+hf2)*-1.0, apos[ia] );
    }

    fapos[ing].add( hf1     );
    fapos[jng].add( hf2     );
    fapos[ia].sub( hf1+hf2 );
    return E;
}

inline double evalSigmaPi( int ia, int ing, int ipi, double K ){
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
    fapos[ing ].add( hf1 );
    fpipos[ipi].add( hf2 );
    fapos[ia].sub( hf1+hf2 );
    return E;
}

inline double evalPiPi( int ia, int ipi, int jpi, double K ){
    //Vec3d h1,h2;
    //Vec3d d  = apos[ing] - apos[jng];   double rij = d .normalize();
    Vec3d h1 = pipos[ipi];   // double ir1  = 1/h1.normalize();
    Vec3d h2 = pipos[jpi];   //double ir2  = 1/h2.normalize();
    double c = h1.dot(h2);
    Vec3d hf1,hf2; // unitary vectors of force - perpendicular to bonds
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    double E = K*c*c;
    double fang = -K*c*2;
    hf1.mul( fang );
    hf2.mul( fang );
    fapos [ipi].add( hf1 );
    fpipos[jpi].add( hf2 );
    fapos[ia].sub( hf1+hf2 );
    return E;
}

double eval_bond(int ib){
    //printf( "bond %i\n", ib );
    Vec2i iat = bond2atom[ib];
    Vec3d f; f.set_sub( apos[iat.y], apos[iat.x] );
    //if(pbcShifts)f.add( pbcShifts[ib] );
    //printf( "bond[%i|%i,%i] (%g,%g,%g) (%g,%g,%g) \n", ib,iat.a,iat.b, apos[iat.x].x, apos[iat.x].y, apos[iat.x].z, apos[iat.y].x, apos[iat.y].y, apos[iat.y].z );
    double l = f.normalize();
    //printf( "bond[%i|%i,%i] (%g,%g,%g) %g \n", ib,iat.a,iat.b, f.x, f.y, f.z, l );
    //lbond [ib] = l;
    //hbond [ib] = f;
    //if(bondMasked)if(bondMasked[ib])return 0;
    const double k = bond_k[ib];
    double dl = (l-bond_l0[ib]);
    f.mul( dl*k*2 );

    //glColor3f(1.,0.,0.);
    //Draw3D::drawVecInPos(  f, apos[iat.x] );
    //Draw3D::drawVecInPos(  f, apos[iat.y] );

    fapos[iat.x].add( f );
    fapos[iat.y].sub( f );
    double E = k*dl*dl;
    /*
    // --- Pi Bonds
    double Epi = 0;
    if(doPi){ // interaction between pi-bonds of given atom
        int i0=(iat.a+1)*nneigh_max;
        int j0=(iat.b+1)*nneigh_max;
        for(int i=-1;i>=-2;i--){
            int ipi=aneighs[i0+i];
            if(ipi<ipi0) continue;
            for(int j=-1;j>=-2;j--){
                int jpi=aneighs[j0+j];
                if(jpi<ipi0) continue;
                Epi += evalPi( iat, ipi0-ipi,ipi0-jpi, Kneighs[i0+i]*Kneighs[j0+j] );
            }
        }
    }
    E += Epi;
    */
    return E;
}

double eval_neighs(int ia){
    double E=0;
    int ioff = ia*nneigh_max;
    //bool bDEBUG=(ia==7);
    //if(bDEBUG)printf( "#-atom[%i] [%i,%i,%i,%i] \n", ia, aneighs[ioff+0],aneighs[ioff+1],aneighs[ioff+2],aneighs[ioff+3] );
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
            K = 1.0;
            if( bipi ){
                if(bjpi){ return evalPiPi     ( ia, -ing-1, -jng-1, K ); }
                else    { return evalSigmaPi  ( ia, jng, -ing-1, K );    }  
            }else{
                if(bjpi){ return evalSigmaPi  ( ia, ing, -jng-1, K ); }
                else    { return evalAngle_cos( ia, ing, jng, K );    }  
            }
        }
    }
    return E;
}

double eval_bonds(){
    double E=0;
    for(int ib=0; ib<nbonds; ib++){  E+=eval_bond(ib); }
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
    if(bClean){
        cleanAtomForce();   //printf( "DEBUG MMFFsp3.eval() 2 \n" );
        cleanPiForce();     //printf( "DEBUG MMFFsp3.eval() 3 \n" );
    }
    normalizePi();
    Eb = eval_bonds();              //printf( "DEBUG MMFFsp3.eval() 4 \n" );
    Ea = eval_neighs();           //printf( "DEBUG MMFFsp3.eval() 4 \n" );
    return Eb+Ea+Ep;
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

void printBonds(){
    for(int i=0;i<nbonds;i++){
        printf( "bond[%i|%i,%i] l0 %g k %g \n", i, bond2atom[i].i,bond2atom[i].j, bond_l0[i], bond_k[i] );
    }
}


};


#endif
