/*
Non-Bonded Force-Field
 should be easily plugged with any molecular dynamics, by either sharing pointer to same data buffer, or by copying data
*/

#ifndef NBFF_h
#define NBFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

#include "Forces.h"

bool checkPairsSorted( int n, Vec2i* pairs ){
    int ia=-1,ja=-1;
    for(int i=0;i<n; i++){
        const Vec2i& b = pairs[i];
        //printf( "pair[%i] %i,%i | %i %i  | %i %i %i ", i, b.i, b.j,   ia,ja ,   b.i>=b.j,  b.i<ia, b.j<=ja );
        if(b.i>=b.j){ return false; }
        if(b.i<ia)  { return false; }
        else if (b.i>ia){ia=b.i; ja=-1; };
        if(b.j<=ja){ return false; }
        ja=b.j;
    }
    return true;
}

// =========== SR repulsion functions

//   good SR repulsion function should have spike in the center ... this is the case of Lorenz
//   (R2-r2)^2
//   (R2-r2)*((R2/r2)-1)^2
//   ((R2/r2)-1)^2
//   (((R2+w2)/(r2+w2))-1)^2

class AtomicSystem{ public:
    int      n;
    int*     types; // non-equivalent atoms, pointer to particular atom in fitting parameters
    Vec3d*   ps;
    //Vec3d* fs;
    //double E;
    AtomicSystem (      ){ n=0;  types=0;            ps=0;               }
    AtomicSystem (int n_){ n=n_; types=new int  [n]; ps   =new Vec3d[n]; }
    ~AtomicSystem(      ){ if(types)delete [] types; if(ps)delete [] ps; }

    void atomsToXYZ(FILE* fout){
        for(int i=0; i<n; i++){
            fprintf( fout, "%i %20.10f %20.10f %20.10f\n", types[i], ps[i].x,ps[i].y,ps[i].z );
        }
    }
};

//class NBsystem : AtomicSystem { public: // Can be Child of 
class NBsystem{ public: // Can be Child of 
    int n;
    int    *atypes=0; // Not necessarily used
    Vec3d  *REQs=0;
    Vec3d  *ps=0;
    Vec3d  *fs=0;
    Vec3d  *PLQs=0;    // used only in combination with GridFF
    Quat4i *neighs=0;  //

    void evalPLQs(double K){
        for(int i=0; i<n; i++){
            //printf( "makePLQs[%i] \n", i );
            //printf( "makePLQs[%i] REQ(%g,%g,%g) \n", i, REQs[i].x,REQs[i].y,REQs[i].z);
            PLQs[i]=REQ2PLQ( REQs[i], K );
            //printf( "makePLQs[%i] REQ(%g,%g,%g) PLQ(%g,%g,%g)\n", i, REQs[i].x,REQs[i].y,REQs[i].z,  PLQs[i].x,PLQs[i].y,PLQs[i].z );
        }
        //printf("NBsystem::makePLQs() DONE => exit(0) \,"); exit(0);
    }

    void makePLQs(double K){
        _realloc(PLQs,n);
        evalPLQs(K);
    }

    void fromRigid( Vec3d* ps0, const Vec3d& p0, const Mat3d& rot ){ for(int i=0; i<n; i++){ rot.dot_to_T( ps0[i], ps[i] ); ps[i].add(p0); } }
    void torq     ( const Vec3d& p0, Vec3d& tq ){ for(int i=0; i<n; i++){ Vec3d d; d.set_sub(ps[i],p0); tq.add_cross(fs[i],d); } }

    void shift( Vec3d d ){ for(int i=0; i<n; i++){  ps[i].add(d); } }


    double evalLJQs(){
        const int N=n;
        double E=0;
        for(int i=0; i<N; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi = ps[i];
            const Vec3d& REQi = REQs[i];
            for(int j=i+1; j<N; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
                E += addAtomicForceLJQ( ps[j]-pi, fij, REQij );
                fs[j].sub(fij);
                fi   .add(fij);
            }
            fs[i].add(fi);
        }
        return E;
    }

    double evalLJQs_ng4( const int* neighs ){
        const int N=n;
        double E=0;
        for(int i=0; i<N; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi = ps[i];
            const Vec3d& REQi = REQs[i];
            const int* ngs = neighs+i*4;
            for(int j=i+1; j<N; j++){    // atom-atom (no self interaction, no double-counting)
                if( (ngs[0]==j)||(ngs[1]==j)||(ngs[2]==j)||(ngs[3]==j) ) continue;
                Vec3d fij = Vec3dZero;
                Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
                E += addAtomicForceLJQ( ps[j]-pi, fij, REQij );
                fs[j].sub(fij);
                fi   .add(fij);
            }
            fs[i].add(fi);
        }
        return E;
    }

    double evalLJQs_PBC( const Mat3d& lvec, Vec3i nPBC=Vec3i{1,1,1} ){
        const int N=n;
        double E=0;
        int npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1) -1;
        Vec3d shifts[npbc]; // temporary store for lattice shifts
        int ipbc=0;
        for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){ 
            if((ia==0)&&(ib==0)&&(ic==0))continue; // skipp pbc0
            shifts[ipbc] = (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic);   
            ipbc++; 
        }}}
        for(int i=0; i<N; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi = ps[i];
            const Vec3d& REQi = REQs[i];
            for(int j=i; j<N; j++){    // atom-atom (yes self interaction, no double-counting)
                Vec3d fij = Vec3dZero;
                Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
                for(ipbc=0; ipbc<npbc; ipbc++){
                    E += addAtomicForceLJQ( ps[j]-pi-shifts[ipbc], fij, REQij );
                }
                fs[j].sub(fij);
                fi   .add(fij);
            }
            fs[i].add(fi);
        }
        printf( "npbc %i ipbc %i E %g \n", npbc, ipbc, E );
        return E;
    }

    double evalLJQs_ng4_PBC( Quat4i* neighs, Quat4i* neighCell, const Mat3d& lvec, Vec3i nPBC=Vec3i{1,1,1}, double Rdamp=1.0 ){
        double R2damp = Rdamp*Rdamp;
        double E=0;    
        int        npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
        const bool bPBC = npbc>0;
        Vec3d shifts[npbc]; // temporary store for lattice shifts
        int ipbc=0;
        if(bPBC>1){
            for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){ 
                shifts[ipbc] = (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic);   
                ipbc++; 
            }}}
        }
        //for(int i=0; i<n; i++)printf( "CPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQs[i].x,REQs[i].y,REQs[i].z );
        for (int i=0; i<n; i++ ){
            Vec3d fi = Vec3dZero;
            Vec3d pi = ps[i];
            const Vec3d& REQi = REQs     [i];
            const Quat4i ng   = neighs   [i];
            const Quat4i ngC  = neighCell[i];
            //for (int j=i+1; j<n; j++){
            //if(i==4){ printf( "CPU_LJQ[%i] ng(%i,%i,%i,%i) ngC(%i,%i,%i,%i) npbc=%i\n", i, ng.x,ng.y,ng.z,ng.w,   ngC.x,ngC.y,ngC.z,ngC.w, npbc ); } 
            //printf( "CPU_LJQ[%i] ng(%i,%i,%i,%i) ngC(%i,%i,%i,%i) npbc=%i\n", i, ng.x,ng.y,ng.z,ng.w,   ngC.x,ngC.y,ngC.z,ngC.w, npbc );
            for (int j=0; j<n; j++){  // DO ALL TO ALL (to be consistent with GPU)
                if(i==j)continue;
                const Vec3d dp = ps[j]-pi;
                Vec3d fij=Vec3dZero;
                Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
                const bool bBonded = ((j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w));
                //const bool bBonded = false;
                for(ipbc=0; ipbc<npbc; ipbc++){
                    if(bBonded){
                        if(
                              ((j==ng.x)&&(ipbc==ngC.x))
                            ||((j==ng.y)&&(ipbc==ngC.y))
                            ||((j==ng.z)&&(ipbc==ngC.z))
                            ||((j==ng.w)&&(ipbc==ngC.w))
                        ){
                            //printf("skip[%i,%i]ipbc=%i\n", i, j, ipbc );
                            continue; // skipp pbc0
                        }
                    }
                    //E += addAtomicForceLJQ( dp + shifts[ipbc], fij, REQij );
                    Vec3f fij_; E+=getLJQ( (Vec3f)(dp+shifts[ipbc]), (Vec3f)REQij, R2damp, fij_ );
                    fij.add((Vec3d)fij_);
                    //if(i==4){ printf( "CPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g\n" , i,j, ipbc, fij_.x,fij_.y,fij_.z, R2damp, REQij.x,REQij.y,REQij.z, dp.norm() ); } 
                    //printf( "CPU_LJQ[%i,%i|%i] fj(%g,%g,%g)\n" , i,j, ipbc, fij_.x,fij_.y,fij_.z );
                }
                //if(i==4){ printf( "CPU_LJQ[%i,%i]   fj(%g,%g,%g) bBonded %i \n" , i,j, fij.x,fij.y,fij.z, bBonded ); } 
                //fs[j].sub(fij);
                fi   .add(fij);
            }
            fs[i].add(fi);
        }
        return E;
    }

    double evalLJQ( NBsystem& B, const bool bRecoil ){
        double E=0;
        //printf("DEBUG NBFF_AB.evalLJQ() n,m %i %i \n", n,m);
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            Vec3d Api = ps[i];
            const Vec3d& AREQi = REQs[i];
            for(int j=0; j<B.n; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Vec3d REQij; combineREQ( B.REQs[j], AREQi, REQij );
                E += addAtomicForceLJQ( B.ps[j]-Api, fij, REQij );
                if(bRecoil) B.fs[j].sub(fij);
                fi.add(fij);
            }
            fs[i].add(fi);
        }
        return E;
    }

    double evalMorse( NBsystem& B, const bool bRecoil, double K=-1.0, double RQ=1.0 ){
        double E=0;
        //printf("DEBUG evalMorse() n,B.n %i %i \n", n,B.n );
        double R2Q=RQ*RQ;
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            Vec3d Api = ps[i];
            const Vec3d& AREQi = REQs[i];
            for(int j=0; j<B.n; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Vec3d REQij; combineREQ( B.REQs[j], AREQi, REQij );

                E += addAtomicForceMorseQ( B.ps[j]-Api, fij, REQij.x, REQij.y, REQij.z, K, R2Q );
                if(bRecoil) B.fs[j].sub(fij);
                fi.add(fij);
            }
            //printf( "fs[%i]( %g | %g,%g,%g)  \n", i, fi.norm(),  fi.x,fi.y,fi.z );
            fs[i].add(fi);
            //fs[i].sub(fi);
        }
        //printf( "EvdW %g \n", E );
        return E;
    }

    double evalMorsePBC( NBsystem& B, const Mat3d& cell, Vec3i nPBC=(Vec3i){1,1,0}, double K=-1.0, double RQ=1.0 ){
        //nPBC = {0,0,0};
        //printf( "NBsystem::evalMorsePBC() nPBC(%i,%i,%i) \n", nPBC.x, nPBC.y, nPBC.z );
        //printf( "cell.a(%g,%g,%g) cell.b(%g,%g,%g) cell.c(%g,%g,%g) \n", cell.a.x,cell.a.y,cell.a.z,   cell.b.x,cell.b.y,cell.b.z,   cell.c.x,cell.c.y,cell.c.z );
        double E=0;
        //printf("DEBUG evalMorse() n,B.n %i %i \n", n,B.n );
        double R2Q=RQ*RQ;
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi_ = ps[i];
            const Vec3d& AREQi = REQs[i];
            //if(i==0)printf("pi_(%g,%g,%g) \n",  pi_.x,pi_.y,pi_.z );
            for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                Vec3d  pi = pi_ + (cell.a*ia) + (cell.b*ib) + (cell.c*ic);
                //if(i==0)printf("pi[%2i,%2i,%2i] = (%g,%g,%g) \n", ia,ib,ic,   pi.x,pi.y,pi.z );
                for(int j=0; j<B.n; j++){    // atom-atom
                    Vec3d fij = Vec3dZero;
                    Vec3d REQij; combineREQ( B.REQs[j], AREQi, REQij );
                    E += addAtomicForceMorseQ( B.ps[j]-pi, fij, REQij.x, REQij.y, REQij.z, K, R2Q );
                    //E += addAtomicForceMorseQ( B.ps[j]-pi, fij, REQij.x, REQij.y, 0, K, R2Q );
                    //E += addAtomicForceQ_R2  ( B.ps[j]-pi, fij, REQij.z, K, R2Q );

                    //E=0; fij = B.ps[j]-pi;  // Test - dr
                    //if(i==0)printf("fi(%g,%g,%g) \n",  fij.x,fij.y,fij.z );

                    fi.add(fij);
                }
            }}} // nPBC
            //if(i==0){ printf( "CPU atom[%i]  fe_Cou(%g,%g,%g|%g)  REQKi.z %g \n", i, fi.x,fi.y,fi.z,E, AREQi.z ); }
            fs[i].add(fi);
        }
        return E;
    }

    double evalMorsePLQ( NBsystem& B, Mat3d& cell, Vec3i nPBC, double K=-1.0, double RQ=1.0 ){
        //printf( "NBsystem::evalMorsePLQ() PLQs %li \n", (long)PLQs, K, RQ );
        double E=0;
        //printf( "NBFF nPBC(%i,%i,%i) K %g RQ %g R2Q %g plq.z %g \n", nPBC.x,nPBC.y,nPBC.z, K, sqrt(R2Q), R2Q, plq.z );
        double R2Q=RQ*RQ;
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi = ps[i];
            //const Vec3d& REQi = REQs[i];
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int j=0; j<B.n; j++){
                Vec3d dp0; dp0.set_sub( pi, B.ps[j] );
                Vec3d REQj = B.REQs[j];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + cell.a*ia + cell.b*ib + cell.c*ic;
                    //Vec3d  dp     = dp0;
                    double r2     = dp.norm2();
                    double r      = sqrt(r2);
                    // ----- Morse
                    double e      = exp( K*(r-REQj.x) );
                    double de     = K*e*REQj.y*-2/r;
                    double eM     = e*REQj.y;
                    // ---- Coulomb
                    double ir2    = 1/(r2+R2Q);
                    double ir     = sqrt(ir2);
                    double eQ     = COULOMB_CONST*REQj.z*ir;
                    //if((i==0)&&(j==0))printf("dp(%g,%g,%g) REQj(%g,%g,%g) r %g e %g de %g \n",  dp.x,dp.y,dp.z,  REQj.x,REQj.y,REQj.z, r, e, de );
                    // --- store
                    qp.e+=eM*e; qp.f.add_mul( dp, de*e   ); // repulsive part of Morse
                    ql.e+=eM*2; ql.f.add_mul( dp, de     ); // attractive part of Morse
                    qe.e+=eQ;   qe.f.add_mul( dp, eQ*ir2 ); // Coulomb
                    //printf(  "evalMorsePLQ() k %g r %g e %g E0 %g E %g \n", K, r, e, REQj.y*plq.x/exp( K*(1.487)), qp.e*plq.x );
                    //if(i==0)printf( "[%i] qp(%g,%g,%g)  ql(%g,%g,%g)  qe(%g,%g,%g)\n", j, qp.x,qp.y,qp.z,   ql.x,ql.y,ql.z,   qe.x,qe.y,qe.z );
                }}}
            }
            Vec3d plq = PLQs[i];
            //if(i==0)printf( "plq[0](%g,%g,%g) qp(%g,%g,%g)  ql(%g,%g,%g)  qe(%g,%g,%g)\n", plq.x,plq.y,plq.z,  qp.x,qp.y,qp.z,   ql.x,ql.y,ql.z,   qe.x,qe.y,qe.z );
            Quat4d fe = qp*plq.x + ql*plq.y + qe*plq.z;
            fs[i].add(fe.f);
            E       +=fe.e;
        }
        return E;
    }


    double evalR( NBsystem& B ){
        double E=0;
        for(int i=0; i<n; i++){
            //Vec3d fi   = Vec3dZero;
            Vec3d Api = ps[i];
            for(int j=0; j<B.n; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Vec3d d = B.ps[j]-Api;
                double r = d.norm();
                E += fmax( 0, REQs[j].x-r );
                //fi.add(fij);
            }
            //fs[i].add(fi);
        }
        return E;
    }

    double evalNeighs( double RQ=1.0, double K=-1.5 ){
        double E=0;
        double R2Q=RQ*RQ;
        for(int i=0; i<n; i++){
            Vec3d  fi = Vec3dZero;
            Vec3d  pi = ps[i];
            Quat4i ngi = {-1,-1,-1,-1};
            if(neighs)ngi=neighs[i];
            //pi.add( shift );
            const Vec3d& REQi = REQs[i];
            for(int j=i+1; j<n; j++){    // atom-atom
                if( (j==ngi.x)||(j==ngi.y)||(j==ngi.z)||(j==ngi.w) ) continue;
                Vec3d fij = Vec3dZero;
                Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
                Vec3d dp=ps[j]-pi;
                //if(i==0){ idebug=1; }else{ idebug=0; };
                double ei = addAtomicForceMorseQ( dp, fij, REQij.x, REQij.y, REQij.z, K, R2Q );    E+=ei;
                //double ei = addAtomicForceQ_R2( dp, fij, REQij.z, K, R2Q );    E+=ei;
                //E += addAtomicForceLJQ   ( ps[j]-pi, fij, REQij );
                //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( fij      , pi );
                //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( fij*-1.0f, ps[j] );
                //glColor3f(1.0f,0.0f,1.0f); Draw3D::drawLine( pi, ps[j] );
                //if(i==6){ printf("CPU[%i,%i] dp(%g,%g,%g) fe(%g,%g,%g|%g)\n", i,j,  dp.x,dp.y,dp.z,   fij.x,fij.y,fij.z,ei ); }
                fs[j].sub(fij);
                fi   .add(fij);
            }
            fs[i].add(fi);
        }
        return E;
    }

    void print(){
        printf("NBsystem(n=%i):\n");
        for(int i=0; i<n; i++){
            if(atypes){ printf("nb_atom[%i] REQ(%g,%g,%g) pos(%g,%g,%g) atyp %i \n", i, REQs[i].x,REQs[i].y,REQs[i].z,  ps[i].x,ps[i].y,ps[i].z, atypes[i] ); }
            else      { printf("nb_atom[%i] REQ(%g,%g,%g) pos(%g,%g,%g) \n", i, REQs[i].x,REQs[i].y,REQs[i].z,  ps[i].x,ps[i].y,ps[i].z ); }
        }
    }

    void bindOrRealloc(int n_, Vec3d* ps_, Vec3d* fs_, Vec3d* REQs_ ){
        n=n_;
        _bindOrRealloc(n,ps_  ,ps  );
        _bindOrRealloc(n,fs_  ,fs  );
        _bindOrRealloc(n,REQs_,REQs);
    }

};

class NBFF : public NBsystem{ public:
// non bonded forcefield
    //int n       = 0;
    //Vec3d* REQs = 0;
    //Vec3d* ps   = 0;
    //Vec3d* fs   = 0;

    int nmask   = 0;

    Vec2i* pairMask = 0; // should be ordered ?

    double sr_w      = 0.1;
    double sr_dRmax  = 1.0; // maximum allowed atom movement before neighborlist update
    double sr_K      = 1.0; // collision stiffness
    double sr_Rscale = 0.5; // rescale atoms for short range
    Vec3d* psBack = 0;   // backup positon for neighbor list
    std::vector<Vec2i> neighList; // possibly too slow ?

void realloc(int n_, int nmask_){
    n=n_;
    nmask=nmask_;
    _realloc(REQs,n);
    _realloc(ps  ,n);
    _realloc(fs  ,n);
    _realloc(pairMask ,nmask);
}

void bindOrRealloc(int n_, int nmask_, Vec3d* ps_, Vec3d* fs_, Vec3d* REQs_, Vec2i* pairMask_ ){
    n=n_;
    nmask=nmask_;
    _bindOrRealloc(n,ps_  ,ps  );
    _bindOrRealloc(n,fs_  ,fs  );
    _bindOrRealloc(n,REQs_,REQs);

    _bindOrRealloc(nmask,pairMask_,pairMask);
}

void setREQs(int i0,int i1, const Vec3d& REQ){
    for(int i=i0;i<i1;i++){ REQs[i]=REQ; }
}

void printAtomParams(){
    printf( " # NBFF::printAtomParams() \n" );
    for(int i=0;i<n;i++){
        printf(  "atom[%i] R %g E %g Q %g \n", i, REQs[i].x, REQs[i].y, REQs[i].z );
    }
}

void cleanForce(){
    for(int i=0; i<n; i++){ fs[i].set(0.); }
}

double evalLJQ_sortedMask( const Vec3d& shift=Vec3dZero ){
    //printf(  "evalLJQ_sortedMask \n" );
    int im=0;
    const int N=n;
    double E=0;
    for(int i=0; i<N; i++){
        Vec3d fi = Vec3dZero;
        Vec3d pi = ps[i];
        pi.add( shift );
        const Vec3d& REQi = REQs[i];
        for(int j=i+1; j<N; j++){    // atom-atom
            // --- mask some atom pairs (e.g. those which are bonded), as long as atoms are sorted we can do it efficiently
            if( (im<nmask)&&(i==pairMask[im].i)&&(j==pairMask[im].j) ){
                im++; continue;
            }
            Vec3d fij = Vec3dZero;
            Vec3d REQij; combineREQ( REQs[j], REQi, REQij );
            double E = addAtomicForceLJQ( ps[j]-pi, fij, REQij );

            /*
            if(E>100.0){
                //printf( "%i %i  %g \n", i, j, E );
                glColor3f(0.0,1.0,0.0); Draw3D::drawLine( ps[j],pi);
                //printf("%i %i %g \n", i, j, fij.norm());
                glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fij*-1000,ps[j]);
                glColor3f(1.0,0.0,0.0); Draw3D::drawVecInPos( fij*1000 ,ps[i]);
            }
            */

            fs[j].sub(fij);
            fi   .add(fij);

        }
        fs[i].add(fi);
    }
    //exit(0);
    return E;

}

double evalLJQ_pbc( Mat3d lvec, Vec3i npbc ){
    double E=0;
    for(int ix=-npbc.x;ix<=npbc.x;ix++){
        for(int iy=-npbc.y;iy<=npbc.y;iy++){
            for(int iz=-npbc.z;iz<=npbc.z;iz++){
                //E+=evalLJQ_shifted( lvec.a*ix + lvec.b*iy + lvec.c*iz );
                E+=evalLJQ_sortedMask( lvec.a*ix + lvec.b*iy + lvec.c*iz );
            }
        }
    }
    return E;
}

// ============= Short Range Force using neighbor list

bool checkListValid(){
    double R2=sq(sr_dRmax);
    const int N=n;
    for(int i=0; i<N; i++){
        Vec3d d; d.set_sub(ps[i],psBack[i]);
        if(d.norm2()>R2){ return false; }
    }
    return true;
}

void makeSRList(){
    //Vec2i* m=pairMask;
    int im=0;
    const int N=n;
    //printf( "N %i \n" );
    neighList.clear();
    for(int i=0; i<N; i++){
        const Vec3d& pi = ps[i];
        psBack[i]=pi;
        double Ri = REQs[i].x;
        const Vec3d& REQi = REQs[i];
        for(int j=i+1; j<N; j++){    // atom-atom
            if( (im<nmask)&&(i==pairMask[im].i)&&(j==pairMask[im].j) ){
                im++; continue;
            }
            Vec3d d   = ps[j]-pi;
            double r2 = d.norm2();
            double Rij = (Ri + REQs[j].x)*sr_Rscale + sr_dRmax;
            if( r2<Rij*Rij ){
                neighList.push_back({i,j});  //  TODO: is this fast enough?
            }
        }
    }
}

double evalSRlist(){
    double E=0;
    double w2 = sq(sr_w);
    for(const Vec2i& ij : neighList ){
        Vec3d fij = Vec3dZero;
        double Rij = (REQs[ij.i].x+REQs[ij.i].y)* sr_Rscale;
        //addForceR2   ( ps[ij.j]-ps[ij.i], fij, Rij*Rij, sr_K     );
        //addForceR2inv( ps[ij.j]-ps[ij.i], fij, Rij*Rij, sr_K, w2 );
        E += addForceR2inv( ps[ij.j]-ps[ij.i], fij, Rij*Rij, sr_K, w2 );
        fs[ij.i].sub(fij);
        fs[ij.j].add(fij);
    }
    return E;
}

};

#endif

