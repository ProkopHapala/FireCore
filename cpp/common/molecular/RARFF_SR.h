
#ifndef RARFF_SR_h
#define RARFF_SR_h
// This is modified from RARFFarr.h

#include <vector>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "SMat3.h"
#include "quaternion.h"
#include "Forces.h"

#include "Buckets3D.h"

/*
Rigid Atom Reactive Force-field
===============================

Problem:  binding of several atoms to one center
Solution: bonds repel for opposite side (back-side)

Electron Pairs
 - To prevent many atoms

Passivation atoms
 - except "scafold atoms" (C,O,N) there are termination atoms (-H, Cl) wich r

Optimization:
    - We can separate rotation and movement of atoms, rotation should be much faster
    - This makes only sense if we use forcefield with cheap evaluation of torque for fixed distance
        - we can perhaps save sqrt() calculation


Simple reactive force-field for mix of sp2 and sp3 hybridized atoms. Energy is based on Morse potential where the attractive term is multiplied by product of cosines between oriented dangling-bonds  sticking out  of atoms (the white sticks).

More specifically for two atoms with positions pi,pj (dij=pi-pj, rij=|dij|, hij = dij/rij )

energy is computed as
E = exp( -2a(rij-r0)) - 2 K exp( -a(rij-r0) )

where rij is distance  K = (ci.cj.cj) where ci=dot(hi,hij), cj=dot(hj,hij), cij=dot(hi,hj) and where hi, hj are normalized vector in direction of bonds

*/

#define N_BOND_MAX 4
#define R1SAFE    1e-8

/*
static const double sp3_hs[] = {
-0.57735026919, -0.57735026919, -0.57735026919,
+0.57735026919, +0.57735026919, -0.57735026919,
-0.57735026919, +0.57735026919, +0.57735026919,
+0.57735026919, -0.57735026919, +0.57735026919
};
*/

static const double sp3_hs[] = {
1.000000000000,  0.00000000000,  0.00000000000,
-0.33380685923,  0.94264149109,  0.00000000000,
-0.33380685923, -0.47132074554, -0.81635147794,
-0.33380685923, -0.47132074554, +0.81635147794
};

static const double sp2_hs[] = {
+1.000000000000, -0.00000000000,  0.00000000000,
-0.500000000000, +0.86602540378,  0.00000000000,
-0.500000000000, -0.86602540378,  0.00000000000,
 0.00000000000,   0.00000000000,  1.00000000000     // electron - can be repulsive
};

static const double sp1_hs[] = {
+1.000000000000,  0.00000000000,  0.00000000000,
-1.000000000000,  0.00000000000,  0.00000000000,
 0.000000000000,  1.00000000000,  0.00000000000,    // electron - can be repulsive
 0.000000000000,  0.00000000000,  1.00000000000     // electron - can be repulsive
};


inline void overlapFE(double r, double amp, double beta, double& e, double& fr ){
    //https://www.wolframalpha.com/input/?i=exp(b*x)*(1%2Bb*x%2B(b*x)%5E2%2F3)+derivative+by+x
    double x     = r*beta;
    double expar = amp*exp(-x);
    e  =  expar*(1 +   x + 0.33333333*x*x );
    fr = (expar*(6 + 5*x +            x*x )*beta*0.33333333)/r;
}

struct RigidAtomType{
    int id = -1;
    char name[4]="C";
    int    nbond = N_BOND_MAX;  // number bonds

    //double Rcut;  // total cutoff (also bond cutoff)
    //double Rrep;  // repulsion cutoff
    //double Rbond; // repulsion cutoff
    //double Ebond; // Bond strengh
    //double Wbond; // Wbond - stiffnes of attraction
    //double Arep;  // repulsion strenght

    double Epz    =  0.5;   // pi-bond strenght
    double rbond0 =  0.5;
    double aMorse =  4.0;
    double bMorse = -0.7;
    // TODO : use REQs here ?

    //double c6     = -15.0;
    //double R2vdW  =  8.0;
    Vec3d* bh0s = (Vec3d*)sp3_hs;

    inline void combine(const RigidAtomType& a, const RigidAtomType& b ){
        nbond   = a.nbond;
        aMorse  = a.aMorse * b.aMorse;  // TODO
        bMorse  = a.bMorse + b.bMorse;
        rbond0  = a.rbond0 + b.rbond0;
        //Rcut   = a.Rcut  + b.Rcut;  
        //Rrep   = a.Rrep  + b.Rrep;  
        //Rbond  = a.Rbond  + b.Rbond;
        //Wbond  = a.Wbond  + b.Wbond;
        //Ebond  = a.Ebond * b.Ebond;
        //Erep   = a.Arep  * b.Arep;  
        //printf( "nbond %i Rcut %g Rrep %g Abond %g Arep %g \n", nbond, Rcut, Rrep, Abond, Arep );
    }
    void print(){
        //printf( "nbond  %i Epz %g \n", nbond, Epz );
        //printf( "Bond  Rcut %g Abond %g\n", Rcut, Abond );
        //printf( "Core  Rrep %g Arep  %g\n", Rrep, Arep  );
        printf( "nbond  %i rbond0 %g\n", nbond, rbond0 );
        printf( "aMorse %g bMorse %g Epz %g \n", aMorse, bMorse, Epz );
        //printf( "c6     %g r2vdW  %g\n", c6, R2vdW );
        //exit(0);
    }
    /*
    inline double evalRadial( double r, double fr ){
        if(r>Rb){ // attraction
            const double sc  = 1/( Rcut-Rbond );
            const double iw2 = 1/(Wbond*Wbond);
            double x  = (r-Rbond)*sc;
            double mx = 1-x;
            double  y1 =mx*mx;
            double  yL = 1/(1-x*x*iw2);
            double dy1 =2*(x-1)*sc;
            double dyL =yL*yL*2*x*sc*iw2;
            fr    = y1*dyL + yL*dy1;
            return y1*yL;
        }else{    // contanct - Parabolic appox
            double Rd = Rbond-Rrep;
            double K  = Eb/(Rd*Rd);
            r-=Rbond;
            fr   = K*2*r;
            return K*r*r - Eb;
        }
    }
    */
    //RigidAtomType()=default;
    //RigidAtomType():nbond(nbond_),rbond(rbond_),aMorse(){};

};

struct CapType{
    int    id = -1;
    char name[4]="H";
    double rbond0 =  1.0;
    double RvdW   =  3.0;
    double EvdW   =  0.01;
};

class PairList{ public:
    std::vector<Vec2i> neighList;

    void makeSRList( int na, Vec3d* ps, double R ){
        //printf( "DEBUG makeSRList \n" );
        neighList.clear();
        double R2=R*R;
        for(int i=0; i<na; i++){
            const Vec3d& pi = ps[i];
            for(int j=i+1; j<na; j++){    // atom-atom
                Vec3d d   = ps[j]-pi;
                double r2 = d.norm2();
                if( r2<R2 ){
                    neighList.push_back({i,j});  //  TODO: is this fast enough?
                }
            }
        }
        //printf( "DEBUG makeSRList DONE npair %i \n", neighList.size() );
    }
    void bind( int& n, Vec2i*& buff ){ n=neighList.size(); buff=&neighList[0]; }
};

class RARFF_SR{ public:
    // This is modified from RARFFarr
    double RcutMax = 5.0; // [A]
    double R2cut = RcutMax*RcutMax;
    double RcutCap = 3.0; // [A]
    double R2cutCap = RcutCap*RcutCap;
    Buckets3D map;

    int    npairs = 0;
    Vec2i* pairs  = 0;

    const CapType       * capTypeList = 0;
    const RigidAtomType * typeList    = 0;

    int  AccelType  = 1; // 1=grid, 2=list, 0=brute
    //bool bGridAccel = true;
    bool bRepelCaps = true;
    bool bDonorAcceptorCap = false;

    double lcap       =   1.0;
    double aMorseCap  =   1.0;
    double bMorseCap  =  -0.7;

    double Ecap   = 0.000681; // H
    double Rcap   = 1.4+1.4;  // H
    double r2cap  = Rcap*Rcap;
    double r6cap  = r2cap*r2cap*r2cap;
    double c6cap  = 2*Ecap*(r6cap);
    double c12cap =   Ecap*(r6cap*r6cap);

    double invRotMass = 2.0;

    int natomActive   = 0;
    int natom         = 0;
    //RigidAtomType* types =0;
    //RigidAtom*     atoms =0;

    // atom properties
    const RigidAtomType** types = 0;
    Vec3d*  apos   = 0;
    Quat4d* qrots  = 0;
    Vec3d*  aforce = 0;
    Vec3d*  torqs  = 0;
    Vec3d*  omegas = 0;  // just for MD
    Vec3d*  vels   = 0;  // just for MD

    // aux
    double*        ebonds=0;
    Vec3d*         hbonds=0;
    Vec3d*         fbonds=0;
    int  *         bondCaps = 0;

    bool * fixedAtoms  = 0;
    bool * ignoreAtoms = 0;
    //double F2pos=0;
    //double F2rot=0;

    int n_pairs_tried=0;
    int n_pairs_evaluated=0;

    int countIgnoreAtoms(bool bPrint=false){
        int nignore=0;
        for(int i=0; i<natom; i++){ if(ignoreAtoms[i]) nignore++; }
        if(bPrint)printf( "RARFF_SR::countIgnoreAtoms() natom=%i nignore=%i \n", natom, nignore );
        return nignore;
    }

    int checkAtomsCollapsed( const double Rcut, int bPrint=false ){
        int n=0;
        double R2cut = Rcut*Rcut;
        for(int i=0; i<natom; i++){
            if(ignoreAtoms[i]) continue;
            for(int j=i+1; j<natom; j++){
                if(ignoreAtoms[j]) continue;
                Vec3d d = apos[i]-apos[j];
                double r2 = d.norm2();
                if(r2<R2cut){
                    n++;
                    if(bPrint)printf( "RARFF_SR::checkAtomsCollapsed(%i,%i) |rij|=%g \n", i, j, sqrt(r2) );
                }
            }
        }
        return n;
    }

    int checkAtomsOut( const double Rcut, int bPrint=false){
        int n=0;
        double R2cut = Rcut*Rcut;
        for(int i=0; i<natom; i++){
            if(ignoreAtoms[i]) continue;
            Vec3d d = apos[i];
            double r2 = d.norm2();
            if( (r2>R2cut) || isnan(r2) ){
                n++;
                if(bPrint)printf( "RARFF_SR::checkAtomsCollapsed(%i,%i) |rij|=%g \n", i, sqrt(r2) );
            }
        }
        return n;
    }

    void alloc(int natom_){
        natom=natom_;
        //printf(  "FARFF aloc na %i no %i nDOF %i \n", natom, norb, nDOF );
        _alloc(types  ,natom);
        _alloc(apos   ,natom);
        _alloc(qrots  ,natom);
        _alloc(aforce ,natom);
        _alloc(torqs  ,natom);
        _alloc(omegas ,natom);  // just for MD
        _alloc(vels   ,natom);  // just for MD
        _alloc(ignoreAtoms, natom);
        _alloc(fixedAtoms ,natom);
        _alloc(ebonds ,natom*N_BOND_MAX);
        _alloc(hbonds ,natom*N_BOND_MAX);
        _alloc(fbonds ,natom*N_BOND_MAX);
        _alloc(bondCaps ,natom*N_BOND_MAX);
        natomActive=natom;
        for(int i=0; i<natom; i++){ 
            ignoreAtoms[i]=false; 
            fixedAtoms [i]=false;
        }
    }

    void realloc(int natom_, int nbuff=0 ){
        natomActive=natom_;
        natom=natom_+nbuff;
        //_realloc(atoms,natom   );
        _realloc(types  ,natom);
        _realloc(apos   ,natom);
        _realloc(qrots  ,natom);
        _realloc(aforce ,natom);
        _realloc(torqs  ,natom);
        _realloc(omegas ,natom);
        _realloc(vels   ,natom);
        _realloc(ignoreAtoms, natom);
        _realloc(fixedAtoms,natom);
        _realloc(ebonds ,natom*N_BOND_MAX);
        _realloc(hbonds ,natom*N_BOND_MAX);
        _realloc(fbonds ,natom*N_BOND_MAX);
        _realloc(bondCaps ,natom*N_BOND_MAX);
        for(int i=0;      i<natom_; i++){ ignoreAtoms[i]=false; }
        for(int i=natom_; i<natom;  i++){ ignoreAtoms[i]=true;  }
        for(int i=0; i<natom;       i++){ fixedAtoms [i]=false; }
    }

    void resize( int natom_new ){
        printf( "RARFFarr::resize %i -> %i \n", natom, natom_new );
        int natom_=natom;
        // save old
        const RigidAtomType** types_ = types;
        Vec3d*  apos_   = apos;
        Quat4d* qrots_  = qrots;
        Vec3d*  aforce_ = aforce;
        Vec3d*  torqs_  = torqs;
        Vec3d*  omegas_ = omegas; // just for MD
        Vec3d*  vels_   = vels;   // just for MD
        bool*   ignoreAtoms_ = ignoreAtoms;
        bool*   fixedAtoms_  = fixedAtoms;
        double* ebonds_ = ebonds;
        Vec3d*  hbonds_ = hbonds;
        Vec3d*  fbonds_ = fbonds;
        int  *  bondCaps_  = bondCaps;
        // copy
        alloc(natom_new);
        //printf( "resize natom_ %i natom %i natom_new %i \n", natom_, natom, natom_new );
        int ja=0;
        for(int ia=0; ia<natom_; ia++){
            if(ignoreAtoms_[ia]) continue;
            types [ja] = types_ [ia];
            apos  [ja] = apos_  [ia];
            qrots [ja] = qrots_ [ia];
            aforce[ja] = aforce_[ia];
            torqs [ja] = torqs_ [ia];
            omegas[ja] = omegas_[ia];
            vels  [ja] = vels_  [ia];
            ignoreAtoms[ja] = false;
            fixedAtoms [ja] = fixedAtoms_[ia];
            for(int ib=0; ib<N_BOND_MAX; ib++){
                int i=ia*N_BOND_MAX + ib;
                int j=ja*N_BOND_MAX + ib;
                ebonds  [j] = ebonds_  [i];
                hbonds  [j] = hbonds_  [i];
                fbonds  [j] = fbonds_  [i];
                bondCaps[j] = bondCaps_[i];
            }
            ja++;
        }
        cleanAux();
        natomActive=ja;
        for(int ia=ja; ia<natom; ia++){
            ignoreAtoms[ia] = true;
            fixedAtoms [ia] = false;
        }
        delete[] types_;
        delete[] apos_;
        delete[] qrots_;
        delete[] aforce_;
        delete[] torqs_;
        delete[] omegas_;  // just for MD
        delete[] vels_;    // just for MD
        delete[] ignoreAtoms_;
        delete[] fixedAtoms_;
        delete[] ebonds_;
        delete[] hbonds_;
        delete[] fbonds_;
        delete[] bondCaps_;
    }

    bool tryResize( int nMaskMin=5, int nMaskMax=20, int nMaskGoal=10 ){
        int nmask = 0;
        for(int i=0; i<natom; i++){ if(ignoreAtoms[i]){ nmask++; } }
        if( (nmask<nMaskMin)||(nmask>nMaskMax) ){
            int natom_new = natom-nmask+nMaskGoal;
            printf( "RARFF::resize(%i) from %i nmaxk %i \n", natom_new, natom, nmask );
            resize( natom_new );
            return true;
        }
        return false;
    }

    int inserAtom( RigidAtomType* typ, const int* caps, const Vec3d& p, const Vec3d& dir, const Vec3d& up ){
        Mat3d m;
        //m.c=dir; double r = m.c.normalize();
        m.a=dir; double r = m.a.normalize();
        if(r<1e-3)return -1;
        //m.b.set_cross(up ,m.c); m.b.normalize();
        //m.a.set_cross(m.b,m.c); m.a.normalize();
        m.b=up; m.b.makeOrthoU( m.a ); m.b.normalize();
        m.c.set_cross(m.a,m.b);
        int ia=natomActive;
        printf( "inserAtom ia %i | natom %i \n", ia, natom );
        if( ia>=natom ){ resize(natom+5); }
        natomActive++;
        ignoreAtoms[ia] = false;
        fixedAtoms [ia] = false;
        types[ia] = typ;
        apos [ia] = p;
        qrots[ia].fromMatrixT(m);
        vels [ia]  = Vec3dZero;
        omegas[ia] = Vec3dZero;
        for(int j=0; j<N_BOND_MAX;j++){
            int i = ia*N_BOND_MAX + j;
            bondCaps[i] = caps[j];
            //printf( "caps[%i] %i \n", j, caps[j] );
        }
        return ia;
        map.bResizing=true;
    }

    // ======== Force Evaluation

    inline double pairEF( int ia, int ja, int nbi, int nbj, RigidAtomType& type){

        // subtract atom positions, check if within cutoff
        Vec3d  dij = apos[ja] - apos[ia];
        double r2  = dij.norm2();
        //printf( "pairEF r %g  Rcut %g \n", sqrt(r2), type.Rcut );
        n_pairs_tried++;
        if( r2>R2cut ) return 0;
        n_pairs_evaluated++;
        double rij = sqrt( r2 );
        Vec3d hij  = dij*(1/rij);
        
        // Evaluate radial Potential and Force
        double dexpr;
        double expr = finiteExp( rij - type.rbond0, dexpr, -type.bMorse, RcutMax-type.rbond0 );
        double E     =    type.aMorse*expr* expr;
        double fr    =  2*type.aMorse*expr*dexpr;
        double Eb    = -2*type.aMorse*      expr;
        double frb   = -2*type.aMorse*     dexpr;
        Vec3d force=hij*fr;

        // Find block corresponding to bonds of this atoms in global arrays
        Vec3d* his  = hbonds  +ia*N_BOND_MAX;
        Vec3d* hjs  = hbonds  +ja*N_BOND_MAX;
        Vec3d* fis  = fbonds  +ia*N_BOND_MAX;
        Vec3d* fjs  = fbonds  +ja*N_BOND_MAX;
        double* eis = ebonds  +ia*N_BOND_MAX;
        double* ejs = ebonds  +ja*N_BOND_MAX;
        int* capis  = bondCaps+ia*N_BOND_MAX;
        int* capjs  = bondCaps+ja*N_BOND_MAX;

        SMat3d Dij;  Dij.fromDeriv( dij, 1/(rij*r2) ); // Evaluate angular derivative matrix

        // Interactions between all pairs of bonds on that atom 
        for(int ib=0; ib<nbi; ib++){
            const Vec3d& hi = his[ib];
            Vec3d& fi       = fis[ib];
            double ci       = hij.dot( hi );   // ci = <hi|hij>
            if(ci<0) continue;

            bool capi = ( capis[ib] >= 0 );

            for(int jb=0; jb<nbj; jb++){
                const Vec3d& hj = hjs[jb];
                Vec3d& fj = fjs[jb];
                // Interaction between capped (i.e. passivated) bonds - Only repulsive
                if( bDonorAcceptorCap && capi && (capjs[jb]>=0) ) continue;
                if( bRepelCaps && capi && (capjs[jb]>=0) ){ // repulsion of capping atoms
                    Vec3d pi = apos[ia] + hi*lcap;
                    Vec3d pj = apos[ja] + hj*lcap;
                    Vec3d  dij   = pi-pj;
                    double r2    = dij.norm2() + R2SAFE;
                    if( r2<R2cutCap ){ 
                        // Morse
                        double rij   = sqrt( r2 );
                        //double e     = aMorseCap * exp( bMorseCap*rij );
                        //Vec3d  f     = dij * (-bMorseCap * e / rij);
                        double frc;
                        double e     = aMorseCap * finiteExp( rij, frc, -bMorseCap, RcutCap );   
                        Vec3d  f     = dij * (aMorseCap * frc / rij );     
                        E+=e;
                        force.add(f);
                        f.mul(1.0/lcap);
                        fi.add(f);
                        fj.sub(f);
                        continue;
                    }
                }

                // cosines of angles between the bonds
                double cj       = hij.dot( hj );  // cj  = <hj|hij>
                double cij      = hi .dot( hj );  // cij = <hj|hi>
                if( (cj>0)||(cij>0) ) continue;   // check if bonds oriented toward each other

                // Angular potential     E =  (cos(ai)*cos(aj)*cos(aij)) ^ 4
                double cc  = ci*cj*cij;
                double cc2 = cc*cc;
                double e   = cc2*cc2;   // Angular Energy 
                double de  = 4*cc2*cc;  // Angular derivative

                // Acumulate bond energy
                double eEb = e * Eb;
                eis[ib]+=eEb*0.5;
                ejs[jb]+=eEb*0.5;

                if( isnan(eEb) ){ printf( "ERROR pairEF(ia=%i,ja=%i)(ib=%i,jb=%i): eEb is %i \n", ia,ja, ib,jb, eEb ); exit(0); }

                E      += eEb;

                // Evaluate forces due to bonds
                double deEb  = de*Eb;
                double cicj  = ci*cj;
                double cijcideEb = cij*ci*deEb;
                double cijcjdeEb = cij*cj*deEb;
                double frbe  = frb*e;
                fi.add_mul( hij, -cijcjdeEb ); fi.add_mul( hj, -cicj*deEb );
                fj.add_mul( hij, -cijcideEb ); fj.add_mul( hi, -cicj*deEb );
                Vec3d hv;
                Dij.dot_to( hi, hv ); force.add_mul( hv, cijcjdeEb );
                Dij.dot_to( hj, hv ); force.add_mul( hv, cijcideEb );
                force.add_mul( hij, frb*e );
            }
        }

        // pi-pi interactions ( trying to align pi-orbitals in aromatic molecules, planarize them )
        //printf( "%i,%i nbi %i nbj %i \n", ia, ja, nbi, nbj );
        if( (nbi==3)&&(nbj==3) ){ // align pz-pz in sp2-sp2 pair
            const Vec3d& hi = his[3];
            const Vec3d& hj = hjs[3];
            Vec3d&       fi = fis[3];
            Vec3d&       fj = fjs[3];
            double cdot = hi.dot(hj);
            //double E    = Epz * ( cdot*cdot );
            double de = -2*cdot*type.Epz*Eb;
            //printf( "cdot %g Epz %g de %g \n", cdot, type.Epz, de );
            fi.add_mul(hj,de);
            fj.add_mul(hi,de);
            //force      // TODO: proper derivatives of energy
        }
        // Accumulate atomic forces
        aforce[ja].sub(force);
        aforce[ia].add(force);
        return E;
    }

    double interEF_brute(){   // Naieve algorithm for pairwise interactions between atoms O(n2)
        //Vec3d bhs[N_BOND_MAX];
        double E = 0;
        for(int i=0; i<natom; i++){
            if(ignoreAtoms[i])continue;
            const RigidAtomType& typei = *types[i];
            Vec3d           pi   = apos[i];
            for(int j=i+1; j<natom; j++){
            //for(int j=0; j<natom; j++){
                if(ignoreAtoms[j])continue;
                RigidAtomType pairType;
                pairType.combine( typei, *types[j] );
                E += pairEF( i, j, typei.nbond, types[j]->nbond, pairType );
            }
        }
        return E;
    }

    double interEF_list(){   // Naieve algorithm for pairwise interactions between atoms O(n2)
        //printf( "DEBUG interEF_list() npairs %i pairs %li \n", npairs, (long)pairs );
        double E = 0;
        RigidAtomType pairType;
        for(int i=0; i<npairs; i++){
            //printf( "pair[%i] \n" );
            Vec2i ij = pairs[i];
            //printf( "pair[%i]( %i,%i) \n", i, ij.i, ij.j );
            pairType.combine( *types[ij.i], *types[ij.j] );
            E += pairEF( ij.i, ij.j, types[ij.i]->nbond, types[ij.j]->nbond, pairType );
        }
        //printf( "DEBUG interEF_list() DONE \n" );
        return E;
    }

    double interEF_buckets(){  // Grid-accelerated algorithm for pairwise interactions between atoms O(n), but rather large prefactor
        map.updateNeighsBufferSize(); // make sure neighs has sufficient size
        int* neighs = map.neighs_in;
        double E = 0;
        Vec3i ip;
        RigidAtomType pairType;
        for(int ic=0; ic<map.ncell; ic++){           // go over all cells
            int nic = map.getInCell( ic, neighs );   // list atoms in same cell
            if( nic>0){                              // if cell is not empty
                map.i2ixyz(ic,ip);                   // unfold 1d-cell index to 3d-cell index
                int* neighs_ = neighs+nic;           // rewind neighs to free slots
                int nrest    = map.getForwardNeighbors( ip, neighs_ );    // list atoms in neighboring cells
                for(int i=0; i<nic; i++){            // for all atoms in pivot cell
                    int ia = neighs[i];
                    const RigidAtomType& typei = *types[ia];
                    Vec3d                pi    =  apos [ia];
                    // -- within same cell
                    //printf( "Same Cell\n");
                    for(int j=i+1; j<nic; j++){    // all other atoms in the same cell (prevent double-counting)  
                        int ja = neighs[j];
                        pairType.combine( typei, *types[ja] );                           // evaluate interaction parameters form the two atom types 
                        E += pairEF( ia, ja, typei.nbond, types[ja]->nbond, pairType );  // Interact
                        //glColor3f(0,0,0); Draw3D::drawLine( apos[ia], apos[ja] );
                    }
                    //printf( "Other Cells\n");
                    // -- with neighbor cells
                    for(int j=0; j<nrest; j++){   // all atoms in the other cells
                        int ja = neighs_[j];
                        if( ia==ja ) continue; // this should not happen !!!! ( but for some reason it is happening )
                        pairType.combine( typei, *types[ja] );                           // evaluate interaction parameters form the two atom types 
                        E += pairEF( ia, ja, typei.nbond, types[ja]->nbond, pairType );  // Interact
                        //if( ic==31 ){ glColor3f(0,0,0); Draw3D::drawLine( apos[ia], apos[ja] ); }
                    }
                }
            }
        }
        return E;
    }

    void cleanAux(){
        for(int i=0; i<natom; i++){
            vels  [i]=Vec3dZero;
            omegas[i]=Vec3dZero;
            torqs [i]=Vec3dZero;
            aforce[i]=Vec3dZero;
        }
        int nb   = natom*N_BOND_MAX;
        for(int i=0; i<nb;   i++){
            ebonds  [i]=0;
            //bondCaps[i]=-1;
            fbonds  [i].set(0.0);
        }
    }

    void cleanAtomForce(){
        for(int i=0; i<natom; i++){
            //atoms[i].cleanForceTorq();
            torqs [i]=Vec3dZero;
            aforce[i]=Vec3dZero;
        }
        int nb   = natom*N_BOND_MAX;    for(int i=0; i<nb;   i++){ ebonds[i]=0; }
        //printf("\n"); for(int i=0; i<nb;   i++){ printf("%i ", bondCaps[i] ); }; printf("\n");
        int nval = natom*N_BOND_MAX*3;  for(int i=0; i<nval; i++){ ((double*)fbonds)[i]=0;}
    }

    double evalF2rot()const{ double F2=0; for(int i=0; i<natom; i++){ if(ignoreAtoms[i])continue; F2+=torqs [i].norm2(); }; return F2; }
    double evalF2pos()const{ double F2=0; for(int i=0; i<natom; i++){ if(ignoreAtoms[i])continue; F2+=aforce[i].norm2(); }; return F2; }

    inline void projectAtomBons(int ia){
        qrots[ia].rotateVectors( N_BOND_MAX, types[ia]->bh0s, hbonds+ia*N_BOND_MAX, false );
    }

    void projectBonds(){
        for(int i=0; i<natom; i++){ if(ignoreAtoms[i])continue; projectAtomBons(i); }
    }

    inline Vec3d bondPosBare( int i, double r )const{
        int ia = i/N_BOND_MAX;
        return apos[ia]  + hbonds[i]*r;
    }

    inline Vec3d bondPos( int i, double sc=1.0 )const{
        //int i = ia*N_BOND_MAX + j;
        int ia = i/N_BOND_MAX;
        double R = types[ia]->rbond0;
        return apos[ia]  + hbonds[i]*(R*sc);
    }

    void applyForceHarmonic1D(const Vec3d& h, double x0, double K){
        //printf( "applyForceHarmonic1D %g %g (%g,%g,%g) \n", x0, K, h.x,h.y,h.z  );
        for(int ia=0; ia<natom; ia++){
            if(ignoreAtoms[ia])continue;
            double x = h.dot(apos[ia])-x0;
            aforce[ia].add_mul( h, K*x );
        }
    }

    void applyForceBox(const Vec3d& p0, const Vec3d& p1, double K, double fmax){
        //printf( "applyForceHarmonic1D %g %g (%g,%g,%g) (%g,%g,%g) \n", K, fmax, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z  );
        double f2max = fmax*fmax;
        for(int ia=0; ia<natom; ia++){
            if(ignoreAtoms[ia])continue;
            Vec3d p = apos[ia];
            Vec3d f;
            Vec3d d0,d1;
            d0=p-p0;
            d1=p-p1;
            if( d0.x<0 ){ f.x=K*d0.x; }else if( d1.x>0 ){ f.x=K*d1.x; };
            if( d0.y<0 ){ f.y=K*d0.y; }else if( d1.y>0 ){ f.y=K*d1.y; };
            if( d0.z<0 ){ f.z=K*d0.z; }else if( d1.z>0 ){ f.z=K*d1.z; };
            double f2 = f.norm2();
            if( f2>f2max ){
                f.mul( sqrt(f2max/f2) );
            }
            //printf( "ia %i (%g,%g,%g) (%g,%g,%g) \n", ia, p.x,p.y,p.z,  f.x,f.y,f.z );
            aforce[ia].add(f);
        }
    }

    void evalTorques(){
        for(int ia=0; ia<natom; ia++){
            if(ignoreAtoms[ia])continue;
            //RigidAtom& atomi = atoms[ia];
            //Vec3d torq = Vec3dZero;
            //int nbi =  types[ia]->nbond;
            //for(int ib=0; ib<nbi; ib++){
            for(int ib=0; ib<N_BOND_MAX; ib++){
                int io = 4*ia+ib;
                fbonds[io].makeOrthoU(hbonds[io]);
                //printf( "ia %i ib %i f(%g,%g,%g)\n", ia, ib,  fbonds[io].x,fbonds[io].y,fbonds[io].z );
                torqs[ia].add_cross( hbonds[io], fbonds[io] );
            }
        }
    }

    void eval(){
        cleanAtomForce();
        projectBonds();     
        if(AccelType==1){
            map.pointsToCells( natomActive, apos, ignoreAtoms );
            interEF_buckets();   
        }else if(AccelType==2){
            interEF_list();     
        }else{
            interEF_brute(); 
        }
    }

    void move(double dt){
        evalTorques();
        for(int i=0; i<natom; i++){
            if( ignoreAtoms[i] || fixedAtoms[i] )continue;
            //atoms[i].moveRotGD(dt*invRotMass);
            //atoms[i].movePosGD(dt);
            qrots[i].dRot_exact( dt, torqs[i] );  qrots[i].normalize();          // costly => debug
            //dRot_taylor2( dt, torqs[i] );   qrots[i].normalize_taylor3();  // fast => op
            apos[i].add_mul(aforce[i],dt);
        }
    }

    void moveMDdamp(double dt, double damp){
        evalTorques();
        for(int i=0; i<natom; i++){
            //if( fixedAtoms[i] ){ printf( "fixedAtoms[%i] \n", i ); continue; }
            if(ignoreAtoms[i] || fixedAtoms[i] )continue;
            //atoms[i].moveMDdamp( dt, invRotMass, damp);
            vels  [i].set_lincomb( damp, vels  [i], dt,            aforce[i] );
            omegas[i].set_lincomb( damp, omegas[i], dt*invRotMass, torqs [i] );
            apos  [i].add_mul    ( vels[i], dt    );
            qrots [i].dRot_exact ( dt,  omegas[i] );
            qrots [i].normalize();
        }
    }

    int passivateBonds( double Ecut ){
        printf( " c6cap, c12cap %g %g \n", c6cap, c12cap );
        for(int i=0; i<natom*N_BOND_MAX; i++){ bondCaps[i]=-1; };
        int n=0;
        for(int ia=0; ia<natom; ia++){
            if(ignoreAtoms[ia])continue;
            int nbi =  types[ia]->nbond;
            for(int ib=0; ib<nbi; ib++){
                int i = ia*N_BOND_MAX + ib;
                if(ebonds[i]>Ecut){
                    bondCaps[i]=1;
                    n++;
                };
            }
        }
        return n;
    }

    int saveXYZ(const char* fname )const{
        printf( "RARFFarr::saveXYZ(%s) \n", fname );
        FILE * pFile = fopen(fname,"w");
        int na = 0;
        for(int i=0; i<natom; i++){ if(!ignoreAtoms[i])na++; }
        fprintf( pFile, "        \n" );
        fprintf( pFile, "#comment \n" );
        int n = 0;
        for(int ia=0; ia<natom; ia++){
            if(ignoreAtoms[ia])continue;
            int i0   = ia*N_BOND_MAX;
            fprintf( pFile, "%s %3.6f %3.6f %3.6f \n", types[ia]->name, apos[ia].x,apos[ia].y,apos[ia].z ); n++;
            for(int j=0; j<types[ia]->nbond; j++){
                int icap = bondCaps[i0+j];
                if(icap>=0){
                    const CapType& ct = capTypeList[icap];
                    Vec3d p = apos[ia] + hbonds[i0+j]*ct.rbond0;
                    fprintf( pFile, "%s %3.6f %3.6f %3.6f \n", ct.name, p.x,p.y,p.z ); n++;
                }
            }
        }
        fseek(pFile, 0L, SEEK_SET);
        fprintf( pFile, "%4i ", n ); // write number of atoms at the beggining of file
        fclose(pFile);
        return natom;
    }

    int save(const char* fname )const{
        printf( "RARFFarr::save(%s) \n", fname );
        FILE * pFile = fopen(fname,"w");
        int na = 0;
        for(int i=0; i<natom; i++){ if(!ignoreAtoms[i])na++; }
        fprintf( pFile, "%i %i \n", na, N_BOND_MAX );
        for(int ia=0; ia<natom; ia++){
            if(ignoreAtoms[ia])continue;
            int i0 = ia*N_BOND_MAX;
            int nret = fprintf( pFile, "%i   %3.6f %3.6f %3.6f    %3.6f %3.6f %3.6f %3.6f   %i %i %i %i  \n", types[ia]->id,
                    apos [ia].x,apos [ia].y,apos [ia].z,
                    qrots[ia].x,qrots[ia].y,qrots[ia].z,qrots[ia].w,
                    bondCaps[i0+0],bondCaps[i0+1],bondCaps[i0+2],bondCaps[i0+3]
                );
        }
        fclose(pFile);
        return natom;
    }

    int load(const char* fname ){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){ printf("ERROR RARFFarr::load() cannot find %s\n", fname ); return -1; }
        char buff[1024];
        char * line;
        int nl,nba,na;
        line = fgets( buff, 1024, pFile ); printf("%s",line);
        sscanf( line, "%i %i \n", &na, &nba ); printf( "na %i nba %i \n", na, nba );
        if(nba!=N_BOND_MAX){ printf("ERROR RARFFarr::load() nba(%i)!=N_BOND_MAX(%i)\n", nba, N_BOND_MAX ); fclose(pFile); return -1; }
        //allocate(natoms,nbonds);
        if(natom<na)resize( na );
        //resize( na );
        //line = fgets( buff, 1024, pFile ); // comment
        int ityp;
        for(int ia=0; ia<na; ia++){
            int i0 = ia*N_BOND_MAX;
            line = fgets( buff, 1024, pFile );  printf("%s",line);
            Vec3d  p;
            Quat4d q;
            Quat4i cp;
            int nret = sscanf( line, "%i   %lf %lf %lf    %lf %lf %lf %lf   %i %i %i %i", &ityp, &p.x,&p.y,&p.z, &q.x,&q.y,&q.z,&q.w, &cp.x,&cp.y,&cp.z,&cp.w    );
            apos[ia]=p; qrots[ia]=q;  (*((Quat4i*)(bondCaps+i0)))=cp;
            types[ia] = &typeList[ityp];
            printf( "atom[%i] %i   p(%g,%g,%g)    qrot(%g,%g,%g)  caps(%i,%i,%i,%i)\n", ia, ityp,
                apos [ia].x,apos [ia].y,apos [ia].z,
                qrots[ia].x,qrots[ia].y,qrots[ia].z,
                bondCaps[i0+0],bondCaps[i0+1],bondCaps[i0+2],bondCaps[i0+3] );
            ignoreAtoms[ia]=false;
        }
        if(natomActive<na)natomActive=na;
        cleanAux();
        projectBonds();
        fclose(pFile);
        return natom;
    }


};

#endif
