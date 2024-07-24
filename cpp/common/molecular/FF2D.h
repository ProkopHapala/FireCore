
#ifndef FF2F_h
#define FF2F_h

#include <vector>

#include "testUtils.h"

#include "datatypes.h"
#include "Vec2.h"
#include "geom2D.h"

//              // H He Li Be B C N O F 
static const int atom_ne[9] = { 0,0, 0, 0, 0,0,1,2,3};

class Atom2D{ public:
    int typ=-1;
    Vec2d pos;
    int nb=0,ne=0,npi=0;
    int neighs[4];
    //int remove = 0;

    double getCos0(){
        if       (npi==2){return -1.0;}
        //else if(npi==1){return -0.5;}
        return -0.5;
    }

    void cleanNeighs(){ 
        //printf("Atom2D::cleanNeighs()\n");
        for(int i=0;i<4;i++){neighs[i]=-1;}  };
    //bool addNeigh(int i){ if((nb+ne+npi)>3)return false; neighs[nb]=i; nb++; return true; };

    bool addNeigh(int ja){
        for(int i=0;i<4;i++){ if(neighs[i]<0){ neighs[i]=ja; return true; } }
        return false;
    };

    void repairNeighs(){
        int imax   = 0;
        int nblank = 0;
        int nvoids = 0;
        for(int i=0;i<4;i++){ 
            if(neighs[i]>=0){ imax=i; nvoids=nblank; }else{ nblank++; }
        }
        nb = 4-nblank;
        if(nvoids>0){
            for(int i=0;i<imax;i++){ 
                if(neighs[i]<0){ 
                    _swap(neighs[i],neighs[imax]);
                    //neighs[i]=neighs[imax];
                    //neighs[imax]=-1;
                    while( (imax==0) || (neighs[imax]>=0) ) imax--; 
                }
                if(i>=imax){ break; } 
            }
        }
    }

    void changeNeigh(int ia, int ja){ 
        for(int i=0;i<4;i++){ if(neighs[i]==ia){ neighs[i]=ja; } }
        if(ja<0){ repairNeighs(); }
    }

    void fromString(const char* s ){
        //printf("DEBUG fromString() %s \n", s );
        int n = sscanf( s, "%i %i    %i %i %i %i   %lf %lf ", &typ, &npi, neighs, neighs+1, neighs+2, neighs+3, &(pos.x), &(pos.y) );
        nb=n-2; if(nb>4)nb=4;
        //printf("DEBUG fromString() n=%i %s \n", n, s );
        //for(int i=nb; i<4; i++){ neighs[i]=-1; }
        for(int i=0; i<4; i++){ if(i<nb){neighs[i]--;}else{neighs[i]=-1;} }
        ne=atom_ne[typ-1];
        //print();
    }

    void print(){
        //printf("nb,ne,npi(%i,%i,%i) neighs{%i,%i,%i,%i} pos(%g,%g)\n", nb,ne,npi, neighs[0],neighs[1],neighs[2],neighs[3], pos.x,pos.y );
        printf("neighs{%i,%i,%i,%i} pos(%g,%g)\n", neighs[0],neighs[1],neighs[2],neighs[3], pos.x,pos.y );
    }

};

class FF2D{ public:
    std::vector<Atom2D> atoms;
    std::vector<Vec2i>  bonds;

    int old_size =0;
    //Vec2d* apos  =0;
    Vec2d* vels  =0;
    Vec2d* forces=0;
    bool bRelax = true;

    double Kang=1.,Kbond=1.,L0=1.;
    double Etot=0.,Ebond=0.,Eang=0.;

    void realloc(int n){
        _realloc(vels,n);
        _realloc(forces,n);
        old_size=n;
    }

    void try_realloc(int nnew=-1){
        if(nnew<0) nnew = atoms.size();
        if(old_size==nnew) return;
        realloc(nnew);
        cleanForce   ();
        cleanVelocity();
    }

    void toArrays( int* types, Vec2d* apos, int4* neighs ){
        for(int ia=0; ia<atoms.size(); ia++ ){
            Atom2D& A = atoms[ia];
            if(types )types [ia]=A.typ;
            if(neighs)neighs[ia]=*(int4*)A.neighs;
            if(apos  )apos  [ia]= A.pos;
        }
    }

    void getBonds( Vec2i* bs ){
        for(int ib=0; ib<bonds.size(); ib++){ bs[ib] = bonds[ib]; }
    };

    void repairAtoms(){
        for(int ia=0; ia<atoms.size(); ia++){ atoms[ia].repairNeighs(); }
    }

    int addAtom( Vec2d pos, int type, int* neighs=0 ){
        //printf("addAtom(t=%i,pos(%g,%g)) @neighs=%li \n", type, pos.x, pos.y, (long)neighs );
        Atom2D A;
        A.pos=pos;
        A.typ=type;
        if( neighs ){ for(int i=0; i<4; i++){ A.neighs[i] = neighs[i]; } }else{ A.cleanNeighs(); }
        //printf( "addAtom[%i]", atoms.size() ); A.print();
        atoms.push_back( A );
        return atoms.size()-1;
    }

    int addBond( int ia, int ja ){
        //printf("addBond(%i,%i) \n", ia, ja );
        atoms[ia].addNeigh(ja);
        atoms[ja].addNeigh(ia);
        bonds.push_back( Vec2i{ia,ja} );
        return bonds.size()-1;
    }

    void swapAtoms( int i, int j ){
        Atom2D ai = atoms[i];
        Atom2D aj = atoms[j];
        // re-index neighbors of neighbors
        for(int k=0; k<4; k++){
           int ia = ai.neighs[k]; atoms[ia].changeNeigh(i,j);
           int ja = aj.neighs[k]; atoms[ja].changeNeigh(j,i);
        }
        Atom2D tmp = atoms[i]; atoms[i] = atoms[j]; atoms[j] = tmp; 
    }

    void removeAtomSwap( int i ){
        //printf("removeAtomSwap(%i) \n", i );
        int j = atoms.size()-1;
        Atom2D ai = atoms[i];
        Atom2D aj = atoms[j];
        // re-index neighbors of neighbors
        for(int k=0; k<4; k++){
           int ia = ai.neighs[k]; if(ia>=0){ atoms[ia].changeNeigh(i,-1); }
           int ja = aj.neighs[k]; if(ja>=0){ atoms[ja].changeNeigh(j,i); }
        }
        atoms[i] = atoms[j];    
        atoms.pop_back();
    }

    int bondsFromNeighs( ){
        //printf("bondsFromNeighs() \n" );
        bonds.clear();
        for(int ia=0; ia<atoms.size(); ia++){
            Atom2D& a = atoms[ia];
            //printf("atom[%i] ngs{%i,%i,%i,%i} \n", ia, a.neighs[0],a.neighs[1], a.neighs[2], a.neighs[3] );
            for(int j=0; j<4; j++){
                int ja = a.neighs[j];
                if( (ja>=0) && (ja<ia) ){
                    bonds.push_back( Vec2i{ia,ja} );
                }
            }
        }
        return bonds.size();
    }

    bool removeAtom(int i){ if( (i>0)&&(i<atoms.size() ) ){ atoms[i].typ=-1;              return true; }else{ return false; } }
    bool removeBond(int i){ if( (i>0)&&(i<bonds.size() ) ){ bonds[i].x=-1; bonds[i].y=-1; return true; }else{ return false; } }

    int findBondAt( Vec2d pos, double R ){
        double R2max = R*R;
        double r2min = 1e+300;
        int ifound = -1;
        for(int ib=0; ib<bonds.size(); ib++){
            Vec2i& b = bonds[ib];
            if( (b.x<0) || (b.y<0) )continue;
            Vec2d d = dpLineSegment( pos, atoms[b.x].pos, atoms[b.y].pos );
            double r2 = d.norm2();
            if( (r2 < R2max) && (r2 < r2min) ){
                ifound = ib;
                r2min = r2;
            }
        }
        return ifound;
    } 

    int findAtomAt( Vec2d pos, double R ){
        //printf("findAtomAt(pos=(%g,%g),R=%g) \n", pos.x, pos.y, R );
        double R2max = R*R;
        double r2min = 1e+300;
        int ifound = -1;
        for(int ia=0; ia<atoms.size(); ia++){
            if(atoms[ia].typ<0) continue;
            Vec2d d  = pos-atoms[ia].pos;
            double r2 = d.norm2();
            if( (r2 < R2max) && (r2 < r2min) ){
                ifound = ia;
                r2min = r2;
            }
            //printf("atom[%i] pos(%g,%g) r=%g ifound=%i \n", ia, atoms[ia].pos.x, atoms[ia].pos.y, sqrt(r2), ifound );
        }
        return ifound;
    }

    void print_atoms(){
        printf("FF2D::print_atoms() n=%i \n", atoms.size() );
        for(int ia=0; ia<atoms.size(); ia++){
            printf("[%i]", ia); atoms[ia].print();
            //Atom2D& A = atoms[ia];
            //printf("atom[%i] pos(%g,%g) typ=%i nb=%i neighs{%i,%i,%i,%i} \n", ia, A.pos.x, A.pos.y, A.typ, A.nb, A.neighs[0],A.neighs[1],A.neighs[2],A.neighs[3] );
        }
    }

    double eval(){
        Vec2d   hbs  [4]; 
        double  invls[4];
        Ebond=0;
        Eang =0;
        for(int ia=0; ia<atoms.size(); ia++){
            Atom2D& A = atoms[ia];
            double c0  = A.getCos0();  //printf("atom[%i] c0 %g \n", ia, c0);
            // ------ Bond pass
            for(int j=0; j<A.nb; j++){
                int ja = A.neighs[j];
                Atom2D& B = atoms[ja];
                Vec2d d = A.pos-B.pos;
                double l=d.normalize();
                hbs  [j]=d;
                invls[j]=1./l;
                if(ja>ia){ // double counting prevention
                    //printf( "bond[%i,%i] l=%g\n", ia,ja, l );
                    double dl  = L0-l;
                    double fr  = Kbond*dl;
                    double Eij = fr*0.5*dl;
                    Vec2d f    = d * fr;  //springForce2D( A.pos-B.pos, 1., -1. );
                    forces[ia].add(f);
                    forces[ja].sub(f);
                    Ebond+=Eij;
                }
            }
            // ----- Angle Pass
            for(int j=0; j<A.nb; j++){
                int ja  = A.neighs[j];
                double         il1= invls[j];
                const   Vec2d& h1 = hbs  [j];
                for(int k=j+1; k<A.nb; k++){
                    int ka  = A.neighs[k];
                    double        il2 = invls[k];
                    const  Vec2d& h2  = hbs  [k];
                    // ---- copied from  MMFFsp3::evalSigmaSigma_cos()
                    double c = h1.dot( h2 );
                    Vec2d hf1,hf2; // unitary vectors of force - perpendicular to bonds
                    hf1 = h2 - h1*c;
                    hf2 = h1 - h2*c;
                        // Draw2D::drawVecInPos( (Vec2f)hf1, (Vec2f)atoms[ja].pos );
                        // Draw2D::drawVecInPos( (Vec2f)hf2, (Vec2f)atoms[ka].pos );
                    double c_ = c-c0;
                    //printf("angle[%i|%i,%i] c %g c0 %g \n", ia,ja,ka, c, c0);
                    double Eijk = Kang*c_*c_;
                    double fang = Kang*c_*2;
                    hf1.mul( fang*il1 );
                    hf2.mul( fang*il2 );
                        // Draw2D::drawVecInPos( (Vec2f)hf1, (Vec2f)atoms[ja].pos );
                        // Draw2D::drawVecInPos( (Vec2f)hf2, (Vec2f)atoms[ka].pos );
                    forces[ja].add( hf1     );
                    forces[ka].add( hf2     );
                    forces[ia].sub( hf1+hf2 );
                    Eang+=Eijk;
                }
            }
        }
        //printf("eval() Ebond=%g Eang=%g \n", Ebond, Eang );
        return Ebond+Eang;
    }

    double moveMD(double dt, double damping){
        double cdamp = 1-damping;
        double F2sum = 0;
        double cfv   = 0;
        for(int i=0; i<atoms.size();i++){
            Vec2d        v = vels[i];
            const Vec2d& f = forces[i];
            F2sum += f.norm2();
            cfv   += f.dot(v);
            v.mul(cdamp);
            v.add_mul( f, dt );
            atoms[i].pos.add_mul(v, dt);
            vels[i]=v;
        } 
        if( bRelax && (cfv<0) ){
            for(int i=0; i<atoms.size();i++){ vels[i].set(0.); };
        }
        return F2sum;
    }

    void cleanForce   (){ for(int i=0; i<atoms.size();i++){forces[i].set(0.);}; };
    void cleanVelocity(){ for(int i=0; i<atoms.size();i++){vels  [i].set(0.);}; };
    void setPosRandom (double sz){ for(int i=0; i<atoms.size();i++){ atoms[i].pos.set(randf(-sz,sz),randf(-sz,sz)); }; };

    double step( double dt, double damp){
        cleanForce();
        Etot=eval();
        return moveMD( dt, damp);
    }

    int run( int niter, double dt, double damp, double Fconv, bool bCleanForce){
        printf("FF2D::run() n=%i dt=%g damp=%g Fconv=%g \n", niter, dt, damp , Fconv );
        double F2conv=Fconv*Fconv;
        try_realloc();
        repairAtoms();
        //print_atoms();
        if(bCleanForce)cleanVelocity();
        int itr=0;
        for(itr=0; itr<niter; itr++){
            double F2sum = step( dt, damp);
            //if(verbosity>1)
            printf( "run[%i] E %g |F| %g \n", itr, Etot, sqrt(F2sum) );
            if(F2sum<F2conv)return itr;
        }
        return itr;
    }

    int insertString(const char* str, char sep=';'){
        constexpr const int ntmp=256; 
        char tmp[ntmp];
        const char* pc   = str;
        //const char* pc_a0=pc;
        Atom2D atom;
        int na=0;
        int i0=atoms.size();
        int neigh=-1;
        int itmp=0;
        while(true){
            char  c  = *pc;  //printf("c[%li]: %c  \n", pc-str, c );
            tmp[itmp]=c;
            if( c==sep ){
                //printf("SEP\n");
                tmp[itmp]='\0';
                atom.fromString(tmp);
                atoms.push_back( atom );
                atom.cleanNeighs();
                itmp=0;
                neigh=-1;
                na++;
            }else{
                itmp++;
            }
            pc++;
            if(itmp>=ntmp){ printf("ERROR: FF2D::insertString() char tmp[%i] overflow \n", ntmp ); exit(0); }
            if(c=='\0'){ break; } // null terminated string
        }
        return na;
    }
    
}; // MMFF

#endif
