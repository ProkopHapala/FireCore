
#ifndef FF2F_h
#define FF2F_h

#include <vector>

#include "testUtils.h"

#include "datatypes.h"
#include "Vec2.h"

//              // H He Li Be B C N O F 
static const int atom_ne[9] = { 0,0, 0, 0, 0,0,1,2,3};

class Atom2D{ public:
    int typ;
    Vec2d pos;
    int nb,ne,npi;
    int neighs[4];

    double getCos0(){
        if       (npi==2){return -1.0;}
        //else if(npi==1){return -0.5;}
        return -0.5;
    }

    void cleanNeighs(){ for(int i=0;i<4;i++){neighs[i]=-1;}  };
    bool addNeigh(int i){ if((nb+ne+npi)>3)return false; neighs[nb]=i; nb++; return true; };

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
        printf("nb,ne,npi(%i,%i,%i) neighs{%i,%i,%i,%i} pos(%g,%g)\n", nb,ne,npi, neighs[0],neighs[1],neighs[2],neighs[3], pos.x,pos.y );
    }

};

class FF2D{ public:
    std::vector<Atom2D> atoms;
    //std::vector<Vec2i>  bonds;

    int old_size =0;
    //Vec2d* apos  =0;
    Vec2d* vels  =0;
    Vec2d* forces=0;

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
                    double Eijk =  Kang*c_*c_;
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
        return Ebond+Eang;
    }

    double moveMD(double dt, double damping){
        double cdamp = 1-damping;
        double F2sum=0;
        for(int i=0; i<atoms.size();i++){
            Vec2d        v = vels[i];
            const Vec2d& f = forces[i];
            v.mul(cdamp);
            F2sum += f.norm2();
            v.add_mul( f, dt );
            atoms[i].pos.add_mul(v, dt);
            vels[i]=v;
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

    int run( int n, double dt, double damp, double F2conv, bool bCleanForce){
        try_realloc();
        if(bCleanForce)cleanVelocity();
        for(int itr=0; itr<n; itr++){
            double F2sum = step( dt, damp);
            //printf( "run[%i] E %g |F| %g \n", itr, Etot, sqrt(F2sum) );
            if(F2sum<F2conv)return itr;
        }
    }

    int insertString(const char* str, char sep=';'){
        constexpr const int ntmp=256; 
        char tmp[ntmp];
        DEBUG
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
