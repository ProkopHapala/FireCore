#ifndef MMFFBuilder_h
#define MMFFBuilder_h

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>

#include "macroUtils.h"
#include "testUtils.h"

//#include "Molecule.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "raytrace.h"

#include "molecular_utils.h"

//#include "Molecule.h"
//#include "MMFF.h"
//#include "MMFFmini.h"
#include "MMFFparams.h"

// =============== Structs for Atoms, Bonds etc...

namespace MM{

static const double const_eVA2_Nm = 16.0217662;

// ============================
// ========= Atom
// ============================
struct Atom{
    constexpr const static Vec3d HcapREQ    = Vec3d{ 1.4870, sqrt(0.000681    ), 0 };
    constexpr const static Vec3d defaultREQ = Vec3d{ 1.7,    sqrt(0.0037292524), 0 };

    // this breaks {<brace-enclosed initializer list>} in C++11
    //int type  = -1;
    //int frag  = -1;
    //int iconf = -1;
    //Vec3d pos;
    //Vec3d REQ = defaultREQ;   // constexpr Vec3d{1.7,sqrt(0.0037292524),0}

    int type;
    int frag;
    int iconf;
    Vec3d pos;
    Vec3d REQ;   // constexpr Vec3d{1.7,sqrt(0.0037292524),0}

    //Atom() = default;

    void print()const{ printf( " Atom{ t %i c %i f %i REQ(%g,%g,%g) pos(%g,%g,%g)}", type, iconf, frag, REQ.x, REQ.y, REQ.z, pos.x,pos.y,pos.z ); }

    Atom() = default;
    Atom(const Vec3d& pos_):type(0),frag(-1),iconf(-1),REQ(defaultREQ),pos(pos_){};
    Atom(const Vec3d& pos_,const Vec3d& REQ_):type(0),frag(-1),iconf(-1),REQ(REQ_),pos(pos_){};
    Atom(int type_,int frag_,int iconf_,const Vec3d& pos_,const Vec3d& REQ_):type(type_),frag(frag_),iconf(iconf_),REQ(REQ_),pos(pos_){};
};

#define N_NEIGH_MAX 4
enum class NeighType: int {
    pi    = -2,
    epair = -3,
    H     = -4
};

// ============================
// ========= AtomConf
// ============================

struct AtomConf{

    int iatom=-1;
    uint8_t n     =0;
    uint8_t nbond =0;
    uint8_t npi   =0; // pi bonds
    uint8_t ne    =0; // electron pairs
    uint8_t nH    =0; //
    int neighs[N_NEIGH_MAX]; // neighs  - NOTE: bonds not atoms !!!!

    //AtomConf() = default;

    inline int fillIn_npi_sp3(){ npi=4-nbond-ne-nH; n=4; return npi; };
    inline int clean_pi_sp3  (){ npi=0; n=nbond+ne+nH;   return n;   };

    inline void rebase(int dib){
        for(int i=0; i<N_NEIGH_MAX; i++){ if(neighs[i]>=0)neighs[i]+=dib; };
    }

    inline int findNeigh(int ib)const{
        for(int i=0; i<N_NEIGH_MAX; i++){
            //printf( "MM:AtomConf.findNeigh()[%i] ng(%i)?=ib(%i) %i \n", i, neighs[i], ib, neighs[i]==ib );
            if(neighs[i]==ib) return i;
        }
        return -1;
    }

    inline bool addNeigh(int ib, uint8_t& ninc ){
        if(n>=N_NEIGH_MAX)return false;
        if(ib>=0){ neighs[nbond]=ib; }else{ neighs[N_NEIGH_MAX-(n-nbond)-1]=ib; };
        ninc++;
        n++;
        //printf( "bond.addNeigh n==%i ninc==%i\n", n, ninc );
        return true;
    };

    inline bool replaceNeigh(int ib, int jb){
        for(int i=0; i<N_NEIGH_MAX; i++){
            if(neighs[i]==ib){ neighs[i]=jb; return true; };
        }
        return false;
    }


    inline bool addBond (int i){ return addNeigh(i,nbond); };
    inline bool addH    (     ){ return addNeigh((int)NeighType::H    ,nH ); };
    inline bool addPi   (     ){ return addNeigh((int)NeighType::pi   ,npi); };
    inline bool addEpair(     ){ return addNeigh((int)NeighType::epair,ne ); };

    inline void clearNonBond(){ n=nbond; npi=0;ne=0;nH=0; };
    inline void clearBond   (){ nbond=0; n=npi+ne+nH;     };
    inline void setNonBond(int npi_,int ne_){ npi=npi_; ne=ne_; n=nbond+npi+ne+nH;  }
    inline void init0(){ for(int i=0; i<N_NEIGH_MAX; i++)neighs[i]=-1; nbond=0; clearNonBond(); }

    //void print()const{ printf( " AtomConf{ ia %i, n %i nb %i np %i ne %i nH %i (%i,%i,%i,%i) }", iatom, n, nbond, npi, ne, nH , neighs[0],neighs[1],neighs[2],neighs[3] ); }
    void print()const{ printf( " AtomConf{ia %i, n,nb,np,ne,nH(%i,%i,%i,%i,%i) [%i,%i,%i,%i]}", iatom, n,nbond,npi,ne,nH, neighs[0],neighs[1],neighs[2],neighs[3] ); }


    inline bool checkNeighsRepeat()const{
        for(int i=0; i<N_NEIGH_MAX; i++){
            const int ib=neighs[i];
            if(ib<0)continue; 
            for(int j=i+1; j<N_NEIGH_MAX; j++){
                //printf( "ib,jb %i %i \n", ib,neighs[j]);
                if(ib==neighs[j]) return true;
            }
        }
        return false;
    };


    AtomConf() = default;
    //AtomConf(int iatom_,int npi_)       :iatom(iatom_),npi(npi_),ne(0  ),nH(0),nbond(0),n(npi_    ){};
    AtomConf(int iatom_,int npi_,int ne_):iatom{iatom_},npi{npi_},ne{ne_},nH{0},nbond{0},n{npi_+ne_}{ for(int i=0;i<N_NEIGH_MAX;i++)neighs[i]=-1; };
    //AtomConf(const MMFFAtomConf&) = default;
    //AtomConf(std::initializer_list<MMFFAtomConf>) {};
};

// ============================
// ========= Bond
// ============================
struct Bond{
    // --- this breaks {<brace-enclosed initializer list>} in C++11
    int    type  = -1;
    Vec2i  atoms = (Vec2i){-1,-1};
    double l0=1,k=0,kpp=0;
    Vec3i8 ipbc=Vec3i8{0,0,0}; // for periodic boundary conditions

    //int    type;
    //Vec2i  atoms;
    //double l0,k;

    inline int getNeighborAtom(int ia)const{
        if     (ia==atoms.i){ return atoms.j; }
        else if(ia==atoms.j){ return atoms.i; }
        return -1;
    }

    void print()const{ printf( " Bond{t %i a(%i,%i) l0 %g k %g}", type, atoms.i, atoms.j, l0, k ); };

    Bond()=default;
    Bond(int type_, Vec2i atoms_, double l0_, double k_):type(type_),atoms(atoms_),l0(l0_),k(k_),ipbc{0,0,0}{};
};

// ============================
// ========= Angle
// ============================
struct Angle{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec2i  bonds = (Vec2i){-1,-1};
    //double a0    = 0;
    //double k     = 0;
    //Angle()=default;

    int type;
    Vec2i  bonds;
    double a0;
    double k;

    void print()const{ printf( " Angle{t %i b(%i,%i) a0 %g k %g}", type, bonds.i, bonds.j, a0, k ); }

    Angle()=default;
    Angle( int type_, Vec2i bonds_, double a0_, double k_):type(type_), bonds(bonds_),a0(a0_),k(k_){ };
};

// ============================
// ========= Dihedral
// ============================
struct Dihedral{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec3i  bonds = (Vec3i){-1,-1,-1};
    //int    n=0;
    //double k=0;

    int    type;
    Vec3i  bonds;
    int    n;
    double k;

    //Dihedral()=default;

    void print()const{ printf( " Dihedral{t %i b(%i,%i,%i) n %i k %g}", type, bonds.a, bonds.b,bonds.c, n, k ); }

    Dihedral()=default;
    Dihedral( int type_, Vec3i  bonds_, int n_, double k_ ):type(type_), bonds(bonds_), n(n_), k(k_){};
};


// ============================
// ========= Fragment
// ============================
//#ifdef Molecule_h
struct Fragment{
    int   imolType;
    Vec2i atomRange;
    Vec2i confRange;
    Vec2i bondRange;
    Vec2i angRange;
    Vec2i dihRange;

    Vec3d  pos;
    Quat4d rot;
    //Molecule * mol;     // ToDo : we can replace this by MolID to leave dependence on Molecule_h
    Vec3d    * pos0s;

    void finish(int ia,int ic, int ib,int ig, int id){ atomRange.y=ia; confRange.y=ic; bondRange.y=ib; angRange.y=ig; dihRange.y=id; };

    Fragment()=default;
    Fragment( Vec2i atomRange_, Vec2i confRange_, Vec2i bondRange_,Vec2i angRange_,Vec2i dihRange_):atomRange(atomRange_),confRange(confRange_),bondRange(bondRange_),angRange(angRange_),dihRange(dihRange_){};
    Fragment( int ia,int ic,int ib,int ig, int id):atomRange{ia,ia},confRange{ic,ic},bondRange{ib,ib},angRange{ig,ig},dihRange{id,id}{};
    //Fragment(Molecule* mol_, Vec3d pos_, Quat4d rot_, Vec2i atomRange_ ):
    Fragment(int imolType_, Vec3d pos_, Quat4d rot_, Vec2i atomRange_, Vec3d* pos0s_ ):
        imolType(imolType_),
        pos(pos_),rot(rot_),
        atomRange{atomRange_},bondRange{0,0},angRange{0,0},dihRange{0,0}, //mol{mol_},
        pos0s{pos0s_}{};
};
//#endif // Molecule_h

// ============================
// ========= splitByBond
// ============================

//int splitGraphs( int nb, Vec2i* bonds, int a0, int b0 ){
int splitGraphs( int nb, Vec2i* bonds, int b0, std::unordered_set<int>& innodes ){
    //printf( "# ======== splitGraphs() \n" );
    std::unordered_set<int> exbonds; // excluded bonds
    exbonds.insert(b0);
    int n0;
    do{ // Breadth-first search in bond graph with exclusion
        n0=innodes.size();
        //printf( "#### splitGraphs.n0= %i \n", n0 );
        //printf("inodes: "); for(int i:innodes){ printf("%i ",i); } ;printf("\n");
        for( int ib=0; ib<nb; ib++ ){                         // go over all bonds
            //printf( "ib %i n0 %i \n", ib, n0 );
            if( exbonds.find(ib) != exbonds.end() ) continue; // if bond in excluded bonds, skip
            const Vec2i& b = bonds[ib];
            int ia=-1;
            if     ( innodes.find(b.a) != innodes.end() ){ ia=b.b; }   // if atom.a in nodes, add atom b
            else if( innodes.find(b.b) != innodes.end() ){ ia=b.a; }   // if atom.b in nodes, add atom a
            if(ia>=0){
                //printf( "splitGraphs.add(ib=%i,ia=%i)\n", ib, ia );
                innodes.insert(ia);
                exbonds.insert(ib);
            }
        }
    }while( innodes.size()>n0 ); // as long as new nodes are added
    return innodes.size();
}

int splitByBond( int ib, int nb, Vec2i* bond2atom, Vec3d* apos, int* selection, Vec3d& ax, Vec3d& p0 ){
    const Vec2i& b = bond2atom[ib];
    ax = (apos[b.b]-apos[b.a]).normalized();
    std::unordered_set<int> innodes1; innodes1.insert( b.a );  std::unordered_set<int>* sel1;
    std::unordered_set<int> innodes2; innodes2.insert( b.b );  std::unordered_set<int>* sel2;
    MM::splitGraphs( nb, bond2atom, ib, innodes1 );
    MM::splitGraphs( nb, bond2atom, ib, innodes2 );
    if(innodes1.size()<innodes2.size()){ sel1=&innodes1; sel2=&innodes2; p0=apos[b.a]; }else{ sel1=&innodes2; sel2=&innodes1; ax.mul(-1.); p0=apos[b.b]; } 
    int n=0;
    for( int i:*sel1){ selection[n]=i; n++; };
    for( int i:*sel2){ selection[n]=i; n++; };
    return sel1->size();
}

// ============================
// ========= Builder
// ============================
class Builder{  public:
    int verbosity = 0;
    //bool bDEBUG = false;

    std::unordered_set<int> selection;

    //static int iDebug = 0;
    std::vector<Atom>       atoms;
    std::vector<Bond>       bonds;
    std::vector<Angle>      angles;
    std::vector<Dihedral>   dihedrals;
    std::vector<AtomConf>  confs;
    //std::vector<int>  atom_neighs;

    bool bPBC=false;
    Mat3d lvec = Mat3dIdentity;  // lattice vectors for PBC (periodic boundary conditions)
    //std::vector<Vec3i> bondPBC;

    MMFFparams* params = 0;  // Each function which needs this can take it as parameter
    std::vector       <std::string>*     atomTypeNames = 0;
    std::unordered_map<std::string,int>* atomTypeDict  = 0;

    std::vector<Fragment>   frags;

#ifdef Molecule_h
    //std::vector<Molecule> mols;
    std::vector<Molecule*>              molTypes;
    std::unordered_map<std::string,int> molTypeDict;
    std::unordered_map<size_t,size_t> fragTypes;
    std::unordered_map<size_t,size_t> mol2molType;
#endif // Molecule_h


    //MMFF forcefield parametes
    double Lepair        = 0.5; 
    double Kepair        = 10.0;
    double Ksp_default   = 1.0;
    double Kpp_default   = 0.5;
    double Kpp_e_default = 0.25;

    int ignoreType=-1;
    Vec3d defaultREQ  {  1.5, 0.0, 0.0 };
    Bond  defaultBond { -1, {-1,-1}, 1.5, 10.0 };
    //Angle defaultAngle{ -1, {-1,-1}, 0.0, 0.5 };
    Angle defaultAngle{ -1, {-1,-1}, M_PI, 0.5 };

    Bond bondBrush = defaultBond;

    int itypHcap  =-1;
    int itypEpair =-1;

    Atom capAtom      = Atom{ (int)NeighType::H,     -1,-1, {0,0,0}, Atom::HcapREQ };
    Atom capAtomEpair = Atom{ (int)NeighType::epair, -1,-1, {0,0,0}, {0,0,0} };
    Atom capAtomPi    = Atom{ (int)NeighType::pi,    -1,-1, {0,0,0}, {0,0,0} };
    Bond capBond      = Bond{ -1,  {-1,-1},  1.07, 100/const_eVA2_Nm };
    Vec3d    capUp   = Vec3d{0.0,0.0,1.0};
    bool bDummyPi    = false;
    bool bDummyEpair = false;
    bool bAddCaps    = true;
    bool bAddECaps   = false;

    // =================== Functions =====================


    // =================== Geometry Functions =====================

    Mat3d rotationFromAtoms( int i0, int i1, int i2 ){
        Mat3d rot;
        Vec3d center = atoms[i0].pos;
        rot.fromDirUp( (atoms[i1].pos-center).normalized(), atoms[i2].pos-center);
        return rot;
    }

    Vec3d vecBetweenAtoms(int i, int j){ return atoms[j].pos-atoms[i].pos; }

    int rayBonds( const Vec3d& ro, const Vec3d& rd, double R ){
        int    imin=-1;
        double rmin=1e+300;
        for(int i=0; i<bonds.size();i++){
            const Vec2i& b  = bonds[i].atoms;
            const Vec3d& p0 = atoms[b.a].pos;
            const Vec3d& p1 = atoms[b.b].pos;
            //double r = rayLineDist( ro, rd, p0, p1-p0 );
            double r = capsulaIntersect( ro, rd, p0, p1, R );
            if(r<rmin){ rmin=r; imin=i; }
        }
        return imin;
    }
    void bbox( Vec3d& pmin, Vec3d& pmax, int i0=0, int n=-1, bool bInit=true){
        natom_def(n,i0);
        if(bInit){ pmin=Vec3dmax; pmax=Vec3dmin; }
        for(int i=0; i<n; i++){ 
            const Vec3d& p=atoms[i].pos;
            pmin.setIfLower  (p);
            pmax.setIfGreater(p);
        }
    }
    void move_atoms     ( Vec3d dshift, int i0=0, int n=-1, bool bInit=true){ natom_def(n,i0); for(int i=0; i<n; i++){ atoms[i].pos.add(dshift); } }
    void transform_atoms( Mat3d M, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int n=-1             ){ natom_def(n,i0);  for(int i=0; i<n; i++){ Vec3d p; M.dot_to( atoms[i].pos-orig_old, p); p.add(orig_new); atoms[i].pos=p; } }
    void rotate_atoms   ( double angle, Vec3d axis=Vec3dZ, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int n=-1 ){ Mat3d M; M.fromRotation(angle,axis); transform_atoms( M,orig_old,orig_new,i0,n); }

    void changeCell( const Mat3d& lvs, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int n=-1 ){
        Mat3d M,MM; 
        //lvec.invert_to(M); 
        lvec.invert_T_to(M); 
        //MM.set_mmul_TN(lvec,M);
        MM.set_mmul_TN(lvs,M);
        //MM = M;
        //MM.set_mmul(lvec,M);
        //MM.set_mmul(lvs,M);
        //printf("DEBUG changeCell()  lvec\n"); lvec .print();
        //printf("DEBUG changeCell()  lvs\n"); lvs   .print();
        //printf("DEBUG changeCell()  M (inv(lvs))\n"); M .print();
        //printf("DEBUG changeCell() MM\n"); MM.print(); //exit(0);
        transform_atoms( MM,orig_old,orig_new,i0,n);
        lvec=lvs;
    }


    // =================== Functions =====================


    /*
    void selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
        int ib = rayBonds( ro, rd, 0.3 );
        std::unordered_set<int> innodes1; innodes1.insert( bonds[ib].a );
        std::unordered_set<int> innodes2; innodes2.insert( bonds[ib].b );
        splitGraphs( bonds.size(), &bonds[0], ib, innodes1 );
        splitGraphs( bonds.size(), &bonds[0], ib, innodes2 );
        std::unordered_set<int>* sel;
        if( innodes1.size()<innodes2.size()  ){ sel=innodes1; }else{ sel=innodes2; }
        selection.erase();
        for( int i:*sel){ selection.insert(i); };
    }
    */

    int findNthAtomOfType( int ityp, int n){
        int j=0;
        for(int i=0; i<atoms.size(); i++){
            if( atoms[i].type==ityp ){
                j++; if(j==n) return i;
            }
        }
        return -1;
    }

    void clearBonds(){
        //atoms.clear();
        bonds.clear();
        //bondPBC.clear();
        angles.clear();
        dihedrals.clear();
        confs.clear();
        for(Atom& a:atoms){ a.iconf=-1; };
    }

    void clear(){
        atoms.clear(); //printf("DEBUG a.1 \n");
        confs.clear();
        bonds.clear(); //printf("DEBUG a.2 \n");
        //bondPBC.clear();
        angles.clear();
        dihedrals.clear();
#ifdef Molecule_h
        //mols .clear(); //printf("DEBUG a.3 \n");
        frags.clear(); //printf("DEBUG a.4 \n");
        fragTypes.clear();
#endif // Molecule_h
        printf("MM::Builder::clear() DONE\n");
    }

    void bindParams( MMFFparams* params_ ){
        params = params_;
        atomTypeNames = &params->atomTypeNames;
        atomTypeDict  = &params->atomTypeDict;
        itypHcap  = params->atomTypeDict.at("H");
        itypEpair = params->atomTypeDict.at("E");
    }

    void initDefaultAtomTypeDict(){
        makeDefaultAtomTypeDict( atomTypeNames, atomTypeDict );
    }

    // ============= Add Capping Hydrogens

    const AtomConf* getAtomConf(int ia)const{
        int ic=atoms[ia].iconf;
        if(ic>=0){ return &confs[ic]; }
        return 0;
    }

    int getBondToNeighbor( int ia, int ja )const {
        const AtomConf* conf = getAtomConf(ia);
        if(conf){
            for(int i=0; i<conf->nbond; i++){
                int ib  = conf->neighs[i];
                if(ib<0) continue;
                int jai = bonds[ib].getNeighborAtom(ia);
                if(jai==ja){ return ib; }
            }
        }
        return -1;
    }

    inline int getBondByAtoms(int i, int j)const{
        int ib;
        ib = getBondToNeighbor( i, j ); if( ib>=0 ) return ib;
        ib = getBondToNeighbor( j, i ); if( ib>=0 ) return ib;
        return -1;
    }

    AtomConf* addConfToAtom( int ia, const AtomConf* conf=0 ){
        //printf( "MM::Builder.addConfToAtom ia %i \n", ia );
        int ic = confs.size();
        atoms[ia].iconf = ic;
        if   (conf){ confs.push_back( *conf      );                       }
        else       { confs.push_back(  AtomConf()); confs.back().init0(); }
        confs.back().iatom=ia;
        return &confs.back();
    }
    int tryAddConfsToAtoms( int i0=0, int imax=-1, int ignoreType=1 ){
        if(imax<0){ imax=atoms.size(); }
        int n=0;
        for(int ia=0;ia<imax;ia++){
            if( ignoreType == atoms[ia].type) continue;
            if( atoms[ia].iconf < 0 ){
                addConfToAtom( ia, 0 );
                n++;
            }
        }
        return n;
    }

    AtomConf* insertAtom(const Atom& atom, const AtomConf* conf ){
        atoms.push_back(atom);
        return addConfToAtom( atoms.size()-1, conf );
    }
    int insertAtom(const Atom& atom ){ atoms.push_back(atom); return atoms.size()-1; }

    int insertAtom( int ityp, const Vec3d& pos, const Vec3d* REQ=0, int npi=-1, int ne=0 ){
        Vec3d REQloc;
        if(REQ==0)
            if(params){ 
                AtomType& at = params->atypes[ityp];
                REQloc.x=at.RvdW;
                REQloc.y=at.EvdW;
                REQloc.z=0;
                REQ=&REQloc;
            }else{ REQ=&defaultREQ; }
        int iconf=-1;
        if(npi>=0){
            //printf( "insertAtom npi>0 => make Conf \n" );
            iconf=confs.size();
            confs.push_back( AtomConf(atoms.size(), npi, ne ) );
        }
        atoms.push_back( Atom    ( ityp,-1,iconf, pos, *REQ )  );
        //printf( "insertAtom[%i]", atoms.size() ); atoms.back().print(); puts("");
        return atoms.size()-1;
    }

    int insertAtom( int ityp, const Vec3d* pos=0, const Vec3d* REQ=0, int npi=-1, int ne=0 ){
        if(pos==0)pos=&Vec3dZero;
        int ia = insertAtom( ityp, *pos, REQ, npi, ne );
        //if(bConf) addConfToAtom( ia, 0 );
        return ia;
    }
    int insertAtom( std::string name, const Vec3d* pos=0, const Vec3d* REQ=0, int npi=-1, int ne=0){ int ityp = params->atomTypeDict.at(name); return insertAtom(ityp,pos,REQ,npi,ne); };

    void randomizeAtomPos(double R){
        for(int i=0;i<atoms.size();i++){
            atoms[i].pos.addRandomCube(R);
        }
    }

    void removeAtom(int i, bool checkBrute=true){
        int iend = atoms.size()-1;
        _swap( atoms[i], atoms[iend] );
        atoms.resize(iend);
    };

    bool tryAddBondToAtomConf( int ib, int ia, bool bCheck ){
        int ic = atoms[ia].iconf;
        //printf( "MM::Builder.addBondToAtomConf ia %i ib %i ic %i \n", ia, ib, ic );
        if(ic>=0){
            if(bCheck){
                int ing = confs[ic].findNeigh(ib);
                //printf( "a[%i]ng[%i] ing %i Found? %i \n", ia,ib, ing, 0<=ing );
                if( 0<=ing ){
                    //printf( " ia %i ib %i ing %i RETURN \n", ia, ib, ing );
                    return true; // neighbor already present in conf
                }
            }
            //printf( "MM::Builder.addBondToAtomConf ia %i ib %i ADDED \n", ia, ib );
            bool success = confs[ic].addBond(ib);
            if(!success){printf("ERROR: in confs[%i].addBond(%i) => exit \n", ic, ib); exit(0); }
            //int order = bonds[ib].type;
            //if(order>1){ for(int i=0; i<order-1; i++)confs[ic].addPi(); };
        }
        return false; // some neighbor already present ?
    }

    bool addBondToConfs( int ib, bool bCheck ){
        bool ret=false;
        const Bond& bond = bonds[ib];
        ret|=tryAddBondToAtomConf( ib, bond.atoms.i, bCheck );
        ret|=tryAddBondToAtomConf( ib, bond.atoms.j, bCheck );
        return ret; // some neighbor already present ?
    }
    bool tryAddBondsToConfs( int i0=0, int imax=-1, bool bCheck=true ){
        bool ret=false;
        if(imax<0) imax=bonds.size();
        for(int ib=i0; ib<imax; ib++){
            ret|=addBondToConfs( ib, true );
        }
        return ret;  // some neighbor already present ?
    }

    int insertBond(const Bond& bond ){
        int ib = bonds.size();
        bonds.push_back(bond);
        //printf( "insertBond[%i]", bonds.size() ); bonds.back().print(); puts("");
        tryAddBondToAtomConf( ib, bond.atoms.i, false );
        tryAddBondToAtomConf( ib, bond.atoms.j, false );
        return ib;
    }

    int insertBond( Vec2i ias, int order ){
        double k=0,l0=1.0;
        if(params){
            int iat=atoms[ias.i].type;
            int jat=atoms[ias.j].type;
            params->getBondParams( iat, jat, order, l0, k );
        }
        return insertBond( Bond(order, ias, l0, k) );
    };

    void assignBondParams( int ib ){
        Bond& b = bonds[ib];
        const Atom& ai = atoms[b.atoms.i];
        const Atom& aj = atoms[b.atoms.j];
        int order=1;
        if( (ai.iconf>=0)&&(aj.iconf>=0) ){ 
            const AtomConf& ci = confs[ai.iconf];
            const AtomConf& cj = confs[aj.iconf];
            order+=_min( ci.npi, cj.npi ); 
            //printf("assignBondParams[%i] (%i,%i|%i) pi(%i,%i) \n", ib,  ai.type, aj.type, order, confs[ai.iconf].npi, confs[aj.iconf].npi );

            // Assign pi-pi allignment 
            bool bpi=ci.npi>0;
            bool bpj=cj.npi>0; 
            if     ( bpi && bpj                       ){ b.kpp=Kpp_default;   }  // pi-pi alignement
            else if( bpi&&(cj.ne>0) || bpj&&(ci.ne>0) ){ b.kpp=Kpp_e_default; }  // pi-epair alignment

        }
        //getBondTypeId( ai.type, aj.type, uint8_t order );
        params->getBondParams( ai.type, aj.type, order, b.l0, b.k );
        //printf("assignBondParams[%i] (%i,%i|%i) -> l0 %g k %g \n", ib,  ai.type, aj.type, order,   b.l0, b.k );
    }

    void assignAllBondParams(){ for(int i=0;i<bonds.size();i++){ assignBondParams(i); } };

    //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    void addCap(int ia,const Vec3d& hdir, Atom* atomj ){
        int ja=atoms.size();
        Atom atom_tmp;
        if(atomj==0){
            atom_tmp=capAtom;
            atomj=&atom_tmp;
        }
        atomj->pos = atoms[ia].pos + hdir;
        insertAtom(*atomj);
        capBond.atoms.set(ia,ja);
        insertBond( capBond );
    }

    //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    void addBondedAtom(int ia, int ityp, bool bConf ){
        int ja=atoms.size();
        int npi=-1; if(bConf){npi=0;};
        insertAtom( ityp, 0, 0, npi, 0 );
        insertBond( {ia,ja}, 1 );
    }

    void makeConfGeom(int nb, int npi, Vec3d* hs){
        Mat3d m;
        if(nb==3){ // defined by 3 sigma bonds
            m.b.set_cross( hs[1]-hs[0], hs[2]-hs[0] );
            //printf( "nb,npi(%i,%i) m.b(%g,%g,%g) \n", nb, npi, m.b.x, m.b.y, m.b.z );
            m.b.mul( -1/m.b.norm() );
            if(npi==0){ // sp3 no-pi
                if( 0 < m.b.dot( hs[0]+hs[1]+hs[2] ) ){ m.b.mul(-1.); }
                hs[3]=m.b;
            }else{
                hs[3]=m.b;
            }
        }else if(nb==2){ // defined by 2 sigma bonds
            m.fromCrossSafe( hs[0], hs[1] );
            if      (npi==0){ // -CH2- like sp3 no-pi
                const double cb = 0.81649658092; // sqrt(2/3)
                const double cc = 0.57735026919;  // sqrt(1/3)
                hs[nb  ] = m.c*cc+m.b*cb;
                hs[nb+1] = m.c*cc-m.b*cb;
            }else if(npi==1){ // =CH- like  sp 1-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
                //printf("like =CH- H(%g,%g,%g) pi(%g,%g,%g,) \n", hs[nb].x, hs[nb].y, hs[nb].z, hs[nb+1].x, hs[nb+1].y, hs[nb+1].z );
            }else{            // #C- like sp 2-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
            }
        }else if(nb==1){
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (npi==0){ // -CH3 like sp3 no-pi
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
            }else if(npi==1){ // =CH2 like sp2 1-pi
                const double ca = 0.87758256189;  // 1/2
                const double cc =-0.5;            // sqrt(1/8)
                hs[nb  ] = m.c*cc + m.a*ca;
                hs[nb+1] = m.c*cc - m.a*ca;
                hs[nb+2] = m.b;
            }else{            // #CH sp  2-pi
                hs[nb  ] = m.c*-1;
                hs[nb+1] = m.b;
                hs[nb+2] = m.a;
            }
        }else if(nb==0){
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (npi==0){ //  CH4 like sp3 no-pi
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
                hs[nb+3] = m.c;
            }
        }
        //if( (nb==2) && (npi=1) ){ printf( "hs angles %g %g %g \n", hs[0].getAngle(hs[1])/M_PI, hs[0].getAngle(hs[2])/M_PI, hs[0].getAngle(hs[3])/M_PI ); }
    }


    void assignSp3Params( int ityp, int nb, int npi, int ne, int npi_neigh, Quat4d& par ){
        int iZ = params->atypes[ityp].iZ;
        if(npi==0){
            par.x=par.w=-0.3333333; // 108.5 deg
        }else if(npi==1){
            par.x=par.w=-0.5;       // 120 deg
        }else if(npi==2){
            par.x=par.w=-1.0;       // 180 deg
        }

        if( (npi_neigh>0)&&(npi==0) ){
            par.w=-0.1; 
            //par.w=0.0; 
        }
        //if(iZ=8){ // Oxygen
        //    if(npi==0){ par.x=0.0; } // sp2
        //}
    }

    int assignSp3Type( int ityp_old, int nb, int npi, int ne, int npi_neigh ){
        //  SMILE conf  example   nb   npi  ne   ntot       Ass  Asp   Kpp
        // ----- Carbon
        // -CH2- sp3       CH4        4    0    0    4      108   -     0
        // =CH-  sp2       CH2=CH2    3    1    0    4      120   -     1
        // #CH   sp1       HC#CH      2    2    0    4      180   -     0
        // ----- Nitrogen
        // -NH-  sp3       NH3        3    0    1    4      108   110   1
        // =N-   sp2       CH2=NH     2    1    1    4      120   90    1
        // #N    sp1       HCN        1    2    1    4      180   90    0
        // ----- Oxygen
        // -O-   sp3       H2O        2    0    2    4      108   90    1
        // =O    sp2       CH2=O      1    1    2    4      120   90    1
        // #O    sp1       CO         1    2    1    4      180   -     0
        // ----- Fuorine
        // -F    sp3/sp1   HF         1    0    1    2      180   -     0
        const AtomType& t = params->atypes[ityp_old];
        int iZ            = t.iZ;
        int ityp          = t.subTypes.array[npi];
        return ityp;
    }

    int assignSp3Type_pi( int ityp_old, int npi )const{
        const AtomType& t = params->atypes[ityp_old];
        int ityp          = t.subTypes.array[npi];
        return ityp;
    }

    void assignAllSp3Types(){
        //printf("DEBUG assignAllSp3Types() \n");
        for(int i=0; i<atoms.size(); i++){
            Atom& A = atoms[i];
            if( A.iconf>=0 ){
                int npi = confs[A.iconf].npi;
                A.type = assignSp3Type_pi( A.type, npi );
            }
        }
    }

    int getNeighType( int ia, int j, int* aneighs ){
        int ja = aneighs[ ia*4 + j];
        //printf( "getNeighType() ia %i j %i ja %i\n", ia, j, ja );
        if( (ja<0)||(ja>atoms.size()) ) return -1;
        return atoms[ ja ].type;
    }

    const char* getAtomTypeName(int ia){ return params->atypes[atoms[ia].type].name; };

    bool hasNeighborOfType( int ia, int n, const int* its, bool* bls, int* aneighs ){
        int ic = atoms[ia].iconf;
        for(int j=0; j<n; j++){ bls[j]=false; };
        if(ic<0) return false;
        const AtomConf& c = confs[ic];
        for(int i=0; i<4; i++){
            int ja = aneighs[i];
            int jt = atoms[ja].type;
            for(int j=0; j<n; j++){
                if( jt==its[j] ) bls[j]=true;
            }
        }
        //printf("hasNeighborOfType(%i|%s){", ia, getAtomTypeName(ia) );
        //for(int j=0; j<n; j++){ printf( "%i,", bls[j]); };
        //printf("}\n" );
        return true;
    }

    int assignSpecialTypes( int* aneighs ){
        //printf("#===========assignSpecialTypes() \n");
        // ------ C
        const int it_C_sp3 = params->getAtomType("C_sp3",true);
        const int it_C_sp2 = params->getAtomType("C_sp2",true);
        const int it_C_sp1 = params->getAtomType("C_sp1",true);
        const int it_C_CA  = params->getAtomType("C_CA" ,true);
        const int it_C_ene = params->getAtomType("C_ene",true);
        const int it_C_yne = params->getAtomType("C_yne",true);
        const int it_C_CH3 = params->getAtomType("C_CH3",true);
        const int it_C_ald = params->getAtomType("C_ald",true);
        const int it_C_COO = params->getAtomType("C_COO",true);
        // ------ O
        const int it_O_sp3  = params->getAtomType("O_sp3",true);
        const int it_O_sp2  = params->getAtomType("O_sp2",true);
        const int it_O_sp1  = params->getAtomType("O_sp1",true);
        const int it_O_OH   = params->getAtomType("O_OH" ,true);
        const int it_O_ald  = params->getAtomType("O_ald",true);
        const int it_O_sCOO = params->getAtomType("O_sCOO",true);
        const int it_O_pCOO = params->getAtomType("O_pCOO",true);
        //const int it_O_O  = params->getAtomType("O_");
        // ------ N
        const int it_N_sp3 = params->getAtomType("N_sp3",true);
        const int it_N_sp2 = params->getAtomType("N_sp2",true);
        const int it_N_sp1 = params->getAtomType("N_sp1",true);
        const int it_N_NH2 = params->getAtomType("N_NH2",true);
        // ------- H
        const int it_H_OH  = params->getAtomType("H_OH" ,true);
        const int it_H_COO = params->getAtomType("H_COO",true);
        const int it_H_NH2 = params->getAtomType("H_NH2",true);
        const int it_H_CH3 = params->getAtomType("H_CH3",true);
        const int it_H_ene = params->getAtomType("H_ene",true);
        const int it_H_yne = params->getAtomType("H_yne",true);
        const int it_H_ald = params->getAtomType("H_ald",true);
        //  C                    0        1        2        3
        const static int its_C  [4]{it_O_sp3,it_O_sp2,it_C_sp2,it_C_CA};
        //  N                    0        1
        const static int its_N[2]{it_C_sp2,it_C_CA};
        //  O                    0        1        2
        const static int its_O [3]{it_C_COO,it_C_sp2,it_C_CA};
        int na=atoms.size();
        bool bls[8];
        int nnew=0;
        //for(int i=0; i<na;i++){ printf("aneigs[%i]{%i,%i,%i,%i}\n", i, aneighs[i*4+0],aneighs[i*4+1],aneighs[i*4+2],aneighs[i*4+3]); }
        for(int ia=0; ia<na; ia++){
            int* ngs  = aneighs+ia*4;
            int itnew =-1;
            Atom& A   = atoms[ia];
            const AtomType& t = params->atypes[A.type];
            int iZ            = t.iZ;
            //printf( "atom[%i] %s=%i iZ %i \n", ia, t.name, A.type, iZ );
            switch (iZ){
                case 1: {  // H
                    int ingt =  getNeighType( 0, 0, ngs );       // if(ingt>-1)printf( "H[%i]-(%i|%s)\n", ia, ingt, params->atypes[ingt].name );
                    //printf("H ingt=%i\n", ingt );
                    if(ingt>=0){
                        if     ( ingt==it_O_OH     ){ itnew=it_H_OH;  }
                        else if( ingt==it_O_sCOO   ){ itnew=it_H_COO; }
                        else if( ingt==it_N_NH2    ){ itnew=it_H_NH2; }
                        else if( ingt==it_C_CH3    ){ itnew=it_H_CH3; }
                        else if( ingt==it_C_ene    ){ itnew=it_H_ene; }
                        else if( ingt==it_C_yne    ){ itnew=it_H_yne; }
                        else if( ingt==it_C_ald    ){ itnew=it_H_ald; }
                    }
                }break;
                case 6: { // C
                    //printf("C \n");
                    hasNeighborOfType( ia,4, its_C, bls, ngs  );
                    if( bls[0] && bls[1] && (A.type==it_C_sp2) ){ // COOH
                        itnew = it_C_COO;
                    }else if( bls[2] || bls[3] ){ // Conjugated
                        itnew = it_C_CA;
                    }
                }break;
                case 7: { // N
                    //printf("N \n");
                    //hasNeighborOfType( ia,2, its_N, bls, ngs );
                }break;
                case 8: { // O
                    //printf("O \n");
                    hasNeighborOfType( ia,3, its_O, bls, ngs );
                    if( bls[0] ){ // COOH
                        if      (A.type==it_O_sp3){ itnew=it_O_sCOO; }
                        else if (A.type==it_O_sp2){ itnew=it_O_pCOO; }
                    }
                }break;
            }
            if( (itnew>=0) && (itnew!=A.type) ){ 
                //printf( "atom[%i] type %i -> %i (%s) \n", ia, A.type, itnew, params->atypes[itnew].name );
                A.type = itnew;
                nnew++; 
            }
        }
        return nnew;
    }

    int assignSpecialTypesLoop( int nmax, int* aneighs ){
        int nnew=0;
        int itr=0;
        for(itr=0; itr<nmax; itr++){
            printf( "# --- assignSpecialTypesLoop[%i] \n", itr );
            int ni = assignSpecialTypes( aneighs ); 
            nnew+=ni; 
            if( ni==0 )break; 
        }
        printf("ERROR: assignSpecialTypesLoop not converged in %i iterations\n", itr ); exit(0);
        return nnew;
    }

    void addCaps( int ia, int ncap, int ne, int nb, const Vec3d* hs ){
        bool Hmask[]{1,1,1,1};
        //if(nH!=ncap) Hmask[rand()%ncap]=0;
        //bool breverse = (nH==2)&&(ncap==3);
        bool breverse;
        if(ncap<4){
            if(ne>0) Hmask[rand()%ncap]=0;
            breverse = (ne>1);
        }else{
            for(int i=0;i<ne;i++)Hmask[3-i]=0;
            breverse = 0;
        }
        //printf( "makeSPConf: atom[%i] ncap %i nH %i nb %i npi %i ne %i Hmask{%i,%i,%i,%i}  \n", ia, ncap, nH, nb,npi,ne,  (int)Hmask[0],(int)Hmask[1],(int)Hmask[2],(int)Hmask[3] );
        for(int i=0; i<ncap; i++){
            if     (Hmask[i]!=breverse){ addCap(ia,hs[nb+i],&capAtom     ); }
            else if(bDummyEpair       ){ addCap(ia,hs[nb+i],&capAtomEpair); }
        }
    }

    int addCapTopo(int ia){
        int ic = atoms[ia].iconf;
        if(ic<0) return -1;
        AtomConf& conf = confs[ic];
        int ncap = N_NEIGH_MAX-conf.nbond - conf.npi;
        AtomType&  typ = params->atypes[ atoms[ia].type ];
        int ne   = typ.nepair();
        int nH   = ncap-ne;
        //printf( "addCapTopo[%i] ne,nH(%i,%i) nb,npi(%i,%i) \n", ia, ne, nH, conf.nbond, conf.npi );
        for(int i=0; i<ne; i++){ addBondedAtom(ia,itypEpair,false); };
        for(int i=0; i<nH; i++){ addBondedAtom(ia,itypHcap ,false); };
        return nH;
    }
    int addAllCapTopo(){
        int n=0,na=atoms.size();
        for(int i=0;i<na;i++){
            int nH = addCapTopo(i);
            n+=nH;
        }
        return n;
    }


    void makeSPConf(int ia,int npi,int ne){
        //if(nH==0){ // ToDo : check reasonable limits npi, nh
        int ic = atoms[ia].iconf;
        AtomConf& conf = confs[ic];
        conf.clearNonBond();
        int nb   = conf.nbond;
        int ncap = 4-nb-npi;   // number of possible caps
        int nH   = ncap-ne;
        if(!bAddECaps) ncap=nH;
        //printf("-- "); println(conf);
        if(verbosity>=2)printf( "makeSPConf[ia=%i] nb,npi(%i,%i) ncap,nH,ne(%i,%i,%i)\n", ia, nb,npi, ncap, nH,ne );
        //Mat3d m;
        Vec3d hs[4];
        for(int i=0;i<nb;i++){
            int ib = conf.neighs[i];
            int ja = bonds[ib].getNeighborAtom(ia);
            hs[i]  = atoms[ja].pos - atoms[ia].pos;
            hs[i].normalize();
        }
        makeConfGeom(conf.nbond, npi, hs);
        if(bAddCaps && (ncap>0) ) addCaps( ia, ncap, ne, nb, hs );
        if(bDummyPi){ for(int i=0; i<npi; i++){ addCap(ia,hs[i+ncap+nb],&capAtomPi); } }
        conf.npi=npi;
        conf.ne =ne;
        //printf("-> "); println(conf);
    }

    bool autoConfEPi(int ia){
        int ic=atoms[ia].iconf;
        if(ic<0)return false;
        int ityp=atoms[ia].type;
        //params->assignRE( ityp, REQ );
        AtomConf& conf = confs[ic];
        conf.ne  = params->atypes[ityp].nepair();
        conf.fillIn_npi_sp3();
        //int npi = 4 - conf.nbond - ne - nH;
        //printf( "autoConfEPi[%i] typ %i ne %i npi %i \n", ia, ityp, ne, npi );
        //for(int i=0; i<ne;  i++)conf.addEpair();
        //for(int i=0; i<npi; i++)conf.addPi();
        return true;
    }
    int autoAllConfEPi( ){
        int n=0,na=atoms.size();
        for(int ia=0;ia<na;ia++){
            if( autoConfEPi(ia) ){n++;}
        }
        return n;
    }

    bool tryMakeSPConf(int ia, bool bAutoEPi=false){
        const AtomConf* conf = getAtomConf(ia);
        //printf("tryMakeSPConf %i conf %li\n", ia, (long)conf  );
        if(conf){
            //printf("tryMakeSPConf: proceed !!! \n"  );
            if(bAutoEPi)autoConfEPi(ia);
            makeSPConf(ia,conf->npi,conf->ne);
            return true;
        }
        return false;
    }

    int makeAllConfsSP( bool bAutoEPi=false ){
        int n=0,na=atoms.size();
        for(int i=0;i<na;i++){
            if(tryMakeSPConf(i,bAutoEPi)){n++;}
        }
        return n;
    }

    int cleanPis(){ int n=0; for( AtomConf& c:confs ){ n+=c.npi; c.clean_pi_sp3(); }; return n; }

    int countPiE(int& npi, int i0=0, int n=-1)const{
        nconf_def(n,i0);
        //printf("DEBUG countPiE() n %i i0 %i \n", n,i0);
        int ne=0; 
        npi=0;
        for(int i=0; i<n;i++){
            const AtomConf& c = confs[i0+i]; 
            npi+=c.npi;
            ne +=c.ne;
        }
        return ne;
    }

    int getAtom_npi(int ia)const{
        const Atom& A = atoms[ia];
        if(A.iconf<0)return 0;
        return confs[A.iconf].npi;
    }

    int getAtom_ne(int ia)const{
        const Atom& A = atoms[ia];
        if(A.iconf<0)return 0;
        return confs[A.iconf].ne;
    }

    // int getAtom_nb(int ia){
    //     const Atom& B = atoms[ja];
    //     if(A.iconf<0)return 0;
    //     return confs[A.iconf].npi;
    // }


    int countAtomPiNeighs(int ia)const{
        int npi=0;
        const Atom& A = atoms[ia];
        if(A.iconf<0) return -1;
        const AtomConf& c= confs[A.iconf];
        for(int i=0; i<c.nbond; i++){
            int ib = c.neighs[i];
            int ja = bonds[ib].getNeighborAtom(ia);
            if( getAtom_npi(ja)>0 ) npi++;
        }
        return npi;
    }

    // ============= Angles

    void addAnglesToBond( int ib, int n,const int* neighs, double a0, double k ){
        for(int j=0; j<n; j++){
            angles.push_back( (Angle){-1,  (Vec2i){ neighs[ib], neighs[j]}, a0,k} );
            //printf("[%li]",angles.size()); println(angles.back());
        }
    }

    void addAnglesUpToN( int n, const int* neighs, double a0, double k ){
        for(int i=0; i<n; i++){ addAnglesToBond( i, i, neighs, a0, k ); }
    }

    bool addAnglesToAtom( int ia, double ksigma, double kpi ){
        const AtomConf* conf = getAtomConf(ia);
        if(conf==0) return false;
        int nsigma = conf->nbond;
        //printf("addAnglesToAtom[%i] nsigma %i npi %i \n", ia, nsigma, conf->npi  );
        // ------ Pi Bonds
        if( bDummyPi && (conf->npi>0) ){
            //printf( "addAnglesToAtom[%i] angles to dummy Pi-orbital \n", ia );
            //nsigma -= conf->npi; // ToDo : Check this
            for(int i=0;i<conf->npi;i++){ addAnglesToBond( i+nsigma, i+nsigma, conf->neighs, M_PI_2, kpi ); }
        }
        // ------- Sigma bonds
        static const double a0s[]{ 0.0, 0.0, M_PI, 120*M_PI/180, 109.5*M_PI/180 };
        static const double Ks []{ 0.0, 0.0, 1.4,           1.2,            1.0 };
        if(nsigma>=2){
            int iangType = nsigma+conf->ne;
            double a0 = a0s[iangType];
            ksigma   *=  Ks[iangType];
            //printf( "addAnglesToAtom[%i] ns %i npi %i a0,ks %g %g   {%g,%g,%g,%g} %g \n", ia, nsigma, conf->npi, a0, ksigma, a0s[0],a0s[1],a0s[2],a0s[3] , a0s[nsigma] );
            addAnglesUpToN( nsigma, conf->neighs, a0, ksigma );
            if(verbosity>2){
                printf("addAnglesToAtom[%i] ", ia); printAtomConf(ia);
                printf( " Sigma(%i,a0[%i] %g, k %g)", ia,  nsigma,conf->npi, iangType,  a0*(180/M_PI), ksigma );
                if(bDummyPi && (conf->npi>0) ){ printf( " Pi(n%i,a0 %g, k %g)", conf->npi,  M_PI_2*(180/M_PI), kpi ); }else{ puts(""); };
            }
        }
        return true;
    }

    void autoAngles(double ksigma, double kpi){
        for(int i=0; i<atoms.size(); i++){
            if(atoms[i].iconf>=0){
                addAnglesToAtom( i, ksigma, kpi );
            }
        }
    }

    // =============== Dihedrals

    bool insertDihedralByAtom(const Quat4i& ias, Dihedral& tors ){
        int ib1 = getBondByAtoms(ias.x,ias.y); if(ib1<0) return false;
        int ib2 = getBondByAtoms(ias.y,ias.z); if(ib2<0) return false;
        int ib3 = getBondByAtoms(ias.z,ias.w); if(ib3<0) return false;
        tors.bonds.set(ib1,ib2,ib3);
        dihedrals.push_back(tors);
        return true;
    }

    /*
    void setAtoms( const Atom* brushAtom, const Vec3d* ps=0, int imin=0, int imax=-1,  ){
        if(imax<0){ imax=atoms.size(); }
        for(int i=imin;i<imax;i++){
            if(ps       )brushAtom.pos = ps[i];
            if(brushAtom)builder.insertAtom( *brushAtom, true);
        }
    }
    void setBondTypes( Bond brushBond, const Vec2i* bond2atom=0, int imin=0, int imax=-1 ){
        if(imax<0){ imax=bonds.size(); }
        for(int i=imin;i<imax;i++){
            if(bond2atom)brushBond.atoms=bond2atom[i];
            builder.insertBond(brushBond);
        }
    }
    */
    void insertAtoms( int n, Atom brushAtom, const Vec3d* ps, bool withConf=true ){
        for(int i=0;i<n;i++){
            brushAtom.pos = ps[i];
            if(withConf){ insertAtom( brushAtom, 0 ); }
            else        { insertAtom( brushAtom    ); }
        }
    }
    void insertBonds( int n, Bond brushBond, const Vec2i* bond2atom ){
        for(int i=0;i<n;i++){
            brushBond.atoms=bond2atom[i];
            insertBond( brushBond );
        }
    }
    void setConfs( int npi, int ne, int imin=0, int imax=-1 ){
        if(imax<0){ imax=atoms.size(); }
        for(int i=imin;i<imax;i++){
            makeSPConf( i, npi, ne );
        }
    }

    void touchingAtoms( int i0, int imax, const Vec3d& p, double R0, double Rfac, std::vector<int>& found ){
        for(int i=i0; i<imax; i++){  // for pbc we need all atom pairs
            const Atom& A = atoms[i];
            Vec3d dp = A.pos - p; // pbc here
            double R = (R0 + A.REQ.x)*Rfac;
            if(  dp.norm2() < (R*R) ){
                if(verbosity>2)  printf( "[%i,%i] r %g R %g \n", i0-1, i, dp.norm(), R );
                found.push_back(i);
            }
        }
    }

    void autoBonds( double R=-0.5, int i0=0, int imax=-1 ){
        // ToDo : periodic boundary conditions
        if(imax<0)imax=atoms.size();
        bool byParams = (R<0);
        double Rfac=-R;
        //if( byParams && (params==0) ){ printf("ERROR in MM::Builder.autoBonds() byParams(R<0) but params==NULL \n"); exit(0); }
        for(int i=i0; i<imax; i++){
            const Atom& A = atoms[i];
            for(int j=i+1; j<imax; j++){  // for pbc we need all atom pairs
                const Atom& B = atoms[j];
                Vec3d dp = B.pos - A.pos; // pbc here
                if(byParams){ R = (B.REQ.x + A.REQ.x)*Rfac; }
                if(  dp.norm2() < (R*R) ){
                    bondBrush.atoms={i,j};
                    insertBond( bondBrush );
                }
            }
        }
    }

    inline Vec3d pbcShift( Vec3i G ){ return lvec.a*G.a + lvec.b*G.b + lvec.c*G.c; }

    void autoBondsPBC( double R=-0.5, int i0=0, int imax=-1, Vec3i npbc=Vec3iOne ){
        if(verbosity>0)printf( "autoBondsPBC() \n" );
        if(verbosity>1){ printf( "builder.lvec: \n" ); lvec.print(); };
        bPBC = true;
        // ToDo : periodic boundary conditions
        if(imax<0)imax=atoms.size();
        bool byParams = (R<0);
        double Rfac=-R;
        //if( byParams && (params==0) ){ printf("ERROR in MM::Builder.autoBonds() byParams(R<0) but params==NULL \n"); exit(0); }
        std::vector<int> found;
        for(int i=i0; i<imax; i++){
            const Atom& A = atoms[i];
            R = A.REQ.x;
            int ipbc=0;
            if(verbosity>1)printf( "autoBondsPBC() Atom[%i] R %g \n", i, R );
            for(int ix=-npbc.x;ix<=npbc.x;ix++){
                for(int iy=-npbc.y;iy<=npbc.y;iy++){
                    for(int iz=-npbc.z;iz<=npbc.z;iz++){
                        int   j0=i+1;
                        //Vec3d vpbc = lvec.a*ix + lvec.b*iy + lvec.c*iz;
                        //Vec3d vpbc; lvec.dot_to_T( {(double)ix,(double)iy,(double)iz} );
                        //Vec3d p = A.pos - pbcShift( {ix,iy,iz} );
                        Vec3d p = A.pos - lvec.lincomb( ix, iy, iz );
                        found.clear();
                        touchingAtoms( j0, imax, p, R, Rfac, found );
                        //if(i==12)printf( "# pbc[%i,%i,%i][%i] nfound %i \n", ix,iy,iz, ipbc, found.size() );
                        for(int j:found){
                            //bondPBC.push_back( {ix,iy,iz} );
                            bondBrush.ipbc=Vec3i8{ix,iy,iz};
                            bondBrush.atoms={i,j};
                            insertBond( bondBrush );
                        }
                        ipbc++;
                    }
                }
            }
        }
    }

    // =============== Configurations

    int checkBond2Conf(bool bPrint)const{
        for(int i=0;i<bonds.size(); i++){
            //printf("checkBond2Conf b[%i]\n", i );
            const Bond& b = bonds[i];
            int i_ = getBondByAtoms(b.atoms.i,b.atoms.j);
            if(i_!= i){
                if(bPrint){
                    printf( "MMFFbuilder.checkBond2Conf: getBondByAtoms(bond[%i/%li]) returned %i \n", i,bonds.size(), i_ );
                }
                return i;
            }
        }
        return -1;
    }

    int checkConf2Bond(bool bPrint)const{
        int nb=0;
        std::vector<int> ng(atoms.size(), 0);
        for(const Bond& b: bonds){ ng[b.atoms.i]++; ng[b.atoms.j]++; };
        for(int ia=0;ia<atoms.size(); ia++){
            //printf("checkConf2Bond[%i] \n", ia );
            const AtomConf* conf = getAtomConf(ia); // we need to modify it
            if(conf==0){
                if( nb<bonds.size() ){
                    if(bPrint){
                        printf( "MMFFbuilder.checkConf2Bond: atom[%i/%li].conf==null nb(%i)<bonds.size(%li) \n", ia, atoms.size(), nb,bonds.size()  );
                    }
                    return ia;
                } else continue;
            }
            int nbconf = conf->nbond;
            if(nbconf != ng[ia] ){
                    if(bPrint){
                        printf( "MMFFbuilder.checkConf2Bond: atom[%i/%li].conf.nbond==%i but counted %i bonds \n", ia, atoms.size(), nbconf, ng[ia] );
                        println( (*conf) );
                    }
                    return ia;
            }
            for(int j=0; j<nbconf; j++){
                int ib = conf->neighs[j];
                int ja = bonds[ib].getNeighborAtom(ia);
                if(ja<0){
                    if(bPrint){
                        printf( "MMFFbuilder.checkConf2Bond: atom[%i/%li].neighs[%i/%i]->bonds[%i/%li].getNeighborAtom(%i) returned %i \n", ia,atoms.size(), j,nbconf, ib,bonds.size(), ia, ja );
                        println( (*conf)   );
                        println( bonds[ib] );
                    }
                    return ia;
                }
            }
            //printf("checkConf2Bond[%i] nb %i \n", ia, nb );
            nb+=nbconf;
        }
        return -1;
    }

    bool checkBondOrdered( int ib )const{ const Vec2i& b=bonds[ib].atoms; return b.a<b.b; }
    bool checkBondsOrdered( bool bOrder, bool bPrint ){ 
        bool ret=false;
        for(int ib=0; ib<bonds.size(); ib++){
            if( !checkBondOrdered( ib ) ){
                ret=true;
                if(bPrint) printf( "WARRNING bond[%i] is not ordered atoms(%i,%i) \n", ib, bonds[ib].atoms.a, bonds[ib].atoms.b );
                if(bOrder){ Vec2i& b=bonds[ib].atoms; _swap(b.a,b.b); bonds[ib].ipbc.mul(-1); }
            }
        }
        return ret;
    }

    bool checkAtomHasBond(int ia, int ib, bool bDefault=true)const{
        int ic = atoms[ia].iconf;
        if(ic>=0){
            int i = confs[ic].findNeigh(ib);
            //printf( "DEBUG ia=%i neighs[%i]=%i \n", ia, i, ib );
            return i>=0;
        }
        return bDefault;
    }

    bool checkNeighsRepeat( bool bPrint=true ){
        bool ret = false;
        for(int ia=0; ia<atoms.size(); ia++ ){
            const Atom& A = atoms[ia];
            if(A.iconf>=0){
                const AtomConf& c = confs[A.iconf];
                bool b = c.checkNeighsRepeat();
                if(b&bPrint){ printf("WARRNING repeating neighbors on atom " ); printAtomConf(ia); puts(""); }
                ret |=b;
            }
        }
        return ret;
    }

    bool checkBondInNeighs( int ib )const{ 
        const Vec2i& b=bonds[ib].atoms; 
        bool ba = checkAtomHasBond(b.a,ib);
        bool bb = checkAtomHasBond(b.b,ib);
        //if(! (ba && bb))printf( "DEBUG !bond[%i|%i,%i] %i %i\n", ib, b.a, b.b, ba, bb);
        return ba && bb;
    }

    bool checkBondsInNeighs( bool bPrint=true ){
        bool ret = false;
        for(int ib=0; ib<bonds.size(); ib++ ){
            //printf("\----bond[%i]\n", ib);
            if( ! checkBondInNeighs( ib ) ){
                if(bPrint){ Vec2i b=bonds[ib].atoms; printf("WARRNING bond[%i|%i,%i] is not in neighborhood \n", ib, b.a, b.b); printAtomConf(b.a); puts(""); printAtomConf(b.b); puts(""); }
                ret=true;
            }
        }
        return ret;
    }

    bool checkBondsSorted( int iPrint=0 )const{
        int ia=-1,ja=-1;
        if(iPrint>1)printf("checkBondsSorted %li \n", bonds.size() );
        for(int i=0;i<bonds.size(); i++){
            const Vec2i& b = bonds[i].atoms;
            if(iPrint>1)printf( "pair[%i] %i,%i | %i %i  | %i %i %i \n", i, b.i, b.j,   ia,ja ,   b.i>=b.j,  b.i<ia, b.j<=ja );
            if(b.i>=b.j){ if(iPrint>0){ printf("b.i>=b.j b[%i](%i,%i) ia,ja(%i,%i)\n", i,b.i,b.j,ia,ja); }; return false; }
            if(b.i<ia)  { if(iPrint>0){ printf("b.i<ia   b[%i](%i,%i) ia,ja(%i,%i)\n", i,b.i,b.j,ia,ja); }; return false; }
            else if (b.i>ia){ia=b.i; ja=-1; };
            if(b.j<=ja){  if(iPrint>0){ printf("b.j<=ja  b[%i](%i,%i) ia,ja(%i,%i)\n", i,b.i,b.j,ia,ja); }; return false; }
            ja=b.j;
        }
        if(iPrint>1)printf("... checkBondsSorted DONE !\n");
        return true;
    }

    void permutAtoms(int* permut, bool doBonds=false ){
        for(Bond&     b: bonds){ b.atoms.a=permut[b.atoms.a];  b.atoms.b=permut[b.atoms.b]; }
        // Confs are not re-shufled because they point to bonds, not atoms
        //for(AtomConf& c: confs){ 
        //    for(int j=0; j<N_NEIGH_MAX;j++){
        //        int& ing= c.neighs[j];
        //        if(ing>=0){ ing=permut[ing]; }
        //    }
        //}
        std::vector<Atom> atoms_(atoms);
        for(int ia=0; ia<atoms_.size(); ia++){
            int ja  = permut[ia];
            Atom& A = atoms [ja];
            A       = atoms_[ia];
            if( A.iconf>=0 ){
                confs[A.iconf].iatom = ja;
            }
        }
        if(doBonds){
            for(int i=0; i<bonds.size(); i++){
                Bond& b = bonds[i];
                b.atoms.set( permut[b.atoms.i], permut[b.atoms.j] );
            }
        }
    }

    void sortConfAtomsFirst(){
        int natom=atoms.size();
        int permut[natom];
        int j=0;
        for(int i=0; i<natom; i++){ if( atoms[i].iconf>=0 ){ permut[i] = j; j++; } }
        for(int i=0; i<natom; i++){ if( atoms[i].iconf< 0 ){ permut[i] = j; j++; } }
        //for(int i=0; i<natom; i++){ printf( "atom %i->%i \n", i, permut[i] ); }
        //printf( "natom %i \n", natom );
        permutAtoms(permut);
        //exit(0);
    }

    void sortAtomsOfBonds(){
        for(int i=0; i<bonds.size(); i++){ bonds[i].atoms.order(); }
    }


    bool sortBonds(){
        //printf( "sortBonds \n" );
        // sort bonds so that
        //   1) (b.i<b.j)
        //   1) if(bk.i) (b.i<b.j)
        //   1) (b.i<b.j)

        //int bsort    = new[bonds.size()];
        //Bond * bback = new Bond[bonds.size()];
        //int *   invBsort = new int     [bonds.size()];

        // use smart pointer to solve problems with delete[] when return on fail
        std::unique_ptr<Bond[]> bback   (new Bond[bonds.size()]);
        std::unique_ptr<int []> invBsort(new int [bonds.size()]);

        int nga[N_NEIGH_MAX];
        int ngb[N_NEIGH_MAX];

        sortAtomsOfBonds();
        //printBonds();

        int nb=0;
        for(int ia=0; ia<atoms.size(); ia++ ){
            // assume atoms with conf are first, capping are later
            //const AtomConf* conf = getAtomConf(ia);
            if(nb>=bonds.size())break;
            AtomConf* conf = (AtomConf*)getAtomConf(ia); // we need to modify it
            if(!conf){
                printf( "ERROR in MMFF.sortBonds(): atom[%i/%li] without conf (confs.size(%li)) met before all bonds enumerated nb(%i)<bonds.size(%li) \n", ia, atoms.size(), confs.size(), nb, bonds.size() );
                printf( " This algorithm assumes all atoms with conf precede atoms without confs in the array \n" );
                printf( " => return \n" );
                return false;
            }
            int nbconf=conf->nbond;
            int * neighs = conf->neighs;
            //printf( "ia %i nb %i conf.nb %i\n", ia, nb, nbconf );
            for(int i=0;i<nbconf;i++ ){
                int ib=neighs[i];
                if(ib<0){ printf("ERROR in MMFF.sortBonds(): atom[%i].condf inconsistent nbond=%i neigh[%i]<0 \n", ia, conf->nbond, i ); return false; }
                int ja = bonds[ib].getNeighborAtom(ia);
                //if(ja<ia)continue; // this bond was processed before
                nga[i]=ja;
                ngb[i]=ib;
            }
            int ja=-1;
            for(int i=0;i<nbconf;i++ ){      // take bonds on atom in order
                int ipick = selectMinHigher(ja, nbconf, nga );
                ja=nga[ipick];
                //neighs[i] = ngb[ipick]; // make conf sorted
                //printf( " atom[%i].neigh[%i] %i \n", ia, i, ja  );
                //printf( " atom[i %i -> j %i] ng %i \n", ia, ja, i  );
                if(ja<ia)continue;      // this bond was processed before (Hopefully)
                int ib = ngb[ipick];

                //bsort   [nb]=ib;
                bback[nb]   = bonds[ib];
                invBsort[ib]=nb;
                //printf( " bond[%i] -> bond[%i] \n", ib, nb );
                nb++;
            }
            // clean conf so it can be re-inserted
            conf->nbond=0;
            conf->n-=nbconf;
        }
        bonds.clear();
        for(int i=0; i<nb;i++){
            bback[i].atoms.order();
            //bonds[i]=bback[i];
            insertBond( bback[i] );
            //printf( " bond[%i] (%i,%i) \n", i, bback[i].atoms.i, bback[i].atoms.j );
        }
        for(int i=0; i<angles.size();i++){
            Vec2i& bs = angles[i].bonds;
            bs.a = invBsort[bs.a];
            bs.b = invBsort[bs.b];
        }
        for(int i=0; i<dihedrals.size();i++){
            Vec3i& bs = dihedrals[i].bonds;
            bs.a = invBsort[bs.a];
            bs.b = invBsort[bs.b];
            bs.c = invBsort[bs.c];
        }
        return true;
    }

    bool trySortBonds( int iPrint=0 ){
        bool sorted = checkBondsSorted( iPrint );
        if( !sorted ){
            if( !sortBonds() ){
                printf( " ERROR in builder.sortBonds() => exit \n" );
                exit(0);
            }
        }
        return sorted;
    }

    int findBondsToAtom(int ia, int* is=0, bool bPrint=false)const{
        int i=0;
        for(int ib=0; ib<bonds.size(); ib++){
            if( bonds[ib].atoms.anyEqual(ia) ){ 
                if(is){ is[i]=1; i++; }
                if(bPrint){ printf("bond[%i]",ib); bonds[ib].print(); puts(""); }
            };
        };
        return i;
    }

    void printSizes()const{ printf( "sizes: atoms(%i|%i) bonds(%i) angles(%i) dihedrals(%i) \n", atoms.size(), confs.size(), bonds.size(), angles.size(), dihedrals.size() ); };

    void printAtoms()const{
        printf(" # MM::Builder.printAtoms(na=%i) \n", atoms.size() );
        for(int i=0; i<atoms.size(); i++){
            printf("atom[%i]",i); atoms[i].print(); puts("");
        }
    }
    void printBonds()const{
        printf(" # MM::Builder.printBonds(nb=%i) \n", bonds.size() );
        for(int i=0; i<bonds.size(); i++){
            //printf("bond[%i]",i); bonds[i].print(); if(bPBC)printf(" pbc(%i,%i,%i)",bondPBC[i].x,bondPBC[i].y,bondPBC[i].z); puts("");
            printf("bond[%i]",i); bonds[i].print(); if(bPBC)printf(" pbc(%i,%i,%i)",bonds[i].ipbc.x,bonds[i].ipbc.y,bonds[i].ipbc.z); puts("");
            //printf("typs(%i,%i)", atoms[bonds[i].atoms.i].type, atoms[bonds[i].atoms.j].type ); 
        }
    }
    void printBondParams()const{
        printf(" # MM::Builder.printBondParams(nb=%i) \n", bonds.size() );
        for(int i=0; i<bonds.size(); i++){
            const Bond& b = bonds[i];
            printf("bond[%i]a(%i,%i)iZs(%i,%i)l0,k(%g,%g)\n",i, b.atoms.i,b.atoms.j, params->atypes[atoms[b.atoms.i].type].iZ, params->atypes[atoms[b.atoms.j].type].iZ, b.l0, b.k );
        }
    }
    void printAngles()const{
        printf(" # MM::Builder.printAngles(ng=%i) \n", angles.size() );
        for(int i=0; i<angles.size(); i++){
            printf("angle[%i]", i); angles[i].print(); puts("");
        }
    }
    void printConfs()const{
        printf(" # MM::Builder.printConfs(nc=%i) \n", confs.size());
        for(int i=0; i<confs.size(); i++){
            printf("conf[%i]", i); confs[i].print(); puts("");
        }
    }
    void printAtomConf(int i)const{
        const Atom& A = atoms[i];
        printf("atom[%i] %s=%i ic %i ", i, params->atypes[A.type].name, A.type, A.iconf);
        if(A.iconf>=0){
            const AtomConf& c = confs[A.iconf];
            //printf(" Conf[%i] n %i nb %i npi %i ne %i nH %i ", A.iconf, c.n, c.nbond, c.npi, c.ne, c.nH );
            c.print();
        }
    }

    void printAtomNeighs(int ia)const{
        const Atom& A = atoms[ia];
        printf("atom[%i,t%i,c%i] ", ia, A.type, A.iconf );
        if(A.iconf>=0){
            const AtomConf& c = confs[A.iconf];
            printf("nBPEH(%i|%i,%i,%i,%i) ngs(", c.n,c.nbond,c.npi,c.ne,c.nH);
            for(int i=0; i<N_NEIGH_MAX; i++){
                int ib = c.neighs[i];
                int ja=-2;
                if(ib>=0){
                    ja = bonds[ib].getNeighborAtom(ia);
                }
                printf("%i,", ja );
            }
            printf("]" );
        }
    }

    void printAtomConfs( bool bOmmitCap=true, bool bNeighs=false )const{
        printf(" # MM::Builder.printAtomConfs(na=%i,nc=%i) \n", atoms.size(), confs.size() );
        for(int i=0; i<atoms.size(); i++){ if( bOmmitCap && (atoms[i].iconf==-1) )continue;  if(bNeighs){printAtomNeighs(i);}else{printAtomConf(i);} puts(""); }
    }

    void printAtomGroupType( int ityp )const{
        //printf(" # MM::Builder.printAtomConfs(na=%i,nc=%i) \n", atoms.size(), confs.size() );
        for(int ia=0; ia<atoms.size(); ia++){ 
            const Atom& A=atoms[ia];
            if( (A.type!=ityp) || (A.iconf==-1) )continue; 
            const AtomConf& conf = confs[A.iconf];
            printf("atom(%i|%i):",ia,A.type);
            for(int i=0; i<conf.nbond; i++){
                int ib=conf.neighs[i];
                //Vec2i = bonds[ib].atoms;
                int ja = bonds[ib].getNeighborAtom(ia);
                printf("(%i|%i)",ja,atoms[ja].type);
            }
            puts("");
        }
    }

    void chargeByNeighbors( bool bClean, double factor=0.05, int niters=1, double damp=0.5 ){
        printf("chargeByNeighbors() bClean=%i factor=%g \n", bClean, factor );
        if(bClean)for(int ia=0; ia<atoms.size(); ia++){ atoms[ia].REQ.z = 0; }
        for(int itr=0; itr<niters; itr++){
            for(const Bond& b: bonds){
                Atom& A = atoms[b.atoms.i];
                Atom& B = atoms[b.atoms.j];
                double dQ = ( params->atypes[A.type].Eaff - params->atypes[B.type].Eaff  ) * factor;
                if(itr>0){  dQ += -A.REQ.z + B.REQ.z;  dQ*=damp; }
                A.REQ.z += dQ;
                B.REQ.z -= dQ;
            }
        }
        for(int ia=0; ia<atoms.size(); ia++){ printf("atom[%i] Q=%g \n", ia, atoms[ia].REQ.z ); }; //exit(0);
    }

    int write2xyz( FILE* pfile, const char* comment="#comment" )const{
        //write2xyz( pfile, atoms.size(), int* atypes, Vec3d* pos, const std::unordered_map<std::string,int>& atomTypeDict, const char* comment="#comment" ){
        fprintf(pfile, "%i\n",  atoms.size() );
        fprintf(pfile, "%s \n", comment      );
        for(int i=0; i<atoms.size(); i++){
            int ityp         = atoms[i].type;
            const Vec3d&  pi = atoms[i].pos;
            //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
            fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n", (*atomTypeNames)[ityp].c_str(), pi.x,pi.y,pi.z );
        };
        return atoms.size();
    }

    int save2xyz( const char * fname, const char* comment="#comment"  )const{
        FILE* pfile = fopen(fname, "w");
        if( pfile == NULL ) return -1;
        int n = write2xyz(pfile, comment );
        fclose(pfile);
        return n;
    }

    int saveMol( const char* fname ){
        FILE* pfile = fopen(fname, "w");
        if( pfile == NULL ) return -1;
        fprintf(pfile, "\n" );
        fprintf(pfile, "SimpleSimulationEngine::MMFFBuilder:: saveMol()\n" );
        fprintf(pfile, "\n" );
        fprintf(pfile, "%3i%3i  0  0  0  0  0  0  0  0999 V2000\n", atoms.size(), bonds.size() );
        for(int i=0; i<atoms.size(); i++){
            int ityp         = atoms[i].type;
            const Vec3d&  pi = atoms[i].pos;
            fprintf( pfile, "%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n",  pi.x,pi.y,pi.z, (*atomTypeNames)[ityp].c_str() );
        }
        for(int i=0; i<bonds.size(); i++){
            int ityp          = bonds[i].type;   if(ityp==-1) ityp=1;
            const Vec2i&  ats = bonds[i].atoms;
            fprintf( pfile, "%3i%3i%3i  0  0  0  0\n",  ats.a+1, ats.b+1, ityp );
        }
        fclose(pfile);
        return atoms.size() + bonds.size();
    }

    int load_xyz( const char * fname, bool noH=false, bool bConf=true ){
        if(verbosity>0)printf( "MM::Builder.load_xyz(%s)\n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        int natoms; Vec3d pos,REQ=defaultREQ; char at_name[8]; int npi,ne=0;
        const int nbuf=1024;
        char buff[nbuf]; char* line;
        line = fgets( buff, nbuf, pFile ); // number of atoms
        sscanf( line, "%i", &natoms );
        if(verbosity>1)printf( "natoms %i \n", natoms );
        line = fgets( buff, nbuf, pFile ); // comment, ignore
        int n0 = atoms.size();
        for(int i=0; i<natoms; i++){
            line     = fgets( buff, nbuf, pFile ); // comment, ignore
            int nret = sscanf( line,       "%s %lf %lf %lf %lf %i  ",    at_name, &pos.x, &pos.y, &pos.z, &REQ.z, &npi );
            if(verbosity>1)printf   (  ".xyz[%i] %s %lf %lf %lf %lf %i\n", i, at_name,  pos.x,  pos.y,  pos.z,  REQ.z,  npi  );
            if( nret < 5 ){ REQ.z=0;  };
            if( nret < 6 ){ npi  =-1; };
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                int ityp=it->second;
                if(verbosity>2)printf( " %s -> %i \n", at_name, ityp );
                if(params){
                    params->assignRE( ityp, REQ );
                    ne = params->atypes[ityp].nepair();
                }
                if( noH && (at_name[0]='H') && (at_name[1]='\0')  ) continue;
                insertAtom( it->second, pos, &REQ, npi, ne );
            }
        }
        return atoms.size() - n0;
    }

    inline void natom_def(int& n,int i0)const{ if(n<0){ n=atoms .size()-i0; }; }
    inline void nconf_def(int& n,int i0)const{ if(n<0){ n=confs .size()-i0; }; }
    inline void nbond_def(int& n,int i0)const{ if(n<0){ n=bonds .size()-i0; }; }
    inline void nang_def (int& n,int i0)const{ if(n<0){ n=angles.size()-i0; }; }
    inline void ndih_def (int& n,int i0)const{ if(n<0){ n=dihedrals.size()-i0; }; }

    void export_REQs(Vec3d* REQs, int i0=0, int n=-1)const{
        natom_def(n,i0);
        //_allocIfNull(REQs,n);
        for(int i=0; i<n; i++){ 
            //printf( "export_REQs[%i] REQ(%g,%g,%g)\n", i, atoms[i0+i].REQ.x,atoms[i0+i].REQ.y,atoms[i0+i].REQ.z  );
            REQs[i]= atoms[i0+i].REQ; }
    }

    void export_apos(Vec3d* apos, int i0=0, int n=-1)const{
        natom_def(n,i0);
        //_allocIfNull(apos,n);
        for(int i=0; i<n; i++){ apos[i]= atoms[i0+i].pos; }
    }

    int* export_atypes(int*& atypes, int i0=0, int n=-1)const{
        natom_def(n,i0);
        _allocIfNull(atypes,n);
        for(int i=0; i<n; i++){ atypes[i]= atoms[i0+i].type; }
        return atypes;
    }

    void export_bonds(Vec2i* b2a, double* l0s=0, double* ks=0, double* kPis=0, int i0=0, int n=-1)const{
        nbond_def(n,i0);
        for(int i=0; i<n; i++){
            const Bond& b  = bonds[i0+i];
            b2a[i] = b.atoms;
            if(ks )ks [i] = b.k;
            if(l0s)l0s[i] = b.l0;

            if(kPis){
                bool bi  = getAtom_npi( b.atoms.i )==1;
                bool bei = getAtom_ne ( b.atoms.i )>0;
                bool bj  = getAtom_npi( b.atoms.j )==1;
                bool bej = getAtom_ne ( b.atoms.j )>0;
                printf( "bond[%i|%i,%i]t(%i,%i) ipe(%i,%i) jpe(%i,%i) \n", i, b.atoms.i, b.atoms.j, atoms[b.atoms.i].type, atoms[b.atoms.j].type,    bi,bei,   bj, bej );
                if( bi&&bj ){                   // pi vs pi
                    kPis[i] =  0.25;
                }else if( (bi&&bej) || (bj&&bei) ){  // pi vs e-pair
                    kPis[i] = -0.25;
                }else{
                    kPis[i] = 0;                // other
                }
            }
            //printf( "export_bonds[%i] l0 %g k %g \n", i, b.l0, b.k );
        }
    }

    void export_angles(Vec2i* a2b, double* a0s=0, Vec2d* cs0s=0, double* ks=0, int i0=0, int n=-1)const{
        nang_def(n,i0);
        for(int i=0; i<n; i++){
            const Angle& a = angles[i0+i];
            a2b[i] = a.bonds;
            if(ks )ks [i] = a.k;
            if(a0s)a0s[i] = a.a0;
            if(cs0s)cs0s[i].fromAngle( a.a0 * 0.5 ); // NOTE: we divide angle by 2
        }
    }

    void export_dihedrals(Vec3i* d2b, int* tn=0, double* ks=0, int i0=0, int n=-1)const{
        ndih_def(n,i0);
        for(int i=0; i<n; i++){
            const Dihedral& d = dihedrals[i0+i];
            d2b[i] = d.bonds;
            if(ks)ks[i] = d.k;
            if(tn)tn[i] = d.n;
        }
    }
    
    void assignAtomREQs( const MMFFparams* params ){
        for(int i=0; i<atoms.size(); i++){
            //mmff->aLJq [i]  = atoms[i].type;
            int ityp = atoms[i].type;
            atoms[i].REQ.x = params->atypes[ityp].RvdW;
            atoms[i].REQ.y = params->atypes[ityp].EvdW;
            atoms[i].REQ.z = 0;
            //atomTypes[i]  = atoms[i].type;
        }
    }

    void startFragment (         ){                          frags.push_back( Fragment( atoms.size(), confs.size(), bonds.size(), angles.size(), dihedrals.size() ) ); }
    void finishFragment(int i=-1 ){ if(i<0)i=frags.size()-1; frags[i].finish(           atoms.size(), confs.size(), bonds.size(), angles.size(), dihedrals.size()   ); }

    void insertAtoms( int na, int* atypes, Vec3d* apos, Vec3d* REQs=0, double* qs=0, int* npis=0, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity ){
        //printf( "# MM::Builder::insertFlexibleMolecule  natoms %i nbonds %i \n", mol->natoms, mol->nbonds );
        //startFragment();
        //int natom0  = atoms.size();
        //int nbond0  = bonds.size();
        for(int i=0; i<na; i++){
            Vec3d p,REQ;
            int ne=0,npi=0;
            int ityp = atypes[i];
            if( ityp==ignoreType ) continue;
            double q=0; if(qs){ q=qs[i]; }
            if(REQs  ){ REQ=REQs[i];                        }
            else      { params->assignRE( ityp, REQ,true ); }
            if(params){ ne = params->atypes[ityp].nepair(); npi=params->atypes[ityp].npi(); }
            if(npis  ){npi=npis[i];}
            REQ.z=q;
            //printf( "insert Atom[%i] ityp %i REQ(%g,%g,%g) npi,ne %i %i \n", i, ityp, REQ.x, REQ.y, REQ.z, npi, ne  );
            rot.dot_to( apos[i],p); p.add( pos );
            insertAtom( ityp, p, &REQ, npi, ne );
        }
        //finishFragment();
    }

    void insertBonds( int nb, Vec2i* b2as, int* btypes ){
        int natom0  = atoms.size();
        for(int i=0; i<nb; i++){
            int btyp = defaultBond.type;
            if(btypes){ btyp=btypes[i]; }
            bonds.push_back( Bond( btyp, b2as[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k ) );
        }
    }

    int insertAtomsBonds( int nb, Vec2i* b2as, int* btypes, int na, int* atypes,  Vec3d* apos, Vec3d* REQs=0, double* qs=0, int* npis=0, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity ){
        //printf( "# MM::Builder::insertFlexibleMolecule  natoms %i nbonds %i \n", mol->natoms, mol->nbonds );
        startFragment();
        insertAtoms( na, atypes, apos, REQs, qs, npis, pos,rot);
        insertBonds( nb, b2as, btypes );
        finishFragment();
        return frags.size()-1;
    }

    void insertFrag( int ifrag, const Vec3d& pos0=Vec3dZero, const Vec3d& pos1=Vec3dZero, const Mat3d& rot=Mat3dIdentity ){
        Vec2i atomRange = frags[ifrag].atomRange;
        Vec2i bondRange = frags[ifrag].bondRange;
        int ia0=atoms.size();
        int ib0=bonds.size();
        int dib=ib0-bondRange.a;
        int dia=ia0-atomRange.a;
        for(int i=atomRange.a; i<atomRange.b; i++){
            Atom A = atoms[i];
            rot.dot_to( A.pos-pos0,A.pos); A.pos.add(pos1);
            if(A.iconf>=0){ 
                AtomConf C = confs[A.iconf];
                C.rebase(dib); 
                C.iatom=atoms.size();
                A.iconf=confs.size();
                confs.push_back(C);
            };
            atoms.push_back(A);
        }
        for(int i=bondRange.a; i<bondRange.b; i++){
            Bond B = bonds[i];
            B.atoms.add(dia,dia);
            bonds.push_back(B);
        }
    }

    void multFragPBC( int ifrag, const Vec3i& nPBC, Mat3d lvs ){
        printf( "multFragPBC() nPBC(%i,%i,%i) ifrag=%i ) \n", nPBC.x, nPBC.y, nPBC.z,  ifrag );
        lvs.print();
        printSizes();
        for(int ic=0; ic<nPBC.c; ic++){ for(int ib=0; ib<nPBC.b; ib++){ for(int ia=0; ia<nPBC.a; ia++){
            if( (ia==0)&&(ib==0)&&(ic==0) ) continue;
            Vec3d pos0 = lvs.a*ia + lvs.b*ib + lvs.c*ic;
            startFragment();
            insertFrag( ifrag, Vec3dZero, pos0 );
            finishFragment();
        }}}
        lvec.a.set_mul(lvs.a,nPBC.a);
        lvec.b.set_mul(lvs.b,nPBC.b);
        lvec.c.set_mul(lvs.c,nPBC.c);
        printSizes();
    }

    int findNearestBondPBC( int ib, int& ifrag0, int ifragMin, int ifragMax, Vec3i8& ipbc ){
        int ia0  = frags[ifrag0].atomRange.a;
        const Vec2i & b     = bonds[ib].atoms;
        const Vec3i8& ipbc0 = bonds[ib].ipbc;
        const Vec3d & pi = atoms[b.i].pos;
        //Vec3d pj = atoms[b.j].pos;
        int jfragMin=-1;
        int jaMin   =-1;
        double r20 = (atoms[b.j].pos-pi).norm2();
        for(int jfrag=ifragMin; jfrag<ifragMax; jfrag++ ){
            if(ifrag0==jfrag) continue;
            int ja0  = frags[jfrag].atomRange.a;
            int ja   =  b.j - ia0 + ja0;
            double r2 = (atoms[ja].pos-pi).norm2();
            if(r2<r20){
                jfragMin=jfrag; jaMin=ja;
                ipbc=Vec3i8{0,0,0};
            }else{
                r2 = (atoms[ja].pos-pi + lvec.a*ipbc0.a + lvec.b*ipbc0.b + lvec.c*ipbc0.c ).norm2();
                if(r2<r20){
                    jfragMin=jfrag; jaMin=ja;
                    ipbc=ipbc0;
                }
            }
        }
        ifrag0=jfragMin;
        return jaMin;
    }

    int removeBondFromConfs(int ib){
        int ic,n=0;
        const Vec2i& b  = bonds[ib].atoms;
        //printf( "BEFORE removeBondFromConfs[%i|%i,%i] \n", ib, b.i, b.j ); printAtomConf(b.i); puts(""); printAtomConf(b.j); puts("");
        ic=atoms[b.i].iconf; if(ic>=0) n+=confs[ic].replaceNeigh(ib,-1);
        ic=atoms[b.j].iconf; if(ic>=0) n+=confs[ic].replaceNeigh(ib,-1);
        //b.set(ia,ja);
        //printf( "AFTER removeBondFromConfs[%i|%i,%i] \n", ib, b.i, b.j ); printAtomConf(b.i); puts(""); printAtomConf(b.j); puts("");
        return n;
    }

    int addBondToConfs(int ib){
        int ic,n=0;
        const Vec2i& b  = bonds[ib].atoms;
        //printf( "BEFORE addBondToConfs[%i|%i,%i] \n", ib, b.i, b.j ); printAtomConf(b.i); puts(""); printAtomConf(b.j); puts("");
        ic=atoms[b.i].iconf; if(ic>=0) n+=confs[ic].replaceNeigh( -1,ib);
        ic=atoms[b.j].iconf; if(ic>=0) n+=confs[ic].replaceNeigh( -1,ib);
        //printf( "AFTER addBondToConfs[%i|%i,%i] \n", ib, b.i, b.j ); printAtomConf(b.i); puts(""); printAtomConf(b.j); puts("");
        return n;
    }

    int correctPBCbonds( int ifragMin, int ifragMax ){
        std::vector<int> ibs; ibs.reserve(100);
        //printf( "correctPBCbonds(ifragMin=%i, ifragMax=%i) \n", ifragMin, ifragMax );
        for(int ifrag=ifragMin; ifrag<ifragMax; ifrag++ ){
            //printf( "correctPBCbonds[ifrag=%i] \n", ifrag );
            const Vec2i& bondRange  = frags[ifrag].bondRange;
            for(int ib=bondRange.a; ib<bondRange.b; ib++){
                //printf( "correctPBCbonds[ib=%i] \n", ib );
                Bond& b = bonds[ib ];
                if( b.ipbc.allEqual(0) ) continue;
                //printf( "correctPBCbonds[ib=%i] (%i,%i) ipbc(%i,%i,%i)\n", ib,  b.atoms.i,b.atoms.j, b.ipbc.x,b.ipbc.y,b.ipbc.z );
                int ifrag0 = ifrag; Vec3i8 ipbc;
                int ja = findNearestBondPBC( ib, ifrag0, ifragMin, ifragMax, ipbc);
                if(ja>=0){ // found?
                    //printf( "correctPBCbonds[%i] (%i,%i)->(%i,%i) \n", ib,  b.atoms.i,b.atoms.j,   b.atoms.i,ja );
                    int nremove = removeBondFromConfs( ib );   if(nremove==0){ printf("!!!! bond[%i] not removed from any conf \n", ib); };
                    b.atoms.set(b.atoms.i,ja);
                    b.ipbc = ipbc;
                    ibs.push_back(ib);
                }  
            }
        }
        for(int ib:ibs){
            int nadd = addBondToConfs(ib);
            if(nadd==0){ printf("!!!! bond[%i] not added to any conf \n", ib); };
        }
        return ibs.size();
    }

#ifdef Molecule_h

    int substituteMolecule(Molecule* mol, Vec3d up, int ib, int ipivot=0, bool bSwapBond=false, const Vec3i* axSwap=0, Mat3d* rot_=0 ){
        Bond& b = bonds[ib];
        int ia  = b.atoms.i;
        int ja  = b.atoms.j;
        if(bSwapBond) _swap(ia,ja);
        if(verbosity>0)printf( "substituteMolecule() ib=%i ipivot=%i ja=%i \n", ib, ipivot, ja  );
        Atom& A    = atoms[ja];
        Vec3d dir  = A.pos-atoms[ia].pos; dir.normalize();
        Mat3d rot; rot.fromDirUp( dir, up );
        if( axSwap ){ rot.swap_vecs(*axSwap); }
        if(rot_){ *rot_=rot; }
        //printf( "substituteMolecule() dir=(%g,%g,%g) \n", dir.x,dir.y,dir.z  );
        //printf( "substituteMolecule() up =(%g,%g,%g) \n", up.x,up.y,up.z  );
        //printf( "substituteMolecule() rot= \n" ); rot.print();
        { // set pivot atom
            int atyp = mol->atomType[ipivot];
            A.type = atyp;
            A.REQ  = mol->REQs[ipivot];
            if(A.iconf>=0){ 
                confs[A.iconf].init0(); confs[A.iconf].ne = params->atypes[atyp].nepair(); 
                //printf( "DEBUG pivot has conf ic=%i \n", A.iconf );
            }else{ 
                //printf( "DEBUG pivot has NOT conf ic=%i \n", A.iconf );
                //addConfToAtom(ia); 
            } 
        }
        int natom0  = atoms.size();
        for(int i=1; i<mol->natoms; i++){
            int ne=0,npi=0;
            Vec3d REQ=mol->REQs[i];
            int ityp = mol->atomType[i];
            if(params){
                //printf( "params \n" );
                params->assignRE( ityp, REQ );
                ne = params->atypes[ityp].nepair();
                REQ.z=mol->REQs[i].z;
            }
            if( mol->npis ) npi=mol->npis[i];
            //printf( "insert Atom[%i] ityp %i REQ(%g,%g,%g) npi,ne %i %i \n", i, ityp, REQ.x, REQ.y, REQ.z, mol->npis[i], ne  );
            Vec3d p; rot.dot_to_T(mol->pos[i]-mol->pos[ipivot],p); p.add( A.pos );
            insertAtom( ityp, p, &REQ, npi, ne );
        }
        //printf( "natom0 = %i \n", natom0 );
        for(int i=0; i<mol->nbonds; i++){
            //bonds.push_back( (Bond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k } );
            Vec2i b = mol->bond2atom[i];
            if(b.i==ipivot){ b.i=ja; }else{ b.i+=natom0; if(b.i>ipivot)b.i--; };
            if(b.j==ipivot){ b.j=ja; }else{ b.j+=natom0; if(b.j>ipivot)b.j--; };
            //printf("add bond[%i] a(%i,%i) as a(%i,%i)\n", i,  mol->bond2atom[i].i,  mol->bond2atom[i].j,   b.i, b.j );
            bonds.push_back( Bond(mol->bondType[i], b, defaultBond.l0, defaultBond.k ) );
        }
        return ja;
    }


    void clearMolTypes( bool deep ){
        if(deep){ for(Molecule* mol : molTypes ){ mol->dealloc(); delete mol; } }
        molTypeDict.clear();
        molTypes.clear();
    }

    int loadMolTypeXYZ(const char* fname, MMFFparams* params_=0 ){
        if(params_!=0) params=params_;
        Molecule* mol = new Molecule();
        mol->atomTypeDict = &params->atomTypeDict;
        //printf("mol->atypNames %i %i \n", mol->atypNames, &params->atypNames );
        //int iret = mol->load_xyz( fname ); 
        int iret =  params->loadXYZ( fname, mol->natoms, &mol->pos, &mol->REQs, &mol->atomType, &mol->npis, &lvec );
        //mol->printAtomInfo(); //exit(0);
        if( iret>0 ){ bPBC=true; printf("lvec loaded from %s \n", fname); };
        if(iret<0)return iret;
        if(params) params->assignREs( mol->natoms, mol->atomType, mol->REQs );
        int ityp = molTypes.size();
        mol2molType[(size_t)mol]=ityp;
        molTypes.push_back(mol);
        return molTypes.size()-1;
    }

    int registerRigidMolType( int natoms, Vec3d* pos, Vec3d* REQs, int* atomType ){
        Molecule* mol = new Molecule();
        mol->allocate( natoms, 0 );
        for(int i=0; i<mol->natoms; i++){ mol->pos[i]=pos[i]; mol->REQs[i]=REQs[i]; mol->atomType[i]=atomType[i]; }
        int ityp = molTypes.size();
        mol2molType[(size_t)mol]=ityp;
        molTypes.push_back(mol);
        return molTypes.size()-1;
    }

    int loadMolType(const std::string& fname, const std::string& label, MMFFparams* params=0, bool bCOG0=false ){
        //printf( "fname:`%s` label:`%s` \n", fname.c_str(), label.c_str()  );
        int itype = loadMolTypeXYZ( fname.c_str(), params );
        if(itype<0) return itype;
        //printf("bCOG0 %i \n", bCOG0 ); exit(0);
        if(bCOG0) molTypes[itype]->cogTo0();
        //printf( "fname:`%s` label:`%s` itype %i \n", fname.c_str(), label.c_str(), itype  );
        molTypeDict[label] = itype;
        return itype;
    }

    void insertFlexibleMolecule_ignorH( Molecule * mol, const Vec3d& pos, const Mat3d& rot, int iH = 1 ){
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        std::vector<int> atomInds(mol->natoms);
        std::vector<int> bondInds(mol->nbonds);
        for(int i=0; i<mol->natoms; i++){
            if( mol->atomType[i]==iH ) continue;
            atomInds[i] = atoms.size();
            Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y);
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (Atom){mol->atomType[i], -1, -1, p, REQi } );
        }
        for(int i=0; i<mol->nbonds; i++){
            //bonds.push_back( (Bond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k } );
            bondInds[i] = bonds.size();
            const Vec2i& b = mol->bond2atom[i];
            bonds.push_back( Bond(mol->bondType[i], { atomInds[b.a],atomInds[b.b] }, defaultBond.l0, defaultBond.k ) );
        }
        for(int i=0; i<mol->nang; i++){
            const Vec2i& ang = mol->ang2bond[i];
            angles.push_back( (Angle){ 1, { bondInds[ang.a],bondInds[ang.b] }, defaultAngle.a0, defaultAngle.k } );
        }
    }

    void insertFlexibleMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot, int ignoreType=-1 ){
        //printf( "# MM::Builder::insertFlexibleMolecule  natoms %i nbonds %i \n", mol->natoms, mol->nbonds );
        startFragment();
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        for(int i=0; i<mol->natoms; i++){
            int ne=0,npi=0;
            Vec3d REQ=mol->REQs[i];
            int ityp = mol->atomType[i];
            if( ityp==ignoreType ) continue;
            if(params){
                //printf( "params \n" );
                params->assignRE( ityp, REQ );
                ne = params->atypes[ityp].nepair();
                REQ.z=mol->REQs[i].z;
            }
            if( mol->npis ) npi=mol->npis[i];
            //printf( "insert Atom[%i] ityp %i REQ(%g,%g,%g) npi,ne %i %i \n", i, ityp, REQ.x, REQ.y, REQ.z, mol->npis[i], ne  );
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            insertAtom( ityp, p, &REQ, npi, ne );
        }
        for(int i=0; i<mol->nbonds; i++){
            //bonds.push_back( (Bond){mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k } );
            bonds.push_back( Bond(mol->bondType[i], mol->bond2atom[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k ) );
        }
        for(int i=0; i<mol->nang; i++){
            double alfa0 = defaultAngle.a0;
            if( mol->ang0s ) alfa0 = mol->ang0s[i];
            angles.push_back( (Angle){ 1, mol->ang2bond[i] + ((Vec2i){nbond0,nbond0}), alfa0, defaultAngle.k } );
            //printf( "angle[%i|%i,%i] %g|%g %g \n", i, angles.back().bonds.a, angles.back().bonds.b, angles.back().a0, alfa0, angles.back().k );
        }
        finishFragment();
    }

    int insertRigidMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot ){
        printf( "insertRigidMolecule() \n" );
        int natoms0 = atoms.size();
        Quat4d qrot; qrot.fromMatrix(rot);
        int ifrag = frags.size();
        //printf( "insertMolecule mol->natoms %i \n", mol->natoms );
        for(int i=0; i<mol->natoms; i++){
            //Vec3d REQi = Vec3d{1.0,0.03,mol->}; // TO DO : LJq can be set by type
            //atoms.push_back( (Atom){mol->atomType[i],mol->pos[i], LJq } );
            Vec3d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y); // REQi.z = 0.0;
            Vec3d  p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (Atom){mol->atomType[i], ifrag, -1, p, REQi } );
        }
        //size_t mol_id = (size_t)(mol);
        //frags.push_back( (Fragment){natoms0, atoms.size()-natoms0, pos, qrot, mol}  );
        //frags.push_back( Fragment{ mol, pos, qrot,  {natoms0, atoms.size()-natoms0} } );
        molTypes.push_back(mol);
        size_t mol_id = molTypes.size()-1;
        frags.push_back( Fragment{ mol_id, pos, qrot, {natoms0, atoms.size()-natoms0}, mol->pos } );
        //size_t mol_id = static_cast<size_t>(mol);
        auto got = fragTypes.find(mol_id);
        if ( got == fragTypes.end() ) {
            fragTypes[ mol_id ] = frags.size()-1; // WTF ?
        }else{}
        return ifrag;
    }


    int insertMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot, bool rigid, int noH=-1 ){
        //int natom0  = atoms .size();
        //int nbond0  = bonds .size();
        //int nangle0 = angles.size();
        //mols.push_back( (MMFFmol){mol, (Vec3i){natom0,nbond0,nangle0} } );
        if( rigid ){
            return insertRigidMolecule( mol, pos, rot );
        }else{
            if(noH>0){ insertFlexibleMolecule_ignorH( mol, pos, rot, noH ); }
            else     { insertFlexibleMolecule       ( mol, pos, rot      ); }
            return -1;
        }
    }


    int insertMolecule( int itype                 , const Vec3d& pos, const Mat3d& rot, bool rigid ){ return insertMolecule( molTypes[itype]                 , pos, rot, rigid ); }
    int insertMolecule( const std::string& molName, const Vec3d& pos, const Mat3d& rot, bool rigid ){ return insertMolecule( molTypes[ molTypeDict[molName] ], pos, rot, rigid ); }

    int insertFlexibleMolecule( int itype                 , const Vec3d& pos, const Mat3d& rot, int ignoreType=-1 ){ if(itype<0) return itype; insertFlexibleMolecule( molTypes[itype], pos, rot, ignoreType ); return 0; }
    int insertFlexibleMolecule( const std::string& molName, const Vec3d& pos, const Mat3d& rot, int ignoreType=-1 ){ return insertFlexibleMolecule( molTypeDict[molName], pos, rot, ignoreType ); }
#endif // Molecule_h


#ifdef ForceField_h
    void toForceField( ForceField& ff ){
        if(iDebug>0) printf( " MMFFbuilder.toForceField na %li nb %li nA %li nd %li \n", atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        //mmff->deallocate();
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        for(int i=0; i<atoms.size(); i++){
            ff.apos [i]  = atoms[i].pos;
            if(iDebug>0){ printf("[%i]", i); atoms[i].print(); if( atoms[i].iconf>=0){confs[atoms[i].iconf].print();} puts(""); }
        }
        for(int i=0; i<bonds.size(); i++){
            const Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            ff.bond2atom[i]    = ib;
            //if(params){
            //    params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, ff.bond_l0[i], ff.bond_k[i] );
            //}else{
            //    //printf( "no params \n" );
            //    ff.setBondParam(i, b.l0, b.k );
            //}
            ff.setBondParam(i, b.l0, b.k );
            if(iDebug>0){  printf( "bond[%i] (%i,%i) %g %g | %g %g\n", i, ff.bond2atom[i].i, ff.bond2atom[i].j, ff.bond_l0[i], ff.bond_k[i], b.l0, b.k ); }
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            const Angle& a  = angles[i];
            ff.ang2bond[i] = a.bonds;
            ff.setAngleParam(i, a.a0, a.k );
            if(iDebug>0){  printf( "angle[%i] (%i,%i) (%g,%g) %g\n", i, ff.ang2bond[i].i, ff.ang2bond[i].j, ff.ang_cs0[i].x, ff.ang_cs0[i].y, ff.ang_k[i] ); }
        }
        for(int i=0; i<dihedrals.size(); i++){
            const Dihedral& d  = dihedrals[i];
            ff.tors2bond[i] = d.bonds;
            ff.setTorsParam( i, d.n, d.k );
            if(iDebug>0){ printf( "dihedrals[%i] (%i,%i,%i) %i %g\n", i, ff.tors2bond[i].a, ff.tors2bond[i].b, ff.tors2bond[i].c, ff.tors_n[i], ff.tors_k[i] ); }
        }
        ff.angles_bond2atom();
        ff.torsions_bond2atom();
        //exit(0);
    }
#endif // ForceField_h

void updatePBC( Vec3d* pbcShifts, Mat3d* M=0 ){
    if(M==0) M=&lvec;
    for(int i=0; i<bonds.size(); i++){
        //pbcShifts[i] = pbcShift( bondPBC[i] );
        //pbcShifts[i] = pbcShift( (Vec3i)bonds[i].ipbc );
        const Vec3i8& G = bonds[i].ipbc;
        pbcShifts[i] = M->a*G.a + M->b*G.b + M->c*G.c;
        //if( pbcShifts[i].norm2()>0.1 ){ printf( "PBC-bond[%i]  atoms(%i,%i)  pbcShift(%g,%g,%g) ipb(%i,%i,%i)\n",  bonds[i].atoms.a, bonds[i].atoms.b, pbcShifts[i].x,pbcShifts[i].y,pbcShifts[i].z, bonds[i].ipbc.x,bonds[i].ipbc.y,bonds[i].ipbc.z );  };  // DEBUG
    }
}

#ifdef MMFFsp3_h
    //void toMMFFsp3( MMFFsp3& ff, bool bRealloc=true, double K_sigma=1.0, double K_pi=1.0, double K_ecap=0.75, bool bATypes=true ){
    void toMMFFsp3( MMFFsp3& ff, bool bRealloc=true, bool bEPairs=true ){
        int nAmax = atoms.size();
        int nBmax = bonds.size();
        int nCmax = confs.size();
        // {
        //     printf("!!!! WARRNING DEBUG HACK !!!! Builder::toMMFFsp3(): change array sizes \n");
        //     printf("!!!! Before: nAmax %i nBmax %i \n", nAmax, nBmax);
        //     nAmax = frags[0].atomRange.b;
        //     nBmax = frags[0].bondRange.b;
        //     nCmax = frags[0].confRange.b;
        //     printf("!!!! After: nAmax %i nBmax %i \n", nAmax, nBmax);
        // }
        //printf("toMMFFsp3() verbosity %i \n", verbosity );
        int npi,ne; ne=countPiE( npi, 0,nCmax );
        if(!bEPairs) ne=0;
        int nconf = nCmax;
        int ncap  = nAmax - nconf;
        int nb    = bonds.size();
        if(verbosity>0)printf(  "MM:Builder::toMMFFsp3() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        if(bRealloc)ff.realloc( nconf, nb+ne, npi, ncap+ne );
        export_apos     ( ff.apos  ,0,nAmax);
        export_atypes   ( ff.atype ,0,nAmax);
        export_bonds    ( ff.bond2atom, ff.bond_l0, ff.bond_k, ff.bond_kPi,  0,nBmax);
        if ( ff.nneigh_max != N_NEIGH_MAX  ){ printf( "ERROR in MM:Builder.toMMFFsp3() N_NEIGH_MAX(%i) != ff.nneigh_max(%i) ", N_NEIGH_MAX, ff.nneigh_max ); exit(0); } 
        Vec3d hs[N_NEIGH_MAX];
        int ipi=0;
        //int ja=0;
        int ie0=nconf+ncap;
        int iie = 0;

        for(int i=0; i<ff.nnode*4;  i++){ ff.aneighs[i]=-1; };

        //for(int i=0; i<ff.nnode*ff.nneigh_max; i++){ ff.Kneighs[i]=K_sigma; }
        //for(int ia=0; ia<atoms.size(); ia++ ){  
        for(int ia=0; ia<nAmax; ia++ ){
            //printf( "atom[%i] \n", ia  );
            int ic = atoms[ia].iconf;
            //printf( "atom[%i] iconf %i \n", ia, ic  );
            if(ic>=0){
                AtomConf& conf = confs[ic];
                //if(verbosity>1) printf( "atom[%i] conf[%i] n,nb,npi,ne(%i,%i,%i,%i)[", ia, ic, conf.n,conf.nbond,conf.npi,conf.ne  );
                //printf( "atom[%i] conf[%i] n,nb,npi,ne(%i,%i,%i,%i)[ \n", ia, ic, conf.n,conf.nbond,conf.npi,conf.ne  );
                int*    ngs  = ff.aneighs + ia*N_NEIGH_MAX;
                int*    nbs  = ff.abonds  + ia*N_NEIGH_MAX;
                //double* kngs = ff.Kneighs + ia*N_NEIGH_MAX;
                // -- atoms
                for(int k=0; k<conf.nbond; k++){
                    int ib = conf.neighs[k];
                    int ja = bonds[ib].getNeighborAtom(ia);
                    hs[k]  = atoms[ja].pos - atoms[ia].pos;
                    hs[k].normalize();
                    ngs[k] = ja;
                    nbs[k] = ib;
                    //kngs[k]=K_sigma;
                }
                for(int k=conf.nbond; k<ff.nneigh_max; k++ ) { nbs[k] = -1; }
                makeConfGeom( conf.nbond, conf.npi, hs );
                //printf( "atom[%i] nb,npi,ne(%i,%i,%i) ngs{%i,%i,%i,%i} h0(%g,%g,%g) h1(%g,%g,%g) h2(%g,%g,%g) h3(%g,%g,%g) \n",  ia, conf.nbond, conf.npi, conf.ne,  conf.neighs[0],conf.neighs[1],conf.neighs[2],conf.neighs[3],    hs[0].x,hs[0].y,hs[0].z,    hs[1].x,hs[1].y,hs[1].z,    hs[2].x,hs[2].y,hs[2].z,    hs[3].x,hs[3].y,hs[3].z   );

                int npi_neigh = countAtomPiNeighs(ia);
                assignSp3Params( ff.atype[ia], conf.nbond, conf.npi, conf.ne, npi_neigh, ff.NeighParams[ia] );

                //if( (conf.nbond==2) && (conf.npi==1) ){ printf( "atom[%i](nb=%i,npi=%i,ne=%i) angles(%g,%g,%g)\n", ia, conf.nbond,conf.npi,conf.ne, hs[3].getAngle(hs[0])/M_PI, hs[3].getAngle(hs[1])/M_PI, hs[3].getAngle(hs[2])/M_PI ); }
                //for(int k=0; k<N_NEIGH_MAX; k++){
                //    if((N_NEIGH_MAX-k)<=conf.npi){ glColor3f(1.,0.,0.); }else{ glColor3f(1.,0.,0.); }
                //    Draw3D::drawVecInPos( hs[k], atoms[ia].pos );
                //}
                // pi-bonds
                for(int k=0; k<conf.npi; k++ ){
                    int ik=N_NEIGH_MAX-1-k;
                    ff.pipos[ipi] = hs[ik];
                    ngs[ik] = -ipi-2;
                    //kngs[ik]=K_pi;
                    //printf( "pi[%i] a[%i] ng[%i] h(%g,%g,%g) \n", ipi, ia, ngs[ik], ff.pipos[ipi].x,ff.pipos[ipi].y,ff.pipos[ipi].z   );
                    ipi++;
                }
                if(bEPairs){
                    // e-cap
                    int etyp=-1;  if(params) etyp=params->atomTypeDict["E"];
                    for(int k=conf.nbond; k<conf.nbond+conf.ne; k++ ){
                        int ie=ie0+iie;
                        ff.apos[ie]=atoms[ia].pos + hs[k]*0.5;
                        ff.atype[ie] = etyp; // electon-pair type
                        //kngs[k] = K_ecap;
                        int ib=nb+iie;
                        ff.bond2atom[ib]=(Vec2i){ia,ie};
                        ff.bond_l0  [ib]=0.5;
                        ff.bond_k   [ib]=defaultBond.k;
                        ngs[k] = ie;
                        nbs[k] = ib;
                        //ngs[k]=0;
                        iie++;
                    }
                }
                if(verbosity>1){ for(int k=0; k<N_NEIGH_MAX; k++ ){ printf( " %i,", ngs[k] ); }; printf( "] \n" ); }
            }
        }
        ff.ne =iie;
        ff.ie0=ie0;
        if( bPBC ){ ff.initPBC(); updatePBC( ff.pbcShifts ); }
        //printf( "check number of pi bonds ipi %i npi %i \n", ipi, npi );
        if(verbosity>0)printf(  "... MM:Builder::toMMFFsp3() DONE \n"  );
    }
#endif // MMFFmini_h

void makeNeighs( int*& aneighs, int perAtom ){
    int na = atoms.size();
    _allocIfNull( aneighs, na );
    int ntot= na*perAtom;
    for(int i=0;i<ntot; i++){ aneighs[i]=-1; }; // back neighbors
    for(int ia=0; ia<na; ia++ ){
        const Atom& A =  atoms[ia];
        if(A.iconf>=0){
            const AtomConf& conf = confs[A.iconf];
            for(int k=0; k<conf.nbond; k++){
                int ib        = conf.neighs[k];
                if(ib<0) continue;
                const Bond& B = bonds[ib];
                int ja        = B.getNeighborAtom(ia);
                aneighs[ ia*perAtom + k ] = ja;
                int jc = atoms[ja].iconf;
                //printf( "makeNeighs[%i|%i] ja %i jc %i \n", ia, k, ja, jc );
                if( jc==-1 ){ aneighs[ ja*perAtom ]=ia; }
            }
        }
    }
    //for(int ia=0; ia<na; ia++){ printf( "aneigh[%i](%i,%i,%i,%i)\n", ia, aneighs[ia*perAtom],aneighs[ia*perAtom+1],aneighs[ia*perAtom+2],aneighs[ia*perAtom+3] ); };
}

#ifdef MMFFsp3_loc_h
void toMMFFsp3_loc( MMFFsp3_loc& ff, bool bRealloc=true, bool bEPairs=true ){

        //double c0s[3]{-0.33333,-0.5,-1.0}; // cos(angle)   sp1 sp2 sp3
        double ang0s[3]{ 109.5 *M_PI/180.0, 120.0*M_PI/180.0, 180.0*M_PI/180.0 }; // cos(angle)   sp1 sp2 sp3

        int nAmax = atoms.size();
        int nCmax = confs.size();
        int npi,ne; ne=countPiE( npi, 0,nCmax );
        if(!bEPairs) ne=0;
        int nconf = nCmax;
        int ncap  = nAmax - nconf;
        int nb    = bonds.size();
        //if(verbosity>0)
        printf(  "MM:Builder::toMMFFsp3_loc() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        if(bRealloc)ff.realloc( nconf, ncap+ne );
        Vec3d hs[4];
        int ie0=nconf+ncap;
        int iie = 0;

        //params->printAtomTypeDict();
        int etyp=-1;  if(params) etyp=params->atomTypeDict["E"];
        for(int i=0; i<ff.nnode;  i++){ ff.aneighs[i]=Quat4i{-1,-1,-1,-1};  ff.bLs[i]=Quat4dZero, ff.bKs[i]=Quat4dZero, ff.Ksp[i]=Quat4dZero, ff.Kpp[i]=Quat4dZero; }; // back neighbors
        for(int ia=0; ia<nAmax; ia++ ){
            const Atom& A =  atoms[ia];
            ff.apos  [ia] = A.pos;
            ff.atypes[ia] = A.type;
            if(A.iconf>=0){

                // Prepare params and orientation
                AtomConf& conf = confs[A.iconf];
                int npi_neigh = countAtomPiNeighs(ia);
                //assignSp3Params( A.type, conf.nbond, conf.npi, conf.ne, npi_neigh, ff.NeighParams[ia] );

                // setup atom (onsite)
                // ff.apars[ia].x = c0s[conf.npi];    // ssC0  // cos(angle) for angles (sigma-siamg)
                // ff.apars[ia].y = 1.0;              // ssK   // stiffness  for angles
                // ff.apars[ia].z = 0.0;              // piC0  // stiffness  for orthogonalization sigma-pi 

                double ang0 = ang0s[conf.npi];
                ang0 *= 0.5;
                ff.apars[ia].x = cos(ang0);    // ssC0  // cos(angle) for angles (sigma-siamg)
                ff.apars[ia].y = sin(ang0);
                ff.apars[ia].z = 1.0;              // ssK   // stiffness  for angles
                ff.apars[ia].w = 0.0;              // piC0  // stiffness  for orthogonalization sigma-pi 
                //printf( "atom[%i] npi(%i)=> angle %g cs(%g,%g) \n", ia, conf.npi, ang0*180./M_PI, ff.apars[ia].x, ff.apars[ia].y  ); 

                // setup ff neighbors
                
                int*     ngs  = ff.aneighs[ia].array;
                double*  bL   = ff.bLs[ia].array;
                double*  bK   = ff.bKs[ia].array;
                double*  Ksp  = ff.Ksp[ia].array;
                double*  Kpp  = ff.Kpp[ia].array;
                //printf( "BEFOR atom[%i] ngs{%i,%i,%i,%i}\n", ia, ngs[0],ngs[1],ngs[2],ngs[3] );

                // -- atoms
                //printf( "atom[%i] ne %i \n", ia, conf.ne, conf.nbond );
                // --- Generate Bonds
                for(int k=0; k<conf.nbond; k++){
                    int ib = conf.neighs[k];
                    const Bond& B = bonds[ib];
                    int ja = B.getNeighborAtom(ia);
                    hs[k]  = atoms[ja].pos - A.pos;
                    hs[k].normalize();
                    ngs[k] = ja;
                    bL [k]=B.l0;
                    bK [k]=B.k;
                    if( (conf.npi>0)||(conf.ne>0) ){
                        Ksp[k]=Ksp_default;
                    }
                    int nej  = getAtom_ne (ja);
                    int npij = getAtom_npi(ja);
                    Kpp[k]=B.kpp;
                }
                makeConfGeom( conf.nbond, conf.npi, hs );
                if(bEPairs){ // --- Generate electron pairs
                    int ns = conf.nbond+conf.ne;
                    for(int k=conf.nbond; k<ns; k++){    
                        int ie=ie0+iie;
                        ngs [k] =ie;
                        //printf( "atom[%i|%i] ie %i \n", ia, k, ie );
                        ff.apos  [ie] = atoms[ia].pos + hs[k]*Lepair;
                        ff.atypes[ie] = etyp;
                        bK [k]=Kepair;
                        bL [k]=Lepair;
                        if( conf.npi>0 ) Ksp[k]=Ksp_default;   // only electron on atoms without pi-orbital are conjugted with pi-orbitas on neighboring atoms
                        iie++;
                    }
                }
                ff.pipos[ia] = hs[N_NEIGH_MAX-1]; // Pi orientation
                //if(verbosity>1){ for(int k=0; k<N_NEIGH_MAX; k++ ){ printf( " %i,", ngs[k] ); }; printf( "] \n" ); }
                //printf( "AFTER atom[%i] ngs{%i,%i,%i,%i}\n", ia, ngs[0],ngs[1],ngs[2],ngs[3] );
            } // if(A.iconf>=0){
        }
        ff.makeBackNeighs();
        //if( bPBC ){ ff.initPBC(); updatePBC( ff.pbcShifts ); }
        if(verbosity>0)printf(  "... MM:Builder::toMMFFsp3_loc() DONE \n"  );
    }
#endif // MMFFf4_h

#ifdef MMFFf4_h
void toMMFFf4( MMFFf4& ff,  bool bRealloc=true, bool bEPairs=true ){

        //float c0s[3]{-0.33333,-0.5,-1.0}; // cos(angle)   sp1 sp2 sp3
        float ang0s[3]{ 109.5 *M_PI/180.0, 120.0*M_PI/180.0, 180.0*M_PI/180.0 }; // cos(angle)   sp1 sp2 sp3

        int nAmax = atoms.size();
        int nCmax = confs.size();
        int npi,ne; ne=countPiE( npi, 0,nCmax );
        if(!bEPairs) ne=0;
        int nconf = nCmax;
        int ncap  = nAmax - nconf;
        int nb    = bonds.size();
        if(verbosity>0)printf(  "MM:Builder::toMMFFf4() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        if(bRealloc)ff.realloc( nconf, ncap+ne );
        Vec3d hs[N_NEIGH_MAX];
        int ie0=nconf+ncap;
        int iie = 0;

        for(int i=0; i<ff.nnode;  i++){ ff.aneighs[i]=Quat4i{-1,-1,-1,-1};  ff.bLs[i]=Quat4fZero, ff.bKs[i]=Quat4fZero, ff.Ksp[i]=Quat4fZero, ff.Kpp[i]=Quat4fZero; }; // back neighbors

        for(int ia=0; ia<nAmax; ia++ ){
            const Atom& A =  atoms[ia];
            ff.apos[ia].f = (Vec3f) A.pos;
            ff.apos[ia].e = 0;
            if(A.iconf>=0){

                // Prepare params and orientation
                AtomConf& conf = confs[A.iconf];
                int npi_neigh = countAtomPiNeighs(ia);
                //assignSp3Params( A.type, conf.nbond, conf.npi, conf.ne, npi_neigh, ff.NeighParams[ia] );

                // setup atom (onsite)
                // ff.apars[ia].x = c0s[conf.npi];    // ssC0  // cos(angle) for angles (sigma-siamg)
                // ff.apars[ia].y = 1.0;              // ssK   // stiffness  for angles
                // ff.apars[ia].z = 0.0;              // piC0  // stiffness  for orthogonalization sigma-pi 
                // ff.apars[ia].w = 0.0; 

                double ang0 = ang0s[conf.npi];
                ang0 *= 0.5;
                ff.apars[ia].x = cos(ang0);    // ssCos0  // cos(angle) for angles (sigma-siamg)
                ff.apars[ia].y = sin(ang0);    // ssSin0
                ff.apars[ia].z = 1.0;          // ssK     // stiffness  for angles
                ff.apars[ia].w = 0.0;          // piCos0  // stiffness  for orthogonalization sigma-pi 
                //printf( "atom[%i] npi(%i)=> angle %g cs(%g,%g) \n", ia, conf.npi, ang0*180./M_PI, ff.apars[ia].x, ff.apars[ia].y  ); 

                // setup ff neighbors
                
                int*    ngs  = ff.aneighs[ia].array;
                float*  bL   = ff.bLs[ia].array;
                float*  bK   = ff.bKs[ia].array;
                float*  Ksp  = ff.Ksp[ia].array;
                float*  Kpp  = ff.Kpp[ia].array;
                // -- atoms
                //printf( "atom[%i] ne %i \n", ia, conf.ne, conf.nbond );
                // --- Generate Bonds
                for(int k=0; k<conf.nbond; k++){
                    int ib        = conf.neighs[k];
                    const Bond& B = bonds[ib];
                    int ja = B.getNeighborAtom(ia);
                    hs[k] = atoms[ja].pos - A.pos;
                    hs[k].normalize();
                    ngs[k] = ja;
                    bL [k]=B.l0;
                    bK [k]=B.k;
                    if( (conf.npi>0)||(conf.ne>0) ){
                        Ksp[k]=Ksp_default;
                    }
                    int nej  = getAtom_ne (ja);
                    int npij = getAtom_npi(ja);
                    Kpp[k]=B.kpp;
                }
                makeConfGeom( conf.nbond, conf.npi, hs );
                if(bEPairs){  // --- Generate electron pairs
                    int ns = conf.nbond+conf.ne;
                    for(int k=conf.nbond; k<ns; k++){    
                        int ie=ie0+iie;
                        ngs [k] =ie;
                        //printf( "atom[%i|%i] ie %i \n", ia, k, ie );
                        ff.apos[ie].f = (Vec3f)( atoms[ia].pos + hs[k]*Lepair );
                        ff.apos[ie].e = 0;
                        bK [k]=Kepair;
                        bL [k]=Lepair;
                        if( conf.npi>0 ) Ksp[k]=Ksp_default;   // only electron on atoms without pi-orbital are conjugted with pi-orbitas on neighboring atoms
                        iie++;
                    }
                }
                ff.pipos[ia].f = (Vec3f)hs[N_NEIGH_MAX-1]; // Pi orientation
                ff.pipos[ia].e = 0;
                //if(verbosity>1){ for(int k=0; k<N_NEIGH_MAX; k++ ){ printf( " %i,", ngs[k] ); }; printf( "] \n" ); }
                //printf( "AFTER atom[%i] ngs{%i,%i,%i,%i}\n", ia, ngs[0],ngs[1],ngs[2],ngs[3] );
            } // if(A.iconf>=0){
        }
        ff.makeBackNeighs();
        if(verbosity>0)printf(  "... MM:Builder::toMMFFf4() DONE \n"  );
    }
#endif // MMFFf4_h

#ifdef MMFFmini_h
    void toMMFFmini( MMFFmini& ff, const MMFFparams* params ){
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        export_apos     ( ff.apos );
        export_bonds    ( ff.bond2atom,   ff.bond_l0, ff.bond_k );
        export_angles   ( ff.ang2bond, 0, ff.ang_cs0, ff.ang_k  );
        export_dihedrals( ff.tors2bond,   ff.tors_n,  ff.tors_k );

        if( bPBC ){ ff.initPBC(); updatePBC( ff.pbcShifts ); }

        ff.angles_bond2atom();
        ff.torsions_bond2atom();
    }

    // ----  OLD version ----
    /*
    void toMMFFmini( MMFFmini& ff, const MMFFparams* params ){
        //printf( "na %i nb %i nA %i \n", atoms.size(), bonds.size(), dihedrals.size() );
        //mmff->deallocate();
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        for(int i=0; i<atoms.size(); i++){
            ff.apos [i]  = atoms[i].pos;
            //println(atoms[i]);
            //if( atoms[i].iconf>=0 ) println(confs[atoms[i].iconf]);
        }
        for(int i=0; i<bonds.size(); i++){
            const Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            ff.bond2atom[i]    = ib;
            if(params){
                params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, ff.bond_l0[i], ff.bond_k[i] );
            }else{
                //printf( "no params \n" );
                ff.setBondParam(i, b.l0, b.k );
            }
            //printf( "bond[%i] (%i,%i) %g %g | %g %g\n", i, ff.bond2atom[i].i, ff.bond2atom[i].j, ff.bond_l0[i], ff.bond_k[i], b.l0, b.k );
            //bondTypes[i]       = bonds[i].type;
        }
        //printf( "toMMFFmini . Angles.size() %i \n", angles.size() );
        for(int i=0; i<angles.size(); i++){
            const Angle& a  = angles[i];
            ff.ang2bond[i] = a.bonds;
            ff.setAngleParam(i, a.a0, a.k );
            //printf( "angle[%i] (%i,%i) (%g,%g) %g\n", i, ff.ang2bond[i].i, ff.ang2bond[i].j, ff.ang_cs0[i].x, ff.ang_cs0[i].y, ff.ang_k[i] );
        }
        for(int i=0; i<dihedrals.size(); i++){
            const Dihedral& d  = dihedrals[i];
            ff.tors2bond[i] = d.bonds;
            ff.setTorsParam( i, d.n, d.k );
            //printf( "dihedrals[%i] (%i,%i,%i) %i %g\n", i, ff.tors2bond[i].a, ff.tors2bond[i].b, ff.tors2bond[i].c, ff.tors_n[i], ff.tors_k[i] );
        }
        ff.angles_bond2atom();
        ff.torsions_bond2atom();
        //exit(0);
    }
    */
#endif // MMFFmini_h

#ifdef MMFF_h
    void toMMFF( MMFF * mmff, MMFFparams* params ){
        //mmff->deallocate();
        mmff->allocate( atoms.size(), bonds.size(), angles.size(), 0 );
        //int * atomTypes = new int[atoms.size()];
        //int * bondTypes = new int[bonds.size()];
        for(int i=0; i<atoms.size(); i++){
            mmff->atypes[i] = atoms[i].type;
            mmff->atom2frag[i] = atoms[i].frag;
            mmff->apos [i]  = atoms[i].pos;
            mmff->aREQ [i]  = atoms[i].REQ;
            //atomTypes[i]  = atoms[i].type;
            //printf( "iatom %i atype %i ifrag %i pos (%g,%g,%g) REQ (%g,%g,%g) \n", i, atoms[i].type, atoms[i].frag, atoms[i].pos.x,atoms[i].pos.y,atoms[i].pos.z, atoms[i].REQ.x,atoms[i].REQ.y,atoms[i].REQ.z );
        }
        for(int i=0; i<bonds.size(); i++){
            mmff->bond2atom[i] = bonds[i].atoms;
            Vec2i ib           = bonds[i].atoms;
            params->getBondParams( atoms[ib.x].type, atoms[ib.y].type, bonds[i].type, mmff->bond_0[i], mmff->bond_k[i] );
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            mmff->ang2bond[i] = angles[i].bonds;
            mmff->ang_0[i] = {1.0,0.0}; // TODO FIXME
            mmff->ang_k[i] = 0.5;       // TODO FIXME
        }
        if( frags.size()>0 ){
            mmff->allocFragment( frags.size() );
            for(int i=0; i<frags.size(); i++){
                MM::Fragment& fragi = frags[i];
                mmff->frag2a  [i] = fragi.atomRange.x;
                mmff->fragNa  [i] = fragi.atomRange.y;
                //mmff->fapos0s [i] = fragi.mol->pos;
                double * posi= (mmff->poses + i*8);
                *(Vec3d *)(posi  )= fragi.pos;
                *(Quat4d*)(posi+4)= fragi.rot;
            }
        }
        //params.fillBondParams( mmff->nbonds, mmff->bond2atom, bondTypes, atomTypes, mmff->bond_0, mmff->bond_k );
        //delete [] atomTypes;
        //delete [] bondTypes;
    }
#endif // MMFF_h

#ifdef EFF_h
//#if defined(EFF_h) || defined(EFF_old_h)
    void toEFF( EFF& ff, const EFFAtomType* params, double esize, double dpair ){
        //int ne = bonds.size() * 2; // ToDo
        int ne = 0;
        int na = 0;
        DEBUG
        for(int i=0; i<atoms.size(); i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            na++;
            ne += params[ ityp ].ne;
            printf( "[%i] ityp %i ne %i  %i \n", i, ityp, params[ityp].ne, ne );
        }
        DEBUG
        printf( "na %i ne %i | %i \n", atoms.size(), ne, bonds.size()*2 );
        ff.realloc( na, ne );
        DEBUG
        for(int i=0; i<na; i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            //printf( "[%i] ityp %i \n" );
            ff.apos [i]  = atoms[i].pos;
            //ff.aQ [i]  = params[ ityp ].ne; // ToDo
            ff.aPars[i]  = EFF::default_AtomParams[ityp];
        }
        DEBUG
        for(int i=0; i<bonds.size(); i++){
            const MM::Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            double c1=0.5-dpair;
            double c2=1-c1;
            int i2 = i*2;
            ff.epos [i2] = atoms[ib.a].pos*c1 + atoms[ib.b].pos*c2;
            ff.esize[i2] = esize;
            ff.espin[i2] = 1;
            i2++;
            ff.epos [i2] = atoms[ib.a].pos*c2 + atoms[ib.b].pos*c1;
            ff.esize[i2] = esize;
            ff.espin[i2] = -1;
            //break;
        }
        DEBUG
        // ToDo:  pi-bonds & e-pairs

    }

    /*
    int countValenceElectrons(){
        int ne=0;
        for(int i=0; i<atoms.size(); i++){

        }
    }
    */
#endif  // EFF_h


#ifdef EFF_old_h
    void toEFF_old( EFF& ff, const EFFAtomType* params, double esize, double dpair ){
        //int ne = bonds.size() * 2; // ToDo
        int ne = 0;
        int na = 0;
        DEBUG
        for(int i=0; i<atoms.size(); i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            na++;
            ne += params[ ityp ].ne;
            printf( "[%i] ityp %i ne %i  %i \n", i, ityp, params[ityp].ne, ne );
        }
        DEBUG
        printf( "na %i ne %i | %i \n", atoms.size(), ne, bonds.size()*2 );
        ff.realloc( na, ne );
        DEBUG
        for(int i=0; i<na; i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            //printf( "[%i] ityp %i \n" );
            ff.apos [i]  = atoms[i].pos;
            ff.aQ   [i]  = params[ ityp ].ne; // ToDo
        }
        DEBUG
        for(int i=0; i<bonds.size(); i++){
            const MM::Bond& b  = bonds[i];
            const Vec2i& ib    = b.atoms;
            double c1=0.5-dpair;
            double c2=1-c1;
            int i2 = i*2;
            ff.epos [i2] = atoms[ib.a].pos*c1 + atoms[ib.b].pos*c2;
            ff.esize[i2] = esize;
            ff.espin[i2] = 1;
            i2++;
            ff.epos [i2] = atoms[ib.a].pos*c2 + atoms[ib.b].pos*c1;
            ff.esize[i2] = esize;
            ff.espin[i2] = -1;
            //break;
        }
        DEBUG
        // ToDo:  pi-bonds & e-pairs

    }

    /*
    int countValenceElectrons(){
        int ne=0;
        for(int i=0; i<atoms.size(); i++){

        }
    }
    */
#endif  // EFF_old_h

}; // MMFFBuilder


} // namespace MMFF

#endif // MMFFBuilder_h
