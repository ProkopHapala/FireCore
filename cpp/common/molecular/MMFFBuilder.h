#ifndef MMFFBuilder_h
#define MMFFBuilder_h
/// @file MMFFBuilder.h   @brief Classes for building and editing molecular topology and building instances of force-fields for paticular molecule or system 
/// @ingroup Classical_Molecular_Mechanics

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>

#include  "globals.h"

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


//#include "MolecularGraph.h"
#include "LimitedGraph.h"


#include "Kekule.h"


// =============== Structs for Atoms, Bonds etc...

/// @namespace MM       @brief Molecular-Mechanics namespace
namespace MM{

static const double const_eVA2_Nm = 16.02176634;

// ============================
// ========= Atom
// ============================

/// @brief Atom context in Molecular-Mechanics  structure with 3D coordinates, type, neighborlist etc.
struct Atom{
    constexpr const static Quat4d HcapREQ    = Quat4d{ 1.4870, 0.026095977, 0., 0. }; // sqrt(0.000681)=0.026095977
    constexpr const static Quat4d defaultREQ = Quat4d{ 1.7,    0.061067605, 0., 0. }; // sqrt(0.0037292524)=0.061067605
    int id;
    int type;    ///< atomic type
    int frag;    ///< to which fragment it belongs ?
    int iconf;   ///< index of the configuration object keeping the list of bonded neighbors
    Vec3d  pos;  ///< atomic coordinates
    Quat4d REQ;  /// non-covaloent interaction parameters {RvdW,EvdW,Q}              constexpr Vec3d{1.7,sqrt(0.0037292524),0}
    //Atom() = default;

    void print()const{ printf( " Atom{id %i t %i c %i f %i REQ(%g,%g,%g,%g) pos(%g,%g,%g)}", id, type, iconf, frag, REQ.x, REQ.y, REQ.z,REQ.w, pos.x,pos.y,pos.z ); }

    Atom() = default;
    Atom(const Vec3d& pos_):type{0},frag{-1},iconf{-1},REQ{defaultREQ},pos{pos_}{};
    Atom(const Vec3d& pos_,const Quat4d& REQ_):type{0},frag{-1},iconf{-1},REQ{REQ_},pos{pos_}{};
    Atom(int type_,int frag_,int iconf_,const Vec3d& pos_,const Quat4d& REQ_):id(-1),type{type_},frag{frag_},iconf{iconf_},REQ{REQ_},pos{pos_}{};
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

/// @brief Atom configurations contains information about bonding topology of the atoms (list of neighbors, number of sigma-bonds, pi-orbitas and free electron pairs)
struct AtomConf{

    int iatom=-1;
    uint8_t n     =0; ///< total sum, should be equal to valency of the atom
    uint8_t nbond =0; ///< number of sigma bonds
    uint8_t npi   =0; ///< number of pi-orbitals (bonds)
    uint8_t ne    =0; ///< number of electron pairs
    uint8_t nH    =0; ///< number of capping atoms (e.g. hydrogen)
    int neighs[N_NEIGH_MAX]; // neighs  - NOTE: bonds not atoms !!!!

    Vec3d  pi_dir{0.,0.,0.};

    //AtomConf() = default;

    // compute the number of pi bonds
    inline int fillIn_npi_sp3(){ npi=4-nbond-ne-nH; n=4; return npi; };

    // cleanup
    inline int clean_pi_sp3  (){ npi=0; n=nbond+ne+nH;   return n;   };

    // shifting all bond indexes by dib
    inline void rebase(int dib){
        for(int i=0; i<N_NEIGH_MAX; i++){ if(neighs[i]>=0)neighs[i]+=dib; };
    }

    // given a bond, it finds if the atom is its neighbor
    inline int findNeigh(int ib)const{
        for(int i=0; i<N_NEIGH_MAX; i++){ if(neighs[i]==ib) return i;
        }
        return -1;
    }

    // add a bond (or pi bond, or electron pair, or capping atom) into an atom neighbor list
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

    inline int countBonds(){ nbond=0; for(int i=0; i<N_NEIGH_MAX;i++){ if(neighs[i]>=0)nbond++; } return nbond; }

    inline int findMinBondIndex(int i0){
        int vmin=1000000;
        int imin=-1;
        for(int i=i0; i<N_NEIGH_MAX; i++){
            int ngi=neighs[i]; 
            if( (ngi>=0)&&(ngi<vmin) ){ imin=i; vmin=ngi; };
        }
        return imin;
    }

    inline bool sortBonds(){ // bouble sort
        bool change=false;
        for (int i=0; i<N_NEIGH_MAX-1; i++) {  // Bouble-sort
            int imin = findMinBondIndex(i+1);
            //printf( "sortBonds()[%i] imin=%i ngs{%i,%i,%i,%i}\n", i, imin, neighs[0],neighs[1],neighs[2],neighs[3] );
            if( imin<0 ) break;
            int ngi=neighs[i];
            if( (ngi<0)|| (ngi>neighs[imin]) ){ _swap(neighs[i],neighs[imin]); };
        }
        //printf( "sortBonds() DONE ngs{%i,%i,%i,%i}\n", neighs[0],neighs[1],neighs[2],neighs[3] );
        return change;
    }


    inline bool addBond (int i){ return addNeigh(i,nbond); };
    inline bool addH    (     ){ return addNeigh((int)NeighType::H    ,nH ); };
    inline bool addPi   (     ){ return addNeigh((int)NeighType::pi   ,npi); };
    inline bool addEpair(     ){ return addNeigh((int)NeighType::epair,ne ); };
    inline int  updateNtot  (){ n=nbond+npi+ne+nH; return n; };
    inline void updateNeighs(){ sortBonds(); countBonds(); updateNtot(); };

    inline void clearNonBond(){ n=nbond; npi=0;ne=0;nH=0; };
    inline void clearBond   (){ nbond=0; updateNtot();     };
    inline void setNonBond(int npi_,int ne_){ npi=npi_; ne=ne_; updateNtot();  }
    inline void init0(){ for(int i=0; i<N_NEIGH_MAX; i++)neighs[i]=-1; nbond=0; clearNonBond(); }

    //void print()const{ printf( " AtomConf{ ia %i, n %i nb %i np %i ne %i nH %i (%i,%i,%i,%i) }", iatom, n, nbond, npi, ne, nH , neighs[0],neighs[1],neighs[2],neighs[3] ); }
    inline void print()const{ printf( " AtomConf{ia %i, n,nb,np,ne,nH(%i,%i,%i,%i,%i) [%i,%i,%i,%i]}", iatom, n,nbond,npi,ne,nH, neighs[0],neighs[1],neighs[2],neighs[3] ); }

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
    AtomConf(int iatom_,int npi_,int ne_):iatom{iatom_},npi{npi_},ne{ne_},nH{0},nbond{0},n{npi_+ne_},pi_dir{Vec3dZero}{ for(int i=0;i<N_NEIGH_MAX;i++)neighs[i]=-1; };
    //AtomConf(const MMFFAtomConf&) = default;
    //AtomConf(std::initializer_list<MMFFAtomConf>) {};
};

// ============================
// ========= Bond
// ============================
/// @brief Bond store information about atom and bond indexes linked by the bond and bond parameters (equlibirum lengh, stiffness, pi-stiffness), and index of replica in periodic-boundary conditions
struct Bond{
    // --- this breaks {<brace-enclosed initializer list>} in C++11
    int    type  = -1;             ///< type of the bond ( e.g. single, double, triple )
    Vec2i  atoms = (Vec2i){-1,-1}; ///< indexes of atoms linked by the bond
    double l0=1.0; ///< equlibirum bond length
    double k=0.0;  ///< bond stiffness
    double kpp=0;  ///< pi-pi stiffness
    Vec3i8 ipbc=Vec3i8{0,0,0}; ///< index of cell image in periodic boundary conditions
    double order=0.0; ///< actual bond order (unlike type can be non-integer, e.g. for aromatic bonds)
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
    Bond(int type_, Vec2i atoms_, double l0_, double k_, double kpp_=0.0):type(type_),atoms(atoms_),l0(l0_),k(k_),ipbc{0,0,0},kpp{kpp_}{};
};

// ============================
// ========= Angle
// ============================
/// @brief Angle store information about atom indexes linked by the angle and angle parameters (equlibirum angle, stiffness)
struct Angle{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec2i  bonds = (Vec2i){-1,-1};
    //double a0    = 0;
    //double k     = 0;
    //Angle()=default;

    int type;
    Vec2i  bonds;        ///< indexes of the bonds linked by the angle
    double a0;           ///< equeilibrium angle
    double k;            ///< stiffness 
    double C0,C1,C2,C3;  ///< and coeffiicients of expansion of the angle potential
    Vec3i  atoms;        ///< indexes of the atoms involved in the angle

    void print()const{ printf( " Angle{t %i b(%i,%i) a0 %g k %g}", type, bonds.i, bonds.j, a0, k ); }

    Angle()=default;
    Angle( int type_, Vec2i bonds_, double a0_, double k_):type(type_), bonds(bonds_),a0(a0_),k(k_){ };
};

// ============================
// ========= Dihedral
// ============================
/// @brief Dihedral angle store information about atom and bond indexes linked by the dihedral angle and dihedral angle parameters (equlibirum angle, stiffness)
struct Dihedral{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec3i  bonds = (Vec3i){-1,-1,-1};
    //int    n=0;
    //double k=0;

    int    type;
    Vec3i  bonds;   ///< indexes of the bonds linked by the angle
    Quat4i atoms;   ///< indexes of the atoms involved in the angle
    int    d,n;    /// multiplicity and number of terms in the dihedral potential
    double k;    ///< stiffness 
    double a0;   ///< equeilibrium angle

    //Dihedral()=default;

    void print()const{ printf( " Dihedral{t %i b(%i,%i,%i) k %g d %i n %i}", type, bonds.a, bonds.b,bonds.c, k, d, n ); }

    Dihedral()=default;
    Dihedral( int type_, Vec3i  bonds_, int n_, double k_, int d_):type(type_), bonds(bonds_), n(n_), k(k_), d(d_){};
    Dihedral( int type_, Quat4i atoms_, int n_, double k_, int d_):type(type_), atoms(atoms_), n(n_), k(k_), d(d_){};
};

// ============================
// ========= Inversion
// ============================
/// @brief Inproper dihedral angle store information about atom and bond indexes linked by the dihedral angle and dihedral angle parameters (equlibirum angle, stiffness)
struct Inversion{

    // --- this breaks {<brace-enclosed initializer list>} in C++11
    //int type     = -1;
    //Vec3i  bonds = (Vec3i){-1,-1,-1};
    //int    n=0;
    //double k=0;

    int    type;
    Vec3i  bonds;   ///< indexes of the bonds linked by the angle
    Quat4i atoms;   ///< indexes of the atoms involved in the angle
    double k;        ///< stiffness 
    double C0,C1,C2; ///< and coeffiicients of expansion of the angle potential

    void print()const{ printf( " Inversion{t %i b(%i,%i,%i) k %g C0 %g C1 %g C2 %g}", type, bonds.a, bonds.b,bonds.c, k, C0, C1, C2 ); }

    Inversion()=default;
    Inversion( int type_, Vec3i  bonds_, double k_, double C0_, double C1_, double C2_):type(type_), bonds(bonds_), k(k_), C0(C0_), C1(C1_), C2(C2_){};
    Inversion( int type_, Quat4i atoms_, double k_, double C0_, double C1_, double C2_):type(type_), atoms(atoms_), k(k_), C0(C0_), C1(C1_), C2(C2_){};
};

// ============================
// ========= Fragment
// ============================
/// @brief Fragment groups range of atoms and bonds (and angles, dihedrals, inversions) that belong e.g. to the same molecule, or to the same residue. Useful especially for rigid-body dynamics and Bound-Boxes acceleration of non-covalent interactions
struct Fragment{
    int   imolType;
    Vec2i atomRange;
    Vec2i confRange;
    Vec2i bondRange;
    Vec2i angRange;
    Vec2i dihRange;

    Vec3d  pos;
    Quat4d rot;
    Vec3d  * pos0s;
    //Molecule * mol;     // ToDo : we can replace this by MolID to leave dependence on Molecule_h
    uint32_t color;


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
#include  "globals.h"
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
/// @brief Comprehensive builder for for building and editing molecular topology and building instances of force-fields for paticular molecule or system.
/// @details It implements many algorithms for automatic finding of bonding topology, assignment of atomic types, adding/removing explicit free electron paris and capping atoms, substitution, segmentation to groups, etc.
class Builder{  public:
    //int verbosity = 0;
    //bool bDebug = false;

    std::unordered_set<int> selection;
    std::vector<int>        atom_permut;
    std::vector<int>        atom2group;
    int ngroups = 0;


    //static int iDebug = 0;
    std::vector<Atom>       atoms;
    std::vector<Bond>       bonds;
    std::vector<Angle>      angles;
    std::vector<Dihedral>   dihedrals;
    std::vector<Inversion>  inversions;
    std::vector<AtomConf>   confs;
    //std::vector<int>  atom_neighs;

    bool bPBC  = false;
    Mat3d lvec = Mat3dIdentity;
    //Mat3d lvec = Mat3dIdentity*1000.0;  // lattice vectors for PBC (periodic boundary conditions)
    //std::vector<Vec3i> bondPBC;

    MMFFparams* params = 0;  // Each function which needs this can take it as parameter
    std::vector       <std::string>*     atomTypeNames = 0;
    std::unordered_map<std::string,int>* atomTypeDict  = 0;

    std::unordered_set<int> capping_types;
    std::unordered_set<int> sp3types;

    std::vector<Fragment>   frags;

    // ToDo: we should rather use /home/prokop/git/FireCore/cpp/common/molecular/MolecularGraph.h
    //std::vector<int> graph_color; 
    //std::vector<int> graph_dist;

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
    Quat4d defaultREQ  {  1.5, 0.0, 0.0, 0.0    };
    Bond   defaultBond { -1, {-1,-1}, 1.5, 10.0 };
    Angle defaultAngle { -1, {-1,-1}, M_PI, 0.5 };
    //Angle defaultAngle{ -1, {-1,-1}, 0.0, 0.5 };

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
    bool bAutoTypes  = true;
    bool bAddCaps    = true;
    bool bAddECaps   = false;

    // =================== Functions =====================

    void randomFragmentCollors(){
        printf("Builder::randomFragmentCollors()\n");
        srand( 1234 );
        for(int i=0; i<frags.size(); i++){
            //frags[i].color = hash_Wang( 55 + i )&0xFF  + ((hash_Wang( 4487 + i*12 )&0xFF)<<8)+ ((hash_Wang( 15455 + i*123 )&0xFF)<<16);
            //int r = rand()&0xFF;
            //int g = rand()&0xFF;
            //int b = rand()&0xFF;
            int r = (rand()>>23)&0xFF;
            int g = (rand()>>23)&0xFF;
            int b = (rand()>>23)&0xFF;
            //frags[i].color = ( rand()&0xFF  + ((rand()&0xFF)<<8)+ ((rand()&0xFF)<<16)) | 0xFF000000;
            frags[i].color = ( r + (g<<8) + (b<<16) ) | 0xFF000000;
            printf("frags[%i].color %X (%i,%i,%i)\n", i, frags[i].color, r, g, b );
        }
    }

    int addCappingTypesByIz( int iZ ){
        //printf( "Builder::addCappingTypesByIz()\n" );
        int n=0; 
        for( int i=0; i<params->atypes.size(); i++ ){ 
            if(params->atypes[i].iZ==iZ){ 
                capping_types.insert(i);
                //printf( "Builder::addCappingTypesByIz()[%i] `%s`\n", i, params->atypes[i].name );
                n++;
            } 
        } 
        return n; 
    }

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
    void move_atoms     ( Vec3d dshift,                                                                        int i0=0, int imax=-1 ){ if(imax<0){imax=atoms.size();} for(int i=i0; i<imax; i++){ atoms[i].pos.add(dshift); } }
    void transform_atoms( Mat3d M,                         Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int imax=-1 ){ if(imax<0){imax=atoms.size();} for(int i=i0; i<imax; i++){ Vec3d p; M.dot_to( atoms[i].pos-orig_old, p); p.add(orig_new); atoms[i].pos=p; } }
    void rotate_atoms   ( double angle, Vec3d axis=Vec3dZ, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int imax=-1 ){ Mat3d M; M.fromRotation(angle,axis);    transform_atoms( M,orig_old,orig_new,i0,imax); }
    void orient_atoms   ( Vec3d fw, Vec3d up,              Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int imax=-1 ){ Mat3d M; M.fromDirUp(fw,up); transform_atoms( M,orig_old,orig_new,i0,imax); }

    void changeCell( const Mat3d& lvs, Vec3d orig_old=Vec3dZero, Vec3d orig_new=Vec3dZero, int i0=0, int n=-1 ){
        Mat3d M,MM; 
        //lvec.invert_to(M); 
        lvec.invert_T_to(M); 
        //MM.set_mmul_TN(lvec,M);
        MM.set_mmul_TN(lvs,M);
        //MM = M;
        //MM.set_mmul(lvec,M);
        //MM.set_mmul(lvs,M);
        //printf("Debug changeCell()  lvec\n"); lvec .print();
        //printf("Debug changeCell()  lvs\n"); lvs   .print();
        //printf("Debug changeCell()  M (inv(lvs))\n"); M .print();
        //printf("Debug changeCell() MM\n"); MM.print(); //exit(0);
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


    int findHighestValence(){
        int imax=-1;
        int nmax=-1;
        for(int i=0; i<atoms.size(); i++){
            //const AtomConf* c = getAtomConf(int ia);
            int ic = atoms[i].iconf;
            if(ic<0)continue;
            int nb = confs[ic].nbond;
            if(nb>nmax){ nmax=nb; imax=i; } 
            
        }
        return imax;
    }

    /*

    // ToDo: we should rather use /home/prokop/git/FireCore/cpp/common/molecular/MolecularGraph.h

    int branch_Graph( int i, int oi=-1 ){}

    int walk_graph( int i, int oi ){
        int ic = atoms[i].iconf;
        if(ic){ ic<; };
    }

    void findBridges( int ia0=-1 ){
        int nc = confs.size();
        int na = atoms.size();
        color   .resize(nc);
        distance.resize(nc);
        for( int& c: color ){ c=-1; }

        if(ia<0){ ia=findHighestValence(); }

        branch_Graph( ia );


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
        atoms.clear();
        confs.clear();
        bonds.clear(); 
        //bondPBC.clear();
        angles.clear();
        dihedrals.clear();
#ifdef Molecule_h
        //mols .clear();
        frags.clear();
        fragTypes.clear();
#endif // Molecule_h
        //printf("MM::Builder::clear() DONE\n");
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

    Vec2i getBondAtomTypes(int ib){
        const MM::Bond& b = bonds[ib];
        return Vec2i{ atoms[b.atoms.a].type, atoms[b.atoms.b].type };
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

    int tryAddConfToAtom( int ia ){
        int t = atoms[ia].type;
        //printf( "Builder::tryAddConfsToAtoms()[ia=%i] t %i capping_types.contains(t)=%i \n", ia, t, capping_types.contains(t) );
        if( capping_types.contains( t ) ) return -1;
        if( atoms[ia].iconf < 0 ){
            addConfToAtom( ia, 0 );
            return confs.size()-1;
        }
        return -1;
    }

    int tryAddConfsToAtoms( int i0=0, int imax=-1 ){
        if(imax<0){ imax=atoms.size(); }
        int n=0;
        for(int ia=0; ia<imax; ia++){
            // int t = atoms[ia].type;
            // //printf( "Builder::tryAddConfsToAtoms()[ia=%i] t %i capping_types.contains(t)=%i \n", ia, t, capping_types.contains(t) );
            // if( capping_types.contains( t ) ) continue;
            // if( atoms[ia].iconf < 0 ){
            //     addConfToAtom( ia, 0 );
            //     n++;
            // }
            tryAddConfToAtom( ia );
        }
        return n;
    }

    AtomConf* insertAtom(const Atom& atom, const AtomConf* conf ){
        atoms.push_back(atom);
        return addConfToAtom( atoms.size()-1, conf );
    }
    int insertAtom(Atom& atom ){ int ia=atoms.size(); atom.id=ia; atoms.push_back(atom); return ia; }

    int insertAtom( int ityp, const Vec3d& pos, const Quat4d* REQ=0, int npi=-1, int ne=0 ){
        //printf( "insertAtom ityp %i pos( %g %g %g ) npi=%i ne=%i \n", ityp, pos.x, pos.y, pos.z, npi, ne );
        Quat4d REQloc;
        if(REQ==0)
            if(params){ 
                AtomType& at = params->atypes[ityp];
                REQloc.x=at.RvdW;
                REQloc.y=at.EvdW;
                REQloc.z=0;
                REQ=&REQloc;
            }else{ REQ=&defaultREQ; }
        int iconf=-1;
        if(npi>=0){ if( !capping_types.contains(ityp)){ 
            //printf( "insertAtom npi>0 => make Conf \n" );
            iconf=confs.size();
            confs.push_back( AtomConf(atoms.size(), npi, ne ) );
        }}
        atoms.push_back( Atom    ( ityp,-1,iconf, pos, *REQ )  );
        //printf( "insertAtom[%i]", atoms.size() ); atoms.back().print(); puts("");
        return atoms.size()-1;
    }

    int insertAtom( int ityp, const Vec3d* pos=0, const Quat4d* REQ=0, int npi=-1, int ne=0 ){
        if(pos==0)pos=&Vec3dZero;
        int ia = insertAtom( ityp, *pos, REQ, npi, ne );
        //if(bConf) addConfToAtom( ia, 0 );
        return ia;
    }
    int insertAtom( std::string name, const Vec3d* pos=0, const Quat4d* REQ=0, int npi=-1, int ne=0){ int ityp = params->atomTypeDict.at(name); return insertAtom(ityp,pos,REQ,npi,ne); };

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
            //printf( "MM::Builder.addBondToAtomConf ia %i ib %i ic %i \n", ia, ib, ic );
            bool success = confs[ic].addBond(ib);
            //printf( "MM::Builder.addBondToAtomConf ia %i ib %i success %i \n", ia, ib, success );
            if(!success){
                printf( "MM::Builder.addBondToAtomConf ia %i confs[%i].addBond(%i) success=%i (failed) \n", ia, ic, ib, success );
                //printf("ERROR: in confs[%i].addBond(%i) => exit \n", ic, ib); 
                printf("confs[%i]: ", ic ); confs[ic].print(); printf("\n");
                printBondsOfAtom( ia );
                int it1 = atoms[bonds[ib].atoms.a].type;
                int it2 = atoms[bonds[ib].atoms.b].type;
                printf( "bond(%i-%i) %s-%s \n", bonds[ib].atoms.a, bonds[ib].atoms.b, params->atypes[it1].name, params->atypes[it2].name );
                printBonds();
                exit(0); 
            }
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

    int insertBond(const Bond& bond, bool bCheck=false ){
        int ib = bonds.size();
        bonds.push_back(bond);
        //printf( "insertBond[%i]", bonds.size() ); bonds.back().print(); puts("");
        tryAddBondToAtomConf( ib, bond.atoms.i, bCheck );
        tryAddBondToAtomConf( ib, bond.atoms.j, bCheck );
        return ib;
    }

    int insertBond( Vec2i ias, int order, bool bCheck=false ){
        double k=0,l0=1.0;
        if(params){
            int iat=atoms[ias.i].type;
            int jat=atoms[ias.j].type;
            params->getBondParams( iat, jat, order, l0, k );
        }
        return insertBond( Bond(order, ias, l0, k), bCheck );
    };

    Vec2d assignBondParamsUFF( int ib ){
        Bond& b = bonds[ib];
        const Atom& ai = atoms[b.atoms.i];
        const Atom& aj = atoms[b.atoms.j];
        int npi=0; if(ai.iconf>=0){ npi=confs[ai.iconf].npi; }
        int npj=0; if(aj.iconf>=0){ npj=confs[aj.iconf].npi; }
        const AtomType& ti    = params->atypes[ai.type];
        const AtomType& tj    = params->atypes[aj.type];
        const ElementType* ei = params->elementOfAtomType(ai.type);
        const ElementType* ej = params->elementOfAtomType(aj.type);
        double Ei = fabs(ei->Eaff);
        double Ej = fabs(ej->Eaff);
        double Qi = ei->Quff;
        double Qj = ej->Quff;

        double BO = 1 + _min( npi, npj );            
        double ri = ti.Ruff;
        double rj = tj.Ruff;

        double rBO = -0.1332*(ri+rj)*log( BO );
        double rEN = ri*rj*sq( sqrt(Ei) - sqrt(Ej) )/( Ei*ri + Ej*rj );
        double rij = ri + rj + rBO - rEN;

        //double kij   = 664.12 * Qi*Qj/( rij*rij*rij ); 
        double kij   = 28.79898 * Qi*Qj/( rij*rij*rij ); 
        
        printf( "bondUFF[%s,%s,%g] r=%g(%g,%g|%g,%g) k=%g(%g,%g) E(%g,%g)   %s %s %i %i \n", ti.name, tj.name, BO, rij,ri,rj,rBO,rEN,      kij,   Qi,Qj,Ei,Ej , ei->name,ej->name, ti.element, tj.element     );
        return { rij, kij };

    }

    void assignBondParams( int ib ){
        //printf( "MMFFBuilder::assignBondParams(ib=%i)\n", ib );
        Bond& b = bonds[ib];
        //printf( "MMFFBuilder::assignBondParams(ib=%i|ia=%i,ja=%i) order=%g \n", ib, b.atoms.i, b.atoms.j, b.order );
        const Atom& ai = atoms[b.atoms.i];
        const Atom& aj = atoms[b.atoms.j];
        printf( "MMFFBuilder::assignBondParams(ib=%i|ia=%i,ja=%i) types: %i=%s %i=%s \n", ib, b.atoms.i, b.atoms.j,  ai.type, params->atypes[ai.type].name, aj.type, params->atypes[aj.type].name );
        int order=1;
        if( (ai.iconf>=0)&&(aj.iconf>=0) ){ 
            const AtomConf& ci = confs[ai.iconf];
            const AtomConf& cj = confs[aj.iconf];
            order+=_min( ci.npi, cj.npi ); 
            //printf("assignBondParams[%i] (%i,%i|%i) pi(%i,%i) \n", ib,  ai.type, aj.type, order, ci.npi, cj.npi );

            // Assign pi-pi allignment 
            bool bpi=ci.npi>0;
            bool bpj=cj.npi>0; 
            if     ( bpi && bpj                       ){ b.kpp=Kpp_default;   }  // pi-pi alignement
            else if( bpi&&(cj.ne>0) || bpj&&(ci.ne>0) ){ b.kpp=Kpp_e_default; }  // pi-epair alignment

        }
        //printf( "MMFFBuilder::assignBondParams(ib=%i|ia=%i,ja=%i) types(%i,%i) order=%i \n", ib, b.atoms.i, b.atoms.j, ai.type, aj.type, order );
        //getBondTypeId( ai.type, aj.type, uint8_t order );
        params->getBondParams( ai.type, aj.type, order, b.l0, b.k );
        //printf("assignBondParams[%i] (%i,%i|%i) -> l0 %g k %g \n", ib,  ai.type, aj.type, order,   b.l0, b.k );
    }

    void assignAllBondParams(){ for(int i=0;i<bonds.size();i++){ assignBondParams(i); } };

    //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    void addCap(int ia,const Vec3d& hdir, Atom* atomj, double l=1.0 ){
        //printf( "addCap(%i)\n", ia );
        int ja=atoms.size();
        Atom atom_tmp;
        if(atomj==0){
            atom_tmp=capAtom;
            atomj=&atom_tmp;
        }
        atomj->pos = atoms[ia].pos + hdir*l;
        insertAtom(*atomj);
        // { // Debug
        //     Atom& a= atoms.back();
        //     printf( "addCap[%i : %i ] pos(%g,%g,%g) pos0(%g,%g,%g) hdir(%g,%g,%g) \n", atoms.size(), ia, a.pos.x,a.pos.y,a.pos.z, atoms[ia].pos.x,atoms[ia].pos.y,atoms[ia].pos.z, hdir.x,hdir.y,hdir.z );
        // }
        Bond B=capBond;
        B.atoms.set(ia,ja);
        int ib = insertBond( B );
        if(params)assignBondParams(ib);
    }

    //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    void addEpair(int ia, const Vec3d& hdir, double l=-0.5 ){
        //printf( "addEpairsByPi[%i, typ=%i] add epair[%i] type %i h(%g,%g,%g)\n", ia, ityp, i, ecap.type,   hs[ib].x,hs[ib].y,hs[ib].z );
        int ja=atoms.size();
        capAtom.type = itypEpair;
        if(params){ 
            int ityp = atoms[ia].type;
            capAtom.type = params->atypes[ityp].ePairType; 
            if(l<0)l=params->atypes[capAtom.type].Ruff;  // NOTE: we use Ruff as default length for epair, this is questionable, but this parameter has no other use for epair
        }else{ l=-l; }
        //printf( "addEpair[%i] type %i |h|=%g l=%g\n", ja, capAtom.type,   hdir.norm(), l );
        capAtom.pos = atoms[ia].pos + hdir*l;
        insertAtom(capAtom);
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
            //printf( "makeConfGeom nb=%i npi=%i \n", 3, npi );
            m.b.set_cross( hs[1]-hs[0], hs[2]-hs[0] );
            m.b.mul( -1/m.b.norm() );
            if(npi==0){ // sp3 no-pi
                if( 0 < m.b.dot( hs[0]+hs[1]+hs[2] ) ){ m.b.mul(-1.); }
                hs[3]=m.b;
            }else{
                hs[3]=m.b;
            }
        }else if(nb==2){ // defined by 2 sigma bonds
            //printf( "makeConfGeom nb=%i npi=%i \n", 2, npi );
            m.fromCrossSafe( hs[0], hs[1] );
            if      (npi==0){ // -CH2- like sp3 no-pi  => 109.5 deg.
                const double cb = 0.81649658092; // sqrt(2/3)
                const double cc = 0.57735026919; // sqrt(1/3)
                hs[nb  ] = m.c*cc+m.b*cb;
                hs[nb+1] = m.c*cc-m.b*cb;
                //printf( "hs:\n" );
            }else if(npi==1){ // =CH- like  sp 1-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
            }else{            // #C- like sp 2-pi
                hs[nb  ] = m.c;
                hs[nb+1] = m.b;
            }
        }else if(nb==1){
            //printf( "makeConfGeom nb=%i npi=%i \n", 1, npi );
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (npi==0){ // -CH3 like sp3 no-pi => 109.5 deg.
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
            }else if(npi==1){ // =CH2 like sp2 1-pi,   => 60 deg.
                //const double ca = 0.87758256189;  // 1/2  // this seems to be wrong,   0.87758256189^2 + 0.5^2 = 1.02015115293, should be 1.0
                const double ca = 0.86602540378;    // 1/2
                const double cc =-0.5;              // sqrt(1/8)
                hs[nb  ] = m.c*cc + m.a*ca;
                hs[nb+1] = m.c*cc - m.a*ca;
                hs[nb+2] = m.b;
            }else{            // #CH sp  2-pi
                hs[nb  ] = m.c*-1;
                hs[nb+1] = m.b;
                hs[nb+2] = m.a;
            }
        }else if(nb==0){
            //printf( "makeConfGeom nb=%i npi=%i \n", 0, npi );
            m.c = hs[0]; m.c.normalize();
            m.c.getSomeOrtho(m.b,m.a);
            if      (npi==0){ //  CH4 like sp3 no-pi  => 109.5 deg.
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


    void makeConfGeomCap(int nb, int npi, Vec3d* hs){
        Mat3d m;
        // --- sp3
        if(npi==0){
            if      (nb==3){
                //printf( "makeConfGeomCap() sp3 (1,1,1,0)\n" );
                m.b.set_cross( hs[1]-hs[0], hs[2]-hs[0] );
                m.b.mul( -1/m.b.norm() );
                if( 0 < m.b.dot( hs[0]+hs[1]+hs[2] ) ){ m.b.mul(-1.); }
                hs[3]=m.b;
            }else if(nb==2){
                //printf( "makeConfGeomCap() sp3 (1,1,0,0)\n" );
                m.fromCrossSafe( hs[0], hs[1] );
                const double cb = 0.81649658092; // sqrt(2/3)
                const double cc = 0.57735026919; // sqrt(1/3)
                hs[nb  ] = m.c*cc+m.b*cb;
                hs[nb+1] = m.c*cc-m.b*cb;
            }else if(nb==1){
                //printf( "makeConfGeomCap() sp3 (1,0,0,0)\n" );
                m.c = hs[0]; m.c.normalize();
                m.c.getSomeOrtho(m.b,m.a);
                const double ca = 0.81649658092;  // sqrt(2/3)
                const double cb = 0.47140452079;  // sqrt(2/9)
                const double cc =-0.33333333333;  // 1/3
                hs[nb  ] = m.c*cc + m.b*(cb*2) ;
                hs[nb+1] = m.c*cc - m.b* cb    + m.a*ca;
                hs[nb+2] = m.c*cc - m.b* cb    - m.a*ca;
            }
        }else
        // --- sp2
        if(npi==1){

        }else
        // --- sp1
        if(npi==2){
        }
        //printf("hs[0] l=%6.3f", hs[0].norm()); printVec(hs[0]);
        //printf("hs[1] l=%6.3f", hs[1].norm()); printVec(hs[1]);
        //printf("hs[2] l=%6.3f", hs[2].norm()); printVec(hs[2]);
        //printf("hs[3] l=%6.3f", hs[3].norm()); printVec(hs[3]);
    }




    bool makeConfGeomPi(int nb, int npi, const Vec3d& pi_dir, Vec3d* hs){
        //Mat3d m;
        if(nb==1){  // e.g. =O
            if(npi==1){
                Vec3d lf; lf.set_cross( pi_dir, hs[0] ); lf.normalize();
                hs[1] = lf*+0.86602540378+hs[0]*-0.5;
                hs[2] = lf*-0.86602540378+hs[0]*-0.5;
                // what was this old bad number: -0.87758256189
                return true;
            }else if(npi==2){
                Vec3d lf; lf.set_cross( pi_dir, hs[0] ); lf.normalize();
                hs[1] = hs[0]*-1.0;
                return true;
            }
        }
        return false;
    }


    Mat3d findMainAxes(int i0=0,int imax=-1, const bool bRemoveCog=true, const bool bRot=true, Vec3i permut=Vec3i{2,1,0} ){
        if(imax<0){ imax=atoms.size();}
        Mat3d M=Mat3dZero;
        Mat3d R;
        Vec3d cog=Vec3dZero;
        for(int i=i0; i<imax; i++){ cog.add(atoms[i].pos); }; cog.mul( 1./(imax-i0) );
        for(int i=i0; i<imax; i++){
            const Vec3d p = atoms[i].pos-cog; 
            M.addOuter(p,p,1.0);
        }
        Vec3d evs;
        M.eigenvals(evs);
        evs.sort();
        //printf( "findMainAxes() evs: %g %g %g \n", evs.x,evs.y,evs.z );
        M.eigenvec( evs.array[permut.x], R.a );
        M.eigenvec( evs.array[permut.y], R.b );
        M.eigenvec( evs.array[permut.z], R.c );
        if(bRot){
            for(int i=i0; i<imax; i++){
                Vec3d p;
                R.dot_to( atoms[i].pos-cog, p ); 
                if(bRemoveCog){ atoms[i].pos=p; }else{ atoms[i].pos=p+cog; }
            }
        }
        return R;
    }

    /*
    bool checkPointSymmetry( int i0=0,int imax=-1, double tol=0.1 ){
        if(imax<0){ imax=atoms.size();}
        for(int i=i0; i<imax; i++){ 
            atoms[i]
        }
    };
    */

    void findSymmetry( int* found, int i0=0,int imax=-1, double tol=0.1 ){
        /// this function should be called after findMainAxes()
        double tol2=tol*tol;
        if(imax<0){ imax=atoms.size();}
        bool bPoint=true;
        bool bX=true;
        bool bY=true;
        bool bZ=true;
        for(int i=i0; i<imax; i++){ 
            const Atom& A = atoms[i]; 
            bool bPoint_i=false;
            bool bX_i=false;
            bool bY_i=false;
            bool bZ_i=false;
            double r2min_x = 1e+300;
            double r2min_y = 1e+300;
            double r2min_z = 1e+300;
            double r2min_p = 1e+300;
            int imin_x=-1;
            int imin_y=-1;
            int imin_z=-1;
            int imin_p=-1;
            for(int j=i0; j<imax; j++){
                const Atom& B = atoms[j]; 
                if( A.type==B.type ){
                    //printf( "#### (a%i:t%i)-(a%i:t%i)\n", i,A.type,  j,B.type );
                    Vec3d d;
                    if( !bPoint_i ){   // point symmetry
                        d = A.pos - B.pos*-1.0;             //printf( "   point d(%g,%g,%g)\n", d.x,d.y,d.z );
                        double r2 = d.norm2();
                        if( r2<tol2 ){ bPoint_i=true; };
                        if( r2<r2min_p ){ r2min_p=r2; imin_p=j; }
                    }
                    if( !bX_i ){   // X-mirror symmetry
                        d=A.pos; d.x*=-1; d.sub(B.pos);     //printf( "   X d(%g,%g,%g) A(%g,%g,%g) B(%g,%g,%g)\n", d.x,d.y,d.z,  A.pos.x,A.pos.y,A.pos.z,  B.pos.x,B.pos.y,B.pos.z );   
                        double r2 = d.norm2();
                        if( r2<tol2 ){ bX_i=true; };
                        if( r2<r2min_x ){ r2min_x=r2; imin_x=j; }
                    }
                    if( !bY_i ){   // Y-mirror symmetry
                        d=A.pos; d.y*=-1; d.sub(B.pos);     //printf( "   Y d(%g,%g,%g)\n", d.x,d.y,d.z ); 
                        double r2 = d.norm2(); 
                        if( r2<tol2 ){ bY_i=true; };
                        if( r2<r2min_y ){ r2min_y=r2; imin_y=j; }
                    }
                    if( !bZ_i ){   // Z-mirror symmetry
                        d=A.pos; d.z*=-1; d.sub(B.pos);     //printf( "   Z d(%g,%g,%g)\n", d.x,d.y,d.z );  
                        double r2 = d.norm2();
                        if( r2<tol2 ){ bZ_i=true; };
                        if( r2<r2min_z ){ r2min_z=r2; imin_z=j; }
                    }
                }
            }
            //printf( "Atom[%i,t=%i] rmin(%g,%g,%g|%g) imin(%i,%i,%i|%i) \n", i, A.type, sqrt(r2min_x),sqrt(r2min_y),sqrt(r2min_z),sqrt(r2min_p), imin_x,imin_y,imin_z,imin_p );
            bPoint&=bPoint_i;   
            bX&=bX_i;
            bY&=bY_i;
            bZ&=bZ_i;
        }
        found[0]=bX;
        found[1]=bY;
        found[2]=bZ;
        found[3]=bPoint;
        printf( "findSymmetry() DONE mirror(x=%i,y=%i,z=%i) point(%i)\n", bX,bY,bZ, bPoint );
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

    int assignResonantType( int ityp_old )const{
        const AtomType& t = params->atypes[ityp_old];
        char tmp_name[8];
        sprintf( tmp_name, "%c_%s", t.name[0], "R" );
        int ityp = params->getAtomType(tmp_name, true);  // search for the subtypes in the atom type dictionary
        return ityp;
    }

    void assignAllSp3Types(){
        //printf("Builder::assignAllSp3Types() \n");
        for(int i=0; i<atoms.size(); i++){
            Atom& A = atoms[i];
            if( params->atypes[A.type].parrent != 0 ) continue; // already assigned
            if( A.iconf>=0 ){
                int npi = confs[A.iconf].npi;
                if( sp3types.count(A.type) == 0 ) continue;
                A.type = assignSp3Type_pi( A.type, npi );
                //printf( "Builder::assignAllSp3Types() [ia=%i] ityp=%i \n", i, A.type );
            }
        }
    }

    int getNeighType( int ia, int j, int* neighs ){
        int ja = neighs[ ia*4 + j];
        //printf( "getNeighType() ia %i j %i ja %i\n", ia, j, ja );
        if( (ja<0)||(ja>atoms.size()) ) return -1;
        return atoms[ ja ].type;
    }

    const char* getAtomTypeName(int ia){ return params->atypes[atoms[ia].type].name; };

    int countNeighPi( int ia, int* neighs ){
        int npi=0;
        int ic = atoms[ia].iconf;
        if(ic<0) return false;
        const AtomConf& c = confs[ic];
        for(int i=0; i<4; i++){
            int ja = neighs[i];
            if(ja<0) continue;
            int jc = atoms[ja].iconf;
            if(jc<0) continue;
            if( confs[jc].npi ){ npi++; }
            //printf( "atom[%i]ng[%i] c.npi=%i npi_sum=%i \n", ia, i, confs[jc].npi, npi );
        }
        return npi;
    }

    bool hasNeighborOfType( int ia, int n, const int* its, int* Ns, int* neighs ){
        int ic = atoms[ia].iconf;
        for(int j=0; j<n; j++){ Ns[j]=0; };
        if(ic<0) return false;
        const AtomConf& c = confs[ic];
        for(int i=0; i<4; i++){
            int ja = neighs[i];
            if(ja<0) continue;
            int jt = atoms[ja].type;
            for(int j=0; j<n; j++){
                if( jt==its[j] ) Ns[j]++;
            }
        }
        //printf("hasNeighborOfType(%i|%s){", ia, getAtomTypeName(ia) );
        //for(int j=0; j<n; j++){ printf( "%i,", bls[j]); };
        //printf("}\n" );
        return true;
    }



#define _Atyp(T)  const int i##T = params->getAtomType(#T);
#define _Btyp(E,T) const int E##_##T = _##E.size(); _##E.push_back(i##T);
#define _TCap(T1,T2) else if( ingt==i##T1 ){ itnew=i##T2; }
#define T(i)   (count[i]>0)

    int assignSpecialTypes( int* neighs ){
        //printf("#===========assignSpecialTypes() \n");
        // ------ C
        const int nt0=10;
        std::vector<int> _C;
        std::vector<int> _N;
        std::vector<int> _O;
        std::vector<int> _H;

        //  ==========  Type Definition

        // ------ C
        _Atyp(C_3) 
        _Atyp(C_2) 
        _Atyp(C_1) 
        _Atyp(C_R) 
        _Atyp(C_ene)  
        _Atyp(C_yne) 
        _Atyp(C_CH3)
        _Atyp(C_ald) 
        _Atyp(C_COO)
        // ------ O
        _Atyp(O)
        _Atyp(O_3)
        _Atyp(O_2)
        _Atyp(O_1)
        _Atyp(O_R)
        _Atyp(O_OH)
        _Atyp(O_ald)
        _Atyp(O_sCOO)
        _Atyp(O_pCOO)
        // ------ N
        _Atyp(N)
        _Atyp(N_3)
        _Atyp(N_2)
        _Atyp(N_1)
        _Atyp(N_R)
        // ------ H
        _Atyp(H)
        _Atyp(H_OH)
        _Atyp(H_COO)
        _Atyp(H_NH2)
        _Atyp(H_CH3)
        _Atyp(H_ene)
        _Atyp(H_yne)
        _Atyp(H_ald)

        //  ==========  Neighbor definition
        // ---- C
        _Btyp(C,C_2)
        _Btyp(C,C_R)
        _Btyp(C,O_3)
        _Btyp(C,O_2)
        _Btyp(C,N_2)
        // ---- O
        _Btyp(O,C_COO)
        _Btyp(O,C_2)
        _Btyp(O,C_R)
        _Btyp(O,C_ald)
        // ---- N
        _Btyp(N,C_2)
        _Btyp(N,C_R)
        _Btyp(N,O_2)

        int  na=atoms.size();
        int  count[100];
        int  nnew=0;
        //for(int i=0; i<na;i++){ printf("aneigs[%i]{%i,%i,%i,%i}\n", i, neighs[i*4+0],neighs[i*4+1],neighs[i*4+2],neighs[i*4+3]); }
        for(int ia=0; ia<na; ia++){
            int* ngs  = neighs+ia*4;
            Atom& A           = atoms[ia];
            const AtomType& t = params->atypes[A.type];
            int iZ            = t.iZ;
            int npineigh = countNeighPi( ia, ngs );
            int itnew =A.type;

            //printf( "atom[%i] %s=%i iZ %i \n", ia, t.name, A.type, iZ );
            switch (iZ){
                case 1: {  // H
                    int ingt =  getNeighType( 0, 0, ngs );       // if(ingt>-1)printf( "H[%i]-(%i|%s)\n", ia, ingt, params->atypes[ingt].name );
                    //printf("H ingt=%i\n", ingt );
                    if(ingt>=0){
                        if (false){}    
                        _TCap(O_OH ,H_OH )
                        _TCap(N_3,H_NH2)
                        _TCap(N_R,H_NH2)
                        _TCap(C_3,H_CH3)
                        _TCap(C_CH3,H_CH3)
                        _TCap(C_2,H_ene)
                        _TCap(C_R,H_ene)
                        _TCap(C_ene,H_ene)
                        _TCap(C_yne,H_yne)
                        _TCap(C_ald,H_ald)
                    }
                }break;
                case 6: { // C
                    //printf("C%i npineigh %i \n", ia, npineigh );
                    hasNeighborOfType( ia, _C.size(), &_C[0], count, ngs  );
                    if  ( (A.type==iC_2) ){
                        if      ( ( T(C_O_3) && T(C_O_2) ) ){ itnew = iC_COO; } // -COOH
                        else if ( ( T(C_O_2)             ) ){ itnew = iC_ald; } // -(C=O)-
                        //else if ( ( count[C_C_2] + count[C_C_R] + count[C_N_2] + count[C_O_2] )>=2 ){ itnew = iC_R;  } // Conjugated sp2 carbond ()
                        else if ( npineigh>=2 ){ itnew = iC_R; };
                    }
                }break;
                case 7: { // N
                    //printf("N \n");
                    //hasNeighborOfType( ia, _N.size(), &_N[0], count, ngs );
                    if( A.type==iN_3 ) if( npineigh>=1 ){ itnew = iN_R; };
                }break;
                case 8: { // O
                    //printf("O \n");
                    hasNeighborOfType( ia, _O.size(), &_O[0], count, ngs );
                    if( T(O_C_COO) ){ // COOH
                        if      (A.type==iO_3){ itnew=iO_sCOO; } // -OH of -COOH
                        else if (A.type==iO_2){ itnew=iO_pCOO; } // =O  of -COOH
                    }else if( T(O_C_ald)     ){ itnew=iO_ald;  } // -(C=O)-

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
#undef _Atyp
#undef _Btyp 
#undef _TCap
#undef T

    int assignSpecialTypesLoop( int nmax, int* neighs ){
        int nnew=0;
        int itr=0;
        for(itr=0; itr<nmax; itr++){
            //printf( "# --- assignSpecialTypesLoop[%i] \n", itr );
            int ni = assignSpecialTypes( neighs ); 
            nnew+=ni; 
            if( ni==0 ){ return nnew; } 
        }
        printf("ERROR: assignSpecialTypesLoop not converged in %i iterations\n", itr ); exit(0);
        return nnew;
    }

    // advanced atom-type assignement
    void assignTypes( int* neighs=0, int niterMax=10, bool bDeallocNeighs=true ){ 
        // assign all like X_3 or X_2 or X_1
        assignAllSp3Types();
        bDeallocNeighs &= (neighs==0);
        makeNeighs            ( neighs, 4        );
        assignSpecialTypesLoop( niterMax, neighs );
        //printAtomConfs(false);
        if(bDeallocNeighs)delete [] neighs;
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
        //printf( "addCaps[%i] ne %i b %i \n", ia, ne, bDummyEpair ); 
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
        int ne   = typ.nepair;
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

    void loadNeighbors(int ia, int nb, const int* neighs, Vec3d* hs ){
        if( nb>N_NEIGH_MAX   ){ printf("ERROR: loadNeighbors: nb %i > N_NEIGH_MAX %i \n", nb, N_NEIGH_MAX ); exit(0); }
        //if( ia>=atoms.size() ){ printf("ERROR: loadNeighbors: ia %i >= atoms.size() %i \n", ia, atoms.size() ); exit(0); }
        for(int i=0;i<nb;i++){
            int ib = neighs[i];
            if(ib<0 || ib>=bonds.size()){ printf("ERROR: loadNeighbors: ib %i >= nbonds.size() %i \n", ib, bonds.size() ); exit(0); }
            int ja = bonds[ib].getNeighborAtom(ia);
            if(ja<0 || ja>=atoms.size()){ printf("ERROR: loadNeighbors: ja %i >= atoms.size() %i \n", ja, atoms.size() ); exit(0); }
            hs[i]  = atoms[ja].pos - atoms[ia].pos;
            hs[i].normalize();
        }
    }

    AtomConf* tryGetNeighDirs( int ia, Vec3d* hs ){
        int ic=atoms[ia].iconf;
        if(ic<0){ return 0; }
        AtomConf& conf = confs[ic];
        int nb = conf.nbond;
        if(nb<2){ return 0; };
        loadNeighbors( ia, nb, conf.neighs, hs );
        makeConfGeom( nb, conf.npi, hs);
        //if(conf.npi>=0) conf.pi_dir = hs[3]; 
        return &conf;
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
        //printf( "makeSPConf[%i] npi=%i ne=%i ncap=%i bDummyEpair=%i bAddCaps=%i \n", ia, npi,ne,ncap, bDummyEpair,bAddCaps  );
        Vec3d hs[4];
        loadNeighbors(ia, nb, conf.neighs, hs );
        makeConfGeom(conf.nbond, npi, hs);
        if(bAddCaps && (ncap>0) )               addCaps( ia, ncap, ne, nb, hs );
        if(bDummyPi){ for(int i=0; i<npi; i++){ addCap(ia,hs[i+ncap+nb],&capAtomPi); } }
        conf.npi=npi;
        conf.ne =ne;
        //printf("-> "); println(conf);
    }

/*
    int addEpairsToAtoms(int ia, double l=0.5 ){
        int ic=atoms[ia].iconf;
        if(ic<0)return false;
        AtomConf& conf = confs[ic];
        int nb   = conf.nbond;
        Vec3d hs[4];
        for(int i=0;i<nb;i++){
            int ib = conf.neighs[i];
            int ja = bonds[ib].getNeighborAtom(ia);
            hs[i]  = atoms[ja].pos - atoms[ia].pos;
            hs[i].normalize();
        }
        makeConfGeom( nb, conf.npi, hs);
        //printf( "addEpairsToAtoms.1(%i) nb=%i npi=%i ne=%i ntot=%i \n", ia, conf.nbond, conf.npi, conf.ne, conf.n  );
        int ne=conf.ne;
        conf.ne=0;
        conf.n-=ne;
        //printf( "addEpairsToAtoms.2(%i) nb=%i npi=%i ne=%i ntot=%i \n", ia, conf.nbond, conf.npi, conf.ne, conf.n  );
        makeConfGeomPi(int nb, int npi, const Vec3d& pi_dir, Vec3d* hs);                
        for( int i=0; i<ne; i++ ){
            int ib=nb+i;
            //printf( "addEpairsToAtoms[%i] i=%i ib=%i h(%g,%g,%g) \n", ia, i, ib, hs[ib].x,hs[ib].y,hs[ib].z );
            addCap(ia,hs[ib],&capAtomEpair, l );
        }
        //printf( "addEpairsToAtoms.3(%i) nb=%i npi=%i ne=%i ntot=%i \n", ia, conf.nbond, conf.npi, conf.ne, conf.n  );
        return conf.ne;
    }

    bool setPiDir(int ia){
        Vec3d hs[4];
        AtomConf* conf = tryGetNeighDirs( int ia, Vec3d* hs );
        if( conf ){
            conf.npi=hs[3];
            return true;
        }
        return false;
    }
*/

    /**
     * Sets the pi direction vector for a given atom configuration according to the directions of pi-vectors of its neighbors.
     * 
     * @param ic The index of the atom configuration.
     * @return True if the pi direction vector was successfully set, false otherwise.
     */
    bool setPiByNeigh(int ic){
        //int ic=atoms[ia].iconf;
        //if(ic<0)return false;
        AtomConf& conf = confs[ic]; 
        Vec3d p = conf.pi_dir;
        double r2 = p.norm2();
        if(r2>0.1){ return false; } // already set
        p = Vec3dZero;
        const int* ngs = conf.neighs;
        for(int i=0; i<conf.nbond; i++){
            int ib = ngs[i];
            int ja = bonds[ib].getNeighborAtom(conf.iatom);
            //printf( "setPiByNeigh[%i] i=%i ib=%i ja=%i \n", ic, i, ib, ja );
            int jc = atoms[ja].iconf;
            if(jc<0) continue;
            Vec3d pi = confs[jc].pi_dir;
            if(i>0){
                double c = p.dot(pi); if(c<0){ pi.mul(-1.); };
            }
            p.add(pi);
        }
        double r = p.norm();
        if(r<1e-6){ return false; }
        p.mul(1./r);
        conf.pi_dir=p;
        return true;
    }

    // TBD what is ia0 here for? 
    int setPiLoop( int ia0=0, int imax=-1, int nMaxIter=10 ){
        //printf( "setPiLoop() confs.size()=%i \n", confs.size() );
        if(imax<0){ imax=atoms.size(); }
        for(int itr=0; itr<nMaxIter; itr++){
            int new_pi=0;
            for(int ic=0; ic<confs.size(); ic++){
                new_pi += setPiByNeigh(ic);
            }
            if(new_pi==0)return itr;
        }
        return nMaxIter;
    }

    bool autoConfEPi(int ia, double l=-0.5 ){
        int ic=atoms[ia].iconf;
        if(ic<0)return false;
        int ityp=atoms[ia].type;
        //params->assignRE( ityp, REQ );
        AtomConf& conf = confs[ic];
        int ne = params->atypes[ityp].nepair;
        if  ( (conf.nbond+conf.nH+ne)>N_NEIGH_MAX  ){ printf( "ERROR int autoConfEPi(ia=%i) ne(%i)+nb(%i)+nH(%i)>N_NEIGH_MAX(%i) => Exit() \n", ia, ne, conf.nbond, conf.nH, N_NEIGH_MAX ); exit(0);}
        else{ conf.ne=ne; }
        conf.fillIn_npi_sp3();

        if(params){ // limit the number of pi-bonds
            const ElementType* el = params->elementOfAtomType( atoms[ia].type );
            //printf( "atom[%i] typ=%i %s %s piMax=%i \n", ia, atoms[ia].type, params->atypes[atoms[ia].type].name, el->name,  el->piMax );
            conf.npi = _min( conf.npi, el->piMax );
        }

        int nb = conf.nbond;
        if(nb>=2){ // for less than 2 neighbors makeConfGeom does not make sense
            Vec3d hs[4];
            //AtomConf* conf = tryGetNeighDirs( ia, hs );
            loadNeighbors( ia, conf.nbond, conf.neighs, hs );
            makeConfGeom (     conf.nbond, conf.npi,    hs );

            // { // Debug
            //     sprintf( tmpstr, "atom%03i_hs.xyz", ia );
            //     FILE* fout=fopen(tmpstr, "w");
            //     fprintf(fout,"5\n");
            //     fprintf(fout,"#\n");
            //     Vec3d p = atoms[ia].pos;
            //     fprintf(fout, "%s %g %g %g\n",  params->atypes[ityp].name, p.x,p.y,p.z  );
            //     for(int i=0;  i<nb; i++){ Vec3d p_=p+hs[i]; fprintf(fout, "H %g %g %g\n", p_.x,p_.y,p_.z  ); }
            //     for(int i=nb; i<4;  i++){ Vec3d p_=p+hs[i]; fprintf(fout, "E %g %g %g\n", p_.x,p_.y,p_.z  ); }
            //     fclose(fout);
            // }

            if(conf.npi>0){ conf.pi_dir=hs[3]; }
            if( bDummyEpair && (conf.ne>0) ){
                conf.ne=0;
                conf.n-=ne;
                for( int i=0; i<ne; i++ ){
                    int ib=nb+i;
                    //printf( "addEpairsToAtoms[%i] i=%i ib=%i |h|=%g \n", ia, i, ib, hs[ib].norm() );
                    addEpair(ia,hs[ib],l);
                }
            }
        }
        return true;
    }
    int autoAllConfEPi( int ia0=0, int imax=-1 ){
        //printf( "autoAllConfEPi() bDummyEpair=%i \n", bDummyEpair  );
        int n=0;
        if(imax<0){ imax=atoms.size(); }
        for(int ia=ia0;ia<imax;ia++){
            if( autoConfEPi(ia) ){n++;}
        }
        return n;
    }

    bool addEpairsByPi(int ia, double l=-0.5 ){
        //printf( "addEpairsByPi[%i] \n", ia  );
        int ic=atoms[ia].iconf;
        if(ic<0)return false;
        //int ityp=atoms[ia].type;
        AtomConf& conf = confs[ic];
        int ne = conf.ne;
        if( (ne<1)||(conf.nbond>1)||(conf.npi<1) )return false;
        conf.ne=0;
        conf.n-=ne;
        int nb = conf.nbond;
        Vec3d hs[4];
        loadNeighbors ( ia, nb,       conf.neighs, hs );
        makeConfGeomPi( nb, conf.npi, conf.pi_dir, hs );
        //if(byPi){ makeConfGeomPi( nb, conf.npi, conf.pi_dir, hs ); } // NOTE: we need to asign pi_dir before calling makeConfGeomPi(), this is however necessary for atoms like =O which do not have other bonds direction of e-pair is not defined if pi-plane is not defined
        //else    { makeConfGeom  ( nb, conf.npi,              hs ); }  
        //printf( "addEpairsByPi[%i, typ=%i=%s] npi=%i hs[0](%6.3f,%6.3f,%6.3f) hs[1](%6.3f,%6.3f,%6.3f) hs[2](%6.3f,%6.3f,%6.3f) hs[3](%6.3f,%6.3f,%6.3f) \n", ia, ityp, params->atypes[ityp].name,  hs[0].x,hs[0].y,hs[0].z,   hs[1].x,hs[1].y,hs[1].z,   hs[2].x,hs[2].y,hs[2].z, hs[3].x,hs[3].y,hs[3].z );

        // { // Debug
        //     sprintf( tmpstr, "atom%03i_hs.xyz", ia );
        //     FILE* fout=fopen(tmpstr, "w");
        //     fprintf(fout,"5\n");
        //     fprintf(fout,"#\n");
        //     Vec3d p = atoms[ia].pos;
        //     fprintf(fout, "%s %g %g %g\n",  params->atypes[ityp].name, p.x,p.y,p.z  );
        //     for(int i=0;  i<nb; i++){ Vec3d p_=p+hs[i]; fprintf(fout, "H %g %g %g\n", p_.x,p_.y,p_.z  ); }
        //     for(int i=nb; i<4;  i++){ Vec3d p_=p+hs[i]; fprintf(fout, "E %g %g %g\n", p_.x,p_.y,p_.z  ); }
        //     fclose(fout);
        // }

        for( int i=0; i<ne; i++ ){
            int ib=nb+i;
            //printf( "addEpairsToAtoms[%i] i=%i ib=%i |h|=%g |hb|=%g |hpi|=%g  l=%g \n", ia, i, ib, hs[ib].norm(), hs[0].norm(), conf.pi_dir.norm(), l );
            addEpair(ia,hs[ib],l);
        }
        return true;
    }
    int addAllEpairsByPi( int ia0=0, int imax=-1, bool byPi=true ){
        int n=0;
        if(imax<0){ imax=atoms.size(); }
        for(int ia=ia0;ia<imax;ia++){
            if( addEpairsByPi(ia ) ){n++;}
        }
        return n;
    }

    bool addCapByPi(int ia, int cap_typ, double l=1.0 ){
        int ic=atoms[ia].iconf;
        if(ic<0)return false;
        //int ityp=atoms[ia].type;
        AtomConf& c = confs[ic];
        c.updateNeighs();
        int nCap = 4 - (c.ne + c.npi + c.nbond);
        if(nCap<=0) return false;
        Atom A = capAtom; A.type = cap_typ; A.iconf=-1;
        printf( "addCapsByPi[%i] nCap=%i conf: ", ia, nCap ); c.print(); //printf("\n");
        Vec3d hs[4];
        loadNeighbors ( ia, c.nbond,    c.neighs, hs );
        //makeConfGeomPi( c.nbond, c.npi, c.pi_dir, hs );
        //makeConfGeom( c.nbond, c.npi, hs);
        makeConfGeomCap( c.nbond, c.npi, hs );
        int nb0 = c.nbond;
        for( int i=0; i<nCap; i++ ){
            int ib=nb0+i;
            printf( "addCapsByPi[%i] i=%i ib=%i |h|=%g |hb|=%g |hpi|=%g  l=%g \n", ia, i, ib, hs[ib].norm(), hs[0].norm(), c.pi_dir.norm(), l );
            //addEpair(ia,hs[ib],l);
            addCap(ia,hs[ib], &A, l );
        }
        return true;
    }
    int addAllCapsByPi( int cap_typ, int ia0=0, int imax=-1, bool byPi=true ){
        int n=0;
        if(imax<0){ imax=atoms.size(); }
        for(int ia=ia0;ia<imax;ia++){
            if( addCapByPi(ia,cap_typ) ){n++;}
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
        //printf("Builder::countPiE() n %i i0 %i \n", n,i0);
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

    int selectBondsBetweenTypes( int imin, int imax, int it1, int it2, bool byZ=false, bool bOnlyFirstNeigh=false ){
        printf( "Builder::selectBondsBetweenTypes[%i,%i] it1 %i it2 %i byZ %i bOnlyFirstNeigh %i \n", imin, imax, it1, it2, byZ, bOnlyFirstNeigh );
        selection.clear();
        std::unordered_set<int> nodes;
        for(int ib=imin; ib<imax; ib++){
            const Vec2i& b = bonds[ib].atoms;
            const Atom&  A = atoms[b.a];
            const Atom&  B = atoms[b.b];
            bool match;
            bool reverse;
            if( byZ ){
                int iz1 = params->atypes[A.type].iZ;
                int iz2 = params->atypes[B.type].iZ;
                //printf( "selectBondsBetweenTypes[%i] zs(%i,%i) ts(%i,%i) \n", ib, iz1,iz2,   it1, it2  );
                match   = (iz1==it1)&&(iz2==it2);
                reverse = (iz1==it2)&&(iz2==it1); 
            }else{
                match   = (A.type==it1)&&(B.type==it2);
                reverse = (A.type==it2)&&(B.type==it1); 
            }
            match |= reverse;
            if(match){
                //printf( "selectBondsBetweenTypes match ib=%i \n", ib, b.a,b.b  );
                if(bOnlyFirstNeigh){
                    int ia;
                    if(reverse){ ia=b.b; }else{ ia=b.a; };
                    int na = nodes.count(ia); 
                    if( na==0 ){
                        selection.insert(ib);
                        nodes.insert(ia);
                    }
                }else{
                    selection.insert(ib);
                }
            }
        }
        return selection.size();
    }


    int pickBond( const Vec3d& ro, const Vec3d& rh, double Rmax ){
        return rayPickBond( ro, rh, bonds.size(), [&](int ib,Vec3d&pa,Vec3d&pb){ 
            Vec2i b = bonds[ib].atoms; pa=atoms[b.i].pos; pb=atoms[b.j].pos;
        }, Rmax, false );
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

    bool addTorsionsToBond( int ibond ){
        Vec2i b = bonds[ibond].atoms;
        const AtomConf* ca = getAtomConf(b.i);
        const AtomConf* cb = getAtomConf(b.j);
        if( (ca==0) || (cb==0)) return false;
        int ityp = atoms[b.i].type;
        int jtyp = atoms[b.j].type;
        for(int i=0; i<ca->nbond; i++ ){
            for(int j=0; j<cb->nbond; j++ ){
                int ib = ca->neighs[i];
                int jb = cb->neighs[j];
                int ia = bonds[ ib ].getNeighborAtom(b.i);  int iat = atoms[ia].type;
                int ja = bonds[ jb ].getNeighborAtom(b.j);  int jat = atoms[ja].type;
                Dihedral tor;
                DihedralType* dtyp = params->getDihedralType( iat, ityp, jtyp, jat, bonds[ibond].type );
                tor.a0    = dtyp->ang0;
                tor.k     = dtyp->ang0;
                tor.n     = dtyp->n;
                tor.atoms = Quat4i{iat,ityp,jtyp,jat};
                tor.bonds = {ibond,ib,jb};
                dihedrals.push_back( tor );
            }
        }
        return true;
    }

    void autoTorsions(){
        for(int ib=0; ib<bonds.size(); ib++){
            addTorsionsToBond( ib );
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


void assignTorsions( bool bNonPi=false, bool bNO=true ){
    printf( "MM::Builder::assignTorsions()\n" );
    int default_order = 2;
    int nb = bonds.size();
    //Quat4i*  tors2atom =0;
    //Quat4d*  torsParams=0;
    //std::vector<Quat4i> tors2atom;
    //std::vector<Quat4d> torsParams;
    for(int ib=0; ib<nb; ib++ ){
        //printf( "assignTorsions[ib=%i]\n", ib );
        Vec2i b = bonds[ib].atoms;
        int ic = atoms[b.i].iconf;
        int jc = atoms[b.j].iconf;
        if( (ic<0)||(jc<0) ) continue;  // must be bond between node atoms
        const AtomConf& ci = confs[ic];
        const AtomConf& cj = confs[jc];
        printf( "#### assignTorsions[ib=%i,%s-%s]   ci.npi=%i cj.npi=%i bNonPi=%i \n", ib, params->atypes[ atoms[b.i].type ].name, params->atypes[ atoms[b.j].type ].name, ci.npi, cj.npi, bNonPi );
        if((ci.npi>1)||(cj.npi>1)) continue;          // No triple-bond torsion
        if((ci.npi==0)||(cj.npi==0)){ // sp3 atoms only if O or N attaced to sp2
            if(!bNonPi)continue;
            int iZ = params->atypes[ atoms[b.i].type ].iZ;
            int jZ = params->atypes[ atoms[b.j].type ].iZ;
            printf( "#### assignTorsions[ib=%i,%s-%s]  iZ=%i jZ=%i ci.npi=%i cj.npi=%i bNonPi=%i \n", ib, params->atypes[ atoms[b.i].type ].name, params->atypes[ atoms[b.j].type ].name,  iZ,jZ, ci.npi, cj.npi, bNonPi );
            if(( ((iZ==7)||(iZ==8)) && ( cj.npi==1 ) )
             ||( ((jZ==7)||(jZ==8)) && ( ci.npi==1 ) )){
            }else{ continue; }
        }

        int order = bonds[ib].type;
        if(order<0){
            printf( "WARRNING: in assignTorsions[ib=%i] order=%i changed to %i \n", ib, order, default_order ); 
            order=default_order;
        } 

        for(int i=0; i<ci.nbond; i++ ){
            if( ci.neighs[i] == ib ) continue;
            int iia = bonds[ci.neighs[i]].getNeighborAtom(b.i);
            for(int j=0; j<cj.nbond; j++ ){
                if( cj.neighs[j] == ib ) continue;
                int jja = bonds[cj.neighs[j]].getNeighborAtom(b.j);
                printf( "dih atoms{%i,%i,%i,%i} types{%i,%i,%i,%i} type_names{%s,%s,%s,%s}\n",   iia,b.i,b.j,jja,    atoms[iia].type,atoms[b.i].type,atoms[b.j].type,atoms[jja].type,  params->atypes[atoms[iia].type].name,params->atypes[atoms[b.i].type].name,params->atypes[atoms[b.j].type].name,params->atypes[atoms[jja].type].name     );
                DihedralType* diht = params->getDihedralType( atoms[iia].type, atoms[b.i].type, atoms[b.j].type, atoms[jja].type, order, true, true );
                if(diht==0){ printf("ERROR in MM::Builder::assignTorsions(ib=%i,%i,%i) cannot find angle type(%i,%i,%i,%i|%i) =>Exit()\n", ib,i,j, atoms[iia].type, atoms[b.i].type, atoms[b.j].type, atoms[jja].type, order ); exit(0); };
                //tors2atom .push_back( Quat4i{ iia,b.i,b.j,jja} );
                //torsParams.push_back( Quat4d{cos(dih->ang0), sin(dih->ang0), dih->k, dih->n} );
                dihedrals.push_back( Dihedral{ diht->bo, Quat4i{iia,b.i,b.j,jja}, diht->n, diht->k, diht->ang0 } );

            }
        }
    }
    printf( "MM::Builder::assignTorsions() DONE\n" );
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
            double Rj = params->atypes[A.type].Ruff;
            //double Rj = (R0 + A.REQ.x)*Rfac;    // Using RvdW
            double R = (Rj+R0)*Rfac;
            //printf( "touchingAtoms[%i] R=%g Ri+j=%g Ri=%g Rj=%g \n", i, R, R0+Rj, R0, Rj );
            if(  dp.norm2() < (R*R) ){
                //if(verbosity>2) 
                //printf( "bond[%i,%i] r %g R %g(%g,%g) \n", i0-1, i,    dp.norm(), R, R0, Rj );
                found.push_back(i);
            }
            //else{printf( "NON bond[%i,%i] r %g R %g \n", i0-1, i, dp.norm(), R );}
        }
    }

    // ToDo:  R=0.5 worsk for hydrocarbons (HCNOF), R=0.65 works for silicon (SiH), we should assign it according to the bond types maybe ?
    int autoBonds( double R=-1.20, int i0=0, int imax=-1 ){
        //printf( "MM::Builder::autoBonds() \n" );
        if(verbosity>0){ printf( "MM::Builder::autoBonds() \n" ); }
        if(imax<0)imax=atoms.size();
        bool byParams = (R<0);
        double Rfac   = -R;
        //if( byParams && (params==0) ){ printf("ERROR in MM::Builder.autoBonds() byParams(R<0) but params==NULL \n"); exit(0); }
        std::vector<int> found;
        int nbond=0;
        for(int i=i0; i<imax; i++){
            //printf( "autoBonds() atom[%i] typ(%i==%s)\n", i, atoms[i].type,  params->atypes[ atoms[i].type ].name );
            const Atom& A = atoms[i];
            bool bCap_i = capping_types.count( A.type ) > 0; 
            //double Ri = A.REQ.x;   // Using RvdW
            double Ri = params->atypes[A.type].Ruff;   
            found.clear();
            touchingAtoms( i+1, imax, A.pos, Ri, Rfac, found );
            for(int j:found){
                if( bCap_i ){ if( capping_types.count( atoms[j].type ) > 0 ) continue ; }  // prevent bonds between two capping atoms
                //printf( "bond[%i] (%i,%i) %s-%s \n", nbond, i,j,  params->atypes[ atoms[i].type ].name, params->atypes[ atoms[j].type ].name    );
                bondBrush.ipbc=Vec3i8{ -1,-1,-1 };
                bondBrush.atoms={i,j};
                insertBond( bondBrush );
                nbond++;
            }
        }
        //printf( "MM::Builder::autoBonds() DONE !!!!!!!!!!!!!!!!\n\n\n" );
        return nbond;
    }

    inline Vec3d pbcShift( Vec3i G ){ return lvec.a*G.a + lvec.b*G.b + lvec.c*G.c; }

    // find 
    int autoBondsPBC( double R=-1.35, int i0=0, int imax=-1, Vec3i npbc=Vec3iOne ){
        //printf( "MM::Builder::autoBondsPBC() \n" );
        //if(verbosity>0){ printf( "MM::Builder::autoBondsPBC() \n" );                             }
        if(verbosity>1){ printf( "MM::Builder::autoBondsPBC() builder.lvec: \n" ); lvec.print(); };
        if(imax<0)imax=atoms.size();
        bool byParams = (R<0);
        double Rfac   = -R;
        //if( byParams && (params==0) ){ printf("ERROR in MM::Builder.autoBonds() byParams(R<0) but params==NULL \n"); exit(0); }
        std::vector<int> found;
        int nbond=0;
        for(int i=i0; i<imax; i++){
            const Atom& A = atoms[i];
            bool bCap_i = capping_types.count( A.type ) > 0; 
            //double Ri = A.REQ.x;  // Using RvdW
            double Ri = params->atypes[A.type].Ruff;
            int ipbc=0;
            //if(verbosity>1)
            //printf( "autoBondsPBC() Atom[%i] R %g \n", i, R );
            for(int ix=-npbc.x;ix<=npbc.x;ix++){
                for(int iy=-npbc.y;iy<=npbc.y;iy++){
                    for(int iz=-npbc.z;iz<=npbc.z;iz++){
                        int   j0=i+1;
                        //Vec3d vpbc = lvec.a*ix + lvec.b*iy + lvec.c*iz;
                        //Vec3d vpbc; lvec.dot_to_T( {(double)ix,(double)iy,(double)iz} );
                        //Vec3d p = A.pos - pbcShift( {ix,iy,iz} );
                        Vec3d p = A.pos - lvec.lincomb( ix, iy, iz );
                        found.clear();
                        // find overlapping atoms
                        touchingAtoms( j0, imax, p, Ri, Rfac, found ); 
                        //if(i==12)printf( "# pbc[%i,%i,%i][%i] nfound %i \n", ix,iy,iz, ipbc, found.size() );
                        for(int j:found){
                            if( bCap_i ){ if( capping_types.count( atoms[j].type ) > 0 ) continue ; }  // prevent bonds between two capping atoms
                            //bondPBC.push_back( {ix,iy,iz} );
                            //bondBrush.ipbc=Vec3i8{ix,iy,iz};
                            bondBrush.ipbc=Vec3i8{ (int8_t)ix, (int8_t)iy, (int8_t)iz };
                            bondBrush.atoms={i,j};
                            insertBond( bondBrush );
                            nbond++;
                        }
                        ipbc++;
                    }
                }
            }
        }
        //printf( "MM::Builder::autoBondsPBC() DONE\n" );
        return nbond;
    }


    Vec3i shortesPBCbond( int ib, Vec3i npbc=Vec3iOne ){
        Bond& b = bonds[ib];
        Vec3d dij = atoms[b.atoms.j].pos - atoms[b.atoms.i].pos;
        Vec3i ipbc_min = {0, 0, 0};
        double r2min = 1e300;
        for (int ix = -npbc.x; ix <= npbc.x; ix++) {
            for (int iy = -npbc.y; iy <= npbc.y; iy++) {
                for (int iz = -npbc.z; iz <= npbc.z; iz++) {
                    Vec3d p = dij + lvec.lincomb(ix, iy, iz);
                    double r2 = p.norm2();
                    if (r2 < r2min) {
                        r2min = r2;
                        ipbc_min = {ix, iy, iz};
                    }
                }
            }
        }
        b.ipbc = Vec3i8{(int8_t)ipbc_min.x, (int8_t)ipbc_min.y, (int8_t)ipbc_min.z};
        return ipbc_min;
    }


    bool checkAllAtomsBonded( bool bPrint=true, bool bExit=true, int nbmax=N_NEIGH_MAX, int nbmin=1 ){
        std::vector<int> nbonds(atoms.size(),0);
        for(int ib=0; ib<bonds.size(); ib++){
            const Bond& b = bonds[ib];
            nbonds[b.atoms.i]++;
            nbonds[b.atoms.j]++;
        }
        bool bRet = false;
        for(int i=0; i<atoms.size(); i++){
            int nb=nbonds[i];
            if( (nb<nbmin)||(nb>nbmax) ){
                printf( "ERROR in Builder::checkAllAtomsBonded(): atom[%i].nbonds is out of range (%i,%i)\n", i, nb, nbmin, nbmax );
                bRet=true;
            }
        }
        return bRet;
    }

    bool checkNumberOfBonds( bool bPrint=true, bool bExitOnError=true, bool bAllAtomsBonded=true ){
        const int na = atoms.size();
        std::vector<int> nbonds(na,0); 
        //for(int ia=0; ia<na; ia++){ printf( "checkNumberOfBonds nbonds[%i]=%i\n", ia, nbonds[ia] ); };
        for(int ib=0; ib<bonds.size(); ib++){
            const Vec2i& b = bonds[ib].atoms;
            //printf( "checkNumberOfBonds[ib=%i] b(%i,%i)\n", ib, b.i, b.j );
            nbonds[ b.i ] ++;
            nbonds[ b.j ] ++; 
        }
        bool err = false;
        for(int ia=0; ia<na; ia++){
            const Atom& A = atoms[ia];
            const AtomType& t = params->atypes[A.type];
            //printf( "checkNumberOfBonds nbonds[%i]=%i\n", ia, nbonds[ia] );
            if(bAllAtomsBonded) { if( nbonds[ia]<=0 ){ err|=true; if( bPrint ){ printf( "WARNING checkNumberOfBonds[%i] `%s` has no bonds nbond(%i)<1 bonds\n", ia, t.name, nbonds[ia] ); } } }
            if(A.iconf>=0){
                const AtomConf& c = confs[A.iconf];
                int nb = nbonds[ia];
                if( nb>N_NEIGH_MAX){ err|=true; if( bPrint ){ printf( "WARNING checkNumberOfBonds[%i] `%s` nbonds(%i)>N_NEIGH_MAX (%i)\n", ia, t.name, nb, N_NEIGH_MAX ); } }
                if( nb>t.valence  ){ err|=true; if( bPrint ){ printf( "WARNING checkNumberOfBonds[%i] `%s` nbonds(%i)>valence     (%i)\n", ia, t.name, nb, t.valence   ); } }
                if( nb!=c.nbond   ){ err|=true; if( bPrint ){ printf( "WARNING checkNumberOfBonds[%i] `%s` nbonds(%i)!=conf.nbond (%i)\n", ia, t.name, nb, c.nbond     ); } } 
            }else{
                if( nbonds[ia]>1  ){ err|=true; if( bPrint ){ printf( "WARNING checkNumberOfBonds[%i] `%s` capping atom has %i>1 bonds\n", ia, t.name, nbonds[ia] ); } }
            }
        }
        if( err && bExitOnError ){ printf("ERROR in checkNumberOfBonds() => Exit()\n"); exit(0); }
        return err;
    }

    // =============== Configurations

    bool checkConfValid( int ia, bool bPrint=true )const{
        int natoms = atoms.size();
        int nbonds = bonds.size();
        int ic     = atoms[ia].iconf;
        if( ic<0 ) return true;
        if( ic>=confs.size() ){  if( bPrint ){ printf( "checkConfValid(ia=%i) ic(%i) >= confs.size(%li)\n", ia, ic, confs.size() ); } return false; }
        const AtomConf& conf = confs[ic];
        int nbond = conf.nbond;
        if( nbond>N_NEIGH_MAX ){  if( bPrint ){ printf( "checkConfValid(ia=%i|ic=%i) nbond(%i) > N_NEIGH_MAX(%i)\n", ia, ic, nbond, N_NEIGH_MAX ); } return false; }
        for( int i=0; i<nbond; i++){
            int ib = conf.neighs[i];
            if( (ib<0) || (ib>=nbonds) ){  if( bPrint ){ printf( "checkConfValid(ia=%i|ic=%i) bond[%i] ib(%i) >= bonds.size(%i)\n", ia, ic, i, ib, nbonds ); } return false; }
            Vec2i b = bonds[ib].atoms;
            if( (b.i>=natoms) || (b.j>=natoms) ){  if( bPrint ){ printf( "checkConfValid(ia=%i|ic=%i) bond[%i] atoms(%i,%i) >= atoms.size(%i)\n", ia, ic, i, ib, b.i, b.j, natoms ); } return false; }
            if( (b.i!=ia) && (b.j!=ia) ){  if( bPrint ){ printf( "checkConfValid(ia=%i|ic=%i) bond[%i] ib(%i)  none of atoms b(i=%i,j=%i) != ia(%i)\n", ia, ic, i, ib, b.i, b.j, ia ); } return false; }
        }
        return true;
    }

    int checkConfsValid( bool bExit=true, bool bPrint=true, bool bPrintArraysExit=true )const{
        int nbad = 0;
        for(int ia=0;ia<atoms.size(); ia++){ if( !checkConfValid( ia, bPrint ) ) nbad++; }
        if( nbad>0) {
            if( bPrint ){ printf( "checkConfsValid() failed: nbad(%i)\n", nbad ); }
            if( bExit  ){ 
                if( bPrintArraysExit ){
                    printAtomConfs();
                    printBonds();
                }
                exit(0); 
            }
        }
        return nbad;
    }

    int checkBond2Conf(bool bPrint=true)const{
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
            //printf( "Builder::checkAtomHasBond ia=%i neighs[%i]=%i \n", ia, i, ib );
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
        //if(! (ba && bb))printf( "Builder::checkBondInNeighs !bond[%i|%i,%i] %i %i\n", ib, b.a, b.b, ba, bb);
        return ba && bb;
    }

    bool checkBondsInNeighs( bool bPrint=true ){
        bool ret = false;
        for(int ib=0; ib<bonds.size(); ib++ ){
            //printf("\----bond[%i]\n", ib);
            if( ! checkBondInNeighs( ib ) ){
                if(bPrint){ Vec2i b=bonds[ib].atoms; printf("WARNING bond[%i|%i,%i] is not in neighborhood \n", ib, b.a, b.b); printAtomConf(b.a); puts(""); printAtomConf(b.b); puts(""); }
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

    void permutAtoms(int* permut, bool doBonds=false, bool bGroups=true ){
        bGroups &= (atom2group.size()>0);
        std::vector<int> atom2group_;
        if(bGroups){
            if(atom2group.size()>0){
            if (atom2group.size()!=atoms.size()){
                printf("ERROR permutAtoms bGroups && atom2group.size(%i)!=atoms.size(%i)\n", atom2group.size(), atoms.size());
                exit(1);
            }
            atom2group_ = atom2group;
            }
        }

        for(Bond&     b: bonds){ b.atoms.a=permut[b.atoms.a];  b.atoms.b=permut[b.atoms.b]; }
        // Confs are not re-shuffled because they point to bonds, not atoms
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
            if(bGroups){
                atom2group[ja] = atom2group_[ia];
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

    int findBondsToAtom(int ia, int* is=0, bool bPrint=false, int nmax=1000000)const{
        int i=0;
        for(int ib=0; ib<bonds.size(); ib++){
            if( bonds[ib].atoms.anyEqual(ia) ){ 
                if(is){ is[i]=1; i++; }
                if(bPrint){ printf("bond[%i]",ib); bonds[ib].print(); puts(""); }
                if(i>=nmax){ printf("WARRNING findBondsToAtom[ia=%i] found %i>nmax(%i) bonds \n", ia, i, nmax); break; }
            };
        };
        return i;
    }

    void numberAtoms(){
        for(int i=0; i<atoms.size(); i++){ atoms[i].id=i; };
    }

    void setup_atom_permut( bool bPrint=false ){
        atom_permut.resize(atoms.size());
        for(int i=0; i<atoms.size(); i++){ 
            atom_permut[ atoms[i].id ]=i; 
            //printf( "setup_atom_permut[%i]-> %i \n", i, atoms[i].id );
        };
        //for(int i=0; i<atoms.size(); i++){ printf( "atom_permut[%i] %i \n", i, atom_permut[i] ); };
    }

    void printSizes()const{ printf( "sizes: atoms(%i|%i) bonds(%i) angles(%i) dihedrals(%i) \n", atoms.size(), confs.size(), bonds.size(), angles.size(), dihedrals.size() ); };

    void printAtoms()const{
        printf(" # MM::Builder.printAtoms(na=%i) \n", atoms.size() );
        for(int i=0; i<atoms.size(); i++){
            printf("atom[%i]",i); atoms[i].print(); puts("");
        }
    }
    void printBonds()const{
        // printf(" # MM::Builder.printBonds(nb=%i) \n", bonds.size() );
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
        if(params){ printf("atom[%i] t%i=`%s` c%i ", ia, A.type, params->atypes[A.type].name, A.iconf );  }
        else      { printf("atom[%i] t%i c%i ", ia, A.type, A.iconf ); }        
        if(A.iconf>=0){
            const AtomConf& c = confs[A.iconf];
            printf("nbpeH(%i|%i,%i,%i,%i) neighs{", c.n,c.nbond,c.npi,c.ne,c.nH);
            for(int i=0; i<N_NEIGH_MAX; i++){
                int ib = c.neighs[i];
                int ja=-2;
                if(ib>=0){
                    ja = bonds[ib].getNeighborAtom(ia);
                }
                printf("%3i,", ja );
            }
            printf("}" );
        }
    }

    inline void printBondsOfAtom( int ia ){
        int ic = atoms[ia].iconf;
        if(ic<0){printf("printBondsOfAtom(%i): atom has no conf (capping) \n", ia); return; };
        const AtomConf& c = confs[ic];
        // printf("printBondsOfAtom(%i): ", ia);
        for(int i=0; i<c.nbond; i++){
            int ib = c.neighs[i];
            Vec2i b = bonds[ib].atoms;
            // printf( "(%i|%i,%i) ", ib, b.i, b.j );
        }
        printf("\n");
    }

    void printAtomConfs( bool bOmmitCap=true, bool bNeighs=false )const{
        // printf(" # MM::Builder.printAtomConfs(na=%i,nc=%i) \n", atoms.size(), confs.size() );
        for(int i=0; i<atoms.size(); i++){ if( bOmmitCap && (atoms[i].iconf==-1) )continue;  if(bNeighs){printAtomNeighs(i);}else{printAtomConf(i);} puts(""); }
    }

    void printAtomTypes( )const{
        // printf(" # MM::Builder.printAtomTypes(na=%i,nc=%i) \n", atoms.size(), confs.size() );
        for(int i=0; i<atoms.size(); i++){ 
            const Atom& A = atoms[i];
            int it=A.type;
            if(A.iconf>=0){
                AtomConf c = confs[A.iconf];
                // printf( "atom[%4i] %ip %ie {%3i,%3i,%3i,%3i} type[%3i]=`%s`   \n", i, c.npi, c.ne, c.neighs[0],c.neighs[1],c.neighs[2],c.neighs[3], it, params->atypes[it].name ); 
            } else{
                // printf( "atom[%4i] type[%3i]=`%s`   \n", i, it, params->atypes[it].name  ); 
            }
            
        }
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
                double dQ = ( params->elementOfAtomType(A.type)->Eaff - params->elementOfAtomType(B.type)->Eaff  ) * factor;
                if(itr>0){  dQ += -A.REQ.z + B.REQ.z;  dQ*=damp; }
                A.REQ.z += dQ;
                B.REQ.z -= dQ;
            }
        }
        for(int ia=0; ia<atoms.size(); ia++){ printf("atom[%i] Q=%g \n", ia, atoms[ia].REQ.z ); }; //exit(0);
    }

    void str2groups(const char* buff, int maxSteps=1024, int _0 = 0){
        // this function reads groups from a string and stores them in atom2group 
        // groups are separated by semicolon ";" and atoms in one group are separated by comma ","
        // there can be any number of groups and each group can have any number of atoms 
        // example: "1,3,5; 2,4,7; 8,9,12"
        
        const char* p = buff;
        int group = 0;
        int atom;
        atom2group.resize(atoms.size(),-1); // -1 means not assigned
        
        printf(" # MM::Builder.str2groups(%s) \n", buff );
        
        char* endptr;
        while (*p && group < maxSteps) { // safety limit on groups
            // Skip any whitespace or non-numeric characters until we find a digit or end
            while (*p && !(*p >= '0' && *p <= '9')) {
                if (*p == ';') {
                    group++;
                    printf("str2groups() moving to group: %i\n", group);
                }
                p++;
            }
            // If we hit the end, break
            if (!*p) break;
            // Try to convert the number
            atom = strtol(p, &endptr, 10);
            atom -= _0;
            if (p != endptr) { // if conversion successful
                printf("str2groups() group: %i atom: %i\n", group, atom);
                if (atom >= 0 && atom < atoms.size()) {
                    atom2group[atom] = group;
                }
                p = endptr; // advance past the number
            } else {
                p++; // shouldn't happen, but just in case
            }
        }
        
        if (group >= maxSteps) {
            printf("Warning: str2groups reached maximum groups (%d). Parsing terminated.\n", maxSteps);
        }
        printf("Finished parsing. Found %d groups\n", group + 1);
    }
        
    void printAtom2Groups()const{
        printf("printAtom2Groups() n=%i na=%i \n", atom2group.size(), atoms.size() );
        for(int i=0; i<atom2group.size(); i++){
            printf("atom[%i] group[%i] \n", i, atom2group[i] );
        }
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
        fprintf(pfile,"M  END\n");
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
        int natoms; Vec3d pos;
        Quat4d REQ=defaultREQ; 
        char at_name[8]; 
        int npi,ne=0;
        const int nbuf=1024;
        char buff[nbuf]; char* line;
        line = fgets( buff, nbuf, pFile ); // number of atoms
        sscanf( line, "%i", &natoms );
        if(verbosity>1)printf( "natoms %i \n", natoms );
        line = fgets( buff, nbuf, pFile ); // comment, ignore
        int n0 = atoms.size();
        for(int i=0; i<natoms; i++){
            line     = fgets( buff, nbuf, pFile ); // comment, ignore
            //int nret     = sscanf( line,     "%s %lf %lf %lf %lf %i  ",    at_name, &pos.x, &pos.y, &pos.z, &REQ.z, &npi );
            //if(verbosity>1)printf(  ".xyz[%i] %s %lf %lf %lf %lf %i\n", i, at_name,  pos.x,  pos.y,  pos.z,  REQ.z,  npi  );
            int nret     = sscanf( line,     "%s %lf %lf %lf %lf %i  ",    at_name, &pos.x, &pos.y, &pos.z, &REQ.z, &REQ.w, &npi );
            if(verbosity>1)printf(  ".xyz[%i] %s %lf %lf %lf %lf %i\n", i, at_name,  pos.x,  pos.y,  pos.z,  REQ.z,  REQ.w,  npi );
            if( nret < 5 ){ REQ.z=0;  };
            if( nret < 6 ){ REQ.w=0;  };
            if( nret < 7 ){ npi  =-1; };
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                int ityp=it->second;
                if(verbosity>2)printf( " %s -> %i \n", at_name, ityp );
                if(params){
                    params->assignRE( ityp, REQ );
                    ne = params->atypes[ityp].nepair;
                }
                if( noH && (at_name[0]='H') && (at_name[1]='\0')  ) continue;
                insertAtom( it->second, pos, &REQ, npi, ne );
            }
        }
        return atoms.size() - n0;
    }

    /// @brief Load a .mol file (MDL molfile format).
    ///
    /// The file is expected to have a header (first 4 lines), followed by a block of
    /// atom lines and then a block of bond lines. The 5th line (index 4) must contain
    /// the number of atoms and bonds.
    /// Load a MDL .mol file.
    /// The file is assumed to contain:
    ///   - A header with at least 5 lines.
    ///   - The 5th line (index 4) provides the number of atoms and bonds.
    ///   - Next, the atom block (one line per atom) with at least 4 tokens:
    ///         x y z element
    ///   - Then the bond block (one line per bond) with at least 3 tokens:
    ///         bond_id atom1 atom2 ...
    int load_mol(const char* fname, const Vec3d& pos0=Vec3dZero, const Mat3d& rot=Mat3dIdentity ) {
        if(verbosity>0)printf( "MM::Builder.load_mol(%s)\n", fname );
        FILE* fp = fopen(fname, "r");
        if(fp == NULL){
            printf("ERROR in MM::Builder::load_mol(): Cannot open %s\n", fname);
            return -1;
        }
        const int nbuf = 1024;
        char buff[nbuf];
        int lineIndex = 0;
        int na = 0, nb = 0;
        int n0 = atoms.size();

        int ifrag=-1;
        startFragment();
        ifrag   = frags.size()-1;  // frags.size()-1
        
        while (fgets(buff, nbuf, fp)) {
            // The 5th line (lineIndex==4) should contain atom and bond counts.
            if (lineIndex == 4) {
                if (sscanf(buff, "%d %d", &na, &nb) < 2) {
                    printf("Error reading atom/bond counts in file %s\n", fname);
                    fclose(fp);
                    return -1;
                }
            }
            // Process atom lines
            else if (lineIndex > 4 && lineIndex <= 4 + na) {
                Vec3d pos;
                char elem[16] = {0};
                int ret = sscanf(buff, "%lf %lf %lf %15s", &pos.x, &pos.y, &pos.z, elem);
                if(ret < 4){
                    printf("Error reading atom at line %d in %s\n", lineIndex, fname);
                    continue;
                }
                
                Quat4d REQ = defaultREQ;
                int npi = -1;
                int ne = 0;
                
                auto it = atomTypeDict->find(elem);
                if(it != atomTypeDict->end()){
                    int ityp = it->second;
                    if(verbosity>2)printf( " %s -> %i \n", elem, ityp );
                    if(params){
                        params->assignRE(ityp, REQ);
                        ne = params->atypes[ityp].nepair;
                    }
                    //if(noH && (elem[0]=='H') && (elem[1]=='\0')) continue;

                    Vec3d p; rot.dot_to(pos,p); p.add( pos0 ); //p.sub(cog);

                    insertAtom(ityp, p, &REQ, npi, ne);
                }
            }
            // Process bond lines
            else if (lineIndex > 4 + na && lineIndex <= 4 + na + nb) {
                int bond_id, a1, a2;
                int ret = sscanf(buff, "%d %d %d", &bond_id, &a1, &a2);
                if(ret < 3){
                    printf("Error reading bond at line %d in %s\n", lineIndex, fname);
                    continue;
                }
                // Convert 1-based atom indices to 0-based
                Bond bond;
                bond.atoms.i = a1 - 1;
                bond.atoms.j = a2 - 1;
                bonds.push_back(bond);
            }
            lineIndex++;
        }
        fclose(fp);
        if(verbosity>0)printf("Builder::load_mol() read atoms[%d] from %s\n", (int)(atoms.size()-n0), fname);
        return ifrag;
    }

    /// @brief Load a .mol2 file (Tripos mol2 format).
    ///
    /// The file must include at least the following sections:
    ///   - @<TRIPOS>MOLECULE (optional lattice vectors via a comment with "lvs")
    ///   - @<TRIPOS>ATOM
    ///   - @<TRIPOS>BOND (if bonds are present)
    /// Load a Tripos .mol2 file.
    /// The file must contain at least:
    ///   - An @<TRIPOS>MOLECULE section (optionally with a lattice comment starting with '#' and containing "lvs")
    ///   - An @<TRIPOS>ATOM section with lines formatted as:
    ///         atom_id atom_name x y z atom_type ... [charge]
    ///   - Optionally an @<TRIPOS>BOND section with lines formatted as:
    ///         bond_id origin_atom_id target_atom_id bond_type
    bool save_mol2(const char* fname, int ifrag=-1) {
        if(verbosity>0)printf( "MM::Builder.save_mol2(%s)\n", fname );
        FILE* fp = fopen(fname, "w");
        if(fp == NULL){
            printf("ERROR in MM::Builder::save_mol2(): Cannot open %s\n", fname);
            return false;
        }

        // Get atom range for the fragment
        int ia0=0,ia1=atoms.size();
        int ib0=0,ib1=bonds.size();
        if(ifrag>=0){
            ia0=frags[ifrag].atomRange.a;
            ia1=frags[ifrag].atomRange.b;
            ib0=frags[ifrag].bondRange.a;
            ib1=frags[ifrag].bondRange.b;
        }
        int na = ia1-ia0;
        int nb = ib1-ib0;

        // Write MOLECULE section
        fprintf(fp, "@<TRIPOS>MOLECULE\n");
        fprintf(fp, "*****\n");
        fprintf(fp, " %d %d 0 0 0\n", na, nb);
        fprintf(fp, "SMALL\n");
        fprintf(fp, "GASTEIGER\n\n");

        // Write ATOM section
        fprintf(fp, "@<TRIPOS>ATOM\n");
        for(int i=ia0; i<ia1; i++){
            const Atom& atom = atoms[i];
            const char* elem = params->atypes[atom.type].name;
            // Format: atom_id atom_name x y z atom_type subst_id subst_name charge
            fprintf(fp, "%6d %-4s %10.4f %10.4f %10.4f %-4s.3 %5d HOH%1d %10.4f\n",
                i+1-ia0,           // atom_id (1-based)
                elem,             // atom_name
                atom.pos.x,       // x
                atom.pos.y,       // y
                atom.pos.z,       // z
                elem,             // atom_type
                1,                // subst_id
                0,                // subst_name number
                atom.REQ.z        // charge
            );
        }

        // Write BOND section
        fprintf(fp, "@<TRIPOS>BOND\n");
        for(int i=ib0; i<ib1; i++){
            const Bond& bond = bonds[i];
            // Format: bond_id origin_atom_id target_atom_id bond_type
            fprintf(fp, "%6d %6d %6d %1d\n",
                i+1-ib0,                    // bond_id (1-based)
                bond.atoms.i+1-ia0,        // origin_atom_id (1-based)
                bond.atoms.j+1-ia0,        // target_atom_id (1-based)
                bond.type                   // bond_type
            );
        }

        fclose(fp);
        if(verbosity>0)printf("Builder::save_mol2() wrote %d atoms and %d bonds to %s\n", na, nb, fname);
        return true;
    }

    int load_mol2(const char* fname, const Vec3d& pos0=Vec3dZero, const Mat3d& rot=Mat3dIdentity ) {
        if(verbosity>0)printf( "MM::Builder.load_mol2(%s)\n", fname );
        FILE* fp = fopen(fname, "r");
        if(fp == NULL){
            printf("ERROR in MM::Builder::load_mol2(): Cannot open %s\n", fname);
            return -1;
        }
        const int nbuf = 1024;
        char buff[nbuf];

        bool bLoadGroups = false;
        char group_str[nbuf];

        bool inMolecule = false;
        bool inAtom = false;
        bool inBond = false;
        int n0 = atoms.size();
        
        int ifrag=-1;
        startFragment();
        ifrag   = frags.size()-1;  // frags.size()-1
        int iline=0;
        while(fgets(buff, nbuf, fp)) {

            if( buff[0]=='@' ){
                // Check section headers
                if (strncmp(buff, "@<TRIPOS>MOLECULE", 17) == 0) {
                    inMolecule = true;
                    inAtom = false;
                    inBond = false;
                }else             
                if (strncmp(buff, "@<TRIPOS>ATOM", 13) == 0) {
                    inMolecule = false;
                    inAtom = true;
                    inBond = false;
                }else 
                if (strncmp(buff, "@<TRIPOS>BOND", 13) == 0) {
                    inAtom = false;
                    inBond = true;
                }else 
                if (strncmp(buff, "@lvs", 4) == 0) {    // expected line like this:    @lvs 20.0 0.0 0.0    0.0 0.5 0.0   20.0 0.0 0.0
                    sscanf( buff+4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &lvec.xx, &lvec.xy, &lvec.xz, &lvec.yx, &lvec.yy, &lvec.yz, &lvec.zx, &lvec.zy, &lvec.zz );
                    //sscanf( buff+4, "%lf %lf %lf %lf %lf %lf %lf %lf %lf", &lvec.xx, &lvec.xy, &lvec.xz, &lvec.yx, &lvec.yy, &lvec.yz, &lvec.zx, &lvec.zy, &lvec.zz );
                    printf("Builder::load_mol2() lvec loaded\n"); printMat(lvec);
                    bPBC = true;
                }else
                if( strncmp(buff, "@groups", 7) == 0 ){
                    //str2groups(buff+7, nbuf-7 );
                    strcpy(group_str, buff+7);
                    bLoadGroups=true;
                    //exit(0);
                }
                continue;
            }
            //printf( "---Builder::load_mol2() [%i] a%i b%i %s", iline, inAtom, inBond, buff);
            
            // Process atom lines
            if (inAtom) {
                int  atom_id;
                char atom_name[16];
                char atom_type[32];
                Vec3d pos;
                double charge = 0.0;
                
                int ret = sscanf(buff, "%d %15s %lf %lf %lf %31s", &atom_id, atom_name, &pos.x, &pos.y, &pos.z, atom_type);
                if(ret < 6) continue;
                
                // Try to read charge if present
                sscanf(buff, "%*d %*s %*lf %*lf %*lf %*s %*s %*s %lf", &charge);
                printf( "Builder::load_mol2() [%i] atom_id=%d atom_name=%s pos=(%lf %lf %lf) atom_type=%s charge=%lf\n", iline, atom_id, atom_name, pos.x, pos.y, pos.z, atom_type, charge );
                
                // substitute character '.' with '_' in atom_type
                char* dot = strchr(atom_type, '.'); if(dot) { *dot = '_'; }

                Quat4d REQ = defaultREQ;
                REQ.z = charge;  // Set charge from mol2 file
                //int npi = -1;
                //int ne = 0;
                
                auto it = atomTypeDict->find(atom_type);

                if(it == atomTypeDict->end() ){ // if we cannot fined specific type we try to find just the element  like 'S_3' we try to fine 'S'
                    printf( "Builder::load_mol2() [%i] cannot find type `%s` find element `%s` instead \n", iline, atom_type, atom_name );
                    if(dot) { *dot = '\0'; }else{
                        char* dot = strchr(atom_type, '_');
                        if(dot) { *dot = '\0'; }
                    }
                    it = atomTypeDict->find(atom_type);
                    printf( "Builder::load_mol2() [%i] atom_type `%s` atom_name: `%s` \n", iline, atom_type, atom_name );
                }

                if(it != atomTypeDict->end()){
                    int ityp = it->second;
                    if(verbosity>2)printf( " %s -> %i \n", atom_type, ityp );
                    if(params){
                        params->assignRE(ityp, REQ);
                        //ne = params->atypes[ityp].nepair;
                    }
                    Vec3d p; rot.dot_to(pos,p); p.add( pos0 ); //p.sub(cog);
                    int ia = insertAtom(ityp, p, &REQ, -1, 0 );
                    tryAddConfToAtom( ia );
                }else{
                    printf( "ERROR in MM::Builder::load_mol2(): unknown atom type `%s` for atom(%i) line(%i): %s \n", atom_type, atom_id, iline, buff );
                    exit(1);
                }
            }
            
            // Process bond lines
            if (inBond) {
                int bond_id, a1, a2, bond_type;
                //char bond_type[8];
                int ret = sscanf(buff, "%d %d %d %d", &bond_id, &a1, &a2, &bond_type);
                printf( "Builder::load_mol2() [%i] bond_id=%d a1=%d a2=%d type=%i\n", iline, bond_id, a1, a2, bond_type );
                if(ret < 3) continue;
                // Convert 1-based atom indices to 0-based
                Bond bond;
                bond.atoms.i = a1 - 1;
                bond.atoms.j = a2 - 1;
                bond.type = bond_type;
                //bonds.push_back(bond);
                int ib = insertBond( bond, true );
                if( bPBC ){ shortesPBCbond( ib ); }
            }
            iline++;
        } // while( fgets() )

        if( bLoadGroups ){ 
            str2groups(group_str, 100, 1 );
            printAtom2Groups();
        }

        fclose(fp);
        if(verbosity>0)printf("Builder::load_mol2() END read natoms[%d] nbonds[%d] from %s\n", (int)(atoms.size()-n0), (int)(bonds.size()-n0), fname);
        return ifrag;
    }













    inline void natom_def(int& n,int i0)const{ if(n<0){ n=atoms .size()-i0; }; }
    inline void nconf_def(int& n,int i0)const{ if(n<0){ n=confs .size()-i0; }; }
    inline void nbond_def(int& n,int i0)const{ if(n<0){ n=bonds .size()-i0; }; }
    inline void nang_def (int& n,int i0)const{ if(n<0){ n=angles.size()-i0; }; }
    inline void ndih_def (int& n,int i0)const{ if(n<0){ n=dihedrals.size()-i0; }; }

    void export_REQs(Quat4d* REQs, int i0=0, int n=-1)const{
        natom_def(n,i0);
        //_allocIfNull(REQs,n);
        //printf( "export_REQs n %i @REQs=%li \n", n, (long)REQs );
        for(int i=0; i<n; i++){ 
            //printf( "export_REQs[%i] REQ(%g,%g,%g,%g)\n", i, atoms[i0+i].REQ.x,atoms[i0+i].REQ.y,atoms[i0+i].REQ.z,atoms[i0+i].REQ.w  );
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
                //printf( "export_bonds[%i|%i,%i]t(%i,%i) ipe(%i,%i) jpe(%i,%i) \n", i, b.atoms.i, b.atoms.j, atoms[b.atoms.i].type, atoms[b.atoms.j].type,    bi,bei,   bj, bej );
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
    int  finishFragment(int i=-1 ){ if(i<0)i=frags.size()-1; frags[i].finish(           atoms.size(), confs.size(), bonds.size(), angles.size(), dihedrals.size()   ); return i; }

    int insertAtoms( int na, int* atypes, Vec3d* apos, Quat4d* REQs=0, double* qs=0, int* npis=0, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity, const Vec3d& pos0=Vec3dZero ){
        //printf( "MM::Builder::insertAtoms  natoms %i atypes=%li apos=%li REQs=%li qs=%li npis=%li \n", na, (long)atypes, (long)apos, (long)REQs, (long)qs, (long)npis );
        //startFragment();
        //int natom0  = atoms.size();
        //int nbond0  = bonds.size();
        //std::vector<int> atomInds(mol->natoms);
        //std::vector<int> bondInds(mol->nbonds);
        int ncap=0;
        for(int i=0; i<na; i++){
            //printf( "insert Atom[%i] ityp %i %s \n", i, atypes[i], params->atypes[atypes[i]].name );
            Vec3d p;
            Quat4d REQ;
            int ne=0,npi=0;
            int ityp = atypes[i];
            if( ityp==ignoreType ) continue;
            //if( capping_types.contains(ityp)  ){ ncap++; continue; }
            double q=0; if(qs){ q=qs[i]; }
            if(REQs  ){ REQ=REQs[i];                        }
            else      { params->assignRE( ityp, REQ,true ); }
            if(params){ ne = params->atypes[ityp].nepair; npi=params->atypes[ityp].npi; }
            if(npis  ){npi=npis[i];}
            REQ.z=q;
            //printf( "insert Atom[%i] ityp %i REQ(%g,%g,%g) npi,ne %i %i \n", i, ityp, REQ.x, REQ.y, REQ.z, npi, ne  );
            rot.dot_to( apos[i]-pos0,p); p.add( pos );
            insertAtom( ityp, p, &REQ, npi, ne );
        }
        //finishFragment();
        return 0;
    }

    int insertAtoms( const Atoms& mol, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity, const Vec3d& pos0=Vec3dZero ){
        return insertAtoms( mol.natoms, mol.atypes, mol.apos, 0, 0, 0, pos, rot, pos0 );
    }

    Atoms* exportAtoms( int i0=0, int n=-1 ){
        natom_def(n,i0);
        Atoms* out = new Atoms( n, bPBC, true );
        if(bPBC){ *(out->lvec)=lvec; };
        for(int i=0; i<n; i++){
            const Atom& a = atoms[i0+i];
            out->atypes[i] = a.type;
            out->apos  [i] = a.pos;
        }
        return out;
    }

    void insertBonds( int nb, Vec2i* b2as, int* btypes ){
        int natom0  = atoms.size();
        for(int i=0; i<nb; i++){
            int btyp = defaultBond.type;
            if(btypes){ btyp=btypes[i]; }
            bonds.push_back( Bond( btyp, b2as[i] + ((Vec2i){natom0,natom0}), defaultBond.l0, defaultBond.k ) );
        }
    }

    int insertAtomsBonds( int nb, Vec2i* b2as, int* btypes, int na, int* atypes,  Vec3d* apos, Quat4d* REQs=0, double* qs=0, int* npis=0, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity ){
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
        //printf( "multFragPBC() nPBC(%i,%i,%i) ifrag=%i ) \n", nPBC.x, nPBC.y, nPBC.z,  ifrag );
        lvs.print();
        //printSizes();
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
        //printSizes();
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

    void reindexConfs(int* ias, int* ibs){
        // ToDo: this would be much simpler if Confs are part of atoms
        int nc = confs.size();
        std::vector<int> new_inds(nc,-1);
        int inew=0;
        for(int ic=0; ic<nc; ic++){
            AtomConf& c = confs[ic];
            c.iatom = ias[c.iatom];
            if(c.iatom<0) continue;
            atoms[c.iatom].iconf = inew;
            int nb=0;
            for(int j=0; j<N_NEIGH_MAX; j++){
                int jb = c.neighs[j];
                if(jb<0) continue;
                c.neighs[j] = ibs[jb];
                if( c.neighs[j] >=0 ) nb++;
            }
            c.nbond = nb;
            c.updateNtot();
            confs[inew]=c;
            inew++;
        }
        confs.resize(inew);
    }

    void reindexBonds(int* ias, bool bUpdateConfs=true ){
        int nb = bonds.size();
        std::vector<int> new_inds(nb,-1);
        int inew=0;
        for(int ib=0; ib<nb; ib++){
            Vec2i& b = bonds[ib].atoms;
            b.i = ias[b.i];
            b.j = ias[b.j];
            bool bbad = (b.i==-1)||(b.j==-1);
            if(bbad){ 
                new_inds[ib]=-1; 
            }else{ 
                new_inds[ib]=inew; 
                bonds[inew]=bonds[ib];
                inew++; 
            }
        }
        reindexConfs( ias, new_inds.data() );
        bonds.resize(inew);
    }

    void selectAll    (){ selection.clear(); for(int i=0; i<atoms.size(); i++)selection.insert(i); };
    void selectInverse(){ std::unordered_set<int> s(selection.begin(),selection.end());selection.clear(); for(int i=0; i<atoms.size();i++) if( !s.contains(i) )selection.insert(i);  };
    void selectCaping (){ selection.clear(); for(int i=0; i<atoms.size();i++){ if(atoms[i].iconf<0)selection.insert(i); } }

    int selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){
        printf( "Builder::selectRect() p0(%g,%g,%g) p1(%g,%g,%g) \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z );
        Vec3d Tp0,Tp1,Tp;
        //Mat3d rot = (Mat3d)cam.rot;
        rot.dot_to(p0,Tp0);
        rot.dot_to(p1,Tp1);
        _order(Tp0.x,Tp1.x);
        _order(Tp0.y,Tp1.y);
        Tp0.z=-1e+300;
        Tp1.z=+1e+300;
        selection.clear();
        for(int i=0; i<atoms.size(); i++ ){
            rot.dot_to( atoms[i].pos,Tp);
            if( Tp.isBetween(Tp0,Tp1) ){
                selection.insert( i );
            }
        }
        return selection.size();
    }

    int deleteAtoms(int n, int* ias, bool bUpdateBonds=true, bool bUpdateConfs=true, bool bCheckError=true ){
        int na = atoms.size();
        std::vector<int> old_inds(na,-1);
        std::vector<int> new_inds(na,-1);
        // mask removed atoms
        for(int i=0; i<n; i++){ 
            int ia = ias[i];
            //printf("deleteAtoms[%i] ia=%i \n", i,ia);
            if(bCheckError)if( (ia<0)||(ia>=na) ){ printf("Builder::deleteAtoms() ias[%i]=%i out of range 0..atoms.size(%i) => exit()\n", i, ia, na ); exit(0); }
            old_inds[ia]=1; 
        }

        printf("confs before ==== \n"); printAtomConfs();
        //printAtoms();
        // reindex remaining atoms
        int inew = 0;
        int ndel=0;
        for(int i=0; i<na; i++){
            if( old_inds[i]==1 ){
                new_inds[i]=-1;
                ndel++;
                //atoms[inew]=atoms[i];  // Debug
            }else{
                new_inds[i]=inew;
                atoms[inew]=atoms[i];
                inew++;
            }
        }
        //for(int i=0; i<na; i++){    printf("new_inds[%i] == %i \n", i, new_inds[i]);  }
        //printAtoms();
        if(bUpdateBonds) reindexBonds( new_inds.data(), bUpdateConfs );
        atoms.resize(inew);
        printf("confs after ==== \n"); printAtomConfs();
        return ndel;
    }

    inline void changeBondInAtomConf( int ia, int ib, int jb ){ 
        //printf("changeBondInAtomConf ia=%i ib=%i jb=%i \n", ia, ib, jb);
        int ic=atoms[ia].iconf; 
        //printf("changeBondInAtomConf ia=%i ib=%i jb=%i ic=%i\n", ia, ib, jb, ic);
        if(ic>=0){ confs[ ic ].replaceNeigh( ib, jb ); 
            if(jb<0){
                confs[ ic ].sortBonds();
                confs[ ic ].countBonds();
                confs[ ic ].updateNtot();
            }
        }
    }

    void deleteBond(int ib){
        int   jb = bonds.size()-1;
        Vec2i bi = bonds[ib].atoms;
        if(ib!=jb){ // not last bond   
            Vec2i bj = bonds[jb].atoms; 
            changeBondInAtomConf( bi.a, ib, -1 );
            changeBondInAtomConf( bi.b, ib, -1 );
            changeBondInAtomConf( bj.a, jb, ib );
            changeBondInAtomConf( bj.b, jb, ib );
            bonds[ib] = bonds[jb];
        }
        bonds.pop_back();
    }


#ifdef LimitedGraph_h
    bool toLimitedGraph( LimitedGraph<N_NEIGH_MAX>& G, bool bNoCap=true, bool bExitOnError=true, bool bRealloc=true ){
        bool err=false;
        int nc = confs.size();
        int na = atoms.size();
        int nb = bonds.size();
        if( bRealloc ){ if(bNoCap){G.realloc(nc);}else{G.realloc(na);} }
        for(int i=0; i<nb; i++){
            Vec2i b = bonds[i].atoms;
            if(bNoCap){
                int ic1 = atoms[b.a].iconf;
                int ic2 = atoms[b.b].iconf;
                if( (ic1<0)||(ic2<0) ) continue;
                err |= G.addEdge( ic1, ic2 );
                if(err && bExitOnError){ printf( "ERROR in MM::Builder::toLimitedGraph() cannot add bond[%i] neighs are filled nng[%i]=%i nng[%i]=%i  \n => Exit(); \n", i, b.a,G.nneighs[b.a], b.b,G.nneighs[b.b] ); exit(0); };
                //G.print();
            }
        }
        return err;
    }
#endif // LimitedGraph_h


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
                confs[A.iconf].init0(); confs[A.iconf].ne = params->atypes[atyp].nepair; 
            }else{ 
                //addConfToAtom(ia); 
            } 
        }
        int natom0  = atoms.size();
        for(int i=1; i<mol->natoms; i++){
            int ne=0,npi=0;
            Quat4d REQ=mol->REQs[i];
            int ityp = mol->atomType[i];
            if(params){
                //printf( "params \n" );
                params->assignRE( ityp, REQ );
                ne = params->atypes[ityp].nepair;
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

/*
    int substituteMolecule(Molecule* mol, int ia_bone, int ia_frag=0, bool bSwapBond=false ){
        int ib;
        int findBondsToAtom( ia_bone, &ib,false,1);
        int findBondsToAtom( ia_bone, &ib,false,1);
        Bond& b = bonds[ib];
        int ia  = b.atoms.i;
        int ja  = b.atoms.j;
        if(bSwapBond) _swap(ia,ja);
        if(verbosity>0)printf( "substituteMolecule() ib=%i ipivot=%i ja=%i \n", ib, ipivot, ja  );
        Atom& A    = atoms[ja];
        { // set pivot atom
            int atyp = mol->atomType[ipivot];
            A.type = atyp;
            A.REQ  = mol->REQs[ipivot];
            if(A.iconf>=0){ 
                confs[A.iconf].init0(); confs[A.iconf].ne = params->atypes[atyp].nepair; 
            }else{ 
                //addConfToAtom(ia); 
            } 
        }
        int natom0  = atoms.size();
        for(int i=1; i<mol->natoms; i++){
            int ne=0,npi=0;
            Quat4d REQ=mol->REQs[i];
            int ityp = mol->atomType[i];
            if(params){
                //printf( "params \n" );
                params->assignRE( ityp, REQ );
                ne = params->atypes[ityp].nepair;
                REQ.z=mol->REQs[i].z;
            }
            if( mol->npis ) npi=mol->npis[i];
            insertAtom( ityp, p, &REQ, npi, ne );
        }
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
*/

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
        // TBD how ret can be -1???
        if     ( iret<0  ){ return iret; }
        else if( iret==0 ){ bPBC=false;  }
        else if( iret>0  ){ bPBC=true;   }
        //printf("MM::Builder::loadMolTypeXYZ(%s) bPBC=%i \n", fname, bPBC );
        if(params) params->assignREs( mol->natoms, mol->atomType, mol->REQs );
        int ityp = molTypes.size();
        mol2molType[(size_t)mol]=ityp;
        molTypes.push_back(mol);
        return molTypes.size()-1;
    }

    int loadXYZ_Atoms(const char* fname, MMFFparams* params_=0, int iH=-1, bool bCOG=false, const Vec3d& pos=Vec3dZero, const Mat3d& rot=Mat3dIdentity ){
        printf( "MM::Builder::loadXYZ_Atoms(%s) \n", fname );
        if(params_!=0) params=params_;
        //  this is a bit stupid - we allocate and deallocate mol just because we need some temporary storage for calling params->loadXYZ()
        Atoms   mol;
        Quat4d* REQs=0;
        int iret =  params->loadXYZ( fname, mol.natoms, &mol.apos, &REQs, &mol.atypes, 0, &lvec );
        //printAtomConfs(false, true );
        int ifrag=-1;
        if( iret>=0  ){ 
            if( iret==0 ){ bPBC=false; }else{ bPBC=true; }
            //printf("MM::Builder::loadXYZ_Atoms(%s) iH=%i bCOG=%i \n", fname, iH, bCOG );
            //if(params) params->assignREs( mol.natoms, mol.atypes, REQs );
            Vec3d cog=Vec3dZero;
            if(bCOG) cog=mol.getBBcog();
            startFragment();
            ifrag   = frags.size()-1;  // frags.size()-1
            for(int i=0; i<mol.natoms; i++){
                const int ityp = mol.atypes[i];
                if( ityp==iH ) continue;
                int ne=0,npi=0; 
                Quat4d REQ=REQs[i];
                if(params_){  //printf( "params \n" );
                    params->assignRE( ityp, REQ, true );
                    // printf( "MM::Builder::loadXYZ_Atoms() assignRE[%i] REQ(%g,%g,%g)\n", i, REQ.x,REQ.y,REQ.z );
                    ne = params->atypes[ityp].nepair;
                    REQ.z=REQs[i].z;
                }
                //Quat4d  REQi = REQs[i]; REQi.y = sqrt(REQi.y);
                Vec3d p; rot.dot_to(mol.apos[i],p); p.add( pos ); p.sub(cog);
                //atoms.push_back( (Atom){mol.atypes[i], -1, -1, p, REQi } );
                int ia = insertAtom( ityp, p, &REQ, npi, ne );
                atoms[ia].frag = ifrag;
            }
        }
        //printAtomConfs(false, true );
        delete [] REQs;
        //printf( "MM::Builder::loadXYZ_Atoms(%s) ifrag=%i DONE \n", fname, ifrag );
        return ifrag;
    }

      


    int registerRigidMolType( int natoms, Vec3d* pos, Quat4d* REQs, int* atomType ){
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

    int insertFlexibleMolecule_ignorH( Molecule * mol, const Vec3d& pos, const Mat3d& rot, int iH = 1 ){
        printf( "# MM::Builder::insertFlexibleMolecule_ignorH()  natoms %i nbonds %i \n", mol->natoms, mol->nbonds );
        startFragment();
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        std::vector<int> atomInds(mol->natoms);
        std::vector<int> bondInds(mol->nbonds);
        for(int i=0; i<mol->natoms; i++){
            if( mol->atomType[i]==iH ) continue;
            atomInds[i] = atoms.size();
            Quat4d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y);
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
        return finishFragment();
    }

    int insertFlexibleMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot, int ignoreType=-1 ){
        printf( "# MM::Builder::insertFlexibleMolecule()  natoms %i nbonds %i \n", mol->natoms, mol->nbonds );
        startFragment();
        int ifrag   = frags.size()-1;  // frags.size()-1
        int natom0  = atoms.size();
        int nbond0  = bonds.size();
        for(int i=0; i<mol->natoms; i++){
            int ne=0,npi=0;
            Quat4d REQ=mol->REQs[i];
            int ityp = mol->atomType[i];
            if( ityp==ignoreType ) continue;
            if(params){
                //printf( "params \n" );
                params->assignRE( ityp, REQ );
                ne = params->atypes[ityp].nepair;
                REQ.z=mol->REQs[i].z;
            }
            if( mol->npis ) npi=mol->npis[i];
            //printf( "insert Atom[%i] ityp %i REQ(%g,%g,%g) npi,ne %i %i \n", i, ityp, REQ.x, REQ.y, REQ.z, mol->npis[i], ne  );
            Vec3d p; rot.dot_to(mol->pos[i],p); p.add( pos );
            int ia = insertAtom( ityp, p, &REQ, npi, ne );
            atoms[ia].frag = ifrag;
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
        return ifrag;
    }

    int insertRigidMolecule( Molecule * mol, const Vec3d& pos, const Mat3d& rot ){
        printf( "insertRigidMolecule() \n" );
        int natoms0 = atoms.size();
        Quat4d qrot; qrot.fromMatrix(rot);
        int ifrag = frags.size();
        //printf( "insertMolecule mol->natoms %i \n", mol->natoms );
        for(int i=0; i<mol->natoms; i++){
            //Quat4d REQi = Vec3d{1.0,0.03,mol->}; // TO DO : LJq can be set by type
            //atoms.push_back( (Atom){mol->atomType[i],mol->pos[i], LJq } );
            Quat4d  REQi = mol->REQs[i];   REQi.y = sqrt(REQi.y); // REQi.z = 0.0;
            Vec3d  p; rot.dot_to(mol->pos[i],p); p.add( pos );
            atoms.push_back( (Atom){mol->atomType[i], ifrag, -1, p, REQi } );
        }
        //size_t mol_id = (size_t)(mol);
        //frags.push_back( (Fragment){natoms0, atoms.size()-natoms0, pos, qrot, mol}  );
        //frags.push_back( Fragment{ mol, pos, qrot,  {natoms0, atoms.size()-natoms0} } );
        molTypes.push_back(mol);
        size_t mol_id = molTypes.size()-1;
        frags.push_back( Fragment{ mol_id, pos, qrot, Vec2i{natoms0, (int)(atoms.size()-natoms0)}, mol->pos } );
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
            if(noH>0){ return insertFlexibleMolecule_ignorH( mol, pos, rot, noH ); }
            else     { return insertFlexibleMolecule       ( mol, pos, rot      ); }
            //return -1;
        }
    }

    int insertMolecule( int itype                 , const Vec3d& pos, const Mat3d& rot, bool rigid ){ return insertMolecule( molTypes[itype]                 , pos, rot, rigid ); }
    int insertMolecule( const std::string& molName, const Vec3d& pos, const Mat3d& rot, bool rigid ){ return insertMolecule( molTypes[ molTypeDict[molName] ], pos, rot, rigid ); }

    int insertFlexibleMolecule( int itype                 , const Vec3d& pos, const Mat3d& rot, int ignoreType=-1 ){ if(itype<0){ return itype; }else{ return insertFlexibleMolecule( molTypes[itype], pos, rot, ignoreType ); }  }
    int insertFlexibleMolecule( const std::string& molName, const Vec3d& pos, const Mat3d& rot, int ignoreType=-1 ){ return insertFlexibleMolecule( molTypeDict[molName], pos, rot, ignoreType ); }
#endif // Molecule_h


#ifdef SimpleForceField_h
    void toSimpleForceField( SimpleForceField& ff ){
        if(idebug>0) printf( " MMFFbuilder.toSimpleForceField na %li nb %li nA %li nd %li \n", atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        //mmff->deallocate();
        ff.realloc( atoms.size(), bonds.size(), angles.size(), dihedrals.size() );
        for(int i=0; i<atoms.size(); i++){
            ff.apos [i]  = atoms[i].pos;
            if(idebug>0){ printf("[%i]", i); atoms[i].print(); if( atoms[i].iconf>=0){confs[atoms[i].iconf].print();} puts(""); }
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
            if(idebug>0){  printf( "bond[%i] (%i,%i) %g %g | %g %g\n", i, ff.bond2atom[i].i, ff.bond2atom[i].j, ff.bond_l0[i], ff.bond_k[i], b.l0, b.k ); }
            //bondTypes[i]       = bonds[i].type;
        }
        for(int i=0; i<angles.size(); i++){
            const Angle& a  = angles[i];
            ff.ang2bond[i] = a.bonds;
            ff.setAngleParam(i, a.a0, a.k );
            if(idebug>0){  printf( "angle[%i] (%i,%i) (%g,%g) %g\n", i, ff.ang2bond[i].i, ff.ang2bond[i].j, ff.ang_cs0[i].x, ff.ang_cs0[i].y, ff.ang_k[i] ); }
        }
        for(int i=0; i<dihedrals.size(); i++){
            const Dihedral& d  = dihedrals[i];
            ff.tors2bond[i] = d.bonds;
            ff.setTorsParam( i, d.n, d.k );
            if(idebug>0){ printf( "dihedrals[%i] (%i,%i,%i) %i %g\n", i, ff.tors2bond[i].a, ff.tors2bond[i].b, ff.tors2bond[i].c, ff.tors_n[i], ff.tors_k[i] ); }
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
        //if( pbcShifts[i].norm2()>0.1 ){ printf( "PBC-bond[%i]  atoms(%i,%i)  pbcShift(%g,%g,%g) ipb(%i,%i,%i)\n",  bonds[i].atoms.a, bonds[i].atoms.b, pbcShifts[i].x,pbcShifts[i].y,pbcShifts[i].z, bonds[i].ipbc.x,bonds[i].ipbc.y,bonds[i].ipbc.z );  };
    }
}

#ifdef MMFFsp3_h
    //void toMMFFsp3( MMFFsp3& ff, bool bRealloc=true, double K_sigma=1.0, double K_pi=1.0, double K_ecap=0.75, bool bATypes=true ){
    void toMMFFsp3( MMFFsp3& ff, bool bRealloc=true, bool bEPairs=true ){
        int nAmax = atoms.size();
        int nBmax = bonds.size();
        int nCmax = confs.size();
        // {
        //     printf("!!!! WARRNING Debug HACK !!!! Builder::toMMFFsp3(): change array sizes \n");
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
        if(verbosity>0)printf(  "MM::Builder::toMMFFsp3() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        if(bRealloc)ff.realloc( nconf, nb+ne, npi, ncap+ne );
        export_apos     ( ff.apos  ,0,nAmax);
        export_atypes   ( ff.atype ,0,nAmax);
        export_bonds    ( ff.bond2atom, ff.bond_l0, ff.bond_k, ff.bond_kPi,  0,nBmax);
        if ( ff.nneigh_max != N_NEIGH_MAX  ){ printf( "ERROR in MM::Builder.toMMFFsp3() N_NEIGH_MAX(%i) != ff.nneigh_max(%i) ", N_NEIGH_MAX, ff.nneigh_max ); exit(0); } 
        Vec3d hs[N_NEIGH_MAX];
        int ipi=0;
        //int ja=0;
        int ie0=nconf+ncap;
        int iie = 0;

        for(int i=0; i<ff.nnode*4;  i++){ ff.neighs[i]=-1; };

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
                int*    ngs  = ff.neighs + ia*N_NEIGH_MAX;
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
        if(verbosity>0)printf(  "MM::Builder::toMMFFsp3() DONE \n"  );
    }
#endif // MMFFmini_h

void makeNeighs( int*& neighs, int perAtom ){
    int na = atoms.size();
    int ntot= na*perAtom;
    _allocIfNull( neighs, ntot );
    //printf( "MM::Builder::makeNeighs() ntot=%i  ^neighs=%li \n", ntot, neighs );
    for(int i=0;i<ntot; i++){ neighs[i]=-1; }; // back neighbors
    for(int ia=0; ia<na; ia++ ){
        //printf( "MM::Builder::makeNeighs()[%i] \n", ia );
        const Atom& A =  atoms[ia];
        if(A.iconf>=0){
            const AtomConf& conf = confs[A.iconf];
            for(int k=0; k<conf.nbond; k++){
                int ib        = conf.neighs[k];
                if(ib<0) continue;
                const Bond& B = bonds[ib];
                int ja        = B.getNeighborAtom(ia);
                neighs[ ia*perAtom + k ] = ja;
                int jc = atoms[ja].iconf;
                //printf( "makeNeighs[%i|%i] ja %i jc %i \n", ia, k, ja, jc );
                if( jc==-1 ){ neighs[ ja*perAtom ]=ia; }
            }
        }
    }
    //for(int ia=0; ia<na; ia++){ printf( "neigh[%i](%i,%i,%i,%i)\n", ia, neighs[ia*perAtom],neighs[ia*perAtom+1],neighs[ia*perAtom+2],neighs[ia*perAtom+3] ); };
}

#ifdef UFF_h
void toUFF( UFF& ff, bool bRealloc=true ){
    printf( "MM::Builder::toUFF() \n" );

    int natoms      = atoms.size();
    int nbonds      = bonds.size();
    int nangles     = angles.size();
    int ndihedrals  = dihedrals.size();
    int ninversions = inversions.size();
    if(verbosity>0) printf(  "MM::Builder::toUFF() natoms %i nbonds %i nangles %i ndihedrals %i ninversions %i\n", natoms, nbonds, nangles, ndihedrals, ninversions );
    if(bRealloc) ff.realloc( natoms, nbonds, nangles, ndihedrals, ninversions );
    
    for(int i=0; i<nbonds; i++ ){
        const Bond& B = bonds[i];
        ff.bonAtoms[i]  = Vec2i{B.atoms.x, B.atoms.y};
        ff.bonParams[i] = Vec2d{B.k, B.l0};
    }
    for(int i=0; i<nangles; i++ ){
        const Angle& A = angles[i];
        ff.angAtoms[i]  = Vec3i{A.atoms.x, A.atoms.y, A.atoms.z};
        ff.angParams[i] = double5{A.k, A.C0, A.C1, A.C2, A.C3};
    }
    for(int i=0; i<ndihedrals; i++ ){
        const Dihedral& D = dihedrals[i];
        ff.dihAtoms[i]  = Quat4i{D.atoms.x, D.atoms.y, D.atoms.z, D.atoms.w};
        ff.dihParams[i] = Vec3d{D.k, D.d, D.n};
    }
    for(int i=0; i<ninversions; i++ ){
        const Inversion& I = inversions[i];
        ff.invAtoms[i]  = Quat4i{I.atoms.x, I.atoms.y, I.atoms.z, I.atoms.w};
        ff.invParams[i] = Quat4d{I.k, I.C0, I.C1, I.C2};
    }

    for(int ia=0; ia<natoms; ia++ ){
        const Atom& A =  atoms[ia];
        ff.atypes[ia] = A.type;
        ff.apos [ia]  = A.pos; 
        ff.neighs[ia]=Quat4i{-1,-1,-1,-1};
        if(A.iconf>=0){
            AtomConf& conf = confs[A.iconf];
            int* ngs  = ff.neighs[ia].array;
            for(int k=0; k<conf.nbond; k++){
                int ib = conf.neighs[k];
                const Bond& B = bonds[ib];
                int ja = B.getNeighborAtom(ia);
                ngs[k] = ja;
            }
        }
    } 

    ff.makeNeighBs();         
    ff.bakeAngleNeighs();     
    ff.bakeDihedralNeighs();  
    ff.bakeInversionNeighs(); 
    ff.mapAtomInteractions(); 
    ff.printSizes();          
 
    ff.bPBC = bPBC;
    if(verbosity>0)printf("MM::Builder::toUFF DONE\n");

}
#endif // UFF_h

#ifdef MMFFsp3_loc_h
void toMMFFsp3_loc( MMFFsp3_loc& ff, bool bRealloc=true, bool bEPairs=true, bool bUFF=false ){

        //double c0s[3]{-0.33333,-0.5,-1.0}; // cos(angle)   sp1 sp2 sp3
        double ang0s[3]{ 109.5 *M_PI/180.0, 120.0*M_PI/180.0, 180.0*M_PI/180.0 }; // cos(angle)   sp1 sp2 sp3

        int nAmax = atoms.size();
        int nCmax = confs.size();
        int npi,ne; ne=countPiE( npi, 0,nCmax );
        if(!bEPairs) ne=0;
        int nconf = nCmax;
        int ncap  = nAmax - nconf;
        int nb    = bonds.size();
        if(verbosity>0)printf(  "MM::Builder::toMMFFsp3_loc() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        int ntors=0; if(ff.bTorsion) ntors=dihedrals.size();
        if(bRealloc)ff.realloc( nconf, ncap+ne, ntors );
        Vec3d hs[4];
        int ie0=nconf+ncap;
        int iie = 0;

        if(ff.bTorsion){
            for(int i=0; i<ff.ntors; i++){
                const Dihedral& dih=dihedrals[i];
                ff.tors2atom [i]=dih.atoms;
                //ff.torsParams[i]=Quat4d{cos(dih.a0),sin(dih.a0),dih.k,dih.n};
                ff.torsParams[i]=Quat4d{cos(dih.a0),sin(dih.a0),dih.k, (double)dih.n}; 
            }
        }

        //params->printAtomTypeDict();
        int etyp=-1;  if(params) etyp=params->atomTypeDict["E"];
        for(int i=0; i<ff.natoms; i++){ ff.neighs[i]=Quat4i{-1,-1,-1,-1}; };
        for(int i=0; i<ff.nnode;  i++){ ff.bLs[i]=Quat4dZero, ff.bKs[i]=Quat4dZero, ff.Ksp[i]=Quat4dZero, ff.Kpp[i]=Quat4dZero; }; // back neighbors
        for(int ia=0; ia<nAmax; ia++ ){
            const Atom& A =  atoms[ia];
            ff.apos  [ia] = A.pos;
            ff.atypes[ia] = A.type;
            AtomType& atyp = params->atypes[A.type];

            if(A.iconf>=0){
                // Prepare params and orientation
                AtomConf& conf = confs[A.iconf];
                //printf( "Builder::toMMFFsp3_loc() [%i] ", ia ); conf.print(); printf("\n");
                int npi_neigh = countAtomPiNeighs(ia);
                //assignSp3Params( A.type, conf.nbond, conf.npi, conf.ne, npi_neigh, ff.NeighParams[ia] );

                // setup atom (onsite)
                // ff.apars[ia].x = c0s[conf.npi];    // ssC0  // cos(angle) for angles (sigma-siamg)
                // ff.apars[ia].y = 1.0;              // ssK   // stiffness  for angles
                // ff.apars[ia].z = 0.0;              // piC0  // stiffness  for orthogonalization sigma-pi 
                //printf( "MM::Builder::toMMFFsp3_loc()[%i] conf.npi=%i \n", ia, conf.npi );

                if( conf.npi>2 ){ printf("ERROR in MM::Builder::toMMFFsp3_loc(): atom[%i].conf.npi(%i)>2 => exit() \n", ia, conf.npi); printAtomConf(ia); exit(0); }
                //double ang0 = ang0s[conf.npi];
                double ang0   = atyp.Ass*deg2rad;
                ang0 *= 0.5;
                ff.apars[ia].x = cos(ang0);    // ssC0    // cos(angle) for angles (sigma-siamg)
                ff.apars[ia].y = sin(ang0);
                //ff.apars[ia].z = 1.0;          // ssK     // stiffness  for angles
                //ff.apars[ia].w = 0;            // piC0    // angle0 for orthogonalization sigma-pi 

                ff.apars[ia].z = atyp.Kss*4.0;   // ssK     // stiffness  for angles    ... ToDo: check if K/4 or K*4
                ff.apars[ia].w = sin(atyp.Asp*deg2rad);  // piC0    // angle0 for orthogonalization sigma-pi 

                //printf( "atom[%i] npi(%i)=> angle %g cs(%g,%g) \n", ia, conf.npi, ang0*180./M_PI, ff.apars[ia].x, ff.apars[ia].y  ); 

                // setup ff neighbors                
                int*     ngs  = ff.neighs[ia].array;
                double*  bL   = ff.bLs[ia].array;
                double*  bK   = ff.bKs[ia].array;
                double*  Ksp  = ff.Ksp[ia].array;
                double*  Kpp  = ff.Kpp[ia].array;
                //printf( "BEFOR atom[%i] ngs{%i,%i,%i,%i}\n", ia, ngs[0],ngs[1],ngs[2],ngs[3] );

                // --- Generate Bonds
                //printf( "ia[%i nbond=%i \n", ia, conf.nbond  );
                for(int k=0; k<conf.nbond; k++){
                    int ib = conf.neighs[k];
                    const Bond& B = bonds[ib];
                    //{ int ti=atoms[B.atoms.a].type; int tj=atoms[B.atoms.b].type;   printf( "ia,ib[%4i,%4i] l0=%7.3f k=%7.2f ts(%3i,%3i) %s-%s \n", ia, ib, B.l0, B.k, ti,tj, params->atypes[ti].name, params->atypes[tj].name ); }
                    int ja = B.getNeighborAtom(ia);
                    //{ printf( "ia=%4i ja=%4i ib=%i \n", ia, ja, ib ); }
                    const Atom& Aj =  atoms[ja];
                    AtomType& jtyp = params->atypes[Aj.type];
                    hs[k]  = atoms[ja].pos - A.pos;
                    hs[k].normalize();
                    ngs[k] = ja;
                    bL [k]=B.l0;
                    bK [k]=B.k;
                    if(bUFF){
                        Vec2d bLK = assignBondParamsUFF( ib );
                        bL [k]=bLK.x;
                        bK [k]=bLK.y;
                    }
                    //Ksp[k]=0;
                    if( (conf.npi>0)||(conf.ne>0) ){ Ksp[k]= atyp.Ksp;}else{ Ksp[k]=0; }
                    int nej  = getAtom_ne (ja);
                    int npij = getAtom_npi(ja);
                    Kpp[k]   = sqrt( atyp.Kpp * jtyp.Kpp );
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
                        //Ksp[k]=0;                angs[iang].x = cos(ang0);
                        if( conf.npi>0 ){ Ksp[k]=atyp.Ksp; }else{ Ksp[k]=0; }  // only electron on atoms without pi-orbital are conjugted with pi-orbitas on neighboring atoms
                        iie++;
                    }
                }
                ff.pipos[ia] = hs[N_NEIGH_MAX-1]; // Pi orientation
                //if(verbosity>1){ for(int k=0; k<N_NEIGH_MAX; k++ ){ printf( " %i,", ngs[k] ); }; printf( "] \n" ); }
                //printf( "AFTER atom[%i] ngs{%i,%i,%i,%i}\n", ia, ngs[0],ngs[1],ngs[2],ngs[3] );
            } // if(A.iconf>=0){
        }

        ff.bPBC = bPBC;
        ff.makeBackNeighs();
        //if( bPBC ){ ff.initPBC(); updatePBC( ff.pbcShifts ); }
        //if(verbosity>0)
        printf(  "MM::Builder::toMMFFsp3_loc() DONE \n"  );
    }

void assignAnglesMMFFsp3( MMFFsp3_loc& ff, bool bUFF=false ){
    printf( "MM::Builder::assignAnglesMMFFsp3()\n" );
    for(int ia=0; ia<ff.nnode; ia++ ){
        printf( "assignAnglesMMFFsp3[ia=%i]\n", ia );
        int iat = ff.atypes[ia];
        AtomType& atyp = params->atypes[iat];
        int*    ngs = ff.neighs[ia].array;
        Vec3d* angs = ff.angles    + ia*6;
        double*  bL   = ff.bLs[ia].array;
        double*  bK   = ff.bKs[ia].array;
        int iang=0;
        for(int i=0; i<3; i++){
            int ing = ngs[i];
            if(ing<0) break;
            int it = ff.atypes[ing];
            for(int j=i+1; j<4; j++){
                int jng  = ngs[j];
                if(jng<0) break; 
                int jt = ff.atypes[jng];
                double ang0;
                double k;
                if(bUFF){
                    double ang0    = atyp.Ass*deg2rad;
                    double Kss_uff = params->assignAngleParamUFF( iat, it, jt, bL[i], bL[j] );
                }else{
                    AngleType* ang = params->getAngleType( it, iat, jt, true, true );
                    if(ang==0){ printf("ERROR in MM::Builder::assignAnglesMMFFsp3(ia=%i,%i,%i) cannot find angle type(%i,%i,%i)(%s,%s,%s) =>Exit()\n", ia,i,j, it, iat, jt,  params->atypes[it].name,params->atypes[iat].name, params->atypes[jt].name ); params->printAngleTypesDict();   exit(0); };
                    ang0 = ang->angle0;
                    k    = ang->stiffness;
                }
                ang0 *= 0.5;
                angs[iang].x = cos(ang0);
                angs[iang].y = sin(ang0);
                angs[iang].z = k;
                iang++; 
            }
        }
    }
}


#endif // MMFFsp3_loc_h




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
        if(verbosity>0)printf(  "MM::Builder::toMMFFf4() nconf %i ncap %i npi %i ne %i \n", nconf, ncap, npi, ne  );
        if(bRealloc)ff.realloc( nconf, ncap+ne );
        Vec3d hs[N_NEIGH_MAX];
        int ie0=nconf+ncap;
        int iie = 0;

        for(int i=0; i<ff.nnode;  i++){ ff.neighs[i]=Quat4i{-1,-1,-1,-1};  ff.bLs[i]=Quat4fZero, ff.bKs[i]=Quat4fZero, ff.Ksp[i]=Quat4fZero, ff.Kpp[i]=Quat4fZero; }; // back neighbors

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
                if( conf.npi>2 ){ printf("ERROR in MM::Builder::toMMFFsp3_loc(): atom[%i].conf.npi(%i)>2 => exit() \n", ia, conf.npi); printAtomConf(ia); exit(0); }
                double ang0 = ang0s[conf.npi];
                ang0 *= 0.5;
                ff.apars[ia].x = cos(ang0);    // ssCos0  // cos(angle) for angles (sigma-siamg)
                ff.apars[ia].y = sin(ang0);    // ssSin0
                ff.apars[ia].z = 1.0;          // ssK     // stiffness  for angles
                ff.apars[ia].w = 0.0;          // piCos0  // stiffness  for orthogonalization sigma-pi 
                //printf( "atom[%i] npi(%i)=> angle %g cs(%g,%g) \n", ia, conf.npi, ang0*180./M_PI, ff.apars[ia].x, ff.apars[ia].y  ); 

                // setup ff neighbors
                
                int*    ngs  = ff.neighs[ia].array;
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
        if(verbosity>0)printf(  "MM::Builder::toMMFFf4() DONE \n"  );
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
            mmff->REQ  [i]  = atoms[i].REQ;
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
        for(int i=0; i<atoms.size(); i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            na++;
            ne += params[ ityp ].ne;
            printf( "[%i] ityp %i ne %i  %i \n", i, ityp, params[ityp].ne, ne );
        }
        printf( "na %i ne %i | %i \n", atoms.size(), ne, bonds.size()*2 );
        ff.realloc( na, ne );
        for(int i=0; i<na; i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            //printf( "[%i] ityp %i \n" );
            ff.apos [i]  = atoms[i].pos;
            //ff.aQ [i]  = params[ ityp ].ne; // ToDo
            ff.aPars[i]  = EFF::default_AtomParams[ityp];
        }
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
        for(int i=0; i<atoms.size(); i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            na++;
            ne += params[ ityp ].ne;
            printf( "[%i] ityp %i ne %i  %i \n", i, ityp, params[ityp].ne, ne );
        }
        printf( "na %i ne %i | %i \n", atoms.size(), ne, bonds.size()*2 );
        ff.realloc( na, ne );
        for(int i=0; i<na; i++){
            int ityp = atoms[i].type;
            if( ityp==capAtomEpair.type || ityp==capAtomPi.type ) continue;
            //printf( "[%i] ityp %i \n" );
            ff.apos [i]  = atoms[i].pos;
            ff.aQ   [i]  = params[ ityp ].ne; // ToDo
        }
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

    void assignUFFtypes_trivial( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& Ai = atoms[ia];
            AtomConf& Ci = confs[Ai.iconf];
            if( params->atypes[Ai.type].name[0] == 'H' ){ // all hydrogens
                Ai.type = params->getAtomType("H_", true);
                set_atom[ia] = true;
                int ja = neighs[ia*4];
                int ib = getBondByAtoms( ia, ja );
                BOs[ib] = 1.0;
                BOs_int[ib] = 1;
                set_bond[ib] = true;
            } else {
                AtomConf& Ci = confs[Ai.iconf];
                if( params->atypes[Ai.type].name[0] == 'C' && Ci.nbond == 4 ){ // all sp3 carbons
                    Ai.type = params->getAtomType("C_3", true);
                    set_atom[ia] = true;
                    for(int in=ia*4; in<ia*4+4; in++){
                        int ja = neighs[in];
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.0;
                        BOs_int[ib] = 1;
                        set_bond[ib] = true;
                    }    
                }else if( params->atypes[Ai.type].name[0] == 'N' && Ci.nbond == 3 ){ // all sp3 nitrogens
                    // not setting ufftype and final bond orders, as some may be resonant
                    for(int in=ia*4; in<ia*4+3; in++){
                        int ja = neighs[in];
                        int ib = getBondByAtoms( ia, ja );
                        BOs_int[ib] = 1;
                    }    
                }else if( params->atypes[Ai.type].name[0] == 'N' && Ci.nbond == 1 ){ // all sp1 nitrogens
                    Ai.type = params->getAtomType("N_1", true);
                    set_atom[ia] = true;
                    int ja = neighs[ia*4];
                    int ib = getBondByAtoms( ia, ja );
                    BOs[ib] = 3.0;
                    BOs_int[ib] = 3;
                    set_bond[ib] = true;
                    Atom& Aj = atoms[ja];
                    // nitrile group
                    if( params->atypes[Aj.type].name[0] == 'C' ){ 
                        Aj.type = params->getAtomType("C_1", true);
                        set_atom[ja] = true;
                        for(int jn=ja*4; jn<ja*4+2; jn++){
                            int ka = neighs[jn];
                            if( ka == ia ) continue;
                            int jb = getBondByAtoms( ja, ka );
                            BOs[jb] = 1.0;
                            BOs_int[jb] = 1;
                            set_bond[jb] = true;
                        }
                    }
                }else if( params->atypes[Ai.type].name[0] == 'O' && Ci.nbond == 2 ){ // all sp3 oxygens
                    // not setting ufftype and final bond orders, as some may be resonant
                    for(int in=ia*4; in<ia*4+2; in++){
                        int ja = neighs[in];
                        int ib = getBondByAtoms( ia, ja );
                        BOs_int[ib] = 1;
                    }    
                }else if( params->atypes[Ai.type].name[0] == 'O' && Ci.nbond == 1 ){ // all sp2 oxygens
                    // not setting ufftype and final bond orders, as some may be resonant
                    int ja = neighs[ia*4];
                    int ib = getBondByAtoms( ia, ja );
                    BOs_int[ib] = 2;
                }
            }
        }

    }

    void assignUFFtypes_nitro( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& Ai = atoms[ia];
            if( params->atypes[Ai.type].name[0] == 'H' ) continue;
            AtomConf& Ci = confs[Ai.iconf];
            if( params->atypes[Ai.type].name[0] == 'N' && Ci.nbond == 3 ){ 
                int n = 0;
                for(int in=ia*4; in<ia*4+3; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'O' ){n++;};
                }
                if( n == 2 ){
                    Ai.type = params->getAtomType("N_R", true);
                    set_atom[ia] = true;
                    for(int in=ia*4; in<ia*4+3; in++){
                        int ja = neighs[in];
                        Atom& Aj = atoms[ja];
                        if ( params->atypes[Aj.type].name[0] == 'O' ){
                            Aj.type = params->getAtomType("O_R", true);
                            set_atom[ja] = true;
                            int ib = getBondByAtoms( ia, ja );
                            BOs[ib] = 1.5;
                            set_bond[ib] = true;
                        }else{
                            int ib = getBondByAtoms( ia, ja );
                            BOs[ib] = 1.0;
                            set_bond[ib] = true;
                        }
                    }
                }
            }
        }

    }

    bool assignUFFtypes_checkall( int* neighs, int* BOs_int ){

        // check saturation
        for (; true; ) {
            bool changed = false;
            for(int ia=0; ia<atoms.size(); ia++){
                bool found = false;
                Atom& A = atoms[ia];
                if( params->atypes[A.type].name[0] == 'H' ) continue;
                AtomConf& C = confs[A.iconf];
                // skip already set atoms
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    int ib = getBondByAtoms( ia, ja );
                    if( BOs_int[ib] < 0 ){ found = true; break; }
                }
                if( !found ) continue;
                // look for atoms that misses only one bond (or only single bonds)
                int nset = 0;
                int val = 0;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    int ib = getBondByAtoms( ia, ja );
                    if( BOs_int[ib] > 0 ){ nset++; val+=BOs_int[ib]; }
                }
                // if found, set it
                if( ( nset == C.nbond - 1 ) || ( C.nbond - nset == params->atypes[A.type].valence - val ) ){
                    changed = true;
                    if( nset == C.nbond-1 ){
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ja = neighs[in];
                            int ib = getBondByAtoms( ia, ja );
                            if( BOs_int[ib] < 0 ){ BOs_int[ib] = params->atypes[A.type].valence - val; break; }
                        }
                    }else if( C.nbond - nset == params->atypes[A.type].valence - val ){
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ja = neighs[in];
                            int ib = getBondByAtoms( ia, ja );
                            if( BOs_int[ib] < 0 ){ BOs_int[ib] = 1; }
                        }
                    }
                }
            }
            if (!changed) break;
        }

        // check if we are done
        bool done = true;
        // check if the bond orders are all set
        for (int ib = 0; ib < bonds.size(); ib++) {
            if ( BOs_int[ib] < 0 ) { done = false; break; }
        }
        return done;

    }

    bool assignUFFtypes_checkatom( int ia, int* neighs, int* BOs_int ){

        bool exclude = true;
        Atom& A = atoms[ia];
        AtomConf& C = confs[A.iconf];
        for(int in=ia*4; in<ia*4+C.nbond; in++){
            int ja = neighs[in];
            int ib = getBondByAtoms( ia, ja );
            if( BOs_int[ib] < 0 ){ exclude = false; break; }
        }
        return exclude;

    }

    void assignUFFtypes_treewalk( int* neighs, int* BOs_int ){

        // see if the problem is trivial
        bool done = assignUFFtypes_checkall ( neighs, &BOs_int[0] );
        if ( done ) return;

        // --- explore the graph
        // NB: at this point, only sp2 C, sp2 N and sp C (both triple bond and 2x double bond) are left to be assigned
        // TBD use arrays instead of vectors
        std::vector<int> BOs_int_tmp0(bonds.size());
        std::vector<int> BOs_int_tmp1(bonds.size());
        std::vector<int> BOs_int_tmp2(bonds.size());
        std::vector<int> BOs_int_tmp3(bonds.size());
        for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp0[ib] = BOs_int[ib]; }

        // pick an atom to attempt assigning a double bond
        for(int ia1=0; ia1<atoms.size(); ia1++){
            Atom& A1 = atoms[ia1];
            if( params->atypes[A1.type].name[0] == 'H' ) continue;
            if ( assignUFFtypes_checkatom( ia1, neighs, &BOs_int[0] ) ) continue;
            // pick a a neighbor
            AtomConf& C1 = confs[A1.iconf];
            for(int in1=ia1*4; in1<ia1*4+C1.nbond; in1++){
                int ib1 = getBondByAtoms( ia1, neighs[in1] );
                if( BOs_int_tmp0[ib1] > 0 ) continue;
                // try to set it up as a double bond
                for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp1[ib] = BOs_int_tmp0[ib]; }
                BOs_int_tmp1[ib1] = 2;
                // see if this solves the problem
                if ( assignUFFtypes_checkall ( neighs, &BOs_int_tmp1[0] ) ){ 
                    for (int ib = 0; ib < bonds.size(); ib++) { BOs_int[ib] = BOs_int_tmp1[ib]; } 
                    return; 
                }

                // otherwise, continue with the second nested attempt
                // pick an atom to attempt assigning a double bond
                for(int ia2=0; ia2<atoms.size(); ia2++){
                    Atom& A2 = atoms[ia2];
                    if( params->atypes[A2.type].name[0] == 'H' ) continue;
                    if ( assignUFFtypes_checkatom( ia2, neighs, &BOs_int[0] ) ) continue;
                    // pick a a neighbor
                    AtomConf& C2 = confs[A2.iconf];
                    for(int in2=ia2*4; in2<ia2*4+C2.nbond; in2++){
                        int ib2 = getBondByAtoms( ia2, neighs[in2] );
                        if( BOs_int_tmp1[ib2] > 0 ) continue;
                        // try to set it up as a double bond
                        for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp2[ib] = BOs_int_tmp1[ib]; }
                        BOs_int_tmp2[ib2] = 2;
                        // see if this solves the problem
                        if ( assignUFFtypes_checkall ( neighs, &BOs_int_tmp2[0] ) ){ 
                            for (int ib = 0; ib < bonds.size(); ib++) { BOs_int[ib] = BOs_int_tmp2[ib]; } 
                            return; 
                        }

                        // otherwise, continue with the third nested attempt
                        // pick an atom to attempt assigning a double bond
                        for(int ia3=0; ia3<atoms.size(); ia3++){
                            Atom& A3 = atoms[ia3];
                            if( params->atypes[A3.type].name[0] == 'H' ) continue;
                            if ( assignUFFtypes_checkatom( ia3, neighs, &BOs_int[0] ) ) continue;
                            // pick a a neighbor
                            AtomConf& C3 = confs[A3.iconf];
                            for(int in3=ia3*4; in3<ia3*4+C3.nbond; in3++){
                                int ib3 = getBondByAtoms( ia3, neighs[in3] );
                                if( BOs_int_tmp2[ib3] > 0 ) continue;
                                // try to set it up as a double bond
                                for (int ib = 0; ib < bonds.size(); ib++) { BOs_int_tmp3[ib] = BOs_int_tmp2[ib]; }
                                BOs_int_tmp3[ib3] = 2;
                                // see if this solves the problem
                                if ( assignUFFtypes_checkall ( neighs, &BOs_int_tmp3[0] ) ){ 
                                    for (int ib = 0; ib < bonds.size(); ib++) { BOs_int[ib] = BOs_int_tmp3[ib]; } 
                                    return; 
                                }
                            }
                        }
                    }
                }
            }
        }

        printf("ERROR TREEWALK: tree walk failed\n");
        printf("STOP\n");
        exit(0);

    }

    void assignUFFtypes_simplerule( double tol, int* neighs, double* BOs, bool* set_atom, bool* set_bond ){

        // simplest-est rule: atoms with two sp2 neighbors (or heteroatoms) are resonant
        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            if( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if( ( params->atypes[A.type].name[0] == 'C' && C.nbond == 3 ) || // C_2
                ( params->atypes[A.type].name[0] == 'N' && C.nbond > 1 ) ||  // N_3 or N_2
                ( params->atypes[A.type].name[0] == 'O' && C.nbond == 2 ) ){ // O_3
                int n = 0;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( ( params->atypes[Aj.type].name[0] == 'C' && Cj.nbond == 3 ) || // C_2
                        ( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond > 1 ) ||  // N_3 or N_2
                        ( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 2 ) ){ n++; }
                }
                if( n > 1 ){
                    // set atoms
                    if ( set_atom[ia] ) {
                        if ( params->atypes[A.type].name[2] != 'R' ) printf("WARNING SIMPLERULE: atom %i would be set to resonant but it has already a type of %s\n", ia+1, params->atypes[A.type].name);
                        continue;
                    }
                    if ( params->atypes[A.type].name[0] == 'C' ){ 
                        A.type = params->getAtomType("C_R", true); 
                    }else if ( params->atypes[A.type].name[0] == 'N' ){ 
                        A.type = params->getAtomType("N_R", true);
                    }else if ( params->atypes[A.type].name[0] == 'O' ){
                        A.type = params->getAtomType("O_R", true);
                    }else{
                        printf("ERROR SIMPLERULE: atom %i type %i not recognized\n", ia, A.type);
                        printf("STOP\n");
                        exit(0);
                    }
                    set_atom[ia] = true;
                    // set bonds
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ja = neighs[in];
                        Atom& Aj = atoms[ja];
                        if( params->atypes[Aj.type].name[0] == 'H' ) continue;
                        AtomConf& Cj = confs[Aj.iconf];
                        if( ( params->atypes[Aj.type].name[0] == 'C' && Cj.nbond == 3 ) || // C_2
                            ( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond > 1 ) ||  // N_3 or N_2
                            ( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 2 ) ){ // O_3
                            int ib = getBondByAtoms( ia, neighs[in] );
                            if ( set_bond[ib] ){
                                if ( abs(BOs[ib]-1.5) > tol ) printf("WARNING SIMPLERULE: bond %i between atoms %i and %i would be set to 1.5 but it has already a bond order of %g\n", ib+1, ia+1, ja+1, BOs[ib]);
                                continue;
                            }
                            BOs[ib] = 1.5;
                            set_bond[ib] = true;
                        }
                    }
                }
            }
        }
    }

    void assignUFFtypes_findrings( double tol, int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){
        bool   explored[atoms.size()];
        for(int ia=0; ia<atoms.size(); ia++){explored[ia]=false;}
        std::vector<int>  nb(8);
        for (; true; ) {
            bool changed = false;
            for(int ia=0; ia<atoms.size(); ia++){
                nb[1] = ia; // 1st atom
                Atom& A1 = atoms[nb[1]];
                if ( params->atypes[A1.type].name[0] == 'H' ) continue;
                AtomConf& C1 = confs[A1.iconf];
                for(int in1=nb[1]*4; in1<nb[1]*4+C1.nbond; in1++){
                    nb[2] = neighs[in1]; // 2nd atom
                    Atom& A2 = atoms[nb[2]];
                    if ( params->atypes[A2.type].name[0] == 'H' ) continue;
                    AtomConf& C2 = confs[A2.iconf];
                    for(int in2=nb[2]*4; in2<nb[2]*4+C2.nbond; in2++){
                        nb[3] = neighs[in2]; // 3rd atom
                        if ( nb[3] == nb[1] ) continue;
                        Atom& A3 = atoms[nb[3]];
                        if ( params->atypes[A3.type].name[0] == 'H' ) continue;
                        AtomConf& C3 = confs[A3.iconf];
                        for(int in3=nb[3]*4; in3<nb[3]*4+C3.nbond; in3++){
                            nb[4] = neighs[in3]; // 4th atom
                            if ( nb[4] == nb[2] ) continue;
                            Atom& A4 = atoms[nb[4]];
                            if ( params->atypes[A4.type].name[0] == 'H' ) continue;
                            AtomConf& C4 = confs[A4.iconf];
                            for(int in4=nb[4]*4; in4<nb[4]*4+C4.nbond; in4++){
                                nb[5] = neighs[in4]; // 5th atom
                                if ( nb[5] == nb[3] ) continue;
                                Atom& A5 = atoms[nb[5]];
                                if ( params->atypes[A5.type].name[0] == 'H' ) continue;
                                AtomConf& C5 = confs[A5.iconf];
                                for(int in5=nb[5]*4; in5<nb[5]*4+C5.nbond; in5++){
                                    nb[6] = neighs[in5]; // 6th atom
                                    if ( nb[6] == nb[4] ) continue;

                                    // found 5-member ring
                                    if ( nb[6] == nb[1] ) {
                                        if ( explored[nb[1]] && explored[nb[2]] && explored[nb[3]] && explored[nb[4]] && explored[nb[5]] ) continue;
                                        nb[0] = nb[5];
                                        // count double bonds & heteroatoms
                                        if ( assignUFFtypes_checkaroma ( 5, &nb[0], &BOs[0], &BOs_int[0], tol ) ) {
                                            changed = true;
                                            for(int i=1; i<6; i++){explored[nb[i]] = true;}
                                            assignUFFtypes_setaroma ( 5, &nb[0], &neighs[0], &BOs[0], &set_atom[0], &set_bond[0], tol );
                                        }
                                    }

                                    Atom& A6 = atoms[nb[6]];
                                    if ( params->atypes[A6.type].name[0] == 'H' ) continue;
                                    AtomConf& C6 = confs[A6.iconf];
                                    for(int in6=nb[6]*4; in6<nb[6]*4+C6.nbond; in6++){
                                        nb[7] = neighs[in6]; // 7th atom
                                        if ( nb[7] == nb[5] ) continue;

                                        // found 6-member ring
                                        if ( nb[7] == nb[1] ) {
                                            if ( explored[nb[1]] && explored[nb[2]] && explored[nb[3]] && explored[nb[4]] && explored[nb[5]] && explored[nb[6]] ) continue;
                                            nb[0] = nb[6];
                                            // count double bonds & heteroatoms
                                            if ( assignUFFtypes_checkaroma ( 6, &nb[0], &BOs[0], &BOs_int[0], tol ) ) {
                                                changed = true;
                                                for(int i=1; i<7; i++){explored[nb[i]] = true;}
                                                assignUFFtypes_setaroma ( 6, &nb[0], &neighs[0], &BOs[0], &set_atom[0], &set_bond[0], tol );
                                            }
                                        }

                                    } // in6
                                } // in5
                            } // in4
                        } // in3
                    } // in2
                } // in1
            } // ia

            if (!changed) break;
        }

    }

    bool assignUFFtypes_checkaroma ( int n, int* nb, double* BOs, int* BOs_int, double tol ){

        // check that all carbons in the ring are sp2
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            AtomConf& C = confs[A.iconf];
            if ( params->atypes[A.type].name[0] == 'C' && C.nbond != 3 ) return false;
        }
                
        int npi = 0;

        // count heteroatoms in the ring that can donate an electron pair        
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            AtomConf& C = confs[A.iconf];
            if ( ( params->atypes[A.type].name[0] == 'N' && C.nbond == 3 ) || 
                 ( params->atypes[A.type].name[0] == 'O' && C.nbond == 2 ) ) npi += 2;
        }

        int npi_save = npi;

        // count pi electrons in the ring considering double bonds in the ring only
        for(int i=1; i<n+1; i++){
            int ib1 = getBondByAtoms( nb[i], nb[i+1] );
            int ib2 = getBondByAtoms( nb[i-1], nb[i] );
            if ( BOs_int[ib1] == 2 || BOs_int[ib2] == 2 ) npi += 1;
        }

        if ( npi == 6 ) return true;

        npi = npi_save;

        // count pi electrons in the ring considering also fused rings
        for(int i=1; i<n+1; i++){
            int ib1 = getBondByAtoms( nb[i], nb[i+1] );
            int ib2 = getBondByAtoms( nb[i-1], nb[i] );
            if ( BOs_int[ib1] == 2 || abs(BOs[ib1]-1.5) < tol || BOs_int[ib2] == 2 || abs(BOs[ib2]-1.5) < tol ) npi += 1;
        }

        if ( npi == 6 ) return true;

        return false;

    }

    void assignUFFtypes_setaroma ( int n, int* nb, int* neighs, double* BOs, bool* set_atom, bool* set_bond, double tol ){

        // set atoms
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            if ( set_atom[nb[i]] ) {
                printf("WARNING SETAROMA: atom %c %i would be set to resonant but it has already a type of %i\n", params->atypes[A.type].name[0], nb[i]+1, A.type);
                continue;
            }
            if ( params->atypes[A.type].name[0] == 'C' ){ 
                A.type = params->getAtomType("C_R", true); 
            }else if ( params->atypes[A.type].name[0] == 'N' ){ 
                A.type = params->getAtomType("N_R", true);
            }else if ( params->atypes[A.type].name[0] == 'O' ){
                A.type = params->getAtomType("O_R", true);
            }else{
                printf("ERROR SIMPLERULE: atom %i type %i not recognized\n", nb[i]+1, A.type);
                printf("STOP\n");
                exit(0);
            }
            set_atom[nb[i]] = true;
        }
        // set ring bonds
        for(int i=1; i<n+1; i++){
            int ib = getBondByAtoms( nb[i], nb[i+1] );
            if ( set_bond[ib] && abs(BOs[ib]-1.5) > tol ) {
                Atom& A1 = atoms[nb[i]];
                Atom& A2 = atoms[nb[i+1]];
                printf("WARNING SETAROMA: bond %i between atoms %c %i and %c %i would be set to 1.5 but it has already a bond order of %g\n", ib+1, params->atypes[A1.type].name[0], nb[i]+1, params->atypes[A2.type].name[0], nb[i+1]+1, BOs[ib]);
                continue;
            }
            BOs[ib] = 1.5;
            set_bond[ib] = true;
        }

        // set other bonds
        for(int i=1; i<n+1; i++){
            Atom& A = atoms[nb[i]];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            for(int in=nb[i]*4; in<nb[i]*4+C.nbond; in++){
                bool ok = true;
                for(int j=1; j<n+1; j++){ if ( neighs[in] == nb[j] ) { ok = false; break; } }
                if ( !ok ) continue;
                Atom& Aj = atoms[neighs[in]];
                if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                AtomConf& Cj = confs[Aj.iconf];
                int ib = getBondByAtoms( nb[i], neighs[in] );
                if ( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 1 ) { // delocalized carbonyl
                    if ( set_atom[neighs[in]] && params->atypes[Aj.type].name[2] != 'R' ) {
                        printf("WARNING SETAROMA: carbonyl atom %c %i would be set to resonant but it has already a type of %i\n", params->atypes[Aj.type].name[0], neighs[in]+1, Aj.type);
                        continue;
                    }
                    if ( set_bond[ib] && abs(BOs[ib]-1.5) > tol ) {
                        printf("WARNING SETAROMA: bond %i between atoms %c %i and %c %i would be set to 1.5 but it has already a bond order of %g\n", ib+1, params->atypes[A.type].name[0], nb[i]+1, params->atypes[Aj.type].name[0], neighs[in]+1, BOs[ib]);
                        continue;
                    }
                    Aj.type = params->getAtomType("O_R", true);
                    set_atom[neighs[in]] = true;
                    BOs[ib] = 1.5;
                    set_bond[ib] = true;
                } //else {
                  //  BOs[ib] = 1.0;
                  //  set_bond[ib] = true; }
            }
        }

    }

    void assignUFFtypes_assignrest( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond ){

        for(int ia=0; ia<atoms.size(); ia++){
            if ( set_atom[ia] ) continue;
            Atom& A = atoms[ia];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if ( params->atypes[A.type].name[0] == 'C' ) {
                if ( C.nbond == 2 ) { 
                    A.type = params->getAtomType("C_1", true); 
                    set_atom[ia] = true; 
                } else if ( C.nbond == 3 ) { 
                    A.type = params->getAtomType("C_2", true); 
                    set_atom[ia] = true; 
                }
            } else if ( params->atypes[A.type].name[0] == 'N' ) {
                if ( C.nbond == 2 ) { 
                    A.type = params->getAtomType("N_2", true); 
                    set_atom[ia] = true; 
                } else if ( C.nbond == 3 ) { 
                    A.type = params->getAtomType("N_3", true); 
                    set_atom[ia] = true;
                    for(int in=ia*4; in<ia*4+3; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        BOs[ib] = 1.0;
                        BOs_int[ib] = 1;
                        set_bond[ib] = true;
                    }
                }
            } else if ( params->atypes[A.type].name[0] == 'O' ) {
                if ( C.nbond == 1 ) { 
                    A.type = params->getAtomType("O_2", true); 
                    set_atom[ia] = true; 
                    int ib = getBondByAtoms( ia, neighs[ia*4] );
                    BOs[ib] = 2.0;
                    BOs_int[ib] = 2;
                    set_bond[ib] = true;
                    Atom& Aj = atoms[neighs[ia*4]];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if ( params->atypes[A.type].name[0] == 'C' && !set_atom[neighs[ia*4]] && Cj.nbond == 3 ) {
                        Aj.type = params->getAtomType("C_2", true);
                        set_atom[neighs[ia*4]] = true;
                    }
                } else if ( C.nbond == 2 ) { 
                    A.type = params->getAtomType("O_3", true); 
                    set_atom[ia] = true; 
                    for(int in=ia*4; in<ia*4+2; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        BOs[ib] = 1.0;
                        BOs_int[ib] = 1;
                        set_bond[ib] = true;
                    }
                }
            }
        }

    }

    void assignUFFtypes_fixsaturation( int* neighs, double* BOs, int* BOs_int, bool* set_bond, double tol ){

        for (; true; ) {
            bool changed = false;

            for(int ia=0; ia<atoms.size(); ia++){
                Atom& A = atoms[ia];
                if ( params->atypes[A.type].name[0] == 'H' ) continue;
                AtomConf& C = confs[A.iconf];
                // if the atom has all bonds set, skip
                bool ok = true;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( BOs[ib] < 0.0 ) { ok = false; break; }
                }
                if ( ok ) continue;
                // compute number of already set bonds and corresponding atom valence
                int nset = 0;
                double val = 0.0;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( BOs[ib] > 0.0 ) { 
                        nset++; 
                        val+=(double)BOs_int[ib]; 
                    }
                }
                val = (double)params->atypes[A.type].valence - val;
                // only one bond to be set
                if ( nset == C.nbond-1 ) {
                    // check for weird valence
                    if ( abs( val - (double)rint(val) ) > tol || abs(val) < tol ) {
                        printf("ERROR FIXSATURATION: atom %c %i would have a valence of %g\n", params->atypes[A.type].name[0], ia+1, val);
                        printf("STOP\n");
                        exit(0);
                    }
                    changed = true;
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        if ( BOs[ib] < 0.0 ) {
                            BOs[ib] = val;
                            BOs_int[ib] = (int)rint(val);
                            set_bond[ib] = true;
                            break;
                        }
                    }
                // or multiple bonds to be set, but they are all single
                } else if ( (int)rint(val) == C.nbond-nset ) {
                        changed = true;
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ib = getBondByAtoms( ia, neighs[in] );
                            if ( BOs[ib] < 0.0 ) {
                                BOs[ib] = 1.0;
                                BOs_int[ib] = 1;
                                set_bond[ib] = true;
                            }
                        }
                }
            }
                
            if( !changed ) break;
        }

    }

    void assignUFFtypes_cumulene( int* neighs, double* BOs, int* BOs_int, bool* set_bond ){

        // ad hoc exception for getting rid of cumulenes in the backbone, i.e.
        /* ... C(ia3)                                           ... C(ia3)                          
               \\                              /                    \\                __            /     
                C(ia4) == C(ia1) == C(ia2) == C(ia3)	    ==>       C(ia4) -- C(ia1) == C(ia2) -- C(ia3) 
               /                              \\                     /                              \\           
                                                C(ia4) ...                                           C(ia4) ... */
        for(int ia1=0; ia1<atoms.size(); ia1++){
            Atom& A1 = atoms[ia1];
            if( params->atypes[A1.type].name[0] != 'C' ) continue; // looking for C
            AtomConf& C1 = confs[A1.iconf];
            if( C1.nbond != 2 ) continue; // looking for sp C
            for(int in1=ia1*4; in1<ia1*4+2; in1++){
                int ia2 = neighs[in1];
                Atom& A2 = atoms[ia2];
                if( params->atypes[A2.type].name[0] != 'C' ) continue; // looking for C
                AtomConf& C2 = confs[A2.iconf];
                if( C2.nbond != 2 ) continue; // looking for sp C
                int ib1 = getBondByAtoms( ia1, ia2 );
                int ia4;
                if(in1==ia1*4){ia4=neighs[ia1*4+1];}else{ia4=neighs[ia1*4];}
                Atom& A4 = atoms[ia4];
                if( params->atypes[A4.type].name[0] != 'C' ) continue; // looking for C
                AtomConf& C4 = confs[A4.iconf];
                if( C4.nbond != 3 ) continue; // looking for sp2 C
                int ib4 = getBondByAtoms( ia1, ia4 );
                for(int in2=ia2*4; in2<ia2*4+2; in2++){
                    int ia3 = neighs[in2];
                    if( ia3 == ia1 ) continue; // looking for the other neighbor of ia2
                    Atom& A3 = atoms[ia3];
                    if( params->atypes[A3.type].name[0] != 'C' ) continue; // looking for C
                    AtomConf& C3 = confs[A3.iconf];
                    if( C3.nbond != 3 ) continue; // looking for sp2 C
                    for(int in3=ia3*4; in3<ia3*4+3; in3++){
                        if( neighs[in3] != ia4 ) continue; // looking for the right ia3
                        int ib2 = getBondByAtoms( ia2, ia3 );
                        int ib3 = getBondByAtoms( ia3, ia4 );
                        BOs[ib1] = 3.0;
                        BOs[ib2] = 1.0;
                        BOs[ib3] = 2.0;
                        BOs[ib4] = 1.0;
                        BOs_int[ib1] = 3;
                        BOs_int[ib2] = 1;
                        BOs_int[ib3] = 2;
                        BOs_int[ib4] = 1;
                        set_bond[ib1] = true;
                        set_bond[ib2] = true;
                        set_bond[ib3] = true;
                        set_bond[ib4] = true;
                    }
                }
            }
        }
    }

    void assignUFFtypes_conjugation( int* neighs, double* BOs ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if( (params->atypes[A.type].name[0]=='N'&&C.nbond==3) || (params->atypes[A.type].name[0]=='O'&&C.nbond==2) ){
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[2]=='R' ) {
                        if ( params->atypes[A.type].name[0] == 'N' ) { A.type = params->getAtomType("N_R", true); } 
                        else if ( params->atypes[A.type].name[0] == 'O' ) { A.type = params->getAtomType("O_R", true); }
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.5;
                    }
                }
            }
        }

    }

    void assignUFFtypes_checks( int* neighs, double* BOs, int* BOs_int, bool* set_atom, bool* set_bond, double tol ){

        // sanity checks
        for(int ia=0; ia<atoms.size(); ia++){
            if(!set_atom[ia]){printf("ERROR CHECKS: atom %i is not set\n", ia+1);printf("STOP\n");exit(0);}
        }
        for(int ib=0; ib<bonds.size(); ib++){
            if(!set_bond[ib]){printf("ERROR CHECKS: bond %i is not set\n", ib+1);printf("STOP\n");exit(0);}
            if(BOs[ib]<0.0){printf("ERROR CHECKS: order for bond %i order is not set\n", ib+1);printf("STOP\n");exit(0);}
            if(BOs_int[ib]<0){printf("ERROR CHECKS: integer order for bond %i order is not set\n", ib+1);printf("STOP\n");exit(0);}
        }

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            // check that resonant atoms have at least one 1.5 bond
            if ( params->atypes[A.type].name[2]=='R' ) { 
                AtomConf& C = confs[A.iconf]; 
                bool found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( abs(BOs[ib]-1.5) < tol ) { found = true; break; }
                }
                if ( !found ) {
                    printf("ERROR CHECKS: atom %c %i is resonant but it has no 1.5 bonds\n", params->atypes[A.type].name[0], ia+1);
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        printf("  bond %i atoms %c %i %c %i bond order %g %i\n", ib+1, params->atypes[A.type].name[0], ia+1, params->atypes[atoms[neighs[in]].type].name[0], neighs[in]+1, BOs[ib], BOs_int[ib]);
                    }
                    printf("STOP\n");
                    exit(0);
                }
            // check that sp2 atoms have one localized double bond
            } else if ( params->atypes[A.type].name[2]=='2' ) { 
                AtomConf& C = confs[A.iconf]; 
                bool found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ib = getBondByAtoms( ia, neighs[in] );
                    if ( abs(BOs[ib]-2.0) < tol ) { 
                        if ( found ) {
                            printf("ERROR CHECKS: atom %c %i is sp2 but it has more than one double bond\n", params->atypes[A.type].name[0], ia+1);
                            for(int in=ia*4; in<ia*4+C.nbond; in++){
                                int ib = getBondByAtoms( ia, neighs[in] );
                                printf("  bond %i atoms %c %i %c %i bond order %g %i\n", ib+1, params->atypes[A.type].name[0], ia+1, params->atypes[atoms[neighs[in]].type].name[0], neighs[in]+1, BOs[ib], BOs_int[ib]);
                            }
                            printf("STOP\n");
                            exit(0);
                        }
                        found = true;
                    }
                }
                if ( !found ) {
                    found = false;
                    for(int in=ia*4; in<ia*4+C.nbond; in++){
                        int ib = getBondByAtoms( ia, neighs[in] );
                        if ( abs(BOs[ib]-1.5) < tol ) { found = true; break; }
                    }
                    if ( found ) {
                        printf("WARNING CHECKS: atom %c %i is sp2 and it has no double bonds, only 1.5 bonds\n", params->atypes[A.type].name[0], ia+1);
                    } else {
                        printf("ERROR CHECKS: atom %c %i is sp2 but it has no double bonds\n", params->atypes[A.type].name[0], ia+1);
                        for(int in=ia*4; in<ia*4+C.nbond; in++){
                            int ib = getBondByAtoms( ia, neighs[in] );
                            printf("  bond %i atoms %c %i %c %i bond order %g %i\n", ib+1, params->atypes[A.type].name[0], ia+1, params->atypes[atoms[neighs[in]].type].name[0], neighs[in]+1, BOs[ib], BOs_int[ib]);
                        }
                        printf("STOP\n");
                        exit(0);
                    }
                }
            }
        }

        // check that 1.5 bonds must be either 1 or 2 in the limit resonance structure
        for(int ib=0; ib<bonds.size(); ib++){
            if ( abs(BOs[ib]-1.5) < tol ) {
                if ( BOs_int[ib] != 1 && BOs_int[ib] != 2 ) {
                    printf("ERROR CHECKS: bond %i is 1.5 but in the limit resonance structure is neither single nor double\n", ib+1);
                    printf("STOP\n");
                    exit(0);
                }
            }
        }
    }

    void assignUFFtypes_amide( int* neighs, double* BOs ){

        for(int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            if ( params->atypes[A.type].name[0] == 'H' ) continue;
            AtomConf& C = confs[A.iconf];
            if( params->atypes[A.type].name[0] == 'C' && C.nbond == 3 ) {
                // look for a carbonyl oxygen
                bool found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 1 ) { found = true; break; }
                }
                if ( !found ) continue;
                // look for an amino nitrogen
                found = false;
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond == 3 ) { found = true; break; }
                }
                if ( !found ) continue;
                A.type = params->getAtomType("C_R", true);
                for(int in=ia*4; in<ia*4+C.nbond; in++){
                    int ja = neighs[in];
                    Atom& Aj = atoms[ja];
                    if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
                    AtomConf& Cj = confs[Aj.iconf];
                    if( params->atypes[Aj.type].name[0] == 'O' && Cj.nbond == 1 ) {
                        Aj.type = params->getAtomType("O_R", true);
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.5;
                    } else if( params->atypes[Aj.type].name[0] == 'N' && Cj.nbond == 3 ) {
                        Aj.type = params->getAtomType("N_R", true);
                        int ib = getBondByAtoms( ia, ja );
                        BOs[ib] = 1.41;
                    }
                }
            }
        }

    }

    // UFF atom-type assignement
    void assignUFFtypes( int* neighs=0, bool bCumulene=false, bool bDeallocNeighs=true, bool b141=true, bool bSimple=true, bool bConj=false ){ 

        // init working arrays
        double tol = 0.05;
        bool   set_atom[atoms.size()];
        bool   set_bond[bonds.size()];
        std::vector<double> BOs(bonds.size());
        std::vector<int>    BOs_int(bonds.size());
        for(int ia=0; ia<atoms.size(); ia++){set_atom[ia]=false;}
        for(int ib=0; ib<bonds.size(); ib++){BOs[ib]=-1.0; BOs_int[ib]=-1; set_bond[ib]=false;}

        bDeallocNeighs &= (neighs==0);
        makeNeighs ( neighs, 4 );

        // assign bond orders and atom types for trivial cases
        // e.g. sp3 carbons, hydrogens, etc...
        assignUFFtypes_trivial( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );

        // assign bond orders and atom types for nitro groups
        assignUFFtypes_nitro( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );

        // find a (hopefully) valid limit resonance structure
        assignUFFtypes_treewalk( neighs, &BOs_int[0] );

        if ( bSimple ) {
            // assign resonant atoms according to the "simple" rule: atoms with two sp2 neighbors (or heteroatoms) are resonant
            assignUFFtypes_simplerule( tol, neighs, &BOs[0], &set_atom[0], &set_bond[0] );
        } else {
            // assign resonant atoms according to the Huckel rule for aromaticity
            assignUFFtypes_findrings( tol, neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );
        }

        // assign the rest of atom types
        assignUFFtypes_assignrest( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0] );

        // try to assign the rest of double bonds
        assignUFFtypes_fixsaturation( neighs, &BOs[0], &BOs_int[0], &set_bond[0], tol );

        // exception to avoid cumulenes
        if( bCumulene ) { assignUFFtypes_cumulene( neighs, &BOs[0], &BOs_int[0], &set_bond[0] ); }
  
        // manually change sp3 nitrogen and oxygen to "resonant" when they are bonded to an sp2 atom (conjugation)
        if( bConj ){ assignUFFtypes_conjugation( neighs, &BOs[0] ); }

        // some check before assigning parameters
        assignUFFtypes_checks( neighs, &BOs[0], &BOs_int[0], &set_atom[0], &set_bond[0], tol );

        // exception for amide groups
        assignUFFtypes_amide( neighs, &BOs[0] );

//DEBUG
for(int ia=0; ia<atoms.size(); ia++){
    Atom& A = atoms[ia];
    // printf("MMFFBuilder::assignUFFtypes(%i) %s\n", ia, params->atypes[A.type].name );
}
//DEBUG
//exit(0);

        // store bond orders
        for(int ib=0; ib<bonds.size(); ib++){
            Bond& B = bonds[ib];
            B.order = BOs[ib];
            //if(verbosity>0)printf( "bondOrder[%i](%i-%i)= %g\n", ib+1, B.atoms.i+1, B.atoms.j+1, B.order );            
        }

        // deallocate neighbor array
        if(bDeallocNeighs)delete [] neighs;

    }

    double assignUFFparams_calcrij ( int ib ){

        Bond& B = bonds[ib];   // reference
        //Bond* B_ = &bonds[ib]; // pointer

        AtomType& Ti = params->atypes[atoms[B.atoms.i].type];
        AtomType& Tj = params->atypes[atoms[B.atoms.j].type];
        const ElementType& Ei = *params->elementOfAtomType(atoms[B.atoms.i].type);
        const ElementType& Ej = *params->elementOfAtomType(atoms[B.atoms.j].type);
        double rBO = -0.1332 * ( Ti.Ruff + Tj.Ruff ) * log( B.order );
        double rEN = Ti.Ruff * Tj.Ruff * sq( sqrt(-Ei.Eaff) - sqrt(-Ej.Eaff) ) / ( -Ei.Eaff*Ti.Ruff -Ej.Eaff*Tj.Ruff );
        double rIJ = Ti.Ruff + Tj.Ruff + rBO - rEN;
        return rIJ;

    }

    void assignUFFparams_vdws( ){

        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            params->assignRE( A.type, A.REQ, true );
        }

    }

    void assignUFFparams_bonds( ){

        for ( int ib=0; ib<bonds.size(); ib++ ) {
            Bond& B = bonds[ib];
            const ElementType& Ei = *params->elementOfAtomType(atoms[B.atoms.i].type);
            const ElementType& Ej = *params->elementOfAtomType(atoms[B.atoms.j].type);
            B.l0 = assignUFFparams_calcrij(ib);
            B.k = 0.5 * 28.7989689090648 * Ei.Quff * Ej.Quff / ( B.l0*sq(B.l0) );
        }

    }

    void assignUFFparams_angles( int* neighs ){

        angles.clear();
        for( int j=0; j<atoms.size(); j++){
            Atom& Aj = atoms[j];
            AtomType& Tj = params->atypes[Aj.type];
            double ct = cos(Tj.Ass*deg2rad);
            double st2 = sq(sin(Tj.Ass*deg2rad));
            if ( params->atypes[Aj.type].name[0] == 'H' ) continue;
            AtomConf& Cj = confs[Aj.iconf];
            for(int in1=j*4; in1<j*4+Cj.nbond-1; in1++){
                int i = neighs[in1];
                const ElementType& Ei = *params->elementOfAtomType(atoms[i].type);
                for(int in2=in1+1; in2<j*4+Cj.nbond; in2++){
                    int k = neighs[in2];
                    const ElementType& Ek = *params->elementOfAtomType(atoms[k].type);
                    Angle a;
                    a.atoms.x = i;
                    a.atoms.y = j;
                    a.atoms.z = k;
                    a.bonds.x = getBondByAtoms( i, j );
                    a.bonds.y = getBondByAtoms( j, k );
                    double rij = assignUFFparams_calcrij(a.bonds.x);
                    double rjk = assignUFFparams_calcrij(a.bonds.y);
                    double rik = sqrt( rij*rij + rjk*rjk - 2.0*rij*rjk*ct );
                    double kappa = 28.7989689090648 * Ei.Quff * Ek.Quff / (sq(rik)*sq(rik)*rik) * ( 3.0*rij*rjk*st2 - sq(rik)*ct );
                    if ( params->atypes[Aj.type].name[2]=='1' || params->atypes[Aj.type].name[2]=='2' || 
                    params->atypes[Aj.type].name[2]=='R' ) { // cosine/periodic
                        if ( params->atypes[Aj.type].name[2]=='1' ) {
                            //a.k = 0.5 * kappa;
                            a.k = kappa;
                            a.C0 = 1.0;
                            a.C1 = 1.0;
                            a.C2 = 0.0;
                            a.C3 = 0.0;
                        } else if ( params->atypes[Aj.type].name[2]=='2' || params->atypes[Aj.type].name[2]=='R' ) {
                            a.k = kappa / 9.0;
                            a.C0 = 1.0;
                            a.C1 = 0.0;
                            a.C2 = 0.0;
                            a.C3 = -1.0;
                        } 
                    } else if ( params->atypes[Aj.type].name[2]=='3' ) { // fourier
                        a.k = kappa;
                        a.C2 = 1.0 / ( 4.0 * st2 );
                        a.C1 = -4.0 * a.C2 * ct;
                        a.C0 = a.C2 * ( 2.0*sq(ct) + 1.0 );
                        a.C3 = 0.0;
                    }
                    angles.push_back(a);
                }
            }
        }

    }

    void assignUFFparams_dihedrals( int* neighs ){

        dihedrals.clear();
        for( int i1=0; i1<atoms.size(); i1++){
            Atom& A1 = atoms[i1];
            int n1;
            if ( params->atypes[A1.type].name[0] == 'H' ){ n1 = 1; } 
            else { AtomConf& C1 = confs[A1.iconf]; n1 = C1.nbond; }
            for(int in1=i1*4; in1<i1*4+n1; in1++){
                int i2 = neighs[in1];
                Atom& A2 = atoms[i2];
                if ( params->atypes[A2.type].name[0] == 'H' ) continue;
                if ( params->atypes[A2.type].name[2] == '1' ) continue; // Torsional potentials for central bonds involving sp-hybridized centers X-1 were assigned a value of zero
                AtomConf& C2 = confs[A2.iconf];
                for(int in2=i2*4; in2<i2*4+C2.nbond; in2++){
                    int i3 = neighs[in2];
                    if ( i3 != i1 ) {
                        Atom& A3 = atoms[i3];
                        if ( params->atypes[A3.type].name[0] == 'H' ) continue;
                        if ( params->atypes[A3.type].name[2] == '1' ) continue; // Torsional potentials for central bonds involving sp-hybridized centers X-1 were assigned a value of zero
                        AtomConf& C3 = confs[A3.iconf];
                        for(int in3=i3*4; in3<i3*4+C3.nbond; in3++){
                            int i4 = neighs[in3];
                            if ( i4 != i2 && i4 > i1 ) { // avoid 3-membered rings and double-counting
                                Atom& A4 = atoms[i4];
                                Dihedral d;
                                d.atoms.x = i1;
                                d.atoms.y = i2;
                                d.atoms.z = i3;
                                d.atoms.w = i4;
                                d.bonds.x = getBondByAtoms( i1, i2 );
                                d.bonds.y = getBondByAtoms( i2, i3 );
                                d.bonds.z = getBondByAtoms( i3, i4 );
                                const ElementType& E2 = *params->elementOfAtomType(atoms[i2].type);
                                const ElementType& E3 = *params->elementOfAtomType(atoms[i3].type);
                                // specific general case (a): * - sp3 - sp3 - *
                                if ( params->atypes[A2.type].name[2] == '3' && params->atypes[A3.type].name[2] == '3' ){
                                    d.k = sqrt( E2.Vuff * E3.Vuff );
                                    d.d = 1;
                                    d.n = 3;
                                    // special case of * - group 16 sp3 - group 16 sp3 - *
                                    if ( (params->atypes[A2.type].name[0]=='O'||params->atypes[A2.type].name[0]=='S') &&
                                        (params->atypes[A3.type].name[0]=='O'||params->atypes[A3.type].name[0]=='S') ){
                                            d.k = 4.1840/60.2214076/1.602176634; // 1 kcal/mol to eV
                                            d.n = 2;
                                            if ( params->atypes[A2.type].name[0] == 'O' ){ d.k *= 2.0; }
                                            else { d.k *= 6.8; }
                                            if ( params->atypes[A3.type].name[0] == 'O' ){ d.k *= 2.0; }
                                            else { d.k *= 6.8; }
                                            d.k = sqrt( d.k );
                                    }
                                } 
                                // specific general case (b): * - sp3 - sp2 - *
                                else if ( (params->atypes[A2.type].name[2]=='3'&&(params->atypes[A3.type].name[2]=='2'||params->atypes[A3.type].name[2]=='R')) ||
                                ((params->atypes[A2.type].name[2]=='2'||params->atypes[A2.type].name[2]=='R')&&params->atypes[A3.type].name[2]=='3') ){
                                    d.k = 4.1840/60.2214076/1.602176634; // 1 kcal/mol to eV
                                    d.d = -1;
                                    d.n = 6;
                                    // special case of * - group 16 sp3 - sp2 - *
                                    if ( (params->atypes[A2.type].name[2]=='3'&&(params->atypes[A2.type].name[0]=='O'||params->atypes[A2.type].name[0]=='S')) ||
                                        (params->atypes[A3.type].name[2]=='3'&&(params->atypes[A3.type].name[0]=='O'||params->atypes[A3.type].name[0]=='S')) ){
                                        int ib = getBondByAtoms( i2, i3 );
                                        Bond& B = bonds[ib];
                                        d.k = 5.0 * sqrt( E2.Uuff * E3.Uuff ) * ( 1.0 + 4.18 * log( B.order ) );
                                        d.d = 1;
                                        d.n = 2;
                                    }
                                    // special case of * - sp3 - sp2 bounded to another sp2 atom
                                    if ( params->atypes[A2.type].name[2]=='3' ){
                                        bool found = false;
                                        for(int in=i3*4; in<i3*4+C3.nbond; in++){
                                            Atom& A = atoms[neighs[in]];
                                            if ( params->atypes[A.type].name[2]=='2' || params->atypes[A.type].name[2]=='R' ) {
                                                found = true;
                                                break;
                                            }
                                        }
                                        if ( found ) {
                                            d.k = 2.0 * 4.1840/60.2214076/1.602176634; // 2 kcal/mol to eV
                                            d.d = 1;
                                            d.n = 3;
                                        }
                                    } else if ( params->atypes[A3.type].name[2]=='3' ){
                                        bool found = false;
                                        for(int in=i2*4; in<i2*4+C2.nbond; in++){
                                            Atom& A = atoms[neighs[in]];
                                            if ( params->atypes[A.type].name[2]=='2' || params->atypes[A.type].name[2]=='R' ) {
                                                found = true;
                                                break;
                                            }
                                        }
                                        if ( found ) {
                                            d.k = 2.0 * 4.1840/60.2214076/1.602176634; // 2 kcal/mol to eV
                                            d.d = 1;
                                            d.n = 3;
                                        }
                                    }
                                }
                                // specific general case (c): * - sp2 - sp2 - *
                                else if ( (params->atypes[A2.type].name[2]=='2'||params->atypes[A2.type].name[2]=='R') &&
                                        (params->atypes[A3.type].name[2]=='2'||params->atypes[A3.type].name[2]=='R') ){
                                    int ib = getBondByAtoms( i2, i3 );
                                    Bond& B = bonds[ib];
                                    d.k = 5.0 * sqrt( E2.Uuff * E3.Uuff ) * ( 1.0 + 4.18 * log( B.order ) );
                                    d.d = -1;
                                    d.n = 2;
                                } else {
                                    printf("ERROR: assignUFFparams_dihedrals: torsion case not found for atoms %s %s %s %s\n", params->atypes[A1.type].name, params->atypes[A2.type].name, params->atypes[A3.type].name, params->atypes[A4.type].name);
                                    printf("STOP\n");
                                    exit(0);
                                }
                                d.k = 0.5 * d.k / ( (double)(C2.nbond-1) * (double)(C3.nbond-1) );
                                dihedrals.push_back(d);
                            }
                        }
                    }
                }
            }
        }

    }       

    void assingUFFparams_assigninversion( int i1, int i2, int i3, int i4 ){

        Atom& A1 = atoms[i1];
        Atom& A2 = atoms[i2];
        Atom& A3 = atoms[i3];
        Atom& A4 = atoms[i4];

        Inversion i;
        i.atoms.x = i1;
        i.atoms.y = i2;
        i.atoms.z = i3;
        i.atoms.w = i4;
        i.bonds.x = getBondByAtoms( i1, i2 );
        i.bonds.y = getBondByAtoms( i1, i3 );
        i.bonds.z = getBondByAtoms( i1, i4 );
        // sp2 carbon
        if ( params->atypes[A1.type].name[0] == 'C' && (params->atypes[A1.type].name[2]=='2'||params->atypes[A1.type].name[2]=='R') ){
            // carbonyl
            if ( ( params->atypes[A2.type].name[0] == 'O' && params->atypes[A2.type].name[2] == '2' ) ||
                 ( params->atypes[A3.type].name[0] == 'O' && params->atypes[A3.type].name[2] == '2' ) ||
                 ( params->atypes[A4.type].name[0] == 'O' && params->atypes[A4.type].name[2] == '2' ) ){
                i.k = 50.0 * 4.1840/60.2214076/1.602176634; } // 50 kcal/mol to eV
            else { i.k = 6.0 * 4.1840/60.2214076/1.602176634; } // 6 kcal/mol to eV
            i.C0 = 1.0;
            i.C1 = -1.0;
            i.C2 = 0.0;
        }
        // sp2 nitrogen
        else if ( params->atypes[A1.type].name[0] == 'N' && (params->atypes[A1.type].name[2]=='2'||params->atypes[A1.type].name[2]=='R') ){
            i.k = 6.0 * 4.1840/60.2214076/1.602176634; // 6 kcal/mol to eV
            i.C0 = 1.0;
            i.C1 = -1.0;
            i.C2 = 0.0;
        }
        // sp3 nitrogen
        else if ( params->atypes[A1.type].name[0] == 'N' && params->atypes[A1.type].name[2]=='3' ){
            i.k = 0.0;
            i.C0 = 0.0;
            i.C1 = 0.0;
            i.C2 = 0.0;
        }
        // sp3 group 15
        else if ( params->atypes[A1.type].name[0] == 'P' || params->atypes[A1.type].name[0] == '3' ){
            double omega0 = 84.4339 * deg2rad;
            i.C0 = 4.0 * sq(cos(omega0)) - cos(2.0*omega0);
            i.C1 = -4.0 * cos(omega0);
            i.C2 = 1.0;
            i.k = 22.0 * 4.1840/60.2214076/1.602176634 / ( i.C0 + i.C1 + i.C2 ); // 22 kcal/mol to eV
        } else {
            printf("ERROR: assignUFFparams_inversions: improper case not found for atoms %s %s %s %s\n", params->atypes[A1.type].name, params->atypes[A2.type].name, params->atypes[A3.type].name, params->atypes[A4.type].name);
            printf("STOP\n");
            exit(0);
        }
        i.k = i.k / 3.0;
        inversions.push_back(i);

    }

    void assignUFFparams_inversions( int* neighs ){

        inversions.clear();
        for( int i1=0; i1<atoms.size(); i1++){
            Atom& A1 = atoms[i1];
            if ( params->atypes[A1.type].name[0] == 'H' ) continue;
            AtomConf& C1 = confs[A1.iconf]; 
            if ( C1.nbond != 3 ) continue;
            int i2 = neighs[i1*4];
            int i3 = neighs[i1*4+1];
            int i4 = neighs[i1*4+2];
            assingUFFparams_assigninversion( i1, i2, i3, i4 );
            assingUFFparams_assigninversion( i1, i4, i2, i3 );
            assingUFFparams_assigninversion( i1, i3, i4, i2 );
        }

    }

    void assignUFFparams( int* neighs=0, bool bDeallocNeighs=true ){

        bDeallocNeighs &= (neighs==0);
        makeNeighs ( neighs, 4 );

        // assign vdw parameters
        assignUFFparams_vdws();   

        // assign bond parameters
        assignUFFparams_bonds();

        // populate angle array & assign angle parameters
        assignUFFparams_angles( neighs );
   
        // assign torsion parameters
        assignUFFparams_dihedrals( neighs );
    
        // assign inversion parameters
        assignUFFparams_inversions( neighs );

        // deallocate neighbor array
        if(bDeallocNeighs)delete [] neighs;

        // write data
        //assignUFFparams_writedata();
          
    }
  
    void assignUFFparams_writedata( ){

        double tokcal = 60.2214076*1.602176634/4.1840;
        //double toeV = 4.1840/60.2214076/1.602176634;
        FILE* f=fopen( "mol.data.firecore", "w");
        fprintf( f, "\n" );
        fprintf( f, "\n" );
        fprintf( f, "%i atoms\n", atoms.size() );
        fprintf( f, "%i atom types\n", atoms.size() );
        fprintf( f, "%i bonds\n", bonds.size() );
        fprintf( f, "%i bond types\n", bonds.size() );
        fprintf( f, "%i angles\n", angles.size() );
        fprintf( f, "%i angle types\n", angles.size() );
        fprintf( f, "%i dihedrals\n", dihedrals.size() );
        fprintf( f, "%i dihedral types\n", dihedrals.size() );
        fprintf( f, "%i impropers\n", inversions.size() );
        fprintf( f, "%i improper types\n", inversions.size() );
        fprintf( f, "\n" );
        fprintf( f, "0.0 64.0 xlo xhi\n" );
        fprintf( f, "0.0 64.0 ylo yhi\n" );
        fprintf( f, "0.0 64.0 zlo zhi\n" );
        fprintf( f, "0.0 0.0 0.0 xy xz yz\n" );
        fprintf( f, "\n" );
        fprintf( f, "Masses\n" );
        fprintf( f, "\n" );
        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            fprintf( f, "%i 1.0 # %s\n", ia+1, params->atypes[A.type].name );
        }
        fprintf( f, "\n" );
        fprintf( f, "Pair Coeffs\n" );
        fprintf( f, "\n" );
        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            fprintf( f, "%i %f %f # %s %s\n", ia+1, sq(A.REQ.y)*tokcal, A.REQ.x*2.0*0.890898718140339, params->atypes[A.type].name, params->atypes[A.type].name );
        }
        if ( bonds.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Bond Coeffs\n" );
            fprintf( f, "\n" );
            for( int ib=0; ib<bonds.size(); ib++){
                Bond& B = bonds[ib];
                fprintf( f, "%i %f %f # %s %s\n", ib+1, 0.5*B.k*tokcal, B.l0, params->atypes[atoms[B.atoms.i].type].name, params->atypes[atoms[B.atoms.j].type].name );
            }
        }
        if ( angles.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Angle Coeffs\n" );
            fprintf( f, "\n" );
            bool bFourier = false;
            bool bCosine = false;
            for( int ia=0; ia<angles.size(); ia++){
                Angle& A = angles[ia];
                int i = A.atoms.j;
                if ( params->atypes[atoms[i].type].name[2] == '1' || params->atypes[atoms[i].type].name[2] == '2' || params->atypes[atoms[i].type].name[2] == 'R' ) { bCosine = true; }
                if ( params->atypes[atoms[i].type].name[2] == '3' ) { bFourier = true; }
            }
            for( int ia=0; ia<angles.size(); ia++){
                Angle& A = angles[ia];
                int i = A.atoms.j;
                if ( bFourier && bCosine ) { 
                    if ( params->atypes[atoms[i].type].name[2] == '1' ) { // cosine/periodic
                        fprintf( f, "%i cosine/periodic %f 1 1 # %s %s %s\n", ia+1, A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '2' || params->atypes[atoms[i].type].name[2] == 'R' ) { // cosine/periodic
                        fprintf( f, "%i cosine/periodic %f -1 3 # %s %s %s\n", ia+1, 9.0*A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '3' ) { // fourier
                        fprintf( f, "%i fourier %f %f %f %f # %s %s %s\n", ia+1, A.k*tokcal, A.C0, A.C1, A.C2, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    }
                } else {
                    if ( params->atypes[atoms[i].type].name[2] == '1' ) { // cosine/periodic
                        fprintf( f, "%i %f 1 1 # %s %s %s\n", ia+1, A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '2' || params->atypes[atoms[i].type].name[2] == 'R' ) { // cosine/periodic
                        fprintf( f, "%i %f -1 3 # %s %s %s\n", ia+1, 9.0*A.k*tokcal, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    } else if ( params->atypes[atoms[i].type].name[2] == '3' ) { // fourier
                        fprintf( f, "%i %f %f %f %f # %s %s %s\n", ia+1, A.k*tokcal, A.C0, A.C1, A.C2, params->atypes[atoms[A.atoms.x].type].name, params->atypes[atoms[A.atoms.y].type].name, params->atypes[atoms[A.atoms.z].type].name );
                    }
                }
            }
        }
        if ( dihedrals.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Dihedral Coeffs\n" );
            fprintf( f, "\n" );
            for( int id=0; id<dihedrals.size(); id++){
                Dihedral& D = dihedrals[id];
                fprintf( f, "%i %f %i %i # %s %s %s %s\n", id+1, D.k*tokcal, D.d, D.n, params->atypes[atoms[D.atoms.x].type].name, params->atypes[atoms[D.atoms.y].type].name, params->atypes[atoms[D.atoms.z].type].name, params->atypes[atoms[D.atoms.w].type].name );
            }
        }
        if ( inversions.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Improper Coeffs\n" );
            fprintf( f, "\n" );
            for( int ii=0; ii<inversions.size(); ii++){
                Inversion& I = inversions[ii];
                fprintf( f, "%i %f %f %f %f 0 # %s %s %s %s\n", ii+1, I.k*tokcal, I.C0, I.C1, I.C2, params->atypes[atoms[I.atoms.x].type].name, params->atypes[atoms[I.atoms.y].type].name, params->atypes[atoms[I.atoms.z].type].name, params->atypes[atoms[I.atoms.w].type].name );
            }
        }
        fprintf( f, "\n" );
        fprintf( f, "Atoms\n" );
        fprintf( f, "\n" );
        for( int ia=0; ia<atoms.size(); ia++){
            Atom& A = atoms[ia];
            fprintf( f, "%i 1 %i 0.0 %g %g %g\n", A.id+1, ia+1, A.pos.x, A.pos.y, A.pos.z );
        }
        if ( bonds.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Bonds\n" );
            fprintf( f, "\n" );
            for( int ib=0; ib<bonds.size(); ib++){
                Bond& B = bonds[ib];
                fprintf( f, "%i %i %i %i\n", ib+1, ib+1, atoms[B.atoms.i].id+1, atoms[B.atoms.j].id+1 );
            }
        }
        if ( angles.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Angles\n" );
            fprintf( f, "\n" );
            for( int ia=0; ia<angles.size(); ia++){
                Angle& A = angles[ia];
                fprintf( f, "%i %i %i %i %i\n", ia+1, ia+1, atoms[A.atoms.i].id+1, atoms[A.atoms.j].id+1, atoms[A.atoms.k].id+1 );
            }
        }
        if ( dihedrals.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Dihedrals\n" );
            fprintf( f, "\n" );
            for( int id=0; id<dihedrals.size(); id++){
                Dihedral& D = dihedrals[id];
                fprintf( f, "%i %i %i %i %i %i\n", id+1, id+1, atoms[D.atoms.x].id+1, atoms[D.atoms.y].id+1, atoms[D.atoms.z].id+1, atoms[D.atoms.w].id+1 );
            }
        }
        if ( inversions.size() > 0 ){
            fprintf( f, "\n" );
            fprintf( f, "Impropers\n" );
            fprintf( f, "\n" );
            for( int ii=0; ii<inversions.size(); ii++){
                Inversion& I = inversions[ii];
                fprintf( f, "%i %i %i %i %i %i\n", ii+1, ii+1, atoms[I.atoms.x].id+1, atoms[I.atoms.y].id+1, atoms[I.atoms.z].id+1, atoms[I.atoms.w].id+1 );
            }
        }
        fclose(f);

    }

}; // MMFFBuilder


} // namespace MMFF

#endif // MMFFBuilder_h
        