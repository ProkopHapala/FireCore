#ifndef MMFFBuilderBase_h
#define MMFFBuilderBase_h
/// @file MMFFBuilder.h   @brief Classes for building and editing molecular topology and building instances of force-fields for paticular molecule or system 
/// @ingroup Classical_Molecular_Mechanics

#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <memory>

#include  "globals.h"

#include "macroUtils.h"
//#include "testUtils.h"

//#include "Molecule.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
//#include "raytrace.h"

//#include "molecular_utils.h"

//#include "Molecule.h"
//#include "MMFF.h"
//#include "MMFFmini.h"
#include "MMFFparams.h"


//#include "MolecularGraph.h"
//#include "LimitedGraph.h"


//#include "Kekule.h"


// =============== Structs for Atoms, Bonds etc...

/// @namespace MM       @brief Molecular-Mechanics namespace
namespace MM{

//static const double const_eVA2_Nm = 16.02176634;

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

    void print()const{ printf( "Atom{id %i t %i c %i f %i REQ(%g,%g,%g,%g) pos(%g,%g,%g)}", id, type, iconf, frag, REQ.x, REQ.y, REQ.z,REQ.w, pos.x,pos.y,pos.z ); }

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

    inline int sum_bonds()const{ return nbond+npi+ne+nH; }

    // compute the number of pi bonds
    inline int fillIn_npi_sp3(){ npi=4-nbond-ne-nH; npi=(npi>2)?2:npi; n=4; return npi; };

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
        n=sum_bonds();
        if(n>=N_NEIGH_MAX){ return false; }
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
// ========= Builder
// ============================
/// @brief Comprehensive builder for for building and editing molecular topology and building instances of force-fields for paticular molecule or system.
/// @details It implements many algorithms for automatic finding of bonding topology, assignment of atomic types, adding/removing explicit free electron paris and capping atoms, substitution, segmentation to groups, etc.
class BuilderBase{  public:
    //static int iDebug = 0;
    std::vector<Atom>       atoms;
    std::vector<Bond>       bonds;
    std::vector<Angle>      angles;
    std::vector<Dihedral>   dihedrals;
    std::vector<Inversion>  inversions;
    std::vector<AtomConf>   confs;
    std::vector<Fragment>   frags;
    //std::vector<int>  atom_neighs;
    // --- other
    bool  bPBC = false;
    Mat3d lvec = Mat3dIdentity;
    MMFFparams* params = 0;  // Each function which needs this can take it as parameter
    std::unordered_set<int> capping_types;
    // --------- Aux Params
    // --- Types
    int  itypHcap      = -1;
    int  itypEpair     = -1;
    int  itypSigmaHole = -1;  // optional dedicated sigma-hole dummy type (e.g., E_h)
    int  itype_min     =  1; // type 0 is * (not valid type)
    int  ignoreType    = -1;
    // --- Default params
    Quat4d defaultREQ   {  1.5, 0.0, 0.0, 0.0    };
    Bond   defaultBond  { -1, {-1,-1}, 1.5, 10.0 };
    Angle  defaultAngle { -1, {-1,-1}, M_PI, 0.5 };
    Bond   bondBrush = defaultBond;
    // --- Caps
    Atom   capAtom      = Atom{ (int)NeighType::H,     -1,-1, {0,0,0}, Atom::HcapREQ };
    Atom   capAtomEpair = Atom{ (int)NeighType::epair, -1,-1, {0,0,0}, {0,0,0} };
    Atom   capAtomPi    = Atom{ (int)NeighType::pi,    -1,-1, {0,0,0}, {0,0,0} };
    Bond   capBond      = Bond{ -1,  {-1,-1},  1.07, 100/const_eVA2_Nm };
    Vec3d  capUp        = Vec3d{0.0,0.0,1.0};
    // --- bools
    bool   bondByRvdW  = false;
    bool   bDummyPi    = false;
    bool   bDummyEpair = false;
    bool   bAutoTypes  = true;
    bool   bAddCaps    = true;
    bool   bAddECaps   = false;

    void cloneFrom( const BuilderBase& o ){
        atoms         = o.atoms;
        bonds         = o.bonds;
        angles        = o.angles;
        dihedrals     = o.dihedrals;
        inversions    = o.inversions;
        confs         = o.confs;
        frags         = o.frags;
        // --- other
        bPBC          = o.bPBC;
        lvec          = o.lvec;
        params        = o.params;
        capping_types = o.capping_types;
        // --------- Aux Params
        // --- Types
        itypHcap      = o.itypHcap;
        itypEpair     = o.itypEpair;
        itypSigmaHole = o.itypSigmaHole;
        itype_min     = o.itype_min;
        ignoreType    = o.ignoreType;
        // --- Default params
        defaultREQ    = o.defaultREQ;
        defaultBond   = o.defaultBond;
        defaultAngle  = o.defaultAngle;
        bondBrush     = o.bondBrush;
        // --- Caps
        capAtom       = o.capAtom;
        capAtomEpair  = o.capAtomEpair;
        capAtomPi     = o.capAtomPi;
        capBond       = o.capBond;
        capUp         = o.capUp;
        // --- bools
        bondByRvdW    = o.bondByRvdW;
        bDummyPi      = o.bDummyPi;
        bDummyEpair   = o.bDummyEpair;
        bAutoTypes    = o.bAutoTypes;
        bAddCaps      = o.bAddCaps;
        bAddECaps     = o.bAddECaps;
    }
    
    // ============== Functions

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
    }

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
                int it  = atoms[ia].type;
                printf("ia %i t %i %-8s ic %i", ia, it, params->atypes[it].name ); confs[ic].print();
                printf("\nbond(%i-%i) %s-%s \n", bonds[ib].atoms.a, bonds[ib].atoms.b, params->atypes[it1].name, params->atypes[it2].name );
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

    void setConfs( int npi, int ne, int imin=0, int imax=-1 ){
        if(imax<0){ imax=atoms.size(); }
        for(int i=imin;i<imax;i++){
            makeSPConf( i, npi, ne );
        }
    }

    inline double getTypeRadius( int type ){
        if(bondByRvdW){ return params->atypes[type].RvdW; }
        else          { return params->atypes[type].Ruff; }
    }

    void touchingAtoms( int i0, int imax, const Vec3d& p, double R0, double Rfac, std::vector<int>& found ){
        for(int i=i0; i<imax; i++){  // for pbc we need all atom pairs
            const Atom& A = atoms[i];
            Vec3d dp = A.pos - p; // pbc here
            double Rj = getTypeRadius(A.type);
            //double Rj = (R0 + A.REQ.x)*Rfac;    // Using RvdW
            double R = (Rj+R0)*Rfac;
            //printf( "touchingAtoms[%i] r=%16.6f <? R(%12.6f) | Ri+j=%16.6f Ri=%16.6f Rj=%16.6f Rfac=%g \n", i, dp.norm(), R, R0+Rj, R0, Rj, Rfac );
            if(  dp.norm2() < (R*R) ){
                //if(verbosity>2) 
                //printf( "bond[%i,%i] r %g R %g(%g,%g) \n", i0-1, i,    dp.norm(), R, R0, Rj );
                found.push_back(i);
            }
            //else{printf( "NON bond[%i,%i] r %g R %g \n", i0-1, i, dp.norm(), R );}
        }
    }

    // ToDo:  R=0.5 worsk for hydrocarbons (HCNOF), R=0.65 works for silicon (SiH), we should assign it according to the bond types maybe ?
    int autoBonds( double R=-1.35, int i0=0, int imax=-1 ){
        //printf( "MM::Builder::autoBonds() R=%g bondByRvdW=%i \n", R, bondByRvdW );
        //if(verbosity>0){ printf( "MM::Builder::autoBonds() \n" ); }
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
            double Ri = getTypeRadius(A.type);   
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


        //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    void addBondedAtom(int ia, int ityp, bool bConf ){
        int ja=atoms.size();
        int npi=-1; if(bConf){npi=0;};
        insertAtom( ityp, 0, 0, npi, 0 );
        insertBond( {ia,ja}, 1 );
    }
    inline Vec3d pbcShift( Vec3i G ){ return lvec.a*G.a + lvec.b*G.b + lvec.c*G.c; }

    // find 
    int autoBondsPBC( double R=-1.35, int i0=0, int imax=-1, Vec3i npbc=Vec3iOne ){
        //printf( "MM::Builder::autoBondsPBC() \n" );
        //if(verbosity>0){ printf( "MM::Builder::autoBondsPBC() \n" );                             }
        //if(verbosity>1){ printf( "MM::Builder::autoBondsPBC() builder.lvec: \n" ); lvec.print(); };
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
            double Ri = getTypeRadius(A.type);
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
        //if(params)assignBondParams(ib);
    }

    //void addCap(int ia,Vec3d& hdir, Atom* atomj, int btype){
    Vec3d addEpair(int ia, const Vec3d& hdir, double l=-0.5, bool bInsertBond=true, bool bInsertAtoms=true, bool bUseCurrentType=false ){
        //printf( "addEpairsByPi[%i] h(%g,%g,%g) bInsertBond=%i bInsertAtoms=%i\n", ia, hdir.x,hdir.y,hdir.z, bInsertBond, bInsertAtoms );
        int ja=atoms.size();
        if(!bUseCurrentType){
            capAtom.type = itypEpair;
            if(params){ 
                int ityp = atoms[ia].type;
                capAtom.type = params->atypes[ityp].ePairType; 
            }
        }else{
            if(capAtom.type<itype_min){ capAtom.type = itypEpair; }
        }
        if(params){ if(l<0) l=params->atypes[capAtom.type].Ruff; }else{ l=-l; }
        
        capAtom.pos = atoms[ia].pos + hdir*l;
        double hn = hdir.norm();
        //if(params) { printf("MMFFBuilder::addEpair ia=%d hostT=%s capT=%s |h|=%g l=%g pos=(%g,%g,%g) host=(%g,%g,%g)\n", ia, params->atypes[atoms[ia].type].name, params->atypes[capAtom.type].name, hn, l, capAtom.pos.x,capAtom.pos.y,capAtom.pos.z, atoms[ia].pos.x,atoms[ia].pos.y,atoms[ia].pos.z); }
        //else       { printf("MMFFBuilder::addEpair ia=%d |h|=%g l=%g\n", ia, hn, l);}
        if(hn<1e-8){ printf("WARNING addEpair(): zero direction for ia=%d -> dummy placed on host!\n", ia); }
        if(bInsertAtoms)insertAtom(capAtom);
        if(bInsertBond){
            capBond.atoms.set(ia,ja);
            insertBond( capBond );
        }
        return capAtom.pos;
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

    void makeNeighs( int*& neighs, int perAtom ){
        int na = atoms.size();
        int ntot= na*perAtom;
        _allocIfNull( neighs, ntot );
        for(int i=0;i<ntot; i++){ neighs[i]=-1; }; // back neighbors
        for(int ia=0; ia<na; ia++ ){
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
                    if( jc==-1 ){ neighs[ ja*perAtom ]=ia; }
                }
            }
        }
    }




    // ======== Printing

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
        printf("printBondsOfAtom(%i): ", ia);
        for(int i=0; i<c.nbond; i++){
            int ib = c.neighs[i];
            Vec2i b = bonds[ib].atoms;
            printf( "(%i|%i,%i) ", ib, b.i, b.j );
        }
        printf("\n");
    }

    void printAtomConfs( bool bOmmitCap=true, bool bNeighs=false )const{
        printf(" # MM::Builder.printAtomConfs(na=%i,nc=%i) \n", atoms.size(), confs.size() );
        for(int i=0; i<atoms.size(); i++){ if( bOmmitCap && (atoms[i].iconf==-1) )continue;  if(bNeighs){printAtomNeighs(i);}else{printAtomConf(i);} puts(""); }
    }

    void printAtomTypes( )const{
        printf(" # MM::Builder.printAtomTypes(na=%i,nc=%i) \n", atoms.size(), confs.size() );
        for(int i=0; i<atoms.size(); i++){ 
            const Atom& A = atoms[i];
            int it=A.type;
            if(A.iconf>=0){
                AtomConf c = confs[A.iconf];
                printf( "atom[%4i] %ip %ie {%3i,%3i,%3i,%3i} type[%3i]=`%s`   \n", i, c.npi, c.ne, c.neighs[0],c.neighs[1],c.neighs[2],c.neighs[3], it, params->atypes[it].name ); 
            } else{
                printf( "atom[%4i] type[%3i]=`%s`   \n", i, it, params->atypes[it].name  ); 
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

    void printCappingTypes()const{
        printf(" # MM::Builder.printCappingTypes(nc=%i) \n", capping_types.size() );
        for(int it: capping_types ){ printf("capping_types[%3i] %8s \n", it, params->atypes[it].name ); }
    }

}; // MMFFBuilder


} // namespace MMFF

#endif // MMFFBuilder_h
        