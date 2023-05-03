#ifndef MMFFparams_h
#define MMFFparams_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "integerOps.h"

#include "molecular_utils.h"

#include <string>
#include <unordered_map>

inline uint64_t getBondTypeId( uint16_t at1, uint16_t at2, uint8_t order ){
    if (at1>at2){ _swap(at1,at2); }
    return pack64( at1, at2, order, 0 );
}

class BondType{ public:
    double length;
    double stiffness;
    uint16_t at1,at2;
    uint8_t  order;
    inline void sort(){ if (at1>at2){ _swap(at1,at2); }  }
    inline uint64_t getId(){ sort(); return pack64( at1, at2, order, 0 ); }
};

class AngleType{ public:
    double angle0;
    double stiffness;
    Vec3i  atoms; // a-b-c  or (lever1,fulcrum,lever2)
    inline void     sort(){ if (atoms.x>atoms.z){ _swap(atoms.x,atoms.z); }  }
    inline uint64_t getId(){ sort(); return pack64( atoms.b,atoms.a,atoms.c, 0 ); }
};

/*
class ElementType{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint32_t  color;
    //double  RvdW;
    //double  EvdW;
}

class AtomType_{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint32_t  color;
    double  RvdW;
    double  EvdW;
    double  Q;
    double Eaff,Ehard,Ra,eta;  // Charge Equilibration Params
}
*/


class ElementType{ public:
    char      name[4];
    uint8_t   iZ;
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint32_t  color;
    double    RvdW;
    double    EvdW;
    double    Quff;

    // QEq
    bool   bQEq;
    double Eaff,Ehard,Ra,eta;

    char* toString( char * str, bool bParams=false ){
        str         +=sprintf( str, "%s %i %i %i %x", name,  iZ, neval, valence, color );
        if(bParams)   sprintf( str, "%g %g %g   %g %g %g %g",  RvdW, EvdW, Quff,  Eaff,Ehard,Ra,eta  );
        return str;
    }

    void print(int i, bool bParams=false ){ 
        printf           ( "AtomType[%i,%s] %i(%i,%i) %x ", i,name,  iZ, neval, valence, color ); 
        if(bParams)printf( "REQuff(%g,%g,%g) QEq(%g,%g,%g)", RvdW, EvdW, Quff,  Eaff,Ehard,Ra,eta ); 
        printf( "\n"  ); 
    }

    inline uint8_t nepair(){ return (neval-valence)/2; };
    
};

class AtomType{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   valence;    // sum of bond orders of all bonds
    uint8_t   nepair;     // number of electron pairs
    uint8_t   npi;        // number of pi orbitals
    uint8_t   sym;        // sp3=0,sp2=1, sp1=2,  ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;
    double    Ruff;       // UFF natural bond radius
    double    RvdW;       // van der Waals potential minimum position
    double    EvdW;       // van der Waals potential minimum energy
    double    Qbase;      // base charge
    double    Hb;         // hydrogen bond correction
    int       parrent;
    int       element;
    int       ePairType;
    // ---- MMFF
    bool bMMFF;
    double  Ass,Asp, Kss,Ksp,Kep,Kpp;

    Vec3i     subTypes=Vec3iZero;  // sp1 sp2 sp3    // Q1 Q2 Q3 (polarized)
    
    //double  eRvdW,eEvdW,eQ;  // electron pair REQ parameters
    // ---- charge equlibration
    //bool   bQEq;
    //double Eaff,Ehard,Ra,eta;
    // ==== additional
    //double    piStiff;
    //double    electroneg;
    //double    polarizability;


    char* toString( char * str, bool bParams=false ){
        str         +=sprintf( str, "%s %i %i %i %i %i", name,  iZ,  valence, nepair,npi,sym );
        if(bParams)   sprintf( str, "%g %g %g %g   %g %g %g %g %g %g",  RvdW,EvdW,Qbase,Hb,  Ass,Asp, Kss,Ksp,Kep,Kpp );
        return str;
    }

    void print(int i, bool bParams=false ){ 
        printf           ( "AtomType[%i,%s] %i(%i,%i,%i,%i) ", i,name,  iZ,  valence, nepair,npi,sym  ); 
        if(bParams)printf( "REQH(%g,%g,%g,%g) MMFF(%g,%g,%g,%g,%g,%g)", RvdW,EvdW,Qbase,Hb,  Ass,Asp, Kss,Ksp,Kep,Kpp ); 
        printf( "\n"  ); 
    }

};

// ======================
// ====   MMFFparams
// ======================
static const int z2typ0[]{ 
    0, //0
    1, //H
    -1, //He
    -1, //Li
    -1, //Be
    -1, //B
    2, //C
    3, //N
    4, //O
    5  //F 
};

class MMFFparams{ public:

    // http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
    
    std::vector       <ElementType>        etypes;
    std::vector       <AtomType>           atypes;
    std::vector       <std::string    >    atomTypeNames;
    std::unordered_map<std::string,int>    elementTypeDict;
    std::unordered_map<std::string,int>    atomTypeDict;
    std::unordered_map<uint64_t,BondType>  bonds;
    std::unordered_map<uint64_t,AngleType> angles;


    double default_bond_length      = 2.0;
    double default_bond_stiffness   = 1.0;
    //Quat4d  default_REQ            = {1.487, 0.0006808, 0.0, 0.};  // Hydrogen
    Quat4d  default_REQ              = {1.500, 0.0005000, 0.0, 0.};  // Hydrogen like

    bool reportIfMissing=true;
    bool exitIfMissing  =true;

    void initDefaultAtomTypeDict(){
        makeDefaultAtomTypeDict( atomTypeNames, atomTypeDict );
    }

    void printAtomTypeDict(){
        for(int i=0; i<atomTypeNames.size(); i++){ printf( "AtomType[%i] %s %i\n", i, atypes[i].name, atomTypeDict[atypes[i].name] );  };
    }
    
    void printElementTypeDict(){
        for(int i=0; i<atomTypeNames.size(); i++){ printf( "ElementType[%i] %s %i\n", i,  etypes[i].name, elementTypeDict[etypes[i].name] );  };
    }

    int getAtomType(const char* s, bool bErr=true){
        //printf( "getAtomType(%s) bErr=%i \n", s, bErr );
        auto found = atomTypeDict.find(s);
        if(found==atomTypeDict.end()){ 
            if(bErr){ printf( "ERROR: MMFFparams::getAtomType(%s) not found !!! => exit() \n", s ); printAtomTypeDict(); exit(0); }
            return -1; 
        }
        return found->second;
    }

    int getElementType(const char* s, bool bErr=true){
        //printf( "getAtomType(%s) bErr=%i \n", s, bErr );
        auto found = elementTypeDict.find(s);
        if(found==elementTypeDict.end()){ 
            if(bErr){ printf( "ERROR: MMFFparams::getElementType(%s) not found !!! => exit() \n", s ); printAtomTypeDict(); exit(0); }
            return -1; 
        }
        return found->second;
    }

    inline ElementType* elementOfAtomType( int it ){ return &etypes[atypes[it].element]; }

    AtomType* getRootParrent(AtomType* t, int nrecur=0){
        if(nrecur>10){ printf("ERROR in MMFFparams.getRootParrent() rootParrent of type(%s) not found in %i recursions => Exit() \n", t->name, nrecur ); exit(0); }
        if( t->parrent==0 ) return t;
        if( (t->parrent<0)||(t->parrent>=atypes.size()) ){ printf("ERROR in MMFFparams.getRootParrent() type(%s).parrent==%i => Exit() \n", t->name,t->parrent ); exit(0); }
        AtomType* par = &atypes[t->parrent];
        return getRootParrent(par,nrecur+1); // recursion
    }

    void string2AtomType(const char * str, AtomType& atyp ){
        char      parent_name[8];
        char      epair_name[8];
        char      element_name[4];
        int iZ_, neval_, valence_, nepair_, npi_, sym_;
        //char sclr[6];           1           2           3        4        5        6        7         8            9           10            11           12        13          14        15         16         17         18        19         20           21          22        23
        //int nret = sscanf( str, " %s         %s           %i      %i       %i        %i      %lf       %lf          %x            %lf          %lf          %lf       %lf         %lf       %lf        %lf       %lf        %lf       %lf         %i           %lf         %lf      %lf ", 
        //                         atyp.name, parent_name, &iZ_, &neval_, &valence_, &sym_,  &atyp.RvdW, &atyp.EvdW, &atyp.color,   &atyp.Eaff, &atyp.Ehard, &atyp.Ra, &atyp.eta,   &atyp.Ass,&atyp.Asp, &atyp.Kss,&atyp.Ksp,&atyp.Kep,&atyp.Kpp, &atyp.ne,&atyp.eRvdW,&atyp.eEvdW,&atyp.eQ );

        // char      name[8];
        // uint8_t   iZ;         // proton number
        // uint8_t   valence;    // sum of bond orders of all bonds
        // uint8_t   ne;         // number of electron pairs
        // uint8_t   npi;        // number of pi orbitals
        // uint8_t   sym;        // sp3=0,sp2=1, sp1=2,  ... tetrahedral, triangle, linear, kink, square, octahedral
        // uint32_t  color;
        // double    RvdW;       // van der Waals potential minimum position
        // double    EvdW;       // van der Waals potential minimum energy
        // double    Qbase;      // base charge
        // double    Hb;         // hydrogen bond correction
        // int       parrent;
        // int       element;
        // int       ePairType;
        // // ---- MMFF
        // bool bMMFF;
        // double  Ass,Asp, Kss,Ksp,Kep,Kpp;
        //char sclr[6];           1           2           3              4             5          6        7       8          9         10        11           12         13            14        15         16         17         18        19      
        int nret = sscanf( str, " %s         %s          %s              %s            %i         %i      %i      %i        %lf         %lf       %lf         %lf        %lf            %lf       %lf        %lf       %lf        %lf       %lf   ", 
                                 atyp.name, parent_name, element_name, epair_name,  &valence_, &nepair_, &npi_, &sym_,    &atyp.Ruff, &atyp.RvdW, &atyp.EvdW, &atyp.Qbase, &atyp.Hb,    &atyp.Ass,&atyp.Asp, &atyp.Kss,&atyp.Ksp,&atyp.Kep,&atyp.Kpp );
        atyp.valence=valence_; atyp.nepair=nepair_; atyp.npi=npi_; atyp.sym=sym_;
        if(atypes.size()!=0){
            int iet  = getElementType(element_name, false); if(iet <0){ printf("ERROR in MMFFparams::string2AtomType(): cannot find elementType (%s) of type(%s) => Exit() \n", element_name, atyp.name ); exit(0);  };
            int ipar = getAtomType   (parent_name,  false); if(ipar<0){ printf("ERROR in MMFFparams::string2AtomType(): cannot find parrent type(%s) of type(%s) => Exit() \n", parent_name,  atyp.name ); exit(0); };
            int iept = getAtomType   (epair_name,   false); if(iept<0){ printf("ERROR in MMFFparams::string2AtomType(): cannot find epair   type(%s) of type(%s) => Exit() \n", epair_name,   atyp.name ); exit(0); };
            const ElementType& et = etypes[iet];
            atyp.iZ    = et.iZ;
            atyp.color = et.color; 
            atyp.parrent  =ipar;
            atyp.ePairType=iept;
        }
        if(nret<19){ atyp.bMMFF = false; atyp.Ass=0;atyp.Asp=0;atyp.Kss=0;atyp.Ksp=0;atyp.Kep=0;atyp.Kpp=0; }else{ atyp.bMMFF = true; }
        atyp.subTypes=Vec3iZero;
        //printf( "AtomType: %s iZ %i ne %i nb %i sym %i RE(%g,%g) %x \n", atyp.name,  atyp.iZ,   atyp.neval, atyp. valence,   atyp.sym,    atyp.RvdW, atyp.EvdW,   atyp.color );
        //char ss[256]; printf("%s\n", toString(ss) );
    }

    void string2ElementType(const char * str, ElementType& atyp ){
        // char      name[4];
        // uint8_t   iZ;
        // uint8_t   neval;      // number of valence electrons
        // uint8_t   valence;    // sum of bond orders of all bonds
        // uint32_t  color;
        // double    RvdW;
        // double    EvdW;
        // // QEq
        // bool   bQEq;
        // double Eaff,Ehard,Ra,eta;
        // // UFF
        // bool   bUFF;
        // double uff_l0,uff_q;
        
        int iZ_, neval_, valence_, sym_;
        //char sclr[6];           1           2        3       4         5              6          7              8              9            10       11       12           
        int nret = sscanf( str, " %s          %i      %i       %i        %x             %lf       %lf            %lf           %lf          %lf        %lf     %lf", 
                                 atyp.name,  &iZ_, &neval_, &valence_,   &atyp.color, &atyp.RvdW, &atyp.EvdW, &atyp.Quff,  &atyp.Eaff, &atyp.Ehard, &atyp.Ra, &atyp.eta  );
        atyp.iZ=iZ_; atyp.neval=neval_; atyp.valence=valence_;
        const int nretmin=8;
        if(nret<nretmin){ printf( "ERROR in MMFFparams::string2ElementType(iZ=%i,%s) not complete (nret(%i)<nretmin(%i)) => Exit() \n", iZ_, atyp.name, nret, nretmin  ); printf("%s\n", str ); exit(0); }
        if(nret<12     ){ atyp.bQEq=false;  atyp.Eaff=0; atyp.Ehard=0; atyp.Ra=0; atyp.eta=0;  }else{ atyp.bQEq=true; }
    }

    int loadAtomTypes(const char * fname, bool exitIfFail=true){
        //printf( "loadAtomTypes %s \n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;

        AtomType atyp;
        int i=0;
        for(i=0; i<10000; i++){
        //for(int i; i<0xFFFF; i++){
            //printf( "loadAtomTypes %i \n", i );
            line = fgets( buff, 1024, pFile );
            if(line==NULL)  break;
            //printf( "loadAtomTypes[%i] line=%s", i, line );
            if(line[0]=='#')continue;
            //atyp.fromString( line );
            string2AtomType( line, atyp );
            atypes.push_back(atyp);
            atomTypeNames.push_back( atyp.name );
            atomTypeDict[atyp.name] = atypes.size()-1;
            //char str[1000];
            //atyp.toString( str );
            //printf( "%i %s %i %s \n", i, atyp.name, atypNames[atyp.name], str );
        }
        return i;
    }

    int loadElementTypes(const char * fname, bool exitIfFail=true){
        //printf( "loadAtomTypes %s \n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        ElementType etyp;
        int i=0;
        for(i=0; i<200; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL)  break;
            if(line[0]=='#')continue;
            string2ElementType( line, etyp );
            etypes.push_back(etyp);
            elementTypeDict[etyp.name] = etypes.size()-1;
        }
        return i;
    }

    inline void assignSubTypes( AtomType& t ){
        //printf( "assignSubTypes %s(iZ=%i)\n", t.name, t.iZ );
        char tmp_name[8];
        const char* ssub[3]{"3","2","1"};
        for(int i=0;i<3;i++){
            sprintf( tmp_name, "%s_%s", t.name, ssub[i] );
            //printf( "assignSubTypes `%s`(iZ=%i)[%i] %s\n", t.name, t.iZ, i, tmp_name );
            int it = getAtomType(tmp_name, false);
            //printf( "assignSubTypes %s(iZ=%i)[%i] %s=%i\n", t.name, t.iZ, i, tmp_name, it );
            if(it<0)continue;
            t.subTypes.array[i] = it;
            //printf( "assignSubTypes %s(iZ=%i)[%i] %s=%i\n", t.name, t.iZ, tmp_name, it );
        }
    }
    inline void assignAllSubTypes(){
        int n=atypes.size();
        std::vector<bool> doIt(256,true); // Warrning : we assume maximum proton number 256
        for(int i=0;i<n;i++){
            AtomType& t = atypes[i];
            //printf( "DEBUG assignAllSubTypes() t.iZ %i doIt.size()= %i \n", t.iZ, doIt.size() );
            //if( t.iZ>=doIt.size() ){ printf("ERROR: atype[%i] t.iZ(%i) > =doIt.size(%i) \n", i, t.iZ, doIt.size()  ); }
            if(doIt[t.iZ]){ assignSubTypes(t); doIt[t.iZ]=false; }
        }
    }

    inline void assignRE( int ityp, Quat4d& REQ, bool bSqrtE=false )const{
        REQ.x    = atypes[ityp].RvdW;
        double e = atypes[ityp].EvdW;
        if(bSqrtE) e=sqrt(e);
        REQ.y = e;
    }

    void assignREs( int n, int * itypes, Quat4d * REQs, bool bSqrtE=false, bool bQ0=false )const{
        //printf( "assignREs(%i) %li \n", n, itypes );
        for(int i=0; i<n; i++){
            //printf( " assignREs[%i] %i \n", i, itypes[i] );
            assignRE( itypes[i], REQs[i], bSqrtE );
            if(bQ0) REQs[i].z=0;
        }
    }

    void assignQEq( int n, int* itypes, double* affins, double* hards ){
        for(int i=0; i<n; i++){
            int ityp = itypes[i];
            int iet  = atypes[ityp].element; 
            affins[i]=etypes[iet].Eaff;
            hards [i]=etypes[iet].Ehard;
        }
    }

    int loadBondTypes(const char * fname, bool exitIfFail=true){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char * line;
        BondType bt;
        //line = fgets( buff, 1024, pFile ); //printf("%s",line);
        //sscanf( line, "%i %i\n", &n );
        int i=0;
        for( i; i<1000; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            //printf("%s",line);
            sscanf(  line, "%i %i %i %lf %lf\n", &bt.at1, &bt.at2, &bt.order, &bt.length, &bt.stiffness );
            //printf(        "%i %i %i %lf %lf\n",  bt.at1,  bt.at2,  bt.order,  bt.length,  bt.stiffness );
            uint64_t id = bt.getId();
            //uint64_t id = getBondTypeId( bt.at1, bt.at2, bt.order );
            //printf( "loadBondTypes[%i] iZ(%i,%i|%i) id=%i \n", i, bt.at1, bt.at2, bt.order, id );
            //bt.at1--; bt.at2--;
            bonds[id]=bt;
        }
        return i;
    }

    void printAngle(const AngleType& at){
        printf( "%s %s %s %g %g\n", atypes[at.atoms.a].name, atypes[at.atoms.b].name, atypes[at.atoms.c].name, at.angle0, at.stiffness );
    }

    void printAtomTypes(bool bParams){
        for(int i=0; i<atypes.size(); i++ ){ atypes[i].print(i, bParams );  }
    }

    int loadAgnleType(const char * fname, bool exitIfFail=true){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char * line;
        AngleType ang;
        int i=0;
        char name[3][8];
        for( i; i<1000; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            sscanf(  line,            "%s %s %s %lf %lf\n",  name[0], name[1], name[2], &ang.angle0, &ang.stiffness );
            //printf( "loadAgnleType[%i] %s %s %s %lf %lf\n",i,name[0], name[1], name[2],  ang.angle0,  ang.stiffness );
            ang.atoms.a = getAtomType(name[0]);
            ang.atoms.b = getAtomType(name[1]);
            ang.atoms.c = getAtomType(name[2]);
            uint64_t id = ang.getId();
            auto found = angles.find(id);
            if( found != angles.end() ){ printf( "WARRNIMG!!! loadAgnleType() angleType[%i] same as ", i); printAngle(found->second); };
            angles[id]=ang;
        }
        return i;
    }

    bool getBondParams( int atyp1, int atyp2, int btyp, double& l0, double& k )const{
        uint64_t id  = getBondTypeId( atypes[atyp1].iZ, atypes[atyp2].iZ, btyp );
        //printf( "(%i,%i,%i) -> %i \n", atypes[atyp1].iZ, atypes[atyp2].iZ, btyp, id );
        //uint64_t id  = getBondTypeId( atyp1, atyp2, btyp );
        auto it = bonds.find(id);
        if   ( it == bonds.end() ){ 
            if(reportIfMissing){ printf("WARRNING!!! getBondParams() missing atyps(%i,%i) iZs(%i,%i) o(%i)id(%i)=>l0 %g k %g \n", atyp1, atyp2, atypes[atyp1].iZ, atypes[atyp2].iZ, btyp, id, default_bond_length, default_bond_stiffness ); };
            if(exitIfMissing){ printf("=> exit(0)\n");exit(0); };
            l0=default_bond_length; k=default_bond_stiffness; return false;
        }else{ 
            l0=it->second.length;   k=it->second.stiffness;   return true; 
        }
    }

    void fillBondParams( int nbonds, Vec2i * bond2atom, int * bondOrder, int * atomType, double * bond_0, double * bond_k ){
        //printf("fillBondParams: %i\n", nbonds);
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            getBondParams( atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
            //printf( "%i (%i %i) %i %g %g \n", i, atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
        }
    }

    void init(const char* felementTypes=0, const char* fatomtypes=0, const char* fbondtypes=0, const char* fagnletypes=0){
        //if(verbosity>0) 
        //printf("MMFFparams::init(%s,%s,%s)\n", fatomtypes, fbondtypes, fagnletypes );
        if(felementTypes ){
            loadElementTypes( felementTypes );
        }
        if(fatomtypes ){
            loadAtomTypes( fatomtypes );
            assignAllSubTypes();
        }
        if(fbondtypes )loadBondTypes( fbondtypes  );
        if(fagnletypes)loadAgnleType( fagnletypes );
    }

    bool cellFromString( char* s, Mat3d& lvec )const{
        char c[3]; Mat3d M;
        int n = sscanf( s, "%c%c%c %lf %lf %lf   %lf %lf %lf   %lf %lf %lf", c,c+1,c+2, &(M.a.x),&(M.a.y),&(M.a.z),   &(M.b.x),&(M.b.y),&(M.b.z),   &(M.c.x),&(M.c.y),&(M.c.z) );
        //printf( "DEBUG cellFromString() n=%i %c%c%c (%g,%g,%g)(%g,%g,%g)(%g,%g,%g) \n", n, c[0],c[1],c[2], M.a.x,M.a.y,M.a.z,   M.b.x,M.b.y,M.b.z,   M.c.x,M.c.y,M.c.z );
        if( (n==12) && (c[0]=='l')&&(c[1]=='v')&&(c[2]=='s') ){
            lvec=M;
            //printf( "DEBUG cellFromString() lvec (%g,%g,%g)(%g,%g,%g)(%g,%g,%g) \n",  lvec.a.x,lvec.a.y,lvec.a.z,   lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
            return true;
        }
        return false;
    }

    int loadXYZ(const char* fname, int& natoms, Vec3d** apos_, Quat4d** REQs_=0, int** atype_=0, int** npis_=0, Mat3d* lvec=0, int verbosity=0 )const{
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        const int nbuf=1024;
        char buff[nbuf];
        char * line;
        int nl;
        line = fgets( buff, nbuf, pFile );
        sscanf( line, "%i \n", &natoms );
        Vec3d* apos  =_allocPointer( apos_,  natoms );
        Quat4d* REQs =_allocPointer( REQs_,  natoms );
        int*   atype =_allocPointer( atype_, natoms );
        int*   npis  =_allocPointer( npis_,  natoms );
        line = fgets( buff, nbuf, pFile ); // comment
        int ret=0;
        if(lvec){ if( cellFromString( buff, *lvec ) ){ ret=1; }else{ printf("WARRNING: lvec not read from %s \n", fname ); } }
        double Q,H;
        for(int i=0; i<natoms; i++){
            char at_name[8];
            double junk; 
            int npi;
            line = fgets( buff, nbuf, pFile );  //printf("%s",line);
            //int nret = sscanf( line, "%s %lf %lf %lf %lf \n", at_name, &apos[i].x, &apos[i].y, &apos[i].z, &Q, &npi );
            int nret = sscanf( line, "%s %lf %lf %lf %lf \n", at_name, &apos[i].x, &apos[i].y, &apos[i].z, &Q, &H, &npi  );
            if( nret < 5 ){ Q=0; };
            if( nret < 6 ){ H=0; };
            if( nret < 7 ){ npi=-1; };
            if(npis){ npis[i] =npi; };
            auto it = atomTypeDict.find( at_name );
            if( it != atomTypeDict.end() ){
                int ityp=it->second;
                if(atype_)atype[i] = ityp;
                if(REQs_){ assignRE( ityp, REQs[i], true ); REQs[i].z=Q; REQs[i].w=H; }
            }else{
                if(atype_)atype[i] = -1;
                if(REQs_)REQs[i]  = default_REQ;
            }
        }
        fclose(pFile);
        return ret;
    }

    void writeXYZ( FILE* pfile, int n, const int* atyps, const Vec3d* apos, const char* comment="#comment", const Quat4d* REQs=0, bool just_Element=true, int npi=0 ){
        //printf( "MMFFparams::writeXYZ() n=%i REQs=%li \n", n, (long)REQs );
        fprintf(pfile, "%i\n", n+npi );
        fprintf(pfile, "%s \n", comment );
        for(int i=0; i<n; i++){
            //printf( "DEBUG writeXYZ()[%i] \n", i );
            int ityp   = atyps[i];
            const Vec3d&  pi = apos[i];
            const char* symbol; 
            bool byName = true;
            if(just_Element){ 
                int iZ = atypes[ityp].iZ;
                if( iZ<=11 ){
                    int it = z2typ0[iZ];
                    symbol = atomTypeNames[it].c_str();
                    //printf( "atom[%i] ityp %i iZ %i it %i `%s` \n", i, ityp, iZ, it, symbol );
                    byName = false;
                }
            }
            if(byName){ symbol = atomTypeNames[ityp].c_str(); }
            //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
            if(REQs){ fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f     %10.6f\n", symbol, pi.x,pi.y,pi.z, REQs[i].z ); }
            else    { fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n"          , symbol, pi.x,pi.y,pi.z            ); }
        }
        for(int i=0; i<npi; i++){
            const Vec3d&  pi = apos[i] + apos[i+n];
            fprintf( pfile, "Pi   %15.10f   %15.10f   %15.10f \n", pi.x,pi.y,pi.z            );

        }
    }

    int saveXYZ( const char * fname, int n, const int* atyps, const Vec3d* apos, const char* comment="#comment", const Quat4d* REQs=0, const char* mode="w", bool just_Element=true ){
        FILE* pfile = fopen(fname, mode );
        //printf( "saveXYZ(%s) \n", fname );
        if( pfile == NULL ) return -1;
        writeXYZ( pfile, n, atyps, apos, comment, REQs, just_Element );
        fclose(pfile);
        return n;
    }

};

#endif
