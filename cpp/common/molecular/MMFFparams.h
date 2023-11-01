#ifndef MMFFparams_h
#define MMFFparams_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "integerOps.h"

#include "constants.h"

#include "molecular_utils.h"

#include <string>
#include <unordered_map>


class BondType{ public:
    double length;
    double stiffness;
    Vec2i   atoms;
    uint8_t order;
    inline bool            sort (){ if(atoms.x>atoms.y){ _swap(atoms.x,atoms.y); return true; } return false;  }
    inline static uint64_t getId( uint16_t at1, uint16_t at2, uint8_t order ){ if (at1>at2){ _swap(at1,at2); } return pack64( at1, at2, order, 0 ); }
    inline uint64_t        id   (){ return getId( atoms.x, atoms.y, order ); }
};

class AngleType{ public:
    double angle0;
    double stiffness;
    Vec3i  atoms; // a-b-c  or (lever1,fulcrum,lever2)
    inline bool            sort (){ if (atoms.x>atoms.z){ _swap(atoms.x,atoms.z); return true; } return false;  }
    inline static uint64_t getId(  uint16_t a, uint16_t b, uint16_t c ){ if (a>c){ _swap(a,c); } return pack64( b,a,c, 0 );  }
    inline uint64_t        id   (){ return getId( atoms.x, atoms.y, atoms.z ); }

};

class DihedralType{ public:
    Quat4i atoms;
    int    bo; // bond order of central atoms
    int    n; 
    double k;
    double ang0;

    inline bool            sort (){ if (atoms.y>atoms.z){ _swap(atoms.y,atoms.z); _swap(atoms.x,atoms.w); return true; } return false; }
    inline static uint64_t getId(  uint16_t a, uint16_t b, uint16_t c, uint16_t d, int order ){ if (b>c){ _swap(b,c); _swap(a,d); } return pack64( b,c,a,d+order);  }
    inline        uint64_t id   (){ return getId(atoms.x,atoms.y,atoms.z,atoms.w,bo); }

};


class ElementType{ public:
    char      name[4];    // symbol
    uint8_t   iZ;         // atomic number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint32_t  color;      // color
    double    RvdW;       // LJ distance parameter
    double    EvdW;       // LJ energy parameter
    double    Quff;       // effective charge in UFF
    bool      bQEq;       // flag to see if QEq params are provided
    double    Eaff;       // electronegativity
    double    Ehard;      // chemical hardness
    double    Ra;         // atomic size
    double    eta;        // valence orbital exponent

    char* toString( char * str, bool bParams=false )const{
        str         +=sprintf( str, "%s %i %i %i %x", name,  iZ, neval, valence, color );
        if(bParams)   sprintf( str, "%g %g %g   %g %g %g %g",  RvdW, EvdW, Quff,  Eaff,Ehard,Ra,eta  );
        return str;
    }

    void print(int i, bool bParams=false )const{ 
        printf           ( "AtomType[%i,%s] %i(%i,%i) %x ", i,name,  iZ, neval, valence, color ); 
        if(bParams)printf( "REQuff(%g,%g,%g) QEq(%g,%g,%g)", RvdW, EvdW, Quff,  Eaff,Ehard,Ra,eta ); 
        printf( "\n"  ); 
    }

    inline uint8_t nepair(){ return (neval-valence)/2; };
    
};

class AtomType{ public:
    char      name[8];    // symbol
    uint8_t   iZ;         // atomic number
    uint8_t   valence;    // sum of bond orders of all bonds
    uint8_t   nepair;     // number of electron pairs
    uint8_t   npi;        // number of pi orbitals
    // TBD to check about this, it can be useful for UFF...
    uint8_t   sym;        // sp3=0,sp2=1, sp1=2,  ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;      // color
    double    Ruff;       // UFF natural bond radius
    double    RvdW;       // LJ distance parameter
    double    EvdW;       // LJ energy parameter
    double    Qbase;      // atomic charge
    double    Hb;         // hydrogen bond correction
    int       parrent;    // parent type
    int       element;    // corresponding element type
    int       ePairType;  // type of lone pair owned
    bool      bMMFF;      // flag to see if MMFF params are provided
    double    Ass;        // equilibrium angle value for sigma-sigma interaction (angle bending)
    double    Asp;        // equilibrium angle value for sigma-pi interaction
    double    Kss;        // force constant for sigma-sigma interaction (angle bending)
    double    Ksp;        // force constant for sigma-pi interaction
    double    Kep;        // force constant for pi-lone pair interaction
    double    Kpp;        // force constant for lone pair-lone pair interaction
    // TBD to get rid of it?
    Vec3i     subTypes=Vec3iZero;  // sp1 sp2 sp3    // Q1 Q2 Q3 (polarized)

    char* toString( char * str, bool bParams=false )const{
        str         +=sprintf( str, "%s %i %i %i %i %i", name,  iZ,  valence, nepair,npi,sym );
        if(bParams)   sprintf( str, "%g %g %g %g   %g %g %g %g %g %g",  RvdW,EvdW,Qbase,Hb,  Ass,Asp, Kss,Ksp,Kep,Kpp );
        return str;
    }

    void print(int i, bool bParams=false )const{ 
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

    // http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html LINK EXPIRED
    
    std::vector       <ElementType>        etypes;
    std::vector       <AtomType>           atypes;
    std::vector       <std::string    >    atomTypeNames;
    std::unordered_map<std::string,int>    elementTypeDict;
    std::unordered_map<std::string,int>    atomTypeDict;

    std::vector<BondType>                      bonds;
    std::unordered_map<uint64_t,BondType*>     bonds_;
    std::unordered_map<std::string,int>        bondDict;

    std::vector<AngleType>                     angles;
    std::unordered_map<uint64_t,AngleType*>    angles_;
    std::unordered_map<std::string,int>        angleDict;

    std::vector<DihedralType>                  dihedrals;
    std::unordered_map<uint64_t,DihedralType*> dihedrals_;
    std::unordered_map<std::string,int>        dihedralDict;


    double default_bond_length      = 2.0;
    double default_bond_stiffness   = 1.0;
    //Quat4d  default_REQ           = {1.487, 0.0006808, 0.0, 0.};  // Hydrogen
    Quat4d  default_REQ             = {1.500, 0.0005000, 0.0, 0.};  // Hydrogen like

    bool echoTry        =true;
    bool reportIfMissing=true;
    bool exitIfMissing  =true;

    //////////////////////////////////////////////////////////////////////////////////
    // LOAD TABLES WITH ELEMENT, ATOM, BOND, ANGLE AND DIHEDRAL TYPE SPECIFICATIONS //
    //////////////////////////////////////////////////////////////////////////////////
    // extract variables from one line of the ElementTypes file
    void string2ElementType(const char * str, ElementType& etyp ){
        //char      name[4];    // symbol
        //uint8_t   iZ;         // atomic number
        //uint8_t   neval;      // number of valence electrons
        //uint8_t   valence;    // sum of bond orders of all bonds
        //uint32_t  color;      // color
        //double    RvdW;       // LJ radius
        //double    EvdW;       // LJ energy well 
        //double    Quff;       // effective charge in UFF
        //bool   bQEq;          // flag to see if QEq params are provided
        //double Eaff;          // electronegativity
        //double Ehard;         // chemical hardness
        //double Ra;            // atomic size
        //double eta;           // valence orbital exponent
        //                       1          2         3            4              5            6           7           8           9           10           11        12           
        int nret = sscanf( str, "%s         %i        %i           %i             %x           %lf         %lf         %lf         %lf         %lf          %lf       %lf", 
                                 etyp.name, &etyp.iZ, &etyp.neval, &etyp.valence, &etyp.color, &etyp.RvdW, &etyp.EvdW, &etyp.Quff, &etyp.Eaff, &etyp.Ehard, &etyp.Ra, &etyp.eta );
        const int nretmin=8;
        if(nret<nretmin){ printf( "ERROR in MMFFparams::string2ElementType(iZ=%i,%s) not complete (nret(%i)<nretmin(%i)) => Exit() \n", etyp.iZ, etyp.name, nret, nretmin  ); printf("%s\n", str ); exit(0); }
        if(nret<12     ){ etyp.bQEq=false; etyp.Eaff=0; etyp.Ehard=0; etyp.Ra=0; etyp.eta=0; }else{ etyp.bQEq=true; }
    }

    // read and store element types
    int loadElementTypes(const char * fname, bool exitIfFail=true){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char * line;
        ElementType etyp;
        int i;
        for(i=0; ; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL)  break;
            if(line[0]=='#')continue;
            string2ElementType( line, etyp );
            etypes.push_back(etyp);
            if( !elementTypeDict.insert({ etyp.name, etypes.size()-1} ).second ){ printf("WARNING: elementType[%i](%s) is duplicated !!! => Ignore \n", i, line ); exit(0); };
            printf("loadElementTypes[%i] name='%s' iZ=%i neval=%i valence=%i\n", etypes.size(), etyp.name, etyp.iZ, etyp.neval, etyp.valence );
        }
        return i;
    }

    // extract variables from one line of the AtomTypes file
    void string2AtomType(const char * str, AtomType& atyp ){
        //char      name[8];    // symbol
        //uint8_t   iZ;         // atomic number
        //uint8_t   valence;    // sum of bond orders of all bonds
        //uint8_t   nepair;     // number of electron pairs
        //uint8_t   npi;        // number of pi orbitals
        //uint8_t   sym;        // sp3=0,sp2=1, sp1=2,  ... tetrahedral, triangle, linear, kink, square, octahedral
        //uint32_t  color;      // color
        //double    Ruff;       // UFF natural bond radius
        //double    RvdW;       // LJ distance parameter
        //double    EvdW;       // LJ energy parameter
        //double    Qbase;      // atomic charge
        //double    Hb;         // hydrogen bond correction
        //int       parrent;    // parent type
        //int       element;    // corresponding element type
        //int       ePairType;  // type of lone pair owned
        //bool      bMMFF;      // flag to see if MMFF params are provided
        //double    Ass;        // equilibrium angle value for sigma-sigma interaction (angle bending)
        //double    Asp;        // equilibrium angle value for sigma-pi interaction
        //double    Kss;        // force constant for sigma-sigma interaction (angle bending)
        //double    Ksp;        // force constant for sigma-pi interaction
        //double    Kep;        // force constant for pi-lone pair interaction
        //double    Kpp;        // force constant for lone pair-lone pair interaction
        //Vec3i     subTypes=Vec3iZero;  // sp1 sp2 sp3    // Q1 Q2 Q3 (polarized)
        char parent_name[8];
        char epair_name[8];
        char element_name[4];
        //                       1          2            3             4           5              6             7          8          9           10          11          12           13        14         15         16         17         18         19
        int nret = sscanf( str, "%s         %s           %s            %s          %i             %i            %i         %i         %lf         %lf         %lf         %lf          %lf       %lf        %lf        %lf        %lf        %lf        %lf", 
                                 atyp.name, parent_name, element_name, epair_name, &atyp.valence, &atyp.nepair, &atyp.npi, &atyp.sym, &atyp.Ruff, &atyp.RvdW, &atyp.EvdW, &atyp.Qbase, &atyp.Hb, &atyp.Ass, &atyp.Asp, &atyp.Kss, &atyp.Ksp, &atyp.Kep, &atyp.Kpp );
        if(atypes.size()!=0){
            int iet  = getElementType(element_name, false); if(iet <0){ printf("ERROR in MMFFparams::string2AtomType(): cannot find elementType (%s) of type(%s) => Exit() \n", element_name, atyp.name ); exit(0);  };
            int ipar = getAtomType   (parent_name,  false); if(ipar<0){ printf("ERROR in MMFFparams::string2AtomType(): cannot find parrent type(%s) of type(%s) => Exit() \n", parent_name,  atyp.name ); exit(0); };
            int iept = getAtomType   (epair_name,   false); if(iept<0){ printf("ERROR in MMFFparams::string2AtomType(): cannot find epair   type(%s) of type(%s) => Exit() \n", epair_name,   atyp.name ); exit(0); };
            const ElementType& et = etypes[iet];
            atyp.element   = iet;
            atyp.parrent   = ipar;
            atyp.ePairType = iept;
            atyp.iZ        = et.iZ;
            atyp.color     = et.color; 
        }
        atyp.subTypes=Vec3iZero;
        const int nretmin=13;
        if(nret<nretmin){ printf( "ERROR in MMFFparams::string2AtomType(%s) not complete (nret(%i)<nretmin(%i)) => Exit() \n", atyp.name, nret, nretmin  ); printf("%s\n", str ); exit(0); }
        if(nret<19){ atyp.bMMFF = false; atyp.Ass=0; atyp.Asp=0; atyp.Kss=0; atyp.Ksp=0; atyp.Kep=0; atyp.Kpp=0; }else{ atyp.bMMFF = true; }
    }

    // read and store atom types
    int loadAtomTypes(const char * fname, bool exitIfFail=true){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char * line;
        AtomType atyp;
        int i;
        for(i=0; ; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL)  break;
            if(line[0]=='#')continue;
            string2AtomType( line, atyp );
            atypes.push_back(atyp);
            atomTypeNames.push_back( atyp.name );
            if( !atomTypeDict.insert({atyp.name, atypes.size()-1}).second ){ printf("WARNING: atomType[%i](%s) is duplicated !!! => Ignore \n", i, line ); exit(0); };
            printf("loadAtomTypes[%i] name='%s' valence=%i nepair=%i npi=%i sym=%i\n", atypes.size(), atyp.name, atyp.valence, atyp.nepair, atyp.npi, atyp.sym );
        }
        return i;
    }



    int loadBondTypes(const char * fname, bool exitIfFail=true, bool bWarnFlip=true){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char names[2][8];
        char name[64];
        char * line;
        BondType bt;
        //line = fgets( buff, 1024, pFile ); //printf("%s",line);
        //sscanf( line, "%i %i\n", &n );
        int i=0;
        for( i; i<1000; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            if(line[0]=='#') continue;
            //printf("%s",line);
            sscanf(  line, "%s %s %i %lf %lf\n", names[0], names[1], &bt.order, &bt.length, &bt.stiffness );

            bt.atoms.x = getAtomType(names[0]);
            bt.atoms.y = getAtomType(names[1]);
            if( bt.sort() && bWarnFlip ){ printf("WARRNING: bondType[%i](%s) is flipped in %s\n", i, line, fname ); };
            sprintf( name, "%s-%s-%i", atypes[bt.atoms.x].name , atypes[bt.atoms.y].name, bt.order );
            bonds.push_back( bt );
            if( !bondDict.insert({ name, bonds.size()-1} ).second ){ printf("WARRNING: bondType[%i](%s) is duplicated !!! => Ignore \n", i, line ); exit(0); };

            //printf(        "%i %i %i %lf %lf\n",  bt.atoms.x,  bt.atoms.y,  bt.order,  bt.length,  bt.stiffness );
            //uint64_t id = bt.getId();
            //uint64_t id = getBondTypeId( bt.at1, bt.at2, bt.order );
            //printf( "loadBondTypes[%i] iZ(%i,%i|%i) id=%i \n", i, bt.at1, bt.at2, bt.order, id );
            //bt.at1--; bt.at2--;
            //bonds[id]=bt;
        }
        return i;
    }

    int loadAngleTypes(const char * fname, bool exitIfFail=true, bool bWarnFlip=true ){
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char names[3][8];
        char name[64];
        char * line;
        AngleType ang;
        int i=0;
        for( i; i<1000; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            if(line[0]=='#') continue;
            sscanf(  line,            "%s %s %s %lf %lf\n",     names[0], names[1], names[2], &ang.angle0, &ang.stiffness );
            //printf( "loadAgnleTypes[%i] %s %s %s %lf %lf\n",i,names[0], names[1], names[2],  ang.angle0,  ang.stiffness );
            
            ang.atoms.x = getAtomType(names[0]);
            ang.atoms.y = getAtomType(names[1]);
            ang.atoms.z = getAtomType(names[2]);
            if( ang.sort() && bWarnFlip ){ printf("WARRNING: angleType[%i](%s) is flipped in %s\n", i, line, fname ); };
            sprintf( name, "%s-%s-%s", atypes[ang.atoms.x].name , atypes[ang.atoms.y].name, atypes[ang.atoms.z].name );
            //printf( "loadAngleTypes[%i](%s) ang0 %g k %g \n", angles.size(), name, ang.angle0, ang.stiffness );
            angles.push_back(ang);
            if( !angleDict.insert({ name, angles.size()-1} ).second ){ printf("WARRNING: angleType[%i](%s) is duplicated !!! => Ignore \n", i, line ); exit(0); };

            // uint64_t id = ang.getId();
            // auto found = angles_.find(id);
            // if( found != angles.end() ){ printf( "WARRNIMG!!! loadAgnleTypes() angleType[%i] same as ", i); printAngle(found->second); };
            //angles_[id]=ang;
        }
        return i;
    }

    int loadDihedralTypes(const char * fname, bool exitIfFail=true, bool bWarnFlip=true){
        printf( "MMFFparams::loadDihedralTypes(%s)\n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            if(exitIfFail)exit(0);
            return -1;
        }
        char buff[1024];
        char names[4][8];
        char name[64];
        char * line;
        DihedralType dih;
        int i=0;
        for( i; i<1000; i++){
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            if(line[0]=='#') continue;
            //printf( "line(%s)\n", line );
            sscanf(  line, "%s %s %s %s %i %lf %lf %i\n", names[0], names[1], names[2], names[3], &dih.bo, &dih.k, &dih.ang0, &dih.n );
            dih.atoms.x = getAtomType(names[0]);
            dih.atoms.y = getAtomType(names[1]);
            dih.atoms.z = getAtomType(names[2]);
            dih.atoms.w = getAtomType(names[3]);
            if( dih.sort() && bWarnFlip ){ printf("WARRNING: dihedralType[%i](%s) is flipped in %s\n", i, line, fname ); }; 
            sprintf( name, "%s-%s-%s-%s-%i", atypes[dih.atoms.x].name , atypes[dih.atoms.y].name, atypes[dih.atoms.z].name, atypes[dih.atoms.w].name, dih.bo );
            //printf( "dihedral[%i](%s) %s-%s-%s-%s-%i k %g ang0 %g n %i \n", dihedrals.size(), atypes[dih.atoms.x].name , atypes[dih.atoms.y].name, atypes[dih.atoms.z].name, atypes[dih.atoms.w].name, dih.bo, dih.k, dih.ang0, dih.n  );
            //printf( "dihedral[%i](%s) k %g ang0 %g n %i \n", dihedrals.size(), name, dih.k, dih.ang0, dih.n  );
            dihedrals.push_back(dih);
            if( !dihedralDict.insert({ name, dihedrals.size()-1} ).second ){ printf("WARRNING: dihedralType[%i](%s) is duplicated !!! => Ignore \n", i, line ); exit(0); };
        }
        return i;
    }

    // TBD - shall we get rid of subtypes...?
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

    // TBD - move dicts into load
    void makeIdDicts(){
        bonds_    .clear(); for(int i=0; i<bonds    .size(); i++){ bonds_    .insert( { bonds    [i].id(), &bonds    [i] } ); }
        angles_   .clear(); for(int i=0; i<angles   .size(); i++){ angles_   .insert( { angles   [i].id(), &angles   [i] } ); }
        dihedrals_.clear(); for(int i=0; i<dihedrals.size(); i++){ dihedrals_.insert( { dihedrals[i].id(), &dihedrals[i] } ); }
    }

    void init(const char* fElementTypes=0, const char* fAtomTypes=0, const char* fBondTypes=0, const char* fAngleTypes=0, const char* fDihedralTypes=0){
        //printf( "MMFFparams::init sElementTypes='%s' sAtomTypes='%s' sBondTypes='%s' sAngleTypes='%s' sDihedralTypes='%s' \n", fElementTypes, fAtomTypes, fBondTypes, fAngleTypes, fDihedralTypes ); exit(0);
        if(fElementTypes )loadElementTypes ( fElementTypes  );
        if(fAtomTypes    ){
            loadAtomTypes( fAtomTypes );
            assignAllSubTypes();
        }
        exit(0);
        // FIN QUA...
        if(fBondTypes    )loadBondTypes    ( fBondTypes     );
        if(fAngleTypes   )loadAngleTypes   ( fAngleTypes    );
        if(fDihedralTypes)loadDihedralTypes( fDihedralTypes );
        makeIdDicts();
    }
    // dictionaries
    void initDefaultAtomTypeDict(){
        makeDefaultAtomTypeDict( atomTypeNames, atomTypeDict );
    }
    ////////////////////////////////////
    // END OF LOADING TABLES OF TYPES //
    ////////////////////////////////////

    //////////////////////////////////////////////////////////
    // ASSIGN ELEMENT, ATOM, BOND, ANGLE AND DIHEDRAL TYPES //
    //////////////////////////////////////////////////////////
    // TBD return pointer?
    int getElementType(const char* s, bool bErr=true)const{
        //printf( "getAtomType(%s) bErr=%i \n", s, bErr );
        auto found = elementTypeDict.find(s);
        if(found==elementTypeDict.end()){ 
            if(bErr){ printf( "ERROR: MMFFparams::getElementType(%s) not found !!! => exit() \n", s ); printAtomTypeDict(); exit(0); }
            return -1; 
        }
        return found->second;
    }

    // TBD return pointer?
    int getAtomType(const char* s, bool bErr=true)const{
        //printf( "getAtomType(%s) bErr=%i \n", s, bErr );
        auto found = atomTypeDict.find(s);
        if(found==atomTypeDict.end()){ 
            if(bErr){ printf( "ERROR: MMFFparams::getAtomType(%s) not found !!! => exit() \n", s ); printAtomTypeDict(); exit(0); }
            return -1; 
        }
        return found->second;
    }

    BondType* getBondType( int ityp, int jtyp, int order, bool bParrents=true, bool bElem=true )const{
        //uint64_t id = BondType::getId( atypes[atyp1].iZ, atypes[atyp2].iZ, btyp );
        uint64_t id   = BondType::getId( ityp, jtyp, order );
        auto it       = bonds_.find(id);
        if( it != bonds_.end() ){ return it->second; } 
        if(bParrents==0){
            if(reportIfMissing){ printf("WARNING getBondParams(ityp=%i,jtyp=%i,order=%i) missing, trying find by parents(%i,%i) \n", ityp, jtyp, order,   atypes[ityp].parrent, atypes[jtyp].parrent ); };
            id  = BondType::getId( atypes[ityp].parrent,        jtyp,          order ); it = bonds_.find(id); if(it!=bonds_.end()){ return it->second; } 
            id  = BondType::getId(        ityp,          atypes[jtyp].parrent, order ); it = bonds_.find(id); if(it!=bonds_.end()){ return it->second; } 
            id  = BondType::getId( atypes[ityp].parrent, atypes[jtyp].parrent, order ); it = bonds_.find(id); if(it!=bonds_.end()){ return it->second; } 
        }
        if(bElem){
            if(reportIfMissing){ printf("WARNING getBondParams(ityp=%i,jtyp=%i,order=%i) missing, trying find by elements(%i,%i) \n", ityp, jtyp, order, atypes[ityp].iZ, atypes[jtyp].iZ  ); };
            int i0 = atomTypeDict.find( elementOfAtomType(ityp)->name )->second;
            int j0 = atomTypeDict.find( elementOfAtomType(jtyp)->name )->second; 
            id  = BondType::getId( i0, j0, order ); 
            it = bonds_.find(id); 
            if( it!=bonds_.end()){ return it->second; }
        }
        if(reportIfMissing){ printf("WARNING getBondParams(ityp=%i,jtyp=%i,order=%i) missing => defaults: l0 %g k %g \n", ityp, jtyp, order, default_bond_length, default_bond_stiffness ); };
        if(exitIfMissing){ printf("=> exit(0)\n");exit(0); };  
        return 0;
    }

    AngleType* getAngleType( int iat, int ityp, int jat, bool bWildcards=true, bool bParrents=true ){
        char tmp[64];
        if(iat>jat){ _swap(iat,jat); }

        sprintf( tmp, "%s-%s-%s", atypes[iat].name,  atypes[ityp].name,  atypes[jat].name );
        if(echoTry)printf("try angleDict[%s]\n", tmp);
        auto found = angleDict.find(tmp);
        if( found != angleDict.end() ){
            printf( "found angleDict[%s]\n", tmp );
            return &angles[found->second];
        }
        
        if(bWildcards){
            if(reportIfMissing){ printf("WARRNING getAngleType(%s-%s-%s) missing, trying find by wildcard(*-%s-*) \n", atypes[iat].name,atypes[ityp].name,atypes[jat].name,  atypes[ityp].name   ); };
            // wildcard
            sprintf( tmp, "*-%s-*", atypes[ityp].name );
            if(echoTry)printf("try angleDict[%s]\n", tmp);
            found = angleDict.find(tmp);
            if( found != angleDict.end() ){
                printf( "found angleDict[%s]\n", tmp );
                return &angles[found->second];
            }
        }

        if(bParrents){
            if(reportIfMissing){ printf("WARRNING getAngleType(%s-%s-%s) missing, trying find by parrents(%s-%s-%s) \n", atypes[iat].name,atypes[ityp].name,atypes[jat].name,   atypes[atypes[iat].parrent].name, atypes[atypes[ityp].parrent].name, atypes[atypes[jat].parrent].name ); };
            int i1,i2,i3;
            for(int i=0; i<8; i++ ){
                if(i&1){ i1=atypes[ityp].parrent; }else{ i1=ityp; }
                if(i&2){ i2=atypes[iat ].parrent; }else{ i2=iat;  }
                if(i&4){ i3=atypes[jat ].parrent; }else{ i3=jat;  }
                sprintf( tmp, "%s-%s-%s", atypes[i2].name, atypes[i1].name, atypes[i3].name );
                if(echoTry)printf("try angleDict[%s]\n", tmp);
                found = angleDict.find(tmp);
                if( found != angleDict.end() ){
                    printf( "found[%i] angleDict[%s]\n", i, tmp );
                    return &angles[found->second];
                }
            }
            if(bWildcards){
                if(reportIfMissing){ printf("WARRNING getAngleType(%s-%s-%s) missing, trying find by wildcard&parents(*-%s-*) \n", atypes[iat].name,atypes[ityp].name,atypes[jat].name,    atypes[atypes[ityp].parrent].name ); };
                int i1;
                for(int i=0; i<2; i++ ){
                    if(i&1){ i1=atypes[ityp].parrent; }else{ i1=ityp; }
                    sprintf( tmp, "*-%s-*", atypes[i1].name );
                    if(echoTry)printf("try angleDict[%s]\n", tmp);
                    found = angleDict.find(tmp);
                    if( found != angleDict.end() ){
                        printf( "found[%i] angleDict[%s]\n", i, tmp );
                        return &angles[found->second];
                    }
                }   
            }
        }

        return 0;
    }

    DihedralType* getDihedralType( int iat, int ityp, int jtyp, int jat, int order, bool bWildcards=true, bool bParrents=true ){
        char tmp[64];
        if(ityp>jtyp){ _swap(ityp,jtyp); _swap(iat,jat); }
        sprintf( tmp, "%s-%s-%s-%s-%i", atypes[iat].name,atypes[ityp].name,atypes[jtyp].name,atypes[jat].name, order );
        if(echoTry)printf("try dihedralDict[%s]\n", tmp);
        auto found = dihedralDict.find(tmp);
        if( found != dihedralDict.end() ){ 
            printf( "found dihedralDict[%s]\n", tmp );
            return &dihedrals[found->second];
        }
        if(bWildcards){
            if(reportIfMissing){ printf("WARRNING getDihedralType(%s-%s-%s-%s-%i) missing, trying find by wildcard(*-%s-%s-*-%i) \n", atypes[iat].name,atypes[ityp].name,atypes[jtyp].name,atypes[jat].name,order,  atypes[ityp].name,atypes[jtyp].name,order  ); };
            sprintf( tmp, "*-%s-%s-*-%i", atypes[ityp].name,atypes[jtyp].name,  order );
            if(echoTry)printf("try dihedralDict[%s]\n", tmp);
            found = dihedralDict.find(tmp);
            if( found != dihedralDict.end() ){ 
                printf( "found dihedralDict[%s]\n", tmp );
                return &dihedrals[found->second];
            }
        }
        if(bParrents){
            if(reportIfMissing){ printf("WARRNING getDihedralType(%s-%s-%s-%s-%i) missing, trying find by parrents(*-%s-%s-*-%i) \n", atypes[iat].name,atypes[ityp].name,atypes[jtyp].name,atypes[jat].name,order,  atypes[ityp].name,atypes[jtyp].name,order   ); };
            int i1,i2,i3,i4;
            for(int i=0; i<16; i++ ){
                if(i&1){ i1=atypes[iat ].parrent; }else{ i1=iat;  }
                if(i&4){ i2=atypes[ityp].parrent; }else{ i2=ityp; }
                if(i&8){ i3=atypes[jtyp].parrent; }else{ i3=jtyp; }
                if(i&2){ i4=atypes[jat ].parrent; }else{ i4=jat;  }
                if(i2>i3){ _swap(i2,i3); _swap(i1,i4); }
                sprintf( tmp, "%s-%s-%s-%s-%i", atypes[i1].name,atypes[i2].name,atypes[i3].name,atypes[i4].name,  order );
                if(echoTry)printf("try dihedralDict[%s]\n", tmp);
                found = dihedralDict.find(tmp);
                if( found != dihedralDict.end() ){ 
                    printf( "found dihedralDict[%s]\n", tmp );
                    return &dihedrals[found->second];
                }
            }
            if(bWildcards){
                if(reportIfMissing){ printf("WARRNING getDihedralType(%s-%s-%s-%s-%i) missing, trying find by wildcard&parrents(*-%s-%s-*-%i) \n", atypes[iat].name,atypes[ityp].name,atypes[jtyp].name,atypes[jat].name,order,  atypes[ityp].name,atypes[jtyp].name,order   ); };
                int i1,i2;
                for(int i=0; i<16; i++ ){
                    if(i&1){ i1=atypes[ityp].parrent; }else{ i1=ityp; }
                    if(i&2){ i2=atypes[jtyp].parrent; }else{ i2=jtyp; }
                    if(i1>i2){ _swap(i1,i2); }
                    sprintf( tmp, "*-%s-%s-*-%i", atypes[i1].name,atypes[i2].name,  order );
                    if(echoTry)printf("try dihedralDict[%s]\n", tmp);
                    found = dihedralDict.find(tmp);
                    if( found != dihedralDict.end() ){
                        printf( "found angleDict[%s]\n", tmp );
                        return &dihedrals[found->second];
                    }
                }   
            }
        }
        return 0;
    }
    ////////////////////////////
    // END OF TYPE ASSIGNMENT //
    ////////////////////////////

    /////////////////////////////
    // ASSIGNING FF PARAMETERS //
    /////////////////////////////
    // TBD should we do it here?
    // assigning bond parameters
    bool getBondParams( int ityp, int jtyp, int order,  double& l0, double& k, bool bParrents=true, bool bElem=true )const{
        BondType* bp = getBondType( ityp, jtyp, order, bParrents, bElem );
        if( bp==0 ){ l0=bp->length; k=bp->stiffness; return false; }else{ l0=bp->length; k=bp->stiffness; return true; }
    }

    void fillBondParams( int nbonds, Vec2i * bond2atom, int * bondOrder, int * atomType, double * bond_0, double * bond_k ){
        //printf("fillBondParams: %i\n", nbonds);
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            getBondParams( atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
            //printf( "%i (%i %i) %i %g %g \n", i, atomType[ib.x], atomType[ib.y], bondOrder[i], bond_0[i], bond_k[i] );
        }
    }

    // assigning angle parameters UFF
    // TBD it returns the force constant only, not the equilibrium angle...
    double assignAngleParamUFF( int ic, int ia, int ib, double ra, double rb )const{  
        const AtomType& tc    = atypes[ic];
        const ElementType* ei = elementOfAtomType(ia);
        const ElementType* ej = elementOfAtomType(ib);
        double Qi = ei->Quff;
        double Qj = ej->Quff;
        double ang0 = tc.Ass*deg2rad;
        double c0 = cos(ang0);
        double s0 = sin(ang0);
        double r = sqrt(  ra*ra + rb*rb - 2*ra*rb*c0 ); // equlibirum distance (bond lenght) between peripheral atoms
        double K =  28.79898 * Qi*Qj * ( 3. *ra*rb*s0*s0 - r*r*c0  ) /( r*r*r*r*r );
        return K;
    }

    // vector with all non-bonded parameters
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

    // assign QEq parameters
    // TBD only electronegativity and hardness are considered...
    void assignQEq( int n, int* itypes, double* affins, double* hards )const{
        for(int i=0; i<n; i++){
            int ityp = itypes[i];
            int iet  = atypes[ityp].element; 
            affins[i]=etypes[iet].Eaff;
            hards [i]=etypes[iet].Ehard;
        }
    }
    ///////////////////////
    // END FF PARAMETERS //
    ///////////////////////

    /////////////////
    // PRINT STUFF //
    /////////////////
    void printAtomTypeDict()const{
        for(int i=0; i<atomTypeNames.size(); i++){ printf( "AtomType[%i] %s %i\n", i, atypes[i].name, atomTypeDict.find(atypes[i].name)->second );  };
    }
    
    void printElementTypeDict()const{
        for(int i=0; i<atomTypeNames.size(); i++){ printf( "ElementType[%i] %s %i\n", i,  etypes[i].name, elementTypeDict.find(etypes[i].name)->second );  };
    }
    void printBond(int i)const{
        const BondType& t = bonds[i];
        printf( "bondType[%3i] %s-%s l0(%7.3f) k(%7.3f)\n", i, atypes[t.atoms.x].name, atypes[t.atoms.y].name, t.length, t.stiffness );
    }
    void printAngle(int i)const{
        const AngleType& t = angles[i];
        printf( "angleType[%3i] %s-%s-%s ang0(%7.3f) k(%7.3f) \n", i, atypes[t.atoms.x].name, atypes[t.atoms.y].name, atypes[t.atoms.z].name, t.angle0, t.stiffness );
    }
    void printDihedral(int i)const{
        const DihedralType& t = dihedrals[i];
        printf( "dihedralType[%3i] %s-%s-%s-%s ang0(%7.3f) k(%7.3f) n(%i)\n", i, atypes[t.atoms.x].name, atypes[t.atoms.y].name, atypes[t.atoms.z].name, atypes[t.atoms.w].name, t.ang0, t.k, t.n );
    }
    void printAtomTypes(bool bParams)const{   for(int i=0; i<atypes.size(); i++ ){ atypes[i].print(i, bParams );  }  }
    void printBondTypes    ()const{ printf("MMFFparams::printBondTypes()\n");     for(int i=0; i<bonds.size();     i++ ){ printBond(i);     } }
    void printAngleTypes   ()const{ printf("MMFFparams::printAngleTypes()\n");    for(int i=0; i<angles.size();    i++ ){ printAngle(i);    } }
    void printDihedralTypes()const{ printf("MMFFparams::printDihedralTypes()\n"); for(int i=0; i<dihedrals.size(); i++ ){ printDihedral(i); } }
    void printAngleTypesDict()const{ printf("MMFFparams::printAngleTypesDict()\n");    for( const auto& it : angleDict ){ printf("angle(%s)[%i]\n", it.first.c_str(), it.second ); } }
    //////////////////
    // END OF PRINT //
    //////////////////

    ///////////////////
    // DEAL WITH XYZ //
    ///////////////////
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
            //printf( "atom[%i] name(%s) it=%i \n", i, at_name, it->second );
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

    void writeXYZ( FILE* pfile, int n, const int* atyps, const Vec3d* apos, const char* comment="#comment", const Quat4d* REQs=0, bool just_Element=true, int npi=0, Vec3i nPBC=Vec3i{1,1,1}, Mat3d lvec=Mat3dIdentity ){
        //printf( "MMFFparams::writeXYZ() n=%i REQs=%li just_Element=%i\n", n, (long)REQs, just_Element );
        int npbc = nPBC.totprod();
        fprintf(pfile, "%i\n", (n+npi)*npbc );
        // TBD print lattice vectors
        fprintf(pfile, "%s \n", comment );
        //printf( "MM::Params::writeXYZ() nPBC={%i,%i,%i}\n", nPBC.x,nPBC.y,nPBC.z );
        for(int ic=0;ic<nPBC.z;ic++){ for(int ib=0;ib<nPBC.y;ib++){ for(int ia=0;ia<nPBC.x;ia++){
            Vec3d shift = lvec.a*ia + lvec.b*ib + lvec.c*ic; 
            //printf( "MM::Params::writeXYZ() iPBC(%i,%i,%i) shift(%g,%g,%g)\n", ia,ib,ic, shift.x,shift.y,shift.z );
            //// //////////////------------------- 
        for(int i=0; i<n; i++){
            //printf( "DEBUG writeXYZ()[%i] \n", i );
            int ityp   = atyps[i];
            const Vec3d&  pi = apos[i] + shift;
            const char* symbol; 
            bool byName = true;
            if(just_Element){ 
                byName = false;
                symbol =  etypes[ atypes[ityp].element ].name;
            }
            if(byName){ symbol = atomTypeNames[ityp].c_str(); }
            //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, atypes[ityp].name );
            if(REQs){ fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f     %10.6f\n", symbol, pi.x,pi.y,pi.z, REQs[i].z ); }
            else    { fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n"          , symbol, pi.x,pi.y,pi.z            ); }
        }
        }}}
        for(int i=0; i<npi; i++){
            const Vec3d&  pi = apos[i] + apos[i+n];
            fprintf( pfile, "Pi   %15.10f   %15.10f   %15.10f \n", pi.x,pi.y,pi.z            );

        }
    }

    int saveXYZ( const char * fname, int n, const int* atyps, const Vec3d* apos, const char* comment="#comment", const Quat4d* REQs=0, const char* mode="w", bool just_Element=true, Vec3i nPBC=Vec3i{1,1,1}, Mat3d lvec=Mat3dIdentity ){
        //printf( "MMFFparams::saveXYZ(%s) \n", fname );
        FILE* pfile = fopen(fname, mode );
        if( pfile == NULL ) return -1;
        writeXYZ( pfile, n, atyps, apos, comment, REQs, just_Element, 0, nPBC, lvec );
        fclose(pfile);
        return n;
    }
    ////////////////
    // END OF XYZ //
    ////////////////

    //////////
    // MISC //
    //////////
    // TBD comment a/o remove
    inline const ElementType* elementOfAtomType( int it )const{ return &etypes[atypes[it].element]; }

    // following the graph and getting the ancestor
    const AtomType* getRootParrent(const AtomType* t, int nrecur=0)const{
        if(nrecur>10){ printf("ERROR in MMFFparams.getRootParrent() rootParrent of type(%s) not found in %i recursions => Exit() \n", t->name, nrecur ); exit(0); }
        if( t->parrent==0 ) return t;
        if( (t->parrent<0)||(t->parrent>=atypes.size()) ){ printf("ERROR in MMFFparams.getRootParrent() type(%s).parrent==%i => Exit() \n", t->name,t->parrent ); exit(0); }
        const AtomType* par = &atypes[t->parrent];
        return getRootParrent(par,nrecur+1); // recursion
    }

    void clear( bool bShring=false ){
        etypes.clear();
        atypes.clear();
        atomTypeNames.clear();
        elementTypeDict.clear();
        atomTypeDict.clear();

        bonds.clear();
        bonds_.clear();
        bondDict.clear();

        angles.clear();
        angles_.clear();
        angleDict.clear();

        dihedrals.clear();
        dihedrals_.clear();
        dihedralDict.clear();
        if(bShring){
            etypes.shrink_to_fit();
            atypes.shrink_to_fit();
            atomTypeNames.shrink_to_fit();
            elementTypeDict.rehash(0);
            atomTypeDict.rehash(0);

            bonds.shrink_to_fit();
            bonds_.rehash(0);
            bondDict.rehash(0);

            angles.shrink_to_fit();
            angles_.rehash(0);
            angleDict.rehash(0);

            dihedrals.shrink_to_fit();
            dihedrals_.rehash(0);
            dihedralDict.rehash(0);
        }
    }
    //////////////
    // END MISC //
    //////////////

};

#endif
