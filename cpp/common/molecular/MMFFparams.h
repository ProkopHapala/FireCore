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

class AtomType{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint8_t   sym;        // sp, sp2, sp3 ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;
    double    RvdW;
    double    EvdW;
    // ==== additional
    double    piStiff;
    double    electroneg;
    double    polarizability;

    char* toString( char * str ){
        sprintf( str, "%s %i %i %i %i %lf %lf %x", name,  iZ,   neval,  valence,   sym,    RvdW, EvdW,   color );
        return str;
    }

    inline uint8_t nepair(){ return (neval-valence)/2; };

    void fromString( char * str ){
        int iZ_, neval_, valence_, sym_;
        //char sclr[6];
        sscanf( str, " %s %i %i %i %i %lf %lf %x \n", name, &iZ_, &neval_, &valence_, &sym_,  &RvdW, &EvdW, &color );
        iZ=iZ_; neval=neval_; valence=valence_; sym=sym_;
        //printf( "AtomType: %s iZ %i ne %i nb %i sym %i RE(%g,%g) %x \n", name,  iZ,   neval_,  valence,   sym,    RvdW, EvdW,   color );
        //char ss[256]; printf("%s\n", toString(ss) );
    }

};

// ======================
// ====   MMFFparams
// ======================

class MMFFparams{ public:

    // http://www.science.uwaterloo.ca/~cchieh/cact/c120/bondel.html
    std::vector       <AtomType>           atypes;
    std::vector       <std::string    >    atomTypeNames;
    std::unordered_map<std::string,int>    atomTypeDict;
    std::unordered_map<uint64_t,BondType>  bonds;
    std::unordered_map<uint64_t,AngleType> angles;

    double default_bond_length      = 2.0;
    double default_bond_stiffness   = 1.0;

    bool reportIfMissing=true;
    bool exitIfMissing  =true;

    void initDefaultAtomTypeDict(){
        makeDefaultAtomTypeDict( atomTypeNames, atomTypeDict );
    }

    void printAtomTypeDict(){
        for(int i=0; i<atomTypeNames.size(); i++){ printf( "AtomType[%i] %s %i\n", i, atomTypeNames[i].c_str(), atomTypeDict[atomTypeNames[i]] ); };
    }

    int getAtomType(const char* s){
        auto found = atomTypeDict.find(s);
        if(found==atomTypeDict.end()){ return -1; }
        return found->second;
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
        for(i=0; i<0xFFFF; i++){
        //for(int i; i<0xFFFF; i++){
            //printf( "loadAtomTypes %i \n", i );
            line = fgets( buff, 1024, pFile );
            if(line==NULL) break;
            atyp.fromString( line );
            atypes.push_back(atyp);
            atomTypeNames.push_back( atyp.name );
            atomTypeDict[atyp.name] = atypes.size()-1;

            //char str[1000];
            //atyp.toString( str );
            //printf( "%i %s %i %s \n", i, atyp.name, atypNames[atyp.name], str );
        }
        return i;
    }

    inline void assignRE( int ityp, Vec3d& REQ )const{
        REQ.x = atypes[ityp].RvdW;
        REQ.y = atypes[ityp].EvdW;
    }

    void assignREs( int n, int * itypes, Vec3d * REQs )const{
        //printf( "assignREs %i   %i %i %i \n", n,  itypes, REQs, atypes );
        for(int i=0; i<n; i++){
            //mmff->aLJq [i]  = atoms[i].type;
            //int ityp = atoms[i].type;
            //printf( "i %i \n", i );
            //int ityp = itypes[i];
            //printf( "ityp %i \n", ityp );
            //REQs[i].x = atypes[ityp].RvdW;
            //REQs[i].y = atypes[ityp].EvdW;
            assignRE( itypes[i], REQs[i] );
            //printf( "assignREs i %i ityp %i RE  %g %g  \n", i, ityp, REQs.x, REQs.y );
            //REQs.z = 0;
            //atomTypes[i]  = atoms[i].type;
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
            if(reportIfMissing){ printf("WARRNING!!! getBondParams() missing atype(%i,%i)iZs(%i,%i)o(%i)id(%i)=>l0 %g k %g \n", atyp1, atyp2, atypes[atyp1].iZ, atypes[atyp2].iZ, btyp, id, default_bond_length, default_bond_stiffness ); };
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

};

#endif
