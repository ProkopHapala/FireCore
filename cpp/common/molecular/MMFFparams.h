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
    uint8_t   sym;        // sp3=0,sp2=1, sp1=2,  ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;
    double    RvdW;
    double    EvdW;

    // ---- charge equlibration
    bool   bQEq;
    double Eaff,Ehard,Ra,eta;
    // ==== additional
    double    piStiff;
    //double    electroneg;
    //double    polarizability;

    char* toString( char * str ){
        sprintf( str, "%s %i %i %i %i %lf %lf %x", name,  iZ,   neval,  valence,   sym,    RvdW, EvdW,   color );
        return str;
    }

    void print(int i){ printf( "AtomType[%i] %s (%i,%i,%i,%i) vdW(R=%lf,E=%lf) clr %x Eaff,Ehard (%g,%g) \n", i, name,  iZ,   neval,  valence,   sym,    RvdW, EvdW,   color, Eaff,Ehard ); }

    inline uint8_t nepair(){ return (neval-valence)/2; };
    inline uint8_t npi   (){ return sym; };

    void fromString( char * str ){
        int iZ_, neval_, valence_, sym_;
        //char sclr[6];                                                            1   2      3         4        5         6     7      8        9        10   11     12
        int nret = sscanf( str, " %s %i %i %i %i %lf %lf %x %lf %lf %lf %lf\n", name, &iZ_, &neval_, &valence_, &sym_,  &RvdW, &EvdW, &color,   &Eaff, &Ehard, &Ra, &eta );
        iZ=iZ_; neval=neval_; valence=valence_; sym=sym_;
        if(nret<10){ bQEq = false; }else{ bQEq=true; }
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
    //Vec3d  default_REQ            = {1.487, 0.0006808, 0.0};  // Hydrogen
    Vec3d  default_REQ              = {1.500, 0.0005000, 0.0};  // Hydrogen like

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
        for(i=0; i<10000; i++){
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

    inline void assignRE( int ityp, Vec3d& REQ, bool bSqrtE=false )const{
        REQ.x = atypes[ityp].RvdW;
        double e=atypes[ityp].EvdW;
        if(bSqrtE) e=sqrt(e);
        REQ.y = e;
    }

    void assignREs( int n, int * itypes, Vec3d * REQs, bool bSqrtE=false, bool bQ0=false )const{
        for(int i=0; i<n; i++){
            //printf( " assignREs[%i] %i \n", i, itypes[i] );
            assignRE( itypes[i], REQs[i], bSqrtE );
            if(bQ0) REQs[i].z=0;
        }
    }

    void assignQEq( int n, int* itypes, double* affins, double* hards ){
        for(int i=0; i<n; i++){
            int ityp = itypes[i];
            affins[i]=atypes[ityp].Eaff;
            hards [i]=atypes[ityp].Ehard;
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

    void printAtomTypes(){
        for(int i=0; i<atypes.size(); i++ ){ atypes[i].print(i); }
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

    void init(const char* fatomtypes=0, const char* fbondtypes=0, const char* fagnletypes=0){
        if(fatomtypes )loadAtomTypes( fatomtypes );
        if(fbondtypes )loadBondTypes( fbondtypes );
        if(fagnletypes)loadAgnleType( fagnletypes );
    }


    int loadXYZ(const char* fname, int& natoms, Vec3d** apos_, Vec3d** REQs_=0, int** atype_=0, int verbosity=0 )const{
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
        Vec3d* apos  =_allocPointer( apos_, natoms );
        Vec3d* REQs  =_allocPointer( REQs_, natoms );
        int*   atype =_allocPointer( atype_, natoms );
        line = fgets( buff, nbuf, pFile ); // comment
        double Q;
        for(int i=0; i<natoms; i++){
            char at_name[8];
            double junk;
            line = fgets( buff, nbuf, pFile );  //printf("%s",line);
            int nret = sscanf( line, "%s %lf %lf %lf %lf \n", at_name, &apos[i].x, &apos[i].y, &apos[i].z, &Q );
            if( nret < 5 ){ Q=0; };

            auto it = atomTypeDict.find( at_name );
            if( it != atomTypeDict.end() ){
                int ityp=it->second;
                if(atype_)atype[i] = ityp;
                if(REQs_){ assignRE( ityp, REQs[i], true ); REQs[i].z=Q; }
            }else{
                if(atype_)atype[i] = -1;
                if(REQs_)REQs[i]  = default_REQ;
            }
        }
        fclose(pFile);
        return natoms;
    }

    void writeXYZ( FILE* pfile, int n, const int* atypes, const Vec3d* apos, const char* comment="#comment", const Vec3d* REQs=0 ){
        printf( "MMFFparams::writeXYZ() n=%i REQs=%li \n", n, (long)REQs );
        fprintf(pfile, "%i\n", n );
        fprintf(pfile, "%s \n", comment );
        for(int i=0; i<n; i++){
            //printf( "DEBUG writeXYZ()[%i] \n", i );
            int ityp   = atypes[i];
            const Vec3d&  pi = apos[i];
            //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
            if(REQs){ fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f     %10.6f\n", atomTypeNames[ityp].c_str(), pi.x,pi.y,pi.z, REQs[i].z ); }
            else    { fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n"    , atomTypeNames[ityp].c_str(), pi.x,pi.y,pi.z            ); }
        };
    }

    int saveXYZ( const char * fname, int n, const int* atypes, const Vec3d* apos, const char* comment="#comment", const Vec3d* REQs=0 ){
        FILE* pfile = fopen(fname, "w");
        if( pfile == NULL ) return -1;
        writeXYZ( pfile, n, atypes, apos, comment, REQs );
        fclose(pfile);
        return n;
    }

};

#endif
