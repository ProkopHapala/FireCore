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


class AtomType{ public:
    char      name[8];
    uint8_t   iZ;         // proton number
    uint8_t   neval;      // number of valence electrons
    uint8_t   valence;    // sum of bond orders of all bonds
    uint8_t   sym;        // sp3=0,sp2=1, sp1=2,  ... tetrahedral, triangle, linear, kink, square, octahedral
    uint32_t  color;
    double    RvdW;
    double    EvdW;
    Vec3i     subTypes=Vec3iZero;  // sp1 sp2 sp3    // Q1 Q2 Q3 (polarized)

    // ---- MMFF
    bool bMMFF;
    double  Ass,Asp, Kss,Ksp,Kep,Kpp;

    // ---- Epairs
    int ne=0;                // number of electron pairs
    double  eRvdW,eEvdW,eQ;  // electron pair REQ parameters
 
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

    //void print(int i){ printf( "AtomType[%i] %s (%i,%i,%i,%i) vdW(R=%lf,E=%lf) clr %x Eaff,Ehard (%g,%g) \n", i, name,  iZ,   neval,  valence,   sym,    RvdW, EvdW,   color, Eaff,Ehard ); }
    void print(int i){ printf( "AtomType[%i] %s (%i,%i,%i,%i) LJ(%g,E=%g) QEq(%g,%g) MMMFF(%g,%g|%g,%g,%g,%g) Epair(%i,%g,%g,%g)\n", i, name,  iZ,   neval,valence,sym,    RvdW,EvdW,   Eaff,Ehard,   Ass,Asp,Kss,Ksp,Kep,Kpp,   ne,eRvdW,eEvdW,eQ  ); }

    inline uint8_t nepair(){ return (neval-valence)/2; };
    inline uint8_t npi   (){ return sym; };

    void fromString( char * str ){
        int iZ_, neval_, valence_, sym_;
        //char sclr[6];           1   2       3         4        5         6     7      8        9        10   11     12     13   14   15    16   17   18   19  20   21    22
        int nret = sscanf( str, " %s    %i    %i       %i        %i      %lf    %lf    %x        %lf    %lf     %lf   %lf    %lf  %lf   %lf   %lf  %lf   %lf  %i   %lf   %lf    %lf ", 
                                 name, &iZ_, &neval_, &valence_, &sym_,  &RvdW, &EvdW, &color,   &Eaff, &Ehard, &Ra, &eta,   &Ass,&Asp, &Kss,&Ksp,&Kep,&Kpp, &ne,&eRvdW,&eEvdW,&eQ );
        iZ=iZ_; neval=neval_; valence=valence_; sym=sym_;
        if(nret<10){ bQEq  = false; Eaff=0;Ehard=0;Ra=0;eta=0;           }else{ bQEq  = true; }
        if(nret<18){ bMMFF = false; Ass=0;Asp=0;Kss=0;Ksp=0;Kep=0;Kpp=0; }else{ bMMFF = true; }
        if(nret<22){ ne=-1,eRvdW=0,eEvdW=0,eQ=0; }
        subTypes=Vec3iZero;
        //printf( "AtomType: %s iZ %i ne %i nb %i sym %i RE(%g,%g) %x \n", name,  iZ,   neval_,  valence,   sym,    RvdW, EvdW,   color );
        //char ss[256]; printf("%s\n", toString(ss) );
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

    int getAtomType(const char* s, bool bErr=true){
        auto found = atomTypeDict.find(s);
        if(found==atomTypeDict.end()){ 
            if(bErr){ printf( "ERROR: MMFFparams::getAtomType(%s) not found !!! => exit() \n", s ); printAtomTypeDict(); exit(0); }
            return -1; 
        }
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
            if(line==NULL)  break;
            //printf( "loadAtomTypes[%i] line=%s", i, line );
            if(line[0]=='#')continue;
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

    inline void assignSubTypes( AtomType& t ){
        //printf( "assignSubTypes %s(iZ=%i)\n", t.name, t.iZ );
        char tmp_name[8];
        const char* ssub[3]{"sp3","sp2","sp1"};
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

    inline void assignRE( int ityp, Vec3d& REQ, bool bSqrtE=false )const{
        REQ.x    = atypes[ityp].RvdW;
        double e = atypes[ityp].EvdW;
        if(bSqrtE) e=sqrt(e);
        REQ.y = e;
    }

    void assignREs( int n, int * itypes, Vec3d * REQs, bool bSqrtE=false, bool bQ0=false )const{
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
        //if(verbosity>0) 
        printf("MMFFparams::init(%s,%s,%s)\n", fatomtypes, fbondtypes, fagnletypes );
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

    int loadXYZ(const char* fname, int& natoms, Vec3d** apos_, Vec3d** REQs_=0, int** atype_=0, int** npis_=0, Mat3d* lvec=0, int verbosity=0 )const{
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
        Vec3d* REQs  =_allocPointer( REQs_,  natoms );
        int*   atype =_allocPointer( atype_, natoms );
        int*   npis  =_allocPointer( npis_,  natoms );
        line = fgets( buff, nbuf, pFile ); // comment
        int ret=0;
        if(lvec){ if( cellFromString( buff, *lvec ) ){ ret=1; }else{ printf("WARRNING: lvec not read from %s \n", fname ); } }
        double Q;
        for(int i=0; i<natoms; i++){
            char at_name[8];
            double junk; 
            int npi;
            line = fgets( buff, nbuf, pFile );  //printf("%s",line);
            int nret = sscanf( line, "%s %lf %lf %lf %lf \n", at_name, &apos[i].x, &apos[i].y, &apos[i].z, &Q, &npi );
            if( nret < 5 ){ Q=0; };
            if( nret < 6 ){ npi=-1; };
            if(npis){ npis[i] =npi; };
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
        return ret;
    }

    void writeXYZ( FILE* pfile, int n, const int* atyps, const Vec3d* apos, const char* comment="#comment", const Vec3d* REQs=0, bool just_Element=true ){
        //printf( "MMFFparams::writeXYZ() n=%i REQs=%li \n", n, (long)REQs );
        fprintf(pfile, "%i\n", n );
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
        };
    }

    int saveXYZ( const char * fname, int n, const int* atyps, const Vec3d* apos, const char* comment="#comment", const Vec3d* REQs=0, const char* mode="w", bool just_Element=true ){
        FILE* pfile = fopen(fname, mode );
        //printf( "saveXYZ(%s) \n", fname );
        if( pfile == NULL ) return -1;
        writeXYZ( pfile, n, atyps, apos, comment, REQs, just_Element );
        fclose(pfile);
        return n;
    }

};

#endif
