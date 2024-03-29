//#pragma once

#ifndef Molecule_h
#define Molecule_h

#include <unordered_map>
#include <vector>
#include <string>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
//#include "IOutils.h"

#include "MMFFparams.h"

inline int atomName2int(char ch ){
    int i = -1;
    switch(ch){
        case 'H': i=1;  break;
        case 'B': i=5;  break;
        case 'C': i=6;  break;
        case 'N': i=7;  break;
        case 'O': i=8;  break;
        case 'F': i=9;  break;
        case 'P': i=15; break;
        case 'S': i=16; break;
    }
    return i;
}

void cpstr( const char* str, char* tmp, int i0, int n ){
    tmp[n]='\0';
    for(int i=0; i<n; i++){
        tmp[i]=str[i0+i];
    }
}

double getDouble(const char* str, int i0, int i1 ){
    const int n=i1-i0;
    char  tmp[n+1];
    cpstr( str, tmp, i0, n );
    //double f; //sscanf( "%lf", &f );
    return atof(tmp);
}

double getInt(const char* str, int i0, int i1 ){
    const int n=i1-i0;
    char  tmp[n+1];
    cpstr( str, tmp, i0, n );
    //printf("str: %s",   str );
    //printf("tmp: %s\n", tmp );
    //int i; sscanf( "%i", &i );
    return atoi(tmp);
}


inline int otherAtom( Vec2i ib, int ia ){ return (ia==ib.x)?ib.y:ib.x; }

#define NBMAX 4

class Molecule{ public:
    int  natoms=0, nbonds=0;
    Vec3d  * pos       = NULL;
    Vec2i  * bond2atom = NULL;
    //double * charges   = NULL;
    Quat4d * REQs      = NULL;
    int    * atomType  = NULL;
    int    * bondType  = NULL;
    int    * npis      = NULL;

    int   * atom_nb    = NULL;
    int   * atom2bond  = NULL;

    int nang = 0;
    Vec2i  * ang2bond   = NULL;
    double * ang0s      = NULL;
    //int * ang2atom = NULL;

    const MMFFparams* params = 0;
    // ToDo: this is redudant can be changed to params
    const std::vector       <std::string    >* atomTypeNames=0;
    const std::unordered_map<std::string,int>* atomTypeDict =0;

    void allocate( int natoms_, int nbonds_ ){
        natoms=natoms_; nbonds=nbonds_;
        _realloc( pos       , natoms );
        //_realloc( charges   , natoms );
        _realloc( REQs      , natoms );
        _realloc( atomType  , natoms );
        _realloc( npis      , natoms );
        _realloc( bondType  , nbonds );
        _realloc( bond2atom , nbonds );

    }

    void bondsOfAtoms(){
        atom_nb   = new int[natoms];
        atom2bond = new int[natoms*NBMAX];
        for(int i=0; i<natoms; i++){ atom_nb[i]=0; }
        for(int i=0; i<nbonds; i++){
            Vec2i ib = bond2atom[i];
            int nb1 = atom_nb[ib.x];
            int nb2 = atom_nb[ib.y];
            if( (nb1>=4)||(nb2>=4) ){
                printf( "ERROR in Molecule.bondsOfAtoms() [ib %i] bond count(%i,%i) > NBMAX(%i) \n", i, nb1, nb2, NBMAX );
                exit(0);
            }
            atom2bond[ib.x*NBMAX+nb1] = i;
            atom2bond[ib.y*NBMAX+nb2] = i; // -i to denote directionality ?
            atom_nb[ib.x]=nb1+1;
            atom_nb[ib.y]=nb2+1;
        }
    }

    /*

    void sortAtomsOfBonds(){

    }

    void sortAtomNeighs(int ia){
        int nb = atom_nb[ia];
        int* ibs = atom2bond + i*NBMAX;
        int nghs[NBMAX];
        int bs  [NBMAX];
        for(int i=0; i<nb i++){

        }
        for(int i=0; i<nb i++){
            int ib = selectMinHigher(ja, nbconf, ibs );
        }
    }

    void remakeBondsSorted(){
        if( (atom_nb==0)||(atom2bond==0) )bondsOfAtoms();
        int ib=0;

        for(int i=0; i<natoms; i++){
            int nb = atom_nb[i];
            for(int ij=0; ij<nb; ij++){
                int j = atom2bond[i*NBMAX+ij];
                bond2atom[ib] = (Vec2i){i,j};
                ib++;
            }
        }
    }
    */

#ifdef MMFFparams_h
    void assignREQs( MMFFparams& params ){
        for(int i=0; i<natoms; i++){
            int ityp=atomType[i];
            const AtomType& typ = params.atypes[ityp];
            //printf( "assignREQs[%i] ityp %i RvdW %g REQs %li params.atypes %li \n", i, ityp, typ.RvdW, (long)REQs, (long)(params.atypes[0]) );
            REQs[i].x=typ.RvdW;
            REQs[i].y=typ.EvdW;
        }
    }
#endif

    void findBonds_brute(double Rmax=1.8, const bool bUseREQs=false){
        double R2max = sq(Rmax);
        std::vector<Vec2i> pairs;
        for(int i=0;i<natoms;i++){
            const Vec3d& pi = pos[i];
            double Ri; if( bUseREQs ){ Ri=REQs[i].x; }
            for(int j=i+1;j<natoms;j++){
                Vec3d d; d.set_sub(pos[j],pi);
                double r2 = d.norm2();
                double R2;
                if( bUseREQs ){
                    R2=(Ri+REQs[j].x)*Rmax;
                    R2*=R2;
                    //printf( "bond[%i : %i] r %g R %g factor %g Bonded=%i\n", i, j, sqrt(r2), sqrt(R2), Rmax, r2<R2 ) ;
                }else{
                    R2=R2max;
                }
                //printf( "%i,%i %g %g  (%g,%g,%g)  (%g,%g,%g) \n", i, j, r2, R2max, pi.x,pi.y,pi.z,    pos[j].x,pos[j].y,pos[j].z );
                if(r2<R2){
                    pairs.push_back({i,j});
                }
            }
        }
        nbonds=pairs.size();
        _realloc( bond2atom , nbonds );
        _realloc( bondType  , nbonds );
        //printf("Molecule::findBonds_brute() bond2atom=%li \n", bond2atom);
        for(int i=0;i<nbonds;i++){  bond2atom[i]=pairs[i]; bondType[i]=-1; }
    }

    void bindParams( const MMFFparams* params_ ){
        params = params_;
        atomTypeNames = &params->atomTypeNames;
        atomTypeDict  = &params->atomTypeDict;
    }

    int countAtomType( int ityp )const{
        int n=0;
        for(int i=0; i<natoms; i++){ if(atomType[i]==ityp)n++; }
        return n;
    }

    Mat3d findAxis( Vec2i ifw, Vec2i ilf ){
        Mat3d M;
        M.fromSideUp( pos[ifw.y]-pos[ifw.x], pos[ilf.y]-pos[ilf.x] );
        return M;
    }

    void orient( int i0, Vec2i ifw, Vec2i ilf ){
        Mat3d M;
        Vec3d p0 = pos[i0];
        M.fromSideUp( pos[ifw.y]-pos[ifw.x], pos[ilf.y]-pos[ilf.x] );
        for(int i=0; i<natoms; i++){
            Vec3d p = pos[i]; p.sub( p0 );
            M.dot_to_T( p, pos[i] );
        }
    }

    inline double getAng0( int ia ){
        int ityp = atomType[ia];
        if(params) ityp = params->atypes[ityp].iZ; // ToDo : later this can be improved
        int nb   = atom_nb [ia];
        printf( "getAng0 iT %i iZ %i nb %i \n", atomType[ia], ityp, nb );
        switch(ityp){
            case 1: return M_PI; break; // H
            case 6: // C
                if     (nb==2){ return M_PI;               }
                else if(nb==3){ return M_PI*0.66666666666; }
                else if(nb==4){ return M_PI*0.60833333333; }
                else          { return M_PI;               }
                break;
            case 7: // N
                if     (nb==2){ return M_PI*0.6666666666;  }
                else if(nb==3){ return M_PI*0.60833333333; }
                else          { return M_PI;               }
                break;
            case 5: return M_PI; break; // B
            case 8:  // O
            case 14: // Si
            case 16: // S
                return M_PI*0.60833333333;
            default: return M_PI;
        }
        return M_PI;
    }

    int autoAngles( bool bAng0=false ){
        // create angle between each pair of bonds of an atom
        nang = 0;
        if( (atom_nb==0)||(atom2bond==0) )bondsOfAtoms();
        for(int i=0; i<natoms; i++){ int nb = atom_nb[i]; nang+=(nb*(nb-1))>>1; }
        ang2bond       = new Vec2i [nang];
        if(bAng0)ang0s = new double[nang];
        int * a2b = atom2bond;
        int iang = 0;
        for(int ia=0; ia<natoms; ia++){
            int nb = atom_nb[ia];
            double a0;
            if((bAng0)&&(nb>1)){
                a0 = getAng0( ia );
                printf( "atom[%i] a0 %g \n ", ia, a0 );
            }
            for(int i=0; i<nb; i++){
                int ib1 = a2b[i];
                for(int j=0; j<i; j++){
                    ang2bond[iang]={ib1,a2b[j]};
                    if(bAng0){
                        ang0s[iang] = a0;
                        printf( "autoAngles[%i][ia%i|%i,%i] %g[rad] | nZ %i nb %i \n", iang, ia,i,j, ang0s[iang], atomType[ia], nb );
                    }
                    iang++;
                }
            }
            a2b+= NBMAX;
        }

        return iang;
    }

    Vec3d getCOG_av()const{
        Vec3d cog =  Vec3dZero;
        for(int i=0; i<natoms; i++){
            cog.add(pos[i]);
        }
        cog.mul(1.0/natoms);
        return cog;
    }

    Vec3d getCOG_minmax()const{
        Vec3d p0 = Vec3dmax;
        Vec3d p1 = Vec3dmin;
        for(int i=0; i<natoms; i++){
            p0.setIfLower  (pos[i]);
            p1.setIfGreater(pos[i]);
        }
        //printf("getCOG_minmax pmin(%g,%g,%g) pmax(%g,%g,%g) \n", p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z );
        return (p0+p1)*0.5;
    }

    void addToPos( Vec3d dp ){ for(int i=0; i<natoms; i++){ pos[i].add(dp); } }

    void cogTo0(){ Vec3d cog = getCOG_minmax(); addToPos( cog*-1.0 ); };

    void FindRotation( Mat3d& rot){
        Mat3d XX=Mat3dZero;
        for(int i=0; i<natoms; i++){
            XX.addOuter( pos[i], pos[i], 1.0 );
        }
        //printf( "XX: " ); printMat(XX);
        Vec3d evs;
        XX.eigenvals(evs);
        evs.sort();
        printf(  "FindRotation evs(%g,%g,%g) \n", evs.x,evs.y,evs.z );
        Vec3d c;
        XX.eigenvec( evs.x, rot.a );
        XX.eigenvec( evs.y, rot.b );
        XX.eigenvec( evs.z, rot.c );
    }

    void orient( Vec3d center, Vec3d dir, Vec3d up ){
        Mat3d rot; rot.fromDirUp( dir, up );
        for(int i=0; i<natoms; i++){
            Vec3d p = pos[i]-center;
            rot.dot_to(p,pos[i]);
        }
    }

    void orient( int i0, int i1, int i2 ){
        Vec3d center = pos[i0];
        orient( center, (pos[i1]-center).normalized(), pos[i2]-center );
    }

    /*
    // DEPRECATED
    int loadMol_old( const char* fname ){
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        sscanf( line, "%i %i\n", &natoms, &nbonds );
        //printf("%i %i \n", natoms, nbonds );

        allocate(natoms,nbonds);
        for(int i=0; i<natoms; i++){
            //char ch;
            char at_name[8];
            double junk;
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%lf %lf %lf %s %lf %lf\n", &pos[i].x, &pos[i].y, &pos[i].z,  at_name, &junk, &REQs[i].z );
            //printf(       "%lf %lf %lf %s %lf %lf\n",  pos[i].x,  pos[i].y,  pos[i].z,  at_name,  junk,  REQs[i].z );
            // atomType[i] = atomChar2int( ch );
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                atomType[i] = it->second;
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
        }
        for(int i=0; i<nbonds; i++){
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line, "%i %i %i\n", &bond2atom[i].x, &bond2atom[i].y, &bondType[i] );
            printf(       "%i %i %i\n",  bond2atom[i].x,  bond2atom[i].y,  bondType[i] );
            bond2atom[i].x--;
            bond2atom[i].y--;
        }
        return natoms + nbonds;
    }
    */

    void assignAtomType(int i, const char * at_name ){
        if(atomTypeDict){
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                atomType[i] = it->second;
                //printf( " %s -> %i \n", at_name,  atomType[i] );
                if(1==it->second)REQs[i].x=1.0; // Hydrogen is smaller
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
            //printf( "at_name[%i] %s -> %i \n", i, at_name, atomType[i] );
        }else{
            atomType[i] = atomName2int( at_name[0] );
        }
    }


    int loadMol( const char* fname, int verbosity=0 ){
        // 0        10         20
        //   -13.0110  -15.2500   -0.0030 N   0  0  0  0  0  0  0  0  0  0  0  0
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }

        char buff[1024];
        char * line;
        int nl;
        try {
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        //sscanf( line, "%i %i\n", &natoms, &nbonds );
        //printf("natoms, nbonds %i %i \n", natoms, nbonds );
        //line = fgets( buff, 1024, pFile );
        natoms = getInt( line, 0, 3 );
        nbonds = getInt( line, 3, 6 );
        //printf("loadMol(%s):  natoms, nbonds %i %i \n", fname, natoms, nbonds );
        //exit(0);
        //return 0;
        allocate(natoms,nbonds);
        int istr_x      =0;
        int istr_y      =10;
        int istr_z      =20;
        int istr_name   =25;
        int istr_charge =30;
        for(int i=0; i<natoms; i++){
            //char ch;
            double junk;
            char at_name[8];
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            sscanf( line,               "%lf %lf %lf %s %lf %lf\n", &pos[i].x, &pos[i].y, &pos[i].z,  at_name, &junk, &(REQs[i].z) );
            if(verbosity>1)printf( ".mol %lf %lf %lf %s %lf %lf\n",  pos[i].x,  pos[i].y,  pos[i].z,  at_name,  junk,   REQs[i].z   );
            REQs[i].x = 1.5; // [A]  van der Waals radius default
            REQs[i].y = 0;   // [eV] van der Waals binding energy default
            assignAtomType( i, at_name );
        }
        for(int i=0; i<nbonds; i++){
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            bond2atom[i].x = getInt( line, 0, 3 );
            bond2atom[i].y = getInt( line, 3, 6 );
            bondType[i]    = getInt( line, 6, 9 );
            //sscanf( line, "%i %i %i\n", &bond2atom[i].x, &bond2atom[i].y, &bondType[i] );
            //printf(       "%i %i %i\n",  bond2atom[i].x,  bond2atom[i].y,  bondType[i] );
            bond2atom[i].x--;
            bond2atom[i].y--;
        }
        }catch(...){
            printf( "ERROR in loadMol(%s) \n", fname );
            fclose(pFile);
            return -1;
        }
        fclose(pFile);
        return natoms + nbonds;
    }

    int load_xyz( const char * fname, int verbosity=0 ){
        if(verbosity>0)printf( "Molecule::load_xyz(%s)\n", fname );
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        int natoms; char at_name[8]; int npi,ne=0;
        const int nbuf=1024;
        char buff[nbuf]; char* line;
        line = fgets( buff, nbuf, pFile ); // number of atoms
        sscanf( line, "%i", &natoms );
        if(verbosity>0)printf( "natoms %i \n", natoms );
        line = fgets( buff, nbuf, pFile ); // comment, ignore
        allocate(natoms,0);
        for(int i=0; i<natoms; i++){
            //printf( "%i \n", i );
            Vec3d& p = pos[i];
            line     = fgets( buff, nbuf, pFile ); // comment, ignore
            //int nret = sscanf( line,            "%s %lf %lf %lf %lf %i  ",    at_name, &p.x, &p.y, &p.z, &REQs[i].z, &npis[i]  );
            //if(verbosity>1)printf   (  ".xyz[%i] %s %lf %lf %lf %lf %i\n", i, at_name,  p.x,  p.y,  p.z,  REQs[i].z,  npis[i]  );
            int nret     = sscanf( line,            "%s %lf %lf %lf %lf %i  ",    at_name, &p.x, &p.y, &p.z, &REQs[i].z, &REQs[i].w, &npis[i] );
            if(verbosity>1)printf(  ".xyz[%i] %s %lf %lf %lf %lf %i\n",        i, at_name,  p.x,  p.y,  p.z,  REQs[i].z,  REQs[i].w,  npis[i] );
            if( nret < 5 ){ REQs[i].z = 0; };
            if( nret < 6 ){ REQs[i].w = 0; };
            if( nret < 7 ){ npis[i]  =-1; };
            assignAtomType(i, at_name );
        }
        return natoms;
    }

    int loadXYZ(const char* fname, int verbosity=0 ){
        if(verbosity>0)printf( "Molecule::loadXYZ(%s)\n", fname );
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        if(params)atomTypeDict=&(params->atomTypeDict);
        const int nbuf=1024;
        char buff[nbuf];
        char * line;
        int nl;
        line = fgets( buff, nbuf, pFile ); //printf("%s",line);
        sscanf( line, "%i \n", &natoms );
        //printf("natoms %i \n", natoms );
        allocate(natoms,0);
        //allocate(natoms,nbonds);
        line = fgets( buff, nbuf, pFile ); // comment
        for(int i=0; i<natoms; i++){
            //char ch;
            char at_name[8];
            double junk;
            line = fgets( buff, nbuf, pFile );  //printf(">>%s<<",line);
            int nret = sscanf( line, "%s %lf %lf %lf %lf\n", at_name, &pos[i].x, &pos[i].y, &pos[i].z, &REQs[i].z, &REQs[i].w, &npis[i] );
            if( nret < 5 ){ REQs[i].z= 0; };
            if( nret < 6 ){ REQs[i].w= 0; };
            if( nret < 7 ){ npis[i]  =-1; };
            //if(verbosity>1)
            //printf(       ".xyz[%i] %s p(%lf,%lf,%lf) Q %lf    |nret(%i)\n", i,  at_name, pos[i].x,  pos[i].y,  pos[i].z,   REQs[i].z, nret );
            // atomType[i] = atomChar2int( ch );
            auto it = atomTypeDict->find( at_name );
            if( it != atomTypeDict->end() ){
                atomType[i] = it->second;
            }else{
                //atomType[i] = atomChar2int( at_name[0] );
                atomType[i] = -1;
            }
            //printf( " i %i name %s ityp %i \n", i, at_name, atomType[i] );
        }
        fclose(pFile);
        //printf( "atypNames.size() %i \n", atypNames->size() );
        //for ( auto element : *atypNames ){
	    //    printf(  "atypNames[%s]-> %i \n", element.first.c_str(), element.second );
        //}
        //printf("loadXYZ DONE \n");
        return natoms;
    }

    int loadXYZ_bas(const char* fname, int verbosity=0 ){
        // xxxxx.xxxxyyyyy.yyyyzzzzz.zzzz aaaddcccssshhhbbbvvvHHHrrriiimmmnnneee
        // http://www.daylight.com/meetings/mug05/Kappler/ctfile.pdf
        FILE * pFile = fopen(fname,"r");
        if( pFile == NULL ){
            printf("cannot find %s\n", fname );
            return -1;
        }
        char buff[1024];
        char * line;
        int nl;
        line = fgets( buff, 1024, pFile ); //printf("%s",line);
        sscanf( line, "%i \n", &natoms );
        //printf("natoms %i \n", natoms );
        //allocate(natoms,0);
        nbonds=0;
        allocate(natoms,nbonds);
        //line = fgets( buff, 1024, pFile ); // comment
        for(int i=0; i<natoms; i++){
            //char ch;
            //char at_name[8];
            //double junk;
            line = fgets( buff, 1024, pFile );  //printf("%s",line);
            double Q;
            int nret = sscanf( line, "%i %lf %lf %lf %lf\n", &atomType[i], &pos[i].x, &pos[i].y, &pos[i].z, &Q );
            if( nret >= 5 ){  REQs[i].z=Q; }else{ REQs[i].z=0; };
            if(verbosity>1)printf(       ".bas[%i] %i %lf %lf %lf  %lf    n", i,  atomType[i], pos[i].x,  pos[i].y,  pos[i].z,   REQs[i].z  );
            // atomType[i] = atomChar2int( ch );
        }
        fclose(pFile);
        return natoms;
    }

    void printAtom2Bond()const{
        if(atom2bond==0){ printf("# Molecule.printAtom2Bond() atom2bond == NULL ; return \n"); return; };
        printf("# Molecule.printAtom2Bond()\n");
        int * a2b = atom2bond;
        for(int ia=0; ia<natoms; ia++){
            int nb = atom_nb[ia];
            printf("%i :", ia);
            for(int i=0; i<nb; i++){
                //printf(" %i", a2b[i] );
                printf(" %i", otherAtom( bond2atom[a2b[i]], ia) );
            }
            printf("\n");
            a2b+= NBMAX;
        }
    }
    
    void printBondsInfo()const{
        if(bond2atom==0){ printf("# Molecule.printBondsInfo() bond2atom == NULL ; return \n"); return; };
        printf("# Molecule.printBondsInfo()\n");
        for(int i=0; i<nbonds; i++){
            printf( "bond[%i] atoms(%i,%i) typs(%i,%i)\n" , i, bond2atom[i].a, bond2atom[i].b, atomType[bond2atom[i].a], atomType[bond2atom[i].b] );
        }
    }

    void printAtomInfo()const{
        printf(" # Molecule.printAtomInfo() : \n" );
        for(int i=0; i<natoms; i++){
            int npi=0;
            if(npis) npi=npis[i];
            printf( "atom[%i] t %i pos(%g,%g,%g) REQs(%g,%g,%g) npi %i \n", i, atomType[i], pos[i].x,pos[i].y,pos[i].z, REQs[i].x, REQs[i].y, REQs[i].z, npi );
        }
    }

    void printAngleInfo()const{
        if(ang2bond==0){ printf("# Molecule.printAngleInfo() ang2bond == NULL ; return \n"); return; };
        printf(" # Molecule.printAngleInfo()\n");
        for(int i=0; i<nang; i++){
            double a0=0; if(ang0s)a0=ang0s[i];
            printf( "angle[%i|%i,%i] a0 %g[rad] \n", i, ang2bond[i].a, ang2bond[i].b, a0 );
        }
    }

    int loadByExt( const std::string&  fname, int verbosity=0 ){
		int idot = fname.find_last_of("."); 
		std::string ext = fname.substr( idot + 1);
        printf("'%s' '%s'\n", fname.c_str(), ext.c_str());
        if     (ext=="bas"){ return loadXYZ_bas(fname.c_str(),verbosity); }
        else if(ext=="xyz"){ return loadXYZ    (fname.c_str(),verbosity); }
        else if(ext=="mol"){ return loadMol    (fname.c_str(),verbosity); }
        return -1;
    }

    void dealloc(){
        _dealloc( pos       );
        _dealloc( bond2atom );
        _dealloc( REQs      );
        _dealloc( atomType  );
        _dealloc( npis      );
        _dealloc( bondType  );
        _dealloc( atom_nb   );
        _dealloc( atom2bond );
        _dealloc( ang2bond  );
    }
    ~Molecule(){
        dealloc();
    };


    void init_xyz( const char* fname, MMFFparams* params_=0, bool bBonds=true ){
        if(params_)bindParams( params_ );
        //atomTypeDict = params.atomTypeDict;
        loadXYZ( fname, 1 );
        if(params_)assignREQs( *params_ );
        //printAtomInfo();
        findBonds_brute(0.5,true);
        //printBondsInfo();
    }

};

#endif






