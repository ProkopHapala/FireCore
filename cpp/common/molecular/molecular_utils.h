
#ifndef molecular_utils_h
#define molecular_utils_h

#include <string>
#include <vector>
#include <unordered_map>

bool isnan(Vec3d&  v){ return (isnan(v.x)||isnan(v.y)||isnan(v.z)); }
bool isnan(Vec3f&  v){ return (isnan(v.x)||isnan(v.y)||isnan(v.z)); }
bool isnan(Quat4d& v){ return (isnan(v.x)||isnan(v.y)||isnan(v.z)||isnan(v.w)); }
bool isnan(Quat4f& v){ return (isnan(v.x)||isnan(v.y)||isnan(v.z)||isnan(v.w)); }

#define _printIfNan(var)    if(isnan(var)){printf("_printIfNan(%s)= ",#var);print(var);puts("");} 

template<typename T,typename Func>
bool ckeckNaN(int n, int m, T* xs, Func func, bool bPrint=true ){
    bool ret = false;
    for(int i=0; i<n;i++){
        bool b=false;
        for(int j=0; j<m;j++){
            int ij=i*m+j;
            b|=isnan( xs[ij] );
            b|=isinf( xs[ij] );
        }
        ret |= b;
        if(b  && bPrint ){
            func();
            printf("[%i](", i );
            for(int j=0; j<m;j++){
                int ij=i*m+j;
                printf("%g,", xs[ij] );
            }
            printf(")\n");
        }
    }
    return ret;
}

template<typename T>
bool checkNumRange( int i, T val, T min, T max, const char* pre, bool bPrint=true, bool bExit=false ){
    if( ((val<min)||(val>max)) ){
        if(bPrint){
            printf("%s[%i]==(%g) out of range (%g:%g)\n", pre, i, val, min, max);
        }
        if(bExit)exit(0);
        return true;
    }
    return false;
}

template<typename T>
bool ckeckRange(int n, int m, T* xs, T min, T max, const char* pre, bool bPrint=true ){
    bool ret = false;
    for(int i=0; i<n;i++){ 
        bool b=false;
        for(int j=0; j<m;j++){  
            int ij=i*m+j;
            T val = xs[ij];
            b|=((val<min)||(val>max)); 
        }
        if(b && bPrint ){
            printf("%s[%i](", pre, i );
            for(int j=0; j<m;j++){  int ij=i*m+j; printf("%g,", xs[ij] );   }
            printf(") outof (%g,%g)\n", (double)min, (double)max );
        }
        ret |= b;
        //ret |= checkNumRange<T>(i,m, xs[ij],min,max,pre,bPrint); 
    }
    return ret;
}

bool ckeckNaN_d(int n, int m, double* xs, const char* pre, bool bPrint=true ){
    bool ret = false;
    for(int i=0; i<n;i++){
        bool b=false;
        for(int j=0; j<m;j++){
            int ij=i*m+j;
            b|=isnan( xs[ij] );
            b|=isinf( xs[ij] );
        }
        ret |= b;
        if(b && bPrint ){
            printf("%s[%i](", pre, i );
            for(int j=0; j<m;j++){
                int ij=i*m+j;
                printf("%g,", xs[ij] );
            }
            printf(")\n");
        }
    }
    return ret;
}

int whereNaN_d(int n, int m, double* xs, const char* pre ){
    for(int i=0; i<n;i++){
        bool b=false;
        for(int j=0; j<m;j++){
            int ij=i*m+j;
            b|=isnan( xs[ij] );
            b|=isinf( xs[ij] );
        }
        if(b){ return i; }
    }
    return -1;
}

bool ckeckNaN_f(int n, int m, float* xs, const char* pre, bool bPrint=true ){
    bool ret = false;
    for(int i=0; i<n;i++){
        bool b=false;
        for(int j=0; j<m;j++){
            int ij=i*m+j;
            b|=isnan( xs[ij] );
            b|=isinf( xs[ij] );
        }
        ret |= b;
        if(b && bPrint ){
            printf("%s[%i](", pre, i );
            for(int j=0; j<m;j++){
                int ij=i*m+j;
                printf("%g,", xs[ij] );
            }
            printf(")\n");
        }
    }
    return ret;
}

void nameList( std::vector<std::string>& names, const std::string& s ){
    size_t pos = 0;
    while ( (pos = s.find(' ', pos) ) != std::string::npos ){
        names.push_back(  s.substr(0, pos) );
    }
}

void listToMap( const std::vector<std::string>& names, std::unordered_map<std::string,int>& dct ){
    for( int i=0; i<names.size(); i++ ){  dct[ names[i] ]=i; }
}

void countingDict( int i0, std::unordered_map<std::string,int>& dct, const std::string& s ){
    size_t pos = 0;
    int i=i0;
    while ( (pos = s.find(' ', pos) ) != std::string::npos ){
        dct[ s.substr(0, pos) ] = i;
        i++;
    }
}

void makeDefaultAtomTypeDict( std::vector<std::string>& names, std::unordered_map<std::string,int>& dct ){
    nameList( names,
        "H He "
        "Li Be B C N O F Ne "
        "Na Mg Al Si P S Cl Ar "
        "K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr "
        "Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
        "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn "
        "Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt "
    );
    listToMap( names, dct );
}
void makeDefaultAtomTypeDict( std::vector<std::string>*& names, std::unordered_map<std::string,int>*& dct ){
    if(names==0) names = new std::vector       <std::string>    ();
    if(dct==0)   dct   = new std::unordered_map<std::string,int>();
    makeDefaultAtomTypeDict( *names, *dct );
}

/*
void writeXYZ( FILE* pfile, int n, const int* atypes, const Vec3d* apos, const std::vector<std::string>& atomTypeNames, const char* comment="#comment" ){
    fprintf(pfile, "%i\n", n );
    fprintf(pfile, "%s \n", comment );
    for(int i=0; i<n; i++){
        //printf( "DEBUG writeXYZ()[%i] \n", i );
        int ityp   = atypes[i];
        const Vec3d&  pi = apos[i];
        //printf( "write2xyz %i %i (%g,%g,%g) %s \n", i, ityp, pi.x,pi.y,pi.z, params->atypes[ityp].name );
        fprintf( pfile, "%s   %15.10f   %15.10f   %15.10f \n", atomTypeNames[ityp].c_str(), pi.x,pi.y,pi.z );
    };
}

int saveXYZ( const char * fname, int n, const int* atypes, const Vec3d* apos, const std::vector<std::string>& atomTypeNames, const char* comment="#comment"  ){
    FILE* pfile = fopen(fname, "w");
    if( pfile == NULL ) return -1;
    writeXYZ( pfile, n, atypes, apos, atomTypeNames, comment );
    fclose(pfile);
    return n;
}
*/

#endif
