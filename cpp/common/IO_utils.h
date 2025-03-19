
#ifndef  IO_utils_h
#define  IO_utils_h

#include <stdlib.h>
#include <stdio.h>
#include <cstdarg>
#include <cstring>

#include <vector>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>

#include <string>
#include <vector>
#include <unordered_map>

#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

const int N_CHAR_TMP = 256;

template <typename T> void toBuff  (  const T& v, void* buff, int& i){ (*(T*)(buff+i))=v; i+=sizeof(T); };
template <typename T> void fromBuff(        T& v, void* buff, int& i){ v=*((T*)(buff+i)); i+=sizeof(T); };

inline int print(       char*   v){ return printf( "%s", v ); };
inline int print(       float   v){ return printf( "%g", v ); };
inline int print(       double  v){ return printf( "%g", v ); };
inline int print(       int     v){ return printf( "%i", v ); };
// Moved to Vec2.h, Vec3.h, quaternion.h
//inline int print( const Vec2f&  v){ return printf( "%lg %g", v.x, v.y ); };
//inline int print( const Vec2d&  v){ return printf( "%lg %g", v.x, v.y ); };
//inline int print( const Vec2i&  v){ return printf( "%li %i", v.x, v.y ); };
//inline int print( const Vec3f&  v){ return printf( "%g %g %g", v.x, v.y, v.z ); };
//inline int print( const Vec3d&  v){ return printf( "%g %g %g", v.x, v.y, v.z ); };
//inline int print( const Vec3i&  v){ return printf( "%i %i %i", v.x, v.y, v.z ); };
//inline int print( const Quat4f& v){ return printf( "%g %g %g %g", v.x, v.y, v.z, v.w ); };
//inline int print( const Quat4d& v){ return printf( "%g %g %g %g", v.x, v.y, v.z, v.w ); };
//inline int print( const Quat4i& v){ return printf( "%i %i %i %i", v.x, v.y, v.z, v.w ); };

inline int fprint(FILE* sbuff,       char*   v){ return fprintf(sbuff, "%s ", v ); };
inline int fprint(FILE* sbuff,       float   v){ return fprintf(sbuff, "%g ", v ); };
inline int fprint(FILE* sbuff,       double  v){ return fprintf(sbuff, "%g ", v ); };
inline int fprint(FILE* sbuff,       int     v){ return fprintf(sbuff, "%i ", v ); };
inline int fprint(FILE* sbuff, const Vec2f&  v){ return fprintf(sbuff, "%g %g ", v.x, v.y ); };
inline int fprint(FILE* sbuff, const Vec2d&  v){ return fprintf(sbuff, "%g %g ", v.x, v.y ); };
inline int fprint(FILE* sbuff, const Vec2i&  v){ return fprintf(sbuff, "%i %i ", v.x, v.y ); };
inline int fprint(FILE* sbuff, const Vec3f&  v){ return fprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int fprint(FILE* sbuff, const Vec3d&  v){ return fprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int fprint(FILE* sbuff, const Vec3i&  v){ return fprintf(sbuff, "%i %i %i ", v.x, v.y, v.z ); };
inline int fprint(FILE* sbuff, const Quat4f& v){ return fprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int fprint(FILE* sbuff, const Quat4d& v){ return fprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int fprint(FILE* sbuff, const Quat4i& v){ return fprintf(sbuff, "%i %i %i %i ", v.x, v.y, v.z, v.w ); };

inline int sprint(char* sbuff,       char*   v){ return sprintf(sbuff, "%s ", v ); };
inline int sprint(char* sbuff,       float   v){ return sprintf(sbuff, "%g ", v ); };
inline int sprint(char* sbuff,       double  v){ return sprintf(sbuff, "%g ", v ); };
inline int sprint(char* sbuff,       int     v){ return sprintf(sbuff, "%i ", v ); };
inline int sprint(char* sbuff, const Vec2f&  v){ return sprintf(sbuff, "%g %g ", v.x, v.y ); };
inline int sprint(char* sbuff, const Vec2d&  v){ return sprintf(sbuff, "%g %g ", v.x, v.y ); };
inline int sprint(char* sbuff, const Vec2i&  v){ return sprintf(sbuff, "%i %i ", v.x, v.y ); };
inline int sprint(char* sbuff, const Vec3f&  v){ return sprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int sprint(char* sbuff, const Vec3d&  v){ return sprintf(sbuff, "%g %g %g ", v.x, v.y, v.z ); };
inline int sprint(char* sbuff, const Vec3i&  v){ return sprintf(sbuff, "%i %i %i ", v.x, v.y, v.z ); };
inline int sprint(char* sbuff, const Quat4f& v){ return sprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int sprint(char* sbuff, const Quat4d& v){ return sprintf(sbuff, "%g %g %g %g ", v.x, v.y, v.z, v.w ); };
inline int sprint(char* sbuff, const Quat4i& v){ return sprintf(sbuff, "%i %i %i %i ", v.x, v.y, v.z, v.w ); };

inline int sscan(char* sbuff, char*&  v){ int n; sscanf(sbuff, "%s%n",  v  , &n );                                return n; };
inline int sscan(char* sbuff, float&  v){ int n; sscanf(sbuff, "%f%n",  &v , &n );                                return n; };
inline int sscan(char* sbuff, double& v){ int n; sscanf(sbuff, "%lf%n", &v , &n );                                return n; };
inline int sscan(char* sbuff, int&    v){ int n; sscanf(sbuff, "%i%n",  &v , &n );                                return n; };
inline int sscan(char* sbuff, Vec2f&  v){ int n; sscanf(sbuff, "%f %f%n", &v.x, &v.y , &n);                       return n; };
inline int sscan(char* sbuff, Vec2d&  v){ int n; sscanf(sbuff, "%lf %lf%n", &v.x, &v.y , &n);                     return n; };
inline int sscan(char* sbuff, Vec2i&  v){ int n; sscanf(sbuff, "%i %i%n", &v.x, &v.y , &n);                       return n; };
inline int sscan(char* sbuff, Vec3f&  v){ int n; sscanf(sbuff, "%f %f %f%n",    &v.x, &v.y, &v.z , &n);           return n; };
inline int sscan(char* sbuff, Vec3d&  v){ int n; sscanf(sbuff, "%lf %lf %lf%n",    &v.x, &v.y, &v.z , &n);        return n; };
inline int sscan(char* sbuff, Vec3i&  v){ int n; sscanf(sbuff, "%i %i %i%n",    &v.x, &v.y, &v.z , &n);           return n; };
inline int sscan(char* sbuff, Quat4f& v){ int n; sscanf(sbuff, "%f %f %f %f%n", &v.x, &v.y, &v.z, &v.w , &n);     return n; };
inline int sscan(char* sbuff, Quat4d& v){ int n; sscanf(sbuff, "%lf %lf %lf %lf%n", &v.x, &v.y, &v.z, &v.w , &n); return n; };
inline int sscan(char* sbuff, Quat4i& v){ int n; sscanf(sbuff, "%i %i %i %i%n",    &v.x, &v.y, &v.z, &v.w , &n);  return n; };

//#define OPT(OP,) OP(tok)
//#define OPT(OP,tok) OP(tok)
//#define DO_10(OP,tok,...) OP(tok) DO_9(__VA_ARGS__)

#define _sprintf( sbuff, format, ... ){ sbuff+=sprintf( sbuff, format, __VA_ARGS__); }
#define _sprint ( sbuff, tok ){ sbuff+=sprint( sbuff, tok ); }
#define _sscan  ( sbuff, tok ){ sbuff+=sscan( sbuff, tok ); }

#define _lprint(         tok        ) {         printf(        " %s: ",#tok);          print(      tok); }
#define _lfprint( fout,  tok        ) {         fprintf( fout,  " %s: ",#tok);        fprint(fout, tok); }
#define _lsprint( sbuff, tok        ) { sbuff+= sprintf( sbuff, " %s: ",#tok); sbuff+=sprint(sbuff,tok); }
#define _Lprint(         tok, s     ) {         printf(         " %s: ",s);          print(      tok); }
#define _Lfprint( fout,  tok, s     ) {         fprintf( fout,  " %s: ",s);        fprint(fout, tok); }
#define _Lsprint( sbuff, tok, s     ) { sbuff+= sprintf( sbuff, " %s: ",s); sbuff+=sprint(sbuff,tok); }
//define _lscan  ( sbuff, tok       ) { sbuff+= sprintf( sbuff, " %s: ",#tok); sbuff+=sprint(sbuff,tok); }
//#define _lsprint_( sbuff, tok, ... ) { _lsprint(sbuff,tok); _lsprint_(sbuff,__VA_OPT__); }
//#define _lsprint_( sbuff, tok, ... ) { _lsprint(sbuff,tok); _lsprint_(sbuff,__VA_ARGS__); } // C-macros are stupid. non recursive, this would not work

#define _toDict( mp, tok ){ mp[#tok] = tok; }
#define _fromDict( mp, tok ){ tok = mp[#tok]; }

inline bool file_exist(const char* fname) { if (FILE *file = fopen(fname, "r")) { fclose(file); return true; } else { return false; } }

inline char * fgetsNonComment(char * str, int num, FILE * stream, char commentChar, int nMaxTry = 100 ){
    char* s = nullptr;
    // no more than 100 comments expected, ensure we don't get stuck in infinite loop
    for (int i = 0; i < nMaxTry; i++) {
        s = fgets(str, num, stream);
        //printf("fgetsNonComment [%i] '%s'", i, s );
        if (s){
            if (s[0] == commentChar)continue;
        }
        break;
    }
    //printf( "fgetsNonComment '%s'", s );
    return s;
}

inline int readMatrix( const char* fname, int nrow, int ncol, double* buff, bool bTranspose=0 ){
    FILE *file = fopen(fname, "r");
    if ( file==0 )[[unlikely]]{
        printf( "ERROR in readMatrix(%s): no such file \n", fname );
        return -1;
    }
    int di=ncol,dj=1;
    if(bTranspose){ di=1; dj=nrow; }
    //char [];
    for(int i=0; i<nrow; i++){
        for(int j=0; j<ncol; j++){
            double val;
            fscanf( file, "%lf \n", &val );
            //printf( "readMatrix[%i,%i] %g \n", i, j, val );
            buff[i*di+j*dj] = val;
        }
    }
    return nrow*ncol;
}

//  TODO:
// Universal data loading idea:
//  - load all data to std::map<string,string> "craft.velocity"->"0.0 1.0 3.0"
//  - pull named tokens from map and parse them to data   "0.0 1.0 3.0" - > Vec3d{0.0,1.0,30.}

// TODO: GUI should have also such macros
//       - bind pointers &float to data item to GUI sliders


/*
// making named items in dictionary

#define _toDict  ( char* mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }
#define _fromDict( char* mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }

template <typename T> void toBuff  (  const T& v, void* buff, int& i){ (*(T*)(buff+i))=v; i+=sizeof(T); };

#define _toSDict  ( mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }
#define _fromSDict( mp, tok ){ sbuff+=sprint("%s",#tok); sbuff+=sprint(sbuff,tok); }
*/


/*
inline char* file2str(char const  *fname ){
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
    char *buff = NULL;
	fseek (fptr, 0, SEEK_END);
    int length = ftell(fptr);
    buffer = new char[length];
    fseek (fptr, 0, SEEK_SET);
    if (buffer){ fread (buff, 1, length, fptr); }
    fclose (fptr);
	return buff;
}
*/

inline int fileExist(const char * fname ){
    FILE *file;
    if ( (file = fopen(fname, "r")) ) {
        fclose(file);
        return 1;
    } else {
        return 0;
    }
}

inline bool checkAllFilesExist( int n, const char** fnames, bool bPrint=true, bool bExit=true ){
    //printf( "checkAllFilesExist() n=%i \n", n );
    bool bExist=true;
    for(int i=0; i<n; i++){
        const char* fname = fnames[i];
        if( 0==fname ){ 
            if(bExit){ printf("ERROR in checkAllFilesExist() fnames[%i]==null => exit() \n", i ); exit(0); }
            continue;
        }
        //printf( "checkAllFilesExist()[%i] '%s' \n", i, fname );
        if(fname==0) continue;
        FILE* f=fopen( fname,"rb"); 
        //printf( "checkAllFilesExist()[%i] '%s' @file=%li \n", i, fname, (long)f );
        if(0==f){ bExist=false; if(bPrint)printf("File(%s) Not Found\n", fname); }else{ fclose(f); };
    };
    //if(fname_Paul){ FILE* f=fopen( fname_Paul,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Paul); bExist=false; }else{ fclose(f); };} 
    //if(fname_Lond){ FILE* f=fopen( fname_Lond,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Lond); bExist=false; }else{ fclose(f); };} 
    //if(fname_Coul){ FILE* f=fopen( fname_Coul,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Coul); bExist=false; }else{ fclose(f); };} 
    return bExist;
}

inline bool tryMakeDir( const char* dirname, bool bExit=true, bool bPrint=true ){
    struct stat statbuf;
    bool bDirReady=false;
    if (stat(dirname, &statbuf) != 0) {   // Check if directory exists
          if ( mkdir(dirname, 0755) == -1 ){ if(bPrint)printf("ERROR in tryMakeDir() cannot mkdir(%s)\n",                       dirname ); if(bExit)exit(0); }else{ bDirReady=true; }
    }else if ( !S_ISDIR(statbuf.st_mode)  ){ if(bPrint)printf("ERROR in tryMakeDir() path(%s) exists but is not a directory\n", dirname ); if(bExit)exit(0); }else{ bDirReady=true; }
    return bDirReady;
}

inline bool tryChangeDir( const char* dirname, bool bExit=true, bool bPrint=true, char* msg="tryChangeDir()" ){
    if (chdir(dirname) == -1) { if(bPrint)printf("ERROR in %s chdir(%s) => exit()\n", msg, dirname ); if(bExit)exit(0); return true; }
    return false;
}

//#include "Tree.h"

inline bool printDir( char* dirName ){
    // from : https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c
    DIR *dir;
    struct dirent *ent;
    if ((dir = opendir (dirName)) != NULL) {
        while ((ent = readdir (dir)) != NULL) {
            printf ("%s\n", ent->d_name);
        }
        closedir(dir);
        return true;
    }else {
        //perror ("");
        return false;
    }
}

// list files in directory
//  https://stackoverflow.com/questions/612097/how-can-i-get-the-list-of-files-in-a-directory-using-c-or-c


inline int listDirContaining( char * dirName, char * fname_contains, std::vector<std::string>& fnames_found ){
    DIR *dir=NULL;
    int n=0;
    struct dirent *ent=NULL;
    int i=0;
    if ( (dir = opendir( dirName )) != NULL) {
        while ( (ent = readdir (dir)) != NULL) {
            char* found = strstr( ent->d_name, fname_contains );
            if( found ){
                printf("%i %s\n", i, ent->d_name);
                fnames_found.push_back( ent->d_name );
                i++;
            }
        }
        n++;
        closedir(dir);
    } else {
        printf("Cannot open directory %s \n", dirName );
        return -1;
    }
    return i;
}

/*
int dir2tree(TreeViewTree& node, char * name, int level ){

    if (niters >100) return -1;
    niters++;

    if((name[0]=='.'))return 0;

    for(int i=0; i<level; i++) printf("_");

    node.content.caption = name;
    DIR *dir=NULL;
    struct dirent *ent=NULL;

    if( chdir(name)==0 ){
    //if( (dir = opendir( name )) != NULL){
        dir = opendir( "." );
        printf("dir '%s' | %i \n", name, level );
        while( (ent = readdir(dir)) != NULL){
            node.branches.push_back( TreeViewTree() );
            dir2tree( node.branches.back(), ent->d_name, level+1 );
        }
        closedir(dir);
        chdir("..");
    }else{
        printf("leaf '%s' | %i \n", name, level );
    }
    return 0;
}
*/

template <typename Func>
int processFileLines( const char * fname, Func func ){
    FILE * pFile;
    const int nbuff = 4096;
    char str[nbuff];
    pFile = fopen ( fname , "r");
    if (pFile == NULL){ printf("file not found: %s \n", fname ); return(-1); };
    int n=0;
    while ( fgets( str , nbuff, pFile) != NULL ){
        if (str[0]=='#') continue;
        func( str );
        n++;
    }
    fclose(pFile);
    return n;
}

inline char* stripWhite( char* s ){
    for(;;s++){
        char c=*s;
        if(c=='\0') return s;
        if(c>=33) break;
    }
    for( char* s_=s;;s_++ ){
        if(*s_<33){
            *s_='\0';
            return s;
        }
    }
}

inline int strcmp_noWhite( const char * s1, const char * s2 ){
    while( true ){
        char c=*s1;
        if(c=='\0') return -1;
        if(c>=33) break;
        s1++;
    }
    while( true ){
        char c2=*s2; if(c2=='\0') break;
        char c1=*s1;
        int d = c1-c2;
        if(d) return d;
        s1++;s2++;
    }
    return 0;
}

inline int str2enum( const char * str, int nnames, const char **names ){
    for(int i=0; i<nnames; i++ ){
        //printf( ">>%s<< >>%s<<\n", str, names[i] );
        if( strcmp_noWhite( str, names[i] ) == 0 ) return i;
    }
    return -1;
}

inline int saveBin( const char *fname, int n, char * data ){
    FILE *ptr_myfile=0;
    ptr_myfile=fopen( fname,"wb");
    if (!ptr_myfile){ printf("saveBin(): Unable to open file `%s` for writing !!!\n", fname ); return -1; }
    int nchar = 1024;
    //printf( "saveBin %s nbyte=%i @data=%li \n", fname, n, data );
    for( int i=1; i<=n; i+=nchar ){
        int len = nchar;
        if( (n-i)<nchar ) len = (n-i);
        fwrite( data+i, sizeof(char), len, ptr_myfile);
    }
    fclose(ptr_myfile);
    return 0;
}

inline int loadBin( const char *fname, int n, char * data ){
    FILE *ptr_myfile=0;
    //printf("loadBin %s \n", fname );
    ptr_myfile=fopen( fname,"rb");
    if (!ptr_myfile){ printf("loadBin(): Unable to open file `%s` for reading !!! \n", fname ); return -1; }
    int nchar = 1024;
    for( int i=1; i<=n; i+=nchar ){
        int len = nchar;
        if( (n-i)<nchar ) len = (n-i);
        fread( data+i, sizeof(char), len, ptr_myfile);
    }
    fclose(ptr_myfile);
    return 0;
}


inline int save_npy(const char *fname, int ndims, int* shape, const char *data, int nBytePerElement=8, const char *dtype="<f8", bool fortran_order = false) {
    FILE *ptr_myfile = fopen(fname, "wb");
    if (!ptr_myfile) {   printf("saveNpy(): Unable to open file `%s` for writing !!!\n", fname); return -1; }
    // Write magic string and version number (1.0)
    const char magic_string[] = "\x93NUMPY";
    fwrite(magic_string, sizeof(char), 6, ptr_myfile);
    char version[2] = {1, 0};  // Version 1.0
    fwrite(version, sizeof(char), 2, ptr_myfile);
    // Create the shape string using sprintf
    char shape_str[256];
    int pos = 0;  // Track position in buffer
    pos += sprintf(shape_str + pos, "(");
    for (int i = 0; i < ndims; ++i) {
        pos += sprintf(shape_str + pos, "%d", shape[i]);  // Add shape dimension
        //printf( "saveNpy()[dim=%i] shape_str=`%s`\n", i, shape_str );
        if ( (ndims==1) || (i<ndims-1)   ) {
            pos += sprintf(shape_str + pos, ", ");  // Add comma and space if not the last dimension
        }
    }
    sprintf(shape_str + pos, ")");  // Close the shape string
    //printf( "saveNpy() shape_str=`%s`\n", shape_str );
    // Create the header
    char header[1024];  // You can adjust size if needed
    snprintf(header, sizeof(header), "{'descr': '%s', 'fortran_order': %s, 'shape': %s, }\n", dtype, fortran_order ? "True" : "False", shape_str );
    printf( "saveNpy(%s) header=%s", fname, header );
    // Calculate padding for alignment to 16 bytes
    int header_len = strlen(header);
    int padding = 16 - ((10 + header_len) % 16);  // 10 bytes for magic, version, header length
    header_len += padding;
    // Write header length (in little endian)
    unsigned short header_size = header_len;  // max header size should fit in unsigned short for version 1.0
    fwrite(&header_size, sizeof(header_size), 1, ptr_myfile);
    // Write the actual header
    fwrite(header, sizeof(char), strlen(header), ptr_myfile);
    for (int i = 0; i < padding; i++) {  fputc(' ', ptr_myfile);   } // Padding with spaces
    //fputc('\n', ptr_myfile);  // Final newline
    int n_elements = 1;  for (int i=0; i<ndims; ++i) { n_elements *= shape[i]; }   // Calculate the total number of elements
    fwrite(data, nBytePerElement, n_elements, ptr_myfile);   // Write the binary data
    fclose(ptr_myfile);
    return 0;
}


class NumpyFile{ public:
    int ndims; 
    int shape[8]; 
    int ntot;
    int nBytePerElement;
    char dtype[8];
    char* data=0;

    NumpyFile( const char* fname ){
        if (load(fname) != 0) {
            printf( "ERROR NumpyFile::NumpyFile(%s): Unable to load file.\n", fname ); exit(0);
            // Handle load error appropriately
            data = 0;
        }
    }

    void save(char* fname, bool fortran_order = false ){
        save_npy(fname, ndims, shape, data, nBytePerElement, dtype, fortran_order );
    }

    void print(){
        printf("NumpyFile() dtype=`%s` nbyte=%i ntot=%i shape(", dtype, nBytePerElement, ntot );
        for(int i=0; i<ndims; i++){ printf( "%i,", shape[i] ); }
        printf(")\n");
    }

    inline int load(const char *fname ) {
        //char cwd[128]; getcwd(cwd, 128); printf( "NumpyFile::load(%s/%s)\n", cwd, fname );
        FILE *ptr_myfile = fopen(fname, "rb");
        if (!ptr_myfile) {   printf("ERROR NumpyFile::load(%s): unable to open file for reading !!!\n", fname);  return -1;  }

        // Read magic string
        char magic_string[7];     if (fread(magic_string, 1, 6, ptr_myfile) != 6) { printf("ERROR NumpyFile::load(%s): Unable to read magic string.\n", fname);              fclose(ptr_myfile);  return -1; }
        magic_string[6] = '\0';   if (strncmp(magic_string, "\x93NUMPY", 6) != 0) { printf("ERROR NumpyFile::load(%s): Invalid magic string (%s).\n", fname, magic_string ); fclose(ptr_myfile);  return -1; }
        //printf("NumpyFile::load() magic_string=%s\n", magic_string);

        // Read version
        unsigned char version[2];   if (fread(version, 1, 2, ptr_myfile) != 2) {     printf("ERROR NumpyFile::load(%s): Unable to read version.\n", fname);  fclose(ptr_myfile);  return -1;   }
        //printf("NumpyFile::load() version=%d.%d\n", version[0], version[1]);

        // Read header length based on version
        uint32_t header_length;
        if (version[0] == 1) { // Version 1.0: 2-byte little endian header length
            uint16_t header_len_16;
            if (fread(&header_len_16, 1, 2, ptr_myfile) != 2) {  printf("ERROR NumpyFile::load(%s): Unable to read 2-byte header length.\n", fname); fclose(ptr_myfile);     return -1;     }
            header_length = header_len_16;
        }else if (version[0] == 2 || version[0] == 3) {   // Version 2.0 and above: 4-byte little endian header length
            if (fread(&header_length, 1, 4, ptr_myfile) != 4) {  printf("ERROR NumpyFile::load(%s): Unable to read 4-byte header length.\n", fname);  fclose(ptr_myfile);  return -1;    }
        } else                                                {  printf("ERROR NumpyFile::load(%s): Unsupported .npy version %d.%d\n", fname, version[0], version[1]);  fclose(ptr_myfile); return -1;  }
        //printf("NumpyFile::load() header_length=%i\n", header_length);

        // Read the header
        char header[header_length + 1];
        if (fread(header, 1, header_length, ptr_myfile) != header_length) { printf("ERROR NumpyFile::load(%s): Unable to read header.\n", fname); fclose(ptr_myfile);   return -1;  }
        header[header_length] = '\0';  // Null-terminate the header for easier parsing
        //printf("NumpyFile::load() header='%s'\n", header);

        // Extract dtype
        char* descr_start = strstr(header, "'descr': '");
        if (!descr_start)                                                 { printf("ERROR NumpyFile::load(%s): 'descr' not found in header.\n", fname); fclose(ptr_myfile); return -1;}
        descr_start += strlen("'descr': '");
        char* descr_end = strchr(descr_start, '\'');
        if (!descr_end) { printf("ERROR NumpyFile::load(%s): 'descr' value not properly terminated.\n", fname); fclose(ptr_myfile); return -1; }
        int descr_len = descr_end - descr_start;
        if (descr_len >= sizeof(dtype)) { printf("ERROR NumpyFile::load(%s): 'descr' value too long.\n", fname); fclose(ptr_myfile); return -1; }
        strncpy(dtype, descr_start, descr_len);
        dtype[descr_len] = '\0';  // Null-terminate dtype string
        //printf("NumpyFile::load() dtype='%s'\n", dtype);

        // Determine bytes per element based on dtype
        if      (strcmp(dtype, "<f8") == 0) { nBytePerElement = 8; } 
        else if (strcmp(dtype, "<f4") == 0) { nBytePerElement = 4; } 
        else if (strcmp(dtype, "<i4") == 0) { nBytePerElement = 4; } 
        else if (strcmp(dtype, "<i8") == 0) { nBytePerElement = 8; } 
        else if (strcmp(dtype, "<u1") == 0) { nBytePerElement = 1; } 
        else { printf("NumpyFile::load(%s): Unsupported dtype `%s`\n", fname, dtype); fclose(ptr_myfile);  return -2;  }

        // Extract shape
        char* shape_start = strstr(header, "'shape': (");
        if (!shape_start) { printf("ERROR NumpyFile::load(%s): 'shape' not found in header.\n", fname);  fclose(ptr_myfile); return -1;}
        shape_start += strlen("'shape': (");
        char* shape_end = strchr(shape_start, ')');
        if (!shape_end) { printf("ERROR NumpyFile::load(%s): 'shape' value not properly terminated.\n", fname); fclose(ptr_myfile); return -1; }
        int shape_str_len = shape_end - shape_start;
        if (shape_str_len >= 256) { printf("ERROR NumpyFile::load(%s): 'shape' string too long.\n", fname);  fclose(ptr_myfile);   return -1;  }
        char shape_str[256];
        strncpy(shape_str, shape_start, shape_str_len);
        shape_str[shape_str_len] = '\0';  // Null-terminate shape string

        // Count dimensions
        ndims = 1;
        for (char* c = shape_str; *c != '\0'; ++c) {
            if (*c == ',') ndims++;
        } if (ndims > 8) { printf("ERROR NumpyFile::load(%s): Too many dimensions (%d).\n", fname, ndims);     fclose(ptr_myfile); return -1;  }

        // Parse shape into array
        char* token = strtok(shape_str, ", ");
        int dim = 0;
        while (token && dim < ndims) {
            shape[dim++] = atoi(token);
            token = strtok(NULL, ", ");
        }

        // Calculate number of elements
        ntot = 1;
        for (int i = 0; i < ndims; ++i) {
            ntot *= shape[i];
        }

        // Allocate buffer for data
        data = new char[ntot * nBytePerElement];
        if (!data) { printf("ERROR NumpyFile::load(%s): Memory allocation failed.\n", fname); fclose(ptr_myfile); return -1; }

        // Read binary data
        size_t items_read = fread(data, nBytePerElement, ntot, ptr_myfile);
        if (items_read != (size_t)ntot) { printf("ERROR NumpyFile::load(%s): Unable to read data. Expected %d elements, got %zu.\n", fname, ntot, items_read);  delete[] data;fclose(ptr_myfile); return -1;  }

        fclose(ptr_myfile);
        return 0;
    }

};



inline char * fgets_comment( char * line, int num, FILE * stream ){
    constexpr int NMaxComment = 10;
    for(int i=0; i<NMaxComment; i++){
        char *str = fgets( line, num, stream );
        printf(">> %s", line );
        if( str[0] != '#' ) return str;
    }
    return NULL;
}

// A simple function that will read a file into an allocated char pointer buffer
inline  char* filetobuf(char const  *fname){
	FILE *fptr;
	long length;
	char *buf;
	fptr = fopen(fname, "rb");			// Open file for reading
	if(fptr==NULL){
        printf("Failed to load %s \n", fname );
	    return NULL;
    }
	fseek(fptr, 0, SEEK_END); 			// Seek to the end of the file
	length = ftell(fptr); 				// Find out how many bytes into the file we are
	//buf = (char*)malloc(length+1); 		// Allocate a buffer for the entire length of the file and a null terminator
	buf = new char[length+1];
	fseek(fptr, 0, SEEK_SET); 			// Go back to the beginning of the file
	fread(buf, length, 1, fptr); 		// Read the contents of the file in to the buffer
	fclose(fptr); 						// Close the file
	buf[length] = 0; 					// Null terminator
	return buf; 						// Return the buffer
}

inline int fileGetNextKey( FILE  *fptr, const char * keystr, char * tmp ){
	int nch = strlen(keystr);
	//printf( "fileGetNextKey [%s] %i \n", keystr, nch );
    while(fgets(tmp, N_CHAR_TMP, fptr)){
        //printf( "fileGetNextKey: %s %i \n", tmp, strncmp(tmp,keystr,nch) );
        if( strncmp(tmp,keystr,nch)==0 ) return ftell(fptr);
        //printf( "fileGetNextKey: NOOOO! \n" );
    };
    return -1;
}

inline char * fileCut( FILE * fptr, int ibeg, int iend ){
    int nch     = iend-ibeg;
    char * buff = new char[nch+1];
    fseek(fptr, ibeg, SEEK_SET);
    fread(buff, nch, 1, fptr);
    buff[nch]='\0';
    return buff;
}

// cut piece out of file
inline char* fileGetSection(const char *fname, const char * keystr, const char * endstr ){
	//printf( "fileGetSection \n" );
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
    char *buff = NULL;
	char      tmp[N_CHAR_TMP];
    int ibeg = fileGetNextKey(fptr, keystr, tmp);
    if( ibeg>=0 ){
        int iend = fileGetNextKey(fptr, endstr, tmp);
        if(iend>=0){           // cut (ibeg,iend)
            buff = fileCut(fptr, ibeg, iend-strlen(endstr)-1 );
        }
    }
    fclose (fptr);
	return buff;
}

inline int whichKey( const char* tmp, int nkey, char const*const* keys ){
    for(int ikey=0; ikey<nkey; ikey++){
        const char * key = keys[ikey];
        //printf( "whichKey[%i] %s %s \n", ikey, key, tmp );
        int i=0;
        bool match=true;
        while(key[i]!='\0'){ if(key[i]!=tmp[i])match=false; i++; }
        if(match) return ikey;
    }
    return -1;
}

inline char ** fileGetSections(const char *fname, int nkey, char const*const* keys, const char* begstr ){
	FILE  *fptr = fopen(fname, "rb");	  // Open file for reading
	if(fptr==NULL){ printf("Failed to load %s \n", fname ); return NULL; }
	char      tmp[N_CHAR_TMP];
	int nb = strlen(begstr);
	char** result = new char*[nkey];
	for(int ikey=0; ikey<nkey; ikey++){ result[ikey] = NULL; }
	int ikey=-1,i0=-1,i1;
	//printf( "fileGetSections\n" );
	while( (i1=fileGetNextKey( fptr, begstr, tmp ))>=0 ){
        //printf( " ikey %i i0  %i i1 %i \n", ikey, i0, i1 );
        //if((ikey>=0)&&(i0>=0)){
        //    //printf(" ikey %i i0  %i i1 %i \n", ikey, i0, i1 );
        //    result[ikey] = fileCut( fptr, i0, i1 );
        //}
        if((ikey>=0)&&(i0>=0)){ result[ikey] = fileCut( fptr, i0, i1 ); }
        ikey = whichKey( tmp+nb, nkey, keys );
        i0=i1;
	};
	fseek(fptr, 0, SEEK_END);
	if((ikey>=0)&&(i0>=0)) result[ikey] = fileCut( fptr, i0, ftell(fptr) ); // seaction at end of file
	fclose(fptr);
	return result;
}

inline int checkNullSections( int n, char ** sections ){
    for(int i=0; i<n; i++){ if(sections[i]==NULL) return i; }
    return 0;
}

inline void saveStr( const char * fname, const char * str ){
    int n = strlen(str);
    FILE  *fptr = fopen(fname, "wb");
    fwrite( str, n, 1, fptr);
    fclose(fptr);
}



/*
int loadColums(char const  *fname, char const  *format, ... ){
    FILE * pFile = fopen (fname,"r");
    char buff[1024];
    char * line;
    while( line = fgets( buff, 1024, pFile ) ){
        sscanf( buff, format, ... );
        printf()
    }
    va_list args;
    va_start(args, fmt);
    scanf( format );
    va_end(args);
    fclose(pFile);
}
*/

inline  int allocateIOBuffs( int nitems, char const *format, void **buffs ){
    int nbuffs = 0;
    //int ibuff  = 0;
    while (*format != '\0') {
        printf( "format %c nbuffs %i \n", *format, nbuffs );
        switch( *format ){
            case 'i': buffs[nbuffs] = new int   [nitems]; nbuffs++; break;
            case 'f': buffs[nbuffs] = new float [nitems]; nbuffs++; break;
            case 'd': buffs[nbuffs] = new double[nitems]; nbuffs++; break;
            case '2': buffs[nbuffs] = new Vec3d [nitems]; nbuffs++; break;  // TODO
            case '3': buffs[nbuffs] = new Vec2d [nitems]; nbuffs++; break;
        }
        format++;
    }
    return nbuffs;
}

inline  int loadColumns( char const  *fname, char const *format, void **buffs ){
    FILE * pFile = fopen(fname,"r");
    if( pFile == NULL ){
        printf("cannot find %s\n", fname );
        return -1;
    }
    char buff[1024];
    char * line;
    int nl;
    line = fgets( buff, 1024, pFile );
    if(line==NULL){
        printf("read nl line NULL \n");
        fclose(pFile);
        return -1;
    }
    sscanf( buff,"%i", &nl );
    printf(" nl = %i \n", nl);
    allocateIOBuffs( nl, format, buffs );
    for(int il=0; il<nl; il++){
        line = fgets( buff, 1024, pFile );
        //while( line = fgets( buff, 1024, pFile ) ){
        //printf("%s \n", line );
        int ib = 0;
        const char *formati = format;
        char *tok = strtok(line, " \t");
        while (tok != NULL) {
            //my_array[i++] = atof(tok);
            //printf( "%s   %c \n", tok, *formati );
            switch( *formati ){
                case 'i': ((int   *)buffs[ib])[il]=atoi(tok); ib++; break;
                case 'f': ((float *)buffs[ib])[il]=atof(tok); ib++; break;
                case 'd': ((double*)buffs[ib])[il]=atof(tok); ib++; break;
                case '2':{
                    Vec2d& v2 = ((Vec2d*)buffs[ib])[il];
                    v2.x = atof(tok); tok = strtok(NULL, " \t");
                    v2.y = atof(tok);
                    }; ib++; break;
                case '3':{
                    //printf(" Vdfsdfsdf455464 \n");
                    Vec3d& v3 = ((Vec3d*)buffs[ib])[il];
                    v3.x = atof(tok); tok = strtok(NULL, " \t"); // printf("tok_y %s \n",tok);
                    v3.y = atof(tok); tok = strtok(NULL, " \t"); // printf("tok_z %s \n",tok);
                    v3.z = atof(tok);
                    }; ib++; break;
            }
            tok = strtok(NULL, " \t");
            formati++;
        }
    }
    return nl;
}


#endif
