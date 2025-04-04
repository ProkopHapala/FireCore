
#ifndef  testUtils_h
#define  testUtils_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"







inline void _endl(){ printf("\n"); };

//#define print(A)    {A.print(); }
//#define println(A)  {A.print(); puts("");}

inline void printArray( int n, double * vec ){	for (int i=0; i<n; i++){ printf( " %f ", vec[i] ); }; printf("\n"); };
inline void printArray( int n, float  * vec ){	for (int i=0; i<n; i++){ printf( " %f ", vec[i] ); }; printf("\n"); };
inline void printArray( int n, int    * vec ){	for (int i=0; i<n; i++){ printf( " %i ", vec[i] ); }; printf("\n"); };

inline void printVec( const Vec2d& v ){ printf( " %f %f \n", v.x, v.y ); };
inline void printVec( const Vec2f& v ){ printf( " %f %f \n", v.x, v.y );};
inline void printVec( const Vec2i& v ){ printf( " %i %i \n", v.x, v.y );};

inline void printVec( const Vec3d& v ){ printf( " %f %f %f \n", v.x, v.y, v.z );};
inline void printVec( const Vec3f& v ){ printf( " %f %f %f \n", v.x, v.y, v.z );};
inline void printVec( const Vec3i& v ){ printf( " %i %i %i \n", v.x, v.y, v.z );};

inline void printQuat( const Quat4d& q ){	printf( " %f %f %f %f \n", q.x, q.y, q.z, q.w ); }
inline void printQuat( const Quat4f& q ){	printf( " %f %f %f %f \n", q.x, q.y, q.z, q.w ); }
inline void printQuat( const Quat4i& q ){ 	printf( " %i %i %i %i \n", q.x, q.y, q.z, q.w ); }

inline void printVecs( int n, const Vec3d * vs ){ 	for(int i=0; i<n;i++)printf( "[%i] (%g,%g,%g)\n",   i, vs[i].x, vs[i].y, vs[i].z ); }
inline void printVecs( int n, const Vec3f * vs ){ 	for(int i=0; i<n;i++)printf( "[%i] (%g,%g,%g)\n",   i, vs[i].x, vs[i].y, vs[i].z ); }
inline void printVecs( int n, const Quat4f* vs ){ 	for(int i=0; i<n;i++)printf( "[%i] (%g,%g,%g,%g)\n",i, vs[i].x, vs[i].y, vs[i].z, vs[i].w ); }

// inline void printMat( const Mat3d& mat  ){
// 	printf( " %f %f %f \n", mat.ax, mat.ay, mat.az );
// 	printf( " %f %f %f \n", mat.bx, mat.by, mat.bz );
// 	printf( " %f %f %f \n", mat.cx, mat.cy, mat.cz );
// }

template<typename T> inline void fprintVec(FILE* fout,int n      ,T* vec, const char* fmt="%20.10f "){ for(int i=0;i<n;i++){ fprintf( fout,fmt, vec[i] ); }; fprintf(fout,"\n"); }
template<typename T> inline void fprintMat(FILE* fout,int n,int m,T* mat, const char* fmt="%20.10f "){ for(int i=0;i<n;i++){ fprintVec(fout,m,mat+i*m, fmt); } }
template<typename T> inline void mat2file(const char* fname,int n, int m, T* vec, const char* fmt="%20.10f ", const char* mode="w" ){ FILE* fout=fopen(fname,mode); fprintMat(fout,n,m,vec,fmt); fclose(fout); }


inline bool compareVecs( int n, const Vec3d* vs, const Quat4f* qs, double Rmax, int verbosity=1 ){
    bool ret=false;
    double R2max=Rmax*Rmax;
    int imax=0;
    double r2max=0; 
    for(int i=0; i<n;i++){
        Vec3f d = ((Vec3f)vs[i]) - qs[i].f;
        double r2 = d.norm2();
        if(r2>r2max){ imax=i; r2max=r2; }
        r2/=( vs[i].norm2() + qs[i].f.norm2() + 1. );
        if(r2>R2max){
            if(verbosity>1)printf( "compareVecs[%i] v(%20.10f,%20.10f,%20.10f) q.f(%20.10f,%20.10f,%20.10f) |r| %f d(%20.10f,%20.10f,%20.10f) \n", i,  vs[i].x,vs[i].y,vs[i].z,   qs[i].x,qs[i].y,qs[i].z, sqrt(r2),  d.x,d.y,d.z   );
            ret=true;
        }
    } 
    if( ret && (verbosity>0)){
        int i=imax; 
        Vec3f  d = ((Vec3f)vs[i]) - qs[i].f;
        double r2 = d.norm2();
        r2/=( vs[i].norm2() + qs[i].f.norm2() + 1. );
        printf( "compareVecs[imax=%i] v(%20.10f,%20.10f,%20.10f) q.f(%20.10f,%20.10f,%20.10f) |r| %f d(%20.10f,%20.10f,%20.10f) \n", i,  vs[i].x,vs[i].y,vs[i].z,   qs[i].x,qs[i].y,qs[i].z, sqrt(r2),  d.x,d.y,d.z   );
    };
    return ret;
}

inline void printMat( const Mat3d& mat  ){
	printf( " %f %f %f \n", mat.ax, mat.ay, mat.az );
	printf( " %f %f %f \n", mat.bx, mat.by, mat.bz );
	printf( " %f %f %f \n", mat.cx, mat.cy, mat.cz );
}

_template_T void println(const T& t){t.print(); puts(""); };
_template_T void println(const char* s,const T& t){ printf("%s",s);t.print(); puts(""); };

// CPU ticks timer
// http://stackoverflow.com/questions/6432669/variance-in-rdtsc-overhead

#define DEBUG   printf( "DEBUG #l %i %s \n",    __LINE__, __FUNCTION__ );
#define DEBUGF  printf( "DEBUG #l %i %s %s \n", __LINE__, __FUNCTION__, __FILE__ );

//void dbg(char* s){ printf("DEBUG (%s) \n", s); };

#define DBG(format,args...) { printf("DEBUG "); printf(format,## args); }

/*
int  dbg(int priority, const char *format, ...){
    va_list args;
    va_start(args, format);
    if(priority & PRIO_LOG) vprintf(format, args);
    va_end(args);
}
*/


inline void save1Dslice( const char* fname, int n, Vec3i i0, Vec3i di, Vec3i ns, int m, double* data, const char* fmt="%g "  ){
    FILE* f=fopen(fname,"w");
    for(int i=0;i<n;i++){
        int ii = (i0.x+i*di.x) + ((i0.y+i*di.y) + (i0.z+i*di.z)*ns.x)*ns.y;
        for(int j=0;j<m;j++){
            fprintf(f,fmt, data[ii*m+j] );
        }
        fprintf(f,"\n");
    }
    fclose(f);
}



template<typename Func>
double checkDeriv(Func getEF, double x, double dx, double& fE, double& f ){
    double f1=0,f2=0;
    double e1 = getEF(x   ,f1);
    double e2 = getEF(x+dx,f2);
    f   = (f1+f2)*0.5;
    fE  = (e2-e1)/dx;
    double err = f - fE;
    printf( " |f-fE|: %g f %g fE %g E %g \n", err, f, fE, (e1+e2)*0.5 );
    return err;
}

template<typename Func>
double checkDeriv3d(Func getEF,const Vec3d p0, double d, Vec3d& fE, Vec3d& f ){
    fE=Vec3dZero;
    f =Vec3dZero;
    getEF(p0,f);
    for(int i=0;i<3;i++){
        double E0,E1;
        Vec3d p=p0;
        p .array[i]-=d;   E0=getEF(p,f);
        p .array[i]+=d*2; E1=getEF(p,f);
        p .array[i]-=d;
        fE.array[i]=(E1-E0)/(2*d);
    }
    double err = (f-fE).norm();
    printf( " |f-fE|: %g   fE(%g,%g,%g)   f(%g,%g,%g) \n", err, fE.x, fE.y, fE.z, f.x,f.y,f.z );
    return err;
}


inline uint64_t getCPUticks(){
    uint32_t lo, hi;
    __asm__ __volatile__ (
      "xorl %%eax, %%eax\n"
      "cpuid\n"
      "rdtsc\n"
      : "=a" (lo), "=d" (hi)
      :
      : "%ebx", "%ecx" );
    return (uint64_t)hi << 32 | lo;
}

class StopWatch{ public:
    long t1;
    double T;
    void   start(){ t1=getCPUticks(); };
    double stop (){ T=getCPUticks()-t1; return T; };
};

static StopWatch stopWatch;


inline void genRandomArray( int n, double * vals, double vmin, double vmax ){
    for( int i=0; i<n; i++ ){ vals[i] = randf( vmin, vmax ); }
}

inline void makeSamples2D( Vec2i ns, Vec3d p0, Vec3d da, Vec3d db, Vec3d *ps ){
    for(int ib=0; ib<ns.y; ib++){
        Vec3d p = p0+db*ib;
        for(int ia=0; ia<ns.x; ia++){
            *ps = p;
            p.add(da);
            ps++;
        }
    }
}

//template<void func(Vec3d pos, Vec3d& fout)>
template <typename Func>
void sampleVecField(Func func, int n, Vec3d *ps, Vec3d *fs ){
    for(int i=0; i<n; i++){ func(ps[i],fs[i]); }
}

//template<void func(Vec3d pos, Vec3d& fout)>
template <typename Func>
int sampleVecField(Func func, Vec2i ns, Vec3d p0, Vec3d a, Vec3d b, Vec3d*&ps, Vec3d*&fs ){
    int ntot = ns.x*ns.y;
    ps = new Vec3d[ntot];
    fs = new Vec3d[ntot];
    makeSamples2D(ns, p0, a, b, ps);
    //sampleFroce<func>( ntot, ps, fs, func );
    sampleForce( func, ntot, ps, fs  );
    return ntot;
}

//template<void func(Vec3d pos, Vec3d& fout)>
template <typename Func>
void sampleScalarField(Func func, int n, Vec3d *ps, double *Es, Vec2d& val_range ){
    val_range={+1e+300,-1e+300};
    for(int i=0; i<n; i++){ double Ei=func(ps[i]); Es[i]=Ei; val_range.enclose( Ei );
       // printf( "sampleScalarField [%i] p(%g,%g,%g) -> %g \n", i, ps[i].x,ps[i].y,ps[i].z,  Ei  );
    }
}

//template<void func(Vec3d pos, Vec3d& fout)>
template <typename Func>
int sampleScalarField(Func func, Vec2i ns, Vec3d p0, Vec3d a, Vec3d b, Vec3d*&ps, double*&Es, Vec2d& val_range ){
    int ntot = ns.x*ns.y;
    ps = new Vec3d[ntot];
    Es = new double[ntot];
    makeSamples2D(ns, p0, a, b, ps);
    //sampleFroce<func>( ntot, ps, fs, func );
    sampleScalarField( func, ntot, ps, Es,  val_range );
    //delete ps;
    return ntot;
}

#define STORE_ERROR(err)                    \
    {                                       \
        err_sum += err*err;                 \
        err_max = fmax(err_max,fabs(err));  \
    }

#define TEST_ERROR_PROC_N( caption, proc, n ) \
    do{                                    \
    double err_sum  = 0.0;                  \
    double err_max  = 0.0;                  \
    for( int i=0; i<n; i++ ){               \
        proc;                               \
    }                                       \
    printf( "%s MaxErr: %g RMSE: %g \n", caption, err_max, sqrt(err_sum/n) );   \
    } while (0) \



#define SPEED_TEST_FUNC( caption, func, xmin, xmax, ncall ) \
    do{                                    \
    double sum  = 0.0;                     \
    double dx   = (xmax-xmin)/(ncall-1);   \
    long tstart = getCPUticks();           \
    double x = xmin;                       \
    for( int i=0; i<ncall; i++ ){          \
        sum += func ( x );                 \
        x   += dx;                         \
        /*printf("%i %f %f %f \n", i, x, func ( x ), sum  );*/  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \


#define SPEED_TEST_FUNC_ARRAY( caption, func, arr, n, m ) \
    do{                                    \
    double sum  = 0.0;                     \
    long tstart = getCPUticks();           \
    for(int j=0; j<m; j++){                \
        for( int i=0; i<n; i++ ){          \
            sum += func ( arr[i] );        \
            /*printf("%i %f %f %f \n", i, arr[i], func ( arr[i] ), sum  );*/  \
        }                                  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    int ncall = n*m;                       \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \


#define SPEED_TEST_FUNC_NM( caption, func, n, m ) \
    do{                                    \
    double sum  = 0.0;                     \
    long tstart = getCPUticks();           \
    for(int j=0; j<m; j++){                \
        for( int i=0; i<n; i++ ){          \
            sum += func;                   \
        }                                  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    int ncall = n*m;                       \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \

#define SPEED_TEST_PROC_NM( caption, proc, n, m ) \
    do{                                    \
    double sum  = 0.0;                     \
    long tstart = getCPUticks();           \
    for(int j=0; j<m; j++){                \
        for( int i=0; i<n; i++ ){          \
            proc;                          \
        }                                  \
    }                                      \
    long time   = getCPUticks() - tstart;  \
    int ncall = n*m;                       \
    printf( "%s : %3.3f ticks/call ( %g %g ) | %g \n", caption, time/double(ncall), (double)time, (double)ncall, sum );   \
    } while (0) \

#endif




