
#ifndef datatypes_utils_h
#define datatypes_utils_h

#include "datatypes.h"
#include "Vec3.h"
#include "quaternion.h"

inline int add_except(int i,int di,int ino){ return (i==ino)?i:i+di; }
inline Quat4i add_except(Quat4i q, int di, int ino=-1){
    return (Quat4i){ 
        add_except(q.x,di,ino), 
        add_except(q.y,di,ino),
        add_except(q.z,di,ino),
        add_except(q.w,di,ino)
    };
};

void     copy_add(int n, Quat4i* from, Quat4i* to, int i0, int ino=-1 ){ for(int i=0; i<n; i++){ to[i]= add_except( from[i],i0,ino); } }

void     copy    (int n, Quat4i* from, Quat4i* to){ for(int i=0; i<n; i++){ to[i]=from[i]; } }
void     copy    (int n, Quat4f* from, Quat4f* to){ for(int i=0; i<n; i++){ to[i]=from[i]; } }

void     set     (int n, Quat4f* qs, Quat4f v=Quat4fZero ){ for(int i=0; i<n; i++){ qs[i]=v; } }

void     pack    (int n, Quat4d* fs, Quat4f* qs, float K=0 ){ for(int i=0; i<n; i++){ qs[i]  =(Quat4f)fs[i];            } }
void     pack    (int n, Vec3d*  fs, Quat4f* qs, float K=0 ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];  qs[i].e=K; } }
void   unpack    (int n, Vec3d*  fs, Quat4f* qs            ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f;           } }
double unpack_add(int n, Vec3d*  fs, Quat4f* qs            ){ double E=0; for(int i=0; i<n; i++){ fs[i].add( (Vec3d)qs[i].f ); E+=qs[i].e; }; return E; }

void pack(int n, Vec3d* fs, double* es, Quat4f* qs ){
    for(int i=0; i<n; i++){ 
        float e; if(es){e=es[i];}else{e=0;}
        //float e; if(es){e=es[i];}else{e=0;}
        qs[i].f=(Vec3f)fs[i];
        qs[i].e=e; 
    } 
}
Quat4f* pack(int n, Vec3d* fs, double* es=0){
    Quat4f* qs = new Quat4f[n];
    pack( n, fs, es, qs );
    return qs;
}

void unpack(int n, Vec3d* fs, double* es, Quat4f* qs ){
    for(int i=0; i<n; i++){ 
        const Quat4f& q= qs[i];
        fs[i]      =(Vec3d)q.f;
        if(es)es[i]=       q.e; 
    } 
}

#endif
