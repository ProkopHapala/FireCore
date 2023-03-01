﻿
#ifndef datatypes_utils_h
#define datatypes_utils_h

#include "datatypes.h"
#include "Vec3.h"
#include "quaternion.h"

void     pack    (int n, Vec3d* fs, Quat4f* qs, float K=0 ){ for(int i=0; i<n; i++){ qs[i].f=(Vec3f)fs[i];  qs[i].e=K; } }
void   unpack    (int n, Vec3d* fs, Quat4f* qs            ){ for(int i=0; i<n; i++){ fs[i]  =(Vec3d)qs[i].f;           } }
double unpack_add(int n, Vec3d* fs, Quat4f* qs            ){ double E=0; for(int i=0; i<n; i++){ fs[i].add( (Vec3d)qs[i].f ); E+=qs[i].e; }; return E; }

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