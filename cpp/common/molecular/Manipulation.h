#ifndef Manipulation_h
#define Manipulation_h

#include <vector>
#include <stdio.h>
#include <math.h>

#include "Vec3.h"
#include "Bspline.h"

struct SplineConstraintManager{

    std::vector<Vec3d> cps;

    int    iAnchor  = -1;
    double Kanchor  = 10.0;

    bool   bActive  = false;

    // Parameterization along spline
    double t        = 0.0;   // [0,1]
    Vec3d  pos      = Vec3dZero;

    // direct manipulation of anchor target position (independent of spline)
    bool   bUseSpline = true;

    void clear(){ cps.clear(); iAnchor=-1; t=0.0; pos=Vec3dZero; bActive=false; bUseSpline=true; }

    void setAnchor(int ia){
        iAnchor=ia;
        printf("SplineConstraintManager::setAnchor(ia=%i)\n", ia);
    }

    void setT(double t_){
        t=t_;
        if(t<0.0)t=0.0; if(t>1.0)t=1.0;
        if(bUseSpline) pos = evalPos(t);
    }

    void moveT(double dt){ setT(t+dt); }

    // ======== spline evaluation (uniform cubic B-spline)

    inline Vec3d evalPos_cubic(int i0, double u)const{
        const Quat4d b = Bspline::basis(u);
        Vec3d p = Vec3dZero;
        p.add_mul(cps[i0  ], b.x);
        p.add_mul(cps[i0+1], b.y);
        p.add_mul(cps[i0+2], b.z);
        p.add_mul(cps[i0+3], b.w);
        return p;
    }

    Vec3d evalPos(double t_)const{
        const int n = (int)cps.size();
        if(n<=0) return pos;
        if(n==1) return cps[0];
        if(n<4){
            // linear polyline fallback
            double s = t_*(n-1);
            int i = (int)s; if(i<0)i=0; if(i>n-2)i=n-2;
            double u = s - i;
            return cps[i]*(1-u) + cps[i+1]*u;
        }
        double s = t_*(n-3);            // segments [0..n-4]
        int i0   = (int)s; if(i0<0)i0=0; if(i0>n-4)i0=n-4;
        double u = s - i0;
        return evalPos_cubic(i0,u);
    }

    void updatePos(){ if(bUseSpline){ pos = evalPos(t); } }

    // ======== I/O

    bool loadDat(const char* fname){
        FILE* f = fopen(fname,"r");
        if(!f){ printf("ERROR: SplineConstraintManager::loadDat(%s) cannot open => exit()\n", fname); exit(0); }
        cps.clear();
        char line[1024];
        int il=0;
        while(fgets(line,1024,f)){
            il++;
            if(line[0]=='#') continue;
            Vec3d p;
            int nret = sscanf(line,"%lf %lf %lf", &p.x,&p.y,&p.z);
            if(nret==3){ cps.push_back(p); }
        }
        fclose(f);
        printf("SplineConstraintManager::loadDat(%s) loaded %zu points\n", fname, cps.size());
        bUseSpline = true;  // ensure we follow spline after loading
        updatePos();
        return cps.size()>0;
    }

    // ======== forces

    inline void addAnchorForce(int ia, const Vec3d& ap, Vec3d& fa)const{
        if(!bActive) return;
        if(ia!=iAnchor) return;
        const Vec3d dp = pos - ap;
        fa.add_mul(dp, Kanchor);
    }

};

#endif
