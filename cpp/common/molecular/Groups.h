
#ifndef Groups_h
#define Groups_h

#include <vector>

#include "Vec3.h"
#include "quaternion.h"
#include "Mat3.h"

struct Group{
    Vec2i i0n;
    Vec3d cog;
    Vec3d fw;
    Vec3d up;

    Mat3d rotMat()const{
        Vec3d lf = cross(fw,up);
        lf.normalize();
        return Mat3d{ fw.x,fw.y,fw.z, up.x,up.y,up.z, lf.x,lf.y,lf.z };
    }
};

class Groups{
    std::vector<Group> groups;
    Vec3d*   apos    = 0;
    Vec3d*  fapos    = 0;
    Quat4f* weights  = 0; // {wcog,wfw,wup}
    Vec2f*  fweights = 0; // {wpush, wtorq}

    Vec3d evalCog( int ig ){
        Group& g = groups[ig];
        Vec3d cog=Vec3dZero;
        double wsum=0;
        for(int i=0; i<g.i0n.y; i++){
            int ia = i + g.i0n.x;
            double w = weights[ia].x;
            cog.add_mul( apos[ia], w );
            //wsum+=w;
        }        
        //cog.mul(1/wsum);  // we assume weights are already normalized
        g.cog = cog;
        return cog;
    }

    void evalRot( int ig ){
        Group& g = groups[ig];
        Vec3d fw=Vec3dZero;
        Vec3d up=Vec3dZero;
        double fsum=0;
        double usum=0;
        for(int i=0; i<g.i0n.y; i++){
            int ia = i + g.i0n.x;
            double wf = weights[ia].y;
            double wu = weights[ia].z;
            const Vec3d& p = apos[ia];
            fw.add_mul( p, wf );
            up.add_mul( p, wu );
            //wsum+=w;
        }        
        //cog.mul(1/wsum);  // we assume weights are already normalized
        fw.normalize();
        up.makeOrthoU(fw);
        up.normalized();
        g.up = up;
        g.fw = fw;
    }

    void applyForce( int ig, Vec3d force, bool bLocal=true ){
        const Group& g = groups[ig];
        for(int i=0; i<g.i0n.y; i++){
            int ia = i + g.i0n.x;
            fapos[ia].add_mul( force,  fweights[ia].x  );
        }
    }

    void applyTorq( int ig, Vec3d dir, double scale, bool bLocal=true ){
        const Group& g = groups[ig];
        Mat3d M  = g.rotMat();
        Vec3d tq;  
        M.dot_to_T( dir, tq );
        for(int i=0; i<g.i0n.y; i++){
            int ia  = i + g.i0n.x;
            Vec3d dp = apos[ia] - g.cog;
            Vec3d ft = cross( tq, dp );
            fapos[ia].add_mul( ft,  fweights[ia].y*scale  );
        }
    }

};

#endif
