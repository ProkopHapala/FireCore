
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

class Groups{ public:
    std::vector<Group> groups;
    int* a2g=0;
    int* g2a=0;
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
        //printf( "Groups::evalRot()[%i] i0=%i n=%i \n", ig, g.i0n.x, g.i0n.y );
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

    void evalAllPoses(){
        printf( "Groups::evalAllPoses()\n" );
        for(int i=0; i<groups.size(); i++){
            //printf( "Groups::evalAllPoses()[%i]\n", i );
            evalCog( i );
            evalRot( i );
            const Group& g = groups[i];
            printf(  "group[%i] cog(%g,%g,%g) fw(%g,%g,%g) up(%g,%g,%g)\n", i,   g.cog.x,g.cog.y,g.cog.z,  g.fw.x,g.fw.y,g.fw.z, g.up.x,g.up.y,g.up.z );
        }
    }

    void initWeights(int natom){
        _realloc0( weights, natom, Quat4fOnes );
        _realloc0( fweights, natom, Vec2fOnes );
    }

    void bindAtoms(Vec3d* apos_, Vec3d* fapos_){ apos=apos_; fapos=fapos_; } 

    void setGroupMapping( int nAtoms, int nGroup, int* a2g_){
        printf( "Groups::setGroupMapping() nAtoms=%i nGroup=%i \n", nAtoms, nGroup );
        a2g     = a2g_;
        g2a     = new int [nAtoms];

        // --- map from a2g to g2a
        // -- count number of atoms in each group
        groups.resize(nGroup);
        for(int i=0; i<nGroup; i++){ groups[i].i0n=(Vec2i){0,0}; }  
        for(int i=0; i<nAtoms;  i++){ // count number of atoms in each group
            int ig = a2g[i];
            if(ig>=0){
                //printf( "ia(%i) -> ig(%i) \n", i, ig );
                groups[ig].i0n.y++;
            } 
        } 
        int nsum=0;
        for(int i=0; i<nGroup; i++){ int ni=groups[i].i0n.y; groups[i].i0n.x=nsum; nsum+=ni; } // accumulate group2atom index offsets
        //for(int i=0; i<nGroupTot; i++){ printf( "granges[%i] i0,n(%i,%i)\n", i, granges[i].x, granges[i].y  ); } 
        for(int i=0; i<nGroup; i++){ groups[i].i0n.y=0; }  // clean atoms in group count
        // -- back mapping
        for(int i=0; i<nAtoms; i++){
            //printf( "[%i] ", i );
            int  ig =  a2g[i];
            if(ig<0) continue;
            Vec2i& g = groups[ig].i0n;
            int j   = g.x + g.y;
            //printf( " i=%i ig=%i j=%i  | nAtomTot %i \n", i, ig, j, nAtomTot );
            g2a[j] = i;
            g.y++;
        }
        //for(int i=0; i<nGroupTot; i++){ granges[i].y=granges[i].y+granges[i].x; }
        // --- print Group->Atom mapping
        for(int ig=0; ig<nGroup; ig++){ 
            Vec2i ni  = groups[ig].i0n;
            printf( "Groups::setGroupMapping() group[%i] grange i0=%i n=%i \n", ig, ni.x, ni.y );
            for(int i=0; i<ni.y; i++){
                int ia = g2a[i + ni.x];
                printf( "group[%i][%i] ia=%i \n", ig, i, ia );
            }
        }
        //for(int i=0; i<nAtoms; i++){    printf( "atom[%i] -> group # %i \n", i, a2g_[i] );}
        //exit(0);
    }

};

#endif
