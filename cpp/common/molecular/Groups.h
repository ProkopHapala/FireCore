
#ifndef Groups_h
#define Groups_h

#include <vector>

#include "Vec3.h"
#include "quaternion.h"
#include "Mat3.h"

#include "Draw3D.h" // debug

struct Group{
    Vec2i i0n;
    Vec3d cog;
    Vec3d fw;
    Vec3d up;
    
    Vec3d force=Vec3dZero;
    Vec3d torq =Vec3dZero;

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
            int ia = g2a[i + g.i0n.x];
            double w = weights[ia].x;
            cog.add_mul( apos[ia], w );
            wsum+=w;
        }        
        cog.mul(1/wsum);  // we assume weights are already normalized
        g.cog = cog;
        return cog;
    }

    void evalRot( int ig, bool bOrthonNorm=true ){
        Group& g = groups[ig];
        Vec3d fw = Vec3dZero;
        Vec3d up = Vec3dZero;
        //double fsum=0;
        //double usum=0;
        //printf( "Groups::evalRot()[%i] i0=%i n=%i \n", ig, g.i0n.x, g.i0n.y );
        for(int i=0; i<g.i0n.y; i++){
            int ia = g2a[i + g.i0n.x];
            const Quat4f w = weights[ia];
            const Vec3d& p = apos[ia];
            fw.add_mul( p, w.y );
            up.add_mul( p, w.z );
            //wsum+=w;
        }        
        //cog.mul(1/wsum);  // we assume weights are already normalized
        if(bOrthonNorm){
            fw.normalize();
            up.makeOrthoU(fw);
            up.normalize();
        }
        g.up = up;
        g.fw = fw;
    }

    void applyForce( int ig, Vec3d force, bool bLocal=true ){
        const Group& g = groups[ig];
        for(int i=0; i<g.i0n.y; i++){
            int ia = g2a[i + g.i0n.x];
            fapos[ia].add_mul( force,  fweights[ia].x  );
        }
    }

    void applyTorq( int ig, Vec3d torq, bool bLocal=true ){
        const Group& g = groups[ig];
        Vec3d tq;
        if(bLocal){
            Mat3d M  = g.rotMat();
            M.dot_to_T( torq, tq );
            {
                glColor3f( 1.0f,0.0f,0.0f ); Draw3D::drawVecInPos( M.a, g.cog );
                glColor3f( 0.0f,1.0f,0.0f ); Draw3D::drawVecInPos( M.b, g.cog );
                glColor3f( 0.0f,0.0f,1.0f ); Draw3D::drawVecInPos( M.c, g.cog );
                glColor3f( 1.0f,1.0f,1.0f ); Draw3D::drawVecInPos( tq*10.0, g.cog );
            }
        }else{
            tq=torq;
        }
        Vec3d fcog = Vec3dZero;
        double wsum =0;
        for(int i=0; i<g.i0n.y; i++){
            int ia  = g2a[i + g.i0n.x];
            Vec3d dp = apos[ia] - g.cog;
            Vec3d ft = cross( tq, dp );
            double w = fweights[ia].y; 
            wsum+=w;
            ft.mul(w);
            fcog.add(ft);
            fapos[ia].add( ft );
            {
                glColor3f( 1.0f,0.0f,0.0f );
                Draw3D::drawVecInPos( ft*10.0, apos[ia] );
            }
        }
        fcog.mul(-1/wsum);
        for(int i=0; i<g.i0n.y; i++){
            int ia  = g2a[i + g.i0n.x];
            fapos[ia].add_mul( fcog, fweights[ia].y );
            // {
            //     glColor3f( 1.0f,0.0f,0.0f );
            //     Draw3D::drawVecInPos( fapos[ia]*10.0, apos[ia] );
            // }
        }
    }

    void applyAllForces( double fsc=1.0, double tsc=1.0, bool bLocalTq=true, bool bLocalF=false ){
        for(int ig=0; ig<groups.size(); ig++){
            const Group& g = groups[ig];
            applyForce( ig, g.force*fsc, bLocalF );
            applyTorq ( ig, g.torq*tsc,  bLocalTq );
        }
    }

    void forceAtom(int ia){
        int ig = a2g[ia];
        if(ig>=0){
            const Group& g = groups[ig];
            const Vec3d dp = apos[ia] - g.cog;
            const Mat3d rot = g.rotMat();
            Vec3d tq; rot.dot_to_T( g.torq, tq ); 
            const Vec3d ft  = cross( tq, dp );
            //fapos[ia].add( ft + g.force );
            {
                glColor3f( 1.0f,0.0f,0.0f );
                Draw3D::drawVecInPos( ft, apos[ia] );
            }
        }
    }

    void defPoseByAtoms( int ig, Vec2i ifw=Vec2i{-1,-1}, Vec2i iup=Vec2i{-1,-1}, int i0=-1 ){
        printf( "Groups::defPoseByAtoms()\n" );
        Group& g = groups[ig];
        // --- cog
        if(i0>=0){
            if( a2g[i0]!=ig ){ printf("ERROR in Groups::defPoseByAtoms() COG atom(i0=%i) is not member of group ig(%i) => exit(0) \n", i0, ig ); exit(0); }
            for(int i=0; i<g.i0n.y; i++){
                int ia  = g2a[i + g.i0n.x];
                weights[ia].x=0;
            }
            weights[i0].x=1;
        }else{
            double wsum=0;
            for(int i=0; i<g.i0n.y; i++){
                int ia  = g2a[i + g.i0n.x];
                wsum   += weights[ia].x;
            }
            double renorm =1./wsum;
            for(int i=0; i<g.i0n.y; i++){
                int ia  = g2a[i + g.i0n.x];
                weights[ia].x *= renorm;
            }
        }
        evalCog( ig );
        //printf( "Groups::defPoseByAtoms() cog(%g,%g,%g) \n", g.cog.x,g.cog.y,g.cog.z );
        // --- fw
        Vec3d fw;
        if( ifw.x>=0 ){
            fw = apos[ifw.x] - apos[ifw.y];
            fw.normalize();
            g.fw = fw;
            //printf( "Groups::defPoseByAtoms() fw(%g,%g,%g) \n", fw.x,fw.y,fw.z );
            double csum=0;
            for(int i=0; i<g.i0n.y; i++){
                int ia   = g2a[i + g.i0n.x];
                Vec3d  d = apos[ia] - g.cog;
                double c = fw.dot( d );
                csum += c;
                weights[ia].y = c;
            }
            csum/=g.i0n.y;
            double csum2=0;
            for(int i=0; i<g.i0n.y; i++){
                int ia   = g2a[i + g.i0n.x];
                weights[ia].y -= csum;
                csum2 += weights[ia].y;
            }
            //printf( "fw csum2 = %g \n", csum2 );
        }
        // --- up
        Vec3d up;
        if( iup.x>=0 ){
            up = apos[iup.x] - apos[iup.y];
            up.makeOrthoU(fw);
            up.normalize();
            g.up = up;
            //printf( "Groups::defPoseByAtoms() up(%g,%g,%g) \n", up.x,up.y,up.z );
            double csum=0;
            for(int i=0; i<g.i0n.y; i++){
                int ia   = g2a[i + g.i0n.x];
                Vec3d  d = apos[ia] - g.cog;
                double c = up.dot( d );
                csum += c;
                weights[ia].z = c;
            }
            csum/=g.i0n.y;
            double csum2=0;
            for(int i=0; i<g.i0n.y; i++){
                int ia   = g2a[i + g.i0n.x];
                weights[ia].z -= csum;
                csum2 += weights[ia].z;
            }
            //printf( "up csum2 = %g \n", csum2 );
        }

        { // correct linear combination
            evalRot( ig, false );
            //crammer();
            Vec2d f{ fw.dot( g.fw ), up.dot( g.fw ) };
            Vec2d u{ fw.dot( g.up ), up.dot( g.up ) };

            double af = 1/(f.x - f.y*u.x/u.y);
            Vec2d cf{ af, -(f.y/u.y) * af };
            
            //Vec2d cu{ -cf.y, cf.x };

            double au = 1/(f.y - f.x*u.y/u.x);
            Vec2d  cu{ au, -(f.x/u.x) * au };

            Vec3d fw_ = g.fw *  cf.x  +    g.up * cf.y;
            Vec3d up_ = g.fw *  cu.x  +    g.up * cu.y;

            if( up.dot(up_)<0 ){ cu.mul(-1.0);  up_.mul(-1); }

            printf(  "|fw|=%g |up|=%g <fw|up>=%g \n", fw.norm(), up.norm(), fw.dot(up) ); 
            printf(  "<fw|fw_>=%g fw(%g,%g,%g)  fw_(%g,%g,%g)  cf(%g,%g)\n",  fw.dot(fw_), fw.x,fw.y,fw.z, fw_.x,fw_.y,fw_.z,  cf.x, cf.y );  
            printf(  "<up|up_>=%g up(%g,%g,%g)  up_(%g,%g,%g)  cu(%g,%g)\n",  up.dot(up_), up.x,up.y,up.z, up_.x,up_.y,up_.z,  cu.x, cu.y );  

            //Vec2d ufw crammer( , {1.0,0.0} );
            //Vec2d cup{ up.dot( g.fw ), up.dot( g.up ) };
            Vec2f ws[g.i0n.y];
            for(int i=0; i<g.i0n.y; i++){
                int ia = g2a[ i + g.i0n.x ];
                ws[i].x=weights[ia].y;
                ws[i].y=weights[ia].z;
            }
            for(int i=0; i<g.i0n.y; i++){
                int ia = g2a[ i + g.i0n.x ];
                Vec2f w = ws[i];
                weights[ia].y = w.x*cf.x + w.y*cf.y;
                weights[ia].z = w.x*cu.x + w.y*cu.y;
            }
            evalRot( ig, false );
        }

    }

    void evalAllPoses(){
        //printf( "Groups::evalAllPoses()\n" );
        for(int i=0; i<groups.size(); i++){
            //printf( "Groups::evalAllPoses()[%i]\n", i );
            evalCog( i );
            evalRot( i );
            const Group& g = groups[i];
            //printf(  "group[%i] cog(%g,%g,%g) fw(%g,%g,%g) up(%g,%g,%g)\n", i,   g.cog.x,g.cog.y,g.cog.z,  g.fw.x,g.fw.y,g.fw.z, g.up.x,g.up.y,g.up.z );
        }
    }

    void realloc(int natom){
        printf( "Groups::realloc(natom=%i)\n", natom );
        _realloc0( weights,  natom, Quat4fZero );
        _realloc0( fweights, natom, Vec2fZero );
        _realloc0( a2g,      natom, -1 );
        _realloc0( g2a,      natom, -1 );
    }

    void fill( int ig, int n, int* ias ){
        printf( "Groups::fill(ig=%i,n=%i)\n", ig, n );
        Group& g = groups[ig];
        int i0=0;
        if(ig>0){ i0 = groups[ig-1].i0n.x+groups[ig-1].i0n.y; }
        g.i0n.y = n;
        g.i0n.x = i0;
        for(int i=0; i<n; i++){
            int ia  = ias[i];
            g2a[i+i0] = ia;
            a2g[ia  ] = ig; 
        }
    }

    int addGroup( int n, int* ias ){
        int ig = groups.size();
        groups.resize(ig+1);
        fill( ig, n, ias );
        return ig;
    }

    void bindAtoms(Vec3d* apos_, Vec3d* fapos_){ apos=apos_; fapos=fapos_; } 

    void setGroupMapping( int nAtoms, int nGroup, int* a2g_, bool bResize=true ){
        printf( "Groups::setGroupMapping() nAtoms=%i nGroup=%i \n", nAtoms, nGroup );
        //a2g     = a2g_;
        //g2a     = new int [nAtoms];

        // --- map from a2g to g2a
        // -- count number of atoms in each group
        if(bResize)groups.resize(nGroup);
        for(int i=0; i<nGroup; i++){ groups[i].i0n=(Vec2i){0,0}; }  
        for(int i=0; i<nAtoms; i++){ // count number of atoms in each group
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


    void print_groups2atoms(){
        printf( "Groups::print_groups2atoms()\n" );
        for(int ig=0; ig<groups.size(); ig++){
            Group& g = groups[ig];
            printf( "group[%i](i0=%i,n=%i) ", ig, g.i0n.x, g.i0n.y  );
            for(int i=0; i<g.i0n.y; i++){
                int ia = g2a[i+g.i0n.x];
                printf( "%i,", ia );
            }
            printf( "\n" );
        }
    }

};

#endif
