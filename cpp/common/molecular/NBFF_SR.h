/*
Non-Bonded Force-Field
 should be easily plugged with any molecular dynamics, by either sharing pointer to same data buffer, or by copying data
*/

#ifndef NBFF_SR_h
#define NBFF_SR_h

#include "fastmath.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Atoms.h"


#include "Forces.h"


int breadthFirst( int nstart, int* istarts, const int* neighs, int* vals, std::vector<int>* found, int ifill, int ifree=-1, int neighmax=4, int nfront0=100, int nmaxiter=100 ){
    std::unordered_set<int> front1,front2;
    std::unordered_set<int> *pfrom=&front1, *pto=&front2;
    //printf("breadthFirst(): front.size(%i)\n", front.size(%) );
    for(int i=0; i<nstart; i++){ int ii=istarts[i]; front1.insert(ii); vals[ii]=ifill; }
    int n=0;
    for( int itr=0; itr<nmaxiter; itr++ ){
        for( int i: *pfrom ){
            //vals[i]=ifill;
            const int* ngs = neighs[i*neighmax];
            for(int j=0; j<neighs; j++){
                ing = ngs[j];
                if(ing<0) continue;
                if(vals[ing]==ifree){
                    vals[ing]==ifill;
                    pto.insert(ing);
                    if(found)found.push_back(ing);
                    n++;
                }
            }
        }
        if(pto.size()==0) break;
        pfrom.clear();
        _swap( pfrom,pto );
    }
    return n;
}



class NBFF_SR: public NBFF{ public:
    // int     n      =0;        // from Atoms
    // int    *atypes =0;        // from Atoms
    // Vec3d  *apos   =0;        // from Atoms
    // Vec3d    *fapos  =0;      // from NBFF
    // Quat4d   *REQs   =0;      // from NBFF
    // Quat4i   *neighs =0;      // from NBFF
    // Quat4i   *neighCell=0;    // from NBFF
    // double alphaMorse = 1.5;  // from NBFF
    // //double  KMorse  = 1.5;  // from NBFF
    // double  Rdamp     = 1.0;  // from NBFF
    // Mat3d   lvec;             // from NBFF
    // Vec3i   nPBC;             // from NBFF
    // bool    bPBC=false;       // from NBFF
    // int    npbc   =0;         // from NBFF
    // Vec3d* shifts =0;         // from NBFF
    // Quat4f *PLQs  =0;         // from NBFF
    // Vec3d  shift0 =Vec3dZero; // from NBFF

    bool bInGroupNB=false;
    int     ngroups     =0;
    Vec2i*  groupin     =0;
    Quat4f* groupPR     =0;
    int*    a2g=0;
    int*    g2a=0;

    void reallocGroups(int ngroups_){
        ngroups=ngroups_;
        _realloc(groups);
    }

    void segmetToGroups( int ngroups, int* istarts, bool bAlloc=true ){
        std::vector<int> found{100};
        int i=0;
        for(int ig=0; ig<ngroups; ig++){
            int istart=istarts[ig];
            found.clear();
            int nf = breadthFirst( 1, &istart, neighs, a2g, &found, int ifill,-1,4);
            groupin[ig] = Vec2i{ i,nf };
            Vec3d pmin=Vec3dmax;
            Vec3d pmax=Vec3dmin;
            for(int ia:found){  g2a[i]=ia; i++; };
        }
    }

    void  updateGroupGeom(int ig){
        Vec2i gin = groupin[ig];
        Vec3d pmin=Vec3dmax;
        Vec3d pmax=Vec3dmin;
        for(int i=gin.a; i<gin.b; i++){
            int ia=g2a[i]; 
            apos[ia].update_bounds(pmin,pmax);
        }
        Vec3d cog = (pmin+pmax)*0.5;
        double r2max=0;
        for(int i=gin.a; i<gin.b; i++){ 
            int ia=g2a[i]; 
            Vec3d d=apos[ia]-cog; 
            r2max=fmax(r2max,d.norm2()); 
        };
        groupPR[ig].f = (Vec3f)cog;
        groupPR[ig].r = (float)sqrt(r2max); 
    }
    void updateGroupsGeom(){ for(int ig=0; ig<ngroups; ig++){ updateGroupGeom(int ig); } }

    double evalInGroupSR( Vec2i gin ){
        double E=0;

    }

    double evalInterGroupSR( Vec2i gin, Quat4f Gi, Vec2i gjn, Quat4f Gj ){
        double E=0;
        for(int i=gin.a; i<gin.b; i++){
            int ia=g2a[i];
            Vec3d pi=apos[ia];
            Vec3f dg = (Vec3f)pi -Gi.f;
            if(dg.norm2()>(Gi.w*Gi.w) )continue;
            for(int j=gjn.a; j<gjn.b; j++){
                int ja=g2a[i];
                Vec3d d = apos[ja]-pi;
                double r2 = d.norm2();
                if(r2<R2cut){
                    E+=getNBSR( REQi, R );
                }
            }
        }
    }

    double evalNonBondSR_groups(){
        double E=0;
        for(int ig=0; ig<ngroups; ig++){
            Vec2i   gin = groupin[ig];
            Quat4f  Gi  = groupPR[ig];
            // self interactin within group
            if(bInGroupNB)E+ = evalInGroupsSR( gin );  // non-covalent interactions inside groups can be often set to zero
            for(int jg=ig+1; jg<=ngroups; jg++){
                cost Quat4f&  Gj  = groupPR[jg];
                Vec3f    dg = Gj - Gi;
                float r2g = dg.norm2();
                float Rij = Gi.w + Gj.w;
                if( r2g>(Rij*Rij) ) continue;
                E+ = evalInterGroupsSR( gin,Gi, groupin[ig],Gj );
            }   
        }
    }


};

#endif

