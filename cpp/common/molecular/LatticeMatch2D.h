
#ifndef LatticeMatch2D_h
#define LatticeMatch2D_h

#include "Vec2.h"
#include  <algorithm>

/*
void findLatticeR( double Rmin, double Rmax, Vec2d a, Vec2d b ){
    double R0 = (Rmin-Rmax)*0.5;
    double ra = a.norm();
    double rb = b.norm();
    double dna = R0/ra; int na=(int)(dna+0.5); dna-=na; 
    double dnb = R0/rb; int nb=(int)(dnb+0.5); dnb-=nb;
    double bonus = (a.x*b.x + a.y*b.y)*2;
}
*/

struct Latmiss{
    int n,ia,ib;
    double d;
    double alpha; 
};

//bool compare_Latmiss_alpha(const Latmiss& a, const Latmiss& b){ return a.alpha>b.alpha; }

class LatticeMatch2D{ public:
    Vec2d lat0[2];
    Vec2d lat1[2];

    std::vector<Latmiss> match_u;
    std::vector<Latmiss> match_v;
    std::vector<Vec2i>   matches;

    void normalizeLatticeRotation(){
        Vec2d rot = lat0[0];
        rot.normalize();
        lat0[0].udiv_cmplx(rot);
        lat0[1].udiv_cmplx(rot);
        lat1[0].udiv_cmplx(rot);
        lat1[1].udiv_cmplx(rot);
    }

    void walk2D( double Rmax, double dmax=0.1 ){
        // Walk over all possible sites (n*a+m*b) of 2d latticle (a,b) and evaluate if lenght rab=|n*a+m*b| of this site is close to some multiple of lenght of the other lattice vectors (u,v), store combinations found into vector match_u, match_v
        //Vec2d a1=;
        //Vec2d b1=;
        double ru = lat1[0].norm();
        double rv = lat1[1].norm();

        double angUV = lat1[1].angle( lat1[0] );

        //double la = lat1[0].norm();
        //double lb = lat1[1].norm();
        double la = lat0[0].x;
        double lb = lat0[1].y;

        double k = lat0[1].x / lat0[0].x;

        int na=(int)(Rmax/la);
        int nb=(int)(Rmax/lb);
        for(int ib=0; ib<nb+2; ib++ ){
            Vec2d vb = lat0[1]*ib;
            int ia0 = (int)(-k*ib);
            for(int ia=ia0-na-1; ia<ia0+na+2; ia++ ){
                Vec2d p = vb + lat0[0]*ia;

                glColor3f(0.5,0.5,0.5); Draw3D::drawPointCross( {p.x,p.y,0.0}, 0.1 );

                double rab = p.normalize();
                if(rab>Rmax) continue;
                double alpha = atan2(p.y,p.x);
                double du = rab/ru; int nu=(int)du; du-=nu; if(du<dmax){ match_u.push_back({nu,ia,ib,du,alpha      }); }else if( (1-du)<dmax ){ match_u.push_back({nu+1,ia,ib,du,alpha      }); };
                double dv = rab/rv; int nv=(int)dv; dv-=nv; if(dv<dmax){ match_v.push_back({nv,ia,ib,dv,alpha+angUV}); }else if( (1-du)<dmax ){ match_v.push_back({nv+1,ia,ib,dv,alpha+angUV}); };
                
            }
        }

        glColor3f(1.0,0.0,0.0); for( const Latmiss& l: match_u ){ Vec2d p = lat0[0]*l.ia + lat0[1]*l.ib; Draw3D::drawPointCross( {p.x,p.y,0.0}, 0.2 ); }
        glColor3f(0.0,0.0,1.0); for( const Latmiss& l: match_v ){ Vec2d p = lat0[0]*l.ia + lat0[1]*l.ib; Draw3D::drawPointCross( {p.x,p.y,0.0}, 0.2 ); }
    };



    void angleToRange( ){
        for(int i=0; i<match_v.size(); i++){ double a = match_v[i].alpha; if(a<0){ a+=2*M_PI; }else if(a>2*M_PI){ a=-2*M_PI; }; match_v[i].alpha = a; }
        for(int i=0; i<match_v.size(); i++){ double a = match_u[i].alpha; if(a<0){ a+=2*M_PI; }else if(a>2*M_PI){ a=-2*M_PI; }; match_u[i].alpha = a; }
    }

    void sort(){
        std::sort( match_u.begin(), match_u.end(), [](Latmiss& a,Latmiss& b){return a.alpha>b.alpha;} );
        std::sort( match_v.begin(), match_v.end(), [](Latmiss& a,Latmiss& b){return a.alpha>b.alpha;} );
    }

    void matchAngles( double dangMax=0.1 ){
        //  Walk 1
        // assumption is that both arrays are ordered by angle in range < 0 .. 2pi >
        int     j  = 0;
        double  aj = match_v[j].alpha;
        double  a0 = 0;
        int nv = match_v.size();
        for(int i=0; i<match_u.size(); i++){ 
            Latmiss& mui = match_u[i];
            double ang = mui.alpha;
            while(  (aj+a0)> (ang-dangMax) ){ j--; if(j<0  ){j+=nv;a0-=2*M_PI;}; aj=match_v[j].alpha; }
            while(  (aj+a0)< (ang+dangMax) ){ j++; if(j>=nv){j-=nv;a0+=2*M_PI;}; aj=match_v[j].alpha; 
                matches.push_back( {i,j} );
            }
        }
    }

    void match( double Rmax, double dmax=0.1, double dangMax=0.1 ){
        walk2D( Rmax, dmax=0.1 );
        angleToRange( );
        sort();
        matchAngles();
    }

};

#endif



