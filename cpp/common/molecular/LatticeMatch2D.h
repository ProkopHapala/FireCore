
#ifndef LatticeMatch2D_h
#define LatticeMatch2D_h

#include  <algorithm>
#include "datatypes.h"
#include "Vec2.h"

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

typedef bool (*CompLatmiss)(const Latmiss& a, const Latmiss& b);

bool isStame(const Latmiss& a, const Latmiss& b){
    if( (a.n == b.n) && ( sq(a.d-b.d)<1e-6 ) ) return true;
    return false;
}



struct IVal{
    int      i;
    double val;
};

//bool compare_Latmiss_alpha(const Latmiss& a, const Latmiss& b){ return a.alpha>b.alpha; }

class LatticeMatch2D{ public:
    Vec2d lat0[2];
    Vec2d lat1[2];
    double angUV = 0;
    double angU  = 0;
    double angV  = 0;
    double lu0;
    double lv0;

    double dmax=0.1;

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

    int checkLengthMatch( std::vector<Latmiss>& out, int ia, int ib, double l0, double rab, double alpha ){
        double d = rab/l0; int n=(int)d; d-=n; 
        if     (    d <dmax ){ out.push_back(Latmiss{n  ,ia,ib,d  ,alpha}); return n; }
        else if( (1-d)<dmax ){ out.push_back(Latmiss{n+1,ia,ib,d-1,alpha}); return n; }
        return -1;
    }

    void walk2D( double Rmax, double dmax_=0.1 ){
        dmax=dmax_;
        //printf( " walk2D lat0 [(%g,%g)  (%g,%g)] \n", lat0[0].x,lat0[0].y,  lat0[1].x,lat0[1].y );
        //printf( " walk2D lat1 [(%g,%g)  (%g,%g)] \n", lat1[0].x,lat1[0].y,  lat1[1].x,lat1[1].y );

        // Walk over all possible sites (n*a+m*b) of 2d latticle (a,b) and evaluate if lenght rab=|n*a+m*b| of this site is close to some multiple of lenght of the other lattice vectors (u,v), store combinations found into vector match_u, match_v
        //Vec2d a1=;
        //Vec2d b1=;
        lu0 = lat1[0].norm();
        lv0 = lat1[1].norm();

        angU = atan2(lat1[0].y,lat1[0].x);
        angV = atan2(lat1[1].y,lat1[1].x);
        angUV = lat1[1].angle( lat1[0] );

        //printf("angU  %g \n", angU);
        //printf("angV  %g \n", angV);
        //printf("angUV %g \n", angUV);

        //double la = lat1[0].norm();
        //double lb = lat1[1].norm();
        double la = lat0[0].x;
        double lb = lat0[1].y;

        double k = lat0[1].x / lat0[0].x;

        //printf( "DEBUG walk2D u(%g,%g) v(%g,%g) \n", lat1[0].x,lat1[0].y,   lat1[1].x,lat1[1].y );
        //printf( "DEBUG walk2D k %g lat0.a(%g,%g), lat0.b(%g,%g) \n", lat0[0].x,lat0[0].y,   lat0[1].x,lat0[1].y );

        match_u.clear();
        match_v.clear();

        int na=(int)(Rmax/la);
        int nb=(int)(Rmax/lb);
        //printf( "DEBUG walk2D na,nb(%i,%i) Rmax,la,lb(%g,%g,%g) \n", na,nb, Rmax, la, lb );
        for(int ib=-nb-2; ib<nb+2; ib++ ){
            Vec2d vb = lat0[1]*ib;
            int ia0 = (int)(-k*ib);
            for(int ia=ia0-na-1; ia<ia0+na+2; ia++ ){
                Vec2d p = vb + lat0[0]*ia;
                //glColor3f(0.5,0.5,0.5); Draw3D::drawPointCross( {p.x,p.y,0.0}, 0.1 );
                double rab = p.normalize();
                if(rab>Rmax)         continue;
                if((ia==0)&&(ib==0)) continue;
                double alpha = atan2(p.y,p.x);
                //printf("angle(%i,%i) %g  angUV %g \n", ib,ia, alpha, angUV );
                int nu = checkLengthMatch( match_u, ia,ib, lu0, rab, alpha       );  //if(nu>0)printf("u[%i]: n,ia,ib(%i|%i,%i) rab %g   p(%g,%g) \n",match_u.size()-1, nu,ia,ib, rab, p.x*rab,p.y*rab  );
                //int nv = checkLengthMatch( match_v, ia,ib, lv0, rab, alpha+angUV );  //if(nv>0)printf("v[%i]: n,ia,ib(%i|%i,%i) rab %g   p(%g,%g) \n",match_u.size()-1, nv,ia,ib, rab, p.x*rab,p.y*rab  );
            }
        }

        //printVecs( match_u, 'u' );
        //printVecs( match_u, 'v' );

        //for( int i=0;i<match_u.size();i++ ){ printf( "match_u[%i] n,a,b(%i|%i,%i) d %g[%] alpha %g[/pi] \n", i, match_u[i].n, match_u[i].ia, match_u[i].ib, match_u[i].d*100, match_u[i].alpha/M_PI ); }
        //for( int i=0;i<match_v.size();i++ ){ printf( "match_v[%i] n,a,b(%i|%i,%i) d %g[%] alpha %g[/pi] \n", i, match_v[i].n, match_v[i].ia, match_v[i].ib, match_v[i].d*100, match_v[i].alpha/M_PI ); }

        //int nu = ((int)Rmax/lu0); glColor3f(1.0,0.0,0.0); 
        //int nv = ((int)Rmax/lv0); glColor3f(0.0,0.0,1.0); 
        //glColor3f(1.0,0.9,0.8); for( const Latmiss& l: match_u ){ Vec2d p = lat0[0]*l.ia + lat0[1]*l.ib; Draw3D::drawPointCross( {p.x,p.y,0.0}, 0.2 ); } for(int i=0;i<(nu+1);i++){ Draw3D::drawCircleAxis(100,Vec3dZero, Vec3dX, Vec3dZ, lu0*i  ); }
        //glColor3f(0.8,0.9,1.0); for( const Latmiss& l: match_v ){ Vec2d p = lat0[0]*l.ia + lat0[1]*l.ib; Draw3D::drawPointCross( {p.x,p.y,0.0}, 0.2 ); } for(int i=0;i<(nv+1);i++){ Draw3D::drawCircleAxis(100,Vec3dZero, Vec3dX, Vec3dZ, lv0*i  ); }
        printf( " walk2D nu %i nv %i Rmax %g dmax %g \n", match_u.size(), match_v.size(), Rmax, dmax );
        printf( " walk2D lat0 [(%g,%g)  (%g,%g)] \n", lat0[0].x,lat0[0].y,  lat0[1].x,lat0[1].y );
        printf( " walk2D lat1 [(%g,%g)  (%g,%g)] \n", lat1[0].x,lat1[0].y,  lat1[1].x,lat1[1].y );

    };

    void printVecs( std::vector<Latmiss>& match, char c ){
        for( int i=0;i<match.size();i++ ){ printf( "match_%c[%i] n,a,b(%i|%i,%i) d %g[0.01] alpha %g[/pi] \n", c, i, match[i].n, match[i].ia, match[i].ib, match[i].d*100, match[i].alpha/M_PI ); }
    }

    void cleansRedudat( std::vector<Latmiss>& match, CompLatmiss func=isStame ){
        int n = match.size()-1;
        //printf( "cleansRedudat nu %i nv %i \n", match_u.size(), match_v.size() );
        //printf( "cleansRedudat() n=%i match.size()=%i \n", n, match.size() );
        int i =1;
        while(i<n){
            const Latmiss& a = match[i];
            bool found=false;
            //printf( "cleansRedudat[%i] a.n=%i \n", i, a.n );
            for(int j=0; j<i; j++){
                bool hit = func(a, match[j]);
                //printf( "hit[%i,%i]=%i  | i.n=%i j.n=%i \n", i,j, hit,  a.n, match[j].n );
                if( hit ){ // are the same
                    _swap( match[i], match[n] );
                    found=true;
                    break;
                }
            }
            //printf( "cleansRedudat[%i<%i] a.n=%i found=%i \n", i, n, a.n, found );
            if(found){
                n--;
            }else{
                i++;
            }
        }
        match.resize(n);
    }

    Vec2d reproduce_vec( const Latmiss& L, Vec2d v, double ang0 ){
        v.rotate(L.alpha-ang0); v.mul(L.n);
        return v;
    }

    Vec2d reproduce_grid( const Latmiss& L ){
        return lat0[0]*L.ia  + lat0[1]*L.ib;
    }

    void makeMatch( Vec2i m, Vec2d& u, Vec2d& v, bool bOnGrid=true ){
        const Latmiss& Lu = match_u[m.i];
        const Latmiss& Lv = match_v[m.j];
        if( bOnGrid ){
            u = reproduce_grid( Lu );
            v = reproduce_grid( Lv );
        }else{
            u = reproduce_vec( Lu, lat1[0], angU       );
            v = reproduce_vec( Lv, lat1[1], angV+angUV );
        }
    }

    void angleToRange( ){
        for(int i=0; i<match_v.size(); i++){ match_v[i].alpha = dangle(  match_v[i].alpha ); }
        for(int i=0; i<match_u.size(); i++){ match_u[i].alpha = dangle(  match_u[i].alpha ); }
    }

    void sort(){
        std::sort( match_u.begin(), match_u.end(), [](Latmiss& a,Latmiss& b){return a.alpha<b.alpha;} );
        std::sort( match_v.begin(), match_v.end(), [](Latmiss& a,Latmiss& b){return a.alpha<b.alpha;} );
    }

    int matchAngles( double dangMax=0.1 ){
        if( (match_u.size()==0) || (match_v.size()==0) ){ printf( "WARRNING LatticeMatch2D::matchAngles() no matches found (match_size(%i,%i)) \n", match_u.size(), match_v.size() ); return 0; }
        //if( (match_u.size()==0) || (match_v.size()==0) ){ printf( "ERROR LatticeMatch2D::matchAngles() no matches found (match_size(%i,%i)) \n", match_u.size(), match_v.size() ); exit(0); }
        //printf( "DEBUG matchAngles() match_size(%i,%i) \n", match_u.size(), match_v.size() );
        matches.clear();
        //  Walk 1
        // assumption is that both arrays are ordered by angle in range < 0 .. 2pi >
        int     j  = 0;
        double  aj = match_v[j].alpha;
        double  a0 = 0;
        int nv = match_v.size();

        //printf( "match_u: n %i min %2.2f max %2.3f \n", match_u.size(),  match_u[0].alpha, match_u.back().alpha );
        //printf( "match_v: n %i min %2.2f max %2.3f \n", match_v.size(),  match_v[0].alpha, match_v.back().alpha );

        for(int i=0; i<match_u.size(); i++){ 
            //printf( "DEBUG matchAngles[%i] \n", i );
            Latmiss& mui = match_u[i];
            double amin=mui.alpha-dangMax;
            double amax=mui.alpha+dangMax;
            //printf( "### ang_u[%i]  arange(%2.2f,%2.2f) aj %2.4f \n", i, amin, amax, aj );
            while(  aj>amin ){ j--; if(j<0  ){j+=nv;a0-=2*M_PI;}; aj=match_v[j].alpha+a0; //printf( "%i ", j );
            //    printf( "ang_v[%i] %2.2f < %2.3f ? \n", j, amin,  aj );
            }
            bool b=true;
            while( b ){ j++; if(j>=nv){j-=nv;a0+=2*M_PI;}; aj=match_v[j].alpha+a0; //printf( "%i ", j );
                b = aj<amax;
                if( b && (aj>amin) ){
                    //printf( "ang_v[%i] %2.2f < %2.3f < %2.2f :-) \n", j, amin,  aj, amax );
                    matches.push_back( {i,j} );
                }else{
                    //printf( "ang_v[%i] %2.2f < %2.3f < %2.2f ? \n", j, amin,  aj, amax );
                }
            }
            //printf( "\n" );
        }
        return matches.size();
    }

    int exportVecMatch( std::vector<Latmiss>& match, int nmax, int2* inds=0, int* ns=0, double* errs=0, int* isorts=0, bool bSort=true, bool bPrint=false ){
        //printf( "exportVecMatch nu %i nv %i \n", match_u.size(), match_v.size() );
        //printf( "exportVecMatch() match.size()=%i bSort=%i \n", match.size(), bSort );
        int n   = match.size();
        if(nmax>n)nmax=n;
        //printf("DEBUG nmax %i \n", nmax);
        std::vector<IVal> ivals{n};
        if( (isorts!=0) || bSort){
            for(int i=0; i<n; i++){ ivals[i].i=i; ivals[i].val=sq(match[i].d); };
            std::sort( ivals.begin(),ivals.end() , [](IVal& a, IVal& b){ return a.val<b.val;} );
            if(isorts) for(int i=0; i<n; i++) isorts[i]=ivals[i].i;
        }
        for(int i=0; i<nmax; i++){
            int j=i;
            if(bSort) j=ivals[i].i;
            Latmiss& L = match[j];
            if(inds){ inds[i]=(int2){L.ia,L.ib}; }
            if(errs){ errs[i]=L.d; }
            if(ns  ){ ns  [i]=L.n; }
            if(bPrint){
                Vec2d u = reproduce_grid( L );
                printf( "[%i] (%i,%i|%i) err %g[0.01] L=%g[A] \n", i, inds[i].x, inds[i].y, ns[i], errs[i]*100.0, u.norm() );
            }
        }
        return nmax;
    }

    void exportMatch( int4* inds, int2* ns, double4* errs=0, int* isorts=0, bool bSort=true, double4 K=(double4){1.,1.,1.,0.} ){
        int n = matches.size();
        int4   * inds_ = inds;
        int2   * ns_   = ns;
        double4* errs_ = errs;
        if(bSort){
            if(inds)inds_= new int4   [n];
            if(ns  )ns_  = new int2   [n];
            if(errs)errs_= new double4[n];
        }
        for(int i=0; i<n; i++){
            Vec2i m = matches[i];
            const Latmiss& Lu = match_u[m.i];
            const Latmiss& Lv = match_v[m.j];
            if(inds){
                inds_[i].x = Lu.ia;
                inds_[i].y = Lu.ib;
                inds_[i].z = Lv.ia;
                inds_[i].w = Lv.ib;
            }
            if(ns){
                ns_  [i].x=Lu.n;
                ns_  [i].y=Lv.n;
            }
            //double ru0 = lat1[0].norm(); 
            //double ru0 = lat1[0].norm();
            if(errs){
                Vec2d u = reproduce_grid( Lu );
                Vec2d v = reproduce_grid( Lv );
                double lu = lu0*Lu.n;
                double lv = lv0*Lv.n;
                double ru  = u.norm();   errs_[i].x=(ru-lu)/lu;          // |u| cost
                double rv  = v.norm();   errs_[i].y=(rv-lv)/lv;          // |v| cost
                double ang = v.angle(u); errs_[i].z=(ang-angUV)/M_PI_2;  // angle distortion cost
                                         errs_[i].w=Lu.n*Lv.n;           // cell volume cost
                //printf(  "err[%i] %g %g %g %g  (%g,%i) (%g,%i) \n", i, errs_[i].x, errs_[i].y, errs_[i].z, errs_[i].w, lu0,Lu.n,  lv0,Lv.n );
            }
        }
        if( (isorts!=0) || bSort){
            std::vector<IVal> ivals{n};
            for(int i=0; i<n; i++){ ivals[i].i=i; double4 E=errs_[i]; ivals[i].val= sq(E.x*K.x)+sq(E.y*K.y)+sq(E.z*K.z)+sq(E.w*K.w); };
            std::sort( ivals.begin(),ivals.end() , [](IVal& a, IVal& b){ return a.val<b.val;} );
            //for(int i=0; i<n; i++){ printf("[%i] %i %g \n ", i, ivals[i].i, ivals[i].val ); }
            if(isorts) for(int i=0; i<n; i++) isorts[i]=ivals[i].i;
            //printf("DEBUG  4 \n");
            if(bSort){
                if(inds){ for(int i=0; i<n; i++){ inds[i]=inds_[ ivals[i].i ]; }  delete [] inds_; }
                if(errs){ for(int i=0; i<n; i++){ errs[i]=errs_[ ivals[i].i ]; }  delete [] errs_; }
                if(ns  ){ for(int i=0; i<n; i++){ ns  [i]=ns_  [ ivals[i].i ]; }  delete [] ns_;   }
            }
        }
        //printf("DEBUG  5 \n");
    }

    int match( double Rmax, double dmax=0.1, double dangMax=0.1 ){
        walk2D( Rmax, dmax=0.1 );
        angleToRange( );
        sort();
        matchAngles();
        return matches.size();
    }

};

#endif



