#ifndef GridFF_h
#define GridFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Forces.h"
#include "MMFFparams.h"

static bool bDebug__ = 0;

class GridFF{ public: 
    GridShape   grid;
    //Vec3d  *FFPauli    = NULL;
    //Vec3d  *FFLondon   = NULL;
    //Vec3d  *FFelec     = NULL;

    Quat4f *FFPauli  = NULL;
    Quat4f *FFLondon = NULL;
    Quat4f *FFelec   = NULL;

    //Vec3d  *FFtot    = NULL; // total FF is not used since each atom-type has different linear combination

    int  natoms     = 0;
    int    * atypes = NULL;
    Vec3d  * apos   = NULL;   // atomic position
    Vec3d  * aREQs  = NULL;
    //Vec3d  * aPLQ = NULL;

    Vec3d shift = Vec3dZero;

    double alpha  = -1.5;
    double Rdamp  =  1.0;

    int iDebugEvalR = 0;
    bool bCellSet = false;

    double findTop(){ double zmax=-1e+300; for(int i=0;i<natoms; i++){ double z=apos[i].z; if(z>zmax)zmax=z; }; return zmax; }

    void bindSystem(int natoms_, int* atypes_, Vec3d* apos_, Vec3d* aREQs_ ){
        natoms=natoms_; atypes=atypes_; apos=apos_; aREQs=aREQs_;
    }

    void allocateFFs(){
        int ntot = grid.getNtot();
        //FFtot    = new Vec3d[ntot];
        //FFPauli  = new Vec3d[ntot];
        //FFLondon = new Vec3d[ntot];
        //FFelec   = new Vec3d[ntot];
        
        FFPauli  = new Quat4f[ntot];
        FFLondon = new Quat4f[ntot];
        FFelec   = new Quat4f[ntot];
    }
    
    void allocateAtoms(int natoms_){
        natoms = natoms_;
        //atypes = new int  [natoms];
        apos   = new Vec3d[natoms];
        aREQs  = new Vec3d[natoms];
    }

    int loadCell(const char * fname ){ return grid.loadCell(fname ); }

    /*
    inline void addForce( const Vec3d& pos, const Vec3d& PLQ, Vec3d& f ) const {
        Vec3d gpos;
        grid.cartesian2grid(pos, gpos);
        //printf( "pos: (%g,%g,%g) PLQ: (%g,%g,%g) \n", pos.x, pos.y, pos.z,  PLQ.x, PLQ.y, PLQ.z );
        f.add_mul( interpolate3DvecWrap( FFPauli,  grid.n, gpos ) , PLQ.x );
        f.add_mul( interpolate3DvecWrap( FFLondon, grid.n, gpos ) , PLQ.y );
        f.add_mul( interpolate3DvecWrap( FFelec,   grid.n, gpos ) , PLQ.z );
        //f = interpolate3DvecWrap( FFLondon,  grid.n, gpos );
        //printf( "p(%5.5e,%5.5e,%5.5e) g(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", pos.x, pos.y, pos.z, gpos.x, gpos.y, gpos.z, f.x,f.y,f.z );
    }

    inline void addForce_surf( Vec3d pos, const Vec3d& PLQ, Vec3d& f ) const {
        pos.add( shift );
        if     ( pos.z > grid.cell.c.z ){ pos.z = grid.dCell.c.z*-0.1 + grid.cell.c.z; }
        else if( pos.z < 0             ){ pos.z = grid.dCell.c.z* 0.1;                 }
        return addForce( pos, PLQ, f );
    }
    */

    inline void addForce( const Vec3d& pos, const Vec3d& PLQ, Quat4f& fe ) const {
        Vec3d gpos;
        grid.cartesian2grid(pos, gpos);
        //printf( "pos: (%g,%g,%g) PLQ: (%g,%g,%g) \n", pos.x, pos.y, pos.z,  PLQ.x, PLQ.y, PLQ.z );
        Quat4f fp=interpolate3DvecWrap( FFPauli,   grid.n, gpos );   fe.add_mul( fp, PLQ.x );
        Quat4f fl=interpolate3DvecWrap( FFLondon,  grid.n, gpos );   fe.add_mul( fl, PLQ.y );
        Quat4f fq=interpolate3DvecWrap( FFelec,    grid.n, gpos );   fe.add_mul( fq, PLQ.z );
        //Quat4f fq=interpolate3DvecWrap( FFelec,    grid.n, gpos );   fe.add_mul( fq, 1.0 );
        //if(bDebug__)printf( "CPU[0] apos(%g,%g,%g)  PLQ(%g,%g,%g)\n", pos.x,pos.y,pos.z, fp.w,fl.w,fq.w );
        //printf( "fp(%g,%g,%g|%g)*(%g) + fl((%g,%g,%g|%g)*(%g) + fq(%g,%g,%g|%g)*(%g) \n", fp.x,fp.y,fp.z,fp.e, PLQ.x,  fl.x,fl.y,fl.z,fl.e, PLQ.y,  fq.x,fq.y,fq.z,fq.e, PLQ.z );
        //printf( "E(%g,%g,%g) PLQ(%g,%g,%g)\n", fp.e,fl.e,fq.e, PLQ.x,PLQ.y,PLQ.z );

        //fe.add_mul( interpolate3DvecWrap( FFPauli,  grid.n, gpos ) , PLQ.x );
        //fe.add_mul( interpolate3DvecWrap( FFLondon, grid.n, gpos ) , PLQ.y );
        //fe.add_mul( interpolate3DvecWrap( FFelec,   grid.n, gpos ) , PLQ.z );

        //f = interpolate3DvecWrap( FFLondon,  grid.n, gpos );
        //printf( "p(%5.5e,%5.5e,%5.5e) g(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", pos.x, pos.y, pos.z, gpos.x, gpos.y, gpos.z, f.x,f.y,f.z );
    }
    inline void addForce_surf( Vec3d pos, const Vec3d& PLQ, Quat4f& f ) const {
        pos.add( shift );
        if     ( pos.z > grid.cell.c.z ){ pos.z = grid.dCell.c.z*-0.1 + grid.cell.c.z; }
        else if( pos.z < 0             ){ pos.z = grid.dCell.c.z* 0.1;                 }
        return addForce( pos, PLQ, f );
    }

    inline double eval( int n, const Vec3d* ps, const Vec3d* PLQs, Vec3d* fs, bool bSurf=false ) const {
        double E=0;
        //printf("GridFF::eval() n %i ps %li PLQs %li \n", n,  (long)ps,  (long)PLQs );
        if(bSurf){ for(int i=0; i<n; i++){ Quat4f fe=Quat4fZero; addForce_surf( ps[i], PLQs[i], fe );  fs[i].add( (Vec3d)fe.f ); E+=fe.e; } }
        else     { for(int i=0; i<n; i++){ Quat4f fe=Quat4fZero; bDebug__=(i==0); addForce     ( ps[i], PLQs[i], fe );  fs[i].add( (Vec3d)fe.f ); E+=fe.e;  
            //if(i==0)printf("CPU[0] apos(%g,%g,%g) PLQs[0](%g,%g,%g|%g) \n", n, ps[i].x,ps[i].y,ps[i].z,  PLQs[i].x,PLQs[i].y,PLQs[i].z,alpha ); 
        } }
        return E;
    }

    void init( Vec3i n, Mat3d cell, Vec3d pos0 ){
        grid.n     = n;
        grid.setCell(cell);
        grid.pos0  = pos0;
        //zmax = cell.c.z;
        allocateFFs();
    }

    void setAtoms( int natoms_, Vec3d * apos_, Vec3d * REQs_ ){
        natoms = natoms_;
        //atypes = new int  [natoms];
        apos   = apos_;
        aREQs  = REQs_;
    }

    void evalGridFFel(int natoms, Vec3d * apos, Vec3d * aREQs, Vec3d * FF ){
        //interateGrid3D( Vec3d{0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p)->void{
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = Vec3dZero;
            for(int ia=0; ia<natoms; ia++){ addAtomicForceQ( p-apos[ia], f, aREQs[ia].z ); }
            FF[ibuff]=f;
        });
    }

    void evalGridFFexp(int natoms, Vec3d * apos, Vec3d * aREQs, double alpha, double A, Vec3d * FF ){
        //interateGrid3D(  Vec3dZero, grid.n, grid.dCell, [=](int ibuff, Vec3d p){
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f =  Vec3dZero;
            for(int ia=0; ia<natoms; ia++){
                //printf( " %i (%g,%g,%g) (%g,%g)\n", ia, apos[ia].x, apos[ia].y, apos[ia].z,  aLJq[ia].x, aLJq[ia].y  );
                addAtomicForceExp( p-apos[ia], f, aREQs[ia].x, aREQs[ia].y,    alpha );
                //addAtomicForceExp( p-apos[ia], f, aLJq[ia].x, aLJq[ia].y,    alpha*2 );
                //addAtomicForceExp( p-apos[ia], f, aLJq[ia].x, aLJq[ia].y*-2, alpha   );
            }
            //printf( " >> %i %i (%g,%g,%g) %g \n", ibuff, natoms, f.x, f.y, f.z, A  );
            FF[ibuff]=f*A;
            //printf( " %i (%g,%g,%g) \n", ibuff, p.x, p.y, p.z );
            //FF[ibuff]=p;
        });
    }

    void evalGridFFs(int natoms, Vec3d * apos, Vec3d * REQs ){
        double R2damp=Rdamp*Rdamp;
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int ia=0; ia<natoms; ia++){
                Vec3d dp; dp.set_sub( p, apos[ia] );
                Vec3d REQi = aREQs[ia];
                double r      = dp.norm();
                double ir     = 1/(r+R2damp);
                double expar  = exp( alpha*(r-REQi.x) );
                double fexp   = alpha*expar*REQi.y*ir*2;
                double eCoul  = COULOMB_CONST*REQi.z*ir;
                qp.e+=expar*expar; qp.f.add_mul( dp, fexp*expar  ); // repulsive part of Morse
                ql.e+=expar*2;     ql.f.add_mul( dp, fexp        ); // attractive part of Morse
                qe.e+=eCoul;       qe.f.add_mul( dp, eCoul*ir*ir ); // Coulomb
            }
            if(FFPauli)  FFPauli [ibuff]=(Quat4f)qp;
            if(FFLondon) FFLondon[ibuff]=(Quat4f)ql;
            if(FFelec)   FFelec  [ibuff]=(Quat4f)qe;
        });
    }

    void evalGridFFs(int natoms_, Vec3d * apos_, Vec3d * REQs_, Vec3i nPBC ){
        printf( "GridFF::evalGridFFs() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, grid.pos0.x,grid.pos0.y,grid.pos0.z );
        //printf( "GridFF nPBC(%i,%i,%i) K %g R %g R2Q %g \n", nPBC.x,nPBC.y,nPBC.z, alpha, Rdamp, Rdamp*Rdamp );
        //for(int i=0; i<natoms_; i++){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", i, apos_[i].x, apos_[i].y, apos_[i].z, REQs_[i].z ); }
        //int Q=0; for(int i=0; i<natoms_; i++){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", i, apos_[i].x, apos_[i].y, apos_[i].z, REQs_[i].z ); Q+=REQs_[i].z; }; printf("evalGridFFs Qtot=%g \n", Q );
        //Rdamp = 0.1; // WARRNIN DEBUG;
        //Rdamp = 1.0; // WARRNIN DEBUG;
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            double R2damp=Rdamp*Rdamp;    
            double K=alpha;
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            //if(ibuff<(grid.n.x*grid.n.y))printf( "evalGridFFs p(%g,%g,%g)\n", p.x,p.y,p.z );
            for(int ia=0; ia<natoms_; ia++){
                Vec3d dp0; dp0.set_sub( p, apos_[ia] );
                Vec3d REQi = REQs_[ia];
                if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    //Vec3d  dp = dp0;
                    double r2     = dp.norm2();
                    double r      = sqrt(r2);
                    // ----- Morse
                    double e      = exp( K*(r-REQi.x) );
                    double eM     = e*REQi.y;
                    double de     = K*eM*-2/r;
                    // ---- Coulomb
                    double ir2    = 1/(r2+R2damp);
                    double ir     = sqrt(ir2);
                    double eQ     = COULOMB_CONST*REQi.z*ir;
                    
                    // --- store
                    qp.e+=eM*e; qp.f.add_mul( dp, de*e   ); // repulsive part of Morse
                    ql.e+=eM*2; ql.f.add_mul( dp, de     ); // attractive part of Morse
                    qe.e+=eQ;   qe.f.add_mul( dp, eQ*ir2 ); // Coulomb

                }}}
            }
            if(FFPauli)  FFPauli [ibuff]=(Quat4f)qp;
            if(FFLondon) FFLondon[ibuff]=(Quat4f)ql;
            if(FFelec)   FFelec  [ibuff]=(Quat4f)qe;
        });
    }

    void evalGridR(int natoms, Vec3d * apos, Vec3d * REQs, Vec3i nPBC ){
        printf( "GridFF::evalGridR() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, grid.pos0.x,grid.pos0.y,grid.pos0.z );
        double R2damp=Rdamp*Rdamp;
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int iat=0; iat<natoms; iat++){
                Vec3d dp0; dp0.set_sub( p, apos[iat] );
                Vec3d REQi = aREQs[iat];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    double r      = dp.norm();
                    qp.e+=fmax(0, REQi.x-r);
                    ql.e+=0;    
                    qe.e+=0;      

                }}}
            }
            if(FFPauli)  FFPauli [ibuff]=(Quat4f)qp;
            if(FFLondon) FFLondon[ibuff]=(Quat4f)ql;
            if(FFelec)   FFelec  [ibuff]=(Quat4f)qe;
        });
    }

    void evalGridFFs( Vec3i nPBC){
        //evalGridFFexp( natoms, apos, aREQs, alpha*2,  1, FFPauli  );
        //evalGridFFexp( natoms, apos, aREQs, alpha  , 1, FFLondon );  // -2.0 coef is in  REQ2PLQ
        //evalGridFFel ( natoms, apos, aREQs,              FFelec   );
        //evalGridFFs( natoms, apos, aREQs );
        if(iDebugEvalR>0){ evalGridR  ( natoms, apos, aREQs, nPBC );}
        else             { evalGridFFs( natoms, apos, aREQs, nPBC ); }
    }

    void evalCombindGridFF( Vec3d REQ, Quat4f * FF ){
        Vec3d PLQ = REQ2PLQ( REQ, alpha );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            //Vec3d f =  Vec3dZero;
            Quat4f f = Quat4fZero;
            if(FFPauli ) f.add_mul( FFPauli [ibuff], PLQ.x );
            if(FFLondon) f.add_mul( FFLondon[ibuff], PLQ.y );
            if(FFelec  ) f.add_mul( FFelec  [ibuff], PLQ.z );
            FF[ibuff] =  f;
        });
    }

    void evalCombindGridFF_CheckInterp( Vec3d REQ, Quat4f * FF ){
        Vec3d PLQ = REQ2PLQ( REQ, alpha );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4f f = Quat4fZero;
            addForce( p+Vec3d{0.1,0.1,0.1}, PLQ, f );
            FF[ibuff] = f;
        });
    }

    /*
    void evalFFline( int n, Vec3d p0, Vec3d p1, Vec3d PLQ, Vec3d * pos, Vec3d * fs ){
        Vec3d dp = p1-p0; dp.mul(1.0/(n-1));
        Vec3d  p = p0;
        for(int i=0; i<n; i++){
            if(fs ){
                Vec3d fp =  Vec3dZero;
                Vec3d fl =  Vec3dZero;
                Vec3d fq =  Vec3dZero;
                for(int ia=0; ia<natoms; ia++){
                    Vec3d d = p-apos[ia];
                    addAtomicForceExp( d, fp, aREQs[ia].x, aREQs[ia].y, 2*alpha );
                    addAtomicForceExp( d, fl, aREQs[ia].x, aREQs[ia].y,   alpha );
                    //addAtomicForceQ  ( d, fq, aLJq[ia].z );
                    //printf( "%i %i %g  (%g,%g)   %g \n", i, ia, d.z, aLJq[ia].x, aLJq[ia].y, alpha  );
                }
                fs[i] = fp*PLQ.x + fl*PLQ.y + fq*PLQ.z;
            }
            if(pos)pos[i]=p;
            //printf("%i %20.10f %20.10f %20.10f    %20.10e %20.10e %20.10e\n" , i, pos[i].x, pos[i].y, pos[i].z, fs[i].x, fs[i].y, fs[i].z );
            p.add(dp);
        }
    }

    void evalFFlineToFile( int n, Vec3d p0, Vec3d p1, Vec3d REQ, const char * fname ){
        Vec3d * pos = new Vec3d[n];
        Vec3d * fs  = new Vec3d[n];
        evalFFline( n, p0, p1, REQ2PLQ(REQ, alpha), pos, fs );
        FILE* pfile;
        pfile = fopen(fname, "w" );
        for(int i=0; i<n; i++){
            fprintf( pfile, "%i %20.10f %20.10f %20.10f    %20.10e %20.10e %20.10e\n", i, pos[i].x, pos[i].y, pos[i].z, fs[i].x, fs[i].y, fs[i].z);
        }
        fclose(pfile);
        delete [] pos; delete [] fs;
    }
    */

/*
void evalGridFFs_symetrized( Vec3i nPBC, double cmax=0.1 ){
    printf( "DEBUG evalGridFFs_symetrized() \n" );
    double cmin=1-cmax;
    std::vector<Vec3d> apos_  ;(%g,%g)
    std::vector<Vec3d> REQs_  ;
    std::vector<int>   atypes_;
    Mat3d M; grid.cell.invert_T_to( M );
    const Vec3d& a = grid.cell.a;
    const Vec3d& b = grid.cell.b;
    for(int i=0; i<natoms; i++){
        Vec3d p_;
        int typ        = atypes[i];
        Vec3d        Q = aREQs[i];
        const Vec3d& p = apos[i];
        M.dot_to( p,p_);
        bool alo  = p_.a < cmax;
        bool ahi  = p_.a > cmin;
        bool blo  = p_.b < cmax;
        bool bhi  = p_.b > cmin;
        bool aa   = (alo||ahi);
        bool bb   = (blo||ahi);
        double w = 1./( (1+aa) * (1+bb) ); // number of replicas ?
        Q.z*=w; // Q
        Q.y*=w; // E0
        apos_.push_back(p);  REQs_.push_back(Q);   atypes_.push_back(typ); 
        if(aa){
            Vec3d p_=p;
            if     ( alo ){ p_.add(a); }
            else          { p_.sub(a); };
            apos_.push_back(p_); REQs_.push_back(Q);   atypes_.push_back(typ); 
        }
        if(bb){
            Vec3d p_=p;
            if     ( alo ){ p_.add(b); }
            else          { p_.sub(b); };
            apos_.push_back(p_); REQs_.push_back(Q);   atypes_.push_back(typ); 
            if(aa){
                if     ( alo ){ p_.add(a); }
                else          { p_.sub(a); };
                apos_.push_back(p_); REQs_.push_back(Q);  atypes_.push_back(typ); 
            }
        }
    }
    
    printf( "na %i | %i %i %i \n", natoms, apos_.size(), REQs_.size(), atypes_.size() );
    params_glob->saveXYZ( "symtrized.xyz", apos_.size() , &atypes_[0] , &apos_[0], "#", &REQs_[0] );

    evalGridFFs( apos_.size(), &apos_[0], &REQs_[0], nPBC );

}
*/

void evalGridFFs_symetrized( Vec3i nPBC, double d=0.1 ){
    printf( "DEBUG evalGridFFs_symetrized() \n" );
    double cmax =-0.5+d;
    double cmin = 0.5-d;
    std::vector<Vec3d> apos_  ;
    std::vector<Vec3d> REQs_  ;
    std::vector<int>   atypes_;
    Mat3d M; grid.cell.invert_T_to( M );
    const Vec3d& a = grid.cell.a;
    const Vec3d& b = grid.cell.b;
    printf( "DEBUG evalGridFFs_symetrized() cmin %g cmax %g \n", cmin, cmax );
    for(int i=0; i<natoms; i++){
        Vec3d p_;
        int typ        = atypes[i];
        Vec3d        Q = aREQs[i];
        const Vec3d& p = apos[i];
        M.dot_to( p,p_);
        bool alo  = p_.a < cmax;
        bool ahi  = p_.a > cmin;
        bool blo  = p_.b < cmax;
        bool bhi  = p_.b > cmin;
        bool aa   = (alo||ahi);
        bool bb   = (blo||ahi);

        printf( "atom[%i](%g,%g) (%g,%g)\n", i, p.y,p.y, p_.a, p_.b );
        double w = 1./( (1+aa) * (1+bb) ); // number of replicas ?
        Q.z*=w; // Q
        Q.y*=w; // E0
        apos_.push_back(p);  REQs_.push_back(Q);   atypes_.push_back(typ); 
        if(aa){
            Vec3d p_=p;
            if     ( alo ){ p_.add(a); }
            else          { p_.sub(a); };
            apos_.push_back(p_); REQs_.push_back(Q);   atypes_.push_back(typ); 
        }
        if(bb){
            Vec3d p_=p;
            if     ( alo ){ p_.add(b); }
            else          { p_.sub(b); };
            apos_.push_back(p_); REQs_.push_back(Q);   atypes_.push_back(typ); 
            if(aa){
                if     ( alo ){ p_.add(a); }
                else          { p_.sub(a); };
                apos_.push_back(p_); REQs_.push_back(Q);  atypes_.push_back(typ); 
            }
        }
    }
    //printf( "na %i | %i %i %i \n", natoms, apos_.size(), REQs_.size(), atypes_.size() );
    //params_glob->saveXYZ( "symtrized.xyz", apos_.size() , &atypes_[0] , &apos_[0], "#", &REQs_[0] );
    evalGridFFs( apos_.size(), &apos_[0], &REQs_[0], nPBC );

}

 #ifdef IO_utils_h
    bool tryLoad( const char* fname_Coulomb, const char* fname_Pauli, const char* fname_London, bool recalcFF=false, Vec3i nPBC={1,1,0}, bool bSaveDebugXSFs=false, bool bSymetrized=false ){
        //printf( "DEBUG GridFF::tryLoad() 0 \n" );
        //printf( "DEBUG GridFF::tryLoad() fname_Pauli >>%s<< fname_London >>%s<< fname_Coulomb >>%s<< \n", fname_Pauli, fname_London, fname_Coulomb );
        //printDir( "../" );
        //printDir( "./" );
        { FILE* f=fopen( fname_Pauli,  "rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Pauli );  recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_London, "rb"); if(0==f){ printf("File(%s) Not Found\n", fname_London);  recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_Coulomb,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Coulomb); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        printf( "DEBUG GridFF::tryLoad() recalcFF %i \n", recalcFF );
        //printf( "fname_Pauli(%s) fname_London(%s) fname_Coulomb(%s) \n", fname_Pauli, fname_London, fname_Coulomb );
        //int nbyte= grid.getNtot()*sizeof(Vec3d);
        int nbyte= grid.getNtot()*sizeof(Quat4f);
        if( recalcFF ){
            printf( "\nBuilding GridFF for substrate ... (please wait... )\n" );
            printf("DEBUG tryLoad() bSymetrized %i \n", bSymetrized );
            if(bSymetrized){
                evalGridFFs_symetrized( nPBC, 0.1 );
                //evalGridFFs( natoms, apos, aREQs, nPBC );
            }else{
                evalGridFFs( nPBC );
            }
            if(bSaveDebugXSFs){
                if(FFPauli)  grid.saveXSF( "FFLond_E.xsf", (float*)FFLondon, 4,3  );
                if(FFLondon) grid.saveXSF( "FFelec_E.xsf", (float*)FFelec,   4,3  );
                if(FFelec )  grid.saveXSF( "FFPaul_E.xsf", (float*)FFPauli,  4,3  );
                if(FFPauli)  grid.saveXSF( "FFLond_z.xsf", (float*)FFLondon, 4,2  );
                if(FFLondon) grid.saveXSF( "FFelec_z.xsf", (float*)FFelec,   4,2  );
                if(FFelec )  grid.saveXSF( "FFPaul_z.xsf", (float*)FFPauli,  4,2  );
            }
            if(FFPauli)  saveBin( fname_Pauli,    nbyte, (char*)FFPauli  );
            if(FFLondon) saveBin( fname_London,   nbyte, (char*)FFLondon );
            if(FFelec )  saveBin( fname_Coulomb,  nbyte, (char*)FFelec   );
        }else{
            if(FFPauli)  loadBin( fname_Pauli,    nbyte, (char*)FFPauli  );
            if(FFLondon) loadBin( fname_London,   nbyte, (char*)FFLondon );
            if(FFelec )  loadBin( fname_Coulomb,  nbyte, (char*)FFelec   );
        }
        return recalcFF;
    }
#endif

}; // RigidSubstrate


#endif
