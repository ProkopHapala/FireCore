#ifndef GridFF_h
#define GridFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Forces.h"
#include "MMFFparams.h"

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

    double findTop(){ double zmax=-1e+300; for(int i=0;i<natoms; i++){ double z=apos[i].z; if(z>zmax)zmax=z; }; printf("findTop() %i %g \n", natoms, zmax); return zmax; }

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
        fe.add_mul( interpolate3DvecWrap( FFPauli,  grid.n, gpos ) , PLQ.x );
        fe.add_mul( interpolate3DvecWrap( FFLondon, grid.n, gpos ) , PLQ.y );
        fe.add_mul( interpolate3DvecWrap( FFelec,   grid.n, gpos ) , PLQ.z );
        //f = interpolate3DvecWrap( FFLondon,  grid.n, gpos );
        //printf( "p(%5.5e,%5.5e,%5.5e) g(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", pos.x, pos.y, pos.z, gpos.x, gpos.y, gpos.z, f.x,f.y,f.z );
    }
    inline void addForce_surf( Vec3d pos, const Vec3d& PLQ, Quat4f& f ) const {
        pos.add( shift );
        if     ( pos.z > grid.cell.c.z ){ pos.z = grid.dCell.c.z*-0.1 + grid.cell.c.z; }
        else if( pos.z < 0             ){ pos.z = grid.dCell.c.z* 0.1;                 }
        return addForce( pos, PLQ, f );
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
        //interateGrid3D( (Vec3d){0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p)->void{
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = (Vec3d){0.0,0.0,0.0};
            for(int ia=0; ia<natoms; ia++){ addAtomicForceQ( p-apos[ia], f, aREQs[ia].z ); }
            FF[ibuff]=f;
        });
    }

    void evalGridFFexp(int natoms, Vec3d * apos, Vec3d * aREQs, double alpha, double A, Vec3d * FF ){
        //interateGrid3D( (Vec3d){0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p){
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = (Vec3d){0.0,0.0,0.0};
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

    void evalGridFFs(int natoms, Vec3d * apos, Vec3d * REQs, Vec3i nPBC ){
        //printf( "GridFF::evalGridFFs() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, grid.pos0.x,grid.pos0.y,grid.pos0.z );
        //printf( "GridFF nPBC(%i,%i,%i) K %g R %g R2Q %g \n", nPBC.x,nPBC.y,nPBC.z, alpha, Rdamp, Rdamp*Rdamp );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            double R2damp=Rdamp*Rdamp;    
            double K=alpha;
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int iat=0; iat<natoms; iat++){
                Vec3d dp0; dp0.set_sub( p, apos[iat] );
                Vec3d REQi = aREQs[iat];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    //Vec3d  dp = dp0;
                    double r2     = dp.norm2();
                    double r      = sqrt(r2);
                    // ----- Morse
                    double e      = exp( K*(r-REQi.x) );
                    double de     = K*e*REQi.y*-2/r;
                    double eM     = e*REQi.y;
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
            //Vec3d f = (Vec3d){0.0,0.0,0.0};
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
            addForce( p+(Vec3d){0.1,0.1,0.1}, PLQ, f );
            FF[ibuff] = f;
        });
    }

    /*
    void evalFFline( int n, Vec3d p0, Vec3d p1, Vec3d PLQ, Vec3d * pos, Vec3d * fs ){
        Vec3d dp = p1-p0; dp.mul(1.0/(n-1));
        Vec3d  p = p0;
        for(int i=0; i<n; i++){
            if(fs ){
                Vec3d fp = (Vec3d){0.0,0.0,0.0};
                Vec3d fl = (Vec3d){0.0,0.0,0.0};
                Vec3d fq = (Vec3d){0.0,0.0,0.0};
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

 #ifdef IO_utils_h
    bool tryLoad( const char* fname_Coulomb, const char* fname_Pauli, const char* fname_London, bool recalcFF=false, Vec3i nPBC={1,1,0}, bool bSaveDebugXSFs=false ){
        //printf( "DEBUG GridFF::tryLoad() 0 \n" );
        { FILE* f=fopen( fname_Pauli,  "rb"); if(0==f){ recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_London, "rb"); if(0==f){ recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_Coulomb,"rb"); if(0==f){ recalcFF=true; }else{ fclose(f); };} // Test if file exist
        //printf( "DEBUG GridFF::tryLoad() recalcFF %i \n", recalcFF );
        //int nbyte= grid.getNtot()*sizeof(Vec3d);
        int nbyte= grid.getNtot()*sizeof(Quat4f);
        if( recalcFF ){
            printf( "\nBuilding GridFF for substrate ... (please wait... )\n" );
            evalGridFFs( nPBC );
            if(bSaveDebugXSFs){
                if(FFPauli)  grid.saveXSF( "FFLond_E.xsf", (float*)FFLondon, 4,3  );
                if(FFLondon) grid.saveXSF( "FFelec_E.xsf", (float*)FFelec,   4,3  );
                if(FFelec )  grid.saveXSF( "FFPaul_E.xsf", (float*)FFPauli,  4,3  );
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
