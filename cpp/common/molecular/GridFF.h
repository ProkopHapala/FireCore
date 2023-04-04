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

inline double evalDipole( int n, Vec3d* ps, Quat4d* REQs, Vec3d& Dout, Vec3d& p0out ){
    // we want to approximate charge distribution by charge `Q` and dipole `D` around center `C`
    // we choose center `c` so that it minimizes size of dipole `|d|^2`
    // D = sum_i{ ( p_i - C ) * q_i }
    // |D|^2       = Dx^2 + Dy^2 + Dz^2
    // minimize   |D|^2 with respocet of position of center C
    // d|D|^2/dCx  = 2*Dx * d{sum_i{(x_i-Cx)*qi}}/dCx   = 0
    //             = -2*Dx *   sum_i{ qi } = 0
    //             = -2*Dx * Q = 0     ( if Q=0 it does not matter)
    //   if Q is not =0
    //  sum_i{ ( p_i - C ) * q_i } = 0
    //  sum_i{ p_i*q_i } = Q*C
    //  C = sum_i{ p_i * q_i } / Q
    Vec3d  D  = Vec3dZero;
    Vec3d  p0 = Vec3dZero;
    double Q  = 0;
    double W  = 0; 
    double Qsafe   = 1e-4;
    double Q2safe  = Qsafe*Qsafe;
    double Wscale  = 1e-4; 
    for(int i=0; i<n; i++){
        const Vec3d& p  = ps  [i];
        double       qi = REQs[i].z;
        //double wi = (Q2safe + qi*qi)*Wscale;  // Q2safe is just for case when all charges are zero
        double wi = Q2safe;
        Q += qi;
        W += wi;
        D .add_mul( p, qi );  // dipome
        p0.add_mul( p, wi );  // this p0 is used only if Q~=0
        //printf( "evalDipole[%i] q %g p(%g,%g,%g)\n", i, qi, p.x,p.y,p.z );
    }
    double denom = 1./( W*W + Q*Q );     // if 
    p0.mul    (      W*denom );   // if Q~=0 this term dominates
    //printf( "denom %g Q %g W %g \n", denom, Q, W );
    //printf( "p0W(%g,%g,%g)\n", p0.x,p0.y,p0.z );
    p0.add_mul(  D,  Q*denom );   // if Q!=0 this term dominates
    D .add_mul( p0, -Q ); 
    //printf( "p0Q(%g,%g,%g)\n", D.x*Q*denom, D.y*Q*denom, D.z*Q*denom );
    //printf( "p0 (%g,%g,%g)\n", p0.x,p0.y,p0.z );
    Dout  = D;
    p0out = p0;
    return Q;
}

class GridFF{ public: 
    GridShape   grid;

    Quat4f *FFPaul = 0;
    Quat4f *FFLond = 0;
    Quat4f *FFelec = 0;

    //Vec3d  *FFtot    = NULL; // total FF is not used since each atom-type has different linear combination

    // ------ ToDo: this should be put inside NBFF
    int  natoms     = 0;
    int    * atypes = 0;
    Vec3d  * apos   = 0;   // atomic position
    Quat4d * REQs  = 0;
    //Quat4d * PLQ = 0;

    std::vector<Vec3d>  apos_  ;
    std::vector<Quat4d> REQs_  ;
    std::vector<int>    atypes_;

    Vec3i nPBC{1,1,0};
    Vec3d shift = Vec3dZero;

    // dipole approximation
    Vec3d dip_p0;
    Vec3d dip;
    double Q;

    double alpha  =  1.5;
    double Rdamp  =  1.0;

    int iDebugEvalR = 0;
    bool bCellSet = false;

    double findTop(){ double zmax=-1e+300; for(int i=0;i<natoms; i++){ double z=apos[i].z; if(z>zmax)zmax=z; }; return zmax; }

    void bindSystem(int natoms_, int* atypes_, Vec3d* apos_, Quat4d* REQs_ ){
        natoms=natoms_; atypes=atypes_; apos=apos_; REQs=REQs_;
    }

    void allocateFFs(){
        int ntot = grid.getNtot();
        _realloc( FFPaul , ntot );
        _realloc( FFLond, ntot );
        _realloc( FFelec  , ntot );
    }
    
    void allocateAtoms(int natoms_){
        natoms = natoms_;
        //_realloc( atypes = new int[natoms];
        _realloc(apos ,natoms);
        _realloc(REQs,natoms);
    }

    int loadCell(const char * fname ){ return grid.loadCell(fname ); }

    void evalCellDipole(){
        Q = evalDipole( natoms, apos, REQs, dip, dip_p0 );
        printf( "GridFF::evalDipole(na=%i): D(%g,%g,%g|%g) p0(%g,%g,%g) \n", natoms, dip.x,dip.y,dip.z,Q,   dip_p0.x,dip_p0.y,dip_p0.z );
    }

inline Quat4f getForce( const Vec3d& p, const Quat4f& PLQ, bool bSurf=true, double off=10. ) const {
    Vec3d u;
    grid.iCell.dot_to( p - grid.pos0, u );
    Vec3i n = grid.n;
    //printf( "pos(%g,%g,%g) u(%g,%g,%g)  p0(%g,%g,%g) \n", p.x,p.y,p.z, u.x,u.y,u.z,    grid.pos0.x,grid.pos0.y,grid.pos0.z );
    if( bSurf ){ double ivnz=1./n.z; double uzmax=1-ivnz*1.5; double uzmin=ivnz+0.5; if(u.z<uzmin){ u.z=uzmin; }else if(u.z>uzmax){u.z=uzmax; } }
    u.add(off);
    u.x=(u.x-(int)u.x)*n.x;
    u.y=(u.y-(int)u.y)*n.y;
    u.z=(u.z-(int)u.z)*n.z;

    //printf( " u(%g,%g,%g) \n", u.x,u.y,u.z );

    // -----
	const int   ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const float tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    const float mx = 1-tx      ,  my = 1-ty      ,  mz = 1-tz      ;
    //------
    int jx = ix+1; jx=(jx<n.x)?jx:0;
    int jy = iy+1; jy=(jy<n.y)?jy:0;
    int jz = iz+1; jz=(jz<n.z)?jz:0;
	//------
	const float f00 = my*mx; 
    const float f01 = my*tx; 
    const float f10 = ty*mx; 
    const float f11 = ty*tx;
    const int   i00 = n.x*(iy + n.y*iz);
    const int   i10 = n.x*(jy + n.y*iz);
    const int   i11 = n.x*(jy + n.y*jz);
    const int   i01 = n.x*(iy + n.y*jz);
    const int 
        i000=ix+i00, i100=jx+i00,
        i010=ix+i10, i110=jx+i10,
        i011=ix+i11, i111=jx+i11,
        i001=ix+i01, i101=jx+i01;
    const float 
        f000=mx*f00, f100=tx*f00,
        f010=mx*f10, f110=tx*f10,
        f011=mx*f11, f111=tx*f11,
        f001=mx*f01, f101=tx*f01;

    { // DEBUG
        Quat4f fDBG = FFPaul[ i000 ];
        /*
         ((FFPaul[ i000 ]*f000) + (FFPaul[ i100 ]*f100)
        + (FFPaul[ i010 ]*f010) + (FFPaul[ i110 ]*f110)  
        + (FFPaul[ i011 ]*f011) + (FFPaul[ i111 ]*f111)
        + (FFPaul[ i001 ]*f001) + (FFPaul[ i101 ]*f101));
        */
        printf( " u(%g,%g,%g) fe(%g,%g,%g|%g) p0(%g,%g,%g) \n", u.x,u.y,u.z, fDBG.x,fDBG.y,fDBG.z,fDBG.w, grid.pos0.x,grid.pos0.y,grid.pos0.z );
    }
	return  // 3 * 8 * 4 = 96 floats   // SIMD optimize ?????
         ((FFPaul[ i000 ]*f000) + (FFPaul[ i100 ]*f100)
        + (FFPaul[ i010 ]*f010) + (FFPaul[ i110 ]*f110)  
        + (FFPaul[ i011 ]*f011) + (FFPaul[ i111 ]*f111)
        + (FFPaul[ i001 ]*f001) + (FFPaul[ i101 ]*f101))*PLQ.x

        +((FFLond[ i000 ]*f000) + (FFLond[ i100 ]*f100)
        + (FFLond[ i010 ]*f010) + (FFLond[ i110 ]*f110)  
        + (FFLond[ i011 ]*f011) + (FFLond[ i111 ]*f111)
        + (FFLond[ i001 ]*f001) + (FFLond[ i101 ]*f101))*PLQ.y

        +((FFelec[ i000 ]*f000) + (FFelec[ i100 ]*f100)
        + (FFelec[ i010 ]*f010) + (FFelec[ i110 ]*f110)  
        + (FFelec[ i011 ]*f011) + (FFelec[ i101 ]*f111)
        + (FFelec[ i001 ]*f001) + (FFelec[ i101 ]*f101))*PLQ.z;
}


    inline void addForce( const Vec3d& pos, const Quat4f& PLQ, Quat4f& fe ) const {
        //printf( "GridFF::addForce() \n" );
        //printf( "GridFF::addForce() pointers(%li,%li,%li)\n", FFPaul, FFLond, FFelec );
        Vec3f gpos; grid.cartesian2grid(pos, gpos);
        //printf( "pos: (%g,%g,%g) PLQ: (%g,%g,%g) pointers(%li,%li,%li)\n", pos.x, pos.y, pos.z,  PLQ.x, PLQ.y, PLQ.z, FFPaul, FFLond, FFelec );
        Quat4f fp=interpolate3DvecWrap( FFPaul,  grid.n, gpos );   fe.add_mul( fp, PLQ.x );
        Quat4f fl=interpolate3DvecWrap( FFLond,  grid.n, gpos );   fe.add_mul( fl, PLQ.y );
        Quat4f fq=interpolate3DvecWrap( FFelec,  grid.n, gpos );   fe.add_mul( fq, PLQ.z );
        //Quat4f fq=interpolate3DvecWrap( FFelec,    grid.n, gpos );   fe.add_mul( fq, 1.0 );
        //if(bDebug__)printf( "CPU[0] apos(%g,%g,%g)  PLQ(%g,%g,%g)\n", pos.x,pos.y,pos.z, fp.w,fl.w,fq.w );
        //printf( "fp(%g,%g,%g|%g)*(%g) + fl((%g,%g,%g|%g)*(%g) + fq(%g,%g,%g|%g)*(%g) \n", fp.x,fp.y,fp.z,fp.e, PLQ.x,  fl.x,fl.y,fl.z,fl.e, PLQ.y,  fq.x,fq.y,fq.z,fq.e, PLQ.z );
        //printf( "E(%g,%g,%g) PLQ(%g,%g,%g)\n", fp.e,fl.e,fq.e, PLQ.x,PLQ.y,PLQ.z );
        //fe.add_mul( interpolate3DvecWrap( FFPaul,  grid.n, gpos ) , PLQ.x );
        //fe.add_mul( interpolate3DvecWrap( FFLond, grid.n, gpos ) , PLQ.y );
        //fe.add_mul( interpolate3DvecWrap( FFelec,   grid.n, gpos ) , PLQ.z );
        //f = interpolate3DvecWrap( FFLond,  grid.n, gpos );
        //printf( "p(%5.5e,%5.5e,%5.5e) g(%5.5e,%5.5e,%5.5e) f(%5.5e,%5.5e,%5.5e) \n", pos.x, pos.y, pos.z, gpos.x, gpos.y, gpos.z, f.x,f.y,f.z );
    }
    inline void addForce_surf( Vec3d pos, const Quat4f PLQ, Quat4f& f ) const {
        pos.add( shift );
        if     ( pos.z > grid.cell.c.z ){ pos.z = grid.dCell.c.z*-0.1 + grid.cell.c.z; }
        else if( pos.z < 0             ){ pos.z = grid.dCell.c.z* 0.1;                 }
        return addForce( pos, PLQ, f );
    }

    inline double eval( int n, const Vec3d* ps, const Quat4f* PLQs, Vec3d* fs, bool bSurf=false ) const {
        double E=0;
        //printf("GridFF::eval() n %i ps %li PLQs %li \n", n,  (long)ps,  (long)PLQs );
        if(bSurf){ for(int i=0; i<n; i++){ Quat4f fe=Quat4fZero;                  addForce_surf( ps[i], (Quat4f)PLQs[i], fe );  fs[i].add( (Vec3d)fe.f ); E+=fe.e; } }
        else     { for(int i=0; i<n; i++){ Quat4f fe=Quat4fZero; bDebug__=(i==0); addForce     ( ps[i], (Quat4f)PLQs[i], fe );  fs[i].add( (Vec3d)fe.f ); E+=fe.e;  
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

    void setAtoms( int natoms_, Vec3d * apos_, Quat4d * REQs_ ){
        natoms = natoms_;
        //atypes = new int  [natoms];
        apos   = apos_;
        REQs  = REQs_;
    }

    void evalGridFFel(int natoms, Vec3d * apos, Quat4d * REQs, Vec3d * FF ){
        //interateGrid3D( Vec3d{0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p)->void{
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = Vec3dZero;
            for(int ia=0; ia<natoms; ia++){ addAtomicForceQ( p-apos[ia], f, REQs[ia].z ); }
            FF[ibuff]=f;
        });
    }

    void evalGridFFexp(int natoms, Vec3d * apos, Quat4d * REQs, double alpha, double A, Vec3d * FF ){
        //interateGrid3D(  Vec3dZero, grid.n, grid.dCell, [=](int ibuff, Vec3d p){
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f =  Vec3dZero;
            for(int ia=0; ia<natoms; ia++){
                //printf( " %i (%g,%g,%g) (%g,%g)\n", ia, apos[ia].x, apos[ia].y, apos[ia].z,  aLJq[ia].x, aLJq[ia].y  );
                addAtomicForceExp( p-apos[ia], f, REQs[ia].x, REQs[ia].y,    alpha );
                //addAtomicForceExp( p-apos[ia], f, aLJq[ia].x, aLJq[ia].y,    alpha*2 );
                //addAtomicForceExp( p-apos[ia], f, aLJq[ia].x, aLJq[ia].y*-2, alpha   );
            }
            //printf( " >> %i %i (%g,%g,%g) %g \n", ibuff, natoms, f.x, f.y, f.z, A  );
            FF[ibuff]=f*A;
            //printf( " %i (%g,%g,%g) \n", ibuff, p.x, p.y, p.z );
            //FF[ibuff]=p;
        });
    }

    void evalGridFFs(int natoms_, Vec3d * apos_, Quat4d * REQs_){
        printf( "GridFF::evalGridFFs() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, grid.pos0.x,grid.pos0.y,grid.pos0.z );
        //printf( "GridFF nPBC(%i,%i,%i) K %g R %g R2Q %g \n", nPBC.x,nPBC.y,nPBC.z, alpha, Rdamp, Rdamp*Rdamp );
        //for(int i=0; i<natoms_; i++){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", i, apos_[i].x, apos_[i].y, apos_[i].z, REQs_[i].z ); }
        //int Q=0; for(int i=0; i<natoms_; i++){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", i, apos_[i].x, apos_[i].y, apos_[i].z, REQs_[i].z ); Q+=REQs_[i].z; }; printf("evalGridFFs Qtot=%g \n", Q );
        //Rdamp = 0.1; // WARRNIN DEBUG;
        //Rdamp = 1.0; // WARRNIN DEBUG;
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            const double R2damp=Rdamp*Rdamp;    
            const double K=-alpha;
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            //if(ibuff<(grid.n.x*grid.n.y))printf( "evalGridFFs p(%g,%g,%g)\n", p.x,p.y,p.z );
            for(int ia=0; ia<natoms_; ia++){
                Vec3d dp0; dp0.set_sub( p, apos_[ia] );
                Quat4d REQi = REQs_[ia];
                //if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    //Vec3d  dp = dp0;
                    double r2     = dp.norm2();
                    double r      = sqrt(r2);
                    // ---- Coulomb
                    double ir2    = 1/(r2+R2damp);
                    double eQ     = COULOMB_CONST*REQi.z*sqrt(ir2);
                    // ----- Morse
                    double e      = exp( K*(r-REQi.x) );
                    double eM     = e*REQi.y;
                    double de     = 2*K*eM/r;                    
                    // --- store
                    qp.e+=eM*e;  qp.f.add_mul( dp,  de*e  ); // repulsive part of Morse
                    ql.e+=eM*-2; ql.f.add_mul( dp, -de    ); // attractive part of Morse
                    qe.e+=eQ;    qe.f.add_mul( dp, eQ*ir2 ); // Coulomb

                }}}
            }
            if(FFPaul)  FFPaul [ibuff]=(Quat4f)qp;
            if(FFLond) FFLond[ibuff]=(Quat4f)ql;
            if(FFelec)   FFelec  [ibuff]=(Quat4f)qe;
        });
    }

    void makeGridFF_omp(int natoms_, Vec3d * apos_, Quat4d * REQs_){
        printf( "GridFF::makeGridFF_omp() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z,  grid.pos0.x,grid.pos0.y,grid.pos0.z );
        const double R2damp=Rdamp*Rdamp;    
        const double K=-alpha;
        // ----- makePBC shifts
        int nshift = (nPBC.a*2+1)*(nPBC.b*2+1)*(nPBC.c*2+1);
        std::vector<Vec3d> shifts{nshift};
        int ipbc=0;
        for(int ia=-nPBC.a; ia<=nPBC.a; ia++){ 
            for(int ib=-nPBC.b; ib<=nPBC.b; ib++){ 
                for(int ic=-nPBC.c; ic<=nPBC.c; ic++){
                    shifts[ipbc] = grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    ipbc++;
                }
            }
        }
        // ----- Project Forcefield on the grid
        #pragma omp parallel
        {
        #pragma omp for simd
        for ( int iz=0; iz<grid.n.z; iz++ ){
            for ( int iy=0; iy<grid.n.y; iy++ ){
                for ( int ix=0; ix<grid.n.x; ix++ ){
                    const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
                    Quat4d qp = Quat4dZero;
                    Quat4d ql = Quat4dZero;
                    Quat4d qe = Quat4dZero;
                    for(int ia=0; ia<natoms_; ia++){
                        const Vec3d dp0   = pos - apos_[ia];
                        const Quat4d REQi = REQs_[ia];
                        //if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }              
                        for(int ipbc=0; ipbc<nshift; ipbc++ ){
                            const Vec3d  dp = dp0 + shifts[ipbc];
                            //Vec3d  dp = dp0;
                            double r2     = dp.norm2();
                            double r      = sqrt(r2);
                            // ---- Coulomb
                            double ir2    = 1/(r2+R2damp);
                            double eQ     = COULOMB_CONST*REQi.z*sqrt(ir2);
                            // ----- Morse
                            double e      = exp( K*(r-REQi.x) );
                            double eM     = e*REQi.y;
                            double de     = 2*K*eM/r;                    
                            // --- store
                            qp.e+=eM*e;   qp.f.add_mul( dp,  de*e  ); // repulsive part of Morse
                            ql.e+=eM*-2.; ql.f.add_mul( dp, -de    ); // attractive part of Morse
                            qe.e+=eQ;     qe.f.add_mul( dp, eQ*ir2 ); // Coulomb
                        }
                    }
                    const int   ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
                    FFPaul [ibuff]=(Quat4f)qp;
                    FFLond[ibuff]=(Quat4f)ql;
                    FFelec  [ibuff]=(Quat4f)qe;
                }
            }
        }
        }
    }

    void evalGridR(int natoms, Vec3d * apos, Quat4d * REQs ){
        printf( "GridFF::evalGridR() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, grid.pos0.x,grid.pos0.y,grid.pos0.z );
        double R2damp=Rdamp*Rdamp;
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int iat=0; iat<natoms; iat++){
                Vec3d dp0; dp0.set_sub( p, apos[iat] );
                Quat4d REQi = REQs[iat];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + grid.cell.a*ia + grid.cell.b*ib + grid.cell.c*ic;
                    double r      = dp.norm();
                    qp.e+=fmax(0, REQi.x-r);
                    ql.e+=0;    
                    qe.e+=0;      
                }}}
            }
            if(FFPaul)  FFPaul [ibuff]=(Quat4f)qp;
            if(FFLond) FFLond[ibuff]=(Quat4f)ql;
            if(FFelec)   FFelec  [ibuff]=(Quat4f)qe;
        });
    }

    void evalGridFFs( Vec3i nPBC_=Vec3i{-1,-1,-1} ){
        if(nPBC_.x>=0) nPBC=nPBC_;
        //evalGridFFexp( natoms, apos, REQs, alpha*2,  1, FFPaul  );
        //evalGridFFexp( natoms, apos, REQs, alpha  , 1, FFLond );  // -2.0 coef is in  REQ2PLQ
        //evalGridFFel ( natoms, apos, REQs,              FFelec   );
        //evalGridFFs( natoms, apos, REQs );
        if(iDebugEvalR>0){ evalGridR  ( natoms, apos, REQs );}
        else             { evalGridFFs( natoms, apos, REQs ); }
    }

    void evalCombindGridFF( Quat4d REQ, Quat4f * FF ){
        Quat4f PLQ = REQ2PLQ( REQ, alpha );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4f f = Quat4fZero;
            if(FFPaul ) f.add_mul( FFPaul [ibuff], PLQ.x );
            if(FFLond) f.add_mul( FFLond[ibuff], PLQ.y );
            if(FFelec  ) f.add_mul( FFelec  [ibuff], PLQ.z );
            FF[ibuff] =  f;
        });
    }

    void evalCombindGridFF_CheckInterp( Quat4d REQ, Quat4f * FF ){
        Quat4f PLQ = REQ2PLQ( REQ, alpha );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4f f = Quat4fZero;
            addForce( p+Vec3d{0.1,0.1,0.1}, (Quat4f)PLQ, f );
            FF[ibuff] = f;
        });
    }

void setAtomsSymetrized( int n, int* atypes, Vec3d* apos, Quat4d* REQs, double d=0.1 ){
    //printf( "evalGridFFs_symetrized() \n" );
    double cmax =-0.5+d;
    double cmin = 0.5-d;
    apos_  .clear();
    REQs_  .clear();
    atypes_.clear();
    Mat3d M; grid.cell.invert_T_to( M );
    const Vec3d& a = grid.cell.a;
    const Vec3d& b = grid.cell.b;
    //printf( "evalGridFFs_symetrized() cmin %g cmax %g \n", cmin, cmax );
    for(int i=0; i<n; i++){
        Vec3d p_;
        int        typ = atypes[i];
        Quat4d     REQ = REQs [i];
        const Vec3d& p = apos  [i];
        M.dot_to( p,p_);
        bool alo  = p_.a < cmax;
        bool ahi  = p_.a > cmin;
        bool blo  = p_.b < cmax;
        bool bhi  = p_.b > cmin;
        bool aa   = (alo||ahi);
        bool bb   = (blo||ahi);
        //printf( "atom[%i](%g,%g) (%g,%g)\n", i, p.y,p.y, p_.a, p_.b );
        double w = 1./( (1+aa) * (1+bb) ); // number of replicas ?
        REQ.z*=w; // Q
        REQ.y*=w; // E0
        apos_.push_back(p);  REQs_.push_back(REQ);   atypes_.push_back(typ); 
        if(aa){
            Vec3d p_=p;
            if     ( alo ){ p_.add(a); }
            else          { p_.sub(a); };
            apos_.push_back(p_); REQs_.push_back(REQ);   atypes_.push_back(typ); 
        }
        if(bb){
            Vec3d p_=p;
            if     ( alo ){ p_.add(b); }
            else          { p_.sub(b); };
            apos_.push_back(p_); REQs_.push_back(REQ);   atypes_.push_back(typ); 
            if(aa){
                if     ( alo ){ p_.add(a); }
                else          { p_.sub(a); };
                apos_.push_back(p_); REQs_.push_back(REQ);  atypes_.push_back(typ); 
            }
        }
    }
    printf( "setAtomsSymetrized() END na_new=%i na_old=%i \n", atypes_.size(), n );
    bindSystem( atypes_.size(), &atypes_[0], &apos_[0], &REQs_[0] );
}


void evalGridFFs_symetrized( double d=0.1, Vec3i nPBC_=Vec3i{-1,-1,-1} ){
    if(nPBC_.x>=0) nPBC=nPBC_;
    setAtomsSymetrized( natoms, atypes, apos, REQs, d );
    //printf( "na %i | %i %i %i \n", natoms, apos_.size(), REQs_.size(), atypes_.size() );
    //params_glob->saveXYZ( "symtrized.xyz", apos_.size() , &atypes_[0] , &apos_[0], "#", &REQs_[0] );
    //evalGridFFs( apos_.size(), &apos_[0], &REQs_[0] );
    makeGridFF_omp( apos_.size(), &apos_[0], &REQs_[0] );

}

 #ifdef IO_utils_h
    bool tryLoad( const char* fname_Coulomb, const char* fname_Pauli, const char* fname_London, bool recalcFF=false, Vec3i nPBC={1,1,0}, bool bSymetrized=false ){
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
                evalGridFFs_symetrized( 0.1, nPBC );
                //evalGridFFs( natoms, apos, REQs, nPBC );
            }else{
                evalGridFFs( nPBC );
            }
            if(FFPaul)  saveBin( fname_Pauli,    nbyte, (char*)FFPaul  );
            if(FFLond) saveBin( fname_London,   nbyte, (char*)FFLond );
            if(FFelec )  saveBin( fname_Coulomb,  nbyte, (char*)FFelec   );
        }else{
            if(FFPaul)  loadBin( fname_Pauli,    nbyte, (char*)FFPaul  );
            if(FFLond) loadBin( fname_London,   nbyte, (char*)FFLond );
            if(FFelec )  loadBin( fname_Coulomb,  nbyte, (char*)FFelec   );
        }
        return recalcFF;
    }
#endif

}; // RigidSubstrate


#endif
