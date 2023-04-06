#ifndef GridFF_h
#define GridFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Forces.h"
#include "MMFFparams.h"

#include "NBFF.h"

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

class GridFF : public NBFF{ public: 
    // -----  From NBFF 
    // int  natoms     = 0;
    // int    * atypes = 0;
    // Vec3d  * apos   = 0;   // atomic position
    // Quat4d * REQs  = 0;
    //Quat4d  * PLQ = 0;
    //Vec3d   * fapos  =0;
    //Quat4d  * REQs   =0;
    //Quat4i  * neighs =0;
    //Quat4i  * neighCell=0;
    //double  Rdamp  = 1.0;
    //Mat3d   lvec;
    //Vec3i   nPBC;
    //bool    bPBC=false;
    //int     npbc=0;
    //Vec3d*  shifts=0;
    //Quat4f  *PLQs   =0;  // used only in combination with GridFF

    // -----------------------
    GridShape   grid;
    Quat4f   *FFPaul = 0;
    Quat4f   *FFLond = 0;
    Quat4f   *FFelec = 0;
    //Quat4f *FFtot    = 0; // total FF is not used since each atom-type has different linear combination

    // ------ ToDo: this should be put inside NBFF
    std::vector<Vec3d>  apos_  ;
    std::vector<Quat4d> REQs_  ;
    std::vector<int>    atypes_;

    //Vec3i nPBC{1,1,0};
    Vec3d shift = Vec3dZero;

    // dipole approximation
    Vec3d dip_p0;
    Vec3d dip;
    double Q;

    //double Rdamp  =  1.0;
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

inline Quat4f getForce( Vec3d p, const Quat4f& PLQ, bool bSurf=true ) const {
    #pragma omp SIMD
    {

    Vec3d u;
    p.add(shift);
    p.sub(grid.pos0);
    grid.iCell.dot_to( p, u );
    Vec3i n = grid.n;
    //printf( "pos(%g,%g,%g) u(%g,%g,%g)  p0(%g,%g,%g) \n", p.x,p.y,p.z, u.x,u.y,u.z,    grid.pos0.x,grid.pos0.y,grid.pos0.z );
    if( bSurf ){ double ivnz=1./n.z; double uzmax=1-ivnz*1.5; double uzmin=ivnz*0.5; if(u.z<uzmin){ u.z=uzmin; }else if(u.z>uzmax){u.z=uzmax; } } // "clamp" boundrary condition in z-direction

    // ---- PBC condiction ( within outer cell bboundaries ? )
    // u.x=(u.x-(int)u.x);
    // u.y=(u.y-(int)u.y);
    // u.z=(u.z-(int)u.z);

    //printf( "GridFF::getForce() u(%g,%g,%g) p(%g,%g,%g) shift(%g,%g,%g) p0(%g,%g,%g) \n", u.x,u.y,u.z,  p.x,p.y,p.z,  shift.x,shift.y,shift.z,   grid.pos0.x,grid.pos0.y,grid.pos0.z );

    // ---- Clamp ( within outer cell bboundaries ? )
    if((u.x<0.)||(u.x>1.)||(u.y<0.)||(u.y>1.)){ return Quat4fZero; };
    //if((u.z<0.)||(u.z>1.)){ return Quat4fZero; };

    u.x*=n.x; 
    u.y*=n.y; 
    u.z*=n.z;

    //printf( " u(%g,%g,%g) \n", u.x,u.y,u.z );

    // -----
	const int   ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const float tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    const float mx = 1-tx      ,  my = 1-ty      ,  mz = 1-tz      ;
    //------
    int kx = ix+1; kx=(kx<n.x)?kx:0;
    int ky = iy+1; ky=(ky<n.y)?ky:0;
    int kz = iz+1; kz=(kz<n.z)?kz:0;
	//------
	const float f00 = my*mz; 
    const float f10 = ty*mz; 
    const float f01 = my*tz; 
    const float f11 = ty*tz;
    const int   i00 = n.x*(iy + n.y*iz);
    const int   i01 = n.x*(iy + n.y*kz);
    const int   i10 = n.x*(ky + n.y*iz);
    const int   i11 = n.x*(ky + n.y*kz);
    const int 
        i000=ix+i00, i100=kx+i00,
        i001=ix+i01, i101=kx+i01,
        i010=ix+i10, i110=kx+i10,
        i011=ix+i11, i111=kx+i11;
        
    const float 
        f000=mx*f00, f100=tx*f00,
        f001=mx*f01, f101=tx*f01,
        f010=mx*f10, f110=tx*f10,
        f011=mx*f11, f111=tx*f11;
    { // DEBUG
        Quat4f fDBG = FFPaul[ i000 ];
        //printf( "GridFF::getForce() u(%g,%g,%g)/[%3i,%3i,%3i] fe(%g,%g,%g|%g) p0(%g,%g,%g) PLQ(%g,%g,%g)\n", u.x,u.y,u.z, grid.n.x,grid.n.y,grid.n.z, fDBG.x,fDBG.y,fDBG.z,fDBG.w, grid.pos0.x,grid.pos0.y,grid.pos0.z, PLQ.x,PLQ.y,PLQ.z );
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
         + (FFelec[ i001 ]*f001) + (FFelec[ i101 ]*f101))*PLQ.z
        ;
    }
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
            const double K=-alphaMorse;
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
                    qp.e+=eM*e;  qp.f.add_mul( dp,  de*e   ); // repulsive part of Morse
                    ql.e+=eM*-2; ql.f.add_mul( dp, -de     ); // attractive part of Morse
                    qe.e+=eQ;    qe.f.add_mul( dp,  eQ*ir2 ); // Coulomb

                }}}
            }
            if(FFPaul)FFPaul[ibuff]=(Quat4f)qp;
            if(FFLond)FFLond[ibuff]=(Quat4f)ql;
            if(FFelec)FFelec[ibuff]=(Quat4f)qe;
        });
    }

    void makeGridFF_omp(int natoms_, Vec3d * apos_, Quat4d * REQs_){
        printf( "GridFF::makeGridFF_omp() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z,  grid.pos0.x,grid.pos0.y,grid.pos0.z );
        if(shifts==0)makePBCshifts( nPBC, lvec );
        const double R2damp=Rdamp*Rdamp;    
        const double K=-alphaMorse;
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
                        for(int ipbc=0; ipbc<npbc; ipbc++ ){
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
                            qp.e+=eM*e;   qp.f.add_mul( dp,  de*e   ); // repulsive part of Morse
                            ql.e+=eM*-2.; ql.f.add_mul( dp, -de     ); // attractive part of Morse
                            qe.e+=eQ;     qe.f.add_mul( dp, -eQ*ir2 ); // Coulomb
                        }
                    }
                    const int ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
                    FFPaul[ibuff]=(Quat4f)qp;
                    FFLond[ibuff]=(Quat4f)ql;
                    FFelec[ibuff]=(Quat4f)qe;
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
        Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4f f = Quat4fZero;
            if(FFPaul ) f.add_mul( FFPaul [ibuff], PLQ.x );
            if(FFLond) f.add_mul( FFLond[ibuff], PLQ.y );
            if(FFelec  ) f.add_mul( FFelec  [ibuff], PLQ.z );
            FF[ibuff] =  f;
        });
    }

    void evalCombindGridFF_CheckInterp( Quat4d REQ, Quat4f * FF ){
        Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );
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


// ============= Debugging and checking

void getEFprofile( int n, Vec3d p0, Vec3d p1, Quat4d REQ, Quat4d* fes, bool bPrint=false){
    if(bPrint){ printf("GridFF::getEFprofile(n=%i,p2{%6.3f,%6.3f,%6.3f},p1{,%6.3f,%6.3f,%6.3f}) \n", n, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z ); };
    Vec3d dp=p1-p0; dp.mul(1./n);
    Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );   //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    for(int i=0; i<n; i++){
        Vec3d  p  = p0 + dp*i;
        Quat4f fe = getForce( p, PLQ, true );
        if(fes)fes[i]=(Quat4d)fe;
        if(bPrint){ printf( "%i %6.3f %6.3f %6.3f %g %g %g\n", i, p.x, p.y, p.z, fe.x,fe.y,fe.z  ); };
    }
}

double checkEFProfileVsNBFF( int n, Vec3d p0, Vec3d p1, const Quat4d& REQ, double tol=1e-2, bool bExit=false, bool bPrint=false, bool bWarn=true ){
    if(bPrint){ printf("GridFF::checkEFProfileVsNBFF(np=%i,natoms=%i,npbc=%i,p2{%6.3f,%6.3f,%6.3f},p1{,%6.3f,%6.3f,%6.3f}) \n", n, natoms,npbc, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z ); };
    if(bPrint){ printf("i   x y z     E  Eref     fx fx_ref      fy fy_ref     fz  fz_ref\n"); }
    double tol2=tol*tol;
    Vec3d dp=p1-p0; dp.mul(1./n);
    Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );   //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    //bool err = false;
    double err2Max=0;
    int imax = -1;
    for(int i=0; i<n; i++){
        Vec3d fref;
        Vec3d  p    = p0 + dp*i;
        Quat4f fe   = getForce          ( p, PLQ, true );
        float  Eref = getMorseQH_PBC_omp( p, REQ, fref );
        float  dE   = fe.e-Eref;
        Vec3f  df   = fe.f - (Vec3f)fref;
        double e2e = dE*dE     /( Eref*Eref + fe.e*fe.e        + 1 );         
        double e2f = df.norm2()/( fe.f.norm2() + fref.norm2() + 1 ); 
        if( e2e>err2Max ){ err2Max=e2e; imax=i; }
        if( e2f>err2Max ){ err2Max=e2f; imax=i; }
        if( (e2e>tol2) || (e2f>tol2) ){
            if(bWarn)printf("WARRNING[%i/%i] dE=%g |dF|=%g p(%6.3f,%6.3f,%6.3f) GridFF(%g,%g,%g|%g)  NBFF(%g.%g,%g|%g)\n", i, n, dE, df.norm(), p.x,p.y,p.z,   fe.x,fe.y,fe.z,fe.w,   fref.x,fref.y,fref.z,Eref  );
            if(bExit){ printf("ERROR in GridFF::checkEFProfileVsNBFF() - GridFF force does not match NBFF reference at test point %i MaxRelativeError=%g => Exit()\n", i, sqrt(err2Max) ); exit(0); }
        } 
        if(bPrint){ printf( "%i    %6.3f %6.3f %6.3f    %g %g   %g %g    %g %g    %g %g\n", i, p.x, p.y, p.z,    fe.e, Eref,    fe.x,fref.x,         fe.y,fref.y,           fe.z,fref.z     ); };
    }
    if(bWarn && (imax>=0) ){
        printf("WARRNING MaxRelativeError=%g at point[%i]\n",  sqrt(err2Max), imax );
        Vec3d fref;
        Vec3d  p    = p0 + dp*imax;
        Quat4f fe   = getForce          ( p, PLQ, true );
        float  Eref = getMorseQH_PBC_omp( p, REQ, fref );
        float  dE   = fe.e-Eref;
        Vec3f  df   = fe.f - (Vec3f)fref;
        double e2e = dE*dE     /(Eref*Eref + fe.e*fe.e + 1);          err2Max=fmax(err2Max,e2e);
        double e2f = df.norm2()/( fe.f.norm2() + fref.norm2() + 1 );  err2Max=fmax(err2Max,e2f);
        printf("WARRNING MaxError at[%i/%i] dE=%g |dF|=%g p(%6.3f,%6.3f,%6.3f) GridFF(%g.%g,%g|%g)  NBFF(%g.%g,%g|%g)\n", imax, n, dE, df.norm(), p.x,p.y,p.z, fe.x,fe.y,fe.z,fe.w,   fref.x,fref.y,fref.z,Eref  );
    }
    return sqrt(err2Max);
}

bool checkZProfilesOverAtom( int ia, int n, double zmin, double zmax, const Quat4d& REQ, double tol=1e-2, bool bExit=true, bool bPrint=true ){
    Vec3d p0=apos[ia]; p0.z=zmin;
    Vec3d p1=p0;       p1.z=zmax;
    double err =  checkEFProfileVsNBFF( n, p0, p1, REQ, tol, false,false );
    if(err>tol){    
        if(bPrint)checkEFProfileVsNBFF( n, p0, p1, REQ, tol,  false, true, false );
        if(bExit){ printf("ERROR in GridFF::checkZProfilesOverAtom(%i) - GridFF force does not match NBFF reference, MaxRelativeError=%g => Exit()\n", ia, err ); exit(0); }
        return true;
    }
    return false;
}

bool evalCheck( int imin=0, int imax=1, bool bExit=true, bool bPrint=true, double tol=1e-2, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.1, 0.}, double dz=0.05 ){
    REQ=Quat4d{ 1.487, 0.02609214441, +0.0, 0.};
    printf( "GridFF::evalCheck() natoms=%i npbc=%i apos=%li REQs=%li shifts=%li \n", natoms, npbc, apos, REQs, shifts );
    _checkNull(shifts)
    _checkNull(REQs)
    _checkNull(apos)
    bool err = false;
    double zmin=grid.pos0.z+1.0;
    double zmax=grid.pos0.z+grid.cell.c.z -1.0;
    zmax=fmin(zmax,10);
    int nz = ((int)((zmax-zmin)/dz)) + 1;
    for( int ia=imin; ia<imax; ia++ ){
        err |= checkZProfilesOverAtom( ia, nz, zmin, zmax, REQ, tol, bExit, bPrint );
    }
    return err;
}

 #ifdef IO_utils_h
    bool tryLoad( const char* fname_Coul, const char* fname_Paul, const char* fname_Lond, bool recalcFF=false ){
        //printf( "GridFF::tryLoad() \n" );
        //printf( "GridFF::tryLoad() fname_Pauli >>%s<< fname_London >>%s<< fname_Coulomb >>%s<< \n", fname_Pauli, fname_London, fname_Coulomb );
        //printDir( "../" );
        //printDir( "./" );
        { FILE* f=fopen( fname_Paul,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Paul); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_Lond,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Lond); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_Coul,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Coul); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        //printf( "GridFF::tryLoad() recalcFF %i \n", recalcFF );
        //printf( "fname_Pauli(%s) fname_London(%s) fname_Coulomb(%s) \n", fname_Pauli, fname_London, fname_Coulomb );
        //int nbyte= grid.getNtot()*sizeof(Vec3d);
        int nbyte= grid.getNtot()*sizeof(Quat4f);
        if( recalcFF ){
            printf( "\nBuilding GridFF for substrate ... (please wait... )\n" );
            makeGridFF_omp( apos_.size(), &apos_[0], &REQs_[0] );
            if(FFPaul) saveBin( fname_Paul,  nbyte, (char*)FFPaul );
            if(FFLond) saveBin( fname_Lond,  nbyte, (char*)FFLond );
            if(FFelec) saveBin( fname_Coul,  nbyte, (char*)FFelec );
        }else{
            if(FFPaul) loadBin( fname_Paul,  nbyte, (char*)FFPaul );
            if(FFLond) loadBin( fname_Lond,  nbyte, (char*)FFLond );
            if(FFelec) loadBin( fname_Coul,  nbyte, (char*)FFelec );
        }
        return recalcFF;
    }
#endif

}; // RigidSubstrate


#endif
