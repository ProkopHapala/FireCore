#ifndef GridFF_h
#define GridFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Forces.h"
#include "MMFFparams.h"

#include "Multipoles.h"

#include "InterpolateTricubic.h"
#include "InterpolateTrilinear.h"
#include "Bspline.h"
#include "Bspline_fit.h"
#include "Bspline_fit_2D.h"
#include "Bspline_fit_3D.h"
#include "VecN.h"

#include "NBFF.h"
#include "EwaldGrid.h"

#include "IO_utils.h"

static bool bDebug__ = 0;

template<typename T>
T sum( int n, T* data, T t ){
    for(int i=0; i<n; i++){ t.add(data[i]); }
    return t;
};


/**
 * Automatically calculates the number of periodic boundary conditions (nPBC) based on minimum length.
 * 
 * @param cell The cell dimensions represented by a 3x3 matrix (Mat3d).
 * @param nPBC The number of periodic boundary conditions in each direction (x, y, z). This parameter will be modified by the function.
 * @param Lmin The minimum length for calculating the number of periodic boundary conditions. Default value is 30.0.
 */
inline void autoNPBC( const Mat3d& cell, Vec3i& nPBC, double Lmin=30.0 ){
    if(nPBC.x!=0){ nPBC.x=(int)Lmin/cell.a.norm(); }
    if(nPBC.y!=0){ nPBC.y=(int)Lmin/cell.b.norm(); }
    if(nPBC.z!=0){ nPBC.z=(int)Lmin/cell.c.norm(); }
    //printf("autoNPBC(): (%i,%i,%i) \n", nPBC.x, nPBC.y, nPBC.z );
}


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

enum class GridFFmod{ Direct=0, LinearFloat=1, LinearDouble=2, HermiteFloat=3, HermiteDouble=4, BsplineFloat=5, BsplineDouble=6 };

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


    Quat4d   *FFPaul_d = 0;
    Quat4d   *FFLond_d = 0;
    Quat4d   *FFelec_d = 0;
    //Quat4f *FFtot    = 0; // total FF is not used since each atom-type has different linear combination

    Quat4d  FEscale{0.0,0.0,0.0,1.0};
    Vec3i   gridN{0,0,0};
    Quat4d   *VPLQH   = 0;
    double   *V_debug = 0;

    double *HHermite_d      = 0; 
    double *Bspline_Pauli   = 0;
    double *Bspline_London  = 0;
    double *Bspline_Coulomb = 0;

    bool bUseEwald=false;
    EwaldGrid* ewald = 0;


    Vec3d  *Bspline_PLQ     = 0;
    Quat4f *Bspline_PLQf    = 0;

    Quat4i cubic_yqis[4];
    Quat4i cubic_xqis[4];


    int iDBG=-1;


    //GridFFmod mode = GridFFmod::LinearFloat;
    GridFFmod mode = GridFFmod::BsplineDouble;
    int perVoxel = 4;

    // ------ ToDo: this should be put inside NBFF
    std::vector<Vec3d>  apos_  ;
    std::vector<Quat4d> REQs_  ;
    std::vector<int>    atypes_;

    //Vec3i nPBC{1,1,0};

    // dipole approximation
    Vec3d dip_p0;
    Vec3d dip;
    double Q;

    double Mpol[10];

    //double Rdamp  =  1.0;
    int iDebugEvalR = 0;
    bool bCellSet    = false;
    bool bSymetrized = false;

    

    double findTop(){ double zmax=-1e+300; for(int i=0;i<natoms; i++){ double z=apos[i].z; if(z>zmax)zmax=z; }; return zmax; }

    void bindSystem(int natoms_, int* atypes_, Vec3d* apos_, Quat4d* REQs_ ){
        natoms=natoms_; atypes=atypes_; apos=apos_; REQs=REQs_;
    }

    void allocateFFs( bool bDouble=false ){
        int ntot = grid.getNtot();
        _realloc( FFPaul, ntot );
        _realloc( FFLond, ntot );
        _realloc( FFelec, ntot );
        if(bDouble){
            //printf( "GridFF::allocateFFs bDouble=%i \n", bDouble );
            _realloc( FFPaul_d, ntot );
            _realloc( FFLond_d, ntot );
            _realloc( FFelec_d, ntot );
        }
    }

    void clear(){
        _dealloc( FFPaul );
        _dealloc( FFLond );
        _dealloc( FFelec );

        _dealloc( FFPaul_d );
        _dealloc( FFLond_d );
        _dealloc( FFelec_d );

        _dealloc( VPLQH   );   
        _dealloc( V_debug );   

        apos_  .clear();
        REQs_  .clear();
        atypes_.clear();

        //NBFF::clear();

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


    // =========================================================
    // ================      Sampling               ============
    // =========================================================

__attribute__((hot))  
inline Quat4f getForce( Vec3d p, const Quat4f& PLQ, bool bSurf=true ) const {
    //printf(  "getForce() p(%g,%g,%g) PLQ(%g,%g,%g,%g) bSurf=%i @FFPaul=%li @FFLond=%li @FFelec=%li \n", p.x,p.y,p.z,  PLQ.x,PLQ.y,PLQ.z,PLQ.w, (long)FFPaul,(long)FFLond,(long)FFelec );
    Vec3d u;
    //p.sub(shift0);
    p.sub(grid.pos0);
    grid.iCell.dot_to( p, u );
    Vec3i n = grid.n;
    //printf( "pos(%g,%g,%g) u(%g,%g,%g)  p0(%g,%g,%g) \n", p.x,p.y,p.z, u.x,u.y,u.z,    grid.pos0.x,grid.pos0.y,grid.pos0.z );
    if( bSurf ){ double ivnz=1./n.z; double uzmax=1-ivnz*1.5; double uzmin=ivnz*0.5; if(u.z<uzmin){ u.z=uzmin; }else if(u.z>uzmax){u.z=uzmax; } } // "clamp" boundrary condition in z-direction
    // ---- PBC condiction ( within outer cell bboundaries ? )
    //u.x=(u.x-(int)u.x)*n.x;
    //u.y=(u.y-(int)u.y)*n.y;
    //u.z=(u.z-(int)u.z)*n.z;
    u.x=(u.x-((int)(u.x+10)-10))*n.x;
    u.y=(u.y-((int)(u.y+10)-10))*n.y;
    u.z=(u.z-((int)(u.z+10)-10))*n.z;
    return Trilinear::fe3f_fe3( (Vec3f)u, n, PLQ.f, FFPaul, FFLond, FFelec );
}
__attribute__((hot))  
inline float addForce( const Vec3d& p, const Quat4f& PLQ, Vec3d& f, bool bSurf=true )const{
    Quat4f fe = getForce( p, PLQ, bSurf ); 
    f.add( (Vec3d)fe.f );   
    return fe.e;
}

__attribute__((hot))  
inline Quat4d getForce_d( Vec3d p, const Quat4d& PLQ, bool bSurf=true ) const {
    Vec3d u;
    //p.sub(shift0);
    p.sub(grid.pos0);
    grid.iCell.dot_to( p, u );
    Vec3i n = grid.n;
    //printf( "pos(%g,%g,%g) u(%g,%g,%g)  p0(%g,%g,%g) \n", p.x,p.y,p.z, u.x,u.y,u.z,    grid.pos0.x,grid.pos0.y,grid.pos0.z );
    if( bSurf ){ double ivnz=1./n.z; double uzmax=1-ivnz*1.5; double uzmin=ivnz*0.5; if(u.z<uzmin){ u.z=uzmin; }else if(u.z>uzmax){u.z=uzmax; } } // "clamp" boundrary condition in z-direction
    // ---- PBC condiction ( within outer cell bboundaries ? )
    //u.x=(u.x-(int)u.x)*n.x;
    //u.y=(u.y-(int)u.y)*n.y;
    //u.z=(u.z-(int)u.z)*n.z;
    u.x=(u.x-((int)(u.x+10)-10))*n.x;
    u.y=(u.y-((int)(u.y+10)-10))*n.y;
    u.z=(u.z-((int)(u.z+10)-10))*n.z;
    return Trilinear::fe3d_fe3( u, n, PLQ.f, FFPaul_d, FFLond_d, FFelec_d );
}
__attribute__((hot))  
inline double addForce_d( const Vec3d& p, const Quat4d& PLQ, Vec3d& f, bool bSurf=true )const{
    Quat4d fe = getForce_d( p, PLQ, bSurf ); 
    f.add( fe.f );   
    return fe.e;
}
__attribute__((hot))  
double addForces_d( int natoms, Vec3d* apos, Quat4d* PLQs, Vec3d* fpos, bool bSurf=true )const{ 
    double E=0;
    for(int ia=0; ia<natoms; ia++){ E+=addForce_d( apos[ia], PLQs[ia], fpos[ia], bSurf ); };
    return E;
}



__attribute__((hot))  
inline Quat4d getForce_Tricubic( Vec3d p, const Quat4d& PLQH, bool bSurf=true ) const {
    Vec3d u;
    //p.sub(shift0);
    p.sub(grid.pos0);
    grid.iCell.dot_to( p, u );
    u.x=(u.x-((int)(u.x+10)-10))*grid.n.x-1;
    u.y=(u.y-((int)(u.y+10)-10))*grid.n.y-1;
    if(u.z<0.0){ u.z=0.0; }else if(u.z>1.0){ u.z=1.0; }; u.z=u.z*grid.n.z-1;
    return Spline_Hermite::fe3d_v4( PLQH, u, gridN, VPLQH ) * FEscale;
}
__attribute__((hot))  
inline float addForce_Tricubic( const Vec3d& p, const Quat4d& PLQ, Vec3d& f, bool bSurf=true )const{
    Quat4d fe = getForce_Tricubic( p, PLQ, bSurf ); 
    f.add( (Vec3d)fe.f );   
    return fe.e;
}


inline int fold_cubic2( int i, int n )const{
    //i+=n*100;
    return (i+(n<<8) )%n;
    //if( i==0 )[[unlikely]]{  return n-1; } else if (i>n)[[unlikely]]{ return i-n; } else [[likely]] {  return i; };
}

__attribute__((hot))  
inline Quat4d getForce_HHermit( Vec3d p, const Quat4d& PLQH, bool bSurf=true ) const {
    //printf( "GridFF::getForce_HHermit() p(%g,%g,%g)\n", p.x,p.y,p.z );
    Vec3d t;
    //p.sub(shift0);
    p.sub(grid.pos0);
    grid.diCell.dot_to( p, t );
    Vec3d inv_dg2{ -grid.diCell.xx, -grid.diCell.yy, -grid.diCell.zz };
    int ix=(int)t.x;
    int iy=(int)t.y;
    int iz=(int)t.z;

    if(t.x<0){ ix=ix-1; };
    if(t.y<0){ iy=iy-1; };
    const int ix_ = fold_cubic2( ix, grid.n.x );
    const int iy_ = fold_cubic2( iy, grid.n.y );
    //double dx=t.x-ix;
    //double dy=t.y-iy;
    //printf( "t.x %g dx %g ix %i \n", t.x, dx, ix );
    //if(dx<0){ dx=dx-1; };
    //if(dy<0){ dy=dy-1; };
    //Quat4d fe = Spline_Hermite::fe3d_comb3( Vec3d{t.x-ix,t.y-iy,t.z-iz}, Vec3i{ix_,iy_,iz}, grid.n, (Vec2d*)HHermite_d, PLQH.f );
    Quat4d fe = Spline_Hermite::fe3d_comb3( Vec3d{t.x-ix,t.y-iy,t.z-iz}, Vec3i{ix_,iy_,iz}, gridN, (Vec2d*)HHermite_d, PLQH.f );
    //{ const int i0 = (iz+grid.n.z*(iy+grid.n.y*ix))*6; printf( "GridFF::getForce_HHermit() p(%g,%g,%g) i(%i,%i,%i) VPLQ(%g,%g,%g) PLQ(%g,%g,%g) \n", p.x,p.y,p.z, ix,iy,iz, HHermite_d[i0+0], HHermite_d[i0+2], HHermite_d[i0+4], PLQH.x,PLQH.y,PLQH.z ); }
    fe.f.mul(inv_dg2);
    return fe;
}
__attribute__((hot))  
inline float addForce_HHermit( const Vec3d& p, const Quat4d& PLQ, Vec3d& f, bool bSurf=true )const{
    Quat4d fe = getForce_HHermit( p, PLQ, bSurf ); 
    f.add( (Vec3d)fe.f );   
    return fe.e;
}

__attribute__((hot))  
inline Quat4d getForce_Bspline( Vec3d p, const Quat4d& PLQH, bool bSurf=true ) const {
    // printf( "GridFF::getForce_Bspline() p(%8.4f,%8.4f,%8.4f) PLQH(%8.4f,%8.4f,%8.4f,%8.4f)\n", p.x,p.y,p.z, PLQH.x,PLQH.y,PLQH.z,PLQH.w );
    
    Vec3d t;
    //p.sub(shift0);
    Vec3d adjust {0.0,0.0,1.5};
    // p.sub(adjust);
    p.sub(grid.pos0);
    // p.sub(grid.pos0+adjust);
    // p.add(adjust);
    // p.add(grid.pos0);
    // p.add(grid.pos0+adjust);
    // p.add(grid.pos0+shift0);
    
    grid.diCell.dot_to( p, t );
    Vec3d inv_dg2{ -grid.diCell.xx, -grid.diCell.yy, -grid.diCell.zz };

    if(verbosity>1)printf("DEBUG: Bspline | Grid pos0=%.3f,%.3f,%.3f | Atom pos=%.3f,%.3f,%.3f \n", grid.pos0.x,grid.pos0.y,grid.pos0.z, p.x,p.y,p.z);
    //Quat4d fe = Quat4dZero;
    
    Quat4d fe = Bspline::fe3d_pbc_comb3( t, grid.n, Bspline_PLQ, PLQH.f, cubic_xqis, cubic_yqis ); 

    //printf( "GridFF::getForce_Bspline() p(%8.4f,%8.4f,%8.4f) g0(%8.4f,%8.4f,%8.4f) PLQH(%8.4f,%8.4f,%8.4f,%8.4f) fe(%g,%g,%g,%g)\n", p.x,p.y,p.z, grid.pos0.x,grid.pos0.y,grid.pos0.z, PLQH.x,PLQH.y,PLQH.z,PLQH.w, fe.x,fe.w,fe.z,fe.w ); 

    // Quat4d fe = Bspline::fe3d( t, grid.n, Bspline_Pauli   )*PLQH.x
    //           + Bspline::fe3d( t, grid.n, Bspline_London  )*PLQH.y  
    //           + Bspline::fe3d( t, grid.n, Bspline_Coulomb )*PLQH.z;  


    //Quat4d fe = Bspline::fe3d( t, gridN, Bspline_Coulomb );
    //printf( "GridFF::getForce_Bspline() p(%g,%g,%g) fe(%g,%g,%g,%g)\n", p.x,p.y,p.z, fe.x,fe.y,fe.z,fe.w );

    //Quat4d fe = Bspline::fe3d( t, gridN, Bspline_Coulomb );
    //Quat4d fe = Bspline::fe3d( Vec3d{t.z,t.y,t.x}, Vec3i{gridN.z,gridN.y,gridN.x}, Bspline_Coulomb );

    // int i = ((int)t.z) + grid.n.z*( ((int)t.y) + ((int)t.x)*grid.n.y );
    // fe.w = Bspline_PLQ[i].dot( PLQH.f );
    // printf( "GridFF::getForce_Bspline() p(%8.4f,%8.4f,%8.4f) ixyz(%4i,%4i,%4i,%8i) BsplinePLQ(%g,%g,%g)\n", p.x,p.y,p.z,  ((int)t.z), ((int)t.y), ((int)t.x), i,    Bspline_PLQ[i].x,Bspline_PLQ[i].y,Bspline_PLQ[i].z  );
    // fe.w = Bspline_Coulomb[i];
    // printf( "GridFF::getForce_Bspline() p(%g,%g,%g) t(%g,%g,%g) Gs[%i]=%g \n", p.x,p.y,p.z,  t.x,t.y,t.z,  i, fe.w );

    fe.f.mul(inv_dg2);
    return fe;
}
__attribute__((hot))  
inline float addForce_Bspline( const Vec3d& p, const Quat4d& PLQ, Vec3d& f, bool bSurf=true )const{
    Quat4d fe = getForce_Bspline( p, PLQ, bSurf ); 
    f.add( (Vec3d)fe.f );   
    return fe.e;
}

__attribute__((hot))  
double addForces( int natoms, Vec3d* apos, Quat4f* PLQs, Vec3d* fpos, bool bSurf=true )const{ 
    double E=0;
    for(int ia=0; ia<natoms; ia++){ E+=addForce( apos[ia], PLQs[ia], fpos[ia], bSurf ); };
    return E;
}
__attribute__((hot))  
inline void addForce( const Vec3d& pos, const Quat4f& PLQ, Quat4f& fe ) const {
    //printf( "GridFF::addForce() \n" );
    //printf( "GridFF::addForce() pointers(%li,%li,%li)\n", FFPaul, FFLond, FFelec );
    Vec3f gpos; grid.cartesian2grid(pos, gpos);
    // printf("DEBUG: Grid pos0.z=%.3f | Atom pos.z=%.3f\n", grid.pos0.z, pos.z);
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




    __attribute__((hot))  
    inline void addForce_surf( Vec3d pos, const Quat4f PLQ, Quat4f& f ) const {
        //pos.add( shift0 );
        if     ( pos.z > grid.cell.c.z ){ pos.z = grid.dCell.c.z*-0.1 + grid.cell.c.z; }
        else if( pos.z < 0             ){ pos.z = grid.dCell.c.z* 0.1;                 }
        return addForce( pos, PLQ, f );
    }

    __attribute__((hot))  
    inline double eval( int n, const Vec3d* ps, const Quat4f* PLQs, Vec3d* fs, bool bSurf=false ) const {
        double E=0;
        //printf("GridFF::eval() n %i ps %li PLQs %li \n", n,  (long)ps,  (long)PLQs );
        if(bSurf){ for(int i=0; i<n; i++){ Quat4f fe=Quat4fZero;                  addForce_surf( ps[i], (Quat4f)PLQs[i], fe );  fs[i].add( (Vec3d)fe.f ); E+=fe.e; } }
        else     { for(int i=0; i<n; i++){ Quat4f fe=Quat4fZero; bDebug__=(i==0); addForce     ( ps[i], (Quat4f)PLQs[i], fe );  fs[i].add( (Vec3d)fe.f ); E+=fe.e;  
            //if(i==0)printf("CPU[0] apos(%g,%g,%g) PLQs[0](%g,%g,%g|%g) \n", n, ps[i].x,ps[i].y,ps[i].z,  PLQs[i].x,PLQs[i].y,PLQs[i].z,alpha ); 
        } }
        return E;
    }

    inline double addAtom( const Vec3d& pos, const Quat4d& PLQ, Vec3d& fout )const{
        Quat4d fed;
        switch( mode ){
            // void evalGridFFPoint( int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3d pos, Quat4d& qp, Quat4d& ql, Quat4d& qe )const{
            // case GridFFmod::Direct       :{ Quat4d qp,ql,qe; evalGridFFPoint( apos.size(), apos.data(), REQs_.data(), pos, qp, ql, qe ); fed = qp*PLQ.x + ql*PLQ.y + qe*PLQ.z;
            case GridFFmod::Direct       :{ Quat4d qp,ql,qe; evalGridFFPoint( apos_.size(), apos_.data(), REQs_.data(), pos, qp, ql, qe ); fed = qp*PLQ.x + ql*PLQ.y + qe*PLQ.z; 
                // printf( "GridFF::addAtom( E(%g,%g,%g|%g)  PLQ(%g,%g,%g)  pos(%g,%g,%g) \n", qp.e, ql.e, qe.e, fed.e,  PLQ.x,PLQ.y,PLQ.z, pos.x,pos.y,pos.z  );
            } break;
            case GridFFmod::LinearFloat  :{ fed=(Quat4d)getForce( pos, (Quat4f)PLQ ); }break;
            case GridFFmod::LinearDouble :{ fed=getForce_d      ( pos, PLQ);          }break;
            case GridFFmod::HermiteDouble:{ fed=getForce_HHermit( pos, PLQ );         }break;
            case GridFFmod::BsplineDouble:{ fed=getForce_Bspline( pos, PLQ );         }break;
                // printf( "@@@@@@@@@@GridFF::addAtom( E(%g,%g,%g|%g)  PLQ(%g,%g,%g)  pos(%g,%g,%g) \n", qp.e, ql.e, qe.e, fed.e,  PLQ.x,PLQ.y,PLQ.z, pos.x,pos.y,pos.z );
            //case GridFFmod::HermiteFloat:  { }break;
            //cast GridFFmod::BSplineFloat:  { }break;
            default: { printf("ERROR GridFF::eval() mode=%i NOT IMPLEMENTED !!! \n", mode); exit(0); }
        }
        //printf("GridFF::addAtom() FE(%10.5f,%10.5f,%10.5f|%10.5f) pos(%7.5f,%7.5f,%7.5f) PLQ(%7.5f,%10.6f,%7.3f) \n",   fed.x,fed.y,fed.z,fed.w,     pos.x,pos.y,pos.z,  PLQ.x,PLQ.y,PLQ.z );
        fout.add( fed.f );
        //fout.sub( fed.f );
        return fed.e;
    }

    Vec3d findIso(double isoval, Vec3d p0, Vec3d p1, const Quat4d PLQ, double xtol = 0.01, bool bPrintError=true )const{
        //printf( "GridFF::findIso() iso=%g p0(%6.3f,%6.3f,%6.3f) p1(%6.3f,%6.3f,%6.3f) PLQ(%g,%g,%g,%g) xtol=%g \n", isoval, p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, PLQ.x,PLQ.y,PLQ.z,PLQ.w, xtol );
        Vec3d fout;  // Force output vector (unused in binary search)
        // Evaluate energy at the endpoints
        double f0 = addAtom(p0, PLQ, fout)-isoval;
        double f1 = addAtom(p1, PLQ, fout)-isoval;
        //printf( "GridFF::findIso() iso=%g p0(%6.3f,%6.3f,%6.3f)f0=%fg p1(%6.3f,%6.3f,%6.3f)f1=%g PLQ(%g,%g,%g,%g) xtol=%g \n", isoval, p0.x,p0.y,p0.z,f0, p1.x,p1.y,p1.z,f1, PLQ.x,PLQ.y,PLQ.z,PLQ.w, xtol );
        if( f0*f1 > 0.0) {
            if(bPrintError) printf("ERROR GridFF::findIso()  E1(%10.5e)*E2(%10.5e) > 0.0 |  p0(%8.4f,%8.4f,%8.4f) p1(%8.4f,%8.4f,%8.4f) PLQ(%10.5e,%10.5e,%10.5e,%10.5e) \n", f0,f1,  p0.x,p0.y,p0.z, p1.x,p1.y,p1.z, PLQ.x,PLQ.y,PLQ.z,PLQ.w );
            //exit(0);
            //return p0;
            return Vec3dNAN;
        }
        double sgn = (f0 > 0.0) ? 1.0 : -1.0;
        double r2tol = xtol*xtol;
        Vec3d pmid;
        int iter=0;

        while (  (p0-p1).norm2() > r2tol ) {
            pmid  = (p0 + p1) * 0.5;
            double fmid = addAtom(pmid, PLQ, fout);
            if ( (fmid-isoval)*sgn < 0.0 ) { p1 = pmid; } 
            else                           { p0 = pmid; }
            //printf( "p0.z=%6.3f p1.z=%6.3f fmid=%g \n", p0.z, p1.z, fmid );
            iter++;
        }
        return pmid;
    }

    // =========================================================
    // ================       Grid Preparation      ============
    // =========================================================

    void init( Vec3i n, Mat3d cell, Vec3d pos0, bool bDouble=false ){
        grid.n     = n;
        grid.setCell(cell);
        grid.pos0  = pos0;
        //zmax = cell.c.z;
        allocateFFs( bDouble );
    }

    void setAtoms( int natoms_, Vec3d * apos_, Quat4d * REQs_ ){
        natoms = natoms_;
        //atypes = new int  [natoms];
        apos   = apos_;
        REQs  = REQs_;
    }

    __attribute__((hot))  
    void evalGridFFel(int natoms, const Vec3d * apos, const Quat4d * REQs, Vec3d * FF )const{
        //interateGrid3D( Vec3d{0.0,0.0,0.0}, grid.n, grid.dCell, [=](int ibuff, Vec3d p)->void{
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Vec3d f = Vec3dZero;
            for(int ia=0; ia<natoms; ia++){ addAtomicForceQ( p-apos[ia], f, REQs[ia].z ); }
            FF[ibuff]=f;
        });
    }


    //__attribute__((pure))
    __attribute__((hot))
    void evalGridFFPoint( int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3d pos, Quat4d& qp, Quat4d& ql, Quat4d& qe )const{
        //const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
        const double R2damp=Rdamp*Rdamp;    
        const double K=-alphaMorse;
        qp = Quat4dZero; 
        ql = Quat4dZero; 
        qe = Quat4dZero;
        //#pragma omp for simd
        for(int ia=0; ia<natoms_; ia++){
            const Vec3d dp0   = pos - apos_[ia]+shift0;
            const Quat4d REQi = REQs_[ia];
            //if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }              
            for(int ipbc=0; ipbc<npbc; ipbc++ ){
                const Vec3d  dp = dp0 + shifts[ipbc];
                //Vec3d  dp = dp0;
                const double r2     = dp.norm2();
                const double r      = sqrt(r2 + 1e-32);
                // ---- Coulomb
                const double ir2    = 1/(r2+R2damp);
                const double eQ     = COULOMB_CONST*REQi.z*sqrt(ir2);
                // ----- Morse
                const double e      = exp( K*(r-REQi.x) );
                const double eM     = e*REQi.y;
                const double de     = 2*K*eM/r;                    
                // --- store
                qp.e+=eM*e;   qp.f.add_mul( dp, -de*e   ); // repulsive part of Morse
                ql.e+=eM*-2.; ql.f.add_mul( dp,  de     ); // attractive part of Morse
                qe.e+=eQ;     qe.f.add_mul( dp,  eQ*ir2 ); // Coulomb
                //qp.e += exp( -r2/0.16 ); // Debug
            }
        }
        //printf( "evalGridFFPoint() E_PLQ(%g,%g,%g) pos(%g,%g,%g) \n", qp.e, ql.e, qe.e, pos.x, pos.y, pos.z );
        //const int ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
        //FFPaul[ibuff]=(Quat4f)qp;
        //FFLond[ibuff]=(Quat4f)ql;
        //FFelec[ibuff]=(Quat4f)qe;
    }

    __attribute__((hot))
    void evalGridFFPoint_Mors( int npbc, const Vec3d* shifts, int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3d pos, Quat4d& qp, Quat4d& ql )const{
        //const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
        //const double R2damp=Rdamp*Rdamp;    
        const double K=-alphaMorse;
        qp = Quat4dZero; 
        ql = Quat4dZero; 

        
        //#pragma omp for simd
        for(int ia=0; ia<natoms_; ia++){
            const Vec3d dp0   = pos - apos_[ia]+shift0; 
            const Quat4d REQi = REQs_[ia];
            // if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }              
            for(int ipbc=0; ipbc<npbc; ipbc++ ){
                const Vec3d  dp = dp0 + shifts[ipbc];
                //Vec3d  dp = dp0;
                const double r2     = dp.norm2();
                const double r      = sqrt(r2 + 1e-32);
                // ----- Morse
                const double e      = exp( K*(r-REQi.x) );
                const double eM     = e*REQi.y;
                const double de     = 2*K*eM/r;          
                //if( (iDBG>=0) && (ia==0) && (ipbc==npbc/2) ){ printf( "evalGridFFPoint_Mors() dp(%8.4f,%8.4f,%8.4f|r=%8.4f) e=%g EPaul=%g ELond=%g alphaMorse=%g R0=%g E0=%g shift(%8.4f,%8.4f,%8.4f) apos(%8.4f,%8.4f,%8.4f) \n", dp.x,dp.y,dp.z, r, e,  eM*e, eM*-2.,  alphaMorse, REQi.x, REQi.y, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z, apos_[ia].x,apos_[ia].y,apos_[ia].z );  };          
                // --- store
                qp.e+=eM*e;   qp.f.add_mul( dp, -de*e   ); // repulsive part of Morse
                ql.e+=eM*-2.; ql.f.add_mul( dp,  de     ); // attractive part of Morse
            }
        }
    }

    __attribute__((pure))
    __attribute__((hot))
    Quat4d evalGridFFPoint_Coul( int npbc, const Vec3d* shifts, int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3d pos ) const {
        //const double R2damp=Rdamp*Rdamp;    
        const double R2damp=1e-32; 
        //const double K=-alphaMorse;
        Quat4d qe     = Quat4dZero;
        //#pragma omp for simd
        for(int ia=0; ia<natoms_; ia++){
            const Vec3d dp0   = pos - apos_[ia]+shift0;
            const Quat4d REQi = REQs_[ia];
            //if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }              
            for(int ipbc=0; ipbc<npbc; ipbc++ ){
                const Vec3d  dp = dp0 + shifts[ipbc];
                const double r2     = dp.norm2();
                const double ir2    = 1/(r2+R2damp);
                const double eQ     = COULOMB_CONST*REQi.z*sqrt(ir2);
                qe.e+=eQ;     qe.f.add_mul( dp,  eQ*ir2 ); // Coulomb
            }
        }
        return qe;
    }

    __attribute__((pure))
    __attribute__((hot))
    Quat4d evalGridFFPoint_REQ( Quat4d REQ, int npbc, const Vec3d* shifts, int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3d pos ) const {
        //const double R2damp=Rdamp*Rdamp;    
        const double R2damp=1e-32; 
        const double K=-alphaMorse;
        Quat4d qe     = Quat4dZero;
        //#pragma omp for simd
        for(int ia=0; ia<natoms_; ia++){
            const Vec3d dp0   = pos - apos_[ia]+shift0;
            Quat4d REQ; combineREQ( REQ, REQs_[ia], REQ );
            //if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }              
            for(int ipbc=0; ipbc<npbc; ipbc++ ){
                const Vec3d  dp = dp0 + shifts[ipbc];
                const double r2     = dp.norm2();
                const double ir2    = 1/(r2+R2damp);
                const double  r     = sqrt( r2   );
                const double eQ     = COULOMB_CONST*REQ.z*sqrt(ir2);

                // Morse:
                const double  e  = exp( K*(r-REQ.x) );
                const double  Ae = REQ.y*e;
                //const double  He  = e * ( (REQ.w<0) ? REQH.w*REQH.y : 0.0 );  // H-bond correction
                const double E =  Ae*(e - 2)            + eQ;
                const double F = (Ae*(e - 1)*2 )*-K/r   - eQ*ir2;
                qe.f.set_mul( dp, F );
                qe.e = E;
            }
        }
        return qe;
    }

    Quat4d evalMorsePBC_PLQ(  Vec3d pi, Quat4d PLQH, int natoms, const Vec3d * apos, const Quat4d * REQs )const{
        //printf( "GridFF::evalMorsePBC() debug fi(%g,%g,%g) REQi(%g,%g,%g)\n",  fi.x,fi.y,fi.z, REQi.x,REQi.y,REQi.z,REQi.w  );
        Quat4d qp,ql,qe;
        evalGridFFPoint( natoms, apos, REQs, pi, qp, ql, qe );
        return qp*PLQH.x + ql*PLQH.y + qe*PLQH.z;
    }
    Quat4d evalMorsePBC_PLQ_sym( Vec3d  pi, Quat4d  PLQH )const{ return evalMorsePBC_PLQ( pi, PLQH, apos_.size(), &apos_[0], &REQs_[0] ); }
    // Quat4d evalMorsePBCatoms_PLQ_sym( int na, Vec3d* ps, Quat4d* REQs, Vec3d* forces ){
    //     double E = 0;
    //     for(int ia=0; ia<na; ia++){ evalMorsePBC_PLQ_sym( ps[ia], REQs[ia], forces[ia] ); };
    //     return E;
    // }

    void evalAtPoints_REQ( int n, Vec3d* ps, Quat4d* FFout, Quat4d REQ, int natoms_, Vec3d * apos_, Quat4d * REQs_, Vec3i* nPBC_=0 ){
        Vec3i nPBC__ = nPBC;
        if(nPBC_) nPBC__=*nPBC_;
        //printf( "GridFF::evalAtPoints_Split() n=%i natoms_=%i bCoul=%i bMors=%i nPBC_Coul(%i,%i,%i) nPBC_Coul(%i,%i,%i)\n", n, natoms_, bCoul,bMors, nPBC_coul.x,nPBC_coul.y,nPBC_coul.z, nPBC_mors.x,nPBC_mors.y,nPBC_mors.z );
        int i=0;
        Vec3d* shifts=0;
        int npbc=makePBCshifts_( nPBC__, lvec, shifts );;
        printf( "GridFF::evalAtPoints_REQ() n=%i natoms_=%i REQ(%g,%g,%g,%g) nPBC(%i,%i,%i|%i)\n", n, natoms_, REQ.x,REQ.y,REQ.z,REQ.w, nPBC.x,nPBC.y,nPBC.z, npbc );
        #pragma omp parallel for shared(i)
        for(int i=0; i<n; i++){
            FFout[i] = evalGridFFPoint_REQ( REQ, npbc, shifts,  natoms_, apos_, REQs_, ps[i] );
        }
        //double T = (getCPUticks()-t0); printf( "evalAtPoints_Split(n=%i) time=%g[MTick] %g[kTick/point]\n", n, T*(1.e-6), (T*1e-3)/n );
        delete[] shifts;
    }
    void evalAtPoints_REQ( int n, Vec3d* ps, Quat4d* FFout, Quat4d REQH, Vec3i* nPBC=0 ){ evalAtPoints_REQ( n, ps, FFout, REQH, natoms, apos, REQs, nPBC ); };

    void evalAtPoints( int n, const Vec3d* ps, Quat4d* FFout, Quat4d PLQH, int natoms_, const Vec3d * apos_, const Quat4d * REQs_, Vec3i* nPBC=0 ){
        //printf( "GridFF::evalAtPoints() n=%i natoms_=%i \n", n, natoms_ );
        int i=0;
        // backup 
        Vec3d* shifts_=shifts;
        int    npbc_ = npbc;
        bool bUserPBC = (nPBC!=0);
        if( bUserPBC ){
            shifts=0;
            npbc = makePBCshifts_( *nPBC, lvec, shifts ); 
        }
        #pragma omp parallel for shared(i)
        for(int i=0; i<n; i++){
            Quat4d qp,ql,qe;
            //printf( "ps[%i](%7.3f,%7.3f,%7.3f)\n", i, ps[i].x,ps[i].y,ps[i].z  );
            evalGridFFPoint( natoms_, apos_, REQs_, ps[i], qp, ql, qe );
            FFout[i] = qp*PLQH.x + ql*PLQH.y + qe*PLQH.z;
        }
        // restore
        if( bUserPBC ){
            npbc  = npbc_;
            shifts = shifts_;
        }
    }
    //void evalAtPoints( int n, Vec3d* ps, Quat4d* FFout, Quat4d PLQH ){ evalAtPoints( n, ps, FFout, PLQH, natoms, apos, REQs ); };
    void evalAtPoints( int n, const Vec3d* ps, Quat4d* FFout, Quat4d PLQH, Vec3i* nPBC=0 ){ evalAtPoints( n, ps, FFout, PLQH, apos_.size(), apos_.data(), REQs_.data(), nPBC ); };

    void evalAtPoints_Split( int n, Vec3d* ps, Quat4d* FFout, Quat4d PLQH, int natoms_, Vec3d * apos_, Quat4d * REQs_, Vec3i* nPBC=0 ){
        Vec3i nPBC_mors{4,4,0};
        Vec3i nPBC_coul{25,25,0}; if(nPBC) nPBC_coul=*nPBC;
        bool bCoul = fabs(PLQH.z) > 1e-32;
        bool bMors = fmax( fabs(PLQH.x), fabs(PLQH.y) ) > 1e-32;
        //printf( "GridFF::evalAtPoints_Split() n=%i natoms_=%i bCoul=%i bMors=%i nPBC_Coul(%i,%i,%i) nPBC_Coul(%i,%i,%i)\n", n, natoms_, bCoul,bMors, nPBC_coul.x,nPBC_coul.y,nPBC_coul.z, nPBC_mors.x,nPBC_mors.y,nPBC_mors.z );
        int i=0;
        Vec3d* shift_coul=0;
        Vec3d* shift_mors=0;
        int npbc_coul=0;
        int npbc_mors=0;
        if(bCoul){ npbc_coul = makePBCshifts_( nPBC_coul, lvec, shift_coul ); } // Coulomb converge much more slowly with nPBC
        if(bMors){ npbc_mors = makePBCshifts_( nPBC_mors, lvec, shift_mors ); }
        //#pragma omp parallel for shared(i)
        //long t0 = getCPUticks();
        //printf( "GridFF::evalAtPoints_Split() npbc_coul=%i npbc_mors=%i \n", npbc_coul, npbc_mors );
        printf( "GridFF::evalAtPoints_Split() n=%i natoms_=%i PLQH(%g,%g,%g,%g) bCoul=%i bMors=%i nPBC_Coul(%i,%i,%i|%i) nPBC_mors(%i,%i,%i|%i)\n", n, natoms_, bCoul,bMors,  PLQH.x,PLQH.y,PLQH.z,PLQH.w, nPBC_coul.x,nPBC_coul.y,nPBC_coul.z,npbc_coul,   nPBC_mors.x,nPBC_mors.y,nPBC_mors.z,npbc_mors );
        #pragma omp parallel for shared(i)
        for(int i=0; i<n; i++){
            Quat4d qp=Quat4dZero,ql=Quat4dZero,qe=Quat4dZero;
            if(bMors){    evalGridFFPoint_Mors( npbc_mors, shift_mors,  natoms_, apos_, REQs_, ps[i], qp, ql ); }
            if(bCoul){ qe=evalGridFFPoint_Coul( npbc_coul, shift_coul,  natoms_, apos_, REQs_, ps[i]         ); }
            FFout[i] = qp*PLQH.x + ql*PLQH.y + qe*PLQH.z;
        }
        //double T = (getCPUticks()-t0); printf( "evalAtPoints_Split(n=%i) time=%g[MTick] %g[kTick/point]\n", n, T*(1.e-6), (T*1e-3)/n );
        delete[] shift_coul;
        delete[] shift_mors;
    }
    //void evalAtPoints_Split( int n, Vec3d* ps, Quat4d* FFout, Quat4d PLQH, Vec3i* nPBC=0 ){ evalAtPoints_Split( n, ps, FFout, PLQH, apos_.size(), apos_.data(), REQs_.data(), nPBC ); };
    void evalAtPoints_Split( int n, Vec3d* ps, Quat4d* FFout, Quat4d PLQH, Vec3i* nPBC=0 ){ evalAtPoints_Split( n, ps, FFout, PLQH, natoms, apos, REQs, nPBC ); };


    __attribute__((hot))  
    void makeGridFF_omp(int natoms_, Vec3d * apos_, Quat4d * REQs_ ){
        printf( "GridFF::makeGridFF_omp()  nPBC(%i,%i,%i) npbc=%i natoms_=%i pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z,  npbc, natoms_,  grid.pos0.x,grid.pos0.y,grid.pos0.z );        
        if(shifts==0)makePBCshifts( nPBC, lvec );
        const double R2damp=Rdamp*Rdamp;    
        const double K=-alphaMorse;

        //for(int i=0; i<npbc; i++ ){ printf( "shifts[%i](%g,%g,%g) \n", i, shifts[i].x, shifts[i].y, shifts[i].z  ); };
        //for(int ia=0; ia<natoms_; ia++)    { printf( "apos[%i] pos(%g,%g,%g) REQ(%g,%g,%g,%g) \n", ia, apos_[ia].x, apos_[ia].y, apos_[ia].z, REQs_[ia].x, REQs_[ia].y, REQs_[ia].z, REQs_[ia].w ); };
        int ix=0,iy=0,iz=0;
        #pragma omp parallel for shared(ix,iy,iz,FFPaul,FFLond,FFelec) collapse(3)
        for ( iz=0; iz<grid.n.z; iz++ ){
            for ( iy=0; iy<grid.n.y; iy++ ){
                for ( ix=0; ix<grid.n.x; ix++ ){
                    Quat4d qp,ql,qe;
                    const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
                    evalGridFFPoint( natoms_, apos_, REQs_, pos, qp, ql, qe );
                    const int ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
                    FFPaul[ibuff]=(Quat4f)qp;
                    FFLond[ibuff]=(Quat4f)ql;
                    FFelec[ibuff]=(Quat4f)qe;                   
                }
            }
        }
    }
    void makeGridFF(){ makeGridFF_omp(natoms,apos,REQs); }

    __attribute__((hot))  
    void makeGridFF_omp_d(int natoms_, Vec3d * apos_, Quat4d * REQs_ ){
        printf( "GridFF::makeGridFF_omp_d() nPBC(%i,%i,%i) pos0(%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z,  grid.pos0.x,grid.pos0.y,grid.pos0.z );        
        if(shifts==0)makePBCshifts( nPBC, lvec );
        const double R2damp=Rdamp*Rdamp;    
        const double K=-alphaMorse;
        int ix=0,iy=0,iz=0;
        #pragma omp parallel for shared(ix,iy,iz,FFPaul,FFLond,FFelec) collapse(3)
        for ( iz=0; iz<grid.n.z; iz++ ){
            for ( iy=0; iy<grid.n.y; iy++ ){
                for ( ix=0; ix<grid.n.x; ix++ ){
                    Quat4d qp,ql,qe;
                    const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
                    evalGridFFPoint( natoms_, apos_, REQs_, pos, qp, ql, qe );
                    const int ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
                    FFPaul_d[ibuff]=qp;
                    FFLond_d[ibuff]=ql;
                    FFelec_d[ibuff]=qe;
                }
            }
        }
    }
    void makeGridFF_d(){ makeGridFF_omp_d(natoms,apos,REQs); }

    __attribute__((hot))  
    void makeGridFF_Hherm_d(int natoms_, Vec3d * apos_, Quat4d * REQs_ ){
        //printf( "GridFF::makeGridFF_Hherm_d() grid.n(%i,%i,%i) nPBC(%i,%i,%i)  npbc=%i natoms=%i pos0(%g,%g,%g) K=%g Rdamp=%g \n",   grid.n.x,grid.n.y,grid.n.z,  nPBC.x,nPBC.y,nPBC.z,  npbc, natoms_,  grid.pos0.x,grid.pos0.y,grid.pos0.z, alphaMorse, Rdamp );   
        //printf( "GridFF::makeGridFF_Hherm_d() grid.n(%i,%i,%i) gridN(%i,%i,%i) \n", grid.n.x,grid.n.y,grid.n.z, gridN.x,gridN.y,gridN.z ); 
        //grid.printCell();    
        
        //if(HHermite_d==0) _realloc( HHermite_d, grid.n.totprod()*6 );
        if(HHermite_d==0) _realloc( HHermite_d, gridN.totprod()*6 );
        if(shifts==0)makePBCshifts( nPBC, lvec );

        double dz = -grid.dCell.c.z;

        //const int ntot = grid.n.totprod();
        //const int nzy  = grid.n.z*grid.n.y;
        const int ntot = gridN.totprod();
        const int nzy  = gridN.z*gridN.y;

        int i=0;
        #pragma omp parallel for shared(i,HHermite_d)
        for(int i=0; i<ntot; i++){
            //int iz =  i%grid.n.z;
            //int iy = (i/grid.n.z)%grid.n.y;
            int iz =  i%gridN.z;
            int iy = (i/gridN.z)%gridN.y;
            int ix =  i/nzy;
                    Quat4d qp,ql,qe;
                    //const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
                    const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*(iy-1) + grid.dCell.a*(ix-1);
                    evalGridFFPoint( natoms_, apos_, REQs_, pos, qp, ql, qe );
                    //const int ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
                    const int ibuff = iz + gridN.z*( iy + gridN.y * ix );
                    int i6 = ibuff*6;
                    HHermite_d[i6+0] = qp.w; HHermite_d[i6+1] = qp.z*dz;
                    HHermite_d[i6+2] = ql.w; HHermite_d[i6+3] = ql.z*dz;
                    HHermite_d[i6+4] = qe.w; HHermite_d[i6+5] = qe.z*dz;     
        }
        printf( "GridFF::makeGridFF_Hherm_d() DONE \n" );
    }
    void makeGridFF_Hherm_d(){ makeGridFF_Hherm_d(natoms,apos,REQs); }

    __attribute__((hot))  
    void evalBsplineRef(int natoms_, Vec3d * apos_, Quat4d * REQs_, double* VPaul, double* VLond, double* VCoul ){
        Vec3i nPBC_coul=Vec3i{25,25,0};
        Vec3i nBPC_mors=nPBC;
        bool bEvalCoulomb = !bUseEwald;
        printf( "GridFF::evalBsplineRef() natoms=%i Rdamp=%g alphaMorse=%g bEvalCoulomb=%i nBPC_mors(%i,%i,%i) nPBC_coul(%i,%i,%i) @VPaul=%li @VLond=%li @VCoul=%li \n", natoms_, Rdamp, alphaMorse, bEvalCoulomb, nBPC_mors.x,nBPC_mors.y,nBPC_mors.z, nPBC_coul.x,nPBC_coul.y,nPBC_coul.z, (long)VPaul, (long)VLond, (long)VCoul );  
        Vec3d* shift_coul=0; int npbc_coul=0;
        Vec3d* shift_mors=0;
        if(bEvalCoulomb) npbc_coul = makePBCshifts_( nPBC_coul, lvec, shift_coul ); // Coulomb converge much more slowly with nPBC
        int npbc_mors = makePBCshifts_( nBPC_mors, lvec, shift_mors );
        double dz = -grid.dCell.c.z;
        Vec3i ns = grid.n;
        //Vec3i ns = gridN;
        const int ntot = ns.totprod();
        const int nxy  = ns.x*ns.y;
        int i=0;
        long t0=getCPUticks();
        //#pragma omp parallel for shared(i,VPaul,VLond,VCoul)
        for(int i=0; i<ntot; i++){
            //int iz =  i%grid.n.z;
            //int iy = (i/grid.n.z)%grid.n.y;
            int iz =  i/nxy;
            int ixy = i - iz*nxy;
            int iy =  ixy/ns.x;
            int ix =  ixy - iy*ns.x;
            const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
            // if((ix==0)&&(iy==0)){ iDBG=iz; 
            //     printf( "CPU pos(%8.4f,%8.4f,%8.4f) iz=%i ", pos.x,pos.y,pos.z );
            // }else{ iDBG=-1; }
            //printf( "GridFF::evalBsplineRef()[%8i] ixyz(%3i,%3i,%3i) pos(%8.4f,%8.4f,%8.4f)\n", i, ix,iy,iz,   pos.x,pos.y,pos.z );
            Quat4d qp,ql;evalGridFFPoint_Mors( npbc_mors, shift_mors,  natoms_, apos_, REQs_, pos, qp, ql );
            const int ibuff = ix + ns.x*( iy + ns.y * iz );
            VPaul[i] = qp.w; 
            VLond[i] = ql.w;
            if(bEvalCoulomb){           
                Quat4d qe  = evalGridFFPoint_Coul( npbc_coul, shift_coul,  natoms_, apos_, REQs_, pos );
                VCoul[i] = qe.w;   
            }
        }
        double T=getCPUticks()-t0;
        if(bEvalCoulomb){  delete [] shift_coul; }
        delete [] shift_mors;
        printf( "GridFF::evalBsplineRef(ntot=%i,natom=%i,nPBC(%i,%i,%i)(coul=%i,morse=%i) DONE in %g [GTick] %g [kTick/point] \n", ntot, natoms_, nPBC.x,nPBC.y,nPBC.z, npbc_coul, shift_mors, T*1e-9, (T*1e-3)/ntot );
    }

    void makeVPLQHeval( int natoms_, Vec3d * apos_, Quat4d * REQs_ ){
        printf( "GridFF::makeVPLQHeval() \n" );
        gridN.x = grid.n.x+3;
        gridN.y = grid.n.y+3;
        gridN.z = grid.n.z+3;
        FEscale.set( -grid.diCell.xx, -grid.diCell.yy, -grid.diCell.zz, 1.0 );  // NOTE: this will not work for non-orthogonal grids
        _realloc0( VPLQH, gridN.totprod(), Quat4dNAN );
        //_realloc0( V_debug, gridN.totprod(), 0.0/0.0 );
        for ( int iz=0; iz<gridN.z; iz++ ){
            int iz_ = fold_cubic( iz, grid.n.z );
            for ( int iy=0; iy<gridN.y; iy++ ){
                int iy_ = fold_cubic( iy, grid.n.y );
                for ( int ix=0; ix<gridN.x; ix++ ){
                    int ix_ = fold_cubic( ix, grid.n.x );
                    Quat4d qp,ql,qe;
                    const Vec3d pos = grid.pos0 + grid.dCell.c*iz_ + grid.dCell.b*iy_ + grid.dCell.a*ix_;
                    evalGridFFPoint( natoms_, apos_, REQs_, pos, qp, ql, qe );
                    const int i = ix + gridN.x* ( iy + gridN.y * iz );
                    VPLQH[ i ] = Quat4d{ qp.w, ql.w, qe.w, 0.0 };
                    //V_debug[ i ] = FFelec_d[j].w; // debug
                }
            }
        }
        printf( "GridFF::makeVPLQHeval() DONE\n" );
    }

    template<typename T>
    void copyPitch( int n, T* dst, int i0dst, int mdst, const T* src, int i0src, int msrc ){
        printf( "copyPitch() n=%i  @src=%li i0src=%i msrc=%i   @dst=%li i0dst=%i mdst=%i \n", n,   (long)src, i0src, msrc,  (long)dst,  i0dst, mdst );
        for ( int i=0; i<n; i++ ){ 
            //if((i%1000)==0){ printf("copyPitch()[%i/%i]\n", i, n); }
            dst[ i*mdst+i0dst ] = src[ i*msrc+i0src ];
        }
    }

    template<typename T>
    void copyPitchTransp( Vec3i ndst, Vec3i transp, T* dst, int i0dst, int mdst, const T* src, int i0src, int msrc ){
        Vec3i nsrc = Vec3i{ ndst.array[transp.x], ndst.array[transp.y], ndst.array[transp.z] };
        printf( "copyPitchTransp() ndst(%i,%i,%i) nsrc(%i,%i,%i)  @src=%li i0src=%i msrc=%i   @dst=%li i0dst=%i mdst=%i \n", ndst.x,ndst.y,ndst.z, nsrc.x,nsrc.y,nsrc.z,   (long)src, i0src, msrc,  (long)dst,  i0dst, mdst );
        Vec3i i;
        for( i.z=0; i.z<ndst.z; i.z++ ){  for( i.y=0; i.y<ndst.y; i.y++ ){  for( i.x=0; i.x<ndst.x; i.x++ ){
            const Vec3i i_ = Vec3i{ i.array[transp.x], i.array[transp.y], i.array[transp.z] };
            const int idst = i.x  + ndst.x*( i.y  + ndst.y * i.z  );
            const int isrc = i_.x + nsrc.x*( i_.y + nsrc.y * i_.z );
            dst[ idst*mdst+i0dst ] = src[ isrc*msrc+i0src ];
        }}}
    }

    template<typename T>
    void copyPBC( Vec3i nsrc, T* src, int i0src, int msrc, Vec3i ndst, T* dst, int i0dst, int mdst ){
        for ( int iz=0; iz<ndst.z; iz++ ){
            int iz_ = fold_cubic( iz, nsrc.z );
            for ( int iy=0; iy<ndst.y; iy++ ){
                int iy_ = fold_cubic( iy, nsrc.y );
                for ( int ix=0; ix<ndst.x; ix++ ){
                    int ix_ = fold_cubic( ix, nsrc.x );
                    const int i = ix  + ndst.x*( iy  + ndst.y * iz  );
                    const int j = ix_ + nsrc.x*( iy_ + nsrc.y * iz_ );
                    src[ i*mdst+i0dst ] =  dst[ j*msrc+i0src ];
                    //V_debug[ i ] = FFelec_d[j].w; // debug
                }
            }
        }
    }

    void makeVPLQH(){
        printf( "GridFF::makeVPLQH() \n" );
        gridN.x = grid.n.x+3;
        gridN.y = grid.n.y+3;
        gridN.z = grid.n.z+3;
        FEscale.set( -grid.diCell.xx, -grid.diCell.yy, -grid.diCell.zz, 1.0 );  // NOTE: this will not work for non-orthogonal grids
        _realloc0( VPLQH, gridN.totprod(), Quat4dNAN );
        //_realloc0( V_debug, gridN.totprod(), 0.0/0.0 );
        for ( int iz=0; iz<gridN.z; iz++ ){
            int iz_ = fold_cubic( iz, grid.n.z );
            for ( int iy=0; iy<gridN.y; iy++ ){
                int iy_ = fold_cubic( iy, grid.n.y );
                for ( int ix=0; ix<gridN.x; ix++ ){
                    int ix_ = fold_cubic( ix, grid.n.x );
                    const int i = ix  + gridN.x* ( iy  + gridN.y  * iz  );
                    const int j = ix_ + grid.n.x*( iy_ + grid.n.y * iz_ );
                    //if( (iz==0)&&(iy==0) ){  printf( "GridFF::makeVPLQH() ix(%i) -> ix(%i)\n", ix, ix_ ); }
                    VPLQH[ i ] = Quat4d{ FFPaul_d[j].w, FFLond_d[j].w, FFelec_d[j].w, 0.0 };
                    //V_debug[ i ] = FFelec_d[j].w; // debug
                }
            }
        }
        printf( "GridFF::makeVPLQH() DONE\n" );
    }

    void FitBsplines( double Ftol=1e-8, int nmaxiter=1000, double dt=0.1 ){
        printf( "GridFF::FitBsplines() \n" );
        int n = gridN.totprod();
        double* Vtemp = new double[ n ];
        double* Ws    = new double[ n ];  for(int i=0; i<n; i++){ Ws[i] = 1.0; }
        copyPBC<double>( gridN, Vtemp,0,1, grid.n, (double*)FFPaul_d, 3,4 ); Bspline::fit3D( gridN, Bspline_Pauli,   Vtemp, Ws, Ftol, nmaxiter, dt );
        copyPBC<double>( gridN, Vtemp,0,1, grid.n, (double*)FFLond_d, 3,4 ); Bspline::fit3D( gridN, Bspline_London,  Vtemp, Ws, Ftol, nmaxiter, dt );
        copyPBC<double>( gridN, Vtemp,0,1, grid.n, (double*)FFelec_d, 3,4 ); Bspline::fit3D( gridN, Bspline_Coulomb, Vtemp, Ws, Ftol, nmaxiter, dt );
        delete [] Vtemp;
        printf( "GridFF::FitBsplines() DONE\n" );
    }

    [[deprecated]]
    void makeGridFF_Bspline_HH_d( double Ftol=1e-8, int nmaxiter=1000, double dt=0.3 ){
        int n = gridN.totprod();
        printf( "GridFF::makeGridFF_Bspline_d() n=%i nmaxiter=%i Ftol=%g dt=%g \n", n );
        _realloc( Bspline_Pauli,    n );
        _realloc( Bspline_London,   n );
        _realloc( Bspline_Coulomb,  n );
        double* Vtemp = new double[ n ];
        //double* Ws    = new double[ n ];  for(int i=0; i<n; i++){ Ws[i] = 1.0; }
        printf( "GridFF::makeGridFF_Bspline_d() START Fitting @Vtemp=%li  @HHermite_d=%li \n", (long)Vtemp, (long)HHermite_d );
        
        //Vec3i nsrc{gridN.z,gridN.y,gridN.x};
        Vec3i transp{2,1,0};
        Vec3i nsrc=gridN; nsrc.swap( transp );
        printf( "GridFF::makeGridFF_Bspline_d() nsrc(%i,%i,%i) transp(%i,%i,%i) \n", nsrc.x, nsrc.y, nsrc.z, transp.x, transp.y, transp.z );
        
        copyPitchTransp<double>( gridN, transp, Vtemp,0,1, (double*)HHermite_d, 0,6 ); Bspline::fit3D_omp( gridN, Bspline_Pauli,   Vtemp, 0, Ftol, nmaxiter, dt, true ); // printf( "GridFF::makeGridFF_Bspline_d() Fit(Bspline_Pauli)   DONE \n" ); //gBS.saveXSF( "debug_VtempPaul.xsf", Vtemp,  1,0  ); gBS.saveXSF( "debug_BsplinePaul.xsf", Bspline_Pauli,   1,0 );
        copyPitchTransp<double>( gridN, transp, Vtemp,0,1, (double*)HHermite_d, 2,6 ); Bspline::fit3D_omp( gridN, Bspline_London,  Vtemp, 0, Ftol, nmaxiter, dt, true ); // printf( "GridFF::makeGridFF_Bspline_d() Fit(Bspline_London)  DONE \n" ); //gBS.saveXSF( "debug_VtempLond.xsf", Vtemp,  1,0  ); gBS.saveXSF( "debug_BsplineLond.xsf", Bspline_London,  1,0 );
        copyPitchTransp<double>( gridN, transp, Vtemp,0,1, (double*)HHermite_d, 4,6 ); Bspline::fit3D_omp( gridN, Bspline_Coulomb, Vtemp, 0, Ftol, nmaxiter, dt, true ); // printf( "GridFF::makeGridFF_Bspline_d() Fit(Bspline_Coulomb) DONE \n" ); //gBS.saveXSF( "debug_VtempCoul.xsf", Vtemp,  1,0  ); gBS.saveXSF( "debug_BsplineCoul.xsf", Bspline_Coulomb, 1,0 );
        printf( "GridFF::makeGridFF_Bspline_d() FINISHED Fitting \n" );
        delete [] Vtemp;
        //delete [] Ws;
        printf( "GridFF::makeGridFF_Bspline_d() DONE\n" );
    }



#ifdef WITH_FFTW

    void makeCoulombEwald(  int natoms_, Vec3d * apos_, Quat4d * REQs_, double* VCoul=0 ){
        std::vector<double> qs(natoms_); 
        std::vector<Vec3d> apos(natoms_);
        for( int i=0; i<natoms_; i++ ){
             qs[i] = REQs_[i].z; 
             apos[i] = apos_[i]+shift0;
        }
        ewald->copy( grid );
        ewald->enlarge( Vec3d{0.0,0.0,20.0} ); // add 20A of Vaccuum to simulated slabe (not periodic in z direction)
        ewald->potential_of_atoms( grid.n.z, VCoul, natoms_, apos_, qs.data() );
    }

#else

    void makeCoulombEwald(  int natoms_, Vec3d * apos_, Quat4d * REQs_, double* VCoul=0 ){
        printf( "GridFF::makeCoulombEwald() ERROR: FFTW is not compiled in this version of GridFF\n" );
        exit(0);
    }

#endif 




    void tryLoadGridFF_potentials( int natoms_, Vec3d * apos_, Quat4d * REQs_, int ntot, double*& VPaul, double*& VLond, double*& VCoul, bool bSaveNPY=false, bool bSaveXSF=false, bool bForceRecalc=false ){
        printf( "GridFF::tryLoadGridFF_potentials() natoms_=%i bUseEwald=%i bForceRecalc=%i bSaveNPY=%i bSaveXSF=%i \n", natoms_, bUseEwald, bForceRecalc, bSaveNPY, bSaveXSF );
        //bool bPBC  = false;
        //Vec3i ns = gridN;
        //Vec3i ns = grid.n;
        //int n = ns.totprod();
        // _allocIfNull( VPaul, ntot );
        // _allocIfNull( VLond, ntot );
        // _allocIfNull( VCoul, ntot );
        bool bPBC    = true;
        GridShape gBS; gBS.copy( grid ); // gBS.n = gridN;
        int nbyte = ntot*sizeof(double);
        //const char* fnames[3] = { "VPaul.bin", "VLond.bin", "VCoul.bin" };
        const char* fnames[3] = { "VPaul.npy", "VLond.npy", "VCoul.npy" };
        if( (!bForceRecalc) && ( checkAllFilesExist( 3, fnames, false ) ) ){
            printf( "found old Bspline referece Energy files %s %s %s\n", fnames[0], fnames[1], fnames[2] );
            // loadBin( fnames[0], nbyte, (char*)VPaul );
            // loadBin( fnames[1], nbyte, (char*)VLond );
            // loadBin( fnames[2], nbyte, (char*)VCoul );
            NumpyFile fPaul( fnames[0] ); VPaul=(double*)fPaul.data;
            NumpyFile fLond( fnames[1] ); VLond=(double*)fLond.data;
            NumpyFile fCoul( fnames[2] ); VCoul=(double*)fCoul.data;
        }else{
            _allocIfNull( VPaul, ntot );
            _allocIfNull( VLond, ntot );
            _allocIfNull( VCoul, ntot );
            if(bUseEwald){
                evalBsplineRef( natoms_, apos_, REQs_, VPaul, VLond, 0 );
                // std::vector<double> qs(natoms_); for( int i=0; i<natoms_; i++ ){ qs[i] = REQs_[i].z; }
                // //std::vector<Vec3d>  ps(natoms_); for( int i=0; i<natoms_; i++ ){  Vec3d p=apos_[i];  p.z-=grid.dCell.c.z; ps[i]=p; } // NOTE: For some unknown reason we have to shift atoms by 1 pixel in z-direction (perhaps something about grid alignement)
                // ewald->copy( grid );
                // ewald->enlarge( Vec3d{0.0,0.0,20.0} ); // add 20A of Vaccuum to simulated slabe (not periodic in z direction)
                // //ewald->potential_of_atoms( grid.n.z, VCoul, natoms_, ps.data(), qs.data() );
                // ewald->potential_of_atoms( grid.n.z, VCoul, natoms_, apos_, qs.data() );
                makeCoulombEwald( natoms_, apos_, REQs_, VCoul );
            }else{
                evalBsplineRef( natoms_, apos_, REQs_, VPaul, VLond, VCoul );
            }
            Vec3i sh{grid.n.z,grid.n.y,grid.n.x};
            save_npy( fnames[0], 3, (int*)&sh, (char*)VPaul );
            save_npy( fnames[1], 3, (int*)&sh, (char*)VLond );
            save_npy( fnames[2], 3, (int*)&sh, (char*)VCoul );
            // saveBin( fnames[0], nbyte, (char*)VPaul );
            // saveBin( fnames[1], nbyte, (char*)VLond );
            // saveBin( fnames[2], nbyte, (char*)VCoul );
            if(bSaveNPY){
                save_npy( "debug_VPaul.npy", 3, (int*)&sh, (char*)VPaul );
                save_npy( "debug_VLond.npy", 3, (int*)&sh, (char*)VLond );
                save_npy( "debug_VCoul.npy", 3, (int*)&sh, (char*)VCoul );
            }
            if(bSaveXSF){
                gBS.saveXSF( "debug_VPaul.xsf",   VPaul, 1,0 );
                gBS.saveXSF( "debug_VLond.xsf",   VLond, 1,0 );
                gBS.saveXSF( "debug_VCoul.xsf",   VCoul, 1,0 );
            }
        }
    }

    bool makeGridFF_Bspline_d( int natoms_, Vec3d * apos_, Quat4d * REQs_, bool bSaveNPY=false, bool bSaveXSF=false, bool bFit=true, bool bRefine=false ){
        printf( "GridFF::makeGridFF_Bspline_d() bUseEwald=%i bSaveNPY=%i bSaveXSF=%i bFit=%i \n", bUseEwald, bSaveNPY, bSaveXSF, bFit );
        //Vec3i ns = gridN;
        Vec3i ns = grid.n;
        int ntot = ns.totprod();
        //int nbyte = ntot*sizeof(double);
        double *VPaul=0, *VLond=0, *VCoul=0;
        tryLoadGridFF_potentials( natoms_, apos_, REQs_, ntot, VPaul, VLond, VCoul, bSaveNPY, bSaveXSF );
        if(!bFit){
            printf( "GridFF::makeGridFF_Bspline_d() bFit=false => exit() \n" );
            _dealloc( VPaul );_dealloc( VLond );_dealloc( VCoul );
            //exit(0);
            return false;
        }
        //save1Dslice( "VPaul_x.log", n, Vec3i i0, Vec3i di, Vec3i ns, int m, double* data, const char* fmt="%g "  );
        // golbal_array_dict.insert( { "Bspline_Pauli",   NDArray{ Bspline_Pauli, Quat4i{gridN.z,gridN.y,gridN.x,-1}} }  );
        //double* Ws    = new double[ n ];  for(int i=0; i<n; i++){ Ws[i] = 1.0; }
        _realloc( Bspline_Pauli,   ntot );
        _realloc( Bspline_London,  ntot );
        _realloc( Bspline_Coulomb, ntot );
        Quat4d testREQ{  1.8, sqrt(0.0091063), 0.6*0.3, 0.0 }; // ~ Oxygen
        Quat4f testPLQ = REQ2PLQ( testREQ, alphaMorse );
        
        double vmax_Paul = VecN::absmax( ntot, VPaul );
        double vmax_Lond = VecN::absmax( ntot, VLond ); 
        double vmax_Coul = VecN::absmax( ntot, VCoul );
        printf( "GridFF::makeGridFF_Bspline_d() testPLQ( Paul=%g, Lond=%g, Coul=%g )  vmaxs( P=%g, L=%g, C=%g )\n", testPLQ.x, testPLQ.y, testPLQ.z, vmax_Paul, vmax_Lond, vmax_Coul );
        
        bool   bPBC = true;
        bool   bInitGE = true;
        double Ftol  = 1e-6;
        double dt    = 0.001;
        int nmaxiter = 50000;
        
        //double Ftol_pre     = 1e-8;
        //double dt_pre       = 0.1;
        //int    nmaxiter_pre = 2000; 

        double Ftol_pre     = 1e-13;
        double dt_pre       = 0.3;
        int    nmaxiter_pre = 5000; 
        double damp         = 0.0; // NOTE / ToDo : It is a bit strange than in GridFF.cl we use damp=0.3 and it is faster, but here the fit is fastes with damp=0.0

        printf("GridFF::makeGridFF_Bspline_d() ========= FITTING start======== \n");
        int niter=0;
        printf( "GridFF::makeGridFF_Bspline_d().fit(Bspline_Pauli  ) : "); 
        niter=           Bspline::fit3D_omp( ns, Bspline_Pauli,   VPaul, 0, Ftol_pre, nmaxiter_pre, dt_pre, damp, bPBC, bInitGE );                   // raw pre-fit 
        //if(bRefine)niter=Bspline::fit3D_omp( ns, Bspline_Pauli,   VPaul, 0, Ftol    , nmaxiter    , dt    , bPBC,  false, vmax_Paul*1e-3 );    // refine fit
        if(niter>=nmaxiter-1){printf("fit3D_omp not cnverged => exit() \n"); exit(0);}; //printf( "GridFF::makeGridFF_Bspline_d() Fit(Bspline_Pauli)   DONE \n" );
        
        printf( "GridFF::makeGridFF_Bspline_d().fit(Bspline_London ) : "); 
        niter=           Bspline::fit3D_omp( ns, Bspline_London,  VLond, 0, Ftol_pre, nmaxiter_pre, dt_pre, damp, bPBC, bInitGE );                   // raw pre-fit 
        //if(bRefine)niter=Bspline::fit3D_omp( ns, Bspline_London,  VLond, 0, Ftol    , nmaxiter,     dt,     bPBC,   false, vmax_Lond*1e-3 );   // refine fit
        if(niter>=nmaxiter-1){printf("fit3D_omp not cnverged => exit() \n"); exit(0);}; //printf( "GridFF::makeGridFF_Bspline_d() Fit(Bspline_London)  DONE \n" );
        
        printf( "GridFF::makeGridFF_Bspline_d().fit(Bspline_Coulomb) : "); 
        niter=           Bspline::fit3D_omp( ns, Bspline_Coulomb, VCoul, 0, Ftol_pre, nmaxiter_pre, dt_pre, damp, bPBC, bInitGE );                  // raw pre-fit 
        //if(bRefine)niter=Bspline::fit3D_omp( ns, Bspline_Coulomb, VCoul, 0, Ftol,     nmaxiter,     dt,     bPBC,   false, vmax_Coul*1e-3 );  // refine fit
        if(niter>=nmaxiter-1){printf("fit3D_omp not cnverged => exit() \n"); exit(0);}; //printf( "GridFF::makeGridFF_Bspline_d() Fit(Bspline_Coulomb) DONE \n" );

        printf("GridFF::makeGridFF_Bspline_d() ========= FITTING finish ======== \n");

        GridShape gBS; gBS.copy( grid ); // gBS.n = gridN;
        if(bPBC){
            if(bSaveNPY){
                Vec3i sh{grid.n.z,grid.n.y,grid.n.x};
                save_npy( "debug_BsplinePaul_pbc.npy", 3, (int*)&sh, (char*)Bspline_Pauli   );
                save_npy( "debug_BsplineLond_pbc.npy", 3, (int*)&sh, (char*)Bspline_London  );
                save_npy( "debug_BsplineCoul_pbc.npy", 3, (int*)&sh, (char*)Bspline_Coulomb );
            }
            if(bSaveXSF){
                gBS.saveXSF( "debug_BsplinePaul_pbc.xsf",   Bspline_Pauli,   1,0 );
                gBS.saveXSF( "debug_BsplineLond_pbc.xsf",   Bspline_London,  1,0 );
                gBS.saveXSF( "debug_BsplineCoul_pbc.xsf",   Bspline_Coulomb, 1,0 );
            }
        }else{
            if(bSaveNPY){
                Vec3i sh{grid.n.z,grid.n.y,grid.n.x};
                save_npy( "debug_BsplinePaul_pbc.npy", 3, (int*)&sh, (char*)Bspline_Pauli   );
                save_npy( "debug_BsplineLond_pbc.npy", 3, (int*)&sh, (char*)Bspline_London  );
                save_npy( "debug_BsplineCoul_pbc.npy", 3, (int*)&sh, (char*)Bspline_Coulomb );
            }
            if(bSaveXSF){
                gBS.saveXSF( "debug_BsplinePaul_nopbc.xsf", Bspline_Pauli,   1,0 );
                gBS.saveXSF( "debug_BsplineLond_nopbc.xsf", Bspline_London,  1,0 );
                gBS.saveXSF( "debug_BsplineCoul_nopbc.xsf", Bspline_Coulomb, 1,0 );
            }
        }
        //printf( "GridFF::makeGridFF_Bspline_d() FINISHED Fitting \n" );
        _dealloc( VPaul );
        _dealloc( VLond );
        _dealloc( VCoul );
        //delete [] Ws;
        //printf( "GridFF::makeGridFF_Bspline_d() DONE\n" );
        return true;
    }

    void pack_Bspline_d( ){
        printf( "GridFF::pack_Bspline_d() grid.n(%i,%i,%i) @Bspline_Pauli=%li @Bspline_London=%li @Bspline_Coulomb=%li \n", grid.n.x,grid.n.y,grid.n.z, (long)Bspline_Pauli, (long)Bspline_London, (long)Bspline_Coulomb );
        int ntot = grid.n.totprod();
        _realloc( Bspline_PLQ, ntot );
        for(int ix=0; ix<grid.n.x; ix++){
            for(int iy=0; iy<grid.n.y; iy++){
                for(int iz=0; iz<grid.n.z; iz++){
                    int j = ix + grid.n.x*( iy + iz*grid.n.y );
                    int i = iz + grid.n.z*( iy + ix*grid.n.y );
                    Bspline_PLQ[i] = Vec3d{  Bspline_Pauli[j], Bspline_London[j], Bspline_Coulomb[j] };  
                }
            }
        }
        printf( "GridFF::pack_Bspline_d() DONE @Bspline_PLQ=%li\n", (long)Bspline_PLQ );
    }

    int Bspline_to_f4( bool bAlloc ){
        int nxyz = grid.n.totprod();
        if(bAlloc) _realloc( Bspline_PLQf, nxyz );
        for(int i=0; i<nxyz; i++){
            Vec3d v = Bspline_PLQ[i];
            Bspline_PLQf[i] = Quat4f{ (float)v.x, (float)v.y, (float)v.z, 0 };
        }
        return nxyz;
    }

    double evalMorsePBC(  Vec3d pi, Quat4d REQi, Vec3d& fi, int natoms, Vec3d * apos, Quat4d * REQs ){
        //printf( "GridFF::evalMorsePBC() debug fi(%g,%g,%g) REQi(%g,%g,%g)\n",  fi.x,fi.y,fi.z, REQi.x,REQi.y,REQi.z,REQi.w  );
        const double R2damp=Rdamp*Rdamp;    
        const double K =-alphaMorse;
        double       E = 0;
        Vec3d f = Vec3dZero;
        //printf("GridFF::evalMorsePBC() npbc=%i natoms=%i bSymetrized=%i \n", npbc, natoms, bSymetrized );
        if(!bSymetrized){ printf("ERROR  GridFF::evalMorsePBC() not symmetrized, call  GridFF::setAtomsSymetrized() first => exit()\n"); exit(0); }
        if( (shifts==0) || (npbc==0) ){ printf("ERROR in GridFF::evalMorsePBC() pbc_shift not intitalized !\n"); };     
        for(int j=0; j<natoms; j ++ ){    // atom-atom
            Vec3d fij = Vec3dZero;
            Vec3d dp0 = pi - apos[j] + shift0;
            Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
            //printf( "GridFF::evalMorsePBC() j %i/%i \n", j,natoms );
            for(int ipbc=0; ipbc<npbc; ipbc++ ){
                //printf( "GridFF::evalMorsePBC() j %i/%i ipbc %i/%i \n", j,natoms, ipbc,npbc );
                const Vec3d  dp = dp0 + shifts[ipbc];
                Vec3d fij=Vec3dZero;
                E += addAtomicForceMorseQ( dp, fij, REQij.x, REQij.y, REQij.z, K, R2damp );
                //E  += exp(-dp.norm2()/0.16 );
                f.sub( fij );
            }
        }
        fi.add( f );
        //printf( "CPU[iG=0,iS=0] fe(%10.6f,%10.6f,%10.6f)\n", fi.x,fi.y,fi.z );
        return E;
    }
    double evalMorsePBC_sym     (         Vec3d  pi, Quat4d  REQi, Vec3d& fi     ){ return evalMorsePBC( pi, REQi, fi, natoms, &apos[0], &REQs[0] ); }
    // double evalMorsePBC_sym     (         Vec3d  pi, Quat4d  REQi, Vec3d& fi     ){ return evalMorsePBC( pi, REQi, fi, apos_.size(), &apos_[0], &REQs_[0] ); }
    double evalMorsePBCatoms_sym( int na, Vec3d* ps, Quat4d* REQs, Vec3d* forces ){
        double E = 0;
        for(int ia=0; ia<na; ia++){ evalMorsePBC_sym( ps[ia], REQs[ia], forces[ia] ); };
        return E;
    }


    __attribute__((hot))  
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
            if(FFPaul) FFPaul[ibuff]=(Quat4f)qp;
            if(FFLond) FFLond[ibuff]=(Quat4f)ql;
            if(FFelec) FFelec[ibuff]=(Quat4f)qe;
        });
    }

    // void evalGridFFs( Vec3i nPBC_=Vec3i{-1,-1,-1} ){
    //     if(nPBC_.x>=0) nPBC=nPBC_;
    //     //evalGridFFexp( natoms, apos, REQs, alpha*2,  1, FFPaul  );
    //     //evalGridFFexp( natoms, apos, REQs, alpha  , 1, FFLond );  // -2.0 coef is in  REQ2PLQ
    //     //evalGridFFel ( natoms, apos, REQs,              FFelec   );
    //     //evalGridFFs( natoms, apos, REQs );
    //     if(iDebugEvalR>0){ evalGridR  ( natoms, apos, REQs );}
    //     else             { evalGridFFs( natoms, apos, REQs ); }
    // }

    void evalCombindGridFF( Quat4d REQ, Quat4f * FF ){
        Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );
        interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
            Quat4f f = Quat4fZero;
            if(FFPaul) f.add_mul( FFPaul [ibuff], PLQ.x );
            if(FFLond) f.add_mul( FFLond[ibuff], PLQ.y );
            if(FFelec) f.add_mul( FFelec  [ibuff], PLQ.z );
            FF[ibuff] =  f;
        });
    }

    // void evalCombindGridFF_CheckInterp( Quat4d REQ, Quat4f * FF ){
    //     Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );
    //     interateGrid3D( grid, [=](int ibuff, Vec3d p)->void{
    //         Quat4f f = Quat4fZero;
    //         addForce( p+Vec3d{0.1,0.1,0.1}, (Quat4f)PLQ, f );
    //         FF[ibuff] = f;
    //     });
    // }

void setAtomsSymetrized( int n, int* atypes, Vec3d* apos, Quat4d* REQs, double d=0.1 ){
    //printf( "GridFF::setAtomsSymetrized() \n" );
    double cmax =-0.5+d;
    double cmin = 0.5-d;
    apos_  .clear();
    REQs_  .clear();
    atypes_.clear();
    Mat3d M; grid.cell.invert_T_to( M );
    const Vec3d& a = grid.cell.a;
    const Vec3d& b = grid.cell.b;
    //printf( "GridFF::setAtomsSymetrized() cmin %g cmax %g \n", cmin, cmax );
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
    //printf( "GridFF::setAtomsSymetrized() END na_new=%i na_old=%i \n", atypes_.size(), n );
    //bindSystem( atypes_.size(), &atypes_[0], &apos_[0], &REQs_[0] );
    bSymetrized = true;
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
    Vec3d dp=p1-p0; dp.mul(1./n);
    //Quat4f PLQf = REQ2PLQ(   REQ, alphaMorse );   //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    Quat4d PLQd = REQ2PLQ_d( REQ, alphaMorse );
    if(bPrint){ printf("GridFF::getEFprofile(n=%i,p2{%6.3f,%6.3f,%6.3f},p1{,%6.3f,%6.3f,%6.3f}) REQ(%g,%g,%g)  PLQ(%g,%g,%g) shift0(%g,%g,%g) grid.pos0(%g,%g,%g) \n", n, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z, REQ.x,REQ.y,REQ.z,  PLQd.x,PLQd.y,PLQd.z, shift0.x,shift0.y,shift0.z, grid.pos0.x,grid.pos0.y,grid.pos0.z  ); };
    for(int i=0; i<n; i++){
        Vec3d  p  = p0 + dp*i;
        //Quat4f fe = getForce( p, PLQ, true );
        //if(fes)fes[i]=(Quat4d)fe;
        Quat4d fed=Quat4dZero;
        fed.e = addAtom( p, PLQd, fed.f );
        if(bPrint){ printf( "%4i %6.3f %6.3f %6.3f    %20.10e %20.10e %20.10e    %20.10e\n", i, p.x, p.y, p.z, fed.x,fed.y,fed.z, fed.w  ); };
    }
}

void getEFprofileToFile( const char* fname, int n, Vec3d p0, Vec3d p1, Quat4d REQ ){
    freopen(fname, "w", stdout);
    getEFprofile( n, p0, p1, REQ, 0, true );
    fclose(stdout);
    freopen("/dev/tty", "a", stdout);
}

double checkEFProfileVsNBFF( int n, Vec3d p0, Vec3d p1, const Quat4d& REQ, double tol=1e-2, bool bExit=false, bool bPrint=false, bool bWarn=true, const char* logfiflename="checkEFProfileVsNBFF.log" ){
    if(bPrint){ printf("GridFF::checkEFProfileVsNBFF(np=%i,natoms=%i,npbc=%i,p2{%6.3f,%6.3f,%6.3f},p1{,%6.3f,%6.3f,%6.3f}) \n", n, natoms,npbc, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z ); };
    FILE * logf=0;
    if(logfiflename){ 
        logf = fopen(logfiflename,"w");
        fprintf( logf, "GridFF::checkEFProfileVsNBFF(np=%i,natoms=%i,npbc=%i,p2{%6.3f,%6.3f,%6.3f},p1{,%6.3f,%6.3f,%6.3f}REQ{%g,%g,%g,%g}) \n", n, natoms,npbc, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z, REQ.x,REQ.y,REQ.z,REQ.w  );
        fprintf( logf, "i   x y z     E  Eref     fx fx_ref      fy fy_ref     fz  fz_ref\n");
    }
    if(bPrint){     printf("i   x y z     E  Eref     fx fx_ref      fy fy_ref     fz  fz_ref\n"); }
    double tol2=tol*tol;
    Vec3d dp=p1-p0; dp.mul(1./n);
    Quat4f PLQ = REQ2PLQ( REQ, alphaMorse );   //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    //bool err = false;
    bool bErr=false;
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
            bErr=true;
            if(bWarn)printf("WARRNING[%i/%i] dE=%g |dF|=%g p(%6.3f,%6.3f,%6.3f) GridFF(%g,%g,%g|%g)  NBFF(%g.%g,%g|%g)\n", i, n, dE, df.norm(), p.x,p.y,p.z,   fe.x,fe.y,fe.z,fe.w,   fref.x,fref.y,fref.z,Eref  );
            if(bExit){ printf("ERROR in GridFF::checkEFProfileVsNBFF() - GridFF force does not match NBFF reference at test point %i MaxRelativeError=%g => Exit()\n", i, sqrt(err2Max) ); exit(0); }
        } 
        if(bPrint){ printf(       "%i    %6.3f %6.3f %6.3f    %g %g   %g %g    %g %g    %g %g\n", i, p.x, p.y, p.z,    fe.e, Eref,    fe.x,fref.x,    fe.y,fref.y,    fe.z,fref.z ); }
        if(logf  ){fprintf( logf, "%3i    %6.3f %6.3f %6.3f    %14.6f %14.6f    %14.6f %14.6f    %14.6f %14.6f %14.6f    %14.6f\n", i, p.x, p.y, p.z,    fe.e, Eref,    fe.x,fref.x,    fe.y,fref.y,    fe.z,fref.z ); }
    }
    if(logf){ fclose(logf); }
    if(bWarn && bErr ){
        //printf("WARRNING GridFF MaxRelativeError=%g at point[%i]\n",  sqrt(err2Max), imax );
        Vec3d fref;
        Vec3d  p    = p0 + dp*imax;
        Quat4f fe   = getForce          ( p, PLQ, true );
        float  Eref = getMorseQH_PBC_omp( p, REQ, fref );
        float  dE   = fe.e-Eref;
        Vec3f  df   = fe.f - (Vec3f)fref;
        double e2e = dE*dE     /(Eref*Eref + fe.e*fe.e + 1);          err2Max=fmax(err2Max,e2e);
        double e2f = df.norm2()/( fe.f.norm2() + fref.norm2() + 1 );  err2Max=fmax(err2Max,e2f);
        printf("WARRNING GridFF MaxError=%g at[%i/%i] dE=%g |dF|=%g p(%6.3f,%6.3f,%6.3f) GridFF(%g.%g,%g|%g)  NBFF(%g.%g,%g|%g)\n",  sqrt(err2Max), imax, n, dE, df.norm(), p.x,p.y,p.z, fe.x,fe.y,fe.z,fe.w,   fref.x,fref.y,fref.z,Eref  );
    }
    return sqrt(err2Max);
}

bool checkZProfilesOverAtom( int ia, int n, double zmin, double zmax, const Quat4d& REQ, double tol=1e-2, bool bExit=true, bool bPrint=true ){
    Vec3d p0=apos[ia]; p0.z=zmin;
    Vec3d p1=p0;       p1.z=zmax;
    double err =  checkEFProfileVsNBFF( n, p0, p1, REQ, tol, false,false, true );
    if(err>tol){    
        if(bPrint)checkEFProfileVsNBFF( n, p0, p1, REQ, tol,  false, true, false );
        if(bExit){ printf("ERROR in GridFF::checkZProfilesOverAtom(%i) - GridFF force does not match NBFF reference, MaxRelativeError=%g => Exit()\n", ia, err ); exit(0); }
        return true;
    }
    //checkEFProfileVsNBFF( n, p0, p1, REQ, tol,  false, true, false );
    return false;
}

bool evalCheck( int imin=0, int imax=1, bool bExit=true, bool bPrint=true, double tol=1e-2, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.1, 0.}, double dz=0.05 ){
    REQ=Quat4d{ 1.487, 0.02609214441, -0.1, 0.};
    printf( "GridFF::evalCheck() natoms=%i npbc=%i apos=%li REQs=%li shifts=%li \n", natoms, npbc, apos, REQs, shifts );
    _checkNull(shifts)
    _checkNull(REQs)
    _checkNull(apos)
    bool err = false;
    double zmin=grid.pos0.z+shift0.z+1.0;
    double zmax=grid.pos0.z+grid.cell.c.z*0.5;
    zmax=fmin(zmax,10);
    int nz = ((int)((zmax-zmin)/dz)) + 1;
    for( int ia=imin; ia<imax; ia++ ){
        err |= checkZProfilesOverAtom( ia, nz, zmin, zmax, REQ, tol, bExit, bPrint );
    }
    return err;
}

void log_z(const char* fname, int ix=0, int iy=0){
    FILE* logf=fopen( fname, "w");
    if(logf==0){ printf("ERROR in GridFF::makeGridFF_omp() cannot open logfile(%s) => Exit()\n", fname ); exit(0); }
    fprintf( logf, "#i   z  Ep_Paul Fz_Paul   Ep_Lond Fz_Lond  E_Coul Fz_Coul \n" );
    for ( int iz=0; iz<grid.n.z; iz++ ){
        const Vec3d pos = grid.pos0 + grid.dCell.c*iz + grid.dCell.b*iy + grid.dCell.a*ix;
        const int ibuff = ix + grid.n.x*( iy + grid.n.y * iz );
        Quat4f qp = FFPaul[ibuff];
        Quat4f ql = FFLond[ibuff];
        Quat4f qe = FFelec[ibuff];
        if(logf && (ix==0) && (iy==0) ){ fprintf( logf,  "%3i %8.3f    %14.6f %14.6f    %14.6f %14.6f    %14.6f %14.6f\n", iz, pos.z,  qp.w,qp.z,   ql.w,ql.z,  qe.w,qe.z ); }
    }
    fclose(logf);
}

void checkSum( bool bDouble ){
    int n = grid.getNtot();
    if( bDouble ){
        Quat4d paul = sum( n, FFPaul_d, Quat4dZero );
        Quat4d lond = sum( n, FFLond_d, Quat4dZero );
        Quat4d elec = sum( n, FFelec_d, Quat4dZero );
        printf( "GridFF::checkSum() bDouble=1 Pauli %g %g %g %g | London %g %g %g %g | Coulomb %g %g %g %g \n", paul.x,paul.y,paul.z,paul.w, lond.x,lond.y,lond.z,lond.w, elec.x,elec.y,elec.z,elec.w );
    }else{
        Quat4f paul = sum( n, FFPaul, Quat4fZero );
        Quat4f lond = sum( n, FFLond, Quat4fZero );
        Quat4f elec = sum( n, FFelec, Quat4fZero );
        printf( "GridFF::checkSum() bDouble=0 Pauli %g %g %g %g | London %g %g %g %g | Coulomb %g %g %g %g \n", paul.x,paul.y,paul.z,paul.w, lond.x,lond.y,lond.z,lond.w, elec.x,elec.y,elec.z,elec.w );
    }
} 

void initGridFF( const char * name, double z0=NAN, bool bAutoNPBC=true, bool bSymetrize=false, double rAutoPBC=20.0 ){
    printf( "GridFF::initGridFF(%s) bSymetrize=%i bAutoNPBC=%i z0=%g \n", name, bSymetrize, bAutoNPBC, z0 );
    if( isnan(z0) ){  
        z0=findTop();   
        //if(verbosity>0) 
        printf("GridFF::findTop() z0=%g \n", z0);  
    }
    grid.pos0.z=z0;
    if(verbosity>1)grid.printCell();
    bool bGridDouble = (mode == GridFFmod::LinearDouble) || (mode == GridFFmod::HermiteDouble) || (mode == GridFFmod::BsplineDouble); 
    allocateFFs( bGridDouble );
    //gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false, {1,1,0}, bSaveDebugXSFs );
//    nPBC=Vec3i{1,1,0};
    if(bAutoNPBC){ autoNPBC( grid.cell, nPBC, 20.0 ); }   printf( "GridFF::initGridFF() nPBC(%i,%i,%i)\n", nPBC.x, nPBC.y, nPBC.z );
    lvec = grid.cell;     // ToDo: We should unify this
    makePBCshifts( nPBC, lvec );
    if(bSymetrize)setAtomsSymetrized( natoms, atypes, apos, REQs, 0.1 );
    //apos_.size(), &apos_[0], &REQs_[0];
    std::vector<double> qs( apos_.size() ); for(int i=0; i<apos_.size(); i++){ qs[i]=REQs_[i].z; }
    Vec3d qcog;
    Multiplole::project( &qcog, apos_.size(), apos_.data(), qs.data(), 2, Mpol, true );
    printf( "GridFF::initGridFF() nPBC(%i,%i,%i|%i) Mpol qcog(%g,%g,%g) Q=%g dipol(%g,%g,%g) quadrupol:xx,yy,zz(%g,%g,%g)yz,xz,xy(%g,%g,%g)\n", nPBC.x,nPBC.x,nPBC.x,npbc, qcog.x,qcog.y,qcog.z, Mpol[0], Mpol[1],Mpol[2],Mpol[3],  Mpol[4],Mpol[5],Mpol[6], Mpol[7],Mpol[8],Mpol[9] );
    //bSaveDebugXSFs=true;
    gridN=grid.n; gridN.x+=3; gridN.y+=3;
}



 #ifdef IO_utils_h
    bool tryLoad( const char* fname_Coul, const char* fname_Paul, const char* fname_Lond, bool recalcFF=false, bool bDouble=false ){
        //printf( "GridFF::tryLoad() \n" );
        //printf( "GridFF::tryLoad() fname_Pauli >>%s<< fname_London >>%s<< fname_Coulomb >>%s<< \n", fname_Pauli, fname_London, fname_Coulomb );
        //printDir( "../" );
        //printDir( "./" );
        { FILE* f=fopen( fname_Paul,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Paul); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_Lond,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Lond); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        { FILE* f=fopen( fname_Coul,"rb"); if(0==f){ printf("File(%s) Not Found\n", fname_Coul); recalcFF=true; }else{ fclose(f); };} // Test if file exist
        printf( "GridFF::tryLoad() recalcFF %i \n", recalcFF );
        //printf( "fname_Pauli(%s) fname_London(%s) fname_Coulomb(%s) \n", fname_Pauli, fname_London, fname_Coulomb );
        //int nbyte= grid.getNtot()*sizeof(Vec3d);
        int nbyte = grid.getNtot()*sizeof(Quat4f);
        char* cFFPaul = (char*)FFPaul;
        char* cFFLond = (char*)FFLond;
        char* cFFelec = (char*)FFelec;
        if(bDouble){
            printf( "GridFF::tryLoad() bDouble %i grid.n(%i,%i,%i)\n", bDouble, grid.n.x, grid.n.y, grid.n.z );
            nbyte = grid.getNtot()*sizeof(Quat4d);
            cFFPaul = (char*)FFPaul_d;
            cFFLond = (char*)FFLond_d;
            cFFelec = (char*)FFelec_d;
        }
        if( recalcFF ){
            printf( "\nBuilding GridFF for substrate (bDouble=%i) grid.n(%i,%i,%i) ... (please wait... )\n", bDouble, grid.n.x, grid.n.y, grid.n.z );
            if(bDouble){ makeGridFF_omp_d( apos_.size(), &apos_[0], &REQs_[0] ); }
            else       { makeGridFF_omp  ( apos_.size(), &apos_[0], &REQs_[0] ); }
            if(cFFPaul) saveBin( fname_Paul,  nbyte, cFFPaul );
            if(cFFLond) saveBin( fname_Lond,  nbyte, cFFLond );
            if(cFFelec) saveBin( fname_Coul,  nbyte, cFFelec );
        }else{
            if(cFFPaul) loadBin( fname_Paul,  nbyte, cFFPaul );
            if(cFFLond) loadBin( fname_Lond,  nbyte, cFFLond );
            if(cFFelec) loadBin( fname_Coul,  nbyte, cFFelec );
        }
        return recalcFF;
    }

    bool tryLoad_new( bool bSymetrize=true, bool bFit=true, bool bRefine=true, bool bPrint=false ){
        printf( "GridFF::tryLoad_new() mode=%i bSymetrize=%i bFit=%i bRefine=%i \n", (int)mode, bSymetrize, bFit, bRefine );
        // print current directory
        char cwd[128]; getcwd(cwd, 128); printf( "Current  directory: %s\n", cwd );
        //const char* fname_Coul=0; 
        //const char* fname_Paul=0; 
        //const char* fname_Lond=0;
        const char* fnames[4];
        int npoint=0;
        int nbyte=0;
        bool bRecalc = false; 
        long t0 = getCPUticks();

        int     na   = natoms; 
        Vec3d*  aps  = apos;
        Quat4d* reqs = REQs;
        if(bSymetrize){
            na   = apos_.size(); 
            aps  = &apos_[0];
            reqs = &REQs_[0];
        }
        //printf( "GridFF::tryLoad_new() na=%i @apos=%li @REQs=%li \n", na, (long)aps, (long)reqs );

        switch (mode){
            
            case GridFFmod::LinearFloat:{
                fnames[0]="FFPaul_f.bin"; fnames[1]="FFLond_f.bin"; fnames[2]="FFelec_f.bin"; 
                npoint = grid.n.totprod();
                nbyte  = npoint*sizeof(Quat4f);
                if( checkAllFilesExist( 3, fnames, bPrint ) ){
                    loadBin( fnames[0], nbyte, (char*)FFPaul );
                    loadBin( fnames[1], nbyte, (char*)FFLond );    
                    loadBin( fnames[2], nbyte, (char*)FFelec );
                }else{
                    bRecalc=true;
                    makeGridFF_omp  ( na, aps, reqs ); 
                    saveBin( fnames[0], nbyte, (char*)FFPaul );
                    saveBin( fnames[1], nbyte, (char*)FFLond );    
                    saveBin( fnames[2], nbyte, (char*)FFelec );
                }
            } break;
            
            case GridFFmod::LinearDouble:{
                fnames[0]="FFPaul_d.bin"; fnames[1]="FFLond_d.bin"; fnames[2]="FFelec_d.bin"; 
                npoint = grid.n.totprod();
                nbyte = npoint*sizeof(Quat4d);
                if( checkAllFilesExist( 3, fnames, bPrint ) ){
                    loadBin( fnames[0], nbyte, (char*)FFPaul_d );
                    loadBin( fnames[1], nbyte, (char*)FFLond_d );    
                    loadBin( fnames[2], nbyte, (char*)FFelec_d );
                }else{
                    bRecalc=true;
                    makeGridFF_omp_d(  na, aps, reqs ); 
                    saveBin( fnames[0], nbyte, (char*)FFPaul_d );
                    saveBin( fnames[1], nbyte, (char*)FFLond_d );    
                    saveBin( fnames[2], nbyte, (char*)FFelec_d );
                }
            } break;
            
            case GridFFmod::HermiteFloat    :{ printf("ERROR in GridFF::tryLoad_new() GridFFmode::HermiteFloat NOT IMPLEMENTED \n"); exit(0); } break;
            case GridFFmod::HermiteDouble   :{ 
                fnames[0] = "GridFF_HH_d.bin";
                npoint = gridN.totprod();
                nbyte  = npoint*sizeof(double)*6;
                perVoxel=6;
                if( checkAllFilesExist( 1, fnames, bPrint ) ){
                    printf("Load HermiteDouble npoint=%i nbyte=%i \n", npoint, nbyte );
                    _realloc( HHermite_d, npoint*6 );
                    loadBin( fnames[0], nbyte, (char*)HHermite_d );
                }else{
                    printf("Recalc & Load HermiteDouble npoint=%i nbyte=%i \n", npoint, nbyte );
                    bRecalc=true;
                    makeGridFF_Hherm_d(  na, aps, reqs ); 
                    saveBin( fnames[0], nbyte, (char*)HHermite_d );
                }
            } break;

            case GridFFmod::BsplineFloat    :{ printf("ERROR in GridFF::tryLoad_new() GridFFmode::BsplineFloat NOT IMPLEMENTED \n"); exit(0); } break;
            case GridFFmod::BsplineDouble   :{
                //fnames[0]="Bspline_Pauli_d.bin"; fnames[1]="Bspline_London_d.bin"; fnames[2]="Bspline_Coulomb_d.bin";
                fnames[0]="Bspline_PLQd.npy";
                fnames[1]="Bspline_PLQ_d.bin";
                //Vec3i ns = gridN;
                Vec3i ns = grid.n;
                npoint = ns.totprod();
                nbyte  = npoint*sizeof(double);
                bool done=false;
                if(       fileExist( fnames[0] ) ){
                    printf( "GridFF::tryLoad_new() load Bspline_PLQ from %s\n", fnames[0] );
                    //_realloc( Bspline_PLQ, npoint*3 );
                    NumpyFile npy(fnames[0]);
                    npy.print();
                    Bspline_PLQ = (Vec3d*)npy.data;
                    done=true;
                } else if( fileExist( fnames[1] ) ){
                    _realloc( Bspline_PLQ, npoint*3 );
                    loadBin( fnames[1], nbyte*3, (char*)Bspline_PLQ );
                    done=true;
                }else{
                    //printf( "GridFF::tryLoad_new() GridFFmod::BsplineDouble no file exist, but we do not want generate new, takes too long, we lazy \n" ); exit(0);
                    bRecalc=true;
                    bool bSaveNPY=false;
                    bool bSaveXSF=false;
                    //printf( "tryLoad_new() mode=BsplineDouble bSaveNPY=%i bSaveXSF=%i bFit=%i bRefine=%i \n", bSaveNPY, bSaveXSF, bFit, bRefine );
                    done = makeGridFF_Bspline_d(  na, aps, reqs, bSaveNPY, bSaveXSF, bFit, bRefine);
                    if(done){
                        pack_Bspline_d();
                        //saveBin( fnames[1], nbyte*3, (char*)Bspline_PLQ );
                        Quat4i sh{grid.n.x,grid.n.y,grid.n.z,3};
                        save_npy( fnames[0], 4,(int*)&sh, (char*)Bspline_PLQ );
                    }
                }

                { // Debug - check content of Bspline_PLQ
                    //printf( "GridFF::tryLoad_new() \n", grid.n.x,grid.n.y,grid.n.z,  grid.dCell.xx, grid.dCell.yy, grid.dCell.zz,  grid.cell.xx, grid.cell.yy, grid.cell.zz, grid.pos0, grid.pos0, grid.pos0 )
                    grid.printCell();
                    Vec3d vmin=Vec3dmax,vmax=Vec3dmin;
                    for(int i=0; i<npoint; i++){ Bspline_PLQ[i].update_bounds(vmin,vmax); }
                    printf( "GridFF::tryLoad_new() Bspline_PLQ Paul(%g,%g) Lond(%g,%g) Coul(%g,%g) \n", vmin.x,vmax.x, vmin.y,vmax.y, vmin.z,vmax.z  );
                    // int ix=18, iy=18, nxy=grid.n.x*grid.n.y;
                    // for(int iz=0; iz<grid.n.z; iz++){
                    //     int i= iz*nxy + iy*grid.n.x + ix; 
                    //     printf( "GridFF::tryLoad_new()  Bspline_PLQ[%4i,%4i,%4i|%8i](%g,%g,%g)\n", ix,iy,iz,i, Bspline_PLQ[i].x,Bspline_PLQ[i].y,Bspline_PLQ[i].z);
                    // }
                }
                
                if(done){
                    //pack_Bspline_d();
                    make_inds_pbc( grid.n.x, cubic_xqis );
                    make_inds_pbc( grid.n.y, cubic_yqis );
                    // printf("GridFF::tryLoad_new() BsplineDouble DONE @Bspline_Pauli=%li  @Bspline_London=%li  @Bspline_Coulomb=%li \n", (long)Bspline_Pauli, (long)Bspline_London, (long)Bspline_Coulomb );
                    // golbal_array_dict.insert( { "Bspline_Pauli",   NDArray{ Bspline_Pauli,   Quat4i{ns.z,ns.y,ns.x,-1}} }  );
                    // golbal_array_dict.insert( { "Bspline_London",  NDArray{ Bspline_London,  Quat4i{ns.z,ns.y,ns.x,-1}} }  );
                    // golbal_array_dict.insert( { "Bspline_Coulomb", NDArray{ Bspline_Coulomb, Quat4i{ns.z,ns.y,ns.x,-1}} }  );
                }else{
                    printf( "ERROR GridFF::tryLoad_new( mode=BsplineDouble ) NOT DONE, cubic_xqis,yqis  not initialized !!! => exit(0); \n"); exit(0);
                }
            } break;
        } // switch( mode )
        double T = getCPUticks()-t0;
        //printf( "GridFF::tryLoad_new(mode=%i) DONE (recalc=%i) time %g[Mticks] %g[tick/point]\n", (int)mode, bRecalc, T*1.0e-6, T/(double)npoint );
        return false;
    }

/**
 * Saves the grid data in XSF format for debugging purposes.
 * 
 * @param bE Whether to save the energy-related grids (default: true)
 * @param bFz Whether to save the force-related grids (default: true)
 * @param bComb Whether to save the combined forcefield grid (default: true)
 * @param testREQ The test Quat4d value for evaluating the combined forcefield (default: Quat4d{ 1.487, 0.02609214441, 0., 0.})
 */
    void saveXsfDebug( bool bE=true, bool bFz=true, bool bComb=true, Quat4d testREQ=Quat4d{ 1.487, 0.02609214441, 0., 0.} ){
        // not testREQ.y [eV^0.5] = sqrt(Eii), 
        // e.g. for Hydrogen 0.02609214441 ev^0.5 = sqrt( 0.0006808 eV )
        // e.g. for Carbon   0.06106717612 ev^0.5 = sqrt( 0.0037292 eV )
        if(bE){
            if(FFPaul) grid.saveXSF( "FFLond_E.xsf", (float*)FFLond, 4,3  );
            if(FFLond) grid.saveXSF( "FFelec_E.xsf", (float*)FFelec, 4,3  );
            if(FFelec) grid.saveXSF( "FFPaul_E.xsf", (float*)FFPaul, 4,3  );
        }
        if(bFz){
            if(FFPaul) grid.saveXSF( "FFLond_z.xsf", (float*)FFLond, 4,2  );
            if(FFLond) grid.saveXSF( "FFelec_z.xsf", (float*)FFelec, 4,2  );
            if(FFelec) grid.saveXSF( "FFPaul_z.xsf", (float*)FFPaul, 4,2  );
        }
        // ---- Save combined forcefield
        if(bComb){
            Quat4f * FFtot = new Quat4f[grid.getNtot()];
            evalCombindGridFF ( testREQ, FFtot );
            grid.saveXSF( "E_PLQ.xsf",  (float*)FFtot, 4, 3, natoms, atypes, apos );
            delete [] FFtot;
        }
    }

#endif


}; // class GridFF


#endif
