#ifndef InterpolateTrilinear_h
#define InterpolateTrilinear_h


#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"

#include "datatypes.h"

namespace Trilinear{

//__attribute__((hot))  
template<typename T8>
inline void basis( const Vec3d u, const Vec3i n, int8& qi, T8& qf ){
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
    qi = int8{
        ix+i00, kx+i00,
        ix+i01, kx+i01,
        ix+i10, kx+i10,
        ix+i11, kx+i11 };
    qf = T8{ 
        mx*f00, tx*f00,
        mx*f01, tx*f01,
        mx*f10, tx*f10,
        mx*f11, tx*f11 };
}  

__attribute__((pure))  
__attribute__((hot))  
inline Quat4f lincomb( const int8 qi, const float8 qf, const Quat4f* FE ){
    return (FE[ qi.x  ]*qf.x  ) + (FE[ qi.y  ]*qf.y)
         + (FE[ qi.z  ]*qf.z  ) + (FE[ qi.w  ]*qf.w)  
         + (FE[ qi.hx ]*qf.hx ) + (FE[ qi.hy ]*qf.hy)
         + (FE[ qi.hz ]*qf.hz ) + (FE[ qi.hw ]*qf.hw);
}

__attribute__((pure))  
__attribute__((hot))
inline Quat4d lincomb( const int8 qi, const double8 qf, const Quat4d* FE ){
    return (FE[ qi.x  ]*qf.x  ) + (FE[ qi.y  ]*qf.y)
         + (FE[ qi.z  ]*qf.z  ) + (FE[ qi.w  ]*qf.w)  
         + (FE[ qi.hx ]*qf.hx ) + (FE[ qi.hy ]*qf.hy)
         + (FE[ qi.hz ]*qf.hz ) + (FE[ qi.hw ]*qf.hw);
}

__attribute__((pure))  
__attribute__((hot))  
inline Quat4f fe3f_3( const Vec3d u, const Vec3i n, const Vec3d C, const Quat4f* FE1, const Quat4f* FE2, const Quat4f* FE3 ){    
    int8 qi; float8 qf;
    basis( u, n, qi, qf );
    return  lincomb( qi, qf, FE1 )*C.x
          + lincomb( qi, qf, FE2 )*C.y
          + lincomb( qi, qf, FE3 )*C.z;
}

__attribute__((pure))  
__attribute__((hot))  
inline Quat4d fe3d_3( const Vec3d u, const Vec3i n, const Vec3d C, const Quat4d* FE1, const Quat4d* FE2, const Quat4d* FE3 ){    
    int8 qi; double8 qf;
    basis( u, n, qi, qf );
    return  lincomb( qi, qf, FE1 )*C.x
          + lincomb( qi, qf, FE2 )*C.y
          + lincomb( qi, qf, FE3 )*C.z;
}

__attribute__((pure))
__attribute__((hot))  
inline Quat4d fe3d_fe3( const Vec3d u, const Vec3i n, const Vec3d C, const Quat4d* FE1, const Quat4d* FE2, const Quat4d* FE3 ){    

	const int    ix = (int)u.x  ,  iy = (int)u.y  ,  iz = (int)u.z  ;
    const double tx = u.x - ix  ,  ty = u.y - iy  ,  tz = u.z - iz  ;
    const double mx = 1-tx      ,  my = 1-ty      ,  mz = 1-tz      ;
    //------
    int kx = ix+1; kx=(kx<n.x)?kx:0;
    int ky = iy+1; ky=(ky<n.y)?ky:0;
    int kz = iz+1; kz=(kz<n.z)?kz:0;
	//------
	const double f00 = my*mz; 
    const double f10 = ty*mz; 
    const double f01 = my*tz; 
    const double f11 = ty*tz;
    const int   i00 = n.x*(iy + n.y*iz);
    const int   i01 = n.x*(iy + n.y*kz);
    const int   i10 = n.x*(ky + n.y*iz);
    const int   i11 = n.x*(ky + n.y*kz);
    const int 
        i000=ix+i00, i100=kx+i00,
        i001=ix+i01, i101=kx+i01,
        i010=ix+i10, i110=kx+i10,
        i011=ix+i11, i111=kx+i11;
        
    const double 
        f000=mx*f00, f100=tx*f00,
        f001=mx*f01, f101=tx*f01,
        f010=mx*f10, f110=tx*f10,
        f011=mx*f11, f111=tx*f11;
    // { // DEBUG
    //     Quat4f fDBG = FFPaul[ i000 ];
    //     //printf( "GridFF::getForce() u(%g,%g,%g)/[%3i,%3i,%3i] fe(%g,%g,%g|%g) p0(%g,%g,%g) PLQ(%g,%g,%g)\n", u.x,u.y,u.z, grid.n.x,grid.n.y,grid.n.z, fDBG.x,fDBG.y,fDBG.z,fDBG.w, grid.pos0.x,grid.pos0.y,grid.pos0.z, PLQ.x,PLQ.y,PLQ.z );
    // }
    //printf( "GridFF::getForce() ixyz(%i,%i,%i) nxyz(%i,%i,%i) ntot=%i  ((%i,%i)(%i,%i))((%i,%i)(%i,%i))\n", ix,iy,iz, n.x, n.y, n.z, grid.getNtot(), i000, i001, i010, i011,  i100, i101, i110, i111 );

	return  // 3 * 8 * 4 = 96 floats   // SIMD optimize ?????
          ((FE1[ i000 ]*f000) + (FE1[ i100 ]*f100)
         + (FE1[ i010 ]*f010) + (FE1[ i110 ]*f110)  
         + (FE1[ i011 ]*f011) + (FE1[ i111 ]*f111)
         + (FE1[ i001 ]*f001) + (FE1[ i101 ]*f101))*C.x

         +((FE2[ i000 ]*f000) + (FE2[ i100 ]*f100)
         + (FE2[ i010 ]*f010) + (FE2[ i110 ]*f110)  
         + (FE2[ i011 ]*f011) + (FE2[ i111 ]*f111)
         + (FE2[ i001 ]*f001) + (FE2[ i101 ]*f101))*C.y

         +((FE3[ i000 ]*f000) + (FE3[ i100 ]*f100)
         + (FE3[ i010 ]*f010) + (FE3[ i110 ]*f110)  
         + (FE3[ i011 ]*f011) + (FE3[ i101 ]*f111)
         + (FE3[ i001 ]*f001) + (FE3[ i101 ]*f101))*C.z
        ;
}

__attribute__((pure))
__attribute__((hot))  
inline Quat4f fe3f_fe3( const Vec3f u, const Vec3i n, const Vec3f C, const Quat4f* FE1, const Quat4f* FE2, const Quat4f* FE3 ){    

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
    // { // DEBUG
    //     Quat4f fDBG = FFPaul[ i000 ];
    //     //printf( "GridFF::getForce() u(%g,%g,%g)/[%3i,%3i,%3i] fe(%g,%g,%g|%g) p0(%g,%g,%g) PLQ(%g,%g,%g)\n", u.x,u.y,u.z, grid.n.x,grid.n.y,grid.n.z, fDBG.x,fDBG.y,fDBG.z,fDBG.w, grid.pos0.x,grid.pos0.y,grid.pos0.z, PLQ.x,PLQ.y,PLQ.z );
    // }
    //printf( "GridFF::getForce() ixyz(%i,%i,%i) nxyz(%i,%i,%i) ntot=%i  ((%i,%i)(%i,%i))((%i,%i)(%i,%i))\n", ix,iy,iz, n.x, n.y, n.z, grid.getNtot(), i000, i001, i010, i011,  i100, i101, i110, i111 );

	return  // 3 * 8 * 4 = 96 floats   // SIMD optimize ?????
          ((FE1[ i000 ]*f000) + (FE1[ i100 ]*f100)
         + (FE1[ i010 ]*f010) + (FE1[ i110 ]*f110)  
         + (FE1[ i011 ]*f011) + (FE1[ i111 ]*f111)
         + (FE1[ i001 ]*f001) + (FE1[ i101 ]*f101))*C.x

         +((FE2[ i000 ]*f000) + (FE2[ i100 ]*f100)
         + (FE2[ i010 ]*f010) + (FE2[ i110 ]*f110)  
         + (FE2[ i011 ]*f011) + (FE2[ i111 ]*f111)
         + (FE2[ i001 ]*f001) + (FE2[ i101 ]*f101))*C.y

         +((FE3[ i000 ]*f000) + (FE3[ i100 ]*f100)
         + (FE3[ i010 ]*f010) + (FE3[ i110 ]*f110)  
         + (FE3[ i011 ]*f011) + (FE3[ i101 ]*f111)
         + (FE3[ i001 ]*f001) + (FE3[ i101 ]*f101))*C.z
        ;
}

}; // namespace Trilinear

#endif
