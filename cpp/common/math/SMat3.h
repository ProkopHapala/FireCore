
#ifndef  SMat3_h
#define  SMat3_h

#include <math.h>
#include <cstdlib>
#include <cstdio>
#include <stdint.h>

#include "fastmath.h"
//#include "fastmath_light.h"
#include "Vec3.h"

//template <class T, class VEC, class MAT>
//template <class T, class VEC>
template <class T>
class SMat3{
	using VEC = Vec3T<T>;
	using MAT = SMat3<T>;
	public:
	union{
		struct{
			T xx,yy,zz;
			//T yz,xz,xy;
			union{
			    struct{ T yz,xz,xy; };
			    struct{ T zy,zx,yx; };
			};
		};
		struct{	VEC diag,offd; };
		T array[6];
		VEC  vecs[2];
	};


  // TODO : This should be good also for rotation matrix (?)

    //T xx,yy,zz;
    //union{
    //    struct{ T yz,xz,xy; };
    //    struct{ T zy,zx,yx; };
    //};

	inline void setOne(     ){ xx=yy=zz=1; xy=xz=yz=0; };
	inline void set   ( T f ){ xx=yy=zz=f; xy=xz=yz=0; };

	inline MAT operator* ( T f   ) const { MAT m; m.diag.set_mul(diag,f); m.offd.set_mul(offd,f); return m; };

    inline void mul ( T f           ){ diag.mul(f); offd.mul(f);  };

    inline void add_mul( MAT m, T f ){ diag.add_mul(m.diag,f); offd.add_mul(m.offd,f); }
   
// ====== dot product with vector

	inline VEC dot( const VEC&  v ) const {
		VEC vout;
		vout.x = xx*v.x + xy*v.y + xz*v.z;
		vout.y = xy*v.x + yy*v.y + yz*v.z;
		vout.z = xz*v.x + yz*v.y + zz*v.z;
		return vout;
	}

	inline void dot_to( const VEC&  v, VEC&  vout ) const {
        T vx=v.x,vy=v.y,vz=v.z; // to make it safe use inplace
		vout.x = xx*vx + xy*vy + xz*vz;
		vout.y = xy*vx + yy*vy + yz*vz;
		vout.z = xz*vx + yz*vy + zz*vz;
	};

    inline void fromDeriv( VEC& dij, T sc ){
        //double ir3 = 1/(rij*r2);
        T xx_  = dij.x*dij.x;
        T yy_  = dij.y*dij.y;
        T zz_  = dij.z*dij.z;
        xy=-dij.x*dij.y*sc;
        xz=-dij.x*dij.z*sc;
        yz=-dij.y*dij.z*sc;
        xx=(yy_+zz_)   *sc;
        yy=(xx_+zz_)   *sc;
        zz=(xx_+yy_)   *sc;
    }

    inline void from_outer(const Vec3T<T>& h){
        xy=h.x*h.y;
        xz=h.x*h.z;
        yz=h.y*h.z;
        xx=h.x*h.x;
        yy=h.y*h.y;
        zz=h.z*h.z;
    }

    inline void add_outer(const Vec3T<T>& h){
        xy+=h.x*h.y;
        xz+=h.x*h.z;
        yz+=h.y*h.z;
        xx+=h.x*h.x;
        yy+=h.y*h.y;
        zz+=h.z*h.z;
    }

    inline void from_dhat(const Vec3T<T>& h){
        // derivatives of normalized vector
        //double ir  = irs[i];
        T hxx = h.x*h.x;
        T hyy = h.y*h.y;
        T hzz = h.z*h.z;
        xy=-h.x*h.y;
        xz=-h.x*h.z;
        yz=-h.y*h.z;
        xx=(hyy+hzz);
        yy=(hxx+hzz);
        zz=(hxx+hyy);
    }

    inline void dhat_dot( const Vec3T<T>& h, Vec3T<T>& f )const{
        f.x += h.x*xx + h.y*xy + h.z*xz;
        f.y += h.x*xy + h.y*yy + h.z*yz;
        f.z += h.x*xz + h.y*yz + h.z*zz;
    }

    inline void mad_ddot( const Vec3T<T>& h, Vec3T<T>& f, T k )const{
        f.x += ( h.x*xx + h.y*xy + h.z*xz )*k;
        f.y += ( h.x*xy + h.y*yy + h.z*yz )*k;
        f.z += ( h.x*xz + h.y*yz + h.z*zz )*k;
    }


};

/*
class Mat3i : public Mat3T< int   , Vec3i, Mat3i >{};
class Mat3f : public Mat3T< float , Vec3f, Mat3f >{};
class MAT : public Mat3T< T, VEC, MAT >{};
*/

using SMat3i = SMat3< int   >;
using SMat3f = SMat3< float >;
using SMat3d = SMat3< double>;

static constexpr SMat3d SMat3dIdentity = (SMat3d){ 1.0,1.0,1.0, 0.0,0.0,0.0 };
static constexpr SMat3d SMat3dZero     = (SMat3d){ 0.0,0.0,0.0, 0.0,0.0,0.0 };

static constexpr SMat3f SMat3fIdentity = (SMat3f){1.0f,1.0f,1.0f, 0.0f,0.0f,0.0f };
static constexpr SMat3f SMat3fZero     = (SMat3f){0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f };

inline void convert( const SMat3f& from, SMat3d& to ){ convert( from.diag, to.diag ); convert( from.offd, to.offd ); };
inline void convert( const SMat3d& from, SMat3f& to ){ convert( from.diag, to.diag ); convert( from.offd, to.offd ); };

inline SMat3f toFloat( const SMat3d& from){ SMat3f to; convert( from.diag, to.diag ); convert( from.offd, to.offd ); return to; }

#endif

