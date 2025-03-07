
#ifndef  Mat3_h
#define  Mat3_h

#include <math.h>
#include <cstdlib>
#include <stdint.h>

#include "fastmath.h"
//#include "gonioApprox.h"
#include "Vec3.h"

//template <class T, class VEC, class MAT>
//template <class T, class VEC>
template <class T>
class Mat3T{
	using VEC = Vec3T<T>;
	using MAT = Mat3T<T>;
	public:
	union{
		struct{
			T xx,xy,xz;
			T yx,yy,yz;
			T zx,zy,zz;
		};
		struct{
			T ax,ay,az;
			T bx,by,bz;
			T cx,cy,cz;
		};
		struct{	VEC a,b,c;    };
		struct{	VEC lf,up,fw; };
		T array[9];
		VEC  vecs [3];
	};


// ====== initialization

	inline explicit operator Mat3T<double>()const{ return Mat3T<double>{ (double)xx,(double)xy,(double)xz, (double)yx,(double)yy,(double)yz, (double)zx,(double)zy,(double)zz }; }
	inline explicit operator Mat3T<float >()const{ return Mat3T<float >{ (float )xx,(float )xy,(float )xz, (float )yx,(float )yy,(float )yz, (float )zx,(float )zy,(float )zz }; }
	inline explicit operator Mat3T<int   >()const{ return Mat3T<int   >{ (int   )xx,(int   )xy,(int   )xz, (int   )yx,(int   )yy,(int   )yz, (int   )zx,(int   )zy,(int   )zz }; }

	inline void setOne(        ){ xx=yy=zz=1; xy=xz=yx=yz=zx=zy=0; };
	inline void set   ( T f ){ xx=yy=zz=f; xy=xz=yx=yz=zx=zy=0; };
	inline void set   ( const T& xx_, const T& xy_, const T& xz_, const T& yx_, const T& yy_, const T& yz_, const T& zx_, const T& zy_, const T& zz_ ){ xx=xx_; xy=xy_; xz=xz_; yx=yx_; yy=yy_; yz=yz_; zx=zx_; zy=zy_; zz=zz_; };

	inline void set  ( const VEC& va, const VEC& vb, const VEC& vc ){ a.set(va); b.set(vb); c.set(vc); }
	inline void set  ( const MAT& M ){
		xx=M.xx; xy=M.xy; xz=M.xz;
		yx=M.yx; yy=M.yy; yz=M.yz;
		zx=M.zx; zy=M.zy; zz=M.zz;
	};

	inline void set_outer  ( const VEC& a, const VEC& b ){
		xx=a.x*b.x; xy=a.x*b.y; xz=a.x*b.z;
		yx=a.y*b.x; yy=a.y*b.y; yz=a.y*b.z;
		zx=a.z*b.x; zy=a.z*b.y; zz=a.z*b.z;
	};

    inline void add_outer( const VEC& a, const VEC& b, T f=1.0 ){
        xx+=a.x*b.x*f; xy+=a.x*b.y*f; xz+=a.x*b.z*f;
        yx+=a.y*b.x*f; yy+=a.y*b.y*f; yz+=a.y*b.z*f;
        zx+=a.z*b.x*f; zy+=a.z*b.y*f; zz+=a.z*b.z*f;
    };


	inline void diag_add( T f ){ xx+=f; yy+=f; zz+=f; };

	inline VEC getColx(){ VEC out; out.x = xx; out.y = yx; out.z = zx; return out; };
    inline VEC getColy(){ VEC out; out.x = xy; out.y = yy; out.z = zy; return out; };
    inline VEC getColz(){ VEC out; out.x = xz; out.y = yz; out.z = zz; return out; };

	inline void  colx_to( VEC& out){ out.x = xx; out.y = yx; out.z = zx; };
    inline void  coly_to( VEC& out){ out.x = xy; out.y = yy; out.z = zy; };
    inline void  colz_to( VEC& out){ out.x = xz; out.y = yz; out.z = zz; };

	inline void  setColx( const VEC v ){ xx = v.x; yx = v.y; zx = v.z; };
	inline void  setColy( const VEC v ){ xy = v.x; yy = v.y; zy = v.z; };
	inline void  setColz( const VEC v ){ xz = v.x; yz = v.y; zz = v.z; };

	inline void swap_vecs( const Vec3i& inds ){ MAT M=*this; a=M.vecs[inds.a]; b=M.vecs[inds.b]; c=M.vecs[inds.c]; }
	inline void swap_rows( const Vec3i& inds ){ a.swap(inds); b.swap(inds); c.swap(inds); }

	// Don't need this, because we use union: use representation a,b,c
	//inline VEC getRowx(){ VEC out; out.x = xx; out.y = xy; out.z = xz; return out; };
	//inline VEC getRowy(){ VEC out; out.x = yx; out.y = yy; out.z = yz; return out; };
	//inline VEC getRowz(){ VEC out; out.x = zx; out.y = zy; out.z = zz; return out; };
	//inline void rowx_to( VEC& out ){ out.x = xx; out.y = xy; out.z = xz; };
	//inline void rowy_to( VEC& out ){ out.x = yx; out.y = yy; out.z = yz; };
	//inline void rowz_to( VEC& out ){ out.x = zx; out.y = zy; out.z = zz; };
	//inline void  setRowx( const VEC& v ){ xx = v.x; xy = v.y; xz = v.z; };
	//inline void  setRowy( const VEC& v ){ yx = v.x; yy = v.y; yz = v.z; };
	//inline void  setRowz( const VEC& v ){ zx = v.x; zy = v.y; zz = v.z; };

// ====== transpose

	inline void makeT(){
		T t1=yx; yx=xy; xy=t1;
		T t2=zx; zx=xz; xz=t2;
		T t3=zy; zy=yz; yz=t3;
	};

	inline void setT  ( const MAT& M ){
		xx=M.xx; xy=M.yx; xz=M.zx;
		yx=M.xy; yy=M.yy; yz=M.zy;
		zx=M.xz; zy=M.yz; zz=M.zz;
	};

	inline MAT transposed(){ MAT t; t.setT(*this); return t; }

	inline void setT  ( const VEC& va, const VEC& vb, const VEC& vc ){
		a.set( va.x, vb.x, vc.x );
		b.set( va.y, vb.y, vc.y );
		c.set( va.z, vb.z, vc.z );
	};

	inline MAT operator* ( T f           ) const { MAT m; m.a.set_mul(a,f);   m.b.set_mul(b,f);   m.c.set_mul(c,f);   return m; };
    inline MAT operator+ ( const MAT& m  ) const { MAT o; o.a.set_add(a,m.a); o.b.set_add(b,m.b); o.c.set_add(c,m.c); return o; };

    inline void add    ( const MAT& m      ){ a.add(m.a);       b.add(m.b);       c.add(m.c);  };
    inline void sub    ( const MAT& m      ){ a.sub(m.a);       b.sub(m.b);       c.sub(m.c);  };
    inline void add_mul( const MAT& m , T f){ a.add_mul(m.a,f); b.add_mul(m.b,f); c.add_mul(m.c,f); };

    inline void mul ( T f        ){ a.mul(f);    b.mul(f);    c.mul(f);    };
    inline void mul ( const VEC& va ){ a.mul(va.a); b.mul(va.b); c.mul(va.c); };

    inline void div ( const VEC& va ){ a.mul(1/va.a); b.mul(1/va.b); c.mul(1/va.c); };

    inline void mulT ( const VEC& va ){
		ax*=va.x; ay*=va.y; az*=va.z;
		bx*=va.x; by*=va.y; bz*=va.z;
		cx*=va.x; cy*=va.y; cz*=va.z;
	};

    inline void divT ( const VEC& va ){
        T fx=1/va.x,fy=1/va.y,fz=1/va.z;
		ax*=fx; ay*=fy; az*=fz;
		bx*=fx; by*=fy; bz*=fz;
		cx*=fx; cy*=fy; cz*=fz;
	};


// ====== dot product with vector

	inline VEC dot( const VEC&  v ) const {
		return VEC{
		    xx*v.x + xy*v.y + xz*v.z,
		    yx*v.x + yy*v.y + yz*v.z,
		    zx*v.x + zy*v.y + zz*v.z  };
	}

    inline VEC dotT( const VEC&  v ) const {
		return VEC{
		    xx*v.x + yx*v.y + zx*v.z,
		    xy*v.x + yy*v.y + zy*v.z,
		    xz*v.x + yz*v.y + zz*v.z };
	}
	inline VEC lincomb( T fx, T fy, T fz )const{  return dotT({fx,fy,fz}); }

	inline void dot_to( const VEC&  v, VEC&  vout ) const {
        T vx=v.x,vy=v.y,vz=v.z; // to make it safe use inplace
		vout.x = xx*vx + xy*vy + xz*vz;
		vout.y = yx*vx + yy*vy + yz*vz;
		vout.z = zx*vx + zy*vy + zz*vz;
	};

	inline void dot_to_T( const VEC&  v, VEC&  vout ) const {
        T vx=v.x,vy=v.y,vz=v.z;
		vout.x = xx*vx + yx*vy + zx*vz;
		vout.y = xy*vx + yy*vy + zy*vz;
		vout.z = xz*vx + yz*vy + zz*vz;
	};

    inline bool tryOrthoNormalize( double errMax, int ia, int ib, int ic ){
        VEC& a = vecs[ia];
        VEC& b = vecs[ib];
        VEC& c = vecs[ic];
        bool res = false;
        res |= a.tryNormalize    ( errMax );
        res |= b.tryOrthogonalize( errMax, a );
        res |= b.tryNormalize    ( errMax );
        res |= c.tryOrthogonalize( errMax, a );
        res |= c.tryOrthogonalize( errMax, b );
        res |= c.tryNormalize    ( errMax );
		return res;
	};


    inline void orthogonalize( int ia, int ib, int ic ){
        VEC& a = vecs[ia];
        VEC& b = vecs[ib];
        VEC& c = vecs[ic];
        a.normalize ();
        b.makeOrthoU(a);
        b.normalize ();
        c.makeOrthoU(a);
        c.makeOrthoU(b);
        c.normalize();
	};

	inline void orthogonalize_taylor3( int ia, int ib, int ic ){
        VEC& a = vecs[ia];
        VEC& b = vecs[ib];
        VEC& c = vecs[ic];
        a.normalize_taylor3();
        b.makeOrthoU(a);
        b.normalize_taylor3();
        c.makeOrthoU(a);
        c.makeOrthoU(b);
        c.normalize_taylor3();
	};


// ====== matrix multiplication

	inline void set_mmul( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.xy*B.yx + A.xz*B.zx;
		xy = A.xx*B.xy + A.xy*B.yy + A.xz*B.zy;
		xz = A.xx*B.xz + A.xy*B.yz + A.xz*B.zz;
		yx = A.yx*B.xx + A.yy*B.yx + A.yz*B.zx;
		yy = A.yx*B.xy + A.yy*B.yy + A.yz*B.zy;
		yz = A.yx*B.xz + A.yy*B.yz + A.yz*B.zz;
		zx = A.zx*B.xx + A.zy*B.yx + A.zz*B.zx;
		zy = A.zx*B.xy + A.zy*B.yy + A.zz*B.zy;
		zz = A.zx*B.xz + A.zy*B.yz + A.zz*B.zz;
	};

	inline void set_mmul_NT( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.xy*B.xy + A.xz*B.xz;
		xy = A.xx*B.yx + A.xy*B.yy + A.xz*B.yz;
		xz = A.xx*B.zx + A.xy*B.zy + A.xz*B.zz;
		yx = A.yx*B.xx + A.yy*B.xy + A.yz*B.xz;
		yy = A.yx*B.yx + A.yy*B.yy + A.yz*B.yz;
		yz = A.yx*B.zx + A.yy*B.zy + A.yz*B.zz;
		zx = A.zx*B.xx + A.zy*B.xy + A.zz*B.xz;
		zy = A.zx*B.yx + A.zy*B.yy + A.zz*B.yz;
		zz = A.zx*B.zx + A.zy*B.zy + A.zz*B.zz;
	};

	inline void set_mmul_TN( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.yx*B.yx + A.zx*B.zx;
		xy = A.xx*B.xy + A.yx*B.yy + A.zx*B.zy;
		xz = A.xx*B.xz + A.yx*B.yz + A.zx*B.zz;
		yx = A.xy*B.xx + A.yy*B.yx + A.zy*B.zx;
		yy = A.xy*B.xy + A.yy*B.yy + A.zy*B.zy;
		yz = A.xy*B.xz + A.yy*B.yz + A.zy*B.zz;
		zx = A.xz*B.xx + A.yz*B.yx + A.zz*B.zx;
		zy = A.xz*B.xy + A.yz*B.yy + A.zz*B.zy;
		zz = A.xz*B.xz + A.yz*B.yz + A.zz*B.zz;
	};

	inline void set_mmul_TT( const MAT& A, const MAT& B ){
		xx = A.xx*B.xx + A.yx*B.xy + A.zx*B.xz;
		xy = A.xx*B.yx + A.yx*B.yy + A.zx*B.yz;
		xz = A.xx*B.zx + A.yx*B.zy + A.zx*B.zz;
		yx = A.xy*B.xx + A.yy*B.xy + A.zy*B.xz;
		yy = A.xy*B.yx + A.yy*B.yy + A.zy*B.yz;
		yz = A.xy*B.zx + A.yy*B.zy + A.zy*B.zz;
		zx = A.xz*B.xx + A.yz*B.xy + A.zz*B.xz;
		zy = A.xz*B.yx + A.yz*B.yy + A.zz*B.yz;
		zz = A.xz*B.zx + A.yz*B.zy + A.zz*B.zz;
	};

// ====== matrix solver

   inline T determinant()const{
        T fCoxx = yy * zz - yz * zy;
        T fCoyx = yz * zx - yx * zz;
        T fCozx = yx * zy - yy * zx;
        T fDet  = xx * fCoxx + xy * fCoyx + xz * fCozx;
        return fDet;
    };

	inline void invert_to( MAT& Mout ) const{
        T idet = 1/determinant(); // we dont check det|M|=0
		//printf("Mat3d::invert_to() idet = %g \n", idet);
        Mout.xx = ( yy * zz - yz * zy ) * idet;
        Mout.xy = ( xz * zy - xy * zz ) * idet;
        Mout.xz = ( xy * yz - xz * yy ) * idet;
        Mout.yx = ( yz * zx - yx * zz ) * idet;
        Mout.yy = ( xx * zz - xz * zx ) * idet;
        Mout.yz = ( xz * yx - xx * yz ) * idet;
        Mout.zx = ( yx * zy - yy * zx ) * idet;
        Mout.zy = ( xy * zx - xx * zy ) * idet;
        Mout.zz = ( xx * yy - xy * yx ) * idet;
    };

    inline void invert_T_to( MAT& Mout ) const{
        T idet = 1/determinant(); // we dont check det|M|=0
        Mout.xx = ( yy * zz - yz * zy ) * idet;
        Mout.yx = ( xz * zy - xy * zz ) * idet;
        Mout.zx = ( xy * yz - xz * yy ) * idet;
        Mout.xy = ( yz * zx - yx * zz ) * idet;
        Mout.yy = ( xx * zz - xz * zx ) * idet;
        Mout.zy = ( xz * yx - xx * yz ) * idet;
        Mout.xz = ( yx * zy - yy * zx ) * idet;
        Mout.yz = ( xy * zx - xx * zy ) * idet;
        Mout.zz = ( xx * yy - xy * yx ) * idet;
    };

    inline void adjoint_to( MAT& Mout ) const{
        Mout.xx = yy * zz - yz * zy;
        Mout.xy = xz * zy - xy * zz;
        Mout.xz = xy * yz - xz * yy;
        Mout.yx = yz * zx - yx * zz;
        Mout.yy = xx * zz - xz * zx;
        Mout.yz = xz * yx - xx * yz;
        Mout.zx = yx * zy - yy * zx;
        Mout.zy = xy * zx - xx * zy;
        Mout.zz = xx * yy - xy * yx;
    };

// ======= Rotation

	inline void rotate( T angle, VEC axis  ){
		//VEC uaxis;
		//uaxis.set( axis * axis.norm() );
		axis.normalize();
		T ca   = cos(angle);
		T sa   = sin(angle);
 		rotate_csa( ca, sa, axis );
	};

	inline void rotate_csa( T ca, T sa, const VEC& uaxis ){
		a.rotate_csa( ca, sa, uaxis );
		b.rotate_csa( ca, sa, uaxis );
		c.rotate_csa( ca, sa, uaxis );
		//a.set(1);
		//b.set(2);
		//c.set(3);
	};

	inline void drotate_omega6( const VEC& w ){
        // consider not-normalized vector omega
        T ca,sa;
        sincosR2_taylor(w.norm2(), sa, ca );
        a.drotate_omega_csa(w,ca,sa);
        b.drotate_omega_csa(w,ca,sa);
        c.drotate_omega_csa(w,ca,sa);
	};

	void dRotateToward( int pivot, const MAT& rot0, T dPhi ){
        int i3 = pivot*3;
        VEC& piv  = *(VEC*)(     array+i3);
        VEC& piv0 = *(VEC*)(rot0.array+i3);
        VEC ax; ax.set_cross(piv,piv0);
        T sa = ax.norm();
        if( sa > dPhi ){
            ax.mul(1.0/sa);
            Vec2d csa; csa.fromAngle( dPhi );
            rotate_csa( csa.x, csa.y, ax );
        }else{
            set(rot0);
        }
    }

	// ==== generation


	inline void fromRotation( T angle, VEC axis ){ setOne(); rotate(angle,axis); };

	inline void fromDirUp( const VEC& dir, const VEC& up ){
		// make orthonormal rotation matrix c=dir; b=(up-<b|c>c)/|b|; a=(c x b)/|a|;
		c.set(dir);
		c.normalize(); // we assume dir is already normalized
		b.set(up);
		b.add_mul( c, -b.dot(c) );   //
		b.normalize();
		a.set_cross(b,c);
		//a.normalize(); // we don't need this since b,c are orthonormal
	};

    inline void fromSideUp( const VEC& side, const VEC&  up ){
		// make orthonormal rotation matrix c=dir; b=(up-<b|c>c)/|b|; a=(c x b)/|a|;
		a.set(side);
		//c.normalize(); // we assume dir is already normalized
		b.set(up);
		b.add_mul( a, -b.dot(a) );   //
		b.normalize();
		c.set_cross(b,a);
		//a.normalize(); // we don't need this since b,c are orthonormal
	};

	inline void fromCrossSafe( const Vec3d& v1, const Vec3d& v2 ){
        b.set_cross( v1, v2 );
        a.set_sub(v2,v1); a.normalize();
        double r2b = b.norm2();
        if( r2b<1e-15 ){
            a.getSomeOrtho(b,c);
        }else{
            b.mul( 1/sqrt(r2b) );
            c.set_cross(b,a);
        }
	}

/**
 * @brief Sets the matrix to a rotation matrix that represents the given Euler angles.
 * 
 * The rotation is applied in the order of phi (roll), theta (pitch), and psi (yaw).
 * @param phi The roll angle in radians.
 * @param theta The pitch angle in radians.
 * @param psi The yaw angle in radians.
 */
	inline void fromEuler( T phi, T theta, T psi ){
        // http://mathworld.wolfram.com/EulerAngles.html
        T ca=1,sa=0, cb=1,sb=0, cc=1,sc=0;
        //if(phi*phi    >1e-16){ ca=cos(phi);   sa=sin(phi); }
        //if(theta*theta>1e-16){ cb=cos(theta); sb=sin(theta); }
        //if(psi*psi    >1e-16){ cc=cos(psi);   sc=sin(psi); }
        ca=cos(phi);   sa=sin(phi);
        cb=cos(theta); sb=sin(theta);
        cc=cos(psi);   sc=sin(psi);
        /*
        xx =  cc*ca-cb*sa*sc;
		xy =  cc*sa+cb*ca*sc;
		xz =  sc*sb;
		yx = -sc*ca-cb*sa*cc;
		yy = -sc*sa+cb*ca*cc;
		yz =  cc*sb;
		zx =  sb*sa;
		zy = -sb*ca;
		zz =  cb;
		*/

        xx =  cc*ca-cb*sa*sc;
		xy =  cc*sa+cb*ca*sc;
		xz =  sc*sb;
		zx = -sc*ca-cb*sa*cc;
		zy = -sc*sa+cb*ca*cc;
		zz =  cc*sb;
		yx =  sb*sa;
		yy = -sb*ca;
		yz =  cb;
	};


    /**
     * @brief Sets the matrix to a rotation matrix that corresponds to the given Euler angles in the order of inclination, longitude ascending node, and argument of periapsis.
     * 
     * @param inc The inclination angle in radians.
     * @param lan The longitude ascending node angle in radians.
     * @param apa The argument of periapsis angle in radians.
     */
    inline void fromEuler_orb( T inc, T lan, T apa ){
        T ci=cos(inc), si=sin(inc); // inc  inclination
        T cl=cos(lan), sl=sin(lan); // lan  longitude ascedning node  capital omega
        T ca=cos(apa), sa=sin(apa); // apa  argument of periapsis     small omega
        xx =  cl*ca - sl*sa*ci;
		xy =  sl*ca + cl*sa*ci;
		xz =  sa*si;
		yx = -cl*sa - sl*ca*ci;
		yy = -sl*sa + cl*ca*ci;
		yz =  ca*si;
		zx =  sl*si;
		zy = -cl*si;
		zz =  ci;
	};


     //  RAND _ ROTATION   Author: Jim Arvo, 1991
     //  This routine maps three values (x[0], x[1], x[2]) in the range [0,1]
     //  into a 3x3 rotation matrix, M.  Uniformly distributed random variables
     //  x0, x1, and x2 create uniformly distributed random rotation matrices.
     //  To create small uniformly distributed "perturbations", supply
     //  samples in the following ranges
     //      x[0] in [ 0, d ]
     //      x[1] in [ 0, 1 ]
     //      x[2] in [ 0, d ]
     // where 0 < d < 1 controls the size of the perturbation.  Any of the
     // random variables may be stratified (or "jittered") for a slightly more
     // even distribution.
     //=========================================================================
/**
 * @brief Generates a random rotation matrix using the given random vector.
 * 
 * The random vector is used to distribute points over the sphere via the reflection I - V Transpose(V).
 * This formulation of V will guarantee that if x[1] and x[2] are uniformly distributed, the reflected points will be uniform on the sphere.
 * 
 *  see
 * 	http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
 *  http://www.realtimerendering.com/resources/GraphicsGems/gemsiii/rand_rotation.c
 *  http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.53.1357&rep=rep1&type=pdf
 * 
 * @param vrand The random vector used to generate the rotation matrix.
 */
	inline void fromRand( const VEC& vrand  ){
        T theta = vrand.x * M_TWO_PI; // Rotation about the pole (Z).
        T phi   = vrand.y * M_TWO_PI; // For direction of pole deflection.
        T z     = vrand.z * 2.0;      // For magnitude of pole deflection.
        // Compute a vector V used for distributing points over the sphere
        // via the reflection I - V Transpose(V).  This formulation of V
        // will guarantee that if x[1] and x[2] are uniformly distributed,
        // the reflected points will be uniform on the sphere.  Note that V
        // has length sqrt(2) to eliminate the 2 in the Householder matrix.
        T r  = sqrt( z );
        T Vx = sin ( phi ) * r;
        T Vy = cos ( phi ) * r;
        T Vz = sqrt( 2.0 - z );
        // Compute the row vector S = Transpose(V) * R, where R is a simple
        // rotation by theta about the z-axis.  No need to compute Sz since
        // it's just Vz.
        T st = sin( theta );
        T ct = cos( theta );
        T Sx = Vx * ct - Vy * st;
        T Sy = Vx * st + Vy * ct;
        // Construct the rotation matrix  ( V Transpose(V) - I ) R, which
        // is equivalent to V S - R.
        xx = Vx * Sx - ct;   xy = Vx * Sy - st;   xz = Vx * Vz;
        yx = Vy * Sx + st;   yy = Vy * Sy - ct;   yz = Vy * Vz;
        zx = Vz * Sx;        zy = Vz * Sy;        zz = 1.0 - z;   // This equals Vz * Vz - 1.0
	}

    // took from here
    // Smith, Oliver K. (April 1961), "Eigenvalues of a symmetric 3 × 3 matrix.", Communications of the ACM 4 (4): 168
    // http://www.geometrictools.com/Documentation/EigenSymmetric3x3.pdf
    // https://www.geometrictools.com/GTEngine/Include/Mathematics/GteSymmetricEigensolver3x3.h
/**
 * @brief Calculates the eigenvalues of a 3x3 matrix and stores them in the given evs vector. 
 * 
 * took from here:
 *   Smith, Oliver K. (April 1961), "Eigenvalues of a symmetric 3 × 3 matrix.", Communications of the ACM 4 (4): 168
 *   http://www.geometrictools.com/Documentation/EigenSymmetric3x3.pdf
 *   https://www.geometrictools.com/GTEngine/Include/Mathematics/GteSymmetricEigensolver3x3.h
 * 
 * @param evs (out) vector to store the eigenvalues in.
 */
	inline void eigenvals( VEC& evs, bool bSort=true ) const {
		const T inv3  = 0.33333333333;
        const T root3 = 1.73205080757;
		T amax = array[0];
		for(int i=1; i<9; i++){ double a=array[i]; if(a>amax)amax=a; }
		T c0 = xx*yy*zz + 2*xy*xz*yz -  xx*yz*yz   - yy*xz*xz   -  zz*xy*xy;
		T c1 = xx*yy - xy*xy + xx*zz - xz*xz + yy*zz - yz*yz;
		T c2 = xx + yy + zz;
		T amax2 = amax*amax; c2/=amax; c1/=amax2; c0/=(amax2*amax);
		T c2Div3 = c2*inv3;
		T aDiv3  = (c1 - c2*c2Div3)*inv3;
		if (aDiv3 > 0.0) aDiv3 = 0.0;
		T mbDiv2 = 0.5*( c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1) );
		T q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
		if (q > 0.0) q = 0.0;
		T magnitude = sqrt(-aDiv3);
		T angle = atan2( sqrt(-q), mbDiv2 ) * inv3;
		T cs    = cos(angle);
		T sn    = sin(angle);
		// are the eigenvalues already sorted ?  - I think they are, from 
		evs.a = amax*( c2Div3 + 2.0*magnitude*cs );
		evs.b = amax*( c2Div3 - magnitude*(cs + root3*sn) );
		evs.c = amax*( c2Div3 - magnitude*(cs - root3*sn) );
		if(bSort){
			_order(evs.a,evs.b);
			_order(evs.b,evs.c);
			_order(evs.a,evs.b);
		}
	}

/**
 * Calculates the eigenvector corresponding to the given eigenvalue of the 3x3 matrix.
 * @param eval (in)  The eigenvalue to find the eigenvector for.
 * @param evec (out) The resulting eigenvector.
 */
	inline void eigenvec( T eval, VEC& evec ) const{
		VEC row0;  row0.set( ax - eval, ay, az );
		VEC row1;  row1.set( bx, by - eval, bz );
		VEC row2;  row2.set( cx, cy,  cz- eval );
		VEC r0xr1; r0xr1.set_cross(row0, row1);
		VEC r0xr2; r0xr2.set_cross(row0, row2);
		VEC r1xr2; r1xr2.set_cross(row1, row2);
		T d0 = r0xr1.dot( r0xr1);
		T d1 = r0xr2.dot( r0xr2);
		T d2 = r1xr2.dot( r1xr2);
		T dmax = d0; int imax = 0;
		if (d1 > dmax) { dmax = d1; imax = 1; }
		if (d2 > dmax) { imax = 2;            }
		if      (imax == 0) { evec.set_mul( r0xr1, 1/sqrt(d0) ); }
		else if (imax == 1) { evec.set_mul( r0xr2, 1/sqrt(d1) ); }
		else                { evec.set_mul( r1xr2, 1/sqrt(d2) ); }
	}

	inline void gaussElimination()
	{
		double tolerance = 1e-10;
		Vec3d temp;
		if (abs(xx) < tolerance)
		{
			temp = a;
			a = b;
			b = temp;
			if (abs(xx) < tolerance)
			{
				temp = a;
				a = c;
				c = temp;
			}
		}
		if (abs(xx) > tolerance)
		{
			a.mul(1 / xx);
			temp.set_mul(a, yx);
			b.sub(temp);
			temp.set_mul(a, zx);
			c.sub(temp);
		}
		else
		{
			temp = a;
			a = b;
			b = c;
			c = temp;
		}
		if (abs(yy) < tolerance)
		{
			temp = b;
			b = c;
			c = temp;
		}
		if (abs(yy) > tolerance)
		{
			b.mul(1 / yy);
			temp.set_mul(b, zy);
			c.sub(temp);
			if(abs(xx) < tolerance){
				temp.set_mul(b, xy);
				a.sub(temp);
				temp = a;
				a = b;
				b = temp;
			}
		}
		else if(abs(zz) > tolerance && abs(yz) > tolerance)
		{
				yz = 1;
				zz = 0;
		}
		if(abs(zz) > tolerance)
		{
			c.mul(1/zz);
		}
	}
	inline int backPropagation(Mat3T& v){
		int count = 1;
		double tolerance = 1e-10;
		v.a.z = 1;
		if(abs(yy) > tolerance){
			v.a.y = -yz/yy;
			a.sub(b.mul(xy/yy));
		}
		else if(abs(yz) > tolerance){
			v.a.z = 0;
			v.a.y = 1;
		}
		else{
			count = 2;
			v.a.y = 1;
			v.a.z = 0;
			v.b.y = 0;
			v.b.z = 1;
		}
		if(abs(xx) > tolerance){
			v.a.x = -v.a.z*xz/xx;
			if(abs(xz) < tolerance){
				v.a.x = -xy/xx;
			}
			if(count > 1){
				v.b.x = -xy/xx;
			}
		}
		else if(abs(xy) > tolerance){
			v.a.x = 1;
			v.b.x = 1;
			v.a.y = -v.a.z*xz/xy;
			v.b.y = -v.b.z*xz/xy;
		}else if(abs(xz) > tolerance){
			v.a.x = 1;
			v.b.x = 1;
			v.a.y = 0;
			v.b.y = 1;
			v.a.z = 0;
			v.b.z = 0;
		}
		else{
			count = 3;
			v.set({1,0,0,0,1,0,0,0,1});
		}
		return count;
	}
/**
 * Calculates the eigenvectors of the 3x3 hermitian matrix and stores them in the given evecs matrix.
 * https://hal.science/hal-01501221/document
*/
	inline void eigenvec_and_eigenvals( VEC& evals, MAT& evecs ) const{
		double tolerance = 1e-10;

		if (abs(xy) < tolerance && abs(xz) < tolerance && abs(yz) < tolerance)
		{
			evals = {xx,yy,zz};
			
			evecs.a = {1,0,0};
			evecs.b = {0,1,0};
			evecs.c = {0,0,1};
			return;
		}

		double x1 = xx*xx+yy*yy+zz*zz-xx*yy-yy*zz-zz*xx+3*(xy*xy+yz*yz+xz*xz);
		double x2 = 0-(2*xx-yy-zz)*(2*yy-zz-xx)*(2*zz-xx-yy)+9*((2*zz-xx-yy)*xy*xy+(2*yy-xx-zz)*xz*xz+(2*xx-yy-zz)*yz*yz)-54*xy*xz*yz;
		double phi;
		if (x2 > tolerance)
		{
			phi = atan(sqrt(4 * x1 * x1 * x1 - x2 * x2) / x2);
		}
		else if (x2 < -tolerance)
		{
			phi = atan(sqrt(4 * x1 * x1 * x1 - x2 * x2) / x2) + M_PI;
		}
		else
		{
			phi = M_PI / 2;
		}
		if(4 * x1 * x1 * x1 - x2 * x2 < tolerance) phi = M_PI;
		evals.x = (xx+yy+zz-2*sqrt(x1)*cos(phi/3))/3;
		evals.y = (xx+yy+zz+2*sqrt(x1)*cos((phi-M_PI)/3))/3;
		evals.z = (xx+yy+zz+2*sqrt(x1)*cos((phi+M_PI)/3))/3;

		Mat3T temp, W;
		temp.set( xx-evals.x, xy, xz, yx, yy-evals.x, yz, zx, zy, zz-evals.x );
		temp.gaussElimination();
		int rank = temp.backPropagation(W);
		evecs.a = W.a;


		temp.set( xx-evals.y, xy, xz, yx, yy-evals.y, yz, zx, zy, zz-evals.y );
		temp.gaussElimination();
		rank = temp.backPropagation(W);

		evecs.b = W.a;


		temp.set( xx-evals.z, xy, xz, yx, yy-evals.z, yz, zx, zy, zz-evals.z );
		temp.gaussElimination();
		rank = temp.backPropagation(W);
		evecs.c = W.a;



		




		//

		// if(abs(xy) < tolerance && abs(xz) < tolerance){
		// 	printf("xy and xz are zero\n");
		// 	evals.x = xx;
		// 	evecs.a = {1,0,0};
		// 	double D = sqrt(4*yz*yz+(yy-zz)*(yy-zz));
		// 	evals.y = (yy+zz-D)/2;
		// 	evals.z = (yy+zz+D)/2;
		// 	evecs.b = {(yy-zz+D)/(2*yz), 1, 0};
		// 	evecs.c = {(yy-zz-D)/(2*yz), 1, 0};
		// 	return;
		// }
		// if(abs(xz) < tolerance && abs(yz) < tolerance){
		// 	printf("xz and yz are zero\n");
		// 	evals.z = zz;
		// 	evecs.c = {0,0,1};
		// 	double D = sqrt(4*xy*xy+(xx-yy)*(xx-yy));
		// 	evals.y = (xx+yy-D)/2;
		// 	evals.x = (xx+yy+D)/2;
		// 	evecs.b = {(xx-yy+D)/(2*xy), 1, 0};
		// 	evecs.a = {(xx-yy-D)/(2*xy), 1, 0};
		// 	return;
		// }
		// if(abs(xz) < tolerance)
		// {
		// 	printf("xz is zero\n");
		// 	eigenvals(evals);
		// 	eigenvec(evals.a, evecs.a);
		// 	eigenvec(evals.b, evecs.b);
		// 	eigenvec(evals.c, evecs.c);
		// 	return;
		// }

		// if(xy < tolerance && xz < tolerance && yz < tolerance)
		// {
		// 	evals.x = xx;
		// 	evals.y = yy;
		// 	evals.z = zz;
		// 	evecs.a = {1,0,0};
		// 	evecs.b = {0,1,0};
		// 	evecs.c = {0,0,1};
		// 	return;
		// }




		// double m1 = (xy*(zz-evals.x)-yz*xz)/(xz*(yy-evals.x)-xy*yz);
		// double m2 = (xy*(zz-evals.y)-yz*xz)/(xz*(yy-evals.y)-xy*yz);
		// double m3 = (xy*(zz-evals.z)-yz*xz)/(xz*(yy-evals.z)-xy*yz);

		// evecs.a = {(evals.x-zz-yz*m1)/xz, m1, 1};
		// evecs.b = {(evals.y-zz-yz*m2)/xz, m2, 1};
		// evecs.c = {(evals.z-zz-yz*m3)/xz, m3, 1};
	}

	inline void print() const {
        printf( " %f %f %f \n", ax, ay, az );
        printf( " %f %f %f \n", bx, by, bz );
        printf( " %f %f %f \n", cx, cy, cz );
    }

    inline void printOrtho() const { printf( " %f %f %f   %e %e %e \n", a.norm2(),b.norm2(),c.norm2(),   a.dot(b),a.dot(c),b.dot(c) ); }
    inline void printOrthoErr() const { printf( " %e %e %e   %e %e %e \n", a.norm()-1,b.norm()-1,c.norm()-1,   a.dot(b),a.dot(c),b.dot(c) ); }

    inline void transformVectors( int n, Vec3T<T>* v0s, Vec3T<T>* vs )const{
        for( int j=0; j<n; j++ ){
            Vec3T<T> v;
            //mrot.dot_to_T( h0s[j], h );
            dot_to( v0s[j], v );
            vs[j] = v;
            //ps[j].set_add_mul( pos, p_, r0 );
        }
    }

    inline void transformPoints0( int n, Vec3T<T>* v0s, Vec3T<T>* ps, const Vec3T<T>& toPos )const{
        for( int j=0; j<n; j++ ){
            Vec3T<T> v;
            //mrot.dot_to_T( apos0[j], v );
            dot_to( v0s[j], v );
            ps[j].set_add( v, toPos );
            //printf( "frag2atoms[%i]  (%g,%g,%g) (%g,%g,%g) \n", j,  apos0[j].x, apos0[j].y, apos0[j].z,   apos[j].x, apos[j].y, apos[j].z  );
            //printf( "%i %i  (%g,%g,%g) (%g,%g,%g) \n", ifrag, j,  m_apos[j].x, m_apos[j].y, m_apos[j].z,   Tp.x, Tp.y, Tp.z  );
        }
    }

    inline void transformPoints( int n, Vec3T<T>* p0s, Vec3T<T>* ps, const Vec3T<T>& pos0 )const{
        for( int j=0; j<n; j++ ){
            Vec3T<T> v0,v;
            v0.set_sub( p0s[j], pos0 );
            dot_to( v0, v );
            ps[j].set_add( pos0, v );
            //printf( "frag2atoms[%i]  (%g,%g,%g) (%g,%g,%g) \n", j,  apos0[j].x, apos0[j].y, apos0[j].z,   apos[j].x, apos[j].y, apos[j].z  );
            //printf( "%i %i  (%g,%g,%g) (%g,%g,%g) \n", ifrag, j,  m_apos[j].x, m_apos[j].y, m_apos[j].z,   Tp.x, Tp.y, Tp.z  );
        }
    }

    inline void scalePoint ( const Vec3T<T>& p0, Vec3T<T>& p, const Vec3T<T>& pos0, const Vec3T<T>& sc )const{
        Vec3T<T> v,v_;
        v.set_sub( p0, pos0 );
        dot_to   ( v, v_  );
        v_.mul   ( sc     );
        dot_to_T ( v_, v  );
        p.set_add( v, pos0);
    };

    inline void scalePoints( int n, Vec3T<T>* p0s, Vec3T<T>* ps, const Vec3T<T>& pos0, const Vec3T<T>& sc )const{
        for( int j=0; j<n; j++ ){ scalePoint( p0s[j], ps[j], pos0, sc ); }
    }

    inline void scalePoints( int n, int* selection, Vec3T<T>* p0s, Vec3T<T>* ps, const Vec3T<T>& pos0, const Vec3T<T>& sc )const{
        for( int j=0; j<n; j++ ){ int i=selection[j]; scalePoint( p0s[i], ps[i], pos0, sc ); }
    }

    inline void addOuter( const VEC& v1, const VEC& v2, T f ){
        xx+=v1.x*v2.x*f; xy+=v1.x*v2.y*f; xz+=v1.x*v2.z*f;
        yx+=v1.y*v2.x*f; yy+=v1.y*v2.y*f; yz+=v1.y*v2.z*f;
        zx+=v1.z*v2.x*f; zy+=v1.z*v2.y*f; zz+=v1.z*v2.z*f;
    }

    inline void setOuter( const VEC& v1, const VEC& v2, T f ){
        xx=v1.x*v2.x*f; xy=v1.x*v2.y*f; xz=v1.x*v2.z*f;
        yx=v1.y*v2.x*f; yy=v1.y*v2.y*f; yz=v1.y*v2.z*f;
        zx=v1.z*v2.x*f; zy=v1.z*v2.y*f; zz=v1.z*v2.z*f;
    }

	inline Vec3T<int> nearestCell( const VEC& d, const int ncellMax=1 ){
		VEC u;
		T off = ncellMax + 0.5;
		dot_to( d, u );
		return Vec3T<int>{
			ncellMax - (int)(u.a+off),
			ncellMax - (int)(u.b+off),
			ncellMax - (int)(u.c+off)
		};
	}

	/*
	template<typename T>
	inline Vec3i nearestCell( const Mat3T<T>& invLvec, const Vec3T<T>& d ){
		Vec3T<T> u;
		invLvec.dot_to( d, u );
		return Vec3i{
			1 -(int)(u.a+1.5),
			1 -(int)(u.b+1.5),
			1 -(int)(u.c+1.5)
		};
	}
	*/

	inline void mat_sqrt( const MAT& A ){
		xx = std::sqrt(A.xx); xy = std::sqrt(A.xy); xz = std::sqrt(A.xz);
		yx = std::sqrt(A.yx); yy = std::sqrt(A.yy); yz = std::sqrt(A.yz);
		zx = std::sqrt(A.zx); zy = std::sqrt(A.zy); zz = std::sqrt(A.zz);
	}

	inline void symetrize_matrix_SVD(MAT& rot){
		MAT A;
		A.set(*this);
		
		double tolerance = 1e-6;
		int iteration_count = 0;
        bool do51 = true; //"This is a way to deal with those nasty gotos in the FORTRAN code"
        int iflag;
        int ix;
		int iy;
		int iz;
        double sigma;
        double gamma;
        double sg;
        Vec3d bb, cc;
        while (true)
        {
            if (do51)
            {
                iflag = 0;
                ix = 0;
            }

            // If the number of iterations exceeds 500, give up
            ++iteration_count;
            if (iteration_count > 100)
            {
                break;
            }

			iy=ix+1;
			if(iy>2) iy=0;
			iz=3-ix-iy;

			sigma = A.vecs[iz].array[iy] - A.vecs[iy].array[iz];
			gamma = A.vecs[iy].array[iy] + A.vecs[iz].array[iz];

            sg = sqrt(sigma * sigma + gamma * gamma);

            if (sg == 0)
            {
                ++ix;
                if (iflag == 0){break;}
                if (ix < 3){do51 = false;}
				else{do51 = true;}
                continue;
            }

            sg = 1.0 / sg;
            if (fabs(sigma) < (tolerance * fabs(gamma)))
			{
                ++ix;
                if (iflag == 0){break;}
                if (ix < 3){do51 = false;}
                else{do51 = true;}
                continue;
            }

			bb.set_add_mul(Vec3dZero, A.vecs[iy], gamma);
			bb.add_mul(A.vecs[iz], sigma);
			bb.mul(sg);
			cc.set_add_mul(Vec3dZero, A.vecs[iz], gamma);
			cc.add_mul(A.vecs[iy], 0-sigma);
			cc.mul(sg);
			A.vecs[iy] = bb;
			A.vecs[iz] = cc;
			
			bb.set_add_mul(Vec3dZero, rot.vecs[iy], gamma);
			bb.add_mul(rot.vecs[iz], sigma);
			bb.mul(sg);
			cc.set_add_mul(Vec3dZero, rot.vecs[iz], gamma);
			cc.add_mul(rot.vecs[iy], 0-sigma);
			cc.mul(sg);
			rot.vecs[iy] = bb;
			rot.vecs[iz] = cc;

            iflag = 1;

            ++ix;
            if (iflag == 0){break;}
            if (ix < 3){do51 = false;}
            else{do51 = true;}
            continue;
        } // End while loop
	}

	inline void SVD(MAT &U, VEC &val, MAT &V)
	{
		MAT A, B;
		A.set(*this);

		B.set_mmul_TN(A, A);
		//B.diag_add(1e-8);
		B.eigenvec_and_eigenvals(val, V);

		V.orthogonalize(0, 1, 2);
		V.c.set_cross(V.a, V.b);

		val.x = sqrt(abs(val.x));
		val.y = sqrt(abs(val.y));
		val.z = sqrt(abs(val.z));
		int zeroValue = -1;
		for (int i = 0; i < 3; i++)
		{
			if (val.array[i] < 1e-10)
			{
				val.array[i] = 1;
				zeroValue = i;
			}
		}

		U.set_mmul_NT(A, V);
		U.divT(val);

		U.makeT();
		int zerovec = -1;
		for (int i = 0; i < 3; i++)
		{
			if (U.vecs[i].norm2() < 1e-10)
			{
				zerovec = i;
			}
		}
		if (zerovec == 0)
		{
			U.vecs[0].set_cross(U.vecs[1], U.vecs[2]);
		}
		else if (zerovec == 1)
		{
			U.vecs[1].set_cross(U.vecs[2], U.vecs[0]);
		}
		else if (zerovec == 2)
		{
			U.vecs[2].set_cross(U.vecs[0], U.vecs[1]);
		}
		U.makeT();
		V.makeT();
		if(zeroValue != -1){
            val.array[zeroValue] = 0;
        }

		// B.set_mmul_NT(A, A);
		// B.eigenvec_and_eigenvals(val, U);
		// U.a.normalize();
		// U.b.normalize();
		// U.c.normalize();
		// // if(U.determinant() < 0)
		// // {
		// // 	int min_val = 0;
		// // 	if (val.x > val.y)
		// // 	{
		// // 		min_val = 1;
		// // 	}
		// // 	if (val.x > val.z)
		// // 	{
		// // 		min_val = 2;
		// // 	}
		// // 	if (val.z > val.y)
		// // 	{
		// // 		min_val = 1;
		// // 	}
		// // 	U.vecs[min_val].mul(-1);
		// // }
		//printf("U.det %f\n", U.determinant());

		// // double max_val = val.x;
		// // VEC tmp;
		// // if(val.x < val.y){ max_val = val.y; val.y = val.x; val.x = max_val; tmp = U.a; U.a = U.b; U.b = tmp; tmp = A.a; A.a = A.b; A.b = tmp; }
		// // if(val.y < val.z){ max_val = val.z; val.z = val.y; val.y = max_val; tmp = U.b; U.b = U.c; U.c = tmp; tmp = A.b; A.b = A.c; A.c = tmp; }
		// // if(val.x < val.y){ max_val = val.y; val.y = val.x; val.x = max_val; tmp = U.a; U.a = U.b; U.b = tmp; tmp = A.a; A.a = A.b; A.b = tmp; }
		// // printf("val %f %f %f\n", val.x, val.y, val.z);
		
		// val.x = sqrt(val.x);
		// val.y = sqrt(val.y);
		// val.z = sqrt(val.z);		
		// // V.a = A.dotT(U.a);
		// // V.b = A.dotT(U.b);
		// // V.c = A.dotT(U.c);
		// V.set_mmul_TN(A, U);
		// V.div(val);
		//printf("V.det %f\n", V.determinant());
		// // B.set_mmul_TN(A, A);
		// // B.eigenvec_and_eigenvals(val, V);
		// // V.a.normalize();
		// // V.b.normalize();
		// // V.c.normalize();
		// // max_val = val.x;
		// // if(val.x < val.y){ max_val = val.y; val.y = val.x; val.x = max_val; tmp = V.a; V.a = V.b; V.b = tmp; }
		// // if(val.y < val.z){ max_val = val.z; val.z = val.y; val.y = max_val; tmp = V.b; V.b = V.c; V.c = tmp; }
		// // if(val.x < val.y){ max_val = val.y; val.y = val.x; val.x = max_val; tmp = V.a; V.a = V.b; V.b = tmp; }


	}
};




/*
class Mat3i : public Mat3T< int   , Vec3i, Mat3i >{};
class Mat3f : public Mat3T< float , Vec3f, Mat3f >{};
class MAT : public Mat3T< T, VEC, MAT >{};
*/

using Mat3i = Mat3T< int   >;
using Mat3f = Mat3T< float >;
using Mat3d = Mat3T< double>;

static constexpr Mat3d Mat3dIdentity = Mat3d{1.0,0.0,0.0, 0.0,1.0,0.0,  0.0,0.0,1.0};
static constexpr Mat3d Mat3dZero     = Mat3d{0.0,0.0,0.0, 0.0,0.0,0.0,  0.0,0.0,0.0};

static constexpr Mat3f Mat3fIdentity = Mat3f{1.0f,0.0f,0.0f, 0.0f,1.0f,0.0f,  0.0f,0.0f,1.0f};
static constexpr Mat3f Mat3fZero     = Mat3f{0.0f,0.0f,0.0f, 0.0f,0.0f,0.0f,  0.0f,0.0f,0.0f};

inline void convert( const Mat3f& from, Mat3d& to ){ convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); };
inline void convert( const Mat3d& from, Mat3f& to ){ convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); };

inline Mat3f toFloat( const Mat3d& from){ Mat3f to; convert( from.a, to.a ); convert( from.b, to.b ); convert( from.c, to.c ); return to; }

template<typename T>
inline void wrapBondVec( Vec3T<T>& d, const Mat3T<T>& lvec, const Mat3T<T>& invLvec ){
    Vec3T<T> u;
    invLvec.dot_to( d, u );
    u.a=u.a+(1-(int)(u.a+1.5));
    u.b=u.b+(1-(int)(u.b+1.5));
    u.c=u.c+(1-(int)(u.c+1.5));
    lvec.dot_to_T( u, d );
}

#endif

