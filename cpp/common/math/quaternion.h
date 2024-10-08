

// read also:
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/quaternionToAngle/

#ifndef  quaternion_h
#define  quaternion_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"

#include "datatypes.h"

template <class T>
inline T project_beam_to_sphere( T r, T x, T y ){
	T z;
	T r2 = r * r;
	T d2 = x*x + y*y;
	if ( d2 < ( 0.5 * r2 ) ) {
		z = sqrt( r2 - d2 );
	} else {
		T t2 = 0.5 * r;
		z = sqrt( t2 / d2 );
	}
	return z;
}

/*
float project_beam_to_sphere( float r, float x, float y ){
	float d, t, z;
	d = sqrt( x*x + y*y );
	if ( d < r * 0.70710678118654752440 ) {
		z = sqrt( r*r - d*d );
	} else {
		t = r * 0.70710678118654752440; // == 1/sqrt(2)
		z = t*t / d;
	}
	return z;
}
*/

//template <class T, class VEC, class MAT, class QUAT>
//template <class T, class VEC>
template <class T>
class Quat4T {
	using VEC  = Vec3T<T>;
    using VEC2 = Vec2T<T>;
	using MAT  = Mat3T<T>;
	using QUAT = Quat4T<T>;
	public:
	union{
		struct{ T x,y,z,w; };
        struct{ VEC2 lo,hi;};
		struct{ VEC f;T e; }; // like force and energy
		struct{ VEC p;T s; }; // like molecular orbital basiset
		T array[4];
	};

    constexpr inline explicit operator Quat4T<double>()const{ return (Quat4T<double>){ (double)x, (double)y, (double)z, (double)w }; }
    constexpr inline explicit operator Quat4T<float> ()const{ return (Quat4T<float>) { (float )x, (float )y, (float )z, (float )w }; }
    constexpr inline explicit operator Quat4T<int>   ()const{ return (Quat4T<int>)   { (int   )x, (int   )y, (int   )z, (int   )w }; }

    constexpr inline explicit operator float4  ()const{ return float4 { (double)x, (double)y, (double)z, (double)w }; }
    constexpr inline explicit operator double4 ()const{ return double4{ (float )x, (float )y, (float )z, (float )w }; }
	constexpr inline explicit operator int4    ()const{ return int4   { (int   )x, (int   )y, (int   )z, (int   )w }; }

    //inline       T& operator[](int i){ return array[i]; }
    //inline const T& operator[](int i){ return array[i]; }

    constexpr inline void set   ( T f                     ){ x=f   ; y=f   ; z=f   ; w=f   ; }
	constexpr inline void set   ( const  QUAT& q          ){ x=q.x ; y=q.y ; z=q.z ; w=q.w ; }
	constexpr inline void set   ( T fx, T fy, T fz, T fw  ){ x=fx  ; y=fy  ; z=fz  ; w=fw  ; }
	constexpr inline void setXYZ( const VEC& v            ){ x=v.x ; y=v.y ; z=v.z ;         }
	constexpr inline void setOne(                         ){ x=y=z=0; w=1;                   }



	// swizzles
	constexpr inline VEC2 xy() const { return {x,y}; };
	constexpr inline VEC2 xz() const { return {x,z}; };
    constexpr inline VEC2 xw() const { return {x,w}; };

	constexpr inline VEC2 yx() const { return {y,x}; };
    constexpr inline VEC2 yz() const { return {y,z}; };
    constexpr inline VEC2 yw() const { return {y,w}; };

	constexpr inline VEC2 zx() const { return {z,x}; };
    constexpr inline VEC2 zy() const { return {z,y}; };
    constexpr inline VEC2 zw() const { return {z,w}; };

    inline VEC2 wx() const { return {w,x}; };
    inline VEC2 wy() const { return {w,y}; };
    inline VEC2 wz() const { return {w,z}; };




    constexpr inline void setView_XY (){set(1.0,0.0,0.0,0.0);}    // top view  x=x, y=-y,
    constexpr inline void setView_XYm(){set(0.0,0.0,0.0,1.0);}  // bottom view
    //inline void set( 0.0,0.0,1.0,0.0 );  // bottom view
    //inline void set( 0.0,1.0,0.0,0.0 );  // top view  x=-x, y=y,
    //inline void setView_(){set( +0.70710678118,0.0,0.0,  0.70710678118 );} // side down
    //inline void setView_(){set( -0.70710678118,0.0,0.0,  0.70710678118 );} // x=x, z=y,  y=-y,
    //inline void setView_(){set( 0.0, -0.70710678118,0.0, 0.70710678118 );} // z=-x, y=y
    //inline void setView_(){set( 0.0, +0.70710678118,0.0, 0.70710678118 );} // z=+x, y=y
    //inline void setView_(){set( 0.0,0.0, +0.70710678118, 0.70710678118 );} // y=-x, x=y
    //inline void setView_(){set( 0.0,0.0, -0.70710678118, 0.70710678118 );} // y=x, x=-y

	constexpr inline void setInverseUnitary( const  QUAT& q){ x=-q.x; y=-q.y; z=-q.z; w=q.w; };
	constexpr inline void setInverse       ( const  QUAT& q){ setInverseUnitary(); mul(1./q.norm2()); };

	constexpr inline QUAT get_inv(){ QUAT q; q.x=-x; q.y=-y; q.z=-z; q.w=w; return q; }

    constexpr inline bool allEqual(T f)const{ return (x==f)&&(y==f)&&(z==f)&&(w==f);};
    constexpr inline bool anyEqual(T f)const{ return (x==f)||(y==f)||(z==f)||(w==f);};

// ====== basic aritmetic

// ============== Basic Math

    constexpr inline void mul        ( T f                                  ){ x*=f;            y*=f;           z*=f;           w*=f;            };
    constexpr inline void add        ( const QUAT& v                        ){ x+=v.x;          y+=v.y;         z+=v.z;         w+=v.w;          };
    constexpr inline void sub        ( const QUAT& v                        ){ x-=v.x;          y-=v.y;         z-=v.z;         w-=v.w;          };
	constexpr inline void set_add    ( const QUAT& a, const QUAT& b         ){ x =a.x+b.x;      y =a.y+b.y;     z =a.z+b.z;     w =a.w+b.w;      };
	constexpr inline void set_sub    ( const QUAT& a, const QUAT& b         ){ x =a.x-b.x;      y =a.y-b.y;     z =a.z-b.z;     w =a.w-b.w;      };
	constexpr inline void set_mul    ( const QUAT& a, T f                   ){ x =a.x*f;        y =a.y*f;       z =a.z*f;       w =a.w*f;        };
	constexpr inline void set_mul    ( const QUAT& a, const QUAT& b         ){ x =a.x*b.x;      y =a.y*b.y;     z =a.z*b.z;     w =a.w*b.w;      };
	constexpr inline void add_mul    ( const QUAT& a, T f                ){ x+=a.x*f;        y+=a.y*f;       z+=a.z*f;       w+=a.w*f;        };
	constexpr inline void set_add_mul( const QUAT& a, const QUAT& b, T f ){ x =a.x + f*b.x;  y =a.y + f*b.y; z =a.z + f*b.z; w =a.w + f*b.w;  };

    constexpr inline T dot  ( QUAT q   ) const {  return       w*q.w + x*q.x + y*q.y + z*q.z;   }
    constexpr inline T norm2(          ) const {  return       w*  w + x*  x + y*  y + z*  z;   }
	constexpr inline T norm (          ) const {  return sqrt( w*  w + x*  x + y*  y + z*  z ); }
    constexpr inline T normalize() {
		T norm  = sqrt( x*x + y*y + z*z + w*w );
		T inorm = 1/norm;
		x *= inorm;    y *= inorm;    z *= inorm;   w *= inorm;
		return norm;
    }
    constexpr inline void normalizeW(){
        double r2 = x*x + y*y + z*z;
        if(r2>1){ printf( "ERROR Quat4<T>::normalizeW() r2>1\n"); exit(0); };
        w = sqrt( 1-r2 );
    }

    constexpr inline void checkNormalized( T D2 ){
        T r2 =  x*x + y*y + z*z + w*w;
        //printf( " (%g,%g,%g,%g) r2 %g \n", x,y,z,w,   r2 );
        T d2 = r2 - 1;
        if( (d2>D2) || (d2<-D2) ){ //printf( "renorm\n" );
            T inorm = 1/sqrt( r2 );
            x *= inorm;    y *= inorm;    z *= inorm;   w *= inorm;
        }
    }

    constexpr inline T normalize_taylor3(){
        // sqrt(1+x) ~= 1 + 0.5*x - 0.125*x*x
        // sqrt(r2) = sqrt((r2-1)+1) ~= 1 + 0.5*(r2-1)
        // 1/sqrt(1+x) ~= 1 - 0.5*x + (3/8)*x^2 - (5/16)*x^3 + (35/128)*x^4 - (63/256)*x^5
        T dr2    = x*x+y*y+z*z+w*w-1;
        T invr = 1 + dr2*( -0.5 + dr2*( 0.375 + dr2*-0.3125 ) );
        x*=invr;
        y*=invr;
        z*=invr;
        w*=invr;
        return dr2;
    }

    constexpr inline T normalize_hybrid(T dr4max){
        T dr2    = x*x+y*y+z*z+w*w-1;
        T dr4    = dr2*dr2;
        T invr;
        if(dr4<dr4max){ invr = 1 + dr2*-0.5 + dr4*( 0.375 + dr2*-0.3125 ); }
        else          { invr = 1/sqrt(dr2+1); };
        x*=invr;
        y*=invr;
        z*=invr;
        w*=invr;
        return dr2;
    }

    constexpr inline T makeOrthoU(const QUAT& q){ T c = dot(q);           add_mul(q,-c); return c; }
    constexpr inline T makeOrtho (const QUAT& q){ T c = dot(q)/q.norm2(); add_mul(q,-c); return c; }

// ====== Quaternion multiplication


	constexpr inline void setQmul( const QUAT& a, const QUAT& b) {
        x =  a.x * b.w + a.y * b.z - a.z * b.y + a.w * b.x;
        y = -a.x * b.z + a.y * b.w + a.z * b.x + a.w * b.y;
        z =  a.x * b.y - a.y * b.x + a.z * b.w + a.w * b.z;
        w = -a.x * b.x - a.y * b.y - a.z * b.z + a.w * b.w;
    };

    constexpr inline void qmul( const QUAT& a) {
        T aw = a.w, ax = a.x, ay = a.y, az = a.z;
        T x_ =  x * aw + y * az - z * ay + w * ax;
        T y_ = -x * az + y * aw + z * ax + w * ay;
        T z_ =  x * ay - y * ax + z * aw + w * az;
                w = -x * ax - y * ay - z * az + w * aw;
        x = x_; y = y_; z = z_;
    };

    constexpr inline void qmul_T( const QUAT& a) {
        T aw = a.w, ax = a.x, ay = a.y, az = a.z;
        T x_ =  ax * w + ay * z - az * y + aw * x;
        T y_ = -ax * z + ay * w + az * x + aw * y;
        T z_ =  ax * y - ay * x + az * w + aw * z;
              w = -ax * x - ay * y - az * z + aw * w;
        x = x_; y = y_; z = z_;
    };

    constexpr inline void qmul_it( QUAT& a) const {
        T aw = a.w, ax = a.x, ay = a.y, az = a.z;
        a.x =  x * aw + y * az - z * ay + w * ax;
        a.y = -x * az + y * aw + z * ax + w * ay;
        a.z =  x * ay - y * ax + z * aw + w * az;
        a.w = -x * ax - y * ay - z * az + w * aw;
    };

    constexpr inline void qmul_it_T( QUAT& a) const {
        T aw = a.w, ax = a.x, ay = a.y, az = a.z;
        a.x = -x * aw - y * az + z * ay + w * ax;
        a.y = +x * az - y * aw - z * ax + w * ay;
        a.z = -x * ay + y * ax - z * aw + w * az;
        a.w = +x * ax + y * ay + z * az + w * aw;
    };

    constexpr inline QUAT operator%( const QUAT& b)const{ 
        return QUAT{
             x * b.w + y * b.z - z * b.y + w * b.x,
            -x * b.z + y * b.w + z * b.x + w * b.y,
             x * b.y - y * b.x + z * b.w + w * b.z,
            -x * b.x - y * b.y - z * b.z + w * b.w
        };
    }

    constexpr inline VEC operator%( const VEC& v)const{  // Hamilton product ( q * v ) vector is rotated by quaternion
        T tx = 2*( z*v.y - y*v.z  + w*v.x);
        T ty = 2*( x*v.z - z*v.x  + w*v.y);
        T tz = 2*( y*v.x - x*v.y  + w*v.z);
        return VEC{
            v.x + ty*z - tz*y,
            v.y + tz*x - tx*z,
            v.z + tx*y - ty*x
        };
    }

    constexpr inline VEC operator|( const VEC& v )const{  // inverse Hamilton product ( v * q ) vector is rotated by quaternion backwards
        T tx = 2*( y*v.z - z*v.y  + w*v.x);
        T ty = 2*( z*v.x - x*v.z  + w*v.y);
        T tz = 2*( x*v.y - y*v.x  + w*v.z);
        return VEC{
            v.x + y*tz - z*ty,
            v.y + z*tx - x*tz,
            v.z + x*ty - y*tx
        };
    }

    constexpr inline void transformVec( const VEC& vec, VEC& out ) const{
        //QUAT qv; qv.set(vec.x,vec.y,vec.z,0.0);
        //qmul_it  (qv);
        //qmul_it_T(qv);
        //out.set(qv.x,qv.y,qv.z);
        // https://blog.molecular-matters.com/2013/05/24/a-faster-quaternion-vector-multiplication/

        // http://stackoverflow.com/questions/22497093/faster-quaternion-vector-multiplication-doesnt-work
        //t = 2 * cross(q.xyz, v); v' = v + q.w * t + cross(q.xyz, t)
        //
        //t = 2 * cross(q.xyz, v)
        //T tx = 2*( y*vec.z - z*vec.y );
        //T ty = 2*( z*vec.x - x*vec.z );
        //T tz = 2*( x*vec.y - y*vec.x );
        // v' = v + q.w * t + cross(q.xyz, t)
        //out.x   = vec.x + w*tx + y*tz - z*ty;
        //out.y   = vec.y + w*ty + z*tx - x*tz;
        //out.z   = vec.z + w*tz + x*ty - y*tx;

        // v' = v + 2.0 * cross(cross(v, q.xyz) + q.w * v, q.xyz);
        T tx = 2*( vec.y*z - vec.z*y  + w*vec.x);
        T ty = 2*( vec.z*x - vec.x*z  + w*vec.y);
        T tz = 2*( vec.x*y - vec.y*x  + w*vec.z);
        out.x   = vec.x + ty*z - tz*y;
        out.y   = vec.y + tz*x - tx*z;
        out.z   = vec.z + tx*y - ty*x;
        // 15 mult, 12 add
    }

    constexpr inline void untransformVec( const VEC& vec, VEC& out ) const{
        T tx = 2*( y*vec.z - z*vec.y  + w*vec.x);
        T ty = 2*( z*vec.x - x*vec.z  + w*vec.y);
        T tz = 2*( x*vec.y - y*vec.x  + w*vec.z);
        out.x   = vec.x + y*tz - z*ty;
        out.y   = vec.y + z*tx - x*tz;
        out.z   = vec.z + x*ty - y*tx;
    }

    constexpr inline void invertUnitary() { x=-x; y=-y; z=-z; }

    constexpr inline void invert() {
		T norm = sqrt( x*x + y*y + z*z + w*w );
		if ( norm > 0.0 ) {
			T invNorm = 1/norm;
			x *= -invNorm; y *= -invNorm;z *= -invNorm;	w *=  invNorm;
		}
    };

    //inline QUAT operator+ ( T f   ) const { return QUAT{ x+f, y+f, z+f, w+f }; };
    constexpr inline QUAT operator* ( T f   ) const { return QUAT{ x*f, y*f, z*f, w*f }; };

    constexpr inline QUAT operator+ ( const QUAT& vi ) const { return QUAT{ x+vi.x, y+vi.y, z+vi.z, w+vi.w }; };
    constexpr inline QUAT operator- ( const QUAT& vi ) const { return QUAT{ x-vi.x, y-vi.y, z-vi.z, w-vi.w }; };
    constexpr inline QUAT operator* ( const QUAT& vi ) const { return QUAT{ x*vi.x, y*vi.y, z*vi.z, w*vi.w }; };
    constexpr inline QUAT operator/ ( const QUAT& vi ) const { return QUAT{ x/vi.x, y/vi.y, z/vi.z, w/vi.w }; };

// ======= metric
    // https://fgiesen.wordpress.com/2013/01/07/small-note-on-quaternion-distance-metrics/
    // metric on quaternions : http://www.cs.cmu.edu/~cga/dynopt/readings/Rmetric.pdf = /home/prokop/Dropbox/KnowDev/quaternions/Rotation_metric_Rmetric.pdf
    //  q and -q denote the same rotation !!!
    //  rho(q,q0)   = |q - q0|
    //  rho(q,q0)   = arccos( dot(q,q) )  ~=   1 - dot(q,q)

    constexpr inline double dist_cos( const QUAT& q0 ) const {
        double cdot = dot( q0 );
        return 1-((cdot>=0)?cdot:-cdot);    // consider q=-q
    }

    constexpr inline double ddist_cos( const QUAT& q0, QUAT& dRdq ) const {
        dRdq.set(q0);
        double cdot = dot( q0 );
        if( cdot<0 ){ dRdq.mul(-1); cdot=-cdot; };  // consider q=-q
        return 1-cdot;
    }

    constexpr inline double sub_paralel_fast( const QUAT& q ) {
        // substract component of *this paralel to q assuming that q is normalized
        double cdot = dot( q );
        add_mul( q, -cdot );
		return cdot;
    }


// ======= Conversion : Angle & Axis

    constexpr inline void fromUniformS3(VEC u){ // u is vec of 3 random numbers from (0..1)
        //  http://planning.cs.uiuc.edu/node198.html
        T a=sqrt(1-u.x);
        T b=sqrt(  u.x);
        u.y *= 2*M_PI;
        u.z *= 2*M_PI;
        x = a*sin(u.y);
        y = a*cos(u.y);
        z = b*sin(u.z);
        w = b*cos(u.z);
    }

    constexpr inline void setRandomRotation(){ fromUniformS3( {randf(), randf(), randf()} ); }


	constexpr inline void fromAngleAxis( T angle, const VEC& axis ){
		T ir   = 1/axis.norm();
		VEC  hat  = axis * ir;
		T a    = 0.5 * angle;
		T sa   = sin(a);
		w =           cos(a);
		x = sa * hat.x;
		y = sa * hat.y;
		z = sa * hat.z;
	};

	constexpr void fromCosAngleAxis( T scal_prod, const VEC& axis ){
	// we assume -phi instead of phi!!!, minus effectively implies sa -> -s
		constexpr T cos_cutoff = 1 - 1e-6;
		T cosphi, sinphi, sa, phi, sgn_sinphi;
		T ir    = 1.0 / axis.norm();
		VEC  hat   = axis * ir;
		cosphi     = scal_prod;
		sgn_sinphi = 1.0; // ?
		if( cosphi > cos_cutoff ){
			sa = 0; w = 1;
		} else if( cosphi < -( cos_cutoff ) ){
			sa = -1; w = 0;
		} else {
			sa = + sqrt( (1 - cosphi) / 2.0 );
			w  = - sqrt( (1 + cosphi) / 2.0 ) * sgn_sinphi;
//			sa = -sa; w = -w;
		}
		x = sa * hat.x;
		y = sa * hat.y;
		z = sa * hat.z;

	}

	#define TRACKBALLSIZE ( 0.8 )

	/*
	void fromTrackball( T p1x, T p1y, T p2x, T p2y ){
		VEC  axis; // axis of rotation
		//T phi;  // angle of rotation
		VEC  p1, p2, d;
		//T t;
		//if( ( sq(p2x-p1x)+sq(p2y-p1y) ) < 1e-8 ){   }
		if( ( p2x == p1x ) && ( p2y == p1y ) ){ setOne(); return; }
		p1.set( p1x, p1y, project_beam_to_sphere<T>( TRACKBALLSIZE, p1x, p1y ) );
		p2.set( p2x, p2y, project_beam_to_sphere<T>( TRACKBALLSIZE, p2x, p2y ) );
		axis.set_cross( p2, p1 );


		T t = d.norm() / ( 2.0 * TRACKBALLSIZE );
		if( t > 1.0 )  t =  1.0;
		if( t < -1.0 ) t = -1.0;
        T phi = 2.0 * asin( t );
        fromAngleAxis( phi, axis );


        //T t = sqrt( 1 - d.norm2() );
        //fromCosAngleAxis( t, axis );

		// T cosphi = ;
		// phi = 2.0 * asin( t );
		// fromCosAngleAxis( cosphi, axis );

	}
	*/

	constexpr void fromTrackball( T p1x, T p1y, T p2x, T p2y ){
		VEC  axis; // axis of rotation
		//T phi;  // angle of rotation
		VEC  p1, p2, d;
		//T t;
		//if( p1x == p2x && p1y == p2y ){	setOne(); return; }
		if( ( sq<T>(p2x-p1x)+sq<T>(p2y-p1y) ) < 1e-8 ){ setOne(); return; }
		p1.set( p1x, p1y, project_beam_to_sphere<T>( TRACKBALLSIZE, p1x, p1y ) );
		p2.set( p2x, p2y, project_beam_to_sphere<T>( TRACKBALLSIZE, p2x, p2y ) );
		axis.set_cross( p2, p1 );
		d.set_sub( p1, p2 );

        /*
		T t = d.norm() / ( 2.0 * TRACKBALLSIZE );
		if( t > 1.0 )  t =  1.0;
		if( t < -1.0 ) t = -1.0;
		T phi = 2.0 * asin( t );
		fromAngleAxis( phi, axis );
		*/

        T t = sqrt( 1 - d.norm2() );
        fromCosAngleAxis( t, axis );

	}


/*
    void fromTrackball_0( T px, T py ){
		VEC  axis; // axis of rotation
		T phi;  // angle of rotation
		VEC  p1, p2, d;
		T t;
		if( ( px == 0 ) && ( p == py ) ){
			setOne();
			return;
		}
		p2.set( px, py, project_beam_to_sphere<T>( TRACKBALLSIZE, px, py ) );
		p1.set( 0, 0, TRACKBALLSIZE );
		axis.set_cross( p2, p1 );
		d.set_sub( p1, p2 );
		t = d.norm() / ( 2.0 * TRACKBALLSIZE );
		if( t > 1.0 )  t =  1.0;
		if( t < -1.0 ) t = -1.0;

		phi = 2.0 * asin( t );
		fromAngleAxis( phi, axis );

	}
*/

// =======  pitch, yaw, roll

	constexpr inline void dpitch( T angle ){ T ca,sa; sincos_taylor2(angle*(T)0.5,sa,ca); pitch( ca, sa );  };
	constexpr inline void pitch ( T angle ){ angle*=(T)0.5; pitch( cos(angle), sin(angle) );  };
    constexpr inline void pitch ( T ca, T sa ) {
        T x_ =  x * ca + w * sa;
        T y_ =  y * ca + z * sa;
        T z_ = -y * sa + z * ca;
             w = -x * sa + w * ca;
        x = x_; y = y_; z = z_;
    };

	constexpr inline void dyaw( T angle ){ T ca,sa; sincos_taylor2(angle*(T)0.5,sa,ca); yaw( ca, sa );  };
	constexpr inline void yaw ( T angle ){ angle*=(T)0.5; yaw( cos(angle), sin(angle) );  };
    constexpr inline void yaw ( T ca, T sa ) {
        T x_ =  x * ca - z * sa;
        T y_ =  y * ca + w * sa;
        T z_ =  x * sa + z * ca;
             w = -y * sa + w * ca;
        x = x_; y = y_; z = z_;
    };

	constexpr inline void droll( T angle ){ T ca,sa; sincos_taylor2(angle*(T)0.5,sa,ca); roll( ca, sa );  };
	constexpr inline void roll ( T angle ){ angle*=(T)0.5; roll( cos(angle), sin(angle) );  };
    constexpr inline void roll ( T ca, T sa ) {
        T x_ =  x * ca + y * sa;
        T y_ = -x * sa + y * ca;
        T z_ =  z * ca + w * sa;
             w = -z * sa + w * ca;
        x = x_; y = y_; z = z_;
    };


	constexpr inline void dpitch2( T angle ){ T ca,sa; sincos_taylor2(angle*(T)0.5,sa,ca); pitch2( ca, sa );  };
	constexpr inline void pitch2 ( T angle ){ angle*=(T)0.5; pitch2( cos(angle), sin(angle) );  };
    constexpr inline void pitch2 ( T ca, T sa ) {
        T x_ =  sa * w + ca * x;
        T y_ = -sa * z + ca * y;
        T z_ =  sa * y + ca * z;
             w  = -sa * x + ca * w;
        x = x_; y = y_; z = z_;
    };

	constexpr inline void dyaw2( T angle ){ T ca,sa; sincos_taylor2(angle*(T)0.5,sa,ca); yaw2( ca, sa );  };
	constexpr inline void yaw2 ( T angle ){ angle*=(T)0.5; yaw2( cos(angle), sin(angle) );  };
    constexpr inline void yaw2 ( T ca, T sa ) {
        T x_ = + sa * z  + ca * x;
        T y_ = + sa * w  + ca * y;
        T z_ = - sa * x  + ca * z;
             w  = - sa * y  + ca * w;
        x = x_; y = y_; z = z_;
    };

	constexpr inline void droll2( T angle ){ T ca,sa; sincos_taylor2(angle*(T)0.5,sa,ca); roll2( ca, sa );  };
	constexpr inline void roll2 ( T angle ){ angle*=(T)0.5; roll2( cos(angle), sin(angle) );  };
    constexpr inline void roll2 ( T ca, T sa ) {
        //ca *=0.5; sa *=0.5; // seems that should be just half
        T x_ = - sa * y + ca * x;
        T y_ = + sa * x + ca * y;
        T z_ = + sa * w + ca * z;
             w  = - sa * z + ca * w;
        x = x_; y = y_; z = z_;
    };

// ====== Differential rotation

	constexpr inline void dRot_exact ( T dt, const VEC& omega ) {
		T hx   = omega.x;
		T hy   = omega.y;
		T hz   = omega.z;
		T r2   = hx*hx + hy*hy + hz*hz;
		if(r2>0){
			T norm = sqrt( r2 );
			T a    = dt * norm * 0.5;
			T sa   = sin( a )/norm;  // we normalize it here to save multiplications
			T ca   = cos( a );
			hx*=sa; hy*=sa; hz*=sa;            // hat * sin(a)
			T x_ = x, y_ = y, z_ = z, w_ = w;
			x =  hx*w_ + hy*z_ - hz*y_ + ca*x_;
			y = -hx*z_ + hy*w_ + hz*x_ + ca*y_;
			z =  hx*y_ - hy*x_ + hz*w_ + ca*z_;
			w = -hx*x_ - hy*y_ - hz*z_ + ca*w_;
		}
	};


	constexpr inline void dRot_taylor2 ( T dt, VEC& omega ) {
		T hx   = omega.x;
		T hy   = omega.y;
		T hz   = omega.z;
		T r2   = hx*hx + hy*hy + hz*hz;
		T b2   = dt*dt*r2;
		const T c2 = 1.0/8;    // 4  *2
		const T c3 = 1.0/48;   // 8  *2*3
		const T c4 = 1.0/384;  // 16 *2*3*4
		const T c5 = 1.0/3840; // 32 *2*3*4*5
		T sa   = dt * ( 0.5 - b2*( c3 - c5*b2 ) );
		T ca   =      ( 1   - b2*( c2 - c4*b2 ) );
		hx*=sa; hy*=sa; hz*=sa;  // hat * sin(a)
		T x_ = x, y_ = y, z_ = z, w_ = w;
		x =  hx*w_ + hy*z_ - hz*y_ + ca*x_;
		y = -hx*z_ + hy*w_ + hz*x_ + ca*y_;
		z =  hx*y_ - hy*x_ + hz*w_ + ca*z_;
		w = -hx*x_ - hy*y_ - hz*z_ + ca*w_;
	};


	constexpr inline void toMatrix( MAT& result) const {
		    T r2 = w*w + x*x + y*y + z*z;
		    //T s  = (r2 > 0) ? 2d / r2 : 0;
			T s  = 2 / r2;
		    // compute xs/ys/zs first to save 6 multiplications, since xs/ys/zs
		    // will be used 2-4 times each.
		    T xs = x * s;  T ys = y * s;  T zs = z * s;
		    T xx = x * xs; T xy = x * ys; T xz = x * zs;
		    T xw = w * xs; T yy = y * ys; T yz = y * zs;
		    T yw = w * ys; T zz = z * zs; T zw = w * zs;
		    // using s=2/norm (instead of 1/norm) saves 9 multiplications by 2 here
		    result.xx = 1 - (yy + zz);
		    result.xy =     (xy - zw);
		    result.xz =     (xz + yw);
		    result.yx =     (xy + zw);
		    result.yy = 1 - (xx + zz);
		    result.yz =     (yz - xw);
		    result.zx =     (xz - yw);
		    result.zy =     (yz + xw);
		    result.zz = 1 - (xx + yy);
	};


    constexpr inline void toMatrix_T( MAT& result) const {
		    T r2 = w*w + x*x + y*y + z*z;
		    //T s  = (r2 > 0) ? 2d / r2 : 0;
			T s  = 2 / r2;
		    // compute xs/ys/zs first to save 6 multiplications, since xs/ys/zs
		    // will be used 2-4 times each.
		    T xs = x * s;  T ys = y * s;  T zs = z * s;
		    T xx = x * xs; T xy = x * ys; T xz = x * zs;
		    T xw = w * xs; T yy = y * ys; T yz = y * zs;
		    T yw = w * ys; T zz = z * zs; T zw = w * zs;
		    // using s=2/norm (instead of 1/norm) saves 9 multiplications by 2 here
		    result.xx = 1 - (yy + zz);
		    result.yx =     (xy - zw);
		    result.zx =     (xz + yw);
		    result.xy =     (xy + zw);
		    result.yy = 1 - (xx + zz);
		    result.zy =     (yz - xw);
		    result.xz =     (xz - yw);
		    result.yz =     (yz + xw);
		    result.zz = 1 - (xx + yy);
	};

/*
	inline void toMatrix_unitary( MAT& result)  const  {
		T xx = x * x;
		T xy = x * y;
		T xz = x * z;
		T xw = x * w;
		T yy = y * y;
		T yz = y * z;
		T yw = y * w;
		T zz = z * z;
		T zw = z * w;
		result.xx = 1 - 2 * ( yy + zz );
		result.xy =     2 * ( xy - zw );
		result.xz =     2 * ( xz + yw );
		result.yx =     2 * ( xy + zw );
		result.yy = 1 - 2 * ( xx + zz );
		result.yz =     2 * ( yz - xw );
		result.zx =     2 * ( xz - yw );
		result.zy =     2 * ( yz + xw );
		result.zz = 1 - 2 * ( xx + yy );
	};
*/

	constexpr inline void toMatrix_unitary( MAT& result)  const  {
		T x2 = 2*x;
		T y2 = 2*y;
		T z2 = 2*z;
		T xx = x2 * x;
		T xy = x2 * y;
		T xz = x2 * z;
		T xw = x2 * w;
		T yy = y2 * y;
		T yz = y2 * z;
		T yw = y2 * w;
		T zz = z2 * z;
		T zw = z2 * w;
		result.xx = 1 - ( yy + zz );
		result.xy =     ( xy - zw );
		result.xz =     ( xz + yw );
		result.yx =     ( xy + zw );
		result.yy = 1 - ( xx + zz );
		result.yz =     ( yz - xw );
		result.zx =     ( xz - yw );
		result.zy =     ( yz + xw );
		result.zz = 1 - ( xx + yy );
	};

    constexpr inline void toMatrix_unitary_T( MAT& result)  const  {
		T x2 = 2*x;
		T y2 = 2*y;
		T z2 = 2*z;
		T xx = x2 * x;
		T xy = x2 * y;
		T xz = x2 * z;
		T xw = x2 * w;
		T yy = y2 * y;
		T yz = y2 * z;
		T yw = y2 * w;
		T zz = z2 * z;
		T zw = z2 * w;
		result.xx = 1 - ( yy + zz );
		result.yx =     ( xy - zw );
		result.zx =     ( xz + yw );
		result.xy =     ( xy + zw );
		result.yy = 1 - ( xx + zz );
		result.zy =     ( yz - xw );
		result.xz =     ( xz - yw );
		result.yz =     ( yz + xw );
		result.zz = 1 - ( xx + yy );
	};

    // This allos passing Quad to functions accepting Mat3f (e.g. to plotting functions)
    //inline explicit operator Mat3T<T>()const{ Mat3T<T> mat; toMatrix_unitary(mat); return mat; }
    constexpr inline Mat3T<T> toMat () const {  Mat3T<T> mat; toMatrix_unitary  (mat); return mat; };
    constexpr inline Mat3T<T> toMatT() const {  Mat3T<T> mat; toMatrix_unitary_T(mat); return mat; };

	// this will compute force on quaternion from force on some point "p" in coordinate system of the quaternion
	//   EXAMPLE : if "p" is atom in molecule, it should be local coordinate in molecular local space, not global coordante after quaternion rotation is applied
	constexpr inline void addForceFromPoint( const VEC& p, const VEC& fp, QUAT& fq ) const {
		// dE/dx = dE/dpx * dpx/dx
		// dE/dx = fx
		T px_x =    p.b*y +             p.c*z;
		T py_x =    p.a*y - 2*p.b*x +   p.c*w;
		T pz_x =    p.a*z -   p.b*w - 2*p.c*x;
		fq.x += 2*( fp.x * px_x  +  fp.y * py_x  + fp.z * pz_x );

		T px_y = -2*p.a*y +   p.b*x -   p.c*w;
		T py_y =    p.a*x +             p.c*z;
		T pz_y =    p.a*w +   p.b*z - 2*p.c*y;
		fq.y += 2*( fp.x * px_y  +  fp.y * py_y  + fp.z * pz_y );

		T px_z = -2*p.a*z +   p.b*w +   p.c*x;
		T py_z =   -p.a*w - 2*p.b*z +   p.c*y;
		T pz_z =    p.a*x +   p.b*y;
		fq.z += 2*( fp.x * px_z  +  fp.y * py_z  + fp.z * pz_z );

		T px_w =    p.b*z -   p.c*y;
		T py_w =   -p.a*z +   p.c*x;
		T pz_w =    p.a*y -   p.b*x;
		fq.w += 2*( fp.x * px_w  +  fp.y * py_w  + fp.z * pz_w );

	}

    //T makeOrthoU( const VEC& a ){ T c = dot(a);          add_mul(a, -c); return c; }
	constexpr inline T outproject( const QUAT& q ){ T cdot = dot(q); add_mul( q, -cdot ); return cdot; };

	constexpr inline void fromMatrix ( const VEC& a, const VEC& b, const VEC& c ) { fromMatrix( a.x,  a.y,  a.z,  b.x,  b.y,  b.z,  c.x,  c.y,  c.z  );  }
	constexpr inline void fromMatrix ( const MAT& M                             ) { fromMatrix( M.ax, M.ay, M.az, M.bx, M.by, M.bz, M.cx, M.cy, M.cz );  }
	constexpr inline void fromMatrixT( const MAT& M                             ) { fromMatrix( M.ax, M.bx, M.cx, M.ay, M.by, M.cy, M.az, M.bz, M.cz );  }
	constexpr inline void fromMatrix ( T m00, T m01, T m02,    T m10, T m11, T m12,        T m20, T m21, T m22) {
        // Use the Graphics Gems code, from
        // ftp://ftp.cis.upenn.edu/pub/graphics/shoemake/quatut.ps.Z
        T t = m00 + m11 + m22;
        // we protect the division by s by ensuring that s>=1
        if (t >= 0) { // by w
            T s = sqrt(t + 1);
            w = 0.5 * s;
            s = 0.5 / s;
            x = (m21 - m12) * s;
            y = (m02 - m20) * s;
            z = (m10 - m01) * s;
        } else if ((m00 > m11) && (m00 > m22)) { // by x
            T s = sqrt(1 + m00 - m11 - m22);
            x = s * 0.5;
            s = 0.5 / s;
            y = (m10 + m01) * s;
            z = (m02 + m20) * s;
            w = (m21 - m12) * s;
        } else if (m11 > m22) { // by y
            T s = sqrt(1 + m11 - m00 - m22);
            y = s * 0.5;
            s = 0.5 / s;
            x = (m10 + m01) * s;
            z = (m21 + m12) * s;
            w = (m02 - m20) * s;
        } else { // by z
            T s = sqrt(1 + m22 - m00 - m11);
            z = s * 0.5;
            s = 0.5 / s;
            x = (m02 + m20) * s;
            y = (m21 + m12) * s;
            w = (m10 - m01) * s;
        }
	}

    //template<typename T>
    constexpr inline Mat3T<T> rotateVectors( int n, Vec3T<T>* v0s, Vec3T<T>* vs, const bool bT )const{
        Mat3T<T> mrot;
        if(bT){ toMatrix_T(mrot); }else{ toMatrix(mrot); };
        //if(bool bT){ toMatrix_unitary_T(mrot); }else{ toMatrix_unitary(mrot); };
        mrot.transformVectors(n,v0s,vs);
        return mrot;
    }

    //template<typename T>
    constexpr inline Mat3T<T> rotatePoints0( int n, Vec3T<T>* p0s, Vec3T<T>* ps, const Vec3T<T>& toPos, const bool bT )const{
        Mat3T<T> mrot;
        if(bT){ toMatrix_T(mrot); }else{ toMatrix(mrot); };
        //if(bool bT){ toMatrix_unitary_T(mrot); }else{ toMatrix_unitary(mrot); };
        mrot.transformPoints0(n,p0s,ps,toPos);
        return mrot;
    }

    constexpr inline Mat3T<T> rotatePoints( int n, Vec3T<T>* p0s, Vec3T<T>* ps, const Vec3T<T>& pos0, const bool bT )const{
        Mat3T<T> mrot;
        if(bT){ toMatrix_T(mrot); }else{ toMatrix(mrot); };
        //if(bool bT){ toMatrix_unitary_T(mrot); }else{ toMatrix_unitary(mrot); };
        mrot.transformPoints(n,p0s,ps,pos0);
        return mrot;
    }

    constexpr inline QUAT& setIfLower  (const QUAT& a){ if(a.x<x)x=a.x;if(a.y<y)y=a.y;if(a.z<z)z=a.z;if(a.w<w)w=a.w; return *this; }
    constexpr inline QUAT& setIfGreater(const QUAT& a){ if(a.x>x)x=a.x;if(a.y>y)y=a.y;if(a.z>z)z=a.z;if(a.w>w)w=a.w; return *this; }
    constexpr inline void update_bounds(QUAT& pmin,QUAT& pmax)const{ 
        if ( x < pmin.x ){ pmin.x=x; } else if ( x > pmax.x ){ pmax.x=x; };
        if ( y < pmin.y ){ pmin.y=y; } else if ( y > pmax.y ){ pmax.y=y; };
        if ( z < pmin.z ){ pmin.z=z; } else if ( z > pmax.z ){ pmax.z=z; };
        if ( z < pmin.w ){ pmin.w=w; } else if ( w > pmax.w ){ pmax.w=w; };
    }

};

/*
class Quat4i : public Quat4T< int,    Vec3i, Mat3i, Quat4i >{  };
class Quat4f : public Quat4T< float,  Vec3f, Mat3f, Quat4f >{  };
class QUAT : public Quat4T< T, VEC, MAT, QUAT >{  };
*/

template<typename T1,typename T2>
constexpr inline Quat4T<T2> cast(const Quat4T<T1>& i){ return Quat4T<T2>{ (T2)i.x, (T2)i.y, (T2)i.z, (T2)i.w }; };

//inline float4 cast(const Quat4T<T1>& i){ return float4{ (float)i.x; (float)i.y; (float)i.z; (float)i.w }; };
//inline float4 cast(const Quat4T<T1>& i){ Quat4T<T2> o; o.x=(T2)i.x; o.y=(T2)i.y; o.z=(T2)i.z; o.w=(T2)i.w; return o; };

using Quat4i = Quat4T< int>;
using Quat4f = Quat4T< float>;
using Quat4d = Quat4T< double >;

static constexpr Quat4i Quat4iZero      {0,0,0,0};
static constexpr Quat4i Quat4iOnes      {1,1,1,1};
static constexpr Quat4i Quat4iMinusOnes {-1,-1,-1,-1};

static constexpr Quat4d Quat4dNAN  {NAN,NAN,NAN,NAN};
static constexpr Quat4d Quat4dZero {0.0,0.0,0.0,0.0};
static constexpr Quat4d Quat4dOnes {1.0,1.0,1.0,1.0};
static constexpr Quat4d Quat4dW    {0.0,0.0,0.0,1.0};
static constexpr Quat4d Quat4dX    {1.0,0.0,0.0,0.0};
static constexpr Quat4d Quat4dY    {0.0,1.0,0.0,0.0};
static constexpr Quat4d Quat4dZ    {0.0,0.0,1.0,0.0};
static constexpr Quat4d Quat4dmin  {-1.e-300,-1.e-300,-1.e-300,-1.e-300};
static constexpr Quat4d Quat4dmax  {+1.e-300,+1.e-300,+1.e-300,+1.e-300};

static constexpr Quat4f Quat4fNAN  {NAN,NAN,NAN,NAN};
static constexpr Quat4f Quat4fZero {0.0f,0.0f,0.0f,0.0f};
static constexpr Quat4f Quat4fOnes {1.0f,1.0f,1.0f,1.0f};
static constexpr Quat4f Quat4fW    {0.0f,0.0f,0.0f,1.0f};
static constexpr Quat4f Quat4fX    {1.0f,0.0f,0.0f,0.0f};
static constexpr Quat4f Quat4fY    {0.0f,1.0f,0.0f,0.0f};  
static constexpr Quat4f Quat4fZ    {0.0f,0.0f,1.0f,0.0f};
static constexpr Quat4f Quat4fmin  {-1.e+37,-1.e+37,-1.e+37,-1.e+37};
static constexpr Quat4f Quat4fmax  {+1.e+37,+1.e+37,+1.e+37,+1.e+37};


// default quaternion poses
// http://www.euclideanspace.com/maths/geometry/rotations/conversions/eulerToQuaternion/steps/index.htm
// https://www.quantstart.com/articles/Mathematical-Constants-in-C
//qCamera.fromMatrixT( Mat3f{ -1.0,0.0,0.0,  0.0,0.0, 1.0,  0.0, 1.0,0.0 } );  _Lprint( qCamera, "\nmXZY" );
//qCamera.fromMatrixT( Mat3f{  1.0,0.0,0.0,  0.0,0.0,-1.0,  0.0, 1.0,0.0 } );  _Lprint( qCamera, "\nXmZY" );
//qCamera.fromMatrixT( Mat3f{  1.0,0.0,0.0,  0.0,0.0, 1.0,  0.0,-1.0,0.0 } );  _Lprint( qCamera, "\nXZmY" );
//qCamera.fromMatrixT( Mat3f{ 0.0,-1.0,0.0,   1.0,0.0,0.0,   0.0,0.0, 1.0 } ); _Lprint( qCamera, "\nmYXZ" );
//qCamera.fromMatrixT( Mat3f{ 0.0, 1.0,0.0,  -1.0,0.0,0.0,   0.0,0.0, 1.0 } ); _Lprint( qCamera, "\nYmXZ" );
//qCamera.fromMatrixT( Mat3f{ 0.0, 1.0,0.0,   1.0,0.0,0.0,   0.0,0.0,-1.0 } ); _Lprint( qCamera, "\nYXmZ" );
//qCamera.fromMatrixT( Mat3f{ 0.0,-1.0,0.0,  0.0,0.0, 1.0,   1.0,0.0,0.0  } ); _Lprint( qCamera, "\nmYZX" );
//qCamera.fromMatrixT( Mat3f{ 0.0, 1.0,0.0,  0.0,0.0,-1.0,   1.0,0.0,0.0  } ); _Lprint( qCamera, "\nYmZX" );
//qCamera.fromMatrixT( Mat3f{ 0.0, 1.0,0.0,  0.0,0.0, 1.0,  -1.0,0.0,0.0  } ); _Lprint( qCamera, "\nYZmX" );
//qCamera.fromMatrixT( Mat3f{ 0.0,0.0,-1.0,   0.0, 1.0,0.0,   1.0,0.0,0.0  } ); _Lprint( qCamera, "\nmZYX" );
//qCamera.fromMatrixT( Mat3f{ 0.0,0.0, 1.0,   0.0,-1.0,0.0,   1.0,0.0,0.0  } ); _Lprint( qCamera, "\nZmYX" );
//qCamera.fromMatrixT( Mat3f{ 0.0,0.0, 1.0,   0.0, 1.0,0.0,  -1.0,0.0,0.0  } ); _Lprint( qCamera, "\nZYmX" );
//static constexpr Quat4f Quat4fXYZ      = Quat4f{       0.0,        0.0,       0.0,       1.0 };
//static constexpr Quat4f Quat4fmXYmZ    = Quat4f{       0.0,        1.0,       0.0,       0.0 };
//static constexpr Quat4f Quat4fmXZY     = Quat4f{       0.0,  M_SQRT1_2, M_SQRT1_2,       0.0 };
//static constexpr Quat4f Quat4fXmZY     = Quat4f{-M_SQRT1_2,        0.0,       0.0, M_SQRT1_2 };
//static constexpr Quat4f Quat4fXZmY     = Quat4f{ M_SQRT1_2,        0.0,       0.0, M_SQRT1_2 };
//static constexpr Quat4f Quat4fmZYX     = Quat4f{       0.0,  M_SQRT1_2,       0.0, M_SQRT1_2 };
//static constexpr Quat4f Quat4fZmYX     = Quat4f{ M_SQRT1_2,        1.0, M_SQRT1_2,       0.0 };
//static constexpr Quat4f Quat4fZYmX     = Quat4f{       0.0, -M_SQRT1_2,       0.0, M_SQRT1_2 };
//static constexpr Quat4f Quat4fFront    = Quat4fmXYmZ;
//static constexpr Quat4f Quat4fBack     = Quat4fXYZ;
//static constexpr Quat4f Quat4fTop      = Quat4fXmZY;
//static constexpr Quat4f Quat4fBotton   = Quat4fXZmY;
//static constexpr Quat4f Quat4fLeft     = Quat4fmZYX;
//static constexpr Quat4f Quat4fRight    = Quat4fZYmX;

static constexpr Quat4f Quat4fIdentity {       0.f,       0.f,       0.f,       1.f };
static constexpr Quat4f Quat4fBack     {       0.f,       0.f,       0.f,       1.f };
static constexpr Quat4f Quat4fFront    {       0.f,       1.f,       0.f,       0.f };
static constexpr Quat4f Quat4fTop      {-M_SQRT1_2,       0.f,       0.f, M_SQRT1_2 };
static constexpr Quat4f Quat4fBotton   { M_SQRT1_2,       0.f,       0.f, M_SQRT1_2 };
static constexpr Quat4f Quat4fLeft     {       0.f,  M_SQRT1_2,       0.f, M_SQRT1_2 };
static constexpr Quat4f Quat4fRight    {       0.f, -M_SQRT1_2,       0.f, M_SQRT1_2 };

static constexpr Quat4d Quat4dIdentity {       0.,        0.,       0.,       1.  };
static constexpr Quat4d Quat4dBack     {       0.,        0.,       0.,       1.  };
static constexpr Quat4d Quat4dFront    {       0.,        1.,       0.,       0.  };
static constexpr Quat4d Quat4dTop      {-M_SQRT1_2,       0.,       0., M_SQRT1_2 };
static constexpr Quat4d Quat4dBotton   { M_SQRT1_2,       0.,       0., M_SQRT1_2 };
static constexpr Quat4d Quat4dLeft     {       0.,  M_SQRT1_2,      0., M_SQRT1_2 };
static constexpr Quat4d Quat4dRight    {       0., -M_SQRT1_2,      0., M_SQRT1_2 };

// -------- this is used for camera
static constexpr Quat4f qFront {-M_SQRT1_2,      0.f,     0.f,M_SQRT1_2 };
static constexpr Quat4f qBack  {      0.0f, M_SQRT1_2,M_SQRT1_2,     0.0f };
static constexpr Quat4f qTop   { 1.0f, 0.0f, 0.0f, 0.0f};
static constexpr Quat4f qBottom{ 0.0f, 0.0f, 0.0f, 1.0f};
static constexpr Quat4f qLeft  {-0.5f, 0.5f, 0.5f, 0.5f};
static constexpr Quat4f qRight { 0.5f, 0.5f, 0.5f,-0.5f};


inline void convert( const Quat4f& from, Quat4d& to ){ to.x=from.x;        to.y=from.y;        to.z=from.z;        to.w=from.w;        };
inline void convert( const Quat4d& from, Quat4f& to ){ to.x=(float)from.x; to.y=(float)from.y; to.z=(float)from.z; to.w=(float)from.w; };

inline int print( const Quat4f& v){ return printf( "%g %g %g %g", v.x, v.y, v.z, v.w ); };
inline int print( const Quat4d& v){ return printf( "%g %g %g %g", v.x, v.y, v.z, v.w ); };
inline int print( const Quat4i& v){ return printf( "%i %i %i %i", v.x, v.y, v.z, v.w ); };


template <class T>
class Quat8T { public:
	union{
        struct{ Quat4T<T> lo,hi; };
		T array[8];
	};
};
using Quat8i = Quat8T< int>;
using Quat8f = Quat8T< float>;
using Quat8d = Quat8T< double >;

#endif

