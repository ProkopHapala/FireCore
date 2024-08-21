
#ifndef  Vec3_h
#define  Vec3_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
//#include "gonioApprox.h"
#include "Vec2.h"

//template <class T,class VEC>
template <class T>
class Vec3T{
	using VEC  = Vec3T<T>;
	using VEC2 = Vec2T<T>;
	public:
	union{
		struct{ T x,y,z; };
		struct{ T a,b,c; };
		struct{ T i,j,k; };
		T array[3];
	};

	// Constructors would prevent us from making Unions etc. so don't do it
	// https://stackoverflow.com/questions/4178175/what-are-aggregates-and-pods-and-how-why-are-they-special
	// but here it seems to work  https://www.youtube.com/watch?v=14Cyfz_tE20&index=10&list=PLlrATfBNZ98fqE45g3jZA_hLGUrD4bo6_
	//Vec3T() = default;
	//constexpr Vec3T(T x_, T y_, T z_ ): x(x_),y(y_),z(z_){};
	//constexpr Vec3T() = default;
	//constexpr Vec3T(T x_, T y_, T z_ ): x(x_),y(y_),z(z_){};

	// ===== methods

	// Automatic conversion (works) but would be problematic
	//inline operator Vec3T<float >()const{ return (Vec3T<float >){(float)x,(float)y,(float)z}; }
	//inline operator Vec3T<double>()const{ return (Vec3T<double>){(double)x,(double)y,(double)z}; }
	//inline operator Vec3T<int   >()const{ return (Vec3T<int   >){(int)x,(int)y,(int)z}; }

	// Explicit conversion
	inline explicit operator Vec3T<double>()const{ return Vec3T<double>{(double)x,(double)y,(double)z}; }
    inline explicit operator Vec3T<float >()const{ return Vec3T<float >{(float )x,(float )y,(float )z}; }
	inline explicit operator Vec3T<int   >()const{ return Vec3T<int   >{(int   )x,(int   )y,(int   )z}; }
    inline explicit operator Vec3T<int8_t>()const{ return Vec3T<int8_t>{(int8_t)x,(int8_t)y,(int8_t)z}; }

	//inline operator (const char*)()const{ return (; }

	//inline Vec3T<double> toDouble()const{ return (Vec3T<double>){ (double)x,(double)y,(double)z}; }
	//inline Vec3T<float > toFloat ()const{ return (Vec3T<float >){ (float)x, (double)y,(double)z}; }
	//inline Vec3T<int >   toInt   ()const{ return (Vec3T<int   >){ (int)x,      (int)y,   (int)z}; }

	// swizzles
	inline VEC2 xy() const { return {x,y}; };
	inline VEC2 xz() const { return {x,z}; };
	inline VEC2 yz() const { return {y,z}; };
    inline VEC2 yx() const { return {y,x}; };
	inline VEC2 zx() const { return {z,x}; };
	inline VEC2 zy() const { return {z,y}; };
    inline VEC xzy() const { return {x,z,y}; };
	inline VEC yxz() const { return {y,x,z}; };
	inline VEC yzx() const { return {y,z,x}; };
	inline VEC zxy() const { return {z,x,y}; };
	inline VEC zyx() const { return {z,y,x}; };

    
    inline VEC  swaped (const Vec3T<int>& inds          ) const{ return VEC{array[inds.x],array[inds.y],array[inds.z]}; };
    inline void swap   (const Vec3T<int>& inds          ){ *this=swaped(inds); };
    //inline void swap_to(const Vec3T<int>& inds, VEC& out) const{ out.x=array[inds.x]; out.y=array[inds.y]; out.z=array[inds.z]; };

	inline VEC& set( T f              ) { x=f;   y=f;   z=f;   return *this; };
    inline VEC& set( T fx, T fy, T fz ) { x=fx;  y=fy;  z=fz;  return *this; };
    inline VEC& set( const VEC& v     ) { x=v.x; y=v.y; z=v.z; return *this; };
	inline VEC& set( T* arr           ) { x=arr[0]; y=arr[1]; z=arr[2]; return *this; };

    inline VEC& get( T& fx, T& fy, T& fz ) { fx=x;  fy=y;  fz=z;           return *this; };
	inline VEC& get( T* arr              ) { arr[0]=x; arr[1]=y; arr[2]=z; return *this; };

    inline VEC& add( T f ) { x+=f; y+=f; z+=f; return *this;};
    inline VEC& mul( T f ) { x*=f; y*=f; z*=f; return *this;};

    inline VEC& add( const VEC&  v ) { x+=v.x; y+=v.y; z+=v.z; return *this;};
    inline VEC& sub( const VEC&  v ) { x-=v.x; y-=v.y; z-=v.z; return *this;};
    inline VEC& mul( const VEC&  v ) { x*=v.x; y*=v.y; z*=v.z; return *this;};
    inline VEC& div( const VEC&  v ) { x/=v.x; y/=v.y; z/=v.z; return *this;};

    inline VEC& set_inv( const VEC&  v ) { x=1/v.x; y=1/v.y; z=1/v.z; return *this; };
    inline VEC  get_inv()                { VEC o; o.x=1/x; o.y=1/y; o.z=1/z; return o; };

    inline VEC& add( T fx, T fy, T fz ) { x+=fx; y+=fy; z+=fz; return *this;};
    inline VEC& sub( T fx, T fy, T fz ) { x-=fx; y-=fy; z-=fz; return *this;};
    inline VEC& mul( T fx, T fy, T fz ) { x*=fx; y*=fy; z*=fz; return *this;};
    inline VEC& div( T fx, T fy, T fz ) { x/=fx; y/=fy; z/=fz; return *this;};

	inline VEC& set_add( const VEC& a, T f ){ x=a.x+f; y=a.y+f; z=a.z+f; return *this;};
	inline VEC& set_mul( const VEC& a, T f ){ x=a.x*f; y=a.y*f; z=a.z*f; return *this;};
	inline VEC& set_mul( const VEC& a, const VEC& b, T f ){ x=a.x*b.x*f; y=a.y*b.y*f; z=a.z*b.z*f; return *this; };

	inline VEC& set_add( const VEC& a, const VEC& b ){ x=a.x+b.x; y=a.y+b.y; z=a.z+b.z; return *this; };
	inline VEC& set_sub( const VEC& a, const VEC& b ){ x=a.x-b.x; y=a.y-b.y; z=a.z-b.z; return *this; };
	inline VEC& set_mul( const VEC& a, const VEC& b ){ x=a.x*b.x; y=a.y*b.y; z=a.z*b.z; return *this; };
	inline VEC& set_div( const VEC& a, const VEC& b ){ x=a.x/b.x; y=a.y/b.y; z=a.z/b.z; return *this; };

	inline VEC& add_mul( const VEC& a, T f                ){ x+=a.x*f;     y+=a.y*f;     z+=a.z*f;   return *this;};
	inline VEC& add_mul( const VEC& a, const VEC& b       ){ x+=a.x*b.x;   y+=a.y*b.y;   z+=a.z*b.z; return *this;};
	inline VEC& sub_mul( const VEC& a, const VEC& b       ){ x-=a.x*b.x;   y-=a.y*b.y;   z-=a.z*b.z; return *this;};
	inline VEC& add_mul( const VEC& a, const VEC& b, T f  ){ x+=a.x*b.x*f; y+=a.y*b.y*f; z+=a.z*b.z*f;   return *this;};


	inline VEC& set_add_mul( const VEC& a, const VEC& b, T f ){ x= a.x + f*b.x;     y= a.y + f*b.y;     z= a.z + f*b.z;  return *this;};


	inline VEC& set_lincomb( T fa, const VEC& a, T fb, const VEC& b ){ x = fa*a.x + fb*b.x;  y = fa*a.y + fb*b.y;  z = fa*a.z + fb*b.z; return *this;};
	inline VEC& add_lincomb( T fa, const VEC& a, T fb, const VEC& b ){ x+= fa*a.x + fb*b.x;  y+= fa*a.y + fb*b.y;  z+= fa*a.z + fb*b.z; return *this;};

	inline VEC& set_lincomb( T fa, T fb, T fc, const VEC& a, const VEC& b, const VEC& c ){ x = fa*a.x + fb*b.x + fc*c.x;  y = fa*a.y + fb*b.y + fc*c.y;  z = fa*a.z + fb*b.z + fc*c.z; return *this;};
	inline VEC& add_lincomb( T fa, T fb, T fc, const VEC& a, const VEC& b, const VEC& c ){ x+= fa*a.x + fb*b.x + fc*c.x;  y+= fa*a.y + fb*b.y + fc*c.y;  z+= fa*a.z + fb*b.z + fc*c.z; return *this;};

    inline VEC& set_lincomb( const VEC& fs, const VEC& a, const VEC& b, const VEC& c ){ x = fs.a*a.x + fs.b*b.x + fs.c*c.x;  y = fs.a*a.y + fs.b*b.y + fs.c*c.y;  z = fs.a*a.z + fs.b*b.z + fs.c*c.z; return *this;};
	inline VEC& add_lincomb( const VEC& fs, const VEC& a, const VEC& b, const VEC& c ){ x+= fs.a*a.x + fs.b*b.x + fs.c*c.x;  y+= fs.a*a.y + fs.b*b.y + fs.c*c.y;  z+= fs.a*a.z + fs.b*b.z + fs.c*c.z; return *this;};

    inline VEC& set_cross( const VEC& a, const VEC& b ){ x =a.y*b.z-a.z*b.y; y =a.z*b.x-a.x*b.z; z =a.x*b.y-a.y*b.x; return *this;};
	inline VEC& add_cross( const VEC& a, const VEC& b ){ x+=a.y*b.z-a.z*b.y; y+=a.z*b.x-a.x*b.z; z+=a.x*b.y-a.y*b.x; return *this;};
	inline VEC& sub_cross( const VEC& a, const VEC& b ){ x-=a.y*b.z-a.z*b.y; y-=a.z*b.x-a.x*b.z; z-=a.x*b.y-a.y*b.x; return *this;};

	T makeOrthoU( const VEC& a ){ T c = dot(a);           add_mul(a, -c); return c; }
	T makeOrtho ( const VEC& a ){ T c = dot(a)/a.norm2(); add_mul(a, -c); return c; }

	inline void operator+=( const VEC& v ){ x+=v.x; y+=v.y; z=z+v.z; };
    inline void operator*=( const VEC& v ){ x*=v.x; y*=v.y; z=z*v.z; };

    //inline VEC operator+ ( T f   ) const { return VEC{ x+f, y+f, z+f }; };
    inline VEC operator* ( T f   ) const { return VEC{ x*f, y*f, z*f }; };

    inline VEC operator+ ( const VEC& vi ) const { return VEC{ x+vi.x, y+vi.y, z+vi.z }; };
    inline VEC operator- ( const VEC& vi ) const { return VEC{ x-vi.x, y-vi.y, z-vi.z }; };
    inline VEC operator* ( const VEC& vi ) const { return VEC{ x*vi.x, y*vi.y, z*vi.z }; };
    inline VEC operator/ ( const VEC& vi ) const { return VEC{ x/vi.x, y/vi.y, z/vi.z }; };

    inline T bidot  ( const VEC& a, const VEC& b ) const { return x*a.x*b.x + y*a.y*b.y + z*a.z*b.z;  };
    inline T antidot( const VEC& a, const VEC& b ) const { return x*a.y*b.z + y*a.z*b.x + z*a.x*b.y;  };

    inline T triple_product( const VEC& a, const VEC& b ) const { return x*(a.y*b.z-a.z*b.y) + y*(a.z*b.x-a.x*b.z) + z*(a.x*b.y-a.y*b.x);  };  // https://en.wikipedia.org/wiki/Triple_product

	inline T dot  ( const VEC& a ) const { return x*a.x + y*a.y + z*a.z;  };
	inline T norm2(              ) const { return x*x + y*y + z*z;        };
	inline T norm ( ) const { return  sqrt( x*x + y*y + z*z ); };
    inline T normalize() {
		T norm  = sqrt( x*x + y*y + z*z );
		T inVnorm = 1.0/norm;
		x *= inVnorm;    y *= inVnorm;    z *= inVnorm;
		return norm;
    }
    inline VEC normalized()const{
        VEC v; v.set(*this);
        v.normalize();
        return v;
    }

    inline T cos_v( const VEC& b ){ return dot(b)/sqrt( norm2() * b.norm2() ); }


    inline bool tryNormalize(double errMax){
        double r2 = norm2();
        if( fabs(r2-1.0)>errMax ){
            mul( 1/sqrt(r2) );
            return true;
        }
        return false;
    }

    inline bool tryOrthogonalize( double errMax, const VEC& u ){
        double c = dot(u);
        if( fabs(c)>errMax ){
            add_mul( u, -c );
            return true;
        }
        return false;
    }

    inline VEC getOrtho( VEC& up ) const {
        up.makeOrthoU(*this); up.normalize();
        VEC out; out.set_cross(*this,up);
        return out;
	}

    inline void sort(){
        if(x>y) _swap(x,y);
        if(y>z) _swap(y,z);
        if(x>y) _swap(x,y);
    }

	inline T normalize_taylor3(){
        // sqrt(1+x) ~= 1 + 0.5*x - 0.125*x*x
        // sqrt(r2) = sqrt((r2-1)+1) ~= 1 + 0.5*(r2-1)
        // 1/sqrt(1+x) ~= 1 - 0.5*x + (3/8)*x^2 - (5/16)*x^3 + (35/128)*x^4 - (63/256)*x^5
        T dr2  = x*x+y*y+z*z-1;
        T invr = 1 + dr2*( -0.5 + dr2*( 0.375 + dr2*-0.3125 ) );
        x*=invr;
        y*=invr;
        z*=invr;
        //return *this;
        return invr;
	}

	inline T fixSphere( const VEC& pc, T r){ sub(pc); T l=norm(); mul(r/l); add(pc); return l; }

	inline T fixSphere_taylor3( const VEC& pc, T r){
        x-=pc.x; y-=pc.y; z-=pc.z;
        T dr2  = x*x+y*y+z*z-1;
        T invr = r*(1 + dr2*( -0.5 + dr2*( 0.375 + dr2*-0.3125 ) ));
        x=x*invr+pc.x;
        y=y*invr+pc.y;
        z=z*invr+pc.z;
        return invr;
    }


	inline void getSomeOrtho( VEC& v1, VEC& v2 ) const {
        T xx = x*x;
        T yy = y*y;
		if(xx<yy){
//			x : y*vz - z*vy;
//			y : z*vx - x*vz;
//			z : x*vy - y*vx;
//			x : y*0 - z*0 ;
//			y : z*1 - x*0 ;
//			z : x*0 - y*1 ;
//			float vx = 0; float vy = z; float vz =-y;
			v1.x =  -yy -z*z;
			v1.y =  x*y;
			v1.z =  x*z;
		}else{
//			x : y*0 - z*1;
//			y : z*0 - x*0;
//			z : x*1 - y*0;
//			float vx = -z; float vy = 0; float vz = x;
			v1.x =  y*x;
			v1.y =  -z*z -xx;
			v1.z =  y*z;
		}
		v2.x = y*v1.z - z*v1.y;
		v2.y = z*v1.x - x*v1.z;
		v2.z = x*v1.y - y*v1.x;
	}

    inline VEC& drotate_omega(const VEC& w){
        T dx =y*w.z-z*w.y;
        T dy =z*w.x-x*w.z;
        T dz =x*w.y-y*w.x;
        //x+=dx; y+=dy; z+=dz;
        x-=dx; y-=dy; z-=dz;
        return *this;
    }

    inline VEC& drotate_omega2(const VEC& w){
        T dx  = y*w.z- z*w.y;
        T dy  = z*w.x- x*w.z;
        T dz  = x*w.y- y*w.x;
        T ddx =dy*w.z-dz*w.y;
        T ddy =dz*w.x-dx*w.z;
        T ddz =dx*w.y-dy*w.x;
        //x+=dx - ddx*0.5;
        //y+=dy - ddy*0.5;
        //z+=dz - ddz*0.5;
        x-=dx - ddx*0.5;
        y-=dy - ddy*0.5;
        z-=dz - ddz*0.5;
        return *this;
    }

    inline VEC& drotate_omega_csa(const VEC& w, T ca, T sa){
        T dx  =  y*w.z -  z*w.y;
        T dy  =  z*w.x -  x*w.z;
        T dz  =  x*w.y -  y*w.x;
        T ddx = dy*w.z - dz*w.y;
        T ddy = dz*w.x - dx*w.z;
        T ddz = dx*w.y - dy*w.x;
        //x+=dx*sa + ddx*ca;
        //y+=dy*sa + ddy*ca;
        //z+=dz*sa + ddz*ca;
        x-=dx*sa + ddx*ca;
        y-=dy*sa + ddy*ca;
        z-=dz*sa + ddz*ca;
        return *this;
    }

    inline VEC& drotate_omega6(const VEC& w){
        /*
        constexpr T c2 = -1.0/2;
        constexpr T c3 = -1.0/6;
        constexpr T c4 =  1.0/24;
        constexpr T c5 =  1.0/120;
        constexpr T c6 = -1.0/720;
        T r2  = w.x*w.x + w.y*w.y + w.z*w.z;
        T sa  =   1 + r2*( c3 + c5*r2 );
        T ca  =  c2 + r2*( c4 + c6*r2 );
        */
        T ca,sa;
        sincosR2_taylor(w.norm2(), sa, ca );
        drotate_omega_csa(w,ca,sa);
        return *this;
    }

    inline T getAngle(const VEC& b){
        T c = dot(b)/sqrt( norm2()*b.norm2() );
        return acos(c); 
    }

    inline T getAngle_unitary(const VEC& b){ return acos( dot(b) ); }

	// Rodrigues rotation formula: v' = cosa*v + sina*(uaxis X v) + (1-cosa)*(uaxis . v)*uaxis
	inline VEC& rotate( T angle, const VEC& axis  ){
		VEC uaxis;
		uaxis.set_mul( axis, 1/axis.norm() );
		T ca   = cos(angle);
		T sa   = sin(angle);
 		rotate_csa( ca, sa, uaxis );
 		return *this;
	};

	inline VEC& rotate_csa( T ca, T sa, const VEC& uaxis ){
		T cu = (1-ca)*dot(uaxis);
		T utx  = uaxis.y*z - uaxis.z*y;
		T uty  = uaxis.z*x - uaxis.x*z;
		T utz  = uaxis.x*y - uaxis.y*x;
		T x_ = ca*x + sa*utx + cu*uaxis.x;
		T y_ = ca*y + sa*uty + cu*uaxis.y;
		       z  = ca*z + sa*utz + cu*uaxis.z;
		x = x_; y = y_;
		return *this;
	};

    inline VEC& rotate_csa( T ca, T sa, const VEC& uaxis, const VEC& p0        ){ sub(p0); rotate_csa( ca, sa, uaxis ); add(p0); return *this; };
    inline VEC& rotate    ( T angle,    const VEC& axis , const VEC& p0        ){ sub(p0); rotate    ( angle,   axis ); add(p0); return *this; };
    inline VEC& scale     (             const VEC& sc   , const VEC& p0        ){ sub(p0); mul(sc);                     add(p0); return *this; };

	inline VEC& rotateTo( const VEC& rot0, double coef ){
        //rot.add_mul( rot0, coef ); rot.normalize();
        VEC ax; ax.set_cross( *this, rot0 );
        double sa2 = ax.norm2();
        if( sa2 < coef*coef ){
            ax.mul( 1/sqrt(sa2) ); // this is safe if coef is large enough
            double ca = sqrt( 1-coef*coef );
            rotate_csa( ca, coef, ax );
        }else{
            set(rot0);
        }
        return *this;
    }

    inline void getInPlaneRotation( const VEC& rot0, const VEC& xhat, const VEC& yhat, double& ca, double& sa ){
        double x0 = rot0.dot(xhat);
        double y0 = rot0.dot(yhat);
        double x_ = dot(xhat);
        double y_ = dot(yhat);
        // http://mathworld.wolfram.com/ComplexDivision.html
        double renorm = 1.0/sqrt( (x0*x0 + y0*y0)*(x_*x_ + y_*y_) );
        ca = ( x0*x_ + y0*y_ ) * renorm;
        sa = ( y0*x_ - x0*y_ ) * renorm;
    }

	inline T along_hat( const VEC& hat, const VEC& p ){ VEC ap; ap.set( p.x-x, p.y-y ); return hat.dot( ap ); }
	inline T along    ( const VEC& b,   const VEC& p ){
		VEC ab,ap;
		ab.set( b.x - x, b.y - y, b.z - z );
		ap.set( p.x - x, p.y - y, b.z - z );
		return ab.dot(ap) / ab.norm(ab);
	}

    inline void abs(){ x=_abs(x);y=_abs(y);z=_abs(z); };
	inline int maxComponent(){ return (x>y)?((x>z)?0:2):((y>z)?1:2); };
	inline int minComponent(){ return (x<y)?((x<z)?0:2):((y<z)?1:2); };

	inline void minmaxComponent(int& imin, int& imax){
        if(x>y){ imax=(x>z)?0:2; imin=(y<z)?1:2; }else{ imax=(y>z)?1:2; imin=(x<z)?0:2; };
    }

	//inline int maxComponentAbs(){
    //    Vec3d  tmp; tmp.x=fabs(x); tmp.y=fabs(y); tmp.z=fabs(z);
    //    return tmp.maxComponent();
	//}

    inline bool isLower  ( const VEC& vmax ) const { return (x<vmax.x)&&(y<vmax.y)&&(x<vmax.z); }
    inline bool isGreater( const VEC& vmin ) const { return (x>=vmin.x)&&(y>=vmin.y)&&(x>=vmin.z); }
    inline bool isBetween( const VEC& vmin, const VEC& vmax ) const { return (x>=vmin.x)&&(x<vmax.x)&&(y>=vmin.y)&&(y<vmax.y)&&(z>=vmin.z)&&(z<vmax.z); }

    inline VEC& setIfLower  (const VEC& a){ if(a.x<x)x=a.x;if(a.y<y)y=a.y;if(a.z<z)z=a.z; return *this; }
    inline VEC& setIfGreater(const VEC& a){ if(a.x>x)x=a.x;if(a.y>y)y=a.y;if(a.z>z)z=a.z; return *this; }
    inline void update_bounds(VEC& pmin,VEC& pmax)const{ 
        if ( x < pmin.x ){ pmin.x=x; } else if ( x > pmax.x ){ pmax.x=x; };
        if ( y < pmin.y ){ pmin.y=y; } else if ( y > pmax.y ){ pmax.y=y; };
        if ( z < pmin.z ){ pmin.z=z; } else if ( z > pmax.z ){ pmax.z=z; };
        //return *this;
    }

    //inline VEC min(VEC a){ return {fmin(x,a.x),fmin(y,a.y),fmin(z,a.z)}; };
    //inline VEC max(VEC a){ return {fmax(x,a.x),fmax(y,a.y),fmax(z,a.z)}; };
    //inline VEC set_min(VEC a,VEC b){ return {fmin(x,a.x),fmin(y,a.y),fmin(z,a.z)}; };
    //inline VEC set_max(VEC a,VEC b){ return {fmax(x,a.x),fmax(y,a.y),fmax(z,a.z)}; };

    inline T dist2( const VEC& a ) const { VEC d; d.set( x-a.x, y-a.y, z-a.z ); return d.norm2(); }
    inline T dist ( const VEC& a ) const { VEC d; d.set( x-a.x, y-a.y, z-a.z ); return d.norm (); }

    inline T totprod    ()const{ return x*y*z; };
    inline T sum        ()const{ return x+y+z; };
    inline bool allEqual(T f)const{ return (x==f)&&(y==f)&&(z==f);};
    inline bool anyEqual(T f)const{ return (x==f)||(y==f)||(z==f);};

    inline void octDir( int iface ){
        // mix-weights for inverse mapping of octahedron
        if(iface&1) a=-a;
        if(iface&2) b=-b;
        if(iface&4) c=-c;
        //double invr = fastInv1sqrt(  );
        //return {a*invr,b*invr,c*invr};
    }

    inline int octFace()const{ return (x>0)|((y>0)<<1)|((z>0)<<2); }

    inline int octCoord( VEC& d )const{
        int    i=0;
        //double r=0;
        if(x>0){ i|=1; d.x=x; }else{ d.x=-x; };
        if(y>0){ i|=2; d.y=y; }else{ d.y=-y; };
        if(z>0){ i|=4; d.z=z; }else{ d.z=-z; };
        double invr=1/(d.x+d.y+d.z);
        d.x*=invr; d.y*=invr; d.z*=invr;
        return i;
    }

    T angleInPlane( const VEC& a, const VEC& b ){
        T x = dot(a);
        T y = dot(b);
        return atan2( y, x );
    }

    inline VEC& setHomogenousSphericalSample( T u, T v ){
        T  r = sqrt(1-u*u);
        T  c = cos(v);
        T  s = sin(v);
        //printf( "%f %f  %f %f %f \n", u,v,  r, c, s );
        x = r*c;
        y = r*s;
        z = u;
        return *this;
    }

    inline VEC& fromRandomSphereSample(){
        setHomogenousSphericalSample( (randf()*2)-1, randf()*2*M_PI );
        return *this;
    }

    inline VEC& addRandomCube ( T d ){ x+=randf(-d,d); y+=randf(-d,d); z+=randf(-d,d);  return *this; }
    inline VEC& fromRandomCube( T d ){ x =randf(-d,d); y =randf(-d,d); z =randf(-d,d);  return *this; }
    inline VEC& fromRandomBox ( const VEC& vmin, const VEC& vmax ){ x=randf(vmin.x,vmax.x); y=randf(vmin.y,vmax.y); z=randf(vmin.z,vmax.z);  return *this; }

    inline VEC& fromLinearSolution( const VEC& va, const VEC& vb, const VEC& vc, const VEC& p ){
        // https://en.wikipedia.org/wiki/Cramer%27s_rule
        // 30 multiplications
        T Dax = vb.y*vc.z - vb.z*vc.y;
        T Day = vb.x*vc.z - vb.z*vc.x;
        T Daz = vb.x*vc.y - vb.y*vc.x;
        T idet = 1/( va.x*Dax - va.y*Day + va.z*Daz );
        x =  idet*( p.x*Dax - p.y*Day + p.z*Daz );
        y = -idet*( p.x*(va.y*vc.z - va.z*vc.y) - p.y*(va.x*vc.z - va.z*vc.x) + p.z*(va.x*vc.y - va.y*vc.x) );
        z =  idet*( p.x*(va.y*vb.z - va.z*vb.y) - p.y*(va.x*vb.z - va.z*vb.x) + p.z*(va.x*vb.y - va.y*vb.x) );
        return *this;
    }

};

template<typename VEC> inline VEC cross( VEC a, VEC b ){ return (VEC){ a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x }; }
template<typename VEC> inline VEC add  ( VEC a, VEC b ){ return (VEC){ a.x+b.x, a.z+b.z, a.z+b.z }; }

using Vec3i = Vec3T<int>;
using Vec3f = Vec3T<float>;
using Vec3d = Vec3T<double>;
using Vec3i8  = Vec3T<int8_t>;
using Vec3ui8 = Vec3T<uint8_t>;
using Vec3b   = Vec3T<bool>;

static constexpr Vec3d Vec3dNAN {NAN,NAN,NAN};
static constexpr Vec3d Vec3dZero{0.0,0.0,0.0};
static constexpr Vec3d Vec3dOne {1.0,1.0,1.0};
static constexpr Vec3d Vec3dX   {1.0,0.0,0.0};
static constexpr Vec3d Vec3dY   {0.0,1.0,0.0};
static constexpr Vec3d Vec3dZ   {0.0,0.0,1.0};
static constexpr Vec3d Vec3dmin {-1e+300,-1e+300,-1e+300};
static constexpr Vec3d Vec3dmax {+1e+300,+1e+300,+1e+300};

static constexpr Vec3f Vec3fNAN {NAN,NAN,NAN};
static constexpr Vec3f Vec3fZero{0.0f,0.0f,0.0f};
static constexpr Vec3f Vec3fOne {1.0f,1.0f,1.0f};
static constexpr Vec3f Vec3fX   {1.0f,0.0f,0.0f};
static constexpr Vec3f Vec3fY   {0.0f,1.0f,0.0f};
static constexpr Vec3f Vec3fZ   {0.0f,0.0f,1.0f};
static constexpr Vec3f Vec3fmin {-1e+37,-1e+37,-1e+37};
static constexpr Vec3f Vec3fmax {+1e+37,+1e+37,+1e+37};

static constexpr Vec3i Vec3iZero {0,0,0};
static constexpr Vec3i Vec3iOne  {1,1,1};
static constexpr Vec3i Vec3iX    {1,0,0};
static constexpr Vec3i Vec3iY    {0,1,0};
static constexpr Vec3i Vec3iZ    {0,0,1};
static constexpr Vec3i Vec3imin  {-2147483647,-2147483647,-2147483647};
static constexpr Vec3i Vec3imax  {+2147483647,+2147483647,+2147483647};


inline uint64_t scalar_id  ( const Vec3i& v){ return ( v.x | (((uint64_t)v.y)<<16) | (((uint64_t)v.z)<<32) ); }
inline Vec3i    from_id    ( uint64_t id   ){
    Vec3i vi;
    vi.x=( id & 0xFFFF ); id=id>>16;
    vi.y=( id & 0xFFFF ); id=id>>16;
    vi.z=( id & 0xFFFF );
    return vi;
}

template<typename T1,typename T2>
inline void convert(const Vec3T<T1>& i, Vec3T<T2>& o){  o.x=(T2)i.x; o.y=(T2)i.y; o.z=(T2)i.z; };

template<typename T1,typename T2>
inline Vec3T<T2> cast(const Vec3T<T1>& i){ Vec3T<T2> o; o.x=(T2)i.x; o.y=(T2)i.y; o.z=(T2)i.z; return o; };


//inline void convert( const Vec3f& from, Vec3d& to ){ to.x=from.x;        to.y=from.y;        to.z=from.z; };
//inline void convert( const Vec3d& from, Vec3f& to ){ to.x=(float)from.x; to.y=(float)from.y; to.z=(float)from.z; };
//inline Vec3f toFloat( const Vec3d& from){ return Vec3f{(float)from.x,(float)from.y,(float)from.z}; }

//inline void print(Vec3d p){printf("(%.16g,%.16g,%.16g)", p.x,p.y,p.z);};
//inline void print(Vec3f p){printf("(%.8g,%.8g,%.8g)", p.x,p.y,p.z);};
//inline void print(Vec3d p){printf("(%lg,%lg,%lg)", p.x,p.y,p.z);};
//inline void print(Vec3f p){printf("(%g,%g,%g)", p.x,p.y,p.z);};
//inline void print(Vec3i p){printf("(%i,%i,%i)", p.x,p.y,p.z);};

inline int print( const Vec3f&  v){ return printf( "%g %g %g", v.x, v.y, v.z ); };
inline int print( const Vec3d&  v){ return printf( "%g %g %g", v.x, v.y, v.z ); };
inline int print( const Vec3i&  v){ return printf( "%i %i %i", v.x, v.y, v.z ); };


template <class T>
class Vec6T { public:
	union{
        struct{ Vec3T<T> lo,hi; };
        struct{ T xx,yy,zz, yz,xz,xy; };
		T array[6];
	};
};
using Vec6i = Vec6T< int>;
using Vec6f = Vec6T< float>;
using Vec6d = Vec6T< double >;



// Structured binding support for Vec3T
// allows to do this :
// Vec3f v{1,2,3};
// auto [x,y,z] = v;

// #include <tuple>

// namespace std {
//     template <class T>
//     struct tuple_size<Vec3T<T>> : std::integral_constant<size_t, 3> {};

//     template <class T, std::size_t Index>
//     struct tuple_element<Index, Vec3T<T>> {
//         using type = T;
//     };
// }

// // Non-const get function
// template <std::size_t Index, class T>
// T& get(Vec3T<T>& v) {
//     if constexpr (Index == 0) return v.x;
//     else if constexpr (Index == 1) return v.y;
//     else if constexpr (Index == 2) return v.z;
//     else static_assert(Index < 3, "Index out of bounds for Vec3T");
// }

// // Const get function
// template <std::size_t Index, class T>
// const T& get(const Vec3T<T>& v) {
//     if constexpr (Index == 0) return v.x;
//     else if constexpr (Index == 1) return v.y;
//     else if constexpr (Index == 2) return v.z;
//     else static_assert(Index < 3, "Index out of bounds for Vec3T");
// }

// // Rvalue get function (for temporary objects)
// template <std::size_t Index, class T>
// T&& get(Vec3T<T>&& v) {
//     if constexpr (Index == 0) return std::move(v.x);
//     else if constexpr (Index == 1) return std::move(v.y);
//     else if constexpr (Index == 2) return std::move(v.z);
//     else static_assert(Index < 3, "Index out of bounds for Vec3T");
// }

#endif



