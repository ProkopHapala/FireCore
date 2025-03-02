
#ifndef  Vec2_h
#define  Vec2_h

#include "fastmath.h"

//template <class T,class VEC>
template <class T>
class Vec2T {
	using VEC = Vec2T<T>;
	public:
	union{
		struct{ T x,y; };
		struct{ T a,b; };
		struct{ T i,j; };
		T array[2];
	};

	// ===== methods

	inline bool order       (){ if(a>b){ SWAP(a,b,T); return true; }; return false; }
	inline bool orderReverse(){ if(a<b){ SWAP(a,b,T); return true; }; return false; }
    inline VEC  swaped(){return VEC{b,a}; };

    inline explicit operator Vec2T<float >()const{ return (Vec2T<float >){(float)x,(float)y}; }
	inline explicit operator Vec2T<double>()const{ return (Vec2T<double>){(double)x,(double)y}; }
	inline explicit operator Vec2T<int   >()const{ return (Vec2T<int   >){(int)x,(int)y}; }


	inline void set( T f            ) { x=f;   y=f;   };
    inline void set( T fx, T fy  ) { x=fx;  y=fy;  };
    inline void set( const VEC& v      ) { x=v.x; y=v.y; };
	inline void set( T* arr         ) { x=arr[0]; y=arr[1]; };

	//inline Vec2T(){};
	//inline Vec2T( T f            ){ x=f;   y=f; };
	//inline Vec2T( T fx, T fy  ){ x=fx;  y=fy; };
	//inline Vec2T( const VEC& v      ){ x=v.x; y=v.y; };
	//inline Vec2T( T* arr         ){ x=arr[0]; y=arr[1]; };

    inline void get( T& fx, T& fy ) { fx=x;  fy=y;        };
	inline void get( T* arr                    ) { arr[0]=x; arr[1]=y; };

    inline void add( T f ) { x+=f; y+=f; };
    inline void mul( T f ) { x*=f; y*=f; };

    inline void add( const VEC&  v ) { x+=v.x; y+=v.y; };
    inline void sub( const VEC&  v ) { x-=v.x; y-=v.y; };
    inline void mul( const VEC&  v ) { x*=v.x; y*=v.y; };
    inline void div( const VEC&  v ) { x/=v.x; y/=v.y; };

    inline void add( T fx, T fy) { x+=fx; y+=fy; };
    inline void sub( T fx, T fy) { x-=fx; y-=fy; };
    inline void mul( T fx, T fy ) { x*=fx; y*=fy; };
    inline void div( T fx, T fy ) { x/=fx; y/=fy; };

    inline void set_inv( const VEC&  v ) { x=1/v.x; y=1/v.y; };

	inline void set_add( const VEC& a, T f ){ x=a.x+f; y=a.y+f; };
	inline void set_mul( const VEC& a, T f ){ x=a.x*f; y=a.y*f; };
	inline void set_mul( const VEC& a, const VEC& b, T f ){ x=a.x*b.x*f; y=a.y*b.y*f; };

	inline void set_add( const VEC& a, const VEC& b ){ x=a.x+b.x; y=a.y+b.y; };
	inline void set_sub( const VEC& a, const VEC& b ){ x=a.x-b.x; y=a.y-b.y; };
	inline void set_mul( const VEC& a, const VEC& b ){ x=a.x*b.x; y=a.y*b.y; };
	inline void set_div( const VEC& a, const VEC& b ){ x=a.x/b.x; y=a.y/b.y; };

	inline void add_mul( const VEC& a, T f                ){ x+=a.x*f;     y+=a.y*f;     };
	inline void add_mul( const VEC& a, const VEC& b          ){ x+=a.x*b.x;   y+=a.y*b.y;   };
	inline void sub_mul( const VEC& a, const VEC& b          ){ x-=a.x*b.x;   y-=a.y*b.y;   };
	inline void add_mul( const VEC& a, const VEC& b, T f  ){ x+=a.x*b.x*f; y+=a.y*b.y*f; };

	inline void set_add_mul( const VEC& a, const VEC& b, T f ){ x= a.x + f*b.x;     y= a.y + f*b.y; };

	inline void set_lincomb( T fa, const VEC& a, T fb, const VEC& b ){ x = fa*a.x + fb*b.x;  y = fa*a.y + fb*b.y; };
	inline void add_lincomb( T fa, const VEC& a, T fb, const VEC& b ){ x+= fa*a.x + fb*b.x;  y+= fa*a.y + fb*b.y; };

	inline void set_lincomb( T fa, T fb, T fc, const VEC& a, const VEC& b, const VEC& c ){ x = fa*a.x + fb*b.x + fc*c.x;  y = fa*a.y + fb*b.y + fc*c.y; };
	inline void add_lincomb( T fa, T fb, T fc, const VEC& a, const VEC& b, const VEC& c ){ x+= fa*a.x + fb*b.x + fc*c.x;  y+= fa*a.y + fb*b.y + fc*c.y; };

	//inline void set_lincomb_perp( const VEC& v, const VEC& cs ){ x=v.x*cs.x - v.y*cs.y; y=v.y*cs.x + v.x*cs.y; }; //  NOT REQUIRED: wequivalent to  set_mul_cmplx

    inline VEC operator+ ( T f   ) const { VEC vo; vo.x=x+f; vo.y=y+f; return vo; };
    inline VEC operator* ( T f   ) const { VEC vo; vo.x=x*f; vo.y=y*f; return vo; };

    inline VEC operator+ ( const VEC& vi ) const { VEC vo; vo.x=x+vi.x; vo.y=y+vi.y; return vo; };
    inline VEC operator- ( const VEC& vi ) const { VEC vo; vo.x=x-vi.x; vo.y=y-vi.y; return vo; };
    inline VEC operator* ( const VEC& vi ) const { VEC vo; vo.x=x*vi.x; vo.y=y*vi.y; return vo; };
    inline VEC operator/ ( const VEC& vi ) const { VEC vo; vo.x=x/vi.x; vo.y=y/vi.y; return vo; };

	inline T dot      ( const VEC& a ) const { return x*a.x + y*a.y; };
	inline T dot_perp ( const VEC& a ) const { return y*a.x - x*a.y; };
	inline T norm2(              ) const { return x*x + y*y;     };

	inline T norm ( ) const { return  sqrt( x*x + y*y ); };
    inline T normalize() {
		T norm  = sqrt( x*x + y*y );
		T inVnorm = 1.0/norm;
		x *= inVnorm;    y *= inVnorm;
		return norm;
    };

    inline T dist2( const VEC& a) const { VEC d; T dx = x-a.x; T dy = y-a.y; return dx*dx + dy*dy; }
    inline T dist ( const VEC& a) const { return sqrt( dist2(a) ); }

    inline double makePerpUni( const VEC& a ) { double cdot=x*a.x+y*a.y; x-=a.x*cdot; y-=a.y*cdot; return cdot; }

    inline void     perp(              )  { T x_=x; x=-y;   y=x_;  }
	inline void set_perp( const VEC& a )  {         x=-a.y; y=a.x; }
	inline T cross ( const VEC& a ) const { return x*a.y - y*a.x;  }

    inline T crammer( const VEC& v, const VEC& a, const VEC& b ){
        T det   = a.x * b.y - a.y * b.x;
        T idet  = 1/det;
        x = (v.x * b.y  - v.y * b.x)*idet;
        y = (a.x * v.y  - a.y * v.x)*idet;
        return det;
    }

	bool isBetweenRotations( const VEC& a, const VEC& b ){ return (cross(a)<0)&&(cross(b)>0);  }

    // what about using bitwize operators for this ^ (multiplication) and | (division) ?
    //inline VEC operator%(const VEC& a, const VEC& b)const{ return VEC{ a.x*b.x - a.y*b.y, a.y*b.x + a.x*b.y; };  // mul_cmplx
    //inline VEC operator|(const VEC& a, const VEC& b)const{ return VEC{ a.x*b.x + a.y*b.y, a.y*b.x - a.x*b.y; };  // udiv_cmplx

    inline VEC operator%(const VEC& b)const{ return Vec2T<T>{ x*b.x - y*b.y, y*b.x + x*b.y }; }  // mul_cmplx
    inline VEC operator|(const VEC& b)const{ return Vec2T<T>{ x*b.x + y*b.y, y*b.x - x*b.y }; }  // udiv_cmplx

	inline void     mul_cmplx (               const VEC& b ){                            T x_ =    x*b.x -   y*b.y;         y =    y*b.x +   x*b.y;       x=x_;  }
	inline void     udiv_cmplx(               const VEC& b ){                            T x_ =    x*b.x +   y*b.y;         y =    y*b.x -   x*b.y;       x=x_;  }
	inline void pre_mul_cmplx ( const VEC& a               ){                            T x_ =  a.x*  x - a.y*  y;         y =  a.y*  x + a.x*  y;       x=x_;  }
	inline void set_mul_cmplx ( const VEC& a, const VEC& b ){                            T x_ =  a.x*b.x - a.y*b.y;         y =  a.y*b.x + a.x*b.y;       x=x_;  }
	inline void set_udiv_cmplx( const VEC& a, const VEC& b ){                            T x_ =  a.x*b.x + a.y*b.y;         y =  a.y*b.x - a.x*b.y;       x=x_;  }
	inline void set_div_cmplx ( const VEC& a, const VEC& b ){ T ir2 = 1/b.norm2();  T x_ = (a.x*b.x + a.y*b.y)*ir2;    y = (a.y*b.x - a.x*b.y)*ir2;    x=x_;  }

	inline void fromAngle        ( T phi ){	x = cos( phi );	y = sin( phi );	    }
	inline void fromAngle_taylor2( T phi ){	sincos_taylor2<T>( phi, y, x );	}
	inline void fromCos          ( T ca  ){  x=ca; y=sqrt(1-ca*ca); }
	inline void fromSin          ( T sa  ){  y=sa; x=sqrt(1-sa*sa); }
	inline T    toAngle()const{ return atan2(y,x); }

	inline void rotate( T phi ){
		T bx = cos( phi );   		  T by = sin( phi );
		T x_ =    x*bx -   y*by;         y =    y*bx +   x*by;       x=x_;
	}

	inline void rotate_taylor2( T phi ){
		T bx,by;  sincos_taylor2<T>( phi, by, bx );
		T x_ =    x*bx -   y*by;         y =    y*bx +   x*by;       x=x_;
	}

	inline T along_hat( const VEC& hat, const VEC& p ){ VEC ap; ap.set( p.x-x, p.y-y ); return hat.dot( ap ); }
	inline T along    ( const VEC& b,   const VEC& p ){
		VEC ab,ap;
		ab.set( b.x - x, b.y - y );
		ap.set( p.x - x, p.y - y );
		return ab.dot(ap) / ab.norm(ab);
	}

	inline T angle( const VEC& a ){
		T d =  dot ( a );
		T c = cross( a );
		//printf( "angle: dot %g cross %g \n", d, c );
		return atan2( c, d );
	}

	inline T totprod(){ return x*y; }
	inline bool allEqual(T f)const{ return (x==f)&&(y==f);};
    inline bool anyEqual(T f)const{ return (x==f)||(y==f);};

	inline void fromLinearSolution( const VEC& va, const VEC& vb, const VEC& rhs ){
        T invD = 1/va .cross(vb);
        T Da   =  rhs.cross(vb);
        T Db   = -rhs.cross(va);
        x = Da*invD;
        y = Db*invD;
	};

	inline void enclose(T v){ x=fmin(x,v); y=fmax(y,v); }


    inline bool isLower  ( const VEC& vmax ) const { return (x<vmax.x)&&(y<vmax.y); }
    inline bool isGreater( const VEC& vmin ) const { return (x>=vmin.x)&&(y>=vmin.y); }
    inline bool isBetween( const VEC& vmin, const VEC& vmax ) const { return (x>=vmin.x)&&(x<vmax.x)&&(y>=vmin.y)&&(y<vmax.y); }

    inline VEC& setIfLower  (const VEC& a){ if(a.x<x)x=a.x;if(a.y<y)y=a.y; return *this; }
    inline VEC& setIfGreater(const VEC& a){ if(a.x>x)x=a.x;if(a.y>y)y=a.y; return *this; }
    inline void update_bounds(VEC& pmin,VEC& pmax)const{ 
        if ( x < pmin.x ){ pmin.x=x; } else if ( x > pmax.x ){ pmax.x=x; };
        if ( y < pmin.y ){ pmin.y=y; } else if ( y > pmax.y ){ pmax.y=y; };
        //return *this;
    }

	// === static functions:

	static inline T dot      (VEC& a, VEC& b){ return a.x*b.y - a.y*b.x; };
	static inline T cross    (VEC& a, VEC& b){ return a.x*b.x + a.y*b.y ; };
	static inline VEC  mul_cmplx(VEC& a, VEC& b){ return (VEC){ a.x*b.x-a.y*b.y, a.y*b.x+a.x*b.y };  }

	static inline VEC newFromAngle( T phi ){ VEC v;  v.x = cos( phi ); v.y = sin( phi ); return v; }

};

using Vec2i = Vec2T<int>;
using Vec2f = Vec2T<float>;
using Vec2d = Vec2T<double>;

static constexpr Vec2d Vec2dNAN  {NAN,NAN};
static constexpr Vec2d Vec2dZero {0.0,0.0};
static constexpr Vec2d Vec2dOnes {1.0,1.0};
static constexpr Vec2d Vec2dX    {1.0,0.0};
static constexpr Vec2d Vec2dY    {0.0,1.0};
static constexpr Vec2d Vec2dmin  {-1.e-300,-1.e-300};
static constexpr Vec2d Vec2dmax  {+1.e-300,+1.e-300};

static constexpr Vec2f Vec2fNAN  {NAN,NAN};
static constexpr Vec2f Vec2fZero {0.0f,0.0f};
static constexpr Vec2f Vec2fOnes {1.0f,1.0f};
static constexpr Vec2f Vec2fX    {1.0f,0.0f};
static constexpr Vec2f Vec2fY    {0.0f,1.0f};
static constexpr Vec2f Vec2fmin  {-1.e+37,-1.e+37};
static constexpr Vec2f Vec2fmax  {+1.e+37,+1.e+37};

static constexpr Vec2i Vec2iZero {0,0};
static constexpr Vec2i Vec2iOnes {1,1};
static constexpr Vec2i Vec2iX    {1,0};
static constexpr Vec2i Vec2iY    {0,1};

inline uint64_t scalar_id  ( const Vec2i& v){ return (((uint64_t)v.a)<<32)|v.b; }
inline uint64_t symetric_id( const Vec2i& v){ if( v.a>v.b ){ return (((uint64_t)v.b)<<32)|v.a; }else{ return (((uint64_t)v.a)<<32)|v.b; }}

inline void convert( const Vec2f& from, Vec2d& to ){ to.x=from.x;        to.y=from.y;        };
inline void convert( const Vec2d& from, Vec2f& to ){ to.x=(float)from.x; to.y=(float)from.y; };

inline Vec2f toFloat( const Vec2d& from){ return (Vec2f){(float)from.x,(float)from.y}; }

inline int print( const Vec2f&  v){ return printf( "%g %g", v.x, v.y ); };
inline int print( const Vec2d&  v){ return printf( "%lg %lg", v.x, v.y ); };
inline int print( const Vec2i&  v){ return printf( "%i %i", v.x, v.y ); };

#endif


