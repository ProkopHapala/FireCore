
#ifndef  VecN_h
#define  VecN_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"

template<typename Func>
void evalFunc1D(Func func, int n, double * xs, double * ys){
    for(int i=0; i<n; i++){ ys[i]=func(xs[i]); }
};

namespace VecN{

    // =======  iterration over array
	inline double norm2 (int n,const double* a                 ){ double sum =0; for (int i=0; i<n; i++ ){ double ai=a[i]; sum+=ai*ai;          } return sum; }
	inline double wnorm2(int n,const double* a,const double* w ){ double sum =0; for (int i=0; i<n; i++ ){ double ai=a[i]; sum+= (ai*ai)*w[i];  } return sum; }
	inline double wdot  (int n,const double* a,const double* b,const double* w){ double sum =0; for (int i=0; i<n; i++ ){ sum+= a[i]*b[i]*w[i];  } return sum; }

	inline double dot    (int n, const double* a, const double* b ){ double sum =0; for (int i=0; i<n; i++ ){ sum+= a[i]*b[i];             } return sum; }
	inline double sum    (int n, const double* a            ){ double sum =0; for (int i=0; i<n; i++ ){ sum+= a[i];                  } return sum; };
    inline double sum2   (int n, const double* a            ){ double sum =0; for (int i=0; i<n; i++ ){ double ai=a[i]; sum+=ai*ai;  } return sum; };
    inline double min    (int n, const double* a            ){ double amax=-1e+300; for (int i=0; i<n; i++ ){ double ai=a[i]; amax=fmax(amax,ai); } return amax; };
    inline double max    (int n, const double* a            ){ double amin=+1e+300; for (int i=0; i<n; i++ ){ double ai=a[i]; amin=fmin(amin,ai); } return amin; };
    inline void   bounds (int n, const double* a, double amin, double amax){ amin=+1e+300; amax=-1e+300; for (int i=0; i<n; i++ ){ double ai=a[i]; amin=fmin(amin,ai); amin=fmax(amax,ai); }; };
    inline double absmax (int n, const double* a            ){ double amax=0;       for (int i=0; i<n; i++ ){ double ai=a[i]; amax=fmax(amax,fabs(ai)); } return amax; };

    inline void minmax(int n, const double* a, double& vmin, double& vmax ){ vmin=+1e+300; vmax=-1e+300; for (int i=0; i<n; i++ ){ double ai=a[i]; vmin=_min(vmin,ai); vmax=_max(vmax,ai); } };

    inline int    imin(int n, const double* a ){ double amin=+1e+300; int im=-1; for (int i=0; i<n; i++ ){ double ai=a[i]; if(ai<amin){amin=ai;im=i;} } return im;  }
    inline int    imax(int n, const double* a ){ double amax=-1e+300; int im=-1; for (int i=0; i<n; i++ ){ double ai=a[i]; if(ai>amax){amax=ai;im=i;} } return im;  }

    inline double err2     (int n, const double* y1s, double* y2s ){ double sum=0;        for (int i=0; i<n; i++ ){ double d=(y1s[i]-y2s[i]); sum+=d*d;          } return sum;  }
    inline double errAbsMax(int n, const double* y1s, double* y2s ){ double amax=-1e+300; for (int i=0; i<n; i++ ){ double d=(y1s[i]-y2s[i]); amax=_max(d,amax); } return amax; }

    // =======  basic vector arritmetics

	inline void set( int n, double  f,            double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f;	      } }
	inline void add( int n, double  f, const double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f+b[i];    } }
	inline void mul( int n, double  f, const double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = f*b[i];    } }

	inline void set( int n, const double* a,            double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i];      } }
	inline void add( int n, const double* a, const double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]+b[i]; } }
	inline void sub( int n, const double* a, const double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]-b[i]; } }
	inline void mul( int n, const double* a, const double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]*b[i]; } }
	inline void div( int n, const double* a, const double* b, double* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]/b[i]; } }
	inline void fma( int n, const double* a, const double* b, double f, double* out ){ for(int i=0; i<n; i++) { out[i]=a[i]+f*b[i]; }  }

	// =======  function-array operations

    inline double dot      (int n, double* xs, double* ys, Func1d func ){ double sum=0;        for (int i=0; i<n; i++ ){ sum+= ys[i]*func(xs[i]); }                         return sum;  }
    inline double err2     (int n, double* xs, double* ys, Func1d func ){ double sum=0;        for (int i=0; i<n; i++ ){ double d=(func(xs[i])-ys[i]); sum+=d*d;          } return sum;  }
    inline double errAbsMax(int n, double* xs, double* ys, Func1d func ){ double amax=-1e+300; for (int i=0; i<n; i++ ){ double d=(func(xs[i])-ys[i]); amax=_max(d,amax); } return amax; }
    inline double min (int n, double* a, Func1d func          ){ double amax=-1e+300; for (int i=0; i<n; i++ ){ double ai=a[i]; amax=_max(amax,ai); } return amax; };

    inline void set( int n, double* xs,             Func1d func, double* out ){ for (int i=0; i<n; i++ ){ out[i] = func(xs[i]);       } };
    //inline void add  ( int n, double* xs, double* ys, Func1d func, double* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]+func(xs[i]); } }
	//inline void sub  ( int n, double* xs, double* ys, Func1d func, double* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]-func(xs[i]); } }
	//inline void mul  ( int n, double* xs, double* ys, Func1d func, double* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]*func(xs[i]); } }
	//inline void div  ( int n, double* xs, double* ys, Func1d func, double* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]/func(xs[i]); } }

    // =======  function-function operations

    inline double dot( int n, double* xs, Func1d funcA, Func1d funcB ){ double sum=0; for (int i=0; i<n; i++ ){ double xi=xs[i]; sum+=funcA(xi)*funcB(xi); } return sum; }
    inline void mdot( int n, int m, const double* A, const double* x, double* Ax ){
        for(int i=0; i<m; i++){
            Ax[i] = VecN::dot( n, A, x );
            A+=n;
        }
    }
	// =======  initialization and I/O

	inline void arange ( int n, double xmin, double dx  , double* out ){ double x=xmin; for(int i=0; i<n; i++){out[i]=x; x+=dx; } }
	inline void linspan( int n, double xmin, double xmax, double* out ){ double dx=(xmax-xmin)/n; arange( n, xmin, dx,out );    }

	inline void random_vector ( int n, double xmin, double xmax, double * out ){
		double xrange = xmax - xmin;
		for (int i=0; i<n; i++ ){		out[i] = xmin + xrange*randf();	}
	}

	inline void print_vector( int n, double * a ){
		for (int i=0; i<n; i++ ){	printf( "%f ", a[i] );	}
		printf( "\n" );
	}



 // =======  iterration over array
	template<typename T> inline T norm2 (int n,const T* a                 ){ T sum =0; for (int i=0; i<n; i++ ){ T ai=a[i]; sum+=ai*ai;          } return sum; }
	template<typename T>inline T wnorm2(int n,const T* a,const T* w ){ T sum =0; for (int i=0; i<n; i++ ){ T ai=a[i]; sum+= (ai*ai)*w[i];  } return sum; }
	template<typename T>inline T wdot  (int n,const T* a,const T* b,const T* w){ T sum =0; for (int i=0; i<n; i++ ){ sum+= a[i]*b[i]*w[i];  } return sum; }

	template<typename T>inline T dot    (int n, const T* a, const T* b ){ T sum =0; for (int i=0; i<n; i++ ){ sum+= a[i]*b[i];             } return sum; }
	template<typename T>inline T sum    (int n, const T* a            ){ T sum =0; for (int i=0; i<n; i++ ){ sum+= a[i];                  } return sum; };
    template<typename T>inline T sum2   (int n, const T* a            ){ T sum =0; for (int i=0; i<n; i++ ){ T ai=a[i]; sum+=ai*ai;  } return sum; };
    template<typename T>inline T min    (int n, const T* a            ){ T amax=-1e+300; for (int i=0; i<n; i++ ){ T ai=a[i]; amax=fmax(amax,ai); } return amax; };
    template<typename T>inline T max    (int n, const T* a            ){ T amin=+1e+300; for (int i=0; i<n; i++ ){ T ai=a[i]; amin=fmin(amin,ai); } return amin; };
    template<typename T>inline void   bounds (int n, const T* a, T amin, T amax){ amin=+1e+300; amax=-1e+300; for (int i=0; i<n; i++ ){ T ai=a[i]; amin=fmin(amin,ai); amin=fmax(amax,ai); }; };
    template<typename T>inline T absmax (int n, const T* a            ){ T amax=0;       for (int i=0; i<n; i++ ){ T ai=a[i]; amax=fmax(amax,fabs(ai)); } return amax; };

    template<typename T>inline void minmax(int n, const T* a, T& vmin, T& vmax ){ vmin=+1e+300; vmax=-1e+300; for (int i=0; i<n; i++ ){ T ai=a[i]; vmin=_min(vmin,ai); vmax=_max(vmax,ai); } };

    template<typename T>inline int    imin(int n, const T* a ){ T amin=+1e+300; int im=-1; for (int i=0; i<n; i++ ){ T ai=a[i]; if(ai<amin){amin=ai;im=i;} } return im;  }
    template<typename T>inline int    imax(int n, const T* a ){ T amax=-1e+300; int im=-1; for (int i=0; i<n; i++ ){ T ai=a[i]; if(ai>amax){amax=ai;im=i;} } return im;  }

    template<typename T>inline T err2     (int n, const T* y1s, T* y2s ){ T sum=0;        for (int i=0; i<n; i++ ){ T d=(y1s[i]-y2s[i]); sum+=d*d;          } return sum;  }
    template<typename T>inline T errAbsMax(int n, const T* y1s, T* y2s ){ T amax=-1e+300; for (int i=0; i<n; i++ ){ T d=(y1s[i]-y2s[i]); amax=_max(d,amax); } return amax; }

    // =======  basic vector arritmetics

	template<typename T>inline void set( int n, T  f,            T* out ){  	for (int i=0; i<n; i++ ){ out[i] = f;	      } }
	template<typename T>inline void add( int n, T  f, const T* b, T* out ){  	for (int i=0; i<n; i++ ){ out[i] = f+b[i];    } }
	template<typename T>inline void mul( int n, T  f, const T* b, T* out ){  	for (int i=0; i<n; i++ ){ out[i] = f*b[i];    } }

	template<typename T>inline void set( int n, const T* a,            T* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i];      } }
	template<typename T>inline void add( int n, const T* a, const T* b, T* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]+b[i]; } }
	template<typename T>inline void sub( int n, const T* a, const T* b, T* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]-b[i]; } }
	template<typename T>inline void mul( int n, const T* a, const T* b, T* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]*b[i]; } }
	template<typename T>inline void div( int n, const T* a, const T* b, T* out ){  	for (int i=0; i<n; i++ ){ out[i] = a[i]/b[i]; } }
	template<typename T>inline void fma( int n, const T* a, const T* b, T f, T* out ){ for(int i=0; i<n; i++) { out[i]=a[i]+f*b[i]; }  }

	// =======  function-array operations

    template<typename T>inline T dot      (int n, T* xs, T* ys, Func1d func ){ T sum=0;        for (int i=0; i<n; i++ ){ sum+= ys[i]*func(xs[i]); }                         return sum;  }
    template<typename T>inline T err2     (int n, T* xs, T* ys, Func1d func ){ T sum=0;        for (int i=0; i<n; i++ ){ T d=(func(xs[i])-ys[i]); sum+=d*d;          } return sum;  }
    template<typename T>inline T errAbsMax(int n, T* xs, T* ys, Func1d func ){ T amax=-1e+300; for (int i=0; i<n; i++ ){ T d=(func(xs[i])-ys[i]); amax=_max(d,amax); } return amax; }
    template<typename T>inline T min (int n, T* a, Func1d func          ){ T amax=-1e+300; for (int i=0; i<n; i++ ){ T ai=a[i]; amax=_max(amax,ai); } return amax; };

    template<typename T>inline void set( int n, T* xs,             Func1d func, T* out ){ for (int i=0; i<n; i++ ){ out[i] = func(xs[i]);       } };
    //template<typename T>inline void add  ( int n, T* xs, T* ys, Func1d func, T* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]+func(xs[i]); } }
	//template<typename T>inline void sub  ( int n, T* xs, T* ys, Func1d func, T* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]-func(xs[i]); } }
	//template<typename T>inline void mul  ( int n, T* xs, T* ys, Func1d func, T* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]*func(xs[i]); } }
	//template<typename T>inline void div  ( int n, T* xs, T* ys, Func1d func, T* out ){ for (int i=0; i<n; i++ ){ out[i] = ys[i]/func(xs[i]); } }

    // =======  function-function operations

    template<typename T>inline T dot( int n, T* xs, Func1d funcA, Func1d funcB ){ T sum=0; for (int i=0; i<n; i++ ){ T xi=xs[i]; sum+=funcA(xi)*funcB(xi); } return sum; }
    template<typename T>inline void mdot( int n, int m, const T* A, const T* x, T* Ax ){
        for(int i=0; i<m; i++){
            Ax[i] = VecN::dot( n, A, x );
            A+=n;
        }
    }
	// =======  initialization and I/O

	template<typename T>inline void arange ( int n, T xmin, T dx  , T* out ){ T x=xmin; for(int i=0; i<n; i++){out[i]=x; x+=dx; } }
	template<typename T>inline void linspan( int n, T xmin, T xmax, T* out ){ T dx=(xmax-xmin)/n; arange( n, xmin, dx,out );    }

	template<typename T>inline void random_vector ( int n, T xmin, T xmax, T * out ){
		T xrange = xmax - xmin;
		for (int i=0; i<n; i++ ){		out[i] = xmin + xrange*randf();	}
	}

	template<typename T>inline void print_vector( int n, T * a ){
		for (int i=0; i<n; i++ ){	printf( "%f ", a[i] );	}
		printf( "\n" );
	}


}



#endif
