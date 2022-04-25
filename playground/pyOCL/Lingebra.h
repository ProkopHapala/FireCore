
#ifndef  Lingebra_h
#define  Lingebra_h

#include <math.h>
#include <cstdlib>
#include <stdio.h>

#include "fastmath.h"
#include "VecN.h"


void (dotFunc)(int nx, int ny, double * x, double * y );


//#define IJ(i,j) n*i+j

namespace Lingebra{

	// ==== function declarations

    //double ** from_continuous( int m, int n, double *p );
    double ** from_continuous( int m, int n, double *p, double ** A=0 );
	double ** new_matrix     ( int m, int n );
	void      delete_matrix  ( int m,        double** A );
	void transpose           ( int m, int n, double** A, double** TA );
	void dot                 ( int m, int n, double** A, double* x, double* out );
	void dotT                ( int m, int n, double** A, double* x, double* out );

    void symCopy             ( int m, double* A, bool toLower );


	void mmul_ik_kj( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ik_jk( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ki_kj( int ni, int nj, int nk, double** A, double** B, double** out );
	void mmul_ki_jk( int ni, int nj, int nk, double** A, double** B, double** out );

	void random_matrix( int m, int n, double xmin, double xmax, double** out );
	void print_matrix( int m, int n, double ** A );

	void   makeQuadricFormMatrix( int m, int n, double * ks, double ** A, double ** Q );
	double evalQudraticForm( int n, double* x, double** Q );
	double evalQudraticFormDirs( int m, int n, double* x, double* k, double** A );

	void GaussElimination    ( int n, double ** A, double * c, int * index );
	void linSolve_gauss      ( int n, double ** A, double * b, int * index, double * x );
	void linSolve_CG         ( int n, double ** A, double * b, double * x );
	void linSolve_BCG        ( int n, int m, double ** A, double * b, double * x );
	void leastSquareFit_Gauss( int n, int m, double ** A, double * b, double * x );

	int eig_Jacobi( int n, double* A, double* V, double* es, double tol, int nMaxIter  );
	double eig_Jacobi_init( int n, double* A, double* V, int* mjs, int& imax, int& jmax );
	void eig_Jacobi_step  ( int n, double* A, double* V, int* mjs, int& imax, int& jmax, double& vmax );

	// ==== inline function

	inline void set( int m, int n, double   f,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, f,       out[i] );  } }
	inline void add( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, f, B[i], out[i] );  } }
	inline void mul( int m, int n, double   f, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, f, B[i], out[i] );  } }

	inline void set( int m, int n, double** A,             double** out ){  for (int i=0; i<m; i++ ){ VecN::set( n, A[i],       out[i] );   } }
	inline void add( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::add( n, A[i], B[i], out[i] );   } }
	inline void sub( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::sub( n, A[i], B[i], out[i] );   } }
	inline void mul( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::mul( n, A[i], B[i], out[i] );   } }
	inline void div( int m, int n, double** A, double** B, double** out ){  for (int i=0; i<m; i++ ){ VecN::div( n, A[i], B[i], out[i] );   } }

	//template<void dotFunc(int nx, int ny, double* x, double* y)>
	template<typename F>
    void genLinSolve_CG( int n, double * b, double * x , F dotFunc ){
        const int    maxIters   = 10;
        const double    maxErr2 = 1e-5;
        double *  r     = new double[n];
        double *  r2    = new double[n];
        double *  p     = new double[n];
        double *  Ap    = new double[n];
        dotFunc( n, n, x, r );
        VecN::sub( n, b, r, r );
        VecN::set( n, r, p );
        double rho = VecN::dot(n, r,r);
        double alpha = 0;
        for ( int i =0; i<maxIters; i++) {
            dotFunc( n, n, p, Ap);
            alpha = rho / VecN::dot(n, p, Ap);
            VecN::fma( n, x, p ,  alpha,   x );
            VecN::fma( n, r, Ap, -alpha,   r2 );
            double err2 = VecN::dot(n, r2,r2);
            //printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
            printf( " iter: %i  err2: %f \n", i, err2 );
            if (err2 < maxErr2 ) break;
            double rho2 = VecN::dot(n, r2,r2);
            double beta = rho2 / rho;
            VecN::fma( n, r2, p, beta, p );
            rho = rho2;
            double * swap = r; r = r2; r2 = swap;
        }
        delete [] r;
        delete [] r2;
        delete [] p;
        delete [] Ap;
    }

};


//template<void dotFunc(int nx, int ny, double* x, double* y)>
class LinSolver{ public:
    //static constexpr int       maxIters = 10;
    //static constexpr double    maxErr2  = 1e-5;
    int n;
    // temp
    int    istep=0;
    double *  r=0;
    double *  r2=0;
    double *  p=0;
    double *  Ap=0;
    double rho = 0;
    double alpha = 0;

    // to solve
    double* x = 0;
    double* b = 0;
    double* M = 0;

    virtual void dotFunc( int n, double * x, double * Ax )=0;

    void realloc(int n_){
        if(n_!=n){
            n = n_;
            _realloc(r,n);
            _realloc(r2,n);
            _realloc(p,n);
            _realloc(Ap,n);
        }
    };

    void dealloc(){
        delete [] r;
        delete [] r2;
        delete [] p;
        delete [] Ap;
    };

    void setLinearProblem(int n_, double* x_, double* b_, double* M_ = 0 ){
        realloc(n_);
        x=x_; b=b_; M=M_;
    }

    double step_GD(double dt){
        dotFunc  ( n, x, r );
        VecN::sub( n, b, r, r );
        //VecN::add( n, b, r, r );
        VecN::fma( n, x, r, dt, x );
        return VecN::dot(n, r,r);
    }

    double step_CG(){
        // see https://en.wikipedia.org/wiki/Conjugate_gradient_method
        //printf( "LinSolver::step_CG %i \n", istep );
        if(istep==0){
            dotFunc  ( n, x, r );
            //printf("r   "); VecN::print_vector(n, r);
            VecN::sub( n, b, r, r ); // r = b - A*x
            //printf("r_  "); VecN::print_vector(n, r);
            VecN::set( n, r, p );    // p = r
            rho = VecN::dot(n, r,r);
            alpha = 0;
            //printf( "rho %f alpha %f \n", rho, alpha );
        }else{
            double rho2 = VecN::dot(n, r2,r2);
            double beta = rho2 / rho;
            VecN::fma( n, r2, p, beta, p );
            rho = rho2;
            double * tmp = r; r = r2; r2 = tmp;
        }
        // NOTE : BCQ can be done if (A.T()*A) is applied instead of A in dotFunc
        //printf("p  "); VecN::print_vector(n, p);
        dotFunc( n, p, Ap);
        //printf("Ap "); VecN::print_vector(n, Ap);
        alpha = rho / VecN::dot(n, p, Ap);    // a  = <r|r>/<p|A|p>
        //printf( "rho %f alpha %f \n", rho, alpha );
        VecN::fma( n, x, p ,  alpha,   x );   // x  = x - a*p
        VecN::fma( n, r, Ap, -alpha,   r2 );  // r2 = r - a*A|p>
        double err2 = VecN::dot(n, r2,r2);
        istep++;
        return err2;
        //printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
        //printf( " iter: %i  err2: %f \n", i, err2 );
        //if (err2 < maxErr2 ) break;
    }

    void solve_CG( int maxIters, double maxErr2 ){
        istep = 0;
        for ( int i =0; i<maxIters; i++) {
            if ( step_CG() > maxErr2 ) break;
        }
    }

};





// ======================= IMPLEMENTATION
// ======================= IMPLEMENTATION
// ======================= IMPLEMENTATION
// ======================= IMPLEMENTATION
// ======================= IMPLEMENTATION





#include "Lingebra.h" // THE HEADER

/*
Discussion:

Should matrix be stored as array of pointers A[i][j] or as plain array A[i+j*n] ?
 - A[i][j] is more general (we may easily select arbitrary sub-matrix without copy)
 - A[i+j*n] could be faster (?) but if n is large cache coherency is bad anyway
*/

namespace Lingebra {

// creates double** from any continuous memory block folowing *p
//double ** from_continuous( int m, int n, double *p, double ** A=0 ){
double ** from_continuous( int m, int n, double *p, double ** A ){
    if(A==0) A = new double*[m];
    for (int i=0; i<m; i++ ){ A[i] = &(p[i*n]); }
    return A;
}

/*
void from_continuous( int m, int n, double *p, double** A ){
    //double ** A = new double*[m];
    for (int i=0; i<m; i++ ){ A[i] = &(p[i*n]); }
    return A;
}
*/

double ** new_matrix( int m, int n ){
    double ** A = new double*[m];
    for (int i=0; i<m; i++ ){ A[i] = new double[n]; }
    return A;
}

void delete_matrix( int m, double** A ){
    for (int i=0; i<m; i++ ){ delete A[i]; }
    delete [] A;
}

void transpose( int m, int n, double** A, double** TA ){
    for (int i=0; i<m; i++ ){
        for (int j=0; j<n; j++ ){
            TA[i][j] = A[j][i];
        }
    }
}

void symCopy ( int n, double* A, bool toLower ){
    if ( toLower ){ for (int i=0;i<n;i++){ for (int j=0;j<i;j++){ A[i*n+j] = A[j*n+i]; } } }
    else          { for (int i=0;i<n;i++){ for (int j=0;j<i;j++){ A[j*n+i] = A[i*n+j]; } } }
};

void dot( int m, int n, double** A, double* x, double* out ){
    for (int i=0; i<m; i++ ){
        out[i] = VecN::dot( n, A[i], x );
    }
}

void dotT( int m, int n, double** A, double* x, double* out ){
    for (int i=0; i<m; i++ ){
        double doti = 0;
        for (int j=0; j<n; j++ ){ doti += A[j][i] * x[j];	}
        out[i] = doti;
    }
}


void mmul_ik_kj( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[i][k] * B[k][j];	}
			out[i][j] = dotij;
		}
	}
}

void mmul_ik_jk( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[i][k] * B[j][k];	}
			out[i][j] = dotij;
		}
	}
}

void mmul_ki_kj( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[k][i] * B[k][j];	}
			out[i][j] = dotij;
		}
	}
}

void mmul_ki_jk( int ni, int nj, int nk, double** A, double** B, double** out ){
	for(int i=0; i<ni; i++){
		for(int j=0; j<nj; j++){
			double dotij = 0;
			for(int k=0; k<nk; k++){	dotij += A[k][i] * B[j][k];	}
			out[i][j] = dotij;
		}
	}
}

void random_matrix( int m, int n, double xmin, double xmax, double** out ){
	//double xrange = xmax - xmin;
	for (int i=0; i<m; i++ ){ VecN::random_vector ( n, xmin, xmax, out[i] ); }
}

void print_matrix( int m, int n, double ** A ){
	for (int i=0; i<m; i++ ){ VecN::print_vector( n, A[i] );	}
}

// ===============================================
// ======                                  =======
// ======              Qadric              =======
// ======                                  =======
// ===============================================

//  Q =  Sum_i |a_i> k_i <a_i|
void makeQuadricFormMatrix( int m, int n, double * ks, double ** A, double ** Q ){
	for (int im=0; im<m; im++ ){
		double kim   = ks[im];
		double * Aim =  A[im];
		for (int i=0; i<n; i++ ){
			for (int j=0; j<i; j++ ){
				double Qij = kim*Aim[i]*Aim[j];
				Q[i][j] += Qij;
				Q[j][i] += Qij;
			}
			Q[i][i] += kim*Aim[i]*Aim[i];
		}
	}
}

//  y =  < x | Q.x >
double evalQudraticForm( int n, double* x, double** Q ){
	double y =  0;
	for (int i=0; i<n; i++ ){	y += x[i] * VecN::dot( n, Q[i], x );   }
	return y;
}

//  y = Sum_i k_i < a_i | x >^2
double evalQudraticFormDirs( int m, int n, double* x, double* k, double** A ){
	double y =  0;
	for (int i=0; i<m; i++ ){
		double ri = VecN::dot( n, A[i], x );
		y += k[i]*ri*ri;
	}
	return y;
}


// ===============================================
// ======                                  =======
// ======       Gauss Elimination          =======
// ======                                  =======
// ======     Linsolve, invertMatrix       =======
// ======                                  =======
// ===============================================


// GaussElimination
// Method to carry out the partial-pivoting Gaussian
// elimination.  Here index[] stores pivoting order.

void GaussElimination( int n, double ** A, double * c, int * index ) {

	// Initialize the index
	for (int i=0; i<n; ++i) index[i] = i;

	// Find the rescaling factors, one from each row
	for (int i=0; i<n; ++i) {
	  double c1 = 0;
	  for (int j=0; j<n; ++j) {
		double c0 = fabs( A[i][j] );
		if (c0 > c1) c1 = c0;
	  }
	  c[i] = c1;
	}

	// Search the pivoting element from each column
	int k = 0;
	for (int j=0; j<n-1; ++j) {
	  double pi1 = 0;
	  for (int i=j; i<n; ++i) {
		double pi0 = fabs( A[ index[i] ][j] );
		pi0 /= c[ index[i] ];
		if (pi0 > pi1) {
		  pi1 = pi0;
		  k = i;
		}
	  }

	  // Interchange rows according to the pivoting order
	  int itmp = index[j];
	  index[j] = index[k];
	  index[k] = itmp;
	  for (int i=j+1; i<n; ++i) {
		double pj = A[ index[i] ][ j ]/A[ index[j] ][j];

	   // Record pivoting ratios below the diagonal
		A[ index[i] ][j] = pj;

	   // Modify other elements accordingly
		for (int l=j+1; l<n; ++l)
		  A[ index[i] ][l] -= pj*A[ index[j] ][l];
	  }
	}
}

void linSolve_gauss( int n, double ** A, double * b, int * index, double * x ) {

	// Transform the matrix into an upper triangle
	GaussElimination( n, A, x, index);

	// Update the array b[i] with the ratios stored
	for(int i=0; i<n-1; ++i) {
		for(int j =i+1; j<n; ++j) {
			b[index[j]] -= A[index[j]][i]*b[index[i]];
		}
	}

	// Perform backward substitutions
	x[n-1] = b[index[n-1]]/A[index[n-1]][n-1];
	for (int i=n-2; i>=0; --i) {
		x[i] = b[index[i]];
		for (int j=i+1; j<n; ++j) {
			x[i] -= A[index[i]][j]*x[j];
		}
		x[i] /= A[index[i]][i];
	}
}

void linSolve_CG( int n, double ** A, double * b, double * x ){
	const int maxIters   = 10;
	const double maxErr2 = 1e-5;
	double *  r     = new double[n];
	double *  r2    = new double[n];
	double *  p     = new double[n];
	double *  Ap    = new double[n];
	dot( n, n, A, x, r );
	VecN::sub( n, b, r, r );
	VecN::set( n, r, p );
	double rho = VecN::dot(n, r,r);
	double alpha = 0;
	for ( int i =0; i<maxIters; i++) {
		dot( n, n, A, p, Ap);
		alpha = rho / VecN::dot(n, p, Ap);
		VecN::fma( n, x, p ,  alpha,   x );
		VecN::fma( n, r, Ap, -alpha,   r2 );
		double err2 = VecN::dot(n, r2,r2);
		//printf( " iter: %i  err2: %f |  alpha %f \n", i, err2,     alpha );
		printf( " iter: %i  err2: %f \n", i, err2 );
		if (err2 < maxErr2 ) break;
		double rho2 = VecN::dot(n, r2,r2);
		double beta = rho2 / rho;
		VecN::fma( n, r2, p, beta, p );
		rho = rho2;
		double * swap = r; r = r2; r2 = swap;
	}
	delete r;
	delete r2;
	delete p;
	delete Ap;
}

void linSolve_BCG( int n, int m, double ** A, double * b, double * x ){
	const int maxIters   = 10;
	const double maxErr2 = 1e-5;
	double *  r     = new double[n];
	double *  r2    = new double[n];
	double *  p     = new double[n];
	double *  Ap    = new double[n];
	dot ( n, m, A, x, r );
	dotT( m, n, A, x, r );
	VecN::sub( n, b, r, r );
	VecN::set( n, r, p );
	double rho = VecN::dot(n, r,r);
	double alpha = 0;
	for ( int i =0; i<maxIters; i++) {
		dot ( n, m, A, p, Ap);
		dotT( m, n, A, p, Ap);
		alpha = rho / VecN::dot(n, p, Ap);
		VecN::fma( n, x, p ,  alpha,   x );
		VecN::fma( n, r, Ap, -alpha,   r2 );
		double err2 = VecN::dot(n, r2,r2);
		printf( " iter: %i  err2: %f \n", i, err2 );
		if (err2 < maxErr2 ) break;
		double rho2 = VecN::dot(n, r2,r2);
		double beta = rho2 / rho;
		VecN::fma( n, r2, p, beta, p );
		rho = rho2;
		double * swap = r; r = r2; r2 = swap;
	}
	delete r;
	delete r2;
	delete p;
	delete Ap;
}

void leastSquareFit_Gauss( int n, int m, double ** A, double * b, double * x ){
	double ** AA    = new_matrix( n, n );
	double *  Ab    = new double[n];
	int    *  index = new int   [n];
	dotT( n, m, A, b, Ab );
	mmul_ki_kj( n, n, m, A, A, AA );
	//linSolve_gauss( n, AA, Ab, index, x );
	linSolve_CG( n, AA, Ab, x );
	VecN::print_vector( n, x );
	delete AA;
	delete Ab;
	delete index;
}



// ===============================================
// ======       EigenValueSolver          =======
// ===============================================


void get_diag_vector( int n, double* a, double* v ){
  for ( int i = 0; i < n; i++ ){ v[i] = a[i*n+i]; }
}

void fill_identity ( int n, double* a ){
  int k = 0;
  for (int j = 0; j < n; j++ ){
    for (int i = 0; i < n; i++ ){
      if ( i == j ){ a[k] = 1.0; }else{ a[k] = 0.0; }
      k = k + 1;
    }
  }
}

inline void jacobi_rot( double* g_, double* h_, double s, double tau ){
    double g=*g_; double h=*h_;
    *g_ = g - s * ( h + g * tau );
    *h_ = h + s * ( g - h * tau );
}


inline void jacobi_rot_cs( double* g_, double* h_, double c, double s ){
    double g=*g_; double h=*h_;
    *g_ = c*g - s*h;
    *h_ = s*g + c*h;
}

void jacobi_rotation( int n, double* A, double* V, int k, int l ){
    //  https://en.wikipedia.org/wiki/Jacobi_rotation
    double aDiff = A[l*n+l] - A[k*n+k];
    double akl   = A[k*n+l];
    double t;
    if ( fabs(akl) < fabs(aDiff)*1.0e-36d ){
        t = akl/aDiff;
    }else{
        double phi = aDiff/(2*akl);
        t = 1/(fabs(phi) + sqrt(phi*phi + 1));
        if (phi < 0.0) t = -t;
    }
    double c   = 1/sqrt(t*t + 1);
    double s   = t*c;
    double tau = s/(1+c);
    A[k*n+l] = 0;
    A[k*n+k] = A[k*n+k] - t*akl;
    A[l*n+l] = A[l*n+l] + t*akl;
    for (int i = 0    ; i < k; i++ ){ jacobi_rot( A+(i*n+k), A+(i*n+l), s, tau ); }
    for (int i = k + 1; i < l; i++ ){ jacobi_rot( A+(k*n+i), A+(i*n+l), s, tau ); }
    for (int i = l + 1; i < n; i++ ){ jacobi_rot( A+(k*n+i), A+(l*n+i), s, tau ); }
    for (int i = 0    ; i < n; i++ ){ jacobi_rot( V+(i*n+k), V+(i*n+l), s, tau ); }
}

void jacobi_rotation_small( int n, double* A, double* V, int k, int l ){
    //  https://en.wikipedia.org/wiki/Jacobi_rotation
    // this version is optimized for small matrices (like 3x3, 4x4 etc. )
    const int kn=k*n;
    const int ln=l*n;
    double * Ak = A+kn;
    double * Al = A+ln;
    double aDiff = Al[l] - Ak[k];
    double akl   = Ak[l];
    double t;
    if ( fabs(akl) < fabs(aDiff)*1.0e-36d ){
        t = akl/aDiff;
    }else{
        double phi = aDiff/(2*akl);
        t = 1/(fabs(phi) + sqrt(phi*phi + 1));
        if (phi < 0.0) t = -t;
    }
    double c    = 1/sqrt(t*t + 1);
    double s    = t*c;
    double tau  = s/(1+c);
    double takl = t*akl;
    Ak[l] = 0;
    Ak[k] -= takl;
    Al[l] += takl;
    for (int in = 0  ; in<kn; in+=n ){ jacobi_rot( A+(in+k),  A+(in +l), s, tau ); }
    for (int i  = k+1; i <l;  i++   ){ jacobi_rot( Ak+i,      A+(i*n+l), s, tau ); } // is the same as  Al+i
    for (int i  = l+1; i <n;  i++   ){ jacobi_rot( Ak+i,      Al+i,      s, tau ); }
    //for (int i = 0    ; i < n; i++ ){ jacobi_rot( V+(i*n+k), V+(i*n+l), s, tau ); }
    double * Vk = V+kn;
    double * Vl = V+ln;
    for (int i = 0    ; i < n; i++ ){ jacobi_rot( Vk+i, Vl+i, s, tau ); }   // Much faster if the matrix is transposed

}









template<double func(double x)>
inline int getMaxIndex(int n, double* v ){
    int imax    = 0;
    double vmax = func(v[0]);
    for(int i=1;i<n;i++){
        double vi = func(v[i]);
        //printf(">>>>> %i %2.5f %2.5f \n", i, vi, vmax);
        if(vi>vmax){ vmax=vi; imax=i; }
    }
    return imax;
}

double updateRowMax(int n, double* A, int* mjs, int& imax, int& jmax){
    double vmax=0;
    int    imax_=imax;
    int    jmax_=jmax;
    for(int i=0; i<(n-1); i++){
        int     mj = mjs[i];
        double vmj = fabs( A[i*n+mj] );

        if ((i==imax)||(i==jmax)||(mj==imax)||(mj==jmax)){
            int j=i+1;
            mj   = getMaxIndex<fabs>( (n-j), A+(i*n+j) )+j;
            vmj  = fabs( A[i*n+mj] );
            //printf("%i case row \n", i, mj, vmj);
        }else{
            if (imax>i){
                double v=fabs( A[i*n+imax] );
                if(v>vmj){ mj=imax; vmj=v; }
                //printf("%i case imax %i %f\n", i, mj, vmj );
            }
            if (jmax>i){
                double v=fabs(A[i*n+jmax]);
                if(v>vmj){ mj=jmax; vmj=v; }   // //vmj=v is not required
                //printf("%i case jmax %i %f\n", i, mj, vmj );
            }
        }
        mjs[i]=mj;

        if(vmj>vmax){
            vmax=vmj;
            imax_=i;
            jmax_=mj;
        }
        //printf("mj = %i \n", mj );
    }
    imax=imax_; jmax=jmax_;
    //printf("imax,jmax (%i,%i)\n", imax,jmax );
    return vmax;
}

void printmatrix( int ni, int nj, double* A, char* format ){
    for(int i=0; i<ni; i++){
        for(int j=0; j<nj; j++){
            printf( format, A[i*nj+j]);
        }
        printf("\n");
    }
}

double eig_Jacobi_init( int n, double* A, double* V, int* mjs, int& imax, int& jmax ){
    fill_identity( n, V );
    double vmax=0;
    for(int i=0;i<(n-1);i++){
        int j=i+1;
        int mj = getMaxIndex<fabs>( (n-j), A+(i*n+j) )+j;
        mjs[i] = mj;
        double v=fabs(A[i*n+mj]);
        //printf(">> %i %i %f\n", i, mj, v );
        if(v>vmax){ vmax=v; imax=i; jmax=mj; }
    }
    return vmax;
}

void eig_Jacobi_step( int n, double* A, double* V, int* mjs, int& imax, int& jmax, double& vmax ){
    jacobi_rotation    ( n, A, V,   imax, jmax );
    vmax = updateRowMax( n, A, mjs, imax, jmax );
    /*
    vmax=0;
    for(int i=0;i<(n-1);i++){
        int j=i+1;
        int mj = getMaxIndex<fabs>( (n-j), A+(i*n+j) )+j;
        mjs[i] = mj;
        double v=fabs(A[i*n+mj]);
        //printf(">> %i %i %f\n", i, mj, v );
        if(v>vmax){ vmax=v; imax=i; jmax=mj; }
    }
    */
}

int eig_Jacobi( int n, double* A, double* V, double* es, double tol, int nMaxIter  ){
    fill_identity( n, V );
    int* mjs = new int[n];
    double vmax=0;
    int    imax,jmax;
    for(int i=0;i<(n-1);i++){
        int j=i+1;
        int mj = getMaxIndex<fabs>( (n-j), A+(i*n+j) )+j;
        mjs[i] = mj;
        double v=fabs(A[i*n+mj]);
        //printf(">> %i %i %f\n", i, mj, v );
        if(v>vmax){ vmax=v; imax=i; jmax=mj; }
    }
    int iter=0;
    for(int iter=0; iter<nMaxIter; iter++){
        //printf("%i (%i,%i) %f\n", iter, imax,jmax,vmax);
        //printmatrix( n,n, A, " %2.3f" );
        jacobi_rotation    ( n, A, V,   imax, jmax );
        vmax = updateRowMax( n, A, mjs, imax, jmax );
        if( vmax < tol ) break;
    }

    if( vmax < tol ){  printf( "converged by %i rotations; error < %e \n", iter, tol ); }
    else            {  printf( "not converged in %i rotations; A[%i,%i] = %e \n", iter, imax,jmax, vmax  );  }

    for(int i=0;i<n;i++){ es[i] = A[i*n+i]; };
    delete mjs;
    return iter;
}



void aproxOrthoNormStep(int nv, int m, double * M ){
    for(int i=0;i<nv;i++){
        double* vi = M+i*m;
        for(int j=0;j<i;j++){
            double* vj = M+j*m;
            double cij = 0;
            for(int k=0;k<m;k++){ cij += vi[k]*vj[k]; }
            // approx half angle
            // cos(a/2) = sqrt((cos(a)+1)/2) = sqrt( 1 +   (cos(a)-1)/2 )  =  sqrt( 1 +   x )
            double x  = (cij-1)*0.5;
            // taylor for sqrt( 1 +   x )
            double ch = -(1 + x*(0.5 + x*(-0.125 + x*0.0625)));
            for(int k=0;k<m;k++){
                double vik = vi[k];
                vi[k]+=vj[k]*ch;
                vj[k]+=vik  *ch;
            }
        }
    }
}


void orthtoForce(int nv, int m, double * P, double * F ){
    for(int i=0;i<nv;i++){
        double* pi = P+i*m;
        double* fi = F+i*m;
        double  cii = 0;
        for(int k=0;k<m;k++){ cii += pi[k]*pi[k]; }
        // approx half angle
        // cos(a/2) = sqrt((cos(a)+1)/2) = sqrt( 1 +   (cos(a)-1)/2 )  =  sqrt( 1 +   x )
        double x  = (cii-1);
        // taylor for sqrt( 1 +   x )
        double ch = -x*( -0.5d + x*( 0.375d + x*-0.3125d ) );
        for(int k=0;k<m;k++){ fi[k]+=pi[k]*ch; }
        for(int j=0;j<i;j++){
            double* pj = P+j*m;
            double* fj = F+j*m;
            double cij = 0;
            for(int k=0;k<m;k++){ cij += pi[k]*pj[k]; }
            // approx half angle
            // cos(a/2) = sqrt((cos(a)+1)/2) = sqrt( 1 +   (cos(a)-1)/2 )  =  sqrt( 1 +   x )
            double x  = (cij-1)*0.5;
            // taylor for sqrt( 1 +   x )
            double ch = -(1 + x*(0.5 + x*(-0.125 + x*0.0625)));
            for(int k=0;k<m;k++){ fi[k]+=pj[k]*ch; }
            for(int k=0;k<m;k++){ fj[k]+=pi[k]*ch; }
        }
    }
}

};















#endif

