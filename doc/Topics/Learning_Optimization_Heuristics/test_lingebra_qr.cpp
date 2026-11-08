// Unit test for Lingebra::leastSquareFit_QR zero/duplicate/rank-deficient column handling (§11.2C).
// Build: g++ -O2 -I../../../cpp/common test_lingebra_qr.cpp -o test_lingebra_qr && ./test_lingebra_qr
#include <cstdio>
#include <cmath>
#include "../../../cpp/common/math/Lingebra.h"

int nFail = 0;
void check( bool ok, const char* msg ){ if(!ok){ nFail++; printf("FAIL %s\n", msg); } else printf("ok   %s\n", msg); }
bool finite( const double* v, int n ){ for(int i=0;i<n;i++) if(!std::isfinite(v[i])) return false; return true; }

int main(){
    const int n=8, m=4;
    double A[n*m], b[n], x[m], Q[n*m], R[m*m], vw[n]; int piv[m];
    for(int i=0;i<n;i++)b[i]=0;

    // --- test 1: one all-zero column (the reported NaN case) ---
    for(int j=0;j<m;j++)for(int i=0;i<n;i++)A[j*n+i]=(j==1)?0.0:((i==j)?2.0:0.1*(i+j));
    for(int i=0;i<n;i++)b[i]=2.0*A[i]+1.0;
    int kc = Lingebra::leastSquareFit_QR(n,m,A,b,x,Q,R,piv,vw,1e-8);
    check( finite(x,m),          "zero column: x finite" );
    check( finite(Q,n*m) && finite(R,m*m), "zero column: Q,R finite" );
    check( x[1]==0.0,            "zero column: its coefficient stays 0" );
    printf("      rank=%i x=(%g %g %g %g)\n", kc, x[0],x[1],x[2],x[3]);

    // --- test 2: duplicate column (rank deficient but nonzero) ---
    for(int j=0;j<m;j++)for(int i=0;i<n;i++)A[j*n+i]=(j==2)?A[0*n+i]:((i==j)?2.0:0.1*(i+j));
    for(int i=0;i<n;i++)b[i]=2.0*A[i]+1.0;
    kc = Lingebra::leastSquareFit_QR(n,m,A,b,x,Q,R,piv,vw,1e-8);
    check( finite(x,m),          "duplicate column: x finite" );
    printf("      rank=%i (expect <=3 of 4) x=(%g %g %g %g)\n", kc, x[0],x[1],x[2],x[3]);
    check( kc<=3,                "duplicate column: dropped" );

    // --- test 3: all columns zero (extreme degenerate) ---
    for(int j=0;j<n*m;j++)A[j]=0;
    kc = Lingebra::leastSquareFit_QR(n,m,A,b,x,Q,R,piv,vw,1e-8);
    check( finite(x,m) && kc==0, "all-zero matrix: rank 0, x=0 finite" );

    // --- test 4: well-conditioned solve parity vs known answer ---
    for(int j=0;j<m;j++)for(int i=0;i<n;i++)A[j*n+i]=(i==j)?1.0:0.05*(i+j+1);
    for(int i=0;i<n;i++)b[i]=A[i]+2.0*A[n+i];      // b = col0 + 2*col1 exactly
    kc = Lingebra::leastSquareFit_QR(n,m,A,b,x,Q,R,piv,vw,1e-8);
    double res=0; for(int i=0;i<n;i++){ double s=b[i]; for(int j=0;j<m;j++)s-=A[j*n+i]*x[j]; res+=s*s; }
    printf("      rank=%i residual=%.3e x=(%g %g %g %g)\n", kc, sqrt(res), x[0],x[1],x[2],x[3]);
    check( sqrt(res)<1e-10,      "exact system: residual ~0" );

    printf("%s (%i failures)\n", nFail?"FAILED":"ALL PASS", nFail);
    return nFail;
}
