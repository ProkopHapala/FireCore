
#ifndef  ProjectiveDynamics_d_cpp
#define  ProjectiveDynamics_d_cpp

#include <stdio.h>
#include <string.h>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
//#include "datatypes.h"  
#include "VecN.h"

#include "Buckets.h"

#include "raytrace.h"
#include "geom3D.h"
#include "Interfaces.h"

#include "Cholesky.h"
#include "CG.h"
#include "arrayAlgs.h"

#include "SparseMatrix.h"
#include "SparseMatrix2.h"

#include "ProjectiveDynamics_d.h"

#include "IO_utils.h"
#include "testUtils.h"

double checkDist(int n, const Vec3d* vec, const Vec3d* ref, int verb, double tol ){
    double err=0;
    for(int i=0; i<n; i++){ 
        double ei = (vec[i]-ref[i]).norm();
        err = fmax( err, ei );
        if(verb>1){
           if(ei>tol) printf( "[%i] vec(%20.10f,%20.10f,%20.10f) ref(%20.10f,%20.10f,%20.10f)\n", i, vec[i].x,vec[i].y,vec[i].z,  ref[i].x,ref[i].y,ref[i].z ); 
        }
    }; 
    if( verb>0 )printf( "Err| v - ref |=%g\n", err );
    return err;
}

void print_vector( int n, double * a, int pitch, int j0, int j1 ){
    for(int i=0; i<n; i++){
        double* ai = a + i*pitch;	
        for (int j=j0; j<j1; j++ ){	printf( "%f ", ai[j] );	}
    }    
    printf( "\n" );
}


// ============ Functions of ProjectiveDynamics_d ============

    void ProjectiveDynamics_d::recalloc( int nPoint_, int nNeighMax_, int nBonds_){
        nPoint = nPoint_; nNeighMax = nNeighMax_;
        nNeighTot = nPoint*nNeighMax;
        _realloc0( points, nPoint   , Quat4dZero );
        _realloc0( forces, nPoint   , Quat4dZero );
        _realloc0( vel,    nPoint   , Quat4dZero );
        _realloc0( vel0,    nPoint   , Quat4dZero );
        _realloc0( bvec,   nPoint   , Quat4dZero );
        _realloc0( bvec0,  nPoint   , Vec3dZero  );
        _realloc0( ps_0,   nPoint   , Vec3dZero  );
        _realloc0( ps_cor       ,nPoint, Vec3dZero );
        _realloc0( ps_pred      ,nPoint, Vec3dZero );
        _realloc0( linsolve_yy , nPoint, Vec3dZero ); // used also to store momentum in MomentumDiff solver
        //_realloc( impuls, nPoint  , Quat4dZero ); );
        _realloc0( params,  nNeighTot , Quat4dZero );
        _realloc0( neighs,  nNeighTot , 0 );
        _realloc0( neighBs, nNeighTot , Vec2iZero );
        _realloc0( neighB2s,nNeighTot , 0 );

        if(nBonds_>0){
            nBonds=nBonds_;
            _realloc0( bonds,     nBonds, Vec2iZero );
            _realloc0( strain,    nBonds, 0.0 );
            //_realloc( l0s,       nBonds );
            _realloc0( maxStrain, nBonds, Vec2dZero );
            _realloc0( bparams,   nBonds, Quat4dZero );
        }
    }

    void ProjectiveDynamics_d::realloc_LinearSystem( bool bCholesky, int nNeighMaxLDLT_, bool bDens ){
        printf( "ProjectiveDynamics_d_d::realloc_LinearSystem() nPoint=%i nNeighMaxLDLT=%i nNeighMax=%i\n", nPoint, nNeighMaxLDLT_, nNeighMax );
        nNeighMaxLDLT = nNeighMaxLDLT_;  if(nNeighMaxLDLT<nNeighMax){ printf("ERROR in ProjectiveDynamics_d::prepare_Cholesky(): nNeighMaxLDLT(%i)<nNeighMax(%i) \n", nNeighMaxLDLT, nNeighMax); exit(0); }
        int n2 = nPoint*nPoint;
        // --- sparse Linear system Matrix and its Cholesky L*D*L^T decomposition
        if(bDens){ _realloc0( PDmat,  n2,     0.0 ); }
        _realloc0( linsolve_b ,nPoint, Vec3dZero );
        if(bCholesky){
            _realloc0( LDLT_L,      n2,     0.0             );
            _realloc0( LDLT_D,      nPoint, 0.0             );
            _realloc0( neighsLDLT, nPoint*nNeighMaxLDLT, -1 );
        }
    }

    double ProjectiveDynamics_d::norm_butFixed( Vec3d* ps ){
        double r2 = 0.0;
        for (int i=0; i<nPoint; i++){ 
            if(kFix){ if(kFix[i]>1e-8) continue; }
            r2 += ps[i].norm2();
        }
        return sqrt(r2);;
    }
    double ProjectiveDynamics_d::norm_butFixed( Quat4d* ps ){
        double r2 = 0.0;
        for (int i=0; i<nPoint; i++){ 
            if(kFix){ if(kFix[i]>1e-8) continue; }
            r2 += ps[i].f.norm2();
        }
        return sqrt(r2);;
    }


    // =================== Solver Using Projective Dynamics and Cholesky Decomposition

    __attribute__((hot)) 
    void ProjectiveDynamics_d::make_PD_Matrix( double* A, double dt ) {
        printf( "ProjectiveDynamics_d_d::make_PD_Matrix() dt=%g\n", dt );
        //for(int i=0; i<nPoint; i++){ printf("mass[%i] %g\n", points[i].w ); }; // debug
        //for(int i=0; i<nBonds; i++){ printf("ks[%i] %g\n",   params[i].y ); }; // debug
        int n = nPoint;
        int n2 = n*n;
        //for(int i=0; i<np2; i++){ PDmat[i]=0.0; }
        double idt2 = 1.0 / (dt * dt);
        for (int i = 0; i < n; ++i) { // point i
            double   Aii  = points[i].w * idt2; // Assuming points[i].w stores the mass
            if(kFix) Aii += kFix[i];            // fixed point
            //printf( "==i,i %i %g %g %g \n", i, Aii, points[i].w, idt2  );
            Vec2i* ngi = neighBs + (i*nNeighMax);
            for( int ing=0; ing<nNeighMax; ing++ ){
                Vec2i  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                //if(ib<0){ printf("ERROR int ProjectiveDynamics_d_d::make_PD_Matrix() invalid neighbor [i=%i,ing=%i] ng(%i,%i) neigh=%i\n", i, ing, ib, ng.x, neighs[i*nNeighMax+ing] ); exit(0); }  
                //printf( "make_PD_Matrix()[i=%i,ing=%i] ng(%i,%i)\n", i, ing, ib, ng.x );
                Vec2i   b = bonds[ib];
                //int j  = neighs[b.x];
                int j    = b.y;
                double k = getLinearBondStiffness( ib );  

                //printf( "i,ib %i %i %g \n", i, ib+1, k  );

                Aii += k;
                // if (j > i) {
                //     A[i+j*n] = -k;
                //     A[j+i*n] = -k;
                // }

                if     ( b.y > i ){
                    A[i  *n+b.y] = -k;
                    A[b.y*n+i  ] = -k;
                }else if ( b.x > i ){
                    A[i  *n+b.x] = -k;
                    A[b.x*n+i  ] = -k;
                }
            }
            A[i+i*n] += Aii;
        }
    }

    __attribute__((hot)) 
    void ProjectiveDynamics_d::make_PDmat_sparse( SparseMatrix<double>& A, double dt, bool bRealloc ) {
        printf( "ProjectiveDynamics_d_d::make_PDmat_sparse() dt=%g nNeighMax=%i \n", dt, nNeighMax );
        double idt2 = 1.0 / (dt * dt);
        // --- count neighbors
        if( bRealloc ){ A.realloc( nPoint, nNeighMax+1 ); }
        
        for (int i=0; i<nPoint; i++) {
            double   Aii  = points[i].w * idt2; 
            if(kFix) Aii += kFix[i];  // fixed point
            Vec2i* ngi = neighBs + (i*nNeighMax);
            for( int ing=0; ing<nNeighMax; ing++ ){
                Vec2i  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                Vec2i   b = bonds[ib];
                int j    = b.y;
                double k = getLinearBondStiffness( ib );  
                //printf( "i,ib %i %i %g \n", i, ib+1, k  );
                Aii += k;
                if       ( b.y > i ){
                    A.set( i,   b.y, -k );
                    A.set( b.y, i  , -k );
                }else if ( b.x > i ){
                    A.set( i,   b.x, -k );
                    A.set( b.x, i  , -k );
                }
            }
            A.set(i,i, Aii );
        }
    }

    __attribute__((hot)) 
    void ProjectiveDynamics_d::rhs_ProjectiveDynamics(const Vec3d* pnew, Vec3d* b) {
        double idt2 = 1.0 / (dt * dt);
        for (int i=0; i<nPoint; i++) {
            double   Ii  = points[i].w*idt2; 
            if(kFix) Ii += kFix[i];
            Vec3d bi; bi.set_mul(pnew[i], Ii );  // points[i].w is the mass
            Vec2i* ngi = neighBs + (i*nNeighMax);
            int ni = 0;
            for( int ing=0; ing<nNeighMax; ing++ ){
                Vec2i  ng = ngi[ing];
                int   ib = ng.y;
                if(ib<0) break;
                double k  = getLinearBondStiffness( ib );  
                double l0 = bparams[ib].x;
                Vec2i   b  = bonds[ib];
                int j     = (b.x == i) ? b.y : b.x;
                //printf( "rhs[i=%2i,ib=%2i,j=%i] k=%g l0=%g\n", i, ib, j, k, l0 );
                Vec3d d = pnew[i] -  pnew[j];
                bi.add_mul(d, k * l0 / d.norm());  // params[ib].x is l0
                ni++;
            }
            //printf( "rhs[i=%2i] ni=%i bi(%g,%g,%g)\n", i, ni, bi.x,bi.y,bi.z );
            b[i] = bi;
        }
    }

    __attribute__((hot)) 
    Vec3d ProjectiveDynamics_d::rhs_ProjectiveDynamics_i(int i, const Vec3d* pnew) {
        double idt2 = 1.0 / (dt * dt);
        double   Ii  = points[i].w*idt2; 
        if(kFix) Ii += kFix[i];
        Vec3d bi; bi.set_mul(pnew[i], Ii );  // points[i].w is the mass
        Vec2i* ngi = neighBs + (i*nNeighMax);
        int ni = 0;
        for( int ing=0; ing<nNeighMax; ing++ ){
            Vec2i  ng = ngi[ing];
            int   ib = ng.y;
            if(ib<0) break;
            double k  = getLinearBondStiffness( ib );  
            double l0 = bparams[ib].x;
            Vec2i   b  = bonds[ib];
            int j     = (b.x == i) ? b.y : b.x;
            //printf( "rhs[i=%2i,ib=%2i,j=%i] k=%g l0=%g\n", i, ib, j, k, l0 );
            Vec3d d = pnew[i] -  pnew[j];
            bi.add_mul(d, k * l0 / d.norm());  // params[ib].x is l0
            ni++;
        }
        //printf( "rhs[i=%2i] ni=%i bi(%g,%g,%g)\n", i, ni, bi.x,bi.y,bi.z );
        return bi;
    }

    __attribute__((hot)) 
    void ProjectiveDynamics_d::rhs_ProjectiveDynamics_(const Vec3d* pnew, Vec3d* b) {
        double idt2 = 1.0 / (dt * dt);
        for (int i=0; i<nPoint; i++) {
            b[i] = rhs_ProjectiveDynamics_i(i,pnew);
        }
    }

    void ProjectiveDynamics_d::updatePD_RHS( const Vec3d* pnew, Quat4d* bvec ){  // RHS for projective dynamics
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = pnew[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d bi; bi.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                 // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                const Quat4d par = params[j];
                const Vec3d  pj  = pnew[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                bi   .add_mul( dij, k*par.x/l );   // RHS      bi  = \sum_j K_{ij} d_{ij} (+ inertial term M_i/dt^2 p'_i )
                Aii  += k;                         // diagonal Aii =  \sum_j k_{ij}       (+ inertial term M_i/dt^2      )
            }
            bvec[i].f = bi;  // use forces array to store RHS vector
            bvec[i].w = Aii;
        }
    }

    // evaluate f = A p
    void ProjectiveDynamics_d::dotPD( const Vec3d* p, Vec3d* f ){  // RHS for projective dynamics
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = p[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d fi           = Vec3dZero;   // sum_j  =  \sum_j ( K_{ij} p_j} )   + Aii * pi        
            double Aii  = Ii;                 // A_{ii} =  \sum_j   K_{ij}          + M_i/dt^2 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                const Quat4d par = params[j];
                const Vec3d  pj  = p[jG];
                const double k   = par.z;
                fi.add_mul( pj, -k ); 
                Aii  += k;                  
            }
            fi .add_mul( pi, Aii );
            f[i] = fi;
        }
    }

    void ProjectiveDynamics_d::updatePD_dRHS( const Vec3d* pnew, Quat4d* bvec ){  // RHS for projective dynamics
        // Calculates db = b - A p0 ( i.e. combines dotPD and updatePD_RHS into a single loop) 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = pnew[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d bi = Vec3dZero; 
            //bi.set_mul(pi,  Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )   
            //bi.add_mul(pi, -Ii ); 
            double Aii = Ii;                 // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                const Quat4d par = params[j];
                const Vec3d  pj  = pnew[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                //bi   .add_mul( pi, -k ); 
                //bi   .add_mul( pj,  k ); 
                //bi   .add_mul( dij,-k );
                //bi   .add_mul( dij, k*par.x/l );   // RHS      bi  = \sum_j K_{ij} d_{ij} (+ inertial term M_i/dt^2 p'_i )
                bi   .add_mul( dij, k*(par.x/l - 1) ); 
                Aii  += k;                         // diagonal Aii =  \sum_j k_{ij}       (+ inertial term M_i/dt^2      )
            }
            //bi .add_mul( pi, -Aii );
            bvec[i].f = bi;  // use forces array to store RHS vector
            bvec[i].w = Aii;
        }
    }

    void ProjectiveDynamics_d::updateJacobi_lin( Vec3d* ps_in, Vec3d* ps_out, Quat4d* bvec, Vec3d* rs ){
        for(int i=0; i<nPoint; i++){
            const int j0 = i * nNeighMax;
            Vec3d sum_j  = Vec3dZero;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                sum_j.add_mul( ps_in[jG],  params[j].z );   // kPull
            }
            const Quat4d bi     = bvec[i];  // forces array stores RHS vector
            if(ps_out ) ps_out[i] = (bi.f + sum_j       )*(1/bi.w);
            if(rs     ) rs    [i] = (bi.f + sum_j - ps_in[i]*bi.w );
        }
    }

    void ProjectiveDynamics_d::updateGaussSeidel_lin( Vec3d* ps, Quat4d* bvec ){
        for(int i=0; i<nPoint; i++){
            const int j0 = i * nNeighMax;
            Vec3d sum_j  = Vec3dZero;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                sum_j.add_mul( ps[jG],  params[j].z );   // kPull
            }
            const Quat4d bi     = bvec[i];  // forces array stores RHS vector
            ps[i]  = (bi.f + sum_j)*(1/bi.w);
        }
    }

    void ProjectiveDynamics_d::evalTrussForce( Vec3d* ps_in, Vec3d* force ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps_in[i];
            Vec3d          fi = Vec3dZero;
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                const Quat4d par = params[j];
                const Vec3d  pj  = ps_in[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double dl  = l - par.x; 
                const double k   = par.z;
                fi.add_mul( dij, k*dl/l ); 
            }
            force[i]  = fi;
        }
    }

    void ProjectiveDynamics_d::updateJacobi_fly( Vec3d* ps_in, Vec3d* ps_out ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps_in[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d sum_j; sum_j.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                       // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                //   update RHS bi
                const Quat4d par = params[j];
                const Vec3d  pj  = ps_in[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                sum_j.add_mul( dij, k*par.x/l );   //  b_i  +=  \sum_j ( K_{ij} d_{ij} )   
                sum_j.add_mul( pj , k         );   //  s_j  +=  \sum_j ( K_{ij} p_j    )
                Aii  += k;

            }
            sum_j.mul( 1/Aii );
            ps_out[i]  = sum_j;
        }
    }

    void ProjectiveDynamics_d::updateGaussSeidel_fly( Vec3d* ps ){
        // This solver calculates bi and Aii on-the fly, preventing updatePD_RHS call and possibly improving numerical stability (?)
        // NOTE: this is goes beyond linar solution of Ap=b because we update b every iteration, and use non-linear operations like d_{ij}=(pi-pj)/|pi-pj| 
        for(int i=0; i<nPoint; i++){
            const Vec3d    pi = ps[i];
            const double idt2 = 1.0 / (dt * dt);
            double         Ii = points[i].w*idt2; 
            if(kFix)       Ii += kFix[i];
            Vec3d sum_j; sum_j.set_mul(pi, Ii );   // b_i    =  M_i/dt^2 p'_i   +  \sum_j ( K_{ij} d_{ij} )         
            double Aii = Ii;                       // A_{ii} =  M_i/dt^2        +  \sum_j   K_{ij} 
            const int j0 = i * nNeighMax;
            for(int jj=0; jj<nNeighMax; jj++){
                const int j  = j0 + jj;
                const int jG = neighs[j];
                if(jG == -1) break;
                
                //   update RHS bi
                const Quat4d par = params[j];
                const Vec3d  pj  = ps[jG];
                const Vec3d  dij = pi - pj;
                const double l   = dij.norm();
                const double k   = par.z;
                sum_j.add_mul( dij, k*par.x/l );   //  b_i  +=  \sum_j ( K_{ij} d_{ij} )   
                sum_j.add_mul( pj , k         );   //  s_j  +=  \sum_j ( K_{ij} p_j    )
                Aii  += k;

            }
            sum_j.mul( 1/Aii );
            ps[i]  = sum_j;
        }
    }

    void ProjectiveDynamics_d::updateIterativeMomentum( Vec3d* psa, Vec3d* psb ){
        //Vec3d* psa = ps_in;
        //Vec3d* psb = ps_out;
        //if(LinSolveMethod::JacobiMomentum == (LinSolveMethod)linSolveMethod  ){
        updatePD_RHS(psa, bvec ); // TODO: we need to calculate this only whem we use linear solvers, not if we use _fly solvers
        //}
        for (int i = 0; i < nSolverIters; i++) {  
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::JacobiMomentum:    { updateJacobi_lin( psa, psb, bvec ); } break;
                case LinSolveMethod::JacobiFlyMomentum: { updateJacobi_fly( psa, psb );       } break;
                case LinSolveMethod::GSMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_lin( psb, bvec ); 
                } break; 
                case LinSolveMethod::GSFlyMomentum:     { 
                    for (int j=0; j<nPoint; j++){ psb[j]=psa[j]; }
                    updateGaussSeidel_fly( psb ); 
                } break; 
            }
            //updateJacobi_lin( psa, psb, bvec );      //    p'_{k+1}    = Solve(A,p_k,b)
            double bmix = mixer.get_bmix( i );       //    bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; // make sure we stop momentum mixing for the last iteration      
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( linsolve_yy[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                linsolve_yy[j] = p - psa[j];         //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }   // we use linsolve_yy to store momentum
            //for (int i=0; i<nPoint; i++){ ps_pred[i]=ps_cor[i]; }
        }

        if( bApplyResudualForce ){
            double dt2 = dt * dt;
            evalTrussForce( psa, linsolve_b );
            //updateJacobi_lin( psa, 0, bvec, linsolve_b );
            // for (int i=0; i<nPoint; i++){ 
            //     if(kFix){ if(kFix[i]>1e-8) continue; }
            //     psa[i].add_mul( linsolve_b[i], residualForceFactor*dt2/points[i].w ); 
            // }
        }

        for (int i=0; i<nPoint; i++){ psb[i]=psa[i]; }
    }

    void ProjectiveDynamics_d::updateIterativeJacobi( Vec3d* psa, Vec3d* psb ){
        updatePD_RHS(psa, bvec );
        for (int i = 0; i < nSolverIters; i++) {  
            updateJacobi_lin( psa, psb, bvec ); 
            for (int i=0; i<nPoint; i++){ psa[i]=psb[i]; }
        }
    }

    // Here we replace Ap=b   by    A(p-p0) = b-Ap0
    void ProjectiveDynamics_d::updateIterativeJacobiDiff( Vec3d* psa, Vec3d* psb ){
        //printf( "updateIterativeJacobiDiff() \n" );
        const bool bUse_dRHS = true; 
        //const bool bUse_dRHS = false; 
        if(bUse_dRHS){
            updatePD_dRHS( psa, bvec  );    // b = Kd - Ap' = \sum_j{ k_{ij} ( l0_{ij}/l_{ij} - 1 ) ( p_i - p_j )
        }else{
            dotPD       ( psa, bvec0 );     // b0 = Ap0
            updatePD_RHS( psa, bvec  );     // b = Kd + Ip'
            for (int i=0; i<nPoint; i++){ bvec[i].f.sub(bvec0[i]);  }       //  db = b - b0 
        }
    //    updatePD_RHS( psa, bvec  );    // b = Kd + Ip'
    //    dotPD       ( psa, bvec0 );     // b0 = Ap0
    //    for (int i=0; i<nPoint; i++){ bvec0[i]= bvec[i].f - bvec0[i];  } 
    //    updatePD_dRHS( psa, bvec  ); 
    //    for (int i=0; i<nPoint; i++){ printf( "bvec %i |dr|: %12.6e b(%12.6e,%12.6e,%12.6e) b0(%12.6e,%12.6e,%12.6e)\n", i, (bvec[i].f-bvec0[i]).norm(), bvec[i].x, bvec[i].y, bvec[i].z, bvec0[i].x, bvec0[i].y, bvec0[i].z ); } 
    //    exit(0);
        for (int i=0; i<nPoint; i++){ ps_0[i]=psa[i]; psa[i]=Vec3dZero; }   //  dp = p - p0      
        for (int i = 0; i < nSolverIters; i++) {  
            updateJacobi_lin( psa, psb, bvec );
            for (int i=0; i<nPoint; i++){ psa[i]=psb[i]; }
        }
        for (int i=0; i<nPoint; i++){ psb[i].add(ps_0[i]); }   //  p = p0 + dp 
        // Debug:
        //    double db_norm  = norm_butFixed( bvec );
        //    double dp_norm  = norm_butFixed( psa );
        //    printf( "updateIterativeJacobiDiff() p0_norm: %g b0_norm: %g b_norm: %g db_norm: %g dp_norm: %g \n", p0_norm, b0_norm, b_norm, db_norm, dp_norm );
    }

    void ProjectiveDynamics_d::updateIterativeMomentumDiff( Vec3d* psa, Vec3d* psb ){
        const bool bUse_dRHS = true; 
        updatePD_dRHS( psa, bvec  );    // b = Kd - Ap' = \sum_j{ k_{ij} ( l0_{ij}/l_{ij} - 1 ) ( p_i - p_j )
        for (int i=0; i<nPoint; i++){ ps_0[i]=psa[i]; psa[i]=Vec3dZero; }   //  dp = p - p0 , start with zero      
        for (int i = 0; i < nSolverIters; i++) {  
            updateJacobi_lin( psa, psb, bvec );
            double bmix = mixer.get_bmix( i );       //    bmix_k = bmix(k)
            if( (i==0) || (i>=(nSolverIters-1) ) ) bmix = 0.0; // make sure we stop momentum mixing for the last iteration      
            for(int j=0; j<nPoint; j++){ 
                Vec3d p = psb[j]; 
                p.add_mul( linsolve_yy[j], bmix );   //    p_{k+1} = p'_k + bmix d_k
                linsolve_yy[j] = p - psa[j];         //    d_{k+1} = p_{k+1} - p_k
                psa[j]         = p;
            }
        }
        for (int i=0; i<nPoint; i++){ psb[i].add(ps_0[i]); }   //  p = p0 + dp 
    }

    // Here we replace Ap=b   by    A(p-p0) = b-Ap0
    void ProjectiveDynamics_d::updateIterativeExternDiff( Vec3d* psa, Vec3d* psb ){
        dotPD       ( psa, bvec0 );     // b0 = Ap0
        updatePD_RHS( psa, bvec  );     // b = Kd + Ip'
        for (int i=0; i<nPoint; i++){ Quat4d b=bvec[i]; b.f.sub(bvec0[i]); extern_b[i]=(Quat4f)b; }   //  db = b - b0 
        for (int i=0; i<nPoint; i++){ extern_x[i] = Quat4fZero;           }       //  dp = p - p0   
        extern_solve();
        for (int i=0; i<nPoint; i++){ psb[i] = (Vec3d)extern_x[i].f + psa[i]; }   //  p = p0 + dp 
    }

    /*    double run_updateJacobi_smart( int psa, int psb, int itr ) {
        //run_updateJacobi_neighs(psa, psb);
        // I want some function of form     A/( B + itr ) which interpolates between  bmix_end and bmix_start when itr transitions from nitr_end to nitr_start
        if(itr<nitr_start){ run_updateJacobi_neighs(psa, psb);  }
        else              { 
            if ( itr > nitr_end ) { bmix.y = bmix_end; }
            else                  { bmix.y = bmix_start + (bmix_end - bmix_start ) * (itr - nitr_start) / (nitr_end - nitr_start); }
            run_updateJacobi_mix(psa, psb); 
        }
        
       return bmix.y;
    }
    */

    void ProjectiveDynamics_d::prepare_LinearSystem( bool bRealloc, bool bCholesky, int nNeighMaxLDLT_, bool bDens ){ 
        double dt = this->dt;
        printf( "ProjectiveDynamics_d_d::prepare_LinearSystem() nPoint=%i dt=%g nNeighMaxLDLT_=%i\n", nPoint, dt, nNeighMaxLDLT_ );
        //nNeighMaxLDLT=nNeighMaxLDLT_;
        if(bRealloc)realloc_LinearSystem( bCholesky, nNeighMaxLDLT_, bDens );

        this->dt = dt;
        int n2 = nPoint*nPoint;

        //long t0;
        //t0=getCPUticks(); make_PD_Matrix(    PDmat,    dt );         printf("@time make_PD_Matrix()    t= %g [MTicks] \n", (getCPUticks()-t0)*1.0e-6 );
        //t0=getCPUticks(); make_PDmat_sparse( PDsparse, dt, true );   printf("@time make_PD_Matrix()    t= %g [MTicks] \n", (getCPUticks()-t0)*1.0e-6 );

        //timeit( "TIME make_PDmat_sparse() t= %g [MTicks]\n", 1e-6, [&](){ make_PDmat_sparse( PDsparse, dt, true ); });
        mat2file<int>   ( "PDsparse_inds.log",  nPoint, nNeighMax+1, (int*)   PDsparse.inds, "%5i " );
        mat2file<double>( "PDsparse_vals.log",  nPoint, nNeighMax+1, (double*)PDsparse.vals );

        if(bDens){
            //timeit( "TIME make_PD_Matrix()    t= %g [MTicks]\n", 1e-6, [&](){ make_PD_Matrix(    PDmat,    dt );       });
            if( PDsparse.checkDens( PDmat, 1e-16, true ) ){ printf("ERROR in ProjectiveDynamics_d_d::prepare_LinearSystem() PDsparse does not match PDmat => exit \n"); exit(0); }
        }

        // { // Check M.dot(vec) sparse
        //     for(int i=0; i<nPoint; i++){ ps_pred[i] = points[i].f *1.1; }; 
        //     dotM_ax                   ( nPoint,3, PDmat, (double*)ps_pred, (double*)linsolve_b );
        //     PDsparse.dot_dens_vector_m(   3,             (double*)ps_pred, (double*)ps_cor     );
        //     for(int i=0; i<nPoint; i++){ printf( "[%i] PDmat*v(%20.10f,%20.10f,%20.10f) PDsparse(%20.10f,%20.10f,%20.10f)\n", i, linsolve_b[i].x,linsolve_b[i].y,linsolve_b[i].z,  ps_cor[i].x,ps_cor[i].y,ps_cor[i].z ); }; 
        //     exit(0);
        // }
        
        if( bCholesky ){
            for(int i=0; i<nPoint; i++){  for(int j=0; j<nNeighMaxLDLT; j++){ if(j>=nNeighMax){neighsLDLT[i*nNeighMaxLDLT+j]=-1;}else{ neighsLDLT[i*nNeighMaxLDLT+j]=neighs[i*nNeighMax+j]; };  } }
            
            bool bSparse = false;
            //bool bSparse = true;
            if(bSparse){
                Lingebra::CholeskyDecomp_LDLT_sparse( PDmat, LDLT_L, LDLT_D, neighsLDLT, nPoint, nNeighMaxLDLT );
                sortNeighs( nPoint, nNeighMaxLDLT, neighsLDLT);
            }else{
                //Lingebra::CholeskyDecomp_LDLT( PDmat, LDLT_L, LDLT_D, nPoint, true );
                //LsparseT.fromDense( nPoint, LDLT_L, 1.e-16, true );
                //for(int i=0; i<n2; i++){ LDLT_L[i]=0; }
                Lingebra::CholeskyDecomp_LDLT( PDmat, LDLT_L, LDLT_D, nPoint );
                //Lsparse.fromDense( nPoint, LDLT_L, 1.e-16 );
                Lsparse.fromDense( nPoint, LDLT_L, 1e-300 );
                LsparseT.fromFwdSubT( nPoint, LDLT_L);
                //mat2file<double>( "PDmat.log",  nPoint,nPoint, PDmat  );
                //mat2file<double>( "LDLT_L.log", nPoint,nPoint, LDLT_L );
            }
           

            // {   printf( "# ---------- Test Cholesky Lsparse \n");
            //     for(int i=0; i<nPoint; i++){ ps_pred[i] = points[i].f*1.1; }; 
            //     rhs_ProjectiveDynamics( ps_pred, linsolve_b );
            //     Lingebra::forward_substitution_m( LDLT_L, (double*)linsolve_b,  (double*)linsolve_yy, nPoint, 3 );
            //     Lsparse.fwd_subs_m( 3,  (double*)linsolve_b,  (double*)ps_cor );
            //     checkDist( nPoint, ps_cor, linsolve_yy, 2 );
            // }

            //mat2file<double>( "LDLT_L.log", nPoint,nPoint, LDLT_L );
            Lsparse.fprint_inds("Lsparse_inds.log");
            //Lsparse.fprint_vals("Lsparse_vals.log");

            LsparseT.fprint_inds("LsparseT_inds.log");
            //LsparseT.fprint_vals("LsparseT_vals.log");
            //exit(0);

            double* L_reconstructed  = Lsparse .reconstruct_dense( );
            double* LT_reconstructed = LsparseT.reconstruct_dense( );
            //writeMatrix("LDLT_L_dense.log",         nPoint, nPoint,  LDLT_L, false);
            //writeMatrix("LDLT_LT_reconstructed.log", nPoint, nPoint, LT_reconstructed, false);
            //writeMatrix("LDLT_L_reconstructed.log", nPoint, nPoint,  L_reconstructed, false);
            mat2file<double>( "LDLT_L_dense.log", nPoint,nPoint, LDLT_L );
            mat2file<double>( "LDLT_LT_reconstructed.log", nPoint,nPoint, LT_reconstructed );
            mat2file<double>( "LDLT_L_reconstructed.log", nPoint,nPoint,  L_reconstructed );
            delete[] L_reconstructed;
            delete[] LT_reconstructed;
        }


    }

    // void ProjectiveDynamics_d::dampPoints( double Vdamping ){
    //     for (int i : damped_points){
    //         Vec3d vi  = vel [i].f;
    //         Vec3d vi0 = vel0[i].f;
    //         Vec3d acc = (vi0-vi) * Vdamping;
    //         vel0[i].f.add_mul( acc, -dt );
    //         forces[i].f.add_mul( acc, points[i].w );
    //     }
    // }

    __attribute__((hot)) 
    void ProjectiveDynamics_d::run_LinSolve(int niter) {
        //printf( "ProjectiveDynamics_d::run_LinSolve()  linSolveMethod=%i nSolverIters=%i ps_cor=%p points=%p \n", linSolveMethod, nSolverIters, ps_cor, points  );
        const int m=3;
        memcpy(ps_cor, points, nPoint * sizeof(Vec3d));
        double dt2 = dt * dt;
        double inv_dt = 1/dt;
        cleanForce();

        for (int iter = 0; iter < niter; iter++) {
            // Evaluate forces (assuming you have a method for this)
            
            // ---- Update Forces
            cleanForce();
            //evalForces();
            if (user_update){ user_update(dt);}
            //evalEdgeVerts();
            //dampPoints( 2.5 );   // Damping - this is not rigorous - it does not perfectly preserve momentum / angular momentum, resp. points have like double the mass

            //  Predictor step
            for (int i=0;i<nPoint;i++){ 
                forces[i].f.add( getPointForce( i ) );
                //forces[i].f.add_mul( Gravity, points[i].w ); // Grafivity
                //forces[i].f.add_mul( vel[i].f, Cdrag ); // Drag force
                ps_pred[i] = points[i].f + vel[i].f*dt + forces[i].f*(dt2/points[i].w); 
                //printf( "ps_pred[%3i](%10.6f,%10.6f,%10.6f) v(%10.6f,%10.6f,%10.6f) p(%10.6f,%10.6f,%10.6f) dt=%g \n", i, ps_pred[i].x,ps_pred[i].y,ps_pred[i].z, vel[i].x,vel[i].y,vel[i].z, points[i].x,points[i].y,points[i].z, dt );
            }

            // Apply fixed constraints
            if(kFix) for (int i = 0; i < nPoint; i++) { if (kFix[i] > 1e-8 ) { ps_pred[i] = points[i].f; } }
            // Compute right-hand side
            
            long t0 = getCPUticks();
            switch( (LinSolveMethod)linSolveMethod ){
                case LinSolveMethod::Jacobi:{
                    //printf( "Jacobi nSolverIters=%i \n", nSolverIters  );
                    updateIterativeJacobi( ps_pred, ps_cor  );
                } break;
                case LinSolveMethod::JacobiDiff:{
                    updateIterativeJacobiDiff( ps_pred, ps_cor );
                } break;
                case LinSolveMethod::MomentumDiff:{
                    updateIterativeMomentumDiff( ps_pred, ps_cor );
                }break;
                case LinSolveMethod::ExternDiff:{
                    updateIterativeExternDiff( ps_pred, ps_cor );
                } break;
                case LinSolveMethod::GaussSeidel:{
                    //printf( "Jacobi nSolverIters=%i \n", nSolverIters  );
                    for (int i=0; i<nPoint; i++){ ps_cor[i]=ps_pred[i]; }
                    updatePD_RHS(ps_cor, bvec );
                    for (int i = 0; i < nSolverIters; i++) {  
                        updateGaussSeidel_lin( ps_cor, bvec ); 
                    }
                } break;
                case LinSolveMethod::JacobiMomentum:
                case LinSolveMethod::GSMomentum: 
                case LinSolveMethod::JacobiFlyMomentum:
                case LinSolveMethod::GSFlyMomentum:
                    updateIterativeMomentum( ps_pred, ps_cor );
                    break;
                case LinSolveMethod::CholeskySparse:{
                    //printf("ProjectiveDynamics_d::run_LinSolve()  LinSolveMethod::CholeskySparse \n");
                    // Solve using LDLT decomposition (assuming you have this method)
                    //solve_LDLT_sparse(b, ps_cor);
                    rhs_ProjectiveDynamics(ps_pred, linsolve_b );
                    Lingebra::forward_substitution_sparse           ( nPoint,m,  LDLT_L, (double*)linsolve_b,  (double*)linsolve_yy, neighsLDLT, nNeighMaxLDLT );
                    for (int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
                    Lingebra::forward_substitution_transposed_sparse( nPoint,m,  LDLT_L, (double*)linsolve_yy, (double*)ps_cor,      neighsLDLT, nNeighMaxLDLT );
                } break;
                case LinSolveMethod::Cholesky:{
                    //printf("ProjectiveDynamics_d::run_LinSolve()  LinSolveMethod::Cholesky \n");
                    //Lingebra::forward_substitution_m( LDLT_L, (double*)linsolve_b,  (double*)linsolve_yy, nPoint,m );
                    //Lsparse.fwd_subs_m( m,  (double*)linsolve_b,  (double*)ps_cor );
                    //if( checkDist( nPoint, ps_cor, linsolve_yy, 1 ) ){ printf("ERROR run_LinSolve.checkDist() => exit()"); exit(0); };
                    rhs_ProjectiveDynamics(ps_pred, linsolve_b );
                    Lsparse.fwd_subs_m( m,  (double*)linsolve_b,  (double*)linsolve_yy );
                    for (int i=0; i<nPoint; i++){ linsolve_yy[i].mul(1/LDLT_D[i]); } // Diagonal 
                    LsparseT.fwd_subs_T_m( m,  (double*)linsolve_yy,  (double*)ps_cor );
                    //Lingebra::forward_substitution_T_m( LDLT_L, (double*)linsolve_yy, (double*)ps_cor,      nPoint,m );
                    //if( checkDist( nPoint, ps_pred, ps_cor, 2 ) ){ printf("ERROR run_LinSolve.checkDist() => exit()"); exit(0); };

                } break;
            }
            //time_LinSolver += (getCPUticks()-t0)*1e-6;
            //mat2file<double>( "points.log",  nPoint,4, (double*)points      );
            //mat2file<double>( "vel.log",     nPoint,4, (double*)vel         );
            //mat2file<double>( "ps_pred.log", nPoint,3, (double*)ps_pred     );
            //mat2file<double>( "b.log",       nPoint,3, (double*)linsolve_b  );
            //mat2file<double>( "yy.log",      nPoint,3, (double*)linsolve_yy );
            //mat2file<double>( "ps_cor.log",  nPoint,3, (double*)ps_cor      );

            //exit(0);

            // Compute residual
            // double res = 0.0;
            // for (int i=0;i<nPoint;i++) {
            //     Vec3d  d = ps_cor[i] - points[i].f;
            //     double l = d.norm();
            //     if (l > res) res = l;
            // }
            // printf("residual[%d] : %f\n", iter, res);
            // Update velocity and points
            
            //   Corrector step
            for (int i=0;i<nPoint;i++) {
                //Vec3d dv = ps_cor[i] - ps_pred[i];
                //dv.mul(1.0 / dt);
                //vel   [i].f.add( dv );

                // To-Do : We need more rigorous approach how to update velocity

                // position-based velocity update
                double vr2 = vel[i].norm2();
                Vec3d v    = ps_cor[i] - points[i].f;
                v.mul(inv_dt);
                //v.mul( sqrt(  vel[i].norm2()/( v.norm2() + 1e-8 ) ) ); // This is unphysicsl
                if(kFix){ if( kFix[i] > 1e-8){ v=Vec3dZero; } }
                vel[i].f = v;
                // update positions
                points[i].f = ps_cor[i];
            }            
            time += dt;
            //printf( "STEP: %6i time: %16.8f p.y: %16.8f v.y: %16.8f f.y: %16.8f mass: %16.8f\n", iter, time, points[0].y, vel[0].y,  forces[0].y, points[0].w );
        } // for iter ( time stepping )
        //exit(0);
    }

    // =================== Truss Simulation

    void ProjectiveDynamics_d::updateInveriants(bool bPrint){
        mass=0;
        cog =Vec3dZero;
        vcog=Vec3dZero;
        for(int i=0; i<nPoint; i++){ double mi=points[i].w; mass+=mi;  cog.add_mul( points[i].f, mi ); vcog.add_mul( vel[i].f, mi ); }
        cog .mul( 1.0/mass );
        vcog.mul( 1.0/mass );
        I=Mat3dZero;
        L=Vec3dZero;
        torq=Vec3dZero;
        for(int i=0; i<nPoint; i++){ 
            double mi=points[i].w; 
            Vec3d d; 
            d   .set_sub   ( points[i].f, cog    );
            L   .add_crossw( d, vel[i]   .f, mi  );
            torq.add_cross ( d, forces[i].f      );
            I   .add_outer ( d, d, mi            );
        }
        if(bPrint){
            printf( "ProjectiveDynamics_d::updateInveriants() mass %g cog(%g,%g,%g) vcog(%g,%g,%g)  L(%g,%g,%g) torq(%g,%g,%g) \n", mass, cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, L.x,L.y,L.z, torq.x,torq.y,torq.z );
            printf( "ProjectiveDynamics_d::updateInveriants() I \n" ); I.print();
        }
    }

    double ProjectiveDynamics_d::evalBondTension(){
        double E = 0.0;
        for(int i=0;i<nBonds; i++ ){
            Vec2i   b  = bonds[i];
            double l0 = bparams[i].x;
            //double l0 = l0s[i];
            double l   = (points[b.y].f-points[b.x].f).norm();
            double dl  = l - l0;
            double s   = dl/l0;
            //if( fabs(s)>0.5 ){ printf( "evalBondTension[%i] strain=%g l=%g l0=%g\n", i, s, l, l0 ); }
            //if(i==6272){ printf( "evalBondTension[%i](%i,%i) strain=%g l=%g l0=%g\n", i, b.x,b.y, s, l, l0 ); }
            strain[i] = s;
            E+= 0.5*dl*dl*bparams[i].z;
            // ToDo: break the bond if strain > maxStrain;
        }
        //exit(0);
        return E;
    }

    void ProjectiveDynamics_d::printNeighs(int i){
        int j0 = i*nNeighMax;
        for(int jj=0;jj<nNeighMax;jj++){
            int j=j0+jj;
            int ing = neighs[j];
            if(ing<0) break;
            Quat4d par = params[j];
            printf( "ng[%i,%i|%i] l0,kP,kT,damp(%g,%g,%g,%g)\n", i, ing, jj, par.x,par.y,par.z,par.w );
        }
    }
    void ProjectiveDynamics_d::printAllNeighs(){ printf("ProjectiveDynamics_d::printAllNeighs(nPoint=%i,nNeighMax=%i)\n",nPoint,nNeighMax); for(int i=0;i<nPoint;i++){ printNeighs(i); }; };

    double ProjectiveDynamics_d::getFmax(){ 
        double fmax=0;
        for(int i=0; i<nPoint; i++){ double f=forces[i].norm(); fmax=fmax>f?fmax:f; }   
        //printf( "|fmax|=%g\n", fmax );
        return fmax;
    }

    void ProjectiveDynamics_d::setFixPoints( int n, int* fixPoints, double Kfix, bool bRealloc ){
        if(bRealloc) reallocFixed();
        for(int i=0;i<n;i++){
            int ip = fixPoints[i]; 
            if(ip>nPoint){ printf("ERROR in ProjectiveDynamics_d::setFixPoints fixing point %i > sim.nPoint(%i) \n", ip, nPoint); exit(0); }
            printf("ProjectiveDynamics_d_d::setFixPoints()[%3i]: %6i \n", i, ip ); 
            kFix[ip] = Kfix ; 
        }
    }


// inline void ProjectiveDynamics_d::setOpt( double dt_, double damp_ ){
//     dt      = dt_max   = dt_;  dt_min=0.1*dt_max;
//     damping = damp_max = damp_;
//     cleanForce( );
//     cleanVel  ( );
// }

//void ProjectiveDynamics_d::print_points(){}

#endif
