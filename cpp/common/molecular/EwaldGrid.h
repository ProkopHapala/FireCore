#ifndef EwaldGrid_h
#define EwaldGrid_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Bspline.h"

#ifdef WITH_FFTW
#include <fftw3.h>

__attribute__((hot))
void array2fftc( int n, const double* in, fftw_complex* out){
    for (int i=0; i<n; i++ ) {
        out[i][0] = in[i];
        out[i][1] = 0.0;
    }
}

__attribute__((hot))
void fftc2array( int n,  const fftw_complex* in, double* out) {
    for (int i=0; i<n; i++ ) { out[i] = in[i][0];}
}


/*

void array2fftc( int ntot, const double* in, fftw_complex* out) const {
    for (int ix = 0; ix < n.x; ix++) {
        for (int iy = 0; iy < n.y; iy++) {
            for (int iz = 0; iz < n.z; ++iz ) {
                int index = ix * n.y * n.z + iy * n.z + iz;
                out[index][0] = in[index];
                out[index][1] = 0.0;
            }
        }
    }
}

void fftc2array( const fftw_complex* in, double* out) const {
    for (int ix = 0; ix < n.x; ix++) {
        for (int iy = 0; iy < n.y; iy++) {
            for (int iz = 0; iz < n.z; ++iz ) {
                int index = ix * n.y * n.z + iy * n.z + iz;
                out[index] = in[index][0];
            }
        }
    }
}
*/

#endif


class EwaldGrid : public GridShape { public: 

__attribute__((hot)) 
void project_atom_on_grid( const Vec3d pi, const double qi, double* dens ) const {
    //printf("project_atom_on_grid() pi(%g,%g,%g) q=%g \n", pi.x, pi.y, pi.z, qi );
    const Vec3d gp = diCell.dot( pi-pos0 );
    const int ix = (int) gp.x;
    const int iy = (int) gp.y;
    const int iz = (int) gp.z;
    const double tx = gp.x - ix;
    const double ty = gp.y - iy;
    const double tz = gp.z - iz;
    // ToDo: Periodic Boundary Conditions for points which are close to the boundary
    const Quat4d bx = Bspline::basis(tx);
    const Quat4d by = Bspline::basis(ty);
    const Quat4d bz = Bspline::basis(tz);

    //printf("project_atom_on_grid() pi(%g,%g,%g) q=%g \n", pi.x, pi.y, pi.z, qi );
    const int nxy = n.x * n.y;
    //int ii=0;
    for (int dz = 0; dz < 4; dz++) {
        const int gz  = iz + dz - 1;
        const int iiz = gz*nxy;
        for (int dy = 0; dy < 4; dy++) {
            const int gy  = iy + dy - 1;
            const int iiy = iiz + gy*n.x;
            const double qbyz = qi * by.array[dy] * bz.array[dz];
            for (int dx = 0; dx < 4; dx++) {
                const int gx = ix + dx - 1;
                const int ig = gx + iiy;

                //printf("project_atom_on_grid()[%i] dxyz(%i,%i,%i) igxyz(%i,%i,%i) ig=%i /%i \n", ii, dx,dy,dz,   gx,gy,gz, ig, n.totprod() );
                // ToDo: Periodic Boundary Conditions for points which are close to the boundary
                //if (gx >= 0 && gx < nx && gy >= 0 && gy < ny && gz >= 0 && gz < nz) {
                    dens[ig] += qbyz * bx.array[dx]; 
                //}

                //ii++;
            }
        }
    }
}

__attribute__((hot)) 
void project_atoms_on_grid( int na, const Vec3d* apos, const double* qs, double* dens ) const {
    for (int ia=0; ia<na; ia++){
        project_atom_on_grid( apos[ia], qs[ia], dens );
    }
}

#ifdef WITH_FFTW

fftw_plan    fft_plan;
fftw_plan    ifft_plan;
fftw_complex *Vw=0,*V=0;

__attribute__((hot))
void laplace_reciprocal_kernel( fftw_complex* VV ){
    int nx=n.x;
    int ny=n.y;
    int nz=n.z;
    double freq_x = (2.0 * M_PI) / cell.a.norm();
    double freq_y = (2.0 * M_PI) / cell.b.norm();
    double freq_z = (2.0 * M_PI) / cell.c.norm();

    //printf("laplace_reciprocal_kernel() nxyz(%i,%i,%i) freq(%g,%g,%g) \n", nx, ny, nz, freq_x, freq_y, freq_z );

    for (int ix = 0; ix < nx; ix++) {
        const double kx  =  ( (ix <= ny / 2) ? ix : ix - nx ) * freq_x;
        for (int iy = 0; iy < ny; iy++) {
            const double ky   = ( (iy <= ny / 2) ? iy : iy - ny ) * freq_y;
            const double k2xy = ky*ky + kx*kx;
            const int    ixy  = (ix*ny + iy)*nz;
            for (int iz = 0; iz < nz; ++iz ) {
                const int i     = iz + ixy;
                const double kz = ( ( iz <= nz / 2) ? iz : iz - nz ) * freq_z;
                const double k2 = k2xy + kz*kz;
                if ( k2 > 1e-32 ){
                    const double invk2 = 1.0/k2;
                    VV[i][0] *= invk2; // Real part
                    VV[i][1] *= invk2; // Imaginary part
                } else {
                    VV[i][0] = 0.0; // Avoid division by zero (DC component)
                    VV[i][1] = 0.0;
                }
            }
        }
    }
}

    void prepare_laplace( int flags=-1 ){
        int ntot = n.totprod();
        V  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot );
        Vw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot );
        //fft_plan =            fftw_plan_dft_3d( nx,ny,nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        if(flags<0){ flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT; }

        // FFTW_PRESERVE_INPUT  // 1<<4
        // FFTW_DESTROY_INPUT   // 1<<0

        // FFTW_ESTIMATE    // 1<<6
        // FFTW_MEASURE     // 0
        // FFTW_PATIENT     // 1<<5
        // FFTW_EXHAUSTIVE  // 1<<3

        fft_plan  = fftw_plan_dft_3d( n.x,n.y,n.z, V,  Vw, FFTW_FORWARD,  flags );
        ifft_plan = fftw_plan_dft_3d( n.x,n.y,n.z, Vw, V,  FFTW_BACKWARD, flags );
        /*
        https://www.fftw.org/fftw3_doc/Planner-Flags.html
        flags: This parameter provides hints to FFTW about how to optimize the planning process. The flags can include:
           * FFTW_ESTIMATE: This tells FFTW to create a plan quickly without running and measuring actual transforms. The plan might not be optimal, but it’s created quickly, which is useful in scenarios where you need to create the plan fast, or the optimal performance isn't critical.
           * FFTW_MEASURE: This instructs FFTW to run and measure the performance of several different FFT algorithms and choose the best one. This can take some time but usually results in a faster plan.
           * FFTW_PATIENT and FFTW_EXHAUSTIVE are more time-consuming options that further refine the plan by exploring even more possible strategies.
        */
    }

    void prepare_laplace_omp( int flags=-1 ){
        int ntot = n.totprod();
        fftw_init_threads();
        V  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot );
        Vw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntot );
        
        if(flags<0){ flags = FFTW_ESTIMATE | FFTW_PRESERVE_INPUT; }
        
        fftw_plan_with_nthreads(4); 
        fft_plan  = fftw_plan_dft_3d( n.x,n.y,n.z, V,  Vw, FFTW_FORWARD,  flags );
        ifft_plan = fftw_plan_dft_3d( n.x,n.y,n.z, Vw, V,  FFTW_BACKWARD, flags );

    }

    __attribute__((hot))
    void solve_laplace( const double* dens, double* Vout=0 ){
        int ntot =  n.totprod();
        array2fftc( ntot, dens, V );
        fftw_execute(fft_plan);
        long t0 = getCPUticks();
        laplace_reciprocal_kernel( Vw );
        double t = (getCPUticks()-t0)*1e-6; printf( "solve_laplace() ng(%i,%i,%i) T(laplace_reciprocal_kernel)=%g [Mticks] \n", n.x,n.y,n.z, t );
        fftw_execute(ifft_plan);
        if(Vout) fftc2array( ntot, V, Vout );        
    }

    void destroy_laplace(){
        fftw_destroy_plan(fft_plan);
        fftw_destroy_plan(ifft_plan);
        fftw_free(V);
        fftw_free(Vw);
    }

    void destroy_laplace_omp(){
        fftw_destroy_plan(fft_plan);
        fftw_destroy_plan(ifft_plan);
        fftw_free(V);
        fftw_free(Vw);
        fftw_cleanup_threads();

    }

#else

    void prepare_laplace(){
        printf("ERROR: you invoke prepare_laplace() while WITH_FFTW=false => exit()\n" ); exit(0);
    }

    void solve_laplace(){
        printf("ERROR: you invoke solve_laplace() while WITH_FFTW=false => exit()\n" ); exit(0);
    };

    void destroy_laplace(){
        printf("ERROR: you invoke destroy_laplace() while WITH_FFTW=false => exit()\n" ); exit(0);
    }

#endif 

}; // class GridFF

#endif
