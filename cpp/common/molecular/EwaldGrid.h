#ifndef EwaldGrid_h
#define EwaldGrid_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Bspline.h"

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

}; // class GridFF




#ifdef WITH_FFTW

fftw_plan fft_plan;
fftw_plan ifft_plan;
fftw_complex *V=0, *rho=0;


void laplace_reciprocal_kernel(){
    int nx=grid.n.x;
    int ny=grid.n.y;
    int nz=grid.n.z;
    // Solve the Laplace equation in Fourier space
    double kConst  = (2.0 * M_PI / L) * (2.0 * M_PI / L); 
    for (int ix = 0; ix < nx; ix++) {
        double kx = (ix <= ny / 2) ? ix : ix - nx;
        for (int iy = 0; iy < ny; iy++) {
            double ky = (iy <= ny / 2) ? iy : iy - ny;
            for (int iz = 0; iz < nzx; ++iz ) {
                int index = ix * ny * nz + iy * nz + iz;
                double kz = ( iz <= nz / 2) ? iz : iz - nz;

                double k_squared = 1/( (kx * kx + ky * ky + kz * kz) * kConst );

                if (k_squared != 0) {
                    phi_hat[index][0] *= -k_squared; // Real part
                    phi_hat[index][1] *= -k_squared; // Imaginary part
                } else {
                    phi_hat[index][0] = 0.0; // Avoid division by zero (DC component)
                    phi_hat[index][1] = 0.0;
                }
            }
        }
    }
}

    void prepare_laplace(){
        int nx=grid.n.x;
        int ny=grid.n.y;
        int nz=grid.n.z;
        fftw_complex *in, *out;
        V  = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*ny*nz );
        Vw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nx*ny*nz );
        //fft_plan =            fftw_plan_dft_3d( nx,ny,nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
        fft_plan = fftw_plan_dft_3d( nx,ny,nz, V, Vw,   FFTW_FORWARD,  FFTW_ESTIMATE );
        ifft_plan = fftw_plan_dft_3d( nx,ny,nz, Vwt, V, FFTW_BACKWARD, FFTW_ESTIMATE );


    }

    void solve_laplace(){
        fftw_execute(fft_plan);
        laplace_reciprocal_kernel
        fftw_execute(backward_plan);
    }

    void destroy_laplace(){
        fftw_destroy_plan(forward_plan);
        fftw_destroy_plan(backward_plan);
        fftw_free(fft_plan);
        fftw_free(ifft_plan);
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

#endif
