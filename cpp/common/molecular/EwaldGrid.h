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


#endif
