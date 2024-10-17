#ifndef EwaldGrid_h
#define EwaldGrid_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Grid.h"
#include "Bspline.h"
//#include "VecN.h"
#include "IO_utils.h"

#ifdef WITH_OMP
#include <omp.h>
#endif

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

__attribute__((hot))
void fftc2array_mul( int n, const fftw_complex* in, double* out, double f) {
    for (int i=0; i<n; i++ ) { out[i] = in[i][0] * f;}
}

#endif



inline int pbc_ifw(int i, int n){ i++; return (i<n )?  i :  i-n; };
inline int pbc_ibk(int i, int n){ i--; return (i>=0)?  i :  i+n; };

class EwaldGrid : public GridShape { public: 

double Qtot=0.0;
double Qabs=0.0;
double Qtot_g=0.0;
double Qabs_g=0.0;
Vec3d dipole = Vec3dZero;
double* V_work  = 0; 
double* vV_work = 0;  

bool bCubicPBCIndexesDone = false;
Quat4i xqs_o3[4];
Quat4i yqs_o3[4];
Quat4i zqs_o3[4];

bool bQuinticPBCIndexesDone = false;
Vec6i xqs_o5[6];
Vec6i yqs_o5[6];
Vec6i zqs_o5[6];

__attribute__((hot)) 
inline void project_atom_on_grid_linear( const Vec3d pi, const double qi, double* dens ) const {
    //printf("project_atom_on_grid() pi(%g,%g,%g) q=%g \n", pi.x, pi.y, pi.z, qi );
    const Vec3d gp = diCell.dot( pi-pos0 );
    const int ix = (int) gp.x; 
    const int iy = (int) gp.y;
    const int iz = (int) gp.z;
    const double tx = gp.x - ix;
    const double ty = gp.y - iy;
    const double tz = gp.z - iz;
    const double mx = 1-tx;
    const double my = 1-ty;
    const double mz = 1-tz;
    //printf("project_atom_on_grid() pi(%g,%g,%g) q=%g \n", pi.x, pi.y, pi.z, qi );
    const int nxy = n.x * n.y;
    //int ii=0;

    int ig = iz*nxy + iy*n.x + ix;

    const double f00 = qi*my*mx;
    const double f01 = qi*my*tx;
    const double f10 = qi*ty*mx;
    const double f11 = qi*ty*tx;

    dens[ig          ] += f00*mz;
    dens[ig+1        ] += f01*mz;
    dens[ig  +n.x    ] += f10*mz;
    dens[ig+1+n.x    ] += f11*mz;
    dens[ig      +nxy] += f00*tz;
    dens[ig+1+   +nxy] += f01*tz;
    dens[ig  +n.x+nxy] += f10*tz;
    dens[ig+1+n.x+nxy] += f11*tz;
    
}

inline void split_atoms_parallel( int na, Vec3d* apos, double Rcut ){
    // -- split atoms so that we are sure they do no overlap within Rcut => we can project them on grid without fear of memory write collision
}

__attribute__((hot)) 
double project_atom_on_grid_cubic( const Vec3d pi, const double qi, double* dens ) const {
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
    double Qsum=0;
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
                double qi = qbyz * bx.array[dx]; 
                dens[ig] += qi;
                Qsum += qi;
            }
        }
    }
    return Qsum;
}

__attribute__((hot)) 
double project_atom_on_grid_cubic_pbc(const Vec3d pi, const double qi, double* dens) const {
    // Convert atomic position to grid position
    const Vec3d gp = diCell.dot(pi - pos0);
    int ix = (int) gp.x;     if(gp.x<0) ix--;
    int iy = (int) gp.y;     if(gp.y<0) iy--;
    int iz = (int) gp.z;     if(gp.z<0) iz--;
    const double tx = gp.x - ix;
    const double ty = gp.y - iy;
    const double tz = gp.z - iz;

    // Calculate the B-spline weights
    const Quat4d bx = Bspline::basis(tx);
    const Quat4d by = Bspline::basis(ty);
    const Quat4d bz = Bspline::basis(tz);

    // Retrieve the number of grid points in each direction
    const int nxy = n.x * n.y;

    // Pre-calculate periodic boundary condition indices for each dimension
    ix=modulo(ix-1,n.x); const Quat4i xqs = choose_inds_pbc_3(ix, n.x, xqs_o3 );  // Assuming you pre-calculate xqs, yqs, zqs
    iy=modulo(iy-1,n.y); const Quat4i yqs = choose_inds_pbc_3(iy, n.y, yqs_o3 );
    iz=modulo(iz-1,n.z); const Quat4i zqs = choose_inds_pbc_3(iz, n.z, zqs_o3 );
    //printf( "project_atom_on_grid_cubic_pbc() ixyz(%i,%i,%i) xqs(%i,%i,%i,%i) yqs(%i,%i,%i,%i) nxyz(%i,%i,%i)\n", ix,iy,iz,  xqs.x,xqs.y,xqs.z,xqs.w,   yqs.x,yqs.y,yqs.z,yqs.w,   n.x,n.y,n.z );
    // Loop over the B-spline grid contributions
    double Qsum=0;
    for (int dz = 0; dz < 4; dz++) {
        const int gz = zqs.array[dz];
        const int iiz = gz * nxy;
        for (int dy = 0; dy < 4; dy++) {
            const int gy = yqs.array[dy];
            const int iiy = iiz + gy * n.x;
            const double qbyz = qi * by.array[dy] * bz.array[dz];
            for (int dx = 0; dx < 4; dx++) {
                const int gx = xqs.array[dx];
                const int ig = gx + iiy;
                double qi = qbyz * bx.array[dx];
                dens[ig] += qi;
                Qsum += qi;
            }
        }
    }
    return Qsum;
}



__attribute__((hot)) 
double project_atom_on_grid_quintic( const Vec3d pi, const double qi, double* dens ) const {
    //printf("project_atom_on_grid() pi(%g,%g,%g) q=%g \n", pi.x, pi.y, pi.z, qi );
    const Vec3d gp = diCell.dot( pi-pos0 );

    //printf("project_atom_on_grid() pi(%7.3f,%7.3f,%7.3f) gp(%7.3f,%7.3f,%7.3f) \n", pi.x, pi.y, pi.z, gp.x, gp.y, gp.z );
    const int ix = (int) gp.x;
    const int iy = (int) gp.y;
    const int iz = (int) gp.z;
    const double tx = gp.x - ix;
    const double ty = gp.y - iy;
    const double tz = gp.z - iz;
    // ToDo: Periodic Boundary Conditions for points which are close to the boundary
    const Vec6d bx = Bspline::basis5(tx);
    const Vec6d by = Bspline::basis5(ty);
    const Vec6d bz = Bspline::basis5(tz);

    //printf("project_atom_on_grid() pi(%g,%g,%g) q=%g \n", pi.x, pi.y, pi.z, qi );
    const int nxy = n.x * n.y;
    //int ii=0;
    double Qsum=0;
    for (int dz = 0; dz < 6; dz++) {
        const int gz  = iz + dz - 2;
        const int iiz = gz*nxy;
        for (int dy = 0; dy < 6; dy++) {
            const int gy  = iy + dy - 2;
            const int iiy = iiz + gy*n.x;
            const double qbyz = qi * by.array[dy] * bz.array[dz];
            for (int dx = 0; dx < 6; dx++) {
                const int gx = ix + dx - 2;
                const int ig = gx + iiy;
                double qi = qbyz * bx.array[dx];
                dens[ig] += qi;
                Qsum += qi;
            }
        }
    }
    return Qsum;
}

__attribute__((hot)) 
double project_atom_on_grid_quintic_pbc(const Vec3d pi, const double qi, double* dens) const {
    // Convert atomic position to grid position
    const Vec3d gp = diCell.dot(pi - pos0);
    int ix = (int) gp.x;    if(gp.x<0) ix--;
    int iy = (int) gp.y;    if(gp.y<0) iy--;
    int iz = (int) gp.z;    if(gp.z<0) iz--;
    const double tx = gp.x - ix;
    const double ty = gp.y - iy;
    const double tz = gp.z - iz;

    // Calculate the quintic B-spline weights
    const Vec6d bx = Bspline::basis5(tx);
    const Vec6d by = Bspline::basis5(ty);
    const Vec6d bz = Bspline::basis5(tz);

    // Retrieve the number of grid points in each direction
    const int nxy = n.x * n.y;

    // Pre-calculate periodic boundary condition indices for each dimension
    ix=modulo(ix-2,n.x); const Vec6i xqs = choose_inds_pbc_5(ix, n.x, xqs_o5 );  // Assuming you pre-calculate xqs, yqs, zqs
    iy=modulo(iy-2,n.y); const Vec6i yqs = choose_inds_pbc_5(iy, n.y, yqs_o5 );
    iz=modulo(iz-2,n.z); const Vec6i zqs = choose_inds_pbc_5(iz, n.z, zqs_o5 );

    // Loop over the B-spline grid contributions
    double Qsum=0;
    for (int dz = 0; dz < 6; dz++) {
        const int gz = zqs.array[dz];
        const int iiz = gz * nxy;
        for (int dy = 0; dy < 6; dy++) {
            const int gy = yqs.array[dy];
            const int iiy = iiz + gy * n.x;
            const double qbyz = qi * by.array[dy] * bz.array[dz];
            for (int dx = 0; dx < 6; dx++) {
                const int gx = xqs.array[dx];
                const int ig = gx + iiy;
                double qi = qbyz * bx.array[dx];
                dens[ig] += qi;
                Qsum += qi;
            }
        }
    }
    return Qsum;
}


__attribute__((hot)) 
void project_atoms_on_grid_linear( int na, const Vec3d* apos, const double* qs, double* dens, bool bPBC=true ) {
    //printf("project_atoms_on_grid_linear() na=%i ns(%i,%i,%i) pos0(%g,%g,%g)\n", na, n.x,n.y,n.z, pos0.x,pos0.y,pos0.z );
    for (int ia=0; ia<na; ia++){
        project_atom_on_grid_linear( apos[ia], qs[ia], dens );
    }
}

__attribute__((hot)) 
void project_atoms_on_grid_cubic( int na, const Vec3d* apos, const double* qs, double* dens, bool bPBC=true ) {
    if( (!bCubicPBCIndexesDone) && bPBC ){ 
        make_inds_pbc(n.x, xqs_o3); 
        make_inds_pbc(n.y, yqs_o3); 
        make_inds_pbc(n.z, zqs_o3);
        //for(int i=0; i<4; i++){  printf( "xqs_o3[%i]{%i,%i,%i,%i}\n",       i, xqs_o3[i].x,xqs_o3[i].y,xqs_o3[i].z,xqs_o3[i].w ); };  
        bCubicPBCIndexesDone=true;
    }
    Qtot_g=0;
    Qabs_g=0;
    //printf("project_atoms_on_grid_cubic() na=%i ns(%i,%i,%i) pos0(%g,%g,%g)\n", na, n.x,n.y,n.z, pos0.x,pos0.y,pos0.z );
    for (int ia=0; ia<na; ia++){ 
        double qg;
        if(bPBC){ qg=project_atom_on_grid_cubic_pbc( apos[ia], qs[ia], dens ); } 
        else    { qg=project_atom_on_grid_cubic    ( apos[ia], qs[ia], dens ); }
        //printf( "project_atoms_on_grid_quintic() qs[%i]=%g   qg=%g |qi-qg|=%20.10e \n", ia, qs[ia], qg, qs[ia]-qg );
        Qtot_g += qg;
        Qabs_g += fabs(qg);    
    }
}

__attribute__((hot)) 
void project_atoms_on_grid_quintic( int na, const Vec3d* apos, const double* qs, double* dens, bool bPBC=true ) {
    //printf("project_atoms_on_grid_quintic() na=%i ns(%i,%i,%i) pos0(%g,%g,%g)\n", na, n.x,n.y,n.z, pos0.x,pos0.y,pos0.z );
    if( (!bQuinticPBCIndexesDone) && bPBC ){ 
        make_inds_pbc_5(n.x, xqs_o5); 
        make_inds_pbc_5(n.y, yqs_o5); 
        make_inds_pbc_5(n.z, zqs_o5);
        bQuinticPBCIndexesDone = true;
    }
    Qtot_g=0;
    Qabs_g=0;
    for (int ia=0; ia<na; ia++){
        double qg;
        if(bPBC){ qg=project_atom_on_grid_quintic_pbc( apos[ia], qs[ia], dens ); }
        else    { qg=project_atom_on_grid_quintic    ( apos[ia], qs[ia], dens ); }
        //printf( "project_atoms_on_grid_quintic() qs[%i]=%g   qg=%g  |qi-qg|=%20.10e  \n", ia, qs[ia], qg, qs[ia]-qg );
        Qtot_g += qg;
        Qabs_g += fabs(qg);
    }
    //if(bPBC){ for (int ia=0; ia<na; ia++){ project_atom_on_grid_quintic_pbc( apos[ia], qs[ia], dens ); } }
    //else    { for (int ia=0; ia<na; ia++){ project_atom_on_grid_quintic    ( apos[ia], qs[ia], dens ); } }  
}

int setup( Vec3d pos0_, Mat3d dCell_, Vec3i ns_, bool bPrint=false ){
    n     = ns_;
    pos0  = pos0_;
    dCell = dCell_;
    updateCell_2();
    if(bPrint){printCell();}
    return n.totprod();
}

void projectAtoms( int na, Vec3d* apos, double* qs, double* dens, int order ){
    long t0 = getCPUticks();
    dipole=Vec3dZero;
    Qabs=0;
    Qtot=0;
    for (int ia=0; ia<na; ia++){ dipole.add_mul(apos[ia],qs[ia]); Qtot+=qs[ia]; Qabs+=fabs(qs[ia]); };  // dipole for slab correction and debugging
    //printf( "EwaldGrid::projectAtoms() na=%i dipole(%g,%g,%g) \n", na, dipole.x,dipole.y,dipole.z   );
    switch(order){
        case 1: project_atoms_on_grid_linear ( na, apos, qs, dens ); break;
        case 2: project_atoms_on_grid_cubic  ( na, apos, qs, dens ); break;
        case 3: project_atoms_on_grid_quintic( na, apos, qs, dens ); break;
        default: printf("ERROR in projectAtoms() order=%i NOT IMPLEMETED !!! \n", order ); exit(0); break;
    }
    //if( bQuintic ){  }
    //else          { W.gewald.project_atoms_on_grid        ( na, (Vec3d*)apos, qs, dens ); }
    double t = (getCPUticks()-t0)*1e-6; printf( "EwaldGrid::projectAtoms(order=%i) na=%i ng(%i,%i,%i) T(project_atoms_on_grid)=%g [Mticks] \n", order, na, n.x,n.y,n.z, t );
}

__attribute__((hot))
double laplace_real( double* Vin, double* Vout, double cSOR ){
    int nxy = n.x * n.y;
    const double fac = 1/6.0;
    double err2 = 0.0; 
    for (int iz = 1; iz < n.z-1; ++iz ) {
        for (int iy = 1; iy < n.y-1; iy++) {
            for (int ix = 1; ix < n.x-1; ix++) {
                const int i = iz*nxy + iy*n.z + ix;
                double vi = 
                    Vin[ i-1   ] + Vin[ i+1   ] + 
                    Vin[ i-n.x ] + Vin[ i+n.x ] + 
                    Vin[ i-nxy ] + Vin[ i+nxy ];
                vi*=fac;
                const double vo = Vin[i];
                vi += (vi-vo)*cSOR; 
                const double dv = vi - vo;
                err2 += dv*dv;
                Vout[i] = vi;
            }
        }
    }
    return err2;
}

__attribute__((hot))
double laplace_real_pbc( double* Vin, double* Vout, double cSOR=0.0 ){
    int nxy = n.x * n.y;
    const double fac = 1/6.0;
    double err2 = 0.0; 
    
    for (int iz = 0; iz < n.z; ++iz ) {
        const int iiz =          iz      *nxy;
        const int ifz =  pbc_ifw(iz, n.z)*nxy;
        const int ibz =  pbc_ibk(iz, n.z)*nxy;
        for (int iy = 0; iy < n.y; iy++) {
            const int iiy =          iy      *n.x;
            const int ify =  pbc_ifw(iy, n.y)*n.x;
            const int iby =  pbc_ibk(iy, n.y)*n.x;
            for (int ix = 0; ix < n.x; ix++) {
                const int ifx =  pbc_ifw(ix, n.x);
                const int ibx =  pbc_ibk(ix, n.x);
                double vi = 
                    Vin[ ibx + iiy + iiz ] + Vin[ ifx + iiy + iiz ] + 
                    Vin[ ix  + iby + iiz ] + Vin[ ix  + ify + iiz ] + 
                    Vin[ ix  + iiy + ibz ] + Vin[ ix  + iiy + ifz ];
                vi*=fac;
                const int i = ix + iiy + iiz;
                const double vo = Vin[ i ];
                vi += (vi-vo)*cSOR; 
                const double dv = vi - vo;
                err2 += dv*dv;
                Vout[i] = vi;
            }
        }
    }
    return err2;
}

__attribute__((hot))
int laplace_real_loop( double* V, int nmaxiter=1000, double tol=1e-6, bool bPBC=true, double cSOR=0.0 ){
    //bPBC = false;
    int ntot = n.totprod();
    double* V_ =0;
    if(V_work ){ V_ = V_work;  }else{ double* V_ = new double[ntot]; };
    //printf("laplace_real_loop(bPBC=%i) nmaxiter=%i tol=%g @V=%li @V_=%li \n", bPBC,  nmaxiter, tol, (long)V, (long)V_ );
    int iter=0;
    for(iter=0; iter<nmaxiter; iter++){ 
        if(bPBC){ laplace_real_pbc( V, V_, cSOR ); }
        else    { laplace_real    ( V, V_, cSOR ); }
        _swap( V, V_ );
    }
    if(iter%2==1){ for(int i=0; i<ntot; i++){ V_[i] = V[i]; } _swap( V, V_ ); }
    //printf("laplace_real_loop(bPBC=%i) DONE  iter=%i @V=%li @V_=%li \n", bPBC, iter, tol, (long)V, (long)V_ );
    if(V_work ==0) delete[] V_;
    return iter;
}

__attribute__((hot))
int laplace_real_loop_inert( double* V, int nmaxiter=1000, double tol=1e-6, bool bPBC=true, double cSOR=0.0, double cV=0.5 ){
    //bPBC = false;
    int ntot = n.totprod();
    double* V_ =0;
    double* vV =0; 
    if(V_work ){ V_ = V_work;  }else{ double* V_ = new double[ntot]; };
    if(vV_work){ vV = vV_work; }else{ double* vV = new double[ntot]; };
    //printf("laplace_real_loop_inert(bPBC=%i) nmaxiter=%i tol=%g @V=%li @V_=%li \n", bPBC,  nmaxiter, tol, (long)V, (long)V_ );
    int iter=0;
    for(iter=0; iter<nmaxiter; iter++){ 
        //printf("laplace_real_loop_inert()[iter=%i]\n", iter );
        if(bPBC){ laplace_real_pbc( V, V_, cSOR ); }
        else    { laplace_real    ( V, V_, cSOR ); }
        for(int i=0; i<ntot; i++){ 
            double v = V_[i]-V[i];
            if(iter>0){ v = v*cV + vV[i]*(1-cV); }
            vV[i] = v; 
            V_[i] = V[i] + v;
        }
        _swap( V, V_ );
    }
    if(iter%2==1){ for(int i=0; i<ntot; i++){ V_[i] = V[i]; } _swap( V, V_ ); }
    //printf("laplace_real_loop(bPBC=%i) DONE  iter=%i @V=%li @V_=%li \n", bPBC, iter, tol, (long)V, (long)V_ );
    if(V_work ==0) delete[] V_;
    if(vV_work==0) delete[] vV;
    return iter;
}

void slabPotential( int nz_slab, double* Vin, double* Vout ){ // 
    // NOTE: the cell must be sufficiently large in z-direction, recomanded at lease 2 x as big as in x,y , i.e. L_z > 2 * max( L_x, L_y ) 
    // ToDo: We assume the cell is neutral, if cell as non-zero net charge there should be another term, described in the paper  
    //  See article: https://pubs.aip.org/aip/jcp/article/111/7/3155/294442/Ewald-summation-for-systems-with-slab-geometry
    //double Vol       = getVolume() * ( nz_slab /(n.z) );
    double Vol       = getVolume();
    double dz = dCell.c.norm();   // lenght of grid step along z-axis
    Vec3d  hz = dCell.c *(1/dz);  // normalized direction of z-axis 
    //double dVcor  =  16.0 * COULOMB_CONST * hz.dot(dipole)/Vol;
    double dVcor  = 4.0 *M_PI * COULOMB_CONST * hz.dot(dipole)/Vol;
    double Vcor0 = -dVcor * cell.c.norm()/2; // 0.5*Lz
    printf( "EwaldGrid::slabPotential() ns(%i,%i,%i) Ls(%g,%g,%g) nz_slab=%i \n", n.x,n.y,n.z, cell.a.norm(), cell.b.norm(), cell.c.norm(), nz_slab );
    printf( "EwaldGrid::slabPotential() dipole(%g,%g,%g) dz=%g [A] Vol=%g[A^3] dVcor=%g [eV/A] \n", dipole.x,dipole.y,dipole.z, dz, Vol, dVcor );
    for (int iz=0; iz<nz_slab; iz++ ) {
        double Vcor_z = Vcor0 + dVcor * (iz*dz);
        for (int iy=0; iy<n.y; iy++) { 
            const int    iyz=(iz*n.y + iy)*n.x;
            for (int ix=0; ix<n.x; ix++) {
                const int i      = ix + iyz;
                Vout[i] = Vin[i] + Vcor_z;
                //Vout[i] = Vin[i] + Vcor_z*0.0000000000; // Debug
            }
        }
    }
}

#ifdef WITH_FFTW

fftw_plan    fft_plan;
fftw_plan    ifft_plan;
fftw_complex *Vw=0,*V=0;

__attribute__((hot))
void laplace_reciprocal_kernel( fftw_complex* VV ){
    const int nx=n.x;
    const int ny=n.y;
    const int nz=n.z;
    const double freq_x = (2.0 * M_PI) / cell.a.norm();
    const double freq_y = (2.0 * M_PI) / cell.b.norm();
    const double freq_z = (2.0 * M_PI) / cell.c.norm();
    //printf("laplace_reciprocal_kernel() nxyz(%i,%i,%i) Ls(%g,%g,%g) freq(%g,%g,%g) \n", nx,ny,nz,  cell.a.norm(),cell.b.norm(),cell.c.norm(), freq_x,freq_y,freq_z );
    const int nx2=n.x/2;
    const int ny2=n.y/2;
    const int nz2=n.z/2;
    for (int iz=0; iz<nz; iz++ ) {
        const double kz = ( (iz<=nz2) ? iz : iz-nz ) * freq_z;
        for (int iy=0; iy<ny; iy++) { 
            const double ky   = ( (iy<=ny2) ? iy : iy-ny ) * freq_y;
            const double k2yz = ky*ky + kz*kz;
            const int    iyz=(iz*ny + iy)*nx;
            for (int ix=0; ix<nx; ix++) {
                const double kx  =  ( (ix<=nx2) ? ix : ix-nx ) * freq_x;
                const double k2  = kx*kx + k2yz;
                const int i      = ix + iyz;
                if ( k2 > 1e-32 ){
                    const double invk2 = 1.0/k2;
                    VV[i][0] *= invk2; // Real part
                    VV[i][1] *= invk2; // Imaginary part
                    //VV[i][0] *= k2; // Real part
                    //VV[i][1] *= 0; // Imaginary part
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

        //printf("prepare_laplace() nxyz(%i,%i,%i) flags=%i \n", n.x, n.y, n.z, flags );

        //fft_plan  = fftw_plan_dft_3d( n.x,n.y,n.z, V,  Vw, FFTW_FORWARD,  flags );
        //ifft_plan = fftw_plan_dft_3d( n.x,n.y,n.z, Vw, V,  FFTW_BACKWARD, flags );

        fft_plan  = fftw_plan_dft_3d( n.z,n.y,n.x, V,  Vw, FFTW_FORWARD,  flags );
        ifft_plan = fftw_plan_dft_3d( n.z,n.y,n.x, Vw, V,  FFTW_BACKWARD, flags );
        /*
        https://www.fftw.org/fftw3_doc/Planner-Flags.html
        flags: This parameter provides hints to FFTW about how to optimize the planning process. The flags can include:
           * FFTW_ESTIMATE: This tells FFTW to create a plan quickly without running and measuring actual transforms. The plan might not be optimal, but it’s created quickly, which is useful in scenarios where you need to create the plan fast, or the optimal performance isn't critical.
           * FFTW_MEASURE: This instructs FFTW to run and measure the performance of several different FFT algorithms and choose the best one. This can take some time but usually results in a faster plan.
           * FFTW_PATIENT and FFTW_EXHAUSTIVE are more time-consuming options that further refine the plan by exploring even more possible strategies.
        */
    }

/*
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
*/
    __attribute__((hot))
    void solve_laplace( const double* dens, double* Vout=0 ){
        int ntot =  n.totprod();
        array2fftc( ntot, dens, V );
        fftw_execute(fft_plan);
        long t0 = getCPUticks();
        laplace_reciprocal_kernel( Vw );
        double t = (getCPUticks()-t0)*1e-6; 
        //printf( "solve_laplace() ng(%i,%i,%i) T(laplace_reciprocal_kernel)=%g [Mticks] \n", n.x,n.y,n.z, t );
        fftw_execute(ifft_plan);
        //if(Vout) fftc2array( ntot, V, Vout );

        double dV      = dCell.a.norm() * dCell.b.norm() * dCell.c.norm();
        double scEwald = COULOMB_CONST * 4.0 * M_PI / (ntot*dV);

        if(Vout) fftc2array_mul( ntot, V, Vout, scEwald );
        //if(Vout) VecN::mul( ntot, scEwald, const double* b, double* out ){
    }

    void destroy_laplace(){
        fftw_destroy_plan(fft_plan);
        fftw_destroy_plan(ifft_plan);
        fftw_free(V);
        fftw_free(Vw);
    }

    void solve_laplace_macro( double* dens, int nz_slab, double* Vout, bool bPrepare=true, bool bDestroy=true, int flags=-1, bool bOMP=false, int nBlur=0, double cSOR=0, double cV=0.95 ){
        long t0 = getCPUticks();
        //if(bPrepare){ if(bOMP){ prepare_laplace_omp( flags );} else { prepare_laplace( flags );} }
        if(bPrepare){ prepare_laplace( flags ); }

        bool bSlab = (nz_slab>0);
        double* Vtmp = Vout;
        if(bSlab){ Vtmp = new double[n.totprod()]; }

        long t1 = getCPUticks();
        solve_laplace( dens, Vtmp );
        long t2 = getCPUticks();
        if(bDestroy){ destroy_laplace( ); }
        long t3=0,t4=0;
        if(nBlur>0){
            int ntot = n.totprod();
            _allocIfNull( V_work,  ntot );
            _allocIfNull( vV_work, ntot );
            t3 = getCPUticks();
            if( cV<-1.0 ){ laplace_real_loop      ( Vtmp, nBlur, 1e-32, true, cSOR     ); }
            else         { laplace_real_loop_inert( Vtmp, nBlur, 1e-32, true, cSOR, cV ); }
            t4 = getCPUticks();
        }

        { // save debug numpy
            bool bSaveNPY = true;
            bool bSaveXSF = false;
            if(bSaveNPY){
                Vec3i sh{n.z,n.y,n.x};
                save_npy( "debug_Ewald_V.npy",    3, (int*)&sh, (char*)Vtmp );
                save_npy( "debug_Ewald_dens.npy", 3, (int*)&sh, (char*)dens );
            }
            // if(bSaveXSF){
            //     gBS.saveXSF( "debug_VCoul.xsf",   VCoul, 1,0 );
            // }
        }

        if( bSlab ){ 
            slabPotential( nz_slab, Vtmp, Vout ); 
            delete[] Vtmp;
        }

        int nCPUs=1;
        #ifdef WITH_OMP
            nCPUs=omp_get_max_threads();
        #endif
    
        printf( "EwaldGrid::solve_laplace_macro() DONE flags=%i nCPUs=%i n(%i,%i,%i) T(prepare_laplace)= %g [Mticks] T(solve_laplace)= %g [Mticks] T(laplace_real_loop)= %g [Mticks]\n", flags, nCPUs, n.x,n.y,n.z, (t1-t0)*1e-6, (t2-t1)*1e-6, (t4-t3)*1e-6 );
        // if(bDestroy){ if(bOMP){ destroy_laplace_omp( ); } else { destroy_laplace    ( ); } }
    }

    void potential_of_atoms( int nz, double* VCoul, int natoms, Vec3d* apos, double* qs, int order=3, int nBlur=4, double cV=0.95, bool bClearCharge=true ){
        printf("EwaldGrid::potential_of_atoms() natoms=%i order=%i nBlur=%i cV=%g \n", natoms, order, nBlur, cV );
        int ntot = n.totprod();
        double* dens  = new double[ ntot ];
        if(bClearCharge){ for(int i=0; i<ntot; i++){ dens[i]=0; }; }
        projectAtoms( natoms, apos, qs, dens, order );
        solve_laplace_macro( dens, nz, VCoul, true, true, -1, false, nBlur, 0, cV );
        delete [] dens;
        //delete [] qs;
    }


#else

    void prepare_laplace( int flags=-1 ){
        printf("ERROR: you invoke prepare_laplace() while WITH_FFTW=false => exit()\n" ); exit(0);
    }

    void solve_laplace( const double* dens, double* Vout=0 ){
        printf("ERROR: you invoke solve_laplace() while WITH_FFTW=false => exit()\n" ); exit(0);
    };

    void destroy_laplace(){
        printf("ERROR: you invoke destroy_laplace() while WITH_FFTW=false => exit()\n" ); exit(0);
    }

#endif 

}; // class GridFF

#endif
