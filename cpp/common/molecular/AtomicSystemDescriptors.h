
#ifndef AtomicSystemDescriptors_h
#define AtomicSystemDescriptors_h

#include <vector>
#include "Vec3.h"
#include "quaternion.h"
#include "Atoms.h"
//#include "NBFF.h"
//#include "MMFFparams.h"
//#include "Forces.h"


/*
Something like SOAP but using linear (or exponential) projection of bond-lengths and angles on a grid
*/


/**
 * @brief a 1D linear function onto an array of values.
 * @param x The position of the linear function.
 * @param w The weight of the linear function.
 * @param n The size of the target array.
 * @param target The array to project the linear function onto.
 */
void project_linear_1d( double x, double w, int n, double* target ){
    int i  = (int)x;
    double dx = x-i;
    double mx = 1-dx;
    target[i  ]= mx*w;
    target[i+1]= dx*w;
}

/**
 * @brief Projects a 2D point onto a 2D grid using linear interpolation.
 * @param x The x-coordinate of the point to project.
 * @param y The y-coordinate of the point to project.
 * @param w The weight to assign to the projected point.
 * @param nx The number of grid points in the x-direction.
 * @param ny The number of grid points in the y-direction.
 * @param target A pointer to an array of size [nx*ny] that will hold the projected values.
 */
void project_linear_2d( double x, double y, double w, int nx, int ny, double* target ){
    int ix  = (int)x;
    double dx = x-ix;
    double mx = 1-dx;
    int iy  = (int)y;
    double dy = y-iy;
    double my = 1-dy;
    int ix0 = iy*nx;
    target[ix  + ix0   ]= mx*my*w;
    target[ix+1+ ix0   ]= dx*my*w;
    target[ix  + ix0+nx]= mx*dy*w;
    target[ix+1+ ix0+nx]= dx*dy*w;
}


class AtomicSystemDescriptors : public Atoms { public:

    int nmaxneighs; // maximum number of neighbors
    int ntypes;     // number of atom types
    int ntype2;     // number of remapped atom types
    int ntype_comp; // number of atom type combinations

    int * nneigh=0;     // [natoms ] number of neighbors
    int * neighs=0;     // [natoms*nmaxneighs] list of neighbors
    int * type_remap=0;   // [ntypes] maps atom types to a smaller set of types so that we reduce number of descriptors
    int * type_compose=0; // 
    int * ntype_comp=0;   // [ntype_comp] number of atom type combinations

    DescriptorSettings* settings=0;
    
    float * bond_coeffs=0; // coefficient for bond-lengths projection (for each type)
    float * angle_coeffs=0; // coefficient for bond-angles angular projection (for each type combination)


    int nang=0, nbond=0;
    double Rmin=1.0, Rmax=2.718281828;

    //Mat3d* typeRotations=0;   // [ntype_comp] rotation matrices for each atomic type combination
    //Vec3d* typeRotEigVals=0;  // [ntype_comp] eigenvalues of rotation matrices for each atomic type combination

    //void project_bonds(){} 
    //void project_angles(){}

    inline int bondTypeIndex( int i, int j ){
        if(i<j) _swap(i,j);
        return return i*(i+1)/2 + j;
    }

    /**
     * @brief Projects angles and bonds for a given atomic system onto grid. It can decompose the projections by atom type. It can also project 2d map of angles and bonds simultaneously.
     * @param sys The atomic system to project angles and bonds for.
     * @param neighs An array of neighbor indices for each atom in the system.
     * @param bonds output array of bond distances for each atom in the system.
     * @param angs output array of angles between each pair of bonded atoms in the system.
     * @param data2d output array of 2D projections of angles and bond distances for each pair of bonded atoms in the system.
     * @param split_types A boolean indicating whether to split the projections by atom type.
     */
    void project_angles_and_bonds( const Atoms& sys, int* neighs, double* bonds, double* angs, double* data2d, bool split_types=false ){ 
        double Rmax2 = Rmax*Rmax;
        Quat4d hs[nmaxneighs];
        int    ts[nmaxneighs];
        double udang  = nang/(2*M_PI);
        double udbond = nbond/(Rmax-Rmin);
        for( int ia=0; ia<natoms; ia++){
            Vec3d pi = sys.apos[ia];
            int ti   = sys.atypes[ia];
            ti=type_remap[ti];
            int*  ngs = neighs[ia*nmaxneighs];
            // evaluate direction vectors and distances
            for(int ing=0; ing<nmaxneighs; ing++){
                int ja = ngs[ing];
                if(ja<0) break;
                Vec3d dj = sys.apos[ja] - pi;
                int tj = sys.atypes[ja];
                tj=type_remap[tj];
                ts[ing]=tj;
                hs[ing].e=dj.normalize();
                hs[ing].f=dj;

                if(bonds){ 
                    project_linear_1d( hi.e, nbond, bonds ); 
                    if(split_types){
                        btyp = bondTypeIndex( ti, tj );
                        project_linear_1d( hi.e, nbond, bonds + btyp*nbond ); 
                    }
                }
                //double r2 = dj.norm2();
                //if( r2<Rmax2 ){    double r = sqrt(r2);}  
            }
            // project angles and bonds
            for(int ing=0; ing<nmaxneighs; ing++){
                Quat4d hi = hs[ing];
                for(int jng=0; jng<ing; jng++){
                    Quat4d hj = hs[jng];
                    double ang = hi.f.getAngle_unitary( hj.f );
                    ang = (ang+M_PI)     *udang;
                    double b1 = (hi.e-Rmin)*udbond;
                    double b2 = (hj.e-Rmin)*udbond;
                    if(angs){
                        project_linear_1d( ang, nbond, angs+angtyp*nang ); 
                        if(split_types){
                            int angtyp = ti*ntype2 + bondTypeIndex( ts[ing], ts[jng] );
                            project_linear_1d( ang, nbond, angs+angtyp*nang ); 
                        }
                    }
                    if(data2d){
                        project_linear_2d( ang, b1, 1, nbond, nang, data2d );
                        project_linear_2d( ang, b2, 1, nbond, nang, data2d );
                    }
                }
            }
        }

    }


    /**
     * @brief Calculates the rotation matrices which characterize angular distribution of atoms of each type around the reference point (center of mass). This is done by calculating the quarupole matrices for each atom and summing them for each type. The quarupole matrices are then diagonalized to obtain the rotation matrices.
     * @param sys The atomic system for which to calculate the rotation matrices
     * @param p0 The reference point for the quarupole calculation.
     * @param typeRotations The rotation matrices for each type are stored here. 
     * @param typeRotEigVals The eigenvalues of the rotation matrices for each type are stored here (i.e. length of the principal axes of the ellipsoid ).
     * @param bClear If true, the quarupole matrices for each type are set to zero before calculation
     * @param bDiagonalize If true, diagonalizes the quarupole matrices for each type after the summing.
     */
    void evalTypeOrientations( const Atoms& sys, Vec3d p0, Mat3d typeRotations, Vec3d* typeRotEigVals, bool bClear=true, bool bDiagonalize=true ){
        if(bClear) for(int i=0; i<(ntype2+ntype_comp); i++){ typeRotations[i].set(0.0); }
        // sum quarupole matrices for each type
        for( int ia=0; ia<natoms; ia++){
            int ityp  = sys.atypes[ia];
            int ityp2 = type_remap[ityp]; 
            Vec3d p  = sys.apos[ia] - p0;
            Mat3d dq; dq.setOuter( p, p );
            typeRotations[ityp2].add( dq );
        }
        // sum comoosed types
        int i=0;
        for( int ityp=0; ityp<ntype_comp; ityp++){
            int nt = ntype_comp[ityp];
            int iout = type_compose[i];
            for( int j=0; j<nt; j++ ){
                int iin = type_compose[i];
                typeRotations[iout].add( typeRotations[iin] );
                i++;
            }
        }
        // diagonalize
        if(bDiagonalize){
            for( int ityp=0; ityp<ntype_comp; ityp++){
                Vec3d eig; Mat3d rot;
                typeRotations[ityp].eigenvals( eig );
                typeRotEigVals[ityp] = eig;
                eigenvec( eig.a, rot.a );
                eigenvec( eig.c, rot.c );
                rot.b.set_cross(rot.c,rot.a);
            }
        }

    }


}; // class FitREQ


#endif
