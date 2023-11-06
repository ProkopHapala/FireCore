
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

class DescriptorSettings{


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


    double Rmin=1.0, Rmax=2.718281828;

    //Mat3d* typeRotations=0;   // [ntype_comp] rotation matrices for each atomic type combination
    //Vec3d* typeRotEigVals=0;  // [ntype_comp] eigenvalues of rotation matrices for each atomic type combination

    //void project_bonds(){} 
    //void project_angles(){}

    void project_angles_and_bonds( const Atoms& sys, int* neighs, int nang, int nb, double* bonds, double* angs,  double* data2d ){
        double Rmax2 = Rmax*Rmax;
        Quat4d hs[nmaxneighs];
        for( int ia=0; ia<natoms; ia++){
            Vec3d pi = sys.apos[ia];
            int*  ngs = neighs[ia*nmaxneighs];
            // evaluate direction vectors and distances
            for(int ing=0; ing<nmaxneighs; ing++){
                int ja = ngs[ing];
                if(ja<0) break;
                Vec3d dj = sys.apos[ja] - pi;

                hs[ing].e=dj.normalize();
                hs[ing].f=dj;
                //double r2 = dj.norm2();
                //if( r2<Rmax2 ){    double r = sqrt(r2);}  
            }
            // project angles and bonds
            for(int ing=0; ing<nmaxneighs; ing++){
                Quat4d hi = hs[ing];
                for(int jng=0; jng<ing; jng++){
                    Quat4d hj = hs[jng];
                    double ang = hi.f.getAngle_unitary( hj.f );

                }
            }
        }

    }


    /**
     * Calculates the rotation matrices which characterize angular distribution of atoms of each type around the reference point (center of mass). This is done by calculating the quarupole matrices for each atom and summing them for each type. The quarupole matrices are then diagonalized to obtain the rotation matrices.
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
