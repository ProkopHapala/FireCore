
#ifndef UFF_h
#define UFF_h

#include <omp.h>
#include "fastmath.h"   // fast math operations
#include "Vec2.h"       // 2D vector
#include "Vec3.h"       // 3D vector
#include "quaternion.h" // quaternions
#include "Forces.h"     // various physical interactions
#include "SMat3.h"             // Symmetric Matrix
#include "molecular_utils.h"   // various molecular utilities
#include "NBFF.h" // Non-Bonded Force Field

#include "Draw3D.h"  // just for debug

// ========================
// ====   UFF          ====
// ========================
//
//    This is an implementation of Universal Force Field. It was made using MMFFsp3_loc as a template.
//

bool checkVec3Match( Vec3d f, Vec3d f_, const char* label, int iPrint=1 ){
    double l  = f.norm();
    double l_ = f_.norm();
    double r = l/l_;
    double c = f.dot(f_)/(l*l_);
    if(iPrint>0){
        printf( "%s : cos(f,f_)=%g |f|/|f_|=%g \n", label,  c, r );
    }else if(iPrint>1){
        printf( "%s : cos(f,f_)=%g |f|/|f_|=%g | f(%g,%g,%g) f_ref(%g,%g,%g) \n", label,  c, r,  f.x, f.y, f.z, f_.x, f_.y, f_.z );
    }
    return  ( fabs(c-1)<1e-6 ) && ( fabs(r-1)<1e-6 );
}

bool checkVec3Matches( int n, Vec3d* v, Vec3d* v_, const char* label, int iPrint=1 ){
    char strbuf[256];
    bool bMatch = true;
    for(int i=0; i<n; i++){ 
        sprintf( strbuf, "%s_%i ", label, i );
        checkVec3Match( v[i], v_[i], strbuf, iPrint ); 
    }
    return bMatch;
}

class UFF : public NBFF { public:

    // === inherited
    //Vec3d *   apos  =0;    // [natoms] // from Atoms
    //Vec3d *  fapos  =0;    // [natoms] // from NBFF
    //Quat4i*  neighs =0;    // [natoms] // from NBFF
    //Quat4i*  neighCell=0;  // [natoms] // from NBFF

    // dimensions of the system
    double Etot, Eb, Ea, Ed, Ei;                          // total, bond, angle, dihedral, inversion energies
    int natoms, nbonds, nangles, ndihedrals, ninversions; // number of bonds, angles, dihedrals, inversions
    int nDOFs;                                            // total number of degrees of freedom
    double*  DOFs __attribute__((aligned(64))) = 0;                                  // degrees of freedom
    double* fDOFs __attribute__((aligned(64))) = 0;                                  // forces                               
    Vec3d * vapos __attribute__((aligned(64))) = 0;      // [natoms] velocities of atoms

    Mat3d   invLvec;    // inverse lattice vectors

    // Auxiliary Variables
    
    Quat4d* hneigh __attribute__((aligned(64))) = 0;  // [natoms*4]     bond vectors (normalized in .xyz=f ) and their inverse length in .w=e
                        //                for each atom and each neighbor (the array in already unrolled)
    //Vec3d * fbon __attribute__((aligned(64))) = 0;  // [nbonds*2]     temporary store of forces on atoms from bonds (before the assembling step)
    Vec3d * fang __attribute__((aligned(64))) = 0;  // [nangles*3]    temporary store of forces on atoms from bonds (before the assembling step)
    Vec3d * fdih __attribute__((aligned(64))) = 0;  // [ndihedrals*4] temporary store of forces on atoms from bonds (before the assembling step)
    Vec3d * finv __attribute__((aligned(64))) = 0;  // [nimpropers*4] temporary store of forces on atoms from bonds (before the assembling step)

    // Params
    Quat4i *  neighBs   __attribute__((aligned(64))) = 0; // [natoms]      bond indices for each neighbor
    Vec2i  *  bonAtoms  __attribute__((aligned(64))) = 0; // [nbonds]      bonds atoms
    Vec2d  *  bonParams __attribute__((aligned(64))) = 0; // [nbonds]      bonds parameters
    Vec3i  *  angAtoms  __attribute__((aligned(64))) = 0; // [nangles]     angles atoms
    double5*  angParams __attribute__((aligned(64))) = 0; // [nangles]     angles parameters
    Quat4i *  dihAtoms  __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals atoms
    Vec3d  *  dihParams __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals parameters
    Quat4i *  invAtoms  __attribute__((aligned(64))) = 0; // [ninversions] inversions atoms
    Quat4d *  invParams __attribute__((aligned(64))) = 0; // [ninversions] inversions parameters

    Vec3i * dihNgs __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals neighbor index
    Vec2i * angNgs __attribute__((aligned(64))) = 0; // [nangles]     angles neighbor index
    Vec3i * invNgs __attribute__((aligned(64))) = 0; // [ninversions] inversions neighbor index

    /*

    // TBD this can be done on interactions e.g. f from bonds, f from angles...
    Vec3d * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Vec3d * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)
    */

    // =========================== Functions

    // reallocate UFF
    void realloc( int natoms_, int nbonds_, int nangles_, int ndihedrals_, int ninversions_ ){

        natoms=natoms_; 
        nbonds=nbonds_; 
        nangles=nangles_; 
        ndihedrals=ndihedrals_; 
        ninversions=ninversions_;
        nDOFs = natoms*3;
        _realloc0(  DOFs, nDOFs, (double)NAN );
        _realloc0( fDOFs, nDOFs, (double)NAN );
        apos   = (Vec3d*) DOFs ;
        fapos  = (Vec3d*)fDOFs;
        _realloc0( neighs    , natoms,   Quat4iMinusOnes );  // neighbor indices for each atom
        _realloc0( neighBs   , natoms,   Quat4iMinusOnes );  // bond indices for each neighbor
        _realloc0( neighCell , natoms,   Quat4iMinusOnes );  // cell indices for each neighbor
        _realloc0( hneigh    , natoms*4, Quat4dNAN );
        _realloc0( atypes    , natoms,   -1        );
        // ---- Aux
        //_realloc0( fbon  , nbonds*2,      Vec3dNAN );
        _realloc0( fang  , nangles*3,     Vec3dNAN );
        _realloc0( fdih  , ndihedrals*4,  Vec3dNAN );
        _realloc0( finv  , ninversions*4, Vec3dNAN );
        // ----- Params 
        _realloc0( bonAtoms  , nbonds,    Vec2iZero );
        _realloc0( bonParams , nbonds,    Vec2dNAN  );
        _realloc0( angAtoms  , nangles,   Vec3iZero );
        _realloc0( angParams , nangles,   (double5){(double)NAN,(double)NAN,(double)NAN,(double)NAN,(double)NAN}  );
        _realloc0( dihAtoms  , ndihedrals,  Quat4iZero );
        _realloc0( dihNgs    , ndihedrals,  Vec3iZero  );
        _realloc0( dihParams , ndihedrals,  Vec3dNAN   );
        _realloc0( invAtoms  , ninversions, Quat4iZero );
        _realloc0( invParams , ninversions, Quat4dNAN  );

        _realloc0( invNgs    , nangles,     Vec3iZero  );
        _realloc0( angNgs    , ninversions, Vec2iZero  );

    }

    // deallocate UFF
    void dealloc(){

        natoms=0; 
        nbonds=0; 
        nangles=0; 
        ndihedrals=0; 
        ninversions=0;
        nDOFs= 0;
        _dealloc(DOFs );
        _dealloc(fDOFs);
        apos   = 0;
        fapos  = 0;
        _dealloc(neighs);
        _dealloc(neighCell);
        _dealloc(hneigh);
        //_dealloc(fbon);
        _dealloc(fang);
        _dealloc(fdih);
        _dealloc(finv);
        _dealloc(bonAtoms);
        _dealloc(bonParams);
        _dealloc(angAtoms);
        _dealloc(angParams);
        _dealloc(dihAtoms);
        _dealloc(dihParams);
        _dealloc(invAtoms);
        _dealloc(invParams);
        
        _dealloc(angNgs);
        _dealloc(dihNgs);
        _dealloc(invNgs);

    }

    // set lattice vectors
    void setLvec(const Mat3d& lvec_){ lvec=lvec_; lvec.invert_T_to( invLvec ); }


    // ============== Evaluation

    void makeNeighBs(){
        for(int i=0; i<nbonds; i++){  
            const Vec2i& b = bonAtoms[i];
            int* ngi = neighBs[b.x].array;
            int* ngj = neighBs[b.y].array;
            for(int j=0; j<4; j++){ if(ngi[j]<0){ ngi[j]=i; break; } }
            for(int j=0; j<4; j++){ if(ngj[j]<0){ ngj[j]=i; break; } }
        }
    };

    void bakeDihedralNeighs(){
        for( int id=0; id<ndihedrals; id++){
            int i = dihAtoms[id].x;
            int j = dihAtoms[id].y;
            int k = dihAtoms[id].z;
            int l = dihAtoms[id].w;
            const int* ingsj = neighs[j].array; // neighbors
            const int* ingsk = neighs[k].array; // neighbors
            //Vec3d  r12, r32;
            //double l12, l32;
            for(int in=0; in<4; in++){
                int ing = ingsj[in];
                if(ing<0) { break; }
                if     (ing==i) { 
                    //r12 = hneigh[j*4+in].f;  
                    dihNgs[id].x = j*4+in;   // j-i
                }   
                else if(ing==k) { 
                    //r32 = hneigh[j*4+in].f;  
                    dihNgs[id].y = j*4+in;  // j-k
                } 
            }
            Vec3d r43;
            double l43;
            for(int in=0; in<4; in++){
                int ing = ingsk[in];
                if(ing<0) { break; }
                if     (ing==l) { 
                    //r43 = hneigh[k*4+in].f; 
                    dihNgs[id].z = k*4+in; // k-l
                }   
            }
        }
    }

    void bakeAngleNeighs(){
        for( int ia=0; ia<nangles; ia++){
            int i = angAtoms[ia].x;
            int j = angAtoms[ia].y;
            int k = angAtoms[ia].z;
            const int*  ings = neighs   [j].array; // neighbors
            //Vec3d  rij, rkj;
            //double lij, lkj;
            for(int in=0; in<4; in++){
                int ing = ings[in];
                if(ing<0) { break; }
                if     (ing==i) { 
                    //rij = hneigh[j*4+in].f; 
                    //lij = hneigh[j*4+in].e; 
                    angNgs[ia].x = j*4+in; // j-i

                }else if(ing==k) { 
                    //rkj = hneigh[j*4+in].f; 
                    //lkj = hneigh[j*4+in].e; 
                    angNgs[ia].y = j*4+in; // j-k
                } 
            }
        }
    }

    void bakeInversionNeighs(){
        for(int ii = 0; ii<ninversions; ii++){
            int i = invAtoms[ii].x;
            int j = invAtoms[ii].y;
            int k = invAtoms[ii].z;
            int l = invAtoms[ii].w;
            const int*    ings = neighs   [i].array; // neighbors
            //Vec3d  r21, r31, r41;
            //double l21, l31, l41;
            for(int in=0; in<3; in++){
                int ing = ings[in];
                if     (ing==j) { 
                    //r21 = hneigh[i*4+in].f; l21 = 1.0/hneigh[i*4+in].e; 
                    invNgs[ii].x = i*4+in; // i-j
                } else if(ing==k) { 
                    //r31 = hneigh[i*4+in].f; l31 = 1.0/hneigh[i*4+in].e; 
                    invNgs[ii].y = i*4+in; // i-k
                } else if(ing==l) { 
                    //r41 = hneigh[i*4+in].f; l41 = 1.0/hneigh[i*4+in].e; 
                    invNgs[ii].z = i*4+in; // i-l
                } 
            }
        }
    }

    // clear forces on all atoms and other DOFs
    //void cleanForce(){ Etot=0.0; for(int i=0; i<nDOFs; i++){ fDOFs[i]=0.0; } }
    void cleanForce(){ 
        Etot=0.0; 
        for(int i=0; i<natoms;        i++){ fapos[i]=Vec3d{0.0,0.0,0.0}; }
        //for(int i=0; i<nbonds*2;      i++){ fbon [i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<nangles*3;     i++){ fang [i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<ndihedrals*4;  i++){ fdih [i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<ninversions*4; i++){ finv [i]=Vec3d{0.0,0.0,0.0}; }
    }

    // make list of neighbors cell index (in periodic boundary conditions), by going through all periodic images
    void makeNeighCells( const Vec3i nPBC_ ){ 
        nPBC=nPBC_;
        for(int ia=0; ia<natoms; ia++){
            Quat4i ngC = Quat4i{-1,-1,-1,-1};
            for(int j=0; j<4; j++){
                const int ja = neighs[ia].array[j];
                if( ja<0 )continue;
                const Vec3d d0 = apos[ja] - apos[ia];
                int ipbc =  0;
                int imin = -1;
                double r2min = 1e+300;
                // go through all periodic images and find nearest distance
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){   
                    Vec3d d = d0 + (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz); 
                    double r2 = d.norm2();
                    if(r2<r2min){   // find nearest distance
                        r2min = r2;
                        imin  = ipbc;
                    }
                    ipbc++; 
                }}}
                ngC.array[j] = imin;
            }
            neighCell[ia]=ngC; // set neighbor cell index
        }
    }

    // make list of neighbors cell index (in periodic boundary conditions) using precomputed pbc_shifts
    void makeNeighCells( int npbc, Vec3d* pbc_shifts ){ 
        for(int ia=0; ia<natoms; ia++){
            for(int j=0; j<4; j++){
                int ja = neighs[ia].array[j];
                if( ja<0 )continue;
                const Vec3d d = apos[ja] - apos[ia];
                // ------- Brute Force method
                int imin=-1;
                float r2min = 1.e+300;
                for( int ipbc=0; ipbc<npbc; ipbc++ ){   
                    Vec3d shift = pbc_shifts[ipbc]; 
                    shift.add(d);
                    float r2 = shift.norm2();
                    if(r2<r2min){   // find nearest distance
                        r2min=r2;
                        imin=ipbc;
                    }
                }
                neighCell[ia].array[j] = imin;
            }
        }
    }

    // Full evaluation of UFF intramolecular force-field
    double eval( bool bClean=true ){

        if(bClean)cleanForce();

        Eb = evalBonds();
        Ea = evalAngles();
        Ed = evalDihedrals();
        Ei = evalInversions();
        
        //Etot = Eb + Ea + Ed + Ei;
        //assembleForces();

// DEBUG 
//Etot = Eb; 
//assembleForcesDEBUG(true,false,false,false);
//Etot = Ea; 
//assembleForcesDEBUG(false,true,false,false);
Etot = Ed; 
assembleForcesDEBUG(false,false,true,false);
//Etot = Ei; 
//assembleForcesDEBUG(false,false,false,true);
double tokcal = 60.2214076*1.602176634/4.1840;
FILE *file = fopen("out","w");
fprintf( file, "%g\n", Etot*tokcal );
fprintf( file, "%i\n", natoms );
for(int ia=0; ia<natoms; ia++){
    fprintf( file, "%i %g %g %g %g %g %g\n", ia+1, apos[ia].x, apos[ia].y, apos[ia].z, fapos[ia].x*tokcal, fapos[ia].y*tokcal, fapos[ia].z*tokcal );
}
fclose(file);
//printf("ADES SON ARIVA' FIN QUA -> UFF.h::eval()\n");exit(0);  

        return Etot;

    }

    void assembleForcesDEBUG(bool bbonds, bool bangles, bool bdihedrals, bool binversions){
        printf("assembleForcesDEBUG(bonds(%i|%i) angles(%i|%i) dihedrals(%i|%i) inversions(%i|%i) )\n", bbonds,nbonds, bangles, nangles, bdihedrals, ndihedrals, binversions, ninversions ); 
        // if(bbonds){
        //     // bonds
        //     for(int i=0; i<nbonds; i++){
        //         printf("bond[%i]\n",i);
        //         int ia1 = bonAtoms[i].x;
        //         int ia2 = bonAtoms[i].y;
        //         fapos[ia1].add( fbon[i*2] );
        //         fapos[ia2].add( fbon[i*2+1] );
        //     }
        // }
        if(bangles){
            // angles
            for(int i=0; i<nangles; i++){
                printf("angle[%i]\n",i);
                int ia1 = angAtoms[i].x;
                int ia2 = angAtoms[i].y;
                int ia3 = angAtoms[i].z;
                fapos[ia1].add( fang[i*3] );
                fapos[ia2].add( fang[i*3+1] );
                fapos[ia3].add( fang[i*3+2] );
            }
        }
        if(bdihedrals){
            // dihedrals
            for(int i=0; i<ndihedrals; i++){
                printf("dihedral[%i]\n",i);
                int ia1 = dihAtoms[i].x;
                int ia2 = dihAtoms[i].y;
                int ia3 = dihAtoms[i].z;
                int ia4 = dihAtoms[i].w;
                fapos[ia1].add( fdih[i*4] );
                fapos[ia2].add( fdih[i*4+1] );
                fapos[ia3].add( fdih[i*4+2] );
                fapos[ia4].add( fdih[i*4+3] );
            }
        }
        if(binversions){
            // inversions
            for(int i=0; i<ninversions; i++){
                printf("inversion[%i]\n",i);
                int ia1 = invAtoms[i].x;
                int ia2 = invAtoms[i].y;
                int ia3 = invAtoms[i].z;
                int ia4 = invAtoms[i].w;
                fapos[ia1].add( finv[i*4] );
                fapos[ia2].add( finv[i*4+1] );
                fapos[ia3].add( finv[i*4+2] );
                fapos[ia4].add( finv[i*4+3] );
            }
        }
    }

    void assembleForces(){
        // NOTE: this is not parallelized ( wee need somethig which loops over atoms otherwise we would need atomic add )
        // bonds
        // for(int i=0; i<nbonds; i++){
        //     int ia1 = bonAtoms[i].x;
        //     int ia2 = bonAtoms[i].y;
        //     fapos[ia1].add( fbon[i*2] );
        //     fapos[ia2].add( fbon[i*2+1] );
        // }
        // angles
        for(int i=0; i<nangles; i++){
            const Vec3i ii = angAtoms[i];
            //int ia1 = angAtoms[i].x;
            //int ia2 = angAtoms[i].y;
            //int ia3 = angAtoms[i].z;
            fapos[ii.x].add( fang[i*3] );
            fapos[ii.y].add( fang[i*3+1] );
            fapos[ii.z].add( fang[i*3+2] );
        }
        // dihedrals
        for(int i=0; i<ndihedrals; i++){
            const Quat4i ii = dihAtoms[i];
            // int ia1 = dihAtoms[i].x;
            // int ia2 = dihAtoms[i].y;
            // int ia3 = dihAtoms[i].z;
            // int ia4 = dihAtoms[i].w;
            fapos[ii.x].add( fdih[i*4] );
            fapos[ii.y].add( fdih[i*4+1] );
            fapos[ii.z].add( fdih[i*4+2] );
            fapos[ii.w].add( fdih[i*4+3] );
        }
        // inversions
        for(int i=0; i<ninversions; i++){
            // int ia1 = invAtoms[i].x;
            // int ia2 = invAtoms[i].y;
            // int ia3 = invAtoms[i].z;
            // int ia4 = invAtoms[i].w;
            const Quat4i ii = dihAtoms[i];
            fapos[ii.x].add( finv[i*4] );
            fapos[ii.y].add( finv[i*4+1] );
            fapos[ii.z].add( finv[i*4+2] );
            fapos[ii.w].add( finv[i*4+3] );
        }
    }

    void assembleAtomForce(const int ia){
        fapos[ia].add( fang[ia*3] );
    }

    void assembleAtomsForces(){
        for(int ia=0; ia<natoms; ia++){ assembleAtomForce(ia); }
    }

    inline double evalAtomBonds(const int ia){
        double E=0.0;
        const Vec3d   pa   = apos     [ia]; 
        const int*    ings = neighs   [ia].array; // neighbors
        const int*    ingC = neighCell[ia].array; // neighbors cell index
        const int*    inbs = neighBs  [ia].array; // neighbors bond index
        for(int in=0; in<4; in++){
            const int ing = ings[in];
            if(ing<0) break;
            // --- Bond vectors
            const int inn=ia*4+in;
            const Vec3d pi = apos[ing]; 
            Vec3d dp;               
            dp.set_sub( pi, pa );
            // Periodic Boundary Conditions
            if(bPBC){ 
                if(shifts){ // if we have bond shifts vectors we use them
                    int ipbc = ingC[in]; 
                    dp.add( shifts[ipbc] );
                }else{ // if we don't have bond shifts vectors we use lattice vectors
                    Vec3i g  = invLvec.nearestCell( dp );
                    Vec3d sh = lvec.a*(double)g.x + lvec.b*(double)g.y + lvec.c*(double)g.z;
                    dp.add( sh );
                }
            }
            const double l = dp.norm();
            hneigh[inn].f.set_mul( dp, 1.0/l ); 
            hneigh[inn].e = 1.0/l;
            // --- Bond Energy
            //if(ing<ia) continue; // avoid double computing - NOTE: but we do double computing (in order to make it better parallelizable)
            //int ib;
            // ToDo: this should be optimized !!!!!!!!!!!!!!!!
            // for(int i=0; i<nbonds; i++){  
            //     if( ( bonAtoms[i].x == ia && bonAtoms[i].y == ing ) || ( bonAtoms[i].y == ia && bonAtoms[i].x == ing ) ) { ib = i; break; }
            // }
            const int   ib  = inbs[in];
            const Vec2d par = bonParams[ib];
            const double dl = l-par.y;
            E += par.x*dl*dl;
            Vec3d f; f.set_mul( dp, 2.0*par.x*dl*hneigh[inn].e );
            //fbon[ib*2  ]=f; // force on atom i
            //f.mul(-1.0);
            //fbon[ib*2+1]=f;      
            fapos[ia].add(f);
            //fapos[ing].sub(f);       // NOTE: this will cause problems in parallelization !!!!!!!! ( we need to use atomic add )   
            // TBD exclude non-bonded interactions between 1-2 neighbors
        }
        return E;
    }

    double evalBonds(){
        double E=0.0;
        for(int ia=0; ia<natoms; ia++){ E+= evalAtomBonds(ia); }
        return E;
    }

    inline double evalAngle_Prokop( const int ia ){
        const Vec2i  ngs = angNgs[id];  
        const Quat4d qij = hneigh[ngs.x];  // ji
        const Quat4d qkj = hneigh[ngs.y];  // jk
        // ---- Angle ( cos, sin )
        const Vec3d  h = qij.f + qkj.f;
        const double c = 0.5*(h.norm2()-2.0);
        const Vec2d  cs{ c, sqrt(1.0-c*c+1e-14) };
        const double inv_sin = 1.0/cs.y;
        const double5 par    = angParams[ia];
        // ---- Energy & Force Fourier
        //double E = par.k * ( par.c0 + par.c1*cs.x +     par.c2*cs2.x         +     par.c3*cs3.x         );
        //double f = par.k * (          par.c1      + 2.0*par.c2*cs2.y*inv_sin + 3.0*par.c3*cs3.y*inv_sin );
        double E  = par.c0;
        double f  = par.c1;
        Vec2d csn = cs;
        csn.mul_cmplx(cs);
        E += par.c2*csn.x;
        f += par.c2*csn.y*inv_sin*2.0;
        csn.mul_cmplx(cs);
        E += par.c3*csn.x;
        f += par.c3*csn.y*inv_sin*3.0;
        E *= par.k;
        f *= par.k;
        // --- Force Vectors
        const double fi  = f*qij.w;
        const double fk  = f*qkj.w;
        const double fic = fi*c;
        const double fkc = fk*c;
        Vec3d fpi; fpi.set_lincomb(  fic,    qij.f, -fi,     qkj.f );
        Vec3d fpk; fpk.set_lincomb( -fk,     qij.f,  fkc,    qkj.f );
        Vec3d fpj; fpj.set_lincomb(  fk-fic, qij.f,  fi-fkc, qkj.f );
        const int i3=ia*3;
        fang[i3  ]=fpi;
        fang[i3+1]=fpj;
        fang[i3+2]=fpk;
        // TBD exclude non-bonded interactions between 1-3 neighbors
        { // Debug Draw
            glColor3f(1.0,0.0,1.0);
            const Vec3i ijk = angAtoms[id];
            const Vec3d pi = apos[ijk.x]; 
            const Vec3d pj = apos[ijk.y];
            const Vec3d pk = apos[ijk.z];
            Draw3D::drawArrow( pi, pi+fpi, 0.03 );
            Draw3D::drawArrow( pj, pj+fpj, 0.03 );
            Draw3D::drawArrow( pk, pk+fpk, 0.03 );
            //glColor3f(0.0,0.0,1.0); Draw3D::drawArrow( pj, pj+qij.f*(1/qij.w), 0.03 );
            //glColor3f(1.0,0.0,0.0); Draw3D::drawArrow( pj, pj+qkj.f*(1/qkj.w), 0.03 );
            //Draw3D::drawArrow( pk, pk+fpk, 0.03 );
        }
        return E;
    }

    inline double evalAngle_Paolo( const int ia ){
        int i = angAtoms[ia].x;
        int j = angAtoms[ia].y;
        int k = angAtoms[ia].z;
        const int*    ings = neighs   [j].array; // neighbors
        Vec3d  rij, rkj;
        double lij, lkj;
        for(int in=0; in<4; in++){
            int ing = ings[in];
            if(ing<0) { break; }
            if     (ing==i) { rij = hneigh[j*4+in].f; lij = hneigh[j*4+in].e; }   
            else if(ing==k) { rkj = hneigh[j*4+in].f; lkj = hneigh[j*4+in].e; } 
        }
        Vec3d h;
        h.set_add( rij, rkj );
        double cos = 0.5*(h.norm2()-2.0);
        double sin = sqrt(1.0-cos*cos+1e-14); // must be positive number !!!
        Vec2d cs{cos,sin};
        Vec2d cs2{cos,sin};
        cs2.mul_cmplx(cs);
        Vec2d cs3{cos,sin};
        cs3.mul_cmplx(cs2);
        double5 par = angParams[ia];
        double E    = par.k * ( par.c0 + par.c1*cs.x +     par.c2*cs2.x      +     par.c3*cs3.x      );
        double fact = par.k * (          par.c1      + 2.0*par.c2*cs2.y/cs.y + 3.0*par.c3*cs3.y/cs.y );
        Vec3d vec_i, vec_k;
        vec_i.set_mul(rij,cs.x);
        vec_i.set_sub(vec_i,rkj);
        vec_k.set_mul(rkj,cs.x);
        vec_k.set_sub(vec_k,rij);
        Vec3d fpi = vec_i*(fact*lij);
        Vec3d fpk = vec_k*(fact*lkj);
        Vec3d fpj = (fpi + fpk)*(-1.0);
        fang[ia*3]  =fpi;
        fang[ia*3+2]=fpk;
        fang[ia*3+1]=fpj;
        { // Debug Draw
            glColor3f(0.0,1.0,0.0);
            const Vec3i ijk = angAtoms[id];
            const Vec3d pi = apos[ijk.x]; 
            const Vec3d pj = apos[ijk.y];
            const Vec3d pk = apos[ijk.z];
            Draw3D::drawArrow( pi, pi+fpi, 0.03 );
            Draw3D::drawArrow( pj, pj+fpj, 0.03 );
            Draw3D::drawArrow( pk, pk+fpk, 0.03 );
        }
        // TBD exclude non-bonded interactions between 1-3 neighbors
        return E;
    }

    double evalAngles(){
        double E=0.0;
        for( int ia=0; ia<nangles; ia++){ 
            E+= evalAngle_Prokop(ia); 
        }
        return E;
    }

    // ====================== Dihedrals

    double evalDihedral_Prokop( const int id ){
        //double E=0.0;
        const Vec3i ngs = dihNgs[id];   // {ji, jk, kl}
        const Quat4d q12 =    hneigh[ngs.x];  // ji
        const Quat4d q32 =    hneigh[ngs.y];  // jk
        const Quat4d q43 =    hneigh[ngs.z];  // kl
        // --- normals 
        Vec3d n123; n123.set_cross( q12.f, q32.f );   //  |n123| = sin( r12, r32 ) 
        Vec3d n234; n234.set_cross( q43.f, q32.f );   //  |n234| = sin( r43, r32 )
        const double il2_123 = 1/n123.norm2();        //  il2_123 =  1/ ( sin( r12, r32 ) )^2
        const double il2_234 = 1/n234.norm2();        //  il2_234 =  1/ ( sin( r43, r32 ) )^2
        const double inv_n12 = sqrt(il2_123*il2_234); //  inv_n12 =  1/ ( sin( r12, r32 ) * sin( r43, r32 ) )
        // --- Energy
        const Vec2d cs{
             n123.dot(n234 )*inv_n12,
            -n123.dot(q43.f)*inv_n12
        };
        Vec2d csn = cs;
        const Vec3d par = dihParams[id];
        const int n = (int)par.z;
        for(int i=1; i<n; i++){ csn.mul_cmplx(cs); }
        const double E  =  par.x * ( 1.0 + par.y * csn.x );
        // --- Force on end atoms
        const double f = -par.x * par.y * par.z * csn.y; 
        Vec3d fp1; fp1.set_mul(n123,-f*il2_123*q12.w );
        Vec3d fp4; fp4.set_mul(n234, f*il2_234*q43.w );
        // --- Recoil forces on axis atoms
        const double c123   = q32.f.dot(q12.f)*(q32.w/q12.w);
        const double c432   = q32.f.dot(q43.f)*(q32.w/q43.w);
        Vec3d fp3; fp3.set_lincomb( -c123,   fp1, -c432-1., fp4 );   // from condition torq_p2=0  ( conservation of angular momentum )
        Vec3d fp2; fp2.set_lincomb( +c123-1, fp1, +c432   , fp4 );   // from condition torq_p3=0  ( conservation of angular momentum )
        const int i4=id*4;
        fdih[i4  ]=fp1;
        fdih[i4+1]=fp2;
        fdih[i4+2]=fp3;
        fdih[i4+3]=fp4;

        { // Debug Draw
            glColor3f(1.0,0.0,0.0);
            const Quat4i ijkl = dihAtoms[id];
            const Vec3d p1 = apos[ijkl.x]; 
            const Vec3d p2 = apos[ijkl.y];
            const Vec3d p3 = apos[ijkl.z];
            const Vec3d p4 = apos[ijkl.w];
            Draw3D::drawArrow( p1, p1+fp1, 0.1 );
            Draw3D::drawArrow( p2, p2+fp2, 0.1 );
            Draw3D::drawArrow( p3, p3+fp3, 0.1 );
            Draw3D::drawArrow( p4, p4+fp4, 0.1 );
        }
        return E;
    }

    double evalDihedral_Prokop_Old( const int id ){
        //double E=0.0;
        const Quat4i ijkl = dihAtoms[id];
        const Vec3d p2  = apos[ijkl.y];
        const Vec3d p3  = apos[ijkl.z];
        const Vec3d r32 = p3-p2;
        const Vec3d r12 = apos[ijkl.x]-p2; 
        const Vec3d r43 = apos[ijkl.w]-p3;
        Vec3d n123; n123.set_cross( r12, r32 );  //  |n123| = |r12| |r32| * sin( r12, r32 ) 
        Vec3d n234; n234.set_cross( r43, r32 );  //  |n234| = |r43| |r32| * sin( r43, r32 )
        
        // ===== Prokop's
        const double l32     = r32   .norm ();   // we can avoid this sqrt() if we read it from hneigh
        const double il2_123 = 1/n123.norm2();   //  il2_123 =  1/ (   |r12| |r32| * sin( r12, r32 ) )^2 
        const double il2_234 = 1/n234.norm2();   //  il2_234 =  1/ (   |r43| |r32| * sin( r43, r32 ) )^2
        const double inv_n12 = sqrt(il2_123*il2_234); // inv_n12 = 1/ ( |r12| |r32| * sin( r12, r32 ) * |r43| |r32| * sin( r43, r32 ) )
        // --- Energy
        const Vec2d cs{
             n123.dot(n234)*inv_n12     ,
            -n123.dot(r43 )*inv_n12*l32
        };
        Vec2d csn = cs;
        Vec3d par = dihParams[id];
        const int n = (int)par.z;
        for(int i=1; i<n; i++){ csn.mul_cmplx(cs); }
        double E  =  par.x * ( 1.0 + par.y * csn.x );
        // --- Force on end atoms
        double f = -par.x * par.y * par.z * csn.y; 
        f*=l32;
        Vec3d fp1; fp1.set_mul(n123,-f*il2_123 );
        Vec3d fp4; fp4.set_mul(n234, f*il2_234 );
        // --- Recoil forces on axis atoms
        double il2_32 = -1/(l32*l32);
        double c123   = r32.dot(r12)*il2_32;
        double c432   = r32.dot(r43)*il2_32;
        Vec3d fp3; fp3.set_lincomb(  c123,   fp1,  c432-1., fp4 );   // from condition torq_p2=0  ( conservation of angular momentum )
        Vec3d fp2; fp2.set_lincomb( -c123-1, fp1, -c432   , fp4 );   // from condition torq_p3=0  ( conservation of angular momentum )
        //Vec3d fp2_ = (fp1_ + fp4_ + fp3_ )*-1.0;                   // from condition ftot=0     ( conservation of linear  momentum )
        //Vec3d fp3_ = (fp1_ + fp4_ + fp2_ )*-1.0;                   // from condition ftot=0     ( conservation of linear  momentum )
        const int i4=id*4;
        fdih[i4  ]=fp1;
        fdih[i4+1]=fp2;
        fdih[i4+2]=fp3;
        fdih[i4+3]=fp4;

        { // Debug Draw
            glColor3f(1.0,0.0,1.0);
            const Quat4i ijkl = dihAtoms[id];
            const Vec3d p1 = apos[ijkl.x]; 
            const Vec3d p2 = apos[ijkl.y];
            const Vec3d p3 = apos[ijkl.z];
            const Vec3d p4 = apos[ijkl.w];
            Draw3D::drawArrow( p1, p1+fp1, 0.01 );
            Draw3D::drawArrow( p2, p2+fp2, 0.01 );
            Draw3D::drawArrow( p3, p3+fp3, 0.01 );
            Draw3D::drawArrow( p4, p4+fp4, 0.01 );
        }
        return E;
    }

    double evalDihedral_Paolo( const int id ){
        // int i = dihAtoms[id].x;
        // int j = dihAtoms[id].y;
        // int k = dihAtoms[id].z;
        // int l = dihAtoms[id].w;
        // const int*    ingsj = neighs   [j].array; // neighbors
        // const int*    ingsk = neighs   [k].array; // neighbors
        // Vec3d  r12, r32;
        // double l12, l32;
        // for(int in=0; in<4; in++){
        //     int ing = ingsj[in];
        //     if(ing<0) { break; }
        //     if     (ing==i) { r12 = hneigh[j*4+in].f; l12 = 1.0/hneigh[j*4+in].e; }   
        //     else if(ing==k) { r32 = hneigh[j*4+in].f; l32 = 1.0/hneigh[j*4+in].e; } 
        // }
        // Vec3d r43;
        // double l43;
        // for(int in=0; in<4; in++){
        //     int ing = ingsk[in];
        //     if(ing<0) { break; }
        //     if     (ing==l) { r43 = hneigh[k*4+in].f; l43 = 1.0/hneigh[k*4+in].e; }   
        // }

        //{ // we need to read the normalized vectros for hneigh because of PBC
        //printf( "evalDihedral_Paolo() id %i \n", id );
        const Vec3i ngs = dihNgs[id];   // {ji, jk, kl}
        printf( "evalDihedral_Paolo() ngs %i %i %i \n", ngs.x, ngs.y, ngs.z );
        const Vec3d  r32 =    hneigh[ngs.y].f;  // jk
        const double l32 = 1./hneigh[ngs.y].e; 
        const Vec3d  r12 =    hneigh[ngs.x].f;  // ji
        const double l12 = 1./hneigh[ngs.x].e;
        const Vec3d  r43 =    hneigh[ngs.z].f;  // kl
        const double l43 = 1./hneigh[ngs.z].e;
        //}

        //printf( "evalDihedral_Paolo() l12 %g l32 %g l43 %g \n", l12, l32, l43 );

        Vec3d r12abs; r12abs.set_mul( r12, l12 );
        Vec3d r32abs; r32abs.set_mul( r32, l32 );
        Vec3d r43abs; r43abs.set_mul( r43, l43 );
        Vec3d n123; n123.set_cross( r12abs, r32abs );
        Vec3d n234; n234.set_cross( r43abs, r32abs );
        double l123 = n123.normalize();
        double l234 = n234.normalize();
        double cos = n123.dot(n234);
        double sin = sqrt(1.0-cos*cos+1e-14); // must be positive number !!!
        Vec3d par = dihParams[id];
        int n = (int)par.z;
        Vec2d cs{cos,sin};
        Vec2d csn{cos,sin};
        for(int i=1; i<n; i++){
            csn.mul_cmplx(cs);
        }
        double E = par.x * ( 1.0 + par.y * csn.x );
        Vec3d scaled_123; scaled_123.set_mul( n123, cos );
        Vec3d scaled_234; scaled_234.set_mul( n234, cos );
        Vec3d tmp_123; tmp_123.set_sub( n123, scaled_234 );
        Vec3d tmp_234; tmp_234.set_sub( n234, scaled_123 );
        Vec3d f_12; f_12.set_cross( r32, tmp_234 );
        Vec3d f_43; f_43.set_cross( r32, tmp_123 );
        Vec3d tmp1_32; tmp1_32.set_mul( tmp_234, l12/l123 );
        Vec3d tmp2_32; tmp2_32.set_mul( tmp_123, l43/l234 );
        Vec3d vec1_32; vec1_32.set_cross( tmp1_32, r12 );
        Vec3d vec2_32; vec2_32.set_cross( tmp2_32, r43 );
        Vec3d vec_32; vec_32.set_add( vec1_32, vec2_32 );
        double fact = -par.x * par.y * par.z * csn.y / sin ;
        Vec3d f_32; f_32.set_mul(vec_32,fact);

        Vec3d fp1 = f_12 * ( fact*l32/l123 );
        Vec3d fp4 = f_43 * ( fact*l32/l234 );
        Vec3d fp2 = ( f_32 + fp1 )*-1.0; 
        Vec3d fp3 = ( f_32 - fp4 );
        const int i4=id*4;
        fdih[i4  ]=fp1;
        fdih[i4+1]=fp2;
        fdih[i4+2]=fp3;
        fdih[i4+3]=fp4;
        // fdih[id*4  ].set_mul(f_12,fact*l32/l123);
        // fdih[id*4+3].set_mul(f_43,fact*l32/l234);
        // fdih[id*4+1].set_add(fdih[id*4],f_32);
        // fdih[id*4+1].set_mul(fdih[id*4+1],-1.0);
        // fdih[id*4+2].set_sub(f_32,fdih[id*4+3]);
        { // Debug Draw
            glColor3f(0.0,0.8,0.0);
            const Quat4i ijkl = dihAtoms[id];
            const Vec3d p1 = apos[ijkl.x]; 
            const Vec3d p2 = apos[ijkl.y];
            const Vec3d p3 = apos[ijkl.z];
            const Vec3d p4 = apos[ijkl.w];
            Draw3D::drawArrow( p1, p1+fp1, 0.02 );
            Draw3D::drawArrow( p2, p2+fp2, 0.02 );
            Draw3D::drawArrow( p3, p3+fp3, 0.02 );
            Draw3D::drawArrow( p4, p4+fp4, 0.02 );

            // glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( r12abs, p2 );
            // glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos( r32abs, p2 );
            // glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos( r43abs, p3 );
            //Draw3D::drawVecInPos( r12abs, p2 )
        }
        return E;
    }

    double evalDihedrals(){
        double E=0.0;
        for( int id=0; id<ndihedrals; id++){  
            E+= evalDihedral_Prokop(id);
            //E+= evalDihedral_Paolo(id);
        }
        return E;
    }

    inline double evalInversions( const int ii ){
        // int i = invAtoms[ii].x;
        // int j = invAtoms[ii].y;
        // int k = invAtoms[ii].z;
        // int l = invAtoms[ii].w;
        //const int*    ings = neighs   [i].array; // neighbors
        // Vec3d  r21, r31, r41;
        // double l21, l31, l41;
        // for(int in=0; in<3; in++){
        //     int ing = ings[in];
        //     if     (ing==j) { r21 = hneigh[i*4+in].f; l21 = 1.0/hneigh[i*4+in].e; }   
        //     else if(ing==k) { r31 = hneigh[i*4+in].f; l31 = 1.0/hneigh[i*4+in].e; } 
        //     else if(ing==l) { r41 = hneigh[i*4+in].f; l41 = 1.0/hneigh[i*4+in].e; } 
        // }

        const Vec3i ngs  = invNgs[ii];   // {ji, ki, li}
        const Vec3d  r21 =    hneigh[ngs.y].f;  // ji
        const double l21 = 1./hneigh[ngs.y].e; 
        const Vec3d  r31 =    hneigh[ngs.x].f;  // ki
        const double l31 = 1./hneigh[ngs.x].e;
        const Vec3d  r41 =    hneigh[ngs.z].f;  // li
        const double l41 = 1./hneigh[ngs.z].e;


        Vec3d r21abs; r21abs.set_mul( r21, l21 );
        Vec3d r31abs; r31abs.set_mul( r31, l31 );
        Vec3d n123; n123.set_cross( r21abs, r31abs );
        double l123 = n123.normalize();
        double sin = -n123.dot(r41);
        double cos = sqrt(1.0-sin*sin+1e-14); // must be positive number !!!
        Quat4d par = invParams[ii];
        Vec2d cs{cos,sin};
        Vec2d cs2{cos,sin};
        cs2.mul_cmplx(cs);
        double E = par.x * ( par.y + par.z * cos + par.w * cs2.x );
        Vec3d scaled_123; scaled_123.set_mul( n123, -sin);
        Vec3d scaled_41;  scaled_41.set_mul (  r41, -sin);
        Vec3d tmp_123; tmp_123.set_sub( n123, scaled_41  );
        Vec3d tmp_41;  tmp_41.set_sub (  r41, scaled_123 );
        Vec3d f_21; f_21.set_cross( r31, tmp_41 );
        Vec3d f_31; f_31.set_cross( tmp_41, r21 );
        double fact = -par.x * ( par.z * sin + 2.0 * par.w * cs2.y ) / cos ;
        finv[ii*4+1].set_mul(f_21,   fact*l31/l123);
        finv[ii*4+2].set_mul(f_31,   fact*l21/l123);
        finv[ii*4+3].set_mul(tmp_123,fact/l41     );
        finv[ii*4  ].set_lincomb(-1.0,-1.0,-1.0,finv[ii*4+1],finv[ii*4+2],finv[ii*4+3]);
        return E;
    }

    double evalInversions(){
        double E=0.0;
        for( int ii=0; ii<ninversions; ii++){ E+=evalInversions(ii); }
        return E;
    }


    // constrain atom to fixed position
    /*
    void constrainAtom( int ia, double Kfix=1.0 ){
        printf( "constrainAtom(i=%i,K=%g)\n", ia, Kfix );
        constr[ia].f=apos[ia];
        constr[ia].w=Kfix;
    };
    */


    // ================== Print functions  

    
    void printSizes     (      ){ printf( "MMFFf4::printSizes(): nDOFs(%i) natoms(%i) nbonds(%i) nangles(%i) ndihedrals(%i) ninversions(%i) npbc(%i)\n", nDOFs,natoms,nbonds,nangles,ndihedrals,ninversions,npbc); }
    /*
    void printAtomParams(int ia){ printf("atom[%i] t%i ngs{%3i,%3i,%3i,%3i} par(%5.3f,%5.3f,%5.3f,%5.3f)  bL(%5.3f,%5.3f,%5.3f,%5.3f) bK(%6.3f,%6.3f,%6.3f,%6.3f)  Ksp(%5.3f,%5.3f,%5.3f,%5.3f) Kpp(%5.3f,%5.3f,%5.3f,%5.3f) \n", ia, atypes[ia], neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,    apars[ia].x,apars[ia].y,apars[ia].z,apars[ia].w,    bLs[ia].x,bLs[ia].y,bLs[ia].z,bLs[ia].w,   bKs[ia].x,bKs[ia].y,bKs[ia].z,bKs[ia].w,     Ksp[ia].x,Ksp[ia].y,Ksp[ia].z,Ksp[ia].w,   Kpp[ia].x,Kpp[ia].y,Kpp[ia].z,Kpp[ia].w  ); };
    void printNeighs    (int ia){ printf("atom[%i] neigh{%3i,%3i,%3i,%3i} neighCell{%3i,%3i,%3i,%3i} \n", ia, neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,   neighCell[ia].x,neighCell[ia].y,neighCell[ia].z,neighCell[ia].w ); }
    void printBKneighs  (int ia){ printf("atom[%i] bkngs{%3i,%3i,%3i,%3i} \n", ia, bkneighs[ia].x,bkneighs[ia].y,bkneighs[ia].z,bkneighs[ia].w ); }
    void printAtomParams(      ){ printf("MMFFsp3_loc::printAtomParams()\n" ); for(int i=0; i<nnode;  i++){ printAtomParams(i); }; }
    void printNeighs    (      ){ printf("MMFFsp3_loc::printNeighs()\n"     ); for(int i=0; i<natoms; i++){ printNeighs    (i);     }; }
    void printBKneighs  (      ){ printf("MMFFsp3_loc::printBKneighs()\n"   ); for(int i=0; i<natoms; i++){ printBKneighs  (i);   }; }
    void print_pipos    (      ){ printf("MMFFsp3_loc::print_pipos()\n"     ); for(int i=0; i<nnode;  i++){ printf( "pipos[%i](%g,%g,%g) r=%g\n", i, pipos[i].x,pipos[i].y,pipos[i].z, pipos[i].norm() ); } }
    void print_apos     (      ){ printf("MMFFsp3_loc::print_apos()\n"      ); for(int i=0; i<natoms; i++){ printf( "apos [%i](%g,%g,%g)\n",      i, apos[i].x ,apos[i].y ,apos[i].z                   ); } }
    void print_pbc_shifts(     ){ printf("MMFFsp3_loc::print_pbc_shifts()\n"); for(int i=0; i<npbc;   i++){ printf( "pbc_shifts[%i](%g,%g,%g)\n", i, shifts[i].x,shifts[i].y,shifts[i].z                   ); } }
    void printAngles(int ia){
        Vec3d* angles_i = angles+(ia*6);
        int* ings       = neighs[ia].array; 
        int iang=0;
        for(int i=0; i<3; i++){
            int ing = ings[i];
            if(ing<0) break;
            for(int j=i+1; j<4; j++){
                int jng  = ings[j];
                if(jng<0) break;
                if(bEachAngle){
                    Vec2d cs0_ss = angles_i[iang].xy();  
                    double ssK    = angles_i[iang].z;
                    printf( "atom[%i|%i]types{%i,%i,%i} ssK %g cs0(%g,%g) \n", ia,iang,  atypes[ing], atypes[ia], atypes[jng], ssK, cs0_ss.x, cs0_ss.y  );
                    iang++; 
                }
            }
        }
    }
    void printAngles(      ){ printf("MMFFsp3_loc::printAngles()\n"); for(int i=0; i<nnode; i++){ printAngles(i); } }
    void printTorsions(      ){ printf("MMFFsp3_loc::printTorsions()\n"); for(int i=0; i<ntors; i++){ printf( "torsion[%i]{%i,%i,%i,%i} {%i,%i,%i,%i}\n", i, tors2atom[i].x,tors2atom[i].y,tors2atom[i].z,tors2atom[i].w,   torsParams[i].x,torsParams[i].y,torsParams[i].z,torsParams[i].w  ); } }
    void printAtomsConstrains( bool bWithOff=false ){ printf("MMFFsp3_loc::printAtomsConstrains()\n"); for(int i=0; i<natoms; i++){ if(bWithOff || (constr[i].w>0.0f) )printf( "consrt[%i](%g,%g,%g|K=%g)\n", i, constr[i].x,constr[i].y,constr[i].z,constr[i].w ); } }
    */

    // check if there are NaNs in the arrays (and exit with error-message if there are)
    /*
    bool checkNans( bool bExit=true, bool bNg=true, bool bPi=true, bool bA=true ){
        bool ret = false;
        if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*) apos,   "apos"  );
        if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*)fapos,  "fapos"  );
        if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*) pipos,  "pipos" );
        if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*)fpipos, "fpipos" );
        if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneigh,  "fneigh"   );
        if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneighpi,"fneighpi" );
        if(bExit&&ret){ printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); 
            printSizes();
            printAtomParams();
            printNeighs();
            print_pbc_shifts();
            printDEBUG(  false, false );
            eval_atoms(true,true);
            printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); 
            exit(0); 
        }
        return ret;
    }
    */

    // rotate selected atoms including their caps
    /*
    void rotateNodes(int n, int* sel, Vec3d p0, Vec3d ax, double phi ){
        ax.normalize();
        double ca=cos(phi);
        double sa=sin(phi);
        for(int i=0;i<n; i++){
            int ia = sel[i];
            apos [ia].rotate_csa( ca, sa, ax, p0 );
            if(ia>=nnode)continue;
            pipos[ia].rotate_csa( ca, sa, ax     );
            int* ngs=neighs[ia].array; 
            for(int j=0;j<4;j++){
                int ja = ngs[j];
                if(ja>=0){ if(ja>nnode) apos[ ja  ].rotate_csa( ca, sa, ax, p0 ); }
            }
        }
    }
    */

   // ============== Move atoms in order to minimize energy
   
    // Loop which iteratively evaluate MMFFsp3 intramolecular force-field and move atoms to minimize energy
    /*
    int run( int niter, double dt, double Fconv, double Flim ){
        double F2conv = Fconv*Fconv;
        double E,F2;
        int    itr;
        for(itr=0; itr<niter; itr++){
            E=0;
            // ------ eval MMFF
            for(int ia=0; ia<natoms; ia++){ 
                if(ia<nnode)E += eval_atom(ia);
                E += evalLJQs_ng4_PBC_atom( ia ); 
            }
            // ---- assemble (we need to wait when all atoms are evaluated)
            for(int ia=0; ia<natoms; ia++){
                assemble_atom( ia );
            }
            // ------ move
            F2=0;
            for(int i=0; i<nvecs; i++){
                F2 += move_atom_MD( i, dt, Flim, 0.99 );
            }
            if(F2<F2conv)break;
        }
        return itr;
    }
    */

    // Loop which iteratively evaluate MMFFsp3 intramolecular force-field and move atoms to minimize energy, wih OpenMP parallelization
    /*
    int run_omp( int niter, double dt, double Fconv, double Flim ){
        double F2conv = Fconv*Fconv;
        double E=0,F2=0;
        int    itr=0;
        #pragma omp parallel shared(E,F2) private(itr)
        for(itr=0; itr<niter; itr++){
            // This {} should be done just by one of the processors
            #pragma omp single
            {E=0;F2=0;}
            // ------ eval MMFF
            #pragma omp for reduction(+:E)
            for(int ia=0; ia<natoms; ia++){ 
                if(verbosity>3)printf( "atom[%i]@cpu[%i/%i]\n", ia, omp_get_thread_num(), omp_get_num_threads()  );
                if(ia<nnode)E += eval_atom(ia);
                E += evalLJQs_ng4_PBC_atom_omp( ia ); 
            }
            // ---- assemble (we need to wait when all atoms are evaluated)
            #pragma omp for
            for(int ia=0; ia<natoms; ia++){
                assemble_atom( ia );
            }
            // ------ move
            #pragma omp for reduction(+:F2)
            for(int i=0; i<nvecs; i++){
                F2 += move_atom_MD( i, dt, Flim, 0.99 );
            }
            for(int i=natoms; i<nvecs; i++){
                F2 += move_atom_MD( i, dt, Flim, 0.99 );
            }
            #pragma omp single
            if(verbosity>2){printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, F2, omp_get_num_threads() );}
        }
        return itr;
    }
    */

    // update atom positions using gradient descent
    /*
    inline double move_atom_GD(int i, float dt, double Flim){
        Vec3d  f   = fapos[i];
        Vec3d  p = apos [i];
        double fr2 = f.norm2();
        if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
        const bool bPi = i>=natoms;
        if(bPi){ f.add_mul( p, -p.dot(f) ); }
        p.add_mul( f, dt );
        if(bPi)p.normalize();
        apos [i] = p;
        return fr2;
    }
    */

    // update atom positions using molecular dynamics (damped leap-frog)
    /*
    inline double move_atom_MD( int i, const float dt, const double Flim, const double cdamp=0.9 ){
        Vec3d  f = fapos[i];
        Vec3d  v = vapos[i];
        Vec3d  p = apos [i];
        const double fr2 = f.norm2();
        const bool bPi = i>=natoms;
        if(bPi)f.add_mul( p, -p.dot(f) );           //b|=ckeckNaN_d( 1,3, (double*)&f, "f.2" );
        v.mul    ( cdamp );
        v.add_mul( f, dt );                         //b|=ckeckNaN_d( 1,3, (double*)&v, "v.2" );
        if(bPi)v.add_mul( p, -p.dot(v) );           //b|=ckeckNaN_d( 1,3, (double*)&v, "v.3" );
        p.add_mul( v, dt );
        if(bPi)p.normalize();
        apos [i] = p;
        vapos[i] = v;
        return fr2;
    }
    */

    // update atom positions using FIRE (Fast Inertial Relaxation Engine)
    /*
    inline double move_atom_FIRE( int i, float dt, double Flim, double cv, double cf ){
        Vec3d  f   = fapos[i];
        Vec3d  v   = vapos[i];
        Vec3d  p   = apos [i];
        double fr2 = f.norm2();
        const bool bPi = i>=natoms;
        if(bPi)f.add_mul( p, -p.dot(f) ); 
        v.mul(             cv );
        v.add_mul( f, dt + cf );
        if(bPi) v.add_mul( p, -p.dot(v) );
        p.add_mul( v, dt );
        if(bPi)  p.normalize();
        apos [i] = p;
        vapos[i] = v;
        return fr2;
    }
    */

    // update atom positions using modified FIRE algorithm
    /*
    inline double move_atom_kvaziFIRE( int i, float dt, double Flim ){
        Vec3d  f   = fapos[i];
        Vec3d        v   = vapos[i];
        Vec3d        p   = apos [i];
        double fr2 = f.norm2();
        if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); }
        const bool bPi = i>=natoms;
        if(bPi)f.add_mul( p, -p.dot(f) ); 
        double vv  = v.norm2();
        double ff  = f.norm2();
        double vf  = v.dot(f);
        double c   = vf/sqrt( vv*ff + 1.e-16    );
        double v_f =    sqrt( vv/( ff + 1.e-8)  );
        double cv;
        double cf = cos_damp_lin( c, cv, 0.01, -0.7,0.0 );
        v.mul(                 cv );
        v.add_mul( f, dt + v_f*cf );
        if(bPi) v.add_mul( p, -p.dot(v) );
        p.add_mul( v, dt );
        if(bPi)  p.normalize();
        apos [i] = p;
        vapos[i] = v;
        return fr2;
    }
    */

    // update atom positions using gradient descent
    /*
    double move_GD(float dt, double Flim=100.0 ){
        double F2sum=0;
        for(int i=0; i<nvecs; i++){
            F2sum += move_atom_GD(i, dt, Flim);
        }
        return F2sum;
    }
    */

    // find optimal time-step for dy FIRE optimization algorithm 
    /*
    double optimalTimeStep(double m=1.0){
        double Kmax = 1.0;
        for(int i=0; i<nnode; i++){ 
            Kmax=fmax(Kmax, bKs[i].x ); 
            Kmax=fmax(Kmax, bKs[i].y ); 
            Kmax=fmax(Kmax, bKs[i].z ); 
            Kmax=fmax(Kmax, bKs[i].w ); 
        }
        return M_PI*2.0*sqrt(m/Kmax)/10.0;  // dt=T/10;   T = 2*pi/omega = 2*pi*sqrt(m/k)
    }
    */

};

#endif

