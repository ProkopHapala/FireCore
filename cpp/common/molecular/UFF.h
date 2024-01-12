
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

// ========================
// ====   UFF          ====
// ========================
//
//    This is an implementation of Universal Force Field. It was made using MMFFsp3_loc as a template.
//


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
    double*  DOFs   = 0;                                  // degrees of freedom
    double* fDOFs   = 0;                                  // forces
    Mat3d   invLvec;                                      // inverse lattice vectors
    Vec3d * vapos  __attribute__((aligned(64))) = 0;      // [natoms] velocities of atoms

    // Auxiliary Variables
    Quat4d* hneigh= 0;  // [natoms*4]     bond vectors (normalized in .xyz=f ) and their inverse length in .w=e
                        //                for each atom and each neighbor (the array in already unrolled)
    Vec3d * fbon  = 0;  // [nbonds*2]     temporary store of forces on atoms from bonds (before the assembling step)
    Vec3d * fang  = 0;  // [nangles*3]    temporary store of forces on atoms from bonds (before the assembling step)
    Vec3d * fdih  = 0;  // [ndihedrals*4] temporary store of forces on atoms from bonds (before the assembling step)
    Vec3d * finv  = 0;  // [nimpropers*4] temporary store of forces on atoms from bonds (before the assembling step)

    // Params
    Vec2i  *  bonAtoms  __attribute__((aligned(64))) = 0; // [nbonds]      bonds atoms
    Vec2d  *  bonParams __attribute__((aligned(64))) = 0; // [nbonds]      bonds parameters
    Vec3i  *  angAtoms  __attribute__((aligned(64))) = 0; // [nangles]     angles atoms
    double5*  angParams __attribute__((aligned(64))) = 0; // [nangles]     angles parameters
    Quat4i *  dihAtoms  __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals atoms
    Vec3d  *  dihParams __attribute__((aligned(64))) = 0; // [ndihedrals]  dihedrals parameters
    Quat4i *  invAtoms  __attribute__((aligned(64))) = 0; // [ninversions] inversions atoms
    Quat4d *  invParams __attribute__((aligned(64))) = 0; // [ninversions] inversions parameters

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
        _realloc0( neighs    , natoms, Quat4iMinusOnes );
        _realloc0( neighCell , natoms, Quat4iMinusOnes );
        _realloc0( hneigh    , natoms*4, Quat4dNAN );
        _realloc0( atypes    , natoms, -1 );
        // ---- Aux
        _realloc0( fbon  , nbonds*2, Vec3dNAN );
        _realloc0( fang  , nangles*3, Vec3dNAN );
        _realloc0( fdih  , ndihedrals*4, Vec3dNAN );
        _realloc0( finv  , ninversions*4, Vec3dNAN );
        // ----- Params 
        _realloc0( bonAtoms  , nbonds, Vec2iZero );
        _realloc0( bonParams , nbonds, Vec2dNAN  );
        _realloc0( angAtoms  , nangles, Vec3iZero );
        _realloc0( angParams , nangles, (double5){(double)NAN,(double)NAN,(double)NAN,(double)NAN,(double)NAN}  );
        _realloc0( dihAtoms  , ndihedrals, Quat4iZero );
        _realloc0( dihParams , ndihedrals, Vec3dNAN  );
        _realloc0( invAtoms  , ninversions, Quat4iZero );
        _realloc0( invParams , ninversions, Quat4dNAN  );

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
        _dealloc(fbon);
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

    }

    // set lattice vectors
    void setLvec(const Mat3d& lvec_){ lvec=lvec_; lvec.invert_T_to( invLvec ); }


    // ============== Evaluation

    // clear forces on all atoms and other DOFs
    //void cleanForce(){ Etot=0.0; for(int i=0; i<nDOFs; i++){ fDOFs[i]=0.0; } }
    void cleanForce(){ 
        Etot=0.0; 
        for(int i=0; i<natoms; i++){       fapos[i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<nbonds*2; i++){      fbon[i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<nangles*3; i++){     fang[i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<ndihedrals*4; i++){  fdih[i]=Vec3d{0.0,0.0,0.0}; }
        for(int i=0; i<ninversions*4; i++){ finv[i]=Vec3d{0.0,0.0,0.0}; }
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
        //Ei = evalInversions();
        
        //Etot = Eb + Ea + Ed + Ei;
        //assembleForces();

// DEBUG 
//Etot = Eb; 
//assembleForcesDEBUG(true,false,false,false);
//Etot = Ea; 
//assembleForcesDEBUG(false,true,false,false);
//Etot = Ed; 
//assembleForcesDEBUG(false,false,true,false);
Etot = Ei; 
assembleForcesDEBUG(false,false,false,true);
double tokcal = 60.2214076*1.602176634/4.1840;
FILE *file = fopen("out","w");
fprintf( file, "%g\n", Etot*tokcal );
fprintf( file, "%i\n", natoms );
for(int ia=0; ia<natoms; ia++){
    fprintf( file, "%i %g %g %g %g %g %g\n", ia+1, apos[ia].x, apos[ia].y, apos[ia].z, fapos[ia].x*tokcal, fapos[ia].y*tokcal, fapos[ia].z*tokcal );
//printf("%i %g %g %g %g %g %g\n", ia+1, apos[ia].x, apos[ia].y, apos[ia].z, fapos[ia].x*tokcal, fapos[ia].y*tokcal, fapos[ia].z*tokcal );
}
fclose(file);
printf("ADES SON ARIVA' FIN QUA -> UFF.h::eval()\n");exit(0);  

        return Etot;

    }

    void assembleForcesDEBUG(bool bbonds, bool bangles, bool bdihedrals, bool binversions){
        if(bbonds){
            // bonds
            for(int i=0; i<nbonds; i++){
                int ia1 = bonAtoms[i].x;
                int ia2 = bonAtoms[i].y;
                fapos[ia1].add( fbon[i*2] );
                fapos[ia2].add( fbon[i*2+1] );
            }
        }
        if(bangles){
            // angles
            for(int i=0; i<nangles; i++){
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

        // bonds
        for(int i=0; i<nbonds; i++){
            int ia1 = bonAtoms[i].x;
            int ia2 = bonAtoms[i].y;
            fapos[ia1].add( fbon[i*2] );
            fapos[ia2].add( fbon[i*2+1] );
        }
        // angles
        for(int i=0; i<nangles; i++){
            int ia1 = angAtoms[i].x;
            int ia2 = angAtoms[i].y;
            int ia3 = angAtoms[i].z;
            fapos[ia1].add( fang[i*3] );
            fapos[ia2].add( fang[i*3+1] );
            fapos[ia3].add( fang[i*3+2] );
        }
        // dihedrals
        for(int i=0; i<ndihedrals; i++){
            int ia1 = dihAtoms[i].x;
            int ia2 = dihAtoms[i].y;
            int ia3 = dihAtoms[i].z;
            int ia4 = dihAtoms[i].w;
            fapos[ia1].add( fdih[i*4] );
            fapos[ia2].add( fdih[i*4+1] );
            fapos[ia3].add( fdih[i*4+2] );
            fapos[ia4].add( fdih[i*4+3] );
        }
        // inversions
        for(int i=0; i<ninversions; i++){
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

    double evalBonds(){

        double E=0.0;
        for(int ia=0; ia<natoms; ia++){
            const Vec3d   pa   = apos     [ia]; 
            const int*    ings = neighs   [ia].array; // neighbors
            const int*    ingC = neighCell[ia].array; // neighbors cell index
            for(int in=0; in<4; in++){
                int ing = ings[in];
                if(ing<0) break;
                // --- Bond vectors
                int inn=ia*4+in;
                Vec3d  pi = apos[ing]; 
                Vec3d  dp;               
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
                double l = dp.norm();
                hneigh[inn].f.set_mul( dp, 1.0/l ); 
                hneigh[inn].e    = 1.0/l;

                // --- Bond Energy
                if(ing<ia) continue; // avoid double computing
                int ib;
                for(int i=0; i<nbonds; i++){
                    if( ( bonAtoms[i].x == ia && bonAtoms[i].y == ing ) || 
                    ( bonAtoms[i].y == ia && bonAtoms[i].x == ing ) ) { ib = i; break; }
                }

                Vec2d par= bonParams[ib];
                double dl= l-par.y;
                E += par.x*dl*dl;
                fbon[ib*2].set_mul( dp, 2.0*par.x*dl*hneigh[inn].e ); // force on atom i
                fbon[ib*2+1].set_mul(fbon[ib*2],-1.0);                         // force on atom j
            }
        }
        return E;

    }

    double evalAngles(){

        double E=0.0;
        for( int ia=0; ia<nangles; ia++){
            int i = angAtoms[ia].x;
            int j = angAtoms[ia].y;
            int k = angAtoms[ia].z;
            const int*    ings = neighs   [j].array; // neighbors
            const int*    ingC = neighCell[j].array; // neighbors cell index
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
            E += par.k * ( par.c0 + par.c1*cs.x + par.c2*cs2.x + par.c3*cs3.x );
            double fact = par.k * ( par.c1 + 2.0*par.c2*cs2.y/cs.y + 3.0*par.c3*cs3.y/cs.y );
            Vec3d vec_i, vec_k;
            vec_i.set_mul(rij,cs.x);
            vec_i.set_sub(vec_i,rkj);
            vec_k.set_mul(rkj,cs.x);
            vec_k.set_sub(vec_k,rij);
            fang[ia*3].  set_mul(vec_i,fact*lij);
            fang[ia*3+2].set_mul(vec_k,fact*lkj);
            fang[ia*3+1].set_add(fang[ia*3],fang[ia*3+2]);
            fang[ia*3+1].set_mul(fang[ia*3+1],-1.0);
        }

        return E;
    }

    double evalDihedrals(){

        double E=0.0;
        for( int id=0; id<ndihedrals; id++){

            int i = dihAtoms[id].x;
            int j = dihAtoms[id].y;
            int k = dihAtoms[id].z;
            int l = dihAtoms[id].w;
            const int*    ingsj = neighs   [j].array; // neighbors
            const int*    ingCj = neighCell[j].array; // neighbors cell index
            const int*    ingsk = neighs   [k].array; // neighbors
            const int*    ingCk = neighCell[k].array; // neighbors cell index
            Vec3d  r12, r32;
            double l12, l32;
            for(int in=0; in<4; in++){
                int ing = ingsj[in];
                if(ing<0) { break; }
                if     (ing==i) { r12 = hneigh[j*4+in].f; l12 = 1.0/hneigh[j*4+in].e; }   
                else if(ing==k) { r32 = hneigh[j*4+in].f; l32 = 1.0/hneigh[j*4+in].e; } 
            }
            Vec3d r43;
            double l43;
            for(int in=0; in<4; in++){
                int ing = ingsk[in];
                if(ing<0) { break; }
                if     (ing==l) { r43 = hneigh[k*4+in].f; l43 = 1.0/hneigh[k*4+in].e; }   
            }
            Vec3d r12abs; r12abs.set_mul( r12, l12 );
            Vec3d r32abs; r32abs.set_mul( r32, l32 );
            Vec3d r43abs; r43abs.set_mul( r43, l43 );
            Vec3d n123; n123.set_cross( r12abs, r32abs );
            Vec3d n234; n234.set_cross( r43abs, r32abs );
            double l123 = n123.normalize();
            double l234 = n234.normalize();
            double cos = n123.dot(n234);
            //double cos = abs(n123.dot(n234))/(l123*l234);
            double sin = sqrt(1.0-cos*cos+1e-14); // must be positive number !!!
            Vec3d par = dihParams[id];
            int n = (int)par.z;
            Vec2d cs{cos,sin};
            Vec2d csn{cos,sin};
            for(int i=1; i<n; i++){
                csn.mul_cmplx(cs);
            }
            E += par.x * ( 1.0 + par.y * csn.x );
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
            fdih[id*4  ].set_mul(f_12,fact*l32/l123);
            fdih[id*4+3].set_mul(f_43,fact*l32/l234);
            fdih[id*4+1].set_add(fdih[id*4],f_32);
            fdih[id*4+1].set_mul(fdih[id*4+1],-1.0);
            fdih[id*4+2].set_sub(f_32,fdih[id*4+3]);
// OLD
/*
            r12.set_mul( r12, l12 );
            r32.set_mul( r32, l32 );
            r43.set_mul( r43, l43 );
            Vec3d n123; n123.set_cross( r12, r32 );
            Vec3d n234; n234.set_cross( r43, r32 );
            double l123 = n123.norm();
            double l234 = n234.norm();
            double cos = n123.dot(n234)/(l123*l234);
            //double cos = abs(n123.dot(n234))/(l123*l234);
            double sin = sqrt(1.0-cos*cos+1e-14); // must be positive number !!!
            Vec3d par = dihParams[id];
            int n = (int)par.z;
            Vec2d cs{cos,sin};
            Vec2d csn{cos,sin};
            for(int i=1; i<n; i++){
                csn.mul_cmplx(cs);
            }
            E += par.x * ( 1.0 + par.y * csn.x );
            Vec3d tmp1_123; tmp1_123.set_mul( n123, 1.0/(l123*l234) );
            Vec3d tmp2_123; tmp2_123.set_mul( n123, cos/(l123*l123) );
            Vec3d tmp1_234; tmp1_234.set_mul( n234, 1.0/(l123*l234) );
            Vec3d tmp2_234; tmp2_234.set_mul( n234, cos/(l234*l234) );
            Vec3d tmp_12; tmp_12.set_sub( tmp1_234, tmp2_123 );
            Vec3d vec_12; vec_12.set_cross( r32, tmp_12 );
            Vec3d tmp_43; tmp_43.set_sub( tmp1_123, tmp2_234 );
            Vec3d vec_43; vec_43.set_cross( r32, tmp_43 );
            Vec3d vec1_32; vec1_32.set_cross( tmp_12, r12 );
            Vec3d vec2_32; vec2_32.set_cross( tmp_43, r43 );
            Vec3d vec_32; vec_32.set_add( vec1_32, vec2_32 );
            double fact = -par.x * par.y * par.z * csn.y / sin ;
            Vec3d f32; f32.set_mul(vec_32,fact);
            fdih[id*4  ].set_mul(vec_12,fact);
            fdih[id*4+3].set_mul(vec_43,fact);
            fdih[id*4+1].set_add(fdih[id*4],f32);
            fdih[id*4+1].set_mul(fdih[id*4+1],-1.0);
            fdih[id*4+2].set_sub(f32,fdih[id*4+3]);
*/            
//printf("dihedral=%i energy=%g i=%i j=%i k=%i l=%i fdih_i=(%g,%g,%g) fdih_j=(%g,%g,%g) fdih_k=(%g,%g,%g) fdih_l=(%g,%g,%g)\n",id+1,i+1,j+1,k+1,l+1,par.x*(1.0+par.y*csn.x)*23.060547831,
//fdih[id*4].x*23.060547831,fdih[id*4].y*23.060547831,fdih[id*4].z*23.060547831,fdih[id*4+1].x*23.060547831,fdih[id*4+1].y*23.060547831,fdih[id*4+1].z*23.060547831,fdih[id*4+2].x*23.060547831,fdih[id*4+2].y*23.060547831,fdih[id*4+2].z*23.060547831,fdih[id*4+3].x*23.060547831,fdih[id*4+3].y*23.060547831,fdih[id*4+3].z*23.060547831);
        }

        return E;
    }


/*        // ======= Angle Step : we compute forces due to angles between bonds, and also due to angles between bonds and pi-vectors
        const double R2damp=Rdamp*Rdamp;
        if(doAngles){
        double  ssK,ssC0;
        Vec2d   cs0_ss;
        Vec3d*  angles_i;
        if(bEachAngle){ // if we compute angle energy for each angle separately, we need to store angle parameters for each angle
            angles_i = angles+(ia*6);
        }else{          // otherwise we use common parameters for all angles
            ssK    = apar.z;
            cs0_ss = Vec2d{apar.x,apar.y};
            ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
        }
        int iang=0;
        for(int i=0; i<3; i++){
            int ing = ings[i];
            if(ing<0) break;
            const Quat4d& hi = hs[i];
            for(int j=i+1; j<4; j++){
                int jng  = ings[j];
                if(jng<0) break;
                const Quat4d& hj = hs[j];    
                if(bEachAngle){ //
                    // 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
                    //  0    1    2    3    4    5 
                    cs0_ss = angles_i[iang].xy();  
                    ssK    = angles_i[iang].z;
                    iang++; 
                };
                if( bAngleCosHalf ){
                    E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
                }else{ 
                    E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
                }
                fa    .sub( f1+f2  );  // apply force on the central atom
                // ----- Error is HERE
                if(bSubtractAngleNonBond){ // subtract non-bonded interactions between atoms which have common neighbor
                    Vec3d fij=Vec3dZero;
                    Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);  // combine van der Waals parameters for the pair of atoms
                    Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f ); // compute vector between neighbors i,j
                    E -= getLJQH( dp, fij, REQij, R2damp ); // subtract non-bonded interactions 
                    f1.sub(fij); // 
                    f2.add(fij);
                }
                // apply forces on neighbors
                fbs[i].add( f1     ); 
                fbs[j].add( f2     );
            }
        }
        }

        fapos [ia]=fa; 
        fpipos[ia]=fpi;
        return E;
    }
*/
/*
        for(int ia=0; ia<natoms; ia++){
            const Vec3d   pa   = apos     [ia]; 
            const int*    ings = neighs   [ia].array; // neighbors
            const int*    ingC = neighCell[ia].array; // neighbors cell index
            for(int in=0; in<4; in++){
                int ing = ings[in];
                if(ing<0) break;
                // --- Bond vectors
                int inn=ia*4+in;
                Vec3d  pi = apos[ing]; 
                Vec3d  dp;               
                dp.set_sub( pi, pa );
                if(bPBC){ // Periodic Boundary Conditions
                    if(shifts){ // if we have bond shifts vectors we use them
                        int ipbc = ingC[in]; 
                        dp.add( shifts[ipbc] );
                    }else{ // if we don't have bond shifts vectors we use lattice vectors
                        Vec3i g  = invLvec.nearestCell( dp );
                        Vec3d sh = lvec.a*(double)g.x + lvec.b*(double)g.y + lvec.c*(double)g.z;
                        dp.add( sh );
                    }
                }
                double l = dp.norm();
                hneigh[inn].f.set_mul( dp, 1.0/l ); 
                hneigh[inn].e    = 1.0/l;
//printf("atom=%i neighbor=%i hb=(%g,%g,%g,%g)\n",ia,ing,hneigh[inn].x,hneigh[inn].y,hneigh[inn].z,hneigh[inn].w);

                // --- Bond Energy
                if(ing<ia) continue; // avoid double computing
                int ib;
                for(int i=0; i<nbonds; i++){
                    if( ( bonAtoms[i].x == ia && bonAtoms[i].y == ing ) || 
                    ( bonAtoms[i].y == ia && bonAtoms[i].x == ing ) ) { ib = i; break; }
                }
                Vec2i iva= bonAtoms [ib];
                Vec2d par= bonParams[ib];
                double dl= l-par.y;
                E += par.x*dl*dl;
                fbon[ib].set_mul( dp, 2.0*par.x*dl*hneigh[inn].e );
//printf("bond=%i energy=%g fbon=(%g,%g,%g)\n",ib,par.y*dl*dl,fbon[ib].x,fbon[ib].y,fbon[ib].z);
            }
        }
*/

    // evaluate energy and forces for single atom (ia) depending on its neighbors
    /*
    double eval_atom(const int ia){
        double E=0;
        const Vec3d pa  = apos [ia]; 
        const Vec3d hpi = pipos[ia]; 
        Vec3d fa   = Vec3dZero;
        Vec3d fpi  = Vec3dZero; 
        //--- array aliases
        const int*    ings = neighs   [ia].array; // neighbors
        const int*    ingC = neighCell[ia].array; // neighbors cell index
        const double* bK   = bKs      [ia].array; // bond stiffness
        const double* bL   = bLs      [ia].array; // bond length
        const double* Kspi = Ksp      [ia].array; // pi-sigma stiffness
        const double* Kppi = Kpp      [ia].array; // pi-pi stiffness
        Vec3d* fbs  = fneigh   +ia*4;             // forces on bonds
        Vec3d* fps  = fneighpi +ia*4;             // forces on pi vectors
        const Quat4d& apar  = apars[ia]; // [c0, Kss, Ksp, c0_e] c0 is cos of equilibrium angle, Kss is bond stiffness, Ksp is pi-sigma stiffness, c0_e is cos of equilibrium angle for pi-electron interaction
        const double  piC0 = apar.w;     // cos of equilibrium angle for pi-electron interaction
        //--- Aux Variables 
        Quat4d  hs[4]; // bond vectors (normalized in .xyz ) and their inverse length in .w
        Vec3d   f1,f2; // working forces
        for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break in the loop
        // --------- Bonds Step
        for(int i=0; i<4; i++){ // loop over bonds
            int ing = ings[i];
            if(ing<0) break;
            Vec3d  pi = apos[ing];
            Quat4d h; 
            h.f.set_sub( pi, pa );
            // Periodic Boundary Conditions
            if(bPBC){    
                if(shifts){ // if we have bond shifts vectors we use them
                    int ipbc = ingC[i]; 
                    h.f.add( shifts[ipbc] );
                }else{  // if we don't have bond shifts vectors we use lattice vectors
                    Vec3i g  = invLvec.nearestCell( h.f );
                    Vec3d sh = lvec.a*g.x + lvec.b*g.y + lvec.c*g.z;
                    h.f.add( sh );
                }
            }
            // initial bond vectors

            // TBD actual length between ia and neighbor (in the right cell)
            double l = h.f.normalize(); 
            h.e    = 1/l;

            // TBD this is an array with normalized directions and magnitude for all neighbors
            // TBD we will compute and store for all bonds
            hs [i] = h;
            if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
                if(doBonds){
                    // bond length force
                    E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1); 
                }

                double kpp = Kppi[i];
                if( (doPiPiI) && (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                    // pi-pi interaction (make them parallel)
                    E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                }
            } 
            // pi-sigma 
            double ksp = Kspi[i];
            if( doPiSigma && (ksp>1e-6) ){  
                // pi-sigma interaction (make them orthogonal)
                E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)
            }
        }



        // ======= Angle Step : we compute forces due to angles between bonds, and also due to angles between bonds and pi-vectors
        const double R2damp=Rdamp*Rdamp;
        if(doAngles){
        double  ssK,ssC0;
        Vec2d   cs0_ss;
        Vec3d*  angles_i;
        if(bEachAngle){ // if we compute angle energy for each angle separately, we need to store angle parameters for each angle
            angles_i = angles+(ia*6);
        }else{          // otherwise we use common parameters for all angles
            ssK    = apar.z;
            cs0_ss = Vec2d{apar.x,apar.y};
            ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
        }
        int iang=0;
        for(int i=0; i<3; i++){
            int ing = ings[i];
            if(ing<0) break;
            const Quat4d& hi = hs[i];
            for(int j=i+1; j<4; j++){
                int jng  = ings[j];
                if(jng<0) break;
                const Quat4d& hj = hs[j];    
                if(bEachAngle){ //
                    // 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
                    //  0    1    2    3    4    5 
                    cs0_ss = angles_i[iang].xy();  
                    ssK    = angles_i[iang].z;
                    iang++; 
                };
                if( bAngleCosHalf ){
                    E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
                }else{ 
                    E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
                }
                fa    .sub( f1+f2  );  // apply force on the central atom
                // ----- Error is HERE
                if(bSubtractAngleNonBond){ // subtract non-bonded interactions between atoms which have common neighbor
                    Vec3d fij=Vec3dZero;
                    Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);  // combine van der Waals parameters for the pair of atoms
                    Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f ); // compute vector between neighbors i,j
                    E -= getLJQH( dp, fij, REQij, R2damp ); // subtract non-bonded interactions 
                    f1.sub(fij); // 
                    f2.add(fij);
                }
                // apply forces on neighbors
                fbs[i].add( f1     ); 
                fbs[j].add( f2     );
            }
        }
        }

        fapos [ia]=fa; 
        fpipos[ia]=fpi;
        return E;
    }
    */

   /*
    double eval_atoms( bool bDebug=false, bool bPrint=false ){
        //FILE *file = fopen("out","w");
        //Etot,
        //Eb=0;Ea=0;Eps=0;EppT=0;EppI=0;
        double E=0;
        if  (bDebug){ for(int ia=0; ia<nnode; ia++){  E+=eval_atom_debug(ia,bPrint); }}
        else        { for(int ia=0; ia<nnode; ia++){  E+=eval_atom(ia);              }}
        //printf( "MYOUTPUT Ebond= %g Eangle= %g Edihed= %g Eimpr = %g Etot=%g\n", Eb, Ea, EppI, Eps, E );
        //fprintf( file, "Ebond= %g Eangle= %g Edihed= %g Eimpr= %g Etot=%g \n", Eb, Ea, EppI, Eps, E );
        //fclose(file);
        return E;
    }
    */

    // compute torsion energy and forces for a given torsion angle
    /*
    double eval_torsion(int it){ 
        printf( "MMFFsp3_loc::eval_torsion(%i)\n", it );
        Quat4i ias = tors2atom [it];
        Quat4d par = torsParams[it];
        // vectors between atoms x--y--z--w 
        Vec3d ha    = apos[ ias.x ] - apos[ ias.y ]; // ha  = x-y
        Vec3d hb    = apos[ ias.w ] - apos[ ias.z ]; // hb  = w-z
        Vec3d hab   = apos[ ias.z ] - apos[ ias.y ]; // hab = z-y
        double ila  = 1/ha.normalize();
        double ilb  = 1/hb.normalize();
        double ilab = 1/hab.normalize();
        double ca   = hab.dot(ha); // ca  = cos(alpha) = <  ha | hab > 
        double cb   = hab.dot(hb); // cb  = cos(beta ) = <  hb | hab >
        double cab  = ha .dot(hb); // cab = <  ha | hb  >
        double sa2  = (1-ca*ca);   // sa2 = sin^2( alpha ) = 1 - cos^2( alpha )
        double sb2  = (1-cb*cb);   // sb2 = sin^2( beta  ) = 1 - cos^2( beta  )
        double invs = 1/sqrt( sa2*sb2 );
        Vec2d cs,csn;
        cs.x = ( cab - ca*cb )*invs;   
        cs.y = sqrt(1-cs.x*cs.x);      // can we avoid this sqrt ?
        cs.udiv_cmplx( par.xy() );     // cs = cs0 * exp( i*phi )
        const int n = (int)par.w; // I know it is stupid store n as double, but I don't make another integer arrays just for it
        for(int i=0; i<n-1; i++){ // cs = cs0^n
            csn.mul_cmplx(cs);
        }
        // check here : https://www.wolframalpha.com/input/?i=(x+%2B+isqrt(1-x%5E2))%5En+derivative+by+x
        // energy
        const double k = par.z;
        double E       = k  *(1-csn.x);
        double dcn     = k*n*   csn.x;
        // derivatives to get forces
        double invs2 = invs*invs;
        dcn *= invs;
        double dcab  = dcn;                          // dc/dcab = dc/d<ha|hb>
        double dca   = (1-cb*cb)*(ca*cab - cb)*dcn;  // dc/dca  = dc/d<ha|hab>
        double dcb   = (1-ca*ca)*(cb*cab - ca)*dcn;  // dc/dca  = dc/d<hb|hab>
        Vec3d fa,fb,fab;
        fa =Vec3dZero;
        fb =Vec3dZero;
        fab=Vec3dZero;
        // Jacobina: derivative of <ha|hb> = <ha|hab> + <hb|hab> 
        SMat3d J;
        // derivative by ha
        J.from_dhat(ha);    // -- by ha
        J.mad_ddot(hab,fa, dca ); // dca /dha = d<ha|hab>/dha
        J.mad_ddot(hb ,fa, dcab); // dcab/dha = d<ha|hb> /dha
        J.from_dhat(hb);    // -- by hb
        J.mad_ddot(hab,fb, dcb ); // dcb /dhb = d<hb|hab>/dha
        J.mad_ddot(ha ,fb, dcab); // dcab/dhb = d<hb|ha> /dha
        J.from_dhat(hab);         // -- by hab
        J.mad_ddot(ha,fab, dca);  // dca/dhab = d<ha|hab>/dhab
        J.mad_ddot(hb,fab, dcb);  // dcb/dhab = d<hb|hab>/dhab
        // multiply by inverse lengths
        fa .mul( ila  );
        fb .mul( ilb  );
        fab.mul( ilab );
        // apply forces and recoils to atoms
        fapos[ias.x].sub( fa );
        fapos[ias.y].add( fa  - fab );
        fapos[ias.z].add( fab - fb  );
        fapos[ias.w].add( fb );
        return E;
    }
    */

    //  compute torsion energy and forces for all torsion angles
    /*
    double eval_torsions(){
        printf( "MMFFsp3_loc::eval_torsions(ntors=%i)\n", ntors );
        double E=0;
        for(int it=0; it<ntors; it++){ 
            E+=eval_torsion(it); 
        }
        return E;
    }
    */

    // constrain atom to fixed position
    /*
    void constrainAtom( int ia, double Kfix=1.0 ){
        printf( "constrainAtom(i=%i,K=%g)\n", ia, Kfix );
        constr[ia].f=apos[ia];
        constr[ia].w=Kfix;
    };
    */

    // add neighbor recoil forces to an atom (ia)
    /*
    void assemble_atom(int ia){
        Vec3d fa=Vec3dZero,fp=Vec3dZero;
        const int* ings = bkneighs[ia].array;
        bool bpi = ia<nnode; // if this is node atom, it has pi-orbital
        for(int i=0; i<4; i++){
            int j = ings[i];
            if(j<0) break;
            fa.add(fneigh  [j]);
            if(bpi){
                fp.add(fneighpi[j]);
            }
        }
        fapos [ia].add( fa ); 
        if(bpi){ // if this is node atom, apply recoil to pi-orbital as well
            fpipos[ia].add( fp );
            fpipos[ia].makeOrthoU( pipos[ia] );  // subtract force component which change pi-vector size
        }
    }
    */

    // add neighbor recoil forces to all atoms
    /*
    void asseble_forces(){
        for(int ia=0; ia<natoms; ia++){
            assemble_atom(ia);
        }
    }
    */

    // make list of back-neighbors
    /*
    void makeBackNeighs( bool bCapNeighs=true ){
        for(int i=0; i<natoms; i++){ bkneighs[i]=Quat4i{-1,-1,-1,-1}; } // clear back-neighbors to -1
        for(int ia=0; ia<nnode; ia++){
            for(int j=0; j<4; j++){        // 4 neighbors
                int ja = neighs[ia].array[j];
                if( ja<0 )continue;
                bool ret = addFirstEmpty( bkneighs[ja].array, 4, ia*4+j, -1 ); // add neighbor to back-neighbor-list on the first empty place
                if(!ret){ printf("ERROR in MMFFsp3_loc::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
            }
        }
        if(bCapNeighs){   // set neighbors for capping atoms
            for(int ia=nnode; ia<natoms; ia++){ neighs[ia]=Quat4i{-1,-1,-1,-1};  neighs[ia].x = bkneighs[ia].x/4;  }
        }
    }
    */

    // ================== Print functions  

    /*
    void printSizes     (      ){ printf( "MMFFf4::printSizes(): nDOFs(%i) natoms(%i) nnode(%i) ncap(%i) nvecs(%i) npbc(%i)\n", nDOFs,natoms,nnode,ncap,nvecs,npbc ); }
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

