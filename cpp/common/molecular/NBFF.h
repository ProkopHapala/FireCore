/*
Non-Bonded Force-Field

    should be easily plugged with any molecular dynamics, by either sharing pointer to same data buffer, or by copying data
    
*/

#ifndef NBFF_h
#define NBFF_h

#include "fastmath.h"
//#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Atoms.h"

#include "Buckets.h"
#include "Forces.h"
#include "ForceField.h"

#include "simd.h"

void fitAABB( Vec6d& bb, int n, int* c2o, Vec3d* ps ){
    //Quat8d bb;
    //bb.lo = bb.lo = ps[c2o[0]];
    for(int i=0; i<n; i++){ 
        //printf( "fitAABB() i %i \n", i );
        int ip = c2o[i];
        //printf( "fitAABB() i=%i ip=%i \n", i, ip );
        Vec3d p = ps[ip];
        bb.lo.setIfLower  ( p );
        bb.hi.setIfGreater( p );
    }; 
    //return bb;
}

int makePBCshifts_( Vec3i nPBC, const Mat3d& lvec, Vec3d*& shifts ){
    const int npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
    //printf( "makePBCshifts_() npbc=%i nPBC{%i,%i,%i}\n", npbc, nPBC.x,nPBC.y,nPBC.z );
    if(shifts==0)_realloc(shifts,npbc);
    int ipbc=0;
    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ 
        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ 
            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
                shifts[ipbc] = (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz);   
                ipbc++; 
            }
        }
    }
    return npbc;
}

__attribute__((pure))
__attribute__((hot))
Quat4d evalPointCoulPBC( Vec3d pos, int npbc, const Vec3d* shifts, int natoms, const Vec3d * apos, const double* Qs, double Rdamp ){
    const double R2damp=Rdamp*Rdamp;    
    //const double K=-alphaMorse;
    Quat4d qe     = Quat4dZero;
    //#pragma omp for simd
    for(int ia=0; ia<natoms; ia++){
        const Vec3d dp0   = pos - apos[ia];
        const double Qi = Qs[ia];
        //if( (ibuff==0) ){ printf( "DEBUG a[%i] p(%g,%g,%g) Q %g \n", ia,apos_[ia].x, apos_[ia].y, apos[ia].z, REQi.z ); }              
        for(int ipbc=0; ipbc<npbc; ipbc++ ){
            const Vec3d  dp = dp0 + shifts[ipbc];
            const double r2     = dp.norm2();
            const double ir2    = 1/(r2+R2damp);
            const double eQ     = COULOMB_CONST*Qi*sqrt(ir2);
            qe.e+=eQ;     qe.f.add_mul( dp,  eQ*ir2 ); // Coulomb
        }
    }
    return qe;
}


void sampleCoulombPBC( int nps, const Vec3d* ps, Quat4d* fe, int natom,  Vec3d* apos, double* Qs, Mat3d lvec, Vec3i nPBC, double Rdamp ){
    Vec3d * shifts =0;
    int npbc = makePBCshifts_( nPBC, lvec, shifts );
    printf( "sampleCoulombPBC() npbc=%i nPBC{%i,%i,%i} nps=%i natom=%i \n", npbc, nPBC.x,nPBC.y,nPBC.z, nps, natom );
    for(int ip=0; ip<nps; ip++){
        Vec3d pos = ps[ip];
        fe[ip] = evalPointCoulPBC( pos, npbc, shifts, natom, apos, Qs, Rdamp );
    };
    delete[] shifts;
}






// Force-Field for Non-Bonded Interactions
class NBFF: public ForceField{ public:
    
    // ---  inherited from Atoms
    //int     natoms =0; // from Atoms
    //int    *atypes =0; // from Atoms
    //Vec3d  *apos   =0; // from Atoms
    //Vec3d  *fapos  =0; // forces on atomic positions
    //Vec3d  *vapos  = 0; // [natom]  velocities of atoms

    Quat4d   *REQs  __attribute__((aligned(64))) =0; // non-bonding interaction paramenters (R: van dew Waals radius, E: van dew Waals energy of minimum, Q: Charge, H: Hydrogen Bond pseudo-charge )
    Quat4i   *neighs   =0; // list of neighbors (4 per atom)
    Quat4i   *neighCell=0; // list of neighbors (4 per atom)

    //  --- Use this to speed-up short range interaction (just repulsion no L-J or Coulomb)
    //       * there can be additional hydrogen bonds just for selected pairs of atoms
    //       * use function repulsion_R4() to calculate repulsion without need of square root ( it is defined in Forces.h )
    // Bounding Boxes
    int      nBBs=0;
    Vec6d*   BBs=0; // bounding boxes (can be either AABB, or cylinder, capsula) 
    Buckets  pointBBs;    // buckets for collision detection

    // --- Parameters
    //bool    bClampNonBonded =  0.0; // if >0 then we clamp non-bonded forces to this value
    //double  FmaxNonBonded   = 10.0; // if bClampNonBonded>0 then we clamp non-bonded forces to this value

    double drSR  = 0.5;     // R_SR = R_cut - drSR
    double ampSR = 0.15;   // Amplitude of short-range repulsion  = ampSR * EvdW

    double alphaMorse = 1.5; // alpha parameter for Morse potential
    //double  KMorse  = 1.5; // spring constant for Morse potential
    double  Rdamp     = 1.0; // damping radius for LJQ and MorseQ
    Mat3d   lvec __attribute__((aligned(64)));  // lattice vectors
    Vec3i   nPBC;  // number of periodic images in each direction 
    bool    bPBC=false; // periodic boundary conditions ?

    int    npbc   =0;  // total number of periodic images
    Vec3d* shifts __attribute__((aligned(64))) =0;  // array of bond vectors shifts in periodic boundary conditions
    Quat4f *PLQs  __attribute__((aligned(64))) =0;  // non-bonding interaction paramenters in PLQ format form (P: Pauli strenght, L: London strenght, Q: Charge ), for faster evaluation in factorized form, especially when using grid

    Quat4d *PLQd  __attribute__((aligned(64))) =0; 

    Vec3d  shift0 __attribute__((aligned(64))) =Vec3dZero; 

#ifdef WITH_AVX
    // ========== Try SIMD
    Vec3sd* apos_simd __attribute__((aligned(64))) =0; // forces on atomic positions
    Vec4sd* REQs_simd __attribute__((aligned(64))) =0; // non-bonding interaction paramenters (R: van dew Waals radius, E: van dew Waals energy of minimum, Q: Charge, H: Hydrogen Bond pseudo-charge )
#endif //WITH_AVX
    // ==================== Functions


    // calculate total torque on the molecule (with respect to point p0) 
    void torq     ( const Vec3d& p0,  Vec3d& tq                   ){ for(int i=0; i<natoms; i++){ Vec3d d; d.set_sub(apos[i],p0); tq.add_cross(fapos[i],d); } }

    // bind PBC-shift vectors 
    void bindShifts(int npbc_, Vec3d* shifts_ ){ npbc=npbc_; shifts=shifts_; }

    // make PBC-shift vectors
    int makePBCshifts( Vec3i nPBC_, bool bRealloc=true ){
        bPBC=true;
        nPBC=nPBC_;
        if(bRealloc){ _dealloc(shifts); }
        return makePBCshifts_(nPBC,lvec,shifts);
        /*
        npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
        if(bRealloc) _realloc(shifts,npbc);
        int ipbc=0;
        for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
            shifts[ipbc] = (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz);   
            ipbc++; 
        }}}
        if(npbc!=ipbc){ printf( "ERROR in MMFFsp3_loc::makePBCshifts() final ipbc(%i)!=nbpc(%i) => Exit()\n", ipbc,npbc ); exit(0); }
        return npbc;
        */
    }

    // pre-calculates PLQs from REQs (for faster evaluation in factorized form, especially when using grid)
    void evalPLQs(double K){
        //printf( "NBFF::evalPLQs() \n" );
        if(PLQs==0){ _realloc(PLQs,natoms);  }
        for(int i=0; i<natoms; i++){
            //printf( "makePLQs[%i] \n", i );
            //printf( "makePLQs[%i] REQ(%g,%g,%g) \n", i, REQs[i].x,REQs[i].y,REQs[i].z);
            PLQs[i]=REQ2PLQ( REQs[i], K );
            //printf( "makePLQs[%i] REQ(%g,%g,%g) PLQ(%g,%g,%g)\n", i, REQs[i].x,REQs[i].y,REQs[i].z,  PLQs[i].x,PLQs[i].y,PLQs[i].z );
        }
        //printf("NBFF::makePLQs() DONE => exit(0) \,"); exit(0);
    }
    void makePLQs(double K){
        _realloc(PLQs,natoms);
        evalPLQs(K);
    }


    // pre-calculates PLQs from REQs (for faster evaluation in factorized form, especially when using grid)
    void evalPLQd(double K){
        if(PLQd==0){ _realloc(PLQd,natoms);  }
        for(int i=0; i<natoms; i++){
            PLQd[i]=REQ2PLQ_d( REQs[i], K );
        }
    }
    void makePLQd(double K){
        _realloc(PLQd,natoms);
        evalPLQd(K);
    }

    __attribute__((hot))  
    inline void updatePointBBs( bool bInit=true){
        const Buckets& buckets = pointBBs;
        //printf( "updatePointBBs() START \n" );
        for(int ib=0; ib<buckets.ncell; ib++){
            //printf( "updatePointBBs() ib %i \n", ib );
            if(bInit){ BBs[ib].lo = Vec3dmax; BBs[ib].hi = Vec3dmin; }
            int n = buckets.cellNs[ib];
            if(n>0){
                int i0 = buckets.cellI0s[ib];
                //printf( "updatePointBBs() ib %i n %i i0 %i \n", ib, n, i0 );
                fitAABB( BBs[ib], n, buckets.cell2obj+i0, vapos );
            }
        }
        //printf( "updatePointBBs() DONE \n" );
    }

    __attribute__((hot))  
    inline int selectInBox( const Vec6d& bb, const int ib, Vec3d* ps, Quat4d* paras, int* inds ){
        const int  npi = pointBBs.cellNs[ib];
        const int* ips = pointBBs.cell2obj +pointBBs.cellI0s[ib];
        int  n   = 0;
        for(int i=0; i<npi; i++){
            int ia = ips[i];
            const Vec3d& p = vapos[ ia ];
            if( (p.x>bb.lo.x)&&(p.x<bb.hi.x)&&
                (p.y>bb.lo.y)&&(p.y<bb.hi.y)&&
                (p.z>bb.lo.z)&&(p.z<bb.hi.z) 
            ){
                ps[i]    = p;
                paras[i] = REQs[ ips[i] ];
                inds[i]  = ia;
                n++;
            }
        }
        return n;
    }

    __attribute__((hot))  
    double evalSortRange_BBs( double Rcut, int ngmax ){
        //printf( "evalSortRange_BBs() START \n" );
        double E=0;
        const double R2cut = Rcut*Rcut;
        Quat4d REQis[ngmax];
        Quat4d REQjs[ngmax];
        Vec3d  pis  [ngmax];
        Vec3d  pjs  [ngmax];
        int   nis   [ngmax];
        int   njs   [ngmax];

        for(int ib=0; ib<nBBs; ib++){

            // --- Within the same bucket
            const int  ni =  selectInBox( BBs[ib], ib, pis, REQis, nis );  // this is maybe not optimal (we just select all atoms in the bucket)
            for(int ip=0; ip<ni; ip++){
                const Vec3d&  pi   = pis[ip];
                const Quat4d& REQi = REQis[ip];
                const int     ia   = nis[ip];
                for(int jp=ip+1; jp<ni; jp++){
                    const Vec3d& d  = pjs[jp] - pi;
                    int          ja = nis[jp];
                    Quat4d REQij; combineREQ( REQis[ip], REQis[jp], REQij );
                    Vec3d f; repulsion_R4( d, f, REQij.x-drSR, REQij.x, REQij.y*ampSR );
                    fapos[ia].sub(f);
                    fapos[ja].add(f);
                }
            }

            // --- between different buckets 
            for(int jb=0; jb<ib; jb++){
                Vec6d bb;
                bb = BBs[jb]; bb.lo.add(-Rcut); bb.hi.add(Rcut); const int ni =  selectInBox( bb, ib, pis, REQis, nis );  // atoms in bucket i overlapping with bucket j
                bb = BBs[ib]; bb.lo.add(-Rcut); bb.hi.add(Rcut); const int nj =  selectInBox( bb, jb, pjs, REQjs, njs );  // atoms in bucket j overlapping with bucket i
                for(int jp=0; jp<nj; jp++){
                    const Vec3d&  pj   = pjs[jp];
                    const Quat4d& REQj = REQjs[jp];
                    const int     ja   = nis[jp];
                    for(int ip=0; ip<ni; ip++){
                        const Vec3d& d = pis[jp] - pj;
                        const int ia   = nis[ip];
                        Quat4d REQij; combineREQ( REQis[ip], REQj, REQij );
                        Vec3d f; repulsion_R4( d, f, REQij.x-drSR, REQij.x, REQij.y*ampSR );
                        fapos[ia].sub(f);
                        fapos[ja].add(f);
                    }
                }
            }

        }
        //printf( "evalSortRange_BBs() DONE \n" );
        return E;

    }

    // evaluate non-bonding interaction using Lenard-Jones potential and Coulomb potential
    __attribute__((hot))  
    double evalLJQs( double Rdamp=1.0 ){
        //printf( "NBFF::evalLJQs() \n" );
        double R2damp = Rdamp*Rdamp;
        const int N=natoms;
        double E=0;
        for(int i=0; i<N; i++){
            Vec3d fi = Vec3dZero;
            const Vec3d pi = apos[i];
            const Quat4d& REQi = REQs[i];
            for(int j=i+1; j<N; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Quat4d REQij; combineREQ( REQs[j], REQi, REQij ); // combine non-bonding interaction parameters of atoms i and j
                //E += addAtomicForceLJQ( apos[j]-pi, fij, REQij );
                E += getLJQH( apos[j]-pi, fij, REQij, R2damp ); // calculate non-bonding interaction energy and force using Lenard-Jones potential and Coulomb potential and 
                fapos[j].sub(fij);
                fi   .add(fij);
            }
            fapos[i].add(fi);
        }
        return E;
    }

    __attribute__((hot))  
    double evalLJQs_PBC( const Mat3d& lvec, Vec3i nPBC=Vec3i{1,1,1}, double Rdamp=1.0 ){
        //printf( "NBFF::evalLJQs_PBC() \n" );
        double R2damp = Rdamp*Rdamp;
        const int N   = natoms;
        double E=0;
        int npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1) -1;
        Vec3d shifts[npbc]; // temporary store for lattice shifts
        int ipbc=0;
        // initialize shifts
        for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){ 
            if((ia==0)&&(ib==0)&&(ic==0))[[unlikely]]{  continue; } // skipp pbc0
            shifts[ipbc] = (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic);   
            ipbc++; 
        }}}
        // calculate non-bonding interaction
        for(int i=0; i<N; i++){
            Vec3d fi = Vec3dZero;
            const Vec3d     pi = apos[i];
            const Quat4d& REQi = REQs[i];
            for(int j=i; j<N; j++){    // atom-atom (yes self interaction, no double-counting)
                Vec3d fij = Vec3dZero;
                Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
                for(ipbc=0; ipbc<npbc; ipbc++){
                    //E += addAtomicForceLJQ( apos[j]-pi-shifts[ipbc], fij, REQij );
                    E +=getLJQH( apos[j]-pi-shifts[ipbc], fij, REQij, R2damp );
                }
                fapos[j].sub(fij);
                fi   .add(fij);
            }
            fapos[i].add(fi);
        }
        //printf( "npbc %i ipbc %i E %g \n", npbc, ipbc, E );
        return E;
    }

    // evaluate non-bonding interaction for given atom (ia) excluding bonded atoms, assume max 4 bonds per atom
    __attribute__((hot))  
    double evalLJQs_ng4_atom( int ia ){
        //printf( "NBFF::evalLJQs_ng4_atom() %li %li \n" );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi   = apos  [ia];           // global read   apos  [ia]
        const Quat4d  REQi = REQs  [ia];           // global read   REQs  [ia]
        const Quat4i ng   = neighs[ia];           // global read   neighs[ia]
        double       E    = 0;
        Vec3d        fi   = Vec3dZero;
        for(int j=0; j<natoms; j++){
            if(ia==j) continue;
            if( (ng.x==j)||(ng.y==j)||(ng.z==j)||(ng.w==j) ) [[unlikely]] { continue; }
            Vec3d fij   = Vec3dZero;
            const Vec3d pj     = apos[j];                       // global read   apos[j]
            const Quat4d REQj  = REQs[j];                       // global read   REQs[j]
            Quat4d REQij; combineREQ( REQj, REQi, REQij );      // pure function (no globals)
            E          += getLJQH( pj-pi, fij, REQij, R2damp );  // pure function (no globals)
            fi.add(fij);
        }
        fapos[ia].add(fi);           // global write fapos[ia]
        return E;
    }

    // evaluate all non-bonding interactions excluding bonded atoms, with OpenMP parallelization
    __attribute__((hot))  
    double evalLJQs_ng4_omp( ){
        //printf( "NBFF::evalLJQs_ng4_omp() \n" );
        double E =0;
        //#pragma omp parallel for reduction(+:E) shared(neighs,R2damp) schedule(dynamic,1)
        for(int ia=0; ia<natoms; ia++){
           //printf("omp_get_num_threads=%i\n", omp_get_num_threads() );
           //printf("eval_atoms(%i) @cpu[%d/%d]\n", i, omp_get_thread_num(), omp_get_num_threads() );
           E += evalLJQs_ng4_atom( ia );
        }

        return E;
    }

    // evaluate all non-bonding interactions excluding bonded atoms (assume max 4 bonds per atom), single-threaded version
    __attribute__((hot))  
    double evalLJQs_ng4( const Quat4i* neighs, double Rdamp=1.0 ){
        //printf( "NBFF::evalLJQs_ng4() \n" );
        double R2damp = Rdamp*Rdamp;
        const int N=natoms;
        double E   =0;
        int i=0;
        //printf( "NBFF::evalLJQs_ng4() neighs=%li \n", neighs );
        for(i=0; i<N; i++){
            Vec3d fi = Vec3dZero;
            const Vec3d  pi = apos[i];
            const Quat4d& REQi = REQs[i];
            //const int*   ngs = neighs+i*4;
            const Quat4i& ngs  = neighs[i];
            //printf( "NBFF::evalLJQs_ng4()[%i] ngs(%i,%i,%i,%i) \n", i, ngs.x,ngs.y,ngs.z,ngs.w );
            for(int j=i+1; j<N; j++){    // atom-atom (no self interaction, no double-counting)
                //printf( "NBFF::evalLJQs_ng4()[%i,%j] neighs=%li \n", neighs );
                if( (ngs.x==j)||(ngs.y==j)||(ngs.z==j)||(ngs.w==j) ) [[unlikely]]  {continue; }
                Vec3d fij = Vec3dZero;
                Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
                //E += addAtomicForceLJQ( apos[j]-pi, fij, REQij );
                E +=getLJQH( apos[j]-pi, fij, REQij, R2damp );
                fapos[j].sub(fij);
                fi   .add(fij);
            }
            fapos[i].add(fi);
        }
        return E;
    }


    // evaluate all non-bonding interactions using Morse potential and Coulomb potential in periodic boundary conditions with OpenMP parallelization
    __attribute__((hot))  
    inline double addMorseQH_PBC_omp( Vec3d pi, const Quat4d&  REQi, Vec3d& fout ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        pi.sub(shift0);
        const double R2damp = Rdamp*Rdamp;
        double E=0,fx=0,fy=0,fz=0;
        #pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            //if(ia==j)continue;    ToDo: Maybe we can keep there some ignore list ?
            const Quat4d& REQj  = REQs[j];
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp      = apos[j]-pi;
            Vec3d fij           = Vec3dZero;
            for(int ipbc=0; ipbc<npbc; ipbc++){
                const Vec3d dpc = dp + shifts[ipbc];
                double eij      = getMorseQH( dpc, fij, REQij, alphaMorse, R2damp );
                E +=eij;
                fx+=fij.x;
                fy+=fij.y;
                fz+=fij.z;
            }
        }
        fout.add(fx,fy,fz);
        return E;
    }
    inline double getMorseQH_PBC_omp( Vec3d pi, const Quat4d&  REQi, Vec3d& fout ){ fout=Vec3dZero; return addMorseQH_PBC_omp( pi, REQi, fout ); }

    // evaluate all non-bonding interactions using Lenard-Jones potential and Coulomb potential in periodic boundary conditions with OpenMP parallelization
    __attribute__((hot))  
    double getLJQs_PBC_omp( const Vec3d& pi, const Quat4d&  REQi, Vec3d& fout ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        double E=0,fx=0,fy=0,fz=0;
        #pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            //if(ia==j)continue;   ToDo: Maybe we can keep there some ignore list ?
            const Quat4d& REQj  = REQs[j];
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp      = apos[j]-pi;
            Vec3d fij           = Vec3dZero;
            for(int ipbc=0; ipbc<npbc; ipbc++){
                // --- We calculate non-bonding interaction every time (most atom pairs are not bonded)
                const Vec3d dpc = dp + shifts[ipbc];    //   dp = pj - pi + pbc_shift = (pj + pbc_shift) - pi 
                double eij      = getLJQH( dpc, fij, REQij, R2damp );
                //printf( "getLJQs_PBC_omp[%i] dp(%6.3f,%6.3f,%6.3f) REQ(%g,%g,%g,%g) \n", eij, dp.x,dp.y,dp.z, REQij.x,REQij.y,REQij.z,REQij.w );
                E +=eij;
                fx+=fij.x;
                fy+=fij.y;
                fz+=fij.z;
            }
        }
        fout=Vec3d{fx,fy,fz};
        return E;
    }

    // evaluate all non-bonding interactions using Lenard-Jones potential and Coulomb potential in periodic boundary conditions with OpenMP SIMD parallelizati
    __attribute__((hot))  
    double evalLJQs_PBC_atom_omp( const int ia, const double Fmax2 ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const Vec3d   pi      = apos[ia];
        const Quat4d  REQi    = REQs[ia];
        const double  R2damp = Rdamp*Rdamp;
        double E=0,fx=0,fy=0,fz=0;
        #pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            //if(ia==j)continue;   ToDo: Maybe we can keep there some ignore list ?
            if(ia==j)[[unlikely]]{continue;}
            const Quat4d& REQj  = REQs[j];
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp      = apos[j]-pi;
            Vec3d fij           = Vec3dZero;
            for(int ipbc=0; ipbc<npbc; ipbc++){
                // --- We calculate non-bonding interaction every time (most atom pairs are not bonded)
                const Vec3d dpc = dp + shifts[ipbc];    //   dp = pj - pi + pbc_shift = (pj + pbc_shift) - pi 
                double eij      = getLJQH( dpc, fij, REQij, R2damp );
                if(bClampNonBonded)[[likely]]{ clampForce( fij, Fmax2 ); }
                //printf( "getLJQs_PBC_omp[%i] dp(%6.3f,%6.3f,%6.3f) REQ(%g,%g,%g,%g) \n", eij, dp.x,dp.y,dp.z, REQij.x,REQij.y,REQij.z,REQij.w );
                E +=eij;
                fx+=fij.x;
                fy+=fij.y;
                fz+=fij.z;
            }
        }
        fapos[ia].add( Vec3d{fx,fy,fz} );
        return E;
    }
    __attribute__((hot))  
    double evalLJQs_PBC_simd(){
        //printf("NBFF::evalLJQs_PBC_simd()\n" );
        double E=0;
        const double Fmax2 = FmaxNonBonded*FmaxNonBonded;
        for(int ia=0; ia<natoms; ia++){  
            //printf("ffls[%i].evalLJQs_PBC_simd(%i)\n", id, ia ); 
            E+=evalLJQs_PBC_atom_omp( ia, Fmax2 ); 
        }
        return E;
    }

    // evaluate all non-bonding interactions using Lenard-Jones potential and Coulomb potential in periodic boundary conditions with OpenMP SIMD parallelization
    __attribute__((hot))  
    double evalLJQs_ng4_PBC_atom_omp(const int ia ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        //printf("DEBUG 1 id=%i ia=%i REQs=%li \n", id, ia, REQs );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi   = apos     [ia];
        const Quat4d REQi = REQs    [ia];
        const Quat4i ng   = neighs   [ia];
        const Quat4i ngC  = neighCell[ia];
        double E=0,fx=0,fy=0,fz=0;
        //printf("DEBUG 1 id=%i ia=%i \n", id, ia );
        //#pragma omp simd collapse(2) reduction(+:E,fx,fy,fz)
        #pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            if(ia==j)continue;
            const Quat4d& REQj  = REQs[j];
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp     = apos[j]-pi;
            Vec3d fij          = Vec3dZero;
            const bool bBonded = ((j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w));
            //printf("DEBUG 2 id=%i ia=%i j=%i \n", id, ia, j );
            for(int ipbc=0; ipbc<npbc; ipbc++){
                // --- We calculate non-bonding interaction every time (most atom pairs are not bonded)
                const Vec3d dpc = dp + shifts[ipbc];    //   dp = pj - pi + pbc_shift = (pj + pbc_shift) - pi 
                double eij      = getLJQH( dpc, fij, REQij, R2damp );
                // --- If atoms are bonded we don't use the computed non-bonding interaction energy and force
                if(bBonded) [[unlikely]]  { 
                    if(   ((j==ng.x)&&(ipbc==ngC.x))
                        ||((j==ng.y)&&(ipbc==ngC.y))
                        ||((j==ng.z)&&(ipbc==ngC.z))
                        ||((j==ng.w)&&(ipbc==ngC.w))
                    ) [[unlikely]]  { 
                        continue;
                    }
                }
                E +=eij;
                fx+=fij.x;
                fy+=fij.y;
                fz+=fij.z;
                //fi+=fij;
            }
        }
        //printf("DEBUG 3 id=%i ia=%i \n", id, ia );
        fapos[ia].add( Vec3d{fx,fy,fz} );
        return E;
    }
    __attribute__((hot))  
    double evalLJQs_ng4_PBC_simd(){
        //printf("NBFF::evalLJQs_ng4_PBC_simd()\n" );
        double E=0;
        for(int ia=0; ia<natoms; ia++){  
            //printf("ffls[%i].evalLJQs_ng4_PBC_atom_omp(%i)\n", id, ia ); 
            E+=evalLJQs_ng4_PBC_atom_omp(ia); 
        }
        return E;
    }
    __attribute__((hot))  
    double evalLJQs_atom_omp( const int ia, const double Fmax2 ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi   = apos     [ia];
        const Quat4d REQi = REQs     [ia];
        Vec3d fi = Vec3dZero;
        double E=0,fx=0,fy=0,fz=0;
        //#pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            if(ia==j)[[unlikely]]{continue;}
            const Quat4d& REQj  = REQs[j];
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp      = apos[j]-pi;
            Vec3d fij           = Vec3dZero;
            double eij = getLJQH( dp, fij, REQij, R2damp );
            if(bClampNonBonded)[[likely]]{ clampForce( fij, Fmax2 ); }
            E +=eij;
            fx+=fij.x;
            fy+=fij.y;
            fz+=fij.z;
            //fi+=fij;
            // {
            //     glLineWidth(3.0);
            //     glColor3d(0.0,0.0,1.0);  Draw3D::drawVecInPos( dp,  apos[ia] );
            //     glColor3d(1.0,0.0,0.0);  Draw3D::drawVecInPos( fij, apos[ia] );
            // } 
            //printf( "evalLJQs_atom_omp(%i) E %g f(%g,%g,%g) \n", ia, E, fx, fy, fz );
        }
        //printf( "evalLJQs_atom_omp(%i) E %g f(%g,%g,%g) \n", ia, E, fx, fy, fz );
        fapos[ia].add( Vec3d{fx,fy,fz} );
        return E;
    }
    __attribute__((hot))  
    double evalLJQs_simd(){
        //printf("NBFF::evalLJQs_simd()\n" );
        double E=0;
        const double Fmax2 = FmaxNonBonded*FmaxNonBonded;
        for(int ia=0; ia<natoms; ia++){ E+=evalLJQs_atom_omp(ia, Fmax2 ); }
        return E;
    }

    __attribute__((hot))  
    double evalLJQs_ng4_atom_omp( const int ia ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi   = apos     [ia];
        const Quat4d  REQi = REQs     [ia];
        const Quat4i ng   = neighs   [ia];
        Vec3d fi = Vec3dZero;
        double E=0,fx=0,fy=0,fz=0;

        #pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            if( (ia==j)  || (j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w) ) [[unlikely]]  { continue; }
            const Quat4d& REQj  = REQs[j];
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp      = apos[j]-pi;
            Vec3d fij           = Vec3dZero;
            double eij = getLJQH( dp, fij, REQij, R2damp );
            E +=eij;
            fx+=fij.x;
            fy+=fij.y;
            fz+=fij.z;
            //fi+=fij; 
        }
        fapos[ia].add( Vec3d{fx,fy,fz} );
        return E;
    }
    __attribute__((hot))  
    double evalLJQs_ng4_simd(){
        //printf("NBFF::evalLJQs_ng4_simd()\n" );
        double E=0;
        for(int ia=0; ia<natoms; ia++){ E+=evalLJQs_ng4_atom_omp(ia); }
        return E;
    }

    __attribute__((hot))  
    Quat4d evalLJQs( Vec3d pi, Quat4d REQi, double Rdamp )const{
        const double R2damp = Rdamp*Rdamp;
        Quat4d fe = Quat4dZero;
        for(int i=0; i<natoms; i++){
            const Quat4d  REQij = _mixREQ(REQi,(REQs[i])); 
            Vec3d dp = pi - apos[i];
            Vec3d fij;
            fe.e += getLJQH( dp, fij, REQij, R2damp );
            fe.f.add(fij);
        }
        return fe;
    }

#ifdef WITH_AVX
    __attribute__((hot))  
    double evalLJQs_atom_avx( const int ia, const double Fmax2 ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi_   = apos[ia];
        
        const Vec3sd pi = Vec3sd{ 
            _mm256_set_pd( pi_.x, pi_.x, pi_.x, pi_.x ), 
            _mm256_set_pd( pi_.y, pi_.y, pi_.y, pi_.y ), 
            _mm256_set_pd( pi_.z, pi_.z, pi_.z, pi_.z ) 
        };
        const Quat4d REQi_ = REQs[ia];
        const Vec4sd REQi = Vec4sd{ 
            _mm256_set_pd( REQi_.x, REQi_.x, REQi_.x, REQi_.x ), 
            _mm256_set_pd( REQi_.y, REQi_.y, REQi_.y, REQi_.y ), 
            _mm256_set_pd( REQi_.z, REQi_.z, REQi_.z, REQi_.z ), 
            _mm256_set_pd( REQi_.w, REQi_.w, REQi_.w, REQi_.w ) 
        };
        Vec3sd fi = Vec3sd{ _mm256_set_pd( 0.,0.,0.,0. ), _mm256_set_pd( pi_.y, pi_.y, pi_.y, pi_.y ), _mm256_set_pd( pi_.z, pi_.z, pi_.z, pi_.z ) };
        __m256d E = _mm256_set_pd( 0.,0.,0.,0. );
        for (int j=0; j<natoms; j++){ 
            const Vec4sd& REQj  = REQs_simd[j];

            const __m256d Rij = REQi.x+REQj.x;
            const __m256d Eij = REQi.y*REQj.y;
            const __m256d Qij = REQi.z*REQj.z;
            const __m256d Hij = REQi.w*REQj.w;
            const Vec3sd dp      = apos_simd[j]-pi;

            const __m256d  r2  = dp.norm2();
            __m256d F;
            // ---- Electrostatic
            const __m256d ir2_ = 1./( r2 + R2damp  );
            E +=  COULOMB_CONST* ( Qij*_mm256_sqrt_pd( ir2_ ) );
            F  =  E*ir2_ ;
            // --- LJ 
            const __m256d  ir2 = 1./r2;
            const __m256d  u2  = Rij*Rij*ir2;
            const __m256d  u6  = u2*u2*u2;
            const __m256d vdW  = u6*Eij;
            const __m256d   H  = u6*u6* ((Hij<0) ? Hij : 0.0);  // H-bond correction
            E   +=  (u6-2.)*vdW + H             ;
            F   += ((u6-1.)*vdW + H )*ir2*12 ;
            fi.add_mul( dp, -F );
        }
        fapos[ia].add( 
            hsum_double_avx( fi.x),
            hsum_double_avx( fi.y),
            hsum_double_avx( fi.z)       
        );
        return hsum_double_avx( E );
    }
#endif // WITH_AVX

    double evalCollisionDamp_atom_omp( const int ia, double damp_rate, double dRcut1=-0.2, double dRcut2=0.3 ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi     = apos [ia];
        const Vec3d  vi     = vapos[ia];
        const double Ri     = REQs [ia].x;
        Vec3d fi = Vec3dZero;
        double fx=0,fy=0,fz=0;

        //#pragma omp simd reduction(+:fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            double R = Ri + REQs[j].x;
            const Vec3d d  = apos[j]-pi;
            double r2 = d.norm2();
            if(j == ia) continue;
            // reduced_mass = mi*mj/(mi+mj) = mi*mj/m_cog
            // delta_v      = reduced_mass * dv
            // fi           = mi * delta_v/dt
            // for masses = 1.0 we have reduced_mass = 1*1/(1+1) = 0.5
            // fi           = 0.5 * dv/dt = 0.5 * damp_rate .... because damp_rate = 1/(dt*ndampstep)

            double w    = damp_rate * 0.5 * smoothstep_down(sqrt(r2), R+dRcut1, R+dRcut2 ); // ToDo : we can optimize this by using some other cutoff function which depends only on r2 (no sqrt)
            //double w    = damp_rate * 0.5 * R8down         (r2,      R, R+dRcut );
            double fcol = w * d.dot( vapos[j]-vi );                                // collisionDamping ~ 1/(dt*ndampstep);     f = m*a = m*dv/dt
            Vec3d fij; fij.set_mul( d, fcol/r2 ); //  vII = d*d.fot(v)/|d|^2 
            fx+=fij.x;
            fy+=fij.y;
            fz+=fij.z;
            //fi+=fij; 
        }
        fapos[ia].add( Vec3d{fx,fy,fz} );
        return 0;
    }
    double evalCollisionDamp_omp( double damp_rate, double dRcut=0.5 ){
        //printf("NBFF::evalCollisionDamp_omp()\n" );
        double E=0;
        for(int ia=0; ia<natoms; ia++){ E+=evalCollisionDamp_atom_omp(ia, damp_rate, dRcut ); }
        return E;
    }


    double evalLJQs_ng4_PBC_atom(const int ia ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        const bool   bPBC = npbc>1;
        const Vec3d  pi   = apos     [ia];
        const Quat4d& REQi = REQs     [ia];
        const Quat4i ng   = neighs   [ia];
        const Quat4i ngC  = neighCell[ia];
        Vec3d fi = Vec3dZero;
        double E =0 ;
        for (int j=0; j<natoms; j++){ if(ia==j)continue;  // all-to-all makes it easier to paralelize 
            const Vec3d dp = apos[j]-pi;
            Vec3d fij      = Vec3dZero;
            Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
            const bool bBonded = ((j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w));
            if(bPBC){
                for(int ipbc=0; ipbc<npbc; ipbc++){
                    if(bBonded){
                        if(   ((j==ng.x)&&(ipbc==ngC.x))
                            ||((j==ng.y)&&(ipbc==ngC.y))
                            ||((j==ng.z)&&(ipbc==ngC.z))
                            ||((j==ng.w)&&(ipbc==ngC.w))
                        ){ continue;}
                    }
                    const Vec3d dpc = dp + shifts[ipbc];    //   dp = pj - pi + pbc_shift = (pj + pbc_shift) - pi 
                    double eij = getLJQH( dpc, fij, REQij, R2damp );
                    E+=eij;
                    fi.add(fij);
                }
            }else{
                if(bBonded) continue;  // Bonded ?
                E+=getLJQH( dp, fij, REQij, R2damp );
                fi.add(fij);
            }
        }
        fapos[ia].add(fi);
        return E;
    }

    double evalLJQs_ng4_PBC_omp(){
        double E=0;
        for (int ia=0; ia<natoms; ia++ ){ 
            evalLJQs_ng4_PBC_atom( ia ); 
        }
        return E;
    }

    //double evalLJQs_ng4_PBC( Quat4i* neighs, Quat4i* neighCell, const Mat3d& lvec, Vec3i nPBC=Vec3i{1,1,1}, double Rdamp=1.0 ){
    double evalLJQs_ng4_PBC( Quat4i* neighs, Quat4i* neighCell, int npbc, const Vec3d* shifts, double Rdamp=1.0 ){
        //printf( "NBFF::evalLJQs_ng4_PBC() \n" );
        //printf( "evalLJQs_ng4_PBC() nPBC(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
        //int ia_DBG = 0;
        // printf( "evalLJQs_ng4_PBC() n=%i nPBC(%i,%i,%i) Rdamp %g \n lvec\n", n, nPBC.x,nPBC.y,nPBC.z, Rdamp );
        // printMat(lvec);
        // for(int i=0; i<n; i++){
        //     printf("evalLJQs[%i] ", i);
        //     printf("neighs(%i,%i,%i,%i) "   , neighs[i].x   ,neighs[i].y   ,neighs[i].z   ,neighs[i].w    );
        //     printf("neighCell(%i,%i,%i,%i) ", neighCell[i].x,neighCell[i].y,neighCell[i].z,neighCell[i].w );
        //     printf("apos(%6.3f,%6.3f,%6.3f) ",  apos[i].x ,apos[i].y ,apos[i].z  );
        //     printf("fapos(%6.3f,%6.3f,%6.3f) ", fapos[i].x,fapos[i].y,fapos[i].z );
        //     printf("REQ(%6.3f,%6.3f,%6.3f) ",   REQs[i].x ,REQs[i].y ,REQs[i].z  );
        //     printf("\n"); 
        // }
        
        const double R2damp = Rdamp*Rdamp;
        const int    n      = natoms;
        const bool   bPBC   = npbc>1;

        double E=0;
        //for(int i=0; i<n; i++)printf( "CPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQs[i].x,REQs[i].y,REQs[i].z );
        for (int i=0; i<n; i++ ){
            Vec3d fi = Vec3dZero;
            const Vec3d pi = apos[i];
            const Quat4d& REQi = REQs     [i];
            const Quat4i ng   = neighs   [i];
            const Quat4i ngC  = neighCell[i];
            //for (int j=i+1; j<n; j++){
            //if(i==4){ printf( "CPU_LJQ[%i] ng(%i,%i,%i,%i) ngC(%i,%i,%i,%i) npbc=%i\n", i, ng.x,ng.y,ng.z,ng.w,   ngC.x,ngC.y,ngC.z,ngC.w, npbc ); } 
            //printf( "CPU_LJQ[%i] ng(%i,%i,%i,%i) ngC(%i,%i,%i,%i) npbc=%i\n", i, ng.x,ng.y,ng.z,ng.w,   ngC.x,ngC.y,ngC.z,ngC.w, npbc );
            for (int j=0; j<n; j++){ if(i==j)continue;  // all-to-all makes it easier to paralelize 
            //for (int j=i+1; j<n; j++){
                const Vec3d dp = apos[j]-pi;
                Vec3d fij      = Vec3dZero;
                Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
                const bool bBonded = ((j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w));
                //const bool bBonded = false;
                if(bPBC){
                    for(int ipbc=0; ipbc<npbc; ipbc++){
                        if(bBonded){
                            if(   ((j==ng.x)&&(ipbc==ngC.x))
                                ||((j==ng.y)&&(ipbc==ngC.y))
                                ||((j==ng.z)&&(ipbc==ngC.z))
                                ||((j==ng.w)&&(ipbc==ngC.w))
                            ){
                                //printf("skip[%i,%i]ipbc=%i\n", i, j, ipbc );
                                continue; // skipp pbc0
                            }
                        }
                        const Vec3d dpc = dp + shifts[ipbc];    //   dp = pj - pi + pbc_shift = (pj + pbc_shift) - pi 

                        // if( (i==3)&&(j==1) ){ 
                        //     Vec3d shpi = apos[i] - shifts[ipbc];  
                        //     printf( "LJ(%i,%i)  ic %i shpi(%g,%g,%g) \n", i,j, ipbc,  shpi.x,shpi.y,shpi.z  );  
                        // }

                        double eij = getLJQH( dpc, fij, REQij, R2damp );
                        //if( (i==36)&&(j==35) )
                        //if( (i==36)&&(j==18) )
                        // if( eij<-0.3 )
                        // { 
                        //     double r = dpc.norm();
                        //     printf( "CPU_LJQ[%i,%i|%i] r,e,fr(%6.3f,%10.7f,%10.7f)    REQ(%6.3f,%10.7f,%4.2f) pbc_shift(%6.3f,%6.3f,%6.3f)\n" , i,j,ipbc, r, eij, dpc.dot(fij)/r,    REQij.x,REQij.y,REQij.z, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z );
                        //     //printf( "CPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g pbc_shift(%g,%g,%g)\n" , i,j, ipbc, fij.x,fij.y,fij.z, R2damp, REQij.x,REQij.y,REQij.z, dpc.norm(), shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z ); 
                            
                        // } 
                        // if( eij<-0.2 ){ Draw3D::drawVecInPos( dpc, pi ); };
                        E+=eij;
                        fi.add(fij);
                    }
                }else{
                    if(bBonded) continue;  // Bonded ?
                    E+=getLJQH( dp, fij, REQij, R2damp );
                    //if(i==ia_DBG){ printf( "CPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n" , i,j, fij.x,fij.y,fij.z, R2damp, REQij.x,REQij.y,REQij.z, dp.norm() ); } 
                    fi.add(fij);
                }
                //fapos[j].sub(fij);
            }
            //if(i==ia_DBG){ printf( "CPU_LJQ[%i] fapos(%g,%g,%g) += fi(%g,%g,%g) \n" , i, fapos[i].x,fapos[i].y,fapos[i].z, fi.x,fi.y,fi.z ); } 
            fapos[i].add(fi);
        }
        return E;
    }






    double evalLJQ( NBFF& B, const bool bRecoil, double Rdamp=1.0 ){
        //printf( "evalLJQ() \n" );
        const int n  = natoms;
        const int nb = B.natoms;
        const double R2damp = Rdamp*Rdamp;
        double E=0;
        //printf("DEBUG NBFF_AB.evalLJQ() n,m %i %i \n", n,m);
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            const Vec3d Api = apos[i];
            const Quat4d& REQi = REQs[i];
            for(int j=0; j<nb; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Quat4d REQij = _mixREQ( B.REQs[j], REQi );
                //E += addAtomicForceLJQ( B.apos[j]-Api, fij, REQij );
                E+=getLJQH( B.apos [j]-Api, fij, REQij, R2damp );
                if(bRecoil) B.fapos[j].sub(fij);
                fi.add(fij);
            }
            fapos[i].add(fi);
        }
        return E;
    }

    double evalMorse( NBFF& B, const bool bRecoil, double K=-1.0, double RQ=1.0 ){
        //printf( "NBFF::evalMorse() \n" );
        const int    n   = natoms;
        const int    nb  = B.natoms;
        const double R2Q = RQ*RQ;
        double E=0;
        //printf("evalMorse() n,B.n %i %i \n", n,B.n );
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            const Vec3d   pi   = apos[i];
            const Quat4d& REQi = REQs[i];
            for(int j=0; j<nb; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                Quat4d REQij; combineREQ( B.REQs[j], REQi, REQij );
                E += addAtomicForceMorseQ( B.apos[j]-pi, fij, REQij.x, REQij.y, REQij.z, K, R2Q );
                if(bRecoil) B.fapos[j].sub(fij);
                fi.add(fij);
            }
            //printf( "fapos[%i]( %g | %g,%g,%g)  \n", i, fi.norm(),  fi.x,fi.y,fi.z );
            fapos[i].add(fi);
            //fapos[i].sub(fi);
        }
        //printf( "EvdW %g \n", E );
        return E;
    }

    double evalMorsePBC( NBFF& B, const Mat3d& lvec, Vec3i nPBC=(Vec3i){1,1,0}, double K=-1.0, double RQ=1.0 ){
        //printf( "NBFF::evalMorsePBC() \n" );
        //nPBC = {0,0,0};
        //printf( "NBFF::evalMorsePBC() nPBC(%i,%i,%i) \n", nPBC.x, nPBC.y, nPBC.z );
        //printf( "cell.a(%g,%g,%g) cell.b(%g,%g,%g) cell.c(%g,%g,%g) \n", cell.a.x,cell.a.y,cell.a.z,   cell.b.x,cell.b.y,cell.b.z,   cell.c.x,cell.c.y,cell.c.z );
        //printf("evalMorse() n,B.n %i %i \n", n,B.n );
        const int    n   = natoms;
        const int    nb  = B.natoms;
        const double R2Q = RQ*RQ;
        double E=0;
        for(int i=0; i<n; i++){
            Vec3d fi  = Vec3dZero;
            const Vec3d pi_ = apos[i];
            const Quat4d& REQi = REQs[i];
            //if(i==0)printf("pi_(%g,%g,%g) \n",  pi_.x,pi_.y,pi_.z );
            for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                const Vec3d  pi = pi_ + (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic);
                //if(i==0)printf("pi[%2i,%2i,%2i] = (%g,%g,%g) \n", ia,ib,ic,   pi.x,pi.y,pi.z );
                for(int j=0; j<nb; j++){    // atom-atom
                    Vec3d fij = Vec3dZero;
                    Quat4d REQij; combineREQ( B.REQs[j], REQi, REQij );
                    E += addAtomicForceMorseQ( B.apos[j]-pi, fij, REQij.x, REQij.y, REQij.z, K, R2Q );
                    //E += addAtomicForceMorseQ( B.apos[j]-pi, fij, REQij.x, REQij.y, 0, K, R2Q );
                    //E += addAtomicForceQ_R2  ( B.apos[j]-pi, fij, REQij.z, K, R2Q );
                    //E=0; fij = B.apos[j]-pi;  // Test - dr
                    //if(i==0)printf("fi(%g,%g,%g) \n",  fij.x,fij.y,fij.z );
                    fi.add(fij);

                    // if((ia==0)&&(ib==0)&&(ic==0)){ // debug draw
                    //     const Vec3d shift = (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic);
                    //     Draw3D::drawLine( B.apos[j] - shift , pi_  );
                    // }

                }
            }}} // nPBC
            //if(i==0){ printf( "CPU atom[%i]  fe_Cou(%g,%g,%g|%g)  REQKi.z %g \n", i, fi.x,fi.y,fi.z,E, REQi.z ); }
            //fapos[i].add(fi);
        }
        return E;
    }

    double evalMorsePLQ( NBFF& B, Mat3d& cell, Vec3i nPBC, double K=-1.0, double RQ=1.0 ){
        //printf( "NBFF::evalMorsePLQ() PLQs %li \n", (long)PLQs, K, RQ );
        double E=0;
        //printf( "NBFF nPBC(%i,%i,%i) K %g RQ %g R2Q %g plq.z %g \n", nPBC.x,nPBC.y,nPBC.z, K, sqrt(R2Q), R2Q, plq.z );
        double R2Q=RQ*RQ;
        for(int i=0; i<natoms; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi = apos[i];
            //const Quat4d& REQi = REQs[i];
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int j=0; j<B.natoms; j++){
                Vec3d dp0; dp0.set_sub( pi, B.apos[j] );
                Quat4d REQj = B.REQs[j];
                for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){
                    Vec3d  dp = dp0 + cell.a*ia + cell.b*ib + cell.c*ic;
                    //Vec3d  dp     = dp0;
                    double r2     = dp.norm2();
                    double r      = sqrt(r2);
                    // ----- Morse
                    double e      = exp( K*(r-REQj.x) );
                    double de     = K*e*REQj.y*-2/r;
                    double eM     = e*REQj.y;
                    // ---- Coulomb
                    double ir2    = 1/(r2+R2Q);
                    double ir     = sqrt(ir2);
                    double eQ     = COULOMB_CONST*REQj.z*ir;
                    //if((i==0)&&(j==0))printf("dp(%g,%g,%g) REQj(%g,%g,%g) r %g e %g de %g \n",  dp.x,dp.y,dp.z,  REQj.x,REQj.y,REQj.z, r, e, de );
                    // --- store
                    qp.e+=eM*e; qp.f.add_mul( dp, de*e   ); // repulsive part of Morse
                    ql.e+=eM*2; ql.f.add_mul( dp, de     ); // attractive part of Morse
                    qe.e+=eQ;   qe.f.add_mul( dp, eQ*ir2 ); // Coulomb
                    //printf(  "evalMorsePLQ() k %g r %g e %g E0 %g E %g \n", K, r, e, REQj.y*plq.x/exp( K*(1.487)), qp.e*plq.x );
                    //if(i==0)printf( "[%i] qp(%g,%g,%g)  ql(%g,%g,%g)  qe(%g,%g,%g)\n", j, qp.x,qp.y,qp.z,   ql.x,ql.y,ql.z,   qe.x,qe.y,qe.z );
                }}}
            }
            Quat4f plq = PLQs[i];
            //if(i==0)printf( "plq[0](%g,%g,%g) qp(%g,%g,%g)  ql(%g,%g,%g)  qe(%g,%g,%g)\n", plq.x,plq.y,plq.z,  qp.x,qp.y,qp.z,   ql.x,ql.y,ql.z,   qe.x,qe.y,qe.z );
            Quat4d fe = qp*plq.x + ql*plq.y + qe*plq.z;
            fapos[i].add(fe.f);
            E       +=fe.e;
        }
        return E;
    }

    double evalR( NBFF& B ){
        //printf( "NBFF::evalR() \n" );
        const int n  = natoms;
        const int nb = B.natoms;
        double E=0;
        for(int i=0; i<n; i++){
            //Vec3d fi   = Vec3dZero;
            const Vec3d Api = apos[i];
            for(int j=0; j<nb; j++){    // atom-atom
                Vec3d fij = Vec3dZero;
                const Vec3d d = B.apos[j]-Api;
                double r = d.norm();
                E += fmax( 0, REQs[j].x-r );
                //fi.add(fij);
            }
            //fapos[i].add(fi);
        }
        return E;
    }
/*
    double evalNeighs( double RQ=1.0, double K=-1.5 ){
        //printf( "NBFF::evalNeighs() \n" );
        double   E=0;
        const double R2Q=RQ*RQ;
        const int    n=natoms;
        for(int i=0; i<n; i++){
            Vec3d  fi = Vec3dZero;
            const Vec3d  pi = apos[i];
            Quat4i ngi = {-1,-1,-1,-1};
            if(neighs)ngi=neighs[i];
            //pi.add( shift );
            const Quat4d& REQi = REQs[i];
            for(int j=i+1; j<n; j++){    // atom-atom
                if( (j==ngi.x)||(j==ngi.y)||(j==ngi.z)||(j==ngi.w) ) continue;
                Vec3d fij = Vec3dZero;
                Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
                const Vec3d dp=apos[j]-pi;
                //if(i==0){ idebug=1; }else{ idebug=0; };
                double ei = addAtomicForceMorseQ( dp, fij, REQij.x, REQij.y, REQij.z, K, R2Q );    E+=ei;
                //double ei = addAtomicForceQ_R2( dp, fij, REQij.z, K, R2Q );    E+=ei;
                //E += addAtomicForceLJQ   ( apos[j]-pi, fij, REQij );
                //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( fij      , pi );
                //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos( fij*-1.0f, apos[j] );
                //glColor3f(1.0f,0.0f,1.0f); Draw3D::drawLine( pi, apos[j] );
                //if(i==6){ printf("CPU[%i,%i] dp(%g,%g,%g) fe(%g,%g,%g|%g)\n", i,j,  dp.x,dp.y,dp.z,   fij.x,fij.y,fij.z,ei ); }
                fapos[j].sub(fij);
                fi     .add(fij);
            }
            fapos[i].add(fi);
        }
        return E;
    }
*/

    int makePBCshifts( Vec3i nPBC, const Mat3d& lvec ){
        npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
        printf( "NBFF::makePBCshifts() npbc=%i nPBC{%i,%i,%i}\n", npbc, nPBC.x,nPBC.y,nPBC.z );
        _realloc(shifts,npbc);
        int ipbc=0;
        for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ 
            for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ 
                for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
                    shifts[ipbc] = (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz);   
                    ipbc++; 
                }
            }
        }
        return npbc;
    }

    void print_nonbonded(){
        printf("NBFF::print_nonbonded(n=%i)\n", natoms );
        for(int i=0; i<natoms; i++){
            if(atypes){ printf("nb_atom[%i] REQ(%7.3f,%g,%g,%g) pos(%7.3f,%7.3f,%7.3f) atyp %i \n", i, REQs[i].x,REQs[i].y,REQs[i].z,REQs[i].w,   apos[i].x,apos[i].y,apos[i].z, atypes[i] ); }
            else      { printf("nb_atom[%i] REQ(%7.3f,%g,%g,%g) pos(%7.3f,%7.3f,%7.3f) \n",         i, REQs[i].x,REQs[i].y,REQs[i].z,REQs[i].w,   apos[i].x,apos[i].y,apos[i].z            ); }
        }
    }


    void checkREQlimits( const Quat4d vmin=Quat4d{ 0.2,0.0,-1.0,-1.0}, const Quat4d vmax=Quat4d{ 3.0,0.2,+1.0,+1.0}  ){
        if( REQs==0 ){ printf( "NBFF::checkREQlimits() REQs is NULL=> Exit()\n" ); exit(0); }
        if( checkLimits( natoms, 4, (double*)REQs, (double*)&vmin, (double*)&vmax, "REQs" ) ){
            printf("ERROR NBFF::checkREQlimits(): REQs are out of range (%g,%g,%g,%g) .. (%g,%g,%g,%g) => Exit() \n", vmin.x,vmin.y,vmin.z,vmin.w,  vmax.x,vmax.y,vmax.z,vmax.w ); print_nonbonded(); exit(0);
        }
    }

    void bindOrRealloc(int n_, Vec3d* apos_, Vec3d* fapos_, Quat4d* REQs_, int* atypes_ ){
        natoms=n_;
        //printf( "NBFF::bindOrRealloc(natoms=%i) @apos_=%li @fapos_=%li @REQs_=%li @atypes_=%li \n", natoms, (long)apos_, (long)fapos_, (long)REQs_, (long)atypes_ );
        _bindOrRealloc(natoms, apos_  , apos   );
        _bindOrRealloc(natoms, fapos_ , fapos  );
        _bindOrRealloc(natoms, REQs_  , REQs   );
        _bindOrRealloc(natoms, atypes_, atypes );
    }

    void dealloc(){
        natoms=0;
        _dealloc(apos   );
        _dealloc(fapos  );
        _dealloc(REQs   );
        _dealloc(atypes );
        _dealloc(neighs );
    }

};

#endif

