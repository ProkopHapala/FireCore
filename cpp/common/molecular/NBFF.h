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

#include "Forces.h"


bool checkLimits( int n, int m, const double* vals, const double* vmin, const double* vmax, const char* message, bool bPrint=true ){
    //for(int j=0; j<m; j++){ printf( "checkLimits[%i] [%g,%g]\n", j, vmin[j], vmax[j] ); }
    bool b=false;
    for(int i=0; i<n; i++){
        const double* vali = vals+i*m;
        for(int j=0; j<m; j++){
            double v = vali[j];
            if( v<vmin[j] || v>vmax[j] || isnan(v) ){  
                b=true; 
                if(bPrint){
                    printf( "%s[%i/%i,%i/%i] %g out of limits [%g,%g] \n", message, i,n, j,m,   v, vmin[j], vmax[j]  );
                }
            }
        }
    }
    return b;
}



class NBFF: public Atoms{ public:
    //int     n      =0; // from Atoms
    //int    *atypes =0; // from Atoms
    //Vec3d  *apos   =0; // from Atoms
    Vec3d    *fapos  =0;
    Quat4d   *REQs   =0;
    Quat4i   *neighs =0;
    Quat4i   *neighCell=0;

    double alphaMorse = 1.5;
    //double  KMorse  = 1.5;
    double  Rdamp     = 1.0;
    Mat3d   lvec;
    Vec3i   nPBC;
    bool    bPBC=false;

    int    npbc   =0;
    Vec3d* shifts =0;
    Quat4f *PLQs  =0;  // used only in combination with GridFF
    Vec3d  shift0 =Vec3dZero;

    // ==================== Functions

    void torq     ( const Vec3d& p0,  Vec3d& tq                   ){ for(int i=0; i<natoms; i++){ Vec3d d; d.set_sub(apos[i],p0); tq.add_cross(fapos[i],d); } }

    void bindShifts(int npbc_, Vec3d* shifts_ ){ npbc=npbc_; shifts=shifts_; }

    int makePBCshifts( Vec3i nPBC_, bool bRealloc=true ){
        bPBC=true;
        nPBC=nPBC_;
        npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
        if(bRealloc) _realloc(shifts,npbc);
        int ipbc=0;
        for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
            shifts[ipbc] = (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz);   
            ipbc++; 
        }}}
        if(npbc!=ipbc){ printf( "ERROR in MMFFsp3_loc::makePBCshifts() final ipbc(%i)!=nbpc(%i) => Exit()\n", ipbc,npbc ); exit(0); }
        return npbc;
    }

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
                Quat4d REQij; combineREQ( REQs[j], REQi, REQij );
                //E += addAtomicForceLJQ( apos[j]-pi, fij, REQij );
                E += getLJQH( apos[j]-pi, fij, REQij, R2damp );
                fapos[j].sub(fij);
                fi   .add(fij);
            }
            fapos[i].add(fi);
        }
        return E;
    }

    double evalLJQs_PBC( const Mat3d& lvec, Vec3i nPBC=Vec3i{1,1,1}, double Rdamp=1.0 ){
        //printf( "NBFF::evalLJQs_PBC() \n" );
        double R2damp = Rdamp*Rdamp;
        const int N   = natoms;
        double E=0;
        int npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1) -1;
        Vec3d shifts[npbc]; // temporary store for lattice shifts
        int ipbc=0;
        for(int ia=-nPBC.a; ia<(nPBC.a+1); ia++){ for(int ib=-nPBC.b; ib<(nPBC.b+1); ib++){ for(int ic=-nPBC.c; ic<(nPBC.c+1); ic++){ 
            if((ia==0)&&(ib==0)&&(ic==0))continue; // skipp pbc0
            shifts[ipbc] = (lvec.a*ia) + (lvec.b*ib) + (lvec.c*ic);   
            ipbc++; 
        }}}
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
            if( (ng.x==j)||(ng.y==j)||(ng.z==j)||(ng.w==j) ) continue;
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
                if( (ngs.x==j)||(ngs.y==j)||(ngs.z==j)||(ngs.w==j) ) continue;
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
                if(bBonded){
                    if(   ((j==ng.x)&&(ipbc==ngC.x))
                        ||((j==ng.y)&&(ipbc==ngC.y))
                        ||((j==ng.z)&&(ipbc==ngC.z))
                        ||((j==ng.w)&&(ipbc==ngC.w))
                    ){ continue;}
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
    double evalLJQs_ng4_PBC_simd(){
        //printf("NBFF::evalLJQs_ng4_PBC_simd()\n" );
        double E=0;
        for(int ia=0; ia<natoms; ia++){  
            //printf("ffls[%i].evalLJQs_ng4_PBC_atom_omp(%i)\n", id, ia ); 
            E+=evalLJQs_ng4_PBC_atom_omp(ia); 
        }
        return E;
    }

    double evalLJQs_ng4_atom_omp( const int ia ){
        //printf( "NBFF::evalLJQs_ng4_PBC_atom(%i)   apos %li REQs %li neighs %li neighCell %li \n", ia,  apos, REQs, neighs, neighCell );
        const double R2damp = Rdamp*Rdamp;
        const Vec3d  pi   = apos     [ia];
        const Quat4d  REQi = REQs     [ia];
        const Quat4i ng   = neighs   [ia];
        const Quat4i ngC  = neighCell[ia];
        Vec3d fi = Vec3dZero;
        double E=0,fx=0,fy=0,fz=0;

        #pragma omp simd reduction(+:E,fx,fy,fz)
        for (int j=0; j<natoms; j++){ 
            if( (ia==j)  || (j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w) ) continue;
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
    double evalLJQs_ng4_simd(){
        //printf("NBFF::evalLJQs_ng4_simd()\n" );
        double E=0;
        for(int ia=0; ia<natoms; ia++){ E+=evalLJQs_ng4_atom_omp(ia); }
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
                }
            }}} // nPBC
            //if(i==0){ printf( "CPU atom[%i]  fe_Cou(%g,%g,%g|%g)  REQKi.z %g \n", i, fi.x,fi.y,fi.z,E, REQi.z ); }
            fapos[i].add(fi);
        }
        return E;
    }

/*
    double evalMorsePLQ( NBFF& B, Mat3d& cell, Vec3i nPBC, double K=-1.0, double RQ=1.0 ){
        //printf( "NBFF::evalMorsePLQ() PLQs %li \n", (long)PLQs, K, RQ );
        double E=0;
        //printf( "NBFF nPBC(%i,%i,%i) K %g RQ %g R2Q %g plq.z %g \n", nPBC.x,nPBC.y,nPBC.z, K, sqrt(R2Q), R2Q, plq.z );
        double R2Q=RQ*RQ;
        for(int i=0; i<n; i++){
            Vec3d fi = Vec3dZero;
            Vec3d pi = apos[i];
            //const Quat4d& REQi = REQs[i];
            Quat4d qp = Quat4dZero;
            Quat4d ql = Quat4dZero;
            Quat4d qe = Quat4dZero;
            for(int j=0; j<B.n; j++){
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
*/

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


    void checkREQlimits( const Quat4d vmin=Quat4d{ 0.2,0.0,-1.,-1e-8}, const Quat4d vmax=Quat4d{ 3.0,0.2,+1.,+1e-8}  ){
        if( REQs==0 ){ printf( "NBFF::checkREQlimits() REQs is NULL=> Exit()\n" ); exit(0); }
        if( checkLimits( natoms, 4, (double*)REQs, (double*)&vmin, (double*)&vmax, "REQs" ) ){
            printf("ERROR NBFF::checkREQlimits(): REQs are out of range (%g,%g,%g,%g) .. (%g,%g,%g,%g) => Exit() \n", vmin.x,vmin.y,vmin.z,vmin.w,  vmax.x,vmax.y,vmax.z,vmax.w ); print_nonbonded(); exit(0);
        }
    }

    void bindOrRealloc(int n_, Vec3d* apos_, Vec3d* fapos_, Quat4d* REQs_, int* atypes_ ){
        natoms=n_;
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

