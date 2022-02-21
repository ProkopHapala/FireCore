#ifndef FireCore_h
#define FireCore_h

#include <dlfcn.h>

#include "MMFFparams.h"

namespace FireCore{

// subroutine firecore_evalForce( nmax_scf, forces_ )  bind(c, name='firecore_evalForce')
//  subroutine firecore_init( natoms_, atomTypes, atomsPos ) bind(c, name='firecore_init')
//typedef void (*Pprocedure)();
//typedef void (*Pfirecore_evalForce)(int,double*);

typedef void (*P_evalForce  )(int,double*,double*);
typedef void (*P_init       )(int,int*,double*);
//typedef void (*P_getCharges )(double*);
typedef void (*P_1d         )(double*);
typedef void (*P_proc       )();

class Lib{ public:
    void        *lib_handle = 0;
    P_evalForce  evalForce =0;
    P_init       init      =0;
    P_1d getCharges=0;
    P_1d set_lvs=0;

    P_proc preinit=0;

    //P_1d getp_ratom=0;
    //P_1d getp_ftot =0;


    inline int loadLib( const char* fullname ){
        char *error;
        if(lib_handle){ dlclose(lib_handle); lib_handle=0; };
        //char fullname[1024] = "/home/prokop/git/FireCore/build/libFireCore.so";
        printf("loading %s :\n", fullname );
        void* plib = dlopen( fullname, RTLD_LAZY | RTLD_GLOBAL );
        if (plib){
            lib_handle=plib;
        }else{ printf( "%s\n", dlerror()); return -1; }
        init      = (P_init)     dlsym(lib_handle, "firecore_init");
        if ((error = dlerror())){ printf( "%s\n", error); init    =0; exit(0); }
        evalForce = (P_evalForce)dlsym(lib_handle, "firecore_evalForce");
        if ((error = dlerror())){ printf("%s\n", error); evalForce=0; exit(0); }
        getCharges = (P_1d)dlsym(lib_handle, "firecore_getCharges");
        if ((error = dlerror())){ printf("%s\n", error); getCharges=0; exit(0); }
        set_lvs    = (P_1d)dlsym(lib_handle, "firecore_set_lvs");
        if ((error = dlerror())){ printf("%s\n", error); set_lvs=0; exit(0); }

        preinit    = (P_proc)dlsym(lib_handle, "firecore_preinit");
        if ((error = dlerror())){ printf("%s\n", error); preinit=0; exit(0); }

        //getp_ratom = (P_1d)dlsym(lib_handle, "firecore_getPointer_ratom");
        //!if ((error = dlerror())){ printf("%s\n", error); getp_ratom=0; exit(0); }
        //!getp_ftot = (P_1d)dlsym(lib_handle, "firecore_getPointer_ftot");
        //!if ((error = dlerror())){ printf("%s\n", error); getp_ftot=0; exit(0); }

        return 0;
    }

}; // class Lib

class QMMM{ public:
    int   nqm;      // number of atoms in QM
    int*  imms=0;   // map QM to those MM indexes
    bool* isCap=0;
    //int* iMMs=0;
    int*    atypeZ=0;
    int*    atype=0;
    Vec3d*  apos=0;
    Vec3d*  aforce=0;
    double* charges=0;
    MMFFparams* params=0;


    P_evalForce  p_evalForce=0;
    P_1d         p_getCharges=0;
    int nmax_scf = 100;

    void bindFireCoreLib(const Lib& lib){
        p_evalForce  = lib.evalForce;
        p_getCharges = lib.getCharges;
    }

    void init(int nqm_){
        nqm=nqm_;
        imms  = new int [nqm];
        isCap = new bool[nqm];
        atypeZ  = new int  [nqm];
        atype   = new int  [nqm];
        apos    = new Vec3d[nqm];
        aforce  = new Vec3d[nqm];
        charges = new double[nqm];
    }

    void setAtypes(int* atypes_){
        for(int i=0; i<nqm; i++){
            int im = imms[i];
            int Z=6;
            int ityp = atypes_[im];
            if(params)Z=params->atypes[ ityp ].iZ;
            atypeZ[i]= isCap[i] ? 1 : Z;   
            atype [i]= isCap[i] ? 0 : ityp;
            printf( "DEBUG atom %i, type %i -> iZ = %i \n", i, atypeZ[i] );
        }  // NOTE : This is just temporary hack
    }

    void qm2mm( int nj, const double* qm, double* mm)const{
        for(int iq=0; iq<nqm; iq++){
            int im = imms[iq];
            const double* qi=qm+iq*nj;
            double*       mi=mm+im*nj;
            for(int j=0; j<nj; j++){
                mi[j] = qi[j];
            }
        }
    }

    void mm2qm( int nj, double* qm, const double* mm)const{
        for(int iq=0; iq<nqm; iq++){
            int im = imms[iq];
            double*       qi=qm+iq*nj;
            const double* mi=mm+im*nj;
            for(int j=0; j<nj; j++){
                qi[j] = mi[j];
            }
            //for(int j=0; j<nj; j++){ qi[nj+j] = mi[j];  }
        }
    }

    void applyCharges( Vec3d* aREQs, bool bDiff=false ){
        p_getCharges( charges );
        if( bDiff ){
            if(params==0){ printf( "ERROR: QMMM::applyCharges() substracting valence charge but *params==null \n"); exit(0); }
            for(int iq=0; iq<nqm; iq++){ 
                const AtomType& at = params->atypes[ atype[iq] ];
                printf( "charges[%i|ityp=%i] Q %g Q0 %i iZ %i name %s\n", iq, atype[iq], charges[iq], at.neval, at.iZ, at.name  );
                charges[iq]-=params->atypes[ atype[iq] ].neval; 
            }
        }
        for(int iq=0; iq<nqm; iq++){
            int im      = imms[iq];
            if(params)
            aREQs[im].z = charges[iq];
        }

    }

    void load_apos  ( const Vec3d* mmpos   )     { 
        mm2qm(3,(      double*)apos  ,(const double*)mmpos   ); 
        //for(int i=0; i<nqm; i++){ int im=imms[i]; printf( "QM apos[%i->%i] cap %i Z %i %g %g %g |  %g %g %g \n",i,im, isCap[i], atypeZ[i], apos[i].x,apos[i].y,apos[i].z,   mmpos[im].x,mmpos[im].y,mmpos[im].z  ); }
    };
    void save_aforce(       Vec3d* mmforce )const{ qm2mm(3,(const double*)aforce,(      double*)mmforce ); };
    void evalQM(const Vec3d* mmpos, Vec3d* mmforce){
        load_apos( mmpos );
        p_evalForce ( nmax_scf, (double*)apos, (double*)aforce );
        save_aforce( mmforce );
    }

#ifdef MMFFmini_h
    void maskMMFF( MMFFmini& ff )const{
        std::unordered_set<int> imm_set;
        for(int i=0; i<nqm; i++){ imm_set.insert( imms[i] ); }
        for(int i:imm_set){  };
        for(int i=0; i<ff.nbonds; i++ ){
            printf( "DEBUG mask bond?  %i/%i \n", i, ff.nbonds );
            Vec2i iat = ff.bond2atom[i];
            if( imm_set.end() == imm_set.find(iat.a) ) continue;
            if( imm_set.end() == imm_set.find(iat.b) ) continue;
            ff.bondMasked[i] = true;
            printf( "DEBUG bond %i, (%i,%i) masked \n", i, iat.a, iat.b );
        }
        printf("DEBUG bondMasked DONE \n");
        for(int i=0; i<ff.nang; i++ ){
            Vec3i iat = ff.ang2atom[i];
            if( imm_set.end() == imm_set.find(iat.a) ) continue;
            if( imm_set.end() == imm_set.find(iat.b) ) continue;
            if( imm_set.end() == imm_set.find(iat.c) ) continue;
            ff.angMasked[i] = true;
            printf( "DEBUG angle %i, (%i,%i,%i) masked \n", i, iat.a, iat.b, iat.b );
        }
        printf("DEBUG angMasked DONE \n");
        for(int i=0; i<ff.nang; i++ ){
            Quat4i iat = ff.tors2atom[i];
            if( imm_set.end() == imm_set.find(iat.x) ) continue;
            if( imm_set.end() == imm_set.find(iat.y) ) continue;
            if( imm_set.end() == imm_set.find(iat.z) ) continue;
            if( imm_set.end() == imm_set.find(iat.w) ) continue;
            ff.torsMasked[i] = true;
            printf( "DEBUG torsion %i, (%i,%i,%i,%i) masked \n", i, iat.x, iat.y, iat.z, iat.w );
        }
        printf("DEBUG torsMasked DONE \n");
    }
#endif

}; // class QMMM

} // namespace FireCore

#endif
