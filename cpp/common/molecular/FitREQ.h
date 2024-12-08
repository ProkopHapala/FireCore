#ifndef FitREQ_h
#define FitREQ_h

#include <vector>
#include "Vec3.h"
#include "quaternion.h"
//#include "NBFF.h"
#include "Atoms.h"
#include "MMFFparams.h"
#include "Forces.h"

#include "MMFFBuilder.h"
#include "functions.h"

/**
 * Saves the coordinates of composite system comprising atoms stored in two systems A and B to a file in XYZ format. 
 * 
 * @param fname The name of the file to save the coordinates to.
 * @param A A pointer to the Atoms object containing the coordinates of the first set of atoms.
 * @param B A pointer to the Atoms object containing the coordinates of the second set of atoms.
 * @param comment An optional comment to include in the file.
 * @param mode The mode in which to open the file. Defaults to "w".
 */
void savexyz(const char* fname, Atoms* A, Atoms* B, const char* comment=0, const char* mode="w" ){ 
    FILE* fout = fopen( fname, mode);
    if(fout==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    int n=0;
    if(A)n+=A->natoms;
    if(B)n+=B->natoms;
    fprintf(fout,"%i\n", n);
    if(comment){ fprintf(fout,"%s\n", comment ); }else{ fprintf(fout, "comment\n", n); };
    if(A)A->atomsToXYZ(fout);
    if(B)B->atomsToXYZ(fout);
    fclose(fout);
}

/**
 * Applies a rigid transformation to an array of 3D points.
 * @param shift The translation vector to be applied to each point after the rotation.
 * @param unshift The optional vector to be subtracted from each point before the rotation (i.e. the center of rotation, and the inverse of the shift).
 * @param dir The direction vector for cosntructing the rotation matrix.
 * @param up The up vector for cosntructing the rotation matrix.
 * @param n The number of points in the input array.
 * @param pin The input array of 3D points.
 * @param pout The output array of 3D points after the transformation.
 */
void rigid_transform( Vec3d shift, Vec3d* unshift, Vec3d dir, Vec3d up, int n, Vec3d* pin, Vec3d* pout ){ 
    Mat3d M;
    M.c=dir;
    M.b=up;
    M.a.set_cross( dir, up );
    for(int i=0; i<n; i++){
        Vec3d p,p0=pin[i];
        if(unshift)p0.sub(*unshift);
        M.dot_to_T( p0, p );
        p.add(shift);
        pout[i]=p;
    }
}

inline double getEpairAtom( double r, double Hij, double k, double& dEdH ){
    dEdH = exp( -k * r );
    return Hij * dEdH;
}

struct AddedData{
    int    nep  =0;      // number of electron pairs
    Vec2i* bs   =0;      // bonds between electron pairs and host atoms  (host_atom_index, epair_index)
    Vec3d* dirs =0;      // directions of electron pairs

    // New members for H-bond optimization
    double  Emodel0     = NAN;     // Static model energy without H-bond corrections
    int     HBna        = 0;       // Number of atoms with significant H-bond corrections
    int     HBn0        = -1;      // Split point between fragment 1 and 2 in HB-filtered arrays
    int*    HBatomsInd  = 0;       // Index of host atom for electron pairs (-1 for non-electron-pair atoms)
    int*    HBatomsHost = 0;       // Index of host atom for electron pairs (-1 for non-electron-pair atoms)
    int*    HBatomsType = 0;       // Positions of H-bond corrected atoms
    Vec3d*  HBatomsPos  = 0;       // Positions of H-bond corrected atoms
    Quat4d* HBatomsREQs = 0;       // REQ parameters of H-bond corrected atoms

    AddedData() = default;
    AddedData(int natoms, int nhb) {
        realloc(natoms, nhb);
    }

    void realloc(int natoms, int nhb) {
        HBna = nhb;
        _realloc0( bs, natoms,   Vec2i{-1,-1} );
        _realloc0( dirs, natoms, Vec3dNAN );
        if(nhb==0) return;
        _realloc0( HBatomsInd,   nhb, -1 );
        _realloc0( HBatomsHost,  nhb, -1 );
        _realloc0( HBatomsType,  nhb, -1 );
        _realloc0( HBatomsPos,   nhb, Vec3dNAN );
        _realloc0( HBatomsREQs,  nhb, Quat4dNAN );
    }

    void dealloc() {
        _dealloc( bs );
        _dealloc( dirs );
        if(HBna==0) return;
        _dealloc( HBatomsInd );
        _dealloc( HBatomsHost );
        _dealloc( HBatomsType );
        _dealloc( HBatomsPos );
        _dealloc( HBatomsREQs );
    }

    ~AddedData() {
        dealloc();
    }

    int initFromREQs( int natoms, int n0, Quat4d* REQs, int* atypes, Vec3d* apos, bool bWrite=false ) {
        // Select only hydrogen bonded atoms
        int nfound  =  0;
        int nfound0 = -1;
        for(int ia=0; ia<natoms; ia++ ) {
            Quat4d REQ = REQs[ia];
            if( fabs(REQ.w) < 1e-300 ) continue;
            if( bWrite ) {
                HBatomsInd [nfound] = ia;
                HBatomsHost[nfound] = bs[ia].x;
                HBatomsType[nfound] = atypes[ia];
                HBatomsPos [nfound] = apos[ia];
                HBatomsREQs[nfound] = REQ;
            }
            if( ia < n0 ) nfound0 = nfound;
            nfound++;
        }
        HBna = nfound;
        HBn0 = nfound0+1;
        return nfound;
    }

};

// =================================
// ====   class    FitREQ       ====
// =================================

/**
 * @class FitREQ
 * @brief Class for fitting non-colvalent interaction parameters (like Lenard-Jonnes, Charge, Hydrogen Bond-correction) of a system of atoms to a set of training examples.
 * 
 * This class contains various arrays and variables for fitting parameters of a molecular system to a set of training examples.
 * It includes arrays for storing the parameters of each atom type, regularization parameters, and stiffness parameters.
 * It also includes arrays for storing the degrees of freedom (DOFs) of the system, as well as temporary arrays for accumulation of derivatives.
 * Additionally, it includes arrays for decomposition of energy components and arrays for storing the poses of the system.
 * 
 * The class provides methods for initializing the types of atoms, setting the system, setting rigid samples, and evaluating the derivatives of the system.
 */
class FitREQ{ public:
    
    //NBFF* nbff;
    int nDOFs=0,ntype=0,nbatch=0,n0=0,n1=0;
    int imodel=1;
    alignas(32) Quat4d*    typeREQs      =0; // [ntype] parameters for each type
    alignas(32) Quat4d*    typeREQsMin   =0; // [ntype] equlibirum value of parameters for regularization 
    alignas(32) Quat4d*    typeREQsMax   =0; // [ntype] equlibirum value of parameters for regularization 
    
    alignas(32) Quat4d*    typeREQs0     =0; // [ntype] equlibirum value of parameters for regularization
    alignas(32) Quat4d*    typeREQs0_low =0; // [ntype] equlibirum value of parameters for regularization (lower wall)
    alignas(32) Quat4d*    typeREQs0_high=0; // [ntype] equlibirum value of parameters for regularization (upper wall)

    alignas(32) Quat4d*    typeKreg      =0; // [ntype] regulatization stiffness
    alignas(32) Quat4d*    typeKreg_low  =0; // [ntype] regulatization stiffness (lower wall)
    alignas(32) Quat4d*    typeKreg_high =0; // [ntype] regulatization stiffness (upper wall)

    alignas(32) Quat4i*    typToREQ      =0; // [ntype] map each unique atom type to place in DOFs;
    // alignas(32)  Vec2i* REQtoTyp; // Maps DOF index to (type_index, component)
    // alignas(32)  Vec3d* DOFregX;  // regularization positions (xmin,x0,xmax) for each DOF
    // alignas(32)  Vec3d* DOFregK;  // regularization stiffness (Kmin,K0,Kmax) for each DOF

    std::vector<Vec2i> REQtoTyp; // Maps DOF index to (type_index, component)
    std::vector<Vec3d> DOFregX;  // regularization positions (xmin,x0,xmax) for each DOF
    std::vector<Vec3d> DOFregK;  // regularization stiffness (Kmin,K0,Kmax) for each DOF

    alignas(32) double*   DOFs =0;       // [nDOFs]
    alignas(32) double*   fDOFs=0;       // [nDOFs]
    alignas(32) double*   vDOFs=0;       // [nDOFs]
    alignas(32) double*   DOFs_old =0;   // [nDOFs]
    alignas(32) double*   fDOFs_old=0;   // [nDOFs]

    alignas(32) double*   sample_fdofs = 0; // [nDOFs*nsamples] - arrays for parallelization to avoid atomic-write conflicts

    bool  bEvalJ          = false;    // Should we evaluate variational derivatives on Fregment J 
    bool  bWriteJ         = false;    // Should we write variational derivatives to Fregment J ( inner loop over j )
    bool  bCheckRepulsion = false;    // Should we check maximum repulsion (EijMax) inside inner loop over j for each sample atoms ?
    bool  bRegularize     = true;     // Should we apply additional regularization forces to otimizer ( beside the true variational forces from inter-atomic forcefield ? )
    bool  bEpairs         = true;     // Should we add electron pairs to the molecule ?
    //bool  bOptEpR = false;          // Should we optimize electron pair distance (from host atom) ?
    bool bBroadcastFDOFs = false;

    // parameters
    double EmaxSample = 100.0; // maximum energy for sampling
    double Kneutral      = 1.0;
    //double max_step      = 0.1345;
    //double max_step      = 0.05; // maximal variation (in percentage) of any parameters allowed in one step during optimization
    double Rfac_tooClose = 0.0;
    double kMorse        = 1.6;

    double EijMax    = 5.0;
    int isamp_debug  = 0;
    int iBadFound    = 0; 
    int nBadFoundMax = 10;

    double*  weights = 0; // [nbatch] scaling importaince of parameters

    std::vector<Atoms*> samples;    // ToDo: would be more convenient to store Atoms* rather than Atoms

    MM::Builder builder;
    MMFFparams* params=0; 

    // New members for H-bond optimization
    bool bEvalOnlyCorrections = false;

// =================================
// =========== Functions ===========
// =================================

/**
 * @brief Reallocates the DOFs (degrees of freedom) array to the given size. affects the size of the fDOFs and vDOFs arrays as well.
 * 
 * @param nDOFs_ The new size of the DOFs array.
 */
void realloc( int nDOFs_ ){
    //nDOFs=nR+nE+nQ;
    nDOFs=nDOFs_;
    _realloc(  DOFs, nDOFs );
    _realloc( fDOFs, nDOFs );
    _realloc( vDOFs, nDOFs ); for(int i=0;i<nDOFs;i++){vDOFs[i]=0;}
    _realloc(  DOFs_old, nDOFs );
    _realloc( fDOFs_old, nDOFs );
    //Rs=DOFs;Es=DOFs+nR;Qs=DOFs+nR+nE;
    //fRs=fDOFs;fEs=fDOFs+nR;fQs=fDOFs+nR+nE;
}

/**
 * This function reallocates the memory for the fs array if the number of atoms in any batch is greater than the current maximum number of atoms.
 * It first determines the maximum number of atoms in any batch and then checks if it is greater than the current maximum number of atoms.
 * If it is, then it reallocates the memory for the fs array to accommodate the new maximum number of atoms.
 */
// void tryRealocSamp(){
//     int n=nmax;
//     for(int i=0; i<samples.size(); i++){ int ni = samples[i]->natoms; if(ni>n){ n=ni;} }
//     if(n>nmax){ _realloc( fs, n );  nmax=n; };
// }

void realloc_sample_fdofs(){
    _realloc0( sample_fdofs, nDOFs*samples.size(), 0.0 );
}

void reduce_sample_fdofs(){
    int nsamples = samples.size();
    for(int j=0; j<nDOFs; j++){ fDOFs[j] = 0.0; }
    for(int i=0; i<nsamples; i++){
        double* fs = sample_fdofs + i*nDOFs;
        for(int j=0; j<nDOFs; j++){ fDOFs[j] += fs[j]; }
    }
}

/**
 * Load a file of types involved in the parameter fitting. The file should contain lines with the following format:
 * atom_name mask_RvdW mask_EvdW mask_Q mask_Hb
 * where mask_RvdW, mask_EvdW, mask_Q, and mask_Hb are integers representing the indices of the fitting parameters for the RvdW, EvdW, Q, and Hb parameters, respectively.
 * @param fname The name of the file to load.
 * @return The number of types loaded.
*/
int loadTypeSelection( const char* fname ){
    if(verbosity>0)printf( "FitREQ::loadTypeSelection(fname=%s) \n", fname );
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    int ntypesel = 0;
    int nwExpected = 21;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        Quat4i tm;
        Quat4d tx0l, tx0h;
        Quat4d tkl, tkh;
        int nw = sscanf( line, "%s %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
        at_name, &tm.x,&tm.y,&tm.z,&tm.w, &tx0l.x,&tkl.x,&tx0h.x,&tkh.x, &tx0l.y,&tkl.y,&tx0h.y,&tkh.y, &tx0l.z,&tkl.z,&tx0h.z,&tkh.z, &tx0l.w,&tkl.w,&tx0h.w,&tkh.w ); // 20
        if(nw!=nwExpected){
            printf("ERROR: FitREQ::loadTypeSelection(fname=%s) line[%i] has %i words (expected %i) \n", fname, ntypesel, nw, nwExpected );
            printf("line[%i]: %s", ntypesel, line );
            exit(0);
        }
        ntypesel++;
    }
    fseek( fin, 0, SEEK_SET );
    std::vector<int> t;
    std::vector<Quat4i> typeMask;
    std::vector<Quat4d> ts0_high, ts0_low;
    std::vector<Quat4d> tk_high, tk_low;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        Quat4i tm;               // on/off regularization constrain {}   Quat because of REQH componets
        Quat4d tx0l, tx0h, tx0m; // position  of regularization constrain  low, high, medium
        Quat4d tkl,  tkh,  tkm;  // stiffness of regularization constrain  low, high, medium
        //printf( "line[%i]: %s", ntypesel, line );
        int nw = sscanf( line, "%s %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
        at_name, &tm.x,&tm.y,&tm.z,&tm.w, &tx0l.x,&tkl.x,&tx0h.x,&tkh.x, &tx0l.y,&tkl.y,&tx0h.y,&tkh.y, &tx0l.z,&tkl.z,&tx0h.z,&tkh.z, &tx0l.w,&tkl.w,&tx0h.w,&tkh.w ); // 20
        if(nw!=nwExpected){
            printf("ERROR: FitREQ::loadTypeSelection(fname=%s) line[%i] has %i words (expected %i) \n", fname, ntypesel, nw, nwExpected );
            printf("line[%i]: %s", ntypesel, line );
            exit(0);
        }
        //if(tm.z!=0){printf("ERROR: FitREQ::loadTypeSelection() mask_Q should be 0 for now\n"); exit(0); }
        typeMask.push_back( tm );
        ts0_low .push_back( tx0l );
        ts0_high.push_back( tx0h );
        tk_low  .push_back( tkl );
        tk_high .push_back( tkh );
        int ityp = params->getAtomType(at_name);
        t.push_back( ityp );
        AtomType& t  = params->atypes[ityp];  
        if(verbosity>0)printf( "ityp(%i): %s tm(%i,%i,%i,%i)   R(xl=%lf,Kl=%lf|xh=%lf,Kh=%lf) E(xl=%lf,Kl=%lf|xh=%lf,Kh=%lf)  Q(xl=%lf,Kl=%lf|xh=%lf,Kh=%lf)  H(xl=%lf,Kl=%lf|xh=%lf,Kh=%lf)  \n", 
        ityp, at_name,  tm.x,tm.y,tm.z,tm.w,    tx0l.x,tkl.x,tx0h.x,tkh.x,   tx0l.y,tkl.y,tx0h.y,tkh.y,    tx0l.z,tkl.z,tx0h.z,tkh.z,    tx0l.w,tkl.w,tx0h.w,tkh.w );
        //tREQs.push_back( Quat4d{ t.RvdW, t.EvdW, t.Qbase, t.Hb } );
    }
    fclose(fin);
    init_types( ntypesel, &t[0], &typeMask[0], &ts0_low[0], &tk_low[0], &ts0_high[0], &tk_high[0]);
    return ntypesel;
}

int loadWeights( const char* fname ){
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    int n = 0;
    while( fgets(line, nline, fin) ){
        double w;
        sscanf( line, "%g\n", &w ); 
        n++;
    }
    fseek( fin, 0, SEEK_SET );
    _realloc( weights, n );
    int i = 0;
    while( fgets(line, nline, fin) ){
        double w;
        sscanf( line, "%g\n", &w ); 
        weights[i]=w;
        i++;
    }
    fclose(fin);
    return n;
}

void initDOFregs(){
    printf( "FitREQ::initDOFregs() nDOFs=%i \n", nDOFs );
    DOFregX.resize(nDOFs);
    DOFregK.resize(nDOFs);
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = REQtoTyp[i];        // which type and which component (REQH)
        DOFregX[i] = Vec3d{
            typeREQs0_low [rt.x].array[rt.y],  // xmin
            typeREQs0     [rt.x].array[rt.y],  // x0
            typeREQs0_high[rt.x].array[rt.y]   // xmax
        };
        DOFregK[i] = Vec3d{
            typeKreg_low  [rt.x].array[rt.y],  // Kmin
            typeKreg      [rt.x].array[rt.y],  // K0
            typeKreg_high [rt.x].array[rt.y]   // Kmax
        };
    }
}

/**
 * Initializes the types of the FitREQ object. This function calculates the number of degrees of freedom (nDOFs) and initializes the typToREQ array, which maps each atom type to its corresponding REQ values.
 * @param ntype_ The number of types.
 * @param ntypesel The number of selected types.
 * @param tsel An array of integers representing the selected types.
 * @param typeMask An array of Quat4i indicating which of the 4 parameters (Rvdw,Evdw,Q,Hb) are free to be fitted.
 * @param typeREQs An array of Quat4d objects representing the non-colvalent interaction parameters (Rvdw,Evdw,Q,Hb) for each type.
 * @param typeREQsMin An array of Quat4d objects representing the minimum values of the non-colvalent interaction parameters (Rvdw,Evdw,Q,Hb) for each type.
 * @param typeREQsMax An array of Quat4d objects representing the maximum values of the non-colvalent interaction parameters (Rvdw,Evdw,Q,Hb) for each type.
 * @param typeREQs0 An array of Quat4d objects representing the equilibrium values of the non-colvalent interaction parameters (Rvdw,Evdw,Q,Hb) for each type.
 * @param typeKreg An array of Quat4d objects representing the regularization stiffness for each type.
 * @return The number of degrees of freedom.
 */
int init_types( int ntypesel, int* tsel, Quat4i* typeMask, Quat4d* ts0_low, Quat4d* tk_low, Quat4d* ts0_high, Quat4d* tk_high ){
    int ntype_ = params->atypes.size();
    //printf( "FitREQ::init_types() ntypesel=%i ntypes=%i \n", ntypesel, ntype_ );
    int nDOFs=0;
    typToREQ = new Quat4i[ntype_];
    for(int i=0; i<ntype_; i++){
        bool bFound=false;
        for(int j=0; j<ntypesel; j++){
            if(tsel[j]==i){
                const Quat4i& tm=typeMask[j];
                Quat4i&       tt=typToREQ[i];
                for(int k=0; k<4; k++){
                    if(tm.array[k]){
                        tt.array[k]=nDOFs;
                        nDOFs++;
                    }else{
                        tt.array[k]=-1;
                    }
                }
                //printf(  "init_types() [%i] typeMask(%i,%i,%i,%i) typToREQ(%i,%i,%i,%i) nDOF %i\n", i, tm.x,tm.y,tm.z,tm.w, tt.x,tt.y,tt.z,tt.w, nDOFs );
                bFound=true;
                break;
            }
        }
        if(!bFound){
            Quat4i&       tt=typToREQ[i];
            for(int j=0; j<4; j++){
                tt.array[j]=-1;
            }
        }
    }
    realloc(nDOFs);
    ntype = ntype_; 
    //typeREQs       = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs      [i] = Quat4d{ params->atypes[i].RvdW, params->atypes[i].EvdW, params->atypes[i].Qbase, params->atypes[i].Hb }; }
    typeREQs       = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs      [i] = Quat4d{ params->atypes[i].RvdW, sqrt(params->atypes[i].EvdW), params->atypes[i].Qbase, params->atypes[i].Hb }; }  // We do square root of EvdW, no need to do it in evalExampleDerivs_*
    typeREQsMin    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin   [i] = Quat4d{ -1e300, -1e300, -1e300, -1e300 }; }
    typeREQsMax    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax   [i] = Quat4d{ 1e+300, 1e+300, 1e+300, 1e+300 }; }
    typeREQs0_low  = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0_low [i] = Quat4dZero;  }
    typeREQs0      = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0     [i] = typeREQs[i]; } // Initialize equilibrium point to current (initial) values
    typeREQs0_high = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0_high[i] = Quat4dZero;  }
    typeKreg_low   = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg_low  [i] = Quat4dZero;  }
    typeKreg       = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg      [i] = Quat4dZero;  }
    typeKreg_high  = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg_high [i] = Quat4dZero;  }
    for(int i=0; i<ntype; i++){ 
        for(int j=0; j<ntypesel; j++){
            if(tsel[j]==i){
                typeREQs0_low [i] = ts0_low [j];
                typeKreg_low  [i] = tk_low  [j];
                typeREQs0_high[i] = ts0_high[j];
                typeKreg_high [i] = tk_high [j];
                break;
            }
        }
    }
    //DOFsFromTypes(); 
    typesToDOFs();
    // Build inverse mapping from DOFs to types
    REQtoTyp.resize(nDOFs);
    for(int i=0; i<ntype; i++){
        const Quat4i& tt = typToREQ[i];
        if(tt.x>=0) REQtoTyp[tt.x] = Vec2i{i,0};
        if(tt.y>=0) REQtoTyp[tt.y] = Vec2i{i,1};
        if(tt.z>=0) REQtoTyp[tt.z] = Vec2i{i,2};
        if(tt.w>=0) REQtoTyp[tt.w] = Vec2i{i,3};
    }
    initDOFregs();
    if(verbosity>0){
        printDOFsToTypes();
        printTypesToDOFs();
        printDOFregularization();
    }    
    return nDOFs;
}


// add electron pairs
Atoms* addEpairs( Atoms* mol ){
    //MM::Builder builder;
    builder.params = params;
    builder.clear();
    if(mol->lvec){ builder.lvec = *(mol->lvec); builder.bPBC=true; }
    builder.insertAtoms(*mol);
    //int ia0=builder.frags[ifrag].atomRange.a;
    //int ic0=builder.frags[ifrag].confRange.a;
    //builder.printBonds();
    //builder.printAtomConfs(true, false );
    //if(iret<0){ printf("!!! exit(0) in MolWorld_sp3::loadGeom(%s)\n", name); exit(0); }
    //builder.addCappingTypesByIz(1);  // insert H caps
    //for( int it : builder.capping_types ){ printf( "capping_type[%i] iZ=%i name=`%s`  \n", it, builder.params->atypes[it].iZ, builder.params->atypes[it].name ); };
    builder.tryAddConfsToAtoms( 0, -1 );
    builder.cleanPis();
    //if(verbosity>2)
    //builder.printAtomConfs(false);
    //builder.export_atypes(atypes);
    // ------- Load lattice vectros
    // NOTE: ERROR IS HERE:  autoBonds() -> insertBond() -> tryAddBondToAtomConf( int ib, int ia, bool bCheck )
    //       probably we try to add multiple bonds for hydrogen ?
    if( builder.bPBC ){ 
        builder.autoBondsPBC( -.5,  0      , mol->n0     ); // we should not add bonds between rigid and flexible parts
        builder.autoBondsPBC( -.5,  mol->n0, mol->natoms ); 
    }else{ 
        builder.autoBonds( -.5,  0      , mol->n0     ); // we should not add bonds between rigid and flexible parts
        builder.autoBonds( -.5,  mol->n0, mol->natoms ); 
    }
    builder.checkNumberOfBonds( true, true );
    //if(verbosity>2)
    //builder.printBonds ();
    //if( fAutoCharges>0 )builder.chargeByNeighbors( true, fAutoCharges, 10, 0.5 );
    //if(substitute_name) substituteMolecule( substitute_name, isubs, Vec3dZ );
    //if( builder.checkNeighsRepeat( true ) ){ printf( "ERROR: some atoms has repating neighbors => exit() \n"); exit(0); };
    // --- Add Epairs
    builder.bDummyEpair = true;
    builder.autoAllConfEPi( ); 
    // --- add Epairs to atoms like =O (i.e. with 2 just one neighbor)
    builder.setPiLoop       ( 0, -1, 10 );
    builder.addAllEpairsByPi( 0, -1 );    
    //builder.printAtomConfs( false, true );
    //builder.assignAllBondParams();    //if(verbosity>1)
    //builder.finishFragment(ifrag);    
    return builder.exportAtoms(); 
}

/**
 * @brief Process a molecular system by adding electron pairs and reordering atoms to maintain a specific structure
 * 
 * This function performs the following steps:
 * 1. Adds electron pairs to the molecular system using the builder
 * 2. Identifies root atoms and directions for electron pairs
 * 3. Reorders atoms to maintain the structure: Atoms(mol1), Epairs(mol1), Atoms(mol2), Epairs(mol2)
 * 4. Stores the electron pair data in the atoms' userData
 * 
 * @param atoms Input molecular system to process
 * @param fout Optional file pointer to write XYZ output (can be null)
 * @return void, but modifies atoms in place
 */
void addAndReorderEpairs(Atoms*& atoms, FILE* fout=nullptr) {
    //printf( "addAndReorderEpairs()\n" );
    // Store original system information
    Atoms* bak = atoms;
    int n0bak  = bak->n0;        // Number of atoms in first molecule
    int natbak = bak->natoms;    // Total number of atoms before adding epairs
    int n1bak  = natbak - n0bak; // Number of atoms in second molecule

    // Add electron pairs to the system
    atoms = addEpairs(atoms);      // This modifies the builder's state
    atoms->n0 = bak->n0;
    atoms->Energy = bak->Energy;
    for(int i=0; i<natbak; i++){ atoms->charge[i] = bak->charge[i]; }
    delete bak;

    // Get electron pair information from builder
    Vec2i* bs   =nullptr;    // electron pair bond { ia= Host_atom_index,  ja=Electron_Pair_Index }
    Vec3d* dirs =nullptr;    // direction of electron pair bond  (Builder will allocate bs,dirs internally )
    int nep_found = builder.listEpairBonds(bs, dirs);  // This uses builder's state from addEpairs
    if(nep_found <= 0) { printf("Warning: No electron pairs found\n"); return;  }

    // Store original bonds and directions before any reordering
    std::vector<Vec2i> bsbak;
    std::vector<Vec3d> dirsbak;
    bsbak  .reserve(nep_found);
    dirsbak.reserve(nep_found);
    for(int i=0; i < nep_found; i++){
        dirs[i].normalize();
        bsbak  .push_back(bs  [i]);
        dirsbak.push_back(dirs[i]);
    }

    // Initialize array marking which atoms are electron pairs
    int* isep = new int[atoms->natoms];
    for(int i=0; i < natbak; i++){ isep[i]=0; }         // Original atoms
    for(int i=natbak; i < atoms->natoms; i++){ isep[i]=1; } // New epairs

    // Make temporary copy for reordering
    Atoms bak2;
    bak2.copyOf(*atoms);

    // Process electron pairs for first molecule (mol1)
    int nE0 = 0;  // Number of epairs in first molecule
    for(int i=0; i<nep_found; i++) {
        if(bsbak[i].x < n0bak) {  // Check if root atom is in first molecule
            int j = n0bak + nE0;   // Target position for this epair
            int k = bsbak[i].y;    // Source position of this epair
            
            if(j < atoms->natoms && k < bak2.natoms) {  // Bounds check
                // Copy epair data to its new position
                atoms->apos[j]   = bak2.apos[k];
                atoms->atypes[j] = bak2.atypes[k];
                atoms->charge[j] = bak2.charge[k];
                atoms->n0++;
                
                // Update bond and direction information
                if(nE0 < nep_found) {  // Bounds check
                    bs[nE0].x = bsbak[i].x;  // Root atom index stays same
                    bs[nE0].y = j;           // Update epair's new position
                    dirs[nE0] = dirsbak[i];
                    isep[j] = 1;
                    nE0++;
                }
            }
        }
    }

    int nE1 = atoms->natoms - natbak - nE0;  // Number of epairs in second molecule

    // Copy atoms from second molecule (mol2)
    for(int i=0; i<n1bak && (atoms->n0 + i) < atoms->natoms; i++) {
        int j = atoms->n0 + i;     // Target position
        int k = n0bak + i;         // Source position
        
        if(k < bak2.natoms) {  // Bounds check
            atoms->apos[j]   = bak2.apos[k];
            atoms->atypes[j] = bak2.atypes[k];
            atoms->charge[j] = bak2.charge[k];
            isep[j] = 0;
        }
    }

    // Process electron pairs for second molecule (mol2)
    int iE = -1;
    for(int i=0; i<nep_found; i++) {
        if(bsbak[i].x >= n0bak) {  // Check if root atom is in second molecule
            iE++;
            int j = n0bak + nE0 + n1bak + iE;  // Target position for this epair
            int k = bsbak[i].y;                // Source position of this epair
            
            if(j < atoms->natoms && k < bak2.natoms && (nE0+iE) < nep_found) {  // Bounds check
                // Copy epair data to its new position
                atoms->apos[j]   = bak2.apos[k];
                atoms->atypes[j] = bak2.atypes[k];
                atoms->charge[j] = bak2.charge[k];
                
                // Update bond and direction information
                bs  [nE0+iE].x = bsbak[i].x + nE0;  // Adjust root atom index
                bs  [nE0+iE].y = j;
                dirs[nE0+iE]   = dirsbak[i];
                isep[j] = 1;
            }
        }
    }

    // Store electron pair relationships in AddedData
    AddedData* data = new AddedData();
    data->nep       = nep_found;
    data->bs        = bs;       //  bs[i].x is epair index, bs[i].y is host atom index
    data->dirs      = dirs;   
    atoms->userData = data;


    // Write output if requested
    if(fout) {
        char line[1024];
        sprintf(line, "#	n0 %i E_tot %g", atoms->n0, atoms->Energy);
        params->writeXYZ(fout, atoms, line, 0, true);
    }
    //printf( "addAndReorderEpairs() DONE \n" );
}

/**
 * Loads XYZ file and creates Atoms objects for each molecule inside it. It saves the molecules to the "samples" vector.
 * Optionally adds electron pairs and outputs XYZ file with epairs.
 *
 * @param fname The name of the XYZ file to load.
 * @param bAddEpairs Flag indicating whether to add epairs to the loaded atoms.
 * @param bOutXYZ Flag indicating whether to output XYZ file with epairs.
 * @return The number of batches created.
 */
int loadXYZ( const char* fname, bool bAddEpairs=false, bool bOutXYZ=false ){
    printf( "FitREQ::loadXYZ()\n" );
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    // --- Open output file
    FILE* fout=0;
    if(bAddEpairs && bOutXYZ){
        sprintf(line,"%s_Epairs.xyz", fname );
        fout = fopen(line,"w");
    }
    int il   = 0;
    nbatch=0;
    Atoms* atoms=0;
    while( fgets(line, nline, fin) ){
        if      ( il==0 ){               // --- Read number of atoms
            int na=-1;
            sscanf( line, "%i", &na );
            if( (na>0)&&(na<10000) ){
                atoms = new Atoms(na);
            }else{ printf( "ERROR in FitREQ::loadXYZ() Suspicious number of atoms (%i) while reading `%s`  => Exit() \n", na, fname ); exit(0); }
        }else if( il==1 ){               // --- Read comment line ( read reference energy )
            sscanf( line, "%*s %*s %i %*s %lf ", &(atoms->n0), &(atoms->Energy) );
        }else if( il<atoms->natoms+2 ){  // --- Road atom line (type, position, charge)
            double x,y,z,q;
            int nret = sscanf( line, "%s %lf %lf %lf %lf", at_name, &x, &y, &z, &q );
            if(nret<5){q=0;}
            int i=il-2;
            atoms->apos[i].set(x,y,z);
            atoms->atypes[i]=params->getAtomType(at_name);
            atoms->charge[i]=q;
        }
        il++;
        if( il >= atoms->natoms+2 ){    // ---- store sample atoms to batch
            if(bAddEpairs){ addAndReorderEpairs(atoms, fout); }
            samples.push_back( atoms );
            il=0; nbatch++;
        }
    }
    if(fout)fclose(fout);
    fclose(fin);
    //init_types();
    return nbatch;
}

int export_Erefs( double* Erefs ){
    if(Erefs){ for(int i=0; i<nbatch; i++){ Erefs[i] = samples[i]->Energy; } }
    return nbatch;
}

void DOFsToTypes(){
    //printf( "DOFsToTypes() \n" );
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = REQtoTyp[i];
        //printf("%i -> %i|%i  = %s.%c \n", i, rt.x, rt.y,  params->atypes[rt.x].name, comp);
        //printf( "DOFsToType()[%3i] rt(%3i|%i) %8s.%c  %g \n", i, rt.x,rt.y,   params->atypes[rt.x].name,"REQH"[rt.y], DOFs[i] );
        typeREQs[rt.x].array[rt.y] = DOFs[i];
    }
}

void DOFsFromTypes(){
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = REQtoTyp[i];
        DOFs[i] = typeREQs[rt.x].array[rt.y];
    }
}

/**
 * @brief reads non-colvalent interaction parameter REQH(Rvdw,Evdw,Q,Hb) of given atom-type from aproprieate degrees of freedom (DOF) according to index stored in typToREQ[ityp]. If index of DOF is negative, the parameter is not read.  
 * @param ityp Index of the type
 */
inline void typeFromDOFs(int ityp){
    const Quat4i& tt = typToREQ[ityp];
    Quat4d& REQ      = typeREQs[ityp];
    if(tt.x>=0)REQ.x = DOFs[tt.x];
    if(tt.y>=0)REQ.y = DOFs[tt.y];
    if(tt.z>=0)REQ.z = DOFs[tt.z];
    if(tt.w>=0)REQ.w = DOFs[tt.w];
    //printf( "DOFsToType(%i) %g %g %g %g \n", ityp,   tt.x,tt.y,tt.z,tt.w,  REQ.x,REQ.y,REQ.z,REQ.w ); // Debug
}
void typesFromDOFs(){ for(int i=0; i<ntype; i++ ){ typeFromDOFs(i); } }
void getTypeFromDOFs(int i, Quat4d& REQ ){ typeREQs[i]=REQ; typeFromDOFs(i); }

/**
 * @brief writes non-colvalent interaction parameter REQH(Rvdw,Evdw,Q,Hb) of given atom-type from aproprieate degrees of freedom (DOF) according to index stored in typToREQ[ityp]. If index of DOF is negative, the parameter is not writen.  
 * 
 * @param ityp The index of the type.
 */
inline void typeToDOFs(int ityp){
    const Quat4i& tt  = typToREQ[ityp];
    const Quat4d& REQ = typeREQs[ityp];
    if(tt.x>=0)DOFs[tt.x] = REQ.x;
    if(tt.y>=0)DOFs[tt.y] = REQ.y;
    if(tt.z>=0)DOFs[tt.z] = REQ.z;
    if(tt.w>=0)DOFs[tt.w] = REQ.w;
}
void typesToDOFs(){ for(int i=0; i<ntype; i++ ){ typeToDOFs(i); } }
void setTypeToDOFs(int i, Quat4d REQ ){ typeREQs[i]=REQ; typeToDOFs(i); }

// ======================================
// =========  EVAL DERIVS  ==============
// ======================================

//void clean_fs(int n){ for(int i=0; i<n; i++){fs[i]=Quat4dZero;} }

__attribute__((hot)) 
void fillTempArrays( const Atoms* atoms, Vec3d* apos, double* Qs  )const{
    //printf( "FillTempArrays() bEpairs=%i \n", bEpairs );
    for(int j=0; j<atoms->natoms; j++){
        Qs[j]   = atoms->charge[j];
        apos[j] = atoms->apos [j];
        //isEp[j] = 0;
    }
    if( bEpairs ){  // electron pairs
        AddedData * ad = (AddedData*)atoms->userData;
        int   nep = ad->nep;
        Vec2i* bs = ad->bs;
        //printf( "FillTempArrays() 2  ad->nep=%i \n", ad->nep );
        for(int j=0; j<ad->nep; j++){
            int    iX  = bs[j].x; // index of the host atom
            int    iE  = bs[j].y; // index of the electron pair
            Quat4d REQ = typeREQs[atoms->atypes[iE]];
            //printf( "FillTempArrays()[iap=%i] iE=%i iX=%i  name=%8s REQ(%12.3e,%12.3e,%12.3e,%12.3e) \n", j, iE, iX, params->atypes[atoms->atypes[iE]].name, REQ.x,REQ.y,REQ.z,REQ.w );
            double Qep = REQ.z; // charge of the electron pair
            Qs[iE]     = Qep;
            Qs[iX]    -= Qep;
            apos[iE] = apos[iX] + ad->dirs[j] * typeREQs[atoms->atypes[iE]].w;   // We move the electron pair to proper distance from the atom
            //isEp[iE] = 1;
        }
    }
}

__attribute__((hot)) 
double evalSample( int isamp, const Atoms* atoms, double wi, Quat4d* fs ) const {
    //double wi   = (weights)? weights[isamp] : 1.0; 
    alignas(32) double Qs  [atoms->natoms];
    alignas(32) Vec3d  apos[atoms->natoms];   // atomic positions
    fillTempArrays( atoms, apos, Qs );
    int     nj = atoms->n0;
    int     j0 = 0; 
    int     ni = atoms->natoms - atoms->n0;
    int     i0 = atoms->n0;
    for(int i=0; i<atoms->natoms; i++){ fs[i]=Quat4dZero; }
    double E=0;
    bool bJ = bEvalJ && ( !bWriteJ );
    switch (imodel){
        case 0:{ 
            E =   evalExampleDerivs( funcVar_LJQH2, i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs( funcVar_LJQH2, j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fs );    // variational derivatives on molecule 2
        }break;
        case 1:{ 
            E =   evalExampleDerivs_LJQH2   ( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs_LJQH2   ( j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fs );    // variational derivatives on molecule 2
        }break;
        case 2:{ 
            E =   evalExampleDerivs_MorseQH2( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs_MorseQH2( j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fs );    // variational derivatives on molecule 2
        }break;
    }
    if( E>EmaxSample ){
        if(verbosity>0) printf( "skipped sample [%i] E(%g)>EmaxSample(%g) atoms too close \n", isamp, E, EmaxSample );
        return E;
    } 
    return E;
}

double evalSample_corr( int isamp, const AddedData* adata, double wi, Quat4d* fs ) const {
    int     nj = adata->HBn0;
    int     j0 = 0; 
    int     ni = adata->HBna - adata->HBn0;
    int     i0 = adata->HBn0;
    for(int i=0; i<adata->HBna; i++){ fs[i]=Quat4dZero; }
    double E=0;
    bool bJ = bEvalJ && ( !bWriteJ );
    switch (imodel){
        case 1:{ 
            // double evalExampleDerivs_LJQH2_corr( int i0, int ni, int j0, int nj, Vec3d* ps, Quat4d* REQs, int* host, Quat4d* dEdREQs ) const {
            E =   evalExampleDerivs_LJQH2_corr   ( i0, ni, j0, nj, adata->HBatomsPos, adata->HBatomsREQs, adata->HBatomsHost, fs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs_LJQH2_corr   ( j0, nj, i0, ni, adata->HBatomsPos, adata->HBatomsREQs, adata->HBatomsHost, fs );    // variational derivatives on molecule 2
        }break;
    }
    if( E>EmaxSample ){
        if(verbosity>0) printf( "skipped sample [%i] E(%g)>EmaxSample(%g) atoms too close \n", isamp, E, EmaxSample );
        return E;
    } 
    return E;
}

__attribute__((hot)) 
double evalSampleError( int isamp, double& E ){
    //isamp_debug = i;
    Atoms* atoms  = samples[isamp];
    double wi     = (weights)? weights[isamp] : 1.0;
    if(wi<1e-300) return 0;
    alignas(32) Quat4d fs  [atoms->natoms];
    E = evalSample( isamp, atoms, wi, fs );
    //printf( "evalDerivsSamp() isamp: %3i E: %20.10f Eref: %20.10f \n", i, E, Eref );
    double Eref    = atoms->Energy;
    double dE      = E - Eref;
    double dEw     = 2.0*dE*wi;
    double Error   = dE*dE*wi;
    double* fDOFs_ = fDOFs;
    if(bBroadcastFDOFs){ 
        fDOFs_ = sample_fdofs + isamp*nDOFs; 
        for(int k=0; k<nDOFs; k++){ fDOFs_[k]=0; }
    }
    acumDerivs( atoms->natoms, atoms->atypes, dEw, fs, fDOFs_ );
    if( bEpairs ){
        AddedData * ad = (AddedData*)atoms->userData;
        acumHostDerivs( ad->nep, ad->bs, atoms->atypes, dEw, fs, fDOFs_ );
    }
    return Error;
}

__attribute__((hot)) 
double evalSampleError_corr( int isamp, double& E ){
    //isamp_debug = i;
    Atoms* atoms  = samples[isamp];
    const AddedData* adata = (const AddedData*)(atoms->userData);
    double wi     = (weights)? weights[isamp] : 1.0;
    if(wi<1e-300) return 0;
    alignas(32) Quat4d fs  [adata->HBna];
    E = evalSample_corr( isamp, adata, wi, fs );
    //printf( "evalDerivsSamp() isamp: %3i E: %20.10f Eref: %20.10f \n", i, E, Eref );
    double Emodel0 = adata->Emodel0;
    double Eref    = atoms->Energy;
    double dE      = E - Emodel0 - Eref;
    double dEw     = 2.0*dE*wi;
    double Error   = dE*dE*wi;
    double* fDOFs_ = fDOFs;
    if(bBroadcastFDOFs){ 
        fDOFs_ = sample_fdofs + isamp*nDOFs; 
        for(int k=0; k<nDOFs; k++){ fDOFs_[k]=0; }
    }
    acumDerivs_corr( adata->HBna, adata->HBatomsType, adata->HBatomsHost, dEw, fs, fDOFs_ );
    return Error;
}

// =========== Correction Only optimized versions of evalSample and evalSampleError



/**
 * @brief Evaluates the variational derivatives of the fitting error (sumed batch training samples) with respect to all fitting parameters (i.e. non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) for each atom type). The atomic systems are not assumed rigid (i.e. the atomic positions can vary freely from sample to sample).
 * 
 * @param Eout array to store the non-covalent interaction energy values of each atomic system in the batch. if Eout==null, the function will not store the energy values.
 * @return double, returns the total fitting error.
 */
__attribute__((hot)) 
double evalSamples_noOmp( double* Eout=0 ){ 
    //printf( "evalSamples_noOmp() \n" );
    double Error = 0.0;
    int nsamp = samples.size();
    for(int i=0; i<nsamp; i++){
       double E; 
       Error+=evalSampleError( i, E );
       if(bEvalOnlyCorrections){ ((AddedData*)samples[i]->userData)->Emodel0 = E; }
       if(Eout)Eout[i]=E;
    }
    return Error;
}

__attribute__((hot)) 
double evalSamples_omp( double* Eout=0 ){ 
    //printf( "evalSamples_omp() \n" );
    double Error = 0.0;
    int nsamp = samples.size();
    #pragma omp parallel for reduction(+:Error)
    for(int i=0; i<nsamp; i++){
       double E; 
       Error+=evalSampleError( i, E );
       if(Eout)Eout[i]=E;
    }
    return Error;
}

void printDebugArrays(int na, int* atypes, Vec3d* apos, double* aq, int* aisep, int nj, int* jtyp, Vec3d* jpos, double* jq, int* jisep, int nep, Vec2i* bs, const Atoms* atoms, bool bEpairs) {
    //printf("=== DEBUG: Arrays in evalDerivsSamp ===\n");
    //printf("Number of atoms in first molecule (nj):  %d\n", nj);
    //printf("Number of atoms in second molecule (na): %d\n", na);
    //printf("Number of electron pairs (nep): %d\n", nep);
    AddedData * ad = (AddedData*)atoms->userData;
    printf("idx  type  charge  isep   position(x,y,z)  \n");
    printf("Molecule  #1:\n"); 
    for(int j=0; j<nj; j++){  printf("[%-3d] type %-3d %8s  Q=%7.3f ep?=%i  apos(%8.3f,%8.3f,%8.3f)\n",  j, jtyp[j], params->atypes[jtyp[j]].name, jq[j], jisep[j], jpos[j].x, jpos[j].y, jpos[j].z ); }
    printf("Molecule  #2:\n"); 
    for(int i=0; i<na; i++){  printf("[%-3d] type %-3d %8s  Q=%7.3f ep?=%i  apos(%8.3f,%8.3f,%8.3f)\n",  i, atypes[i], params->atypes[atypes[i]].name, aq[i], aisep[i], apos[i].x, apos[i].y, apos[i].z );}
    // Print electron pair bonds if present
    if(bEpairs && nep>0 && bs){
        printf("Electron pairs\n");
        //printf("idx  host_atom  e-pair\n");
        //printf("----------------------\n");
        for(int i=0; i<nep; i++){ printf("[%-3d]  ihost,iep(%-3d,%-3d) dir(%6.3f,%6.3f,%6.3f)\n", i, bs[i].x, bs[i].y, ad->dirs[i].x,ad->dirs[i].y,ad->dirs[i].z ); }
    }
    //printf("\n=== End Debug Print ===\n\n");
}

__attribute__((hot)) 
void acumDerivs_omp_atomic( int n, int* types, double dEw, Quat4d* fs ){
    for(int i=0; i<n; i++){
        int t            = types[i];     // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];  // get index of degrees of freedom for atom type t
        const Quat4d  f  = fs[i];
        if(tt.x>=0){ 
            #pragma omp atomic
            fDOFs[tt.x]+=f.x*dEw; }
        if(tt.y>=0){ 
            #pragma omp atomic
            fDOFs[tt.y]+=f.y*dEw; }
        if(tt.z>=0){ 
            #pragma omp atomic
             fDOFs[tt.z]+=f.z*dEw; }
        if(tt.w>=0){ 
            #pragma omp atomic
            fDOFs[tt.w]+=f.w*dEw; }
        //if((tt.y>=0)&&(i==0))printf( "acumDerivs i= %i f= %g f*dEw= %g fDOFs= %g\n", i, f.y, f.y*dEw, fDOFs[tt.y] );
        //if((tt.x>=0)||(tt.y>=0))printf( "acumDerivs i= %i dE= %g f.x= %g fDOFs= %g f.y= %g fDOFs= %g\n", i, dE, f.x, fDOFs[tt.x], f.y, fDOFs[tt.y] );
    }
    //exit(0);
}

__attribute__((hot)) 
void acumHostDerivs_omp_atomic( int nepair, Vec2i* epairAndHostIndex, int* types, double dEw, Quat4d* fs  ){
    for(int i=0; i<nepair; i++){
        Vec2i ab         = epairAndHostIndex[i];
        int t            = types[ab.i];      // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];      // get index of degrees of freedom for atom type t
        const Quat4d&  f = fs[ab.j];         // get variation from the host atom
        if(tt.z>=0){ 
            #pragma omp atomic
            fDOFs[tt.z]  -= f.z*dEw; } // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
    }
}

__attribute__((hot)) 
void acumDerivs( int n, int* types, double dEw, Quat4d* fs, double* fDOFs ){
    for(int i=0; i<n; i++){
        int t            = types[i];     // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];  // get index of degrees of freedom for atom type t
        const Quat4d  f  = fs[i];
        if(tt.x>=0){ fDOFs[tt.x]+=f.x*dEw; }
        if(tt.y>=0){ fDOFs[tt.y]+=f.y*dEw; }
        if(tt.z>=0){ fDOFs[tt.z]+=f.z*dEw; }
        if(tt.w>=0){ fDOFs[tt.w]+=f.w*dEw; }
        //if((tt.y>=0)&&(i==0))printf( "acumDerivs i= %i f= %g f*dEw= %g fDOFs= %g\n", i, f.y, f.y*dEw, fDOFs[tt.y] );
        //if((tt.x>=0)||(tt.y>=0))printf( "acumDerivs i= %i dE= %g f.x= %g fDOFs= %g f.y= %g fDOFs= %g\n", i, dE, f.x, fDOFs[tt.x], f.y, fDOFs[tt.y] );
    }
    //exit(0);
}

__attribute__((hot)) 
void acumHostDerivs( int nepair, Vec2i* epairAndHostIndex, int* types, double dEw, Quat4d* fs, double* fDOFs  ){
    for(int i=0; i<nepair; i++){
        Vec2i ab         = epairAndHostIndex[i];
        int t            = types[ab.i];      // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];      // get index of degrees of freedom for atom type t
        const Quat4d&  f = fs[ab.j];         // get variation from the host atom
        if(tt.z>=0){ fDOFs[tt.z]  -= f.z*dEw; } // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
    }
}


__attribute__((hot)) 
void acumDerivs_corr( int n, int* types, int* host, double dEw, Quat4d* fs, double* fDOFs ){
    for(int i=0; i<n; i++){
        int t            = types[i];     // map atom index i to atom type t
        int ih           = host[i];
        const Quat4i& tt = typToREQ[t];  // get index of degrees of freedom for atom type t
        const Quat4d  f  = fs[i];
        if(tt.x>=0){ fDOFs[tt.x]+=f.x*dEw; }
        if(tt.y>=0){ fDOFs[tt.y]+=f.y*dEw; }
        if(tt.z>=0){ fDOFs[tt.z]+=f.z*dEw; }
        if(tt.w>=0){ fDOFs[tt.w]+=f.w*dEw; }
        if(ih>=0){  // electron pair
            const Quat4i& tth = typToREQ[ih];       // get index of degrees of freedom for atom type ih
            if(tth.z>=0){ fDOFs[tth.z]-=f.z*dEw; }  // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
        }
    }
}



/**
 * Calculates the correction to the electrostatic energy and its derivative with respect to the charge.
 */
__attribute__((hot)) 
double corr_elec( double ir, double ir2, double Q, Vec3d d, Vec3d* dirs, int i, int nep, Vec2i* bs, int nj, Vec3d* pos, int j, Vec3d* ps, double &dE_dQ){
    double dE_dr = ir * ir2 * COULOMB_CONST * Q * d.dot(dirs[i]);
    for(int k=0; k<nep; k++){
        if(bs[k].y==nj+i){
            Vec3d dd = pos[j] - ps[bs[k].x-nj];
            dE_dQ -= COULOMB_CONST / dd.norm();
            break;
        }
    }
    return dE_dr;
}


void printAtomsParams( int i0, int n, int* types, Vec3d* ps, Quat4d* typeREQs, double* Qs ){
    printf("printAtomsParams()\n");
    //printf("# i  ti   pi(xyz)   REQi(REQH)   Qi\n");
    for(int ii=0; ii<n; ii++){
        const int i=i0+ii;
        const Quat4d& REQi = typeREQs[types[i]];
        printf("Atom[%3i] %3i=%-8s pos(%7.3f,%7.3f,%7.3f) REQH(%7.3f,%7.3f,%7.3f,%7.3f) Q: %7.3f \n",  i, types[i], params->atypes[types[i]].name, ps[i].x, ps[i].y, ps[i].z,REQi.x, REQi.y, REQi.z, REQi.w, Qs[i] );
    }
}

void printAtomParamDerivs( int na, Quat4d* dEdREQs, int isamp ){
    printf( "printAtomParamDerivs() isamp: %i natoms: %i Forces: \n", isamp, na );
    for(int i=0; i<na; i++){
        printf( "dEdREQs[%3i] %12.3f %12.3f %12.3f %12.3f\n", i, dEdREQs[i].x, dEdREQs[i].y, dEdREQs[i].z, dEdREQs[i].w );
    }
}

void saveDebugXYZ( int i0, int n, int* types, Vec3d* ps, double* Qs, const char* fname="debug.xyz", const char* comment="" ){
    FILE* fout = fopen(fname,"a");
    fprintf(fout, "%i\n", n );
    fprintf(fout, "%s\n", comment );
    for(int ii=0; ii<n; ii++){ // loop over all atoms[i] in system
        const int      i    = i0+ii;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        fprintf(fout, "%c %7.3f %7.3f %7.3f    %7.3f \n", params->atypes[ti].name[0], pi.x, pi.y, pi.z, Qi );
    }
    fclose(fout);
}

double findRmin( Atoms* atoms, Vec2i* inds=0 )const{
    double rmin = 1e+300;
    int ni = atoms->n0;
    int n0 = atoms->n0;
    int nj = atoms->natoms - n0;
    for(int i=0; i<ni; i++){ // loop over all atoms[i] in system
        int itype        = atoms->atypes[i];
        if(params->atypes[itype].name[0]=='E'){ continue; } // ignore electron pairs ( ToDo: should be done in more robuts way )
        const Vec3d&  pi = atoms->apos[i ]; 
        for(int jj=0; jj<nj; jj++){ 
            const int   j        = n0+jj;
            const int   jtype    = atoms->atypes[j];
            if(params->atypes[jtype].name[0]=='E'){ continue; } // ignore electron pairs ( ToDo: should be done in more robuts way )
            const Vec3d      dij = atoms->apos[j] - pi;
            const double r       = dij.norm();
            if(r<rmin){
                rmin=r;
                if(inds){ inds->i=i; inds->j=j; }
            }
        }
    }
    return rmin;
}

bool checkSampleRepulsion( double Eij, int i, int j, int ti, int tj, double r, bool bPrint=true, bool bSave=true )const{
    if(Eij<EijMax) return false;
    char comment[256];
    //char comment[256];
    double rmin = findRmin( samples[isamp_debug] );
    sprintf(comment, "checkSampleRepulsion(isamp=%3i)[%3i,%3i] (%8s,%8s) r: %20.10f rmin: %20.10f ELJ:  %20.10f Eref: %20.10f", isamp_debug, i,j, params->atypes[ti].name, params->atypes[tj].name , r, rmin, Eij, samples[isamp_debug]->Energy  );
    //printf( "evalExampleDerivs_LJQH2(isamp=%3i)[%3i,%3i] (%8s,%8s) r: %20.10f ELJ:  %20.10f \n", isamp_debug, i,j, params->atypes[ti].name, params->atypes[tj].name , r, ELJ  ); 
    if(bPrint)printf( "%s\n", comment );
    if(bSave)samples[isamp_debug]->saveXYZ( "debug.xyz", "a", true, true, Vec3i{1,1,1}, comment );
    //saveDebugXYZ( 0, ni+nj, types, ps, typeREQs, Qs, "debug.xyz", comment );
    //iBadFound++; if(iBadFound>=nBadFoundMax){ printf("ERROR in evalExampleDerivs_LJQH2(): too many bad pairs (%i) \n", iBadFound); exit(0); }
    return true;
}
    
// Optimized versions of evaluation functions
__attribute__((hot)) 
double evalExampleDerivs_LJQH2_corr( int i0, int ni, int j0, int nj, Vec3d* ps, Quat4d* REQs, int* host, Quat4d* dEdREQs ) const {
    double Etot = 0.0;
    for(int ii=0; ii<ni; ii++) {
        const int     i     = i0+ii;
        const Vec3d&  pi    = ps[i];
        const int     ih    = host[i];
        const bool    bEpi  = ih>=0;  // is this atom an electron pair?
        const Quat4d& REQi  = REQs[ii];
        Quat4d        fREQi = Quat4dZero;
        for(int jj=0; jj<nj; jj++) {
            const int j = j0+jj;
            const Quat4d& REQj = REQs[j];
            double H         = REQi.w * REQj.w;   
            if( H>-1e+300 ){ continue; } // only atractive H-bonds
            int jh           = host[j];
            const double R0  = REQi.x + REQj.x;
            const double E0  = REQi.y * REQj.y;
            const double Q   = REQi.z + REQj.z;
            const Vec3d  dij = ps[j] - pi;
            const double r       = dij.norm();
            const double ir      = 1/r;
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;
            Etot += Eel;
            if       (bEpi ){ // i is an electron pair
                double dE_dH;
                Etot     += getEpairAtom( r, H, REQi.x, dE_dH );
                fREQi.z  += dE_dQ * REQj.z;  // dEtot/dH2i
                fREQi.w  += dE_dH * REQj.w;  // dEtot/dH2i
            }else if (jh>=0){ // j is an electron pair
                double dE_dH;
                Etot     += getEpairAtom( r, H, REQj.x, dE_dH ); // Note - radius of electron REQi.x is used for decay constant
                fREQi.w  += dE_dH * REQi.w;  // dEtot/dH2i
            }else{           // both i and j are real atoms
                const double u       = R0/r;
                const double u3      = u*u*u;
                const double u6      = u3*u3;
                const double u6p     = (1.0 + H) * u6;
                const double dE_dE0  = u6 * (u6p - 2.0);
                const double dE_dH   = -E0 * u6 * u6;
                Etot                +=  E0 * dE_dE0;
                fREQi.w             +=  dE_dH * REQj.w;     // dEtot/dH2i
            }
        }
        dEdREQs[i].add(fREQi);
    }
    return Etot;
}
    
__attribute__((hot)) 
//double evalExampleDerivs_LJQH2( int i0, int ni, int j0, int nj, int*  types, Vec3d* ps, Quat4d*  typeREQs, double* Qs, Quat4d* dEdREQs )const{
double evalExampleDerivs_LJQH2( int i0, int ni, int j0, int nj, int* __restrict__ types, Vec3d* __restrict__ ps, Quat4d* __restrict__ typeREQs, double* __restrict__ Qs, Quat4d* __restrict__ dEdREQs )const{
    double Etot = 0.0;
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system
        const int      i    = i0+ii;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;
        for(int jj=0; jj<nj; jj++){ 
            const int   j        = j0+jj;
            const double     Qj  = Qs[j];
            const Vec3d      dij = ps[j] - pi;
            const int        tj  = types[j];
            const Quat4d& REQj   = typeREQs[tj];
            const double R0      = REQi.x + REQj.x;
            const double E0      = REQi.y * REQj.y; 
            const double Q       = Qi     * Qj    ;
            double       H2      = REQi.w * REQj.w;      // expected H2<0
            const double sH      = (H2<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H2 *= sH;
            // --- Electrostatic
            const double r       = dij.norm();
            const double ir      = 1/r;
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;
            const double u       = R0/r;
            const double u3      = u*u*u;
            const double u6      = u3*u3;
            const double u6p     = ( 1.0 + H2 ) * u6;                     
            const double dE_dE0  = u6 * ( u6p - 2.0 );
            const double dE_dR0  = 12.0 * (E0/R0) * u6 * ( u6p - 1.0 );
            const double dE_dH2  = -E0 * u6 * u6;
            const double ELJ     =  E0 * dE_dE0;
            //if(bCheckRepulsion)[[unlikely]]{ checkSampleRepulsion( ELJ, i,j, ti,tj, r, true, true ); }
            // --- Energy and forces
            Etot    +=  ELJ + Eel;
            //printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) ELJ:  %20.10f   Eel: %20.10f \n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel  );
            //{ int itypPrint=4; if( (ti==itypPrint) || (tj==itypPrint) ){ printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) ELJ,Eel: %12.3e,%12.3e Q(%12.3e|%12.3e,%12.3e) dEdREQH(%12.3e,%12.3e,%12.3e,%12.3e)\n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel, Q,Qi,Qj,  dE_dR0, dE_deps, dE_dQ, dE_dH2  ); } }
            fREQi.x +=  dE_dR0;                    // dEtot/dR0_i
            fREQi.y +=  dE_dE0  * 0.5 * REQj.y;    // dEtot/dE0_i
            fREQi.z +=  dE_dQ   * Qj;              // dEtot/dQ_i
            fREQi.w +=  dE_dH2  * REQj.w * sH;     // dEtot/dH2i
            if(bWriteJ){ dEdREQs[j].add( Quat4d{
                        dE_dR0,                    // dEtot/dR0_j
                        dE_dE0  * 0.5 * REQi.y,    // dEtot/dE0_j
                        dE_dQ   * Qi,              // dEtot/dQ_j
                        dE_dH2  * REQi.w * sH,     // dEtot/dH2j
            }); }
        }
        dEdREQs[i].add(fREQi);
    }
    // printAtomParamDerivs( ni+nj, dEdREQs, isamp_debug );
    //printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

__attribute__((hot)) 
double evalExampleDerivs_MorseQH2( int i0, int ni, int j0, int nj, int* __restrict__ types, Vec3d* __restrict__ ps, Quat4d* __restrict__ typeREQs, double* __restrict__ Qs, Quat4d* __restrict__ dEdREQs )const{
    double Etot = 0.0;
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system
        const int      i    = i0+ii;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;
        for(int jj=0; jj<nj; jj++){ 
            const int   j        = j0+jj;
            const double     Qj  = Qs[j];
            const Vec3d      dij = ps[j] - pi;
            const int        tj  = types[j];
            const Quat4d& REQj   = typeREQs[tj];
            const double R0      = REQi.x + REQj.x;
            const double E0      = REQi.y * REQj.y; 
            const double Q       = Qi     * Qj    ;
            double       H2      = REQi.w * REQj.w;      // expected H2<0
            const double sH      = (H2<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H2 *= sH;
            // --- Electrostatic
            double r      = dij.norm();
            double ir     = 1/r;    
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;
            // --- Morse
            //double alpha   = 6.0 / R0;
            double alpha   = 1.6; 
            double e       = exp( -alpha * ( r - R0 ) );
            double e2      = e * e;
            double e2p     = ( 1.0 - H2 ) * e2;
            double dE_dE0  = e2p - 2.0 * e;
            double dE_dR0  = 2.0 * alpha * E0 * ( e2p - e );
            double dE_dH2  = - E0 * e2;
            double ELJ     = E0 * dE_dE0;
            //if(bCheckRepulsion)[[unlikely]]{ checkSampleRepulsion( ELJ, i,j, ti,tj, r, true, true ); }
            // --- Energy and forces
            Etot    +=  ELJ + Eel;
            //printf( "evalExampleDerivs_MorseQH2()[%3i,%3i] (%8s,%8s) ELJ:  %20.10f   Eel: %20.10f \n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel  );
            //{ int itypPrint=4; if( (ti==itypPrint) || (tj==itypPrint) ){ printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) ELJ,Eel: %12.3e,%12.3e Q(%12.3e|%12.3e,%12.3e) dEdREQH(%12.3e,%12.3e,%12.3e,%12.3e)\n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel, Q,Qi,Qj,  dE_dR0, dE_deps, dE_dQ, dE_dH2  ); } }
            fREQi.x +=  dE_dR0;                    // dEtot/dR0_i
            fREQi.y +=  dE_dE0  * 0.5 * REQj.y;    // dEtot/dE0_i
            fREQi.z +=  dE_dQ   * Qj;              // dEtot/dQ_i
            fREQi.w +=  dE_dH2  * REQj.w * sH;     // dEtot/dH2i
            if(bWriteJ){ dEdREQs[j].add( Quat4d{
                        dE_dR0,                    // dEtot/dR0_j
                        dE_dE0  * 0.5 * REQi.y,    // dEtot/dE0_j
                        dE_dQ   * Qi,              // dEtot/dQ_j
                        dE_dH2  * REQi.w * sH,     // dEtot/dH2j
            }); }
            //printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
        dEdREQs[i].add(fREQi);
    }
    //printAtomParamDerivs( ni+nj, dEdREQs, isamp_debug );
    //printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

__attribute__((hot)) 
inline static double funcVar_LJQH2( double r, double R0, double E0, double Q, double H, Quat4d& dEdREQH ){
    //const double ir2     = 1.0 / dij.norm2();
    //const double ir      = sqrt( ir2 );
    const double ir      = 1/r;
    const double ir2     = ir*ir;
    const double u2      = ir2 * ( R0 * R0 );
    const double u4      = u2 * u2;
    const double u6      = u4 * u2;
    const double u6p     = ( 1.0 + H ) * u6;                     
    dEdREQH.x            = 12.0 * (E0/R0) * u6 * ( u6p - 1.0 );
    dEdREQH.y            = u6 * ( u6p - 2.0 );
    dEdREQH.z            = ir * COULOMB_CONST;
    dEdREQH.w            = -E0 * u6 * u6;
    return E0*dEdREQH.y + Q*dEdREQH.z;
}

template <typename Func>
double evalExampleDerivs( Func func, int i0, int ni, int j0, int nj, int* types, Vec3d* ps, Quat4d* typeREQs, double* Qs, Quat4d* dEdREQs )const{
    double Etot = 0.0;
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system
        const int      i    = i0+ii;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;
        for(int jj=0; jj<nj; jj++){
            const int   j        = j0+jj;
            const double     Qj  = Qs[j];
            const Vec3d      dij = ps[j] - pi;
            const int        tj  = types[j];
            const Quat4d& REQj   = typeREQs[tj];
            double        H      = REQi.w * REQj.w;      // expected H2<0
            const double  sH     = (H<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H *= sH;
            // --- Electrostatic
            Quat4d fREQH;
            Etot    += func( dij.norm(), REQi.x+REQj.x, REQi.y*REQj.y, Qi*Qj, H, fREQi );
            // --- Energy and forces
            fREQi.x +=  fREQH.x;                   // dEtot/dR0_i
            fREQi.y +=  fREQH.y * 0.5 * REQj.y;    // dEtot/dE0_i
            fREQi.z +=  fREQH.z   * Qj;            // dEtot/dQ_i
            fREQi.w +=  fREQH.w  * REQj.w * sH;    // dEtot/dH2i
            if(bWriteJ){ dEdREQs[j].add( Quat4d{
                        fREQH.x,                   // dEtot/dR0_j
                        fREQH.y * 0.5 * REQi.y,    // dEtot/dE0_j
                        fREQH.z   * Qi,            // dEtot/dQ_j
                        fREQH.w  * REQi.w * sH,    // dEtot/dH2j
            }); }
        }
        dEdREQs[i].add(fREQi);
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}


// ======================================
// =========  OPTIMIZE  =================
// ======================================

__attribute__((hot)) 
double run( int nstep, double Fmax, double dt, int imodel_, int ialg, bool bOMP, bool bClamp, double max_step ){
    imodel=imodel_;
    double Err=0;
    if( verbosity>1){ printf( "FitREQ::run() nstep %i Fmax %g dt %g isamp %i \n", nstep, Fmax, dt ); }
    if(bOMP){ bBroadcastFDOFs=true; realloc_sample_fdofs();  }
    double F2max=Fmax*Fmax;
    double F2;
    for(int i=0; i<nstep; i++){
        //printf("[%i]  DOFs=", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
        DOFsToTypes(); 
        clean_fDOFs();
        //printf("[%i]  DOFs=", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
        if(bOMP){ 
            Err = evalSamples_omp(); 
            reduce_sample_fdofs();
        }
        else    { Err = evalSamples_noOmp(); }
        if( verbosity>0)printStepDOFinfo( i, Err, "FitREQ::run() BEFORE REGULARIZATION" );
        if(bRegularize){ regularizeDOFs(); }      
        //if( verbosity>0)printStepDOFinfo( i, Err, "FitREQ::run() AFTER  REGULARIZATION" );
        //exit(0);        
        switch(ialg){
            case 0: F2 = move_GD( dt, max_step ); break;
            case 1: F2 = move_MD( dt, max_step ); break;
            case 2: F2 = move_GD_BB_short( i, dt, max_step ); break;
            case 3: F2 = move_GD_BB_long( i, dt, max_step ); break;
            case 4: F2 = move_MD_nodamp( dt, max_step ); break;
        }
        // regularization must be done before evaluation of derivatives
        if(bClamp     ){ limit_params();  }
        //printf("step= %i dt= %g\n", i, dt );
        //printStepDOFinfo( i, Err, "FitREQ::run() AFTER MOVE" );
        //if( F2<0.0   ){ printf("DYNAMICS STOPPED after %i iterations \n", i); printf("VERY FINAL DOFs= ");for(int j=0;j<nDOFs;j++){ printf("%.15g ",DOFs[j]); };printf("\n"); return Err; }
        if( F2<F2max ){ printf("CONVERGED in %i iterations \n", i); }
    }
    printf("VERY FINAL  DOFs= "); for(int j=0;j<nDOFs;j++){ printf("%.15g ", DOFs[j]); };printf("\n");
    printf("VERY FINAL fDOFs= "); for(int j=0;j<nDOFs;j++){ printf("%.15g ",fDOFs[j]); };printf("\n");
    return Err;
}

__attribute__((hot)) 
double run_omp( int nstep, double Fmax, double dt, int imodel_, int ialg, bool bClamp, double max_step ){
    imodel=imodel_;
    double Err=0;
    if( verbosity>1){ printf( "FitREQ::run() nstep %i Fmax %g dt %g isamp %i \n", nstep, Fmax, dt ); }
    double F2max=Fmax*Fmax;
    double F2;
    int nsamp = samples.size();
    double Error=0;
    int itr=0;
    { bBroadcastFDOFs=true; realloc_sample_fdofs();  }
    #pragma omp parallel shared(itr, Error, nsamp, ialg, nstep, dt, max_step, F2, bRegularize, bClamp, verbosity )
    while(itr<nstep){
        #pragma omp single
        {
            DOFsToTypes(); 
            clean_fDOFs();
            Error = 0.0;
        }
        #pragma omp for reduction(+:Error)
        for(int i=0; i<nsamp; i++){
            double E; 
            Error+=evalSampleError( i, E );
        }
        #pragma omp single
        {
            reduce_sample_fdofs();
            if( verbosity>0)printStepDOFinfo( itr, Err, "FitREQ::run() BEFORE REGULARIZATION" );
            if(bRegularize){ regularizeDOFs(); }      
            switch(ialg){
                case 0: F2 = move_GD         (      dt, max_step ); break;
                case 1: F2 = move_MD         (      dt, max_step ); break;
                case 2: F2 = move_GD_BB_short( itr, dt, max_step ); break;
                case 3: F2 = move_GD_BB_long ( itr, dt, max_step ); break;
                case 4: F2 = move_MD_nodamp  (      dt, max_step ); break;
            }
            if(bClamp     ){ limit_params();  }
            if( F2<F2max ){ 
                printf("CONVERGED in %i iterations \n", itr); 
                itr = nstep;  // to terminate the while-loop
            }
            itr++;
        }
    } // while(itr<nstep){
    printf("VERY FINAL  DOFs= "); for(int j=0;j<nDOFs;j++){ printf("%.15g ", DOFs[j]); };printf("\n");
    printf("VERY FINAL fDOFs= "); for(int j=0;j<nDOFs;j++){ printf("%.15g ",fDOFs[j]); };printf("\n");
    return Error;
}

/**
 * Calculates the regularization force for each degree of freedom (DOF) based on the difference between the current and target values of the fitted non-colvalent interaction parameters REQH(Rvdw,Evdw,Q,Hb).
 * The regularization force tries to minimize the difference between the current (REQ) and the default value (REQ0) of the parameters. Its strength is controlled by the regularization stiffness (Kreg).
 * The resulting forces are stored in the fDOFs array.
 */
__attribute__((hot)) 
void regularization_force(){
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt       = typToREQ [i];
        const Quat4d& REQ      = typeREQs [i];
        const Quat4d& REQ0     = typeREQs0[i];
        const Quat4d& K        = typeKreg [i];
        if(tt.x>=0)fDOFs[tt.x] = (REQ0.x-REQ.x)*K.x;
        if(tt.y>=0)fDOFs[tt.y] = (REQ0.y-REQ.y)*K.y;
        if(tt.z>=0)fDOFs[tt.z] = (REQ0.z-REQ.z)*K.z;
        if(tt.w>=0)fDOFs[tt.w] = (REQ0.w-REQ.w)*K.w;
    }
}

__attribute__((hot)) 
inline double constrain( double x, double xmin, double xmax, double Kmin, double Kmax, double& f ){
    double E=0;
    if(x<xmin){
        double d = x-xmin; 
        f =      -d*Kmin;
        E = 0.5*d*d*Kmin;
    }else if(x>xmax){
        double d = x-xmax; 
        f =      -d*Kmax;
        E = 0.5*d*d*Kmax;
    }
    return E;
}

__attribute__((hot)) 
void regularization_force_walls(){
    double E = 0;
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt    = typToREQ      [i];
        const Quat4d& REQ   = typeREQs      [i];
        const Quat4d& REQ0l = typeREQs0_low [i];
        const Quat4d& Kl    = typeKreg_low  [i];
        const Quat4d& REQ0h = typeREQs0_high[i];
        const Quat4d& Kh    = typeKreg_high [i];
        if(tt.x>=0){ E+= constrain( REQ.x, REQ0l.x, REQ0h.x, Kl.x,Kh.x, fDOFs[tt.x] ); }
        if(tt.x>=0){ E+= constrain( REQ.y, REQ0l.y, REQ0h.y, Kl.y,Kh.y, fDOFs[tt.y] ); }
        if(tt.x>=0){ E+= constrain( REQ.z, REQ0l.z, REQ0h.z, Kl.z,Kh.z, fDOFs[tt.z] ); }
        if(tt.x>=0){ E+= constrain( REQ.w, REQ0l.w, REQ0h.w, Kl.w,Kh.w, fDOFs[tt.w] ); }
    }
}

// New regularization function operating per-DOF
__attribute__((hot)) 
double regularizeDOFs(){
    double Etot = 0;
    for(int i=0; i<nDOFs; i++){
        const Vec3d& regX = DOFregX[i];
        const Vec3d& regK = DOFregK[i];
        // Hard walls outside [xmin,xmax]
        Etot += constrain(DOFs[i], regX.x, regX.z, regK.x, regK.z, fDOFs[i]);
        // Soft harmonic around x0
        if(regK.y > 0){
            double d  = DOFs[i] - regX.y;
            fDOFs[i] -= d * regK.y;
            Etot     += 0.5 * d * d * regK.y;
        }
    }
    return Etot;
}

/**
 * Limits the fitted non-colvalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) to be btween minimum and maximum (REQmin and REQmax).
 * 
 * @param none
 * @return void
 */
__attribute__((hot)) 
void limit_params(){
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt     = typToREQ[i];
        const Quat4d& REQ    = typeREQs [i];
        const Quat4d& REQmin = typeREQsMin[i];
        const Quat4d& REQmax = typeREQsMax[i];
        // Note: should we limit derivatives or the value itself?
        if(tt.x>=0){ DOFs[tt.x]=_clamp(DOFs[tt.x],REQmin.x,REQmax.x ); }
        if(tt.y>=0){ DOFs[tt.y]=_clamp(DOFs[tt.y],REQmin.y,REQmax.y ); }
        if(tt.z>=0){ DOFs[tt.z]=_clamp(DOFs[tt.z],REQmin.z,REQmax.z ); }
        if(tt.w>=0){ DOFs[tt.w]=_clamp(DOFs[tt.w],REQmin.w,REQmax.w ); }
        //if(tt.x>=0){ fDOFs[tt.x]=_clamp(fDOFs[tt.x],REQmin.x,REQmax.x ); }
        //if(tt.y>=0){ fDOFs[tt.y]=_clamp(fDOFs[tt.y],REQmin.y,REQmax.y ); }
        //if(tt.z>=0){ fDOFs[tt.z]=_clamp(fDOFs[tt.z],REQmin.z,REQmax.z ); }
        //if(tt.w>=0){ fDOFs[tt.w]=_clamp(fDOFs[tt.w],REQmin.w,REQmax.w ); }
    }
}

/**
 * This function limits the time step based on the maximum step size and the actual maximum force magnitude in the fDOFs array. It is used in the gradient descent and molecular dynamics optimization algorithms to make sure it is stable.
 * @param dt The time step to be limited.
 * @return The limited time step.
 */
__attribute__((hot)) 
double limit_dt(double dt, double max_step){
    double fm=0;
    int ifm;
    //for(int i=0; i<nDOFs; i++){fm=_max(fm,fabs(fDOFs[i]));}
    for(int i=0; i<nDOFs; i++){double f=fabs(fDOFs[i]/DOFs[i]); if(f>fm){fm=f; ifm=i;}}
    if(dt*fm>max_step){
        dt = max_step/fm;
        //printf( "limit_dt %g max variation[%i] %g old %g -> new %g \n", dt, ifm+1, -fDOFs[ifm]/DOFs[ifm]*dt, DOFs[ifm], DOFs[ifm]-fDOFs[ifm]*dt );
        printf( "limit_dt %g\n", dt );
    }
    return dt;
}

__attribute__((hot)) 
double limit_dt_MD(double dt, double max_step, double cdamp){
    double sm = 0.0;
    int ism;
    for(int i=0; i<nDOFs; i++){double s=fabs((vDOFs[i]*cdamp-fDOFs[i]*dt)*dt); if(s>sm){sm=s;ism=i;}}
    if(sm>max_step){
        double v = vDOFs[ism] * cdamp;
        double finv = 1.0/fDOFs[ism];
        dt = sqrt( 0.25*v*v*finv*finv + max_step*fabs(finv) ) - 0.5*fabs(v*finv);
//sm=0.0;for(int i=0; i<nDOFs; i++){double s=fabs((vDOFs[i]*cdamp-fDOFs[i]*dt)*dt); if(s>sm){sm=s;}}
        printf( "limit_dt %g\n", dt );
    }
//printf( "check_dt dt=%g max_step=%g\n", dt, sm );
    return dt;
}

__attribute__((hot)) 
double limit_dt_MD_nodamp(double dt, double max_step){
    double sm = 0.0;
    int ism;
    for(int i=0; i<nDOFs; i++){double s=fabs((vDOFs[i]-fDOFs[i]*dt)*dt); if(s>sm){sm=s;ism=i;}}
    if(sm>max_step){
        double v = vDOFs[ism];
        double finv = 1.0/fDOFs[ism];
        dt = sqrt( 0.25*v*v*finv*finv + max_step*fabs(finv) ) - 0.5*fabs(v*finv);
//sm=0.0;for(int i=0; i<nDOFs; i++){double s=fabs((vDOFs[i]-fDOFs[i]*dt)*dt); if(s>sm){sm=s;}}
        printf( "limit_dt %g\n", dt );
    }
//printf( "check_dt dt=%g max_step=%g\n", dt, sm );
    return dt;
}

/**
 * Moves fitting parameters using the gradient descent (GD) method in order to minimize the fitting error.
 * @param dt The time step ( the higher the value, the faster the convergence, but the less stable the algorithm is).
 * @return Sum of squares of the variatinal derivatives of the fitting error with respect to all fitting parameters.
 */
__attribute__((hot)) 
double move_GD( double dt, double max_step ){
    //printf("now in move_GD\n");
    double F2 = 0;
    dt=limit_dt(dt,max_step);
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        DOFs[i] -= f*dt;
        F2 += f*f;
    }
    return F2;
}

// compute optimal dt according to the Barzilai-Borwein method
__attribute__((hot)) 
double move_GD_BB_short( int step, double dt, double max_step ){
    double F2 = 0;
    if(step>0){
        double dxdg = 0.0;
        double dgdg = 0.0;
        for(int i=0; i<nDOFs; i++){
            double dx = DOFs[i] - DOFs_old[i];
            double dg = fDOFs[i] - fDOFs_old[i];
            dxdg += dx*dg;
            dgdg += dg*dg;
        }
        dt = fabs(dxdg)/dgdg;
        //dt = dxdg/dgdg;
        if(fabs(dxdg)<=1e-10*dgdg){dt=1e-11;}
    }
    dt=limit_dt(dt,max_step);
    printf("step= %i dt= %g\n", step, dt );
    for(int i=0; i<nDOFs; i++){
        DOFs_old[i] = DOFs[i];
        fDOFs_old[i] = fDOFs[i];
        double f = fDOFs[i];
        DOFs[i] -= f*dt;
        F2 += f*f;
    }
    // stop the algorithm if the step is too small
    //if(dt<1e-10){ return -F2; }
    return F2;
}

// compute optimal dt according to the Barzilai-Borwein method
__attribute__((hot)) 
double move_GD_BB_long( int step, double dt, double max_step ){
    double F2 = 0;
    if(step>0){
        double dxdx = 0.0;
        double dxdg = 0.0;
        for(int i=0; i<nDOFs; i++){
            double dx = DOFs[i] - DOFs_old[i];
            double dg = fDOFs[i] - fDOFs_old[i];
            dxdx += dx*dx;
            dxdg += dx*dg;
        }
        dt = dxdx/fabs(dxdg); // long BB step
        //dt = dxdx/dxdg; // long BB step
        if(dxdx<=1e-10*fabs(dxdg)){dt=1e-11;}
    }
    dt=limit_dt(dt,max_step);
    printf("step= %i dt= %g\n", step, dt );
    for(int i=0; i<nDOFs; i++){
        DOFs_old[i] = DOFs[i];
        fDOFs_old[i] = fDOFs[i];
        double f = fDOFs[i];
        DOFs[i] -= f*dt;
        F2 += f*f;
    }
    // stop the algorithm if the step is too small
    //if(dt<1e-10){ return -F2; }
    return F2;
}

/**
 * Moves fitting parameters using the damped molecular dynamics (leap-frog algorithm) in order to minimize the fitting error.
 * @param dt The time step ( the higher the value, the faster the convergence, but the less stable the algorithm is).
 * @param damp The damping factor to apply to the velocity of each degree of freedom. Default value is 0.1.
 * @return Sum of squares of the variatinal derivatives of the fitting error with respect to all fitting parameters.
 */
__attribute__((hot)) 
double move_MD( double dt, double max_step=-0.1, double damp=0.1 ){
    double cdamp = 1.0-damp;
    double F2 = 0.0;
    if(max_step>0.0){ dt=limit_dt_MD(dt,max_step,cdamp); };
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        double v = vDOFs[i];
        v       *= cdamp;
        v       -= f*dt;
        vDOFs[i] = v;
        DOFs[i] += v*dt;
        F2      += f*f;
    }
    return F2;
}

__attribute__((hot)) 
double move_MD_nodamp( double dt, double max_step=-0.1 ){
    double fv = 0.0; for(int i=0; i<nDOFs; i++){ fv += fDOFs[i]*vDOFs[i]; }
    double F2 = 0.0;
    if(fv<0.0){ 
        for(int i=0; i<nDOFs; i++){ 
            vDOFs[i] = 0.0; 
        } 
    }
    if(max_step>0.0){ dt=limit_dt_MD_nodamp(dt,max_step); };
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        double v = vDOFs[i];
        //DOFs[i] += v*dt - 0.5*f*dt*dt;
        v       -= f*dt;
        vDOFs[i] = v;
        DOFs[i] += v*dt;
        F2      += f*f;
    }
    return F2;
}

/// @brief Cleans the derivative array by setting all its elements to zero.
void clean_fDOFs(){ for(int i=0; i<nDOFs; i++){fDOFs[i]=0;} }

Vec2d getMinMax( int n, double* vec ){
    Vec2d bd = Vec2d{+1e+300,-1e+300};
    for(int i=0; i<n; i++){
        double v = vec[i];
        if(v<bd.x){ bd.x=v; }
        if(v>bd.y){ bd.y=v; }
    }
    return bd;
} 

void printStepDOFinfo( int istep, double Err, const char* label="" ){
    Vec2d bd = getMinMax( nDOFs, fDOFs );
    if( verbosity>2){ 
        printf( "%s step: %i Err= %g  min,max= %g %g \n", label, istep, Err, bd.x, bd.y );
        printDOFs(); 
    }else if( verbosity>1){
        printf("step= %i  DOFs= ", istep);for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
        printf("step= %i fDOFs= ", istep);for(int j=0;j<nDOFs;j++){ printf("%g ",fDOFs[j]); };printf("\n");
    }
    if( isnan(Err)                   ){ printf( "ERROR in %s step: %i Err= %g \n"        , label, istep, Err ); exit(0); }
    if ( bd.x < -1e+8 || bd.y > 1e+8 ){ printf( "ERROR in %s step: %i Fmin,max= %g %g \n", label, istep, bd.x, bd.y ); exit(0); }
}

const char* REQcomponentToStr(int i){
    switch(i){
        case 0: return "R";
        case 1: return "E";
        case 2: return "Q";
        case 3: return "H";
        default: return "?";
    }
}

void printDOFmapping() {
    printf("=== DOF Mapping ===\n");
    for(int i=0; i<ntype; i++){
        const Quat4i& tt = typToREQ[i];
        if(tt.x>=0 || tt.y>=0 || tt.z>=0 || tt.w>=0){
            printf("Type %3d (%8s):\n", i, params->atypes[i].name);
            if(tt.x>=0) printf("  DOF[%3d] -> %s (vdW radius)\n",     tt.x, REQcomponentToStr(0));
            if(tt.y>=0) printf("  DOF[%3d] -> %s (vdW energy)\n",     tt.y, REQcomponentToStr(1));
            if(tt.z>=0) printf("  DOF[%3d] -> %s (charge)\n",         tt.z, REQcomponentToStr(2));
            if(tt.w>=0) printf("  DOF[%3d] -> %s (H-bond param)\n",   tt.w, REQcomponentToStr(3));
        }
    }
    printf("=================\n");
}

void printDOFsToTypes() const {
    printf("printDOFsToTypes() nDOFs=%i\n", nDOFs);
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = REQtoTyp[i];
        const char comp = "REQH"[rt.y];
        printf("%i -> %i|%i  = %s.%c \n", i, rt.x, rt.y,  params->atypes[rt.x].name, comp);
    }
    //printf("\n");
}

void printTypesToDOFs() const {
    printf("printTypesToDOFs() ntype=%i\n", ntype);
    for(int i=0; i<ntype; i++){
        const Quat4i& tt = typToREQ[i];
        if(tt.x>=0 || tt.y>=0 || tt.z>=0 || tt.w>=0){
            printf("%d(%s): ", i, params->atypes[i].name);
            if(tt.x>=0) printf("R=%d,", tt.x);
            if(tt.y>=0) printf("E=%d,", tt.y);
            if(tt.z>=0) printf("Q=%d,", tt.z);
            if(tt.w>=0) printf("H=%d,", tt.w);
            printf("\n");
        }
    }
    //printf("\n");
}

void printDOFs() const {
    printf("printDOFvalues() nDOFs=%i\n", nDOFs);
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = REQtoTyp[i];
        char* tname = params->atypes[rt.x].name;
        const char comp = "REQH"[rt.y];
        printf("%3i ->(%3i|%i) %8s.%c : %20.10f dE/d%c: %20.10f \n", i, rt.x,rt.y, tname,  comp, DOFs[i], comp, fDOFs[i]);
    }
    //printf("\n");
}

void printDOFregularization(){
    printf( "\nprintDOFregularization()\n" );
    //printf("DOF  Type      Comp  Current       Position(min,x0,max)           Stiffness(min,k0,max)\n");
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt   = REQtoTyp[i];
        const Vec3d& regX = DOFregX[i];
        const Vec3d& regK = DOFregK[i];
        printf("[%2i] %8s.%c(%8.3f) reg: x(%8.3f,%8.3f,%8.3f) K(%8.3f,%8.3f,%8.3f)\n",   i, params->atypes[rt.x].name, "REQH"[rt.y],  DOFs[i],    regX.x, regX.y, regX.z,      regK.x, regK.y, regK.z        );
    }
}

}; // class FitREQ



#endif
