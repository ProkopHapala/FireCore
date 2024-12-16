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

#include "IO_utils.h"

//#include "OptRandomWalk.h"

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

inline double getEpairAtom( double r, double Hij, double w, double& dEdH, double& dEdw ){
    
    double iw = 1.0/w;
    double u = r*iw; 
    dEdH = exp( -u );
    dEdw = dEdH*u*iw;
    //dEdH = exp( -u*u );
    return Hij * dEdH;
    

    // if( r<w ){
    //     return Hij;
    // }else{
    //     return 0.0;
    // }
}

struct AddedData{
    int    nep  =0;      // number of electron pairs
    Vec2i* bs   =0;      // bonds between electron pairs and host atoms  (host_atom_index, epair_index)
    Vec3d* dirs =0;      // directions of electron pairs
    int*   host = 0;       // [natoms]  Index of host atom for electron pairs (-1 for non-electron-pair atoms)

    // New members for H-bond optimization
    double  Emodel0     = NAN;     // Static model energy without H-bond corrections
    int     HBna        = 0;       // Number of atoms with significant H-bond corrections
    int     HBn0        = -1;      // Split point between fragment 1 and 2 in HB-filtered arrays
    int*    HBtoFit     = 0;       // [natom] for each old atom, the index of atom in the fitted array [HBna], -1 if not fitted
    int*    HBatomsInd  = 0;       // [HBna]  original index of this fitted atom 
    //int*    HBatomsHost = 0;       // [HBna]  Index of host atom for electron pairs (-1 for non-electron-pair atoms)
    //int*    HBatomsType = 0;       // [HBna]  Positions of H-bond corrected atoms
    //Vec3d*  HBatomsPos  = 0;       // [HBna]  Positions of H-bond corrected atoms
    //Quat4d* HBatomsREQs = 0;       // [HBna]  REQ parameters of H-bond corrected atoms
    //int*    host = 0;       // [natoms]  Index of host atom for electron pairs (-1 for non-electron-pair atoms)
    //int*    type = 0;       // [natoms]  Positions of H-bond corrected atoms
    Quat4d* REQs = 0;       // [natoms]  REQ parameters of H-bond corrected atoms

    //AddedData() = default;
    //AddedData(int natoms, int nhb) {    realloc(natoms, nhb);}

    void reallocHB(int natoms, int nhb) {
        //printf( "reallocHB(%i, %i) \n", natoms, nhb );
        HBna = nhb;
        // _realloc0( bs,   natoms, Vec2i{-1,-1} );
        // _realloc0( dirs, natoms, Vec3dNAN );
        if(nhb==0) return;
        _realloc0( HBtoFit, natoms, -1 );
        //_realloc0( host,    natoms, -1 );
        _realloc0( REQs,    natoms, Quat4dNAN );
        _realloc0( HBatomsInd, nhb, -1 );
        // _realloc0( HBatomsHost, nhb, -1 );
        // _realloc0( HBatomsType, nhb, -1 );
        // _realloc0( HBatomsREQs, nhb, Quat4dNAN );
        //_realloc0( HBatomsPos,  nhb, Vec3dNAN );
    }

    void deallocHB() {
        // _dealloc( bs   );
        // _dealloc( dirs );
        if(HBna==0) return;
        _dealloc( HBtoFit     );
        //_dealloc( host );
        _dealloc( REQs );
        _dealloc( HBatomsInd  );
        //_dealloc( HBatomsHost );
        //_dealloc( HBatomsType );
        //_dealloc( HBatomsPos  );
        //_dealloc( HBatomsREQs );

    }

    //~AddedData() {dealloc();}

};

// ================================
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
    int nDOFs=0,ntype=0,nbatch=0;
    int imodel   =1;
    int iparallel=0;
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
    // alignas(32)  Vec2i* DOFtoTyp; // Maps DOF index to (type_index, component)
    // alignas(32)  Vec3d* DOFregX;  // regularization positions (xmin,x0,xmax) for each DOF
    // alignas(32)  Vec3d* DOFregK;  // regularization stiffness (Kmin,K0,Kmax) for each DOF

    std::vector<Vec2i> DOFtoTyp;  // Maps DOF index to (type_index, component)
    std::vector<Vec3d> DOFregX;   // regularization positions (xmin,x0,xmax) for each DOF
    std::vector<Vec3d> DOFregK;   // regularization stiffness (Kmin,K0,Kmax) for each DOF
    std::vector<Vec2d> DOFlimits;   // limits (xmin,xmax) for each DOF

    alignas(32) double*   DOFs =0;       // [nDOFs]
    alignas(32) double*   fDOFs=0;       // [nDOFs]
    alignas(32) double*   vDOFs=0;       // [nDOFs]
    alignas(32) double*   DOFs_old =0;   // [nDOFs]
    alignas(32) double*   fDOFs_old=0;   // [nDOFs]
    alignas(32) Vec2d*    fDOFbounds=0; // [nDOFs] minimum and maximum value of fDOFs over all samples

    alignas(32) double*   sample_fdofs = 0; // [nDOFs*nsamples] - arrays for parallelization to avoid atomic-write conflicts

    bool  bEvalJ          = false;    // Should we evaluate variational derivatives on Fregment J 
    bool  bWriteJ         = false;    // Should we write variational derivatives to Fregment J ( inner loop over j )
    bool  bCheckRepulsion = false;    // Should we check maximum repulsion (EijMax) inside inner loop over j for each sample atoms ?
    bool  bRegularize     = true;     // Should we apply additional regularization forces to otimizer ( beside the true variational forces from inter-atomic forcefield ? )
    bool  bAddRegError    = true;     // Should we add regularization error to total error ?
    bool  bEpairs         = true;     // Should we add electron pairs to the molecule ?
    bool  bEpairDistByType = false;   // Should we use different electron pair distances for different atom types ?
    //bool  bOptEpR = false;          // Should we optimize electron pair distance (from host atom) ?
    bool  bBroadcastFDOFs = false;    // Should we broadcast fDOFs (each sample to its own chunk of memory) to prevent atomic-write conflicts ?
    bool  bUdateDOFbounds = true;     // Should we update fDOFbounds after each sample ?
    bool  bClearDOFboundsEpoch = false; // Should we clear fDOFbounds after each epoch ?
    bool  bEvalOnlyCorrections = false;  // Split evaluation and optimization to Emodel0 and Ecorrection (where only Ecorrection is updated every iteration)

    // what to do with samples with E>EmodelCut ?
    bool bListOverRepulsive    = true;   // Should we list overrepulsive samples? 
    bool bSaveOverRepulsive    = false;  // Should we save overrepulsive samples to .xyz file?
    bool bPrintOverRepulsive   = true;   // Should we print overrepulsive samples? 
    bool bDiscardOverRepulsive = true;   // Should we discard overrepulsive samples? ( i.e. ignore them as training examples )
    bool bWeightByEmodel       = true;   // Should we weight samples by their model energy ?
    std::vector<int> overRepulsiveList;  // List of overrepulsive samples indexes (for recent epoch)
    char* fname_overRepulsive  = "overRepulsive.xyz";

    bool bPrintDOFs     = false;
    bool bPrintfDOFs    = true;
    bool bPrintBeforReg = true;
    bool bPrintAfterReg = false;

    bool bSaveSampleToXYZ = false;
    char* xyz_out         = "out.xyz";

    std::vector<int> typesPresent;
    std::vector<int> fittedTypes;

    // parameters
    int    iWeightModel    = 1;    // weight of model energy 1=linear, 2=cubic_smooth_step  
    double EmodelCut       = 10.0; // sample model energy when we consider it too repulsive and ignore it during fitting
    double EmodelCutStart  = 5.0;  // sample model energy when we start to decrease weight in the fitting  
    double invWsum         = 1.0;

    //double kMorse          = 1.4;
    //double kMorse          = 1.5;
    //double kMorse          = 1.6;
    //double kMorse          = 1.7;
    double kMorse            = 1.8;
    double Lepairs           = 1.0;
    

    // check pairwise repulsion betwee atoms within one sample
    double EijMax    = 5.0;
    int isamp_debug  = 0;
    int iBadFound    = 0; 
    int nBadFoundMax = 10;

    double*  weights = 0; // [nbatch] scaling importaince of parameters

    std::vector<Atoms*> samples;    // ToDo: would be more convenient to store Atoms* rather than Atoms

    MM::Builder builder;
    MMFFparams* params=0; 

    // New members for H-bond optimization
    

// =================================
// =========== Functions ===========
// =================================

/**
 * @brief Reallocates the DOFs (degrees of freedom) array to the given size. affects the size of the fDOFs and vDOFs arrays as well.
 * 
 * @param nDOFs_ The new size of the DOFs array.
 */
void realloc_DOFs( int nDOFs_ ){
    //nDOFs=nR+nE+nQ;
    nDOFs=nDOFs_;
    _realloc0(  DOFs, nDOFs, 0.0 );
    _realloc0( fDOFs, nDOFs, 0.0 );
    _realloc0( vDOFs, nDOFs, 0.0 ); 
    _realloc0(  DOFs_old,  nDOFs, 0.0 );
    _realloc0( fDOFs_old,  nDOFs, 0.0 );
    _realloc0( fDOFbounds, nDOFs, Vec2dZero );
    //Rs=DOFs;Es=DOFs+nR;Qs=DOFs+nR+nE;
    //fRs=fDOFs;fEs=fDOFs+nR;fQs=fDOFs+nR+nE;
}

void initDOFs(){
    for (int iDOF=0; iDOF<nDOFs; iDOF++) {
        double xstart = DOFregX[iDOF].y;
        DOFs[iDOF] = xstart;
        fDOFs[iDOF] = 0.0;
        vDOFs[iDOF] = 0.0;
    }
}

void realloc_sample_fdofs(){
    _realloc0( sample_fdofs, nDOFs*samples.size(), 0.0 );
}


void initTypeParamsFromDOFs() {
    for (int iDOF=0; iDOF<nDOFs; iDOF++) {
        Vec2i rt = DOFtoTyp[iDOF];
        int ityp = rt.x;
        int comp = rt.y;
        typToREQ[ityp].array[comp] = iDOF;
        double xstart = DOFregX[iDOF].y;
        printf("ityp=%i comp=%i iDOF=%i xstart=%g\n", ityp, comp, iDOF, xstart);
        typeREQs      [ityp].array[comp] = xstart;  // Initialize with xstart
        typeREQs0     [ityp].array[comp] = xstart;
        typeREQsMin   [ityp].array[comp] = DOFlimits[iDOF].x;
        typeREQsMax   [ityp].array[comp] = DOFlimits[iDOF].y;
        typeREQs0_low [ityp].array[comp] = DOFregX  [iDOF].x;
        typeREQs0_high[ityp].array[comp] = DOFregX  [iDOF].z;
        typeKreg_low  [ityp].array[comp] = DOFregK  [iDOF].x;
        typeKreg_high [ityp].array[comp] = DOFregK  [iDOF].z;
    }
}

int initAllTypes(){
    ntype = params->atypes.size();
    printf( "FitREQ::initAllTypes() ntype=%i \n", ntype );
    reallocTypeParams(ntype);
    for(int i=0; i<ntype; i++){
        Quat4d tREQH = Quat4d{ params->atypes[i].RvdW, sqrt(params->atypes[i].EvdW), params->atypes[i].Qbase, params->atypes[i].Hb };
        typeREQs [i] = tREQH;
        typeREQs0[i] = tREQH;
    }
    fittedTypes.resize(ntype);  for(int i=0; i<ntype; i++){ fittedTypes[i]=0; }  // initialize all types as not fitted
    return ntype;
}

void reallocTypeParams(int ntype_) {
    ntype = ntype_;
    _realloc0( typToREQ,       ntype_, Quat4i{-1,-1,-1,-1} );
    _realloc0( typeREQs,       ntype_, Quat4dZero );
    _realloc0( typeREQs0,      ntype_, Quat4dZero );
    _realloc0( typeREQsMin,    ntype_, Quat4d{ -1e+300, -1e+300, -1e+300, -1e+300 } );
    _realloc0( typeREQsMax,    ntype_, Quat4d{  1e+300,  1e+300,  1e+300,  1e+300 } );
    _realloc0( typeREQs0_low,  ntype_, Quat4dZero );
    _realloc0( typeREQs0_high, ntype_, Quat4dZero );
    _realloc0( typeKreg_low,   ntype_, Quat4dZero );
    _realloc0( typeKreg,       ntype_, Quat4dZero );
    _realloc0( typeKreg_high,  ntype_, Quat4dZero );
}

void reduce_sample_fdofs(){
    int nsamples = samples.size();
    for(int j=0; j<nDOFs; j++){ fDOFs[j] = 0.0; }
    for(int i=0; i<nsamples; i++){
        double* fs = sample_fdofs + i*nDOFs;
        for(int j=0; j<nDOFs; j++){ fDOFs[j] += fs[j]; }
    }
}

void printTypeParams( bool bOnlyPresent=true ){
    printf("printTypeParams():\n");
    int ntype = params->atypes.size();
    bool bCounted = typesPresent.size() > 0;
    for(int i=0; i<ntype; i++){
        int ncount = -1;
        if(bCounted){
            ncount = typesPresent[i];
            if( bOnlyPresent && (ncount==0) ){ continue; }
        }
        Quat4d tREQH = typeREQs[i]; 
        printf("type %3i %-8s count: %6i REQH: %10.3f %10.3f %10.3f %10.3f \n", i, params->atypes[i].name, ncount, tREQH.x, tREQH.y, tREQH.z, tREQH.w );
    }
}

void countTypesPresent( ){
    typesPresent.resize( params->atypes.size() );
    for( int it=0; it<typesPresent.size(); it++ ){ typesPresent[it]=0; }
    int nsamples = samples.size();
    for( int isamp=0; isamp<nsamples; isamp++ ){
        const Atoms* atoms = samples[isamp];
        for(int i=0; i<atoms->natoms; i++){
            int ityp = atoms->atypes[i];
            typesPresent[ityp]+=1;
        }
    }
}

double updateWeightsSum(){
    double wsum = 0.0;
    int nsamples = samples.size();
    if(weights==0){
        wsum=nsamples;
    }else{
        for(int i=0; i<samples.size(); i++){ wsum += weights[i]; }
    }
    invWsum = 1.0/wsum;
    printf("updateWeightsSum() nsamples=%i sum=%g invWsum=%g \n", samples.size(), wsum, invWsum );
    return wsum;
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

/**
 * @brief Load DOF selection from new format file that specifies individual degrees of freedom per type
 * Format: typename comp Min Max xlo xhi Klo Khi xstart
 * where:
 * - typename: Atom type name (e.g. "N_3", "H_O")
 * - comp: Which component (0=R, 1=E, 2=Q, 3=H) 
 * - Min/Max: Hard limits on parameter value
 * - xlo/xhi: Regularization targets
 * - Klo/Khi: Regularization constants
 * - xstart: Initial value
 * @return Number of DOFs loaded
 */
int loadDOFSelection( const char* fname ){
    initAllTypes();
    printf("loadDOFSelection() fname=%s\n", fname);
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    int nw_min    = 1 + 1 + 2 + 2 +2 + 1;
    int nw_xstart = nw_min + 1;
    DOFtoTyp .clear();
    DOFregX  .clear();
    DOFregK  .clear();
    DOFlimits.clear();
    int iDOF = 0;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        int comp;
        double xmin, xmax, xlo, xhi, klo, khi, K0, xstart;
        int nw = sscanf( line, "%s %i %lf %lf %lf %lf %lf %lf %lf %lf",  at_name, &comp, &xmin, &xmax, &xlo, &xhi, &klo, &khi, &K0, &xstart );
        if(nw < nw_min){ printf("ERROR in loadDOFSelection(): expected min. %i words, got %i in line '%s'\n", nw_min, nw, line); exit(0); }
        int ityp = params->getAtomType(at_name);
        if(ityp < 0){  printf("ERROR in loadDOFSelection(): unknown atom type %s\n", at_name); exit(0); }
        // Add new DOF
        fittedTypes[ityp] = 1;
        if(nw<nw_xstart){ 
            printf("iDOF: %i %8s.%i xstart missing (nw=%i) => xstart=%g \n", iDOF, at_name, comp, nw, typeREQs0[ityp].array[comp]  );
            xstart = typeREQs0[ityp].array[comp]; 
        }
        DOFtoTyp .push_back( {ityp, comp}       );
        DOFregX  .push_back( {xlo, xstart, xhi} );
        DOFregK  .push_back( {klo, K0,     khi} );
        DOFlimits.push_back( {xmin, xmax}       );
        //if(verbosity>0)
            printf( "DOF[%3i] %3i|%i %-8s %c  range(%g,%g) reg(x0=(%g,%g),K=(%g,%g)) xstart=%g\n",  iDOF, ityp,comp,  at_name, "REQH"[comp], xmin, xmax, xlo, xhi, klo, khi, xstart );
        iDOF++;
    }
    fclose(fin);
    realloc_DOFs( DOFtoTyp.size() );
    initDOFs();
    initTypeParamsFromDOFs();
    //DOFsToTypes();
    if(verbosity>0){
        printDOFsToTypes();
        printTypesToDOFs();
        printDOFregularization();
    }
    return nDOFs;
}

Atoms* addEpairs( Atoms* mol ){
    builder.params = params;
    builder.clear();
    if(mol->lvec){ builder.lvec = *(mol->lvec); builder.bPBC=true; }
    builder.insertAtoms(*mol);
    //builder.printCappingTypes();
    //builder.printAtomConfs(false, true );
    builder.tryAddConfsToAtoms( 0, -1 );
    builder.cleanPis();
    if( builder.bPBC ){ 
        builder.autoBondsPBC( -.5,  0      , mol->n0     ); // we should not add bonds between rigid and flexible parts
        builder.autoBondsPBC( -.5,  mol->n0, mol->natoms ); 
    }else{ 
        builder.autoBonds( -.5,  0      , mol->n0     ); // we should not add bonds between rigid and flexible parts
        builder.autoBonds( -.5,  mol->n0, mol->natoms ); 
    }
    //builder.printAtomConfs(false, true );
    builder.checkNumberOfBonds( true, true );
    builder.bDummyEpair = true;
    builder.autoAllConfEPi( ); 
    //builder.printAtomConfs(false, true );
    builder.setPiLoop       ( 0, -1, 10 );
    builder.addAllEpairsByPi( 0, -1 );  
    builder.addSigmaHoles();     
    //builder.printAtomConfs(false, true );
    //exit(0);
    return builder.exportAtoms(); 
}

int initFittedAdata( Atoms* atoms, AddedData* adata, bool byType=true ) {
    //printf("initFittedAdata() isamp=%3i: ntypes=%3i natoms=%3i \n", samples.size(), fittedTypes.size(), atoms->natoms );
    int nfound  =  0;
    int nfound0 = -1;
    //for(int i=0; i<fittedTypes.size(); i++){ printf("initFittedAdata(): fittedTypes[%i]: %i %s \n", i, fittedTypes[i], params->atypes[i].name ); }
    if( adata ){
        //printf("initFittedAdata(): nepairs=%i \n", adata->nep );
        for(int i=0; i<adata->nep; i++){ 
            //printf("initFittedAdata(): adata->bs[%i]: %i %i \n", i, adata->bs[i].x, adata->bs[i].y );
            Vec2i b=adata->bs[i]; adata->host[b.y] = b.x; 
        } // (host_atom_index, epair_index)
    }
    for(int ia=0; ia<atoms->natoms; ia++ ) {
        int ityp   = atoms->atypes[ia];
        Quat4d REQ = typeREQs[ ityp ];
        if( adata ){
            //printf( "initFittedAdata(): ia=%3i ityp=%3i %-8s fitted: %i \n", ia, ityp, params->atypes[ityp].name, fittedTypes[ityp] );
            //adata->host[ia] = adata->bs[ia].x;
            //adata->type[ia] = ityp;
            adata->REQs[ia] = REQ;
        }
        if( byType ){
            if( !fittedTypes[ityp] ) continue; 
        }else{
            if( fabs(REQ.w) < 1e-300 ) continue;
        }
        if( adata ){
           // DEBUG
            adata->HBtoFit    [ia    ] = nfound;
            adata->HBatomsInd [nfound] = ia;
            // adata->HBtoFit    [ia    ] = nfound;
            // adata->HBatomsInd [nfound] = ia;
            // adata->HBatomsHost[nfound] = adata->bs[ia].x;
            // adata->HBatomsType[nfound] = ityp;
            // adata->HBatomsPos [nfound] = atoms->apos[ia];
            // adata->HBatomsREQs[nfound] = REQ;
        }
        if( ia < atoms->n0 ) nfound0 = nfound;
        //printf( "initFittedAdata(): ia=%3i nfound=%i ityp=%3i %-8s  fitted: %i \n", ia, nfound, ityp, params->atypes[ityp].name, fittedTypes[ityp] );
        nfound++;
    }
    if( adata ){
        adata->HBna = nfound;
        adata->HBn0 = nfound0+1;
    }
    return nfound;
}

void printSampleFittedAtoms(int isamp){
    // prints only the fitted atoms from the sample
    //printf("printFittedAdata():\n");
    Atoms* atoms     = samples[isamp];
    AddedData* adata = (AddedData*) atoms->userData;
    int nFit = adata->HBna;
    int n0   = adata->HBn0;
    printf("printFittedAdata(isamp=%i): (natoms=%i,n0=%i) (HBna=%i,HBn0=%i) \n", isamp,   atoms->natoms, atoms->n0, nFit, n0 );
    for(int i=0; i<nFit; i++){
        int ia = adata->HBatomsInd[i];
        int it = atoms->atypes[ia];
        Quat4d REQ = adata->REQs[ia];
        //int it = adata->HBatomsType[i];
        //Quat4d REQ = adata->HBatomsREQs[i];
        printf("fit_atom: %3i ia: %3i frag: %i it: %3i %-8s REQ: %10.3e %10.3e %10.3e %10.3e \n", i, ia,  i>=n0, it, params->atypes[it].name, REQ.x, REQ.y, REQ.z, REQ.w );
    }
}

void printSampleFitSplit(int isamp){
    // print sample atoms showing which atoms are fitted
    //printf("printFittedAdata():\n");
    Atoms* atoms     = samples[isamp];
    AddedData* adata = (AddedData*) atoms->userData;
    int na   = atoms->natoms; 
    int n0   = atoms->n0;
    printf("printSampleFitSplit(isamp=%i): (natoms=%i,n0=%i) (HBna=%i,HBn0=%i) \n", isamp,   na, n0, adata->HBna, adata->HBn0 );
    for(int ia=0; ia<na; ia++){
        int it = atoms->atypes[ia];
        Quat4d REQ = adata->REQs[ia];
        printf("atom: ia: %3i frag: %i fitted: %3i it: %3i %-8s REQ: %10.3e %10.3e %10.3e %10.3e \n", ia, ia>=n0, adata->HBtoFit[ia], it, params->atypes[it].name, REQ.x, REQ.y, REQ.z, REQ.w );
    }
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
    //printf( "addAndReorderEpairs() samples.size()=%i \n", samples.size() );
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

    //int ntypesTot = params->atypes.size();
    //bool fittedTypes[ ntypesTot ];

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

    // // Store electron pair relationships in AddedData
    AddedData* data = new AddedData();
    atoms->userData = data;
    data->nep       = nep_found;
    data->bs        = bs;       //  bs[i].x is epair index, bs[i].y is host atom index
    data->dirs      = dirs;   
    _realloc0( data->host, atoms->natoms, -1 );
    for(int i=0; i<nep_found; i++){ data->host[bs[i].y] = bs[i].x; } // bs[i].x is host atom index, bs[i].y is epair index
    //data->isep      = isep;

    if(bEvalOnlyCorrections){
        int naFit = initFittedAdata( atoms, 0 );  // first pass just count the number of fitted atoms
        data->reallocHB( atoms->natoms, naFit );  // then we can allocate the arrays
        initFittedAdata( atoms, data );           // second pass store the fitted atoms
    }
    
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
            //if(bEvalOnlyCorrections){ printFittedAdata( samples.size()-1 ); }
            il=0; nbatch++;
        }
    }
    if(fout)fclose(fout);
    fclose(fin);
    //init_types();
    countTypesPresent( );
    printTypeParams( true );
    return nbatch;
}

int export_Erefs( double* Erefs ){
    if(Erefs){ for(int i=0; i<nbatch; i++){ Erefs[i] = samples[i]->Energy; } }
    return nbatch;
}

void DOFsToTypes(){
    //printf( "DOFsToTypes() nDOFs: %i @DOFtoTyp=%p \n", nDOFs, DOFtoTyp );
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = DOFtoTyp[i];
        //printf( "DOFsToType()[%3i] rt(%3i|%i) %8s.%c  %g \n", i, rt.x,rt.y,   params->atypes[rt.x].name,"REQH"[rt.y], DOFs[i] );
        typeREQs[rt.x].array[rt.y] = DOFs[i];
    }
}

void DOFsFromTypes(){
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = DOFtoTyp[i];
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
            double lep = Lepairs;
            if( bEpairDistByType ){ typeREQs[atoms->atypes[iE]].w; }
            apos[iE] = apos[iX] + ad->dirs[j] * lep;  // We move the electron pair to proper distance from the atom
            //isEp[iE] = 1;
        }
    }
}

__attribute__((hot)) 
double evalSample( int isamp, const Atoms* atoms, double wi, Quat4d* fREQs ) const {
    //double wi   = (weights)? weights[isamp] : 1.0; 
    const AddedData* adata = (const AddedData*)(atoms->userData);
    alignas(32) double Qs  [atoms->natoms];
    alignas(32) Vec3d  apos[atoms->natoms];   // atomic positions
    fillTempArrays( atoms, apos, Qs );
    int     nj = atoms->n0;
    int     j0 = 0; 
    int     ni = atoms->natoms - atoms->n0;
    int     i0 = atoms->n0;
    for(int i=0; i<atoms->natoms; i++){ fREQs[i]=Quat4dZero; }
    double E=0;
    bool bJ = bEvalJ && ( !bWriteJ );
    switch (imodel){
        case 0:{ 
            E =   evalExampleDerivs( funcVar_LJQH2, i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs( funcVar_LJQH2, j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 2
        }break;
        case 1:{ 
            E =   evalExampleDerivs_LJQH2   ( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs_LJQH2   ( j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 2
        }break;
        case 2:{ 
            E =   evalExampleDerivs_MorseQH2( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs_MorseQH2( j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 2
        }break;
        case 3:{ 
            E =   evalExampleDerivs_LJQH2_SR   ( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, adata->host, fREQs );    // variational derivatives on molecule 1
            if(bJ)evalExampleDerivs_LJQH2_SR   ( j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, adata->host, fREQs );    // variational derivatives on molecule 2
        }break;
        //case 4:{ 
        //    E =   evalExampleDerivs_MorseQH2_SR( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 1
        //    if(bJ)evalExampleDerivs_MorseQH2_SR( j0, nj, i0, ni, atoms->atypes, apos, typeREQs, Qs, fREQs );    // variational derivatives on molecule 2
        //}break;
    }
    return E;
}

double evalSample_corr( int isamp, Quat4d* fREQs ) const {
    printf( "evalSample_corr() \n" );
    Atoms*           atoms = samples[isamp];
    const AddedData* adata = (const AddedData*)(atoms->userData);
    int     nj = adata->HBn0;
    int     j0 = 0; 
    int     ni = adata->HBna - adata->HBn0;
    int     i0 = adata->HBn0;

    int     nj_ = atoms->n0;
    int     j0_ = 0; 
    int     ni_ = atoms->natoms - atoms->n0;
    int     i0_ = atoms->n0;

    for(int i=0; i<adata->HBna; i++){ fREQs[i]=Quat4dZero; }
    double E=0;
    bool bJ = bEvalJ && ( !bWriteJ );
    switch (imodel){
        case 1:{ 
            // double evalExampleDerivs_LJQH2_corr( int i0, int ni, int j0, int nj, Vec3d* ps, Quat4d* REQs, int* host, Quat4d* dEdREQs ) const {
            // NOTE: there may be problem with double-counting of Energy, if both i and j atom are fitted
            E = evalExampleDerivs_LJQH2_corr( i0, ni, j0_, nj_, adata->HBatomsInd, atoms->apos, adata->REQs, adata->host, adata->HBtoFit, fREQs, false );    // variational derivatives on molecule 1
            +   evalExampleDerivs_LJQH2_corr( j0, nj, i0_, ni_, adata->HBatomsInd, atoms->apos, adata->REQs, adata->host, adata->HBtoFit, fREQs, true  );    // variational derivatives on molecule 2
        }break;
    }
    return E;
}

double evalSample_uncorr( int isamp ) const {
    printf( "evalSample_uncorr() \n" );
    Atoms* atoms     = samples[isamp];
    AddedData* adata = (AddedData*)atoms->userData;
    int     nj = atoms->n0;
    int     j0 = 0; 
    int     ni = atoms->natoms - atoms->n0;
    int     i0 = atoms->n0;
    double  E=0;
    switch (imodel){
        case 1:{ 
            E = evalExampleDerivs_LJQH2_uncorr( i0, ni, j0, nj, atoms->apos, adata->REQs, adata->host, adata->HBtoFit ); 
        }break;
    }
    adata->Emodel0 = E;
    return E;
}

void evalSamples_uncorr( ){ 
    for(int isamp=0; isamp<samples.size(); isamp++){evalSample_uncorr( isamp );}
}

__attribute__((hot)) 
void smoothWeight( double E, double& wi )const{
    switch (iWeightModel){
        case 1:{ wi = linstep_down   ( E, EmodelCutStart, EmodelCut);} break;
        case 2:{ wi = smoothstep_down( E, EmodelCutStart, EmodelCut);} break;
    }
}

__attribute__((hot)) 
void handleOverRepulsive( int isamp, double E, const Atoms* atoms, double wi ){
    if(bListOverRepulsive   ){ if(E>EmodelCut) overRepulsiveList.push_back( isamp  );  };
    if(bPrintOverRepulsive || bSaveOverRepulsive ){
        char tmp[256];
        sprintf(tmp, "sample[%4i] Eref: %16.6e Emodel: %16.6e wi: %16.6e EmodelCutStart,EmodelCut: %16.6e %16.6e", isamp, atoms->Energy, E, wi, EmodelCutStart, EmodelCut );
        //if(bSaveOverRepulsive ){ atoms->saveXYZ( fname_overRepulsive, "a", true, true, {1,1,1}, tmp, true ); };
        //if(bSaveOverRepulsive ){ params->saveXYZ( fname_overRepulsive, "a", true, true, {1,1,1}, tmp, true ); };
        if(bSaveOverRepulsive ){ saveDebugXYZ( 0, atoms->natoms, atoms->atypes, atoms->apos, fname_overRepulsive, tmp );}
        if(bPrintOverRepulsive){ printf( "handleOverRepulsive() skipped %s \n", tmp ); }
    }    
}

void printOverRepulsiveList(){
    printf( "printOverRepulsiveList() nfound=%i: { ", overRepulsiveList.size() );
    for(auto i : overRepulsiveList){ printf( "%i,", i ); }  
    printf( "}\n" );  
}

__attribute__((hot)) 
double evalSampleError( int isamp, double& E ){
    //isamp_debug = i;
    Atoms* atoms  = samples[isamp];
    double wi     = (weights)? weights[isamp] : 1.0;
    if(wi<-1e-300) return 0;
    alignas(32) Quat4d fREQs [atoms->natoms];
    alignas(32) double fDOFs_[nDOFs];
    E = evalSample( isamp, atoms, wi, fREQs );
    //if(verbosity>3)
    //printf( "evalSampleError() isamp: %3i Emodel: %20.6f Eref: %20.6f bBroadcastFDOFs=%i @sample_fdofs=%p \n", isamp, E, atoms->Energy, bBroadcastFDOFs, sample_fdofs );
    // if( E>EmodelCutStart ){ 
    //     if(bWeightByEmodel){ smoothWeight( E, wi ); };
    //     handleOverRepulsive( isamp, E, atoms, wi );  
    //     if( weights ) weights[isamp] = wi;   // store to weights so that we know 
    //     if( bDiscardOverRepulsive && (E>EmodelCut) ){ E=NAN; return 0; }  
    // }
    if( bSaveSampleToXYZ ){
        char comment[256];
        sprintf(comment, "# %4i Eref: %16.6e Emodel: %16.6e wi: %16.6e", isamp, atoms->Energy, E, wi );
        //printf( "evalSampleError() saving %s comment: %s \n", xyz_out, comment );
        saveDebugXYZ( 0, atoms->natoms, atoms->atypes, atoms->apos, xyz_out, comment );
    }
    double Eref    = atoms->Energy;
    double dE      = E - Eref;
    wi*= invWsum;
    double dEw     = 2.0*dE*wi;
    double Error   =  dE*dE*wi;
    //printf( "evalSampleError() isamp: %3i Emodel: %20.6f Eref: %20.6f bBroadcastFDOFs=%i @sample_fdofs=%p \n", isamp, E, atoms->Energy, bBroadcastFDOFs, sample_fdofs );
    double* fDOFs__ = bBroadcastFDOFs ? sample_fdofs + isamp*nDOFs : fDOFs_;   // broadcast fDOFs ?
    for(int k=0; k<nDOFs; k++){ fDOFs__[k]=0; }                                // clean fDOFs
    acumDerivs( atoms->natoms, atoms->atypes, dEw, fREQs, fDOFs_ );            // accumulate fDOFs from fREQs
    if( bEpairs ){
        AddedData * ad = (AddedData*)atoms->userData;
        acumHostDerivs( ad->nep, ad->bs, atoms->atypes, dEw, fREQs, fDOFs_ );
    }
    if( bUdateDOFbounds  ){ updateDOFbounds( fDOFs__ ); }
    
    if( !bBroadcastFDOFs ){ 
        double F2=0.0;
        for(int k=0; k<nDOFs; k++){ double fi=fDOFs__[k]; fDOFs[k] += fi; F2+=fi*fi;  } 
        //printf( "evalSampleError() isamp: %3i Eref: %20.6f Emodel: %20.6f |F|= f: %20.6f  wi: %20.6f dEw: %20.6f \n", isamp, atoms->Energy, E, sqrt(F2), wi, dEw );
    }
    return Error;
}

__attribute__((hot)) 
double evalSampleError_corr( int isamp, double& E ){
    printf( "evalSampleError_corr() \n" );
    //isamp_debug = i;
    Atoms* atoms  = samples[isamp];
    const AddedData* adata = (const AddedData*)(atoms->userData);
    double wi     = (weights)? weights[isamp] : 1.0;
    //if(wi<1e-300) return 0;
    alignas(32) Quat4d fREQs  [adata->HBna];
    E = evalSample_corr( isamp, fREQs );
    //if(bWeightByEmodel){ if(E>EmodelCut){ return 0; } smoothWeight(E,wi); }
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
    //acumDerivs_corr( adata->HBna, adata->HBatomsType, adata->HBatomsHost, dEw, fREQs, fDOFs_ );
    acumDerivs_corr( adata->HBna, adata->HBatomsInd, atoms->atypes, adata->host, dEw, fREQs, fDOFs_ );
    return Error;
}

// =========== Correction Only optimized versions of evalSample and evalSampleError

void clear_fDOFbounds(){
    for(int i=0; i<nDOFs; i++){ fDOFbounds[i] = Vec2d{1.e+300,-1.e+300}; }
}

/**
 * @brief Evaluates the variational derivatives of the fitting error (sumed batch training samples) with respect to all fitting parameters (i.e. non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) for each atom type). The atomic systems are not assumed rigid (i.e. the atomic positions can vary freely from sample to sample).
 * 
 * @param Eout array to store the non-covalent interaction energy values of each atomic system in the batch. if Eout==null, the function will not store the energy values.
 * @return double, returns the total fitting error.
 */
__attribute__((hot)) 
double evalSamples_noOmp( double* Eout=0 ){ 
    if( bListOverRepulsive   ) overRepulsiveList.clear();
    if( bSaveOverRepulsive   ) clearFile( fname_overRepulsive );
    if( bClearDOFboundsEpoch ) clear_fDOFbounds(); // clean fDOFbounds
    bBroadcastFDOFs    = false; 
    double Error = 0.0;
    int nsamp = samples.size();
    //printf( "evalSamples_noOmp() nsamp=%i\n", nsamp );
    for(int isamp=0; isamp<nsamp; isamp++){
        double E; 
        Error+=evalSampleError( isamp, E );
        //if(bEvalOnlyCorrections){ ((AddedData*)samples[isamp]->userData)->Emodel0 = E; }
        if(Eout){
            //if( bListOverRepulsive ){ if(overRepulsiveList.back()==isamp) E=NAN; }
            //printf( "evalSamples_noOmp() isamp: %3i Emodel: %20.6f Eref: %20.6f @Eout=%p \n", isamp, E, samples[isamp]->Energy, Eout );
            Eout[isamp]=E;
        }
    }
    if(bListOverRepulsive){ printOverRepulsiveList(); }
    return Error;
}

__attribute__((hot)) 
double evalSamples_omp( double* Eout=0 ){ 
    bListOverRepulsive = false;
    bBroadcastFDOFs    = true; 
    double Error       = 0.0;
    int nsamp = samples.size();
    printf( "evalSamples_omp() nsamp=%i\n", nsamp );
    #pragma omp parallel for reduction(+:Error)
    for(int isamp=0; isamp<nsamp; isamp++){
       double E; 
       Error+=evalSampleError( isamp, E );
       if(Eout)Eout[isamp]=E;
    }
    if(bBroadcastFDOFs){ reduce_sample_fdofs(); }
    return Error;
}

__attribute__((hot)) 
double evalSamples_corr( double* Eout=0 ){ 
    bBroadcastFDOFs = false; 
    double Error    = 0.0;
    int nsamp       = samples.size();
    for(int isamp=0; isamp<nsamp; isamp++){
        double E; 
        Error+=evalSampleError_corr( isamp, E );
        if(Eout){ Eout[isamp]=E;}
    }
    if(bListOverRepulsive){ printOverRepulsiveList(); }
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
void acumDerivs_omp_atomic( int n, int* types, double dEw, Quat4d* fREQs ){
    for(int i=0; i<n; i++){
        int t            = types[i];     // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];  // get index of degrees of freedom for atom type t
        const Quat4d  f  = fREQs[i];
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
void acumHostDerivs_omp_atomic( int nepair, Vec2i* epairAndHostIndex, int* types, double dEw, Quat4d* fREQs  ){
    for(int i=0; i<nepair; i++){
        Vec2i ab         = epairAndHostIndex[i];
        int t            = types[ab.i];      // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];      // get index of degrees of freedom for atom type t
        const Quat4d&  f = fREQs[ab.j];         // get variation from the host atom
        if(tt.z>=0){ 
            #pragma omp atomic
            fDOFs[tt.z]  -= f.z*dEw; } // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
    }
}

__attribute__((hot)) 
void acumDerivs( int n, int* types, double dEw, Quat4d* fREQs, double* fDOFs ){
    for(int ia=0; ia<n; ia++){           // loop over all atoms in the sample
        int t            = types[ia];    // [natom] map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];  // [ntype] get DOF index for type t
        const Quat4d  f  = fREQs[ia];    // [natom] variational derivative of parameters of atom i
        if(tt.x>=0){ fDOFs[tt.x]+=f.x*dEw; }
        if(tt.y>=0){ fDOFs[tt.y]+=f.y*dEw; }
        if(tt.z>=0){ fDOFs[tt.z]+=f.z*dEw; }
        if(tt.w>=0){ fDOFs[tt.w]+=f.w*dEw; }
        //printf( "acumDerivs() ia: %3i  %2i|%-8s  dEw= %g fREQ( %10.2e %10.2e %10.2e %10.2e ) tt(%2i,%2i,%2i,%2i)\n", ia, t, params->atypes[t].name, dEw, f.x,f.y,f.z,f.w, tt.x,tt.y,tt.z,tt.w );
        //if((tt.y>=0)&&(i==0))printf( "acumDerivs i= %i f= %g f*dEw= %g fDOFs= %g\n", i, f.y, f.y*dEw, fDOFs[tt.y] );
        //if((tt.x>=0)||(tt.y>=0))printf( "acumDerivs i= %i dE= %g f.x= %g fDOFs= %g f.y= %g fDOFs= %g\n", i, dE, f.x, fDOFs[tt.x], f.y, fDOFs[tt.y] );
    }
    //exit(0);
}

__attribute__((hot)) 
void acumHostDerivs( int nepair, Vec2i* epairAndHostIndex, int* types, double dEw, Quat4d* fREQs, double* fDOFs  ){
    for(int i=0; i<nepair; i++){
        Vec2i ab         = epairAndHostIndex[i];
        int t            = types[ab.i];      // map atom index i to atom type t
        const Quat4i& tt = typToREQ[t];      // get index of degrees of freedom for atom type t
        const Quat4d&  f = fREQs[ab.j];         // get variation from the host atom
        if(tt.z>=0){ fDOFs[tt.z]  -= f.z*dEw; } // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
    }
}


__attribute__((hot)) 
void acumDerivs_corr( int n, int* aInd, int* types, int* host, double dEw, Quat4d* fREQs, double* fDOFs ){
    for(int i=0; i<n; i++){
        int ia           = aInd[i];
        int t            = types[ia];    // map atom index i to atom type t
        int ih           = host [ia];
        const Quat4i& tt = typToREQ[t];  // get index of degrees of freedom for atom type t
        const Quat4d  f  = fREQs[i];
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

void saveDebugXYZ( int i0, int n, const int* types, const Vec3d* ps, const char* fname="debug.xyz", const char* comment="" )const{
    FILE* fout = fopen(fname,"a");
    fprintf(fout, "%i\n", n );
    fprintf(fout, "%s\n", comment );
    for(int ii=0; ii<n; ii++){ // loop over all atoms[i] in system
        const int      i    = i0+ii;
        const Vec3d&  pi    = ps      [i ]; 
        //const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        //printf( "atom[%i] ti=%i\n", ii, ti );
        fprintf(fout, "%c %7.3f %7.3f %7.3f \n", params->atypes[ti].name[0], pi.x, pi.y, pi.z );
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
    //printf( "checkSampleRepulsion(isamp=%3i)[%3i,%3i] (%8s,%8s) r: %20.10f ELJ:  %20.10f \n", isamp_debug, i,j, params->atypes[ti].name, params->atypes[tj].name , r, ELJ  ); 
    if(bPrint)printf( "%s\n", comment );
    if(bSave)samples[isamp_debug]->saveXYZ( "debug.xyz", "a", true, true, Vec3i{1,1,1}, comment );
    //saveDebugXYZ( 0, ni+nj, types, ps, typeREQs, Qs, "debug.xyz", comment );
    //iBadFound++; if(iBadFound>=nBadFoundMax){ printf("ERROR in evalExampleDerivs_LJQH2(): too many bad pairs (%i) \n", iBadFound); exit(0); }
    return true;
}

// Optimized versions of evaluation functions
__attribute__((hot)) 
static double evalExampleDerivs_LJQH2_uncorr( int i0, int ni, int j0, int nj, Vec3d* ps, Quat4d* REQs, int* host, int* fitted ){
    // this function compute energy of pairs of atoms which are not fitted, and  not variational derivatives, i.e. they are not included in the optimization
    printf("evalExampleDerivs_LJQH2_uncorr() i0: %i ni: %i j0: %i nj: %i \n", i0, ni, j0, nj );
    double Etot = 0.0;
    for(int ii=0; ii<ni; ii++) {
        const int     i     = i0+ii;
        const Vec3d&  pi    = ps[i];
        const int     ih    = host[i];
        const bool    bEpi  = ih>=0;  // is this atom an electron pair?
        const Quat4d& REQi  = REQs[ii];
        //Quat4d        fREQi = Quat4dZero;
        const int     ommit_i  = (fitted[i]>=0);
        for(int jj=0; jj<nj; jj++) {
            const int j = j0+jj;
            //printf( "evalExampleDerivs_LJQH2_uncorr() i: %3i j: %3i  ommit_i: %i fitted[j]: %i\n", i, j, ommit_i, fitted[j] );
            if( ommit_i || ( fitted[j]>=0) ) continue;  // if atom i or j is fitted, skip
            const Quat4d& REQj = REQs[j];
            double H         = REQi.w * REQj.w;   
            const double sH      = (H<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H *= sH;
            int jh           = host[j];
            const double R0  = REQi.x + REQj.x;
            const double E0  = REQi.y * REQj.y;
            const double Q   = REQi.z + REQj.z;
            const Vec3d  dij = ps[j] - pi;
            const double r       = dij.norm();
            const double ir      = 1/r;
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;
            double Eij = Eel;
            double dE_dH;
            double dE_dw;
            if       (bEpi ){ // i is an electron pair
                Eij  += getEpairAtom( r, H, REQi.x, dE_dH, dE_dw );
            }else if (jh>=0){ // j is an electron pair
                double dE_dH;
                Eij  += getEpairAtom( r, H, REQj.x, dE_dH, dE_dw ); // Note - radius of electron REQi.x is used for decay constant
            }else{           // both i and j are real atoms
                const double u       = R0/r;
                const double u3      = u*u*u;
                const double u6      = u3*u3;
                const double u6p     = (1.0 + H) * u6;                     
                const double dE_dE0  = u6 * (u6p - 2.0);
                Eij                 +=  E0 * dE_dE0;
            }
            Etot += Eij;
            printf( "evalExampleDerivs_LJQH2_uncorr() i: %3i j: %3i Eij: %10.3e     fitted[i]: %3i fitted[j]: %3i \n", i, j, Eij, fitted[i], fitted[j] );
            //printf( "evalExampleDerivs_LJQH2_uncorr() i: %3i j: %3i H: %10.3e R0: %10.3e E0: %10.3e Q: %10.3e Eel: %10.3e Eij: %10.3e\n", i, j, R0,E0,Q,H,  Eij );
        }
    }
    //exit(0);
    return Etot;
}


// Optimized versions of evaluation functions
__attribute__((hot)) 
static double evalExampleDerivs_LJQH2_corr( int i0, int ni, int j0, int nj, int* aInd, Vec3d* ps, Quat4d* REQs, int* host, int* fitted, Quat4d* dEdREQs, bool bCheckAddJ ){
    printf("evalExampleDerivs_LJQH2_corr() i0: %i ni: %i j0: %i nj: %i \n", i0, ni, j0, nj );
    double Etot = 0.0;
    for(int ii=0; ii<ni; ii++) {
        int           i      = i0+ii;
        const int     ia     = aInd[i];
        //printf( "evalExampleDerivs_LJQH2_corr() i: %3i ia: %3i \n", ii, ia );
        const Vec3d&  pi    = ps[ia];
        const int     ih    = host[ia];
        const bool    bEpi  = ih>=0;  // is this atom an electron pair?
        const Quat4d& REQi  = REQs[ia];
        Quat4d        fREQi = Quat4dZero;
        for(int jj=0; jj<nj; jj++) {
            const int j = j0+jj;
            //if( fitted[j]<0 ) continue;     // evaluate only fitted atoms 
            const Quat4d& REQj = REQs[j];
            double H         = REQi.w * REQj.w;   
            const double sH  = (H<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H *= sH;
            int jh           = host[j];
            const double R0  = REQi.x + REQj.x;
            const double E0  = REQi.y * REQj.y;
            const double Q   = REQi.z + REQj.z;
            const Vec3d  dij = ps[j] - pi;
            const double r       = dij.norm();
            const double ir      = 1/r;
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;
            double Eij = Eel;
            double dE_dw;
            if       (bEpi ){ // i is an electron pair
                double dE_dH;    
                Eij     += getEpairAtom( r, H, REQi.x, dE_dH, dE_dw );
                fREQi.z  += dE_dQ * REQj.z;  // dEtot/dH2i
                fREQi.w  += dE_dH * REQj.w;  // dEtot/dH2i
            }else if (jh>=0){ // j is an electron pair
                double dE_dH;
                Eij     += getEpairAtom( r, H, REQj.x, dE_dH, dE_dw ); // Note - radius of electron REQi.x is used for decay constant
                fREQi.w  += dE_dH * REQi.w;  // dEtot/dH2i
            }else{           // both i and j are real atoms
                const double u       = R0/r;
                const double u3      = u*u*u;
                const double u6      = u3*u3;
                const double u6p     = (1.0 + H) * u6;                     
                const double dE_dE0  = u6 * (u6p - 2.0);
                const double dE_dH   = -E0 * u6 * u6;
                Eij                 +=  E0 * dE_dE0;
                fREQi.w             +=  dE_dH * REQj.w * sH;     // dEtot/dH2i
            }
            if( bCheckAddJ ){  if( fitted[j]>=0 ){ Eij=0.0; } }  // We must check if j is also fitted atom to avoid double counting
            printf( "evalExampleDerivs_LJQH2_corr() i: %3i j: %3i  Eij: %10.3e       fitted[i]: %3i fitted[j]: %3i \n", ia, j,   Eij,  fitted[ia], fitted[j]  );
            //printf( "evalExampleDerivs_LJQH2_corr() i: %3i j: %3i H: %10.3e R0: %10.3e E0: %10.3e Q: %10.3e Eel: %10.3e Eij: %10.3e\n", i, j, R0,E0,Q,H,  Eij );
            Etot += Eij;
        }
        if(dEdREQs)dEdREQs[i].add(fREQi);
    }
    return Etot;
}
    
__attribute__((hot)) 
double evalExampleDerivs_LJQH2_SR( int i0, int ni, int j0, int nj, int*  types, Vec3d* ps, Quat4d* typeREQs, double* Qs, int* host, Quat4d* dEdREQs )const{
    //printf("evalExampleDerivs_LJQH2() i0: %i ni: %i j0: %i nj: %i \n", i0, ni, j0, nj );
    // this differes form evalExampleDerivs_LJQH2() in that it uses short range corrections for electron pairs
    double Etot = 0.0;
    const bool bWJ = bWriteJ&&dEdREQs;
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system
        const int      i    = i0+ii;
        const int     ih    = host[i];
        const bool    bEpi  = ih>=0;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;

        //if( bEpi ){ Vec3d dr = ps[ih]-pi; printf( "i ih |dij|=%g \n", i, ih, dr.norm() );  }
        for(int jj=0; jj<nj; jj++){ 
            const int   j        = j0+jj;
            const int   jh       = host[j];
            const double     Qj  = Qs[j];
            const Vec3d      dij = ps[j] - pi;
            const int        tj  = types[j];
            const Quat4d& REQj   = typeREQs[tj];
            const double R0      = REQi.x + REQj.x;
            const double E0      = REQi.y * REQj.y; 
            const double Q       = Qi     * Qj    ;
            double       H       = REQi.w * REQj.w;     // expected H2<0
            const double sH      = (H<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H *= sH;
            // --- Electrostatic
            const double r       = dij.norm();
            const double ir      = 1/r;
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;

            double Eij    = Eel;
            double dE_dR0 = 0.0;
            double dE_dE0 = 0.0;
            double dE_dH  = 0.0;
            double dE_dw  = 0.0;
            const bool bEpj = jh>=0;
            if( bEpi||bEpj ){     // i or j is an electron pair
                // double w=0; 
                // if( bEpj       ){ w += REQj.x; }
                // if( bEpi       ){ w += REQi.x; }
                // if( bEpi&&bEpi ){ w *= 0.5;    }
                // Eij += getEpairAtom( r, H, w, dE_dH, dE_dw );
            }else{
                const double u   = R0/r;
                const double u3  = u*u*u;
                const double u6  = u3*u3;
                const double u6p = ( 1.0 + H ) * u6;                     
                dE_dE0           = u6 * ( u6p - 2.0 );
                dE_dR0           = 12.0 * (E0/R0) * u6 * ( u6p - 1.0 );
                dE_dH            = -E0 * u6 * u6;
                Eij              += E0 * dE_dE0;
            }

            if( bWJ ){ dEdREQs[j].add( Quat4d{
                        dE_dR0,                    // dEtot/dR0_j
                       -dE_dE0  * 0.5 * REQi.y,    // dEtot/dE0_j
                        dE_dQ   * Qi,              // dEtot/dQ_j
                        dE_dH   * REQi.w * sH,     // dEtot/dH2j
            }); }

            //if(bCheckRepulsion)[[unlikely]]{ checkSampleRepulsion( ELJ, i,j, ti,tj, r, true, true ); }
            // --- Energy and forces
            Etot    +=  Eij;
            fREQi.x +=  dE_dR0;                    // dEtot/dR0_i
            fREQi.y +=  dE_dE0  * 0.5 * REQj.y;    // dEtot/dE0_i
            fREQi.z += -dE_dQ   * Qj;              // dEtot/dQ_i
            fREQi.w +=  dE_dH   * REQj.w * sH;     // dEtot/dH2i
        }
        if(dEdREQs)dEdREQs[i].add(fREQi);
    }
    // printAtomParamDerivs( ni+nj, dEdREQs, isamp_debug );
    //printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

__attribute__((hot)) 
//double evalExampleDerivs_LJQH2( int i0, int ni, int j0, int nj, int*  types, Vec3d* ps, Quat4d*  typeREQs, double* Qs, Quat4d* dEdREQs )const{
double evalExampleDerivs_LJQH2( int i0, int ni, int j0, int nj, int* __restrict__ types, Vec3d* __restrict__ ps, Quat4d* __restrict__ typeREQs, double* __restrict__ Qs, Quat4d* __restrict__ dEdREQs )const{
    //printf("evalExampleDerivs_LJQH2() i0: %i ni: %i j0: %i nj: %i \n", i0, ni, j0, nj );
    double Etot = 0.0;
    const bool bWJ = bWriteJ&&dEdREQs;
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
            double       H       = REQi.w * REQj.w;      // expected H2<0
            const double sH      = (H<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H *= sH;
            // --- Electrostatic
            const double r       = dij.norm();
            const double ir      = 1/r;
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;

            const double u       = R0/r;
            const double u3      = u*u*u;
            const double u6      = u3*u3;
            const double u6p     = ( 1.0 + H ) * u6;                     
            const double dE_dE0  = u6 * ( u6p - 2.0 );
            const double dE_dR0  = 12.0 * (E0/R0) * u6 * ( u6p - 1.0 );
            const double dE_dH  = -E0 * u6 * u6;
            const double ELJ     = E0 * dE_dE0;
            //if(bCheckRepulsion)[[unlikely]]{ checkSampleRepulsion( ELJ, i,j, ti,tj, r, true, true ); }
            // --- Energy and forces
            Etot    +=  ELJ + Eel;
            //printf( "evalExampleDerivs_LJQH2() i: %3i j: %3i  Eij: %10.3e\n", i, j, ELJ + Eel );
            //printf( "evalExampleDerivs_LJQH2() i: %3i j: %3i H: %10.3e R0: %10.3e E0: %10.3e Q: %10.3e Eel: %10.3e Eij: %10.3e\n", i, j, R0,E0,Q,H,  ELJ + Eel );
            //printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) ELJ:  %20.10f   Eel: %20.10f \n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel  );
            //{ int itypPrint=4; if( (ti==itypPrint) || (tj==itypPrint) ){ printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) ELJ,Eel: %12.3e,%12.3e Q(%12.3e|%12.3e,%12.3e) dEdREQH(%12.3e,%12.3e,%12.3e,%12.3e)\n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel, Q,Qi,Qj,  dE_dR0, dE_deps, dE_dQ, dE_dH2  ); } }
            fREQi.x +=  dE_dR0;                    // dEtot/dR0_i
            fREQi.y +=  dE_dE0  * 0.5 * REQj.y;    // dEtot/dE0_i
            fREQi.z +=  -dE_dQ   * Qj;              // dEtot/dQ_i
            fREQi.w +=  dE_dH   * REQj.w * sH;     // dEtot/dH2i
            //printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) dE_dH:  %20.10f   ELJ: %20.10f \n", i,j, params->atypes[ti].name, params->atypes[tj].name, dE_dH, ELJ, sH, REQj.w, REQi.w );
            if( bWJ ){ dEdREQs[j].add( Quat4d{
                        dE_dR0,                    // dEtot/dR0_j
                        -dE_dE0  * 0.5 * REQi.y,    // dEtot/dE0_j
                        dE_dQ   * Qi,              // dEtot/dQ_j
                        dE_dH   * REQi.w * sH,     // dEtot/dH2j
            }); }
        }
        if(dEdREQs)dEdREQs[i].add(fREQi);
    }
    // printAtomParamDerivs( ni+nj, dEdREQs, isamp_debug );
    //printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}


__attribute__((hot)) 
double evalExampleDerivs_MorseQH2( int i0, int ni, int j0, int nj, int* __restrict__ types, Vec3d* __restrict__ ps, Quat4d* __restrict__ typeREQs, double* __restrict__ Qs, Quat4d* __restrict__ dEdREQs )const{
    double Etot = 0.0;
    const bool bWJ = bWriteJ&&dEdREQs;
    const double alpha   = kMorse; 
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
            double       H       = REQi.w * REQj.w;      // expected H2<0
            const double sH      = (H<0.0) ? 1.0 : 0.0; // sH=1.0 if H2<0
            H *= sH;
            // --- Electrostatic
            double r             = dij.norm();
            double ir            = 1/r;    
            const double dE_dQ   = ir * COULOMB_CONST;
            const double Eel     = Q * dE_dQ;

            // --- Morse
            //double alpha   = 6.0 / R0;
            
            double e       = exp( -alpha * ( r - R0 ) );
            double e2      = e * e;
            double e2p     = ( 1.0 + H ) * e2;
            double dE_dE0  = e2p - 2.0 * e;
            double dE_dR0  = 2.0 * alpha * E0 * ( e2p - e );
            double dE_dH   = - E0 * e2;
            double ELJ     =   E0 * dE_dE0;
            //if(bCheckRepulsion)[[unlikely]]{ checkSampleRepulsion( ELJ, i,j, ti,tj, r, true, true ); }
            // --- Energy and forces
            Etot    +=  ELJ + Eel;
            //printf( "evalExampleDerivs_MorseQH2()[%3i,%3i] (%8s,%8s) ELJ:  %20.10f   Eel: %20.10f \n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel  );
            //{ int itypPrint=4; if( (ti==itypPrint) || (tj==itypPrint) ){ printf( "evalExampleDerivs_LJQH2()[%3i,%3i] (%8s,%8s) ELJ,Eel: %12.3e,%12.3e Q(%12.3e|%12.3e,%12.3e) dEdREQH(%12.3e,%12.3e,%12.3e,%12.3e)\n", i,j, params->atypes[ti].name, params->atypes[tj].name , ELJ,Eel, Q,Qi,Qj,  dE_dR0, dE_deps, dE_dQ, dE_dH2  ); } }
            fREQi.x +=  dE_dR0;                    // dEtot/dR0_i
            fREQi.y +=  dE_dE0  * 0.5 * REQj.y;    // dEtot/dE0_i
            fREQi.z +=  -dE_dQ   * Qj;              // dEtot/dQ_i
            fREQi.w +=  dE_dH   * REQj.w * sH;     // dEtot/dH2i
            if( bWJ ){ dEdREQs[j].add( Quat4d{
                        dE_dR0,                    // dEtot/dR0_j
                        -dE_dE0 * 0.5 * REQi.y,    // dEtot/dE0_j
                        dE_dQ   * Qi,              // dEtot/dQ_j
                        dE_dH   * REQi.w * sH,     // dEtot/dH2j
            }); }
            //printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
        if(dEdREQs)dEdREQs[i].add(fREQi);
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
    const bool bWJ = bWriteJ&&dEdREQs;
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
            if(bWJ){ dEdREQs[j].add( Quat4d{
                        fREQH.x,                   // dEtot/dR0_j
                        fREQH.y * 0.5 * REQi.y,    // dEtot/dE0_j
                        fREQH.z   * Qi,            // dEtot/dQ_j
                        fREQH.w  * REQi.w * sH,    // dEtot/dH2j
            }); }
        }
        if(dEdREQs)dEdREQs[i].add(fREQi);
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}


// ======================================
// =========  OPTIMIZE  =================
// ======================================

__attribute__((hot)) 
double evalFitError(int itr, bool bOMP=true, bool bEvalSamples=true){
    double Err=0.0;
    DOFsToTypes(); 
    clean_fDOFs();
    if(bEvalSamples)[[likely]]{
        if(bOMP){  Err = evalSamples_omp();   }
        else    {  Err = evalSamples_noOmp(); }
    }
    if(bPrintBeforReg)printStepDOFinfo( itr, Err, "BEFOR_REG: " );
    if(bRegularize){ 
        double Ereg = regularizeDOFs(); 
        if(bAddRegError)Err += Ereg;
    }   
    if(bPrintAfterReg)printStepDOFinfo( itr, Err, "AFTER_REG: " );
    return Err;
}

__attribute__((hot)) 
double run( int ialg, int nstep, double Fmax, double dt, double max_step, double damping, bool bClamp, bool bOMP ){
    if( verbosity>1){ printf( "FitREQ::run() imodel %i ialg %i nstep %i Fmax %g dt %g max_step %g \n", imodel, ialg, nstep, Fmax, dt, max_step ); }
    if(weights){updateWeightsSum();}
    double Err=0;
    //noiseToDOFs(0.2);
    clear_fDOFbounds();
    if(bOMP){ bBroadcastFDOFs=true; realloc_sample_fdofs();  }
    double F2max=Fmax*Fmax;
    double F2;
    int nsamp = samples.size();
    for(int itr=0; itr<nstep; itr++){
        Err = evalFitError( itr, bOMP );
        switch(ialg){
            case 0: F2 = move_GD         (      dt, max_step          ); break;
            case 1: F2 = move_MD         (      dt, max_step, damping ); break;
            case 2: F2 = move_GD_BB_short( itr, dt, max_step          ); break;
            case 3: F2 = move_GD_BB_long ( itr, dt, max_step          ); break;
        }
        if( bClamp ){ limitDOFs(); } // TODO: should we put this before evalFitError() ?
        if( F2<F2max ){ printf("CONVERGED in %i iterations \n", itr); break; }
    }
    printf("VERY FINAL |E|=%.15g  DOFs= ",Err     ); for(int j=0;j<nDOFs;j++){ printf("%.15g ", DOFs[j]); };printf("\n");
    printf("VERY FINAL |F|=%.15g fDOFs= ",sqrt(F2)); for(int j=0;j<nDOFs;j++){ printf("%.15g ",fDOFs[j]); };printf("\n");
    return Err;
}

__attribute__((hot)) 
double run_omp( int ialg, int nstep, double Fmax, double dt, double max_step, double damping, bool bClamp ){
    double Err=0;
    if( verbosity>1){ printf( "FitREQ::run() nstep %i Fmax %g dt %g isamp %i \n", nstep, Fmax, dt ); }
    if(weights){updateWeightsSum();}
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
            //if(bPrintAfterReg)printStepDOFinfo( itr, Err, "FitREQ::run() BEFORE REGULARIZATION" );
            if(bRegularize){ 
                double Ereg = regularizeDOFs(); 
                if(bAddRegError)Err += Ereg;
            }   
            switch(ialg){
                case 0: F2 = move_GD         (      dt, max_step          ); break;
                case 1: F2 = move_MD         (      dt, max_step, damping ); break;
                case 2: F2 = move_GD_BB_short( itr, dt, max_step          ); break;
                case 3: F2 = move_GD_BB_long ( itr, dt, max_step          ); break;
            }
            if(bClamp     ){ limitDOFs();  }
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

__attribute__((hot)) 
__attribute__((pure))
static inline double constrain( double x, const Vec3d& regX, const Vec3d& regK, double& fout ){
    double E = 0;
    double f = 0;
    if(x<regX.x){
        double d = x      -regX.x; 
        f        =      -d*regK.x;
        E        = 0.5*d*d*regK.x;
    }else if(x>regX.z){
        double d = x      -regX.z; 
        f        =      -d*regK.z;
        E        = 0.5*d*d*regK.z;
    }
    if( regK.y > 0){
        double d  = x-regX.y;
        f        -=       d*regK.y;
        E        += 0.5*d*d*regK.y;
    }
    fout += f;
    return E;
}

void noiseToDOFs(double d){
    for(int i=0; i<nDOFs; i++){ DOFs[i]+=randf(-d,d); }
}


// New regularization function operating per-DOF
__attribute__((hot)) 
double regularizeDOFs(){
    //printf("regularizeDOFs() nDOFs=%i @DOFregX=%p @DOFregX=%p @DOFs=%p @fDOFs=%p @DOFtoTyp=%p \n",  nDOFs, DOFregX, DOFregK, DOFs, fDOFs, DOFtoTyp); 
    double E = 0;
    for(int i=0; i<nDOFs; i++){   
        E += constrain( DOFs[i], DOFregX[i], DOFregK[i], fDOFs[i] );
        //{ Vec2i tc    = DOFtoTyp[i]; printf( "regularizeDOFs() i: %3i   %8s|%c x: %10.3e f: %10.3e E: %10.3e\n", i,  params->atypes[tc.x].name, "REQH"[tc.y],  DOFs[i], fDOFs[i], E );   }
    }
    //printf( "regularizeDOFs() END E: %10.3e\n", E );
    return E;
}

__attribute__((hot)) 
void limitDOFs(){
    for(int i=0; i<nDOFs; i++){   DOFs[i] = _clamp(DOFs[i], DOFlimits[i].x, DOFlimits[i].y); }
    //printf("limitDOFs() DOFs: "); for(int i=0;i<nDOFs;i++){ printf("%8.2e ", DOFs[i]); };printf("\n");
    //printf("limitDOFs() DOFs: "); for(int i=0;i<nDOFs;i++){ printf("%8.2e ", DOFlimits[i].x); };printf("\n");
    //printf("limitDOFs() DOFs: "); for(int i=0;i<nDOFs;i++){  printf("%8.2e ", DOFlimits[i].y); };printf("\n");
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
    for(int i=0; i<nDOFs; i++){ double f=fabs(fDOFs[i]);  if(f>fm){fm=f; ifm=i;}  }
    //for(int i=0; i<nDOFs; i++){double f=fabs(fDOFs[i]/DOFs[i]); if(f>fm){fm=f; ifm=i;}}
    if(dt*fm>max_step){
        double new_dt = max_step/fm;
        //printf( "limit_dt %10.3e -> %10.3e max_fDOF[%3i]=%10.3e  max_fDOF*dt: %10.3e max_step: %8.3e \n", dt, new_dt, ifm, fDOFs[ifm], -fDOFs[ifm]*dt, max_step );
        //printf( "limit_dt %g\n", dt );
        dt=new_dt;
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
double move_GD( double dt, double max_step=-0.1 ){
    //printf("now in move_GD\n");
    double F2 = 0;
    if(max_step>0)dt=limit_dt(dt,max_step);
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        DOFs[i] += f*dt;
        F2      += f*f;
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
    if(max_step>0)dt=limit_dt(dt,max_step);
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
    if(max_step>0)dt=limit_dt(dt,max_step);
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
double move_MD( double dt, double max_step=-0.1, double damp=0.1, bool bClimbBreak=true ){
    double cdamp = 1.0-damp;
    if(bClimbBreak){
        double fv = 0.0; for(int i=0; i<nDOFs; i++){ fv += vDOFs[i]*fDOFs[i]; }
        printf( "move_MD fv= %g\n", fv );
        if(fv<0.0){      for(int i=0; i<nDOFs; i++){       vDOFs[i] = 0.0;    }; printf( "ClimbBreak\n" ); } 
    }
    double F2 = 0.0;
    if(max_step>0.0){ dt=limit_dt_MD(dt,max_step,cdamp); };
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        double v = vDOFs[i];
        v       *= cdamp;
        v       += f*dt;
        vDOFs[i] = v;
        DOFs[i] += v*dt;
        F2      += f*f;
    }
    return F2;
}

/// @brief Cleans the derivative array by setting all its elements to zero.
void clean_fDOFs    ()                { for(int i=0; i<nDOFs; i++){ fDOFs[i]=0;                          }}
void updateDOFbounds( double* fDOFs_ ){ for(int i=0; i<nDOFs; i++){ fDOFbounds[i].enclose( fDOFs_[i] );  }}

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
        //int nsamp = samples.size();
        double F2 = 0.0;                           for(int i=0; i<nDOFs; i++){ F2 += fDOFs[i]*fDOFs[i]; }
        double F = sqrt(F2);
        if(bPrintDOFs ){ printf("%s step: %5i |E|: %8.3e |F|: %8.3e  DOFs: ",label, istep, Err, F); for(int j=0; j<nDOFs; j++){ printf("%10.2e ", DOFs[j]); }; printf("\n"); }
        if(bPrintfDOFs){ printf("%s step: %5i |E|: %8.3e |F|: %8.3e fDOFs: ",label, istep, Err, F); for(int j=0; j<nDOFs; j++){ printf("%10.2e ",fDOFs[j]); }; printf("\n"); }
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
    printf("printDOFmapping()\n");
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
}

void printDOFsToTypes() const {
    printf("printDOFsToTypes() nDOFs=%i\n", nDOFs);
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = DOFtoTyp[i];
        const char comp = "REQH"[rt.y];
        printf("%3i -> %3i|%i  %-8s %c \n", i, rt.x, rt.y,  params->atypes[rt.x].name, comp);
    }
    //printf("\n");
}

void printTypesToDOFs() const {
    printf("printTypesToDOFs() ntype=%i\n", ntype);
    for(int i=0; i<ntype; i++){
        const Quat4i& tt = typToREQ[i];
        if(tt.x>=0 || tt.y>=0 || tt.z>=0 || tt.w>=0){
            printf("%-8s %3i { ", params->atypes[i].name, i );
            if(tt.x>=0) printf("R=%d ", tt.x);
            if(tt.y>=0) printf("E=%d ", tt.y);
            if(tt.z>=0) printf("Q=%d ", tt.z);
            if(tt.w>=0) printf("H=%d ", tt.w);
            printf("}\n");
        }
    }
    //printf("\n");
}

void printDOFs() const {
    printf("printDOFvalues() nDOFs=%i\n", nDOFs);
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = DOFtoTyp[i];
        char* tname = params->atypes[rt.x].name;
        const char comp = "REQH"[rt.y];
        printf("%3i ->(%3i|%i) %8s.%c : %20.10f dE/d%c: %20.10f \n", i, rt.x,rt.y, tname,  comp, DOFs[i], comp, fDOFs[i]);
    }
    //printf("\n");
}

void printDOFregularization(){
    printf( "printDOFregularization()\n" );
    //printf("DOF  Type      Comp  Current       Position(min,x0,max)           Stiffness(min,k0,max)\n");
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt   = DOFtoTyp[i];
        const Vec2d& lim  = DOFlimits[i];
        const Vec3d& regX = DOFregX[i];
        const Vec3d& regK = DOFregK[i];
        //printf("[%2i] %8s.%c(%8.3f) limits( %10.2e %10.2e ) reg: x(%10.2e,%10.2e,%10.2e ) K(%10.2e,%10.2e,%10.2e )\n",   i, params->atypes[rt.x].name, "REQH"[rt.y],  DOFs[i], lim.x, lim.y,   regX.x, regX.y, regX.z,      regK.x, regK.y, regK.z        );
        printf("DOF: %2i %-8s %c x: %10.3f limits( %10.2e %10.2e ) reg: x( %10.3f %10.3f %10.3f ) K( %10.3f %10.3f %10.3f )\n",   i, params->atypes[rt.x].name, "REQH"[rt.y],  DOFs[i], lim.x, lim.y,   regX.x, regX.y, regX.z,      regK.x, regK.y, regK.z        );
    }
}

}; // class FitREQ

#endif

