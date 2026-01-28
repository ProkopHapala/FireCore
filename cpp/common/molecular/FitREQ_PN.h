#ifndef FitREQ_PN_h
#define FitREQ_PN_h

#include <vector>
#include <math.h>
#include "Vec3.h"
#include "quaternion.h"
#include "Atoms.h"
#include "MMFFparams.h"
#include "Forces.h"
#include "MMFFBuilder.h"
#include "functions.h"
#include "IO_utils.h"

inline double soft_clamp(double y, double y1, double y2, double& dy_new){
    if (y > y1){
        double y12 = y2 - y1;
        double z   = (y - y1) / y12;
        dy_new = 1.0 / ((1.0 + z)*(1.0 + z));
        return y1 + y12 * (1.0 - 1.0 / (1.0 + z));
    } else {
        dy_new = 1.0;
        return y;
    }
}

constexpr double ARGEXP_LIMIT = 700.0;
constexpr double UPEXP_LIMIT  = 1.0142320547350045e304; // exp(700)
constexpr double DOWNEXP_LIMIT = 9.85967654375977e-305; // exp(-700)
inline double safe_exp(double x) {
    if (x >  ARGEXP_LIMIT) return UPEXP_LIMIT;
    if (x < -ARGEXP_LIMIT) return DOWNEXP_LIMIT;
    return exp(x);
}

constexpr double ARGEXP2_LIMIT  = 350.0;
constexpr double UPEXP2_LIMIT   = 1.004769916070824e151; // exp(350)
constexpr double DOWNEXP2_LIMIT = 9.95245916070824e-152; // exp(-350)
inline double safe_exp2(double x) {
    if (x >  ARGEXP2_LIMIT) return UPEXP2_LIMIT;
    if (x < -ARGEXP2_LIMIT) return DOWNEXP2_LIMIT;
    return exp(x);
}

constexpr double ARGEXP3_LIMIT = 233.0;
constexpr double UPEXP3_LIMIT   = 1.531774938317012e101; // exp(233)
constexpr double DOWNEXP3_LIMIT = 6.52879859425664e-102; // exp(-233)
inline double safe_exp3(double x) {
    if (x >  ARGEXP3_LIMIT) return UPEXP3_LIMIT;
    if (x < -ARGEXP3_LIMIT) return DOWNEXP3_LIMIT;
    return exp(x);
}

inline double getSR_PN( double r, double H, double R0, double& dEdH, double& dEdR0 ){
    const double iR0 = 1.0 / R0;
    const double u  = r * iR0; 
    dEdH = safe_exp( -u );
    dEdR0 = H * u * iR0 * dEdH;
    return H * dEdH;
}

inline double getSR2_PN( double r, double H, double R0, double& dEdH, double& dEdR0 ){
    const double iR0 = 1.0 / R0;
    const double u  = r * iR0; 
    const double u2  = u * u; 
    dEdH = safe_exp( -u2 );
    dEdR0 = 2.0 * H * u2 * iR0 * dEdH;
    return H * dEdH;
}

inline double getSR3_PN( double r, double H, double R0, double& dEdH, double& dEdR0 ){
    const double iR0 = 1.0 / R0;
    const double u  = r * iR0;
    const double ep  = safe_exp3( u );
    const double em  = 1.0 / ep;
    const double s = 1.0 / ( ep + em );
    dEdH = 2.0 * s;
    dEdR0 = H * u * iR0 * ( ep - em ) * dEdH * s;
//if (!std::isfinite(dEdR0)||!std::isfinite(dEdH)) { printf("getSR3_PN: r=%f H=%f R0=%f u=%f dEdH=%f dEdR0=%f\n", r, H, R0, u, dEdH, dEdR0); exit(0); }
    return H * dEdH;
}

inline double coulomb_potential_bare(double r){ return 1.0/r; }

// Soft-clamp on the bare Coulomb potential value y=1/r (no charges, no COULOMB_CONST)
inline double dampCoulomb_SoftClamp(double r, double y1, double y2){
    double y = coulomb_potential_bare(r);
    if(y<=y1) return y;
    double y12 = y2 - y1; if(y12<=0) return y1; // guard
    double z   = (y - y1)/y12;
    return y1 + y12 * (1.0 - 1.0/(1.0+z));
}

// Exact Boys-like damping: erf(r)/r, stable at r->0 where it tends to 2/sqrt(pi)
inline double boys_potential_exact(double r){
    if(r<1e-12) return 1.1283791670955126; // 2/sqrt(pi)
    return ::erf(r)/r;
}

// Polynomial approximations for r < r_min (coefficients correspond to r_min=1.5 from Python derivation)
// NOTE: These are hard-coded for r_min ~ 1.5. If you change boys_rmin, recompute coefficients.
inline double boys_poly_c1_hermite_cubic(double r){
    const double c3=0.07607654346400593, c2=-0.3193203709421615, c0=1.12837916709551;
    double r2=r*r; return (c3*r + c2)*r2 + c0;
}

inline double boys_poly_c2_hermite_quintic(double r){
    const double c5=0.0978009599229947, c4=-0.3186499989199512, c3=0.37187006074364537, c2=-0.3761263890318375, c0=1.12837916709551;
    double r2=r*r; return (((c5*r + c4)*r + c3)*r + c2)*r2 + c0;
}

inline double boys_poly_c1_quartic_even(double r){
    const double c4=0.02535884782133531, c2=-0.26226296334415705, c0=1.12837916709551;
    double r2=r*r; return (c4*r2 + c2)*r2 + c0;
}

inline double boys_poly_c2_sextic_even(double r){
    const double c6=0.010677274768021069, c4=-0.022688888634759506, c2=-0.20820925983105038, c0=1.12837916709551;
    double r2=r*r; return ((c6*r2 + c4)*r2 + c2)*r2 + c0;
}

// Piecewise damped Coulomb using Boys: for r<rmin use selected poly/erf, for r>=rmin fall back to 1/r.
// mode: 0=exact erf/r, 1=C1 Hermite Cubic, 2=C2 Hermite Quintic, 3=C1 Quartic Even, 4=C2 Sextic Even
inline double dampCoulomb_Boys(double r, double rmin, int mode){
    if(r<rmin){
        switch(mode){
            case 0: return boys_potential_exact(r);
            case 1: return boys_poly_c1_hermite_cubic(r);
            case 2: return boys_poly_c2_hermite_quintic(r);
            case 3: return boys_poly_c1_quartic_even(r);
            default:return boys_poly_c2_sextic_even(r);
        }
    }else{
        return coulomb_potential_bare(r);
    }
}

// === Damped Coulomb helpers (free functions) ===
// Bare 1/r potential used as the base (without COULOMB_CONST and without charges)
inline void boys_potential_exact_val_deriv(double r, double& y, double& dy_dr){
    if(r<1e-12){ y=1.1283791670955126; dy_dr=0.0; return; }
    const double exp_r2 = exp(-r*r);
    y     = ::erf(r)/r;
    dy_dr = ( 1.1283791670955126 * r * exp_r2 - ::erf(r) ) / (r*r);
}

inline void boys_poly_c1_hermite_cubic_val_deriv(double r, double& y, double& dy_dr){
    const double c3=0.07607654346400593, c2=-0.3193203709421615, c0=1.12837916709551;
    const double r2=r*r;
    y     = (c3*r + c2)*r2 + c0;
    dy_dr = (3.0*c3*r + 2.0*c2)*r;
}

inline void boys_poly_c2_hermite_quintic_val_deriv(double r, double& y, double& dy_dr){
    const double c5=0.0978009599229947, c4=-0.3186499989199512, c3=0.37187006074364537, c2=-0.3761263890318375, c0=1.12837916709551;
    const double r2=r*r;
    y     = (((c5*r + c4)*r + c3)*r + c2)*r2 + c0;
    dy_dr = (((5.0*c5*r + 4.0*c4)*r + 3.0*c3)*r + 2.0*c2)*r;
}

inline void boys_poly_c1_quartic_even_val_deriv(double r, double& y, double& dy_dr){
    const double c4=0.02535884782133531, c2=-0.26226296334415705, c0=1.12837916709551;
    const double r2=r*r;
    y     = (c4*r2 + c2)*r2 + c0;
    dy_dr = (4.0*c4*r2 + 2.0*c2)*r;
}

inline void boys_poly_c2_sextic_even_val_deriv(double r, double& y, double& dy_dr){
    const double c6=0.010677274768021069, c4=-0.022688888634759506, c2=-0.20820925983105038, c0=1.12837916709551;
    const double r2=r*r;
    y     = ((c6*r2 + c4)*r2 + c2)*r2 + c0;
    dy_dr = ((6.0*c6*r2 + 4.0*c4)*r2 + 2.0*c2)*r;
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
    Quat4d* REQs = 0;       // [natoms]  REQ parameters of H-bond corrected atoms

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

};

// ================================
// ====   class    FitREQ       ====
// =================================

/**
 * @class FitREQ_PN
 * @brief Class for fitting non-colvalent interaction parameters (like Lenard-Jonnes, Charge, Hydrogen Bond-correction) of a system of atoms to a set of training examples.
 * 
 * This class contains various arrays and variables for fitting parameters of a molecular system to a set of training examples.
 * It includes arrays for storing the parameters of each atom type, regularization parameters, and stiffness parameters.
 * It also includes arrays for storing the degrees of freedom (DOFs) of the system, as well as temporary arrays for accumulation of derivatives.
 * Additionally, it includes arrays for decomposition of energy components and arrays for storing the poses of the system.
 * 
 * The class provides methods for initializing the types of atoms, setting the system, setting rigid samples, and evaluating the derivatives of the system.
 */
class FitREQ_PN{ public:
    
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
    alignas(32) double*    typeKregCount =0; // [ntype] number of samples for each type
    alignas(32) Quat4i*    typToREQ      =0; // [ntype] map each unique atom type to place in DOFs;
    std::vector<Vec2i> DOFtoTyp;    // Maps DOF index to (type_index, component)
    std::vector<Vec3d> DOFregX;     // regularization positions (xmin,x0,xmax) for each DOF
    std::vector<Vec3d> DOFregK;     // regularization stiffness (Kmin,K0,Kmax) for each DOF
    std::vector<Vec2d> DOFlimits;   // limits (xmin,xmax) for each DOF
    std::vector<double> DOFinvMass; // inverse mass for each DOF ( for dynamical relaxation )
    std::vector<double> DOFcount;   // number of samples for each DOF
    alignas(32) double*   DOFs =0;          // [nDOFs]
    alignas(32) double*   fDOFs=0;          // [nDOFs]
    alignas(32) double*   vDOFs=0;          // [nDOFs]
    alignas(32) double*   DOFs_old =0;      // [nDOFs]
    alignas(32) double*   fDOFs_old=0;      // [nDOFs]
    alignas(32) Vec2d*    fDOFbounds=0;     // [nDOFs] minimum and maximum value of fDOFs over all samples
    alignas(32) double*   sample_fdofs = 0; // [nDOFs*nsamples] - arrays for parallelization to avoid atomic-write conflicts

    bool bEvalJ                = false; // Should we evaluate variational derivatives on Fregment J 
    bool bWriteJ               = false; // Should we write variational derivatives to Fregment J ( inner loop over j )
    bool bCheckRepulsion       = false; // Should we check maximum repulsion (EijMax) inside inner loop over j for each sample atoms ?
    bool bRegularize           = true;  // Should we apply additional regularization forces to otimizer ( beside the true variational forces from inter-atomic forcefield ? )
    bool bAddRegError          = true;  // Should we add regularization error to total error ?
    bool bEpairs               = true;  // Should we add electron pairs to the molecule ?
    bool bEpairDistByType      = false; // Should we use different electron pair distances for different atom types ?
    bool bBroadcastFDOFs       = false; // Should we broadcast fDOFs (each sample to its own chunk of memory) to prevent atomic-write conflicts ?
    bool bUpdateDOFbounds      = true;  // Should we update fDOFbounds after each sample ?
    bool bClearDOFboundsEpoch  = false; // Should we clear fDOFbounds after each epoch ?
    bool bEvalOnlyCorrections  = false; // Split evaluation and optimization to Emodel0 and Ecorrection (where only Ecorrection is updated every iteration)
    bool bSaveJustElementXYZ   = false; // Should we save just element names in the output .xyz file ?
    bool bRegCountWeight       = false; // Should we weight regularization by number of samples for each type ?  
    bool bClamp                = false; // Should we hardly restrain the parameter values?
    bool bListOverRepulsive    = true;  // Should we list overrepulsive samples? 
    bool bSaveOverRepulsive    = false; // Should we save overrepulsive samples to .xyz file?
    bool bPrintOverRepulsive   = true;  // Should we print overrepulsive samples? 
    bool bDiscardOverRepulsive = true;  // Should we discard overrepulsive samples? ( i.e. ignore them as training examples )
    bool bWeightByEmodel       = true;  // Should we weight samples by their model energy ?
    std::vector<int> overRepulsiveList; // List of overrepulsive samples indexes (for recent epoch)
    char* fname_overRepulsive  = "overRepulsive.xyz";
    bool bPrintDOFs            = false;
    bool bPrintfDOFs           = true;
    bool bPrintBeforReg        = true;
    bool bPrintAfterReg        = false;
    bool bUpdateHostCharge     = true;
    bool bSaveSampleToXYZ      = false;
    char* xyz_out              = "out.xyz";
    std::vector<int> typesPresent;
    std::vector<int> fittedTypes;

    // parameters
    int    iWeightModel    = 1;     // weight of model energy 1=linear, 2=cubic_smooth_step  
    double EmodelCut       = 10.0;  // sample model energy when we consider it too repulsive and ignore it during fitting
    double EmodelCutStart  = 5.0;   // sample model energy when we start to decrease weight in the fitting  
    bool   bSoftClamp      = false; // Should we apply soft clamp to the model energy ?
    double softClamp_start = 4.0;   // if (dE=(Ei_model-Ei_ref) > softClamp_start then apply soft clamp 
    double softClamp_max   = 6.0;   // dE will never be larger than softClamp_max
    double invWsum         = 1.0;
    double kMorse          = 1.8;

    // check pairwise repulsion betwee atoms within one sample
    double EijMax    = 5.0;
    int nsamp_debug  = 10;
    int isamp_debug  = 0;
    int iBadFound    = 0; 
    int nBadFoundMax = 10;
    double*  weights = 0;           // [nbatch] scaling importaince of parameters
    std::vector<Atoms*> samples;    // ToDo: would be more convenient to store Atoms* rather than Atoms
    MM::Builder builder;
    MMFFparams* params = 0; 

    // export optimization trajectory 
    double* trj_E     = 0;
    double* trj_F     = 0;
    double* trj_DOFs  = 0;
    double* trj_fDOFs = 0;

    // === Electrostatics damping/clamping parameters ===
    // Boys-based damping
    double boys_rmin   = 1.5;   // matching radius for piecewise Boys approximations
    int    boys_mode   = 4;     // 0=exact erf/r, 1=cubic C1, 2=quintic C2, 3=quartic even C1, 4=sextic even C2
    // Soft-clamp thresholds on bare Coulomb potential y=1/r (before multiplying by COULOMB_CONST*Qi*Qj)
    double clamp_y1    = 3.0;   // start clamping when y>y1
    double clamp_y2    = 6.0;   // asymptote towards y2 as y->infty
    int    clamp_mode  = 1;     // 1=soft, 2=smooth, 3=soft_neg, 4=smooth_neg
    
    // PN
    int    ivdW      = 4;     // 0=no vdW, 1=LJ, 2=LJr8, 3=LJr9, 4=Morse, 5=Buck
    int    iCoul     = 1;     // 0=no Coul, 1=point charges, 2=soft clamping, 10-14=Boys clamping, 10=exact erf/r, 11=cubic C1, 12=quintic C2, 13=quartic even C1, 14=sextic even C2
    int    iHbond    = 0;     // 0=no HBond correction, 1=H1 correction, 2=H2 correction
    int    iEpairs   = 0;     // 0=no interaction, 1=SR interaction, 2=SR2 interaction, 3=SR3 interaction
    double Lepairs   = 1.0;   // Ang, distance host atom-Epair
    bool   bPN       = true;  // use PN model for vdW and Coulomb interactions

// =================================
// =========== Functions ===========
// =================================

// loadDOFSelection
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
    //              at_name,comp, xmin,xmax, xlo,xhi, klo,khi, K0,       xstart     invMass     
    int nw_min    = 1 + 1 + 2 + 2 + 2 +1;
    int nw_xstart = nw_min + 1;
    DOFtoTyp  .clear();
    DOFregX   .clear();
    DOFregK   .clear();
    DOFlimits .clear();
    DOFinvMass.clear();
    int iDOF = 0;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        int comp;
        double xmin, xmax, xlo, xhi, klo, khi, K0, xstart, invMass;
        int nw = sscanf( line, "%s %i %lf %lf %lf %lf %lf %lf %lf %lf %lf",  at_name, &comp, &xmin, &xmax, &xlo, &xhi, &klo, &khi, &K0, &xstart, &invMass );
        if(nw < nw_min){ printf("ERROR in loadDOFSelection(): expected min. %i words, got %i in line '%s'\n", nw_min, nw, line); exit(0); }
        int ityp = params->getAtomType(at_name);
        if(ityp < 0){  printf("ERROR in loadDOFSelection(): unknown atom type %s\n", at_name); exit(0); }
        // Add new DOF
        fittedTypes[ityp] = 1;
        if(nw<nw_xstart){ 
            printf("iDOF: %3i %8s.%i xstart missing (nw=%3i) => xstart=%16.8f \n", iDOF, at_name, comp, nw, typeREQs0[ityp].array[comp]  );
            xstart = typeREQs0[ityp].array[comp]; 
        }
        if (nw<nw_xstart+1){ printf("iDOF: %i %8s.%i invMass missing (nw=%i) => invMass=%g \n", iDOF, at_name, comp, nw, 1.0 ); invMass = 1.0; }
        DOFtoTyp  .push_back( {ityp, comp}       );
        DOFregX   .push_back( {xlo, xstart, xhi} );
        DOFregK   .push_back( {klo, K0,     khi} );
        DOFlimits .push_back( {xmin, xmax}       );
        DOFinvMass.push_back( invMass           );
        DOFcount  .push_back( 0 );
        //if(verbosity>0)
            printf( "DOF[%3i] %3i|%i %-8s %c | range: %+10.4e ,%+10.4e | reg{ x0: %+10.4e , %+10.4e | K: %+10.4e , %+10.4e } xstart: %+10.4e invMass: %+10.4e \n",  iDOF, ityp,comp,  at_name, "REQH"[comp], xmin, xmax, xlo, xhi, klo, khi, xstart, invMass );
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

int initAllTypes(){
    ntype = params->atypes.size();
    printf( "FitREQ_PN::initAllTypes() ntype=%3i \n", ntype );
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
    _realloc0( typeKregCount,  ntype_, 0.0 );
}

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

void initTypeParamsFromDOFs() {
    //printf("FitREQ_PN::initTypeParamsFromDOFs() typesPresent.size()=%i @DOFcount=%p\n", typesPresent.size(), DOFcount);
    for (int iDOF=0; iDOF<nDOFs; iDOF++) {
        Vec2i rt = DOFtoTyp[iDOF];
        int ityp = rt.x;
        int comp = rt.y;
        typToREQ[ityp].array[comp] = iDOF;
        double xstart = DOFregX[iDOF].y;
        //printf("FitREQ_PN::initTypeParamsFromDOFs() ityp: %3i comp: %i iDOF: %3i xstart: %16.8f\n", ityp, comp, iDOF, xstart);
        typeREQs      [ityp].array[comp] = xstart;  // Initialize with xstart
        typeREQs0     [ityp].array[comp] = xstart;
        typeREQsMin   [ityp].array[comp] = DOFlimits[iDOF].x;
        typeREQsMax   [ityp].array[comp] = DOFlimits[iDOF].y;
        typeREQs0_low [ityp].array[comp] = DOFregX  [iDOF].x;
        typeREQs0_high[ityp].array[comp] = DOFregX  [iDOF].z;
        typeKreg_low  [ityp].array[comp] = DOFregK  [iDOF].x;
        typeKreg_high [ityp].array[comp] = DOFregK  [iDOF].z;
        //DOFcount[iDOF] = typesPresent[ityp];
    }
    //printf("FitREQ_PN::initTypeParamsFromDOFs() DONE\n");
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

void printDOFregularization(){
    printf( "printDOFregularization()\n" );
    //printf("DOF  Type      Comp  Current       Position(min,x0,max)           Stiffness(min,k0,max)\n");
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt   = DOFtoTyp[i];
        const Vec2d& lim  = DOFlimits[i];
        const Vec3d& regX = DOFregX[i];
        const Vec3d& regK = DOFregK[i];
        //printf("[%2i] %8s.%c(%8.3f) limits( %10.2e %10.2e ) reg: x(%10.2e,%10.2e,%10.2e ) K(%10.2e,%10.2e,%10.2e )\n",   i, params->atypes[rt.x].name, "REQH"[rt.y],  DOFs[i], lim.x, lim.y,   regX.x, regX.y, regX.z,      regK.x, regK.y, regK.z        );
        printf("DOF: %2i %-8s %c x: %10.3f limits( %10.2e %10.2e ) reg: x( %10.3f %10.3f %10.3f ) K( %10.3f %10.3f %10.3f ) invMass: %10.3f\n",   i, params->atypes[rt.x].name, "REQH"[rt.y],  DOFs[i], lim.x, lim.y,   regX.x, regX.y, regX.z,      regK.x, regK.y, regK.z , DOFinvMass[i] );
    }
}

// loadXYZ
/**
 * @brief Load an XYZ file into `Atoms` samples; optionally add epairs and export.
 *
 * Behavior:
 *  - When `bAddEpairs=true`, calls `addAndReorderEpairs()` which adds epairs and sets their positions
 *    to `Lepairs` along stored directions. This ensures the exported `.xyz` reflects the geometry used
 *    during fitting (same as in `fillTempArrays()`).
 *  - When `bOutXYZ=true`, also writes `<fname>_Epairs.xyz` with the augmented system.
 *
 * @param fname Path to input .xyz
 * @param bAddEpairs Whether to add epairs
 * @param bOutXYZ Whether to save an augmented .xyz
 * @param bAppend Whether to append to the output file
 * @return Number of batches created
 */
int loadXYZ( const char* fname, bool bAddEpairs=false, bool bOutXYZ=false, char* OutXYZ_fname="", bool bAppend=false ){
    printf( "FitREQ_PN::loadXYZ(%s) bAddEpairs=%i bOutXYZ=%i OutXYZ_fname=%s bAppend=%i\n", fname, bAddEpairs, bOutXYZ, OutXYZ_fname, bAppend );
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    if(!bAppend){
        for(Atoms* a : samples) delete a;
        samples.clear();
        nbatch = 0;
    }
    // --- Open output file
    FILE* fout=0;
    if(bAddEpairs && bOutXYZ){
        if (strlen(OutXYZ_fname)>0){
            sprintf(line,"%s", OutXYZ_fname );
        }else{
            sprintf(line,"%s_Epairs.xyz", fname );
        }
        printf("FitREQ_PN::loadXYZ() opened output file %s\n", line);  
        fout = fopen(line, bAppend ? "a" : "w");
    }
    int il   = 0;
    Atoms* atoms=0;
    while( fgets(line, nline, fin) ){
        if      ( il==0 ){               // --- Read number of atoms
            int na=-1;
            sscanf( line, "%i", &na );
            if( (na>0)&&(na<10000) ){
                atoms = new Atoms(na);
            }else{ printf( "ERROR in FitREQ_PN::loadXYZ() Suspicious number of atoms (%i) while reading `%s`  => Exit() \n", na, fname ); exit(0); }
        }else if( il==1 ){               // --- Read comment line ( read reference energy )
            sscanf( line, "%*s %*s %i %*s %lf ", &(atoms->n0), &(atoms->Energy) );
            //printf("FitREQ_PN::loadXYZ() nbatch[%i] Energy %lf\n", nbatch, atoms->Energy );             
        }else if( il<atoms->natoms+2 ){  // --- Read atom line (type, position, charge)
            double x,y,z,q;
            int nret = sscanf( line, "%s %lf %lf %lf %lf", at_name, &x, &y, &z, &q );
            //printf( ".xyz[%i] %s %lf %lf %lf %lf\n", il, at_name, x, y, z, q );
            if(nret<5){q=0;}
            int i=il-2;
            atoms->apos[i].set(x,y,z);
            atoms->atypes[i]=params->getAtomType(at_name);
            atoms->charge[i]=q;
        }
        il++;
        if( il >= atoms->natoms+2 ){    // ---- store sample atoms to batch
            if(bAddEpairs){ 
                addAndReorderEpairs(atoms, fout); 
            }else{
                AddedData* data = new AddedData();
                _realloc0( data->host, atoms->natoms, -1 );
                atoms->userData = data;
            }
            int nbad = atoms->checkTypeInRange( params->atypes.size()-1, 1, true );
            if( nbad>0 ){ 
                printf("ERROR in FitREQ_PN::loadXYZ() samples[%5i]->checkTypeInRange() return nbad=%i => exit()\n", nbatch, nbad ); 
                params->printTypesOfAtoms(atoms->natoms, atoms->atypes);
                exit(0); 
            }
            samples.push_back( atoms );
            //if(bEvalOnlyCorrections){ printFittedAdata( samples.size()-1 ); }
            il=0; nbatch++;
        }
    }
    if(fout){
        //printf("FitREQ_PN::loadXYZ() closed output file %s\n", );  
        fclose(fout);
    }
    fclose(fin);
    //init_types();
    countTypesPresent( );
    printTypeParams( true );
    return nbatch;
}

/**
 * @brief Add dummy electron pairs and reorder atoms to a canonical layout.
 *
 * Principle:
 *  - Uses `MM::Builder::buildBondsAndEpairs()` to add dummy atoms (electron pairs, sigma holes).
 *  - Reorders atoms to the layout: Atoms(mol1), Epairs(mol1), Atoms(mol2), Epairs(mol2).
 *  - Stores epair bonds and unit directions into `atoms->userData` (see `AddedData`).
 *  - After reordering, it positions epairs at distance `Lepairs` from their host along `dirs`,
 *    so that the saved .xyz is consistent with the geometry used later in `fillTempArrays()` during fitting.
 *
 * Notes:
 *  - Historical default epair distance (~0.5 Å) originates in `MMFFBuilder::addEpair()` where,
 *    if no parameters are bound, a default `l=-0.5` becomes 0.5 Å; otherwise the builder uses
 *    `params->atypes[epairType].Ruff` as the initial length. We override that here to `Lepairs` for consistency.
 *
 * @param atoms Input molecular system to process
 * @param fout Optional file pointer to write XYZ output (can be null)
 */
void addAndReorderEpairs(Atoms*& atoms, FILE* fout=nullptr) {
    //printf( "addAndReorderEpairs() samples.size()=%i \n", samples.size() );
    // Store original system information
    Atoms* bak = atoms;
    int n0bak  = bak->n0;        // Number of atoms in first molecule
    int natbak = bak->natoms;    // Total number of atoms before adding epairs
    int n1bak  = natbak - n0bak; // Number of atoms in second molecule

    // Add electron pairs to the system
    //atoms = addEpairs(atoms);      // This modifies the builder's state
    builder.params = params;
    atoms = builder.buildBondsAndEpairs(atoms);

    atoms->n0 = bak->n0;
    atoms->Energy = bak->Energy;
    for(int i=0; i<natbak; i++){ atoms->charge[i] = bak->charge[i]; }
    

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

    int nbad = atoms->checkTypeInRange( params->atypes.size()-1, 1, true );
    if( nbad>0 ){ 
        printf("ERROR in FitREQ_PN::addAndReorderEpairs() samples[%5i]->checkTypeInRange() return nbad=%i => exit()\n", nbatch, nbad ); 
        params->printTypesOfAtoms(atoms->natoms, atoms->atypes);
        params->printTypesOfAtoms(bak->natoms,   bak->atypes);
        exit(0); 
    }

    delete bak;

    // // Store electron pair relationships in AddedData
    AddedData* data = new AddedData();
    atoms->userData = data;
    data->nep       = nep_found;
    data->bs        = bs;       //  bs[i].x is epair index, bs[i].y is host atom index
    data->dirs      = dirs;   
    _realloc0( data->host, atoms->natoms, -1 );
    for(int i=0; i<nep_found; i++){ data->host[bs[i].y] = bs[i].x; } // bs[i].x is host atom index, bs[i].y is epair index
    //data->isep      = isep;

    // Position epairs at the fitting distance (Lepairs) so that exported .xyz matches fitting geometry
    // This mirrors the placement done later in fillTempArrays(): apos[iE] = apos[iX] + dirs[i]*Lepairs
    if(bEpairs && nep_found>0){
        for(int i=0; i<nep_found; i++){
            const int iX = bs[i].x;  // host atom
            const int iE = bs[i].y;  // epair atom (after reordering)
            atoms->apos[iE] = atoms->apos[iX] + dirs[i]*Lepairs;
        }
    }

    if(bEvalOnlyCorrections){
        int naFit = initFittedAdata( atoms, 0 );  // first pass just count the number of fitted atoms
        data->reallocHB( atoms->natoms, naFit );  // then we can allocate the arrays
        initFittedAdata( atoms, data );           // second pass store the fitted atoms
    }
    
    // Write output if requested
    if(fout) {
        char line[1024];
        sprintf(line, "#	n0 %i E_tot %g", atoms->n0, atoms->Energy);
        params->writeXYZ(fout, atoms, line, 0, bSaveJustElementXYZ);
    }
    //printf( "addAndReorderEpairs() DONE \n" );
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
    updateTypeStats();
}

void updateTypeStats(){
    printf("FitREQ_PN::updateTypeStats(): nDOFs=%i DOFcount.size()=%i typesPresent.size()=%i\n", nDOFs, DOFcount.size(), typesPresent.size() );
    for (int iDOF=0; iDOF<nDOFs; iDOF++) {
        Vec2i rt = DOFtoTyp[iDOF];
        int ityp = rt.x; // int comp = rt.y;
        printf("DOF %3i ityp: %3i typesPresent[ityp]=%i\n", iDOF, ityp, typesPresent[ityp] );
        DOFcount[iDOF] = typesPresent[ityp];
    }
    printf("FitREQ_PN::updateTypeStats(): DONE\n");
}

void printTypeParams( bool bOnlyPresent=true ){
    printf("FitREQ_PN::printTypeParams():\n");
    int ntype = params->atypes.size();
    bool bCounted = typesPresent.size() > 0;
    for(int i=0; i<ntype; i++){
        int ncount = -1;
        if(bCounted){
            ncount = typesPresent[i];
            if( bOnlyPresent && (ncount==0) ){ continue; }
        }
        Quat4d tREQH = typeREQs[i]; 
        //typeKregCount[i] = ncount;
        printf("type %3i %-8s count: %6i REQH: %10.3f %10.3f %10.3f %10.3f \n", i, params->atypes[i].name, ncount, tREQH.x, tREQH.y, tREQH.z, tREQH.w );
    }
}

// run_PN()
//  |-- evalFitError()
//  |    |-- evalSamples_[serial,omp]()
//  |    |    |-- evalSampleError()
//  |    |    |    |-- evalSample_PN()
//  |    |    |    |    | -- evalEnergyDerivs()
//  |    |    |    |-- apply weights and softclamp, then accumulate atom-based-gradients (fREQs) into type-based-gradients (fDOFs)
//  |    |    |-- loop over samples
//  |    |-- apply regularization to fDOFs
//  |-- move_[GD,MD,GD_BB_short,GD_BB_long]()
//  |-- limitDOFs()
__attribute__((hot)) 
double run_PN( int ialg, int nstep, double Fmax, double dt, double max_step, double damping, bool bOMP ){
    bSaveSampleToXYZ=false; 
    //bSaveSampleToXYZ=true; 
    if(verbosity>1){ printf( "FitREQ_PN::run() imodel %i ialg %i nstep %i Fmax %g dt %g max_step %g \n", imodel, ialg, nstep, Fmax, dt, max_step ); }
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
        saveTrajectory(itr, Err, F2);
        if( F2<F2max ){ printf("CONVERGED in %i iterations (|F|=%g < F2max= %g) \n", itr, sqrt(F2), Fmax ); break; }
    }
    printf("VERY FINAL |E|=%.15g  DOFs= ",Err     ); for(int j=0;j<nDOFs;j++){ printf("%.15g ", DOFs[j]); };printf("\n");
    printf("VERY FINAL |F|=%.15g fDOFs= ",sqrt(F2)); for(int j=0;j<nDOFs;j++){ printf("%.15g ",fDOFs[j]); };printf("\n");
    printDOFvalues();
    return Err;
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

void clear_fDOFbounds(){
    for(int i=0; i<nDOFs; i++){ fDOFbounds[i] = Vec2d{1.e+300,-1.e+300}; }
}

void realloc_sample_fdofs(){
    _realloc0( sample_fdofs, nDOFs*samples.size(), 0.0 );
}

// evalFitError
__attribute__((hot)) 
double evalFitError(int itr, bool bOMP=true, bool bEvalSamples=true){
    double Err=0.0;
    DOFsToTypes(); 
    clean_fDOFs();
    if(bEvalSamples)[[likely]]{
        if(bOMP){  Err = evalSamples_omp();   }
        else    {  Err = evalSamples_serial(); }
    }
    if(bPrintBeforReg)printStepDOFinfo( itr, Err, "evalFitError() BEFOR_REG: " );
    if(bRegularize){ 
        double Ereg = regularizeDOFs(); 
        if(bAddRegError)Err += Ereg;
    }   
    if(bPrintAfterReg)printStepDOFinfo( itr, Err, "evalFitError() AFTER_REG: " );
    return Err;
}

void DOFsToTypes(){
    //printf( "DOFsToTypes() nDOFs: %i @DOFtoTyp=%p \n", nDOFs, DOFtoTyp );
    for(int i=0; i<nDOFs; i++){
        const Vec2i& rt = DOFtoTyp[i];
        //printf( "DOFsToType()[%3i] rt(%3i|%i) %8s.%c  %g \n", i, rt.x,rt.y,   params->atypes[rt.x].name,"REQH"[rt.y], DOFs[i] );
        typeREQs[rt.x].array[rt.y] = DOFs[i];
    }
}

/// @brief Cleans the derivative array by setting all its elements to zero.
void clean_fDOFs    ()                { for(int i=0; i<nDOFs; i++){ fDOFs[i]=0;                          }}

/**
 * @brief Evaluates the variational derivatives of the fitting error (sumed batch training samples) with respect to all fitting parameters (i.e. non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) for each atom type). The atomic systems are not assumed rigid (i.e. the atomic positions can vary freely from sample to sample).
 * 
 * @param Eout array to store the non-covalent interaction energy values of each atomic system in the batch. if Eout==null, the function will not store the energy values.
 * @return double, returns the total fitting error.
 */
__attribute__((hot)) 
double evalSamples_serial( double* Eout=0 ){ 
    if( bListOverRepulsive   ) overRepulsiveList.clear();
    if( bSaveOverRepulsive   ) clearFile( fname_overRepulsive );
    if( bClearDOFboundsEpoch ) clear_fDOFbounds(); // clean fDOFbounds
    bBroadcastFDOFs    = false; 
    double Error = 0.0;
    int nsamp = samples.size();
    //printf( "evalSamples_serial() nsamp=%i\n", nsamp );
    for(int isamp=0; isamp<nsamp; isamp++){
        double E; 
        Error+=evalSampleError( isamp, E );
//printf( "evalSamples_serial() isamp=%i E=%18.8e\n", isamp, E );
        if(Eout)Eout[isamp]=E;
    }
    if(bPrintOverRepulsive){ printOverRepulsiveList(); }
//exit(0); // JAMME
    return Error;
}

__attribute__((hot)) 
double evalSamples_omp( double* Eout=0 ){ 
    bListOverRepulsive = false;
    bBroadcastFDOFs    = true; 
    double Error       = 0.0;
    int nsamp = samples.size();
    //printf( "evalSamples_omp() nsamp=%i\n", nsamp );
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
double evalSampleError( int isamp, double& E ){
    //printf( "evalSampleError() isamp = %i \n", isamp );
    isamp_debug = isamp;
    Atoms* atoms  = samples[isamp];
    double wi     = (weights)? weights[isamp] : 1.0;
    if(wi<-1e-300) return 0;
    alignas(32) Quat4d fREQs [atoms->natoms];
    alignas(32) double fDOFs_[nDOFs];
    //if(bPN){ E = evalSample_PN( isamp, atoms, wi, fREQs ); }
    //else{ E = evalSample( isamp, atoms, wi, fREQs ); }
    E = evalSample_PN( isamp, atoms, wi, fREQs );

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
        sprintf(comment, "# %4i Eref: %+10.4e Emodel: %+10.4e wi: %+10.4e", isamp, atoms->Energy, E, wi );
        //printf( "evalSampleError() saving %s comment: %s \n", xyz_out, comment );
        saveDebugXYZ( 0, atoms->natoms, atoms->atypes, atoms->apos, xyz_out, comment );
    }
    double Eref       = atoms->Energy;
    double dE_        = E - Eref;
    double dClamp_dE  = 1.0;
    double dE = (bSoftClamp) ? soft_clamp(dE_, softClamp_start, softClamp_max, dClamp_dE ) : dE_;
    ///printf( "evalSampleError() isamp: %3i Emodel: %20.6f Eref: %20.6f bBroadcastFDOFs=%i @sample_fdofs=%p \n", isamp, E, atoms->Energy, bBroadcastFDOFs, sample_fdofs );
    wi*= invWsum;
    double dEw        = 2.0*dE*wi*dClamp_dE; // chain derivatives   d(wi*Delta_E^2)/dE  = 2*wi*Delta_E * (Delta_E/dE)
    double Error      =  dE*dE*wi;
    //if(isamp_debug<nsamp_debug){ printf( "evalSampleError() isamp: %3i Emodel: %20.6f Eref: %20.6f bBroadcastFDOFs=%i @sample_fdofs=%p \n", isamp, E, atoms->Energy, bBroadcastFDOFs, sample_fdofs ); }
    //if(isamp_debug<nsamp_debug){ 
    //    if( abs(dE_)>1 )
    //    printf( "evalSampleError() isamp: %3i  dEw: %+10.4e wi: %+10.4e dE: %+10.4e dE_: %+10.4e Emodel: %+10.4e Eref: %+10.4e \n", isamp, dEw, wi, dE, dE_, E, atoms->Energy ); 
    //}
    double* fDOFs__ = bBroadcastFDOFs ? sample_fdofs + isamp*nDOFs : fDOFs_;   // broadcast fDOFs ?
    for(int k=0; k<nDOFs; k++){ fDOFs__[k]=0; }                                // clean fDOFs
    acumDerivs( atoms->natoms, atoms->atypes, dEw, fREQs, fDOFs_ );            // accumulate fDOFs from fREQs
    if( bEpairs && bUpdateHostCharge ){
        AddedData * ad = (AddedData*)atoms->userData;
        acumHostDerivs( ad->nep, ad->bs, atoms->atypes, dEw, fREQs, fDOFs_ );
    }
    if( bUpdateDOFbounds  ){ updateDOFbounds( fDOFs__ ); }
    if( !bBroadcastFDOFs ){ 
        double F2=0.0;
        for(int k=0; k<nDOFs; k++){ double fi=fDOFs__[k]; fDOFs[k] += fi; F2+=fi*fi;  } 
        //printf( "evalSampleError() isamp: %3i Eref: %20.6f Emodel: %20.6f |F|: %20.6f  wi: %20.6f dEw: %20.6f \n", isamp, atoms->Energy, E, sqrt(F2), wi, dEw );
    }
    //printf("isamp=%i E=%g Eref=%g dE_=%g dE=%g wi=%g Error=%g\n", isamp, E, Eref, dE_, dE, wi, Error);
    return Error;
}

__attribute__((hot)) 
double evalSample_PN( int isamp, const Atoms* atoms, double wi, Quat4d* fREQs ) const {
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
    double E=0.0;
    E = evalEnergyDerivs    ( i0, ni, j0, nj, atoms->atypes, apos, typeREQs, Qs, adata->host, fREQs );
    return E;
}

__attribute__((hot)) 
void fillTempArrays( const Atoms* atoms, Vec3d* apos, double* Qs  )const{
    //printf( "FillTempArrays() bEpairs=%i \n", bEpairs );
    for(int j=0; j<atoms->natoms; j++){
        Qs[j]   = atoms->charge[j];
        apos[j] = atoms->apos [j];
        
        //{ // for types with fitted charges update charge from DOFs
            int ityp = atoms->atypes[j];
            Quat4i tt = typToREQ[ityp];
            if( tt.z>=0 ){ Qs[j] = DOFs[tt.z]; };
        //}
        //isEp[j] = 0;
        //printf( "FillTempArrays()[ia=%3i] t %3i %-8s iDOF=%3i Q=%12.3e \n", j, ityp, params->atypes[ityp].name,  tt.z, Qs[j] );
    }
    if( bEpairs ){  // electron pairs
        AddedData * ad = (AddedData*)atoms->userData;
        int   nep = ad->nep;
        Vec2i* bs = ad->bs;
        //printf( "FillTempArrays() 2  ad->nep=%i \n", ad->nep );
        for(int j=0; j<ad->nep; j++){
            int    iX  = bs[j].x; // index of the host atom
            int    iE  = bs[j].y; // index of the electron pair
            //Quat4d REQ = typeREQs[atoms->atypes[iE]];
            //printf( "FillTempArrays()[iap=%i] iE=%i iX=%i  name=%8s REQ(%12.3e,%12.3e,%12.3e,%12.3e) \n", j, iE, iX, params->atypes[atoms->atypes[iE]].name, REQ.x,REQ.y,REQ.z,REQ.w );
            //double Qep = REQ.z; // charge of the electron pair
            
            int ityp = atoms->atypes[iE];
            Quat4i tt = typToREQ[ityp];
            double Qep = DOFs[tt.z];
            //printf( "FillTempArrays()[iap=%3i] iE=%3i iX=%3i  t %3i %-8s iDOF=%3i Qep=%12.3e \n", j, iE, iX, ityp, params->atypes[ityp].name,  tt.z, Qep );

            Qs[iE]     = Qep;
            if(bUpdateHostCharge){
                Qs[iX]    -= Qep;
            }
            // #if IF_DEBUG
            //     if(  (fabs(Qs[iX])>1e+10) || (fabs(Qs[iE])>1e+10) ){ printf( "fillTempArrays() j=%i Qs[iX]=%12.3e Qs[iE]=%12.3e Qep=%12.3e \n", j, Qs[iX], Qs[iE], Qep ); exit(0); }
            // #endif
            double lep = Lepairs;
            //if(isamp_debug<1){ printf( "FillTempArrays() isamp %3i [iap=%3i] iE=%3i iX=%3i  t %3i %-8s iDOF=%3i Qep=%12.3e lep=%12.3e \n", isamp_debug, j, iE, iX, ityp, params->atypes[ityp].name,  tt.z, Qep, lep ); }
            if( bEpairDistByType ){ typeREQs[atoms->atypes[iE]].w; }
            apos[iE] = apos[iX] + ad->dirs[j] * lep;  // We move the electron pair to proper distance from the atom
            //isEp[iE] = 1;
        }
    }
}

__attribute__((hot)) 
double evalEnergyDerivs ( int i0, int ni, int j0, int nj, int* types, Vec3d* ps, Quat4d* typeREQs, double* Qs, int* host, Quat4d* dEdREQs )const{
    double E_tot = 0.0;
//double E_Coul = 0.0, E_vdW = 0.0, E_Hcorr = 0.0, E_Epairs = 0.0;    
    //const double alpha   = kMorse;
//printf("JAMME     i     -dE_dR0_i       -dE_dRi         -dE_deps_i      -dE_dQ_i        -dE_dH_i        j     -dE_dR0_j       -dE_dRj         -dE_deps_j      -dE_dQ_j        -dE_dH_j        E_vdW           E_Hcorr         E_Coul          E_Epair\n");
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system2
        const int     i     = i0+ii;
        const int     ih    = host[i];
        const bool    bEpi  = ih>=0;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;

        
        for(int jj=0; jj<nj; jj++){ // loop over all atoms[j] in system1
            const int     j    = j0+jj;
            const int     jh   = host[j];
            const bool    bEpj = jh>=0;
            const double  Qj   = Qs[j];
            const Vec3d   dij  = ps[j] - pi;
            const int     tj   = types[j];
            const Quat4d& REQj = typeREQs[tj];
            const double  R0   = REQi.x + REQj.x;
            const double  eps  = REQi.y * REQj.y; 
            const double  Q    = Qi     * Qj    ;
            double        H    = REQi.w * REQj.w;
            const double  sH   = (H<0.0) ? 1.0 : 0.0; 
            const double  r    = dij.norm();

            double Eij_Coul = 0.0, Eij_vdW = 0.0, Eij_Hcorr = 0.0, Eij_Epairs = 0.0;
            double dE_dR0 = 0.0, dE_deps = 0.0, dE_dQ = 0.0, dE_dH = 0.0;
            double fR = 0.0, fA = 0.0, fH1 = 0.0, fH2 = 0.0, fR0 = 0.0, fB = 1.0;
            double alpha = kMorse;
//double dE_dRi = 0.0, dE_dRj = 0.0; // JAMME

            if( bEpi ){  
                if(bEpj) continue; // dummy atoms should not interact with each other
/*
if(bEpj){
    double dE_dR = 0.0;
    if(iEpairs==1){ 
        Eij_Epairs = getSR_PN( r, H, REQj.x, dE_dH, dE_dR );
    }else if(iEpairs==2){ 
        Eij_Epairs = getSR2_PN( r, H, REQj.x, dE_dH, dE_dR );
    }
    dEdREQs[j].x -= dE_dR;
}
*/
                // --- Electron pair interaction
                double dE_dR = 0.0;
                if(iEpairs==1){ 
                    Eij_Epairs = getSR_PN( r, H, REQi.x, dE_dH, dE_dR );
                }else if(iEpairs==2){ 
                    Eij_Epairs = getSR2_PN( r, H, REQi.x, dE_dH, dE_dR );
                }else if(iEpairs==3){ 
                    Eij_Epairs = getSR3_PN( r, H, REQi.x, dE_dH, dE_dR );
                }
                dEdREQs[i].x -= dE_dR;
//dE_dRi = dE_dR; // JAMME
            }else if( bEpj ){
                // --- Electron pair interaction
                double dE_dR = 0.0;
                if(iEpairs==1){ 
                    Eij_Epairs = getSR_PN( r, H, REQj.x, dE_dH, dE_dR );
                }else if(iEpairs==2){ 
                    Eij_Epairs = getSR2_PN( r, H, REQj.x, dE_dH, dE_dR );
                }else if(iEpairs==3){ 
                    Eij_Epairs = getSR3_PN( r, H, REQj.x, dE_dH, dE_dR );
                }
                dEdREQs[j].x -= dE_dR;
//dE_dRj = dE_dR; // JAMME
            }else{
                // --- Electrostatic interaction
                if(iCoul==1){ // point charges
                    dE_dQ = COULOMB_CONST / r ;
                }else if(iCoul==2){ // point charges with softclamp
                    dE_dQ = dampCoulomb_SoftClamp(r, clamp_y1, clamp_y2) * COULOMB_CONST;
                }else if(iCoul>9){ // Boys clamping with different approximations
                    dE_dQ = dampCoulomb_Boys(r, boys_rmin, iCoul-10) * COULOMB_CONST;
                }
                Eij_Coul = Q * dE_dQ;                
                // --- Van der Waals interaction
                if(ivdW==1){ // Lennard-Jones 12-6
                    const double u  = R0 / r;
                    const double u3 = u * u * u;
                    fA              = u3 * u3;
                    fR              = fA * fA;
                    fH1             = 1.0;
                    fH2             = 2.0;
                    fR0             = 12.0 / R0;
                }else if(ivdW==2){ // Lennard-Jones 8-6
                    const double u  = R0 / r;
                    const double u2 = u * u;
                    fA              = u2 * u2 * u2;
                    fR              = fA * u2;
                    fH1             = 3.0;
                    fH2             = 4.0;
                    fR0             = 24.0 / R0;
                }else if(ivdW==3){ // Lennard-Jones 9-6
                    const double u  = R0 / r;
                    const double u3 = u * u * u;
                    fA              = u3 * u3;
                    fR              = fA * u3;
                    fH1             = 2.0;
                    fH2             = 3.0;
                    fR0             = 18.0 / R0;
                }else if(ivdW==4){ // Morse
                    if( kMorse < 0.0 ){ 
                        alpha = 6.0 / R0; 
                        fR0   = 2.0 * alpha * ( 1.0 + (r - R0) / R0);
                    }else{
                        fR0   = 2.0 * alpha;
                    }
                    fA                 = safe_exp2( -alpha * ( r - R0 ) );
                    fR                 = fA * fA;
                    fH1                = 1.0;
                    fH2                = 2.0;
                }else if(ivdW==5){ // Buckingham
                    if( kMorse < 0.0 ){ 
                        alpha = 6.0 / R0; 
                        fB    = 1.0 + (r - R0) / R0;
                    }
                    const double u     = R0 / r;
                    const double u3    = u * u * u;
                    const double e     = safe_exp2( -alpha * ( r - R0 ) );
                    fA                 = u3 * u3;
                    fR                 = e * e;
                    fH1                = 1.0;
                    fH2                = 2.0;
                    fR0                = 2.0 * alpha;
                }
                dE_deps = fH1 * fR - fH2 * fA;
                dE_dR0  = eps * fR0 * ( fB * fR - fA );
                Eij_vdW   = eps * dE_deps;
                // --- Hydrogen-bond corrections
                if(sH>0.0){
                    if(iHbond==1||iHbond==3){
                        const double f         = fH1 * fR;
                        const double dE_deps_H = H * f;
                        dE_deps       += dE_deps_H;
                        dE_dR0        += H * eps * fB * fR0 * fR;
                        dE_dH          = eps * f;
                        Eij_Hcorr      = eps * dE_deps_H;
                    }
                    if(iHbond==2||iHbond==3){
                        const double f         = fH2 * fA;
                        const double dE_deps_H = H * f;
                        dE_deps       += dE_deps_H;
                        dE_dR0        += H * eps * fR0 * fA;
                        dE_dH         += eps * f;
                        Eij_Hcorr     += eps * dE_deps_H;
                    }
                }
            }

            double Eij = Eij_Coul + Eij_vdW + Eij_Hcorr + Eij_Epairs;
//E_Coul += Eij_Coul; E_vdW += Eij_vdW; E_Hcorr += Eij_Hcorr; E_Epairs += Eij_Epairs;            
/*
printf("JAMME %5i %15.9f %15.9f %15.9f %15.9f %15.9f %5i %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n",
j+1, -dE_dR0,-dE_dRj,-dE_deps*REQi.y,-dE_dQ*Qi,-dE_dH*REQi.w, 
i+1, -dE_dR0,-dE_dRi,-dE_deps*REQj.y,-dE_dQ*Qj,-dE_dH*REQj.w, Eij_vdW, Eij_Hcorr, Eij_Coul, Eij_Epairs);
*/
            dEdREQs[j].sub( Quat4d{ 
                dE_dR0, 
                dE_deps * REQi.y, 
                dE_dQ   * Qi, 
                dE_dH   * REQi.w 
            } );

            // --- Energy and forces
            E_tot    +=  Eij;
            fREQi.x -= dE_dR0;              // dEtot/dR0_i
            fREQi.y -= dE_deps * REQj.y;    // dEtot/dE0_i
            fREQi.z -= dE_dQ   * Qj;        // dEtot/dQ_i
            fREQi.w -= dE_dH   * REQj.w;    // dEtot/dH_i
        }  // end loop over all atoms[j] in system1

        dEdREQs[i].add(fREQi);
    }  // end loop over all atoms[i] in system2
/*
printf("JAMMETOT     Etot\n");
printf("JAMMETOT %15.9f\n", E_tot );
printf("JAMMETOT     i     x_i             y_i             z_i             R0_i            eps_i           q_i             H_i             F_i[1]          F_i[2]          F_i[3]          F_i[4]\n");
for(int jj=0; jj<nj; jj++){ // loop over all atoms[j] in system1
    const int     j     = j0+jj;
    const int     tj    = types   [j ];
    const Quat4d& REQj  = typeREQs[tj];
    printf("JAMMETOT %5i %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
        j+1, ps[j].x,ps[j].y,ps[j].z, REQj.x,REQj.y,Qs[j],REQj.w, dEdREQs[j].x,dEdREQs[j].y,dEdREQs[j].z,dEdREQs[j].w );
}
for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system2
    const int     i     = i0+ii;
    const int     ti    = types   [i ];
    const Quat4d& REQi  = typeREQs[ti];
    printf("JAMMETOT %5i %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f\n", 
        i+1, ps[i].x,ps[i].y,ps[i].z, REQi.x,REQi.y,Qs[i],REQi.w, dEdREQs[i].x,dEdREQs[i].y,dEdREQs[i].z,dEdREQs[i].w );
}
*/
//exit(0);
    return E_tot;
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
        //fprintf(fout, "%c %7.3f %7.3f %7.3f \n", params->atypes[ti].name[0], pi.x, pi.y, pi.z );
        fprintf(fout, "%s %23.15f %23.15f %23.15f \n", params->atypes[ti].name, pi.x, pi.y, pi.z );
    }
    fclose(fout);
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
        //if( ( tt.x>=0 ) || ( tt.y>=0 ) || ( tt.z>=0 ) || ( tt.w>=0 ) )printf( "acumDerivs() ia: %3i  %2i|%-8s  dEw= %g fREQ( %10.2e %10.2e %10.2e %10.2e ) tt(%2i,%2i,%2i,%2i)\n", ia, t, params->atypes[t].name, dEw, f.x,f.y,f.z,f.w, tt.x,tt.y,tt.z,tt.w );
        //if( ( tt.x>=0 ) || ( tt.y>=0 ) || ( tt.z>=0 ) || ( tt.w>=0 ) )printf( "acumDerivs() ia: %3i  %2i|%-8s f.z %10.2e f.z*dEw %10.2e dEw %10.2e \n", ia, t, params->atypes[t].name, f.z, f.z*dEw, dEw );
        //if((tt.y>=0)&&(i==0))printf( "acumDerivs i= %i f= %g f*dEw= %g fDOFs= %g\n", i, f.y, f.y*dEw, fDOFs[tt.y] );
        //if((tt.x>=0)||(tt.y>=0))printf( "acumDerivs i= %i dE= %g f.x= %g fDOFs= %g f.y= %g fDOFs= %g\n", i, dE, f.x, fDOFs[tt.x], f.y, fDOFs[tt.y] );
    }
    // if( isamp_debug < nsamp_debug ){
    //     printf( "acumDerivs() isamp=%i  dEw %+10.6e\n", isamp_debug, dEw );
    //     for(int i=0; i<nDOFs; i++){  printf( "fDOFs[%3i] %+10.6e\n", i, fDOFs[i] );}
    // }
    //exit(0);
}

__attribute__((hot)) 
void acumHostDerivs( int nepair, Vec2i* epairAndHostIndex, int* types, double dEw, Quat4d* fREQs, double* fDOFs  ){
    //printf( "acumHostDerivs() nepair %i\n", nepair );
    for(int i=0; i<nepair; i++){
        Vec2i b             = epairAndHostIndex[i];
        int  tE             = types[b.y];        // type of the electron pair
        const Quat4i&  ttE  = typToREQ[tE];      // get index of degrees of freedom for atom type t
        //printf( "acumHostDerivs() iEpair %i t %i %-8s tt.z %i f.z %g f.z*dEw %g dEw %g\n", i, t, params->atypes[t].name, tt.z, f.z, f.z*dEw, dEw );
        if(ttE.z>=0){ 
            //int tX              = types[b.x];        // type of the host atom
            fDOFs[ttE.z]       -= fREQs[b.x].z*dEw; 
            //printf( "acumHostDerivs() iEpair %3i t %3i %-8s tt.z %2i f.z %10.2e f.z*dEw %10.2e dEw %10.2e\n", i, tE, params->atypes[tE].name, ttE.z, fREQs[b.x].z, fREQs[b.x].z*dEw, dEw );
        } // we subtract the variational derivative of the host atom charge because the charge is transfered from host to the electron pair
    }
}

void updateDOFbounds( double* fDOFs_ ){ for(int i=0; i<nDOFs; i++){ fDOFbounds[i].enclose( fDOFs_[i] );  }}

void printOverRepulsiveList(){
    printf( "printOverRepulsiveList() nfound=%i: { ", overRepulsiveList.size() );
    for(auto i : overRepulsiveList){ printf( "%i,", i ); }  
    printf( "}\n" );  
}

void reduce_sample_fdofs(){
    int nsamples = samples.size();
    for(int j=0; j<nDOFs; j++){ fDOFs[j] = 0.0; }
    for(int i=0; i<nsamples; i++){
        double* fs = sample_fdofs + i*nDOFs;
        for(int j=0; j<nDOFs; j++){ fDOFs[j] += fs[j]; }
    }
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

// New regularization function operating per-DOF
__attribute__((hot)) 
double regularizeDOFs(){
    //printf("regularizeDOFs() nDOFs=%i @DOFregX=%p @DOFregX=%p @DOFs=%p @fDOFs=%p @DOFtoTyp=%p \n",  nDOFs, DOFregX, DOFregK, DOFs, fDOFs, DOFtoTyp); 
    double E = 0;
    for(int i=0; i<nDOFs; i++){   
        double w = bRegCountWeight ? 1.0/DOFcount[i] : 1.0; // weight by number of samples for each type
        E += constrain( DOFs[i], DOFregX[i], DOFregK[i]*w, fDOFs[i] );
        //{ Vec2i tc    = DOFtoTyp[i]; printf( "regularizeDOFs() i: %3i   %8s|%c x: %10.3e f: %10.3e E: %10.3e\n", i,  params->atypes[tc.x].name, "REQH"[tc.y],  DOFs[i], fDOFs[i], E );   }
    }
    //printf( "regularizeDOFs() END E: %10.3e\n", E );
    return E;
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

// move
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
        DOFs[i] += f*dt;
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
        DOFs[i] += f*dt;
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
        //printf( "move_MD fv= %g\n", fv );
        if(fv<0.0){      for(int i=0; i<nDOFs; i++){       vDOFs[i] = 0.0;    }; printf( "ClimbBreak <f|v>=%g\n", fv ); } 
    }
    double F2 = 0.0;
    if(max_step>0.0){ dt=limit_dt_MD(dt,max_step,cdamp); };
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        double v = vDOFs[i];
        v       *= cdamp;
        v       += f*dt * DOFinvMass[i];
        vDOFs[i] = v;
        DOFs[i] += v*dt;
        F2      += f*f;
    }
    return F2;
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
void limitDOFs(){
    for(int i=0; i<nDOFs; i++){   DOFs[i] = _clamp(DOFs[i], DOFlimits[i].x, DOFlimits[i].y); }
    //printf("limitDOFs() DOFs: "); for(int i=0;i<nDOFs;i++){ printf("%8.2e ", DOFs[i]); };printf("\n");
    //printf("limitDOFs() DOFs: "); for(int i=0;i<nDOFs;i++){ printf("%8.2e ", DOFlimits[i].x); };printf("\n");
    //printf("limitDOFs() DOFs: "); for(int i=0;i<nDOFs;i++){  printf("%8.2e ", DOFlimits[i].y); };printf("\n");
}

void saveTrajectory(int itr, double Err, double F2){
    if(trj_E    ){ trj_E[itr]=Err; }
    if(trj_F    ){ trj_F[itr]=F2; }
    int i0=itr*nDOFs;
    if(trj_DOFs ){ for(int i=0;i<nDOFs;i++){ trj_DOFs [i0+i]=DOFs [i]; } }
    if(trj_fDOFs){ for(int i=0;i<nDOFs;i++){ trj_fDOFs[i0+i]=fDOFs[i]; } }
}

void printDOFvalues(){
    printf( "FitREQ_PN::printDOFvalues() \n" );
    printf( "# i          ityp    type_name                x                   dE/dx \n" );
    for(int i=0;i<nDOFs;i++){
        int ityp = DOFtoTyp[i].x;
        printf( "DOF %3i  t: %3i  %-8s %c :  %30.15f   %10.2e \n", i, ityp, params->atypes[ityp].name, "REQH"[DOFtoTyp[i].y],  DOFs[i], fDOFs[i] );
    }
}

void evalEnergyComponents ( int i0, int ni, int j0, int nj, int* types, Vec3d* ps, Quat4d* typeREQs, double* Qs, int* host, double& E_Coul, double& E_vdW, double& E_Hcorr, double& E_Epairs )const{
    E_Coul = 0.0, E_vdW = 0.0, E_Hcorr = 0.0, E_Epairs = 0.0;    
    for(int ii=0; ii<ni; ii++){ // loop over all atoms[i] in system2
        const int     i     = i0+ii;
        const int     ih    = host[i];
        const bool    bEpi  = ih>=0;
        const Vec3d&  pi    = ps      [i ]; 
        const double  Qi    = Qs      [i ]; 
        const int     ti    = types   [i ];
        const Quat4d& REQi  = typeREQs[ti];
        Quat4d        fREQi = Quat4dZero;
        for(int jj=0; jj<nj; jj++){ // loop over all atoms[j] in system1
            const int     j    = j0+jj;
            const int     jh   = host[j];
            const bool    bEpj = jh>=0;
            const double  Qj   = Qs[j];
            const Vec3d   dij  = ps[j] - pi;
            const int     tj   = types[j];
            const Quat4d& REQj = typeREQs[tj];
            const double  R0   = REQi.x + REQj.x;
            const double  eps  = REQi.y * REQj.y; 
            const double  Q    = Qi     * Qj    ;
//const double  Q    = Qi     * Qj    * 1.;
            double        H    = REQi.w * REQj.w;
            const double  sH   = (H<0.0) ? 1.0 : 0.0; 
            const double  r    = dij.norm();

            double Eij_Coul = 0.0, Eij_vdW = 0.0, Eij_Hcorr = 0.0, Eij_Epairs = 0.0;
            double dE_dR0 = 0.0, dE_deps = 0.0, dE_dQ = 0.0, dE_dH = 0.0;
            double fR = 0.0, fA = 0.0, fH1 = 0.0, fH2 = 0.0, fR0 = 0.0, fB = 1.0;
            double alpha = kMorse;

            if( bEpi ){  
                if(bEpj) continue; // dummy atoms should not interact with each other
                // --- Electron pair interaction
                double dE_dR = 0.0;
                if(iEpairs==1){ 
                    Eij_Epairs = getSR_PN( r, H, REQi.x, dE_dH, dE_dR );
                }else if(iEpairs==2){ 
                    Eij_Epairs = getSR2_PN( r, H, REQi.x, dE_dH, dE_dR );
                }
            }else if( bEpj ){
                // --- Electron pair interaction
                double dE_dR = 0.0;
                if(iEpairs==1){ 
                    Eij_Epairs = getSR_PN( r, H, REQj.x, dE_dH, dE_dR );
                }else if(iEpairs==2){ 
                    Eij_Epairs = getSR2_PN( r, H, REQj.x, dE_dH, dE_dR );
                }
            }else{
                // --- Electrostatic interaction
                if(iCoul==1){ // point charges
                    dE_dQ = COULOMB_CONST / r ;
                }else if(iCoul==2){ // point charges with softclamp
                    dE_dQ = dampCoulomb_SoftClamp(r, clamp_y1, clamp_y2) * COULOMB_CONST;
                }else if(iCoul>9){ // Boys clamping with different approximations
                    dE_dQ = dampCoulomb_Boys(r, boys_rmin, iCoul-10) * COULOMB_CONST;
                }
                Eij_Coul = Q * dE_dQ;                
                // --- Van der Waals interaction
                if(ivdW==1){ // Lennard-Jones 12-6
                    const double u  = R0 / r;
                    const double u3 = u * u * u;
                    fA              = u3 * u3;
                    fR              = fA * fA;
                    fH1             = 1.0;
                    fH2             = 2.0;
                    fR0             = 12.0 / R0;
                }else if(ivdW==2){ // Lennard-Jones 8-6
                    const double u  = R0 / r;
                    const double u2 = u * u;
                    fA              = u2 * u2 * u2;
                    fR              = fA * u2;
                    fH1             = 3.0;
                    fH2             = 4.0;
                    fR0             = 24.0 / R0;
                }else if(ivdW==3){ // Lennard-Jones 9-6
                    const double u  = R0 / r;
                    const double u3 = u * u * u;
                    fA              = u3 * u3;
                    fR              = fA * u3;
                    fH1             = 2.0;
                    fH2             = 3.0;
                    fR0             = 18.0 / R0;
                }else if(ivdW==4){ // Morse
                    if( kMorse < 0.0 ){ 
                        alpha = 6.0 / R0; 
                        fR0   = 2.0 * alpha * ( 1.0 + (r - R0) / R0);
                    }else{
                        fR0   = 2.0 * alpha;
                    }
                    fA                 = exp( -alpha * ( r - R0 ) );
                    fR                 = fA * fA;
                    fH1                = 1.0;
                    fH2                = 2.0;
                }else if(ivdW==5){ // Buckingham
                    if( kMorse < 0.0 ){ 
                        alpha = 6.0 / R0; 
                        fB    = 1.0 + (r - R0) / R0;
                    }
                    const double u     = R0 / r;
                    const double u3    = u * u * u;
                    const double e     = exp( -alpha * ( r - R0 ) );
                    fA                 = u3 * u3;
                    fR                 = e * e;
                    fH1                = 1.0;
                    fH2                = 2.0;
                    fR0                = 2.0 * alpha;
                }
                dE_deps = fH1 * fR - fH2 * fA;
                dE_dR0  = eps * fR0 * ( fB * fR - fA );
                Eij_vdW   = eps * dE_deps;
                // --- Hydrogen-bond corrections
                if(sH>0.0){
                    if(iHbond==1||iHbond==3){
                        const double f         = fH1 * fR;
                        const double dE_deps_H = H * f;
                        dE_deps       += dE_deps_H;
                        dE_dR0        += H * eps * fB * fR0 * fR;
                        dE_dH          = eps * f;
                        Eij_Hcorr        = eps * dE_deps_H;
                    }
                    if(iHbond==2||iHbond==3){
                        const double f         = fH2 * fA;
                        const double dE_deps_H = H * f;
                        dE_deps       += dE_deps_H;
                        dE_dR0        += H * eps * fR0 * fA;
                        dE_dH         += eps * f;
                        Eij_Hcorr       += eps * dE_deps_H;
                    }
                }
            }

            E_Coul += Eij_Coul; 
            E_vdW += Eij_vdW; 
            E_Hcorr += Eij_Hcorr; 
            E_Epairs += Eij_Epairs;            

        }  // end loop over all atoms[j] in system1

    }  // end loop over all atoms[i] in system2
}

Vec2d getMinMax( int n, double* vec ){
    Vec2d bd = Vec2d{+1e+300,-1e+300};
    for(int i=0; i<n; i++){
        double v = vec[i];
        if(v<bd.x){ bd.x=v; }
        if(v>bd.y){ bd.y=v; }
    }
    return bd;
} 

}; // class FitREQ_PN

#endif
