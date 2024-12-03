
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

struct AddedData{
    int    nep=0;
    Vec2i* bs=0;
    Vec3d* dirs=0;
    int*  isep=0;
};

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
    Quat4d*    typeREQs =0;   // [ntype] parameters for each type
    Quat4d*    typeREQsMin=0; // [ntype] equlibirum value of parameters for regularization 
    Quat4d*    typeREQsMax=0; // [ntype] equlibirum value of parameters for regularization 
    Quat4d*    typeREQs0=0;   // [ntype] equlibirum value of parameters for regularization
    Quat4d*    typeKreg =0;   // [ntype] regulatization stiffness
    Quat4d*    typeREQs0_low =0;   // [ntype] equlibirum value of parameters for regularization (lower wall)
    Quat4d*    typeKreg_low  =0;   // [ntype] regulatization stiffness (lower wall)
    Quat4d*    typeREQs0_high=0;   // [ntype] equlibirum value of parameters for regularization (upper wall)
    Quat4d*    typeKreg_high =0;   // [ntype] regulatization stiffness (upper wall)
    Quat4i*    typToREQ =0;   // [ntype] map each unique atom type to place in DOFs;
    
    double*   DOFs =0;       // [nDOFs]
    double*   fDOFs=0;       // [nDOFs]
    double*   vDOFs=0;       // [nDOFs]
    double*   DOFs_old =0;       // [nDOFs]
    double*   fDOFs_old=0;       // [nDOFs]
    bool      bEpairs= false;

    // parameters
    double Kneutral      = 1.0;
    //double max_step      = 0.1345;
    //double max_step      = 0.05; // maximal variation (in percentage) of any parameters allowed in one step during optimization
    double Rfac_tooClose = 0.0;
    double kMorse        = 1.6;

    //double* Rs =0,Es =0,Qs =0;
    //double* fRs=0,fEs=0,fQs=0;

    //std::vector(Atoms) examples; // Training examples
    Atoms*   batch=0;     // [nbatch]  // ToDo: would be more convenient to store Atoms* rather than Atoms
    double*  weights = 0; // [nbatch] scaling importaince of parameters
    double*  Es      = 0; 
    Mat3d*   poses   = 0; // [nbatch]
   
    // for rigid fitting
    Atoms* systemTest0=0; //[1]
    Atoms* systemTest =0; //[1]
    
    // system 0
    Atoms* system0=0;   // [1]

    // Temporary array for accumulation of derivs
    int     nmax = 0;
    Quat4d* fs   = 0; //[nmax] 
    //std::vector<Vec3d> fs;

    std::vector<Atoms> batch_vec;    // ToDo: would be more convenient to store Atoms* rather than Atoms
    std::vector<Atoms*> samples;    // ToDo: would be more convenient to store Atoms* rather than Atoms


    MM::Builder builder;

    MMFFparams* params=0; 

    // ------- Arrays for decomposition of energy components
    bool bDecomp = false;
    Quat4d** Elines  = 0;
    Quat4d** Eilines = 0;
    Quat4d** Ejlines = 0;

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
void tryRealocSamp(){
    int n=nmax;
    for(int i=0; i<samples.size(); i++){ int ni = samples[i]->natoms; if(ni>n){ n=ni;} }
    if(n>nmax){ _realloc( fs, n );  nmax=n; };
}

/**
 * Load a file of types involved in the parameter fitting. The file should contain lines with the following format:
 * atom_name mask_RvdW mask_EvdW mask_Q mask_Hb
 * where mask_RvdW, mask_EvdW, mask_Q, and mask_Hb are integers representing the indices of the fitting parameters for the RvdW, EvdW, Q, and Hb parameters, respectively.
 * @param fname The name of the file to load.
 * @return The number of types loaded.
*/
int loadTypeSelection( const char* fname, int imodel ){
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    int ntypesel = 0;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        Quat4i tm;
        int nw = sscanf( line, "%s %i %i %i %i\n", at_name, &tm.x,&tm.y,&tm.z,&tm.w ); 
        if(nw==5){ ntypesel++; }
    }
    fseek( fin, 0, SEEK_SET );
    std::vector<int> t;
    std::vector<Quat4i> typeMask;
    //std::vector<Quat4d> tREQs;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        Quat4i tm;
        int nw = sscanf( line, "%s %i %i %i %i\n", at_name, &tm.x,&tm.y,&tm.z,&tm.w ); 
        if(nw!=5)continue;
        //if(tm.z!=0){printf("ERROR: FitREQ::loadTypeSelection() mask_Q should be 0 for now\n"); exit(0); }
        typeMask.push_back( tm );
        int ityp = params->getAtomType(at_name);
        t.push_back( ityp );
        AtomType& t  = params->atypes[ityp];  
        //tREQs.push_back( Quat4d{ t.RvdW, t.EvdW, t.Qbase, t.Hb } );
    }
    fclose(fin);
    //init_types( ntypesel, &typeMask[0], &tREQs[0], true );
    init_types_new( ntypesel, &t[0], &typeMask[0], imodel );
    return ntypesel;
}

int loadTypeSelection_walls( const char* fname ){
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    int ntypesel = 0;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        Quat4i tm;
        Quat4d tx0l, tx0h;
        Quat4d tkl, tkh;
        int nw = sscanf( line, "%s %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
        at_name, &tm.x,&tm.y,&tm.z,&tm.w, &tx0l.x,&tkl.x,&tx0h.x,&tkh.x, &tx0l.y,&tkl.y,&tx0h.y,&tkh.y, &tx0l.z,&tkl.z,&tx0h.z,&tkh.z, &tx0l.w,&tkl.w,&tx0h.w,&tkh.w ); 
        if(nw==21){ ntypesel++; }
    }
    fseek( fin, 0, SEEK_SET );
    std::vector<int> t;
    std::vector<Quat4i> typeMask;
    std::vector<Quat4d> ts0_high, ts0_low;
    std::vector<Quat4d> tk_high, tk_low;
    while( fgets(line, nline, fin) ){
        if(line[0]=='#')continue;
        Quat4i tm;
        Quat4d tx0l, tx0h, tx0m;
        Quat4d tkl, tkh, tkm;
        int nw = sscanf( line, "%s %i %i %i %i %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n", 
        at_name, &tm.x,&tm.y,&tm.z,&tm.w, &tx0l.x,&tkl.x,&tx0h.x,&tkh.x, &tx0l.y,&tkl.y,&tx0h.y,&tkh.y, &tx0l.z,&tkl.z,&tx0h.z,&tkh.z, &tx0l.w,&tkl.w,&tx0h.w,&tkh.w ); 
        if(nw!=21)continue;
        //if(tm.z!=0){printf("ERROR: FitREQ::loadTypeSelection() mask_Q should be 0 for now\n"); exit(0); }
        typeMask.push_back( tm );
        ts0_low.push_back( tx0l );
        ts0_high.push_back( tx0h );
        tk_low.push_back( tkl );
        tk_high.push_back( tkh );
        int ityp = params->getAtomType(at_name);
        t.push_back( ityp );
        AtomType& t  = params->atypes[ityp];  
        //tREQs.push_back( Quat4d{ t.RvdW, t.EvdW, t.Qbase, t.Hb } );
    }
    fclose(fin);
    //init_types( ntypesel, &typeMask[0], &tREQs[0], true );
    //init_types_new( ntypesel, &t[0], &typeMask[0], imodel );
    init_types_walls( ntypesel, &t[0], &typeMask[0], &ts0_low[0], &tk_low[0], &ts0_high[0], &tk_high[0]);
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
int init_types_new( int ntypesel, int* tsel, Quat4i* typeMask, int imodel ){
    int ntype_ = params->atypes.size();
    printf( "FitREQ::init_types_new() ntypesel=%i ntypes=%i \n", ntypesel, ntype_ );
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
                printf(  "init_types_new() [%i] typeMask(%i,%i,%i,%i) typToREQ(%i,%i,%i,%i) nDOF %i\n", i, tm.x,tm.y,tm.z,tm.w, tt.x,tt.y,tt.z,tt.w, nDOFs );
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
    typeREQs    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs[i]    = Quat4d{ params->atypes[i].RvdW, params->atypes[i].EvdW, params->atypes[i].Qbase, params->atypes[i].Hb }; }
    if((imodel==2)||(imodel==5)||(imodel==8)){
        typeREQsMin = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Quat4d{ 0.0,    0.0,   -1e+300, -1.0 }; }
        typeREQsMax = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Quat4d{ 1e+300, 1e+300, 1e+300,  1.0 }; }
    }else{
        typeREQsMin = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Quat4d{ 0.0,    0.0,   -1e+300, -1e+300 }; }
        typeREQsMax = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Quat4d{ 1e+300, 1e+300, 1e+300,  1e+300 }; }
    }
    typeREQs0   = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0[i]   = typeREQs[i]; }
    typeKreg    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg[i]    = Quat4dZero; }
    DOFsFromTypes(); 
    return nDOFs;
}

int init_types_walls( int ntypesel, int* tsel, Quat4i* typeMask, Quat4d* ts0_low, Quat4d* tk_low, Quat4d* ts0_high, Quat4d* tk_high ){
    int ntype_ = params->atypes.size();
    printf( "FitREQ::init_types_new() ntypesel=%i ntypes=%i \n", ntypesel, ntype_ );
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
                printf(  "init_types_new() [%i] typeMask(%i,%i,%i,%i) typToREQ(%i,%i,%i,%i) nDOF %i\n", i, tm.x,tm.y,tm.z,tm.w, tt.x,tt.y,tt.z,tt.w, nDOFs );
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
    typeREQs    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs[i]    = Quat4d{ params->atypes[i].RvdW, params->atypes[i].EvdW, params->atypes[i].Qbase, params->atypes[i].Hb }; }
    typeREQsMin = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Quat4d{ -1e300, -1e300, -1e300, -1e300 }; }
    typeREQsMax = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Quat4d{ 1e+300, 1e+300, 1e+300, 1e+300 }; }
    typeREQs0_low  = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0_low [i] = Quat4dZero; }
    typeKreg_low   = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg_low  [i] = Quat4dZero; }
    typeREQs0_high = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0_high[i] = Quat4dZero; }
    typeKreg_high  = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg_high [i] = Quat4dZero; }
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
    DOFsFromTypes(); 
    return nDOFs;
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
int loadXYZ_new( const char* fname, bool bAddEpairs=false, bool bOutXYZ=false ){
    printf( "FitREQ::loadXYZ_new()\n" );
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
            if(bAddEpairs){
                // add Epairs
                Atoms* bak = atoms; 
                int n0bak  = bak->n0;
                int natbak = bak->natoms;
                int n1bak  = natbak - n0bak;
//for(int i=0; i<atoms->natoms; i++){printf( "000 atoms[%i] %s pos=%g %g %g q=%g\n", i, params->atypes[atoms->atypes[i]].name, atoms->apos[i].x, atoms->apos[i].y, atoms->apos[i].z, atoms->charge[i] );};printf("\n");
                atoms = addEpairs( atoms );
                atoms->n0 = bak->n0;
                atoms->Energy = bak->Energy;
                for(int i=0; i<natbak; i++){ atoms->charge[i] = bak->charge[i]; }
                delete bak;
                // find root atoms and directions
                Vec2i* bs=0;
                Vec3d* dirs=0;
                int nep_found = builder.listEpairBonds( bs, dirs );
                for(int i=0; i < nep_found; i++){dirs[i].normalize();};
                int*  isep=0;
                isep = new int[atoms->natoms];
                for(int i=0; i < natbak; i++){ isep[i]=0; }
                for(int i=natbak; i < atoms->natoms; i++){ isep[i]=1; } // Epairs are after atoms
//for(int i=0; i<atoms->natoms; i++){printf( "BEFORE atoms[%i] %s pos=%g %g %g q=%g\n", i, params->atypes[atoms->atypes[i]].name, atoms->apos[i].x, atoms->apos[i].y, atoms->apos[i].z, atoms->charge[i] );};
//for(int i=0; i<nep_found; i++){printf( "BEFORE Epair[%i] %i %i %g %g %g\n", i, bs[i].x, bs[i].y, dirs[i].x, dirs[i].y, dirs[i].z );};
//printf("\n");
                // shuffle atoms so that they are ordered as Atoms(mol1), Epairs(mol1), Atoms(mol2), Epairs(mol2)
                Atoms bak2; bak2.copyOf( *atoms ); 
                std::vector<Vec2i> bsbak( bs, bs+nep_found );
                std::vector<Vec3d> dirsbak( dirs, dirs+nep_found );
                // Atoms(mol1): nothing to do
                // Epairs(mol1)
                int iE=-1;
                for(int i=0; i<nep_found; i++){ 
                    if(bsbak[i].x<n0bak){
                        iE++;
                        int j=n0bak+iE;
                        int k=bsbak[i].y;
                        atoms->apos[j]   = bak2.apos[k];
                        atoms->atypes[j] = bak2.atypes[k];
                        atoms->charge[j] = bak2.charge[k];
                        atoms->n0++;
                        bs[iE].x = bsbak[i].x;
                        bs[iE].y = j;
                        dirs[iE] = dirsbak[i];
                        isep[j] = 1;
                    }
                }
                int nE0 = iE+1;
                int nE1 = atoms->natoms - natbak - nE0;
                // Atoms(mol2)
                for(int i=0; i<n1bak; i++){
                    int j=atoms->n0+i;
                    int k=n0bak+i;
                    atoms->apos[j]   = bak2.apos[k];
                    atoms->atypes[j] = bak2.atypes[k];
                    atoms->charge[j] = bak2.charge[k];
                    isep[j] = 0;
                }
                // Epairs(mol2)
                iE=-1;
                for(int i=0; i<nep_found; i++){ 
                    if(bsbak[i].x>=n0bak){
                        iE++;
                        int j=n0bak+nE0+n1bak+iE;
                        int k=bsbak[i].y;
                        atoms->apos[j]   = bak2.apos[k];
                        atoms->atypes[j] = bak2.atypes[k];
                        atoms->charge[j] = bak2.charge[k];
                        bs[nE0+iE].x = bsbak[i].x+nE0;
                        bs[nE0+iE].y = j;
                        dirs[nE0+iE] = dirsbak[i];
                        isep[j] = 1;
                    }
                }
                //delete bak2; 
//for(int i=0; i<atoms->natoms; i++){printf( "AFTER atoms[%i] %s pos=%g %g %g q=%g isep=%d\n", i, params->atypes[atoms->atypes[i]].name, atoms->apos[i].x, atoms->apos[i].y, atoms->apos[i].z, atoms->charge[i], isep[i] );};
//for(int i=0; i<nep_found; i++){printf( "AFTER Epair[%i] %i %i %g %g %g\n", i, bs[i].x, bs[i].y, dirs[i].x, dirs[i].y, dirs[i].z );};
//printf("FIN QUA\n");exit(0);
//if(nbatch==0){for(int i=0; i<nep_found; i++){printf( "AFTER Epair[%i] %i %i %g %g %g\n", i, bs[i].x, bs[i].y, dirs[i].x, dirs[i].y, dirs[i].z );};}
                // store root atoms and directions
                AddedData * ad = new AddedData();
                ad->nep  = nep_found;
                ad->bs   = bs;
                ad->dirs = dirs;
                ad->isep = isep;
                atoms->userData = ad;
                if(fout){
                    sprintf(line,"#	n0 %i E_tot %g", atoms->n0, atoms->Energy ); 
                    params->writeXYZ( fout, atoms, line, 0, true );
                }
            }
            samples.push_back( atoms );
            il=0; nbatch++;
        }
    }
    if(fout)fclose(fout);
    fclose(fin);
    //init_types_new();
    return nbatch;
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
 * @brief reads non-colvalent interaction parameter REQH(Rvdw,Evdw,Q,Hb) of given atom-type from aproprieate degrees of freedom (DOF) according to index stored in typToREQ[ityp]. If index of DOF is negative, the parameter is not read.  
 * @param ityp Index of the type
 */
inline void DOFsToType(int ityp){
    const Quat4i& tt = typToREQ[ityp];
    Quat4d& REQ      = typeREQs[ityp];
    if(tt.x>=0)REQ.x = DOFs[tt.x];
    if(tt.y>=0)REQ.y = DOFs[tt.y];
    if(tt.z>=0)REQ.z = DOFs[tt.z];
    if(tt.w>=0)REQ.w = DOFs[tt.w];
}
void DOFsToTypes(){ for(int i=0; i<ntype; i++ ){ DOFsToType(i); } }
void getType(int i, Quat4d& REQ ){ typeREQs[i]=REQ; DOFsToType(i); }

/**
 * @brief writes non-colvalent interaction parameter REQH(Rvdw,Evdw,Q,Hb) of given atom-type from aproprieate degrees of freedom (DOF) according to index stored in typToREQ[ityp]. If index of DOF is negative, the parameter is not writen.  
 * 
 * @param ityp The index of the type.
 */
inline void DOFsFromType(int ityp){
    const Quat4i& tt  = typToREQ[ityp];
    const Quat4d& REQ = typeREQs[ityp];
    if(tt.x>=0)DOFs[tt.x] = REQ.x;
    if(tt.y>=0)DOFs[tt.y] = REQ.y;
    if(tt.z>=0)DOFs[tt.z] = REQ.z;
    if(tt.w>=0)DOFs[tt.w] = REQ.w;
}
void DOFsFromTypes(){ for(int i=0; i<ntype; i++ ){ DOFsFromType(i); } }
void setType(int i, Quat4d REQ ){ typeREQs[i]=REQ; DOFsFromType(i); }

// ======================================
// =========  EVAL DERIVS  ==============
// ======================================

void clean_fs(int n){ for(int i=0; i<n; i++){fs[i]=Quat4dZero;} }

/**
 * @brief Evaluates the variational derivatives of the fitting error (sumed batch training samples) with respect to all fitting parameters (i.e. non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) for each atom type). The atomic systems are not assumed rigid (i.e. the atomic positions can vary freely from sample to sample).
 * 
 * @param Eout array to store the non-covalent interaction energy values of each atomic system in the batch. if Eout==null, the function will not store the energy values.
 * @return double, returns the total fitting error.
 */
double evalDerivsSamp( double* Eout=0 ){ 
    //printf( "FitREQ::evalDerivsSamp() nbatch %i imodel %i verbosity %i \n", samples.size(), imodel, verbosity );
    tryRealocSamp();
    double Error = 0.0;
    std::vector<double> qs_;
    std::vector<Vec3d>  apos_; 
    std::vector<int>    isep_;
    std::vector<Vec3d>  dirs_; 
    for(int i=0; i<samples.size(); i++){
        const Atoms* atoms = samples[i];
        qs_.resize(atoms->natoms);
        apos_.resize(atoms->natoms);
        isep_.resize(atoms->natoms);
        dirs_.resize(atoms->natoms);
        for(int j=0; j<atoms->natoms; j++){
            qs_[j]   = atoms->charge[j];
            apos_[j] = atoms->apos[j];
            isep_[j] = 0;
            dirs_[j] = Vec3dZero;
        }
        // electron pairs
        int nep=0;
        Vec2i* bs   = 0;
        if( bEpairs ){
            int*  isep = 0;
            Vec3d* dirs = 0;
            AddedData * ad = (AddedData*)atoms->userData;
            nep  = ad->nep;
            bs   = ad->bs;
            dirs = ad->dirs;
            isep = ad->isep;
            // set values of the Epairs' and root atoms' charges according to the current values of the parameters
            for(int j=0; j<nep; j++){
                int iX=bs[j].x;
                int iE=bs[j].y;
                double q  = typeREQs[atoms->atypes[iE]].z;
                qs_[iE]   = q;
                qs_[iX]   = qs_[iX] - q;
                apos_[iE] = apos_[iX] + dirs[j] * typeREQs[atoms->atypes[iE]].w;
                dirs_[iE] = dirs[j];
                isep_[iE] = 1;
            }
//for(int j=0; j<nep; j++){printf( "dirs[%i] %g %g %g\n", j, dirs[j].x, dirs[j].y, dirs[j].z );};
//for(int j=0; j<atoms->natoms; j++){printf( "atom[%i] %g %g %g\n", j, dirs_[j].x, dirs_[j].y, dirs_[j].z );};
//exit(0);
        }
        int     nj      = atoms->n0;
        int*    jtyp    = atoms->atypes;
        //Vec3d*  jpos  = atoms->apos;
        Vec3d*  jpos    = apos_.data();
        double* jq      = qs_.data();
        int*    jisep   = isep_.data();
        Vec3d*  jdirs   = dirs_.data();
        int     na      = atoms->natoms - atoms->n0;
        int*    atypes  = atoms->atypes + atoms->n0;
        //Vec3d* apos   = atoms->apos   + atoms->n0;
        Vec3d*  apos    = apos_.data()  + atoms->n0;
        double* aq      = qs_.data()    + atoms->n0;
        int*    aisep   = isep_.data()  + atoms->n0;
        Vec3d*  adirs   = dirs_.data()  + atoms->n0;
        double Eref=atoms->Energy;
        double wi = 1.0; 
        if(weights) wi = weights[i];
        // ToDo: we need to initialize fs according to DOFs before calling clean_fs()
        //clean_fs(na);
        clean_fs(atoms->natoms);
        double E;

//for(int i=0; i<na; i++){printf( "atoms[%i] %s pos=(%g %g %g) q=%g isep=%i dirs=(%g %g %g)\n", i, params->atypes[atypes[i]].name, apos[i].x, apos[i].y, apos[i].z, aq[i], aisep[i], adirs[i].x, adirs[i].y, adirs[i].z );};
//printf("\n");
//for(int j=0; j<nj; j++){printf( "atoms[%i] %s pos=(%g %g %g) q=%g isep=%i dirs=(%g %g %g)\n", j, params->atypes[jtyp[j]].name, jpos[j].x, jpos[j].y, jpos[j].z, jq[j], jisep[j], jdirs[j].x, jdirs[j].y, jdirs[j].z );};
//printf("\n");
//for(int i=0; i<nep; i++){printf( "Epair[%i] %i %i\n", i, bs[i].x, bs[i].y );};
//exit(0);
        switch (imodel){
            case 0:  E = evalExampleDerivs_LJQ        (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 1:  E = evalExampleDerivs_LJQH1      (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 2:  E = evalExampleDerivs_LJQH2      (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 3:  E = evalExampleDerivs_LJQH1H2    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 4:  E = evalExampleDerivs_BuckQ      (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 5:  E = evalExampleDerivs_BuckQH1    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 6:  E = evalExampleDerivs_BuckQH2    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 7:  E = evalExampleDerivs_BuckQH1H2  (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 8:  E = evalExampleDerivs_MorseQ     (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 9:  E = evalExampleDerivs_MorseQH1   (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 10: E = evalExampleDerivs_MorseQH2   (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 11: E = evalExampleDerivs_MorseQH1H2 (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 12: E = evalExampleDerivs_LJx2Q      (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 13: E = evalExampleDerivs_LJx2QH1    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 14: E = evalExampleDerivs_LJx2QH2    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 15: E = evalExampleDerivs_LJx2QH1H2  (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 16: E = evalExampleDerivs_LJr8Q      (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 17: E = evalExampleDerivs_LJr8QH1    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 18: E = evalExampleDerivs_LJr8QH2    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 19: E = evalExampleDerivs_LJr8QH1H2  (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 20: E = evalExampleDerivs_LJr9Q      (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 21: E = evalExampleDerivs_LJr9QH1    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break; 
            case 22: E = evalExampleDerivs_LJr9QH2    (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
            case 23: E = evalExampleDerivs_LJr9QH1H2  (na, atypes, apos, aq, aisep, adirs, nj, jtyp, jpos, jq, jisep, jdirs, nep, bs); break;
        }
        if( E>1e+300 ){
            if(verbosity>0) printf( "skipped sample [%i] atoms too close \n", i );
            continue;
        } 
        if(Eout){ Eout[i]=E; };
        double dE = (E - Eref);
        Error += dE*dE*wi;
        acumDerivs( atoms->natoms, jtyp, 2.0*dE*wi );
    }
    return Error;
}

/**
 * Accumulates the derivatives of the non-covalent interaction energy with respect to fitting parameters in the array fDOFs.
 * @param n The number of atoms in the system.
 * @param types An array of integers representing the types of atoms in the system.
 * @param dE The difference between the energy of the system and the reference energy (scalled by the weight for this sample (system)
*/
void acumDerivs( int n, int* types, double dEw){
    
    for(int i=0; i<n; i++){
        int t            = types[i];
        const Quat4i& tt = typToREQ[t];
        const Quat4d  f  = fs[i];
        if(tt.x>=0)fDOFs[tt.x]+=f.x*dEw;
        if(tt.y>=0)fDOFs[tt.y]+=f.y*dEw;
        if(tt.z>=0)fDOFs[tt.z]+=f.z*dEw;
        if(tt.w>=0)fDOFs[tt.w]+=f.w*dEw;
        //if((tt.y>=0)&&(i==0))printf( "acumDerivs i= %i f= %g f*dEw= %g fDOFs= %g\n", i, f.y, f.y*dEw, fDOFs[tt.y] );
        //if((tt.x>=0)||(tt.y>=0))printf( "acumDerivs i= %i dE= %g f.x= %g fDOFs= %g f.y= %g fDOFs= %g\n", i, dE, f.x, fDOFs[tt.x], f.y, fDOFs[tt.y] );
    }
    //exit(0);
}

/**
 * @brief Calculates the correction to the electrostatic energy and its derivative with respect to the charge.
 */
double corr_elec( double ir, double ir2, double Q, Vec3d d, Vec3d* dirs, int i, int nep, Vec2i* bs, int nj, Vec3d* pos, int j, Vec3d* ps, double &dE_dQ){
    double dE_dr = ir * ir2 * COULOMB_CONST * Q * d.dot(dirs[i]);
    for(int k=0; k<nep; k++){
        if(bs[k].y==nj+i){
            Vec3d dd = pos[j] - ps[bs[k].x-nj];
            dE_dQ -= COULOMB_CONST / sqrt( dd.norm2() );
            break;
        }
    }
    return dE_dr;
}

/**
 * Calculates the total energy and forces of a system of atoms interacting through the Lennard-Jones potential and electrostatic interactions.
 * 
 * @param n The number of atoms in system(=second molecule).
 * @param types An array of integers representing the type of each atom in system(=second molecule).
 * @param ps An array of Vec3d representing the position of each atom in system(=second molecule).
 * @param aq An array of doubles representing the charge of each atom in system(=second molecule).
 * @param nj The number of atoms in the system0(=first molecule).
 * @param jtyp An array of integers representing the type of each atom in system0(=first molecule).
 * @param jpos An array of Vec3d representing the position of each atom in system0(=first molecule).
 * @param jq An array of doubles representing the charge of each atom in system0(=first molecule).
 * 
 * @return The total energy of the system.
 */
double evalExampleDerivs_LJQ( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJQ (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double dE_deps = u6 * ( u6 - 2.0 );
            double dE_dR0  = 12.0 * eps / (R0+1e-300) * u6 * ( u6 - 1.0 );
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y; // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                     // dEtot/dQi
                fsi.w = dE_dri;                            // dEtot/dri
            }
            fsj.x = dE_dR0;                                // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y; // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                     // dEtot/dQj
                fsj.w = dE_drj;                            // dEtot/drj
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJQH1( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJQH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double h1p     = 1.0 + H1;
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double dE_deps = u6 * ( u6 - 2.0 * h1p );
            double dE_dR0  = 12.0 * eps / (R0+1e-300) * u6 * (  u6 - h1p );
            double dE_dH1  = -2.0 * eps * u6;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJQH2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJQH2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double u6p     = ( 1.0 - H2 ) * u6;                     
            double dE_deps = u6 * ( u6p - 2.0 );
            double dE_dR0  = 12.0 * eps / (R0+1e-300) * u6 * ( u6p - 1.0 );
            double dE_dH2  = -eps * u6 * u6;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJQH1H2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJQH1H2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double h1p     = 1.0 + H1;
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double u6p     = ( 1.0 - H2 ) * u6;                     
            double dE_deps = u6 * ( u6p - 2.0 * h1p );
            double dE_dR0  = 12.0 * eps / (R0+1e-300) * u6 * ( u6p - h1p );
            double dE_dH1  = -2.0 * eps * u6;
            double dE_dH2  = -eps * u6 * u6;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

/* The Buckingham potential can be written in the form:
   E = epsilon * { 6/(alpha-6) * exp[alpha*(1-r/R0)] - alpha/(alpha-6) * (R0/r)^6 }
   alpha is a free parameter; it can be chosen such as:
   - alpha = 12, which corresponds to the same dispersion part of the LJ potential
   - alpha = 19/2 + sqrt(73)/2, which corresponds to the same curvature at the minimum of the LJ potential
   - alpha is fitted to reference data (not implemented yet)*/
double evalExampleDerivs_BuckQ( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot  = 0.0;
    double alpha = 12.0; //double alpha = 19.0/2.0 + sqrt(73.0)/2.0;
    double ia6   = 1.0 / ( alpha - 6.0 );
    //printf( "evalExampleDerivs_BuckQ (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u1 = ir * R0;
            double u2 = u1 * u1;
            double u4 = u2 * u2;
            double u6 = u4 * u2;
            // --- Buckingham
            double iu1     = 1.0 / u1;
            double e       = exp( alpha * ( 1.0 - iu1 ) );
            double dE_deps = ia6 * ( 6.0 * e - alpha * u6 );
            double dE_dR0  = 6.0 * eps * alpha / (R0+1e-300) * ia6 * ( iu1 * e - u6 );
            double EBuck   = eps * dE_deps;
            // --- Energy and forces
            Etot  += EBuck + Eel;
            fsi.x = dE_dR0;                                // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y; // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                     // dEtot/dQi
                fsi.w = dE_dri;                            // dEtot/dri
            }
            fsj.x = dE_dR0;                                // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y; // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                     // dEtot/dQj
                fsj.w = dE_drj;                            // dEtot/drj
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_BuckQH1( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot  = 0.0;
    double alpha = 12.0; //double alpha = 19.0/2.0 + sqrt(73.0)/2.0;
    double ia6   = 1.0 / ( alpha - 6.0 );
    //printf( "evalExampleDerivs_BuckQH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u1  = ir * R0;
            double u2  = u1 * u1;
            double u4  = u2 * u2;
            double u6  = u4 * u2;
            double u6p = ( 1.0 + H1 ) * u6;                     
            // --- Buckingham
            double iu1     = 1.0 / u1;
            double e       = exp( alpha * ( 1.0 - iu1 ) );
            double dE_deps = ia6 * ( 6.0 * e - alpha * u6p );
            double dE_dR0  = 6.0 * eps * alpha * ia6 / (R0+1e-300) * ( iu1 * e - u6p );
            double dE_dH1  = -eps * alpha * ia6 * u6;
            double EBuck   = eps * dE_deps;
            // --- Energy and forces
            Etot  += EBuck + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_BuckQH2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot  = 0.0;
    double alpha = 12.0; //double alpha = 19.0/2.0 + sqrt(73.0)/2.0;
    double ia6   = 1.0 / ( alpha - 6.0 );
    //printf( "evalExampleDerivs_BuckQH2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u1 = ir * R0;
            double u2 = u1 * u1;
            double u4 = u2 * u2;
            double u6 = u4 * u2;
            // --- Buckingham
            double iu1     = 1.0 / u1;
            double e       = exp( alpha * ( 1.0 - iu1 ) );
            double ep      = e * ( 1.0 - H2 );
            double dE_deps = ia6 * ( 6.0 * ep - alpha * u6 );
            double dE_dR0  = 6.0 * eps * alpha * ia6 / (R0+1e-300) * ( iu1 * ep - u6 );
            double dE_dH2  = -6.0 * eps * ia6 * e;
            double EBuck   = eps * dE_deps;
            // --- Energy and forces
            Etot  += EBuck + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_BuckQH1H2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot  = 0.0;
    double alpha = 12.0; //double alpha = 19.0/2.0 + sqrt(73.0)/2.0;
    double ia6   = 1.0 / ( alpha - 6.0 );
    //printf( "evalExampleDerivs_BuckQH1H2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u1  = ir * R0;
            double u2  = u1 * u1;
            double u4  = u2 * u2;
            double u6  = u4 * u2;
            double u6p = ( 1.0 + H1 ) * u6;                     
            // --- Buckingham
            double iu1     = 1.0 / u1;
            double e       = exp( alpha * ( 1.0 - iu1 ) );
            double ep      = e * ( 1.0 - H2 );
            double dE_deps = ia6 * ( 6.0 * ep - alpha * u6p );
            double dE_dR0  = 6.0 * eps * alpha * ia6 / (R0+1e-300) * ( iu1 * ep - u6p );
            double dE_dH1  = -eps * alpha * ia6 * u6;
            double dE_dH2  = -6.0 * eps * ia6 * e;
            double EBuck   = eps * dE_deps;
            // --- Energy and forces
            Etot  += EBuck + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

/* The Morse potential can be written in the form:
   E = epsilon * { exp[-2*alpha*(r-R0)] - 2 * exp[-alpha*(r-R0)] }
   so that it goes to zero as r goes to infinity
   alpha is a free parameter; it can be chosen such as:
   - alpha = 6/R0, which corresponds to the same curvature at the minimum of the LJ potential
   - have the dispersion part as close as possible to LJ (not derived yet)
   - alpha is fitted to reference data (not implemented yet)*/
double evalExampleDerivs_MorseQ( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_MorseQ (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Morse
            double alpha   = 6.0 / R0;
            double r       = 1.0 / ir;
            double e       = exp( -alpha * ( r - R0 ) );
            double e2      = e * e;
            double dE_deps = e2 - 2.0 * e;
            double dE_dR0  = 2.0 * alpha * eps * ( e2 - e );
            double EMorse  = eps * dE_deps;
            // --- Energy and forces
            Etot  += EMorse + Eel;
            fsi.x = dE_dR0;                                // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y; // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                     // dEtot/dQi
                fsi.w = dE_dri;                            // dEtot/dri
            }
            fsj.x = dE_dR0;                                // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y; // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                     // dEtot/dQj
                fsj.w = dE_drj;                            // dEtot/drj
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_MorseQH1( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_MorseQH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Morse
            double alpha   = 6.0 / R0;
            double r       = 1.0 / ir;
            double e       = exp( -alpha * ( r - R0 ) );
            double e2      = e * e;
            double ep      = ( 1.0 + H1 ) * e;
            double dE_deps = e2 - 2.0 * ep;
            double dE_dR0  = 2.0 * alpha * eps * ( e2 - ep );
            double dE_dH1  = -2.0 * eps * e;
            double EMorse  = eps * dE_deps;
            // --- Energy and forces
            Etot  += EMorse + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_MorseQH2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_MorseQH2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Morse
            double alpha   = 6.0 / R0;
            double r       = 1.0 / ir;
            double e       = exp( -alpha * ( r - R0 ) );
            double e2      = e * e;
            double e2p     = ( 1.0 - H2 ) * e2;
            double dE_deps = e2p - 2.0 * e;
            double dE_dR0  = 2.0 * alpha * eps * ( e2p - e );
            double dE_dH2  = - eps * e2;
            double EMorse  = eps * dE_deps;
            // --- Energy and forces
            Etot  += EMorse + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_MorseQH1H2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_MorseQH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Morse
            double alpha   = 6.0 / R0;
            double r       = 1.0 / ir;
            double e       = exp( -alpha * ( r - R0 ) );
            double e2      = e * e;
            double ep      = ( 1.0 + H1 ) * e;
            double e2p     = ( 1.0 - H2 ) * e2;
            double dE_deps = e2p - 2.0 * ep;
            double dE_dR0  = 2.0 * alpha * eps * ( e2p - ep );
            double dE_dH1  = -2.0 * eps * e;
            double dE_dH2  = -eps * e2;
            double EMorse  = eps * dE_deps;
            // --- Energy and forces
            Etot  += EMorse + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "Etot= %g\n", Etot );exit(0);    
    return Etot;
}

// Modified LJ potential with a quadratic term to model repulsion
double evalExampleDerivs_LJx2Q( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJx2Q (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones or parabola
            double r = 1.0 / ir;
            double dE_deps;
            double dE_dR0;
            if(r<R0){
                double s1 = r / R0 - 1.0;
                dE_deps   = 36.0 * s1 * s1 - 1.0;
                dE_dR0    = -72.0 * eps * r / ( R0 * R0 ) * s1;
            }else{
                double u2 = ir2 * ( R0 * R0 );
                double u4 = u2 * u2;
                double u6 = u4 * u2;
                dE_deps   = u6 * ( u6 - 2.0 );
                dE_dR0    = 12.0 * eps / (R0+1e-300) * u6 * ( u6 - 1.0 );
            }
            double ELJ = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y; // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                     // dEtot/dQi
                fsi.w = dE_dri;                            // dEtot/dri
            }
            fsj.x = dE_dR0;                                // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y; // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                     // dEtot/dQj
                fsj.w = dE_drj;                            // dEtot/drj
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJx2QH1( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJx2QH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones or parabola
            double h1p    = 1.0 + H1;
            double sqh1p  = sqrt( h1p );
            double sih1p  = cbrt( sqh1p );
            double isih1p = 1.0 / sih1p;
            double R0p    = R0 * isih1p;
            double r      = 1.0 / ir;
            double dE_deps;
            double dE_dR0;
            double dE_dH1;
            if(r<R0p){
                double h1p2 = h1p * h1p;
                double s    = r / R0p;
                double s1   = s - 1.0;
                double s14  = s - 13.0 / 14.0;
                dE_deps     = 36.0 * ( s1 * s1 - 1.0 ) * h1p2;
                dE_dR0      = -72.0 * eps * h1p2 * r / ( R0p * R0p ) * s1 * isih1p;
                dE_dH1      = 2.0 * eps * h1p * ( 42.0 * s14 * s14 - 17.0 / 14.0 );
            }else{
                double u2 = ir2 * ( R0 * R0 );
                double u4 = u2 * u2;
                double u6 = u4 * u2;
                dE_deps   = u6 * ( u6 - 2.0 * h1p );
                dE_dR0    = 12.0 * eps / (R0+1e-300) * u6 * ( u6 - h1p );
                dE_dH1    = -2.0 * eps * u6;
            }
            double ELJ = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJx2QH2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJx2QH2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones or parabola
            double h2p   = 1.0 - H2;
            double sqh2p = sqrt( h2p );
            double sih2p = cbrt( sqh2p );
            double R0p   = R0 * sih2p;
            double r     = 1.0 / ir;
            double dE_deps;
            double dE_dR0;
            double dE_dH2;
            if(r<R0p){
                double ih2p = 1.0 / h2p;
                double epsp = eps * ih2p;
                double s    = r / R0p;
                double s1   = s - 1.0;
                double s8   = s - 7.0 / 8.0;
                dE_deps     = 36.0 * ( s1 * s1 - 1.0 ) * ih2p;
                dE_dR0      = -72.0 * epsp * r / ( R0p * R0p ) * s1 * sih2p;
                dE_dH2      = epsp * ih2p * ( 48.0 * s8 * s8 - 7.0 / 4.0 );
            }else{
                double u2  = ir2 * ( R0 * R0 );
                double u4  = u2 * u2;
                double u6  = u4 * u2;
                double u6p = h2p * u6;                     
                dE_deps    = u6 * ( u6p - 2.0 );
                dE_dR0     = 12.0 * eps / (R0+1e-300) * u6 * ( u6p - 1.0 );
                dE_dH2     = -eps * u6 * u6;
            }
            double ELJ = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJx2QH1H2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJx2QH1H2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones or parabola
            double h1p  = 1.0 + H1;
            double h2p  = 1.0 - H2;
            double sqhp = sqrt( h2p / h1p );
            double sihp = cbrt( sqhp );
            double R0p  = R0 * sihp;
            double r    = 1.0 / ir;
            double dE_deps;
            double dE_dR0;
            double dE_dH1;
            double dE_dH2;
            if(r<R0p){
                double s    = r / R0p;
                double s1   = s - 1.0;
                double s14  = s - 13.0 / 14.0;
                double s8   = s - 7.0 / 8.0;
                double hp11 = h1p / h2p;
                double hp21 = hp11 * h1p;
                dE_deps     = ( 36.0 * s1 * s1 - 1.0 ) * hp21;
                dE_dR0      = -72.0 * eps * hp21 * r / ( R0p * R0p ) * s1 * sihp;
                dE_dH1      = 2.0 * eps * hp11 * ( 42.0 * s14 * s14 - 17.0 / 14.0 );
                dE_dH2      = eps * hp11 * hp11 * ( 48.0 * s8 * s8 - 7.0 / 4.0 );
            }else{
                double u2  = ir2 * ( R0 * R0 );
                double u4  = u2 * u2;
                double u6  = u4 * u2;
                double u6p = h2p * u6;                     
                dE_deps    = u6 * ( u6p - 2.0 * h1p );
                dE_dR0     = 12.0 * eps / (R0+1e-300) * u6 * ( u6p - h1p );
                dE_dH1     = -2.0 * eps * u6;
                dE_dH2     = -eps * u6 * u6;
            }
            double ELJ = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

// 8-6 Lennard-Jones potential
double evalExampleDerivs_LJr8Q( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr8Q (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double dE_deps = u6 * ( 3.0 * u2 - 4.0 );
            double dE_dR0  = 24.0 * eps / (R0+1e-300) * u6 * ( u2 - 1.0 );
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y; // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                     // dEtot/dQi
                fsi.w = dE_dri;                            // dEtot/dri
            }
            fsj.x = dE_dR0;                                // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y; // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                     // dEtot/dQj
                fsj.w = dE_drj;                            // dEtot/drj
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJr8QH1( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr8QH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double h1p     = 1.0 + H1;
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double dE_deps = u6 * ( 3.0 * u2 - 4.0 * h1p );
            double dE_dR0  = 24.0 * eps / (R0+1e-300) * u6 * (  u2 - h1p );
            double dE_dH1  = -4.0 * eps * u6;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJr8QH2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr8QH2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double u2p     = ( 1.0 - H2 ) * u2;                     
            double dE_deps = u6 * ( 3.0 * u2p - 4.0 );
            double dE_dR0  = 24.0 * eps / (R0+1e-300) * u6 * ( u2p - 1.0 );
            double dE_dH2  = -3.0 * eps * u6 * u2;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJr8QH1H2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr8QH1H2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double h1p     = 1.0 + H1;
            double u2      = ir2 * ( R0 * R0 );
            double u4      = u2 * u2;
            double u6      = u4 * u2;
            double u2p     = ( 1.0 - H2 ) * u2;                     
            double dE_deps = u6 * ( 3.0 * u2p - 4.0 * h1p );
            double dE_dR0  = 24.0 * eps / (R0+1e-300) * u6 * ( u2p - h1p );
            double dE_dH1  = -4.0 * eps * u6;
            double dE_dH2  = -3.0 * eps * u6 * u2;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

// 9-6 Lennard-Jones potential
double evalExampleDerivs_LJr9Q( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr9Q (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u1      = ir * R0;
            double u3      = u1 * u1 * u1;
            double u6      = u3 * u3;
            double dE_deps = u6 * ( 2.0 * u3 - 3.0 );
            double dE_dR0  = 18.0 * eps / (R0+1e-300) * u6 * ( u3 - 1.0 );
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y; // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                     // dEtot/dQi
                fsi.w = dE_dri;                            // dEtot/dri
            }
            fsj.x = dE_dR0;                                // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y; // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                     // dEtot/dQj
                fsj.w = dE_drj;                            // dEtot/drj
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJr9QH1( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr9QH1 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x; 
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j]; 
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double h1p     = 1.0 + H1;
            double u1      = ir * R0;
            double u3      = u1 * u1 * u1;
            double u6      = u3 * u3;
            double dE_deps = u6 * ( 2.0 * u3 - 3.0 * h1p );
            double dE_dR0  = 18.0 * eps / (R0+1e-300) * u6 * (  u3 - h1p );
            double dE_dH1  = -3.0 * eps * u6;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJr9QH2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr9QH2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double u1      = ir * R0;
            double u3      = u1 * u1 * u1;
            double u6      = u3 * u3;
            double u3p     = ( 1.0 - H2 ) * u3;                     
            double dE_deps = u6 * ( 2.0 * u3p - 3.0 );
            double dE_dR0  = 18.0 * eps / (R0+1e-300) * u6 * ( u3p - 1.0 );
            double dE_dH2  = -2.0 * eps * u6 * u3;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

double evalExampleDerivs_LJr9QH1H2( int n, int* types, Vec3d* ps, double* aq, int* aisep, Vec3d* adirs, int nj=0, int* jtyp=0, Vec3d* jpos=0, double* jq=0, int* jisep=0, Vec3d* jdirs=0, int nep=0, Vec2i* bs=0){
    double Etot = 0.0;
    //printf( "evalExampleDerivs_LJr9QH1H2 (n=%i,nj=%i)\n", n,nj );
    for(int i=0; i<n; i++){ // loop over all atoms[i] in system
        int   ti           = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        for(int j=0; j<nj; j++){ // loop over all atoms[j] in system0
            int tj             = jtyp[j];
            const Quat4d& REQj = typeREQs[tj];
            Vec3d d            = jpos[j] - pi;
            Quat4d fsi         = Quat4dZero;
            Quat4d fsj         = Quat4dZero;
            double R0          = REQi.x + REQj.x;
            double eps         = sqrt( REQi.y * REQj.y );
            double Q           = aq[i] * jq[j];
            double H1          = fmax( 0.0, -REQi.z * REQj.z );
            double H2          = fmax( 0.0, -REQi.w * REQj.w );
            // --- Electrostatic
            double ir2    = 1.0 / ( d.norm2() );
            double ir     = sqrt( ir2 );
            double dE_dQ  = ir * COULOMB_CONST;
            double Eel    = Q * dE_dQ;
            double dE_dri = 0.0;
            double dE_drj = 0.0;
            if(aisep[i]==1){dE_dri+=corr_elec(ir,ir2,Q,d,adirs,i,nep,bs,nj,jpos,j,ps,dE_dQ);}
            if(jisep[j]==1){dE_drj-=corr_elec(ir,ir2,Q,d,jdirs,j,nep,bs,0, ps,i,jpos,dE_dQ);}          
            // --- Lennard-Jones
            double h1p     = 1.0 + H1;
            double u1      = ir * R0;
            double u3      = u1 * u1 * u1;
            double u6      = u3 * u3;
            double u3p     = ( 1.0 - H2 ) * u3;                     
            double dE_deps = u6 * ( 2.0 * u3p - 3.0 * h1p );
            double dE_dR0  = 18.0 * eps / (R0+1e-300) * u6 * ( u3p - h1p );
            double dE_dH1  = -3.0 * eps * u6;
            double dE_dH2  = -2.0 * eps * u6 * u3;
            double ELJ     = eps * dE_deps;
            // --- Energy and forces
            Etot  += ELJ + Eel;
            fsi.x = dE_dR0;                                         // dEtot/dR0i
            fsi.y = dE_deps * 0.5 / (eps+1e-300) * REQj.y;          // dEtot/depsi
            if(aisep[i]==1){
                fsi.z = dE_dQ * jq[j];                              // dEtot/dQi
                fsi.w = dE_dri;                                     // dEtot/dri
            }else{
                fsi.z = dE_dH1 * H1 / (REQi.z+sign(REQi.z)*1e-300); // dEtot/dH1i
                fsi.w = dE_dH2 * H2 / (REQi.w+sign(REQi.w)*1e-300); // dEtot/dH2i
            }
            fsj.x = dE_dR0;                                         // dEtot/dR0j
            fsj.y = dE_deps * 0.5 / (eps+1e-300) * REQi.y;          // dEtot/depsj
            if(jisep[j]==1){
                fsj.z = dE_dQ * aq[i];                              // dEtot/dQj
                fsj.w = dE_drj;                                     // dEtot/drj
            }else{
                fsj.z = dE_dH1 * H1 / (REQj.z+sign(REQj.z)*1e-300); // dEtot/dH1j
                fsj.w = dE_dH2 * H2 / (REQj.w+sign(REQj.w)*1e-300); // dEtot/dH2j
            }
            fs[nj+i].add(fsi);
            fs[j].add(fsj);
//printf( "debug i= %i j= %i fsi.x= %g fsi.y= %g fsi.z= %g fsi.w= %g fsj.x= %g fsj.y= %g fsj.z= %g fsj.w= %g\n", i, j, fsi.x, fsi.y, fsi.z, fsi.w, fsj.x, fsj.y, fsj.z, fsj.w );
        }
    }
//printf( "debug Etot= %g\n", Etot );exit(0);    
    return Etot;
}

// ======================================
// =========  OPTIMIZE  =================
// ======================================

/**
 * Calculates the regularization force for each degree of freedom (DOF) based on the difference between the current and target values of the fitted non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb).
 * The regularization force tries to minimize the difference between the current (REQ) and the default value (REQ0) of the parameters. Its strength is controlled by the regularization stiffness (Kreg).
 * The resulting forces are stored in the fDOFs array.
 */
void regularization_force(){
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt   = typToREQ [i];
        const Quat4d& REQ  = typeREQs [i];
        const Quat4d& REQ0 = typeREQs0[i];
        const Quat4d& K    = typeKreg [i];
        if(tt.x>=0)fDOFs[tt.x] = (REQ0.x-REQ.x)*K.x;
        if(tt.y>=0)fDOFs[tt.y] = (REQ0.y-REQ.y)*K.y;
        if(tt.z>=0)fDOFs[tt.z] = (REQ0.z-REQ.z)*K.z;
        if(tt.w>=0)fDOFs[tt.w] = (REQ0.w-REQ.w)*K.w;
    }
}

void regularization_force_walls(){
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt   = typToREQ [i];
        const Quat4d& REQ  = typeREQs [i];
        const Quat4d& REQ0l = typeREQs0_low [i];
        const Quat4d& Kl    = typeKreg_low  [i];
        const Quat4d& REQ0h = typeREQs0_high[i];
        const Quat4d& Kh    = typeKreg_high [i];
        if(tt.x>=0){
            if(REQ.x<REQ0l.x){fDOFs[tt.x] -= Kl.x*(REQ0l.x-REQ.x);}
            if(REQ.x>REQ0h.x){fDOFs[tt.x] -= Kh.x*(REQ0h.x-REQ.x);}
        }
        if(tt.y>=0){
            if(REQ.y<REQ0l.y){fDOFs[tt.y] -= Kl.y*(REQ0l.y-REQ.y);}
            if(REQ.y>REQ0h.y){fDOFs[tt.y] -= Kh.y*(REQ0h.y-REQ.y);}
        }
        if(tt.z>=0){
            if(REQ.z<REQ0l.z){fDOFs[tt.z] -= Kl.z*(REQ0l.z-REQ.z);}
            if(REQ.z>REQ0h.z){fDOFs[tt.z] -= Kh.z*(REQ0h.z-REQ.z);}
        }
        if(tt.w>=0){
            if(REQ.w<REQ0l.w){fDOFs[tt.w] -= Kl.w*(REQ0l.w-REQ.w);}
            if(REQ.w>REQ0h.w){fDOFs[tt.w] -= Kh.w*(REQ0h.w-REQ.w);}
        }
    }
}

/**
 * Limits the fitted non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) to be btween minimum and maximum (REQmin and REQmax).
 * 
 * @param none
 * @return void
 */
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

double run( int nstep, double Fmax, double dt, int imodel_, int isampmode, int ialg, bool bRegularize, bool bClamp, double max_step, bool bEpairs_ ){
    imodel=imodel_;
    bEpairs=bEpairs_;
    double Err=0;
    if( verbosity>1){ printf( "run( nstep %i Fmax %g dt %g isamp %i )\n", nstep, Fmax, dt, isampmode  ); }
    double F2max=Fmax*Fmax;
    double F2;
    for(int i=0; i<nstep; i++){
        //printf("[%i]  DOFs=", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
        DOFsToTypes(); 
        clean_derivs();
        //printf("[%i]  DOFs=", i);for(int j=0;j<W.nDOFs;j++){ printf("%g ",W. DOFs[j]); };printf("\n");
        switch(isampmode){
            case 0: Err = evalDerivsRigid(); break;
            case 1: Err = evalDerivs     (); break;
            case 2: Err = evalDerivsSamp (); break;
        }   
        if( verbosity>1){
            printf("step= %i DOFs= ", i);for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
            //if(bRegularize){ W.regularization_force(); }
            printf("step= %i fDOFs= ", i);for(int j=0;j<nDOFs;j++){ printf("%g ",fDOFs[j]); };printf("\n");
        }
        if(bRegularize){ regularization_force_walls(); }
        if( verbosity>1){ printf("step= %i after_reg fDOFs= ", i);for(int j=0;j<nDOFs;j++){ printf("%g ",fDOFs[j]); };printf("\n"); }
//exit(0);        
        switch(ialg){
            case 0: F2 = move_GD( dt, max_step ); break;
            case 1: F2 = move_MD( dt, max_step ); break;
            case 2: F2 = move_GD_BB_short( i, dt, max_step ); break;
            case 3: F2 = move_GD_BB_long( i, dt, max_step ); break;
            case 4: F2 = move_MD_nodamp( dt, max_step ); break;
        }
        // regularization must be done before evaluation of derivatives
        if(bClamp     ){ limit_params();         }
        //printf("step= %i dt= %g\n", i, dt );
        if( verbosity>1){
            printf("step= %i RMSE= %g |F|= %g\n", i, sqrt(Err), sqrt(abs(F2)) );
            printf("[%i]\n", i );
        }
        if( F2<0.0   ){ printf("DYNAMICS STOPPED after %i iterations \n", i); printf("VERY FINAL DOFs= ");for(int j=0;j<nDOFs;j++){ printf("%.15g ",DOFs[j]); };printf("\n"); return Err; }
        if( F2<F2max ){ printf("CONVERGED in %i iterations \n", i);           printf("VERY FINAL DOFs= ");for(int j=0;j<nDOFs;j++){ printf("%.15g ",DOFs[j]); };printf("\n"); return Err; }
    }
    printf("step= %i DOFs= ", nstep); for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
    printf("VERY FINAL DOFs= ");      for(int j=0;j<nDOFs;j++){ printf("%.15g ",DOFs[j]); };printf("\n");
    return Err;
}





// ------------------------------------------------------------------------------------------------
// --------------------------- NOT USED (ANYMORE OR YET) ------------------------------------------
// ------------------------------------------------------------------------------------------------

/**
 * @brief Cleans the derivative array by setting all its elements to zero.
 * 
 */
void clean_derivs(){ for(int i=0; i<nDOFs; i++){fDOFs[i]=0;} }

/**
 * @brief Renormalizes the weights of a batch of data points. The sum of weights is set to R.
 * 
 * @param R The target sum of weights after renormalization.
 */
void renormWeights(double R=1.0){ 
    double s=0; 
    for(int i=0; i<nbatch; i++){ s+=weights[i]; }
    //printf( "renormWeights() s=%g R=%g R/s=%g \n", s,R, R/s );
    s=R/s; 
    for(int i=0; i<nbatch; i++){ weights[i]*=s; }
    //for(int i=0; i<nbatch; i++){ printf("W.weights[%i]=%g\n",i,weights[i]); }       
}

/**
 * This function reallocates the memory for the fs array if the number of atoms in any batch is greater than the current maximum number of atoms.
 * It first determines the maximum number of atoms in any batch and then checks if it is greater than the current maximum number of atoms.
 * If it is, then it reallocates the memory for the fs array to accommodate the new maximum number of atoms.
 */
void tryRealocTemp(){
    int n=nmax;
    for(int i=0; i<nbatch; i++){ int ni = batch[i].natoms; if(ni>n){ n=ni;} }
    if(n>nmax){ _realloc( fs, n );  nmax=n; };
}

/**
 * This function checks if the current number of atoms in the system exceeds the maximum number of atoms allowed. If so, it reallocates memory for the force array to accommodate the new number of atoms.
 */
void tryRealocTemp_rigid(){
    if(nmax<systemTest0->natoms){
        nmax=systemTest0->natoms;
        _realloc( fs, nmax );
    }
}

/**
 * Initializes the types of the FitREQ object.
 * @param ntype_ The number of types.
 * @param typeMask An array of Quat4i indicating which of the 4 parameters (Rvdw,Evdw,Q,Hb) are free to be fitted.
 * @param tREQs An array of Quat4d objects representing the non-colvalent interaction parameters (Rvdw,Evdw,Q,Hb) for each type.
 * @param bCopy A boolean indicating whether to copy the input arrays to newly allocated internal arrays.
 * @return The number of degrees of freedom.
 */
int init_types( int ntype_, Quat4i* typeMask, Quat4d* tREQs=0, bool bCopy=false ){
    
    printf( "FitREQ::init_types() ntype_=%i ntype=%i bCopy=%i tREQs=%li \n", ntype_, ntype, bCopy, (long)tREQs );
    int nDOFs=0;
    typToREQ = new Quat4i[ntype_];
    for(int i=0; i<ntype_; i++){
        const Quat4i& tm=typeMask[i];
        Quat4i&       tt=typToREQ[i];
        for(int j=0; j<4; j++){
            if(tm.array[j]){
                tt.array[j]=nDOFs;
                nDOFs++;
            }else{
                tt.array[j]=-1;
            }
        }
        printf(  "init_types() [%i] typeMask(%i,%i,%i,%i) typToREQ(%i,%i,%i,%i)  nDOF %i  tREQs(%g,%g,%g,%g) \n", i, tm.x,tm.y,tm.z,tm.w, tt.x,tt.y,tt.z,tt.w, nDOFs, tREQs[i].x,tREQs[i].y,tREQs[i].z,tREQs[i].w );
    }
    realloc(nDOFs);
    //poses = new Mat3d[nbatch];
    //if( bRigid){
    //    poses = new Mat3d[nbatch];
    //}
    if(tREQs){
        ntype = ntype_; 
        if(bCopy){
            typeREQs    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs[i]    = tREQs[i]; }
            typeREQsMin = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Quat4dmin; }
            typeREQsMax = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Quat4dmax; }
            typeREQs0   = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0[i]   = tREQs[i]; }
            typeKreg    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg[i]    = Quat4dZero; }
        }else{ 
            typeREQs = tREQs; 
        }
        DOFsFromTypes(); 
    }
    //for(int i=0; i<ntype; i++){ printf( "init_types()[%i] typeREQs(%g,%g,%g,%g) \n", i,  ); }
    //printf(" DOFs=");for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
    return nDOFs;
}

/**
 * Initializes the types and parameters for fitting the REQ (Quadrupole) values.
 * This function calculates the number of degrees of freedom (nDOFs) and initializes the typToREQ array, which maps each atom type to its corresponding REQ values.
 * The params->atypes array is used to initialize the typeREQs array. t->sym controls which of the 4 parameters (Rvdw,Evdw,Q,Hb) are free to be fitted.
 * It also initializes other arrays related to typeREQs.
 *
 * @return The number of degrees of freedom (nDOFs).
 */
int init_types_par(){
    printf( "FitREQ::init_types_par()\n" );
    int ntype_ = params->atypes.size();
    int nDOFs=0;
    typToREQ = new Quat4i[ntype_];
    for(int i=0; i<params->atypes.size(); i++){
        AtomType& t  = params->atypes[i];        
        Quat4i&   tt = typToREQ[i];
        int bitmask = 1;
        for(int j=0; j<4; j++){
            if( bitmask & t.sym ){
                tt.array[j]=nDOFs;
                nDOFs++;
            }else{
                tt.array[j]=-1;
            }
            bitmask<<=1;
        }
        //printf(  "init_types() [%i] typeMask(%i,%i,%i) typToREQ(%i,%i,%i)  nDOF %i  tREQs(%g,%g,%g,%g) \n", i, tm.x,tm.y,tm.z, tt.x,tt.y,tt.z, nDOFs, tREQs[i].x,tREQs[i].y,tREQs[i].z,tREQs[i].w );
    }
    realloc(nDOFs);
    ntype = ntype_; 
    typeREQs    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs[i]    = Quat4d{ params->atypes[i].RvdW, params->atypes[i].EvdW, params->atypes[i].Qbase, params->atypes[i].Hb }; }
    typeREQs0   = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0[i]   = typeREQs[i]; }
    typeREQsMin = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Quat4dmin;   }
    typeREQsMax = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Quat4dmax;   }
    typeKreg    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg[i]    = Quat4dZero;  }
    DOFsFromTypes(); 
    //for(int i=0; i<ntype; i++){ printf( "init_types()[%i] typeREQs(%g,%g,%g,%g) \n", i,  ); }
    //printf(" DOFs=");for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
    return nDOFs;
}

/**
 * Sets one instace of atomic system
 * @param isys The index of the training system to set. If isys is negative, we set systemTest(-3), systemTest0(-2), or system0(-1).
 * @param na The number of atoms in the system.
 * @param types An array of integers representing the atom types.
 * @param ps An array of Vec3d representing the positions of the atoms.
 * @param bCopy A boolean indicating whether to copy the input arrays to newly allocated internal arrays. default=false
 */
void setSystem( int isys, int na, int* types, Vec3d* ps, bool bCopy=false ){
    
    printf( "FitREQ::setSystem(%i) \n", isys );
    Atoms** pS;
    Atoms*  S;
    if(isys<0){
        if     (isys==-3){ pS=&systemTest;  n1=na; }
        else if(isys==-2){ pS=&systemTest0; n1=na; }
        else if(isys==-1){ pS=&system0;     n0=na; }
        if(*pS==0   ){ *pS = new Atoms(na); }
        S=*pS;
    }else             { S=&batch[isys]; }
    if(bCopy){
        if(S->natoms!=na){ printf("ERORRO FitREQ.setSystem() na(%i)!=S.n(%i) \n", na, S->natoms ); exit(0); }
        for(int i=0; i<na; i++){ S->atypes[i]=types[i]; S->apos[i]=ps[i]; }
    }else{
        S->natoms =na;
        S->atypes =types;
        S->apos   =ps;
    }
}

/**
 * Sets the rigid samples (i.e. poses=translation,rotation) of the atomic system to be fitted.
 * @param n The number of samples.
 * @param Es_     array of doubles containing the energies of each sample.
 * @param poses_  array of Mat3d   containing the rotation and translation of each sample (a=translation, b,c=rotation)
 * @param bCopy  If true, the function will allocate memory and copy the input arrays. If false, the function will only store the pointers to the input arrays.
 * @param bAlloc If true, the function will allocate memory for the arrays. If false, the function will assume that the arrays have already been allocated.
 */
void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy=false, bool bAlloc=false ){

    printf( "FitREQ::setRigidSamples() \n" );
    nbatch = n; 
    if(bCopy){
        if(Es_   ){Es    = new double[nbatch]; }
        if(poses_){poses = new Mat3d [nbatch]; }
        if(Es_   )for(int i=0;i<nbatch; i++){ Es   [i]=Es_   [i]; }
        if(poses_)for(int i=0;i<nbatch; i++){ poses[i]=poses_[i]; }
    }else{
        if(Es_   )Es = Es_; 
        if(poses_)poses = poses_;
    }
}

/*
bool checkToClose(double Rfac){
    int    nj   =system0->natoms;
    int*   jtyp =system0->atypes;
    Vec3d* jpos =system0->apos;
    for(int i=0; i<n; i++){
        int   ti          = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        Qtot+=REQi.z;
        for(int j=0; j<nj; j++){
            int tj              = jtyp[j];
            const Quat4d& REQj  = typeREQs[tj];
            Vec3d d             = jpos[j] - pi;
            double R  = REQi.x+REQj.x;

            double r2 = d.norm2();
            if( r2 < (R*R*Rfac2) ){ return true; }
           
        }
        fs[i].add(fsi);
    }
    return false
}
*/

/**
 * Calculates the energy components of a Lenard-Jones potential, electrostatic interactions and hydrogen bond correction H1 between two systems of atoms
 * @param ni The number of particles in the input set of atoms.
 * @param types An array of integers representing the types of atoms in the input system.
 * @param ps An array of Vec3d representing the positions of atoms in the input system.
 * @param isamp An integer representing the index of current system (sample).
 * @return A Quat4d representing the sum of non-covalent interaction energy components (Pauli,London, Coulomb, Hbond) of the input system with the system0.
 */
Quat4d evalExampleEnergyComponents_LJQH2(int ni, int* types, Vec3d* ps, int isamp,  int nj=0, int* jtyp=0, Vec3d* jpos=0   ){
    
    if( jpos==0 ){   
        nj   =system0->natoms;
        jtyp =system0->atypes;
        jpos =system0->apos;
    }
    Quat4d Esum = Quat4dZero;

    if(Eilines)for(int i=0; i<nj; i++){ Quat4d* Eiline=Eilines[i]; if(Eiline){ Eiline[isamp]=Quat4dZero; } }
    if(Ejlines)for(int j=0; j<nj; j++){ Quat4d* Ejline=Ejlines[j]; if(Ejline){ Ejline[isamp]=Quat4dZero; } }

    for(int i=0; i<ni; i++){
        const int     ti   = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        Quat4d* Eiline=0; if(Eilines){ Eiline=Eilines[i];  }
        for(int j=0; j<nj; j++){
            const int tj        = jtyp     [j];
            const Quat4d& REQj  = typeREQs [tj];
            //const Quat4d  REQij = _mixREQ(REQi,REQj);
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.z*REQj.z;
            double H  = REQi.w*REQj.w;

            const Vec3d d       = jpos     [j] - pi;

            if(H>0) H=0;
            //printf( "ij[%i=%i,%i=%i] REQij(%6.3f,%10.7f,%6.3f,%6.3f) test:REQi(%6.3f,%10.7f,%6.3f,%6.3f) sys0:REQj(%6.3f,%10.7f,%6.3f,%6.3f)\n", i,ti, j,tj,  R,E0,Q,H,   REQi.x,REQi.y,REQi.z,REQi.w,   REQj.x,REQj.y,REQj.z,REQj.w );
            // --- Eectrostatic
            //double ir2     = 1/( d.norm2() + 1e-4 );
            double ir2     = 1/( d.norm2() );
            // --- Lenard-Jones
            double u2   = ir2*(R*R);
            double u6   = u2*u2*u2;
            double u12  = u6*u6;
            Quat4d Eij = Quat4d{
                E0   *u12,                   // Pauli
                E0*-2*u6 ,                   // London
                Q *COULOMB_CONST*sqrt(ir2),  // Coulomb
                H*u12,                       // Hbond 
            };
            Esum.add(Eij);
            if( bDecomp ){
                {                                       if(Eiline){ Eiline[isamp].add(Eij); } }
                if(Ejlines){ Quat4d* Ejline=Ejlines[j]; if(Ejline){ Ejline[isamp].add(Eij); } }
                if(Elines){
                    int ij = i + j*ni;
                    Quat4d* Eline = Elines[ij];
                    if( Eline ){    Eline[isamp]=Eij;  }
                }
            }
        }
    }
    return Esum;
}

Quat4d evalExampleEnergyComponents_BuckinghamQH2(int ni, int* types, Vec3d* ps, int isamp,  int nj=0, int* jtyp=0, Vec3d* jpos=0   ){
    
    if( jpos==0 ){   
        nj   =system0->natoms;
        jtyp =system0->atypes;
        jpos =system0->apos;
    }
    Quat4d Esum = Quat4dZero;

    if(Eilines)for(int i=0; i<nj; i++){ Quat4d* Eiline=Eilines[i]; if(Eiline){ Eiline[isamp]=Quat4dZero; } }
    if(Ejlines)for(int j=0; j<nj; j++){ Quat4d* Ejline=Ejlines[j]; if(Ejline){ Ejline[isamp]=Quat4dZero; } }

    for(int i=0; i<ni; i++){
        const int     ti   = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        Quat4d* Eiline=0; if(Eilines){ Eiline=Eilines[i];  }
        for(int j=0; j<nj; j++){
            const int tj        = jtyp     [j];
            const Quat4d& REQj  = typeREQs [tj];
            //const Quat4d  REQij = _mixREQ(REQi,REQj);
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.z*REQj.z;
            double H  = REQi.w*REQj.w;

            const Vec3d d       = jpos     [j] - pi;

            if(H>0) H=0;
            double r = sqrt( d.norm2() );
            //double B  = exp  ( R *  kMorse );
            //double EB = B*exp( r * -kMorse );
            double e = exp( (R-r) * kMorse );

            double ir      = 1/r;
            double u       = ir*R;
            // --- Buckingam
            double u2   = u*u;
            double u6   = u2*u2*u2;
            Quat4d Eij = Quat4d{
                E0*e,                 // Pauli
                E0*-2*u6 ,            // London
                Q *COULOMB_CONST*ir,  // Coulomb
                H*e,                  // Hbond 
            };
            Esum.add(Eij);
            if( bDecomp ){
                {                                       if(Eiline){ Eiline[isamp].add(Eij); } }
                if(Ejlines){ Quat4d* Ejline=Ejlines[j]; if(Ejline){ Ejline[isamp].add(Eij); } }
                if(Elines){
                    int ij = i + j*ni;
                    Quat4d* Eline = Elines[ij];
                    if( Eline ){    Eline[isamp]=Eij;  }
                }
            }
        }
    }
    return Esum;
}

/**
 * Calculates the derivative of the energy with respect to the total charge in order to force the system toward the neutral charge.
 * @param n The number of atoms in the system.
 * @param types An array of integers representing the types of atoms in the system.
 * @param ps An array of Vec3d representing the positions of atoms in the system.
 * @param Qtot The total charge of the system.
 * @return Energy penalty due to the deviation of the total charge from the neutral charge.
 */
double evalExampleDerivs_Qneutral(int n, int* types, Vec3d* ps, double Qtot ){
    
    // Qtot = Sum_i * Q_i
    // E = K*Qtot^2
    // d E^2 / dqi = K*Qtot*(  )
    double E = 0;
    double fQ = -Kneutral*Qtot;
    for(int i=0; i<n; i++){
        //int   ti          = types[i];
        //const Quat4d& REQi = typeREQs[ti]; 
        fs[i].z += fQ;
        E       += Kneutral*Qtot*Qtot;
    }
    return  E;
}

/**
 * @brief Evaluates the variational derivatives of the fitting error (sumed batch training samples) with respect to all fitting parameters (i.e. non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) for each atom type). The atomic systems are not assumed rigid (i.e. the atomic positions can vary freely from sample to sample).
 * 
 * @param Eout array to store the non-covalent interaction energy values of each atomic system in the batch. if Eout==null, the function will not store the energy values.
 * @return double, returns the total fitting error.
 */
double evalDerivs( double* Eout=0 ){ 
    printf( "FitREQ::evalDerivs() nbatch %i imodel %i verbosity %i \n", nbatch, imodel, verbosity );
    tryRealocTemp();
    double Error = 0;
    for(int i=0; i<nbatch; i++){
        //printf("evalDerivs[%i]\n", i );
        const Atoms& C = batch[i];
        //C.print();
        double Eref=Es[i];
        double wi = 1; 
        if(weights) wi = weights[i];
        // switch(ikind){ case ikind_LJQ:
        clean_fs(C.natoms);
        double E=0;
        switch (imodel){
            // TBD
            //case 1: E = evalExampleDerivs_LJQH         (C.natoms, C.atypes, C.apos ); break;
            //case 2: E = evalExampleDerivs_LJQH2        (C.natoms, C.atypes, C.apos ); break;
            //case 3: E = evalExampleDerivs_MorseQH2     (C.natoms, C.atypes, C.apos ); break;
            //case 4: E = evalExampleDerivs_BuckinghamQH2(C.natoms, C.atypes, C.apos ); break;
        }
        if( E>1e+300 ){
            if(verbosity>0) printf( "skipped sample [%i] atoms too close \n", i );
            continue;
        } 
        if(Eout){ Eout[i]=E; };
        double dE = (E - Eref)*wi;
        Error += dE;
        acumDerivs(C.natoms, C.atypes, dE );
        //double dE = C.E - E;
        //double E_ = evalExampleDerivs_LJQ(C.n, C.atypes, C.apos, dE*wi );  // Backward pass
        //exit(0);
    }
    return Error;
}

/**
 * @brief Evaluates the variational derivatives of the fitting error (sumed batch training samples) with respect to all fitting parameters (i.e. non-covalent interaction parameters REQH(Rvdw,Evdw,Q,Hb) for each atom type). The atomic systems are assumed rigid (i.e. the system is translated and rotated rigidly from sample to sample).
 * 
 * @param Eout array to store the non-covalent interaction energy values of each atomic system in the batch. if Eout==null, the function will not store the energy values.
 * @return double, eturns the total fitting error.
 */
double evalDerivsRigid( double* Eout=0 ){
    printf( "FitREQ::evalDerivsRigid() nbatch %i imodel %i verbosity %i \n", nbatch, imodel, verbosity );
    tryRealocTemp_rigid();
    const Atoms& C0 = *systemTest0;
    const Atoms& C  = *systemTest;
    double Error = 0; 
    for(int i=0; i<nbatch; i++){
        double Eref=Es[i]; 
        rigid_transform( poses[i].a, 0, poses[i].c, poses[i].b, C.natoms, C0.apos, C.apos );
        //savexyz("FitREQ_debug.xyz", system0, systemTest,0, "a" );
        clean_fs(C.natoms);
        double E =0;
        switch (imodel){
            // TBD
            //case 1: E = evalExampleDerivs_LJQH (C.natoms, C.atypes, C.apos ); break;
            //case 2: E = evalExampleDerivs_LJQH2(C.natoms, C.atypes, C.apos ); break;
            //case 3: E = evalExampleDerivs_MorseQH2     (C.natoms, C.atypes, C.apos ); break;
            //case 4: E = evalExampleDerivs_BuckinghamQH2(C.natoms, C.atypes, C.apos ); break;
        }
        if( E>1e+300 ){

        } 
        if(Eout){ Eout[i]=E; };
        double wi = 1;
        if(weights) wi = weights[i];
        //if(bLimitForce){ };
        double dE = -(E - Eref)*wi;
        Error += dE*dE;
        acumDerivs(C.natoms, C.atypes, dE );
        //printf( "[%i] x %g E %g Eref %g ", i, poses[i].a.x, E, Eref );printf("}\n");
        //printf("fs={");for(int j=0;j<C.natoms;j++){ printf("(%g,%g,%g)",fs[j].x,fs[j].y,fs[j].z); };//printf("}\n");
        //printf("fDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",fDOFs[j]); };printf("}\n");
    }
    //printf("fDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",fDOFs[j]); };printf("}\n");
    //printf("vDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",vDOFs[j]); };printf("}\n");
    return Error;
}

/**
 * @brief Loads atomic coordinates and types from an XYZ file and split them into two sets of atoms (rigid fragment system0 and moving fragment stored to training batch). The selection of atoms for system0 is specified by the indices in the i0s array. The selection of atoms for the training batch is specified by the indices in the itests array.
 * 
 * @param fname The name of the XYZ file to load.
 * @param n0 The number of atoms in the system0 vector.
 * @param i0s An array of indices indicating which atoms in the XYZ file correspond to the atoms in the system0 vector (i.e. the rigid fragment).
 * @param ntest The number of atoms in the batch vector.
 * @param itests An array of indices indicating which atoms in the XYZ file correspond to the atoms in the traning batch vector (i.e. the moving fragment).
 * @param types0 An optional array of atom types for the atoms in the system0 vector. if types0==0, the function will read the atom types from the XYZ file.
 * @param testtypes An optional array of atom types for the atoms in the batch vector. if testtypes==0, the function will read the atom types from the XYZ file.
 * @return The number of batches loaded.
 */
int loadXYZ( const char* fname, int n0, int* i0s, int ntest, int* itests, int* types0=0, int* testtypes=0 ){
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    int isys = 0;
    int il   = 0;
    int na   = 0; 
    Atoms atoms;
    if(!system0){ system0=new Atoms(n0); }else { system0->realloc( n0); }
    Atoms S         ( ntest );
    bool bReadTypes = !(types0 && testtypes);
    printf( "FitREQ::loadXYZ() n0 %i ntest %i \n", n0, ntest );
    while( fgets(line, nline, fin) ){
        //printf( ">>%s<<\n", line );
        //printf( "il %i isys %i na %i \n", il, isys, na, ntest,  );
        if      ( il==0 ){
            sscanf( line, "%i", &na );
            if( (na>0)&&(na<10000) ){
                if(isys==0){
                    atoms.allocNew(na);
                }
                S.allocNew(ntest);
            }else{ printf( "ERROR in FitREQ::loadXYZ() Suspicious number of atoms (%i) while reading `%s`  => Exit() \n", na, fname ); exit(0); }
        }else if( il==1 ){
            // comment - ToDo : Here we can read the reference Energy directly from .xyz 
        }else if( il<na+2 ){
            double x,y,z,q;
            //printf( ">>%s<<\n", line );
            int nret = sscanf( line, "%s %lf %lf %lf %lf", at_name, &x, &y, &z, &q );
            if(nret<5){q=0;}
            int i=il-2;
            //printf( "[%i] `%s` (%g,%g,%g) q %g \n", i, at_name, x, y, z, q );
            atoms.apos  [i].set(x,y,z);
            if(bReadTypes){ atoms.atypes[i]=params->getAtomType(at_name); }else{ atoms.atypes[i]=-1; }
        }
        il++;
        if( il >= na+2 ){ 
            //printf( "===== FitREQ::loadXYZ() isys %i na %i n0 %i ntest %i \n", isys, na, system0->natoms, S.natoms  );
            //printf( "atoms " ); atoms.print();
            if(isys==0)for(int i=0; i<n0;    i++ ){ int ia=i0s   [i];  system0->apos[i]=atoms.apos[ia]; int t; if(types0   ){ t=types0   [i]; }else{ t=atoms.atypes[ia]; }; system0->atypes[i]=t; } // store to system0
            {          for(int i=0; i<ntest; i++ ){ int ia=itests[i];         S.apos[i]=atoms.apos[ia]; int t; if(testtypes){ t=testtypes[i]; }else{ t=atoms.atypes[ia]; };  S      .atypes[i]=t; } // store to batch[isys]
            //printf( "system0 " ); system0->print();
            //printf( "stests  " ); S       .print();
            batch_vec.push_back( S ); }
            il=0; isys++; 
        }
    }
    nbatch =  batch_vec.size();
    batch  = &batch_vec[0];
    _realloc(Es,nbatch);
    atoms.dealloc();
    fclose(fin);
    return nbatch;
}

// void atomsToXYZ(FILE* fout, int n, int* types, Vec3d* ps){
//     for(int i=0; i<n; i++){
//         fprintf( fout, "%i %20.10f %20.10f %20.10f\n", types[i], ps[i].x,ps[i].y,ps[i].z );
//     }
// }

}; // class FitREQ


#endif