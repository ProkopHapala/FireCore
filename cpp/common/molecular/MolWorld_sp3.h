
#ifndef MolWorld_sp3_h
#define MolWorld_sp3_h
/// @file MolWorld_sp3.h @brief contains MolWorld_sp3 class, which is a comprehensive class storing the state of a molecular simulation including bonding,non-bodning of molecules and molecules with substrate
/// @ingroup Classical_Molecular_Mechanics

#include <stdlib.h>
#include <stdio.h>
#include <sys/stat.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <omp.h>

#include "IO_utils.h"

//#include "testUtils.h"
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"

#include "MMFFparams.h"
static MMFFparams* params_glob;

//#include "raytrace.h"
#include "constants.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "MMFFsp3_loc.h"
#include "MMFFf4.h"
#include "UFF.h"

#include "NBFF.h"
#include "GridFF.h"
#include "RigidBodyFF.h"
#include "QEq.h"
#include "constrains.h"
#include "molecular_utils.h"

#include "LimitedGraph.h"

#include "Molecule.h"
#include "MMFFBuilder.h"
#include "SMILESparser.h"
#include "DynamicOpt.h"

#include "MultiSolverInterface.h"
#include "GlobalOptimizer.h"
#include "GOpt.h"
#include "Groups.h"

#include "datatypes_utils.h"

#include "arrayAlgs.h"
#include "SVG_render.h"

#include "EwaldGrid.h"



enum class MolWorldVersion{ BASIC=0, QM=1, GPU=2 };
//static inline MolWorldVersion operator|(MolWorldVersion a, MolWorldVersion b) { return (MolWorldVersion)( ((int)a) | ((int)b)  );  };
//static inline MolWorldVersion operator&(MolWorldVersion a, MolWorldVersion b) { return (MolWorldVersion)( ((int)a) & ((int)b)  );  };
//static inline int operator|(MolWorldVersion a, MolWorldVersion b) { return ( ((int)a) | ((int)b)  );  };
//static inline int operator&(MolWorldVersion a, MolWorldVersion b) { return ( ((int)a) & ((int)b)  );  };
//bool hasFlag(MolWorldVersion value, MolWorldVersion flag) { return (bool)( ((int)value) & ((int)flag)  );  }

/// @brief Comprehensive class storing the state of a molecular simulation including bonding,non-bodning of molecules and molecules with substrate
/// @details It impolements SolverInterface, various methods of initialization and various similation loops (serial or OpenMP parallelized)   
class MolWorld_sp3 : public SolverInterface { public:
    bool isInitialized=false;
    //const char* data_dir     = "common_resources";
    const char* data_dir     = 0;
    const char* xyz_name     = "input";
    //const char* surf_name    = "surf";
    const char* surf_name       = 0;
    const char* substitute_name = 0;    int isubs;
    //const char* lvs_name     ="input.lvs";
    //const char* surflvs_name ="surf.lvs";
    const char* smile_name   = 0;
    const char* constr_name  = 0;
    Vec3i nMulPBC  = Vec3iZero; 

    //const char* trj_fname    = "trj.xyz";
    const char* trj_fname    = 0;
    int savePerNsteps = 1;
    OptLog opt_log;
    Vec3i nPBC_save{1,1,1};

    // ---  Parameters for Molecular Dynamics at non-zero temperature   -- Moved to GOpt go;
    //bool bThermalSampling = false;   // if >0 then we do thermal sampling
    //double T_target         = 300.0;   // target temperature for thermal sampling (if bThermalSampling>0)
    //double T_current        = 0.0;   // current temperature for thermal sampling (if bThermalSampling>0)
    //double gamma_damp       = 0.01;  // damping factor for thermal sampling (if bThermalSampling>0)

    double fAutoCharges=-1;
    bool bEpairs = false;
    bool bToCOG=false;
    bool bCellBySurf=false;
    int   bySurf_ia0 =0;
    Vec2d bySurf_c0=Vec2dZero;
    Vec2d bySurf_lat[2];

    Mat3d new_lvec=Mat3dIdentity;

    Mat3d debug_rot; // This was used for debuging molecular orientation

	// Building
	MMFFparams    params;
    //MMFFparams*   params = null;
	MM::Builder  builder;
	SMILESparser smiles;

    Vec3d* apos_bak=0;

	// Force-Fields & Dynamics
	MMFFsp3      ff;
    MMFFsp3_loc  ffl;
    //MMFFf4     ff4;
    UFF          ffu;
	Constrains   constrs;
	//NBFF_old   nff;
    NBFF         surf, nbmol;
	GridFF       gridFF;
    bool bGridDouble = true;

    EwaldGrid gewald;

    std::vector<int> atom2group;
    Groups groups;
    Mat3d bbox{   -1e+8,-1e+8,-1e+8,  +1e+8,+1e+8,+1e+8,    -1,-1,-1 };



    RigidBodyFF  rbff;
    QEq          qeq;
	DynamicOpt   opt;
    DynamicOpt   optRB;  // rigid body optimizer

    LimitedGraph<N_NEIGH_MAX> graph; // used to find bridges in the molecule, and other topology-related algorithms

    GlobalOptimizer gopt;
    GOpt            go;
    bool   bGopt =false;
    int gopt_ifound=0;
    int gopt_nfoundMax=1000000;

    bool bCheckStuck=false;
    double RStuck=0.2;
    int nStuckMax=100;
    int nStuckTrj=10;
    int nStuck=0;

    GridShape MOgrid;

    SVG_render svg;

    int  icurIter = 0;
    int  iterPerFrame=50;
    int  iParalel   =  -101; 
    int  iParalelMax=1;
    int  iParalelMin=0;
    int iParalel_default=1;
    bool bOcl=false; // used only in Ocl version

    double gridStep = 0.1; 
    //double gridStep = 0.05; 
    //double gridStep = 0.2; 
    //Vec3i nPBC{0,0,0};   // just debug
    Vec3i nPBC{1,1,0};
    //Vec3i nPBC{1,3,0};
    int    npbc       = 0;
    Vec3d* pbc_shifts = 0;
    int    ipbc0=0;

	// state
	bool bConverged = false;
	double  Etot=0;
	double  maxVcog = 1e-9;
	double  maxFcog = 1e-9;
	double  maxTg   = 1e-1;
	double  Kmorse = -1.0;
    double  Ftol_default = 1e-4;
    double  dt_default   = 0.05;

    double time_per_iter = 0;

    int  iAtomNotConv = -1;
    FILE* atomTrjFile=0;

	// force-eval-switchefs
    int  imethod=0;

	bool doBonded          = false; // 1
	bool bNonBonded        = true;  // 2
    bool bGroups           = false; // 3
    bool bConstrains       = false; // 4
	bool bSurfAtoms        = false; // 5
    bool bGridFF           = false; // 6
    bool bTricubic         = false; // 7
	bool bPlaneSurfForce   = false; // 7
    bool bMMFF             = true;  // 8
    bool bUFF              = false; // 9
    bool b141              = true;  // 10 // seems not to be used in assignUFFtypes()
    bool bSimple           = false; // 11 // use assignUFFtypes_simplerule() or assignUFFtypes_findrings() in assignUFFtypes()
    //bool bSimple           = true;   // use assignUFFtypes_simplerule() or assignUFFtypes_findrings() in assignUFFtypes()
    bool bConj             = true; // 12  // manually change sp3 nitrogen and oxygen to "resonant" when they are bonded to an sp2 atom (conjugation) in assignUFFtypes()
    bool bCumulene         = true; // 13  // exception to avoid cumulenes in assignUFFtypes()
    bool bRigid            = false; // 14
	bool bOptimizer        = true;  // 15
	bool bPBC              = false; // 16
	bool bCheckInvariants  = true;  // 17
    bool bRelaxPi          = false; // 18
    bool bChargeUpdated    = false; // 19
    bool bAnimManipulation = false; // 20
    bool bNonBondNeighs    = false; // 21
    bool bWhichAtomNotConv = false; // 22
    bool bCheckInit        = false; // 23

    Vec3d anim_vec;
    float anim_speed;
	
    // ToDo: pivot can be later replaced by gizmo
    Vec3d pivotPoint = Vec3dZero;
    Mat3d pivotRot   = Mat3dIdentity;
    Vec3d cog,vcog,fcog,tqcog;
    int nloop=0;

	// Selecteion & Manipulasion
	std::vector<int>            selection;
    //std::unordered_map<int,int> selection_map;
    std::unordered_set<int>     selection_set;
	Vec3d manipulation_p0=Vec3dZero; 
	Vec3d manipulation_ax=Vec3dZ;
	int*  manipulation_sel=0;
	int   manipulation_nsel=0;

    std::vector<int> constrain_list;
    double Kfix=1.0;

    int ipicked    = -1; // picket atom 
    int ibpicked   = -1; // picket bond
    int iangPicked = -1; // picket angle
    Vec3d* picked_lvec = 0;
    Vec3d pick_hray, pick_ray0;

    bool   bConstrZ=false;
    double ConstrZ_xmin=0.0;
    double ConstrZ_l=0.0;
    double ConstrZ_k=1.0;

    std::vector<Vec3i> Hbonds;

    Mat3d* dlvec = 0;

    // lattice scan
    bool   bLatScan      = false;
    Mat3d* latscan_dlvec = 0;
    Vec2i  latscan_n{0,0};

	// IO temp & aux
	FILE* xyz_file=0;
	//char* tmpstr;

    double Kpick  = -2.0;
    double QEpair = -0.2;
    //int nEp=0;    // number of electron pairs
    //int etyp=0;   // electron pair type

    int itest = 0;

    int nSystems    = 1;
    int iSystemCur  = 0;    // currently selected system replica

    bool bRelax=false;


    // ========== from python interface

    virtual int getMolWorldVersion() const { return (int)MolWorldVersion::BASIC; };

    virtual int getGroupPose( Quat4f*& gpos, Quat4f*& gfw, Quat4f*& gup ){ gpos=0; gfw=0; gup=0; return 0; };
    virtual void stopExploring (){ go.bExploring=false; };
    virtual void startExploring(){ go.startExploring(); };
    virtual int getMultiConf( float* Fconvs , bool* bExplors ){ return 0; };

    virtual int getTitle(char* s){
        int nstr=0;
        if(xyz_name ){ nstr += sprintf( s     , "%s",   xyz_name  ); }
        if(surf_name){ nstr += sprintf( s+nstr, " @%s", surf_name ); }
        return nstr;
    }

    virtual void init(){
        printf( "MolWorld_sp3::init() verbosity=%i\n", verbosity );
        //params.verbosity=verbosity;
        //printf(  "MolWorld_sp3:init() params.verbosity = %i \n", params.verbosity );
        printf("params.atypes.size() %i\n", params.atypes.size() );
        if( params.atypes.size() == 0 ){
            initParams( "common_resources/ElementTypes.dat", "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat", "common_resources/DihedralTypes.dat" );
        }
        // initialization global optimizer
        gopt.solver = this;
        params_glob = &params;
        if(verbosity>0){
            printf("\n#### MolWorld_sp3::init()\n");
            if(smile_name   )printf("smile_name  (%s)\n", smile_name );
            if(data_dir     )printf("data_dir    (%s)\n", data_dir );
            if(xyz_name     )printf("xyz_name    (%s)\n", xyz_name );
            if(surf_name    )printf("surf_name   (%s)\n", surf_name );
            if(substitute_name)printf("substitute_name  (%s)\n", substitute_name );
            // TBD we should also print if we use UFF or not...
            printf( "MolWorld_sp3::init() bMMFF %i bUFF %i bRigid %i\n", bMMFF, bUFF, bRigid );
        }
        if(surf_name ){
            bGridFF = true;
            //double z0 = 0.0;   // This is how we have it in python API i.e. MMFF.py
            double z0 = NAN;   // This makes inconsistency with python API i.e. MMFF.py
            loadSurf( surf_name, bGridFF, idebug>0, z0 );
        }
        if ( smile_name ){               
            insertSMILES( smile_name );    
            builder.addAllCapTopo();       
            builder.randomizeAtomPos(1.0); 
            bMMFF=true;
        }else if ( xyz_name ){
            if( bMMFF ){ 
                printf("buildMolecule_xyz( %s )\n", xyz_name);
                buildMolecule_xyz( xyz_name );
            }else{
                printf("MolWorld_sp3::init() loading %s\n", xyz_name);
                loadNBmol( xyz_name ); 
                if(bRigid)initRigid();
            }
        }
        builder.randomFragmentCollors();
        if(bMMFF){     
            makeFFs();
            if(bCheckStuck)apos_bak = new Vec3d[ffl.natoms];
        }
        if(!bUFF){ builder.setup_atom_permut( true ); }
        if(constr_name ){ constrs.loadBonds( constr_name, &builder.atom_permut[0], 0 );  }
        if(dlvec       ){ add_to_lvec(*dlvec);    }  // modify lattice after initialization - it helps to build constrained systems 
        //builder.printAtoms();
        //printf( "MolWorld_sp3::init() ffl.neighs=%li ffl.neighCell-%li \n", ffl.neighs, ffl.neighCell );
        //ffl.printNeighs();
        if(verbosity>0) 
        printf( "#### MolWorld_sp3::init() DONE\n\n");
    }

    virtual void pre_loop(){
        groups.print_groups2atoms(); //exit(0);
        // int ngroup = 0;
        // printf("atom2group.size()==%i\n", atom2group.size() );
        // for(int i=0; i<atom2group.size(); i++){ 
        //     //printf("atom2group[%i]==%i\n", i, atom2group[i]);
        //     ngroup=_max(ngroup,atom2group[i]); 
        // } 
        // ngroup++;
        // if( ngroup>0 ){
        //     groups.setGroupMapping( ffl.natoms, ngroup, &atom2group[0] );
        //     groups.bindAtoms(ffl.apos, ffl.fapos);
        //     //groups.initWeights(ffl.natoms);
        //     groups.evalAllPoses();
        // }
    }

    // ========== Render to SVG

    void renderSVG( const char* fname, Vec3d nPBC=Vec3d{0,0,0}, Mat3d rotMat=Mat3dIdentity, bool bAtoms=true, bool bBonds=true, bool bCaps=true, bool bAtomIndex=false, float Rsc=0.25, float Rsub=0.5 ){
    
        svg.rot = rotMat;
        Vec3d cog =  cog_bbox( ffl.natoms, ffl.apos );  svg.cog = cog;

        char str[256];

        int na = ffl.nnode;
        if(bCaps) na = ffl.natoms; 
        svg.findViewport( na, ffl.apos );
        svg.open(fname);

        //bool bOrig = true;

        // order of rendering:    arrange atoms from back to front (using z-order in the current view)
        std::vector<int>   z_order(na);
        std::vector<float> zs     (na);
        for(int ia=0; ia<na; ia++){ 
            Vec3d p; rotMat.dot_to( ffl.apos[ia], p );
            zs     [ia] = p.z; 
            z_order[ia] = ia;
        }
        //sort( z_order.begin(), z_order.end(), [&zs](int i1, int i2){ return zs[i1]<zs[i2]; } ); // sort indexes by z
        quickSort<float>( &zs[0], &z_order[0], 0, na-1 );

        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
            for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    svg.cog = cog + ffl.lvec.a*ix + ffl.lvec.b*iy + ffl.lvec.c*iz;
                    printf( "ix,iy,iz %i %i %i cog(%g,%g,%g) \n", ix,iy,iz, svg.cog.x,svg.cog.y,svg.cog.z );
                    bool bOrig = (ix==0)&&(iy==0)&&(iz==0);

                    // --- bonds
                    if(bBonds){
                        svg.stroke_opacity = 1.0;
                        svg.stroke_width   = 3.0;
                        svg.beginPath();
                        for(int ia=0; ia<na;    ia++){
                            for(int j=0; j<4; j++){
                                int ib = ffl.neighs[ia].array[j];
                                if (ib<0)                         continue; 
                                if ( (!bCaps) && (ia>ffl.nnode) ) continue;
                                svg.path_move( ffl.apos[ia] );
                                Vec3d pj = ffl.apos[ib];
                                if(ffl.shifts) pj.add( ffl.shifts[ffl.neighCell[ia].array[j]]);
                                svg.path_line( pj );
                            }
                        }
                        svg.endPath();
                    }

                    // --- atoms
                    svg.stroke_width   = 1.0;
                    svg.stroke_opacity = 0.5;
                    const char* atom_style = "atom_style";
                    svg.writeCurrentStyle( atom_style );
                    if(bAtoms){
                        for(int i=0; i<na; i++){
                            int ia = z_order[i]; 
                            //int ia = i;
                            //printf( "zorder[i=%i] ia=%i \n", i, ia );
                            int it         = ffl.atypes[ia];
                            svg.color_fill = params_glob->atypes[it].color;
                            float r        = (params_glob->atypes[it].RvdW-Rsub)*Rsc; 
                            svg.drawCircle( ffl.apos[ia], r, atom_style ); 
                        }
                    }
                    
                    if(bOrig & bAtomIndex){
                        const char* label_style = "label_style";
                        svg.beginStyle( label_style );
                        svg.print("font-family: \"DejaVu Mono, monospace\"\n");
                        svg.print("font-size: 20\n");
                        svg.print("fill: #000000\n");
                        svg.endStyle();
                        svg.color_fill = 0x000000;
                        svg.font_size = 15;
                        for(int ia=0; ia<ffl.natoms;    ia++){
                            sprintf(str,"%i",ia);
                            if(bAtomIndex) svg.drawText( str, ffl.apos[ia], label_style );
                        }
                    }
        
                }
            }
        }
        svg.close();
    }


    // ===============================================
    //       Implement    SolverInterface
    // ===============================================

    virtual int getHOMO(){ return 0; };
    virtual int projectOrbital(int iMO, double*& ewfaux ){ ewfaux=0; return 0; };  
    virtual int projectDensity(         double*& ewfaux ){ ewfaux=0; return 0; };



/**
 * Adds a distance constraint between two atoms to the molecular world.
 *
 * @param i0 The index of the first atom.
 * @param i1 The index of the second atom.
 * @param lmin The minimum allowed distance between the atoms (default: 1.0).
 * @param lmax The maximum allowed distance between the atoms (default: 2.0).
 * @param kmin The minimum force constant for the constraint (default: 0.0).
 * @param kmax The maximum force constant for the constraint (default: 1.0).
 * @param flim The maximum force limit for the constraint (default: 10.0).
 * @param shift The shift vector to be applied to the atoms (default: 0).
 * @param bOldIndex Flag indicating whether to use old atom indices (default: false).
 */
void addDistConstrain( int i0,int i1, double lmin=1.0,double lmax=2.0,double kmin=0.0,double kmax=1.0,double flim=10.0, Vec3d shift=Vec3dZero, bool bOldIndex=false ){
    if(bOldIndex){
        i0 = builder.atom_permut[i0];
        i1 = builder.atom_permut[i1];
    }
    constrs.bonds.push_back( DistConstr{ {i0,i1}, {lmax,lmin}, {kmax,kmin}, flim, shift } );
}


/**
 * Sets the constraints for the molecular world.
 * 
 * @param bClear Flag indicating whether to clear existing constraints. Default is true.
 * @param Kfix_ The fixed constraint strength. Default is 1.0.
 */
    virtual void setConstrains(bool bClear=true, double Kfix_=1.0 ){
        double Kfix=Kfix_;
        for(int i=0; i<ffl.natoms; i++){ ffl.constr[i].w=-1; }
        for(int i: constrain_list     ){ 
            //printf( "setConstrains %i \n", i );
            ffl.constr[i].w=Kfix; ffl.constr[i].f=ffl.apos[i]; 
        }
    }

/**
 * @brief Changes the lattice vector of the molecular world.
 * 
 * This function updates the lattice vector of the molecular world
 * and performs necessary calculations and bindings based on the new lattice vector.
 * 
 * @param lvec The new lattice vector to be set.
 */
virtual void change_lvec( const Mat3d& lvec ){
    ffl.setLvec( lvec );
    //npbc = makePBCshifts( nPBC, lvec );
    evalPBCshifts( nPBC, ffl.lvec, pbc_shifts );
    ffl.bindShifts(npbc,pbc_shifts);
    builder.lvec = lvec;
}

/**
 * Adds the given displacement vector to the lattice vector of the molecular world.
 * This function updates the lattice vector, evaluates the periodic boundary conditions (PBC) shifts,
 * and binds the shifts to the molecular world.
 * 
 * @param dlvec The displacement vector to be added to the lattice vector.
 */
virtual void add_to_lvec( const Mat3d& dlvec ){
    ///printf("MolWold_sp3::add_to_lvec()\n");
    //printf(  "BEFORE ffl.lvec " ); printMat(ffl.lvec);
    ffl.setLvec( ffl.lvec+dlvec );
    //npbc = makePBCshifts( nPBC, lvec );
    evalPBCshifts( nPBC, ffl.lvec, pbc_shifts );
    ffl.bindShifts(npbc,pbc_shifts);
    builder.lvec = ffl.lvec;
    //printf(  "AFTER ffl.lvec " ); printMat(ffl.lvec);
}


/**
 * Changes the lattice vector of the molecular world using relaxation method.
 * 
 * @param nstep     The number of relaxation steps to perform.
 * @param nMaxIter  The maximum number of iterations for the relaxation solver.
 * @param tol       The tolerance for convergence of the relaxation solver.
 * @param dlvec     The change in lattice vector to be added at each step.
 */
virtual void change_lvec_relax( int nstep, int nMaxIter, double tol, const Mat3d& dlvec ){
    printf( "MolWorld_sp3::change_lvec_relax() \n" );
    for(int i=0; i<nstep; i++){
        add_to_lvec( dlvec );
        printf( "change_lvec_relax()[%i] lvec(%6.2f,%6.2f,%6.2f) \n", i, ffl.lvec.a.x,ffl.lvec.a.x,ffl.lvec.a.x );
        solve( nMaxIter, tol );
    }
}

/**
 * Solves the molecular system using the specified parameters.
 * 
 * @param nmax The maximum number of iterations.
 * @param tol The tolerance for convergence.
 * @return The total energy of the system.
 */
virtual double solve( int nmax, double tol )override{
    if(nmax>1){
        bRelax=true;
    }else{
        bRelax=false;
    }

    long t0=getCPUticks();
    int nitr = run_omp_Milan( nmax, opt.dt_max, tol, 1000.0, -1. );
    long t=(getCPUticks()-t0); if(verbosity>1)printf( "time run_omp[%i] %g[Mtick] %g[ktick/iter]  %g[s] %g[ms/iter]\n", nitr, t*1e-6, t*1.e-3/nitr, t*tick2second, t*tick2second*1000/nitr  );
    if(std::isnan(Etot)){ printf( "ERROR: Etot is NaN\n" ); }
    return Etot;
}

/**
 * Sets the geometry of the atoms
 * 
 * @param ps Pointer to an array of Vec3d representing the positions of the atoms.
 * @param lvec Pointer to a Mat3d representing the lattice vectors.
 */
virtual void setGeom( Vec3d* ps, Mat3d *lvec )override{
    //printf( "MolWorld::setGeom()\n" );
    //printf("ffl.lvec\n"    ); printMat( ffl.lvec );
    //printf("   *lvec\n"    ); printMat(    *lvec );
    change_lvec( *lvec );
    //printMat( ffl.lvec );
    //printPBCshifts();
    for(int i=0; i<ffl.natoms; i++){
        //printf( "setGeom[%i] ffl.apos(%6.3f,%6.3f,%6.3f) ps(%6.3f,%6.3f,%6.3f) \n", i, ffl.apos[i].x, ffl.apos[i].y, ffl.apos[i].z,   ps[i].x,ps[i].y,ps[i].z );
        ffl.apos[i] = ps[i];
        ffl.vapos[i] = Vec3dZero;
    }
    //ffl.initPi( pbc_shifts );
}

/**
 * @brief Retrieves the geometry of the atoms
 * 
 * @param ps Pointer to an array of Vec3d objects to store the atom positions.
 * @param lvec Pointer to a Mat3d object to store the lattice vectors (optional).
 * @return The total energy of the molecular system.
 */
virtual double getGeom( Vec3d* ps, Mat3d *lvec )override{
    //printf( "MolWorld::getGeom()\n" );
    //printf("getGeom ffl.lvec\n"    ); printMat( ffl.lvec );
    if(lvec){ *lvec=ffl.lvec; }
    //for(int i=0; i<ffl.nvecs; i++){
    for(int i=0; i<ffl.natoms; i++){
        ps[i]=ffl.apos[i];
    }
    //printf( "MolWorld_sp3::getGeom() Etot=%g \n ", Etot );
    return Etot;
}




/**
 * @brief Optimizes the lattice in one dimension.
 * 
 * This function optimizes the lattice in one dimension by performing a lattice scan.
 * It adjusts the lattice vectors and atom positions to minimize the energy of the system.
 * 
 * @param n1 The number of lattice scans in the negative direction.
 * @param n2 The number of lattice scans in the positive direction.
 * @param dlvec The displacement vector used for the lattice scan.
 */
virtual void optimizeLattice_1d( int n1, int n2, Mat3d dlvec ){
    printf("\n\n\n######### MolWorld_sp3::optimizeLattice_1d(%i.%i) \n", n1, n2 );
    //printMat( ffl.lvec );
    //printPBCshifts();
    //ffl.print_apos();
    //printf("ffl.lvec\n"    ); printMat( ffl.lvec    );
    //printf("ffl.invLvec\n" ); printMat( ffl.invLvec );
    //gopt.reallocPop( n1+n2, ffl.nvecs );
    //gopt.atypes = ffl.atypes;

    
    gopt.reallocPop( n1+n2, ffl.natoms, true );

    //gopt.tolerance = 0.02;
    gopt.tolerance = 0.01;

    //ffl.constrainAtom(10);
    
    for(int i=0; i<ffl.natoms; i++ ){ gopt.atypes[i]= params.atypes[ffl.atypes[i]].iZ; }
    //Mat3d lvec0 = builder.lvec;
    Mat3d lvec0 = ffl.lvec;
    //printf("optimizeLattice_1d lvec0\n"    ); printMat( lvec0    );
    if(n1>0){
        //gopt.lattice_scan_1d( n1, lvec0, dlvec*-1, "lattice_scan_1d_bk.xyz", n1-1,-1 );
        gopt.lattice_scan_1d( n1, lvec0, dlvec*-1, 0, n1-1,-1 );
        setGeom( gopt.population[n1-1]->apos, &lvec0 );
    }
    if(n2>0){
        //gopt.lattice_scan_1d( n2, lvec0, dlvec   ,initMode, "lattice_scan_1d_fw.xyz", n1,1 );
        gopt.lattice_scan_1d( n2, lvec0+dlvec, dlvec   , 0, n1,1 );
        setGeom( gopt.population[n1-1]->apos, &lvec0 );
    }
    gopt.popToXYZ( "lattice_scan_1d_all.xyz");
    gopt.popToXYZ( "lattice_scan_1d_all_2x2.xyz",0,-1,{2,2,1});
    
}




virtual void upload_pop( const char* fname ){
    printf("MolWorld_sp3::upload_pop(%s) : We do lattice constant relaxation here \n", fname );
    gopt.loadPopXYZ( fname );

}

virtual void evalAFMscan( GridShape& scan, Quat4f*& OutFE, Quat4f*& OutPos, Quat4f** ps=0, bool bSaveDebug=false ){ printf( "MolWorld_sp3::evalAFMscan() NOT IMPLEMENTED, use GPU accelerated class MolWorld_sp3_multi instead! \n" ); }
virtual void evalAFM_FF ( GridShape& grid, Quat4f*  data=0,                                bool bSaveDebug=false ){ printf( "MolWorld_sp3::evalAFM_FF()  NOT IMPLEMENTED, use GPU accelerated class MolWorld_sp3_multi instead! \n" ); }

virtual void setSystemReplica (int i){ int nsys=countSystemReplica(); if(nsys<1)return; iSystemCur = i; printf( "MolWorld_sp3::setSystemReplica(%i/%i)\n", iSystemCur, nsys ); gopt.setGeom( iSystemCur ); };
virtual int countSystemReplica(     ){ return gopt.population.size(); }
void nextSystemReplica(){ int nsys=countSystemReplica(); int i=iSystemCur+1; if(i>=nsys)i=0;      setSystemReplica( i ); };
void prevSystemReplica(){ int nsys=countSystemReplica(); int i=iSystemCur-1; if(i<0    )i=nsys-1; setSystemReplica( i ); };

virtual char* getStatusString( char* s, int nmax ){
    s += sprintf(s, "iSystemCur %i\n",  iSystemCur );
    //if(bMMFF) s += sprintf(s, "eval:Ang,ps,ppT,ppI(%i|%i,%i,%i)\n",  ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0 );
    s += sprintf(s, "cog (%g,%g,%g)\n", cog .x,  cog .y,  cog .z );
    s += sprintf(s, "vcog(%15.5e,%15.5e,%15.5e)\n",  vcog.x,  vcog.y,  vcog.z);
    s += sprintf(s, "fcog(%15.5e,%15.5e,%15.5e)\n",  fcog.x,  fcog.y,  fcog.z);
    s += sprintf(s, "torq(%15.5e,%15.5e,%15.5e)\n", tqcog.x, tqcog.y, tqcog.z);
    return s;
}

// =================== Functions


virtual void swith_method(){ bGridFF=!bGridFF; };
virtual char* info_str   ( char* str=0 ){ if(str==0)str=tmpstr; sprintf(str,"bGridFF %i ffl.bAngleCosHalf %i \n", bGridFF, ffl.bAngleCosHalf ); return str; }

    // =================== PBC Shifts ( Should go to ForceField class ? )

/**
 * @brief Evaluates the periodic boundary condition (PBC) shifts for a given set of PBC and lattice vectors.
 * 
 * @param nPBC The dimensions of the PBC in each direction (x, y, z).
 * @param lvec The lattice vectors defining the unit cell.
 * @param shifts An array to store the calculated PBC shifts.
 * @return The total number of PBC shifts calculated.
 */
int evalPBCshifts( Vec3i nPBC, const Mat3d& lvec, Quat4f* shifts ){
    int ipbc=0;
    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
        if( (ix==0) && (iy==0) && (iz==0) ) ipbc0 = ipbc;
        shifts[ipbc].f = (Vec3f)( (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz) );   
        //printf( "shifts[%3i=%2i,%2i,%2i] (%7.3f,%7.3f,%7.3f)\n",  ipbc, ix,iy,iz, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z );
        ipbc++; 
    }}}
    return ipbc;
}

int evalPBCshifts( Vec3i nPBC, const Mat3d& lvec, Vec3d* shifts ){
    int ipbc=0;
    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){  
        shifts[ipbc] = (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz);   
        if( (ix==0) && (iy==0) && (iz==0) ) ipbc0 = ipbc;
        //printf( "shifts[%3i=%2i,%2i,%2i] (%7.3f,%7.3f,%7.3f)\n",  ipbc, ix,iy,iz, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z );
        ipbc++; 
    }}}
    return ipbc;
}

/**
 * @brief Calculates the periodic boundary condition (PBC) shifts for a given set of PBC and lattice vectors.
 * 
 * @param nPBC The dimensions of the PBC in each direction (x, y, z).
 * @param lvec The lattice vectors defining the unit cell.
 * @param shifts An array to store the calculated PBC shifts.
 * @return The total number of PBC shifts calculated.
 */
int makePBCshifts( Vec3i nPBC, const Mat3d& lvec ){
    npbc = (nPBC.x*2+1)*(nPBC.y*2+1)*(nPBC.z*2+1);
    //pbc_shifts = new Vec3d[npbc];
    _realloc(pbc_shifts,npbc);
    int npbc_eval = evalPBCshifts( nPBC, lvec, pbc_shifts );
    if(npbc!=npbc_eval){ printf( "ERROR in MolWorld_sp3::makePBCshifts() final ipbc(%i)!=nbpc(%i) => Exit()\n", npbc_eval,npbc ); exit(0); }
    return npbc;
}

void printPBCshifts(){
    printf("printPBCshifts():\n");
    for(int i=0; i<npbc; i++){ printf("pbc_shift[%i] (%6.3f,%6.3f,%6.3f)\n", i, pbc_shifts[i].x,pbc_shifts[i].y,pbc_shifts[i].z ); }
}

    // =================== Building / Modifying Molecule or ForceField


  /**

  /**
 * Substitutes a molecule in the builder.
 * 
 * @param fname The filename of the molecule to be substituted.
 * @param ib The index of the atom in the builder where the substitution will occur.
 * @param up The up vector for the substitution.
 * @param ipivot The index of the pivot atom for the substitution.
 * @param bSwapBond Flag indicating whether to swap the bond during substitution.
 * @param axSwap The axis swap vector for the substitution.
 * 
 * @return The index of the substituted molecule in the builder.
 */
    int substituteMolecule( const char* fname,  int ib, Vec3d up, int ipivot=0, bool bSwapBond=false, const Vec3i* axSwap=0 ){
        //builder.printAtomConfs(false);
        //builder.printBonds();
        printf( " ===================== Substitute molecule START !!! \n");
        Molecule* mol = new Molecule(); mol->init_xyz( fname, &params, true );
        //Vec3i axSwap={2,1,0};
        //Vec3i axSwap={2,0,1}; // THIS ONE
        //Vec3i axSwap={0,1,2};
        //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false );
        //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &axSwap, &debug_rot );
        int ja = builder.substituteMolecule( mol, Vec3dZ, ib, ipivot, false, 0, &debug_rot );
        //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &(Vec3i{2,1,0}), &debug_rot );
        //builder.addCappingTypesByIz(1);
        builder.tryAddConfsToAtoms( 0, -1 ); 
        builder.sortConfAtomsFirst();              
        builder.tryAddBondsToConfs( );      
        builder.finishFragment();       
        //builder.printAtomConfs(false);
        //builder.printBonds();
        //builder.printBondParams();
        delete mol;
        printf( "====================== Substitute molecule DONE !!! \n");
        return ja;
    }


/**
 * Finds the bridge bonds in the molecular world.
 * 
 * @return The number of bridge bonds found.
 */
    int findBridgeBonds(){
        //LimitedGraph graph;
        //graph.init( builder.na, builder.bonds.size(), builder.bonds.data() );
        builder.toLimitedGraph( graph );
        if(verbosity>2)graph.bPrint = true; // Debug
        graph.bridge();
        if(verbosity>0)printf( "MolWorld_sp3::findBridgeBonds()  graph.found.size() %i \n", graph.found.size() );
        return graph.found.size();
    }

/**
 * Multiplies the periodic boundary conditions (PBC) for a given fragment.
 * 
 * @param nMulPBC_ The vector representing the number of times to multiply the PBC in each direction (a, b, c)
 * @param ifrag The index of the fragment.
 */
    void PBC_multiply( Vec3i& nMulPBC_, int ifrag ){
        if(verbosity>0) printf( "PBC_multiply n(%i,%i,%i) ifrag=%i \n", nMulPBC_.x,nMulPBC_.y,nMulPBC_.z, ifrag );
        //printf("surface  lattice:\n"); gridFF .grid.cell.print();
        //printf("molecule lattice:\n"); builder.lvec.print();
        builder.multFragPBC( ifrag, nMulPBC_, builder.lvec );
        //printf("molecule lattice:\n"); builder.lvec.print();
        //builder.printAtoms();
        //new_lvec.ax=builder.lvec.a.norm(); new_lvec.by=builder.lvec.b.norm(); new_lvec.cz=builder.lvec.c.norm();
        builder.correctPBCbonds( ifrag, builder.frags.size() ); // correct bonds for newly added fragments
        builder.checkBondsInNeighs(true);
        builder.sortConfAtomsFirst(); 
        //printf("molecule lattice:\n"); builder.lvec.print();
        //builder.printAtomConfs();
        //builder.printBonds();    
    }

/**
 * Changes the cell by specifying two vectors, a and b.
 * Optionally, the function can also specify the index of an atom (ia0) and a shift vector (c0) to move the atoms.
 * If ia0 is provided, the atoms are shifted by the specified vector relative to the position of atom ia0.
 * 
 * @param a The first vector used to change the cell.
 * @param b The second vector used to change the cell.
 * @param ia0 (Optional) The index of the atom to use as a reference for shifting the atoms.
 * @param c0 (Optional) The shift vector to move the atoms.
 */
    void changeCellBySurf( Vec2d a, Vec2d b, int ia0=-1, Vec2d c0=Vec2dZero ){
        //printf( "changeCellBySurf() a(%g,%g) b(%g,%g) \n", a.x,a.y, b.x,b.y  );
        double la0=builder.lvec.a.norm();
        double lb0=builder.lvec.b.norm();
        Mat3d lvs;
        lvs.a=gridFF.grid.cell.a*a.a + gridFF.grid.cell.b*a.b;
        lvs.b=gridFF.grid.cell.a*b.a + gridFF.grid.cell.b*b.b; 
        lvs.c=builder.lvec.c;
        builder.changeCell( lvs );
        //Vec3d pmin,pmax; builder.bbox(pmin,pmax); printf( "BBOX pmin(%g,%g,%g) pmax(%g,%g,%g)\n", pmin.x,pmin.y,pmin.z,  pmax.x,pmax.y,pmax.z ); builder.move_atoms(pmin*-1);
        //builder.move_atoms( builder.atoms[0].pos*-1.);
        if(ia0>=0){
            Vec3d shift =  builder.atoms[ia0].pos*-1 + gridFF .grid.cell.a*c0.a + gridFF .grid.cell.b*c0.b;
            builder.move_atoms( shift );
        }
        printf( "changeCellBySurf() DONE, |a,b|=%g,%g (old |a,b|=%g,%g) \n", builder.lvec.a.norm(), builder.lvec.b.norm(), la0, lb0 );
        //builder.lvec = lvs;
        //builder.printAtomConfs();
        //builder.printBonds();    
    }

    int hideEPairs(){
        //printf( "plotNonBondGrid() removing EPairs %i \n" );
        int etyp = params.getAtomType("E");
        ffl.chargeToEpairs( -QEpair, etyp );
        selectByType( params.getElementType("E"), true );
        int nEp = selection.size();
        nbmol.natoms -= nEp;
        return nEp;
    }

    void unHideEPairs(){
        //if(nEp)
        //nbmol.natoms += nEp;
        nbmol.natoms = ffl.natoms;
        int etyp = params.getAtomType("E");
        ffl.chargeToEpairs( QEpair, etyp );
        //nEp = -1;
    }

/**
 * Calculates the auto charges for the molecules in the molecular world.
 * This function assigns charges to the atoms based on the given force field parameters.
 * It also applies constraints to certain atom types and performs charge relaxation using molecular dynamics.
 * After the charges are calculated, they are copied to the corresponding molecules in the molecular world.
 */
    void autoCharges(bool bVerbose=false, bool bFromScratch=false ){
        //if(verbosity>0)
        printf("MolWorld_sp3::autoCharges() ff.natoms=%i @nbmol.REQs=%li @ff.REQs=%li @ffu.REQs=%li @ffl.REQs=%li \n", ff.natoms, (long)nbmol.REQs, (long)ff.REQs, (long)ffu.REQs, (long)ffl.REQs  );
        //printf("MolWorld_sp3::autoCharges() START REQ.q[-1] \n", nbmol.REQs[nbmol.natoms-1].z );
        //bVerbose=true;

        ForceField* ff;
        if(bUFF){ ff = &ffu; }else{ ff = &ffl; }

        qeq.realloc( ff->natoms );
        params.assignQEq ( ff->natoms, ff->atypes, qeq.affins, qeq.hards );
        int etyp = params.getAtomType("E");    //Constrain electron pairs
        qeq.constrainTypes( ff->atypes, etyp );
        
        if(!bFromScratch){ copy( qeq.n, 4,2,(double*)nbmol.REQs, 1,0,(double*)qeq.qs ); }  // Initial charges 
        //for(int i=0; i<qeq.n; i++){ printf( "qeq.qs[%i]=%g  REQ.z=%g  \n", i, qeq.qs[i], nbmol.REQs[i].z );}
        qeq.relaxChargeMD ( ff->apos, 1000, 1e-2, 0.1, 0.0, bVerbose, bFromScratch );
        copy( qeq.n, 1,0,(double*)qeq.qs, 4,2,(double*)nbmol.REQs );
        if(bFromScratch){ ffl.chargeToEpairs( QEpair, etyp ); }
        if(bUFF){
            for(int i=0; i<ffu.natoms; i++){
                printf( "uff[ia=%i] atype[%i] REQs(%g,%g,%g)\n", i, ffu.atypes[i], nbmol.REQs[i].x, nbmol.REQs[i].y, nbmol.REQs[i].z );
            }
        }
        nbmol.evalPLQs(gridFF.alphaMorse);
        printf("MolWorld_sp3::autoCharges() END REQ.q[-1] \n", nbmol.REQs[nbmol.natoms-1].z );
        bChargeUpdated=true;
    }
                                             
                                            
    // =================== Initialization of different parts of the system ( different force-fields )


/**
 * Initializes the grid-based force field (GridFF)
 * 
 * @param name The name of the grid file.
 * @param bGrid Flag indicating whether to enable grid-based force field.
 * @param bSaveDebugXSFs Flag indicating whether to save debug XSF files.
 * @param z0 The z-coordinate of the cell origin.
 * @param cel0 The initial cell position.
 * @param bAutoNPBC Flag indicating whether to automatically set non-periodic boundary conditions.
 */
    virtual void initGridFF( const char * name, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bSymetrize=true, bool bAutoNPBC=true, bool bCheckEval=true, bool bUseEwald=true, bool bFit=true, bool bRefine=true ){
        //if(verbosity>0)
        printf("MolWorld_sp3::initGridFF(%s,bSymetrize=%i,bAutoNPBC=%i bCheckEval=%i bUseEwald=%i bFit=%i bRefine=%i bGridDouble=%i gridStep=%g,z0=%g,cel0={%g,%g,%g} )\n",  name, bSymetrize, bAutoNPBC,bCheckEval,bUseEwald,bFit,bRefine, bGridDouble, gridStep, z0, cel0.x,cel0.y,cel0.z  );
        sprintf(tmpstr, "%s.lvs", name );
        if( file_exist(tmpstr) ){  gridFF.grid.loadCell( tmpstr, gridStep );  gridFF.bCellSet=true; }
        if( !gridFF.bCellSet ){
            bGridFF=false; 
            printf( "WARRNING!!! GridFF not initialized because %s not found\n", tmpstr );
            return;
        }
        //double* ffgrid = 0;
        gridFF.grid.center_cell( cel0 );
        bGridFF=true;
        gridFF.bindSystem(surf.natoms, surf.atypes, surf.apos, surf.REQs );
        gridFF.initGridFF( name, z0, bAutoNPBC, bSymetrize );
        char wd0[1024]; getcwd(wd0,1024); printf( "MolWorld_sp3::initGridFF() 1 wd0=`%s`\n", wd0 );
        const char* last_slash = strrchr(name, '/');
        const char* result = (last_slash) ? last_slash + 1 : name;
        char wd[128]; sprintf( wd, "data/%s", result); printf( "MolWorld_sp3::initGridFF() 1 wd=`%s`\n", wd );
        tryMakeDir  ( wd );
        tryChangeDir( wd );
        gridFF.bUseEwald = bUseEwald;
        gridFF.ewald     = &gewald;
        gridFF.tryLoad_new( bSymetrize, bFit, bRefine );
        //ffgrid = gridFF.HHermite_d;
        tryChangeDir( wd0 );
        //getcwd(tmpstr, 1024 ); printf( "initGridFF() 3 WD=`%s`\n", tmpstr );
        gridFF.shift0 = Vec3d{0.,0.,-2.0};
        //gridFF.shift0 = Vec3d{0.,0.,0.0};
        //if(bCheckEval)gridFF.evalCheck();    // WARRNING:  CHECK FOR gridFF TURNED OFF !!!!!!!!!!!!!!!!!!!!!!!!!
        //return ffgrid;
    }

    void initNBmol( NBFF* ff, bool bCleanCharge=true ){
        if(verbosity>0)printf( "MolWorld_sp3::initNBmol() na %i \n", ff->natoms  );
        //void bindOrRealloc(int n_, Vec3d* apos_, Vec3d* fapos_, Quat4d* REQs_, int* atypes_ ){
        nbmol.bindOrRealloc( ff->natoms, ff->apos, ff->fapos, 0, ff->atypes );    
        //nbmol.bindOrRealloc( na, apos, fapos, 0, 0 );   
        //builder.export_atypes( nbmol.atypes );     
        builder.export_REQs( nbmol.REQs   );       ff->REQs=nbmol.REQs;
        nbmol  .makePLQs   ( gridFF.alphaMorse );  ff->PLQs=nbmol.PLQs; 
        nbmol  .makePLQd   ( gridFF.alphaMorse );  ff->PLQd=nbmol.PLQd; 
        //nbmol.print_nonbonded();
        if(bCleanCharge)for(int i=builder.atoms.size(); i<ff->natoms; i++){ nbmol.REQs[i].z=0; }  // Make sure that atoms not present in Builder has well-defined chanrge        
        params.assignREs( ff->natoms, nbmol.atypes, nbmol.REQs, true, false  );
        if(verbosity>1)nbmol.print();                              
    }


/**
 * Loads a molecular structure from an XYZ file.
 * 
 * @param name The name of the XYZ file.
 */
    void loadNBmol( const char* name){
        if(verbosity>0)printf( "MolWorld_sp3::loadNBmol() \n" );
        sprintf(tmpstr, "%s.xyz", name );
        params.loadXYZ( tmpstr, nbmol.natoms, &nbmol.apos, &nbmol.REQs, &nbmol.atypes );
        _realloc(nbmol.fapos,nbmol.natoms);
        nbmol  .makePLQs     ( gridFF.alphaMorse );  
        ffl.PLQs=nbmol.PLQs; 
        if(verbosity>1)nbmol.print();                              
    }

/**
 * @brief Loads a surface from a file and initializes the necessary data structures.
 * 
 * @param name The name of the file to load the surface from.
 * @param bGrid Flag indicating whether to generate a grid.
 * @param bSaveDebugXSFs Flag indicating whether to save debug XSF files.
 * @param z0 The z-coordinate of the surface.
 * @param cel0 The initial cell coordinates.
 * @return True if the surface is successfully loaded, false otherwise.
 */
    bool loadSurf(const char* name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0} ){
        char fname[256];
        sprintf(fname, "%s.xyz", name );
        int ret = params.loadXYZ( fname, surf.natoms, &surf.apos, &surf.REQs, &surf.atypes, 0, &gridFF.grid.cell );
        if     ( ret<0 ){ getcwd(tmpstr,1024); printf("ERROR in MolWorld_sp3::loadSurf() file(%s) not found in path(%s)=> Exit() \n", fname, tmpstr ); exit(0); }
        if     ( ret==0){ printf("ERROR in MolWorld_sp3::loadSurf() no lattice vectors in (%s) => Exit() \n", fname ); exit(0); }
        else if( ret>0 ){ gridFF.grid.updateCell(gridStep); gridFF.bCellSet=true;  }
        //gridFF.grid.printCell(); 
        if(verbosity>0)printf("MolWorld_sp3::loadSurf(%s) 1 natoms %i apos %li atyps %li \n", name, surf.natoms, (long)surf.apos, (long)surf.atypes  );
        //surf.print_nonbonded();
        bSurfAtoms=true;
        if(bGrid){
            bool bSymmetrize=true;
            initGridFF( name,z0,cel0, bSymmetrize );
        }
        return true;
    }

    int insertMolecule( const char* fname, const char* name=0, Vec3d pos={0,0,0}, Mat3d rot=Mat3dIdentity ){
        // ToDo: we should be able to insert molecule without actually creating molecule-type
        //sprintf(tmpstr, "%s.xyz", name );
        int iret=-1;
        // {   // Previous version using Molecule class
        //     if(name==0){ name=fname; }
        //     int imol  = builder.loadMolType( fname, name );
        //     if (imol<0){ printf("ERROR MolWorld_sp3::builder.loadMolType(%s) imol=%i \n", name, imol ); exit(0); }
        //     iret  = builder.insertFlexibleMolecule( imol, pos, Mat3dIdentity, -1 );
        //     if (iret<0){ printf("ERROR MolWorld_sp3::insertMolecule(%s) iret=%i \n", name, iret ); exit(0); }
        // }
        {   // New version using Atoms class
            iret = builder.loadXYZ_Atoms( fname, &params, -1, false, pos, rot );
        }
        //int ifrag = builder.frags.size()-1;
        //return ifrag;
        return iret;
    }

    void makeMoleculeTopology(){
        // ToDo: this will make topology for the whole system (all molecules), in future we should think how to make topology for individual molecules separately
        //builder.addCappingTypesByIz(1);
        builder.tryAddConfsToAtoms( 0, -1 );
        builder.cleanPis();
        //if(verbosity>2)builder.printAtomConfs(false);
        // --- we assume the lattice vectors are already loaded, and bPBC is set
        if( bPBC ){ builder.autoBondsPBC(); }
        else      { builder.autoBonds();    }
        //checkAllAtomsBonded( bool bPrint=true, bool bExit=true, int nbmax=N_NEIGH_MAX, int nbmin=1 );
        if(bCheckInit)builder.checkNumberOfBonds( true, true );
        //if(verbosity>2)builder.printBonds ();
    }

    void assingMoleculeTopoTypes(){
        if( fAutoCharges>0 )builder.chargeByNeighbors( true, fAutoCharges, 10, 0.5 );
        //if(substitute_name) substituteMolecule( substitute_name, isubs, Vec3dZ );
        //if( builder.checkNeighsRepeat( true ) ){ printf( "ERROR: some atoms has repating neighbors => exit() \n"); exit(0); };
        builder.autoAllConfEPi  (           ); 
        builder.setPiLoop       ( 0, -1, 10 ); // setup pi-orbitals
        if(bEpairs)builder.addAllEpairsByPi( );    
        //builder.printAtomConfs(false, false );
        //builder.printAtomConfs(false, true );
        // TBD here FF params are assigned already, but types are not yet found out...
        builder.assignAllBondParams();    //if(verbosity>1) 
    }

/**
 * Loads the geometry from a file and inserts it into the molecular world.
 * 
 * @param name The name of the file containing the geometry.
 * @return The index of the inserted fragment in the builder.
 */
    int loadGeom( const char* name ){ // TODO : overlaps with buildFF()
        if(verbosity>0)printf("MolWorld_sp3::loadGeom(%s)\n", name );
        // ------ Load geometry
        sprintf(tmpstr, "%s.xyz", name );
        /*
        int imol  = builder.loadMolType( tmpstr, name );
        int iret  = builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
        int ifrag = builder.frags.size()-1;
        */
        int ifrag = insertMolecule( tmpstr, name, pivotPoint, pivotRot );
        //builder.printAtomConfs(false, true );
        builder.addCappingTypesByIz(1);   // Find all hydrogen cappings
        builder.tryAddConfsToAtoms( 0, -1 );
        //builder.printAtomConfs(false, true );
        builder.cleanPis();
        if(verbosity>2)builder.printAtomConfs(false);
        // ------- Load lattice vectros
        sprintf(tmpstr, "%s.lvs", name );
        if( file_exist(tmpstr) ){
            if ( builder.bPBC ) printf("WARNING: Lattice vectors were already read from XYZ. Now they will be overwritten by the content of the file %s.\n", tmpstr);
            builder.bPBC=true;
            readMatrix( tmpstr, 3, 3, (double*)&builder.lvec );
        }
        bPBC=builder.bPBC;  //printf( "builder.bPBC %i \n", builder.bPBC );
        if( bPBC ){ builder.autoBondsPBC(); }
        else      { builder.autoBonds();    }
        if(bCheckInit)builder.checkNumberOfBonds( true, true );
        if(verbosity>2)builder.printBonds ();
        return ifrag;
    }

    // int loadmol(const char* fname_mol ){
    //     int imol = builder.loadMolType( fname_mol, "molecule" );
    //     builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
    //     return imol;
    // }

  /**
  /**
 * Inserts a SMILES string into the molecular world.
 * 
 * @param s The SMILES string to be inserted.
 */
    void insertSMILES(const char* s){
        smiles.builder=&builder;
        smiles.parseString( 10000, s );
    }


/**
 * Sets the optimizer for the molecular world.
 * 
 * @param n The number of degrees of freedom.
 * @param ps Pointer to the array of positions.
 * @param fs Pointer to the array of forces.
 */
    void setOptimizer( int n, double* ps, double* fs, double* vs=0 ){
        //opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
        opt.bindOrAlloc( n, ps, vs, fs, 0 );
        if(dt_default<0){ dt_default=ffl.optimalTimeStep(); } 
        if(verbosity>0)printf("MolWorld_sp3::setOptimizer(): optimal time_step = %g \n", dt_default);
        opt.initOpt( dt_default );
        opt.cleanVel();
        //opt.verbosity=2;
    }
    void setOptimizer(){ setOptimizer( ff.nDOFs, ff.DOFs, ff.fDOFs ); };

/**
 * Initializes the rigid bodies
 */
    void initRigid(){
        int nrb = builder.frags.size();
        //printf("# --- initRigid() nrb=%i \n", nrb);
        int n0rot=nrb*3;
        optRB.bindOrAlloc( n0rot + nrb*4, 0, 0, 0, 0);
        rbff.realloc( nrb, (Vec3d*)optRB.pos, (Quat4d*)(optRB.pos+n0rot), (Vec3d*)(optRB.force+n0rot), (Vec3d*)(optRB.vel+n0rot), 0, 0 );
        int natom=0;
        //printf("# --- initRigid() rbff.n=%i \n", rbff.n );
        for(int i=0; i<nrb; i++){
            const MM::Fragment& frag = builder.frags[i]; // problem - some atoms are not in builder - e.g. Epair
            int i0 = frag.atomRange.x;
            int ni = frag.atomRange.y - i0;
            //printf("# initRigid[%i] i0 %i ni %i \n", i, i0, ni );
            nbmol.apos + i0;
            rbff.mols[i].bindOrRealloc(ni, nbmol.apos+i0, nbmol.fapos+i0, nbmol.REQs+i0, nbmol.atypes+i0 );
            natom+=ni;
        }
        rbff.makePos0s();
        //printf("# --- initRigid() END \n");
    }
                                             
                                             
/**
 * @brief Initializes the MolWorld_sp3 object with a SMILES string.
 * 
 * @param s The SMILES string to initialize the object with.
 * @param bPrint Flag indicating whether to print debug information.
 * @param bCap Flag indicating whether to add all-cap topology.
 * @param bNonBonded_ Flag indicating whether to initialize non-bonded interactions. #spare
 * @param bOptimizer_ Flag indicating whether to set up an optimizer.
 */
    void initWithSMILES(const char* s, bool bPrint=false, bool bCap=true, bool bNonBonded_=false, bool bOptimizer_=true ){
        params.init("common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
        //params.printAtomTypeDict();
        //params.printAtomTypes();
        //params.printBond();
        builder.bindParams(&params);
        insertSMILES( s );
        if(bCap)builder.addAllCapTopo();
        //builder.autoAngles( 10.0, 10.0 );
        builder.randomizeAtomPos(1.0);
        builder.toMMFFsp3( ff );
        if(bPrint){   
            printf("=============\n"); printf("%s\n", s);
            ff.printBonds();
            ff.printNeighs();
        }
        //if(bNonBonded)init_nonbond();
        if(bOptimizer){ setOptimizer(); }
        _realloc( manipulation_sel, ff.natoms );
        printf( "... MolWorld_sp3::initWithSMILES() DONE\n" );
    }




// void ini_in_dir(){
//     params.init( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
//     builder.bindParams(&params);
//     int nheavy = 0;  // ---- Load Atomic Type Parameters
//     if( file_exist("cel.lvs") ){ 
//         loadGeom( "mm" ); 
//         if(bGridFF)makeGridFF();
//         // ----- Optimizer setup
//         //opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.fapos, 0 );
//         setOptimizer();
//         //double E = ff.eval(true);
//     }else{
// 		printf("WARNING: cel.lvs not found => molecular system not initialized in [MolWorld_sp3::ini_in_dir()] \n" );
// 	}
// }

/**
 * Initializes the parameters for the MolWorld_sp3 class.
 * 
 * @param sElemTypes    The string containing the element types.
 * @param sAtomTypes    The string containing the atom types.
 * @param sBondTypes    The string containing the bond types.
 * @param sAngleTypes   The string containing the angle types.
 * @param sDihedralTypes    The string containing the dihedral types (optional).
 */
    void initParams( const char* sElemTypes, const char* sAtomTypes, const char* sBondTypes, const char* sAngleTypes, const char* sDihedralTypes=0 ){
        printf( "MolWorld_sp3::initParams():\n\tsElemTypes(%s)\n\tsAtomTypes(%s)\n\tsBondTypes(%s)\n\tsAngleTypes(%s)\n", sElemTypes, sAtomTypes, sBondTypes, sAngleTypes );
        params.init( sElemTypes, sAtomTypes, sBondTypes, sAngleTypes, sDihedralTypes );
        builder.bindParams(&params);
        params_glob = &params;
        builder.capAtomEpair.type = params.getAtomType("E");
        builder.addCappingTypesByIz(1);   // hydrogens
        builder.addCappingTypesByIz(200); // electron pairs
        //params.printAtomTypeDict();
        //params.printAtomTypes();
        //params.printBond();
        //params.printAngleTypes();
        //params.printDihedralTypes();
        builder.sp3types.insert( params.getAtomType("C") );
        builder.sp3types.insert( params.getAtomType("N") );
        builder.sp3types.insert( params.getAtomType("O") );
    }


/**
 * Builds a molecule from an XYZ file.
 * 
 * @param xyz_name The name of the XYZ file.
 * @return The index of the built molecule.
 */
    int buildMolecule_xyz( const char* xyz_name ){
        int ifrag = loadGeom( xyz_name );
        printf( "MolWorld_sp3::buildMolecule_xyz(%s) ifrag=%i \n", xyz_name, ifrag );
        int ia0=builder.frags[ifrag].atomRange.a;
        int ic0=builder.frags[ifrag].confRange.a;
        // TBD not sure that I got how charges are assigned in here...
        if( fAutoCharges>0 )builder.chargeByNeighbors( true, fAutoCharges, 10, 0.5 );
        if(substitute_name) substituteMolecule( substitute_name, isubs, Vec3dZ );
        if( builder.checkNeighsRepeat( true ) ){ printf( "ERROR: some atoms has repating neighbors => exit() \n"); exit(0); };
        builder.autoAllConfEPi  ( ia0 );
        builder.setPiLoop       ( ic0, -1, 10 );
        if(bEpairs)builder.addAllEpairsByPi( ia0=0 ); 
        //builder.printAtomConfs(false, false );
        //builder.printAtomConfs(false, true );
        // TBD here FF params are assigned already, but types are not yet found out...
        builder.assignAllBondParams();    //if(verbosity>1)
        builder.finishFragment(ifrag);    
        //builder.printAtoms();
        //printf( "buildMolecule_xyz: nMulPBC(%i,%i,%i) \n",nMulPBC.x,nMulPBC.y,nMulPBC.z  );
        //if( nMulPBC    .totprod()>1 ){ PBC_multiply    ( nMulPBC, ifrag ); };
        //if( bCellBySurf             ){ changeCellBySurf( bySurf_lat[0], bySurf_lat[1], bySurf_ia0, bySurf_c0 ); };
        //printf("builder.lvec\n");builder.lvec.print();
        return ifrag;
    }

/**
 * Prepare the molecular world for MMFF calculations.
 * This includes checking bonds, assigning atom types, assigning torsions, finding bridge bonds,
 * converting to MMFF sp3 format, assigning angles, converting to MMFF f4 format, and setting up
 * periodic boundary conditions if enabled.
 */
    void makeMMFFs(){
        print("MolWorld_sp3::makeMMFFs()\n" );
        // check if all bonds are in atom neighbors
        if( builder.checkBondsInNeighs(true) ) { 
            printf("ERROR some bonds are not in atom neighbors => exit"); 
            exit(0); 
        };
        // reshuffling atoms in order to have non-capping first
        builder.numberAtoms();
        builder.sortConfAtomsFirst();
        builder.checkBondsOrdered( true, false );

        // make assignement of atom types and force field parameters
        if( bUFF ){  // according to UFF
            builder.assignUFFtypes( 0, bCumulene, true, b141, bSimple, bConj); 
            builder.assignUFFparams( 0, true );
        }else{      // according to MMFF
            builder.assignTypes();
        }

        // passing them to FFs
        if ( bUFF ){
            builder.toUFF( ffu, true );
        }else{
            if( ffl.bTorsion ){ builder.assignTorsions( true, true ); }  //exit(0);

            builder.printAtomConfs();
            builder.printBonds();

            builder.toMMFFsp3_loc( ffl, true, bEpairs, bUFF );   
            //ffl.printAtomParams();
            if(ffl.bTorsion){  ffl.printTorsions(); } // without electron pairs
            if(ffl.bEachAngle){ builder.assignAnglesMMFFsp3  ( ffl, false      ); ffl.printAngles();   }  //exit(0);
            //builder.toMMFFf4     ( ff4, true, bEpairs );  //ff4.printAtomParams(); ff4.printBKneighs(); 
            builder.toMMFFsp3    ( ff , true, bEpairs );
            ffl.flipPis( Vec3dOne );

            ffl.printNeighs();
            //ff4.flipPis( Vec3fOne );

            //ffl.printAtomParams();
            //ffl.print_nonbonded();
        }

        // setting up PBC
        if(bPBC){
            if ( bUFF ){
                ffu.setLvec( builder.lvec);   
                npbc = makePBCshifts( nPBC, builder.lvec );
                ffu.bindShifts(npbc,pbc_shifts);
                //ffu.makeNeighCells( nPBC );      
                ffu.makeNeighCells( npbc, pbc_shifts ); 
            }else{
                ff.bPBCbyLvec = true;
                ff .setLvec( builder.lvec);
                ffl.setLvec( builder.lvec);   
                //ff4.setLvec((Mat3f)builder.lvec);
                npbc = makePBCshifts( nPBC, builder.lvec );
                ffl.bindShifts(npbc,pbc_shifts);
                //ff4.makeNeighCells  ( nPBC );       
                //ffl.makeNeighCells( nPBC );      
                ffl.makeNeighCells( npbc, pbc_shifts ); 
            }
        }

    }

/**
 * This function is responsible for performing various operations to generate force fields (FFs).
 * It calls the makeMMFFs() function, initializes the non-bonded molecules, sets the non-bonded interactions,
 * and performs additional calculations and checks.
 * If the bChargeToEpair flag is set to true, it converts charges to electron pairs.
 * Finally, if the bOptimizer flag is set to true, it sets the optimizer and performs relaxation if bRelaxPi is also true.
 */
    virtual void makeFFs(){
        print("MolWorld_sp3::makeFFs()\n" );
        makeMMFFs();
        if ( bUFF ){
            initNBmol( &ffu );
            setNonBond( bNonBonded );
            ffu.go = &go;
            nbmol.evalPLQs(gridFF.alphaMorse);
            ffu.atomForceFunc = [&](int ia,const Vec3d p,Vec3d& f)->double{    
                //printf( "ffu.atomForceFunc() ia=%i \n", ia  );
                double E=0;
                if   (bGridFF){ 
                    E += gridFF.addAtom( p, nbmol.PLQd[ia], f );
                    //Vec3d fi=Vec3dZero;
                    //E+= gridFF.addForce( p, nbmol.PLQs[ia], fi, true  ); 
                    //E += gridFF.addAtom( p, nbmol.PLQd[ia], fi );
                    //printf("MolWorld_sp3::ffu.atomForceFunc(ia=%i,gridFF.mode=%i) p(%g,%g,%g) fi(%g,%g,%g) PLQ(%g,%g,%g) @gridFF.Bspline_PLQ=%li \n", ia, (int)gridFF.mode, p.x, p.y, p.z, fi.x,fi.y,fi.z, nbmol.PLQs[ia].x, nbmol.PLQs[ia].y, nbmol.PLQs[ia].z, (long)gridFF.Bspline_PLQ );
                    //f.add( fi );
                }  // GridFF
                if(bConstrZ){
                    springbound( p.z-ConstrZ_xmin, ConstrZ_l, ConstrZ_k, f.z );
                }
                if(ipicked==ia)[[unlikely]]{ 
                    const Vec3d fs = getForceSpringRay( p, pick_hray, pick_ray0,  Kpick ); 
                    f.add( fs );
                }
                
                return E;
            };

            if(bOptimizer){ 
                //setOptimizer( ffu.nDOFs, ffu.DOFs, ffu.fDOFs );
                setOptimizer( ffu.natoms*3, (double*)ffu.apos, (double*)ffu.fapos );
                ffu.vapos = (Vec3d*)opt.vel;
            }                         
        }else{
            initNBmol( &ffl );
            //ffl.printAtomParams();
            setNonBond( bNonBonded );
            //ffl.print_nonbonded();
            bool bChargeToEpair=true;
            //bool bChargeToEpair=false;
            if(bChargeToEpair){
                int etyp=-1; etyp=params.atomTypeDict["E"];
                //ff.chargeToEpairs( nbmol.REQs, QEpair, etyp );  
                ffl.chargeToEpairs( QEpair, etyp ); 
            }
            nbmol.evalPLQs(gridFF.alphaMorse);

            if(bCheckInit){
                idebug=1;
                ffl.checkREQlimits();
                ffl.eval_check();
                //ff4.eval_check();
                //ff .eval_check();
                idebug=0;
            }

            //ffl.print_nonbonded(); exit(0);
            if(bOptimizer){ 
                //setOptimizer(); 
                //setOptimizer( ff.nDOFs, ff .DOFs,  ff.fDOFs );
                setOptimizer( ffl.nDOFs, ffl.DOFs, ffl.fDOFs );
                if(bRelaxPi) ffl.relax_pi( 1000, 0.1, 1e-4 );
                ffl.vapos = (Vec3d*)opt.vel;
            }                         
            _realloc( manipulation_sel, ff.natoms );
        }
        //ffl.print_nonbonded();
    }

virtual void updateFromBuilder(){
    printf("MolWorld_sp3::updateFromBuilder()\n");
    clearFFs();
    makeFFs();
};

void clearFFs(){
    printf("MolWorld_sp3::clearFFs()\n");
    ffl.dealloc();
    ff.dealloc();
    //ff4.dealloc();
    _dealloc( apos_bak );
    // --- nbmol
    nbmol.neighs=0;   // NOTE : if we set pointer to zero it does not try to deallocate it !!!
    nbmol.apos=0;  
    nbmol.fapos=0;  
    nbmol.atypes=0;
    nbmol.dealloc();
    // --- opt
    opt.pos = 0;
    opt.force = 0;
    opt.dealloc();
}


/**
 * @brief Clears the MolecularWorld object.
 * 
 * This function clears the MolecularWorld object by deallocating memory and resetting variables.
 * 
 * @param bParams Flag indicating whether to clear the parameters as well. Default is true.
 */
virtual void clear( bool bParams=true, bool bSurf=false ){
    //printf("MolWorld_sp3::clear() \n");
    builder.clear();
    selection.clear();
    selection_set.clear();
    clearFFs();
    constrs.clear();
    if(bSurf){  // ToDo: deallocate gridFF and surf ?
        //gridFF.dealloc();
        //gridFF.bCellSet=false;
        //gridFF.bGridFF=false;
    }
    if(bParams){
        params.clear();
    }
    // if(database){
    //     printf("MolWorld_sp3::clear() database->clear(); \n");
    //     database->dealloc();
    // }
}

    virtual int getMultiSystemPointers( int*& M_neighs,  int*& M_neighCell, Quat4f*& M_apos, int& nvec ){
        // int nsys=0,nvec=0;
        // int    * M_neighs    =0;
        // int    * M_neighCell =0;
        // Quat4f * M_apos     =0;
        return 0;
    }


/**
 * Scans the surface with force field, changing forces.
 * 
 * @param n The number of particles.
 * @param ps An array of Quat4f representing the particle positions.
 * @param REQs An array of Quat4f representing non-bonded interactions.
 * @param fs An array of Quat4f to store the calculated forces.
 */
    virtual void scanSurfFF( int n, Quat4f* ps, Quat4f* REQs, Quat4f* fs ){
        for(int i=0; i<n; i++){
            Quat4f PLQ = REQ2PLQ        ( (Quat4d)REQs[i], gridFF.alphaMorse );
            fs[i]      = gridFF.getForce( (Vec3d)ps[i].f, PLQ, true );
        }
    }


/**
 * Checks the invariants of the molecular world.
 * 
 * @param maxVcog The maximum value for the center of gravity velocity.
 * @param maxFcog The maximum value for the center of gravity force.
 * @param maxTg   The maximum value for the torque of the center of gravity.
 * @return True if any of the invariants exceed the maximum values, false otherwise.
 */
    bool checkInvariants( double maxVcog, double maxFcog, double maxTg ){
        cog   = average( ff.natoms, ff.apos  );
        vcog  = sum    ( ff.natoms, (Vec3d*)opt.vel  );
        fcog  = sum    ( ff.natoms, ff.fapos );
        tqcog = torq   ( ff.natoms, ff.apos, ff.fapos, cog );
        //tqcog.add( ff.evalPiTorq() );
        return ( vcog.norm()>maxVcog ) || ( fcog.norm()>maxFcog ) || ( tqcog.norm() );
    }

    //void open_xyzFile (const char* fname){ xyz_file=fopen( fname,"w" ); };
    //void close_xyzFile(){fclose(xyz_file)};

    // double eval_f4(){
    //     pack( ff4.natoms, ffl.apos , ff4.apos  );
    //     pack( ff4.nnode,  ffl.pipos, ff4.pipos );
    //     double E = ff4.eval();
    //     //ff4.move_GD( 0.01);
    //     unpack( ff4.natoms, ffl. apos, ff4. apos  );
    //     unpack( ff4.natoms, ffl.fapos, ff4.fapos  );
    //     unpack( ff4.nnode,  ffl. pipos,ff4. pipos );
    //     unpack( ff4.nnode,  ffl.fpipos,ff4.fpipos );
    //     //for(int i=0; i<ff4.nnode; i++) printf("pi[%i] <fpi,pi> %g |pi| %g \n", i, ffl.fpipos[i].dot( ffl.pipos[i] ), ffl.pipos[i].norm() );
    //     return E;   
    // };


/**
 * Sets the non-bonded flag for the molecular world.
 * 
 * @param bNonBonded A boolean value indicating whether non-bonded interactions should be considered.
 */
    void setNonBond( bool bNonBonded ){
        ffl.bSubtractAngleNonBond = bNonBonded;
        //ff4.bSubtractAngleNonBond = bNonBonded;
        // if(bNonBonded){
        //     ffl.REQs = nbmol.REQs;
        //     ff .REQs = nbmol.REQs;
        //     // if(ff4.REQs==0){
        //     //     ff4.REQs = new Quat4f[nbmol.natoms];
        //     //     for(int i=0; i<nbmol.natoms; i++ ){ ff4.REQs[i] = (Quat4f)nbmol.REQs[i]; };
        //     // }
        // }
    }



    int updateBuilderFromFF(bool bPos=true, bool bQ=true){
        if( ffl.nnode  != builder.confs.size() ){printf( "ERROR: MolWorld_sp3::updateBuilderFromFF() ffl.nnode(%i)  != builder->confs.size(%i) \n", ffl.nnode,  builder.confs.size() ); exit(0); }
        //if( ffl.natoms != builder.atoms.size() ){printf( "ERROR: MolWorld_sp3::updateBuilderFromFF() ffl.natoms(%i) != builder->atoms.size(%i) \n", ffl.natoms, builder.atoms.size() ); exit(0); }
        //printf( "MolWorld_sp3::updateBuilderFromFF(nnode=%i,ncap=%i) \n", ffl.nnode, ffl.natoms-ffl.nnode );
        int na = builder.atoms.size();  if(na>ffl.natoms)na=ffl.natoms;
        for(int i=0; i<na; i++){
            if(bPos){ builder.atoms[i].pos   = ffl.apos[i];   }
            if(bQ  ){ builder.atoms[i].REQ.z = ffl.REQs[i].z; }
        }
        return 0;
    }

    double evalEkTemp(){
        double Ek = ffl.evalKineticEnergy();
        int nDOFs = ffl.natoms*3;
        //nDOFs += ffl.nnode*2; // pi-rotations
        go.T_current = ( 2*Ek/ (const_kB*nDOFs) );
        //printf( "MolWorld_sp3::evalEkTemp() T=%g[K] T_target=%g[K] gamma=%g[1/dtu] Ek=%g[eV] nDOFs=%i \n", T_current, T_target, gamma_damp, Ek, nDOFs );
        return go.T_current;
    }

    void update_GOpt(){

        go.nExplore = 0;
    }

    double whichAtomNotConv( int& imax ){
        double F2max=0;
        imax =-1;
        for(int ia=0; ia<ffl.natoms; ia++){
            double r2 = ffl.fapos[ia].norm2();
            if(r2>F2max){ F2max=r2; imax=ia; }
        }
        return sqrt(F2max);
    }

    double getMostDisplacedAtom( Vec3d* ps, int& imax ){
        imax=-1;
        double R2max =  0;
        for(int ia=0; ia<ffl.natoms; ia++){
            Vec3d d=ffl.apos[ia]-ps[ia];
            double r2 = d.norm2();
            if(r2>R2max){ R2max=r2; imax=ia; }
        }
        return sqrt(R2max);
   }

    bool checkStuck( double R ){
        double R2 = R*R;
        int imax  = -1;
        double R2max = getMostDisplacedAtom( apos_bak, imax );
        if(R2max>R2)[[unlikely]]{ 
            for(int ia=0; ia<ffl.natoms; ia++){ apos_bak[ia]=ffl.apos[ia]; }
            nStuck=0;
           if(atomTrjFile)fclose(atomTrjFile); 
            atomTrjFile=0;
            return true;
        }
        if( (atomTrjFile==0) && ( nStuck>=(nStuckMax-nStuckTrj) ) ){
            atomTrjFile=fopen( "StuckAtomTrj.log", "w" );
        }
        if( nStuck>nStuckMax )[[unlikely]]{
            printf( "MolWorld_sp3::checkStuck() nStuck(%i)>nStuckMax(%i) => exit(0) \n", nStuck, nStuckMax );
            //saveXYZ( "stuck.xyz", ffl.natoms, ffl.apos, ffl.atypes );
            saveXYZ( "stuck.xyz", tmpstr, false, "a", nPBC_save );
            if(atomTrjFile)fclose(atomTrjFile);
            exit(0);
        }
        nStuck++;
        return false;
    }

    void handleStuckAtom(int itr, Vec3d cvf ){
        if( nStuck>=(nStuckMax-nStuckTrj) ){
            int imax; double Fmax = whichAtomNotConv( imax );
            Vec3d pi = ffl.apos [imax];
            Vec3d fi = ffl.fapos[imax];
            Vec3d vi = ffl.vapos[imax];
            if(atomTrjFile){
                double f = sqrt(cvf.z);
                double v = sqrt(cvf.y);
                fprintf( atomTrjFile, "%4i %4i   %20.15f %20.15f %20.15f   %20.15f %20.15f %20.15f    %20.15f %20.15f %20.15f    %6.4f %20.15f %20.15f \n", itr, imax, pi.x,pi.y,pi.z, vi.x,vi.y,vi.z, fi.x,fi.y,fi.z,    cvf.x/(f*v), v, f );
            }
        }
        //if( nStuck>(nStuckMax-2) ){
        // printf( "MolWorld_sp3::run_no_omp() Stuck due to atom[%i] |F[%i]|=%g fapos[%i](%g,%g,%g) apos[%i](%g,%g,%g) \n", imax, imax, Fmax, imax, fi.x,fi.y,fi.z, imax, pi.x,pi.y,pi.z );
        // gridFF.getEFprofileToFile( "Stuck_Fx.txt", 100, pi-Vec3dX*0.1, pi+Vec3dX*0.1, ffl.REQs[imax] );
        // gridFF.getEFprofileToFile( "Stuck_Fy.txt", 100, pi-Vec3dY*0.1, pi+Vec3dY*0.1, ffl.REQs[imax] );
        // gridFF.getEFprofileToFile( "Stuck_Fz.txt", 100, pi-Vec3dZ*0.1, pi+Vec3dZ*0.1, ffl.REQs[imax] );
        // exit(0);
        //}
    }

    __attribute__((hot))  
    double eval( ){
        if(verbosity>0)[[unlikely]]{ printf( "#### MolWorld_sp3::eval()\n"); }
        //ffl.doBonds       = false;
        //ffl.doPiPiI       = false;
        //ffl.doPiSigma     = false;
        //ffl.doAngles      = false;
        //ffl.bAngleCosHalf = false;
        //ffl.bEachAngle    = true;
        //printf("MolWorld_sp3::eval() bConstrains %i bNonBonded %i ffl.bSubtractAngleNonBond %i  ffl.bPBC %i ffl.doBonds %i ffl.doPiPiI %i ffl.doPiSigma %i ffl.doAngles %i ffl.bAngleCosHalf %i ffl.bEachAngle %i \n", bConstrains, bNonBonded, ffl.bSubtractAngleNonBond,   ffl.bPBC, ffl.doBonds, ffl.doPiPiI, ffl.doPiSigma, ffl.doAngles, ffl.bAngleCosHalf, ffl.bEachAngle );
        double E=0;
        //setNonBond( bNonBonded );  // Make sure ffl subtracts non-covalent interction for angles
        //ffl.print_nonbonded();
        //ffl.printAtomParams();
        //ffl.print_pbc_shifts();
        //printf("lvec: ");printMat(builder.lvec);
        if(bMMFF){ 
            if(bUFF){ E += ffu.eval(); }
            else{ E += ffl.eval(); }
            
        }else{ VecN::set( nbmol.natoms*3, 0.0, (double*)nbmol.fapos );  }      
        //bPBC=false;
        if(bNonBonded){
            //E += nbmol.evalLJQs_ng4_PBC_omp( );
            E += ffl  .evalLJQs_ng4_PBC_omp( );
            /*
            if(bMMFF){    
                if  (bPBC){ E += nbmol.evalLJQs_ng4_PBC( ffl.neighs, ffl.neighCell, npbc, pbc_shifts, gridFF.Rdamp ); }   // atoms outside cell
                else      { E += nbmol.evalLJQs_ng4    ( ffl.neighs );                                   }   // atoms in cell ignoring bondede neighbors       
                //else      { E += nbmol.evalLJQs_ng4_omp( ffl.neighs );                                   }   // atoms in cell ignoring bondede neighbors  
            }else{
                if  (bPBC){ E += nbmol.evalLJQs_PBC    ( ff.lvec, {1,1,0} ); }   // atoms outside cell
                else      { E += nbmol.evalLJQs        ( );                  }   // atoms in cell ignoring bondede neighbors    
            }
            */
        }
        //printf( "bConstrains=%i constrs.bonds.size()=%i \n", bConstrains, constrs.bonds.size() );
        if(bConstrains)constrs.apply( nbmol.apos, nbmol.fapos, &ffl.lvec );
        /*
        if(bSurfAtoms){ 
            if   (bGridFF){ E+= gridFF.eval(nbmol.natoms, nbmol.apos, nbmol.PLQs, nbmol.fapos ); }
            //else        { E+= nbmol .evalMorse   ( surf, false,                  gridFF.alphaMorse, gridFF.Rdamp );  }
            else          { E+= nbmol .evalMorsePBC( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );  }
        }
        */
        //printf( "eval() bSurfAtoms %i bGridFF %i \n", bSurfAtoms, bGridFF );
        //for(int i=0; i<nbmol.natoms; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fapos[i].x,nbmol.fapos[i].y,nbmol.fapos[i].z ); }    
        //ffl.printDebug(  false, false );
        //exit(0);
        if(verbosity>0)[[unlikely]]{ printf( "#### MolWorld_sp3::eval() DONE\n\n"); }
        return E;
    }







   __attribute__((hot))  
    int run_omp_Milan( int niter_max, double dt, double Fconv=1e-6, double Flim=1000, double timeLimit=0.02, double* outE=0, double* outF=0, double* outV=0, double* outVF=0 ){
        if(dt>0){ opt.setTimeSteps(dt); }else{ dt=opt.dt; }
        //printf( "run_omp() niter_max %i dt %g Fconv %g Flim %g timeLimit %g outE %li outF %li \n", niter_max, dt, Fconv, Flim, timeLimit, (long)outE, (long)outF );
        long T0 = getCPUticks();
        double E=0,F2=0,F2conv=Fconv*Fconv;
        double ff=0,vv=0,vf=0;
        int itr=0,niter=niter_max;
        bConverged = false;
        if(bToCOG && (!bGridFF) ){ 
            Vec3d cog=average( ffl.natoms, ffl.apos );  
            move( ffl.natoms, ffl.apos, cog*-1.0 ); 
        }
        if(bCheckStuck){ checkStuck( RStuck ); }
        //#pragma omp parallel shared(E,F2,ff,vv,vf,ffl) private(itr)
        #pragma omp parallel shared(niter,itr,E,F2,ff,vv,vf,ffl,T0,bConstrains,bConverged)
        while(itr<niter || bConverged) {
            if(itr<niter){
            //#pragma omp barrier
            #pragma omp single
            {E=0;F2=0;ff=0;vv=0;vf=0;
            }

            //------ eval forces
            //#pragma omp barrier
            #pragma omp for reduction(+:E)
            for(int ia=0; ia<ffl.natoms; ia++){ 
                {                 ffl.fapos[ia           ] = Vec3dZero; } // atom pos force
                if(ia<ffl.nnode){ ffl.fapos[ia+ffl.natoms] = Vec3dZero; } // atom pi  force

if(std::isnan(E)){printf("Before eval_atom\n");exit(1);}
                if(ia<ffl.nnode){ E+=ffl.eval_atom(ia); }


if(std::isnan(E)){printf("After eval_atom\n");exit(1);}
                // ----- Error is HERE
                if(bPBC){ E+=ffl.evalLJQs_ng4_PBC_atom_omp( ia ); }
                else    { E+=ffl.evalLJQs_ng4_atom_omp    ( ia ); } 
if(std::isnan(E)){printf("After evalLJQs_ng4_atom_omp\n");exit(1);}
bGridFF=false;
                if   (bGridFF){ E+= gridFF.addForce          ( ffl.apos[ia], ffl.PLQs[ia], ffl.fapos[ia], true ); }        // GridFF

                if(ipicked==ia)[[unlikely]]{ 
                    const Vec3d f = getForceSpringRay( ffl.apos[ia], pick_hray, pick_ray0,  Kpick ); 
                    ffl.fapos[ia].add( f );
                }

                if(bConstrZ){
                    springbound( ffl.apos[ia].z-ConstrZ_xmin, ConstrZ_l, ConstrZ_k, ffl.fapos[ia].z );
                }

            }
            // if(ffl.bTorsion){
            //     #pragma omp for reduction(+:E)
            //     for(int it=0; it<ffl.ntors; it++){ 
            //         E+=ffl.eval_torsion(it); 
            //     }
            // }
            #pragma omp single
            {
                if(bConstrains  ){ /*std::cout<<"E before constrains :" << E;*/ E+=constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec ); /*std::cout<<" E after constrains :" << E <<std::endl;*/ }
                if(!bRelax){ gopt.constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec ); }
            }
            // ---- assemble (we need to wait when all atoms are evaluated)
            //#pragma omp barrier
            #pragma omp for
            for(int ia=0; ia<ffl.natoms; ia++){
                ffl.assemble_atom( ia );
            }
            
            //#pragma omp barrier
            { //  ==== FIRE
                #pragma omp for reduction(+:vf,vv,ff)
                for(int i=0; i<opt.n; i++){
                    double v=opt.vel  [i];
                    double f=opt.force[i];
                    vv+=v*v; ff+=f*f; vf+=v*f;
                }
                #pragma omp single
                { opt.vv=vv; opt.ff=ff; opt.vf=vf; F2=ff; opt.FIRE_update_params(); }
                // ------ move
                #pragma omp for
                for(int i=0; i<ffl.nvecs; i++){
                    if( bRelax ){
                        ffl.move_atom_FIRE( i, opt.dt, 10000.0, opt.cv, opt.renorm_vf*opt.cf );
                        //ffl.move_atom_FIRE( i, dt, 10000.0, 0.9, 0 ); // Equivalent to MDdamp
                    }
                }
                sprintf(tmpstr,"# %i E %g |F| %g istep=%i", gopt_ifound, Etot, sqrt(ffl.cvf.z), go.istep );
            }
            
            } 
            
            //#pragma omp barrier
            #pragma omp single
            { 
                Etot=E;
                itr++; 
                if(timeLimit>0)[[unlikely]]{
                    double t = (getCPUticks() - T0)*tick2second;
                    if(t>0.02){ 
                        niter=0; 
                        if(verbosity>1) [[unlikely]] { printf( "run_omp() ended due to time limit after %i nsteps ( %6.3f [s]) \n", itr, t ); }
                    }
                }
                if(F2<F2conv)[[unlikely]]{ 
                    niter=0; 
                    bConverged = true;
                    double t = (getCPUticks() - T0)*tick2second;
                    if(verbosity>1) [[unlikely]] { printf( "run_omp() CONVERGED in %i/%i nsteps E=%g |F|=%g time= %g [ms]( %g [us/%i iter])\n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr ); }

                    if(verbosity) printf( "run_omp() CONVERGED in %i/%i nsteps E=%g |F|=%g time= %g [ms]( %g [us/%i iter])\n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr );
                    if(bGopt  && bRelax ){
                        gopt_ifound++;
                        sprintf(tmpstr,"# %i E %g |F| %g istep=%i", gopt_ifound, Etot, sqrt(ffl.cvf.z), go.istep );
                        saveXYZ( "gopt.xyz", tmpstr, false, "a", nPBC_save );printf( "run_omp().save %s \n", tmpstr );//}
                    }
                }
            }
           
        }{
        double t = (getCPUticks() - T0)*tick2second;
        if( (itr>=niter_max)&&(verbosity>1)) [[unlikely]] {printf( "run_omp() NOT CONVERGED in %i/%i dt=%g E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, opt.dt, E,  sqrt(F2), t*1e+3, t*1e+6/itr, itr ); }
        }
        return itr;
    }




































  
/**
 * Performs relaxation of the molecular system with FIRE algorithm.
 * 
 * @param niter The number of relaxation iterations to perform.
 * @param Ftol The force tolerance for convergence. Defaults to 1e-6.
 * @param bWriteTrj Flag indicating whether to write trajectory information. Defaults to false.
 * @return True if relaxation converged, false otherwise.
 */
 __attribute__((hot))  
bool relax( int niter, double Ftol = 1e-6, bool bWriteTrj=false ){
    printf( "MolWorld_sp3::relax() niter %i Ftol %g bWriteTrj %i \n", niter, Ftol, bWriteTrj );
    Etot=0.0;
    double f2tol=Ftol*Ftol;
    bConverged=false; 
    if(bWriteTrj){ xyz_file=fopen( "relax_trj.xyz","w" ); }
    for(int itr=0; itr<niter; itr++){
        Etot=eval();                                                  
        if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
        double f2 = opt.move_FIRE();
        //if(bWriteTrj){ toXYZ(); ;printf("DEBUB[%i] 4 \n", itr); };
        if(bWriteTrj){  sprintf(tmpstr,"# relax[%i] E=%g f2=%g", itr, Etot, sqrt(f2) );  toXYZ(tmpstr); };
        printf( "relax[%i] |F| %g (Ftol=%g)  Etot %g \n", itr, sqrt(f2), Ftol, Etot );
        if(f2<f2tol){ bConverged=true; break; }
    }
    if(bWriteTrj){ fclose(xyz_file); }
    return bConverged;
}

/**
 * Runs the simulation for a specified number of steps.
 * 
 * @param nstepMax The maximum number of steps to run.
 * @param dt The time step size. Default value is -1.
 * @param Fconv The convergence threshold for the force. Default value is 1e-6.
 * @param ialg The algorithm to use for optimization. Default value is 2.
 * @param outE Pointer to an array to store the total energy at each step. Default value is 0.
 * @param outF Pointer to an array to store the force at each step. Default value is 0.
 * @return The number of iterations performed.
 */
//int run( int nstepMax, double dt, double Fconv=1e-6, int ialg=0, double* outE, double* outF ){ 
    __attribute__((hot))  
    virtual int run( int nstepMax, double dt=-1, double Fconv=1e-6, int ialg=2, double* outE=0, double* outF=0, double* outV=0, double* outVF=0 ){ 
        //printf( "MolWorld_sp3::run(%i) \n", nstepMax );
        //printf( "MolWorld_sp3::run() nstepMax %i double dt %g Fconv %g ialg %g \n", nstepMax, dt, Fconv, ialg );
        //printf( "opt.damp_max %g opt.damping %g \n", opt.damp_max, opt.damping );
        double F2conv=Fconv*Fconv;
        double F2 = 1.0;
        double Etot=0;
        int itr=0;
        //if( (ialg!=0)&(!opt_initialized) ){ printf("ERROR ialg(%i)>0 but optimizer not initialized => call initOpt() first !"); exit(0); };
        if(dt>0){ opt.setTimeSteps(dt); }
        //if(ialg>0){ opt.cleanVel( ); }
        for(itr=0; itr<nstepMax; itr++ ){        
            //ff.clearForce();
            Etot = eval();
            switch(ialg){
                case  0: ffl.move_GD      (opt.dt);      break;
                case -1: opt.move_LeapFrog(opt.dt);      break;
                case  1: F2 = opt.move_MD (opt.dt,opt.damping); break;
                case  2: F2 = opt.move_FIRE();          break;
                case  3: F2 = opt.move_FIRE_smooth();   break;
            }
            opt_log.set(itr, opt.cos_vf, opt.f_len, opt.v_len, opt.dt, opt.damping );
            if(outE){ outE[itr]=Etot; }
            if(outF){ outF[itr]=F2;   }
            if( (trj_fname) && (itr%savePerNsteps==0) )[[unlikely]]{
                sprintf(tmpstr,"# %i E %g |F| %g", itr, Etot, sqrt(F2) );
                saveXYZ( trj_fname, tmpstr, false, "a", nPBC_save );
            }
            if(verbosity>1)[[unlikely]]{ printf("[%i] Etot %g[eV] |F| %g [eV/A] \n", itr, Etot, sqrt(F2) ); };
            if(F2<F2conv)[[unlikely]]{
                bConverged=true;
                if(verbosity>0)[[unlikely]]{ printf("MolWorld_sp3::run() Converged in %i iteration Etot %g[eV] |F| %g[eV/A] <(Fconv=%g) \n", itr, Etot, sqrt(F2), Fconv ); };
                if( trj_fname )[[unlikely]]{
                    sprintf(tmpstr,"# %i E %g |F| %g", itr, Etot, sqrt(F2) );
                    saveXYZ( trj_fname, tmpstr, false, "a", nPBC_save );
                }
                break;
            }
        }
        //printShortestBondLengths();
        return itr;
    }


/**
 * Pulls an atom towards a target position using a spring force.
 * 
 * @param ia The index of the atom to be pulled.
 * @param apos The array of atom positions.
 * @param fapos The array of final atom positions after applying the force.
 * @param K The spring constant (default value: -2.0).
 */
void pullAtom( int ia, Vec3d* apos, Vec3d* fapos, float K=-2.0 ){ 
    Vec3d f = getForceSpringRay( apos[ia], pick_hray, pick_ray0, K ); fapos[ia].add( f );
}

    virtual void MDloop( int nIter, double Ftol=-1 ){
        if(iParalel<-100){ iParalel=iParalel_default; };
        if(Ftol<0)Ftol=Ftol_default;
        ffu.bNonBonded     = bNonBonded;
        ffu.bNonBondNeighs = bNonBondNeighs;
        long T0 = getCPUticks();
        int nitr=0;
        if(bUFF){
            switch(iParalel){
                case  0: nitr=ffu.run    ( nIter, dt_default, Ftol, 1000.0 ); break;
                case  1: nitr=ffu.run_omp( nIter, dt_default, Ftol, 1000.0 ); break;
                default: [[unlikely]] {
                    printf( "ERROR: MolWorld_sp3::MDloop() iParalel(%i) not implemented (use 0=run_no_omp(), 1=run_omp()) \n", iParalel );
                    exit(0);
                }
            }
        }else{
            switch(iParalel){
                case  0: nitr=run_no_omp( nIter, dt_default, Ftol, 1000.0 ); break;
                case  1: nitr=run_omp   ( nIter, dt_default, Ftol, 1000.0 ); break;
                //case  0: nitr=run_no_omp( nIter, 0.02, Ftol, 1000.0 ); break;
                //case  1: nitr=run_omp   ( nIter, 0.02, Ftol, 1000.0 ); break;
                default: [[unlikely]] {
                    printf( "ERROR: MolWorld_sp3::MDloop() iParalel(%i) not implemented (use 0=run_no_omp(), 1=run_omp()) \n", iParalel );
                    exit(0);
                }
            }
        }

        { // Measure time 
            double ticks = (getCPUticks() - T0);
            double t = ticks*tick2second;
            double c_smooth = 0.1;
            
            double titer =  ( t*1e+6/nitr );
            if( (time_per_iter>(titer*2))||(time_per_iter<(titer*0.5)) ) time_per_iter=titer;
            time_per_iter = time_per_iter*(1-c_smooth) + titer*c_smooth;
            printf( "MolWorld_sp3::MDloop()  (bUFF=%i,iParalel=%i,bSurfAtoms=%i,bGridFF=%i,GridFF.mode=%i,GridFF.npbc=%i,npbc=%i,bPBC=%i,bNonBonded=%ibNonBondNeighs=%i,go.bExploring=%i,dt=%g,niter=%i) time=%g[ms/%i](%g[us/iter] tick2second=%g)\n", bUFF,iParalel,bSurfAtoms,bGridFF,(int)gridFF.mode,gridFF.npbc,npbc,bPBC,bNonBonded,bNonBondNeighs,go.bExploring,dt_default,nitr, t*1e+3,nitr, time_per_iter, tick2second );
        }

        //run( nIter );
        icurIter+=nitr;
        bChargeUpdated=false;
    }

double eval_no_omp(){
    double E=0;
    double F2max = ffl.FmaxNonBonded*ffl.FmaxNonBonded;
    ffl.bNonBonded=bNonBonded; ffl.setNonBondStrategy( bNonBondNeighs*2-1 );
    for(int i=0; i<ffl.natoms; i++){ ffl.fapos[i]=Vec3dZero; }
    for(int ia=0; ia<ffl.natoms; ia++){ 
        {                 ffl.fapos[ia           ] = Vec3dZero; } // atom pos force
        if(ia<ffl.nnode){ ffl.fapos[ia+ffl.natoms] = Vec3dZero; } // atom pi  force
        if(bMMFF)[[likely]]{
            if(ia<ffl.nnode){ E+=ffl.eval_atom(ia); }
        }
        //printf( "debug.1 E[%i]=%g\n", ia, E );
        // ----- Error is HERE
        if(bNonBonded){
            if(bNonBondNeighs)[[likely]]{
                if(bPBC)[[likely]]{ E+=ffl.evalLJQs_ng4_PBC_atom_omp( ia ); }
                else              { E+=ffl.evalLJQs_ng4_atom_omp    ( ia ); } 
            }else{
                if(bPBC)[[likely]]{ E+=ffl.evalLJQs_PBC_atom_omp( ia, F2max ); }
                else              { E+=ffl.evalLJQs_atom_omp    ( ia, F2max ); } 
            }
        }
        //printf( "debug.2 E[%i]=%g\n", ia, E );
        if(bSurfAtoms)[[likely]]{ 
            if(bGridFF)[[likely]]{  // with gridFF
                E += gridFF.addAtom( ffl.apos[ia], ffl.PLQd[ia], ffl.fapos[ia] );
            }else{ // Without gridFF (Direct pairwise atoms)
                //{ E+= nbmol .evalMorse   ( surf, false,                  gridFF.alphaMorse, gridFF.Rdamp );  }
                //{ E+= nbmol .evalMorsePBC    ( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );  }
                { E+= gridFF.evalMorsePBC_sym( ffl.apos[ia], ffl.REQs[ia],  ffl.fapos[ia] );   }
            }
        }
        //printf( "debug.3 E[%i]=%g\n", ia, E );
        if(bConstrZ){
            E+=springbound( ffl.apos[ia].z-ConstrZ_xmin, ConstrZ_l, ConstrZ_k, ffl.fapos[ia].z );
        }
        //printf( "debug.4 E[%i]=%g\n", ia, E );
    }
    // ---- assembling
    for(int ia=0; ia<ffl.natoms; ia++){
        ffl.assemble_atom( ia );
    } 
    if(bConstrains){
        E += constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec );
    }
    //printf( "debug.5 E=%g\n", E );
    if( go.bExploring){
        E += go.constrs.apply( ffl.apos, ffl.fapos, &(ffl.lvec) );
    }
    //printf( "debug.6 E=%g\n", E );
    if(bGroups){ groups.applyAllForces(0.0, 0.2*sin(nloop*0.02) ); }
    return E;
}


    __attribute__((hot))  
    int run_no_omp( int niter_max, double dt, double Fconv=1e-6, double Flim=1000, double damping=-1.0, double* outE=0, double* outF=0, double* outV=0, double* outVF=0 ){
        //printf( "MolWorld_sp3::run_no_omp() niter_max %i dt %g Fconv %g Flim %g damping %g out{E,vv,ff,vf}(%li,%li,%li,%li) \n", niter_max, dt, Fconv, Flim, damping, (long)outE, (long)outF, (long)outV, (long)outVF );
        printf( "MolWorld_sp3::run_no_omp() ffl.natoms=%i \n", ffl.natoms );
        nloop++;
        if(dt>0){ opt.setTimeSteps(dt); }else{ dt=opt.dt; }
        //if(verbosity>1)[[unlikely]]{ printf( "MolWorld_sp3::run_no_omp() niter_max %i dt %g Fconv %g Flim %g damping %g out{E,vv,ff,vf}(%li,%li,%li,%li) \n", niter_max, dt, Fconv, Flim, damping, (long)outE, (long)outF, (long)outV, (long)outVF ); }
        long T0 = getCPUticks();
        double E=0,F2=0,F2conv=Fconv*Fconv;
        double ff=0,vv=0,vf=0;
        int itr=0,niter=niter_max;
        bConverged = false;
        bool bFIRE = true;

        double cdamp = ffl.colDamp.update( dt );  if(cdamp>1)cdamp=1;
        // if(damping>0){ cdamp = 1-damping; if(cdamp<0)cdamp=0;}
        double F2max = ffl.FmaxNonBonded*ffl.FmaxNonBonded;

        ffl.bNonBonded=bNonBonded; ffl.setNonBondStrategy( bNonBondNeighs*2-1 );

        //printf( "MolWorld_sp3::run_no_omp(itr=%i/%i) gridFF.mode=%i bGridFF=%i bSurfAtoms=%i   gridFF.Bspline_PLQ=%li  FFPaul_d=%li FFLond_d=%li FFelec_d=%li \n", itr,niter_max, gridFF.mode, bGridFF, bSurfAtoms, (long)gridFF.Bspline_PLQ, (long)gridFF.FFPaul_d, (long)gridFF.FFLond_d, (long)gridFF.FFelec_d );

        //printf( "MolWorld_sp3::run_no_omp() bNonBonded=%i bNonBondNeighs=%i bSubtractBondNonBond=%i bSubtractAngleNonBond=%i bClampNonBonded=%i\n", bNonBonded, bNonBondNeighs, ffl.bSubtractBondNonBond, ffl.bSubtractAngleNonBond, ffl.bClampNonBonded );

        //if(bToCOG){ printf("bToCOG=%i \n", bToCOG ); center(true); }
        //if(bToCOG){ Vec3d cog=average( ffl.natoms, ffl.apos );  move( ffl.natoms, ffl.apos, cog*-1.0 ); }
        if(bCheckStuck){ checkStuck( RStuck ); }
        //if(verbosity>0)printf( "MolWorld_sp3::run_no_omp(niter=%i,bColB=%i,bColNB=%i) dt %g damping %g colB %g colNB %g \n", niter_max, ffl.bCollisionDamping, ffl.bCollisionDampingNonBond, dt, 1-cdamp, ffl.col_damp*dt, ffl.col_damp_NB*dt );
        //if(verbosity>1)[[unlikely]]{printf( "MolWorld_sp3::run_no_omp(niter=%i,bCol(B=%i,A=%i,NB=%i)) dt %g damp(cM=%g,cB=%g,cA=%g,cNB=%g)\n", niter, ffl.colDamp.bBond, ffl.colDamp.bAng, ffl.colDamp.bNonB, dt, 1-cdamp, ffl.colDamp.cdampB*dt, ffl.colDamp.cdampAng*dt, ffl.colDamp.cdampNB*dt ); }
        for(itr=0; itr<niter; itr++){
            //double ff=0,vv=0,vf=0;
            if(bGopt){
                if( bGopt         ) go.update();
                if( go.bExploring ) bConverged = false;
            }

            if(bGroups){
                groups.evalAllPoses();
            }

            Vec3d cvf_bak = ffl.cvf;
            E=0; ffl.cvf = Vec3dZero;
            //------ eval forces
            //long t1 = getCPUticks();

            for(int i=0; i<ffl.natoms; i++){ ffl.fapos[i]=Vec3dZero; }
            for(int ia=0; ia<ffl.natoms; ia++){ 
                {                 ffl.fapos[ia           ] = Vec3dZero; } // atom pos force
                if(ia<ffl.nnode){ ffl.fapos[ia+ffl.natoms] = Vec3dZero; } // atom pi  force
                if(bMMFF)[[likely]]{
                    if(ia<ffl.nnode){ E+=ffl.eval_atom(ia); }
                }
                // ----- Error is HERE
                if(bNonBonded){
                    if(bNonBondNeighs)[[likely]]{
                        if(bPBC)[[likely]]{ E+=ffl.evalLJQs_ng4_PBC_atom_omp( ia ); }
                        else              { E+=ffl.evalLJQs_ng4_atom_omp    ( ia ); } 
                    }else{
                        if(bPBC)[[likely]]{ E+=ffl.evalLJQs_PBC_atom_omp( ia, F2max ); }
                        else              { E+=ffl.evalLJQs_atom_omp    ( ia, F2max ); } 
                    }
                }
                if(bSurfAtoms)[[likely]]{ 
                    if(bGridFF)[[likely]]{  // with gridFF
                        gridFF.addAtom( ffl.apos[ia], ffl.PLQd[ia], ffl.fapos[ia] );
                        //Vec3d fi=Vec3dZero;
                        //gridFF.addAtom( ffl.apos[ia], ffl.PLQd[ia], fi );
                        //printf("MolWorld_sp3::ffu.atomForceFunc(ia=%i,gridFF.mode=%i) p(%g,%g,%g) fi(%g,%g,%g) PLQ(%g,%g,%g)\n", ia, (int)gridFF.mode, ffl.apos[ia].x, ffl.apos[ia].y, ffl.apos[ia].z, fi.x,fi.y,fi.z, ffl.PLQd[ia].x, ffl.PLQd[ia].y, ffl.PLQd[ia].z );
                        //ffl.fapos[ia].add( fi );
                    }else{ // Without gridFF (Direct pairwise atoms)
                        //{ E+= nbmol .evalMorse   ( surf, false,                  gridFF.alphaMorse, gridFF.Rdamp );  }
                        //{ E+= nbmol .evalMorsePBC    ( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );  }
                        { E+= gridFF.evalMorsePBC_sym( ffl.apos[ia], ffl.REQs[ia],  ffl.fapos[ia] );   }
                    }
                }
                // if   (bGridFF){ 
                //     if  (bTricubic){ E+= gridFF.addForce_Tricubic( ffl.apos[ia], ffl.PLQd[ia], ffl.fapos[ia], true  ); }
                //     else           { E+= gridFF.addForce         ( ffl.apos[ia], ffl.PLQs[ia], ffl.fapos[ia], true  ); }
                // }  // GridFF
                //if( ffl.colDamp.bNonB ){ ffl.evalCollisionDamp_atom_omp( ia, ffl.colDamp.cdampNB, ffl.colDamp.dRcut1, ffl.colDamp.dRcut2 ); }
                //if( ffl.bCollisionDampingNonBond ){ ffl.evalCollisionDamp_atom_omp( ia, ffl.col_damp_NB, ffl.col_damp_dRcut1, ffl.col_damp_dRcut2 ); }
                if(bConstrZ){
                    springbound( ffl.apos[ia].z-ConstrZ_xmin, ConstrZ_l, ConstrZ_k, ffl.fapos[ia].z );
                }
                //if(bGroups){ groups.forceAtom(ia); }
                if(ipicked==ia)[[unlikely]]{ 
                    const Vec3d f = getForceSpringRay( ffl.apos[ia], pick_hray, pick_ray0,  Kpick ); 
                    ffl.fapos[ia].add( f );
                }
            }
            // ---- assembling
            for(int ia=0; ia<ffl.natoms; ia++){
                ffl.assemble_atom( ia );
            } 

            //double t_eval = (getCPUticks()-t1);
            //printf( "MolWorld_sp3::run_no_omp() (bPBC=%i,bGridFF=%i,bNonBondNeighs=%i,|Fmax|=%g,dt=%g,niter=%i) %g[tick]\n", bPBC,bGridFF,bNonBondNeighs,sqrt(F2max),opt.dt,niter, t_eval );
            //for(int i=0; i<ffl.natoms; i++){ printf( "ffl.fapos[%i] (%g,%g,%g)\n", i, ffl.fapos[i].x, ffl.fapos[i].y, ffl.fapos[i].z ); }
            if(bConstrains){
                //printf( "run_no_omp() constrs[%li].apply() bThermalSampling=%i \n", constrs.bonds.size(), bThermalSampling );
                //ffl.cleanForce();
                //printf( "run_omp() constrs[%i].apply()\n", constrs.bonds.size() );
                E+= constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec );
            }
            if( go.bExploring){
                ffl.Etot += go.constrs.apply( ffl.apos, ffl.fapos, &(ffl.lvec) );
            }
            if(bGroups){ groups.applyAllForces(0.0, 0.2*sin(nloop*0.02) ); }
            if(bFIRE){
                for(int i=0; i<opt.n; i++){
                    double v=opt.vel  [i];
                    double f=opt.force[i];
                    //vv+=v*v; ff+=f*f; vf+=v*f;
                    ffl.cvf.x += v*f;
                    ffl.cvf.y += v*v;
                    ffl.cvf.z += f*f;
                }
                opt.vf=ffl.cvf.x;
                opt.vv=ffl.cvf.y; 
                opt.ff=ffl.cvf.z; 
                opt.FIRE_update_params();
            }else{
                // ----- Dynamics
                if(ffl.colDamp.medium<0){ cdamp=ffl.colDamp.tryAccel(); }
                //if(ffl.colDamp.medium<0){  if( ffl.colDamp.canAccelerate() ){ cdamp=1-ffl.colDamp.medium;  }else{ cdamp=1+ffl.colDamp.medium*0.2; } }
                //if(ffl.colDamp.medium<0){  if( ffl.colDamp.canAccelerate() ){ cdamp=1-ffl.colDamp.medium; }else{ cdamp=1+ffl.colDamp.medium; } }
            }
            for(int i=0; i<ffl.nvecs; i++){
                if( go.bExploring ){
                    //if(i==0)printf( "run_no_omp() move_atom_Langevin() gamma_damp%g T_target=%g \n", go.gamma_damp, go.T_target );
                    ffl.move_atom_Langevin( i, dt, 10000.0, go.gamma_damp, go.T_target );
                }else if(bFIRE){
                    ffl.move_atom_FIRE( i, opt.dt, 10000.0, opt.cv, opt.renorm_vf*opt.cf );
                }else{
                    ffl.cvf.add( ffl.move_atom_MD( i, dt, Flim, cdamp ) );
                }
                //ffl.cvf.add( ffl.move_atom_MD( i, dt, Flim, 1.0 ) );
                //ffl.cvf.add( ffl.move_atom_MD( i, dt, Flim, 1.01 ) );
                //ffl.move_atom_MD( i, 0.05, Flim, 0.9 );
                //ffl.move_atom_MD( i, 0.05, 1000.0, 0.9 );
                //ffl.move_atom_FIRE( i, dt, 10000.0, 0.9, 0 ); // Equivalent to MDdamp
            }
            if( (!bFIRE) && (!go.bExploring) ){
                //printf( "Acceleration: cdamp=%g(%g) nstepOK=%i \n", cdamp_, cdamp, ffl.colDamp.nstepOK );
                //printf( "Kinetic energy decrease: factor=%g v2_new=%g v2_old=%g \n", ffl.cvf.y/cvf_bak.y,  ffl.cvf.y, cvf_bak.y );
                double cos_vf = ffl.colDamp.update_acceleration( ffl.cvf );
                if(ffl.cvf.x<0)[[unlikely]]{ ffl.cleanVelocity(); }else{
                    //double renorm_vf  = sqrt( ffl.cvf.y / ( ffl.cvf.z  + 1e-32 ) );
                    //for(int i=0; i<ffl.nvecs; i++){ ffl.vapos[i] = ffl.vapos[i]*(1-cdamp) + ffl.fapos[i]*( renorm_vf * cdamp ); }
                }
            }
            Etot=E;
            if(outE )outE [itr]=Etot;
            if(outF )outF [itr]=sqrt(ffl.cvf.z);
            if(outV )outV [itr]=sqrt(ffl.cvf.y);
            if(outVF)outVF[itr]=ffl.cvf.x/sqrt(ffl.cvf.z*ffl.cvf.y + 1e-32);
            //if( isnan(Etot) || isnan(ffl.cvf.z) || isnan(ffl.cvf.y) )[[unlikely]]{  printf( "MolWorld_sp3::run_no_omp(itr=%i) ERROR NaNs : E=%f cvf(%g,%g,%g)\n", itr, E, ffl.cvf.x, sqrt(ffl.cvf.y), sqrt(ffl.cvf.z) ); return -itr; }
            if( (trj_fname) && ( (itr%savePerNsteps==0) ||(niter==0) ) )[[unlikely]]{
                sprintf(tmpstr,"# %i E %g |F| %g", itr, Etot, sqrt(ffl.cvf.z) );
                //printf( "run_no_omp::save() %s \n", tmpstr );
                saveXYZ( trj_fname, tmpstr, false, "a", nPBC_save );
            }
            if(ffl.cvf.z<F2conv)[[unlikely]]{ 
                //niter=0; 
                bConverged=true;
                double t = (getCPUticks() - T0)*tick2second;
                if(verbosity>1)[[unlikely]]{printf( "MolWorld_sp3::run_no_omp() CONVERGED in %i/%i nsteps E=%g |F|=%g time=%g[ms/%i](%g[us/iter])\n", itr,niter_max, E, sqrt(ffl.cvf.z), itr, t*1e+3, t*1e+6/itr ); }
                if( bGopt && !go.bExploring )[[unlikely]]{
                    gopt_ifound++;
                    sprintf(tmpstr,"# %i E %g |F| %g istep=%i ", gopt_ifound, Etot, sqrt(ffl.cvf.z), go.istep );
                    printf( "run_no_omp::save() %s \n", tmpstr );
                    saveXYZ( "gopt.xyz", tmpstr, false, "a", nPBC_save );
                    go.startExploring();
                    //go.apply_kick( ffl.natoms, ffl.apos, ffl.vapos );
                    bConverged=false;
                }
                break;
            }else{
                //printf( "nStuck %i \n", nStuck );
                if( bCheckStuck )[[unlikely]] { handleStuckAtom(itr, ffl.cvf ); }
                if(verbosity>3)  [[unlikely]] { printf( "MolWorld_sp3::run_no_omp(itr=%i/%i) E=%g |F|=%g |v|=%g cos(v,f)=%g dt=%g cdamp=%g\n", itr,niter_max, E, sqrt(ffl.cvf.z), sqrt(ffl.cvf.y), ffl.cvf.x/sqrt(ffl.cvf.z*ffl.cvf.y+1e-32), dt, cdamp ); }
            }
        }
        double ticks = (getCPUticks() - T0);
        double t = ticks*tick2second;
        if( (itr>=(niter_max-1)) && (verbosity>1) ) [[unlikely]] {
            double c_smooth = 0.1;
            time_per_iter = time_per_iter*(1-c_smooth) + ( t*1e+6/itr )*c_smooth;
            //printf( "MolWorld_sp3::run_no_omp() NOT CONVERGED in %i/%i dt=%g E=%g |F|=%g time=%g[ms/%i](%g[us/iter])\n", itr,niter_max, opt.dt, E,  sqrt(ffl.cvf.z), t*1e+3,itr,t*1e+6/itr ); 
            printf( "MolWorld_sp3::run_no_omp() NOT CONVERGED (bPBC=%i,bGridFF=%i,bNonBondNeighs=%i,go.bExploring=%i,|Fmax|=%g,dt=%g,niter=%i) time=%g[ms/%i](%g[us/iter]) | tick2second=%g ticks=%g \n", bPBC,bGridFF,bNonBondNeighs,go.bExploring,sqrt(F2max),opt.dt,niter, t*1e+3,itr, time_per_iter,  tick2second,  ticks );
        }
        return itr;
    }
  
/**
 * Runs the simulation using OpenMP parallelization.
 *
 * @param niter_max The maximum number of iterations.
 * @param dt The time step.
 * @param Fconv The convergence threshold for the force.
 * @param Flim The force limit.
 * @param timeLimit The time limit for the optimization.
 * @param outE Pointer to an array to store the total energy at each iteration (optional).
 * @param outF Pointer to an array to store the squared force at each iteration (optional).
 * @return The number of iterations performed.
 */
    __attribute__((hot))  
    int run_omp( int niter_max, double dt, double Fconv=1e-6, double Flim=1000, double timeLimit=0.02, double* outE=0, double* outF=0, double* outV=0, double* outVF=0 ){
        nloop++;
        if(dt>0){ opt.setTimeSteps(dt); }else{ dt=opt.dt; }
        //int ncpu = omp_get_num_threads(); printf("MolWorld_sp3::run_omp() ncpu=%i \n", ncpu );
        //printf( "run_omp() niter_max %i dt %g Fconv %g Flim %g timeLimit %g outE %li outF %li \n", niter_max, dt, Fconv, Flim, timeLimit, (long)outE, (long)outF );
        long T0 = getCPUticks();
        double E=0,F2=0,F2conv=Fconv*Fconv;
        double ff=0,vv=0,vf=0;
        int itr=0,niter=niter_max;
        bConverged = false;
        bool bFIRE = true;

        double cdamp = ffl.colDamp.update( dt );  if(cdamp>1)cdamp=1;
        // if(damping>0){ cdamp = 1-damping; if(cdamp<0)cdamp=0;}
        double F2max = ffl.FmaxNonBonded*ffl.FmaxNonBonded;

        ffl.bNonBonded=bNonBonded; ffl.setNonBondStrategy( bNonBondNeighs*2-1 );
        //printf( "MolWorld_sp3::run_no_omp() bNonBonded=%i bNonBondNeighs=%i bSubtractBondNonBond=%i bSubtractAngleNonBond=%i bClampNonBonded=%i\n", bNonBonded, bNonBondNeighs, ffl.bSubtractBondNonBond, ffl.bSubtractAngleNonBond, ffl.bClampNonBonded );


        //if( bGridDouble ){ printf( "run_omp() bGridDouble %i @ffl.PLQd=%li @FFPaul_d=%li @FFLond_d=%li @FFPaul_d=%li \n", bGridDouble, (long)ffl.PLQd, (long)gridFF.FFPaul_d, (long)gridFF.FFLond_d, (long)gridFF.FFPaul_d );    }

        // if(bToCOG && (!bGridFF) ){ 
        //     Vec3d cog=average( ffl.natoms, ffl.apos );  
        //     move( ffl.natoms, ffl.apos, cog*-1.0 ); 
        // }
        if(bCheckStuck){ checkStuck( RStuck ); }
        //#pragma omp parallel shared(E,F2,ff,vv,vf,ffl) private(itr)
        #pragma omp parallel shared(niter,itr,E,F2,ff,vv,vf,ffl,T0,bConstrains,bConverged)
        while(itr<niter){
            if(itr<niter){
            //#pragma omp barrier
            #pragma omp single
            {E=0;F2=0;ff=0;vv=0;vf=0;
                if( bGopt         ) go.update();
                if( go.bExploring ) bConverged = false;
            }
            //------ eval forces
            //#pragma omp barrier
            #pragma omp for reduction(+:E)
            for(int ia=0; ia<ffl.natoms; ia++){ 
                {                 ffl.fapos[ia           ] = Vec3dZero; } // atom pos force
                if(ia<ffl.nnode){ ffl.fapos[ia+ffl.natoms] = Vec3dZero; } // atom pi  force
                //if(verbosity>3)
                //printf( "atom[%i]@cpu[%i/%i]\n", ia, omp_get_thread_num(), omp_get_num_threads()  );
                if(ia<ffl.nnode){ E+=ffl.eval_atom(ia); }
                //if(ia<ffl.nnode){ E+=ffl.eval_atom_opt(ia); }

                // ----- Error is HERE
                if(bNonBonded){
                    if(bNonBondNeighs)[[likely]]{
                        if(bPBC)[[likely]]{ E+=ffl.evalLJQs_ng4_PBC_atom_omp( ia ); }
                        else              { E+=ffl.evalLJQs_ng4_atom_omp    ( ia ); } 
                    }else{
                        if(bPBC)[[likely]]{ E+=ffl.evalLJQs_PBC_atom_omp( ia, F2max ); }
                        else              { E+=ffl.evalLJQs_atom_omp    ( ia, F2max ); } 
                    }
                }

                if(bSurfAtoms)[[likely]]{ 
                    if(bGridFF)[[likely]]{  // with gridFF
                        gridFF.addAtom( ffl.apos[ia], ffl.PLQd[ia], ffl.fapos[ia] );
                    }else{ // Without gridFF (Direct pairwise atoms)
                        //{ E+= nbmol .evalMorse   ( surf, false,                  gridFF.alphaMorse, gridFF.Rdamp );  }
                        //{ E+= nbmol .evalMorsePBC    ( surf, gridFF.grid.cell, nPBC, gridFF.alphaMorse, gridFF.Rdamp );  }
                        { E+= gridFF.evalMorsePBC_sym( ffl.apos[ia], ffl.REQs[ia],  ffl.fapos[ia] );   }
                    }
                }

                //if     (bGridFF){ E+= gridFF.addMorseQH_PBC_omp( ffl.apos[ia], ffl.REQs[ia], ffl.fapos[ia] ); }  // NBFF
                if(bConstrZ){
                    springbound( ffl.apos[ia].z-ConstrZ_xmin, ConstrZ_l, ConstrZ_k, ffl.fapos[ia].z );
                }
                if(ipicked==ia)[[unlikely]]{ 
                    const Vec3d f = getForceSpringRay( ffl.apos[ia], pick_hray, pick_ray0,  Kpick ); 
                    ffl.fapos[ia].add( f );
                }
            }
            // if(ffl.bTorsion){
            //     #pragma omp for reduction(+:E)
            //     for(int it=0; it<ffl.ntors; it++){ 
            //         E+=ffl.eval_torsion(it); 
            //     }
            // }
            #pragma omp single
            {
                if(bConstrains  ){ E+=constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec ); }
                if(go.bExploring){ go.constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec ); }
            }
            // ---- assemble (we need to wait when all atoms are evaluated)
            //#pragma omp barrier
            #pragma omp for
            for(int ia=0; ia<ffl.natoms; ia++){
                ffl.assemble_atom( ia );
            }
            // #pragma omp barrier
            // #pragma omp for reduction(+:F2)
            // for(int i=0; i<ffl.nvecs; i++){
            //     F2 += ffl.move_atom_MD     ( i, dt, Flim, 0.95 );
            //     //F2 += ffl.move_atom_kvaziFIRE( i, dt, Flim );
            // }
            //#pragma omp barrier
            
            //#pragma omp barrier
            // --- FIRE pre
                #pragma omp for reduction(+:vf,vv,ff)
                for(int i=0; i<opt.n; i++){
                    double v=opt.vel  [i];
                    double f=opt.force[i];
                    vv+=v*v; ff+=f*f; vf+=v*f;
                }
                #pragma omp single
                { opt.vv=vv; opt.ff=ff; opt.vf=vf; F2=ff; opt.FIRE_update_params(); }

                // ------ move
                #pragma omp for
                for(int i=0; i<ffl.nvecs; i++){
                    if( go.bExploring ){
                        ffl.move_atom_Langevin( i, dt, 10000.0, go.gamma_damp, go.T_target ); 
                    }else{
                        //ffl.move_atom_MD( i, opt.dt, Flim, 0.9 );
                        //ffl.move_atom_MD( i, 0.05, Flim, 0.9 );
                        //ffl.move_atom_MD( i, 0.05, 1000.0, 0.9 );
                        ffl.move_atom_FIRE( i, opt.dt, 10000.0, opt.cv, opt.renorm_vf*opt.cf );
                        //ffl.move_atom_FIRE( i, dt, 10000.0, 0.9, 0 ); // Equivalent to MDdamp
                    }
            } 
            
            //#pragma omp barrier
            #pragma omp single
            { 
                Etot=E;
                itr++; 
                // if(timeLimit>0)[[unlikely]]{
                //     double t = (getCPUticks() - T0)*tick2second;
                //     if(t>0.02){ 
                //         niter=0; 
                //         if(verbosity>1) [[unlikely]] { printf( "run_omp() ended due to time limit after %i nsteps ( %6.3f [s]) \n", itr, t ); }
                //     }
                // }
                if(F2<F2conv)[[unlikely]]{ 
                    niter=0; 
                    bConverged = true;
                    double t = (getCPUticks() - T0)*tick2second;
                    if(verbosity>1) [[unlikely]] { printf( "MolWorld_sp3::run_omp() CONVERGED in %i/%i nsteps E=%g |F|=%g time= %g [ms]( %g [us/%i iter])\n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr ); }

                    //printf( "MolWorld_sp3::run_omp() CONVERGED in %i/%i nsteps E=%g |F|=%g time= %g [ms]( %g [us/%i iter])\n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr );
                    if(bGopt  && (!go.bExploring) ){
                        gopt_ifound++;
                        sprintf(tmpstr,"# %i E %g |F| %g istep=%i", gopt_ifound, Etot, sqrt(F2), go.istep );
                        printf( "MolWorld_sp3::run_omp().save %s \n", tmpstr );
                        saveXYZ( "gopt.xyz", tmpstr, false, "a", nPBC_save );
                        go.startExploring();
                        //go.apply_kick( ffl.natoms, ffl.apos, ffl.vapos );
                        bConverged=false;
                        // if(bToCOG && bGridFF ){ 
                        //     Vec3d cog=average( ffl.natoms, ffl.apos );   cog.z=0;  
                        //     move( ffl.natoms, ffl.apos, cog*-1.0 ); 
                        // }
                    }
                }else{
                    if( bCheckStuck  )[[unlikely]]{ handleStuckAtom(itr, Vec3d{opt.vf,opt.vv,opt.ff} ); }
                }
                //printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(F2), omp_get_num_threads() ); 
                //{printf( "step[%i] dt %g(%g) cv %g cf %g cos_vf %g \n", itr, opt.dt, opt.dt_min, opt.cv, opt.cf, opt.cos_vf );}
                //if(verbosity>2){printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(F2), omp_get_num_threads() );}
            }
            } // if(itr<niter)
        } // while(itr<niter)
           
        {
        double t = (getCPUticks() - T0)*tick2second;
        if( (itr>=niter_max)&&(verbosity>1)) [[unlikely]] {printf( "MolWorld_sp3::run_omp() NOT CONVERGED in %i/%i dt=%g E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, opt.dt, E,  sqrt(F2), t*1e+3, t*1e+6/itr, itr ); }
        }
        return itr;
    }

void scan_rigid( int nconf, Vec3d* poss, Mat3d* rots, double* Es, Vec3d* aforces, Vec3d* aposs, bool omp ){
    printf("MolWorld_sp3::scan_rigid(nconf=%i,omp=%i) @poss=%li @rots=%li @Es=%li @aforces=%li @aposs=%li \n", nconf, omp, (long)poss, (long)rots, (long)Es, (long)aforces, (long)aposs);
    printf("MolWorld_sp3::scan_rigid() bNonBonded=%i bNonBondNeighs=%i bPBC=%i bSurfAtoms=%i bGridFF=%i gridFF.mode=%i \n", bNonBonded, bNonBondNeighs, bPBC, bSurfAtoms, bGridFF, gridFF.mode );

    for(int ia=0; ia<ffl.natoms; ia++){  printf( "MolWorld_sp3::scan_rigid()[ia=%i] pos(%8.4f,%8.4f,%8.4f) REQ(%8.4f,%16.8f,%8.4f,%8.4f) PLQd(%16.8f,%16.8f,%16.8f,%16.8f) \n", ia, ffl.apos[ia].x, ffl.apos[ia].y, ffl.apos[ia].z, ffl.REQs[ia].x, ffl.REQs[ia].y, ffl.REQs[ia].z, ffl.REQs[ia].w, ffl.PLQd[ia].x, ffl.PLQd[ia].y, ffl.PLQd[ia].z, ffl.PLQd[ia].w );  }

    Atoms atoms;
    atoms.copyOf( ffl );
    std::vector<Vec3d> pipos(ffl.nnode); for(int i=0; i<ffl.nnode; i++){ pipos[i]=ffl.pipos[i]; }
    for(int i=0; i<nconf; i++){
        Vec3d pos; if(poss){ pos=poss[i]; }else{ pos=Vec3dZero; }
        Mat3d rot; if(rots){ rot=rots[i]; }else{ rot=Mat3dIdentity; }
        ffl.setFromRef( atoms.apos, pipos.data(), pos, rot );
        double E = eval_no_omp();
        //printf( "scan_rigid[%i] E=%g \n", i, E );
        if(Es){ Es[i]=E; }
        if(aforces){ ffl.copyForcesTo( aforces + i*ffl.natoms ); }
        if(aposs  ){ ffl.copyPosTo   ( aposs   + i*ffl.natoms ); }
    }
}

void scan_relaxed( int nconf, Vec3d* poss, Mat3d* rots, double* Es, Vec3d* aforces, Vec3d* aposs, bool omp, int niter_max, double dt, double Fconv=1e-6, double Flim=1000 ){
    printf("MolWorld_sp3::scan_relaxed(nconf=%i,omp=%i) @poss=%li @rots=%li @Es=%li @aforces=%li @aposs=%li \n", nconf, omp, (long)poss, (long)rots, (long)Es, (long)aforces, (long)aposs);
    Atoms atoms;
    atoms.copyOf( ffl );
    std::vector<Vec3d> pipos(ffl.nnode); for(int i=0; i<ffl.nnode; i++){ pipos[i]=ffl.pipos[i]; }
    for(int i=0; i<nconf; i++){
        Vec3d pos; if(poss){ pos=poss[i]; }else{ pos=Vec3dZero; }
        Mat3d rot; if(rots){ rot=rots[i]; }else{ rot=Mat3dIdentity; }
        ffl.setFromRef( atoms.apos, pipos.data(), pos, rot );
        double E = eval_no_omp();
        int niterdone = run_no_omp( niter_max, dt, Fconv, Flim, 1000.0);
        if(Es){ Es[i]=E; }
        if(aforces){ ffl.copyForcesTo( aforces + i*ffl.natoms ); }
        if(aposs  ){ ffl.copyPosTo   ( aposs   + i*ffl.natoms ); }
    }
}

// void makeGridFF( bool bSaveDebugXSFs=false, Vec3i nPBC={1,1,0} ) {
//     gridFF.bindSystem(surf.natoms, 0, surf.apos, surf.REQs );
//     //if(verbosity>1)
//     gridFF.grid.printCell();
//     gridFF.allocateFFs();
//     //double x0= ( gridFF.grid.cell.a.x + gridFF.grid.cell.b.x )*-0.5;
//     //double y0= ( gridFF.grid.cell.a.y + gridFF.grid.cell.b.y )*-0.5;
//     //gridFF.grid.pos0 = Vec3d{ x0,y0,-8.0};
//     //gridFF.shift   = Vec3d{0.0,0.0,-8.0};
//     //gridFF.tryLoad( "data/FFelec.bin", "data/FFPaul.bin", "data/FFLond.bin", true, {0,0,0} );
//     //gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false, {1,1,1} );
//     gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false, nPBC, bSaveDebugXSFs );
// }


//inline int pickParticle( const Vec3d& ray0, const Vec3d& hRay, double R, int n, Vec3d * ps, bool* ignore=0 ){
//int pickParticle( Vec3d ray0, Vec3d hray, double R=0.5 ){ return pickParticle( ray0, hray, R, ff.natoms, ff.apos ); }

// int toXYZ(const char* comment="#comment", bool bNodeOnly=false){
//     if(xyz_file==0){ printf("ERROR no xyz file is open \n"); return -1; }
//     params.writeXYZ( xyz_file, (bNodeOnly ? ffl.nnode : ffl.natoms) , nbmol.atypes, nbmol.apos, comment );
//     return 0;
// }

/**
 * Prepares the MolecularWorld object to be written to an XYZ file by function writeXYZ.
 *
 * @param comment The comment to be included in the XYZ file (default: "#comment").
 * @param bNodeOnly Flag indicating whether to convert only the nodes or the entire structure (default: false).
 * @param file The file pointer to write the XYZ data to (default: 0).
 * @param bPi Flag indicating whether to include pi electrons in the XYZ file (default: false).
 * @param just_Element Flag indicating whether to include only the element symbols in the XYZ file (default: true).
 * @return 0 if successful, -1 if an error occurs.
 */
int toXYZ(const char* comment="#comment", bool bNodeOnly=false, FILE* file=0, bool bPi=false, bool just_Element=true ){
    if(file==0){ file=xyz_file; };
    if(file==0){ printf("ERROR no xyz file is open \n"); return -1; }
    int n=ffl.natoms; if(bNodeOnly){ n=ffl.nnode; }
    int npi=0; if(bPi)npi=ffl.nnode;
    params.writeXYZ( file, n, nbmol.atypes, nbmol.apos, comment, 0,just_Element, npi );
    return 0;
}

/**
 * Saves the molecular structure in XYZ format to a file.
 *
 * @param fname The name of the file to save the structure to.
 * @param comment The comment to include in the XYZ file. Default is "#comment".
 * @param bNodeOnly Flag indicating whether to save only the nodes or all atoms. Default is false.
 * @param mode The file mode to open the file with. Default is "w".
 * @param nPBC The number of periodic boundary conditions in each direction. Default is {1, 1, 1}.
 * @return The number of atoms saved to the file.
 */
int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false, const char* mode="w", Vec3i nPBC=Vec3i{1,1,1} ){ 
    char str_tmp[1024];
    if(bPBC){ sprintf( str_tmp, "lvs %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f %s", ffl.lvec.a.x, ffl.lvec.a.y, ffl.lvec.a.z, ffl.lvec.b.x, ffl.lvec.b.y, ffl.lvec.b.z, ffl.lvec.c.x, ffl.lvec.c.y, ffl.lvec.c.z, comment ); }
    else    { sprintf( str_tmp, "%s", comment ); }
    return params.saveXYZ( fname, (bNodeOnly ? ffl.nnode : ffl.natoms) , nbmol.atypes, nbmol.apos, str_tmp, nbmol.REQs, mode, true, nPBC, ffl.lvec ); 
}
    //int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false){ return params.saveXYZ( fname, (bNodeOnly ? ff.nnode : ff.natoms) , ff.atype, ff.apos, comment, nbmol.REQs ); }
    //int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false){ return params.saveXYZ( fname, (bNodeOnly ? ff.nnode : ff.natoms) , nbmol.atypes, nbmol.apos, comment, nbmol.REQs ); }

    // ========= Analysis

    double findHbonds_PBC( int ia, double Rcut, double Hcut, double cosMin, Vec3d dir, std::vector<Vec3i>& out ){
        //printf( "MolWorld_sp3::findHbonds_PBC(%i) Rcut=%g Hcut=%g\n", ia, Rcut, Hcut );
        const double R2cut = Rcut*Rcut;
        const NBFF&  ff   = ffl;
        const Vec3d  pi   = ff.apos     [ia];
        const Quat4d REQi = ff.REQs     [ia];
        if( REQi.w < 0 ) return 0;
        const Quat4i ng   = ff.neighs   [ia];
        const Quat4i ngC  = ff.neighCell[ia];
        double E=0,fx=0,fy=0,fz=0;
        for (int j=0; j<ff.natoms; j++){ 
            if(ia==j)continue;
            const Quat4d& REQj  = ff.REQs[j];
            if( REQj.w > 0 ) continue;
            const Quat4d  REQij = _mixREQ(REQi,REQj); 
            const Vec3d dp     = ff.apos[j]-pi;
            Vec3d fij          = Vec3dZero;
            const bool bBonded = ((j==ng.x)||(j==ng.y)||(j==ng.z)||(j==ng.w));
            for(int ipbc=0; ipbc<npbc; ipbc++){
                //printf( "[ia=%i,j=%i,ipbc=%i]\n", ia, j, ipbc );
                // --- We calculate non-bonding interaction every time (most atom pairs are not bonded)
                const Vec3d dpc = dp + ff.shifts[ipbc];    //   dp = pj - pi + pbc_shift = (pj + pbc_shift) - pi 
                //double eij      = getLJQH( dpc, fij, REQij, R2damp );
                // --- If atoms are bonded we don't use the computed non-bonding interaction energy and force
                double r2 = dpc.norm2();
                if( (r2>R2cut) || (REQij.w>-Hcut) )[[likely]] continue;
                //if( (r2>R2cut) )[[likely]] continue;
                if(bBonded) [[unlikely]]  { 
                    if(   ((j==ng.x)&&(ipbc==ngC.x))
                        ||((j==ng.y)&&(ipbc==ngC.y))
                        ||((j==ng.z)&&(ipbc==ngC.z))
                        ||((j==ng.w)&&(ipbc==ngC.w))
                    ) [[unlikely]]  { 
                        continue;
                    }
                }
                double c = dir.dot( dpc )/sqrt(r2);
                //printf( "[ia=%i,j=%i,ipbc=%i] r=%g/Rcut(%g) cos=%g/cosMin(%g)\n", ia, j, ipbc,  sqrt(r2), Rcut, c, cosMin );
                if( c<cosMin ) continue;
                //printf( "[ia=%i,j=%i,ipbc=%i] found  r=%g/Rcut(%g) \n", ia, j, ipbc,  sqrt(r2), Rcut );
                out.push_back( Vec3i{ ia, j, ipbc } );
            }
        }
        return E;
    }
    double findHbonds_PBC( double Rcut, double Hcut, double angMax, std::vector<Vec3i>* out =0 ){
        //printf( "MolWorld_sp3::findHbonds_PBC()\n" );
        double cosMin = cos(angMax);
        if(out==0){ out = &Hbonds; }
        double E=0;
        const NBFF&  ff   = ffl;
        for(int ia=0; ia<ffl.natoms; ia++){ 

            if( ff.REQs[ia].w < Hcut ) continue;

            // --- find vector of hydrogen bond to base
            const Quat4i ng  = ff.neighs   [ia];
            const Quat4i ngC = ff.neighCell[ia];
            Vec3d dir = ffl.apos[ia] - ffl.apos[ng.x] + ff.shifts[ngC.x]; 
            dir.normalize();
            //printf( "[ia=%i,j=%i,ipbc=%i] dir(%g,%g,%g) \n", ia, ng.x, ngC.x, dir.x,dir.y,dir.z );
            //out->push_back( Vec3i{ ia, ng.x, ngC.x } );

            E+=findHbonds_PBC( ia, Rcut, Hcut, cosMin, dir, *out );
        }
        //printf( "MolWorld_sp3::findHbonds_PBC() DONE\n" );
        return E;
    }
// ========= Manipulation with the molecule

/**
 * Shifts the positions of atoms in the given selection by the specified displacement vector d.
 *
 * @param n The number of atoms in the selection.
 * @param selection An array of indices representing the atoms in the selection.
 * @param d The displacement vector to shift the atoms by.
 */
void shift_atoms ( int n, int* selection, Vec3d d                          ){ move  ( n, selection, ff.apos, d           ); };
/**
 * Rotates the atoms specified by the given selection around the specified axis ax by the specified angle phi.
 *
 * @param n The number of atoms in the selection.
 * @param selection An array of integers representing the indices of the atoms to be rotated.
 * @param p0 The pivot point around which the atoms will be rotated.
 * @param ax The axis of rotation.
 * @param phi The angle of rotation in radians.
 */
void rotate_atoms( int n, int* selection, Vec3d p0, Vec3d ax, double phi   ){ rotate( n, selection, ff.apos, p0, ax, phi ); };
/**
 * Shifts the atoms in the given selection by a specified distance l along the line connecting the atoms with indices ia0 and ia1.
 *
 * @param n         The number of atoms in the selection.
 * @param selection An array of indices representing the atoms in the selection.
 * @param ia0       The index of the first atom
 * @param ia1       The index of the second atom
 * @param l         The distance to shift the atoms by.
 */
void shift_atoms ( int n, int* selection, int ia0, int ia1, double l              ){ Vec3d d=(ff.apos[ia1]-ff.apos[ia0]).normalized()*l; move( n, selection, ff.apos, d); };
/**
 * Rotates the specified atoms in the selection around an axis defined by the atoms with indices iax0 and iax1 by the specified angle phi.
 *
 * @param n The number of atoms in the selection.
 * @param selection An array of indices representing the atoms to be rotated.
 * @param ia0 The index of the reference atom.
 * @param iax0 The index of the starting atom of the rotation axis.
 * @param iax1 The index of the ending atom of the rotation axis.
 * @param phi The angle of rotation in radians.
 */
void rotate_atoms( int n, int* selection, int ia0, int iax0, int iax1, double phi ){ rotate( n, selection, ff.apos, ff.apos[ia0], (ff.apos[iax1]-ff.apos[iax0]).normalized(), phi ); };

/**
 * Splits the selection at a specified bond index.
 * 
 * @param ib The bond index to split at.
 * @param selection The selection array to modify. If set to 0, the global manipulation_sel array will be used.
 * @return The number of atoms in the resulting selection.
 */
int splitAtBond( int ib, int* selection ){
    bool bGlob=(selection==0); 
    if(bGlob){ selection=manipulation_sel; }
    int n = MM::splitByBond( ib, ff.nbonds, ff.bond2atom, ff.apos, selection, manipulation_ax, manipulation_p0 );
    if(bGlob){ manipulation_nsel=n; }
    return n;
}

int selectByType( int itype, bool bByElement=false ){
    selection.clear();
    for(int i=0; i<ffl.natoms; i++){ 
        int it = ffl.atypes[i];
        if( bByElement ){ it = params.atypes[it].element; }
        if( it == itype ){ selection.push_back( i ); }
    }
    return selection.size();
}

/**
 * Selects atoms within a rectangular region defined by two points in 3D space.
 * 
 * @param p0 The first point defining the rectangular region.
 * @param p1 The second point defining the rectangular region.
 * @param rot The rotation matrix to transform the points to the desired coordinate system.
 */
int selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){
    printf( "MolWorld_sp3::selectRect() p0(%g,%g,%g) p1(%g,%g,%g) \n", p0.x,p0.y,p0.z, p1.x,p1.y,p1.z );
    Vec3d Tp0,Tp1,Tp;
    //Mat3d rot = (Mat3d)cam.rot;
    rot.dot_to(p0,Tp0);
    rot.dot_to(p1,Tp1);
    _order(Tp0.x,Tp1.x);
    _order(Tp0.y,Tp1.y);
    Tp0.z=-1e+300;
    Tp1.z=+1e+300;
    selection.clear();
    for(int i=0; i<ff.natoms; i++ ){
        rot.dot_to(ff.apos[i],Tp);
        if( Tp.isBetween(Tp0,Tp1) ){
            selection.push_back( i );
        }
    }
    return selection.size();
}

int selectFragment( int ifrag ){
    //printf( "MolWorld_sp3::selectFragment() ifrag %i \n", ifrag );
    selection.clear();
    if(ifrag<0)selection.size();
    // ToDo: this is not most efficient way to do it, but for know we don't have reliable way to get all atoms which belongs to the fragment when atoms are reordered
    for(int i=0; i<builder.atoms.size(); i++ ){
        const MM::Atom& a = builder.atoms[i];
        if( a.frag == ifrag ){
            selection.push_back( i );
        }
    }
    return selection.size();
}

int selectAllBonded( int ia ){
    selection_set.clear();
    selection.clear();
    if( (ia<0)||(ia>=builder.atoms.size()) ){  printf( "ERROR : MolWorld_sp3::selectAllBonded(ia=%i) but builder.atoms.size(%i) \n" );  exit(0); };
    if( builder.atoms[ia].iconf == -1 ){  // if this atom is capping
        for(int ib=0; ib<builder.bonds.size(); ib++){
            const Vec2i b = builder.bonds[ib].atoms;
            if(b.x==ia){ selection_set.insert( b.y ); }else if(b.y==ia){ selection_set.insert( b.x ); }
        }
    }else{
        selection_set.insert( ia );
    }
    int nfound = 1;
    while( nfound>0){
        int osz = selection_set.size();
        for( int ia : selection_set ){
            int ic = builder.atoms[ia].iconf;
            if( ic < 0 )continue;  
            MM::AtomConf& conf = builder.confs[ic];
            for( int j = 0; j<conf.nbond; j++ ){
                int ib = conf.neighs[j];
                int ja = builder.bonds[ib].getNeighborAtom(ia);
                selection_set.insert( ja );
            }
        }
        nfound = selection_set.size() - osz;
    }
    for( int ia : selection_set ){ selection.push_back( ia ); }
    return selection.size();
}

// ========= Geometry operations on selected atoms

void selectionFromBuilder(){ for(int i: builder.selection){ selection.push_back(i); }; }


bool trySel( int*& sel, int& n ){ 
    if(sel==0){
        sel=selection.data(); 
        n  =selection.size();
        return true; 
    }
    return false;
}

Vec3d center( bool dotIt=false, int* sel=0, int n=-1 ){
    trySel( sel, n );
    Vec3d c = Vec3dZero;
    for(int i=0; i<n; i++){  // ToDo: is it better to do it with min/max ?
        c.add( ffl.apos[selection[i]] ); 
    }
    c.mul( 1./n );
    printf( "MolWorld_sp3::center() cog(%g,%g,%g) selection.size()=%i dotIt=%i \n", c.x,c.y,c.z, selection.size(), dotIt );
    if(dotIt){ for(int i=0; i<n; i++){ ffl.apos[selection[i]].sub(c); } }
    return c;
}

Mat3d getInertiaTensor( int* sel=0, int n=-1, Vec3d* cog=0 ){
    trySel( sel, n );
    Vec3d c;
    if(cog){ c=*cog; }else{ c=center(0, sel, n ); }
    Mat3d I = Mat3dZero;
    for(int i=0; i<n; i++){
        int ia = selection[i];
        Vec3d p = ffl.apos[ia];
        Vec3d d = p - c;
        //printf( "getInertiaTensor()[%i->%i] d(%g,%g,%g) p(%g,%g,%g) \n", i, ia,  d.x,d.y,d.z, p.x,p.y,p.z );
        I.addOuter(d,d, 1.0 );
    }
    //printf( "getInertiaTensor() I(%g,%g,%g) (%g,%g,%g) (%g,%g,%g) \n", I.a.x,I.a.y,I.a.z, I.b.x,I.b.y,I.b.z, I.c.x,I.c.y,I.c.z );
    return I;
}

Mat3d alignToAxis( Vec3i ax={2,1,0}, Mat3d* I_=0, Vec3d* cog=0, bool doIt=true, int* sel=0, int n=-1 ){
    trySel( sel, n );
    Vec3d c; if(cog){ c=*cog;  }else{ c = center(0, sel, n ); }
    Mat3d I; if(I_ ){ I = *I_; }else{ I = getInertiaTensor( sel, n, &c ); }
    Vec3d ls; I.eigenvals(ls);
    printf( "alignToAxis() evals(%g,%g,%g) \n", ls.x,ls.y,ls.z );
    Mat3d rot;
    I.eigenvec( ls.array[ax.x], rot.a );
    I.eigenvec( ls.array[ax.y], rot.b );
    I.eigenvec( ls.array[ax.z], rot.c );
    if(doIt){
        for(int i=0; i<n; i++){
            int ia = selection[i];
            Vec3d p; rot.dot_to( ffl.apos[ia]-c, p );
            ffl.apos[ia] = p + c; 
        }
    }
    return rot;
}

virtual int deleteAtomSelection(){
    return builder.deleteAtoms( selection.size(), selection.data() );   
}

void clearSelections(){
    selection.clear();  
    selection_set.clear();    
    builder.selection.clear();
}

void selectAll    (){ selection.clear(); for(int i=0; i<nbmol.natoms; i++)selection.push_back(i); };
void selectInverse(){ std::unordered_set<int> s(selection.begin(),selection.end());selection.clear(); for(int i=0; i<nbmol.natoms;i++) if( !s.contains(i) )selection.push_back(i);  };

int fragmentsByBonds(){
    builder.frags.clear();
    for( MM::Atom& a : builder.atoms ){ a.frag = -1; }
    //int nfound = 1;
    for( int ia=0; ia<builder.atoms.size(); ia++ ){
        if( builder.atoms[ia].frag == -1 ){
            int nfound = selectAllBonded( ia );
            int ifrag = builder.frags.size();
            //builder.startFragment();
            builder.frags.push_back( MM::Fragment() );
            printf( "fragmentsByBonds() ifrag %i nfound %i \n", ifrag, nfound );
            for( int i : selection ){
                builder.atoms[i].frag = ifrag;
                //builder.frags[ifrag].atoms.push_back( i );
            }
        }    
    }
    selection_set.clear();
    selection.clear();
    builder.randomFragmentCollors();
    return builder.frags.size();
}

/**
 * Performs a translation scan along the specified direction d for a given number of steps.
 *
 * @param n The number of atoms in the selection.
 * @param selection An array of indices representing the atoms in the selection.
 * @param d The direction of translation.
 * @param nstep The number of steps for the translation scan.
 * @param Es An array to store the energy values at each step (optional).
 * @param trjName The name of the trajectory file to write the coordinates (optional).
 * @param bAddjustCaps Flag indicating whether to adjust the caps (default: false).
 */
void scanTranslation_ax( int n, int* selection, Vec3d d, int nstep, double* Es,const char* trjName, bool bAddjustCaps=false ){
    //if(selection==0){ selection=manipulation_sel; n=manipulation_nsel; }
    //Vec3d d=(*(Vec3d*)(vec)); 
    d.mul(1./nstep);
    FILE* file=0;
    if(trjName){ file=fopen( trjName, "w" ); }
    for(int i=0; i<nstep; i++){
        if(file){ toXYZ(tmpstr,false,file); };
        double E = eval();
        if(Es)Es[i]=E;
        move( n, selection, nbmol.apos, d);
    }
    if(file){ fclose(file); }
}

/**
 * Performs a translation scan along the line connecting the atoms with indices ia0 and ia1 for a given number of steps.
 *
 * @param n The number of atoms in the system.
 * @param selection An array of indices representing the atoms to be translated.
 * @param ia0 The index of the first atom involved in the translation.
 * @param ia1 The index of the second atom involved in the translation.
 * @param l The length of the translation vector.
 * @param nstep The number of steps in the translation scan.
 * @param Es An array to store the calculated energies.
 * @param trjName The name of the trajectory file to save the scan results.
 * @param bAddjustCaps Flag indicating whether to adjust the caps of the molecular system.
 */
void scanTranslation( int n, int* selection, int ia0, int ia1, double l, int nstep, double* Es, const char* trjName, bool bAddjustCaps=false ){ Vec3d d=(nbmol.apos[ia1]-nbmol.apos[ia0]).normalized()*l; scanTranslation_ax(n,selection, d, nstep, Es, trjName , bAddjustCaps); };

/**
 * Performs a rotational scan of the selected atoms around a given axis ax for a given number of steps.
 *
 * @param n The number of atoms in the selection.
 * @param selection An array of integers representing the indices of the atoms in the selection.
 * @param p0 The reference point for rotation.
 * @param ax The axis of rotation.
 * @param phi The total rotation angle.
 * @param nstep The number of steps in the rotational scan.
 * @param Es An array to store the energy values at each step (optional).
 * @param trjName The name of the trajectory file to write the rotational scan data (optional).
 */
void scanRotation_ax( int n, int* selection, Vec3d p0, Vec3d ax, double phi, int nstep, double* Es, const char* trjName ){
    //if(p0==0) p0=(double*)&manipulation_p0;
    //if(ax==0) ax=(double*)&manipulation_ax;
    //if(selection==0){selection=manipulation_sel; n=manipulation_nsel; }
    double dphi=phi/nstep;
    FILE* file=0;
    if(trjName){ file=fopen( trjName, "w" ); }
    //printf( "MolWorld_sp3_simple::scanRotation_ax() nstep=%i phi=%g nsel=%i ax(%g,%g,%g) p0(%g,%g,%g) \n", nstep, phi, n, ax.x,ax.y,ax.z,  p0.x,p0.y,p0.z );
    for(int i=0; i<nstep; i++){
        double E = eval();
        Vec3d tq = torq( n, nbmol.apos, nbmol.fapos, p0, selection );
        if(file){  sprintf(tmpstr,"# rotScan[%i] E=%g tq=(%g,%g,%g)", i, E, tq.x,tq.y,tq.z );  toXYZ(tmpstr, false, file, true ); };
        if(Es)Es[i]=E;
        //printf("scanRotation_ax[%i] phi %g E %g \n", i, phi*i, E );
        ffl.rotateNodes(n, selection, p0, ax, dphi );
    }
    if(file){ fclose(file); }
}

/**
 * Scans the rotation of the selected atoms around the axis defined by the atoms with indices iax0 and iax1 for a given number of steps.
 *
 * @param n The number of atoms in the system.
 * @param selection An array of indices representing the atoms to be rotated.
 * @param ia0 The index of the central atom.
 * @param iax0 The index of the starting atom of the rotation axis.
 * @param iax1 The index of the ending atom of the rotation axis.
 * @param phi The rotation angle in radians.
 * @param nstep The number of steps in the rotation.
 * @param Es An array to store the calculated energies.
 * @param trjName The name of the trajectory file to save the rotated configurations.
 */
void scanRotation( int n, int* selection,int ia0, int iax0, int iax1, double phi, int nstep, double* Es, const char* trjName ){ Vec3d ax=(nbmol.apos[iax1]-nbmol.apos[iax0]).normalized(); scanRotation_ax(n,selection, nbmol.apos[ia0], ax, phi, nstep, Es, trjName ); };


void scanAngleToAxis_ax( int n, int* selection, double r, double R, Vec3d p0, Vec3d ax, int nstep, double* angs, double* Es, const char* trjName ){
    //printf( "scanAngleToAxis_ax()\n" );
    FILE* file=0;
    if(trjName){ file=fopen( trjName, "w" ); }
    for(int i=0; i<nstep; i++){
        double ang = angs[i];
        Vec2d cs; cs.fromAngle(ang);
        //printf( "scanAngleToAxis_ax[%i] ang=%g cs(%g,%g)\n", i, ang, cs.x, cs.y );
        for(int j=0; j<n; j++){
            int ia = selection[j];
            Vec3d d = ffl.apos[ia]-p0;
            double l2=d.norm2();
            double c = ax.dot(d);
            //double sign=(c>0)1:-1;
            d.add_mul(ax, -c );                 // remove axial component
            //d.mul( 1.0/sqrt(l2-c*c) );
            d.mul( (r*cs.x + R)/sqrt(l2-c*c) ); // renormalize radial compent 
            d.add_mul(ax, r*cs.y  );            // add back new axial component 
            ffl.apos[ia] = p0 + d;
        }
        double E = eval();
        if(file){ sprintf(tmpstr,"# scanAngleToAxis[%i] E=%g ", i, E);  toXYZ(tmpstr, false, file, true );  };
        if(Es)Es[i]=E;
    }
    if(file){ fclose(file); }
}


/**
 * Calculates automaticaly the charges for a molecular system.
 *
 * @param natoms The number of atoms in the system.
 * @param atypes An array of atom types.
 * @param REQs An array of quaternion values representing the atomic charges.
 * @param neighs An array of quaternion values representing the neighboring atoms.
 * @param nMaxIter The maximum number of iterations for the charge calculation (default: 10).
 * @param K The charge force constant (default: 1.0).
 * @param K0 The charge force constant for the total charge (default: 1.0).
 * @param Q0 The desired total charge (default: 0.0).
 * @param dt The time step for the charge calculation (default: 0.1).
 * @param damping The damping factor for the charge calculation (default: 0.1).
 * @param Fconv The convergence criterion for the charge calculation (default: 1e-6).
 */
__attribute__((hot))  
void autoCharges(int natoms, int* atypes, Quat4d* REQs, Quat4i* neighs, int nMaxIter=10, double K=1.0, double K0=1.0, double Q0=0.0, double dt=0.1, double damping=0.1, double Fconv=1e-6 ){
    std::vector<double> fs(natoms);
    std::vector<double> vs(natoms,0.);
    for(int i=0; i<nMaxIter; i++){
        // --- eval charge forces
        double Qtot  = 0;
        for(int ia=0; ia<natoms; ia++ ){
            int* ng            = neighs[ia].array;
            int  it            = atypes[ia];
            const ElementType* et = params.elementOfAtomType(it);
            double qi   = REQs[ia].z;
            double Qtot = qi;
            double f    = et->Eaff + et->Ehard*qi;
            for(int j=0;j<4;j++){
                int ja             = ng[j];
                if(ja<0) continue;
                f           = REQs[ja].z * K;
            }
            //nsum++;
            fs[ia] =f;
        }
        double fQ = K0*(Q0-Qtot);
        // --- move
        double cdamp = 1-damping;
        double f2=0;
        for(int ia=0; ia<natoms; ia++ ){
            double qi = REQs[ia].z;
            float f   = fs[ia] + fQ*qi;
            float v   = vs[ia];
            f2+=f*f;
            v*=cdamp;
            v+=f;
            qi        += v*dt;
            vs  [ia]   = v;
            REQs[ia].z = qi;

        }
        if(f2<(Fconv*Fconv)){ break; }
    }
}


/**
 * Prints the current state of the switches used in MolWorld_sp3_simple.
 * The switches include bCheckInvariants, bPBC, bNonBonded, bMMFF, ffl.doAngles,
 * ffl.doPiSigma, ffl.doPiPiI, and ffl.bSubtractAngleNonBond.
 */
virtual void printSwitches(){
    printf( "MolWorld_sp3_simple::printSwitches() bCheckInvariants=%i bPBC=%i bNonBonded=%i bMMFF=%i ffl.doAngles=%i ffl.doPiSigma=%i ffl.doPiPiI=%i ffl.bSubtractAngleNonBond=%i \n", bCheckInvariants, bPBC, bNonBonded, bMMFF, ffl.doAngles, ffl.doPiSigma, ffl.doPiPiI, ffl.bSubtractAngleNonBond );
}



bool addSnapshot(bool ifNew = false, char* fname = 0)
{
    if(!gopt.database){
        gopt.database = new MolecularDatabase();
        gopt.database->setDescriptors();
    }
//     loadNBmol("butan");
//     int nMembers = database->getNMembers();

//     database->addMember(&nbmol);

//     std::string comment = "#" + std::to_string(nMembers);
//     const char* cstr = comment.c_str();
//     std::string trjName1 = "trj" + std::to_string(nMembers) + "_1.xyz";        
//     const char* cstr1 = trjName1.c_str();  
//     FILE* file1=0;
//     file1=fopen( cstr1, "w" );
//     params.writeXYZ(file1, &nbmol, cstr);
//     fclose(file1);
// srand(457);


//             double move_x = 0;//randf() * 0.5;
//             double move_y = 0;//randf() * 0.5;
//             double move_z = 0;//randf() * 0.5;
//     double angle = 0;
// if(nMembers==0)angle = 0;//randf() * 360;
// //printf("angle: %g\n", angle);
//     //if(nMembers==1)move_x = 0.2;//;
//         Vec3d axis ={1,1,1};// {randf(), randf(), randf()};
//         Vec3d p0 = Vec3dZero;//{randf() * 5 - 2.5, randf() * 5 - 2.5, randf() * 5 - 2.5};
//         for (int i = 0; i < nbmol.natoms; i++)
//         {            
            

//             //if(i==0 && nMembers == 0) move_x = 5;//randf() * 0.5;
//             nbmol.apos[i].add({move_x, move_y, move_z});
//             nbmol.apos[i].rotate(2 * 3.14159 / 360 * angle, axis, p0);
//          }
// //if(nMembers>1)    printf("database->compareAtoms(&nbmol, 0) %g\n", database->compareAtoms(&nbmol, 0));
// //printf("ff.nbonds %i\n", ff.nbonds);
// loadNBmol("changed_butan");
// //buildMolecule_xyz("changed_butan");

// int nbFixed = 2;
// int* fixed = new int[nbFixed];
//     while(nMembers < 100){
//     nMembers = database->getNMembers();

//     database->addMember(&nbmol);

//     std::string comment = "#" + std::to_string(nMembers);
//     const char* cstr = comment.c_str();
//     std::string trjName1 = "trj" + std::to_string(nMembers) + "_1.xyz";        
//     const char* cstr1 = trjName1.c_str();  
//     FILE* file1=0;
//     file1=fopen( cstr1, "w" );
//     params.writeXYZ(file1, &nbmol, cstr);
//     fclose(file1);


        
//         fixed[0] = -1;
//         fixed[1] = -1;
//         database->as_rigid_as_possible(&nbmol, 0, ff.nbonds, ff.bond2atom, nbFixed, fixed);
        
//         //printf("nbmol.neighs[0]: %i %i %i %i\n", nbmol.neighs[0].x, nbmol.neighs[0].y, nbmol.neighs[0].z, nbmol.neighs[0].w);
//     }
// delete[] fixed;



//    return true;

    if (fname)
    {
        loadNBmol(fname);
        //buildMolecule_xyz(fname);
    }
    if (ifNew)
    {
        int ID = gopt.database->addIfNewDescriptor(&nbmol);
        if (ID != -1)
        {
            printf("Same as %d\n", ID);
            return false;
        }
    }
    else
    {
        gopt.database->addMember(&nbmol);
    }
    return true;
}

void printDatabase(){
    //if(database) database->print();

    int Findex = 0;
    std::vector<double> a;
    std::vector<double> b;
    std::vector<int> boundaryRules;
    for (int i = 0; i < nbmol.natoms * 3; i++)
    {
        a.push_back(-20);
        b.push_back(20);
        boundaryRules.push_back(1);
    }
    double Fstar = 0;
    int maxeval = 15000;
    int nRelax = 500;
    int nExplore = __INT_MAX__;
    int bShow = 1;
    int bSave = 0;
    int bDatabase = 0;
    runGlobalOptimization(Findex, &a, &b, &boundaryRules, Fstar, maxeval, nRelax, nExplore, bShow, bSave, bDatabase);

}

double computeDistance(int i, int j){
    if(gopt.database) return gopt.database->computeDistance(i, j);
    printf("Error, database does not exist"); return -1;
}

void runGlobalOptimization(int Findex, std::vector<double>* a, std::vector<double>* b, std::vector<int>* boundaryRules, 
    double Fstar, int maxeval, int nRelax, int nExploring, int bShow, int bSave, int bDatabase){
    //gopt.init_heur(nbmol, params, Findex, a, b, boundaryRules, Fstar, maxeval, nRelax, nExploring, bShow, bSave, bDatabase)


    std::vector<double> par_mut = {0.01};
    std::vector<double> par_alg = {1e-5, 100, maxeval, 100};    
    
    constrs.loadBonds("hexan-dicarboxylic.cons");
    printf("\nbConstrains: %i\n", bConstrains);
    bConstrains = true;

    
    gopt.init_heur(&nbmol, &params, Findex, a, b, boundaryRules, 0, maxeval, nRelax, nExploring, bShow, 2, 0);
    gopt.SPSA(0, &par_mut, &par_alg);
    //gopt.randomBrutal(&nbmol, &params, 0, 0, 0, &boundaryRules, 0, maxeval, bShow, 3, &par_mut, &par_alg);
    //nbmol.print();
printf("after optimization\n");
}








};

#endif