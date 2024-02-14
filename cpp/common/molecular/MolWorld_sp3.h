﻿
#ifndef MolWorld_sp3_h
#define MolWorld_sp3_h

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <omp.h>
#include <memory>
#include "IO_utils.h"

//#include "testUtils.h"
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"

#include "MMFFparams.h"
static MMFFparams* params_glob;

//#include "raytrace.h"
#include "Forces.h"
#include "MMFFsp3.h"
#include "MMFFsp3_loc.h"
#include "MMFFf4.h"

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

#include "datatypes_utils.h"

#include "arrayAlgs.h"
#include "SVG_render.h"



#include "MolecularDatabase.h"



/**
 * @class MolWorld_sp3
 * @brief A class representing a molecular world in a 3D space.
 *
 * This class inherits from SolverInterface and provides functionality for manipulating and simulating molecular systems.
 * It contains various member variables and methods for handling molecular data, force fields, dynamics, and more.
 */
class MolWorld_sp3 : public SolverInterface { public:
    const char* data_dir     = "common_resources";
    const char* xyz_name     = "input";
    //const char* surf_name    = "surf";
    const char* surf_name       = 0;
    const char* substitute_name = 0;    int isubs;
    //const char* lvs_name     ="input.lvs";
    //const char* surflvs_name ="surf.lvs";
    const char* smile_name   = 0;
    const char* constr_name  = 0;
    Vec3i nMulPBC  = Vec3iZero;                    // Vector for multiple Periodic Boundary Conditions

    //const char* trj_fname    = "trj.xyz";
    const char* trj_fname    = 0;
    int savePerNsteps = 1;
    OptLog opt_log;
    Vec3i nPBC_save{1,1,1};                         // Vector for saving Periodic Boundary Conditions


    double fAutoCharges=-1;
    bool bEpairs = false;

    bool bCellBySurf=false;
    int   bySurf_ia0 =0;
    Vec2d bySurf_c0=Vec2dZero;                      // Vector for surface coordinates
    Vec2d bySurf_lat[2];                            // Array for surface lattice

    Mat3d new_lvec=Mat3dIdentity;                   // Matrix for new lattice vectors

    Mat3d debug_rot; // This was used for debuging molecular orientation

	// Building
	MMFFparams   params;
	MM::Builder  builder;
	SMILESparser smiles;

	// Force-Fields & Dynamics
	MMFFsp3      ff;
    MMFFsp3_loc  ffl;
    MMFFf4       ff4;
	Constrains   constrs;
	//NBFF_old   nff;
    NBFF surf, nbmol; // Non-bonded force fields for surface and molecule
    GridFF gridFF;    // Grid force field
    RigidBodyFF rbff; // Rigid body force field
    QEq qeq;          // Charge equilibration
    DynamicOpt opt;
    DynamicOpt optRB; // rigid body optimizer

    LimitedGraph<N_NEIGH_MAX> graph; // used to find bridges in the molecule, and other topology-related algorithms

    GlobalOptimizer gopt;

    GridShape MOgrid;

    SVG_render svg;

    int  iterPerFrame=50;
    int  iParalel=0; 
    int  iParalelMax=1;
    int  iParalelMin=0;
    bool bOcl=false; // used only in Ocl version

    double gridStep = 0.1; 
    //double gridStep = 0.2; 
    //Vec3i nPBC{0,0,0};   // just debug
    //Vec3i nPBC{1,1,0};
    Vec3i nPBC{1,3,0};
    int    npbc       = 0;
    Vec3d* pbc_shifts = 0;

    bool bConverged = false; // Boolean flag indicating if the simulation has converged
    double Etot = 0;         // Total energy of the system
    double maxVcog = 1e-9;   // Maximum velocity of the center of gravity
    double maxFcog = 1e-9;   // Maximum force on the center of gravity
    double maxTg = 1e-1;     // Maximum torque on the center of gravity
    double Kmorse = -1.0;    // Morse potential constant

    // force-eval-switchefs
    int  imethod=0;

    bool doBonded          = false; // If true, the program considers bonded interactions.
    bool bNonBonded        = true;  // If true, the program considers non-bonded interactions.
    bool bConstrains       = false; // If true, the program applies certain constraints (like bond lengths or angles).
    bool bSurfAtoms        = false; // If true, the program considers surface atoms for some calculations.
    bool bGridFF           = false; // If true, the program uses a grid-based force field.
    bool bPlaneSurfForce   = false; // If true, the program applies a force related to a plane surface.
    bool bMMFF             = true;  // If true, the program uses the Merck Molecular Force Field for calculations.
    bool bRigid            = false; // If true, the program treats the molecule as rigid (no internal degrees of freedom).
    bool bOptimizer        = true;  // If true, the program uses an optimizer to find the lowest energy state.
    bool bPBC              = false; // If true, the program uses Periodic Boundary Conditions.
    bool bCheckInvariants  = true;  // If true, the program checks for invariants (properties that don't change during the simulation).
    bool bRelaxPi          = false; // If true, the program allows relaxation of pi systems (double and triple bonds).
    bool bChargeUpdated    = false; // If true, the program updates the charges on atoms during the simulation.
    bool bAnimManipulation = false; // If true, the program allows for manipulation of the animation of the molecule.

    Vec3d cog,vcog,fcog,tqcog;
    int nloop=0;

	// Selecteion & Manipulasion
	std::vector<int> selection;
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

    Mat3d* dlvec = 0;

    // lattice scan
    bool   bLatScan      = false;
    Mat3d* latscan_dlvec = 0;
    Vec2i  latscan_n{0,0};

	// IO temp & aux
	FILE* xyz_file=0;
	char* tmpstr;

    double Kpick = -2.0;

    int itest = 0;

    //int nSystems    = 1;
    int iSystemCur  = 0;    // currently selected system replica


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
    
    }}}
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
        printf( "setConstrains %i \n", i );
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
    //printf( "MolWorld::solve(nmax=%i,tol=%g)\n", nmax, tol );
    //ffl.print_apos();
    //printf("ffl.lvec\n"    ); printMat( ffl.lvec    );
    //printf("ffl.invLvec\n" ); printMat( ffl.invLvec );
    //printf("npbc %i nPBC(%i,%i,%i) \n", npbc, nPBC.x,nPBC.y,nPBC.z );
    //nmax = 10;
    long t0=getCPUticks();
    int nitr = run_omp( nmax, opt.dt_max, tol, 1000.0, -1. );
    long t=(getCPUticks()-t0); printf( "time run_omp[%i] %g[Mtick] %g[ktick/iter]  %g[s] %g[ms/iter]\n", nitr, t*1e-6, t*1.e-3/nitr, t*tick2second, t*tick2second*1000/nitr  );
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
    ffl.initPi( pbc_shifts );
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
    if(npbc!=npbc_eval){ printf( "ERORRO in MolWorld_sp3::makePBCshifts() final ipbc(%i)!=nbpc(%i) => Exit()\n", npbc_eval,npbc ); exit(0); }
    return npbc;
}

void printPBCshifts(){
    printf("printPBCshifts():\n");
    for(int i=0; i<npbc; i++){ printf("pbc_shift[%i] (%6.3f,%6.3f,%6.3f)\n", i, pbc_shifts[i].x,pbc_shifts[i].y,pbc_shifts[i].z ); }
}

/**
 * Calculates the auto charges for the molecules in the molecular world.
 * This function assigns charges to the atoms based on the given force field parameters.
 * It also applies constraints to certain atom types and performs charge relaxation using molecular dynamics.
 * After the charges are calculated, they are copied to the corresponding molecules in the molecular world.
 */
void autoCharges(){
    if(verbosity>0)printf("MolWorld_sp3::autoCharges() \n");
    qeq.realloc( ff.natoms );
    params.assignQEq ( ff.natoms, ff.atype, qeq.affins, qeq.hards );
    int iconstr = params.getAtomType("E");    //printf("constrain type %i \n", iconstr );
    qeq.constrainTypes( ff.atype, iconstr );
    qeq.relaxChargeMD ( ff.apos, 1000, 1e-2, 0.1, 0.1 );
    copy( qeq.n, 1, 0, (double*)qeq.qs, 3, 2, (double*)nbmol.REQs );
    bChargeUpdated=true;
}

/**
 * Automatically calculates the number of periodic boundary conditions (nPBC) based on minimum length.
 * 
 * @param cell The cell dimensions represented by a 3x3 matrix (Mat3d).
 * @param nPBC The number of periodic boundary conditions in each direction (x, y, z). This parameter will be modified by the function.
 * @param Lmin The minimum length for calculating the number of periodic boundary conditions. Default value is 30.0.
 */
static void autoNPBC( const Mat3d& cell, Vec3i& nPBC, double Lmin=30.0 ){
    if(nPBC.x!=0){ nPBC.x=(int)Lmin/cell.a.norm(); }
    if(nPBC.y!=0){ nPBC.y=(int)Lmin/cell.b.norm(); }
    if(nPBC.z!=0){ nPBC.z=(int)Lmin/cell.c.norm(); }
    printf("autoNPBC(): (%i,%i,%i) \n", nPBC.x, nPBC.y, nPBC.z );
}

/**
 * Saves the grid data in XSF format for debugging purposes.
 * 
 * @param bE Whether to save the energy-related grids (default: true)
 * @param bFz Whether to save the force-related grids (default: true)
 * @param bComb Whether to save the combined forcefield grid (default: true)
 * @param testREQ The test Quat4d value for evaluating the combined forcefield (default: Quat4d{ 1.487, 0.02609214441, 0., 0.})
 */
void saveGridXsfDebug( bool bE=true, bool bFz=true, bool bComb=true, Quat4d testREQ=Quat4d{ 1.487, 0.02609214441, 0., 0.} ){
    // not testREQ.y [eV^0.5] = sqrt(Eii), 
    // e.g. for Hydrogen 0.02609214441 ev^0.5 = sqrt( 0.0006808 eV )
    // e.g. for Carbon   0.06106717612 ev^0.5 = sqrt( 0.0037292 eV )
    if(bE){
        if(gridFF.FFPaul) gridFF.grid.saveXSF( "FFLond_E.xsf", (float*)gridFF.FFLond, 4,3  );
        if(gridFF.FFLond) gridFF.grid.saveXSF( "FFelec_E.xsf", (float*)gridFF.FFelec, 4,3  );
        if(gridFF.FFelec) gridFF.grid.saveXSF( "FFPaul_E.xsf", (float*)gridFF.FFPaul, 4,3  );
    }
    if(bFz){
        if(gridFF.FFPaul) gridFF.grid.saveXSF( "FFLond_z.xsf", (float*)gridFF.FFLond, 4,2  );
        if(gridFF.FFLond) gridFF.grid.saveXSF( "FFelec_z.xsf", (float*)gridFF.FFelec, 4,2  );
        if(gridFF.FFelec) gridFF.grid.saveXSF( "FFPaul_z.xsf", (float*)gridFF.FFPaul, 4,2  );
    }
    // ---- Save combined forcefield
    if(bComb){
        Quat4f * FFtot = new Quat4f[gridFF.grid.getNtot()];
        gridFF.evalCombindGridFF ( testREQ, FFtot );
        gridFF.grid.saveXSF( "E_PLQ.xsf",  (float*)FFtot, 4, 3, gridFF.natoms, gridFF.atypes, gridFF.apos );
        delete [] FFtot;
    }
}

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
virtual void initGridFF( const char * name, bool bGrid=true, bool bSaveDebugXSFs=false, double z0=NAN, Vec3d cel0={-0.5,-0.5,0.0}, bool bAutoNPBC=true ){
    if(verbosity>0)printf("MolWorld_sp3::initGridFF(%s,bGrid=%i,z0=%g,cel0={%g,%g,%g})\n",  name, bGrid, z0, cel0.x,cel0.y,cel0.z  );
    sprintf(tmpstr, "%s.lvs", name );
    if( file_exist(tmpstr) ){  gridFF.grid.loadCell( tmpstr, gridStep );  gridFF.bCellSet=true; }
    if( !gridFF.bCellSet ){
        bGridFF=false; 
        printf( "WARRNING!!! GridFF not initialized because %s not found\n", tmpstr );
        return;
    }
    if(bGrid){
        gridFF.grid.center_cell( cel0 );
        bGridFF=true;
        gridFF.bindSystem(surf.natoms, surf.atypes, surf.apos, surf.REQs );
        if( isnan(z0) ){  z0=gridFF.findTop();   if(verbosity>0) printf("GridFF::findTop() %g \n", z0);  };
        gridFF.grid.pos0.z=z0;
        //gridFF.grid.pos0.z=-5;
        if(verbosity>1)gridFF.grid.printCell();
        gridFF.allocateFFs();
        //gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false, {1,1,0}, bSaveDebugXSFs );
        gridFF.nPBC=Vec3i{1,1,0};
        if(bAutoNPBC){ autoNPBC( gridFF.grid.cell, gridFF.nPBC, 20.0 ); }
        //gridFF.nPBC = (Vec3i){0,0,0};
        //gridFF.nPBC = (Vec3i){1,1,0};
        //gridFF.nPBC = (Vec3i){10,10,0};
        gridFF.lvec = gridFF.grid.cell;     // ToDo: We should unify this
        gridFF.makePBCshifts     ( gridFF.nPBC, gridFF.lvec );
        gridFF.setAtomsSymetrized( gridFF.natoms, gridFF.atypes, gridFF.apos, gridFF.REQs, 0.1 );
        //bSaveDebugXSFs=true;
        gridFF.tryLoad( "FFelec.bin", "FFPaul.bin", "FFLond.bin", false );
        gridFF.log_z( "initGridFF_iz_ix0_iy0.log" ,0,0);
        if(bSaveDebugXSFs)saveGridXsfDebug();
        bGridFF   =true; 
        //bSurfAtoms=false;
    }
    gridFF.shift0 = Vec3d{0.,0.,-2.0};
    //gridFF.shift0 = Vec3d{0.,0.,0.0};
    gridFF.evalCheck();
}

/**
 * Initializes the non-bonded molecules
 * 
 * @param na The number of atoms.
 * @param apos An array of Vec3d representing the positions of the atoms.
 * @param fapos An array of Vec3d representing the final positions of the atoms.
 * @param atypes An array of integers representing the atom types.
 * @param bCleanCharge A boolean indicating whether to clean the charges of atoms not present in the builder.
 */
void initNBmol( int na, Vec3d* apos, Vec3d* fapos, int* atypes, bool bCleanCharge=true ){
    if(verbosity>0)printf( "MolWorld_sp3::initNBmol() na %i \n", na  );
	nbmol.bindOrRealloc( na, apos, fapos, 0, atypes );    
    //nbmol.bindOrRealloc( na, apos, fapos, 0, 0 );   
    //builder.export_atypes( nbmol.atypes );     
	builder.export_REQs( nbmol.REQs   );    
    nbmol  .makePLQs   ( gridFF.alphaMorse );  
    ffl.PLQs=nbmol.PLQs; 
    if(bCleanCharge)for(int i=builder.atoms.size(); i<na; i++){ nbmol.REQs[i].z=0; }  // Make sure that atoms not present in Builder has well-defined chanrge                       
    params.assignREs( na, nbmol.atypes, nbmol.REQs, true, false  );
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
    sprintf(tmpstr, "%s.xyz", name );
	int ret = params.loadXYZ( tmpstr, surf.natoms, &surf.apos, &surf.REQs, &surf.atypes, 0, &gridFF.grid.cell );
	if     ( ret<0 ){ printf("ERROR in MolWorld_sp3::loadSurf() file(%s) not found => Exit() \n",         tmpstr ); exit(0); }
    if     ( ret==0){ printf("ERROR in MolWorld_sp3::loadSurf() no lattice vectors in (%s) => Exit() \n", tmpstr ); exit(0); }
    else if( ret>0 ){ gridFF.grid.updateCell(gridStep); gridFF.bCellSet=true;  }
    //gridFF.grid.printCell(); 
    if(verbosity>0)printf("MolWorld_sp3::loadSurf(%s) 1 natoms %i apos %li atyps %li \n", name, surf.natoms, (long)surf.apos, (long)surf.atypes  );
    //surf.print();
	bSurfAtoms=true;
    initGridFF( name,bGrid,bSaveDebugXSFs,z0,cel0 );
	return true;
}

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
    builder.printAtomConfs(false);
    builder.printBonds();
    printf( " ===================== Substitute molecule START !!! \n");
    Molecule* mol = new Molecule(); mol->init_xyz( fname, &params, true );
    //Vec3i axSwap={2,1,0};
    //Vec3i axSwap={2,0,1}; // THIS ONE
    //Vec3i axSwap={0,1,2};
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false );
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &axSwap, &debug_rot );
    int ja = builder.substituteMolecule( mol, Vec3dZ, ib, ipivot, false, 0, &debug_rot );
    //builder.substituteMolecule( mol, Vec3dZ, 4, 0, false, &(Vec3i{2,1,0}), &debug_rot );
    builder.addCappingTypesByIz(1);
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
    graph.bPrint = true;
    graph.bridge();
    printf( "graph.found.size() %i \n", graph.found.size() );
    return graph.found.size();
}



/**
 * Loads the geometry from a file and inserts it into the molecular world.
 * 
 * @param name The name of the file containing the geometry.
 * @return The index of the inserted fragment in the builder.
 */
int loadGeom( const char* name ){ // TODO : overlaps with buildFF()
    if(verbosity>0)printf("MolWorld_sp3::loadGeom(%s)\n",  name );
    // ------ Load geometry
    sprintf(tmpstr, "%s.xyz", name ); //printf("tmpstr=`%s`\n", tmpstr);
    int imol  = builder.loadMolType( tmpstr, name );
    int iret  = builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
    int ifrag = builder.frags.size()-1;
    if(iret<0){ printf("!!! exit(0) in MolWorld_sp3::loadGeom(%s)\n", name); exit(0); }
    builder.addCappingTypesByIz(1);
    //for( int it : builder.capping_types ){ printf( "capping_type[%i] iZ=%i name=`%s`  \n", it, builder.params->atypes[it].iZ, builder.params->atypes[it].name ); };
    builder.tryAddConfsToAtoms( 0, -1 );
    builder.cleanPis();
    if(verbosity>2)builder.printAtomConfs(false);
    //builder.export_atypes(atypes);
    // ------- Load lattice vectros
    sprintf(tmpstr, "%s.lvs", name );
    //builder.printAtomConfs(true);
    if( file_exist(tmpstr) ){
        builder.bPBC=true;
        readMatrix( tmpstr, 3, 3, (double*)&builder.lvec );
    }
    bPBC=builder.bPBC;  //printf( "builder.bPBC %i \n", builder.bPBC );
    if( bPBC ){ builder.autoBondsPBC(); }
    else      { builder.autoBonds();    }
    builder.checkNumberOfBonds( true, true );
    if(verbosity>2)builder.printBonds ();
    return ifrag;
}

// int loadmol(const char* fname_mol ){
//     int imol = builder.loadMolType( fname_mol, "molecule" );
//     builder.insertFlexibleMolecule( imol, {0,0,0}, Mat3dIdentity, -1 );
//     return imol;
// }

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
void setOptimizer( int n, double* ps, double* fs ){
    //opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
    opt.bindOrAlloc( n, ps, 0, fs, 0 );
    double dtopt=ff.optimalTimeStep(); 
    if(verbosity>0)printf("MolWorld_sp3::setOptimizer(): optimnal time step = %g \n", dtopt);
    opt.initOpt( dtopt );
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
    params.printAngleTypes();
    params.printDihedralTypes();
}

/**
 * Builds a molecule from an XYZ file.
 * 
 * @param xyz_name The name of the XYZ file.
 * @return The index of the built molecule.
 */
int buildMolecule_xyz( const char* xyz_name ){
    int ifrag = loadGeom( xyz_name );
    int ia0=builder.frags[ifrag].atomRange.a;
    int ic0=builder.frags[ifrag].confRange.a;
    //builder.printBonds();
    //builder.printAtomConfs(true, false );
    if( fAutoCharges>0 )builder.chargeByNeighbors( true, fAutoCharges, 10, 0.5 );
    if(substitute_name) substituteMolecule( substitute_name, isubs, Vec3dZ );
    if( builder.checkNeighsRepeat( true ) ){ printf( "ERROR: some atoms has repating neighbors => exit() \n"); exit(0); };
    builder.autoAllConfEPi  ( ia0 ); 
    builder.setPiLoop       ( ic0, -1, 10 );
    if(bEpairs)builder.addAllEpairsByPi( ia0=0 );    
    //builder.printAtomConfs(false, false );
    //builder.printAtomConfs(false, true );
    builder.assignAllBondParams();    //if(verbosity>1)
    builder.finishFragment(ifrag);    
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
    //builder.printAtoms();          
    //if( builder.checkBondsOrdered( false, true ) ) { printf("ERROR Bonds are not ordered => exit"); exit(0); };
    if( builder.checkBondsInNeighs(true) ) { 
        printf("ERROR some bonds are not in atom neighbors => exit"); 
        exit(0); 
    };
    builder.numberAtoms();
    builder.sortConfAtomsFirst();
    //builder.printAtomConfs(false,true);
    builder.checkBondsOrdered( true, false );
    builder.assignTypes();
    //builder.printAtomTypes();
    if( ffl.bTorsion ){ builder.assignTorsions( true, true ); }  //exit(0);

    // --- find bridge bonds
    findBridgeBonds();
    
    builder.toMMFFsp3_loc( ffl, true, bEpairs );   if(ffl.bTorsion){  ffl.printTorsions(); } // without electron pairs
    if(ffl.bEachAngle){ builder.assignAnglesMMFFsp3  ( ffl, false      ); ffl.printAngles();   }  //exit(0);
    builder.toMMFFf4     ( ff4, true, bEpairs );  //ff4.printAtomParams(); ff4.printBKneighs(); 
    builder.toMMFFsp3    ( ff , true, bEpairs );
    
    ffl.flipPis( Vec3dOne );
    ff4.flipPis( Vec3fOne );
    if(bPBC){  
        //ff.printAtomParams();
        ff.bPBCbyLvec = true;
        ff .setLvec( builder.lvec);
        ffl.setLvec( builder.lvec);   
        ff4.setLvec((Mat3f)builder.lvec);
        npbc = makePBCshifts( nPBC, builder.lvec );
        ffl.bindShifts(npbc,pbc_shifts);
        ff4.makeNeighCells  ( nPBC );       
        //ffl.makeNeighCells( nPBC );      
        ffl.makeNeighCells( npbc, pbc_shifts ); 
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
    makeMMFFs();
    initNBmol( ffl.natoms, ffl.apos, ffl.fapos, ffl.atypes ); 
    setNonBond( bNonBonded );
    bool bChargeToEpair=true;
    //bool bChargeToEpair=false;
    if(bChargeToEpair){
        int etyp=-1; etyp=params.atomTypeDict["E"];
        ff.chargeToEpairs( nbmol.REQs, -0.2, etyp );  
    }
    nbmol.evalPLQs(gridFF.alphaMorse);
    { // check FFS
        //ffl.printAtomParams();
        //printf("npbc %i\n", npbc ); ffl.printNeighs();
        //builder.printBonds();
        //printf("!!!!! builder.toMMFFsp3() DONE \n");
        idebug=1;
        //printf( "!!!!!!!!!!!!!!!! ffl.checkREQlimits(); \n" );
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

/**
 * Initializes the MolWorld_sp3 object.
 * 
 * @param bGrid A boolean indicating whether to initialize the grid.
 */
virtual void init( bool bGrid ){
    // params.init("common_resources/ElementTypes.dat", "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat" );
	// builder.bindParams(&params);
    // params.printAtomTypes(true);
    // //params.printAtomTypeDict();
    // //params.printBond();
    if( params.atypes.size() == 0 ){
        initParams( "common_resources/ElementTypes.dat", "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat", "common_resources/AngleTypes.dat", "common_resources/DihedralTypes.dat" );
    }
    gopt.solver = this;
    params_glob = &params;
    builder.verbosity=verbosity;
    if(verbosity>0){
        printf("\n#### MolWorld_sp3::init()\n");
        if(smile_name   )printf("smile_name  (%s)\n", smile_name );
        if(data_dir     )printf("data_dir    (%s)\n", data_dir );
        if(xyz_name     )printf("xyz_name    (%s)\n", xyz_name );
        if(surf_name    )printf("surf_name   (%s)\n", surf_name );
        if(substitute_name)printf("substitute_name  (%s)\n", substitute_name );
        //if(lvs_name     )printf("lvs_name    (%s)\n", lvs_name );
        //if(surflvs_name )printf("surflvs_name(%s)\n", surflvs_name );
        printf( "MolWorld_sp3::init() bMMFF %i bRigid %i \n", bMMFF, bRigid );
        //for(int i=0; i<10; i++){ float x = -1.0+i*0.2; printf( "x %g ix %i wx %g \n", x, (int)x, x+1-(int)(x+1.5) ); }; exit(0);
    }
    if(surf_name )loadSurf( surf_name, bGrid, idebug>0 );
    if ( smile_name ){               
        insertSMILES( smile_name );    
        builder.addAllCapTopo();       
        builder.randomizeAtomPos(1.0); 
        bMMFF=true;
    }else if ( xyz_name ){
        if( bMMFF ){ 
            buildMolecule_xyz( xyz_name );
        }else{
            loadNBmol( xyz_name ); 
            if(bRigid)initRigid();
        }
    }
    if(bMMFF){      
        makeFFs();
    }
    builder.setup_atom_permut( true );
    if(constr_name ){ constrs.loadBonds( constr_name, &builder.atom_permut[0], 0 );  }
    if(dlvec       ){ add_to_lvec(*dlvec);    }  // modify lattice after initialization - it helps to build constrained systems 

    database = new MolecularDatabase;

    //builder.printAtoms();
    //printf( "MolWorld_sp3::init() ffl.neighs=%li ffl.neighCell-%li \n", ffl.neighs, ffl.neighCell );
    //ffl.printNeighs();
    if(verbosity>0) printf( "... MolWorld_sp3::init() DONE \n");
}

/**
 * @brief Clears the MolecularWorld object.
 * 
 * This function clears the MolecularWorld object by deallocating memory and resetting variables.
 * 
 * @param bParams Flag indicating whether to clear the parameters as well. Default is true.
 */
virtual void clear( bool bParams=true ){
    //printf("MolWorld_sp3.clear() \n");
    builder.clear();
    ffl.dealloc();
    ff.dealloc();
    ff4.dealloc();
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

    constrs.clear();

    if(bParams){
        params.clear();
    }
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

/**
 * Evaluates the energy using the ff4 potential and returns the result.
 * 
 * @return The energy evaluated using the ff4 potential.
 */
double eval_f4(){
    pack( ff4.natoms, ffl.apos , ff4.apos  );
    pack( ff4.nnode,  ffl.pipos, ff4.pipos );
    double E = ff4.eval();
    //ff4.move_GD( 0.01);
    unpack( ff4.natoms, ffl. apos, ff4. apos  );
    unpack( ff4.natoms, ffl.fapos, ff4.fapos  );
    unpack( ff4.nnode,  ffl. pipos,ff4. pipos );
    unpack( ff4.nnode,  ffl.fpipos,ff4.fpipos );
    //for(int i=0; i<ff4.nnode; i++) printf("pi[%i] <fpi,pi> %g |pi| %g \n", i, ffl.fpipos[i].dot( ffl.pipos[i] ), ffl.pipos[i].norm() );
    return E;   
};

/**
 * Sets the non-bonded flag for the molecular world.
 * 
 * @param bNonBonded A boolean value indicating whether non-bonded interactions should be considered.
 */
void setNonBond( bool bNonBonded ){
    ffl.bSubtractAngleNonBond = bNonBonded;
    ff4.bSubtractAngleNonBond = bNonBonded;
    if(bNonBonded){
        ffl.REQs = nbmol.REQs;
        ff .REQs = nbmol.REQs;
        if(ff4.REQs==0){
            ff4.REQs = new Quat4f[nbmol.natoms];
            for(int i=0; i<nbmol.natoms; i++ ){ ff4.REQs[i] = (Quat4f)nbmol.REQs[i]; };
        }
    }
}

/**
 * Calculates the energy of the molecular system.
 *
 * @return The calculated energy of the molecular system.
 */
double eval( ){

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
        //E += ff .eval();
        E += ffl.eval(); 
        //ffl.printDEBUG(  false, false );
        //for(int i=0; i<nbmol.natoms; i++){ printf("atom[%i] f(%g,%g,%g)\n", i, nbmol.fapos[i].x,nbmol.fapos[i].y,nbmol.fapos[i].z ); }   
        
        //printf( "ffl.lvec\n" );    printMat( ffl.lvec    );
        //printf( "ffl.invLvec\n" ); printMat( ffl.invLvec );
        //exit(0);
        //E += eval_f4();
        //printf( "atom[0] nbmol(%g,%g,%g) ff(%g,%g,%g) ffl(%g,%g,%g) \n", nbmol.apos[0].x,nbmol.apos[0].y,nbmol.apos[0].z,  ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,  ffl.apos[0].x,ffl.apos[0].y,ffl.apos[0].z );
        
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
    //ffl.printDEBUG(  false, false );
    //exit(0);


    return E;
}


/**
 * Performs relaxation of the molecular system with FIRE algorithm.
 * 
 * @param niter The number of relaxation iterations to perform.
 * @param Ftol The force tolerance for convergence. Defaults to 1e-6.
 * @param bWriteTrj Flag indicating whether to write trajectory information. Defaults to false.
 * @return True if relaxation converged, false otherwise.
 */
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
virtual int run( int nstepMax, double dt=-1, double Fconv=1e-6, int ialg=2, double* outE=0, double* outF=0 ){ 
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
        if( (trj_fname) && (itr%savePerNsteps==0) ){
            sprintf(tmpstr,"# %i E %g |F| %g", itr, Etot, sqrt(F2) );
            saveXYZ( trj_fname, tmpstr, false, "a", nPBC_save );
        }
        if(verbosity>1){ printf("[%i] Etot %g[eV] |F| %g [eV/A] \n", itr, Etot, sqrt(F2) ); };
        if(F2<F2conv){
            bConverged=true;
            if(verbosity>0){ printf("Converged in %i iteration Etot %g[eV] |F| %g[eV/A] <(Fconv=%g) \n", itr, Etot, sqrt(F2), Fconv ); };
            if( trj_fname ){
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

virtual void MDloop( int nIter, double Ftol = 1e-6 ){
    //printf( "MolWorld_sp3::MDloop() \n" );
    //ff.doPiPiI  =false;
    //ff.doPiPiT  =false;
    //ff.doPiSigma=false;
    //ff.doAngles =false;

    /*
    //ff.cleanAll();
    for(int itr=0; itr<nIter; itr++){
        double E = eval();
        //ckeckNaN_d( nbmol.natoms, 3, (double*)nbmol.fapos, "nbmol.fapos" );
		//if( bPlaneSurfForce )for(int i=0; i<ff.natoms; i++){ ff.fapos[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) ); }
        //printf( "apos(%g,%g,%g) f(%g,%g,%g)\n", ff.apos[0].x,ff.apos[0].y,ff.apos[0].z,   ff.fapos[0].x,ff.fapos[0].y,ff.fapos[0].z );
        //if(bCheckInvariants){ checkInvariants(maxVcog,maxFcog,maxTg); }
        if(ipicked>=0){ pullAtom( ipicked );  }; // printf( "pullAtom(%i) E=%g\n", ipicked, E ); };
        //ff.fapos[  10 ].set(0.0); // This is Hack to stop molecule from moving
        //opt.move_GD(0.001);
        //opt.move_LeapFrog(0.01);
        //double f2=opt.move_MDquench();

        double f2=0;
        //opt.damp_max = 0.0;
        //opt.cv_kill  = 0.5;
        //double f2=opt.move_FIRE();
        for(int i=0; i<ffl.nvecs; i++){ f2+=ffl.move_atom_MD( i, 0.05, 1000.0, 0.9 ); } 
        //double f2=opt.move_VSpread( 0.1,  0.01, 0.05 ); 
        //printf( "[%i] E= %g [eV] |F|= %g [eV/A]\n", nloop, E, sqrt(f2) );
        //double f2=1;
        if(f2<sq(Ftol)){
            bConverged=true;
        }
        nloop++;
    }
    */
    
    verbosity = 1;
    

    


    //run_omp( 500, opt.dt, 1e-6, 1000.0 );
    //ffl.run_omp( 10, 0.05, 1e-6, 1000.0 );
    run_omp( nIter, 0.05, 1e-6, 1000.0 );
    //run_omp( 100, 0.05, 1e-6, 1000.0 );
    //run_omp( 1, opt.dt, 1e-6, 1000.0 );
    //run_omp( 2, opt.dt, 1e-6, 1000.0 );
    //run_omp( 100, opt.dt, 1e-6, 1000.0 );
    //run_omp( 500, 0.05, 1e-6, 1000.0 );
    //run_omp( 500, 0.05, 1e-6, 1000.0 );
    
    //run( nIter );
    
    bChargeUpdated=false;
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
int run_omp( int niter_max, double dt, double Fconv=1e-6, double Flim=1000, double timeLimit=0.02, double* outE=0, double* outF=0 ){
    //printf( "run_omp() niter_max %i dt %g Fconv %g Flim %g timeLimit %g outE %li outF %li \n", niter_max, dt, Fconv, Flim, timeLimit, (long)outE, (long)outF );
    if(dt>0){ opt.setTimeSteps(dt); }
    long T0 = getCPUticks();
    double E=0,F2=0,F2conv=Fconv*Fconv;
    double ff=0,vv=0,vf=0;
    int itr=0,niter=niter_max;

    //#pragma omp parallel shared(E,F2,ff,vv,vf,ffl) private(itr)
    #pragma omp parallel shared(niter,itr,E,F2,ff,vv,vf,ffl,T0)
    while(itr<niter){
        if(itr<niter){
        //#pragma omp barrier
        #pragma omp single
        {E=0;F2=0;ff=0;vv=0;vf=0;}
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
            if(bPBC){ E+=ffl.evalLJQs_ng4_PBC_atom_omp( ia ); }
            else    { E+=ffl.evalLJQs_ng4_atom_omp    ( ia ); } 
            if   (bGridFF){ E+= gridFF.addForce          ( ffl.apos[ia], ffl.PLQs[ia], ffl.fapos[ia], true ); }        // GridFF
            //if     (bGridFF){ E+= gridFF.addMorseQH_PBC_omp( ffl.apos[ia], ffl.REQs[ia], ffl.fapos[ia]       ); }    // NBFF
            
            if(ipicked==ia){ 
                const Vec3d f = getForceSpringRay( ffl.apos[ia], pick_hray, pick_ray0,  Kpick ); 
                ffl.fapos[ia].add( f );
            }

        }
        if(ffl.bTorsion){
            #pragma omp for reduction(+:E)
            for(int it=0; it<ffl.ntors; it++){ 
                E+=ffl.eval_torsion(it); 
            }
        }
        if(bConstrains){
            #pragma omp single
            {
                //printf( "run_omp() constrs[%i].apply()\n", constrs.bonds.size() );
                E+= constrs.apply( ffl.apos, ffl.fapos, &ffl.lvec );
            }
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
                ffl.move_atom_FIRE( i, opt.dt, 10000.0, opt.cv, opt.renorm_vf*opt.cf );
                //ffl.move_atom_FIRE( i, dt, 10000.0, 0.9, 0 ); // Equivalent to MDdamp
            }
        }
        
        }
        //#pragma omp barrier
        #pragma omp single
        { 
            Etot=E;
            itr++; 
            if(timeLimit>0){
                double t = (getCPUticks() - T0)*tick2second;
                if(t>0.02){ 
                    niter=0; 
                    if(verbosity>0)printf( "run_omp() ended due to time limit after %i nsteps ( %6.3f [s]) \n", itr, t ); 
                }
            }
            if(F2<F2conv){ 
                niter=0; 
                double t = (getCPUticks() - T0)*tick2second;
                if(verbosity>0)printf( "run_omp() CONVERGED in %i/%i nsteps E=%g |F|=%g time= %g [ms]( %g [us/%i iter])\n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr );
            }

            if(outE){ outE[itr]=Etot; }
            if(outF){ outF[itr]=F2;   }

            if( (trj_fname) && ( (itr%savePerNsteps==0) ||(niter==0) ) ){
                sprintf(tmpstr,"# %i E %g |F| %g", itr, Etot, sqrt(F2) );
                //printf( "run_omp::save() %s \n", tmpstr );
                saveXYZ( trj_fname, tmpstr, false, "a", nPBC_save );
            }

            //printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(F2), omp_get_num_threads() ); 
            //{printf( "step[%i] dt %g(%g) cv %g cf %g cos_vf %g \n", itr, opt.dt, opt.dt_min, opt.cv, opt.cf, opt.cos_vf );}
            //if(verbosity>2){printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, sqrt(F2), omp_get_num_threads() );}
        }
    }{
    double t = (getCPUticks() - T0)*tick2second;
    if(itr>=niter_max)if(verbosity>0)printf( "run_omp() NOT CONVERGED in %i/%i E=%g |F|=%g time= %g [ms]( %g [us/%i iter]) \n", itr,niter_max, E, sqrt(F2), t*1e+3, t*1e+6/itr, itr );
    }
    return itr;
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
    sprintf( str_tmp, "lvs %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f    %8.3f %8.3f %8.3f %s", ffl.lvec.a.x, ffl.lvec.a.y, ffl.lvec.a.z, ffl.lvec.b.x, ffl.lvec.b.y, ffl.lvec.b.z, ffl.lvec.c.x, ffl.lvec.c.y, ffl.lvec.c.z, comment );
    return params.saveXYZ( fname, (bNodeOnly ? ffl.nnode : ffl.natoms) , nbmol.atypes, nbmol.apos, str_tmp, nbmol.REQs, mode, true, nPBC, ffl.lvec ); 
}
//int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false){ return params.saveXYZ( fname, (bNodeOnly ? ff.nnode : ff.natoms) , ff.atype, ff.apos, comment, nbmol.REQs ); }
//int saveXYZ(const char* fname, const char* comment="#comment", bool bNodeOnly=false){ return params.saveXYZ( fname, (bNodeOnly ? ff.nnode : ff.natoms) , nbmol.atypes, nbmol.apos, comment, nbmol.REQs ); }

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

/**
 * Selects atoms within a rectangular region defined by two points in 3D space.
 * 
 * @param p0 The first point defining the rectangular region.
 * @param p1 The second point defining the rectangular region.
 * @param rot The rotation matrix to transform the points to the desired coordinate system.
 */
void selectRect( const Vec3d& p0, const Vec3d& p1, const Mat3d& rot ){
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


/**
 * Scans angles to axis.
 *
 * This function scans a set of angles and modifies the positions of atoms based on the given parameters.
 * The positions of atoms are adjusted by removing the axial component, renormalizing the radial component,
 * and adding back a new axial component. The energy of the system is evaluated at each step.
 *
 * @param n         The number of atoms in the selection.
 * @param selection An array of atom indices representing the selection.
 * @param r         The radial component parameter.
 * @param R         The renormalization parameter.
 * @param p0        The reference point for the axial component.
 * @param ax        The axial vector.
 * @param nstep     The number of steps to perform.
 * @param angs      An array of angles to scan.
 * @param Es        An array to store the evaluated energies.
 * @param trjName   The name of the trajectory file to write the positions to (optional).
 */
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


MolecularDatabase* database;


void addSnapshot(){

    int nMembers = database->getNMembers();
    database->setDescriptors(&params, &nbmol);



srand(time(0));

    double angle = randf() * 360;
    double move_x = randf() * 0.5;
    double move_y = randf() * 0.5;
    double move_z = randf() * 0.5;
    Vec3d axis = {randf(), randf(), randf()};
    Vec3d p0 = {randf()*5-2.5, randf()*5-2.5, randf()*5-2.5};
    for (int i = 0; i < nbmol.natoms; i++)
    {
        nbmol.apos[i].rotate(2 * 3.14159 / 360 * angle, axis, p0);
        nbmol.apos[i].add({move_x, move_y, move_z});    
    }



    

    for (int i = 0; i<nMembers; i++)
    {
        printf("%d <-> %d: %lf\n", nMembers, i, database->compareAtoms(&nbmol, i));
    }
    
    database->testHash(&nbmol);    
    
    std::string trjName1 = "trj" + std::to_string(nMembers) + "_0.xyz";        
    const char* cstr1 = trjName1.c_str();  
    FILE* file1=0;
    file1=fopen( cstr1, "w" );
    params.writeXYZ(file1, &nbmol, "#comment");
    fclose(file1);
    // for (int i = 0; i < nbmol.natoms; i++)
    // {
    //     nbmol.apos[i].rotate(-2 * 3.14159 / 360 * 20, Vec3d{0, 0, 1}, Vec3d{0, 0, 0});
    // }
    //params.writeXYZ(file2, &nbmol, "#comment");
    
    //fclose(file2);
    //if(nMembers>0) printf("database->compareDescriptors(nMembers-1, 0): %f\n", database->compareDescriptors(nMembers, 0));


    
    // if(nMembers >=1){
    //     if(database->compareDescriptors(nMembers, 0)){
    //         printf( "MolWorld_sp3_simple::addSnapshot() WARNING: duplicate descriptor \n" );
    //     }
    //     else{
    //         printf( "MolWorld_sp3_simple::addSnapshot() new descriptor \n" );
    //     }
    // }
    //nbmol.natoms = 6;
}

void printDatabase(){
    database->print();
}








};

#endif
