#ifndef MolGUI_h
#define MolGUI_h

#include "Camera.h"
#include "Draw.h"
#include "GLMesh.h"
#include "Renderer.h"
#include "Vec3.h"
#include "globals.h"

#include "macroUtils.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "quaternion.h"
#include "testUtils.h"
#include "IO_utils.h"

#include <SDL2/SDL.h>

#include "Draw3D.h"
#include "SDL_utils.h"
//#include "Solids.h"

#include "MolWorld_sp3.h"

#include "raytrace.h"
#include "Draw3D_Molecular.h"  // it needs to know MMFFparams
#include "MolecularDraw.h"
#include "MarchingCubes.h"
#include "GUI.h"
#include "Console.h"
#include "EditorGizmo.h"
#include "SimplexRuler.h"
#include "AppSDL2OGL_3D.h"

#include "AtomsInGrid.h"

#include <chrono>

#include "MolGUI_tests.h"

#include "FullscreenShader.h"

#include "MethodDict.h"


#include "DipoleMap.h"

//using Action  = std::function<void(double val)>; 
//using CommandDict = std::unordered_map<std::string,>;

void plotNonBondLine( const NBFF& ff, Quat4d REQi, double Rdamp, Vec3d p1, Vec3d p2, int n, Vec3d up=Vec3dZ, bool bForce=false ){
    Vec3d d = (p2-p1)*(1.0/n);
    Vec3d p = p1;
    double fsc = up.norm();
    //opengl1renderer.begin(GL_LINE_STRIP);
    opengl1renderer.begin(GL_TRIANGLE_STRIP);
    for(int i=0; i<=n; i++){
        Quat4d fe = ff.evalLJQs( p, REQi, Rdamp );
        //Draw3D::drawVecInPos( p, f*0.1, 0.1 );
        Vec3d pf;
        if(bForce){ pf = p + fe.f*fsc; } // Force 
        else      { pf = p + up*fe.e;  } // Energy
        //printf( "plotNonBondLine[%i] p(%g,%g,%g) E=%g pf(%g,%g,%g) \n", i, p.x,p.y,p.z, fe.w, pf.x,pf.y,pf.z );
        opengl1renderer.vertex3f( pf.x, pf.y, pf.z );
        opengl1renderer.vertex3f( p.x, p.y, p.z );
        p.add(d);
    }
    opengl1renderer.end();
}

Vec2d evalNonBondGrid2D( const NBFF& ff, Quat4d REQi, double Rdamp, Vec2i ns, double*& Egrid, Vec3d p0, Vec3d a, Vec3d b, Vec3d up=Vec3dZ, bool bForce=false ){
    //printf( "evalNonBondGrid2D().1 ns(%i,%i) @Egrid=%li \n", ns.x, ns.y, (long)Egrid );
    if(Egrid==0)Egrid = new double[(ns.a)*(ns.b)];
    //printf( "evalNonBondGrid2D().2 ns(%i,%i) @Egrid=%li \n", ns.x, ns.y, (long)Egrid );
    a.mul( 1.0/ns.a );
    b.mul( 1.0/ns.b );
    up.normalize();
    Vec2d Erange = Vec2d{1e+300,-1e+300};
    for(int iy=0; iy<ns.b; iy++){
        for(int ix=0; ix<ns.a; ix++){
            Vec3d p = p0 + a*ix + b*iy;
            Quat4d fe = ff.evalLJQs( p, REQi, Rdamp );
            double val;
            if(bForce){ val = fe.f.dot(up); } // Force 
            else      { val = fe.e;         } // Energy
            Erange.enclose( val );
            Egrid[ix*ns.a+iy] = val;
        }
    }
    return Erange;
}

// namespace Draw3D{
// } // namespace Draw3D


void drawDipoleMapGrid( DipoleMap& dipoleMap, Vec2d sc=Vec2d{1.,1.}, bool radial=true, bool azimuthal=true ){
    int nphi = dipoleMap.nphi;
    int nr   = dipoleMap.particles.size()/nphi;
    if(radial)
    for(int ip=0; ip<nphi; ip++ ){
        opengl1renderer.begin(GL_LINE_STRIP);
        for(int ir=0; ir<nr; ir++ ){
            int i = ip + ir*nphi;
            //printf( "drawDipoleMap()[%i|ip=%i,ir=%i]\n", i, ip, ir );
            Vec3d p = dipoleMap.particles[i];
            p.z    +=  dipoleMap.FE[i].e*sc.x + dipoleMap.FE2[i].e*sc.y;
            opengl1renderer.vertex3f( p.x, p.y, p.z );
        }
        opengl1renderer.end();
    }
    if(azimuthal)
    for(int ir=0; ir<nr; ir++ ){
       opengl1renderer.begin(GL_LINE_STRIP);
        for(int ip=0; ip<nphi; ip++ ){
            int i = ip + ir*nphi;
            Vec3d p = dipoleMap.particles[i];
            p.z    +=  dipoleMap.FE[i].e*sc.x + dipoleMap.FE2[i].e*sc.y;
            opengl1renderer.vertex3f( p.x, p.y, p.z );
        }
        opengl1renderer.end();
    }
}

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

//void command_example(double x, void* caller);

static const char* testFragmentShader = R"(
#version 100

uniform sampler2D uTexture;
varying mediump vec2 fUV;

void main(){
    gl_FragColor = texture2D(uTexture, fUV);
}

)";

class MolGUI : public AppSDL2OGL_3D { public:

    //using Action  = std::function<void(double val)>;
    using Action = std::function<void(GUIAbstractPanel* panel)>; 

    Vec3d  rotation_center = Vec3dZero;
    Vec3d  rotation_axis   = Vec3dZ;
    double rotation_step   = 15.0 * (M_PI/180.0);

    //int ipicked    = -1; // picket atom 
    //int ibpicked   = -1; // picket bond
    //int iangPicked = -1; // picket angle
    //Vec3d* picked_lvec = 0;
    //int perFrame =  1;
    int perFrame =  100;

    bool bDoMM=true;//,bDoQM=true;
    bool bConverged = false;
    bool bRunRelax  = false;
    double cameraMoveSpeed = 1.0;
    //bool useGizmo=true;
    bool useGizmo=false;
    bool bDrawHexGrid=true;
    bool bHexDrawing=false; 

    bool bConsole=false;


    bool bWriteOptimizerState = true;
    //bool bPrepared_mm = false;
    //bool bPrepared_qm = false;

    double z0_scan = 0.0;

    double subs_iso = 0.05;

    MolWorld_sp3* W=0;

    bool bDipoleMap = false;
    DipoleMap dipoleMap;

    // ---- GUI

    Console console;
    GUI gui;
    GUIPanel*    Qpanel=0;
    EditorGizmo  gizmo;
    SimplexRuler ruler; // Helps paiting organic molecules
    enum class Gui_Mode { base, edit, scan };
    Gui_Mode  gui_mode = Gui_Mode::base;
    //int gui_mode = Gui_Mode::edit;
    DropDownList* panel_Frags=0;
    GUIPanel*     panel_iMO  =0;
    GUIPanel*     panel_AFM  =0;

    GUIPanel*  BondLengh_types=0;
    GUIPanel*  BondLengh_min=0;
    GUIPanel*  BondLengh_max=0;

    // --- NonBond plot
    //bool bDrawNonBond = false;
    bool bDrawNonBondLines=false;
    bool bDrawNonBondGrid =false;
    bool bDrawParticles = false;
    bool hideEp           =false;
    CheckBoxList* panel_NBPlot=0;
    MultiPanel* panel_Edit = 0;
    MultiPanel* panel_NonBondPlot=0;
    MultiPanel* panel_PickedType=0;
    MultiPanel* panel_TestType=0;
    MultiPanel* panel_GridXY=0;

    // https://docs.lammps.org/Howto_tip5p.html
    // TIP5P - Transferable Intermolecular Potential 5 Point   - 5 particles   (  O(0.0e), H1(+0.241e), H1(+0.241), E1(-0.241e), E2(-0.241e) ), E ~ 0.7A from Oxygen
    Quat4d particle_REQ { 1.187, sqrt(0.0006808    ), -0.24, 0.0 };
    Quat4d particle_REQ2{ 1.750, sqrt(0.00260184625), +0.24, 0.0 };
    double particle_Lbond = 1.0;
    std::vector<Quat4d> particles;
    std::vector<Quat4d> particles2;
    std::vector<int>    particlePivots; // index of atom to which the particle is attached


    Dict<Action> panelActions;   // used for binding GUI actions using Lua scripts 

    // this is used for binding to Lua and to GUI
    std::unordered_map<std::string,int> member_offsets_double;
    std::unordered_map<std::string,int> member_offsets_int;

    MethodDict<MolGUI> actions;

    // ---- Visualization params

    int iSystemCur = 0;
    int which_MO  = 0; 


    //float textSize        = 0.025;
    float textSize        = 0.020;
    //float textSize        = 0.015;
    double ForceViewScale = 1.0;

    // ---- nice balls and sticks for school
    // double mm_Rsc         = 0.20;
    // double mm_Rsub        = 0.0;

    // double ForceViewScale = 1.0;
    // double mm_Rsc         = 0.25;
    // double mm_Rsub        = 0.5;

    // ---- small balls and sticks for debugging
    // double ForceViewScale = 100.0;
    double mm_Rsc            = 0.1;
    double mm_Rsub           = 0.0;

    bool bBuilderChanged     = false;

    bool   bViewBuilder      = false;
    //bool   bViewBuilder      = true;
    bool   bViewAxis         = false;
    bool   bViewCell         = true;

    bool   mm_bAtoms         = true;
    bool   bViewMolCharges   = false;
    bool   bViewHBondCharges = false;
    bool   bViewAtomLabels   = true;
    bool   bViewAtomTypes    = false;
    bool   bViewColorFrag    = false;
    bool   bViewBondLabels   = false;
    bool   bViewAtomSpheres  = true;
    bool   bViewAtomForces   = true;
    bool   bViewBondLenghts  = false;
    bool   bViewBonds        = true;
    bool   bViewPis          = false;
    bool   bViewSubstrate    = true;
    bool   isoSurfRenderType = 1;
    bool   bDebug_scanSurfFF = false;
    Quat4d testREQ;
    Quat4f testPLQ;

    double myAngle=0.0;

    Mat3d dlvec { 0.1,0.0,0.0,   0.0,0.0,0.0,  0.0,0.0,0.0 };
    Mat3d dlvec2{ 0.0,0.1,0.0,   0.0,0.0,0.0,  0.0,0.0,0.0 };
    Vec2f mouse_pix;

    // ----- Visualization Arrays - allows to switch between forcefields, and make it forcefield independnet
    int    natoms=0,nnode=0,nbonds=0;
    int*   atypes;
    Vec2i* bond2atom = 0; 
    Vec3d* pbcShifts = 0; 
    Vec3d* apos      = 0;
    Vec3d* fapos     = 0;
    Vec3d* pipos     = 0;
    Vec3d* fpipos    = 0;
    Quat4d* REQs     = 0;

    Quat4i* neighs    = 0;
    Quat4i* neighCell = 0;

    int nsys=0,nvec=0;
    int    * M_neighs    =0;
    int    * M_neighCell =0;
    Quat4f * M_apos      =0;

    Constrains* constrs = 0;

    double* bL0s = 0;

    // ---- Graphics objects    // ToDO: maybe move to globals.h
    int  fontTex=-1,fontTex3D=-1;

    int  ogl_afm       = 0;
    int  ogl_afm_trj   = 0;
    int  ogl_esp       = 0;
    int  ogl_mol       = 0;
    GLMesh_NC ogl_isosurf = GLMesh_NC(GL_TRIANGLES);
    int  ogl_MO        = 0;
    int  ogl_nonBond   = 0;
    int  ogl_Hbonds    = 0;
    int  ogl_trj       = 0;
    int  ogl_surf_scan = 0;

    FullscreenShader TESTfullscreenShader = FullscreenShader(testFragmentShader);

    std::vector<Quat4f> debug_ps;
    std::vector<Quat4f> debug_fs;
    std::vector<Quat4f> debug_REQs;

    std::vector<Vec2i> bondsToShow;
    Vec3d * bondsToShow_shifts = 0; 

    // ----- AFM scan
    GridShape afm_scan_grid{ Vec3d{-10.,-10.,0.0}, Vec3d{10.,10.,8.0}, 0.1 };
    GridShape afm_ff_grid;  //  { Mat3d{}, 0.1 };
    Quat4f    *afm_ff=0,*afm_Fout=0,*afm_PPpos=0; 
    Quat4f    *afm_ps0=0;
    //Vec3d     *afm_ps=0;
    int  afm_iz    = 25;
    int  afm_nconv = 10;
    bool afm_bDf   = true;

    // ----- Atoms in grid
    AtomsInGrid atomsInGrid;

    // ======================= Functions 

    void setPanelAction( int ipanel, const char* name ){ 
        //gui.panels[ipanel]->setCommand( panelActions.get(name) );
        int id = panelActions.getId(name);
        if(id<0){
            printf( "ERROR: setPanelAction(%i,%s) failed\n", ipanel, name );
            exit(0);
        }else{ 
            gui.panels[ipanel]->setCommand( panelActions.vec[id] );
        }
        // ToDo:       cannot convert ‘void (MolGUI::*)()’ to ‘const std::function<void(GUIAbstractPanel*)>&’
    };

	virtual void draw() override;
	virtual void drawHUD() override;
	//virtual void mouseHandling( ) override;
	virtual void eventHandling   ( const SDL_Event& event  ) override;
	virtual void keyStateHandling( const Uint8 *keys )       override;

    void eventMode_default( const SDL_Event& event );
    void eventMode_scan   ( const SDL_Event& event );
    void eventMode_edit   ( const SDL_Event& event );
    void mouse_default( const SDL_Event& event );

	MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_=0 );
    
    //void initMol(MolWorld_sp3* W_);
    //void InitQMMM();
    bool visual_FF_test();
	//int  loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int  loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );
    void tryLoadGridFF();
    //void makeGridFF   (bool recalcFF=false, bool bRenderGridFF=true);
    //void renderGridFF( double isoVal=0.001, int isoSurfRenderType=0, double colorScale = 50. );
    void renderSurfAtoms( Vec3i nPBC, bool bPointCross=false, float qsc=0.5, float Rsc=1, float Rsub=0 );
    //void renderGridFF    ( double isoVal=0.1, int isoSurfRenderType=0, double colorScale = 50. );
    void renderGridFF_new( double isoVal=0.1, int isoSurfRenderType=0, double colorScale = 0.25, Quat4d REQ=Quat4d{ 1.487, sqrt(0.0006808), 0., 0.} );
    void renderESP( Quat4d REQ=Quat4d{ 1.487, 0.02609214441, 1., 0.} );
    void renderAFM( int iz, int offset );
    void renderAFM_trjs( int di );
    void Fz2df( int nxy, int izmin, int izmax, const Quat4f* afm_Fout, float* dfout );
    void makeAFM( int iz=-1 );
    void lattice_scan( int n1, int n2, const Mat3d& dlvec );
    void makeBondLengths0();
    void makeBondColoring( Vec2i typs, Vec2d lrange, double*& clr,bool bNew=true);

    void bindMolecule(int natoms_, int nnode_, int nbonds_, int* atypes_,Vec3d* apos_,Vec3d* fapos_,Quat4d* REQs_, Vec3d* pipos_, Vec3d* fpipos_, Vec2i* bond2atom_, Vec3d* pbcShifts_);
    void bindMolecule(const MolWorld_sp3* W );
    void bindMolWorld( MolWorld_sp3* W );
    void unBindMolecule();
    void drawSystemShifts( int n, const Vec3d* shifts, int i0 );
    inline void drawSystemSingle() { drawSystemShifts(1, &Vec3dZero, 0); }
    void drawBuilder( Vec3i ixyz=Vec3iZero );
    
    void drawPi0s( float sc );
    void  showAtomGrid( char* s, int ia, bool bDraw=true );
    Vec3d showNonBond ( char* s, Vec2i b, bool bDraw=true );
    void tryPlotNonBond();
    void plotNonBondLines();
    void plotNonBondGrid();
    void plotNonBondGridAxis();
    void relaxNonBondParticles( double dt = 0.2, double Fconv = 1e-6, int niter = 1000);
    void drawParticles();
    void drawDipoleMap();
    //void drawParticleNonBonds();
    void showBonds();
    void printMSystem( int isys, int perAtom, int na, int nvec, bool bNg=true, bool bNgC=true, bool bPos=true );
    //void flipPis( Vec3d ax );
    //void drawSystemQMMM();
    void renderOrbital(int i, double iso=0.1);
    void renderDensity(       double iso=0.1);
	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );
    void scanSurfFF      ( int n, Vec3d p0, Vec3d p1, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.2, 0.}, int evalMode=0, int viewMode=0,  double sc=1.0, Vec3d hat=Vec3dX );
    void debug_scanSurfFF( int n, Vec3d p0, Vec3d p1, Quat4d REQ=Quat4d{ 1.487, 0.02609214441, +0.2, 0.}, double sc=NAN );

    void initGUI();
    void updateGUI();
    int  clearGUI(int n);
    void initMemberOffsets();
    void initWiggets();
    //void initWiggets_Run();
    void initCommands();
    void nonBondGUI();
    void drawingHex(double z0);

    // ======================= Functions for AtomsInGrid
    // What is this Non-uniform grid? I don't remember when and why I implemented this.
    // Class representing quadrature mesh created by embeding spherical atoms in 3D rectangular grid. 

    void makeNonuniformGrid(int niter=5, double dt=0.1 ){
        printf( "makeNonuniformGrid() \n" );
        //atomsInGrid.setup( Vec3d pmin, Vec3d pmax_, double step );
        atomsInGrid.snap_corners( W->nbmol );
        //for(int iter=0; iter<niter; iter++){
        //    atomsInGrid.eval_forces();
        //    atomsInGrid.move_points( dt );
        //}
    };

    void relaxNonuniformGrid(int niter=10, double dt=0.1 ){
        printf( "relaxNonuniformGrid() niter=%i double dt=%g \n", niter, dt );
        for(int iter=0; iter<niter; iter++){
            atomsInGrid.eval_forces();
            double f2 = atomsInGrid.move_points( dt );
            printf( "relaxNonuniformGrid[%i] |F|=%g \n", iter, sqrt(f2) );
        }
    };

    void plotNonuniformGrid( bool bFixedOnly=false ){
        for (int i=0; i<atomsInGrid.gpoints.size(); i++){
            GridPointDynamics& gp = atomsInGrid.gpoints[i];

            if( bFixedOnly && !gp.fixed ) continue;
            Vec3i ipos; atomsInGrid.i2ixyz( gp.ic, ipos );
            Vec3d pj; 

            if( gp.fixed ){ 
                Vec3f col = COLOR_RED;
                pj = atomsInGrid.get_gpos({ipos.x,   ipos.y,   ipos.z-1}); Draw3D::drawLine( gp.pos, pj, col );
                pj = atomsInGrid.get_gpos({ipos.x,   ipos.y,   ipos.z+1}); Draw3D::drawLine( gp.pos, pj, col );
                pj = atomsInGrid.get_gpos({ipos.x,   ipos.y-1, ipos.z  }); Draw3D::drawLine( gp.pos, pj, col );
                pj = atomsInGrid.get_gpos({ipos.x,   ipos.y+1, ipos.z  }); Draw3D::drawLine( gp.pos, pj, col );
                pj = atomsInGrid.get_gpos({ipos.x-1, ipos.y,   ipos.z  }); Draw3D::drawLine( gp.pos, pj, col );
                pj = atomsInGrid.get_gpos({ipos.x+1, ipos.y,   ipos.z  }); Draw3D::drawLine( gp.pos, pj, col );
            } else { 
                Vec3f col = COLOR_BLUE;
                Vec3d pj;
                bool bfix; 
                bfix = atomsInGrid.get_gpos2( {ipos.x,   ipos.y,   ipos.z-1}, pj ); if(!bfix) Draw3D::drawLine( gp.pos, pj, col );
                bfix = atomsInGrid.get_gpos2( {ipos.x,   ipos.y,   ipos.z+1}, pj ); if(!bfix) Draw3D::drawLine( gp.pos, pj, col );
                bfix = atomsInGrid.get_gpos2( {ipos.x,   ipos.y-1, ipos.z  }, pj ); if(!bfix) Draw3D::drawLine( gp.pos, pj, col );
                bfix = atomsInGrid.get_gpos2( {ipos.x,   ipos.y+1, ipos.z  }, pj ); if(!bfix) Draw3D::drawLine( gp.pos, pj, col );
                bfix = atomsInGrid.get_gpos2( {ipos.x-1, ipos.y,   ipos.z  }, pj ); if(!bfix) Draw3D::drawLine( gp.pos, pj, col );
                bfix = atomsInGrid.get_gpos2( {ipos.x+1, ipos.y,   ipos.z  }, pj ); if(!bfix) Draw3D::drawLine( gp.pos, pj, col );
            }
            
            /*
            
            Vec3d pj; 
            Vec3i ip;
            ip={ipos.x,   ipos.y,   ipos.z-1};  if( gp.fixed || !atomsInGrid.gpoints[ atomsInGrid.ixyz2i(ip) ].fixed  ) Draw3D::drawLine( gp.pos, atomsInGrid.box2pos(ip) );
            ip={ipos.x,   ipos.y,   ipos.z+1};  if( gp.fixed || !atomsInGrid.gpoints[ atomsInGrid.ixyz2i(ip) ].fixed  ) Draw3D::drawLine( gp.pos, atomsInGrid.box2pos(ip) );
            ip={ipos.x,   ipos.y-1, ipos.z  };  if( gp.fixed || !atomsInGrid.gpoints[ atomsInGrid.ixyz2i(ip) ].fixed  ) Draw3D::drawLine( gp.pos, atomsInGrid.box2pos(ip) );
            ip={ipos.x,   ipos.y+1, ipos.z  };  if( gp.fixed || !atomsInGrid.gpoints[ atomsInGrid.ixyz2i(ip) ].fixed  ) Draw3D::drawLine( gp.pos, atomsInGrid.box2pos(ip) );
            ip={ipos.x-1, ipos.y,   ipos.z  };  if( gp.fixed || !atomsInGrid.gpoints[ atomsInGrid.ixyz2i(ip) ].fixed  ) Draw3D::drawLine( gp.pos, atomsInGrid.box2pos(ip) );
            ip={ipos.x+1, ipos.y,   ipos.z  };  if( gp.fixed || !atomsInGrid.gpoints[ atomsInGrid.ixyz2i(ip) ].fixed  ) Draw3D::drawLine( gp.pos, atomsInGrid.box2pos(ip) );
            */
            // pj = atomsInGrid.get_gpos({ipos.x,   ipos.y,   ipos.z+1}); Draw3D::drawLine( gp.pos, pj );
            // pj = atomsInGrid.get_gpos({ipos.x,   ipos.y-1, ipos.z  }); Draw3D::drawLine( gp.pos, pj );
            // pj = atomsInGrid.get_gpos({ipos.x,   ipos.y+1, ipos.z  }); Draw3D::drawLine( gp.pos, pj );
            // pj = atomsInGrid.get_gpos({ipos.x-1, ipos.y,   ipos.z  }); Draw3D::drawLine( gp.pos, pj );
            // pj = atomsInGrid.get_gpos({ipos.x+1, ipos.y,   ipos.z  }); Draw3D::drawLine( gp.pos, pj );
            
            //Draw3D::drawPointCross( gp.pos, 0.1 );
        }
    }

};

//=================================================
//                   INIT()
//=================================================


void MolGUI::initMemberOffsets(){
    {
    auto& m = member_offsets_double;
    mapByteOffIn(m, z0_scan );
    }
    {
    auto& m = member_offsets_int;
    mapByteOffIn(m, perFrame );
    }
}

void MolGUI::initCommands(){
    panelActions.add( "perFrame", [&](GUIAbstractPanel* p){ perFrame = ((GUIPanel*)p)->getIntVal(); printf( "GUIPanel(%s) set perFrame=%i \n", ((GUIPanel*)p)->caption.c_str(), perFrame ); return 0; } );
}

void MolGUI::initWiggets(){
    printf( "MolGUI::initWiggets() \n" );

    // TODO: adding GUI widgets would be better witth LUA for fast experimentation
    GUI_stepper ylay(1,2 );
    GUI_stepper gx  (1,16);

    // ylay.step(2);
    // Qpanel = new GUIPanel( "Q_pick: ", 5,ylay.x0,5+100,ylay.x1, true, true ); 
    // Qpanel->setRange(-1.0,1.0)
    //       ->setValue(0.0)
    //     //->command = [&](GUIAbstractPanel* p){ zoom = ((GUIPanel *)p)->value; return 0; };
    //       ->setCommand( [&](GUIAbstractPanel* p){ W->nbmol.REQs[W->ipicked].z = ((GUIPanel *)p)->value; return 0; } );
    // (GUIPanel*)gui.addPanel( Qpanel );

    // ------ Table(   Lattice Vectros )

    ylay.step(6);
    Table* tab1 = new Table( 9, sizeof(W->builder.lvec.a), (char*)&W->builder.lvec );
    tab1->addColum( &(W->builder.lvec.a.x), 1, DataType::Double    );
    tab1->addColum( &(W->builder.lvec.a.y), 1, DataType::Double    );
    tab1->addColum( &(W->builder.lvec.a.z), 1, DataType::Double    );
    ((TableView*)gui.addPanel( new TableView( tab1, "lattice", 5, ylay.x0,  0, 0, 3, 3 ) ))->input = new GUITextInput();

    // ------ GUIPanel(   "Zoom: " )

    ylay.step(3); 
    ((GUIPanel*)gui.addPanel( new GUIPanel( "Zoom: ", 5,ylay.x0,5+100,ylay.x1, true, true ) ) )
        ->setRange(5.0,50.0)
        ->setValue(zoom)
        //->command = [&](GUIAbstractPanel* p){ zoom = ((GUIPanel *)p)->value; return 0; };
        ->setCommand( [&](GUIAbstractPanel* p){ zoom = ((GUIPanel *)p)->value; return 0; } );

    // ------ DropDownList(   "Pick Mode:"  )

    ylay.step(3); 
    ((DropDownList*)gui.addPanel( new DropDownList("Pick Mode:",5,ylay.x0,5+100, 3 ) ) )
        ->addItem("pick_atoms")
        ->addItem("pick_bonds")
        ->addItem("pick_angles");


    // ------ DropDownList(   "Fragments:" )

    ylay.step(3); 
    panel_Frags = ((DropDownList*)gui.addPanel( new DropDownList("Fragments:",5,ylay.x0,5+100, 3 ) ) );
    panel_Frags->setCommand( [&](GUIAbstractPanel* me_){ int i=((DropDownList*)me_)->iSelected; printf( "panel_Frags %02i \n", i );  W->selectFragment(i); return 0; } );   

    // ------ DropDownList(   "View Side"  )

    ylay.step(6); 
    ((DropDownList*)gui.addPanel( new DropDownList("View Side",5,ylay.x0,5+100, 3 ) ) )
        ->addItem("Top")
        ->addItem("Botton")
        ->addItem("Front")
        ->addItem("Back")
        ->addItem("Left")
        ->addItem("Right")
        ->setCommand( [&](GUIAbstractPanel* me_){ 
            DropDownList& me = *(DropDownList*)me_;
            printf( "old cam.qrot(%g,%g,%g,%g) -> %s \n", cam.qrot().x,cam.qrot().y,cam.qrot().z,cam.qrot().w, me.labels[me.iSelected].c_str()  );
            switch(me.iSelected){
                case 0: cam.setQrot(qTop);    break;
                case 1: cam.setQrot(qBottom); break;
                case 2: cam.setQrot(qFront);  break;
                case 3: cam.setQrot(qBack);   break;
                case 4: cam.setQrot(qLeft);   break;
                case 5: cam.setQrot(qRight);  break;
            }
            printf( "->new cam.qrot(%g,%g,%g,%g) \n", cam.qrot().x,cam.qrot().y,cam.qrot().z,cam.qrot().w );
            printf( "cam: aspect %g zoom %g \n", cam.aspect(), cam.zoom());
            printMat((Mat3d)cam.rotMat());
            }
        );

    GUIPanel*     p   =0;
    MultiPanel*   mp  =0;
    CheckBoxList* chk =0;

    // ------ MultiPanel(    "Edit"   )

    ylay.step( 1 ); ylay.step( 2 );
    mp= new MultiPanel( "Edit", gx.x0, ylay.x0, gx.x1, 0,-13); gui.addPanel( mp ); panel_Edit=mp;
    //GUIPanel* addPanel( const std::string& caption, Vec3d vals{min,max,val}, bool isSlider, bool isButton, bool isInt, bool viewVal, bool bCmdOnSlider );
    // mp->addPanel( "Sel.All", {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ 
    //     if(bViewBuilder){ W->builder.selection.clear(); for(int i=0; i<W->builder.atoms.size(); i++)W->builder.selection.insert(i); return 0; }
    //     else            { W->        selection.clear(); for(int i=0; i<W->nbmol.natoms; i++)W->selection.push_back(i);              return 0; }
    // };
    // mp->addPanel( "Sel.Inv", {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ 
    //     if(bViewBuilder){ std::unordered_set<int> s(W->builder.selection.begin(),W->builder.selection.end()); W->builder.selection.clear(); for(int i=0; i<W->builder.atoms.size(); i++) if( !s.contains(i) )W->builder.selection.insert(i);  return 0;  }
    //     else            { std::unordered_set<int> s(W->selection.begin(),        W->selection.end());         W->selection.clear();         for(int i=0; i<W->nbmol.natoms;         i++) if( !s.contains(i) )W->selection.push_back(i);       return 0;  }
    // };

    mp->addPanel( "print.nonB",  {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->ffl.print_nonbonded();   return 0; };   // 1
    mp->addPanel( "print.Aconf", {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->builder.printAtomConfs(); return 0; };  // 2
    mp->addPanel( "Sel.All", {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ if(bViewBuilder){ W->builder.selectAll();     }else{ W->selectAll();    } return 0; };  // 3
    mp->addPanel( "Sel.Inv", {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ if(bViewBuilder){ W->builder.selectInverse(); }else{ W->selectInverse();} return 0; };  // 4
    mp->addPanel( "Sel.Cap", {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->builder.selectCaping(); for(int ia: W->builder.selection) W->selection.push_back(ia); return 0; }; // 5
    mp->addPanel( "Add.CapHs",{0.0,1.0, 0.0}, 0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){  // 6
        //printf("====== AtomConf before Add.CapHs \n"); W->builder.printAtomConfs();
        bBuilderChanged = W->builder.addAllCapsByPi( W->params.getAtomType("H") ) > 0; 
        //printf("====== AtomConf After Add.CapHs \n"); 
        //W->builder.printAtomConfs();
        return 0; }; 
    //mp->addPanel( "rot3a"  , {0.0,1.0, 0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->rot3a();            return 0; };
    mp->addPanel( "toCOG"  , {-3.0,3.0,0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ if(bViewBuilder){W->selectionFromBuilder();} W->center(true);         if(bViewBuilder){W->updateBuilderFromFF();} return 0; }; // 7
    mp->addPanel( "toPCAxy", {-3.0,3.0,0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ if(bViewBuilder){W->selectionFromBuilder();} W->alignToAxis({2,1,0}); if(bViewBuilder){W->updateBuilderFromFF();} return 0; }; // 5
    p=mp->addPanel( "save.xyz",{-3.0,3.0,0.0},  0,1,0,0,0 );p->command = [&](GUIAbstractPanel* p){ const char* fname = ((GUIPanel*)p)->inputText.c_str(); if(bViewBuilder){ W->builder.save2xyz(fname);}else{W->saveXYZ(fname);} return 0; }; p->inputText="out.xyz";  // 9
    p=mp->addPanel( "save.mol:",{-3.0,3.0,0.0},  0,1,0,0,0 );p->command = [&](GUIAbstractPanel* p){ const char* fname = ((GUIPanel*)p)->inputText.c_str(); W->updateBuilderFromFF(); W->builder.saveMol(fname); return 0; }; p->inputText="out.mol";  // 10

    //mp->addPanel( "VdwRim", {1.0,3.0,1.5},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){  double R=((GUIPanel*)p)->value; dipoleMap.points_along_rim( R, {5.0,0.0,0.0}, Vec2d{0.0,0.1} );  bDipoleMap=true;  return 0; };
    //mp->addPanel( "VdwRim", {1.0,3.0,1.5},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){  double R=((GUIPanel*)p)->value; dipoleMap.prepareRim( 5, {1.0,2.0}, {0.4,0.6}, {0.0,-7.0,0.0}, {1.0,0.0} );  bDipoleMap=true;  return 0; };
    //mp->addPanel( "VdwRim", {1.0,3.0,1.5},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){  double R=((GUIPanel*)p)->value;   double Rs[]{R,R+0.1,R+0.2,R+0.3,R+0.5};  dipoleMap.prepareRim2( 5, Rs, 0.3, {0.0,-7.0,0.0}, {1.0,0.0} );  bDipoleMap=true;  return 0; };
    mp->addPanel( "VdwRim", {1.0,3.0,1.5},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){  // 11
        //double Rs[]{0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.5,3.0,3.5,4.5,6.0,8.0,10.0};  
        double Rs[20]{0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 2.0, 2.5, 3.0, 3.5, 4.5, 6.0, 8.0, 10.0, 14.0, 20.0 };  
        //dipoleMap.prepareRim2( 15, Rs, 0.3, {0.0,-7.0,0.0}, {1.0,0.0} ); 
        W->hideEPairs();
        dipoleMap.genAndSample( 20, Rs, 0.3, {0.0,-7.0,0.0}, {1.0,0.0} );
        //dipoleMap.genAndSample( 1, Rs, 0.05, {0.0,-7.0,0.0}, {1.0,0.0} );
        dipoleMap.saveToXyz( "dipoleMap.xyz" );
        bDipoleMap=true;  return 0;      
    };
    //int npan = mp->subs.size();

    auto lamb_scanSurf = [&](GUIAbstractPanel* p){ 
        if(W->ipicked<0)return; 
        //particle_REQ = 0.0;
        int npan =6;
        double sc = pow( 10., panel_Edit->subs[npan+0]->value );
        int iview = panel_Edit->subs[npan+1]->value;
        int iscan = panel_Edit->subs[npan+2]->value;
        printf( "lamb_scanSurf() npan=%i iview=%i iscan=%i sc=%g SurfFF_view.range(%g,%g)\n", npan, iview, iscan, sc,  panel_Edit->subs[npan+1]->vmin, panel_Edit->subs[npan+1]->vmax  );
        scanSurfFF( 100, apos[W->ipicked]+Vec3d{0.0,0.0,-5.0}, apos[W->ipicked]+Vec3d{0.0,0.0,+5.0}, W->ffl.REQs[W->ipicked], iscan, iview, sc, Vec3dX );
        //scanSurfFF( 100, apos[W->ipicked]-Vec3d{0.0,0.0,5.0}, apos[W->ipicked]-Vec3d{0.0,0.0,1.0}, particle_REQ, iscan, iview, sc );
    };
    mp->addPanel( "scanSurfFF",  {-3.0,3.0,0.0}, 1,1,0,1,1 )->command = lamb_scanSurf;   // 12
    mp->addPanel( "SurfFF_view", {-0.01,1.01,0}, 1,1,1,1,1 )->command = lamb_scanSurf;   // 13
    mp->addPanel( "SurfFF_scan", {-0.01,1.01,0}, 1,1,1,1,1 )->command = lamb_scanSurf;   // 14
    ylay.step( (mp->nsubs+1)*2 ); ylay.step( 2 );

    // mp= new MultiPanel( "Run", gx.x0, ylay.x0, gx.x1, 0,-2); gui.addPanel( mp ); //panel_NonBondPlot=mp;
    // mp->addPanel( "NonBond"  , {-3.0,3.0,0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->bNonBonded=!W->bNonBonded;         return 0; };
    // mp->addPanel( "NonBondNG", {-3.0,3.0,0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->bNonBondNeighs=!W->bNonBondNeighs; return 0; };
    // mp->addPanel( "Grid",      {-3.0,3.0,0.0},  0,1,0,0,0 )->command = [&](GUIAbstractPanel* p){ W->bGridFF=!W->bGridFF;               return 0; };
    
    // ------ CheckBoxList( "Run" )

    chk = new CheckBoxList( gx.x0, ylay.x0, gx.x1 );   gui.addPanel( chk );  chk->caption ="Run"; //chk->bgColor = 0xFFE0E0E0;
    chk->addBox( "NonBond"   , &W->bNonBonded     );
    chk->addBox( "NonBondNG" , &W->bNonBondNeighs );
    chk->addBox( "GridFF"    , &W->bGridFF        );
    chk->addBox( "tricubic"  , &W->bTricubic      );
    ylay.step( (chk->boxes.size()+1)*2 ); ylay.step( 2 );

    // ------ MultiPanel(   BondLenghs   )

    // ==== bond length calculation
    auto blFunc = [&](GUIAbstractPanel* p){
        //printf("MolGUI::[]blFunc() %s\n", BondLengh_types->inputText.c_str() );
        Vec2d lrange{ BondLengh_min->value, BondLengh_max->value };
        Vec2i typs = W->params.parseBondAtomTypes( BondLengh_types->inputText.c_str(), true );
        makeBondColoring( typs, lrange, bL0s, true );
        bViewBuilder     = false;
        bViewBondLenghts = true;
        bViewBondLabels  = false;
        bViewAtomLabels  = false;
    };

    mp= new MultiPanel( "BondLenghs", gx.x0, ylay.x0, gx.x1, 0,-3); gui.addPanel( mp );
    p=mp->addPanel( "types:",  {1.0,3.0,1.5 },  0,1,0,0,0 );p->command=blFunc; BondLengh_types=p; p->inputText="Si-Si";
    //p=mp->addPanel( "min.BL:", {2.2,2.6,2.30},  1,1,0,1,1 );p->command=blFunc; BondLengh_min=p;
    //p=mp->addPanel( "max.BL:", {2.2,2.6,2.40},  1,1,0,1,1 );p->command=blFunc; BondLengh_max=p;
    p=mp->addPanel( "min.BL:", {2.2,2.6,2.32},  1,1,0,1,1 );p->command=blFunc; BondLengh_min=p;
    p=mp->addPanel( "max.BL:", {2.2,2.6,2.43},  1,1,0,1,1 );p->command=blFunc; BondLengh_max=p;
    printf( "MolGUI::initWiggets() WorldVersion=%i \n", W->getMolWorldVersion() );
     ylay.step((mp->nsubs+1)*2);
    //exit(0);

    // ------ GUIPanel(   Mol. Orb.   )

    if( W->getMolWorldVersion() & (int)MolWorldVersion::QM ){ 
        // --- Selection of orbital to plot
        ylay.step( 2 );
        panel_iMO = ((GUIPanel*)gui.addPanel( new GUIPanel( "Mol. Orb.", 5,ylay.x0,5+100,ylay.x1, true, true, true ) ) );
        panel_iMO->setRange(-5.0,5.0);
        panel_iMO->setValue(0.0);
        panel_iMO->setCommand( [&](GUIAbstractPanel* p){ 
            which_MO = ((GUIPanel *)p)->value;
            int iHOMO = W->getHOMO(); printf( "plot HOMO+%i (HOMO=eig#%i) \n", iHOMO+which_MO, iHOMO );
            renderOrbital( iHOMO + which_MO );
        return 0; });
        ylay.step( 2 );
    }

    // ------ GUIPanel(   Mol. Orb.   )

    if( W->getMolWorldVersion() & (int)MolWorldVersion::GPU ){ 
        // --- Selection of orbital to plot
        ylay.step( 2 );
        panel_AFM = ((GUIPanel*)gui.addPanel( new GUIPanel( "AFM iz", 5,ylay.x0,5+100,ylay.x1, true, true, true ) ) );
        panel_AFM->setRange(0,20.0);
        panel_AFM->setValue(2);
        panel_AFM->setCommand( [&](GUIAbstractPanel* p){ 
            afm_iz = ((GUIPanel *)p)->value;
            makeAFM( afm_iz );
        return 0; });
    }

    MolGUI::nonBondGUI();
}

void MolGUI::nonBondGUI(){
    GUI_stepper gx(100,6);
    //bDrawNonBond = true;
    // ---- NonBond plot Options
    CheckBoxList* chk = new CheckBoxList( gx.x0, 10, gx.x1 );   panel_NBPlot=chk;
    gui.addPanel( chk );
    //panel_Plot = chk; 
    chk->caption = "NBPlot"; chk->bgColor = 0xFFE0E0E0;
    chk->addBox( "minima", &bDrawParticles    );
    chk->addBox( "lines ", &bDrawNonBondLines );
    chk->addBox( "gridXY", &bDrawNonBondGrid  );
    chk->addBox( "hideEp", &hideEp            );
    //chk->addBox( "findHb", &findHb            );

    
    //chk->addBox( "grid", &plot1.bGrid );
    //chk->addBox( "axes", &plot1.bAxes );

    GUIPanel* p=0;
    MultiPanel* mp=0;

    // ----- 1D plot option
    gx.step( 2 ); gx.step( 8+10 );
    mp= new MultiPanel( "PlotNonBond", gx.x0, 10, gx.x1, 0,-7); gui.addPanel( mp ); panel_NonBondPlot=mp;
    //GUIPanel* addPanel( const std::string& caption, Vec3d vals{min,max,val}, bool isSlider, bool isButton, bool isInt, bool viewVal, bool bCmdOnSlider );
    mp->addPanel( "Mode  : ", {0.0,2.0, 0.0},  1,0,1,1,0 );   // ToDo: perhaps it would be better to use bit-mask
    mp->addPanel( "Ezoom : ", {-3.0,3.0,-1.0}, 1,0,0,1,0 );
    mp->addPanel( "Rplot : ", {0.0,10.0,5.0},  1,0,0,1,0 );
    mp->addPanel( "dstep : ", {0.02,0.5,0.1},  1,0,0,1,0 );
    mp->addPanel( "Rdamp : ", {-2.0,2.0,0.1},  1,0,0,1,0 );
    mp->addPanel( "Rcut  : ", {-2.0,2.0,10.1}, 1,0,0,1,0 );
    mp->addPanel( "findHb:",  {1.5,5.0,4.0},   1,1,0,1,1 )->command = [&](GUIAbstractPanel* p){ double Rc = ((GUIPanel*)p)->value; W->Hbonds.clear(); W->findHbonds_PBC( Rc, 0.01, 30.0*deg2rad ); };

    // ----- TestAtom Params
    gx.step( 2 ); gx.step( 8+10 );
    mp = new MultiPanel( "GridXY", gx.x0, 10, gx.x1, 0,-4);   gui.addPanel( mp );    panel_GridXY=mp;
    mp->addPanel( "n     : ", {  10,200, 150  }, 1,0,1,1,0 );
    mp->addPanel( "size  : ", { 2.0,20.0,15.0  }, 1,0,0,1,0 );
    mp->addPanel( "vmin  : ", {-6.0,6.0,-0.0  }, 1,0,0,1,0 );
    mp->addPanel( "z_cut : ", {-5.0,5.0, 0.0  }, 1,0,0,1,0 );

    // ----- TestAtom Params
    gx.step( 2 ); gx.step( 8+10 );
    mp = new MultiPanel( "TestType", gx.x0, 10, gx.x1, 0,-4);   gui.addPanel( mp );    panel_TestType=mp;
    mp->addPanel( "RvdW  : ", { 0.0,2.50,1.5    },  1,0,0,1,0 );
    mp->addPanel( "EvdW  : ", { 0.0,0.02,0.0006808},1,0,0,1,0 );
    mp->addPanel( "Charge: ", {-0.5,0.5,+0.2    },  1,0,0,1,0 );
    mp->addPanel( "Hbond : ", {-1.0,1.0, 0.0    },  1,0,0,1,0 );

    // // ----- PickedAtom params
    gx.step( 2 ); gx.step( 8+10 );
    mp = new MultiPanel( "PickedType", gx.x0, 10, gx.x1, 0,-4 ); gui.addPanel( mp );  panel_PickedType=mp;
    mp->addPanel( "RvdW  : ", { 0.0,2.50,1.5    },  1,0,0,1,0 );
    mp->addPanel( "EvdW  : ", { 0.0,0.02,0.0006808},1,0,0,1,0 );
    mp->addPanel( "Charge: ", {-0.5,0.5, 0.0    },  1,0,0,1,0 );
    mp->addPanel( "Hbond : ", {-1.0,1.0, 0.0    },  1,0,0,1,0 );

}

void MolGUI::plotNonBondLines(){
    //GUIPanel* p=0;
    MultiPanel* mp=0;
    mp=panel_TestType;    
    Quat4d REQtest{ 
        mp->subs[0]->value, 
        sqrt(mp->subs[1]->value), 
        mp->subs[2]->value, 
        mp->subs[3]->value 
    };
    REQtest.w*=REQtest.y;  // Hbond = %H * sqrt(EvdW_ii)
    double dstep, Rplot, Ezoom, Rdamp, Rcut;
    mp=panel_NonBondPlot; Ezoom=pow(10.,mp->subs[1]->value); Rplot=mp->subs[2]->value; dstep=mp->subs[3]->value;  Rdamp = mp->subs[4]->value; Rcut = mp->subs[5]->value;
    
    // --- plot lines along each e-pair in the molecule
    // -- display list - delete if exist
    int epair_element = W->params.getElementType("E");
    //double R0 = 1.5;

    for(int i=0; i<W->nbmol.natoms; i++){
        int ityp = W->nbmol.atypes[i];
        AtomType& t = W->params.atypes[ityp];
        //printf( "plotNonBond[%i] %s t.element=%i epair_element=%i \n", i, t.name, t.element,  epair_element );
        if(t.element!=epair_element) continue;
        //printf( "plotNonBond[%i] %s | %s \n", i, t.name, "IS E-PAIR" );

        int j = W->ffl.neighs[i].x;
        Vec3d pe    = W->nbmol.apos[i];
        Vec3d pa    = W->nbmol.apos[j];
        Quat4d REQa = W->nbmol.REQs[j];
        Vec3d d = pe-pa;
        double r = d.normalize();

        // double plotNonBondLine( const NBFF& ff, Quat4d REQi, double Rdamp, Vec3d p1, Vec3d p2, int n, Vec3d up=Vec3dZ, bool bForce ){
        double R0ij = REQa.x + REQtest.x;
        Vec3d p0 = pa + d*R0ij*0.5;
        Vec3d p1 = p0 + d*Rplot;
        int n = (int)(Rplot/dstep);

        int jtyp = W->nbmol.atypes[j];
        //printf( "plotNonBond[%i]: jtyp=%i j=%i n=%i \n", i, jtyp,j, n );
        //printf( "plotNonBond[%i]: %s(%i) n=%i \n", i, W->params.atypes[jtyp].name,j, n );

        plotNonBondLine( W->nbmol, REQtest, Rdamp, p0, p1, n, Vec3dZ*Ezoom, false );
        
    }
}

void MolGUI::plotNonBondGrid(){
    //GUIPanel* p=0;
    MultiPanel* mp=0;
    mp=panel_TestType;    
    Quat4d REQtest{ 
        mp->subs[0]->value, 
        sqrt(mp->subs[1]->value), 
        mp->subs[2]->value, 
        mp->subs[3]->value 
    };
    REQtest.w*=REQtest.y;  // Hbond = %H * sqrt(EvdW_ii)
    mp=panel_NonBondPlot; 
    double Ezoom = pow(10.,mp->subs[1]->value); 
    double Rplot = mp->subs[2]->value; 
    double dstep = mp->subs[3]->value;  
    double Rdamp = mp->subs[4]->value; 
    double Rcut  = mp->subs[5]->value;

    mp=panel_GridXY; 
    int    n    = mp->subs[0]->getIntVal();
    double sz   = mp->subs[1]->value; 
    double vmax = pow(10.,mp->subs[2]->value);
    double z_cut= mp->subs[3]->value;
    
    Vec2i ns{n,n};
    Vec3d p0{-sz    ,-sz   ,0.0};
    Vec3d a { sz*2.0,   0.0,0.0};
    Vec3d b { 0.0   ,sz*2.0,0.0};
    double* Egrid=0;

    int nEp=0, etyp=0;
    //printf( "plotNonBondGrid() hideEp %i \n", hideEp );
    if(hideEp){
        // //printf( "plotNonBondGrid() removing EPairs %i \n" );
        // etyp = W->params.getAtomType("E");
        // W->ffl.chargeToEpairs( -W->QEpair, etyp );
        // W->selectByType( W->params.getElementType("E"), true );
        // nEp = W->selection.size();
        // W->nbmol.natoms -= nEp;
        W->hideEPairs();
    }
    Vec2d Erange = evalNonBondGrid2D( W->nbmol, REQtest, Rdamp, ns, Egrid, p0, a,b );
    if(hideEp){
        //printf( "plotNonBondGrid() adding EPairs %i \n" );
        // W->nbmol.natoms += nEp;
        // W->ffl.chargeToEpairs( W->QEpair, etyp );
        W->unHideEPairs();
    }
    
    double lvmax = 10.0;

    //printf( "Erange(%g,%g) vlim(%g,%g)\n", Erange.x, Erange.y, -vmax, vmax );
    //Draw3D::drawScalarGrid( ns, p0,b*(1./ns.b),a*(1./ns.a), Egrid, -vmax, vmax, Draw::colors_RWB ); //  const uint32_t * colors, int ncol );
    opengl1renderer.lineWidth(0.25);
    Draw3D::drawScalarGridLines( ns, p0, b*(1./ns.b), a*(1./ns.a), Vec3dZ, Egrid, lvmax/vmax, Vec2d{-vmax,vmax} );
    opengl1renderer.lineWidth(1.0);
    // saddly text drawing in drawAxis3D is not working when baked in display list, so we have to draw it on the fly
    //opengl1renderer.color3f( 0.0, 0.0, 0.0 );
    //int nEtick =10;
    //Draw3D::drawAxis3D( {ns.x,ns.y,nEtick}, p0+(Vec3dZero*(-lvmax/nEtick)), {sz/ns.a,sz/ns.b,}, {sz*-0.5/ns.a,sz*-0.5/ns.b,nEtick*-0.5*vmax}, {sz/ns.a,sz/ns.b,nEtick*vmax}, fontTex, textSize, "%6.6f" );
    
    delete [] Egrid;
}

void MolGUI::plotNonBondGridAxis(){
    MultiPanel* mp=0;
    mp=panel_GridXY; 
    int    n    = mp->subs[0]->getIntVal();
    double sz   = mp->subs[1]->value; 
    double vmax = pow(10.,mp->subs[2]->value);
    double lvmax = 10.0;
    int nEtick =10;
    n/=10;
    Vec2i ns{n,n};
    Vec3d p0{-sz    ,-sz   ,-lvmax};
    Vec3d a { sz*2.0,   0.0,0.0};
    Vec3d b { 0.0   ,sz*2.0,0.0};

    //double to_meV = 1000.0;
    double to_meV = 0.0;
    Draw3D::drawAxis3D( {ns.x,ns.y,nEtick}, p0, {sz*2/ns.a,sz*2/ns.b,2*lvmax/nEtick}, {-sz,-sz,-vmax*to_meV }, {2*sz/ns.a,2*sz/ns.b,2*vmax*to_meV/nEtick }, fontTex3D, 0.1, textSize*0.7, "%.2f" );
}

void MolGUI::tryPlotNonBond(){
    // --- check if parameters changed
    bool bChanged =  (panel_GridXY->clearChanged()>=0) || (panel_NonBondPlot->clearChanged()>=0) || (panel_TestType->clearChanged()>=0) || (panel_PickedType->clearChanged()>=0) || (panel_NBPlot->clearChanged()>=0);
    //printf( " MolGUI::plotNonBond() bChanged=%i  ogl_nonBond=%i \n", bChanged, ogl_nonBond );
    bool ogl0     = (ogl_nonBond<=0);
    if( !( bChanged || ogl0 ) ){ return; };
    if( !ogl0 ){ opengl1renderer.deleteLists(ogl_nonBond,1); }
    //printf( " MolGUI::plotNonBond() bChanged=%i  ogl_nonBond=%i \n", bChanged, ogl_nonBond );
    ogl_nonBond = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_nonBond, GL_COMPILE);
    if( bDrawNonBondLines ){ plotNonBondLines(); };
    if( bDrawNonBondGrid  ){ plotNonBondGrid();  };
    opengl1renderer.endList();
}

void MolGUI::relaxNonBondParticles( double dt, double Fconv, int niter){
    if( panel_NBPlot->clearChanged()<0 ) return;
    //printf( "relaxNonBondParticles() dt=%g Fconv=%g niter=%i \n", dt, Fconv, niter );
    double F2conv = Fconv*Fconv;
    MultiPanel* mp=0;
    mp=panel_TestType;    
    Quat4d REQtest{ 
        mp->subs[0]->value, 
        sqrt(mp->subs[1]->value), 
        mp->subs[2]->value, 
        mp->subs[3]->value 
    };
    // ----------- SELECT E-PAIRS and add hydrogen particles next to it
    particles.clear();
    particlePivots.clear();
    int epair_element = W->params.getElementType("E");
    for(int i=0; i<W->nbmol.natoms; i++){
        int ityp = W->nbmol.atypes[i];
        AtomType& t = W->params.atypes[ityp];
        if(t.element!=epair_element) continue;
        int j = W->ffl.neighs[i].x;
        Vec3d pe    = W->nbmol.apos[i];
        Vec3d pa    = W->nbmol.apos[j];
        Quat4d REQa = W->nbmol.REQs[j];
        Vec3d d = pe-pa;
        double r = d.normalize();
        double R0ij = REQa.x + REQtest.x;
        Quat4d p; p.f = pa + d*R0ij;
        particles.push_back( p );
        particlePivots.push_back(j); // atom to which the particle is attached
    }
    // ------------ RELAXATION
    if(hideEp){ W->hideEPairs(); }
    //std::vector<Vec3d> vpos( particles.size() );
    bool bTrj = true;
    //int ogl_trj = 0;
    if(bTrj){
        if(ogl_trj>0) opengl1renderer.deleteLists(ogl_trj,1);
        ogl_trj = opengl1renderer.genLists(1);
        opengl1renderer.newList(ogl_trj, GL_COMPILE);
    }
    for( int i=0; i<particles.size(); i++ ){
        Vec3d p    = particles[i].f;
        //Vec3d v    = vpos[i];
        Vec3d v    = Vec3dZero;
        Quat4d fe;
        if(bTrj){ Draw3D::drawPointCross( p, 0.2 ); opengl1renderer.begin(GL_LINE_STRIP); }
        for(int iter=0; iter<niter; iter++){
            if(bTrj){ opengl1renderer.vertex3f(p.x,p.y,p.z); }
            fe = W->nbmol.evalLJQs( p, REQtest, W->ffl.Rdamp );
            fe.f.mul( -1.0 );
            if( fe.norm2() < F2conv ) break;
            double cvf = v.dot(fe.f);
            //printf( "relaxNonBondParticles[%i] iter %i E=%g |f|=%g cvf=%4.2f \n", i, iter, fe.e, fe.norm2(), cvf/sqrt(fe.norm2()*v.norm2()) );
            if( cvf<0 ){ v.set(0.0); }
            v += fe.f * dt;
            p += v    * dt;
        }
        if(bTrj){ opengl1renderer.end(); }
        //vpos[i]        = v;
        particles[i].f = p;
        particles[i].e = fe.e;
    }
    if(bTrj){ opengl1renderer.endList(); }
    if(hideEp){ W->unHideEPairs(); }
}

void MolGUI::drawParticles(){
    //printf( "drawParticles() particles.size()=%i \n", particles.size() );
    // int ityp       = W->params.getAtomType("H");
    // AtomType& atyp = W->params.atypes[ityp];
    // double r = ((atyp.RvdW-mm_Rsub)*mm_Rsc);
    Mat3d m  = Mat3dIdentity*0.2;
    opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.shadeModel(GL_SMOOTH);
    for( int i=0; i<particles.size(); i++ ){
        Quat4d& p = particles[i];
        Draw3D::drawSphere((Vec3f)p.f, 1, COLOR_GREEN);
    }
    opengl1renderer.color3f( 0.0, 0.0, 0.0 );
    opengl1renderer.begin(GL_TRIANGLES);
    double to_meV = 1000.0;
    for( int i=0; i<particles.size(); i++ ){
        Quat4d& p = particles[i];    Draw3D::drawDouble( p.f, p.e*to_meV, fontTex3D, textSize, "%.3fmeV" );
        
    }
    // ---- Pivots
    for( int i=0; i<particles.size(); i++ ){
        Vec3d& p = particles[i].f;
        int    j = particlePivots[i];
        Vec3d& a = W->nbmol.apos [j];
        Draw3D::drawLine( p, a, COLOR_GREEN );
    }
    opengl1renderer.color3f( 0.0, 0.0, 0.0 );
    for( int i=0; i<particles.size(); i++ ){
        Vec3d& p = particles[i].f;  
        int    j = particlePivots[i];
        Vec3d& a = W->nbmol.apos [j];
        double l = (a-p).norm();
        Draw3D::drawDouble( (a+p)*0.5, l, fontTex3D, textSize, "%.2fA" );   
    }
}

void MolGUI::drawDipoleMap(){
    MultiPanel* mp=0;
    //mp=panel_TestType;    
    // Quat4d REQtest{ 
    //     mp->subs[0]->value, 
    //     sqrt(mp->subs[1]->value), 
    //     mp->subs[2]->value, 
    //     mp->subs[3]->value 
    // };
    //REQtest.w*=REQtest.y;  // Hbond = %H * sqrt(EvdW_ii)
    // double dstep, Rplot, Ezoom, Rdamp, Rcut;
    mp=panel_NonBondPlot; //Ezoom=pow(10.,mp->subs[1]->value);
    double Ezoom=1/pow(10.,mp->subs[1]->value);
    //printf( "drawDipoleMap() Ezoom=%g  val=%g \n", Ezoom, mp->subs[1]->value ); 


    Draw3D::drawPointCross( {5.0,0.0,0.0}, 0.2 );
    Draw3D::drawPointCross( dipoleMap.particles[0], 0.1 );
    for(int i=0; i<dipoleMap.particles.size(); i++ ){
        Draw3D::drawPointCross( dipoleMap.particles[i], 0.05, COLOR_BLACK );
        Draw3D::drawLine      ( dipoleMap.particles[i], dipoleMap.particles2[i], COLOR_BLUE );
    }
    
    //int nphi = dipoleMap.nphi;
    //int nr   = dipoleMap.FE.size() / nphi;
    //printf( "drawDipoleMap() nphi=%i nr=%i FE.size(%i) particles.size(%i) \n", nphi, nr, dipoleMap.FE.size(), dipoleMap.particles.size(), nr*nphi );

    // for(int ip=0; ip<nphi; ip++ ){
    //     opengl1renderer.begin(GL_LINE_STRIP);
    //     for(int ir=0; ir<nr; ir++ ){
    //         int i = ip + ir*nphi;
    //         //printf( "drawDipoleMap()[%i|ip=%i,ir=%i]\n", i, ip, ir );
    //         Vec3d p = dipoleMap.particles[i];
    //         p.z    += dipoleMap.FE[i].e/Ezoom;
    //         opengl1renderer.vertex3f( p.x, p.y, p.z );
    //     }
    //     opengl1renderer.end();
    // }    
    // for(int ir=0; ir<nr; ir++ ){
    //    opengl1renderer.begin(GL_LINE_STRIP);
    //     for(int ip=0; ip<nphi; ip++ ){
    //         int i = ip + ir*nphi;
    //         Vec3d p = dipoleMap.particles[i];
    //         p.z    += dipoleMap.FE[i].e/Ezoom;
    //         opengl1renderer.vertex3f( p.x, p.y, p.z );
    //     }
    //     opengl1renderer.end();
    // }

    opengl1renderer.lineWidth(1.0);
    opengl1renderer.color3f( 0.0, 0.7, 0.0 ); drawDipoleMapGrid( dipoleMap, Vec2d{1.0,1.0}*Ezoom, true, true );
    opengl1renderer.lineWidth(0.25);
    opengl1renderer.color3f( 1.0, 0.5, 0.0 ); drawDipoleMapGrid( dipoleMap, Vec2d{1.0,0.0}*Ezoom, true, true );
    opengl1renderer.color3f( 0.0, 0.7, 1.0 ); drawDipoleMapGrid( dipoleMap, Vec2d{0.0,1.0}*Ezoom, true, true );
    opengl1renderer.lineWidth(1.0);
}

MolGUI::MolGUI( int& id, int WIDTH_, int HEIGHT_, MolWorld_sp3* W_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    printf( "MolGUI::MolGUI() \n" );
    if( W_ ) W=W_;
    long T0=getCPUticks();
    float nseconds = 0.1;
    SDL_Delay( (int)(1000*0.1) );
    tick2second = nseconds/(getCPUticks()-T0);
    printf( "CPU speed calibration: tick=%g [s] ( %g GHz)\n", tick2second, 1.0e-9/tick2second );
    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" ); GUI_fontTex = fontTex;
    fontTex3D = makeTexture    ( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );
    actions.vec.resize( 256 );
    initGUI();
}

void MolGUI::initGUI(){

    initMemberOffsets();

    initCommands();

    // ---- Graphics setup
    //Console::init( int lineLength=256, SDL_Window* window_=0 ){
    console.init( 256, window );
    console.callback = [&](const char* s){ printf( "console.callback(%s)\n", s ); return 0; };
    console.fontTex = fontTex;

    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //opengl1renderer.lightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//opengl1renderer.lightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	opengl1renderer.lightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );
    // ---- Gizmo
    cam.setPersp(false);
    gizmo.cam = &cam;
    //gizmo.bindPoints(W->ff.natoms, W->ff.apos      );
    //gizmo.bindEdges (W->ff.nbonds, W->ff.bond2atom );
    gizmo.pointSize = 0.5;
    //gizmo.iDebug    = 2;
    ruler.setStep( 1.5 * sqrt(3) );

}

void MolGUI::updateGUI(){
    gizmo.bindPoints( natoms, apos      );
    gizmo.bindEdges ( nbonds, bond2atom );
    if(panel_Frags){
        panel_Frags->labels.clear();
        for(int i=0; i<W->builder.frags.size(); i++){
            char s[16]; sprintf(s,"Frag_%02i", i );
            panel_Frags->addItem( s );
        }
    }
}

int MolGUI::clearGUI(int n){
    printf( "MolGUI::clearGUI(%i) \n", n );
    Qpanel      = 0;
    panel_Frags = 0;
    panel_iMO   = 0;
    return gui.clear( n );
}

void MolGUI::bindMolWorld( MolWorld_sp3* W_ ){
    //if(verbosity>0)
    printf("MolGUI::bindMolWorld() \n");
    W = W_;
    // if(W_!=0){ W = W_;}
    // if(W==0){
    //     W = new MolWorld_sp3(); 
    //     W->init();
    // }else{ 
    //     if( W->isInitialized==false){ W->init(); } 
    // }
    MolGUI::bindMolecule( W );

    W->getTitle(tmpstr); SDL_SetWindowTitle( window, tmpstr );
    //initGUI();
    initWiggets();
    updateGUI();
    //dipoleMap.ff = &(W->ffl);
    dipoleMap.ff     = &(W->nbmol);
    dipoleMap.params = &(W->params);
    //if(verbosity>0)printf("... MolGUI::bindMolWorld() DONE\n");
}

    bool MolGUI::visual_FF_test(){ // ==== MolGUI TESTS   Torsions ( Paolo vs Prokop optimized )
    
        if( frameCount==0){
            Vec3d ax = apos[1]-apos[0];
            apos[2].add_mul( ax, -0.5 );
            apos[3].add_mul( ax, -0.5 );
            for(int i=4; i<6; i++){ apos[i].mul(0.8); }
            for(int i=0; i<natoms; i++){ apos[i].mul(1.5); }

            apos[0].z += 0.5;
            apos[1].z -= 0.5;
        }

        double angle = 0.0;
        //if( keys[ SDL_SCANCODE_LEFTBRACKET  ] ){ angle=0.1; }
        //if( keys[ SDL_SCANCODE_RIGHTBRACKET ] ){ angle=-0.1; }

        if( keys[ SDL_SCANCODE_LEFTBRACKET  ] ){ apos[3].z -= 0.1; }
        if( keys[ SDL_SCANCODE_RIGHTBRACKET ] ){ apos[3].z += 0.1; }

        int sel[3]{1,4,5};
       
        Vec3d ax = apos[1]-apos[0];
        Vec3d p0 = apos[0];
        ax.normalize();
        Vec2d cs; cs.fromAngle( angle );
        for(int i=0; i<3; i++){
            apos[sel[i]].rotate_csa( cs.x, cs.y, ax, p0 );
        }

        // { // Check Dihedrals
        //     W->ffu.evalBonds();
        //     Vec3d fbak[4];
        //     double E,E_;
        //     E_=W->ffu.evalDihedral_Paolo( 0 );
        //     //E_=W->ffu.evalDihedral_Prokop_Old( 0 );
        //     for(int i=0; i<4; i++){ fbak[i]=W->ffu.fdih[i]; }
        //     //E=W->ffu.evalDihedral_Prokop_Old( 0 );
        //     E=W->ffu.evalDihedral_Prokop( 0 );
        //     printf( " Eerr %g |   E_ref %g E %g \n", E-E_, E_, E );
        //     checkVec3Matches( 4, W->ffu.fdih, fbak, "dih_fp", 1 );
        // }

        // { // Check Angles
        //     W->ffu.evalBonds();
        //     Vec3d fbak[4];
        //     double E,E_;
        //     E_=W->ffu.evalAngle_Paolo( 0 );
        //     for(int i=0; i<3; i++){ fbak[i]=W->ffu.fang[i]; }
        //     E=W->ffu.evalAngle_Prokop( 0 );
        //     printf( " Eerr %g |   E_ref %g E %g \n", E-E_, E_, E );
        //     checkVec3Matches( 3, W->ffu.fang, fbak, "dih_fp", 1 );
        // }

        // { // Check Inversions
        //     //W->ffu.printSizes();
        //     W->ffu.evalBonds();
        //     Vec3d fbak[4];
        //     double E,E_;
        //     E_=W->ffu.evalInversions_Paolo( 0 );
        //     for(int i=0; i<4; i++){ fbak[i]=W->ffu.fang[i]; }
        //     E=W->ffu.evalInversions_Prokop( 0 );
        //     printf( " Eerr %g |   E_ref %g E %g \n", E-E_, E_, E );
        //     //checkVec3Matches( 4, W->ffu.fang, fbak, "dih_fp", 1 );
        // }

        // { // check OMP
        //     Vec3d fapos[ W->ffu.natoms ];
        //     double E_ = W->ffu.eval();
        //     for(int i=0; i<W->ffu.natoms; i++){ fapos[i]=W->ffu.fapos[i]; }
        //     double E  = W->ffu.eval_omp();
        //     printf( " Eerr %g | E_ref %g E %g \n", E-E_, E_, E);
        //     bool b = checkVec3Matches( W->ffu.natoms, W->ffu.fapos, fapos, "fpos", 1 );
        //     if(b){ printf( "ffu.eval_omp() OK \n" ); }else{ printf( "ffu.eval_omp() FAILED \n" ); }
        //     exit(0);
        // }

        return false;
    }

//=================================================
//                   DRAW()
//=================================================

void MolGUI::draw(){
    GLES::active_camera = &cam;

    opengl1renderer.clearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    // Smooth lines : https://vitaliburkov.wordpress.com/2016/09/17/simple-and-fast-high-quality-antialiased-lines-with-opengl/
    //opengl1renderer.enable(GL_LINE_SMOOTH);
    opengl1renderer.enable(GL_BLEND);
    opengl1renderer.enable(GL_LIGHTING );
    opengl1renderer.enable(GL_DEPTH_TEST);

    //printf( "MolGUI::draw()[frameCount=%i] \n", frameCount );
    if(W->bLatScan){ lattice_scan( W->latscan_n.x, W->latscan_n.y, *W->latscan_dlvec ); quit(); }

    //if( (ogl_isosurf==0) && W->bGridFF ){ renderGridFF( subs_iso ); }
    //if( ogl_esp==0 ){ renderESP(); }

    if(frameCount==0){ cam.dpitch( M_PI ); }

    //debug_scanSurfFF( 100, {0.,0.,z0_scan}, {0.0,3.0,z0_scan}, 10.0 );

    // --- Mouse Interaction / Visualization
	//ray0 = (Vec3d)(  cam.rotMat().a*mouse_begin_x  +  cam.rotMat().b*mouse_begin_y  +  cam.pos );
    ray0 = mouseRay0();

    W->pick_hray = (Vec3d)cam.rotMat().c;
    W->pick_ray0 = (Vec3d)ray0;

    // if( (frameCount==0) && (W->bGridFF) ){ char fname[128]; 
    //     sprintf(fname,"gridFF_EFprofile_mod%i.log",   (int)W->gridFF.mode ); W->gridFF.getEFprofileToFile( fname, 1000, Vec3d{0.0,0.0,-10.0}, Vec3d{0.0,0.0,10.0}, W->ffl.REQs[0] ); 
    //     sprintf(fname,"gridFF_EFprofile_mod%i_x.log", (int)W->gridFF.mode ); W->gridFF.getEFprofileToFile( fname, 1000, Vec3d{-10.0,0.0,0.0}, Vec3d{10.0,0.0,0.0}, W->ffl.REQs[0] ); 
    //     sprintf(fname,"gridFF_EFprofile_mod%i_y.log", (int)W->gridFF.mode ); W->gridFF.getEFprofileToFile( fname, 1000, Vec3d{0.0,-10.0,0.0}, Vec3d{0.0,10.0,0.0}, W->ffl.REQs[0] );  
    //     //exit(0);
    // }  

    if(bRunRelax){ 
        bool bRelaxOld = W->bConverged;
        //printf( "MolGUI::draw().W->MDloop(%i) bUFF %i \n", perFrame, W->bUFF );    
        W->MDloop(perFrame); 
        // if( W->bConverged && !bRelaxOld ){  // it relaxed just now
        //     if(ogl_MO>0){ int iHOMO = W->getHOMO(); renderOrbital( iHOMO + which_MO );  }
        //     //updateGUI(); 
        // }
    }
    //if( bViewBuilder ){  W->updateBuilderFromFF(); }
    //if(bRunRelax){ W->relax( perFrame ); }

    Draw3D::drawPointCross( ray0, 0.1 );        // Mouse Cursor 
    //if(W->ipicked>=0) Draw3D::drawLine( W->ff.apos[W->ipicked], ray0); // Mouse Dragging Visualization
    if(W->ipicked>=0) Draw3D::drawLine( apos[W->ipicked], (Vec3d)ray0, COLOR_BLACK); // Mouse Dragging Visualization
    
    {   // draw mouse selection box;   ToDo:   for some reason the screen is upside-down
        //Vec3d ray0_ = ray0;            ray0_.y=-ray0_.y;
        //Vec3d ray0_start_=ray0_start;  ray0_start_.y=-ray0_start_.y;
        if(bDragging){
            drawMuseSelectionBox();
        }
    }

    //printf( "bViewSubstrate %i ogl_isosurf %i W->bGridFF %i \n", bViewSubstrate, ogl_isosurf, W->bGridFF );

    if(bViewCell){ 
        //Draw3D::drawTriclinicBox( W->builder.lvec, Vec3d{0.0,0.0,0.0}, Vec3d{1.0,1.0,1.0} ); 
        Draw3D::drawTriclinicBoxT( W->builder.lvec, Vec3d{0.0,0.0,0.0}, Vec3d{1.0,1.0,1.0} ); 
    }

    TESTfullscreenShader.begin();

    if( bViewSubstrate ){
        if( ( W->bGridFF )&&( ((int)(W->gridFF.mode))!=0) ){
            if( (ogl_isosurf.vertexCount()==0) ){ renderGridFF_new( subs_iso ); }
            viewSubstrate( {-5,10}, {-5,10}, &ogl_isosurf, W->gridFF.grid.cell.a, W->gridFF.grid.cell.b );
        }else{
            renderSurfAtoms(  Vec3i{1,1,0}, false );  
        }
    }

    // ----- Visualization of the Groups of Atoms
    if( W->bGroups ){  
        Quat4f *gpos=0,*gfw=0,*gup=0; int ng=W->getGroupPose( gpos, gfw, gup );
        if(gpos){
            for(int ig=0; ig<ng; ig++){
                Draw3D::drawVecInPos( gfw[ig].f*5., gpos[ig].f, COLOR_RED );
                Draw3D::drawVecInPos( gup[ig].f*5., gpos[ig].f, COLOR_GREEN );
                Draw3D::drawVecInPos( cross(gfw[ig].f,gup[ig].f), gpos[ig].f, COLOR_BLUE );
            }
        }
    }

    if( ogl_esp     ){ opengl1renderer.callList(ogl_esp);      }
    if( ogl_afm_trj ){ opengl1renderer.callList(ogl_afm_trj);  }
    if( ogl_afm     ){ opengl1renderer.callList(ogl_afm);      }

    if(bDrawNonBondGrid || bDrawNonBondLines){  tryPlotNonBond();  if( ogl_nonBond){ opengl1renderer.lineWidth(0.25); opengl1renderer.callList(ogl_nonBond); opengl1renderer.lineWidth(1.00); opengl1renderer.color3f(.0f,.0f,.0f); plotNonBondGridAxis(); } }
    if(bDrawParticles){ relaxNonBondParticles();  opengl1renderer.color3f(.0f,1.0f,.5f); if(ogl_trj){ opengl1renderer.callList(ogl_trj); } drawParticles(); }; 

    //Draw3D::drawMatInPos( W->debug_rot, W->ff.apos[0] ); // Debug

    //if(bDoQM)drawSystemQMMM();

    //if( ogl_MO && W->bConverged ){ 
    if( ogl_MO && (!bRunRelax) ){ 
        opengl1renderer.pushMatrix();
        //Vec3d c = W->builder.lvec.a*-0.5 + W->builder.lvec.b*-0.5 + W->builder.lvec.c*-0.5;
        //Vec3d c = Vec3dZero;
        //Vec3d c = W->gridFF.shift0;
        //W->cog = average( W->ffl.natoms, W->ffl.apos  );
        Vec3d pmin,pmax; bbox( pmin, pmax, W->ffl.natoms, W->ffl.apos, 0 ); W->cog=(pmin+pmax)*0.5;
        Vec3d c = W->cog + W->builder.lvec.a*-0.5 + W->builder.lvec.b*-0.5 + W->builder.lvec.c*-0.5;
        //printf( "ogl_MO c (%g,%g,%g) cog (%g,%g,%g) \n", c.x, c.y, c.z, W->cog.x, W->cog.y, W->cog.z );
        opengl1renderer.translatef( c.x, c.y, c.z );
            opengl1renderer.color3f(1.0,1.0,1.0); 
            opengl1renderer.callList(ogl_MO); 
        opengl1renderer.popMatrix();
    }

    // Draw the actual system ( molecules : atoms, bonds etc. )
    if(bDoMM){
        if( bViewBuilder ){ drawBuilder(); }   // Draw Builder 
        else if(W->builder.bPBC){              // Draw System with PBC 
            drawSystemShifts( W->npbc, W->pbc_shifts, W->ipbc0 );
            Draw3D::drawBBox( W->bbox.a, W->bbox.b );
        }else{ drawSystemSingle(); } // Draw System without PBC
    }

    TESTfullscreenShader.end();

    plotNonuniformGrid();

    if( bDipoleMap )drawDipoleMap();

    // Draw hydorgen bonds (if any)
    if( W->Hbonds.size() > 0 ){
        //int nfound = W->Hbonds.size();
        //printf( "findHb() Rc=%g nfound=%i \n", Rc, nfound );
        //if(ogl_Hbonds>0){ opengl1renderer.deleteLists(ogl_Hbonds,0); ogl_Hbonds=0; }
        //if( nfound>0 ){
        // opengl1renderer.newList( ogl_Hbonds, GL_COMPILE )
        opengl1renderer.lineWidth(5.0);
        //printf( "W->Hbonds.size(%i)\n", W->Hbonds.size() ); 
        opengl1renderer.begin(GL_LINES);
        opengl1renderer.color3f( 0.0,1.0,0.0 );
        for( Vec3i b: W->Hbonds ){
            //printf( "Hbonds(%i,%i,%i)\n", b.x,b.y,b.z ); 
            Draw3D::vertex( W->ffl.apos[b.a]                      );
            Draw3D::vertex( W->ffl.apos[b.b] + W->ffl.shifts[b.c] );
            Draw3D::vertex( W->ffl.apos[b.a] - W->ffl.shifts[b.c] );
            Draw3D::vertex( W->ffl.apos[b.b]  );
            //Draw3D::vertex( W->ffl.apos[b.b] + Vec3d{5.0,5.0,5.0}  );
        }
        opengl1renderer.end();
        opengl1renderer.lineWidth(1.0);
        // opengl1renderer.endList();
        //}
    }

    // Draw constraints (if any)
    if(constrs){
        // bond constrains
        for( DistConstr con : constrs->bonds ){ 
            Vec3d sh; W->builder.lvec.dot_to_T( con.shift, sh );
            Draw3D::drawLine( apos[con.ias.a],    apos[con.ias.b] + sh, {0, 0.7, 0} ); 
            Draw3D::drawLine( apos[con.ias.a]-sh, apos[con.ias.b]     , {0, 0.7, 0} ); 
        }
        // angle constrains
        opengl1renderer.color3f(0.0f,0.8f,0.8f);
        for( AngleConstr con : constrs->angles ){ 
            const Mat3d& lvec = W->builder.lvec;
            Vec3d ash = lvec.a*con.acell.a + lvec.b*con.acell.b + lvec.c*con.acell.c;
            Vec3d bsh = lvec.a*con.bcell.a + lvec.b*con.bcell.b + lvec.c*con.bcell.c;
            //Draw3D::drawTriangle( apos[con.ias.b] + ash,   apos[con.ias.a],   apos[con.ias.c] + bsh );
            //Draw3D::drawTriangle( apos[con.ias.b],   apos[con.ias.a],   apos[con.ias.c] );
            opengl1renderer.begin(GL_TRIANGLES);
            opengl1renderer.color3f(0.0f,1.0f,0.5f); Draw3D::vertex(apos[con.ias.b] + ash);
            opengl1renderer.color3f(0.0f,0.7f,0.7f); Draw3D::vertex(apos[con.ias.a]      );
            opengl1renderer.color3f(0.0f,0.5f,1.0f); Draw3D::vertex(apos[con.ias.c] + bsh);
            opengl1renderer.end();
        }
    }

    opengl1renderer.color3f(0.0f,0.5f,0.0f); showBonds();

    //visual_FF_test();

    if(W->ipicked>-1){ 
        Vec3d p,f;
        if(W->bUFF){
            p = W->ffu.apos[W->ipicked];
            f = getForceSpringRay( W->ffu.apos[W->ipicked], W->pick_hray, W->pick_ray0, W->Kpick ); 
        }else{
            p = W->ffl.apos[W->ipicked];
            f = getForceSpringRay( W->ffl.apos[W->ipicked], W->pick_hray, W->pick_ray0, W->Kpick ); 
        }
        Draw3D::drawVecInPos( f*-ForceViewScale, p, COLOR_RED );
    }

    if( ogl_surf_scan>0 ){ opengl1renderer.callList(ogl_surf_scan); }

    if(bDebug_scanSurfFF){ // --- GridFF debug_scanSurfFF()
        // ------ Deprecated ?
        Vec3d p0=W->gridFF.grid.pos0; 
        if(W->ipicked>-1){ p0.z = W->ffl.apos[W->ipicked].z; }{
            p0.z = W->ffl.apos[0].z;
        } 
        //printf( "p0.x %g \n", p0.z );
        int nx=W->gridFF.grid.n.x*3*5;
        int ny=W->gridFF.grid.n.y*3*5;
        Vec3d a = W->gridFF.grid.cell.a;
        Vec3d b = W->gridFF.grid.cell.b;
        //debug_scanSurfFF( nx, p0-a,       p0+a*2.       );
        //debug_scanSurfFF( ny, p0-b,       p0+b*2.       );
        debug_scanSurfFF( nx, p0-a+b*0.5, p0+a*2.+b*0.5 );
        debug_scanSurfFF( ny, p0-b+a*0.5, p0+b*2.+a*0.5 );
    }

    if(bViewBuilder){ opengl1renderer.color3f( 0.f,1.f,0.f ); for(int ia : W->builder.selection ){ Draw3D::drawSphereOctLines( 8, 0.5, W->builder.atoms[ia].pos ); } }
    else            { opengl1renderer.color3f( 0.f,1.f,0.f ); for(int ia : W->selection         ){ Draw3D::drawSphereOctLines( 8, 0.5, W->nbmol.apos[ia]        ); } }

    // --- Drawing Population of geometies overlay
    if(frameCount>=1){ 
        if(frameCount==1){
            nsys=W->getMultiSystemPointers( M_neighs, M_neighCell, M_apos, nvec); 
            //for(int i=0;i<nsys; i++){ printMSystem(i, 4, natoms, nvec ); }
        }
        //printf( "nsys %i \n", nsys );
        if(nsys>0){ 
            float dt = 2.*M_PI/nsys;
            float shift = 2.0*M_PI/3.;
            for(int isys=0; isys<nsys; isys++){
            float r = cos(dt*isys           )*0.5 + 0.5;
            float g = cos(dt*isys + shift   )*0.5 + 0.5;
            float b = cos(dt*isys + shift*2 )*0.5 + 0.5;
            opengl1renderer.color3f( r, g, b );
            //printf( "#=========== isys= %i \n", isys );
            Draw3D::neighs_multi(natoms,4,M_neighs,M_neighCell,M_apos, W->pbc_shifts, isys, nvec ); 
        } } 
    }

    //if(iangPicked>=0){
    //    opengl1renderer.color3f(0.,1.,0.);      Draw3D::angle( W->ff.ang2atom[iangPicked], W->ff.ang_cs0[iangPicked], W->ff.apos, fontTex3D );
    //}
    if(useGizmo){ gizmo.draw(); }
    if(bHexDrawing)drawingHex(5.0);
    if(bViewAxis){ opengl1renderer.lineWidth(3);  Draw3D::drawAxis(1.0); opengl1renderer.lineWidth(1); }
};

void MolGUI::printMSystem( int isys, int perAtom, int na, int nvec, bool bNg, bool bNgC, bool bPos ){
    int i0a = (isys * na  )*perAtom;
    int i0v = (isys * nvec); //*perAtom;
    printf( "### MolGUI::MSystem[%i] \n", isys );
    //if(blvec ){printf("lvec\n");  printMat(lvecs [isys]); }
    //if(bilvec){printf("ilvec\n"); printMat(ilvecs[isys]); }
    for(int i=0; i<na; i++){
        int ii=i*perAtom;
        printf("[%i]", i );
        if(bNg )printf("ng (%3i,%3i,%3i,%3i)",  M_neighs   [i0a+ii+0], M_neighs   [i0a+ii+1], M_neighs   [i0a+ii+2], M_neighs   [i0a+ii+3] );
        if(bNgC)printf("ngC(%2i,%2i,%2i,%2i)",  M_neighCell[i0a+ii+0], M_neighCell[i0a+ii+1], M_neighCell[i0a+ii+2], M_neighCell[i0a+ii+3] );
        if(bPos)printf("pa(%6.3f,%6.3f,%6.3f)", M_apos     [i0v+i].x , M_apos     [i0v+i].y , M_apos     [i0v+i].z );
        printf("\n" );
    }
}

void MolGUI::showBonds(  ){
    Vec3d* ps = W->ffl.apos;
    opengl1renderer.begin(GL_LINES);
    for(int i=0; i<bondsToShow.size(); i++ ){
        Vec2i b  = bondsToShow       [i];
        Vec3d d  = bondsToShow_shifts[i];
        Vec3d pi = ps[b.i];
        Vec3d pj = ps[b.j];
        printf( "b[%i,%i] pi(%7.3f,%7.3f,%7.3f) pj(%7.3f,%7.3f,%7.3f) d(%7.3f,%7.3f,%7.3f) \n", pi.x,pi.y,pi.z,    pj.x,pj.y,pj.z,  d.x,d.y,d.z );
        Draw3D::vertex(pi-d);  Draw3D::vertex(pj  );
        Draw3D::vertex(pi  );  Draw3D::vertex(pj+d);
    }
    opengl1renderer.end();
}

void MolGUI::showAtomGrid( char* s, int ia, bool bDraw ){
    if( (ia>=W->ffl.natoms) ){ printf( "ERROR showAtomGrid(%i) out of atom range [0 .. %i] \n", ia, W->ffl.natoms ); }
    const GridShape& grid = W->gridFF.grid;

    Vec3d  pi  = W->ffl.apos[ia];
    Quat4f PLQ = W->ffl.PLQs[ia];

    Quat4f fp = W->gridFF.getForce( pi, Quat4f{ PLQ.x , 0.    , 0.    , 0. }, true );
    Quat4f fl = W->gridFF.getForce( pi, Quat4f{ 0.0   , PLQ.y , 0.    , 0. }, true );
    Quat4f fq = W->gridFF.getForce( pi, Quat4f{ 0.0   , 0.    , PLQ.z , 0. }, true );
    Quat4f fe = fp+fl+fq;

    s += sprintf(s, "GridFF(ia=%i) Etot %15.10f Epaul %15.10f; EvdW %15.10f Eel %15.10f p(%7.3f,%7.3f,%7.3f) \n", ia,   fe.e, fp.e, fl.e, fq.e,  pi.x,pi.y,pi.z  );
    //if( fabs(fe.e-fe_ref.e)>1e-8 ){ s += sprintf(s, "ERROR: getLJQH(%15.10f) Differs !!! \n", fe_ref.e ); }
    if( fe.e>0 ){ opengl1renderer.color3f(0.7f,0.f,0.f); }else{ opengl1renderer.color3f(0.f,0.f,1.f); }
    Draw::drawText( tmpstr, {0, 0}, fontSizeDef, {150,20} );
    opengl1renderer.translatef( 0.0,fontSizeDef*2,0.0 );
}

Vec3d MolGUI::showNonBond( char* s, Vec2i b, bool bDraw ){
    // function to evaluate non-bonded interaction between two atoms and plotting results to screen
    int na = W->ffl.natoms ;
    if( (b.i>=na)||(b.j>=na) ){ printf( "ERROR showNonBond(%i,%i) out of atom range [0 .. %i] \n", b.i,b.j, W->ffl.natoms ); }
    // --- mix non-bonded parameters
    Quat4d* REQs = W->ffl.REQs;
    Vec3d * ps   = W->ffl.apos;
    Mat3d& lvec  = W->ffl.lvec;
    Quat4d REQH  = _mixREQ( REQs[b.i], REQs[b.j] );
    Vec3d pi     = ps[b.i];
    Vec3d pj     = ps[b.j];
    Vec3d d      = pj-pi;
    
    Vec3i g      = W->ffl.invLvec.nearestCell( d );
    Vec3d shift  = lvec.a*g.x + lvec.b*g.y + lvec.c*g.z;
    d.add( shift );

    // --- Reference
    double R2damp = W->gridFF.Rdamp;
    Vec3d f;
    double Eref  = getLJQH( d, f, REQH, R2damp );
    // --- Decomposed
    double  r2  = d.norm2();
    // ---- Electrostatic
    double ir2_ = 1/( r2 + R2damp  );
    double Eel =  COULOMB_CONST*REQH.z*sqrt( ir2_ );
    double Fel =  Eel*ir2_ ;
    // --- LJ 
    double  ir2   = 1/r2;
    double  u2    = REQH.x*REQH.x*ir2;
    double  u6    = u2*u2*u2;
    double vdW    = u6*REQH.y;
    double EvdW   = -2.  *vdW;
    double Epaul  = u6   *vdW;
    double EH     = u6*u6* ((REQH.w<0) ? REQH.w : 0.0);  // H-bond correction
    double Etot   = Epaul + EvdW + EH + Eel;
    //F      +=  (12.*(u6-1.)*vdW + H*6.)*ir2;
    // --- print Energy components to string
    s += sprintf(s, "PLQH[%i,%i] r=%6.3f Etot %15.10f Epaul %15.10f; EvdW %15.10f EH %15.10f Eel %15.10f\n", b.i,b.j, sqrt(r2), Etot, Epaul, EvdW, EH, Eel );
    if( fabs(Etot-Eref)>1e-8 ){ s += sprintf(s, "ERROR: getLJQH(%15.10f) Differs !!! \n", Eref ); }
    if( Etot>0 ){ opengl1renderer.color3f(0.7f,0.f,0.f); }else{ opengl1renderer.color3f(0.f,0.f,1.f); }
    // --- draw to screen
    Draw::drawText( tmpstr, {0, 0}, fontSizeDef, {150,20} );
    opengl1renderer.translatef( 0.0,fontSizeDef*2,0.0 );
    return shift;
}

void MolGUI::drawHUD(){
    opengl1renderer.disable ( GL_LIGHTING );
    gui.draw();

    opengl1renderer.pushMatrix();
    Vec3f textPos = {0, 0, 0};
    if(W->bCheckInvariants){
        opengl1renderer.translatef( 10.0,HEIGHT-20.0,0.0 );
        textPos += {10.0,HEIGHT-20.0,0.0};
        opengl1renderer.color3f(0.5,0.0,0.3);
        //char* s=tmpstr;
        //printf( "(%i|%i,%i,%i) cog(%g,%g,%g) vcog(%g,%g,%g) fcog(%g,%g,%g) torq (%g,%g,%g)\n", ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0,  cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, fcog.x,fcog.y,fcog.z, tq.x,tq.y,tq.z );
        //printf( "neval Ang %i nevalPiSigma %i PiPiT %i PiPiI %i v_av %g \n", ff.nevalAngles, ff.nevalPiSigma, ff.nevalPiPiT, ff.nevalPiPiI, v_av );
        // s += sprintf(s, "iSystemCur %i\n",  W->iSystemCur );
        // if(W->bMMFF) s += sprintf(s, "eval:Ang,ps,ppT,ppI(%i|%i,%i,%i)\n",  W->ff.nevalAngles>0, W->ff.nevalPiSigma>0, W->ff.nevalPiPiT>0, W->ff.nevalPiPiI>0 );
        // s += sprintf(s, "cog (%g,%g,%g)\n", W->cog .x,W->cog .y,W->cog .z);
        // s += sprintf(s, "vcog(%15.5e,%15.5e,%15.5e)\n", W->vcog.x,W->vcog.y,W->vcog.z);
        // s += sprintf(s, "fcog(%15.5e,%15.5e,%15.5e)\n", W->fcog.x,W->fcog.y,W->fcog.z);
        // s += sprintf(s, "torq(%15.5e,%15.5e,%15.5e)\n", W->tqcog.x,W->tqcog.y,W->tqcog.z);
        W->getStatusString( tmpstr, ntmpstr );
        Draw::drawText( tmpstr, textPos, fontSizeDef, {100,20} );
    }
    if(bWriteOptimizerState){
        double T = W->evalEkTemp();   //printf( "T_kinetic=%g[K] \n", T );
        opengl1renderer.translatef( 0.0,fontSizeDef*-5*2,0.0 );
        textPos += {0.0,fontSizeDef*-5*2,0.0};
        opengl1renderer.color3f(0.0,0.5,0.0);
        char* s=tmpstr;
        //dt 0.132482 damp 3.12175e-17 n+ 164 | cfv 0.501563 |f| 3.58409e-10 |v| 6.23391e-09
        double v=sqrt(W->opt.vv);
        double f=sqrt(W->opt.ff);
        //s += sprintf(s,"time/iter=%6.2f[us] bGopt=%i bExploring=%i go.istep=%i T= %7.5f[K] dt=%7.5f damp=%7.5f n+ %4i | cfv=%7.5f |f|=%12.5e |v|=%12.5e \n", time_per_iter, W->bGopt, W->go.bExploring, W->go.istep, T, W->opt.dt, W->opt.damping, W->opt.lastNeg, W->opt.vf/(v*f), f, v );
        s += sprintf(s,"time/iter=%6.2f[us] T= %7.5f[K] dt=%7.5f damp=%7.5f n+ %4i | cfv=%7.5f |f|=%12.5e |v|=%12.5e \n", W->time_per_iter, T, W->opt.dt, W->opt.damping, W->opt.lastNeg, W->opt.vf/(v*f), f, v );
        Draw::drawText( tmpstr, textPos, fontSizeDef, {200,20} );
        opengl1renderer.translatef( 0.0,fontSizeDef*-5*2,0.0 );
        textPos += {0, fontSizeDef*-5*2, 0};
        Draw::drawText( W->info_str(tmpstr), textPos, fontSizeDef, {100,20} );
    }
    opengl1renderer.popMatrix();
    textPos = {0, 0, 0};

    if( W->getMolWorldVersion() == (int)MolWorldVersion::GPU ){
        opengl1renderer.pushMatrix();
        bool  bExplors[W->nSystems];
        float Fconvs  [W->nSystems];
        W->getMultiConf( Fconvs , bExplors );
        opengl1renderer.translatef( 0.0,fontSizeDef*2*30,0.0 );
        textPos += {0.0,fontSizeDef*2*30,0.0};
        opengl1renderer.color3f(0.5,0.,1.);
        for( int i=0; i<W->nSystems; i++ ){  
            sprintf( tmpstr, "SYS[%3i][%i]|F|=%g \n", i, bExplors[i], Fconvs[i] );
            Draw::drawText( tmpstr, textPos, fontSizeDef, {200,20} ); 
            opengl1renderer.translatef( 0.0,fontSizeDef*2,0.0 );
            textPos += {0.0,fontSizeDef*2,0.0};
        };
        opengl1renderer.popMatrix();
    }

    /*
    opengl1renderer.translatef( 0.0,fontSizeDef*-2*2,0.0 );
    if( !bondsToShow_shifts ){
        bondsToShow.push_back( {18,36} );
        bondsToShow.push_back( {35,7}  );
        _alloc( bondsToShow_shifts, bondsToShow.size() );
        // ---- Setup H-bond corrections
        Quat4d* REQs = W->ffl.REQs;
        REQs[18].w = -0.7 * REQs[18].y;
        REQs[36].w = +0.7 * REQs[36].y;
        REQs[35].w = +0.7 * REQs[35].y;
        REQs[7 ].w = -0.7 * REQs[7 ].y;
    }
    for(int i=0; i<bondsToShow.size(); i++ ){
        bondsToShow_shifts[i] =  MolGUI::showNonBond( str, bondsToShow[i] );
    }
    */

    mouse_pix = ((Vec2f){ 2*mouseX/float(HEIGHT) - ASPECT_RATIO,
                          2*mouseY/float(HEIGHT) - 1      });// *(1/zoom);

    if(bConsole) console.draw();
}

static GLMesh<GLMESH_FLAG_COLOR> drawingHexMesh;
void MolGUI::drawingHex(double z0){
    Vec2i ip; Vec2d dp;
    Vec3d p3 = rayPlane_hit( (Vec3d)ray0, (Vec3d)cam.rotMat().c, {0.0,0.0,1.0}, {0.0,0.0,z0} );
    Vec2d p{p3.x,p3.y};
    double off=1000.0;
    bool s = ruler.simplexIndex( p+(Vec2d){off,off}, ip, dp );
    //ruler.nodePoint( ip, p );    opengl1renderer.color3f(1.,1.,1.); Draw3D::drawPointCross( {p.x,p.y, 5.0}, 0.5 );
    Vec3f col = s ? (Vec3f){1, 0.2, 1} : (Vec3f){0.2,1,1};
    ruler.tilePoint( ip, s, p ); Draw3D::drawPointCross( {p.x-off,p.y-off, z0}, 0.2, col );
    
    bool bLine=true;
    if(bDrawHexGrid){
        drawingHexMesh.drawMode = bLine ? GL_LINES : GL_POINTS;
        drawingHexMesh.clear();
        ruler.simplexIndex( (Vec2d){off,off}, ip, dp );
        double sc = ruler.step/sqrt(3.0);
        for(int ix=0;ix<10;ix++ ){
            for(int iy=0;iy<10;iy++ ){
                Vec2i ip_{ip.x+ix,ip.y+iy};
                ruler.tilePoint( ip_, true,  p ); 
                p.sub(off,off);
                if(bLine){
                    col = {1.0,0.2,1.0};
                    Vec2d p2;
                    drawingHexMesh.addVertex(Vec3f{p.x,p.y,z0}, Vec3fZero, col); p2=p+ruler.lvecs[0]*sc; drawingHexMesh.addVertex(Vec3f{p2.x,p2.y,z0}, Vec3fZero, col);
                    drawingHexMesh.addVertex(Vec3f{p.x,p.y,z0}, Vec3fZero, col); p2=p+ruler.lvecs[1]*sc; drawingHexMesh.addVertex(Vec3f{p2.x,p2.y,z0}, Vec3fZero, col);
                    drawingHexMesh.addVertex(Vec3f{p.x,p.y,z0}, Vec3fZero, col); p2=p+ruler.lvecs[2]*sc; drawingHexMesh.addVertex(Vec3f{p2.x,p2.y,z0}, Vec3fZero, col);
                }else{
                    drawingHexMesh.addVertex(Vec3f{p.x,p.y,z0}, Vec3fZero, {1.0, 0.2, 1.0}); ruler.tilePoint( ip_, false, p );  p.add(off,off);
                    drawingHexMesh.addVertex(Vec3f{p.x,p.y,z0}, Vec3fZero, {0.2, 1.0, 1.0});
                }
            }
        }
        drawingHexMesh.draw();
    }
}

void MolGUI::selectRect( const Vec3d& p0, const Vec3d& p1 ){ 
    if(bViewBuilder){
        W->builder.selectRect( p0, p1, (Mat3d)cam.rotMat() ); 
        W->selection.clear();
        for( int i : W->builder.selection){ W->selection.push_back(i); }
        { // Debug
           printf("selection: "); for( int i : W->selection){ printf("%i ", i); } printf("\n");
        }
    }else{ W->selectRect( p0, p1, (Mat3d)cam.rotMat() ); }
    
}

void  MolGUI::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
    int ib = pickBondCenter( W->ff.nbonds, W->ff.bond2atom, W->ff.apos, ro, rd, 0.5 );
    printf( "picked bond  %i \n", ib );
    W->selection.reserve(W->ff.natoms);
    W->splitAtBond( ib, &(W->selection[0]) );
}

/*void MolGUI::renderGridFF( double isoVal, int isoSurfRenderType, double colorSclae ){
    if(verbosity>0) printf( "MolGUI::renderGridFF()\n" );
    //int iatom = 11;
    testREQ = Quat4d{ 1.487, sqrt(0.0006808), 0., 0.}; // H
    testPLQ = REQ2PLQ( testREQ, W->gridFF.alphaMorse );
    Quat4f * FFtot = new Quat4f[ W->gridFF.grid.getNtot() ];
    W->gridFF.evalCombindGridFF ( testREQ, FFtot );
    //W->gridFF.grid.saveXSF( "E_renderGridFF.xsf",  (float*)FFtot, 4, 3, W->gridFF.natoms, W->gridFF.atypes, W->gridFF.apos );
    //if(idebug>0) W->gridFF.grid.saveXSF( "FFtot_z.xsf",  (float*)FFtot, 4, 2, W->gridFF.natoms, W->gridFF.atypes, W->gridFF.apos );
    ogl_isosurf = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_isosurf, GL_COMPILE);
    opengl1renderer.shadeModel( GL_SMOOTH );
    opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.enable(GL_DEPTH_TEST);
    //bool sign=false;
    bool sign=true;
    int nvert = renderSubstrate_( W->gridFF.grid, FFtot, W->gridFF.FFelec, +isoVal, sign, colorSclae );   //printf("Debug: renderGridFF() renderSubstrate() -> nvert= %i ", nvert );
    // ---- This seems still not work properly
    //int ntris=0;
    //opengl1renderer.color3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( W->gridFF.grid,  isoVal, (double*)FFtot, isoSurfRenderType,  3,2 );
    //opengl1renderer.color3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( W->gridFF.grid, -isoVal, (double*)FFtot, isoSurfRenderType,  3,2 );
    //opengl1renderer.color3f(0.,0.,1.); Draw3D::drawTriclinicBox( W->gridFF.grid.cell, Vec3d{0.0, 0.0, 0.0}, Vec3d{1.0, 1.0, 1.0} );
    //Draw3D::drawAxis(1.0);
    opengl1renderer.endList();
    delete [] FFtot;
    if(verbosity>0) printf( "... MolGUI::renderGridFF() DONE\n" );
}*/

//void MolGUI::renderGridFF_new( double isoVal, int isoSurfRenderType, double colorScale, Quat4d REQ = Quat4d{ 1.487, sqrt(0.0006808), 0., 0.} ){
void MolGUI::renderGridFF_new( double isoVal, int isoSurfRenderType, double colorScale, Quat4d REQ ){
    Quat4d PLQ = REQ2PLQ_d( REQ, W->gridFF.alphaMorse );
    printf( "MolGUI::renderGridFF_new() isoVal=%g REQ{%g,%g,%g,%g} PLQ{%g,%g,%g,%g}\n", isoVal, REQ.x, REQ.y, REQ.z,  REQ.z, PLQ.x, PLQ.y, PLQ.z, PLQ.w );

    opengl1renderer.shadeModel( GL_SMOOTH );
    opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.enable(GL_DEPTH_TEST);

    Vec2d zrange{-5.0,5.0};
    renderSubstrate_new( &ogl_isosurf, W->gridFF, Vec2d{zrange.x,zrange.y}, isoVal, PLQ, colorScale );
    if(verbosity>0) printf( "... MolGUI::renderGridFF_new() DONE\n" );
}

void MolGUI::renderSurfAtoms( Vec3i nPBC, bool bPointCross, float qsc, float Rsc, float Rsub ){
    //if(verbosity>0) printf( "MolGUI::renderSurfAtoms() nPBC(%i,%i,%i) qsc=%g Rsc=%g Rsub=%g W->gridFF.apos_.size()=%li bPointCross=%i\n", nPBC.x, nPBC.y, nPBC.z, qsc, Rsc, Rsub, W->gridFF.apos_.size(), bPointCross );
    
    opengl1renderer.shadeModel( GL_SMOOTH );
    opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.enable(GL_DEPTH_TEST);

    Mat3d& lvec =  W->gridFF.grid.cell;
    int icell=0;
    for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
            for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                Vec3d shift = lvec.a*ix + lvec.b*iy + lvec.c*iz;
                Draw3D::atomsREQ( W->gridFF.natoms, W->gridFF.apos, W->gridFF.REQs, qsc, Rsc, Rsub, bPointCross, W->gridFF.shift0 + shift );
                icell++;
            }
        } 
    }
    //if(verbosity>0) printf( "... MolGUI::renderSurfAtoms() DONE\n" );
}

void MolGUI::renderESP( Quat4d REQ ){
    printf( "MolGUI::renderESP() %li \n", ogl_esp );
    ogl_esp = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_esp, GL_COMPILE);
    opengl1renderer.shadeModel( GL_SMOOTH );
    opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.enable(GL_DEPTH_TEST);
    const NBFF& nbmol = W->nbmol;
    int nvert = Draw3D::drawESP( nbmol.natoms, nbmol.apos, nbmol.REQs, REQ );
    opengl1renderer.endList();
};

void MolGUI::renderOrbital(int iMO, double iso ){
    printf( "MolGUI::renderOrbital() \n" );
    double * ewfaux=0;
    W->projectOrbital( iMO, ewfaux );
    if(ewfaux==0)return;
    if(ogl_MO){ opengl1renderer.deleteLists(ogl_MO,1); }
    ogl_MO  = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_MO, GL_COMPILE);
    int ntris=0;  
    opengl1renderer.color3f(0.0,0.0,1.0); ntris += Draw3D::MarchingCubesCross( W->MOgrid,  iso, ewfaux, isoSurfRenderType);
    opengl1renderer.color3f(1.0,0.0,0.0); ntris += Draw3D::MarchingCubesCross( W->MOgrid, -iso, ewfaux, isoSurfRenderType);
    opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::drawTriclinicBox(W->MOgrid.cell.transposed(), Vec3dZero, Vec3dOne );
    opengl1renderer.endList();
    W->bConverged=true;
    bRunRelax=false;
    delete [] ewfaux;
}

void MolGUI::renderDensity(double iso){
    printf( "MolGUI::renderDensity() \n" );
    double * ewfaux=0;
    W->projectDensity( ewfaux );
    if(ewfaux==0)return;
    if(ogl_MO){ opengl1renderer.deleteLists(ogl_MO,1); }
    ogl_MO  = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_MO, GL_COMPILE);
    int ntris = Draw3D::MarchingCubesCross( W->MOgrid, iso, ewfaux, isoSurfRenderType  );
    //printf( "renderOrbital() ntris %i \n", ntris );
    opengl1renderer.endList();
    W->bConverged=true;
    bRunRelax=false;
    delete [] ewfaux;
}


void MolGUI::Fz2df( int nxy, int izmin, int izmax, const Quat4f* afm_Fout, float* dfout ){
    // conversion of vertical force Fz to frequency shift
    // according to:
    // Giessibl, F. J. A direct method to calculate tip-sample forces from frequency shifts in frequency-modulation atomic force microscopy Appl. Phys. Lett. 78, 123 (2001)
    // oscialltion amplitude of cantilever is A = n * dz
    // from Python:
    // x  = np.linspace(-1,1,n+1)
    // y  = np.sqrt(1-x*x)
    // dy =  ( y[1:] - y[:-1] )/(dz*n)
    // fpi    = (n-2)**2
    // prefactor = -1 * ( 1 + fpi*(2/np.pi) ) / (fpi+1) # correction for small n
    // return dy*prefactor,

    // ---- Build the convolution mask
    const int nz = izmax-izmin;
    float dz = (float)afm_scan_grid.dCell.c.z;
    float ws[nz];
    float x  = -1.0; 
    float oy = 0;
    float fpi = sq(nz-2);
    float prefactor = -1 * ( 1 + fpi*(2/M_PI) ) / (fpi+1);
    float  W = prefactor/(nz*dz);
    for(int i=0; i<nz; i++){
        x+=dz;
        float y  = sqrt(1-x*x);
        float dy = y-oy;
        ws[i]    = W*dy;
        oy=y;
    }
    // ---- Fz->df using the convolution mask
    for(int ixy=0; ixy<nxy; ixy++){
        float df = 0.0;
        for(int iz=0; iz<nz; iz++){
            df += afm_Fout[ (iz+izmin)*nxy+ixy ].z * ws[iz];
        }
        dfout[ixy] = df;
    }
}

void MolGUI::renderAFM( int iz, int offset ){
    printf( "MolGUI::renderAFM( iz=%i, offset=%i afm_scan_grid.n.z=%i afm_nconv=%i)\n", iz, offset, afm_scan_grid.n.z, afm_nconv );
    if(afm_Fout==0){ printf("WARRNING: MolGUI::renderAFM() but afm_Fout not allocated \n"); return; };
    if(iz<0 )iz=0;
    if(iz>=afm_scan_grid.n.z)iz=afm_scan_grid.n.z-1;
    int pitch = 4;
    int nxy =afm_scan_grid.n.x * afm_scan_grid.n.y;
    float* data_iz =  (float*)(afm_Fout + iz*nxy);
    float* dfdata=0;
    if(afm_bDf){
        dfdata= new float[nxy];
        float* data_iz = &dfdata[0];
        Fz2df( nxy, iz, iz+afm_nconv, afm_Fout, data_iz );
    }
    printf( "MolGUI::renderAFM() ogl_afm=%i \n", ogl_afm ); //exit(0);
    if(ogl_afm>0)opengl1renderer.deleteLists(ogl_afm,1);
    ogl_afm = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_afm, GL_COMPILE);
    opengl1renderer.shadeModel( GL_SMOOTH );
    //opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.disable(GL_LIGHTING);
    opengl1renderer.enable(GL_DEPTH_TEST);
    double vmin=+1e+300;
    double vmax=-1e+300;
    for(int i=0; i<nxy; i++){ float f=data_iz[i*pitch+offset]; vmin=fmin(f,vmin); vmax=fmax(f,vmax);  }
    //_realloc(afm_ps, nxy);
    //for(int i=0; i<nxy; i++){ afm_ps[i]=(Vec3d)afm_ps0[i].f; }

    opengl1renderer.pushMatrix();
    opengl1renderer.translatef( 0.0,0.0,-5.0 - iz * afm_scan_grid.dCell.c.z );
    printf( "MolGUI::renderAFM() vmin=%g vmax=%g @data_iz=%li @colors_afmhot=%li\n", vmin, vmax, (long)data_iz, (long)Draw::colors_afmhot );
    Draw3D::drawScalarField( {afm_scan_grid.n.x,afm_scan_grid.n.y}, afm_ps0,                                  data_iz, pitch,offset, vmin,vmax, Draw::colors_afmhot );

    //Draw3D::drawColorScale( 10, Vec3dZero, Vec3dY, Vec3dX, Draw::colors_afmhot );
    //Draw3D::drawScalarGrid ( {afm_scan_grid.n.x,afm_scan_grid.n.y}, afm_scan_grid.pos0, afm_scan_grid.dCell.a, afm_scan_grid.dCell.a, data_iz, pitch,offset, vmin, vmax );
    opengl1renderer.popMatrix();
    opengl1renderer.endList();
    if(dfdata){ delete [] dfdata; }
};

void MolGUI::renderAFM_trjs( int di ){
    if(afm_Fout==0){ printf("WARRNING: MolGUI::renderAFM_trjs() but afm_PPpos not allocated \n"); return; };
    printf( "MolGUI::renderAFM_trjs(di=%i) afm_PPpos=%li \n", di, afm_PPpos ); //exit(0);
    if(ogl_afm_trj>0)opengl1renderer.deleteLists(ogl_afm_trj,1);
    ogl_afm_trj = opengl1renderer.genLists(1);
    opengl1renderer.newList(ogl_afm_trj, GL_COMPILE);
    opengl1renderer.shadeModel( GL_SMOOTH );
    //opengl1renderer.enable(GL_LIGHTING);
    opengl1renderer.disable(GL_LIGHTING);
    opengl1renderer.enable(GL_DEPTH_TEST);
    Vec3i ns = afm_scan_grid.n;
    int nxy=ns.x*ns.y;
    for(int iy=0; iy<ns.y; iy+=di ){
        for(int ix=0; ix<ns.x; ix+=di ){
            //afm_ps0
            opengl1renderer.begin(GL_LINE_STRIP);
            for(int iz=0; iz<ns.z; iz++ ){
                int i = afm_scan_grid.n.y;
                Vec3f p = afm_PPpos[ iz*nxy + iy*ns.x + ix ].f;
                opengl1renderer.vertex3f( p.x,p.y,p.z ); 
            }
            opengl1renderer.end();
        }
    }
    opengl1renderer.endList();
};


void MolGUI::makeAFM( int iz ){
    printf( "\n ==== MolGUI::makeAFM() %li \n", ogl_afm ); //exit(0);
    afm_ff_grid.cell = W->builder.lvec;
    afm_ff_grid.updateCell(0.1);
    W->evalAFM_FF ( afm_ff_grid,   afm_ff,                        false );
    W->evalAFMscan( afm_scan_grid, afm_Fout, afm_PPpos, &afm_ps0, false );
    MolGUI::renderAFM_trjs( 5 );
    //MolGUI::renderAFM(  afm_scan_grid.n.z-10, 2 );
    if(iz>0){ afm_iz=iz; }
    MolGUI::renderAFM(  afm_iz, 2 );
    printf( "\n ==== MolGUI::makeAFM() DONE \n\n", ogl_afm ); //exit(0);
};

void MolGUI::drawPi0s( float sc=1.0 ){
    const MMFFsp3& ff = W->ff;
    opengl1renderer.begin(GL_LINES);
    for(int ia=0; ia<ff.nnode; ia++){
        int* ngs = ff.neighs + ia*ff.nneigh_max;
        for(int j=0; j<ff.nneigh_max; j++){
            int ing = ngs[j]; 
            if(ing<0){
                int ipi = -ing-2;
                Vec3f p=(Vec3f)ff.apos[ia];      opengl1renderer.vertex3f(p.x,p.y,p.z);
                p.add_mul( ff.pi0s[ipi].f, sc);  opengl1renderer.vertex3f(p.x,p.y,p.z);
                //printf("drawPi0s[%i,%i|%i] (%g,%g,%g)\n", ia, j, ipi, ff.pi0s[ipi].f.z, ff.pi0s[ipi].f.y, ff.pi0s[ipi].f.z );
            }
        }
    }
    opengl1renderer.end();
}

/*
void MolGUI::flipPis( Vec3d ax ){
    for(int i=0; i<nnode; i++){
        double c = pipos[i].dot(ax);
        if( c<0 ){ pipos[i].mul(-1); } 
    }
}
*/

void MolGUI::bindMolecule( int natoms_, int nnode_, int nbonds_, int* atypes_,Vec3d* apos_,Vec3d* fapos_,Quat4d* REQs_, Vec3d* pipos_, Vec3d* fpipos_, Vec2i* bond2atom_, Vec3d* pbcShifts_ ){
    natoms=natoms_; nnode=nnode_; nbonds=nbonds_;
    if(atypes_)atypes=atypes_;
    if(apos_  )apos  =apos_;
    if(fapos_ )fapos =fapos_;
    if(REQs_  )REQs  =REQs_;
    if(pipos_ )pipos =pipos_;
    if(fpipos_)fpipos=fpipos_;
    if(bond2atom_)bond2atom=bond2atom_;
    if(pbcShifts_)pbcShifts=pbcShifts_;
}

void MolGUI::bindMolecule(const MolWorld_sp3* W ){
    printf( "MolGUI::bindMolecule() \n" );
    //natoms=ff.natoms; nnode=ff.nnode; nbonds=ff.nbonds;
    //atypes=ff.atypes; apos=ff.apos; fapos=ff.fapos; REQs=ff.REQs; pipos=ff.pipos; fpipos=ff.fpipos; bond2atom=ff.bond2atom; pbcShifts=ff.pbcShifts;
    if(W->bUFF){
        bindMolecule( W->ffu.natoms, 0, W->ffu.nbonds, W->ffu.atypes, W->ffu.apos, W->ffu.fapos, W->ffu.REQs, 0, 0, W->ffu.bonAtoms, W->ffu.shifts  );
        neighs    = W->ffu.neighs;
        neighCell = W->ffu.neighCell;
    }else{    
        bindMolecule( W->ffl.natoms, W->ffl.nnode, W->ff.nbonds, W->nbmol.atypes, W->nbmol.apos, W->nbmol.fapos, W->nbmol.REQs, W->ffl.pipos, W->ffl.fpipos, W->ff.bond2atom, W->ff.pbcShifts );
        neighs    = W->ffl.neighs;
        neighCell = W->ffl.neighCell;
    }
    constrs   = (Constrains*)&W->constrs;
}

void MolGUI::unBindMolecule(){
    natoms=0; nnode=0; nbonds=0; atypes=0; apos=0; fapos=0; REQs=0; pipos=0; fpipos=0; bond2atom=0; pbcShifts=0;
    neighs=0; neighCell=0;
}

void MolGUI::makeBondLengths0(){
    int nb = W->builder.bonds.size();
    _realloc0( bL0s, nb, -1. );
    for( int i=0; i<nb; i++ ){
        const MM::Bond& b = W->builder.bonds[i];
        //printf( "MolGUI::makeBondLengths0() [%i] l0=%g \n", i, b.l0 );
        bL0s[i] = b.l0; 
    }
}

void MolGUI::makeBondColoring( Vec2i typs, Vec2d lrange, double*& clr,bool bNew){
    //printf( "MolGUI::makeBondColoring() \n" ); //exit(0);
    int nb = W->builder.bonds.size();
    if(bNew){ _realloc0( clr, nb, -1. ); }
    double invr = 1./(lrange.y-lrange.x);
    for( int i=0; i<nb; i++ ){
        //const MM::Bond& b = W->builder.bonds[i];
        Vec2i b = W->builder.bonds[i].atoms;
        int ti  = W->builder.atoms[b.a].type;
        int tj  = W->builder.atoms[b.b].type;
        double l = (W->ffl.apos[b.b]-W->ffl.apos[b.a]).norm();
        if( (ti==typs.x) && (tj==typs.y) ){
            double c = (l-lrange.x)*invr;
            c = _clamp( c, -1.+1.e-6, 1.-1.e-6 );
            clr[i] = c;
        }else{
            clr[i] = -1;
        }
        //printf( "MolGUI::makeBondLengths0() [%i] l0=%g \n", i, b.l0 );
    }
}

void MolGUI::drawBuilder( Vec3i ixyz ){
    //printf( "MolGUI::drawBuilder() ixyz(%i,%i,%i)\n", ixyz.x,ixyz.y,ixyz.z );
    //float textSize=0.015;
    //if(bViewColorFrag){ printf( " B.frags.size(%i) \n", B.frags.size() ); }
    opengl1renderer.enable(GL_DEPTH_TEST);
    MM::Builder& B = W->builder;
    bool bOrig = (ixyz.x==0)&&(ixyz.y==0)&&(ixyz.z==0);
    if(bViewAtomSpheres&&mm_bAtoms ){       
        opengl1renderer.enable(GL_LIGHTING);
        opengl1renderer.shadeModel(GL_SMOOTH);                     
        for(int ia=0; ia<B.atoms.size(); ia++){
            const MM::Atom& a = B.atoms[ia];
            //const MM::AtomType& atyp = W->params.atypes[a.type];
            const AtomType& atyp = W->params.atypes[a.type];
            Vec3f atomColor;
            if(bViewColorFrag){ 
                if(a.frag<0) continue;
                if( (a.frag<0)||(a.frag>=B.frags.size()) ){ printf( "ERROR MolGUI::drawBuilder() a.frag(%i) out of range 0 .. B.frags.size(%i) \n", a.frag, B.frags.size() ); exit(0); }
                Draw::setRGB( B.frags[a.frag].color );
                atomColor = COL2VEC(B.frags[a.frag].color);
            }else{ Draw::setRGB( W->params.atypes[a.type].color ); atomColor = COL2VEC(W->params.atypes[a.type].color); }
            float sz = (atyp.RvdW-mm_Rsub)*mm_Rsc;
            Draw3D::drawSphere((Vec3f)a.pos, sz, atomColor);
        }    
    }
    if(mm_bAtoms&&bViewAtomLabels ){ 
        opengl1renderer.color3f(0.0f,0.0f,0.0f); 
        if(bViewMolCharges){
            //void atomPropertyLabel( int n, double* data, Vec3d* ps, int pitch, int offset, int fontTex, float sz=0.02, const char* format="%4.2f\0" ){
            for(int i=0; i<B.atoms.size(); i++){
                Draw3D::drawDouble( B.atoms[i].pos, B.atoms[i].REQ.z, fontTex, textSize, "%4.2f" );
                //drawInt( ps[i], (int)data[i*pitch+offset], fontTex, sz );
            }
        }else
        for(int ia=0; ia<B.atoms.size(); ia++){ 
            int ii=ia;
            if(bViewColorFrag){ ii = B.atoms[ia].frag; }
            Draw3D::drawInt( B.atoms[ia].pos, ii, fontTex, textSize );
        } 
    }
    if( bViewBonds ){
        opengl1renderer.disable(GL_LIGHTING);
        opengl1renderer.color3f(0.0f,0.0f,0.0f);
        opengl1renderer.lineWidth(3.0f);
        opengl1renderer.begin(GL_LINES);
        for(int ib=0; ib<B.bonds.size(); ib++){
            Vec2i b  = B.bonds[ib].atoms;
            Draw3D::vertex(B.atoms[b.a].pos); 
            Draw3D::vertex(B.atoms[b.b].pos);
        }
        opengl1renderer.end();
        opengl1renderer.lineWidth(1.0);
    }
}

void MolGUI::drawSystemShifts( int n, const Vec3d* shifts, int i0 ){
    opengl1renderer.enable(GL_DEPTH_TEST);

    bool bViewBL = bViewBondLenghts &&  (bL0s!=0);

    // bonds
    if( neighs && (!bViewBL) ){
        opengl1renderer.color3f(0.0f,0.0f,0.0f); 
        opengl1renderer.lineWidth(1.0);
        GLMesh<0>* neighMesh = Draw3D::makeNeighsMesh(  natoms, 4, (int*)neighs, (int*)neighCell, apos, W->pbc_shifts );

        for (int i=0; i<n; i++) neighMesh->draw( (Vec3f)shifts[i] );
    }

    // atom spheres
    if(bViewAtomSpheres && mm_bAtoms && (!bViewBondLenghts)){ // TODO
        for (int i=0; i<n; i++) Draw3D::atoms( natoms, apos, atypes, W->params, 1.0, mm_Rsc, mm_Rsub, shifts[i] );
    }

    // == i0 system gets special rendering ==
    if(bViewAtomForces    &&  fapos           ){ opengl1renderer.color3f(1.0f,0.0f,0.0f); Draw3D::drawVectorArray  ( natoms, apos, fapos, ForceViewScale, 10000.0 );   }
    if(mm_bAtoms&&bViewAtomLabels&&(!bViewBL) ){ opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::atomLabels       ( natoms, apos,                                    fontTex3D, textSize );  }
    if(mm_bAtoms&&bViewAtomTypes              ){ opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::atomTypes        ( natoms, apos, atypes, &(params_glob->atypes[0]), fontTex3D, textSize );  }
    if(bViewMolCharges && (W->nbmol.REQs!=0)  ){ opengl1renderer.color3f(0.0,0.0,0.0);    Draw3D::atomPropertyLabel( natoms,  (double*)REQs,  apos, 4, 2,             fontTex3D, textSize ); }
    if(bViewHBondCharges && (W->nbmol.REQs!=0)){ opengl1renderer.color3f(0.0,0.0,0.0);    Draw3D::atomPropertyLabel( natoms,  (double*)REQs,  apos, 4, 3,             fontTex3D, textSize ); }
    opengl1renderer.enable( GL_DEPTH_TEST );
    
    {// Graph

        // --- draw whole molecule skeleton stored in W->graph
        //opengl1renderer.color3f(1.0,0.0,1.0);
        //for(int i=0; i<W->graph.n; i++){ 
        //    for(int j=0; j<W->graph.nneighs[i]; j++ ) Draw3D::drawLine( apos[i], apos[W->graph.neighs[i][j]] ); 
        //};

        // --- draw only the bridge bonds stored in W->graph.found
        for(int i=0; i<W->graph.found.size(); i++){ Vec2i b = W->graph.found[i]; Draw3D::drawLine( apos[b.i], apos[b.j], {1, 0, 1} );  };
    }
    if( bond2atom ){
        if(bViewPis &&  fpipos ){ opengl1renderer.color3f(0.0f,1.0f,1.0f); Draw3D::drawVectorArray( nnode, apos, pipos, 1.0, 100.0 );          }
        if( bViewBondLabels    ){ opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::bondLabels( nbonds, bond2atom, apos, fontTex3D,  textSize  ); }
        opengl1renderer.enable( GL_DEPTH_TEST ); 
        
        if(bViewBondLenghts){ 
            if(bViewAtomLabels ){ opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw3D::bondsLengths( nbonds, bond2atom, apos, fontTex3D, textSize ); }
            opengl1renderer.enable( GL_DEPTH_TEST );
            //if(bL0s==0){ makeBondLengths0(); }
            opengl1renderer.lineWidth( 10.0 );
            //Draw3D::bondLengthColorMap( nbonds, bond2atom, apos, bL0s, 0.01 );
            //Draw3D::bondLengthColorMap(nbonds, bond2atom, apos, Vec2d{2.33,2.35} );
            //printf( "@bL0s=%li\n", (long)bL0s );
            if(bL0s) Draw3D::bondLengthColorMap(nbonds, bond2atom, apos, bL0s );
            opengl1renderer.lineWidth( 1.0 );
        }
    }
}

void MolGUI::saveScreenshot( int i, const char* fname ){
    //char str[64];
    sprintf( tmpstr, fname, i );               
    printf( "save to %s \n", tmpstr );
    unsigned int *screenPixels = new unsigned int[WIDTH*HEIGHT*4];  
    opengl1renderer.flush();                                                      
    opengl1renderer.finish();                                                     
    //opengl1renderer.readPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_INT, screenPixels);
    opengl1renderer.readPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, screenPixels);  
    //SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff );   
    SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000 );  
    SDL_SaveBMP(bitmap, tmpstr);   
    SDL_FreeSurface(bitmap);
    delete[] screenPixels;
}

void MolGUI::scanSurfFF( int n, Vec3d p0, Vec3d p1, Quat4d REQ, int evalMode, int viewMode, double sc, Vec3d hat ){
    printf( " MolGUI::scanSurfFF() evalMode=%i viewMode=%i n=%i p0(%g,%g,%g) p1(%g,%g,%g) REQ(%g,%g,%g,%g) ogl_surf_scan=%i sc=%g hat(%g,%g,%g)\n", evalMode,viewMode,  n, p0.x,p0.y,p0.z,  p1.x,p1.y,p1.z,    REQ.x,REQ.y,REQ.z,REQ.w,    ogl_surf_scan, sc, hat.x,hat.y,hat.z );
    if( !ogl_surf_scan ){ opengl1renderer.deleteLists(ogl_surf_scan,1); }
    ogl_surf_scan = opengl1renderer.genLists(1);
    Vec3d dp=p1-p0; dp.mul(1./n);
    opengl1renderer.newList(ogl_surf_scan, GL_COMPILE );
    opengl1renderer.begin(GL_LINES);
    Quat4d PLQ = REQ2PLQ_d( REQ, W->gridFF.alphaMorse );
    Vec3d op2;
    for(int i=0; i<n; i++){
        Vec3d  p = p0 + dp*i;
        Quat4d fe;
        if     ( evalMode==0 ){ fe = W->gridFF.getForce_d       ( p, PLQ, true ); }
        else if( evalMode==1 ){ fe = W->gridFF.getForce_Tricubic( p, PLQ, true ); }
        //printf(   "MolGUI::scanSurfFF()[%i] p(%g,%g,%g) fe(%g,%g,%g|%g)\n", i, p.x,p.y,p.z,   fe.x,fe.y,fe.z,fe.w  );
        Vec3d p2;
        if     ( viewMode==0 ){ p2 = p + hat*fe.e*sc; }  // Energy
        else if( viewMode==1 ){ p2 = p +     fe.f*sc; }  // force
        Draw3D::vertex( p   ); Draw3D::vertex(p2);
        if(i>0){ 
            Draw3D::vertex( p-dp ); Draw3D::vertex(p);
            Draw3D::vertex( op2  ); Draw3D::vertex(p2); 
        }
        op2=p2;
    }
    opengl1renderer.end();
    opengl1renderer.endList();
}

void MolGUI::debug_scanSurfFF( int n, Vec3d p0, Vec3d p1, Quat4d REQ, double sc ){
    if(isnan(sc)){ sc=ForceViewScale; }
    //printf("========= MolGUI::scanFF\n");
    //sc=10.0;
    //p0=Vec3d{0.0,0.0,z0_scan}; p1=Vec3d{5.0,5.0,z0_scan}; // Horizontal scan
    //p0=Vec3d{0.0,0.0,z0_scan}; p1=Vec3d{0.0,5.0,z0_scan}; // Horizontal scan
    //p0=Vec3d{0.0,z0_scan,0.0}; p1=Vec3d{0.0,z0_scan,10.0,}; // Vertical scan
    Vec3d dp=p1-p0; dp.mul(1./n);
    //Quat4f PLQ = REQ2PLQ( REQ, W->gridFF.alphaMorse );
    //printf( "PLQ %6.3f %10.7f %6.3f \n", PLQ.x,PLQ.y,PLQ.z   );
    debug_ps  .resize(n);
    debug_fs  .resize(n);
    debug_REQs.resize(n);
    for(int i=0; i<n; i++){ 
        debug_REQs[i]  =(Quat4f)REQ;
        debug_ps  [i].f=(Vec3f )(p0+dp*i);
    }
    W->scanSurfFF( n, &debug_ps[0], &debug_REQs[0], &debug_fs[0] );
    opengl1renderer.begin(GL_LINES);
    for(int i=0; i<n; i++){
        //Vec3d  p  = p0 + dp*i;
        //Quat4f fe = Quat4fZero; W->gridFF.addForce_surf( p, PLQ, fe );
        //Quat4f fe = W->gridFF.getForce( p, PLQ, true );
        //printf( "debug_scanSurfFF[%i] p(%6.3f,%6.3f,%6.3f) fe(%g,%g,%g|%g)\n", i,  p.x,p.y,p.z, fe.x,fe.y,fe.z,fe.w );
        Vec3f f = debug_fs[i].f;
        Vec3f p = debug_ps[i].f;

        Draw3D::vertex( p ); Draw3D::vertex(p + (Vec3f)dp   );
        Draw3D::vertex( p ); Draw3D::vertex(p + f*sc        );

        //double E=0;
        //if   (W->bGridFF){ E+= W->gridFF.eval        (1, W->nbmol.apos, W->nbmol.PLQs,  W->nbmol.fapos );                        }
        //else           { E+= W->nbmol .evalMorse   (W->surf, false,                        W->gridFF.alphaMorse, W->gridFF.Rdamp ); }
        //else             { E+= W->nbmol .evalMorsePBC(W->surf, W->gridFF.grid.cell, W->nPBC, W->gridFF.alphaMorse, W->gridFF.Rdamp ); }
        //Draw3D::vertex( W->nbmol.apos[0] ); Draw3D::vertex( W->nbmol.apos[0]+W->nbmol.fapos[0]*sc );       // Force Vectro
        //Draw3D::vertex( W->nbmol.apos[0] ); Draw3D::vertex( W->nbmol.apos[0]+(Vec3d{0.0,0.0,E})*sc  );  // Energy -> z
        //Draw3D::vertex( W->nbmol.apos[0] ); Draw3D::vertex( W->nbmol.apos[0]+(Vec3d{0.0,E,0.0})*sc  );  // Energy -> x
    }
    opengl1renderer.end();
}

void MolGUI::lattice_scan( int n1, int n2, const Mat3d& dlvec ){
    printf( "MolGUI::lattice_scan(%i,%i, dvec", n1,n2  ); printMat(dlvec);
    long T0 = getCPUticks();
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    W->optimizeLattice_1d( n1,n2, dlvec ); 
    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    double time_s     = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() * 1e-9;
    double time_GTick = (getCPUticks()-T0)*1e-9;
    printf( "Time{W->optimizeLattice_1d(20,20)} %g[s] %g[GTick]  %g[GTick/s] \n", time_s,time_GTick, time_GTick/time_s  );
}

void MolGUI::mouse_default( const SDL_Event& event ){
    //printf( "MolGUI::mouse_default() \n" );
    switch( event.type ){
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                //case SDL_BUTTON_LEFT: { ray0_start = ray0;  bDragging = true; }break;
                case SDL_BUTTON_LEFT: mouseStartSelectionBox(); break;
                case SDL_BUTTON_RIGHT:{ }break;
            } break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    if( ray0.dist2(ray0_start)<0.1 ){ // too small for selection box 
                        int ipick = pickParticle( (Vec3d)ray0, (Vec3d)cam.rotMat().c, 0.5, W->nbmol.natoms, W->nbmol.apos );
                        if( ipick>=0 ){ 
                            printf( "MolGUI::mouse_default() ipick=%i \n", ipick ); 
                            W->selection.clear();  
                        };
                        if( ipick == W->ipicked ){  W->ipicked=-1; }else{ W->ipicked = ipick; };
                        if(W->ipicked>=0){ 
                            //printf( "picked atom %i \n", W->ipicked ); 
                            //printf( "MolGUI::mouse_default() SDL_BUTTON_LEFT W->selection.clear(); (B) ipick=%i W->ipicked=%i\n", ipick, W->ipicked );
                            //W->selection.clear();
                            W->selection.push_back(W->ipicked); 
                            //if( Qpanel==0 ){ printf( "MolGUI::mouse_default() Qpanel==0 \n" ); exit(0); }
                            // printf( "MolGUI::mouse_default() @Qpanel=%li @W=%li W->nbmol.REQs=%li W->ipicked=%i \n", (long)Qpanel, (long)W, (long)W->nbmol.REQs, W->ipicked );
                            // int ip = W->ipicked;             
                            // double z = W->nbmol.REQs[ip].z;  
                            // Qpanel->value = z;               
                            if( Qpanel ){ 
                                Qpanel->value = W->nbmol.REQs[W->ipicked].z;
                            }
                            
                        };
                        //printf( "picked atom %i \n", W->ipicked );
                    }else{
                        selectRect( (Vec3d)ray0_start, (Vec3d)ray0 );
                    }
                    bDragging=false;
                    break;
                case SDL_BUTTON_RIGHT:{ 
                    // int ib = W->builder.pickBond( (Vec3d)ray0, (Vec3d)cam.rotMat().c, 0.3 );
                    // //printf( "MolGUI::pickBond: %i \n", ib  );
                    // if(ib>=0){ printf( "MolGUI::delete bond: %i \n", ib  );  W->builder.deleteBond(ib); bBuilderChanged=true; }
                    W->ipicked=-1; 
                } break;
            }break;
    }
}

void MolGUI::eventMode_edit( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN : { if(gui.bKeyEvents)switch( event.key.keysym.sym ){
                case SDLK_t:{
                    affineTransform( W->ff.natoms, W->ff.apos, W->ff.apos, W->builder.lvec, W->new_lvec );
                    W->builder.updatePBC( W->ff.pbcShifts, &(W->new_lvec) );
                    _swap( W->builder.lvec, W->new_lvec );
                }break;
                case SDLK_i:
                    //selectShorterSegment( (Vec3d)(cam.rotMat().a*mouse_begin_x + cam.rotMat().b*mouse_begin_y + cam.rotMat().c*-1000.0), (Vec3d)cam.rotMat().c );
                    selectShorterSegment( (Vec3d)ray0, (Vec3d)cam.rotMat().c );
                    //selection.erase();
                    //for(int i:builder.selection){ selection.insert(i); };
                    break;
            }}break;
        case SDL_WINDOWEVENT:{switch (event.window.event) {case SDL_WINDOWEVENT_CLOSE:{ quit(); }break;} } break;
    }
}

void MolGUI::eventMode_scan( const SDL_Event& event  ){
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN : { if(gui.bKeyEvents)switch( event.key.keysym.sym ){
            //case SDLK_KP_1: picked_lvec = &W->builder.lvec.a; break;
            //case SDLK_KP_2: picked_lvec = &W->builder.lvec.b; break;
            //case SDLK_KP_3: picked_lvec = &W->builder.lvec.c; break;
            //case SDLK_KP_7: picked_lvec->x+=xstep; break;
            //case SDLK_KP_4: picked_lvec->x-=xstep; break;
            //case SDLK_KP_8: picked_lvec->y+=xstep; break;
            //case SDLK_KP_5: picked_lvec->y-=xstep; break;
            //case SDLK_KP_9: picked_lvec->z+=xstep; break;
            //case SDLK_KP_6: picked_lvec->z-=xstep; break;

            //case SDLK_COMMA:  W->change_lvec( W->ffl.lvec+dlvec    ); break;
            //case SDLK_PERIOD: W->change_lvec( W->ffl.lvec+dlvec*-1 ); break;

            case SDLK_COMMA:  W->add_to_lvec( dlvec    ); break;
            case SDLK_PERIOD: W->add_to_lvec( dlvec*-1 ); break;

            case SDLK_a: apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
            case SDLK_d: apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

            case SDLK_w: apos[1].mul( 1.1 ); break;
            case SDLK_s: apos[1].mul( 0.9 ); break;

            case SDLK_KP_PLUS:  z0_scan = z0_scan+0.1; break;
            case SDLK_KP_MINUS: z0_scan = z0_scan-0.1; break;

            case SDLK_SPACE: 
                bRunRelax=!bRunRelax; 
                // printf( "bRunRelax %i \n", bRunRelax );
                // if(!bRunRelax){ if(ogl_MO>0){ int iHOMO = W->getHOMO(); renderOrbital( iHOMO + which_MO );  } }
                break;

        }  }break;
        case SDL_MOUSEBUTTONDOWN:
        case SDL_MOUSEBUTTONUP:
            mouse_default( event );
            break;
        case SDL_WINDOWEVENT:{switch (event.window.event) {case SDL_WINDOWEVENT_CLOSE:{ quit(); }break;} } break;
    }
}

void MolGUI::eventMode_default( const SDL_Event& event ){
    if(useGizmo)gizmo.onEvent( mouse_pix, event );
    //printf( "MolGUI::eventMode_default() bConsole=%i \n", bConsole );
    switch( event.type ){
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_KEYDOWN : 
                if (bConsole){ bConsole=console.keyDown( event.key.keysym.sym ); }
                else 
                if(gui.bKeyEvents) switch( event.key.keysym.sym ){
                case SDLK_KP_0: cam.setQrot(Quat4fIdentity); break;

                //case SDLK_COMMA:  which_MO--; printf("which_MO %i \n", which_MO ); break;
                //case SDLK_PERIOD: which_MO++; printf("which_MO %i \n", which_MO ); break;

                //case SDLK_INSERT:   break;
                case SDLK_DELETE:   {  bBuilderChanged = W->deleteAtomSelection()>0;      W->clearSelections();     } break;

                //case SDLK_HOME:     break;
                //case SDLK_END:      break;

                //case SDLK_PAGEUP  : W->add_to_lvec( dlvec2     ); break;
                //case SDLK_PAGEDOWN: W->add_to_lvec( dlvec2*-1  ); break;

                case SDLK_BACKQUOTE:{ bConsole = !bConsole;}break;   // ` SDLK_ for key '`' 

                case SDLK_0:      W->add_to_lvec( dlvec2     ); break;
                case SDLK_9:      W->add_to_lvec( dlvec2*-1  ); break;

                case SDLK_COMMA:  W->add_to_lvec( dlvec    ); break;
                case SDLK_PERIOD: W->add_to_lvec( dlvec*-1 ); break;

                case SDLK_SEMICOLON:    afm_iz++; if(afm_iz>=afm_scan_grid.n.z-afm_nconv)afm_iz=0;  renderAFM(afm_iz,2); break;
                case SDLK_QUOTE:        afm_iz--; if(afm_iz<0)afm_iz=afm_scan_grid.n.z-1-afm_nconv; renderAFM(afm_iz,2);  break;

                //case SDLK_LESS:    which_MO--; printf("which_MO %i \n"); break;
                //case SDLK_GREATER: which_MO++; printf("which_MO %i \n"); break;

                //case SDLK_m: renderOrbital( which_MO ); break;
                //case SDLK_r: renderDensity(          ); break;
                //case SDLK_s: W->saveXYZ( "out.xyz", "#comment", false ); break;
                case SDLK_s: { 
                    ((Atoms*)(&W->ffl))->lvec = &(W->ffl.lvec);
                    FILE*file=fopen("snapshot.xyz","w");
                    std::vector<int> iZs(W->ffl.natoms,-1); for(int i=0; i<W->ffl.natoms; i++){ iZs[i]=W->params.atypes[W->ffl.atypes[i]].iZ; } 
                    int* atypes_bak=W->ffl.atypes;
                    W->ffl.atypes=&iZs[0];
                    W->ffl.atomsToXYZ(file, true,true, Vec3i{5,7,1}, "", false );
                    W->ffl.atypes=atypes_bak;
                    fclose(file);
                    } break;

                //case SDLK_p: saveScreenshot( frameCount ); break;
                case SDLK_p:{ 
                    saveScreenshot( frameCount );
                    W->renderSVG( "screenshot.svg",     {0,0,0}, Mat3dIdentity );
                    W->renderSVG( "screenshot_2x2.svg", {1,1,0}, Mat3dIdentity );
                    W->saveXYZ( "screenshot.xyz",     "#comment", false, "w", {1,1,1} );
                    W->saveXYZ( "screenshot_3x3.xyz", "#comment", false, "w", {3,3,1} );
                }break;

                // case SDLK_h:if( ( W->getMolWorldVersion() & MolWorldVersion::QM ) ){ printf( "makeAFM(): is supported only in GPU version of MolWorld \n" ); }else{
                //     //int iMO = which_MO;
                //     int iHOMO = W->getHOMO(); printf( "plot HOMO+%i (HOMO=eig#%i) \n", iHOMO+which_MO, iHOMO );
                //     renderOrbital( iHOMO + which_MO ); break;
                // }break;

                case SDLK_h:{
                    bViewHBondCharges  ^= 1;
                    //printf( "view Hydrogen Bonds \n" );
                }break;
                
                //case SDLK_o: W->optimizeLattice_1d( 10,40, Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  } ); break;
                case SDLK_o:{
                    MolGUI::lattice_scan( 20,20, Mat3d{   0.0,0.5,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0 } );    
                    // long T0 = getCPUticks();
                    // std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                    // W->optimizeLattice_1d( 20,20, Mat3d{   0.0,0.5,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  } ); 
                    // std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                    // double time_s     = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() * 1e-9;
                    // double time_GTick = (getCPUticks()-T0)*1e-9;
                    // printf( "Time{W->optimizeLattice_1d(20,20)} %g[s] %g[GTick]  %g[GTick/s] \n", time_s,time_GTick, time_GTick/time_s  );
                    }break;
                
                case SDLK_u:{ 
                        long T0 = getCPUticks();
                        std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
                        W->upload_pop        ( "population.xyz" ); 
                        std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
                        double time_s     = std::chrono::duration_cast<std::chrono::nanoseconds>(end-begin).count() * 1e-9;
                        double time_GTick = (getCPUticks()-T0)*1e-9;
                        printf( "Time{W->upload_pop(population.xyz)} %g[s] %g[GTick]  %g[GTick/s] \n", time_s,time_GTick, time_GTick/time_s  );
                    }break;
                //case SDLK_u: W->upload_pop        ( "population.xyz" ); break;
                //case SDLK_o: W->optimizeLattice_1d( 0,2, Mat3d{   0.2,0.0,0.0,    0.0,0.0,0.0,    0.0,0.0,0.0  } ); break;

                case SDLK_LEFTBRACKET:  W->prevSystemReplica(); break;
                case SDLK_RIGHTBRACKET: W->nextSystemReplica(); break;

                //case SDLK_LEFTBRACKET:  myAngle-=0.1; printf( "myAngle %g \n", myAngle ); break;
                //case SDLK_RIGHTBRACKET: myAngle+=0.1; printf( "myAngle %g \n", myAngle );  break;

                //case SDLK_g: useGizmo=!useGizmo; break;
                //case SDLK_g: W->bGridFF=!W->bGridFF; break;
                //case SDLK_g: W->swith_gridFF(); break;
                case SDLK_c: W->autoCharges(); break;
                
                case SDLK_v:{ 
                    if( W->getMolWorldVersion() & (int)MolWorldVersion::GPU ){ makeAFM(); }else{ printf( "makeAFM(): is supported only in GPU version of MolWorld \n" ); } 
                    } break;
                case SDLK_KP_MULTIPLY:  afm_iz++; if(afm_iz>=afm_scan_grid.n.z-afm_nconv)afm_iz=0;  renderAFM(afm_iz,2); break;
                case SDLK_KP_DIVIDE:    afm_iz--; if(afm_iz<0)afm_iz=afm_scan_grid.n.z-1-afm_nconv; renderAFM(afm_iz,2);  break;
                
                //case SDLK_SEMICOLON:    afm_iz++; if(afm_iz>=afm_scan_grid.n.z-afm_nconv)afm_iz=0;  renderAFM(afm_iz,2); break;
                //case SDLK_QUOTE:        afm_iz--; if(afm_iz<0)afm_iz=afm_scan_grid.n.z-1-afm_nconv; renderAFM(afm_iz,2);  break;

                //case SDLK_0:            afm_iz++; if(afm_iz>=afm_scan_grid.n.z-afm_nconv)afm_iz=0;  renderAFM(afm_iz,2); break;
                //case SDLK_9:            afm_iz--; if(afm_iz<0)afm_iz=afm_scan_grid.n.z-1-afm_nconv; renderAFM(afm_iz,2);  break;

                case SDLK_LEFTPAREN:    afm_iz++; if(afm_iz>=afm_scan_grid.n.z-afm_nconv)afm_iz=0;  renderAFM(afm_iz,2); break;
                case SDLK_RIGHTPAREN:   afm_iz--; if(afm_iz<0)afm_iz=afm_scan_grid.n.z-1-afm_nconv; renderAFM(afm_iz,2);  break;
                
                //case SDLK_LESS:      afm_iz++; if(afm_iz>=afm_scan_grid.n.z-afm_nconv)afm_iz=0;  renderAFM(afm_iz,2); break;
                //case SDLK_GREATER:   afm_iz--; if(afm_iz<0)afm_iz=afm_scan_grid.n.z-1-afm_nconv; renderAFM(afm_iz,2);  break;

                //case SDLK_j: makeNonuniformGrid(); break;
                //case SDLK_k: relaxNonuniformGrid(); break;

                case SDLK_k: bViewColorFrag ^= 1; break;
                case SDLK_j: bViewBuilder   ^= 1; break;

                case SDLK_g: W->bGridFF=!W->bGridFF; break;
                //case SDLK_c: W->bOcl=!W->bOcl;       break;
                case SDLK_m: W->swith_method();      break;
                //case SDLK_h: W->ff4.bAngleCosHalf = W->ffl.bAngleCosHalf = !W->ffl.bAngleCosHalf; break;
                //case SDLK_k: bDebug_scanSurfFF ^=1; break;
                //case SDLK_q: W->autoCharges(); break;
                case SDLK_a: bViewAtomSpheres ^= 1; break;
                case SDLK_l: bViewAtomLabels  ^= 1; break;
                case SDLK_t: bViewAtomTypes   ^= 1; break;
                case SDLK_b: bViewBondLabels  ^= 1; break;  
                case SDLK_q: bViewMolCharges  ^= 1; break;
                case SDLK_f: bViewAtomForces  ^= 1; break;
                case SDLK_w: bViewSubstrate   ^= 1; break;
                case SDLK_i: bViewBondLenghts ^= 1; break;
                
                //case SDLK_LEFTBRACKET: rotate( W->selection.size(), &W->selection[0], W->ff.apos, rotation_center, rotation_axis, +rotation_step ); break;
                //case SDLK_RIGHTBRACKET: rotate( W->selection.size(), &W->selection[0], W->ff.apos, rotation_center, rotation_axis, -rotation_step );  break;
                case SDLK_SPACE: 
                    bRunRelax=!bRunRelax;  

                    if( bRunRelax ){ // Update ffl from builder
                        if (bBuilderChanged){
                            W->updateFromBuilder();
                            bindMolecule(W);
                            bBuilderChanged = false;
                        }
                        //bViewBuilder    = false;
                    }else{
                        //bViewBuilder = true;
                        W->updateBuilderFromFF();
                    }

                    //if( bRunRelax ){ if (W->go.bExploring){ W->stopExploring(); }else{ W->startExploring(); }; }

                    // printf( "bRunRelax %i \n", bRunRelax );
                    //if(bRunRelax)W->setConstrains();                  
                    if(!bRunRelax){ if(ogl_MO>0){ int iHOMO = W->getHOMO(); renderOrbital( iHOMO + which_MO );  } }
                    break;

                default:{
                    uint8_t k = event.key.keysym.sym & 0xFF;
                    if( k != event.key.keysym.sym ){ k+=128; } 
                    printf( "free key: 0xFF&(%i)  hex=%x \n", k, event.key.keysym.sym );
                    actions.actionDispatch( this, k );
                } break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
        case SDL_MOUSEBUTTONUP: mouse_default( event ); break;
        case SDL_WINDOWEVENT:{switch (event.window.event) {case SDL_WINDOWEVENT_CLOSE:{ quit(); }break;} } break;
    } // switch( event.type ){
    //printf( "MolGUI::eventMode_default() END bConsole=%i \n", bConsole );
}

void MolGUI::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    if( gui.onEvent( mouseX, mouseY, event ) )return;

    switch (gui_mode){
        case  Gui_Mode::edit: eventMode_edit(event); break;
        case  Gui_Mode::scan: eventMode_scan(event); break;
        case  Gui_Mode::base: 
        default:              eventMode_default(event); break;   
    }
    //AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void MolGUI::keyStateHandling( const Uint8 *keys ){
    if(bConsole){ return; }
    if(gui.bTextEvents){return;}
    double dstep=0.025;
    switch (gui_mode){
        case Gui_Mode::edit: 
        case Gui_Mode::scan: 
        case Gui_Mode::base:
        default:{
            //if( keys[ SDL_SCANCODE_X ] ){ cam.pos.z +=0.1; }
            //if( keys[ SDL_SCANCODE_Z ] ){ cam.pos.z -=0.1; }
            //if( keys[ SDL_SCANCODE_LEFT  ] ){ qCamera.dyaw  (  keyRotSpeed ); }
            //if( keys[ SDL_SCANCODE_RIGHT ] ){ qCamera.dyaw  ( -keyRotSpeed ); }
            //if( keys[ SDL_SCANCODE_UP    ] ){ qCamera.dpitch(  keyRotSpeed ); }
            //if( keys[ SDL_SCANCODE_DOWN  ] ){ qCamera.dpitch( -keyRotSpeed ); }
            //if( keys[ SDL_SCANCODE_A ] ){ cam.pos.add_mul( cam.rotMat().a, -cameraMoveSpeed ); }
            //if( keys[ SDL_SCANCODE_D ] ){ cam.pos.add_mul( cam.rotMat().a,  cameraMoveSpeed ); }
            //if( keys[ SDL_SCANCODE_W ] ){ cam.pos.add_mul( cam.rotMat().b,  cameraMoveSpeed ); }
            //if( keys[ SDL_SCANCODE_S ] ){ cam.pos.add_mul( cam.rotMat().b, -cameraMoveSpeed ); }
            //if( keys[ SDL_SCANCODE_Q ] ){ cam.pos.add_mul( cam.rotMat().c, -cameraMoveSpeed ); }
            //if( keys[ SDL_SCANCODE_E ] ){ cam.pos.add_mul( cam.rotMat().c,  cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_KP_4 ] ){ W->nbmol.shift( {-0.1,0.,0.} ); }
            if( keys[ SDL_SCANCODE_KP_6 ] ){ W->nbmol.shift( {+0.1,0.,0.} ); }
            if( keys[ SDL_SCANCODE_KP_8 ] ){ W->nbmol.shift( {0.,+0.1,0.} ); }
            if( keys[ SDL_SCANCODE_KP_2 ] ){ W->nbmol.shift( {0.,-0.1,0.} ); }
            if( keys[ SDL_SCANCODE_KP_7 ] ){ W->nbmol.shift( {0.,0.,+0.1} ); }
            if( keys[ SDL_SCANCODE_KP_9 ] ){ W->nbmol.shift( {0.,0.,-0.1} ); }
            if( keys[ SDL_SCANCODE_LEFT  ] ){ cam.shift( cam.rotMat().a* -cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_RIGHT ] ){ cam.shift( cam.rotMat().a*  cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_UP    ] ){ cam.shift( cam.rotMat().b*  cameraMoveSpeed ); }
            if( keys[ SDL_SCANCODE_DOWN  ] ){ cam.shift( cam.rotMat().b* -cameraMoveSpeed ); }
            //AppSDL2OGL_3D::keyStateHandling( keys );
        } break;   
    }
};


#endif