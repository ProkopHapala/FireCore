

int verbosity = 0;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include "testUtils.h"

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw3D.h"
#include "SDL_utils.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"
#include "Vec3Utils.h"

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFsp3.h"
#include "NBFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"
#include "SMILESparser.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams

//#include "NBSRFF.h"
#include "IO_utils.h"

#include "AppSDL2OGL_3D.h"

bool bPrint = true;
constexpr int  ntmp_str=2048;
char tmp_str[ntmp_str];

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppMMFFsp3 : public AppSDL2OGL_3D { public:

    int verbosity = 2;

	Molecule    mol;
	MMFFparams  params;
    MMFFsp3    ff;
    NBFF       nff;
    SMILESparser smiles;
    MM::Builder builder;
    DynamicOpt  opt;

    int* atypes = 0;

    bool bNonBonded = true;

    Vec3d cog,vcog,fcog,tq;
    double* bondLenghts=0;

    std::vector<int> selection;
    bool bDragging = false;
    Vec3d  ray0_start;
    Vec3d  rotation_center = Vec3dZero;
    Vec3d  rotation_axis   = Vec3dZ;
    double rotation_step   = 15.0 * (M_PI/180.0);

    Vec3d lvec_a0;
    int icell = 0;
    int frameCountPrev=0;

    bool bConverged = false;
    bool bRunRelax  = false;

    int     fontTex;
    int     ogl_sph=0;
    int     ogl_mol=0;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1, iangPicked = -1;
    int perFrame =  50;

    bool bAtomsSpheres=true;
    double mm_Rsc =  0.25;
    double mm_Rsub = 1.0;
    bool   mm_bAtoms = false;

    double drndv =  10.0;
    double drndp =  0.5;

	virtual void draw   ();
	virtual void drawHUD();
	virtual void eventHandling   ( const SDL_Event& event  );

	TestAppMMFFsp3( int& id, int WIDTH_, int HEIGHT_ );

	//void drawSystem( );
    void drawSystem( bool bAtoms=true, bool bBonds=true, bool bForces=false, float texSize=0.007 );

	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

    void initFromXYZ   (const char* fname);
    void initFromSMILES(const char* s, bool bCap=true);

};

void TestAppMMFFsp3::selectRect( const Vec3d& p0, const Vec3d& p1 ){
    Vec3d Tp0,Tp1,Tp;
    Mat3d rot = (Mat3d)cam.rot;
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

void  TestAppMMFFsp3::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
    //int ib         = builder.rayBonds( ro, rd, 0.3 );
    int ib = pickBondCenter( ff.nbonds, ff.bond2atom, ff.apos, ro, rd, 0.5 );
    printf( "picked bond  %i \n", ib );
    const Vec2i& b = builder.bonds[ib].atoms;
    rotation_axis = (ff.apos[b.b]-ff.apos[b.a]).normalized();
    std::unordered_set<int> innodes1; innodes1.insert( b.a );
    std::unordered_set<int> innodes2; innodes2.insert( b.b );
    MM::splitGraphs( ff.nbonds, ff.bond2atom, ib, innodes1 );
    MM::splitGraphs( ff.nbonds, ff.bond2atom, ib, innodes2 );
    std::unordered_set<int>* sel;
    if( innodes1.size()<innodes2.size()  ){ sel=&innodes1;  rotation_center = ff.apos[b.a]; }else{ sel=&innodes2; rotation_center = ff.apos[b.b]; rotation_axis.mul(-1); }
    selection.clear();
    for( int i:*sel){ selection.push_back(i); };
}

void TestAppMMFFsp3::initFromXYZ(const char* fname){
    int iret=builder.insertFlexibleMolecule( builder.loadMolType( fname, "mol1",0, true ), {0,0,0}, Mat3dIdentity, -1 );
    if(iret<0){ printf("ERROR: Molecule not loaded!!!\n" ); exit(0); }
    //if(verbosity>1)builder.printAtoms ();
    //if(verbosity>1)builder.printConfs ();
    if(verbosity>1)builder.printAtomConfs();
    builder.export_atypes(atypes);
    builder.autoBonds ();               if(verbosity>1)builder.printBonds ();
    //builder.autoBondsPBC();            printf("//======== autoBondsPBC()\n"); //if(verbosity>1)builder.printBonds ();  // exit(0);
    //if(verbosity>1)builder.printAtomConfs();
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});   if(verbosity>1)builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     if(verbosity>1)builder.printAngles();
    //builder.autoAngles( 10.0, 10.0 );  printf("//======== autoAngles()\n"); //if(verbosity>1)builder.printAngles();
    //if(verbosity>1)builder.printAtomConfs();
    //builder.sortConfAtomsFirst();      //printf("//======== sortConfAtomsFirst()\n");
    //if(verbosity>1)builder.printAtomConfs();
    //builder.makeAllConfsSP();          printf("//======== makeAllConfsSP()\n");
    //if(verbosity>1)builder.printAtomConfs();
    builder.assignAllBondParams( );   builder.printBondParams();
}

void TestAppMMFFsp3::initFromSMILES(const char* s, bool bCap){
    smiles.builder=&builder;
    smiles.parseString( 10000, s );     //printf("DONE: smiles.parseString()     \n");
    if(bCap) builder.addAllCapTopo();   //printf("DONE: builder.addAllCapTopo(); \n");
    if(verbosity>1)builder.printBonds();
    if(verbosity>1)builder.printAtomConfs();
    builder.export_atypes(atypes);
    builder.randomizeAtomPos(1.0);
}

TestAppMMFFsp3::TestAppMMFFsp3( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    verbosity = 0;

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    // ------ using molecule from mol-file does not seem to work now - There is some problem with AtomConfs
    // >> This algorithm assumes all atoms with conf precede atoms without confs in the array
    // >>   ERROR in builder.sortBonds() => exit

    int nheavy = 0;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);
    builder.verbosity = verbosity;

    readMatrix( "common_resources/polymer-2.lvs", 3, 3, (double*)&builder.lvec );
    printf( "builder.lvec \n" ); builder.lvec.print();
    //initFromXYZ("common_resources/Benzene_deriv.xyz");
    //initFromSMILES("C=C");
    //initFromSMILES("C=N");
    //initFromSMILES("C=O");
    initFromSMILES("C=C1NC=CC1CO");

    builder.saveMol( "check.mol" );
    builder.toMMFFsp3( ff );
    
    if(verbosity>1)ff.printNeighs();
    if(verbosity>1)ff.printBonds();

    //ff.doPiPiI  =false;
    //ff.doPiPiT  =false;
    //ff.doPiSigma=false;
    //ff.doAngles =false;

    //ff.printNeighs();
    //ff.printBonds();

    bNonBonded = false;
    if(bNonBonded){   
        nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.fapos, 0, ff.bond2atom );
        builder.export_REQs( nff.REQs );
        if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
            printf( "ERROR: nff.pairMask is not sorted => exit \n" );
            exit(0);
        };
    }else{
        printf( "WARRNING : we ignore non-bonded interactions !!!! \n" );
    }

    opt.bindOrAlloc( ff.nDOFs, ff.DOFs,0, ff.fDOFs, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    // ======== Test before we run
    if(verbosity>1)nff.printAtomParams();
    //ckeckNaN_d( ff.natoms, ff.nneigh_max, ff.Kneighs, "ff.Kneighs" );
    ckeckNaN_d( ff.nbonds,             1, ff.bond_k,  "ff.bond_k"  );

    //ff.doPi = 0;
    double E = ff.eval(true); printf( "ff.eval() E = %g \n", E );
    ff.bDEBUG_plot=true;

    //exit(0);
    //Draw3D::makeSphereOgl( ogl_sph, 3, 1.0 );
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );

    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    //selection.insert( selection.end(), {12, 16, 14, 6, 2, 3,   20,18,31,25,26} );
    //selection.insert( selection.end(), {13,29,30} );
    //splitGraphs( ff.nbonds, ff.bond2atom, 12, 22 );
    //std::unordered_set<int> innodes;
    //innodes.insert(12);
    //innodes.insert(11);
    //MM::splitGraphs( ff.nbonds, ff.bond2atom, 22, innodes );
    //for( int i : innodes ){ selection.push_back(i); }
    printf("TestAppMMFFsp3() DONE \n");
}

void TestAppMMFFsp3::drawSystem( bool bAtoms, bool bBonds, bool bForces, float texSize ){
    if(bBonds){
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLines ( ff.nbonds, (int*)ff.bond2atom, ff.apos );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec );
        //glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, ff.pbcShifts );
        glColor3f(0.0f,0.0f,0.0f);Draw3D::bonds( ff.nbonds, ff.bond2atom, ff.apos ); 
        //glColor3f(0.0f,0.0f,1.0f); Draw3D::bondLabels( ff.nbonds, ff.bond2atom, ff.apos, fontTex, 0.02 );                     
        if(bondLenghts) glColor3f(0.0f,0.0f,1.0f); Draw3D::bondPropertyLabel( ff.nbonds, bondLenghts, ff.bond2atom, ff.apos, 1,0, fontTex, texSize, "%4.2f\0" );
    }
    if(bForces){ glColor3f(1.0f,0.0f,0.0f); Draw3D::vecsInPoss( ff.natoms, ff.fapos, ff.apos, 300.0              ); }
    if(bAtoms){
        //glColor3f(0.0f,0.0f,1.0f); Draw3D::atomPropertyLabel( ff.natoms, (double*)nff.REQs, ff.apos, 3,2, fontTex, 0.02, "%4.2f\0" );
        glColor3f(0.5f,0.0f,0.0f); Draw3D::atomLabels( ff.natoms, ff.apos, fontTex, texSize                     );                     
        //if(bSpheres)Draw3D::atomsREQ  ( ff.natoms, ff.apos,   nff.REQs, ogl_sph, 1.0, 0.25, 1.0 );
        if(bAtomsSpheres)Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, mm_Rsc, mm_Rsub );  
        //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 1.0, 1.0 );      
        //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.5, 1.0 );       
        //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.25, 1.0 );       
    }
}

void TestAppMMFFsp3::draw(){
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( "====== Frame # %i \n", frameCount );
    //cam.qrot.rotate(M_PI*0.01,Vec3fX);

    //glColor3f(0.,0.,0.); drawBonds( builder );
    //glColor3f(0.,0.,0.); drawNeighs( builder );
    //drawNeighs( ff );
    //builder.toMMFFsp3( ff, false );

    if(frameCount==1){
        qCamera.pitch( M_PI );
        //ff.printAtomPos();
        //ff.printBondParams();
        //ff.printAngleParams();
        //ff.printTorsionParams();
        lvec_a0 = builder.lvec.a;
        //printf( "lvec_a0  (%g %g,%g) \n", lvec_a0.x, lvec_a0.y, lvec_a0.z );
    }
    { // bondlengths
        //printf(  "ff.nbonds %i ff.natoms %i \n", ff.nbonds, ff.natoms );
        _allocIfNull( bondLenghts, ff.nbonds );
        evalLenghs(ff.nbonds,ff.bond2atom,ff.apos,bondLenghts);
        //for(int i=0; i<ff.nbonds; i++){ Vec2i b=ff.bond2atom[i]; bondLenghts[i]=ff.apos[b.i].dist(ff.apos[b.j]); printf("bondLength[%i](%i,%i):%g\n", i, b.i,b.j, bondLenghts[i] ); };
    }

    //Draw3D::drawAxis(  10. );

	//ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);

    //perFrame = 1;
    //perFrame = 100;
    //perFrame = 20;
    //bRunRelax = false;
    double Ftol = 1e-6;
    //double Ftol = 1e-2;
    double v_av =0;
    if(bRunRelax){
        builder.lvec.a    = lvec_a0 + Vec3d{-1.0,0.0,0.0};
        for(int itr=0; itr<perFrame; itr++){
            double E=0;
            // --- Eval Forces
            ff.cleanAtomForce();
            E += ff.eval(true);

            //if(bNonBonded){ E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} ); }
            //for(int i=0; i<ff.natoms; i++){ ff.fapos[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) ); }
            if(ipicked>=0){ Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 ); ff.fapos[ipicked].add( f ); };
            
            cog  = average( ff.natoms, ff.apos  );
            vcog = sum    ( ff.natoms, (Vec3d*)opt.vel  );
            fcog = sum    ( ff.natoms, ff.fapos );
            //tq   = torq   ( ff.natoms, ff.apos, ff.fapos, cog );
            tq   = torq   ( ff.natoms, ff.apos, ff.fapos, cog );
            tq.add( ff.evalPiTorq() );

            //printf(  "cog(%g,%g,%g) vcog(%g,%g,%g) fcog(%g,%g,%g) torq (%g,%g,%g)\n", cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, fcog.x,fcog.y,fcog.z, tq.x,tq.y,tq.z );

            // --- Move
            v_av = sqrt( VecN::norm2( opt.n, opt.vel )/opt.n*3 ); 
            double f2;
            //opt.move_MD( 0.05, 0.1*v_av );
            opt.move_FIRE();
            //opt.move_GD( 0.01 );
            //printf( "E %g |F| %g |Ftol %g \n", E, sqrt(f2), Ftol );
            if(f2<sq(Ftol)){
                bConverged=true;
            }

        }
    }
    //printf( "(%i|%i,%i,%i) cog(%g,%g,%g) vcog(%g,%g,%g) fcog(%g,%g,%g) torq (%g,%g,%g)\n", ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0,  cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, fcog.x,fcog.y,fcog.z, tq.x,tq.y,tq.z );
    //printf( "neval Ang %i nevalPiSigma %i PiPiT %i PiPiI %i v_av %g \n", ff.nevalAngles, ff.nevalPiSigma, ff.nevalPiPiT, ff.nevalPiPiI, v_av );

    //drawSystem();
    //glColor3f(0.,0.,0.); drawBonds( ff );
    drawSystem(true,true,false);
    Draw3D::drawNeighs( ff, 0.0 );

    //glColor3f(0.,0.,0.); drawBonds( builder );
    //glColor3f(0.,0.,0.); drawNeighs( builder );
    //drawNeighs( ff );

    /*
    if(bDragging)Draw3D::drawTriclinicBox(cam.rot.transposed(), (Vec3f)ray0_start, (Vec3f)ray0 );
    //Draw3D::drawTriclinicBox(builder.lvec, Vec3dZero, Vec3dOne );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::drawTriclinicBox(builder.lvec.transposed(), Vec3dZero, Vec3dOne );
    //glColor3f(0.6f,0.6f,0.6f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    //glColor3f(0.95f,0.95f,0.95f); Draw3D::plotSurfPlane( (Vec3d){0.0,0.0,1.0}, -3.0, {3.0,3.0}, {20,20} );
    if(builder.bPBC){
        //printf( "draw PBC \n" );
        //Draw3D::drawPBC( (Vec3i){1,1,0}, builder.lvec, [&](){drawSystem();} );
        Draw3D::drawPBC( (Vec3i){2,2,0}, builder.lvec, [&](){drawSystem();} );
    }else{
        drawSystem();
    }
    for(int i=0; i<selection.size(); i++){
        int ia = selection[i];
        glColor3f( 0.f,1.f,0.f );
        Draw3D::drawSphereOctLines( 8, 0.3, ff.apos[ia] );
    }
    if(makeScreenshot){
        saveScreenshot( icell );
        icell++;
    }
    */
};


void TestAppMMFFsp3::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_KP_7: builder.lvec.a.x+=xstep; break;
                //case SDLK_KP_4: builder.lvec.a.x-=xstep; break;
                //case SDLK_KP_8: builder.lvec.a.y+=xstep; break;
                //case SDLK_KP_5: builder.lvec.a.y-=xstep; break;
                //case SDLK_KP_9: builder.lvec.a.z+=xstep; break;
                //case SDLK_KP_6: builder.lvec.a.z-=xstep; break;
                
                case SDLK_o:    bAtomsSpheres=!bAtomsSpheres; break;

                case SDLK_p:    ff.bDEBUG_plot=!ff.bDEBUG_plot; break;
                case SDLK_KP_7: ff.iDEBUG_pick++; if(ff.iDEBUG_pick>=ff.iDEBUG_n)ff.iDEBUG_pick=0; break;
                case SDLK_KP_4: ff.iDEBUG_pick--; if(ff.iDEBUG_pick<0           )ff.iDEBUG_pick=ff.iDEBUG_n-1; break;

                case SDLK_f:
                    //selectShorterSegment( (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y + cam.rot.c*-1000.0), (Vec3d)cam.rot.c );
                    selectShorterSegment( ray0, (Vec3d)cam.rot.c );
                    //selection.erase();
                    //for(int i:builder.selection){ selection.insert(i); };
                    break;

                case SDLK_LEFTBRACKET:
                    rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, +rotation_step );
                    break;
                case SDLK_RIGHTBRACKET:
                    rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, -rotation_step );
                    break;
                case SDLK_SPACE: bRunRelax=!bRunRelax; break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                    ff.iDEBUG_pick=ipicked;
                    printf( "ipicked %i \n", ipicked );
                    /*
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                    selection.clear();
                    if(ipicked>=0){ selection.push_back(ipicked); };
                    printf( "picked atom %i \n", ipicked );
                    */
                    ray0_start = ray0;
                    bDragging = true;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = ff.pickBond( ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf("ibpicked %i \n", ibpicked);
                    break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = -1;
                    //ray0_start
                    if( ray0.dist2(ray0_start)<0.1 ){
                        ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natoms, ff.apos );
                        selection.clear();
                        if(ipicked>=0){ selection.push_back(ipicked); };
                        printf( "picked atom %i \n", ipicked );
                    }else{
                        selectRect( ray0_start, ray0 );
                    }
                    bDragging=false;
                    break;
                case SDL_BUTTON_RIGHT:
                    //ibpicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
    //STOP = false;
}

void TestAppMMFFsp3::drawHUD(){
    glDisable ( GL_LIGHTING );

    glTranslatef( 10.0,HEIGHT-20.0,0.0 );
	glColor3f(0.5,0.0,0.3);
	char* s=str;
    //printf( "(%i|%i,%i,%i) cog(%g,%g,%g) vcog(%g,%g,%g) fcog(%g,%g,%g) torq (%g,%g,%g)\n", ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0,  cog.x,cog.y,cog.z, vcog.x,vcog.y,vcog.z, fcog.x,fcog.y,fcog.z, tq.x,tq.y,tq.z );
    //printf( "neval Ang %i nevalPiSigma %i PiPiT %i PiPiI %i v_av %g \n", ff.nevalAngles, ff.nevalPiSigma, ff.nevalPiPiT, ff.nevalPiPiI, v_av );
    s += sprintf(s, "eval:Ang,ps,ppT,ppI(%i|%i,%i,%i)\n",  ff.nevalAngles>0, ff.nevalPiSigma>0, ff.nevalPiPiT>0, ff.nevalPiPiI>0 );
    s += sprintf(s, "cog (%g,%g,%g)\n", cog .x,cog .y,cog .z);
    s += sprintf(s, "vcog(%15.5e,%15.5e,%15.5e)\n", vcog.x,vcog.y,vcog.z);
    s += sprintf(s, "fcog(%15.5e,%15.5e,%15.5e)\n", fcog.x,fcog.y,fcog.z);
    s += sprintf(s, "torq(%15.5e,%15.5e,%15.5e)\n", tq  .x,tq  .y,tq  .z);
    Draw::drawText( str, fontTex, fontSizeDef, {100,20} );

}

// ===================== MAIN

TestAppMMFFsp3 * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	//thisApp = new TestAppMMFFsp3( junk , 800, 600 );
    thisApp = new TestAppMMFFsp3( junk , 1600, 1000 );
	//thisApp = new TestAppMMFFsp3( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















