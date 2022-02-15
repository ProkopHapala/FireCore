
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

#include "raytrace.h"
#include "Forces.h"

#include "Molecule.h"
#include "MMFFmini.h"
#include "NBFF.h"
#include "MMFFparams.h"
#include "MMFFBuilder.h"
#include "DynamicOpt.h"
#include "QEq.h"

#include "Draw3D_Molecular.h"  // it needs to know MMFFparams

#include "FireCoreAPI.h"

//#include "NBSRFF.h"
#include "IO_utils.h"

#include "AppSDL2OGL_3D.h"

// ===========================================
// ================= MAIN CLASS ==============
// ===========================================

class TestAppMMFFmini : public AppSDL2OGL_3D { public:
	Molecule    mol;
	MMFFparams  params;
    MMFFmini    ff;
    NBFF       nff;
    MM::Builder builder;
    DynamicOpt  opt;
    FireCore::Lib fireCore;

    int* atypes = 0;
    int* atypeZ = 0;

    bool bNonBonded = true;

    std::vector<Molecule*> molecules;

    std::vector<int> selection;
    bool bDragging = false;
    Vec3d  ray0_start;
    Vec3d  rotation_center = Vec3dZero;
    Vec3d  rotation_axis   = Vec3dZ;
    double rotation_step   = 15.0 * (M_PI/180.0);

    Vec3d lvec_a0;
    int icell = 0;
    int frameCountPrev=0;

    //std::unordered_map<std::string,int> atomTypeDict;

    //Mat3d lvec;

    bool bConverged = false;
    bool bRunRelax  = false;

    int  fontTex;
    int  ogl_sph=0;
    int  ogl_mol=0;

    char str[256];

    Vec3d ray0;
    int ipicked  = -1, ibpicked = -1, iangPicked = -1;
    Vec3d* picked_lvec = 0;
    int perFrame =  50;

    double drndv =  10.0;
    double drndp =  0.5;

    

	virtual void draw   ();
	virtual void drawHUD();
	//virtual void mouseHandling( );
	virtual void eventHandling   ( const SDL_Event& event  );
	//virtual void keyStateHandling( const Uint8 *keys );

	TestAppMMFFmini( int& id, int WIDTH_, int HEIGHT_ );

	//int makeMoleculeInline();
	//int makeMoleculeInlineBuilder( bool bPBC );
	int loadMoleculeMol( const char* fname, bool bAutoH, bool loadTypes );
	//int loadMoleculeXYZ( const char* fname, bool bAutoH );
	int loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH=false );

	void drawSystem( );

	void selectShorterSegment( const Vec3d& ro, const Vec3d& rd );
	void selectRect( const Vec3d& p0, const Vec3d& p1 );

	void saveScreenshot( int i=0, const char* fname="data/screenshot_%04i.bmp" );

};

//=================================================
//                   INIT()
//=================================================

TestAppMMFFmini::TestAppMMFFmini( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

    fontTex = makeTexture( "common_resources/dejvu_sans_mono_RGBA_inv.bmp" );

    int nheavy = 0;
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);

    readMatrix( "common_resources/polymer-2.lvs", 3, 3, (double*)&builder.lvec );
    molecules.push_back( new Molecule() ); molecules[0]->atomTypeDict = builder.atomTypeDict; molecules[0]->load_xyz("common_resources/polymer-2.xyz", true);
    //molecules.push_back( new Molecule() ); molecules[1]->atomTypeDict = builder.atomTypeDict; molecules[1]->load_xyz("common_resources/polymer-2-monomer.xyz", true);
    builder.insertFlexibleMolecule(  molecules[0], {0,0,0}           , Mat3dIdentity, -1 );
    //builder.insertFlexibleMolecule(  molecules[1], builder.lvec.a*1.2, Mat3dIdentity, -1 );

    builder.lvec.a.x *= 2.3;

    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    builder.export_atypes(atypes);
    builder.bDEBUG = true;
    //builder.autoBonds ();             builder.printBonds ();
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     builder.printAngles();
    builder.toMMFFmini( ff, &params );
    builder.saveMol( "data/polymer.mol" );

    //builder.lvec.a.x *= 2.0;

    nff.bindOrRealloc( ff.natoms, ff.nbonds, ff.apos, ff.aforce, 0, ff.bond2atom );
    builder.export_REQs( nff.REQs );

    //bNonBonded = false;
    if(bNonBonded){
        if( !checkPairsSorted( nff.nmask, nff.pairMask ) ){
            printf( "ERROR: nff.pairMask is not sorted => exit \n" );
            exit(0);
        };
    }else{
        printf( "WARRNING : we ignore non-bonded interactions !!!! \n" );
    }

    fireCore.loadLib( "/home/prokop/git/FireCore/build/libFireCore.so" );
    atypeZ = new int[ff.natoms];
    for(int i=0; i<ff.natoms; i++){ 
        atypeZ[i]= (atypes[i]==0)?  1 : 6 ;
        printf( "DEBUG atom %i, type %i -> iZ = %i \n", i, atypes[i], atypeZ[i] ); 
    }  // NOTE : This is just temporary hack 
    fireCore.init   ( ff.natoms, atypeZ, (double*)ff.apos );

    opt.bindOrAlloc( 3*ff.natoms, (double*)ff.apos, 0, (double*)ff.aforce, 0 );
    //opt.setInvMass( 1.0 );
    opt.cleanVel( );

    // ======== Test before we run
    nff.printAtomParams();
    ff.printAtomPos();
    ff.printBondParams();
    ff.printAngleParams();
    ff.printTorsionParams();
    double E = ff.eval(true);
    printf( "iter0 E = %g \n", E );
    printf("TestAppMMFFmini.init() DONE \n");
    //exit(0);

    //Draw3D::makeSphereOgl( ogl_sph, 3, 1.0 );
    Draw3D::makeSphereOgl( ogl_sph, 5, 1.0 );

    //float l_diffuse  []{ 0.9f, 0.85f, 0.8f,  1.0f };
	float l_specular []{ 0.0f, 0.0f,  0.0f,  1.0f };
    //glLightfv    ( GL_LIGHT0, GL_AMBIENT,   l_ambient  );
	//glLightfv    ( GL_LIGHT0, GL_DIFFUSE,   l_diffuse  );
	glLightfv    ( GL_LIGHT0, GL_SPECULAR,  l_specular );

    std::unordered_set<int> innodes;
    //innodes.insert(12);
    innodes.insert(11);
    MM::splitGraphs( ff.nbonds, ff.bond2atom, 22, innodes );
    for( int i : innodes ){ selection.push_back(i); }

    picked_lvec = &builder.lvec.a;
}

//=================================================
//                   DRAW()
//=================================================

void TestAppMMFFmini::draw(){
    //glClearColor( 0.5f, 0.5f, 0.5f, 1.0f );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    //printf( "====== Frame # %i \n", frameCount );
    //cam.qrot.rotate(M_PI*0.01,Vec3fX);
    //Draw3D::drawAxis(  10. );
    //if( ogl_mol ){ glCallList( ogl_mol ); return; }
    //printf( "builder.lvec: " ); builder.lvec.print();

    if(frameCount==1){
        qCamera.pitch( M_PI );
        ff.printAtomPos();
        ff.printBondParams();
        ff.printAngleParams();
        ff.printTorsionParams();
        //lvec_a0 = builder.lvec.a;
        //printf( "lvec_a0  (%g %g,%g) \n", lvec_a0.x, lvec_a0.y, lvec_a0.z );
    }

	//ibpicked = world.pickBond( ray0, camMat.c , 0.5 );
    ray0 = (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    Draw3D::drawPointCross( ray0, 0.1 );
    //Draw3D::drawVecInPos( camMat.c, ray0 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);

    double Ftol = 1e-6;
    //double Ftol = 1e-2;

    //ff.apos[0].set(-2.0,0.0,0.0);
    //perFrame = 1;
    //perFrame = 100;
    perFrame = 20;
    //bRunRelax = false;

    bool makeScreenshot = false;
    if(bRunRelax){
        //builder.lvec.a    = lvec_a0 + Vec3d{-1.0,0.0,0.0};

        for(int itr=0; itr<perFrame; itr++){
            double E=0;
            ff.cleanAtomForce();
            E += ff.eval(false);
            if(bNonBonded){
                //E += nff.evalLJQ_sortedMask();   // This is fast but does not work in PBC
                E += nff.evalLJQ_pbc( builder.lvec, {1,1,1} );
            }
            if(ipicked>=0){
                Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 );
                //printf( "f (%g,%g,%g)\n", f.x, f.y, f.z );
                ff.aforce[ipicked].add( f );
            };
            float K = -0.01;
            for(int i=0; i<ff.natoms; i++){
                //ff.aforce[i].add( getForceHamakerPlane( ff.apos[i], {0.0,0.0,1.0}, -3.0, 0.3, 2.0 ) );
                ff.aforce[i].add( getForceMorsePlane( ff.apos[i], {0.0,0.0,1.0}, -5.0, 0.0, 0.01 ) );
                //ff.aforce[i].z += ff.apos[i].z * K;

                //printf( "%g %g %g\n",  world.aforce[i].x, world.aforce[i].y, world.aforce[i].z );
            }

            ff.aforce[  10 ].set(0.0); // This is Hack to stop molecule from moving

            //opt.move_LeapFrog(0.01);
            //opt.move_MDquench();
            //opt.move_GD(0.001);
            double f2 = opt.move_FIRE();

            // =========== Molecular Dynamics with random velicity
            //float d = 0.05;
            //for(int i=0; i<ff.natoms; i++){ ff.aforce[i].add({randf(-d,d),randf(-d,d),randf(-d,d)});  };
            //double f2; opt.move_MD( 0.1, 0.005 );

            //printf( "E %g |F| %g |Ftol %g \n", E, sqrt(f2), Ftol );
            if(f2<sq(Ftol)){
                bConverged=true;
            }

        }
    }

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
        glColor3f( 0.f,1.f,0.f ); Draw3D::drawSphereOctLines( 8, 0.3, ff.apos[ia] );
    }

    if(iangPicked>=0){
        glColor3f(0.,1.,0.);
        Draw3D::angle( ff.ang2atom[iangPicked], ff.ang_cs0[iangPicked], ff.apos, fontTex );
    }

    if(makeScreenshot){ saveScreenshot( icell ); icell++; }

};

void TestAppMMFFmini::selectRect( const Vec3d& p0, const Vec3d& p1 ){
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

void  TestAppMMFFmini::selectShorterSegment( const Vec3d& ro, const Vec3d& rd ){
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

int TestAppMMFFmini::loadMoleculeXYZ( const char* fname, const char* fnameLvs, bool bAutoH ){
    params.loadAtomTypes( "common_resources/AtomTypes.dat" );
    params.loadBondTypes( "common_resources/BondTypes.dat" );
    builder.bindParams(&params);
    int nheavy = builder.load_xyz( fname, bAutoH, true, true );
    readMatrix( fnameLvs, 3, 3, (double*)&builder.lvec );

    //builder.printAtoms ();
    //builder.printConfs ();
    builder.printAtomConfs();
    //exit(0);
    builder.export_atypes(atypes); // NOTE : these are not proton numbers !!!!

    builder.bDEBUG = true;
    //builder.autoBonds ();             builder.printBonds ();
    builder.autoBondsPBC();             builder.printBonds ();  // exit(0);
    //builder.autoBondsPBC(-0.5, 0, -1, {0,0,0});             builder.printBonds ();  // exit(0);
    //builder.autoAngles( 0.5, 0.5 );     builder.printAngles();
    builder.autoAngles( 10.0, 10.0 );     builder.printAngles();

    //bNonBonded = false;
    //exit(0);
    builder.toMMFFmini( ff, &params );

    builder.saveMol( "data/polymer.mol" );

    return nheavy;
}

int TestAppMMFFmini::loadMoleculeMol( const char* fname, bool bAutoH, bool bLoadTypes ){

    /// should look in   test_SoftMolecularDynamics.cpp

    if(bLoadTypes){
        printf( "bLoadTypes==True : load atom and bond types from file \n" );
        params.loadAtomTypes( "common_resources/AtomTypes.dat" );
        params.loadBondTypes( "common_resources/BondTypes.dat");
        //builder.params = &params;
    }else{
        printf( "bLoadTypes==False : make default Atom names dictionary \n" );
        params.initDefaultAtomTypeDict( );
    }
    //mol.atomTypeDict  = &params.atomTypeDict;
    //mol.atomTypeNames = &params.atomTypeNames;
    mol.bindParams(&params);
    mol.loadMol( fname );
    //mol.loadMol_const( "common_resources/propylacid.mol");
    //mol.loadMol_const( "/home/prokop/Dropbox/TEMP/ERC2021/Molecules/chain--frag-4---N-----.mol" );
    //exit(0);
    //int iH = 1;
    int iH = params.atomTypeDict["H"];
    int nh     = mol.countAtomType(iH); printf( "nh %i\n", nh );
    int nheavy = mol.natoms - nh;
    if(bAutoH){
        printf( "bAutoH==True : MMFFBuilder will re-calculate hydrogens, pi-orbitals and free electron pairs \n" );
        builder.insertFlexibleMolecule_ignorH( &mol, Vec3dZero, Mat3dIdentity, iH );
        //builder.setConfs( 0, 0, 0, mol.natoms-nh );
        //for(int i=0;i<(mol.natoms-nh);i++){ builder.makeSPConf(i,0,0); }
        //for(int i=0;i<mol.natoms;i++)     { builder.makeSPConf(i,0,0); }
    }else{
        printf( "bAutoH==False : Angles assigned by simple algorithm Molecule::autoAngles \n" );
        //mol.bondsOfAtoms();
        params.assignREs( mol.natoms, mol.atomType, mol.REQs );
        mol.autoAngles(true);
        Vec3d cog = mol.getCOG_av();
        mol.addToPos( cog*-1.0d );
        builder.insertMolecule(&mol, Vec3dZero, Mat3dIdentity, false );
        builder.toMMFFmini( ff, &params );
    }

    //builder.sortAtomsOfBonds();
    builder.tryAddConfsToAtoms(0, nh);
    builder.tryAddBondsToConfs();
    //for(int i=0; i<nh; i++){ builder.addConfToAtom(i); }
    //builder.tryAddBondsToConfs();

    //mol.printAtomInfo();
    //mol.printAtom2Bond();
    //mol.printAngleInfo();
    builder.printAtoms();
    //builder.printBonds();
    //builder.printAngles();
    //builder.printConfs();

    //bNonBonded = false;      // ToDo : WARRNING : this is just hack, because builder.sortBonds() does not seem to work, we have to switch off Non-Bonding interactions
    builder.trySortBonds();
    //builder.sortBonds();
    builder.printBonds();
    builder.printAngles();
    builder.printConfs();
    builder.toMMFFmini( ff, &params );

    //Draw3D::shapeInPoss( ogl_sph, ff.natoms, ff.apos, 0 );

    /*
    ogl_mol = glGenLists(1);
    glNewList( ogl_mol, GL_COMPILE );
        Draw3D::drawLines( mol.nbonds, (int*)mol.bond2atom, mol.pos );
    glEndList();
    */

    return nheavy;
}

void TestAppMMFFmini::drawSystem( ){
    //glColor3f(0.0f,0.0f,0.0f); Draw3D::drawLines ( ff.nbonds, (int*)ff.bond2atom, ff.apos );
    glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsPBC  ( ff.nbonds, ff.bond2atom, ff.apos, &builder.bondPBC[0], builder.lvec ); // DEBUG
    glColor3f(0.5f,0.0f,0.0f); Draw3D::atomLabels( ff.natoms, ff.apos, fontTex                     );                     //DEBUG
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::bondLabels( ff.nbonds, ff.bond2atom, ff.apos, fontTex, 0.02 );                     //DEBUG
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::atomPropertyLabel( ff.natoms, (double*)nff.REQs, ff.apos, 3,2, fontTex, 0.02, "%4.2f\0" );

    //glColor3f(1.0f,0.0f,0.0f); Draw3D::vecsInPoss( ff.natoms, ff.aforce, ff.apos, 300.0              );
    //Draw3D::atomsREQ  ( ff.natoms, ff.apos,   nff.REQs, ogl_sph, 1.0, 0.25, 1.0 );
    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 1.0, 1.0 );       //DEBUG
    //Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.5, 1.0 );       //DEBUG
    Draw3D::atoms( ff.natoms, ff.apos, atypes, params, ogl_sph, 1.0, 0.25, 1.0 );       //DEBUG
}

void TestAppMMFFmini::saveScreenshot( int i, const char* fname ){
    //if(makeScreenshot){
        char str[64];
        sprintf( str, fname, i );               // DEBUG
        printf( "save to %s \n", str );
        unsigned int *screenPixels = new unsigned int[WIDTH*HEIGHT*4];  //DEBUG
        glFlush();                                                      //DEBUG
        glFinish();                                                     //DEBUG
        //glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_INT, screenPixels);
        glReadPixels(0, 0, WIDTH, HEIGHT, GL_RGBA, GL_UNSIGNED_BYTE, screenPixels);   //DEBUG
        //SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0xff000000, 0x00ff0000, 0x0000ff00, 0x000000ff );   //DEBUG
        SDL_Surface *bitmap = SDL_CreateRGBSurfaceFrom(screenPixels, WIDTH, HEIGHT, 32, WIDTH*4, 0x000000ff, 0x0000ff00, 0x00ff0000, 0xff000000 );   //DEBUG
        SDL_SaveBMP(bitmap, str);    //DEBUG
        SDL_FreeSurface(bitmap);
        delete[] screenPixels;
}

void TestAppMMFFmini::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    float xstep = 0.2;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                //case SDLK_p:  first_person = !first_person; break;
                //case SDLK_o:  perspective  = !perspective; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;

                //case SDLK_LEFTBRACKET:  ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;
                //case SDLK_RIGHTBRACKET: ff.i_DEBUG=(ff.i_DEBUG+1)%ff.nang; break;

                //case SDLK_v: for(int i=0; i<world.natoms; i++){ ((Vec3d*)opt.vel)[i].add(randf(-drndv,drndv),randf(-drndv,drndv),randf(-drndv,drndv)); } break;
                //case SDLK_p: for(int i=0; i<world.natoms; i++){ world.apos[i].add(randf(-drndp,drndp),randf(-drndp,drndp),randf(-drndp,drndp)); } break;

                //case SDLK_LEFTBRACKET:  if(ibpicked>=0) world.bond_0[ibpicked] += 0.1; break;
                //case SDLK_RIGHTBRACKET: if(ibpicked>=0) world.bond_0[ibpicked] -= 0.1; break;

                //case SDLK_a: world.apos[1].rotate(  0.1, {0.0,0.0,1.0} ); break;
                //case SDLK_d: world.apos[1].rotate( -0.1, {0.0,0.0,1.0} ); break;

                //case SDLK_w: world.apos[1].mul( 1.1 ); break;
                //case SDLK_s: world.apos[1].mul( 0.9 ); break;

                //case SDLK_KP_7: builder.lvec.a.mul(   1.01); break;
                //case SDLK_KP_4: builder.lvec.a.mul( 1/1.01); break;

                case SDLK_KP_1: picked_lvec = &builder.lvec.a; break;
                case SDLK_KP_2: picked_lvec = &builder.lvec.b; break;
                case SDLK_KP_3: picked_lvec = &builder.lvec.c; break;

                case SDLK_KP_7: picked_lvec->x+=xstep; break;
                case SDLK_KP_4: picked_lvec->x-=xstep; break;

                case SDLK_KP_8: picked_lvec->y+=xstep; break;
                case SDLK_KP_5: picked_lvec->y-=xstep; break;

                case SDLK_KP_9: picked_lvec->z+=xstep; break;
                case SDLK_KP_6: picked_lvec->z-=xstep; break;

                case SDLK_f:
                    //selectShorterSegment( (Vec3d)(cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y + cam.rot.c*-1000.0), (Vec3d)cam.rot.c );
                    selectShorterSegment( ray0, (Vec3d)cam.rot.c );
                    //selection.erase();
                    //for(int i:builder.selection){ selection.insert(i); };
                    break;

                case SDLK_LEFTBRACKET:
                    Vec3d::rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, +rotation_step );
                    break;
                case SDLK_RIGHTBRACKET:
                    Vec3d::rotate( selection.size(), &selection[0], ff.apos, rotation_center, rotation_axis, -rotation_step );
                    break;

                case SDLK_SPACE: bRunRelax=!bRunRelax; break;

                case SDLK_g: iangPicked=(iangPicked+1)%ff.nang;
                    printf( "ang[%i] cs(%g,%g) k %g (%i,%i,%i)\n", iangPicked, ff.ang_cs0[iangPicked].x, ff.ang_cs0[iangPicked].y, ff.ang_k[iangPicked],
                        ff.ang2atom[iangPicked].a,ff.ang2atom[iangPicked].b,ff.ang2atom[iangPicked].c );
                    break;

            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
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

void TestAppMMFFmini::drawHUD(){
    glDisable ( GL_LIGHTING );

}

// ===================== MAIN

TestAppMMFFmini * thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppMMFFmini( junk , 800, 600 );
	//thisApp = new TestAppMMFFmini( junk , 800, 400 );
	thisApp->loop( 1000000 );
	return 0;
}
















