
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "VecN.h"

/*
#include "Multipoles.h"
#include "PotentialFlow.h"
#include "grids3D.h"
#include "MultipoleGrid.h"
*/

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"
#include "Draw3D_Molecular.h"

//#include "MMFF.h"
//#include "RARFF.h"
//#include "RARFF2.h"
//#include "RARFFarr.h"
#include "RARFF_SR.h"
#include "QEq.h"

#include "argparse.h"
#include "Console.h"
LambdaDict funcs;
#ifdef WITH_LUA
//#include "LuaUtils.h"
//#include "LuaClass.h"
#include "LuaHelpers.h"
lua_State  * theLua=0;
#endif // WITH_LUA

#define R2SAFE  1.0e-8f


// ============= This is From LIB ReactiveFF.cpp
// /home/prokop/git/ProbeParticleModel_OclPy3/cpp/ReactiveFF.cpp


// ============= Global Variables

std::vector<RigidAtomType*> atomTypes;
RARFF_SR ff;
PairList pairList;
QEq      qeq;

// ============= Functions

void generate_atoms( RARFF_SR& ff, int natom, double xspan, double step, int ntyps, RigidAtomType* types, Quat4i capsBrush ){
    //Vec3d pmin={-5.,-5.,-1.0};
    //Vec3d pmin={-5.,-5.,-1.0};
    ff.map.setup_Buckets3D( Vec3d{-xspan,-xspan,-step}, Vec3d{xspan,xspan,step}, step );
    int nat=natom;
    //int nat=2;
    ff.realloc( nat, 100 );
    for(int i=0; i<nat; i++){
        //if(randf()>0.5){ ff.types[i]=&type1;  }else{ ff.types[i]=&type2; }
        //ff.types[i]=curType;
        ff.types[i]=&types[0];
        ff.apos [i].fromRandomBox( ff.map.pos0 , ff.map.pmax );
        ff.qrots[i].setRandomRotation();
        ((Quat4i*)ff.bondCaps)[i]=capsBrush;
    }
    ff.apos [0]=Vec3dZero;
    ff.qrots[0]=Quat4dIdentity;
    ff.cleanAux();
}

void generate_quad_row(RARFF_SR& ff, int ntypes, RigidAtomType* types, Quat4i capsBrush){
    //Quat4d rots  [] = { Quat4dFront, Quat4dBack, Quat4dFront, Quat4dBack }; 
    //Quat4d y90{0.0,1.0,0.0,0.0};
    // ---- Rotation by x-axis
    // Quat4d q0 {  0.00000000000, 0.0,0.0, 1.00000000000  };
    // Quat4d x30{ -0.25881904510, 0.0,0.0, 0.96592582628  };
    // Quat4d x45{ -0.38268343236, 0.0,0.0, 0.92387953251  };
    // Quat4d x60{ -0.50000000000, 0.0,0.0, 0.86602540378  };
    // Quat4d x90{ -M_SQRT1_2    , 0.0,0.0, M_SQRT1_2      };
    // Quat4d q0 {  0.00000000000, 0.0,0.0, 1.00000000000  };
    // Quat4d x30{ +0.25881904510, 0.0,0.0, 0.96592582628  };
    // Quat4d x45{ +0.38268343236, 0.0,0.0, 0.92387953251  };
    // Quat4d x60{ +0.50000000000, 0.0,0.0, 0.86602540378  };
    // Quat4d x90{ +M_SQRT1_2    , 0.0,0.0, M_SQRT1_2      };
    // ----- rotation by y-axis
    Quat4d q0 { 0.0,  0.00000000000, 0.0, 1.00000000000  };
    Quat4d y30{ 0.0, +0.25881904510, 0.0, 0.96592582628  };
    Quat4d y45{ 0.0, +0.38268343236, 0.0, 0.92387953251  };
    Quat4d y60{ 0.0, +0.50000000000, 0.0, 0.86602540378  };
    Quat4d y90{ 0.0, +M_SQRT1_2    , 0.0, M_SQRT1_2      };
    //Quat4d y90{0.0,1.0,0.0,0.0};
    int np = 5;
    Quat4d rots  [np] = { q0, y30, y45, y60, y90 }; 
    ff.realloc( np, 100 );
    for (int i=0; i<np; i++){
        ff.types[i] = &types[0];
        ff.apos [i] = Vec3d{i*2.0,0.0,0.0};
        ff.qrots[i] = rots[i];
        ((Quat4i*)ff.bondCaps)[i] = capsBrush;
    }
}

void generateDiamond111Surface( RARFF_SR& ff, double alat, Vec3i n, int ntypes, RigidAtomType* types, Quat4i capsBrush ) {
    printf( "generateDiamond111Surface() \n" );
    //n=Vec3i{3,3,4};
    n=Vec3i{5,5,4};
    //double a = 5.43; // Silicon lattice constant, use 3.57 for diamond
    alat = 3.57; // Diamond lattice constant
    double sqrt3 = sqrt(3.0);
    // Hexagonal close-packed layer spacing
    //double layerSpacing = alat * sqrt3 / 2.0;
    double bL             = alat * sqrt3 / 4.0; // Bond length
    printf( "generateDiamond111Surface() alat %g bL %g \n", alat, bL );
    int size = n.x; // Number of layers

    Vec3d avec = Vec3d{  3.0/4.0, sqrt3/4.0, 0.0 }*alat;
    Vec3d bvec = Vec3d{  0.0    , sqrt3/2.0, 0.0 }*alat;
    Vec3d zvec = Vec3d{  0.0    , 0.0      , 1.0 };

    Quat4d q0 {  0.00000000000, 0.0,0.0, 1.00000000000  };
    Quat4d x30 {   0.25881904510, 0.0,0.0, 0.96592582628  };
    Quat4d x30m{ -0.25881904510, 0.0,0.0, 0.96592582628  };
    //Quat4d y60{ 0.0, +0.50000000000, 0.0, 0.86602540378  };
    Quat4d y90 { 0.0, +M_SQRT1_2    , 0.0, M_SQRT1_2      };
    Quat4d y90m{ 0.0, -M_SQRT1_2    , 0.0, M_SQRT1_2      };

    Quat4d rots  [] = { y90%x30, y90m%x30, y90%x30m, y90m%x30m }; 
    Vec3d  shifts[] = { {0.0,0.0,0.0}, (avec+bvec)*(1./3.)  + zvec*(bL/3.0), (avec+bvec)*(1./3.) + zvec*(bL*(4/3.)), Vec3d{0.0,0.0,bL*(5/3.)} };

    bool fixed[] = { false, false, true, true };
    //bool fixed[] = { false, false, false, false };

    Vec3d pmin = Vec3dmax;
    Vec3d pmax = Vec3dmin;
    int i = 0;
    ff.realloc( n.totprod(), 100 );
    // Generate points in a hexagonal pattern
    for (int iz=0; iz<n.z; iz++){
        Vec3d  shift = shifts[iz];
        Quat4d rot   = rots  [iz];
        for (int iy=0; iy<n.y; iy++){
            for (int ix=0; ix<n.x; ix++){
                Vec3d p = avec*ix + bvec*iy + shift;
                pmin = pmin.setIfLower(p);
                pmax = pmax.setIfGreater(p);
                ff.types[i] = &types[0];
                ff.apos [i] = avec*ix + bvec*iy + shift;
                ff.qrots[i] = rot;
                ff.fixedAtoms[i] = fixed[iz];
                ((Quat4i*)ff.bondCaps)[i] = capsBrush;
                i++;
            }
        }
    }

    double step = 2*bL;
    pmin.add( -step, -step, -step );
    pmax.add(  step,  step,  step );

    ff.map.setup_Buckets3D( pmin, pmax, step );
    printf( "generateDiamond111Surface() pmin(%g,%g,%g) pmax(%g,%g,%g) ntot=%i n(%i,%i,%i) \n", pmin.x,pmin.y,pmin.z, pmax.x,pmax.y,pmax.z, ff.map.ntot, ff.map.n.x,ff.map.n.y,ff.map.n.z );
    printf( "generateDiamond111Surface() DONE\n" );
}

template< typename FF>
int pickBond( FF ff, Vec3d& ray0, const Vec3d& hRay, double R ){
    double tmin =  1e+300;
    int imin    = -1;
    for(int ia=0; ia<ff.natom; ia++){
        if(ff.ignoreAtoms)if(ff.ignoreAtoms[ia])continue;
        int nb = ff.types[ia]->nbond;
        for(int j=0; j<nb; j++){
            int i = ia*N_BOND_MAX + j;
            Vec3d pbond = ff.bondPos(i,1.0);
            double ti = raySphere( ray0, hRay, R, pbond );
            //printf( "atom %i t %g %g hRay(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, tmin, hRay.x,hRay.y,hRay.z,  ps[i].x, ps[i].y, ps[i].z );
            //printf( "atom %i t %g %g %g ray0(%g,%g,%g) ps(%g,%g,%g) \n", i, ti, R, tmin, ray0.x,ray0.y,ray0.z,  ps[i].x, ps[i].y, ps[i].z );
            if(ti<tmin){ imin=i; tmin=ti; }
        }
    }
    return imin;
}

void testEF( RARFF_SR& ff, double rmin, double rmax, int n, double* Eout=0, double* Fout=0, double* xs=0, double* Fnum=0 ){
    printf( "RFFapp::testEF() \n" );
    int i1=0;
    int i2=1;
    ff.apos[i1]=Vec3dZero;
    double dx = (rmax-rmin)/n;
    RigidAtomType pairType;
    //printf( "RFFapp::testEF() 2 \n" );
    for(int i=0; i<n; i++){
        double x = rmin+dx*i;
        ff.apos[i2]=Vec3d{ x, 0.0,0.0 };
        pairType.combine( *ff.types[i1], *ff.types[i2] );
        ff.aforce[i2] = Vec3dZero;
        //printf( "RFFapp::testEF()[%i] x %g -> ", i, x );
        double E = ff.pairEF( i1, i2, ff.types[i1]->nbond, ff.types[i1]->nbond, pairType);
        //printf( "RFFapp::testEF() E %g F(%g,%g,%g) \n", E, ff.aforce[i2].x,ff.aforce[i2].y,ff.aforce[i2].z );
        if(Eout) Eout[i]=E;
        if(Fout) Fout[i]=ff.aforce[i2].x;
        if(xs) xs[i]=x;
        if(Fnum){ if(i<2) continue;
            Fnum[i-1] = (Eout[i]-Eout[i-2])/(-2*dx);
            //Fnum[i-1] -= Fout[i-1];  Fnum[i-1]*=1000; 
        }
    }
    if(Fnum){ Fnum[0]=NAN; Fnum[n-1]=NAN; }
    //exit(0);
}
class TestAppRARFF: public AppSDL2OGL_3D { public:

    //RigidAtom     atom1;
    //  type1,type2, typH2O;
    RigidAtomType  typeList [4];
    //RigidAtomType* typeList_[3];
    RigidAtomType* curType=0;

    CapType capTypes[2];

    Quat4i capsBrush;

    int renderMode=0;
    int nrenderModes=2;
    bool bViewBox =true;
    bool bViewGrid=false;
    //bool bViewGrid=true;
    bool bBlockAddAtom=false;
    int ipicked = -1;
    Vec3d ray0;
    Vec3d mouse_p0;

    int ogl_sph=0;

    const char* workFileName="data/work.rff";

    Plot2D plot1;
    //bool bRun = true;
    bool bRun = false;
    //int  perFrame = 100;
    //int  perFrame = 50;
    //int  perFrame = 20;
    int  perFrame = 10;

    int      fontTex;
    bool bConsole=false;
    Console console;

    int nAtomCount= 0;
    int nIgnore   = 0;

    // ========== Functions

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

    void simulation();
    void makePotentialPlot();
    void visualize_cells(const bool bDrawPoints=false);
    void visualize_atoms();
    //void generate_atoms( int natom, double xspan, double step );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {
    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    {RigidAtomType& typ=typeList[1]; // C-sp3
        typ.id = 1;
        typ.nbond  = 4;  // number bonds
        typ.rbond0 =  1.5/2;
        typ.aMorse =  1.0;
        typ.bMorse = -0.7;
        typ.bh0s   = (Vec3d*)sp3_hs;
        typ.print();
    }
    curType = &typeList[1];
    capsBrush.set(-1,-1,-1,-1);
    ff.capTypeList = capTypes;
    ff.typeList    = typeList;
    ff.bRepelCaps  =  false; // ignore caps .... good e.g. for water
    ff.bDonorAcceptorCap = true;

    double alat_Si      = 5.43;
    //double alat_Diamond = 3.57;
    //generate_atoms(   1000, 20.0, ff.RcutMax );
    //generate_atoms( ff, 1000, 20.0, ff.RcutMax, 1, &typeList[1], capsBrush );
    //generate_diamond( ff, alat_Diamond, {3,3,1}, 1, &typeList[1], capsBrush );
    generateDiamond111Surface( ff, alat_Si, {3,3,1}, 1, &typeList[1], capsBrush );
    //ff.qrots[1]=Quat4dFront;
    //ff.qrots[0]=Quat4dBack;
    ff.projectBonds();
    //makePotentialPlot();
    Draw3D::makeSphereOgl( ogl_sph, 3, 0.25 );

    ff.map.pointsToCells( ff.natomActive, ff.apos, ff.ignoreAtoms );

    ff.map.bResizing=true;
    console.init( 256, window );
    console.callback = [&](const char* s){ printf( "console.callback(%s)\n", s ); return 0; };
    console.fontTex = fontTex;
#ifdef WITH_LUA
    console.callback = [&](const char* str)->bool{
       lua_State* L=theLua;
        printf( "console.lua_callback: %s\n", str );
        if (luaL_dostring(L, str) != LUA_OK) [[unlikely]]{
            // If there's an error, print it
            //fprintf(stderr, "Error: %s\n", lua_tostring(L, -1));
            printf( "Error: %s\n", lua_tostring(L, -1) );
            lua_pop(L, 1);  // Remove error message from the stack
            return false;
        }
        return true;
    };
#endif // WITH_LUA
}

void TestAppRARFF::simulation(){
    //printf( "TestAppRARFF::simulation() frame=%i AccelType=%i \n", frameCount, ff.AccelType );
    if(ff.AccelType==2){ 
        pairList.makeSRList( ff.natomActive, ff.apos, ff.RcutMax ); 
        pairList.bind( ff.npairs, ff.pairs );
    }
    for(int i=0; i<perFrame; i++){     
        ff.eval();
        //ff.applyForceHarmonic1D( Vec3dZ, 0.0, -1.0); // Press atoms together in z-diraction (like on substrate) 
        if(ipicked>=0){                                // Drag atoms by Mouse
            Vec3d f = getForceSpringRay( ff.apos[ipicked], (Vec3d)cam.rot.c, ray0, -1.0 );
            ff.aforce[ipicked].add( f );
        }
        //ff.moveMDdamp(0.05, 0.9);
        //ff.moveMDdamp(0.05, 0.98);
        ff.moveMDdamp(0.03, 0.98);
    }
}

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    opengl1renderer.clearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    opengl1renderer.enable(GL_DEPTH_TEST);
    opengl1renderer.disable( GL_LIGHTING );
    //printf("frame %i \n", frameCount);
    //if( ff.tryResize( 5, 100, 10) );
    //return;
    // ---------- Simulate
    //bRun = false;
    //perFrame = 10;
    //ff.bGridAccel=false;
    ray0 = (Vec3d)( cam.pos + cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y);
    if(bRun){
        long T=getCPUticks();
        simulation();
        T=getCPUticks()-T;
        //printf( "T %g[MTicks/Step] %g[tick/pair] pairs_eval(%g) pairs_try(%g) n^2(%i)\n", T*1.0e-6/perFrame, T/(ff.n_pairs_evaluated*1.), ff.n_pairs_tried/(perFrame*1.), ff.n_pairs_evaluated/(perFrame*1.), ff.natom*ff.natom );
        ff.n_pairs_tried    =0;
        ff.n_pairs_evaluated=0;
    }else{
        if(ipicked>=0){
            Vec3d dpos = ray0 - ff.apos[ipicked];
            dpos.makeOrthoU( (Vec3d)cam.rot.c );
            ff.apos[ipicked].add(dpos);
        }
    }
    //ff.checkAtomsCollapsed( 0.2, true );
    ff.checkAtomsOut( 100.0, true );
    if(bViewBox && (ff.AccelType==1) ){
        opengl1renderer.color3f(0.7f,0.7f,0.7f); Draw3D::drawBBox( (Vec3f)ff.map.pos0, (Vec3f)ff.map.pmax );
    }
    if(bViewGrid && (ff.AccelType==1) ){
        //if(ff.AccelType==1)
        visualize_cells();
    }
    Draw3D::drawPointCross( ray0, 0.1 );
    if(ipicked>=0) Draw3D::drawLine( ff.apos[ipicked], ray0);
    Draw3D::drawAxis( 1.0);

    // switch(renderMode){
    //     case 0: opengl1renderer.disable(GL_LIGHTING); break;
    //     case 1: opengl1renderer.enable (GL_LIGHTING); break;
    // }
    opengl1renderer.enable (GL_LIGHTING);
    visualize_atoms();
};

void TestAppRARFF::drawHUD(){
    opengl1renderer.pushMatrix();
	opengl1renderer.translatef( 400.0,400.0,0.0 );
	opengl1renderer.scalef    ( 40.0,40.0,1.0  );
	plot1.view();
    opengl1renderer.popMatrix();

    if(bConsole) console.draw();
}

void TestAppRARFF::keyStateHandling( const Uint8 *keys ){
    if(bConsole){ return; }
    if( keys[ SDL_SCANCODE_R ] ){
        if(ipicked>=0){
            //Mat3d rot;
            Quat4d qrot;  qrot.f=(Vec3d)cam.rot.c*0.01; qrot.normalizeW();
            ff.qrots[ipicked].qmul(qrot);
            ff.projectAtomBons(ipicked);
        }
    }
	AppSDL2OGL_3D::keyStateHandling( keys );
};

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );

    Quat4i* caps = &capsBrush;
    if(ipicked>=0) caps=((Quat4i*)ff.bondCaps)+ipicked;

    switch( event.type ){
        case SDL_KEYDOWN :
            if (bConsole){ bConsole=console.keyDown( event.key.keysym.sym ); }
            else switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;
                case SDLK_BACKQUOTE:{ bConsole = !bConsole;}break;   // ` SDLK_ for key '`' 
                case SDLK_g:  ff.AccelType = (ff.AccelType+1)%3; break;
                case SDLK_h:  bViewGrid = !bViewGrid; break;
                case SDLK_r:  renderMode=(renderMode+1)%nrenderModes; break;
                
                case SDLK_KP_0: caps->array[0]*=-1; break;
                case SDLK_KP_1: caps->array[1]*=-1; break;
                case SDLK_KP_2: caps->array[2]*=-1; break;
                case SDLK_KP_3: caps->array[3]*=-1; break;
                case SDLK_KP_MULTIPLY: caps->mul(-1); printf("capsBrush(%i,%i,%i,%i)\n",  caps->x,caps->y,caps->z,caps->w ); break;

                case SDLK_l: ff.load(workFileName); bRun=false; break;
                case SDLK_k: ff.save(workFileName);  break;
                case SDLK_x: ff.saveXYZ("data/RARFF_out.xyz");  break;
                case SDLK_SPACE: bRun = !bRun;  break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEWHEEL:{
            if     (event.wheel.y > 0){ zoom/=1.2; }
            else if(event.wheel.y < 0){ zoom*=1.2; }}break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:{
                    int ip = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    bBlockAddAtom=false;
                    if( ip>=0){
                        if(ipicked==ip){ ipicked=-1; bBlockAddAtom=true; printf("inv\n"); }
                        else           { ipicked=ip;                     printf("set\n"); }
                    }else{ ipicked=-1; }
                    mouse_p0 = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                    //printf( "LMB DOWN picked %i/%i bblock %i \n", ipicked,ip, bBlockAddAtom );
                    }break;
                case SDL_BUTTON_RIGHT:{
                    /*
                    ipicked = pickParticle( ray0, (Vec3d)cam.rot.c, 0.5, ff.natom, ff.apos, ff.ignoreAtoms );
                    printf( "remove atom %i \n", ipicked );
                    ff.ignoreAtoms[ ipicked ] = true;
                    */
                    }break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //printf( "LMB UP picked %i bblock %i \n", ipicked, bBlockAddAtom );
                    if( (ipicked==-1)&&(!bBlockAddAtom) ){
                        Vec3d mouse_p = (Vec3d)( cam.rot.a*mouse_begin_x + cam.rot.b*mouse_begin_y );
                        Vec3d dir = mouse_p-mouse_p0;
                        if(dir.norm2()<1e-6){ dir=(Vec3d)cam.rot.b; };
                        int ib = pickBond( ff, mouse_p0, (Vec3d)cam.rot.c, 0.5 );
                        if(ib>=0){
                            printf( "add atom to bond %i of atom %i \n", ib%N_BOND_MAX, ib/N_BOND_MAX );
                            Vec3d p0 = ff.bondPos(ib, 2.0 );
                            //ff.inserAtom( {nbBrush,4,4}, mouse_p0, dir, (Vec3d)cam.rot.b );
                            //ff.inserAtom( &type1, (const int[]){0,0,0,0}, p0, dir, (Vec3d)cam.rot.b  );
                            int ia = ff.inserAtom( curType, (int*)&capsBrush, p0, ff.hbonds[ib]*-1, dir  );
                            ff.projectAtomBons(ia);
                            bRun = 0;
                        }
                        //ipicked = -1;
                    }
                    break;
                case SDL_BUTTON_RIGHT:
                    ipicked = -1;
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

void TestAppRARFF::makePotentialPlot(){
    DataLine2D * line_Er = new DataLine2D(100); line_Er->clr = 0xFF0000FF; 
    DataLine2D * line_Fr = new DataLine2D(100); line_Fr->clr = 0xFFFF0000;   line_Fr->replace_xs( line_Er->xs );
    DataLine2D * line_Fn = new DataLine2D(100); line_Fn->clr = 0xFF008000;   line_Fn->replace_xs( line_Er->xs );
    testEF( ff, 0.0, 6.0, 60, line_Er->ys, line_Fr->ys,  line_Er->xs, line_Fn->ys );
    plot1.init();
    plot1.scaling.y = 1.0;
    plot1.fontTex = fontTex;
    plot1.clrGrid = 0xFF404040;
    //plot1.clrBg   = 0xFF408080;
    //plot1.lines.push_back( line1  );
    plot1.lines.push_back( line_Er  );
    plot1.lines.push_back( line_Fr  );
    plot1.lines.push_back( line_Fn  );
    plot1.render();
}

void TestAppRARFF::visualize_cells( bool bDrawPoints ){
    //print( "visualize_cells() \n" );
    bDrawPoints=true;
    for(int ic=0;ic<ff.map.ncell;ic++){
        Vec3i ip;
        Draw  ::color_of_hash( 464+645*ic );
        ff.map.i2ixyz( ic, ip );
        Draw3D::drawBBox( ff.map.box2pos2( ip, {0.,0.,0.} ),  ff.map.box2pos2( ip, {1.,1.,1.} ) );
        if(bDrawPoints){
            int i0=ff.map.cellI0s[ic];
            int ni=ff.map.cellNs [ic];
            //printf( "ic %i i0 %i ni %i \n", ic, i0, ni );
            if( (ni<0) || (i0<0) || ( (ni!=0) && (i0+ni>ff.map.nobjSize) ) )[[unlikely]]{ printf( "ERROR TestAppRARFF::visualize_cells() ic %i i0 %i ni %i ff.map.nobjSize=%i \n", ic, i0, ni, ff.map.nobjSize ); exit(0); }
            for(int j=i0; j<i0+ni;j++){
                int io  = ff.map.cell2obj[j];
                Vec3d p = ff.apos[io];
                //printf( "j %i io %i p(%g,%g,%g) \n", j, io, p.x,p.y,p.z );
                //Draw  ::color_of_hash( 464+645*ic );
                Draw3D::drawPointCross( p, 0.2 );            
            }
        }
    }
}

void TestAppRARFF::visualize_atoms(){
    int na = 0;
    // ---------- Draw
    opengl1renderer.lineWidth(3.0f); 
    opengl1renderer.color3f(0.0,0.0,0.0);
    double fsc = 0.1;
    double tsc = 0.1;
    //printf( "ff.natom %i \n", ff.natom );
    uint32_t clrs[]={0xFF000000,0xFF0000FF,0xFF00FF00,0xFFFF0000};
    for(int ia=0; ia<ff.natom; ia++){
        if(ff.ignoreAtoms[ia])continue;
        na++;
        //opengl1renderer.color3f(0.3,0.3,0.3);
        opengl1renderer.color3f(0.5,0.5,0.5);
        Draw3D::drawShape( ogl_sph, ff.apos[ia], Mat3dIdentity*0.7 );
        for(int j=0; j<ff.types[ia]->nbond; j++){
            int i=ia*N_BOND_MAX+j;
            Vec3d pb = ff.bondPos( i );
            //printf( "bondCaps[%i] %i\n", i, ff.bondCaps[i] );
            if( ff.bondCaps[i]>=0 ){ opengl1renderer.color3f(1.0,0.0,0.0); } else{ opengl1renderer.color3f(0.0,0.0,0.0); }
            //Draw::setRGB( clrs[j] );
            opengl1renderer.color3f(0.5,0.5,0.5);;
            Draw3D::drawLine( ff.apos[ia] , pb );
            opengl1renderer.color3f(0.0,1.0,0.0); Draw3D::drawVecInPos( ff.fbonds[i]*fsc, pb );
            //opengl1renderer.color3f(0.0,0.0,0.0); Draw3D::drawVecInPos( ff.hbonds[i], ff.apos[i] );
            //opengl1renderer.color3f(0.0,1.0,0.0); Draw3D::drawVecInPos( ff.fbonds[io]*fsc, ff.apos[i]+ff.hbonds[io] );
        }
    };
    opengl1renderer.lineWidth(1.0f); 
    //if(na<nAtomCount){ printf( "TestAppRARFF::visualize_atoms() number of atoms decreased( na=%i nAtomCount=%i )\n", na, nAtomCount ); nAtomCount=na; }
    //printf( "TestAppRARFF::visualize_atoms() na=%i \n", na );
    //int nig = ff.countIgnoreAtoms(true);
    //int nig = ff.countIgnoreAtoms(false);
    //if( nig>nIgnore) { printf( "TestAppRARFF::visualize_atoms() number of atoms decreased( na=%i nAtomCount=%i )\n", nig,  nIgnore ); nIgnore=nig; }
}

// ===================== MAIN

TestAppRARFF* app=0;

#ifdef WITH_LUA
int l_insertQuickCommand(lua_State *L){
    const char* s = Lua::getString(L,1);
    printf( "l_insertQuickCommand `%s`\n", s );
    app->console.quick_tab.table.push_back( s );
    app->console.quick_tab.sort();
    printf( "l_insertQuickCommand table.size() %i \n", app->console.quick_tab.table.size() );
    for(int i=0; i<app->console.quick_tab.table.size(); i++){
        printf( "l_insertQuickCommand[%i] `%s`\n", i, app->console.quick_tab.table[i].c_str() );
    }
    return 1; // number of return values to Lua environment
}

int initMyLua(){
    printf( "initMyLua()\n" );
    theLua         = luaL_newstate();
    lua_State  * L = theLua;
    luaL_openlibs(L);
    // lua_register(L, "fix",   l_fixAtom      );
    // lua_register(L, "natom", l_getAtomCount );
    // lua_register(L, "apos",  l_getAtomPos   );
    // lua_register(L, "run",   l_toggleStop   );
    lua_register(L, "command", l_insertQuickCommand  );
    printf( "initMyLua() DONE\n" );
    return 1;
}
#endif // WITH_LUA


int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
    SDL_DisplayMode DM;
    SDL_GetCurrentDisplayMode(0, &DM);
    app = new TestAppRARFF( junk , DM.w-100, DM.h-100 );
	//app = new TestAppRARFF( junk , 1400, 1000 );

    #ifdef WITH_LUA
        initMyLua();
        funcs["-lua"]={1,[&](const char** ss){ if( Lua::dofile(theLua,ss[0]) ){ printf( "ERROR in funcs[-lua] no such file %s => exit()\n", ss[0] ); exit(0); }; }};
    #endif // WITH_LUA
    //char str[40];  sprintf(str,  );
	//SDL_SetWindowTitle( app->child_windows[1]->window, " test_eFF " );
    process_args( argc, argv, funcs );

	app->loop( 1000000 );
	return 0;
}



