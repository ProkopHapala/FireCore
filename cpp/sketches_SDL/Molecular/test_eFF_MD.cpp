

#include <globals.h>
//int verbosity = 0;

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include  <functional>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>
#include "Draw.h"
#include "Draw3D.h"
#include "Solids.h"

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Vec3Utils.h"
#include "VecN.h"

#include "AppSDL2OGL_3D.h"
#include "testUtils.h"
#include "SDL_utils.h"
#include "Plot2D.h"

//#include "MMFF.h"

#define R2SAFE  1.0e-8f

#include "DynamicOpt.h"

#include "InteractionsGauss.h"
#include "eFF.h"
#include "Forces.h"
#include "eFF_plots.h"

// ======= THE CLASS

class TestAppRARFF: public AppSDL2OGL_3D { public:
    bool bRun = false;
    EFF  ff;
    DynamicOpt opt;
    double dt_MD            = 1e-4;
    double* invMasses_Relax =0;
    double* invMasses_MD    =0;

    int perFrame = 1;
    double F2conv = 1e-6;
    bool bRelaxed=false;
    bool bMD=true;

    std::vector<int> ifixe;

    Vec2i field_ns;
    Vec2d Erange;
    double E0,Espread;
    Vec3d  * field_ps=0;
    double * field_Es=0;
    std::function<void   (const Vec3d& p, Vec3d& f)>  FFfunc;
    std::function<double (const Vec3d& p)          >  Efunc ;

    int ipicked  = 0;

    // DEBUG STUFF
    GLint ogl_fs = 0;
    GLint oglSph = 0;

    Plot2D plot1;

    int  fontTex;

    // ---- Functions

    virtual void draw   ();
    virtual void drawHUD();
    //virtual void mouseHandling( );
    virtual void eventHandling   ( const SDL_Event& event  );
    //virtual void keyStateHandling( const Uint8 *keys );

    TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ );

};

TestAppRARFF::TestAppRARFF( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_, " test_eFF " ) {

    fontTex   = makeTextureHard( "common_resources/dejvu_sans_mono_RGBA_pix.bmp" );
    plot1.fontTex=fontTex;

    DEBUG
    // ===== SETUP GEOM
    //char* fname = "data/H_eFF.xyz";
    //ff.loadFromFile_fgo( "data/H2O.fgo" );    ff.writeTo_fgo( "test.fgo", false ); //exit(0);
    ff.loadFromFile_fgo( "data/H2O_shoot.fgo", true );  ifixe.push_back(ff.ne-1);      bMD=true; bRelaxed=true;
    //ff.loadFromFile_fgo( "data/C2H4.fgo" );
    //ff.loadFromFile_fgo( "data/C2H2.fgo" );

    //ff.bEvalAECoulomb = 0;
    //ff.bEvalAEPauli   = 0;
    //ff.bEvalCoulomb   = 0;
    //ff.bEvalPauli     = 0;
    //ff.bEvalKinetic   = 0;
    //ff.bEvalAA        = 0;

    // ==== Bind Optimizer
    //opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, 0, ff.fDOFs, 0 );
    opt.bindOrAlloc( ff.nDOFs, ff.pDOFs, ff.vDOFs, ff.fDOFs, 0 );
    ff.makeMasses(invMasses_MD);
    if(!bRelaxed){
        opt.cleanVel( );
        opt.initOpt( 0.01, 0.2 );
        opt.f_limit = 1000.0;
    }else{
        opt.invMasses = invMasses_MD;
    }    

    //ff.iPauliModel = 0; // dens overlap
    ff.iPauliModel = 1; // addPauliGauss   from the article using  KRSrho
    //ff.iPauliModel = 2; // addPauliGaussVB valence bons
    ff.info();

    double E = ff.eval();
     ff.printEnergies();

    oglSph=Draw::list(oglSph);
    Draw3D::drawSphere_oct(3,1.0d,Vec3d{0.,0.,0.});
    glEndList();

}

/*
void TestAppRARFF::step_MD(int nStep){
    for(int itr=0;itr<nStep;itr++){
        ff.clearForce();
        Etot = ff.eval();
        opt.move_LeapFrog(double dt_loc);
    }
};
double TestAppRARFF::step_Relax( int nMax=100, double Fconv=1e-4 ){
    double F2conv = Fconv*Fconv;
    for(int itr=0;itr<perFrame;itr++){
        ff.clearForce();
        cleanDebugForce( ff.ne, ff.na);
        Etot = ff.eval();
        F2 = opt.move_FIRE();
        if(F2<F2F2conv) break;
    }
    return F2; 
};
*/

void TestAppRARFF::draw(){
    //printf( " ==== frame %i \n", frameCount );
    glClearColor( 1.0f, 1.0f, 1.0f, 1.0f );
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);

    //return;
    double vminOK = 1e-6;
    double vmaxOK = 1e+3;
    perFrame=10; // ToDo : why it does not work properly for  perFrame>1 ?
    //perFrame = 1;
    double sum = 0;
    if(bRun){
        double F2 = 1.0;
        double Etot;
        for(int itr=0;itr<perFrame;itr++){
            ff.clearForce();
            Etot = ff.eval();
            if(bRelaxed){
                if(bMD){ 
                    opt.move_LeapFrog(dt_MD); // MD-step
                    double vemax = sqrt(maxR2( ff.ne, ((Vec3d*)opt.vel)+ff.na ));
                    printf( "MD E= %g ve_max= %g  \n", Etot, vemax );
                }
            }else{
                for(int ie: ifixe ){  ff.fixElectron(ie, opt.vel); }; // fix fixed electrons
                F2 = opt.move_FIRE();                                 // relaxation step
                printf( "Relax E= %g |F| = %g \n", Etot, sqrt(F2) );
                if(F2<F2conv){                                        // converged?
                    ff.writeTo_fgo( "relaxed.fgo", false );
                    bRelaxed=true; 
                    invMasses_Relax= opt.invMasses;
                    opt.invMasses  = invMasses_MD ;
                    break; 
                }
            }
        }
    }

    drawEFF( ff, oglSph, 1.0, 0.1, 0.1, 1.5 );

    glCallList(ogl_fs);

    char strtmp[256];
    double Qsz = 0.05;
    double fsc = 1.0;
};


void TestAppRARFF::drawHUD(){
	//glTranslatef( 100.0, 250.0, 0.0 );
	//glScalef    ( 100.0, 100.0, 1.0 );
	//plot1.view();
    glTranslatef( 10.0,HEIGHT-20.0,0.0 );
	glColor3f(0.5,0.0,0.3);
    //Draw::drawText( "AHOJ ", fontTex, fontSizeDef, {100,20} );
    int nstr=2048;
	char str[nstr];
	char* s=str;
	s+=ff.Eterms2str(s);
	ff.orbs2str(s);
    Draw::drawText( str, fontTex, fontSizeDef, {120,20} );
}

/*
void TestAppRARFF::mouseHandling( ){
    int mx,my; Uint32 buttons = SDL_GetRelativeMouseState( &mx, &my);
    if ( buttons & SDL_BUTTON(SDL_BUTTON_LEFT)) {
    }
    AppSDL2OGL_3D::mouseHandling( );
};
*/

void TestAppRARFF::eventHandling ( const SDL_Event& event  ){
    //printf( "NonInert_seats::eventHandling() \n" );
    Vec3d pa0;
    switch( event.type ){
        case SDL_KEYDOWN :
            switch( event.key.keysym.sym ){
                case SDLK_p:  first_person = !first_person; break;
                case SDLK_o:  perspective  = !perspective; break;

                case SDLK_i: ff.info(); break;
                //case SDLK_LEFTBRACKET :  Espread *= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                //case SDLK_RIGHTBRACKET:  Espread /= 1.2; ogl_fs=genFieldMap(ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread ); break;
                //case SDLK_e: bMapElectron=!bMapElectron; break;
                // case SDLK_f:{
                //     pa0 = ff.apos[ipicked];
                //     sampleScalarField( Efunc, field_ns, {-5.0,-5.0,+0.1}, {0.1,0.0,0.0}, {0.0,0.1,0.0}, field_ps, field_Es, Erange );
                //     E0 = field_Es[0];
                //     ogl_fs = genFieldMap( ogl_fs, field_ns, field_ps, field_Es, E0-Espread, E0+Espread );
                //     ff.apos[ipicked]= pa0;
                //     }break;
                case SDLK_SPACE: bRun = !bRun; break;
                //case SDLK_r:  world.fireProjectile( warrior1 ); break;
            }
            break;
        case SDL_MOUSEBUTTONDOWN:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ipicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                break;
            }
            break;
        case SDL_MOUSEBUTTONUP:
            switch( event.button.button ){
                case SDL_BUTTON_LEFT:
                    //ibpicked = pickParticle( ff.natoms, ff.apos, ray0, (Vec3d)cam.rot.c , 0.5 );
                    //printf( "dist %i %i = ", ipicked, ibpicked );
                    break;
            }
            break;
    };
    AppSDL2OGL::eventHandling( event );
}

// ===================== MAIN

TestAppRARFF* thisApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	//SDL_SetRelativeMouseMode( SDL_TRUE );
	int junk;
	thisApp = new TestAppRARFF( junk , 800, 600 );
	thisApp->loop( 1000000 );

    //char str[40];  sprintf(str,  );
	//SDL_SetWindowTitle( thisApp->child_windows[1]->window, " test_eFF " );

	return 0;
}
















