
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>

#include <SDL2/SDL.h>


#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "DynamicOpt.h"
#include "RBodyConfDyn.h"

#include "Draw3D.h"
#include "AppSDL2OGL_3D.h"

constexpr int nsubstrate = 8;
double substrate_[8*3] = {
    +1.0,-1.0,-1.0,
    -1.0,-1.0,-1.0,
    +1.0,+1.0,-1.0,
    -1.0,+1.0,-1.0,
    +1.0,-1.0,+1.0,
    -1.0,-1.0,+1.0,
    +1.0,+1.0,+1.0,
    -1.0,+1.0,+1.0
 };
Vec3d * substrate = (Vec3d*)substrate_;

double externalForces(int n, double* pose, double* force){
    //Vec3d pos = *(Vec3d*)pos_;
    //force.set_mul( )
    constexpr double k=-0.2;
    //for(int i=0; i<3; i++){ force[i]+=k*pos[i]; }   // parabolic well on positions

    double E = 0;
    Vec3d   pos = *(Vec3d*)pose;
    Vec3d* fpos =  (Vec3d*)force;
    for(int i=0; i<nsubstrate; i++){
        Vec3d dpos; dpos.set_sub(pos, substrate[i]);
        double ir2 = 1/( 0.25*dpos.norm2() + 0.1);
        fpos->add_mul( dpos, (ir2*ir2 - ir2)  );
    }
    return E;
}


void plotFFplan( int na, int nb, const Vec3d& da, const Vec3d& db, Vec3d pos0, double fscale ){
    for(int ia=0; ia<na; ia++){
        for(int ib=0; ib<nb; ib++){
            Vec3d p = pos0 + (da*ia)+(db*ib);
            Vec3d f; f.set(0.0);
            externalForces( 3, (double*)&p, (double*)&f);
            Draw3D::drawLine( p, p+f*(fscale), COLOR_BLACK );
            Draw3D::drawPointCross( p, 0.01 );
        }
    }
}




// ======================  TestApp

class TestConfDynamics : public AppSDL2OGL_3D {
	public:

	RBodyConfDyn confWorld;
	int vobConf,vobSubstrate,vobFF;

	// ---- function declarations

	virtual void draw();

	TestConfDynamics( int& id, int WIDTH_, int HEIGHT_ );

};

TestConfDynamics::TestConfDynamics( int& id, int WIDTH_, int HEIGHT_ ) : AppSDL2OGL_3D( id, WIDTH_, HEIGHT_ ) {

	vobConf=opengl1renderer.genLists(1);
	opengl1renderer.newList( vobConf, GL_COMPILE );
        opengl1renderer.disable ( GL_LIGHTING );
        Draw3D::drawAxis(0.2);
	opengl1renderer.endList();

    vobSubstrate=opengl1renderer.genLists(1);
	opengl1renderer.newList( vobSubstrate, GL_COMPILE );
        opengl1renderer.enable( GL_LIGHTING );
        for(int i=0; i<nsubstrate;i++){
            Draw3D::drawSphere_oct(2,0.5, substrate[i]);
        }
	opengl1renderer.endList();

    vobFF=opengl1renderer.genLists(1);
	opengl1renderer.newList( vobFF, GL_COMPILE );
        opengl1renderer.disable ( GL_LIGHTING );
        opengl1renderer.color3f(0.7f,0.7f,0.7f);
        plotFFplan( 40, 40, {0.0,0.1,0.0}, {0.1,0.0,0.0}, {-2.0,-2.0,0.0}, 0.05 );
	opengl1renderer.endList();

	confWorld.init( 200 );
	confWorld.pos->set(2.0,2.0,0.0);

	//confWorld.storeThisConf();
    //confWorld.confs[0].print();

    confWorld.objectiveFuncDerivs = externalForces;

	confWorld.mutate_near();
    ((RBodyPose*)confWorld.optimizer.pos)->print();

}

void TestConfDynamics::draw   (){
    opengl1renderer.clearColor( 0.5f, 0.5f, 0.5f, 1.0f );
	opengl1renderer.clear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	Draw3D::drawShape( vobConf, *confWorld.pos, *confWorld.rot );


    double Frest  = 1e+300;
    //delay = 100;
	Frest = confWorld.optStep( );
	//((RBodyPose*)confWorld.optimizer.pos)->print();
	//exit(0);
	if( Frest < 1e-3 ){
        confWorld.storeThisConf();
        confWorld.mutate_near();
	}

	for(int i=0; i<confWorld.nConfs; i++){
        //confWorld.confs[0].print();
        Draw3D::drawShape( vobConf, confWorld.confs[i].pos, confWorld.confs[i].rot );
	}


    //Draw3D::drawShape( {0.0,0.0,0.0},  {0.0,0.0,0.0,1.0}, vobSubstrate );
    opengl1renderer.callList(vobSubstrate);
    opengl1renderer.callList(vobFF);

	//opengl1renderer.disable ( GL_LIGHTING );
	//opengl1renderer.callList( point_cloud );
	//Draw3D::drawAxis ( 3.0f );

};

// ===================== MAIN

TestConfDynamics * testApp;

int main(int argc, char *argv[]){
	SDL_Init(SDL_INIT_VIDEO);
	SDL_GL_SetAttribute(SDL_GL_SHARE_WITH_CURRENT_CONTEXT, 1);
	int junk;
	testApp = new TestConfDynamics( junk , 800, 600 );
	testApp->loop( 1000000 );
	return 0;
}
















