
#ifndef  Draw3D_Molecular_h
#define  Draw3D_Molecular_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <stdint.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"

#include "Draw.h"
//#include "Mesh.h"
#include "Draw3D.h"

#include "Forces.h"

//#include "MMFFParams.h"

//#include <SDL2/SDL_opengl.h>

namespace Draw3D{

void plotSurfPlane( Vec3d normal, double c0, Vec2d d, Vec2i n ){
    Vec3d da,db;
    normal.getSomeOrtho( da,db );
    da.mul( d.a/da.norm() );
    db.mul( d.b/db.norm() );
    //glColor3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(normal, {0.0,0.0,0.0} );
    //glColor3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos(da*10, {0.0,0.0,0.0} );
    //glColor3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(db*10, {0.0,0.0,0.0} );
    Draw3D::drawRectGridLines( n*2, (da*-n.a)+(db*-n.b) + normal*c0, da, db );
}

void torsion( Quat4i t, const Vec3d* apos ){
    float sz=0.1;
    //Vec2i b=ff.ang2bond[i];
    //Quat4i t=ff.tors2atom[i];
    Draw3D::drawArrow(apos[t.x],apos[t.y],sz);
    Draw3D::drawArrow(apos[t.y],apos[t.z],sz);
    Draw3D::drawArrow(apos[t.z],apos[t.w],sz);
    //if(b.i&SIGN_MASK){ Draw3D::drawArrow(apos[a.x],apos[a.y]); }else{ Draw3D::drawArrow(apos[a.x],apos[a.y]); };
    //if(b.j&SIGN_MASK){   };
}

void makeSphereOgl( int& ogl, int nsub, float sz ){
    ogl = Draw::list(ogl);
    //glNewList( ogl, GL_COMPILE );
        //glEnable( GL_LIGHTING );
        //glColor3f( 0.8f, 0.8f, 0.8f );
        //Draw3D::drawSphere_oct(3, 0.5, {0.0,0.0,0.0} );
        Draw3D::drawSphere_oct( nsub, sz, {0.0,0.0,0.0} );
    glEndList();
}

void atomsREQ( int n, Vec3d* ps, Vec3d* REQs, int ogl_sph, float qsc=1, float Rsc=1, float Rsub=0, bool bPointCross=false ){
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    for(int i=0; i<n; i++){
        float q = (float)REQs[i].z*qsc;
        glColor3f(1-fmax(0,-q),1-fmax(q,-q),1-fmax(0,+q));
        if(bPointCross){
            Draw3D::drawPointCross( ps[i], (REQs[i].x-Rsub)*Rsc );
        }else{
            //Draw3D::drawShape( ogl_sph, ps[i], Mat3dIdentity*(REQs[i].x*Rsc) );
            Draw3D::drawShape( ogl_sph, ps[i], Mat3dIdentity*((REQs[i].x-Rsub)*Rsc) );
        }
    }
}

void bondLabels( int n, const Vec2i* b2a, const Vec3d* apos, int fontTex, float sz=0.02 ){
    for(int i=0; i<n; i++){
        Vec2i ib = b2a[i];
        drawInt( (apos[ib.x]+apos[ib.y])*0.5, i, fontTex, sz );
    }
}

void atomLabels( int n, const Vec3d* apos, int fontTex, float sz=0.02 ){
    for(int i=0; i<n; i++){
        drawInt( apos[i], i, fontTex, sz );
    }
}

void vecsInPos( int n, const Vec3d* vecs, const Vec3d* pos, double sz=1.0 ){
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        Vec3d p=pos[i];         glVertex3f(p.x,p.y,p.z);
        p.add_mul(vecs[i], sz); glVertex3f(p.x,p.y,p.z);
    }
    glEnd();
}

void atomPropertyLabel( int n, double* data, Vec3d* ps, int pitch, int offset, int fontTex, float sz=0.02, const char* format="%4.2f\0" ){
    for(int i=0; i<n; i++){
        drawDouble( ps[i], data[i*pitch+offset], fontTex, sz, format );
        //drawInt( ps[i], (int)data[i*pitch+offset], fontTex, sz );
    }
}

void bondPropertyLabel( int n, double* data, const Vec2i* b2a,  Vec3d* ps, int pitch, int offset, int fontTex, float sz=0.01, const char* format="%4.2f\0" ){
    for(int i=0; i<n; i++){
        Vec2i b = b2a[i];
        Vec3d p = (ps[b.i]+ps[b.j])*0.5;
        drawDouble( p, data[i*pitch+offset], fontTex, sz, format );
        //drawInt( ps[i], (int)data[i*pitch+offset], fontTex, sz );
    }
}

void bonds( int n, const Vec2i* b2a, const Vec3d* apos){
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        Vec2i b = b2a[i];
        //Draw3D::drawLine( apos[b.b], apos[b.a] );
        Draw3D::vertex( apos[b.b] );
        Draw3D::vertex( apos[b.a] );
    }
    glEnd();
}

void bondsLengths( int n, const Vec2i* b2a, const Vec3d* apos, int fontTex, float sz=0.01, const char* format="%4.2f\0" ){
    glBegin(GL_LINES);
    for(int i=0; i<n; i++){
        Vec2i b = b2a[i];
        const Vec3d& pi=apos[b.i];
        const Vec3d& pj=apos[b.j];
        double r = (pi-pj).norm();
        drawDouble( (pi+pj)*0.5, r, fontTex, sz, format );
    }
    glEnd();
}


void pbcBondNeighLabels( int n, const Vec2i* b2a, const Vec3d* apos, const Vec3d* pbcShifts, int fontTex, float sz=0.01 ){
    for(int ib=0; ib<n; ib++){
        if (pbcShifts[ib].norm2()<0.1) continue;
        Vec2i b = b2a[ib];
        //Vec3d pa=apos[b.a]-pbcShifts[ib];
        //Vec3d pb=apos[b.b]+pbcShifts[ib];
        //printf( "Draw pbc_bond[%i](%i,%i) pa(%g,%g,%g) pb(%g,%g,%g) \n", ib, b.a, b.b, pa.x,pa.y,pa.z,   pb.x,pb.y,pb.z );
        drawInt( apos[b.a]-pbcShifts[ib], b.a, fontTex, sz );
        drawInt( apos[b.b]+pbcShifts[ib], b.b, fontTex, sz );
    }
}


void bondsPBC( int n, const Vec2i* b2a, const Vec3d* apos, const Vec3d* pbc_shifts ){
    //printf( "bondsPBC &b2a=%li &apos=%li &pbc_shifts=%li \n", (long)b2a, (long)apos, (long)pbc_shifts );
    for(int i=0; i<n; i++){
        Vec2i b = b2a[i];
        Draw3D::drawLine( apos[b.b], apos[b.a]- pbc_shifts[i] );
        Draw3D::drawLine( apos[b.a], apos[b.b]+ pbc_shifts[i] );
    }
}

void bondsPBC( int n, const Vec2i* b2a, const Vec3d* apos, const Vec3i* pbc, const Mat3d& lvec ){
    for(int i=0; i<n; i++){
        Vec2i b = b2a[i];
        Vec3i G = pbc[i];
        if((G.a!=0)||(G.b!=0)||(G.c!=0)){
        Draw3D::drawLine( apos[b.b], apos[b.a]+ lvec.a*-G.a + lvec.b*-G.b + lvec.c*-G.c );}
        Draw3D::drawLine( apos[b.a], apos[b.b]+ lvec.a*G.a + lvec.b*G.b + lvec.c*G.c );
    }
}

void angle( const Vec3i& ang, const Vec2d& cs0, const Vec3d* apos, int fontTex ){
    Draw3D::drawTriangle( apos[ang.a], apos[ang.b], apos[ang.c], true );
    Draw3D::drawDouble( (apos[ang.a]+apos[ang.c])*0.5, atan2( cs0.y, cs0.x )*2*180/M_PI, fontTex );
}



int drawESP( int na, Vec3d* apos, Vec3d* REQs, Vec3d REQ ){
    const int nsph=300;
    Vec3d sph[nsph];

    std::vector<int> selection;

    // Fibonachi sampling : see here: http://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
    double goldenRatio = (1. + sqrt(5.))/2.;
    for(int i=0; i<nsph; i++){
        double theta = i*(2*M_PI/ goldenRatio);
        double phi   = acos(1 - 2*(i+0.5)/nsph);
        Vec2d cst; cst.fromAngle(theta); 
        Vec2d csp; csp.fromAngle(phi); 
        sph[i]=Vec3d{ 
            cst.x * csp.y,
            cst.y * csp.y,
            csp.x           
        };
        //printf("sph[%i] p(%g,%g,%g) tp(%g,%g) cst(%g,%g) csp(%g,%g)\n", i, sph[i].x,sph[i].y,sph[i].z, theta, phi, cst.x,cst.y, csp.x,csp.y ); 
    }
    int np=0;
    double fIn = 1.0;
    //printf("na %i npsh =%i \n", na, nsph );  
    for(int ia=0; ia<na; ia++){   
        //printf("ia=%i\n", ia );   
        glPointSize(10.0);
        glBegin(GL_POINTS);  
        //glBegin(GL_LINE_STRIP); 
        

        double Rsph = REQs[ia].x + REQ.x;
        Vec3d pi    = apos[ia];
        Vec3d REQi  = REQs[ia];

        /*

        /// NOTE : This optimization would work only for generating the points, we still need to calculate contribution form all atoms !!!!

        // we select only those atoms which can possibly overlap - Optimization for large systems
        selection.clear();
        for(int ja=0; ja<na; ja++){ 
            Vec3d dp=pi-apos[ja];
            double R  = Rsph + (REQ.x +  REQs[ja].x)*fIn;
            double r2 = dp.norm2(); 
            if( r2<(R*R) ){
                selection.push_back( ja );
            }
        }
        */

        // generate points on surface of atom sphere and filter-out those inside other atoms
        for(int i=0; i<nsph; i++){
            Vec3d p = apos[ia] + sph[i]*Rsph;
            bool bValid=true;
            double e=0;
            Vec3d f;
            for(int ja=0; ja<na; ja++){
            //for(int ja : selection ){ 
                Vec3d dp=p-apos[ja];
                Vec3d REQij;
                combineREQ( REQ, REQs[ja], REQij );
                double r2 = dp.norm2();
                if( (ja!=ia) && ( r2 < sq(REQij.x*fIn) ) ){ bValid=false; break; } // exclude point which is inside other atom
                e += addAtomicForceLJQ( dp, f, REQij );
            }
            if(bValid){
                glColor3f(1-fmax(0,-e),1-fmax(e,-e),1-fmax(0,+e));
                glVertex3f(p.x,p.y,p.z);
                np++;
            }
        }
        glEnd();
    }
    //printf("drawESP np=%i \n", np); //exit(0);
    return np;
}






#ifdef MMFFparams_h
void atoms( int n, Vec3d* ps, int* atypes, const MMFFparams& params, int ogl_sph, float qsc=1, float Rsc=1, float Rsub=0 ){
    glEnable(GL_LIGHTING);
    glEnable(GL_DEPTH_TEST);
    glShadeModel(GL_SMOOTH);
    for(int i=0; i<n; i++){
        const AtomType& atyp = params.atypes[atypes[i]];
        Draw::setRGB( atyp.color );
        Draw3D::drawShape( ogl_sph, ps[i], Mat3dIdentity*((atyp.RvdW-Rsub)*Rsc) );
    }
}
#endif

#ifdef MMFFBuilder_h
void drawBonds( const MM::Builder& builder ){
    //drawSystem( false, true, false );
    glBegin( GL_LINES );
    for(int ib=0; ib<builder.bonds.size(); ib++ ){
        const MM::Bond& b = builder.bonds[ib]; 
        //printf( "bond[%i] (%i,%i) \n)", ib, b.atoms.a, b.atoms.b );
        //Draw3D::drawLine( builder.atoms[b.atoms.a].pos, builder.atoms[b.atoms.b].pos );
        Draw3D::vertex(builder.atoms[b.atoms.a].pos);
        Draw3D::vertex(builder.atoms[b.atoms.b].pos);
    }
    glEnd();
}
#endif
#ifdef MMFFBuilder_h
void drawNeighs( const MM::Builder& builder ){
    glBegin( GL_LINES );
    //printf( "DEBUG drawNeighs() \n" );
    for(int ia=0; ia<builder.atoms.size(); ia++ ){
        const MM::Atom& a = builder.atoms[ia];
        if( a.iconf<0 ) continue;
        //printf( "a.iconf %i \n", a.iconf );
        Vec3d pa = builder.atoms[ia].pos;
        MM::AtomConf c = builder.confs[a.iconf];
        for(int j=0; j<N_NEIGH_MAX; j++ ){
            int ib = c.neighs[j];
            //printf( "neigh[%i] = %i \n", j, ib );
            if( ib>=0 ){
                int ja = builder.bonds[ib].getNeighborAtom(ia);
                //Draw3D::drawLine( pa, builder.atoms[c.neighs[j]].pos  );
                Draw3D::vertex(pa);
                Draw3D::vertex(builder.atoms[ja].pos);
            }
        }
    }
    glEnd();
}
#endif

#ifdef MMFFsp3_h
void drawBonds( const MMFFsp3& ff, double Fsc=0.0 ){
    //drawSystem( false, true, false );
    glBegin( GL_LINES );
    for(int ib=0; ib<ff.nbonds; ib++ ){
        const Vec2i& b = ff.bond2atom[ib]; 
        Draw3D::vertex(ff.apos[b.i]);
        Draw3D::vertex(ff.apos[b.j]);
    }
    glEnd();
}
#endif
#ifdef MMFFsp3_h
void drawNeighs( const MMFFsp3& ff, double Fsc=0.0 ){
    //drawSystem( false, true, false );
    for(int ia=0; ia<ff.nnode; ia++ ){
        //printf( "atom[%i]\n", ia );
        int* ngs = ff.aneighs + ia*ff.nneigh_max;
        for(int j=0; j<ff.nneigh_max; j++ ){
            //printf( "atom[%i]neigh[%i]=%i \n", ia, j, ngs[j] );
            if(ngs[j]>=0){
                glColor3f(0.,0.,0.); Draw3D::drawLine( ff.apos[ia], ff.apos[ngs[j]] );
                if(Fsc>0.0){ glColor3f(1.,0.,0.); Draw3D::drawVecInPos( ff.fapos[ia]*Fsc, ff.apos[ia] ); }
            }else{
                int ipi = -ngs[j]-1;
                glColor3f(0.,0.5,0.); Draw3D::drawVecInPos( ff.pipos[ipi], ff.apos[ia] );
                if(Fsc>0.0){ glColor3f(1.,0.5,0.); Draw3D::drawVecInPos( ff.fpipos[ipi]*Fsc, ff.apos[ia]+ff.pipos[ipi] ); }
            }
        }
    }
}
#endif

}; // namespace Draw3D

#endif

