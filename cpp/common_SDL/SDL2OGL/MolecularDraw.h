#ifndef MolecularDraw_h
#define MolecularDraw_h

#include "Draw3D.h"
#include "AtomicConfiguration.h"
#include "GridFF.h"
#include "Renderer.h"
#include "molecular_utils.h"

void colorRB( float f ){ opengl1renderer.color3f( 0.5+f, 0.5, 0.5-f ); }
//void colorRBH( float f, float h ){ opengl1renderer.color3f( 0.5+f, 0.5+h, 0.5-f ); }
Vec3f colorRBH( float f, float h ){ return {0.5+f+h, 0.5+h, 0.5-f+h}; }
//void colorRB( float f ){ opengl1renderer.color3f( 0.5+f, 0.5+f, 0.5+f ); }
void colorBW( float f ){ opengl1renderer.color3f( 0.5-f, 0.5-f, 0.5-f ); }

void printPoses( int n, double * poses ){
    for( int i=0; i<n; i++ ){
        int i8 = i*8;
        //printf( "force[%04i] %g,%g,%g,%g|%g,%g,%g,%g\n",i, opt.force[i8+0], opt.force[i8+1], opt.force[i8+2], opt.force[i8+3],    opt.force[i8+4], opt.force[i8+5], opt.force[i8+6], opt.force[i8+7]  );
        printf( "[%04i] %g,%g,%g,%g | %g,%g,%g,%g \n",i, poses[i8+0], poses[i8+1], poses[i8+2], poses[i8+3],    poses[i8+4], poses[i8+5], poses[i8+6], poses[i8+7]  );
    }
}

void drawMapedPoints( Renderer* r, const FastAtomicMetric& D, int itest ){
    //atomdist.natoms=1;
    //atomdist.pos[0]=cursor3D;
    //atomdist.toCells(atomdist.ruler.step*0.5-0.01);
    Draw3D::drawBBox( D.ruler.pos0, D.ruler.pmax );
    int j=0;
    for(int i=0; i<D.natoms; i++){
        //Draw3D::drawPointCross( renderer, atomdist.pos[i], atomdist.Rcut );
        //Draw3D::drawPointCross( renderer, atomdist.pos[i], 0.1 );
        bool b = ( i == (itest%D.natoms));
        if(b){ Draw3D::drawSphereOctLines( 16, D.Rcut, D.pos[i] ); }
        else { Draw3D::drawPointCross( r, D.pos[i], 0.1 ); }
        //printf("%i %i \n", i, D.atomNs[i] );
        for(int jj=0; jj<D.atomNs[i];jj++){
            if(b){
                int ic = D.atom2cells[j];
                Vec3i ip;  D.ruler.i2ixyz ( ic, ip );
                Vec3d  p = D.ruler.box2pos( ip, {0.0,0.0,0.0} );
                double d = D.ruler.step;
                Draw3D::drawBBox( p, p+Vec3d{d,d,d} );
            }
            j++;
        }
    }
}

void drawNeighs( Renderer* r, const FastAtomicMetric& D, Vec3d pos ){
    Draw3D::drawBBox( D.ruler.pos0, D.ruler.pmax );
    Draw3D::drawSphereOctLines(16,D.Rcut,pos);
    {
        //ip = atomdist.ruler.i cursor3D
        //Vec3i ip;  atomdist.ruler.i2ixyz ( icell, ip );
        Vec3i ip = D.ruler.ipcell( pos );
        Vec3d  p = D.ruler.box2pos( ip, {0.0,0.0,0.0} );
        double d = D.ruler.step;
        Draw3D::drawBBox( p, p+Vec3d{d,d,d} );
    }
    //printf( "DEBUG 2 \n" );
    if( Box::pointIn( pos, D.ruler.pos0, D.ruler.pmax) ){
        int tmpIs[D.natoms];
        int nfound = D.findNeighs( pos, D.Rcut, tmpIs );
        //printf( "DEBUG 3 \n" );
        //printf( "nfound %i \n", nfound );
        for(int i=0; i<nfound; i++){
            Draw3D::drawLine( pos, D.pos[tmpIs[i]] );
        }
    }
    for(int i=0; i<D.natoms; i++){
        //Draw3D::drawPointCross( renderer, atomdist.pos[i], atomdist.Rcut );
        Draw3D::drawPointCross( r, D.pos[i], 0.1 );
    }
}


void drawPPRelaxTrj( int n, double dt, double damp, GridFF& gff, Vec3d pos, Quat4f PRQ ){
    Vec3d vel = Vec3dZero;
    opengl1renderer.begin(GL_LINE_STRIP);
    for(int i=0; i<n; i++){
        //Vec3d f = Vec3dZero;
        Quat4f fe = Quat4fZero;
        gff.addForce( pos, PRQ, fe );
        vel.mul(damp);
        vel.add_mul( (Vec3d)fe.f, dt);
        pos.add_mul( vel        , dt );
        opengl1renderer.vertex3f(pos.x,pos.y,pos.z);
        //printf( " %i (%g,%g,%g) (%g,%g,%g) \n", i, pos.x,pos.y,pos.z,  f.x,f.y,f.z );
    }
    opengl1renderer.end();
    //exit(0);
}

void drawGridForceAlongLine( Renderer* r, int n, GridFF& gff, Vec3d pos0, Vec3d dpos, Quat4f PRQ, double fsc ){
    Vec3d pos = pos0;
	for( int i=0; i<n; i++ ){
        //Vec3d f = Vec3dZero;
        Quat4f fe = Quat4fZero;
        gff.addForce( pos, PRQ, fe);
        //printf( " %i (%g,%g,%g) (%g,%g,%g) \n", i, pos.x,pos.y,pos.z,  f.x,f.y,f.z );
        Draw3D::drawVecInPos( (Vec3d)fe.f *fsc, pos );
        Draw3D::drawPointCross( r, pos, 0.1 );
        pos.add(dpos);
	}
}

void plotSurfPlane( Vec3d normal, double c0, Vec2d d, Vec2i n ){
    Vec3d da,db;
    normal.getSomeOrtho( da,db );
    da.mul( d.a/da.norm() );
    db.mul( d.b/db.norm() );
    //opengl1renderer.color3f(1.0f,0.0f,0.0f); Draw3D::drawVecInPos(normal, {0.0,0.0,0.0} );
    //opengl1renderer.color3f(0.0f,1.0f,0.0f); Draw3D::drawVecInPos(da*10, {0.0,0.0,0.0} );
    //opengl1renderer.color3f(0.0f,0.0f,1.0f); Draw3D::drawVecInPos(db*10, {0.0,0.0,0.0} );
    Draw3D::drawRectGridLines( n*2, (da*-n.a)+(db*-n.b) + normal*c0, da, db );
}

/*
void renderSubstrate( int n, Vec3d * points, GLenum mode ){
    //printf( "iso_points.size() %i \n", iso_points.size() );
    if( mode == GL_POINTS ){
        opengl1renderer.begin(GL_POINTS);
        for(int i=0; i<iso_points.size(); i++){ opengl1renderer.vertex3f( iso_points[i].x, iso_points[i].y, iso_points[i].z      ); }
        opengl1renderer.end();
    }
}
*/

/*
// =========== OLD ?
void renderSubstrate_( const GridShape& grid, Vec3d * FF, double isoval, bool sign ){
    //printf( "iso_points.size() %i \n", iso_points.size() );
    int nxy = grid.n.x * grid.n.y;
    printf("nxy %i \n", nxy );
    Vec3d * pos     = new Vec3d[nxy];
    Vec3d * normals = new Vec3d[nxy];
    //printf( " -- DEBUG 1 \n" );
    //DEBUG
    getIsoSurfZ( grid, isoval, sign, FF, pos, normals );
    //printf( " -- DEBUG 2 \n" );
    //opengl1renderer.enable(GL_LIGHTING);
    //DEBUG
    for ( int ib=1; ib<grid.n.y; ib++ ){
        opengl1renderer.begin(GL_TRIANGLE_STRIP);
        for ( int ia=0; ia<grid.n.x; ia++ ){
            int ip1 = (ib-1)*grid.n.x + ia;
            int ip2 = (ib  )*grid.n.x + ia;
            //printf( "iba (%i,%i) pos (%g,%g,%g)\n", ib,ia, pos[ip1].x,pos[ip1].y,pos[ip1].z );
            //opengl1renderer.color3f(pos[ip1].z*5-2,1.0f,1.0f); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //opengl1renderer.color3f(pos[ip2].z*5-2,1.0f,1.0f); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
            opengl1renderer.color3f(0.7f,0.7f,0.7f); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            opengl1renderer.color3f(0.8f,0.7f,0.7f); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
        }
        opengl1renderer.end();
    }
    DEBUG
    //printf( " -- DEBUG 3 \n" );
    delete [] pos;
    delete [] normals;
    //exit(0);
}
*/

int renderSubstrate_( const GridShape& grid, Quat4f * FF, Quat4f * FFel, double isoval, bool sign, float sclr=1.0 ){
    //printf( "iso_points.size() %i \n", iso_points.size() );   
    //Vec3d * pos     = new Vec3d[grid.n.x * grid.n.y];
    //Vec3d * pos     = new Vec3d[grid.n.x * grid.n.y];
    //Vec3d * normals = new Vec3d[grid.n.x * grid.n.y];
    double * Zs = new double[grid.n.x * grid.n.y];
    //printf( " -- DEBUG 1 \n" );
    //getIsoSurfZ( grid, isoval, sign, FF, pos, normals );
    getIsoSurfZ( grid, isoval, sign, FF, Zs );
    //opengl1renderer.enable(GL_LIGHTING);
    int nvert = 0;
    //opengl1renderer.disable(GL_LIGHTING);
    for ( int ib=1; ib<=grid.n.y; ib++ ){
        opengl1renderer.begin(GL_TRIANGLE_STRIP);
        //opengl1renderer.begin(GL_LINES);
        for ( int ia=0; ia<=grid.n.x; ia++ ){
            int ip1 = ((ib-1)%grid.n.y)*grid.n.x + (ia%grid.n.x);
            int ip2 = ((ib  )%grid.n.y)*grid.n.x + (ia%grid.n.x);
            //printf( "iba (%i,%i) pos (%g,%g,%g)\n", ib,ia, pos[ip1].x,pos[ip1].y,pos[ip1].z );
            //opengl1renderer.color3f(pos[ip1].z*5-2,1.0f,1.0f); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //opengl1renderer.color3f(pos[ip2].z*5-2,1.0f,1.0f); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
            Vec3f gpos; Quat4f fel1,fel2,  f1,f2;

            Vec3d p1 = grid.dCell.a*ia + grid.dCell.b*(ib-1);
            Vec3d p2 = grid.dCell.a*ia + grid.dCell.b*(ib  );
            p1.z=Zs[ip1];
            p2.z=Zs[ip2];

            grid.cartesian2grid( p1, gpos); fel1 = interpolate3DvecWrap( FFel, grid.n, gpos );
            grid.cartesian2grid( p2, gpos); fel2 = interpolate3DvecWrap( FFel, grid.n, gpos );

            Vec3d nr1,nr2; double invr;
            grid.cartesian2grid( p1, gpos); gpos.z-=3.5; f1 = interpolate3DvecWrap( FF, grid.n, gpos );  f1.z*=0.1; f1.normalize(); f1.mul(-1);
            grid.cartesian2grid( p2, gpos); gpos.z-=3.5; f2 = interpolate3DvecWrap( FF, grid.n, gpos );  f2.z*=0.1; f2.normalize(); f2.mul(-1);
            //opengl1renderer.color3f(0.7f,0.7f,0.7f); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //opengl1renderer.color3f(0.8f,0.7f,0.7f); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
            //opengl1renderer.color3f( fel1.x, fel1.y, fel1.z ); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z);
            //opengl1renderer.color3f( fel2.x, fel2.y, fel2.z ); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z);
            //printf( "[%i,%i]fel1.z %g sclr %g \n", ib,ia, fel1.z, sclr );
            //if( ckeckNaN_d(1, 3, (double*)(pos+ip1), "p2" ) || ckeckNaN_d(1, 3, (double*)(pos+ip2), "p1" ) ){ printf("ERROR in renderSubstrate()[%i,%i]: ip(%i,%i) NaNs Found !!! => Exit() \n",ia,ib, ip1,ip2  );  exit(0); };
            //printf( "renderSubstrate[%i,%i] p1(%g,%g,%g) p2(%g,%g,%g) \n", ia,ib, pos[ip1].x,pos[ip1].y,pos[ip1].z, pos[ip2].x,pos[ip2].y,pos[ip2].z );
            
            
            //colorRB( fel1.z*-sclr ); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z); nvert++;
            //colorRB( fel2.z*-sclr ); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z); nvert++;

            //colorRB( fel1.z*-sclr ); opengl1renderer.normal3f(f1.x,f1.y,f1.z); opengl1renderer.vertex3f(p1.x,p1.y,p1.z); nvert++;
            //colorRB( fel2.z*-sclr ); opengl1renderer.normal3f(f2.x,f2.y,f2.z); opengl1renderer.vertex3f(p2.x,p2.y,p2.z); nvert++;

            colorRB( fel1.z*sclr ); opengl1renderer.normal3f(f1.x,f1.y,f1.z); opengl1renderer.vertex3f(p1.x,p1.y,p1.z); nvert++;
            colorRB( fel2.z*sclr ); opengl1renderer.normal3f(f2.x,f2.y,f2.z); opengl1renderer.vertex3f(p2.x,p2.y,p2.z); nvert++;

            //colorRB( fel1.z*-sclr ); 
            //opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z); opengl1renderer.vertex3f(pos[ip1].x+normals[ip1].x*0.1, pos[ip1].y+normals[ip1].y*0.1, pos[ip1].z+normals[ip1].z*0.1); nvert++;
            
            //printf( "[%i,%i]fel1.e %g sclr %g \n", ib,ia, fel1.e, sclr );
            //colorRB( fel1.e*-sclr ); opengl1renderer.normal3f(normals[ip1].x,normals[ip1].y,normals[ip1].z); opengl1renderer.vertex3f(pos[ip1].x,pos[ip1].y,pos[ip1].z); nvert++;
            //colorRB( fel2.e*-sclr ); opengl1renderer.normal3f(normals[ip2].x,normals[ip2].y,normals[ip2].z); opengl1renderer.vertex3f(pos[ip2].x,pos[ip2].y,pos[ip2].z); nvert++;
        }
        opengl1renderer.end();
    }
    //printf( " -- DEBUG 3 \n" );
    delete [] Zs;
    //delete [] pos;
    //delete [] normals;
    //exit(0);
    return nvert;
}

void renderSubstrate_new( GLMesh* outMesh, const GridFF& gff, Vec2d zrange, double isoval, Quat4d PLQ, double sclr, bool bErrNan=false ){
    printf( "renderSubstrate_new() gff.mode=%i @gff.Bspline_PLQ=%li \n", gff.mode, (long)gff.Bspline_PLQ );
    Quat4d PL{PLQ.x,PLQ.y,0.0,0.0};
    Quat4d Q {0.0,0.0,PLQ.y,0.0};
    Vec3i gn = gff.grid.n;
    Mat3d dCell = gff.grid.dCell;
    int nvert = 0;

    for ( int ib=1; ib<=gn.y; ib++ ){
        outMesh->drawMode = GL_TRIANGLE_STRIP;
        for ( int ia=0; ia<=gn.x; ia++ ){

            Vec3d p1a = dCell.a*ia + dCell.b*(ib-1); p1a.z=zrange.x;
            Vec3d p2a = dCell.a*ia + dCell.b*(ib  ); p2a.z=zrange.x;
            const Vec3d p1b=p1a; p1a.z=zrange.y;
            const Vec3d p2b=p2a; p2a.z=zrange.y;

            Vec3d p1 = gff.findIso( isoval, p1a, p1b, PL, 0.02, bErrNan );
            Vec3d p2 = gff.findIso( isoval, p2a, p2b, PL, 0.02, bErrNan );

            if( isnan(p1.z) ){
                if(bErrNan){
                    printf( "renderSubstrate_new() failed to find isovalue %g at p(%g,%g) zrange(%g,%g) => z-scan: \n", isoval, p1a.x,p1a.y, zrange.x, zrange.y );
                    Vec3d fout;
                    int n=100;
                    Vec3d dp=(p1b-p1a)*(1./n);
                    for(int i=0;i<n;i++){
                        Vec3d p = p1a + dp*(i*1.);
                        double e = gff.addAtom( p, PL, fout);
                        printf( "%8.4f %g\n", p.z, e );
                    }
                    exit(0);
                }else{
                    p1.z = zrange.y;
                }
            }
            if( isnan(p2.z) ){ p2.z = zrange.y; }

            Vec3d f1=Vec3dZ,f2=Vec3dZ;
            double el1 = gff.addAtom( p1, Q, f1 );
            double el2 = gff.addAtom( p2, Q, f2 );

            gff.addAtom( p1, PL, f1 );
            gff.addAtom( p2, PL, f2 );

            
            Vec3f color = colorRBH( el1*sclr, sin(p1.z*1.0)*0.1 );
            Vec3f normal = {f1.x,f1.y,f1.z};
            Vec3f pos = {p1.x, p1.y, p1.z};
            outMesh->addVertex(pos, normal, color);

            color = colorRBH( el2*sclr, sin(p2.z*1.0)*0.1 );
            normal = {f2.x,f2.y,f2.z};
            pos = {p2.x, p2.y, p2.z};
            outMesh->addVertex(pos, normal, color);
        }
    }
}


void viewSubstrate( int nx, int ny, int isoOgl, Vec3d a, Vec3d b, Vec3d pos0=Vec3dZero ){
    opengl1renderer.pushMatrix();
    for( int ix = -nx; ix<=nx; ix++ ){
        for( int iy = -ny; iy<=ny; iy++ ){
            Vec3d pos = a*ix + b*iy + pos0;
            opengl1renderer.translatef(pos.x, pos.y, pos.z);
            opengl1renderer.callList(isoOgl);
            opengl1renderer.translatef(-pos.x, -pos.y, -pos.z);
        }
    }
    opengl1renderer.popMatrix();
}

void viewSubstrate( Renderer* r, Vec2i nxs, Vec2i nys, GLMesh* isoOgl, Vec3d a, Vec3d b, Vec3d pos0=Vec3dZero ){
    for( int ix = nxs.x; ix<=nxs.y; ix++ ){
        for( int iy = nys.x; iy<=nys.y; iy++ ){
            Vec3d pos = a*ix + b*iy + pos0;
            r->drawMesh(isoOgl, (Vec3f)pos);
        }
    }
}

#endif
