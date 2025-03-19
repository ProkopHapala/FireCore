
#include "GLMesh.h"

#include "Draw.h"
#include "Renderer.h"
#include "Vec3.h"
#include <GLES2/gl2.h>
#include "Draw2D.h"  // THE HEADER

//namespace Draw2D{

float Draw2D::z_layer = 0.0f; // should be initialized like this http://stackoverflow.com/questions/19136059/namespace-global-variable-losing-value-c

/*void Draw2D::drawPoint( const Vec2f& vec ){
	//opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_POINTS);
		opengl1renderer.vertex3f( vec.x, vec.y, z_layer );
	opengl1renderer.end();
};*/

void Draw2D::drawPointCross( const Vec2f& vec, float d ){
	//opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);
		opengl1renderer.vertex3f( vec.x-d, vec.y,   z_layer   );  opengl1renderer.vertex3f( vec.x+d, vec.y,   z_layer );
		opengl1renderer.vertex3f( vec.x,   vec.y-d, z_layer   );  opengl1renderer.vertex3f( vec.x,   vec.y+d, z_layer  );
	opengl1renderer.end();
};

void Draw2D::drawVec( const Vec2f& vec ){
	//opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);
		opengl1renderer.vertex3f( 0, 0, z_layer ); opengl1renderer.vertex3f( vec.x, vec.y, z_layer );
	opengl1renderer.end();
};

void Draw2D::drawVecInPos( const Vec2f& v, const Vec2f& pos ){
	//opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);
		opengl1renderer.vertex3f( pos.x, pos.y, z_layer ); opengl1renderer.vertex3f( pos.x+v.x, pos.y+v.y, z_layer );
	opengl1renderer.end();
};

void Draw2D::drawBody2d( const Vec2f& rot, const Vec2f& pos, float l1, float l2 ){
    //printf( "(%3.3f,%3.3f) (%3.3f,%3.3f)\n",  rot.x, rot.y, pos.x, pos.y );
	//opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);
		opengl1renderer.vertex3f( pos.x, pos.y, z_layer ); opengl1renderer.vertex3f( pos.x+rot.x*l1, pos.y+rot.y*l1, z_layer );
		opengl1renderer.vertex3f( pos.x, pos.y, z_layer ); opengl1renderer.vertex3f( pos.x+rot.y*l2, pos.y-rot.x*l2, z_layer );
	opengl1renderer.end();
};

void Draw2D::drawLine( const Vec2f& p1, const Vec2f& p2 ){
	//opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);
		opengl1renderer.vertex3f( p1.x, p1.y, z_layer ); opengl1renderer.vertex3f( p2.x, p2.y, z_layer );
	opengl1renderer.end();
};

void Draw2D::drawTriangle( const Vec2f& p1, const Vec2f& p2, const Vec2f& p3 ){
	opengl1renderer.begin   (GL_TRIANGLES);
		opengl1renderer.normal3f( 0, 0, 1 );
		opengl1renderer.vertex3f( p1.x, p1.y, z_layer );
		opengl1renderer.vertex3f( p2.x, p2.y, z_layer );
		opengl1renderer.vertex3f( p3.x, p3.y, z_layer );
	opengl1renderer.end();
};

static GLMesh makeRectMesh(){
    GLMesh m;
    m.addVertex( {0, 0, 0} );
    m.addVertex( {0, 1, 0} );
    m.addVertex( {1, 1, 0} );
    m.addVertex( {1, 0, 0} );
    return m;
}
static GLMesh rectMesh = makeRectMesh();
void Draw2D::drawRectangle( float p1x, float p1y, float p2x, float p2y, Vec3f color, bool filled){ // TODO: create a list of drawn rects and them draw them at once using instancing
    rectMesh.drawMode = filled ? GL_QUADS : GL_LINE_LOOP;
    rectMesh.color = color;

    const float WIDTH = 1820; // TODO: make these not constant
    const float HEIGHT = 980;

    p1x = 2*p1x/WIDTH - 1;
    p2x = 2*p2x/WIDTH - 1;

    p1y = 2*p1y/HEIGHT - 1;
    p2y = 2*p2y/HEIGHT - 1;

    Mat4f mvp;
    mvp.array[ 0] = p2x-p1x;  mvp.array[ 4] = 0;        mvp.array[ 8] = 0;  mvp.array[12] = p1x;
    mvp.array[ 1] = 0;        mvp.array[ 5] = p2y-p1y;  mvp.array[ 9] = 0;  mvp.array[13] = p1y;
    mvp.array[ 2] = 0;        mvp.array[ 6] = 0;        mvp.array[10] = 0;  mvp.array[14] = z_layer;
    mvp.array[ 3] = 0;        mvp.array[ 7] = 0;        mvp.array[11] = 0;  mvp.array[15] = 1;

    opengl1renderer.disable(GL_DEPTH_TEST);
    rectMesh.drawMVP(mvp);
};

void Draw2D::drawRectangle( const Vec2f& p1, const Vec2f& p2, Vec3f color, bool filled ){
	drawRectangle( p1.x, p1.y, p2.x, p2.y, color, filled );
};
/*
void Draw2D::drawPoint_d( const Vec2d& vec ){
	Vec2f vec_;    convert( vec, vec_ ); drawPoint     ( vec_ );
};

void Draw2D::drawPointCross_d ( const Vec2d& vec, double d         ){
	Vec2f vec_;    convert( vec, vec_ ); drawPointCross( vec_, (float)d);
};*/

void Draw2D::drawVec_d( const Vec2d& vec ){
	Vec2f vec_;
	convert( vec, vec_ );
	drawVec( vec_ );
};

void Draw2D::drawVecInPos_d( const Vec2d& v,   const Vec2d& pos ){
	Vec2f v_,pos_; convert( v, v_ );     convert( pos, pos_ ); drawVecInPos( v_, pos_);
};

void Draw2D::drawBody2d_d( const Vec2d& rot,   const Vec2d& pos, float l1, float l2 ){
	Vec2f rot_,pos_; convert( rot, rot_ );     convert( pos, pos_ ); drawBody2d( rot_, pos_, l1, l2 );
};

void Draw2D::drawLine_d( const Vec2d& p1,  const Vec2d& p2  ){
	Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawLine( p1_, p2_);
};

void Draw2D::drawRectangle_d( const Vec2d& p1,  const Vec2d& p2, Vec3f color, bool filled ){
	Vec2f p1_,p2_; convert( p1, p1_ );   convert( p2, p2_ );   drawRectangle( p1_, p2_, color, filled );
};

void Draw2D::drawTriangle_d( const Vec2d& p1,  const Vec2d& p2, const Vec2d& p3 ){
	Vec2f p1_,p2_,p3_;  convert( p1, p1_ );  convert( p2, p2_ ); convert( p3, p3_ ); drawTriangle( p1_, p2_, p3_ );
};

void Draw2D::drawCircle( const Vec2f& center, float radius, int n, bool filled ){
    //printf( " z_layer %3.3f \n", z_layer );
	if( filled){ opengl1renderer.begin(GL_TRIANGLE_FAN); }else{ opengl1renderer.begin(GL_LINE_LOOP); };
	float dphi =  6.28318530718f / n;
	Vec2f drot; drot.fromAngle( dphi );
	Vec2f v;    v.set( radius, 0.0f );
	for ( int i=0; i<n; i++ ){
		opengl1renderer.vertex3f( center.x + v.x, center.y + v.y, z_layer );
		v.mul_cmplx( drot );
	}
	opengl1renderer.end();
};


void Draw2D::drawCircle_d( const Vec2d& center_, float radius, int n, bool filled ){
	Vec2f center; convert( center_, center );
	drawCircle( center, radius, n , filled );
};

void Draw2D::drawArc( const Vec2f& center, float radius, float a0, float phi, float dang, bool filled ){
    //printf( " z_layer %3.3f \n", z_layer );
	if( filled){ opengl1renderer.begin(GL_TRIANGLE_FAN); }else{ opengl1renderer.begin(GL_LINE_STRIP); };
	//float phi = a1-a0;
	int n =  floor( phi/dang + 1.0 );
	if(n<0) n=-n;
	dang  = phi/n;
	//Vec2f rot;  rot.fromAngle( a0   );
	//printf( "\n", phi, dang );
	Vec2f drot; drot.fromAngle( dang );
	//drot.mul(1.002); // DEBUG
	Vec2f v;    v   .fromAngle( a0   ); v.mul(radius);
	for ( int i=0; i<n+1; i++ ){
		opengl1renderer.vertex3f( center.x + v.x, center.y + v.y, z_layer );
		v.mul_cmplx( drot );
	}
	opengl1renderer.end();
};

void Draw2D::drawArc_d( const Vec2d& center_, float radius, float a0, float a1, float dang, bool filled ){
	Vec2f center; convert( center_, center );
	drawArc( center, radius,  a0, a1, dang,  filled );
};


/*
void Draw2D::drawRotRect( Vec2d pos, Vec2d rot, Vec2d sz ){
    Vec2d rotT; rotT.set_perp(rot);
    Vec2d p =  pos + rot*(sz.a*-0.5) + rotT*(sz.b*-0.5);
    opengl1renderer.begin(GL_LINE_LOOP);
        opengl1renderer.vertex2d(p.x,p.y); p.add_mul(rot , sz.a);
        opengl1renderer.vertex2d(p.x,p.y); p.add_mul(rotT, sz.b);
        opengl1renderer.vertex2d(p.x,p.y); p.add_mul(rot ,-sz.a);
        opengl1renderer.vertex2d(p.x,p.y);
    opengl1renderer.end();
};

void Draw2D::drawRotT   ( Vec2d pos, Vec2d rot, Vec2d sz ){
    Vec2d rotT; rotT.set_perp(rot);
    Vec2d p;
    opengl1renderer.begin(GL_LINES);
        p = pos;                    opengl1renderer.vertex2d(p.x,p.y); p.add_mul(rot ,sz.a); opengl1renderer.vertex2d(p.x,p.y);
        p = pos + rotT*(-sz.b*0.5); opengl1renderer.vertex2d(p.x,p.y); p.add_mul(rotT,sz.b); opengl1renderer.vertex2d(p.x,p.y);
    opengl1renderer.end();
};

void Draw2D::drawRotTriangle( Vec2d pos, Vec2d rot, Vec2d sz ){
    Vec2d rotT; rotT.set_perp(rot);
    Vec2d p =  pos + rotT*(sz.b*-0.5);
    opengl1renderer.begin(GL_LINE_LOOP);
        p = pos + rotT*(sz.b*-0.5); opengl1renderer.vertex2d(p.x,p.y);
        p = pos + rot *(sz.a     ); opengl1renderer.vertex2d(p.x,p.y);
        p = pos + rotT*(sz.b* 0.5); opengl1renderer.vertex2d(p.x,p.y);
    opengl1renderer.end();
}


void Draw2D::drawPoints( int npoints, Vec2d * points ){
	opengl1renderer.begin   (GL_POINTS);
	for( int i=0; i<npoints; i++ ){
		Vec2f p; convert( points[i], p );
		opengl1renderer.vertex3f( p.x, p.y, z_layer );
	}
	opengl1renderer.end();
};

void Draw2D::drawPoints( int npoints, Vec2d * points, float sc  ){
	opengl1renderer.begin   (GL_LINES);
	for( int i=0; i<npoints; i++ ){
		Vec2f p; convert( points[i], p );
		opengl1renderer.vertex3f( p.x-sc, p.y, z_layer );
		opengl1renderer.vertex3f( p.x+sc, p.y, z_layer );
        opengl1renderer.vertex3f( p.x, p.y-sc, z_layer );
		opengl1renderer.vertex3f( p.x, p.y+sc, z_layer );
	}
	opengl1renderer.end();
};



void Draw2D::drawPlot2D( int np, double* xs, double* ys, Vec2d sc, Vec2d p0 ){ // , Vec2d sc=Vec2dOnes, Vec2d p0=Vec2dZero
	opengl1renderer.begin   (GL_LINE_STRIP);
	for( int i=0; i<np; i++ ){
        //printf( "[%i] %g %g | %g %g \n", i, xs[i], ys[i],  p0.x+xs[i]*sc.x, p0.y+ys[i]*sc.y );
		opengl1renderer.vertex3f( p0.x+xs[i]*sc.x, p0.y+ys[i]*sc.y, z_layer );
	}
	opengl1renderer.end();
};


void Draw2D::drawLines( int n, Vec2d * points ){
	opengl1renderer.begin   (GL_LINE_STRIP);
	for( int i=0; i<n; i++ ){
		Vec2f p1;
		convert( points[i], p1 );
		opengl1renderer.vertex3f( p1.x, p1.y, z_layer );
		//drawLine_d( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
	opengl1renderer.end();
};

void Draw2D::drawLines( int n, Vec2d * points, Vec2d * vecs, float sz ){
	opengl1renderer.begin   (GL_LINES);
	for( int i=0; i<n; i++ ){
		Vec2f p1,p2;
		convert( points[i], p1 );
		convert( vecs  [i], p2 ); p2.mul(sz); p2.add(p1);
		opengl1renderer.vertex3f( p1.x, p1.y, z_layer );
		opengl1renderer.vertex3f( p2.x, p2.y, z_layer );
		//drawLine_d( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
		//printf ( "plot %i (%f,%f)(%f,%f) \n", i, points[i].x,points[i].y,    vecs[i].x,vecs[i].y );
	}
	opengl1renderer.end();
};

void Draw2D::drawLines( int nlinks, int * links, Vec2d * points ){
	int n2 = nlinks<<1;
	opengl1renderer.begin   (GL_LINES);
	for( int i=0; i<n2; i+=2 ){
		Vec2f p1, p2;
		convert( points[links[i]],   p1 );  opengl1renderer.vertex3f( p1.x, p1.y, z_layer );
		convert( points[links[i+1]], p2 );  opengl1renderer.vertex3f( p2.x, p2.y, z_layer );
		//drawLine_d( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
	opengl1renderer.end();
};

void Draw2D::drawConvexPolygon( int n, Vec2d * points, bool filled ){
	if( filled ){opengl1renderer.begin   ( GL_TRIANGLE_FAN );	}else{opengl1renderer.begin   ( GL_LINE_LOOP   );	}
	for( int i=0; i<n; i++ ){
		opengl1renderer.vertex3f( (float)points[i].x, (float)points[i].y, z_layer );
	}
	opengl1renderer.end();
};

void Draw2D::drawPolarFunc( double x0, double y0, double fscale, int n, double phi0, double * data ){
		double dphi = M_PI_2/n;
		opengl1renderer.begin(GL_LINE_STRIP);
		for( int i=-1; i<n; i++ ){
			int ii = i;	if( i<0 ) ii=n-1;
			double cd  = data[ii];
			double phi = dphi*i + phi0;
			double x   = cos( phi );
			double y   = sin( phi );
			opengl1renderer.vertex3f( (float)( x0 + fscale*x ), (float)(y0 + fscale*y), z_layer );
		}
		opengl1renderer.end();
};*/


void Draw2D::plot( int n, float dx, double * ys ){
    opengl1renderer.begin(GL_LINE_STRIP);
    for( int i=0; i<n; i++ ){
        //printf("Draw2D::plot i,x,y %i %f %f\n", i, xs[i], ys[i] );
        opengl1renderer.vertex3f( i*dx, (float)ys[i], z_layer );
    }
    opengl1renderer.end();
    //exit(0);
};


void Draw2D::plot( int n, double * xs, double * ys ){
    opengl1renderer.begin(GL_LINE_STRIP);
    //double DEBUG_sum = 0.0;
    for( int i=0; i<n; i++ ){
        //printf("Draw2D::plot i,x,y %i %f %f\n", i, xs[i], ys[i] );
        //DEBUG_sum += ys[i];
        opengl1renderer.vertex3f( (float)xs[i], (float)ys[i], z_layer );
    }
    //printf( "DEBUG_sum %g \n", DEBUG_sum );
    opengl1renderer.end();
    //exit(0);
};

void Draw2D::plot_dots( int n, double * xs, double * ys ){
    opengl1renderer.begin   (GL_POINTS);
	for( int i=0; i<n; i++ ){ opengl1renderer.vertex3f( (float)xs[i], (float)ys[i], z_layer );}
	opengl1renderer.end();
};

void Draw2D::plot_cross( int n, double * xs, double * ys, double sz ){
    opengl1renderer.begin   (GL_LINES);
	for( int i=0; i<n; i++ ){
        float x = (float)xs[i]; float y = (float)ys[i];
		opengl1renderer.vertex3f( x-sz, y   , z_layer );   opengl1renderer.vertex3f( x+sz, y   , z_layer );
        opengl1renderer.vertex3f( x   , y-sz, z_layer );   opengl1renderer.vertex3f( x   , y+sz, z_layer );
	}
	opengl1renderer.end();
};

void Draw2D::plot_X( int n, double * xs, double * ys, double sz ){
    opengl1renderer.begin   (GL_LINES);
	for( int i=0; i<n; i++ ){
        float x = (float)xs[i]; float y = (float)ys[i];
		opengl1renderer.vertex3f( x-sz, y-sz, z_layer );   opengl1renderer.vertex3f( x+sz, y+sz, z_layer );
        opengl1renderer.vertex3f( x+sz, y-sz, z_layer );   opengl1renderer.vertex3f( x-sz, y+sz, z_layer );
	}
	opengl1renderer.end();
};

void Draw2D::plot_O( int n, double * xs, double * ys, double sz, int ncirc ){
    opengl1renderer.begin   (GL_LINE_LOOP);
    float dphi =  6.28318530718f / ncirc;
    Vec2f drot; drot.fromAngle( dphi );
	for( int i=0; i<n; i++ ){
        float x = (float)xs[i]; float y = (float)ys[i];
        Vec2f v;    v.set( sz, 0.0f );
        for( int j=0; j<ncirc; j++ ){
            opengl1renderer.vertex3f( x + v.x, y + v.y, z_layer );
            v.mul_cmplx( drot );
        }
	}
	opengl1renderer.end();
};
/*
void Draw2D::drawFunc( float xmin, float xmax, int n, Func1d func ){
    opengl1renderer.begin(GL_LINE_STRIP);
    float dx = (xmax-xmin)/n;
    for( float x=xmin; x<=xmax; x+=dx ){
        float y = (float) func( x );
        opengl1renderer.vertex3f( x, y, z_layer );
    }
    opengl1renderer.end();
};

void Draw2D::drawCurve( float tmin, float tmax, int n, Func1d2 func ){
    opengl1renderer.begin(GL_LINE_STRIP);
    float dt = (tmax-tmin)/n;
    for( float t=tmin; t<=tmax; t+=dt ){
        double x,y;
        func( t, x, y );
        opengl1renderer.vertex3f( (float)x, (float)y, z_layer );
    }
    opengl1renderer.end();
};

void Draw2D::drawFuncDeriv( float xmin, float xmax, float d, int n, Func1d func ){
    opengl1renderer.begin(GL_LINE_STRIP);
    float dx = (xmax-xmin)/n;
    for( float x=xmin; x<=xmax; x+=dx ){
        float y  = (float) func( x );
        float y_ = (float) func( x + d );
        float dy = (y_ - y)/d;
        opengl1renderer.vertex3f( x, dy, z_layer );
    }
    opengl1renderer.end();
};*/

void Draw2D::drawGrid( float xmin, float ymin, float xmax, float ymax, float dx, float dy ){
    opengl1renderer.begin(GL_LINES);
    // X-grid
    int nmin,nmax;
    if( xmin>0 ){ nmin=(int)(xmin/dx)+1; }else{ nmin=(int)(xmin/dx);   }
    if( xmax>0 ){ nmax=(int)(xmax/dx);   }else{ nmax=(int)(xmax/dx)-1; }
    for( int i=nmin; i<=nmax; i++ ){
        float x = i*dx;
        opengl1renderer.vertex3f( x, ymin, z_layer );
        opengl1renderer.vertex3f( x, ymax, z_layer );
    }
    // Y-grid
    if( ymin>0 ){ nmin=(int)(ymin/dy)+1; }else{ nmin=(int)(ymin/dy); }
    if( ymax>0 ){ nmax=(int)(ymax/dy);   }else{ nmax=(int)(ymax/dy)-1; }
    for( int i=nmin; i<=nmax; i++ ){
        float y = i*dy;
        opengl1renderer.vertex3f( xmin, y, z_layer );
        opengl1renderer.vertex3f( xmax, y, z_layer );
    }
    opengl1renderer.end();
};

void Draw2D::drawGrid( int n, double * ticks, double lmin, double lmax, bool XorY ){
    opengl1renderer.begin(GL_LINES);
    //printf( "Draw2D::drawGrid() n %i \n", n );
    if( XorY ){ // X-grid
        //for( int i=0; i<=n; i++ ){ float x = ticks[i]; opengl1renderer.vertex3f( x,    (float)lmin, z_layer ); opengl1renderer.vertex3f( x,    (float)lmax, z_layer );  }   // <=n ( this caused segfault when Plot2D::drawAxes()  )
        for( int i=0; i<n; i++ ){ float x = ticks[i]; opengl1renderer.vertex3f( x,    (float)lmin, z_layer ); opengl1renderer.vertex3f( x,    (float)lmax, z_layer );  }
    }else{ // Y-grid
        //for( int i=0; i<=n; i++ ){ float y = ticks[i];   opengl1renderer.vertex3f( lmin, (float)y,    z_layer ); opengl1renderer.vertex3f( lmax, (float)y,    z_layer ); }  // <=n ( this caused segfault when Plot2D::drawAxes()  )
        for( int i=0; i<n; i++ ){ float y = ticks[i];   opengl1renderer.vertex3f( lmin, (float)y,    z_layer ); opengl1renderer.vertex3f( lmax, (float)y,    z_layer ); }
    }
    opengl1renderer.end();
};

/*
void Draw2D::drawTicketAxis( int n, double * ticks, double l0, double lsz, bool XorY ){
    opengl1renderer.begin(GL_LINES);
    // X-grid
    float lplus  = l0+lsz;
    float lminus = l0-lsz;
    doub
    if( XorY ){
        opengl1renderer.vertex3f( ticks[0], (float)l0, z_layer );    opengl1renderer.vertex3f( ticks[n-1],    (float)l0, z_layer );
        for( int i=0; i<=n; i++ ){ float x = ticks[i]; opengl1renderer.vertex3f( x, lminus, z_layer ); opengl1renderer.vertex3f( x, lplus, z_layer ); }
    }else{
        opengl1renderer.vertex3f( (float)l0, ticks[0], z_layer );    opengl1renderer.vertex3f(     (float)l0, ticks[n-1], z_layer );
        for( int i=0; i<=n; i++ ){ float y = ticks[i]; opengl1renderer.vertex3f( lminus, y,    z_layer ); opengl1renderer.vertex3f( lplus, y,    z_layer ); }
    }
    opengl1renderer.end();
};
*/

void Draw2D::drawSimplex( float x, float y, bool s, float step ){
    opengl1renderer.begin   ( GL_TRIANGLES );
    if( s ){
        opengl1renderer.vertex3f( (float)(x+0.5*step ), (float)(y+0.86602540378*step), 0.0f );
        opengl1renderer.vertex3f( (float)(x+1.5*step ), (float)(y+0.86602540378*step), 0.0f );
        opengl1renderer.vertex3f( (float)(x+1.0*step ), (float) y               , 0.0f );
    }else{
        opengl1renderer.vertex3f( (float) x,           (float)y,                 0.0f );
        opengl1renderer.vertex3f( (float)(x+step),     (float)y,                 0.0f );
        opengl1renderer.vertex3f( (float)(x+0.5*step), (float)(y+0.86602540378*step), 0.0f );
    };
    opengl1renderer.end( );
};
/*
void Draw2D::drawSimplexGrid( int n, float step ){
    //opengl1renderer.color3f(0.1f,0.1f,0.1f);
    opengl1renderer.begin( GL_LINES );
    int n2 = 2 * n;
    float stepy = step*0.86602540378f;
    for( int i=-n; i<=n; i++ ){
        opengl1renderer.vertex3f( -n*step,  i*stepy,  0.0f);
        opengl1renderer.vertex3f(  n*step,  i*stepy,  0.0f);
    }
    for( int i=-n/2; i<=n/2; i++ ){
        opengl1renderer.vertex3f( (i-n*0.5f)*step, -n*stepy,  0.0f);
        opengl1renderer.vertex3f( (i+n*0.5f)*step,  n*stepy,  0.0f);

        opengl1renderer.vertex3f( (i+n*0.5f)*step, -n*stepy,  0.0f);
        opengl1renderer.vertex3f( (i-n*0.5f)*step,  n*stepy,  0.0f);
    }

    for( int i=-n/2; i<=n/2; i++ ){
        opengl1renderer.vertex3f( -n*step,        -2*i*stepy,  0.0f);
        opengl1renderer.vertex3f( (i-n*0.5f)*step,   n*stepy,  0.0f);

        opengl1renderer.vertex3f(  n*step,        -2*i*stepy,  0.0f);
        opengl1renderer.vertex3f( (i+n*0.5f)*step,  -n*stepy,  0.0f);

        opengl1renderer.vertex3f(  n*step,         2*i*stepy,  0.0f);
        opengl1renderer.vertex3f( (i+n*0.5f)*step,   n*stepy,  0.0f);

        opengl1renderer.vertex3f( -n*step,         2*i*stepy,  0.0f);
        opengl1renderer.vertex3f( (i-n*0.5f)*step,  -n*stepy,  0.0f);
    }

    opengl1renderer.end();
}*/

void Draw2D::drawShape( const Vec2d& pos, const Vec2d& rot, int shape ){
	opengl1renderer.pushMatrix();
	//opengl1renderer.translatef( pos.x, pos.y , 0 );
	//opengl1renderer.rotatef( phi*(180/M_PI), 0, 0, 1 );
	float glMat[16];
	toGLMat( pos, rot, glMat );

    //Draw::printGLmat( glMat );
	opengl1renderer.multMatrixf( glMat );
	opengl1renderer.callList( shape );
	opengl1renderer.popMatrix();
};


// ===== image and sprite-text

/*
void Draw2D::renderImage( int itex, const Rect2d& rec ){
    opengl1renderer.enable( GL_TEXTURE_2D );
    opengl1renderer.bindTexture( GL_TEXTURE_2D, itex );
    opengl1renderer.color3f(1.0f,1.0f,1.0f);
    //printf( " itex %i \n", itex );
    opengl1renderer.begin(GL_QUADS);
        opengl1renderer.texCoord2f( 0.0f, 1.0f ); opengl1renderer.vertex3f( rec.a.x, rec.a.y, 3.0f );
        opengl1renderer.texCoord2f( 1.0f, 1.0f ); opengl1renderer.vertex3f( rec.b.x, rec.a.y, 3.0f );
        opengl1renderer.texCoord2f( 1.0f, 0.0f ); opengl1renderer.vertex3f( rec.b.x, rec.b.y, 3.0f );
        opengl1renderer.texCoord2f( 0.0f, 0.0f ); opengl1renderer.vertex3f( rec.a.x, rec.b.y, 3.0f );
    opengl1renderer.end();
    opengl1renderer.bindTexture(GL_TEXTURE_2D, 0); // this is not most efficient but safe
};*/

/*
void Draw2D::drawString( const char * str, int imin, int imax, float x, float y, float sz, int itex ){
    const int nchars = 95;
    float persprite = 1.0f/nchars;
    opengl1renderer.enable     ( GL_TEXTURE_2D );
    opengl1renderer.bindTexture( GL_TEXTURE_2D, itex );
    //opengl1renderer.color4f(0.0f,0.5f,0.0f,1.0f);
    opengl1renderer.enable(GL_BLEND);
    opengl1renderer.enable(GL_ALPHA_TEST);
    //opengl1renderer.blendFunc( GL_ONE, GL_ZERO );
    opengl1renderer.blendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    //glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
    //opengl1renderer.blendFunc(GL_ONE_MINUS_DST_ALPHA,GL_DST_ALPHA);
    opengl1renderer.begin(GL_QUADS);
    for(int i=imin; i<imax; i++){
        int isprite = str[i] - 33;
        float offset  = isprite*persprite+(persprite*0.57);
        float xi = i*sz + x;
        opengl1renderer.texCoord2f( offset          , 1.0f ); opengl1renderer.vertex3f( xi,    y,      3.0f );
        opengl1renderer.texCoord2f( offset+persprite, 1.0f ); opengl1renderer.vertex3f( xi+sz, y,      3.0f );
        opengl1renderer.texCoord2f( offset+persprite, 0.0f ); opengl1renderer.vertex3f( xi+sz, y+sz*2, 3.0f );
        opengl1renderer.texCoord2f( offset          , 0.0f ); opengl1renderer.vertex3f( xi,    y+sz*2, 3.0f );
    }
    opengl1renderer.end();
    opengl1renderer.disable  ( GL_BLEND );
    opengl1renderer.disable  ( GL_ALPHA_TEST );
    opengl1renderer.disable  ( GL_TEXTURE_2D );
    opengl1renderer.blendFunc( GL_ONE, GL_ZERO );

}

void Draw2D::drawString( const  char * str, float x, float y, float sz, int itex ){
    const int nchars = 95;
    float persprite = 1.0f/nchars;
    opengl1renderer.enable( GL_TEXTURE_2D );
    opengl1renderer.bindTexture( GL_TEXTURE_2D, itex );
    opengl1renderer.enable(GL_BLEND);
    opengl1renderer.enable(GL_ALPHA_TEST);
    opengl1renderer.blendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    //opengl1renderer.color4f(0.0f,1.0f,0.0f,1.0f);
    opengl1renderer.begin(GL_QUADS);
    for(int i=0; i<65536; i++){
        if( str[i] == 0 ) break; // 0-terminated string
        int isprite = str[i] - 33;
        float offset  = isprite*persprite+(persprite*0.57);
        float xi = i*sz + x;
        opengl1renderer.texCoord2f( offset          , 1.0f ); opengl1renderer.vertex3f( xi,    y,    3.0f );
        opengl1renderer.texCoord2f( offset+persprite, 1.0f ); opengl1renderer.vertex3f( xi+sz, y,    3.0f );
        opengl1renderer.texCoord2f( offset+persprite, 0.0f ); opengl1renderer.vertex3f( xi+sz, y+sz*2, 3.0f );
        opengl1renderer.texCoord2f( offset          , 0.0f ); opengl1renderer.vertex3f( xi,    y+sz*2, 3.0f );
    }
    opengl1renderer.end();
    opengl1renderer.disable( GL_TEXTURE_2D );
    opengl1renderer.disable  ( GL_BLEND );
    opengl1renderer.disable  ( GL_ALPHA_TEST );
    opengl1renderer.disable  ( GL_TEXTURE_2D );
    opengl1renderer.blendFunc( GL_ONE, GL_ZERO );
}
*/
/*
void Draw2D::drawText( const char * str, int nchar, Vec2d pos, int fontTex, float textSize ){
    opengl1renderer.disable    ( GL_LIGHTING   );
    opengl1renderer.disable    ( GL_DEPTH_TEST );
    opengl1renderer.shadeModel ( GL_FLAT       );
    opengl1renderer.pushMatrix();
        opengl1renderer.translatef( pos.x, pos.y, z_layer );
        Draw::drawText( str, fontTex, textSize, 0, nchar );
    opengl1renderer.popMatrix();
};
*/

void Draw2D::drawText( const char * str, int nchar, Vec2d pos, float angle, int fontTex, float textSize ){
    opengl1renderer.disable    ( GL_LIGHTING   );
    opengl1renderer.disable    ( GL_DEPTH_TEST );
    opengl1renderer.shadeModel ( GL_FLAT       );
    opengl1renderer.pushMatrix();
        opengl1renderer.translatef( pos.x, pos.y, z_layer );
        opengl1renderer.rotatef( angle, 0,0,1 );
        Draw::drawText( str, fontTex, textSize, nchar );
    opengl1renderer.popMatrix();
};

void Draw2D::drawText( const char * str, Vec2d pos, Vec2d sz, int fontTex, float textSize ){
    Vec2i block_size = {(int) sz.x/textSize, (int)sz.y/(2*textSize) };
    opengl1renderer.disable    ( GL_LIGHTING   );
    opengl1renderer.disable    ( GL_DEPTH_TEST );
    opengl1renderer.shadeModel ( GL_FLAT       );
    opengl1renderer.pushMatrix();
        opengl1renderer.translatef( pos.x, pos.y, z_layer );
        //Draw::drawText( str, fontTex, textSize, 0, nchar );
        Draw::drawText ( str, fontTex, textSize, block_size );
    opengl1renderer.popMatrix();
};
/*
void Draw2D::drawTextBillboard( const char * str, int nchar, Vec2d pos, int fontTex, float textSize ){
    opengl1renderer.disable    ( GL_LIGHTING   );
    opengl1renderer.disable    ( GL_DEPTH_TEST );
    opengl1renderer.shadeModel ( GL_FLAT       );
    opengl1renderer.pushMatrix();
        opengl1renderer.matrixMode(GL_MODELVIEW);
        // DOES NOT WORK !!!!
        float M[16];
        opengl1renderer.getFloatv(GL_PROJECTION_MATRIX, M );
        printf( "MODEL VIEW: \n" );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[ 0],M[ 1],M[ 2],M[ 3] );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[ 4],M[ 5],M[ 6],M[ 7] );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[ 8],M[ 9],M[10],M[11] );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[12],M[13],M[14],M[15] );
        printf( "PROJECTION: \n" );
        opengl1renderer.getFloatv(GL_PROJECTION_MATRIX, M );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[ 0],M[ 1],M[ 2],M[ 3] );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[ 4],M[ 5],M[ 6],M[ 7] );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[ 8],M[ 9],M[10],M[11] );
        printf( "%3.3g %3.3g %3.3g %3.3g \n", M[12],M[13],M[14],M[15] );
        opengl1renderer.scalef(1/M[5],1/M[10],1);
        opengl1renderer.translatef( pos.x, pos.y, z_layer );
        //Draw::drawText( str, fontTex, textSize, 0, nchar );
        Draw::drawText ( str, fontTex, textSize, nchar );
    opengl1renderer.popMatrix();
};


void Draw2D::drawText( const char * str, const Vec2d& pos, float angle, int fontTex, float textSize, int istart, int iend ){
    opengl1renderer.disable    ( GL_LIGHTING   );
    opengl1renderer.disable    ( GL_DEPTH_TEST );
    opengl1renderer.shadeModel ( GL_FLAT       );
    opengl1renderer.pushMatrix();
        //opengl1renderer.matrixMode(GL_MODELVIEW);
        //opengl1renderer.matrixMode(GL_PROJECTION);
        opengl1renderer.translatef( pos.x, pos.y, z_layer );
        opengl1renderer.rotatef( angle, 0,0,1 );
        //Draw::billboardCam( );
        //Draw::billboardCamProj( );
        //Draw2D::drawString( inputText.c_str(), 0, 0, textSize, fontTex );
        Draw::drawText( str, fontTex, textSize, istart, iend );
    opengl1renderer.popMatrix();
};


void Draw2D::draw_attached_vec( const Vec2d& pos, const Vec2d& rot, const Vec2d& pos0, const Vec2d& rot0, const Vec2d& lengths ){
	Vec2d gpos, grot;
	grot  .set_mul_cmplx( rot0, rot );
	gpos  .set_mul_cmplx( rot0, pos );
	gpos.add( pos0 );
	//printf( " platform.rot %f %f grot %f %f  gpos %f %f \n",  platform.rot.x, platform.rot.y, grot.x, grot.y, gpos.x, gpos.y );
	//drawPointCross( gpos, 0.1 ); drawVecInPos( grot, gpos );

	//float lperp = 0.1;
	//float llong = 0.5;
	//lengths
	double lperp = lengths.x;
	double llong = lengths.y;
	opengl1renderer.begin(GL_LINES);
		opengl1renderer.vertex3f( (float)( gpos.x-grot.x*lperp), (float)(gpos.y-grot.y*lperp), 1 );   opengl1renderer.vertex3f( (float)(gpos.x+grot.x*lperp), (gpos.y+grot.y*lperp), 1 );
		opengl1renderer.vertex3f( (float)( gpos.x-grot.y*llong), (float)(gpos.y+grot.x*llong), 1 );   opengl1renderer.vertex3f( (float)(gpos.x+grot.y*llong), (gpos.y-grot.x*llong), 1 );
	opengl1renderer.end();
};


void Draw2D::drawTriaglePatchBas( Vec2i i0, Vec2i n, int NX, int* basins, double vmin, double vmax ){
    Vec2f a,b,p;
    a.set( 1.0d, 0.0d           ); //a.mul(scale);
    b.set( 0.5d, 0.86602540378d ); //b.mul(scale);
    //opengl1renderer.disable(GL_SMOOTH);
    //int ii = 0;
    double renorm=1.0d/(vmax-vmin);
    for (int iy=0; iy<n.y-1; iy++){
        opengl1renderer.begin( GL_TRIANGLE_STRIP );
        int ii = (i0.y+iy)*NX + i0.x;
        for (int ix=0; ix<n.x; ix++){
            p.set( ix*a.x+iy*b.x, ix*a.y+iy*b.y );
            Draw::color_of_hash( 5454+basins[ii]*14787979 );
            opengl1renderer.vertex3f( p.x    , p.y    , 0 );
            opengl1renderer.vertex3f( p.x+b.x, p.y+b.y, 0 );
            ii++;
        }
        opengl1renderer.end();
    }
}*/

//}; // namespace Draw2D
