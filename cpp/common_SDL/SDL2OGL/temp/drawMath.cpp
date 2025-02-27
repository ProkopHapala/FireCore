
//#include <SDL2/SDL.h>


#include "drawMath.h" // THE HEADER

void drawPoint( const Vec3d& vec ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_POINTS);	          	     
		opengl1renderer.vertex3d( vec.x, vec.y, vec.z );
	opengl1renderer.end();
};

void drawVec( const Vec3d& vec ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);	          	     
		opengl1renderer.vertex3d( 0, 0, 0 ); opengl1renderer.vertex3d( vec.x, vec.y, vec.z );
	opengl1renderer.end();
};

void drawVecInPos( const Vec3d& v, const Vec3d& pos ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);	          	     
		opengl1renderer.vertex3d( pos.x, pos.y, pos.z ); opengl1renderer.vertex3d( pos.x+v.x, pos.y+v.y, pos.z+v.z );
	opengl1renderer.end();
};

void drawLine( const Vec3d& p1, const Vec3d& p2 ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);	          	     
		opengl1renderer.vertex3d( p1.x, p1.y, p1.z ); opengl1renderer.vertex3d( p2.x, p2.y, p2.z );
	opengl1renderer.end();
};

void drawVecInPos( const Vec3f& v, const Vec3f& pos ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);	          	     
		opengl1renderer.vertex3f( pos.x, pos.y, pos.z ); opengl1renderer.vertex3d( pos.x+v.x, pos.y+v.y, pos.z+v.z );
	opengl1renderer.end();
};

void drawTriangle( const Vec3d& p1, const Vec3d& p2, const Vec3d& p3 ){
	Vec3d d1,d2,normal;
	d1.set( p2 - p1 );
	d2.set( p3 - p1 );
	normal.set_cross(d1,d2);
	normal.normalize();
	opengl1renderer.begin   (GL_TRIANGLES);
		opengl1renderer.normal3d( normal.x, normal.y, normal.z );	          	     
		opengl1renderer.vertex3d( p1.x, p1.y, p1.z ); 
		opengl1renderer.vertex3d( p2.x, p2.y, p2.z );
		opengl1renderer.vertex3d( p2.x, p2.y, p2.z );
	opengl1renderer.end();
};

void drawMatInPos( const Mat3d& mat, const Vec3d& pos ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);	          	     
		opengl1renderer.color3f( 1, 0, 0 ); opengl1renderer.vertex3d( pos.x, pos.y, pos.z ); opengl1renderer.vertex3d( pos.x+mat.xx, pos.y+mat.xy, pos.z+mat.xz );
		opengl1renderer.color3f( 0, 1, 0 ); opengl1renderer.vertex3d( pos.x, pos.y, pos.z ); opengl1renderer.vertex3d( pos.x+mat.yx, pos.y+mat.yy, pos.z+mat.yz );
		opengl1renderer.color3f( 0, 0, 1 ); opengl1renderer.vertex3d( pos.x, pos.y, pos.z ); opengl1renderer.vertex3d( pos.x+mat.zx, pos.y+mat.zy, pos.z+mat.zz );
	opengl1renderer.end();
};

void drawShape( const Vec3d& pos, const Mat3d& rot, int shape ){ 
	opengl1renderer.pushMatrix();
	float glMat[16];
	toGLMat( pos, rot, glMat );
	opengl1renderer.multMatrixf( glMat );
	opengl1renderer.callList( shape ); 
	opengl1renderer.popMatrix();
};

int drawConeFan( int n, float r, const Vec3f& base, const Vec3f& tip ){
	int nvert=0;
	Vec3f a,b,c,c_hat; 
	c.set_sub( tip, base ); 
	c_hat.set_mul( c, 1/c.norm() );
	c_hat.getSomeOrtho( a, b );
	a.normalize();
	b.normalize();
/*
	printf( "a %f %f %f |a| %f \n", a.x, a.y, a.z, a.norm() );
	printf( "b %f %f %f |b| %f \n", b.x, b.y, b.z, b.norm() );
	printf( "c %f %f %f |c| %f \n", c.x, c.y, c.z, c.norm() );
	printf( "c_hat %f %f %f |c_hat| %f \n", c_hat.x, c_hat.y, c_hat.z, c_hat.norm() );
	printf( "a.b %f a.c %f b.c %f \n", a.dot(b), a.dot(c), b.dot(c) );

	opengl1renderer.color3f(1,0,0); drawVecInPos( a,     base );
	opengl1renderer.color3f(0,1,0); drawVecInPos( b,     base );
	opengl1renderer.color3f(0,0,1); drawVecInPos( c_hat, base );
	opengl1renderer.color3f(1,0,0); drawVecInPos( a,     tip );
	opengl1renderer.color3f(0,1,0); drawVecInPos( b,     tip );
	opengl1renderer.color3f(0,0,1); drawVecInPos( c_hat, tip );
	opengl1renderer.enable (GL_LIGHTING);
	opengl1renderer.color3f( 0.8f, 0.8f, 0.8f );
*/
	float ca, sa;
	float pa=1, pb=0;
	//sincos_taylor2<float>( 2*M_PI/n, sa, ca );
	float alfa = 2*M_PI/n;
	sa = sin( alfa );
	ca = cos( alfa );

	Vec3f q; q.set(c); q.add_mul( a, -r );
	float pnab =  c_hat.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	opengl1renderer.begin   ( GL_TRIANGLE_FAN );
	//opengl1renderer.begin   ( GL_LINES );
	//opengl1renderer.begin   ( GL_LINE_STRIP );

	opengl1renderer.normal3f( c_hat.x, c_hat.z, c_hat.z );
	//printf( "pn0 %f %f %f \n", c_hat.x, c_hat.z, c_hat.z );
	opengl1renderer.vertex3f( tip.x, tip.y, tip.z ); nvert++;
	for(int i=0; i<=n; i++ ){
		float pa_ = ca*pa - sa*pb;
		      pb  = sa*pa + ca*pb;     
	          pa  = pa_;
		Vec3f p,pn;
		p .set(   pa*a.x +  pb*b.x,   pa*a.y +  pb*b.y,   pa*a.z +  pb*b.z );
		pn.set( pnab*p.x + pnc*c_hat.x, pnab*p.y + pnc*c_hat.y, pnab*p.z + pnc*c_hat.z  );
		//pn.set( pnab*p.x - pnc*c_hat.x, pnab*p.y - pnc*c_hat.y, pnab*p.z - pnc*c_hat.z  );
		//printf( "p %f %f %f   pn %f %f %f |pn| %f \n", p.x, p.y, p.z,   pn.x, pn.y, pn.z, pn.norm() );
		//printf( " |a| %f |b| %f |p| %f |pn| %f |a+b| %f |ab+c| %f \n", a.norm(), b.norm(),  p.norm(), pn.norm(), sqrt(pa*pa+pb*pb), sqrt(pnc*pnc + pnab*pnab ) );
		opengl1renderer.normal3f( pn.x, pn.y, pn.z );
		opengl1renderer.vertex3f( base.x + r*p.x, base.y + r*p.y, base.z + r*p.z ); nvert++;
/*
		opengl1renderer.vertex3f( tip.x, tip.y, tip.z );
		opengl1renderer.vertex3f( base.x + r*p.x, base.y + r*p.y, base.z + r*p.z );
		opengl1renderer.vertex3f( base.x + r*p.x, base.y + r*p.y, base.z + r*p.z );
		opengl1renderer.vertex3f( base.x + r*p.x+pn.x, base.y + r*p.y+pn.y, base.z + r*p.z+pn.z );
*/
	}
	opengl1renderer.end();
	return nvert;
};


int drawCylinderStrip( int n, float r1, float r2, const Vec3f& base, const Vec3f& tip ){
	int nvert=0;

	Vec3f a,b,c,c_hat; 
	c.set_sub( tip, base ); 
	c_hat.set_mul( c, 1/c.norm() );
	c_hat.getSomeOrtho( a, b );
	a.normalize();
	b.normalize();

	float ca, sa;
	float pa=1, pb=0;
	//sincos_taylor2<float>( M_PI/n, ca, sa );
	float alfa = 2*M_PI/n;
	sa = sin( alfa );
	ca = cos( alfa );

	Vec3f q; q.set(c); q.add_mul( a, -(r1-r2) );
	float pnab =  c_hat.dot( q )/q.norm();
	float pnc  =  sqrt( 1 - pnab*pnab );

	opengl1renderer.begin   ( GL_TRIANGLE_STRIP );
	//opengl1renderer.begin   ( GL_LINES );
	for(int i=0; i<=n; i++ ){
		float pa_ = ca*pa - sa*pb;
		      pb  = sa*pa + ca*pb;     
	          pa  = pa_;
		Vec3f p,pn;
		p .set(   pa*a.x +  pb*b.x,   pa*a.y +  pb*b.y,   pa*a.z +  pb*b.z );
		pn.set( pnab*p.x - pnc*c_hat.x, pnab*p.y - pnc*c_hat.y, pnab*p.z - pnc*c_hat.z  );
		//printf( "p %f %f %f   pn %f %f %f |pn| %f \n", p.x, p.y, p.z,   pn.x, pn.y, pn.z, pn.norm() );
		opengl1renderer.normal3f( pn.x, pn.y, pn.z );
		opengl1renderer.vertex3f( base.x + r1*p.x, base.y + r1*p.y, base.z + r1*p.z ); nvert++;
		opengl1renderer.vertex3f( tip .x + r2*p.x, tip .y + r2*p.y, tip .z + r2*p.z ); nvert++;
/*
		opengl1renderer.vertex3f( base.x + r1*p.x,      base.y + r1*p.y,      base.z + r1*p.z      );
		opengl1renderer.vertex3f( tip .x + r2*p.x,      tip .y + r2*p.y,      tip .z + r2*p.z      );
		opengl1renderer.vertex3f( base.x + r1*p.x,      base.y + r1*p.y,      base.z + r1*p.z      );
		opengl1renderer.vertex3f( base.x + r1*p.x+pn.x, base.y + r1*p.y+pn.y, base.z + r1*p.z+pn.z );
		opengl1renderer.vertex3f( tip .x + r2*p.x,      tip .y + r2*p.y,      tip .z + r2*p.z      );
		opengl1renderer.vertex3f( tip .x + r2*p.x+pn.x, tip .y + r2*p.y+pn.y, tip .z + r2*p.z+pn.z );
*/
	}
	opengl1renderer.end();
	return nvert;	
};


int drawSphereTriangle( int n, float r, const Vec3f& pos, const Vec3f& a, const Vec3f& b, const Vec3f& c ){
	int nvert=0;
	float d = 1.0f/n;
	Vec3f da,db;
	da.set_sub( a, c ); da.mul( d );
	db.set_sub( b, c ); db.mul( d );
	for( int ia=0; ia<n; ia++ ){
		Vec3f p0,p; p0.set( c );
		p0.add_mul( da, ia );
		p.set_mul( p0, 1.0f/p0.norm() );
		opengl1renderer.begin   (GL_TRIANGLE_STRIP);
		//opengl1renderer.begin   (GL_LINES);
		//opengl1renderer.color3f( d*ia, 0, 0 );
		opengl1renderer.normal3f( p.x, p.y, p.z );
		opengl1renderer.vertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
		//opengl1renderer.vertex3f( r*p.x+pos.x+p.x, r*p.y+pos.y+p.y, r*p.z+pos.z+p.z );
		for( int ib=0; ib<(n-ia); ib++ ){
			Vec3f p; 
			p.set_add( p0, da );
			p.normalize();
			//opengl1renderer.color3f( 0, 1, 0 );
			opengl1renderer.normal3f( p.x, p.y, p.z );
			opengl1renderer.vertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			//opengl1renderer.vertex3f( r*p.x+pos.x+p.x, r*p.y+pos.y+p.y, r*p.z+pos.z+p.z );
			p.set_add( p0, db );
			p.normalize();
			//opengl1renderer.color3f( 0, 0, 1 );
			opengl1renderer.normal3f( p.x, p.y, p.z );
			opengl1renderer.vertex3f( r*p.x+pos.x, r*p.y+pos.y, r*p.z+pos.z );   nvert++;
			//opengl1renderer.vertex3f( r*p.x+pos.x+p.x, r*p.y+pos.y+p.y, r*p.z+pos.z+p.z );
			p0.add( db );
			//printf(" %f %f %f %f \n", p.x, p.y, p.z, p.norm() );
		}
		opengl1renderer.end();
	}
	return nvert;	
};

int drawSphere_oct( int n, double r_, const Vec3d& pos_ ){	
	int nvert=0;
	Vec3f pos,px,mx,py,my,pz,mz;
	convert( pos_, pos ); float r = (float)r_;
	px.set( 1,0,0); py.set(0, 1,0); pz.set(0,0, 1);
	mx.set(-1,0,0); my.set(0,-1,0); mz.set(0,0,-1);
	nvert += drawSphereTriangle( n, r, pos, mz, mx, my );
	nvert += drawSphereTriangle( n, r, pos, mz, my, px );
	nvert += drawSphereTriangle( n, r, pos, mz, px, py );
	nvert += drawSphereTriangle( n, r, pos, mz, py, mx );
	nvert += drawSphereTriangle( n, r, pos, pz, mx, my );
	nvert += drawSphereTriangle( n, r, pos, pz, my, px );
	nvert += drawSphereTriangle( n, r, pos, pz, px, py );
	nvert += drawSphereTriangle( n, r, pos, pz, py, mx );
	return nvert;
};

void drawLines( int nlinks, int * links, Vec3d * points ){
	int n2 = nlinks<<1;
	for( int i=0; i<n2; i+=2 ){
		drawLine( points[links[i]], points[links[i+1]] );
		//printf ( " %i %i %i %f %f \n", i, links[i], links[i+1], points[links[i]].x, points[links[i+1]].x );
	}
};


// =================
// from drawUtils.h
// =================


void drawBox( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b ){
	opengl1renderer.begin(GL_QUADS);
		opengl1renderer.color3f( r, g, b );		          	     
		opengl1renderer.normal3f(0,0,-1); opengl1renderer.vertex3f( x0, y0, z0 ); opengl1renderer.vertex3f( x1, y0, z0 ); opengl1renderer.vertex3f( x1, y1, z0 ); opengl1renderer.vertex3f( x0, y1, z0 ); 
		opengl1renderer.normal3f(0,-1,0); opengl1renderer.vertex3f( x0, y0, z0 ); opengl1renderer.vertex3f( x1, y0, z0 ); opengl1renderer.vertex3f( x1, y0, z1 ); opengl1renderer.vertex3f( x0, y0, z1 ); 
		opengl1renderer.normal3f(-1,0,0); opengl1renderer.vertex3f( x0, y0, z0 ); opengl1renderer.vertex3f( x0, y1, z0 ); opengl1renderer.vertex3f( x0, y1, z1 ); opengl1renderer.vertex3f( x0, y0, z1 );	
		opengl1renderer.normal3f(0,0,+1); opengl1renderer.vertex3f( x1, y1, z1 ); opengl1renderer.vertex3f( x0, y1, z1 ); opengl1renderer.vertex3f( x0, y0, z1 ); opengl1renderer.vertex3f( x1, y0, z1 ); 
		opengl1renderer.normal3f(0,+1,1); opengl1renderer.vertex3f( x1, y1, z1 ); opengl1renderer.vertex3f( x0, y1, z1 ); opengl1renderer.vertex3f( x0, y1, z0 ); opengl1renderer.vertex3f( x1, y1, z0 ); 
		opengl1renderer.normal3f(+1,0,0); opengl1renderer.vertex3f( x1, y1, z1 ); opengl1renderer.vertex3f( x1, y0, z1 ); opengl1renderer.vertex3f( x1, y0, z0 ); opengl1renderer.vertex3f( x1, y1, z0 );		
	opengl1renderer.end();
};

int makeBoxList( float x0, float x1, float y0, float y1, float z0, float z1, float r, float g, float b  ){
	int ilist=opengl1renderer.genLists(1);
	opengl1renderer.newList( ilist, GL_COMPILE );
		drawBox( x0, x1, y0, y1, z0, z1, r, g, b );
	opengl1renderer.endList();
	return( ilist );
	// don't forget use opengl1renderer.deleteLists( ilist ,1); later 
}

void drawAxis( float sc ){
	opengl1renderer.disable (GL_LIGHTING);
	opengl1renderer.begin   (GL_LINES);	          	     
		opengl1renderer.color3f( 1, 0, 0 ); opengl1renderer.vertex3f( 0, 0, 0 ); opengl1renderer.vertex3f( 1*sc, 0, 0 );
		opengl1renderer.color3f( 0, 1, 0 ); opengl1renderer.vertex3f( 0, 0, 0 ); opengl1renderer.vertex3f( 0, 1*sc, 0 );
		opengl1renderer.color3f( 0, 0, 1 ); opengl1renderer.vertex3f( 0, 0, 0 ); opengl1renderer.vertex3f( 0, 0, 1*sc );	
	opengl1renderer.end();
};















