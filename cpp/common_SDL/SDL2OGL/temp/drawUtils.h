
///////////////////////////////////////
//    Geometry Primitives and Shapes
///////////////////////////////////////

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

