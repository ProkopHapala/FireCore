#include "GLMesh.h"
#include "MeshLibrary.h"

#include "Draw.h"
#include "Renderer.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Draw2D.h"  // THE HEADER

GLMesh<MPOS> lineStrip;
GLMesh<MPOS> points;


float Draw2D::z_layer = 0.0f; // should be initialized like this http://stackoverflow.com/questions/19136059/namespace-global-variable-losing-value-c

void Draw2D::drawPointCross( const Vec2f& vec, float d ){
	MeshLibrary::cross.setUniform3f("uColor", opengl1renderer.color);
	MeshLibrary::cross.draw2D({vec.x, vec.y, z_layer}, d);
};

void Draw2D::drawLine( const Vec2f& p1, const Vec2f& p2 ){
	MeshLibrary::line2D.setUniform3f("uColor", opengl1renderer.color);
	MeshLibrary::line2D.draw2D({p1.x, p1.y, z_layer}, p2-p1);
};

void Draw2D::drawRectangle( float p1x, float p1y, float p2x, float p2y, Vec3f color, bool filled){ // TODO: create a list of drawn rects and them draw them at once using instancing?
	MeshLibrary::rect.drawMode = filled ? GL_TRIANGLE_FAN : GL_LINE_LOOP;
	MeshLibrary::rect.setUniform3f("uColor", color);

	glDisable(GL_DEPTH_TEST);
	MeshLibrary::rect.draw2D({p1x, p1y, z_layer}, {p2x-p1x, p2y-p1y});
};

void Draw2D::drawRectangle( const Vec2f& p1, const Vec2f& p2, Vec3f color, bool filled ){
	drawRectangle( p1.x, p1.y, p2.x, p2.y, color, filled );
};

void Draw2D::drawCircle( const Vec2f& center, float radius, bool filled ){
	MeshLibrary::circle.drawMode = filled ? GL_TRIANGLE_FAN : GL_LINE_LOOP;
	MeshLibrary::circle.setUniform3f("uColor", opengl1renderer.color);
	MeshLibrary::circle.draw2D({center.x, center.y, z_layer}, radius);
};

void Draw2D::plot( int n, float dx, double * ys ){
	lineStrip.clear();
	for(int i=0; i<n; i++){
		lineStrip.addVertex({i*dx, (float)ys[i], z_layer});
	}
	lineStrip.setUniform3f("uColor", opengl1renderer.color);
	lineStrip.draw();
};

void Draw2D::plot( int n, double * xs, double * ys ){
	lineStrip.clear();
	for(int i=0; i<n; i++){
		lineStrip.addVertex({(float)xs[i], (float)ys[i], z_layer});
	}
	lineStrip.setUniform3f("uColor", opengl1renderer.color);
	lineStrip.draw();
};

void Draw2D::plot_dots( int n, double * xs, double * ys ){
	points.clear();
	for(int i=0; i<n; i++){
		points.addVertex({(float)xs[i], (float)ys[i], z_layer});
	}
	points.setUniform3f("uColor", opengl1renderer.color);
	points.draw();
};

void Draw2D::plot_cross( int n, double * xs, double * ys, double sz ){
	MeshLibrary::cross.setUniform3f("uColor", opengl1renderer.color);
	for(int i=0; i<n; i++){
		MeshLibrary::cross.draw2D({(float)xs[i], (float)ys[i], z_layer}, sz);
	}
};

void Draw2D::plot_X( int n, double * xs, double * ys, double sz ){
	MeshLibrary::xmark.setUniform3f("uColor", opengl1renderer.color);
	for(int i=0; i<n; i++){
		MeshLibrary::xmark.draw2D({(float)xs[i], (float)ys[i], z_layer}, sz);
	}
};

void Draw2D::plot_O( int n, double * xs, double * ys, double sz ){
	MeshLibrary::circle.drawMode = GL_LINE_LOOP;
	MeshLibrary::circle.setUniform3f("uColor", opengl1renderer.color);
	for(int i=0; i<n; i++){
		MeshLibrary::circle.draw2D({(float)xs[i], (float)ys[i], z_layer}, sz);
	}
};

void Draw2D::drawGrid( int n, double * ticks, double lmin, double lmax, bool XorY ){
	MeshLibrary::line2D.setUniform3f("uColor", opengl1renderer.color);
	if(XorY){ // X-grid
		for(int i=0; i<n; i++){
			float x = ticks[i];
			MeshLibrary::line2D.draw2D({x, (float)lmin, z_layer}, {0, (float)(lmax-lmin)});
		}
	}else{ // Y-grid
		for(int i=0; i<n; i++){
			float y = ticks[i];
			MeshLibrary::line2D.draw2D({(float)lmin, y, z_layer}, {(float)(lmax-lmin), 0});
		}
	}
};

// ===== image and sprite-text

void Draw2D::drawText( const char * str, int nchar, Vec2d pos, float angle, int fontTex, float textSize ){
	glDisable(GL_DEPTH_TEST);
	Draw::drawText( str, {pos.x, pos.y, z_layer}, textSize, nchar );
};

void Draw2D::drawText( const char * str, Vec2d pos, Vec2d sz, int fontTex, float textSize ){
	Vec2i block_size = {(int) sz.x/textSize, (int)sz.y/(2*textSize) };
	glDisable(GL_DEPTH_TEST);
	Draw::drawText(str, {pos.x, pos.y, z_layer}, textSize, block_size );
};
