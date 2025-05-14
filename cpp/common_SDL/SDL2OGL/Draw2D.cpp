#include "GLMesh.h"
#include "MeshLibrary.h"

#include "Draw.h"
#include "Renderer.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Draw2D.h"  // THE HEADER


float Draw2D::z_layer = 0.0f; // should be initialized like this http://stackoverflow.com/questions/19136059/namespace-global-variable-losing-value-c

void Draw2D::drawPointCross( const Vec2f& vec, float d ){
	drawLine( {vec.x-d, vec.y  }, {vec.x+d, vec.y  } );
	drawLine( {vec.x  , vec.y-d}, {vec.x  , vec.y+d} );
};

void Draw2D::drawLine( const Vec2f& p1, const Vec2f& p2 ){
    MeshLibrary::line2D.setUniform3f("uColor", opengl1renderer.color);
    MeshLibrary::line2D.draw2D( {p1.x, p1.y, z_layer}, p2-p1);
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
	MeshLibrary::circle.draw2D( {center.x, center.y, z_layer}, radius );
};

void Draw2D::plot( int n, float dx, double * ys ){
    opengl1renderer.begin(GL_LINE_STRIP);
    for( int i=0; i<n; i++ ){
        opengl1renderer.vertex3f( i*dx, (float)ys[i], z_layer );
    }
    opengl1renderer.end();
};


void Draw2D::plot( int n, double * xs, double * ys ){
    opengl1renderer.begin(GL_LINE_STRIP);
    for( int i=0; i<n; i++ ){
        opengl1renderer.vertex3f( (float)xs[i], (float)ys[i], z_layer );
    }
    opengl1renderer.end();
};

void Draw2D::plot_dots( int n, double * xs, double * ys ){
    opengl1renderer.begin   (GL_POINTS);
	for( int i=0; i<n; i++ ){ opengl1renderer.vertex3f( (float)xs[i], (float)ys[i], z_layer );}
	opengl1renderer.end();
};

void Draw2D::plot_cross( int n, double * xs, double * ys, double sz ){
	for( int i=0; i<n; i++ ){
		float x = (float)xs[i]; float y = (float)ys[i];
		drawLine( {x-sz, y}, {x+sz, y} );
		drawLine( {x, y-sz}, {x, y+sz} );
	}
};

void Draw2D::plot_X( int n, double * xs, double * ys, double sz ){
	for( int i=0; i<n; i++ ){
        float x = (float)xs[i]; float y = (float)ys[i];
        drawLine( {x-sz, y-sz}, {x+sz, y+sz} );
        drawLine( {x+sz, y-sz}, {x-sz, y+sz} );
	}
};

void Draw2D::plot_O( int n, double * xs, double * ys, double sz ){
	for( int i=0; i<n; i++ ){
        drawCircle( {xs[i], ys[i]}, sz, false );
	}
};

void Draw2D::drawGrid( int n, double * ticks, double lmin, double lmax, bool XorY ){
    if( XorY ){ // X-grid
        for( int i=0; i<n; i++ ){ float x = ticks[i]; drawLine( {x, (float)lmin}, {x, (float)lmax}); }
    }else{ // Y-grid
        for( int i=0; i<n; i++ ){ float y = ticks[i]; drawLine({lmin, (float)y}, {lmax, (float) y}); }
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
