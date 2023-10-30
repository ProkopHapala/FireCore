
#ifndef SVG_render_h
#define SVG_render_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"

class SVG_render{ public:

    Vec3d cog = Vec3dZero;
    Mat3d rot = Mat3dIdentity;
    FILE* fout=0;

    //float zoom = 100.0;
    float zoom = 50.0;
    uint32_t color_stroke = 0x000000;
    uint32_t color_fill   = 0x808080;
    float fill_opacity   = 1.0;
    float stroke_opacity = 1.0;
    float stroke_width   = 2.0;
    char* linecap="butt";
    bool bStroke = true;
    bool bFill   = true;
    int font_size = 10;

    Vec2d pmin,pmax;
    Vec2d pcenter;
    

    void open(const char* fname){
        fout = fopen(fname,"w");
        if(fout==0){ printf("ERROR in SVG_render::open(%s):cannot open file\n", fname ); exit(1); }
        printf( "SVG_render::open(%s) pmin(%g,%g) pmax(%g,%g) pc(%g,%g)\n", fname, pmin.x, pmin.y, pmax.x, pmax.y, pcenter.x, pcenter.y );
        //fprintf(fout,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"200\" height=\"200\" viewBox=\"-6 -7 13 13\">\n");
        fprintf(fout,"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%i\" height=\"%i\" >\n", (int)((pmax.x-pmin.x)), (int)((pmax.y-pmin.y))  );
        //fprintf(fout,"<svg xmlns=\"http://www.w3.org/2000/svg\" >\n");
    }
    void close(){ fprintf(fout,"</svg>\n"); fclose(fout); };

    //void update_center(){ pcenter = (pmin+pmax)*0.5; }
    void update_center(){ pcenter = (pmax-pmin)*0.5; }
    void findViewport(int n, Vec3d* ps, float margin=0.5 ){
        pmin=(Vec2d){+1e+37,+1e+37};
        pmax=(Vec2d){-1e+37,-1e+37};
        for(int i=0; i<n; i++){ 
            Vec3d p_; rot.dot_to( ps[i]-cog, p_ );
            Vec2d p__ = p_.xy()*zoom;
            pmin.setIfLower  ( p__ );
            pmax.setIfGreater( p__ );
        }
        margin*=zoom;
        pmin.sub(margin,margin);
        pmax.add(margin,margin);
        update_center();
    }

    // ---- Line
    void drawLine( Vec2d p1, Vec2d p2 ){
        p1.add(pcenter);
        p2.add(pcenter);
        fprintf(fout,"<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" stroke=\"%06x\" stroke-opacity=\"%4.2f\" stroke-width=\"%4.2f\" />\n", p1.x, p1.y, p2.x, p2.y, color_stroke, stroke_opacity, stroke_width );
    }
    void drawLine( const Vec3d& p1, const Vec3d& p2 ){
        Vec3d p1_,p2_;
        rot.dot_to( p1-cog, p1_ );
        rot.dot_to( p2-cog, p2_ );
        drawLine( p1_.xy()*zoom, p2_.xy()*zoom );
    }

    void print(const char* s){ fprintf(fout,"%s",s); };

    // ---- Polyline
    void beginPath(){
        fprintf(fout,"<path d=\"\n");
    }
    void endPath( bool bClose=false, bool bFill=false ){
        if(bClose)fprintf(fout,"Z");
        if(bFill){
            fprintf(fout,"\" stroke=\"#%06x\" stroke-opacity=\"%4.2f\" stroke-width=\"%4.2f\" fill=\"#%06x\" fill-opacity=\"%4.2f\" />\n", color_stroke, stroke_opacity, stroke_width, color_fill, fill_opacity );
        }else{
            fprintf(fout,"\" stroke=\"#%06x\" stroke-opacity=\"%4.2f\" stroke-width=\"%4.2f\" />\n",                                     color_stroke, stroke_opacity, stroke_width );
        }   
    }
    void path_move( Vec2d p ){
        p.add(pcenter);
        fprintf(fout,"M %.0f %.0f \n", p.x, p.y );
    }
    void path_line( Vec2d p ){
        p.add(pcenter);
        fprintf(fout,"L %.0f %.0f \n", p.x, p.y );
    }
    void path_move( const Vec3d& p ){
        Vec3d p_;
        rot.dot_to( p-cog, p_ );
        path_move( p_.xy()*zoom );
    }
    void path_line( const Vec3d& p ){
        Vec3d p_;
        rot.dot_to( p-cog, p_ );
        path_line( p_.xy()*zoom );
    }


    void beginStyle( const char* name ){
        fprintf(fout,"<style>\n .%s {\n", name );
    }
    void writeCurrentStyle(){
        // stroke: #000000;
        // stroke-opacity: 1.00;
        // stroke-width: 1.0;
        // fill: #021440214402144;
        // fill-opacity: 1.00;
        fprintf(fout,"stroke: #%06x ;\n", color_stroke );
        fprintf(fout,"stroke-opacity: %4.1f ;\n", stroke_opacity );
        fprintf(fout,"stroke-width: %4.1f ;\n", stroke_width );
        //fprintf(fout,"fill: \n", color_fill&0xFF, (color_fill>>8)&0xFF, (color_fill>>16)&0xFF, fill_opacity );
        fprintf(fout,"fill-opacity: %4.2f ;\n", fill_opacity );
    }
    void endStyle(){
        fprintf(fout,"}\n</style>\n");
    }
    void writeCurrentStyle(const char* name ){  
        beginStyle( name );
        writeCurrentStyle();
        endStyle();
    }
    // ---- Circle
    void drawCircle( Vec2d p, double r, const char* style=0 ){
        //printf( "SVG_render::drawCircle() p(%g,%g) r=%g \n", p.x, p.y, r );
        r*=zoom;
        p.add(pcenter);
        if(style){
            fprintf(fout,"<circle cx=\"%.0f\" cy=\"%.0f\" r=\"%.0f\" fill=\"#%06x\" class=\"%s\" />\n", p.x, p.y, r, color_fill, style );
        }else{
            if (bStroke){
                if (bFill){
                    fprintf(fout,"<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"#%06x\" stroke-opacity=\"%4.2f\" stroke-width=\"%4.1f\" fill=\"#%06x\" fill-opacity=\"%4.2f\" />\n", p.x, p.y, r, color_stroke, stroke_opacity, stroke_width, color_fill, fill_opacity );
                }else{
                    fprintf(fout,"<circle cx=\"%f\" cy=\"%f\" r=\"%f\" stroke=\"#%06x\" stroke-opacity=\"%4.2f\" stroke-width=\"%4.1f\" />\n", p.x, p.y, r, color_stroke, stroke_opacity, stroke_width );
                }
            }else{
                fprintf(fout,"<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"#%06x\" fill-opacity=\"%4.2f\" />\n", p.x, p.y, r, color_fill, fill_opacity );
            }
        }
    }
    void drawCircle( const Vec3d& p, double r, const char* style=0 ){
        Vec3d p_;
        rot.dot_to( p-cog, p_ );
        drawCircle( p_.xy()*zoom, r, style  );
    }

    // ---- Text
    void drawText( const char* txt, Vec2d p, const char* style=0 ){
        p.add(pcenter);
        if(style){
            fprintf(fout,"<text x=\"%f\" y=\"%f\" class=\"%s\" >%s</text>\n", p.x, p.y, style, txt );
        }else{
            fprintf(fout,"<text x=\"%f\" y=\"%f\" fill=\"#%06x\" fill-opacity=\"%4.2f\" font-size=\"%i\" >%s</text>\n", p.x, p.y, color_fill, font_size, txt );
        }
    }
    void drawText( const char* txt, const Vec3d& p, const char* style=0 ){
        Vec3d p_;
        rot.dot_to( p-cog, p_ );
        drawText( txt, p_.xy()*zoom, style );
    }


};

#endif
