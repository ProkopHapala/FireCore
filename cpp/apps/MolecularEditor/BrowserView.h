/// @file BrowserView.h
/// @brief ACDsee-style thumbnail browser for molecular files (.xyz, .mol, .mol2).
///
/// BrowserView manages directory navigation, molecule loading, and thumbnail
/// rendering for a grid-based file browser. Features:
/// - Directory navigation: enter subdirs, go back to parent (..), absolute paths
/// - Molecule loading from .xyz/.mol/.mol2 files with auto-centering and bond finding
/// - Thumbnail rendering: each molecule rendered to an OpenGL texture via offscreen viewport
/// - Auto-fit zoom per molecule based on bounding box
/// - 2D HUD: path bar, clickable subdirectory buttons, thumbnail grid
/// - Click hit-testing for directory navigation and molecule selection
///
/// Role in repo: Reusable browser component used by MolecularBrowser.cpp.
/// Inspired by BrowserSDL (browser_sdl.h) but uses OpenGL for 3D molecule thumbnails
/// instead of SDL surfaces for 2D images. Does not own the GL window — delegates
/// camera/viewport to the host AppSDL2OGL_3D. Inherits from Browser (browser.h)
/// for directory reading and file filtering.

#ifndef BrowserView_h
#define BrowserView_h

#include <string>
#include <vector>
#include <map>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw3D.h"
#include "Draw2D.h"
#include "Draw.h"
#include "SDL_utils.h"
#include "Solids.h"
#include "Molecule.h"
#include "MMFFparams.h"
#include "Draw3D_Molecular.h"
#include "browser.h"
#include "AppSDL2OGL_3D.h"

inline void makeTextureEmpty( GLuint& tx, int sz ){
    glGenTextures(1, &tx);
    glBindTexture(GL_TEXTURE_2D, tx);
    glTexImage2D(GL_TEXTURE_2D, 0, 4, sz, sz, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}

class BrowserView : public Browser { public:
    // ---- Shared resources (set by host)
    MMFFparams*  params  = nullptr;
    int fontTex  = 0;
    int ogl_sph  = 0;
    int texture_size = 256;
    int thumb_size   = 256;

    // ---- State
    std::vector<Molecule*>  molecules;
    std::vector<GLuint>     thumbnails;

    bool bNeedRerender = true;
    int  pathBarHeight = fontSizeDef*2 + 4;
    int  dirBarHeight  = fontSizeDef*2 + 8;
    int  labelHeight   = fontSizeDef*2 + 2;  // space for filename under each thumbnail
    std::vector<Quat4d> dirRects;
    std::vector<Quat4d> thumbRects;

    int nrow = 0, ncol = 0;
    int selectedThumb = 0;  // keyboard cursor position

    // ---- Methods

    BrowserView( const std::string& work_dir_ ) : Browser( work_dir_ ) {}

    void init( MMFFparams* params_, int fontTex_, int ogl_sph_ ){
        params  = params_;
        fontTex = fontTex_;
        ogl_sph = ogl_sph_;
        extensions.insert( {"xyz",1} );
        extensions.insert( {"mol",1} );
        extensions.insert( {"mol2",1} );
    }

    void clearThumbnails(){
        for(size_t i=0; i<thumbnails.size(); i++){ glDeleteTextures(1, &thumbnails[i]); }
        thumbnails.clear();
        for(size_t i=0; i<molecules.size(); i++){ delete molecules[i]; }
        molecules.clear();
    }

    std::string parentDir( const std::string& dir ){
        if(dir=="/" || dir.empty()) return "/";
        size_t pos = dir.find_last_of('/');
        if(pos==std::string::npos) return ".";
        if(pos==0) return "/";
        return dir.substr(0, pos);
    }

    void navigateToDir( const std::string& newDir ){
        if(work_dir=="."){
            char cwd[4096]; if(getcwd(cwd,sizeof(cwd))) work_dir = cwd;
        }
        std::string target;
        if(newDir==".."){
            target = parentDir(work_dir);
        }else if(newDir[0]=='/'){
            target = newDir;
        }else{
            target = work_dir + "/" + newDir;
        }
        printf("BrowserView::navigateToDir: '%s' -> '%s'\n", work_dir.c_str(), target.c_str());
        work_dir = target;
        clearThumbnails();
        thumbRects.clear();
        dirRects.clear();
        readDir( work_dir );
        readMoleculess();
        bNeedRerender = true;
    }

    void readMoleculess( bool bOrientFlat=true ){
        Mat3d rot;
        for(int i=0; i<(int)molecules.size(); i++){ delete molecules[i]; }
        molecules.clear();
        for(int i=0; i<(int)fileNames.size(); i++){
            printf("==================\n");
            printf("readMoleculess[%i]\n", i);
            Molecule* mol = new Molecule();
            std::string path = work_dir+"/"+fileNames[i];
            mol->bindParams( params );
            int ret = mol->loadByExt( path, 0 );
            if(ret>=0){
                mol->addToPos( mol->getCOG_minmax()*-1.0 );
                Vec3d pmin=Vec3dmax, pmax=Vec3dmin;
                for(int j=0; j<mol->natoms; j++){ pmin.setIfLower(mol->pos[j]); pmax.setIfGreater(mol->pos[j]); }
                Vec3d span = pmax - pmin;
                printf("  bbox: span(%.3f,%.3f,%.3f)\n", span.x, span.y, span.z);
                if(bOrientFlat){
                    mol->FindRotation( rot );
                    mol->orient( {0,0,0}, rot.a, rot.b );
                }
                mol->assignREQs( *params );
                if(mol->nbonds==0) mol->findBonds_brute( 0.5, true );
                molecules.push_back(mol);
                printf("  loaded mol[%i] '%s' natom %i nbond %i\n", i, fileNames[i].c_str(), mol->natoms, mol->nbonds);
            }else{
                delete mol;
            }
        }
    }

    void renderThumbnails( AppSDL2OGL_3D* app ){
        glViewport(0,0,texture_size,texture_size);
        float savedZoom = app->zoom;
        float savedAspect = app->ASPECT_RATIO;
        app->ASPECT_RATIO = 1;
        for(int i=0; i<(int)molecules.size(); i++){
            GLuint tx;
            makeTextureEmpty( tx, texture_size );
            thumbnails.push_back(tx);
            Molecule& mol = *molecules[i];
            Vec3d pmin=Vec3dmax, pmax=Vec3dmin;
            for(int j=0; j<mol.natoms; j++){ pmin.setIfLower(mol.pos[j]); pmax.setIfGreater(mol.pos[j]); }
            Vec3d span = pmax - pmin;
            double maxspan = fmax( span.x, fmax( span.y, span.z ) );
            double zoom_fit = 5.0;
            if(maxspan > 1e-6) zoom_fit = maxspan * 0.7;
            app->zoom = zoom_fit;
            app->camera();
            glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            Draw3D::atoms ( mol.natoms, mol.pos, mol.atomType, *params, ogl_sph, 1.0, 0.5, 1.0 );
            glColor3f(0.0f,0.0f,0.0f);
            Draw3D::bonds ( mol.nbonds, mol.bond2atom, mol.pos);
            glBindTexture(GL_TEXTURE_2D, tx);
            glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, texture_size, texture_size, 0);
        }
        app->zoom = savedZoom;
        app->ASPECT_RATIO = savedAspect;
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glViewport(0, 0, app->WIDTH, app->HEIGHT);
    }

    void drawThumbnail( int itex, Vec2d p0, Vec2d p1, float sz ){
        glColor4f(1.0f,1.0f,1.0f,1.0f);
        glEnable     ( GL_TEXTURE_2D );
        glBindTexture( GL_TEXTURE_2D, itex );
        glBegin(GL_QUADS);
        glTexCoord2f( 0, 0 ); glVertex3f( p0.x, p1.y, 0.0f );
        glTexCoord2f( 1, 0 ); glVertex3f( p1.x, p1.y, 0.0f );
        glTexCoord2f( 1, 1 ); glVertex3f( p1.x, p0.y, 0.0f );
        glTexCoord2f( 0, 1 ); glVertex3f( p0.x, p0.y, 0.0f );
        glEnd();
        glDisable  ( GL_TEXTURE_2D );
    }

    void drawHUD( AppSDL2OGL_3D* app ){
        int WIDTH = app->WIDTH, HEIGHT = app->HEIGHT;
        glDisable(GL_LIGHTING);

        // ---- Path bar
        Draw::setRGB( 0x000000 );
        Draw2D::drawRectangle( 0, HEIGHT-pathBarHeight, WIDTH, HEIGHT, true );
        Draw::setRGB( 0xFFFFFF );
        char str[512];
        sprintf( str, "Dir: %s", work_dir.c_str() );
        Draw2D::drawText( str, 0, {4, HEIGHT-pathBarHeight+2}, 0.0, fontTex, fontSizeDef );

        // ---- Subdirectory buttons
        dirRects.clear();
        float dirY0 = HEIGHT - pathBarHeight - dirBarHeight;
        float dirY1 = HEIGHT - pathBarHeight;
        float dirX  = 4;
        for(int i=0; i<(int)subDirNames.size(); i++){
            const std::string& name = subDirNames[i];
            float w = name.length()*fontSizeDef + 8;
            if(i==0){ Draw::setRGB( 0xCCCCDD ); }else{ Draw::setRGB( 0xDDDDDD ); }
            Draw2D::drawRectangle( dirX, dirY0, dirX+w, dirY1, true );
            Draw::setRGB( 0x000000 );
            Draw2D::drawText( name.c_str(), name.length(), {dirX+4, dirY0+2}, 0.0, fontTex, fontSizeDef );
            dirRects.push_back( {dirX, dirY0, dirX+w, dirY1} );
            dirX += w + 4;
            if(dirX > WIDTH - 50){ dirY0 -= dirBarHeight; dirY1 -= dirBarHeight; dirX = 4; }
        }

        // ---- Thumbnails (cell = thumb image + label below)
        thumbRects.clear();
        ncol = (int)( WIDTH / thumb_size );
        float cellW = thumb_size;
        float cellH = thumb_size + labelHeight;
        float y = dirY0;  // top of current row (Y increases upward)
        int imol = 0;
        while(y > 0 && imol < (int)thumbnails.size()){
            for(int ix=0; ix<ncol; ix++){
                if(imol >= (int)thumbnails.size()) break;
                float x0 = ix*cellW;
                float y1 = y;              // top
                float y0 = y - thumb_size; // bottom of image
                // draw thumbnail image
                drawThumbnail( thumbnails[imol], {x0,y0}, {x0+cellW-1,y1}, texture_size );
                // draw filename below thumbnail
                Draw::setRGB( 0x000000 );
                Draw2D::drawText( fileNames[imol].c_str(), fileNames[imol].length(), {x0, y0-labelHeight}, 0.0, fontTex, fontSizeDef );
                // selection cursor (green border)
                if(imol == selectedThumb){
                    glLineWidth(3.0f);
                    Draw::setRGB( 0x0000FF );
                    Draw2D::drawRectangle( x0-1, y0-1, x0+cellW, y1+1, false );
                    glLineWidth(1.0f);
                }
                thumbRects.push_back( {x0, y0-labelHeight, x0+cellW-1, y1} );
                imol++;
            }
            y -= cellH;
        }

        if(thumbnails.empty()){
            Draw::setRGB( 0x888888 );
            Draw2D::drawText( "No molecule files (.xyz, .mol, .mol2) in this directory", 0, {10, 10}, 0.0, fontTex, fontSizeDef );
        }
    }

    // Returns: -1 = nothing clicked, -2 = parent dir clicked, >=0 = molecule index clicked
    int handleClick( int mx, int my, int HEIGHT ){
        int my_flipped = HEIGHT - my;
        for(int i=0; i<(int)dirRects.size(); i++){
            const Quat4d& r = dirRects[i];
            if( mx>=r.x && mx<=r.z && my_flipped>=r.y && my_flipped<=r.w ){
                printf("click on dir[%i] '%s'\n", i, subDirNames[i].c_str());
                navigateToDir( subDirNames[i] );
                return -2;
            }
        }
        for(int i=0; i<(int)thumbRects.size(); i++){
            const Quat4d& r = thumbRects[i];
            if( mx>=r.x && mx<=r.z && my_flipped>=r.y && my_flipped<=r.w ){
                printf("click on mol[%i]\n", i);
                selectedThumb = i;
                return i;
            }
        }
        return -1;
    }

    // Keyboard navigation for thumbnail selection
    void moveCursor( int dx, int dy ){
        int n = (int)thumbnails.size();
        if(n==0) return;
        int col = selectedThumb % ncol;
        int row = selectedThumb / ncol;
        col += dx; row += dy;
        if(col < 0) col = 0; if(col >= ncol) col = ncol-1;
        int newRow = row;
        if(newRow < 0) newRow = 0;
        int maxRow = (n - 1) / ncol;
        if(newRow > maxRow) newRow = maxRow;
        int newIdx = newRow * ncol + col;
        if(newIdx >= n) newIdx = n - 1;
        if(newIdx < 0) newIdx = 0;
        selectedThumb = newIdx;
    }

    ~BrowserView(){ clearThumbnails(); }
};

#endif
