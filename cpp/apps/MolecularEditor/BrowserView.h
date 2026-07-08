/// @file BrowserView.h
/// @brief ACDsee-style thumbnail browser for molecular files (.xyz, .mol, .mol2).
///
/// BrowserView manages directory navigation, molecule loading, and thumbnail
/// rendering for a grid-based file browser. Features:
/// - Directory navigation: enter subdirs, go back to parent (..), absolute paths
/// - Molecule loading from .xyz/.mol/.mol2 files with auto-centering and bond finding
/// - Thumbnail rendering: each molecule rendered to an OpenGL texture via offscreen viewport
/// - Auto-fit zoom per molecule based on bounding box
/// - Scroll-preserving grid selection (VIEW ↔ BROWSE)
/// - Filters: hide hessian/spectrum NPZ, npy companions, dot-directories
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
#include "CrystalNpz.h"
#include "TopologyNpz.h"

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
    std::vector<TopologyBboxOverlay> topologyOverlays;

    bool bNeedRerender = true;
    int  pathBarHeight = fontSizeDef*2 + 4;
    int  dirBarHeight  = fontSizeDef*2 + 8;
    int  labelHeight   = fontSizeDef*2 + 2;  // space for filename under each thumbnail
    std::vector<Quat4d> dirRects;
    std::vector<Quat4d> thumbRects;
    std::vector<int>    thumbTileIndex; // global tile index per visible thumbRects entry

    int nrow = 0, ncol = 0;
    int selectedThumb = 0;
    int scrollRow = 0; // first visible grid row (preserved across VIEW ↔ BROWSE)

    int totalTiles() const { return (int)subDirNames.size() + (int)thumbnails.size(); }
    int totalRows() const { int n=totalTiles(); return ncol>0 ? (n+ncol-1)/ncol : 0; }
    int rowsVisible( int HEIGHT ) const {
        int avail = HEIGHT - pathBarHeight;
        float cellH = thumb_size + labelHeight;
        int rv = (int)(avail / cellH);
        return rv < 1 ? 1 : rv;
    }
    void ensureSelectionVisible( int HEIGHT ){
        if(ncol<=0) return;
        int selRow = selectedThumb / ncol;
        int rv = rowsVisible(HEIGHT), tr = totalRows();
        if(tr <= rv){ scrollRow = 0; return; }
        if(selRow < scrollRow) scrollRow = selRow;
        else if(selRow >= scrollRow + rv) scrollRow = selRow - rv + 1;
        if(scrollRow < 0) scrollRow = 0;
        int maxScroll = tr - rv;
        if(scrollRow > maxScroll) scrollRow = maxScroll;
    }

    bool isDirTile( int i ) const { return i >= 0 && i < (int)subDirNames.size(); }
    int molIndexFromTile( int i ) const { return i - (int)subDirNames.size(); }

    /// Enter on selection: -2 = navigated dir, -1 = invalid, >=0 = molecule index
    int activateSelectedTile(){
        if(selectedThumb < 0 || selectedThumb >= totalTiles()) return -1;
        if(isDirTile(selectedThumb)){
            navigateToDir(subDirNames[selectedThumb]);
            return -2;
        }
        int imol = molIndexFromTile(selectedThumb);
        if(imol < 0 || imol >= (int)molecules.size()) return -1;
        return imol;
    }

    void drawFolderTile( const std::string& name, Vec2d p0, Vec2d p1, bool selected ){
        Draw::setRGB( selected ? 0xC8D8FF : 0xE8E8F0 );
        Draw2D::drawRectangle( p0.x, p0.y, p1.x, p1.y, true );
        Draw::setRGB( 0x404060 );
        Draw2D::drawRectangle( p0.x+8, p1.y-28, p0.x+48, p1.y-8, true );
        Draw::setRGB( 0x000000 );
        std::string label = std::string("[") + name + "]";
        Draw2D::drawText( label.c_str(), 0, {p0.x+4, p0.y+4}, 0.0, fontTex, fontSizeDef );
        if(selected){
            glLineWidth(3.0f);
            Draw::setRGB( 0x0000FF );
            Draw2D::drawRectangle( p0.x-1, p0.y-1, p1.x+1, p1.y+1, false );
            glLineWidth(1.0f);
        }
    }

    // ---- Methods

    BrowserView( const std::string& work_dir_ ) : Browser( work_dir_ ) {}

    void init( MMFFparams* params_, int fontTex_, int ogl_sph_ ){
        params  = params_;
        fontTex = fontTex_;
        ogl_sph = ogl_sph_;
        extensions.insert( {"xyz",1} );
        extensions.insert( {"mol",1} );
        extensions.insert( {"mol2",1} );
        extensions.insert( {"npz",1} );
        extensions.insert( {"npy",1} );
    }

    void clearThumbnails(){
        for(size_t i=0; i<thumbnails.size(); i++){ glDeleteTextures(1, &thumbnails[i]); }
        thumbnails.clear();
        for(size_t i=0; i<molecules.size(); i++){ delete molecules[i]; }
        molecules.clear();
        topologyOverlays.clear();
    }

    static std::string fileExt( const std::string& fname ){
        size_t idot = fname.find_last_of('.');
        if(idot==std::string::npos) return "";
        std::string ext = fname.substr(idot+1);
        for(char& c : ext) if(c>='A'&&c<='Z') c += 'a'-'A';
        return ext;
    }

    int loadMoleculeFile( Molecule* mol, const std::string& path, TopologyBboxOverlay& overlay ){
        std::string ext = fileExt(path);
        overlay.clear();
        if(ext=="npz"){
            NpzFile npz = npz_read_file(path.c_str());
            if(topology_has_bbox_keys(npz)){
                loadTopologyNpzOnce(npz, *mol, params, overlay);
            }else{
                CrystalData c = crystal_from_npz(npz);
                crystal_fill_molecule(*mol, c, params);
            }
            return mol->natoms;
        }
        if(ext=="npy") return molecule_loadNpy(*mol, path.c_str(), params, 0);
        return mol->loadByExt(path, 0);
    }

    void filterNpyCompanionFiles(){
        bool hasPos = false;
        for(int i=0; i<(int)fileNames.size(); i++) if(fileNames[i]=="pos.npy") hasPos = true;
        if(!hasPos) return;
        std::vector<std::string> keep;
        for(int i=0; i<(int)fileNames.size(); i++){
            const std::string& f = fileNames[i];
            if(f=="Z.npy" || f=="bonds_ij.npy") continue;
            keep.push_back(f);
        }
        fileNames = keep;
    }

    void filterNonViewerNpzFiles(){
        std::vector<std::string> keep;
        for(int i=0; i<(int)fileNames.size(); i++){
            const std::string& f = fileNames[i];
            if(fileExt(f) != "npz"){ keep.push_back(f); continue; }
            std::string stem = f;
            size_t dot = stem.rfind('.');
            if(dot != std::string::npos) stem = stem.substr(0, dot);
            for(char& c : stem) if(c>='A'&&c<='Z') c += 'a'-'A';
            if(stem.find("hessian") != std::string::npos) continue;
            if(stem.find("spectrum") != std::string::npos) continue;
            keep.push_back(f);
        }
        fileNames = keep;
    }

    void refreshDir(){
        readDir(work_dir);
        filterNpyCompanionFiles();
        filterNonViewerNpzFiles();
        readMoleculess();
        bNeedRerender = true;
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
        selectedThumb = 0;
        scrollRow = 0;
        refreshDir();
    }

    void readMoleculess( bool bOrientFlat=true ){
        Mat3d rot;
        for(int i=0; i<(int)molecules.size(); i++){ delete molecules[i]; }
        molecules.clear();
        topologyOverlays.clear();
        for(int i=0; i<(int)fileNames.size(); i++){
            printf("==================\n");
            printf("readMoleculess[%i]\n", i);
            Molecule* mol = new Molecule();
            std::string path = work_dir+"/"+fileNames[i];
            mol->bindParams( params );
            TopologyBboxOverlay overlay;
            int ret = loadMoleculeFile(mol, path, overlay);
            if(ret>=0){
                mol->addToPos( mol->getCOG_minmax()*-1.0 );
                Vec3d pmin=Vec3dmax, pmax=Vec3dmin;
                for(int j=0; j<mol->natoms; j++){ pmin.setIfLower(mol->pos[j]); pmax.setIfGreater(mol->pos[j]); }
                Vec3d span = pmax - pmin;
                printf("  bbox: span(%.3f,%.3f,%.3f)\n", span.x, span.y, span.z);
                if(bOrientFlat && mol->natoms >= 3){
                    mol->FindRotation( rot );
                    mol->orient( {0,0,0}, rot.a, rot.b );
                }
                mol->assignREQs( *params );
                if(mol->nbonds==0) mol->findBonds_brute( 0.5, true );
                molecules.push_back(mol);
                topologyOverlays.push_back(overlay);
                printf("  loaded mol[%i] '%s' natom %i nbond %i n_groups %i\n", i, fileNames[i].c_str(), mol->natoms, mol->nbonds, overlay.n_groups);
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
        sprintf( str, "Dir: %s  |  dirs+files: %i", work_dir.c_str(), totalTiles() );
        Draw2D::drawText( str, 0, {4, HEIGHT-pathBarHeight+2}, 0.0, fontTex, fontSizeDef );

        // ---- ACDsee-style grid: folder tiles first, then molecule thumbnails
        thumbRects.clear();
        dirRects.clear();
        thumbTileIndex.clear();
        ncol = (int)( WIDTH / thumb_size );
        if(ncol < 1) ncol = 1;
        ensureSelectionVisible( HEIGHT );
        float cellW = thumb_size;
        float cellH = thumb_size + labelHeight;
        float y = HEIGHT - pathBarHeight;
        int itile = scrollRow * ncol;
        int ntiles = totalTiles();
        int maxRow = rowsVisible( HEIGHT );
        for(int irow=0; irow<maxRow && y>0 && itile<ntiles; irow++){
            for(int ix=0; ix<ncol; ix++){
                if(itile >= ntiles) break;
                float x0 = ix*cellW;
                float y1 = y;
                float y0 = y - thumb_size;
                bool selected = (itile == selectedThumb);
                if(isDirTile(itile)){
                    drawFolderTile( subDirNames[itile], {x0,y0}, {x0+cellW-1,y1}, selected );
                }else{
                    int imol = molIndexFromTile(itile);
                    drawThumbnail( thumbnails[imol], {x0,y0}, {x0+cellW-1,y1}, texture_size );
                    Draw::setRGB( 0x000000 );
                    Draw2D::drawText( fileNames[imol].c_str(), fileNames[imol].length(), {x0, y0-labelHeight}, 0.0, fontTex, fontSizeDef );
                    if(selected){
                        glLineWidth(3.0f);
                        Draw::setRGB( 0x0000FF );
                        Draw2D::drawRectangle( x0-1, y0-1, x0+cellW, y1+1, false );
                        glLineWidth(1.0f);
                    }
                }
                thumbRects.push_back( {x0, y0-labelHeight, x0+cellW-1, y1} );
                thumbTileIndex.push_back( itile );
                itile++;
            }
            y -= cellH;
        }

        if(ntiles==0){
            Draw::setRGB( 0x888888 );
            Draw2D::drawText( "No subdirs or molecule files (.xyz, .mol, .mol2, .npz, .npy)", 0, {10, 10}, 0.0, fontTex, fontSizeDef );
        }
    }

    // Returns: -1 = nothing clicked, -2 = dir navigated, >=0 = molecule index
    int handleClick( int mx, int my, int HEIGHT ){
        int my_flipped = HEIGHT - my;
        for(int i=0; i<(int)thumbRects.size(); i++){
            const Quat4d& r = thumbRects[i];
            if( mx>=r.x && mx<=r.z && my_flipped>=r.y && my_flipped<=r.w ){
                selectedThumb = thumbTileIndex[i];
                if(isDirTile(selectedThumb)){
                    printf("click on dir[%i] '%s'\n", selectedThumb, subDirNames[selectedThumb].c_str());
                    navigateToDir(subDirNames[selectedThumb]);
                    return -2;
                }
                int imol = molIndexFromTile(selectedThumb);
                printf("click on mol[%i]\n", imol);
                return imol;
            }
        }
        return -1;
    }

    void moveCursor( int dx, int dy ){
        int n = totalTiles();
        if(n==0 || ncol<=0) return;
        int col = selectedThumb % ncol;
        int row = selectedThumb / ncol;
        col += dx; row += dy;
        if(col < 0) col = 0; if(col >= ncol) col = ncol-1;
        if(row < 0) row = 0;
        int maxRow = (n - 1) / ncol;
        if(row > maxRow) row = maxRow;
        int newIdx = row * ncol + col;
        if(newIdx >= n) newIdx = n - 1;
        if(newIdx < 0) newIdx = 0;
        selectedThumb = newIdx;
    }

    ~BrowserView(){ clearThumbnails(); }
};

#endif
