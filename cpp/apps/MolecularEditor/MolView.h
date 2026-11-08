/// @file MolView.h
/// @brief Lightweight 3D molecule viewer with auto-fit zoom and info HUD.
///
/// MolView displays a single molecule in 3D with interactive rotation (via the
/// host app's camera). Features:
/// - Renders atoms (colored spheres) and bonds (lines) using Draw3D
/// - Optional MMFF bond-length color map per bond type (±5% of l0) with 2D legend
/// - Topology AABB wireframe overlay from TopologyNpz
/// - Auto-fit zoom based on molecule bounding box
/// - Info HUD: molecule name, atom count, bond count
/// - Does not own the molecule — takes a pointer from the browser
/// - Does not own the GL window — delegates camera/viewport to host AppSDL2OGL_3D
///
/// Role in repo: Reusable lightweight molecule viewer component used by
/// MolecularBrowser.cpp when in VIEW mode. Lighter than MolGUI (no editing,
/// no force-field, no simulation). Inspired by Vis3D and GLView patterns
/// but specialized for molecular rendering.

#ifndef MolView_h
#define MolView_h

#include <string>
#include <vector>
#include <map>
#include <stdio.h>

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw3D.h"
#include "Draw2D.h"
#include "Draw.h"
#include "Molecule.h"
#include "MMFFparams.h"
#include "Draw3D_Molecular.h"
#include "AppSDL2OGL_3D.h"
#include "GUI.h"
#include "TopologyNpz.h"

class MolView { public:
    // ---- Shared resources (set by host)
    MMFFparams* params = nullptr;
    int fontTex = 0;    // 2D HUD font
    int fontTex3D = 0;  // 3D billboard font for labels
    int ogl_sph = 0;

    // ---- View toggles
    bool bViewAtomLabels  = false;
    bool bViewAtomTypes   = false;
    bool bViewBondLabels  = false;
    bool bViewBondLenghts = false;
    bool bViewBondLengthColor = false;
    bool bViewAabbOverlay = true;
    float textSize3D = 0.02f;
    static constexpr double bondLengthFrac = 0.05; // ±5% symmetric range per MMFF bond type

    struct BondLegendRow { int ityp, jtyp, order; double l0; };
    std::vector<double> bondL0s;
    std::vector<BondLegendRow> bondLegend;
    bool bondLengthCacheValid = false;

    void rebuildBondLengthCache(){
        bondL0s.clear();
        bondLegend.clear();
        bondLengthCacheValid = false;
        if(!mol || !params || mol->nbonds <= 0) return;
        bondL0s.resize(mol->nbonds);
        std::map<uint64_t, BondLegendRow> seen;
        for(int i=0; i<mol->nbonds; i++){
            Vec2i ib = mol->bond2atom[i];
            int order = mol->bondType[i];
            if(order < 1) order = 1;
            double l0, k;
            params->getBondParams(mol->atomType[ib.x], mol->atomType[ib.y], order, l0, k);
            bondL0s[i] = l0;
            int it = mol->atomType[ib.x], jt = mol->atomType[ib.y];
            if(it > jt){ int t=it; it=jt; jt=t; }
            uint64_t key = ((uint64_t)(uint16_t)it<<32) | ((uint64_t)(uint16_t)jt<<16) | (uint16_t)order;
            if(!seen.count(key)) seen[key] = BondLegendRow{it, jt, order, l0};
        }
        for(const auto& kv : seen) bondLegend.push_back(kv.second);
        for(size_t i=0; i<bondLegend.size(); i++)
            for(size_t j=i+1; j<bondLegend.size(); j++)
                if(bondLegend[j].ityp < bondLegend[i].ityp || (bondLegend[j].ityp==bondLegend[i].ityp && bondLegend[j].jtyp < bondLegend[i].jtyp))
                    { BondLegendRow t=bondLegend[i]; bondLegend[i]=bondLegend[j]; bondLegend[j]=t; }
        bondLengthCacheValid = true;
    }

    void drawBondLengthLegend2D( int WIDTH ){
        if(!bViewBondLengthColor || bondLegend.empty()) return;
        const int barW = 120, barH = 8, rowH = fontSizeDef*3 + 4;
        int x0 = WIDTH - barW - 12;
        int y = 8 + (int)checkBoxes->boxes.size() * (fontSizeDef*2) + 8;
        Draw::setRGB( 0xF8F8F8 );
        Draw2D::drawRectangle( x0-4, y-4, x0+barW+4, y + (int)bondLegend.size()*rowH + 4, true );
        for(size_t r=0; r<bondLegend.size(); r++){
            const BondLegendRow& e = bondLegend[r];
            char title[64];
            sprintf(title, "%s-%s", params->atypes[e.ityp].name, params->atypes[e.jtyp].name);
            Draw::setRGB( 0x000000 );
            Draw2D::drawText( title, 0, {(float)x0, (float)(y + r*rowH + barH + 2)}, 0.0, fontTex, fontSizeDef );
            double vmin = e.l0 * (1.0 - bondLengthFrac);
            double vmax = e.l0 * (1.0 + bondLengthFrac);
            for(int ix=0; ix<barW; ix++){
                double t = (barW > 1) ? (double)ix / (barW - 1) : 0.5;
                double l = vmin + t * (vmax - vmin);
                double s = (e.l0 > 1e-8) ? (l - e.l0) / (e.l0 * bondLengthFrac) : 0.0;
                if(s < -1.0) s = -1.0; else if(s > 1.0) s = 1.0;
                double abss = fabs(s);
                float amp = (float)fmin(1.0, abss / 0.12);
                float red = (s > 0.0) ? (float)s * amp : 0.f;
                float blu = (s < 0.0) ? (float)(-s) * amp : 0.f;
                glColor3f(red, 0.f, blu);
                Draw2D::drawRectangle( x0+ix, y + r*rowH, x0+ix+1, y + r*rowH + barH, true );
            }
            char lbl[96];
            sprintf(lbl, "%.3f", vmin);
            Draw::setRGB( 0x0000FF );
            Draw2D::drawText( lbl, 0, {(float)x0, (float)(y + r*rowH - fontSizeDef)}, 0.0, fontTex, fontSizeDef*0.85f );
            sprintf(lbl, "%.3f", e.l0);
            Draw::setRGB( 0x000000 );
            Draw2D::drawText( lbl, 0, {(float)(x0 + barW*0.5f - fontSizeDef*1.2f), (float)(y + r*rowH - fontSizeDef)}, 0.0, fontTex, fontSizeDef*0.85f );
            sprintf(lbl, "%.3f", vmax);
            Draw::setRGB( 0xFF0000 );
            Draw2D::drawText( lbl, 0, {(float)(x0 + barW - fontSizeDef*3.5f), (float)(y + r*rowH - fontSizeDef)}, 0.0, fontTex, fontSizeDef*0.85f );
        }
    }

    // ---- GUI checkboxes for view toggles
    GUI gui;
    CheckBoxList* checkBoxes = nullptr;
    bool bGUIInited = false;

    void initGUI( int WIDTH, int HEIGHT ){
        checkBoxes = new CheckBoxList( 4, 4, 180, fontSizeDef*2 );
        checkBoxes->addBox( "[l] Atom Labels",  &bViewAtomLabels  );
        checkBoxes->addBox( "[t] Atom Types",   &bViewAtomTypes   );
        checkBoxes->addBox( "[b] Bond Labels",  &bViewBondLabels  );
        checkBoxes->addBox( "[i] Bond Lengths", &bViewBondLenghts );
        checkBoxes->addBox( "[c] Bond L Color", &bViewBondLengthColor );
        checkBoxes->update();
        gui.addPanel( checkBoxes );
        bGUIInited = true;
    }

    // ---- State
    Molecule* mol = nullptr;  // not owned — borrowed from BrowserView
    const TopologyBboxOverlay* bboxOverlay = nullptr;
    std::string molName;
    int  molIndex = -1;
    bool bNeedFit  = true;
    float savedZoom = 10.0f;

    void init( MMFFparams* params_, int fontTex_, int fontTex3D_, int ogl_sph_ ){
        params   = params_;
        fontTex  = fontTex_;
        fontTex3D = fontTex3D_;
        ogl_sph  = ogl_sph_;
    }

    void setMolecule( Molecule* mol_, const std::string& name_, int index_, const TopologyBboxOverlay* bbox_=nullptr ){
        mol      = mol_;
        molName  = name_;
        molIndex = index_;
        bboxOverlay = bbox_;
        bNeedFit = true;
        bondLengthCacheValid = false;
    }

    void drawAabbWireframes(){
        if(!bViewAabbOverlay || !bboxOverlay || !bboxOverlay->loaded) return;
        const TopologyBboxOverlay& o = *bboxOverlay;
        glDisable(GL_LIGHTING);
        glLineWidth(1.5f);
        for(int g=0; g<o.n_groups; g++){
            float r,gc,b; topology_wgColor(g, r, gc, b);
            glColor3f(r, gc, b);
            double xmin=o.group_bbox_min[g*3+0], ymin=o.group_bbox_min[g*3+1], zmin=o.group_bbox_min[g*3+2];
            double xmax=o.group_bbox_max[g*3+0], ymax=o.group_bbox_max[g*3+1], zmax=o.group_bbox_max[g*3+2];
            glBegin(GL_LINES);
            glVertex3d(xmin,ymin,zmin); glVertex3d(xmax,ymin,zmin);
            glVertex3d(xmin,ymin,zmin); glVertex3d(xmin,ymax,zmin);
            glVertex3d(xmin,ymin,zmin); glVertex3d(xmin,ymin,zmax);
            glVertex3d(xmax,ymax,zmax); glVertex3d(xmin,ymax,zmax);
            glVertex3d(xmax,ymax,zmax); glVertex3d(xmax,ymin,zmax);
            glVertex3d(xmax,ymax,zmax); glVertex3d(xmax,ymax,zmin);
            glVertex3d(xmin,ymax,zmin); glVertex3d(xmax,ymax,zmin);
            glVertex3d(xmin,ymax,zmin); glVertex3d(xmin,ymax,zmax);
            glVertex3d(xmax,ymin,zmin); glVertex3d(xmax,ymax,zmin);
            glVertex3d(xmax,ymin,zmin); glVertex3d(xmax,ymin,zmax);
            glVertex3d(xmin,ymin,zmax); glVertex3d(xmax,ymin,zmax);
            glVertex3d(xmin,ymin,zmax); glVertex3d(xmin,ymax,zmax);
            glEnd();
        }
        glLineWidth(1.0f);
        glEnable(GL_LIGHTING);
    }

    void autoFitZoom( AppSDL2OGL_3D* app ){
        if(!mol || mol->natoms==0) return;
        Vec3d pmin=Vec3dmax, pmax=Vec3dmin;
        for(int j=0; j<mol->natoms; j++){ pmin.setIfLower(mol->pos[j]); pmax.setIfGreater(mol->pos[j]); }
        Vec3d span = pmax - pmin;
        double maxspan = fmax( span.x, fmax( span.y, span.z ) );
        double zoom_fit = 5.0;
        if(maxspan > 1e-6 && maxspan < 1e300) zoom_fit = maxspan * 0.6;
        app->zoom = zoom_fit;
        printf("MolView::autoFitZoom: '%s' maxspan=%.3f zoom=%.3f\n", molName.c_str(), maxspan, zoom_fit);
    }

    void draw3D( AppSDL2OGL_3D* app ){
        if(!mol) return;
        if(bNeedFit){
            autoFitZoom(app);
            bNeedFit = false;
        }
        app->camera();
        glEnable(GL_LIGHTING);
        glEnable(GL_DEPTH_TEST);
        Draw3D::atoms ( mol->natoms, mol->pos, mol->atomType, *params, ogl_sph, 1.0, 0.5, 1.0 );
        if(bViewBondLengthColor && mol->nbonds > 0){
            if(!bondLengthCacheValid) rebuildBondLengthCache();
            glDisable(GL_LIGHTING);
            glLineWidth(2.5f);
            Draw3D::bondLengthColorMapMMFF( mol->nbonds, mol->bond2atom, mol->pos, bondL0s.data(), bondLengthFrac );
            glLineWidth(1.0f);
        }else{
            glColor3f(0.0f,0.0f,0.0f);
            Draw3D::bonds ( mol->nbonds, mol->bond2atom, mol->pos );
        }
        drawAabbWireframes();

        // ---- Labels (using Draw3D_Molecular.h functions)
        glDisable(GL_LIGHTING);
        if(bViewAtomLabels ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomLabels ( mol->natoms, mol->pos, fontTex3D, textSize3D ); }
        if(bViewAtomTypes  ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::atomTypes  ( mol->natoms, mol->pos, mol->atomType, &params->atypes[0], fontTex3D, textSize3D ); }
        if(bViewBondLabels ){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondLabels ( mol->nbonds, mol->bond2atom, mol->pos, fontTex3D, textSize3D ); }
        if(bViewBondLenghts){ glColor3f(0.0f,0.0f,0.0f); Draw3D::bondsLengths( mol->nbonds, mol->bond2atom, mol->pos, fontTex3D, textSize3D ); }
        glEnable(GL_LIGHTING);
    }

    void drawHUD( AppSDL2OGL_3D* app ){
        int WIDTH = app->WIDTH, HEIGHT = app->HEIGHT;
        if(!bGUIInited) initGUI(WIDTH, HEIGHT);
        glDisable(GL_LIGHTING);
        // ---- Info bar at top
        Draw::setRGB( 0x000000 );
        Draw2D::drawRectangle( 0, HEIGHT-fontSizeDef*2-4, WIDTH, HEIGHT, true );
        Draw::setRGB( 0xFFFFFF );
        char str[512];
        if(mol){
            sprintf( str, "Mol: %s  |  atoms: %i  bonds: %i  |  [l] labels [t] types [b] bond# [i] len [c] bond color [g] AABB %s  |  [Esc] back  [Ctrl+Q] quit", molName.c_str(), mol->natoms, mol->nbonds, bViewAabbOverlay?"on":"off" );
        }else{
            sprintf( str, "No molecule loaded  |  [Enter] back to browser" );
        }
        Draw2D::drawText( str, 0, {4, HEIGHT-fontSizeDef*2-2}, 0.0, fontTex, fontSizeDef );
        // ---- Checkbox panel (bottom-left) via GUI
        checkBoxes->syncRead();
        gui.draw();
        if(bViewBondLengthColor){
            if(!bondLengthCacheValid) rebuildBondLengthCache();
            drawBondLengthLegend2D( WIDTH );
        }
    }

    // Returns true if the click was consumed by the GUI
    bool handleMouse( int mx, int my, const SDL_Event& event ){
        if(!bGUIInited) return false;
        return gui.onEvent( mx, my, event ) != nullptr;
    }
};

#endif
