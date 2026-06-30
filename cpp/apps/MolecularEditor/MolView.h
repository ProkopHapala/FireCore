/// @file MolView.h
/// @brief Lightweight 3D molecule viewer with auto-fit zoom and info HUD.
///
/// MolView displays a single molecule in 3D with interactive rotation (via the
/// host app's camera). Features:
/// - Renders atoms (colored spheres) and bonds (lines) using Draw3D
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
    float textSize3D = 0.02f;

    // ---- GUI checkboxes for view toggles
    GUI gui;
    CheckBoxList* checkBoxes = nullptr;
    bool bGUIInited = false;

    void initGUI( int WIDTH, int HEIGHT ){
        checkBoxes = new CheckBoxList( 4, 4, 120, fontSizeDef*2 );
        checkBoxes->addBox( "[l] Atom Labels",  &bViewAtomLabels  );
        checkBoxes->addBox( "[t] Atom Types",   &bViewAtomTypes   );
        checkBoxes->addBox( "[b] Bond Labels",  &bViewBondLabels  );
        checkBoxes->addBox( "[i] Bond Lengths", &bViewBondLenghts );
        checkBoxes->update();
        gui.addPanel( checkBoxes );
        bGUIInited = true;
    }

    // ---- State
    Molecule* mol = nullptr;  // not owned — borrowed from BrowserView
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

    void setMolecule( Molecule* mol_, const std::string& name_, int index_ ){
        mol      = mol_;
        molName  = name_;
        molIndex = index_;
        bNeedFit = true;
    }

    void autoFitZoom( AppSDL2OGL_3D* app ){
        if(!mol || mol->natoms==0) return;
        Vec3d pmin=Vec3dmax, pmax=Vec3dmin;
        for(int j=0; j<mol->natoms; j++){ pmin.setIfLower(mol->pos[j]); pmax.setIfGreater(mol->pos[j]); }
        Vec3d span = pmax - pmin;
        double maxspan = fmax( span.x, fmax( span.y, span.z ) );
        double zoom_fit = 5.0;
        if(maxspan > 1e-6) zoom_fit = maxspan * 0.6;
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
        glColor3f(0.0f,0.0f,0.0f);
        Draw3D::bonds ( mol->nbonds, mol->bond2atom, mol->pos );

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
            sprintf( str, "Mol: %s  |  atoms: %i  bonds: %i  |  [l] atom labels [t] atom types [b] bond labels [i] bond lengths  |  [Enter] back  [Esc] quit", molName.c_str(), mol->natoms, mol->nbonds );
        }else{
            sprintf( str, "No molecule loaded  |  [Enter] back to browser" );
        }
        Draw2D::drawText( str, 0, {4, HEIGHT-fontSizeDef*2-2}, 0.0, fontTex, fontSizeDef );
        // ---- Checkbox panel (bottom-left) via GUI
        checkBoxes->syncRead();
        gui.draw();
    }

    // Returns true if the click was consumed by the GUI
    bool handleMouse( int mx, int my, const SDL_Event& event ){
        if(!bGUIInited) return false;
        return gui.onEvent( mx, my, event ) != nullptr;
    }
};

#endif
