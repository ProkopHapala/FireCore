#ifndef SplineGUI_h
#define SplineGUI_h

#include <stdio.h>
#include <math.h>

#include <SDL2/SDL.h>

#include "Vec3.h"
#include "Draw3D.h"
#include "Manipulation.h"


inline void drawSplineSegment(const Vec3d& p0, const Vec3d& p1, const Vec3d& p2, const Vec3d& p3, int nsub){
    glBegin(GL_LINE_STRIP);
    for(int i=0; i<=nsub; i++){
        double u = (double)i / nsub;
        const Quat4d b = Bspline::basis(u);
        Vec3d p = Vec3dZero;
        p.add_mul(p0, b.x);
        p.add_mul(p1, b.y);
        p.add_mul(p2, b.z);
        p.add_mul(p3, b.w);
        glVertex3f(p.x,p.y,p.z);
    }
    glEnd();
}

struct SplineGUI{

    SplineConstraintManager* man = 0;

    bool bEdit    = false;   // enable keyboard control
    bool bEditCps = false;   // edit spline control points

    int  iCp      = 0;

    double step   = 0.1;
    double dt     = 0.01;

    void bind(SplineConstraintManager* m){ man=m; }

    inline void draw(){
        if(!man) return;
        if(man->cps.size()==0) return;

        // draw cps
        glDisable(GL_LIGHTING);
        glLineWidth(2.0);

        // draw smooth spline curve
        glColor3f(0.2f,0.8f,0.2f);
        const int nsub = (int)(1.0/dt);  // subdivisions matching +/- step
        const int n = (int)man->cps.size();
        if(n>=4){
            for(int i=0; i<n-3; i++){
                drawSplineSegment(man->cps[i], man->cps[i+1], man->cps[i+2], man->cps[i+3], nsub);
            }
        }else if(n>=2){
            // fallback to polyline for insufficient points
            glBegin(GL_LINE_STRIP);
            for(const Vec3d& p: man->cps){ glVertex3f(p.x,p.y,p.z); }
            glEnd();
        }

        // draw control points
        glLineWidth(1.0);
        for(int i=0;i<(int)man->cps.size();i++){
            const Vec3d& p = man->cps[i];
            if(i==iCp && bEditCps){ glColor3f(1.0f,0.5f,0.0f); }else{ glColor3f(0.0f,0.5f,1.0f); }
            Draw3D::drawPointCross(p, 0.2);
        }

        // draw anchor target
        glColor3f(1.0f,0.0f,0.0f);
        Draw3D::drawPointCross(man->pos, 0.3);

        glEnable(GL_LIGHTING);
    }

    inline void onKeyDown(SDL_Keycode sym){
        printf("SplineGUI.onKeyDown(): sym=%s @man=%p bEdit=%i\n", SDL_GetKeyName(sym), man, bEdit);
        if(!man) return;
        if(!bEdit) return;

        Vec3d d=Vec3dZero;

        // ---- numpad moves
        bool valid=0;
        switch(sym){
            case SDLK_KP_4: d.x=-step; valid=1; break;
            case SDLK_KP_6: d.x=+step; valid=1; break;
            case SDLK_KP_2: d.y=-step; valid=1; break;
            case SDLK_KP_8: d.y=+step; valid=1; break;
            case SDLK_KP_1: d.z=-step; valid=1; break;
            case SDLK_KP_3: d.z=+step; valid=1; break;

            case SDLK_KP_PLUS : man->moveT(+dt); printf("SplineGUI: t=%g pos(%g,%g,%g)\n", man->t, man->pos.x,man->pos.y,man->pos.z ); valid=1; break;
            case SDLK_KP_MINUS: man->moveT(-dt); printf("SplineGUI: t=%g pos(%g,%g,%g)\n", man->t, man->pos.x,man->pos.y,man->pos.z ); valid=1; break;

            case SDLK_KP_5: bEditCps=!bEditCps; printf("SplineGUI: bEditCps=%i\n", (int)bEditCps ); valid=1; break;

            case SDLK_KP_ENTER: 
                man->bUseSpline=true; 
                man->updatePos(); 
                printf("SplineGUI: re-attach to spline t=%g pos(%g,%g,%g)\n", man->t, man->pos.x,man->pos.y,man->pos.z ); 
                valid=1; 
                break;

            default: break;
        }

        if(!valid) return;

        printf("SplineGUI.onKeyDown(): valid=%i sym=%d\n", valid, (int)sym); 

        if(bEditCps){
            if( (int)man->cps.size()>0 ){
                iCp = _clamp(iCp, 0, (int)man->cps.size()-1);
                man->cps[iCp].add(d);
                printf("SplineGUI: move Cp[%i] to (%g,%g,%g)\n", iCp, man->cps[iCp].x,man->cps[iCp].y,man->cps[iCp].z );
                man->updatePos();
            }
        }
        // When not editing control points, do NOT allow direct anchor movement
        // Only +/- should move along the spline
    }

};

#endif
