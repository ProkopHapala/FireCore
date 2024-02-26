
#ifndef GOpt_h
#define GOpt_h

#include <stdlib.h>
#include <stdio.h>
//#include <string.h>
#include <vector>

#include "IO_utils.h"

//#include "testUtils.h"
#include "fastmath.h"
#include "Vec3.h"
//#include "Mat3.h"
//#include "Vec3Utils.h"
#include "constants.h"
#include "Forces.h"
//#include "MMFFsp3_loc.h"
//#include "UFF.h"

#include "constrains.h"
//#include "molecular_utils.h"

struct GOpt{
    bool bExploring = false;
    int  istep   =0;
    int  nExplore=0;
    int  nRelax  =0;
    double vel_kick = 1.0;
    double pos_kick = 0.25;

    double bThermalSampling = 0.0;   // if >0 then we do thermal sampling
    double T_target         = 300.0;   // target temperature for thermal sampling (if bThermalSampling>0)
    double T_current        = 0.0;   // current temperature for thermal sampling (if bThermalSampling>0)
    double gamma_damp       = 0.01;  // damping factor for thermal sampling (if bThermalSampling>0)

    Constrains constrs;

    void startExploring(){
        printf( "GOpt::startExploring()\n" );
        bExploring = true;
        istep=0;
        constrs.update_drives();
        constrs.update_drives();
    }

    bool update(){
        istep++;
        if(bExploring){
            if(istep>=nExplore){ 
                printf( "GOpt::update() stop exploring istep(%i)>nExplore(%i) \n" );
                bExploring=false; istep=0; return true; 
            }
        }else{
           // if(istep>=nRelax  ){ bExploring=true; istep=0; return true; }
        }
        //bExploring       = go.bExploring; 
        //bThermalSampling = bExploring;
        //bConstrains      = bExploring;
        //if(bExploring) bConverged = false;
        return false;
    }

    void apply_kick( int n, Vec3d* pos=0, Vec3d* vel=0 ){
        //printf( "GOpt::apply_kick() n=%i |pos|=%g |vel=%g|\n", n, pos_kick, vel_kick );
        Vec3d vcog = Vec3dZero;
        for(int i=0; i<n; i++){
            if(pos) pos[i].add( randf(-pos_kick,pos_kick), randf(-pos_kick,pos_kick), randf(-pos_kick,pos_kick) );
            if(vel){
                Vec3d dv = Vec3d{ randf(-vel_kick,vel_kick), randf(-vel_kick,vel_kick), randf(-vel_kick,vel_kick) };
                vel[i].add( dv );
                vcog  .add( dv );
            }
        }
        vcog.mul( 1.0/n );
        if(vel)for(int i=0; i<n; i++){
            vel[i].sub( vcog );
        }
    }

    void clear(){
        bExploring = false;
        istep=0;
        constrs.clear();
    }

    void copy(const GOpt& go){
        //printf( "GOpt::copy() \n" );
        nExplore = go.nExplore;
        nRelax   = go.nRelax;
        vel_kick = go.vel_kick;
        pos_kick = go.pos_kick;
        T_target = go.T_target;
        gamma_damp = go.gamma_damp;
        constrs.copy( go.constrs );
        printf( "GOpt::copy() nExplore=%i nRelax=%i vel_kick=%g pos_kick=%g T_target=%g gamma_damp=%g \n", nExplore, nRelax, vel_kick, pos_kick, T_target, gamma_damp );
    }


    void print(){ printf( "GOpt::print() nExplore=%i nRelax=%i vel_kick=%g pos_kick=%g T_target=%g gamma_damp=%g \n", nExplore, nRelax, vel_kick, pos_kick, T_target, gamma_damp ); }

};

#endif
