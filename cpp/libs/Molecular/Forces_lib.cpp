
#include <functional>
#include <unordered_map>
#include <string>


#include "globals.h"
#include "fastmath.h"
#include "Forces.h"
#include "testUtils.h"


using ForceEvalFunc = std::function<void(int, double*, Vec2d*, double*)>;

std::unordered_map<std::string, ForceEvalFunc> evaluators;

//============================

extern "C"{

// A function to set up our map of evaluators
int init() {
    printf("initialize_evaluators() \n");

    // --- Lennard-Jones + Coulomb (double precision) ---
    printf("initialize_evaluators() 'getLJQ' params: [R, E, Q, R2damp] \n");
    evaluators["getLJQ"] = 
        [](int n, double* xs, Vec2d* FEs, double* params) {
        // Unpack parameters from the array
        Quat4d REQH;
        REQH.x = params[0]; // R (equilibrium distance)
        REQH.y = params[1]; // E (well depth)
        REQH.z = params[2]; // Q (charge)
        double R2damp = params[3];

        for (int i = 0; i < n; ++i) {
            Vec3d dp = {xs[i], 0.0, 0.0}; // Assuming 1D evaluation along x-axis
            Vec3d force_vec;
            double energy = getLJQ(dp, force_vec, REQH, R2damp);
            FEs[i] = {energy, force_vec.x}; // Store energy and the x-component of the force
        }
    };

    // --- Morse + Coulomb + Hydrogen Bond ---
    printf("initialize_evaluators() 'getMorseQH' params: [R, E, Q, H, K, R2damp] \n");
    evaluators["getMorseQH"] = 
        [](int n, double* xs, Vec2d* FEs, double* params) {
        // Unpack parameters
        Quat4d REQH;
        REQH.  x      = params[0]; // R
        REQH.  y      = params[1]; // E
        REQH.  z      = params[2]; // Q
        REQH.  w      = params[3]; // H-bond parameter
        double K      = params[4];      // Morse stiffness
        double R2damp = params[5];
        for (int i = 0; i < n; ++i) {
            Vec3d dp = {xs[i], 0.0, 0.0};
            Vec3d force_vec;
            double energy = getMorseQH(dp, force_vec, REQH, K, R2damp);
            FEs[i] = {energy, force_vec.x};
        }
    };
    

    // --- Morse + Coulomb + Hydrogen Bond ---
    printf("initialize_evaluators() 'getMorseP4' params: [R, E, K] \n");
    evaluators["getMorseP4"] = 
        [](int n, double* xs, Vec2d* FEs, double* params) {
        // Unpack parameters
        double R     = params[0]; // R
        double E     = params[1]; // E
        double K     = params[2];      // Morse stiffness
        double R2damp = params[3];
        for (int i = 0; i < n; ++i) {
            Vec3d dp = {xs[i], 0.0, 0.0};
            Vec3d force_vec;
            double energy = getMorseP4(dp, force_vec, R, E, K);
            FEs[i] = {energy, force_vec.x};
        }
    };


    return  evaluators.size();
}

void evalForce( char* funcName, int n, double* xs, Vec2d* FEs, double* params ){
    //initialize_evaluators();
    long T0 = getCPUticks();
    if( evaluators.find(funcName) == evaluators.end() ){ printf("ERROR: evalForce() funcName=%s NOT FOUND!!! =\u003e exit\n", funcName ); };
    evaluators[funcName](n,xs,FEs,params);
    long T1 = getCPUticks();
    printf("evalForce(%s,%i) %g [tick/eva] \n", funcName, n, ((double)(T1-T0))/n );
}

} // extern "C"
