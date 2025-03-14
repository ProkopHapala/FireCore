#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <math.h>
#include <functional>

#include "fastmath.h"
#include "Vec3.h"
#include "VecN.h"
#include "DynamicOpt.h"
#include "InteractionsGauss.h"
#include "eFF.h"
#include "argparse.h"

#define R2SAFE 1.0e-8f
//#define ExportFunc __declspec(dllexport)

#ifdef _WIN32
#define EXPORT_API __declspec(dllexport)
#else
#define EXPORT_API
#endif



class ConsoleEFF {
public:
    EFF ff;
    DynamicOpt opt;
    int perFrame = 10;
    bool bRun = true;

    void initEFFsystem(const char* fname, bool bTestEval = true) {
        ff.loadFromFile_fgo(fname);
        opt.bindOrAlloc(ff.nDOFs, ff.pDOFs, nullptr, ff.fDOFs, nullptr);
        opt.initOpt(0.0015, 0.01);
        opt.f_limit = 1000.0;
        ff.iPauliModel = 1;
        
        if(bTestEval) {
            double E = ff.eval();
            ff.printEnergies();
        }
    }

    void runSimulation() {

        double F2 = 1.0;
        for(int itr=0; itr<perFrame; itr++) {
            ff.clearForce();
            double Etot = ff.eval();
            F2 = opt.move_FIRE();

            printf("Step %d: E_total = %.6f, |F| = %.6f\n", itr, Etot, sqrt(F2));
            ff.printEnergies();
            
            if(F2 < 1e-6) {
                printf("Convergence reached\n");
                bRun = false;
                break;
            }
        }
        //printSystemState();
    
    }

    void printSystemState() {
        printf("\n=== Atomic Positions ===\n");
        for(int i=0; i<ff.na; i++) {
            printf("Atom %d: (%.3f, %.3f, %.3f)\n", i, ff.apos[i].x, ff.apos[i].y, ff.apos[i].z);
        }

        printf("\n=== Electron Positions/Sizes ===\n");
        for(int i=0; i<ff.ne; i++) {
            printf("Electron %d: Pos(%.3f, %.3f, %.3f) Size=%.3f\n",
                i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z, ff.esize[i]);
        }
    }
};

ConsoleEFF eff;

int main(int argc, char *argv[]) {
    
    if(argc > 1) {
        eff.initEFFsystem(argv[1]);
        while (eff.bRun) {
            eff.runSimulation();
        }
    } else {
        printf("Usage: %s <input.fgo>\n", argv[0]);
    }
    
    return 0;
}

extern "C" EXPORT_API int* unityInit(const char* fileName) {
    eff.initEFFsystem(fileName);
    return new int[2]{eff.ff.ne, eff.ff.na};
}

// Export the frame update function returning float array pointer
extern "C" EXPORT_API float* unityNextFrame(int* size) {
    eff.runSimulation();
    
    // Calculate total size needed (3 coordinates per position)
    *size = ((eff.ff.ne + eff.ff.na) * 3) + eff.ff.ne;
    float* positions = new float[*size];
    
    int idx = 0;
    // Copy electron positions
    for(int i=0; i < eff.ff.ne; i++) {
        positions[idx++] = eff.ff.epos[i].x;
        positions[idx++] = eff.ff.epos[i].y;
        positions[idx++] = eff.ff.epos[i].z;
    }
    
    // Copy atom positions
    for(int i=0; i < eff.ff.na; i++) {
        positions[idx++] = eff.ff.apos[i].x;
        positions[idx++] = eff.ff.apos[i].y;
        positions[idx++] = eff.ff.apos[i].z;
    }

    for (int i = 0; i < eff.ff.ne; i++) {
        positions[idx++] = eff.ff.esize[i];
    }
    
    return positions;
}

// Add cleanup function to prevent memory leaks
extern "C" EXPORT_API void cleanupPositions(float* positions) {
    delete[] positions;
}
