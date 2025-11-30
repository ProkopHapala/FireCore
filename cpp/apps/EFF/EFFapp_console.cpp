#include <memory>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <string>
#include <vector>
#include <math.h>
#include <functional>
#include <list>

#include "fastmath.h"
#include "Vec3.h"
#include "VecN.h"
#include "DynamicOpt.h"
#include "InteractionsGauss.h"
#include "eFF.h"
#include "argparse.h"
#include "Forces.h"

#define R2SAFE 1.0e-8f
//#define ExportFunc __declspec(dllexport)

#ifdef _WIN32
#define EXPORT_API __declspec(dllexport)
#else
#define EXPORT_API __attribute__((visibility("default")))
#endif



class ConsoleEFF {
public:
    EFF ff;
    DynamicOpt opt;
    int perFrame = 10;
    bool bRun = true;

    std::vector<int> fixedAtoms = {};
    std::vector<int> fixedElectrons = {};

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
    
    // Set position of a single atom
    bool setAtomPosition(int atomIndex, float x, float y, float z) {
        if (atomIndex < 0 || atomIndex >= ff.na) {
            printf("Error: Atom index %d out of range (0-%d)\n", atomIndex, ff.na-1);
            return false;
        }
        
        ff.apos[atomIndex].x = x;
        ff.apos[atomIndex].y = y;
        ff.apos[atomIndex].z = z;
        
        // No need to update DOFs as apos directly points to the DOFs array
        return true;
    }
    
    // Set position and size of a single electron
    bool setElectronPosition(int electronIndex, float x, float y, float z, float size) {
        if(setElectronPosition(electronIndex, x, y, z)) {
            ff.esize[electronIndex] = size;
            return true;
        }
        return false;
    }
    bool setElectronPosition(int electronIndex, float x, float y, float z) {
        if (electronIndex < 0 || electronIndex >= ff.ne) {
            printf("Error: Electron index %d out of range (0-%d)\n", electronIndex, ff.ne-1);
            return false;
        }
        
        ff.epos[electronIndex].x = x;
        ff.epos[electronIndex].y = y;
        ff.epos[electronIndex].z = z;
        
        
        // No need to update DOFs as epos and esize directly point to the DOFs array
        return true;
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
    
    // Calculate the size of the return array: 2 (for ne and na) + ne (for espin array)
    int returnSize = 2 + eff.ff.ne;
    int* returnArray = new int[returnSize];
    
    // Set the first two elements (ne and na)
    returnArray[0] = eff.ff.ne;
    returnArray[1] = eff.ff.na;
    
    // Append the espin array
    for (int i = 0; i < eff.ff.ne; i++) {
        returnArray[2 + i] = eff.ff.espin[i];
    }
    
    return returnArray;
}

// Export the frame update function returning float array pointer
extern "C" EXPORT_API float* unityNextFrame(int* size) {
    // save previous positions for fixed particles before the simulation step
    int nfa = eff.fixedAtoms.size();
    int nfe = eff.fixedElectrons.size();
    
    std::vector<Vec3d> fixedAPositions;
    std::vector<Vec3d> fixedEPositions;

    for(int i=0; i < nfe; i++) {
        fixedEPositions.push_back(eff.ff.epos[eff.fixedElectrons[i]]);
    }

    for(int i=0; i < nfa; i++) {
        fixedAPositions.push_back(eff.ff.apos[eff.fixedAtoms[i]]);
    }


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

    // set the fixed original positions of fixed paricles
    for(int i=0; i < nfe; i++) {
        eff.setElectronPosition(eff.fixedElectrons[i], fixedEPositions[i].x, fixedEPositions[i].y, fixedEPositions[i].z, eff.ff.esize[eff.fixedElectrons[i]]);
    }

    for(int i=0; i < nfa; i++) {
        eff.setAtomPosition(eff.fixedAtoms[i], fixedAPositions[i].x, fixedAPositions[i].y, fixedAPositions[i].z);
    }
    
    return positions;
}

// Add cleanup function to prevent memory leaks
extern "C" EXPORT_API void cleanupPositions(float* positions) {
    delete[] positions;
}

// Export function to set a single atom position
extern "C" EXPORT_API void unitySetAtomPosition(int atomIndex, float x, float y, float z) {
    eff.setAtomPosition(atomIndex, x, y, z);
}

extern "C" EXPORT_API void fixAtom(int atomIndex) {
    eff.fixedAtoms.push_back(atomIndex);
}

extern "C" EXPORT_API void unfixAtom(int atomIndex) {
    for (int i=0; i<eff.fixedAtoms.size(); i++) {
        if (eff.fixedAtoms[i] == atomIndex) {
            eff.fixedAtoms.erase(eff.fixedAtoms.begin()+i);
        }
    }
}

extern "C" EXPORT_API void fixElectron(int electronIndex) {
    eff.fixedElectrons.push_back(electronIndex);
}

extern "C" EXPORT_API void unfixElectron(int electronIndex) {
    for (int i=0; i<eff.fixedElectrons.size(); i++) {
        if (eff.fixedElectrons[i] == electronIndex) {
            eff.fixedElectrons.erase(eff.fixedElectrons.begin()+i);
        }
    }
}

// Export function to set a single electron position and size
extern "C" EXPORT_API void unitySetElectronPosition(int electronIndex, float x, float y, float z, float size) {
    eff.setElectronPosition(electronIndex, x, y, z, size);
}

extern "C" EXPORT_API int* unityGetProtonNumbers() {
    int* protonNumbers = new int[eff.ff.na];

    for(int i = 0; i < eff.ff.na; i++) {
        protonNumbers[i] = static_cast<int>(eff.ff.aPars[i].x);
    }

    return protonNumbers;
}

// Add cleanup function for proton numbers to prevent memory leaks
extern "C" EXPORT_API void cleanupProtonNumbers(int* protonNumbers) {
    delete[] protonNumbers;
}

// Add cleanup function for unityInit return value to prevent memory leaks
extern "C" EXPORT_API void cleanupInitData(int* data) {
    delete[] data;
}
