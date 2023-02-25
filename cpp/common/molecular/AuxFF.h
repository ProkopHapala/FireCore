
#ifndef AuxFF_h
#define AuxFF_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"

// ======================
// ====   AuxFF
// ======================
/*

Couples MMFFf4 with Non-Bonded interactions
There can be additional centers ( e.g. bond-charges, Pi-orbitals, Free-electron pairs attached to atoms )

*/





class AuxFF{ public:

    // Bond Centers
    Vec2i* bond;   // center on bond
    float* bondf;  // position of center on bond

    // Pi Center (there are two atoms )
    int*   bond;   // center on bond
    float* bondf;  // position of center on bond

    Quat4f*  REQs =0;  // [natom] parameters of non-covalent interactions



    void projectBond(){
    }

    void projectPis(){
    }

};


#endif
