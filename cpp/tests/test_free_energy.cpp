#include <gtest/gtest.h>

#include "testUtils.h"


#include "MolWorld_sp3.h"

// def init(
//     xyz_name  =None, 
//     surf_name =None, 
//     smile_name=None, 
//     constr_name=None,
//     sElementTypes  = "data/ElementTypes.dat",
//     sAtomTypes     = "data/AtomTypes.dat", 
//     sBondTypes     = "data/BondTypes.dat", 
//     sAngleTypes    = "data/AngleTypes.dat",
//     sDihedralTypes = "data/DihedralTypes.dat",
//     bMMFF=True, 
//     bEpairs=False,  
//     nPBC=(1,1,0), 
//     gridStep=0.1,
//     bUFF=False,
//     b141=True,
//     bSimple=False,
//     bConj=True,
//     bCumulene=True
// ):
// global glob_bMMFF
// glob_bMMFF = bMMFF
// nPBC=np.array(nPBC,dtype=np.int32)
// return lib.init( cstr(xyz_name), cstr(surf_name), cstr(smile_name), cstr(constr_name), bMMFF, bEpairs, bUFF, b141, bSimple, bConj, bCumulene, nPBC, gridStep, cstr(sElementTypes), cstr(sAtomTypes), cstr(sBondTypes), cstr(sAngleTypes), cstr(sDihedralTypes) )
int main(int argc, char **argv) {
    // ::testing::InitGoogleTest(&argc, argv);
    MolWorld_sp3 W;
    W.init();
    char* current_dir = getcwd(NULL, 0);
    printf("%s", current_dir);

    return 0;//RUN_ALL_TESTS();
}

