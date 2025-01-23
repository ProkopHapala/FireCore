#include <gtest/gtest.h>

#include "globals.h"
#include "testUtils.h"


#include "MolWorld_sp3.h"



char dir_cpp[1024];
void init(MolWorld_sp3& W, const char* xyz_name){
	W.smile_name = nullptr;
	W.xyz_name = xyz_name;
	W.surf_name  = nullptr;
    W.constr_name= nullptr;
	W.bMMFF      = true;
    W.bEpairs    = false;
    W.gridStep   = 0.01;
    W.nPBC       = Vec3i({1,1,1});
    W.bUFF       = false; 
    W.b141       = true;
    W.bSimple    = false;
    W.bConj      = true;
    W.bCumulene  = true;
    W.bGridFF    = false;



    char path_ElementTypes[1024];
    char path_AtomTypes[1024];
    char path_BondTypes[1024];
    char path_AngleTypes[1024];
    char path_DihedralTypes[1024];
    snprintf(path_ElementTypes, sizeof(path_ElementTypes), "%s/common_resources/ElementTypes.dat", dir_cpp);
    snprintf(path_AtomTypes, sizeof(path_AtomTypes), "%s/common_resources/AtomTypes.dat", dir_cpp);
    snprintf(path_BondTypes, sizeof(path_BondTypes), "%s/common_resources/BondTypes.dat", dir_cpp);
    snprintf(path_AngleTypes, sizeof(path_AngleTypes), "%s/common_resources/AngleTypes.dat", dir_cpp);
    snprintf(path_DihedralTypes, sizeof(path_DihedralTypes), "%s/common_resources/DihedralTypes.dat", dir_cpp);

    W.initParams( path_ElementTypes, path_AtomTypes, path_BondTypes, path_AngleTypes, path_DihedralTypes);


    printf("INITTIALIZING\n");
    W.init();
}

TEST(MolWorld_sp3, initParams) {
    MolWorld_sp3 W;
    char path_ElementTypes[1024];
    char path_AtomTypes[1024];
    char path_BondTypes[1024];
    char path_AngleTypes[1024];
    char path_DihedralTypes[1024];
    snprintf(path_ElementTypes, sizeof(path_ElementTypes), "%s/common_resources/ElementTypes.dat", dir_cpp);
    snprintf(path_AtomTypes, sizeof(path_AtomTypes), "%s/common_resources/AtomTypes.dat", dir_cpp);
    snprintf(path_BondTypes, sizeof(path_BondTypes), "%s/common_resources/BondTypes.dat", dir_cpp);
    snprintf(path_AngleTypes, sizeof(path_AngleTypes), "%s/common_resources/AngleTypes.dat", dir_cpp);
    snprintf(path_DihedralTypes, sizeof(path_DihedralTypes), "%s/common_resources/DihedralTypes.dat", dir_cpp);

    W.initParams( path_ElementTypes, path_AtomTypes, path_BondTypes, path_AngleTypes, path_DihedralTypes);
    ASSERT_GT(W.params.atypes.size(), 0);
}


TEST(MolWorld_sp3, mexican_hat_TI ){
    MolWorld_sp3 W;
    char path_pyridine[1024];
    snprintf(path_pyridine, sizeof(path_pyridine), "%s/common_resources/xyz/pyridine", dir_cpp);
    init(W, path_pyridine);
    double initial = 0.5;
    double final   = 4.0;

    double E = W.mexican_hat_TI(initial, final, 100, 100000, 10000, 100, 100, 0.5);
    printf("E = %f\n", E);
    EXPECT_NEAR(E, -0.018800922191230408, 0.001);
}

TEST(MolWorld_sp3, three_atoms_problem_TI){
    MolWorld_sp3 W;
    char path_three_atoms[1024];
    snprintf(path_three_atoms, sizeof(path_three_atoms), "%s/common_resources/three_atoms", dir_cpp);
    init(W, path_three_atoms);
    double initial = 0.0;
    double final   = 5.0;
    double eps = 0.1;
    for(int i=0; i<W.ffl.natoms; i++){
        W.ffl.REQs[i].y = 0.316228; // this number corresponds to write 0.1 in AtomTypes Hydrogen EvdW
    }

    double E = W.three_atoms_problem_TI(initial, final, 100, 100000, 10000, 100, 300.0, 0.5);
    printf("E = %f\n", E);
    EXPECT_NEAR(E, 0.026256787647647592, 0.001);
}

int main(int argc, char **argv) {    
    ::testing::InitGoogleTest(&argc, argv);

    char* current_dir = getcwd(NULL, 0);
    char* cpp_pos = strstr(current_dir, "/cpp");
    if(cpp_pos) {*(cpp_pos+4) = '\0';}
    else { "Error: cannot find /cpp in current directory -> Cannot run tests"; exit(1); }
    snprintf(dir_cpp, sizeof(dir_cpp), "%s", current_dir);
    free(current_dir);

    return RUN_ALL_TESTS();
}

