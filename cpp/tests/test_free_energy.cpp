#include <gtest/gtest.h>

#include "globals.h"
#include "testUtils.h"


#include "MolWorld_sp3.h"



char dir_cpp[1024];
void init(MolWorld_sp3& W){
	W.smile_name = nullptr;
	W.xyz_name = "../common_resources/xyz/pyridine";
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
    init(W);
    double initial = 0.5;
    double final   = 4.0;
    int N = 1;
    int colectiveVariable[N] = {0};
    double E = W.compute_Free_energy(initial, final, N, colectiveVariable, 100, 100000, 10000, 100, 100, 0.5);
    printf("E = %f\n", E);
    EXPECT_DOUBLE_EQ(E, 0.08241649205729196);
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

