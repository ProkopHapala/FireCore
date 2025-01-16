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

    W.init();
}

// test that the temperature is enough close to the desired value
TEST(MMFFsp3_loc, test_move_atom_Langevin){
    MolWorld_sp3 W;
    verbosity = -1;
    char path_H2O[1024];
    snprintf(path_H2O, sizeof(path_H2O), "%s/common_resources/xyz/H2O", dir_cpp);
    init(W, path_H2O);
    W.ffl.apos[0] = {1.0, 0.0, 0.0};

    int N = 10000;
    int Neq = 100000;

    // store desired data for temperature
    double T_desired = 300.0;
    double T_measured = 0.0;
    double dt = 0.5;
    double tdamp = 20;
    double gamma = 1/(dt*tdamp);

    for(int i = 0; i < N+Neq; i++){
        W.ffl.fapos[0] = {0.0, 0.0, 0.0};

        Vec3d d = W.ffl.apos[0];
        double l = d.norm();
        double f;
        Vec2d Ks = {0.005, 0.005};
        spring(l, Vec2dZero, Ks, 100000, f);
        W.ffl.fapos[0].add_mul(d, f/l);        
        W.ffl.move_atom_Langevin(0, dt, 10000, gamma, T_desired);

        if(i>Neq)T_measured += W.ffl.vapos[0].norm2() / (3.0 * const_kB * 1);
    }
    T_measured /= N;
    EXPECT_NEAR(T_desired, T_measured, 10.0);
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