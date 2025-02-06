#include "globals.h"
#include "testUtils.h"

#include "MolWorld_sp3.h"



char dir_cpp[1024];
void init(MolWorld_sp3& W, const char* xyz_name, const char* constr_name=nullptr){
	W.smile_name = nullptr;
	W.xyz_name   = xyz_name;
	W.surf_name  = nullptr;
    W.constr_name= constr_name;
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
    
    verbosity = -1;
    W.initParams( path_ElementTypes, path_AtomTypes, path_BondTypes, path_AngleTypes, path_DihedralTypes);

    W.init();
}

void mexican_hat(double lamda1, double lamda2, int nbStep, int nMDsteps, int nEQsteps, double tdamp, double T, double dt, int nbProdSteps, int nrealization, int nSampleSteps){
    char xyz_path [1024];
    snprintf(xyz_path, sizeof(xyz_path), "%s/common_resources/xyz/H2O", dir_cpp);


    MolWorld_sp3 W;
    init(W, xyz_path);

    // run mexican hat with TI
    W.mexican_hat_TI(lamda1, lamda2, nbStep, nMDsteps, nEQsteps, tdamp, T, dt);

    // run mexican hat with JE
    W.mexican_hat_JE(lamda1, lamda2, nbProdSteps, nrealization, nEQsteps, tdamp, T, dt,  nSampleSteps);
}

void three_atoms_problem(double lamda1, double lamda2, int nbStep, int nMDsteps, int nEQsteps, double tdamp, double T, double dt, int nbProdSteps, int nrealization, int nSampleSteps){
    char xyz_path [1024];
    snprintf(xyz_path, sizeof(xyz_path), "%s/common_resources/xyz/H2O", dir_cpp);

    MolWorld_sp3 W;
    init(W, xyz_path);

    // run three atom problem with TI
    W.three_atoms_problem_TI(lamda1, lamda2, nbStep, nMDsteps, nEQsteps, tdamp, T, dt);
    // run three atom problem with JE
    W.three_atoms_problem_JE(lamda1, lamda2, nbProdSteps, nrealization, nEQsteps, tdamp, T, dt,  nSampleSteps);
}

void entropic_spring(double lamda1, double lamda2, int n, int *dc, int nbStep, int nMDsteps, int nEQsteps, double tdamp, double T, double dt, int nbProdSteps, int nrealization, int nSampleSteps, int natoms){
    if(!dc){
        dc = new int[n];
        for(int i=0; i<n; i++){
            dc[i] = i;
            printf("%d\n", dc[i]);
        }
    }

    
    char xyz_path [1024];
    snprintf(xyz_path, sizeof(xyz_path), "%s/common_resources/entropic_spring_%d", dir_cpp, natoms);

    char constr_path [1024];
    snprintf(constr_path, sizeof(constr_path), "%s/common_resources/entropic_spring_%d.cons", dir_cpp, natoms);

    MolWorld_sp3 W;
    init(W, xyz_path, constr_path);

    // run entropic spring with TI
    W.entropic_spring_TI(lamda1, lamda2, n, dc, nbStep, nMDsteps, nEQsteps, tdamp, T, dt);
    // run entropic spring with JE
    W.entropic_spring_JE(lamda1, lamda2, n, dc, nbProdSteps, nrealization, nEQsteps, tdamp, T, dt,  nSampleSteps);
}



int main(int argc, char **argv)
{
    printf("generate data for toy models for free energy\n");

    char* current_dir = getcwd(NULL, 0);
    char* cpp_pos = strstr(current_dir, "/cpp");
    if(cpp_pos) {*(cpp_pos+4) = '\0';}
    else { "Error: cannot find /cpp in current directory path -> Cannot run program"; exit(1); }
    snprintf(dir_cpp, sizeof(dir_cpp), "%s", current_dir);
    free(current_dir);



    //mexican_hat(std::stod(argv[1]), std::stod(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]), std::stoi(argv[11]));
    //three_atoms_problem(std::stod(argv[1]), std::stod(argv[2]), std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]), std::stoi(argv[11]));
    entropic_spring(std::stod(argv[1]), std::stod(argv[2]), 1, nullptr, std::stoi(argv[3]), std::stoi(argv[4]), std::stoi(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]), std::stoi(argv[9]), std::stoi(argv[10]), std::stoi(argv[11]), std::stoi(argv[12]));
    return 0;
}