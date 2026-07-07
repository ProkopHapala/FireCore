/// @file CrystalNpz.h
/// @brief Crystal geometry NPZ schema v1: pos, Z, bonds_ij.
/// Parity reference: web/common_js/npzIO.js readCrystalNpz / molToCrystalArrays
#ifndef CrystalNpz_h
#define CrystalNpz_h

#include "NpzIO.h"
#include "Molecule.h"
#include "MMFFparams.h"
#include <string>
#include <unistd.h>

struct CrystalData { public:
    int natoms = 0;
    std::vector<double> pos;   // 3N
    std::vector<int> Z;          // N
    int nbonds = 0;
    std::vector<int> bonds_ij; // 2*nbonds
};

inline int crystal_atomTypeFromZ(const MMFFparams& p, int Z){
    for(int i=0; i<(int)p.etypes.size(); i++){
        if(p.etypes[i].iZ == (uint8_t)Z) return p.getAtomType(p.etypes[i].name);
    }
    printf("CrystalNpz ERROR: unknown atomic number Z=%i\n", Z); exit(1);
    return -1;
}

inline void crystal_fill_molecule(Molecule& mol, const CrystalData& c, MMFFparams* params){
    if(!params){ printf("CrystalNpz ERROR: params required\n"); exit(1); }
    mol.bindParams(params);
    int nb = c.nbonds;
    mol.allocate(c.natoms, nb);
    for(int i=0; i<c.natoms; i++){
        mol.pos[i].set(c.pos[i*3], c.pos[i*3+1], c.pos[i*3+2]);
        mol.atomType[i] = crystal_atomTypeFromZ(*params, c.Z[i]);
        mol.REQs[i] = Quat4dZero;
        mol.npis[i] = -1;
    }
    for(int k=0; k<nb; k++){
        mol.bond2atom[k].set(c.bonds_ij[k*2], c.bonds_ij[k*2+1]);
        mol.bondType[k] = 0;
    }
}

inline CrystalData crystal_from_npz(const NpzFile& npz){
    CrystalData c;
    const NpyArray* ppos = npz.get("pos");
    const NpyArray* pZ   = npz.get("Z");
    if(!ppos || !pZ){ printf("CrystalNpz ERROR: missing pos/Z\n"); exit(1); }
    const double* pos = npy_as_f8(*ppos);
    const int32_t* Z  = npy_as_i4(*pZ);
    int n = 0;
    if(const NpyArray* pn = npz.get("natoms")){
        n = npy_as_i4(*pn)[0];
    }else{
        n = (int)pZ->count();
    }
    if(ppos->ndim != 2 || ppos->shape[0]!=n || ppos->shape[1]!=3){
        printf("CrystalNpz ERROR: pos shape mismatch natoms=%i\n", n); exit(1);
    }
    if((int)pZ->count() != n){ printf("CrystalNpz ERROR: Z length %zu != natoms %i\n", pZ->count(), n); exit(1); }
    c.natoms = n;
    c.pos.assign(pos, pos + n*3);
    c.Z.assign(Z, Z + n);
    if(const NpyArray* pb = npz.get("bonds_ij")){
        const int32_t* bij = npy_as_i4(*pb);
        int nb = (pb->ndim==2) ? pb->shape[0] : (int)(pb->count()/2);
        c.nbonds = nb;
        c.bonds_ij.assign(bij, bij + nb*2);
    }
    return c;
}

inline CrystalData loadCrystalNpz(const char* path){
    NpzFile npz = npz_read_file(path);
    return crystal_from_npz(npz);
}

inline int molecule_loadCrystalNpz(Molecule& mol, const char* path, MMFFparams* params, int verbosity=0){
    if(verbosity>0) printf("Molecule::loadCrystalNpz(%s)\n", path);
    CrystalData c = loadCrystalNpz(path);
    crystal_fill_molecule(mol, c, params);
    return c.natoms;
}

/// Standalone .npy fast path: stem must be pos, Z, or bonds_ij; caller merges via sidecar cache in BrowserView.
inline bool crystal_load_npy_stem(const char* path, const std::string& stem, CrystalData& acc){
    NpyArray arr = npy_load_file(path, true);
    if(stem=="pos"){
        const double* p = npy_as_f8(arr);
        int n = (arr.ndim==2) ? arr.shape[0] : (int)(arr.count()/3);
        acc.pos.assign(p, p + n*3);
        if(acc.natoms==0) acc.natoms = n;
        else if(acc.natoms != n){ printf("CrystalNpz ERROR: pos natoms mismatch\n"); exit(1); }
        return true;
    }
    if(stem=="Z"){
        const int32_t* z = npy_as_i4(arr);
        int n = (int)arr.count();
        acc.Z.assign(z, z + n);
        if(acc.natoms==0) acc.natoms = n;
        else if(acc.natoms != n){ printf("CrystalNpz ERROR: Z natoms mismatch\n"); exit(1); }
        return true;
    }
    if(stem=="bonds_ij"){
        const int32_t* b = npy_as_i4(arr);
        int nb = (arr.ndim==2) ? arr.shape[0] : (int)(arr.count()/2);
        acc.nbonds = nb;
        acc.bonds_ij.assign(b, b + nb*2);
        return true;
    }
    printf("CrystalNpz ERROR: unsupported .npy stem '%s' in '%s'\n", stem.c_str(), path);
    exit(1);
    return false;
}

inline int molecule_loadNpy(Molecule& mol, const char* path, MMFFparams* params, int verbosity=0){
    if(verbosity>0) printf("Molecule::loadNpy(%s)\n", path);
    std::string p(path);
    size_t slash = p.find_last_of("/\\");
    std::string dir = (slash==std::string::npos) ? "" : p.substr(0, slash+1);
    size_t dot   = p.find_last_of('.');
    if(dot==std::string::npos || dot < slash){ printf("CrystalNpz ERROR: bad .npy path '%s'\n", path); exit(1); }
    size_t name0 = (slash==std::string::npos) ? 0 : slash+1;
    std::string stem = p.substr(name0, dot-name0);
    CrystalData c;
    crystal_load_npy_stem(path, stem, c);
    if(c.pos.empty() || c.Z.empty()){
        std::string sib = dir + (stem=="pos" ? "Z.npy" : (stem=="Z" ? "pos.npy" : ""));
        if(!sib.empty()){
            std::string sstem = (stem=="pos") ? "Z" : "pos";
            crystal_load_npy_stem(sib.c_str(), sstem, c);
        }else if(stem=="bonds_ij"){
            crystal_load_npy_stem((dir+"pos.npy").c_str(), "pos", c);
            crystal_load_npy_stem((dir+"Z.npy").c_str(), "Z", c);
        }
    }
    if(c.pos.empty() || c.Z.empty()){
        printf("CrystalNpz ERROR: need pos.npy + Z.npy (got stem '%s' in '%s')\n", stem.c_str(), path);
        exit(1);
    }
    std::string bondsPath = dir + "bonds_ij.npy";
    if(c.bonds_ij.empty() && access(bondsPath.c_str(), R_OK)==0)
        crystal_load_npy_stem(bondsPath.c_str(), "bonds_ij", c);
    crystal_fill_molecule(mol, c, params);
    return c.natoms;
}

#endif
