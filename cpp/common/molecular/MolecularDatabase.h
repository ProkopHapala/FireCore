#ifndef MolecularDatabase_h
#define MolecularDatabase_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Mat3.h"
#include "Atoms.h"
#include "Vec3Utils.h"


class MolecularDatabase : public MetaData{
private:
    Vec3d* descriptors;
    Atoms* atoms;
    int nMembers;
    int nMaxCell;
public:
    MolecularDatabase(){
        nMaxCell = 100;
        alloc(nMaxCell);
    };
    MolecularDatabase(int n){
        alloc(n);
    };
    void alloc(int n){
        this->nMembers = 0;
        this->nMaxCell = n;
        this->atoms = new Atoms[n];
        this->descriptors = new Vec3d[n];
    };
    void realloc(int n){
        this->nMaxCell = n;
        Atoms* temp = new Atoms[n];
        memcpy(temp, this->atoms, sizeof(Atoms)*n);
        delete[] this->atoms;
        this->atoms = temp;
        Vec3d* temp2 = new Vec3d[n];
        memcpy(temp2, this->descriptors, sizeof(Vec3d)*n);
        delete[] this->descriptors;
        this->descriptors = temp2;        
    };
    void addMember(Atoms structure){
        if(nMembers >= nMaxCell){
            realloc(nMaxCell*2);
        }
        atoms[nMembers].copyOf(structure);
        calculatePrincipalAxes(&structure);
        this->nMembers++;
    };
    void  print(){
        printf("nMembers: %d\n", nMembers);
        printf("%-10s %-10s %-10s %-10s %-10s\n", "Member", "Ax_1", "Ax_2", "Ax_3", "nAtoms");
        for(int i = 0; i < nMembers; i++){
            printf("%-10d %-10f %-10f %-10f %-10d\n", i, descriptors[i].x, descriptors[i].y, descriptors[i].z, atoms[i].natoms);
        }
    };
    //virtual char* tostr(int id){};    
    void calculatePrincipalAxes(Atoms* a){
        Mat3d I;
        I.set(0.0);
        int natoms = a->natoms;
        double* ws = new double[natoms];
        for(int i=0; i<natoms; i++){
            //ws[i] = atypes[i];
            ws[i] = 1.0;
        }
        Vec3d cog = Vec3dZero;
        getCOG(natoms, a->apos, ws, cog);
        for(int i=0; i<natoms; i++){
            Vec3d r;
            Mat3d E=Mat3dIdentity;
            r.set_sub(a->apos[i], cog);
            I.addOuter(r, r, -ws[i]);
            I.add_mul(E,r.norm2()*ws[i]);
        }
        I.eigenvals(descriptors[nMembers]);
        descriptors[nMembers].sort();
        delete [] ws;
        printf("%lf %lf %lf\n", descriptors[nMembers].x, descriptors[nMembers].y, descriptors[nMembers].z);
    }


};











#endif