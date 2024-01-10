#ifndef MolecularDatabase_h
#define MolecularDatabase_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Mat3.h"
#include "Atoms.h"
#include "Vec3Utils.h"

#include <map>




class Descriptor : public MetaData{
    public:
    double* descs=0;

    Descriptor(){};
    Descriptor(int n_){
        alloc(n_);
    };
    ~Descriptor(){
        delete[] descs;
    };
    void alloc(int n_){
        dimensionOfDescriptorSpace=n_;
        descs = new double[dimensionOfDescriptorSpace];
    };
    void realloc(int n_){
        dimensionOfDescriptorSpace=n_;
        double* temp = new double[dimensionOfDescriptorSpace];
        memcpy(temp, this->descs, sizeof(double)*dimensionOfDescriptorSpace);
        delete[] descs;
        descs = temp;
    };

    void print(){
        printf("n: %d\n", dimensionOfDescriptorSpace);
        for(int i = 0; i < dimensionOfDescriptorSpace; i++){
            printf("%lf ", descs[i]);
        }
        printf("\n");
    };

    void calculatePrincipalAxes(Atoms* a, int atom_type){
        Vec3d* pos;
        int nbAtomsToPos=0;
        if(atom_type > 0){
            for(int i = 0; i < a->natoms; i++){
                if(a->atypes[i] == atom_type){
                    nbAtomsToPos++;
                }
            }
            pos = new Vec3d[nbAtomsToPos];
            int j = 0;
            for(int i =0; i< a->natoms; i++){
                if(a->atypes[i] == atom_type){
                    pos[j] = a->apos[i];
                    j++;
                }
            }
        }
        else{
            pos = a->apos;
            nbAtomsToPos = a->natoms;
        }
        if(nbAtomsToPos < 3){
            printf("Not enough atoms to calculate principal axes\n");
            delete [] pos;
            return;
            
        }
        this->dimensionOfDescriptorSpace += 3;
        if(!descs){
             alloc(dimensionOfDescriptorSpace);
        }else{
            realloc(dimensionOfDescriptorSpace);
        }
        Mat3d I;
        I.set(0.0);
        double* ws = new double[nbAtomsToPos];
        for(int i=0; i<nbAtomsToPos; i++){
            //ws[i] = a->atypes[i];
            ws[i] = 1.0;
        }
        Vec3d cog = Vec3dZero;
        getCOG(nbAtomsToPos, pos, ws, cog);
        for(int i=0; i<nbAtomsToPos; i++){
            Vec3d r;
            Mat3d E=Mat3dIdentity;
            r.set_sub(pos[i], cog);
            I.addOuter(r, r, -ws[i]);
            I.add_mul(E,r.norm2()*ws[i]);
        }
        Vec3d e;
        I.eigenvals(e);
        e.sort();
        descs[dimensionOfDescriptorSpace-3] = e.x;
        descs[dimensionOfDescriptorSpace-2] = e.y;
        descs[dimensionOfDescriptorSpace-1] = e.z;
        delete [] ws, pos;
    }

    void calculateNatoms(Atoms* a, int atom_type){
        this->dimensionOfDescriptorSpace += 1;
        if(!descs){
            alloc(dimensionOfDescriptorSpace);
        }else{
            realloc(dimensionOfDescriptorSpace);
        }
        if(atom_type > 0){
            int nbAtomsOfTheType=0;
            for(int i = 0; i < a->natoms; i++){
                if(a->atypes[i] == atom_type){
                    nbAtomsOfTheType++;
                }
            }
            descs[dimensionOfDescriptorSpace-1] = nbAtomsOfTheType;
        }
        else{
            descs[dimensionOfDescriptorSpace-1] = a->natoms;
        }
    }



    void chooseDescriptor(MolecularDatabaseDescriptor md, Atoms* a){
        switch (md.type)
        {
        case principal_axes:
            // for(int i=0; i<a->natoms; i++){
            //     printf("a->atypes[%d] = %d\n", i, a->atypes[i]);
            // }
            calculatePrincipalAxes(a, md.atom_type);
            break;
        case number_of_atoms:
            calculateNatoms(a, md.atom_type);
            break;
        default:
            break;
        }
    }

    bool compare(Descriptor* d){
        if(this->dimensionOfDescriptorSpace != d->dimensionOfDescriptorSpace){
            return false;
        }
        for(int i = 0; i < dimensionOfDescriptorSpace; i++){
            if(((this->descs[i] - d->descs[i])*(this->descs[i] - d->descs[i])) > 0.0001){
                return false;
            }
        }
        return true;
    }

    void copyOf(Descriptor* d){
        if(this->dimensionOfDescriptorSpace != d->dimensionOfDescriptorSpace){
            alloc(d->dimensionOfDescriptorSpace);
        }
        memcpy(this->descs, d->descs, sizeof(double)*d->dimensionOfDescriptorSpace);
    }

};





class MolecularDatabase : public MetaData{
private:
    Descriptor* descriptors = nullptr;
    Atoms* atoms = nullptr;
    int nMembers = 0;
    int nMaxCell = 0;
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
        this->descriptors = new Descriptor[n];
    };
    void realloc(int n){
        this->nMaxCell = n;
        Atoms* temp = new Atoms[n];
        memcpy(temp, this->atoms, sizeof(Atoms)*n);
        delete[] this->atoms;
        this->atoms = temp;
        Descriptor* temp2 = new Descriptor[n];
        memcpy(temp2, this->descriptors, sizeof(Descriptor)*n);
        delete[] this->descriptors;
        this->descriptors = temp2;        
    };
    void addMember(Atoms* structure, Descriptor* descriptor){
        if(nMembers >= nMaxCell){
            realloc(nMaxCell*2);
        }
        atoms[nMembers].copyOf(*structure);
        // for(int i = 0; i < nbOfusedDescriptors; i++){
        //     descriptors[nMembers].chooseDescriptor(usedDescriptors[i], &atoms[nMembers]);
        // }
        descriptors[nMembers].copyOf(descriptor);
        if(nMembers == 0)   dimensionOfDescriptorSpace = descriptors[nMembers].dimensionOfDescriptorSpace;
        this->nMembers++;
    };
    void  print(){
        printf("nMembers: %d\n", nMembers);
        printf("%-10s", "Member");
        for(int j = 0; j < dimensionOfDescriptorSpace; j++){
            printf(" %-10s", ("Desc" + std::to_string(j)).c_str());
        }
        printf("\n");

        for(int i = 0; i < nMembers; i++){
            printf("%-10d", i);
            for(int j = 0; j < dimensionOfDescriptorSpace; j++){
                printf(" %-10lf", descriptors[i].descs[j]);
            }
            printf("\n");
        }
    };
    int getNMembers(){   return nMembers;   };
    bool compareDescriptors(int i, int j){
        return descriptors[i].compare(&descriptors[j]);
    };

    void testHash(Atoms* a){
        Descriptor d;
        for(int i = 0; i < nbOfusedDescriptors; i++){
            d.chooseDescriptor(usedDescriptors[i], a);
        }
        hashDescriptors(&d);
        if(1){
            addMember(a, &d);
        }

    };

    int hashDescriptors(Descriptor* d){
        return 0;
    };
    //virtual char* tostr(int id){};    



};











#endif