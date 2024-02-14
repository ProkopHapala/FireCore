#ifndef MolecularDatabase_h
#define MolecularDatabase_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Mat3.h"
#include "Atoms.h"
#include "Vec3Utils.h"
#include "MMFFparams.h"
#include "AtomicConfiguration.h"

#include <map>

enum DescriptorType {principal_axes, number_of_atoms};

struct MolecularDatabaseDescriptor
{
    DescriptorType type;
    int atom_type = -1;
};


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

    void calculatePrincipalAxes(Atoms* a, int atom_type, MMFFparams* params=0){
        Vec3d* pos;
        int nbAtomsToPos=0;
        if(atom_type > 0){
            for(int i = 0; i < a->natoms; i++){
                int ityp = a->atypes[i];
                if(params) ityp = params->atypes[ityp].element;
                if( ityp == atom_type){
                    nbAtomsToPos++;
                }
            }
            pos = new Vec3d[nbAtomsToPos];
            int j = 0;
            for(int i =0; i< a->natoms; i++){
                int ityp = a->atypes[i];
                if(params) ityp = params->atypes[ityp].element;
                if(ityp == atom_type){
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

    void calculateNatoms(Atoms* a, int atom_type, MMFFparams* params=0){
        this->dimensionOfDescriptorSpace += 1;
        if(!descs){
            alloc(dimensionOfDescriptorSpace);
        }else{
            realloc(dimensionOfDescriptorSpace);
        }
        if(atom_type > 0){
            int nbAtomsOfTheType=0;
            for(int i = 0; i < a->natoms; i++){
                int ityp = a->atypes[i];
                if(params) ityp = params->atypes[ityp].element;
                if( ityp == atom_type){
                    nbAtomsOfTheType++;
                }
            }
            descs[dimensionOfDescriptorSpace-1] = nbAtomsOfTheType;
        }
        else{
            descs[dimensionOfDescriptorSpace-1] = a->natoms;
        }
    }



    void chooseDescriptor(MolecularDatabaseDescriptor md, Atoms* a, MMFFparams* params=0){
        switch (md.type)
        {
        case principal_axes:
            calculatePrincipalAxes(a, md.atom_type, params);
            break;
        case number_of_atoms:
            calculateNatoms(a, md.atom_type, params);
            break;
        default:
            break;
        }
    }

    double compare(Descriptor* d){
        if(this->dimensionOfDescriptorSpace != d->dimensionOfDescriptorSpace){
            return -1;
        }
        double sum = 0;
        for(int i = 0; i < dimensionOfDescriptorSpace; i++){
            sum += (this->descs[i] - d->descs[i])*(this->descs[i] - d->descs[i]);
        }
        return sum;
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
    MolecularDatabaseDescriptor* usedDescriptors=nullptr;
    int nbOfusedDescriptors = 0;
public:    
// User can choose desriptors to use by writing it in the usedDescriptor array (chosen order is important):
// first parameter is the type of descriptor (principal_axes, number_of_atoms)
// second parameter is the atom type (-1 means all atoms)
    MolecularDatabase(){

        nMaxCell = 100;
        alloc(nMaxCell);
    };
    MolecularDatabase(int n){
        alloc(n);
    };
    ~MolecularDatabase(){
        delete[] atoms;
        delete[] descriptors;
        delete[] usedDescriptors;
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

    Atoms loadAtoms(int i){
        if(!atoms) return nullptr;
        return atoms[i];
    };


    void testHash(Atoms* a, MMFFparams* params=0){
        Descriptor d;
        for(int i = 0; i < nbOfusedDescriptors; i++){
            d.chooseDescriptor(usedDescriptors[i], a, params);
        }
        hashDescriptors(&d);
        if(1){
            addMember(a, &d);
        }

    };

    int hashDescriptors(Descriptor* d){
        return 0;
    };

    void setDescriptors( MMFFparams* params=0, Atoms* atoms=0, MolecularDatabaseDescriptor* md = 0, int nbOfUsedDescriptors_ = 0){
        if(nbOfUsedDescriptors_ == 0){
            nbOfusedDescriptors = 2;
            usedDescriptors = new MolecularDatabaseDescriptor[nbOfusedDescriptors];
            int ityp = atoms->atypes[1];
            usedDescriptors[0] = {principal_axes, -1}; 
            usedDescriptors[1] = {number_of_atoms, -1};
        }
        else{
            nbOfusedDescriptors = nbOfUsedDescriptors_;
            usedDescriptors = new MolecularDatabaseDescriptor[nbOfusedDescriptors];
            for(int i = 0; i < nbOfusedDescriptors; i++){
                usedDescriptors[i] = md[i];
            }
        }
    };

    void FindRotation( Mat3d& rot, Atoms* a){
        Mat3d XX=Mat3dZero;
        for(int i=0; i<a->natoms; i++){
            XX.addOuter( a->apos[i], a->apos[i], 1.0 );
        }
        //printf( "XX: " ); printMat(XX);
        Vec3d evs;
        XX.eigenvals(evs);
        evs.sort();        
        //printf(  "FindRotation evs(%g,%g,%g) \n", evs.x,evs.y,evs.z );

        XX.eigenvec( evs.x, rot.a );
        XX.eigenvec( evs.y, rot.b );
        XX.eigenvec( evs.z, rot.c );
    }

    void orient( Vec3d center, Vec3d dir, Vec3d up, Atoms* a ){
        Mat3d rot; rot.fromDirUp( dir, up );
        for(int i=0; i< a->natoms; i++){
            Vec3d p = a->apos[i]-center;
            rot.dot_to(p,a->apos[i]);
        }
    }

    double compareAtoms(Atoms* a, int i){
        if(!nMembers){
            printf("No members in the database\n");
            return -1;
        }

        double dist = 0;

        Mat3d rot = Mat3dIdentity;
        Mat3d rot2 = Mat3dIdentity;
        Mat3d rot_orig, rot2_orig;
        Vec3d cog, cog2;
        double* ws = new double[a->natoms];
        for(int j=0; j<a->natoms; j++){
            //ws[j] = a->atypes[j];
            ws[j] = 1.0;
        }

        getCOG(a->natoms, a->apos, ws, cog);

        orient(cog, rot.c, rot.b, a);
        FindRotation(rot, a);
  
        orient(Vec3dZero, rot.c, rot.b, a);

        getCOG(atoms[i].natoms, atoms[i].apos, ws, cog2);
        orient(cog2, rot2.c, rot2.b, &atoms[i]);
        FindRotation(rot2, &atoms[i]);
        orient(Vec3dZero, rot2.c, rot2.b, &atoms[i]);

        bool mirror_x = false;
        bool mirror_y = false;
        bool mirror_z = false;
        double threshold = 0.001;
        for (int k = 0; k < a->natoms; k++)
        {
            if ((a->apos[k].x) * (a->apos[k].x) > threshold)
            {
                if ((a->apos[k].x + atoms[i].apos[k].x) * (a->apos[k].x + atoms[i].apos[k].x) < 0.01)
                {
                    for (int j = 0; j < a->natoms; j++)
                        a->apos[j].x = -a->apos[j].x;
                    mirror_x = true;
                }
                break;
            }
        }
        for (int k = 0; k < a->natoms; k++)
        {
            if ((a->apos[k].y) * (a->apos[k].y) > threshold)
            {
                if ((a->apos[k].y + atoms[i].apos[k].y) * (a->apos[k].y + atoms[i].apos[k].y) < 0.01)
                {
                    for (int j = 0; j < a->natoms; j++)
                        a->apos[j].y = -a->apos[j].y;
                    mirror_y = true;
                }
                break;
            }
        }
        for (int k = 0; k < a->natoms; k++)
        {
            if ((a->apos[k].z) * (a->apos[k].z) > threshold)
            {
                if ((a->apos[k].z + atoms[i].apos[k].z) * (a->apos[k].z + atoms[i].apos[k].z) < 0.01)
                {
                    for (int j = 0; j < a->natoms; j++)
                        a->apos[j].z = -a->apos[j].z;
                    mirror_z = true;
                }
                break;
            }
        }

        AtomicConfiguration ac, ac2;
        ac.natoms = a->natoms;
        ac.types = a->atypes;
        ac.pos = a->apos;

        ac2.natoms = atoms[i].natoms;
        ac2.types = atoms[i].atypes;
        ac2.pos = atoms[i].apos;

        dist = ac.dist(ac2);

        rot2.fromDirUp(rot2.c, rot2.b);
        rot2 = rot2.transposed();
        for (int j = 0; j < atoms[i].natoms; j++)
        {
            Vec3d p = atoms[i].apos[j];
            rot2.dot_to(p, atoms[i].apos[j]);
            atoms[i].apos[j].add(cog2);
        }

        if (mirror_x)
            for (int j = 0; j < a->natoms; j++)
                a->apos[j].x = -a->apos[j].x;
        if (mirror_y)
            for (int j = 0; j < a->natoms; j++)
                a->apos[j].y = -a->apos[j].y;
        if (mirror_z)
            for (int j = 0; j < a->natoms; j++)
                a->apos[j].z = -a->apos[j].z;
        

        rot.fromDirUp(rot.c, rot.b);
        rot = rot.transposed();
        for (int j = 0; j < a->natoms; j++)
        {
            Vec3d p = a->apos[j];
            rot.dot_to(p, a->apos[j]);
            a->apos[j].add(cog);
        }
        delete[] ws;
        return dist;
    }

    double compareDescriptors(int i, int j)
    {
        return descriptors[i].compare(&descriptors[j]);
    };
    // virtual char* tostr(int id){};
};

#endif