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
#include "macroUtils.h"

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
        Mat3d I=Mat3dZero;
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
        descs[dimensionOfDescriptorSpace-3] = sqrt(e.x);
        descs[dimensionOfDescriptorSpace-2] = sqrt(e.y);
        descs[dimensionOfDescriptorSpace-1] = sqrt(e.z);
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

    double compareDescriptors(Descriptor* d){
        if(this->dimensionOfDescriptorSpace != d->dimensionOfDescriptorSpace){
        //     printf("Descriptors have different dimensions\n");
        //     printf("This: %d, Other: %d\n", this->dimensionOfDescriptorSpace, d->dimensionOfDescriptorSpace);
            return -1;
        }
        double sum = 0;
        for(int i = 0; i < dimensionOfDescriptorSpace; i++){
            sum += (this->descs[i] - d->descs[i])*(this->descs[i] - d->descs[i]);
        }
        return sqrt(sum);
    }

    void copyOf(const Descriptor& d){
        if(this->dimensionOfDescriptorSpace != d.dimensionOfDescriptorSpace){
            alloc(d.dimensionOfDescriptorSpace);
        }
        memcpy(this->descs, d.descs, sizeof(double)*d.dimensionOfDescriptorSpace);
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
        
        Atoms* temp = new Atoms[n];
        for(int i = 0; i < nMembers; i++){
            temp[i].copyOf(this->atoms[i]);
        }
        delete[] this->atoms;
        this->atoms = temp;



        Descriptor *temp2 = new Descriptor[n];

        for(int i = 0; i < nMembers; i++){
            temp2[i].copyOf(this->descriptors[i]);
        }

        delete[] this->descriptors;
        this->descriptors = temp2;

        this->nMaxCell = n;

    };
    void addMember(Atoms* structure, Descriptor* descriptor){
        if(nMembers >= nMaxCell){
            printf("Reallocating to maximaly store %d molecules\n", nMaxCell*2);
            realloc(nMaxCell*2);
        }
        
        atoms[nMembers].copyOf(*structure);
        descriptors[nMembers].copyOf(*descriptor);
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
        if(!atoms || i > nMembers) return nullptr;
        return atoms[i];
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
                
        //printf(  "FindRotation evs(%g,%g,%g) \n", evs.x,evs.y,evs.z );
        evs.sort();
        
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

    double compareAtoms(int h, int i){
        if(!nMembers){
            printf("No members in the database\n");
            return -1;
        }else if (h >= nMembers || i >= nMembers){
            printf("Index out of bounds\n");
            return -1;
        }
        
        double dist = 0;

        Mat3d rot = Mat3dIdentity;
        Mat3d rot2 = Mat3dIdentity;
        Vec3d cog, cog2;
        double* ws = new double[atoms[h].natoms];
        for(int j=0; j<atoms[h].natoms; j++){
            //ws[j] = atoms[h].atypes[j];
            ws[j] = 1.0;
        }

        getCOG(atoms[h].natoms, atoms[h].apos, ws, cog);
        orient(cog, rot.c, rot.b, &atoms[h]);
        FindRotation(rot, &atoms[h]);

        getCOG(atoms[i].natoms, atoms[i].apos, ws, cog2);
        orient(cog2, rot2.c, rot2.b, &atoms[i]);
        FindRotation(rot2, &atoms[i]);
        
        
        orient(Vec3dZero, rot.c, rot.b, &atoms[h]);
        orient(Vec3dZero, rot2.c, rot2.b, &atoms[i]);

        bool mirror_x_h = false;
        bool mirror_y_h = false;
        bool mirror_z_h = false;
        bool mirror_x_i = false;
        bool mirror_y_i = false;
        bool mirror_z_i = false;
        double threshold = 0.001;

        if((atoms[h].apos[0].x) * (atoms[h].apos[0].x) > threshold && atoms[h].apos[0].x<0){
                for (int j = 0; j < atoms[h].natoms; j++)
                    atoms[h].apos[j].x = -atoms[h].apos[j].x;
                mirror_x_h = true;
        }
        if((atoms[h].apos[0].y) * (atoms[h].apos[0].y) > threshold && atoms[h].apos[0].y<0){
                for (int j = 0; j < atoms[h].natoms; j++)
                    atoms[h].apos[j].y = -atoms[h].apos[j].y;
                mirror_y_h = true;
        }
        if((atoms[h].apos[0].z) * (atoms[h].apos[0].z) > threshold && atoms[h].apos[0].z<0){
                for (int j = 0; j < atoms[h].natoms; j++)
                    atoms[h].apos[j].z = -atoms[h].apos[j].z;
                mirror_z_h = true;
        }

        if((atoms[i].apos[0].x) * (atoms[i].apos[0].x) > threshold && atoms[i].apos[0].x<0){
                for (int j = 0; j < atoms[i].natoms; j++)
                    atoms[i].apos[j].x = -atoms[i].apos[j].x;
                mirror_x_i = true;
        }
        if((atoms[i].apos[0].y) * (atoms[i].apos[0].y) > threshold  && atoms[i].apos[0].y<0){
                for (int j = 0; j < atoms[i].natoms; j++)
                    atoms[i].apos[j].y = -atoms[i].apos[j].y;
                mirror_y_i = true;
        }
        if((atoms[i].apos[0].z) * (atoms[i].apos[0].z) > threshold && atoms[i].apos[0].z<0){
                for (int j = 0; j < atoms[i].natoms; j++)
                    atoms[i].apos[j].z = -atoms[i].apos[j].z;
                mirror_z_i = true;
        }

        AtomicConfiguration ac, ac2;
        ac.natoms = atoms[h].natoms;
        ac.types = atoms[h].atypes;
        ac.pos = atoms[h].apos;

        ac2.natoms = atoms[i].natoms;
        ac2.types = atoms[i].atypes;
        ac2.pos = atoms[i].apos;

        dist = ac.dist(ac2);




        if (mirror_x_h)
            for (int j = 0; j < atoms[h].natoms; j++)
                atoms[h].apos[j].x = -atoms[h].apos[j].x;
        if (mirror_y_h)
            for (int j = 0; j < atoms[h].natoms; j++)
                atoms[h].apos[j].y = -atoms[h].apos[j].y;
        if (mirror_z_h)
            for (int j = 0; j < atoms[h].natoms; j++)
                atoms[h].apos[j].z = -atoms[h].apos[j].z;
        if (mirror_x_i)
            for (int j = 0; j < atoms[i].natoms; j++)
                atoms[i].apos[j].x = -atoms[i].apos[j].x;
        if (mirror_y_i)
            for (int j = 0; j < atoms[i].natoms; j++)
                atoms[i].apos[j].y = -atoms[i].apos[j].y;
        if (mirror_z_i)
            for (int j = 0; j < atoms[i].natoms; j++)
                atoms[i].apos[j].z = -atoms[i].apos[j].z;
        
        rot2.fromDirUp(rot2.c, rot2.b);
        rot2 = rot2.transposed();
        for (int j = 0; j < atoms[i].natoms; j++)
        {
            Vec3d p = atoms[i].apos[j];
            rot2.dot_to(p, atoms[i].apos[j]);
            atoms[i].apos[j].add(cog2);
        }

        rot.fromDirUp(rot.c, rot.b);
        rot = rot.transposed();
        for (int j = 0; j < atoms[h].natoms; j++)
        {
            Vec3d p = atoms[h].apos[j];
            rot.dot_to(p, atoms[h].apos[j]);
            atoms[h].apos[j].add(cog);
        }
        
        delete[] ws;
        return dist;
    }

    double compareDescriptors(int i, int j)
    {
        return descriptors[i].compareDescriptors(&descriptors[j]);
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
    // virtual char* tostr(int id){};
};

#endif