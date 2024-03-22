#ifndef MolecularDatabase_h
#define MolecularDatabase_h

#include <cmath>

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "VecN.h"
#include "Mat3.h"
#include "quaternion.h"
#include "NBFF.h"
#include "Vec3Utils.h"
#include "MMFFparams.h"
#include "AtomicConfiguration.h"
#include "macroUtils.h"

#include "MMFFBuilder.h"

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
            if(verbosity) printf("Not enough atoms to calculate principal axes\n");
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



    void setDescriptor(MolecularDatabaseDescriptor md, Atoms* a, MMFFparams* params=0){
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
    inline int getNMembers(){   return nMembers;   };

    Atoms loadAtoms(int i){
        if(!atoms || i > nMembers) return nullptr;
        return atoms[i];
    };


    void setDescriptors( Atoms* atoms=0, MolecularDatabaseDescriptor* md = 0, int nbOfUsedDescriptors_ = 0){
        if(verbosity) printf("MolecularDatabase::setDescriptors\n");
        if(nbOfUsedDescriptors_ == 0){
            nbOfusedDescriptors = 2;
            usedDescriptors = new MolecularDatabaseDescriptor[nbOfusedDescriptors];
            int ityp0 = atoms->atypes[0];
            int ityp1 = atoms->atypes[1];
            int ityp2 = atoms->atypes[2];
            usedDescriptors[0] = {principal_axes, -1};
            usedDescriptors[1] = {number_of_atoms, -1};
            //usedDescriptors[2] = {principal_axes, ityp1};
            //usedDescriptors[3] = {principal_axes, ityp2};
        }
        else
        {
            nbOfusedDescriptors = nbOfUsedDescriptors_;
            usedDescriptors = new MolecularDatabaseDescriptor[nbOfusedDescriptors];
            for (int i = 0; i < nbOfusedDescriptors; i++)
            {
                usedDescriptors[i] = md[i];
            }
        }
    };

    void FindRotation(Mat3d &rot, Atoms *a)
    {
        Mat3d XX = Mat3dZero;
        for (int i = 0; i < a->natoms; i++)
        {
            XX.addOuter(a->apos[i], a->apos[i], 1.0);
        }
        // printf( "XX: " ); printMat(XX);
        Vec3d evs;
        XX.eigenvals(evs);

        // printf(  "FindRotation evs(%g,%g,%g) \n", evs.x,evs.y,evs.z );
        evs.sort();

        XX.eigenvec(evs.x, rot.a);
        XX.eigenvec(evs.y, rot.b);
        XX.eigenvec(evs.z, rot.c);
    }

    void orient(Vec3d center, Vec3d dir, Vec3d up, Atoms *a)
    {
        Mat3d rot;
        rot.fromDirUp(dir, up);
        for (int i = 0; i < a->natoms; i++)
        {
            Vec3d p = a->apos[i] - center;
            rot.dot_to(p, a->apos[i]);
        }
    }

    void replace(Atoms* a, int i, MMFFparams* params=0){
        printf("replace %d\n", i);
        if(i < nMembers && i >= 0){
            atoms[i]=*a;
            Descriptor* d=new Descriptor();
            for(int j = 0; j < nbOfusedDescriptors; j++){
                d->setDescriptor(usedDescriptors[j], a, params);
            }
            descriptors[i]=*d;
        }
        else{
            printf("Cannot replace: Index %d out of range\n", i);
        }
    }

    void remove(int i){
        if(i < nMembers && i >= 0){
            for(int j = i; j < nMembers-1; j++){
                atoms[j]=atoms[j+1];
                descriptors[j]=descriptors[j+1];
            }
            nMembers--;
        }
        else{
            printf("Cannot remove: Index %d out of range\n", i);
        }
    }

    // double compareAtoms(Atoms *atoms_h, Atoms *atoms_i)
    // {
    //     if (atoms_h == nullptr || atoms_i == nullptr)
    //     {
    //         printf("Null pointer passed\n");
    //         return -1;
    //     }

    //     double dist = 0;

    //     Mat3d rot = Mat3dIdentity;
    //     Mat3d rot2 = Mat3dIdentity;
    //     Vec3d cog, cog2;
    //     double *ws = new double[atoms_h->natoms];
    //     for (int j = 0; j < atoms_h->natoms; j++)
    //     {
    //         // ws[j] = atom_h->atypes[j];
    //         ws[j] = 1.0;
    //     }

    //     getCOG(atoms_h->natoms, atoms_h->apos, ws, cog);
    //     orient(cog, rot.c, rot.b, atoms_h);
    //     FindRotation(rot, atoms_h);

    //     getCOG(atoms_i->natoms, atoms_i->apos, ws, cog2);
    //     orient(cog2, rot2.c, rot2.b, atoms_i);
    //     FindRotation(rot2, atoms_i);

    //     orient(Vec3dZero, rot.c, rot.b, atoms_h);
    //     orient(Vec3dZero, rot2.c, rot2.b, atoms_i);

    //     Vec3i direction_h =Vec3iOne, direction_i =Vec3iOne;

    //     bool rot_x_h = false;
    //     bool rot_y_h = false;
    //     bool rot_z_h = false;
    //     bool rot_x_i = false;
    //     bool rot_y_i = false;
    //     bool rot_z_i = false;
    //     double threshold = 0.0;

    //     if ((atoms_h->apos[0].x) * (atoms_h->apos[0].x) > threshold && atoms_h->apos[0].x < 0)
    //     {
    //         direction_h.x = std::copysign(1.0, atoms_h->apos[0].x);
    //     }
    //     if ((atoms_h->apos[0].y) * (atoms_h->apos[0].y) > threshold && atoms_h->apos[0].y < 0)
    //     {
    //         direction_h.y = std::copysign(1.0, atoms_h->apos[0].y);
    //     }
    //     if ((atoms_h->apos[0].z) * (atoms_h->apos[0].z) > threshold && atoms_h->apos[0].z < 0)
    //     {
    //         direction_h.z = std::copysign(1.0, atoms_h->apos[0].z);
    //     }
    //     if ((atoms_i->apos[0].x) * (atoms_i->apos[0].x) > threshold && atoms_i->apos[0].x < 0)
    //     {
    //         direction_i.x = std::copysign(1.0, atoms_i->apos[0].x);
    //     }
    //     if ((atoms_i->apos[0].y) * (atoms_i->apos[0].y) > threshold && atoms_i->apos[0].y < 0)
    //     {
    //         direction_i.y = std::copysign(1.0, atoms_i->apos[0].y);
    //     }
    //     if ((atoms_i->apos[0].z) * (atoms_i->apos[0].z) > threshold && atoms_i->apos[0].z < 0)
    //     {
    //         direction_i.z = std::copysign(1.0, atoms_i->apos[0].z);
    //     }

    //     //printf("direction_h: %d %d %d\n", direction_h.x, direction_h.y, direction_h.z);
    //     //printf("direction_i: %d %d %d\n", direction_i.x, direction_i.y, direction_i.z);


    //     if (direction_h.totprod() * direction_i.totprod() > 0)
    //     {
    //         switch (direction_h.x)
    //         {
    //         case 1:
    //             switch (direction_h.y)
    //             {
    //             case 1:
    //                 switch (direction_h.z)
    //                 {
    //                 case 1:
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_h->natoms; j++)
    //                         atoms_h->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);
    //                     rot_y_h = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_h.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms_h->natoms; j++)
    //                         atoms_h->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
    //                     rot_z_h = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_h->natoms; j++)
    //                         atoms_h->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
    //                     rot_x_h = true;
    //                     break;
    //                 }
    //                 break;
    //             }
    //             break;
    //         case -1:
    //             switch (direction_h.y)
    //             {
    //             case 1:
    //                 switch (direction_h.z)
    //                 {
    //                 case 1:
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_h->natoms; j++)
    //                         atoms_h->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);
    //                     rot_y_h = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_h.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms_h->natoms; j++)
    //                         atoms_h->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
    //                     rot_z_h = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_h->natoms; j++)
    //                         atoms_h->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
    //                     rot_x_h = true;
    //                     break;
    //                 }
    //                 break;
    //             }
    //             break;
    //         }
    //         switch (direction_i.x)
    //         {
    //         case 1:
    //             switch (direction_i.y)
    //             {
    //             case 1:
    //                 switch (direction_i.z)
    //                 {
    //                 case 1:
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_i->natoms; j++)
    //                         atoms_i->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
    //                     rot_y_i = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_i.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms_i->natoms; j++)
    //                         atoms_i->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
    //                     rot_z_i = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_i->natoms; j++)
    //                         atoms_i->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
    //                     rot_x_i = true;
    //                     break;
    //                 }
    //                 break;
    //             }
    //             break;
    //         case -1:
    //             switch (direction_i.y)
    //             {
    //             case 1:
    //                 switch (direction_i.z)
    //                 {
    //                 case 1:
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_i->natoms; j++)
    //                         atoms_i->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
    //                     rot_y_i = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_i.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms_i->natoms; j++)
    //                         atoms_i->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
    //                     rot_z_i = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms_i->natoms; j++)
    //                         atoms_i->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
    //                     rot_x_i = true;
    //                     break;
    //                 }
    //                 break;
    //             }
    //             break;
    //         }
    //     }

    //     AtomicConfiguration ac, ac2;
    //     ac.natoms = atoms_h->natoms;
    //     ac.types = atoms_h->atypes;
    //     ac.pos = atoms_h->apos;

    //     ac2.natoms = atoms_i->natoms;
    //     ac2.types = atoms_i->atypes;
    //     ac2.pos = atoms_i->apos;

    //     dist = sqrt(ac.dist(ac2));
    //     // for (int j =0; j< atoms_i->natoms; j++)
    //     // {
    //     //     dist = dist + (atoms_i->apos[j]-atoms_h->apos[j]).norm2();
    //     // }
    //     dist = sqrt(dist);

    //     if (rot_z_i)
    //         for (int j = 0; j < atoms_i->natoms; j++)
    //             atoms_i->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
    //     if (rot_y_i)
    //         for (int j = 0; j < atoms_i->natoms; j++)
    //             atoms_i->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
    //     if (rot_x_i)
    //         for (int j = 0; j < atoms_i->natoms; j++)
    //             atoms_i->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
    //     if (rot_z_h)
    //         for (int j = 0; j < atoms_h->natoms; j++)
    //             atoms_h->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
    //     if (rot_y_h)
    //         for (int j = 0; j < atoms_h->natoms; j++)
    //             atoms_h->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
    //     if (rot_x_h)
    //         for (int j = 0; j < atoms_h->natoms; j++)
    //             atoms_h->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);

    //     rot2.fromDirUp(rot2.c, rot2.b);
    //     rot2 = rot2.transposed();
    //     for (int j = 0; j < atoms_i->natoms; j++)
    //     {
    //         Vec3d p = atoms_i->apos[j];
    //         rot2.dot_to(p, atoms_i->apos[j]);
    //         atoms_i->apos[j].add(cog2);
    //     }

    //     rot.fromDirUp(rot.c, rot.b);
    //     rot = rot.transposed();
    //     for (int j = 0; j < atoms_h->natoms; j++)
    //     {
    //         Vec3d p = atoms_h->apos[j];
    //         rot.dot_to(p, atoms_h->apos[j]);
    //         atoms_h->apos[j].add(cog);
    //     }

    //     delete[] ws;
    //     return dist;
    // }

    double compareAtoms(int h, int i)
    {
        return compareAtoms(&atoms[h], &atoms[i]);
    }

    double compareAtoms(Atoms *atoms_h, int i)
    {
        return compareAtoms(atoms_h, &atoms[i]);
    }

    double compareDescriptors(int i, int j)
    {
        return descriptors[i].compareDescriptors(&descriptors[j]);
    };

    double compareDescriptors(Atoms* atoms, int j)
    {
        Descriptor d;
        for(int i = 0; i < nbOfusedDescriptors; i++)
            d.setDescriptor(usedDescriptors[i], atoms);
        return d.compareDescriptors(&descriptors[j]);
    };

    void addMember(Atoms* structure, Descriptor* descriptor=0){
        if(nMembers >= nMaxCell){
            if(verbosity) printf("Reallocating to maximaly store %d molecules\n", nMaxCell*2);
            realloc(nMaxCell*2);
        }

        if(!descriptor){
            descriptor = new Descriptor();
            for(int i = 0; i < nbOfusedDescriptors; i++)
                descriptor->setDescriptor(usedDescriptors[i], structure);
        }
        
        atoms[nMembers].copyOf(*structure);
        descriptors[nMembers].copyOf(*descriptor);
        if(nMembers == 0)   dimensionOfDescriptorSpace = descriptors[nMembers].dimensionOfDescriptorSpace;
        this->nMembers++;
    };

    int addIfNewDescriptor(Atoms* a, MMFFparams* params=0){
        Descriptor d;
        int sameDescriptor = -1;
        //if(verbosity) printf("testHash\n");
        for(int i = 0; i < nbOfusedDescriptors; i++){
            d.setDescriptor(usedDescriptors[i], a, params);
        }
        for(int i = 0; i < nMembers; i++){
            if(d.compareDescriptors(&descriptors[i]) < 0.2){
                //addMember(a, &d);
                printf("compareAtoms: %lf\n", compareAtoms(a, &atoms[i]) ); 
                if (compareAtoms(a, &atoms[i]) > 1.0)
                {
                    nMembers--;
                }
                else
                {
                    sameDescriptor = i;
                    break;
                }
            }
        }
        if(sameDescriptor == -1){
            addMember(a, &d);
        }
        //hashDescriptors(&d);

        return sameDescriptor;
    };

    int hashDescriptors(Descriptor* d){
        return 0;
    };
    // virtual char* tostr(int id){};












































    double compareAtoms(Atoms *atoms_h, Atoms *atoms_i)
    {
        if (atoms_h == nullptr || atoms_i == nullptr)
        {
            printf("Null pointer passed\n");
            return -1;
        }

        int natm = atoms_h->natoms;

        double dist = 0;

        printf("Comparing atoms\n");

        //---Center on origin---//
        Vec3d x0_center = Vec3dZero;
        Vec3d x1_center = Vec3dZero;

        double *ws = new double[natm];
        for (int j = 0; j < natm; j++)
        {
            // ws[j] = coord0_->atypes[j];
            ws[j] = 1.0;
        }

        getCOG(natm, atoms_h->apos, ws, x0_center);
        getCOG(natm, atoms_i->apos, ws, x1_center);

        for (int j = 0; j < natm; j++)
        {
            atoms_h->apos[j] = atoms_h->apos[j] - x0_center;
            atoms_i->apos[j] = atoms_i->apos[j] - x1_center;
        }

        //---End Center on origin---//


        Mat3d rot;
        rot = superimposer(atoms_h, atoms_i, atoms_h->natoms);

        orient(Vec3dZero, rot.c, rot.b, atoms_h);


        AtomicConfiguration ac, ac2;
        ac.natoms = atoms_h->natoms;
        ac.types = atoms_h->atypes;
        ac.pos = atoms_h->apos;

        ac2.natoms = atoms_i->natoms;
        ac2.types = atoms_i->atypes;
        ac2.pos = atoms_i->apos;

        dist = sqrt(ac.dist(ac2));

        rot.makeT();
        orient(Vec3dZero, rot.c, rot.b, atoms_h);
        for(int j=0 ; j<natm; j++){
            atoms_h->apos[j] = atoms_h->apos[j] + x0_center;
            atoms_i->apos[j] = atoms_i->apos[j] + x1_center;
        }

 
        return dist;
    }


/*
Transformed code from https://github.com/willhooper/superpose/tree/master
Works similar as Kabsch algorithm
Aligns two sets of coordinates by symetrizing matrix from Kabsch algorithm
Symmetric matrix is already diagonal, but in diferent orthogonal basis
*/
    Mat3d superimposer(Atoms* coord0_, Atoms* coord1_, unsigned int natm)
    {
        Vec3d* coord0 = coord0_->apos;
        Vec3d* coord1 = coord1_->apos;

        Mat3d mtx = Mat3dZero;
        float tolerance = 0.0001;

        // Error checking
        unsigned int MAXATOM = 200;
        if (natm > MAXATOM)
        {
            printf( "Error! Increase maxatom\n");
            throw 0;
        }
        if (natm < 4)
        {
            printf( "Error! Too few atoms (must be 4 or greater)\n");
            throw 0;
        }

        double** x0_ = new double*[3];
        double** x1_ = new double*[3];
        for (int i = 0; i < 3; i++)
        {
            x0_[i] = new double[natm];
            x1_[i] = new double[natm];
        }
        for (unsigned int i = 0; i < natm; ++i)
        {
            x0_[0][i] = coord0[i].x;
            x1_[0][i] = coord1[i].x;
            x0_[1][i] = coord0[i].y;
            x1_[1][i] = coord1[i].y;
            x0_[2][i] = coord0[i].z;
            x1_[2][i] = coord1[i].z;
        }

        Mat3d aa = Mat3dZero;
        aa.xx = VecN::dot(natm, x0_[0], x1_[0]);
        aa.xy = VecN::dot(natm, x0_[0], x1_[1]);
        aa.xz = VecN::dot(natm, x0_[0], x1_[2]);
        aa.yx = VecN::dot(natm, x0_[1], x1_[0]);
        aa.yy = VecN::dot(natm, x0_[1], x1_[1]);
        aa.yz = VecN::dot(natm, x0_[1], x1_[2]);
        aa.zx = VecN::dot(natm, x0_[2], x1_[0]);
        aa.zy = VecN::dot(natm, x0_[2], x1_[1]);
        aa.zz = VecN::dot(natm, x0_[2], x1_[2]);



        // Initialize rotation matrix
        Mat3d rotation = Mat3dIdentity;

        //---Iterative rotation scheme---//
        int iteration_count = 0;
        bool do51 = true; //"This is a way to deal with those nasty gotos in the FORTRAN code"
        int iflag;
        int ix;
        double sigma;
        double gamma;
        double sg;
        Vec3d bb = Vec3dZero;
        Vec3d cc = Vec3dZero;
        Vec3i ind = {0,1,2};
        while (true)
        {
            if (do51)
            {
                iflag = 0;
                ix = 0;
            }

            // If the number of iterations exceeds 500, give up
            ++iteration_count;
            if (iteration_count > 100)
            {
                break;
            }

            aa.swap_vecs(ind);
            rotation.swap_vecs(ind);
            ind = {1, 2, 0};

            switch (ix)
            {
            case 0:
                sigma = aa.c.y - aa.b.z;
                gamma = aa.b.y + aa.c.z;
                break;
            case 1:
                sigma = aa.c.z - aa.b.x;
                gamma = aa.b.z + aa.c.x;
                break;
            case 2:
                sigma = aa.c.x - aa.b.y;
                gamma = aa.b.x + aa.c.y;
                break;

            default:
                break;
            }

            sg = sqrt(sigma * sigma + gamma * gamma);

            if (sg == 0)
            {
                ++ix;
                if (iflag == 0)
                {
                    break;
                }
                if (ix < 3)
                {
                    do51 = false;
                }
                else
                {
                    do51 = true;
                }
                continue;
            }

            sg = 1.0 / sg;
            if (fabs(sigma) < (tolerance * fabs(gamma)))
            {
                ++ix;
                if (iflag == 0)
                {
                    break;
                }
                if (ix < 3)
                {
                    do51 = false;
                }
                else
                {
                    do51 = true;
                }

                continue;
            }



            bb.set_add_mul(Vec3dZero, aa.b, gamma);
            bb.add_mul(aa.c, sigma);
            bb.mul(sg);
            cc.set_add_mul(Vec3dZero, aa.c, gamma);
            cc.add_mul(aa.b, 0-sigma);
            cc.mul(sg);
            aa.b = bb;
            aa.c = cc;

            bb.set_add_mul(Vec3dZero, rotation.b, gamma);
            bb.add_mul(rotation.c, sigma);
            bb.mul(sg);
            cc.set_add_mul(Vec3dZero, rotation.c, gamma);
            cc.add_mul(rotation.b, -sigma);
            cc.mul(sg);
            rotation.b = bb;
            rotation.c = cc;


            iflag = 1;

            ++ix;
            if (iflag == 0)
            {
                break;
            }
            if (ix < 3)
            {
                do51 = false;
            }
            else
            {
                do51 = true;
            }

            continue;

        } // End while loop


        return rotation;
    }
















    void as_rigid_as_possible(Atoms* atoms_h_, int index, int* neighs=0, int nNeighs=4){
        printf("as_rigid_as_possible::Aligning new structure %d to atoms[%d]\n", nMembers-1, index);
        if(!neighs){
            printf("No neighbours given\n");
            return;
        }
        Atoms* atoms_h = atoms_h_;
        //atoms_h.copyOf(atoms_h_);
        Atoms* atoms_i = &atoms[index];
        int natm = atoms_h->natoms;
        double E = 0;
        Mat3d S = Mat3dZero;
        Mat3d* R = new Mat3d[natm];
        Mat3d U, V;
        Vec3d* F = new Vec3d[natm];
        Vec3d temp;
        Quat4d* w_mat = new Quat4d[natm];
        double* w_vec = new double[natm];
        double alpha = 0.01;
        for (int j = 0; j < natm; j++){R[j]=Mat3dZero; w_mat[j] = Quat4dOnes; w_vec[j] = 1;
            // for(int k = 0; k < nNeighs; k++){
            //     if(neighs[nNeighs*j+k] == -1) break;
            //     w_mat[j].array[k] = atoms_h;
            // }
        
        }    //ToDo: set weights

        const int qmax =1;
        int q = 0;
        while(qmax > q){
            E = 0;
            for(int j = 0; j < natm; j++){
                S = Mat3dZero;
                bool capping = false;
                for(int k = 0; k < nNeighs; k++){
                    if(neighs[nNeighs*j+1] == -1){ capping = true; break;}
                    if(neighs[nNeighs*j+k] == -1) break;
                    Vec3d e_h = atoms_h->apos[j] - atoms_h->apos[neighs[nNeighs*j+k]];
                    Vec3d e_i = atoms_i->apos[j] - atoms_i->apos[neighs[nNeighs*j+k]];
                    // printf("e_h: %lf %lf %lf\n", e_h.x, e_h.y, e_h.z);
                    // printf("e_i: %lf %lf %lf\n", e_i.x, e_i.y, e_i.z);
                    S.addOuter(e_h, e_i, w_mat[j].array[k]);
                }
                //S.print();
                if(capping)continue;


                S.SVD2(U, temp, V);

                R[j].set_mmul_NT(V, U);

                double E_j = 0;
                for(int k = 0; k < nNeighs; k++){
                    Vec3d e_h = atoms_h->apos[j] - atoms_h->apos[neighs[nNeighs * j + k]];
                    Vec3d e_i = atoms_i->apos[j] - atoms_i->apos[neighs[nNeighs * j + k]];
                    double d = (e_h - R[j].dot(e_i)).norm2();
                    E_j += w_mat[j].array[k]*d;
                }
                E += w_vec[j]*E_j;
            }
            for(int j = 0; j < natm; j++){
                Mat3d Rot;
                F[j] = Vec3dZero;
                for(int k = 0; k < nNeighs; k++){
                    if(neighs[nNeighs*j+1] == -1){break;}
                    if(neighs[nNeighs*j+k] == -1) break;
                    Vec3d e_h = atoms_h->apos[j] - atoms_h->apos[neighs[nNeighs*j+k]];
                    Vec3d e_i = atoms_i->apos[j] - atoms_i->apos[neighs[nNeighs * j + k]];
                    Rot = Mat3dZero;
                    Rot.add(R[j]);
                    Rot.add(R[neighs[nNeighs * j + k]]);
                    F[j] = F[j] + (Rot.dot(e_h).mul(-0.5) + e_i).mul(4 * w_mat[j].array[k]);
                    
                }
                printf("%d F: %lf %lf %lf\n",j, F[j].x, F[j].y, F[j].z);
            }
            for(int j = 0; j < natm; j++){
                if(neighs[nNeighs*j+1] == -1){ atoms_h->apos[j] = atoms_h->apos[j] - F[neighs[nNeighs*j]].mul(alpha);}
                else{ atoms_h->apos[j] = atoms_h->apos[j] - F[j].mul(alpha); }
            }
            printf("%d: %lf\n", q, E);
            q++;
            atoms_h->print();
        }
        delete[] R;
        delete[] w_mat;
        delete[] w_vec;
        delete[] F;
        //printf("S: %lf %lf %lf\n", S.x, S.y, S.z);
    }
};

#endif