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

    void addMember(Atoms* structure, Descriptor* descriptor){
        if(nMembers >= nMaxCell){
            if(verbosity) printf("Reallocating to maximaly store %d molecules\n", nMaxCell*2);
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
    inline int getNMembers(){   return nMembers;   };

    Atoms loadAtoms(int i){
        if(!atoms || i > nMembers) return nullptr;
        return atoms[i];
    };


    void setDescriptors( Atoms* atoms=0, MolecularDatabaseDescriptor* md = 0, int nbOfUsedDescriptors_ = 0){
        if(verbosity) printf("MolecularDatabase::setDescriptors\n");
        if(nbOfUsedDescriptors_ == 0){
            nbOfusedDescriptors = 4;
            usedDescriptors = new MolecularDatabaseDescriptor[nbOfusedDescriptors];
            int ityp0 = atoms->atypes[0];
            int ityp1 = atoms->atypes[1];
            int ityp2 = atoms->atypes[2];
            usedDescriptors[0] = {principal_axes, -1};
            usedDescriptors[1] = {principal_axes, ityp0};
            usedDescriptors[2] = {principal_axes, ityp1};
            usedDescriptors[3] = {principal_axes, ityp2};
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

    // double compareAtoms(int h, int i)
    // {
    //     if (!nMembers)
    //     {
    //         printf("No members in the database\n");
    //         return -1;
    //     }
    //     else if (h >= nMembers || i >= nMembers)
    //     {
    //         printf("Index out of bounds\n");
    //         return -1;
    //     }

    //     double dist = 0;

    //     Mat3d rot = Mat3dIdentity;
    //     Mat3d rot2 = Mat3dIdentity;
    //     Vec3d cog, cog2;
    //     double *ws = new double[atoms[h].natoms];
    //     for (int j = 0; j < atoms[h].natoms; j++)
    //     {
    //         // ws[j] = atoms[h].atypes[j];
    //         ws[j] = 1.0;
    //     }

    //     getCOG(atoms[h].natoms, atoms[h].apos, ws, cog);
    //     orient(cog, rot.c, rot.b, &atoms[h]);
    //     FindRotation(rot, &atoms[h]);

    //     getCOG(atoms[i].natoms, atoms[i].apos, ws, cog2);
    //     orient(cog2, rot2.c, rot2.b, &atoms[i]);
    //     FindRotation(rot2, &atoms[i]);

    //     orient(Vec3dZero, rot.c, rot.b, &atoms[h]);
    //     orient(Vec3dZero, rot2.c, rot2.b, &atoms[i]);

    //     Vec3i direction_h, direction_i;

    //     bool rot_x_h = false;
    //     bool rot_y_h = false;
    //     bool rot_z_h = false;
    //     bool rot_x_i = false;
    //     bool rot_y_i = false;
    //     bool rot_z_i = false;
    //     double threshold = 0.001;

    //     if ((atoms[h].apos[0].x) * (atoms[h].apos[0].x) > threshold && atoms[h].apos[0].x < 0)
    //     {
    //         direction_h.x = std::copysign(1.0, atoms[h].apos[0].x);
    //     }
    //     if ((atoms[h].apos[0].y) * (atoms[h].apos[0].y) > threshold && atoms[h].apos[0].y < 0)
    //     {
    //         direction_h.y = std::copysign(1.0, atoms[h].apos[0].y);
    //     }
    //     if ((atoms[h].apos[0].z) * (atoms[h].apos[0].z) > threshold && atoms[h].apos[0].z < 0)
    //     {
    //         direction_h.z = std::copysign(1.0, atoms[h].apos[0].z);
    //     }
    //     if ((atoms[i].apos[0].x) * (atoms[i].apos[0].x) > threshold && atoms[i].apos[0].x < 0)
    //     {
    //         direction_i.x = std::copysign(1.0, atoms[i].apos[0].x);
    //     }
    //     if ((atoms[i].apos[0].y) * (atoms[i].apos[0].y) > threshold && atoms[i].apos[0].y < 0)
    //     {
    //         direction_i.y = std::copysign(1.0, atoms[i].apos[0].y);
    //     }
    //     if ((atoms[i].apos[0].z) * (atoms[i].apos[0].z) > threshold && atoms[i].apos[0].z < 0)
    //     {
    //         direction_i.z = std::copysign(1.0, atoms[i].apos[0].z);
    //     }

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
    //                     for (int j = 0; j < atoms[h].natoms; j++)
    //                         atoms[h].apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);
    //                     rot_y_h = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_h.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms[h].natoms; j++)
    //                         atoms[h].apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
    //                     rot_z_h = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms[h].natoms; j++)
    //                         atoms[h].apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
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
    //                     for (int j = 0; j < atoms[h].natoms; j++)
    //                         atoms[h].apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);
    //                     rot_y_h = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_h.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms[h].natoms; j++)
    //                         atoms[h].apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
    //                     rot_z_h = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms[h].natoms; j++)
    //                         atoms[h].apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
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
    //                     for (int j = 0; j < atoms[i].natoms; j++)
    //                         atoms[i].apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
    //                     rot_y_i = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_i.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms[i].natoms; j++)
    //                         atoms[i].apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
    //                     rot_z_i = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms[i].natoms; j++)
    //                         atoms[i].apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
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
    //                     for (int j = 0; j < atoms[i].natoms; j++)
    //                         atoms[i].apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
    //                     rot_y_i = true;
    //                     break;
    //                 }
    //                 break;
    //             case -1:
    //                 switch (direction_i.z)
    //                 {
    //                 case 1:
    //                     for (int j = 0; j < atoms[i].natoms; j++)
    //                         atoms[i].apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
    //                     rot_z_i = true;
    //                     break;
    //                 case -1:
    //                     for (int j = 0; j < atoms[i].natoms; j++)
    //                         atoms[i].apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
    //                     rot_x_i = true;
    //                     break;
    //                 }
    //                 break;
    //             }
    //             break;
    //         }
    //     }

    //     AtomicConfiguration ac, ac2;
    //     ac.natoms = atoms[h].natoms;
    //     ac.types = atoms[h].atypes;
    //     ac.pos = atoms[h].apos;

    //     ac2.natoms = atoms[i].natoms;
    //     ac2.types = atoms[i].atypes;
    //     ac2.pos = atoms[i].apos;

    //     dist = ac.dist(ac2);


    //     if(rot_z_i)
    //         for (int j = 0; j < atoms[i].natoms; j++)
    //             atoms[i].apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
    //     if(rot_y_i)
    //         for (int j = 0; j < atoms[i].natoms; j++)
    //             atoms[i].apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
    //     if(rot_x_i)
    //         for (int j = 0; j < atoms[i].natoms; j++)
    //             atoms[i].apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
    //     if(rot_z_h)
    //         for (int j = 0; j < atoms[h].natoms; j++)
    //             atoms[h].apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
    //     if(rot_y_h)
    //         for (int j = 0; j < atoms[h].natoms; j++)
    //             atoms[h].apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
    //     if(rot_x_h)
    //         for (int j = 0; j < atoms[h].natoms; j++)
    //             atoms[h].apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);

        
    //     rot2.fromDirUp(rot2.c, rot2.b);
    //     rot2 = rot2.transposed();
    //     for (int j = 0; j < atoms[i].natoms; j++)
    //     {
    //         Vec3d p = atoms[i].apos[j];
    //         rot2.dot_to(p, atoms[i].apos[j]);
    //         atoms[i].apos[j].add(cog2);
    //     }

    //     rot.fromDirUp(rot.c, rot.b);
    //     rot = rot.transposed();
    //     for (int j = 0; j < atoms[h].natoms; j++)
    //     {
    //         Vec3d p = atoms[h].apos[j];
    //         rot.dot_to(p, atoms[h].apos[j]);
    //         atoms[h].apos[j].add(cog);
    //     }
        
    //     delete[] ws;
    //     return dist;
    // }


double compareAtoms(Atoms* atoms_h, Atoms* atoms_i)
{
    if (atoms_h == nullptr || atoms_i == nullptr)
    {
        printf("Null pointer passed\n");
        return -1;
    }

    double dist = 0;

    Mat3d rot = Mat3dIdentity;
    Mat3d rot2 = Mat3dIdentity;
    Vec3d cog, cog2;
    double *ws = new double[atoms_h->natoms];
    for (int j = 0; j < atoms_h->natoms; j++)
    {
        // ws[j] = atom_h->atypes[j];
        ws[j] = 1.0;
    }

    getCOG(atoms_h->natoms, atoms_h->apos, ws, cog);
    orient(cog, rot.c, rot.b, atoms_h);
    FindRotation(rot, atoms_h);

    getCOG(atoms_i->natoms, atoms_i->apos, ws, cog2);
    orient(cog2, rot2.c, rot2.b, atoms_i);
    FindRotation(rot2, atoms_i);

    orient(Vec3dZero, rot.c, rot.b, atoms_h);
    orient(Vec3dZero, rot2.c, rot2.b, atoms_i);

    Vec3i direction_h, direction_i;

    bool rot_x_h = false;
    bool rot_y_h = false;
    bool rot_z_h = false;
    bool rot_x_i = false;
    bool rot_y_i = false;
    bool rot_z_i = false;
    double threshold = 0.001;

        if ((atoms_h->apos[0].x) * (atoms_h->apos[0].x) > threshold && atoms_h->apos[0].x < 0)
        {
            direction_h.x = std::copysign(1.0, atoms_h->apos[0].x);
        }
        if ((atoms_h->apos[0].y) * (atoms_h->apos[0].y) > threshold && atoms_h->apos[0].y < 0)
        {
            direction_h.y = std::copysign(1.0, atoms_h->apos[0].y);
        }
        if ((atoms_h->apos[0].z) * (atoms_h->apos[0].z) > threshold && atoms_h->apos[0].z < 0)
        {
            direction_h.z = std::copysign(1.0, atoms_h->apos[0].z);
        }
        if ((atoms_i->apos[0].x) * (atoms_i->apos[0].x) > threshold && atoms_i->apos[0].x < 0)
        {
            direction_i.x = std::copysign(1.0, atoms_i->apos[0].x);
        }
        if ((atoms_i->apos[0].y) * (atoms_i->apos[0].y) > threshold && atoms_i->apos[0].y < 0)
        {
            direction_i.y = std::copysign(1.0, atoms_i->apos[0].y);
        }
        if ((atoms_i->apos[0].z) * (atoms_i->apos[0].z) > threshold && atoms_i->apos[0].z < 0)
        {
            direction_i.z = std::copysign(1.0, atoms_i->apos[0].z);
        }

        if (direction_h.totprod() * direction_i.totprod() > 0)
        {
            switch (direction_h.x)
            {
            case 1:
                switch (direction_h.y)
                {
                case 1:
                    switch (direction_h.z)
                    {
                    case 1:
                        break;
                    case -1:
                        for (int j = 0; j < atoms_h->natoms; j++)
                            atoms_h->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);
                        rot_y_h = true;
                        break;
                    }
                    break;
                case -1:
                    switch (direction_h.z)
                    {
                    case 1:
                        for (int j = 0; j < atoms_h->natoms; j++)
                            atoms_h->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
                        rot_z_h = true;
                        break;
                    case -1:
                        for (int j = 0; j < atoms_h->natoms; j++)
                            atoms_h->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
                        rot_x_h = true;
                        break;
                    }
                    break;
                }
                break;
            case -1:
                switch (direction_h.y)
                {
                case 1:
                    switch (direction_h.z)
                    {
                    case 1:
                        break;
                    case -1:
                        for (int j = 0; j < atoms_h->natoms; j++)
                            atoms_h->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);
                        rot_y_h = true;
                        break;
                    }
                    break;
                case -1:
                    switch (direction_h.z)
                    {
                    case 1:
                        for (int j = 0; j < atoms_h->natoms; j++)
                            atoms_h->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
                        rot_z_h = true;
                        break;
                    case -1:
                        for (int j = 0; j < atoms_h->natoms; j++)
                            atoms_h->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
                        rot_x_h = true;
                        break;
                    }
                    break;
                }
                break;
            }
            switch (direction_i.x)
            {
            case 1:
                switch (direction_i.y)
                {
                case 1:
                    switch (direction_i.z)
                    {
                    case 1:
                        break;
                    case -1:
                        for (int j = 0; j < atoms_i->natoms; j++)
                            atoms_i->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
                        rot_y_i = true;
                        break;
                    }
                    break;
                case -1:
                    switch (direction_i.z)
                    {
                    case 1:
                        for (int j = 0; j < atoms_i->natoms; j++)
                            atoms_i->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
                        rot_z_i = true;
                        break;
                    case -1:
                        for (int j = 0; j < atoms_i->natoms; j++)
                            atoms_i->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
                        rot_x_i = true;
                        break;
                    }
                    break;
                }
                break;
            case -1:
                switch (direction_i.y)
                {
                case 1:
                    switch (direction_i.z)
                    {
                    case 1:
                        break;
                    case -1:
                        for (int j = 0; j < atoms_i->natoms; j++)
                            atoms_i->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
                        rot_y_i = true;
                        break;
                    }
                    break;
                case -1:
                    switch (direction_i.z)
                    {
                    case 1:
                        for (int j = 0; j < atoms_i->natoms; j++)
                            atoms_i->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
                        rot_z_i = true;
                        break;
                    case -1:
                        for (int j = 0; j < atoms_i->natoms; j++)
                            atoms_i->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
                        rot_x_i = true;
                        break;
                    }
                    break;
                }
                break;
            }
        }

        AtomicConfiguration ac, ac2;
        ac.natoms = atoms_h->natoms;
        ac.types = atoms_h->atypes;
        ac.pos = atoms_h->apos;

        ac2.natoms = atoms_i->natoms;
        ac2.types = atoms_i->atypes;
        ac2.pos = atoms_i->apos;

        dist = ac.dist(ac2);


        if(rot_z_i)
            for (int j = 0; j < atoms_i->natoms; j++)
                atoms_i->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog2);
        if(rot_y_i)
            for (int j = 0; j < atoms_i->natoms; j++)
                atoms_i->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog2);
        if(rot_x_i)
            for (int j = 0; j < atoms_i->natoms; j++)
                atoms_i->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog2);
        if(rot_z_h)
            for (int j = 0; j < atoms_h->natoms; j++)
                atoms_h->apos[j].rotate(M_PI, Vec3d{0, 0, 1}, cog);
        if(rot_y_h)
            for (int j = 0; j < atoms_h->natoms; j++)
                atoms_h->apos[j].rotate(M_PI, Vec3d{1, 0, 0}, cog);
        if(rot_x_h)
            for (int j = 0; j < atoms_h->natoms; j++)
                atoms_h->apos[j].rotate(M_PI, Vec3d{0, 1, 0}, cog);

        
        rot2.fromDirUp(rot2.c, rot2.b);
        rot2 = rot2.transposed();
        for (int j = 0; j < atoms_i->natoms; j++)
        {
            Vec3d p = atoms_i->apos[j];
            rot2.dot_to(p, atoms_i->apos[j]);
            atoms_i->apos[j].add(cog2);
        }

        rot.fromDirUp(rot.c, rot.b);
        rot = rot.transposed();
        for (int j = 0; j < atoms_h->natoms; j++)
        {
            Vec3d p = atoms_h->apos[j];
            rot.dot_to(p, atoms_h->apos[j]);
            atoms_h->apos[j].add(cog);
        }
        
        delete[] ws;
        return dist;
    }

double compareAtoms(int h, int i){
    return compareAtoms(&atoms[h], &atoms[i]);
}

double compareAtoms(Atoms* atoms_h, int i)
{
    return compareAtoms(atoms_h, &atoms[i]);
}

    double compareDescriptors(int i, int j)
    {
        return descriptors[i].compareDescriptors(&descriptors[j]);
    };

    double compareDescriptors(Atoms* atoms, int i)
    {
        Descriptor d;
        for(int i = 0; i < nbOfusedDescriptors; i++)
            d.setDescriptor(usedDescriptors[i], atoms);
        return d.compareDescriptors(&descriptors[i]);
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
        printf("nMembers: %d\n", nMembers);
        printf("sameDescriptor: %d\n", sameDescriptor);
        //hashDescriptors(&d);

        return sameDescriptor;
    };

    int hashDescriptors(Descriptor* d){
        
        return 0;
    };
    // virtual char* tostr(int id){};
};

#endif