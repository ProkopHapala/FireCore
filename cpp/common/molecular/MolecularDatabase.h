#ifndef MolecularDatabase_h
#define MolecularDatabase_h

#include <cmath>
#include <chrono>  // for std::chrono


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
#include "LimitedGraph.h"

#include "MMFFBuilder.h"

enum DescriptorType
{
    principal_axes,
    number_of_atoms
};

struct MolecularDatabaseDescriptor
{
    DescriptorType type;
    int atom_type = -1;
};

class Descriptor : public MetaData
{
public:
    double *descs = nullptr;

    Descriptor(){};
    Descriptor(int n_)
    {
        alloc(n_);
    };
    ~Descriptor()
    {
        dealloc();
        delete[] descs;
    };
    void alloc(int n_)
    {
        dimensionOfDescriptorSpace = n_;
        descs = new double[dimensionOfDescriptorSpace];
    };
    void realloc(int n_)
    {
        if (descs == nullptr)
        {
            alloc(n_);
            return;
        }
        double *temp = new double[dimensionOfDescriptorSpace];
        delete[] descs;
        descs = temp;
    };
    void dealloc()
    {
        dimensionOfDescriptorSpace = 0;
    };

    void print()
    {
        printf("n: %d\n", dimensionOfDescriptorSpace);
        for (int i = 0; i < dimensionOfDescriptorSpace; i++)
        {
            printf("%lf ", descs[i]);
        }
        printf("\n");
    };

    void calculatePrincipalAxes(Atoms *a, int atom_type, MMFFparams *params = 0)
    {
        Vec3d *pos;
        int nbAtomsToPos = 0;
        if (atom_type > 0)
        {
            for (int i = 0; i < a->natoms; i++)
            {
                int ityp = a->atypes[i];
                if (params)
                    ityp = params->atypes[ityp].element;
                if (ityp == atom_type)
                {
                    nbAtomsToPos++;
                }
            }
            pos = new Vec3d[nbAtomsToPos];
            int j = 0;
            for (int i = 0; i < a->natoms; i++)
            {
                int ityp = a->atypes[i];
                if (params)
                    ityp = params->atypes[ityp].element;
                if (ityp == atom_type)
                {
                    pos[j] = a->apos[i];
                    j++;
                }
            }
        }
        else
        {
            pos = a->apos;
            nbAtomsToPos = a->natoms;
        }
        if (nbAtomsToPos < 3)
        {
            if (verbosity)
                printf("Not enough atoms to calculate principal axes\n");
            delete[] pos;
            return;
        }
        this->dimensionOfDescriptorSpace += 3;
        if (descs == nullptr)
        {
            alloc(dimensionOfDescriptorSpace);
        }
        else
        {
            double *temp = new double[dimensionOfDescriptorSpace - 3];
            std::copy(descs, descs + dimensionOfDescriptorSpace - 3, temp);
            realloc(dimensionOfDescriptorSpace);
            std::copy(temp, temp + dimensionOfDescriptorSpace - 3, descs);
            delete[] temp;
        }
        Mat3d I = Mat3dZero;
        double *ws = new double[nbAtomsToPos];
        for (int i = 0; i < nbAtomsToPos; i++)
        {
            // ws[i] = a->atypes[i];
            ws[i] = 1.0;
        }
        Vec3d cog = Vec3dZero;
        getCOG(nbAtomsToPos, pos, ws, cog);
        for (int i = 0; i < nbAtomsToPos; i++)
        {
            Vec3d r;
            Mat3d E = Mat3dIdentity;
            r.set_sub(pos[i], cog);
            I.addOuter(r, r, -ws[i]);
            I.add_mul(E, r.norm2() * ws[i]);
        }
        Vec3d e;
        I.eigenvals(e);
        e.sort();
        descs[dimensionOfDescriptorSpace - 3] = sqrt(e.x);
        descs[dimensionOfDescriptorSpace - 2] = sqrt(e.y);
        descs[dimensionOfDescriptorSpace - 1] = sqrt(e.z);
        delete[] ws, pos;
    }

    void calculateNatoms(Atoms *a, int atom_type, MMFFparams *params = 0)
    {
        this->dimensionOfDescriptorSpace += 1;
        if (descs == nullptr)
        {
            alloc(dimensionOfDescriptorSpace);
        }
        else
        {
            double *temp = new double[dimensionOfDescriptorSpace - 1];
            std::copy(descs, descs + dimensionOfDescriptorSpace - 1, temp);
            realloc(dimensionOfDescriptorSpace);
            std::copy(temp, temp + dimensionOfDescriptorSpace - 1, descs);
            delete[] temp;
        }
        if (atom_type > 0)
        {
            int nbAtomsOfTheType = 0;
            for (int i = 0; i < a->natoms; i++)
            {
                int ityp = a->atypes[i];
                if (params)
                    ityp = params->atypes[ityp].element;
                if (ityp == atom_type)
                {
                    nbAtomsOfTheType++;
                }
            }
            descs[dimensionOfDescriptorSpace - 1] = nbAtomsOfTheType;
        }
        else
        {
            descs[dimensionOfDescriptorSpace - 1] = a->natoms;
        }
    }

    void setDescriptor(MolecularDatabaseDescriptor md, Atoms *a, MMFFparams *params = 0)
    {
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

    double compareDescriptors(Descriptor *d)
    {
        if (this->dimensionOfDescriptorSpace != d->dimensionOfDescriptorSpace)
        {
            //     printf("Descriptors have different dimensions\n");
            //     printf("This: %d, Other: %d\n", this->dimensionOfDescriptorSpace, d->dimensionOfDescriptorSpace);
            return -1;
        }
        double sum = 0;
        for (int i = 0; i < dimensionOfDescriptorSpace; i++)
        {
            sum += (this->descs[i] - d->descs[i]) * (this->descs[i] - d->descs[i]);
        }
        return sqrt(sum);
    }

    void copyOf(const Descriptor &d)
    {
        realloc(d.dimensionOfDescriptorSpace);
        for (int i = 0; i < dimensionOfDescriptorSpace; i++)
        {
            this->descs[i] = d.descs[i];
        }
    }
};

class MolecularDatabase : public MetaData
{
private:
    Descriptor *descriptors = nullptr;
    Atoms *atoms = nullptr;
    int nMembers = 0;
    int nMaxCell = 0;
    MolecularDatabaseDescriptor *usedDescriptors = nullptr;
    int nbOfusedDescriptors = 0;
    

public:
    // ================= Constructor and Destructor ================= //
    std::vector<bool> convergedStructure;
    std::vector<int> structureOccurence;
    std::vector<std::chrono::high_resolution_clock::time_point> structureTimestamps;
    int totalEntries = 0;
    MolecularDatabase()
    {
        nMaxCell = 100;
        alloc(nMaxCell);
    };
    MolecularDatabase(int n)
    {
        alloc(n);
    };
    ~MolecularDatabase()
    {
        dealloc();
        delete[] usedDescriptors;
        delete[] descriptors;
        delete[] atoms;
    };
    void alloc(int n)
    {
        this->nMembers = 0;
        this->nMaxCell = n;
        this->atoms = new Atoms[n];
        this->descriptors = new Descriptor[n];
    };
    void realloc(int n)
    {

        Atoms *temp = new Atoms[n];
        for (int i = 0; i < nMembers; i++)
        {
            temp[i].copyOf(this->atoms[i]);
        }
        delete[] this->atoms;
        this->atoms = temp;

        Descriptor *temp2 = new Descriptor[n];
        for (int i = 0; i < nMembers; i++)
        {
            temp2[i].copyOf(this->descriptors[i]);
        }
        delete[] this->descriptors;
        this->descriptors = temp2;

        this->nMaxCell = n;
    };
    void dealloc()
    {
        for (int i = 0; i < nMembers; i++)
        {
            descriptors[i].dealloc();
            atoms[i].apos = 0;
            atoms[i].atypes = 0;
            atoms[i].dealloc();
        }
        atoms = nullptr;
        descriptors = nullptr;
        usedDescriptors = nullptr;
        nMembers = 0;
    };

    // ================= Accessors ================= //

    Atoms *GetAtoms(int i)
    {
        if (i < nMembers && i >= 0)
        {
            return &atoms[i];
        }
        else
        {
            printf("Cannot get atoms: Index %d out of range\n", i);
            return nullptr;
        }
    }


    

    void print()
    {
        printf("nMembers: %d\n", nMembers);
        printf("convergedStructure: %lu\n", convergedStructure.size());
        printf("%10s %10s %13s %13s %13s", "Member", "Occurence", "Age (sec)", "Occ/sec", "Converged?");
        for (int j = 0; j < dimensionOfDescriptorSpace; j++)
        {
            printf("     %-12s", ("Desc" + std::to_string(j)).c_str());
        }
        printf("\n");

        auto currentTime = std::chrono::high_resolution_clock::now();

        for (int i = 0; i < nMembers; i++)
        {
            double age = std::chrono::duration<double, std::ratio<1, 1>>(currentTime - structureTimestamps[i]).count();
            double occurencePerTime = structureOccurence[i] / age;
            printf("%10d %10d %13.1f %13.1f %13s", i, structureOccurence[i], age, occurencePerTime, convergedStructure[i] ? "true" : "false");
            for (int j = 0; j < dimensionOfDescriptorSpace; j++)
            {
                printf("     %-12lf", descriptors[i].descs[j]);
            }
            printf("\n");
        }
    };

    inline int getNMembers() { return nMembers; };

    // ================= Editing members of database ================= //

    void loadAtoms(int i, Atoms *a)
    {
        if (!atoms || i > nMembers || i < 0)
            return;
        if(atoms[i].natoms != a->natoms)
            a->realloc(atoms[i].natoms);
        a->copyOf(atoms[i]);
    };

    // User can choose desriptors to use by writing it in the usedDescriptor array (chosen order is important):
    // first parameter is the type of descriptor (principal_axes, number_of_atoms)
    // second parameter is the atom type (-1 means all atoms)
    void setDescriptors(Atoms *atoms = 0, MolecularDatabaseDescriptor *md = 0, int nbOfUsedDescriptors_ = 0)
    {
        if (verbosity > 1)
            printf("MolecularDatabase::setDescriptors\n");
        if (nbOfUsedDescriptors_ == 0)
        {
            nbOfusedDescriptors = 2;
            usedDescriptors = new MolecularDatabaseDescriptor[nbOfusedDescriptors];
            int ityp0 = atoms->atypes[0];
            int ityp1 = atoms->atypes[1];
            int ityp2 = atoms->atypes[2];
            usedDescriptors[0] = {principal_axes, -1};
            usedDescriptors[1] = {number_of_atoms, -1};
            // usedDescriptors[2] = {principal_axes, ityp1};
            // usedDescriptors[3] = {principal_axes, ityp2};
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

    void addMember(Atoms *a, Descriptor *descriptor = 0, MMFFparams *params = 0)
    {
        if (nMembers >= nMaxCell)
        {
            if (verbosity > 1)
                printf("Reallocating to maximaly store %d molecules\n", nMaxCell * 2);
            realloc(nMaxCell * 2);
        }
        if (!descriptor)
            for (int i = 0; i < nbOfusedDescriptors; i++)
                descriptors[nMembers].setDescriptor(usedDescriptors[i], a, params);
        else
            descriptors[nMembers].copyOf(*descriptor);

        atoms[nMembers].copyOf(*a);
        if (nMembers == 0)
            dimensionOfDescriptorSpace = descriptors[nMembers].dimensionOfDescriptorSpace;
        this->nMembers++;
        structureTimestamps.push_back(std::chrono::high_resolution_clock::now());
    };


    int addIfNewDescriptor(Atoms *a, Mat3d* lat_vec = 0, MMFFparams *params = 0)
    {
        totalEntries++;
        Descriptor d;
        int sameDescriptor = -1;
        // if(verbosity) printf("testHash\n");
        for (int i = 0; i < nbOfusedDescriptors; i++)
        {
            d.setDescriptor(usedDescriptors[i], a, params);
        }
        for (int i = 0; i < nMembers; i++)
        {
            if (d.compareDescriptors(&descriptors[i]) < 5)
            {
                if(!lat_vec){
                    if (computeDistance(a, &atoms[i])/a->natoms < 0.1)
                    {
                        //printf("computeDistance: %lf\n", computeDistance(a, &atoms[i]));
                        sameDescriptor = i;
                        break;
                    }
                }
                else{
                    if(computeDistanceOnSurf(a, &atoms[i], lat_vec)/a->natoms < 0.1)
                    {
                        sameDescriptor = i;
                        break;
                    }
                }
            }
        }
        if (sameDescriptor == -1)
        {
            addMember(a, &d);
            structureOccurence.push_back(1);
        }
        else
        {
            structureOccurence[sameDescriptor]++;
        }
        // hashDescriptors(&d);

        return sameDescriptor;
    };

    int hashDescriptors(Descriptor *d)
    {
        return 0;
    };

    void replace(Atoms *a, int i, MMFFparams *params = 0)
    {
        if (i < nMembers && i >= 0)
        {
            atoms[i].realloc(a->natoms);
            atoms[i].copyOf(*a);

            descriptors[i].dealloc();
            for (int j = 0; j < nbOfusedDescriptors; j++)
                descriptors[i].setDescriptor(usedDescriptors[j], &atoms[i], params);
        }
        else
        {
            printf("Cannot replace: Index %d out of range\n", i);
        }
    }


    // ================= Comparing structures ================= //

    double compareDescriptors(int i, int j) { return descriptors[i].compareDescriptors(&descriptors[j]); }

    double compareDescriptors(Atoms *atoms, int j)
    {
        Descriptor d;
        for (int i = 0; i < nbOfusedDescriptors; i++)
            d.setDescriptor(usedDescriptors[i], atoms);
        return d.compareDescriptors(&descriptors[j]);
    };

    double computeDistance(int h, int i) { return computeDistance(&atoms[h], &atoms[i]); }

    double computeDistance(Atoms *atoms_h, int i) { return computeDistance(atoms_h, &atoms[i]); }

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

    double computeDistance(Atoms *atoms_h, Atoms *atoms_i)
    {
        if (atoms_h == nullptr || atoms_i == nullptr)
        {
            printf("Null pointer passed\n");
            return -1;
        }
        int natm = atoms_h->natoms;
        if (abs(natm-atoms_i->natoms)>0.5)
        {
            printf("Different number of atoms\n");
            return -1;
        }
        double dist = 0;

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

        Mat3d rot;
        rot = superimposer(atoms_h->apos, atoms_i->apos, atoms_h->natoms);
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
        for (int j = 0; j < natm; j++)
        {
            atoms_h->apos[j] = atoms_h->apos[j] + x0_center;
            atoms_i->apos[j] = atoms_i->apos[j] + x1_center;
        }
        delete[] ws;
        return dist;
    }

    Vec3d shift_to_elementary_cell(Atoms *atoms, Mat3d *lat_vec, bool Surf=0)
    {
        Mat3d A;
        Vec3d shift = Vec3dZero;
        Vec3d sh;
        lat_vec->invert_to(A);
        shift = A.dot(atoms->apos[0]);
        for (int k = 0; k < 2; k++)
        {
            if (floor(abs(shift.array[k]) < 1))
            {
                continue;
            }
            for (int j = 0; j < atoms->natoms; j++)
            {
                sh.set_mul(lat_vec->vecs[k], floor(shift.array[k]));
                // printf("floor(shift.array[k] + 1e-4): %g\n", floor(shift.array[k] + 1e-4));
                atoms->apos[j].sub(sh);
            }
        }
        return shift;
    }

    void symmetry4mm(Atoms *atoms, int symmetry_ind){
        if(symmetry_ind < 0 || symmetry_ind > 7){
            printf("Symmetry index out of range\n");
            return;
        }
        Mat3d A[8] = {
            Mat3dIdentity,

            // mirror (-x,-y)
            {-1,0,0,
            0,-1,0,
            0,0,1},

            // rotation 90 degrees
            {0,-1,0,
            1,0,0,
            0,0,1},

            {0,1,0,
            -1,0,0,
            0,0,1},

            {-1,0,0,
            0,1,0,
            0,0,1},

            {1,0,0,
            0,-1,0,
            0,0,1},

            {0,1,0,
            1,0,0,
            0,0,1},

            {0,-1,0,
            -1,0,0,
            0,0,1}
        };
        for (int j = 0; j < atoms->natoms; j++)
        {
            A[symmetry_ind].dot_to(atoms->apos[j], atoms->apos[j]);
        }
    }

    double computeDistanceOnSurf(int h, int i, Mat3d* lat_vec){
        return computeDistanceOnSurf(&atoms[h], &atoms[i], lat_vec);
    }

    double computeDistanceOnSurf(Atoms *atoms_h, Atoms *atoms_i, Mat3d* lat_vec){
        // move both apos[0] into elementary cell
        Mat3d A;
        int n = atoms_h->natoms;
        Vec3d shift_h = shift_to_elementary_cell(atoms_h, lat_vec);
        Vec3d shift_i = shift_to_elementary_cell(atoms_i, lat_vec);

        // try all symmetries to find the smallest RMSD
        AtomicConfiguration ac, ac2;
        double dist = __DBL_MAX__;
        Atoms a_h, a_i;
        for(int sym = 0; sym < 8; sym++){
            a_h.copyOf(*atoms_h);
            a_i.copyOf(*atoms_i);
            symmetry4mm(&a_h, sym); //TODO: all point groups should be implemented
            ac.natoms = a_h.natoms;
            ac.types = a_h.atypes;
            ac.pos = a_h.apos;

            ac2.natoms = a_i.natoms;
            ac2.types = a_i.atypes;
            ac2.pos = a_i.apos;

            if(ac.dist(ac2) < dist){
                dist = ac.dist(ac2);
            }
        }
        dist = sqrt(dist);

        //revert back to original positions
        for (int k = 0; k < 2; k++)
        {
            if (floor(abs(shift_h.array[k]) < 1))
            {
                continue;
            }
            for (int j = 0; j < n; j++)
            {
                Vec3d sh;
                sh.set_mul(lat_vec->vecs[k], floor(shift_h.array[k]));
                atoms_h->apos[j].add(sh);
            }
        }
        for (int k = 0; k < 2; k++)
        {
            if (floor(abs(shift_i.array[k]) < 1))
            {
                continue;
            }
            for (int j = 0; j < n; j++)
            {
                Vec3d sh;
                sh.set_mul(lat_vec->vecs[k], floor(shift_i.array[k]));
                atoms_i->apos[j].add(sh);
            }
        }
        return dist;
    }

    /*
    Transformed code from https://github.com/willhooper/superpose/tree/master
    Works similar as Kabsch algorithm
    Aligns two sets of coordinates by symetrizing matrix from Kabsch algorithm
    Symmetric matrix is already diagonal, but in diferent orthogonal basis
    */
    Mat3d superimposer(Vec3d *coord0, Vec3d *coord1, unsigned int natm)
    {

        float tolerance = 0.0001;

        // Error checking
        unsigned int MAXATOM = 200;
        if (natm > MAXATOM)
        {
            printf("Error! Increase maxatom\n");
            throw 0;
        }
        if (natm < 4)
        {
            printf("Error! Too few atoms (must be 4 or greater)\n");
            throw 0;
        }

        double **x0_ = new double *[3];
        double **x1_ = new double *[3];
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
        aa.symetrize_matrix_SVD(rotation);
        for (int i = 0; i < 3; i++)
        {
            delete[] x0_[i];
            delete[] x1_[i];
        }
        delete[] x0_;
        delete[] x1_;
        return rotation;
    }

    // ================= Development ================= //

    void as_rigid_as_possible(Atoms *atoms_h_, int index, int nBonds, Vec2i *bonds, int nbFixed=0, int* fixed =0, int *neighs = 0)
    {
        const int nNeighs = 4;
        //printf("as_rigid_as_possible::Aligning new structure %d to atoms[%d]\n", nMembers - 1, index);
        // for(int i = 0; i < nBonds; i++)
        //     printf("Bond %d: %d %d\n", i, bonds[i].a, bonds[i].b);
        if (!neighs)
        {
            neighs = new int[nNeighs * atoms_h_->natoms];
            for (int i = 0; i < atoms_h_->natoms; i++)
            {
                for (int j = 0; j < nNeighs; j++)
                {
                    neighs[nNeighs * i + j] = -1;
                }
            }
            for (int i = 0; i < nBonds; i++)
            {
                int currentNeighs_a, currentNeighs_b;
                for (int j = 0; j < nNeighs; j++)
                {
                    if (neighs[nNeighs * bonds[i].a + j] == -1)
                    {
                        currentNeighs_a = j;
                        break;
                    }
                }
                for (int j = 0; j < nNeighs; j++)
                {
                    if (neighs[nNeighs * bonds[i].b + j] == -1)
                    {
                        currentNeighs_b = j;
                        break;
                    }
                }
                neighs[nNeighs * bonds[i].a + currentNeighs_a] = bonds[i].b;
                neighs[nNeighs * bonds[i].b + currentNeighs_b] = bonds[i].a;
            }
        }
        Atoms *atoms_h = atoms_h_;
        // atoms_h.copyOf(atoms_h_);
        Atoms *atoms_i = &atoms[index];
        int natm = atoms_h->natoms;
        double E = 0;
        Mat3d S = Mat3dZero;
        Mat3d *R = new Mat3d[natm];
        Mat3d U, V;
        Vec3d *F = new Vec3d[natm];
        Vec3d temp;
        Mat3d Rot;
        Quat4d *w_mat = new Quat4d[natm];
        double *w_vec = new double[natm];
        double alpha = 0.001 ;
        Vec3d e_h, e_i;
        for (int j = 0; j < natm; j++)
        {
            R[j] = Mat3dZero;
            double nbNeighs = 0;
            for (int k = 0; k < nNeighs; k++)
            {
                if (neighs[nNeighs * j + k] == -1)
                    break;
                nbNeighs += 1;
            }
            w_mat[j] = Quat4dOnes;
            //w_mat[j].mul(4/nbNeighs);
            w_mat[j].mul(exp(1/nbNeighs));
            w_vec[j] = 1;
        } // ToDo: set weights

        //atoms_h->print();
        //atoms_i->print();

        const int qmax = 1;
        int q = 0;
        while (qmax > q)
        {
            E = 0;
            for (int j = 0; j < natm; j++)
            {
                S = Mat3dZero;
                bool capping = false;
                for (int k = 0; k < nNeighs; k++)
                {
                    if (k == 0 && neighs[nNeighs * j + 1] == -1)
                    {
                        capping = true;
                        break;
                    }
                    if (neighs[nNeighs * j + k] == -1)
                        break;

                    e_h = atoms_h->apos[j] - atoms_h->apos[neighs[nNeighs * j + k]];
                    e_i = atoms_i->apos[j] - atoms_i->apos[neighs[nNeighs * j + k]];
                    //printf("e_h: %lf %lf %lf\n", e_h.x, e_h.y, e_h.z);
                    //printf("e_i: %lf %lf %lf\n", e_i.x, e_i.y, e_i.z);
                    S.addOuter(e_h, e_i, w_mat[j].array[k]);
                    // printf("S: %d(%d)\n", j, k);
                    // S.print();
                }
                //if (capping){continue;}
                //  printf("U: %d\n", j);
                //  U.print();
                //  printf("V: %d\n", j);
                //  V.print();
                S.SVD(U, temp, V);
                // if(capping){
                // Mat3d StS;
                // StS.set_mmul_TN(S, S);
                // printf("SSt: %d\n", j);
                // StS.print();
                //  printf("U: %d\n", j);
                //  U.print();
                //  printf("V: %d\n", j);
                //  V.print();
                // }
                //  printf("temp: %lf %lf %lf\n", temp.x, temp.y, temp.z);
                //  printf("R: %d\n", j);
                
    
                R[j].set_mmul_NT(U,V);
                //printf("R[%d].determinant(): %lf\n", j, R[j].determinant());
                //R[j].print();


                double E_j = 0;
                for (int k = 0; k < nNeighs; k++)
                {   
                    if(capping)
                        break;
                    if (neighs[nNeighs * j + k] == -1)
                        break;
                    e_h = atoms_h->apos[j] - atoms_h->apos[neighs[nNeighs * j + k]];
                    e_i = atoms_i->apos[j] - atoms_i->apos[neighs[nNeighs * j + k]];
                    double d = (e_h - R[j].dot(e_i)).norm2();
                    E_j += w_mat[j].array[k] * d;
                }
                E += w_vec[j] * E_j;
            }
            for (int j = 0; j < natm; j++)
            {
                //printf("R[%d]:\n", j);
                //R[j].print();
                F[j] = Vec3dZero;
                for (int k = 0; k < nNeighs; k++)
                {
                    Rot = Mat3dZero;
                    if (neighs[nNeighs * j + 1] == -1)
                        Rot.add(R[neighs[nNeighs * j + k]]);
                    if (neighs[nNeighs * j + k] == -1)
                        break;
                    e_h = atoms_h->apos[j] - atoms_h->apos[neighs[nNeighs * j + k]];
                    e_i = atoms_i->apos[j] - atoms_i->apos[neighs[nNeighs * j + k]];
                    // printf("e_h: %lf %lf %lf\n", e_h.x, e_h.y, e_h.z);
                    // printf("e_i: %lf %lf %lf\n", e_i.x, e_i.y, e_i.z);
                    
                    Rot.add(R[j]);
                    // printf("Rot: %d\n", j);
                    // Rot.print();
                    // printf("\n");

                    if(neighs[nNeighs*neighs[nNeighs * j + k]+1] != -1)
                        Rot.add(R[neighs[nNeighs * j + k]]);
                    else
                        Rot.add(R[j]);
                    
                    Rot.mul(0.5);
                    // printf("Rot:\n");
                    // Rot.print();

                    //printf("Rot.dot(e_h).mul(-0.5): %lf %lf %lf\n", Rot.dot(e_h).mul(-0.5).x, Rot.dot(e_h).mul(-0.5).y, Rot.dot(e_h).mul(-0.5).z);
                    //printf("Rot.dot(e_h).mul(-0.5) + e_i: %lf %lf %lf\n", (Rot.dot(e_h).mul(-0.5) + e_i).x, (Rot.dot(e_h).mul(-0.5) + e_i).y, (Rot.dot(e_h).mul(-0.5) + e_i).z);
                    //printf("Rot.dot(e_h).mul(0.5): %lf %lf %lf\n", Rot.dot(e_h).mul(0.5).x, Rot.dot(e_h).mul(0.5).y, Rot.dot(e_h).mul(0.5).z);
                    F[j].add(e_h - Rot.dot(e_i)).mul(4 * w_mat[j].array[k]);
                }
            }
            //atoms_h->print();
            for (int j = 0; j < natm; j++)
            {//printf("%d F: %lf %lf %lf\n", j, F[j].x, F[j].y, F[j].z);
                for(int f = 0; f < nbFixed; f++){
                    if(j == fixed[f]){
                        F[j] = Vec3dZero;
                        atoms_h->apos[j] = atoms_i->apos[j];
                    }
                }
                // if(neighs[nNeighs * j + 1] == -1)
                // {
                //     continue;
                // }
                // else
                // {
                    atoms_h->apos[j] = atoms_h->apos[j] - F[j].mul(alpha);
                //}
            }
            
            // for (int j = 0; j < natm; j++)
            // {
            //     if(neighs[nNeighs * j + 1] == -1)
            //     {
            //         atoms_h->apos[j] = atoms_h->apos[neighs[nNeighs * j + 0]]  - R[neighs[nNeighs * j + 0]].dot(atoms_i->apos[neighs[nNeighs * j + 0]] - atoms_i->apos[j]);
            //     }
            // }
            q++;
            //atoms_h->print();
        }
        //atoms_h->print();
        //atoms_i->print();
        //for(int i = 0; i < natm; i++)
        // {
        //     printf("F[%d]: %lf %lf %lf\n", i, F[i].x, F[i].y, F[i].z);
        // }
        Vec3d F_tot = Vec3dZero;
        for(int i = 0; i < natm; i++)
        {
            F_tot.add(F[i]);
        }
        printf("%d: E = %lf, F = %lf\n", q, E, F_tot.norm());
        delete[] R;
        //printf("delete R\n");
        delete[] w_mat;
        //printf("delete w_mat\n");
        delete[] w_vec;
        //printf("delete w_vec\n");
        delete[] F;
        //printf("as_rigid_as_possible::End\n");
        // printf("S: %lf %lf %lf\n", S.x, S.y, S.z);
    }
};

#endif

// void FindRotation(Mat3d &rot, Atoms *a)
// {
//     Mat3d XX = Mat3dZero;
//     for (int i = 0; i < a->natoms; i++)
//     {
//         XX.addOuter(a->apos[i], a->apos[i], 1.0);
//     }
//     Vec3d evs;
//     XX.eigenvals(evs);

//     evs.sort();

//     XX.eigenvec(evs.x, rot.a);
//     XX.eigenvec(evs.y, rot.b);
//     XX.eigenvec(evs.z, rot.c);
// }
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