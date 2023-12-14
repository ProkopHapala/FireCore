#ifndef MolecularDatabase_h
#define MolecularDatabase_h

#include <string.h>
#include <stdio.h>
#include "Vec3.h"
#include "Mat3.h"
#include "Atoms.h"
#include "Vec3Utils.h"

class MolecularDatabase : public Atoms
{
private:
    int* key;
    DatabaseMember* members;
    int nMembers;

public:
    MolecularDatabase();
    ~MolecularDatabase();
    void search();
    void duplicateCheck();
    void addMember(double E, Vec3d& F, Atoms& atoms);
    void addMember(DatabaseMember& member);
    int getNMembers(){ return this->nMembers; };
    DatabaseMember* getMember(int i){ return &this->members[i]; };

};




class DatabaseMember : public MolecularDatabase{
    public:
    int    key;
    double E;
    Vec3d  F;
    Atoms* atoms;
    Vec3d principalAxes;
    DatabaseMember(double E, Vec3d& F, Atoms& atoms){
        this->E = E;
        this->F = F;
        this->atoms = &atoms;
        calculatePrincipalAxes();
    };
    //DatabaseMember(int natoms){ realloc(natoms); };
    //DatabaseMember(const DatabaseMember& p){ copyOf(p);  };
    //~DatabaseMember(){};

    /**
     * Calculates the principal axes of the molecular database.
     * This function calculates the principal axes of the molecular database
     * using the moment of inertia tensor and stores the result in the
     * `principalAxes` variable.
     * @f[I=(r_i)^2*E-r_i⊗r_i@f] 
     *  
     * 
     */
    void calculatePrincipalAxes(){
        int natoms = atoms->natoms;
        Mat3d I;
        I.set(0.0);
        Vec3d ws = Vec3dOne;
        Vec3d cog = Vec3dZero;
        getCOG(natoms, atoms->apos, ws, cog);
        for(int i=0; i<natoms; i++){
            I.xx += (apos[i].y-cog.y)*(apos[i].y-cog.y) + (apos[i].z-cog.z)*(apos[i].z-cog.z);
            I.xy += -(apos[i].x-cog.x)*(apos[i].y-cog.y);
            I.xz += -(apos[i].x-cog.x)*(apos[i].z-cog.z);
            I.yy += (apos[i].x-cog.x)*(apos[i].x-cog.x) + (apos[i].z-cog.z)*(apos[i].z-cog.z);
            I.yz += -(apos[i].y-cog.y)*(apos[i].z-cog.z);
            I.zz += (apos[i].x-cog.x)*(apos[i].x-cog.x) + (apos[i].y-cog.y)*(apos[i].y-cog.y);
        }
        I.eigenValues(principalAxes);
    }

    Vec3d getPrincipalAxes(){ return this->principalAxes;};


};




#endif