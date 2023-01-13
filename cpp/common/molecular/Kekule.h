
#ifndef Kekule_h
#define Kekule_h

#include <stdint.h>

/*
#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
*/

#include "Vec2.h"

class Kekule{ public:
    //int ialg=1;
    double Katom=1.0,Kbond=1.0,KatomInt=0.1,KbondInt=0.1;
    double damping=0.1;
    double Etot=0,Eb=0,Ea=0;
    int natom=0;
    int nbond=0;

    // Params (targets)
    Vec2i*  bond2atom      =0; // [nbond] atoms coresponding to the bond;
    double* bondOrderMin   =0; // [nbond] wanted obnd order of that bond
    double* bondOrderMax   =0;
    //double* bondOrderK   =0; // [nbond] stiffness of that bond order
    double* atomValenceMin =0; // [natom] how many bonds atom can have?
    double* atomValenceMax =0; // [natom] how many bonds atom can have?
    //double* atomValenceK =0; // [natom] stiffness of that bond order

    // dynamical variables
    double* bondOrder    =0; // [natom]
    double* bondOrderF   =0; // [nbond]
    double* bondOrderV   =0; // [nbond]
    double* atomValence  =0; // [natom]
    double* atomValenceF =0; // [natom]

    void realloc( int natom_, int nbond_, Vec2i* bond2atom_=0, double* bondOrderMin_=0, double* bondOrderMax_=0,  double* atomValenceMin_=0, double* atomValenceMax_=0 ){
        natom=natom_; nbond=nbond_;

        _bindOrRealloc( nbond, bond2atom_, bond2atom );

        _bindOrRealloc( natom, atomValenceMin_, atomValenceMin );
        _bindOrRealloc( natom, atomValenceMax_, atomValenceMax );
        //_bindOrRealloc( natom, atomValence0_, atomValence0   );
        //_bindOrRealloc( natom, atomValenceK_, atomValenceK   );
        
        _bindOrRealloc( nbond, bondOrderMin_,   bondOrderMin );
        _bindOrRealloc( nbond, bondOrderMax_,   bondOrderMax );
        //_bindOrRealloc( nbond, bondOrder0_,   bondOrder0   );
        //_bindOrRealloc( nbond, bondOrderK_,   bondOrderK   );

        _realloc(atomValence ,natom);
        _realloc(atomValenceF,natom);
        //_realloc(atomValence0,natom);
        //_realloc(atomValenceK,natom);
        //_realloc(bond2atom ,nbond);

        _realloc(bondOrder ,nbond);
        _realloc(bondOrderF,nbond);
        _realloc(bondOrderV,nbond);
        //_realloc(bondOrderK,nbond);
        //_realloc(bondOrder0,nbond);

    }

    void clearAtomValence(){ for(int i=0; i<natom; i++){ atomValence[i]=0; }; }
    void clearVelocity   (){ for(int i=0; i<nbond; i++){ bondOrderV [i]=0; }; }
    void setRandomAtoms  (){ for(int i=0; i<natom; i++){ atomValence[i]=randf(atomValenceMin[i],atomValenceMax[i]); } }
    void setRandomBonds  (){ for(int i=0; i<nbond; i++){ bondOrder  [i]=randf(bondOrderMin[i],bondOrderMax[i]);     } }
    void setDefaultBondOrders(double min, double max){  for(int i=0; i<nbond; i++){ bondOrderMax[i]=max; bondOrderMin[i]=min; }}

    void pinBondOrders(int n, int* ibonds, int* target ){
        //if(bSetZeros)for(int i=0; i<nbond; i++){ bondOrder0[i]=1; bondOrderK[i]=0; };
        for(int i=0; i<n; i++){
            int ib = ibonds[i];
            bondOrderMin[ib] = target[i];
            bondOrderMax[ib] = target[i];
            //bondOrderK[ib] = K;
        }
    }

    void evalAtomForces(){
        Ea=0;
        for(int i=0; i<natom; i++){
            double val  = atomValence[i];
            double dint = val - (int)(val+0.5);
            double dmin = val - atomValenceMin[i]; if(dmin>0)dmin=0;
            double dmax = val - atomValenceMax[i]; if(dmax<0)dmax=0;
            atomValenceF[i] = -((dmin+dmax)*Katom + dint*KatomInt);
            double e=0.5*( Katom*( dmin*dmin + dmax*dmax )  +  KatomInt*dint*dint   );
            Ea+=e;
        }
    }
    
    void evalBondForces(){
        Eb=0;
        for(int i=0; i<nbond; i++){
            double val = bondOrder[i];
            double dint = val - (int)(val+0.5);
            double dmin = val - bondOrderMin[i];   if(dmin>0)dmin=0;
            double dmax = val - bondOrderMax[i];   if(dmax<0)dmax=0;
            bondOrderF[i] = -(( dmin + dmax )*Kbond + dint*KbondInt);
            double e=0.5*( Kbond*( (dmin*dmin) + (dmax*dmax) )  +  KbondInt*dint*dint   );
            if(verbosity>2){  printf( "val %g dmin %g dmax %g dint %g E %g\n", val, dmin, dmax, dint, e ); }
            Eb+=e;
        }
    }

    void projectValence(){
        //for(int i=0; i<natom; i++){ atomValence[i]=0; };
        clearAtomValence();
        for(int i=0; i<nbond; i++){
            const Vec2i& b = bond2atom[i];
            double bo      = bondOrder[i]; 
            atomValence[b.a] += bo;
            atomValence[b.b] += bo;
        }
    }

    void projectAtomForces(){
        for(int i=0; i<nbond; i++){
            const Vec2i& b = bond2atom[i];
            bondOrderF[i] += atomValenceF[b.a] + atomValenceF[b.b];
        }
    }

    double move_GD( double dt ){
        double F2sum = 0;
        for(int i=0; i<nbond; i++){
            double f = bondOrderF[i];
            bondOrder[i] += f * dt;
            //printf( "bf[%i] %g \n", i, f );
            F2sum+=f*f;
        }
        return F2sum;
    }

    double move_MDdamp( double dt, double damping ){
        double F2sum = 0;
        double damp = 1-damping;
        for(int i=0; i<nbond; i++){
            double f = bondOrderF[i];
            double v = bondOrderV[i];
            v *= damp;
            v += f * dt;
            bondOrderV[i]  = v;
            bondOrder [i] += v * dt;
            //printf( "bf[%i] %g \n", i, f );
            F2sum+=f*f;
        }
        return F2sum;
    }

    double eval(){
        evalBondForces();
        projectValence();
        evalAtomForces();
        projectAtomForces();
        Etot=Ea+Eb;
        //printf("Etot, Ea,Eb %g %g %g \n", Etot, Ea,Eb);
        return Etot;
    }

    double update( double dt, int ialg ){
        eval();
        double F2;
        switch(ialg){
            case 0: F2=move_GD    ( dt          ); break;
            case 1: F2=move_MDdamp( dt, damping ); break;
        }
        return F2;
    }

    double relax( double dt, double F2conv=1e-6, int maxIter=1000, bool bRandStart=true, int ialg=1 ){
        if(bondOrderV)clearVelocity();
        if(bRandStart){
            setRandomAtoms();
            setRandomBonds();
        }
        double F2sum=0;
        for(int i=0; i<maxIter; i++){
            F2sum = update( dt, ialg );
            if(verbosity>0){ printf("Iter %i Ea %g Eb %g F2sum %g ", i, Ea, Eb, F2sum); if(verbosity>1){printBondOrders();}else{ printf("\n"); } }
            if(F2sum<F2conv) break;
        }
        //printBondOrders();
        return F2sum;
    }


    void printBondOrders(){
        printf("BOs[");
        for(int i=0; i<nbond; i++){
            //printf( "bondOrder[%i]= %g \n", i, bondOrder[i] );
            printf( "%1.2f ",  bondOrder[i] );
        }
        printf("]\n");
    }

}; // Kekule

#endif



