
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
    double Katom=1.0;
    double Kbond=1.0;
    double KatomInt=0.1;
    double KbondInt=0.1;
    //double damping=0.1;
    double Etot=0.0;
    double Eb=0.0;
    double Ea=0.0;
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
    double* bondOrder    =0; // [nbond]
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

    void clear( bool bB2A=false, bool bAO=true, bool bBO=true ){
        if(bB2A) _dealloc( bond2atom  );
        if(bAO){
            _dealloc( atomValenceMin );
            _dealloc( atomValenceMax);
            //_dealloc( atomValence0);
            //_dealloc( atomValenceK);
        }
        if(bBO){
            _dealloc( bondOrderMin );
            _dealloc( bondOrderMax );
        }
        _dealloc( bondOrder  );
        _dealloc( bondOrderF );
        _dealloc( bondOrderV );
        //_dealloc( bondOrder0 );
        //_dealloc( bondOrderK );
        _dealloc( atomValence  );
        _dealloc( atomValenceF );
        //_dealloc( atomValence0 );
        //_dealloc( atomValenceK );
    }
    
    void clearAtomValence(){ for(int i=0; i<natom; i++){ atomValence[i]=0; }; }
    void clearVelocity   (){ for(int i=0; i<nbond; i++){ bondOrderV [i]=0; }; }
    void setRandomAtoms  (){ for(int i=0; i<natom; i++){ atomValence[i]=randf(atomValenceMin[i],atomValenceMax[i]); } }
    void setRandomBonds  (){ for(int i=0; i<nbond; i++){ bondOrder  [i]=randf(bondOrderMin[i],bondOrderMax[i]);     } }
    
    void pinBondOrders(int n, int* ibonds, int* target ){
        //if(bSetZeros)for(int i=0; i<nbond; i++){ bondOrder0[i]=1; bondOrderK[i]=0; };
        for(int i=0; i<n; i++){
            int ib = ibonds[i];
            bondOrderMin[ib] = target[i];
            bondOrderMax[ib] = target[i];
            //bondOrderK[ib] = K;
        }
    }

    // --- init
    void setDefaultBondOrders(double min, double max, double BO0=-1. ){  
        for(int ib=0; ib<nbond; ib++){ 
            bondOrderMax[ib]=max; 
            bondOrderMin[ib]=min;  
            if(BO0>0){ bondOrder[ib]=BO0; };
        }
    }
    void setValences(){
        for(int ia=0; ia<natom; ia++){ 
            atomValence[ia]=0.0; 
        }
        for(int ib=0; ib<nbond; ib++){
            const Vec2i& b = bond2atom[ib];
            atomValence[b.a] += bondOrder[ib];
            atomValence[b.b] += bondOrder[ib];
        }
    }
    // --- end init

    // --- eval
   void evalBondForces( bool hack15 ){
        Eb=0;
        for(int ib=0; ib<nbond; ib++){
            double val = bondOrder[ib];
            double dint;
            if(hack15 && val>1.25 && val<1.75) {
                dint = val - 1.5;
            }else{
                dint = val - (int)(val+0.5);
            }
            double dmin = val - bondOrderMin[ib];   
            if(dmin>0.0)dmin=0.0;
            double dmax = val - bondOrderMax[ib];   
            if(dmax<0.0)dmax=0.0;
            bondOrderF[ib] = -(( dmin + dmax )*Kbond + dint*KbondInt);
            double e=0.5*( Kbond*( (dmin*dmin) + (dmax*dmax) ) + KbondInt*dint*dint );
            Eb+=e;
        }
    }

    void projectValence(){
        clearAtomValence();
        for(int ib=0; ib<nbond; ib++){
            const Vec2i& b = bond2atom[ib];
            double bo      = bondOrder[ib]; 
            atomValence[b.a] += bo;
            atomValence[b.b] += bo;
        }
    }

    void evalAtomForces(){
        Ea=0.0;
        for(int ia=0; ia<natom; ia++){
            double val  = atomValence[ia];
            double dint = val - (int)(val+0.5);
            double dmin = val - atomValenceMin[ia]; 
            if(dmin>0)dmin=0;
            double dmax = val - atomValenceMax[ia]; 
            if(dmax<0)dmax=0;
            atomValenceF[ia] = -((dmin+dmax)*Katom + dint*KatomInt);
            double e=0.5*( Katom*( dmin*dmin + dmax*dmax ) + KatomInt*dint*dint );
            Ea+=e;
        }
    }

    void projectAtomForces(){
        for(int ib=0; ib<nbond; ib++){
            const Vec2i& b = bond2atom[ib];
            bondOrderF[ib] += atomValenceF[b.a] + atomValenceF[b.b];
        }
    }

    double eval( bool hack15=false ){
        evalBondForces( hack15 );
        projectValence();
        evalAtomForces();
        projectAtomForces();
        return Ea+Eb;
    }
    // --- end eval

    // --- move
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

    double move_MDdamp( double dt, double damp ){
        double F2sum = 0;
        for(int ib=0; ib<nbond; ib++){
            double f = bondOrderF[ib];
            double v = bondOrderV[ib];
            v *= damp;
            v += f * dt;
            bondOrderV[ib]  = v;
            bondOrder [ib] += v * dt;
            F2sum+=f*f;
        }
        return F2sum;
    }
    // --- end move

    double update( double dt, double damp, int ialg, bool hack15 ){
        eval( hack15 );
        double F2;
        switch(ialg){
            case 0: F2=move_GD    ( dt          ); break;
            case 1: F2=move_MDdamp( dt, damp ); break;
        }
        return F2;
    }

    double relax( double dt, double damping=0.1, double F2conv=1e-6, int maxIter=1000, bool bRandStart=true, int ialg=1, bool hack15=false ){
        if(bondOrderV)clearVelocity();
        if(bRandStart){
            setRandomAtoms();
            setRandomBonds();
        }
        double F2sum=0.0;
        double damp = 1.0-damping;
        for(int i=0; i<maxIter; i++){
            F2sum = update( dt, damp, ialg, hack15 );
            if(verbosity>2){ printf("Iter %i Ea %g Eb %g F2sum %g ", i, Ea, Eb, F2sum); if(verbosity>1){printBondOrders();}else{ printf("\n"); } }
            if(F2sum<F2conv) break;
        }
        return F2sum;
    }

    void roundBondOrders( double tol=0.01 ){
        bool ok=true;
        for(int ib=0; ib<nbond; ib++){
            if( fabs(bondOrder[ib]-1.0) < tol ){
                bondOrder[ib] = 1.0;
            }else if( fabs(bondOrder[ib]-1.5) < tol ){
                bondOrder[ib] = 1.5;
            }else if( fabs(bondOrder[ib]-2.0) < tol ){
                bondOrder[ib] = 2.0;
            }else if( fabs(bondOrder[ib]-3.0) < tol ){
                bondOrder[ib] = 3.0;
            }else{
                ok=false;
                printf( "ERROR roundBondOrders: bondOrder[%i]= %g \n", ib, bondOrder[ib]);
            }
        }
        if(!ok){exit(0);}
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



