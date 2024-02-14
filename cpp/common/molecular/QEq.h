
#ifndef QEq_h
#define QEq_h

#include "fastmath.h"
#include "Vec3.h"

// implementation of "ChargeEquilibrationforMolecularDynamicsSimulations" 
// http://www.sklogwiki.org/SklogWiki/index.php/Charge_equilibration_for_molecular_dynamics_simulations
// ChargeEquilibrationforMolecularDynamicsSimulations
// https://pubs.acs.org/doi/pdf/10.1021/j100161a070
// AnthonyK.Rappe, William A.Goddard, Phys.Chem.1991,95,3358-3363

void makeCoulombMatrix(int n, Vec3d* ps, double* J){
    //double Jmax = 0.0;
    for(int i=0; i<n; i++){
        Vec3d pi = ps[i];
        J[i*n+i]=0;
        for(int j=i+1; j<n; j++){
            Vec3d  d = ps[j] - pi;
            double Jij = 14.3996448915/(1+d.norm());        // ToDo: there should be other (gaussian, erf) damping
            //double Jij = 14.3996448915/sqrt(1+d.norm2());
            //printf("Jij[%i,%i] %g \n", i,j, Jij);
            //Jmax = fmax(Jmax,Jij);
            J[i*n+j] =Jij;
            J[j*n+i] =Jij;
        }
    }
    //printf("Jmax %g \n", Jmax);
    //exit(0);
}

class QEq{ public:
    int     n    = 0;
    //Vec3d*  ps = 0;
    double* J    = 0;
    double* qs   = 0;
    double* fqs  = 0;
    double* vqs  = 0;

    double* affins = 0;
    double* hards  = 0;
    bool*   constrain = 0; 

    double Qtarget = 0.0;
    double Qtot    = 0.0;

    void realloc(int n_){
        n=n_;
        //_realloc( J, n*n );
        _realloc( qs , n );
        _realloc( fqs, n );
        _realloc( vqs, n );
        _realloc( hards,  n );
        _realloc( affins, n );
        _realloc( constrain, n ); for(int i=0; i<n;i++){ constrain[i]=false; }
    }

    void init(){ for(int i=0; i<n; i++){ qs[i]=0; fqs[i]=0; vqs[i]=0; } }
    void constrainTypes(const int* atypes,int iconst ){ for(int i=0; i<n; i++){ if(atypes[i]==iconst){constrain[i]=true;}; }  }

    double getQvars(){
        double err2=0;
        //Qtot         = 0.0;
        double fqtot = 0.0;
        int nsum=0;
        for(int i=0; i<n; i++){
            if(constrain[i]) continue;
            double qi = qs[i];
            //Qtot += qi;
            double fq = affins[i] + hards[i]*qi;
            for(int j=0; j<n; j++){
                fq += J[i*n+j]*qs[j];
            }
            //fq *= -1;   // invert sign of force 
            //printf( "QEq.getQvars[%i] fq %g \n", i, fq );
            fqtot +=fq; 
            fqs[i]=fq;
            nsum++;
            //err2 += fq*fq;
        }
        // constrain TOTAL CHARGE
        //fqtot*=(Qtot-Qtarget)/n;
        //double invn = 1/n;
        double dfqtot= fqtot/nsum;
        for(int i=0; i<n; i++){
            if(constrain[i]) continue;
            fqs[i] -= dfqtot; 
            err2   += fqs[i]*fqs[i]; 
        };
        //printf( "Qtot %g F2 %g fqtot \n", Qtot, err2, fqtot );
        //for(int i=0; i<n; i++){ qs[i]-=qi; };
        return err2;
    }

    double moveMDdamp(double dt, double damping){
        Qtot = 0.0;
        double damp=1-damping;
        int nsum=0;
        double cvf=0;
        int nconstr=0;
        for(int i=0; i<n; i++){
            double qi  = qs[i];
            if(constrain[i]){ Qtot+=qi; continue; }
            double fqi = fqs[i];
            double vqi = vqs[i];
            vqi*=damp;
            vqi-=fqi*dt;
            qi +=vqi*dt;
            //vqs[i]  = vqs[i]*damp - fqs[i]*dt;
            // qs[i] += vqs[i]*dt;
            vqs[i] = vqi;
            qs[i]  = qi;
            Qtot  += qi;
            cvf   += fqi*vqi;
            nsum++;
        }
        //if(cvf<0){ for(int i=0; i<n; i++){ vqs[i]=0; }; }
        // force Qtarget
        //printf( "Qtot %g \n" );
        double dQ    = (Qtarget-Qtot)/(nsum);
        for(int i=0; i<n; i++){ 
            if(constrain[i]){ continue; }
            qs[i] += dQ; 
        }
        return Qtot;
    }

    double relaxChargeMD( Vec3d* ps, int nsteps=1000, double Fconf=1e-2, double dt=0.1, double damp=0.0, bool bVerbose=false, bool bInit=true ){
        //printf( "QEq.relaxChargeMD() \n" );
        double F2conf=Fconf*Fconf;
        J = new double[n*n];
        makeCoulombMatrix( n, ps, J );
        double F2=1.0;
        if(bInit){
            init();
        }else{
            for(int i=0; i<n; i++){ qs[i]*=-1.0; }; // invert charges (e-) convention
        }
        //printf( "q[-1](start) %g\n", qs[n-1] );
        for(int itr=0; itr<nsteps; itr++){
            F2 = getQvars();
            if(F2<F2conf) break;
            if(bVerbose)printf( "QEq.relaxChargeMD()[%i] |F|=%g \n", itr, sqrt(F2) );
            //printf( "F2 %g \n", F2 );
            double cvf = moveMDdamp(dt, damp);
            //printf( "q[-1](itr=%i) %g\n", qs[n-1], itr );
            if(cvf<0){ 
                if(bVerbose)printf( "QEq.relaxChargeMD()[%i](cvf(%g)<0) => v[:]=0 \n", itr, cvf );
                for(int i=0; i<n; i++){ vqs[i]=0; }; 
            }
        }
        double Qtot=0;
        for(int i=0; i<n; i++){ qs[i]*=-1.0; Qtot+=qs[i]; }; // invert charges (e-) convention
        if( fabs(Qtot-Qtarget) > 1e-6 ){ printf( "QEq.relaxChargeMD() WARNING: Qtot(%g)-Qtarget(%g)=%g \n", Qtot, Qtarget, Qtot-Qtarget ); }
        // for(int i=0; i<n; i++){ 
        //     if(bVerbose)printf( "QEq.relaxChargeMD() qs[%i] %g \n", i, qs[i] );
        // };   // reverse charges (e-)
        //printf( "q[-1](end) %g\n", qs[n-1] );
        delete [] J; J=0;
        return F2;
    }
};


#endif
