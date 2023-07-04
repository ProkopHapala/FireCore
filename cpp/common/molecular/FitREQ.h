
#ifndef FitREQ_h
#define FitREQ_h

#include <vector>
#include "Vec3.h"
#include "quaternion.h"
//#include "NBFF.h"
#include "Atoms.h"
#include "MMFFparams.h"
#include "Forces.h"


void savexyz(const char* fname, Atoms* A, Atoms* B, const char* comment=0, const char* mode="w" ){
    FILE* fout = fopen( fname, mode);
    if(fout==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    int n=0;
    if(A)n+=A->natoms;
    if(B)n+=B->natoms;
    fprintf(fout,"%i\n", n);
    if(comment){ fprintf(fout,"%s\n", comment ); }else{ fprintf(fout, "comment\n", n); };
    if(A)A->atomsToXYZ(fout);
    if(B)B->atomsToXYZ(fout);
    fclose(fout);
}

void rigid_transform( Vec3d shift, Vec3d* unshift, Vec3d dir, Vec3d up, int n, Vec3d* pin, Vec3d* pout ){
    Mat3d M;
    M.c=dir;
    M.b=up;
    M.a.set_cross( dir, up );
    for(int i=0; i<n; i++){
        Vec3d p,p0=pin[i];
        if(unshift)p0.sub(*unshift);
        M.dot_to_T( p0, p );
        p.add(shift);
        pout[i]=p;
    }
}

class FitREQ{ public:
    //NBFF* nbff;
    int nDOFs=0,ntype=0,nbatch=0,n0=0,n1=0;
    int imodel=1;
    Quat4d*    typeREQs =0;   // [ntype] parameters for each type
    Quat4d*    typeREQsMin=0; // [ntype] equlibirum value of parameters for regularization 
    Quat4d*    typeREQsMax=0; // [ntype] equlibirum value of parameters for regularization 
    Quat4d*    typeREQs0=0;   // [ntype] equlibirum value of parameters for regularization 
    Quat4d*    typeKreg =0;   // [ntype] regulatization stiffness
    Quat4i*    typToREQ =0;   // [ntype] map each unique atom type to place in DOFs;
    
    double*   DOFs =0;       // [nDOFs]
    double*   fDOFs=0;       // [nDOFs]
    double*   vDOFs=0;       // [nDOFs]

    double  Kneutral = 1.0;
    //double max_step=-1.;
    double  max_step=0.01345;

    //double* Rs =0,Es =0,Qs =0;
    //double* fRs=0,fEs=0,fQs=0;

    //std::vector(Atoms) examples; // Training examples
    Atoms*   batch=0;     // [nbatch]  // ToDo: would be more convenient to store Atoms* rather than Atoms
    double*  weights = 0; // [nbatch] scaling importaince of parameters
    double*  Es      = 0; 
    Mat3d*   poses   = 0; // [nbatch]
   
    // for rigid fitting
    Atoms* systemTest0=0; //[1]
    Atoms* systemTest =0; //[1]
    
    // system 0
    Atoms* system0=0;   // [1]

    // Temporary array for accumulation of derivs
    int     nmax = 0;
    Quat4d* fs   = 0; //[nmax]
    //std::vector<Vec3d> fs;

    std::vector<Atoms> batch_vec;    // ToDo: would be more convenient to store Atoms* rather than Atoms

    MMFFparams* params=0; 

    // ------- Arrays for decomposition of energy components
    bool bDecomp = false;
    Quat4d** Elines  = 0;
    Quat4d** Eilines = 0;
    Quat4d** Ejlines = 0;


//void realoc( int nR=0, int nE=0, int nQ=0 ){
void realloc( int nDOFs_ ){
    //nDOFs=nR+nE+nQ;
    nDOFs=nDOFs_;
    _realloc(  DOFs , nDOFs );
    _realloc( fDOFs , nDOFs );
    _realloc( vDOFs , nDOFs ); for(int i=0;i<nDOFs;i++){vDOFs[i]=0;}
    //Rs=DOFs;Es=DOFs+nR;Qs=DOFs+nR+nE;
    //fRs=fDOFs;fEs=fDOFs+nR;fQs=fDOFs+nR+nE;
}

void tryRealocTemp(){
    int n=nmax;
    for(int i=0; i<nbatch; i++){ int ni = batch[i].natoms; if(ni>n){ n=ni;} }
    if(n>nmax){ _realloc( fs, n );  nmax=n; };
}

void tryRealocTemp_rigid(){
    if(nmax<systemTest0->natoms){
        nmax=systemTest0->natoms;
        _realloc( fs, nmax );
    }
}

int init_types( int ntype_, Quat4i* typeMask, Quat4d* tREQs=0, bool bCopy=false ){
    printf( "FitREQ::init_types() ntype_=%i ntype=%i \n", ntype_, ntype );
    int nDOFs=0;
    typToREQ = new Quat4i[ntype_];
    for(int i=0; i<ntype_; i++){
        const Quat4i& tm=typeMask[i];
        Quat4i&       tt=typToREQ[i];
        for(int j=0; j<4; j++){
            if(tm.array[j]){
                tt.array[j]=nDOFs;
                nDOFs++;
            }else{
                tt.array[j]=-1;
            }
        }
        //printf(  "tt[%i] tm(%i,%i,%i) tt(%i,%i,%i)  nDOF %i \n", i, tm.x,tm.y,tm.z, tt.x,tt.y,tt.z, nDOFs);
    }
    realloc(nDOFs);
    //poses = new Mat3d[nbatch];
    //if( bRigid){
    //    poses = new Mat3d[nbatch];
    //}
    if(tREQs){
        ntype = ntype_; 
        if(bCopy){
            typeREQs    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs[i]    = tREQs[i]; }
            typeREQsMin = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Quat4dmin; }
            typeREQsMax = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Quat4dmax; }
            typeREQs0   = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0[i]   = tREQs[i]; }
            typeKreg    = new Quat4d[ntype]; for(int i=0; i<ntype; i++){ typeKreg[i]    = Quat4dZero; }
        }else{ typeREQs = tREQs; }
        DOFsFromTypes(); 
    }
    //printf(" DOFs=");for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
    return nDOFs;
}

void setSystem( int isys, int na, int* types, Vec3d* ps, bool bCopy=false ){
    printf( "FitREQ::setSystem(%i) \n", isys );
    Atoms** pS;
    Atoms*  S;
    if(isys<0){
        if     (isys==-3){ pS=&systemTest;  n1=na; }
        else if(isys==-2){ pS=&systemTest0; n1=na; }
        else if(isys==-1){ pS=&system0;     n0=na; }
        if(*pS==0   ){ *pS = new Atoms(na); }
        S=*pS;
    }else             { S=&batch[isys]; }
    if(bCopy){
        if(S->natoms!=na){ printf("ERORRO FitREQ.setSystem() na(%i)!=S.n(%i) \n", na, S->natoms ); exit(0); }
        for(int i=0; i<na; i++){ S->atypes[i]=types[i]; S->apos[i]=ps[i]; }
    }else{
        S->natoms =na;
        S->atypes =types;
        S->apos   =ps;
    }
}

void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy=false, bool bAlloc=false ){
    printf( "FitREQ::setRigidSamples() \n" );
    nbatch = n; 
    if(bCopy){
        if(Es_   ){Es    = new double[nbatch]; }
        if(poses_){poses = new Mat3d [nbatch]; }
        if(Es_   )for(int i=0;i<nbatch; i++){ Es   [i]=Es_   [i]; }
        if(poses_)for(int i=0;i<nbatch; i++){ poses[i]=poses_[i]; }
    }else{
        if(Es_   )Es = Es_; 
        if(poses_)poses = poses_;
    }
}

double evalExampleDerivs_LJQH(int n, int* types, Vec3d* ps ){
    int    nj   =system0->natoms;
    int*   jtyp =system0->atypes;
    Vec3d* jpos =system0->apos;
    double Etot=0;
    double Qtot=0;
    for(int i=0; i<n; i++){
        int   ti          = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        Quat4d fsi         = Quat4dZero;
        Qtot+=REQi.z;
        for(int j=0; j<nj; j++){
            int tj              = jtyp[j];
            const Quat4d& REQj  = typeREQs[tj];
            //const Quat4d& REQj = REQ0[j]; // optimization
            Vec3d d             = jpos[j] - pi;
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.z*REQj.z;
            double H  = REQi.w*REQj.w;

            //printf( "ij[%i=%i,%i=%i] REQij(%6.3f,%10.7f,%6.3f,%6.3f) test:REQi(%6.3f,%10.7f,%6.3f,%6.3f) sys0:REQj(%6.3f,%10.7f,%6.3f,%6.3f)\n", i,ti, j,tj,  R,E0,Q,H,   REQi.x,REQi.y,REQi.z,REQi.w,   REQj.x,REQj.y,REQj.z,REQj.w );

            // --- Eectrostatic
            //double ir2     = 1/( d.norm2() + 1e-4 );
            double ir2     = 1/( d.norm2() );
            double ir      = sqrt(ir2);
            double dE_dQ   = ir * COULOMB_CONST;
            double Eel     = Q*dE_dQ;
            // --- Lenard-Jones
            double u2   = ir2*(R*R);
            double u4   = u2*u2;
            double u6   = u4*u2;

            // ELJ      = E0*( (R/r)^12 - 2*(R/r)^6 )
            // dELJ/dR  = E0*( 12*(R/r)^11/r - 12*(R/r)^5/r    )

            double dE_dE = u6   *( u6 - 2 );
            double dE_dH = u6;
            double dE_dR = E0*12*ir*( u6*u4 - u4 );
            double ELJ   = E0*dE_dE + H*dE_dH;

            Etot  += ELJ + Eel;

            fsi.x += dE_dR;        // dEtot/dRi
            fsi.y += dE_dE*REQj.y; // dEtot/dEi
            fsi.z += dE_dQ*REQj.z; // dEtot/dQi
            fsi.w += dE_dH*REQj.w; // dEtot/dHi
        }
        fs[i].add(fsi);
    }

    return Etot;
}

double evalExampleDerivs_LJQH2(int n, int* types, Vec3d* ps ){
    //printf( "evalExampleDerivs_LJQH2()\n" );
    int    nj   =system0->natoms;
    int*   jtyp =system0->atypes;
    Vec3d* jpos =system0->apos;
    double Etot=0;
    double Qtot=0;
    for(int i=0; i<n; i++){
        int   ti          = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        Quat4d fsi         = Quat4dZero;
        Qtot+=REQi.z;
        for(int j=0; j<nj; j++){
            int tj              = jtyp[j];
            const Quat4d& REQj  = typeREQs[tj];
            //const Quat4d& REQj = REQ0[j]; // optimization
            Vec3d d             = jpos[j] - pi;
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.z*REQj.z;
            double H  = REQi.w*REQj.w;

            if(H>0) H=0;

            //printf( "ij[%i=%i,%i=%i] REQij(%6.3f,%10.7f,%6.3f,%6.3f) test:REQi(%6.3f,%10.7f,%6.3f,%6.3f) sys0:REQj(%6.3f,%10.7f,%6.3f,%6.3f)\n", i,ti, j,tj,  R,E0,Q,H,   REQi.x,REQi.y,REQi.z,REQi.w,   REQj.x,REQj.y,REQj.z,REQj.w );

            // --- Eectrostatic
            //double ir2     = 1/( d.norm2() + 1e-4 );
            double ir2     = 1/( d.norm2() );
            double ir      = sqrt(ir2);
            double dE_dQ   = ir * COULOMB_CONST;
            double Eel     = Q*dE_dQ;
            // --- Lenard-Jones
            double u2   = ir2*(R*R);
            double u4   = u2*u2;
            double u6   = u4*u2;

            //if( (i==0) && (j==nj) ){ printf( "r %g \n", 1/ir ); }
            //if( i==0 ){ printf( "[%i,%i] r %g \n", i,j, 1/ir ); }


            // ELJ      = E0*( (R/r)^12 - 2*(R/r)^6 )
            // dELJ/dR  = E0*( 12*(R/r)^11/r - 12*(R/r)^5/r    )

            double dE_dE = u6   *( u6 - 2 );
            double dE_dH = u6*u6;
            double dE_dR = E0*12*ir*( u6*u4 - u4 );
            double ELJ   = E0*dE_dE + H*dE_dH;

            //printf( "ij[%i=%i,%i=%i] EH %g EPaul %g EvdW %g Eel %g REQij(%6.3f,%10.7f,%6.3f,%6.3f) \n", i,ti, j,tj,  H*u6*u6, E0*u6*u6, E0*u6*-2, Q*dE_dQ,  R,E0,Q,H  );

            Etot  += ELJ + Eel;

            fsi.x += dE_dR;        // dEtot/dRi
            fsi.y += dE_dE*REQj.y; // dEtot/dEi
            fsi.z += dE_dQ*REQj.z; // dEtot/dQi
            fsi.w += dE_dH*REQj.w; // dEtot/dHi
        }
        fs[i].add(fsi);
    }

    return Etot;
}

Quat4d evalExampleEnergyComponents_LJQH2(int ni, int* types, Vec3d* ps, int isamp ){
    const int    nj   = system0->natoms;
    const int*   jtyp = system0->atypes;
    const Vec3d* jpos = system0->apos;
    Quat4d Esum = Quat4dZero;

    if(Eilines)for(int i=0; i<nj; i++){ Quat4d* Eiline=Eilines[i]; if(Eiline){ Eiline[isamp]=Quat4dZero; } }
    if(Ejlines)for(int j=0; j<nj; j++){ Quat4d* Ejline=Ejlines[j]; if(Ejline){ Ejline[isamp]=Quat4dZero; } }

    for(int i=0; i<ni; i++){
        const int     ti   = types[i];
        const Vec3d&  pi   = ps[i]; 
        const Quat4d& REQi = typeREQs[ti];
        Quat4d* Eiline=0; if(Eilines){ Eiline=Eilines[i];  }
        for(int j=0; j<nj; j++){
            const int tj        = jtyp     [j];
            const Quat4d& REQj  = typeREQs [tj];
            //const Quat4d  REQij = _mixREQ(REQi,REQj);
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.z*REQj.z;
            double H  = REQi.w*REQj.w;

            const Vec3d d       = jpos     [j] - pi;

            if(H>0) H=0;
            //printf( "ij[%i=%i,%i=%i] REQij(%6.3f,%10.7f,%6.3f,%6.3f) test:REQi(%6.3f,%10.7f,%6.3f,%6.3f) sys0:REQj(%6.3f,%10.7f,%6.3f,%6.3f)\n", i,ti, j,tj,  R,E0,Q,H,   REQi.x,REQi.y,REQi.z,REQi.w,   REQj.x,REQj.y,REQj.z,REQj.w );
            // --- Eectrostatic
            //double ir2     = 1/( d.norm2() + 1e-4 );
            double ir2     = 1/( d.norm2() );
            double ir      = sqrt(ir2);
            // --- Lenard-Jones
            double u2   = ir2*(R*R);
            double u6   = u2*u2*u2;
            double u12  = u6*u6;
            Quat4d Eij = Quat4d{
                E0   *u12,                   // Pauli
                E0*-2*u6 ,                   // London
                Q *COULOMB_CONST*sqrt(ir2), // Coulomb
                E0*u12,                      // Hbond 
            };
            Esum.add(Eij);
            if( bDecomp ){
                {                                       if(Eiline){ Eiline[isamp].add(Eij); } }
                if(Ejlines){ Quat4d* Ejline=Ejlines[j]; if(Ejline){ Ejline[isamp].add(Eij); } }
                if(Elines){
                    int ij = i + j*ni;
                    Quat4d* Eline = Elines[ij];
                    if( Eline ){    Eline[isamp]=Eij;  }
                }
            }
        }
    }
    return Esum;
}

double evalExampleDerivs_Qneutral(int n, int* types, Vec3d* ps, double Qtot ){
    // Qtot = Sum_i * Q_i
    // E = K*Qtot^2
    // d E^2 / dqi = K*Qtot*(  )
    double E = 0;
    double fQ = -Kneutral*Qtot;
    for(int i=0; i<n; i++){
        //int   ti          = types[i];
        //const Quat4d& REQi = typeREQs[ti]; 
        fs[i].z += fQ;
        E       += Kneutral*Qtot*Qtot;
    }
    return  E;
}

inline void DOFsToType(int i){
    const Quat4i& tt = typToREQ[i];
    Quat4d& REQ      = typeREQs[i];
    if(tt.x>=0)REQ.x = DOFs[tt.x];
    if(tt.y>=0)REQ.y = DOFs[tt.y];
    if(tt.z>=0)REQ.z = DOFs[tt.z];
    if(tt.w>=0)REQ.w = DOFs[tt.w];
}
void DOFsToTypes(){ for(int i=0; i<ntype; i++ ){ DOFsToType(i); } }
void getType(int i, Quat4d& REQ ){ typeREQs[i]=REQ; DOFsToType(i); }


inline void DOFsFromType(int i){
    const Quat4i& tt  = typToREQ[i];
    const Quat4d& REQ = typeREQs[i];
    if(tt.x>=0)DOFs[tt.x] = REQ.x;
    if(tt.y>=0)DOFs[tt.y] = REQ.y;
    if(tt.z>=0)DOFs[tt.z] = REQ.z;
    if(tt.w>=0)DOFs[tt.w] = REQ.w;
}
void DOFsFromTypes(){ for(int i=0; i<ntype; i++ ){ DOFsFromType(i); } }
void setType(int i, Quat4d REQ ){ typeREQs[i]=REQ; DOFsFromType(i); }

void acumDerivs( int n, int* types, double dE){
    for(int i=0; i<n; i++){
        int t            = types[i];
        const Quat4i& tt = typToREQ[t];
        const Quat4d  f  = fs[i];
        //printf( "acumDerivs[%i] t %i tt(%i,%i,%i) f(%g,%g,%g)\n", i, t, tt.x,tt.y,tt.z,  f.x,f.y,f.z );
        if(tt.x>=0)fDOFs[tt.x]+=f.x*dE;
        if(tt.y>=0)fDOFs[tt.y]+=f.y*dE;
        if(tt.z>=0)fDOFs[tt.z]+=f.z*dE;
        if(tt.w>=0)fDOFs[tt.w]+=f.w*dE;
    }
}

void regularization_force(){
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt   = typToREQ [i];
        const Quat4d& REQ  = typeREQs [i];
        const Quat4d& REQ0 = typeREQs0[i];
        const Quat4d& K    = typeKreg [i];
        if(tt.x>=0)fDOFs[tt.x] = (REQ0.x-REQ.x)*K.x;
        if(tt.y>=0)fDOFs[tt.y] = (REQ0.y-REQ.y)*K.y;
        if(tt.z>=0)fDOFs[tt.z] = (REQ0.z-REQ.z)*K.z;
        if(tt.w>=0)fDOFs[tt.w] = (REQ0.w-REQ.w)*K.w;
    }
}

void limit_params(){
    for(int i=0; i<ntype; i++ ){
        const Quat4i& tt     = typToREQ[i];
        const Quat4d& REQ    = typeREQs [i];
        const Quat4d& REQmin = typeREQsMin[i];
        const Quat4d& REQmax = typeREQsMax[i];
        if(tt.x>=0){ fDOFs[tt.x]=_clamp(fDOFs[tt.x],REQmin.x,REQmax.x ); }
        if(tt.y>=0){ fDOFs[tt.y]=_clamp(fDOFs[tt.y],REQmin.y,REQmax.y ); }
        if(tt.z>=0){ fDOFs[tt.z]=_clamp(fDOFs[tt.z],REQmin.z,REQmax.z ); }
        if(tt.w>=0){ fDOFs[tt.w]=_clamp(fDOFs[tt.w],REQmin.w,REQmax.w ); }
    }
}

void renormWeights(double R){
    double s=0; 
    for(int i=0; i<nbatch; i++){ s+=weights[i]; }
    s=R/s; 
    for(int i=0; i<nbatch; i++){ weights[i]*=s; }
}

void clean_fs(int n){
    for(int i=0; i<n; i++){ fs[i]=Quat4dZero; }
}

double evalDerivs( double* Eout=0 ){
    printf( "FitREQ::evalDerivs() \n" );
    tryRealocTemp();
    double Error = 0;
    for(int i=0; i<nbatch; i++){
        //printf("evalDerivs[%i]\n", i );
        const Atoms& C = batch[i];
        //C.print();
        double Eref=Es[i];
        double wi = 1; 
        if(weights) wi = weights[i];
        // switch(ikind){ case ikind_LJQ:
        clean_fs(C.natoms);
        double E=0;
        switch (imodel){
            case 1: E = evalExampleDerivs_LJQH (C.natoms, C.atypes, C.apos ); break;
            case 2: E = evalExampleDerivs_LJQH2(C.natoms, C.atypes, C.apos ); break;
        }
        if(Eout){ Eout[i]=E; };
        double dE = (E - Eref)*wi;
        Error += dE;
        acumDerivs(C.natoms, C.atypes, dE );
        //double dE = C.E - E;
        //double E_ = evalExampleDerivs_LJQ(C.n, C.atypes, C.apos, dE*wi );  // Backward pass
    }
    return Error;
}

double evalDerivsRigid( double* Eout=0 ){
    //printf( "FitREQ::evalDerivsRigid() \n" );
    tryRealocTemp_rigid();
    const Atoms& C0 = *systemTest0;
    const Atoms& C  = *systemTest;
    double Error = 0; 
    for(int i=0; i<nbatch; i++){
        double Eref=Es[i]; 
        rigid_transform( poses[i].a, 0, poses[i].c, poses[i].b, C.natoms, C0.apos, C.apos );
        //savexyz("FitREQ_debug.xyz", system0, systemTest,0, "a" );
        clean_fs(C.natoms);
        double E =0;
        switch (imodel){
            case 1: E = evalExampleDerivs_LJQH (C.natoms, C.atypes, C.apos ); break;
            case 2: E = evalExampleDerivs_LJQH2(C.natoms, C.atypes, C.apos ); break;
        }
        if(Eout){ Eout[i]=E; };
        double wi = 1;
        if(weights) wi = weights[i];
        //if(bLimitForce){ };
        double dE = -(E - Eref)*wi;
        Error += dE*dE;
        acumDerivs(C.natoms, C.atypes, dE );
        //printf( "[%i] x %g E %g Eref %g ", i, poses[i].a.x, E, Eref );printf("}\n");
        //printf("fs={");for(int j=0;j<C.natoms;j++){ printf("(%g,%g,%g)",fs[j].x,fs[j].y,fs[j].z); };//printf("}\n");
        //printf("fDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",fDOFs[j]); };printf("}\n");
    }
    //printf("fDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",fDOFs[j]); };printf("}\n");
    //printf("vDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",vDOFs[j]); };printf("}\n");
    return Error;
}

int loadXYZ( const char* fname, int n0, int* i0s, int ntest, int* itests, int* types0=0, int* testtypes=0 ){
    FILE* fin = fopen( fname, "r" );
    if(fin==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    const int nline=1024;
    char line[1024];
    char at_name[8];
    int isys = 0;
    int il   = 0;
    int na   = 0; 
    Atoms atoms;
    if(!system0){ system0=new Atoms(n0); }else { system0->realloc( n0); }
    Atoms S         ( ntest );
    bool bReadTypes = !(types0 && testtypes);
    printf( "FitREQ::loadXYZ() n0 %i ntest %i \n", n0, ntest );
    while( fgets(line, nline, fin) ){
        //printf( ">>%s<<\n", line );
        //printf( "il %i isys %i na %i \n", il, isys, na, ntest,  );
        if      ( il==0 ){
            sscanf( line, "%i", &na );
            if( (na>0)&&(na<10000) ){
                if(isys==0){
                    atoms.allocNew(na);
                }
                S.allocNew(ntest);
            }else{ printf( "ERROR in FitREQ::loadXYZ() Suspicious number of atoms (%i) while reading `%s`  => Exit() \n", na, fname ); exit(0); }
        }else if( il==1 ){
            // comment - ToDo : Here we can read the reference Energy directly from .xyz 
        }else if( il<na+2 ){
            double x,y,z,q;
            //printf( ">>%s<<\n", line );
            int nret = sscanf( line, "%s %lf %lf %lf %lf", at_name, &x, &y, &z, &q );
            if(nret<5){q=0;}
            int i=il-2;
            //printf( "[%i] `%s` (%g,%g,%g) q %g \n", i, at_name, x, y, z, q );
            atoms.apos  [i].set(x,y,z);
            if(bReadTypes){ atoms.atypes[i]=params->getAtomType(at_name); }else{ atoms.atypes[i]=-1; }
        }
        il++;
        if( il >= na+2 ){ 
            //printf( "===== FitREQ::loadXYZ() isys %i na %i n0 %i ntest %i \n", isys, na, system0->natoms, S.natoms  );
            //printf( "atoms " ); atoms.print();
            if(isys==0)for(int i=0; i<n0;    i++ ){ int ia=i0s   [i];  system0->apos[i]=atoms.apos[ia]; int t; if(types0   ){ t=types0   [i]; }else{ t=atoms.atypes[ia]; }; system0->atypes[i]=t; } // store to system0
            {          for(int i=0; i<ntest; i++ ){ int ia=itests[i];         S.apos[i]=atoms.apos[ia]; int t; if(testtypes){ t=testtypes[i]; }else{ t=atoms.atypes[ia]; };  S      .atypes[i]=t; } // store to batch[isys]
            //printf( "system0 " ); system0->print();
            //printf( "stests  " ); S       .print();
            batch_vec.push_back( S ); }
            il=0; isys++; 
        }
    }
    nbatch =  batch_vec.size();
    batch  = &batch_vec[0];
    _realloc(Es,nbatch);
    atoms.dealloc();
    fclose(fin);
    return nbatch;
}


void clean_derivs(){
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0; }
}

double limit_dt(double dt){
    double fm=0;
    for(int i=0; i<nDOFs; i++){fm=_max(fm,fabs(fDOFs[i]));}
    if(dt*fm>max_step){
        dt = max_step/fm;
        printf( "limit_dt %g \n", dt );
    }
    return dt;
}

double move_GD( double dt ){
    double F2 = 0;
    if(max_step>0){ dt=limit_dt(dt);};
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        DOFs[i] += f*dt;
        F2 += f*f;
    }
    return F2;
}

double move_MD( double dt, double damp=0.1 ){
    double cdamp = 1-damp;
    double F2 = 0;
    //if(max_step>0){ dt=sqrt(limit_dt(dt));};
    for(int i=0; i<nDOFs; i++){
        double f = fDOFs[i];
        double v = vDOFs[i];
        v*=cdamp;
        v+=f*dt;
        vDOFs[i] = v;
        DOFs[i] += v*dt;
        F2      += f*f;
    }
    return F2;
}

// void atomsToXYZ(FILE* fout, int n, int* types, Vec3d* ps){
//     for(int i=0; i<n; i++){
//         fprintf( fout, "%i %20.10f %20.10f %20.10f\n", types[i], ps[i].x,ps[i].y,ps[i].z );
//     }
// }

}; // class FitREQ


#endif
