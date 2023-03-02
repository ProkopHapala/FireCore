
#ifndef FitREQ_h
#define FitREQ_h

//#include <vector>
#include "Vec3.h"
#include "NBFF.h"

void savexyz(const char* fname, AtomicSystem* A, AtomicSystem* B, const char* comment=0, const char* mode="w" ){
    FILE* fout = fopen( fname, mode);
    if(fout==0){ printf("cannot open '%s' \n", fname ); exit(0);}
    int n=0;
    if(A)n+=A->n;
    if(B)n+=B->n;
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
    //NBsystem* nbff;
    int nDOFs=0,ntype=0,nbatch=0,n0=0,n1=0;
    Vec3d*    typeREQs =0;   // [ntype] parameters for each type
    Vec3d*    typeREQsMin=0; // [ntype] equlibirum value of parameters for regularization 
    Vec3d*    typeREQsMax=0; // [ntype] equlibirum value of parameters for regularization 
    Vec3d*    typeREQs0=0;   // [ntype] equlibirum value of parameters for regularization 
    Vec3d*    typeKreg =0;   // [ntype] regulatization stiffness
    Vec3i*    typToREQ =0;   // [ntype] map each unique atom type to place in DOFs;
    
    double*   DOFs =0;       // [nDOFs]
    double*   fDOFs=0;       // [nDOFs]
    double*   vDOFs=0;       // [nDOFs]

    double  Kneutral = 1.0;
    //double max_step=-1.;
    double  max_step=0.01345;

    //double* Rs =0,Es =0,Qs =0;
    //double* fRs=0,fEs=0,fQs=0;

    //std::vector(AtomicSystem) examples; // Training examples
    AtomicSystem* batch=0;     // [nbatch] 
    double*       weights = 0; // [nbatch] scaling importaince of parameters
    double*       Es      = 0; 
    Mat3d*        poses   = 0; // [nbatch]
   
    // for rigid fitting
    AtomicSystem* systemTest0=0; //[1]
    AtomicSystem* systemTest =0; //[1]
    
    // system 0
    AtomicSystem* system0=0;   // [1]
    //Vec3d* REQ0=0;    // [system0.n] REQparams for system0
    //Vec3i* jnds =0;   // [system0.n] typToREQ  for system0
    //Vec3d* fs0  =0;   // [system0.n] typToREQ  for system0

    // Temporary array for accumulation of derivs
    int     nmax = 0;
    Vec3d * fs   = 0; //[nmax]
    //std::vector<Vec3d> fs;

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
    for(int i=0; i<nbatch; i++){ int ni = batch[i].n; if(ni>n){ n=ni;} }
    if(n>nmax){ _realloc( fs, n );  nmax=n; };
}

void tryRealocTemp_rigid(){
    if(nmax<systemTest0->n){
        nmax=systemTest0->n;
        _realloc( fs, nmax );
    }
}

int init_types( int ntype_, Vec3i* typeMask, Vec3d* tREQs=0, bool bCopy=false ){
    printf( "init_types(%i) \n", ntype_ );
    int nDOFs=0;
    typToREQ = new Vec3i[ntype];
    //tREQs    = new Vec3d[ntyp];
    for(int i=0; i<ntype_; i++){
        const Vec3i& tm=typeMask[i];
        Vec3i&       tt=typToREQ[i];
        for(int j=0; j<3; j++){
            if(tm.array[j]){
                tt.array[j]=nDOFs;
                nDOFs++;
            }
            else{
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
            typeREQs    = new Vec3d[ntype]; for(int i=0; i<ntype; i++){ typeREQs[i]    = tREQs[i]; }
            typeREQsMin = new Vec3d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMin[i] = Vec3dmin; }
            typeREQsMax = new Vec3d[ntype]; for(int i=0; i<ntype; i++){ typeREQsMax[i] = Vec3dmax; }
            typeREQs0   = new Vec3d[ntype]; for(int i=0; i<ntype; i++){ typeREQs0[i]   = tREQs[i]; }
            typeKreg    = new Vec3d[ntype];  for(int i=0; i<ntype; i++){ typeKreg[i]   = Vec3dZero; }
        }else{ typeREQs = tREQs; }

        DOFsFromTypes(); 
    }
    //printf(" DOFs=");for(int j=0;j<nDOFs;j++){ printf("%g ",DOFs[j]); };printf("\n");
    return nDOFs;
}

void setSystem( int isys, int na, int* types, Vec3d* ps, bool bCopy=false ){
    printf( "setSystem(%i) \n", isys );
    AtomicSystem** pS;
    AtomicSystem*  S;
    if(isys<0){
        if     (isys==-3){ pS=&systemTest;  n1=na; }
        else if(isys==-2){ pS=&systemTest0; n1=na; }
        else if(isys==-1){ pS=&system0;     n0=na; }
        if(*pS==0   ){ *pS = new AtomicSystem(na); }
        S=*pS;
    }else             { S=&batch[isys]; }
    if(bCopy){
        if(S->n!=na){ printf("ERORRO FitREQ.setSystem() na(%i)!=S.n(%i) \n", na, S->n ); exit(0); }
        for(int i=0; i<na; i++){ S->types[i]=types[i]; S->ps[i]=ps[i]; }
    }else{
        S->n     =na;
        S->types =types;
        S->ps    =ps;
    }
}

void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy=false, bool bAlloc=false ){
    printf( "setRigidSamples() \n" );
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

double evalExampleDerivs_LJQ(int n, int* types, Vec3d* ps ){
    const double COULOMB_CONST_ = 14.3996448915;  //  [V*A/e] = [ (eV/A) * A^2 /e^2]
    int    nj   =system0->n;
    int*   jtyp =system0->types;
    Vec3d* jpos =system0->ps;
    double Etot=0;
    double Qtot=0;
    for(int i=0; i<n; i++){
        int   ti          = types[i];
        const Vec3d& REQi = typeREQs[ti];
        const Vec3d& pi   = ps[i]; 
        Vec3d fsi         = Vec3dZero;
        Qtot+=REQi.z;
        for(int j=0; j<nj; j++){
            int tj              = jtyp[i];
            const Vec3d& REQj   = typeREQs[tj];
            //const Vec3d& REQj = REQ0[j]; // optimization
            Vec3d d             = jpos[j] - pi;
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.z*REQj.z;

            //printf( "ij[%i,%i] REQij(%g,%g,%g) REQi(%g,%g,%g) REQj(%g,%g,%g)\n", i,j,  R,E0,Q,   REQi.x,REQi.y,REQi.z,   REQj.x,REQj.y,REQj.z );

            // --- Eectrostatic
            double ir2     = 1/( d.norm2() + 1e-4 );
            double ir      = sqrt(ir2);
            double dEel_dQ = ir * COULOMB_CONST_;
            double Eel     = Q*dEel_dQ;
            // --- Lenard-Jones
            double u2   = ir2*(R*R);
            double u4   = u2*u2;
            double u6   = u4*u2;

            // ELJ      = E0*( (R/r)^12 - 2*(R/r)^6 )
            // dELJ/dR  = E0*( 12*(R/r)^11/r - 12*(R/r)^5/r    )
            double dELJ_dE = u6   *( u6 - 2 );
            double ELJ     = E0*dELJ_dE;
            // ---- Total E
            Etot  += ELJ + Eel;

            double dELJ_dR =  E0*12*ir*( u6*u4 - u4 );
            fsi.x += dELJ_dR;         // dEtot/dRi
            fsi.y += dELJ_dE*REQj.y;  // dEtot/dEi
            fsi.z += dEel_dQ*REQj.z;  // dEtot/dQi
        }
        fs[i].add(fsi);
    }

    return Etot;
}

double evalExampleDerivs_Qneutral(int n, int* types, Vec3d* ps, double Qtot ){
    // Qtot = Sum_i * Q_i
    // E = K*Qtot^2
    // d E^2 / dqi = K*Qtot*(  )
    double E = 0;
    double fQ = -Kneutral*Qtot;
    for(int i=0; i<n; i++){
        //int   ti          = types[i];
        //const Vec3d& REQi = typeREQs[ti]; 
        fs[i].z += fQ;
        E       += Kneutral*Qtot*Qtot;
    }
    return  E;
}

void DOFsToTypes(){
    for(int i=0; i<ntype; i++ ){
        const Vec3i& tt= typToREQ[i];
        Vec3d& REQ     = typeREQs[i];
        if(tt.x>=0)REQ.x = DOFs[tt.x];
        if(tt.y>=0)REQ.y = DOFs[tt.y];
        if(tt.z>=0)REQ.z = DOFs[tt.z];
    }
}

void DOFsFromTypes(){
    for(int i=0; i<ntype; i++ ){
        const Vec3i& tt  = typToREQ[i];
        const Vec3d& REQ = typeREQs[i];
        if(tt.x>=0)DOFs[tt.x] = REQ.x;
        if(tt.y>=0)DOFs[tt.y] = REQ.y;
        if(tt.z>=0)DOFs[tt.z] = REQ.z;
    }
}

void acumDerivs( int n, int* types, double dE){
    for(int i=0; i<n; i++){
        int t           = types[i];
        const Vec3i& tt = typToREQ[t];
        const Vec3d  f  = fs[i];
        //printf( "acumDerivs[%i] t %i tt(%i,%i,%i) f(%g,%g,%g)\n", i, t, tt.x,tt.y,tt.z,  f.x,f.y,f.z );
        if(tt.x>=0)fDOFs[tt.x]+=f.x*dE;
        if(tt.y>=0)fDOFs[tt.y]+=f.y*dE;
        if(tt.z>=0)fDOFs[tt.z]+=f.z*dE;
    }
}

void regularization_force(){
    for(int i=0; i<ntype; i++ ){
        const Vec3i& tt   = typToREQ[i];
        const Vec3d& REQ  = typeREQs [i];
        const Vec3d& REQ0 = typeREQs0[i];
        const Vec3d& K    = typeKreg [i];
        if(tt.x>=0)fDOFs[tt.x] = (REQ0.x-REQ.x)*K.x;
        if(tt.y>=0)fDOFs[tt.y] = (REQ0.y-REQ.y)*K.y;
        if(tt.z>=0)fDOFs[tt.z] = (REQ0.z-REQ.z)*K.z;
    }
}

void limit_params(){
    for(int i=0; i<ntype; i++ ){
        const Vec3i& tt     = typToREQ[i];
        const Vec3d& REQ    = typeREQs [i];
        const Vec3d& REQmin = typeREQsMin[i];
        const Vec3d& REQmax = typeREQsMax[i];
        if(tt.x>=0){ fDOFs[tt.x]=_clamp(fDOFs[tt.x],REQmin.x,REQmax.x ); }
        if(tt.y>=0){ fDOFs[tt.y]=_clamp(fDOFs[tt.y],REQmin.y,REQmax.y ); }
        if(tt.z>=0){ fDOFs[tt.z]=_clamp(fDOFs[tt.z],REQmin.z,REQmax.z ); }
    }
}

void renormWeights(double R){
    double s=0; 
    for(int i=0; i<nbatch; i++){ s+=weights[i]; }
    s=R/s; 
    for(int i=0; i<nbatch; i++){ weights[i]*=s; }
}

void clean_fs(int n){
    for(int i=0; i<n; i++){ fs[i]=Vec3dZero; }
}

double evalDerivs( double* Eout=0 ){
    tryRealocTemp();
    double Error = 0;
    for(int i=0; i<nbatch; i++){
        const AtomicSystem& C = batch[i];
        double Eref=Es[i];
        double wi = 1; 
        if(weights) wi = weights[i];
        // switch(ikind){ case ikind_LJQ:
        clean_fs(C.n);
        double E  = evalExampleDerivs_LJQ(C.n, C.types, C.ps );  // Forward pass
        if(Eout){ Eout[i]=E; };
        double dE = (E - Eref)*wi;
        Error += dE;
        acumDerivs(C.n, C.types, dE );
        //double dE = C.E - E;
        //double E_ = evalExampleDerivs_LJQ(C.n, C.types, C.ps, dE*wi );  // Backward pass
    }
    return Error;
}

double evalDerivsRigid( double* Eout=0 ){
    //printf( "evalDerivsRigid() \n" );
    tryRealocTemp_rigid();
    const AtomicSystem& C0 = *systemTest0;
    const AtomicSystem& C  = *systemTest;
    double Error = 0; 
    for(int i=0; i<nbatch; i++){
        double Eref=Es[i]; 
        rigid_transform( poses[i].a, 0, poses[i].c, poses[i].b, C.n, C0.ps, C.ps );
        //savexyz("FitREQ_debug.xyz", system0, systemTest,0, "a" );
        clean_fs(C.n);
        double E  = evalExampleDerivs_LJQ( C.n, C.types, C.ps );  // Forward pass
        if(Eout){ Eout[i]=E; };
        double wi = 1;
        if(weights) wi = weights[i];
        //if(bLimitForce){ };
        double dE = -(E - Eref)*wi;
        Error += dE*dE;
        acumDerivs(C.n, C.types, dE );
        //printf( "[%i] x %g E %g Eref %g ", i, poses[i].a.x, E, Eref );printf("}\n");
        //printf("fs={");for(int j=0;j<C.n;j++){ printf("(%g,%g,%g)",fs[j].x,fs[j].y,fs[j].z); };//printf("}\n");
        //printf("fDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",fDOFs[j]); };printf("}\n");
    }
    //printf("fDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",fDOFs[j]); };printf("}\n");
    //printf("vDOFs={");for(int j=0;j<nDOFs;j++){ printf("%g,",vDOFs[j]); };printf("}\n");
    return Error;
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
