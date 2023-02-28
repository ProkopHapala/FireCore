
#ifndef FitREQ_h
#define FitREQ_h

//#include <vector>
#include "Vec3.h"
#include "NBFF.h"


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
    int ntype=0;
    Vec3d*    typeREQs=0;  // parameters for each type
    Vec3i*    typToREQ=0;  //  map each unique atom type to place in DOFs;
    int       nDOFs=0,nR=0,nE=0,nQ=0;
    double*   DOFs =0;
    double*   fDOFs=0;

    //double* Rs =0,Es =0,Qs =0;
    //double* fRs=0,fEs=0,fQs=0;

    //std::vector(AtomicSystem) examples; // Training examples
    int nbatch=0;
    AtomicSystem* batch=0;     // [nbatch] 
    double* weights = 0; // [nbatch] scaling importaince of parameters
    double* Es      = 0; 
   
    // for rigid fitting
    AtomicSystem* systemTest0=0; //[1]
    AtomicSystem* systemTest =0; //[1]
    Mat3d*        poses      =0; //[nbatch]

    // system 0
    AtomicSystem* system0=0;   // [1]
    Vec3d* REQ0 =0;      // [system0.n] REQparams for system0
    //Vec3i* jnds =0; // [system0.n] typToREQ  for system0
    //Vec3d* fs0  =0; // [system0.n] typToREQ  for system0

    // Temporary array for accumulation of derivs
    int nmax=0;
    Vec3d * fs=0;
    //std::vector<Vec3d> fs;

//void realoc( int nR=0, int nE=0, int nQ=0 ){
void realloc( int nDOFs_ ){
    //nDOFs=nR+nE+nQ;
    nDOFs=nDOFs_;
    _realloc(  DOFs , nDOFs );
    _realloc( fDOFs , nDOFs );
    //Rs=DOFs;Es=DOFs+nR;Qs=DOFs+nR+nE;
    //fRs=fDOFs;fEs=fDOFs+nR;fQs=fDOFs+nR+nE;
}

void tryRealocTemp(){
    int n=nmax;
    for(int i=0; i<nbatch; i++){ int ni = batch[i].n; if(ni>n){ n=ni;} }
    if(n>nmax){ _realloc( fs, n );  nmax=n; };
}

int init_types( int ntype_, Vec3i* typeMask, Vec3d* tREQs=0 ){
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
        }
    }
    realloc(nDOFs);
    //poses = new Mat3d[nbatch];
    //if( bRigid){
    //    poses = new Mat3d[nbatch];
    //}
    if(tREQs){
        ntype    = ntype_; 
        typeREQs = tREQs;
        DOFsFromTypes(); 
    }
    return nDOFs;
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
        if(tt.x>=0)DOFs[tt.y] = REQ.y;
        if(tt.x>=0)DOFs[tt.z] = REQ.z;
    }
}

void setSystem( int isys, int na, int* types, Vec3d* ps, bool bCopy=false ){
    AtomicSystem* S;
    if     (isys==-3){ S=systemTest;     }
    else if(isys==-2){ S=systemTest0;    }
    else if(isys==-1){ S=system0;        }
    else             { S=&batch[isys];   }
    if(bCopy){
        if(S==0   ){ S = new AtomicSystem(na); }
        if(S->n!=na){ printf("ERORRO FitREQ.setSystem() na(%i)!=S.n(%i) \n", na, S->n ); exit(0); }
        for(int i=0; i<na; i++){ S->types[i]=types[i]; S->ps[i]=ps[i]; }
    }else{
        if(S==0){ S = new AtomicSystem(); }
        S->n     =na;
        S->types =types;
        S->ps    =ps;
    }
}

void setRigidSamples( int n, double* Es_, Mat3d* poses_, bool bCopy=false, bool bAlloc=false ){
    nbatch = n; 
    if(bAlloc){
        if(Es_   )Es    = new double[nbatch];
        if(poses_)poses = new Mat3d [nbatch];
    }
    if(bCopy){
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
    for(int i=0; i<n; i++){
        int t = types[i];
        //const Vec3i& inds = typToREQ[t];
        const Vec3d& pi  = ps[i]; 
        const Vec3d  REQi = typeREQs[t];
        Vec3d fsi   = Vec3dZero;
        for(int j=0; j<nj; j++){
            Vec3d d           = jpos[j] - pi;
            const Vec3d& REQj = REQ0[j];
            double R  = REQi.x+REQj.x;
            double E0 = REQi.y*REQj.y;
            double Q  = REQi.x*REQj.x;

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

            //if(bStore){
                double dELJ_dR =  E0*12*ir*( u6*u4 - u4 );
                //dQ_dQi  = Qj
                //dE0_dEi = Ej
                //dR_Ri   = 1;
                // minimize U = Sum_i w_i * DE_i^2 = Sum_i w_i ( E_i,train - E_i,test )^2
                // dU_dxi = 2*DE_i* ( dE_i/dx_i )
                // fDOFs[inds.x] +=    dE*dELJ_dR       ; // dEtot/dRi
                // fDOFs[inds.y] +=    dE*dELJ_dE*REQj.y; // dEtot/dEi
                // fDOFs[inds.z] +=    dE*dEel_dQ*REQj.z; // dEtot/dQi
                // fDOFs[jnds.x] += wj*dE*dELJ_dR       ; // dEtot/dRj
                // fDOFs[jnds.y] += wj*dE*dELJ_dE*REQi.y; // dEtot/dEj
                // fDOFs[jnds.z] += wj*dE*dEel_dQ*REQi.z; // dEtot/dQj

                fsi.x += dELJ_dR;         // dEtot/dRi
                fsi.y += dELJ_dE*REQj.y;  // dEtot/dEi
                fsi.z += dEel_dQ*REQj.z;  // dEtot/dQi
                
                /*
                // --- Update fs0 ??
                fs0[j].add( Vec3d{
                    dE*dELJ_dR,        // dEtot/dRi
                    dE*dELJ_dE*REQi.y, // dEtot/dEi
                    dE*dEel_dQ*REQi.z  // dEtot/dQi
                });
                */
            //}
            //ei =  addAtomicForceLJQ( ps[j]-pi, fij, REQij );
        }
        fs[i].add(fsi);
    }
    return Etot;
}

void acumDerivs( int n, int* types, double dE){
    for(int i=0; i<n; i++){
        int t = types[i];
        const Vec3i& inds = typToREQ[t];
        const Vec3d  f = fs[i];
        fDOFs[inds.x]+=f.x;
        fDOFs[inds.y]+=f.y;
        fDOFs[inds.z]+=f.z;
    }
}

double evalDerivs(){
    tryRealocTemp();
    double Error = 0;
    for(int i=0; i<nbatch; i++){
        const AtomicSystem& C = batch[i];
        double Eref=Es[i];
        double wi = 1; 
        if(weights) wi = weights[i];
        // switch(ikind){ case ikind_LJQ:
        double E  = evalExampleDerivs_LJQ(C.n, C.types, C.ps );  // Forward pass
        double dE = (E - Eref)*wi;
        Error += dE;
        acumDerivs(C.n, C.types, dE );
        //double dE = C.E - E;
        //double E_ = evalExampleDerivs_LJQ(C.n, C.types, C.ps, dE*wi );  // Backward pass
    }
    return Error;
}

double evalDerivsRigid(){
    tryRealocTemp();
    const AtomicSystem& C0 = *systemTest0;
    const AtomicSystem& C  = *systemTest;
    double Error = 0; 
    for(int i=0; i<nbatch; i++){
        double Eref=Es[i];
        double wi = 1; 
        if(weights) wi = weights[i];
        rigid_transform( poses[i].a, 0, poses[i].c, poses[i].b, C.n, C0.ps, C.ps );
        double E  = evalExampleDerivs_LJQ( C.n, C.types, C.ps );  // Forward pass
        double dE = (E - Eref)*wi;
        Error += dE;
        acumDerivs(C.n, C.types, dE );
    }
    return Error ;
}

void move_GD( double dt ){
    for(int i=0; i<nDOFs; i++){
        DOFs[i] += fDOFs[i]*dt;
    }
}

}; // class FitREQ


#endif
