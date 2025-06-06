#ifndef EFF_h
#define EFF_h
/// @file EFF.h @brief Implements electron force-field solver based on Floating Gaussian Orbitals
/// @defgroup Electron_Forcefield  Electron_Forcefield
/// @addtogroup Electron_Forcefield
/// @{


#include "datatypes.h"
#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "Forces.h"


#include "InteractionsGauss.h"

#include "eFF_LAMMPS_funcs.h"


/*

eFF : Electron Force Field
---------------------------

[1] http://aip.scitation.org/doi/10.1063/1.3272671
The dynamics of highly excited electronic systems: Applications of the electron force field
Julius T. Su, William A. Goddard

[2] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
Non-adiabatic dynamics modeling framework for materials in extreme conditions
Hai Xiao, Andrés Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard,

NOTES:
-------

1) It seems that decrease of kinetic energy by sharing electron between two atoms is dominant contribution which makes formation of bonds fabourable
    * H2 Molecule perhaps cannot be stable without this contribution ( i.e. with fixed radius of electron blobs )

Params
            a             b             c           d           e           s-core
Al     0.486000       1.049000      0.207000                              1.660000
Si     0.320852       2.283269      0.814857                              1.691398
C     22.721015       0.728733      1.103199     17.695345    6.693621    0.621427
O     25.080199       0.331574      1.276183     12.910142    3.189333    0.167813

*/

/*
Erf approximation:
# Gaussian:    F = (x2-1)**2 / sqrtPi
# Erf          E = x*(1 + x2 * ( -0.66666666666 + 0.2*x2 ) ) * (2/(16.0/15.0))
*/

// ============= Constnts

//double wae = 1.0;
//double bAE = -5.0;
//double aAE = 20.0;

//#define QE -2.0
#define QE -1.0

constexpr static const double Hartree_to_eV = 27.211386245988;

// ToDo : Later properly
// ToDo : Perhaps we should use differenw size(width) for Pauli and Coulomb similarly to CLCFGO
constexpr static const Vec3d default_eAbWs[] = {
// amp[eV]  exp[1/A]  size[A]
{ 0.0,  0.0, 0.0},  // Q = 0 //
{ 0.0,  0.0, 0.01}, // Q = 1 // H
{ 2.0, -3.0, 0.1},  // Q = 2 // Be
{ 2.0, -3.0, 0.1},  // Q = 3 // B
{ 2.0, -3.0, 0.1},  // Q = 4 // C
{ 2.0, -3.0, 0.1},  // Q = 5 // N
{ 2.0, -3.0, 0.1},  // Q = 6 // O
{ 2.0, -3.0, 0.1},  // Q = 6 // F
};

constexpr static const Vec3d default_aAbWs[] = {
{ 0.0,  0.0, 0.1},  // Q = 0 //
{ 0.0,  0.0, 0.01}, // Q = 1 // H
{ 1.0, -5.0, 0.1},  // Q = 2 // Be
{ 1.0, -5.0, 0.1},  // Q = 3 // B
{ 1.0, -5.0, 0.1},  // Q = 4 // C
{ 1.0, -5.0, 0.1},  // Q = 5 // N
{ 1.0, -5.0, 0.1},  // Q = 6 // O
{ 1.0, -5.0, 0.1},  // Q = 6 // F
};
// SEE
// [2] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
// Non-adiabatic dynamics modeling framework for materials in extreme conditions
// Hai Xiao, Andrés Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard,

constexpr static const  double default_EPCs_sp[] = {
// s-core       a                b               c             d             e
 0.621427,   22.721015,       0.728733,      1.103199,     17.695345,    6.693621,  // C
 0.000000,   0.0000000,       0.000000,      0.000000,     0.0000000,    0.000000,  // N  // ToDo: Maybe we can try interpolate C and O ?
 0.167813,   25.080199,       0.331574,      1.276183,     12.910142,    3.189333,  // O
};
constexpr static const  double default_EPCs_ss[] = {
 1.660000,    0.486000,       1.049000,      0.207000,     -1       ,   -1       ,  // Al
 1.691398,    0.320852,       2.283269,      0.814857,     -1       ,  - 1       ,  // Si
};

// Reformulate    d=b  ,   e=c
/*
constexpr static const default_EPCs[] = {
// s-core        a               b              c             d               e
 0.621427,   22.721015,      17.695345,      6.693621,     ,    ,  // C
 0.167813,   25.080199,      12.910142,      3.189333,     ,    ,  // O
 1.660000,    0.486000,       1.049000,      0.207000,     -1       ,   -1       ,  // Al
 1.691398,    0.320852,       2.283269,      0.814857,     -1       ,  - 1       ,  // Si
};
*/

struct EFFAtomType{
    int ne;       // number of valence electrons
    double Rcore; // core electron radius
    double Acore; // core electron repulsion (depend on number of core electrons)
    // ToDo : maybe distinct core radius for pi-electrons (?)
    //double RcoreP; // core electron radius
    //double AcoreP; // core electron repulsion (depend on number of core electrons)
};


EFFAtomType EFFparams[11]{
    {0,0.,0.}, // void ... (or electron?)
    {1,0.,0.}, // H
    {2,0.,0.}, // He // only core ?

    {1,0.1,2.0}, // Li
    {2,0.1,2.0}, // Be
    {3,0.1,2.0}, // B
    {4,0.1,2.0}, // C
    {5,0.1,2.0}, // N
    {6,0.1,2.0}, // O
    {7,0.1,2.0}, // F
    {8,0.1,2.0}  // Ne // only core ?
};





inline double frEFF_r2( double r2, double w1, double w2 ){
    return ( 1/(1 + w1*r2) - 1/(1 + w2*r2) )/( 1/w1 + 1/w2 );
}

inline void combineAbW( const Vec3d& abwi, const Vec3d& abwj, Vec3d& abw ){
    abw.x = abwi.x*abwj.x;
    abw.y = (abwi.y+abwj.y)*0.5;
    abw.z = abwi.z+abwj.z;
};

inline double addPairEF_expQ( const Vec3d& d, Vec3d& f, double w2, double qq, double bExp, double aExp ){
    double r2     = d.norm2();
    double invrw2 = 1./( r2 + w2 );
    double invrw  = sqrt(invrw2);
    double E      =  const_El_eVA*qq*invrw;
    double fr     = -E*invrw2;
    if( bExp<0 ){
        double r      = sqrt( r2+R2SAFE );
        double Eexp  = aExp*exp( bExp*r );
        fr          += bExp*Eexp/r;
        E           += Eexp;
    }
    f.set_mul( d, fr );
    return E;
}

inline double interp_gx4(double r2, double y1, double y2 ){
    double c = (1-r2);
    c*=c;
    return c*y1 + (1-c)*y2;
}

/// EFF solver
class EFF{ public:

    int nAtomParams = 10;
    constexpr static const Quat4d default_AtomParams[] = {
    //  Q   sQ   sP   cP
    { 0.,  1.0, 1.0, 0.0 }, // 0
    { 1.,  0.0, 0.0, 0.0 }, // 1 H
    { 0.,  1.0, 1.0, 1.0 }, // 2 He
    { 1.,  0.0, 0.1, 1.0 }, // 3 Li
    { 2.,  0.0, 0.1, 1.0 }, // 4 Be
    { 3.,  0.0, 0.1, 1.0 }, // 5 B
    { 4.,  0.0, 0.1, 1.0 }, // 6 C
    { 5.,  0.0, 0.1, 1.0 }, // 7 N
    { 6.,  0.0, 0.1, 1.0 }, // 8 O
    { 7.,  0.0, 0.1, 1.0 }, // 9 F
    };


    constexpr static const double8 default_AtomParams2[] = {
    // Z_nuc, R_eff, Zcore_eff,   PA,        PB,        PC,        PD,        PE
{ 0.0,   1.0,      0.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 0: Dummy
{ 1.0,   0.0,      0.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 1: H (Bare nucleus)
{ 2.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 2: He (Simple core: 1 pair (2e), radius 0.3. sQ=0.3, sP=0.3, cP=1.0)
{ 3.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 3: Li (Simple core: Z=3, sQ=0.5, sP=0.5, cP=1.0 for 1s2)
{ 4.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 4: Be (Simple core: Z=4, sQ=0.4, sP=0.4, cP=1.0 for 1s2)
{ 5.0,   0.1,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 5: B  (Simple core: Z=5, sQ=0.35,sP=0.35,cP=1.0 for 1s2)
{ 6.0,   0.621427, 2.0,     22.721015, 0.728733,  1.103199,  17.695345, 6.693621  }, // 6: C (ECP: Z_nuc=6, R_core=0.621, Z_core=2. p-type)
{ 7.0,   0.0,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      }, // 7: N (ECP: Z_nuc=7, R_core=0.0,   Z_core=2. p-type)
{ 8.0,   0.167813, 2.0,     25.080199, 0.331574,  1.276183,  12.910142, 3.189333  }, // 8: O (ECP: Z_nuc=8, R_core=0.167, Z_core=2. p-type)
{ 9.0,   0.3,      2.0,      0.0,       0.0,       0.0,       0.0,       0.0      }  // 9: F (Simple core: Z=9, sQ=0.3, sP=0.3, cP=1.0 for 1s2)
// Add Al, Si etc. as needed
};

    const Quat4d*  atom_params   = default_AtomParams;
    const double8* atom_params2 = default_AtomParams2;

//                                              H   He  Li    Be      B      C     N     O      F
constexpr static const double aMasses[9] = {  1.0, 4.0, 7.0, 9.0,  11.0,  12.0,  14.0, 16.0,  19.0 };

    bool bDealoc = false;

    bool bNegativeSizes=false;

    int iPauliModel = 1;
    bool bCoreCoul  = true;

    double KPauliOverlap = 50.0; // ToDo : This is just "bulgarian constant" for now
    double KPauliKin     = 50.0; // ToDo : Not sure if we should use this - perhaps this model of pauli energy should work "ab inition"

    constexpr static const double default_esize = 0.5;
    constexpr static const double min_esize     = 0.1;
    constexpr static const Vec3d KRSrho = { 1.125, 0.9, -0.2 }; ///< eFF universal parameters
    //constexpr static const Vec3d KRSrho = { 1.125, 0.9, -0.3 }; ///< eFF universal parameters // If it is sufficiently strong electron pairs are formed
    //constexpr static const Vec3d KRSrho = { 1.125, 0.9, -1.0 }; ///< eFF universal parameters // If it is sufficiently strong electron pairs are formed
    //Vec3d KRSrho = { 1.125, 0.9, 0.2 };

    bool bEvalKinetic   = true;
    bool bEvalEE        = true;
    bool bEvalCoulomb   = true;
    bool bEvalPauli     = true;
    bool bEvalAE        = true;
    bool bEvalAECoulomb = true;
    bool bEvalAEPauli   = true;
    bool bEvalAA        = true;
    bool bEvalCoreCorect = true;

    int ne=0,na=0,nDOFs=0; ///< number of electrons, atoms, degrees of freedom
    //int*   atype  =0;

    Vec3d  * apos   =0; ///< atomic positions
    Vec3d  * aforce =0; ///< atomic forces
    Vec3d  * avel   =0; ///< atomic velocities

    Quat4d * aPars  = 0;   /// electron params { x=Q,y=sQ,z=sP,w=cP }
    Quat4d * aPars2 = 0;   /// electron params { x=Q,y=sQ,z=sP,w=cP }

    bool  bUseECPs = false;
    int   iECPmodel = 0;

    //double * espin  =0;
    int    * espin  =0; ///< electron spins
    double * echarge =0; ///< electron charge (QE for individual, 2*QE for pair)
    Vec3d  * epos   =0; ///< electron positions
    Vec3d  * eforce =0; ///< electron forces
    Vec3d  * evel   =0; ///< electron velocities
    double * esize  =0; ///< electron size
    double * fsize  =0; ///< electron force on size
    double * vsize  =0; ///< electron velocity on size
    double * eE = 0;

    double* pDOFs =0;  ///< buffer of degrees of freedom
    double* fDOFs =0;  ///< buffer of forces on degrees of freedom
    double* vDOFs =0;  ///< buffer of velocities on degrees of freedom
    double* invMasses=0; ///< buffer of inverse masses [nDOFs] (both electrons and nuclei)

    int nfix=0;
    Quat4d* fixed_poss=0; // [nfix], {x,y,z,w} position of fixed particles
    Vec2i*  fixed_inds=0; // [nfix]  {ia/-ie, bitmask{x|y|z|w}}, ia/-ie is index of atom/ or negative index of electron, bitmax indicate which component is fixed

    double Etot=0,Ek=0, Eee=0,EeePaul=0,EeeExch=0,  Eae=0,EaePaul=0,  Eaa=0; ///< different kinds of energy

void realloc(int na_, int ne_, bool bVel=false){
    bDealoc = true;
    na=na_; ne=ne_;
    nDOFs=na*3+ne*3 + ne;
    _realloc( pDOFs, nDOFs);
    _realloc( fDOFs, nDOFs);
    _realloc( aPars, na);
    _realloc( aPars2, na);
    _realloc( espin, ne);
    _realloc( eE, ne );
    _realloc( echarge, ne);  for(int i=0; i<ne; i++){ echarge[i]=1; }

    apos   = (Vec3d*)pDOFs;
    aforce = (Vec3d*)fDOFs;

    epos   = (Vec3d*)(pDOFs + na*3);
    eforce = (Vec3d*)(fDOFs + na*3);
    esize  =          pDOFs + na*3 + ne*3;
    fsize  =          fDOFs + na*3 + ne*3;

    if(bVel){
        _realloc( vDOFs, nDOFs      ); // velocities
        _realloc( invMasses, nDOFs  ); // velocities
        avel  = (Vec3d*)vDOFs;
        evel  = (Vec3d*)(vDOFs + na*3);
        vsize =          vDOFs + na*3 + ne*3;
        for(int i=0; i<nDOFs; i++){ vDOFs[i]=0;     }
        for(int i=0; i<nDOFs; i++){ invMasses[i]=1; }
    }

}

void realloc_fixed(int nfix_){
    nfix=nfix_;
    _realloc0( fixed_poss, nfix, Quat4dNAN );
    _realloc0( fixed_inds, nfix, Vec2iZero );
}

void apply_hard_fix(){
    //printf( "fFF:apply_fixed() nfix %i \n", nfix );
    for(int i=0; i<nfix; i++){
        const int     ia   = fixed_inds[i].x;
        const int     mask = fixed_inds[i].y;
        const Quat4d& p0   = fixed_poss[i];
        // print bit mask as binary
        //printf( "fFF:apply_fixed() i %3i ia %3i bit %08b pos0 %16.8f %16.8f %16.8f %16.8f \n", i, ia, mask, p0.x, p0.y, p0.z, p0.w );
        if(ia<0){
            int ie = -ia-1;
            // if(mask&1){ epos[ie].x = p0.x; eforce[ie].x = 0; evel[ie].x = 0; }
            // if(mask&2){ epos[ie].y = p0.y; eforce[ie].y = 0; evel[ie].y = 0; }
            // if(mask&4){ epos[ie].z = p0.z; eforce[ie].z = 0; evel[ie].z = 0; }
            // if(mask&8){ esize[ie]  = p0.w; fsize[ie]    = 0; vsize[ie]  = 0; }

            if(mask&1){ epos[ie].x = p0.x; }
            if(mask&2){ epos[ie].y = p0.y; }
            if(mask&4){ epos[ie].z = p0.z; }
            if(mask&8){ esize[ie]  = p0.w; }
        }else{
            //if(mask&1){ apos[ia].x = p0.x; aforce[ia].x = 0; avel[ia].x = 0; }
            //if(mask&2){ apos[ia].y = p0.y; aforce[ia].y = 0; avel[ia].y = 0; }
            //if(mask&4){ apos[ia].z = p0.z; aforce[ia].z = 0; avel[ia].z = 0; }

            if(mask&1){ apos[ia].x = p0.x;  }
            if(mask&2){ apos[ia].y = p0.y; }
            if(mask&4){ apos[ia].z = p0.z;  }
        }
    }
}

void clear_fixed_force(){
    //printf( "fFF:apply_fixed() nfix %i \n", nfix );
    for(int i=0; i<nfix; i++){
        const int     ia   = fixed_inds[i].x;
        const int     mask = fixed_inds[i].y;
        const Quat4d& p0   = fixed_poss[i];
        // print bit mask as binary
        //printf( "fFF:apply_fixed() i %3i ia %3i bit %08b pos0 %16.8f %16.8f %16.8f %16.8f \n", i, ia, mask, p0.x, p0.y, p0.z, p0.w );
        if(ia<0){
            int ie = -ia-1;
            if(mask&1){ eforce[ie].x = 0; }
            if(mask&2){ eforce[ie].y = 0; }
            if(mask&4){ eforce[ie].z = 0; }
            if(mask&8){ fsize[ie]    = 0; }
        }else{
            if(mask&1){ aforce[ia].x = 0; }
            if(mask&2){ aforce[ia].y = 0; }
            if(mask&4){ aforce[ia].z = 0; }
        }
    }
}

void clear_fixed_dynamics(){
    //printf( "fFF:apply_fixed() nfix %i avel=%p evel=%p vsize=%p \n", nfix, avel, evel, vsize );
    for(int i=0; i<nfix; i++){
        const int     ia   = fixed_inds[i].x;
        const int     mask = fixed_inds[i].y;
        const Quat4d& p0   = fixed_poss[i];
        // print bit mask as binary
        //printf( "fFF:apply_fixed() i %3i ia %3i bit %08b pos0 %16.8f %16.8f %16.8f %16.8f \n", i, ia, mask, p0.x, p0.y, p0.z, p0.w );
        if(ia<0){
            int ie = -ia-1;
            if(mask&1){ eforce[ie].x = 0; evel[ie].x = 0; }
            if(mask&2){ eforce[ie].y = 0; evel[ie].y = 0; }
            if(mask&4){ eforce[ie].z = 0; evel[ie].z = 0; }
            if(mask&8){ fsize[ie]    = 0; vsize[ie]  = 0; }
        }else{
            if(mask&1){ aforce[ia].x = 0; avel[ia].x = 0; }
            if(mask&2){ aforce[ia].y = 0; avel[ia].y = 0; }
            if(mask&4){ aforce[ia].z = 0; avel[ia].z = 0; }
        }
    }
}

void makeMasses(double*& invMasses, double m_const=-1){
    //double atomic_mass   = 1.6605402e-27;
    //double electron_mass = 9.1093837e-31; 
    double au_Me           = 1822.88973073;
    double eV_MeAfs        = 17.5882001106;   // [ Me*A^2/fs^2]  
    if(invMasses==0) invMasses = new double[nDOFs];
    if(m_const>0){
        for(int i=0; i<nDOFs; i++){ invMasses[i] = m_const; }  
    }else{
        int na3 =na*3;
        int ne3 =ne*3; 
        double* buff=invMasses; for(int i=0; i<na;  i++){ double m = eV_MeAfs/( aMasses[ (int)(aPars[i].x-0.5) ] * au_Me ); int i3=i*3; buff[i3]=m;buff[i3+1]=m;buff[i3+2]=m;  } // assign atomic masses   [Me] i.e. in units of electron mass
        buff+=na3;              for(int i=0; i<ne3; i++){ buff[i]  = eV_MeAfs;         }                                      // assign electron masses [Me] i.e. in units of electron mass
        buff+=ne3;              for(int i=0; i<ne;  i++){ buff[i]  = 0.01*eV_MeAfs/( 0.5 )*echarge[i]; }                                      // assign electron size massses ????  ToDo:  What should be this mass ?????
    }
    /*
    // Force units:  [eV/A]
    // Mass  units:  [eV/]
    time [fs]:   
        Me = 9.1093837e-31;    // [kg]
        eV = 1.602176634e-19;  // [J]
        A  = 1e-10;            // [m]
        fs = 1e-15;            // [s]

        a        = F/m = (E/l)/m 
        [m/s^2]  = [J/m]/[kg]
        [A/fs^2] = [eV/A]/[C*Me]
        [C*Me]/9.1093837e-31 = [eV*fs^2/A^2]/((1.602176634e-19)*(1.e-15/1.e-10)^2)
        [C*Me]/9.1093837e-31 = [eV*fs^2/A^2]/1.602176634e-29
        // C  = 17.5882001106 = 1/0.05685630102
        // eV = 17.5882001106 [Me*A^2/fs^2]   

        // Force of [eV/A] acting per one 1[fs] accelerates one electron (or other particle with mass 1[Me]) to speed 17.5882001106 [A/fs]  

    */
}

void dealloc(){
    delete [] pDOFs;
    delete [] fDOFs;
    //delete [] aQ   ;
    //delete [] aAbWs;
    //delete [] eAbWs;
    delete [] aPars;
    delete [] aPars2;
    delete [] espin;
    delete [] echarge;
    delete [] eE;
}
~EFF(){ if(bDealoc)dealloc(); }

void fixElectron(int ie, double* vs=0){
    eforce[ie].set(0.);
    fsize [ie]=0;
    if(vs){  
        int io=(na+ie)*3;
        vs[io  ]=0; // epos.x velocity
        vs[io+1]=0; // epos.y velocity
        vs[io+2]=0; // epos.z velocity
        vs[(na+ne)*3+ie]=0;  // size velocity
    };
}

void set_electron(int ie, const Vec3d& pos, double size, int spin, double charge=1.0){
    epos   [ie] = pos;
    esize  [ie] = size;
    espin  [ie] = spin;
    echarge[ie] = charge;
}

/// evaluate kinetic energy of each electron
double evalKinetic(){
    Ek=0;
    bNegativeSizes=false;
    for(int i=0; i<ne; i++){
        //double fs=0;
        //double dEk = addKineticGauss( esize[i]*M_SQRT1_2, fs );
        //fsize[i]+=fs*M_SQRT1_2;
        //if( i_DEBUG>0 ) printf( "evalKinetic[%i] s %g -> f %g Ek %g \n", i, esize[i], fsize[i], Ek );
        double s = esize[i];
        if(s<min_esize){ s=min_esize; esize[i]=s; bNegativeSizes=true; };
        double dEk = addKineticGauss_eFF( esize[i], fsize[i], echarge[i] );
        eE[i] =dEk;
        Ek   +=dEk;
        //if(verbosity>2){ printf("%s e%i Ke %5.20f \n",prefix, i,dEk); }
    }
    //if( bNegativeSizes & (verbosity>0) ){ printf( "negative electron sizes => perhaps decrease relaxation time step? \n" ); }
    return Ek;
}

/// evaluate Electron-Electron forces
double evalEE(){
    Eee    =0;
    EeePaul=0;
    //double w2ee = wee*wee;
    //const double qq = QE*QE;
    for(int i=0; i<ne; i++){
        const Vec3d    pi  = epos[i];
        //Vec3d&   fi  = eforce[i];
        const int8_t spini = espin[i];
        //const double   si  = esize[i];
        const double   si  = esize[i];
        double&       fsi  = fsize[i];
        const double  qi   = echarge[i];
        for(int j=0; j<i; j++){
            Vec3d  f  = Vec3dZero;
            const Vec3d  dR = epos [j] - pi;
            //const double sj = esize[j];
            const double sj = esize[j];
            double&     fsj = fsize[j];
            const double qj = echarge[j];
            const double qij = qi*qj;

            const int spinj  = espin[j];
            const int spinij = spini*spinj;
            const double q2  = (spini==0)&&(spinj==0)? 2 : 1.0;

            double dEee=0,dEpaul=0;
            if(bEvalCoulomb){
                dEee = addCoulombGauss( dR, si, sj, f, fsi, fsj, qij );
                //dEee = addCoulombGauss( dR, si*M_SQRT2, sj*M_SQRT2, f, fsi, fsj, qq );
                //dEee = addCoulombGauss( dR, si*2, sj*2, f, fsi, fsj, qq );
            }
            if(bEvalPauli){
                //printf( "Eee[%i,%i]= %g(%g) r %g s(%g,%g) \n", i, j, dEee, Eee, dR.norm(), si,sj );
                if( iPauliModel == 1 ){ // Pauli repulsion form this eFF paper http://aip.scitation.org/doi/10.1063/1.3272671
                    //if( spini==espin[j] ){
                        //printf( "EeePaul_1[%i,%i]  ", i, j );
                        //printf( "evalEE() r %g pi (%g,%g,%g) pj (%g,%g,%g) \n", dR.norm(), epos[j].x,epos[j].y,epos[j].z, pi.x,pi.y,pi.z  );
                        //dEpaul = addPauliGauss  ( dR, si, sj, f, fsi, fsj, spini!=espin[j], KRSrho );
                        //dEpaul = addPauliGauss_New  ( dR, si, sj, f, fsi, fsj, spini!=espin[j], KRSrho );
                        dEpaul = addPauliGauss_New  ( dR, si, sj, f, fsi, fsj, spinij, KRSrho, q2 );
                        //printf( "EeePaul[%i,%i]= %g \n", i, j, dEpaul );
                    //}
                }else if( iPauliModel == 2 ){ // iPauliModel==0 Pauli repulasion from Valence-Bond theory
                    if(spini==espin[j]){ // Pauli repulsion only for electrons with same spin
                        //printf( "EeePaul[%i,%i] >> ", i, j );
                        //double dEpaul = addPauliGaussVB( dR, si*M_SQRT2, sj*M_SQRT2, f, fsi, fsj ); EeePaul+=dEpaul;
                        dEpaul = addPauliGaussVB( dR, si, sj, f, fsi, fsj );
                        //printf( "EeePaul[%i,%i]= %g \n", i, j, dEpaul );
                        //printf( "<< dEpaul %g \n", dEpaul );
                    }
                }else{   // iPauliModel==0 Pauli repulasion as overlap between same spin orbitals 
                    if( spini==espin[j] ){
                        //printf( "EeePaul[%i,%i] ", i, j );
                        //i_DEBUG=1;
                        dEpaul = addDensOverlapGauss_S( dR, si*M_SQRT2, sj*M_SQRT2, KPauliOverlap, f, fsi, fsj );
                        //double dEpaul = addPauliGauss  ( dR, si, sj, f, fsi, fsj, false, KRSrho ); EeePaul+=dEpaul;
                        //i_DEBUG=0;
                        //printf( "EeePaul[%i,%i]= %g \n", i, j, dEpaul );
                    }
                }
            }
            //if(verbosity>2){ printf("%s e%i-e%i Coul %5.20f \n",prefix,i,j,dEee); printf("%s e%i-e%i Paul %5.20f \n",prefix,i,j,dEpaul); }
            Eee    += dEee;
            EeePaul+= dEpaul;
            double dE = 0.5*( dEee + dEpaul );
            eE[i]+=dE;
            eE[j]+=dE;
            //if( spini==espin[j] ) EeePaul += addDensOverlapGauss_S( dR,si,sj, 1, f, fsi, fsj );
            //Eee += addPairEF_expQ( epos[j]-pi, f, w2ee, +1, 0, 0 );
            //if( i_DEBUG>0 ) printf( "evalEE[%i,%i] dR(%g,%g,%g) s(%g,%g) q %g  ->   f(%g,%g,%g) fs(%g,%g) \n", i,j, dR.x,dR.y,dR.z, si,sj, qq,   f.x,f.y,f.z, fsi,fsj );
            eforce[j].sub(f);
            eforce[i].add(f);

            //DEBUG_fa_aa[j].sub(f);
            //DEBUG_fa_aa[i].add(f);
            //if( (i==DEBUG_i)&&(j==DEBUG_j) ){
            //    glColor3f(1.0,0.0,0.0);
            //    Draw3D::drawVecInPos( f*-1., epos[j] );
            //    Draw3D::drawVecInPos( f   , pi      );
            //}
        }
    }
    //printf( "Eee %g EeePaul %g \n", Eee, EeePaul );
    //if( i_DEBUG>0 )  for(int j=0; j<ne; j++){  printf( "evalEE: esize[%i] %g f %g \n", j, esize[j], fsize[j] ); }
    return Eee+EeePaul;
}

/// evaluate Atom-Electron forces
double evalAE(){
    Eae    =0;
    EaePaul=0;
    double Eee_=0;
    //double see2 = see*see;
    //double saa2 = saa*saa;
    //double invSae = 1/( see*see + saa*saa );
    //double w2ae = wae*wae;
    for(int i=0; i<na; i++){
        const Vec3d  pi   = apos[i];
        //const double qqi  = aQ[i]*QE;
        //const Vec3d  abwi = eAbWs[i];
        const Quat4d aPar = aPars[i]; // { x=Q,y=sQ,z=sP,w=cP }
        const double qq  = aPar.x*QE;
        for(int j=0; j<ne; j++){
            Vec3d f=Vec3dZero;
            const Vec3d   dR  = epos [j] - pi;
            const double  sj  = esize[j];
            const double  qj  = echarge[j];
            double& fsj = fsize[j];
            double  fs_junk=0;
            //Eae += addPairEF_expQ( epos[j]-pi, f, abwi.z, qi*QE, abwi.y, abwi.x );
            //Eae  += addCoulombGauss( dR,sj,      f, fsj,      qqi );     // correct
            double dEae=0,dEaePaul=0,dEee=0;
            if(bEvalAECoulomb){
                dEae  = addCoulombGauss( dR, aPar.y, sj, f, fs_junk, fsj, qq );
            }
            if( bEvalAEPauli && (aPar.w>1e-8) ){
                //if(qqi<-1.00001) EaePaul += addDensOverlapGauss_S( dR,sj, abwi.z, abwi.a, f, fsj, fs_junk );     // correct
                //double dEaePaul = addPauliGauss      ( dR, sj, abwi.z, f, fsj, fs_junk, false, KRSrho );     // correct
                //double dEaePaul = addDensOverlapGauss_S( dR, sj, aPar.z, aPar.w, f, fsj, fs_junk );     // correct
                dEaePaul            = addPauliGauss_New( dR, sj, aPar.z, f, fsj, fs_junk, 0, KRSrho, qj      ); // spin=0 means both -1 and +1  
                if(bCoreCoul){ dEee = addCoulombGauss  ( dR, sj, aPar.z, f, fsj, fs_junk,            qj*2.0  ); }
                //printf( "EaePaul[%i,%i] E %g r %g s %g abw(%g,%g) \n", i, j, dEaePaul, dR.norm(), sj, abwi.z, abwi.a );
            }
            //if( i_DEBUG>0 ) printf( "evalAE[%i,%i] dR(%g,%g,%g) s %g q %g  ->   f(%g,%g,%g) fs %g \n", i,j, dR.x,dR.y,dR.z, sj, qqi,   f.x,f.y,f.z, fsj );
            //printf( "evalAE[%i,%i] E %g r %g s(%g,%g) \n", i,j, Eae, dR.norm(), aPar.y, sj );
            if(verbosity>2){ 
                //printf("%s a%i-e%i Coul %5.20f \n",prefix,i,j,dEae); 
                //printf("%s a%i-e%i Paul %5.20f \n",prefix,i,j,dEaePaul); 
            }
            //printf( "evalAE[%i,%i] f(%g,%g,%g) r %g s %g E %g EPaul %g Eee %g \n", i,j, f.x,f.y,f.z, dR.norm(), aPar.y, sj, dEae,dEaePaul,dEee );
            Eae    +=dEae;
            EaePaul+=dEaePaul;
            Eee_   +=dEee;
            eE[j]  +=dEae+dEaePaul+dEee;
            eforce[j].sub(f);
            aforce[i].add(f);

            //DEBUG_fe_ae[j].sub(f);
            //DEBUG_fa_ae[i].add(f);
            //if( (i==DEBUG_i)&&(j==DEBUG_j) ){
            //    glColor3f(1.0,1.0,1.0); Draw3D::drawLine    ( pi, epos[j]    );
            //    glColor3f(0.0,1.0,0.0); Draw3D::drawVecInPos( f*-1., epos[j] );
            //    glColor3f(1.0,0.0,1.0); Draw3D::drawVecInPos( f    , pi      );
            //}
        }
    }
    Eee+=Eee_;
    //if( i_DEBUG>0 )  for(int j=0; j<ne; j++){  printf( "evalAE: esize[%i] %g f %g \n", j, esize[j], fsize[j] ); }
    return Eae+EaePaul+Eee_;
}


double evalECP( Vec3d dR, Quat4d aPar, Quat4d BCDE, double sj, Vec3d& f, double& fs ){
    double dEae = 0;
    if(bEvalAECoulomb){
        double fs_junk;
        dEae += addCoulombGauss(dR, aPar.z, sj, f, fs_junk, fs, (aPar.x - aPar.z) );   // ( Z - Zeff )
    }
    // Pauli: pseudo-atom (i) - electron (j)
    if(bEvalAEPauli){
        const double rc = dR.norm();
        double fe_ = 0.0;
        double dE_ = 0.0;
        double fr_ = 0; 
        double fs_ = 0;
        if (BCDE.w>0){ LAMMPS_NS::PauliCoreElec (rc, sj, &dE_, &fr_, &fe_, aPar.w, BCDE.x, BCDE.y); } 
        else         { LAMMPS_NS::PauliCorePElec(rc, sj, &dE_, &fr_, &fe_, aPar.w, BCDE.x, BCDE.y, BCDE.z, BCDE.w);  }
        dEae += dE_ * Hartree_to_eV;
        fs   -= fe_ * Hartree_to_eV;
        f.add_mul(dR, (fr_ * Hartree_to_eV) / rc);
    }
    return dEae;
}

double evalAE_ECP(){
    Eae    =0;
    EaePaul=0;
    double Eee_=0;
    for(int i=0; i<na; i++){
        const Vec3d  pi   = apos[i];
        const Quat4d aPar = aPars[i]; // {x=Z,y=size,z=Zeff,w=A}
        const double qq  = aPar.x*QE;
        const Quat4d& BCDE = aPars2[i];
        const bool is_eCP  = ( aPar.w  > 0.0 );

        for(int j=0; j<ne; j++){
            Vec3d f=Vec3dZero;
            const Vec3d   dR  = epos [j] - pi;
            const double  sj  = esize[j];
            const double  qj  = echarge[j];
            double& fsj = fsize[j];
            double  fs_junk=0;
            double dEae=0,dEaePaul=0,dEee=0;

            if(bEvalAECoulomb){
                dEae  = addCoulombGauss(dR, aPar.z, sj, f, fs_junk, fsj, (aPar.x - aPar.z) );   // ( Z - Zeff )
            }
            if(bEvalAEPauli){
                if(is_eCP){
                    const double rc = dR.norm();
                    double fe_ = 0.0;
                    double dE_ = 0.0;
                    double fr_ = 0; 
                    double fs_ = 0;
                    switch( iECPmodel ){
                        case 0:
                            if (BCDE.w>0){ LAMMPS_NS::PauliCoreElec (rc, sj, &dE_, &fr_, &fe_, aPar.w, BCDE.x, BCDE.y); } 
                            else         { LAMMPS_NS::PauliCorePElec(rc, sj, &dE_, &fr_, &fe_, aPar.w, BCDE.x, BCDE.y, BCDE.z, BCDE.w);  }
                            break;
                        case 1:
                            if (BCDE.w>0){ LAMMPS_NS::PauliCoreElec_optGemini (rc, sj, dE_, fr_, fe_, aPar.w, BCDE.x, BCDE.y); }
                            else         { LAMMPS_NS::PauliCorePElec_optGemini(rc, sj, dE_, fr_, fe_, aPar.w, BCDE.x, BCDE.y, BCDE.z, BCDE.w); }
                            break;
                        case 2:
                            if (BCDE.w>0){ LAMMPS_NS::PauliCoreElec_optDeepSeek (rc, sj, dE_, fr_, fe_, aPar.w, BCDE.x, BCDE.y); }
                            else         { LAMMPS_NS::PauliCorePElec_optDeepSeek(rc, sj, dE_, fr_, fe_, aPar.w, BCDE.x, BCDE.y, BCDE.z, BCDE.w); }
                            break;
                        default:
                            break;
                    }
                    dEaePaul    += dE_*qj * Hartree_to_eV;
                    fsj         -= fe_*qj * Hartree_to_eV;
                    f.add_mul(dR, (fr_*qj * Hartree_to_eV) / rc);
                }else{
                    dEaePaul   = addPauliGauss_New( dR, sj, aPar.z, f, fsj, fs_junk, 0, KRSrho, qj ); // spin=0 means both -1 and +1  
                    //dEaePaul   = addPauliGauss_New( dR, sj, aPar.z, f, fsj, fs_junk, 0, KRSrho, aPar.z*0.5 ); // spin=0 means both -1 and +1 
                }
            }
            //printf( "evalAE[%i,%i] f(%g,%g,%g) r %g s %g E %g EPaul %g Eee %g \n", i,j, f.x,f.y,f.z, dR.norm(), aPar.y, sj, dEae,dEaePaul,dEee );
            Eae    +=dEae;
            EaePaul+=dEaePaul;
            Eee_   +=dEee;
            eE[j]  +=dEae+dEaePaul+dEee;
            eforce[j].sub(f);
            aforce[i].add(f);
        }
    }
    Eee+=Eee_;
    //if( i_DEBUG>0 )  for(int j=0; j<ne; j++){  printf( "evalAE: esize[%i] %g f %g \n", j, esize[j], fsize[j] ); }
    return Eae+EaePaul+Eee_;
}


/// evaluate Atom-Atom forces
double evalAA(){
    //if( i_DEBUG>0 ) printf( "evalAA \n" );
    Eaa=0;
    //double invSaa = 1/(saa*saa);
    //double w2aa = waa*waa;
    for(int i=0; i<na; i++){
        const Vec3d  pi   = apos[i];
        //const double qi   = aQ[i];
        //const Vec3d  abwi = aAbWs[i];
        const Quat4d aPari = aPars[i];
        for(int j=0; j<i; j++){
            Vec3d f = Vec3dZero; // HERE WAS THE ERROR (missing initialization !!!! )
            //Vec3d  abw;
            //combineAbW( abwi, aAbWs[j], abw );
            const Quat4d& aParj = aPars[j];
            const Vec3d   dR    = apos[j] - pi;
            //if( (i_DEBUG>0) && (1==qi==aQ[j]) ){ printf( " abw(H-H): %i,%i A %g B %g w %g \n", i,j, abw.x, abw.y, abw.z ); }
            //Eaa += addPairEF_expQ( apos[j]-pi, f, abw.z, qi*aQ[j], abw.y, abw.x );
            //Eaa += addPairEF_expQ( apos[j]-pi, f, abw.z, qi*aQ[j], abw.y, abw.x );
            double qq;
            if(bCoreCoul){
                qq = (aPari.x-aPari.w*2)*(aParj.x-aParj.w*2);
                //printf( "evalAA()[%i,%i] qq %g qi(%g,%g) qi(%g,%g)\n", i,j, qq, aPari.x, aPari.w, aParj.x, aParj.w );
            }else{
                qq = aPari.x*aParj.x;
            }
            double dEaa =  addAtomicForceQ( dR, f, qq );
            //if(verbosity>2){ printf("%s a%i-a%i Coul %5.20f \n",prefix,i,j,dEaa); }
            Eaa+=dEaa;
            //printf( " Eaa[%i,%i]  Q %g(%g,%g) \n", i, j, aPari.x*aParj.x, aPari.x, aParj.x  );
            //   ToDo : Pauli Repulsion of core electrons ?????
            aforce[j].sub(f);
            aforce[i].add(f);

            //DEBUG_fa_aa[j].sub(f);
            //DEBUG_fa_aa[i].add(f);
            //if( (i==DEBUG_i)&&(j==DEBUG_j) ){
            //    glColor3f(1.0,0.0,0.0);
            //    //Draw3D::drawVecInPos( dR*-1., apos[j] );
            //    //Draw3D::drawVecInPos( dR    , pi      );
            //    Draw3D::drawVecInPos( f*-1., apos[j] );
            //   Draw3D::drawVecInPos( f    , pi      );
            //}
        }
    }
    return Eaa;
}

double evalCoreCorrection(){
    double Ek_  =0;
    double Eee_ =0;
    double Eae_ =0;
    for(int i=0; i<na; i++){
        double fsi=0,fsj=0;
        Vec3d  f=Vec3dZero;
        const Quat4d aPar = aPars[i]; // { x=Q,y=sQ,z=sP,w=cP }
        if(aPar.w<1e-8) continue;
        double s = aPar.z;
        double w = aPar.w;
        Ek_       += addKineticGauss_eFF( s, fsi, w*2.0 );                                        // kineatic energy of core electrons
      //EaePaul_  += addPauliGauss_New( Vec3dZero, s, s, f, fsi, fsj, -1, KRSrho, aPar.w );             // Pauli repulsion of core electrons
        Eee_      += addCoulombGauss  ( Vec3dZero, s, s, f, fsi, fsj,            aPar.w );             // core electrons with each other
        Eae_      += addCoulombGauss  ( Vec3dZero, aPar.y, s, f, fsi, fsj,       aPar.w*-2.0*aPar.x ); // nuclie with core electrons
        // WARRNING :  This does not apply any forces ????
    }
    //printf( "DEBUG evalCoreCorrection() Ek0 %g Eee0 %g Eae0 %g \n", Ek , Eee , Eae  );
    //printf( "DEBUG evalCoreCorrection() Ek_ %g Eee_ %g Eae_ %g \n", Ek_, Eee_, Eae_ );
    Ek+=Ek_; Eee+=Eee_; Eae+=Eae_;
    return Ek_+Eee_+Eae_;
}


/// evaluate full Electron Forcefild
double eval(){
    //printf( "DEBUG eval() verbosity %i na %i ne %i \n", verbosity, na, ne );
    //clearForce();
    clearForce_noAlias();
    Etot = 0;
    if(bEvalKinetic    ) Etot+= evalKinetic();
    if(bEvalEE         ) Etot+= evalEE();
    if(bEvalAE){
        if(bUseECPs){ Etot+= evalAE_ECP(); }
        else        { Etot+= evalAE();     }
    }
    if(bEvalAA         ) Etot+= evalAA();
    if(bEvalCoreCorect ) Etot+=evalCoreCorrection();
    //printf( "eval() Etot %g epos[0](%g,%g,%g) \n", Etot, epos[0].x, epos[0].y, epos[0].z );
    return Etot;
}

void clearForce(){
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0; }
}

void clearForce_noAlias(){
    for(int i=0;i<na;i++){
        aforce[i]=Vec3dZero;
    }
    for(int i=0;i<ne;i++){
        eforce[i]=Vec3dZero;
        fsize[i]=0;
        eE[i] = 0;
    }
}

double move_GD(double dt){
    double F2sum = 0;
    for(int i=0;i<nDOFs;i++){
        double f  = fDOFs[i];
        pDOFs[i] += f * dt;
        F2sum+=f*f;
    }
    return F2sum;
    //printf( "move_GD sum %g ", sum );
}

void move_GD_noAlias(double dt){
    //Vec3d fe=Vec3dZero,fa=Vec3dZero;
    //double fs=0;
    for(int i=0;i<na;i++){
        apos[i].add_mul( aforce[i], dt );
        //fa.add( aforce[i] );
    }
    for(int i=0;i<ne;i++){
        epos [i].add_mul( eforce[i], dt );
        //fe.add( aforce[i] );
        esize[i] += fsize[i] * dt;
        //fs += fsize[i];
    }
    //printf( "fs %g fe(%g,%g,%g) fa(%g,%g,%g)\n", fs, fe.x,fe.y,fe.z, fa.x,fa.y,fa.z );
}

/*
void autoAbWs( const Vec3d * AAs, const Vec3d * AEs ){
    for(int i=0;i<na;i++){
        int ityp =(int)round(aQ[i]);
        aAbWs[i] = AAs[ityp];
        eAbWs[i] = AEs[ityp];
        //printf( "atom[%i] typ %i aAbW(%g,%g,%g) eAbW(%g,%g,%g) \n", i, ityp, aAbWs[i].x, aAbWs[i].y, aAbWs[i].z,  eAbWs[i].x, eAbWs[i].y, eAbWs[i].z );
    }
}
*/

double electronPotAtPoint( const Vec3d& pi, double si, double Q, int spini=0, bool bEvalCoulomb=true )const{
    double EeeCoul=0;
    double EeePaul=0;
    Vec3d fp; double fsi;
    bool bEvalPauli = (spini!=0);
    // ToDo: perhaps we should make put this loop into a separate function
    for(int j=0; j<ne; j++){
        Vec3d  f  = Vec3dZero;
        const Vec3d  dR = epos [j] - pi;
        //const double sj = esize[j];
        const double sj = esize[j];
        double&     fsj = fsize[j];
        double dEee=0,dEpaul=0;
        if(bEvalCoulomb){
            dEee = addCoulombGauss( dR, si, sj, f, fsi, fsj, Q );
        }
        if(bEvalPauli){
            if( iPauliModel == 1 ){ // Pauli repulsion form this eFF paper http://aip.scitation.org/doi/10.1063/1.3272671
                dEpaul = addPauliGauss_New  ( dR, si, sj, f, fsi, fsj, spini*espin[j], KRSrho );
                //printf( "EeePaul[j=%i] dEpaul=%g spins(%i*%i)\n", j, dEpaul , spini, espin[j] ); // it seems it gives almost ZERO for some reason => investigate
            }else if( iPauliModel == 2 ){ 
                if(spini==espin[j]){ 
                    dEpaul = addPauliGaussVB( dR, si, sj, f, fsi, fsj );
                }
            }else{  
                if( spini==espin[j] ){
                    dEpaul = addDensOverlapGauss_S( dR, si*M_SQRT2, sj*M_SQRT2, KPauliOverlap, f, fsi, fsj );
                }
            }
        }
        EeeCoul+= dEee;
        EeePaul+= dEpaul;
        //printf( "EeePaul[j=%i] EeePaul=%g EeeCoul=%g spins(%i*%i)\n", j, EeePaul, EeeCoul, spini, espin[j] );
    }
    //printf( "electronPotAtPoint() END EeePaul=%g EeeCoul=%g \n", EeePaul, EeeCoul );
    return EeeCoul + EeePaul;
}

double atomsPotAtPoint( const Vec3d& pos, double s, double Q )const{
    double Eae    =0;
    double EaePaul=0;
    Vec3d fp; double fs;
    // ToDo: perhaps we should make put this loop into a separate function
    for(int i=0; i<na; i++){
        const Vec3d  dR   = pos-apos[i];
        const Quat4d aPar = aPars[i]; // { x=Q,y=sQ,z=sP,w=cP }
        const double qq   = aPar.x*Q;
        Vec3d  f=Vec3dZero;
        double fsi=0,fsj=0;
        //addCoulombGauss      ( const Vec3d& dR, double si, double sj,             Vec3d& f, double& fsi, double& fsj, double qq ){
        //addDensOverlapGauss_S( const Vec3d& dR, double si, double sj, double amp, Vec3d& f, double& fsi, double& fsj            ){
        Eae  += addCoulombGauss              ( dR, s, aPar.y,         f, fsi, fsj, qq );
        if( aPar.w>1e-8 ){
            EaePaul+= addDensOverlapGauss_S( dR, s, aPar.z, aPar.w, f, fsi, fsj         );
        }
    }
    return Eae + EaePaul;
}

double* evalPotAtPoints( int n, Vec3d* ps, double* out=0, double s=0.0, double Q=1.0, int spin=0, bool bAtom=true, bool bElectron=false )const{
    if(out==0){ out = new double[n]; };
    for(int i=0; i<n; i++){
        double E = 0;
        if(bAtom    ){ E += atomsPotAtPoint   ( ps[i], s, Q       ); }
        if(bElectron){ E += electronPotAtPoint( ps[i], s, Q, spin ); }
        out[i] = E;
    }
    return out;
}

void printEnergies(){
    printf( "Etot %16.8f | Ek %16.8f Eee,p(%16.8f,%16.8f) Eae,p(%16.8f,%16.8f) Eaa %g \n", Etot, Ek, Eee,EeePaul, Eae,EaePaul, Eaa );
}

void printAtoms(){
    //printf( "Etot %g Ek %g Eel %g(ee %g, ea %g aa %g)  EPaul %g(ee %g, ae %g) \n", Etot, Ek, Eel, Eee,Eae,Eaa,   EPaul, EeePaul, EaePaul );
    for(int i=0; i<na; i++){
        //printf( "a[%i] p(%g,%g,%g) q %g eAbW(%g,%g,%g) aAbW(%g,%g,%g) \n", i, apos[i].x, apos[i].y, apos[i].z, aQ[i], eAbWs[i].z,eAbWs[i].z,eAbWs[i].z, aAbWs[i].z,aAbWs[i].z,aAbWs[i].z );
        printf( "a[%3i ] p(%16.8f ,%16.8f ,%16.8f ) Par(Q,sQ,sP,P)(%g,%g,%g,%g)  \n", i, apos[i].x, apos[i].y, apos[i].z, aPars[i].x,aPars[i].y,aPars[i].z,aPars[i].w );
        //printf( "a[%i] xyzs(%g,%g,%g) fxyzs(%g,%g,%g) \n", i, apos[i].x, apos[i].y, apos[i].z, aforce[i].x, aforce[i].y, aforce[i].z );
    }
}

void printElectrons(){
    for(int i=0; i<ne; i++){
        printf( "e[%3i ] p(%16.8f ,%16.8f ,%16.8f ) sz %16.8f s %i \n", i, epos[i].x, epos[i].y, epos[i].z, esize[i], espin[i] );
        //printf( "e[%i] xyzs(%g,%g,%g,%g) fxyzs(%g,%g,%g,%g) \n", i, ff.epos[i].x, ff.epos[i].y, ff.epos[i].z, ff.esize[i], ff.eforce[i].x, ff.eforce[i].y, ff.eforce[i].z, ff.fsize[i] );
    }
}

void printSwitches(){
    printf( "iPauliModel %i KPauliOverlap %g bCoreCoul %i \n", iPauliModel, KPauliOverlap, bCoreCoul );
    printf( "bEvalKinetic %i bEvalEE %i bEvalCoulomb %i bEvalPauli %i bEvalAE %i bEvalAECoulomb %i bEvalAEPauli %i bEvalAA %i bEvalCoreCorect %i  \n", bEvalKinetic, bEvalEE, bEvalCoulomb, bEvalPauli, bEvalAE, bEvalAECoulomb, bEvalAEPauli, bEvalAA, bEvalCoreCorect  );
    printf( "KRSrho %g %g %g \n", KRSrho.x, KRSrho.y, KRSrho.z );
}

void info(){
    printSwitches();
    printAtoms();
    printElectrons();
}

int Eterms2str(char* str){
    // Ek=0,Eee=0,EeePaul=0,EeeExch=0,Eae=0,EaePaul=0,Eaa=0, Etot=0;
    double Eorbs = 0;
    for(int i=0; i<ne; i++){ Eorbs+=eE[i]; }
    return sprintf( str, "Etot %3.3f|%3.3f Ek %3.3f Eee,P(%3.3f,%3.3f) Eae,P(%3.3f,%3.3f) Eaa %g )\n", Etot, Eorbs+Eaa, Ek, Eee, EeePaul, Eae, EaePaul, Eaa );
}

int orb2str(char* str, int ie){
    if(vDOFs){
        return sprintf( str, "e[%i|%i] E %7.3f p(%5.2f,%5.2f,%5.2f|%5.2f) f(%5.2f,%5.2f,%5.2f|%5.2f) v(%5.2f,%5.2f,%5.2f|%5.2f)\n", ie, espin[ie], eE[ie], epos[ie].x,epos[ie].y,epos[ie].z,esize[ie],    eforce[ie].x,eforce[ie].y,eforce[ie].z,fsize[ie],  evel[ie].x,evel[ie].y,evel[ie].z,vsize[ie]  );
    }else{
        return sprintf( str, "e[%i|%i] E %7.3f p(%5.2f,%5.2f,%5.2f|%5.2f) f(%5.2f,%5.2f,%5.2f|%5.2f)\n", ie, espin[ie], eE[ie],   epos[ie].x,epos[ie].y,epos[ie].z,esize[ie],    eforce[ie].x,eforce[ie].y,eforce[ie].z,fsize[ie]  );
    }
    //return sprintf( str, "e[%i] E %7.3f s %5.2f  p(%5.2f,%5.2f,%5.2f) \n", ie, eE[ie], esize[ie], epos[ie].x,epos[ie].y,epos[ie].z );
}

char* orbs2str(char* str0){
    char* str=str0;
    for(int i=0; i<ne; i++){
        //str+=sprintf(str,"\norb_%i E %g \n", i, oEs[i] );
        str+=orb2str(str, i);
        
    }
    return str;
}

void to_xyz( FILE* pFile, const char* comment=0 ){
    fprintf( pFile, " %i \n", na+ne );
    if(comment!=0){ 
        fprintf( pFile, "na,ne %i %i Etot(%g)=T(%g)+ee(%g)+ea(%g)+aa(%g) | %s \n", na,ne, Etot, Ek, Eee, Eae, Eaa, comment ); 
    }else{ 
        fprintf( pFile, "na,ne %i %i Etot(%g)=T(%g)+ee(%g)+ea(%g)+aa(%g) \n", na,ne, Etot, Ek, Eee, Eae, Eaa ); 
    }
    for (int i=0; i<na; i++){
        int iZ = (int)(aPars[i].x+0.5);
        //if(iZ>1)iZ+=2;
        fprintf( pFile, "%3i %10.6f %10.6f %10.6f \n", iZ, apos[i].x, apos[i].y, apos[i].z );
    }
    for (int i=0; i<ne; i++){
        int e = espin[i];
        if(e==1 ){e=92;}
        if(e==-1){e=109;} // see Jmol colors https://jmol.sourceforge.net/jscolors/
        fprintf( pFile, "%3i %10.6f %10.6f %10.6f 0 %10.6f \n", e, epos[i].x, epos[i].y, epos[i].z,  esize[i]*4.0 );
    }
}

void save_xyz( const char* filename, const char* mode="w", const char* comment=0 ){
    //printf( "EFF::save_xyz(%s)\n", filename );
    FILE * pFile; 
    pFile = fopen (filename,mode);
    if(pFile==0){ printf("ERROR file >>%s<< not found \n", filename ); return; }
    to_xyz( pFile, comment );
    fclose(pFile);
}

bool loadFromFile_xyz( const char* filename ){
    //printf(" filename: >>%s<< \n", filename );
    FILE * pFile;
    pFile = fopen (filename,"r");
    if (pFile==0){ printf("ERROR file >>%s<< not found \n", filename ); return true; }
    int ntot;
    fscanf (pFile, " %i \n", &ntot );
    fscanf (pFile, " %i %i\n", &na, &ne );
    if((ne+na)!=ntot){ printf("ERROR ne(%i)+na(%i) != ntot(%i) >>%s<< \n", ne, na, ntot, filename ); fclose(pFile); return true; }
    //printf("na %i ne %i ntot %i \n", na, ne, ntot );
    realloc( na, ne );
    char buf[1024];
    //int ntot=na+ne;

    double Qasum = 0.0;

    int ia=0,ie=0;
    //int ia=0,ie=nDOFs-1;
    //int ia=0;
    //int ie=nDOFs-1;
    for (int i=0; i<ntot; i++){
        double x,y,z,s;
        int e;
        fgets( buf, 256, pFile); //printf( ">%s<\n", buf );
        //printf( "buf[%i]  >>%s<< \n", ie, buf );
        int nw = sscanf (buf, " %i %lf %lf %lf %lf", &e, &x, &y, &z, &s );
        if( e<0){
            epos[ie]=Vec3d{x,y,z};
            if     (e==-1){ espin[ie]= 1; }
            else if(e==-2){ espin[ie]=-1; }
            if( nw>4 ){ esize[ie]=s; }else{ esize[ie]=default_esize; }
            //esize[ie] = 1.0;
            //esize[ie] = 0.5;
            //esize[ie] = 0.25;
            //printf( " e[%i] pos(%g,%g,%g) spin %i size %g | nw %i \n", ie, epos[ie].x, epos[ie].y, epos[ie].z, espin[ie], esize[ie], nw );
            ie++;
        }else{
            apos[ia]=Vec3d{x,y,z};
            //aQ  [ia]=e;  // change this later
            //aAbws[ia] = default_aAbWs[e];
            //eAbws[ia] = default_eAbWs[e];
            aPars [ia] = default_AtomParams[e]; ;
            Qasum += e;
            ia++;
            //printf( " a[%i] ", ia );
        };
        //printf( " %i %f %f %f  \n", e, x,y,z );
    }
    clearForce();
    //printf( "Qtot = %g (%g - 2*%i) \n",  Qasum - ne, Qasum, ne );
    fclose (pFile);
    return 0;
}


int readSingleFGO(FILE* pFile, bool bVel=true, double fUnits=1., bool bRealloc=true){
    //printf("readSingleFGO() bVel %i fUnits %g bRealloc %i pFile %p\n",bVel, fUnits, bRealloc, pFile);
    int ntot;
    const int nbuff = 1024;
    char buff[nbuff]; char* line;
    //fscanf (pFile, " %i \n", &ntot );
    int natom_=0, nOrb_=0, perOrb_=0; bool bClosedShell=0;
    int iline=0;
    line=fgets(buff,nbuff,pFile); iline++;
    if(line==0){ printf("ERROR in eFF::readSingleFGO : No line found \n"); return -1;}
    sscanf (line, "%i %i %i\n", &natom_, &nOrb_, &perOrb_, &bClosedShell );
    //printf("readSingleFGO() natom_ %i nOrb_ %i perOrb_ %i bClosedShell %i\n", natom_, nOrb_, perOrb_, bClosedShell );
    if(perOrb_!=1){ printf("ERROR in eFF::readSingleFGO : perOrb must =1 ( found %i instead) !!! \n", perOrb_ ); return -1;};
    if(bClosedShell) nOrb_*=2;
    if(bRealloc)realloc( natom_, nOrb_, bVel );
    //printf("readSingleFGO() natom_ %i nOrb_ %i na %i ne %i\n", natom_, nOrb_, na, ne );
    double Qasum = 0.0;
    for(int i=0; i<na; i++){
        double x,y,z;
        double vx,vy,vz;
        double Q,sQ,sP,cP;
        //int bfix;
        fgets( buff, nbuff, pFile); iline++;
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &x, &y, &z, &Q, &sQ, &sP, &cP, &vx, &vy, &vz );
        Q=-Q;
        apos  [i]=Vec3d{x*fUnits,y*fUnits,z*fUnits};
        aPars[i].set(Q,sQ*fUnits,sP*fUnits,cP);
        if( (bVel)&&(nw>=10) ){
            avel[i] = Vec3d{vx*fUnits,vy*fUnits,vz*fUnits};
        }
        Qasum += Q;
    }
    int nBasRead = ne;
    if( bClosedShell ) nBasRead/=2;
    for(int i=0; i<nBasRead; i++){
        double x,y,z;
        double vx,vy,vz,vs,vc;
        double s,c;
        int spin=0;
        fgets( buff, nbuff, pFile); iline++;
        int nw = sscanf (buff, "%lf %lf %lf %lf %lf %i %lf %lf %lf %lf", &x, &y, &z, &s, &c, &spin, &vx, &vy, &vz, &vs );
        epos [i]=Vec3d{x*fUnits,y*fUnits,z*fUnits};
        esize[i]=s*fUnits;
        if(bVel){ 
            if(nw>=9 )evel [i] = Vec3d{vx*fUnits,vy*fUnits,vz*fUnits};
            if(nw>=10)vsize[i] = vs*fUnits;
        }
        if( !bClosedShell ){ if(nw>5)espin[i]=spin; }else{ espin[i]=1; };
    }
    if( bClosedShell ){
        for(int i=0; i<nBasRead; i++){
            int j = i+nBasRead;
            epos [j]=epos[i];
            esize[j]=esize[i];
            espin[j]=-1;
        }
    }
    return iline;
}

bool loadFromFile_fgo( char const* filename, bool bVel=true, double fUnits=1. ){
    FILE * pFile = fopen(filename,"r");
    if( pFile == NULL ){
        printf("ERROR in eFF::loadFromFile_fgo(%s) : No such file !!! \n", filename );
        return -1;
    }
    int iline = readSingleFGO(pFile, bVel, fUnits);
    fclose(pFile);
    return iline < 2;
}

void writeTo_fgo( char const* filename, bool bVel=false, const char* fmode="w", int iconf=-1 ){
    //printf(" writeTo_fgo(%s, %s) \n", filename, fmode);
    FILE * pFile = fopen (filename, fmode );
    if(iconf>=0){
        fprintf( pFile, "%i %i %i #iconf %i\n", na, ne, 1, iconf );
    }else{
        fprintf( pFile, "%i %i %i\n", na, ne, 1 );
    }
    for(int i=0; i<na; i++){               
                fprintf(pFile, "%f %f %f   %f %f %f %f ", apos[i].x, apos[i].y, apos[i].z,   -aPars[i].x, aPars[i].y, aPars[i].z, aPars[i].w );
        if(bVel)fprintf(pFile, "%f %f %f", avel[i].x, avel[i].y, avel[i].z  );
        fprintf(pFile,"\n");
    }
    for(int i=0; i<ne; i++){
                fprintf(pFile, "%f %f %f %f %f  %i ", epos[i].x, epos[i].y, epos[i].z, esize[i], 1.0,  espin[i] );
        if(bVel)fprintf(pFile, "%f %f %f",  evel[i].x, evel[i].y, evel[i].z, vsize[i], 0.0 );
        fprintf(pFile,"\n");
    }
    fclose (pFile);
}

inline double analyse_distance(int ia,int ib)const{
    Vec3d pi=(ia<0)?epos[-ia-1]:apos[ia];
    Vec3d pj=(ib<0)?epos[-ib-1]:apos[ib];
    Vec3d d=pj-pi;
    return sqrt(d.norm2());
}

inline double analyse_angle(int i,int j,int k)const{
    Vec3d pi=(i<0)?epos[-i-1]:apos[i];
    Vec3d pj=(j<0)?epos[-j-1]:apos[j];
    Vec3d pk=(k<0)?epos[-k-1]:apos[k];
    Vec3d u=pi-pj;
    Vec3d v=pk-pj;
    double nu=sqrt(u.norm2());
    double nv=sqrt(v.norm2());
    double ca=(nu>0&&nv>0)?(u.dot(v)/(nu*nv)):1.;
    if(ca>1)ca=1;
    if(ca<-1)ca=-1;
    return acos(ca);
}

inline void analyse_distances(const Vec2i* pairs,int nPairs,double* out){
    for(int i=0;i<nPairs;i++)out[i]=analyse_distance(pairs[i].x,pairs[i].y);
}

inline void analyse_angles(const Vec3i* triples,int nTriples,double* out){
    for(int i=0;i<nTriples;i++)out[i]=analyse_angle(triples[i].x,triples[i].y,triples[i].z);
}


    // ============ Output Handling Methods ==========
    void copyAtomPositions(Vec3d* apos_out, int iconf=0) const {
        if(apos_out){
            Vec3d* apos = apos_out + iconf*na;
            for(int j=0; j<na; j++){ apos[j] = this->apos[j]; }
        }
    }

    void copyElectronPositions(Quat4d* epos_out, int iconf=0) const {
        if(epos_out){
            Quat4d* epos = epos_out + iconf*ne;
            for(int j=0; j<ne; j++){ 
                epos[j].f = this->epos[j]; 
                epos[j].w = this->esize[j]; 
            }
        }
    }

    void copyEnergies(double* outEs, int iconf=0, const int nEperConf=5) const {
        if(outEs){
            double* outEi = outEs + iconf*nEperConf;
            outEi[0] = Etot;
            outEi[1] = Ek;
            outEi[2] = Eee;
            outEi[3] = Eae;
            outEi[4] = Eaa;
        }
    }
};

/// @}
#endif
