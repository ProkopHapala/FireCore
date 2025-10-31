/*

// ============= For automatic generation of interfaces

int nnode
int natom
int nvec = nnode+natom

//------- Dynamical

_RW float4*  apos      [ nSys* nvec ]
_RW float4*  aforce    [ nSys* nvec ]
_RW float4*  avel      [ nSys* nvec ]
_RW float4*  fneigh    [ nSys* nnode*2*4 ]

//------- parameters

_R  int4*    neighs    [ nSys* natom ]
_R  int4*    bkNeighs  [ nSys* natom ]
_R  int4*    neighCell [ nSys* natom ]
_R  float4*  REQKs     [ nSys* natom ]
_R  float4*  apars,    [ nSys* nnode ]
_R  float4*  bLs,      [ nSys* nnode ]
_R  float4*  bKs,      [ nSys* nnode ]
_R  float4*  Ksp,      [ nSys* nnode ]
_R  float4*  Kpp,      [ nSys* nnode ]
_R  cl_Mat3* lvecs,    [ nSys ]
_R  cl_Mat3* ilvecs    [ nSys ]

_R float4*  apos_surf   [ natom_surf ]
_R float4*  aforce_surf [ natom_surf ]

_R float4*  ps,         [ ndipol ]
_R float4*  dipols,     [ ndipol ]

_RW image3d_t  FE_Paul[ ng.x, ng.y, ng.z ]
_RW image3d_t  FE_Lond[ ng.x, ng.y, ng.z ]
_RW image3d_t  FE_Coul[ ng.x, ng.y, ng.z ]

const float4   MDpars
const int4     nGrid
const cl_Mat3  dGrid
const float4   grid_p0

*/


/*
    relax_multi.cl -  OpenCL kernel source code for multi-system relaxation

    The file contains various OpenCL kernels icluding:


    This file is part of FireCode project.

*/



// ======================================================================
// ======================================================================
//                           FUNCTIONS
// ======================================================================
// ======================================================================

#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable

typedef struct __attribute__ ((packed)){
    float4 a;
    float4 b;
    float4 c;
} cl_Mat3;

#define  float4Zero  (float4){0.f,0.f,0.f,0.f}
#define  float3Zero  (float3){0.f,0.f,0.f}
#define  float2Zero  (float3){0.f,0.f,0.f}

#define R2SAFE          1e-4f
#define COULOMB_CONST   14.3996448915f       // [ eV*Ang/e^2 ]
#define const_kB        8.617333262145e-5f   // [ eV/K ]



inline float2 udiv_cmplx( float2 a, float2 b ){ return (float2){  a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y }; }     // divison of unitary complex numbers (i.e. rotation backwards)
//inline void     udiv_cmplx(               const VEC& b ){                            T x_ =    x*b.x +   y*b.y;         y =    y*b.x -   x*b.y;       x=x_;  }

inline float3 rotMat ( float3 v, float3 a, float3 b, float3 c ){ return (float3)(dot(v,a),dot(v,b),dot(v,c)); }  // rotate vector v by matrix (a,b,c)
inline float3 rotMatT( float3 v, float3 a, float3 b, float3 c ){ return a*v.x + b*v.y + c*v.z; }                 // rotate vector v by matrix (a,b,c) transposed

// evaluate angular force and energy using cos(angle) formulation,    - faster, but not good for angles > 90 deg
inline float evalAngCos( const float4 hr1, const float4 hr2, float K, float c0, __private float3* f1, __private float3* f2 ){
    float  c = dot(hr1.xyz,hr2.xyz);
    float3 hf1,hf2;
    hf1 = hr2.xyz - hr1.xyz*c;
    hf2 = hr1.xyz - hr2.xyz*c;
    float c_   = c-c0;
    float E    = K*c_*c_;
    float fang = -K*c_*2;
    hf1 *= fang*hr1.w;
    hf2 *= fang*hr2.w;
    *f1=hf1;
    *f2=hf2;
    return E;
}

// evaluate angular force and energy using cos(angle/2) formulation - a bit slower, but not good for angles > 90 deg
inline float evalAngleCosHalf( const float4 hr1, const float4 hr2, const float2 cs0, float k, __private float3* f1, __private float3* f2 ){
    // This is much better angular function than evalAngleCos() with just a little higher computational cost ( 2x sqrt )
    // the main advantage is that it is quasi-harmonic beyond angles > 90 deg
    float3 h  = hr1.xyz + hr2.xyz;  // h = a+b
    float  c2 = dot(h,h)*0.25f;     // cos(a/2) = |ha+hb|  (after normalization)
    float  s2 = 1.f-c2 + 1e-7;      // sin(a/2) = sqrt(1-cos(a/2)^2) ;  s^2 must be positive (otherwise we get NaNs)
    float2 cso = (float2){ sqrt(c2), sqrt(s2) }; // cso = cos(a/2) + i*sin(a/2)
    float2 cs = udiv_cmplx( cs0, cso );          // rotate back by equilibrium angle
    float  E         =  k*( 1 - cs.x );          // E = k*( 1 - cos(a/2) )  ; Do we need Energy? Just for debugging ?
    float  fr        = -k*(     cs.y );          // fr = k*( sin(a/2) )     ; force magnitude
    c2 *= -2.f;
    fr /=  4.f*cso.x*cso.y;   //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    float  fr1    = fr*hr1.w; // magnitude of force on atom a
    float  fr2    = fr*hr2.w; // magnitude of force on atom b
    *f1 =  h*fr1  + hr1.xyz*(fr1*c2);  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| ); force on atom a
    *f2 =  h*fr2  + hr2.xyz*(fr2*c2);  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| ); force on atom b
    return E;
}

// evaluate angular force and energy for pi-pi alignment interaction
inline float evalPiAling( const float3 h1, const float3 h2,  float K, __private float3* f1, __private float3* f2 ){  // interaction between two pi-bonds
    float  c = dot(h1,h2); // cos(a) (assumes that h1 and h2 are normalized)
    float3 hf1,hf2;        // working forces or direction vectors
    hf1 = h2 - h1*c;       // component of h2 perpendicular to h1
    hf2 = h1 - h2*c;       // component of h1 perpendicular to h2
    bool sign = c<0; if(sign) c=-c; // if angle is > 90 deg we need to flip the sign of force
    float E    = -K*c;     // energy is -K*cos(a)
    float fang =  K;       // force magnitude
    if(sign)fang=-fang;    // flip the sign of force if angle is > 90 deg
    hf1 *= fang;           // force on atom a
    hf2 *= fang;           // force on atom b
    *f1=hf1;
    *f2=hf2;
    return E;
}

// evaluate bond force and energy for harmonic bond stretching
inline float evalBond( float3 h, float dl, float k, __private float3* f ){
    float fr = dl*k;   // force magnitude
    *f = h * fr;       // force on atom a
    return fr*dl*0.5;  // energy
}

// evaluate non-covalent interaction force and energy for Lennard-Jones (Q) and Coulomb interactions of charges (Q) and hydrogen bond correction (pseudo-charges H), damping R2damp is used to avoid singularity at r=0
inline float4 getLJQH( float3 dp, float4 REQ, float R2damp ){
    // ---- Electrostatic (damped Coulomb potential)
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 +  R2damp);              // inverse distance squared and damped
    float   Ec    =  COULOMB_CONST*REQ.z*sqrt( ir2_ ); // Ec = Q1*Q2/sqrt(r^2+R2damp)
    // --- Lennard-Jones and Hydrogen bond correction
    float  ir2 = 1.f/r2;          // inverse distance squared
    float  u2  = REQ.x*REQ.x*ir2; // u2 = (R0/r)^2
    float  u6  = u2*u2*u2;        // u6 = (R0/r)^6
    float vdW  = u6*REQ.y;        // vdW = E0*(R0/r)^6
    float E    =       (u6-2.f)*vdW     + Ec  ;     // E = E0*(R0/r)^6 - E0*(R0/r)^12 + Q1*Q2/sqrt(r^2+R2damp)
    float fr   = -12.f*(u6-1.f)*vdW*ir2 - Ec*ir2_;  // fr = -12*E0*( (R0/r)^8/r + 12*E0*(R0/r)^14) - Q1*Q2/(r^2+R2damp)^1.5
    return  (float4){ dp*fr, E };
}

inline float4 getMorseQH( float3 dp,  float4 REQH, float K, float R2damp ){
    float r2    = dot(dp,dp);
    float ir2_  = 1/(r2+R2damp);
    float r     = sqrt( r2   );
    float ir_   = sqrt( ir2_ );     // ToDo: we can save some cost if we approximate r^2 = r^2 + R2damp;
    float e     = exp ( K*(r-REQH.x));
    //double e2    = e*e;
    //double fMors =  E0*  2*K*( e2 -   e ); // Morse
    //double EMors =  E0*      ( e2 - 2*e );
    float   Ae  = REQH.y*e;
    float fMors = Ae*  2*K*(e - 1); // Morse
    float EMors = Ae*      (e - 2);
    float Eel   = COULOMB_CONST*REQH.z*ir_;
    float fr    = fMors/r - Eel*ir2_ ;
    return  (float4){ dp*fr, EMors+Eel };
}

// evaluate damped Coulomb potential and force
inline float4 getCoulomb( float3 dp, float R2damp ){
    // ---- Electrostatic
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 + R2damp);
    float   E    = COULOMB_CONST*sqrt( ir2_ );
    return  (float4){ dp*-E*ir2_, E };
}

// limit force magnitude to fmax
float3 limnitForce( float3 f, float fmax ){
    float fr2 = dot(f,f);                         // force magnitude squared
    if( fr2>(fmax*fmax) ){ f*=(fmax/sqrt(fr2)); } // if force magnitude is larger than fmax we scale it down to fmax
    return f;
}

float4 getR4repulsion( float3 d, float R, float Rcut, float A ){
    // we use R4blob(r) = A * (1-r^2)^2
    // such that at distance r=R we have force f = fmax
    // f = -dR4blob/dr = 4*A*r*(1-r^2) = fmax
    // A = fmax/(4*R*(1-R^2))
    float R2    = R*R;
    float R2cut = Rcut*Rcut;
    float r2 = dot(d,d);
    if( r2>R2cut ){
        return (float4){0.0f,0.0f,0.0f,0.0f};
    }else if( r2>R2 ){
        float mr2 = R2cut-r2;
        float fr = A*mr2;
        return (float4){ d*(-4*fr), fr*mr2 };
    }else{
        float mr2 = R2cut-R2;
        float fr  = A*mr2;
        float r    = sqrt(r2);
        float fmax = 4*R*fr;
        return (float4){ d* (-fmax/r), fmax*(R-r) + fr*mr2 };
    }
}


// ======================================================================
// ======================================================================
//                           MMFF kernells
// ======================================================================
// ======================================================================


// ======================================================================
//                          getMMFFf4()
// ======================================================================

// 1.  getMMFFf4() - computes bonding interactions between atoms and nodes and its neighbors (max. 4 neighbors allowed), the resulting forces on atoms are stored "fapos" array and recoil forces on neighbors are stored in "fneigh" array
//                   kernel run over all atoms and all systems in parallel to exploit GPU parallelism
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void getMMFFf4(
    const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2  [natoms]     positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  fapos,        // 3  [natoms]     forces on    atoms (just node atoms are evaluated)
    __global float4*  fneigh,       // 4  [nnode*4*2]  recoil forces on neighbors (and pi-orbitals)
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atom  cell index
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} i.e. R0: van der Waals radii, E0: well depth and partial charge, Q: partial charge
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}, i.e. c0ss: cos(equlibrium angle/2) for sigma-sigma; Kss: stiffness of sigma-sigma angle; c0sp: is cos(equlibrium angle) for sigma-pi
    __global float4*  bLs,          // 8  [nnode]  bond length    between node and each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness between node and each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor     (only node atoms have pi-pi alignemnt interaction)
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor (only node atoms have pi-pi alignemnt interaction)
    __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
    __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
    __global float4*  pbc_shifts,   // 14 pbc shifts for each system
    const int npbc,                 // 15 number of pbc shifts
    const int bSubtractVdW          // 16 subtract vdW energy
){

    const int iG = get_global_id (0);   // intex of atom   (iG<nAtoms)
    const int iS = get_global_id (1);   // index of system (iS<nS)
    //const int nG = get_global_size(0);
    //const int nS = get_global_size(1);  // number of systems
    //const int iL = get_local_id  (0);
    //const int nL = get_local_size(0);
    const int nAtoms=nDOFs.x;  // number of atoms in the system
    const int nnode =nDOFs.y;  // number of nodes in the system
    //const int nvec  = nAtoms+nnode;

    if(iG>=nnode) return;

    const int i0a   = iS*nAtoms;         // index of first atom      in the system
    const int i0n   = iS*nnode;          // index of first node atom in the system
    const int i0v   = iS*(nAtoms+nnode); // index of first vector    in the system ( either atom or pi-orbital )

    const int iaa = iG + i0a;  // index of current atom (either node or capping atom)
    const int ian = iG + i0n;  // index of current node atom
    const int iav = iG + i0v;  // index of current vector ( either atom or pi-orbital )

    #define NNEIGH 4

    // ---- Dynamical
    float4  hs [4];              // direction vectors of bonds (h.xyz) and inverse bond lengths (h.w)
    float3  fbs[4];              // force on neighbor sigma    (fbs[i] is sigma recoil force on i-th neighbor)
    float3  fps[4];              // force on neighbor pi       (fps[i] is pi    recoil force on i-th neighbor)
    float3  fa  = float3Zero;    // force on center atom positon

    float E=0;                   // Total Energy of this atom
    // ---- Params
    const int4   ng  = neighs[iaa];    // neighboring atoms
    const float3 pa  = apos[iav].xyz;  // position of current atom
    const float4 par = apars[ian];     // (xy=s0_ss,z=ssK,w=piC0 ) forcefield parameters for current atom


    // Temp Arrays
    const int*   ings  = (int*  )&ng; // neighboring atoms, we cast it to int[] to be index it in for loop


    const float   ssC0   = par.x*par.x - par.y*par.y;                      // cos(2) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf , where ang0 is equilibrium angle
    for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }   // clear recoil forces on neighbors

    float3 f1,f2;         // working forces

    if((iG==0)&&(iS==0)){
        printf( "getMMFFf4() iG %i, iS %i, iaa %i bSubtractVdW %i\n", iG, iS, iaa, bSubtractVdW );
    }

    { // ========= BONDS - here we evaluate pairwise interactions of node atoms with its 4 neighbors

        float3  fpi = float3Zero;                // force on pi-orbital
        const int4   ngC = neighCell[iaa];       // neighboring atom cell index
        const float3 hpi = apos[iav+nAtoms].xyz; // direction of pi-orbital
        const float4 vbL = bLs[ian];             // bond lengths
        const float4 vbK = bKs[ian];             // bond stiffness
        const float4 vKs = Ksp[ian];             // stiffness of sigma-pi othogonalization
        const float4 vKp = Kpp[ian];             // stiffness of pi-pi    alignment

        const int*   ingC  = (int*  )&ngC;   // neighboring atom cell index (we cast it to int[] to be index it in for loop)
        const float* bL    = (float*)&vbL;   // bond lengths
        const float* bK    = (float*)&vbK;   // bond stiffness
        const float* Kspi  = (float*)&vKs;   // stiffness of sigma-pi othogonalization
        const float* Kppi  = (float*)&vKp;   // stiffness of pi-pi    alignment

        const int ipbc0 = iS*npbc;  // index of first PBC shift for current system

        for(int i=0; i<NNEIGH; i++){  // loop over 4 neighbors
            float4 h;                 // direction vector of bond
            const int ing  = ings[i]; // index of i-th neighbor node atom
            const int ingv = ing+i0v; // index of i-th neighbor vector
            const int inga = ing+i0a; // index of i-th neighbor atom
            if(ing<0) break;

            // --- Compute bond direction vector and inverse bond length
            h.xyz    = apos[ingv].xyz - pa;  // direction vector of bond
            { // shift bond to the proper PBC cell
                int ic  = ingC[i];                  // index of i-th neighbor cell
                h.xyz  += pbc_shifts[ipbc0+ic].xyz; // shift bond to the proper PBC cell
            }
            float  l = length(h.xyz);  // compute bond length
            h.w      = 1./l;           // store ivnerse bond length
            h.xyz   *= h.w;            // normalize bond direction vector
            hs[i]    = h;              // store bond direction vector and inverse bond length

            float epp = 0; // pi-pi    energy
            float esp = 0; // pi-sigma energy

            // --- Evaluate bond-length stretching energy and forces
            if(iG<ing){
                E+= evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;   // harmonic bond stretching, fa is force on center atom, fbs[i] is recoil force on i-th neighbor,

                // pi-pi alignment interaction
                float kpp = Kppi[i];
                if( (ing<nnode) && (kpp>1.e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                    epp += evalPiAling( hpi, apos[ingv+nAtoms].xyz, kpp,  &f1, &f2 );   fpi+=f1;  fps[i]+=f2;    //   pi-alignment(konjugation), fpi is force on pi-orbital, fps[i] is recoil force on i-th neighbor's pi-orbital
                    E+=epp;
                }
            }

            // pi-sigma othogonalization interaction
            float ksp = Kspi[i];
            if(ksp>1.e-6){
                esp += evalAngCos( (float4){hpi,1.}, h, ksp, par.w, &f1, &f2 );   fpi+=f1; fa-=f2;  fbs[i]+=f2;    //   pi-planarization (orthogonality), fpi is force on pi-orbital, fbs[i] is recoil force on i-th neighbor
                E+=epp;
            }
        }

        // --- Store Pi-forces                      we store pi-forces here because we don't use them in the angular force evaluation
        const int i4p=(iG + iS*nnode*2 )*4 + nnode*4; // index of first pi-force for current atom
        for(int i=0; i<NNEIGH; i++){
            fneigh[i4p+i] = (float4){fps[i],0}; // store recoil pi-force on i-th neighbor
        }
        fapos[iav+nAtoms]  = (float4){fpi,0};  // store pi-force on pi-orbital of current atom

    }

    { //  ============== Angles   - here we evaluate angular interactions between pair of sigma-bonds of node atoms with its 4 neighbors

        for(int i=0; i<NNEIGH; i++){ // loop over first bond
            int ing = ings[i];
            if(ing<0) break;         // if there is no i-th neighbor we break the loop
            const float4 hi = hs[i];
            const int ingv = ing+i0v;
            const int inga = ing+i0a;
            for(int j=i+1; j<NNEIGH; j++){ // loop over second bond
                int jng  = ings[j];
                if(jng<0) break;           // if there is no j-th neighbor we break the loop
                const int jngv = jng+i0v;
                const int jnga = jng+i0a;
                const float4 hj = hs[j];

                E += evalAngleCosHalf( hi, hj, par.xy, par.z, &f1, &f2 );    // evaluate angular force and energy using cos(angle/2) formulation
                fa    -= f1+f2;

                if(bSubtractVdW)
                { // Remove non-bonded interactions from atoms that are bonded to common neighbor
                    float4 REQi=REQKs[inga];   // non-bonding parameters of i-th neighbor
                    float4 REQj=REQKs[jnga];   // non-bonding parameters of j-th neighbor
                    // combine non-bonding parameters of i-th and j-th neighbors using mixing rules
                    float4 REQij;
                    REQij.x  = REQi.x  + REQj.x;
                    REQij.yz = REQi.yz * REQj.yz;

                    float3 dp = (hj.xyz/hj.w) - (hi.xyz/hi.w);   // recover vector between i-th and j-th neighbors using stored vectos and inverse bond lengths, this should be faster than dp=apos[jngv].xyz-apos[ingv].xyz; from global memory
                    float4 fij = getLJQH( dp, REQij, 1.0f );     // compute non-bonded interaction between i-th and j-th neighbors using Lennard-Jones and Coulomb interactions and Hydrogen bond correction
                    f1 -=  fij.xyz;
                    f2 +=  fij.xyz;
                }

                fbs[i]+= f1;
                fbs[j]+= f2;
            }
        }

    }

    // ========= Save results - store forces on atoms and recoil on its neighbors  (pi-forces are already done)
    const int i4 =(iG + iS*nnode*2 )*4;
    //const int i4p=i4+nnode*4;
    for(int i=0; i<NNEIGH; i++){
        fneigh[i4 +i] = (float4){fbs[i],0};
        //fneigh[i4p+i] = (float4){fps[i],0};
    }
    //fapos[iav     ] = (float4){fa ,0}; // If we do  run it as first forcefield
    fapos[iav       ] += (float4){fa ,0};  // If we not run it as first forcefield
    //fapos[iav+nAtoms]  = (float4){fpi,0};

}


//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void getMMFFf4_bak(
    const int4 nDOFs,               // 1   (nAtoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2  [natoms]    positions of atoms (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  fapos,        // 3  [natoms]    forces on    atoms (just node atoms are evaluated)
    __global float4*  fneigh,       // 4  [nnode*4*2] recoil forces on neighbors (and pi-orbitals)
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atoms
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q}
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 8  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvecs,        // 12 lattice vectors         for each system
    __global cl_Mat3* ilvecs,       // 13 inverse lattice vectors for each system
    __global float4*  pbc_shifts,
    const int npbc,
    const int bSubtractVdW
){

    const int iG = get_global_id (0);   // intex of atom
    const int iS = get_global_id (1);   // index of system
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);  // number of systems
    //const int iL = get_local_id  (0);
    //const int nL = get_local_size(0);
    const int nAtoms=nDOFs.x;
    const int nnode =nDOFs.y;
    const int nvec  = nAtoms+nnode;

    const int i0a   = iS*nAtoms;
    const int i0n   = iS*nnode;
    const int i0v   = iS*nvec;
    const int ipbc0 = iS*npbc;

    const int iaa = iG + i0a;
    const int ian = iG + i0n;
    const int iav = iG + i0v;


    #define NNEIGH 4

    // if(iG==0){
    //     printf( "GPU::getMMFFf4() npbc=%i \n", npbc );
    //     for(int i=0; i<npbc; i++){
    //         printf( "pbcshift[%i](%6.3f,%6.3f,%6.3f)\n", i, pbc_shifts[i].x,pbc_shifts[i].y,pbc_shifts[i].z );
    //     }
    // }

    // if(iav==0){ printf( "GPU::getMMFFf4() nnode=%i nAtoms=%i iS %i nG %i nS %i \n", nnode, nAtoms, iS, nG, nS ); }
    // if(ia==0)for(int i=0; i<nnode; i++){
    //     printf( "GPU[%i] ", i );
    //     printf( "ngs{%2i,%2i,%2i,%2i} ", neighs[i].x, neighs[i].y, neighs[i].z, neighs[i].w );
    //     printf( "apar{%6.3f,%6.3f,%6.3f,%6.3f} ", apars[i].x, apars[i].y, apars[i].z, apars[i].w );
    //     printf(  "BL{%6.3f,%6.3f,%6.3f,%6.3f} ", bLs[i].x, bLs[i].y, bLs[i].z, bLs[i].w );
    //     printf(  "BK{%6.3f,%6.3f,%6.3f,%6.3f} ", bKs[i].x, bKs[i].y, bKs[i].z, bKs[i].w );
    //     printf( "Ksp{%6.3f,%6.3f,%6.3f,%6.3f} ", Ksp[i].x, Ksp[i].y, Ksp[i].z, Ksp[i].w );
    //     printf( "Kpp{%6.3f,%6.3f,%6.3f,%6.3f} ", Kpp[i].x, Kpp[i].y, Kpp[i].z, Kpp[i].w );
    //     printf( "\n" );
    // }

    // ========= Private Memory
    //const cl_Mat3 lvec    = lvecs [iS];
    //const cl_Mat3 invLvec = ilvecs[iS];

    // ---- Dynamical
    float4  hs [4];              // direction vectors of bonds
    float3  fbs[4];              // force on neighbor sigma
    float3  fps[4];              // force on neighbor pi
    float3  fa  = float3Zero;    // force on center position
    float3  fpi = float3Zero;    // force on pi orbital

    float E=0;

    // ---- Params
    const int4   ng  = neighs[iaa];
    const int4   ngC = neighCell[iaa];
    const float3 pa  = apos[iav].xyz;
    const float3 hpi = apos[iav+nAtoms].xyz;
    const float4 par = apars[ian];    //     (xy=s0_ss,z=ssK,w=piC0 )
    const float4 vbL = bLs[ian];
    const float4 vbK = bKs[ian];
    const float4 vKs = Ksp[ian];
    const float4 vKp = Kpp[ian];

    // Temp Arrays
    const int*   ings  = (int*  )&ng;
    const int*   ingC  = (int*  )&ngC;
    const float* bL    = (float*)&vbL;
    const float* bK    = (float*)&vbK;
    const float* Kspi  = (float*)&vKs;
    const float* Kppi  = (float*)&vKp;

    const float   ssC0   = par.x*par.x - par.y*par.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

    //const int iS_DBG = 5;
    // //const int iG_DBG = 0;
    //const int iG_DBG = 0;
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     //for(int i=0; i<nAtoms; i++){ int4 ng_ = neighs[i+iS*nAtoms];  printf( "GPU[%i|%i] neighs(%i,%i,%i,%i) \n", i, iS, ng_.x,ng_.y,ng_.z,ng_.w ); };
    //     //for(int i=0; i<nAtoms*nS; i++){ int iv=i%nAtoms; int4 ng_ = neighs[i]; if(iv==0)printf("----\n"); printf( "%3i [%2i,%2i,%i] neighs(%i,%i,%i,%i) \n", i, i/nAtoms, iv, iv<=nAtoms, ng_.x,ng_.y,ng_.z,ng_.w );  };
    //     //for(int i=0; i<nvec*nS; i++){ int iv=i%nvec; float4 p = apos[i]; if(iv==0)printf("----\n"); printf( "%3i [%2i,%2i,%i,%i] p(%10.5f,%10.5f,%10.5f)  \n", i, i/nvec, iv,iv<=nnode,iv<=nAtoms, p.x,p.y,p.z );  };
    // }
    //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU[%i=%i|%i] neighs(%i,%i,%i,%i) \n", iG, iaa, iS, ings[0],ings[1],ings[2],ings[3] );

    for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }

    // ========= Evaluate Bonds
    float3 f1,f2;         // temporary forces
    for(int i=0; i<NNEIGH; i++){
        float4 h;
        const int ing  = ings[i];
        const int ingv = ing+i0v;
        const int inga = ing+i0a;
        if(ing<0) break;
        h.xyz    = apos[ingv].xyz - pa;    //printf( "[%i|%i] ing=%i h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", ia,i,ing, h.x,h.y,h.z, apos[ing].x,apos[ing].y,apos[ing].z,  pa.x,pa.y,pa.z );

        // { // PBC bond vector correction
        //     float3 u  = (float3){ dot( invLvec.a.xyz, h.xyz ), dot( invLvec.b.xyz, h.xyz ), dot( invLvec.c.xyz, h.xyz ) };
        //     h.xyz   += lvec.a.xyz*(1.f-(int)(u.x+1.5f))
        //              + lvec.b.xyz*(1.f-(int)(u.y+1.5f))
        //              + lvec.c.xyz*(1.f-(int)(u.z+1.5f));
        //     // if((iG==iG_DBG)&&(iS==iS_DBG)){
        //     //     float3 shi =  (float3){(1.f-(int)(u.x+1.5f)),  (1.f-(int)(u.y+1.5f)), (1.f-(int)(u.z+1.5f)) };
        //     //     printf( "GPU:bond[%i,%i] u(%6.3f,%6.3f,%6.3f) shi(%6.3f,%6.3f,%6.3f) \n", iG, ing, u.x,u.y,u.z,   shi.x,shi.y,shi.z );
        //     // }
        // }

        { // PBC shifts
            int ic  = ingC[i];
            h.xyz  += pbc_shifts[ipbc0+ic].xyz;
        }

        float  l = length(h.xyz);
        h.w      = 1./l;
        h.xyz   *= h.w;
        hs[i]    = h;

        float epp = 0;
        float esp = 0;

        //if((iG==iG_DBG)&&(iS==iS_DBG)) printf( "GPU:h[%i|%i=%i] l %g h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", iaa,i,ing, l, h.x,h.y,h.z, apos[ingv].x,apos[ingv].y,apos[ingv].z, pa.x,pa.y,pa.z );
        if(iG<ing){   // we should avoid  2-counting because otherwise node atoms would be computed 2x, but capping only once
            E+= evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;
            //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU bond[%i=%i|%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", iG,iaa,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );

            float kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1.e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                //E += evalPiAling( hpi, pipos[ing].xyz, kpp,  &f1, &f2 );   fpi+=f1;  fps[i]+=f2;        //   pi-alignment     (konjugation)
                epp += evalPiAling( hpi, apos[ingv+nAtoms].xyz, kpp,  &f1, &f2 );   fpi+=f1;  fps[i]+=f2;    //   pi-alignment     (konjugation)
                //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", iaa,ing, kpp, dot(hpi,apos[ingv+nAtoms].xyz), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
                //printf( "GPU:pipi[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, kpp, dot(hpi,apos[ing+nAtoms].xyz), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
                E+=epp;
                //printf( "GPU[%i|%i] hpi(%g,%g,%g) hpj(%g,%g,%g) \n", ia,ing, hpi.x,hpi.y,hpi.z, apos[ing+nAtoms].x,apos[ing+nAtoms].y,apos[ing+nAtoms].z );
            }

            // ToDo: triple bonds ?
        }


        // DEBUG: ERROR: uncomenting this couse drift
        // pi-sigma
        float ksp = Kspi[i];
        if(ksp>1.e-6){
            esp += evalAngCos( (float4){hpi,1.}, h, ksp, par.w, &f1, &f2 );   fpi+=f1; fa-=f2;  fbs[i]+=f2;    //   pi-planarization (orthogonality)
            //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU:sp:h[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g)\n", iaa,ing, ksp,par.z, dot(hpi,h.xyz), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z  );
            //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU:sp[%i|%i] ksp=%g piC0=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", iaa,ing, ksp,par.z, dot(hpi,h.xyz), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            E+=epp;
        }

        //printf( "GPU[%i|%i] esp=%g epp=%g \n", esp, epp );

    }


    //  ============== Angles
    for(int i=0; i<NNEIGH; i++){
        int ing = ings[i];
        if(ing<0) break;
        const float4 hi = hs[i];
        const int ingv = ing+i0v;
        const int inga = ing+i0a;
        for(int j=i+1; j<NNEIGH; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const int jngv = jng+i0v;
            const int jnga = jng+i0a;
            const float4 hj = hs[j];
            //printf( "[%i|%i,%i] hi(%g,%g,%g) hj(%g,%g,%g)\n", ia, i,j, hi.x,hi.y,hi.z,   hj.x,hj.y,hj.z );

            E += evalAngleCosHalf( hi, hj, par.xy, par.z, &f1, &f2 );
            //E += evalAngCos      ( hi, hj, par.z,  ssC0,  &f1, &f2 );     // angles between sigma bonds

            //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", iG,ing,jng, par.z, par.x,par.y, dot(hi.xyz,hj.xyz),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );

            fa    -= f1+f2;

            //if(bSubtractVdW)
            { // Remove vdW
                float4 REQi=REQKs[inga];   // ToDo: can be optimized
                float4 REQj=REQKs[jnga];
                float4 REQij;
                REQij.x  = REQi.x  + REQj.x;
                REQij.yz = REQi.yz * REQj.yz;
                float3 dp = (hj.xyz/hj.w) - (hi.xyz/hi.w);
                float4 fij = getLJQH( dp, REQij, 1.0f );
                //float4 fij = getLJQH( apos[ingv].xyz-apos[jngv].xyz, REQKij, 1.0f );
                f1 -=  fij.xyz;
                f2 +=  fij.xyz;
                //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", iG,ing,jng, length(dp), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
            }

            fbs[i]+= f1;
            fbs[j]+= f2;
            //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU:ANG[%i|%i,%i] fa(%g,%g,%g) fbs[%i](%g,%g,%g) fbs[%i](%g,%g,%g)\n", iG,ing,jng, fa.x,fa.y,fa.z, i,fbs[i].x,fbs[i].y,fbs[i].z,   j,fbs[j].x,fbs[j].y,fbs[j].z  );
            //if(ia==0)printf( "GPU:fa[%i](%g,%g,%g)\n", ia, fa.x,fa.y,fa.z  );
            // ToDo: subtract non-covalent interactions
        }
    }


    // ========= Save results

    const int i4 =(iG + iS*nnode*2 )*4;
    const int i4p=i4+nnode*4;
    //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "GPU i4_max %i %i \n",  ( (nG-1) + (nS-1)*nnode*2 )*4, ( (nG-1) + (nS-1)*nnode*2 )*4 + nnode*4+4 );
    //printf( "[%i,%i] i4 %i i4p %i \n", iG,iS, i4, i4p );
    for(int i=0; i<NNEIGH; i++){
        fneigh[i4 +i] = (float4){fbs[i],0};
        fneigh[i4p+i] = (float4){fps[i],0};
    }
    //fapos[iav     ] = (float4){fa ,0}; // If we do  run it as first forcefield
    fapos[iav       ] += (float4){fa ,0};  // If we not run it as first forcefield
    fapos[iav+nAtoms]  = (float4){fpi,0};


    //printf( "GPU[%i] fa(%g,%g,%g) fpi(%g,%g,%g)\n", ia, fa.x,fa.y,fa.z, fpi.x,fpi.y,fpi.z );

    //fapos[ia]=(float4){1,2,3,ia};

    //if(ia==0){ printf( "GPU::getMMFFf4() DONE\n" ); }

}

// ======================================================================
//                     updateGroups()
// ======================================================================

//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void updateGroups(
    int               ngroup,      // 1 // number of groups (total, for all systems)
    __global int2*    granges,     // 2 // (i0,n) range of indexes specifying the group
    __global int*     g2a,         // 3 // indexes of atoms corresponding to groups defined by granges
    __global float4*  apos,        // 4 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  gcenters,    // 5 // centers of each groups (CoGs)
    __global float4*  gfws,        // 6 // forwad  orietantian vector for each group
    __global float4*  gups,        // 7 // up      orietantian vector for each group
    __global float4*  gweights     // 8 // up      orietantian vector for each group
){
    const int iG = get_global_id  (0); // index of atom
    if(iG>=ngroup) return; // make sure we are not out of bounds of current system

    // if(iG==0){
    //     printf( "GPU ngroup=%i \n", ngroup );
    //     for(int i=0; i<ngroup; i++){
    //         const int2 grange = granges[i];
    //         printf("GPU granges[%i] i0=%i n=%i \n", i, grange.x, grange.y  );
    //         for(int j=0; j<grange.y; j++){
    //             int ia = g2a[ grange.x + j ];
    //             //printf( "[%i] %i \n", j, ia );
    //             printf( "GPU gweights[%i](%g,%g,%g,%g)\n", ia, gweights[ia].x,gweights[ia].y,gweights[ia].z,gweights[ia].w );
    //         }
    //         printf("\n");
    //     }
    // }

    const int2 grange = granges[iG];

    float3 cog = (float3){0.0f,0.0f,0.0f};

    float wsum = 0.f;
    for(int i=0; i<grange.y; i++){
        int ia = g2a[ grange.x + i ];
        //const float4 pe = apos[ia];
        const float4 w = gweights[ia];
        cog    += apos[ia].xyz * w.x;
        wsum   += w.x;
    }
    cog *= ( 1.f/wsum );
    gcenters[iG] = (float4){cog,0.0f};

    float3 up  = (float3){0.0f,0.0f,0.0f};
    float3 fw  = (float3){0.0f,0.0f,0.0f};
    for(int i=0; i<grange.y; i++){
        int ia = g2a[ grange.x + i ];
        //const float4 pe = apos[ia];
        const float4 w = gweights[ia];
        const float3 d = apos[ia].xyz - cog.xyz;
        fw.xyz += d * w.y;
        up.xyz += d * w.z;
    }
    {  // Orthonormalize
        fw  = normalize( fw );
        up += fw * -dot( fw, up );
        up  = normalize( up );
    }

    //printf( "GPU[iG=%i] cog(%g,%g,%g) fw(%g,%g,%g) up(%g,%g,%g) \n", iG, cog.x,cog.y,cog.z,   fw.x,fw.y,fw.z,  up.x,up.y,up.z );
    gfws[iG] = (float4){fw,0.0f};
    gups[iG] = (float4){up,0.0f};
}

// ======================================================================
//                     groupForce()
// ======================================================================

//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void groupForce(
    const int4        n,            // 1 // (natoms,nnode) dimensions of the system
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,       // 3 // forces on atoms
    __global int*     a2g,          // 4 // atom to group maping (index)
    __global float4*  gforces,      // 5 // linar forces appliaed to atoms of the group
    __global float4*  gtorqs,       // 6 // {hx,hy,hz,t} torques applied to atoms of the group
    __global float4*  gcenters,     // 7 // centers of rotation (for evaluation of the torque
    __global float4*  gfws,         // 8 // forward vector of group orientation
    __global float4*  gups,         // 9 // up      vector of group orientation
    __global float2*  gfweights    // 10 // weights for application of forces on atoms
){
    const int natoms = n.x;           // number of atoms
    const int nnode  = n.y;           // number of node atoms
    const int nGrpup = n.w;           // number of node atoms
    const int nvec   = natoms+nnode; // number of vectors (atoms+node atoms)
    const int iG = get_global_id  (0); // index of atom

    if(iG>=natoms) return; // make sure we are not out of bounds of current system

    const int iS = get_global_id  (1); // index of system
    const int nG = get_global_size(0); // number of atoms
    const int nS = get_global_size(1); // number of systems

    // if( (iG==0) && (iS==0) ){
    //     printf( "GPU::groupForce() natom=%i nnode=%i nvec=%i \n", natoms, nnode, nvec );
    // //     int ig_sel = 0;
    // //     int is = 0;
    // //     // for(int ia=0; ia<natoms; ia++){
    // //     //      int iav = ia + is*nvec;
    // //     //     printf( "%i ", a2g[iav] );
    // //     // }
    // //     // printf("\n");

    //     for(int is=0; is<nS; is++){
    //         // printf( "sys[%i] ", is );
    //         // for(int ia=0; ia<natoms; ia++){
    //         //     int iav = ia + is*nvec;
    //         //     printf( "%i ", a2g[iav] );
    //         // }
    //         // printf("\n");
    //         for(int ia=0; ia<natoms; ia++){
    //             int iav = ia + is*nvec;
    //             const int ig = a2g[iav];
    //             if(ig>=0){
    //                 //printf( "GPU:atom[%i|%i,%i] ig=%i(%i/%i) gforces(%10.6f,%10.6f,%10.6f)\n", is, ia, iav, ig, ig-is*nGrpup,nGrpup, gforces[ig].x, gforces[ig].y, gforces[ig].z  );
    //                 printf( "GPU:atom[isys=%i|ia=%i] gfweights[iav=%i](%10.6f,%10.6f) gtorqs[ig=%i](%10.6f,%10.6f,%10.6f,%10.6f)\n", is, ia,     iav,  gfweights[iav].x,gfweights[iav].y,    ig, gtorqs[ig].x, gtorqs[ig].y, gtorqs[ig].z, gtorqs[ig].w  );
    //             }
    //         }
    //     }
    // }

    //const int ian = iG + iS*nnode;
    const int iaa = iG + iS*natoms;  // index of atom in atoms array
    const int iav = iG + iS*nvec;    // index of atom in vectors array

    float4 fe    = aforce[iav]; // position of atom or pi-orbital
    const int ig = a2g[iav];  // index of the group to which this atom belongs

    float2  w = gfweights[ig];

    // --- apply linear forece from the group
    fe.xyz += gforces[ig].xyz * w.x;

    // ToDo: group vectors may be stored in Local Memory ?
    const float3 torq = gtorqs[ig].xyz;
    const float3 fw   = gfws  [ig].xyz;
    const float3 up   = gups  [ig].xyz;
    const float3 lf   = normalize( cross(fw,up) );
    const float3 tq   = fw * torq.x   +  up * torq.y    +   lf * torq.z;

    // --- apply torque from the group
    const float3 dp  = apos[iav].xyz - gcenters[ig].xyz;
    fe.xyz          += cross( dp, tq.xyz ) * w.x;

    // --- store results
    aforce[iav] = fe;

}

// ======================================================================
//                     updateAtomsMMFFf4()
// ======================================================================

/*
float2 KvaziFIREdamp( float c, float2 damp_lims, float2 clim ){
    float2 cvf;
    if      (c < clim.x ){   //-- force against veloctiy
        cvf.x = damp_lims.x; // v    // like 0.5 (strong damping)
        cvf.y = 0;           // f
    }else if(c > clim.y ){   //-- force alingned to velocity
        cvf.x = 1-damping;   // v    // like 0.99 (weak dampong damping)
        cvf.y =   damping;   // f
    }else{                   // -- force ~ perpendicular to velocity
        float f = (c-clim.x )/( clim.y - clim.x  );
        cvf.x = (1.-damping)*f;
        cvf.y =     damping *f;
    }
    return cvf;
}
*/

/*
def KvaziFIREdamp( c, clim, damps ):
    # ----- velocity & force ~ perpendicular
    t = (c-clim[0] )/( clim[1] - clim[0]  )
    cv = damps[0] + (damps[1]-damps[0])*t
    #cf =     damps[1] *t*(1-t)*4
    cf =     damps[1]*t*(1-t)*2
    # ----- velocity & force ~ against each other
    mask_lo     =  c < clim[0]
    cv[mask_lo] = damps[0]  # v    // like 0.5 (strong damping)
    cf[mask_lo] = 0             # f
    # ----- velocity & force ~ alligned
    mask_hi     =  c > clim[1]
    cv[mask_hi] = damps[1]  # v    // like 0.99 (weak dampong damping)
    cf[mask_hi] = 0           # f
    return cv,cf
*/


// Damping function for FIRE algorithm, modified to reduction of forece and velocity arrays to make it more suitable for parallelization
float2 KvaziFIREdamp( float c, float2 clim, float2 damps ){
    float2 cvf;
    if      (c < clim.x ){   //-- force against veloctiy
        cvf.x = damps.x;     // v    // like 0.5 (strong damping)
        cvf.y = 0;           // f
    }else if(c > clim.y ){   //-- force alingned to velocity
        cvf.x = damps.y;     // v    // like 0.99 (weak dampong damping)
        cvf.y = 0;           // f
    }else{                   // -- force ~ perpendicular to velocity
        float t = (c-clim.x )/( clim.y - clim.x );
        cvf.x = damps.x + (damps.y-damps.x)*t;
        cvf.y = damps.y*t*(1.f-t)*2.f;
    }
    return cvf;
}

unsigned int hash_wang(unsigned int bits) {
    //unsigned int bits = __float_as_int(value);
    bits = (bits ^ 61) ^ (bits >> 16);
    bits *= 9;
    bits = bits ^ (bits >> 4);
    bits *= 0x27d4eb2d;
    bits = bits ^ (bits >> 15);
    return bits;
}

float hashf_wang( float val, float xmin, float xmax) {
    //return ( (float)(bits)*(2147483647.0f );
    return (((float)( hash_wang(  __float_as_int(val) ) )) * 4.6566129e-10 )  *(xmax-xmin)+ xmin;
}

// Assemble recoil forces from neighbors and  update atoms positions and velocities
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void updateAtomsMMFFf4(
    const int4        n,            // 1 // (natoms,nnode) dimensions of the system
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  avel,         // 3 // velocities of atoms
    __global float4*  aforce,       // 4 // forces on atoms
    __global float4*  cvf,          // 5 // damping coefficients for velocity and force
    __global float4*  fneigh,       // 6 // recoil forces on neighbors (and pi-orbitals)
    __global int4*    bkNeighs,     // 7 // back neighbors indices (for recoil forces)
    __global float4*  constr,       // 8 // constraints (x,y,z,K) for each atom
    __global float4*  constrK,      // 9 // constraints stiffness (kx,ky,kz,?) for each atom
    __global float4*  MDparams,     // 10 // MD parameters (dt,damp,Flimit)
    __global float4*  TDrives,      // 11 // Thermal driving (T,gamma_damp,seed,?)
    __global cl_Mat3* bboxes,       // 12 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
    __global int*     sysneighs,    // 13 // // for each system contains array int[nMaxSysNeighs] of nearby other systems
    __global float4*  sysbonds      // 14 // // contains parameters of bonds (constrains) with neighbor systems   {Lmin,Lmax,Kpres,Ktens}
){
    const int natoms=n.x;           // number of atoms
    const int nnode =n.y;           // number of node atoms
    const int nMaxSysNeighs = n.w;  // max number of inter-system interactions; if <0 shwitch inter system interactions off
    const int nvec  = natoms+nnode; // number of vectors (atoms+node atoms)
    const int iG = get_global_id  (0); // index of atom

    if(iG>=nvec) return;

    const int iS = get_global_id  (1); // index of system
    const int nG = get_global_size(0); // number of atoms
    const int nS = get_global_size(1); // number of systems

    //const int ian = iG + iS*nnode;
    const int iaa = iG + iS*natoms;  // index of atom in atoms array
    const int iav = iG + iS*nvec;    // index of atom in vectors array

    const float4 MDpars  = MDparams[iS]; // (dt,damp,Flimit)
    const float4 TDrive = TDrives[iS];

    // if((iS==0)&&(iG==0)){
    //     //printf("MDpars[%i] (%g,%g,%g,%g) \n", iS, MDpars.x,MDpars.y,MDpars.z,MDpars.w);
    //     for(int is=0; is<nS; is++){
    //         //printf( "GPU::TDrives[%i](%g,%g,%g,%g)\n", i, TDrives[i].x,TDrives[i].y,TDrives[i].z,TDrives[i].w );
    //         //printf( "GPU::bboxes[%i](%g,%g,%g)(%g,%g,%g)(%g,%g,%g)\n", is, bboxes[is].a.x,bboxes[is].a.y,bboxes[is].a.z,   bboxes[is].b.x,bboxes[is].b.y,bboxes[is].b.z,   bboxes[is].c.x,bboxes[is].c.y,bboxes[is].c.z );
    //         for(int ia=0; ia<natoms; ia++){
    //             int ic = ia+is*natoms;
    //             if(constr[ia+is*natoms].w>0) printf( "GPU:sys[%i]atom[%i] constr(%g,%g,%g|%g) constrK(%g,%g,%g|%g)\n", is, ia, constr[ic].x,constr[ic].y,constr[ic].z,constr[ic].w,   constrK[ic].x,constrK[ic].y,constrK[ic].z,constrK[ic].w  );
    //         }
    //     }
    // }

    const int iS_DBG = 5; // debug system
    //const int iG_DBG = 0;
    const int iG_DBG = 1; // debug atom

    //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "updateAtomsMMFFf4() natoms=%i nnode=%i nvec=%i nG %i iS %i/%i  dt=%g damp=%g Flimit=%g \n", natoms,nnode, nvec, iS, nG, nS, MDpars.x, MDpars.y, MDpars.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     int i0a = iS*natoms;
    //     for(int i=0; i<natoms; i++){
    //         printf( "GPU:constr[%i](%7.3f,%7.3f,%7.3f |K= %7.3f) \n", i, constr[i0a+i].x,constr[i0a+i].y,constr[i0a+i].z,  constr[i0a+i].w   );
    //     }
    // }

    if(iG>=(natoms+nnode)) return; // make sure we are not out of bounds of current system

    //aforce[iav] = float4Zero;

    float4 fe      = aforce[iav]; // force on atom or pi-orbital
    const bool bPi = iG>=natoms;  // is it pi-orbital ?

    // ------ Gather Forces from back-neighbors

    int4 ngs = bkNeighs[ iav ]; // back neighbors indices

    //if(iS==5)printf( "iG,iS %i %i ngs %i,%i,%i,%i \n", iG, iS, ngs.x,ngs.y,ngs.z,ngs.w );
    //if( (iS==0)&&(iG==0) ){ printf( "GPU:fe.1[iS=%i,iG=%i](%g,%g,%g,%g) \n", fe.x,fe.y,fe.z,fe.w ); }

    // sum all recoil forces from back neighbors   - WARRNING : bkNeighs must be properly shifted on CPU by adding offset of system iS*nvec*4
    if(ngs.x>=0){ fe += fneigh[ngs.x]; } // if neighbor index is negative it means that there is no neighbor, so we skip it
    if(ngs.y>=0){ fe += fneigh[ngs.y]; }
    if(ngs.z>=0){ fe += fneigh[ngs.z]; }
    if(ngs.w>=0){ fe += fneigh[ngs.w]; }

    // ---- Limit Forces - WARRNING : Github_Copilot says: this is not the best way to limit forces, because it can lead to drift, better is to limit forces in the first forcefield run (best is NBFF)
    float Flimit = 10.0;
    float fr2 = dot(fe.xyz,fe.xyz);  // squared force
    if( fr2 > (Flimit*Flimit) ){  fe.xyz*=(Flimit/sqrt(fr2)); }  // if force is too big, we scale it down to Flimit

    // =============== FORCE DONE
    aforce[iav] = fe;             // store force before limit
    //aforce[iav] = float4Zero;   // clean force   : This can be done in the first forcefield run (best is NBFF)

    // =============== DYNAMICS

    float4 ve = avel[iav]; // velocity of atom or pi-orbital
    float4 pe = apos[iav]; // position of atom or pi-orbital

    // -------- Fixed Atoms and Bounding Box
    if(iG<natoms){                  // only atoms have constraints, not pi-orbitals
        // ------- bboxes
        const cl_Mat3 B = bboxes[iS];
        // if(B.c.x>0.0f){ if(pe.x<B.a.x){ fe.x+=(B.a.x-pe.x)*B.c.x; }else if(pe.x>B.b.x){ fe.x+=(B.b.x-pe.x)*B.c.x; }; }
        // if(B.c.y>0.0f){ if(pe.y<B.a.y){ fe.y+=(B.a.y-pe.y)*B.c.y; }else if(pe.y>B.b.y){ fe.y+=(B.b.y-pe.y)*B.c.y; }; }
        if(B.c.z>0.0f){ if(pe.z<B.a.z){ fe.z+=(B.a.z-pe.z)*B.c.z; }else if(pe.z>B.b.z){ fe.z+=(B.b.z-pe.z)*B.c.z; }; }
        // ------- constrains
        float4 cons = constr[ iaa ]; // constraints (x,y,z,K)
        if( cons.w>0.f ){            // if stiffness is positive, we have constraint
            float4 cK = constrK[ iaa ];
            cK = max( cK, (float4){0.0f,0.0f,0.0f,0.0f} );
            const float3 fc = (cons.xyz - pe.xyz)*cK.xyz;
            fe.xyz += fc; // add constraint force
            if(iS==0){printf( "GPU::constr[ia=%i|iS=%i] (%g,%g,%g|K=%g) fc(%g,%g,%g) cK(%g,%g,%g)\n", iG, iS, cons.x,cons.y,cons.z,cons.w, fc.x,fc.y,fc.z , cK.x, cK.y, cK.z ); }
        }
    }

    // -------- Inter system interactions
    if( nMaxSysNeighs>0 ){
        for(int i=0; i<nMaxSysNeighs; i++){
            const int j     = iS*nMaxSysNeighs + i;
            const int    jS = sysneighs[j];
            const float4 bj = sysbonds [j];
            const float4 pj = apos[jS*nvec + iG];
            float3 d        = pj.xyz - pe.xyz;
            float  l = length( d );
            if      (l<bj.x){
                d*=(l-bj.x)*bj.z/l;  // f = dx*kPress
            }else if(l>bj.y){
                d*=(bj.y-l)*bj.w/l;  // f = dx*kTens
            }
            fe.xyz += d;
        }
    }

    /*
    // ------ Move (kvazi-FIRE)    - proper FIRE need to reduce dot(f,v),|f|,|v| over whole system (3*N dimensions), this complicates paralell implementaion, therefore here we do it only over individual particles (3 dimensions)
    if(bPi){
        fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // subtract forces  component which change pi-orbital lenght
        ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // subtract veocity component which change pi-orbital lenght
    }
    float ff = dot(fe.xyz,fe.xyz);
	float vv = dot(ve.xyz,ve.xyz);
    float vf = dot(ve.xyz,fe.xyz);
    #define ff_safety 1e-8
    float  c          = vf/sqrt( ff*vv + 1e-8 );
    float  renorm_vf  = sqrt( vv/(ff + 1e-8) );
    //float2 cvf = KvaziFIREdamp( c, MDpars.y*0.f, (float2){-0.2f,0.2f} );
    float2 cvf = KvaziFIREdamp( c, (float2){-0.8f,-0.3f}, (float2){0.8f,0.999f} );
    //if((iG==iG_DBG)&&(iS==iS_DBG)){ printf( "GPU cos=%g cv=%g \n", c, cvf.x ); }
    ve.xyz *= cvf.x;
    ve.xyz += fe.xyz*cvf.y*renorm_vf;
    ve.xyz += fe.xyz*MDpars.x;
    pe.xyz += ve.xyz*MDpars.x;
    if(bPi){
        pe.xyz=normalize(pe.xyz);                   // normalize pi-orobitals
    }
    pe.w=0;ve.w=0;  // This seems to be needed, not sure why ?????
    avel[iav] = ve;__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond(
    const int4 ns,                  // 1
    // Dynamical
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    // Parameters
    __global float4*  REQKs,        // 4
    __global int4*    neighs,       // 5
    __global int4*    neighCell,    // 6
    __global cl_Mat3* lvecs,        // 7
    const int4        nPBC,         // 8
    const float4      GFFParams     // 9
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];

    //const int iG = get_global_id  (0);
    //const int nG = get_global_size(0);
    //const int iS = get_global_id  (1);
    //const int nS = get_global_size(1);
    //const int iL = get_local_id   (0);
    //const int nL = get_local_size (0);

    //const int natoms = ns.x;
    //const int nnode  = ns.y;
    //const int nvec   = ns.x+ns.y;

    //const int i0n = iS*nnode;
    //const int i0a   = iS*natoms;
    //const int i0v   = iS*nvec;
    //const int    iaa  = iG + i0a;
    //const int    iav  = iG + i0v;
    const bool   bPBC = (nPBC.x+nPBC.y+nPBC.z)>0;

    const int4   ng    = neighs   [ get_global_id(0) + get_global_id(1)*ns.x ];
    const int4   ngC   = neighCell[ get_global_id(0) + get_global_id(1)*ns.x ];
    const float4 REQKi = REQKs    [ get_global_id(0) + get_global_id(1)*ns.x ];
    const float3 posi  = atoms    [ get_global_id(0) + get_global_id(1)*(ns.x+ns.y) ].xyz;
    const float  R2damp = GFFParams.x*GFFParams.x;
    float4 fe          = float4Zero;

    const cl_Mat3 lvec = lvecs[ get_global_id(1) ];

    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);

    // ========= Atom-to-Atom interaction ( N-body problem )
    for (int j0=0; j0<get_global_size(0); j0+=get_local_size(0) ){
        const int i=j0+get_local_id(0);
        if(i<ns.x){
            LATOMS[get_local_id(0)] = atoms [ i + get_global_id(1)*(ns.x+ns.y) ];
            LCLJS [get_local_id(0)] = REQKs [ i + get_global_id(1)*ns.x ];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<get_local_size(0); jl++){
            const int ja=j0+jl;
            if( (ja!= get_global_id(0) ) && (ja<ns.x) ){   // ToDo: Should interact withhimself in PBC ?
                const float4 aj = LATOMS[jl];
                float4 REQK     = LCLJS [jl];
                float3 dp       = aj.xyz - posi;
                REQK.x  +=REQKi.x;
                REQK.yz *=REQKi.yz;
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));

                if(bPBC){
                    int ipbc=0;
                    dp += shift0;
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            if( !( bBonded &&(
                                    ((ja==ng.x)&&(ipbc==ngC.x))
                                    ||((ja==ng.y)&&(ipbc==ngC.y))
                                    ||((ja==ng.z)&&(ipbc==ngC.z))
                                    ||((ja==ng.w)&&(ipbc==ngC.w))
                            ))){
                                float4 fij = getLJQH( dp, REQK, R2damp );
                                fe += fij;
                            }
                            ipbc++;
                            dp    += lvec.a.xyz;
                        }
                        dp    += shift_a;
                    }
                }else
                if( !bBonded ){
                    fe += getLJQH( dp, REQK, R2damp );
                }
            }
        }
    }
    if( get_global_id(0) < ns.x ){
        forces[ get_global_id(0) + get_global_id(1)*(ns.x+ns.y) ] = fe;           // If we do  run it as first forcefield
    }
}

    apos[iav] = pe;
    */

    const bool bDrive = TDrive.y > 0.0f;

    // ------ Move (Leap-Frog)
    if(bPi){ // if pi-orbital, we need to make sure that it has unit length
        fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // subtract forces  component which change pi-orbital lenght,
        ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // subtract veocity component which change pi-orbital lenght
    }else{
        // Thermal driving  - Langevin thermostat, see C++ MMFFsp3_loc::move_atom_Langevin()
        if( bDrive ){ // if gamma>0
            fe.xyz    += ve.xyz * -TDrive.y ;  // damping,  check the untis  ... cdamp/dt = gamma
            //const float3 rnd = (float3){ hashf_wang(ve.x+TDrive.w,-1.0,1.0),hashf_wang(ve.y+TDrive.w,-1.0,1.0),hashf_wang(ve.z+TDrive.w,-1.0,1.0)};
            __private float3 ix;
            // + (float3){TDrive.w,TDrive.w,TDrive.w}
            //const float4 rnd = fract( (ve*541547.1547987f + TDrive.wwww), &ix )*2.f - (float4){1.0,1.0,1.0,1.0};  // changes every frame
            const float3 rvec = (float3){  // random vector depending on the index
                (((iG+136  + (int)(1000.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f,
                (((iG+778  + (int)(1013.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f,
                (((iG+4578 + (int)( 998.f*TDrive.w) ) * 2654435761 >> 16)&0xFF) * 0.00390625f
            };
            //const float3 rnd = fract( ( rvec + TDrive.www)*12.4565f, &ix )*2.f - (float3){1.0,1.0,1.0};
            const float3 rnd = sin( ( rvec + TDrive.www )*124.4565f );
            //if(iS==3){  printf( "atom[%i] seed=%g rvec(%g,%g,%g) rnd(%g,%g,%g) \n", iG, TDrive.w, rvec.x,rvec.y,rvec.z, rnd.x,rnd.y,rnd.z ); }
            fe.xyz    += rnd.xyz * sqrt( 2*const_kB*TDrive.x*TDrive.y/MDpars.x );
        }
    }
    cvf[iav] += (float4){ dot(fe.xyz,fe.xyz),dot(ve.xyz,ve.xyz),dot(fe.xyz,ve.xyz), 0.0f };    // accumulate |f|^2 , |v|^2  and  <f|v>  to calculate damping coefficients for FIRE algorithm outside of this kernel
    //if(!bDrive){ ve.xyz *= MDpars.z; } // friction, velocity damping
    ve.xyz *= MDpars.z;             // friction, velocity damping
    ve.xyz += fe.xyz*MDpars.x;      // acceleration
    pe.xyz += ve.xyz*MDpars.x;      // move
    //ve     *= 0.99f;              // friction, velocity damping
    //ve.xyz += fe.xyz*0.1f;        // acceleration
    //pe.xyz += ve.xyz*0.1f;        // move
    if(bPi){        // if pi-orbital, we need to make sure that it has unit length
        pe.xyz=normalize(pe.xyz);                   // normalize pi-orobitals
    }
    pe.w=0;ve.w=0;    // This seems to be needed, not sure why ?????
    avel[iav] = ve;   // store velocity
    apos[iav] = pe;   // store position


    /*
    //------ Move Gradient-Descent
    if(bPi){ fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz ); } // subtract forces  component which change pi-orbital lenght
    //pe.xyz += fe.xyz*MDpars.x*0.01f;
    //pe.xyz += fe.xyz*MDpars.x*0.01f;
    pe.xyz += fe.xyz*0.02f;
    //if(bPi){ pe.xyz=normalize(pe.xyz); }
    pe.w=0;  // This seems to be needed, not sure why ?????
    apos[iav] = pe;
    */

    //if(iG==0){ printf( "GPU::updateAtomsMMFFf4() END\n" ); }

}



// ======================================================================
//                     printOnGPU()
// ======================================================================
// Print atoms and forces on GPU
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void printOnGPU(
    const int4        n,            // 1
    const int4        mask,         // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr        // 8
){
    const int natoms=n.x;
    const int nnode =n.y;
    const int isys  =n.z;
    const int nvec  = natoms+nnode;
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);
    //const int ian = iG + iS*nnode;
    const int iaa = iG + iS*natoms;
    const int iav = iG + iS*nvec;
    //const int iG_DBG = 0;
    const int iG_DBG = 1;

    printf( "#### GPU::printOnGPU(isys=%i) natoms=%i nnode=%i nG,nS(%i,%i) \n", isys,  natoms,nnode,   nS,nG );
    if(mask.x){
        for(int i=0; i<natoms; i++){
            int ia=i + isys*nvec;
            printf( "GPU[%i=%i] ", i, ia );
            //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[ia].x, bkNeighs[ia].y, bkNeighs[ia].z, bkNeighs[ia].w );
            printf( "fapos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[ia].x, aforce[ia].y, aforce[ia].z, aforce[ia].w );
            //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[ia].x, avel[ia].y, avel[ia].z, avel[ia].w );
            printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[ia].x, apos[ia].y, apos[ia].z, apos[ia].w );
            printf(  "constr{%6.3f,%6.3f,%6.3f,%6.3f} ", constr[ia].x, constr[ia].y, constr[ia].z, constr[ia].w );
            printf( "\n" );
        }
    }
    if(mask.y){
        for(int i=0; i<nnode; i++){
            int i1=i+natoms + isys*nvec;
            printf( "GPU[%i=%i] ", i, i1 );
            printf(  "fpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[i1].x, aforce[i1].y, aforce[i1].z, aforce[i1].w );
            //printf(  "vpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z, avel[i1].w );
            printf(   "pipos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z, apos[i1].w );
            printf( "\n" );
        }
    }
    if(mask.z){
        for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
            int i0 = isys*nvec*8;
            int i1=i0 + i*4+j;
            int i2=i0 + i1 + nnode*4;
            printf( "GPU[%i,%i] ", i, j );
            printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i1].x, fneigh[i1].y, fneigh[i1].z, fneigh[i1].w );
            printf( "fneighpi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i2].x, fneigh[i2].y, fneigh[i2].z, fneigh[i2].w );
            printf( "\n" );
        }}
    }

}

// ======================================================================
//                     cleanForceMMFFf4()
// ======================================================================
// Clean forces on atoms and neighbors to prepare for next forcefield evaluation
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void cleanForceMMFFf4(
    const int4        n,           // 2
    __global float4*  aforce,      // 5
    __global float4*  fneigh       // 6
){
    const int natoms=n.x;
    const int nnode =n.y;
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);
    const int nvec = natoms+nnode;

    const int iav = iG + iS*nvec;
    const int ian = iG + iS*nnode;

    aforce[iav]=float4Zero;
    //aforce[iav]=(float4){iG,iS,iav,0.0};

    //if(iav==0){ printf("GPU::cleanForceMMFFf4() iS %i nG %i nS %i \n", iS, nG, nS );}
    //if(iG==0){ for(int i=0;i<(natoms+nnode);i++ ){printf("cleanForceMMFFf4[%i](%g,%g,%g)\n",i,aforce[i].x,aforce[i].y,aforce[i].z);} }
    if(iG<nnode){
        const int i4 = ian*4;
        fneigh[i4+0]=float4Zero;
        fneigh[i4+1]=float4Zero;
        fneigh[i4+2]=float4Zero;
        fneigh[i4+3]=float4Zero;
    }
    //if(iG==0){ printf( "GPU::updateAtomsMMFFf4() END\n" ); }
}



// ======================================================================
//                           getShortRangeBuckets()
// ======================================================================

/*
Algorithm:
 * We do pairwise interactions only for particles which are in overlaping groups
 * each group has neighborlist listing which groups overlap with it (or we can just check overlap of groups, for small systems thi may be faster)
 * work_group is identical with group_i, therefore we should make sure number of atoms in each group is ~ work_group_size
 * for all particles in group_j we check if they are enclosed by bounding box of group_i
    * if that is the case we add them to local memory
 *
*/
__kernel void getShortRangeBuckets(
    const int4 ns,                  // 1
    // Dynamical
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global int2*    buckets,      // 4 // i0,n for bucket i
    __global float8*  BBs,          // 6 // bounding boxes (xmin,xmax,ymin,0,  ymax,zmin,zmax,0 )
    // Parameters
    __global float4*  REQKs,        // 4
    const float Rcut,
    const float SRdR,
    const float SRamp
    //const int4 nPBC,              // 7
    //const cl_Mat3 lvec,           // 8
    //const float Rdamp             // 9
){
    // local size should be equal to maximum size of one bucket (i.e. maximum number of atoms in one bucket)
    __local float4 POS[16];  // atom positions
    __local float4 PAR[16];  // REQKs parameters
    __local bool   mask[16];

    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int  ib = get_group_id(0);
    const int2 bi = buckets[ib];
    if(iL>bi.y) return; // check if atom within group range
    //if(iG>=natoms) return;

    const int nb     = ns.w;
    const int natoms = ns.x;
    //const int nnode =ns.y;

    // ========= Atom-to-Atom interaction ( N-body problem )
    float4 posi  = atoms[iG];
    float4 REQKi = REQKs[iG];
    float4 fe    = float4Zero;
    float8 bbi   = BBs[ib];
    float8 bbi2  = bbi; bbi2.lo.xyz+=Rcut;  bbi2.hi.xyz-=Rcut;


    for(int jb=0; jb<nb; jb++){
        int2 b = buckets[jb];

        // --- PBC replicas ?

        // We do not do this if we make neighborlist for groups
        { // check if bbj overlaps with bbi?
            float8 bbj = BBs[jb];
            if (bbi2.hi.x < bbj.lo.x || bbi2.lo.x > bbj.hi.x ||  // Separated along x-axis?
                bbi2.hi.y < bbj.lo.y || bbi2.lo.y > bbj.hi.y ||  // Separated along y-axis?
                bbi2.hi.z < bbj.lo.z || bbi2.lo.z > bbj.hi.z)    // Separated along z-axis?
            { // No overlap
                continue; // skip this group_j
            }
        }

        int ia = b.x + iL;


        // copy atoms to local memory
        //   * we copy only those atoms which are within the bounding box
        //   * we need to know which atoms were copied, therefore we use mask[]
        mask[iL] = false;
        if( iL < b.y ){
            float4 p = atoms[ia];
            // check if the particle from group_j is inside BBox of group_i
            if( (p.x<bbi.lo.x) && (p.x>bbi.hi.x) &&
                (p.y<bbi.lo.y) && (p.y>bbi.hi.y) &&
                (p.z<bbi.lo.z) && (p.z>bbi.hi.z)
            ){
                POS[iL]  = p;
                PAR[iL]  = REQKs[ia];
                mask[iL] = true;       // we need to know if the atom is in local memory or not
            }
        }
        //mask[iL] = bIn;
        barrier(CLK_LOCAL_MEM_FENCE);

        for (int j=0; j<nL; j++){
            if( mask[j] ){
                const float4 aj = POS[j];
                const float3 dp = aj.xyz - posi.xyz;
                float4 REQK = PAR[j];
                REQK.x +=REQKi.x;
                REQK.yz*=REQKi.yz;
                float4 fij = getR4repulsion( dp, REQK.x-SRdR, REQK.x, REQK.y*SRamp );
                fe += fij;
                //if(iG==4){ printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,0, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    forces[iG] = fe;
    //forces[iG] = fe*(-1.f);

}

// ======================================================================
//                           getShortRangeBuckets2()
// ======================================================================
// This algorithm assumes all atoms which overlap with group_i were already added to group_i list


__kernel void getShortRangeBuckets2(
    const int4 ns,                  // 1
    // Dynamical
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global int2*    buckets,      // 4 // {i0,n} particles which belong to group_i
    __global int2*    bucketsJs,    // 5 // {i0,n} particles which overlap with bounding box of group_i (i.e. can be from any group_j)
    __global int*     overIndex,    //6 indexes of atoms in overlap split by  bucketsJs
    __global int*     overCell,    //6 indexes of atoms in overlap split by  bucketsJs
    // Parameters
    __global float4*  REQKs,        // 8
    __global cl_Mat3* lvecs,
    const float Rcut,
    const float SRdR,
    const float SRamp,
    const int bPBC
    //const int4 nPBC,              // 7
    //const cl_Mat3 lvec,           // 8
    //const float Rdamp             // 9
){
    // local size should be equal to maximum size of one bucket (i.e. maximum number of atoms in one bucket)
    __local float4 POS[16];  // atom positions
    __local float4 PAR[16];  // REQKs parameters
    __local int    Js [16];  // atom index
    __local int    JCs[16];  // cell index

    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int  ib = get_group_id(0);
    const int2 bi = buckets[ib];
    if(iL>bi.y) return; // check if atom within group range
    //if(iG>=natoms) return;

    const int natoms=ns.x;
    const int nnode =ns.y;
    const int nb = ns.w;             // number of buckets

    // only if bPBC=true
    const int iS = get_global_id  (1); // index of system
    const cl_Mat3 lvec = lvecs[iS];

    // ========= Atom-to-Atom interaction ( N-body problem )
    float4 posi  = atoms[iG];
    float4 REQKi = REQKs[iG];
    float4 fe    = float4Zero;
    const int2 bj = bucketsJs[ib];
    for (int j0=0; j0<bj.y; j0+=nL){
        const int j=j0+iL;
        if(j<bj.y){  // copy to local memory
            int ja  = overIndex[bj.x+j];
            POS[iL] = atoms[ja];
            PAR[iL] = REQKs[ja];
            Js [iL] = ja;
            JCs[iL] = overCell[ja];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            int ja = Js[j];
            float4 aj = POS[j];
            if(bPBC){
                const int ilvec = JCs[j];
                const int ilveca = ((ilvec&0xF0)>>4)-8;
                const int ilvecb = ((ilvec&0x0F)   )-8;
                aj += lvec.a*ilveca + lvec.a*ilvecb;
            }
            const float3 dp = aj.xyz - posi.xyz;
            float4 REQK = PAR[j];
            REQK.x +=REQKi.x;
            REQK.yz*=REQKi.yz;
            float4 fij = getR4repulsion( dp, REQK.x-SRdR, REQK.x, REQK.y*SRamp );
            fe += fij;
            //if(iG==4){ printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,0, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    };
    forces[iG] = fe;
    //forces[iG] = fe*(-1.f);
}

// ======================================================================
//                           sortAtomsToBucketOverlaps()
// ======================================================================
// This function will project all atoms into interaction bounding-box of each group
// * each iG (thread) is one group
// * groups share particles in local memory
// NOTE: sortAtomsToBucketOverlaps() does not have to be run every cycle if atoms does not move too far (similar to neighor-list)
__kernel void sortAtomsToBucketOverlaps(
    const int4 ns,                  // 1
    // Dynamical
    __global float4*  atoms,        // 2
    __global float4*  shifts,
    __global int2*    buckets,      // 4 // i0,n for bucket i
    __global float8*  BBs,          // 6 // bounding boxes (xmin,xmax,ymin,0,  ymax,zmin,zmax,0 )
    __global int2*    bucketsJs,    // 5 // {i0,n} particles which overlap with bounding box of group_i (i.e. can be from any group_j)
    __global int*     overIndex,    //6   indexes of atoms in overlap split by  bucketsJs
    __global int*     overCell,     //    index of PBC cell of the atoms in the overlpa
    //__global float4*  overParams,   //6 indexes of atoms in overlap split by  bucketsJs
    //__global float4*  overPos,      //6 indexes of atoms in overlap split by  bucketsJs
    __global cl_Mat3* lvecs,
    const int4 nPBC,
    const float Rcut
){
    // local size should be equal to maximum size of one bucket (i.e. maximum number of atoms in one bucket)
    __local int    IND[16];
    __local float3 POS[16];  // atom positions
    //__local float4 PAR[16];  // REQKs parameters

    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int  ib = get_group_id(0);
    const int2 bi = buckets[ib];
    if(iL>bi.y) return; // check if atom within group range
    //if(iG>=natoms) return;

    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int natoms=ns.x;  // number of atoms
    const int nnode =ns.y;  // number of node atoms
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)

    const int nb     = ns.w;
    const int i0v = iS*nvec;    // index of first atom in vectors array

    // ========= Atom-to-Atom interaction ( N-body problem )
    float8 bbi   = BBs[ib];

    int iB0 = bucketsJs[iG].x;

    const cl_Mat3 lvec = lvecs[iS];
    const bool bPBC = (nPBC.x+nPBC.y+nPBC.z)>0;

    // For simplicity we go over all atoms - ignoring buckets
    int nfound = 0;
    for (int j0=0; j0<nG; j0+=nL){
        const int i=j0+iL;
        if(i<natoms){
            int ja  = i+i0v;
            POS[iL] = atoms[ja].xyz;
            IND[iL] = ja;
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( ja<natoms){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float3 p = POS[jl];    // read atom position   from local memory

                if(bPBC){
                    //int ipbc=0;
                    //dp += shift0;
                    // Fixed PBC size
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        float3 dp = p + lvec.b.xyz*iy - lvec.a.xyz*nPBC.x;
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                            if( (dp.x<bbi.lo.x) && (dp.x>bbi.hi.x) &&
                                (dp.y<bbi.lo.y) && (dp.y>bbi.hi.y) &&
                                (dp.z<bbi.lo.z) && (dp.z>bbi.hi.z)
                            ){
                                int isave = iB0 + nfound;
                                overIndex[isave] = IND[jl];
                                overCell [isave] = (ix+8) + (iy+8)*16;
                                nfound++;
                            }
                            //ipbc++;
                            dp += lvec.a.xyz;
                        }
                        //dp    += lvec.a.xyz;
                    }
                }else{
                    if( (p.x<bbi.lo.x) && (p.x>bbi.hi.x) &&
                        (p.y<bbi.lo.y) && (p.y>bbi.hi.y) &&
                        (p.z<bbi.lo.z) && (p.z>bbi.hi.z)
                    ){
                        int isave = iB0 + nfound;
                        overIndex[isave] = IND[jl];
                        //= overCell [];
                        nfound++;
                    }
                }
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
}

// ======================================================================
//                           getNonBond()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond(
    const int4 ns,                  // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int4*    neighs,       // 5 // neighbors indices      ( to ignore interactions between bonded atoms )
    __global int4*    neighCell,    // 6 // neighbors cell indices ( to know which PBC image should be ignored  due to bond )
    __global cl_Mat3* lvecs,        // 7 // lattice vectors for each system
    const int4        nPBC,         // 8 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams     // 9 // Grid-Force-Field parameters
){

    // we use local memory to store atomic position and parameters to speed up calculation, the size of local buffers should be equal to local workgroup size
    //__local float4 LATOMS[2];
    //__local float4 LCLJS [2];
    //__local float4 LATOMS[4];
    //__local float4 LCLJS [4];
    //__local float4 LATOMS[8];
    //__local float4 LCLJS [8];
    //__local float4 LATOMS[16];
    //__local float4 LCLJS [16];
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LCLJS [32];   // local buffer for atom parameters
    //__local float4 LATOMS[64];
    //__local float4 LCLJS [64];
    //__local float4 LATOMS[128];
    //__local float4 LCLJS [128];

    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory

    const int natoms=ns.x;  // number of atoms
    const int nnode =ns.y;  // number of node atoms
    //const int nAtomCeil =ns.w;
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)

    //const int i0n = iS*nnode;
    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int i0v = iS*nvec;    // index of first atom in vectors array
    //const int ian = iG + i0n;
    const int iaa = iG + i0a; // index of atom in atoms array
    const int iav = iG + i0v; // index of atom in vectors array

    //const int iS_DBG = 0;
    //const int iG_DBG = 0;

    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); }
    //if((iG==iG_DBG)&&(iS==iS_DBG)){
    //    printf( "GPU::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL );
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );
    //     }
    //}


    //if( iL==0 ){ for(int i=0; i<nL; i++){  LATOMS[i]=(float4){ 10000.0, (float)i,(float)iG,(float)iS }; LCLJS[i]=(float4){ 20000.0, (float)i,(float)iG,(float)iS }; } }


    // NOTE: if(iG>=natoms) we are reading from invalid adress => last few processors produce crap, but that is not a problem
    //       importaint is that we do not write this crap to invalid address, so we put   if(iG<natoms){forces[iav]+=fe;} at the end
    //       we may also put these if(iG<natoms){ .. } around more things, but that will unnecessarily slow down other processors
    //       we need these processors with (iG>=natoms) to read remaining atoms to the local memory.

    //if(iG<natoms){
    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;  // PBC is used if any of the PBC dimensions is >0
    //const bool bPBC=false;

    const int4   ng    = neighs   [iaa];  // neighbors indices
    const int4   ngC   = neighCell[iaa];  // neighbors cell indices
    const float4 REQKi = REQKs    [iaa];  // non-bonded parameters
    const float3 posi  = atoms    [iav].xyz; // position of atom
    const float  R2damp = GFFParams.x*GFFParams.x; // squared damping radius
    float4 fe          = float4Zero;  // force on atom

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors for this system

    //if(iG==0){ printf("GPU[iS=%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f} \n", iS, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );  }

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;   // shift of PBC image 0
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                      // shift of PBC image in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                      // shift of PBC image in the outer loop
    //}

    /*
    if((iG==iG_DBG)&&(iS==iS_DBG)){
        printf( "GPU::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) bPBC=%i nPBC(%i,%i,%i)\n", natoms,nnode,nvec, nS,nG,nL, bPBC, nPBC.x,nPBC.y,nPBC.z );
        for(int i=0; i<natoms; i++){
            printf( "GPU a[%i] ", i);
            printf( "p{%6.3f,%6.3f,%6.3f} ", atoms[i0v+i].x,atoms[i0v+i].y,atoms[i0v+i].z  );
            printf( "ng{%i,%i,%i,%i} ", neighs[i0a+i].x,neighs[i0a+i].y,neighs[i0a+i].z,neighs[i0a+i].w );
            printf( "ngC{%i,%i,%i,%i} ", neighCell[i0a+i].x,neighCell[i0a+i].y,neighCell[i0a+i].z,neighCell[i0a+i].w );
            printf( "\n");
        }
    }
    */
    if((iG==0)&&(iS==0)){
        printf( "getNonBond() iG %i, iS %i, iaa %i\n", iG, iS, iaa );
    }
    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){      // loop over all atoms in the system, by chunks of size of local memory
        const int i=j0+iL;              // index of atom in local memory
        if(i<natoms){                   // j0*nL may be larger than natoms, so we need to check if we are not reading from invalid address
            LATOMS[iL] = atoms [i+i0v]; // read atom position to local memory
            LCLJS [iL] = REQKs [i+i0a]; // read atom parameters to local memory
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( (ja!=iG) && (ja<natoms) ){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float4 aj = LATOMS[jl];    // read atom position   from local memory
                float4 REQK     = LCLJS [jl];    // read atom parameters from local memory
                float3 dp       = aj.xyz - posi; // vector between atoms
                //if((iG==44)&&(iS==0))printf( "[i=%i,ja=%i/%i,j0=%i,jl=%i/%i][iG/nG/na %i/%i/%i] aj(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,ja,nG,j0,jl,nL,   iG,nG,natoms,   aj.x,aj.y,aj.z,aj.w,  REQK.x,REQK.y,REQK.z,REQK.w  );
                REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
                //if( (j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                //if( (ji==1)&&(iG==0) )printf( "2 non-bond[%i,%i] bBonded %i\n",iG,ji,bBonded );

                if(bPBC){         // ===== if PBC is used, we need to loop over all PBC images of the atom
                    int ipbc=0;   // index of PBC image
                    dp += shift0; // shift to PBC image 0
                    // Fixed PBC size
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "GPU ipbc %i(%i,%i) shift(%7.3g,%7.3g,%7.3g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z );
                            // Without these IF conditions if(bBonded) time of evaluation reduced from 61 [ms] to 51[ms]
                            if( !( bBonded &&(                     // if atoms are bonded, we do not want to calculate non-covalent interaction between them
                                    ((ja==ng.x)&&(ipbc==ngC.x))    // check if this PBC image is not the same as one of the bonded atoms
                                    ||((ja==ng.y)&&(ipbc==ngC.y))  // i.e. if ja is neighbor of iG, and ipbc is its neighbor cell index then we skip this interaction
                                    ||((ja==ng.z)&&(ipbc==ngC.z))
                                    ||((ja==ng.w)&&(ipbc==ngC.w))
                            ))){
                                //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "ipbc %i(%i,%i) shift(%g,%g,%g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z );
                                float4 fij = getLJQH( dp, REQK, R2damp );  // calculate non-bonded force between atoms using LJQH potential
                                //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); }
                                fe += fij;
                            }
                            ipbc++;
                            dp    += lvec.a.xyz;
                        }
                        dp    += shift_a;
                    }
                }else                                                  // ===== if PBC is not used, it is much simpler
                if( !bBonded ){  fe += getLJQH( dp, REQK, R2damp ); }  // if atoms are not bonded, we calculate non-bonded interaction between them
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(iG<natoms){
        //if(iS==0){ printf( "GPU::getNonBond(iG=%i) fe(%g,%g,%g,%g)\n", iG, fe.x,fe.y,fe.z,fe.w ); }
        forces[iav] = fe;           // If we do    run it as first forcefield, we can just store force (non need to clean it before in that case)
        //forces[iav] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //forces[iav] = fe*(-1.f);
    }
}




#define EXCL_MAX 16

// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
//void func_getNonBond(
__kernel void getNonBond_ex2(
    const int4        nDOFs,        // 1 // (natoms,nnode) dimensions of the system
    // Dynamical
    __global float4*  apos,         // 2 // positions of atoms  (including node atoms [0:nnode] and capping atoms [nnode:natoms] and pi-orbitals [natoms:natoms+nnode] )
    __global float4*  aforce,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQs,         // 4 // non-bonded parameters (RvdW,EvdW,QvdW,Hbond)
    __global int*     excl,         // 5 // packed sorted exclusion list ()   
    __global cl_Mat3* lvecs,        // 6 // lattice vectors for each system
    const int4        nPBC,         // 7 // number of PBC images in each direction (x,y,z)
    const float4      GFFParams
){
    __local float4 LATOMS[32];   // local buffer for atom positions
    __local float4 LCLJS [32];   // local buffer for atom parameters

    const int iG = get_global_id  (0); // index of atom
    const int nG = get_global_size(0); // number of atoms
    const int iS = get_global_id  (1); // index of system
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0); // index of atom in local memory
    const int nL = get_local_size (0); // number of atoms in local memory

    const int natoms=nDOFs.x;  // number of atoms
    const int nnode =nDOFs.y;  // number of node atoms
    const int nvec  =natoms+nnode; // number of vectors (atoms+node atoms)
    const int i0a = iS*natoms;  // index of first atom in atoms array
    const int i0v = iS*nvec;    // index of first atom in vectors array
    const int iaa = iG + i0a; // index of atom in atoms array
    const int iav = iG + i0v; // index of atom in vectors array
    
    //if(iG<natoms){
    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;  // PBC is used if any of the PBC dimensions is >0
    //const bool bPBC=false;

    const float4 REQKi  = REQs     [iaa];  // non-bonded parameters
    const float3 posi   = apos     [iav].xyz; // position of atom
    const float  R2damp = GFFParams.x*GFFParams.x; // squared damping radius
    float4 fe           = float4Zero;  // force on atom

    const cl_Mat3 lvec   = lvecs[iS]; // lattice vectors for this system
    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;   // shift of PBC image 0
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                      // shift of PBC image in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                      // shift of PBC image in the outer loop

    //const int excl_base = iaa*EXCL_MAX;
    int iex             = iaa*EXCL_MAX;
    const int iex_end   = iex + EXCL_MAX-1;
    int jex             = excl[iex];

    if((iG==0)&&(iS==0)){
        printf( "getNonBond_ex2() iG %i, iS %i, iaa %i, iex %i, iex_end %i, jex %i\n", iG, iS, iaa, iex, iex_end, jex );
    }

    // ========= Atom-to-Atom interaction ( N-body problem ), we do it in chunks of size of local memory, in order to reuse data and reduce number of reads from global memory  
    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){     // loop over all atoms in the system, by chunks of size of local memory
        const int i=j0+iL;             // index of atom in local memory
        if(i<natoms){                  // j0*nL may be larger than natoms, so we need to check if we are not reading from invalid address
            LATOMS[iL] = apos [i+i0v]; // read atom position to local memory 
            LCLJS [iL] = REQs [i+i0a]; // read atom parameters to local memory
        }
        barrier(CLK_LOCAL_MEM_FENCE);   // wait until all atoms are read to local memory
        for (int jl=0; jl<nL; jl++){    // loop over all atoms in local memory (like 32 atoms)
            const int ja=j0+jl;         // index of atom in global memory
            if( (ja!=iG) && (ja<natoms) ){   // if atom is not the same as current atom and it is not out of range,  // ToDo: Should atom interact with himself in PBC ?
                const float4 aj = LATOMS[jl];    // read atom position   from local memory
                float4 REQK     = LCLJS [jl];    // read atom parameters from local memory
                float3 dp       = aj.xyz - posi; // vector between atoms
                //if((iG==44)&&(iS==0))printf( "[i=%i,ja=%i/%i,j0=%i,jl=%i/%i][iG/nG/na %i/%i/%i] aj(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,ja,nG,j0,jl,nL,   iG,nG,natoms,   aj.x,aj.y,aj.z,aj.w,  REQK.x,REQK.y,REQK.z,REQK.w  );
                REQK.x  +=REQKi.x;   // mixing rules for vdW Radius
                REQK.yz *=REQKi.yz;  // mixing rules for vdW Energy

                if(jex!=-1){
                   if( (iex<iex_end) && ((jex&0xFFFFFF)<ja) ){ iex++; }
                   jex = excl[iex]; 
                }

                if(bPBC){         // ===== if PBC is used, we need to loop over all PBC images of the atom
                    int ipbc=0;   // index of PBC image
                    dp += shift0; // shift to PBC image 0
                    // Fixed PBC size
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            int jac = (ipbc<<24) | ja;
                            if(jex!=jac){
                                float4 fij = getLJQH( dp, REQK, R2damp );  // calculate non-bonded force between atoms using LJQH potential
                                fe += fij;
                            }
                            ipbc++; 
                            dp    += lvec.a.xyz; 
                        }
                        dp    += shift_a;
                    }
                }else {
                    if(jex!=ja){                                              // ===== if PBC is not used, it is much simpler
                        float4 fij = getLJQH( dp, REQK, R2damp ); 
                        fe += fij;
                    }
                }
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    if(iG<natoms){
        //if(iS==0){ printf( "OCL::getNonBond(iG=%i) fe(%g,%g,%g,%g)\n", iG, fe.x,fe.y,fe.z,fe.w ); }
        aforce[iav] = fe;           // If we do    run it as first forcefield, we can just store force (non need to clean it before in that case)
        //aforce[iav] += fe;        // If we don't run it as first forcefield, we need to add force to existing force
        //aforce[iav] = fe*(-1.f);
    }
}


// ======================================================================
// ======================================================================
//                           GridFF
// ======================================================================
// ======================================================================
// Calculate Grid-Force-Field (GFF) forces on atoms between rigid surface and atoms in the system, this is done by texture interpolation

// NOTE: https://registry.khronos.org/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html
// CLK_ADDRESS_REPEAT - out-of-range image coordinates are wrapped to the valid range. This address mode can only be used with normalized coordinates. If normalized coordinates are not used, this addressing mode may generate image coordinates that are undefined.

// various samplers
__constant sampler_t samp0 =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_1       =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_2       =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_nearest =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_nearest_norm =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_gff      =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_gff_norm =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_gff   =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_gff   =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

// tri-linear interpolation of texture without using hardware interpolation ( because it is inaccurate )
float4 read_imagef_trilin( __read_only image3d_t imgIn, float4 coord ){
    float4 d = (float4)(0.00666666666,0.00666666666,0.00666666666,1.0);
    float4 icoord;
    float4 fc     =  fract( coord/d, &icoord );
    icoord*=d;
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    // NOTE AMD-GPU seems to not accept CLK_NORMALIZED_COORDS_FALSE
    //return read_imagef( imgIn, sampler_2, icoord );
    //return read_imagef( imgIn, sampler_1, coord );
    return
     (( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,0.0,0.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,d.y,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,d.y,0.0,0.0) ) * fc.x )*fc.y )*mc.z
    +(( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,0.0,d.z,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,0.0,d.z,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_2, icoord+(float4)(0.0,d.y,d.z,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_2, icoord+(float4)(d.x,d.y,d.z,0.0) ) * fc.x )*fc.y )*fc.z;
};

// tri-linear interpolation of texture without using hardware interpolation ( because it is inaccurate ), and with not-normalized coordinates
float4 read_imagef_trilin_( __read_only image3d_t imgIn, float4 coord ){
    float4 icoord;
    float4 fc     =  fract( coord, &icoord );
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    // NOTE AMD-GPU seems to not accept CLK_NORMALIZED_COORDS_FALSE
    //return read_imagef( imgIn, sampler_2, icoord );
    //return read_imagef( imgIn, sampler_1, coord );
    return
     (( read_imagef( imgIn, sampler_nearest, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest, icoord+(float4)(1.0f,0.0,0.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_nearest, icoord+(float4)(0.0,1.0f,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest, icoord+(float4)(1.0f,1.0f,0.0,0.0) ) * fc.x )*fc.y )*mc.z
    +(( read_imagef( imgIn, sampler_nearest, icoord+(float4)(0.0,0.0,1.f,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest, icoord+(float4)(1.f,0.0,1.f,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_nearest, icoord+(float4)(0.0,1.f,1.f,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest, icoord+(float4)(1.f,1.f,1.f,0.0) ) * fc.x )*fc.y )*fc.z;
};

// tri-linear interpolation of texture without using hardware interpolation ( because it is inaccurate ), but using normalized coordinates
float4 read_imagef_trilin_norm( __read_only image3d_t imgIn, float4 coord ){
    float4 icoord;
    float4 fc     =  fract( coord, &icoord );
    float4 mc     = (float4)(1.0,1.0,1.0,1.0) - fc;
    // NOTE AMD-GPU seems to not accept CLK_NORMALIZED_COORDS_FALSE
    //return read_imagef( imgIn, sampler_2, icoord );
    //return read_imagef( imgIn, sampler_1, coord );
    return
     (( read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(0.0,0.0,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(1.0f,0.0,0.0,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(0.0,1.0f,0.0,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(1.0f,1.0f,0.0,0.0) ) * fc.x )*fc.y )*mc.z
    +(( read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(0.0,0.0,1.f,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(1.f,0.0,1.f,0.0) ) * fc.x )*mc.y
     +( read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(0.0,1.f,1.f,0.0) ) * mc.x
      + read_imagef( imgIn, sampler_nearest_norm, icoord+(float4)(1.f,1.f,1.f,0.0) ) * fc.x )*fc.y )*fc.z;
};

// hermite cubic spline interpolation in 1D, returns value and derivative, ys.x,y1,y2,y3 are values at node points x=0,1,2,3
float2 spline_hermite( float x, float4 ys ){
    //float4  float y0, float y1, float dy0, float dy1
    float dy0 = (ys.z-ys.x)*0.5;
    float dy1 = (ys.w-ys.y)*0.5;
    float y01 = ys.x-ys.y;
    float p2  = (-3.f*y01 -2.f*dy0 - dy1)*x;
    float p3  = ( 2.f*y01 +    dy0 + dy1)*x*x;
    return (float2){   ys.x + x*(dy0 +   p2 +   p3 ),  // value
	                  dy0 + 2.f*p2 + 3.f*p3        };  // derivative
}

// evaluate hermite cubic spline basis functions and their derivatives in 1D, spline value can be obtained as dot(ys,hspline_basis(x))
float8 hspline_basis( float x ){
	float x2   = x*x;
	float K    =  x2*(x - 1);
	float d0   =  K - x2 + x;   //      x3 -   x2
    float d1   =  K         ;   //      x3 - 2*x2 + x
    return (float8){
    // ------ values
                     d0*-0.5,
      2*K - x2 + 1 - d1*-0.5,   //    2*x3 - 3*x2 + 1
	 -2*K + x2     + d0* 0.5,   //   -2*x3 + 3*x2
                     d0*0.5,
    // ----- derivatives
                     d0*-0.5,
	     2*K - d1*-0.5,   //    6*x2 - 6*x
	    -2*K + d0* 0.5,   //   -6*x2 + 6*x
               d0*0.5
                     };
}

// tri-cubic interpolation of texture without using hardware interpolation ( because it is inaccurate ), but using normalized coordinates
float4 interpolate_tricubic( __read_only image3d_t im, float4 p0 ){
    float4 dx=(float4){1.f,0.f,0.f,0.f};
    float4 dy=(float4){0.f,1.f,0.f,0.f};
    float4 dz=(float4){0.f,0.f,1.f,0.f};
    float4 iu; float4 u = fract(p0,&iu);
    float8 cx = hspline_basis(u.x);
    float8 cy = hspline_basis(u.y);
    float4 p;
    p=iu   -dz   -dy; float2 S00= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu   -dz      ; float2 S01= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu   -dz   +dy; float2 S02= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu   -dz+dy+dy; float2 S03= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    float3 S0  = S00.xyx*cy.s04.xxy + S01.xyx*cy.s15.xxy + S02.xyx*cy.s26.xxy + S03.xyx*cy.s37.xxy;
    p=iu         -dy; float2 S10= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu            ; float2 S11= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu         +dy; float2 S12= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu      +dy+dy; float2 S13= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    float3 S1  = S10.xyx*cy.s04.xxy + S11.xyx*cy.s15.xxy + S12.xyx*cy.s26.xxy + S13.xyx*cy.s37.xxy;
    p=iu   +dz   -dy; float2 S20= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu   +dz      ; float2 S21= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu   +dz   -dy; float2 S22= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu   +dz+dy+dy; float2 S23= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    float3 S2  = S20.xyx*cy.s04.xxy + S21.xyx*cy.s15.xxy + S22.xyx*cy.s26.xxy + S23.xyx*cy.s37.xxy;
    p=iu+dz+dz   -dy; float2 S30= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu+dz+dz      ; float2 S31= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu+dz+dz   +dy; float2 S32= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    p=iu+dz+dz+dy+dy; float2 S33= read_imagef(im,samp0,p-dx ).x*cx.s04 + read_imagef(im,samp0,p).x*cx.s15 + read_imagef(im,samp0,p+dx).x*cx.s26 + read_imagef(im,samp0,p+dx+dx).x*cx.s37;
    float3 S3  = S30.xyx*cy.s04.xxy + S31.xyx*cy.s15.xxy + S32.xyx*cy.s26.xxy + S33.xyx*cy.s37.xxy;
    float8 cz = hspline_basis(u.z);
    return S0.xyzx*cz.s04.xxxy + S1.xyzx*cz.s15.xxxy + S2.xyzx*cz.s26.xxxy + S3.xyzx*cz.s37.xxxy;
}

// interpolation of texture using hardware interpolation, this is inaccurate, but fast
float4 interpFE( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef( imgIn, sampler_1, coord );
    //return coord;
}

// interpolation of texture using software interpolation, this is accurate, but slightly slower
float4 interpFE_prec( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef_trilin( imgIn, coord );
    // read_imagef( imgIn, sampler_1, coord );
    //return coord;
}






// ======================================================================
//                           getNonBond_GridFF()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system, and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by texture interpolation
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
// NOTE: this modified version of getNonBond() just by including Grid-Force-Field (GFF) forces
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF(
    const int4 ns,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW,EvdW,Q,H)
    // GridFF
    __read_only image3d_t  FE_Paul, // 10 // Grid-Force-Field (GFF) for Pauli repulsion
    __read_only image3d_t  FE_Lond, // 11 // Grid-Force-Field (GFF) for London dispersion
    __read_only image3d_t  FE_Coul, // 12 // Grid-Force-Field (GFF) for Coulomb interaction
    const cl_Mat3  diGrid,          // 13 // inverse of grid spacing
    const float4   grid_p0          // 14 // origin of the grid
    //__global cl_Mat3* bboxes      // 15 // bounding box (xmin,ymin,zmin)(xmax,ymax,zmax)(kx,ky,kz)
){
    __local float4 LATOMS[32];         // local memory chumk of positions of atoms
    __local float4 LCLJS [32];         // local memory chumk of atom parameters
    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms=ns.x;         // number of atoms in the system
    const int nnode =ns.y;         // number of nodes in the system
    const int nvec  =natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors of the system

    const int iS_DBG = 0;
    const int iG_DBG = 0;

    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); }
    //if((iG==iG_DBG)&&(iS==iS_DBG)) printf( "getNonBond_GridFF() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     printf( "GPU::getNonBond_GridFF() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL );
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );
    //     }
    // }

    //if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0; // Periodic boundary conditions if any of nPBC.x,nPBC.y,nPBC.z is non-zero
    const int4   ng    = neighs   [iaa];           // indexes of neighbors of the atom
    const int4   ngC   = neighCell[iaa];           // indexes of cells of neighbors of the atom
    const float4 REQKi = REQKs    [iaa];           // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H) of the atom
    const float3 posi  = atoms    [iav].xyz;       // position of the atom
    const float  R2damp = GFFParams.x*GFFParams.x; // damping radius for Lenard-Jones potential
    const float  alphaMorse = GFFParams.y;         // alpha parameter for Morse potential used in Grid-Force-Field (GFF)
    float4 fe           = float4Zero;              // forces on the atom

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;  // shift of the first PBC image
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                  // shift of lattice vector in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                  // shift of lattice vector in the outer loop

    // ========= Atom-to-Atom interaction ( N-body problem )     - we do it by chunks of nL atoms in order to reuse data and reduce number of global memory reads
    for (int j0=0; j0<natoms; j0+= nL ){ // loop over atoms in the system by chunks of nL atoms which fit into local memory
        const int i = j0 + iL;           // global index of atom in the system
        LATOMS[iL] = atoms [i+i0v];      // load positions  of atoms into local memory
        LCLJS [iL] = REQKs [i+i0a];      // load parameters of atoms into local memory
        barrier(CLK_LOCAL_MEM_FENCE);    // wait until all atoms are loaded into local memory
        for (int jl=0; jl<nL; jl++){     // loop over atoms in the local memory chunk
            const int ja=jl+j0;          // global index of atom in the system
            if( (ja!=iG) && (ja<natoms) ){ // atom should not interact with himself, and should be in the system ( j0*nL+iL may be out of range of natoms )
                const float4 aj   = LATOMS[jl]; // position of the atom
                float4       REQK = LCLJS [jl]; // parameters of the atom
                float3 dp   = aj.xyz - posi;    // vector between atoms
                REQK.x  +=REQKi.x;              // mixing of RvdW radii
                REQK.yz *=REQKi.yz;             // mixing of EvdW and Q
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));   // atom is bonded if it is one of the neighbors
                if(bPBC){       // ==== with periodic boundary conditions we need to consider all PBC images of the atom
                    int ipbc=0; // index of PBC image
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                    dp -= shift0;  // shift to the first PBC image
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(  // if bonded in any of PBC images, then we have to check both index of atom and index of PBC image to decide if we should skip this interaction
                                          ((ja==ng.x)&&(ipbc==ngC.x)) // check 1. neighbor and its PBC cell
                                        ||((ja==ng.y)&&(ipbc==ngC.y)) // check 2. neighbor and its PBC cell
                                        ||((ja==ng.z)&&(ipbc==ngC.z)) // ...
                                        ||((ja==ng.w)&&(ipbc==ngC.w)) // ...
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); }
                                    fe += fij; // accumulate forces
                                }
                                ipbc++;         // increment index of PBC image
                                dp+=lvec.a.xyz; // shift to the next PBC image
                            }
                            dp+=shift_a;        // shift to the next PBC image
                        }
                        dp+=shift_b;            // shift to the next PBC image
                    }
                }else{ //  ==== without periodic boundary it is much simpler, not need to care about PBC images
                    if(bBonded) continue;  // Bonded ?
                    float4 fij = getLJQH( dp, REQK, R2damp ); // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // wait until all atoms are processed, ToDo: not sure if it is needed here ?
    }

    if(iG>=natoms) return; // natoms <= nG, because nG must be multiple of nL (loccal kernel size). We cannot put this check at the beginning of the kernel, because it will break reading of atoms to local memory

    // ========== Interaction with Grid-Force-Field (GFF) ==================
    const float3 posg  = posi - grid_p0.xyz;                                                                                   // position of the atom with respect to the origin of the grid

    float4 coord = (float4)( dot(posg, diGrid.a.xyz),   dot(posg,diGrid.b.xyz), dot(posg,diGrid.c.xyz), 0.0f );          // normalized grid coordinates of the atom
    coord.z = clamp( coord.z, 0.0001f, 0.99f );  // prevent periodic boundary in z-direction
    // #if 0
        //coord +=(float4){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
        const float4 fe_Paul = read_imagef( FE_Paul, sampler_gff_norm, coord );    // interpolate Grid-Force-Field (GFF) for Pauli repulsion
        const float4 fe_Lond = read_imagef( FE_Lond, sampler_gff_norm, coord );    // interpolate Grid-Force-Field (GFF) for London dispersion
        const float4 fe_Coul = read_imagef( FE_Coul, sampler_gff_norm, coord );    // interpolate Grid-Force-Field (GFF) for Coulomb interaction
    // #else
        // const float4 fe_Paul = read_imagef_trilin_norm( FE_Paul, coord );
        // const float4 fe_Lond = read_imagef_trilin_norm( FE_Lond, coord );
        // const float4 fe_Coul = read_imagef_trilin_norm( FE_Coul, coord );
    //#endif
    //read_imagef_trilin( imgIn, coord );  // This is for higher accuracy (not using GPU hw texture interpolation)
    const float ej   = exp( alphaMorse* REQKi.x );   // exp(-alphaMorse*RvdW) pre-factor for factorized Morse potential
    const float cL   = ej*REQKi.y;                   // prefactor London dispersion (attractive part of Morse potential)
    const float cP   = ej*cL;                        // prefactor Pauli repulsion   (repulsive part of Morse potential)
    fe  += fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQKi.z;  // total GridFF force and energy on the atom (including Morse(London+Pauli) and Coulomb )

    /*
    if((iG==0)&&(iS==0)){
        printf("GPU:getNonBond_GridFF(natoms=%i)\n", natoms);
        for(int i=0; i<natoms; i++){ printf("GPU:atom[%i] apos(%6.3f,%6.3f,%6.3f|%g) rekq(%6.3f,%10.7f,%6.3f|%g) \n", i, atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   REQKs[i].x,REQKs[i].y,REQKs[i].z,REQKs[i].w ); }

        int ia0 = 0;
        float4 REQK = REQKs[ia0];
        int    np   = 50;
        float  dx   = 0.1;

        const float ej   = exp( REQK.w * REQK.x );
        const float cL   = ej*REQK.y;
        const float cP   = ej*cL;

        printf(  "GPU: grid_p0(%g,%g,%g) \n", grid_p0.x,grid_p0.y,grid_p0.z );
        printf(  "GPU: REQK Ri %g Ei %g Q %g beta %g \n", REQK.x, REQK.y, REQK.z, REQK.w );
        printf(  "GPU: ej %g  cP %g cL %g \n", ej, cP, cL );
        printf(  "#i  r.x  E_LJ E_Paul E_Lond, E_Coul  Fx_LJ Fx_Paul Fx_Lond, Fx_Coul  \n" );
        for(int i=0; i<np; i++){

            float3 pos  = (float3){  0.f, 0.f, 0.f + dx*i  };  // NOTE: p=(0.0f, 0.0f, 0.0f) is lateral center of            at the bottom of the cell
            //float3 pos  = (float3){  2.f, 2.f, 0.f + dx*i  };  // NOTE: p=(2.0f, 2.0f, 0.0f) is lateral conrer of 4x4 A cell at the bottom of the cell

            //const float3 posg  = pos - grid_p0.xyz;
            const float3 posg  = pos;
            const float4 coord = (float4)( dot(posg, diGrid.a.xyz),   dot(posg,diGrid.b.xyz), dot(posg,diGrid.c.xyz), 0.0f );
            const float4 fe_Paul = read_imagef_trilin_norm( FE_Paul, coord );
            const float4 fe_Lond = read_imagef_trilin_norm( FE_Lond, coord );
            const float4 fe_Coul = read_imagef_trilin_norm( FE_Coul, coord );
            const float4 fetot   = fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQKi.z;

            printf(  "%i %8.3f  %g %g %g %g   %g %g %g %g \n", ia0, pos.z,   fetot.w, fe_Paul.w*cP, fe_Lond.w*cL, fe_Coul.w*REQK.z,      fetot.z, fe_Paul.z*cP, fe_Lond.z*cL, fe_Coul.z*REQK.z  );

        }
    }
    */

    forces[iav] = fe;        // If we do    run it as first forcefield, in this case we do not need to clear forces before running this forcefield

}




inline int modulo(const int i, const int m) {
    int result = i % m;
    if (result < 0) {
        result += m;
    }
    return result;
}

inline int4 make_inds_pbc(const int n, const int iG) {
    switch( iG ){
        case 0: { return (int4)(0, 1,   2,   3  ); }
        case 1: { return (int4)(0, 1,   2,   3-n); }
        case 2: { return (int4)(0, 1,   2-n, 3-n); }
        case 3: { return (int4)(0, 1-n, 2-n, 3-n); }
    }
    return (int4)(-100, -100, -100, -100);
    // iqs[0] = (int4)(0, 1,   2,   3  );
    // iqs[1] = (int4)(0, 1,   2,   3-n);
    // iqs[2] = (int4)(0, 1,   2-n, 3-n);
    // iqs[3] = (int4)(0, 1-n, 2-n, 3-n);
}

inline int4 choose_inds_pbc(const int i, const int n, const int4* iqs) {
    if (i >= (n-3)) {
        const int ii = i + 4 - n;
        return iqs[ii];
    }
    return (int4)(0, +1, +2, +3);
}

inline int4 choose_inds_pbc_3( const int i, const int n, const int4* iqs ){
    if(i>=(n-3)){
        const int ii = i+4-n;
        //printf( "choose_inds_pbc() ii=%i i=%i n=%i \n", ii, i, n );
        const int4 d = iqs[ii];
        return (int4){ i+d.x, i+d.y, i+d.z, i+d.w };
    }
    return (int4){ i, i+1, i+2, i+3 };
}


inline float4 basis(float u) {
    const float inv6 = 1.0f / 6.0f;
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        inv6 * t * t * t,
        inv6 * (3.0f * u2 * (u - 2.0f) + 4.0f),
        inv6 * (3.0f * u * (1.0f + u - u2) + 1.0f),
        inv6 * u2 * u
    );
}

inline float4 dbasis(float u) {
    const float u2 = u * u;
    const float t = 1.0f - u;
    return (float4)(
        -0.5f * t * t,
        0.5f * (3.0f * u2 - 4.0f * u),
        0.5f * (-3.0f * u2 + 2.0f * u + 1.0f),
        0.5f * u2
    );
}

inline float2 fe1Dcomb(__global const float4* E, const float4 C, const float4 p, const float4 d) {
    const float4 cs = (float4)(dot(C, E[0]), dot(C, E[1]), dot(C, E[2]), dot(C, E[3]));
    return (float2)(dot(p, cs), dot(d, cs));
}

inline float3 fe2d_comb(int nz, __global const float4* E, int4 di, const float4 C, const float4 pz, const float4 dz, const float4 by, const float4 dy) {
    const float2 fe0 = fe1Dcomb(E + di.x, C, pz, dz);
    const float2 fe1 = fe1Dcomb(E + di.y, C, pz, dz);
    const float2 fe2 = fe1Dcomb(E + di.z, C, pz, dz);
    const float2 fe3 = fe1Dcomb(E + di.w, C, pz, dz);

    return (float3)(
        fe0.x * dy.x + fe1.x * dy.y + fe2.x * dy.z + fe3.x * dy.w,
        fe0.y * by.x + fe1.y * by.y + fe2.y * by.z + fe3.y * by.w,
        fe0.x * by.x + fe1.x * by.y + fe2.x * by.z + fe3.x * by.w
    );
}

inline float4 fe3d_pbc_comb(const float3 u, const int3 n, __global const float4* Es, const float4 PLQH, __local const int4* xqis, __local int4* yqis) {
    int ix = (int)u.x;
    int iy = (int)u.y;
    int iz = (int)u.z;
    if (u.x < 0) ix--;
    if (u.y < 0) iy--;
    const float tx = u.x - ix;
    const float ty = u.y - iy;
    const float tz = u.z - iz;

    if ((iz < 1) || (iz >= n.z - 2)) {
        return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }

    ix = modulo(ix-1, n.x);
    iy = modulo(iy-1, n.y);

    const int nyz = n.z * n.y;
    // int4 qx = xqis[ix%4] * nyz;
    // int4 qy = yqis[iy%4] * n.z;

    int4 qx = choose_inds_pbc( ix, n.x, xqis );
    //const int4 qx = choose_inds_pbc( ix, n.x, xqis )*nyz;
    const int4 qy = choose_inds_pbc( iy, n.y, yqis )*n.z;

    const float4 bz = basis(tz);
    const float4 dz = dbasis(tz);
    const float4 by = basis(ty);
    const float4 dy = dbasis(ty);

    const int i0 = (iz - 1) + n.z * (iy + n.y * ix);

    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) n(%i,%i,%i) \n", u.x,u.y,u.z, ix,iy,iz, n.x,n.y,n.z );
    //printf( "GPU fe3d_pbc_comb() u(%8.4f,%8.4f,%8.4f) ixyz(%i,%i,%i) qx(%i,%i,%i,%i) nyz=%i\n", u.x,u.y,u.z, ix,iy,iz, qx.x,qx.y,qx.z,qx.w, nyz );
    qx*=nyz;

    //return (float4){ 0.0f, 0.0f, 0.0f, dot(PLQH, Es[ i0 ])  };

    float3 E1 = fe2d_comb(n.z, Es + (i0 + qx.x), qy, PLQH, bz, dz, by, dy);
    float3 E2 = fe2d_comb(n.z, Es + (i0 + qx.y), qy, PLQH, bz, dz, by, dy);
    float3 E3 = fe2d_comb(n.z, Es + (i0 + qx.z), qy, PLQH, bz, dz, by, dy);
    float3 E4 = fe2d_comb(n.z, Es + (i0 + qx.w), qy, PLQH, bz, dz, by, dy);

    const float4 bx = basis(tx);
    const float4 dx = dbasis(tx);

    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))
    );
}


// ======================================================================
//                           getNonBond_GridFF_Bspline()
// ======================================================================
// Calculate non-bonded forces on atoms (icluding both node atoms and capping atoms), cosidering periodic boundary conditions
// It calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system, and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by texture interpolation
// it can be run in parallel for multiple systems, in order to efficiently use number of GPU cores (for small systems with <100 this is essential to get good performance)
// This is the most time consuming part of the forcefield evaluation, especially for large systems when nPBC>1
// NOTE: this modified version of getNonBond() just by including Grid-Force-Field (GFF) forces
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF_Bspline(
    const int4 ns,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW,EvdW,Q,H)
    // GridFF
    __global float4*  BsplinePLQ,   // 10 // Grid-Force-Field (GFF) for Pauli repulsion
    const int4     grid_ns,         // 11 // origin of the grid
    const float4   grid_invStep,    // 12 // origin of the grid
    const float4   grid_p0          // 13 // origin of the grid
){
    __local float4 LATOMS[32];         // local memory chumk of positions of atoms
    __local float4 LCLJS [32];         // local memory chumk of atom parameters
    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms=ns.x;         // number of atoms in the system
    const int nnode =ns.y;         // number of nodes in the system
    const int nvec  =natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system

    const float4 REQKi = REQKs    [iaa];           // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H) of the atom
    const float3 posi  = atoms    [iav].xyz;       // position of the atom
    float4 fe          = float4Zero;              // forces on the atom

    const int iS_DBG = 0;
    const int iG_DBG = 0;

    // =================== Non-Bonded interaction ( molecule-molecule )

    { // insulate nbff

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors of the system

    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); }
    //if((iG==iG_DBG)&&(iS==iS_DBG)) printf( "GPU::getNonBond_GridFF_Bspline() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL );
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );
    //     }
    // }

    //if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0; // Periodic boundary conditions if any of nPBC.x,nPBC.y,nPBC.z is non-zero
    const int4   ng    = neighs   [iaa];           // indexes of neighbors of the atom
    const int4   ngC   = neighCell[iaa];           // indexes of cells of neighbors of the atom

    const float  R2damp = GFFParams.x*GFFParams.x; // damping radius for Lenard-Jones potential

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;  // shift of the first PBC image
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                  // shift of lattice vector in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                  // shift of lattice vector in the outer loop

    // ========= Atom-to-Atom interaction ( N-body problem )     - we do it by chunks of nL atoms in order to reuse data and reduce number of global memory reads
    for (int j0=0; j0<natoms; j0+= nL ){ // loop over atoms in the system by chunks of nL atoms which fit into local memory
        const int i = j0 + iL;           // global index of atom in the system
        LATOMS[iL] = atoms [i+i0v];      // load positions  of atoms into local memory
        LCLJS [iL] = REQKs [i+i0a];      // load parameters of atoms into local memory
        barrier(CLK_LOCAL_MEM_FENCE);    // wait until all atoms are loaded into local memory
        for (int jl=0; jl<nL; jl++){     // loop over atoms in the local memory chunk
            const int ja=jl+j0;          // global index of atom in the system
            if( (ja!=iG) && (ja<natoms) ){ // atom should not interact with himself, and should be in the system ( j0*nL+iL may be out of range of natoms )
                const float4 aj   = LATOMS[jl]; // position of the atom
                float4       REQK = LCLJS [jl]; // parameters of the atom
                float3 dp   = aj.xyz - posi;    // vector between atoms
                REQK.x  +=REQKi.x;              // mixing of RvdW radii
                REQK.yz *=REQKi.yz;             // mixing of EvdW and Q
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));   // atom is bonded if it is one of the neighbors
                if(bPBC){       // ==== with periodic boundary conditions we need to consider all PBC images of the atom
                    int ipbc=0; // index of PBC image
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                    dp -= shift0;  // shift to the first PBC image
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(  // if bonded in any of PBC images, then we have to check both index of atom and index of PBC image to decide if we should skip this interaction
                                          ((ja==ng.x)&&(ipbc==ngC.x)) // check 1. neighbor and its PBC cell
                                        ||((ja==ng.y)&&(ipbc==ngC.y)) // check 2. neighbor and its PBC cell
                                        ||((ja==ng.z)&&(ipbc==ngC.z)) // ...
                                        ||((ja==ng.w)&&(ipbc==ngC.w)) // ...
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); }
                                    fe += fij; // accumulate forces
                                }
                                ipbc++;         // increment index of PBC image
                                dp+=lvec.a.xyz; // shift to the next PBC image
                            }
                            dp+=shift_a;        // shift to the next PBC image
                        }
                        dp+=shift_b;            // shift to the next PBC image
                    }
                }else{ //  ==== without periodic boundary it is much simpler, not need to care about PBC images
                    if(bBonded) continue;  // Bonded ?
                    float4 fij = getLJQH( dp, REQK, R2damp ); // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // wait until all atoms are processed, ToDo: not sure if it is needed here ?
    }

    } // insulate nbff

    if(iG>=natoms) return; // natoms <= nG, because nG must be multiple of nL (loccal kernel size). We cannot put this check at the beginning of the kernel, because it will break reading of atoms to local memory

    // ========== Molecule-Grid interaction with GridFF using tricubic Bspline ================== (see. kernel sample3D_comb() in GridFF.cl

    __local int4 xqs[4];
    __local int4 yqs[4];
    { // insulate gridff
        if      (iL<4){             xqs[iL]=make_inds_pbc(grid_ns.x,iL); }
        else if (iL<8){ int i=iL-4; yqs[i ]=make_inds_pbc(grid_ns.y,i ); };
        //const float3 inv_dg = 1.0f / grid_d.xyz;
        barrier(CLK_LOCAL_MEM_FENCE);

        const float ej = exp( GFFParams.y * REQKi.x ); // exp(-alphaMorse*RvdW) pre-factor for factorized Morse potential
        const float4 PLQH = (float4){
            ej*ej*REQKi.y,                   // prefactor London dispersion (attractive part of Morse potential)
            ej*   REQKi.y,                   // prefactor Pauli repulsion   (repulsive part of Morse potential)
            REQKi.z,
            0.0f
        };
        //const float3 p = ps[iG].xyz;
        const float3 u = (posi - grid_p0.xyz) * grid_invStep.xyz;

        float4 fg = fe3d_pbc_comb(u, grid_ns.xyz, BsplinePLQ, PLQH, xqs, yqs);

        //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() fg(%g,%g,%g|%g) u(%g,%g,%g) posi(%g,%g,%g) grid_invStep(%g,%g,%g)\n", fg.x,fg.y,fg.z,fg.w,  u.x,u.y,u.z, posi.x,posi.y,posi.z, grid_invStep.x, grid_invStep.y, grid_invStep.z  ); }

        fg.xyz *= -grid_invStep.xyz;
        fe += fg;
        //fes[iG] = fe;
    }  // insulate gridff

    forces[iav] = fe;        // If we do    run it as first forcefield, in this case we do not need to clear forces before running this forcefield
    //forces[iav] += fe;     // If we don't run it as first forcefield, we need to add forces to the forces calculated by previous forcefields
    //forces[iav] = fe*(-1.f);


}






// ======================================================================
//         B-spline Interpolation Functions using 3D Texture
// ======================================================================

__constant sampler_t sampler_bspline = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;


// 1D B-spline interpolation along Z, combining 4 potential components Reads 4 float4 values from texture at (ix, iy, qz[0..3]),  Combines R,G,B,A channels using C_coeffs (PLQH)m  Returns (Energy, dEnergy/duz)
inline float2 fe1Dcomb_tex(__read_only image3d_t img,
                           sampler_t smp,
                           const float4 C,             // Coefficients for combining P,L,Q,H
                           const float4 pz,            // B-spline basis functions for z (B0, B1, B2, B3)
                           const float4 dz,            // B-spline derivative basis functions for z (B'0, B'1, B'2, B'3)
                           const int ix, const int iy, // X and Y integer coordinates
                           const int4 qz)    // 4 integer Z coordinates
{
    // Combine Pauli, London, Coulomb, H-bond components using C = (P,L,Q,H) for each grid point
    const float4 cs = (float4)(
        dot(C, read_imagef(img, smp, (int4)(ix, iy, qz.x, 0))),
        dot(C, read_imagef(img, smp, (int4)(ix, iy, qz.y, 0))),
        dot(C, read_imagef(img, smp, (int4)(ix, iy, qz.z, 0))),
        dot(C, read_imagef(img, smp, (int4)(ix, iy, qz.w, 0)))
    );
    // Interpolate energy and its derivative w.r.t. u_z using B-spline basis
    return (float2)(
        dot(pz, cs), // Energy
        dot(dz, cs)  // dEnergy/du_z
    );
}

// 2D B-spline interpolation in YZ plane, for a given X-coordinate slice
// Calls fe1Dcomb_tex four times.
// Returns (dEnergy/duy, dEnergy/duz, Energy)
inline float3 fe2d_comb_tex(__read_only image3d_t img,
                            sampler_t smp,
                            const int  ix,         // Single X integer coordinate for this slice
                            const int4 qy,         // 4 integer Y coordinates
                            const int4 qz,         // 4 integer Z coordinates
                            const float4 C,        // Coefficients for combining P,L,Q,H
                            const float4 pz, const float4 dz,  // Basis for Z
                            const float4 py, const float4 dy)  // Basis for Y
{
    // Interpolate along Z for 4 different Y lines (at the given ix)
    const float2 fe0 = fe1Dcomb_tex(img, smp, C, pz, dz, ix, qy.x, qz);
    const float2 fe1 = fe1Dcomb_tex(img, smp, C, pz, dz, ix, qy.y, qz);
    const float2 fe2 = fe1Dcomb_tex(img, smp, C, pz, dz, ix, qy.z, qz);
    const float2 fe3 = fe1Dcomb_tex(img, smp, C, pz, dz, ix, qy.w, qz);

    // feN.x is Energy(yN,      uz_interp)
    // feN.y is dEnergy/duz(yN, uz_interp)

    // Interpolate along Y using results from 1D Z-interpolation
    return (float3)(
        dot(dy, (float4)(fe0.x, fe1.x, fe2.x, fe3.x)),     // dEnergy/du_y = sum_j (B'_j(u_y) * Energy(y_j, u_z_interp))
        dot(py, (float4)(fe0.y, fe1.y, fe2.y, fe3.y)),     // dEnergy/du_z = sum_j (B_j(u_y) * dEnergy/du_z(y_j, u_z_interp))
        dot(py, (float4)(fe0.x, fe1.x, fe2.x, fe3.x))      // Energy       = sum_j (B_j(u_y) * Energy(y_j, u_z_interp))
    );
}

// 3D B-spline interpolation for force and energy
// u: normalized coordinates (fractional cell coordinates)
// n: dimensions of the B-spline grid (texture dimensions)
// img: 3D texture storing (Pauli, London, Coulomb, H-bond_correction) potential values
// PLQH: coefficients to combine the 4 potential components
// xqis, yqis: precomputed PBC index patterns for X and Y dimensions
// Returns (dEnergy/dux, dEnergy/duy, dEnergy/duz, Energy)
inline float4 fe3d_pbc_comb_tex(const float3 u,
                                const int3 n,
                                __read_only image3d_t img,
                                sampler_t smp,
                                const float4 PLQH,
                                __local const int4* xqis, // Patterns from make_inds_pbc for x-dim
                                __local const int4* yqis) // Patterns from make_inds_pbc for y-dim
{
    // Integer part of u (knot index preceding the point)
    // Matches original code's floor logic for ix, iy, iz
    int ix = (int)u.x;
    int iy = (int)u.y;
    int iz = (int)u.z;
    if (u.x < 0) ix--;
    if (u.y < 0) iy--;
    if (u.z < 0) iz--; // Also apply floor logic to z

    // Fractional part of u (position within the cell defined by knot ix, iy, iz)
    const float tx = u.x - ix;
    const float ty = u.y - iy;
    const float tz = u.z - iz;

    // B-spline interpolation requires 4 knots starting from index (i-1).
    // The indices needed are (i-1, i, i+1, i+2).
    const int ix_knot_start = ix - 1;
    const int iy_knot_start = iy - 1;
    const int iz_knot_start = iz - 1;

    // Boundary condition for Z: if iz_raw_knot is too close to edge, return zero.
    // iz_raw_knot must be in [1, n.z - 3] for full 4-knot support (iz_knot_start must be >=0, iz_knot_start+3 must be < n.z).
    if ((iz < 1) || (iz >= n.z - 2)) {
        return (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    }
    // Absolute Z indices (no PBC for Z based on this check and index range)
    const int4 qz = (int4)(iz_knot_start, iz_knot_start + 1, iz_knot_start + 2, iz_knot_start + 3);

    // Apply PBC for X and Y dimensions to get base knot index for choose_inds_pbc_3
    // The base knot index for choose_inds_pbc_3 should be the starting index *after* modulo,
    // i.e., (ix-1) % n.x.
    const int ix_pbc_base = modulo(ix_knot_start, n.x);
    const int iy_pbc_base = modulo(iy_knot_start, n.y);

    // Get the 4 absolute integer grid indices for X and Y using PBC logic
    // These indices will be used directly in read_imagef
    const int4 qx = choose_inds_pbc_3(ix_pbc_base, n.x, xqis);
    const int4 qy = choose_inds_pbc_3(iy_pbc_base, n.y, yqis);

    // Calculate B-spline basis functions and their derivatives
    const float4 bz = basis(tz);  const float4 dz = dbasis(tz);
    const float4 by = basis(ty);  const float4 dy = dbasis(ty);
    const float4 bx = basis(tx);  const float4 dx = dbasis(tx);

    // Interpolate along YZ for 4 different X planes
    // E#.x = dE/duy, E#.y = dE/duz, E#.z = E, all at (qx.#, u_y_interp, u_z_interp)
    const float3 E1 = fe2d_comb_tex(img, smp, qx.x, qy, qz, PLQH, bz, dz, by, dy);
    const float3 E2 = fe2d_comb_tex(img, smp, qx.y, qy, qz, PLQH, bz, dz, by, dy);
    const float3 E3 = fe2d_comb_tex(img, smp, qx.z, qy, qz, PLQH, bz, dz, by, dy);
    const float3 E4 = fe2d_comb_tex(img, smp, qx.w, qy, qz, PLQH, bz, dz, by, dy);

    // Interpolate along X using results from 2D YZ-interpolation
    // Result is (dE/dux, dE/duy, dE/duz, E_total)
    return (float4)(
        dot(dx, (float4)(E1.z, E2.z, E3.z, E4.z)),     // dEnergy/du_x = sum_i (B'_i(u_x) * Energy(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.x, E2.x, E3.x, E4.x)),     // dEnergy/du_y = sum_i (B_i(u_x) * dEnergy/du_y(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.y, E2.y, E3.y, E4.y)),     // dEnergy/du_z = sum_i (B_i(u_x) * dEnergy/du_z(x_i, u_y_interp, u_z_interp))
        dot(bx, (float4)(E1.z, E2.z, E3.z, E4.z))      // Energy       = sum_i (B_i(u_x) * Energy(x_i, u_y_interp, u_z_interp))
    );
}


// ======================================================================
//                           getNonBond_GridFF_Bspline_tex()
// ======================================================================
// Calculate non-bonded forces on atoms (including both node atoms and capping atoms), considering periodic boundary conditions
// It calculates Lenard-Jones, Coulomb and Hydrogen-bond forces between all atoms in the system,
// and interactions of these atoms rigid surface described by Grid-Force-Field (GFF) done by tricubic B-spline texture interpolation
// It can be run in parallel for multiple systems.
// NOTE: This version uses a 3D texture for the B-spline GridFF.
__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF_Bspline_tex( // Renamed kernel to distinguish from buffer version
    const int4 ns,                  // 1 // dimensions of the system (natoms,nnode,nvec)
    // Dynamical
    __global float4*  atoms,        // 2 // positions of atoms
    __global float4*  forces,       // 3 // forces on atoms
    // Parameters
    __global float4*  REQKs,        // 4 // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H)
    __global int4*    neighs,       // 5 // indexes neighbors of atoms
    __global int4*    neighCell,    // 6 // indexes of cells of neighbor atoms
    __global cl_Mat3* lvecs,        // 7 // lattice vectors of the system
    const int4 nPBC,                // 8 // number of PBC images in each direction
    const float4  GFFParams,        // 9 // parameters of Grid-Force-Field (GFF) (RvdW_cutoff_factor_for_LJ, alphaMorse, Q_atom, H_bond_params_unused)
    // GridFF - Using Texture
    __read_only image3d_t BsplinePLQH_tex, // 10 // Grid-Force-Field (GFF) data (Pauli,London,Coulomb,HBond) in a 3D texture (Renamed texture)
    const int4     grid_ns,         // 11 // dimensions of the grid (matches buffer code name)
    const float4   grid_invStep,    // 12 // inverse of grid cell dimensions
    const float4   grid_p0          // 13 // origin of the grid
){
    __local float4 LATOMS[32];         // local memory chumk of positions of atoms
    __local float4 LCLJS [32];         // local memory chumk of atom parameters
    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms=ns.x;         // number of atoms in the system
    const int nnode =ns.y;         // number of nodes in the system
    const int nvec  =natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system

    const float4 REQKi = REQKs    [iaa];           // parameters of Lenard-Jones potential, Coulomb and Hydrogen Bond (RvdW,EvdW,Q,H) of the atom
    const float3 posi  = atoms    [iav].xyz;       // position of the atom
    float4 fe          = float4Zero;              // forces on the atom

    const int iS_DBG = 0;
    const int iG_DBG = 0;

    // =================== Non-Bonded interaction ( molecule-molecule )

    { // insulate nbff

    const cl_Mat3 lvec = lvecs[iS]; // lattice vectors of the system

    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); }
    //if((iG==iG_DBG)&&(iS==iS_DBG)) printf( "GPU::getNonBond_GridFF_Bspline() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     printf( "GPU::getNonBond_GridFF_Bspline() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL );
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );
    //     }
    // }

    //if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0; // Periodic boundary conditions if any of nPBC.x,nPBC.y,nPBC.z is non-zero
    const int4   ng    = neighs   [iaa];           // indexes of neighbors of the atom
    const int4   ngC   = neighCell[iaa];           // indexes of cells of neighbors of the atom

    const float  R2damp = GFFParams.x*GFFParams.x; // damping radius for Lenard-Jones potential

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;  // shift of the first PBC image
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);                  // shift of lattice vector in the inner loop
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);                  // shift of lattice vector in the outer loop

    // ========= Atom-to-Atom interaction ( N-body problem )     - we do it by chunks of nL atoms in order to reuse data and reduce number of global memory reads
    for (int j0=0; j0<natoms; j0+= nL ){ // loop over atoms in the system by chunks of nL atoms which fit into local memory
        const int i = j0 + iL;           // global index of atom in the system
        LATOMS[iL] = atoms [i+i0v];      // load positions  of atoms into local memory
        LCLJS [iL] = REQKs [i+i0a];      // load parameters of atoms into local memory
        barrier(CLK_LOCAL_MEM_FENCE);    // wait until all atoms are loaded into local memory
        for (int jl=0; jl<nL; jl++){     // loop over atoms in the local memory chunk
            const int ja=jl+j0;          // global index of atom in the system
            if( (ja!=iG) && (ja<natoms) ){ // atom should not interact with himself, and should be in the system ( j0*nL+iL may be out of range of natoms )
                const float4 aj   = LATOMS[jl]; // position of the atom
                float4       REQK = LCLJS [jl]; // parameters of the atom
                float3 dp   = aj.xyz - posi;    // vector between atoms
                REQK.x  +=REQKi.x;              // mixing of RvdW radii
                REQK.yz *=REQKi.yz;             // mixing of EvdW and Q
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));   // atom is bonded if it is one of the neighbors
                if(bPBC){       // ==== with periodic boundary conditions we need to consider all PBC images of the atom
                    int ipbc=0; // index of PBC image
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                    dp -= shift0;  // shift to the first PBC image
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(  // if bonded in any of PBC images, then we have to check both index of atom and index of PBC image to decide if we should skip this interaction
                                          ((ja==ng.x)&&(ipbc==ngC.x)) // check 1. neighbor and its PBC cell
                                        ||((ja==ng.y)&&(ipbc==ngC.y)) // check 2. neighbor and its PBC cell
                                        ||((ja==ng.z)&&(ipbc==ngC.z)) // ...
                                        ||((ja==ng.w)&&(ipbc==ngC.w)) // ...
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );  // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); }
                                    fe += fij; // accumulate forces
                                }
                                ipbc++;         // increment index of PBC image
                                dp+=lvec.a.xyz; // shift to the next PBC image
                            }
                            dp+=shift_a;        // shift to the next PBC image
                        }
                        dp+=shift_b;            // shift to the next PBC image
                    }
                }else{ //  ==== without periodic boundary it is much simpler, not need to care about PBC images
                    if(bBonded) continue;  // Bonded ?
                    float4 fij = getLJQH( dp, REQK, R2damp ); // calculate Lenard-Jones, Coulomb and Hydrogen-bond forces between atoms
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); }
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE); // wait until all atoms are processed, ToDo: not sure if it is needed here ?
    }

    } // insulate nbff

    if(iG>=natoms) return; // natoms <= nG, because nG must be multiple of nL (loccal kernel size). We cannot put this check at the beginning of the kernel, because it will break reading of atoms to local memory

    // ========== Molecule-Grid interaction with GridFF using tricubic Bspline (Texture based) ==================

    __local int4 xqs[4]; // Local memory for PBC index patterns for X
    __local int4 yqs[4]; // Local memory for PBC index patterns for Y

    { // insulate gridff
        // Initialize local memory for PBC index patterns. Only first 8 work-items do this.
        if      (iL<4){ xqs[iL]=make_inds_pbc(grid_ns.x,iL); }
        else if (iL<8){ yqs[iL-4]=make_inds_pbc(grid_ns.y,iL-4); };
        barrier(CLK_LOCAL_MEM_FENCE); // Ensure local memory is populated

        // Coefficients for combining Pauli, London, Coulomb, H-bond from the grid field
        // Matches original calculation using GFFParams.y (alphaMorse) and REQKi
        const float alphaMorse = GFFParams.y;
        const float ej = exp( alphaMorse * REQKi.x ); // REQKi.x is RvdW of atom_i
        const float4 PLQH = (float4){
            ej*ej*REQKi.y,                   // Pauli coeff: EvdW_i * exp(2 * alphaMorse * RvdW_i)
            ej*   REQKi.y,                   // London coeff: EvdW_i * exp(alphaMorse * RvdW_i)
            REQKi.z,                         // Coulomb coeff: Q_i
            0.0f                             // H-bond coeff (assuming zeroed out)
        };

        // Calculate normalized coordinates 'u' for B-spline interpolation
        const float3 u = (posi - grid_p0.xyz) * grid_invStep.xyz;

        // Perform 3D B-spline interpolation using texture
        // fg contains (dE/dux, dE/duy, dE/duz, Energy)
        float4 fg = fe3d_pbc_comb_tex(u, grid_ns.xyz, BsplinePLQH_tex, sampler_bspline, PLQH, xqs, yqs);

        // Convert derivatives from du space to real space: Fx = -dE/dx = - (dE/dux) * (dux/dx)
        // dux/dx = grid_invStep.x, etc.
        fg.xyz *= -grid_invStep.xyz; // grid_invStep components are positive

        fe += fg; // Add GridFF force and energy to atom's total
        // fes[iG] = fe; // If you have a separate energy buffer
    }  // insulate gridff

    // Store the total force and energy for this atom
    forces[iav] = fe;
    // Use forces[iav] += fe; if forces buffer accumulates from multiple kernels
}




// ======================================================================
//                           sampleGridFF()
// ======================================================================
// this is just to test interpolation of Grid-Force-Field (GFF) on GPU
//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void sampleGridFF(
    const int4 ns,                  // 1
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global float4*  REQs,         // 4
    const float4  GFFParams,        // 5
    __read_only image3d_t  FE_Paul, // 6
    __read_only image3d_t  FE_Lond, // 7
    __read_only image3d_t  FE_Coul, // 8
    const cl_Mat3  diGrid,          // 9
    const float4   grid_p0          // 10
){
    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int np = ns.x;

    float3 dz = (float3){ 0.0f, 0.0f, 0.1f };

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const float4 REQ        = REQs[iG];
    const float3 posi       = atoms[iG].xyz;
    const float  R2damp     = GFFParams.x*GFFParams.x;
    const float  alphaMorse = GFFParams.y;

    const float ej   = exp( alphaMorse* REQ.x );
    const float cL   = ej*REQ.y;
    const float cP   = ej*cL;

    /*
    if(iG==0){ printf( "GPU::sampleGridFF() np=%i R2damp=%g aMorse=%g p(%g,%g,%g) REQ(%g,%g,%g)  cP=%g cL=%g ej=%g \n", np, R2damp, alphaMorse, posi.x,posi.y,posi.z, REQ.x,REQ.y,REQ.z, cP,cL,ej  ); }
    if(iG==0){
        printf( "GPU_sGFF #i  z  E_Paul Fz_Paul   E_Lond Fz_Lond   E_Coul Fz_Coul  \n" );
        for(int i=0; i<np; i++){
            const float4 REQ  = REQs[i];
            //const float3 posi = atoms[i].xyz;
            const float3 posi = grid_p0.xyz + dz*i;
            const float ej   = exp( alphaMorse* REQ.x );
            const float cL   = ej*REQ.y;
            const float cP   = ej*cL;

            float4 fe          = float4Zero;
            const float3 posg  = posi - grid_p0.xyz;
            const float4 coord = (float4)( dot(posg, diGrid.a.xyz),   dot(posg,diGrid.b.xyz), dot(posg,diGrid.c.xyz), 0.0f );
            #if 0
                //coord +=(float4){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
                const float4 fe_Paul = read_imagef( FE_Paul, sampler_gff_norm, coord );
                const float4 fe_Lond = read_imagef( FE_Lond, sampler_gff_norm, coord );
                const float4 fe_Coul = read_imagef( FE_Coul, sampler_gff_norm, coord );
            #else
                const float4 fe_Paul = read_imagef_trilin_norm( FE_Paul, coord );
                const float4 fe_Lond = read_imagef_trilin_norm( FE_Lond, coord );
                const float4 fe_Coul = read_imagef_trilin_norm( FE_Coul, coord );
            #endif
            //read_imagef_trilin( imgIn, coord );  // This is for higher accuracy (not using GPU hw texture interpolation)
            fe  += fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQ.z;
            //printf( "GPU[%i] z(%g) E,fz(%g,%g)  PLQ(%g,%g,%g) REQ(%g,%g) \n", i, posi.z,  fe.w,fe.z,  cP,cL,REQ.z,  REQ.x,REQ.y  );
            printf(  "GPU_sGFF %3i %8.3f    %14.6f %14.6f    %14.6f %14.6f    %14.6f %14.6f\n", i, posi.z, fe_Paul.w,fe_Paul.z, fe_Lond.w,fe_Lond.z,  fe_Coul.w,fe_Coul.z  );
        }
    }
    */


// NOTE: https://registry.khronos.org/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html
// CLK_ADDRESS_REPEAT - out-of-range image coordinates are wrapped to the valid range. This address mode can only be used with normalized coordinates. If normalized coordinates are not used, this addressing mode may generate image coordinates that are undefined.

    // ========== Interaction with grid
    float4 fe               = float4Zero;
    const float3 posg  = posi - grid_p0.xyz;
    float4 coord = (float4)( dot(posg, diGrid.a.xyz),   dot(posg,diGrid.b.xyz), dot(posg,diGrid.c.xyz), 0.0f );
    if(iG==0){ printf( "coord(%g,%g,%g) pos(%g,%g,%g) diGrid.a(%g,%g,%g)\n", coord.x,coord.y,coord.z,  posi.x,posi.y,posi.z, diGrid.a.x,diGrid.a.y,diGrid.a.z ); }
    //#if 0
        //coord +=(float4){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
        const float4 fe_Paul = read_imagef( FE_Paul, sampler_gff_norm, coord );
        const float4 fe_Lond = read_imagef( FE_Lond, sampler_gff_norm, coord );
        const float4 fe_Coul = read_imagef( FE_Coul, sampler_gff_norm, coord );
    // #else
    //     const float4 fe_Paul = read_imagef_trilin_norm( FE_Paul, coord );
    //     const float4 fe_Lond = read_imagef_trilin_norm( FE_Lond, coord );
    //     const float4 fe_Coul = read_imagef_trilin_norm( FE_Coul, coord );
    //#endif
    //read_imagef_trilin( imgIn, coord );  // This is for higher accuracy (not using GPU hw texture interpolation)
    forces[iG] = fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQ.z;
}


// ======================================================================
//                           make_GridFF()
// ======================================================================

__attribute__((reqd_work_group_size(32,1,1)))
__kernel void make_GridFF(
    const int nAtoms,                // 1
    __global float4*  atoms,         // 2
    __global float4*  REQs,          // 3
    __write_only image3d_t  FE_Paul, // 4
    __write_only image3d_t  FE_Lond, // 5
    __write_only image3d_t  FE_Coul, // 6
    const int4     nPBC,             // 7
    const int4     nGrid,            // 8
    const cl_Mat3  lvec,             // 9
    const float4   grid_p0,          // 10
    const float4   GFFParams         // 11
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x;
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab;

    const float  alphaMorse = GFFParams.y;
    const float  R2damp     = GFFParams.x*GFFParams.x;
    const float3 dGrid_a = lvec.a.xyz*(1.f/(float)nGrid.x);
    const float3 dGrid_b = lvec.b.xyz*(1.f/(float)nGrid.y);
    const float3 dGrid_c = lvec.c.xyz*(1.f/(float)nGrid.z);
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz)

    /*
    if(iG==0){printf("GPU:make_GridFF() nL=%i,nG=%i,nAtoms=%i,nPBC(%i,%i,%i) Rdamp %g alphaMorse %g \n", nL, nG, nAtoms, nPBC.x,nPBC.y,nPBC.z, GFFParams.x, alphaMorse );}
    if(iG==0){printf("GPU:make_GridFF() p0{%6.3f,%6.3f,%6.3f} lvec{{%6.3f,%6.3f,%6.3f},{%6.3f,%6.3f,%6.3f},{%6.3f,%6.3f,%6.3f}} \n", grid_p0.x,grid_p0.y,grid_p0.z,  lvec.a.x,lvec.a.y,lvec.a.z, lvec.b.x,lvec.b.y,lvec.b.z, lvec.c.x,lvec.c.y,lvec.c.z );}
    //if(iG==0){printf("GPU::make_GridFF(nAtoms=%i) \n", nAtoms );}
    if(iG==0){
        printf( "GPU_GFF_z #i   z  Ep_Paul Fz_Paul   Ep_Lond Fz_Lond  E_Coul Fz_Coul\n");
        for(int ic=0; ic<nGrid.z; ic++){
            const float3 pos_    = grid_p0.xyz  + dGrid_a.xyz*ia      + dGrid_b.xyz*ib      + dGrid_c.xyz*ic;  // grid point within cell
            const float3 pos     = pos_ + lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;       // most negative PBC-cell
            //const float3  shift0 = lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
            float4 fe_Paul = float4Zero;
            float4 fe_Lond = float4Zero;
            float4 fe_Coul = float4Zero;
            for (int ja=0; ja<nAtoms; ja++ ){
                const float4 REQ  =       REQs[ja];
                float3       dp   = pos - atoms[ja].xyz;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                            float  r2  = dot(dp,dp);
                            float  r   = sqrt(r2 + 1e-32 );
                            float ir2  = 1.f/(r2+R2damp  );
                            // ---- Electrostatic
                            float   E  = COULOMB_CONST*REQ.z*sqrt(ir2);
                            fe_Coul   += (float4)(dp*(E*ir2), E );
                            // ---- Morse ( Pauli + Dispersion )
                            float    e = exp( -alphaMorse*(r-REQ.x) );
                            float   eM = REQ.y*e;
                            float   de = 2.f*alphaMorse*eM/r;
                            float4  fe = (float4)( dp*de, eM );
                            fe_Paul += fe * e;
                            fe_Lond += fe * (float4)( -1.0f,-1.0f,-1.0f, -2.0f );
                            dp   +=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                    }
                    dp   +=shift_c;
                }
            }
            //printf(  "FE(RvdW[%i]) Paul(%g,%g,%g|%g) Lond(%g,%g,%g|%g) Coul(%g,%g,%g|%g)  \n", ia0, fe_Paul.x,fe_Paul.y,fe_Paul.z,fe_Paul.w,   fe_Lond.x,fe_Lond.y,fe_Lond.z,fe_Lond.w,    fe_Coul.x,fe_Coul.y,fe_Coul.z,fe_Coul.w  );
            //printf(  "%i %8.3f  %g %g    %g %g    %g %g  \n", ia0, dp.x, fe_Paul.x,fe_Paul.w,   fe_Lond.x,fe_Lond.w,    fe_Coul.x,fe_Coul.w  );
            //printf(  "%i %8.3f  %g %g %g %g %g   %g %g %g %g %g  \n", ia0, dp.x,  ELJ, fetot.w, fe_Paul.w,fe_Lond.w,fe_Coul.w*REQK.z,   FLJ, fetot.x, fe_Paul.x,fe_Lond.x,fe_Coul.x*REQK.z  );
            printf(  "GPU_GFF_z %3i %8.3f    %14.6f %14.6f    %14.6f %14.6f    %14.6f %14.6f\n", ic, pos.z, fe_Paul.w,fe_Paul.z, fe_Lond.w,fe_Lond.z,  fe_Coul.w,fe_Coul.z  );
        }
    }
    */

    const int nMax = nab*nGrid.z;
    if(iG>=nMax) return;

    const float3 pos    = grid_p0.xyz  + dGrid_a.xyz*ia      + dGrid_b.xyz*ib      + dGrid_c.xyz*ic       // grid point within cell
                                       +  lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell

    //const float3  shift0 = lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
    float4 fe_Paul = float4Zero;
    float4 fe_Lond = float4Zero;
    float4 fe_Coul = float4Zero;
    for (int j0=0; j0<nAtoms; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = REQs [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( ja<nAtoms ){
                const float4 REQK =       LCLJS [jl];
                float3       dp   = pos - LATOMS[jl].xyz;

                //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                //dp+=lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;

                //float3 shift=shift0;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){

                            //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc[%i,%i,%i] dp(%g,%g,%g)\n", ix,iy,iz, dp.x,dp.y,dp.z );
                            float  r2  = dot(dp,dp);
                            float  r   = sqrt(r2+1e-32 );
                            float ir2  = 1.f/(r2+R2damp);
                            // ---- Electrostatic
                            float   E  = COULOMB_CONST*REQK.z*sqrt(ir2);
                            fe_Coul   += (float4)(dp*(E*ir2), E );
                            // ---- Morse ( Pauli + Dispersion )
                            float    e = exp( -alphaMorse*(r-REQK.x) );
                            float   eM = REQK.y*e;
                            float   de = 2.f*alphaMorse*eM/r;
                            float4  fe = (float4)( dp*de, eM );
                            fe_Paul += fe * e;
                            fe_Lond += fe * (float4)( -1.0f,-1.0f,-1.0f, -2.0f );

                            // if((iG==0)&&(j==0)){
                            //     //float3 sh = dp - pos + LCLJS[j].xyz + lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
                            //     float3 sh = shift;
                            //     printf( "GPU(%2i,%2i,%2i) sh(%7.3f,%7.3f,%7.3f)\n", ix,iy,iz, sh.x,sh.y,sh.z  );
                            // }
                            //ipbc++;

                            dp   +=lvec.a.xyz;
                            //shift+=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                        //shift+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp   +=shift_c;
                    //shift+=shift_c;
                    //dp+=lvec.b.xyz*(nPBC.y*-2.f-1.f);
                    //dp+=lvec.c.xyz;
                }

            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    if(iG>=nMax) return;
    int4 coord = (int4){ia,ib,ic,0};
    //write_imagef( FE_Paul, coord, (float4){pos,(float)iG} );
    write_imagef( FE_Paul, coord, fe_Paul );
    write_imagef( FE_Lond, coord, fe_Lond );
    write_imagef( FE_Coul, coord, fe_Coul );
}

// ======================================================================
//                           getSurfMorse()
// ======================================================================
// This is brute-force alternative to getNonBond_GridFF_Bspline - it describe interactin of molecule with substrate by pairwise interactins with multiple replicas

__kernel void getSurfMorse(
    const int4 ns,                // 1
    __global float4*  atoms,      // 2
    __global float4*  REQs,       // 3
    __global float4*  forces,     // 4
    __global float4*  atoms_s,    // 5
    __global float4*  REQ_s,      // 6
    const int4     nPBC,          // 7
    const cl_Mat3  lvec,          // 8
    const float4   pos0,          // 9
    const float4   GFFParams      // 10
){

    __local float4 LATOMS[32];
    __local float4 LCLJS [32];

    const int nAtoms  = ns.x;

    const int iG = get_global_id  (0); // index of atom in the system
    const int iS = get_global_id  (1); // index of system
    const int iL = get_local_id   (0); // index of atom in the local memory chunk
    const int nG = get_global_size(0); // total number of atoms in the system
    const int nS = get_global_size(1); // total number of systems
    const int nL = get_local_size (0); // number of atoms in the local memory chunk

    const int natoms  = ns.x;         // number of atoms in the system
    const int nnode   = ns.y;         // number of nodes in the system
    const int nvec    = natoms+nnode; // number of vectos (atoms and pi-orbitals) in the system
    const int na_surf = ns.z;         //

    //const int i0n = iS*nnode;    // index of the first node in the system
    const int i0a = iS*natoms;     // index of the first atom in the system
    const int i0v = iS*nvec;       // index of the first vector (atom or pi-orbital) in the system
    //const int ian = iG + i0n;    // index of the atom in the system
    const int iaa = iG + i0a;      // index of the atom in the system
    const int iav = iG + i0v;      // index of the vector (atom or pi-orbital) in the system

    // if( (iG==0) && (iS==0) ){
    //     printf("GPU::getSurfMorse() nglob(%i,%i) nloc(%i) ns(%i,%i,%i) \n", nG,nS, nL, ns.x,ns.y,ns.z  );
    //     //printf("GPU::getSurfMorse() nglob(%i,%i) nloc(%i) ns(%i,%i,%i) nPBC(%i,%i,%i)\n", nG,nS, nL, ns.x,ns.y,ns.z,  nPBC.x,nPBC.y,nPBC.z  );
    //     //for(int i=0; i<ns.x; i++){ printf( "forces[%i] (%g,%g,%g) \n", i, forces[i].x,forces[i].y,forces[i].z );   }
    //     for(int i=0; i<na_surf; i++){ printf( "GPU.atoms_s[%i](%g,%g,%g) \n", i, atoms_s[i].x,atoms_s[i].y,atoms_s[i].z );   }
    // }
    float4 fe   = (float4){0.0f,0.0f,0.0f,0.0f};

    if(iG>=nAtoms) return;

    const float  K          = -GFFParams.y;
    const float  R2damp     =  GFFParams.x*GFFParams.x;
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz)

    const float3 pos  = atoms[iav].xyz - pos0.xyz +  lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell
    const float4 REQi = REQs [iaa];

    // ToDo: perhaps it is efficient to share surface along isys direction ( all system operate with the same surface atoms )

    for (int j0=0; j0<na_surf; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms_s[i];
        LCLJS [iL] = REQ_s  [i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( ja<na_surf ){
                float4 REQH =       LCLJS [jl];
                float3 dp   = pos - LATOMS[jl].xyz;
                REQH.x   += REQi.x;
                REQH.yzw *= REQi.yzw;
                //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
                //dp+=lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
                //float3 shift=shift0;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                            const float4 fej = getMorseQH( dp,  REQH, K, R2damp );
                            fe -= fej;
                            //if( (iG==0) && (iS==0) && (iz==0)&&(iy==0)&&(ix==0)){
                            //if( (iG==0) && (iS==0) && (jl==0) ){
                            // if( (iG==0) && (jl==0)   && (iz==0)&&(iy==0)&&(ix==0)  ){
                            //    printf( "GPU[%3i|%3i/%3i] dp(%10.6f,%10.6f,%10.6f) LATOMS[%i](%10.6f,%10.6f,%10.6f) pos(%10.6f,%10.6f,%10.6f) pos0(%10.6f,%10.6f,%10.6f)  \n", ja,0,0,   dp.x,dp.y,dp.z,  jl,  LATOMS[jl].x,LATOMS[jl].y,LATOMS[jl].z,   pos.x,pos.y,pos.z,  pos0.x,pos0.y,pos0.z );
                            //    //printf( "GPU[%3i(%3i,%3i,%3i)] K,R2damp(%10.6f,%10.6f) l=%10.6f dp(%10.6f,%10.6f,%10.6f) REQij(%10.6f,%10.6f,%10.6f,%10.6f) fij(%10.6f,%10.6f,%10.6f) \n", ja, ix,iy,iz,  K,R2damp, length(dp), dp.x,dp.y,dp.z,  REQH.x, REQH.y, REQH.z, REQH.w,  fej.x,fej.y,fej.z );
                            // }
                            dp   +=lvec.a.xyz;
                            //shift+=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                        //shift+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp   +=shift_c;
                    //shift+=shift_c;
                    //dp+=lvec.b.xyz*(nPBC.y*-2.f-1.f);
                    //dp+=lvec.c.xyz;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    // if( (iG==0) && (iS==0) ){
    //      printf( "GPU[iG=%i,iS=%i] fe(%10.6f,%10.6f,%10.6f)\n", iG,iS, fe.x,fe.y,fe.z );
    // }

    // if( (iG==0) ){
    //       printf( "GPU[iG=%i,iS=%i] fe(%10.6f,%10.6f,%10.6f)\n", iG,iS, fe.x,fe.y,fe.z );
    // }
    forces[iav] += fe;

}



// ======================================
// ======================================
//               PP-AFM
// ======================================
// ======================================

// ======================================
//              PPAFM_makeFF()
// ======================================

__attribute__((reqd_work_group_size(32,1,1)))
__kernel void PPAFM_makeFF(
    const int4 ns,                   // 1
    __global float4*  atoms,         // 2
    __global float4*  REQs,          // 3
    __write_only image3d_t  imgOut,  // 4
    const int4     nPBC,             // 5
    const int4     nGrid,            // 6
    const cl_Mat3  lvec,             // 7
    const float4   grid_p0,          // 8
    const float4 tipParams,          // 9
    const float4 Qs,                 // 10
    const float4 QZs                 // 11
){
    __local float4 LPOS[32];
    __local float4 LREQ[32];
    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x;
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab;
    const int nMax = nab*nGrid.z;

    const int natoms = ns.x;
    const int nnode  = ns.y;
    const int iSys   = ns.z;
    const int i0a    = iSys*natoms;
    const int i0v    = iSys*(natoms+nnode);
    if(iG>=nMax) return;

    if(iG==0){printf("GPU:PPAFM_makeFF() nL=%i,nG=%i,nMax=%i,na=%i,nnd=%i,nGrid(%i,%i,%i), nPBC(%i,%i,%i) \n", nL, nG,nMax, natoms,nnode, nGrid.x,nGrid.y,nGrid.z, nPBC.x,nPBC.y,nPBC.z );}
    if(iG==0){printf("GPU:make_GridFF() p0{%6.3f,%6.3f,%6.3f} lvec{{%6.3f,%6.3f,%6.3f},{%6.3f,%6.3f,%6.3f},{%6.3f,%6.3f,%6.3f}} grid_p0{%6.3f,%6.3f,%6.3f} \n", grid_p0.x,grid_p0.y,grid_p0.z,  lvec.a.x,lvec.a.y,lvec.a.z, lvec.b.x,lvec.b.y,lvec.b.z, lvec.c.x,lvec.c.y,lvec.c.z, grid_p0.x,grid_p0.y,grid_p0.z  );}

    const float3 dGrid_a = lvec.a.xyz*(1.f/(float)nGrid.x);
    const float3 dGrid_b = lvec.b.xyz*(1.f/(float)nGrid.y);
    const float3 dGrid_c = lvec.c.xyz*(1.f/(float)nGrid.z);
    const float3 pos     = grid_p0.xyz  + dGrid_a.xyz*ia     + dGrid_b.xyz*ib      + dGrid_c.xyz*ic       // grid point within cell
                                        + lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell
    //float3 pos = grid_p0.xyz + dlvec.a.xyz*ia + dlvec.b.xyz*ib  + dlvec.c.xyz*ic;
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz)

    float4 fe  = float4Zero;

    //Qs *= COULOMB_CONST;
    for (int i0=0; i0<natoms; i0+= nL ){
        int i = i0 + iL;
        //if(i>=natoms) break;  // wrong !!!!
        LPOS[iL] = atoms[i0v+i];
        LREQ[iL] = REQs [i0a+i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<natoms ){
                //const float4 pi  = LPOS[j];
                float3       dp   = pos - LPOS[j].xyz;
                float4       REQ = LREQ[j];
                float4 Qs_ = Qs * REQ.z;
                if(iG==0){ printf( "atom[%i] Qs_=(%g,%g,%g,%g) Qs=(%g,%g,%g,%g) REQ.z=%g \n", i0+j,  Qs_.x,Qs_.y,Qs_.z,Qs_.w,  Qs.x,Qs.y,Qs.z,Qs.w,  REQ.z ); };

                REQ.x   += tipParams.x;
                REQ.yzw *= tipParams.yzw;

                //float3 shift=shift0;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){

                            fe += getLJQH   ( dp, REQ, 1.0f );
                            fe += getCoulomb( dp+(float3)(0.0f,0.0f,QZs.x) ,0.1f) * Qs_.x;
                            fe += getCoulomb( dp+(float3)(0.0f,0.0f,QZs.y) ,0.1f) * Qs_.y;
                            fe += getCoulomb( dp+(float3)(0.0f,0.0f,QZs.z) ,0.1f) * Qs_.z;
                            fe += getCoulomb( dp+(float3)(0.0f,0.0f,QZs.w) ,0.1f) * Qs_.w;

                            dp   +=lvec.a.xyz;
                            //shift+=lvec.a.xyz;
                        }
                        dp   +=shift_b;
                        //shift+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp   +=shift_c;
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // --- limit maximum force
    float renorm = 100.0/fabs(fe.w);
    if( renorm<1.f ){ fe*=renorm; }
    // --- OUTPUT

    write_imagef( imgOut, (int4){ia,ib,ic,0}, fe );

}

// ======================================
//    PPAFM_relaxStrokesTilted_convZ
// ======================================

#define OPT_FIRE 1
#define FIRE_FTDEC  0.5f
#define FIRE_FTINC  1.1f
#define FIRE_FDAMP  0.99f
#define F2SAFE      1e-8f

float3 tipForce( float3 dpos, float4 stiffness, float4 dpos0 ){
    float r = sqrt( dot( dpos,dpos) );
    return  (dpos-dpos0.xyz) * stiffness.xyz        // harmonic 3D
         + dpos * ( stiffness.w * (r-dpos0.w)/r );  // radial
}

float3 update_FIRE( float3 f, float3 v, float* dt, float* damp,    float dtmin, float dtmax, float damp0 ){
    // Bitzek, E., Koskinen, P., Gähler, F., Moseler, M., & Gumbsch, P. (2006). Structural Relaxation Made Simple. Physical Review Letters, 97(17), 170201.
    // https://doi.org/10.1103/PhysRevLett.97.170201
    // http://users.jyu.fi/~pekkosk/resources/pdf/FIRE.pdf
    float ff = dot(f,f);
    float vv = dot(v,v);
    float vf = dot(v,f);
    if( vf < 0 ){ // if velocity along direction of force
        v      *= 0;
        (*dt)   = fmax( dtmin, (*dt) * FIRE_FTDEC );
        (*damp) = damp0;
    }else{       // if velocity against direction of force
        // v = cV * v  + cF * F
        v       *= (1 - (*damp));
        v       +=  f * ( (*damp) * sqrt( vv / (ff + F2SAFE ) ) );
        (*dt)    = fmin( dtmax, (*dt) * FIRE_FTINC );
        (*damp) *= FIRE_FDAMP;
    }
    return v;
    //v  += f * dt;
    //p  += v * dt;
}

//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void PPAFM_scan(
    __read_only image3d_t  imgIn,   // 1
    __global  float4*      points,  // 2
    __global  float4*      FEs,     // 3
    __global  float4*      PPpos,   // 4
    const cl_Mat3  diGrid,          // 5
    const cl_Mat3  tipRot,          // 6
    float4 stiffness,               // 7
    float4 dpos0,                   // 8
    float4 relax_params,            // 9
    const int nz,                   // 10
    const int nMaxItr               // 11
){
    const float3 dTip   = tipRot.c.xyz * tipRot.c.w;
    float4 dpos0_=dpos0; dpos0_.xyz= rotMatT( dpos0_.xyz, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz );

    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0_.xyz;

    float dt      = relax_params.x;
    float damp    = relax_params.y;

    float dtmax = dt;
    float dtmin = dtmax*0.1f;
    float damp0 = damp;


    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    //const int iL = get_local_id   (0);
    //const int nL = get_local_size (0);


    //if(iG==0){printf("GPU:PPAFM_scan() nG=%i nz=%i relax_params(%g,%g,%g,%g) stiffness(%g,%g,%g,%g) dTip(%g,%g,%g) dpos0(%g,%g,%g)\n", nG, nz, relax_params.x,relax_params.y,relax_params.z,relax_params.w,  stiffness.x,stiffness.y,stiffness.z, stiffness.w, dTip.x,dTip.y,dTip.z, dpos0_.x,dpos0_.y,dpos0_.z );}
    //if(iG==0){printf("GPU:PPAFM_scan() ilvec{{%6.3f,%6.3f,%6.3f},{%6.3f,%6.3f,%6.3f},{%6.3f,%6.3f,%6.3f}}\n", diGrid.a.x,diGrid.a.y,diGrid.a.z, diGrid.b.x,diGrid.b.y,diGrid.b.z, diGrid.c.x,diGrid.c.y,diGrid.c.z  );}

    //if( (get_global_id(0)==0) ){ float4 fe = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );  printf( " pos (%g,%g,%g) feImg(%g,%g,%g,%g) \n", pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w );}
    //if( (get_global_id(0)==0) ){ printf( "dt %g damp %g \n", dt, damp ); }; return;

    int itr_tot = 0;

    for(int iz=0; iz<nz; iz++){
        float4 fe = float4Zero;
        float3 v  = float3Zero;
        int itr=0;
        for(itr=0; itr<nMaxItr; itr++){
        //for(int i=0; i<4; i++){
        //for(int i=0; i<1; i++){ // DEBUG
            //fe            = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            //fe            = interpFE_prec( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            const float4 coord = (float4)( dot(pos, diGrid.a.xyz),   dot(pos,diGrid.b.xyz), dot(pos,diGrid.c.xyz), 0.0f );
            // #if 0
                //coord +=(float4){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
                fe = read_imagef( imgIn, sampler_gff_norm, coord );
            // #else
                // fe = read_imagef_trilin_norm( imgIn, coord );
            //#endif

            float3 f      = fe.xyz * -1.0f;
            float3 dpos   = pos-tipPos;


            //float3 dpos_  = rotMat  ( dpos, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz );    // to tip-coordinates
            //float3 ftip   = tipForce( dpos_, stiffness, dpos0 );
            //f            += rotMatT ( ftip, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz );      // from tip-coordinates
            //f            += tipRot.c.xyz * surfFF.x;                                            // TODO: more sophisticated model of surface potential? Like Hamaker ?

            f +=  tipForce( dpos, stiffness, dpos0_ );  // Not rotated

            #if 1
                v = update_FIRE( f, v, &dt, &damp, dtmin, dtmax, damp0 );
                //if(get_global_id(0)==(64*128+64)){ printf( "itr,iz,i %i %i %i  |F| %g |v| %g <f,v> %g , (%g,%g,%g) (%g,%g,%g) damp %g dt %g \n", itr_tot, iz,i,  sqrt(dot(f,f)), sqrt(dot(v,v)),  dot(f,v),  fe.x,fe.y,fe.z, pos.x, pos.y, pos.z, damp, dt ); }
            #else
                v        *=    (1.f - damp);
                //if(get_global_id(0)==(64*128+64)){ printf( "itr,iz,i %i %i %i  |F| %g |v| %g <f,v> %g , (%g,%g,%g) (%g,%g,%g) damp %g dt %g \n", itr_tot, iz,i,  sqrt(dot(f,f)), sqrt(dot(v,v)),  dot(f,v),  fe.x,fe.y,fe.z, pos.x, pos.y, pos.z, damp, dt ); }
            #endif
            v        += f * dt;
            pos.xyz  += v * dt;
            itr++;
            if(dot(f,f)<relax_params.z) break;
        }
        //itr_tot+itr;

        // if(1){ // output tip-rotated force
        //     fe.xyz = rotMat( fe.xyz, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz);
        // }

        // DEBUG
        // //float3 pos    = tipPos.xyz;
        // const float4 coord   = (float4)( dot(pos, diGrid.a.xyz), dot(pos,diGrid.b.xyz), dot(pos,diGrid.c.xyz), 0.0f );
        // //const float4 coord   = (float4)( tipPos.x*0.1f, tipPos.y*0.1f, tipPos.z*0.1, 0.0f );
        // fe = read_imagef( imgIn, sampler_gff_norm, coord );
        // //fe = (float4){ tipPos.x, tipPos.y, tipPos.z, get_global_id(0)*nz+iz };

        //const float4 coord = (float4)( dot(pos, diGrid.a.xyz),   dot(pos,diGrid.b.xyz), dot(pos,diGrid.c.xyz), 0.0f );
        //fe = read_imagef( imgIn, sampler_gff_norm, coord );

        PPpos[ iz*nG + iG ] = (float4){pos,(float)itr}; // Store Result
        FEs  [ iz*nG + iG ] = fe; // Store Result
        //FEs[ iz*nG + iG ] = (float4){pos,fe.w}; // Store Result
        //FEs[ iz*nG + iG ] = (float4){pos,(float)itr}; // Store Result

        tipPos += dTip.xyz;
        pos    += dTip.xyz;
    }

}


//__attribute__((reqd_work_group_size(1,1,1)))
__kernel void PPAFM_scan_df(
    __read_only image3d_t  imgIn,   // 1
    __global  float4*      points,  // 2
    __constant  float*     weighs,  // 3
    __global  float4*      FEs,     // 4
    const cl_Mat3  diGrid,          // 5
    const cl_Mat3  tipRot,          // 6
    float4 stiffness,               // 7
    float4 dpos0,                   // 8
    float4 relax_params,            // 9
    float4 surfFF,                  // 11
    const int nz,                   // 12
    const int nzout                 // 13
){

    __local float  WEIGHTS[64];

    const float3 dTip   = tipRot.c.xyz * tipRot.c.w;
    float4 dpos0_=dpos0; dpos0_.xyz= rotMatT( dpos0_.xyz, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz );

    float3 tipPos = points[get_global_id(0)].xyz;
    float3 pos    = tipPos.xyz + dpos0_.xyz;

    float dt      = relax_params.x;
    float damp    = relax_params.y;

    float dtmax = dt;
    float dtmin = dtmax*0.1f;
    float damp0 = damp;

    //if( (get_global_id(0)==0) ){     float4 fe = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );  printf( " pos (%g,%g,%g) feImg(%g,%g,%g,%g) \n", pos.x, pos.y, pos.z, fe.x,fe.y,fe.z,fe.w );}
    //if( (get_global_id(0)==0) ){ printf( "dt %g damp %g \n", dt, damp ); }; return;

    const int ioff = get_global_id(0)*nzout;
    const int nzw   = nz-nzout;
    const int iL=get_local_id(0);
    const int nL=get_local_size(0);
    for (int i=iL; i<nzw; i+=nL ){
        WEIGHTS[i] = weighs[i];
    }
    for (int iz=0; iz<nzout; iz++ ){
        FEs[ioff+iz] = 0.0f;
    }
    barrier(CLK_LOCAL_MEM_FENCE);

    int itr_tot = 0;

    for(int iz=0; iz<nz; iz++){
        float4 fe;
        float3 v   = 0.0f;
        for(int i=0; i<128; i++){
        //for(int i=0; i<1; i++){ // DEBUG
            //fe            = interpFE( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            //fe            = interpFE_prec( pos, dinvA.xyz, dinvB.xyz, dinvC.xyz, imgIn );
            const float4 coord = (float4)( dot(pos, diGrid.a.xyz),   dot(pos,diGrid.b.xyz), dot(pos,diGrid.c.xyz), 0.0f );
            // #if 0
                //coord +=(float4){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
                const float4 fe_Paul = read_imagef( imgIn, sampler_gff_norm, coord );
            // #else
                // const float4 fe_Paul = read_imagef_trilin_norm( imgIn, coord );
            //#endif

            float3 f      = fe.xyz;
            float3 dpos   = pos-tipPos;
            float3 dpos_  = rotMat  ( dpos, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz );    // to tip-coordinates
            float3 ftip   = tipForce( dpos_, stiffness, dpos0 );

            f            += rotMatT ( ftip, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz );      // from tip-coordinates
            f            += tipRot.c.xyz * surfFF.x;                                            // TODO: more sophisticated model of surface potential? Like Hamaker ?

            //f      +=  tipForce( dpos, stiffness, dpos0_ );  // Not rotated

            #if 1
                v = update_FIRE( f, v, &dt, &damp, dtmin, dtmax, damp0 );
                //if(get_global_id(0)==(64*128+64)){ printf( "itr,iz,i %i %i %i  |F| %g |v| %g <f,v> %g , (%g,%g,%g) (%g,%g,%g) damp %g dt %g \n", itr_tot, iz,i,  sqrt(dot(f,f)), sqrt(dot(v,v)),  dot(f,v),  fe.x,fe.y,fe.z, pos.x, pos.y, pos.z, damp, dt ); }
            #else
                v        *=    (1 - damp);
                //if(get_global_id(0)==(64*128+64)){ printf( "itr,iz,i %i %i %i  |F| %g |v| %g <f,v> %g , (%g,%g,%g) (%g,%g,%g) damp %g dt %g \n", itr_tot, iz,i,  sqrt(dot(f,f)), sqrt(dot(v,v)),  dot(f,v),  fe.x,fe.y,fe.z, pos.x, pos.y, pos.z, damp, dt ); }
            #endif
            v        += f * dt;
            pos.xyz  += v * dt;

            itr_tot++;
            if(dot(f,f)<relax_params.z) break;
        }

        if(1){ // output tip-rotated force
            fe.xyz = rotMat( fe.xyz, tipRot.a.xyz, tipRot.b.xyz, tipRot.c.xyz);
        }

        // do the convolution
        for(int izout=0;izout<nzout;izout++){
            int jzw = iz - izout;
            if((jzw<nzw)&&(jzw>0)){
                FEs[ ioff + izout] += fe * WEIGHTS[jzw];
            }
        }
        //if( iz<nzout ) FEs[ioff+iz] = fe;
        tipPos += dTip.xyz;
        pos    += dTip.xyz;
    }

}


// ======================================================================
//                           add_DipoleField()
// ======================================================================

__attribute__((reqd_work_group_size(32,1,1)))
__kernel void addDipoleField(
    const int n,                     // 1
    __global float4*  ps,            // 2
    __global float4*  dipols,        // 3
    __write_only image3d_t  FE_Coul, // 4
    const int4     nGrid,            // 5
    const cl_Mat3  dGrid,            // 6
    const float4   grid_p0           // 7
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id  (0);
    const int nL = get_local_size(0);
    const int nab = nGrid.x*nGrid.y;
    const int ia  = iG%nGrid.x;
    const int ib  = (iG%nab)/nGrid.x;
    const int ic  = iG/nab;

    const int nMax = nab*nGrid.z;
    if(iG>nMax) return;

    //if(iG==0){printf("GPU::addDipoleField(nL=%i,nG=%i,nAtoms=%i,nPBC(%i,%i,%i))\n", nL, nG, n  );}

    float3 pos     = grid_p0.xyz + dGrid.a.xyz*ia + dGrid.b.xyz*ib  + dGrid.c.xyz*ic;
    float4 fe  = float4Zero;
    for (int i0=0; i0<n; i0+= nL ){
        int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = ps    [i];
        LCLJS [iL] = dipols[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<n ){
                float4 P     = LCLJS [j];
                float4 atom  = LATOMS[j];
                float3 d     = pos - atom.xyz;
                float  invr2 = 1.f / dot(d,d);
                float  invr  = sqrt(invr2);
                float  invr3 = invr*invr2;
                // https://en.wikipedia.org/wiki/Electric_dipole_moment#Potential_and_field_of_an_electric_dipole
                // Efield(R) = const *(    R*(Q/|R|^3) + R*3*<p|R>/|R|^5 - p/|R|^3

                float  VP  =  dot( P.xyz, d )*invr2;
                float4 fei = (float4){
                    (d*( P.w + 3*VP ) - P.xyz )*invr3,   // Force  (E-filed )
                       ( P.w +   VP           )*invr     // Energy (Potential)
                }*COULOMB_CONST;
                fe    += fei;

            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    int4 coord = (int4){ia,ib,ic,0};
    write_imagef( FE_Coul, coord, fe );
}


// =====================================================================
// =====================================================================
// =====================================================================
// =====================================================================
//
//                 LOCAL      SUPER       KERNELL
//
// =====================================================================
// =====================================================================
// =====================================================================
// =====================================================================

__attribute__((work_group_size_hint(64,1,1)))
__kernel void evalMMFFf4_local_test(
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,         // 2  [natoms]
    __global float4*  avel,         // 3
    //__global float4*  fapos,      //   [natoms] // Only Debug output
    //__global float4*  fneigh,     //   [nnode*4*2]
    __global float4*  constr,       // 4
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 6
    __global int4*    bkneighs,     // 7
    __global float4*  REQs,         // 8  [natoms] non-boding parametes {R0,E0,Q}
    __global float4*  apars,        // 9  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 10  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 11  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 12 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 13 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvecs,        // 14
    __global cl_Mat3* ilvecs,       // 15
    const int4        nPBC,         // 16
    const float4      GFFParams,    // 17
    const float4      MDpars,       // 18
    const int         niter         // 19
){

    const int iG  = get_global_id  (0); // intex of atom
    const int iS  = get_global_id  (1); // index of system
    const int nG  = get_global_size(0);
    const int nS  = get_global_size(1); // number of systems
    const int iL  = get_local_id   (0);
    const int iSL = get_local_id   (1); // number of systems
    const int nL  = get_local_size (0);
    const int nSL = get_local_size (1);

    const int natom = nDOFs.x;
    const int nnode = nDOFs.y;
    const int nvec  = natom+nnode;
    const int ncap  = natom - nnode;
    const int iLc   = nL+nnode;
    const int nL4   = nL*4;

    //if((iG==0)&&(iS==0))printf( "GPU::evalMMFFf4_local_test() G=%i/%i S=%i/%i  L=%i/%i sL=%i/%i natom=%i nnode=%i ncap=%i nvec=%i niter=%i \n", iG,nG, iS,nS, iL,nL, iSL,nSL,     natom, nnode, ncap, nvec, niter );

    const int i0a = iS*natom;
    const int i0n = iS*nnode;
    const int i0v = iS*nvec;
    // --- node atomm global indexes
    const int iaa = iG + i0a;
    const int ian = iG + i0n;
    const int iav = iG + i0v;
    // --- capping atom global indexes
    const int ica = iG + i0a-nnode;
    const int icn = iG + i0n-nnode;
    const int icv = iG + i0v-nnode;

    #define NNEIGH    4
    #define NATOM_LOC 128
    #define NNODE_LOC 64

    // ========= Local Memory

    __local float4  _apos    [NATOM_LOC  ];   // 2  [natoms] // atom position local
    __local float4  _fneigh  [NNODE_LOC*4];   // 4  [nnode*4*2]
    //__local float4  _fneigh  [NATOM_LOC];   // 4  [nnode*4*2]

    if(iG>=nnode){return; }

    // --- node
    float4 pa = apos[iav];
    float4 va = avel[iav];
    float3 fa = float3Zero;

    const float4  vbL  = bLs      [ian]; // bond lengths
    const float4  vbK  = bKs      [ian]; // bond stiffness
    const int4    ng   = neighs   [iaa]; // neighbors
    const int4    bk   = bkneighs [iaa]; // back-neighbors

    //if( (bk.x>=nL4)||(bk.y>=nL4)||(bk.z>=nL4)||(bk.w>=nL4) ){ printf("GPU::ERROR bk>nL4 iG,S[%i,%i] bk{%3i,%3i,%3i,%3i} nL4=%i \n", iG,iS, bk.x,bk.y,bk.z,bk.w, nL4 ); };

    // ---- array aliases
    const int*   ings  = (int*  )&ng;
    const float* bL    = (float*)&vbL;
    const float* bK    = (float*)&vbK;


    // ========== BIG LOOP  (Molecular Dynamics Loop)
    //for(int itr=0; itr<niter; itr++ ){

        float4  hs [4];              // direction vectors of bonds
        float3  fbs[4];              // force on neighbor sigma
        float3  f1,f2;               // force on pair of atoms (private temporary)

        barrier(CLK_LOCAL_MEM_FENCE);
        for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; }  // ---- Clean neighbor forces

        // ========= Bonds evaluation

        for(int i=0; i<NNEIGH; i++){
            float4 h;
            const int ing  = ings[i];
            const int ingv = ing+i0v;
            const int inga = ing+i0a;
            if(ing<0) break;
            if(ing>nnode)continue;
            h.xyz    = _apos[ing].xyz - pa.xyz;
            float  l = length(h.xyz);
            h.w      = 1./l;
            h.xyz   *= h.w;
            hs[i]    = h;
            if(iG<ing){
                evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;
            }
        }

        // ===========  Store locals ( Neighbor Forces )

        if(iG<NNODE_LOC){
            const int i4=iG*4;
            for(int i=0; i<NNEIGH; i++){  _fneigh[i4+i] = (float4){fbs[i],0.f}; }
            //for(int i=0; i<NNEIGH; i++){  _fneigh[i] = float4Zero; }
            //_fneigh[iG]=float4Zero;
        }

        // if(iS==0){
        //     printf( "[iG=%i]fbs1(%g,%g,%g) fbs2(%g,%g,%g) fbs3(%g,%g,%g) fbs4(%g,%g,%g) \n", iG, fbs[0].x,fbs[0].y,fbs[0].z,    fbs[1].x,fbs[1].y,fbs[1].z,    fbs[2].x,fbs[2].y,fbs[2].z,    fbs[3].x,fbs[3].y,fbs[3].z );
        //     printf( "[iG=%i]fng[%i](%g,%g,%g) fng[%i](%g,%g,%g) fng3[%i](%g,%g,%g) fng4[%i](%g,%g,%g) \n", iG, i4,_fneigh[i4].x,_fneigh[i4].y,_fneigh[i4].z,   i4+1, _fneigh[i4+1].x,_fneigh[i4+1].y,_fneigh[i4+1].z,    i4+2,_fneigh[i4+2].x,_fneigh[i4+2].y,_fneigh[i4+2].z,  i4+3,_fneigh[i4+3].x,_fneigh[i4+3].y,_fneigh[i4+3].z );
        // }

        // ===========  Assemble Forces

        barrier(CLK_LOCAL_MEM_FENCE);    // Make sure _fneigh is uptodate for all atoms

        // if((bk.x>=0)&&(bk.x<nL4)){ fa.xyz += _fneigh[bk.x].xyz; }
        // if((bk.y>=0)&&(bk.y<nL4)){ fa.xyz += _fneigh[bk.y].xyz; }
        // if((bk.z>=0)&&(bk.z<nL4)){ fa.xyz += _fneigh[bk.z].xyz; }
        // if((bk.w>=0)&&(bk.w<nL4)){ fa.xyz += _fneigh[bk.w].xyz; }

        if(bk.x>=0){ fa.xyz += _fneigh[bk.x].xyz; }
        if(bk.y>=0){ fa.xyz += _fneigh[bk.y].xyz; }
        if(bk.z>=0){ fa.xyz += _fneigh[bk.z].xyz; }
        if(bk.w>=0){ fa.xyz += _fneigh[bk.w].xyz; }

        //if(iS==0){  printf( "[iG=%i] bk{%3i,%3i,%3i,%3i}/%i fng(%g,%g,%g) \n", iG,  bk.x,bk.y,bk.z,bk.w, NNODE_LOC*4,  fng.x,fng.y,fng.z );}

        // ===========  Move Atoms
        //if(iS==0)printf( "[iG=%i] fa(%g,%g,%g)  bk{%3i,%3i,%3i,%3i} \n", iG, fa.x,fa.y,fa.z,   bk.x,bk.y,bk.z,bk.w );
        // --- move node atom
        va     *= MDpars.y;
        va.xyz += fa.xyz*MDpars.x;
        pa.xyz += va.xyz*MDpars.x;
        pa.w=0;va.w=0;

        //if(iS==0)printf( "[iG=%i] fa(%g,%g,%g) va(%g,%g,%g) pa(%g,%g,%g) \n", iG,  fa.x,fa.y,fa.z,  va.x,va.y,va.z,  pa.x,pa.y,pa.z );

        // ===========  Store locals ( Atomic Positions )

        _apos[iL] = pa;
        fa = float3Zero;

    //} // END OF THE BIG LOOP

    // ---- Store to global arrays
    apos[iav      ]=pa;
    avel[iav      ]=va;

};

__attribute__((work_group_size_hint(64,1,1)))
__kernel void evalMMFFf4_local2(
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,         // 2  [natoms]
    __global float4*  avel,         // 3
    //__global float4*  fapos,      //   [natoms] // Only Debug output
    //__global float4*  fneigh,     //   [nnode*4*2]
    __global float4*  constr,       // 4
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 6
    __global int4*    bkneighs,     // 7
    __global float4*  REQs,         // 8  [natoms] non-boding parametes {R0,E0,Q}
    __global float4*  apars,        // 9  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 10  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 11  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 12 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 13 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvecs,        // 14
    __global cl_Mat3* ilvecs,       // 15
    const int4        nPBC,         // 16
    const float4      GFFParams,    // 17
    const float4      MDpars,       // 18
    const int         niter         // 19
){

    const int iG = get_global_id  (0); // intex of atom
    const int iS = get_global_id  (1); // index of system
    const int nG = get_global_size(0);
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int natom = nDOFs.x;
    const int nnode = nDOFs.y;
    const int nvec  = natom+nnode;
    const int ncap  = natom-nnode;
    const int iLc   = iL+nnode;
    const int nfneigh = nnode*4;

    //if((iG==0)&&(iS==0))printf( "GPU::evalMMFFf4_local2() nL=%i nG=%i nS=%i natom=%i nnode=%i ncap=%i nvec=%i \n", nL,nG,nS, natom, nnode, ncap, nvec );
    //if((iG==0)&&(iS==0)){ if(NNODE_LOC<nnode)printf( "GPU::ERRROR NNODE_LOC(%i) < nnode(%i) \n", NNODE_LOC, nnode );}

    const int i0a = iS*natom;
    const int i0n = iS*nnode;
    const int i0v = iS*nvec;
    // --- node atomm global indexes
    const int iaa = iL + i0a;
    const int ian = iL + i0n;
    const int iav = iL + i0v;
    // --- capping atom global indexes
    const int ica = iL + i0a+nnode;
    const int icn = iL + i0n+nnode;
    const int icv = iL + i0v+nnode;
    const float R2damp = GFFParams.x*GFFParams.x;


    #define NNEIGH    4
    //#define NATOM_LOC 128
    //#define NNODE_LOC 64

    #define NATOM_LOC 64
    #define NNODE_LOC 32
    // ========= Local Memory

    // We want to fit only NODE atoms into local memory, the rest we do later

    __local float4  _apos    [NATOM_LOC  ];   // 2  [natoms] // atom position local
    __local float4  _ppos    [NATOM_LOC  ];   // 2  [natoms] // atom position local
    __local float4  _REQs    [NATOM_LOC  ];   // 6  [natoms] non-boding parametes {R0,E0,Q}
    __local float4  _fneigh  [NNODE_LOC*4];
    __local float4  _fneighpi[NNODE_LOC*4];

    // NOTE: We can avoid using  __local _fneigh[] if we write __private fbs[] and fps[] into __local aforce[] in synchronized manner in 4 steps (4 slots) using barrier(CLK_LOCAL_MEM_FENCE);

    const bool bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const bool bNode = iG<nnode;
    const bool bCap  = iG<ncap;

    // --- node
    float4  pa,va;                // position and velocity
    float4  REQa,cons;            // vdW parametes and constrains
    float4  par,vbL,vbK,vKs,vKp;  // force-field parameters
    float4  hp,vp;                // pi
    int4    ng,ngC,bk;            // neighbors
    if(bNode){
        va        = avel[iav];
        pa        = apos[iav];
        REQa      = REQs[iaa];
        hp        = apos[iav+natom];
        vp        = avel[iav+natom];
        _ppos[iL] = hp;
        _apos[iL] = pa;
        _REQs[iL] = REQa;
        // --- params
        ng   = neighs   [iaa];
        ngC  = neighCell[iaa];
        bk   = bkneighs [iaa];
        cons = constr   [iaa];
        par  = apars    [ian];    //     (xy=s0_ss,z=ssK,w=piC0 )
        vbL  = bLs      [ian];
        vbK  = bKs      [ian];
        vKs  = Ksp      [ian];
        vKp  = Kpp      [ian];
        //if( (bk.x>=nfneigh)||(bk.y>=nfneigh)||(bk.z>=nfneigh)||(bk.w>=nfneigh) ){ printf("GPU::ERROR node[%i|%i].bk{%3i,%3i,%3i,%3i}>nfneigh(%i)\n", iG,iS,  bk.x,bk.y,bk.z,bk.w, nfneigh ); };
    }

    // --- cap
    float4 vc,pc,REQc,c_cons;
    int4   c_ng,c_ngC,c_bk;
    if(bCap){
        vc         = avel     [icv]; // caping atom velocity
        pc         = apos     [icv]; // capping atom postion
        REQc       = REQs     [ica];
        c_ng       = neighs   [ica];
        c_ngC      = neighCell[ica];
        c_bk       = bkneighs [ica];
        c_cons     = constr   [ica];
        _apos[iLc] = pc;
        _REQs[iLc] = REQc;
        //if(iS==0)printf( "iL %i icv %i pc{%g,%g,%g}\n", iL, icv, pc.x,pc.y,pc.z );
        //if( (c_bk.x>=nfneigh)||(c_bk.y>=nfneigh)||(c_bk.z>=nfneigh)||(c_bk.w>=nfneigh) ){ printf("GPU::ERROR node[%i|%i].bk{%3i,%3i,%3i,%3i}>nfneigh(%i)\n", iG,iS,  c_bk.x,c_bk.y,c_bk.z,c_bk.w, nfneigh ); };
    }

    // ---- Parameters
    const cl_Mat3 lvec    = lvecs [iS];
    const cl_Mat3 invLvec = ilvecs[iS];
    const float3 shift0   = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;
    const float3 shift_a  = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b  = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);

    // =================================
    //            BIG LOOP
    // =================================

    for(int itr=0; itr<niter; itr++ ){   // BIG LOOP  (Molecular Dynamics Loop)
    //for(int itr=0; itr<1; itr++ ){   // BIG LOOP  (Molecular Dynamics Loop)

    float3  fa = float3Zero;    // force on atom position
    float3  fc = float3Zero;    // force on capping atom
    float3  fp = float3Zero;    // force on pi orbital

    barrier(CLK_LOCAL_MEM_FENCE);


    // ======================================
    //            COVALENT INTERACTIONS
    // ======================================
    if(bNode){

        const int*   ings  = (int*  )&ng;
        const float* bL    = (float*)&vbL;
        const float* bK    = (float*)&vbK;
        const float* Kspi  = (float*)&vKs;
        const float* Kppi  = (float*)&vKp;
        const float  ssC0  = par.x*par.x - par.y*par.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

        float4  hs [NNEIGH];         // direction vectors of bonds
        float3  fbs[NNEIGH];         // force on neighbor sigma
        float3  fps[NNEIGH];         // force on neighbor pi
        float3  f1,f2;               // force on pair of atoms (private temporary)

        for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }

        // ========= Bonds evaluation

        for(int i=0; i<NNEIGH; i++){
            float4 h;
            const int ing  = ings[i];
            const int ingv = ing+i0v;
            const int inga = ing+i0a;

            if(ing<0) break;
            h.xyz    = _apos[ing].xyz - pa.xyz;

            { // PBC bond vector correction
                float3 u  = (float3){ dot( invLvec.a.xyz, h.xyz ), dot( invLvec.b.xyz, h.xyz ), dot( invLvec.c.xyz, h.xyz ) };
                h.xyz   += lvec.a.xyz*(1.f-(int)(u.x+1.5f))
                        + lvec.b.xyz*(1.f-(int)(u.y+1.5f))
                        + lvec.c.xyz*(1.f-(int)(u.z+1.5f));
            }

            float  l = length(h.xyz);
            h.w      = 1./l;
            h.xyz   *= h.w;
            hs[i]    = h;

            //fbs[i]=(float3){(float)iG,(float)i,0.0f}; // DEBUG
            //evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fa+=f1;
            if(iG<ing){
                evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;

                // pi-pi
                float kpp = Kppi[i];
                if( (ing<nnode) && (kpp>1.e-6) ){
                    evalPiAling( hp.xyz, _ppos[ing].xyz, kpp,  &f1, &f2 );   fp+=f1;  fps[i]+=f2;
                }

            }

            // pi-sigma
            float ksp = Kspi[i];
            if(ksp>1.e-6){
                evalAngCos( (float4){hp.xyz,1.}, h, ksp, par.w, &f1, &f2 );   fp+=f1; fa-=f2;  fbs[i]+=f2;    //   pi-planarization (orthogonality)
            }

        }

        // ========= Angles evaluation

        for(int i=0; i<NNEIGH; i++){
            int ing = ings[i];
            if(ing<0) break;
            const float4 hi = hs[i];
            float4 REQi =_REQs[ing];   // ToDo: can be optimized
            //const int ingv  = ing+i0v;
            //const int inga  = ing+i0a;
            for(int j=i+1; j<NNEIGH; j++){
                int jng  = ings[j];
                if(jng<0) break;
                //const int jngv  = jng+i0v;
                //const int jnga  = jng+i0a;
                const float4 hj = hs[j];
                evalAngleCosHalf( hi, hj, par.xy, par.z, &f1, &f2 );
                fa    -= f1+f2;

                { // Remove vdW
                    float4 REQj=_REQs[jng];
                    float4 REQij;
                    REQij.x    = REQi.x   + REQj.x;
                    REQij.yzw  = REQi.yzw * REQj.yzw;
                    float3 dp  = (hj.xyz/hj.w) - (hi.xyz/hi.w);
                    float4 fij = getLJQH( dp, REQij, 1.0f );
                    f1 -=  fij.xyz;
                    f2 +=  fij.xyz;
                }

                fbs[i]+= f1;
                fbs[j]+= f2;
            }
        }

        // ========= Write neighbor forces (=recoils) to local memory

        const int i4=iG*4;
        for(int i=0; i<NNEIGH; i++){
            _fneigh  [i4+i] = (float4){fbs[i],0.f};
            _fneighpi[i4+i] = (float4){fps[i],0.f};
        }
        // if(iS==0){
        //     printf( "[iG=%i]fbs1(%g,%g,%g) fbs2(%g,%g,%g) fbs3(%g,%g,%g) fbs4(%g,%g,%g) \n", iG, fbs[0].x,fbs[0].y,fbs[0].z,    fbs[1].x,fbs[1].y,fbs[1].z,    fbs[2].x,fbs[2].y,fbs[2].z,    fbs[3].x,fbs[3].y,fbs[3].z );
        //     printf( "[iG=%i]fng[%i](%g,%g,%g) fng[%i](%g,%g,%g) fng3[%i](%g,%g,%g) fng4[%i](%g,%g,%g) \n", iG, i4,_fneigh[i4].x,_fneigh[i4].y,_fneigh[i4].z,   i4+1, _fneigh[i4+1].x,_fneigh[i4+1].y,_fneigh[i4+1].z,    i4+2,_fneigh[i4+2].x,_fneigh[i4+2].y,_fneigh[i4+2].z,  i4+3,_fneigh[i4+3].x,_fneigh[i4+3].y,_fneigh[i4+3].z );
        // }

    }  // if(bNode){


    // ======================================
    //          NON-COVALENT INTERACTIONS
    // ======================================

    // ========= Atom-to-Atom interaction ( N-body problem )

    for(int ja=0; ja<natom; ja++ ){
        const bool bBonded = ((ja==  ng.x)||(ja==  ng.y)||(ja==  ng.z)||(ja==  ng.w));
        const bool cBonded = ((ja==c_ng.x)||(ja==c_ng.y)||(ja==c_ng.z)||(ja==c_ng.w));
        const float4 aj    = _apos[ja];
        float4       REQ   = _REQs[ja];
        float4       cREQ  = REQ;

        REQ.x         += REQa.x;
        REQ.yzw       *= REQa.yzw;

        cREQ.x        += REQc.x;
        cREQ.yzw      *= REQc.yzw;

        float3 dp      = aj.xyz - pa.xyz;
        float3 dpc     = aj.xyz - pc.xyz;

        if(bPBC){

            dp     -= shift0;
            dpc    -= shift0;
            int ipbc=0;
            for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                    for(int ix=-nPBC.x; ix<=nPBC.x; ix++){

                        if(bNode && (ja!=iL) ){          // ---- Node Atom
                            if( !( bBonded && (
                                      ((ja==ng.x)&&(ipbc==ngC.x))
                                    ||((ja==ng.y)&&(ipbc==ngC.y))
                                    ||((ja==ng.z)&&(ipbc==ngC.z))
                                    ||((ja==ng.w)&&(ipbc==ngC.w))
                            ))){
                                fa += getLJQH( dp, REQ, R2damp ).xyz;
                            }
                        }

                        if(bCap && (ja!=iLc) ){       // ---- Capping atom
                            if(!(cBonded && (
                                      ((ja==c_ng.x)&&(ipbc==c_ngC.x))
                                    ||((ja==c_ng.y)&&(ipbc==c_ngC.y))
                                    ||((ja==c_ng.z)&&(ipbc==c_ngC.z))
                                    ||((ja==c_ng.w)&&(ipbc==c_ngC.w))
                            ))){
                                fc += getLJQH( dpc, cREQ, R2damp ).xyz;
                            }
                        }

                        ipbc++;
                        dp+=lvec.a.xyz;
                    }
                    dp+=shift_a;
                }
                dp+=shift_b;
            }

        }else{
            if(bNode && (ja!=iL ) && !bBonded ) fa += getLJQH( dp,   REQ, R2damp ).xyz;
            if(bCap  && (ja!=iLc) && !cBonded ) fc += getLJQH( dpc, cREQ, R2damp ).xyz;
            //if(bNode&& (ja!=iL )){ float3 fa_ = getLJQH( dp,   REQ, R2damp ).xyz; fa_ = limnitForce( fa_, 1.0f ); if(iS==0)printf( "nb[%i,%i] fa_(%g,%g,%g)\n", iL , ja, fa_.x,fa_.y,fa_.z ); }
            //if(bCap && (ja!=iLc)){ float3 fc_ = getLJQH( dpc, cREQ, R2damp ).xyz; fc_ = limnitForce( fc_, 1.0f ); if(iS==0)printf( "nb[%i,%i] fc_(%g,%g,%g)\n", iLc, ja, fc_.x,fc_.y,fc_.z ); }
            //if(bNode&&(ja!=iG )){ float3 fa_ = getLJQH( dp,   REQ, R2damp ).xyz; fa += limnitForce( fa_, 0.01f ); }
            //if(bCap &&(ja!=iLc)){ float3 fc_ = getLJQH( dpc, cREQ, R2damp ).xyz; fc += limnitForce( fc_, 0.01f ); }
        }

    } // for(ja)

    // ======================================
    //                MOVE
    // ======================================

    barrier(CLK_LOCAL_MEM_FENCE);    // Make sure _fneigh is uptodate for all atoms

    if(bNode){

        if(bk.x>=0){ fa += _fneigh[bk.x].xyz;  fp += _fneighpi[bk.x].xyz; }
        if(bk.y>=0){ fa += _fneigh[bk.y].xyz;  fp += _fneighpi[bk.y].xyz; }
        if(bk.z>=0){ fa += _fneigh[bk.z].xyz;  fp += _fneighpi[bk.z].xyz; }
        if(bk.w>=0){ fa += _fneigh[bk.w].xyz;  fp += _fneighpi[bk.w].xyz; }

        if( cons.w>0 ){  fa.xyz += (pa.xyz - cons.xyz)*-cons.w; }
        va     *= MDpars.y;
        va.xyz += fa.xyz*MDpars.x;
        pa.xyz += va.xyz*MDpars.x;
        pa.w=0;va.w=0;
        _apos[iL] = pa;

        // --- move pi    ...   This would simplify if we consider pi- is just normal cap atom (dummy atom)
        fp.xyz += hp.xyz * -dot( hp.xyz, fp.xyz );   // subtract forces  component which change pi-orbital lenght
        vp.xyz += hp.xyz * -dot( hp.xyz, vp.xyz );   // subtract veocity component which change pi-orbital lenght
        vp     *= MDpars.y;
        vp.xyz += fp.xyz*MDpars.x;
        hp.xyz += hp.xyz*MDpars.x;
        hp.xyz=normalize(hp.xyz);                    // normalize pi-orobitals
        hp.w=0;vp.w=0;
        _ppos[iL] = hp;

    }

    if(bCap){

        if(c_bk.x>=0){ fc += _fneigh[c_bk.x].xyz; }
        if(c_bk.y>=0){ fc += _fneigh[c_bk.y].xyz; }
        if(c_bk.z>=0){ fc += _fneigh[c_bk.z].xyz; }
        if(c_bk.w>=0){ fc += _fneigh[c_bk.w].xyz; }

        // --- move cap atom
        if( c_cons.w>0 ){fc.xyz += (pc.xyz - c_cons.xyz)*-c_cons.w; }
        vc     *= MDpars.y;
        vc.xyz += fc.xyz*MDpars.x;
        pc.xyz += vc.xyz*MDpars.x;
        pc.w=0;vc.w=0;

        _apos[iLc] = pc;

    }

    } // END OF THE BIG LOOP


    if(bNode){
        apos[iav      ]=pa;
        avel[iav      ]=va;
        apos[iav+natom]=hp;
        avel[iav+natom]=vp;
    }

    if(bCap){
        apos[ica]=pc;
        avel[ica]=vc;
    }


};

// =====================================================================
// =====================================================================
// =====================================================================
// =====================================================================


__attribute__((work_group_size_hint(32,1,1)))
__kernel void evalMMFFf4_local1(
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,         // 2  [natoms]
    __global float4*  avel,         // 3
    //__global float4*  fapos,      //   [natoms] // Only Debug output
    //__global float4*  fneigh,     //   [nnode*4*2]
    __global float4*  constr,       // 4
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 6
    __global int4*    bkneighs,     // 7
    __global float4*  REQs,         // 8  [natoms] non-boding parametes {R0,E0,Q}
    __global float4*  apars,        // 9  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 10  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 11  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 12 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 13 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvecs,        // 14
    __global cl_Mat3* ilvecs,       // 15
    const int4        nPBC,         // 16
    const float4      GFFParams,    // 17
    const float4      MDpars,       // 18
    const int         niter         // 19
){

    const int iG = get_global_id  (0); // intex of atom
    const int iS = get_global_id  (1); // index of system
    const int nG = get_global_size(0);
    const int nS = get_global_size(1); // number of systems
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int natom = nDOFs.x;
    const int nnode = nDOFs.y;
    const int nvec  = natom+nnode;
    const int ncap  = natom-nnode;
    const int iLc   = iL+nnode;
    const int nfneigh = nnode*4;

    //if((iG==0)&&(iS==0))printf( "GPU::evalMMFFf4_local1() nL=%i nG=%i nS=%i natom=%i nnode=%i ncap=%i nvec=%i \n", nL,nG,nS, natom, nnode, ncap, nvec );
    //if((iG==0)&&(iS==0)){ if(NNODE_LOC<nnode)printf( "GPU::ERRROR NNODE_LOC(%i) < nnode(%i) \n", NNODE_LOC, nnode );}

    const int i0a = iS*natom;
    const int i0n = iS*nnode;
    const int i0v = iS*nvec;
    // --- node atomm global indexes
    const int iaa = iL + i0a;
    const int ian = iL + i0n;
    const int iav = iL + i0v;
    // --- capping atom global indexes
    const int ica = iL + i0a+nnode;
    const int icn = iL + i0n+nnode;
    const int icv = iL + i0v+nnode;
    const float R2damp = GFFParams.x*GFFParams.x;


    #define NNEIGH    4
    //#define NATOM_LOC 128
    //#define NNODE_LOC 64

    #define NATOM_LOC 64
    #define NNODE_LOC 32
    // ========= Local Memory

    // We want to fit only NODE atoms into local memory, the rest we do later

    __local float4  _apos    [NATOM_LOC  ];   // 2  [natoms] // atom position local
    __local float4  _ppos    [NATOM_LOC  ];   // 2  [natoms] // atom position local
    __local float4  _REQs    [NATOM_LOC  ];   // 6  [natoms] non-boding parametes {R0,E0,Q}
    __local float4  _fneigh  [NNODE_LOC*4];
    __local float4  _fneighpi[NNODE_LOC*4];

    // NOTE: We can avoid using  __local _fneigh[] if we write __private fbs[] and fps[] into __local aforce[] in synchronized manner in 4 steps (4 slots) using barrier(CLK_LOCAL_MEM_FENCE);

    const bool bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const bool bNode = iG<nnode;

    float4       va   = avel     [iav];
    float4       pa   = apos     [iav];
    const float4 cons = constr   [iaa];
    const float4 REQa = REQs     [iaa];
    const int4 ng     = neighs   [iaa];
    const int4 ngC    = neighCell[iaa];
    const int4 bk     = bkneighs [iaa];
    _apos[iL] = pa;
    _REQs[iL] = REQa;
    float4  par,vbL,vbK,vKs,vKp;  // force-field parameters
    float4  hp,vp;                // pi
    if(bNode){
        hp        = apos[iav+natom];
        vp        = avel[iav+natom];
        _ppos[iL] = hp;
        // --- params

        par  = apars    [ian];    //     (xy=s0_ss,z=ssK,w=piC0 )
        vbL  = bLs      [ian];
        vbK  = bKs      [ian];
        vKs  = Ksp      [ian];
        vKp  = Kpp      [ian];
        //if( (bk.x>=nfneigh)||(bk.y>=nfneigh)||(bk.z>=nfneigh)||(bk.w>=nfneigh) ){ printf("GPU::ERROR node[%i|%i].bk{%3i,%3i,%3i,%3i}>nfneigh(%i)\n", iG,iS,  bk.x,bk.y,bk.z,bk.w, nfneigh ); };

        // ======= MMFF aliases
    }

    // ---- Parameters
    const cl_Mat3 lvec    = lvecs [iS];
    const cl_Mat3 invLvec = ilvecs[iS];
    const float3 shift0   = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;
    const float3 shift_a  = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b  = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);

    // =================================
    //            BIG LOOP
    // =================================

    for(int itr=0; itr<niter; itr++ ){   // BIG LOOP  (Molecular Dynamics Loop)
    //for(int itr=0; itr<1; itr++ ){   // BIG LOOP  (Molecular Dynamics Loop)

    float3  fa = float3Zero;    // force on atom position
    float3  fp = float3Zero;    // force on pi orbital

    barrier(CLK_LOCAL_MEM_FENCE);


    // ======================================
    //            COVALENT INTERACTIONS
    // ======================================
    if(bNode){

       // ---- Array aliases
        const int*   ings  = (int*  )&ng;
        const float* bL    = (float*)&vbL;
        const float* bK    = (float*)&vbK;
        const float* Kspi  = (float*)&vKs;
        const float* Kppi  = (float*)&vKp;
        const float  ssC0   = par.x*par.x - par.y*par.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

        float4  hs [NNEIGH];         // direction vectors of bonds
        float3  fbs[NNEIGH];         // force on neighbor sigma
        float3  fps[NNEIGH];         // force on neighbor pi
        float3  f1,f2;               // force on pair of atoms (private temporary)

        for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }

        // ========= Bonds evaluation

        for(int i=0; i<NNEIGH; i++){
            float4 h;
            const int ing  = ings[i];
            const int ingv = ing+i0v;
            const int inga = ing+i0a;

            if(ing<0) break;
            h.xyz    = _apos[ing].xyz - pa.xyz;

            if(bPBC){ // PBC bond vector correction
                float3 u  = (float3){ dot( invLvec.a.xyz, h.xyz ), dot( invLvec.b.xyz, h.xyz ), dot( invLvec.c.xyz, h.xyz ) };
                h.xyz   += lvec.a.xyz*(1.f-(int)(u.x+1.5f))
                        +  lvec.b.xyz*(1.f-(int)(u.y+1.5f))
                        +  lvec.c.xyz*(1.f-(int)(u.z+1.5f));
            }

            float  l = length(h.xyz);
            h.w      = 1./l;
            h.xyz   *= h.w;
            hs[i]    = h;

            //fbs[i]=(float3){(float)iG,(float)i,0.0f}; // DEBUG
            //evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fa+=f1;
            if(iL<ing){
                evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;

                // pi-pi
                float kpp = Kppi[i];
                if( (ing<nnode) && (kpp>1.e-6) ){
                    evalPiAling( hp.xyz, _ppos[ing].xyz, kpp,  &f1, &f2 );   fp+=f1;  fps[i]+=f2;
                }

            }

            // pi-sigma
            float ksp = Kspi[i];
            if(ksp>1.e-6){
                evalAngCos( (float4){hp.xyz,1.}, h, ksp, par.w, &f1, &f2 );   fp+=f1; fa-=f2;  fbs[i]+=f2;    //   pi-planarization (orthogonality)
            }

        }


        // ========= Angles evaluation
        for(int i=0; i<NNEIGH; i++){
            int ing = ings[i];
            if(ing<0) break;
            const float4 hi = hs[i];
            float4 REQi =_REQs[ing];   // ToDo: can be optimized
            //const int ingv  = ing+i0v;
            //const int inga  = ing+i0a;
            for(int j=i+1; j<NNEIGH; j++){
                int jng  = ings[j];
                if(jng<0) break;
                //const int jngv  = jng+i0v;
                //const int jnga  = jng+i0a;
                const float4 hj = hs[j];
                evalAngleCosHalf( hi, hj, par.xy, par.z, &f1, &f2 );
                fa    -= f1+f2;

                { // Remove vdW
                    float4 REQj=_REQs[jng];
                    float4 REQij;
                    REQij.x    = REQi.x   + REQj.x;
                    REQij.yzw  = REQi.yzw * REQj.yzw;
                    float3 dp  = (hj.xyz/hj.w) - (hi.xyz/hi.w);
                    float4 fij = getLJQH( dp, REQij, 1.0f );
                    f1 -=  fij.xyz;
                    f2 +=  fij.xyz;
                }

                fbs[i]+= f1;
                fbs[j]+= f2;
            }
        }

        // ========= Write neighbor forces (=recoils) to local memory
        const int i4=iL*4;
        for(int i=0; i<NNEIGH; i++){
            _fneigh  [i4+i] = (float4){fbs[i],0.f};
            _fneighpi[i4+i] = (float4){fps[i],0.f};
        }

    } // if(bNode){


    // ======================================
    //          NON-COVALENT INTERACTIONS
    // ======================================

    // ========= Atom-to-Atom interaction ( N-body problem )

    for(int ja=0; ja<natom; ja++ ){
        if( ja==iL ) continue;
        const bool bBonded = ((ja==  ng.x)||(ja==  ng.y)||(ja==  ng.z)||(ja==  ng.w));
        const float4 aj    = _apos[ja];
        float4       REQ   = _REQs[ja];
        REQ.x         += REQa.x;
        REQ.yzw       *= REQa.yzw;
        float3 dp      = aj.xyz - pa.xyz;

        if(bPBC){
            dp -= shift0;
            int ipbc=0;
            //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z );
            for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                    for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                        if( !( bBonded && (
                                    ((ja==ng.x)&&(ipbc==ngC.x))
                                ||((ja==ng.y)&&(ipbc==ngC.y))
                                ||((ja==ng.z)&&(ipbc==ngC.z))
                                ||((ja==ng.w)&&(ipbc==ngC.w))
                        ))){
                            fa += getLJQH( dp, REQ, R2damp ).xyz;
                        }

                        ipbc++;
                        dp+=lvec.a.xyz;
                    }
                    dp+=shift_a;
                }
                dp+=shift_b;
            }
        }else{
            if(!bBonded) fa += getLJQH( dp,   REQ, R2damp ).xyz;
        }

    } // for(ja)

    // ======================================
    //                MOVE
    // ======================================

    barrier(CLK_LOCAL_MEM_FENCE);    // Make sure _fneigh is uptodate for all atoms

    if(bk.x>=0){ fa += _fneigh[bk.x].xyz; fp += _fneighpi[bk.x].xyz; }
    if(bk.y>=0){ fa += _fneigh[bk.y].xyz; fp += _fneighpi[bk.y].xyz; }
    if(bk.z>=0){ fa += _fneigh[bk.z].xyz; fp += _fneighpi[bk.z].xyz; }
    if(bk.w>=0){ fa += _fneigh[bk.w].xyz; fp += _fneighpi[bk.w].xyz; }
    if( cons.w>0 ){  fa.xyz += (pa.xyz - cons.xyz)*-cons.w; }

    //if(iS==0)printf("iL[%i] pa(%g,%g,%g) fa(%g,%g,%g) bNode %i\n", iL,  pa.x,pa.y,pa.z,  fa.x,fa.y,fa.z, bNode );
    //if((iS==0)&&(iL==1))printf("iL[%i] bNode %i ng{%3i,%3i,%3i,%3i} fa{%g,%g,%g}\n", iL, bNode, ng.x,ng.y,ng.z,ng.w, fa.x,fa.y,fa.z );

    //pa.xyz += fa.xyz*0.01f;

    va     *= MDpars.y;
    va.xyz += fa.xyz*MDpars.x;
    pa.xyz += va.xyz*MDpars.x;
    pa.w=0;va.w=0;
    _apos[iL] = pa;


    if(bNode){
        // --- move pi    ...   This would simplify if we consider pi- is just normal cap atom (dummy atom)
        fp.xyz += hp.xyz * -dot( hp.xyz, fp.xyz );   // subtract forces  component which change pi-orbital lenght
        vp.xyz += hp.xyz * -dot( hp.xyz, vp.xyz );   // subtract veocity component which change pi-orbital lenght
        vp     *= MDpars.y;
        vp.xyz += fp.xyz*MDpars.x;
        hp.xyz += hp.xyz*MDpars.x;
        hp.xyz=normalize(hp.xyz);                    // normalize pi-orobitals
        hp.w=0;vp.w=0;
        _ppos[iL] = hp;
    }

    } // END OF THE BIG LOOP

    apos[iav]=pa;
    avel[iav]=va;
    if(bNode){
        apos[iav+natom]=hp;
        avel[iav+natom]=vp;
    }

};
