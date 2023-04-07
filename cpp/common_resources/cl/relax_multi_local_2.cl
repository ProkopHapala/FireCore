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
#define COULOMB_CONST   14.3996448915f  // [eV*Ang/e^2]

inline float2 udiv_cmplx( float2 a, float2 b ){ return (float2){  a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y }; }
//inline void     udiv_cmplx(               const VEC& b ){                            T x_ =    x*b.x +   y*b.y;         y =    y*b.x -   x*b.y;       x=x_;  }

inline double evalAngleCosHalf( const float4 hr1, const float4 hr2, const float2 cs0, double k, __private float3* f1, __private float3* f2 ){
    // This is much better angular function than evalAngleCos() with just a little higher computational cost ( 2x sqrt )
    float3 h  = hr1.xyz + hr2.xyz;
    float  c2 = dot(h,h)*0.25f;     // cos(a/2) = |ha+hb|
    float  s2 = 1.f-c2 + 1e-7;      // s2 must be positive
    float2 cso = (float2){ sqrt(c2), sqrt(s2) };
    float2 cs = udiv_cmplx( cs0, cso );
    float  E         =  k*( 1 - cs.x );  // just for debug ?
    float  fr        = -k*(     cs.y );
    c2 *= -2.f;
    fr /=  4.f*cso.x*cso.y;   //    |h - 2*c2*a| =  1/(2*s*c) = 1/sin(a)
    float  fr1    = fr*hr1.w;
    float  fr2    = fr*hr2.w;
    *f1 =  h*fr1  + hr1.xyz*(fr1*c2);  //fa = (h - 2*c2*a)*fr / ( la* |h - 2*c2*a| );
    *f2 =  h*fr2  + hr2.xyz*(fr2*c2);  //fb = (h - 2*c2*b)*fr / ( lb* |h - 2*c2*b| );
    return E;
}

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

inline float evalPiAling( const float3 h1, const float3 h2,  float K, __private float3* f1, __private float3* f2 ){  // interaction between two pi-bonds
    float  c = dot(h1,h2);
    float3 hf1,hf2;
    hf1 = h2 - h1*c;
    hf2 = h1 - h2*c;
    bool sign = c<0; if(sign) c=-c;
    float E    = -K*c;
    float fang =  K;
    if(sign)fang=-fang;
    hf1 *= fang;
    hf2 *= fang;
    *f1=hf1;
    *f2=hf2;
    return E;
}

inline float evalBond( float3 h, float dl, float k, __private float3* f ){
    float fr = dl*k;
    *f = h * fr;
    return fr*dl*0.5;
}

inline float4 getLJQ( float3 dp, float3 REQ, float R2damp ){
    // ---- Electrostatic
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 +  R2damp);
    float   Ec    =  COULOMB_CONST*REQ.z*sqrt( ir2_ );
    // --- LJ 
    float  ir2 = 1.f/r2;
    float  u2  = REQ.x*REQ.x*ir2;
    float  u6  = u2*u2*u2;
    float vdW  = u6*REQ.y;
    float E    =       (u6-2.f)*vdW     + Ec  ;
    float fr   = -12.f*(u6-1.f)*vdW*ir2 - Ec*ir2_;
    return  (float4){ dp*fr, E };
}

// ======================================================================
//                          getMMFFf4()
// ======================================================================

// In this kernell we assume that Node atoms and Capping atoms go in pairs => it is efficient to put them in pairs ( node & cap ) 

__kernel void getMMFFf4_local(
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    const int niter,
    // Dynamical
    __global float4*  apos,         // 2  [natoms]
    __global float4*  avel,
    //__global float4*  fapos,      // 3  [natoms] // Only Debug output    
    //__global float4*  fneigh,     // 4  [nnode*4*2]
    __global float4*  constr,       // 8
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,
    __global int4*    bkneighs,
    __global float4*  REQs,         // 6  [natoms] non-boding parametes {R0,E0,Q} 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    
    __global float4*  bLs,          // 8  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
    
    __global float16* bpars,        // 8  [nnode]  bond lengths  for each neighbor
    __global cl_Mat3* lvecs,        // 12
    __global cl_Mat3* ilvecs,       // 13
    const int4 nPBC,                // 8
    const float Rdamp,               // 9
    const float4      MDpars       // 1
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
    const int iLc = nL+nnode;

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

    // We want to fit only NODE atoms into local memory, the rest we do later

    __local float4  _apos    [NATOM_LOC  ];   // 2  [natoms] // atom position local
    __local float4  _REQs    [NATOM_LOC  ];   // 6  [natoms] non-boding parametes {R0,E0,Q}
    __local float4  _ppos    [NATOM_LOC  ];   // 2  [nnode]  // pi-postion local
    __local float4  _fneigh  [NNODE_LOC*4];   // 4  [nnode*4*2]
    __local float4  _fneighpi[NNODE_LOC*4];

    // NOTE: We can avoid using  __local _fneigh[] if we write __private fbs[] and fps[] into __local aforce[] in synchronized manner in 4 steps (4 slots) using barrier(CLK_LOCAL_MEM_FENCE);   

    const bool bNonde = iG<nnode;
    const bool bCap   = iG<(natoms-nnode);

    // --- node
    float4  pa   = apos[iav];
    float4  va   = avel[iav];
    float4  REQa = REQs[iaa];
    _apos[iL]    = pa;
    _REQs[iL]    = REQa;
    // --- cap
    float4  vc   = avel[icv]; // caping atom velocity
    float4  pc   = apos[icv]; // capping atom postion
    float4  REQc = REQs[ica];
    _apos[iLc]   = pc;
    _REQs[iLc]   = REQc;
    // --- pi
    float4  hp  = apos[iav+natom];
    float4  vp  = avel[iav+natom];
    if(bNonde){
        _ppos [iL]       = hp;
    }
    
    barrier(CLK_LOCAL_MEM_FENCE);   // Make sure globals are read into memory

    // ========= Private Memory

    // ---- Parameters
    const cl_Mat3 lvec    = lvecs [iS];
    const cl_Mat3 invLvec = ilvecs[iS];
    const float  R2damp = Rdamp*Rdamp;
    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);
 
    //  ---- Node
    const int4   ng   = neighs   [iaa];
    const int4   ngC  = neighCell[iaa];
    const int4   bk   = bkneighs [iaa];
    const float4 REQi = REQs     [iaa];
    const float4 cons = constr   [iaa];

    //  ---- cap
    const int4   c_ng   = neighs   [ica];
    const int4   c_ngC  = neighCell[ica];
    const int4   c_bk   = bkneighs [ica];
    const float4 c_REQi = REQs     [ica];
    const float4 c_cons = constr   [ica];

    const float4 par  = apars[ian];    //     (xy=s0_ss,z=ssK,w=piC0 )
    const float4 vbL  = bLs  [ian];       
    const float4 vbK  = bKs  [ian];       
    const float4 vKs  = Ksp  [ian];       
    const float4 vKp  = Kpp  [ian];

    // ====== This makes sence only for None-Atoms
    //  => Maybe it would be better to run one thread for 2 atoms: 
    //      (1) one Node atoms and 
    //      (2) one Capping atom (if would affect just non-bonded interactions)
    
    float3  fa = float3Zero;    // force on atom position
    float3  fc = float3Zero;    // force on capping atom 
    float3  fp = float3Zero;    // force on pi orbital

    float   E=0;
    float3  f1,f2;               // force on pair of atoms (private temporary)
    //float3  f1_,f2_;

    // ======= MMFF aliases
    // ---- Array aliases
    const int*   ings  = (int*  )&ng; 
    const float* bL    = (float*)&vbL; 
    const float* bK    = (float*)&vbK;
    const float* Kspi  = (float*)&vKs;  
    const float* Kppi  = (float*)&vKp; 
    const float  ssC0   = par.x*par.x - par.y*par.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
 
    // ---- Dynamical
    float4  hs [4];              // direction vectors of bonds
    float3  fbs[4];              // force on neighbor sigma
    float3  fps[4];              // force on neighbor pi

    // =================================
    //            BIG LOOP
    // =================================

    for(int itr=0; itr<niter; itr++ ){   // BIG LOOP  (Molecular Dynamics Loop) 
    barrier(CLK_LOCAL_MEM_FENCE);

    // ======================================
    //            COVALENT INTERACTIONS
    // ======================================

    // ========= Evaluate Bonds

    for(int i=0; i<NNEIGH; i++){ fbs[i]=float3Zero; fps[i]=float3Zero; }
    
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

        float epp = 0;
        float esp = 0;        
        if(iG<ing){
            E+= evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;                   
            float kpp = Kppi[i];
            // pi-pi
            if( (ing<nnode) && (kpp>1.e-6) ){
                epp += evalPiAling( hp.xyz, _ppos[ing].xyz, kpp,  &f1, &f2 );   fp+=f1;  fps[i]+=f2;
                E+=epp;
            }
        } 
        // pi-sigma 
        float ksp = Kspi[i];
        if(ksp>1.e-6){  
            esp += evalAngCos( (float4){hp.xyz,1.}, h, ksp, par.w, &f1, &f2 );   fp+=f1; fa-=f2;  fbs[i]+=f2;    //   pi-planarization (orthogonality)
            E+=epp;
        }
    }
    
    //  ============== Angles 
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
            E += evalAngleCosHalf( hi, hj, par.xy, par.z, &f1, &f2 );            
            fa    -= f1+f2;
            { // Remove vdW
                float4 REQj=_REQs[jng];
                float4 REQij;
                REQij.x    = REQi.x   + REQj.x;
                REQij.yzw  = REQi.yzw * REQj.yzw; 
                float3 dp  = (hj.xyz/hj.w) - (hi.xyz/hi.w); 
                float4 fij = getLJQ( dp, REQij.xyz, 1.0f );
                f1 -=  fij.xyz;
                f2 +=  fij.xyz;
            }
            fbs[i]+= f1;
            fbs[j]+= f2;
        }
    }

    // ========= Write fneigh to local memory

    const int i4=iG*4;
    for(int i=0; i<NNEIGH; i++){
        _fneigh  [i4+i] = (float4){fbs[i],0};
        _fneighpi[i4+i] = (float4){fps[i],0};
    }

    // ======================================
    //          NON-COVALENT INTERACTIONS
    // ======================================

    // ========= Atom-to-Atom interaction ( N-body problem )    
    for(int ja=0; ja<natom; ja+= nL ){
        if( ja==iG )continue;
        const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
        const bool cBonded = ((ja==c_ng.x)||(ja==c_ng.y)||(ja==c_ng.z)||(ja==c_ng.w));
        const float4 aj    = _apos[ja];
        const float4  REQ  = _REQs[ja];
        const float4 cREQ  = REQ;

        REQ.x         += REQ.x;
        REQ.yzw       *= REQ.yzw;
        float3 dp      = aj.xyz - shift0;
        float3 dpc     = dp - pc.xyz;
        dp            -=      pa.xyz;
 
        const float4 REQK = _REQs[ja];

        int ipbc=0; 
        for(int iy=0; iy<3; iy++){
            for(int ix=0; ix<3; ix++){
                
                // ---- Node Atom
                if(bBonded){
                    if(
                          ((ja==ng.x)&&(ipbc==ngC.x))
                        ||((ja==ng.y)&&(ipbc==ngC.y))
                        ||((ja==ng.z)&&(ipbc==ngC.z))
                        ||((ja==ng.w)&&(ipbc==ngC.w))
                    )continue; // skipp pbc0
                }
                float4 f = getLJQ( dp, REQK.xyz, R2damp );
                fa += f.xyz;

                // ---- Capping atom
                if(bBonded){
                    if(
                          ((ja==c_ng.x)&&(ipbc==c_ngC.x))
                        ||((ja==c_ng.y)&&(ipbc==c_ngC.y))
                        ||((ja==c_ng.z)&&(ipbc==c_ngC.z))
                        ||((ja==c_ng.w)&&(ipbc==c_ngC.w))
                    )continue; // skipp pbc0
                }
                f = getLJQ( dp, REQK.xyz, R2damp );
                fc += f.xyz;
                
                ipbc++; 
                dp+=lvec.a.xyz;
            }
            dp+=shift_a;
        }
    }

    // ======================================
    //                MOVE 
    // ======================================

    // ========= Assemble forces from neighbors
    barrier(CLK_LOCAL_MEM_FENCE);    // Make sure _fneigh is uptodate for all atoms
     
    if(bk.x>=0){ fa += _fneigh[bk.x].xyz; }
    if(bk.y>=0){ fa += _fneigh[bk.y].xyz; }
    if(bk.z>=0){ fa += _fneigh[bk.z].xyz; }
    if(bk.w>=0){ fa += _fneigh[bk.w].xyz; }

    if(bNonde){
        if(bk.x>=0){ fp += _fneighpi[bk.x].xyz; }
        if(bk.y>=0){ fp += _fneighpi[bk.y].xyz; }
        if(bk.z>=0){ fp += _fneighpi[bk.z].xyz; }
        if(bk.w>=0){ fp += _fneighpi[bk.w].xyz; }
    }

    // =============== DYNAMICS   ( Leap Frog)

    // --- move node atom
    if( cons.w>0 ){  fa.xyz += (pa.xyz - cons.xyz)*-cons.w; }
    va     *= MDpars.y;
    va.xyz += fa.xyz*MDpars.x;
    pa.xyz += va.xyz*MDpars.x;
    pa.w=0;va.w=0;
    _apos[iL] = pa;

    // --- move cap atom
    if( c_cons.w>0 ){fc.xyz += (pc.xyz - c_cons.xyz)*-c_cons.w; }
    va     *= MDpars.y;
    vc.xyz += fc.xyz*MDpars.x;
    pc.xyz += vc.xyz*MDpars.x;
    pc.w=0;vc.w=0;
    _apos[iLc] = pc;

    // --- move pi    ...   This would simplify if we consider pi- is just normal cap atom (dummy atom)
    fp.xyz += hp.xyz * -dot( hp.xyz, fp.xyz );   // subtract forces  component which change pi-orbital lenght
    vp.xyz += hp.xyz * -dot( hp.xyz, vp.xyz );   // subtract veocity component which change pi-orbital lenght
    vp     *= MDpars.y;
    vp.xyz += fp.xyz*MDpars.x;
    hp.xyz += hp.xyz*MDpars.x; 
    hp.xyz=normalize(hp.xyz);                    // normalize pi-orobitals
    hp.w=0;vp.w=0;
    _apos[iL+natom] = hp;

    // =============== FORCE DONE
    //aforce[iav] = fe;           // store force before limit
    //aforce[iav] = float4Zero;     // clean force

    } // END OF THE BIG LOOP

    apos[iav]=pa;
    avel[iav]=va;
    if(bNode){
        apos[iav+natom]=hp;
        avel[iav+natom]=vp;
    }

};
