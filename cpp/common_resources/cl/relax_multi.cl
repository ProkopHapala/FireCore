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
#define COULOMB_CONST   14.3996448915f  // [eV*Ang/e^2]

inline float2 udiv_cmplx( float2 a, float2 b ){ return (float2){  a.x*b.x + a.y*b.y,  a.y*b.x - a.x*b.y }; }
//inline void     udiv_cmplx(               const VEC& b ){                            T x_ =    x*b.x +   y*b.y;         y =    y*b.x -   x*b.y;       x=x_;  }

inline float3 rotMat ( float3 v, float3 a, float3 b, float3 c ){ return (float3)(dot(v,a),dot(v,b),dot(v,c)); }
inline float3 rotMatT( float3 v, float3 a, float3 b, float3 c ){ return a*v.x + b*v.y + c*v.z; }

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

inline float4 getLJQH( float3 dp, float4 REQ, float R2damp ){
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

inline float4 getCoulomb( float3 dp, float R2damp ){
    // ---- Electrostatic
    float   r2    = dot(dp,dp);
    float   ir2_  = 1.f/(  r2 + R2damp);
    float   E    = COULOMB_CONST*sqrt( ir2_ );
    return  (float4){ dp*-E*ir2_, E };
}

float3 limnitForce( float3 f, float fmax ){
    float fr2 = dot(f,f);
    if( fr2>(fmax*fmax) ){ f*=(fmax/sqrt(fr2)); }
    return f;
}


// ======================================================================
// ======================================================================
//                           MMFF kernells
// ======================================================================
// ======================================================================


// ======================================================================
//                          getMMFFf4()
// ======================================================================

__attribute__((reqd_work_group_size(1,1,1)))
__kernel void getMMFFf4(
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,         // 2  [natoms]
    __global float4*  fapos,        // 3  [natoms]     
    __global float4*  fneigh,       // 4  [nnode*4*2]
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global int4*    neighCell,    // 5  [nnode]  neighboring atoms
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 8  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvecs,        // 12
    __global cl_Mat3* ilvecs,       // 13
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
        if(iG<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
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
            /*
            if(bSubtractVdW){ // Remove vdW
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
            */

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
    fapos[iav       ] = (float4){fa ,0}; // If we do  run it as first forcefield
    //fapos[iav       ] += (float4){fa ,0};  // If we not run it as first forcefield
    fapos[iav+nAtoms]  = (float4){fpi,0}; 


    //printf( "GPU[%i] fa(%g,%g,%g) fpi(%g,%g,%g)\n", ia, fa.x,fa.y,fa.z, fpi.x,fpi.y,fpi.z );

    //fapos[ia]=(float4){1,2,3,ia};

    //if(ia==0){ printf( "GPU::getMMFFf4() DONE\n" ); }
    
}

// ======================================================================
//                     updateAtomsMMFFf4()
// ======================================================================

/*
float2 KvaziFIREdamp( double c, float2 damp_lims, float2 clim ){
    float2 cvf;
    if      (c < clim.x ){   //-- force against veloctiy
        cvf.x = damp_lims.x; // v    // like 0.5 (strong damping)
        cvf.y = 0;           // f
    }else if(c > clim.y ){   //-- force alingned to velocity
        cvf.x = 1-damping;   // v    // like 0.99 (weak dampong damping)
        cvf.y =   damping;   // f
    }else{                   // -- force ~ perpendicular to velocity
        double f = (c-clim.x )/( clim.y - clim.x  );
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

float2 KvaziFIREdamp( double c, float2 clim, float2 damps ){
    float2 cvf;
    if      (c < clim.x ){   //-- force against veloctiy
        cvf.x = damps.x;     // v    // like 0.5 (strong damping)
        cvf.y = 0;           // f
    }else if(c > clim.y ){   //-- force alingned to velocity
        cvf.x = damps.y;     // v    // like 0.99 (weak dampong damping)
        cvf.y = 0;           // f
    }else{                   // -- force ~ perpendicular to velocity
        double t = (c-clim.x )/( clim.y - clim.x );
        cvf.x = damps.x + (damps.y-damps.x)*t;
        cvf.y = damps.y*t*(1.f-t)*2.f;
    }
    return cvf;
}

__attribute__((reqd_work_group_size(1,1,1)))
__kernel void updateAtomsMMFFf4(
    //const float4      MDpars,       // 1
    const int4        n,            // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr,       // 8
    __global float4*  MDparams     // 8
){
    const int natoms=n.x;
    const int nnode =n.y;
    const int nvec  = natoms+nnode;
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);

    //const int ian = iG + iS*nnode; 
    const int iaa = iG + iS*natoms; 
    const int iav = iG + iS*nvec;

    const float4 MDpars = MDparams[iS];

    //if((iS==0)&&(iG==0)){ printf("MDpars[%i] (%g,%g,%g,%g) \n", iS, MDpars.x,MDpars.y,MDpars.z,MDpars.w);  }

    const int iS_DBG = 5;
    //const int iG_DBG = 0;
    const int iG_DBG = 1;

    //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "updateAtomsMMFFf4() natoms=%i nnode=%i nvec=%i nG %i iS %i/%i  dt=%g damp=%g Flimit=%g \n", natoms,nnode, nvec, iS, nG, nS, MDpars.x, MDpars.y, MDpars.z );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     int i0a = iS*natoms;
    //     for(int i=0; i<natoms; i++){
    //         printf( "GPU:constr[%i](%7.3f,%7.3f,%7.3f |K= %7.3f) \n", i, constr[i0a+i].x,constr[i0a+i].y,constr[i0a+i].z,  constr[i0a+i].w   );
    //     }
    // }
    
    if(iG>=(natoms+nnode)) return;

    //aforce[iav] = float4Zero;

    float4 fe      = aforce[iav]; 
    const bool bPi = iG>=natoms;
    
    // ------ Gather Forces from back-neighbors

    int4 ngs = bkNeighs[ iav ];

    //if(iS==5)printf( "iG,iS %i %i ngs %i,%i,%i,%i \n", iG, iS, ngs.x,ngs.y,ngs.z,ngs.w );

    // WARRNING : bkNeighs must be properly shifted on CPU !!!!!!!!!!!!!!
    if(ngs.x>=0){ fe += fneigh[ngs.x]; }
    if(ngs.y>=0){ fe += fneigh[ngs.y]; }
    if(ngs.z>=0){ fe += fneigh[ngs.z]; }
    if(ngs.w>=0){ fe += fneigh[ngs.w]; }

    // ---- Limit Forces
    float Flimit = 10.0;
    float fr2 = dot(fe.xyz,fe.xyz);
    if( fr2 > (Flimit*Flimit) ){
        fe.xyz*=(Flimit/sqrt(fr2));
    }
    
    // =============== FORCE DONE
    aforce[iav] = fe;           // store force before limit
    //aforce[iav] = float4Zero;   // clean force   : This can be done in the first forcefield run (best is NBFF)
    
    // =============== DYNAMICS

    float4 ve = avel[iav];
    float4 pe = apos[iav];


    // -------constrains
    if(iG<natoms){ 
       float4 cons = constr[ iaa ];
       if( cons.w>0 ){
            fe.xyz += (pe.xyz - cons.xyz)*-cons.w;
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
    avel[iav] = ve;
    apos[iav] = pe;
    */
    
    
    // ------ Move (Leap-Frog)
    if(bPi){ 
        fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // subtract forces  component which change pi-orbital lenght
        ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // subtract veocity component which change pi-orbital lenght
    }
    // ve     *= MDpars.y;
    // ve.xyz += fe.xyz*MDpars.x;
    // pe.xyz += ve.xyz*MDpars.x;

    ve     *= 0.99f;
    ve.xyz += fe.xyz*0.01f;
    pe.xyz += ve.xyz*0.01f;
    if(bPi){ 
        pe.xyz=normalize(pe.xyz);                   // normalize pi-orobitals
    }
    pe.w=0;ve.w=0;  // This seems to be needed, not sure why ?????
    avel[iav] = ve;
    apos[iav] = pe;
    

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

__attribute__((reqd_work_group_size(1,1,1)))
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

__attribute__((reqd_work_group_size(1,1,1)))
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
//                           getNonBond()
// ======================================================================

__attribute__((reqd_work_group_size(32,1,1)))
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
    //__local float4 LATOMS[64];
    //__local float4 LCLJS [64];
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int iL = get_local_id   (0);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);
    const int nL = get_local_size (0);

    const int natoms=ns.x;
    const int nnode =ns.y;
    //const int nAtomCeil =ns.w;
    const int nvec  =natoms+nnode;

    //const int i0n = iS*nnode; 
    const int i0a = iS*natoms; 
    const int i0v = iS*nvec;
    //const int ian = iG + i0n;
    const int iaa = iG + i0a;
    const int iav = iG + i0v;
    
    const int iS_DBG = 0;
    const int iG_DBG = 0;

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
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    //const bool bPBC=false;

    const int4   ng    = neighs   [iaa];
    const int4   ngC   = neighCell[iaa];
    const float4 REQKi = REQKs    [iaa];
    const float3 posi  = atoms    [iav].xyz;
    const float  R2damp = GFFParams.x*GFFParams.x;
    float4 fe          = float4Zero;

    const cl_Mat3 lvec = lvecs[iS];

    //if(iG==0){ printf("GPU[iS=%i] lvec{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f}{%6.3f,%6.3f,%6.3f} \n", iS, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );  }

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);
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

    // ========= Atom-to-Atom interaction ( N-body problem )  

    //barrier(CLK_LOCAL_MEM_FENCE);
    for (int j0=0; j0<nG; j0+=nL){
        const int i=j0+iL;
        if(i<natoms){
            LATOMS[iL] = atoms [i+i0v];
            LCLJS [iL] = REQKs [i+i0a];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=j0+jl;
            if( (ja!=iG) && (ja<natoms) ){   // ToDo: Should interact withhimself in PBC ?
                const float4 aj = LATOMS[jl];
                float4 REQK     = LCLJS [jl];
                float3 dp       = aj.xyz - posi;
                //if((iG==44)&&(iS==0))printf( "[i=%i,ja=%i/%i,j0=%i,jl=%i/%i][iG/nG/na %i/%i/%i] aj(%g,%g,%g,%g) REQ(%g,%g,%g,%g)\n", i,ja,nG,j0,jl,nL,   iG,nG,natoms,   aj.x,aj.y,aj.z,aj.w,  REQK.x,REQK.y,REQK.z,REQK.w  );
                REQK.x  +=REQKi.x;
                REQK.yz *=REQKi.yz;
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
                //if( (j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                //if( (ji==1)&&(iG==0) )printf( "2 non-bond[%i,%i] bBonded %i\n",iG,ji,bBonded );

                /*
                if(bPBC){    
                    int ipbc=0;
                    dp += shift0;
                    // Fixed PBC size
                    for(int iy=0; iy<3; iy++){
                        for(int ix=0; ix<3; ix++){
                            //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "GPU ipbc %i(%i,%i) shift(%7.3g,%7.3g,%7.3g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z ); 
                            // Without these IF conditions if(bBonded) time of evaluation reduced from 61 [ms] to 51[ms]
                            if( !( bBonded &&(
                                    ((ja==ng.x)&&(ipbc==ngC.x))
                                    ||((ja==ng.y)&&(ipbc==ngC.y))
                                    ||((ja==ng.z)&&(ipbc==ngC.z))
                                    ||((ja==ng.w)&&(ipbc==ngC.w))
                            ))){
                                //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                //if( (ji==1)&&(iG==0)&&(iS==0) )printf( "ipbc %i(%i,%i) shift(%g,%g,%g)\n", ipbc,ix,iy, shift.x,shift.y,shift.z ); 
                                float4 fij = getLJQH( dp, REQK, R2damp );
                                //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                fe += fij;
                            }
                            ipbc++; 
                            dp    += lvec.a.xyz; 
                        }
                        dp    += shift_a;
                    }
                }else 
                */
                
                if( !bBonded ){
                    fe += getLJQH( dp, REQK, R2damp );
                }
            }
        }
        //barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    if(iG<natoms){
        forces[iav] = fe;           // If we do  run it as first forcefield 
        //forces[iav] += fe;        // If we not run it as first forcefield
        //forces[iav] = fe*(-1.f);
    }
}



// ======================================================================
// ======================================================================
//                           GridFF
// ======================================================================
// ======================================================================

// NOTE: https://registry.khronos.org/OpenCL/sdk/1.1/docs/man/xhtml/sampler_t.html
// CLK_ADDRESS_REPEAT - out-of-range image coordinates are wrapped to the valid range. This address mode can only be used with normalized coordinates. If normalized coordinates are not used, this addressing mode may generate image coordinates that are undefined.

__constant sampler_t samp0 =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

__constant sampler_t sampler_1       =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_2       =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_nearest =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_nearest_norm =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_gff      =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_gff_norm =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_gff   =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
//__constant sampler_t sampler_gff   =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

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

float2 spline_hermite( float x, float4 ys ){
    //float4  float y0, float y1, float dy0, float dy1 
    float dy0 = (ys.z-ys.x)*0.5;
    float dy1 = (ys.w-ys.y)*0.5;
    float y01 = ys.x-ys.y;
    float p2  = (-3.f*y01 -2.f*dy0 - dy1)*x;
    float p3  = ( 2.f*y01 +    dy0 + dy1)*x*x;
    return (float2){   ys.x + x*(dy0 +   p2 +   p3 ),  // value
	                  dy0 + 2.f*p2 + 3.f*p3        };        // derivative
}


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




float4 interpFE( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef( imgIn, sampler_1, coord );
    //return coord;
}

float4 interpFE_prec( float3 pos, float3 dinvA, float3 dinvB, float3 dinvC, __read_only image3d_t imgIn ){
    const float4 coord = (float4)( dot(pos,dinvA),dot(pos,dinvB),dot(pos,dinvC), 0.0f );
    return read_imagef_trilin( imgIn, coord ); 
    // read_imagef( imgIn, sampler_1, coord );
    //return coord;
}

// ======================================================================
//                           getNonBond_GridFF()
// ======================================================================

__attribute__((reqd_work_group_size(32,1,1)))
__kernel void getNonBond_GridFF(
    const int4 ns,                  // 1
    // Dynamical
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    // Parameters
    __global float4*  REQKs,        // 4
    __global int4*    neighs,       // 5
    __global int4*    neighCell,    // 6
    __global cl_Mat3* lvecs,        // 7
    const int4 nPBC,                // 8
    const float4  GFFParams,        // 9
    // GridFF
    __read_only image3d_t  FE_Paul, // 10
    __read_only image3d_t  FE_Lond, // 11
    __read_only image3d_t  FE_Coul, // 12
    const cl_Mat3  diGrid,          // 13
    const float4   grid_p0          // 14
    
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id  (0);
    const int iS = get_global_id  (1);
    const int iL = get_local_id   (0);
    const int nG = get_global_size(0);
    const int nS = get_global_size(1);
    const int nL = get_local_size (0);

    const int natoms=ns.x;
    const int nnode =ns.y;
    const int nvec  =natoms+nnode;

    //const int i0n = iS*nnode; 
    const int i0a = iS*natoms; 
    const int i0v = iS*nvec;
    //const int ian = iG + i0n;
    const int iaa = iG + i0a;
    const int iav = iG + i0v;
    
    const cl_Mat3 lvec = lvecs[iS];

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
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const int4   ng    = neighs   [iaa];
    const int4   ngC   = neighCell[iaa];
    const float4 REQKi = REQKs    [iaa];
    const float3 posi  = atoms    [iav].xyz;
    const float  R2damp = GFFParams.x*GFFParams.x;
    const float  alphaMorse = GFFParams.y;
    float4 fe           = float4Zero;

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);

    // ========= Atom-to-Atom interaction ( N-body problem )    
    for (int j0=0; j0<natoms; j0+= nL ){
        const int i = j0 + iL;
        LATOMS[iL] = atoms [i+i0v];
        LCLJS [iL] = REQKs [i+i0a];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int jl=0; jl<nL; jl++){
            const int ja=jl+j0;
            if( (ja!=iG) && (ja<natoms) ){   // ToDo: Should interact withhimself in PBC ?
                const float4 aj   = LATOMS[jl];
                float4       REQK = LCLJS [jl];
                float3 dp   = aj.xyz - posi;
                REQK.x  +=REQKi.x;
                REQK.yz *=REQKi.yz;
                const bool bBonded = ((ja==ng.x)||(ja==ng.y)||(ja==ng.z)||(ja==ng.w));
                if(bPBC){
                    int ipbc=0;
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                    dp -= shift0;
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if( !( bBonded &&(
                                          ((ja==ng.x)&&(ipbc==ngC.x))
                                        ||((ja==ng.y)&&(ipbc==ngC.y))
                                        ||((ja==ng.z)&&(ipbc==ngC.z))
                                        ||((ja==ng.w)&&(ipbc==ngC.w))
                                ))){
                                    //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                    float4 fij = getLJQH( dp, REQK, R2damp );
                                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                    fe += fij;
                                }
                                ipbc++; 
                                dp+=lvec.a.xyz;
                            }
                            dp+=shift_a;
                        }
                        dp+=shift_b;
                    }
                }else{
                    if(bBonded) continue;  // Bonded ?
                    float4 fij = getLJQH( dp, REQK, R2damp );
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); } 
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    if(iG>=natoms) return;

    // ========== Interaction with grid
    const float3 posg  = posi - grid_p0.xyz;
    const float4 coord = (float4)( dot(posg, diGrid.a.xyz),   dot(posg,diGrid.b.xyz), dot(posg,diGrid.c.xyz), 0.0f );
    // #if 0
        //coord +=(float4){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
        const float4 fe_Paul = read_imagef( FE_Paul, sampler_gff_norm, coord );
        const float4 fe_Lond = read_imagef( FE_Lond, sampler_gff_norm, coord );
        const float4 fe_Coul = read_imagef( FE_Coul, sampler_gff_norm, coord );
    // #else
        // const float4 fe_Paul = read_imagef_trilin_norm( FE_Paul, coord );
        // const float4 fe_Lond = read_imagef_trilin_norm( FE_Lond, coord );
        // const float4 fe_Coul = read_imagef_trilin_norm( FE_Coul, coord );
    //#endif
    //read_imagef_trilin( imgIn, coord );  // This is for higher accuracy (not using GPU hw texture interpolation)
    const float ej   = exp( alphaMorse* REQKi.x );
    const float cL   = ej*REQKi.y;
    const float cP   = ej*cL;
    fe  += fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQKi.z;
    
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

    forces[iav] = fe;        // If we do  run it as first forcefield
    //forces[iav] += fe;     // If we not run it as first forcefield
    //forces[iav] = fe*(-1.f);
    
}

// ======================================================================
//                           sampleGridFF()
// ======================================================================

__attribute__((reqd_work_group_size(1,1,1)))
__kernel void sampleGridFF(
    const int4 ns,                  // 1
    __global float4*  atoms,        // 2
    __global float4*  forces,       // 3
    __global float4*  REQs,        // 4
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

__attribute__((reqd_work_group_size(1,1,1)))
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


__attribute__((reqd_work_group_size(1,1,1)))
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