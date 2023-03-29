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
    float  c2 = dot(h,h)*0.25f;              // cos(a/2) = |ha+hb|
    float  s2 = 1.f-c2;
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

__kernel void getMMFFf4(
    const int4 nDOFs,               // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,         // 2  [natoms]
    __global float4*  fapos,        // 3  [natoms]     
    __global float4*  fneigh,       // 4  [nnode*4*2]
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
    __global float4*  apars,        // 7  [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 8  [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 9  [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvecs,        // 12
    __global cl_Mat3* ilvecs        // 13
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

    const int i0a = iS*nAtoms; 
    const int i0n = iS*nnode; 
    const int i0v = iS*nvec;

    const int iaa = iG + i0a; 
    const int ian = iG + i0n; 
    const int iav = iG + i0v;

    #define NNEIGH 4

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
    const cl_Mat3 lvec    = lvecs [iS];
    const cl_Mat3 invLvec = ilvecs[iS];

    // ---- Dynamical
    float4  hs [4];              // direction vectors of bonds
    float3  fbs[4];              // force on neighbor sigma
    float3  fps[4];              // force on neighbor pi
    float3  fa  = float3Zero;    // force on center position 
    float3  fpi = float3Zero;    // force on pi orbital
    float E=0;

    // ---- Params
    const int4   ng  = neighs[iaa];    
    const float3 pa  = apos[iav].xyz;
    const float3 hpi = apos[iav+nAtoms].xyz; 
    const float4 par = apars[ian];    //     (xy=s0_ss,z=ssK,w=piC0 )
    const float4 vbL = bLs[ian];       
    const float4 vbK = bKs[ian];       
    const float4 vKs = Ksp[ian];       
    const float4 vKp = Kpp[ian];       

    // Temp Arrays
    const int*   ings  = (int*  )&ng; 
    const float* bL    = (float*)&vbL; 
    const float* bK    = (float*)&vbK;
    const float* Kspi  = (float*)&vKs;  
    const float* Kppi  = (float*)&vKp; 

    const float   ssC0   = par.x*par.x - par.y*par.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

    // const int iS_DBG = 5;
    // //const int iG_DBG = 0;
    // const int iG_DBG = 0;
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
            
            { // Remove vdW
                float4 REQi=REQKs[ing];   // ToDo: can be optimized
                float4 REQj=REQKs[jng];
                float4 REQij;
                REQij.x  = REQi.x  + REQj.x;
                REQij.yz = REQi.yz * REQj.yz; 
                float3 dp = (hj.xyz/hj.w) - (hi.xyz/hi.w); 
                float4 fij = getLJQ( dp, REQij.xyz, 1.0f );
                //float4 fij = getLJQ( apos[ing].xyz-apos[jng].xyz, REQKij.xyz, 1.0f );
                f1 -=  fij.xyz;
                f2 +=  fij.xyz;
                //f1.sub(fij);
                //f2.add(fij);
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
    fapos[iav       ] = (float4){fa ,0};
    fapos[iav+nAtoms] = (float4){fpi,0};
    //fpipos[ia] = (float4){fpi,0};
    

    //printf( "GPU[%i] fa(%g,%g,%g) fpi(%g,%g,%g)\n", ia, fa.x,fa.y,fa.z, fpi.x,fpi.y,fpi.z );

    //fapos[ia]=(float4){1,2,3,ia};

    //if(ia==0){ printf( "GPU::getMMFFf4() DONE\n" ); }
    
}

// ======================================================================
//                     updateAtomsMMFFf4()
// ======================================================================

__kernel void updateAtomsMMFFf4(
    const float4      MDpars,       // 1
    const int4        n,            // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr        // 8
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

    const int iS_DBG = 5;
    //const int iG_DBG = 0;
    const int iG_DBG = 1;

    //if((iG==iG_DBG)&&(iS==iS_DBG))printf( "updateAtomsMMFFf4() natoms=%i nnode=%i nvec=%i iS %i nG %i nS %i  iav_max %i  dt=%g damp=%g Flimit=%g \n", natoms,nnode, nvec, iS, nG, nS,  (nG-1)+(nS-1)*nvec, MDpars.x, MDpars.y, MDpars.z );
    
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

    // constrains
    float4 pe = apos[iav];
    if(iG<natoms){ 
       float4 cons = constr[ iaa ];
       if( cons.w>0 ){
            fe.xyz += (pe.xyz - cons.xyz)*-cons.w;
       }
    }

    // !!!!! Error code was "CL_INVALID_COMMAND_QUEUE" (-36)  is HERE !!!!!!
    aforce[iav] = fe; // store force before limit

    // ---- Limit Forces
    float fr2 = dot(fe.xyz,fe.xyz);
    if( fr2 > (MDpars.z*MDpars.z) ){
        fe.xyz*=(MDpars.z/sqrt(fr2));
    } 

    /*
    // ------ Move (kvazi-FIRE)    - proper FIRE need to reduce dot(f,v),|f|,|v| over whole system (3*N dimensions), this complicates paralell implementaion, therefore here we do it only over individual particles (3 dimensions)
    
    float4 ve = avel[iav];
    if(bPi){ 
        fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // subtract forces  component which change pi-orbital lenght
        ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // subtract veocity component which change pi-orbital lenght
    }
    ve.xyz += fe.xyz*MDpars.x;                       // according to LAMMPS implementation we should update the velocity by force first ... see https://github.com/lammps/lammps/blob/730e5d2e64106f3e5357fd739b44c7eec19c7d2a/src/min_fire.cpp#L393
    float ff = dot(fe.xyz,fe.xyz);
	float vv = dot(ve.xyz,ve.xyz);
    float vf = dot(ve.xyz,fe.xyz);
    #define ff_safety 1e-8
	if( vf < 0.0f ){
		ve.xyz=float3Zero;
	}else{
		float cf  =    MDpars.y * sqrt(vv/(ff+ff_safety));
		float cv  = 1.f - MDpars.y;
		ve.xyz    = cv * ve.xyz  + cf * fe.xyz;
	}
    pe.xyz += ve.xyz*MDpars.x;
    if(bPi){ 
        pe.xyz=normalize(pe.xyz);                   // normalize pi-orobitals
    }
    pe.w=0;ve.w=0;  // This seems to be needed, not sure why ?????
    avel[iav] = ve;
    apos[iav] = pe;
    */

    // // ------ Move (Leap-Frog)
    //float4 pe = apos[iav];
    float4 ve = avel[iav];
    if(bPi){ 
        fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz );   // subtract forces  component which change pi-orbital lenght
        ve.xyz += pe.xyz * -dot( pe.xyz, ve.xyz );   // subtract veocity component which change pi-orbital lenght
    }
    ve     *= MDpars.y;
    ve.xyz += fe.xyz*MDpars.x;
    pe.xyz += ve.xyz*MDpars.x;
    if(bPi){ 
        pe.xyz=normalize(pe.xyz);                   // normalize pi-orobitals
    }
    pe.w=0;ve.w=0;  // This seems to be needed, not sure why ?????
    avel[iav] = ve;
    apos[iav] = pe;
    
    // ------ Move Gradient-Descent
    // float4 pe = apos[iav];
    // //if(bPi){ fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz ); } // subtract forces  component which change pi-orbital lenght
    // pe.xyz += fe.xyz*MDpars.x*0.1f;
    // //if(bPi){ pe.xyz=normalize(pe.xyz); }
    // pe.w=0;  // This seems to be needed, not sure why ?????
    // apos[iav] = pe;
    
    // ------ Store Force DEBUG
    //aforce[iav] = fe; // DEBUG - we do not have to save it, just to print it out on CPU
    //aforce[iav] = 0;  // ToDo:  this allows to ommit  updateAtomsMMFFf4() !!!! 
    //if(iG==0){ printf( "GPU::updateAtomsMMFFf4() END\n" ); }
    
}



// ======================================================================
//                     updateAtomsMMFFf4()
// ======================================================================

__kernel void printOnGPU(
    const int4        n,            // 2
    __global float4*  apos,         // 3
    __global float4*  avel,         // 4
    __global float4*  aforce,       // 5
    __global float4*  fneigh,       // 6
    __global int4*    bkNeighs,     // 7
    __global float4*  constr        // 8
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

    const int iS_DBG = 5;
    //const int iG_DBG = 0;
    const int iG_DBG = 1;

    // if((iG==iG_DBG)&&(iS==iS_DBG))printf( "printOnGPU() natoms=%i nnode=%i nvec=%i iS %i nG %i nS %i  iav_max %i  dt=%g damp=%g Flimit=%g \n", natoms,nnode, nvec, iS, nG, nS,  (nG-1)+(nS-1)*nvec );
    // if((iG==iG_DBG)&&(iS==iS_DBG)){
    //     for(int i=0; i<natoms; i++){
    //         int ia=i + iS*nvec;
    //         printf( "GPU[%i=%i] ", i, ia );
    //         //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[ia].x, bkNeighs[ia].y, bkNeighs[ia].z, bkNeighs[ia].w );
    //         printf( "fapos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[ia].x, aforce[ia].y, aforce[ia].z, aforce[ia].w );
    //         //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[ia].x, avel[ia].y, avel[ia].z, avel[ia].w );
    //         printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[ia].x, apos[ia].y, apos[ia].z, apos[ia].w );
    //         printf(  "constr{%6.3f,%6.3f,%6.3f,%6.3f} ", constr[ia].x, constr[ia].y, constr[ia].z, constr[ia].w );
    //         printf( "\n" );
    //     }
    //     for(int i=0; i<nnode; i++){
    //         int i1=i+natoms + iS*nvec;
    //         printf( "GPU[%i=%i] ", i, i1 );
    //         printf(  "fpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[i1].x, aforce[i1].y, aforce[i1].z, aforce[i1].w );
    //         //printf(  "vpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z, avel[i1].w );
    //         printf(   "pipos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z, apos[i1].w );
    //         printf( "\n" );
    //     }
    //     for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
    //         int i1=i*4+j;
    //         int i2=i1+nnode*4;
    //         printf( "GPU[%i,%i] ", i, j );
    //         printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i1].x, fneigh[i1].y, fneigh[i1].z, fneigh[i1].w );
    //         printf( "fneighpi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i2].x, fneigh[i2].y, fneigh[i2].z, fneigh[i2].w );
    //         printf( "\n" );
    //     }}
    //     for(int i=0; i<nnode*4*2; i++){
    //         printf( "GPU[%i,%i,%i] ", i/4, i%4, i>=(nnode*4) );
    //         printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i].x, fneigh[i].y, fneigh[i].z, fneigh[i].w );
    //         printf( "\n" );
    //     }
    // }
}

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
    const int4 nPBC,                // 8
    const float Rdamp               // 9
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
    
    const int iS_DBG = 5;
    const int iG_DBG = 0;

    // if((iG==iG_DBG)&&(iS==iS_DBG)){ 
    //     printf( "GPU::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); 
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );        
    //     }
    // }

    if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const int4   ng    = neighs   [iaa];
    const int4   ngC   = neighCell[iaa];
    const float4 REQKi = REQKs    [iaa];
    const float3 posi  = atoms    [iav].xyz;
    const float  R2damp = Rdamp*Rdamp;
    float4 fe          = float4Zero;

    const cl_Mat3 lvec = lvecs[iS];

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    // ========= Atom-to-Atom interaction ( N-body problem )    
    for (int i0=0; i0<natoms; i0+= nL ){
        const int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms [i+i0v];
        LCLJS [iL] = REQKs [i+i0a];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            const int ji=j+i0;
            if( (ji!=iG) && (ji<natoms) ){   // ToDo: Should interact withhimself in PBC ?
                const float4 aj = LATOMS[j];
                const float3 dp = aj.xyz - posi;
                float4 REQK = LCLJS[j];
                REQK.x+=REQKi.x;
                REQK.yz*=REQKi.yz;

                const bool bBonded = ((ji==ng.x)||(ji==ng.y)||(ji==ng.z)||(ji==ng.w));

                if(bPBC){
                    // ToDo: it may be more effcient not construct pbc_shift on-the-fly
                    int ipbc=0;
                    float3 shift=lvec.a.xyz*-nPBC.x;
                    for(int ia=-nPBC.x; ia<=nPBC.x; ia++){
                        shift+=lvec.b.xyz*-nPBC.y;
                        for(int ib=-nPBC.y; ib<=nPBC.y; ib++){
                            shift+=lvec.c.xyz*-nPBC.z;
                            for(int ic=-nPBC.z; ic<=nPBC.z; ic++){     
                                if(bBonded){
                                    // Maybe We can avoid this by using Damped LJ or Buckingham potential which we can safely subtract in bond-evaluation ?
                                    if(
                                         ((ji==ng.x)&&(ipbc==ngC.x))
                                       ||((ji==ng.y)&&(ipbc==ngC.y))
                                       ||((ji==ng.z)&&(ipbc==ngC.z))
                                       ||((ji==ng.w)&&(ipbc==ngC.w))
                                    )continue; // skipp pbc0
                                }
                                //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                float4 fij = getLJQ( dp+shift, REQK.xyz, R2damp );
                                //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                fe += fij;
                                ipbc++; 
                                shift+=lvec.c.xyz;
                            }
                            shift+=lvec.b.xyz;
                        }
                        shift+=lvec.a.xyz;
                    }
                }else{
                    if(bBonded) continue;  // Bonded ?
                    // ToDo : what if bond is not within this cell ?????
                    //fe += getMorseQ( dp, REQK, R2damp );
                    //fe += getLJQ( dp, REQK.xyz, R2damp );
                    float4 fij = getLJQ( dp, REQK.xyz, R2damp );
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); } 
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    //forces[iav] = fe;
    forces[iav] += fe;
    //forces[iav] = fe*(-1.f);
    
}


// ======================================================================
// ======================================================================
//                           GridFF
// ======================================================================
// ======================================================================

__constant sampler_t sampler_1       =  CLK_NORMALIZED_COORDS_TRUE  | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;
__constant sampler_t sampler_2       =  CLK_NORMALIZED_COORDS_TRUE | CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_nearest =  CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
__constant sampler_t sampler_gff     =  CLK_NORMALIZED_COORDS_FALSE  | CLK_ADDRESS_REPEAT | CLK_FILTER_LINEAR;
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
//                           getNonBond()
// ======================================================================

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
    const float Rdamp,              // 9
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
    
    // const int iS_DBG = 5;
    // const int iG_DBG = 0;
    // if((iG==iG_DBG)&&(iS==iS_DBG)){ 
    //     printf( "GPU::getNonBond() natoms,nnode,nvec(%i,%i,%i) nS,nG,nL(%i,%i,%i) \n", natoms,nnode,nvec, nS,nG,nL ); 
    //     for(int i=0; i<nS*nG; i++){
    //         int ia = i%nS;
    //         int is = i/nS;
    //         if(ia==0){ cl_Mat3 lvec = lvecs[is];  printf( "GPU[%i] lvec(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f)(%6.3f,%6.3f,%6.3f) \n", is, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z  ); }
    //         //printf( "GPU[%i,%i] \n", is,ia,  );        
    //     }
    // }

    if(iG>=natoms) return;

    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const int4   ng    = neighs   [iaa];
    const int4   ngC   = neighCell[iaa];
    const float4 REQKi = REQKs    [iaa];
    const float3 posi  = atoms    [iav].xyz;
    const float  R2damp = Rdamp*Rdamp;
    float4 fe          = float4Zero;

    const cl_Mat3 lvec = lvecs[iS];

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    const float3 shift0  = lvec.a.xyz*nPBC.x + lvec.b.xyz*nPBC.y + lvec.c.xyz*nPBC.z;
    const float3 shift_a = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);
    const float3 shift_b = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);

    // ========= Atom-to-Atom interaction ( N-body problem )    
    for (int i0=0; i0<natoms; i0+= nL ){
        const int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms [i+i0v];
        LCLJS [iL] = REQKs [i+i0a];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            const int ji=j+i0;
            if( (ji!=iG) && (ji<natoms) ){   // ToDo: Should interact withhimself in PBC ?
                const float4 aj = LATOMS[j];
                float3 dp   = aj.xyz - posi;
                float4 REQK = LCLJS[j];
                REQK.x  +=REQKi.x;
                REQK.yz *=REQKi.yz;
                const bool bBonded = ((ji==ng.x)||(ji==ng.y)||(ji==ng.z)||(ji==ng.w));
                if(bPBC){
                    int ipbc=0;
                    //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                    dp -= shift0;
                    for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                        for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                            for(int ix=-nPBC.x; ix<=nPBC.x; ix++){
                                if(bBonded){
                                    // Maybe We can avoid this by using Damped LJ or Buckingham potential which we can safely subtract in bond-evaluation ?
                                    if(
                                         ((ji==ng.x)&&(ipbc==ngC.x))
                                       ||((ji==ng.y)&&(ipbc==ngC.y))
                                       ||((ji==ng.z)&&(ipbc==ngC.z))
                                       ||((ji==ng.w)&&(ipbc==ngC.w))
                                    )continue; // skipp pbc0
                                }
                                //fe += getMorseQ( dp+shifts, REQK, R2damp );
                                float4 fij = getLJQ( dp, REQK.xyz, R2damp );
                                //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,ipbc, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp+shift)  ); } 
                                fe += fij;
                                ipbc++; 
                                dp+=lvec.a.xyz;
                            }
                            dp+=shift_a;
                            //dp+=lvec.b.xyz;
                        }
                        dp+=shift_b;
                        //dp+=lvec.c.xyz;
                    }
                }else{
                    if(bBonded) continue;  // Bonded ?
                    // ToDo : what if bond is not within this cell ?????
                    //fe += getMorseQ( dp, REQK, R2damp );
                    //fe += getLJQ( dp, REQK.xyz, R2damp );
                    float4 fij = getLJQ( dp, REQK.xyz, R2damp );
                    fe += fij;
                    //if((iG==iG_DBG)&&(iS==iS_DBG)){  printf( "GPU_LJQ[%i,%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); } 
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // ========== Interaction with grid
    const float3 posg  = posi - grid_p0.xyz;
    const float4 coord = (float4)( dot(posg, diGrid.a.xyz),   dot(posg,diGrid.b.xyz), dot(posg,diGrid.c.xyz), 0.0f );
    #if 0
        coord +=(float3){0.5f,0.5f,0.5f,0.0f}; // shift 0.5 voxel when using native texture interpolation
        const float4 fe_Paul = read_imagef( FE_Paul, sampler_gff, coord );
        const float4 fe_Lond = read_imagef( FE_Lond, sampler_gff, coord );
        const float4 fe_Coul = read_imagef( FE_Coul, sampler_gff, coord );
    #else
        const float4 fe_Paul = read_imagef_trilin_( FE_Paul, coord );
        const float4 fe_Lond = read_imagef_trilin_( FE_Lond, coord );
        const float4 fe_Coul = read_imagef_trilin_( FE_Coul, coord );
    #endif
    //read_imagef_trilin( imgIn, coord );  // This is for higher accuracy (not using GPU hw texture interpolation)
    const float ej   = exp( REQKi.w * REQKi.x );
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
            const float4 fe_Paul = read_imagef_trilin_( FE_Paul, coord );
            const float4 fe_Lond = read_imagef_trilin_( FE_Lond, coord );
            const float4 fe_Coul = read_imagef_trilin_( FE_Coul, coord );
            const float4 fetot   = fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQKi.z;

            printf(  "%i %8.3f  %g %g %g %g   %g %g %g %g \n", ia0, pos.z,   fetot.w, fe_Paul.w*cP, fe_Lond.w*cL, fe_Coul.w*REQK.z,      fetot.z, fe_Paul.z*cP, fe_Lond.z*cL, fe_Coul.z*REQK.z  );

        }
    }
    */

    //forces[iav] = fe;
    forces[iav] += fe;
    //forces[iav] = fe*(-1.f);
    
}


// ======================================================================
//                           make_GridFF()
// ======================================================================

__kernel void make_GridFF(
    const int nAtoms,                // 1
    __global float4*  atoms,         // 2
    __global float4*  REQKs,         // 3
    __write_only image3d_t  FE_Paul, // 4
    __write_only image3d_t  FE_Lond, // 5
    __write_only image3d_t  FE_Coul, // 6
    const int4     nPBC,             // 7
    const int4     nGrid,            // 8
    const cl_Mat3  lvec,             // 9
    const float4   grid_p0           // 10
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


    float  R2damp = 1.0;
    /*
    if(iG==0){printf("GPU:make_GridFF(nL=%i,nG=%i,nAtoms=%i,nPBC(%i,%i,%i))\n", nL, nG, nAtoms, nPBC.x,nPBC.y,nPBC.z );}
    //if(iG==0){printf("GPU::make_GridFF(nAtoms=%i) \n", nAtoms );}
    if(iG==0){
        printf("GPU:make_GridFF(natoms=%i)\n", nAtoms);
        for(int i=0; i<nAtoms; i++){ printf("GPU:atom[%i] apos(%6.3f,%6.3f,%6.3f|%g) rekq(%6.3f,%10.7f,%6.3f|%g) \n", i, atoms[i].x,atoms[i].y,atoms[i].z,atoms[i].w,   REQKs[i].x,REQKs[i].y,REQKs[i].z,REQKs[i].w ); }
        int ia0 = 0;
        float4 REQK = REQKs[ia0];
        float Rii = REQK.x*2.f;
        float Eii = REQK.y*REQK.y;
        int np  = 20;
        float dx      = 0.1;
        float x0      = REQK.x*2.f - np*dx/2;
        
        //         expr = np.exp( b*Rj )
        // cP    = Ej*expr*expr
        // cL    = Ej*expr

        const float ej   = exp( REQK.w * REQK.x );
        const float cL   = ej*REQK.y;
        const float cP   = ej*cL;

        printf(  "GPU: REQK Ri %g Ei %g Q %g beta %g \n", REQK.x, REQK.y, REQK.z, REQK.w );
        printf(  "GPU: ej %g  cP %g cL %g \n", ej, cP, cL );
        printf(  "#i  r.x  E_LJ E_Paul E_Lond, E_Coul  Fx_LJ Fx_Paul Fx_Lond, Fx_Coul  \n" );
        for(int i=0; i<np; i++){
            float4 fe_Paul = float4Zero;
            float4 fe_Lond = float4Zero;
            float4 fe_Coul = float4Zero;
            
            float3 dp  = (float3){ x0 + dx*i, 0.f, 0.f  };
            float  r2  = dot(dp,dp);
            float  r   = sqrt(r2);
            float ir2  = 1/(r2+R2damp); 
            // ---- Electrostatic
            float   E  = COULOMB_CONST*REQK.z*sqrt(ir2);
            fe_Coul   += (float4)(dp*(E*ir2), E );
            // ---- Morse ( Pauli + Dispersion )
            float    e = exp( -REQK.w*(r-REQK.x) );
            float   eM = REQK.y*e;
            float   de = 2.f*REQK.w*eM/r;
            float4  fe = (float4)( dp*de, eM );
            fe_Paul += fe * e;
            fe_Lond += fe * (float4)( -1.0f,-1.0f,-1.0f, -2.0f );

            float ex  =      exp( -REQK.w * (r - Rii ) );
            float ELJ =      Eii*( ex*ex - 2.f*ex );
            float FLJ =  2.f*Eii*( ex*ex - ex     )*REQK.w;

            float4 fetot  = fe_Paul*cP  + fe_Lond*cL  +  fe_Coul*REQK.z*0.f;

            //printf(  "FE(RvdW[%i]) Paul(%g,%g,%g|%g) Lond(%g,%g,%g|%g) Coul(%g,%g,%g|%g)  \n", ia0, fe_Paul.x,fe_Paul.y,fe_Paul.z,fe_Paul.w,   fe_Lond.x,fe_Lond.y,fe_Lond.z,fe_Lond.w,    fe_Coul.x,fe_Coul.y,fe_Coul.z,fe_Coul.w  );
            //printf(  "%i %8.3f  %g %g    %g %g    %g %g  \n", ia0, dp.x, fe_Paul.x,fe_Paul.w,   fe_Lond.x,fe_Lond.w,    fe_Coul.x,fe_Coul.w  );
            printf(  "%i %8.3f  %g %g %g %g %g   %g %g %g %g %g  \n", ia0, dp.x,  ELJ, fetot.w, fe_Paul.w,fe_Lond.w,fe_Coul.w*REQK.z,   FLJ, fetot.x, fe_Paul.x,fe_Lond.x,fe_Coul.x*REQK.z  );

        }
    }
    */
    
    const int nMax = nab*nGrid.z;
    if(iG>nMax) return;

    const float3 dGrid_a = lvec.a.xyz*(1.f/(float)nGrid.x);
    const float3 dGrid_b = lvec.b.xyz*(1.f/(float)nGrid.y);
    const float3 dGrid_c = lvec.c.xyz*(1.f/(float)nGrid.z); 

    const float3 pos    = grid_p0.xyz  + dGrid_a.xyz*ia      + dGrid_b.xyz*ib      + dGrid_c.xyz*ic       // grid point within cell
                                       +  lvec.a.xyz*-nPBC.x + lvec .b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;  // most negative PBC-cell
    const float3 shift_b = lvec.b.xyz + lvec.a.xyz*(nPBC.x*-2.f-1.f);      //  shift in scan(iy)
    const float3 shift_c = lvec.c.xyz + lvec.b.xyz*(nPBC.y*-2.f-1.f);      //  shift in scan(iz) 

    float4 fe_Paul = float4Zero;
    float4 fe_Lond = float4Zero;
    float4 fe_Coul = float4Zero;
    for (int i0=0; i0<nAtoms; i0+= nL ){
        const int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms[i];
        LCLJS [iL] = REQKs[i];
        barrier(CLK_LOCAL_MEM_FENCE);
        for (int j=0; j<nL; j++){
            if( (j+i0)<nAtoms ){ 
                const float4 REQK   = LCLJS [j];
                float3 dp    = pos - LATOMS[j].xyz;
            
                //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc NONE dp(%g,%g,%g)\n", dp.x,dp.y,dp.z ); 
                //dp+=lvec.a.xyz*-nPBC.x + lvec.b.xyz*-nPBC.y + lvec.c.xyz*-nPBC.z;
                for(int iz=-nPBC.z; iz<=nPBC.z; iz++){
                    for(int iy=-nPBC.y; iy<=nPBC.y; iy++){
                        for(int ix=-nPBC.x; ix<=nPBC.x; ix++){

                            //if( (i0==0)&&(j==0)&&(iG==0) )printf( "pbc[%i,%i,%i] dp(%g,%g,%g)\n", ix,iy,iz, dp.x,dp.y,dp.z );   
                            float  r2  = dot(dp,dp);
                            float  r   = sqrt(r2);
                            float ir2  = 1/(r2+R2damp); 
                            // ---- Electrostatic
                            float   E  = COULOMB_CONST*REQK.z*sqrt(ir2);
                            fe_Coul   += (float4)(dp*(E*ir2), E );
                            // ---- Morse ( Pauli + Dispersion )
                            float    e = exp( -REQK.w*(r-REQK.x) );
                            float   eM = REQK.y*e;
                            float   de = 2.f*REQK.w*eM/r;
                            float4  fe = (float4)( dp*de, eM );
                            fe_Paul += fe * e;
                            fe_Lond += fe * (float4)( -1.0f,-1.0f,-1.0f, -2.0f );

                            //ipbc++; 
                            
                            dp+=lvec.a.xyz;
                        }
                        dp+=shift_b;
                        //dp+=lvec.a.xyz*(nPBC.x*-2.f-1.f);
                        //dp+=lvec.b.xyz;
                    }
                    dp+=shift_c;
                    //dp+=lvec.b.xyz*(nPBC.y*-2.f-1.f);
                    //dp+=lvec.c.xyz;
                }

            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    int4 coord = (int4){ia,ib,ic,0};
    write_imagef( FE_Paul, coord, fe_Paul );
    write_imagef( FE_Lond, coord, fe_Lond );
    write_imagef( FE_Coul, coord, fe_Coul );
}

// ======================================================================
//                           add_DipoleField()
// ======================================================================

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

    if(iG==0){printf("GPU::addDipoleField(nL=%i,nG=%i,nAtoms=%i,nPBC(%i,%i,%i))\n", nL, nG, n  );}

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