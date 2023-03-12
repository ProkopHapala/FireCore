
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
#define COULOMB_CONST   14.399644f  // [eV*Ang/e^2]

float evalAngCos( const float4 hr1, const float4 hr2, float K, float c0, __private float3* f1, __private float3* f2 ){
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

float4 getLJQ( float3 dp, float3 REQ, float R2damp ){
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
    const int4 nDOFs,              // 1   (nAtoms,nnode)
    // Dynamical
    __global float4*  apos,        // 2    [natoms]
    __global float4*  fapos,       // 3    [natoms]     
    __global float4*  fneigh,      // 4    [nnode*4]
    //__global float4*  pipos,       //  [nnode]    - Stored in apos
    //__global float4*  fpipos,      //  [nnode]    - Stored in fapos
    //__global float4*  fneighpi,    //  [nnode*4]  - stored fneigh
    // parameters
    __global int4*    neighs,       // 5  [nnode]  neighboring atoms
    __global float4*  REQKs,        // 6  [natoms] non-boding parametes {R0,E0,Q} 
    __global float4*  apars,        // 7 [nnode]  per atom forcefield parametrs {c0ss,Kss,c0sp}
    __global float4*  bLs,          // 8 [nnode]  bond lengths  for each neighbor
    __global float4*  bKs,          // 9 [nnode]  bond stiffness for each neighbor
    __global float4*  Ksp,          // 10 [nnode]  stiffness of pi-alignment for each neighbor
    __global float4*  Kpp,          // 11 [nnode]  stiffness of pi-planarization for each neighbor
    __global cl_Mat3* lvec,         // 12
    __global cl_Mat3* invLvec       // 13
){

    const int iG = get_global_id (0);   // intex of atom
    const int iS = get_global_id (1);   // index of system
    //const int nG = get_global_size(0);
    //const int nS = get_global_size(1);  // number of systems
    //const int iL = get_local_id  (0);
    //const int nL = get_local_size(0);
    const int nAtoms=nDOFs.x;
    const int nnode =nDOFs.y;
    const int nvec  = nAtoms+nnode;

    const int iaa = iG + iS*nAtoms; 
    const int ian = iG + iS*nnode; 
    const int iav = iG + iS*nvec;

    #define NNEIGH 4

    /*
    if(ia==0){ printf( "GPU::getMMFFf4() nnode=%i nAtoms=%i size=%i \n", nnode, nAtoms, nG ); }
    if(ia==0)for(int i=0; i<nnode; i++){
        printf( "GPU[%i] ", i );
        printf( "ngs{%2i,%2i,%2i,%2i} ", neighs[i].x, neighs[i].y, neighs[i].z, neighs[i].w );
        printf( "apar{%6.3f,%6.3f,%6.3f,%6.3f} ", apars[i].x, apars[i].y, apars[i].z, apars[i].w );
        printf(  "BL{%6.3f,%6.3f,%6.3f,%6.3f} ", bLs[i].x, bLs[i].y, bLs[i].z, bLs[i].w );
        printf(  "BK{%6.3f,%6.3f,%6.3f,%6.3f} ", bKs[i].x, bKs[i].y, bKs[i].z, bKs[i].w );
        printf( "Ksp{%6.3f,%6.3f,%6.3f,%6.3f} ", Ksp[i].x, Ksp[i].y, Ksp[i].z, Ksp[i].w );
        printf( "Kpp{%6.3f,%6.3f,%6.3f,%6.3f} ", Kpp[i].x, Kpp[i].y, Kpp[i].z, Kpp[i].w );
        printf( "\n" );
    }
    */
    

    // ========= Private Memory

    // ---- Dynamical
    float4  hs [4];              // direction vectors of bonds
    //float4  fbs[4];            // force on neighbor sigma
    //float4  fps[4];            // force on neighbor pi
    float3  fbs[4];              // force on neighbor sigma
    float3  fps[4];              // force on neighbor pi
    float3  fa  = float3Zero;    // force on center position 
    float3  fpi = float3Zero;    // force on pi orbital
    float E=0;

    // ---- Params
    const int4   ng  = neighs[iaa];    
    const float3 pa  = apos[iav].xyz;
    const float3 hpi = apos[iav+nAtoms].xyz; 
    const float4 par = apars[ian];     
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

    // ========= Evaluate Bonds

    float3 f1,f2;         // temporary forces
    for(int i=0; i<NNEIGH; i++){
        float4 h;
        fbs[i]=float3Zero;
        fps[i]=float3Zero;
        //fbs[i]=(float3){1,2,3};
        //fps[i]=(float3){4,5,6};
        int ing = ings[i];
        if(ing<0) break;
        h.xyz    = apos[ing].xyz - pa;    //printf( "[%i|%i] ing=%i h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", ia,i,ing, h.x,h.y,h.z, apos[ing].x,apos[ing].y,apos[ing].z,  pa.x,pa.y,pa.z ); 
        float  l = length(h.xyz); 
        h.w      = 1./l;
        h.xyz   *= h.w;
        hs[i]    = h;


        float epp = 0;
        float esp = 0;

        //printf( "[%i|%i] l %g h(%g,%g,%g) \n", ia,i, l, h.x,h.y,h.z ); 
        if(iG<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            E+= evalBond( h.xyz, l-bL[i], bK[i], &f1 );  fbs[i]-=f1;  fa+=f1;   
            //if(ia==0)printf( "GPU bond[%i|%i] kpp=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
            
            float kpp = Kppi[i];
            if( (ing<nnode) && (kpp>1.e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                //E += evalPiAling( hpi, pipos[ing].xyz, kpp,  &f1, &f2 );   fpi+=f1;  fps[i]+=f2;        //   pi-alignment     (konjugation)
                epp += evalPiAling( hpi, apos[ing+nAtoms].xyz, kpp,  &f1, &f2 );   fpi+=f1;  fps[i]+=f2;    //   pi-alignment     (konjugation)
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
            esp += evalAngCos( (float4){hpi,1.}, h, ksp, par.z, &f1, &f2 );   fpi+=f1; fa-=f2;  fbs[i]+=f2;    //   pi-planarization (orthogonality)
            E+=epp;
        }
        
        //printf( "GPU[%i|%i] esp=%g epp=%g \n", esp, epp );
        
    }
    
    //  ============== Angles 
    for(int i=0; i<NNEIGH; i++){
        int ing = ings[i];
        if(ing<0) break;
        const float4 hi = hs[i];

        for(int j=i+1; j<NNEIGH; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const float4 hj = hs[j];
            //printf( "[%i|%i,%i] hi(%g,%g,%g) hj(%g,%g,%g)\n", ia, i,j, hi.x,hi.y,hi.z,   hj.x,hj.y,hj.z );
            E += evalAngCos( hi, hj, par.y, par.x, &f1, &f2 );     // angles between sigma bonds
            //if(ia==0)printf( "GPU:ang[%i|%i,%i] kss=%g c0=%g c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, par.y, par.x, dot(hi.xyz,hj.xyz),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            { // Remove vdW
                float4 REQKi=REQKs[ing];   // ToDo: can be optimized
                float4 REQKj=REQKs[jng];
                float4 REQKij;
                REQKij.x  = REQKi.x  + REQKj.x;
                REQKij.yz = REQKi.yz * REQKj.yz; 
                float3 dp = (hi.xyz/hi.w) - (hj.xyz/hj.w); 
                float4 fij = getLJQ( dp, REQKij.xyz, 1.0f );
                //float4 fij = getLJQ( apos[ing].xyz-apos[jng].xyz, REQKij.xyz, 1.0f );
                f1 -=  fij.xyz;
                f2 +=  fij.xyz;
            }
            fbs[i]+= f1;
            fbs[j]+= f2;
            fa    -= f1+f2;
            //if(ia==0)printf( "GPU:fa[%i](%g,%g,%g)\n", ia, fa.x,fa.y,fa.z  );
            // ToDo: subtract non-covalent interactions
        }
    }
    
    // ========= Save results

    const int i4 =iG*4;
    const int i4p=i4+nnode*4;
    for(int i=0; i<NNEIGH; i++){
        fneigh[i4 +i] = (float4){fbs[i],0};
        fneigh[i4p+i] = (float4){fps[i],0};
        //fneigh[i4 +i] = (float4){fbs[i],ia};
        //fneigh[i4p+i] = (float4){fps[i],ia};
        //fneighpi[i4+i] = (float4){fps[i],0};
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
    __global int4*    bkNeighs      // 7
){
    const int natoms=n.x;
    const int nnode =n.y;
    const int iG = get_global_id (0);
    const int nG = get_global_size(0);

    //if(iG==0)printf( "updateAtomsMMFFf4() natoms=%i nnode=%i natoms+nnode=%i size=%i dt=%g damp=%g Flimit=%g \n", natoms,nnode, natoms+nnode, nG, MDpars.x, MDpars.y, MDpars.z );
    /*
    if(iG==0){
    for(int i=0; i<natoms; i++){
        printf( "GPU[%i] ", i );
        //printf( "bkngs{%2i,%2i,%2i,%2i} ",         bkNeighs[i].x, bkNeighs[i].y, bkNeighs[i].z, bkNeighs[i].w );
        printf( "fapos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[i].x, aforce[i].y, aforce[i].z, aforce[i].w );
        //printf(  "avel{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i].x, avel[i].y, avel[i].z, avel[i].w );
        printf(  "apos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i].x, apos[i].y, apos[i].z, apos[i].w );
        printf( "\n" );
    }
    for(int i=0; i<nnode; i++){
        int i1=i+natoms;
        printf( "GPU[%i] ", i1 );
        printf(  "fpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", aforce[i1].x, aforce[i1].y, aforce[i1].z, aforce[i1].w );
        //printf(  "vpipos{%6.3f,%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z, avel[i1].w );
        printf(   "pipos{%6.3f,%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z, apos[i1].w );
        printf( "\n" );
    }
    for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
        int i1=i*4+j;
        int i2=i1+nnode*4;
        printf( "GPU[%i,%i] ", i, j );
        printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i1].x, fneigh[i1].y, fneigh[i1].z, fneigh[i1].w );
        printf( "fneighpi{%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i2].x, fneigh[i2].y, fneigh[i2].z, fneigh[i2].w );
        printf( "\n" );
    }}
    for(int i=0; i<nnode*4*2; i++){
        printf( "GPU[%i,%i,%i] ", i/4, i%4, i>=(nnode*4) );
        printf( "fneigh  {%6.3f,%6.3f,%6.3f,%6.3f} ", fneigh[i].x, fneigh[i].y, fneigh[i].z, fneigh[i].w );
        printf( "\n" );
    }
    }
    */

    if(iG>=(natoms+nnode)) return;
    float4 fe      = aforce[iG]; 
    const bool bPi = iG>=natoms;
    
    // ------ Gather Forces from back-neighbors
    int4 ngs;  
    if( bPi ){  // pis 
        const int ip0  = nnode*4;
        ngs            = bkNeighs[iG+natoms];
        ngs           += (int4){ip0,ip0,ip0,ip0};
    }else{     // atoms  
        ngs            = bkNeighs[iG];
    }
    if(ngs.x>=0){ fe += fneigh[ngs.x]; }
    if(ngs.y>=0){ fe += fneigh[ngs.y]; }
    if(ngs.z>=0){ fe += fneigh[ngs.z]; }
    if(ngs.w>=0){ fe += fneigh[ngs.w]; }

    aforce[iG] = fe; // store force before limit
    // ---- Limit Forces
    float fr2 = dot(fe.xyz,fe.xyz);
    if( fr2 > (MDpars.z*MDpars.z) ){
        fe.xyz*=(MDpars.z/sqrt(fr2));
    } 

    // ------ Move (Leap-Frog)
    float4 pe = apos[iG];
    float4 ve = avel[iG];
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
    avel[iG] = ve;
    apos[iG] = pe;
    
    // ------ Move Gradient-Descent
    /*
    float4 pe = apos[iG];
    //if(bPi){ fe.xyz += pe.xyz * -dot( pe.xyz, fe.xyz ); } // subtract forces  component which change pi-orbital lenght
    pe.xyz += fe.xyz*MDpars.x*0.1f;
    //if(bPi){ pe.xyz=normalize(pe.xyz); }
    pe.w=0;  // This seems to be needed, not sure why ?????
    apos[iG] = pe;
    */

    // ------ Store Force DEBUG
    //aforce[iG] = fe; // DEBUG - we do not have to save it, just to print it out on CPU
    //aforce[iG] = 0;  // ToDo:  this allows to ommit  updateAtomsMMFFf4() !!!! 
    //if(iG==0){ printf( "GPU::updateAtomsMMFFf4() END\n" ); }

}

__kernel void cleanForceMMFFf4(
    const int4        n,           // 2
    __global float4*  aforce,      // 5
    __global float4*  fneigh       // 6
){
    const int natoms=n.x;
    const int nnode =n.y;
    const int iG = get_global_id (0);
    if(iG>=(natoms+nnode)) return;
    aforce[iG]=float4Zero;
    //if(iG==0){ for(int i=0;i<(natoms+nnode);i++ ){printf("cleanForceMMFFf4[%i](%g,%g,%g)\n",i,aforce[i].x,aforce[i].y,aforce[i].z);} }
    if(iG<nnode){ 
        fneigh[iG*4+0]=float4Zero;
        fneigh[iG*4+1]=float4Zero;
        fneigh[iG*4+2]=float4Zero;
        fneigh[iG*4+3]=float4Zero;
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
    const int4 nPBC,                // 7
    const cl_Mat3 lvec,             // 8
    const float Rdamp               // 9
){
    __local float4 LATOMS[32];
    __local float4 LCLJS [32];
    const int iG = get_global_id  (0);
    const int nG = get_global_size(0);
    const int iL = get_local_id   (0);
    const int nL = get_local_size (0);

    const int natoms=ns.x;
    const int nnode =ns.y;

    //if(iG==0){ printf( "GPU::getNonBond() natoms %i nnode %i nG %i nL %i \n", natoms,nnode,nG,nL ); }

    if(iG>=natoms) return;
    //const bool   bNode = iG<nnode;   // All atoms need to have neighbors !!!!
    const bool   bPBC  = (nPBC.x+nPBC.y+nPBC.z)>0;
    const int4   ng    = neighs   [iG];
    const int4   ngC   = neighCell[iG];
    const float4 REQKi = REQKs [iG];
    const float3 posi  = atoms [iG].xyz;
    const float  R2damp = Rdamp*Rdamp;
    float4 fe          = float4Zero;

    //if(iG==0){ for(int i=0; i<natoms; i++)printf( "GPU[%i] ng(%i,%i,%i,%i) REQ(%g,%g,%g) \n", i, neighs[i].x,neighs[i].y,neighs[i].z,neighs[i].w, REQKs[i].x,REQKs[i].y,REQKs[i].z ); }

    // ========= Atom-to-Atom interaction ( N-body problem )    
    for (int i0=0; i0<natoms; i0+= nL ){
        const int i = i0 + iL;
        //if(i>=nAtoms) break;  // wrong !!!!
        LATOMS[iL] = atoms [i];
        LCLJS [iL] = REQKs [i];
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
                    float3 shifts=lvec.a.xyz*-nPBC.x;
                    for(int ia=-nPBC.x; ia<=nPBC.x; ia++){
                        shifts+=lvec.b.xyz*-nPBC.y;
                        for(int ib=-nPBC.y; ib<=nPBC.y; ib++){
                            shifts+=lvec.c.xyz*-nPBC.z;
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
                                //fe += getMorseQ( dp, REQK, R2damp );
                                fe += getLJQ( dp, REQK.xyz, R2damp );
                                ipbc++; 
                                shifts+=lvec.c.xyz;
                            }
                            shifts+=lvec.b.xyz;
                        }
                        shifts+=lvec.a.xyz;
                    }
                }else{
                    if(bBonded) continue;  // Bonded ?
                    // ToDo : what if bond is not within this cell ?????
                    //fe += getMorseQ( dp, REQK, R2damp );
                    //fe += getLJQ( dp, REQK.xyz, R2damp );
                    float4 fij = getLJQ( dp, REQK.xyz, R2damp );
                    fe += fij;
                    //if(iG==4){ printf( "GPU_LJQ[%i,%i|%i] fj(%g,%g,%g) R2damp %g REQ(%g,%g,%g) r %g \n", iG,ji,0, fij.x,fij.y,fij.z, R2damp, REQK.x,REQK.y,REQK.z, length(dp)  ); } 
                }
            }
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }
    
    forces[iG] = fe;
    //forces[iG] = fe*(-1.f);
    
}