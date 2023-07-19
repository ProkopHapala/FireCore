
#ifndef MMFFsp3_loc_h
#define MMFFsp3_loc_h

#define ANG_HALF_COS  1

#include <omp.h>

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Forces.h"
#include "quaternion.h"
#include "SMat3.h"
#include "molecular_utils.h"

#include "NBFF.h"

//#include <cstdio>

inline double cos_damp_lin( double c, double& cv, double D, double cmin, double cmax  ){
    double cf;
    if      (c < cmin){
        cv = 0.;
        cf = 0.;
    }else if(c > cmax){
        cv = 1-D;
        cf =   D;
    }else{  // cos(v,f) from [ cmin .. cmax ]
        double u = (c-cmin)/(cmax-cmin);
        cv = (1.-D)*u;
        cf =     D *u;
    }
    return cf;
}

// ======================
// ====   MMFFsp3
// ======================

//class MMFFsp3_loc: public NBFF { public:

class MMFFsp3_loc : public NBFF { public:
    static constexpr const int nneigh_max = 4;
    // int natoms=0;         // [natoms] // from Atoms
    //Vec3d *   apos  =0;    // [natoms] // from Atoms
    //Vec3d *  fapos  =0;    // [natoms] // from NBFF
    //Quat4i*  neighs =0;    // [natoms] // from NBFF
    //Quat4i*  neighCell=0; // [natoms] // from NBFF
    //Quat4d*  REQs =0;       // [nnode]  // from NBFF
    //bool bPBC=false;       // from NBFF
    //Vec3i   nPBC;          // from NBFF 
    //Mat3d   lvec;          // from NBFF
    //double  Rdamp  = 1.0;  // from NBFF

    int  nDOFs=0,nnode=0,ncap=0,nvecs=0,ntors=0;
    double Etot,Eb,Ea, Eps,EppT,EppI;

    double *  DOFs = 0;   // degrees of freedom
    double * fDOFs = 0;   // forces

    bool doBonds  =true;
    bool doNeighs =true;
    bool doPiPiI  =true;
    bool doPiPiT  =true;
    bool doPiSigma=true;
    bool doAngles =true;
    bool doEpi    =true; 

    bool bEachAngle = false;
    bool bTorsion   = false;
    
    //                           c0     Kss    Ksp    c0_e
    Quat4d default_NeighParams{ -1.0,   1.0,   1.0,   -1.0 };

    // Dynamical Varaibles;
    //Vec3d *   apos=0;   // [natom]
    //Vec3d *  fapos=0;   // [natom]
    Vec3d *  pipos=0;   // [nnode]
    Vec3d * fpipos=0;   // [nnode]

    // Aux Dynamil
    Vec3d * fneigh  =0;  // [nnode*4]     temporary store of forces on atoms form neighbors (before assembling step)
    Vec3d * fneighpi=0;  // [nnode*4]     temporary store of forces on pi    form neighbors (before assembling step)

    // Params
    int   *  atypes  =0;

    //Quat4i*  neighs =0;   // [natoms] // from NBFF
    Quat4i*  bkneighs=0;   // [natoms]  inverse neighbors
    
    Quat4d*  apars=0;  // [nnode] per atom forcefield parametrs
    Quat4d*  bLs  =0;  // [nnode] bond lengths
    Quat4d*  bKs  =0;  // [nnode] bond stiffness
    Quat4d*  Ksp  =0;  // [nnode] stiffness of pi-alignment
    Quat4d*  Kpp  =0;  // [nnode] stiffness of pi-planarization

    Vec3d*   angles=0;

    Quat4i*  tors2atom =0;
    Quat4d*  torsParams=0;

    Quat4d*  constr=0;
    Vec3d * vapos = 0;

    Mat3d   invLvec;

    bool    bAngleCosHalf         = true;
    bool    bSubtractAngleNonBond = false;

    //int itr_DBG=0;

// =========================== Functions

void realloc( int nnode_, int ncap_, int ntors_=0 ){
    nnode=nnode_; ncap=ncap_; ntors=ntors_;
    natoms= nnode + ncap; 
    nvecs = natoms+nnode;  // each atom as also pi-orientiation (like up-vector)
    nDOFs = nvecs*3;
    //printf( "MMFFsp3::realloc() natom(%i,nnode=%i,ncap=%i), npi=%i, nbond=%i \n", natoms, nnode, ncap, npi, nbonds );
    int ipi0=natoms;
    
    _realloc0(  DOFs    , nDOFs , (double)NAN );
    _realloc0( fDOFs    , nDOFs , (double)NAN );
    apos   = (Vec3d*) DOFs ;
    fapos  = (Vec3d*)fDOFs;
    pipos  = apos  + ipi0;
    fpipos = fapos + ipi0;
    // ---- Aux
    _realloc0( fneigh  , nnode*4, Vec3dNAN );
    _realloc0( fneighpi, nnode*4, Vec3dNAN );
    // ----- Params [natom]
    _realloc0( atypes    , natoms, -1 );
    _realloc0( neighs    , natoms, Quat4iMinusOnes );
    _realloc0( neighCell , natoms, Quat4iMinusOnes );
    _realloc0( bkneighs  , natoms, Quat4iMinusOnes);
    _realloc0( apars     , nnode, Quat4dNAN );
    _realloc0( bLs       , nnode, Quat4dNAN );
    _realloc0( bKs       , nnode, Quat4dNAN );
    _realloc0( Ksp       , nnode, Quat4dNAN );
    _realloc0( Kpp       , nnode, Quat4dNAN );

    // Additional:
    // Angles
    _realloc0( angles, nnode*6, Vec3dNAN );   // 6=4*3/2
    // Torsions
    _realloc0( tors2atom,  ntors, Quat4iZero );
    _realloc0( torsParams, ntors, Quat4dNAN  ); 

    _realloc0( constr    , natoms, Quat4dOnes*-1. );
}

void clone( MMFFsp3_loc& from, bool bRealloc, bool bREQsDeep=true ){
    realloc( from.nnode, from.ncap  );
    lvec   =from.lvec;
    invLvec=from.invLvec;
    for(int i=0; i<nDOFs; i++){
        DOFs[i]=from. DOFs[i];
     //fDOFs[i]=from.fDOFs[i];
    }
    for(int i=0; i<natoms; i++){
        atypes   [i]=from.atypes   [i];
        neighs   [i]=from.neighs   [i];
        neighCell[i]=from.neighCell[i];
        bkneighs [i]=from.bkneighs [i];
        constr   [i]=from.constr   [i];
    }
    for(int i=0; i<nnode; i++){
        apars[i]=from.apars[i];
        bLs  [i]=from.bLs  [i];    
        bKs  [i]=from.bKs  [i];    
        Ksp  [i]=from.Ksp  [i]; 
        Kpp  [i]=from.Kpp  [i];    
    }
    if(from.REQs){
        if(bREQsDeep){ _realloc(REQs, natoms ); for(int i=0; i<natoms; i++){  REQs[i]=from.REQs[i]; } }
        else         {          REQs=from.REQs;                                                       }
    }
}

void dealloc(){
    _dealloc(DOFs );
    _dealloc(fDOFs);
    apos   = 0;
    fapos  = 0;
    pipos  = 0;
    fpipos = 0;
    _dealloc(atypes);
    _dealloc(neighs);
    _dealloc(neighCell);
    _dealloc(bkneighs);
    _dealloc(apars);
    _dealloc(bLs);
    _dealloc(bKs);
    _dealloc(Ksp);
    _dealloc(Kpp);
    _dealloc(angles);

    _dealloc(tors2atom  );
    _dealloc(torsParams );
    nnode=0; ncap=0; ntors=0; natoms=0; nvecs =0; nDOFs =0; int ipi0=natoms;
}

void setLvec(const Mat3d& lvec_){ lvec=lvec_; lvec.invert_T_to( invLvec ); }

double optimalTimeStep(double m=1.0){
    double Kmax = 1.0;
    for(int i=0; i<nnode; i++){ 
        Kmax=fmax(Kmax, bKs[i].x ); 
        Kmax=fmax(Kmax, bKs[i].y ); 
        Kmax=fmax(Kmax, bKs[i].z ); 
        Kmax=fmax(Kmax, bKs[i].w ); 
    }
    return M_PI*2.0*sqrt(m/Kmax)/10.0;  // dt=T/10;   T = 2*pi/omega = 2*pi*sqrt(m/k)
}

// ============== Evaluation

double eval_atom(const int ia){
    //printf( "MMFFsp3_loc::eval_atom(%i)\n", ia );
    double E=0;
    const Vec3d pa  = apos [ia]; 
    const Vec3d hpi = pipos[ia]; 

    //printf( "apos[%i](%g,%g,%g)\n",ia,apos[ia].x,apos[ia].y,apos[ia].z );
    //return E;

    //Vec3d& fa  = fapos [ia]; 
    //Vec3d& fpi = fpipos[ia];
    Vec3d fa   = Vec3dZero;
    Vec3d fpi  = Vec3dZero; 
    
    //--- array aliases
    const int*    ings = neighs   [ia].array;
    const int*    ingC = neighCell[ia].array;
    const double* bK   = bKs      [ia].array;
    const double* bL   = bLs      [ia].array;
    const double* Kspi = Ksp      [ia].array;
    const double* Kppi = Kpp      [ia].array;
    Vec3d* fbs  = fneigh   +ia*4;
    Vec3d* fps  = fneighpi +ia*4;

    // // --- settings
    // double  ssC0 = apars[ia].x;
    // double  ssK  = apars[ia].y;
    // double  piC0 = apars[ia].z;
    // //bool    bPi  = ings[3]<0;   we distinguish this by Ksp, otherwise it would be difficult for electron pairs e.g. (-O-C=)

    const Quat4d& apar  = apars[ia];
    const double  piC0 = apar.w;

    //printf( "ang0 %g cs0(%g,%g)\n", atan2(cs0_ss.y,cs0_ss.x)*180/M_PI, cs0_ss.x,cs0_ss.x );

    //--- Aux Variables 
    Quat4d  hs[4];
    Vec3d   f1,f2;
    
    //bool bErr=0;
    //const int ia_DBG = 0;
    //if(ia==ia_DBG)printf( "ffl[%i] neighs(%i,%i,%i,%i) \n", ia, ings[0],ings[1],ings[2],ings[3] );
    //printf("lvec %i %i \n", ia, ia_DBG ); // printMat(lvec);

    //if( ia==5 ){ printf( "ffls[%2i] atom[%2i] ng(%3i,%3i,%3i,%3i) ngC(%3i,%3i,%3i,%3i) shifts=%li bPBC=%i\n", id, ia,   ings[0],ings[1],ings[2],ings[3],   ingC[0],ingC[1],ingC[2],ingC[3], (long)shifts, bPBC ); };

    for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break

    // double Eb=0;
    // double Ea=0;
    // double EppI=0;
    // double Eps=0;

    // --------- Bonds Step
    for(int i=0; i<4; i++){
        int ing = ings[i];
        //printf( "bond[%i|%i=%i]\n", ia,i,ing );
        //fbs[i]=Vec3dOne; fps[i]=Vec3dOne;
        //fbs[i]=Vec3dZero; fps[i]=Vec3dZero; // NOTE: wee need to initialize it before, because of the break
        if(ing<0) break;

        //printf("ia %i ing %i \n", ia, ing ); 
        Vec3d  pi = apos[ing];
        Quat4d h; 
        h.f.set_sub( pi, pa );
        //if(idebug)printf( "bond[%i|%i=%i] l=%g pj[%i](%g,%g,%g) pi[%i](%g,%g,%g)\n", ia,i,ing, h.f.norm(), ing,apos[ing].x,apos[ing].y,apos[ing].z, ia,pa.x,pa.y,pa.z  );
        
        //Vec3d h_bak = h.f;    
        //shifts=0; 
        if(bPBC){   
            if(shifts){
                int ipbc = ingC[i]; 
                //Vec3d sh = shifts[ipbc]; //apbc[i]  = pi + sh;
                h.f.add( shifts[ipbc] );
                //if( (ipbc!=4) ){ printf("ffls[%i] atom[%i,%i=%i] ipbc %i shifts(%g,%g,%g)\n", id, ia,i,ing, ipbc, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z); };
                //if( (ipbc!=4) ){ printf("ffls[%i] atom[%i,%i=%i] ipbc %i shifts(%g,%g,%g)\n", id, ia,i,ing, ipbc, shifts[ipbc].x,shifts[ipbc].y,shifts[ipbc].z); };
            }else{
                Vec3i g  = invLvec.nearestCell( h.f );
                // if(ia==ia_DBG){
                //     Vec3d u; invLvec.dot_to(h.f,u);
                //     printf( "CPU:bond[%i,%i] u(%6.3f,%6.3f,%6.3f) shi(%6.3f,%6.3f,%6.3f) \n", ia, ing, u.x,u.y,u.z,   (float)g.x,(float)g.y,(float)g.z );
                // }
                Vec3d sh = lvec.a*g.x + lvec.b*g.y + lvec.c*g.z;
                h.f.add( sh );
                //apbc[i] = pi + sh;
            }
        }

        //wrapBondVec( h.f );
        //printf( "h[%i,%i] r_old %g r_new %g \n", ia, ing, h_bak.norm(), h.f.norm() );
        double l = h.f.normalize();

        h.e    = 1/l;
        hs [i] = h;

        //bErr|=ckeckNaN( 1,4, (double*)(hs+i), [&]{ printf("atom[%i]hs[%i]",ia,i); } );
        // bond length force
        //continue; 

        //if(ia==ia_DBG) printf( "ffl:h[%i|%i=%i] l %g h(%g,%g,%g) pj(%g,%g,%g) pa(%g,%g,%g) \n", ia,i,ing, l, h.x,h.y,h.z, apos[ing].x,apos[ing].y,apos[ing].z, pa.x,pa.y,pa.z );

        if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
            if(doBonds){
                E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);

                //double Ebi = evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);
                //E +=Ebi; 
                //Eb+=Ebi; 

                //printf( "ffl:bond[%i|%i=%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,i,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
                //if(ia==ia_DBG)printf( "ffl:bond[%i|%i=%i] kb=%g l0=%g l=%g h(%g,%g,%g) f(%g,%g,%g) \n", ia,i,ing, bK[i],bL[i], l, h.x,h.y,h.z,  f1.x,f1.y,f1.z  );
                //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fbond[%i]",ia,i); } );
            }

            double kpp = Kppi[i];
            if( (doPiPiI) && (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                
                //double EppIi = evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                //E +=EppIi; 
                //EppI+=EppIi;
                
                //if(verbosity>0)printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g), EppI=%g\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z, EppIi  );
                //printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g), EppI=%g, Epp=%g\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z, EppIi, EppI  );
                //if(ia==ia_DBG)printf( "ffl:pp[%i|%i] kpp=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, kpp, hpi.dot(pipos[ing]), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
                //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fpp1[%i]",ia,i); } );
                //bErr|=ckeckNaN( 1,3, (double*)&f2, [&]{ printf("atom[%i]fpp2[%i]",ia,i); } );
            }
            // ToDo: triple bonds ?
            
        } 
        
        // pi-sigma 
        //if(bPi){    
        double ksp = Kspi[i];
        if( doPiSigma && (ksp>1e-6) ){  
            E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)

            //double Epsi = evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)
            //E +=Epsi; 
            //Eps+=Epsi;

            //printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z  );
            //if(kpp<-1e-6){  E += evalPiAling( hpi, pipos[ing], 1., 1., -kpp, f1, f2 );   fpi.add(f1);  fps[i].add(f2);  }    //   align pi-electron pair (e.g. in Nitrogen)
            //if(ia==ia_DBG)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g hp(%g,%g,%g) h(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), hpi.x,hpi.y,hpi.z,  h.x,h.y,h.z  );
            //if(ia==ia_DBG)printf( "ffl:sp[%i|%i] ksp=%g piC0=%g c=%g f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing, ksp,piC0, hpi.dot(h.f), f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fsp1[%i]",ia,i); } );
            //bErr|=ckeckNaN( 1,3, (double*)&f2, [&]{ printf("atom[%i]fsp2[%i]",ia,i); } );
        }
        //}
        
    }

    //printf( "MMFF_atom[%i] cs(%6.3f,%6.3f) ang=%g [deg]\n", ia, cs0_ss.x, cs0_ss.y, atan2(cs0_ss.y,cs0_ss.x)*180./M_PI );
    // --------- Angle Step
    const double R2damp=Rdamp*Rdamp;
    if(doAngles){

    double  ssK,ssC0;
    Vec2d   cs0_ss;
    Vec3d*  angles_i;
    if(bEachAngle){
        angles_i = angles+(ia*6);
    }else{
        ssK    = apar.z;
        cs0_ss = Vec2d{apar.x,apar.y};
        ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
    }

    int iang=0;
    for(int i=0; i<3; i++){
        int ing = ings[i];
        if(ing<0) break;
        const Quat4d& hi = hs[i];
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            const Quat4d& hj = hs[j];    
            if(bEachAngle){
                // 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
                //  0    1    2    3    4    5 
                cs0_ss = angles_i[iang].xy();  
                ssK    = angles_i[iang].z;
                //printf( "types{%i,%i,%i} ssK %g cs0(%g,%g) \n", atypes[ing], atypes[ia], atypes[jng], ssK, cs0_ss.x, cs0_ss.y  );
                iang++; 
            };
            //bAngleCosHalf = false;
            //double Eai;
            //printf( "bAngleCosHalf= %i\n", bAngleCosHalf);
            if( bAngleCosHalf ){
                E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
                //Eai = evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
            }else{ 
                E += evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
                //Eai = evalAngleCos( hi.f, hj.f, hi.e, hj.e, ssK, ssC0, f1, f2 );     // angles between sigma bonds
            }

            //E +=Eai; 
            //Ea+=Eai;
            
            //printf( "ffl:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //if(ia==ia_DBG)printf( "ffl:ang[%i|%i,%i] kss=%g cs0(%g,%g) c=%g l(%g,%g) f1(%g,%g,%g) f2(%g,%g,%g)\n", ia,ing,jng, ssK, cs0_ss.x,cs0_ss.y, hi.f.dot(hj.f),hi.w,hj.w, f1.x,f1.y,f1.z,  f2.x,f2.y,f2.z  );
            //bErr|=ckeckNaN( 1,3, (double*)&f1, [&]{ printf("atom[%i]fss1[%i,%i]",ia,i,j); } );
            //bErr|=ckeckNaN( 1,3, (double*)&f2, [&]{ printf("atom[%i]fss2[%i,%i]",ia,i,j); } );
            fa    .sub( f1+f2  );
            // ----- Error is HERE
            if(bSubtractAngleNonBond){
                Vec3d fij=Vec3dZero;
                //Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);
                Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                //Vec3d dp   = hj.f*(1./hj.w) - hi.f*(1./hi.w);
                //Vec3d dp   = apbc[j] - apbc[i];
                E -= getLJQH( dp, fij, REQij, R2damp );
                //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                //bErr|=ckeckNaN( 1,3, (double*)&fij, [&]{ printf("atom[%i]fLJ2[%i,%i]",ia,i,j); } );
                f1.sub(fij);
                f2.add(fij);
            }
            fbs[i].add( f1     );
            fbs[j].add( f2     );
            //if(ia==ia_DBG)printf( "ffl:ANG[%i|%i,%i] fa(%g,%g,%g) fbs[%i](%g,%g,%g) fbs[%i](%g,%g,%g)\n", ia,ing,jng, fa.x,fa.y,fa.z, i,fbs[i].x,fbs[i].y,fbs[i].z,   j,fbs[j].x,fbs[j].y,fbs[j].z  );
            // ToDo: subtract non-covalent interactions
        }
    }
    }    
    //if(bErr){ printf("ERROR in ffl.eval_atom[%i] => Exit() \n", ia ); exit(0); }
    
    /*
    double Kfix = constr[ia].w;
    if(Kfix>0){  
        //printf( "applyConstrain(i=%i,K=%g)\n", ia, Kfix );
        Vec3d d = constr[ia].f-pa;
        d.mul( Kfix );
        double fr2 = d.norm2();
        double F2max = 1.0;
        if( fr2>F2max ){
            d.mul( sqrt(F2max/fr2) );
        }
        fa.add( d );
        E += d.norm()*Kfix*0.5;
    }
    */
    
    //fapos [ia].add(fa ); 
    //fpipos[ia].add(fpi);
    fapos [ia]=fa; 
    fpipos[ia]=fpi;

    //printf( "MYOUTPUT atom %i Ebond= %g Eangle= %g Edihed= %g Eimpr = %g Etot=%g\n", ia+1, Eb, Ea, EppI, Eps, E );
    //fprintf( file, "atom %i Ebond= %g Eangle= %g Edihed= %g Eimpr= %g Etot=%g \n", ia+1, Eb, Ea, EppI, Eps, E );
    

    return E;
}



    double eval_atom_opt(const int ia){

        double E=0;
        const Vec3d pa   = apos [ia]; 
        Vec3d       fa   = Vec3dZero;
        const int*  ings = neighs[ia].array;

        Vec3d* fbs  = fneigh   +ia*4;

        const Quat4d& apar  = apars[ia];
        const double  piC0  = apar.w;

        //--- Aux Variables 
        Quat4d  hs[4];
        Vec3d   f1,f2;

        // ========================= Bonds
        { // Bonds

            const Vec3d hpi = pipos[ia]; 
            Vec3d fpi       = Vec3dZero; 
            const int* ingC = neighCell[ia].array;

            const double* bK   = bKs      [ia].array;
            const double* bL   = bLs      [ia].array;
            const double* Kspi = Ksp      [ia].array;
            const double* Kppi = Kpp      [ia].array;

            Vec3d* fps  = fneighpi +ia*4;
            for(int i=0; i<4; i++){ fbs[i]=Vec3dZero; fps[i]=Vec3dZero; } // we initialize it here because of the break

            for(int i=0; i<4; i++){
                const int ing = ings[i];
                if(ing<0) break;

                const Vec3d  pi = apos[ing];
                Quat4d h; 
                h.f.set_sub( pi, pa );
                
                const int ipbc = ingC[i]; 
                //Vec3d sh = shifts[ipbc]; //apbc[i]  = pi + sh;
                h.f.add( shifts[ipbc] );
                const double l = h.f.normalize();
                h.e    = 1/l;
                hs [i] = h;

                if(ia<ing){   // we should avoid double counting because otherwise node atoms would be computed 2x, but capping only once
                    E+= evalBond( h.f, l-bL[i], bK[i], f1 ); fbs[i].sub(f1);  fa.add(f1);

                    const double kpp = Kppi[i];
                    if( (ing<nnode) && (kpp>1e-6) ){   // Only node atoms have pi-pi alignemnt interaction
                        E += evalPiAling( hpi, pipos[ing], 1., 1.,   kpp,       f1, f2 );   fpi.add(f1);  fps[i].add(f2);    //   pi-alignment     (konjugation)
                    }
                } 

                const double ksp = Kspi[i];
                if( (ksp>1e-6) ){  
                    E += evalAngleCos( hpi, h.f      , 1., h.e, ksp, piC0, f1, f2 );   fpi.add(f1); fa.sub(f2);  fbs[i].add(f2);       //   pi-planarization (orthogonality)
                }            
            }
            fpipos[ia]=fpi;

        } // Bonds


        // ========================= Angles
        const double R2damp=Rdamp*Rdamp;
        if(doAngles){

            /*
            double  ssK,ssC0;
            Vec2d   cs0_ss;
            Vec3d*  angles_i;
            if(bEachAngle){
                Vec3d*   angles_i = angles+(ia*6);
            }else{
                ssK    = apar.z;
                cs0_ss = Vec2d{apar.x,apar.y};
                ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf
            }
            */

            const double ssK    = apar.z;
            const Vec2d  cs0_ss = Vec2d{apar.x,apar.y};
            const double ssC0   = cs0_ss.x*cs0_ss.x - cs0_ss.y*cs0_ss.y;   // cos(2x) = cos(x)^2 - sin(x)^2, because we store cos(ang0/2) to use in  evalAngleCosHalf

            int iang=0;
            for(int i=0; i<3; i++){
                int ing = ings[i];
                if(ing<0) break;
                const Quat4d& hi = hs[i];
                for(int j=i+1; j<4; j++){
                    int jng  = ings[j];
                    if(jng<0) break;
                    const Quat4d& hj = hs[j];    

                    /*
                    if(bEachAngle){
                        // 0-1, 0-2, 0-3, 1-2, 1-3, 2-3
                        //  0    1    2    3    4    5 
                        cs0_ss = angles_i[iang].xy();  
                        ssK    = angles_i[iang].z;
                        iang++; 
                    };
                    */
                    E += evalAngleCosHalf( hi.f, hj.f,  hi.e, hj.e,  cs0_ss,  ssK, f1, f2 );
                    fa    .sub( f1+f2  );

                    /*
                    // ----- Error is HERE
                    if(bSubtractAngleNonBond){
                        Vec3d fij=Vec3dZero;
                        //Quat4d REQij; combineREQ( REQs[ing],REQs[jng], REQij );
                        Quat4d REQij = _mixREQ(REQs[ing],REQs[jng]);
                        Vec3d dp; dp.set_lincomb( 1./hj.w, hj.f,  -1./hi.w, hi.f );
                        //Vec3d dp   = hj.f*(1./hj.w) - hi.f*(1./hi.w);
                        //Vec3d dp   = apbc[j] - apbc[i];
                        E -= getLJQH( dp, fij, REQij, R2damp );
                        //if(ia==ia_DBG)printf( "ffl:LJQ[%i|%i,%i] r=%g REQ(%g,%g,%g) fij(%g,%g,%g)\n", ia,ing,jng, dp.norm(), REQij.x,REQij.y,REQij.z, fij.x,fij.y,fij.z );
                        //bErr|=ckeckNaN( 1,3, (double*)&fij, [&]{ printf("atom[%i]fLJ2[%i,%i]",ia,i,j); } );
                        f1.sub(fij);
                        f2.add(fij);
                    }
                    */

                    fbs[i].add( f1     );
                    fbs[j].add( f2     );

                }
            }
        }    // doAngles
        fapos [ia]=fa; 
        return E;   
    }





double eval_atoms(){
    //FILE *file = fopen("out","w");
    //Etot,
    //Eb=0;Ea=0;Eps=0;EppT=0;EppI=0;
    double E=0;
    for(int ia=0; ia<nnode; ia++){ 
        E+=eval_atom(ia); 
    }
    //printf( "MYOUTPUT Ebond= %g Eangle= %g Edihed= %g Eimpr = %g Etot=%g\n", Eb, Ea, EppI, Eps, E );
    //fprintf( file, "Ebond= %g Eangle= %g Edihed= %g Eimpr= %g Etot=%g \n", Eb, Ea, EppI, Eps, E );
    //fclose(file);
    return E;
}

double eval_torsion(int it){
    printf( "MMFFsp3_loc::eval_torsion(%i)\n", it );
    Quat4i ias = tors2atom [it];
    Quat4d par = torsParams[it];

    Vec3d ha    = apos[ ias.x ] - apos[ ias.y ];
    Vec3d hb    = apos[ ias.w ] - apos[ ias.z ];
    Vec3d hab   = apos[ ias.z ] - apos[ ias.y ];

    double ila  = 1/ha.normalize();
    double ilb  = 1/hb.normalize();
    double ilab = 1/hab.normalize();

    double ca   = hab.dot(ha);
    double cb   = hab.dot(hb);
    double cab  = ha .dot(hb);
    double sa2  = (1-ca*ca);
    double sb2  = (1-cb*cb);
    double invs = 1/sqrt( sa2*sb2 );
    //double c    = ;  //  c = <  ha - <ha|hab>hab   | hb - <hb|hab>hab    >

    Vec2d cs,csn;
    cs.x = ( cab - ca*cb )*invs;
    cs.y = sqrt(1-cs.x*cs.x); // can we avoid this sqrt ?
    cs.udiv_cmplx( par.xy() ); 

    const int n = (int)par.w; // I know it is stupid store n as double, but I don't make another integer arrays just for it
    for(int i=0; i<n-1; i++){
        csn.mul_cmplx(cs);
    }

    // check here : https://www.wolframalpha.com/input/?i=(x+%2B+isqrt(1-x%5E2))%5En+derivative+by+x

    const double k = par.z;
    double E       = k  *(1-csn.x);
    double dcn     = k*n*   csn.x;
    //double fr  =  k*n*    csn.y;

    //double c   = cos_func(ca,cb,cab);

    //printf( "<fa|fb> %g cT %g cS %g \n", cs.x, cT, cS );

    // derivatives to get forces

    double invs2 = invs*invs;
    dcn *= invs;
    double dcab  = dcn;                          // dc/dcab = dc/d<ha|hb>
    double dca   = (1-cb*cb)*(ca*cab - cb)*dcn;  // dc/dca  = dc/d<ha|hab>
    double dcb   = (1-ca*ca)*(cb*cab - ca)*dcn;  // dc/dca  = dc/d<hb|hab>

    Vec3d fa,fb,fab;

    fa =Vec3dZero;
    fb =Vec3dZero;
    fab=Vec3dZero;

    //Mat3Sd J;
    SMat3d J;

    J.from_dhat(ha);    // -- by ha
    J.mad_ddot(hab,fa, dca ); // dca /dha = d<ha|hab>/dha
    J.mad_ddot(hb ,fa, dcab); // dcab/dha = d<ha|hb> /dha

    J.from_dhat(hb);    // -- by hb
    J.mad_ddot(hab,fb, dcb ); // dcb /dhb = d<hb|hab>/dha
    J.mad_ddot(ha ,fb, dcab); // dcab/dhb = d<hb|ha> /dha

    J.from_dhat(hab);         // -- by hab
    J.mad_ddot(ha,fab, dca);  // dca/dhab = d<ha|hab>/dhab
    J.mad_ddot(hb,fab, dcb);  // dcb/dhab = d<hb|hab>/dhab
    // derivative cab = <ha|hb>

    fa .mul( ila  );
    fb .mul( ilb  );
    fab.mul( ilab );

    // ToDo : Which order ?
    fapos[ias.x].sub( fa );
    fapos[ias.y].add( fa  - fab );
    fapos[ias.z].add( fab - fb  );
    fapos[ias.w].add( fb );

    /*
    {
        double fsc = 100.0;
        glColor3f(1.0,0.0,1.0);  Draw3D::drawVecInPos(hab.normalized()*E, apos[ia.x]);
        glColor3f(1.0,1.0,1.0);  Draw3D::drawVecInPos(ha.normalized()*E, apos[ia.x]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa*fsc, apos[ia.x]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fa*fsc*-1, apos[ia.y]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fb*fsc*-1, apos[ia.z]);
        glColor3f(1.0,0.0,0.0);  Draw3D::drawVecInPos(fb*fsc, apos[ia.w]);

    }
    */

    return E;
}

double eval_torsions(){
    printf( "MMFFsp3_loc::eval_torsions(ntors=%i)\n", ntors );
    double E=0;
    for(int it=0; it<ntors; it++){ 
        E+=eval_torsion(it); 
    }
    return E;
}


void addjustAtomCapLenghs(int ia){
    const Vec3d pa  = apos [ia]; 
    const int*    ings = neighs   [ia].array;
    //const int*    ingC = neighCell[ia].array;
    const double* bL   = bLs      [ia].array;
    for(int i=0; i<4; i++){
        int ing = ings[i];
        if(ing<nnode)continue;
        Vec3d  pi = apos[ing];
        Vec3d h; 
        h.set_sub( pi, pa );
        h.mul(bL[ing]/h.norm());
        apos[ing].set_add( pa, h ); 
    } 
}

void addjustCapLenghs(){
    for(int ia=0; ia<nnode; ia++){ 
        addjustAtomCapLenghs(ia);
    }
}

void initPi( Vec3d* pbc_shifts, double Kmin=0.0001, double r2min=1e-4, bool bCheck=true ){
    //printf( "MMFFsp3_loc::initPi()\n" );
    for(int ia=0; ia<nnode; ia++){ 
        if(vapos)vapos[natoms+ia]=Vec3dZero;
        const int*    ngs = neighs   [ia].array;
        const int*    ngC = neighCell[ia].array;
        const double* ks  = Ksp      [ia].array;
        Vec3d u,v,p;
        int nfound=0; 
        int j=0;
        while(j<4){ if(ks[j]>Kmin){ u=apos[ngs[j]]+pbc_shifts[ngC[j]]; nfound++; j++;break; }; j++; }
        while(j<4){ if(ks[j]>Kmin){ v=apos[ngs[j]]+pbc_shifts[ngC[j]]; nfound++; j++;break; }; j++; }
        while(j<4){ if(ks[j]>Kmin){ p=apos[ngs[j]]+pbc_shifts[ngC[j]]; nfound++; j++;break; }; j++; }
        //printf("...ia=%i nfound=%i\n",ia,nfound);
        if      ( nfound>=2 ){
            if( nfound==2 ){  p=apos[ia]; }
            //printf( "::[%i] u(%6.3f,%6.3f,%6.3f) v(%6.3f,%6.3f,%6.3f) p(%6.3f,%6.3f,%6.3f) nf=%1i \n", ia, u.x,u.y,u.z,  v.x,v.y,v.z, p.x,p.y,p.z, nfound );
            u.sub(p); v.sub(p);
            Vec3d pi; pi.set_cross( u,v );
            double r2 = pi.norm2();
            //printf( "::pi[%i](%6.3f,%6.3f,%6.3f) u(%6.3f,%6.3f,%6.3f) v(%6.3f,%6.3f,%6.3f) nf=%1i r2 %g \n", ia, u.x,u.y,u.z,  v.x,v.y,v.z, pi.x,pi.y,pi.z, nfound, sqrt(r2) );
            if( r2>r2min ){
                pi.mul(1/sqrt(r2));
                pipos[ia]=pi;
                //printf( "pi[%i]=pi(%5.3f,%5.3f,%5.3f) nfound %1i r2 %g \n", ia, pi.x,pi.y,pi.z, nfound, sqrt(r2) );
                continue;
            }
        }
        //printf( "pi[%i]=Vec3dZero nf=%1i \n", ia, nfound );
        //pipos[ia]=Vec3dZ; // pi cannot be define
        pipos[ia]=Vec3dZero; // pi cannot be define
    }
    for(int ia=0; ia<nnode; ia++){ 
        const int*    ngs = neighs   [ia].array;
        const double* ks  = Kpp      [ia].array;
        int imax=-1;
        double kmax=0.0;
        double r2 = pipos[ia].norm2();
        if( r2>0.1 )continue;
        for(int i=0; i<4; i++){
            double k = ks[i];
            if(k>kmax){ 
                if(bCheck){
                    int ing=ngs[i];
                    if     ( ing<0     ){ printf("ERROR in MMFFsp3_loc::initPi() atom[%i] Kpp=%g but neihgbor[%i] is undefined(%i) => Exit() \n", ia, k,i, ing    ); exit(0); }    
                    else if( ing>nnode ){ printf("ERROR in MMFFsp3_loc::initPi() atom[%i] Kpp=%g but neihgbor[%i] is cap(nnode=%i) => Exit() \n", ia, k,i, nnode  ); exit(0); }
                }
                imax=i; kmax=k; 
            };
        }
        if( kmax>Kmin ){ 
            pipos[ia]=pipos[ ngs[imax] ]; 
            //printf("pi[%i] set by neigh[%i==%i] with Kpp=%g\n", ia, imax, ngs[imax], kmax ); 
        }else{
            pipos[ia]=Vec3dZ; // pi cannot be define
        }
        
    }
    //for(int ia=0; ia<nnode; ia++){  pipos[ia]=Vec3dZ ; } // DEBUG
}

int relax_pi( int niter, double dt, double Fconv, double Flim=1000.0 ){
    double F2conv = Fconv*Fconv;
    double E=0,F2=0;
    int    itr=0;
    for(itr=0; itr<niter; itr++){
        {E=0;F2=0;}
        cleanForce();
        normalizePis();
        {Eb=0;Ea=0;Eps=0;EppT=0;EppI=0;}
        for(int ia=0; ia<nnode; ia++){ 
            E += eval_atom(ia);
        }
        for(int ia=0; ia<nnode; ia++){
            assemble_atom( ia );
        }
        for(int i=natoms; i<nvecs; i++){
            F2 +=move_atom_MD( i, dt, Flim, 0.9 );
        }
        if(F2<F2conv){ return itr; }
    }
    return itr;
}

void normalizePis(){ 
    for(int i=0; i<nnode; i++){ pipos[i].normalize(); } 
}

void constrainAtom( int ia, double Kfix=1.0 ){
    printf( "constrainAtom(i=%i,K=%g)\n", ia, Kfix );
    constr[ia].f=apos[ia];
    constr[ia].w=Kfix;
};

void cleanForce(){ 
    Etot=0;
    //for(int i=0; i<natoms; i++){ fapos [i].set(0.0);  } 
    //for(int i=0; i<nnode;  i++){ fpipos[i].set(0.0);  } 
    for(int i=0; i<nDOFs; i++){ fDOFs[i]=0;  } 
    // NOTE: We do not need clean fneigh,fneighpi because they are set in eval_atoms 
}

void assemble_atom(int ia){
    Vec3d fa=Vec3dZero,fp=Vec3dZero;
    const int* ings = bkneighs[ia].array;
    bool bpi = ia<nnode;
    for(int i=0; i<4; i++){
        int j = ings[i];
        if(j<0) break;
        //if(j>=(nnode*4)){ printf("ERROR bkngs[%i|%i] %i>=4*nnode(%i)\n", ia, i, j, nnode*4 ); exit(0); }
        fa.add(fneigh  [j]);
        if(bpi){
            //printf( "assemble[%i,%i|%i] pi(%g,%g,%g) fp(%g,%g,%g) fpng(%g,%g,%g) \n", ia,i,j, pipos[ia].x,pipos[ia].y,pipos[ia].z, fpipos[ia].x,fpipos[ia].y,fpipos[ia].z, fneighpi[j].x,fneighpi[j].y,fneighpi[j].z );
            fp.add(fneighpi[j]);
            //ckeckNaN( 1,3, (double*)(fneighpi+j), [&]{ printf("assemble.fneighpi[%i]",j); } );
        }
    }
    fapos [ia].add( fa ); 

    if(bpi){
        //if( pipos[ia].norm2()>1.5 ){ printf("pipos[%i](%g,%g,%g) not normalized !!! (assemble.iteration=%i) => Exit() \n", ia, pipos[ia].x,pipos[ia].y,pipos[ia].z, itr_DBG ); exit(0); };
        fpipos[ia].add( fp );
        fpipos[ia].makeOrthoU( pipos[ia] );  // subtract force component which change pi-vector size
        //ckeckNaN( 1,3, (double*)(fpipos+ia), [&]{ printf("assemble.makeOrthoU[%i] c %g pipos(%g,%g,%g) fpipos(%g,%g,%g)",ia,  pi.dot(fi),   pi.x,pi.y,pi.z, fi.x,fi.y,fi.z  ); } );
    }
}

void asseble_forces(){
    for(int ia=0; ia<natoms; ia++){
        assemble_atom(ia);
    }
}

double eval( bool bClean=true, bool bCheck=true ){
    //if(bClean){ cleanAll(); }
    //printf( "print_apos() BEFORE\n" );print_apos();
    if(bClean)cleanForce();
    normalizePis();
    //printf( "print_apos() AFTER \n" ); print_apos();
    Etot += eval_atoms();
    //if(idebug){printf("CPU BEFORE assemble() \n"); printDEBUG();} 
    asseble_forces();
    if(bTorsion) EppI +=eval_torsions();
    //Etot = Eb + Ea + Eps + EppT + EppI;
    return Etot;
}

double eval_check(){
    if(verbosity>0){
        printf(" ============ check MMFFsp3_loc START\n " );
        printSizes();
    }
    //print_pipos();
    eval();
    checkNans();
    if(verbosity>0)printf(" ============ check MMFFsp3_loc DONE\n " );
    return Etot;
}

// ToDo: OpenMP paraelization atempt
int run( int niter, double dt, double Fconv, double Flim ){
    double F2conv = Fconv*Fconv;
    double E,F2;
    int    itr;
    //if(itr_DBG==0)print_pipos();
    //bool bErr=0;
    for(itr=0; itr<niter; itr++){
        E=0;
        // ------ eval MMFF
        for(int ia=0; ia<natoms; ia++){ 
            if(ia<nnode)E += eval_atom(ia);
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+ia), [&]{ printf("eval.MMFF[%i]",ia); } );
            E += evalLJQs_ng4_PBC_atom( ia ); 
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+ia), [&]{ printf("eval.NBFF[%i]",ia); } );
        }
        // ---- assemble (we need to wait when all atoms are evaluated)
        for(int ia=0; ia<natoms; ia++){
            assemble_atom( ia );
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+ia), [&]{ printf("assemble[%i]",ia); } );
        }
        // ------ move
        F2=0;
        for(int i=0; i<nvecs; i++){
            //F2 += move_atom_GD( i, dt, Flim );
            //bErr|=ckeckNaN( 1,3, (double*)(fapos+i), [&]{ printf("move[%i]",i); } );
            F2 += move_atom_MD( i, dt, Flim, 0.99 );
            //F2 += move_atom_kvaziFIRE( i, dt, Flim );
        }
        if(F2<F2conv)break;
        //itr_DBG++;
    }
    return itr;
}

int run_omp( int niter, double dt, double Fconv, double Flim ){
    double F2conv = Fconv*Fconv;
    double E=0,F2=0;
    int    itr=0;
    #pragma omp parallel shared(E,F2) private(itr)
    for(itr=0; itr<niter; itr++){
        // This {} should be done just by one of the processors
        #pragma omp single
        {E=0;F2=0;}
        // ------ eval MMFF
        #pragma omp for reduction(+:E)
        for(int ia=0; ia<natoms; ia++){ 
            if(verbosity>3)printf( "atom[%i]@cpu[%i/%i]\n", ia, omp_get_thread_num(), omp_get_num_threads()  );
            if(ia<nnode)E += eval_atom(ia);
            //E += evalLJQs_ng4_PBC_atom( ia ); 
            E += evalLJQs_ng4_PBC_atom_omp( ia ); 
        }
        // ---- assemble (we need to wait when all atoms are evaluated)
        #pragma omp for
        for(int ia=0; ia<natoms; ia++){
            assemble_atom( ia );
        }
        // ------ move
        #pragma omp for reduction(+:F2)
        for(int i=0; i<nvecs; i++){
            F2 += move_atom_MD( i, dt, Flim, 0.99 );
        }
        for(int i=natoms; i<nvecs; i++){
             F2 += move_atom_MD( i, dt, Flim, 0.99 );
        }
        #pragma omp single
        if(verbosity>2){printf( "step[%i] E %g |F| %g ncpu[%i] \n", itr, E, F2, omp_get_num_threads() );}
    }
    return itr;
}

void flipPis( Vec3d ax ){
    for(int i=0; i<nnode; i++){
        double c = pipos[i].dot(ax);
        if( c<0 ){ pipos[i].mul(-1); } 
    }
}

inline double move_atom_GD(int i, float dt, double Flim){
    Vec3d  f   = fapos[i];
    Vec3d  p = apos [i];
    double fr2 = f.norm2();
    if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    if(bPi){ f.add_mul( p, -p.dot(f) ); }
    p.add_mul( f, dt );
    if(bPi)p.normalize();
    apos [i] = p;
    //fapos[i] = Vec3dZero;
    return fr2;
}

inline double move_atom_MD( int i, const float dt, const double Flim, const double cdamp=0.9 ){
    Vec3d  f = fapos[i];
    Vec3d  v = vapos[i];
    Vec3d  p = apos [i];
    const double fr2 = f.norm2();
    //if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    // bool b=false;
    // b|=ckeckNaN_d( 1,3, (double*)&p, "p.1" );
    // b|=ckeckNaN_d( 1,3, (double*)&v, "v.1" );
    // b|=ckeckNaN_d( 1,3, (double*)&f, "f.1" );
    if(bPi)f.add_mul( p, -p.dot(f) );           //b|=ckeckNaN_d( 1,3, (double*)&f, "f.2" );
    v.mul    ( cdamp );
    v.add_mul( f, dt );                         //b|=ckeckNaN_d( 1,3, (double*)&v, "v.2" );
    if(bPi)v.add_mul( p, -p.dot(v) );           //b|=ckeckNaN_d( 1,3, (double*)&v, "v.3" );
    p.add_mul( v, dt );
    if(bPi)p.normalize();
    // if( bPi &&( p.norm2()>1.5 )){ printf("pipos[%i/%i](%g,%g,%g) not normalized !!! (move.iteration=%i) => Exit() \n", i, natoms, p.x,p.y,p.z, itr_DBG ); exit(0); };
    // b|=ckeckNaN_d( 1,3, (double*)&p, "p.4" );
    // b|=ckeckNaN_d( 1,3, (double*)&v, "v.4" );
    // b|=ckeckNaN_d( 1,3, (double*)&f, "f.4" );
    // if(b){ printf("ERROR NaNs in move_atom_MD[%i] => Exit() \n", i ); exit(0);}
    apos [i] = p;
    vapos[i] = v;
    //fapos[i] = Vec3dZero;
    return fr2;
}

inline double move_atom_FIRE( int i, float dt, double Flim, double cv, double cf ){
    Vec3d  f   = fapos[i];
    Vec3d  v   = vapos[i];
    Vec3d  p   = apos [i];
    double fr2 = f.norm2();
    //if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    if(bPi)f.add_mul( p, -p.dot(f) ); 
    v.mul(             cv );
    v.add_mul( f, dt + cf );
    if(bPi) v.add_mul( p, -p.dot(v) );
    p.add_mul( v, dt );
    if(bPi)  p.normalize();
    apos [i] = p;
    vapos[i] = v;
    //fapos[i] = Vec3dZero;
    return fr2;
}


inline double move_atom_kvaziFIRE( int i, float dt, double Flim ){
    Vec3d  f   = fapos[i];
    Vec3d        v   = vapos[i];
    Vec3d        p   = apos [i];
    double fr2 = f.norm2();
    if(fr2>(Flim*Flim)){ f.mul(Flim/sqrt(fr2)); };
    const bool bPi = i>=natoms;
    if(bPi)f.add_mul( p, -p.dot(f) ); 
    double vv  = v.norm2();
    double ff  = f.norm2();
    double vf  = v.dot(f);
    double c   = vf/sqrt( vv*ff + 1.e-16    );
    double v_f =    sqrt( vv/( ff + 1.e-8)  );
    double cv;
    double cf = cos_damp_lin( c, cv, 0.01, -0.7,0.0 );
    v.mul(                 cv );
    v.add_mul( f, dt + v_f*cf );
    if(bPi) v.add_mul( p, -p.dot(v) );
    p.add_mul( v, dt );
    if(bPi)  p.normalize();
    apos [i] = p;
    vapos[i] = v;
    //fapos[i] = Vec3dZero;
    return fr2;
}




double move_GD(float dt, double Flim=100.0 ){
    double F2sum=0;
    for(int i=0; i<nvecs; i++){
        F2sum += move_atom_GD(i, dt, Flim);
    }
    return F2sum;
}


void makeBackNeighs( bool bCapNeighs=true ){
    for(int i=0; i<natoms; i++){ bkneighs[i]=Quat4i{-1,-1,-1,-1}; };
    for(int ia=0; ia<nnode; ia++){
        for(int j=0; j<4; j++){        // 4 neighbors
            int ja = neighs[ia].array[j];
            if( ja<0 )continue;
            //NOTE: We deliberately ignore back-neighbors from caping atoms 
            bool ret = addFirstEmpty( bkneighs[ja].array, 4, ia*4+j, -1 );
            if(!ret){ printf("ERROR in MMFFsp3_loc::makeBackNeighs(): Atom #%i has >4 back-Neighbors (while adding atom #%i) \n", ja, ia ); exit(0); }
        };
    }
    //for(int i=0; i<natoms; i++){printf( "bkneigh[%i] (%i,%i,%i,%i) \n", i, bkneighs[i].x, bkneighs[i].y, bkneighs[i].z, bkneighs[i].w );}
    //checkBkNeighCPU();
    if(bCapNeighs){   // set neighbors for capping atoms
        for(int ia=nnode; ia<natoms; ia++){ neighs[ia]=Quat4i{-1,-1,-1,-1};  neighs[ia].x = bkneighs[ia].x/4;  }
    }
}

void makeNeighCells( const Vec3i nPBC_ ){ 
    nPBC=nPBC_;
    //printf( "makeNeighCells() nPBC_(%i,%i,%i) lvec (%g,%g,%g) (%g,%g,%g) (%g,%g,%g)\n", nPBC.x,nPBC.y,nPBC.z, lvec.a.x,lvec.a.y,lvec.a.z,  lvec.b.x,lvec.b.y,lvec.b.z,   lvec.c.x,lvec.c.y,lvec.c.z );
    for(int ia=0; ia<natoms; ia++){
        Quat4i ngC = Quat4i{-1,-1,-1,-1};
        for(int j=0; j<4; j++){
            const int ja = neighs[ia].array[j];
            if( ja<0 )continue;
            const Vec3d d0 = apos[ja] - apos[ia];
            int ipbc =  0;
            int imin = -1;
            double r2min = 1e+300;
            for(int iz=-nPBC.z; iz<=nPBC.z; iz++){ for(int iy=-nPBC.y; iy<=nPBC.y; iy++){ for(int ix=-nPBC.x; ix<=nPBC.x; ix++){   
                Vec3d d = d0 + (lvec.a*ix) + (lvec.b*iy) + (lvec.c*iz); 
                double r2 = d.norm2();
                //printf( "[%i,%i][%i,%i,%i] %g   (%g,%g,%g)\n", ia,ja, ix,iy,iz, sqrt(r2), d.x,d.y,d.z );
                if(r2<r2min){   // find nearest distance
                    r2min = r2;
                    imin  = ipbc;
                }
                ipbc++; 
            }}}
            ngC.array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
        }
        //printf("\n", ngC.x,ngC.y,ngC.z,ngC.w);
        neighCell[ia]=ngC;
    }
    //printNeighs();
}

void makeNeighCells( int npbc, Vec3d* pbc_shifts ){ 
    for(int ia=0; ia<natoms; ia++){
        for(int j=0; j<4; j++){
            //printf("ngcell[%i,j=%i] \n", ia, j);
            int ja = neighs[ia].array[j];
            //printf("ngcell[%i,ja=%i] \n", ia, ja);
            if( ja<0 )continue;
            const Vec3d d = apos[ja] - apos[ia];

            // ------- Brute Force method
            int imin=-1;
            float r2min = 1.e+300;
            for( int ipbc=0; ipbc<npbc; ipbc++ ){   
                Vec3d shift = pbc_shifts[ipbc]; 
                shift.add(d);
                float r2 = shift.norm2();
                if(r2<r2min){   // find nearest distance
                    r2min=r2;
                    imin=ipbc;
                }
            }

            /*
            // -------- Fast method
            { //DEBUG
                Vec3d u;
                invLvec.dot_to( d, u );
                int ix = 1-(int)(u.x+1.5);
                int iy = 1-(int)(u.y+1.5);
                int iz = 1-(int)(u.z+1.5);
                Vec3d dpbc = lvec.a*ix + lvec.b*iy + lvec.c*iz;
                printf( "NeighCell[%i,%i] ipbc %i pbc_shift(%6.3f,%6.3f,%6.3f) dpbc(%6.3f,%6.3f,%6.3f) iabc(%2i,%2i,%2i) u(%6.3f,%6.3f,%6.3f)\n", ia, ja, imin, pbc_shifts[imin].x,pbc_shifts[imin].y,pbc_shifts[imin].z,  dpbc.x,dpbc.y,dpbc.z, ix,iy,iz, u.x,u.y,u.z  );
            }
            */

            //printf("ngcell[%i,%i] imin=%i \n", ia, ja, imin);
            neighCell[ia].array[j] = imin;
            //printf("ngcell[%i,%i] imin=%i ---- \n", ia, ja, imin);
        }
    }
}

void printSizes     (      ){ printf( "MMFFf4::printSizes(): nDOFs(%i) natoms(%i) nnode(%i) ncap(%i) nvecs(%i) npbc(%i)\n", nDOFs,natoms,nnode,ncap,nvecs,npbc ); };
void printAtomParams(int ia){ printf("atom[%i] t%i ngs{%3i,%3i,%3i,%3i} par(%5.3f,%5.3f,%5.3f,%5.3f)  bL(%5.3f,%5.3f,%5.3f,%5.3f) bK(%6.3f,%6.3f,%6.3f,%6.3f)  Ksp(%5.3f,%5.3f,%5.3f,%5.3f) Kpp(%5.3f,%5.3f,%5.3f,%5.3f) \n", ia, atypes[ia], neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,    apars[ia].x,apars[ia].y,apars[ia].z,apars[ia].w,    bLs[ia].x,bLs[ia].y,bLs[ia].z,bLs[ia].w,   bKs[ia].x,bKs[ia].y,bKs[ia].z,bKs[ia].w,     Ksp[ia].x,Ksp[ia].y,Ksp[ia].z,Ksp[ia].w,   Kpp[ia].x,Kpp[ia].y,Kpp[ia].z,Kpp[ia].w  ); };
void printNeighs    (int ia){ printf("atom[%i] neigh{%3i,%3i,%3i,%3i} neighCell{%3i,%3i,%3i,%3i} \n", ia, neighs[ia].x,neighs[ia].y,neighs[ia].z,neighs[ia].w,   neighCell[ia].x,neighCell[ia].y,neighCell[ia].z,neighCell[ia].w ); };
void printBKneighs  (int ia){ printf("atom[%i] bkngs{%3i,%3i,%3i,%3i} \n", ia, bkneighs[ia].x,bkneighs[ia].y,bkneighs[ia].z,bkneighs[ia].w ); };
void printAtomParams(      ){ printf("MMFFsp3_loc::printAtomParams()\n" ); for(int i=0; i<nnode;  i++){ printAtomParams(i); }; };
void printNeighs    (      ){ printf("MMFFsp3_loc::printNeighs()\n"     ); for(int i=0; i<natoms; i++){ printNeighs    (i);     }; };
void printBKneighs  (      ){ printf("MMFFsp3_loc::printBKneighs()\n"   ); for(int i=0; i<natoms; i++){ printBKneighs  (i);   }; };
void print_pipos    (      ){ printf("MMFFsp3_loc::print_pipos()\n"     ); for(int i=0; i<nnode;  i++){ printf( "pipos[%i](%g,%g,%g) r=%g\n", i, pipos[i].x,pipos[i].y,pipos[i].z, pipos[i].norm() ); } }
void print_apos     (      ){ printf("MMFFsp3_loc::print_apos()\n"      ); for(int i=0; i<natoms; i++){ printf( "apos [%i](%g,%g,%g)\n",      i, apos[i].x ,apos[i].y ,apos[i].z                   ); } }


void printAngles(int ia){
    Vec3d* angles_i = angles+(ia*6);
    int* ings       = neighs[ia].array; 
    int iang=0;
    for(int i=0; i<3; i++){
        int ing = ings[i];
        if(ing<0) break;
        for(int j=i+1; j<4; j++){
            int jng  = ings[j];
            if(jng<0) break;
            if(bEachAngle){
                Vec2d cs0_ss = angles_i[iang].xy();  
                double ssK    = angles_i[iang].z;
                printf( "atom[%i|%i]types{%i,%i,%i} ssK %g cs0(%g,%g) \n", ia,iang,  atypes[ing], atypes[ia], atypes[jng], ssK, cs0_ss.x, cs0_ss.y  );
                iang++; 
            }
        }
    }
}
void printAngles(      ){ printf("MMFFsp3_loc::printAngles()\n"); for(int i=0; i<nnode; i++){ printAngles(i); } }

void printTorsions(      ){ printf("MMFFsp3_loc::printTorsions()\n"); for(int i=0; i<ntors; i++){ printf( "torsion[%i]{%i,%i,%i,%i} {%i,%i,%i,%i}\n", i, tors2atom[i].x,tors2atom[i].y,tors2atom[i].z,tors2atom[i].w,   torsParams[i].x,torsParams[i].y,torsParams[i].z,torsParams[i].w  ); } }

void printAtomsConstrains( bool bWithOff=false ){ printf("MMFFsp3_loc::printAtomsConstrains()\n"); for(int i=0; i<natoms; i++){ if(bWithOff || (constr[i].w>0.0f) )printf( "consrt[%i](%g,%g,%g|K=%g)\n", i, constr[i].x,constr[i].y,constr[i].z,constr[i].w ); } }

void printDEBUG(  bool bNg=true, bool bPi=true, bool bA=true ){
    printf( "MMFFsp3_loc::printDEBUG()\n" );
    if(bA)for(int i=0; i<natoms; i++){
        printf( "CPU[%i] ", i );
        //printf( "bkngs{%2i,%2i,%2i,%2i,%2i} ",         bkNeighs[i].x, bkNeighs[i].y, bkNeighs[i].z, bkNeighs[i].w );
        printf( "fapos{%6.3f,%6.3f,%6.3f} ", fapos[i].x, fapos[i].y, fapos[i].z );
        //printf(  "avel{%6.3f,%6.3f,%6.3f} ", avel[i].x, avel[i].y, avel[i].z);
        printf(  "apos{%6.3f,%6.3f,%6.3f} ", apos[i].x, apos[i].y, apos[i].z );
        printf( "\n" );
    }
    if(bPi)for(int i=0; i<nnode; i++){
        int i1=i+natoms;
        printf( "CPU[%i] ", i1 );
        printf(  "fpipos{%6.3f,%6.3f,%6.3f} ", fapos[i1].x, fapos[i1].y, fapos[i1].z );
        //printf(  "vpipos{%6.3f,%6.3f,%6.3f} ", avel[i1].x, avel[i1].y, avel[i1].z );
        printf(   "pipos{%6.3f,%6.3f,%6.3f} ", apos[i1].x, apos[i1].y, apos[i1].z );
        printf( "\n" );
    }
    if(bNg)for(int i=0; i<nnode; i++){ for(int j=0; j<4; j++){
        int i1=i*4+j;
        //int i2=(i+natoms)*4+j;
        printf( "CPU[%i,%i] ", i, j );
        printf( "fneigh  {%6.3f,%6.3f,%6.3f} ", fneigh  [i1].x, fneigh  [i1].y, fneigh  [i1].z );
        printf( "fneighpi{%6.3f,%6.3f,%6.3f} ", fneighpi[i1].x, fneighpi[i1].y, fneighpi[i1].z );
        printf( "\n" );
    }}
}

bool checkNans( bool bExit=true, bool bNg=true, bool bPi=true, bool bA=true ){
    bool ret = false;
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*) apos,   "apos"  );
    if(bA)  ret |= ckeckNaN_d(natoms,  3, (double*)fapos,  "fapos"  );
    if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*) pipos,  "pipos" );
    if(bPi) ret |= ckeckNaN_d(nnode,   3, (double*)fpipos, "fpipos" );
    if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneigh,  "fneigh"   );
    if(bNg) ret |= ckeckNaN_d(nnode*4, 3, (double*)fneighpi,"fneighpi" );
    if(bExit&&ret){ printf("ERROR: NaNs detected in %s in %s => exit(0)\n", __FUNCTION__, __FILE__ ); exit(0); };
    return ret;
}

void rotateNodes(int n, int* sel, Vec3d p0, Vec3d ax, double phi ){
    //printf( "MMFFsp3_loc::rotateNodes() nsel=%i phi=%g ax(%g,%g,%g) p0(%g,%g,%g) \n", n,  phi, ax.x,ax.y,ax.z,   p0.x,p0.y,p0.z );
    ax.normalize();
    double ca=cos(phi);
    double sa=sin(phi);
    for(int i=0;i<n; i++){
        int ia = sel[i];
        //printf( "MMFFsp3_loc::rotateNodes() atom[%i](%g,%g,%g) cs(%g,%g) \n", ia, apos[ia].x,apos[ia].y,apos[ia].z, ca,sa );
        apos [ia].rotate_csa( ca, sa, ax, p0 );
        if(ia>=nnode)continue;
        pipos[ia].rotate_csa( ca, sa, ax     );
        int* ngs=neighs[ia].array; 
        for(int j=0;j<4;j++){
            int ja = ngs[j];
            if(ja>=0){ if(ja>nnode) apos[ ja  ].rotate_csa( ca, sa, ax, p0 ); }
        }
    }
}

void chargeToEpairs( Quat4d* REQs, int* atypes, double cQ=-0.2, int etyp=-1 ){
    for( int ia=0; ia<nnode; ia++ ){
        int* ngs=neighs[ia].array; 
        for( int j=0; j<4; j++ ){
            int ja = ngs[j]; 
            if(ja<0) continue;
            if( atypes[ja]==etyp ){ REQs[ja].z+=cQ; REQs[ia].z-=cQ; };
        }
    }
}


inline double measureCosPiPi(int ia, int ib, bool bRenorm=true){
    double c = pipos[ia].dot(pipos[ib]);
    if(bRenorm){ c/=sqrt( pipos[ia].norm2()* pipos[ib].norm2() ); }
    return c;
}
inline double measureAnglePiPi(int ia, int ib, bool bRenorm=true ){ return acos( measureCosPiPi(ia, ib, bRenorm ) ); }
inline double measureCosSigmaPi(int ipi, int ia, int ib){
    Vec3d b = apos[ib]-apos[ia];
    return pipos[ipi].dot(b)/sqrt( pipos[ipi].norm2()*b.norm2() );
}
inline double measureAngleSigmaPi(int ipi, int ia, int ib){ return acos( measureCosSigmaPi(ipi, ia, ib ) ); }

};


#endif
