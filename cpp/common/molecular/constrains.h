
#ifndef constrains_h
#define constrains_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "quaternion.h"
#include "molecular_utils.h"
#include "Forces.h"


struct DistConstr{
    Vec2i ias;
    Vec2d  ls;
    Vec2d  ks;
    Vec3d  shift;   // like PBC_shift
    double flim;
    bool active;

    DistConstr()=default;
    DistConstr( Vec2i ias_, Vec2d ls_, Vec2d ks_, double flim_=1e+300, Vec3d shift_=Vec3dZero ):ias(ias_),ls(ls_),ks(ks_),flim(flim_),shift(shift_),active(true){ };

    __attribute__((hot))  
    inline double apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec =0, Mat3d* dlvec =0 )const{
        Vec3d sh;
        if(lvec){ lvec->dot_to_T( shift, sh ); }else{ sh=shift; }
        Vec3d d   = ps[ias.b] -ps[ias.a] + sh;
        double l  = d.norm(); 
        double f,E;
        E = spring( l, ls, ks, flim, f );
        //f = ks.x*l; E= 0.5*ks.x*l*l;
        d.mul(f/l);
        fs[ias.b].sub(d);
        fs[ias.a].add(d);

        if( dlvec ){  // Forces on Lattice Vector
            //  dE/dsh = (dE/dr) * ( dr/dsh ) = f * ( dr/dsh ) = f * d|a-b+sh|/dsh = f/|a-b+sh| * sh = sh*(f/l) 
            sh.mul( f/l );
            dlvec->a.add_mul( sh, shift.a );
            dlvec->b.add_mul( sh, shift.b );
            dlvec->c.add_mul( sh, shift.c );
        }

        //fs[ias.b].add(d);
        //fs[ias.a].sub(d);
        //printf( "DistConstr:apply(%i,%i) l %g E %g f %g | ls(%g,%g) ks(%g,%g) flim %g lvec %li |sh| %g \n", ias.b, ias.a, l, E,f, ls.x,ls.y, ks.x,ks.y, flim, (long)lvec, sh.norm() );
        return E;
    }

    void print(){ printf( "bond_constr ias(%i,%i) ls(%f,%f) ks(%lf,%lf) shift(%lf,%lf,%lf) flim=%lf \n",   ias.a,ias.b,  ls.a,ls.b,   ks.a,ks.b,    shift.a,shift.b,shift.c,    flim ); };

};


struct AngleConstr{
    Vec3i  ias;     // (ia,ib,ic), ia is the center, ib and ic are the ends
    Vec2d  cs0;     // unitary complex number (cos,sin) of the equilibrium angle
    Vec3i  acell;   // indexes of PBC cell shift of bond a
    Vec3i  bcell;   // indexes of PBC cell shift of bond b
    double k;
    double flim;
    bool active;

    AngleConstr()=default;
    AngleConstr( Vec3i ias_, Vec2d cs0_, double k_ ):ias(ias_),cs0(cs0_),k(k_),active(true){ };

    __attribute__((hot))  
    inline double  apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec =0, Mat3d* dlvec =0 )const{
        Vec3d ashift,bshift;
        if(lvec){
            ashift=lvec->a*acell.a + lvec->b*acell.b + lvec->c*acell.c;
            bshift=lvec->a*bcell.a + lvec->b*bcell.b + lvec->c*bcell.c;
        }else{ ashift=Vec3dZero; bshift=Vec3dZero;}
        Vec3d f1,f2,h1,h2;
        h1.set_sub(ps[ias.x],ps[ias.y]+ashift); double ir1 = 1./h1.normalize();
        h2.set_sub(ps[ias.x],ps[ias.z]+bshift); double ir2 = 1./h2.normalize();
        double E =  evalAngleCosHalf( h1, h2, ir1,ir2, cs0, k, f1,f2 );
        //printf(  "AngleConstr::apply(%i,%i,%i) E= %g   cos(a)= %g \n", ias.x, ias.y, ias.z, E,  f1.dot(f2)/sqrt(f1.norm2()*f2.norm2()) );
        //fs[ias.x].sub(f1); fs[ias.y].add(f1);
        //fs[ias.x].sub(f2); fs[ias.z].add(f2);
        fs[ias.y].sub(f1); fs[ias.x].add(f1);
        fs[ias.z].sub(f2); fs[ias.x].add(f2);
        return E;
    }

    void print(){ printf( "angle_constr ias(%i,%i,%i) cs0(%f,%f) k(%lf) flim=%lf acell(%i,%i,%i) bcell(%i,%i,%i)\n",   ias.a,ias.b,ias.c,   cs0.x,cs0.y,      k,       flim, acell.a, acell.b, acell.c, bcell.a, bcell.b, bcell.c); };

};


struct TorsionConstr{
    Quat4i  ijkl;     // (ia,ib,ic), ia is the center, ib and ic are the ends
     Vec3d par;      
    // TODO: currently we cannot define dihedral constrains crossing PBC boundaries
    bool active;
    bool driving;
    double fDrive;
    Vec2d  fDrives; // force for liear interpolation



    TorsionConstr()=default;
    TorsionConstr( Quat4i ijkl_, Vec3d par_, double k_ ):ijkl(ijkl_),par(par_),active(true),driving(false),fDrive(0.0){ };

    inline void   push_fDrive(double f){    fDrives.x=fDrives.y; fDrives.y=f; }
    inline void   generate   (){ push_fDrive( randf( -fDrive, fDrive ) ); }
    inline double interpolate_fDrive(double t)const{ return fDrives.x*(1-t) + fDrives.y*t; }

    __attribute__((hot))  
    inline double  apply( Vec3d* apos, Vec3d* fapos, double t, Mat3d* lvec=0, Mat3d* dlvec =0 )const{
        const Vec3d p2  = apos[ijkl.y];
        const Vec3d p3  = apos[ijkl.z];
        const Vec3d r32 = p3-p2;
        const Vec3d r12 = apos[ijkl.x]-p2; 
        const Vec3d r43 = apos[ijkl.w]-p3;

        Vec3d n123; n123.set_cross( r12, r32 );  //  |n123| = |r12| |r32| * sin( r12, r32 ) 
        Vec3d n234; n234.set_cross( r43, r32 );  //  |n234| = |r43| |r32| * sin( r43, r32 )
        
        // ===== Prokop's
        const double l32     = r32   .norm ();   // we can avoid this sqrt() if we read it from hneigh
        const double il2_123 = 1/n123.norm2();   //  il2_123 =  1/ (   |r12| |r32| * sin( r12, r32 ) )^2 
        const double il2_234 = 1/n234.norm2();   //  il2_234 =  1/ (   |r43| |r32| * sin( r43, r32 ) )^2
        const double inv_n12 = sqrt(il2_123*il2_234); // inv_n12 = 1/ ( |r12| |r32| * sin( r12, r32 ) * |r43| |r32| * sin( r43, r32 ) )
        // --- Energy
        const Vec2d cs{
             n123.dot(n234)*inv_n12     ,
            -n123.dot(r43 )*inv_n12*l32
        };
        Vec2d csn = cs;
        const int n = (int)par.z;
        for(int i=1; i<n; i++){ csn.mul_cmplx(cs); }
        
        double E  =  par.x * ( 1.0 + par.y * csn.x );
        // --- Force on end atoms
        double f;
        if(driving){ //double f = fDrive;
            f = interpolate_fDrive(t);
            //printf( "TorsionConstr::apply[%i] driving f=%g fDrives(%g,%g) t=%g fDrive %g \n", (int)ijkl.x, f, fDrives.x, fDrives.y, t, fDrive );
        } else  { f = -par.x * par.y * par.z * csn.y; }

        f*=l32;
        Vec3d fp1; fp1.set_mul(n123,-f*il2_123 );
        Vec3d fp4; fp4.set_mul(n234, f*il2_234 );
        // --- Recoil forces on axis atoms
        double il2_32 = -1/(l32*l32);
        double c123   = r32.dot(r12)*il2_32;
        double c432   = r32.dot(r43)*il2_32;
        Vec3d fp3; fp3.set_lincomb(  c123,   fp1,  c432-1., fp4 );   // from condition torq_p2=0  ( conservation of angular momentum )
        Vec3d fp2; fp2.set_lincomb( -c123-1, fp1, -c432   , fp4 );   // from condition torq_p3=0  ( conservation of angular momentum )
        //Vec3d fp2_ = (fp1_ + fp4_ + fp3_ )*-1.0;                   // from condition ftot=0     ( conservation of linear  momentum )
        //Vec3d fp3_ = (fp1_ + fp4_ + fp2_ )*-1.0;                   // from condition ftot=0     ( conservation of linear  momentum )
        // fapos[ijkl.x]=fp1;
        // fapos[ijkl.y]=fp2;
        // fapos[ijkl.z]=fp3;
        // fapos[ijkl.w]=fp4;
        {
            const Vec3d p1 = apos[ijkl.x]; 
            const Vec3d p4 = apos[ijkl.w];
            glColor3f(0.0,0.0,1.0);
            Draw3D::drawVecInPos( n123, p2 );
            Draw3D::drawVecInPos( n234, p3 );
            glColor3f(1.0,0.0,0.0);
            Draw3D::drawVecInPos( fp1, p1 );
            Draw3D::drawVecInPos( fp4, p4 );
            Draw3D::drawVecInPos( fp2, p2 );
            Draw3D::drawVecInPos( fp3, p3 );
            //Draw3D::drawVecInPos( fp1, p2, 0.1, 0xFF0000 );

        }
        return E;
    }

    void print(){ printf( "TorsionConstr ijkl(%i,%i,%i,%i) par(c0=%g,k=%g,|n=%f) fDrive=%g driving=%i \n",   ijkl.x,ijkl.y,ijkl.z,ijkl.w,   par.x,par.y,par.z,  fDrive, driving ); };

};



class Constrains{ public:
    std::vector<DistConstr>  bonds;
    std::vector<AngleConstr> angles;
    std::vector<TorsionConstr> torsions;

    double dt = 0.01;
    int nDriveUpdate = 100;
    int iDriveUpdate = 0;

    int addBond( const char* line, int* atom_permut=0, int _0=1 ){
        int i = bonds.size();
        DistConstr cons; cons.active=true;
        int nret = sscanf( line, "b %i %i   %lf %lf   %lf %lf   %lf    %lf %lf %lf",   &cons.ias.a,&cons.ias.b,  &cons.ls.a,&cons.ls.b,  &cons.ks.a,&cons.ks.b,  &cons.flim, &cons.shift.a,&cons.shift.b,&cons.shift.c );
        cons.ias.a-=_0;
        cons.ias.b-=_0;
        if( atom_permut ){
            printf( "permut bond (%i->%i)-(%i->%i) \n", cons.ias.a, atom_permut[cons.ias.a], cons.ias.b, atom_permut[cons.ias.b] ); 
            cons.ias.a=atom_permut[cons.ias.a];
            cons.ias.b=atom_permut[cons.ias.b];
        }
        if(nret<10){cons.shift=Vec3dZero;}
        if(nret<7 ){ printf("WARRNING : Constrains::loadBonds[%i] bond nret(%i)<10 line=%s", i, nret, line ); }
        cons.print();
        bonds.push_back( cons );
        return i;
    }

    int addAngle( const char* line, int* atom_permut=0, int _0=1 ){
        int i = angles.size();
        double ang;
        AngleConstr cons; cons.active=true;
        int nret = sscanf( line, "g %i %i %i   %lf  %lf %lf   %i %i %i  %i %i %i",    &cons.ias.a,&cons.ias.b,&cons.ias.c,   &ang,   &cons.k,   &cons.flim,   &cons.acell.a,&cons.acell.b,&cons.acell.c,   &cons.bcell.a,&cons.bcell.b,&cons.bcell.c );
        cons.cs0.fromAngle( (ang/180.0)*M_PI*0.5 );
        cons.ias.a-=_0;
        cons.ias.b-=_0;
        cons.ias.c-=_0;
        if( atom_permut ){
            printf( "permut angle (%i->%i)-(%i->%i)-(%i->%i) \n", cons.ias.b, atom_permut[cons.ias.b],   cons.ias.a, atom_permut[cons.ias.a],    cons.ias.c, atom_permut[cons.ias.c] ); 
            cons.ias.a=atom_permut[cons.ias.a];
            cons.ias.b=atom_permut[cons.ias.b];
            cons.ias.c=atom_permut[cons.ias.c];
        }
        if(nret<12){cons.bcell=Vec3iZero;}
        if(nret<9 ){cons.acell=Vec3iZero;}
        if(nret<6){ printf("WARRNING : Constrains::loadBonds[%i] angle nret(%i)<6 line=%s", i, nret, line ); }
        cons  .print();
        angles.push_back( cons );
        return i;
    }

    int addTorsion( const char* line, int* atom_permut=0, int _0=1 ){
        int i = torsions.size();
        Quat4i ijkl;
        Vec3d par;
        double fDrive;
        int nret = sscanf( line, "t %i %i %i %i   %lf %lf %lf    %lf ",    &ijkl.x,&ijkl.y,&ijkl.z,&ijkl.w,   &par.x,&par.y,&par.z,  &fDrive );
        ijkl.x-=_0;
        ijkl.y-=_0;
        ijkl.z-=_0;
        ijkl.w-=_0;
        if( atom_permut ){
            printf( "permut torsion (%i->%i)-(%i->%i)-(%i->%i)-(%i->%i) \n", ijkl.x, atom_permut[ijkl.x],   ijkl.y, atom_permut[ijkl.y],    ijkl.z, atom_permut[ijkl.z],    ijkl.w, atom_permut[ijkl.w] ); 
            ijkl.x=atom_permut[ijkl.x];
            ijkl.y=atom_permut[ijkl.y];
            ijkl.z=atom_permut[ijkl.z];
            ijkl.w=atom_permut[ijkl.w];
        }
        if(nret<8){printf("WARRNING : Constrains::loadBonds[%i] torsion nret(%i)<8 line=%s", i, nret, line ); }
        TorsionConstr cons; cons.active=true;
        cons.ijkl=ijkl;
        cons.par=par;
        if(nret<7){par.z=1;}
        if(nret<8){cons.fDrive=0;}
        if( fabs(fDrive)>1e-16 ){ cons.driving=true; }else{ cons.driving=false; }
        cons.fDrive=fDrive;
        cons.print();
        torsions.push_back( cons );
        return i;
    }

    __attribute__((hot))  
    double apply( Vec3d* ps, Vec3d* fs, Mat3d* lvec=0, Mat3d* dlvec=0 ){
        double E=0;  
        int i=0;
        if(iDriveUpdate>=nDriveUpdate){ iDriveUpdate=0; }
        double t = iDriveUpdate/(double)nDriveUpdate;
        for( const DistConstr&  c : bonds  ){ E+= c.apply(ps,fs, lvec, dlvec ); }
        for( const AngleConstr& c : angles ){ E+= c.apply(ps,fs, lvec, dlvec ); }
        i=0;
        for( TorsionConstr& c : torsions ){ 
            if(iDriveUpdate==0){ c.generate(); }
            E+= c.apply(ps,fs, t, lvec, dlvec ); 
            printf( "TorsionConstr::apply[%i] iDriveUpdate %i t %g \n", i, iDriveUpdate, t );
            i++;
        }
        iDriveUpdate++;
        return E;
    }

    int loadBonds( const char* fname, int* atom_permut=0, int _0=1 ){
        //printf("Constrains::loadBonds(%s) atom_permut=%li _0=%i\n", fname, (long)atom_permut, _0 );
        FILE* pFile = fopen( fname, "r" );
        if(pFile==0){ printf("ERROR in Constrains::loadBonds(%s) - No Such File \n", fname ); return -1; }
        else{
            char buff[1024];            
            int i=0;
            for(i=0; i<10000; i++){
                char* line = fgets( buff, 1024, pFile );
                printf( "Constrains::loadBonds[i=%4i] %s", i, line );
                if(line==NULL)  break;
                if     (line[0]=='#'){continue;}
                else if(line[0]=='b'){
                    addBond( line, atom_permut, _0=_0 );
                }else if(line[0]=='g'){
                    addAngle( line, atom_permut, _0=_0 );
                }else if(line[0]=='t'){
                    addTorsion( line, atom_permut, _0=_0 );
                }
            }
            return i;
            fclose( pFile );
        }
    }

    void clear( bool bShring=false){
        bonds.clear();
        angles.clear();
        if(bShring){
            bonds.shrink_to_fit();
            angles.shrink_to_fit();
        }
    }
    
};

#endif
