#ifndef MolGUI_tests_h
#define MolGUI_tests_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Forces.h"

#include <SDL2/SDL.h>

#include "Draw3D.h"

double torsion_Paolo( Vec3d p1, Vec3d p2, Vec3d p3, Vec3d p4, Vec3d par ){

    // Vec3d  r12, r32;
    // double l12, l32;
    // for(int in=0; in<4; in++){
    //     int ing = ingsj[in];
    //     if(ing<0) { break; }
    //     if     (ing==i) { r12 = hneigh[j*4+in].f; l12 = 1.0/hneigh[j*4+in].e; }   
    //     else if(ing==k) { r32 = hneigh[j*4+in].f; l32 = 1.0/hneigh[j*4+in].e; } 
    // }
    // Vec3d r43;
    // double l43;
    // for(int in=0; in<4; in++){
    //     int ing = ingsk[in];
    //     if(ing<0) { break; }
    //     if     (ing==l) { r43 = hneigh[k*4+in].f; l43 = 1.0/hneigh[k*4+in].e; }   
    // }
    // Vec3d r12abs; r12abs.set_mul( r12, l12 );
    // Vec3d r32abs; r32abs.set_mul( r32, l32 );
    // Vec3d r43abs; r43abs.set_mul( r43, l43 );
    // Vec3d n123; n123.set_cross( r12abs, r32abs );
    // Vec3d n234; n234.set_cross( r43abs, r32abs );

    Vec3d r12 = p1-p2; 
    Vec3d r32 = p3-p2;
    Vec3d r43 = p4-p3;
    Vec3d n123; n123.set_cross( r12, r32 );
    Vec3d n234; n234.set_cross( r43, r32 );

    double l123 = n123.normalize();
    double l234 = n234.normalize();
    double l32  = r32.normalize();
    double l12  = r12.normalize();
    double l43  = r43.normalize();
    double cos = n123.dot(n234);
    //double cos = abs(n123.dot(n234))/(l123*l234);
    double sin = sqrt(1.0-cos*cos+1e-14); // must be positive number !!!
    int n = (int)par.z;
    Vec2d cs{cos,sin};
    Vec2d csn{cos,sin};
    for(int i=1; i<n; i++){
        csn.mul_cmplx(cs);
    }
    double E = par.x * ( 1.0 + par.y * csn.x );
    Vec3d scaled_123; scaled_123.set_mul( n123, cos );
    Vec3d scaled_234; scaled_234.set_mul( n234, cos );
    Vec3d tmp_123; tmp_123.set_sub( n123, scaled_234 );
    Vec3d tmp_234; tmp_234.set_sub( n234, scaled_123 );
    Vec3d f_12; f_12.set_cross( r32, tmp_234 );
    Vec3d f_43; f_43.set_cross( r32, tmp_123 );
    Vec3d tmp1_32; tmp1_32.set_mul( tmp_234, l12/l123 );
    Vec3d tmp2_32; tmp2_32.set_mul( tmp_123, l43/l234 );
    Vec3d vec1_32; vec1_32.set_cross( tmp1_32, r12 );
    Vec3d vec2_32; vec2_32.set_cross( tmp2_32, r43 );
    Vec3d vec_32; vec_32.set_add( vec1_32, vec2_32 );
    double fact = -par.x * par.y * par.z * csn.y / sin ;
    Vec3d f_32; f_32.set_mul(vec_32,fact);

    // ---- Final Force on vectors
    Vec3d fp1  = f_12 * (fact*l32/l123 );
    Vec3d fp4  = f_43 * (fact*l32/l234 );
    Vec3d fp2  = ( fp1 + f_32 )*-1.;
    Vec3d fp3  = (  fp4 - f_32)*-1.; 
    //fdih[id*4  ].set_mul(f_12,fact*l32/l123);
    //fdih[id*4+3].set_mul(f_43,fact*l32/l234);
    //fdih[id*4+1].set_add(fdih[id*4],f_32);
    //fdih[id*4+1].set_mul(fdih[id*4+1],-1.0);
    //fdih[id*4+2].set_sub(f_32,fdih[id*4+3]);    

    {  // CHECK
        //printf(  "|f1|/|f1_| %g  |f3|/|f3_| %g cos(f1,f1_) %g cos(f3,f3_) %g |cs| %g \n", f1.norm()/f1_.norm(), f3.norm2()/f3_.norm2(),  f1.dot(f1_)/sqrt(f1.norm2()*f1_.norm2()),  f3.dot(f3_)/sqrt(f3.norm2()*f3_.norm2()),  c*c+s*s );

        opengl1renderer.lineWidth(2.0);
        Draw3D::drawVecInPos( fp1, p1, COLOR_RED );
        Draw3D::drawVecInPos( fp2, p2, COLOR_RED );
        Draw3D::drawVecInPos( fp3, p3, COLOR_RED );
        Draw3D::drawVecInPos( fp4, p4, COLOR_RED );
        Vec3d ax = r32; ax.normalize();
        Vec3d ftot=Vec3dZero;
        Vec3d torq=Vec3dZero;
        
        torq.set_cross( p1-p2, fp1 );
        torq.add_cross( p4-p2, fp4 );
        //torq.add_cross( p3-p2, (f_32 + f1)*-1 );    
        torq.add_cross( p3-p2, fp3 );   

        //Vec3d fp3_ = ( fp1*ax.dot(p1-p2) + fp4*ax.dot(p4-p2) )*( -1./ax.dot(p3-p2) ); 
        Vec3d fp3_; fp3_.set_lincomb( ax.dot(p1-p2)/-l32, fp1, ax.dot(p4-p2)/-l32, fp4 );
        Vec3d fp2_ = (fp3_ + fp1 + fp4)*-1.0; 

        Draw3D::drawVecInPos( fp3_, p3, {1, 0, 1} );
        Draw3D::drawVecInPos( fp2_, p2, {1, 0, 1} );

        Vec3d ftot2=Vec3dZero;
        Vec3d torq2=Vec3dZero;
        //torq2.set_cross( p1-p3, fp1 );
        //torq2.add_cross( p4-p3, fp4 );
        //torq2.add_cross( p2-p3, fp2 );    
        //torq2.add_cross( p3-p3, fp3 );   

        torq2.set_cross( p1-p2, fp1 );
        torq2.add_cross( p4-p2, fp4 );
        //torq2.add_cross( p3-p2, (f_32 + f1)*-1 );    
        torq2.add_cross( p3-p2, fp3_ );

        ftot2 = fp1 + fp2_ + fp3_ + fp4;
        //printf( "Paolo |ftot|=%g |torq|=%g |torq2|=%g  <torq|ax>=%g \n", ftot.norm(), torq.norm(), torq2.norm(), ax.dot(torq) );

        printf( "Paolo |ftot|=%g |ftot2|=%g |torq|=%g |torq2|=%g  \n", ftot.norm(), ftot2.norm(), torq.norm(), torq2.norm() );

    }
    return E;
}


Vec3d zeroTorqCondition( const Vec3d& d32, double invL2, const Vec3d& d12, const Vec3d& d43, const Vec3d& f1, const Vec3d& f4 ){
    //double il2 = -1/(l32*l32);
    Vec3d f3; f3.set_lincomb( d32.dot(d12)*invL2, f1, d32.dot(d43)*invL2 - 1., f4 ); 
    return f3;
}


double torsion_Prokop( Vec3d p1, Vec3d p2, Vec3d p3, Vec3d p4, Vec3d par ){

    Vec2d cs,csn;
    int   n = (int)par.z;

    // ===== Common
    Vec3d r12 = p1-p2; 
    Vec3d r32 = p3-p2;
    Vec3d r43 = p4-p3;
    Vec3d n123; n123.set_cross( r12, r32 );  //    |n123| = |r12| |r32| * sin( r12, r32 ) 
    Vec3d n234; n234.set_cross( r43, r32 );  //    |n234| = |r43| |r32| * sin( r43, r32 )
    
    // ===== Prokop's
    double l32     = r32   .norm();
    double il2_123 = 1/n123.norm2();
    double il2_234 = 1/n234.norm2();
    double inv_n12 = sqrt(il2_123*il2_234);
    // --- 
    cs.x =  n123.dot(n234)*inv_n12;
    cs.y = -n123.dot(r43 )*inv_n12*l32;
    double c = cs.x;
    double s = cs.y;
    csn = cs;
    for(int i=1; i<n; i++){ csn.mul_cmplx(cs); }
    // --- Force on end atoms
    double f = -par.x * par.y * par.z * csn.y; 
    f*=l32;
    Vec3d fp1_; fp1_.set_mul(n123,-f*il2_123 );
    Vec3d fp4_; fp4_.set_mul(n234, f*il2_234 );
    // --- Recoil forces on axis atoms
    double il2_32 = -1/(l32*l32);
    double c123   = r32.dot(r12)*il2_32;
    double c432   = r32.dot(r43)*il2_32;
    Vec3d fp3_; fp3_.set_lincomb(  c123,   fp1_,  c432-1., fp4_ );   // from condition torq_p2=0  ( conservation of angular momentum )
    Vec3d fp2_; fp2_.set_lincomb( -c123-1, fp1_, -c432   , fp4_ );   // from condition torq_p3=0  ( conservation of angular momentum )
    //Vec3d fp2_ = (fp1_ + fp4_ + fp3_ )*-1.0;                       // from condition ftot=0     ( conservation of linear  momentum )
    //Vec3d fp3_ = (fp1_ + fp4_ + fp2_ )*-1.0;                       // from condition ftot=0     ( conservation of linear  momentum )

    // ===== Paolo's
    double l123 = n123.normalize();
    double l234 = n234.normalize();
    double l32_ = r32.normalize();
    double l12  = r12.normalize();
    double l43  = r43.normalize();
    // --- Energy
    double cos = n123.dot(n234);
    double sin = sqrt(1.0-cos*cos+1e-14); // must be positive number !!!
    cs = Vec2d{cos,sin};
    csn=cs;
    for(int i=1; i<n; i++){ csn.mul_cmplx(cs); }
    double E = par.x * ( 1.0 + par.y * csn.x );
    double fact = -par.x * par.y * par.z * csn.y / sin ;
    // --- Forces of end atoms
    Vec3d scaled_123; scaled_123.set_mul( n123, cos );
    Vec3d scaled_234; scaled_234.set_mul( n234, cos );
    Vec3d tmp_123; tmp_123.set_sub( n123, scaled_234 );
    Vec3d tmp_234; tmp_234.set_sub( n234, scaled_123 );
    Vec3d f_12; f_12.set_cross( r32, tmp_234 );
    Vec3d f_43; f_43.set_cross( r32, tmp_123 );
    // --- Recoil forces on axis atoms
    Vec3d tmp1_32; tmp1_32.set_mul( tmp_234, l12/l123 );
    Vec3d tmp2_32; tmp2_32.set_mul( tmp_123, l43/l234 );
    Vec3d vec1_32; vec1_32.set_cross( tmp1_32, r12 );
    Vec3d vec2_32; vec2_32.set_cross( tmp2_32, r43 );
    Vec3d vec_32; vec_32.set_add( vec1_32, vec2_32 );
    Vec3d f_32; f_32.set_mul(vec_32,fact);
    // ---- Final Force on vectors
    Vec3d fp1  = f_12 * (fact*l32/l123 );
    Vec3d fp4  = f_43 * (fact*l32/l234 );
    Vec3d fp2  = ( fp1 + f_32 )*-1.;
    Vec3d fp3  = ( fp4 - f_32 )*-1.; 

    {  
        Vec3d ftot,tq_p2,tq_p3;
        ftot = fp1 + fp2 + fp3 + fp4;
        tq_p2=Vec3dZero;
        tq_p2.set_cross( p1-p2, fp1 );
        tq_p2.add_cross( p4-p2, fp4 ); 
        tq_p2.add_cross( p3-p2, fp3 ); 
        tq_p3=Vec3dZero;
        tq_p3.set_cross( p1-p3, fp1 );
        tq_p3.add_cross( p4-p3, fp4 );
        tq_p3.add_cross( p2-p3, fp2 );   
        printf( "Paolo |ftot|=%g |torq_p2|=%g |torq_p3|=%g \n", ftot.norm(), tq_p2.norm(), tq_p3.norm() );     

        // plot Paolo's forces
        opengl1renderer.lineWidth(1.0);
        opengl1renderer.color3f(1.0,0.0,0.0); 
        // Draw3D::drawVecInPos( fp1, p1 );
        // Draw3D::drawVecInPos( fp2, p2 );
        // Draw3D::drawVecInPos( fp3, p3 );
        // Draw3D::drawVecInPos( fp4, p4 );
        Draw3D::drawArrow( p1, p1+fp1_, 0.1 );
        Draw3D::drawArrow( p2, p2+fp2_, 0.1 );
        Draw3D::drawArrow( p3, p3+fp3_, 0.1 );
        Draw3D::drawArrow( p4, p4+fp4_, 0.1 );

        // plot Prokop's forces
        opengl1renderer.lineWidth(2.0);
        Draw3D::drawVecInPos( fp1, p1, {1, 0, 1} );
        Draw3D::drawVecInPos( fp2, p2, {1, 0, 1} );
        Draw3D::drawVecInPos( fp3, p3, {1, 0, 1} );
        Draw3D::drawVecInPos( fp4, p4, {1, 0, 1} );
        opengl1renderer.lineWidth(1.0);

        ftot = fp1_ + fp2_ + fp3_ + fp4_;
        tq_p2=Vec3dZero;
        tq_p2.set_cross( p1-p2, fp1_ );
        tq_p2.add_cross( p4-p2, fp4_ ); 
        tq_p2.add_cross( p3-p2, fp3_ ); 
        tq_p3=Vec3dZero;
        tq_p3.set_cross( p1-p3, fp1_ );
        tq_p3.add_cross( p4-p3, fp4_ );
        tq_p3.add_cross( p2-p3, fp2_ );   
        printf( "Prokop |ftot|=%g |torq_p2|=%g |torq_p3|=%g \n", ftot.norm(), tq_p2.norm(), tq_p3.norm() );    

    }

    return E;
}

#endif