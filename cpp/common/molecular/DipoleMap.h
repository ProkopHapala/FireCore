
#ifndef DipoleMap_h
#define DipoleMap_h

#include "NBFF.h"
#include "MMFFparams.h"

class DipoleMap { public:
    // https://docs.lammps.org/Howto_tip5p.html
    // TIP5P - Transferable Intermolecular Potential 5 Point   - 5 particles   (  O(0.0e), H1(+0.241e), H1(+0.241), E1(-0.241e), E2(-0.241e) ), E ~ 0.7A from Oxygen
    
    const MMFFparams* params=0;
    const NBFF* ff=0;
    const char* name1="H";
    const char* name2="F";
    Quat4d REQ { 1.187, sqrt(0.0006808    ), +0.24, 0.0 };
    Quat4d REQ2{ 1.750, sqrt(0.00260184625), -0.24, 0.0 };
    double Lbond = 1.0;
    std::vector<Vec3d> particles;
    std::vector<Vec3d> particles2;
    std::vector<int>   pivots; // index of atom to which the particle is attached

    int nphi = 0;
    std::vector<Quat4d> FE;;
    std::vector<Quat4d> FE2;

    void distribute0( int nphi, int nr, double r0, double dr, double z0 ){
        // ----------- SELECT E-PAIRS and add hydrogen particles next to it
        particles.clear();
        particles2.clear();
        double dphi = 2*M_PI/nphi;
        for(int ix=0; ix<nr; ix++){
            for(int iy=0; iy<nphi; iy++){
                double r = r0 + dr*ix;
                double phi = dphi*iy;
                Vec3d p = Vec3d{ r*cos(phi), r*sin(phi), z0 };
                Vec3d d = p; d.normalize();
                particles.push_back( p );
                particles2.push_back( p + d*Lbond );

            }
        } 
    }

    double findNearestAtom( Vec3d p, int& imin, double R0=0.0 ){
        imin=-1;
        double r2min  = 1e+300;
        //double r2min2 = 1e+300;
        for(int j=0; j<ff->natoms; j++){
            //printf( "findNearestAtom[%i] typ=%i REQ(R=%g,Q=%g)\n", j, ff->atypes[j], ff->REQs[j].x,ff->REQs[j].z, ff->apos[j].x, ff->apos[j].y, ff->apos[j].z );
            Vec3d d = p - ff->apos[j];
            double R = ff->REQs[j].x + R0;
            double r2 = d.norm2() - R*R;
            if(r2<r2min){ r2min=r2; imin=j; }
            //if(imin2){
            //    if((r2<r2min2)&&(r2>r2min)){ r2min2=r2; *imin2=j; }
            // }
        }
        return r2min;
    }

    void orient( double cCog, double cEfield, double cPivot, Vec3d cog=Vec3dZero, bool bClampField=false ){
        particles2.clear();
        bool bEfield = (cEfield>0);  if(!bEfield) cEfield = 0;
        bool bPivot  = (cPivot >0);  if(!bPivot ) cPivot  = 0;
        bool bCog    = (cCog   >0);  if(!bCog   ) cCog    = 0;
        for(int i=0; i<particles.size(); i++){
            //printf( "orient[%i] bEfield=%i bPivot=%i bCog=%i \n", i, bEfield, bPivot, bCog );
            Vec3d p  = particles[i];
            //Vec3d p2 = particles2[i];
            Vec3d d = Vec3dZero;
            double csum = 0;
            if(bCog){
                Vec3d di = p - cog; 
                di.normalize();
                d.add_mul( di, cCog );
                csum += cCog;
            }
            if(bEfield){
                Quat4d fe = ff->evalLJQs( p, { 1.0, 0.0, 1.0 }, 1.0 );
                
                double c = cEfield;
                if(bClampField){ c *= fmax( 0.0, -fe.e ); }

                fe.f.normalize();
                d.add_mul( fe.f , c );
                csum += c;
                printf( "orient[%i] fe.e=%g cEfield=%g c=%g cEfield=%g cCog=%g csum=%g \n", i, fe.e, c, cEfield, cCog, csum  );
            }
            if(bPivot){
                int    imin=-1;
                findNearestAtom( p, imin );
                Vec3d d_ = (p - ff->apos[imin]); 
                d.normalize();
                d.add_mul( d_, cPivot);
                pivots.push_back( imin );
                csum += cPivot;
            }
            d.normalize();
            particles2.push_back( p + d*Lbond );
        }
    }

    void orient_cos( Vec2d cs12=Vec2d{0.0,0.5}, Vec3d cog=Vec3dZero ){
        // don't forget to tur off electron pairs before calling this function
        particles2.clear();
        for(int i=0; i<particles.size(); i++){
            //printf( "orient[%i] bEfield=%i bPivot=%i bCog=%i \n", i, bEfield, bPivot, bCog );
            Vec3d p  = particles[i]; 
            Vec3d d   = p - cog;           
            int    imin=-1;
            findNearestAtom( p, imin );
            //printf(  "orient_cos[%i] imin=%i REQ(R=%g,Q=%g) apos(%g,%g,%g)\n", i, imin, ff->REQs[imin].x,ff->REQs[imin].z );
            if( ff->REQs[imin].z < -0.1 ){ // only eletro-negative atoms can be pivot
                REQ2 = ff->REQs[imin];
                Vec3d piv = ff->apos[imin];
                Vec3d a   = cog - piv;
                Vec3d b   = p   - piv;
                double c = -a.dot(b)/sqrt(a.norm2()*b.norm2());
                //printf(  "orient_cos[%i] imin=%i cos=%g\n", i, imin, c );
                if( c>cs12.x ){    
                    //printf(  "orient_cos[%i] imin=%i cos(%g) > cs12.x(%g) \n", i, imin, c, cs12.x );
                    if( c>cs12.y ){ // take pivot as a reference
                        d = b;
                    }else{  // linear interpolation based on cos of angle
                        double cb = (c-cs12.x)/(cs12.y-cs12.x);
                        double cd = 1-cb;
                        printf(  "orient_cos[%i] imin=%i REQ(R=%g,Q=%g) cos=%g lerp(cd=%g,cb=%g)\n", i, imin, ff->REQs[imin].x,ff->REQs[imin].z, c, cd, cb );
                        //d = d*( f/d.norm() ) + b*( (1-f)/b.norm() );
                        d = d*( cd/d.norm() ) + b*( cb/b.norm() );
                        //d.mul( -1.0 );
                    }
                }
            }
            d.normalize();
            particles2.push_back( p + d*Lbond );
        }
    }


    void orient_Rgrad( double w=0.1, double cut=-0.8 ){
        // don't forget to tur off electron pairs before calling this function
        particles2.clear();
        for(int i=0; i<particles.size(); i++){
            //printf( "orient[%i] bEfield=%i bPivot=%i bCog=%i \n", i, bEfield, bPivot, bCog );
            Vec3d p  = particles[i];        
            Vec3d d = Vec3dZero;
            double c = 0;
            for(int j=0; j<ff->natoms; j++){
                Vec3d d_ = p - ff->apos[j];
                double r2 = d_.norm();
                double ci = 1/( fmax( 0.0, r2/( sq(ff->REQs[j].x) - ff->REQs[j].z ) + cut  ) + w*w );
                //double ci = 1/( fmax( 0.0, r2 + w2 ) + 0.0001 );
                //printf( "orient_2piv[%i,%i] r2(%g) ci(%g) r2/R2(%g) REQ(R=%g)\n", i, j, sqrt(r2), ci,  r2/sq(ff->REQs[j].x), ff->REQs[j].x );
                d.add_mul( d_, ci );
                c += ci;
            }
            //printf( "orient_2piv[%i] ci(%g) \n", i, c );
            d.mul( 1/c );
            d.normalize();
            particles2.push_back( p + d*Lbond );
        }
    }

    void relax_orientation( double dt, double Fconv, int niter ){
        for(int i=0; i<particles.size(); i++){
            Vec3d p  = particles [i];
            Vec3d p2 = particles2[i];
            Vec3d v = Vec3dZero;
            Vec3d d = p2 - p;
            d.normalize();
            for(int j=0; j<niter; j++){
                Vec3d f   = Vec3dZero;
                Quat4d fe = ff->evalLJQs( p2, REQ2, 1.0 );

                double cvf = fe.f.dot( v ); if(cvf<0) v=Vec3dZero; // stop if force and velocity are in opposite directions

                double cf = d.dot( fe.f ); fe.f.add_mul( d, -cf );  // remove component of force parallel to d

                v .add_mul( fe.f, dt );

                double cv = d.dot( fe.f ); fe.f.add_mul( d, -cv );    // remove component of velocity parallel to d

                // --- project to sphere with radius Lbond
                p2.add_mul( v, dt );
                Vec3d d2 = p2 - p;
                d2.mul( Lbond / d2.norm() );
                p2 = p + d2;
                
                if( f.norm2() < Fconv ) break;
            }
            particles2[i] = p + d;
        }
    }

    void points_along_rim( double R, Vec3d p0, Vec2d d0=Vec2d{0.0,0.1}, int niter=10000, double Fconv=1e-4, double dt=0.05, double K=1.0, int npmax = 1000, bool bStopByPhi=false ){
        printf( "DipoleMap::points_along_rim() R %g p0 (%g,%g,%g) d0 (%g,%g) niter %i Fconv %g dt %g K %g npmax %i \n", R, p0.x,p0.y,p0.z, d0.x,d0.y, niter, Fconv, dt, K, npmax );
        double F2conv = Fconv*Fconv;
        double dr = d0.normalize();
        Vec3d  d  = Vec3d{d0.x,d0.y,0.0};
        Vec3d  p  = p0;
        Vec3d  po = p0+(d*-dr);

        Vec2d rot0{p.x,p.y}; rot0.normalize();
        double rotSign; 
        // if( bStopByPhi ){
        //     Vec2d rot{po.x,po.y}; rot.normalize();
        //     rot.udiv_cmplx( rot0 );
        //     rotSign = (rot.y>0)?1.0:-1.0;
        // }
        //printf( "points_along_rim()[i] d(%g,%g) p(%g,%g) po(%g,%g) p0(%g,%g) \n", -1, d.x,d.y, p.x,p.y, po.x,po.y, p0.x,p0.y );
        for(int i=0; i<npmax; i++ ){
            //printf( "points_along_rim()[%i] (%g,%g,%g) \n", i, p.x,p.y,p.z );
            Vec3d v = Vec3dZero;
            for( int j=0; j<niter; j++ ){  // relaxation of particle along vdW-rim in distance dr from p0
                Quat4d fe  = ff->evalLJQs( p, Quat4d{R,-1.0,0.0,0.0}, 1.0 );
                // // -- constrain force
                double c = fe.f.dot(d); fe.f.add_mul( d, -c );  // remove component of force parallel to d
                fe.f.z = 0;
                fe.f.add_mul( d, (d.dot(p-po)-dr) *-K );        // add spring force along d
                if( fe.f.norm2() < F2conv ) break;
                // -- move
                double cvf = fe.f.dot( v ); if(cvf<0) v=Vec3dZero; // stop if force and velocity are in opposite directions
                v.add_mul( fe.f, dt );
                v.z=0;

                // if( i==1 ){
                //     //printf( "points_along_rim()[%i,%i] p(%g,%g,%g) fe.f(%g,%g,%g) v(%g,%g,%g) \n", i,j, p.x,p.y,p.z, fe.f.x,fe.f.y,fe.f.z, v.x,v.y,v.z );
                //     //printf( "points_along_rim()[%i,%i] p(%g,%g) p0(%g,%g) d(%g,%g) fe.f(%g,%g,%g) v(%g,%g,%g) \n", i,j, p.x,p.y, po.x,po.y,  d.x,d.y, fe.f.x,fe.f.y, v.x,v.y );
                //     //printf( "points_along_rim()[%i,%i] d(p,po)=%g constr(%g) dr=%g p(%g,%g) po(%g,%g) d(%g,%g) \n", i,j, (p-p0).norm(), d.dot(p-po), dr, p.x,p.y, po.x,po.y, d.x, d.y  );
                // }

                p.add_mul( v, dt );
                //particles.push_back( p );
            }
            if(i==0){
                po = p + d*-dr;
            }else if(i>0){
                //printf( "points_along_rim()[%i] d(%g,%g) p(%g,%g) po(%g,%g) p0(%g,%g) \n", i, d.x,d.y, p.x,p.y, po.x,po.y, p0.x,p0.y );
                particles.push_back( p );
                if     ( i==1 ){ p0=p; }
                else if( i>5  ){ if( (p-p0).norm2()<(dr*dr) ) break; }
            }
            d  = p-po; d.normalize();
            //printf( "points_along_rim()[%i] d(%g,%g) p(%g,%g) po(%g,%g) p0(%g,%g) \n", i, d.x,d.y, p.x,p.y, po.x,po.y, p0.x,p0.y );
            po = p;
            p.add_mul( d, dr );
            // if( bStopByPhi ){
            //     Vec2d rot{p.x,p.y}; rot.normalize();
            //     rot.udiv_cmplx( rot0 );
            //     //printf( "points_along_rim()[%i] p(%g,%g) p0(%g,%g) rot(%g,%g) rot0(%g,%g) rotSign=%g \n", i, p.x,p.y, p0.x, p0.y, rot.x,rot.y, rot0.x,rot0.y, rotSign  ) ;
            //     if( (i>0.1*npmax) && (rot.x>0.9) && (rotSign*(rot.y)<0) ) break;
            // }
        }
    }

    void place_point_offset( double dR,  Vec2i range, Vec3d cog=Vec3dZero ){
        for( int i=range.x; i<range.y; i++ ){
            Vec3d  d = particles[i] - cog;
            double R = d.norm();
            d.mul( (R+dR)/R );
            particles.push_back( d + cog );
        }
    }

    int prepareRim( int nr, Vec2d Rs, Vec2d dLs, Vec3d p0, Vec2d d0=Vec2d{0.0,0.1}, int niter=10000, double Fconv=1e-4, double dt=0.05, double K=1.0, int npmax=1000 ){
        for( int i=0; i<nr; i++ ){
            double R  = Rs.x  + (Rs.y -Rs.x )*i/(nr-1);
            double dL = dLs.x + (dLs.y-dLs.x)*i/(nr-1);
            points_along_rim( R, p0, d0*dL, niter, Fconv, dt, K, npmax );
        }
        return particles.size();
    }

    int prepareRim2( int nr, double* Rs, double dL0, Vec3d p0, Vec2d d0=Vec2d{0.0,0.1}, Vec3d cog=Vec3dZero, int niter=10000, double Fconv=1e-4, double dt=0.05, double K=1.0, int npmax = 1000 ){
        double R0 = Rs[0];
        Vec2i range; range.x=particles.size();
        points_along_rim( R0, p0, d0*dL0, niter, Fconv, dt, K, npmax );
        range.y=particles.size();
        for( int i=1; i<nr; i++ ){
            place_point_offset( Rs[i]-R0, range, cog );
        }
        return range.y-range.x;
    }

    void eval( double dt, double Fconv, double* Es=0, Quat4d* fes1=0, Quat4d* fes2=0 ){
        printf( "DipoleMap::eval() REQ(%g,%g,%g,%g) REQ2(%g,%g,%g,%g) \n", REQ.x,REQ.y,REQ.z,REQ.w, REQ2.x,REQ2.y,REQ2.z,REQ2.w );
        for(int i=0; i<particles.size(); i++){
            Vec3d p  = particles [i];
            Vec3d p2 = particles2[i];
            Quat4d fe  = ff->evalLJQs( p,  REQ,  1.0 );
            Quat4d fe2 = ff->evalLJQs( p2, REQ2, 1.0 );
            if( Es   ) Es[i]   = fe.e + fe2.e;
            if( fes1 ) fes1[i] = fe;
            if( fes2 ) fes2[i] = fe2;   
        }
    }

    int genAndSample( int nr, double* Rs, double dL0, Vec3d p0, Vec2d d0=Vec2d{0.0,0.1}, Vec3d cog=Vec3dZero, int niter=10000, double Fconv=1e-4, double dt=0.05, double K=1.0, int npmax = 1000 ){
        nphi = prepareRim2( nr, Rs, dL0, p0, d0, cog, niter, Fconv, dt, K,npmax );
        //orient( 0.5, -1.0, 0.5, cog );
        //orient_cos( {0.0,0.5},  cog );
        orient_Rgrad(  );
        //orient( 0.000, 1.0, -1, cog );
        FE .resize( particles.size() );
        FE2.resize( particles.size() );
        eval( dt, Fconv, 0, &FE[0], &FE2[0] );
        return nphi;
    }

    void saveToXyz( const char* fname, bool bCharges=false )const{
        FILE* fout = fopen( fname, "w" );
        for(int i=0; i<particles.size(); i++){
            fprintf( fout, "%i \n",  ff->natoms+2   );
            fprintf( fout, "# %i FE %g %g %g %g  FE2 %g %g %g %g \n", i, FE[i].x, FE[i].y, FE[i].z, FE[i].w,   FE2[i].x, FE2[i].y, FE2[i].z, FE2[i].w );
            if(bCharges){
                fprintf( fout, "%s %15.10f   %15.10f   %15.10f  %5.2f\n",  name1, particles[i].x,  particles[i].y,  particles[i].z  , REQ.z );
                fprintf( fout, "%s %15.10f   %15.10f   %15.10f  %5.2f\n",  name2, particles2[i].x, particles2[i].y, particles2[i].z , REQ2.z );
                params->writeXYZ( fout, ff->natoms, ff->atypes, ff->apos, 0, ff->REQs, true, 0, Vec3i{1,1,1}, Mat3dIdentity, false );
            }else{
                fprintf( fout, "%s %15.10f   %15.10f   %15.10f\n",  name1, particles[i].x,  particles[i].y,  particles[i].z  , REQ.z );
                fprintf( fout, "%s %15.10f   %15.10f   %15.10f\n",  name2, particles2[i].x, particles2[i].y, particles2[i].z , REQ2.z );
                params->writeXYZ( fout, ff->natoms, ff->atypes, ff->apos, 0, 0, true, 0, Vec3i{1,1,1}, Mat3dIdentity, false );
            }
        }
        fclose( fout );
    }


};

#endif
