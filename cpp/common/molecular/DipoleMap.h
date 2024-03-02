
#ifndef DipoleMap_h
#define DipoleMap_h

#include "NBFF.h"

class DipoleMap {
    // https://docs.lammps.org/Howto_tip5p.html
    // TIP5P - Transferable Intermolecular Potential 5 Point   - 5 particles   (  O(0.0e), H1(+0.241e), H1(+0.241), E1(-0.241e), E2(-0.241e) ), E ~ 0.7A from Oxygen
    
    const NBFF* ff=0;
    Quat4d REQ { 1.187, sqrt(0.0006808    ), -0.24, 0.0 };
    Quat4d REQ2{ 1.750, sqrt(0.00260184625), +0.24, 0.0 };
    double Lbond = 1.0;
    std::vector<Vec3d> particles;
    std::vector<Vec3d> particles2;
    std::vector<int>   pivots; // index of atom to which the particle is attached

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

    void orient( double cEfield, double cPivot ){
        bool bEfield = (cEfield>0);  if(!bEfield) cEfield = 0;
        bool bPivot  = (cPivot >0);  if(!bPivot ) cPivot  = 0;
        for(int i=0; i<particles.size(); i++){
            Vec3d p  = particles[i];
            Vec3d p2 = particles2[i];

            Vec3d d = Vec3dZero;
            if(bEfield){
                Quat4d fe = ff->evalLJQs( p, { 1.0, 0.0, 1.0 }, 1.0 );
                fe.f.normalize();
                d.add_mul( fe.f , cEfield );
            }
            if(bPivot){
                int    imin=-1;
                double r2min = 1e+300;
                for(int j=0; j<ff->natoms; j++){
                    Vec3d d = p - ff->apos[j];
                    double r2 = d.norm2();
                    if(r2<r2min){ r2min=r2; imin=j; }
                }
                pivots[i] = imin;
                Vec3d d_ = (p - ff->apos[imin]); d.normalize();
                d.add_mul( d_, cPivot);
            }   
            particles2.push_back( p + d* (Lbond / (cEfield+cPivot) ) );
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


    void eval( double dt, double Fconv, double* Es=0, Quat4d* fes1=0, Quat4d* fes2=0 ){
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


};

#endif
