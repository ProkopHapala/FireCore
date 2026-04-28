#ifndef Ewald2D_h
#define Ewald2D_h

#include "fastmath.h"
#include "Vec2.h"
#include "Vec3.h"

class Ewald2D { public: 

    Vec2d* Vk=0;

    Vec2d lat_a; 
    Vec2d lat_b;
    int   kmax = 10;

    double alpha = 0.5;

    void project( int n, Vec3d* ps, double* Qs ){        
        for (int ia=0; ia<=kmax; ia++) {
            for (int ib=0; ib<=kmax; ib++) {
                if ( (ia==0) && (ib==0) ) continue;
                Vec2d k = lat_a*ia + lat_b*ib;
                double k2 = k.dot(k);
                if (k2 == 0.0) continue;
                double k_prefactor = exp(-k2/(4.0*alpha*alpha))/k2;
                const int k = ia+ib*kmax;
                Vk[k] = Vec2dZero;
                for (int i=0; i<n; i++) {
                    Vec3d pi   = ps[i];
                    double phi = k.x*pi.x + k.y*pi.y;
                    Vec2d u{ cos(phi), sin(phi) };
                    Vk[k] = u * k_prefactor * Qs[i];
                }
                // double rho2 = rho_real * rho_real + rho_imag * rho_imag;
                // energy += k_prefactor * rho2 / volume;
                // for (int i = 0; i < natoms; i++) {
                //     double force_prefactor = 2.0 * Qs[i] * k_prefactor / volume;
                //     Vec3d force = k * (force_prefactor * (rho_imag * cos_k_dot_r[i] - rho_real * sin_k_dot_r[i]));
                //     forces[i] += force;
                // }
            }
        }
    }

    for (int i = 0; i < natoms; i++) {
            double force_prefactor = 2.0 * Qs[i] * k_prefactor / volume;
            Vec3d force = k * (force_prefactor * (rho_imag * cos_k_dot_r[i] - rho_real * sin_k_dot_r[i]));
            forces[i] += force;
    }



    void apply( int n, Vec3d* ps, double* Qs, Vec3d* fs ){
        for (int ia=0; ia<=kmax; ia++) {
            for (int ib=0; ib<=kmax; ib++) {
                if ( (ia==0) && (ib==0) ) continue;
                Vec2d k = lat_a*ia + lat_b*ib;
                double k2 = k.dot(k);
                if (k2 == 0.0) continue;
                double k_prefactor = exp(-k2/(4.0*alpha*alpha))/k2;
                const int k = ia+ib*kmax;
                Vec2d u;
                for (int i=0; i<n; i++) {
                    Vec3d pi   = ps[i];
                    double phi = k.x*pi.x + k.y*pi.y;
                    u.fromAngle(phi);
                    Vk[i]= u * k_prefactor * Qs[i];
                }
            }
        }
    }


    // void project_fast( Vec2d& p, double Q ){
    //     double phi_a = lat_a.dot(p);
    //     double phi_b = lat_b.dot(p);
    //     Vec2d dua; dua.fromAngle(phi_a);
    //     Vec2d dub; dub.fromAngle(phi_b);
    //     Vec2d ua=Vec2d{1.0,0.0};
    //     for (int na=0; na<=kmax; na++) {
    //         Vec2d ub=Vec2d{1.0,0.0};
    //         for (int nb=0; nb<=kmax; nb++) {
    //             Vec2d u; u.set_mul_cmplx(ua,ub);
    //             Q = Q + u.dot(p);
    //             ub.mul_cmplx(dub);
    //         }
    //         ua.mul_cmplx(dua);
    //     }

    // }
    

}; // class GridFF


#endif
