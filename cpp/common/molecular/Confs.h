
#ifndef Confs_h
#define Confs_h

#include "quaternion.h"


/*

F_i = d_( L - L0 ) / d_xi

L^2 = Sum_ig{   Kc*| cg1_ig - cg2_ig |^2  +    Ku* ( 1-cos(u1_ig,u2_ig) )  + Kv*( 1 - cos(u1_ig,u2_ig) }


cos(u1,u2) =  <u1|u2>/(|u1||u2|)


cg1 = Sum_i{ p_i * wc_i  } / Sum_i{ wc_i }

d{c1_ig}/d{p_i} =  wc_i / Sum_i{ wc_i  }

u1 = Sum_i{ p_i * wu_i  } - cg1 * Sum_i{wu_i}

d{u1}/d{x1i} =    wu_i -   wc_i  Sum_i{wu_i}/Sum_i{wc_i}   =   Wu_i     .... just constant (renormalized)
d{u1}/d{x1i} =    wu_i -   wc_i   ... if vec{wu} and vec{wc} are normalized 

d{ cos(u1,u2) }/dx1i = d{ <u1|u2>/(|u1|*|u2|) }/dx1i =     1/(|u1||u2|)  *  d{ <u1|u2> }/dx1i     +   <u1|u2>/|u2|  *  d{ 1/|u1| }/dx1i 
d{ <u1|u2> }/dx1i = x2i
d{ 1/|u1|  }/dx1i = -1/|u1|^3 * x1i
d{ cos(u1,u2) }/dx1_i =  x2_i /(|u1||u2|)  -  x1_i / <u1|u2>/( |u1|^3 * |u2| )

d{ cos(u1,u2) }/dx1_i =  x2_i * ilu1*ilu2  -  x1_i * cos(u1,u2) * ilu1*ilu1 =  Wui * ilu1*( x2_i*ilu2  -  x1_i*ilu1* cos(u1,u2) )   


=> need to store 
ilu1 = 1/|u1| 
ilu2 = 1/|u2|
ilv1 = 1s/|v1|
ilv2 = 1/|v2| 

*/


class Group{ public:
    int     n = 0;   // number of atoms
    int*    iatoms;  // atom indexes belonging to the group
    Quat4f* cs;      // projection coefficients

    inline float distance2( Quat4f* ps1, Quat4f* ps2 )const{
        float lr2 = 0.f; 
        for(int i=0; i<n; i++){
            Vec3f d = ps1[i].f - ps2[i].f;
            //lr2 += d.norm2()*w;  // weighted?
            lr2 += d.norm2();
        }
        return lr2;
    }

    inline Vec3f cog( Quat4f* ps )const{
        Vec3f cog = Vec3fZero; 
        for(int i=0; i<n; i++){
            int ia = iatoms[i];
            //cog.add_mul( ps[ia].f, w );       // weighted?
            cog.add( ps[ia].f ); 
        }
        return cog;
    }

    void pose( Quat4f* ps, Vec3f& cog, Vec3f& u, Vec3f& v )const{
        cog = Vec3fZero;
        u   = Vec3fZero;
        v   = Vec3fZero;
        Quat4f Csum = Quat4fZero;
        for(int i=0; i<n; i++){
            int ia          = iatoms[i];
            const Quat4f& C = cs[i];
            const Vec3f&  p = ps[ia].f;
            Csum.add(C);
            cog.add_mul( p, C.w );
            u.add_mul  ( p, C.x );
            v.add_mul  ( p, C.y );
        }
        cog.mul( 1/Csum.w );
        u.add_mul( cog, -Csum.x );
        v.add_mul( cog, -Csum.y );
    }

    void apply_pose_force( Quat4f* fs, const Vec3f& fcog, const Vec3f& fu, const Vec3f& fv )const{
        for(int i=0; i<n; i++){
            int ia          = iatoms[i];
            const Quat4f& C = cs    [i];
            Vec3f fi = Vec3fZero;
            fi.add_mul( fcog, C.w );
            fi.add_mul( fu,   C.x ); // what about subtracted cog ?
            fi.add_mul( fv,   C.y );
            fs[i].f.add(fi);
        }
    }

}; // class Group{

class Confs{ public:
    int nSys  = 0;  // number of systems
    int natom = 0;  // number of atoms
    int nvec  = 0;  // number of quaternion per system

    Vec3f KL{ 1.0, 1.0, 1.0 };  // relative weight of translation and rotation distance

    Quat4f* atoms      =0;
    Quat4f* aforces    =0;

    std::vector<Group> groups;

    float dist2( int isys1, int isys2 ){   // distance between two systems
        int ng = groups.size();
        int i1 = isys1*nvec;
        int i2 = isys2*nvec;
        float l2 = 0; 
        for(int ig=0; ig<ng; ig++ ){
            const Group& G = groups[ig];
            //Vec3f L = pose_distance( ps1, ps2 );
            Vec3f cog1, u1, v1; G.pose( atoms+i1, cog1, u1, v1 );
            Vec3f cog2, u2, v2; G.pose( atoms+i2, cog2, u2, v2 );
            Vec3f L2 = Vec3f{ cog1.dist2(cog2), 1-u1.cos_v(u2), 1-v1.cos_v(v2) };
            l2     += L2.dot( KL );
        }
        return l2;
    }

    float collision_force( int isys1, int isys2, float l0, float k ){   // distance between two systems
        int ng = groups.size();
        Vec3f Ls[ng];
        int i1 = isys1*nvec;
        int i2 = isys2*nvec;
        float l2 = 0; 
        // ---- evaluate configuration distances (reduction)
        for(int ig=0; ig<ng; ig++ ){
            const Group& G = groups[ig];
            //Vec3f L = pose_distance( ps1, ps2 );
            Vec3f cog1, u1, v1; G.pose( atoms+i1, cog1, u1, v1 );
            Vec3f cog2, u2, v2; G.pose( atoms+i2, cog2, u2, v2 );
            Vec3f L2 = Vec3f{ cog1.dist2(cog2), 1-u1.cos_v(u2), 1-v1.cos_v(v2) };
            Ls[ig].set_mul(L2,KL); 
            l2      += L2.dot( KL );

        }
        // ---- distance to force
        float l   = sqrt(l2);
        float dl  = l - l0;
        float f_r = k*dl/l;
        // ---- apply force
        for(int ig=0; ig<ng; ig++ ){
            const Group& G = groups[ig];
            Vec3f fcog1;  // TO BE DONE 
            Vec3f fu1;    // TO BE DONE
            Vec3f fv1;    // TO BE DONE
            G.apply_pose_force( aforces+i1, fcog1, fu1, fv1 );
            Vec3f fcog2;  // TO BE DONE
            Vec3f fu2;    // TO BE DONE
            Vec3f fv2;    // TO BE DONE
            G.apply_pose_force( aforces+i2, fcog2, fu2, fv2 );
        }
        return l;
    }

};   // class Confs{ public: 


#endif
