
#ifndef  Camera_h
#define  Camera_h

#include <math.h>
#include <cstdlib>
#include <stdint.h>

#include "fastmath.h"
#include "Vec3.h"
#include "Mat3.h"
#include "Mat4.h"
#include "quaternion.h"

template<typename T>
class CameraT{ public:
    Vec3T<T>  pos    = (Vec3T<T>){0.0f,0.0f,-50.0f};
    Quat4T<T> qrot   = (Quat4T<T>)Quat4dIdentity;
    T  zoom   = 10.0;
    T  aspect = 1.0;
    const T  zmin   = -10000.0;
    const T  zmax   = 10000.0;
    bool   persp  = true;

    inline const Mat3T<T> rotMat() const { Mat3T<T>m; qrot.toMatrix(m); return m; }
    
    inline void lookAt( Vec3T<T> p, T R ){ pos = p + rotMat().c*-R; }
    //inline void lookAt( Vec3T<T> p, T R ){ Vec3T<T> p_; convert(p,p_); lookAt(p_,R); }

    inline T getTgX()const{ return 1.0/(zoom*aspect); }
    inline T getTgY()const{ return 1.0/(zoom);            }

    inline void word2screenOrtho( const Vec3T<T>& pWord, Vec3T<T>& pScreen ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rotMat().dot_to( p, p );
        pScreen.x = p.x/(2*zoom*aspect);
        pScreen.y = p.y/(2*zoom);
        pScreen.z = (p.z-zmin)/(zmax-zmin);
    }

    inline Vec2T<T> word2pixOrtho( const Vec3T<T>& pWord, const Vec2f& resolution ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rotMat().dot_to( p, p );
        return Vec2f{ resolution.x*(0.5+p.x/(2*zoom*aspect)),
                      resolution.y*(0.5+p.y/(2*zoom)) };
    }

    inline void word2screenPersp( const Vec3T<T>& pWord, Vec3T<T>& pScreen ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rotMat().dot_to( p, p );
        T  resc = zmin/(2*p.z*zoom);
        pScreen.x = p.x*resc/aspect;
        pScreen.y = p.y*resc;
        //pScreen.z = p.z/zmin;        // cz
        //(2*zmin)/w      0       0              0            0
        //0          (2*zmin)/h   0              0            0
        //0               0       (zmax+zmin)/(zmax-zmin)    (2*zmin*zmax)/(zmax-zmin)
        //0               0       0               0           -1
        //------
        //x_  =  ((2*zmin)/w)  * x
        //y_  =  ((2*zmin)/h ) * y
    }

    inline Vec2T<T> word2pixPersp( const Vec3T<T>& pWord, const Vec2f& resolution ) const {
        Vec3T<T> p; p.set_sub(pWord,pos);
        rotMat().dot_to( p, p );
        T  resc = zmin/(2*p.z*zoom);
        return (Vec2f){
            resolution.x*( 0.5 + p.x*resc/aspect ),
            resolution.y*( 0.5 + p.y*resc        ) };
    }

    inline void pix2rayOrtho( const Vec2f& pix, Vec3T<T>& ro ) const {
        //T  resc = 1/zoom;
        T  resc = zoom;
        //printf( "Camera::pix2rayOrtho() pix(%g,%g) rsc %g b(%g,%g,%g) a(%g,%g,%g)\n", pix.a,pix.b, resc,  rot.a.x,rot.a.y,rot.a.z,  rot.b.x,rot.b.y,rot.b.z );
        ro = rotMat().a*(pix.a*resc) + rotMat().b*(pix.b*resc);
    }

    inline void pix2rayPersp( const Vec2f& pix, Vec3T<T>& rd ) const {
        T  resc = 1/zoom;
        rd = rotMat().a*(pix.a*resc) + rotMat().b*(pix.b*resc);
    }

    //inline Vec3T<T> pix2ray( const Vec2f& pix, Vec3T<T>& rd, Vec3T<T>& ro ){
    inline void pix2ray( const Vec2f& pix, Vec3T<T>& rd, Vec3T<T>& ro ){
        //printf(  "Camera::pix2ray() persp %i \n", persp);
        if(persp){
            ro = pos;
            pix2rayPersp( pix, rd );
        }else{
            rd = rotMat().c;
            pix2rayOrtho( pix, ro );
        }
    }

    inline bool pointInFrustrum( Vec3T<T> p ) const {
        p.sub(pos);
        Vec3T<T> c;
        rotMat().dot_to( p, c );
        T  tgx = c.x*zoom*aspect;
        T  tgy = c.y*zoom;
        T  cz  = c.z*zmin;
        return (tgx>-cz)&&(tgx<cz) && (tgy>-cz)&&(tgy<cz) && (c.z>zmin)&&(c.z<zmax);
    }

    inline bool sphereInFrustrum( Vec3T<T> p, T  R ) const {
        p.sub(pos);
        Vec3T<T> c;
        rotMat().dot_to( p, c );
        T  my = c.z*zmin/zoom;
        T  mx = my/aspect + R;  my+=R;
        return (c.x>-mx)&&(c.x<mx) && (c.y>-my)&&(c.y<my) && ((c.z+R)>zmin)&&((c.z-R)<zmax);
    }

    Mat4T<T> projectionMatrix() const {
        Mat4T<T> m;
        m.setOne();
        if(persp){
            m.setPerspective( -zoom*aspect, zoom*aspect, -zoom, zoom, zmin, zmax );
        }else{
            m.setOrthographic( -zoom*aspect, zoom*aspect, -zoom, zoom, zmin, zmax );
        }
        return m;
    }

    Mat4T<T> viewMatrix() const {
        Mat4T<T> m;
        m.setOne();
        Mat3T<T> mrot;
        qrot.toMatrix_T(mrot);
        m.setRot( mrot );
        m.setPos( -pos );
        
        return m;
    }

    Mat4T<T> viewProjectionMatrix() const {
        Mat4T<T> v = viewMatrix();
        Mat4T<T> vp = v.mmulR(projectionMatrix());
        return vp;
    }

};

using Camera   = CameraT<float>;
using Camera_d = CameraT<double>;

#endif

