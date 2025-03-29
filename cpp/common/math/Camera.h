
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
class CameraT{ 
private:
    Vec3T<T>  _pos    = (Vec3T<T>){0.0f,0.0f,-50.0f};
    Quat4T<T> _qrot   = (Quat4T<T>)Quat4dIdentity;
    T  _zoom   = 10.0;
    T  _aspect = 1.0;
    bool   _persp  = true;

    bool update_vp_mat = true;
    Mat4T<T> _viewMatrix;
    Mat4T<T> _projMat;
    Mat4T<T> _viewProjMatrix;

public:
    const T  zmin   = -10000.0;
    const T  zmax   = 10000.0;

    inline Vec3T<T> pos() const {return _pos;};
    inline void setPos( Vec3T<T> p ){ _pos = p; update_vp_mat = true; }
    inline void shift( Vec3T<T> p ){ _pos += p; update_vp_mat = true; }

    inline Quat4T<T> qrot() const { return _qrot; }
    inline void setQrot( Quat4T<T> q ){ _qrot = q; update_vp_mat = true; }
    inline void dyaw( T a ){ _qrot.dyaw(a); update_vp_mat = true; }
    inline void dpitch( T a ){ _qrot.dpitch(a); update_vp_mat = true; }
    inline void qrotQmul_T( Quat4T<T> q ){ _qrot.qmul_T(q); update_vp_mat = true; }

    inline T zoom() const { return _zoom; }
    inline void setZoom( T z ){ _zoom = z; update_vp_mat = true; }

    inline T aspect() const { return _aspect; }
    inline void setAspect( T a ){ _aspect = a; update_vp_mat = true; }

    inline bool persp() const { return _persp; }
    inline void setPersp( bool p ){ _persp = p; update_vp_mat = true; }

    inline const Mat3T<T> rotMat() const { Mat3T<T>m; _qrot.toMatrix(m); return m; }
    
    inline void lookAt( Vec3T<T> p, T R ){ setPos(p + rotMat().c*-R); }
    //inline void lookAt( Vec3T<T> p, T R ){ Vec3T<T> p_; convert(p,p_); lookAt(p_,R); }

    inline T getTgX()const{ return 1.0/(_zoom*_aspect); }
    inline T getTgY()const{ return 1.0/(_zoom);            }

    inline void word2screenOrtho( const Vec3T<T>& pWord, Vec3T<T>& pScreen ) const {
        Vec3T<T> p; p.set_sub(pWord,_pos);
        rotMat().dot_to( p, p );
        pScreen.x = p.x/(2*_zoom*_aspect);
        pScreen.y = p.y/(2*_zoom);
        pScreen.z = (p.z-zmin)/(zmax-zmin);
    }

    inline Vec2T<T> word2pixOrtho( const Vec3T<T>& pWord, const Vec2f& resolution ) const {
        Vec3T<T> p; p.set_sub(pWord,_pos);
        rotMat().dot_to( p, p );
        return Vec2f{ resolution.x*(0.5+p.x/(2*_zoom*_aspect)),
                      resolution.y*(0.5+p.y/(2*_zoom)) };
    }

    inline void word2screenPersp( const Vec3T<T>& pWord, Vec3T<T>& pScreen ) const {
        Vec3T<T> p; p.set_sub(pWord,_pos);
        rotMat().dot_to( p, p );
        T  resc = zmin/(2*p.z*_zoom);
        pScreen.x = p.x*resc/_aspect;
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
        Vec3T<T> p; p.set_sub(pWord,_pos);
        rotMat().dot_to( p, p );
        T  resc = zmin/(2*p.z*_zoom);
        return (Vec2f){
            resolution.x*( 0.5 + p.x*resc/_aspect ),
            resolution.y*( 0.5 + p.y*resc        ) };
    }

    inline void pix2rayOrtho( const Vec2f& pix, Vec3T<T>& ro ) const {
        //T  resc = 1/zoom;
        T  resc = _zoom;
        //printf( "Camera::pix2rayOrtho() pix(%g,%g) rsc %g b(%g,%g,%g) a(%g,%g,%g)\n", pix.a,pix.b, resc,  rot.a.x,rot.a.y,rot.a.z,  rot.b.x,rot.b.y,rot.b.z );
        ro = rotMat().a*(pix.a*resc) + rotMat().b*(pix.b*resc);
    }

    inline void pix2rayPersp( const Vec2f& pix, Vec3T<T>& rd ) const {
        T  resc = 1/_zoom;
        rd = rotMat().a*(pix.a*resc) + rotMat().b*(pix.b*resc);
    }

    //inline Vec3T<T> pix2ray( const Vec2f& pix, Vec3T<T>& rd, Vec3T<T>& ro ){
    inline void pix2ray( const Vec2f& pix, Vec3T<T>& rd, Vec3T<T>& ro ){
        //printf(  "Camera::pix2ray() persp %i \n", persp);
        if(_persp){
            ro = _pos;
            pix2rayPersp( pix, rd );
        }else{
            rd = rotMat().c;
            pix2rayOrtho( pix, ro );
        }
    }

    inline bool pointInFrustrum( Vec3T<T> p ) const {
        p.sub(_pos);
        Vec3T<T> c;
        rotMat().dot_to( p, c );
        T  tgx = c.x*_zoom*_aspect;
        T  tgy = c.y*_zoom;
        T  cz  = c.z*zmin;
        return (tgx>-cz)&&(tgx<cz) && (tgy>-cz)&&(tgy<cz) && (c.z>zmin)&&(c.z<zmax);
    }

    inline bool sphereInFrustrum( Vec3T<T> p, T  R ) const {
        p.sub(_pos);
        Vec3T<T> c;
        rotMat().dot_to( p, c );
        T  my = c.z*zmin/_zoom;
        T  mx = my/_aspect + R;  my+=R;
        return (c.x>-mx)&&(c.x<mx) && (c.y>-my)&&(c.y<my) && ((c.z+R)>zmin)&&((c.z-R)<zmax);
    }

    void recalculate_vpMat() {
        // view matrix
        _viewMatrix.setOne();
        Mat3T<T> mrot; _qrot.toMatrix_T(mrot);
        _viewMatrix.setRot( mrot );
        _viewMatrix.setPos( -_pos );

        // projection matrix
        _projMat.setOne();
        if(_persp){
            _projMat.setPerspective( -_zoom*_aspect, _zoom*_aspect, -_zoom, _zoom, zmin, zmax );
        }else{
            _projMat.setOrthographic( -_zoom*_aspect, _zoom*_aspect, -_zoom, _zoom, zmin, zmax );
        }

        // view projection matrix
        _viewProjMatrix = _viewMatrix;
        _viewProjMatrix.mmulL(_projMat);

        update_vp_mat = false;
    }

    Mat4T<T> projectionMatrix() {
        if (update_vp_mat) recalculate_vpMat();
        return _projMat;
    }

    Mat4T<T> viewMatrix() {
        if (update_vp_mat) recalculate_vpMat();
        return _viewMatrix;
    }

    Mat4T<T> viewProjectionMatrix() {
        if (update_vp_mat) recalculate_vpMat();
        return _viewProjMatrix;
    }

};

using Camera   = CameraT<float>;
using Camera_d = CameraT<double>;

#endif

