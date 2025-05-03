#ifndef _GLuniform_H_
#define _GLuniform_H_

#include "GLES.h"
#include "Mat4.h"
#include "GLTexture.h"

struct GLuniform{
    enum {NONE, f1, f2, f3, f4, i1, i2, i3, i4, ui1, ui2, ui3, ui4, m2, m3, m4, tex} type;
    union data{
        GLfloat f1;
        Vec2T<GLfloat> f2;
        Vec3T<GLfloat> f3;
        Vec4T<GLfloat> f4;

        GLint i1;
        Vec2T<GLint> i2;
        Vec3T<GLint> i3;
        Vec4T<GLint> i4;

        GLuint ui1;
        Vec2T<GLuint> ui2;
        Vec3T<GLuint> ui3;
        Vec4T<GLuint> ui4;

        //Mat2T<GLfloat> m2;
        Mat3T<GLfloat> m3;
        Mat4T<GLfloat> m4;
        GLTexture* tex;
    }data;

    bool operator==(const GLuniform& other) const {
        if (type != other.type) return false;
        switch (type) {
            case NONE: return true;
            case f1: return data.f1 == other.data.f1;
            case f2: return data.f2 == other.data.f2;
            case f3: return data.f3 == other.data.f3;
            case f4: return data.f4 == other.data.f4;
            case i1: return data.i1 == other.data.i1;
            case i2: return data.i2 == other.data.i2;
            case i3: return data.i3 == other.data.i3;
            case i4: return data.i4 == other.data.i4;
            case ui1: return data.ui1 == other.data.ui1;
            case ui2: return data.ui2 == other.data.ui2;
            case ui3: return data.ui3 == other.data.ui3;
            case ui4: return data.ui4 == other.data.ui4;
            case m3: return data.m3 == other.data.m3;
            case m4: return data.m4 == other.data.m4;
            case tex: return data.tex == other.data.tex;
            default: return false;
        }
    };
};

#endif // _GLuniform_H_
