#ifndef _GLuniform_H_
#define _GLuniform_H_

#include "GLES.h"
#include "Mat4.h"
#include "GLTexture.h"
#include <unordered_map>
#include <string>
#include "constexprString.h"

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

    inline bool operator==(const GLuniform& other) const {
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
            //case m3: return data.m3 == other.data.m3; // comparing matrices is too slow, and they rarely are the same
            //case m4: return data.m4 == other.data.m4;
            case tex: return data.tex == other.data.tex;
            default: return false;
        }
    };

    inline void apply(GLuint location) const {
        switch (type) {
            case GLuniform::f1:   glUniform1f (location, data.f1); break;
            case GLuniform::f2:   glUniform2f (location, data.f2.x, data.f2.y); break;
            case GLuniform::f3:   glUniform3f (location, data.f3.x, data.f3.y, data.f3.z); break;
            case GLuniform::f4:   glUniform4f (location, data.f4.x, data.f4.y, data.f4.z, data.f4.w); break;
            case GLuniform::i1:   glUniform1i (location, data.i1); break;
            case GLuniform::i2:   glUniform2i (location, data.i2.x, data.i2.y); break;
            case GLuniform::i3:   glUniform3i (location, data.i3.x, data.i3.y, data.i3.z); break;
            case GLuniform::i4:   glUniform4i (location, data.i4.x, data.i4.y, data.i4.z, data.i4.w); break;
            case GLuniform::ui1:  glUniform1ui(location, data.ui1); break;
            case GLuniform::ui2:  glUniform2ui(location, data.ui2.x, data.ui2.y); break;
            case GLuniform::ui3:  glUniform3ui(location, data.ui3.x, data.ui3.y, data.ui3.z); break;
            case GLuniform::ui4:  glUniform4ui(location, data.ui4.x, data.ui4.y, data.ui4.z, data.ui4.w); break;
            case GLuniform::m3:   glUniformMatrix3fv(location, 1, GL_FALSE, data.m3.array); break;
            case GLuniform::m4:   glUniformMatrix4fv(location, 1, GL_FALSE, data.m4.array); break;
            case GLuniform::tex:  printf("ERROR: setting texture uniforms in shader is not yet implemented. Please set it in GLMesh instead.\n"); break;
            default: printf("ERROR: invalid uniform type\n");
        }
    }
};



class GLuniformSet{
public:
    using common_uniforms_enum = constexprStringList<
        "uMVPMatrix",
        "uColor",
        "uTexture",
        "uDepth",
        "uPos",
        "uRadius",
        "uTexture1",
        "uTexture2",
        "uDepth1",
        "uDepth2",
        "uMVPinv"
    >;

    GLuniform enum_uniforms[common_uniforms_enum::size()];
    std::unordered_map<std::string, GLuniform> named_uniforms;

    GLuniformSet(){
        for (int i=0; i<common_uniforms_enum::size(); i++) enum_uniforms[i].type = GLuniform::NONE;
    }

    template <constexprString name> inline void set(GLuniform u){
        if constexpr(common_uniforms_enum::template contains<name>()){
            enum_uniforms[common_uniforms_enum::template get_idx<name>()] = u;
        }else{
            named_uniforms[name] = u;
        }
    }

    template <constexprString name> inline void set1f (GLfloat value)       { set<name>({.type=GLuniform::f1,  .data={.f1 =value}}); }
    template <constexprString name> inline void set2f (Vec2T<GLfloat> value){ set<name>({.type=GLuniform::f2,  .data={.f2 =value}}); }
    template <constexprString name> inline void set3f (Vec3T<GLfloat> value){ set<name>({.type=GLuniform::f3,  .data={.f3 =value}}); }
    template <constexprString name> inline void set4f (Vec4T<GLfloat> value){ set<name>({.type=GLuniform::f4,  .data={.f4 =value}}); }
    template <constexprString name> inline void set1i (GLint value)         { set<name>({.type=GLuniform::i1,  .data={.i1 =value}}); }
    template <constexprString name> inline void set2i (Vec2T<GLint> value)  { set<name>({.type=GLuniform::i2,  .data={.i2 =value}}); }
    template <constexprString name> inline void set3i (Vec3T<GLint> value)  { set<name>({.type=GLuniform::i3,  .data={.i3 =value}}); }
    template <constexprString name> inline void set4i (Vec4T<GLint> value)  { set<name>({.type=GLuniform::i4,  .data={.i4 =value}}); }
    template <constexprString name> inline void set1ui(GLuint value)        { set<name>({.type=GLuniform::ui1, .data={.ui1=value}}); }
    template <constexprString name> inline void set2ui(Vec2T<GLuint> value) { set<name>({.type=GLuniform::ui2, .data={.ui2=value}}); }
    template <constexprString name> inline void set3ui(Vec3T<GLuint> value) { set<name>({.type=GLuniform::ui3, .data={.ui3=value}}); }
    template <constexprString name> inline void set4ui(Vec4T<GLuint> value) { set<name>({.type=GLuniform::ui4, .data={.ui4=value}}); }
    template <constexprString name> inline void set3m (Mat3T<GLfloat> value){ set<name>({.type=GLuniform::m3,  .data={.m3 =value}}); }
    template <constexprString name> inline void set4m (Mat4T<GLfloat> value){ set<name>({.type=GLuniform::m4,  .data={.m4 =value}}); }
    template <constexprString name> inline void setTex(GLTexture*     value){ set<name>({.type=GLuniform::tex, .data={.tex=value}}); }
};

#endif // _GLuniform_H_
