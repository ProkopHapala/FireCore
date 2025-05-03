#ifndef _GLATTRIBS_H_
#define _GLATTRIBS_H_

#include "GLES.h"
#include "quaternion.h"
#include <utility>

namespace GLattrib{
    enum Name {Position, Normal, Color, UV};
}


template<typename T>
struct attrib{
    typedef T type;
    GLattrib::Name name;
};


namespace GLattrib{
    template<typename T> constexpr GLenum type2GLenum();
    template<typename T> constexpr GLint type2count();
    template<Name> constexpr const char* name2str();


    template<attrib...attribs, size_t...i>
    static constexpr bool _check_attribs_impl(std::index_sequence<i...>){
        Name name;
        size_t idx;
        bool invalid = false;
        ((
            name = attribs.name,
            idx = i,
            invalid = invalid || ((name == attribs.name && idx != i) || ...)
        ),...);
        return !invalid;
    }
    template<attrib...attribs>
    static constexpr bool check_attribs(){
        return _check_attribs_impl<attribs...>(std::make_index_sequence<sizeof...(attribs)>{});
    }


    template<> constexpr GLenum type2GLenum<float>() { return GL_FLOAT; }
    template<> constexpr GLenum type2GLenum<Vec2f>() { return GL_FLOAT; }
    template<> constexpr GLenum type2GLenum<Vec3f>() { return GL_FLOAT; }
    template<> constexpr GLenum type2GLenum<Quat4f>(){ return GL_FLOAT; }
    template<> constexpr GLenum type2GLenum<int>()  { return GL_INT; }
    template<> constexpr GLenum type2GLenum<Vec2i>(){ return GL_INT; }
    template<> constexpr GLenum type2GLenum<Vec3i>(){ return GL_INT; }
    template<> constexpr GLenum type2GLenum<unsigned int>(){ return GL_UNSIGNED_INT; }

    template<> constexpr GLint type2count<float>() { return 1; }
    template<> constexpr GLint type2count<Vec2f>() { return 2; }
    template<> constexpr GLint type2count<Vec3f>() { return 3; }
    template<> constexpr GLint type2count<Quat4f>(){ return 4; }
    template<> constexpr GLint type2count<int>()  { return 1; }
    template<> constexpr GLint type2count<Vec2i>(){ return 2; }
    template<> constexpr GLint type2count<Vec3i>(){ return 3; }
    template<> constexpr GLint type2count<unsigned int>(){ return 1; }

    template<> constexpr const char* name2str<Position>(){ return "vPosition"; }
    template<> constexpr const char* name2str<Normal>()  { return "vNormal"; }
    template<> constexpr const char* name2str<Color>()   { return "vColor"; }
    template<> constexpr const char* name2str<UV>()      { return "vUV"; }
};

#define MPOS    attrib<Vec3f>(GLattrib::Position)
#define MNORMAL attrib<Vec3f>(GLattrib::Normal)
#define MCOLOR  attrib<Vec3f>(GLattrib::Color)
#define MUV     attrib<Vec2f>(GLattrib::UV)

#endif // _GLATTRIBS_H_
