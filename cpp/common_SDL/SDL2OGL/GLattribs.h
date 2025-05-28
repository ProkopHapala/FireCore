#ifndef _GLATTRIBS_H_
#define _GLATTRIBS_H_

#include "GLES.h"
#include "quaternion.h"
#include <stdexcept>
#include <utility>
#include <cstddef>

namespace GLattrib{
    enum Name {
        Position = 0,
        Normal   = 1,
        Color    = 2,
        UV       = 3,
    ATTRIB_NAME_MAX}; // ATTRIB_NAME_MAX must be last
}


template<typename T>
struct attrib{
    typedef T type;
    GLattrib::Name name;
};

template<attrib...attribs>
struct attribs_monostate{};


template <typename T, typename ... Ts>
struct GLvertex{
    T first;
    GLvertex<Ts...> next;

    GLvertex(T first, Ts... next) : first(first), next(next...) {}

    template<size_t i> auto& get(){
        if constexpr(i==0) return first;
        else return next.template get<i-1>();
    }

    template<size_t i> static constexpr size_t get_offset() {
        if constexpr(i==0) return __builtin_offsetof(GLvertex<T, Ts...>, first); // TODO: using just offsetof() throws error for some reason, so we use __builtin_offsetof() instead - but it isn't portable
        else return __builtin_offsetof(GLvertex<T, Ts...>, next) + decltype(next)::template get_offset<i-1>();
    }
};
template<typename T>
struct GLvertex<T>{
    T first;

    template <size_t i> auto& get(){
        static_assert(i==0, "Error: index out of bounds"); // TODO: is a static assert the correct way to do this?
        return first;
    }

    template<size_t i> static constexpr size_t get_offset() {
        static_assert(i==0, "Error: index out of bounds"); // note: i should be know at compile time (because this is a constexpr), so static assert is ok
        if constexpr(i==0) return offsetof(GLvertex<T>, first);
    }
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

    static inline const char* name2str(Name name){
        switch(name){
            case Position: return "vPosition";
            case Normal:   return "vNormal";
            case Color:    return "vColor";
            case UV:       return "vUV";
            default:       throw std::runtime_error("Error: invalid name");
        }
    }
};

#define MPOS    attrib<Vec3f>(GLattrib::Position)
#define MNORMAL attrib<Vec3f>(GLattrib::Normal)
#define MCOLOR  attrib<Vec3f>(GLattrib::Color)
#define MUV     attrib<Vec2f>(GLattrib::UV)

#endif // _GLATTRIBS_H_
