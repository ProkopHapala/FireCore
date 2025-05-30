#ifndef GLINSTANCEDMESH_H
#define GLINSTANCEDMESH_H

#include "GLES.h"
#include "GLattribs.h"
#include "GLvbo.h"
#include "Shader.h"
#include "Camera.h"

template <class T, template<attrib...> class Template>
struct is_specialization : std::false_type {};

template <template<attrib...> class Template, attrib... Args>
struct is_specialization<Template<Args...>, Template> : std::true_type {};

template <class T, template<attrib...> class Template>
constexpr bool is_specialization_v = is_specialization<T, Template>::value;








template<class VertVboType, attrib...instAttribs>
class GLInstancedMeshBase {private:
    static_assert(is_specialization_v<VertVboType, GLvbo>, "MeshType must be specialization of GLMeshBase");
    static_assert(GLattrib::check_attribs<instAttribs...>(), "Attribute list cannot contain duplicate names.");
    
    template<attrib...vertAttribs>
    static constexpr bool check_attribs_merged(attribs_monostate<vertAttribs...>){
        return GLattrib::check_attribs<vertAttribs..., instAttribs...>();
    }
    static_assert(check_attribs_merged(typename VertVboType::attr_monostate{}), "MeshType and instAttribs must have disjunct attribute names.");

public:
    using vertex = GLvertex<typename decltype(instAttribs)::type ...>;
    using attrIdxSeq = std::make_index_sequence<sizeof...(instAttribs)>;
    using attr_monostate = attribs_monostate<instAttribs...>;

private:
    template<attrib...vertAttribs>
    static constexpr Shader* get_default_shader(attribs_monostate<vertAttribs...>){
        return defaultShader<vertAttribs..., instAttribs...>;
    }

public:
    GLenum drawMode = GL_TRIANGLES;
    VertVboType* verts;
    GLvbo<instAttribs...>* instances;
    Shader* shader;
    GLuniformSet uniforms;

    GLInstancedMeshBase(VertVboType* vbo, GLenum draw_mode=GL_TRIANGLES, Shader* shader=get_default_shader(typename VertVboType::attr_monostate{}))
        : shader(shader), verts(vbo), drawMode(draw_mode), instances(new GLvbo<instAttribs...>()) {}
    
    void addInstance(typename decltype(instAttribs)::type ... args){
        instances->push_back(vertex(args...));
    }

    inline void bind_sync_vbos(){
        verts->bind( [this](GLattrib::Name name){return shader->attrName2Loc(name);} );
        instances->bind( [this](GLattrib::Name name){return shader->attrName2Loc(name);}, 1 );
    }

    inline void draw(GLenum drawMode=0){
        if (drawMode == 0) drawMode = this->drawMode;
        if (GLES::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }

        bind_sync_vbos();
        shader->setUniforms(uniforms);
        shader->setUniform4m("uMVPMatrix", GLES::active_camera->viewProjectionMatrix());
        shader->use();
        glDrawArraysInstanced(drawMode, 0, verts->size(), instances->size());
    }
};



#endif // GLINSTANCEDMESH_H
