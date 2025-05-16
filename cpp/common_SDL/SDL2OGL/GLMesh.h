#ifndef _GLMesh_H_
#define _GLMesh_H_


#include <cstddef>
#include <unordered_map>
#include <utility>
#include "GLattribs.h"
#include "GLES.h"
#include "GLuniform.h"
#include "Shader.h"
#include "Camera.h"
#include "quaternion.h"

template<attrib ... attribs>
class GLMeshBase{
private:
static_assert(GLattrib::check_attribs<attribs...>(), "ERROR: attribute list cannot contain duplicate names.");

    template <typename T, typename ... Ts>
    struct vert{
        T first;
        vert<Ts...> rest;

        vert(T first, Ts... rest) : first(first), rest(rest...) {}

        template<size_t i> auto& get(){
            if constexpr(i==0) return first;
            else return rest.template get<i-1>();
        }

        template<size_t i> static constexpr size_t get_offset() {
            if constexpr(i==0) return __builtin_offsetof(vert<T, Ts...>, first); // TODO: using just offsetof() throws error for some reason, so we use __builtin_offsetof() instead - but it isn't portable
            else return __builtin_offsetof(vert<T, Ts...>, rest) + decltype(rest)::template get_offset<i-1>();
        }
    };

    template<typename T>
    struct vert<T>{
        T first;

        template <size_t i> auto& get(){
            static_assert(i==0, "Error: index out of bounds"); // TODO: is a static assert the correct way to do this?
            return first;
        }

        template<size_t i> static constexpr size_t get_offset() {
            static_assert(i==0, "Error: index out of bounds"); // note: i should be know at compile time (because this is a constexpr), so static assert is ok
            if constexpr(i==0) return offsetof(vert<T>, first);
        }
    };

public:
    using vertex = vert<typename decltype(attribs)::type ...>;
    using attrIdxSeq = std::make_index_sequence<sizeof...(attribs)>;

private:
    std::vector<vertex> vertices;
    std::vector<GLuniform> uniforms;

    std::unordered_map<std::string, GLuniform> lazy_uniforms; // TODO: figure out a better system for lazy initialization

    GLenum usage = GL_STATIC_DRAW;
    GLuint vbo = 0;
    bool vbo_sync = true;

    inline void ensure_vbo(){
        if (vbo) return;
        glGenBuffers(1, &vbo);
        if constexpr (sizeof...(attribs) > 16){
            printf("WARNING: mesh has more than 16 vertex attributes, which may not be supported by OpenGL ES 3.0.\n");
        }
        GLint max_attribs = -1;
        glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &max_attribs);
        if (sizeof...(attribs) > max_attribs){
            printf("ERROR: mesh has more than %d vertex attributes, which is more than the maximum supported by OpenGL ES 3.0.\n", max_attribs);
        }

        for (auto i : lazy_uniforms) {
            setUniformName(i.first.c_str(), i.second);
        }
        lazy_uniforms.clear();
    }

    inline void bind_vbo(){ _bind_vbo_impl(attrIdxSeq{}); }
    template<size_t ... attrIdx>
    inline void _bind_vbo_impl(std::index_sequence<attrIdx...> seq){
        ensure_vbo();

        if (GLES::currentGL_ARRAY_BUFFER == vbo) return;

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        
        (glVertexAttribPointer(
            shader->getAttribLoc(attrIdx), // TODO: should we instead use glSetAttribLocation() in shader to set (attr_loc = i) ?
            GLattrib::type2count <typename decltype(attribs)::type>(),
            GLattrib::type2GLenum<typename decltype(attribs)::type>(),
            GL_FALSE,
            sizeof(vertex),
            (void*)vertex::template get_offset<attrIdx>())
        ,...);
    }

public:
    GLenum drawMode;
    Shader<attribs...>* shader;
    // TODO: texture? color?

    GLMeshBase(GLenum drawMode=GL_TRIANGLES, GLenum usage=GL_STATIC_DRAW, Shader<attribs...>* shader=defaultcolorShader<attribs...>)
        : drawMode(drawMode), usage(usage), shader(shader) {}
    ~GLMeshBase(){ if (vbo) glDeleteBuffers(1, &vbo); }

    void clear(){
        vertices.clear();
        vbo_sync = false;
    }

    void addVertex(typename decltype(attribs)::type ... args){
        vertices.push_back(vertex(args...));
        vbo_sync = false;
    }
    void addVertex_strip(typename decltype(attribs)::type ... args){
        vertex v = vertex(args...);

        if (vertexCount() >= 3){
            vertex v1 = vertices[vertices.size()-1];
            vertex v2 = vertices[vertices.size()-2];

            vertices.push_back(v2);
            vertices.push_back(v1);
        }

        vertices.push_back(v);
        vbo_sync = false;
    }
    inline int vertexCount() const { return vertices.size(); }


    void setUniformName(const char* name, GLuniform u){
        if (vbo == 0){
            lazy_uniforms[std::string(name)] = u;
            return;
        }
        setUniformLoc(shader->getUniformLocation(name), u);
    }

    void setUniformLoc(GLuint loc, GLuniform u){
        if (loc == -1) return;
        while (uniforms.size() <= loc) uniforms.push_back({.type=GLuniform::NONE});
        uniforms[loc] = u;
    }

    void setUniform1f(const char* name, GLfloat v)        {setUniformName(name, {.type=GLuniform::f1, .data={.f1=v}}); }
    void setUniform2f(const char* name, Vec2T<GLfloat> v) {setUniformName(name, {.type=GLuniform::f2, .data={.f2=v}}); }
    void setUniform3f(const char* name, Vec3T<GLfloat> v) {setUniformName(name, {.type=GLuniform::f3, .data={.f3=v}}); }
    void setUniform4f(const char* name, Vec4T<GLfloat> v) {setUniformName(name, {.type=GLuniform::f4, .data={.f4=v}}); }
    void setUniform1i(const char* name, GLint v)          {setUniformName(name, {.type=GLuniform::i1, .data={.i1=v}}); }
    void setUniform2i(const char* name, Vec2T<GLint> v)   {setUniformName(name, {.type=GLuniform::i2, .data={.i2=v}}); }
    void setUniform3i(const char* name, Vec3T<GLint> v)   {setUniformName(name, {.type=GLuniform::i3, .data={.i3=v}}); }
    void setUniform4i(const char* name, Vec4T<GLint> v)   {setUniformName(name, {.type=GLuniform::i4, .data={.i4=v}}); }
    void setUniformMatrix3f(const char* name, Mat3T<GLfloat> v) {setUniformName(name, {.type=GLuniform::m3, .data={.m3=v}}); }
    void setUniformMatrix4f(const char* name, Mat4T<GLfloat> v) {setUniformName(name, {.type=GLuniform::m4, .data={.m4=v}}); }  
    void setUniformTex(const char* name, GLTexture* tex) {setUniformName(name, {.type=GLuniform::tex, .data={.tex=tex}});} 

    inline void bind_sync_vbo(){
        bind_vbo();
        if (vbo_sync) return;

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), usage);
        vbo_sync = true;
    }

    inline void draw(GLenum drawMode=0){
        if (drawMode == 0) drawMode = this->drawMode;
        bind_sync_vbo();

        int texi = 0;
        for (int i=0; i<uniforms.size(); i++){
            if (uniforms[i].type == GLuniform::NONE) continue;
            if (uniforms[i].type == GLuniform::tex){
                glActiveTexture(GL_TEXTURE0 + texi);
                uniforms[i].data.tex->bind();
                shader->setUniformLoc(i, {.type=GLuniform::i1, .data={.i1=texi}});
                texi++;
                continue;
            }
            shader->setUniformLoc(i, uniforms[i]);
        }
        shader->use();

        glDrawArrays(drawMode, 0, vertexCount());
    }

    int x;
};


template<attrib...attribs>
class GLMesh : public GLMeshBase<attribs...>{
public:
    using base = GLMeshBase<attribs...>;
    GLMesh(GLenum drawMode=GL_TRIANGLES, GLenum usage=GL_STATIC_DRAW, Shader<attribs...>* shader=defaultcolorShader<attribs...>)
        : base(drawMode, usage, shader) {}


    void drawMVP(Mat4f mvp, GLenum drawMode=0){
        base::shader->setuMVPMatrix(mvp);
        base::draw(drawMode);
    }

    inline void draw(GLenum drawMode){draw(Vec3fZero, Vec3fOne, drawMode);}
    inline void draw(Vec3f position, float scale, GLenum drawMode=0){draw(position, (Vec3f){scale, scale, scale}, drawMode);}
    void draw(Vec3f position=Vec3fZero, Vec3f scale={1, 1, 1}, GLenum drawMode=0){
        if (GLES::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }
    
        // scale + translation
        Mat4f mvpMatrix;
        mvpMatrix.setDiag(scale.x, scale.y, scale.z, 1);
        mvpMatrix.setPos(position);

        // view + projection
        mvpMatrix.mmulL(GLES::active_camera->viewProjectionMatrix());
    
        drawMVP(mvpMatrix, drawMode);
    }

    void draw(Vec3f position, Quat4f rotation, Vec3f scale={1, 1, 1}, GLenum drawMode=0){
        if (GLES::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }

        // TODO: optimize this
    
        // scale + translation
        Mat4f ScaleTranslationMatrix;
        ScaleTranslationMatrix.setDiag(scale.x, scale.y, scale.z, 1);
        ScaleTranslationMatrix.setPos(position);
        
        // rotation
        Mat3f modelRotMatrix3;
        rotation.toMatrix(modelRotMatrix3);
        Mat4f modelRotMatrix4 = Mat4fIdentity;
        modelRotMatrix4.setRot(modelRotMatrix3);
    
        // model matrix
        Mat4f modelMatrix = modelRotMatrix4;
        modelMatrix.mmulL(ScaleTranslationMatrix);

        // MVP matrix
        Mat4f mvpMatrix = modelMatrix;
        mvpMatrix.mmulL(GLES::active_camera->viewMatrix());
        mvpMatrix.mmulL(GLES::active_camera->projectionMatrix());
    
        drawMVP(mvpMatrix, drawMode);
    }

    inline void draw2D(Vec3f pos, float scale, GLenum drawMode=0){draw2D(pos, (Vec2f){scale, scale}, drawMode);}
    void draw2D(Vec3f pos=Vec3fZero, Vec2f scale={1, 1}, GLenum drawMode=0){
        // convert from screen space ((0, 0)  to (WIDHT, HEIGHT)) to NDC ((-1, -1) to (1, 1))

        const int WIDTH = GLES::screen_size.x;
        const int HEIGHT = GLES::screen_size.y;

        pos.x = pos.x*2/WIDTH - 1;
        pos.y = pos.y*2/HEIGHT - 1;
    
        scale.x = scale.x*2/WIDTH;
        scale.y = scale.y*2/HEIGHT;

        draw2D_NDC(pos, scale, drawMode);
    }

    void draw2D_NDC(Vec3f pos=Vec3fZero, Vec2f scale={1, 1}, GLenum drawMode=0){
        Mat4f mvp;
        mvp.array[ 0] = scale.x; mvp.array[ 4] = 0;       mvp.array[ 8] = 0; mvp.array[12] = pos.x;
        mvp.array[ 1] = 0;       mvp.array[ 5] = scale.y; mvp.array[ 9] = 0; mvp.array[13] = pos.y;
        mvp.array[ 2] = 0;       mvp.array[ 6] = 0;       mvp.array[10] = 1; mvp.array[14] = pos.z;
        mvp.array[ 3] = 0;       mvp.array[ 7] = 0;       mvp.array[11] = 0; mvp.array[15] = 1;

        drawMVP(mvp, drawMode);
    }
};


#endif // _GLMesh_H_

