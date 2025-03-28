#ifndef _GLMESH_H_
#define _GLMESH_H_

#include "GLES2.h"
#include "GLTexture.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Shader.h"
#include <variant>
#include <vector>
#include <cstddef>
#include <type_traits>

template<unsigned int attrib_flags>
class GLMesh{
public:
    struct vertex{
        Vec3f position;
        std::conditional_t<attrib_flags&GLMESH_FLAG_NORMAL, Vec3f, std::monostate> normal;
        std::conditional_t<attrib_flags&GLMESH_FLAG_COLOR , Vec3f, std::monostate> color;
        std::conditional_t<attrib_flags&GLMESH_FLAG_UV    , Vec2f, std::monostate> uv;
    };

private:
    std::vector<vertex> vertices;

    GLenum usage = GL_STATIC_DRAW;
    GLuint vbo = 0;
    bool vbo_sync = true;

    inline void init_vbo(){
        if (vbo) return;
        glGenBuffers(1, &vbo);
    }

public:
    GLenum drawMode;
    Vec3f color = {1.0f, 1.0f, 1.0f};
    Shader<attrib_flags>* shader;
    GLTexture* texture = nullptr;

    GLMesh(GLenum drawMode=GL_TRIANGLES, GLenum usage=GL_STATIC_DRAW, Shader<attrib_flags>* shader=defaultShader<attrib_flags>, GLTexture* texture=nullptr): drawMode(drawMode), usage(usage), shader(shader), texture(texture)
    {
        if (attrib_flags&GLMESH_FLAG_TEX && !texture) printf("Warning: GLMesh created with GLMESH_FLAG_TEX but no texture provided!\n");
        if (attrib_flags&GLMESH_FLAG_TEX && !attrib_flags&GLMESH_FLAG_UV) printf("Warning: GLMesh created with GLMESH_FLAG_TEX but not GLMESH_FLAG_UV!\n");
    };

    void clear(){
        vertices.clear();
        vbo_sync = false;
    }

    void addVertex(Vec3f pos, Vec3f normal=Vec3fZero, Vec3f color=COLOR_WHITE, Vec2f uv=Vec2fZero){
        vertex v;
        v.position = pos;
        if constexpr (attrib_flags&GLMESH_FLAG_NORMAL) v.normal = normal;
        if constexpr (attrib_flags&GLMESH_FLAG_COLOR ) v.color  = color;
        if constexpr (attrib_flags&GLMESH_FLAG_UV    ) v.uv     = uv;
        vertices.push_back(v);

        vbo_sync = false;
    }

    void updateVertex(int i, Vec3f pos, Vec3f normal=Vec3fZero, Vec3f color=COLOR_WHITE, Vec2f uv=Vec2fZero){
        vertex v;
        v.position = pos;
        if constexpr (attrib_flags&GLMESH_FLAG_NORMAL) v.normal = normal;
        if constexpr (attrib_flags&GLMESH_FLAG_COLOR ) v.color  = color;
        if constexpr (attrib_flags&GLMESH_FLAG_UV    ) v.uv     = uv;
        vertices[i] = v;

        vbo_sync = false;
    }

    inline int vertexCount() const { return vertices.size(); }

    void bind_sync_vbo(){
        init_vbo();
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glVertexAttribPointer(SHADER_ATTRIB_POSITION, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, position));
        if constexpr (attrib_flags&GLMESH_FLAG_NORMAL) glVertexAttribPointer(SHADER_ATTRIB_NORMAL, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, normal));
        if constexpr (attrib_flags&GLMESH_FLAG_COLOR ) glVertexAttribPointer(SHADER_ATTRIB_COLOR , 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, color ));
        if constexpr (attrib_flags&GLMESH_FLAG_UV    ) glVertexAttribPointer(SHADER_ATTRIB_UV    , 2, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, uv    ));
        if (vbo_sync) return;

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), GL_STATIC_DRAW);
        vbo_sync = true;
    }

    void drawMVP(Mat4f mvp){
        bind_sync_vbo();
    
        shader->use();
        shader->setUniformMat4f("uMVPMatrix", mvp);
        shader->setUniform3f("uColor", color);
        if constexpr (attrib_flags&GLMESH_FLAG_TEX) {
            glActiveTexture(GL_TEXTURE0);
            shader->setUniform1i("uTexture", 0);
            if (texture) texture->bind();
            else printf("Warning: texture = nullptr!\n");
        }
    
        glDrawArrays(drawMode, 0, vertexCount());
    }

    void draw(Vec3f position=Vec3fZero, Quat4f rotation=Quat4fIdentity, Vec3f scale={1, 1, 1}){
        if (GLES2::active_camera == nullptr){
            printf("Warning: No active camera - skipping rendering.\n");
            return;
        }
    
        // scale
        Mat4f modelScaleMatrix = Mat4fIdentity;
        modelScaleMatrix.mul({scale.x, scale.y, scale.z, 1.0});
        
        // rotation
        Mat3f modelRotMatrix3;
        rotation.toMatrix(modelRotMatrix3);
        Mat4f modelRotMatrix4 = Mat4fIdentity;
        modelRotMatrix4.setRot(modelRotMatrix3);
    
        // translation
        Mat4f modelTranslationMatrix = Mat4fIdentity;
        modelTranslationMatrix.setPos(position);
    
        // model matrix
        Mat4f modelMatrix = Mat4fIdentity;
        modelMatrix.mmulL(modelScaleMatrix);
        modelMatrix.mmulL(modelRotMatrix4);
        modelMatrix.mmulL(modelTranslationMatrix);
    
        // MVP matrix
        Mat4f mvpMatrix = modelMatrix;
        mvpMatrix.mmulL(GLES2::active_camera->viewMatrix());
        mvpMatrix.mmulL(GLES2::active_camera->projectionMatrix());
    
        drawMVP(mvpMatrix);
    }

    ~GLMesh(){
        if (vbo) glDeleteBuffers(1, &vbo);
    }
};

using GLMesh_Texture = GLMesh<GLMESH_FLAG_UV | GLMESH_FLAG_TEX>;
using GLMesh_Normal = GLMesh<GLMESH_FLAG_NORMAL>;
using GLMesh_Color = GLMesh<GLMESH_FLAG_COLOR>;
using GLMesh_NC = GLMesh<GLMESH_FLAG_NORMAL | GLMESH_FLAG_COLOR>;

#endif // _GLMESH_H_
