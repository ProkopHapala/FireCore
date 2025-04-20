#ifndef _GLMESH_H_
#define _GLMESH_H_

#include "GLES.h"
#include "GLTexture.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Shader.h"
#include "Camera.h"
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

    inline void bind_vbo(){
        init_vbo();

        if (GLES::currentGL_ARRAY_BUFFER == vbo) return;

        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glVertexAttribPointer(SHADER_ATTRIB_POSITION, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, position));
        if constexpr (attrib_flags&GLMESH_FLAG_NORMAL) glVertexAttribPointer(SHADER_ATTRIB_NORMAL, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, normal));
        if constexpr (attrib_flags&GLMESH_FLAG_COLOR ) glVertexAttribPointer(SHADER_ATTRIB_COLOR , 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, color ));
        if constexpr (attrib_flags&GLMESH_FLAG_UV    ) glVertexAttribPointer(SHADER_ATTRIB_UV    , 2, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, uv    ));
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

    void addVertex_strip(Vec3f pos, Vec3f normal=Vec3fZero, Vec3f color=COLOR_WHITE, Vec2f uv=Vec2fZero){
        vertex v;
        v.position = pos;
        if constexpr (attrib_flags&GLMESH_FLAG_NORMAL) v.normal = normal;
        if constexpr (attrib_flags&GLMESH_FLAG_COLOR ) v.color  = color;
        if constexpr (attrib_flags&GLMESH_FLAG_UV    ) v.uv     = uv;

        if (vertexCount() >= 3){
            vertex v1 = vertices[vertices.size()-1];
            vertex v2 = vertices[vertices.size()-2];

            vertices.push_back(v2);
            vertices.push_back(v1);
        }

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

    inline void bind_sync_vbo(){
        bind_vbo();

        if (vbo_sync) return;

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), usage);
        vbo_sync = true;
    }

    void drawMVP(Mat4f mvp){
        bind_sync_vbo();
    
        shader->use();
        shader->setuMVPMatrix(mvp);
        shader->setuColor(color);
        if constexpr (attrib_flags&GLMESH_FLAG_TEX) {
            glActiveTexture(GL_TEXTURE0);
            shader->setuTexture(0);
            if (texture) texture->bind();
            else printf("Warning: texture = nullptr!\n");
        }
    
        glDrawArrays(drawMode, 0, vertexCount());
    }

    inline void draw(Vec3f position, float scale){draw(position, (Vec3f){scale, scale, scale});}
    void draw(Vec3f position=Vec3fZero, Vec3f scale={1, 1, 1}){
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
    
        drawMVP(mvpMatrix);
    }

    void draw(Vec3f position, Quat4f rotation, Vec3f scale={1, 1, 1}){
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
    
        drawMVP(mvpMatrix);
    }

    void draw2D(Vec3f pos=Vec3fZero, Vec2f scale={1, 1}){
        // convert from screen space ((0, 0)  to (WIDHT, HEIGHT)) to NDC ((-1, -1) to (1, 1))

        const int WIDTH = GLES::screen_size.x;
        const int HEIGHT = GLES::screen_size.y;

        pos.x = pos.x*2/WIDTH - 1;
        pos.y = pos.y*2/HEIGHT - 1;
    
        scale.x = scale.x*2/WIDTH;
        scale.y = scale.y*2/HEIGHT;

        draw2D_NDC(pos, scale);
    }

    void draw2D_NDC(Vec3f pos=Vec3fZero, Vec2f scale={1, 1}){
        Mat4f mvp;
        mvp.array[ 0] = scale.x; mvp.array[ 4] = 0;       mvp.array[ 8] = 0; mvp.array[12] = pos.x;
        mvp.array[ 1] = 0;       mvp.array[ 5] = scale.y; mvp.array[ 9] = 0; mvp.array[13] = pos.y;
        mvp.array[ 2] = 0;       mvp.array[ 6] = 0;       mvp.array[10] = 1; mvp.array[14] = pos.z;
        mvp.array[ 3] = 0;       mvp.array[ 7] = 0;       mvp.array[11] = 0; mvp.array[15] = 1;

        //glDisable(GL_DEPTH_TEST);
        drawMVP(mvp);
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
