#ifndef _MESH_H_
#define _MESH_H_

#include "GLES2.h"
#include "Vec3.h"
#include "Shader.h"
#include <vector>
#include <cstddef>

class GLMesh {
public:
    struct vertex{
        Vec3f position;
        Vec3f normal;
        Vec3f color;
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
    Shader* shader;

    GLMesh(GLenum drawMode=GL_TRIANGLES, GLenum usage=GL_STATIC_DRAW, Shader* shader=getDefaultShader()): drawMode(drawMode), usage(usage), shader(shader) {};

    void clear(){
        vertices.clear();
        vbo_sync = false;
    }

    void addVertex(Vec3f position, Vec3f normal={0, 0, 0}, Vec3f color={1, 1, 1}){
        vertex v;
        v.position = position;
        v.normal = normal;
        v.color = color;
        vertices.push_back(v);

        vbo_sync = false;
    }

    void updateVertex(int i, Vec3f position, Vec3f normal={0, 0, 0}, Vec3f color={1, 1, 1}){
        vertex v;
        v.position = position;
        v.normal = normal;
        v.color = color;
        vertices[i] = v;

        vbo_sync = false;
    }

    inline int vertexCount() const { return vertices.size(); }

    void bind_sync_vbo(){
        init_vbo();
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glVertexAttribPointer(SHADER_ATTRIB_POSITION, 3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, position));
        glVertexAttribPointer(SHADER_ATTRIB_NORMAL,   3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, normal));
        glVertexAttribPointer(SHADER_ATTRIB_COLOR,    3, GL_FLOAT, GL_FALSE, sizeof(vertex), (void*)offsetof(vertex, color));
        if (vbo_sync) return;

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), GL_STATIC_DRAW);
        vbo_sync = true;
    }

    void drawMVP(Mat4f mvp){
        bind_sync_vbo();
    
        shader->use();
        shader->setUniformMat4f("uMVPMatrix", mvp);
        shader->setUniform3f("uColor", color);
    
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

#endif // _MESH_H_
