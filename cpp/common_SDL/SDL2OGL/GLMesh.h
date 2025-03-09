#ifndef _MESH_H_
#define _MESH_H_

#include "GLES2.h"
#include <vector>
#include "Renderer.h"
#include "Vec3.h"

class GLMesh {
public:
    struct vertex{
        Vec3f position;
        Vec3f normal;
        Vec3f color;
    };

private:
    std::vector<vertex> vertices;

    GLuint vbo = 0;
    bool vbo_sync = true;

    inline void init_vbo(){
        if (vbo) return;
        glGenBuffers(1, &vbo);
    }

public:
    GLenum drawMode = GL_TRIANGLES;
    Vec3f color = {1.0f, 1.0f, 1.0f};

    GLMesh(GLenum drawMode=GL_TRIANGLES): drawMode(drawMode) {};

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

    void bind_sync_vbo( Renderer* r ){
        init_vbo();
        r->bindBuffer(GL_ARRAY_BUFFER, vbo);
        if (vbo_sync) return;

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), GL_STATIC_DRAW);
        vbo_sync = true;
    }

    ~GLMesh(){
        if (vbo) glDeleteBuffers(1, &vbo);
    }
};

#endif // _MESH_H_
