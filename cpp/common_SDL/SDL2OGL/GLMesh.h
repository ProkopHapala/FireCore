#ifndef _MESH_H_
#define _MESH_H_

#include <GLES2/gl2.h>
#include <vector>
#include <SDL2/SDL_opengles2.h>
#include "Vec3.h"

class GLMesh {
public:
    struct vertex{
        Vec3f position;
        Vec3f normal;
        Vec3f color;
    };

    std::vector<vertex> vertices;

    //std::vector<unsigned int> indices;
    int drawMode = GL_TRIANGLES;

    Vec3f color = {1.0f, 1.0f, 1.0f};
    
    GLuint vbo = 0;
    bool vbo_sync = true;

    //GLuint ibo = 0;

    GLMesh(int drawMode): drawMode(drawMode) {
        glGenBuffers(1, &vbo);
    };

    void addVertex(Vec3f position, Vec3f normal={0, 0, 0}, Vec3f color={1, 1, 1}){
        vertex v;
        v.position = position;
        v.normal = normal;
        v.color = color;
        vertices.push_back(v);

        vbo_sync = false;
    }

    void bind_sync_vbo(){
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        if (vbo_sync) return;

        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(vertex), vertices.data(), GL_STATIC_DRAW);
        vbo_sync = true;
    }

};

#endif // _MESH_H_
