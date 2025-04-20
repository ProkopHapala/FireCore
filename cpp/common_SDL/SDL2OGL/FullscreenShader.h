#include "GLES.h"
#include "Shader.h"
#include "GLMesh.h"
#include "GLTexture.h"
#include "Framebuffer.h"

const static char* vertexShaderSource = R"(

attribute vec4 vPosition;

varying vec2 fUV;

void main(){
    gl_Position = vPosition;
    fUV = (vPosition.xy + vec2(1, 1)) / 2.0;
}

)";

class FullscreenShader{
private:
    Shader<GLMESH_FLAG_TEX> shader;
    GLMesh<GLMESH_FLAG_TEX> mesh;
    GLFramebuffer framebuffer;

public:
    FullscreenShader(const char* fragShaderSource) : shader(Shader<GLMESH_FLAG_TEX>(vertexShaderSource, fragShaderSource))
    {
        mesh = GLMesh<GLMESH_FLAG_TEX>(GL_TRIANGLES, GL_STATIC_DRAW, &shader, &framebuffer.colorBuffer);
        mesh.addVertex({-1, -1, -1});
        mesh.addVertex({ 3, -1, -1});
        mesh.addVertex({-1,  3, -1});
    }

    // shader will get applied to everything rendered between begin() and end()
    void begin(){
        framebuffer.begin();
    }

    void end(){
        framebuffer.end();

        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
        glDisable(GL_DEPTH_TEST);
        mesh.draw2D_NDC();
        glEnable(GL_DEPTH_TEST);
    }

    inline void pause(){ framebuffer.pause(); }
    inline void unpause(){ framebuffer.unpause(); }
};
