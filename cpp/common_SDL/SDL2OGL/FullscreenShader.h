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
    Shader<0> shader;
    GLMesh<0> mesh;
    GLFramebuffer framebuffer_back; // everything rendered between "begin()" and "end()" will be rendered to this framebuffer
    
public:
    GLFramebuffer out_framebuffer; // after "end()", the result of the shader is stored in this framebuffer
    
    FullscreenShader(const char* fragShaderSource) : shader(Shader<0>(vertexShaderSource, fragShaderSource))
    {
        printf("%s\n", fragShaderSource);
        mesh = GLMesh<0>(GL_TRIANGLES, GL_STATIC_DRAW, &shader);
        mesh.addVertex({-1, -1, -1});
        mesh.addVertex({ 3, -1, -1});
        mesh.addVertex({-1,  3, -1});
    }

    void begin(){
        framebuffer_back.begin();
    }

    void end(){
        framebuffer_back.end();
        out_framebuffer.begin();

        glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, framebuffer_back.colorBuffer.getHandle());
        glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, framebuffer_back.depthBuffer.getHandle());
        shader.setUniform1i("uTexture", 0);
        shader.setUniform1i("uDepth"  , 1);
        
        glDisable(GL_DEPTH_TEST);
        mesh.draw2D_NDC();

        out_framebuffer.end();
    }

    inline void pause(){ framebuffer_back.pause(); }
    inline void unpause(){ framebuffer_back.unpause(); }
};
