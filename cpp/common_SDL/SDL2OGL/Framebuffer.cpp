#include "Framebuffer.h"
#include "GLES.h"
#include "GLMesh.h"
#include "GLTexture.h"
#include "Shader.h"
#include <GLES3/gl3.h>

const static char* vertexShaderSource = R"(
#version 300 es

in mediump vec4 vPosition;
out mediump vec2 fUV;

void main(){
    gl_Position = vPosition;
    fUV = (vPosition.xy + vec2(1, 1)) / 2.0;
}

)";

static const char* fragmentShaderSource = R"(
#version 300 es

uniform sampler2D uTexture1;
uniform sampler2D uDepth1;

uniform sampler2D uTexture2;
uniform sampler2D uDepth2;

in mediump vec2 fUV;
layout(location=0) out mediump vec4 FragColor;

void main(){
    highp float depth1 = texture(uDepth1, fUV).x;
    highp float depth2 = texture(uDepth2, fUV).x;
    if (depth1 < depth2){
        FragColor = texture(uTexture1, fUV);
        gl_FragDepth = depth1;
    }else{
        FragColor = texture(uTexture2, fUV);
        gl_FragDepth = depth2;
    }
}

)";

static Shader<0> shader = Shader<0>(vertexShaderSource, fragmentShaderSource);

static inline GLMesh<0>makeMergeMesh(){
    GLMesh<0> m = GLMesh<0>(GL_TRIANGLES, GL_STATIC_DRAW, &shader);
    m.addVertex({-1, -1, -1});
    m.addVertex({ 3, -1, -1});
    m.addVertex({-1,  3, -1});
    return m;
}
static GLMesh<0> mergeMesh = makeMergeMesh();


void mergeFramebuffers(GLFramebuffer& fb1, GLFramebuffer& fb2){
    shader.use();
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, fb1.colorBuffer.getHandle());
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, fb2.colorBuffer.getHandle());
    glActiveTexture(GL_TEXTURE2); glBindTexture(GL_TEXTURE_2D, fb1.depthBuffer.getHandle());
    glActiveTexture(GL_TEXTURE3); glBindTexture(GL_TEXTURE_2D, fb2.depthBuffer.getHandle());
    shader.setUniform1i("uTexture1", 0);
    shader.setUniform1i("uTexture2", 1);
    shader.setUniform1i("uDepth1"  , 2);
    shader.setUniform1i("uDepth2"  , 3);
    
    glDepthMask(GL_TRUE);
    glDisable(GL_DEPTH_TEST);
    mergeMesh.draw2D_NDC();
}

