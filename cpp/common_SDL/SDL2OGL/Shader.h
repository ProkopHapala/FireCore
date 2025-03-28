
#ifndef SHADER_H
#define SHADER_H

#include "GLES2.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat4.h"
#include <string>

#define SHADER_ATTRIB_POSITION 0
#define SHADER_ATTRIB_NORMAL 1
#define SHADER_ATTRIB_COLOR 2
#define SHADER_ATTRIB_UV 3

template<unsigned int attrib_flags>
class Shader {
private:
    GLuint programId        = 0;
    GLuint vertexShaderId   = 0;
    GLuint fragmentShaderId = 0;

    char* __vertexShaderSource   = nullptr;
    char* __fragmentShaderSource = nullptr;

public:
    Shader(const char* vertexShaderSource, const char* fragmentShaderSource);
    ~Shader();
    
    inline GLuint getProgramId() const { return programId; }

    void use();
    
    // Uniform setters
    GLint getUniformLocation(const char* name) {
        GLint loc = glGetUniformLocation(programId, name);
        if (loc == -1) {
            printf("Uniform '%s' not found\n", name);
            exit(-1);
        }
        return loc;
    }

    void setUniform1f(const char* name, float value)                        { glUniform1f(getUniformLocation(name), value); } // TODO: don't use glGetUniformLocation
    void setUniform2f(const char* name, float x, float y)                   { glUniform2f(getUniformLocation(name), x, y); }
    void setUniform3f(const char* name, float x, float y, float z)          { glUniform3f(getUniformLocation(name), x, y, z); }
    void setUniform4f(const char* name, float x, float y, float z, float w) { glUniform4f(getUniformLocation(name), x, y, z, w); }
    void setUniformMat4f(const char* name, const float* matrix)      { glUniformMatrix4fv(getUniformLocation(name), 1, GL_FALSE, matrix); }

    inline void setUniform2f(const char* name, Vec2f x) {setUniform2f(name, x.x, x.y);}
    inline void setUniform3f(const char* name, Vec3f x) {setUniform3f(name, x.x, x.y, x.z);}
    inline void setUniformMat4f(const char* name, Mat4f matrix) {setUniformMat4f(name, matrix.array);}

private:
    GLuint compileShader(GLenum shaderType, const char* source);
    GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader);
};







// default shaders
constexpr const std::string buildDefaultVertexShaderSource(unsigned int attrib_flags){
    std::string source = std::string("");
    //source += "#version 110\n";

    // uniforms
    source += "uniform mat4 uMVPMatrix;\n";

    // attributes
    source += "attribute vec4 vPosition;\n";
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "attribute vec3 vNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "attribute vec3 vColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "attribute vec3 vUV    ;\n";

    // varyings
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "varying vec3 fNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "varying vec3 fColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "varying vec3 fUV    ;\n";

    // void main()
    source += "void main() {\n";
    source += "gl_Position = uMVPMatrix * vPosition;\n";
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "fNormal = vNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "fColor = vColor;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "fUV    = vUV;\n";
    source += "}\n";

    return source;
}

constexpr const std::string buildDefaultFragmentShaderSource(unsigned int attrib_flags){
    std::string source = std::string("");
    //source += "#version 110\n";

    // uniforms
    source += "uniform vec3 uColor;\n";

    // varyings
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "varying vec3 fNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "varying vec3 fColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "varying vec3 fUV    ;\n";

    // void main()
    source += "void main() {\n";
    if      (attrib_flags & GLMESH_FLAG_NORMAL) source += "gl_FragColor = vec4(fNormal*uColor, 1.0);\n";
    else if (attrib_flags & GLMESH_FLAG_COLOR ) source += "gl_FragColor = vec4(fColor*uColor, 1.0);\n";
    else                                        source += "gl_FragColor = vec4(uColor, 1.0);\n";
    source += "}\n";

    return source;
}

template <unsigned int attrib_flags>
Shader<attrib_flags>* defaultShader = new Shader<attrib_flags>(buildDefaultVertexShaderSource(attrib_flags).c_str(), buildDefaultFragmentShaderSource(attrib_flags).c_str());

#endif // SHADER_H
