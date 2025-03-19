
#ifndef SHADER_H
#define SHADER_H

#include "GLES2.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat4.h"

#define SHADER_ATTRIB_POSITION 0
#define SHADER_ATTRIB_NORMAL 1
#define SHADER_ATTRIB_COLOR 2
//#define SHADER_ATTRIB_UV 3

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

Shader* getDefaultShader();

#endif // SHADER_H
