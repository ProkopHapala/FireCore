
#include "Shader.h"
#include "GLES.h"
#include <cstring>

GLuint compileShader(GLenum shaderType, const char* source) {
    GLuint shader = glCreateShader(shaderType);
    if (shader == 0) {
        printf("Error creating shader\n");
        exit(-1);
        return 0;
    }
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint infoLen = 0;
    glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLen);
    if (infoLen > 1) {
        char* infoLog = (char*)malloc(sizeof(char) * infoLen);
        glGetShaderInfoLog(shader, infoLen, nullptr, infoLog);
        printf("Error/Warning compiling shader:\n %s\n", infoLog);
        free(infoLog);

        printf("Shader source:\n%s\n", source);
    }
    
    GLint compiled;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (!compiled) {
        glDeleteShader(shader);
        exit(-1);
    }

    return shader;
}

GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader) {
    GLuint program = glCreateProgram();
    if (program == 0) {
        printf("Error creating program\n");
        exit(-1);
        return 0;
    }
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    glLinkProgram(program);
    
    GLint infoLen = 0;
    glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLen);
    if (infoLen > 1) {
        char* infoLog = (char*)malloc(sizeof(char) * infoLen);
        glGetProgramInfoLog(program, infoLen, nullptr, infoLog);
        printf("Error/Warning linking program:\n %s\n", infoLog);
        free(infoLog);
    }
    
    GLint linked;
    glGetProgramiv(program, GL_LINK_STATUS, &linked);
    if (!linked) {
        glDeleteProgram(program);
        exit(-1);
    }

    return program;
}
