
#include "Shader.h"
#include "GLES2.h"
#include <cstring>

template<unsigned int attrib_flags>
GLuint Shader<attrib_flags>::compileShader(GLenum shaderType, const char* source) {
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
    }
    
    GLint compiled;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (!compiled) {
        glDeleteShader(shader);
        exit(-1);
    }

    return shader;
}

template<unsigned int attrib_flags>
GLuint Shader<attrib_flags>::linkProgram(GLuint vertexShader, GLuint fragmentShader) {
    GLuint program = glCreateProgram();
    if (program == 0) {
        printf("Error creating program\n");
        exit(-1);
        return 0;
    }
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    glBindAttribLocation(program, SHADER_ATTRIB_POSITION, "vPosition");
    if constexpr (attrib_flags & GLMESH_FLAG_NORMAL) glBindAttribLocation(program, SHADER_ATTRIB_NORMAL  , "vNormal");
    if constexpr (attrib_flags & GLMESH_FLAG_COLOR ) glBindAttribLocation(program, SHADER_ATTRIB_COLOR   , "vColor" );
    if constexpr (attrib_flags & GLMESH_FLAG_UV    ) glBindAttribLocation(program, SHADER_ATTRIB_UV      , "vUV"    );

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

template<unsigned int attrib_flags>
void Shader<attrib_flags>::use(){
    if (!programId) {
        vertexShaderId   = compileShader(GL_VERTEX_SHADER,   __vertexShaderSource  );
        fragmentShaderId = compileShader(GL_FRAGMENT_SHADER, __fragmentShaderSource);
        programId = linkProgram(vertexShaderId, fragmentShaderId);

        free(__vertexShaderSource  ); __vertexShaderSource   = nullptr;
        free(__fragmentShaderSource); __fragmentShaderSource = nullptr;

        // TODO: this only needs to happen once, not every compilation
        glEnableVertexAttribArray(SHADER_ATTRIB_POSITION);
        glEnableVertexAttribArray(SHADER_ATTRIB_NORMAL);
        glEnableVertexAttribArray(SHADER_ATTRIB_COLOR);
        glEnableVertexAttribArray(SHADER_ATTRIB_UV);

        uMVPloc     = glGetUniformLocation(programId, "uMVPMatrix");
        uColorloc   = glGetUniformLocation(programId, "uColor"    );
        uTextureloc = glGetUniformLocation(programId, "uTexture"  );
    }

    glUseProgram(programId);
}

template<unsigned int attrib_flags>
Shader<attrib_flags>::Shader(const char* vertexShaderSource, const char* fragmentShaderSource){
    __vertexShaderSource = (char*)malloc(strlen(vertexShaderSource) + 1);
    __fragmentShaderSource = (char*)malloc(strlen(fragmentShaderSource) + 1);
    strcpy(__vertexShaderSource, vertexShaderSource);
    strcpy(__fragmentShaderSource, fragmentShaderSource);
}

template<unsigned int attrib_flags>
Shader<attrib_flags>::~Shader(){
    if (programId)glDeleteProgram(programId);
    if (vertexShaderId)glDeleteShader(vertexShaderId);
    if (fragmentShaderId)glDeleteShader(fragmentShaderId);
}

// TODO: is there a cleaner way to do this?
template class Shader< 0>;
template class Shader< 1>;
template class Shader< 2>;
template class Shader< 3>;
template class Shader< 4>;
template class Shader< 5>;
template class Shader< 6>;
template class Shader< 7>;
template class Shader< 8>;
template class Shader< 9>;
template class Shader<10>;
template class Shader<11>;
template class Shader<12>;
template class Shader<13>;
template class Shader<14>;
template class Shader<15>;
