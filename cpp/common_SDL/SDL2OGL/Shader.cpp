
#include "Shader.h"
#include "GLES.h"
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
void Shader<attrib_flags>::ensure_handle(){
    if (programId) return;

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
}

template<unsigned int attrib_flags>
void Shader<attrib_flags>::use(){
    ensure_handle();
    glUseProgram(programId);

    // update uniforms
    while (update_uniforms){
        int i = __builtin_ctz(update_uniforms); // count trailing zeros
        uniform& u = uniforms[i];

        switch (u.type) {
            case uniform::f1:   glUniform1f(i, u.data.f1); break;
            case uniform::f2:   glUniform2f(i, u.data.f2.x, u.data.f2.y); break;
            case uniform::f3:   glUniform3f(i, u.data.f3.x, u.data.f3.y, u.data.f3.z); break;
            case uniform::f4:   glUniform4f(i, u.data.f4.x, u.data.f4.y, u.data.f4.z, u.data.f4.w); break;
            case uniform::i1:   glUniform1i(i, u.data.i1); break;
            case uniform::i2:   glUniform2i(i, u.data.i2.x, u.data.i2.y); break;
            case uniform::i3:   glUniform3i(i, u.data.i3.x, u.data.i3.y, u.data.i3.z); break;
            case uniform::i4:   glUniform4i(i, u.data.i4.x, u.data.i4.y, u.data.i4.z, u.data.i4.w); break;
            case uniform::ui1:  glUniform1ui(i, u.data.ui1); break;
            case uniform::ui2:  glUniform2ui(i, u.data.ui2.x, u.data.ui2.y); break;
            case uniform::ui3:  glUniform3ui(i, u.data.ui3.x, u.data.ui3.y, u.data.ui3.z); break;
            case uniform::ui4:  glUniform4ui(i, u.data.ui4.x, u.data.ui4.y, u.data.ui4.z, u.data.ui4.w); break;
            case uniform::m3:   glUniformMatrix3fv(i, 1, GL_FALSE, u.data.m3.array); break;
            case uniform::m4:   glUniformMatrix4fv(i, 1, GL_FALSE, u.data.m4.array); break;
        }

        update_uniforms &= ~(1 << i);
    }
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
