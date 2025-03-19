
#include "Shader.h"
#include "GLES2.h"
#include <cstring>

GLuint Shader::compileShader(GLenum shaderType, const char* source) {
    GLuint shader = glCreateShader(shaderType);
    if (shader == 0) {
        printf("Error creating shader\n");
        exit(-1);
        return 0;
    }
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);

    GLint compiled;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &compiled);
    if (!compiled) {
        GLint infoLen = 0;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &infoLen);
        if (infoLen > 1) {
            char* infoLog = (char*)malloc(sizeof(char) * infoLen);
            glGetShaderInfoLog(shader, infoLen, nullptr, infoLog);
            printf("Error compiling shader:\n %s\n", infoLog);
            free(infoLog);
            exit(-1);
        }
        glDeleteShader(shader);
        return 0;
    }
    return shader;
}

GLuint Shader::linkProgram(GLuint vertexShader, GLuint fragmentShader) {
    GLuint program = glCreateProgram();
    if (program == 0) {
        printf("Error creating program\n");
        exit(-1);
        return 0;
    }
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    glBindAttribLocation(program, SHADER_ATTRIB_POSITION, "vPosition");
    glBindAttribLocation(program, SHADER_ATTRIB_NORMAL  , "vNormal");
    glBindAttribLocation(program, SHADER_ATTRIB_COLOR   , "vColor");

    glLinkProgram(program);

    GLint linked;
    glGetProgramiv(program, GL_LINK_STATUS, &linked);
    if (!linked) {
        GLint infoLen = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLen);
        if (infoLen > 1) {
            char* infoLog = (char*)malloc(sizeof(char) * infoLen);
            glGetProgramInfoLog(program, infoLen, nullptr, infoLog);
            printf("Error linking program:\n %s\n", infoLog);
            exit(-1);
            free(infoLog);
        }
        glDeleteProgram(program);
        return 0;
    }
    return program;
}

void Shader::use(){
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
    }

    glUseProgram(programId);
}

Shader::Shader(const char* vertexShaderSource, const char* fragmentShaderSource){
    __vertexShaderSource = (char*)malloc(strlen(vertexShaderSource) + 1);
    __fragmentShaderSource = (char*)malloc(strlen(fragmentShaderSource) + 1);
    strcpy(__vertexShaderSource, vertexShaderSource);
    strcpy(__fragmentShaderSource, fragmentShaderSource);
}

Shader::~Shader(){
    if (programId)glDeleteProgram(programId);
    if (vertexShaderId)glDeleteShader(vertexShaderId);
    if (fragmentShaderId)glDeleteShader(fragmentShaderId);
}




// Vertex shader
static const char* defaultVertexShaderSource = R"(
    uniform mat4 uMVPMatrix;

    attribute vec4 vPosition;
    attribute vec3 vNormal;
    attribute vec3 vColor;
    
    varying vec3 fColor;
    varying vec3 fNormal;

    void main() {
        gl_Position = uMVPMatrix * vPosition;
        fNormal = vNormal;
        fColor = vColor;
    }
)";

// Fragment shader
static const char* defaultFragmentShaderSource = R"(
    varying vec3 fColor;
    varying vec3 fNormal;

    uniform vec3 uColor;

    void main() {
        if (fNormal == vec3(0.0, 0.0, 0.0)){
            gl_FragColor = vec4(fColor*uColor, 1.0);
        }else{
            gl_FragColor = vec4(fNormal, 1.0);
        }
    }
)";

static Shader* defaultShader = nullptr;
Shader* getDefaultShader(){
    if (defaultShader != nullptr) return defaultShader;
    defaultShader = new Shader(defaultVertexShaderSource, defaultFragmentShaderSource);
    return defaultShader;
}


