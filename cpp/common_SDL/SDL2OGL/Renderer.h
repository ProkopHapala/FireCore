#ifndef _Renderer_H_
#define _Renderer_H_

//#include <SDL2/SDL.h>

#include <SDL2/SDL_opengles2.h>


#include <math.h>
#include <cstdlib>
#include <stdint.h>
#include <vector>

#include "Vec3.h"
#include "GLMesh.h"
#include "Mat4.h"
#include "quaternion.h"
#include "Camera.h"

#include <functional>

#ifndef __SKIP_DEFS__
#include "SDL_opengl_defs.h" // type and constants for old OpenGL
#endif

#define SHADER_ATTRIB_POSITION 0
#define SHADER_ATTRIB_NORMAL 1
#define SHADER_ATTRIB_COLOR 2

class Renderer {
    private:
        const char* vertexShaderSource = R"(
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
        const char* fragmentShaderSource = R"(
            varying vec3 fColor;
            varying vec3 fNormal;

            void main() {
                gl_FragColor = vec4(fColor, 1.0);
            }
        )";

        GLuint compileShader(GLenum shaderType, const char* source) {
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

        GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader) {
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

        GLuint program;
        GLint mvpMatrixLocation;

    public:
        Camera* active_camera = nullptr;

        Renderer(){
            // Compile shaders
            GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
            GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
            if (vertexShader == 0 || fragmentShader == 0){
                printf("Error compiling shaders\n");
                exit(-1);
                return;
            }
            program = linkProgram(vertexShader, fragmentShader);
            if (program == 0){
                printf("Error linking program\n");
                exit(-1);
                return;
            }
            glUseProgram(program);

            glEnableVertexAttribArray(SHADER_ATTRIB_POSITION);
            glVertexAttribPointer(SHADER_ATTRIB_POSITION, 3, GL_FLOAT, GL_FALSE, sizeof(GLMesh::vertex), 0);

            glEnableVertexAttribArray(SHADER_ATTRIB_NORMAL);
            glVertexAttribPointer(SHADER_ATTRIB_NORMAL, 3, GL_FLOAT, GL_FALSE, sizeof(GLMesh::vertex), (const void*)(3*sizeof(float)));

            glEnableVertexAttribArray(SHADER_ATTRIB_COLOR);
            glVertexAttribPointer(SHADER_ATTRIB_COLOR, 3, GL_FLOAT, GL_FALSE, sizeof(GLMesh::vertex), (const void*)(6*sizeof(float)));
          
            mvpMatrixLocation = glGetUniformLocation(program, "uMVPMatrix");
        }

        void drawMeshMVP(GLMesh* mesh, Mat4f mvp){
            mesh->bind_sync_vbo();

            //printf("MVP Matrix:\n");
            //mvp.print();

            glUseProgram(program);

            GLuint positionAttrib = 0;
            glEnableVertexAttribArray(positionAttrib);
            glVertexAttribPointer(positionAttrib, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), 0);

            GLuint normalAttrib = 1;
            glEnableVertexAttribArray(normalAttrib);
            glVertexAttribPointer(normalAttrib, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (const void*)(3*sizeof(float)));

            GLuint colorAttrib = 2;
            glEnableVertexAttribArray(colorAttrib);
            glVertexAttribPointer(colorAttrib, 3, GL_FLOAT, GL_FALSE, 9*sizeof(float), (const void*)(3*sizeof(float)));
            
            glUniformMatrix4fv(mvpMatrixLocation, 1, GL_FALSE, mvp.array);
            
            glDrawArrays(mesh->drawMode, 0, mesh->vertexCount());
        }

        void drawMesh(GLMesh* mesh, Vec3f position, Quat4f rotation=Quat4fIdentity, Vec3f scale={1, 1, 1});
};

class OpenGL1Renderer {
    public:
        Renderer* renderer;

        Vec3f color = {1, 1, 1};
        float color_alpha = 1;

        Vec3f normal = {0, 0, 0};

        bool begun = false;
        GLenum mode = -1;

        GLMesh* current_mesh;

        GLenum MatrixMode = GL_MODELVIEW;


        std::vector<Mat4f> mvMatStack;
        std::vector<Mat4f> projMatStack;
        std::vector<Mat4f> textureMatStack;
        std::vector<Mat4f> colorMatStack;

        std::vector<std::vector<std::function<void()>>> callLists;

        std::vector<std::function<void(void)>> current_call_list_builder;
        int current_callList = 0;
        GLenum current_callListMode = GL_COMPILE;

    public:
        OpenGL1Renderer();

        void bind_renderer(Renderer* r){
            this->renderer = r;
        }

        // ===== OpenGL 1.x functions =====
        void begin(GLenum mode);
        void end();
        void enable(GLenum cap);
        void disable(GLenum cap);
        void flush();
        void finish();

        void color3f(float r, float g, float b);
        void color3d(double r, double g, double b);
        void color4f(float r, float g, float b, float a);

        void vertex2d(double x, double y);
        void vertex3f(float x, float y, float z);
        void vertex3d(double x, double y, double z);
        
        void normal3f(float x, float y, float z);
        void normal3d(double x, double y, double z);

        void newList(GLuint list, GLenum mode);
        GLuint genLists(GLsizei range);
        void deleteLists(GLuint list, GLsizei range);
        void endList();
        void callList(GLuint list);
        
        void clear(GLbitfield mask);
        void clearColor(float r, float g, float b, float a);

        void pushMatrix();
        void popMatrix();
        void loadMatrixf(const GLfloat *m);
        void matrixMode(GLenum mode);
        void multMatrixf(const GLfloat *m);

        void depthFunc(GLenum func);
        void blendFunc(GLenum sfactor, GLenum dfactor);

        void translatef(float x, float y, float z);
        void readPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels);
        void genTextures(GLsizei n, GLuint *textures);
        void shadeModel(GLenum mode);
        void viewport(GLint x, GLint y, GLsizei width, GLsizei height);
        void bindTexture(GLenum target, GLuint texture);
        void texParameteri(GLenum target, GLenum pname, GLint param);
        void lineWidth(GLfloat width);
        void pixelStorei(GLenum pname, GLint param);
        void loadIdentity();
        void frustum(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar);
        void ortho(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar);
        void getFloatv(GLuint pname, GLfloat *params);
        void rotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z);
        void scalef(GLfloat x, GLfloat y, GLfloat z);
        void pointSize(GLfloat size);
        void scissor(GLint x, GLint y, GLsizei width, GLsizei height);
        void frontFace(GLenum mode);
        void materialfv(GLenum face, GLenum pname, const GLfloat *params);
        void lightModeli(GLenum pname, GLint param);
        void polygonMode(GLenum face, GLenum mode);
        void lightfv(GLenum light, GLenum pname, const GLfloat *params);
        void texImage2D(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid *pixels);
        void texCoord2f(GLfloat s, GLfloat t);
        void hint(GLenum target, GLenum mode);
        void lightModelfv(GLenum pname, const GLfloat *params);

        void TEST();

        const GLubyte* getString(GLenum name);

        void copyTexImage2D(GLenum target, GLint level,
            GLenum internalformat,
            GLint x, GLint y,
            GLsizei width, GLsizei height,
            GLint border);

};

extern OpenGL1Renderer opengl1renderer;




#endif // _Renderer_H_
