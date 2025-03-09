
#include "Renderer.h"
#include "Camera.h"
#include "GLMesh.h"
#include "Mat4.h"
#include <GLES2/gl2.h>
#include <cstdio>
#include <functional>
#include <vector>


// Vertex shader
static const char* vertexShaderSource = R"(
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
static const char* fragmentShaderSource = R"(
    varying vec3 fColor;
    varying vec3 fNormal;

    void main() {
        if (fColor == vec3(1.0, 1.0, 1.0)){
            gl_FragColor = vec4(fNormal, 1.0);
        }else{
            gl_FragColor = vec4(fColor, 1.0);
        }
    }
)";

GLuint Renderer::compileShader(GLenum shaderType, const char* source) {
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

GLuint Renderer::linkProgram(GLuint vertexShader, GLuint fragmentShader) {
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

Renderer::Renderer(){
    // Compile shaders
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
    if (vertexShader == 0 || fragmentShader == 0){
        printf("Error compiling shaders\n");
        exit(-1);
        return;
    }
    defualtProgram = linkProgram(vertexShader, fragmentShader);
    if (defualtProgram == 0){
        printf("Error linking program\n");
        exit(-1);
        return;
    }

    glEnableVertexAttribArray(SHADER_ATTRIB_POSITION);
    glEnableVertexAttribArray(SHADER_ATTRIB_NORMAL);
    glEnableVertexAttribArray(SHADER_ATTRIB_COLOR);

    loadProgram(defualtProgram);
}

void Renderer::loadProgram(GLuint program){
    if (program == current_program) return;
    glUseProgram(program);

    mvpMatrixLocation = glGetUniformLocation(program, "uMVPMatrix");
    current_program = program;
}

void Renderer::bindBuffer(GLenum target, GLuint buffer){
    if (target != GL_ARRAY_BUFFER && target != GL_ELEMENT_ARRAY_BUFFER){
        printf("ERROR: bindBuffer(): invalid target\n");
        exit(-1);
        return;
    }
    switch (target){
        case GL_ARRAY_BUFFER:
            if (buffer == current_gl_array_buffer) return;
            current_gl_array_buffer = buffer;
            break;
        case GL_ELEMENT_ARRAY_BUFFER:
            if (buffer == current_gl_element_array_buffer) return;
            current_gl_element_array_buffer = buffer;
            break;
    }

    glBindBuffer(target, buffer);
    glVertexAttribPointer(SHADER_ATTRIB_POSITION, 3, GL_FLOAT, GL_FALSE, sizeof(GLMesh::vertex), (void*)offsetof(GLMesh::vertex, position));
    glVertexAttribPointer(SHADER_ATTRIB_NORMAL,   3, GL_FLOAT, GL_FALSE, sizeof(GLMesh::vertex), (void*)offsetof(GLMesh::vertex, normal));
    glVertexAttribPointer(SHADER_ATTRIB_COLOR,    3, GL_FLOAT, GL_FALSE, sizeof(GLMesh::vertex), (void*)offsetof(GLMesh::vertex, color));
}

void Renderer::drawMeshMVP(GLMesh* mesh, Mat4f mvp){
    mesh->bind_sync_vbo(this);
    
    loadProgram(defualtProgram);
    glUniformMatrix4fv(mvpMatrixLocation, 1, GL_FALSE, mvp.array);

    glDrawArrays(mesh->drawMode, 0, mesh->vertexCount());
}







// the global variable for the renderer
OpenGL1Renderer opengl1renderer = OpenGL1Renderer();


// OpenGL 1.x functions
OpenGL1Renderer::OpenGL1Renderer(){
    Mat4f m; m.setOne();
    mvMatStack.push_back(m);
    projMatStack.push_back(m);
    textureMatStack.push_back(m);
    colorMatStack.push_back(m);
    callLists.push_back(std::vector<std::function<void()>>());
}


void OpenGL1Renderer::begin(GLenum m){
    if (current_callList != 0){
        current_call_list_builder.push_back(std::function<void()>([this, m](){begin(m);}));
        return;
    }

    if (begun) {
        printf("Error: OpenGL1Renderer::begin() called while already begun.\n");
        *(int*)0 = 15;
        return;
    }
    mode = m; 
    begun = true;

    current_mesh = new GLMesh(mode);
}
void OpenGL1Renderer::end(){
    if (current_callList != 0){
        current_call_list_builder.push_back(std::function<void()>([this](){end();}));
        return;
    }

    if (!begun) {
        printf("Error: OpenGL1Renderer::end() called while not begun.\n");
        return;
    }
    begun = false;

    Mat4f mvpMatrix = mvMatStack.back();
    mvpMatrix.mmulL(projMatStack.back());

    renderer->drawMeshMVP(current_mesh, mvpMatrix);
}

void OpenGL1Renderer::normal3f(float x, float y, float z){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, x, y, z](){normal3f(x, y, z);});
        return;
    }
    normal = {x, y, z};
}
void OpenGL1Renderer::normal3d(double x, double y, double z){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, x, y, z](){normal3d(x, y, z);});
        return;
    }
    normal = {x, y, z};
}

void OpenGL1Renderer::color3d(double r, double g, double b){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, r, g, b](){color3d(r, g, b);});
        return;
    }
    color = {r, g, b}; color_alpha = 1;
}
void OpenGL1Renderer::color3f(float r, float g, float b){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, r, g, b](){color3f(r, g, b);});
        return;
    }
    color = {r, g, b}; color_alpha = 1;
}
void OpenGL1Renderer::color4f(float r, float g, float b, float a){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, r, g, b, a](){color4f(r, g, b, a);});
        return;
    }
    color = {r, g, b}; color_alpha = a;
}

void OpenGL1Renderer::vertex2d(double x, double y){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, x, y](){vertex2d(x, y);});
        return;
    }
    current_mesh->addVertex({x, y, 0}, normal, color);
}
void OpenGL1Renderer::vertex3f(float x, float y, float z)       {
    if (current_callList != 0){
        current_call_list_builder.push_back([this, x, y, z](){vertex3f(x, y, z);});
        return;
    }
    current_mesh->addVertex({x, y, z}, normal, color); }
void OpenGL1Renderer::vertex3d(double x, double y, double z)    {
    if (current_callList != 0){
        current_call_list_builder.push_back([this, x, y, z](){vertex3d(x, y, z);});
        return;
    }
    current_mesh->addVertex({x, y, z}, normal, color); }


void OpenGL1Renderer::flush()           { glFlush(); }
void OpenGL1Renderer::finish()          { glFinish(); }

void OpenGL1Renderer::clear(GLbitfield mask)    {
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. clear");return;}
    glClear(mask); }
void OpenGL1Renderer::enable(GLenum cap)        {
    if (current_callList != 0){
        current_call_list_builder.push_back([this, cap](){enable(cap);});
        return;
    }
    glEnable(cap); }
void OpenGL1Renderer::disable(GLenum cap)       {
    if (current_callList != 0){
        current_call_list_builder.push_back([this, cap](){disable(cap);});
        return;
    }
    glDisable(cap); }
void OpenGL1Renderer::depthFunc(GLenum func)    {
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. depthfunc ");return;}
    glDepthFunc(func); }
void OpenGL1Renderer::shadeModel(GLenum mode)   { /*glShadeModel(mode);*/ }
void OpenGL1Renderer::lineWidth(GLfloat width)  {
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. linewidth ");return;}
    glLineWidth(width); }
void OpenGL1Renderer::pointSize(GLfloat size)   { /*glPointSize(size);*/ }
void OpenGL1Renderer::frontFace(GLenum mode)    {
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. frontface ");return;}
    glFrontFace(mode); }


void OpenGL1Renderer::translatef(float x, float y, float z)     {
    if (current_callList != 0){
        current_call_list_builder.push_back([this, x, y, z](){translatef(x, y, z);});
        return;
    }
    Mat4f m; m.setOne(); m.setPos({x, y, z}); multMatrixf(m.array); }
void OpenGL1Renderer::scalef(GLfloat x, GLfloat y, GLfloat z)   {
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. scalef ");return;}
    Mat4f m; m.setOne(); m.xx=x; m.yy=y; m.zz=z; multMatrixf(m.array); }
void OpenGL1Renderer::rotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z){
    if (current_callList != 0){
        current_call_list_builder.push_back([this, angle, x, y, z](){rotatef(angle, x, y, z);});
        return;
    }
    
    angle = angle * M_PI/180;
    float c = cos(angle);
    float s = sin(angle);
    
    // normalise
    float r = sqrt(x*x + y*y + z*z);
    x /= r;
    y /= r;
    z /= r;
    
    Mat4f m; m.setOne();
    m.xx = x*x*(1-c) + c;
    m.yy = y*y*(1-c) + c;
    m.zz = z*z*(1-c) + c;
    
    m.xy = x*y*(1-c) + z*s;
    m.xz = x*z*(1-c) - y*s;
    
    m.yx = y*x*(1-c) - z*s;
    m.yz = y*z*(1-c) + x*s;

    m.zx = z*x*(1-c) + y*s;
    m.zy = z*y*(1-c) - x*s;
    
    multMatrixf(m.array);
}


void OpenGL1Renderer::loadIdentity(){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. loadIdentity ");return;}
    Mat4f m; m.setOne(); loadMatrixf(m.array); }
void OpenGL1Renderer::pushMatrix(){
    if (current_callList != 0){
        current_call_list_builder.push_back([this](){pushMatrix();});
        return;
    }

    switch (MatrixMode){
        case GL_MODELVIEW:
            mvMatStack.push_back(mvMatStack.back());
            break;
        case GL_PROJECTION:
            projMatStack.push_back(projMatStack.back());
            break;
        case GL_TEXTURE:
            textureMatStack.push_back(textureMatStack.back());
            break;
        case GL_COLOR:
            colorMatStack.push_back(colorMatStack.back());
            break;
        default:
        printf("Error: OpenGL1Renderer::pushMatrix() called with invalid MatrixMode.\n");
    }
}
void OpenGL1Renderer::popMatrix(){
    if (current_callList != 0){
        current_call_list_builder.push_back([this](){popMatrix();});
        return;
    }

    switch (MatrixMode) {
        case GL_MODELVIEW:
            mvMatStack.pop_back();
            break;
        case GL_PROJECTION:
            projMatStack.pop_back();
            break;
        case GL_TEXTURE:
            textureMatStack.pop_back();
            break;
        case GL_COLOR:
            colorMatStack.pop_back();
            break;
        default:
            printf("Error: OpenGL1Renderer::popMatrix() called with invalid MatrixMode.\n");
    }
}
void OpenGL1Renderer::loadMatrixf(const GLfloat *mat){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. loadMatrixf ");return;}

    Mat4f m; for (int i=0; i<16; i++) m.array[i] = mat[i];
    switch (MatrixMode) {
        case GL_MODELVIEW:
            mvMatStack.back() = m;
            break;
        case GL_PROJECTION:
            projMatStack.back() = m;
            break;
        case GL_TEXTURE:
            textureMatStack.back() = m;
            break;
        case GL_COLOR:
            colorMatStack.back() = m;
            break;
        default:
            printf("Error: OpenGL1Renderer::loadMatrixf() called with invalid MatrixMode.\n");
        }
    }
void OpenGL1Renderer::matrixMode(GLenum mode){ if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. matrixMode ");return;} MatrixMode = mode; }
void OpenGL1Renderer::multMatrixf(const GLfloat *mat){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. multMatrixf");return;}
    Mat4f m; for (int i=0; i<16; i++) m.array[i] = mat[i];

    switch (MatrixMode) {
        case GL_MODELVIEW:
        mvMatStack.back().mmulR(m);
        break;
        case GL_PROJECTION:
        projMatStack.back().mmulR(m);
        break;
        case GL_TEXTURE:
        textureMatStack.back().mmulR(m);
        break;
        case GL_COLOR:
        colorMatStack.back().mmulR(m);
        break;
        default:
        printf("Error: OpenGL1Renderer::multMatrixf() called with invalid MatrixMode.\n");
    }
}
void OpenGL1Renderer::frustum(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. frustum ");return;}

    Mat4f m; m.setOne();
    m.xx = 2*zNear/(right-left);
    m.yy = 2*zNear/(top-bottom);
    m.zz = -(zFar+zNear)/(zFar-zNear);
    m.zx = (right+left)/(right-left);
    m.zy = (top+bottom)/(top-bottom);
    m.zw = -1;
    m.wz = -2*(zFar*zNear)/(zFar-zNear);
    m.ww = 0;
    
    multMatrixf(m.array);
}
void OpenGL1Renderer::ortho(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zmin, GLdouble zmax){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. ortho ");return;}

    Mat4f m;

    m.array[ 0] = 2/(right-left);  m.array[ 4] = 0;               m.array[ 8] = 0;               m.array[12] = -(right+left)/(right-left);
    m.array[ 1] = 0;               m.array[ 5] = 2/(top-bottom);  m.array[ 9] = 0;               m.array[13] = -(top+bottom)/(top-bottom);
    m.array[ 2] = 0;               m.array[ 6] = 0;               m.array[10] = -2/(zmax-zmin);  m.array[14] = -(zmax+zmin) /(zmax-zmin) ;
    m.array[ 3] = 0;               m.array[ 7] = 0;               m.array[11] = 0;               m.array[15] = 1;
    
    multMatrixf(m.array);
}




GLuint OpenGL1Renderer::genLists(GLsizei range){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. genLists ");return 0;}
    for (int i=0; i<range; i++) {
        callLists.push_back(std::vector<std::function<void(void)>>());
    }
    printf("======================================= gen List: %i\n", callLists.size() - range);
    return callLists.size() - range;
}
void OpenGL1Renderer::newList(GLuint list, GLenum mode){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. NewList ");return;}
    if (list > callLists.size()) {
        printf("Error: OpenGL1Renderer::newList() called with invalid list.\n");
        return;
    }
    if (mode == GL_COMPILE_AND_EXECUTE){
        printf("Error: OpenGL1Renderer::newList() called with mode GL_COMPILE_EXECUTE - not supported.\n");
        return;
    }
    current_callList = list;
    current_callListMode = mode;
    current_call_list_builder = std::vector<std::function<void(void)>>();

    printf("======================================= newList: %i\n", list);
    //*(int*)0 = 15;
}
void OpenGL1Renderer::endList(){
    printf("======================================= endList: %i\n", current_callList);
    callLists[current_callList] = current_call_list_builder;
    current_callList = 0;
}
void OpenGL1Renderer::deleteLists(GLuint list, GLsizei range){ /*glDeleteLists(list, range);*/ }
void OpenGL1Renderer::callList(GLuint list){
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList. 'CallList' ");return;}
    if (list > callLists.size()) {
        printf("Error: OpenGL1Renderer::callList() called with invalid list.\n");
        return;
    }
    //return; // call lists are *very* slow
    for (int i=0; i<callLists[list].size(); i++) {
        callLists[list][i]();
    }
}




void OpenGL1Renderer::blendFunc(GLenum sfactor, GLenum dfactor) {
    if (current_callList != 0){
        current_call_list_builder.push_back([this, sfactor, dfactor](){blendFunc(sfactor, dfactor);});
        return;
    }
    glBlendFunc(sfactor,dfactor); }
void OpenGL1Renderer::pixelStorei(GLenum pname, GLint param)    { /*glPixelStorei(pname, param);*/ }
void OpenGL1Renderer::genTextures(GLsizei n, GLuint *textures)  { /*glGenTextures(n, textures);*/ }

void OpenGL1Renderer::getFloatv(GLuint pname, GLfloat* params)  {
    if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList.");return;}
    switch (pname) {
        case GL_MODELVIEW_MATRIX:
            for (int i=0; i<16; i++) params[i] = mvMatStack.back().array[i];
            break;
        case GL_PROJECTION_MATRIX:
            for (int i=0; i<16; i++) params[i] = projMatStack.back().array[i];
            break;
        default:
            printf("Warning: OpenGL1Renderer::getFloatv() called with unknown pname.\n");
    }
}



void OpenGL1Renderer::lightModeli(GLenum pname, GLint param)    { /*glLightModeli(pname, param);*/ }
void OpenGL1Renderer::polygonMode(GLenum face, GLenum mode)     { /*glPolygonMode(face, mode);*/ }
void OpenGL1Renderer::texCoord2f(GLfloat s, GLfloat t)          { /*glTexCoord2f(s,t);*/ }
void OpenGL1Renderer::hint(GLenum target, GLenum mode)          { if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList.");return;}glHint(target, mode); }
const GLubyte* OpenGL1Renderer::getString(GLenum name)          { return glGetString(name); }

void OpenGL1Renderer::viewport(GLint x, GLint y, GLsizei width, GLsizei height) { if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList.");return;}glViewport(x,y,width,height); }
void OpenGL1Renderer::bindTexture(GLenum target, GLuint texture)                { /*glBindTexture(target,texture);*/ }
void OpenGL1Renderer::clearColor(float r, float g, float b, float a)            { if (current_callList != 0){printf("Error: Unsuported operation while compiling CallList.");return;}glClearColor(r,g,b,a); }
void OpenGL1Renderer::texParameteri(GLenum target, GLenum pname, GLint param)   { /*glTexParameteri(target,pname,param);*/ }
void OpenGL1Renderer::scissor(GLint x, GLint y, GLsizei width, GLsizei height)  { /*glScissor(x, y, width, height);*/ }
void OpenGL1Renderer::lightModelfv(GLenum pname, const GLfloat *params)         { /*glLightModelfv(pname, params);*/ }
void OpenGL1Renderer::lightfv(GLenum light, GLenum pname, const GLfloat *params){ /*glLightfv(light, pname, params);*/ }

void OpenGL1Renderer::readPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels)   { /*glReadPixels(x,y,width,height,format,type,pixels); */ }
void OpenGL1Renderer::materialfv(GLenum face, GLenum pname, const GLfloat *params)                                              { /*glMaterialfv(face, pname, params);*/ }

void OpenGL1Renderer::texImage2D(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid *pixels){ /*glTexImage2D(target, level, internalformat, width, height, border, format, type, pixels);*/ }

void OpenGL1Renderer::copyTexImage2D(GLenum target, GLint level,
    GLenum internalformat,
    GLint x, GLint y,
    GLsizei width, GLsizei height,
    GLint border){ /*glCopyTexImage2D(target, level, internalformat, x, y, width, height, border);*/ }




void Renderer::drawMesh(GLMesh* mesh, Vec3f position, Quat4f rotation, Vec3f scale){
    if (active_camera == nullptr){
        printf("Warning: No active camera - skipping rendering.\n");
        return;
    }

    // scale
    Mat4f modelScaleMatrix = Mat4fIdentity;
    modelScaleMatrix.mul({scale.x, scale.y, scale.z, 1.0});
    
    // rotation
    Mat3f modelRotMatrix3;
    rotation.toMatrix(modelRotMatrix3);
    Mat4f modelRotMatrix4 = Mat4fIdentity;
    modelRotMatrix4.setRot(modelRotMatrix3);

    // translation
    Mat4f modelTranslationMatrix = Mat4fIdentity;
    modelTranslationMatrix.setPos(position);

    // model matrix
    Mat4f modelMatrix = Mat4fIdentity;
    modelMatrix.mmulL(modelScaleMatrix);
    modelMatrix.mmulL(modelRotMatrix4);
    modelMatrix.mmulL(modelTranslationMatrix);

    // MVP matrix
    Mat4f mvpMatrix = modelMatrix;
    mvpMatrix.mmulL(active_camera->viewMatrix());
    mvpMatrix.mmulL(active_camera->projectionMatrix());

    drawMeshMVP(mesh, mvpMatrix);
}
