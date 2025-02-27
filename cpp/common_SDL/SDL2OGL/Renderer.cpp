#include <SDL2/SDL_opengl.h>

#define __SKIP_DEFS__
#include "Renderer.h"


// the global variable for the renderer
OpenGL1Renderer opengl1renderer = OpenGL1Renderer();


// OpenGL 1.x functions
void OpenGL1Renderer::begin(GLenum mode) { glBegin(mode); }
void OpenGL1Renderer::end() { glEnd(); }

GLuint OpenGL1Renderer::genLists(GLsizei range) { return glGenLists(range); }
void OpenGL1Renderer::deleteLists(GLuint list, GLsizei range) { glDeleteLists(list, range); }

void OpenGL1Renderer::endList() { glEndList(); }
void OpenGL1Renderer::vertex(const Vec3f& v ){ glVertex3f(v.x,v.y,v.z); }
void OpenGL1Renderer::vertex(const Vec3d& v ){ glVertex3f(v.x,v.y,v.z); }
void OpenGL1Renderer::color (const Vec3f& v ){ glColor3f (v.x,v.y,v.z); }
void OpenGL1Renderer::color (const Vec3d& v ){ glColor3f (v.x,v.y,v.z); }
void OpenGL1Renderer::normal(const Vec3f& v ){ glNormal3f(v.x,v.y,v.z); }
void OpenGL1Renderer::normal(const Vec3d& v ){ glNormal3f(v.x,v.y,v.z); }

void OpenGL1Renderer::newList(GLuint list, GLenum mode){ glNewList(list, mode); }
void OpenGL1Renderer::clearColor(float r, float g, float b, float a){ glClearColor(r,g,b,a); }
void OpenGL1Renderer::clear(GLbitfield mask){ glClear(mask); }
void OpenGL1Renderer::enable(GLenum cap){ glEnable(cap); }
void OpenGL1Renderer::disable(GLenum cap){ glDisable(cap); }
void OpenGL1Renderer::callList(GLuint list){ glCallList(list); }
void OpenGL1Renderer::pushMatrix(){ glPushMatrix(); }
void OpenGL1Renderer::popMatrix(){ glPopMatrix(); }
void OpenGL1Renderer::depthFunc(GLenum func){ glDepthFunc(func); }
void OpenGL1Renderer::blendFunc(GLenum sfactor, GLenum dfactor){ glBlendFunc(sfactor,dfactor); }
void OpenGL1Renderer::translatef(float x, float y, float z){ glTranslatef(x,y,z); }
void OpenGL1Renderer::flush(){ glFlush(); }
void OpenGL1Renderer::finish(){ glFinish(); }
void OpenGL1Renderer::readPixels(GLint x, GLint y, GLsizei width, GLsizei height, GLenum format, GLenum type, GLvoid *pixels){ glReadPixels(x,y,width,height,format,type,pixels); }
void OpenGL1Renderer::shadeModel(GLenum mode){ glShadeModel(mode); }
void OpenGL1Renderer::viewport(GLint x, GLint y, GLsizei width, GLsizei height){ glViewport(x,y,width,height); }
void OpenGL1Renderer::bindTexture(GLenum target, GLuint texture){ glBindTexture(target,texture); }
void OpenGL1Renderer::texParameteri(GLenum target, GLenum pname, GLint param){ glTexParameteri(target,pname,param); }
void OpenGL1Renderer::matrixMode(GLenum mode){ glMatrixMode(mode); }
void OpenGL1Renderer::multMatrixf(const GLfloat *m){ glMultMatrixf(m); }
void OpenGL1Renderer::lineWidth(GLfloat width){ glLineWidth(width); }

void OpenGL1Renderer::genTextures(GLsizei n, GLuint *textures){glGenTextures(n, textures);}
void OpenGL1Renderer::pixelStorei(GLenum pname, GLint param){glPixelStorei(pname, param);}
void OpenGL1Renderer::loadIdentity(){glLoadIdentity();}
void OpenGL1Renderer::frustum(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar){glFrustum(left, right, bottom, top, zNear, zFar);}
void OpenGL1Renderer::ortho(GLdouble left, GLdouble right, GLdouble bottom, GLdouble top, GLdouble zNear, GLdouble zFar){glOrtho(left, right, bottom, top, zNear, zFar);}
void OpenGL1Renderer::getFloatv(GLuint pname, GLfloat *params){glGetFloatv(pname, params);}
void OpenGL1Renderer::rotatef(GLfloat angle, GLfloat x, GLfloat y, GLfloat z){glRotatef(angle, x, y, z);}
void OpenGL1Renderer::scalef(GLfloat x, GLfloat y, GLfloat z){glScalef(x, y, z);}
void OpenGL1Renderer::pointSize(GLfloat size){glPointSize(size);}
void OpenGL1Renderer::scissor(GLint x, GLint y, GLsizei width, GLsizei height){glScissor(x, y, width, height);}
void OpenGL1Renderer::frontFace(GLenum mode){glFrontFace(mode);}
void OpenGL1Renderer::materialfv(GLenum face, GLenum pname, const GLfloat *params){glMaterialfv(face, pname, params);}
void OpenGL1Renderer::lightModeli(GLenum pname, GLint param){glLightModeli(pname, param);}
void OpenGL1Renderer::polygonMode(GLenum face, GLenum mode){glPolygonMode(face, mode);}

void OpenGL1Renderer::color3f(float r, float g, float b){glColor3f(r,g,b);}
void OpenGL1Renderer::vertex3f(float x, float y, float z){glVertex3f(x,y,z);}
void OpenGL1Renderer::color4f(float r, float g, float b, float a){glColor4f(r,g,b,a);}

void OpenGL1Renderer::lightfv(GLenum light, GLenum pname, const GLfloat *params){glLightfv(light, pname, params);}
void OpenGL1Renderer::texImage2D(GLenum target, GLint level, GLint internalformat, GLsizei width, GLsizei height, GLint border, GLenum format, GLenum type, const GLvoid *pixels){glTexImage2D(target, level, internalformat, width, height, border, format, type, pixels);}
void OpenGL1Renderer::texCoord2f(GLfloat s, GLfloat t){glTexCoord2f(s,t);}
void OpenGL1Renderer::hint(GLenum target, GLenum mode){glHint(target, mode);}
void OpenGL1Renderer::loadMatrixf(const GLfloat *m){glLoadMatrixf(m);}
void OpenGL1Renderer::lightModelfv(GLenum pname, const GLfloat *params){glLightModelfv(pname, params);}

void OpenGL1Renderer::vertex3d(double x, double y, double z){glVertex3d(x, y, z);}
void OpenGL1Renderer::normal3f(float x, float y, float z){glNormal3f(x, y, z);}
void OpenGL1Renderer::normal3d(double x, double y, double z){glNormal3d(x, y, z);}
void OpenGL1Renderer::vertex2d(double x, double y){glVertex2d(x, y);}
void OpenGL1Renderer::color3d(double r, double g, double b){glColor3d(r, g, b);}

const GLubyte* OpenGL1Renderer::getString(GLenum name){return glGetString(name);}

void OpenGL1Renderer::copyTexImage2D(GLenum target, GLint level,
    GLenum internalformat,
    GLint x, GLint y,
    GLsizei width, GLsizei height,
    GLint border){glCopyTexImage2D(target, level, internalformat, x, y, width, height, border);}
