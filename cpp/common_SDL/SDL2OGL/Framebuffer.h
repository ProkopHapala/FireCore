
#include "GLES2.h"
#include "GLTexture.h"
#include <cstdlib>

class GLFramebuffer{
public:
    GLTexture colorBuffer = GLTexture(GLES2::screen_size, GL_RGBA);

private:
    GLuint handle = 0;
    GLuint depthBufferHandle = 0;
    bool active = false;
    bool paused = false;

    void ensure_handle(){
        if (handle) return;

        // adapted from https://www.opengl-tutorial.org/intermediate-tutorials/tutorial-14-render-to-texture
        glGenFramebuffers(1, &handle);
        glBindFramebuffer(GL_FRAMEBUFFER, handle);

        colorBuffer.setMagFilter(GL_NEAREST);
        colorBuffer.setMinFilter(GL_NEAREST);


        GL_CHECK_ERROR();

        glGenRenderbuffers(1, &depthBufferHandle);
        glBindRenderbuffer(GL_RENDERBUFFER, depthBufferHandle);
        // TODO: GL_DEPTH_COMPONENT16 is more portable, but requires a smaller zmin/zmax in Camera.h
        glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24_OES, GLES2::screen_size.x, GLES2::screen_size.y);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufferHandle);

        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, colorBuffer.getHandle(), 0); GL_CHECK_ERROR();

        //Check framebuffer status
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER); GL_CHECK_ERROR();
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            char* status_str;
            switch(status) {
                case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:         status_str = "GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT"; break;
                case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT: status_str = "GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT"; break;
                case GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS: status_str = "GL_FRAMEBUFFER_INCOMPLETE_DIMENSIONS"; break;
                case GL_FRAMEBUFFER_UNSUPPORTED: status_str = "GL_FRAMEBUFFER_UNSUPPORTED"; break;
                default: status_str = "Unknown error"; break;
            }

            printf("Framebuffer is not complete: %s\n", status_str);
            exit(0);
        }
    }

public:
    GLFramebuffer(){};
    ~GLFramebuffer(){ if (handle) glDeleteFramebuffers(1, &handle); }

    void begin(){
        if (active) {
            printf("ERROR: framebuffer is already active\n");
            exit(1);
        }

        ensure_handle();
        GLES2::pushFramebuffer(handle);
        glViewport(0, 0, GLES2::screen_size.x, GLES2::screen_size.y);
        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        active = true;
    }

    void end(){
        if (paused){
            printf("ERROR: framebuffer is paused\n");
            exit(1);
        }
        if (!active){
            printf("ERROR: framebuffer is not active\n");
            exit(1);
        }
        GLES2::popFramebuffer(handle);
        active = false;
    }

    void pause(){
        if (!active){
            printf("ERROR: framebuffer is not active\n");
            exit(1);
        }
        if (paused){
            printf("ERROR: framebuffer is already paused\n");
            exit(1);
        }
        GLES2::popFramebuffer(handle);
        paused = true;
    }

    void unpause(){
        if (!active){
            printf("ERROR: framebuffer is not active\n");
            exit(1);
        }
        if (!paused){
            printf("ERROR: framebuffer is not paused\n");
            exit(1);
        }
        GLES2::pushFramebuffer(handle);
        paused = false;
    }
};
