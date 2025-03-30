#ifndef _GL_TEXTURE_H_
#define _GL_TEXTURE_H_

#include "GLES2.h"
#include "Vec3.h"
#include <GLES2/gl2.h>
#include <SDL2/SDL_surface.h>
#include <cstring>


class GLTexture{
    private:
        char* path = nullptr;
        Vec2i size = {0, 0};
        GLenum format = 0;

        GLuint handle = 0;

        void ensure_handle(){
            if (!handle){
                glGenTextures(1, &handle);
                if (path){
                    loadBMP(path);
                    free(path);
                    path = nullptr;                    
                }else{
                    glBindTexture(GL_TEXTURE_2D, handle);
                    glTexImage2D(GL_TEXTURE_2D, 0, format, size.x, size.y, 0, format, GL_UNSIGNED_BYTE, 0);
                }
            }
        }
    
    public:

        GLTexture() {};
        GLTexture(const char* path){ 
            this->path = (char*)malloc(strlen(path)+1);
            strcpy(this->path, path); // lazy initialisation
        }
        GLTexture(Vec2i size, GLenum format=GL_RGBA) : size(size), format(format){} // lazy initialisation
        ~GLTexture(){
            if (handle) glDeleteTextures(1, &handle);
        }

        inline GLuint getHandle() { ensure_handle(); return handle;};

        void bind(GLenum target=GL_TEXTURE_2D){
            ensure_handle();
            glBindTexture(target, handle);
        }

        void setMinFilter(GLenum filter)    { bind(); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter ); } // default is GL_NEAREST_MIPMAP_LINEAR
        void setMagFilter(GLenum filter)    { bind(); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter ); } // default is GL_LINER
        void setWrapHorizontal(GLenum wrap) { bind(); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S    , wrap   ); } // default is GL_REPEAT
        void setWrapVertical  (GLenum wrap) { bind(); glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T    , wrap   ); } // default is GL_REPEAT

        void loadImage(int width, int height, GLenum format, GLenum type, const void* data){
            bind();
            glTexImage2D(GL_TEXTURE_2D, 0, format, width, height, 0, format, type, data);
            size = {width, height};
            this->format = format;
        }

        void loadBMP(const char* filename){
            SDL_Surface* surf = SDL_LoadBMP( filename );
            if (!surf){
                printf( "ERROR: could not load image '%s'\n", filename );
                return;
            }
            bind();
            if (surf->format->BytesPerPixel == 1) {
                glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
                loadImage(surf->w, surf->h, GL_ALPHA, GL_UNSIGNED_BYTE, surf->pixels);
            }
            else if (surf->format->BytesPerPixel == 3) { loadImage(surf->w, surf->h, GL_RGB,  GL_UNSIGNED_BYTE, surf->pixels); }
            else if (surf->format->BytesPerPixel == 4) { loadImage(surf->w, surf->h, GL_RGBA, GL_UNSIGNED_BYTE, surf->pixels); }
            else {
                printf("ERROR: could not load image '%s', unknown bytes per pixel (%i)\n", filename, surf->format->BytesPerPixel);
                SDL_FreeSurface(surf);
                return;
            }
            SDL_FreeSurface( surf );
        }
};

#endif // _GL_TEXTURE_H_