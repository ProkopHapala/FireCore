
#ifndef  SDL_utils_h
#define  SDL_utils_h

#include <stdlib.h>
#include <stdio.h>




float* loadDataImageFloat( char * fname, float hsc, int& nx, int& ny, int& BytesPerPixel ){
    //SDL_Surface * surf = IMG_Load( fname );
    SDL_Surface * surf = SDL_LoadBMP( fname );
    if ( surf ){
        nx = surf->w;
        ny = surf->h;
        BytesPerPixel = surf->format->BytesPerPixel;
        int ntot      = nx * ny * BytesPerPixel ;
        float * ret   = new float[ntot];
        float renorm  = hsc / 256.0;
        for(int i=0; i<ntot; i++){
            ret[i] = renorm * ((uint8_t*)surf->pixels)[i];
        }
        SDL_FreeSurface( surf );
        return ret;
    }else{ printf( "cannot load %s\n", fname  ); }
    return 0;
};









GLuint makeTexture( char * fname ){

    //SDL_Surface * surf = IMG_Load( fname );
    SDL_Surface * surf = SDL_LoadBMP( fname );
    if ( surf ){
        GLuint itex=0;
        opengl1renderer.genTextures  ( 1, &itex );
        opengl1renderer.bindTexture  ( GL_TEXTURE_2D, itex );
        //if      (surf->format->BytesPerPixel == 1) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, 1,  surf->w,  surf->h, 0, GL_LUMINANCE, GL_UNSIGNED_BYTE, surf->pixels ); }
        //if      (surf->format->BytesPerPixel == 1) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, 1,   surf->w,  surf->h, 0, GL_INTENSITY, GL_UNSIGNED_BYTE, surf->pixels ); }
        //if      (surf->format->BytesPerPixel == 1) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, 1,       surf->w,  surf->h, 0, GL_ALPHA,     GL_UNSIGNED_BYTE, surf->pixels ); }
        if      (surf->format->BytesPerPixel == 1) {
            opengl1renderer.pixelStorei(GL_UNPACK_ALIGNMENT, 1);
            opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
            opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
            opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, 1,       surf->w,  surf->h, 0, GL_RED,       GL_UNSIGNED_BYTE, surf->pixels );
        }
        else if (surf->format->BytesPerPixel == 3) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, GL_RGB,  surf->w,  surf->h, 0, GL_BGR,       GL_UNSIGNED_BYTE, surf->pixels ); }
        else if (surf->format->BytesPerPixel == 4) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, GL_RGBA, surf->w,  surf->h, 0, 0x8000,       GL_UNSIGNED_BYTE, surf->pixels ); }
        else return 0;
        //printf( "surface->format->Rmask : %i itex %i \n", surf->format->BytesPerPixel, itex  ) ;// surface->format->Rmask/ == 0x000000ff;

        //opengl1renderer.texImage2D   ( GL_TEXTURE_2D, 0, 3, surf->w,  surf->h, 0, GL_BGR, GL_UNSIGNED_BYTE, surf->pixels );
        //opengl1renderer.texImage2D   ( GL_TEXTURE_2D, 0, 3, surf->w,  surf->h, 0, GL_BGRA, GL_UNSIGNED_BYTE, surf->pixels );
        //opengl1renderer.texImage2D(GL_TEXTURE_2D, 0, GL_RGBA, surf->w,  surf->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, surf->pixels );
        //opengl1renderer.texImage2D(GL_TEXTURE_2D, 0, GL_BGRA, surf->w,  surf->h, 0, GL_BGRA, GL_UNSIGNED_BYTE, surf->pixels );
        //opengl1renderer.texImage2D(GL_TEXTURE_2D, 0,  GL_RGBA, surf->w,  surf->h, 0, GL_RGBA, GL_UNSIGNED_BYTE, surf->pixels );

        opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
        opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
        //opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
        //opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        SDL_FreeSurface( surf );
        return itex;
    }else{
        printf( "cannot load %s\n", fname  );
    }
    return 0;
};


GLuint makeTextureHard( char * fname ){

    //SDL_Surface * surf = IMG_Load( fname );
    SDL_Surface * surf = SDL_LoadBMP( fname );
    if ( surf ){
        GLuint itex=0;
        opengl1renderer.genTextures  ( 1, &itex );
        opengl1renderer.bindTexture  ( GL_TEXTURE_2D, itex );
        if      (surf->format->BytesPerPixel == 1) {
            opengl1renderer.pixelStorei(GL_UNPACK_ALIGNMENT, 1);
            opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
            opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
            opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, 1,       surf->w,  surf->h, 0, GL_RED,       GL_UNSIGNED_BYTE, surf->pixels );
        }
        else if (surf->format->BytesPerPixel == 3) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, GL_RGB,  surf->w,  surf->h, 0, GL_BGR,       GL_UNSIGNED_BYTE, surf->pixels ); }
        else if (surf->format->BytesPerPixel == 4) { opengl1renderer.texImage2D( GL_TEXTURE_2D, 0, GL_RGBA, surf->w,  surf->h, 0, 0x8000,       GL_UNSIGNED_BYTE, surf->pixels ); }
        else return 0;
        //printf( "surface->format->Rmask : %i itex %i \n", surf->format->BytesPerPixel, itex  ) ;// surface->format->Rmask/ == 0x000000ff;

        opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
        opengl1renderer.texParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
        SDL_FreeSurface( surf );
        return itex;
    }else{
        printf( "cannot load %s\n", fname  );
    }
    return 0;
};

#endif
