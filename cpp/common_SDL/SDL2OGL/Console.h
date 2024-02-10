
#ifndef  Console_h
#define  Console_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
//#include "Draw3D.h"

//#include "Table.h"
#include <string>
//#include <functional>

    // SDL_SCANCODE_1 = 30,
    // SDL_SCANCODE_2 = 31,
    // SDL_SCANCODE_3 = 32,
    // SDL_SCANCODE_4 = 33,
    // SDL_SCANCODE_5 = 34,
    // SDL_SCANCODE_6 = 35,
    // SDL_SCANCODE_7 = 36,
    // SDL_SCANCODE_8 = 37,
    // SDL_SCANCODE_9 = 38,
    // SDL_SCANCODE_0 = 39,

    // SDL_SCANCODE_KP_DIVIDE   = 84,    1
    // SDL_SCANCODE_KP_MULTIPLY = 85,    2
    // SDL_SCANCODE_KP_MINUS    = 86,    3
    // SDL_SCANCODE_KP_PLUS     = 87,    4
    // SDL_SCANCODE_KP_ENTER    = 88,    5
    // SDL_SCANCODE_KP_1 = 89,           6
    // SDL_SCANCODE_KP_2 = 90,           7
    // SDL_SCANCODE_KP_3 = 91,           8
    // SDL_SCANCODE_KP_4 = 92,           9
    // SDL_SCANCODE_KP_5 = 93,          10
    // SDL_SCANCODE_KP_6 = 94,          11
    // SDL_SCANCODE_KP_7 = 95,          12
    // SDL_SCANCODE_KP_8 = 96,          13
    // SDL_SCANCODE_KP_9 = 97,          14
    // SDL_SCANCODE_KP_0 = 98,          15
    // SDL_SCANCODE_KP_PERIOD = 99,     16

constexpr char printableChars     [] = " !\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
constexpr char printableCharsShift[] = " !\"#$%&\"()*+<_>?)!@#$%^&*(::<+>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ{|}^_`ABCDEFGHIJKLMNOPQRSTUVWXYZ{|}~";
constexpr char kp_scancode        [] = "/*-+ 1234567890.";

class Console{ public:
    bool bShift = false;
    int fontTex = 0;
    SDL_Window* window=0;
    int lineLength=256;
    char * line=0;
    int  ncur=0;
    int  nend=0;
    // lambda function for callback
    std::function<void(const char*)> callback;
    //std::vector(string[]) prev_cmd;   // previous commands (history)
    //std::string(string[]) quick_tab;  // quick tab completion

    // ===== Functions

    void init( int lineLength_=1024, SDL_Window* window_=0 ){
        if(lineLength_>0){ 
            lineLength = lineLength_;
            line = new char[lineLength]; 
        }
        if(window_){ window=window_; }
        for(int i=0;i<lineLength;i++){ line[i]='\0'; }
        //printf( "Console::init() %i %i %i \n", lineLength, window, window_ );
    }

    // void fwd_copy( int n, const char* src, char* dst ){
    //     //printf( "Console::fwd_copy() n %i\n", n );
    //     for(int i=0;   i<n; i++ ){ 
    //         //printf("%c", src[i] );
    //         dst[i] = src[i]; 
    //     }
    //     //printf( "Console::fwd_copy() end\n" );
    // }
    void fwd_copy( int n, const char* src, char* dst ){ for(int i=0;   i<n; i++ ){ dst[i] = src[i]; } }
    void bwd_copy( int n, const char* src, char* dst ){ for(int i=n-1; i>=0; i--){ dst[i] = src[i]; } }

    // console.keyDown( event.key.keysym.sym )
    bool keyDown( const SDL_Keycode key ){ 
        SDL_Keymod modState = SDL_GetModState();  // Get the current state of modifier keys
        bool bShift = modState & KMOD_SHIFT;
        //printf( "Console::keyDown() bShift=%i \n", bShift );
        //printf( "Console::keyDown(key=%i) ncur=%i\n", key, ncur );
        switch( key ){
            case SDLK_BACKQUOTE: return false; break; 
            case SDLK_KP_ENTER: //  [[fallthrough]]
            case SDLK_RETURN: {
                //printf( "Console::run(%s) \n", line ); 
                callback( line );
            } break;
            case SDLK_BACKSPACE: 
            if(ncur>0   ){ 
                //if(ncur<=nend)
                fwd_copy(nend-ncur, line+ncur,   line+ncur-1  ); 
                line[nend-1]='\0'; 
                nend--; 
                ncur--; 
            } break;
            case SDLK_DELETE:    
            if(nend>ncur){ 
                if(ncur< nend)
                fwd_copy(nend-ncur-1, line+ncur+1, line+ncur    ); 
                line[nend-1]='\0';
                nend--;         
            } break;
            case SDLK_LEFT:      if(ncur>0   ){ ncur--; } break;
            case SDLK_RIGHT:     if(ncur<nend){ ncur++; } break;

            //case SDLK_LSHIFT: //  [[fallthrough]]
            //case SDLK_RSHIFT: {bShift=!bShift;}  break;
            // // other non-printable keys
            // case SDLK_ESCAPE:
            // case SDLK_TAB:
            // case SDLK_CAPSLOCK:
            // case SDLK_F1 ... SDLK_F12:  // Range for function keys
            // case SDLK_PRINTSCREEN:
            // case SDLK_SCROLLLOCK:
            // case SDLK_PAUSE:
            // case SDLK_INSERT:
            // //case SDLK_DELETE:
            // case SDLK_HOME:
            // case SDLK_END:
            // case SDLK_PAGEUP:
            // case SDLK_PAGEDOWN:
            // //case SDLK_RIGHT:
            // //case SDLK_LEFT:
            // case SDLK_DOWN:
            // case SDLK_UP:
            // case SDLK_NUMLOCKCLEAR:
            // //case SDLK_LSHIFT:
            // //case SDLK_RSHIFT:
            // case SDLK_RALT:
            // case SDLK_LALT:
            // case SDLK_LCTRL:
            // case SDLK_RCTRL:
            // case SDLK_LGUI:
            // case SDLK_RGUI:
            // case SDLK_MENU:
            // case SDLK_SYSREQ:

            default:             
            if(nend<lineLength){ 
                char c=0;
                if( (key>=31)&&(key<128) ){
                    //printf(  "key(%c|%c|%c) %i \n", key, printableChars[key-32], printableCharsShift[key-32], key );
                    //c = key-32;
                    if(bShift){ 
                        c = printableCharsShift[key-32]; 
                    }else{
                        c = key;
                    }
                }else{
                    int k = (key&0xff)-84;
                    if( (k>=0) && (k<16) ){
                        //printf(  "key(%c) %i   %i \n", kp_scancode[k], k, key );
                        c = kp_scancode[c];
                    }
                }
                if(c!=0){
                    //printf( "key=%i c=%c \n", key, c );
                    if(ncur<nend){ bwd_copy(nend-ncur, line+ncur, line+ncur+1); }
                    line[ncur]=c; 
                    ncur++; 
                    nend++; 
                }
            } break;
        }
        
        // if(ncur<nend){
        //     line[ncur]=key;
        //     ncur++;
        // }
        // if(key==SDLK_BACKQUOTE){ return false; }
        // else if(key==SDLK_RETURN){
        //     printf( "Console::run(%s)\n", key, line );
        //     callback( line );
        //     //if( line[0] != 0 ){ printf( "Console::keyDown(%i) line %s\n", key, line ); }
        //     return true;
        // }else if(key==SDLK_BACKSPACE){
        //     int n = strlen(line);
        //     if(n>0){ line[n-1] = 0; }
        //     return true;
        // }
        return true;
    }

    void draw(){
        int w,h;
        SDL_GetWindowSize(window, &w, &h);
        glColor3f(0.0f,0.0f,0.0f);
        Draw2D::drawRectangle( 0,h-fontSizeDef*2, w,h );
        glColor3f(0.0f,1.0f,0.0f);
        //Draw2D::drawRectangle( w/2,h*2, w,h, false );
        float xcur = fontSizeDef*ncur;
        Draw2D::drawRectangle( xcur,h, xcur+fontSizeDef,h-fontSizeDef*2, false );
        glColor3f(1.0f,1.0f,1.0f);
        //Draw2D::drawString( line, {0,0}, 0.1, 0 );
        //Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2}, 0.0, GUI_fontTex, fontSizeDef );
        Draw2D::drawText( line           , nend            , {0   , h   -fontSizeDef*2}, 0.0, fontTex, fontSizeDef );

    }

};

#endif
