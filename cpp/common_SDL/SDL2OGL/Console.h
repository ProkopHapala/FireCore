
#ifndef  Console_h
#define  Console_h

#include <SDL2/SDL.h>


#include "Draw.h"
#include "Draw2D.h"
//#include "Draw3D.h"

//#include "Table.h"
#include <string>
#include <functional>
#include "CircularArray.h"
#include "SortedStrings.h"

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
    int  nminMatch=2;

    int imatch=-1,imatchn=0;
    //bool historyOpen=false;
    int iHistory=-1;

    // lambda function for callback
    std::function<bool(const char*)> callback;
    //std::vector<std::string> history;   // previous commands (history)
    //std::string(string[]) quick_tab;  // quick tab completion
    CircularArray<char*> history{ 16, true, true }; 

    SortedStrings quick_tab; // quick tab completion

    // ===== Functions

    void close_history(){
        //historyOpen=false;
        iHistory=-1;
    }
    void accept_history(){
        if(iHistory>=0){
            const char* s = history.get(iHistory);
            nend = strlen(s);
            if(ncur>nend){ ncur=nend; }
            strcpy( line, s );
            close_history();
        }
    }

    void init( int lineLength_=1024, SDL_Window* window_=0 ){
        if(lineLength_>0){ 
            lineLength = lineLength_;
            line = new char[lineLength]; 
        }
        if(window_){ window=window_; }
        //for(int i=0;i<lineLength;i++){ line[i]='\0'; }
        set(lineLength,'\0');
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
    void set(int n,char c){ for(int i=0; i<n; i++){ line[i]=c; } }
    void fwd_copy( int n, const char* src, char* dst ){ for(int i=0;   i<n; i++ ){ dst[i] = src[i]; } }
    void bwd_copy( int n, const char* src, char* dst ){ for(int i=n-1; i>=0; i--){ dst[i] = src[i]; } }

    // console.keyDown( event.key.keysym.sym )
    bool keyDown( const SDL_Keycode key ){ 
        if(iHistory<0){ iHistory=-1; }
        SDL_Keymod modState = SDL_GetModState();  // Get the current state of modifier keys
        bool bShift = modState & KMOD_SHIFT;
        //printf( "Console::keyDown() bShift=%i \n", bShift );
        //printf( "Console::keyDown(key=%i|%c) ncur=%i\n", key, key, ncur );
        switch( key ){
            case SDLK_BACKQUOTE: return false; break; 
            case SDLK_KP_ENTER: //  [[fallthrough]]
            case SDLK_RETURN: {
                //printf( "Console::run(%s) \n", line ); 
                char* s = line;
                if(iHistory>=0){ s = history.get(iHistory); } 
                if( callback( s ) ){
                    //printf("Console::SDLK_RETURN() callback() OK\n");
                    if( (nend>0) && (iHistory<0) ){ // save to history
                        char* ss = new char[nend+1];
                        fwd_copy(nend+1, s, ss);
                        //printf( "Console::SDLK_RETURN() history.push(%s) line(%s) \n", ss, line );
                        history.push( ss );
                    }else{
                        //printf( "Console::SDLK_RETURN() not-push() iHistory=%i nend=%i \n", iHistory, nend );
                    }
                    close_history();
                    ncur=0;
                    nend=0;
                    set(lineLength,'\0');
                }
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
            case SDLK_LEFT: { accept_history();    if(ncur>0   ){ ncur--; }} break;
            case SDLK_RIGHT:{ accept_history();    if(ncur<nend){ ncur++; }} break;

            case SDLK_ESCAPE:    close_history(); break;
            case SDLK_TAB: {
                //printf( "Console::keyDown(SDLK_TAB) iHistory=%i \n", iHistory );
                if(iHistory>=0){
                    accept_history();
                }else if(imatch>=0){
                    const std::string& s = quick_tab.table[imatch]; 
                    nend = s.size();
                    strcpy( line, quick_tab.table[imatch].c_str() );
                    printf( "Console::keyDown(SDLK_TAB) quick_tab[%i] nend=%i `%s` line=`%s`\n", imatch, nend, quick_tab.table[imatch].c_str(), line );
                    //close_history();
                }
                //accept_history();
            } break;

            case SDLK_UP:{ 
                //printf("Console::keyDown() history.size() %i \n", history.size() );
                int ih = iHistory+1;
                if(ih<history.size()){
                    iHistory = ih;
                }
            }break;
            case SDLK_DOWN: {
                iHistory-=1; if(iHistory<0){ iHistory=-1; }
            } break;
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
                        c = kp_scancode[k];
                        //printf(  "key(%c) %i   %i \n", c, k, key );
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
        if( (iHistory<0)&&( nend>nminMatch ) ){     
            imatch = quick_tab.findMatch( line );
            //printf( "Console::keyDown() quick_tab.findMatch() imatch(%i) \n", imatch );
            if(imatch>=0){ 
                //printf( "Console::keyDown() quick_tab.findMatch() imatch(%i)>0 \n", imatch );
                imatchn = quick_tab.findMatchEnd( line, imatch );
                //printf( "Console::keyDown() imatch %i \n", imatch );
                //printf( "Console::keyDown() quick_tab[%i] `%s` imatchn=%i \n", imatch, quick_tab.table[imatch].c_str(), imatchn );
            }else{
                imatchn=0;
            }
        }
        return true;
    }

    void draw(){
        int w,h;
        SDL_GetWindowSize(window, &w, &h);
        //printf( "Console::draw() w %i h %i \n", w, h );
        opengl1renderer.color3f(0.0f,0.0f,0.0f); Draw2D::drawRectangle( 0,h-fontSizeDef*2, w,h );
        //opengl1renderer.color3f(0.0f,1.0f,0.0f); Draw2D::drawRectangle( w/2,h/2, w,h );   // just debugging
        opengl1renderer.color3f(0.0f,1.0f,0.0f);
        //Draw2D::drawRectangle( w/2,h*2, w,h, false );
        float xcur = fontSizeDef*ncur;
        Draw2D::drawRectangle( xcur,h, xcur+fontSizeDef,h-fontSizeDef*2, false );
        opengl1renderer.color3f(1.0f,1.0f,1.0f);
        //Draw2D::drawString( line, {0,0}, 0.1, 0 );
        //Draw2D::drawText( caption.c_str(), caption.length(), {xmin, ymax-fontSizeDef*2}, 0.0, GUI_fontTex, fontSizeDef );
        Draw2D::drawText( line           , nend            , {0   , h   -fontSizeDef*2}, 0.0, fontTex, fontSizeDef );
        if(iHistory>=0){
            const char* s = history.get(iHistory);
            int n = strlen(s);
            //printf( "Console::draw() history[%i] history.size(%i) n=%i @%li `%s`\n", iHistory, history.size(), n, (long)s, s );
            //opengl1renderer.color3f(0.3f,0.3f,0.3f); Draw2D::drawRectangle( 0,h+fontSizeDef*2, n*fontSizeDef,h-fontSizeDef*4, true );
            //opengl1renderer.color3f(1.0f,1.0f,1.0f); Draw2D::drawText( s, n, {0, h-fontSizeDef*4}, 0.0, fontTex, fontSizeDef );
            opengl1renderer.color3f(0.0f,0.0f,0.5f); Draw2D::drawRectangle( 0,h, n*fontSizeDef,h-fontSizeDef*2, true );
            opengl1renderer.color3f(1.0f,1.0f,1.0f); Draw2D::drawText( s, n, {0, h-fontSizeDef*2}, 0.0, fontTex, fontSizeDef );
        }else
        if( (imatch>=-1)&&(imatchn>0)){
            //opengl1renderer.color3f(0.0f,0.0f,1.0f); Draw2D::drawRectangle( 0,h-fontSizeDef*2, n*fontSizeDef,h-(1+imatchn)*fontSizeDef*2, true );
            //opengl1renderer.color3f(1.0f,1.0f,1.0f);
            for(int i=0; i<imatchn; i++){
                //printf( "Console::draw(imatch..imatchn) imatch=%i, i=%i, quick_tab.table.size(%i) \n", imatch, i, quick_tab.table.size() );
                const char* s = quick_tab.table[imatch+i].c_str();
                int n = strlen(s);
                //printf( "Console::draw() quick_tab[%i] quick_tab.size(%i) n=%i @%li `%s`\n", imatch, quick_tab.table.size(), n, (long)s, s );
                opengl1renderer.color3f(0.5f,0.0f,0.5f); Draw2D::drawRectangle( 0,h-(1+i)*fontSizeDef*2, n*fontSizeDef,h-(2+i)*fontSizeDef*2, true );
                opengl1renderer.color3f(1.0f,1.0f,1.0f); Draw2D::drawText( s, n, {0, h-(2+i)*fontSizeDef*2}, 0.0, fontTex, fontSizeDef );
            }
        }
    }

};

#endif
