#ifndef _INPUT_MANAGER_H_
#define _INPUT_MANAGER_H_
#include <SDL2/SDL_events.h>
#include <functional>
#include <vector>

class InputManager{
    using EventHandler = std::function<bool(const SDL_Event)>;
    std::vector<EventHandler> handlers;
    Uint8 keystates[SDL_NUM_SCANCODES];

public:
    InputManager(){memset(keystates, false, SDL_NUM_SCANCODES);}

    bool eventHandling( const SDL_Event& event ){
        for (auto handler : handlers){
            if (handler(event)) return true;
        }

        if (event.type == SDL_KEYDOWN || event.type == SDL_KEYUP){
            keystates[event.key.keysym.scancode] = event.key.state == SDL_PRESSED;
        }

        return false;
    }

    void addEventHandler( EventHandler handler ){
        handlers.push_back(handler);
    }

    bool isKeyPressed( SDL_Scancode key ){
        return keystates[key];
    }
};

#endif