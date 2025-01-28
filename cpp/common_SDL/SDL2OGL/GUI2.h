#include <SDL2/SDL_video.h>
#include <cmath>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include "Draw.h"
#include "Vec2.h"


template<class T>
class GUI2Rect2T{ // switch to Rect2d already implemented in geom2D.h ?
    using REC = GUI2Rect2T<T>;
    using VEC = Vec2T<T>;

    public:
    T xmin, ymin, xmax, ymax;

    GUI2Rect2T() : xmin(0), ymin(0), xmax(0), ymax(0) {}
    GUI2Rect2T(const T xmin_, const T ymin_, const T xmax_, const T ymax_) : xmin(xmin_), ymin(ymin_), xmax(xmax_), ymax(ymax_) {}
    GUI2Rect2T(const VEC vmin, const VEC vmax) : xmin(vmin.x), ymin(vmin.y), xmax(vmax.x), ymax(vmax.y) {}

    VEC size(){ return (VEC){xmax-xmin, ymax-ymin}; }
    VEC min() { return (VEC){xmin, ymin}; }
    VEC max() { return (VEC){xmax, ymax}; }

    bool contains( VEC pos ){
        return (pos.x >= xmin && pos.x <= xmax && pos.y >= ymin && pos.y <= ymax);
    }

    bool operator==( const REC& r ){ return (xmin==r.xmin && ymin==r.ymin && xmax==r.xmax && ymax==r.ymax); }
    bool operator!=( const REC& r ){ return (xmin!=r.xmin || ymin!=r.ymin || xmax!=r.xmax || ymax!=r.ymax); }
};

using GUI2Rect2i = GUI2Rect2T<int>;
using GUI2Rect2f = GUI2Rect2T<float>;
#define FULL_RECT GUI2Rect2f(0,0,1,1)


extern int GUI2_fontTex;

class GUI2Node{
    private:
        // Input
        const bool process_input = false;
        bool mouse_over = false;
        bool mouse_down = false;
        bool consumed_mouse_down = false;

        // if any of these change, then `this.rect` and `parent.minSize` need to be updated -> parent->udpate_minSize(), update_rect()
        GUI2Rect2f _anchors;
        Vec2i _pos;
        Vec2i _size;
        Vec2i _minSize;

        // if this changes, then `children[i].rect` needs to be updated -> update_children_rects()
        GUI2Rect2i _rect;

        // if this changes, then `self.minSize` and `newChild.rect` needs to be updated
        std::vector<GUI2Node*> children;
        const bool leaf_node = false; // if leaf_node == true, then this node cannot have children
        GUI2Node* parent = nullptr;

        void set_minSize( Vec2i minSize );
        void update_rect(); // call parent->update_child_rect(this)

    protected:
        void update_minSize();

        std::vector<GUI2Node*> get_children();
        virtual void update_child_rect(GUI2Node* child);
        virtual void update_children_rects();
        virtual Vec2i calculate_minSize();
        virtual void on_rect_updated();
        
        virtual void on_mouse_enter();
        virtual void on_mouse_exit();
        virtual void on_mouse_over();

        virtual void on_mouse_down();
        virtual void on_mouse_up(); // also called when mouse is leaving rather than when released
        virtual void on_mouse_click();
        virtual void on_mouse_drag( const SDL_Event& event );

    public:
        bool is_mouse_over();
        bool is_mouse_down();

        void set_rect( GUI2Rect2i rect ); // this should not be called from the outside

        // constructors
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size );
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool leaf_node );
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool leaf_node, bool process_input );

        // getters
        GUI2Rect2f anchors();
        Vec2i pos();
        Vec2i size();
        Vec2i minSize();
        GUI2Rect2i rect();

        // setters
        void set_anchors( GUI2Rect2f anchors );
        void set_pos( Vec2i pos );
        void set_size( Vec2i size );

        // other
        GUI2Node* addChild( GUI2Node* node );
        void removeChild( GUI2Node* node );
        virtual void draw();
        bool onEvent( const SDL_Event& event );
};

class GUI2Panel : public GUI2Node{
    public:
        uint32_t bgColor;

        GUI2Panel ( GUI2Rect2f anchors, Vec2i pos, Vec2i size, uint32_t bgColor = 0xA0A0A0 );
        GUI2Panel ( uint32_t bgColor = 0xA0A0A0 );

        virtual void draw() override;
};

class GUI2Text : public GUI2Node{
    public:
        enum class Align { TOP_LEFT, TOP_CENTER, TOP_RIGHT, CENTER_LEFT, CENTER, CENTER_RIGHT, BOTTOM_LEFT, BOTTOM_CENTER, BOTTOM_RIGHT };

    private:
        std::string text;
        Align align = Align::TOP_LEFT;
        uint32_t fontSize = fontSizeDef;
        //bool allowOverflow = true; // TODO

        // calculated from the values above
        Vec2i textPos_;
        Vec2i textSize_;

        void recalculate_textPos();
        void recalculate_textSize();

    protected:
        Vec2i calculate_minSize() override;
        void on_rect_updated() override;

    public:
        uint32_t fontColor=0x000000;
        
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text );
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor );
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor, Align align );
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, Align align);
        GUI2Text( std::string text, Align align );

        void setText( std::string text );
        void setFontSize( uint32_t fontSize );
        void setAlign( Align align );

        virtual void draw() override;
};

class GUI2Vlist : public GUI2Node{
    public:
        enum class Align { LEFT, CENTER, RIGHT, STRETCH };

    private:
        unsigned int sepperation;
        Align align;

    protected:
        Vec2i calculate_minSize() override;

        virtual void update_child_rect(GUI2Node* child) override;
        virtual void update_children_rects() override;

    public:
        GUI2Vlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, unsigned int sepperation = 5, Align align = Align::STRETCH );
        GUI2Vlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, Align align );
        GUI2Vlist( unsigned int sepparation = 5, Align align = Align::STRETCH );

        void set_sepperation( unsigned int sepperation );
        unsigned int get_sepperation();

        void set_align(Align align);
        Align get_align();
};

class GUI2Hlist : public GUI2Node{
    public:
        enum class Align { TOP, CENTER, BOTTOM, STRETCH };

    private:
        unsigned int sepperation;
        Align align;

    protected:
        Vec2i calculate_minSize() override;

        virtual void update_child_rect(GUI2Node* child) override;
        virtual void update_children_rects() override;

    public:
        GUI2Hlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, unsigned int sepperation = 5, Align align = Align::STRETCH );
        GUI2Hlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, Align align );
        GUI2Hlist( unsigned int sepperation = 5, Align align = Align::STRETCH );

        void set_sepperation( unsigned int sepperation );
        unsigned int get_sepperation();

        void set_align(Align align);
        Align get_align();
};

class GUI2ButtonBase : public GUI2Node{
    protected:
        void on_mouse_click() override;

    public:
        const std::function<void()> command;

        GUI2ButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command );
};

class GUI2Button : public GUI2ButtonBase{
    private:
        GUI2Panel* panel;

    protected:
        virtual void on_mouse_enter() override;
        virtual void on_mouse_exit() override;
        virtual void on_mouse_down() override;
        virtual void on_mouse_up() override;

    public:
        uint32_t bgColor;
        uint32_t bgColorHover;
        uint32_t bgColorPressed;

        GUI2Button( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command );
        GUI2Button( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed );
};

class GUI2TextButton : public GUI2Button{
    private:
        GUI2Text* text_node;
        std::string text;
    
    public:
        GUI2TextButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, const std::function<void()>& command );
        GUI2TextButton( std::string text, const std::function<void()>& command );
};

class GUI2ToggleButtonBase : public GUI2Node {
    private:
        const std::function<void(bool)> command;
        bool active = false;
        bool* bound_bool = nullptr;

    protected:
        virtual void on_mouse_click() override;

    public:
        bool is_active();
        void bind_bool( bool* bound_bool );

        GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command );
        GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, bool* bound_bool );
        GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool );
};

class GUI2ToggleButton : public GUI2ToggleButtonBase {
    private:
        GUI2Panel* panel;

        bool hovering = false;
        bool pressed = false;

    protected:
        virtual void on_mouse_enter() override;
        virtual void on_mouse_exit() override;
        virtual void on_mouse_down() override;
        virtual void on_mouse_up() override;

        virtual void draw() override;

    public:
        uint32_t bgColor;
        uint32_t bgColorHover;
        uint32_t bgColorPressed;

        uint32_t bgColorActive;
        uint32_t bgColorActiveHover;
        uint32_t bgColorActivePressed;

        GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command );
        GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, bool* bound_bool );
        GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, bool* bound_bool, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed, uint32_t bgColorActive, uint32_t bgColorActiveHover, uint32_t bgColorActivePressed );

        GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool );
        GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed, uint32_t bgColorActive, uint32_t bgColorActiveHover, uint32_t bgColorActivePressed );

};

class GUI2ToggleTextButton : public GUI2ToggleButton {
    private:
        GUI2Text* text_node;
    
    public:
        GUI2ToggleTextButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool, std::string text );
        GUI2ToggleTextButton( bool* bound_bool, std::string text );
};

class GUI2Dragable : public GUI2Node {
    protected:
        virtual void on_mouse_drag( const SDL_Event& event ) override;
    
    public:
        GUI2Dragable( GUI2Rect2f anchors, Vec2i pos, Vec2i size );
};

template <class T>
class GUI2SliderT : public GUI2Node {
    private:
        GUI2Panel* bgPanel;
        GUI2Panel* fillPanel;

        const bool interactive = false;
        const bool limit_bound_value = false;

        T min, max;
        T value;
        T* bound_value;

        const std::function<void(T)> command;

    protected:
        virtual void draw() override;
        virtual void on_mouse_drag( const SDL_Event& event ) override;

    public:
        void set_min(T min);
        void set_max(T max);
        void set_value(T value);
        void bind_value(T* bound_value);

        GUI2SliderT( GUI2Rect2f anchors, Vec2i pos, Vec2i size, T min, T max, T* bound_value, const std::function<void(T)>& command, bool interactive = true );
};
using GUI2Slideri = GUI2SliderT<int>;
using GUI2Sliderf = GUI2SliderT<float>;

template <class T>
class GUI2TextSliderT : public GUI2SliderT<T> {
    private:
        GUI2Text* text;
        std::string format; // the string "$value" in format will be replaced with the std::to_string(value)
                            // for more specific formating, use "$value%.2f" to format the value with "%.2f" (to 2 decimal places)

        const std::function<void(T)> command;
        void on_value_update(T value);

    public:
        GUI2TextSliderT( GUI2Rect2f anchors, Vec2i pos, Vec2i size, T min, T max, T* bound_value, std::string format, const std::function<void(T)>& command, bool interactive = true );
        GUI2TextSliderT( T min, T max, T* bound_value, std::string format, const std::function<void(T)>& command = nullptr, bool interactive = true );
};
using GUI2TextSlideri = GUI2TextSliderT<int>;
using GUI2TextSliderf = GUI2TextSliderT<float>;


class GUI2 {public:
    private:
        GUI2Node root_node;

    public:
        GUI2();

        GUI2Node* addNode( GUI2Node* node );
        void draw( SDL_Window* window );
        bool onEvent( SDL_Event event, SDL_Window* window );
};
