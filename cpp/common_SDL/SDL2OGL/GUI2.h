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
    public:
        enum ExpandMode {UP_RIGHT, DOWN_LEFT, DOWN_RIGHT, UP_LEFT};
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
        ExpandMode _expandMode = UP_RIGHT;

        // if this changes, then `children[i].rect` needs to be updated -> update_children_rects()
        GUI2Rect2i _rect;

        // if this changes, then `self.minSize` and `newChild.rect` needs to be updated
        std::vector<GUI2Node*> children;
        GUI2Node* _parent = nullptr;

        bool _active = true;

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

        virtual void draw_self();

    public:
        bool is_mouse_over();
        bool is_mouse_down();

        void set_rect( GUI2Rect2i rect ); // this should not be called from the outside

        // constructors
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size );
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool process_input );

        // getters
        GUI2Rect2f anchors();
        Vec2i pos();
        Vec2i size();
        Vec2i minSize();
        GUI2Rect2i rect();
        GUI2Node* parent();
        bool active();
        ExpandMode expandMode();

        // setters
        void set_anchors( GUI2Rect2f anchors );
        void set_pos( Vec2i pos );
        void set_size( Vec2i size );
        void set_active( bool active );
        void set_expandMode( ExpandMode expandMode );

        // other
        GUI2Node* addChild( GUI2Node* node );
        void removeChild( GUI2Node* node );
        void draw();
        virtual bool onEvent( const SDL_Event& event );
};

class GUI2Panel : public GUI2Node{
    protected:
        virtual void draw_self() override;

    public:
        uint32_t bgColor;

        GUI2Panel ( GUI2Rect2f anchors, Vec2i pos, Vec2i size, uint32_t bgColor = 0xA0A0A0 );
        GUI2Panel ( uint32_t bgColor = 0xA0A0A0 );
};

class GUI2Text : public GUI2Node{
    public:
        enum class Align { TOP_LEFT, TOP_CENTER, TOP_RIGHT, CENTER_LEFT, CENTER, CENTER_RIGHT, BOTTOM_LEFT, BOTTOM_CENTER, BOTTOM_RIGHT };

    private:
        using GUI2Node::addChild;

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
        void draw_self() override;

    public:
        uint32_t fontColor=0x000000; // TODO: fix - currently does nothing
        
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text );
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor );
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor, Align align );
        GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, Align align);
        GUI2Text( std::string text, Align align );

        void setText( std::string text );
        void setFontSize( uint32_t fontSize );
        void setAlign( Align align );
};

class GUI2Vlist : public GUI2Node{
    public:
        enum class Align { LEFT, CENTER, RIGHT, STRETCH };

    private:
        unsigned int min_sepperation;
        Align align;

    protected:
        Vec2i calculate_minSize() override;

        void update_child_rect(GUI2Node* child) override;
        void update_children_rects() override;

    public:
        GUI2Vlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, unsigned int min_sepperation = 5, Align align = Align::STRETCH );
        GUI2Vlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, Align align );
        GUI2Vlist( unsigned int sepparation = 5, Align align = Align::STRETCH );

        void set_sepperation( unsigned int min_sepperation );
        unsigned int get_sepperation();

        void set_align(Align align);
        Align get_align();
};

class GUI2Hlist : public GUI2Node{
    public:
        enum class Align { TOP, CENTER, BOTTOM, STRETCH };

    private:
        unsigned int min_sepperation;
        Align align;

    protected:
        Vec2i calculate_minSize() override;

        void update_child_rect(GUI2Node* child) override;
        void update_children_rects() override;

    public:
        GUI2Hlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, unsigned int min_sepperation = 5, Align align = Align::STRETCH );
        GUI2Hlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, Align align );
        GUI2Hlist( unsigned int min_sepperation = 5, Align align = Align::STRETCH );

        void set_sepperation( unsigned int min_sepperation );
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
        GUI2ButtonBase( const std::function<void()>& command );
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
        GUI2Button( const std::function<void()>& command );
};

class GUI2TextButton : public GUI2Button{
    private:
        GUI2Text* text_node;
        std::string text;
        int _text_padding = 0;
    
    public:
        int text_padding();
        void set_text_padding(int text_padding);

        GUI2TextButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, const std::function<void()>& command );
        GUI2TextButton( std::string text, const std::function<void()>& command );
};

class GUI2ToggleButtonBase : public GUI2Node {
    private:
        std::function<void(bool)> command;
        bool active = false;
        bool* bound_bool = nullptr;

    protected:
        virtual void on_mouse_click() override;

    public:
        bool is_active();
        void set(bool active);
        void bind_bool( bool* bound_bool );
        void set_command( const std::function<void(bool)>& command );

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

        virtual void draw_self() override;

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

        GUI2ToggleButton( const std::function<void(bool)>& command );
};

class GUI2ToggleTextButton : public GUI2ToggleButton {
    private:
        GUI2Text* text_node;
        int _text_padding = 0;

    public:
        int text_padding();
        void set_text_padding(int text_padding);

        GUI2ToggleTextButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool, std::string text, std::function<void(bool)> command );
        GUI2ToggleTextButton( bool* bound_bool, std::string text );
        GUI2ToggleTextButton( std::string text, std::function<void(bool)> command );
};

class GUI2Dragable : public GUI2Node {
    protected:
        void on_mouse_drag( const SDL_Event& event ) override;
    
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
        void draw_self() override;
        void on_mouse_drag( const SDL_Event& event ) override;

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

class GUI2FloatingMenu : public GUI2Node {
    private:
        using GUI2Node::addChild;
        using GUI2Node::removeChild;

        class MenuButton : public GUI2Button {
            private:
                using GUI2Node::addChild;
                using GUI2Node::removeChild;

                int submenu_idx = -1;
                GUI2FloatingMenu* menu = nullptr;

            public:
                MenuButton( std::string text, std::string keybind_text, const std::function<void()> command );
                MenuButton( std::string text, std::string keybind_text, int submenu_idx, GUI2FloatingMenu* menu );
        };

        class SubmenuButton : public GUI2ToggleButton {
            private:
                using GUI2Node::addChild;
                using GUI2Node::removeChild;

                int submenu_idx;
                GUI2FloatingMenu* menu;

                void on_mouse_enter() override;

            public:
                SubmenuButton( std::string text, int submenu_idx, GUI2FloatingMenu* menu );
        };

        GUI2Vlist* vlist = nullptr;
        GUI2Panel* bgPanel = nullptr;

        std::vector<GUI2FloatingMenu*> submenus;
        std::vector<SubmenuButton*> submenuButtons;
        int current_submenu = -1;

        void update_child_rect(GUI2Node* node) override;
        void update_children_rects() override;
        
    protected:
        bool onEvent(const SDL_Event& event) override;
    
    public:
        void open_submenu(int idx);
        void close_submenu();

        bool collapse_self;
    
        GUI2FloatingMenu( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool collapse_self = true );
        GUI2FloatingMenu( bool collapse_self = true );

        GUI2FloatingMenu* add_submenu(std::string name);
        GUI2Node* add_item(GUI2Node* item);
        GUI2Button* add_button(std::string text, std::string keybind_text, const std::function<void()> command);
};

class GUI2Toolbar : public GUI2Node {
    private:
        using GUI2Node::addChild;
        using GUI2Node::removeChild;

        class CustomTextButton : public GUI2ToggleTextButton{
            private:
                int tabIdx;
                GUI2Toolbar* toolbar;

            protected:
                void on_mouse_enter() override;

            public:
                CustomTextButton( std::string text, const std::function<void(bool)> command, int tabIdx, GUI2Toolbar* toolbar );
        };

        GUI2Panel* bgPanel;
        GUI2Hlist* buttonsHlist;
        
        std::vector<CustomTextButton*> buttons;
        std::vector<GUI2FloatingMenu*> tabs;

        int openTabIdx = -1;

    protected:
        void update_child_rect(GUI2Node* node) override;
        void update_children_rects() override;

    public:
        void open_tab(int idx);
        void close_tab(int idx);

        GUI2Toolbar( );
        GUI2FloatingMenu* addTab( std::string name );

        bool onEvent( const SDL_Event& event ) override;
};

class GUI2 {public:
    private:
        GUI2Node root_node;

    public:
        GUI2();

        GUI2Node* addNode( GUI2Node* node );
        void draw( SDL_Window* window );
        bool onEvent( SDL_Event event, SDL_Window* window );
};
