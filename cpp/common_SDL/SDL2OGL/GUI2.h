#include <SDL2/SDL_video.h>
#include <cmath>
#include <cstdint>
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

    bool operator==( const REC& r ){ return (xmin==r.xmin && ymin==r.ymin && xmax==r.xmax && ymax==r.ymax); }
    bool operator!=( const REC& r ){ return (xmin!=r.xmin || ymin!=r.ymin || xmax!=r.xmax || ymax!=r.ymax); }
};

using GUI2Rect2i = GUI2Rect2T<int>;
using GUI2Rect2f = GUI2Rect2T<float>;

extern int GUI2_fontTex;

class GUI2Node{
    private:
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
    protected:
        void set_rect( GUI2Rect2i rect );
        void set_minSize( Vec2i minSize );
        void update_rect(); // call parent->update_child_rect(this)
        void update_minSize();

        virtual void update_child_rect(GUI2Node* child);
        virtual void update_children_rects();
        virtual Vec2i calculate_minSize();
        virtual void on_rect_updated();
        
    public:
        // constructors
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size );
        GUI2Node( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool leaf_node );

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
        void addChild( GUI2Node* node );
        void removeChild( GUI2Node* node );
        virtual void draw();
};

class GUI2Panel : public GUI2Node{
    public:
        uint32_t bgColor=0xA0A0A0;

        GUI2Panel ( GUI2Rect2f anchors, Vec2i pos, Vec2i size );
        GUI2Panel ( GUI2Rect2f anchors, Vec2i pos, Vec2i size, uint32_t bgColor );

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

        void setText( std::string text );
        void setFontSize( uint32_t fontSize );
        void setAlign( Align align );

        virtual void draw() override;
};

class GUI2 {public:
    private:
        GUI2Node root_node;

    public:
        GUI2() : root_node(GUI2Node( GUI2Rect2f(0,0,1,1), (Vec2i){0,0}, (Vec2i){0,0} )) {}

        void addNode(GUI2Node* node);
        void draw(SDL_Window* window);
};
