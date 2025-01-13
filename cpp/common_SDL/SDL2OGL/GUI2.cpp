#include "GUI2.h"
#include "Draw.h"
#include "Draw2D.h"
#include "Vec2.h"
#include <SDL2/SDL_video.h>
#include <cstdint>
#include <cstdio>
#include <sys/types.h>

int GUI2_fontTex = 0;

// ==============================
//    class GUI2Node
// ==============================

void GUI2Node::draw(){
    for(auto& c:children){
        c->draw();
    }
}

void GUI2Node::update_child_rect(GUI2Node* child){
    Vec2i min = rect.min() + (Vec2i)(child->anchors.min()*(Vec2f)rect.size()) + child->pos;
    Vec2i max = rect.min() + (Vec2i)(child->anchors.max()*(Vec2f)rect.size()) + child->size;
    child->set_rect(GUI2Rect2i(min, max));
}

void GUI2Node::update_children_rects(){
    for(auto& c:children){
        update_child_rect(c);
    }
}

void GUI2Node::set_rect(GUI2Rect2i rect_){
    if (rect == rect_) return; // skip updating the rect if nothing changed
    rect = rect_;
    update_children_rects();
}

void GUI2Node::update_rect(){
    if (parent) { parent->update_child_rect(this); return; }

    // node is a root node, so rect = (pos, pos + size)
    rect = GUI2Rect2i(pos, pos + size);
    update_children_rects();
}

void GUI2Node::set_pos(Vec2i pos_){
    if (pos == pos_) return; // skip updating the rect if nothing changed
    pos = pos_;
    update_rect();
}
void GUI2Node::set_size(Vec2i size_){
    if (size == size_) return; // skip updating the rect if nothing changed
    size = size_;
    update_rect();
}

GUI2Node::GUI2Node(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_){
    anchors = anchors_;
    pos = pos_;
    size = size_;
    update_rect();
}

void GUI2Node::addChild(GUI2Node* node){
    if (node->parent != nullptr){
        printf("ERROR: node already has a parent - ignoring this addChild() call\n");
        return;
    }

    children.push_back(node);
    node->parent = this;
}


// ==============================
//    class GUI2Panel
// =============================

void GUI2Panel::draw(){ // TODO: don't need to redraw every frame
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( rect.xmin, rect.ymin, rect.xmax, rect.ymax, true );

    GUI2Node::draw(); // call super (calls draw() on children)
}

GUI2Panel::GUI2Panel(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_): GUI2Node(anchors_, pos_, size_){}
GUI2Panel::GUI2Panel(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, uint32_t bgColor_): GUI2Node(anchors_, pos_, size_), bgColor(bgColor_){}


// ==============================
//    class GUI2Text
// ==============================

GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text ): GUI2Node(anchors, pos, size), text(text){recalculateTextPos();}
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor ): GUI2Node(anchors, pos, size), text(text), fontSize(fontSize), fontColor(fontColor){recalculateTextPos();}
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor, Align align ): GUI2Node(anchors, pos, size), text(text), fontSize(fontSize), fontColor(fontColor), align(align){recalculateTextPos();}
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, Align align): GUI2Node(anchors, pos, size), text(text), align(align){recalculateTextPos();}

void GUI2Text::recalculateTextPos(){ // TODO: recalculate min_size
    unsigned int lineCount = 1;
    unsigned int lineLength = 0;
    unsigned int maxLineLength = 0;
    for(auto& c:text){
        if (c == '\n') {lineCount++; lineLength=0; }
        else lineLength++;
        if (lineLength > maxLineLength) maxLineLength = lineLength;
    }
    textSize_ = {maxLineLength*fontSize, lineCount*fontSize*2};

    switch (align){
        case Align::TOP_LEFT:      textPos_ = {rect.xmin,                                 rect.ymax}; break;
        case Align::TOP_CENTER:    textPos_ = {rect.xmin + (rect.size().x-textSize_.x)/2, rect.ymax}; break;
        case Align::TOP_RIGHT:     textPos_ = {rect.xmax-textSize_.x,                     rect.ymax}; break;
        case Align::CENTER_LEFT:   textPos_ = {rect.xmin,                                 rect.ymin + (rect.size().y+textSize_.y)/2}; break;
        case Align::CENTER:        textPos_ = {rect.xmin + (rect.size().x-textSize_.x)/2, rect.ymin + (rect.size().y+textSize_.y)/2}; break;
        case Align::CENTER_RIGHT:  textPos_ = {rect.xmax-textSize_.x,                     rect.ymin + (rect.size().y+textSize_.y)/2}; break;
        case Align::BOTTOM_LEFT:   textPos_ = {rect.xmin,                                 rect.ymin+textSize_.y}; break;
        case Align::BOTTOM_CENTER: textPos_ = {rect.xmin + (rect.size().x-textSize_.x)/2, rect.ymin+textSize_.y}; break;
        case Align::BOTTOM_RIGHT:  textPos_ = {rect.xmax-textSize_.x,                     rect.ymin+textSize_.y}; break;
    }

    textPos_.y -= fontSize*2; // Draw2D::drawText() draws first line above textPos, but we want the first line to be below textPos
}

void GUI2Text::setText( std::string text_ ){
    text = text_;
    GUI2Text::recalculateTextPos();
}
void GUI2Text::setFontSize( uint32_t fontSize_ ){
    fontSize = fontSize_;
    GUI2Text::recalculateTextPos();
}
void GUI2Text::setAlign( Align align_ ){
    align = align_;
    GUI2Text::recalculateTextPos();
}

void GUI2Text::set_rect(GUI2Rect2i rect_){
    if (rect_ == rect) return;
    GUI2Node::set_rect(rect_);
    GUI2Text::recalculateTextPos();
}


void GUI2Text::draw(){
    Draw  ::setRGB( fontColor );
    Draw2D::drawText( text.c_str(), textPos_, textSize_, GUI2_fontTex, fontSize );

    GUI2Node::draw(); // call super (calls draw() on children)
}



// ==============================
//    class GUI2
// ==============================

void GUI2::addNode(GUI2Node* node){
    root_node.addChild(node);
}
void GUI2::draw(SDL_Window* window){
    int width, height;
    SDL_GetWindowSize(window, &width, &height);

    root_node.set_pos((Vec2i){0,0});
    root_node.set_size((Vec2i){width, height});

    glDisable   ( GL_LIGHTING    );
    glDisable   ( GL_DEPTH_TEST  );
    glShadeModel( GL_FLAT        );
    root_node.draw();
}
