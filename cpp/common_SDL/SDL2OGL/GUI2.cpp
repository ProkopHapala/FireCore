#include "GUI2.h"
#include "Draw.h"
#include "Draw2D.h"
#include "Vec2.h"
#include <SDL2/SDL_video.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <sys/types.h>

int GUI2_fontTex = 0;

// ==============================
//    class GUI2Node
// ==============================

// getters
GUI2Rect2f GUI2Node::anchors(){return _anchors;}
Vec2i GUI2Node::pos(){return _pos;}
Vec2i GUI2Node::size(){return _size;}
Vec2i GUI2Node::minSize(){return _minSize;}
GUI2Rect2i GUI2Node::rect(){return _rect;}

// setters
void GUI2Node::set_anchors(GUI2Rect2f anchors){
    if (_anchors == anchors) return;
    _anchors = anchors;

    if (parent != nullptr) {parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_pos(Vec2i pos){
    if (_pos == pos) return;
    _pos = pos;

    if (parent != nullptr) {parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_size(Vec2i size){
    if (_size == size) return;
    _size = size;

    if (parent != nullptr) {parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_minSize(Vec2i minSize){
    if (_minSize == minSize) return;
    _minSize = minSize;

    if (parent != nullptr) {parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_rect(GUI2Rect2i rect){
    if (_rect == rect) return;

    _rect = rect;
    on_rect_updated();
    update_children_rects();
}

// update/calculate functions
Vec2i GUI2Node::calculate_minSize(){
    Vec2i minSize = (Vec2i){0,0};
    for (auto c : children){
        Vec2i c_minSize = c->minSize() - c->size();
        if (c->anchors().size().x != 0){ c_minSize.x = c_minSize.x / c->anchors().size().x; }
        if (c->anchors().size().y != 0){ c_minSize.y = c_minSize.y / c->anchors().size().y; }

        minSize.x = std::min(minSize.x, c_minSize.x);
        minSize.y = std::min(minSize.y, c_minSize.y);
    }
    return minSize;
}
void GUI2Node::update_minSize(){
    Vec2i new_minSize = calculate_minSize();
    if (_minSize == new_minSize) return;
    _minSize = new_minSize;
    if (parent != nullptr) parent->update_minSize();
}
void GUI2Node::update_child_rect(GUI2Node* child){
    Vec2i min = _rect.min() + (Vec2i)(child->anchors().min()*(Vec2f)_rect.size()) + child->pos();
    Vec2i max = _rect.min() + (Vec2i)(child->anchors().max()*(Vec2f)_rect.size()) + child->size();
    child->set_rect(GUI2Rect2i(min, max));
}
void GUI2Node::update_children_rects(){
    for(auto& c:children){
        update_child_rect(c);
    }
}
void GUI2Node::update_rect(){
    if (parent != nullptr) { parent->update_child_rect(this); return; }

    // node is a root node, so rect = (pos, pos + size)
    _rect = GUI2Rect2i(_pos, _pos + _size);
    update_children_rects();
}
void GUI2Node::on_rect_updated(){}

// draw
void GUI2Node::draw(){
    for(auto& c:children){
        c->draw();
    }
}

// constructors
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size) : _anchors(anchors), _pos(pos), _size(size){
    update_minSize();
    update_rect();
}
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool leaf_node_): leaf_node(leaf_node_), _anchors(anchors), _pos(pos), _size(size){
    update_minSize();
    update_rect();
}


void GUI2Node::addChild(GUI2Node* node){
    if (leaf_node) {
        throw "ERROR: cannot add a child to leaf node";
    }
    if (node->parent != nullptr){
        printf("ERROR: node already has a parent - ignoring this addChild() call\n");
        return;
    }

    children.push_back(node);
    node->parent = this;
    update_minSize();
    update_child_rect(node);
}
void GUI2Node::removeChild(GUI2Node* node){
    if (leaf_node) {
        throw "ERROR: cannot remove a child from leaf node";
    }
    if (node->parent != this){
        printf("ERROR: node does not have this as a parent - ignoring this removeChild() call\n");
        return;
    }
    node->parent = nullptr;
    children.erase(std::find(children.begin(), children.end(), node));

    update_minSize();
}


// ==============================
//    class GUI2Panel
// =============================

void GUI2Panel::draw(){ // TODO: don't need to redraw every frame ?
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( rect().xmin, rect().ymin, rect().xmax, rect().ymax, true );

    GUI2Node::draw(); // call super (calls draw() on children)
}

GUI2Panel::GUI2Panel(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_): GUI2Node(anchors_, pos_, size_){}
GUI2Panel::GUI2Panel(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, uint32_t bgColor_): GUI2Node(anchors_, pos_, size_), bgColor(bgColor_){}


// ==============================
//    class GUI2Text
// ==============================

GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text ): GUI2Node(anchors, pos, size, true), text(text){recalculate_textSize();}
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor ): GUI2Node(anchors, pos, size, true), text(text), fontSize(fontSize), fontColor(fontColor){recalculate_textSize();}
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor, Align align ): GUI2Node(anchors, pos, size, true), text(text), fontSize(fontSize), fontColor(fontColor), align(align){recalculate_textSize();}
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, Align align): GUI2Node(anchors, pos, size, true), text(text), align(align){recalculate_textSize();}

void GUI2Text::recalculate_textSize(){
    unsigned int lineCount = 1;
    unsigned int lineLength = 0;
    unsigned int maxLineLength = 0;
    for(auto& c:text){
        if (c == '\n') {lineCount++; lineLength=0; }
        else lineLength++;
        if (lineLength > maxLineLength) maxLineLength = lineLength;
    }
    textSize_ = {maxLineLength*fontSize, lineCount*fontSize*2};
    GUI2Text::recalculate_textPos();
}

void GUI2Text::recalculate_textPos(){
    switch (align){
        case Align::TOP_LEFT:      textPos_ = {rect().xmin,                                   rect().ymax}; break;
        case Align::TOP_CENTER:    textPos_ = {rect().xmin + (rect().size().x-textSize_.x)/2, rect().ymax}; break;
        case Align::TOP_RIGHT:     textPos_ = {rect().xmax-textSize_.x,                       rect().ymax}; break;
        case Align::CENTER_LEFT:   textPos_ = {rect().xmin,                                   rect().ymin + (rect().size().y+textSize_.y)/2}; break;
        case Align::CENTER:        textPos_ = {rect().xmin + (rect().size().x-textSize_.x)/2, rect().ymin + (rect().size().y+textSize_.y)/2}; break;
        case Align::CENTER_RIGHT:  textPos_ = {rect().xmax-textSize_.x,                       rect().ymin + (rect().size().y+textSize_.y)/2}; break;
        case Align::BOTTOM_LEFT:   textPos_ = {rect().xmin,                                   rect().ymin+textSize_.y}; break;
        case Align::BOTTOM_CENTER: textPos_ = {rect().xmin + (rect().size().x-textSize_.x)/2, rect().ymin+textSize_.y}; break;
        case Align::BOTTOM_RIGHT:  textPos_ = {rect().xmax-textSize_.x,                       rect().ymin+textSize_.y}; break;
    }

    textPos_.y -= fontSize*2; // Draw2D::drawText() draws first line above textPos, but we want the first line to be below textPos
}

Vec2i GUI2Text::calculate_minSize(){
    GUI2Text::recalculate_textSize();
    return textSize_;
}

void GUI2Text::setText( std::string text_ ){
    text = text_;
    GUI2Text::recalculate_textSize();
}
void GUI2Text::setFontSize( uint32_t fontSize_ ){
    fontSize = fontSize_;
    GUI2Text::recalculate_textSize();
}
void GUI2Text::setAlign( Align align_ ){
    align = align_;
    GUI2Text::recalculate_textPos();
}

void GUI2Text::on_rect_updated(){
    GUI2Text::recalculate_textPos();
}


void GUI2Text::draw(){
    Draw  ::setRGB( fontColor );
    Draw2D::drawText( text.c_str(), (Vec2d)textPos_, (Vec2d)textSize_, GUI2_fontTex, fontSize );

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
