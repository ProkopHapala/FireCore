#include "GUI2.h"
#include "Draw.h"
#include "Draw2D.h"
#include "Vec2.h"
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_video.h>
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <sys/types.h>

int GUI2_fontTex = 0;

// ==============================
//    class GUI2Node
// ==============================

// constructors
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool leaf_node, bool process_input):
    leaf_node(leaf_node),
    process_input(process_input),
    _anchors(anchors),
    _pos(pos),
    _size(size)
    {
        update_minSize();
        update_rect();
    }
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool leaf_node_):
    GUI2Node(anchors, pos, size, leaf_node_, false){}
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size):
    GUI2Node(anchors, pos, size, false, false){}

// getters
GUI2Rect2f GUI2Node::anchors(){return _anchors;}
Vec2i GUI2Node::pos(){return _pos;}
Vec2i GUI2Node::size(){return _size;}
Vec2i GUI2Node::minSize(){return _minSize;}
GUI2Rect2i GUI2Node::rect(){return _rect;}
bool GUI2Node::is_mouse_over(){return mouse_over;}
bool GUI2Node::is_mouse_down(){return mouse_down;}

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
std::vector<GUI2Node*> GUI2Node::get_children(){
    return children;
}

// update/calculate functions
Vec2i GUI2Node::calculate_minSize(){
    Vec2i minSize = (Vec2i){0,0};
    for (auto c : children){
        Vec2i c_minSize = c->minSize() - c->size();
        if (c->anchors().size().x != 0){ c_minSize.x = c_minSize.x / c->anchors().size().x; }
        if (c->anchors().size().y != 0){ c_minSize.y = c_minSize.y / c->anchors().size().y; }

        minSize.x = std::max(minSize.x, c_minSize.x);
        minSize.y = std::max(minSize.y, c_minSize.y);
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


GUI2Node* GUI2Node::addChild(GUI2Node* node){
    if (leaf_node) {
        throw "ERROR: cannot add a child to leaf node";
    }
    if (node->parent != nullptr){
        printf("ERROR: node already has a parent - ignoring this addChild() call\n");
        return nullptr;
    }

    children.push_back(node);
    node->parent = this;
    update_minSize();
    update_child_rect(node);

    return node;
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

bool GUI2Node::onEvent(const SDL_Event& event){
    for(auto& c:children){
        if (c->onEvent(event)) return true;
    }

    if (!process_input) return false;

    switch (event.type){
        case SDL_MOUSEMOTION:
            if ( rect().contains((Vec2i){event.motion.x, event.motion.y}) ){
                if (!mouse_over) on_mouse_enter();
                mouse_over = true;
                on_mouse_over();
            }else{
                if (mouse_down) { mouse_down = false; on_mouse_up(); }
                if (mouse_over) on_mouse_exit();
                mouse_over = false;
            }
            break;

        case SDL_MOUSEBUTTONDOWN:
            switch (event.button.button){
                case SDL_BUTTON_LEFT:
                    if ( rect().contains((Vec2i){event.button.x, event.button.y}) ){
                            on_mouse_down();
                            mouse_down = true;
                            consumed_mouse_down = true;

                            return true;
                        }
                        break;
            }break;
        
        case SDL_MOUSEBUTTONUP:
            switch (event.button.button){
                case SDL_BUTTON_LEFT:
                    if ( rect().contains((Vec2i){event.button.x, event.button.y}) ){
                            if (mouse_down) { on_mouse_click(); on_mouse_up(); }
                            mouse_down = false;
                        }
                    if (consumed_mouse_down) { consumed_mouse_down = false; return true; }
                    break;
            }break;
        default: break;
    }

    return false;
}
void GUI2Node::on_mouse_enter(){}
void GUI2Node::on_mouse_over(){}
void GUI2Node::on_mouse_exit(){}

void GUI2Node::on_mouse_down(){}
void GUI2Node::on_mouse_up(){}
void GUI2Node::on_mouse_click(){}


// ==============================
//    class GUI2Panel
// =============================

GUI2Panel::GUI2Panel(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_):
    GUI2Node(anchors_, pos_, size_){}
GUI2Panel::GUI2Panel(GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, uint32_t bgColor_):
    GUI2Node(anchors_, pos_, size_),
    bgColor(bgColor_){}


void GUI2Panel::draw(){ // TODO: don't need to redraw every frame ?
    Draw  ::setRGB( bgColor );
    Draw2D::drawRectangle ( rect().xmin, rect().ymin, rect().xmax, rect().ymax, true );

    GUI2Node::draw(); // call super (calls draw() on children)
}


// ==============================
//    class GUI2Text
// ==============================

GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text ):
    GUI2Node(anchors, pos, size, true),
    text(text)
    { update_minSize(); }
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor ):
    GUI2Node(anchors, pos, size, true),
    text(text),
    fontSize(fontSize),
    fontColor(fontColor)
    { update_minSize(); }
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor, Align align ):
    GUI2Node(anchors, pos, size, true),
    text(text),
    fontSize(fontSize),
    fontColor(fontColor),
    align(align)
    { update_minSize(); }
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, Align align):
    GUI2Node(anchors, pos, size, true),
    text(text),
    align(align)
    { update_minSize(); }

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
//    class GUI2Vlist
// ==============================

GUI2Vlist::GUI2Vlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_ ):
    GUI2Node(anchors_, pos_, size_){}
GUI2Vlist::GUI2Vlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, unsigned int sepperation_ ):
    GUI2Node(anchors_, pos_, size_),
    sepperation(sepperation_){}
GUI2Vlist::GUI2Vlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, Align align_ ):
    GUI2Node(anchors_, pos_, size_),
    align(align_){}
GUI2Vlist::GUI2Vlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, unsigned int sepperation_, Align align_ ):
    GUI2Node(anchors_, pos_, size_),
    sepperation(sepperation_),
    align(align_){}

void GUI2Vlist::set_sepperation( unsigned int sepperation_ ){
    sepperation = sepperation_;
    update_minSize();
    update_children_rects();
}
unsigned int GUI2Vlist::get_sepperation(){
    return sepperation;
}
void GUI2Vlist::set_align( Align align_ ){
    align = align_;
    update_children_rects();
}
GUI2Vlist::Align GUI2Vlist::get_align(){
    return align;
}

Vec2i GUI2Vlist::calculate_minSize(){
    Vec2i minSize = {0,0};
    for(auto& child:get_children()){
        Vec2i childMinSize = child->minSize();
        minSize.x = std::max(minSize.x, childMinSize.x);
        minSize.y += childMinSize.y + sepperation;
    }
    minSize.y -= sepperation;
    return minSize;
}

void GUI2Vlist::update_child_rect(GUI2Node* child){
    update_children_rects();
}
void GUI2Vlist::update_children_rects(){
    unsigned int ymax = rect().ymax;

    for(auto& child:get_children()){
        GUI2Rect2i R;
        R.ymax = ymax;
        R.ymin = ymax - child->minSize().y;

        switch (align){
            case Align::STRETCH:  R.xmin = rect().xmin;
                                  R.xmax = rect().xmax; break;

            case Align::LEFT:     R.xmin = rect().xmin; 
                                  R.xmax = rect().xmin + child->minSize().x; break;

            case Align::CENTER:   R.xmin = rect().xmin + (rect().size().x-child->minSize().x)/2;
                                  R.xmax = R.xmin + child->minSize().x; break;

            case Align::RIGHT:    R.xmin = rect().xmax - child->minSize().x;
                                  R.xmax = rect().xmax; break;
        }

        child->set_rect(R);
        ymax -= child->minSize().y + sepperation;
    }
}


// ==============================
//    class GUI2Hlist
// ==============================

GUI2Hlist::GUI2Hlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_ ):
    GUI2Node(anchors_, pos_, size_){}
GUI2Hlist::GUI2Hlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, unsigned int sepperation_ ):
    GUI2Node(anchors_, pos_, size_),
    sepperation(sepperation_){}
GUI2Hlist::GUI2Hlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, Align align_ ):
    GUI2Node(anchors_, pos_, size_),
    align(align_){}
GUI2Hlist::GUI2Hlist( GUI2Rect2f anchors_, Vec2i pos_, Vec2i size_, unsigned int sepperation_, Align align_ ):
    GUI2Node(anchors_, pos_, size_),
    sepperation(sepperation_),
    align(align_){}

void GUI2Hlist::set_sepperation( unsigned int sepperation_ ){
    sepperation = sepperation_;
    update_minSize();
    update_children_rects();
}
unsigned int GUI2Hlist::get_sepperation(){
    return sepperation;
}
void GUI2Hlist::set_align( Align align_ ){
    align = align_;
    update_children_rects();
}
GUI2Hlist::Align GUI2Hlist::get_align(){
    return align;
}

Vec2i GUI2Hlist::calculate_minSize(){
    Vec2i minSize = {0,0};
    for(auto& child:get_children()){
        Vec2i childMinSize = child->minSize();
        minSize.y = std::max(minSize.y, childMinSize.y);
        minSize.x += childMinSize.x + sepperation;
    }
    minSize.x -= sepperation;
    return minSize;
}

void GUI2Hlist::update_child_rect(GUI2Node* child){
    update_children_rects();
}
void GUI2Hlist::update_children_rects(){
    unsigned int xmin = rect().xmin;

    for(auto& child:get_children()){
        GUI2Rect2i R;
        R.xmin = xmin;
        R.xmax = xmin + child->minSize().x;

        switch (align){
            case Align::STRETCH:  R.ymin = rect().ymin;
                                  R.ymax = rect().ymax; break;

            case Align::TOP:      R.ymin = rect().ymax - child->minSize().y; 
                                  R.ymax = rect().ymax; break;

            case Align::CENTER:   R.ymin = rect().ymin + (rect().size().y-child->minSize().y)/2;
                                  R.ymax = R.ymin + child->minSize().y; break;

            case Align::BOTTOM:   R.ymin = rect().ymin;
                                  R.ymax = rect().ymin + child->minSize().y; break;
        }

        child->set_rect(R);
        xmin += child->minSize().x + sepperation;
    }
}

// ==============================
//    class GUI2ButtonBase
// ==============================

GUI2ButtonBase::GUI2ButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command ):
    GUI2Node(anchors, pos, size, false, true),
    command(command){}

void GUI2ButtonBase::on_mouse_click(){
    std::invoke(command);
}

// ==============================
//    class GUI2Button
// ==============================

GUI2Button::GUI2Button( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed ):
    GUI2ButtonBase(anchors, pos, size, command),
    bgColor(bgColor),
    bgColorHover(bgColorHover),
    bgColorPressed(bgColorPressed)
    {
        panel = (GUI2Panel*)addChild(new GUI2Panel({0, 0, 1, 1}, {0, 0}, {0, 0}));
        panel->bgColor = bgColor;
    }

GUI2Button::GUI2Button( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command ):
    GUI2Button(anchors, pos, size, command, 0x808080, 0x787878, 0x606060){}

void GUI2Button::on_mouse_enter(){
    panel->bgColor = bgColorHover;
}
void GUI2Button::on_mouse_exit(){
    panel->bgColor = bgColor;
}
void GUI2Button::on_mouse_down(){
    panel->bgColor = bgColorPressed;
}
void GUI2Button::on_mouse_up(){
    panel->bgColor = bgColorHover;
}

// ==============================
//    class GUI2ToggleButtonBase
// ==============================

GUI2ToggleButtonBase::GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command ):
    GUI2Node(anchors, pos, size, false, true),
    command(command){}

void GUI2ToggleButtonBase::on_mouse_click(){
    active = !active;
    command(active);
}

bool GUI2ToggleButtonBase::is_active(){
    return active;
}

// ==============================
//    class GUI2ToggleButton
// ==============================

GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed, uint32_t bgColorActive, uint32_t bgColorActiveHover, uint32_t bgColorActivePressed ):
    GUI2ToggleButtonBase(anchors, pos, size, command),
    bgColor(bgColor),
    bgColorHover(bgColorHover),
    bgColorPressed(bgColorPressed),
    bgColorActive(bgColorActive),
    bgColorActiveHover(bgColorActiveHover),
    bgColorActivePressed(bgColorActivePressed)
    {
        panel = (GUI2Panel*)addChild(new GUI2Panel({0, 0, 1, 1}, {0, 0}, {0, 0}));
        panel->bgColor = bgColor;
    }
GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command ):
    GUI2ToggleButton(anchors, pos, size, command, 0xA0A0A0, 0x989898, 0x808080, 0x00FF00, 0x00E800, 0x00D000){}

void GUI2ToggleButton::on_mouse_enter(){
    panel->bgColor = is_active() ? bgColorActiveHover : bgColorHover;
}
void GUI2ToggleButton::on_mouse_exit(){
    panel->bgColor = is_active() ? bgColorActive : bgColor;
}
void GUI2ToggleButton::on_mouse_down(){
    panel->bgColor = is_active() ? bgColorActivePressed : bgColorPressed;
}
void GUI2ToggleButton::on_mouse_up(){
    panel->bgColor = is_active() ? bgColorActiveHover : bgColorHover;
}



// ==============================
//    class GUI2
// ==============================

GUI2::GUI2():
    root_node(GUI2Node( GUI2Rect2f(0,0,1,1), (Vec2i){0,0}, (Vec2i){0,0} )) {}

GUI2Node* GUI2::addNode(GUI2Node* node){
    return root_node.addChild(node);
}
void GUI2::draw( SDL_Window* window){
    int width, height;
    SDL_GetWindowSize(window, &width, &height);

    root_node.set_pos((Vec2i){0,0});
    root_node.set_size((Vec2i){width, height});

    glDisable   ( GL_LIGHTING    );
    glDisable   ( GL_DEPTH_TEST  );
    glShadeModel( GL_FLAT        );
    root_node.draw();
}

bool GUI2::onEvent( SDL_Event event, SDL_Window* window ){
    int width, height;
    SDL_GetWindowSize(window, &width, &height);

    switch (event.type){
        case (SDL_MOUSEMOTION):
            event.motion.y = height - event.motion.y;
            event.motion.yrel *= -1;
            break;
        case (SDL_MOUSEBUTTONUP):
        case (SDL_MOUSEBUTTONDOWN):
            event.button.y = height - event.button.y;
            break;
    }

    return root_node.onEvent(event);
}
