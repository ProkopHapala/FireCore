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
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool process_input):
    process_input(process_input),
    _anchors(anchors),
    _pos(pos),
    _size(size)
    {
        update_minSize();
        update_rect();
    }
GUI2Node::GUI2Node(GUI2Rect2f anchors, Vec2i pos, Vec2i size):
    GUI2Node(anchors, pos, size, false){}

// getters
GUI2Rect2f GUI2Node::anchors(){return _anchors;}
Vec2i GUI2Node::pos(){return _pos;}
Vec2i GUI2Node::size(){return _size;}
Vec2i GUI2Node::minSize(){return _minSize;}
GUI2Rect2i GUI2Node::rect(){return _rect;}
GUI2Node* GUI2Node::parent(){return _parent;}
GUI2Node::ExpandMode GUI2Node::expandMode(){return _expandMode;}
bool GUI2Node::active(){return _active;}
bool GUI2Node::is_mouse_over(){return mouse_over;}
bool GUI2Node::is_mouse_down(){return mouse_down;}
std::vector<GUI2Node*> GUI2Node::get_children(){return children;}

// setters
void GUI2Node::set_anchors(GUI2Rect2f anchors){
    if (_anchors == anchors) return;
    _anchors = anchors;

    if (_parent != nullptr) {_parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_pos(Vec2i pos){
    if (_pos == pos) return;
    _pos = pos;

    if (_parent != nullptr) {_parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_size(Vec2i size){
    if (_size == size) return;
    _size = size;

    if (_parent != nullptr) {_parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_minSize(Vec2i minSize){
    if (_minSize == minSize) return;
    _minSize = minSize;

    if (_parent != nullptr) {_parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_expandMode(ExpandMode expandMode){
    if (_expandMode == expandMode) return;
    _expandMode = expandMode;

    if (_parent != nullptr) {_parent->update_minSize();}
    update_rect();
}
void GUI2Node::set_active(bool active){
    _active = active;
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
        if (c->anchors().size().x != 0){ c_minSize.x = c_minSize.x / c->anchors().size().x; } else { c_minSize.x = 0; }
        if (c->anchors().size().y != 0){ c_minSize.y = c_minSize.y / c->anchors().size().y; } else { c_minSize.y = 0; }

        minSize.x = std::max(minSize.x, c_minSize.x);
        minSize.y = std::max(minSize.y, c_minSize.y);
    }
    return minSize;
}
void GUI2Node::update_minSize(){
    Vec2i new_minSize = calculate_minSize();
    if (_minSize == new_minSize) return;
    _minSize = new_minSize;
    if (_parent != nullptr) _parent->update_minSize();
    update_rect(); // TODO: does this lead to O(n^2) complexity?
}
void GUI2Node::update_child_rect(GUI2Node* child){
    Vec2i min = _rect.min() + (Vec2i)(child->anchors().min()*(Vec2f)_rect.size()) + child->pos();
    Vec2i max = _rect.min() + (Vec2i)(child->anchors().max()*(Vec2f)_rect.size()) + child->pos() + child->size();

    if ((max - min).x < child->minSize().x){
        if (child->_expandMode == UP_RIGHT || child->_expandMode == DOWN_RIGHT) max.x = min.x + child->minSize().x;
        else min.x = max.x - child->minSize().x;
    }
    if ((max - min).y < child->minSize().y){
        if (child->_expandMode == UP_RIGHT || child->_expandMode == UP_LEFT) max.y = min.y + child->minSize().y;
        else min.y = max.y - child->minSize().y;
    }

    child->set_rect(GUI2Rect2i(min, max));
}
void GUI2Node::update_children_rects(){
    for(auto& c:children){
        update_child_rect(c);
    }
}
void GUI2Node::update_rect(){
    if (_parent != nullptr) { _parent->update_child_rect(this); return; }

    // node is a root node, so rect = (pos, pos + size)
    _rect = GUI2Rect2i(_pos, _pos + _size);
    update_children_rects();
}
void GUI2Node::on_rect_updated(){}

// draw
void GUI2Node::draw(){
    if (!_active) return;
    draw_self();
    for(auto& c:children){
        c->draw();
    }
}
void GUI2Node::draw_self(){}



GUI2Node* GUI2Node::addChild(GUI2Node* node){
    if (node->_parent != nullptr){
        printf("ERROR: node already has a parent - ignoring this addChild() call\n");
        return nullptr;
    }

    children.push_back(node);
    node->_parent = this;
    update_minSize();
    update_child_rect(node);

    return node;
}
void GUI2Node::removeChild(GUI2Node* node){
    if (node->_parent != this){
        printf("ERROR: node does not have this as a parent - ignoring this removeChild() call\n");
        return;
    }
    node->_parent = nullptr;
    children.erase(std::find(children.begin(), children.end(), node));

    update_minSize();
}

bool GUI2Node::onEvent(const SDL_Event& event){
    if (!_active) return false;

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
            if (consumed_mouse_down){
                on_mouse_drag( event );
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
void GUI2Node::on_mouse_drag( const SDL_Event& event){}


// ==============================
//    class GUI2Panel
// =============================

GUI2Panel::GUI2Panel(GUI2Rect2f anchors, Vec2i pos, Vec2i size, uint32_t bgColor):
    GUI2Node(anchors, pos, size),
    bgColor(bgColor){}
GUI2Panel::GUI2Panel(uint32_t bgColor):
    GUI2Panel(FULL_RECT, {0, 0}, {0, 0}, bgColor){}

void GUI2Panel::draw_self(){ // TODO: don't need to redraw every frame ?
    Draw2D::drawRectangle ( (Vec2f)rect().min(), (Vec2f)rect().max(), COL2VEC(bgColor), true );
}


// ==============================
//    class GUI2Text
// ==============================

GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text ):
    GUI2Node(anchors, pos, size),
    text(text),
    textRenderer(text)
    { update_minSize(); }
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor ):
    GUI2Node(anchors, pos, size),
    text(text),
    textRenderer(text),
    fontSize(fontSize),
    fontColor(fontColor)
    { update_minSize(); }
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, uint32_t fontSize, uint32_t fontColor, Align align ):
    GUI2Node(anchors, pos, size),
    text(text),
    textRenderer(text),
    fontSize(fontSize),
    fontColor(fontColor),
    align(align)
    { update_minSize(); }
GUI2Text::GUI2Text( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, Align align):
    GUI2Node(anchors, pos, size),
    text(text),
    textRenderer(text),
    align(align)
    { update_minSize(); }
GUI2Text::GUI2Text( std::string text, Align align ):
    GUI2Text(FULL_RECT, {0,0}, {0,0}, text, align) {}


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
    recalculate_textSize();
    return textSize_;
}

void GUI2Text::setText( std::string text_ ){
    text = text_;
    update_minSize();
    textRenderer.set(text);
}
void GUI2Text::setFontSize( uint32_t fontSize_ ){
    fontSize = fontSize_;
    update_minSize();
}
void GUI2Text::setAlign( Align align_ ){
    align = align_;
    recalculate_textPos();
}

void GUI2Text::on_rect_updated(){
    recalculate_textPos();
}


void GUI2Text::draw_self(){
    textRenderer.draw2D({textPos_.x, textPos_.y, 0}, fontSize, COL2VEC(fontColor));
}

// ==============================
//    class GUI2Vlist
// ==============================

GUI2Vlist::GUI2Vlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, unsigned int min_sepperation, Align align ):
    GUI2Node(anchors, pos, size),
    min_sepperation(min_sepperation),
    align(align){}
GUI2Vlist::GUI2Vlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, Align align ):
    GUI2Vlist(anchors, pos, size, 5, align){}
GUI2Vlist::GUI2Vlist( unsigned int sepparation, Align align ):
    GUI2Vlist(FULL_RECT, {0,0}, {0,0}, sepparation, align){}


void GUI2Vlist::set_sepperation( unsigned int min_sepperation_ ){
    min_sepperation = min_sepperation_;
    update_minSize();
    update_children_rects();
}
unsigned int GUI2Vlist::get_sepperation(){
    return min_sepperation;
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
        minSize.y += childMinSize.y + min_sepperation;
    }
    minSize.y -= min_sepperation;
    return minSize;
}

void GUI2Vlist::update_child_rect(GUI2Node* child){
    update_children_rects();
}
void GUI2Vlist::update_children_rects(){
    unsigned int ymax = rect().ymax;
    unsigned int extra_sepperation = rect().size().y - minSize().y;

    for(unsigned int i=0; i<get_children().size(); i++){
        GUI2Node* child = get_children()[i];
        GUI2Rect2i R;
        unsigned int extra_sep_move = 0;
        if (get_children().size() > 1) extra_sep_move = (extra_sepperation * i) / (get_children().size()-1);

        R.ymax = ymax - extra_sep_move;
        R.ymin = R.ymax - child->minSize().y;

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
        ymax -= child->minSize().y + min_sepperation;
    }
}


// ==============================
//    class GUI2Hlist
// ==============================

GUI2Hlist::GUI2Hlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, unsigned int min_sepperation, Align align ):
    GUI2Node(anchors, pos, size),
    min_sepperation(min_sepperation),
    align(align){}
GUI2Hlist::GUI2Hlist( GUI2Rect2f anchors, Vec2i pos, Vec2i size, Align align ):
    GUI2Hlist(anchors, pos, size, 5, align){}
GUI2Hlist::GUI2Hlist( unsigned int min_sepperation, Align align ):
    GUI2Hlist(FULL_RECT, {0,0}, {0,0}, min_sepperation, align){}

void GUI2Hlist::set_sepperation( unsigned int min_sepperation_ ){
    min_sepperation = min_sepperation_;
    update_minSize();
    update_children_rects();
}
unsigned int GUI2Hlist::get_sepperation(){
    return min_sepperation;
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
        minSize.x += childMinSize.x + min_sepperation;
    }
    minSize.x -= min_sepperation;
    return minSize;
}

void GUI2Hlist::update_child_rect(GUI2Node* child){
    update_children_rects();
}
void GUI2Hlist::update_children_rects(){
    unsigned int xmin = rect().xmin;
    unsigned int extra_sepperation = rect().size().x - minSize().x;

    for(unsigned int i=0; i<get_children().size(); i++){
        GUI2Node* child = get_children()[i];
        GUI2Rect2i R;
        unsigned int extra_sep_move = 0;
        if (get_children().size() > 1) extra_sep_move = (extra_sepperation*i)/(get_children().size()-1);

        R.xmin = xmin + extra_sep_move;
        R.xmax = R.xmin + child->minSize().x;

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
        xmin += child->minSize().x + min_sepperation;
    }
}

// ==============================
//    class GUI2ButtonBase
// ==============================

GUI2ButtonBase::GUI2ButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command ):
    GUI2Node(anchors, pos, size, true),
    command(command){}
GUI2ButtonBase::GUI2ButtonBase( const std::function<void()>& command ):
    GUI2ButtonBase(FULL_RECT, {0,0}, {0,0}, command){}

void GUI2ButtonBase::on_mouse_click(){
    if (command != nullptr) command();
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
        panel = (GUI2Panel*)addChild(new GUI2Panel());
        panel->bgColor = bgColor;
    }
GUI2Button::GUI2Button( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void()>& command ):
    GUI2Button(anchors, pos, size, command, 0xA0A0A0, 0x989898, 0x808080){}
GUI2Button::GUI2Button( const std::function<void()>& command ):
    GUI2Button(FULL_RECT, {0,0}, {0,0}, command){}

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
//    class GUI2TextButton
// ==============================

GUI2TextButton::GUI2TextButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, std::string text, const std::function<void()>& command ):
    GUI2Button(anchors, pos, size, command),
    text(text)
    {
        text_node = (GUI2Text*)addChild(new GUI2Text(text, GUI2Text::Align::CENTER_LEFT));
    }
GUI2TextButton::GUI2TextButton( std::string text, const std::function<void()>& command ):
    GUI2TextButton(FULL_RECT, {0, 0}, {0, 0}, text, command){}

int GUI2TextButton::text_padding(){return _text_padding;}
void GUI2TextButton::set_text_padding(int text_padding_){
    _text_padding = text_padding_;
    text_node->set_pos({_text_padding, _text_padding});
    text_node->set_size({-_text_padding*2, -_text_padding*2});
}

// ==============================
//    class GUI2ToggleButtonBase
// ==============================

GUI2ToggleButtonBase::GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, bool* bound_bool ):
    GUI2Node(anchors, pos, size, true),
    bound_bool(bound_bool),
    command(command){}
GUI2ToggleButtonBase::GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command ):
    GUI2ToggleButtonBase(anchors, pos, size, command, nullptr){}
GUI2ToggleButtonBase::GUI2ToggleButtonBase( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool ):
    GUI2ToggleButtonBase(anchors, pos, size, nullptr, bound_bool){}

void GUI2ToggleButtonBase::on_mouse_click(){
    if (bound_bool != nullptr) active = *bound_bool;
    active = !active;
    if (bound_bool != nullptr) *bound_bool = active;
    if (command != nullptr) command(active);
}

bool GUI2ToggleButtonBase::is_active(){
    if (bound_bool != nullptr) return *bound_bool;
    return active;
}

void GUI2ToggleButtonBase::set(bool active_){
    if (bound_bool != nullptr) active = *bound_bool;
    if (active == active_) return;

    active = active_;
    if (bound_bool != nullptr) *bound_bool = active;

    if (command != nullptr) command(active);
}

void GUI2ToggleButtonBase::bind_bool(bool* bound_bool_){
    bound_bool = bound_bool_;
}

void GUI2ToggleButtonBase::set_command(const std::function<void(bool)>& command_){
    command = command_;
}

// ==============================
//    class GUI2ToggleButton
// ==============================

GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, bool* bound_bool, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed, uint32_t bgColorActive, uint32_t bgColorActiveHover, uint32_t bgColorActivePressed ):
    GUI2ToggleButtonBase(anchors, pos, size, command, bound_bool),
    bgColor(bgColor),
    bgColorHover(bgColorHover),
    bgColorPressed(bgColorPressed),
    bgColorActive(bgColorActive),
    bgColorActiveHover(bgColorActiveHover),
    bgColorActivePressed(bgColorActivePressed)
    {
        panel = (GUI2Panel*)addChild(new GUI2Panel());
        panel->bgColor = bgColor;
    }
GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command, bool* bound_bool ):
    GUI2ToggleButton(anchors, pos, size, command, bound_bool, 0xA0A0A0, 0x989898, 0x808080, 0x00FF00, 0x00E800, 0x00D000){}
GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, const std::function<void(bool)>& command ):
    GUI2ToggleButton(anchors, pos, size, command, nullptr){}

GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool ):
    GUI2ToggleButton(anchors, pos, size, nullptr, bound_bool){}
GUI2ToggleButton::GUI2ToggleButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool, uint32_t bgColor, uint32_t bgColorHover, uint32_t bgColorPressed, uint32_t bgColorActive, uint32_t bgColorActiveHover, uint32_t bgColorActivePressed ):
    GUI2ToggleButton(anchors, pos, size, nullptr, bound_bool, bgColor, bgColorHover, bgColorPressed, bgColorActive, bgColorActiveHover, bgColorActivePressed){}
GUI2ToggleButton::GUI2ToggleButton( const std::function<void(bool)>& command ):
    GUI2ToggleButton(FULL_RECT, {0,0}, {0,0}, command){}

void GUI2ToggleButton::on_mouse_enter(){
    hovering = true;
}
void GUI2ToggleButton::on_mouse_exit(){
    hovering = false;
}
void GUI2ToggleButton::on_mouse_down(){
    pressed = true;
}
void GUI2ToggleButton::on_mouse_up(){
    pressed = false;
}

void GUI2ToggleButton::draw_self(){
    if (pressed){
        panel->bgColor = is_active() ? bgColorActivePressed : bgColorPressed;
    } else if (hovering){
        panel->bgColor = is_active() ? bgColorActiveHover : bgColorHover;
    } else {
        panel->bgColor = is_active() ? bgColorActive : bgColor;
    }
}

// ==============================
//    class GUI2TextToggleButton
// ==============================

GUI2ToggleTextButton::GUI2ToggleTextButton( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool* bound_bool, std::string text, std::function<void(bool)> command ):
    GUI2ToggleButton(anchors, pos, size, command, bound_bool)
    {
        text_node = (GUI2Text*)addChild(new GUI2Text(text, GUI2Text::Align::CENTER_LEFT));
    }
GUI2ToggleTextButton::GUI2ToggleTextButton( bool* bound_bool, std::string text ):
    GUI2ToggleTextButton(FULL_RECT, {0, 0}, {0, 0}, bound_bool, text, nullptr){}
GUI2ToggleTextButton::GUI2ToggleTextButton( std::string text, std::function<void(bool)> command ):
    GUI2ToggleTextButton(FULL_RECT, {0, 0}, {0, 0}, nullptr, text, command){}

int GUI2ToggleTextButton::text_padding(){return _text_padding;}
void GUI2ToggleTextButton::set_text_padding(int padding){
    _text_padding = padding;
    text_node->set_pos({_text_padding, _text_padding});
    text_node->set_size({-_text_padding*2, -_text_padding*2});
}

// ==============================
//    class GUI2Dragable
// ==============================

GUI2Dragable::GUI2Dragable( GUI2Rect2f anchors, Vec2i pos, Vec2i size ):
    GUI2Node(anchors, pos, size, true){}

void GUI2Dragable::on_mouse_drag( const SDL_Event& event ){
    Vec2i delta = {event.motion.xrel, event.motion.yrel};
    set_pos(pos() + delta);
}

// ==============================
//    class GUI2SliderT
// ==============================

template <class T>
GUI2SliderT<T>::GUI2SliderT( GUI2Rect2f anchors, Vec2i pos, Vec2i size, T min, T max, T* bound_value, const std::function<void(T)>& command, bool interactive ):
    GUI2Node(anchors, pos, size, true),
    min(min),
    max(max),
    value(min),
    bound_value(bound_value),
    interactive(interactive),
    command(command)
    {
        bgPanel = (GUI2Panel*)addChild(new GUI2Panel());
        fillPanel = (GUI2Panel*)addChild(new GUI2Panel(0x00FF00));
    }

template <class T>
void GUI2SliderT<T>::bind_value( T* bound_value_ ){
    bound_value = bound_value_;
    if (value == *bound_value) return;
    value = *bound_value;
    if (command != nullptr) command(value);
}

template <class T>
void GUI2SliderT<T>::set_min(T min_) {
    min = min_;
};
template <class T>
void GUI2SliderT<T>::set_max(T max_) {
    max = max_;
}
template <class T>
void GUI2SliderT<T>::set_value(T value_) {
    if (bound_value != nullptr) value = *bound_value;
    if (value == value_) return;

    value = value_;
    *bound_value = value;
    if (command != nullptr) command(value);
}

template <class T>
void GUI2SliderT<T>::on_mouse_drag(const SDL_Event& event ){
    float delta_proggress = 0.5 * (float)event.motion.xrel / rect().size().x;
    T new_value = (delta_proggress * (max - min)) + value;
    if (new_value < min) new_value = min;
    if (new_value > max) new_value = max;
    set_value(new_value);
}


template <class T>
void GUI2SliderT<T>::draw_self(){
    if (bound_value != nullptr && value != *bound_value){
        value = *bound_value;
        if (command != nullptr) command(value);
    }

    float v = value;
    if (v < min) v = min;
    if (v > max) v = max;

    float progress = (float)(v - min) / (float)(max - min);
    fillPanel->set_anchors({0, 0, progress, 1});
}

template class GUI2SliderT<int>;
template class GUI2SliderT<float>;


// ==============================
//    class GUI2TextSliderT
// ==============================

template <class T>
GUI2TextSliderT<T>::GUI2TextSliderT( GUI2Rect2f anchors, Vec2i pos, Vec2i size, T min, T max, T* bound_value, std::string format, const std::function<void(T)>& command, bool interactive ):
    GUI2SliderT<T>(anchors, pos, size, min, max, bound_value, [&](T value){on_value_update(value);}, interactive),
    command(command),
    format(format)
    {
        text = (GUI2Text*)GUI2Node::addChild(new GUI2Text("NULL", GUI2Text::Align::CENTER_LEFT));
    }
template <class T>
GUI2TextSliderT<T>::GUI2TextSliderT( T min, T max, T* bound_value, std::string format, const std::function<void(T)>& command, bool interactive ):
    GUI2SliderT<T>(FULL_RECT, {0, 0}, {0, 0}, min, max, bound_value, [&](T value){on_value_update(value);}, interactive){}

template <class T>
void GUI2TextSliderT<T>::on_value_update( T value ){
    std::string str = format;
    std::string replace = std::string("$value");

    size_t start_pos = str.find(replace);
    if(start_pos != std::string::npos){
        if (start_pos+replace.length() < format.length() && format[start_pos+replace.length()] == '%'){
            // fancier formatting
            size_t end_pos = str.find('%', start_pos+replace.length()+1);
            std::string format_specifier;
            if (end_pos != std::string::npos){
                format_specifier = str.substr(start_pos + replace.length(), end_pos - start_pos - replace.length());
            }else{
                format_specifier = format.substr(start_pos + replace.length());
            }
            char formated[128];
            std::snprintf(formated, sizeof(formated), format_specifier.c_str(), value);
            str.replace(start_pos, end_pos - start_pos, formated);
        }
        else{
            str.replace(start_pos, replace.length(), std::to_string(value));
        }
    }

    text->setText(str);

    if (command != nullptr) command(value);
}
template class GUI2TextSliderT<int>;
template class GUI2TextSliderT<float>;

// ==============================
//    class GUI2FloatingMenu
// ==============================

GUI2FloatingMenu::MenuButton::MenuButton( std::string text, std::string keybind_text, const std::function<void()> command ):
    GUI2Button(command)
    {
        GUI2Hlist* hlist = (GUI2Hlist*)addChild(new GUI2Hlist());
        hlist->addChild(new GUI2Text(text, GUI2Text::Align::CENTER_LEFT));
        hlist->addChild(new GUI2Text("   ", GUI2Text::Align::CENTER)); // spacing TODO: add GUI2Padding
        ((GUI2Text*)hlist->addChild(new GUI2Text(keybind_text, GUI2Text::Align::CENTER_RIGHT)))->fontColor = 0x404040;
    }
GUI2FloatingMenu::MenuButton::MenuButton( std::string text, std::string keybind_text, int submenu_idx, GUI2FloatingMenu* menu ):
    MenuButton(text, keybind_text, nullptr)
    {
        this->submenu_idx = submenu_idx;
        this->menu = menu;
    }

GUI2FloatingMenu::SubmenuButton::SubmenuButton( std::string text, int submenu_idx, GUI2FloatingMenu* menu ):
    submenu_idx(submenu_idx),
    menu(menu),
    GUI2ToggleButton([menu,submenu_idx](bool active){if (active) menu->open_submenu(submenu_idx); else menu->close_submenu(); })
    {
        GUI2Hlist* hlist = (GUI2Hlist*)addChild(new GUI2Hlist());
        hlist->addChild(new GUI2Text(text, GUI2Text::Align::CENTER_LEFT));
        hlist->addChild(new GUI2Text("   ", GUI2Text::Align::CENTER)); // spacing TODO: add GUI2Padding
        ((GUI2Text*)hlist->addChild(new GUI2Text(">", GUI2Text::Align::CENTER_RIGHT)))->fontColor = 0x101010;
    }

void GUI2FloatingMenu::SubmenuButton::on_mouse_enter(){
    GUI2ToggleButton::on_mouse_enter();

    if (menu->current_submenu != submenu_idx && menu->current_submenu != -1){
        menu->open_submenu(submenu_idx);
    }
}

GUI2FloatingMenu::GUI2FloatingMenu( GUI2Rect2f anchors, Vec2i pos, Vec2i size, bool collapse_self ):
    GUI2Node(anchors, pos, size, true),
    collapse_self(collapse_self)
    {
        bgPanel = (GUI2Panel*)addChild(new GUI2Panel());
        vlist = (GUI2Vlist*)addChild(new GUI2Vlist(1));
    }
GUI2FloatingMenu::GUI2FloatingMenu( bool collapse_self ):
    GUI2FloatingMenu({0, 0, 0, 0}, {0, 0}, {0, 0}, collapse_self){}


void GUI2FloatingMenu::open_submenu(int idx){
    if (idx < 0 || idx >= submenus.size()){
        printf("GUI2FloatingMenu::open_submenu(%i): index out of bounds\n", idx);
        return;
    }
    if (current_submenu == idx) return;
    close_submenu();
    current_submenu = idx;

    submenus[idx]->set_active(true);
    submenuButtons[idx]->set(true);
}

void GUI2FloatingMenu::close_submenu(){
    if (current_submenu == -1) return;
    submenus[current_submenu]->close_submenu();
    submenus[current_submenu]->set_active(false);
    submenuButtons[current_submenu]->set(false);
    current_submenu = -1;
}

void GUI2FloatingMenu::update_child_rect(GUI2Node* node){
    update_children_rects();
}
void GUI2FloatingMenu::update_children_rects(){
    if (bgPanel == nullptr) return;
    if (vlist == nullptr) return;
    GUI2Node::update_child_rect(bgPanel);
    GUI2Node::update_child_rect(vlist);

    for (int i=0; i<submenus.size(); i++){
        GUI2Node::update_child_rect(submenus[i]); // calculates the rect size
        GUI2Rect2i rect;
        rect.xmin = submenuButtons[i]->rect().xmax;
        rect.ymax = submenuButtons[i]->rect().ymax;

        rect.xmax = rect.xmin + submenus[i]->rect().size().x;
        rect.ymin = rect.ymax - submenus[i]->rect().size().y;

        submenus[i]->set_rect(rect);
    }
}

GUI2Node* GUI2FloatingMenu::add_item(GUI2Node* item){
    return vlist->addChild(item);
}

GUI2Button* GUI2FloatingMenu::add_button(std::string text, std::string keybind_text, const std::function<void()> command){
    return (GUI2Button*)add_item(new MenuButton(text, keybind_text, command));
}

GUI2FloatingMenu* GUI2FloatingMenu::add_submenu(std::string name){
    int idx = submenus.size();

    SubmenuButton* button = new SubmenuButton(name, idx, this);
    submenuButtons.push_back(button);
    add_item(button);

    GUI2FloatingMenu* submenu = new GUI2FloatingMenu(false);
    submenus.push_back(submenu);
    submenu->set_active(false);
    addChild(submenu);

    return submenu;
}

bool GUI2FloatingMenu::onEvent(const SDL_Event& event){
    if (!active()) return false;

    for (auto c : get_children()){
        if (c->onEvent(event)) return true;
    }

    if (!collapse_self) return false;

    if (event.type != SDL_MOUSEBUTTONUP) return true;
    if (current_submenu != -1) close_submenu();
    set_active(false);

    return true;
}

// ==============================
//    class GUI2Toolbar
// ==============================

GUI2Toolbar::CustomTextButton::CustomTextButton( std::string text, const std::function<void(bool)> command, int tabIdx, GUI2Toolbar* toolbar ):
    tabIdx(tabIdx),
    toolbar(toolbar),
    GUI2ToggleTextButton(text, command){}

void GUI2Toolbar::CustomTextButton::on_mouse_enter(){
    GUI2ToggleTextButton::on_mouse_enter();

    if (toolbar->openTabIdx != tabIdx && toolbar->openTabIdx != -1){
        toolbar->open_tab(tabIdx);
    }
}

GUI2Toolbar::GUI2Toolbar( ):
    GUI2Node( GUI2Rect2f(0,1,1,1), (Vec2i){0,0}, (Vec2i){0, 0} ),
    bgPanel(nullptr),
    buttonsHlist(nullptr)
    {
        set_expandMode(DOWN_RIGHT);
        bgPanel = (GUI2Panel*)addChild(new GUI2Panel());
        buttonsHlist = (GUI2Hlist*)addChild(new GUI2Hlist({0, 0, 0, 1}, {0, 0}, {0, 0}));
    }

void GUI2Toolbar::close_tab(int idx){
    if (idx < 0 || idx >= tabs.size()) {
        printf("GUI2Toolbar::close_tab(%i): index out of bounds, ignoring...\n", idx);
        return;
    }
    openTabIdx = -1;

    tabs[idx]->close_submenu();
    tabs[idx]->set_active(false);
    buttons[idx]->set(false);
}
void GUI2Toolbar::open_tab(int idx){
    if (idx < 0 || idx >= tabs.size()) {
        printf("GUI2Toolbar::open_tab(%i): index out of bounds, ignoring...\n", idx);
        return;
    }
    if (openTabIdx == idx) return;
    if (openTabIdx != -1) close_tab(openTabIdx);
    openTabIdx = idx;

    tabs[idx]->set_active(true);
    buttons[idx]->set(true);

    update_children_rects();
}

GUI2FloatingMenu* GUI2Toolbar::addTab( std::string name ){
    uint idx = tabs.size();

    CustomTextButton* button = (CustomTextButton*)buttonsHlist->addChild(new CustomTextButton(name, [this,idx](bool active){if (active) this->open_tab(idx); else this->close_tab(idx); }, idx, this));
    button->set_text_padding(3);
    buttons.push_back(button);

    GUI2FloatingMenu* tab = (GUI2FloatingMenu*)addChild(new GUI2FloatingMenu(false));
    tab->set_active(false);
    tabs.push_back(tab);

    return tab;
}

bool GUI2Toolbar::onEvent( const SDL_Event& event ){
    for(auto& c:get_children()){
        if (c->onEvent(event)) return true;
    }

    if (openTabIdx == -1) return false;

    if (event.type != SDL_MOUSEBUTTONUP) return true;
    if (openTabIdx != -1) buttons[openTabIdx]->set(false);

    return true;
}

void GUI2Toolbar::update_child_rect(GUI2Node* node){
    update_children_rects();
}
void GUI2Toolbar::update_children_rects(){
    if (buttonsHlist == nullptr) return;
    GUI2Node::update_child_rect(bgPanel);
    GUI2Node::update_child_rect(buttonsHlist);

    for (int i=0; i<tabs.size(); i++){
        GUI2Node::update_child_rect(tabs[i]);
        GUI2Rect2i rect;
        rect.xmin = buttons[i]->rect().xmin;
        rect.ymax = buttons[i]->rect().ymin;

        rect.xmax = rect.xmin + tabs[i]->rect().size().x;
        rect.ymin = rect.ymax - tabs[i]->rect().size().y;

        tabs[i]->set_rect(rect);
    }
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

    glDisable   ( GL_DEPTH_TEST  );
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
