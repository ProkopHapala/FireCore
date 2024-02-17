
#ifndef  GUI_h
#define  GUI_h

#include <SDL2/SDL.h>
#include <SDL2/SDL_opengl.h>

#include "Draw.h"
#include "Draw2D.h"
#include "Draw3D.h"

#include "Table.h"

#include <string>
#include <functional>

/*
TODO:
 - draw small glyph by simple language ... like:
    goto    100.0 100.0
    line    200.0 250.0
    ngon    6 16.5
*/

class GUI_stepper{ public:
    int x0=0,x1=0;
    GUI_stepper()=default;
    GUI_stepper(int n, int l){ x0=n*fontSizeDef; x1=(n+l)*fontSizeDef; }
    inline void step(int n){ x0=x1; x1=x0+n*fontSizeDef; }
};

struct Command{
    int id;
    int key;
    std::string name;
};

class Commander{ public:
    std::vector<Command>        commands;
    std::unordered_map<int,int> keymap;
    //std::unordered_map<std::string,int> byname;

    inline void add(int id, int key, std::string name=""){
        int i=commands.size();
        commands.push_back( (Command){id,key,name} );
        auto got = keymap.find(key);
        if(got!=keymap.end()){ printf("WARRNING KeyBind: key already used %i %c\n", key, (char)key ); return; };
        keymap.insert({key,i});
    }

    inline void rebind( int i, int key ){
        auto got = keymap.find(key);
        if(got!=keymap.end()){ printf("WARRNING KeyBind: key already used %i %c\n", key, (char)key ); return; };
        int old_key  = commands[i].key;
        auto old_got = keymap.find(old_key);
        if( old_got != keymap.end() ){ keymap.erase(old_got); }
        keymap.insert({key,i});
        commands[i].key=key;
    }
};

//class GUIAbstractPanel;

//extern GUIAbstractPanel* GUI_mouse_panel;
//extern Vec2i GUI_mouse_old_pos;
extern int GUI_fontTex;

//void GUI_globalEventHandler(const SDL_Event* event );

class GUI;


class GUIAbstractPanel;

// Alterntively we may use events
//    https://wiki.libsdl.org/SDL_UserEvent
//    https://stackoverflow.com/questions/26707536/event-driven-programming-callback-vs-message-polling
class GUIEventCallback{public:
    virtual int GUIcallback(GUIAbstractPanel* caller)=0;
};

// ==============================
//    class GUITextInput
// ==============================

class GUITextInput{ public:

    int         curPos=0;
	std::string inputText;

    bool     isNumber=true,checkRange=false;
	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0;
	SDL_Keycode num_op = 0;

	bool     modified=true,entered=false;

	// ==== functions

	void         applyVal( float f );

    virtual void view3D ( const Vec3d& pos, int fontTex, float textSize );
    virtual void viewHUD( const Vec2i& pos, int fontTex, bool bBack=true );
    virtual void onKeyDown( SDL_Event e );
	virtual void onText   ( SDL_Event e );

};

// ==============================
//    class GUIAbstractPanel
// ==============================

class GUIAbstractPanel{ public:
    //int textSz = fontSizeDef;
	int  xmin=256,xmax=128,ymin=0,ymax=0;
	bool visible=true, disabled=false;

	uint32_t bgColor=0xA0A0A0, textColor=0x000000;

	bool     redraw=true;
	int      gllist=0;

	//int      fontTex=0;
    //char*    caption=NULL;
    std::string caption;
	//void (*command)(double,void*) = NULL;
    //std::function<void(double)> command =0;
    std::function<void(GUIAbstractPanel*)> command =0; // = [] { return 0; }

	// ==== functions


    void initPanel( const std::string& caption_, int xmin_, int ymin_, int xmax_, int ymax_ );

    GUIAbstractPanel(){};
    GUIAbstractPanel( const std::string& caption, int xmin, int ymin, int xmax, int ymax ){ initPanel(caption, xmin,ymin,xmax,ymax); };
    virtual ~GUIAbstractPanel()=default;

    virtual void open(){};
    virtual void close(){};
	virtual void moveTo(int x, int y);
	virtual void moveBy(int dx,int dy);

	virtual void              onKeyDown(             const SDL_Event& e, GUI& gui );
	virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& e, GUI& gui );
    virtual void              onText(                const SDL_Event& e, GUI& gui );

    virtual void view  ( );
    //virtual void tryRender();
    //void view     ();
    void tryRender();
    virtual void render();

    // inline fnctions

    inline void draw      ( ){ 
        //printf( "GUIAbstractPanel::draw() \n" );
        tryRender(); 
        view(); 
        redraw=false; 
        //printf( "GUIAbstractPanel::draw() END\n" );
    };
    inline bool check     ( int  x, int  y ){
        //printf( "check x %i <%i...%i>   y %i <%i...%i>\n", x, xmin, xmax,   y, ymin, ymax );
        return (x>xmin)&&(x<xmax)&&(y>ymin)&&(y<ymax);
    }
	inline void toRelative ( int& x, int& y ){ x-=xmin; y-=ymin; }
    inline GUIAbstractPanel* setCommand( const std::function<void(GUIAbstractPanel* panel)>& command_ ){ command=command_; return this; }
};

// ==============================
//       class GUIPanel
// ==============================

class GUIPanel : public GUIAbstractPanel { public:
	bool isSlider=true, isButton=false;

	uint32_t barColor=0x00FF00;
    bool viewVal=true;
    bool valBelow=true;

	bool     executed=false;
	int      curPos=0;
	std::string inputText;

	float    vmin=0.0f, vmax=1.0f;
	double   value=0.0;
	bool     isInt = false;
	double* master=0;

    // ==== functions

    GUIPanel()=default;
    GUIPanel( const std::string& caption, int xmin, int ymin, int xmax, int ymax, bool isSlider_=true, bool isButton_=true, bool isInt_=false, bool viewVal_=true ){
        initPanel(caption, xmin,ymin,xmax,ymax); isSlider=isSlider_; isButton=isButton_; isInt=isInt_; viewVal=viewVal_;
        command=0;
    };

    virtual void view() override;
	//virtual void tryRender();
	void render();
    virtual void              onKeyDown( const SDL_Event&  e, GUI& gui )                  override;
    virtual void              onText   ( const SDL_Event&  e, GUI& gui )                  override;
    virtual GUIAbstractPanel* onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ) override;

	// ===== inline functions
    inline int    getIntVal()        { return round(value); };
    inline void   val2text()         { if(isInt){ inputText = std::to_string(getIntVal()); }else{ inputText = std::to_string(value); }; };
	inline double x2val( float  x   ){ return ( x*(vmax-vmin)/(xmax-xmin) )+ vmin; };
	inline float  val2x( double val ){ return (val-vmin)*(xmax-xmin)/(vmax-vmin);  };

	inline int   x2val_int( float  x   ){ return round( x*(vmax-vmin)/(xmax-xmin) )+ vmin; };
	//inline float val2x_int( double val ){ return (val-vmin)*(xmax-xmin)/(vmax-vmin);  };

    inline void syncRead (){ if(master)value=*master; }
    inline void syncWrite(){ if(master)*master=value; }

	inline GUIPanel* setRange(float vmin_, float vmax_){ vmin=vmin_; vmax=vmax_; return this; };
    inline GUIPanel* setValue(float val_){ value=_clamp(val_,vmin,vmax); redraw=true; return this; };

};

// ==============================
//     class  ScisorBox
// ==============================

class ScisorBox : public GUIAbstractPanel { public:

    // ==== functions

    void apply();
    void initScisor( const std::string& caption, int xmin_, int ymin_, int xmax_, int ymax_ );

    ScisorBox(){};
    ScisorBox( const std::string& caption, int xmin_, int ymin_, int xmax_, int ymax_ ){ initScisor( caption, xmin_, ymin_, xmax_, ymax_ ); };

    //virtual void draw     ( );
    //virtual void tryRender( );
    virtual void render( )                                                                override;
    virtual GUIAbstractPanel* onMouse ( int x, int y, const SDL_Event&  event, GUI& gui ) override;

    //virtual void onKeyDown( SDL_Event e ){};
    //virtual void onText   ( SDL_Event e ){};

};

// ==============================
//     class  MultiPanel
// ==============================

class MultiPanel : public GUIAbstractPanel { public:
    int nsubs;
    GUIPanel ** subs;
    bool opened = true;
    int dy;

    // ==== functions

    void initMulti( const std::string& caption, int xmin_, int ymin_, int xmax_, int dy, int nsubs_, bool isSlider=true, bool isButton=true, bool isInt=false, bool viewVal=true );
    MultiPanel(){};
    MultiPanel(const std::string& caption, int xmin, int ymin, int xmax, int dy, int nsubs, bool isSlider=true, bool isButton=true, bool isInt=false, bool viewVal=true ){ if(dy==0){dy=2*fontSizeDef;} initMulti( caption, xmin, ymin, xmax, dy, nsubs, isSlider,isButton,isInt,viewVal); }

    virtual void open()override;
    virtual void close()override;
    virtual void moveBy(int dx, int dy) override;
    virtual void view  ( )override;
    //virtual void tryRender( );
    virtual void render( )override;
    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui )override;

    void toggleOpen();

    //virtual void onKeyDown( SDL_Event e, GUI& gui ){};
    //virtual void onText   ( SDL_Event e, GUI& gui ){};

};

// ==============================
//     class  CheckBoxList
// ==============================

#define _addBox(var)   _.addBox( #var, &var );

struct CheckBox{
    bool  val   =0; //
    bool* master=0; // if not zero, this is the true value
    std::string label;
    inline void read (){ if(master)val=*master; }
    inline void write(){ if(master)*master=val; }
    inline void flip (){ printf("flip  %i %i ", val,*master ); val=!val; write(); printf("-> %i %i \n", val,*master); }
};

class CheckBoxList : public GUIAbstractPanel { public:
    std::vector<CheckBox> boxes;
    int dy;
    uint32_t checkColor=0x00FF00;


    // ==== functions

    void addBox(std::string label, bool* ptr){
        boxes.push_back((CheckBox){false,ptr,label});
    }

    void initCheckBoxList( int xmin_, int ymin_, int xmax_, int dy=fontSizeDef*2 );

    CheckBoxList(){};
    CheckBoxList(int xmin, int ymin, int xmax, int dy=fontSizeDef*2){ initCheckBoxList( xmin, ymin, xmax, dy); }

    //virtual void  open();
    //virtual void close();
    //void    toggleOpen();
    //virtual void moveBy(int dx, int dy);
    void update();

    virtual void view  ( )                                                               override;
    virtual void render( )                                                               override;
    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui )  override;

    inline void syncRead (){ for(CheckBox& b: boxes){ b.read (); } }
    inline void syncWrite(){ for(CheckBox& b: boxes){ b.write(); } }
};


// ==============================
//     class  CommandList
// ==============================

class CommandList : public GUIAbstractPanel { public:
    Commander*                commands;
    //std::function<void(int)>* commandDispatch=0;
    std::function<void(int)> commandDispatch;
    //std::function commandDispatch;
    bool bDispatch=0;
    int dy;
    uint32_t modColor=0x00FFFF;
    int icmdbind = -1;
    //uint32_t checkColor=0x00FF00;

    // ==== functions

    void initCommandList( int xmin_, int ymin_, int xmax_, int dy=fontSizeDef*2 );

    CommandList(){};
    CommandList(int xmin, int ymin, int xmax, int dy=fontSizeDef*2){ initCommandList( xmin, ymin, xmax, dy); }

    bool     getKeyb(int key);
    void         update( );
    virtual void view  ( )                                                              override;
    virtual void render( )                                                              override;
    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui ) override;

};


// ==============================
//     class  KeyTree
// ==============================

/*

#### For dispatching tree of commands
 * BackSpace go one level up
 * Escape restat tree from root

should look like this:
---------------------
* [F]ile [E]dit [V]iew
    * [O]pen [S]ave [L]oad [Q]uit

*/



// ==============================
//     class  DropDownList
// ==============================

//static const char* exampleDropDownListItems[3] = {"Item1","Item2","Item3"};

class DropDownList : public GUIAbstractPanel { public:
    bool bOpened = false;
    int nSlots=5;
    //int iSlot0=0;

    int iSelected=0;
    int iItem0 = 0;
    //int nItems = 0;
    //char** labels = (char**)exampleDropDownListItems;
    std::vector<std::string> labels;
    //GUIEventCallback* onSelect = 0;

    //int nsubs;
    //GUIPanel ** subs;

    // ==== functions

    DropDownList* addItem(const std::string& label);
    void initList( const std::string& caption, int xmin_, int ymin_, int xmax_, int nSlots_ );
    int selectedToStr(char* str);

    DropDownList(){}
    DropDownList( const std::string& caption, int xmin, int ymin, int xmax, int nSlots){ initList(caption,xmin,ymin,xmax,nSlots); }

    virtual void open() override;
    virtual void close() override;

    //virtual void view ( );
    //virtual void tryRender( );
    virtual void render( )                                                                override;
    virtual GUIAbstractPanel* onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ) override;

    //virtual void onKeyDown( SDL_Event e )override{};
    //virtual void onText   ( SDL_Event e )override{};

};



// ==============================
//     class  TreeView
// ==============================

#include "Tree.h"

class TreeViewItem{ public:
    //bool open=false;
    bool open=true;
    int level=0;
    int nth_line=0;
    std::string caption;

    TreeViewItem(){};
    TreeViewItem( std::string caption_){ caption=caption_; };
};
typedef PTree<TreeViewItem> TreeViewTree;

//static const char* exampleDropDownListItems[3] = {"Item1","Item2","Item3"};
class TreeView : public GUIAbstractPanel { public:
    //typedef Tree<TreeViewItem> TreeT;
    int iItem0=0;
    int iSelected=0;
    int nSlots = 5;
    TreeViewTree root;
    std::vector<TreeViewTree*> lines;


    virtual void view()                                                                   override;
    virtual void render()                                                                 override;
    virtual GUIAbstractPanel* onMouse  ( int x, int y, const SDL_Event& event, GUI& gui ) override;

    void updateLines( TreeViewTree& node, int level );
    inline void updateLines(){ lines.clear(); updateLines(root,0); };

    void initTreeView( const std::string& caption, int xmin_, int ymin_, int xmax_, int nSlots_ );

    TreeView(){}
    TreeView( const std::string& caption, int xmin, int ymin, int xmax, int nSlots){
        initTreeView(caption,xmin,ymin,xmax,nSlots);
    }

    //virtual void onKeyDown( SDL_Event e ){};
    //virtual void onText   ( SDL_Event e ){};
};

// ==============================
//     TableView
// ==============================

class TableView : public GUIAbstractPanel { public:
    //typedef Tree<TreeViewItem> TreeT;
    Table* table=0;
    int i0=0,j0=0;
    int imax=0,jmax=0;
    int i=0,j=0;  // cursors
    std::vector<int> nchs;
    std::vector<int> xs;

    GUITextInput* input=0;

    void initTableView( Table* table_, const std::string& caption_, int xmin_, int ymin_, int i0_, int j0_, int imax_, int jmax_ );
    virtual void render()                                                               override;
    virtual void onKeyDown( const SDL_Event& e, GUI& gui )                              override;
    virtual void onText   ( const SDL_Event& e, GUI& gui )                              override;
    virtual GUIAbstractPanel* onMouse( int x, int y, const SDL_Event& event, GUI& gui ) override;

    TableView()=default;
    TableView( Table* table_, const std::string& caption_, int xmin_, int ymin_, int i0_, int j0_, int imax_, int jmax_ ){
        initTableView( table_, caption_, xmin_, ymin_, i0_, j0_, imax_, jmax_ );
    }

    inline void view ( ){
        glCallList( gllist );
        if(input)input->viewHUD( {xmin+xs[j+1]-xs[j0],xmax+((i0-i)*2*fontSizeDef)}, fontSizeDef, true );
    }

};


// ==============================
//    class GUI
// ==============================

class GUI{ public:
    bool bKeyEvents = true;
    GUIAbstractPanel* focused = 0;
    GUIAbstractPanel* dragged = 0;
    std::vector<GUIAbstractPanel*> panels;

    GUIAbstractPanel* addPanel( GUIAbstractPanel* panel );
    GUIAbstractPanel* onEvent( int mouseX, int mouseY, const SDL_Event& event );
    void draw();

    //void layoutRow( int xmin, int xmax );
    void layoutRow( int xmin, int ymin, int xspace=0 );

    int clear(int i0=0){ int n=0; for(int i=i0; i<panels.size(); i++){  if(panels[i]){ n++; delete panels[i];} } panels.resize(i0); return n; }

    inline ~GUI(){ for(GUIAbstractPanel* panel: panels){ delete panel; } }

};

#endif
