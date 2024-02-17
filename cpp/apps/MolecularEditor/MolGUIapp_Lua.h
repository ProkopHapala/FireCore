#ifndef MolGUIapp_Lua_h
#define MolGUIapp_Lua_h

//#include "MolGUI.h"
#include <unordered_map>
#include <functional>

//using Command     = std::function<void(double val)>; 
//using CommandDict = std::unordered_map<std::string,>;


#include "LuaHelpers.h"
lua_State  * theLua=0;

int l_addGUIpanel(lua_State *L){
    const char* label     = Lua::getString(L,1);
    Vec3i xyw             = Lua::getVec3i (L,2);
    Vec3i SliderButtonInt = Lua::getVec3i (L,3);
    Vec3d MinMaxCur       = Lua::getVec3  (L,4);
    const char* command   = Lua::getString(L,5);
    int ipanel = app->gui.panels.size();
    GUIPanel* panel = (GUIPanel*)app->gui.addPanel( new GUIPanel( label, xyw.x,xyw.y, xyw.x+xyw.z*fontSizeDef,xyw.y+fontSizeDef*2, SliderButtonInt.x>0, SliderButtonInt.y>0, SliderButtonInt.z>0 ) );
    //GUIPanel* panel = ((GUIPanel*)app->gui.addPanel( new GUIPanel( "Mol. Orb.", 5,ylay.x0,5+100,ylay.x1, true, true, true ) ) );
    panel->setRange(MinMaxCur.x,MinMaxCur.y);
    panel->setValue(MinMaxCur.z);
    //panel->setCommand( app->actions.get( command ) );
    app->setPanelAction( ipanel, command  );
    lua_pushinteger(L, ipanel);
    return 1;
}

int l_clearGUI   (lua_State *L){
    int n = Lua::getInt(L,1);
    int ret = app->gui.clear( n );
    lua_pushinteger(L, ret );
    return 1;
}

int l_fixAtom(lua_State *L){
    // LuaCall: fixAtom( ia, true )
    int ia = Lua::getInt(L,1);
    int b  = Lua::getInt(L,2);
    printf( "l_fixAtom(ia=%i,b=%i)\n", ia, b ); 
    double Kfix = 10.0; //app->W->ffl.Kfix;
    Quat4d constr;
    constr.f=app->W->ffl.apos[ia];
    constr.e=Kfix; 
    app->W->ffl.constr[ia]=constr; 
    //if     (b>0){ app->W->atomFixed[ia] = true;  }
    //else if(b<0){ app->W->atomFixed[ia] = false; }
    return 0; // number of return values to Lua environment
}

int l_getAtomCount(lua_State *L){
    //Lua::pushInt(L, app->W->ffl.natoms );
    lua_pushinteger(L, app->W->ffl.natoms );
    return 1; // number of return values to Lua environment
}

int l_toggleStop(lua_State *L){
    app->bRunRelax = !app->bRunRelax;
    printf( "l_toggleStop %i\n", app->bRunRelax );
    return 0; // number of return values to Lua environment
}

int l_getAtomPos(lua_State *L){
    int ia = Lua::getInt(L,1);
    Vec3d pos = app->W->ffl.apos[ia];
    Lua::pushVec3(L, pos);
    return 1; // number of return values to Lua environment
}

int l_insertQuickCommand(lua_State *L){
    const char* s = Lua::getString(L,1);
    printf( "l_insertQuickCommand `%s`\n", s );
    app->console.quick_tab.table.push_back( s );
    app->console.quick_tab.sort();
    for(int i=0; i<app->console.quick_tab.table.size(); i++){
        printf( "l_insertQuickCommand[%i] `%s`\n", i, app->console.quick_tab.table[i].c_str() );
    }
    return 1; // number of return values to Lua environment
}

int l_clearMolecules(lua_State *L){
    printf( "l_clearMolecule()\n" );
    //bool bParams = Lua::getInt(L,1);
    //bool bSurf   = Lua::getInt(L,2);
    bool bParams = false;
    bool bSurf   = false;
    printf( "l_clearMolecule(bParams=%i,bSurf=%i)\n", bParams, bSurf );
    app->bRunRelax = false;
    app->W->clear( bParams, bSurf );
    app->W->builder.printSizes();
    app->unBindMolecule();
    app->bViewBuilder = true;
    return 0; // number of return values to Lua environment
}

int l_addMoleculeFile(lua_State *L){
    const char* fname = Lua::getString(L,1);
    Vec3d       pos   = Lua::getVec3  (L,2);
    printf( "l_addMoleculeFile(`%s`,pos(%g,%g,%g))\n", fname, pos.x, pos.y, pos.z );
    //app->W->loadMolecule( s );
    //app->W->buildMolecule_xyz( xyz_name );
    app->bRunRelax = false;
    MolWorld_sp3* W = app->W;
    int ifrag = W->insertMolecule( fname, 0, pos );
    // W->makeMoleculeTopology();
    // W->assingMoleculeTopoTypes();
    // W->makeFFs();
    app->bViewBuilder = true;
    lua_pushinteger(L, ifrag);
    return 1; // number of return values to Lua environment
}

int l_autoCharges(lua_State *L){
    printf( "l_autoCharges()\n" );
    app->W->autoCharges();
    return 0; // number of return values to Lua environment
};

int l_autoFrags(lua_State *L){
    printf( "l_autoFrags()\n" );
    app->W->fragmentsByBonds();
    return 0; // number of return values to Lua environment
};

int l_makeFF(lua_State *L){
    printf( "l_makeFF()\n" );
    MolWorld_sp3* W = app->W;
    printf(" ======== l_makeFF() =========\n" );
    W->builder.printSizes();
    W->builder.printAtoms();
    printf(" ======== l_makeFF() 1 W->makeMoleculeTopology()\n" );
    W->builder.randomFragmentCollors();
    W->makeMoleculeTopology();
    //W->builder.printConfs();
    W->builder.printAtomConfs();
    printf(" ======== l_makeFF() 2 W->assingMoleculeTopoTypes()\n" );
    W->assingMoleculeTopoTypes();
    printf(" ======== l_makeFF() 3 W->assingMoleculeTopoTypes()\n" );
    W->makeFFs();
    app->bindMolecule( W->ffl.natoms, W->ffl.nnode, W->ff.nbonds, W->nbmol.atypes, W->nbmol.apos, W->nbmol.fapos, W->nbmol.REQs, W->ffl.pipos, W->ffl.fpipos, W->ff.bond2atom, W->ff.pbcShifts );
    app->neighs    = W->ffl.neighs;
    app->neighCell = W->ffl.neighCell;
    app->bViewBuilder = false;
    app->updateGUI();
    printf(" ======== l_makeFF() 4 DONE\n" );
    return 0; // number of return values to Lua environment
}

int initMyLua(){
    printf( "initMyLua()\n" );
    theLua         = luaL_newstate();
    lua_State  * L = theLua;
    luaL_openlibs(L);
    lua_register(L, "fix",     l_fixAtom  );
    lua_register(L, "natom",   l_getAtomCount  );
    lua_register(L, "apos",    l_getAtomPos  );
    lua_register(L, "run",     l_toggleStop  );
    lua_register(L, "command", l_insertQuickCommand  );
    lua_register(L, "clear",   l_clearMolecules  );
    lua_register(L, "add",     l_addMoleculeFile  );
    lua_register(L, "make",    l_makeFF  );
    lua_register(L, "autoCharges", l_autoCharges  );
    lua_register(L, "frags",   l_autoFrags  );

    lua_register(L, "button",    l_addGUIpanel );
    lua_register(L, "clear_gui", l_clearGUI    );
    printf( "initMyLua() DONE\n" );
    return 1;
}

#endif // MolGUIapp_Lua



