#ifndef MolGUIapp_multi_Lua_h
#define MolGUIapp_multi_Lua_h

//#include "MolGUI.h"
#include <unordered_map>
#include <functional>

#include "LuaHelpers.h"
#include "Atoms.h"
#include "MolWorld_sp3_multi.h"
#include "MolGUIapp_Lua.h"
//lua_State  * theLua=0;

int l_getReplicaCount(lua_State *L){
    lua_pushinteger(L, W->nSystems );
    return 1; // number of return values to Lua environment
}

int l_shiftReplica(lua_State *L){
    //int   ia = Lua::getInt(L,1);
    int   isys = Lua::getInt (L,1);
    Vec3d d    = Lua::getVec3(L,2);
    //MolWorld_sp3_multi& W = *(app->W);
    //Vec3d* ps = W->ffls[isys].apos;
    //for(int i=0; i<; i++)
    printf( "l_shiftReplica(isys=%i,dpos=(%g,%g,%g))\n", isys, d.x,d.y,d.z );
    W->ffls[isys].shift( d );
    //Lua::pushVec3(L, pos);
    return 0; // number of return values to Lua environment
}

int l_systemsToGPU( lua_State *L ){
    W->pack_systems();
    W->upload();
    return 0; 
}

lua_State* initMyLua_multi( lua_State* L =0 ){
    printf( "initMyLua_multi()\n" );
    //theLua         = luaL_newstate();
    //lua_State  * L = theLua;
    if(L==0){ 
        L =  luaL_newstate(); 
        luaL_openlibs(L);
    }
    lua_register(L, "getReplicaCount",  l_getReplicaCount  );
    lua_register(L, "shiftReplica",     l_shiftReplica     );
    lua_register(L, "systemsToGPU",     l_systemsToGPU     );
    printf( "initMyLua_multi() DONE\n" );
    //return 1;
    return L;
}

#endif // MolGUIapp_Lua



