#ifndef LuaHelpers_h
#define LuaHelpers_h

// http://www.lua.org/manual/5.2/manual.html#lua_load

#include <cstring>
#include <cstdio>
//#include <string>
#include <vector>
//#include <iostream>

#include "Vec2.h"
#include "Vec3.h"
#include "Mat3.h"
#include "quaternion.h"
#include "IO_utils.h"

#ifdef __cplusplus
# include <lua5.2/lua.hpp>
#else
# include <lua5.2/lua.h>
# include <lua5.2/lualib.h>
# include <lua5.2/lauxlib.h>
#endif

namespace Lua{

    void print_error(lua_State* L) { puts( lua_tostring(L, -1)); lua_pop(L, 1); }
    bool checkError( lua_State* L, int ret ){ if(ret!=LUA_OK) print_error(L); return ret==LUA_OK; }
    //bool checkError( lua_State* L, int ret ){ if(ret!=LUA_OK) print_error(L); return ret==LUA_OK; }
    //if (luaL_dofile(L,ss[0]) != LUA_OK){ fprintf(stderr, "Error: %s\n", lua_tostring(L, -1)); lua_pop(L, 1); }

    void getError( int i,  const char * s ){ printf( "LuaERROR @Stack[%i]: %s\n", i, s ); }
    void clean   (lua_State* L) { int n = lua_gettop(L); lua_pop(L, n); }

    void dumpStack(lua_State* L ){
        int n = lua_gettop( L);  printf("LuaStack[%i] = {",n);
        for (int i=1; i<=n; i++){
            int t = lua_type( L, i);
            switch (t) {
                case LUA_TSTRING:  printf(" %i:'%s';", i, lua_tostring ( L, i)); break;
                case LUA_TBOOLEAN: printf(" %i:%s;",   i, lua_toboolean( L, i) ? "true" : "false"); break;
                case LUA_TNUMBER:  printf(" %i:%g;",   i, lua_tonumber ( L, i)); break;
                default:           printf(" %i:%s;",   i, lua_typename ( L, t)); break;
            }
        }
        printf("}\n");
    }

    bool getBool(lua_State* L, int i){ return (bool)lua_toboolean(L, i); }
    bool getBool(lua_State* L  )     { return getBool(L, -1); }

    double getDouble(lua_State* L, int i ) {
        if(!lua_isnumber(L, i)) { getError( i, "Not a double"); }
        return lua_tonumber(L, i);
    }
    double getDouble(lua_State* L ) { return getDouble(L, -1); };

    //int toInt(double f){ long li=floor(f);      }

    long getInt(lua_State* L, int i) {
        long li=0;
        if(lua_isnumber(L, i)) {
            double f = lua_tonumber(L, i);
            //printf( "%f \n", f );
            li = floor(f);
            if( li != f ) getError( i, "Not an integer");
        }else{
            getError( i, "Not a number");
        }
        return li;
    }
    long getInt( lua_State* L ) { return getInt(L, -1); };

    const char * getString(lua_State* L, int i) {
        const char * s = NULL;
        if(lua_isstring(L, i)) {  return lua_tostring(L, i);       }
        else                   {  getError( i, "Not a string"); return NULL; }
    }
    const char * getString(lua_State* L){ return getString(L,-1); };

    double getDoubleField(lua_State* L, const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        //if (!lua_isnumber(L, -1)) getErro( -1, "field not a number" );
        //int result = (int)lua_tonumber(L, -1);
        double f= getDouble(L,-1);
        lua_pop(L, 1);        // remove number
        return f;
    }

    long getIntField(lua_State* L, const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        long li= getDouble(L,-1);
        lua_pop(L, 1);        // remove number
        return li;
    }

    const char * getStringField( lua_State* L, const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        //if (!lua_isnumber(L, -1)) getErro( -1, "field not a number" );
        //int result = (int)lua_tonumber(L, -1);
        const char * s= getString(L,-1);
        lua_pop(L, 1);        // remove number
        return s;
    }

    void getLuaArr( lua_State* L, int n, float * arr, int idx=-1 ){
        luaL_checktype(L, idx, LUA_TTABLE);
        for(int i=0; i<n; i++){
            lua_rawgeti(L, idx, i+1); arr[i] = lua_tonumber(L, -1); lua_pop(L, 1);
        }
    }

    void getLuaArr( lua_State* L,  int n, int * arr, int idx=-1 ){
        luaL_checktype(L, idx, LUA_TTABLE);
        for(int i=0; i<n; i++){
            lua_rawgeti(L, idx, i+1); arr[i] = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        }
    }

    Quat4d getVec4( lua_State* L, int idx){
        Quat4d v;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); v.x = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); v.y = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 3); v.z = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 4); v.z = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_pop(L, 3);
        return v;
    }
    Quat4d getVec4( lua_State* L){ return getVec4(L, -1); }

    Vec3d getVec3( lua_State* L, int idx){
        Vec3d v;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); v.x = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); v.y = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 3); v.z = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_pop(L, 3);
        return v;
    }
    Vec3d getVec3( lua_State* L){  return getVec3(L, -1 ); }

    Vec2d getVec2( lua_State* L, int idx){
        Vec2d v;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); v.x = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); v.y = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_pop(L, 3);
        return v;
    }
    Vec2d getVec2( lua_State* L){ return getVec2(L,-1); }

    Mat3d getMat3( lua_State* L, int idx){
        Mat3d mat;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_pushinteger(L, 1); lua_gettable(L, idx); mat.a=getVec3(L,-1); lua_pop(L, 1);
        lua_pushinteger(L, 2); lua_gettable(L, idx); mat.b=getVec3(L,-1); lua_pop(L, 1);
        lua_pushinteger(L, 3); lua_gettable(L, idx); mat.c=getVec3(L,-1); lua_pop(L, 1);
        return mat;
    }
    Mat3d getMat3(lua_State* L){ return getMat3(L, -1); };

    void getDoubleVector( lua_State* L, std::vector<double>& vec ) {
        luaL_checktype(L, -1, LUA_TTABLE);
        lua_pushnil(L);
        while(lua_next(L, -2)) { vec.push_back( lua_tonumber(L, -1) ); }
        clean(L);
    }

    Quat4i getVec4i( lua_State* L, int idx){
        Quat4i v;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        //double x,y,z,w;
        lua_rawgeti(L, idx, 1); v.x = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); v.y = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        lua_rawgeti(L, idx, 3); v.z = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        lua_rawgeti(L, idx, 4); v.w = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        //lua_rawgeti(L, idx, 1); x = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_rawgeti(L, idx, 2); y = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_rawgeti(L, idx, 3); z = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_rawgeti(L, idx, 4); w = lua_tonumber(L, -1); lua_pop(L, 1);
        //vec.x = floor(x); vec.y = floor(y); vec.z = floor(z); vec.w = floor(w);
        //printf( "getVec4i() (%g,%g,%g,%g) -> (%i,%i,%i,%i)", x,y,z,w, vec.x, vec.y, vec.z, vec.w );
        //lua_pop(L, 3);
        return v;
    }
    Quat4i getVec4i( lua_State* L){ return getVec4i(L, -1); }

    Vec3i getVec3i( lua_State* L, int idx){
        Vec3i v;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); v.x = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); v.y = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        lua_rawgeti(L, idx, 3); v.z = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        //lua_pop(L, 3);
        return v;
    }
    Vec3i getVec3i( lua_State* L){ return getVec3i(L, -1); }

    Vec2i getVec2i( lua_State* L, int idx){
        Vec2i v;
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); v.x = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); v.y = floor(lua_tonumber(L, -1)); lua_pop(L, 1);
        //lua_pop(L, 3);
        return v;
    }
    Vec2i getVec2i( lua_State* L){ return getVec2i(L, -1); }

    void pushVec3( lua_State* L, const Vec3d& vec ){
        lua_createtable(L, 3, 0);
        lua_pushnumber(L, vec.x); lua_rawseti(L, -2, 1); // 
        lua_pushnumber(L, vec.y); lua_rawseti(L, -2, 2);
        lua_pushnumber(L, vec.z); lua_rawseti(L, -2, 3);
    }

    bool getVar(lua_State* L, const char * nameStr) {
        int level = 0;
        int i0 = 0;
        char wstr[64];
        for(int i=0; i<256; i++) {
            char ch = nameStr[i];
            if( (ch=='.') || (ch=='\n') ){
                strncpy(wstr, nameStr+i0, i-i0 );
                if(level == 0)       { lua_getglobal(L,     wstr); }
                else                 { lua_getfield (L, -1, wstr); }
                if(lua_isnil(L, -1)) { printf( "LuaERROR: %s of %s not found \n", wstr, nameStr ); return false; }
                else                 { level++; }
            }
        }
        return true; // success
    }

    bool dofile(lua_State* L, const char* filename){
        bool ret=true;
        if( fileExist(filename ) ){
            int ret = luaL_dofile(L,filename);
            checkError( L, ret );
            if(ret==LUA_OK) return false; //exit(0); 
        }else{
            printf( "Lua::dofile() ERROR: file not found %s\n", filename );
        }
        return ret;
    }

/*


// http://cc.byexamples.com/2008/11/19/lua-stack-dump-for-c/
void dumpStack(lua_State* L){
    int n = lua_gettop( L);  printf("LuaStack[%i] = {",n);
    for (int i=1; i<=n; i++){
        int t = lua_type( L, i);
        switch (t) {
            //case LUA_TSTRING:  printf("%i: str:  '%s';",i, lua_tostring ( L, i)); break;
            //case LUA_TBOOLEAN: printf("%i bool: %s;",  i, lua_toboolean( L, i) ? "true" : "false"); break;
            //case LUA_TNUMBER:  printf("%i num:  %g;",  i, lua_tonumber ( L, i)); break;
            case LUA_TSTRING:  printf(" %i:'%s';", i, lua_tostring ( L, i)); break;
            case LUA_TBOOLEAN: printf(" %i:%s;",   i, lua_toboolean( L, i) ? "true" : "false"); break;
            case LUA_TNUMBER:  printf(" %i:%g;",   i, lua_tonumber ( L, i)); break;
            default:           printf(" %i:%s;",   i, lua_typename ( L, t)); break;
        }
    }
    printf("}\n");
}

// https://www.lua.org/pil/25.1.html
// https://stackoverflow.com/questions/11974806/accessing-a-lua-table-within-a-table-from-c
// " It is generally easier to use functions lua_getfield (for string indexes) or lua_rawgeti (for numerical indexes) than the raw lua_gettable function. "
// https://stackoverflow.com/questions/2705666/lua-how-to-check-if-one-of-the-values-associated-with-the-specified-key-of-a-t
int getfield (lua_State* L, const char *key) {
    lua_pushstring(L, key);
    lua_gettable(L, -2);  // get background[key]
    if (!lua_isnumber(L, -1)) printf(  "Lua ERROR: invalid component in background color\n" );
    int result = (int)lua_tonumber(L, -1);
    lua_pop(L, 1);        // remove number
    return result;
}

// https://stackoverflow.com/questions/41387796/access-nested-tables-from-lua-to-c-to-get-values
// https://stackoverflow.com/questions/25940366/passing-array-to-c-as-argument-in-the-stack

void getVec3(lua_State *L, int idx, Vec3d& vec){
    // universal helper function to get Vec3 function argument from Lua to C++ function
    luaL_checktype(L, idx, LUA_TTABLE);
    lua_rawgeti(L, idx, 1); vec.x = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_rawgeti(L, idx, 2); vec.y = lua_tonumber(L, -1); lua_pop(L, 1);
    lua_rawgeti(L, idx, 3); vec.z = lua_tonumber(L, -1); lua_pop(L, 1);
    //lua_pop(L, 3);
}

void getMat3(lua_State *L, int idx, Mat3d& mat){
    // universal helper function to get Vec3 function argument from Lua to C++ function

    luaL_checktype(L, idx, LUA_TTABLE);
    lua_pushinteger(L, 1); lua_gettable(L, idx); lua_getVec3(L, -1, mat.a ); lua_pop(L, 1);
    lua_pushinteger(L, 2); lua_gettable(L, idx); lua_getVec3(L, -1, mat.b ); lua_pop(L, 1);
    lua_pushinteger(L, 3); lua_gettable(L, idx); lua_getVec3(L, -1, mat.c ); lua_pop(L, 1);

    //luaL_checktype(L, idx, LUA_TTABLE);   printf("DEBUG 0\n");
    //lua_pushinteger(L, 1); lua_gettable(L, idx); printf("DEBUG 0.1\n"); lua_getVec3(L, -1, mat.a ); printf("DEBUG 1\n"); lua_pop(L, 1);
    //lua_pushinteger(L, 2); lua_gettable(L, idx); printf("DEBUG 1.1\n"); lua_getVec3(L, -1, mat.b ); printf("DEBUG 2\n"); lua_pop(L, 1);
    //lua_pushinteger(L, 3); lua_gettable(L, idx); printf("DEBUG 2.1\n"); lua_getVec3(L, -1, mat.c ); printf("DEBUG 3\n"); lua_pop(L, 1);
    //lua_rawgeti(L, idx, 1); lua_gettable(L, idx); printf("DEBUG 0.1\n"); lua_getVec3(L, -1, mat.a ); printf("DEBUG 1\n"); lua_pop(L, 1);
    //lua_rawgeti(L, idx, 2); lua_gettable(L, idx); printf("DEBUG 1.1\n"); lua_getVec3(L, -1, mat.b ); printf("DEBUG 2\n"); lua_pop(L, 1);
    //lua_rawgeti(L, idx, 3); lua_gettable(L, idx); printf("DEBUG 2.1\n"); lua_getVec3(L, -1, mat.c ); printf("DEBUG 3\n"); lua_pop(L, 1);
}
*/


}

#endif

