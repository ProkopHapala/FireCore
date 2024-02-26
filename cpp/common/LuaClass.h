#ifndef LuaClass_h
#define LuaClass_h

#include <lua5.2/lua.hpp>
#include "LuaUtils.h"


class LuaClass{ public:
    static lua_State* L=0;  // don't forget to initialize it
    bool* atomFixed = 0;

    // void l_fixAtom(int ia){
    //     int ia = Lua::getInt(L,1); 
    //     atomFixed[ia] = true;
    // }
    // void l_unfixAtom(){
    //     int ia = Lua::getInt(L,1); 
    //     atomFixed[ia] = false;
    // }

    static void registerLuaMethod(){
        lua_register(L, "fixAtom",this->l_fixAtoms);
    }

    template<typename Lamb>
    void registerLuaMethod(Lamb lambda, char const* name ) {
        *(decltype(lambda)**)lua_newuserdata(L, sizeof(lambda)) = lambda;
        lua_pushcclosure(L, [](lua_State* L) -> int {
            auto lambdaPtr = *(decltype(lambda)**)lua_touserdata(L, lua_upvalueindex(1));
            return (*lambdaPtr)(L);
        }, 1); 
        lua_setglobal(L, name);
    }

    // void registerLuaMethods(){
    //     registerLuaMethod( [this](lua_State* L)->int { this->l_fixAtoms();  return 0; }, "fixAtom"   );
    //     registerLuaMethod( [this](lua_State* L)->int { this->l_unfixAtom(); return 0; }, "unfixAtom" );
    // }

    // =================== Old Lua Functions

    void getError( int i,  const char * s ){ printf( "LuaERROR[%i]: %s\n",s); }
    void clean   () { int n = lua_gettop(L); lua_pop(L, n); }


    void dumpStack(){  // print the stack - for debugging
        int n = lua_gettop( L);     printf("LuaStack[%i] = {",n);
        for (int i=1; i<=n; i++){
            int t = lua_type( L, i);
            switch (t) {
                case LUA_TSTRING:   printf(" %i:'%s';", i, lua_tostring ( L, i)); break;
                case LUA_TBOOLEAN:  printf(" %i:%s;",   i, lua_toboolean( L, i) ? "true" : "false"); break;
                case LUA_TNUMBER:   printf(" %i:%g;",   i, lua_tonumber ( L, i)); break;
                default:            printf(" %i:%s;",   i, lua_typename ( L, t)); break;
            }
        }
        printf("}\n");
    }

    bool getBool(int i){ return (bool)lua_toboolean(L, i); }
    bool getBool()     { return getBool(-1); }

    double getDouble(int i ) {
        if(!lua_isnumber(L, i)) { getError( i, "Not a double"); }
        return lua_tonumber(L, i);
    }
    double getDouble() { return getDouble(-1); };

    long getInt(int i) {
        long li=0;
        if(lua_isnumber(L, i)) {
            double f = lua_tonumber(L, i);
            long   li = floor(f);
            if( li != f ) getError( i, "Not an integer");
        }else{
            getError( i, "Not a number");
        }
        return li;
    }
    long getInt() { return getInt(-1); };

    const char * getString(int i) {
        const char * s = NULL;
        if(lua_isstring(L, i)) {  return lua_tostring(L, i);       }
        else                   {  getError( i, "Not a string"); return NULL; }
    }
    const char * getString(){ return getString(-1); };

    double getDoubleField( const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        //if (!lua_isnumber(L, -1)) getErro( -1, "field not a number" );
        //int result = (int)lua_tonumber(L, -1);
        double f= getDouble(-1);
        lua_pop(L, 1);        // remove number
        return f;
    }

    long getIntField( const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        long li= getDouble(-1);
        lua_pop(L, 1);        // remove number
        return li;
    }

    const char * getStringField( const char *key) {
        lua_pushstring(L, key);
        lua_gettable(L, -2);  // get background[key]
        //if (!lua_isnumber(L, -1)) getErro( -1, "field not a number" );
        //int result = (int)lua_tonumber(L, -1);
        const char * s= getString(-1);
        lua_pop(L, 1);        // remove number
        return s;
    }

    void getVec3( int idx, Vec3d& vec){
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_rawgeti(L, idx, 1); vec.x = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 2); vec.y = lua_tonumber(L, -1); lua_pop(L, 1);
        lua_rawgeti(L, idx, 3); vec.z = lua_tonumber(L, -1); lua_pop(L, 1);
        //lua_pop(L, 3);
    }
    void getVec3( Vec3d& vec){  getVec3( -1, vec); }

    void getMat3( int idx, Mat3d& mat){
        // universal helper function to get Vec3 function argument from Lua to C++ function
        luaL_checktype(L, idx, LUA_TTABLE);
        lua_pushinteger(L, 1); lua_gettable(L, idx); getVec3(-1, mat.a ); lua_pop(L, 1);
        lua_pushinteger(L, 2); lua_gettable(L, idx); getVec3(-1, mat.b ); lua_pop(L, 1);
        lua_pushinteger(L, 3); lua_gettable(L, idx); getVec3(-1, mat.c ); lua_pop(L, 1);
    }
    void getMat3( Mat3d& mat){ getMat3( -1, mat); };

    void getDoubleVector( std::vector<double>& vec ) {
        luaL_checktype(L, -1, LUA_TTABLE);
        lua_pushnil(L);
        while(lua_next(L, -2)) { vec.push_back( lua_tonumber(L, -1) ); }
        clean();
    }

*/

    virtual void init(){ L = luaL_newstate(); };
};

#endif

