print("LUA: makeGUI.lua START")
clear_gui(0)  -- remove all widgets from the GUI


-- const char* label     = Lua::getString(L,1);
-- Vec3i xyw             = Lua::getVec3i (L,2);
-- Vec3i SliderButtonInt = Lua::getVec3i (L,3);
-- Vec3d MinMaxCur       = Lua::getVec3  (L,4);
-- const char* command   = Lua::getString(L,5);

button( "perFrame", {1,2,16}, {1,1,0}, {-12.0,12.0,6.0}, "perFrame" )  -- add a button to the GUI

print("LUA: script.lua END")