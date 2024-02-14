-- 
--dofile("./common_resources/utils.lua")
print("LUA: test_add_mols.lua 1")

function A()
    print("LUA: ============================")
    print("LUA: ======= function A() =======")
    print("LUA: ============================")
    clear()
    add( "./common_resources/H2O.xyz", {0.,0.,0.} )
    add( "./common_resources/NH3.xyz", {10.,0.,0.} )
    print("LUA: ====== make() ==============")
    make()
    autoCharges()
    -- frags()
end

print("LUA: test_add_mols.lua END")