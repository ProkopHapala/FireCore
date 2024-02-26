-- 
dofile("./common_resources/utils.lua")

print("LUA: script.lua 1")

na = natom()

print("LUA: script.lua 2")

print(na)

print("LUA: script.lua 3")

function print_aposs()
    print("LUA: print_aposs()")
    na = natom()
    print("LUA: print_aposs() na=", na)
    for i = 1, na do
        pi = apos(i)
        -- print(i, pi )
        io.write(i, " "); printv( pi )
    end
end

command("print_aposs()")

print("LUA: script.lua 4")