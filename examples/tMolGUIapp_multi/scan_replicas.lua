-- 
-- dofile("./common_resources/utils.lua")

require "math"

print("LUA: scan_replicas.lua.start")

-- shiftReplica( 0, {5.,0.,0.}  );

function move_replicas()
    print("LUA: move_replicas()")
    dx   = 0.5
    dy   = 0.5
    nx   = 5  
    nSys = getReplicaCount()
    print("LUA: move_replicas() nSys=", nSys )
    for i = 1, nSys do
        iy = math.floor(i/nx)
        ix = i - iy*nx
        d  = {ix*dx,iy*dy,0.0}
        io.write( "LUA: ", i, " ", d[1],d[2],d[3] ); 
        shiftReplica( i, d );
    end
    systemsToGPU()
end

lua_button( "move_replicas()", {200,400,16}, {1,1,0},{0.0,1.0,0.0}, "move_replicas()" )

-- move_replicas()

print("LUA: scan_replicas.lua.end")