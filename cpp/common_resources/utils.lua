
function printv( t )
    -- no new line
    io.write("{" ); 
    for k, v in pairs(t) do
        io.write( tostring(v), " " )
    end
    io.write("}\n" );
end

function printTable(t, indent)
    indent = indent or ""
    if type(t) ~= "table" then
        print(indent .. tostring(t))
    else
        for k, v in pairs(t) do
            if type(v) == "table" then
                print(indent .. tostring(k) .. ":")
                printTable(v, indent .. "  ")
            else
                print(indent .. tostring(k) .. ": " .. tostring(v))
            end
        end
    end
end

