#=========================================================
 File reader    
=========================================================#


function readGraph(absolute_path)
    println("\n read graph")
    lines = readlines(absolute_path)
    line1 = split.(lines[1], " ")
    node_num = parse(Int, String(line1[1]))
    edge_num = parse(Int, String(line1[2]))
    edge_fields = Vector{Tuple{Int, Int, Float64}}()
    for line in lines[2:end]
        line = split.(line, " ")
        push!(edge_fields, (parse(Int, String(line[1])), parse(Int, String(line[2])), parse(Float64, String(line[3]))))
    end
    return Graph(node_num, edge_num, edge_fields)
end