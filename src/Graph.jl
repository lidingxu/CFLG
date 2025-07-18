#=========================================================
 Graph
=========================================================#

# Node
struct Node
    # parameters
    node_id::Int # node id

    # construtor
    function Node(node_id::Int)
        new(node_id)
    end
end

# Edge
struct Edge
    # parameters
    edge_id::Int # edge id
    nodes::Dict{Symbol, Int} # end node a id, end node b id : nodea_id < nodeb_id
    length::Float64 # length
    etype::Symbol # {e_normal, e_long}

    # construtor
    function Edge(edge_id::Int, nodea_id::Int, nodeb_id::Int, length::Float64, etype::Symbol=:e_normal)
        ord = lor(nodea_id, nodeb_id)
        edge = new(edge_id, Dict([(:a,ord[1]), (:b,ord[2])]), length, etype)
        return edge
    end
end

# return whether contain a node
function contain(edge::Edge, node_id::Int)
    return edge.nodes[:a] == node_id || edge.nodes[:b] == node_id
end



# Graph
mutable struct Graph
    # parameters
    node_num::Int
    edge_num::Int
    nodes::Vector{Node} # nodes
    edges::Vector{Edge} # edges
    node_ids::Vector{Int} # node ids
    edge_ids::Vector{Int} # edge ids
    adjacent_edges::Vector{Set{Int}} # adjacent edge
    incident_nodes::Vector{Set{Int}} # incident nodes
    node2edge::Dict{Tuple{Int,Int},Int} # map from end node pair to edge
    min_len::Float64 # minimum length of edfes
    max_len::Float64 # maximum edge length
    avg_len::Float64 # average edge length
    total_len::Float64 # total length
    #leaves::Vector{Int} #  leaf nodes

    # construtor
    function Graph(node_num::Int, edge_num::Int, edge_fields)
        graph = new(node_num, edge_num)
        graph.nodes = Vector{Node}() # nodes
        graph.edges = Vector{Edge}() # edges
        graph.node_ids = Vector{Int}(1:node_num) # node ids
        graph.edge_ids = Vector{Int}(1:edge_num) # edge ids
        graph.adjacent_edges = Vector{Set{Int}}() # adjacent edges
        graph.incident_nodes = Vector{Set{Int}}() # incident nodes
        graph.node2edge = Dict{Tuple{Int,Int},Int}() # map from end node pair to edge
        graph.min_len = typemax(Float64) # minimum edge length
        graph.max_len =  typemin(Float64) # maximum edge length
        graph.avg_len =  0 # average edge length
        graph.total_len = 0
        #graph.leaves =  Vector{Int}()  #  leaf nodes

        # add nodes
        for node_id in 1:node_num
            push!(graph.nodes, Node(node_id))
        end

        # add edges, compute node2edge
        for (edge_id, edge_field) in enumerate(edge_fields)
            if length(edge_field) == 3
                push!(graph.edges, Edge(edge_id, lor(edge_field[1], edge_field[2])..., edge_field[3], :e_normal))
            elseif length(edge_field) == 4
                push!(graph.edges, Edge(edge_id, lor(edge_field[1], edge_field[2])..., edge_field[3], edge_field[4]))
            end
            graph.node2edge[lor(edge_field[1], edge_field[2])] = edge_id
            graph.min_len = min(graph.min_len, edge_field[3])
            graph.max_len = max(graph.max_len, edge_field[3])
            graph.avg_len += edge_field[3]
            graph.total_len += edge_field[3]
        end
        graph.avg_len /= graph.edge_num

        # compute adjacent edges
        for node_id in 1: node_num
            push!(graph.adjacent_edges, adjacentEdges(graph, node_id))
        end

        # compute incident edges
        for node_id in 1: node_num
            push!(graph.incident_nodes, incidentNodes(graph, node_id))
        end
        return graph
    end
end

# print graph
function printGraph(graph::Graph)
    println(graph.node_num, " ", graph.edge_num)
    println(graph.node_ids)
    println(graph.edge_ids)
    println(graph.node2edge)
    println(graph.nodes)
    println(graph.edges)
end

# return edges ajdacent to a node
function adjacentEdges(graph::Graph, node_id::Int)
    edge_ids = Set{Int}()
    for edge in graph.edges
        if contain(edge, node_id)
            push!(edge_ids, edge.edge_id)
        end
    end
    return edge_ids
end

# return nodes incident to a node
function incidentNodes(graph::Graph, node_id::Int)
    node_ids = Set{Int}()
    for edge in graph.edges
        if contain(edge, node_id)
            if edge.nodes[:a] == node_id
                push!(node_ids, edge.nodes[:b])
            else
                push!(node_ids, edge.nodes[:a])
            end
        end
    end
    return node_ids
end

# return whether two nodes are incident
function isIncident(graph::Graph, node_id::Int, node_id_::Int)
    return node_id_ in  incidentNodes(graph, node_id)
end

# return wether two edges are adjacent
function isAdjacent(graph::Graph, edge_id::Int, edge_id_::Int)
    return contain(graph.edges[edge_id], graph.edges[edge_id_].nodes[:a]) ||  contain(graph.edges[edge_id], graph.edges[edge_id_].nodes[:b])
end

# nodeCover, mode: :partial: stop nodeCover earlier, :full : no stop, used as shotest distance formulationrithm
function nodeCover!(graph::Graph, s_id::Int, dlt::Float64, d::Dict{Tuple{Int, Int}, Float64}, mode::Symbol)
    Ec = Set{Int}() # complete covered edges
    E = Set{Int}() # affected edges
    V = Set{Int}() # covered nodes
    U = Set{Int}() # uncertain nodes
    Q = PriorityQueue{Int, Float64}() # min heap
    #prev = Dict{Int,Int}() # previous nodes


    for v_id in graph.node_ids
        Q[v_id] = typemax(Float64)
        #prev[v_id] = -1
    end
    # self cover
    #prev[s_id]=s_id
    Q[s_id] = 0
    push!(V, s_id)
    while !isempty(Q)
        # current nearest point
        (u_id, dist) = peek(Q)
        dequeue!(Q)
        if mode == :Partial && dist > dlt # end the formulationrithm
            break
        end
        #@assert(u_id in V)
        d[lor(s_id,u_id)] = min(d[lor(s_id,u_id)], dist) # update distance
        for v_id in graph.incident_nodes[u_id] # incident nodes
            if !haskey(Q, v_id) # only unrelaxed node
                continue
            end
            dist_ = Q[v_id]
            # consider an edge
            e_id = graph.node2edge[lor(u_id, v_id)]
            push!(E, e_id) # affected edge
            #@assert(u_id in V || v_id in V)
            r = dist+ graph.edges[e_id].length # new length
            if r <= dist_ # update distance
                Q[v_id] = r
                #prev[v_id] = u_id
                if r <= dlt
                    push!(Ec,e_id) # edge completed cover
                    push!(V, v_id) # node completed cover
                else
                    push!(U, e_id) # uncertain edge
                end
            end
        end
    end
    for e_id in U
        va_id = graph.edges[e_id].nodes[:a]
        vb_id = graph.edges[e_id].nodes[:b]
        if va_id in V && vb_id in V && 2*dlt - d[lor(s_id, va_id)] -  d[lor(s_id, vb_id)] >= graph.edges[e_id].length
            push!(Ec,e_id) # e is completed covered
        end
    end
    return Ec, E, V
end


# mutual
function mutual(dlt::Float64, e::Edge, ef::Edge, d::Dict{Tuple{Int, Int}, Float64}, cr_tol::Float64, c_tol::Float64)
    va_id = e.nodes[:a]
    vb_id = e.nodes[:b]
    vfa_id = ef.nodes[:a]
    vfb_id = ef.nodes[:b]
    Qlb = typemax(Float64)
    Qub = typemin(Float64)
    len = ef.length
    rs = Dict()
    function r_idf(q::Float64, defa::Float64, defb::Float64, Qlb_id::Float64, Qub_id::Float64)
        if q <= Qlb_id
            return dlt - (defa + q)
        elseif q <= Qub_id
            return 0
        else
            return dlt - (defb + len  - q)
        end
    end
    for v_id in  [va_id, vb_id]
        defa = d[lor(v_id, vfa_id)]
        defb = d[lor(v_id, vfb_id)]
        Q_id = (defb + len- defa) / 2
        Qlb_id = max(min(Q_id, dlt - defa), 0)
        Qub_id = min(max(Q_id, defb + len - dlt),len)
        Qlb = min(Qlb, Qlb_id)
        Qub = max(Qub, Qub_id)
        r_id(q) = r_idf(q, defa, defb, Qlb_id, Qub_id)
        rs[v_id] = r_id
    end
    r(q) = rs[va_id](q) + rs[vb_id](q)
    if r(Qlb) >= len * (1+cr_tol) + c_tol && r(Qub) >= len * (1+cr_tol) + c_tol
        return true
    else
        return false
    end
end


# absorb degree-2 nodes of graph (experimental)
function absorbGraph(graph::Graph)
    v_ids  = Set{Int}(graph.node_ids)
    e_ids = Set{Int}(graph.edge_ids)
    e_lens = Dict{Int, Float64}()
    v2e  = Dict{Tuple{Int,Int}, Float64}()
    inc =  Dict{Int, Set{Int}}()
    degree2  = Set{Int}()
    # init
    for node in graph.nodes
        if length(graph.adjacent_edges[node.node_id]) == 2
            push!(degree2, node.node_id)
        end
        inc[node.node_id] = Set{Int}()
    end
    # graph data
    for edge in graph.edges
        e_lens[edge.edge_id] = edge.length
        v2e[lor(edge.nodes[:a],edge.nodes[:b])] = edge.edge_id
        push!(inc[edge.nodes[:a]],  edge.nodes[:b])
        push!(inc[edge.nodes[:b]],  edge.nodes[:a])
    end

    absorb = true
    enew_id  = length(e_ids)
    while absorb
        absorb = false
        for v_id in degree2
            if !(v_id in v_ids)
                continue
            end
            incv = [v for v in inc[v_id]]
            v1_id = incv[1]
            v2_id = incv[2]
            if v1_id in inc[v2_id] # avoid multi-edge
                continue
            else
                # delete edges
                e1_id = v2e[lor(v_id, v1_id)]
                e2_id = v2e[lor(v_id, v2_id)]
                delete!(e_ids, e1_id)
                delete!(e_ids, e2_id)
                # delete nodes
                delete!(v_ids, v_id)
                delete!(inc[v1_id], v_id)
                delete!(inc[v2_id], v_id)
                # add new edges
                enew_id += 1
                push!(e_ids, enew_id)
                v2e[lor(v1_id, v2_id)] = enew_id
                e12_len = e_lens[e1_id] + e_lens[e2_id]
                e_lens[enew_id] = e12_len
                push!(inc[v1_id], v2_id)
                push!(inc[v2_id], v1_id)
                absorb = true
                break
            end
        end
    end


    node_num = length(v_ids)
    edge_num = length(e_ids)
    edge_fields = Vector{Tuple{Int, Int, Float64}}()

    map = Dict{Int, Int}()

    for (i, v_id) in enumerate(v_ids)
        map[v_id] = i
    end

    for v1_id in v_ids
        for v2_id in inc[v1_id]
            if v1_id < v2_id
                continue
            end
            e_id = v2e[lor(v1_id, v2_id)]
            push!(edge_fields, (map[v1_id], map[v2_id], e_lens[e_id]))
        end
    end

    absorbgraph = Graph(node_num, edge_num, edge_fields)
    println("\nafter absorb:",absorbgraph.node_num, " ", graph.node_num)
    @assert(abs(absorbgraph.total_len - graph.total_len) < 1e-5)
    return absorbgraph
end

# prune a graph to a smaller graph
function projectGraph(graph::Graph, VR, R, required, rho, dlt)

    # Construct the set of nodes for the induced graph
    induced_nodes = Set{Int}()
    for node in graph.nodes
        if !(node.node_id in VR)
            push!(induced_nodes, node.node_id)
        end
    end

    induced_nodes = union(induced_nodes, R)

    # Map old node ids to new node ids (1-based)
    node_map = Dict{Int, Int}()
    for (i, v_id) in enumerate(induced_nodes)
        node_map[v_id] = i
    end

    # Collect edges where both endpoints are in induced_nodes
    edge_fields = Vector{Tuple{Int, Int, Float64, Symbol}}()
    for edge in graph.edges
        a = edge.nodes[:a]
        b = edge.nodes[:b]
        if a in induced_nodes && b in induced_nodes
            push!(edge_fields, (node_map[a], node_map[b], edge.length, edge.etype))
        end
    end

    # Add edges from each v in R to a new node with length rho[v]
    for v in R
        new_node_id = length(node_map) + 1
        node_map[-v] = new_node_id  # Use -v as a unique key for the new node
        push!(edge_fields, (node_map[v], new_node_id,  rho[v], :e_normal ))
    end

    node_num = length(node_map)
    edge_num = length(edge_fields)

    projectgraph = Graph(node_num, edge_num, edge_fields)

    # Construct R_, RE_, R_toR
    R_ = Vector{Int}()
    RE_ = Vector{Int}()
    R_toR = Dict{Int, Int}()
    RtoR_ = Dict{Int, Int}()
    EtoE_ = Dict{Int, Int}()

    for v in R
        newv = node_map[v]
        push!(R_, newv)
        # Find the edge index in projectgraph that connects v to new_node_id
        for eidx in projectgraph.adjacent_edges[newv]
            edge = projectgraph.edges[eidx]
            if (edge.nodes[:a] == newv && edge.nodes[:b] == node_map[-v]) ||
               (edge.nodes[:b] == newv && edge.nodes[:a] == node_map[-v])
                push!(RE_, eidx)
                break
            end
        end
        R_toR[newv] = v
        RtoR_[v] = newv
    end

    for edge in graph.edges
        if edge.nodes[:a] in induced_nodes && edge.nodes[:b] in induced_nodes
            EtoE_[edge.edge_id] = projectgraph.node2edge[lor(node_map[edge.nodes[:a]], node_map[edge.nodes[:b]])]
        end
    end

    @assert( length(RE_) == length(R_) )
    # checking
    num_e = 0
    num_v = 0
    for edge in projectgraph.edges
        num_e = max(num_e, edge.edge_id)
        num_v = max(max(num_v, edge.nodes[:a]), edge.nodes[:b])
    end

    num_v_ = 0

    for node in projectgraph.nodes
        num_v_  =max(num_v_, node.node_id)
    end

    @assert(num_v == num_v_ && num_v_ == projectgraph.node_num && length(projectgraph.nodes) == num_v && num_e == projectgraph.edge_num)
    return projectgraph, R_, RE_, R_toR, RtoR_, EtoE_
end



# break a graph to a smaller graph
function breakGraph(graph::Graph, dlt::Float64, isabsorb::Bool = true, long::Bool = false)
    if dlt >= graph.max_len
        return graph
    end

    if isabsorb
        graph = absorbGraph(graph)
    end


    node_num = graph.node_num
    edge_num = graph.edge_num

    max_piece = -1
    if long # long edge
        edge_fields = Vector{Tuple{Int, Int, Float64, Symbol}}()
        for edge in graph.edges
            if edge.length < dlt
                push!(edge_fields, Tuple{Int, Int, Float64, Symbol}([edge.nodes[:a], edge.nodes[:b], edge.length, :e_normal]))
                continue
            end
            piece_num = Int(ceil(edge.length / (2 * dlt))) # to do int
            if piece_num >= 2
                push!(edge_fields, Tuple{Int, Int, Float64, Symbol}([edge.nodes[:a], edge.nodes[:b], edge.length, :e_long]))
            else # still break
                piece_num = Int(ceil(edge.length / dlt))
                v_ids = Vector{Int}( (node_num + 1) : (node_num + piece_num - 1) )
                len_last = edge.length - (piece_num - 1) * (dlt * (1- 1e-6))
                @assert(len_last < dlt)
                for v_id in v_ids
                    if v_id == node_num + 1
                        push!(edge_fields, Tuple{Int, Int, Float64, Symbol}([edge.nodes[:a], v_id, dlt * (1- 1e-6), :e_normal]))
                    else
                        push!(edge_fields, Tuple{Int, Int, Float64, Symbol}([v_id - 1, v_id, dlt * (1- 1e-6), :e_normal]))
                    end
                end
                push!(edge_fields, Tuple{Int, Int, Float64, Symbol}([node_num + piece_num - 1, edge.nodes[:b], len_last, :e_normal]))
                node_num += piece_num - 1
                edge_num += piece_num - 1
            end
            max_piece = max(piece_num * 2, max_piece)
        end
    else
        edge_fields = Vector{Tuple{Int, Int, Float64}}()
        for edge in graph.edges
            if edge.length < dlt
                push!(edge_fields, Tuple{Int, Int, Float64}([edge.nodes[:a], edge.nodes[:b], edge.length]))
                continue
            end
            piece_num = Int(ceil(edge.length / dlt)) # to do int
            len = edge.length / piece_num
            max_piece = max(piece_num, max_piece)
            v_ids = Vector{Int}( (node_num + 1) : (node_num + piece_num - 1) )
            #len_last = edge.length - (piece_num - 1) * (dlt * (1- 1e-9))
            #print( edge.length, " ", piece_num, " ", len_last, " ", dlt)
            #@assert(len_last < dlt)
            for v_id in v_ids
                if v_id == node_num + 1
                    push!(edge_fields, Tuple{Int, Int, Float64}([edge.nodes[:a], v_id, len]))
                else
                    push!(edge_fields, Tuple{Int, Int, Float64}([v_id - 1, v_id, len]))
                end
            end
            push!(edge_fields, Tuple{Int, Int, Float64}([node_num + piece_num - 1, edge.nodes[:b], len]))
            node_num += piece_num - 1
            edge_num += piece_num - 1
        end
    end

    println("\nmax_piece:", max_piece)

    breakgraph = Graph(node_num, edge_num, edge_fields)

    # checking
    num_e = 0
    num_v = 0
    for edge in breakgraph.edges
        num_e = max(num_e, edge.edge_id)
        num_v  = max(max(num_v, edge.nodes[:a]), edge.nodes[:b])
    end

    num_v_ = 0

    for node in breakgraph.nodes
        num_v_  =max(num_v_, node.node_id)
    end

    @assert(num_v == num_v_ && num_v_ == breakgraph.node_num && length(breakgraph.nodes) == num_v && num_e == breakgraph.edge_num)
    return breakgraph
end


# graph process for vertex formulation, notice dlt should greater than edge length
function processGraph(graph::Graph, dlt::Float64, mode::Symbol, cr_tol::Float64, c_tol::Float64)
    #@assert(dlt >= graph.max_len)
    Ec = Dict{Int, Set{Int}}() # edges can completely cover edges
    Vc = Dict{Int, Set{Int}}() # nodes can completely cover edges
    Ep = Dict{Int, Set{Int}}() # edges can partially cover incident edges of nodes
    Vp = Dict{Int, Set{Int}}() # nodes can partially cover incident edges of nodes
    EIp = Dict{Int, Set{Tuple{Int, Symbol}}}() # edges, nodes can partially cover incident edges of nodes
    EIc = Dict{Int, Set{Tuple{Int, Symbol}}}() # edges, nodes can partially cover incident edges of nodes
    d = Dict{Tuple{Int, Int}, Float64}() # distance record

    # initilization
    for v_id in graph.node_ids
        for v_id_ in graph.node_ids
            if v_id <= v_id_
                d[lor(v_id,v_id_)] = typemax(Float64)
            end
        end
        Ep[v_id] = Set{Int}()
        Vp[v_id] = Set{Int}()
        EIp[v_id] = Set{Tuple{Int,Symbol}}()
        EIc[v_id] = Set{Tuple{Int,Symbol}}()
    end

    for e_id in graph.edge_ids
        if graph.edges[e_id].etype == :e_normal
            Ec[e_id] = Set{Int}([e_id]) # notice, we assume dlt >=  len
        elseif graph.edges[e_id].etype == :e_long
            Ec[e_id] = Set{Int}()
        end
        Vc[e_id] = Set{Int}()
    end

    EcV = Dict{Int, Set{Int}}() # nodes can completely cover edges
    EV = Dict{Int, Set{Int}}() # nodes can affect edges
    VV = Dict{Int, Set{Int}}() # nodes can completely covers nodes

    # compute Vc
    for v_id in graph.node_ids
        (Ec_, E_, V_) = nodeCover!(graph, v_id, dlt, d, mode) # run node cover
        # record
        EcV[v_id] = Ec_
        EV[v_id] = E_
        VV[v_id] = V_
        for e_id in Ec_
            push!(Vc[e_id], v_id) # nodes are completely covered
        end
    end

    # compute Ec
    for e_id in graph.edge_ids
        for ef_id in graph.edge_ids
            if ef_id <= e_id  # notice, we include self cover already
                continue
            end
            e = graph.edges[e_id]
            ef = graph.edges[ef_id]
            if !(ef_id in EcV[e.nodes[:a]] && ef_id in EcV[e.nodes[:b]])
                continue
            end
            if mutual(dlt, e, ef, d, cr_tol, c_tol)
                push!(Ec[e_id], ef_id)
                push!(Ec[ef_id], e_id)
            end
        end
    end


    # compute Vp
    for v_id in graph.node_ids
        for vf_id in VV[v_id]
            for e_id in graph.adjacent_edges[v_id] # adjacent edges of v_id
                if !(e_id in EcV[vf_id]) # not c. covered by vf_id
                    push!(Vp[v_id], vf_id) # vf_id partiall cover adjacent edges of node_id
                    break
                end
            end
        end
    end

    # compute Ep, EIp
    for v_id in graph.node_ids
        for ef_id in EV[v_id] # affected edges
            #println("\n", EV[v_id])
            vfa_id =  graph.edges[ef_id].nodes[:a]
            vfb_id =  graph.edges[ef_id].nodes[:b]
            len = graph.edges[ef_id].length
            #println(d[lor(node_id, nodefa_id)], " ", d[lor(node_id, nodefb_id)]," ", len, " ", nodefa_id, " ", nodefb_id, " ", node_id)
            for e_id in graph.adjacent_edges[v_id]
                if !(ef_id in Ec[e_id]) # not completed cover
                    #println(d[lor(node_id, nodefa_id)], d[lor(node_id, nodefb_id)], len)
                    #if !(v_id in  VV[vfa_id] || v_id in VV[vfb_id])
                    #    println("\n", v_id," ",vfa_id," ",vfb_id," ",d[lor(v_id, vfa_id)], " ", d[lor(v_id, vfb_id)], " ",dlt)
                    #    println(vfa_id in VV[v_id], vfb_id in VV[v_id], ef_id in EV[v_id])
                    #@end
                    #@assert(v_id in VV[vfa_id] || v_id in VV[vfb_id])
                    #@assert(d[lor(v_id, vfa_id)] <= dlt*(1+cr_tol) +c_tol || d[lor(v_id, vfb_id)] <= dlt*(1+cr_tol) +c_tol)
                    if d[lor(v_id, vfa_id)] <= dlt*(1+cr_tol) +c_tol && d[lor(v_id, vfa_id)] <= (d[lor(v_id, vfb_id)] + len)*(1+cr_tol) +c_tol
                        push!(EIp[v_id], (ef_id, :a))
                        push!(Ep[v_id],ef_id)
                    end
                    if d[lor(v_id, vfb_id)] <= dlt*(1+cr_tol) +c_tol && d[lor(v_id, vfb_id)]  <= (d[lor(v_id, vfa_id)] + len)*(1+cr_tol) +c_tol
                        push!(EIp[v_id], (ef_id, :b))
                        push!(Ep[v_id],ef_id)
                    end
                    break
                end
            end
            if d[lor(v_id, vfa_id)] <= dlt*(1+cr_tol) +c_tol && d[lor(v_id, vfa_id)] <= (d[lor(v_id, vfb_id)] + len)*(1+cr_tol) +c_tol
                push!(EIc[v_id], (ef_id, :a))
            end
            if d[lor(v_id, vfb_id)] <= dlt*(1+cr_tol) +c_tol && d[lor(v_id, vfb_id)]  <= (d[lor(v_id, vfa_id)] + len)*(1+cr_tol) +c_tol
                push!(EIc[v_id], (ef_id, :b))
            end
        end
    end

    return Ec, Vc, Ep, Vp, EIp, EIc, d
end


function computeDistance(graph::Graph)
    d = Dict{Tuple{Int, Int}, Float64}() # distance record
    # initilization
    for v_id in graph.node_ids
        for v_id_ in graph.node_ids
            if v_id <= v_id_
                d[lor(v_id,v_id_)] = typemax(Float64)
            end
        end
    end

    # only distance information is used
    for v_id in graph.node_ids
        nodeCover!(graph, v_id, typemax(Float64), d, :full) # run node cover in full mode
    end

    return d
end



# graph process for edge formulation
function processGraphSimple(graph::Graph, dlt::Float64)
    d = computeDistance(graph::Graph)

    # compute big M
    radius = typemin(Float64)
    for v_id in graph.node_ids
        for v_id_ in graph.node_ids
            if v_id < v_id_
                radius = max(radius, d[lor(v_id,v_id_)])
            end
        end
    end
    max_len = typemin(Float64)
    for edge in graph.edges
        max_len = max(max_len, edge.length)
    end
    bigM_EF = 2*(max_len + radius + dlt) + 1
    return bigM_EF, d
end

# Traverse a tree rooted at a given node, applying a function to each node
function traverseTree(graph, ch::Dict{Int, Set{Int}}, pr::Dict{Int, Int}, root::Int, f)
    stack = [root]
    visited = Set{Int}()
    while !isempty(stack)
        node = pop!(stack)
        if node in visited
            continue
        end
        f(graph, ch, pr, node)
        push!(visited, node)
        for child in ch[node]
            push!(stack, child)
        end
    end
end

function findAttachedTrees(graph::Graph)
    # Find all maximal attached trees in the graph
    # Output:
    #   ch: Dict{Int, Set{Int}} - children of each node
    #   pr: Dict{Int, Int}      - parent of each node
    #   R: Set{Int}             - roots of maximal attached trees

    ch = Dict{Int, Set{Int}}()
    pr = Dict{Int, Int}()
    R = Set{Int}()
    L = Set{Int}()

    # Initialize: ch empty, pr empty, R = all degree-one vertices
    for node in graph.nodes
        ch[node.node_id] = Set{Int}()
        # degree-one node
        if length(graph.adjacent_edges[node.node_id]) == 1
            push!(R, node.node_id)
            push!(L, node.node_id)
        end
    end

    changed = true
    while changed
        changed = false
        # Collect vertices to process in this round
        to_add = Set{Int}()
        for v in R
            # Find neighbors of v not already in ch[v] or pr[v]
            neighbors = graph.incident_nodes[v]
            # Exclude children (already in ch[v]) and parent (already in pr)
            unprocessed_neighbors = setdiff(neighbors, ch[v])
            if length(unprocessed_neighbors) == 1
                w = first(unprocessed_neighbors)
                pr[v] = w
                push!(ch[w], v)
                push!(to_add, w)
                changed = true
            end
        end
        # Remove processed v from R, add new roots w to R
        for v in R
            if haskey(pr, v)
                delete!(R, v)
            end
        end
        for w in to_add
            push!(R, w)
        end
    end

    function verify(graph, ch, pr, v)
        for w in ch[v]
            @assert pr[w] == v
        end
        if haskey(pr, v)
            prv = Set([pr[v]])
            @assert graph.incident_nodes[v] == union(ch[v], prv)
            @assert isempty(intersect(ch[v], prv))
        else
            for v_ in graph.incident_nodes[v]
                if ! (v_ in ch[v])
                    @assert( length(graph.adjacent_edges[v_]) != 1 )
                end
            end
        end
    end

    for root in R
        traverseTree(graph, ch, pr, root, verify)
    end

    return ch, pr, R, L

end

