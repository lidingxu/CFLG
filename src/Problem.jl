#=========================================================
 Problem Data
=========================================================#


# problemdata
mutable struct Problem
    # parameters
    graph::Graph # graph
    dlt::Float64 # delta
    instance::String # instance name
    c_tol::Float64 # cover tolerance
    cr_tol::Float64 # relative tolerance
    verbose::Bool # is output model information
    prob_graph::Graph # problem graph: min_len < delta

    # preporcessing data
    d::Dict{Tuple{Int, Int}, Float64} # distance record

    # for vertext formulation
    Ec::Dict{Int, Set{Int}} # edges can completely cover edges
    Vc::Dict{Int, Set{Int}} # nodes can completely cover edges
    Ep::Dict{Int, Set{Int}} # edges can partially cover incident edges of nodes
    Vp::Dict{Int, Set{Int}} # nodes can partially cover incident edges of nodes
    EIp::Dict{Int, Set{Tuple{Int, Symbol}}} # edges, nodes can partially cover incident edges of nodes
   
    # bounds
    Uv::Vector{Float64}
    dltv::Dict{Tuple{Int, Int}, Float64}
    dlte::Dict{Tuple{Int, Int, Symbol}, Float64}
    Mv::Dict{Tuple{Int, Int}, Float64}
    Me::Dict{Tuple{Int, Int, Symbol}, Float64}

    # for edge formulation
    bigM_EF::Float64    


    # construtor
    function Problem(graph::Graph, dlt::Float64, instance::String = "", c_tol::Float64=1e-6, cr_tol::Float64=1e-6, verbose::Bool=true)
        problem = new(graph, dlt, instance, c_tol, cr_tol, verbose)
        return problem
    end
end



# preprocess problem
function preprocess!(prob::Problem, algo::AlgorithmSet)
    if algo == None
        dtf_graph = breakGraph(prob.graph, prob.dlt, true, true)
        sbd_graph = breakGraph(prob.graph, prob.dlt, true, false)
        org_stat = GraphStat(prob.graph.node_num, prob.graph.edge_num, prob.graph.min_len, prob.graph.max_len, prob.graph.avg_len)
        dtf_stat = GraphStat(dtf_graph.node_num, dtf_graph.edge_num, dtf_graph.min_len, dtf_graph.max_len, dtf_graph.avg_len)
        sdb_stat = GraphStat(sbd_graph.node_num, sbd_graph.edge_num, sbd_graph.min_len, sbd_graph.max_len, sbd_graph.avg_len)
        return [org_stat, dtf_stat, sdb_stat]
    end
    if prob.dlt < prob.graph.max_len
        prob.prob_graph = breakGraph(prob.graph, prob.dlt, true, ifelse(algo == RF, true, false))
    else
        prob.prob_graph = prob.graph
    end
    print("problem_graph/original graph:", " node: ", prob.prob_graph.node_num, "/", prob.graph.node_num, " edge: ",
    prob.prob_graph.edge_num, "/", prob.graph.edge_num, " dlt: ", prob.dlt, " break_avg_len: ", prob.prob_graph.avg_len, " break_max_len: ", prob.prob_graph.max_len)
    if algo == EF # edge formulation
        (prob.bigM_EF, prob.d) = processEFGraph(prob.prob_graph, prob.dlt)
        # initilization
    elseif algo == F0
        prob.d = computeDistance(prob.prob_graph)
    else # vertex formulations
        (prob.Ec, prob.Vc, prob.Ep, prob.Vp, prob.EIp, prob.d) = processVFGraph(prob.prob_graph, prob.dlt, :Partial, prob.cr_tol, prob.c_tol)
    end
end



function boundTighten!(prob::Problem)
    # Bound tigntenning
    graph = prob.prob_graph

    Uv = Vector{Float64}(undef, graph.node_num)
    dltv = Dict{Tuple{Int, Int}, Float64}()
    dlte = Dict{Tuple{Int, Int, Symbol}, Float64}()
    Mv = Dict{Tuple{Int, Int}, Float64}()
    Me = Dict{Tuple{Int, Int, Symbol}, Float64}()
    
    # compute Uv
    for v_id in graph.node_ids
        len = typemin(Float64)
        for e_id in graph.adjacent_edges[v_id]
            len = max( ifelse(graph.edges[e_id].etype == :e_long, 2*prob.dlt, graph.edges[e_id].length), len)
        end
        # numerical stable
        #Uv[v_id] = len
        Uv[v_id] = len*(1+prob.cr_tol) + prob.c_tol
    end


    # compute dltv, dlte, Mv, Me
    for v_id in graph.node_ids
        for vf_id in prob.Vp[v_id]
            dlen = prob.d[lor(v_id, vf_id)]
            #dltv[(v_id, vf_id)] = min(Uv[v_id]+ dlen, prob.dlt)
            #Mv[(v_id, vf_id)] = max(0, Uv[v_id]+ dlen-  prob.dlt)
            # numerical stable version
            @assert(prob.dlt - dlen > -1e-6)
            dltv[(v_id, vf_id)] = min(prob.dlt, min(Uv[v_id]+ dlen, prob.dlt) * (1+ prob.cr_tol) + prob.c_tol)
            Mv[(v_id, vf_id)] = max(0, Uv[v_id]+ dlen-  prob.dlt) * (1+ prob.cr_tol) + prob.c_tol
        end
        for (ef_id, end_node) in prob.EIp[v_id]
            ef = graph.edges[ef_id]
            vf_id = ef.nodes[end_node]
            dlen = prob.d[lor(v_id, vf_id)]
            elen =  ifelse(ef.etype == :e_long, 2*prob.dlt, ef.length)
            #dlte[(v_id, ef_id, end_node)] = min(Uv[v_id]+ dlen + elen, prob.dlt)  
            #Me[(v_id, ef_id, end_node)] = max(0, Uv[v_id]+ dlen + elen - prob.dlt)
            # numerical stable version
            dlte[(v_id, ef_id, end_node)] = min(prob.dlt, min(Uv[v_id]+ dlen + elen, prob.dlt) * (1+ prob.cr_tol) + prob.c_tol) # min(Uv[v_id]+ dlen + elen, prob.dlt) + prob.c_tol
            Me[(v_id, ef_id, end_node)] =  max(0, Uv[v_id]+ dlen + elen - prob.dlt) * (1+ prob.cr_tol) + prob.c_tol
            #if(Me[(v_id, ef_id, end_node)] + (dlte[(v_id, ef_id, end_node)] - dlen - elen ) < Uv[v_id])
            #    println(Me[(v_id, ef_id, end_node)] + (dlte[(v_id, ef_id, end_node)] - dlen - elen) - Uv[v_id])
            #end
            @assert(Me[(v_id, ef_id, end_node)] + (dlte[(v_id, ef_id, end_node)] - dlen - elen) >= Uv[v_id])
                
            #if (max(0, Uv[v_id]+ dlen + elen - prob.dlt) * (1+ prob.cr_tol) + prob.c_tol) > prob.dlt
            #    #println(Me[(v_id, ef_id, end_node)] , prob.dlt)
            #end
        end
    end

    prob.Uv = Uv
    prob.dltv = dltv
    prob.dlte = dlte
    prob.Mv = Mv
    prob.Me = Me
end







