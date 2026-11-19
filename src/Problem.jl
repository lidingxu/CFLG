#=========================================================
 Problem Data
=========================================================#

# problem data
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
    EIc::Dict{Int, Set{Tuple{Int, Symbol}}}

    # bounds
    Uv::Vector{Float64}
    dltv::Dict{Tuple{Int, Int}, Float64}
    dlte::Dict{Tuple{Int, Int, Symbol}, Float64}
    Mv::Dict{Tuple{Int, Int}, Float64}
    Me::Dict{Tuple{Int, Int, Symbol}, Float64}
    upper_bd::Int

    # for edge formulation
    bigM_EF::Float64


    # construtor
    function Problem(graph::Graph, dlt::Float64, instance::String = "", c_tol::Float64=1e-6, cr_tol::Float64=1e-6, verbose::Bool=true)
        problem = new(graph, dlt, instance, c_tol, cr_tol, verbose)
        return problem
    end
end



# preprocess problem
function preprocess!(prob::Problem, formulation::FormulationSet)
    if formulation == None
        dtf_graph = breakGraph(prob.graph, prob.dlt, true, true)
        sbd_graph = breakGraph(prob.graph, prob.dlt, true, false)
        org_stat = GraphStat(prob.graph.node_num, prob.graph.edge_num, prob.graph.min_len, prob.graph.max_len, prob.graph.avg_len)
        dtf_stat = GraphStat(dtf_graph.node_num, dtf_graph.edge_num, dtf_graph.min_len, dtf_graph.max_len, dtf_graph.avg_len)
        sdb_stat = GraphStat(sbd_graph.node_num, sbd_graph.edge_num, sbd_graph.min_len, sbd_graph.max_len, sbd_graph.avg_len)
        return [org_stat, dtf_stat, sdb_stat]
    end

    if prob.dlt < prob.graph.max_len
        prob.prob_graph = breakGraph(prob.graph, prob.dlt, true,  mask(formulation, MSK_ISL) )
    else
        prob.prob_graph = prob.graph
    end

    print("problem_graph/original graph:", " node: ", prob.prob_graph.node_num, "/", prob.graph.node_num, " edge: ",
    prob.prob_graph.edge_num, "/", prob.graph.edge_num, " dlt: ", prob.dlt, " break_avg_len: ", prob.prob_graph.avg_len, " break_max_len: ", prob.prob_graph.max_len)

    if formulation == EF # edge formulation needs a simple process
        (prob.bigM_EF, prob.d) = processGraphSimple(prob.prob_graph, prob.dlt)
    else # normal process
        mode = ifelse( mask(formulation, MSK_ISP0), :full, :Partial)
        (prob.Ec, prob.Vc, prob.Ep, prob.Vp, prob.EIp, prob.EIc, prob.d) = processGraph(prob.prob_graph, prob.dlt, :Partial, prob.cr_tol, prob.c_tol)
    end
end



function bounds!(prob::Problem)
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
        Uv[v_id] = min(len*(1+prob.cr_tol), prob.dlt) + prob.c_tol
    end


    # compute dltv, dlte, Mv, Me
    for v_id in graph.node_ids
        for vf_id in prob.Vp[v_id]
            dlen = prob.d[lor(v_id, vf_id)]
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
            # numerical stable version
            dlte[(v_id, ef_id, end_node)] = min(prob.dlt, min(Uv[v_id]+ dlen + elen, prob.dlt) * (1+ prob.cr_tol) + prob.c_tol) # min(Uv[v_id]+ dlen + elen, prob.dlt) + prob.c_tol
            Me[(v_id, ef_id, end_node)] =  max(0, Uv[v_id]+ dlen + elen - prob.dlt) * (1+ prob.cr_tol) + prob.c_tol
            @assert(Me[(v_id, ef_id, end_node)] + (dlte[(v_id, ef_id, end_node)] - dlen - elen) >= Uv[v_id])
        end
    end

    prob.Uv = Uv
    prob.dltv = dltv
    prob.dlte = dlte
    prob.Mv = Mv
    prob.Me = Me
end


"""
    dominanceFilter!(prob::Problem) -> Int

Remove dominated disjunction pairs from `prob.EIp` and return how many were removed.

For a vertex `v` and an edge `f = (u_a, u_b)` of length `l_f` with both `(f, a)` and `(f, b)` in `EIp[v]`,
the pair `(f, b)` is dominated by `(f, a)` if `d(v, u_b) >= d(v, u_a) + l_f` (and symmetrically): then
`tau_vfb(q) = d(v, u_b) + l_f - q >= d(v, u_a) + q = tau_vfa(q)` for every `q` in `[0, l_f]`, so reaching the
facility on `f` through `u_b` never gives more residual coverage than through `u_a`. Both pairs share the same
SOS constraint and the same activation `z <= y_f`, so any solution with `z_vfb = 1` stays feasible with
`z_vfa = 1` instead. By the triangle inequality the test fires exactly when the shortest path from `v` to
`u_b` runs along `f`. For long edges the condition can never hold for two pairs that are both in `EIp[v]`.

Must be called after `processGraph`/`pruning` (which rebuild `EIp`) and before `bounds!`.
"""
function dominanceFilter!(prob::Problem)
    graph = prob.prob_graph
    nrem = 0
    for v_id in graph.node_ids
        kept = Set{Tuple{Int, Symbol}}()
        for (ef_id, end_node) in prob.EIp[v_id]
            other = (ef_id, end_node == :a ? :b : :a)
            ef = graph.edges[ef_id]
            if other in prob.EIp[v_id]
                d_this = prob.d[lor(v_id, ef.nodes[end_node])]
                d_other = prob.d[lor(v_id, ef.nodes[other[2]])]
                if d_this >= d_other + ef.length - 1e-9 * max(1.0, ef.length)
                    nrem += 1
                    continue
                end
            end
            push!(kept, (ef_id, end_node))
        end
        prob.EIp[v_id] = kept
    end
    return nrem
end


function pruning(problem::Problem)
    # Find all maximal attached trees in the graph
    # Output:
    #   ch: Dict{Int, Set{Int}} - children of each node
    #   pr: Dict{Int, Int}      - parent of each node
    #   R: Set{Int}             - roots of maximal attached trees
    graph = problem.prob_graph
    delta = problem.dlt

    ch, pr, R, L = findAttachedTrees(graph)

    VR = Set{Int}()
    function getVR(graph, ch, pr, v)
        push!(VR, v)
    end

    EI = Set{Int}()

    for root in R
        traverseTree(graph, ch, pr, root, getVR)
    end

    function printe(graph, ch, pr, vid)
        for chld in ch[vid]
            eid = graph.node2edge[lor(vid, chld)]
            edge = graph.edges[eid]
            len = edge.length
            print((eid, len, edge.nodes), ",")
        end
        print("\n")
    end
    for root in R
        #print("\n root:", root, ":")
        #traverseTree(graph, ch, pr, root, printe)
    end

    Vbar = Set{Int}(R)
    for node in graph.nodes
        if !(node.node_id in VR)
            push!(Vbar, node.node_id)
        end
    end

    Ebar = Set{Int}()
    for edge in graph.edges
        if edge.nodes[:a] in Vbar && edge.nodes[:b] in Vbar
            push!(Ebar, edge.edge_id)
        end
    end

    r = Dict{Int, Float64}()
    y = Dict{Int, Float64}()
    visited = Set{Int}()
    for v in VR
        r[v] = 0
        y[v] = 0
    end

    for v in L
        push!(visited, v)
    end

    changed = true
    while changed
        changed = false
        for vid in VR
            if !(vid in visited)
                # Check if all children of v are visited
                all_children_visited = all(child -> child in visited, ch[vid])
                rlens = []
                if all_children_visited
                    for chld in ch[vid]
                        eid = graph.node2edge[lor(vid, chld)]
                        edge = graph.edges[eid]
                        len = edge.length
                        glen = len - (delta - r[chld] )
                        rlen = len + r[chld]
                        ye = glen <= 0 ? 0 : floor( glen / (2 * delta)) + 1
                        y[vid] += y[chld] + ye
                        if ye > 0 && edge.etype == :e_normal
                            push!(EI, eid)
                        end
                        push!(rlens, rlen - 2 * delta * ye)
                        #print((chld, y[chld] + ye, rlen - 2*delta * ye))
                    end
                    if !isempty(rlens)
                        idx = findmax(abs, rlens)[2]
                        r[vid] = rlens[idx]
                    else
                        r[vid] = 0.0
                    end
                    push!(visited, vid)
                    changed = true
                end
            end
        end
    end
    #print(r)

    print("Vbar/V, Ebar/E", (length(Vbar), graph.node_num, length(Ebar), graph.edge_num), "\n")
    required = Dict(v => r[v] >= 0 for v in R)
    rho = Dict(v => ( r[v] >= 0 ? r[v] : delta + r[v]) for v in R)
    Ups = Dict(v => y[v] for v in R)

    # pruning statistics (graph before projection = subdivided graph)
    total_len = sum(edge.length for edge in graph.edges; init = 0.0)
    core_len = sum(graph.edges[eid].length for eid in Ebar; init = 0.0)
    RUNINFO["prune_nodes_before"] = graph.node_num
    RUNINFO["prune_edges_before"] = graph.edge_num
    RUNINFO["prune_roots"] = length(R)
    RUNINFO["prune_tree_nodes"] = length(setdiff(VR, R))
    RUNINFO["prune_tree_edges"] = graph.edge_num - length(Ebar)
    RUNINFO["prune_tree_len_share"] = total_len > 0 ? 1 - core_len / total_len : 0.0
    RUNINFO["prune_fixed_facilities"] = sum(values(Ups); init = 0.0)
    RUNINFO["prune_required_roots"] = count(values(required))

    #print((required, "\n", rho, "\n", Ups, "\n", keys(rho)), "\n")

    problem.prob_graph, R_, RE_, R_toR, RtoR_, EtoE_  = projectGraph(graph, VR, R, required, rho, problem.dlt)
    #print(problem.prob_graph.edges, "\n", graph.edges)
    mode = :Partial
    (problem.Ec, problem.Vc, problem.Ep, problem.Vp, problem.EIp, problem.EIc, problem.d) =
        processGraph(problem.prob_graph, problem.dlt, mode, problem.cr_tol, problem.c_tol)
    RUNINFO["prune_nodes_after"] = problem.prob_graph.node_num
    RUNINFO["prune_edges_after"] = problem.prob_graph.edge_num
    Ups_ = Dict(newv => Ups[R_toR[newv]] for newv in R_)
    required_ = Dict(newv => required[R_toR[newv]] for newv in R_)
    return required_, Ups_, R_,  RE_, RtoR_, EtoE_
end