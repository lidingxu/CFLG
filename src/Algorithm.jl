#=========================================================
 Algorithm
=========================================================#

function solve!(problem::Problem, solver_name::String, option::Option, algorithm::String)

    # parse string format algorithm setting
    if algorithm == "EF"
        algo = EF
    elseif algorithm == "EFP"
        algo = EFP
    elseif algorithm == "EFP0"
        algo = EFP0
    elseif algorithm == "EFP"
        algo = EFP
    elseif algorithm == "EFPD"
        algo = EFPD
    elseif algorithm == "EFPLg"
        algo = EFPLg
    elseif algorithm == "EFPV"
        algo = EFPV
    elseif algorithm == "EFPV2"
        algo = EFPV2
    elseif algorithm == "EFPL"
        algo = EFPL
    elseif algorithm == "EVF"
        algo = EVF
    elseif algorithm == "EVFP0"
        algo = EVFP0
    elseif algorithm == "EVFP"
        algo = EVFP
    elseif algorithm == "EVFPV"
        algo = EVFPV
    elseif algorithm == "EVFPL"
        algo = EVFPL
    elseif algorithm == "None"
        algo = None
    else
        print("%s not a formulation!!\n", string(algo))
        @assert(false)
    end
    
    # preprocess the problem
    CPUtic()
    graph_stats = preprocess!(problem, algo)

    if algo == None
        return graph_stats
    end

    # bound tightenning
    if mask(algo, MSK_ISP0) && mask(algo, MSK_ISP1) 
        boundTighten!(problem)
    end

    # enumerate cover patterns
    Pi = nothing
    if algo == EFPV || algo == EFPV2
        K = mask(algo, MSK_ISK2) ? 2 : 1 
        print("\nenumerate cover pattern\n")
        if  K != 0
            Pi = enumCoverPatterns(problem, option, solver_name, K, option.time_limit)
            print("\n find cover pattern: ", length(Pi), "\n")
        end
    end


    preprocess_time = CPUtoc()

    cflg = initModel(solver_name, option, option.time_limit - preprocess_time)
    if algo == EF
        stat, sol = solveEF!(problem, algo, cflg)
    elseif algo == EFPL || algo == EVFPL
        #stat, sol = solveRF!(problem, cflg)
        stat, sol = solveFLs!(problem, algo, cflg)
    elseif algo == EFPLg
        stat, sol = solveEFPLg!(problem, algo, cflg)
    elseif algo == EFPD 
        stat, sol = solveEFPD!(problem, algo, cflg)
    else 
        stat, sol = solveFPVs!(problem, algo, cflg, Pi) 
    end

    stat.preprocess_time = preprocess_time
    #verification(problem, sol)

    #println("Optimal value: ", objective_value(cflg))
    println(stat)

    return stat

end

# algo == EF
function solveEF!(problem::Problem, algo::AlgorithmSet, cflg)
    graph = problem.prob_graph
    edges = graph.edges
    dlt = problem.dlt
    e_pairs = Vector{Tuple{Int, Int}}()
    for ef1_id in graph.edge_ids
        for ef2_id in graph.edge_ids
            if ef1_id == ef2_id
                continue
            end
            push!(e_pairs, (ef1_id, ef2_id))
        end
    end

    nodes =[:a,:b]
    enodes = [:a,:b,:zero]

    # variables
    @variables(cflg, begin 
        y[ef_id in graph.edge_ids], Bin # edge facility
        0 <= q[ef_id in graph.edge_ids] <= edges[ef_id].length # edge coordinate variable
        0 <= r[ef_id in graph.edge_ids, v_id in graph.node_ids]# residual cover variable 
        ri[ef_id in graph.edge_ids, node in nodes, v_id in graph.node_ids]# residual cover variable i
        z[e_id in graph.edge_ids, e_pair in e_pairs], Bin # SOS modelling variable
        zm[ef_id in graph.edge_ids, node in enodes, v_id in graph.node_ids], Bin # big M variables
    end) 

    # constraints
    @constraints(cflg, begin 
        [ef_id in graph.edge_ids, node in nodes, v_id in graph.node_ids], dlt - ifelse(node == :a, q[ef_id], edges[ef_id].length - q[ef_id]) - problem.d[lor(edges[ef_id].nodes[node], v_id)] == ri[ef_id, node, v_id] # cover range
        [e_id in graph.edge_ids, e_pair in e_pairs], r[e_pair[1], edges[e_id].nodes[:a]] + r[e_pair[2], edges[e_id].nodes[:b]] >= z[e_id, e_pair] * (edges[e_id].length * (1+problem.cr_tol) + problem.c_tol) # cover condition
        [e_id in graph.edge_ids], sum(z[e_id, e_pair] for e_pair in e_pairs) +  y[e_id] == 1 # SOS 1
        [e_id in graph.edge_ids,  e_pair in e_pairs], z[e_id, e_pair] <= y[e_pair[1]] # open contion
        [e_id in graph.edge_ids,  e_pair in e_pairs], z[e_id, e_pair] <= y[e_pair[2]] # open contion
    end)


    @constraints(cflg, begin 
        [ef_id in graph.edge_ids, node in nodes, v_id in graph.node_ids], r[ef_id, v_id] >= ri[ef_id, node, v_id] # big M, lb
        [ef_id in graph.edge_ids, node in nodes, v_id in graph.node_ids], r[ef_id, v_id] - (1 - zm[ef_id, node, v_id]) * problem.bigM_EF <= ri[ef_id, node, v_id]# big M, ub
        [ef_id in graph.edge_ids, v_id in graph.node_ids], r[ef_id, v_id] - (1 - zm[ef_id, :zero, v_id]) * problem.bigM_EF <= 0# big M, ub
        [ef_id in graph.edge_ids, v_id in graph.node_ids], sum(zm[ef_id, node, v_id] for node in enodes) == 1 # big M, SOS
    end) 

    # objective
    @objective(cflg, Min, sum(y[ef_id] for ef_id in graph.edge_ids))
    println("\n model loaded\n") 
    optimize!(cflg)

    stat = Stat();
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg).optimizer.model
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end

    stat.algo = algo
    stat.instance = problem.instance
    if has_values(cflg)
        stat.sol_val = objective_value(cflg)
        for ef_id in graph.edge_ids
            if value(y[ef_id]) >= 0.5
                push!(sol, (:e, ef_id,  value(q[ef_id])))
            end
        end
    end

    return stat, sol
end


# 
function solveFPVs!(problem::Problem, algo::AlgorithmSet, cflg, Pi)
    # get graph
    graph = problem.prob_graph

    # get delimitation
    if mask(algo, MSK_ISP0)
        EIp = problem.EIp
        Vp = problem.Vp

    elseif algo == EVF
        EIp = Dict{Int, Set{Tuple{Int, Symbol}}}()
        Vp = Dict{Int, Set{Int}}()
        Eid = Set{Tuple{Int, Symbol}}()
        Vid= Set{Int}(graph.node_ids)
        for e_id in graph.edge_ids
            push!(Eid, Tuple{Int, Symbol}([e_id, :a]))
            push!(Eid, Tuple{Int, Symbol}([e_id, :b]))
        end
        for v_id in graph.node_ids
            Vp[v_id] = Vid
            EIp[v_id] = Eid
        end
    end

    # get bounds
    is_trivial_bound = false
    if !mask(algo, MSK_ISP1)
        is_trivial_bound = true
        bigM = problem.dlt * (problem.cr_tol + 1) + problem.c_tol
    end

    usev = !mask(algo, MSK_ISE)

    # basic variables
    @variables(cflg, begin 
        ye[ef_id in graph.edge_ids], Bin # edge facility
        x[v_id in graph.node_ids], Bin #  node residual indictor cover
        w[e_id in graph.edge_ids], Bin # complete cover indicator variable
        0 <= q[ef_id in graph.edge_ids] <= graph.edges[ef_id].length # edge coordinate variable
        0 <= rv[v_id in graph.node_ids] <= (is_trivial_bound ? bigM : problem.Uv[v_id]) # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # big-M modelling variable on edges
    end) 

    # vertex facility related variables
    if usev
        print("vertex var\n")
        @variables(cflg, begin 
            yv[vf_id in graph.node_ids], Bin # node facility
            zv[v_id in graph.node_ids, vf_id in Vp[v_id]], Bin # big-M modelling variable on nodes
        end) 
    end

    # basic constraints
    @constraints(cflg, begin 
        [e_id in graph.edge_ids, ef_id in problem.Ec[e_id]], w[e_id] >= ye[ef_id] # complete cover by edges: open
        [v_id in graph.node_ids], x[v_id] >= 1- sum( (1 - w[e_id]) for e_id in graph.adjacent_edges[v_id])  # ajdacent non covered 2
        [v_id in graph.node_ids, e_id in graph.adjacent_edges[v_id]], x[v_id] <= w[e_id]  # ajdacent non covered 2
        [e_id in graph.edge_ids], (graph.edges[e_id].length *  (1+problem.cr_tol) + problem.c_tol) * (1 - w[e_id]) <= rv[graph.edges[e_id].nodes[:a]] + rv[graph.edges[e_id].nodes[:b]] # jointly complete cover condition
        [e_id in graph.edge_ids], q[e_id] <= graph.edges[e_id].length * ye[e_id]# redundant bound
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] <= ye[efi[1]] # edge activated constraint
        [v_id in graph.node_ids], rv[v_id] <=  (is_trivial_bound ? bigM : problem.Uv[v_id] ) * (1 - x[v_id]) # big M on x
        [v_id in graph.node_ids, efi in EIp[v_id]], rv[v_id] <= ( is_trivial_bound ? ( bigM +  graph.edges[efi[1]].length*  (1+problem.cr_tol) ) :
         problem.Me[(v_id, efi[1], efi[2])] ) * (1 - ze[v_id, efi]) + ( is_trivial_bound ? problem.dlt : problem.dlte[(v_id, efi[1], efi[2])] ) -
         ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1]], (graph.edges[efi[1]].length - q[efi[1]]) ))  # big M on edges
    end) 

    # vertex facility related constraints
    if usev
        @constraints(cflg, begin 
            [e_id in graph.edge_ids], w[e_id] <= sum(yv[vf_id] for vf_id in problem.Vc[e_id]) + sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
            [e_id in graph.edge_ids, vf_id in problem.Vc[e_id]], w[e_id] >= yv[vf_id] # complete cover by nodes: open
            [v_id in graph.node_ids, vf_id in Vp[v_id]], zv[v_id, vf_id] <= yv[vf_id] # node activated constraint
            [v_id in graph.node_ids, vf_id in Vp[v_id]], rv[v_id] <= (is_trivial_bound ? bigM : problem.Mv[(v_id, vf_id)] ) * (1 - zv[v_id, vf_id]) + (is_trivial_bound ? problem.dlt : problem.dltv[(v_id, vf_id)]) - problem.d[lor(v_id,vf_id)] # big M on nodes
            [v_id in graph.node_ids], x[v_id] + sum(zv[v_id, vf_id] for vf_id in Vp[v_id]) + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
        end) 
    else
        @constraints(cflg, begin 
            [e_id in graph.edge_ids], w[e_id] <= sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
            [v_id in graph.node_ids], x[v_id] + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
        end) 
    end

    root_user_cut = true
    # valid inequalities
    if !mask(algo, MSK_ISE) && mask(algo, MSK_ISV)
        # valid inequalities
        for v_id in graph.node_ids
            # leave fixing
            if length(graph.adjacent_edges[v_id]) == 1
                # constraints
                @constraints(cflg, begin
                    yv[v_id] == 0 # fixing
                    [ef_id in graph.adjacent_edges[v_id]], ye[ef_id] == 0 # fixing
                end)
            else # adjacent edges
                #continue
                @constraint(cflg, 
                    sum(ye[ef_id] for ef_id in graph.adjacent_edges[v_id])  + yv[v_id]<= 1)
            end
        end
    elseif mask(algo, MSK_ISE) && mask(algo, MSK_ISV)
        if Pi !== nothing
            if root_user_cut
                @constraints(cflg, begin
                    [pi in Pi],  sum( (1 - ze[tple[1], tple[2]]) for tple in pi) >= 1
                end)
            end
        end
    end
    #MOI.set(cflg, MOI.UserCutCallback(), user_cut_callback)

    # objective
    @objective(cflg, Min, (usev ? sum(yv[vf_id] for vf_id in graph.node_ids) : 0) +  sum(ye[ef_id] for ef_id in graph.edge_ids))

    println("\n model loaded\n")  
    optimize!(cflg)
    stat = Stat();
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg).optimizer.model
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    stat.algo = algo
    stat.instance = problem.instance
    if has_values(cflg)
        stat.sol_val = objective_value(cflg)
        for ef_id in graph.edge_ids
            if value(ye[ef_id]) >= 0.5
                push!(sol, (:e, ef_id,  value(q[ef_id])))
            end
        end
        if usev
            for v_id in graph.node_ids
                if value(yv[v_id]) >= 0.5
                    push!(sol, (:v, v_id,  0))
                end
            end
        end
    end
    return stat, sol
end


# long edge vertex models
function solveFLs!(problem::Problem, algo::AlgorithmSet, cflg)

    # get graph
    graph = problem.prob_graph

    # find long edges
    edges = graph.edges
    EIp = problem.EIp
    Vp = problem.Vp
    Elong = Set{Int}()
    tail_len = Dict{Int, Float64}()
    fac_num = Dict{Int, Float64}()
    for edge in edges
        if edge.etype == :e_long
            push!(Elong, edge.edge_id)
            tail_len[edge.edge_id] = edge.length - (floor(edge.length / (2 * problem.dlt))) * 2 * problem.dlt 
            fac_num[edge.edge_id] = ceil(edge.length / (2 * problem.dlt))
        end
    end

    usev = !mask(algo, MSK_ISE)

    # variables
    @variables(cflg, begin 
        ye[ef_id in graph.edge_ids], Bin # edge facility
        x[v_id in graph.node_ids], Bin #  node residual indictor cover
        w[e_id in graph.edge_ids], Bin # complete cover indicator variable
        u[e_id in Elong], Bin # indictor variable for phase transition
        0 <= q[ef_id in graph.edge_ids] <= ifelse(edges[ef_id].etype == :e_long, 2*problem.dlt, graph.edges[ef_id].length) # edge coordinate variable
        0 <= rv[v_id in graph.node_ids] <= problem.Uv[v_id] # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # big-M modelling variable on edges
    end) 

    # vertex facility related variables
    if usev
        @variables(cflg, begin 
            yv[vf_id in graph.node_ids], Bin # node facility
            zv[v_id in graph.node_ids, vf_id in Vp[v_id]], Bin # big-M modelling variable on nodes
        end) 
    end


    E = Set{Int}(graph.edge_ids)
    Enormal = setdiff(E, Elong)

    # constraints
    @constraints(cflg, begin 
        [e_id in Enormal, ef_id in problem.Ec[e_id]], w[e_id] >= ye[ef_id] # complete cover by edges: open
        [v_id in graph.node_ids], x[v_id] >= 1- sum( (1 - w[e_id]) for e_id in graph.adjacent_edges[v_id])  # ajdacent non covered 2
        [v_id in graph.node_ids, e_id in graph.adjacent_edges[v_id]], x[v_id] <= w[e_id]  # ajdacent non covered 2
        [e_id in Enormal], (graph.edges[e_id].length *  (1+problem.cr_tol) + problem.c_tol) * (1 - w[e_id]) <= rv[graph.edges[e_id].nodes[:a]] + rv[graph.edges[e_id].nodes[:b]] # jointly complete cover condition
        [e_id in Enormal], q[e_id] <= graph.edges[e_id].length * ye[e_id] # redundant bound
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] <= ye[efi[1]] # edge activated constraint
        [v_id in graph.node_ids], rv[v_id] <= problem.Uv[v_id] * (1 - x[v_id]) # big M on x
        [v_id in graph.node_ids, efi in EIp[v_id]], rv[v_id] <= problem.Me[(v_id, efi[1], efi[2])] * (1 - ze[v_id, efi]) + problem.dlte[(v_id, efi[1], efi[2])] -
         ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1]], edges[efi[1]].etype == :e_normal ? graph.edges[efi[1]].length - q[efi[1]] : 2*problem.dlt * u[efi[1]]+ tail_len[efi[1]]  - q[efi[1]]  )) # # big M on edges
    end) 

    if usev
        @constraints(cflg, begin 
            [e_id in Enormal, vf_id in problem.Vc[e_id]], w[e_id] >= yv[vf_id] # complete cover by nodes: open
            [e_id in Enormal], w[e_id] <= sum(yv[vf_id] for vf_id in problem.Vc[e_id]) + sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
            [ef_id in Enormal, node in [:a,:b]], yv[graph.edges[ef_id].nodes[node]] + ye[ef_id] <= 1 # facility is efither in interior or at end nodes
            [v_id in graph.node_ids], x[v_id] +  sum(zv[v_id, vf_id] for vf_id in Vp[v_id])  + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
            [v_id in graph.node_ids, vf_id in Vp[v_id]], zv[v_id, vf_id] <= yv[vf_id] # node activated constraint
            [v_id in graph.node_ids, vf_id in Vp[v_id]], rv[v_id] <=  problem.Mv[(v_id, vf_id)]  * (1 - zv[v_id, vf_id]) + problem.dltv[(v_id, vf_id)] - problem.d[lor(v_id,vf_id)] # big M on nodes
            [ef_id in Elong], yv[edges[ef_id].nodes[:a]] == 0  #valid inequalities fixing
        end) 
    else
        @constraints(cflg, begin 
            [e_id in Enormal], w[e_id] <= sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
            [v_id in graph.node_ids], x[v_id] + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
        end) 
    end

    # long edge constraints
    @constraints(cflg,begin
        [e_id in Elong], ye[e_id] == 1 # fixing long edges y
        [e_id in Elong], w[e_id] == 0 # fixing long edges w
        [ef_id in Elong], q[ef_id] <=  tail_len[ef_id] * (1 - u[ef_id]) + 2 * problem.dlt *  u[ef_id] # phase transition
        [ef_id in Elong], q[ef_id] >=  tail_len[ef_id] *  u[ef_id] # phase transition
        [ef_id in Elong], rv[edges[ef_id].nodes[:a]] + problem.dlt >=   q[ef_id]  # head cover
        [ef_id in Elong], rv[edges[ef_id].nodes[:b]] +  q[ef_id] - (2*  u[ef_id] - 1) * problem.dlt  >=  tail_len[ef_id]  # tail cover
    end) # to do valid inequalities
    
    # objective
    if usev
        @objective(cflg, Min,  sum(yv[vf_id] for vf_id in graph.node_ids) + sum((edges[ef_id].etype == :e_long ? fac_num[ef_id] - u[ef_id] : ye[ef_id]) for ef_id in graph.edge_ids))
    else
        @objective(cflg, Min,  sum((edges[ef_id].etype == :e_long ? fac_num[ef_id] - u[ef_id] : ye[ef_id]) for ef_id in graph.edge_ids))
    end

    println("\n model loaded\n")
    optimize!(cflg)
    stat = Stat();
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg).optimizer.model
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end

    stat.algo = algo
    stat.instance = problem.instance
    if has_values(cflg)
        stat.sol_val = objective_value(cflg)
        for ef_id in graph.edge_ids
            if value(ye[ef_id]) >= 0.5
                push!(sol, (:e, ef_id,  value(q[ef_id])))
            end
            if edges[ef_id].etype == :e_long 
                #println("\n", fac_num[ef_id])
                if value(u[ef_id]) >= 0.5
                    push!(sol, (:e_long, ef_id,  value(u[ef_id])))
                end
            end
        end
        if usev
            for v_id in graph.node_ids
                if value(yv[v_id]) >= 0.5
                    push!(sol, (:v, v_id,  0))
                end
            end
        end
    end
    #println(graph.edges)
    return stat, sol
end


# EF disjunctive formulation
function solveEFPD!(problem::Problem, algo::AlgorithmSet, cflg)
    graph = problem.prob_graph
    EIp = problem.EIp
    Vp = problem.Vp
    Ep = problem.Ep

    # bound tightenning
    isboundtight = false
    if isboundtight
        Lbs = Dict{Tuple{Int, Tuple{Int, Symbol}}, Float64}()
        Ubs = Dict{Tuple{Int, Tuple{Int, Symbol}}, Float64}()
        for v_id in graph.node_ids
            for efi in EIp[v_id]
                exceed_len = min( max(problem.dlte[(v_id, efi[1], efi[2])] - problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])], 0), graph.edges[efi[1]].length)
                exceed_len = isboundtight * exceed_len
                Ubs[(v_id, efi)] = ifelse(efi[2] == :a, exceed_len, graph.edges[efi[1]].length )  * (1+problem.cr_tol)
                Lbs[(v_id, efi)] = ifelse(efi[2] == :a, 0, (graph.edges[efi[1]].length - exceed_len) * (1 - problem.cr_tol) )
            end
        end
    end

    # variables
    @variables(cflg, begin 
        ye[ef_id in graph.edge_ids], Bin # edge facility
        x[v_id in graph.node_ids], Bin # node residual indictor cover
        w[e_id in graph.edge_ids], Bin # complete cover indicator variable
        0 <= q[ef_id in graph.edge_ids] <= graph.edges[ef_id].length # edge coordinate variable
        0 <= qve[v_id in graph.node_ids, ef_id in Ep[v_id]] <= graph.edges[ef_id].length  # disjunctive edge coordinate variable on nodes
        0<= qvei[v_id in graph.node_ids, efi in EIp[v_id]]
        0 <= rv[v_id in graph.node_ids] <= problem.Uv[v_id] # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # disjunctive indicator variable on edges
    end) 

    # constraints
    @constraints(cflg, begin 
        [e_id in graph.edge_ids, ef_id in problem.Ec[e_id]], w[e_id] >= ye[ef_id] # complete cover by edges: open
        [e_id in graph.edge_ids], w[e_id] <=  sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
        [v_id in graph.node_ids], x[v_id] >= 1- sum( (1 - w[e_id]) for e_id in graph.adjacent_edges[v_id])  # ajdacent non covered 2
        [v_id in graph.node_ids, e_id in graph.adjacent_edges[v_id]], x[v_id] <= w[e_id]  # ajdacent non covered 2
        [e_id in graph.edge_ids], (graph.edges[e_id].length *  (1+problem.cr_tol) + problem.c_tol) * (1 - w[e_id]) <= rv[graph.edges[e_id].nodes[:a]] + rv[graph.edges[e_id].nodes[:b]] # jointly complete cover condition
        [e_id in graph.edge_ids], q[e_id] <= graph.edges[e_id].length * ye[e_id] # redundant bound
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] <= ye[efi[1]] # edge activated constraint        
        [v_id in graph.node_ids], x[v_id] +  sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # disjunctive SOS-1 constraint
        [v_id in graph.node_ids], rv[v_id]  <=  sum( ( problem.dlte[(v_id, efi[1], efi[2])] - problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] - ifelse( efi[2] == :a , 0, (graph.edges[efi[1]].length) ) ) * ze[v_id, efi] + ifelse(efi[2] == :a , -qvei[v_id, efi], qvei[v_id, efi] ) for efi in EIp[v_id]) # disjunctive residual cover aggreagation constraint
        [v_id in graph.node_ids, ef_id in Ep[v_id]], q[ef_id] ==   qve[v_id, ef_id] + sum( ifelse(efi[1] == ef_id, qvei[v_id, efi], 0) for efi in EIp[v_id] )  # disjunctive residual q aggregation constraint
        [v_id in graph.node_ids, efi in EIp[v_id]], qvei[v_id, efi] <=  graph.edges[efi[1]].length * (1+problem.cr_tol)  * ze[v_id, efi] # disjunctive upper bound for q on edge
        #[v_id in graph.node_ids, efi in EIp[v_id]], Lbs[(v_id, efi)] * ze[v_id, efi]  <= qvei[v_id, efi]  # disjunctive lower bound for q on edge
        [v_id in graph.node_ids, ef_id in Ep[v_id]], qve[v_id, ef_id] <= graph.edges[ef_id].length * (1+problem.cr_tol) * (1 -   sum( ifelse( efi[1] == ef_id, ze[v_id, efi], 0) for efi in EIp[v_id])) # disjunctive upper bound for q on v
    end) 



    # objective
    @objective(cflg, Min,  sum(ye[ef_id] for ef_id in graph.edge_ids))

    println("\n model loaded\n")  
    optimize!(cflg)
    stat = Stat();
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg).optimizer.model
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    stat.algo = algo

    stat.instance = problem.instance
    if has_values(cflg)
        stat.sol_val = objective_value(cflg)
        for ef_id in graph.edge_ids
            if value(ye[ef_id]) >= 0.5
                push!(sol, (:e, ef_id,  value(q[ef_id])))
            end
        end
    end
    return stat, sol
end


# EF logical formulation indicator
function solveEFPLg!(problem::Problem, algo::AlgorithmSet, cflg)
    print("Lg\n")
    graph = problem.prob_graph
    EIp = problem.EIp
    Vp = problem.Vp
    Ep = problem.Ep

    # variables
    @variables(cflg, begin 
        ye[ef_id in graph.edge_ids], Bin # edge facility
        x[v_id in graph.node_ids], Bin # node residual indictor cover
        w[e_id in graph.edge_ids], Bin # complete cover indicator variable
        0 <= q[ef_id in graph.edge_ids] <= graph.edges[ef_id].length # edge coordinate variable
        0 <= rv[v_id in graph.node_ids] <= problem.Uv[v_id] # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # indicator modelling variable on edges
    end) 


    #problem.Me[(v_id, efi[1], efi[2])]  +  problem.dlte[(v_id, efi[1], efi[2])]  -
    #( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1]], graph.edges[efi[1]].length - q[efi[1]]))
    # constraints
    @constraints(cflg, begin 
        [e_id in graph.edge_ids, ef_id in problem.Ec[e_id]], w[e_id] >= ye[ef_id] # complete cover by edges: open
        [e_id in graph.edge_ids], w[e_id] <=  sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
        [v_id in graph.node_ids], x[v_id] >= 1- sum( (1 - w[e_id]) for e_id in graph.adjacent_edges[v_id])  # ajdacent non covered 2
        [v_id in graph.node_ids, e_id in graph.adjacent_edges[v_id]], x[v_id] <= w[e_id]  # ajdacent non covered 2
        [e_id in graph.edge_ids], (graph.edges[e_id].length *  (1+problem.cr_tol) + problem.c_tol) * (1 - w[e_id]) <= rv[graph.edges[e_id].nodes[:a]] + rv[graph.edges[e_id].nodes[:b]] # jointly complete cover condition
        [e_id in graph.edge_ids], q[e_id] <= graph.edges[e_id].length * ye[e_id] # redundant bound   
        [v_id in graph.node_ids], x[v_id] + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] <= ye[efi[1]] # edge activated constraint
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] => { rv[v_id] <=  problem.dlte[(v_id, efi[1], efi[2])]  -
        ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1]], graph.edges[efi[1]].length - q[efi[1]])) }
    end) 

    # objective
    @objective(cflg, Min,  sum(ye[ef_id] for ef_id in graph.edge_ids))

    println("\n model loaded\n")  
    optimize!(cflg)
    stat = Stat();
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg).optimizer.model
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    stat.algo = algo

    stat.instance = problem.instance
    if has_values(cflg)
        stat.sol_val = objective_value(cflg)
        for ef_id in graph.edge_ids
            if value(ye[ef_id]) >= 0.5
                push!(sol, (:e, ef_id,  value(q[ef_id])))
            end
        end
    end
    return stat, sol
end


# verify the graph is covered, not for long edge formulation
function verification(problem::Problem, sol::Vector{Tuple{Symbol,Int, Float64}})
    graph = problem.prob_graph
    edges = graph.edges
    dlt = problem.dlt

    d = Dict{Tuple{Int, Int}, Float64}() # distance record

    # initilization
    for v_id in graph.node_ids
        for v_id_ in graph.node_ids
            if v_id <= v_id_
                d[lor(v_id,v_id_)] = typemax(Float64)
            end
        end
    end

    # only distance informulation is used
    for v_id in graph.node_ids
        nodeCover(graph, v_id, dlt, d, :full) # run node cover in full mode
    end

    range = Vector{Dict{Symbol, Float64}}()
    rangect = Vector{Dict{Symbol, Tuple{Symbol,Int, Symbol}}}()
    is_cover = Vector{Bool}()

    for e_id in graph.edge_ids  
        push!(range,  Dict{Symbol, Float64}(:a=>0., :b=>0.))
        push!(rangect,  Dict{Symbol, Tuple{Symbol,Int, Symbol}}(:a=>(:v,-1,:zero), :b=>(:v,-1, :zero)))
        push!(is_cover, false)
    end


    for fac in sol
        if fac[1] == :e
            ef_id = fac[2]
            q = fac[3]
            len = edges[ef_id].length
            is_cover[ef_id] = true
            for e_id in graph.edge_ids
                if is_cover[e_id]
                    continue
                end
                a_symb = :zero
                if d[lor(edges[ef_id].nodes[:a], edges[e_id].nodes[:a])] + q <  d[lor(edges[ef_id].nodes[:b], edges[e_id].nodes[:a])] + len -q
                    a_symb = :a
                else
                    a_symb = :b
                end
                b_symb = :zero
                if d[lor(edges[ef_id].nodes[:a], edges[e_id].nodes[:b])] + q < d[lor(edges[ef_id].nodes[:b], edges[e_id].nodes[:b])] + len - q
                    b_symb = :a
                else
                    b_symb = :b
                end
                range_a = dlt - min(d[lor(edges[ef_id].nodes[:a], edges[e_id].nodes[:a])] + q, d[lor(edges[ef_id].nodes[:b], edges[e_id].nodes[:a])] + len -q)
                range_b = dlt - min(d[lor(edges[ef_id].nodes[:a], edges[e_id].nodes[:b])] + q, d[lor(edges[ef_id].nodes[:b], edges[e_id].nodes[:b])] + len - q)
                if range_a >= range[e_id][:a]
                    rangect[e_id][:a] = (:e, ef_id, a_symb) 
                end
                if range_b >= range[e_id][:b]
                    rangect[e_id][:b] = (:e, ef_id, b_symb) 
                end
                range[e_id][:a] = max(range_a, range[e_id][:a])
                range[e_id][:b] = max(range_b, range[e_id][:b])
                if range[e_id][:a] + range[e_id][:b] >= edges[e_id].length
                    is_cover[e_id] = true
                end
            end
        elseif fac[1] == :v
            v_id = fac[2]
            for e_id in graph.edge_ids
                if is_cover[e_id]
                    continue
                end
                range_a = dlt - d[lor(v_id, edges[e_id].nodes[:a])]
                range_b = dlt - d[lor(v_id, edges[e_id].nodes[:b])]
                if range_a >= range[e_id][:a]
                    rangect[e_id][:a] = (:v, v_id, :zero) 
                end
                if range_b >= range[e_id][:b]
                    rangect[e_id][:b] = (:v, v_id, :zero) 
                end
                range[e_id][:a] = max(range_a, range[e_id][:a])
                range[e_id][:b] = max(range_b, range[e_id][:b])
                if range[e_id][:a] + range[e_id][:b] >= edges[e_id].length
                    is_cover[e_id] = true
                end
            end        
        end
    end

    is_covered = true

    for e_id in graph.edge_ids
        if  !is_cover[e_id]
            is_covered = false
            break
        end
    end

    cover_info = []

    for e_id in graph.edge_ids
        push!(cover_info, (e_id, is_cover[e_id], :a, range[e_id][:a], rangect[e_id][:a], range[e_id][:b], :b, rangect[e_id][:b], edges[e_id].length, range[e_id][:a] + range[e_id][:b] - edges[e_id].length))
        if !is_cover[e_id] && problem.verbose 
            print("\n",cover_info[end])
        end
    end
    if problem.verbose
        print("\n:",sol)
    end

    println("\n is covered:", is_covered, "\n")

end


