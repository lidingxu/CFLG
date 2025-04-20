#=========================================================
 Algorithm
=========================================================#

function solve!(problem::Problem, solver_name::String, option::Option, formulation::String)
    # parse string format formulation setting
    if formulation == "EF"
        formulation = EF
    elseif formulation == "EFP"
        formulation = EFP
    elseif formulation == "EFP0"
        formulation = EFP0
    elseif formulation == "EFP"
        formulation = EFP
    elseif formulation == "EFPB"
        formulation = EFPB
    elseif formulation == "EFPD"
        formulation = EFPD
    elseif formulation == "EFPDB"
        formulation = EFPDB
    elseif formulation == "EFPI"
        formulation = EFPI
    elseif formulation == "EFPV"
        formulation = EFPV
    elseif formulation == "EFPV2"
        formulation = EFPV2
    elseif formulation == "EFPL"
        formulation = EFPL
    elseif formulation == "EVF"
        formulation = EVF
    elseif formulation == "EVFP0"
        formulation = EVFP0
    elseif formulation == "EVFP"
        formulation = EVFP
    elseif formulation == "EVFPV"
        formulation = EVFPV
    elseif formulation == "LEVFP"
        formulation = LEVFP
    elseif formulation == "LEFP"
        formulation = LEFP
    elseif formulation == "LEFP2"
        formulation = LEFP2
    elseif formulation == "LEFPV"
        formulation = LEFPV
    elseif formulation == "LEFPV2"
        formulation = LEFPV2
    elseif formulation == "LEFPB"
        formulation = LEFPB
    elseif formulation == "LEFPI"
        formulation = LEFPI
    elseif formulation == "LEFPD"
        formulation = LEFPD
    elseif formulation == "LEFPDB"
        formulation = LEFPDB
    elseif formulation == "None"
        formulation = None
    else
        print("%s not a formulation!!\n", string(formulation))
        @assert(false)
    end

    # preprocess the problem
    CPUtic()
    graph_stats = preprocess!(problem, formulation)

    if formulation == None
        return graph_stats
    end

    # bound tightenning
    if mask(formulation, MSK_ISP0) && mask(formulation, MSK_ISP1)
        boundTighten!(problem)
    end

    # enumerate cover patterns
    Pi = nothing
    if formulation == EFPV || formulation == EFPV2
        K = mask(formulation, MSK_ISV2) ? 2 : 1
        print("\n enumerate cover pattern\n")
        if  K != 0
            Pi = enumCoverPatterns(problem, option, solver_name, K, option.time_limit)
            print("\n find cover pattern: ", length(Pi), "\n")
        end
    end

    preprocess_time = CPUtoc()

    # initialize a JuMP with the option
    cflg = initModel(solver_name, option, option.time_limit - preprocess_time)

    # solve the formulation
    is_edge = mask(formulation, MSK_ISE)
    is_preprocess = mask(formulation, MSK_ISP0) || mask(formulation, MSK_ISP1)
    is_benders = mask(formulation, MSK_ISB)
    is_long = mask(formulation, MSK_ISL)

    if is_edge
        if is_preprocess
            stat, sol = solveEFP!(problem, formulation, cflg)
        else
            stat, sol = solveEF!(problem, formulation, cflg)
        end
    elseif formulation == LEVFP
        stat, sol = solveLEVF!(problem, formulation, cflg)
    else
        stat, sol = solveFPVs!(problem, formulation, cflg, Pi, is_benders)
    end

    #if formulation == EF
    #    stat, sol = solveEF!(problem, formulation, cflg)
    #elseif formulation == LEVFP
    #    stat, sol = solveLEVF!(problem, formulation, cflg)
    #elseif formulation == EFPI
    #    stat, sol = solveEFPI!(problem, formulation, cflg)
    #elseif formulation == EFPD || formulation == EFPDB
    #    stat, sol = solveEFPD!(problem, formulation, cflg, is_benders)
    #else
    #    stat, sol = solveFPVs!(problem, formulation, cflg, Pi, is_benders)
    #end

    stat.preprocess_time = preprocess_time
    println(stat)

    return stat

end

# edge model disjunctive formulation (preprocessed)
function solveEFP!(problem::Problem, formulation::FormulationSet, cflg)
    is_indicator = mask(formulation, MSK_ISI)
    is_lifted = mask(formulation, MSK_ISD)
    is_benders = mask(formulation, MSK_ISB)
    is_cuts = mask(formulation, MSK_ISV)
    is_morecuts = mask(formulation, MSK_ISV2)
    longedge1 = ! mask(formulation, MSK_ISL2)
    print("\n more cuts:", is_morecuts, "\n")
    graph = problem.prob_graph
    EIp = problem.EIp
    dlt = problem.dlt

    Es = [ef_id for ef_id in graph.edge_ids if graph.edges[ef_id].etype == :e_normal]
    El = [ef_id for ef_id in graph.edge_ids if graph.edges[ef_id].etype == :e_long]
    Vs = [v_id for v_id in graph.node_ids if all(graph.edges[ef_id].etype == :e_normal for ef_id in graph.adjacent_edges[v_id])]
    Vl = [v_id for v_id in graph.node_ids if any(graph.edges[ef_id].etype == :e_long for ef_id in graph.adjacent_edges[v_id])]
    #Vscomp = setdiff(graph.node_ids, Vs)
    ab = (:a, :b)
    abc = (:a, :b, :c)
    abcd = (:a, :b, :c, :d)
    print("formulation:", formulation, length(El), " ", length(Es), " ", length(Vs), "\n")
    # variables
    @variables(cflg, begin
        ye[ef_id in graph.edge_ids], Bin # edge facility
        x[v_id in graph.node_ids], Bin # node residual indictor cover
        w[e_id in Es], Bin # complete cover indicator variable
        0 <= q[ef_id in graph.edge_ids, i in ab] <= graph.edges[ef_id].length # edge coordinate variable
        0 <= rv[v_id in graph.node_ids] <= problem.Uv[v_id] # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # disjunctive indicator variable on edges
    end)

    if longedge1
        @variables(cflg, begin
            yei[ef_id in El, i in ab], Int # other edge facility of long edge
            0 <= qq[ef_id in graph.edge_ids] <= 2 * dlt # edge coordinate variable
        end)
    else
        @variables(cflg, begin
            yei[ef_id in El, i in abcd], Bin # other edge facility of long edge
            0 <= qq[ef_id in graph.edge_ids] <= 2 * dlt
        end)
    end

    qLBs = Dict{Tuple{Int, Symbol}, Float64}()
    qUBs = Dict{Tuple{Int, Symbol}, Float64}()
    for ef_id in graph.edge_ids
        if ef_id in Es
            for i in ab
                qLBs[(ef_id, i)] = 0.0
                qUBs[(ef_id, i)] = graph.edges[ef_id].length
            end
        else
            qLBs[(ef_id, :a)] = 0.0
            qLBs[(ef_id, :b)] = 2 * dlt
            qUBs[(ef_id, :a)] = graph.edges[ef_id].length
            qUBs[(ef_id, :b)] = graph.edges[ef_id].length
        end
    end

    if longedge1
        @constraints(cflg, begin
            [ef_id in El], 0 <= yei[ef_id, :b] # long edge facility
            [ef_id in El], yei[ef_id, :b] <= 1 # long edge facility
            [ef_id in El], -1 + floor(graph.edges[ef_id].length / (2 * dlt) )  <= yei[ef_id, :a]  # long edge facility
            [ef_id in El], yei[ef_id, :a] <= floor(graph.edges[ef_id].length / (2 * dlt) )  # long edge facility
            [ef_id in El], q[ef_id, :b] == q[ef_id, :a] + 2 * dlt * yei[ef_id, :a] + qq[ef_id] # long edge coordinate relation
            [ef_id in El], q[ef_id, :b] >= graph.edges[ef_id].length * yei[ef_id, :b] # long edge coordinate lower bound
            [ef_id in El], qq[ef_id] <= 2 * dlt * yei[ef_id, :b] # long edge coordinate upper bound
        end)
    else
        # x_1 x_2 = 1 => y_12 = 1, x_1 = 0, x_2 = 0 => -1 <= y_12 <= 0, otherwise y_12 = 0
        @variables(cflg, begin
            expres_install[ef_id_ in graph.edge_ids, i in ab], Bin
        end)

        @constraints(cflg, begin
            [ef_id in El], yei[ef_id, :c] <= yei[ef_id, :a]
            [ef_id in El], yei[ef_id, :c] <= yei[ef_id, :b]
            [ef_id in El], yei[ef_id, :c] >= yei[ef_id, :a] + yei[ef_id, :b] - 1
            [ef_id in El], yei[ef_id, :d] <= 1 - yei[ef_id, :c]
            [ef_id in El], q[ef_id, :b] == q[ef_id, :a] + 2 * dlt * (floor(graph.edges[ef_id].length / (2 * dlt) ) - yei[ef_id, :d] ) + qq[ef_id] # long edge coordinate relation
            [ef_id in El], qq[ef_id] <= 2 * dlt * yei[ef_id, :c]
            [ef_id in El], q[ef_id, :a] <= 2 * dlt * (1 - yei[ef_id, :a])  # long edge facility
            [ef_id in El], q[ef_id, :b] >= graph.edges[ef_id].length * yei[ef_id, :b] # long edge coordinate lower bound
            [v_id in Vl], sum( graph.edges[ef_id].etype == :e_long ? (graph.edges[ef_id].nodes[:a] == v_id ? yei[ef_id, :a] : yei[ef_id, :b]) :  0 for ef_id in graph.adjacent_edges[v_id]) <= 1
        end)

        @constraints(cflg, begin
            [ef_id_ in graph.edge_ids, i in ab], expres_install[ef_id_, i] == sum( graph.edges[ef_id].etype == :e_long && ef_id != ef_id_ ? (graph.edges[ef_id].nodes[:a] == graph.edges[ef_id_].nodes[i] ? yei[ef_id, :a] : yei[ef_id, :b]) :  0 for ef_id in graph.adjacent_edges[graph.edges[ef_id_].nodes[i]])
            [ef_id_ in El], q[ef_id_, :a] >= (2 * dlt - problem.c_tol) * expres_install[ef_id_,:a]
            [ef_id_ in El], graph.edges[ef_id_].length - q[ef_id_, :b] >= (2 * dlt - problem.c_tol) * expres_install[ef_id_,:b]
            [ef_id_ in Es], ye[ef_id_] <= 2 - expres_install[ef_id_,:a] - expres_install[ef_id_,:b]
            [ef_id_ in Es], q[ef_id_, :a] >= (graph.edges[ef_id_].length - problem.c_tol) * ( ye[ef_id_] + expres_install[ef_id_,:a] - 1)
            [ef_id_ in Es], graph.edges[ef_id_].length - q[ef_id_, :b] >= (graph.edges[ef_id_].length - problem.c_tol) * ( ye[ef_id_] + expres_install[ef_id_,:b] - 1)
        end)
    end

    # basic constraints
    @constraints(cflg, begin
        # bounds on variables
        [ef_id in El],  ye[ef_id] == 1 # at least one point long edge must be used
        [ef_id in Es], q[ef_id, :a] == q[ef_id, :b]  # same coordinate in short edges
        [ef_id in El], 0 <= q[ef_id, :a]  # leftmost coordinate in long edges
        [ef_id in El], q[ef_id, :a] <= 2 * dlt  # leftmost coordinate in long edges
        [ef_id in El], graph.edges[ef_id].length - 2 * dlt <= q[ef_id, :b] # right coordinate in long edges
        [ef_id in El], q[ef_id, :b] <= graph.edges[ef_id].length # right coordinate in long edges
        # shared constraints on complete covers
        [e_id in Es, ef_id in problem.Ec[e_id]], w[e_id] >= ye[ef_id] # complete cover by edges: open
        [e_id in Es], w[e_id] <=  sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
        [v_id in Vs], x[v_id] >= 1 - sum( (1 - w[e_id]) for e_id in graph.adjacent_edges[v_id])  # ajdacent non covered 2
        # shared constraints on short edge covering
        [e_id in Es], (graph.edges[e_id].length *  (1+problem.cr_tol) + problem.c_tol) * (1 - w[e_id]) <= rv[graph.edges[e_id].nodes[:a]] + rv[graph.edges[e_id].nodes[:b]] # jointly complete cover condition
        # shared constraints on long edge modelling
        [ef_id in El], q[ef_id, :a] <= rv[graph.edges[ef_id].nodes[:a]] + dlt # long edge left cover
        [ef_id in El], graph.edges[ef_id].length - q[ef_id, :b] <= rv[graph.edges[ef_id].nodes[:b]] + dlt # long edge right cover
        # disjunctive activation constraints
        [v_id in graph.node_ids], x[v_id] + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] <= ye[efi[1]] # edge activated constraint
    end)

    if is_indicator
        @constraints(cflg, begin
        [v_id in graph.node_ids], x[v_id] => { rv[v_id] <= 0 } # if x[v]=1 then rv=0
        [v_id in graph.node_ids, efi in EIp[v_id]], ze[v_id, efi] => { rv[v_id] <=  problem.dlte[(v_id, efi[1], efi[2])]  -
        ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1], :a], graph.edges[efi[1]].length - q[efi[1], :b])) }
    end)
    else
        if is_lifted
            @variables(cflg, begin
                0 <= qvei[v_id in graph.node_ids, efi in EIp[v_id]] <= graph.edges[efi[1]].length # disjunctive edge coordinate variable on nodes
            end)
            videfia = [(v_id, efi) for v_id in graph.node_ids for efi in EIp[v_id] if efi[2] == :a]
            videfib = [(v_id, efi) for v_id in graph.node_ids for efi in EIp[v_id] if efi[2] == :b]
            @constraints(cflg, begin
                [v_id in graph.node_ids], rv[v_id]  <=  sum( ( problem.dlte[(v_id, efi[1], efi[2])] - problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] - ifelse( efi[2] == :a , 0, (graph.edges[efi[1]].length) ) ) * ze[v_id, efi] + ifelse(efi[2] == :a , -qvei[v_id, efi], qvei[v_id, efi] ) for efi in EIp[v_id]) # disjunctive residual cover aggreagation constraint
                [(v_id, efi) in videfia], qLBs[efi[1], efi[2]] * ze[v_id, efi] <= qvei[v_id, efi]  # disjunctive upper bound for q on edge
                [(v_id, efi) in videfib], qvei[v_id, efi] <=  qUBs[efi[1], efi[2]]  * ze[v_id, efi] # disjunctive upper bound for q on edge
                [(v_id, efi) in videfib], qLBs[efi[1], efi[2]]  * (1 - ze[v_id, efi]) <= q[efi[1], efi[2]] - qvei[v_id, efi]  # disjunctive upper bound for q on edge
                [(v_id, efi) in videfia], q[efi[1], efi[2]] - qvei[v_id, efi] <=  qUBs[efi[1], efi[2]]  * (1 - ze[v_id, efi]) # disjunctive upper bound for q on edge
                #[v_id in graph.node_ids], rv[v_id] <=  problem.Uv[v_id] * (1 - x[v_id]) # big M on x
                #[v_id in graph.node_ids, efi in EIp[v_id]], rv[v_id] <= problem.Me[(v_id, efi[1], efi[2])]  * (1 - ze[v_id, efi]) + problem.dlte[(v_id, efi[1], efi[2])]  -
                #( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1],:a], (graph.edges[efi[1]].length - q[efi[1],:b]) ))  # big M on edges
            end)
        else
            # bigM constraints
            @constraints(cflg, begin
                [v_id in graph.node_ids], rv[v_id] <=  problem.Uv[v_id] * (1 - x[v_id]) # big M on x
                [v_id in graph.node_ids, efi in EIp[v_id]], rv[v_id] <= problem.Me[(v_id, efi[1], efi[2])]  * (1 - ze[v_id, efi]) + problem.dlte[(v_id, efi[1], efi[2])]  -
                ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1],:a], (graph.edges[efi[1]].length - q[efi[1],:b]) ))  # big M on edges
            end)
        end

        if is_benders
            # Benders decomposition
            if is_lifted
                master_variables = [ye..., yei..., x..., w..., ze..., q..., qq..., rv...]
                sub_variables = [qvei...]
            else
                master_variables = [ye..., yei..., x..., w..., ze..., q..., qq...]
                sub_variables = [rv...]
            end
            setModelAnottion(cflg, master_variables, sub_variables)
            #@constraints(cflg, begin
            #[v_id in graph.node_ids], rv[v_id] <=  problem.Uv[v_id]  * (1 - x[v_id]) # big M on x
            #[v_id in graph.node_ids, efi in EIp[v_id]], rv[v_id] <= (
            # problem.Me[(v_id, efi[1], efi[2])] ) * (1 - ze[v_id, efi]) + ( problem.dlte[(v_id, efi[1], efi[2])] ) -
            # ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1]], (graph.edges[efi[1]].length - q[efi[1]]) ))  # big M on edges
            #end)
        end

    end

    sepatime = 0
    if is_cuts || is_morecuts
        alladded = false
        function user_cut_callback(cb_data, cb_where::Cint)
            if alladded || cb_where != GRB_CB_MIPSOL || cb_where != GRB_CB_MIPNODE
                return
            end
            resultP = Ref{Cdouble}()
            GRBcbget(cb_data, cb_where, GRB_CB_MIPNODE_NODCNT, resultP)
            if resultP[] >= 1.0 && ! is_morecuts
                GRBterminate(backend(model))
                return
            end
            CPUtic()
            Gurobi.load_callback_variable_primal(cb_data, cb_where)
            ze_val = callback_value.(Ref(cb_data), ze)
            q_val = callback_value.(Ref(cb_data), q)
            rv_val = callback_value.(Ref(cb_data), rv)
            for v_id in graph.node_ids
                qvei_expr = Dict()
                qvei_val = Dict()
                for efi in EIp[v_id]
                    if efi[2] == :a
                        # need lower bound
                        lb1 = qLBs[efi[1], efi[2]]  * ze_val[v_id, efi]
                        lb2 = q_val[efi[1], efi[2]] - qUBs[efi[1], efi[2]]  * (1 - ze_val[v_id, efi])
                        if lb1 >= lb2
                            qvei_val[efi] = lb1
                            qvei_expr[efi] = qLBs[efi[1], efi[2]]  * ze[v_id, efi]
                        else
                            qvei_val[efi] = lb2
                            qvei_expr[efi] = q[efi[1], efi[2]] - qUBs[efi[1], efi[2]]  * (1 - ze[v_id, efi])
                        end
                    else
                        # need upper bound
                        ub1 = qUBs[efi[1], efi[2]]  * ze_val[v_id, efi]
                        ub2 = q_val[efi[1], efi[2]] - qLBs[efi[1], efi[2]]  * (1 - ze_val[v_id, efi])
                        if ub1 <= ub2
                            qvei_val[efi] = ub1
                            qvei_expr[efi] = qUBs[efi[1], efi[2]]  * ze[v_id, efi]
                        else
                            qvei_val[efi] = ub2
                            qvei_expr[efi] = q[efi[1], efi[2]] - qLBs[efi[1], efi[2]]  * (1 - ze[v_id, efi])
                        end
                    end
                end
                slack = sum( ( problem.dlte[(v_id, efi[1], efi[2])] - problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] - ifelse( efi[2] == :a , 0, (graph.edges[efi[1]].length) ) ) * ze_val[v_id, efi] + ifelse(efi[2] == :a , -qvei_val[efi], qvei_val[efi] ) for efi in EIp[v_id])
                shouldadd = rv_val[v_id] > slack + 1e-8
                if shouldadd
                    cut = @build_constraint(rv[v_id]  <= 1e-7 + sum( ( problem.dlte[(v_id, efi[1], efi[2])] - problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] - ifelse( efi[2] == :a , 0, (graph.edges[efi[1]].length) ) ) * ze[v_id, efi] + ifelse(efi[2] == :a , -qvei_expr[efi], qvei_expr[efi] ) for efi in EIp[v_id]) )
                    MOI.submit(cflg, MOI.UserCut(cb_data), cut)
                end
            end
            sepatime += CPUtoc()
        end
        set_optimizer_attribute(cflg, "PreCrush", 1)
        MOI.set(cflg, Gurobi.CallbackFunction(), user_cut_callback)
    end

    # objective
    if longedge1
        @objective(cflg, Min,  sum(ye[ef_id] for ef_id in graph.edge_ids) + sum(yei[ef_id, :a] + yei[ef_id, :b] for ef_id in El))
    else
        @objective(cflg, Min,  sum(ye[ef_id] for ef_id in graph.edge_ids) + sum( floor(graph.edges[ef_id].length / (2 * dlt) ) + yei[ef_id, :c] - yei[ef_id, :d] for ef_id in El))
    end

    println("\n model loaded\n")
    optimize!(cflg)
    stat = Stat();
    print("\n sepatime", sepatime, "\n")
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    stat.sepatime = sepatime
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg)
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    if has_values(cflg)
        stat.sol_val = objective_value(cflg)
        for ef_id in graph.edge_ids
            print((ef_id,value(ye[ef_id])))
        end
        print(stat.sol_val)
        for ef_id in El
            print((floor(graph.edges[ef_id].length / (2 * dlt) ), value(yei[ef_id,:a]), " ",value(yei[ef_id,:b]), value(yei[ef_id,:b]), value(q[ef_id, :a]), value(q[ef_id, :b]), value(qq[ef_id]), value(yei[ef_id, :c]) )  )
        end
    end
    stat.formulation = formulation

    stat.instance = problem.instance

    return stat, sol
end

# formulation == EF
function solveEF!(problem::Problem, formulation::FormulationSet, cflg)
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
        cpx = backend(cflg)
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end

    stat.formulation = formulation
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


# edge or edge-vertex model big-M formulations
function solveFPVs!(problem::Problem, formulation::FormulationSet, cflg, Pi, is_benders = false, is_bounding = false)
    # get graph
    graph = problem.prob_graph

    # get delimitation
    if mask(formulation, MSK_ISP0)
        EIp = problem.EIp
        Vp = problem.Vp
    elseif formulation == EVF
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
    if !mask(formulation, MSK_ISP1)
        is_trivial_bound = true
        bigM = problem.dlt * (problem.cr_tol + 1) + problem.c_tol
    end

    usev = !mask(formulation, MSK_ISE)

    # basic variables
    @variables(cflg, begin
        ye[ef_id in graph.edge_ids], Bin # edge facility
        x[v_id in graph.node_ids], Bin #  node residual indictor cover
        w[e_id in graph.edge_ids], Bin # complete cover indicator variable
        0 <= q[ef_id in graph.edge_ids] <= graph.edges[ef_id].length # edge coordinate variable
        0 <= rv[v_id in graph.node_ids] <= (is_trivial_bound ? bigM : problem.Uv[v_id]) # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # big-M modelling variable on edges
    end)

    if is_benders
        master_variables = [ye..., x..., w..., ze..., q...]
        sub_variables = rv
        setModelAnottion(cflg, master_variables, sub_variables)
    end
    # vertex facility related variables
    if usev
        print("vertex var\n")
        @variables(cflg, begin
            yv[vf_id in graph.node_ids], Bin # node facility
            zv[v_id in graph.node_ids, vf_id in Vp[v_id]], Bin # big-M modelling variable on nodes
        end)
    end

    if is_bounding
        @constraints(cflg, begin
            [e_id in graph.edge_ids], ye[e_id] == 0
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
    if !mask(formulation, MSK_ISE) && mask(formulation, MSK_ISV)
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
    elseif mask(formulation, MSK_ISE) && mask(formulation, MSK_ISV)
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
        cpx = backend(cflg)
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    stat.formulation = formulation
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


# long edge vertex model
function solveLEVF!(problem::Problem, formulation::FormulationSet, cflg)

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

    usev = !mask(formulation, MSK_ISE)

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
    @variables(cflg, begin
        yv[vf_id in graph.node_ids], Bin # node facility
        zv[v_id in graph.node_ids, vf_id in Vp[v_id]], Bin # big-M modelling variable on nodes
    end)


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


    @constraints(cflg, begin
        [e_id in Enormal, vf_id in problem.Vc[e_id]], w[e_id] >= yv[vf_id] # complete cover by nodes: open
        [e_id in Enormal], w[e_id] <= sum(yv[vf_id] for vf_id in problem.Vc[e_id]) + sum(ye[ef_id] for ef_id in problem.Ec[e_id]) # complete cover: close
        [ef_id in Enormal, node in [:a,:b]], yv[graph.edges[ef_id].nodes[node]] + ye[ef_id] <= 1 # facility is efither in interior or at end nodes
        [v_id in graph.node_ids], x[v_id] +  sum(zv[v_id, vf_id] for vf_id in Vp[v_id])  + sum(ze[v_id, efi] for efi in EIp[v_id]) == 1 # big-M SOS-1 constraint
        [v_id in graph.node_ids, vf_id in Vp[v_id]], zv[v_id, vf_id] <= yv[vf_id] # node activated constraint
        [v_id in graph.node_ids, vf_id in Vp[v_id]], rv[v_id] <=  problem.Mv[(v_id, vf_id)]  * (1 - zv[v_id, vf_id]) + problem.dltv[(v_id, vf_id)] - problem.d[lor(v_id,vf_id)] # big M on nodes
        [ef_id in Elong], yv[edges[ef_id].nodes[:a]] == 0  #valid inequalities fixing
    end)


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
    @objective(cflg, Min,  sum(yv[vf_id] for vf_id in graph.node_ids) + sum((edges[ef_id].etype == :e_long ? fac_num[ef_id] - u[ef_id] : ye[ef_id]) for ef_id in graph.edge_ids))

    println("\n model loaded\n")
    optimize!(cflg)

    stat = Stat();
    sol = Vector{Tuple{Symbol,Int, Float64}}();
    stat.termintaion_status = termination_status(cflg)
    stat.bound = objective_bound(cflg)
    stat.gap = relative_gap(cflg)
    stat.time = solve_time(cflg)
    if solver_name(cflg) == "CPLEX"
        cpx = backend(cflg)
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end

    stat.formulation = formulation
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
        for v_id in graph.node_ids
            if value(yv[v_id]) >= 0.5
                push!(sol, (:v, v_id,  0))
            end
        end
    end
    #println(graph.edges)
    return stat, sol
end


# edge model disjunctive formulation
function solveEFPD!(problem::Problem, formulation::FormulationSet, cflg, is_benders = false)
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
        0 <= barqvei[v_id in graph.node_ids, efi in EIp[v_id]] <= graph.edges[efi[1]].length  # disjunctive edge coordinate variable on nodes
        0 <= qvei[v_id in graph.node_ids, efi in EIp[v_id]] <= graph.edges[efi[1]].length # disjunctive edge coordinate variable on nodes
        0 <= rv[v_id in graph.node_ids] <= problem.Uv[v_id] # residual cover variable
        ze[v_id in graph.node_ids, efi in EIp[v_id]], Bin # disjunctive indicator variable on edges
    end)

    if is_benders
        master_variables = [ye..., x..., w..., ze..., q...]
        sub_variables = [rv..., barqvei..., qvei...]
        setModelAnottion(cflg, master_variables, sub_variables)
        #@constraints(cflg, begin
        #[v_id in graph.node_ids], rv[v_id] <=  problem.Uv[v_id]  * (1 - x[v_id]) # big M on x
        #[v_id in graph.node_ids, efi in EIp[v_id]], rv[v_id] <= (
        # problem.Me[(v_id, efi[1], efi[2])] ) * (1 - ze[v_id, efi]) + ( problem.dlte[(v_id, efi[1], efi[2])] ) -
        # ( problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] +  ifelse(efi[2] == :a , q[efi[1]], (graph.edges[efi[1]].length - q[efi[1]]) ))  # big M on edges
        #end)
    end

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
        [v_id in graph.node_ids], rv[v_id]  ==  sum( ( problem.dlte[(v_id, efi[1], efi[2])] - problem.d[lor(v_id, graph.edges[efi[1]].nodes[efi[2]])] - ifelse( efi[2] == :a , 0, (graph.edges[efi[1]].length) ) ) * ze[v_id, efi] + ifelse(efi[2] == :a , -qvei[v_id, efi], qvei[v_id, efi] ) for efi in EIp[v_id]) # disjunctive residual cover aggreagation constraint
        [v_id in graph.node_ids, efi in EIp[v_id]], q[efi[1]] ==   barqvei[v_id, efi] + qvei[v_id, efi]  # disjunctive residual q aggregation constraint
        [v_id in graph.node_ids, efi in EIp[v_id]], qvei[v_id, efi] <=  graph.edges[efi[1]].length * (1+problem.cr_tol)  * ze[v_id, efi] # disjunctive upper bound for q on edge
        [v_id in graph.node_ids, efi in EIp[v_id]], barqvei[v_id, efi] <=  graph.edges[efi[1]].length * (1+problem.cr_tol)  * x[v_id] # disjunctive upper bound for q on edge
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
        cpx = backend(cflg)
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    stat.formulation = formulation

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


# Edge model indicator formulation indicator
function solveEFPI!(problem::Problem, formulation::FormulationSet, cflg)
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
        cpx = backend(cflg)
        stat.node = CPLEX.CPXgetnodecnt(cpx.env, cpx.lp)
    else
        stat.node = Int32(node_count(cflg))
    end
    stat.formulation = formulation

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



