#=========================================================
 the main function receive the arguments: #solver_name #time_limit #instance_name #algorithm #raidus
=========================================================#

include("src/CFLG.jl")
using .CFLG


function main(args)
    @assert(length(args) == 5)
    solver_name = args[1]
    time_limit =  parse(Float64, args[2])
    instance = args[3]
    algorithm = args[4]
    dlt= parse(Float64, args[5])

    algorithms = ["F0", "F", "SF", "SFD", "RF", "EF", "None"]
    solver_names = ["Gurobi", "CPLEX", "GLPK"]
    println( solver_name, " ", time_limit, " ", instance, " ", algorithm)
    @assert(algorithm in algorithms)
    @assert(solver_name in solver_names)
    graph = readGraph(instance)
    option = Option(time_limit)
    println("data loaded\n")
    problem = Problem(graph, dlt)
    if algorithm == "F0"
        algo = F0
    elseif algorithm == "F"
        algo = F
    elseif algorithm == "SF"
        algo = SF
    elseif algorithm == "SFD"
        algo = SFD
    elseif algorithm == "RF"
        algo = RF
    elseif algorithm == "EF"
        algo = EF
    elseif algorithm == "None"
        algo = None
    end
    stat = solve!(problem, solver_name, option, algo)
    if algo == None
        stat_info = string("instance: ", instance, "\n", 
        "org_node: ", string(stat[1].num_node), "\n", "org_edge: ", string(stat[1].num_edge),  "\n", "org_min_len: ", string(stat[1].min_len),  "\n", "org_avg_len: ", string(stat[1].avg_len),  "\n", "org_max_len: ", string(stat[1].max_len),  "\n",
        "dtf_node: ", string(stat[2].num_node), "\n", "dtf_edge: ", string(stat[2].num_edge),  "\n", "dtf_min_len: ", string(stat[2].min_len),  "\n", "dtf_avg_len: ", string(stat[2].avg_len),  "\n", "dtf_max_len: ", string(stat[2].max_len),  "\n",
        "sdb_node: ", string(stat[3].num_node), "\n", "sdb_edge: ", string(stat[3].num_edge),  "\n", "sdb_min_len: ", string(stat[3].min_len),  "\n", "sdb_avg_len: ", string(stat[3].avg_len),  "\n", "sdb_max_len: ", string(stat[3].max_len),  "\n", 
        "dlt: ", string(dlt), "\n", "algo: ", string(algo))        
    else
        stat_info = string("instance: ", instance, "\n", "obj: ", string(stat.sol_val), "\n", "bound: ", string(stat.bound), "\n", "gap: ", string(stat.gap), "\n", "time: ", string(stat.time), "\n", 
        "preprocess_time: ", string(stat.preprocess_time), "\n",  "node: ", string(stat.node), "\n", "algo: ", string(stat.algo))
    end
    print(stat_info)
end


main(ARGS)