#=========================================================
 benchmark test
 the main function receive the arguments: #instance_directory #solver_name #time_limit #output_directory #instance_name #algorithm #cover
=========================================================#

include("../src/CFLG.jl")
using .CFLG


function main(args)
    @assert(length(args) == 7)
    instance_dir = args[1]
    solver_name = args[2]
    time_limit =  parse(Float64, args[3])
    output_dir = args[4]
    instance_name = args[5]
    algorithm = args[6]
    cover= args[7]

    abs_instance =  string(instance_dir ,"/", instance_name)
    solver_names = ["Gurobi", "CPLEX", "GLPK", "SCIP"]
    covers = ["Small", "Large"]
    println(instance_dir, " ", solver_name, " ", time_limit, " ", output_dir, " ", instance_name, " ", algorithm, " ", cover)
    @assert(solver_name in solver_names)
    @assert(cover in covers)
    graph = readGraph(abs_instance)
    option = Option(time_limit)
    if cover == "Small"
        dlt = float(graph.avg_len)
    else    cover == "Large"
        dlt = float( 2 * graph.avg_len)
    end
    println("data loaded\n")
    problem = Problem(graph, dlt)
    stat = solve!(problem, solver_name, option, algorithm)
    abs_output =  string(output_dir , "/" , instance_name , "." , algorithm , "." , cover) 
    if algorithm == "None"
        stat_info = string("instance: ", instance_name, "\n", 
        "org_node: ", string(stat[1].num_node), "\n", "org_edge: ", string(stat[1].num_edge),  "\n", "org_min_len: ", string(stat[1].min_len),  "\n", "org_avg_len: ", string(stat[1].avg_len),  "\n", "org_max_len: ", string(stat[1].max_len),  "\n",
        "dtf_node: ", string(stat[2].num_node), "\n", "dtf_edge: ", string(stat[2].num_edge),  "\n", "dtf_min_len: ", string(stat[2].min_len),  "\n", "dtf_avg_len: ", string(stat[2].avg_len),  "\n", "dtf_max_len: ", string(stat[2].max_len),  "\n",
        "sdb_node: ", string(stat[3].num_node), "\n", "sdb_edge: ", string(stat[3].num_edge),  "\n", "sdb_min_len: ", string(stat[3].min_len),  "\n", "sdb_avg_len: ", string(stat[3].avg_len),  "\n", "sdb_max_len: ", string(stat[3].max_len),  "\n", 
        "dlt: ", string(dlt), "\n", "algo: ", algorithm)        
    else
        stat_info = string("instance: ", instance_name, "\n", "obj: ", string(stat.sol_val), "\n", "bound: ", string(stat.bound), "\n", "gap: ", string(stat.gap), "\n", "time: ", string(stat.time), "\n", 
        "preprocess_time: ", string(stat.preprocess_time), "\n",  "node: ", string(stat.node), "\n", "algo: ", algorithm)
    end
    io = open(abs_output, "w")
    print(io, stat_info)
    close(io)
end


main(ARGS)
