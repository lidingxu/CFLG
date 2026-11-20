#=========================================================
 benchmark test
 the main function receives the arguments:
   #instance_directory #instance_name #output_directory #solver_name #time_limit #formulation #cover [#seed]
 with a seed, the result file is <instance>.<formulation>.<cover>.s<seed>, otherwise <instance>.<formulation>.<cover>
 next to the result file: .start (written before solving, so a crash leaves a .start without a result) and
 .log (solver log of the main model)
=========================================================#

include("../src/CFLG.jl")
using .CFLG
using Sockets: gethostname


const WARMUP_INSTANCE = joinpath(@__DIR__, "..", "benchmarks", "test", "r_20_0.1_23.txt")

radius(graph, cover) = cover == "Small" ? float(0.25 * graph.avg_len) : float(2 * graph.avg_len)

# solve a tiny instance with the same formulation first, so that Julia's compilation is not timed
function warmup(solver_name, formulation, cover)
    graph = readGraph(WARMUP_INSTANCE)
    problem = Problem(graph, radius(graph, cover))
    redirect_stdout(devnull) do
        solve!(problem, solver_name, Option(30.0, 1e-4, 0), formulation)
    end
    return
end

function main(args)
    @assert(length(args) == 7 || length(args) == 8)
    # parse arguments
    instance_dir = args[1]
    instance_name = args[2]
    output_dir = args[3]
    solver_name = args[4]
    time_limit =  parse(Float64, args[5])
    formulation = args[6]
    cover= args[7]
    has_seed = length(args) == 8
    seed = has_seed ? parse(Int, args[8]) : 0

    solver_names = ["Gurobi", "GLPK", "SCIP"]
    covers = ["Small", "Large"]
    println(instance_dir, " ", solver_name, " ", time_limit, " ", output_dir, " ", instance_name, " ", formulation, " ", cover, " seed ", seed)
    @assert(solver_name in solver_names)
    @assert(cover in covers)

    abs_output = string(output_dir, "/", instance_name, ".", formulation, ".", cover, has_seed ? ".s$seed" : "")
    mkpath(output_dir)
    host = gethostname()
    cpu = Sys.cpu_info()[1].model
    open(abs_output * ".start", "w") do io
        print(io, "host: ", host, "\n", "cpu: ", cpu, "\n", "started: ", time(), "\n")
    end

    if formulation != "None"
        warmup(solver_name, formulation, cover)
    end

    # read graph and set option
    abs_instance =  string(instance_dir ,"/", instance_name)
    println(abs_instance)
    graph = readGraph(abs_instance)
    log_file = formulation == "None" ? "" : abs_output * ".log"
    log_file != "" && rm(log_file; force = true)
    option = Option(time_limit, 1e-4, 1, 1, true, seed, log_file)
    dlt = radius(graph, cover)
    println("data loaded\n")
    problem = Problem(graph, dlt)
    stat = solve!(problem, solver_name, option, formulation)
    if formulation == "None"
        stat_info = string("instance: ", instance_name, "\n",
        "org_node: ", string(stat[1].num_node), "\n", "org_edge: ", string(stat[1].num_edge),  "\n", "org_min_len: ", string(stat[1].min_len),  "\n", "org_avg_len: ", string(stat[1].avg_len),  "\n", "org_max_len: ", string(stat[1].max_len),  "\n",
        "dtf_node: ", string(stat[2].num_node), "\n", "dtf_edge: ", string(stat[2].num_edge),  "\n", "dtf_min_len: ", string(stat[2].min_len),  "\n", "dtf_avg_len: ", string(stat[2].avg_len),  "\n", "dtf_max_len: ", string(stat[2].max_len),  "\n",
        "sdb_node: ", string(stat[3].num_node), "\n", "sdb_edge: ", string(stat[3].num_edge),  "\n", "sdb_min_len: ", string(stat[3].min_len),  "\n", "sdb_avg_len: ", string(stat[3].avg_len),  "\n", "sdb_max_len: ", string(stat[3].max_len),  "\n",
        "dlt: ", string(dlt), "\n", "formulation: ", formulation)
    else
        extra = stat.extra
        # time: solver time; preprocess_time: preprocessing before the model is built; build_time: JuMP model
        # construction; total_time: wall clock of the whole solve! call (bounded by the time limit)
        stat_info = string("instance: ", instance_name, "\n", "obj: ", string(stat.sol_val), "\n", "bound: ", string(stat.bound), "\n", "gap: ", string(stat.gap), "\n", "time: ", string(stat.time), "\n",
        "preprocess_time: ", string(stat.preprocess_time), "\n", "sepa_time: ", string(stat.sepatime), "\n",  "node: ", string(stat.node), "\n", "formulation: ", formulation, "\n",
        "status: ", string(stat.termintaion_status), "\n", "cover: ", cover, "\n", "seed: ", string(seed), "\n",
        "build_time: ", string(get(extra, "build_time", 0.0)), "\n", "total_time: ", string(get(extra, "total_time", NaN)), "\n",
        "host: ", host, "\n", "cpu: ", cpu, "\n", "solver_version: ", CFLG.solverVersion(solver_name), "\n", "julia_version: ", string(VERSION))
        for key in sort(collect(keys(extra)))
            key in ("build_time", "total_time") && continue
            stat_info *= string("\n", key, ": ", string(extra[key]))
        end
    end
    io = open(abs_output, "w")
    print(io, stat_info)
    close(io)
end


main(ARGS)
