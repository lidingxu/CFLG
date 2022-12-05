#=========================================================
 utility 
=========================================================#

@inline lor(a::Int, b::Int)= min(a,b), max(a,b)


@enum AlgorithmSet begin
    SFD = 0 # discrete vertex formulation, bound tightenning, complete cover, valid inequalities
    RF = 1 # vertex formulation with long edge, bound tightenning, complete cover
    SF = 2 # vertex formulation, bound tightenning, complete cover, valid inequalities
    VF_BTCC = 3 # vertex formulation, bound tightenning, complete cover
    F = 4 # vertex formulation, trivial bound: delta
    F0 = 5 # vertex formulation, full cover
    EF = 6 # edge formulation, from "Covering edges in networks", Fröhlich et al.
    None= 7 # just record statistics of original graph, degree-2-free graph, subdivided graph
end

struct GraphStat
    num_node::Int # node number
    num_edge::Int # node number
    min_len::Float64 # minimum length of edfes
    max_len::Float64 # maximum edge length
    avg_len::Float64 # average edge length

    function GraphStat(num_node::Int, num_edge::Int, min_len::Float64, max_len::Float64, avg_len::Float64)
        new(num_node, num_edge, min_len, max_len, avg_len)
    end
end


struct Option
    time_limit::Float64 # time limit
    rel_gap::Float64 # relative duality gap
    log_level::Int # log level
    thread::Int # thread number
    silent::Bool # silent model

    function Option(time_limit::Float64=40.0, rel_gap::Float64=1e-4, log_level::Int=1, thread::Int=1, silent::Bool=true)
        new(time_limit, rel_gap, log_level, thread, silent)
    end
end

mutable struct Stat
    termintaion_status
    sol_val::Float64
    bound::Float64
    gap::Float64
    preprocess_time::Float64
    time::Float64
    node::Int32
    algo::AlgorithmSet
    instance::String
    solver_name::String
    

    function Stat()
        new()
    end
end

function try_import(name::Symbol)
    try
        @eval import $name
        return true
    catch e
        return false
    end
end