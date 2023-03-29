#=========================================================
 utility 
=========================================================#

@inline lor(a::Int, b::Int)= min(a,b), max(a,b)
@inline non(MSK::UInt8)= ( MSK & UInt8(0b00000000) )


MSK_ZERO =  UInt8(0b00000000)
MSK_ISMOD = UInt8(0b10000000)   # is this a model
MSK_ISE =   UInt8(0b01000000)   # edge or edge-vertex model
MSK_ISP0 =  UInt8(0b00100000)   # is processing by delimitation
MSK_ISP1 =  UInt8(0b00010000)   # is processing by bound tightenning
MSK_ISV =   UInt8(0b00001000)   # is using valid inequalities 
MSK_ISD =   UInt8(0b00000100)   # is disjunctive programming 
MSK_ISL =   UInt8(0b00000010)   # is a long edge model

@inline mask(algo, MSK::UInt8)= ( UInt8(algo) & MSK != MSK_ZERO )


@enum AlgorithmSet::UInt8 begin
    EF =     MSK_ISMOD | MSK_ISE                                  # edge formulation, from "Covering edges in networks", Fröhlich et al.
    EFP0 =   MSK_ISMOD | MSK_ISE | MSK_ISP0                       # edge formulation with simple processing (delimited cover) 
    EFP  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1            # edge formulation with processing (bound tightenning and delimited cover)
    EFPV =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISV  # edge formulation with processing (bound tightenning and delimited cover) and valid inequalities
    EFPD =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISD             # edge disjunctive programming formulation with processing (delimited cover)
    EFPL =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISL  # long edge formulation with processing (bound tightenning and delimited cover)
    EVF =    MSK_ISMOD                                            # edge-vertex formulation 
    EVFP0 =  MSK_ISMOD | MSK_ISP0                                 # edge-vertex formulation with simple processing (delimited cover)
    EVFP  =  MSK_ISMOD | MSK_ISP0 | MSK_ISP1                      # edge-vertex formulation with processing (bound tightenning and delimited cover)
    EVFPV =  MSK_ISMOD | MSK_ISP0 | MSK_ISP1 | MSK_ISV            # edge-vertex formulation with processing (bound tightenning and delimited cover) and valid inequalities
    EVFPL =  MSK_ISMOD | MSK_ISP0 | MSK_ISP1 | MSK_ISL            # long edge-vertex formulation with processing (bound tightenning and delimited cover)
    None =   MSK_ZERO                                             # not a model, just record statistics of original graph, degree-2-free graph, subdivided graph
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

    function Option(time_limit::Float64=360.0, rel_gap::Float64=1e-4, log_level::Int=1, thread::Int=1, silent::Bool=true)
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
