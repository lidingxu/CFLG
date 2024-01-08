#=========================================================
 utility 
=========================================================#

@inline lor(a::Int, b::Int)= min(a,b), max(a,b)
@inline non(MSK::UInt16)= ( MSK & UInt16(0b0000000000000000) )


MSK_ZERO =  UInt16(0b0000000000000000)
MSK_ISMOD = UInt16(0b1000000000000000)   # is this a model
MSK_ISE =   UInt16(0b0100000000000000)   # edge or edge-vertex model
MSK_ISP0 =  UInt16(0b0010000000000000)   # is processing by delimitation
MSK_ISP1 =  UInt16(0b0001000000000000)   # is processing by bound tightenning
MSK_ISV =   UInt16(0b0000100000000000)   # is using valid inequalities
MSK_ISK2 =  UInt16(0b0000010000000000)   # K = 2 
MSK_ISD =   UInt16(0b0000001000000000)   # is disjunctive programming 
MSK_ISL =   UInt16(0b0000000100000000)   # is a long edge model
MSK_ISLg =  UInt16(0b0000000010000000)   # is Logical 
MSK_ISC  =  UInt16(0b0000000001000000)   # is cover preprocessing

@inline mask(algo, MSK::UInt16)= ( UInt16(algo) & MSK != MSK_ZERO )


@enum AlgorithmSet::UInt16 begin
    EF     =   MSK_ISMOD | MSK_ISE                                           # edge formulation, from "Covering edges in networks", Fröhlich et al.
    EFP0   =   MSK_ISMOD | MSK_ISE | MSK_ISP0                                # edge formulation with simple processing (delimited cover) 
    EFP    =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1                     # edge formulation with processing (bound tightenning and delimited cover)
    EFPC   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISC           # edge formulation with processing (bound tightenning and delimited cover) and cover preprocessing
    EFPV   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISV           # edge formulation with processing (bound tightenning and delimited cover) and valid inequalities
    EFPV2  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISV | MSK_ISK2# edge-vertex formulation with processing (bound tightenning and delimited cover) and valid inequalities (K = 2)    
    EFPD   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISD 	     # edge disjunctive programming formulation with processing (delimited cover)
    EFPDC  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISD | MSK_ISC # edge disjunctive programming formulation with processing (delimited cover) and cover preprocessing
    EFPLg  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISLg          # edge logical formulation with processing (delimited cover)
    EFPL   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISL           # long edge formulation with processing (bound tightenning and delimited cover)
    EVF    =   MSK_ISMOD                                                     # edge-vertex formulation 
    EVFP0  =   MSK_ISMOD | MSK_ISP0                                          # edge-vertex formulation with simple processing (delimited cover)
    EVFP   =   MSK_ISMOD | MSK_ISP0 | MSK_ISP1                               # edge-vertex formulation with processing (bound tightenning and delimited cover)
    EVFPV  =   MSK_ISMOD | MSK_ISP0 | MSK_ISP1 | MSK_ISV                     # edge-vertex formulation with processing (bound tightenning and delimited cover) and valid inequalities
    EVFPL  =   MSK_ISMOD | MSK_ISP0 | MSK_ISP1 | MSK_ISL                     # long edge-vertex formulation with processing (bound tightenning and delimited cover)
    None   =   MSK_ZERO                                                      # not a model, just record statistics of original graph, degree-2-free graph, subdivided graph
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

    function Option(time_limit::Float64=360.0,  rel_gap::Float64=1e-4, log_level::Int=1, thread::Int=1, silent::Bool=true)
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



function initModel(solver_name::String,  option::Option, time_limit_sec)
    # set solver
    if solver_name == "Gurobi"  
        #@eval import Gurobi
        #model = optimizer_with_attributes(Gurobi.Optimizer, "Threads" => option.thread, "MIPGap" => option.rel_gap, "MIPGAPABS" => 1,  "TimeLimit" => option.time_limit)
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "Threads", option.thread)
        set_optimizer_attribute(model, "OutputFlag", option.log_level)
        set_optimizer_attribute(model, "MIPGap", option.rel_gap)
        set_optimizer_attribute(model, "MIPGapAbs", 1)
        set_optimizer_attribute(model, "TimeLimit", time_limit_sec) # note Gurobi only supports wall-clock time
    elseif solver_name == "CPLEX"
        #model = optimizer_with_attributes(CPLEX.Optimizer, "CPXPARAM_Threads" => option.thread, "CPXPARAM_MIP_Tolerances_MIPGap" => option.rel_gap, "CPXPARAM_MIP_Tolerances_AbsMIPGap" => 1, "CPXPARAM_ClockType" => 1, "CPXPARAM_TimeLimit" => option.time_limit)
        model =  JuMP.Model(() -> CPLEX.Optimizer())
        set_optimizer_attribute(model, "CPXPARAM_Threads", option.thread)
        set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_MIPGap", option.rel_gap)
        set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_AbsMIPGap", 1)
        set_optimizer_attribute(model, "CPXPARAM_ClockType", 1) # CPU clock time
        set_optimizer_attribute(model, "CPXPARAM_TimeLimit", time_limit_sec) # n
    elseif solver_name == "GLPK"
        model = Model(GLPK.Optimizer)
        set_optimizer_attribute(model, "mip_gap", option.rel_gap)
        set_optimizer_attribute(model, "tm_lim", time_limit_sec) # n
        #to do
    elseif solver_name == "SCIP"
        model = Model(SCIP.Optimizer)
        set_optimizer_attribute(model, "limits/gap", option.rel_gap)
        set_optimizer_attribute(model, "limits/time", time_limit_sec) # n
        set_optimizer_attribute(model, "limits/absgap", 1)
        set_optimizer_attribute(model, "timing/clocktype", 1)
    else
        println("unkown solver name\n")
        @assert(false)
    end
    return model
end
