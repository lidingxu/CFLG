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
MSK_ISK2 =  UInt16(0b0000010000000000)   # if using valid inequalities, whether the rank is 2, otherwise the rank is 1 by default
MSK_ISD =   UInt16(0b0000001000000000)   # is disjunctive programming formulation
MSK_ISL =   UInt16(0b0000000100000000)   # is using long edge modeling
MSK_ISI =   UInt16(0b0000000010000000)   # is using indicator constraint modelling
MSK_ISC  =  UInt16(0b0000000001000000)   # is cover preprocessing
MSK_ISB  =  UInt16(0b0000000000100000)   # is Benders decomposition

@inline mask(formulation, MSK::UInt16)= ( UInt16(formulation) & MSK != MSK_ZERO )

@enum FormulationSet::UInt16 begin
    EF     =   MSK_ISMOD | MSK_ISE                                           # edge model formulation, from "Covering edges in networks", Fröhlich et al.
    EFP0   =   MSK_ISMOD | MSK_ISE | MSK_ISP0                                # edge model formulation with simple processing (delimited cover)
    EFP    =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1                     # edge model big-M formulation with processing (bound tightenning and delimited cover)
    EFPB   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISB           # edge model big-M formulation with processing (bound tightenning and delimited cover) and Benders decomposition
    EFPC   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISC           # edge model big-M formulation with processing (bound tightenning and delimited cover) and cover preprocessing
    EFPV   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISV           # edge model big-M formulation with processing (bound tightenning and delimited cover) and rank-1 valid inequalities
    EFPV2  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISV | MSK_ISK2# edge model big-M formulation with processing (bound tightenning and delimited cover) and rank-2 valid inequalities
    EFPD   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISD 	         # edge model disjunctive programming formulation with processing (delimited cover)
    EFPDB  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISD | MSK_ISB # edge model disjunctive programming formulation with processing (delimited cover) and Benders decomposition
    EFPDC  =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISD | MSK_ISC # edge model disjunctive programming formulation with processing (delimited cover) and cover preprocessing
    EFPI   =   MSK_ISMOD | MSK_ISE | MSK_ISP0 | MSK_ISP1 | MSK_ISI           # edge model indicator constraint formulation with processing (delimited cover)
    EVF    =   MSK_ISMOD                                                     # edge-vertex model big-M formulation
    EVFP0  =   MSK_ISMOD | MSK_ISP0                                          # edge-vertex model big-M formulation  with simple processing (delimited cover)
    EVFP   =   MSK_ISMOD | MSK_ISP0 | MSK_ISP1                               # edge-vertex model big-M formulation  with processing (bound tightenning and delimited cover)
    EVFPV  =   MSK_ISMOD | MSK_ISP0 | MSK_ISP1 | MSK_ISV                     # edge-vertex model big-M formulation  with processing (bound tightenning and delimited cover) and simple valid inequalities
    LEVFP  =   MSK_ISMOD | MSK_ISP0 | MSK_ISP1 | MSK_ISL                     # long edge-vertex model big-M with processing (bound tightenning and delimited cover)
    LEFPI  =   MSK_ISMOD | MSK_ISE  | MSK_ISP0 | MSK_ISP1 | MSK_ISL | MSK_ISI                      # edge model big-M formulation with processing (bound tightenning and delimited cover) and Benders decomposition
    LEFP   =   MSK_ISMOD | MSK_ISE  | MSK_ISP0 | MSK_ISP1 | MSK_ISL                       # edge model big-M formulation with processing (bound tightenning and delimited cover) and Benders decomposition
    LEFPV  =   MSK_ISMOD | MSK_ISE  | MSK_ISP0 | MSK_ISP1 | MSK_ISL | MSK_ISV
    LEFPB  =   MSK_ISMOD | MSK_ISE  | MSK_ISP0 | MSK_ISP1 | MSK_ISB | MSK_ISL             # edge model big-M formulation with processing (bound tightenning and delimited cover) and Benders decomposition
    LEFPD  =   MSK_ISMOD | MSK_ISE  | MSK_ISP0 | MSK_ISP1 | MSK_ISD | MSK_ISL             # edge model disjunctive programming formulation with processing (delimited cover)
    LEFPDB =   MSK_ISMOD | MSK_ISE  | MSK_ISP0 | MSK_ISP1 | MSK_ISD | MSK_ISB | MSK_ISL   # edge model disjunctive programming formulation with processing (delimited cover) and Benders decomposition
    None   =   MSK_ZERO                                                      # not a model, just record statistics of original graph, degree-2-free graph, subdivided graph
end

# statistics of a graph
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


# solver option
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

# result statistics
mutable struct Stat
    termintaion_status # termintaion status
    sol_val::Float64 # best solution found
    bound::Float64 # best bound found
    gap::Float64 # gap obtained
    preprocess_time::Float64 # preprocess time
    time::Float64 # total time
    node::Int32 # nodes of search tree
    formulation::FormulationSet
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
        model = Model(Gurobi.Optimizer)
        set_optimizer_attribute(model, "Threads", option.thread)
        set_optimizer_attribute(model, "OutputFlag", option.log_level)
        set_optimizer_attribute(model, "MIPGap", option.rel_gap)
        set_optimizer_attribute(model, "MIPGapAbs", 1)
        set_optimizer_attribute(model, "TimeLimit", time_limit_sec) # note Gurobi only supports wall-clock time
    elseif solver_name == "CPLEX"
        model = direct_model(CPLEX.Optimizer())
        set_optimizer_attribute(model, "CPXPARAM_Threads", option.thread)
        set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_MIPGap", option.rel_gap)
        set_optimizer_attribute(model, "CPXPARAM_MIP_Tolerances_AbsMIPGap", 1)
        set_optimizer_attribute(model, "CPXPARAM_ClockType", 1) # CPU clock time
        set_optimizer_attribute(model, "CPXPARAM_TimeLimit", time_limit_sec)
    elseif solver_name == "GLPK"
        model = Model(GLPK.Optimizer)
        set_optimizer_attribute(model, "mip_gap", option.rel_gap)
        set_optimizer_attribute(model, "tm_lim", time_limit_sec) #
        #to do
    elseif solver_name == "SCIP"
        model = Model(SCIP.Optimizer)
        set_optimizer_attribute(model, "limits/gap", option.rel_gap)
        set_optimizer_attribute(model, "limits/time", time_limit_sec)
        set_optimizer_attribute(model, "limits/absgap", 1)
        set_optimizer_attribute(model, "timing/clocktype", 1) # CPU clock time
    else
        println("unkown solver name\n")
        @assert(false)
    end
    return model
end



function add_annotation(
    model::JuMP.Model,
    variable_classification::Dict;
    all_variables::Bool = true,
)
    num_variables = sum(length(it) for it in values(variable_classification))
    if all_variables
        @assert num_variables == JuMP.num_variables(model)
    end
    indices, annotations = CPXINT[], CPXLONG[]
    for (key, value) in variable_classification
        for variable_ref in value
            push!(indices, variable_ref.index.value - 1)
            push!(annotations, CPX_BENDERS_MASTERVALUE + key)
        end
    end
    cplex = backend(model)
    index_p = Ref{CPXINT}()
    CPXnewlongannotation(
        cplex.env,
        cplex.lp,
        CPX_BENDERS_ANNOTATION,
        CPX_BENDERS_MASTERVALUE,
    )
    CPXgetlongannotationindex(
        cplex.env,
        cplex.lp,
        CPX_BENDERS_ANNOTATION,
        index_p,
    )
    CPXsetlongannotations(
        cplex.env,
        cplex.lp,
        index_p[],
        CPX_ANNOTATIONOBJ_COL,
        length(indices),
        indices,
        annotations,
    )
    return
end

function setModelAnottion(model::JuMP.Model, master_variables, sub_variables)
    set_optimizer_attribute(model, "CPXPARAM_Benders_Strategy", 1)
    variable_classification = Dict(0 => master_variables, 1 => sub_variables)
    add_annotation(model, variable_classification)
end
