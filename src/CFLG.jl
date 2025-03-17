module CFLG
    using DataStructures
    using JuMP
    using CPUTime
    using MutableNamedTuples
    try
        using CPLEX
    catch e
    end

    try
        using Gurobi
    catch e
    end

    try
        using GLPK
    catch e
    end

    try
        using SCIP
    catch e
    end

    include("utils.jl")
    include("Graph.jl")
    include("Reader.jl")
    include("Problem.jl")
    include("Infeas.jl")
    include("Algorithm.jl")
    export readGraph, Option, FormulationSet, GraphStat, Stat, Problem, solve!
end
