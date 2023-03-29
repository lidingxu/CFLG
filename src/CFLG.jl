module CFLG
    using DataStructures
    using JuMP
    using CPUTime
    using Combinatorics
    try 
        import CPLEX
    catch e
    end

    try 
        import Gurobi
    catch e
    end

    try 
        import GLPK
    catch e
    end

    try 
        import SCIP
    catch e
    end

    include("utils.jl")
    include("Graph.jl")
    include("Reader.jl")
    include("Problem.jl")
    include("Algorithm.jl")
    export readGraph, Option, AlgorithmSet, GraphStat, Stat, Problem, solve!
end
