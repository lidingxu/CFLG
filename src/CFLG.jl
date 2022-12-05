module CFLG
    using DataStructures
    using JuMP
    using ExportAll
    using CPUTime
    include("utils.jl")
    include("Graph.jl")
    include("Reader.jl")
    include("Problem.jl")
    include("Algorithm.jl")

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

    @exportAll()
    #export try_import
    #export Node, Edge, Graph, printGraph, graphProcess
    #export readGraph
    #export Option, Problem
end
