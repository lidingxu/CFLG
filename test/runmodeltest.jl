
#=========================================================
 model test
=========================================================#


@testset "test models" begin    
    solver_names = ["CPLEX"]
    for solver_name in solver_names
        if solver_name == "Gurobi"
            grb = try_import(:Gurobi)
            @test grb
        elseif solver_name == "CPLEX"
            cpx = try_import(:CPLEX)
            @test cpx
        elseif solver_name == "GLPK"
            glpk = try_import(:GLPK)
            @test glpk
        else
            @test false
            println("unkown solver name\n")
        end
        graph = readGraph("../benchmarks/test/city_132.txt")
        #@test graph != nothing
        default_option = Option()
        for algo in [F0]          
            problem = Problem(graph, Float64(graph.avg_len))
            solve!(problem, solver_name, default_option, algo)
        end
    end
end 

