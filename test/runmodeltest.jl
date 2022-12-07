
#=========================================================
 model test
=========================================================#


using CFLG

@testset "test models" begin    
    solver_names = ["CPLEX"]
    for solver_name in solver_names
        if solver_name == "Gurobi"
	    try 
		import Gurobi
	    catch e
	    end
        elseif solver_name == "CPLEX"
	    try 
		import CPLEX
	    catch e
	    end
        elseif solver_name == "GLPK"
	    try 
		import GLPK
	    catch e
	    end
        else
            @test false
            println("unkown solver name\n")
        end
        graph = readGraph("../benchmarks/test/r_40_0.5_388.txt")
        #@test graph != nothing
        default_option = Option()
        for algo in ["DF"]          
            problem = Problem(graph, Float64(graph.avg_len))
            solve!(problem, solver_name, default_option, algo)
        end
    end
end 

