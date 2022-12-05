
#=========================================================
 graph test
=========================================================#


@testset "test graph functions" begin    
    graph = readGraph("../benchmarks/testbenchmark/testa.graph")
    @test graph != nothing
    printGraph(graph)
    (Ec, Vc, Ep, Vp, EIp, d) = graphProcess(graph, 2.1) 
    println("Ec:",Ec)
    println("Vc:",Vc)
    println("Ep:", Ep)
    println("Vp:", Vp)
    println("EIp:", EIp)
    println("d:", d)
end 





