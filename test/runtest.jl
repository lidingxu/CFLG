#=========================================================
 unit tests
=========================================================#

include("../src/CFLG.jl")
using .CFLG
using Test
using Printf


tests = [
    #"rungraphtest.jl",
    "runmodeltest.jl"
    ]

@info("starting all tests")
println()
timings = Dict{String, Float64}()
@testset "all tests" begin
    all_test_time = @elapsed for t in tests
        @info("starting $t tests")
        test_time = @elapsed include(t)
        flush(stdout); flush(stderr)
        @info("finished $t tests in $(@sprintf("%8.2e seconds", test_time))")
        println()
        timings[t] = test_time
    end
    
    @info("finished all tests in $(@sprintf("%8.2e seconds", all_test_time))")
    println("\ntest suite timings (seconds):")
    display(timings)
    println()
    end
    
println()
;
    