#!/bin/bash
timelimit=1800
timebound=2200
formulations=("EF" "EFP" "EFPV" "EFPV2" "EFPI" "EFPD" "EVFP" "LEVFP")
covers=("Small" "Large")
solver="CPLEX" # "Gurobi", "CPLEX", "GLPK", "SCIP"
datapath="benchmarks"
resultpath="results"
gnuparalleltest=1 # enable: 1, disable: 0



runInstance() {
    benchmark_dir=$1
    solver=$2
    timelimit=$3
    result_dir=$4
    instance=$5
    formulation=$6
    cover=$7

    echo "$instance" "$formulation" "$cover"
    
    # example: julia  src/main.jl "benchmarks/test" "city_132.txt" "results/test" "SCIP" "100" "EF" "Small"
    julia  src/main.jl "$benchmark_dir" "$instance" "$result_dir" "$solver" "$timelimit"  "$formulation" "$cover"
                
}
export -f runInstance


benchmarks=$(ls ${datapath})

echo $benchmarks



for benchmark in $benchmarks
do
    if [ $benchmark == "test" ]
    then
        continue
    fi
    
    instances=$(ls $datapath/$benchmark)

    if [ $gnuparalleltest == 0 ]
    then
        for instance in  $instances
        do
            for formulation in ${formulations[@]}
                do
                for cover in ${covers[@]}
                do
                    timeout $timebound runInstance  "$datapath/$benchmark" "$solver" "$timelimit" "$resultpath/$benchmark"  "$instance"  "$formulation"  "$cover"
                done
            done
        done
    else
        parallel --will-cite --jobs 80% --timeout $timebound runInstance  "$datapath/$benchmark" "$solver" "$timelimit" "$resultpath/$benchmark"  ::: "$instances" :::  "${formulations[@]}" :::  "${covers[@]}"
    fi
done

