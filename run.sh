#!/bin/bash
timelimit=1800
timebound=2200
formulations=("EF" "EFP" "EFPV" "EFPV2" "EFPI" "EFPD" "EVFP" "LEVFP" "None")
covers=("Small" "Large")
solver="CPLEX" # "Gurobi", "CPLEX", "GLPK", "SCIP".
datapath="$PWD/benchmarks"
resultpath="$PWD/results"
gnuparalleltest=0 # enable: 1, disable: 0
juliabin="julia" #"/home/lxu/software/julia-1.10.2/bin/julia"


runInstance() {
    juliabin=$1
    benchmark_dir=$2
    solver=$3
    timelimit=$4
    result_dir=$5
    instance=$6
    formulation=$7
    cover=$8
    "
    # example: julia  src/main.jl "benchmarks/test" "K400.10.con.red" "results/test" "CPLEX" "100" "LEFP" "Small"
    "$juliabin"  src/main.jl "$benchmark_dir" "$instance" "$result_dir" "$solver" "$timelimit"  "$formulation" "$cover"
                
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
                    echo  "$juliabin" "src/main.jl" "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "$solver" "$timelimit"  "$formulation" "$cover"
                    timeout $timebound "$juliabin" src/main.jl "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "$solver" "$timelimit"  "$formulation" "$cover"
                    exit
                done
            done
        done
    else
        parallel --will-cite --jobs 80% --timeout $timebound runInstance  "$juliabin" "$datapath/$benchmark" "$solver" "$timelimit" "$resultpath/$benchmark"  ::: "$instances" :::  "${formulations[@]}" :::  "${covers[@]}"
    fi
done

