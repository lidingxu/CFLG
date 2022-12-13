#!/bin/bash
timelimit=1800
timebound=2200
algorithms=("EF" "SF", "RF"  "ESF" "ERF" "EDF" "ESFV" "None")
covers=("Small" "Large")
solver="CPLEX"
datapath=$(cd ../benchmarks;pwd)
resultpath=$(cd ../results;pwd)
gnuparalleltest=1



runInstance() {
    benchmark_dir=$1
    solver=$2
    timelimit=$3
    result_dir=$4
    instance=$5
    algo=$6
    cover=$7

    echo "$instance" "$algo" "$cover"
    python ./checkexec.py  $result_dir $instance $algo $cover
    if [ $? == 1 ]
    then
        return 1
    fi

    julia ./runbenchmark.jl $benchmark_dir $solver "$timelimit" "$result_dir" "$instance" "$algo" "$cover"
                
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
            for algo in ${algorithms[@]}
                do
                for cover in ${covers[@]}
                do
                    #continue
                    #echo "$instance" "$algo" "$cover"
                    timeout $timebound runInstance  "$datapath/$benchmark" "$solver" "$timelimit" "$resultpath/$benchmark"  "$instance"  "$algo"  "$cover"
                done
            done
        done
    else
        parallel --will-cite --jobs 88% --timeout $timebound runInstance  "$datapath/$benchmark" "$solver" "$timelimit" "$resultpath/$benchmark"  ::: "$instances" :::  "${algorithms[@]}" :::  "${covers[@]}"
        #parallel --will-cite --jobs 37% julia ./runbenchmark.jl  "$datapath/$benchmark" "CPLEX" "$timelimit" "$resultpath/$benchmark"  ::: "$instances" :::  "${algorithms[@]}" :::  "${covers[@]}"
        #$instances | parallel --will-cite   --dryrun  "printls {}"
        #parallel --will-cite  printls0 para ::: 1
        #parallel --will-cite  printls "$benchmark" para  ::: "$instances" :::  "${algorithms[@]}"
        #break 
        #parallel --will-cite -j 4 --dryrun  printls ::: "${instances[@]}" :::  "${algorithms[@]}"  
        #parallel --will-cite  runInstance ::: "$algorithms"  "${instances[@]}"  "$benchmark" "para"
    fi
    #find   $datapath/$benchmark -name *.cbp
done

