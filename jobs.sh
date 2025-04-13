#!/bin/bash
# Set variables
timelimit=3600
timebound=4200
formulations=("LEFPD" "LEFPDB" "EF" "EFP" "LEFP" "LEFPI" "LEVFP" "None")
covers=("Small" "Large")
solver="CPLEX" # "Gurobi", "CPLEX", "GLPK", "SCIP".
datapath="$PWD/benchmarks"
resultpath="$PWD/results"
juliabin="julia" #"/home/lxu/software/julia-1.10.2/bin/julia"


# Generate all job combinations and save to file
rm -f job_list.txt
for benchmark in $(ls ${datapath})
do
  if [ $benchmark != "random_B" ]
  then
    continue
  fi

  for instance in $(ls $datapath/$benchmark)
  do
    for formulation in "${formulations[@]}"
    do
      for cover in "${covers[@]}"
      do
        echo "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "$solver" "$timelimit" "$formulation" "$cover" >> job_list.txt
      done
    done
  done
done
