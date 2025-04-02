#!/bin/bash
#SBATCH --job-name=opt_array
#SBATCH --time=02:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=lidingxu.ac@gmail.com

# Set variables
timelimit=4600
timebound=4200
formulations=("EF" "EFP" "LEVFP" "LEFP" "LEFPI" "LEVFP" "LEVFPD" "LEVFPB" "LEVFPDB" "None")
covers=("Small" "Large")
solver="CPLEX" # "Gurobi", "CPLEX", "GLPK", "SCIP".
datapath="$PWD/benchmarks"
resultpath="$PWD/results"
juliabin="julia" #"/home/lxu/software/julia-1.10.2/bin/julia"


# Generate all job combinations and save to file
rm -f job_list.txt
for benchmark in $(ls ${datapath})
do
  if [ $benchmark == "test" ]
  then
    continue
  fi

  for instance in $(ls $datapath/$benchmark)
  do
    for formulation in "${formulations[@]}"
    do
      for cover in "${covers[@]}"
      do
        echo "$benchmark $instance $formulation $cover" >> job_list.txt
      done
    done
  done
done

# Get the specific job from the job list based on SLURM array index
if [ -n "$SLURM_ARRAY_TASK_ID" ]; then
  line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" job_list.txt)
  benchmark=$(echo $line | cut -d' ' -f1)
  instance=$(echo $line | cut -d' ' -f2)
  formulation=$(echo $line | cut -d' ' -f3)
  cover=$(echo $line | cut -d' ' -f4)

  echo "Running benchmark: $benchmark, instance: $instance, formulation: $formulation, cover: $cover"
  "$juliabin" src/main.jl "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "$solver" "$timelimit" "$formulation" "$cover"
else
  # Count total number of jobs
  total_jobs=$(wc -l < job_list.txt)
  echo "Total jobs to run: $total_jobs"

  # Submit the array job
  sbatch --array=1-$total_jobs $0
fi