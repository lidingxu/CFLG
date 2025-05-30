#!/bin/bash
#SBATCH --job-name=distributed_jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=02:08:00
#SBATCH --output=outputs/output_%A_%a.out
#SBATCH --partition=opt_int

PARAM=$(sed -n "$(($SLURM_ARRAY_TASK_ID + 1))p" job_list.txt)

# Split PARAM into variables
read datapath instance resultpath solver timelimit formulation cover <<< "$PARAM"

echo "Data path: $datapath"
echo "Instance: $instance"
echo "Result path: $resultpath"
echo "Solver: $solver"
echo "Time limit: $timelimit"
echo "Formulation: $formulation"
echo "Cover: $cover"

echo "julia" "--project=."  "src/main.jl" "$datapath" "$instance" "$resultpath" "$solver" "$timelimit" "$formulation" "$cover"
srun julia --project=. src/main.jl "$datapath" "$instance" "$resultpath" "$solver" "$timelimit" "$formulation" "$cover"
