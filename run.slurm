#!/bin/bash
#SBATCH --job-name=distributed_jobs
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --mem=8G
#SBATCH --time=02:20:00
#SBATCH --output=outputs/output_%A_%a.out
#SBATCH --partition=big
#SBATCH --hint=nomultithread
#SBATCH --constraint=Gold6342

# JOBLIST: job list file, OFFSET: index of the first line of this array chunk (both set by runjobs.sh)
JOBLIST=${JOBLIST:-job_list.txt}
OFFSET=${OFFSET:-0}
PARAM=$(sed -n "$(($SLURM_ARRAY_TASK_ID + $OFFSET + 1))p" "$JOBLIST")

# Split PARAM into variables (seed is empty for the graph-statistics runs "None")
read datapath instance resultpath solver timelimit formulation cover seed <<< "$PARAM"

echo "Data path: $datapath"
echo "Instance: $instance"
echo "Result path: $resultpath"
echo "Solver: $solver"
echo "Time limit: $timelimit"
echo "Formulation: $formulation"
echo "Cover: $cover"
echo "Seed: $seed"
echo "Host: $(hostname), CPU: $(lscpu | grep 'Model name' | sed 's/Model name: *//')"

# julia --project=. src/main.jl benchmarks/test testg.txt none CPLEX 20 LEFPA Small 0
echo "julia" "--project=."  "src/main.jl" "$datapath" "$instance" "$resultpath" "$solver" "$timelimit" "$formulation" "$cover" $seed
srun julia --project=. src/main.jl "$datapath" "$instance" "$resultpath" "$solver" "$timelimit" "$formulation" "$cover" $seed
