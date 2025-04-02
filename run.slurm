#!/bin/bash
#SBATCH --job-name=distributed_jobs
#SBATCH --array=0-$(($(wc -l < job_list.txt) - 1))  # 任务数 = job_list.txt 行数
#SBATCH --ntasks=1   # 每个任务运行 1 个进程
#SBATCH --cpus-per-task=1  # 每个任务使用 2 个 CPU 核
#SBATCH --nodes=40   # 20 个节点
#SBATCH --output=output_%A_%a.out

# 读取当前任务参数 (基于 SLURM_ARRAY_TASK_ID)
PARAM=$(sed -n "$(($SLURM_ARRAY_TASK_ID + 1))p" job_list.txt)

echo "Running task $SLURM_ARRAY_TASK_ID on $(hostname) with params: $PARAM"

# Split PARAM into variables
read datapath instance resultpath solver timelimit formulation cover <<< "$PARAM"

echo "Data path: $datapath"
echo "Instance: $instance"
echo "Result path: $resultpath"
echo "Solver: $solver"
echo "Time limit: $timelimit"
echo "Formulation: $formulation"
echo "Cover: $cover"

# 运行你的程序 (替换成实际的命令)
julia src/main.jl "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "$solver" "$timelimit" "$formulation" "$cover"
