#!/bin/bash
# Generate the job list for the revision experiments: one line per run
#   datapath instance resultpath solver timelimit formulation cover [seed]
# (no seed for the graph-statistics runs "None"). Runs whose result file already exists are skipped, so the
# same command lists only the missing runs when resubmitting.
#
# Instances (all in benchmarks/):
#   tree instances: one per generator x n in {100,...,350} x p in {0.1, 0.15, 0.2, 0.3, 0.4, 0.5}
#   main group: city, Kgroup, random and the core tree instances (p in {0.1, 0.15, 0.2}, the design of the
#               paper's tree sets) -> all main formulations
#   tree group: the tree instances with p in {0.3, 0.4, 0.5} -> tree-study formulations (additional network
#               study; it also uses the core tree instances, whose runs come from the main group)
#
# usage: bash jobs.sh [job_list_file]
#   JOB_GROUPS="main tree"      which job groups to list (default both)
#   TREE_SIZES="100 200 300"    restrict the tree group to these sizes n (default all)
#   SEEDS="0 1 2"               solver seeds of the main group (seed 0 = Gurobi's default)
#   TREE_SEEDS="0 1 2"          solver seeds of the tree group
joblist=${1:-job_list.txt}
timelimit=7200
solver="Gurobi" # "Gurobi", "CPLEX", "GLPK", "SCIP".
covers=("Small" "Large")
main_seeds=(${SEEDS:-0 1 2})
tree_seeds=(${TREE_SEEDS:-0 1 2})
groups=(${JOB_GROUPS:-main tree})
datapath="$PWD/benchmarks"
resultpath="$PWD/results_rev"

# main: Part I (EVF-DLB, EF, EF-DB, EF-DLB, EF-PB), Part II (EF-PI, EF-PE, EF-PBC),
# EF-PBN (null-callback baseline of EF-PBC) and the dominance-filter variants (EF-PBF, EF-PIF, EF-PEF, EF-PBCF)
main_formulations=("LEVFP" "EF" "EFP" "LEFP" "LEFPA" "LEFPAI" "LEFPAD" "LEFPAV" "LEFPAN" "LEFPAF" "LEFPAIF" "LEFPADF" "LEFPAVF")
# small radius only: EF-PBU (EF-PB without vertex-assignment; at large radius the models are identical)
main_small_formulations=("LEFPAU")
# tree study: EVF-DLB (SOTA), EF-DLB (EF-PB without attached-tree projection), EF-PB, EF-PBF
tree_formulations=("LEVFP" "LEFP" "LEFPA" "LEFPAF")
tree_sizes=" ${TREE_SIZES:-100 150 200 250 300 350} "

# tree instance <prefix>.<n>.<extra>.s<seed> is a core instance if p = extra/n <= 0.2
is_core_tree() {
  local n extra
  IFS=. read -r _ n extra _ <<< "$1"
  [[ $((5 * extra)) -le $n ]]
}

rm -f "$joblist"
total=0
skipped=0
add_job() { # benchmark instance formulation cover [seed]
  local out="$resultpath/$1/$2.$3.$4${5:+.s$5}"
  if [[ -f "$out" ]]; then
    skipped=$((skipped + 1))
    return
  fi
  echo "$datapath/$1" "$2" "$resultpath/$1" "$solver" "$timelimit" "$3" "$4" $5 >> "$joblist"
  total=$((total + 1))
}

add_runs() { # benchmark instance group formulations...
  local benchmark=$1 instance=$2 group=$3; shift 3
  local seeds=("${main_seeds[@]}")
  [[ $group == "tree" ]] && seeds=("${tree_seeds[@]}")
  for cover in "${covers[@]}"; do
    add_job "$benchmark" "$instance" "None" "$cover"
    local extra=()
    [[ $group == "main" && $cover == "Small" ]] && extra=("${main_small_formulations[@]}")
    for formulation in "$@" "${extra[@]}"; do
      for seed in "${seeds[@]}"; do
        add_job "$benchmark" "$instance" "$formulation" "$cover" "$seed"
      done
    done
  done
}

for group in "${groups[@]}"; do
  [[ $group == "main" || $group == "tree" ]] || { echo "unknown group $group"; exit 1; }
done
for benchmark in $(ls "$datapath"); do
  [[ $benchmark == "test" ]] && continue
  for instance in $(ls "$datapath/$benchmark"); do
    if [[ $benchmark != tree_* ]] || is_core_tree "$instance"; then
      [[ " ${groups[*]} " == *" main "* ]] && add_runs "$benchmark" "$instance" main "${main_formulations[@]}"
    else
      n=$(echo "$instance" | cut -d. -f2)
      [[ " ${groups[*]} " == *" tree "* && $tree_sizes == *" $n "* ]] && add_runs "$benchmark" "$instance" tree "${tree_formulations[@]}"
    fi
  done
done
echo "$total jobs listed in $joblist ($skipped already done, skipped)"
