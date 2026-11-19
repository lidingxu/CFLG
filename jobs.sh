#!/bin/bash
# Generate the job list for the revision experiments: one line per run
#   datapath instance resultpath solver timelimit formulation cover [seed]
# (no seed for the graph-statistics runs "None"). Runs whose result file already exists are skipped, so the
# same command lists only the missing runs when resubmitting.
#
# usage: bash jobs.sh [job_list_file]
#   JOB_GROUPS="main tree"      which job groups to list (default both)
#   TREE_SIZES="100 200 300"    restrict the tree study to these sizes n (default all)
#   SEEDS="0 1 2"               solver seeds (seed 0 = Gurobi's default)
joblist=${1:-job_list.txt}
timelimit=7200
solver="Gurobi" # "Gurobi", "CPLEX", "GLPK", "SCIP".
covers=("Small" "Large")
seeds=(${SEEDS:-0 1 2})
groups=(${JOB_GROUPS:-main tree})
resultpath="$PWD/results_rev"

# main benchmark: Part I (EVF-DLB, EF, EF-DB, EF-DLB, EF-PB), Part II (EF-PI, EF-PE, EF-PBC),
# EF-PBN (null-callback baseline of EF-PBC) and the dominance-filter variants (EF-PBF, EF-PIF, EF-PEF, EF-PBCF)
main_formulations=("LEVFP" "EF" "EFP" "LEFP" "LEFPA" "LEFPAI" "LEFPAD" "LEFPAV" "LEFPAN" "LEFPAF" "LEFPAIF" "LEFPADF" "LEFPAVF")
main_datapath="$PWD/benchmarks"
# tree study: EVF-DLB (SOTA), EF-DLB (EF-PB without attached-tree projection), EF-PB, EF-PBF
tree_formulations=("LEVFP" "LEFP" "LEFPA" "LEFPAF")
tree_datapath="$PWD/benchmarks_tree"
tree_sizes=" ${TREE_SIZES:-100 150 200 250 300 350} "

rm -f "$joblist"
total=0
skipped=0
add_job() { # datapath instance resultdir formulation cover [seed]
  local out="$3/$2.$4.$5${6:+.s$6}"
  if [[ -f "$out" ]]; then
    skipped=$((skipped + 1))
    return
  fi
  echo "$1" "$2" "$3" "$solver" "$timelimit" "$4" "$5" $6 >> "$joblist"
  total=$((total + 1))
}

list_group() { # datapath formulations...
  local datapath=$1; shift
  local formulations=("$@")
  for benchmark in $(ls "$datapath"); do
    [[ $benchmark == "test" ]] && continue
    for instance in $(ls "$datapath/$benchmark"); do
      if [[ $datapath == "$tree_datapath" ]]; then
        n=$(echo "$instance" | cut -d. -f2)
        [[ $tree_sizes == *" $n "* ]] || continue
      fi
      for cover in "${covers[@]}"; do
        add_job "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "None" "$cover"
        for formulation in "${formulations[@]}"; do
          for seed in "${seeds[@]}"; do
            add_job "$datapath/$benchmark" "$instance" "$resultpath/$benchmark" "$formulation" "$cover" "$seed"
          done
        done
      done
    done
  done
}

for group in "${groups[@]}"; do
  case $group in
    main) list_group "$main_datapath" "${main_formulations[@]}" ;;
    tree) list_group "$tree_datapath" "${tree_formulations[@]}" ;;
    *) echo "unknown group $group"; exit 1 ;;
  esac
done
echo "$total jobs listed in $joblist ($skipped already done, skipped)"
