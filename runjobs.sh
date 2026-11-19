# Submit the revision experiments: lists the missing runs (jobs.sh) into a timestamped job list and submits it
# as Slurm arrays of at most CHUNK tasks with at most MAXRUN running at once per array.
#   CHUNK=1000 MAXRUN=200 JOB_GROUPS="main tree" TREE_SIZES="100 200 300" bash runjobs.sh
# Rerunning the script resubmits only runs without a result file. Old Slurm outputs are kept.
CHUNK=${CHUNK:-1000}
MAXRUN=${MAXRUN:-200}
mkdir -p jobs outputs
joblist="jobs/job_list_$(date +%Y%m%d_%H%M%S).txt"
bash jobs.sh "$joblist"
lines=$(wc -l < "$joblist")
echo "Number of jobs: $lines"
[ "$lines" -eq 0 ] && exit 0
export GRB_LICENSE_FILE=/software/Gurobi/gurobi1201/gurobi.lic
export LC_ALL=C
export JULIA_DEPOT_PATH=".julia_depot"
julia --project=. -e 'ENV["CPLEX_STUDIO_BINARIES"] = "/software/opt-sw/cplex/bin/x86-64_linux/"; using Pkg; Pkg.instantiate(); Pkg.precompile()'

for ((offset = 0; offset < lines; offset += CHUNK)); do
  n=$(( lines - offset < CHUNK ? lines - offset : CHUNK ))
  sbatch --array=0-$((n - 1))%$MAXRUN --export=ALL,JOBLIST="$joblist",OFFSET=$offset run.slurm
done
