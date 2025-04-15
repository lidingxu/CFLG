bash jobs.sh
#rm -r results/**/*
if [ -d "outputs" ]; then
  rm -r outputs/*
fi
lines=$(wc -l < job_list.txt)
echo "Number of jobs: $lines"
export GRB_LICENSE_FILE=/software/Gurobi/gurobi1201/gurobi.lic
export LC_ALL=C
export JULIA_DEPOT_PATH=".julia_depot"
julia --project=. -e 'ENV["CPLEX_STUDIO_BINARIES"] = "/software/opt-sw/cplex/bin/x86-64_linux/"; using Pkg; Pkg.instantiate(); Pkg.precompile()'

sbatch --array=0-$(($lines - 1)) run.slurm
