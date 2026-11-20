

# cflg: an algorithmic toolkit for continuous set covering on networks

Facility location is an important application in operations research.  Continuous set covering on networks generalizes the classical dicrete set covering problems in graphs, which allows both the demands and facilities continuously locating in edges.
cflg contains several algorithms for solving continuous set covering on networks. cflg is written in Julia based on [JuMP](https://jump.dev/JuMP.jl/stable/installation/).


## Installation
cflg has been developped and tested under *Linux* System.

To use cflg, it requires:
- Julia,
- JuMP,

and one of the following MILP solvers:
- SCIP (default; SCIP.jl 0.12.8 with SCIP 10.0.3, installed with the package environment),
- Gurobi,
- CPLEX,
- GLPK.


## Benchmarks
There are seven benchmarks `city`, `Kgroup_A`, `Kgroup_B`, `random_A`, `random_B`, `tree_A`, `tree_B`, and a small `test` set. Each benchmark contains instances in `.txt` file format.

The almost-tree benchmarks `tree_A` and `tree_B` contain 36 instances each: one per generator (A, B) x number of vertices n in {100, 150, ..., 350} x proportion of extra edges p in {0.1, 0.15, 0.2, 0.3, 0.4, 0.5}. `python genalmosttree.py` regenerates them exactly (see `python genalmosttree.py --help` for other sizes, proportions and seeds).

Each instance file contains data of a graph. Its first line contains the number of nodes and edges of the graph; Edges are represented as a list from the second line: each line indicates an edge's end nodes and weight.

## Usage

Run one instance with
```
julia --project=. src/main.jl data_folder instance result_folder solver time_limit algorithm radius [seed]
```
where `solver` is `SCIP`, `Gurobi`, `CPLEX` or `GLPK`, `radius` is `Small` or `Large`, and `seed` sets the solver's random seed (0 = the solver's default). The result file is `result_folder/instance.algorithm.radius[.sseed]`; the solver log is written next to it with the suffix `.log`.


The algorithm option names and short descriptions (features are inferred from the bitmask flags):


| Algorithm | Model | Edge vs Edge-Vertex | Delimitation (P0) | Bound tighten (P1) | Valid Ineq. (V) | Disjunctive (D) | Long edge (L) | Indicator (I) | Attached pruning (A) | Comment |
|---|:---:|:-------------------:|:-----------------:|:------------------:|:---------------:|:-------:|:-----------:|:-------------:|:--------------------:|---------|
| `EF`     | Yes | Edge        | No  | No  | No  | No  | No  | No  | No  | Edge model (Fröhlich et al.) |
| `EFP0`   | Yes | Edge        | Yes | No  | No  | No  | No  | No  | No  | Edge model + delimitation (simple processing) |
| `EFP`    | Yes | Edge        | Yes | Yes | No  | No  | No  | No  | No  | Edge big‑M + delimitation + bound tightening (EF-DB) |
| `EVF`    | Yes | Edge‑Vertex | No  | No  | No  | No  | No  | No  | No  | Edge‑vertex big‑M formulation |
| `EVFP`   | Yes | Edge‑Vertex | Yes | Yes | No  | No  | No  | No  | No  | EVF + processing (bound tightening + delimitation) |
| `EVFPV`  | Yes | Edge‑Vertex | Yes | Yes | Yes | No  | No  | No  | No  | EVFP + simple valid inequalities |
| `LEVFP`  | Yes | Edge‑Vertex | Yes | Yes | No  | No  | Yes | No  | No  | Long edge‑vertex big‑M + processing (EVF-DLB in the paper) |
| `LEFPI`  | Yes | Edge        | Yes | Yes | No  | No  | Yes | Yes | No  | Edge model with long‑edge + indicator |
| `LEFP`   | Yes | Edge        | Yes | Yes | No  | No  | Yes | No  | No  | Long‑edge edge model + processing (EF-DLB) |
| `LEFPA`  | Yes | Edge        | Yes | Yes | No  | No  | Yes | No  | Yes | LEFP + attached-tree pruning (EF-PB) |
| `LEFPAI` | Yes | Edge        | Yes | Yes | No  | No  | Yes | Yes | Yes | EF-PB with indicator constraints (EF-PI) |
| `LEFPAD` | Yes | Edge        | Yes | Yes | No  | Yes | Yes | No  | Yes | EF-PB with the disjunctive formulation (EF-PE) |
| `LEFPAV` | Yes | Edge        | Yes | Yes | Yes | No  | Yes | No  | Yes | EF-PB + rank-1 lift-and-project cuts separated at the root node by a user-cut callback (EF-PBC) |
| `LEFPAN` | Yes | Edge        | Yes | Yes | No  | No  | Yes | No  | Yes | EF-PB with an empty user-cut callback, the baseline of EF-PBC (EF-PBN) |
| `LEFPAF`, `LEFPAIF`, `LEFPADF`, `LEFPAVF` | Yes | Edge | Yes | Yes | as above | as above | Yes | as above | Yes | EF-PB, EF-PI, EF-PE, EF-PBC with the dominance filter on disjunction pairs (EF-PBF, EF-PIF, EF-PEF, EF-PBCF) |
| `LEFPAU` | Yes | Edge        | Yes | Yes | No  | No  | Yes | No  | Yes | EF-PB without the vertex assignment: every long edge keeps both internal endpoints (EF-PBU) |
| `LEFPD`  | Yes | Edge        | Yes | Yes | No  | Yes | Yes | No  | No  | LEFP + disjunctive formulation |
| `LEFPV`  | Yes | Edge        | Yes | Yes | Yes | No  | Yes | No  | No  | LEFP + valid inequalities |
| `None`   | No  | —           | No  | No  | No  | No  | No  | No  | No  | Not a model — compute statistics (original, degree‑2‑free, subdivided) |

Use the exact algorithm name (e.g., `LEFPA`, `None`) when invoking the program.


## Experiments on a Slurm cluster

The revision experiments run as Slurm job arrays, one run (instance, algorithm, radius, seed) per task:
```
bash runjobs.sh                  # SCIP (default)
SOLVER=Gurobi bash runjobs.sh    # Gurobi
```
- `jobs.sh` lists the runs: every algorithm of the paper on city, Kgroup, random and the tree instances with p <= 0.2, with solver seeds 0, 1, 2 and a 7200 s time limit (group `main`), and EVF-DLB, EF-DLB, EF-PB, EF-PBF on the tree instances with p >= 0.3 (group `tree`). Environment variables `JOB_GROUPS`, `SEEDS`, `TREE_SEEDS` and `TREE_SIZES` restrict the list (see the header of `jobs.sh`).
- Results go to `results_rev/<solver>/<benchmark>/`. Runs whose result file already exists are skipped, so running `runjobs.sh` again submits only the missing runs.
- `runjobs.sh` instantiates the Julia environment in `.julia_depot`, then submits the list in arrays of at most `CHUNK` tasks (default 1000) with at most `MAXRUN` running at once (default 200).
- `run.slurm` gives each task 1 CPU, 16 GB and 2 h 20 min (the 7200 s solver limit plus Julia start-up and preprocessing). It requests partition `big`, one thread per physical core (`--hint=nomultithread`) and nodes with an Intel Xeon Gold 6342 (`--constraint=Gold6342`), so that all runs use the same processor; each Slurm output records the host and CPU model. Partition and feature names are site-specific: adapt them to your cluster (`sinfo -o "%P %N %f"`).
- The time limit is wall-clock time for both SCIP and Gurobi.




## References

If you find cflg useful in your work, we kindly request that you cite the following paper draft ([arXiv preprint](https://arxiv.org/abs/2203.00284)), which is recommended reading for advanced users:

    @misc{pelegrín2022continuous,
        title={Continuous Covering on Networks: Strong Mixed Integer Programming Formulations},
        author={Mercedes Pelegrín and Liding Xu},
        year={2022},
        eprint={2203.00284},
        archivePrefix={arXiv},
        primaryClass={math.OC}
    }

